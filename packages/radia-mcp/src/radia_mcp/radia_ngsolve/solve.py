"""Reusable NGSolve magnetostatic A-formulation solver -- the validated,
high-mu-safe foundation that the force / energy extractors (``force.py``) and the
independent cross-validation tests build on.

#25 safe: put high-permeability bodies DIRECTLY in the surrounding air (no nested
"shell" material) so the interior B is correct; see the ``force_validation`` MCP
tool. Cross-validated against an independent reference solve (sphere field
0.11 %, solenoid B <0.35 %, inductance 0.01 %, ...).
"""
import cmath
import math

from ngsolve import (HCurl, H1, Periodic, NumberSpace, FESpace, BilinearForm,
                     LinearForm, GridFunction, CoefficientFunction, InnerProduct,
                     curl, grad, dx, Integrate, Conj, Variation, Preconditioner,
                     x, y, sqrt, BND)

MU0 = 4.0e-7 * math.pi
NU0 = 1.0 / MU0


def _direct_inverse(matrix, freedofs, inverse, where=""):
    """Direct sparse factorization with an actionable message on the 3D curl-curl
    integer-overflow ceiling.

    NGSolve's sparsecholesky/umfpack factorisation of a 3D HCurl (curl-curl) system
    overflows around ~84k order-2 elements, surfacing as a cryptic
    "bad array new length". That is INTEGER OVERFLOW in the factorisation fill-in,
    NOT out-of-RAM -- it triggers with tens of GB still free. Re-raise with the fix
    rather than letting the opaque message escape. (Bug -> permanent guard, runtime
    layer: see the note in rules.py for why this is not a static lint.)
    """
    try:
        return matrix.Inverse(freedofs, inverse=inverse)
    except Exception as e:
        msg = str(e).lower()
        if "array new length" in msg:
            raise RuntimeError(
                "Direct solver '{}' overflowed the factorization index space{} -- "
                "this is the ~84k order-2 3D HCurl (curl-curl) ceiling: integer "
                "overflow in the sparse fill-in, NOT out-of-RAM. Switch to an "
                "iterative solver -- CG with a BDDC/AMG preconditioner "
                "(Preconditioner(a, 'bddc') registered BEFORE a.Assemble()) -- or "
                "coarsen the mesh / lower the element order.".format(
                    inverse, (" in " + where) if where else "")
            ) from e
        raise


def periodic_h1(mesh, order=2, dirichlet="", antiperiodic=False):
    """H1 space whose DOFs on a netgen-identified PERIODIC boundary pair are tied
    together -- the FEMM "periodic / anti-periodic" BC for motor / machine SECTOR
    models (solve one pole or slot pitch instead of the whole machine).

    The MESH must carry the periodic identification: in 2D netgen
    ``SplineGeometry`` append the SLAVE edge with ``copy=<master_segment>`` (and
    define the boundary loop CCW so the mesher does not stall); in OCC use
    ``edge.Identify(other, name, IdentificationType.PERIODIC, trafo)``. This helper
    only wraps the resulting ``H1`` with :class:`ngsolve.Periodic`.

    antiperiodic=True flips the sign across the seam (A_slave = -A_master) for a
    HALF-period sector (the usual motor case). It is implemented as a phase = -1
    identification, which requires a COMPLEX space, so the returned space is
    complex -- assemble/solve as usual and take ``gfu.real`` for the physical
    field.

    Drop-in for the ``H1(mesh, order, dirichlet=...)`` line in the planar solvers
    (e.g. :func:`solve_planar_magnetostatic`): the identification ties A_z across
    the sector cut, reproducing the cyclic machine from one sector. Validated to
    machine precision against a closed-form periodic strip
    (tests/test_planar_periodic.py: periodic 2e-7, anti-periodic 2e-8 rel L2).
    """
    base = H1(mesh, order=order, dirichlet=dirichlet, complex=antiperiodic)
    return Periodic(base, phase=[-1]) if antiperiodic else Periodic(base)


def reluctivity(mesh, mu_r_by_material):
    """Build a per-material reluctivity CF nu = NU0 / mu_r from a
    ``{material_name: mu_r}`` dict (materials not listed default to mu_r=1)."""
    return CoefficientFunction([NU0 / float(mu_r_by_material.get(m, 1.0))
                                for m in mesh.GetMaterials()])


def azimuthal_coil_current(mesh, total_current, section_area, region="coil"):
    """Divergence-free azimuthal current density J0 * phi_hat over ``region``
    (J0 = total_current / section_area), zero elsewhere. Returns a vector CF."""
    ind = CoefficientFunction([1.0 if m == region else 0.0
                               for m in mesh.GetMaterials()])
    j0 = total_current / section_area
    r_xy = sqrt(x * x + y * y + 1e-12)
    return ind * j0 * CoefficientFunction((-y, x, 0)) / r_xy


def remanent_source(mesh, Br_by_material):
    """Permanent-magnet source for the A-formulation: a HARD magnet (mu_r=1) with
    remanent flux density Br (vector, tesla) in a region contributes the linear
    form  int nu0 Br . curl(v)  there. Returns the vector CF ``nu0 * Br`` (zero
    outside the listed magnets) to hand to the ``curl_source`` argument of
    :func:`solve_magnetostatic_Aform`.

    ``Br_by_material`` maps material name -> (Brx, Bry, Brz) in tesla. Assumes the
    magnet mu_r = 1; for mu_r != 1 build ``nu Br`` and pass it as ``curl_source``.
    """
    zero = (0.0, 0.0, 0.0)
    comps = [Br_by_material.get(m, zero) for m in mesh.GetMaterials()]
    bx = CoefficientFunction([c[0] for c in comps])
    by = CoefficientFunction([c[1] for c in comps])
    bz = CoefficientFunction([c[2] for c in comps])
    return NU0 * CoefficientFunction((bx, by, bz))


def planar_magnet_source(mesh, magnets):
    """2D-planar (Cartesian A_z) permanent-magnet source, FEMM-convention.

    ``magnets`` maps material name -> (Hc, phi_deg): coercivity Hc [A/m] and
    magnetization direction phi measured from +x [deg] (FEMM MagDir). Returns
    the 2D vector CF ``nu*B_rem = Hc*(cos phi, sin phi)`` per magnet material
    (zero elsewhere), to feed :func:`solve_planar_magnetostatic`. This is the
    FEMM ``prob3big.cpp`` magnetization edge-loop written as a continuous source.
    """
    zero = (0.0, 0.0)
    comp = []
    for m in mesh.GetMaterials():
        if m in magnets:
            hc, phi = magnets[m]
            ph = math.radians(phi)
            comp.append((hc * math.cos(ph), hc * math.sin(ph)))
        else:
            comp.append(zero)
    sx = CoefficientFunction([c[0] for c in comp])
    sy = CoefficientFunction([c[1] for c in comp])
    return CoefficientFunction((sx, sy))


def laminated_mu_eff(mu_r, sigma, omega, d_lam, fill=1.0):
    """Complex effective permeability of IN-PLANE laminated steel (FEMM AC
    lamination model -- the ``Lamination & Wire Type`` material).

    A stack of laminations (each steel sheet of thickness ``d_lam``, relative
    permeability ``mu_r``, conductivity ``sigma``, stacking/fill factor
    ``fill``) excited by a field PARALLEL to the sheets develops eddy currents
    that circulate within each sheet -> a 1D skin effect across the lamination
    thickness.  Homogenized, the stack behaves as a single block of complex
    permeability

        mu_eff = mu0 * [ fill * mu_r * tanh(b)/b + (1 - fill) ],
        b      = (d_lam/2) * sqrt(1j*omega*mu0*mu_r*sigma)

    The ``tanh(b)/b`` factor is the lamination flux-exclusion + eddy loss
    (Im(mu_eff) < 0 carries the loss); the ``(1-fill)`` term is the
    non-conducting insulation volume fraction (mu0) in parallel.  ``omega=0`` or
    ``sigma=0`` returns the real static value ``mu0*(fill*mu_r + 1 - fill)``.

    Returns a complex python scalar.  Use ``1/laminated_mu_eff(...)`` as the
    (uniform, complex) reluctivity ``nu`` of the laminated region in
    :func:`solve_planar_eddy`, with ``sigma = 0`` there (the eddy loss is
    already captured by ``Im(mu_eff)`` -- do NOT also mesh-resolve the sheets).
    """
    if omega == 0 or sigma == 0:
        return MU0 * (fill * mu_r + (1.0 - fill))
    b = (d_lam / 2.0) * cmath.sqrt(1j * omega * MU0 * mu_r * sigma)
    factor = cmath.tanh(b) / b
    return MU0 * (fill * mu_r * factor + (1.0 - fill))


def laminated_reluctivity_tensor(mesh, lam_by_material, default_mu_r=1.0,
                                 lam_normal="y"):
    """Static ANISOTROPIC reluctivity matrix CF for laminated regions (FEMM
    lamination, DC / low-frequency).

    Laminations make the steel magnetically anisotropic via the stacking
    (fill) factor: flux PARALLEL to the sheets sees steel and insulation in
    parallel, flux NORMAL to the stack sees them in series:

        mu_par  = mu0 * (fill*mu_r + 1 - fill)        (along the sheets)
        mu_perp = mu0 / (fill/mu_r + 1 - fill)        (across the stack)

    ``lam_by_material`` maps ``material -> (mu_r, fill)``.  ``lam_normal`` is the
    stacking direction (the sheet normal), ``"x"`` or ``"y"``.  Materials not
    listed are isotropic with ``default_mu_r``.

    Returns a 2x2 matrix CoefficientFunction to pass as ``nu`` to
    :func:`solve_planar_magnetostatic` (whose ``nu*grad(u)*grad(v)`` term reads
    a matrix ``nu`` as the anisotropic form ``(nu grad A).grad v``).  Because
    ``B = (dA/dy, -dA/dx)``, the matrix diagonal is the reluctivity felt by
    ``(B from dA/dx, B from dA/dy)`` -- i.e. swapped relative to B, which this
    helper handles: ``lam_normal="y"`` -> diag(1/mu_perp, 1/mu_par).
    """
    mxx, myy = [], []
    for m in mesh.GetMaterials():
        if m in lam_by_material:
            mu_r, fill = lam_by_material[m]
            mu_par = MU0 * (fill * mu_r + (1.0 - fill))
            mu_perp = MU0 / (fill / mu_r + (1.0 - fill))
            nu_par, nu_perp = 1.0 / mu_par, 1.0 / mu_perp
            if lam_normal == "y":      # stack along y: B_y normal, B_x parallel
                mxx.append(nu_perp); myy.append(nu_par)
            else:                       # stack along x: B_x normal, B_y parallel
                mxx.append(nu_par); myy.append(nu_perp)
        else:
            nu_iso = NU0 / float(default_mu_r)
            mxx.append(nu_iso); myy.append(nu_iso)
    Mxx = CoefficientFunction(mxx)
    Myy = CoefficientFunction(myy)
    return CoefficientFunction((Mxx, 0.0, 0.0, Myy), dims=(2, 2))


def stranded_source(mesh, currents_by_region):
    """FEMM "stranded" conductor source: a UNIFORM current density Jz = I/area in
    each region (the conductor is finely stranded / litz, so eddy currents cannot
    redistribute -- the current stays uniform and Rac = Rdc, no skin/proximity in
    the bundle itself).

    ``currents_by_region`` maps material -> total current I [A].  Returns the Jz
    CoefficientFunction [A/m^2] to feed :func:`solve_planar_eddy` (or the static
    solver) WITH ``sigma = 0`` in those regions (the stranded conductor carries
    an imposed current, not an eddy reaction).  The total current per region is
    preserved exactly: int Jz dA = I.

    Contrast with a SOLID conductor (driven_region/total_current and sigma>0),
    whose current redistributes (skin effect) -- see tests/test_stranded.py.
    """
    comp = []
    for m in mesh.GetMaterials():
        if m in currents_by_region:
            area = Integrate(CoefficientFunction(1.0), mesh,
                             definedon=mesh.Materials(m))
            comp.append(currents_by_region[m] / area)
        else:
            comp.append(0.0)
    return CoefficientFunction(comp)


def coil_flux_linkage_2d(A, mesh, pos_region, neg_region, n_turns, depth=1.0):
    """2D flux linkage [Wb] of an ``n_turns`` search coil whose 'go' side lies in
    ``pos_region`` and 'return' side in ``neg_region`` (named mesh materials), from
    the planar A_z solution ``A``:

        lambda = n_turns * depth * ( <A_z>_pos - <A_z>_neg ),

    ``<A_z>`` the AREA-AVERAGE of A_z over each coil side -- exactly FEMM's circuit
    "flux linkage" (``mo_getcircuitproperties``) for an N-turn coil. The motor
    back-EMF follows for rotation at omega: ``e(theta) = -d(lambda)/dt =
    -omega * d(lambda)/d(theta)`` (read lambda(theta) over a rotor-angle sweep and
    differentiate). Validated on a 2-pole PMSM against an independent reference
    (lambda(theta) sinusoid, amplitude + shape)."""
    def mean_over(reg):
        dom = dx(definedon=mesh.Materials(reg))
        area = Integrate(CoefficientFunction(1.0) * dom, mesh)
        return Integrate(A * dom, mesh) / area
    return n_turns * depth * (mean_over(pos_region) - mean_over(neg_region))


def inductance_2d(A, mesh, pos_region, neg_region, n_turns, current, depth=1.0):
    """2D coil INDUCTANCE [H] from a current-driven A_z solution: L = lambda/current,
    lambda = :func:`coil_flux_linkage_2d`.  For synchronous-machine Ld/Lq, drive the
    stator coil along the rotor d- or q-axis and re-solve; the saliency gives Ld != Lq.
    For a SATURATED machine use the frozen-permeability dq method (standard FE technique):
    converge a nonlinear solve at the operating point, FREEZE nu(|B|) to a fixed CF, then
    take these linear dq-perturbation inductances (superposition only holds with frozen nu).

    UNIT caveat (same as coil_flux_linkage_2d): a mm-meshed solve has A_z = 1e3*A_phys,
    so multiply the result by 1e-6 for henries of a 1 mm stack (1e-3 m depth, A /1e3)."""
    return coil_flux_linkage_2d(A, mesh, pos_region, neg_region, n_turns, depth) / current


def saturated_secant_inductance(mesh, nu_of_B, Jz, pos_region, neg_region, n_turns,
                                current, depth=1.0, **solve_kw):
    """SECANT (apparent) winding inductance L = lambda(I)/I at a SATURATED operating
    current: a nonlinear ``solve_planar_magnetostatic_nonlinear`` at the given current,
    then ``coil_flux_linkage_2d / current``. As the iron saturates with current, L falls
    below its unsaturated (linear) value -- the Ld(i) saturation knee that underlies
    saturated dq / MTPA maps. ``nu_of_B`` is the callable reluctivity (smooth bounded so
    Picard converges); ``**solve_kw`` -> relax/max_iter/tol/order of the nonlinear solve.
    Validated examples/comsol_class/saturating_inductance.py."""
    A = solve_planar_magnetostatic_nonlinear(mesh, nu_of_B, Jz=Jz, **solve_kw)
    return coil_flux_linkage_2d(A, mesh, pos_region, neg_region, n_turns, depth) / current


def incremental_inductance(lam_plus, lam_minus, dI):
    """Differential (INCREMENTAL) inductance L_inc = dlambda/dI by central finite difference,
    ``(lam_plus - lam_minus)/(2 dI)``, from nonlinear flux linkages at I0 +/- dI. For a
    SATURATING (concave) lambda(I), L_inc < L_sec = lambda/I and both fall as the iron
    saturates -- the small-signal / control inductance (vs the apparent
    :func:`saturated_secant_inductance`). NB the PROPER frozen-permeability L_inc uses the
    tangent differential-reluctivity TENSOR dH/dB; this FD is the ground truth."""
    return (lam_plus - lam_minus) / (2.0 * dI)


def incremental_inductance_matrix(flux_linkages, currents, di):
    """INCREMENTAL (small-signal) inductance MATRIX L_jk = d lambda_j / d i_k of a multi-winding
    (or dq) machine, by CENTRAL finite difference of the nonlinear flux-linkage map -- the
    multi-port generalisation of :func:`incremental_inductance`.

    ``flux_linkages`` is a callable: current-vector (list) -> list of flux linkages
    ``[lambda_0, lambda_1, ...]`` (e.g. one nonlinear FE solve per call). ``currents`` is the
    operating-point vector; ``di`` the perturbation. Returns the n x n matrix (list of lists),
    column k = d(all lambda)/d i_k from +-di on i_k (so n calls of 2 solves each).

    The matrix is SYMMETRIC -- L_jk = L_kj (RECIPROCITY, from the magnetic co-energy Hessian
    d^2 W'/(di_j di_k)); the OFF-DIAGONALS are the cross-(saturation) inductances (e.g. L_dq,
    which is ZERO without saturation but grows as a shared iron path saturates). The frozen-
    permeability APPARENT/secant inductance (:func:`saturated_secant_inductance`, the freeze of
    nu at the operating |B|) is LARGER than this incremental value in saturation -- use THIS
    tangent matrix for small-signal / control / dq-current-loop design. Validated
    examples/comsol_class/cross_saturation_tensor.py (reciprocity + secant-vs-incremental)."""
    n = len(currents)
    L = [[0.0] * n for _ in range(n)]
    for k in range(n):
        ip, im = list(currents), list(currents)
        ip[k] += di
        im[k] -= di
        lam_p, lam_m = flux_linkages(ip), flux_linkages(im)
        for j in range(n):
            L[j][k] = (lam_p[j] - lam_m[j]) / (2.0 * di)
    return L


def magnetic_coenergy(currents, fluxes):
    """Magnetic CO-ENERGY W'(I) = int_0^I lambda dI' (cumulative trapezoidal) from a sampled
    lambda(I) curve (``currents`` ascending from ~0). Returns the cumulative co-energy at each
    sample. The stored energy is W = lambda*I - W'; for a SATURATING (concave) lambda(I),
    W' > W. Co-energy gives the force F = dW'/dx at constant current (the energy method)."""
    cur = [0.0] + list(currents)
    flx = [0.0] + list(fluxes)
    out, cum = [], 0.0
    for i in range(1, len(cur)):
        cum += 0.5 * (flx[i] + flx[i - 1]) * (cur[i] - cur[i - 1])
        out.append(cum)
    return out


def pm_machine_constant(lam_plus, lam_minus, dtheta):
    """The PM-machine CONSTANT from the PM flux-linkage derivative (central finite
    difference), ``(lam_plus - lam_minus)/(2 dtheta)`` with ``dtheta`` in radians. By energy
    conservation the back-EMF constant equals the torque constant:

        Ke = dlambda_m/dtheta [Wb/rad = V s/rad]  ==  Kt = T/I [Nm/A]   (Kt = Ke),

    because a current I in a winding with PM flux linkage lambda_m(theta) feels the co-energy
    torque T = I dlambda_m/dtheta. Compute Ke here from the no-load back-EMF (PM-only coil
    flux linkage at theta +/- dtheta); the FE torque-per-amp Kt = (T(I) - T(0))/I (eggshell
    torque, current minus cogging) matches it. Validated examples/comsol_class/kt_ke.py."""
    return (lam_plus - lam_minus) / (2.0 * dtheta)


def frozen_reluctivity(A, nu_of_B):
    """FROZEN-PERMEABILITY reluctivity: snapshot nu(|B|) at the field of a CONVERGED
    nonlinear solution ``A`` (planar A_z), returning a FIXED reluctivity CF.

    With nu frozen at the operating point the problem becomes LINEAR, so SUPERPOSITION is
    restored -- the basis of saturated-machine Ld/Lq(id,iq) maps and PM/armature flux
    decomposition. Pass the result as the ``nu`` of linear :func:`solve_planar_magnetostatic`
    solves for each source separately; their flux linkages then sum to the nonlinear total.

    ``A`` : converged H1 A_z from :func:`solve_planar_magnetostatic_nonlinear`.
    ``nu_of_B`` : the SAME callable(B_cf) -> nu_cf used for that nonlinear solve.

    Returns ``nu_of_B(curl A)`` as a DIRECT CF (A is fixed, so it is frozen). Two lessons:
    use the direct CF -- do NOT L2-project a clamped/discontinuous nu (projection smears it
    and the freeze stops reproducing the operating point); and the freeze only reproduces a
    CONVERGED nonlinear solve, so a smooth bounded nu(B) (Picard converges) is required --
    re-solving one linear step with the frozen nu and comparing flux linkage is the
    built-in convergence + correctness check. Validated to 0.5 %
    (examples/comsol_class/frozen_permeability.py)."""
    B = CoefficientFunction((grad(A)[1], -grad(A)[0]))
    return nu_of_B(B)


def dq_torque(lambda_m, Ld, Lq, id_, iq, pole_pairs):
    """Synchronous-machine dq electromagnetic torque [N.m] (motor convention):

        T = (3/2) p [ lambda_m * iq + (Ld - Lq) * id * iq ].

    First term = magnet/alignment torque; second = reluctance torque from the
    saliency (Ld != Lq). ``lambda_m`` is the PM flux linkage (e.g. from
    :func:`coil_flux_linkage_2d` on a PM-only solve), Ld/Lq from
    :func:`inductance_2d` along the rotor d-/q-axis (use the frozen-permeability
    values for a saturated machine). ``p`` = pole pairs.

    Force-anchored against an independent salient-rotor FEA: the (Ld,Lq) extracted by
    projecting the FE phase flux linkages through :func:`park_transform` reproduce the
    MEAN Maxwell-stress-tensor torque over a rotor period; the instantaneous stress-tensor
    torque RIPPLES about this value at the pole harmonic (6*theta_e), and the single-instant
    fundamental :func:`dq_flux_torque` matches exactly (tests/test_dq_torque.py)."""
    return 1.5 * pole_pairs * (lambda_m * iq + (Ld - Lq) * id_ * iq)


def dq_torque_components(lambda_m, Ld, Lq, id_, iq, pole_pairs):
    """Split the dq electromagnetic torque into its magnet and reluctance parts [N.m]:

        T_magnet     = (3/2) p  lambda_m iq             (alignment torque from the PM flux)
        T_reluctance = (3/2) p (Ld - Lq) id iq          (saliency torque, 0 if Ld == Lq)
        T_total      = T_magnet + T_reluctance  ==  :func:`dq_torque`.

    The standard interior-/salient-PM design decomposition: how much of a salient PM machine's
    torque comes from the magnet vs from the rotor saliency. At the MTPA operating point
    (:func:`mtpa_operating_point`) both terms are positive and the reluctance share grows with
    the saliency |Ld - Lq| and the current. Returns ``(T_magnet, T_reluctance, T_total)``.

    Force-anchored: on a salient-pole PM FEA the magnet linkage lambda_m (open-circuit Park-d) and
    the saliency Ld/Lq (current injection) reconstruct, through this split, the Maxwell-stress-tensor
    torque whose argmax over the current angle lands on :func:`mtpa_operating_point` (the magnet term
    alone gives the pure-q torque, the reluctance term the angle-dependent gain)."""
    k = 1.5 * pole_pairs
    t_mag = k * lambda_m * iq
    t_rel = k * (Ld - Lq) * id_ * iq
    return t_mag, t_rel, t_mag + t_rel


def mtpa_operating_point(lambda_m, Ld, Lq, current, pole_pairs):
    """MAXIMUM-TORQUE-PER-AMPERE operating point of a salient PM synchronous
    machine at stator current magnitude ``current``.

    Returns ``(gamma, id, iq, torque)`` where ``gamma`` is the current angle from
    the q-axis [rad] (id = -I sin gamma, iq = I cos gamma) and ``torque`` =
    :func:`dq_torque`. MTPA maximises T at fixed |I|; d/dgamma T = 0 gives

        lambda_m sin g + (Ld - Lq) I cos 2g = 0,

    a quadratic in s = sin g:  2(Ld-Lq) I s^2 - lambda_m s - (Ld-Lq) I = 0. The
    T-maximising real root in [-1, 1] is selected, so this is robust to BOTH
    saliency signs (IPM Ld<Lq -> gamma>0; synchronous reluctance lambda_m=0,
    Ld>Lq -> gamma=-45 deg) and to the non-salient limit (Ld=Lq -> gamma=0, pure
    q-axis). Validated against a numerical argmax (tests/test_motor_mtpa.py)."""
    if current == 0:
        return (0.0, 0.0, 0.0, 0.0)
    cands = [0.0]                                   # pure q-axis is always a candidate
    dL = Ld - Lq
    if abs(dL) > 1e-300:
        a, b, c = 2.0 * dL * current, -lambda_m, -dL * current
        disc = b * b - 4.0 * a * c
        if disc >= 0.0:
            root = math.sqrt(disc)
            for s in ((-b + root) / (2.0 * a), (-b - root) / (2.0 * a)):
                if -1.0 <= s <= 1.0:
                    cands.append(math.asin(s))
    best = None
    for g in cands + [-g for g in cands]:
        id_, iq = -current * math.sin(g), current * math.cos(g)
        T = dq_torque(lambda_m, Ld, Lq, id_, iq, pole_pairs)
        if best is None or T > best[3]:
            best = (g, id_, iq, T)
    return best


def base_speed_electrical(lambda_m, Ld, Lq, Imax, Vmax, pole_pairs):
    """Electrical speed [rad/s] at which the MTPA-at-Imax operating point first
    reaches the inverter voltage limit ``Vmax`` (stator R neglected):

        omega_base = Vmax / |(Ld*id_mtpa + lambda_m, Lq*iq_mtpa)|.

    Below it the machine runs at constant torque (MTPA, current-limited); above
    it field weakening begins. (id,iq)_mtpa come from :func:`mtpa_operating_point`
    at |I| = Imax."""
    _, idm, iqm, _ = mtpa_operating_point(lambda_m, Ld, Lq, Imax, pole_pairs)
    return Vmax / math.hypot(Ld * idm + lambda_m, Lq * iqm)


def field_weakening_operating_point(lambda_m, Ld, Lq, Imax, Vmax, omega_e, pole_pairs):
    """Torque-maximising dq operating point of a salient PM synchronous machine at
    electrical speed ``omega_e`` [rad/s], subject to BOTH the inverter limits

        current:  id^2 + iq^2 <= Imax^2
        voltage:  (Ld*id + lambda_m)^2 + (Lq*iq)^2 <= (Vmax/omega_e)^2     (R neglected)

    Returns ``(id, iq, torque, region)`` with ``region`` one of:

      * ``"MTPA"`` — below base speed; the :func:`mtpa_operating_point` at Imax is
        still voltage-feasible (constant-torque region).
      * ``"FW"``   — field weakening on the current circle: the optimum is the
        current-circle ∩ voltage-ellipse corner (constant ~power; torque falls).
      * ``"MTPV"`` — deep field weakening: the max-torque-per-voltage tangent point
        lies inside the current circle (only when characteristic current
        Ich = lambda_m/Ld < Imax). The optimum leaves the current circle.

    Returns ``None`` when no point satisfies both limits — i.e. ``omega_e`` exceeds
    the machine's maximum speed (a finite max speed exists when Ich >= Imax).

    The optimum always lies on the feasible-region boundary, so it is the best of
    three closed-form candidates (MTPA point, current∩voltage corner, MTPV tangent).
    The MTPV point uses the change of variables X = Ld*id + lambda_m, Y = Lq*iq,
    which turns the voltage ellipse into the circle X^2+Y^2=(Vmax/omega_e)^2 and the
    torque into Y*(lambda_m*Lq + (Ld-Lq)*X); maximising on the circle gives a
    quadratic in cos(theta). Continues :func:`mtpa_operating_point` into the
    constant-power / high-speed regions (drive operating-region / efficiency-map
    foundation). Validated to 0.000 % against an independent numeric constrained
    argmax across all three regions (tests/test_field_weakening.py)."""
    Vlim = Vmax / omega_e

    def _volt_ok(id_, iq):
        return (Ld * id_ + lambda_m) ** 2 + (Lq * iq) ** 2 <= Vlim * Vlim * (1 + 1e-6)

    cand = []
    # 1) MTPA at Imax (constant-torque region) — valid only if voltage-feasible
    _, idm, iqm, Tm = mtpa_operating_point(lambda_m, Ld, Lq, Imax, pole_pairs)
    if _volt_ok(idm, iqm):
        cand.append((idm, iqm, Tm, "MTPA"))

    # 2) current circle ∩ voltage ellipse (substitute iq^2 = Imax^2 - id^2):
    #    (Ld^2-Lq^2) id^2 + 2 Ld lambda_m id + (lambda_m^2 + Lq^2 Imax^2 - Vlim^2) = 0
    A = Ld * Ld - Lq * Lq
    B = 2.0 * Ld * lambda_m
    C = lambda_m * lambda_m + Lq * Lq * Imax * Imax - Vlim * Vlim
    ids = []
    if abs(A) < 1e-300:
        if abs(B) > 1e-300:
            ids = [-C / B]
    else:
        disc = B * B - 4.0 * A * C
        if disc >= 0.0:
            r = math.sqrt(disc)
            ids = [(-B + r) / (2.0 * A), (-B - r) / (2.0 * A)]
    for id_ in ids:
        s = Imax * Imax - id_ * id_
        if s < -1e-12:
            continue
        iq = math.sqrt(max(0.0, s))           # motoring branch iq >= 0
        cand.append((id_, iq, dq_torque(lambda_m, Ld, Lq, id_, iq, pole_pairs), "FW"))

    # 3) MTPV tangent point on the voltage ellipse (X=Ld id+lambda_m, Y=Lq iq):
    #    maximise Y*(lambda_m Lq + (Ld-Lq) X) on X^2+Y^2=Vlim^2
    a, b = Ld - Lq, lambda_m * Lq
    mtpv = None
    if abs(a) < 1e-300:                        # non-salient: max at X=0 (id=-lambda_m/Ld)
        mtpv = (-lambda_m / Ld, Vlim / Lq)
    else:
        aa, bb, cc = 2.0 * a * Vlim, b, -a * Vlim   # 2 a Vlim c^2 + b c - a Vlim = 0
        disc = bb * bb - 4.0 * aa * cc
        if disc >= 0.0:
            r = math.sqrt(disc)
            best = None
            for cth in ((-bb + r) / (2.0 * aa), (-bb - r) / (2.0 * aa)):
                if abs(cth) > 1.0:
                    continue
                sth = math.sqrt(max(0.0, 1.0 - cth * cth))
                for Y in (Vlim * sth, -Vlim * sth):
                    idv, iqv = (Vlim * cth - lambda_m) / Ld, Y / Lq
                    T = dq_torque(lambda_m, Ld, Lq, idv, iqv, pole_pairs)
                    if best is None or T > best[2]:
                        best = (idv, iqv, T)
            if best is not None:
                mtpv = (best[0], best[1])
    if mtpv is not None and mtpv[0] ** 2 + mtpv[1] ** 2 <= Imax * Imax * (1 + 1e-6):
        cand.append((mtpv[0], mtpv[1],
                     dq_torque(lambda_m, Ld, Lq, mtpv[0], mtpv[1], pole_pairs), "MTPV"))

    feasible = [c for c in cand if c[1] >= -1e-9
                and c[0] ** 2 + c[1] ** 2 <= Imax * Imax * (1 + 1e-6)
                and _volt_ok(c[0], c[1])]
    return max(feasible, key=lambda c: c[2]) if feasible else None


def induction_machine_thevenin(V1, R1, X1, Xm):
    """Thevenin source the rotor branch sees in a 3-phase induction-machine single-cage
    equivalent circuit (the shunt magnetizing reactance ``Xm`` moved to the source):

        Vth = V1 * Xm / sqrt(R1^2 + (X1+Xm)^2),
        Zth = Rth + j Xth = (j Xm)(R1 + j X1) / (R1 + j(X1+Xm)).

    Returns ``(Vth, Rth, Xth)``. ``V1`` = stator phase voltage, ``R1``/``X1`` = stator
    resistance/leakage reactance, ``Xm`` = magnetizing reactance [all SI]. The reduction
    that makes the torque-slip curve closed-form (:func:`induction_machine_torque`)."""
    Vth = V1 * Xm / math.hypot(R1, X1 + Xm)
    Zth = (complex(0.0, Xm) * complex(R1, X1)) / complex(R1, X1 + Xm)
    return Vth, Zth.real, Zth.imag


def induction_machine_torque(V1, R1, X1, R2, X2, Xm, omega_s, slip):
    """Air-gap electromagnetic torque [N.m] of a 3-phase induction machine at ``slip``
    from the single-cage equivalent circuit (Thevenin form):

        T(s) = (3/omega_s) * Vth^2 * (R2/s) / [ (Rth + R2/s)^2 + (Xth + X2)^2 ],

    ``omega_s`` = SYNCHRONOUS MECHANICAL speed [rad/s] (= 2 pi f / pole_pairs), ``R2``/``X2``
    = rotor resistance/leakage referred to the stator. T>0 motoring (0<s<1), peaks at the
    breakdown slip (:func:`induction_machine_breakdown`), -> 0 as s -> 0 (synchronism) and
    is the starting torque at s=1. The induction-machine analog of the PM dq torque
    :func:`dq_torque`. Validated against a numeric slip sweep (tests/test_induction_machine.py)."""
    Vth, Rth, Xth = induction_machine_thevenin(V1, R1, X1, Xm)
    r2s = R2 / slip
    return (3.0 / omega_s) * Vth * Vth * r2s / ((Rth + r2s) ** 2 + (Xth + X2) ** 2)


def induction_machine_breakdown(V1, R1, X1, R2, X2, Xm, omega_s):
    """Breakdown (pull-out) point of an induction machine -- the maximum of
    :func:`induction_machine_torque` over slip:

        s_max = R2 / sqrt(Rth^2 + (Xth + X2)^2),
        T_max = 3 Vth^2 / ( 2 omega_s [ Rth + sqrt(Rth^2 + (Xth + X2)^2) ] ).

    Returns ``(s_max, T_max)``. The classic INVARIANT: T_max is INDEPENDENT of the rotor
    resistance R2 -- only the breakdown SLIP scales with R2 (s_max ∝ R2), which is how
    rotor resistance (wound-rotor / deep-bar) shifts peak torque toward standstill for
    starting without changing its height. Closed form == numeric argmax (validated)."""
    Vth, Rth, Xth = induction_machine_thevenin(V1, R1, X1, Xm)
    root = math.hypot(Rth, Xth + X2)
    return R2 / root, 3.0 * Vth * Vth / (2.0 * omega_s * (Rth + root))


def synchronous_power_angle_torque(V, E, Xd, Xq, omega_e, pole_pairs, delta):
    """Electromagnetic torque [N.m] of a SALIENT synchronous machine vs the LOAD ANGLE
    ``delta`` (the angle between the terminal voltage V and the excitation EMF E) -- the
    grid/voltage-frame characteristic (vs the current-frame :func:`dq_torque`/MTPA):

        T(d) = (3 p/omega_e) [ V E/Xd sin d  +  (V^2/2)(1/Xq - 1/Xd) sin 2d ].

    First term = FIELD/excitation torque (~ sin d); second = RELUCTANCE torque (~ sin 2d) from
    the saliency Xd != Xq. ``V``/``E`` per-phase RMS, ``Xd``/``Xq`` synchronous reactances,
    ``omega_e`` electrical speed [rad/s], ``p`` pole pairs. For a PM machine E = omega_e*lambda_m.
    Generating for d>0, motoring for d<0 (mirror)."""
    field = V * E / Xd * math.sin(delta)
    rel = 0.5 * V * V * (1.0 / Xq - 1.0 / Xd) * math.sin(2.0 * delta)
    return (3.0 * pole_pairs / omega_e) * (field + rel)


def synchronous_pullout(V, E, Xd, Xq, omega_e, pole_pairs):
    """Pull-out (steady-state stability limit) of a salient synchronous machine -- the maximum
    of :func:`synchronous_power_angle_torque` over the load angle. dT/dd=0 gives
    ``(V E/Xd) cos d + V^2(1/Xq-1/Xd) cos 2d = 0``, a quadratic in cos d (cos 2d = 2cos^2 d - 1).

    Returns ``(delta_max, T_max)`` [rad, N.m]. LIMITS: NON-salient (Xd=Xq) -> delta_max = 90 deg
    (classic round-rotor pull-out at the quadrature point); RELUCTANCE-only (E=0) -> delta_max =
    45 deg; salient field machine -> 45 < delta_max < 90 deg (the reluctance term advances the
    peak). Closed form == numeric argmax (validated)."""
    sal = V * V * (1.0 / Xq - 1.0 / Xd)
    B = V * E / Xd
    if abs(sal) < 1e-300:                       # non-salient: cos d = 0 -> 90 deg
        d = math.pi / 2.0
        return d, synchronous_power_angle_torque(V, E, Xd, Xq, omega_e, pole_pairs, d)
    A, C = 2.0 * sal, -sal
    disc = B * B - 4.0 * A * C
    best = None
    if disc >= 0.0:
        root = math.sqrt(disc)
        for c in ((-B + root) / (2.0 * A), (-B - root) / (2.0 * A)):
            if -1.0 <= c <= 1.0:
                d = math.acos(c)
                T = synchronous_power_angle_torque(V, E, Xd, Xq, omega_e, pole_pairs, d)
                if best is None or T > best[1]:
                    best = (d, T)
    return best


def dq_voltages(R, Ld, Lq, lambda_m, id_, iq, omega_e):
    """Steady-state dq terminal voltages of a PM synchronous machine (per-phase, motor
    convention):

        v_d = R*id - omega_e*lambda_q ,   v_q = R*iq + omega_e*lambda_d ,
        lambda_d = Ld*id + lambda_m ,     lambda_q = Lq*iq.

    ``omega_e`` electrical speed [rad/s], ``R`` per-phase resistance. Returns ``(vd, vq)``;
    the terminal voltage magnitude is hypot(vd, vq). The voltage side of the dq model that the
    field-weakening limit (:func:`field_weakening_operating_point`) constrains."""
    lam_d = Ld * id_ + lambda_m
    lam_q = Lq * iq
    return R * id_ - omega_e * lam_q, R * iq + omega_e * lam_d


def dq_operating_point(R, Ld, Lq, lambda_m, id_, iq, omega_e, pole_pairs):
    """Full dq steady-state OPERATING POINT of a PM synchronous machine at (id, iq, omega_e):
    the terminal voltage, currents, power split, power factor and efficiency. Where MTPA
    (:func:`mtpa_operating_point`) / field weakening choose the currents, this EVALUATES the
    terminal quantities the inverter sees there. Returns a dict:

        vd, vq, Vmag, Imag, torque, P_in, P_em, P_cu, power_factor, efficiency, omega_mech.

    P_in = (3/2)(vd id + vq iq) (electrical input), P_em = (3/2)omega_e(lambda_d iq -
    lambda_q id) = torque*omega_mech (air-gap power), P_cu = (3/2)R|I|^2 (copper loss);
    P_in = P_cu + P_em exactly, and the torque matches :func:`dq_torque`. power_factor =
    P_in/((3/2)|V||I|); efficiency = P_em/P_in (motoring; copper loss only). The terminal
    (V, I, P, PF) layer of the dq machine model -- inverter/drive sizing."""
    vd, vq = dq_voltages(R, Ld, Lq, lambda_m, id_, iq, omega_e)
    lam_d, lam_q = Ld * id_ + lambda_m, Lq * iq
    Vmag, Imag = math.hypot(vd, vq), math.hypot(id_, iq)
    P_in = 1.5 * (vd * id_ + vq * iq)
    P_em = 1.5 * omega_e * (lam_d * iq - lam_q * id_)
    P_cu = 1.5 * R * Imag * Imag
    T = dq_torque(lambda_m, Ld, Lq, id_, iq, pole_pairs)
    denom = 1.5 * Vmag * Imag
    return {
        "vd": vd, "vq": vq, "Vmag": Vmag, "Imag": Imag, "torque": T,
        "P_in": P_in, "P_em": P_em, "P_cu": P_cu,
        "power_factor": (P_in / denom) if denom > 0 else 0.0,
        "efficiency": (P_em / P_in) if P_in > 0 else float("nan"),
        "omega_mech": omega_e / pole_pairs,
    }


def characteristic_current(lambda_m, Ld):
    """Characteristic current of a PM machine, Ich = lambda_m/Ld [A]. The HIGH-SPEED steady
    3-phase short-circuit current, and the d-axis current that exactly cancels the magnet flux
    (lambda_d = 0). The wide-CPSR / 'infinite max speed' design target is Ich = Imax
    (:func:`field_weakening_operating_point`); Ich < Imax enables the MTPV region. Also the
    worst-case DEMAGNETISING current to check against the magnet knee."""
    return lambda_m / Ld


def short_circuit_dq_currents(R, Ld, Lq, lambda_m, omega_e):
    """Steady-state dq currents of a PM synchronous machine with the 3-phase terminals SHORTED
    (v_d = v_q = 0) at electrical speed ``omega_e``. The dq voltage equations give

        0 = R id - omega_e Lq iq ,   0 = R iq + omega_e (Ld id + lambda_m)
        =>  id = -omega_e^2 Lq lambda_m /(R^2 + omega_e^2 Ld Lq),
            iq = -omega_e R   lambda_m /(R^2 + omega_e^2 Ld Lq).

    Returns ``(id, iq)``. HIGH-SPEED limit (omega_e -> inf, or R -> 0): id -> -lambda_m/Ld =
    -:func:`characteristic_current` (pure d-axis DEMAGNETISING), iq -> 0, so |Isc| -> Ich.
    The BRAKING torque is :func:`dq_torque` of these currents -- it rises from 0, peaks near the
    critical speed (omega_e = R/Ls for a non-salient machine, Ld=Lq=Ls), then -> 0 at high speed.
    The fault / demagnetisation-risk regime (vs the operating-point :func:`dq_operating_point`)."""
    den = R * R + omega_e * omega_e * Ld * Lq
    return (-omega_e * omega_e * Lq * lambda_m / den,
            -omega_e * R * lambda_m / den)


def winding_distribution_factor(harmonic, slots_per_pole_per_phase, slot_angle_elec_deg):
    """DISTRIBUTION factor k_d,n of a distributed AC winding -- the reduction of the n-th
    harmonic EMF because the ``q = slots_per_pole_per_phase`` coils of a phase belt are spread
    over consecutive slots (electrical angle ``slot_angle_elec_deg`` = gamma apart):

        k_d,n = sin(n q gamma/2) / (q sin(n gamma/2))   = |sum_{i<q} exp(j n i gamma)| / q.

    q=1 (concentrated) -> 1. The phasor-sum form is the geometric series of the q slot EMFs."""
    q = slots_per_pole_per_phase
    g = math.radians(slot_angle_elec_deg)
    den = q * math.sin(harmonic * g / 2.0)
    if abs(den) < 1e-14:                       # n*gamma a multiple of 2pi (slot harmonic) -> |kd|=1
        return 1.0
    return math.sin(harmonic * q * g / 2.0) / den


def winding_pitch_factor(harmonic, pitch_fraction):
    """PITCH (chording) factor k_p,n = sin(n * beta * pi/2), ``beta = pitch_fraction`` = coil span
    / pole pitch (<=1). beta=1 (full pitch) -> k_p,1 = 1. SHORT-PITCHING kills harmonics: a 5/6
    pitch gives |k_p,5| = |k_p,7| = sin(5pi/12) = 0.259, suppressing the 5th and 7th."""
    return math.sin(harmonic * pitch_fraction * math.pi / 2.0)


def winding_factor(harmonic, slots_per_pole_per_phase, slot_angle_elec_deg, pitch_fraction=1.0):
    """AC-winding factor k_w,n = k_d,n * k_p,n (:func:`winding_distribution_factor` *
    :func:`winding_pitch_factor`) -- the n-th harmonic EMF is k_w,n times the
    concentrated-full-pitch value. The fundamental k_w,1 (~0.92-0.96 for a good 3-phase winding)
    scales the back-EMF / torque constant (#34, F3); the harmonic k_w,n (small for a
    distributed, short-pitched winding) set the EMF harmonic content and the parasitic torques.
    Validated against the direct phasor sum (tests/test_winding_factor.py)."""
    return (winding_distribution_factor(harmonic, slots_per_pole_per_phase, slot_angle_elec_deg)
            * winding_pitch_factor(harmonic, pitch_fraction))


def single_phase_mmf_harmonic(harmonic, winding_factor_n, turns_per_phase, current, pole_pairs):
    """Peak amplitude [A] of the n-th SPACE harmonic of one phase's air-gap MMF:

        F_phi,n = (4/pi) (k_w,n / n) (N_ph I / 2p),

    the (4/pi)/n square-wave envelope of a full-pitch coil, reduced by the winding factor k_w,n
    (:func:`winding_factor`) of the distributed / short-pitched belt. ``turns_per_phase`` = N_ph
    series turns, ``current`` = peak phase current, ``pole_pairs`` = p. The fundamental
    F_phi,1 = (2/pi) k_w,1 N_ph I / p is the per-phase MMF that drives the magnetizing inductance
    (#69) and the air-gap field; the harmonics drive parasitic torques and rotor losses."""
    return (4.0 / math.pi) * (winding_factor_n / harmonic) * (turns_per_phase * current / (2.0 * pole_pairs))


def rotating_mmf_amplitude(harmonic, phases, single_phase_amplitude):
    """Amplitude of the ROTATING air-gap MMF wave that the n-th space harmonic forms in a balanced
    ``phases``-phase winding: the m phases (spaced 2 pi/m in both space and time) combine to
    (m/2) F_phi,n for the harmonics that form a rotating wave (:func:`mmf_harmonic_direction` != 0)
    and to exactly ZERO for the triplen (n a multiple of m -- e.g. the 3rd / 9th cancel in a 3-phase
    machine). So the rotating fundamental is (3/2) F_phi,1 (:func:`single_phase_mmf_harmonic`); this
    is the (m/2) factor behind the synchronous magnetizing inductance L_m = (m/2) L_mu (#69)."""
    return 0.0 if mmf_harmonic_direction(harmonic, phases) == 0 else 0.5 * phases * single_phase_amplitude


def mmf_harmonic_direction(harmonic, phases):
    """Rotation sense of the n-th MMF space harmonic in a balanced ``phases``-phase winding (the
    rotating-field theorem): +1 FORWARD (with the fundamental) if (n-1) is a multiple of m, -1
    BACKWARD if (n+1) is, and 0 if the harmonic CANCELS (n a multiple of m, or -- for m>=5 -- any n
    that is neither +-1 mod m). For the 3-phase machine the ODD harmonics a symmetric winding
    actually carries are n = 6k +/- 1: forward for 6k+1 (1, 7, 13, ...), backward for 6k-1 (5, 11,
    ...), while the triplen (3, 9, ...) cancel (even harmonics are absent by half-wave symmetry).
    The backward 5th / forward 7th family drive the dominant parasitic asynchronous torques and the
    rotor-surface / magnet-eddy losses."""
    if harmonic % phases == 0:
        return 0
    if (harmonic - 1) % phases == 0:
        return +1
    if (harmonic + 1) % phases == 0:
        return -1
    return 0


def skew_factor(harmonic, skew_angle_elec):
    """Skew factor k_sk,n of a winding/magnet SKEWED continuously over an electrical angle
    ``skew_angle_elec`` (theta_sk):

        k_sk,n = sin(n theta_sk/2)/(n theta_sk/2)        (the sinc / averaging factor).

    Skewing averages the n-th harmonic of the air-gap flux over the skew, so the n-th EMF and the
    n-th torque harmonic are scaled by k_sk,n. The axial twin of the distribution factor (#51,
    which spreads a phase over slots); here the conductor is spread along the STACK. A harmonic
    whose period equals the skew (n theta_sk = 2 pi) is NULLED -- skewing by ONE SLOT PITCH
    (:func:`slot_pitch_skew_angle`) cancels the slot-passing / cogging ripple. The fundamental
    pays a small cost k_sk,1 < 1. Validated against the phasor average (tests/test_skew_factor.py)."""
    x = harmonic * skew_angle_elec / 2.0
    return 1.0 if abs(x) < 1e-12 else math.sin(x) / x


def slot_pitch_skew_angle(slots, pole_pairs):
    """Electrical angle of ONE SLOT PITCH = 2 pi * pole_pairs / slots [rad] -- the usual skew for
    cogging cancellation (the cogging/slot-passing harmonic then sees k_sk = 0,
    :func:`skew_factor`). ``slots`` = stator slots Q, ``pole_pairs`` = p."""
    return 2.0 * math.pi * pole_pairs / slots


def skewed_winding_factor(harmonic, slots_per_pole_per_phase, slot_angle_elec_deg,
                          skew_angle_elec, pitch_fraction=1.0):
    """Total n-th-harmonic winding factor WITH skew: k_w,n * k_sk,n = (:func:`winding_factor`) *
    (:func:`skew_factor`). The complete EMF-scaling of a distributed, chorded AND skewed winding
    -- distribution + pitch suppress harmonics in the slot plane, skew suppresses them along the
    stack (and cancels cogging), at a small fundamental cost."""
    return (winding_factor(harmonic, slots_per_pole_per_phase, slot_angle_elec_deg, pitch_fraction)
            * skew_factor(harmonic, skew_angle_elec))


def skew_average(theta, values, skew_span, n_slices=0):
    """Multi-slice SKEW of a rotor-angle sweep -- the FE-multislice 2D->3D skew
    correction (the ONELAB ``generate_3d.geo`` ``Model3D = Multi-slice`` technique).

    Given a quantity ``values`` (torque, back-EMF, cogging, ...) sampled at UNIFORM,
    PERIODIC rotor angles ``theta`` over exactly one full period, return the curve a
    machine SKEWED by ``skew_span`` (same angle units as ``theta``) would produce --
    the axial average of the rotor field over the skew:

        f_skew(t) = (1/skew_span) int_0^{skew_span} f(t + s) ds          (n_slices=0, exact)
        f_skew(t) = (1/N) sum_{k=0}^{N-1} f(t + (k+0.5) skew_span/N)     (N discrete slices)

    A harmonic of order n is thereby scaled by the sinc factor
    sin(n*skew_span/2)/(n*skew_span/2) = :func:`skew_factor` (EXACT in the continuous
    n_slices=0 mode), and the slot-passing / cogging harmonic is NULLED at one
    slot-pitch skew (:func:`slot_pitch_skew_angle`).  Unlike the analytic
    ``skew_factor``, this operates on the ACTUAL (possibly saturated, non-sinusoidal)
    swept curve -- the FE counterpart that captures every harmonic at once.

    ``values`` must cover EXACTLY one period (so the periodic continuation is valid).
    n_slices=0 uses the spectral (FFT) form (exact continuous skew); n_slices=N
    mimics N physical slices and converges to it as N grows.  Returns an array like
    ``values``.  numpy-only (no FE solve).
    """
    import numpy as np
    v = np.asarray(values, dtype=float)
    th = np.asarray(theta, dtype=float)
    Ns = v.size
    dth = (th[-1] - th[0]) / (Ns - 1)
    if n_slices and n_slices > 0:
        out = np.zeros_like(v)
        idx0 = (th - th[0]) / dth                       # = [0, 1, ..., Ns-1]
        for k in range(int(n_slices)):
            shift = (k + 0.5) * skew_span / n_slices    # midpoint of slice k
            idx = (idx0 + shift / dth) % Ns
            lo = np.floor(idx).astype(int); fr = idx - lo
            out += (1.0 - fr) * v[lo % Ns] + fr * v[(lo + 1) % Ns]
        return out / n_slices
    # spectral (FFT-exact) continuous skew: each mode k -> e^{ik s/2} sinc(k s/2)
    C = np.fft.fft(v)
    k = 2.0 * np.pi * np.fft.fftfreq(Ns, d=dth)          # angular frequency per bin
    x = k * skew_span / 2.0
    sinc = np.where(np.abs(x) < 1e-12, 1.0, np.sin(x) / np.where(x == 0.0, 1.0, x))
    return np.real(np.fft.ifft(C * (np.exp(1j * x) * sinc)))


def slot_leakage_permeance(conductor_height, slot_width, opening_height=0.0, opening_width=None):
    """Specific (per-unit-axial-length) SLOT-LEAKAGE permeance of a rectangular machine slot:

        lambda_s = h_c/(3 w_s)  +  h_o/w_o   [dimensionless],

    the leakage flux that crosses the slot from tooth to tooth and links the conductor without
    reaching the air gap. ``h_c`` = conductor height, ``w_s`` = slot (conductor) width. The
    **1/3 factor** is exact: a uniform-current conductor builds the cross-slot MMF LINEARLY from
    0 at the slot bottom to NI at the conductor top, and integrating B^2 over that triangular
    field gives 1/3 (vs 1 for a current-free region of the same shape). The optional second term
    is the tooth-tip OPENING above the conductor (height ``h_o``, width ``w_o`` <= w_s) -- a
    current-free series air permeance h_o/w_o (defaults to the slot width if not given).

    Use ``slot_leakage_inductance`` for the inductance. The 1/3 conductor term and the
    same-width opening term are EXACT (validated to 0.000 %, examples/comsol_class/
    slot_leakage.py); a NARROWED opening (w_o < w_s) adds a few % of tooth-shoulder fringing the
    simple series form omits, so the closed form is then a slight LOWER bound."""
    lam = conductor_height / (3.0 * slot_width)
    if opening_height > 0.0:
        lam += opening_height / (opening_width if opening_width else slot_width)
    return lam


def slot_leakage_inductance(conductor_height, slot_width, axial_length, turns=1,
                            opening_height=0.0, opening_width=None):
    """Slot-leakage INDUCTANCE [H] of one slot: L = mu0 * N^2 * l_stk * lambda_s, with
    ``lambda_s`` from :func:`slot_leakage_permeance` (``turns`` N, ``axial_length`` l_stk).

    The leakage-reactance contribution X_sigma = omega L that every machine winding carries on
    top of the magnetising reactance -- it limits the short-circuit current and shapes the
    transient response, and (unlike the air-gap inductance) is set by the slot GEOMETRY, not the
    iron. For a single bulk conductor (N=1) per unit length this is exactly mu0 * lambda_s [H/m]."""
    return MU0 * turns * turns * axial_length * slot_leakage_permeance(
        conductor_height, slot_width, opening_height, opening_width)


def deep_bar_skin_ratio(height, f, sigma, mu_r=1.0):
    """Deep-bar skin parameter xi = h/delta of a solid conductor of height ``height`` in a slot,
    delta = sqrt(2/(omega mu sigma)) the skin depth at frequency ``f`` [Hz] (mu = mu_r*mu0):

        xi = height * sqrt(pi f mu0 mu_r sigma).

    The leakage field in an iron slot is purely transverse, so the AC current redistribution is a
    1-D skin effect over the slot DEPTH -- xi sets how deep the deep-bar effect runs (xi<<1 DC,
    xi>>1 strongly skin-limited)."""
    return height * math.sqrt(math.pi * f * MU0 * mu_r * sigma)


def deep_bar_resistance_factor(xi):
    """Deep-bar AC/DC RESISTANCE factor k_R = Rac/Rdc of a solid conductor filling an iron slot
    (Field's formula):

        k_R = xi (sinh 2xi + sin 2xi)/(cosh 2xi - cos 2xi),   xi = h/delta = deep_bar_skin_ratio.

    The slot-leakage field crowds the AC current toward the slot opening, raising the resistance.
    xi->0: k_R->1 (DC); xi>>1: k_R->xi (Rac ~ sqrt(f), current in a one-skin-depth surface layer).
    The DEEP-BAR / double-cage rotor effect that gives induction motors their high STARTING
    resistance (#40), and the AC copper loss of a conductor in a slot at high frequency."""
    if xi < 1e-2:
        return 1.0
    den = math.cosh(2 * xi) - math.cos(2 * xi)
    return xi * (math.sinh(2 * xi) + math.sin(2 * xi)) / den


def deep_bar_reactance_factor(xi):
    """Deep-bar AC/DC slot-leakage REACTANCE (inductance) factor k_X = Lac/Ldc of a solid
    conductor filling an iron slot:

        k_X = (3/2xi)(sinh 2xi - sin 2xi)/(cosh 2xi - cos 2xi),   xi = deep_bar_skin_ratio.

    The same skin redistribution REDUCES the slot-leakage inductance (the lower conductor no
    longer links the full slot flux). xi->0: k_X->1 (the DC slot leakage mu0 h/3w,
    :func:`slot_leakage_permeance`); xi>>1: k_X->3/(2 xi)->0. The reactance companion of
    :func:`deep_bar_resistance_factor` -- together they are the deep-bar circuit Z(slip)."""
    if xi < 1e-2:
        return 1.0
    den = math.cosh(2 * xi) - math.cos(2 * xi)
    return (3.0 / (2.0 * xi)) * (math.sinh(2 * xi) - math.sin(2 * xi)) / den


def dowell_resistance_factor(skin_ratio, layers):
    """Dowell AC/DC resistance factor F_R of an ``layers``-layer (m-layer) winding -- the
    transformer / inductor winding-loss workhorse:

        F_R = Delta [ (sinh2D + sin2D)/(cosh2D - cos2D)
                       + (2/3)(m^2 - 1)(sinhD - sinD)/(coshD + cosD) ],   D = Delta = h/delta,

    ``skin_ratio`` = Delta = foil thickness / skin depth (:func:`deep_bar_skin_ratio`), ``layers`` = m.
    The first (skin) term is the single-conductor deep-bar factor (:func:`deep_bar_resistance_factor`,
    m=1); the second is the PROXIMITY loss, which grows as ~m^2 -- the field from the other layers
    builds up across the stack, so the outer layers see a large field and dissipate heavily. The
    reason multi-layer (and litz vs solid) winding design is dominated by AC loss: halving the foil
    thickness or splitting layers cuts F_R fast. Validated against the m-foil eddy FE to 0.00 %
    (examples/comsol_class/dowell_winding.py)."""
    if skin_ratio < 1e-2:
        return 1.0
    D = skin_ratio
    skin = (math.sinh(2 * D) + math.sin(2 * D)) / (math.cosh(2 * D) - math.cos(2 * D))
    prox = (math.sinh(D) - math.sin(D)) / (math.cosh(D) + math.cos(D))
    return D * (skin + (2.0 / 3.0) * (layers * layers - 1) * prox)


def classical_eddy_loss_density(sigma, thickness, frequency, B_peak):
    """Classical (thin-lamination) eddy-current loss VOLUME density [W/m^3] of a steel sheet of
    thickness ``d`` carrying a sinusoidal flux density of peak ``B_peak`` at ``frequency``:

        P_eddy = pi^2 sigma d^2 f^2 B_peak^2 / 6   (= sigma omega^2 d^2 B_peak^2 / 24).

    The d^2 law is WHY cores are laminated: halving the sheet thickness quarters the eddy loss. It
    is the f^2 ('eddy') term of the Steinmetz core-loss separation P = P_hyst(f) + P_eddy(f^2). Valid
    in the THIN limit d << skin depth (uniform flux across the sheet); the skin-effect roll-off at
    higher f / thicker sheets is :func:`lamination_eddy_skin_factor`. ``sigma`` [S/m], ``thickness``
    [m], ``frequency`` [Hz], ``B_peak`` [T]."""
    return math.pi ** 2 * sigma * thickness ** 2 * frequency ** 2 * B_peak ** 2 / 6.0


def lamination_eddy_skin_factor(thickness, skin_depth):
    """Skin-effect REDUCTION factor F(xi) (<= 1) of the lamination eddy loss vs the classical d^2
    value (:func:`classical_eddy_loss_density`), from the exact 1-D magnetic-diffusion solution of a
    slab in a fixed average AC flux:

        F(xi) = (3/xi) (sinh xi - sin xi)/(cosh xi - cos xi),    xi = d/delta,

    delta = sqrt(2/(omega mu sigma)) the skin depth. F -> 1 as xi -> 0 (thin sheet, the classical
    law) and F -> 3/xi as xi -> inf (heavy skin effect: the flux cannot penetrate, so the loss grows
    only as d, not d^2, and the classical formula over-predicts). The field-driven counterpart of the
    current-driven deep-bar k_R (:func:`deep_bar_resistance_factor`) -- same diffusion, different
    drive. Validated against the numerically-integrated exact slab loss to 1e-6
    (examples/comsol_class/lamination_eddy_loss.py)."""
    xi = thickness / skin_depth
    if xi < 0.1:
        return 1.0 - xi ** 4 / 630.0              # series (the direct form loses precision near 0)
    if xi > 30.0:
        return 3.0 / xi                            # cosh/sinh overflow guard; F -> 3/xi
    return (3.0 / xi) * (math.sinh(xi) - math.sin(xi)) / (math.cosh(xi) - math.cos(xi))


def lamination_eddy_loss_density(sigma, thickness, frequency, B_peak, mu_r=1.0):
    """Skin-corrected lamination eddy-loss VOLUME density [W/m^3]:

        P = (pi^2 sigma d^2 f^2 B_peak^2 / 6) * F(d/delta),

    the classical loss :func:`classical_eddy_loss_density` times the reduction factor
    :func:`lamination_eddy_skin_factor` with delta = sqrt(2/(omega mu0 mu_r sigma)). ``mu_r`` is the
    sheet's relative permeability (large for steel -> thinner skin depth -> bigger xi). In the thin
    limit (typical 0.35 mm sheet at 50/60 Hz) F ~ 1 and this equals the classical value; it rolls off
    for thick sheets or high frequency (PWM harmonics)."""
    omega = 2.0 * math.pi * frequency
    delta = math.sqrt(2.0 / (omega * MU0 * mu_r * sigma))
    return classical_eddy_loss_density(sigma, thickness, frequency, B_peak) \
        * lamination_eddy_skin_factor(thickness, delta)


def carter_coefficient(slot_pitch, gap, slot_opening):
    """Carter coefficient k_C of a slotted air gap -- the factor by which the slot openings make
    the gap behave MAGNETICALLY LARGER:

        gamma = (b_o/g)^2/(5 + b_o/g),   k_C = tau_s/(tau_s - gamma * g)   (>= 1),

    tau_s = slot pitch, g = gap length, b_o = slot opening. A slot opening adds reluctance (the
    flux spreads as it crosses the opening), so the AVERAGE gap flux for a given MMF drops as if
    the gap were k_C*g (:func:`effective_air_gap`). The standard machine-design correction that
    feeds the magnetising inductance, the no-load flux, and the slot-ripple permeance. Carter's
    1901 conformal-mapping result; the (b_o/g)^2/(5+b_o/g) form is the textbook fit (validated
    against the scalar-potential gap permeance to < 0.3 %, examples/comsol_class/carter_coefficient.py)."""
    gamma = (slot_opening / gap) ** 2 / (5.0 + slot_opening / gap)
    return slot_pitch / (slot_pitch - gamma * gap)


def effective_air_gap(slot_pitch, gap, slot_opening):
    """Effective (Carter) air gap g_eff = k_C * g of a slotted machine -- the SMOOTH gap that has
    the same magnetising permeance as the real slotted gap (:func:`carter_coefficient`). Use it in
    place of the physical gap in the air-gap permeance / magnetising-inductance / no-load-flux
    calculations; the slot openings effectively lengthen the gap by a few to ~30 %."""
    return carter_coefficient(slot_pitch, gap, slot_opening) * gap


def magnetizing_inductance_per_phase(gap_diameter, stack_length, effective_gap, pole_pairs,
                                     kw1, turns_per_phase):
    """Per-phase (self) MAGNETISING inductance L_mu of an AC machine -- the air-gap main-flux
    inductance, the dominant part of the dq L_md / L_mq (the remainder is leakage, e.g. slot
    leakage #54):

        L_mu = (2 mu0/pi) (kw1 N_ph)^2 D L_stk / (p^2 g_eff),

    D = air-gap diameter [m], L_stk = stack length [m], g_eff = effective (Carter #57) gap [m],
    p = pole PAIRS, kw1 = fundamental winding factor (#51), N_ph = series turns per phase. Follows
    from the fundamental MMF -> air-gap B -> flux/pole -> flux-linkage chain. SYNTHESISES the
    campaign -- kw1 from :func:`winding_factor`, g_eff from :func:`effective_air_gap` -- and scales
    as 1/g_eff (Carter slotting lowers it), 1/p^2, and (kw1 N_ph)^2. Validated against the explicit
    derivation chain (examples/comsol_class/magnetizing_inductance.py)."""
    return (2.0 * MU0 / math.pi) * (kw1 * turns_per_phase) ** 2 * gap_diameter * stack_length \
        / (pole_pairs ** 2 * effective_gap)


def synchronous_magnetizing_inductance(gap_diameter, stack_length, effective_gap, pole_pairs,
                                       kw1, turns_per_phase, phases=3):
    """Synchronous (per-phase) MAGNETISING inductance L_m = (m/2) L_mu of an m-phase machine -- the
    per-phase self :func:`magnetizing_inductance_per_phase` raised by m/2 (3/2 for three phases) by
    the mutual coupling of the m phases sharing the common rotating air-gap field. This is the L_md
    (= L_mq for a non-salient machine) of the dq model; the terminal Ld = L_m + L_leakage (slot
    leakage #54 + end-winding + harmonic leakage)."""
    return 0.5 * phases * magnetizing_inductance_per_phase(
        gap_diameter, stack_length, effective_gap, pole_pairs, kw1, turns_per_phase)


# Convention: the FIRST component is the MAGNETIZED-axis factor; all three sum to 1.
_DEMAG_FACTORS = {
    "sphere": (1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0),
    "cylinder_transverse": (0.5, 0.5, 0.0),         # magnetized across an infinite cylinder: N=1/2
    "cylinder_axial": (0.0, 0.5, 0.5),              # magnetized along the infinite axis:   N=0
    "thin_slab": (1.0, 0.0, 0.0),                   # magnetized through the thin dimension: N=1
}


def demagnetizing_factor(shape):
    """Demagnetizing factor tuple (N_mag, N_b, N_c) of a uniformly magnetized body -- the geometry
    factor that sets the self-demagnetizing field H = -N M inside it. The FIRST component is the
    factor along the MAGNETIZED axis; for ANY ellipsoid the three sum to one (Nx + Ny + Nz = 1).
    Exact closed forms for the limiting shapes:

        sphere (1/3, 1/3, 1/3); infinite cylinder transverse (1/2, 1/2, 0) / axial (0, 1/2, 1/2);
        thin slab through-thickness (1, 0, 0).

    ``shape`` in {'sphere','cylinder_transverse','cylinder_axial','thin_slab'}. The basis of the PM
    operating point / load line (#42) and of shape-anisotropy design (an axially-magnetized rod
    keeps its remanence; a transversely-magnetized disc loses half)."""
    try:
        return _DEMAG_FACTORS[shape]
    except KeyError:
        raise ValueError("shape must be one of %s" % sorted(_DEMAG_FACTORS))


def magnetized_body_internal_field(remanence, demag_factor):
    """Internal flux density B_in = Br (1 - N) [T] of a uniformly magnetized body (remanence
    Br = mu0 M) along its magnetized axis with demagnetizing factor ``demag_factor`` = N: a sphere
    (N = 1/3) holds 2 Br/3, a transverse-magnetized infinite cylinder (N = 1/2) holds Br/2, a thin
    slab magnetized through-thickness (N = 1) holds ~0. The accompanying internal field is the
    self-demagnetizing :func:`demagnetizing_field`. Validated against an NGSolve transverse-cylinder
    solve to ~0.2 % (open-domain truncation; examples/comsol_class/demagnetizing_factor.py)."""
    return remanence * (1.0 - demag_factor)


def demagnetizing_field(remanence, demag_factor):
    """Internal (self-)demagnetizing field H_in = -N Br/mu0 [A/m] of a uniformly magnetized body
    (remanence Br, demag factor N) -- the field the magnet applies to itself; the larger N (the
    'fatter' the magnet across its magnetization), the deeper into the second quadrant the operating
    point sits, i.e. the closer to the demag knee (#42). Negative (opposes M)."""
    return -demag_factor * remanence / MU0


def solve_planar_magnetostatic(mesh, nu, Jz=None, magnets=None, order=2,
                               dirichlet="outer", magnet_cf=None):
    """2D PLANAR (Cartesian) A_z magnetostatics -- the FEMM ``prob1big`` analog
    on standard NGSolve H1 (no axihenrotte needed; planar has no 1/r singularity).

    Solves  -div(nu grad A_z) = Jz  with optional permanent magnets. The flux
    density is ``B = (grad(A_z)[1], -grad(A_z)[0])`` = (dA/dy, -dA/dx).

    Parameters
    ----------
    mesh      : 2D (x, y) mesh; outer boundary named ``dirichlet`` (A_z = 0).
    nu        : reluctivity CF (e.g. from :func:`reluctivity`).
    Jz        : out-of-plane current density CF [A/m^2] or None.
    magnets   : {material: (Hc, phi_deg)} for permanent magnets, or None.
                See :func:`planar_magnet_source`.
    order     : H1 order (default 2).

    Validated: transverse-magnetized linear cylinder (mu_r=2, Hc=3e5) ->
    interior B = mu0*mu_r*Hc/(mu_r+1) within 0.05 %; external 2D dipole
    B_in*(a/r)^2 within ~1 % (see tests/test_planar_magnet.py).

    Returns the H1 GridFunction ``gfu`` (A_z); flux density ``B = CF((grad(gfu)[1], -grad(gfu)[0]))``.
    """
    fes = H1(mesh, order=order, dirichlet=dirichlet)
    u, v = fes.TnT()
    a = BilinearForm(fes, symmetric=True)
    a += nu * grad(u) * grad(v) * dx
    f = LinearForm(fes)
    if Jz is not None:
        f += Jz * v * dx
    if magnets is not None:
        s = planar_magnet_source(mesh, magnets)
        f += (s[0] * grad(v)[1] - s[1] * grad(v)[0]) * dx
    if magnet_cf is not None:
        # spatially-varying in-plane magnetization s = nu0*Br (e.g. a Halbach array)
        f += (magnet_cf[0] * grad(v)[1] - magnet_cf[1] * grad(v)[0]) * dx
    a.Assemble()
    f.Assemble()
    gfu = GridFunction(fes)
    gfu.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * f.vec
    return gfu


def halbach_dipole_magnetization(mesh, Br, region="magnet"):
    """In-plane magnetization CF for an ideal k=1 (dipole) HALBACH cylinder over
    ``region``: the easy axis rotates as 2*phi, M = (Br/mu0)(cos 2phi, sin 2phi),
    producing a UNIFORM transverse bore field B = Br ln(Ro/Ri). Hand it to the
    ``magnet_cf`` argument of :func:`solve_planar_magnetostatic`; zero outside
    ``region``. (cos 2phi = (x^2-y^2)/r^2, sin 2phi = 2xy/r^2 -- no atan2 / branch cut.)
    """
    ind = mesh.MaterialCF({region: 1.0}, default=0.0)
    r2 = x * x + y * y + 1e-18
    return ind * (Br / MU0) * CoefficientFunction(((x * x - y * y) / r2, 2.0 * x * y / r2))


def halbach_bore_field(Br, Ri, Ro):
    """Exact uniform bore field of an ideal dipole Halbach cylinder: B = Br ln(Ro/Ri)."""
    return Br * math.log(Ro / Ri)


def cylinder_magnet_axial_field(Br, R, L, z):
    """On-axis B_z of a uniformly axially-magnetized CYLINDER magnet (radius R,
    length L, remanence Br, recoil mu_r = 1), with ``z`` measured from the magnet
    CENTRE. Surface-charge (Coulombian) model:

        B_z(z) = (Br/2) [ (z+L/2)/sqrt(R^2+(z+L/2)^2)
                          - (z-L/2)/sqrt(R^2+(z-L/2)^2) ].

    Centre field B_z(0) = Br L / (2 sqrt(R^2 + (L/2)^2)); far field -> the dipole
    Br R^2 L /(2 z^3). Exact for rigid magnetization (mu_r=1); the bread-and-butter
    closed form for PM pull/holding-force and gap-field estimates. ``z`` may be a
    scalar or any iterable (returns a list)."""
    def _bz(zz):
        zp, zm = zz + 0.5 * L, zz - 0.5 * L
        return 0.5 * Br * (zp / math.sqrt(R * R + zp * zp)
                           - zm / math.sqrt(R * R + zm * zm))
    try:
        return [_bz(float(zz)) for zz in z]
    except TypeError:
        return _bz(float(z))


def solenoid_axial_field(n_turns_per_m, current, radius, length, z):
    """On-axis B_z of a finite air-core SOLENOID (current sheet K = n*I), radius a,
    length L, with ``z`` from the CENTRE. Exact:

        B_z(z) = (mu0 n I / 2) [ (L/2 - z)/sqrt((L/2-z)^2+a^2)
                                + (L/2 + z)/sqrt((L/2+z)^2+a^2) ].

    Centre B_z(0) = mu0 n I (L/2)/sqrt((L/2)^2+a^2) = mu0 n I L/sqrt(L^2+D^2) (D=2a);
    long-coil limit (L>>a) -> mu0 n I; the coil end (z=L/2) -> half the centre value.
    ``z`` may be a scalar or any iterable (returns a list). Smooth on axis (no pole
    corner), so the axi L2-projection extraction recovers it cleanly."""
    a, h = radius, 0.5 * length

    def _bz(zz):
        t1 = (h - zz) / math.sqrt((h - zz) ** 2 + a * a)
        t2 = (h + zz) / math.sqrt((h + zz) ** 2 + a * a)
        return 0.5 * MU0 * n_turns_per_m * current * (t1 + t2)
    try:
        return [_bz(float(zz)) for zz in z]
    except TypeError:
        return _bz(float(z))


def solenoid_inductance_long(n_turns_per_m, radius, length):
    """Long-solenoid inductance L = mu0 n^2 (pi a^2) L (the Nagaoka coefficient -> 1
    as 2a/L -> 0). The exact finite value is k_Nagaoka(2a/L) < 1 times this; use this
    as the L>>a reference and read the true finite value from the FEM energy 2W/I^2."""
    return MU0 * n_turns_per_m ** 2 * math.pi * radius ** 2 * length


def magnetic_circuit_gap_field(N, current, gap, iron_path, mu_r):
    """Reluctance-model gap flux density of an iron magnetic circuit (C/H-frame
    dipole, pot core, relay): Ampere's law around the loop with B continuous gives

        B_gap = mu0 N I / (gap + iron_path / mu_r).

    ``gap`` = total air-gap length [m], ``iron_path`` = mean iron path length [m],
    ``N*current`` = amp-turns. For mu_r -> inf the iron carries no MMF and
    B_gap -> mu0 NI / gap. The FEM value sits a few % BELOW this (fringing/leakage
    widen the effective gap) -- that gap is the physical content the reluctance model
    misses, and what a COMSOL cross-check pins down."""
    return MU0 * N * current / (gap + iron_path / mu_r)


def magnetic_circuit_gap_force(N, current, gap, iron_path, mu_r, area):
    """Reluctance HOLDING force of an iron magnetic-circuit air gap (electromagnet /
    relay / solenoid actuator) -- Maxwell stress on the pole face:

        F = B_gap^2 * area / (2 mu0),   B_gap = magnetic_circuit_gap_field(...),

    so  F = mu0 (N I)^2 area / [2 (gap + iron_path/mu_r)^2]  ~ 1/gap^2 (small gap /
    high mu). ``area`` = pole-face area [m^2] (= face height * stack depth in 2D). The
    FE force (B_gap_FE^2 * area / 2mu0 from the solved gap field) tracks this to a few %
    (fringing/leakage), and F*(gap+iron_path/mu_r)^2 is constant across gaps -- the
    1/g^2 reluctance force law. Validated examples/comsol_class/reluctance_actuator.py."""
    B_gap = magnetic_circuit_gap_field(N, current, gap, iron_path, mu_r)
    return B_gap * B_gap * area / (2.0 * MU0)


def magnetic_circuit_inductance(N, area, gap, iron_path, mu_r):
    """Reluctance-model INDUCTANCE of an N-turn coil on a gapped iron magnetic circuit
    (inductor / transformer / actuator winding):

        L = N^2 / R = N^2 mu0 area / (gap + iron_path / mu_r),

    R the magnetic reluctance of the series gap+iron path (``area`` = core/gap
    cross-section [m^2] = leg width * stack depth in 2D). The inductance twin of
    :func:`magnetic_circuit_gap_force` (same effective gap ``gap + iron_path/mu_r``);
    for mu_r -> inf it is gap-limited, L -> N^2 mu0 area / gap.

    IMPORTANT -- this is the LEAKAGE-FREE estimate (a LOWER BOUND). Unlike the gap
    FORCE (which reads only B_gap, so the reluctance model is good to a few %), the
    coil inductance also links the WINDOW LEAKAGE flux, so a 2D FE solve (radia
    ``inductance_2d``) gives L well ABOVE this -- and the excess GROWS as the gap
    opens (more flux short-cuts the window). On the #27 window frame the FE L is
    ~1.27x (g=4 mm) -> ~1.52x (g=9 mm) this model: the gap that leakage adds is
    exactly the physical content the lumped model misses. Validated
    examples/comsol_class/gapped_core_inductor.py (radia FE vs an independent 2D FE
    reference to ~0.4 %; vs this leakage-free model as the documented lower bound)."""
    return N * N * MU0 * area / (gap + iron_path / mu_r)


def magnetic_circuit_bh_operating_point(mmf, iron_path, h_of_b, gap=0.0):
    """Operating flux density B [T] of a NONLINEAR series magnetic circuit -- the saturating
    counterpart of the constant-mu :func:`magnetic_circuit_inductance` / :func:`magnetic_circuit_gap_force`.
    Ampere's law around the loop, with flux continuity (B the same in iron and gap of equal area):

        NI = H_iron(B) * l_fe + (B/mu0) * gap ,

    ``h_of_b`` is the iron B->H curve [A/m] (a callable, e.g. ``lambda B: 51*B + 2.5*B**15``), ``mmf``
    = N I [A-turns], ``iron_path`` = mean iron path length l_fe [m], ``gap`` = total air-gap length
    [m]. Solved for B by bisection (the right side increases monotonically in B). For ``gap=0`` (a
    closed iron loop) every ampere-turn drops in the iron, so B is fixed by the BH KNEE alone; any
    gap adds the large linear term B*gap/mu0 that quickly de-saturates the iron (the leakage-free
    #27/#39/#42 story, now with a real saturating BH). Validated self-consistently (residual -> 0,
    the constant-mu limit reproduces NI/(l_fe/mu0/mu_r + gap/mu0) exactly) and against a live
    nonlinear FE of a closed iron-frame electromagnet (examples/comsol_class/nonlinear_magnetic_circuit.py)."""
    def rhs(B):
        return h_of_b(B) * iron_path + (B / MU0) * gap
    lo, hi = 0.0, 10.0                                   # B bracket [T]; rhs(0)=0 <= mmf <= rhs(10)
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if rhs(mid) > mmf:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


def pm_circuit_loadline_gap_field(Br, magnet_len, gap, iron_path, mu_r, mu_rec=1.0):
    """PM operating point -- the air-gap flux density of a PERMANENT-MAGNET-driven iron
    circuit (magnet length ``magnet_len`` and the gap/iron in series, equal cross-section).
    Ampere's law around the loop (no current) with B continuous and the recoil line
    B = Br + mu0 mu_rec H in the magnet gives the LOAD-LINE gap field

        B_gap = Br * magnet_len / (magnet_len + mu_rec*gap + iron_path/mu_r).

    The PM counterpart of the coil-driven reluctance circuit (:func:`magnetic_circuit_gap_field`,
    NI replaced by the magnet's Br*magnet_len/mu_rec MMF). The magnet operating point lies on
    its demag/recoil line; a larger gap is a steeper load line -> lower B_gap (more demag).

    LEAKAGE -- this is the leakage-FREE estimate (an UPPER bound): some magnet flux leaks
    across the window instead of reaching the gap, so a 2D FE solve gives B_gap BELOW this,
    and the deficit GROWS as the gap opens (the gap reluctance rises, flux prefers the leak
    path). On the #27 window frame the FE B_gap is ~0.86x (g=2 mm) -> ~0.68x (g=8 mm) this
    model. Validated examples/comsol_class/pm_loadline.py (radia FE vs an independent 2D FE
    reference to ~0.3 %; vs this leakage-free load line as the documented upper bound)."""
    return Br * magnet_len / (magnet_len + mu_rec * gap + iron_path / mu_r)


def two_wire_external_inductance(separation, radius):
    """External inductance per unit length of a TWO-WIRE line (parallel wires, radius a,
    centre separation D) -- the magnetic dual of the two-wire capacitance:

        L_ext = (mu0/pi) acosh(D/2a)   [H/m]   (-> (mu0/pi) ln(D/a) for D >> a).

    With the two-wire C = pi eps0/acosh(D/2a): v = 1/sqrt(LC) = c and
    Z0 = sqrt(L/C) = (1/pi) sqrt(mu0/eps0) acosh(D/2a) = 120 acosh(D/2a) [ohm]. radia gets L
    from the FE field energy 2W/I^2 (``force.magnetic_energy_2d`` over the air -- the +-I
    pair's line-dipole field energy CONVERGES, no outer-boundary dependence, unlike a single
    wire). Validated examples/comsol_class/two_wire_inductance.py."""
    return MU0 / math.pi * math.acosh(separation / (2.0 * radius))


def wire_over_ground_inductance(height, radius):
    """External inductance per unit length of a single wire (radius a) at ``height`` h above a
    PERFECT CONDUCTING GROUND PLANE -- the single-ended / microstrip-style line:

        L_ext = (mu0/2 pi) acosh(h/a)   [H/m]   (-> (mu0/2 pi) ln(2h/a) for h >> a).

    The IMAGE METHOD: the ground (A = 0 on the plane) reflects a -I image wire at -h, so the field
    above the ground is that of the +I/-I pair separated by 2h -- but only HALF the energy lives
    above the plane, so L is exactly HALF the two-wire value at separation 2h:
    ``wire_over_ground_inductance(h, a) = 0.5 * two_wire_external_inductance(2 h, a)``. radia reads
    it from the FE energy ABOVE the ground (Dirichlet A=0 plane); with Z0 = c L_ext it is the
    single-conductor transmission line. Validated examples/comsol_class/wire_over_ground.py (FE vs
    closed form to a few %, the open-domain energy truncation, FE below the h->inf-domain value)."""
    return MU0 / (2.0 * math.pi) * math.acosh(height / radius)


def two_wire_loop_inductance(separation, radius):
    """TOTAL (loop) inductance per unit length of a two-wire line at DC -- the EXTERNAL
    field inductance PLUS the INTERNAL inductance of BOTH conductors:

        L_loop = (mu0/pi) acosh(D/2a)  +  2 * mu0/(8 pi)  =  L_ext + mu0/(4 pi).

    This is the inductance a meter reads for the loop (external + both wires' internal),
    vs the external-only :func:`two_wire_external_inductance` (#30). The internal term
    2*:func:`internal_inductance_round_wire` is radius-independent (mu0/4pi) and is the
    LOW-frequency value -- it rolls off via the skin effect (see
    :func:`skin_effect_internal_inductance_ratio`). radia gets it from the FE field energy
    2(W_air + W_wire1 + W_wire2)/I^2. Validated examples/comsol_class/two_wire_loop_inductance.py."""
    return two_wire_external_inductance(separation, radius) + 2.0 * internal_inductance_round_wire()


def two_wire_force_per_length(current1, current2, separation):
    """Force per unit length between two infinitely-long parallel wires carrying ``current1`` and
    ``current2`` a distance ``separation`` d apart -- Ampere's force law:

        F = mu0 I1 I2 / (2 pi d)   [N/m].

    Parallel LIKE currents (I1 I2 > 0) ATTRACT, opposite currents repel; the returned value is the
    signed force-per-length with the convention that a positive product is an attractive force of
    this magnitude. The Lorentz body force int J x B over one wire: only the OTHER wire's field
    contributes the net force (the symmetric self-field integrates to zero). The building block of
    the image-wire force (:func:`image_force_wire_iron` is this with the iron's image at separation
    2h -> mu0 I^2/(4 pi h)). Validated to 0.001 % against an independent 2D magnetostatic Lorentz-
    force solve (a stored regression reference)."""
    return MU0 * current1 * current2 / (2.0 * math.pi * separation)


def image_force_wire_iron(current, height):
    """Attractive force per unit length on a current wire a distance ``height`` h above an
    IRON (mu -> inf) half-plane, by the METHOD OF IMAGES: the iron mirrors the wire as a
    SAME-sign current at -h (field lines enter the iron perpendicularly), so two parallel
    same-direction currents at separation 2h give

        F = mu0 I^2 / (4 pi h)   [N/m]   (toward the iron).

    (A perfect CONDUCTOR plane mirrors an OPPOSITE-sign current -> repulsion of the same
    magnitude.) radia gets F as the Lorentz body force int J x B over the wire (its symmetric
    self-field nets to zero; the iron's image field provides F). With a large-but-FINITE iron
    block F is ~4 % below this (finite block != infinite half-plane); -> exact as the block
    grows. Validated examples/comsol_class/image_force.py."""
    return MU0 * current * current / (4.0 * math.pi * height)


def internal_inductance_round_wire():
    """INTERNAL inductance per unit length of a round wire with a UNIFORM (DC) current,
    from the field energy stored INSIDE the conductor:

        L_int = mu0/(8 pi)   [H/m]   -- independent of the wire radius.

    (B(r)=mu0 I r/(2 pi R^2) inside -> W_in = mu0 I^2/(16 pi) -> L_int = 2 W_in/I^2.) The
    internal companion to the EXTERNAL inductances (coax (mu0/2pi)ln(b/a),
    :func:`two_wire_external_inductance`); it sets the DC/low-frequency inductance and rolls
    OFF as the skin effect expels current at high frequency. radia gets it from the wire-interior
    field energy ``2*force.magnetic_energy_2d(B, mesh, 'wire')/I^2`` (boundary-independent --
    Ampere fixes the interior B from the enclosed current). Validated
    examples/comsol_class/internal_inductance.py."""
    return MU0 / (8.0 * math.pi)


def skin_effect_resistance_ratio(q):
    """AC/DC RESISTANCE ratio Rac/Rdc of a round wire at skin-effect parameter
    ``q = a*sqrt(omega*mu0*sigma) = sqrt(2)*a/delta`` (a = radius, delta = skin depth):

        Rac/Rdc = (q/2) [ber(q) bei'(q) - bei(q) ber'(q)] / [ber'(q)^2 + bei'(q)^2]

    (Kelvin ber/bei). q -> 0 : Rac/Rdc -> 1 (uniform DC current); q >> 1 : -> q/2 + ... (current
    crowds into the ~delta-thick skin). The real part of the round-wire INTERNAL impedance.
    Needs scipy. Validated vs the radia eddy solve (tests/test_planar_eddy.py)."""
    from scipy.special import kelvin
    be, ke, bep, kep = kelvin(q)
    ber, bei, berp, beip = be.real, be.imag, bep.real, bep.imag
    return (q / 2.0) * (ber * beip - bei * berp) / (berp ** 2 + beip ** 2)


def skin_effect_internal_inductance_ratio(q):
    """AC INTERNAL-inductance ROLL-OFF ratio L_int(omega)/L_int(0) of a round wire at
    ``q = a*sqrt(omega*mu0*sigma) = sqrt(2)*a/delta`` (L_int(0) = mu0/(8 pi),
    :func:`internal_inductance_round_wire`):

        L_int/L_int_dc = (4/q) [ber(q) ber'(q) + bei(q) bei'(q)] / [ber'(q)^2 + bei'(q)^2]

    (Kelvin ber/bei). q -> 0 : -> 1 (the DC mu0/8pi, #36); q >> 1 : -> 4/q -> 0 (the skin effect
    expels current from the core, so the interior stores less magnetic energy). The imaginary
    part /omega of the round-wire INTERNAL impedance -- the inductive twin of
    :func:`skin_effect_resistance_ratio`; together Z_internal(omega) = Rdc*Rac_ratio +
    j*omega*(mu0/8pi)*Lint_ratio. Needs scipy. Validated vs the radia eddy solve
    (examples/comsol_class/ac_internal_inductance.py)."""
    from scipy.special import kelvin
    be, ke, bep, kep = kelvin(q)
    ber, bei, berp, beip = be.real, be.imag, bep.real, bep.imag
    return (4.0 / q) * (ber * berp + bei * beip) / (berp ** 2 + beip ** 2)


def skin_effect_slab_attenuation(sigma, frequency, half_thickness, mu_r=1.0):
    r"""Mid-plane field attenuation A_z(0)/A_z(surface) of a conducting slab x in [-a, a]
    in a tangential AC field (skin effect) -- the closed form of the eddy diffusion
    ``lap(A_z) = j*omega*mu0*mu_r*sigma*A_z`` with the surface phasor A_z=A0 on x=+/-a:

        A_z(0)/A0 = 1 / cosh(gamma a),   gamma = (1+j)/delta,
        delta = sqrt(2 / (omega mu0 mu_r sigma))    (skin depth).

    Complex: magnitude = attenuation, phase = the eddy phase lag.  a << delta -> ~1 (the
    field penetrates); a >> delta -> ~2 exp(-a/delta) (expelled to the ~delta-thick skin).
    The slab counterpart of the round-wire Kelvin-function ratios
    (:func:`skin_effect_resistance_ratio`), relevant to lamination / sheet eddy loss.
    Validated against the complex NGSolve eddy solve (tests/test_skin_effect_slab.py) and
    an independent frequency-domain magnetics cross-check."""
    import cmath
    mu0 = 4e-7 * math.pi
    omega = 2.0 * math.pi * frequency
    delta = math.sqrt(2.0 / (omega * mu0 * mu_r * sigma))
    gamma = (1 + 1j) / delta
    return 1.0 / cmath.cosh(gamma * half_thickness)


def cogging_torque_order(slots, poles):
    r"""Spatial order (cycles per MECHANICAL revolution) of the fundamental cogging torque of a
    slotted PM machine: ``N_c = LCM(slots, poles)`` (poles = pole COUNT = 2p).  The cogging
    PERIOD is 360/N_c degrees (:func:`cogging_period_deg`); a higher LCM (fractional-slot, e.g.
    12-slot/10-pole -> 60) gives smaller, higher-frequency cogging.  Design rule: slots != poles
    (equal -> huge cogging).  Confirmed against an independent slotted 4-pole / 6-slot rotor FEA
    -- the cogging waveform is dominated by exactly one cycle per 360/LCM degrees."""
    from math import gcd
    return slots * poles // gcd(slots, poles)


def cogging_period_deg(slots, poles):
    """Mechanical-angle period [deg] of the cogging torque = 360 / LCM(slots, poles)
    (:func:`cogging_torque_order`).  One cogging cycle spans this angle; a one-slot-pitch skew
    (360/slots) nulls the fundamental (sinc(pi)=0, :func:`skew_factor`)."""
    return 360.0 / cogging_torque_order(slots, poles)


def annular_radial_resistance(sigma, r_inner, r_outer, length=1.0):
    r"""DC resistance [ohm] of an annular conductor carrying RADIAL current between the inner
    ring (r_inner) and the outer ring (r_outer), axial length ``length``, conductivity ``sigma``:

        R = ln(r_outer / r_inner) / (2 pi sigma length),

    the closed form of -div(sigma grad V)=0 with V fixed on the two rings (V(r) ~ ln r).  The
    radial-conduction counterpart of the axial L/(sigma A) -- relevant to disc/annulus electrodes
    and radial busbars.  Validated against an NGSolve conduction solve and an independent
    cross-check (tests/test_annular_resistance.py)."""
    import math
    return math.log(r_outer / r_inner) / (2.0 * math.pi * sigma * length)


def demag_operating_field(Br, demag_factor, mu_r=1.0):
    r"""Internal operating field H_op along the magnetization [A/m] of a permanent magnet at its
    no-load operating point: with demagnetizing factor N (=1/2 transverse cylinder, 1/3 sphere),
    remanence Br and recoil mu_r,

        B_in = Br (1-N) / (1 + N (mu_r-1)),    H_op = (B_in - Br) / (mu0 mu_r).

    For mu_r=1 this is H_op = -N Br/mu0 = -N M (e.g. -Br/(2 mu0) for the transverse cylinder).
    In a 2-D field the operating point is DISTRIBUTED -- corners of non-ellipsoidal magnets
    deviate -- so evaluate H.m_hat PER ELEMENT and compare each to the knee
    (:func:`demag_margin`).  Confirmed vs a transverse PM-cylinder FEA (bulk H_x == -Hc/2)."""
    import math
    mu0 = 4e-7 * math.pi
    B_in = Br * (1.0 - demag_factor) / (1.0 + demag_factor * (mu_r - 1.0))
    return (B_in - Br) / (mu0 * mu_r)


def demag_margin(H_operating, H_knee):
    """Irreversible-demagnetization margin [A/m] = H_operating - H_knee (both the m_hat component,
    negative = demagnetizing).  Positive => the operating point sits ABOVE the knee (safe); a
    stator / short-circuit field that drives the WORST element's H.m_hat below H_knee causes
    permanent loss of magnetization.  Pair with the per-element worst H.m_hat, not just the bulk
    :func:`demag_operating_field` (corners reach the knee first)."""
    return H_operating - H_knee


def park_transform(a, b, c, theta_e):
    r"""Amplitude-invariant (peak-preserving) Park transform of three phase quantities a,b,c into
    the rotor (d,q) frame at electrical angle theta_e [rad]:

        d = (2/3)[a cos(th) + b cos(th-2pi/3) + c cos(th+2pi/3)]
        q = -(2/3)[a sin(th) + b sin(th-2pi/3) + c sin(th+2pi/3)]

    theta_e is the rotor d-axis position measured FROM the phase-A MAGNETIC axis (the centre of
    phase A's coil span, NOT the slot where the A conductors sit -- for a full-pitch belt the axis
    leads the belt by 90 elec).  Inverse: :func:`inverse_park`.  Applied to FE phase flux linkages
    it yields lambda_d, lambda_q (-> Ld=lambda_d/i_d at i_q=0, Lq=lambda_q/i_q at i_d=0)."""
    import math
    d = (2.0 / 3.0) * (a * math.cos(theta_e) + b * math.cos(theta_e - 2 * math.pi / 3)
                       + c * math.cos(theta_e + 2 * math.pi / 3))
    q = -(2.0 / 3.0) * (a * math.sin(theta_e) + b * math.sin(theta_e - 2 * math.pi / 3)
                        + c * math.sin(theta_e + 2 * math.pi / 3))
    return d, q


def inverse_park(d, q, theta_e):
    r"""Inverse amplitude-invariant Park: (d,q) at electrical angle theta_e [rad] -> phase a,b,c
    (peak amplitudes).  ``a = d cos(th) - q sin(th)``, etc.  Round-trips with :func:`park_transform`
    (park(inverse_park(d,q,th),th) == (d,q)).  Use to set the three phase currents that realise a
    commanded (i_d,i_q) operating point in an FE saliency/torque extraction."""
    import math
    a = d * math.cos(theta_e) - q * math.sin(theta_e)
    b = d * math.cos(theta_e - 2 * math.pi / 3) - q * math.sin(theta_e - 2 * math.pi / 3)
    c = d * math.cos(theta_e + 2 * math.pi / 3) - q * math.sin(theta_e + 2 * math.pi / 3)
    return a, b, c


def dq_flux_torque(pole_pairs, lam_d, lam_q, i_d, i_q):
    r"""Instantaneous electromagnetic torque [N.m] from the dq flux linkages, general (no linear-
    decomposition assumption):

        T = (3/2) p (lambda_d i_q - lambda_q i_d).

    This is the co-energy route -- it uses the FE flux linkages directly, so it captures the
    FUNDAMENTAL dq interaction exactly.  Compare with the lumped :func:`dq_torque` (which assumes
    lambda_d = lambda_pm + Ld i_d, lambda_q = Lq i_q) and with a Maxwell-stress-tensor torque: for
    a salient machine the stress-tensor torque RIPPLES about this value at the pole-harmonic
    frequency, and its mean over one rotor period equals it."""
    return 1.5 * pole_pairs * (lam_d * i_q - lam_q * i_d)


def coil_pair_axial_field(N, current, radius, separation, z):
    """On-axis B_z of a coaxial PAIR of N-turn circular loops (radius a) on the
    z-axis at z = +/- separation/2, each carrying ``current`` in the SAME sense:

        B_z(z) = (mu0 N I a^2/2)[((z-s/2)^2+a^2)^-3/2 + ((z+s/2)^2+a^2)^-3/2].

    HELMHOLTZ pair = separation == radius (s = a): d^2B/dz^2 = 0 at the centre, so
    the field is maximally uniform, with centre value B0 = (4/5)^(3/2) mu0 N I/a =
    0.7155 mu0 N I/a. (Anti-Helmholtz = opposite senses -> a linear gradient, for a
    quadrupole/magnetic-bottle; pass -current to one loop.) ``z`` from the mid-plane;
    scalar or iterable."""
    a, s = radius, separation
    def _bz(zz):
        zm, zp = zz - 0.5 * s, zz + 0.5 * s
        return 0.5 * MU0 * N * current * a * a * (
            (zm * zm + a * a) ** -1.5 + (zp * zp + a * a) ** -1.5)
    try:
        return [_bz(float(zz)) for zz in z]
    except TypeError:
        return _bz(float(z))


def solve_planar_magnetostatic_nonlinear(mesh, nu_of_B, Jz=None, magnets=None,
                                         order=2, dirichlet="outer", relax=0.5,
                                         max_iter=80, tol=1e-5, min_iter=3):
    """2D PLANAR A_z magnetostatics with SATURATING IRON -- Picard fixed point on
    a field-dependent reluctivity nu(|B|). The planar twin of
    :func:`solve_magnetostatic_nonlinear` (FEMM ``prob1big`` nonlinear).

    ``nu_of_B`` : callable(B_cf) -> nu_cf, where ``B_cf = CF((grad(A)[1], -grad(A)[0]))``
    is the current flux density. Typical saturating iron (Froehlich/Kennelly
    B = H/(alpha+beta H), so nu(B) = alpha/(1-beta B)) clamped below B_sat::

        def nu_of_B(B):
            Bmag = sqrt(InnerProduct(B, B) + 1e-20)
            Bc = IfPos(Bmag - 0.98*Bsat, 0.98*Bsat, Bmag)
            return iron_ind * alpha/(1 - beta*Bc) + (1 - iron_ind) * NU0

    ``Jz`` (source current density) and/or ``magnets`` ({mat:(Hc,phi_deg)}) drive
    the problem. ``relax`` under-relaxes A <- (1-w)A + w A_new; keep w ~ 0.3-0.5
    for saturating curves (too large diverges). Validated by Ampere's law on a
    wire-in-iron-annulus: B(r)=BH(I/2 pi r) to <0.1 % away from the inner
    interface (tests/test_planar_nonlinear.py).

    Returns the converged H1 GridFunction ``A_z``.
    """
    fes = H1(mesh, order=order, dirichlet=dirichlet)
    u, v = fes.TnT()
    src = None
    if magnets is not None:
        src = planar_magnet_source(mesh, magnets)
    gfu = GridFunction(fes)
    gfu.vec[:] = 0.0
    prev = None
    for it in range(max_iter):
        B = CoefficientFunction((grad(gfu)[1], -grad(gfu)[0]))
        a = BilinearForm(fes, symmetric=True)
        a += nu_of_B(B) * grad(u) * grad(v) * dx
        f = LinearForm(fes)
        if Jz is not None:
            f += Jz * v * dx
        if src is not None:
            f += (src[0] * grad(v)[1] - src[1] * grad(v)[0]) * dx
        a.Assemble()
        f.Assemble()
        gnew = GridFunction(fes)
        gnew.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * f.vec
        gfu.vec.data = (1.0 - relax) * gfu.vec + relax * gnew.vec
        Bc = CoefficientFunction((grad(gfu)[1], -grad(gfu)[0]))
        cur = Integrate(InnerProduct(Bc, Bc) * dx, mesh)
        if prev is not None and it + 1 >= min_iter and abs(cur - prev) < tol * max(abs(cur), 1e-30):
            break
        prev = cur
    return gfu


def solve_planar_eddy(mesh, nu, sigma, omega, driven_region=None,
                      total_current=None, applied_Ez=None, Jz=None, order=3,
                      dirichlet="outer", dirichlet_value=0.0):
    """2D PLANAR time-harmonic eddy currents (complex A_z) -- FEMM ``prob2big`` analog.

    Solves  -div(nu grad A_z) + j w sigma A_z = J_src .  The eddy current in
    conductors (sigma>0) is J_eddy = -j w sigma A_z (added automatically). Drive
    modes:

    * ``Jz`` : imposed source current density CF [A/m^2] (e.g. a stranded coil).
    * ``dirichlet_value`` : DRIVEN-SURFACE phasor -- a NON-zero Dirichlet A_z =
      dirichlet_value imposed on the ``dirichlet`` boundary (default 0 = magnetic
      insulation). The slab skin-effect (surface flux A0 driven on the conductor
      faces, the eddy reaction attenuating it inward) is A_z(0)/A0 = 1/cosh(gamma a),
      gamma=(1+j)/delta, delta=sqrt(2/(w mu sigma)).
    * ``driven_region`` + ``total_current`` : CURRENT-DRIVEN solid conductor. A
      NumberSpace scalar Vc (axial driving field -dV/dz) is added and constrained
      so  int_region sigma (-j w A_z + Vc) dA = total_current -- net current is
      fixed while J redistributes (skin/proximity). FEMM "current driven" circuit.
    * ``applied_Ez`` : VOLTAGE-DRIVEN conductor. The axial field Vc = applied_Ez
      [V/m] is PRESCRIBED (a known source  f += sigma*Vc*v, no NumberSpace); the
      net current I = int sigma(-j w A_z + Vc) then follows, giving the per-length
      impedance Z = applied_Ez / I. FEMM "voltage driven" circuit.

    Returns the compound GridFunction (``components = (A_z, Vc)``) in current-driven
    mode, else the H1 GridFunction ``A_z`` (E_z = -j w A_z (+ applied_Ez in the
    conductor), J_z = sigma E_z). Flux density ``B = CF((grad(A_z)[1], -grad(A_z)[0]))``.

    Validated: round-wire skin-effect Rac/Rdc vs the Kelvin-function (ber/bei)
    formula to 0.07 % (current-driven) and voltage-driven Z to 0.07 %
    (tests/test_planar_eddy.py, test_planar_eddy_voltage.py) at q=4.
    """
    if driven_region is not None and total_current is not None:
        fes = H1(mesh, order=order, complex=True, dirichlet=dirichlet) \
            * NumberSpace(mesh, complex=True)
        (Az, Vc), (dA, dV) = fes.TnT()
        a = BilinearForm(fes)
        a += nu * grad(Az) * grad(dA) * dx
        a += 1j * omega * sigma * Az * dA * dx          # eddy reaction
        a += -sigma * Vc * dA * dx                       # Vc drives A in conductor
        a += -1j * omega * sigma * Az * dV * dx          # net-current constraint
        a += sigma * Vc * dV * dx
        f = LinearForm(fes)
        if Jz is not None:
            f += Jz * dA * dx
        a.Assemble()
        f.Assemble()
        f.vec.FV().NumPy()[fes.Range(1).start] += complex(total_current)
        gfu = GridFunction(fes)
        gfu.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="umfpack") * f.vec
        return gfu

    fes = H1(mesh, order=order, complex=True, dirichlet=dirichlet)
    u, v = fes.TnT()
    a = BilinearForm(fes)
    a += nu * grad(u) * grad(v) * dx + 1j * omega * sigma * u * v * dx
    f = LinearForm(fes)
    if Jz is not None:
        f += Jz * v * dx
    if applied_Ez is not None:           # voltage-driven: J_src = sigma * Vc
        f += sigma * applied_Ez * v * dx
    a.Assemble()
    f.Assemble()
    gfu = GridFunction(fes)
    if dirichlet_value:                  # DRIVEN-SURFACE: A_z = dirichlet_value on `dirichlet` (skin slab)
        gfu.Set(CoefficientFunction(complex(dirichlet_value)), BND,
                definedon=mesh.Boundaries(dirichlet))
        r = f.vec - a.mat * gfu.vec
        gfu.vec.data += a.mat.Inverse(fes.FreeDofs(), inverse="umfpack") * r
    else:
        gfu.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="umfpack") * f.vec
    return gfu


def solve_planar_eddy_multi(mesh, nu, sigma, omega, conductors,
                            connection="series", total_current=None,
                            currents=None, Jz=None, order=3, dirichlet="outer"):
    """2D PLANAR time-harmonic eddy currents with a MULTI-CONDUCTOR circuit --
    skin AND proximity effect (FEMM "circuit" grouping several solid conductors).

    The conductors share the SAME complex A_z field, so each one's eddy currents
    are driven by the field of all the others -> proximity effect, on top of each
    one's own skin effect. The net current of each conductor is constrained by an
    axial-E (``Vc = -dV/dz``) NumberSpace unknown, exactly as in the single-
    conductor :func:`solve_planar_eddy`, but one unknown per conductor.

    Parameters
    ----------
    conductors : list of material names, each a solid conductor (sigma > 0 there).
    connection :
        * ``"series"``  -- every conductor carries the SAME net current
          ``total_current`` (one Vc_k per conductor; the per-conductor voltages
          differ). Typical multi-turn winding / adjacent bus-bars.
        * ``"parallel"`` -- all conductors share ONE Vc (one common voltage); the
          net currents redistribute and SUM to ``total_current``. (Equivalent to
          :func:`solve_planar_eddy` with sigma covering every conductor.)
    total_current : circuit current I (series: per conductor; parallel: the sum).
    currents      : optional per-conductor net currents [list, len==conductors],
                    overriding the equal-current series default (independent
                    forced currents, one Vc each).
    Jz            : optional imposed source current density CF [A/m^2].
    order, dirichlet : as in :func:`solve_planar_eddy`.

    Returns the compound GridFunction. Components:
        * series / independent : ``(A_z, Vc_0, Vc_1, ...)`` -- one Vc per conductor;
          conductor k field  E_z = -j w A_z + Vc_k, loss
          ``ohmic_loss_2d(E_z, mesh, sigma, region=conductors[k])``.
        * parallel : ``(A_z, Vc)`` -- shared Vc.
    The circuit AC resistance is ``Rac = 2*sum_k P_k / |I|^2``.

    Validated (tests/test_planar_proximity.py): reduces to the single-wire
    Kelvin Rac; two well-separated wires in series give 2x the isolated Rac; and
    bringing them together raises Rac (proximity), all against the exact ber/bei
    round-wire reference.
    """
    if connection == "parallel":
        fes = H1(mesh, order=order, complex=True, dirichlet=dirichlet) \
            * NumberSpace(mesh, complex=True)
        (Az, Vc), (dA, dV) = fes.TnT()
        cond_cf = dx(definedon=mesh.Materials("|".join(conductors)))
        a = BilinearForm(fes)
        a += nu * grad(Az) * grad(dA) * dx
        a += 1j * omega * sigma * Az * dA * dx
        a += -sigma * Vc * dA * cond_cf
        a += -1j * omega * sigma * Az * dV * cond_cf
        a += sigma * Vc * dV * cond_cf
        f = LinearForm(fes)
        if Jz is not None:
            f += Jz * dA * dx
        a.Assemble()
        f.Assemble()
        f.vec.FV().NumPy()[fes.Range(1).start] += complex(total_current)
        gfu = GridFunction(fes)
        gfu.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="umfpack") * f.vec
        return gfu

    # series / independent: one NumberSpace per conductor
    if currents is None:
        currents = [total_current] * len(conductors)
    h1 = H1(mesh, order=order, complex=True, dirichlet=dirichlet)
    fes = FESpace([h1] + [NumberSpace(mesh, complex=True) for _ in conductors])
    trials = fes.TrialFunction()
    tests = fes.TestFunction()
    Az, dA = trials[0], tests[0]
    a = BilinearForm(fes)
    a += nu * grad(Az) * grad(dA) * dx
    a += 1j * omega * sigma * Az * dA * dx
    for k, reg in enumerate(conductors):
        Vc_k, dV_k = trials[k + 1], tests[k + 1]
        dxr = dx(definedon=mesh.Materials(reg))
        a += -sigma * Vc_k * dA * dxr            # Vc_k drives A only in conductor k
        a += -1j * omega * sigma * Az * dV_k * dxr  # net-current constraint of k
        a += sigma * Vc_k * dV_k * dxr
    f = LinearForm(fes)
    if Jz is not None:
        f += Jz * dA * dx
    a.Assemble()
    f.Assemble()
    fv = f.vec.FV().NumPy()
    for k in range(len(conductors)):
        fv[fes.Range(k + 1).start] += complex(currents[k])
    gfu = GridFunction(fes)
    gfu.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="umfpack") * f.vec
    return gfu


def solve_planar_eddy_nonlinear(mesh, nu_of_B, sigma, omega, Jz=None, order=3,
                                dirichlet="outer", relax=0.5, max_iter=80,
                                tol=1e-5, min_iter=3):
    """2D PLANAR NONLINEAR time-harmonic eddy currents (FEMM nonlinear AC) --
    Picard fixed point with an amplitude-dependent reluctivity nu(|B|).

    Solves  -div(nu(|B|) grad A_z) + j w sigma A_z = Jz  where the saturating
    material reluctivity is evaluated at the PHASOR AMPLITUDE
    ``|B| = sqrt(B . conj(B))`` (the effective-permeability harmonic-balance
    approximation FEMM uses for nonlinear AC -- one complex solve per Picard
    sweep, the standard amplitude-based effective mu).

    ``nu_of_B`` : callable(Bmag_cf) -> nu_cf, taking the SCALAR amplitude
    ``Bmag`` (not the vector, unlike the static nonlinear solver), e.g.::

        def nu_of_B(Bmag):
            Bc = IfPos(Bmag - 0.98*Bsat, 0.98*Bsat, Bmag)
            return iron_ind * alpha/(1 - beta*Bc) + (1 - iron_ind) * NU0

    ``relax`` under-relaxes (0.3-0.5 for hard saturation). Reduces EXACTLY to
    :func:`solve_planar_eddy` (Jz mode) when nu_of_B is constant; saturation
    makes |B| grow sub-linearly with drive (tests/test_planar_eddy_nonlinear.py).

    Returns the converged complex H1 GridFunction ``A_z``.
    """
    fes = H1(mesh, order=order, complex=True, dirichlet=dirichlet)
    u, v = fes.TnT()
    gfu = GridFunction(fes)
    gfu.vec[:] = 0.0
    prev = None
    for it in range(max_iter):
        B = CoefficientFunction((grad(gfu)[1], -grad(gfu)[0]))
        Bmag = sqrt((B[0] * Conj(B[0]) + B[1] * Conj(B[1])).real + 1e-30)
        a = BilinearForm(fes)
        a += nu_of_B(Bmag) * grad(u) * grad(v) * dx
        a += 1j * omega * sigma * u * v * dx
        f = LinearForm(fes)
        if Jz is not None:
            f += Jz * v * dx
        a.Assemble()
        f.Assemble()
        gnew = GridFunction(fes)
        gnew.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="umfpack") * f.vec
        gfu.vec.data = (1.0 - relax) * gfu.vec + relax * gnew.vec
        Bc = CoefficientFunction((grad(gfu)[1], -grad(gfu)[0]))
        cur = Integrate((Bc[0] * Conj(Bc[0]) + Bc[1] * Conj(Bc[1])).real * dx, mesh)
        if prev is not None and it + 1 >= min_iter and abs(cur - prev) < tol * max(abs(cur), 1e-30):
            break
        prev = cur
    return gfu


def solve_axi_magnetostatic(mesh, nu, Jr=None, magnets=None, order=2,
                              dirichlet="axis|outer"):
    """Axisymmetric A_phi magnetostatics via H1Henrotte FESpace (FEMM prob3big axi).

    Weak form in the meridional (r,z) plane (r = NGSolve x, dx = dr dz):
        K += nu*(1/r)*(r*dA/dr+A)*(r*dv/dr+v)*dx + nu*r*dA/dz*dv/dz*dx
    Flux density: B_z = grad(u)[0] + u/r,  B_r = -grad(u)[1]

    Parameters
    ----------
    mesh      : 2D (r,z) mesh; r=0 axis MUST be included in ``dirichlet``.
    nu        : reluctivity CF (NU0/mu_r per material).
    Jr        : phi-direction current density [A/m^2] CF or None.
    magnets   : {material: (Hc, theta_deg)} where theta is from the r-axis [deg].
                theta=90 => axial (+z) magnetization (FEMM MagDir convention).
    order     : H1Henrotte order 1 or 2 (default 2).
    dirichlet : boundary tag for A_phi = 0 (must include the r=0 axis).

    Returns the H1Henrotte GridFunction ``gfu`` (A_phi at DOFs). Validated:
    magnetized sphere B_in = 2 mu0 mu_r Hc/(mu_r+2) to -0.05 %
    (tests/test_axi_magnetostatic.py).
    """
    from radia.radia_axifemm import H1Henrotte
    fes = H1Henrotte(mesh, order=order, dirichlet=dirichlet)
    u, v = fes.TnT()
    r = x
    a = BilinearForm(fes, symmetric=True)
    a += nu * (1.0 / r) * (r * grad(u)[0] + u) * (r * grad(v)[0] + v) * dx
    a += nu * r * grad(u)[1] * grad(v)[1] * dx
    f = LinearForm(fes)
    if Jr is not None:
        f += Jr * r * v * dx
    if magnets is not None:
        for mat, (Hc, theta_deg) in magnets.items():
            th = math.radians(theta_deg)
            reg = mesh.MaterialCF({mat: 1.0}, default=0.0)
            f += reg * Hc * math.sin(th) * (r * grad(v)[0] + v) * dx
            f += -reg * Hc * math.cos(th) * r * grad(v)[1] * dx
    a.Assemble()
    f.Assemble()
    gfu = GridFunction(fes)
    gfu.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * f.vec
    return gfu


def solve_axi_eddy(mesh, nu, sigma, omega, driven_region=None, total_current=None,
                    applied_Vc=None, Jr=None, order=3, dirichlet="axis|outer"):
    """Axisymmetric time-harmonic eddy currents via H1Henrotte (FEMM prob2big axi).

    Solves (K + j*w*M) A_phi = f where K is the H1Henrotte stiffness and
    M = sigma*r (r-weighted sigma-mass for the 3D volume element 2*pi*r dr dz).

    Drive modes:
    * ``Jr`` : imposed phi-direction current density [A/m^2] CF.
    * ``driven_region`` + ``total_current`` : CURRENT-DRIVEN conductor.
      A NumberSpace scalar Vc (= r*E_phi, constant) is added; the constraint
      int sigma*(-j*w*r*A + Vc)*dx = I/(2*pi) fixes the total 3D ring current
      I = 2*pi * int J_phi * r dr dz.
    * ``applied_Vc`` : VOLTAGE-DRIVEN. Vc = r*E_phi is PRESCRIBED (a known
      constant [V/turn/radian]); net current I = 2*pi*int sigma*(-j*w*r*A+Vc)*dx.

    Returns compound gfu (A_phi, Vc) in current-driven mode, else H1Henrotte gfu.
    Flux density: B_z = grad(u)[0] + u/r, B_r = -grad(u)[1]

    Validated: static limit (sigma=0) agrees with solve_axi_magnetostatic;
    time-harmonic via Cu-disk eddy eigenvalue tau_1 = 224.31 us, 0.27 % gap to
    BEM-Foster (see examples/axifemm/research/verification/test_disk_eigenvalue.py).
    """
    from radia.radia_axifemm import H1Henrotte
    r = x

    if driven_region is not None and total_current is not None:
        fes_h1 = H1Henrotte(mesh, order=order, complex=True, dirichlet=dirichlet)
        fes = fes_h1 * NumberSpace(mesh, complex=True)
        (Az, Vc), (dA, dV) = fes.TnT()
        a = BilinearForm(fes)
        a += nu * (1.0 / r) * (r * grad(Az)[0] + Az) * (r * grad(dA)[0] + dA) * dx
        a += nu * r * grad(Az)[1] * grad(dA)[1] * dx
        a += 1j * omega * sigma * r * Az * dA * dx
        # Vc = r*E_phi (constant); J_phi = sigma*(-jw A + Vc/r), so source is sigma*Vc*v*dx
        a += -sigma * Vc * dA * dx
        # current constraint: int sigma*(-jw*r*A + Vc)*dx = I/(2pi)
        a += -1j * omega * sigma * r * Az * dV * dx
        a += sigma * Vc * dV * dx
        f = LinearForm(fes)
        if Jr is not None:
            f += Jr * r * dA * dx
        a.Assemble()
        f.Assemble()
        f.vec.FV().NumPy()[fes.Range(1).start] += complex(total_current) / (2.0 * math.pi)
        gfu = GridFunction(fes)
        gfu.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="umfpack") * f.vec
        return gfu

    fes = H1Henrotte(mesh, order=order, complex=True, dirichlet=dirichlet)
    u, v = fes.TnT()
    a = BilinearForm(fes)
    a += nu * (1.0 / r) * (r * grad(u)[0] + u) * (r * grad(v)[0] + v) * dx
    a += nu * r * grad(u)[1] * grad(v)[1] * dx
    a += 1j * omega * sigma * r * u * v * dx
    f = LinearForm(fes)
    if Jr is not None:
        f += Jr * r * v * dx
    if applied_Vc is not None:
        # Vc = r*E_phi (constant); source: sigma*(Vc/r)*v*r dx = sigma*Vc*v*dx
        f += sigma * applied_Vc * v * dx
    a.Assemble()
    f.Assemble()
    gfu = GridFunction(fes)
    gfu.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="umfpack") * f.vec
    return gfu


def solve_axi_magnetostatic_nonlinear(mesh, nu_of_B, Jr=None, order=2,
                                        dirichlet="axis|outer", relax=0.5,
                                        max_iter=80, tol=1e-5, min_iter=3):
    """Axisymmetric A_phi magnetostatics with SATURATING B-H via Picard.

    Same iteration as ``solve_planar_magnetostatic_nonlinear`` but for the
    axisymmetric H1Henrotte FESpace.

    ``nu_of_B`` : callable(B_cf) -> nu_cf where
        B_cf = CF((grad(u)[0]+u/r, -grad(u)[1]))  is the current (Bz, Br) field.

    Example Froehlich/Kennelly curve (no gradient-recovery issue in axi since we
    sample away from the axis)::

        def nu_of_B(B):
            Bmag = sqrt(InnerProduct(B, B) + 1e-20)
            Bc = IfPos(Bmag - 0.98*Bsat, 0.98*Bsat, Bmag)
            return iron_ind * alpha/(1 - beta*Bc) + (1 - iron_ind) * NU0

    Returns the converged H1Henrotte GridFunction A_phi.
    """
    from radia.radia_axifemm import H1Henrotte
    fes = H1Henrotte(mesh, order=order, dirichlet=dirichlet)
    u_trial, v = fes.TnT()
    r = x
    gfu = GridFunction(fes)
    gfu.vec[:] = 0.0
    prev = None
    for it in range(max_iter):
        B = CoefficientFunction((grad(gfu)[0] + gfu / r, -grad(gfu)[1]))
        nu_cur = nu_of_B(B)
        a = BilinearForm(fes, symmetric=True)
        a += nu_cur * (1.0 / r) * (r * grad(u_trial)[0] + u_trial) * (r * grad(v)[0] + v) * dx
        a += nu_cur * r * grad(u_trial)[1] * grad(v)[1] * dx
        f = LinearForm(fes)
        if Jr is not None:
            f += Jr * r * v * dx
        a.Assemble()
        f.Assemble()
        gnew = GridFunction(fes)
        gnew.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * f.vec
        gfu.vec.data = (1.0 - relax) * gfu.vec + relax * gnew.vec
        B_cur = CoefficientFunction((grad(gfu)[0] + gfu / r, -grad(gfu)[1]))
        cur = Integrate(InnerProduct(B_cur, B_cur) * dx, mesh)
        if prev is not None and it + 1 >= min_iter and abs(cur - prev) < tol * max(abs(cur), 1e-30):
            break
        prev = cur
    return gfu


def solve_axi_eddy_harmonic(mesh, mu_cf, sigma_cf, omega, applied_A,
                            order=1, dirichlet="axis|outer"):
    """Forced time-harmonic axisymmetric eddy currents via the closed-form C++
    AxiHenrotte integrators + a complex scipy solve -- the CORRECT axi eddy path.

    H1Henrotte's symbolic NGSolve weak-form assembly CANNOT build the complex
    ``nu*stiffness + j*w*sigma*mass`` system: the AxiHenrotte DiffOps are
    real-coefficient only, so mixing a real stiffness with the complex j*w*sigma
    mass in one BilinearForm fails (this is why the grad-weak-form solve_axi_eddy
    only works at sigma=0 / via the K,M eigenvalue path).  Instead we assemble the
    REAL closed-form element matrices

        K = AxiHenrotteStiffnessBFI(mu_cf)      (note: PERMEABILITY mu, not nu)
        M = AxiHenrotteSigmaMassBFI(sigma_cf)

    -- the SAME validated integrators behind the Cu-disk tau_1 eigenvalue --
    symmetrise them (BFI stores one triangle), form the complex
    ``S = K + j*w*M`` in scipy and solve  S x = b  with
    ``b_i = int sigma * applied_A * v_i`` (eddy driven by an imposed A_phi field,
    e.g. ``applied_A = B0*x/2`` for a uniform axial B0 -- the Kameari convention).

    ORDER=1 ONLY: the P2 AxiHenrotte element has an axis-singularity NaN.

    Returns ``(gfu, P_eddy)``:
      * ``P_eddy`` -- time-averaged 3-D eddy loss [W] = 0.5*w^2 * Re(x^H M x)
        ( = 0.5 sigma w^2 int |A|^2 2 pi r dr dz ), straight from the matrix M.
        This is the RELIABLE scalar output (matrix-based, no field eval).
      * ``gfu`` -- complex H1Henrotte GridFunction holding the solution DOFs.
        Its pointwise value/gradient eval is now RELIABLE: the AxiHenrotte
        FESpace is complex-capable (iscomplex from the ``complex`` flag) and the
        DiffOps implement the complex ``CalcMatrix`` overloads (src/ext/axifemm),
        validated to machine precision in tests/test_axi_henrotte_complex_eval.py
        -- so |A|/B post-processing of this complex solution works (this is what
        unblocks the nonlinear mu(|B|) Picard layer).  Requires the rebuilt
        radia_axifemm extension.  ``P_eddy`` remains the matrix-based reference.

    Validated (tests/test_axi_eddy_harmonic.py): for a Cu disk in a uniform
    applied B0, the eddy loss P_eddy(w) is positive, rises with w, and is
    w^2-suppressed at low frequency (the eddy onset) -- and the raw conductor-DOF
    solution shows flux expelled (|A| -> 0) as w rises past 1/(2 pi tau_1),
    tau_1 = 224 us.  The K, M are the SAME integrators behind the validated
    Cu-disk tau_1 eigenvalue, now in a forced harmonic solve.
    """
    import numpy as np
    import scipy.sparse as sp
    import scipy.sparse.linalg as spla
    from radia.radia_axifemm import (H1Henrotte, AxiHenrotteStiffnessBFI,
                                     AxiHenrotteSigmaMassBFI)

    fes = H1Henrotte(mesh, order=order, dirichlet=dirichlet)
    aK = BilinearForm(fes, symmetric=True)
    aK += AxiHenrotteStiffnessBFI(mu_cf)
    aK.Assemble()
    aM = BilinearForm(fes, symmetric=True)
    aM += AxiHenrotteSigmaMassBFI(sigma_cf)
    aM.Assemble()
    bf = LinearForm(fes)
    bf += sigma_cf * applied_A * fes.TestFunction() * dx
    bf.Assemble()

    n = fes.ndof

    def _coo(mat):
        rr, cc, vv = mat.COO()
        return sp.csr_matrix((np.asarray(vv, dtype=float),
                              (np.asarray(rr, dtype=int), np.asarray(cc, dtype=int))),
                             shape=(n, n))

    free = np.array([i for i in range(n) if fes.FreeDofs()[i]], dtype=int)
    K = _coo(aK.mat)[free[:, None], free[None, :]]; K = (K + K.T) * 0.5
    M = _coo(aM.mat)[free[:, None], free[None, :]]; M = (M + M.T) * 0.5
    bv = np.asarray(bf.vec, dtype=float)[free].astype(complex)
    xf = spla.spsolve((K + 1j * omega * M).tocsc(), bv)
    P_eddy = 0.5 * omega * omega * float(np.real(np.conj(xf) @ (M @ xf)))

    fes_c = H1Henrotte(mesh, order=order, complex=True, dirichlet=dirichlet)
    gfu = GridFunction(fes_c)
    arr = gfu.vec.FV().NumPy()
    arr[:] = 0.0
    arr[free] = xf
    return gfu, P_eddy


def solve_axi_eddy_harmonic_nonlinear(mesh, mu_of_B, sigma_cf, omega, applied_A,
                                      order=1, dirichlet="axis|outer",
                                      relax=0.5, max_iter=60, tol=1e-4):
    """NONLINEAR (saturating mu(|B|)) forced time-harmonic axisymmetric eddy
    currents -- a Picard fixed point on top of :func:`solve_axi_eddy_harmonic`.

    Each sweep reassembles the closed-form stiffness ``K = AxiHenrotteStiffnessBFI(mu)``
    with the permeability updated from the CURRENT complex A_phi solution, keeps the
    sigma-mass ``M = AxiHenrotteSigmaMassBFI(sigma)`` fixed, and re-solves the complex
    ``S = K + j*w*M`` in scipy.  This is the FEMM "effective-mu" harmonic-balance
    approximation for nonlinear AC -- the axisymmetric analog of
    :func:`solve_planar_eddy_nonlinear`.

    ``mu_of_B`` : callable(Bmag_cf) -> mu_cf, taking the SCALAR amplitude |B| and
    returning a PERMEABILITY CF.  NOTE the convention is PERMEABILITY mu (what
    AxiHenrotteStiffnessBFI expects) -- the RECIPROCAL of the planar solver's
    ``nu_of_B`` reluctivity.  The amplitude is built from the complex A_phi
    (``B_r = -dA/dz``, ``B_z = dA/dr + A/r``, ``|B| = sqrt(|B_r|^2 + |B_z|^2)``)::

        def mu_of_B(Bmag):                       # Frohlich/Kennelly soft saturation
            return MU0 * (1.0 + (mu_r - 1.0) / (1.0 + (Bmag / Bsat) ** 2))

    This RELIES on the complex H1Henrotte field eval (the src/ext/axifemm complex
    CalcMatrix fix): the BFI samples ``mu_of_B(|B(A)|)`` at element centroids, which
    evaluates the complex A_phi value AND gradient there.  ``relax`` under-relaxes
    the DOF vector (0.3-0.5 for hard saturation).  ORDER=1 ONLY (P2 AxiHenrotte has
    an axis-singularity NaN).  Reduces EXACTLY to :func:`solve_axi_eddy_harmonic`
    when ``mu_of_B`` is constant.  Returns ``(gfu, P_eddy)`` with
    ``P_eddy = 0.5 w^2 Re(x^H M x)`` the reliable matrix-based eddy loss [W].
    """
    import numpy as np
    import scipy.sparse as sp
    import scipy.sparse.linalg as spla
    from radia.radia_axifemm import (H1Henrotte, AxiHenrotteStiffnessBFI,
                                     AxiHenrotteSigmaMassBFI)

    fes = H1Henrotte(mesh, order=order, complex=True, dirichlet=dirichlet)
    n = fes.ndof
    free = np.array([i for i in range(n) if fes.FreeDofs()[i]], dtype=int)

    def _coo(mat):
        rr, cc, vv = mat.COO()
        return sp.csr_matrix((np.asarray(vv, dtype=float),
                              (np.asarray(rr, dtype=int), np.asarray(cc, dtype=int))),
                             shape=(n, n))

    # sigma-mass M and the RHS are FIXED across Picard sweeps.
    aM = BilinearForm(fes, symmetric=True)
    aM += AxiHenrotteSigmaMassBFI(sigma_cf)
    aM.Assemble()
    M = _coo(aM.mat)[free[:, None], free[None, :]]; M = (M + M.T) * 0.5

    bf = LinearForm(fes)
    bf += sigma_cf * applied_A * fes.TestFunction() * dx
    bf.Assemble()
    bv = np.asarray(bf.vec, dtype=float)[free].astype(complex)

    gfu = GridFunction(fes)
    arr = gfu.vec.FV().NumPy()
    arr[:] = 0.0

    # |B| amplitude CF from the current complex A_phi (gfu holds the iterate; the
    # CF re-evaluates with gfu's updated values each sweep).  x == r in axisym.
    dA = grad(gfu)
    Br = -dA[1]
    Bz = dA[0] + gfu / x
    Bmag = sqrt((Br * Conj(Br) + Bz * Conj(Bz)).real + 1e-30)
    mu_cf = mu_of_B(Bmag)

    xf = np.zeros(free.shape[0], dtype=complex)
    sweeps = 0
    for it in range(max_iter):
        aK = BilinearForm(fes, symmetric=True)
        aK += AxiHenrotteStiffnessBFI(mu_cf)          # mu sampled from |B(gfu)| at centroids
        aK.Assemble()
        K = _coo(aK.mat)[free[:, None], free[None, :]]; K = (K + K.T) * 0.5
        xnew = spla.spsolve((K + 1j * omega * M).tocsc(), bv)
        dnorm = float(np.linalg.norm(xnew - xf) / (np.linalg.norm(xnew) + 1e-30))
        xf = xnew if it == 0 else (1.0 - relax) * xf + relax * xnew
        arr[:] = 0.0
        arr[free] = xf                                # update gfu for the next |B|
        sweeps = it + 1
        if it > 0 and dnorm < tol:
            break

    P_eddy = 0.5 * omega * omega * float(np.real(np.conj(xf) @ (M @ xf)))
    return gfu, P_eddy


def solve_magnetostatic_Aform(mesh, nu, source=None, curl_source=None, order=2,
                              dirichlet="outer", reg=1e-6):
    """Solve the A-formulation magnetostatic problem  curl(nu curl A) = J.

    Parameters
    ----------
    mesh     : NGSolve Mesh (outer boundary named ``dirichlet`` for n x A = 0).
    nu       : reluctivity CoefficientFunction (e.g. from ``reluctivity``).
    source   : current-density vector CF (RHS int J.v) or None (no free current).
    curl_source : vector CF entering as int (curl_source . curl v) -- the
               permanent-magnet term int nu Br . curl v (Br = remanent flux
               density). Build it with ``remanent_source`` (mu_r=1 hard magnets)
               or as nu*Br yourself; None = no magnet.
    order    : HCurl order (default 2).
    reg      : gradient-nullspace mass regularization coefficient (x NU0); B is
               gauge-invariant so this only fixes A, not the physical field.

    Returns the HCurl GridFunction ``gfu``; the flux density is ``B = curl(gfu)``.
    """
    fes = HCurl(mesh, order=order, dirichlet=dirichlet)
    u, v = fes.TnT()
    a = BilinearForm(fes)
    a += nu * curl(u) * curl(v) * dx + reg * NU0 * u * v * dx
    f = LinearForm(fes)
    if source is not None:
        f += InnerProduct(source, v) * dx
    if curl_source is not None:
        f += InnerProduct(curl_source, curl(v)) * dx
    a.Assemble()
    f.Assemble()
    gfu = GridFunction(fes)
    gfu.vec.data = _direct_inverse(a.mat, fes.FreeDofs(), "sparsecholesky",
                                   "solve_magnetostatic_Aform") * f.vec
    return gfu


def solve_scattered_uniform_field(mesh, mu_r_by_material, B0_vec, order=2,
                                  dirichlet="outer", reg=1e-8):
    """Scattered-field A-form magnetostatics for permeable bodies in a UNIFORM
    applied field ``B0`` (no coil source). The source is the permeability jump::

        curl(nu curl A~) = -curl(nu B0) = curl(nu_pert B0),
        nu_pert = nu0 (1 - 1/mu_r)   inside each permeable body, 0 elsewhere

    and the TOTAL field is ``B = curl(A~) + B0``. This is the validated machinery
    behind the permeable-sphere / shell cross-checks (demagnetisation N=1/3 sphere
    to <6 %, spherical-shell shielding factor to ~1 %).

    ``mu_r_by_material`` maps material name -> relative permeability (default 1 =
    vacuum for any other region). ``B0_vec`` is a 3-tuple or CoefficientFunction.
    Returns the SCATTERED GridFunction ``A~`` (read total B as ``curl(A~)+B0``).
    """
    B0cf = B0_vec if isinstance(B0_vec, CoefficientFunction) else CoefficientFunction(tuple(B0_vec))
    nu_cf = mesh.MaterialCF({m: NU0 / mr for m, mr in mu_r_by_material.items()}, default=NU0)
    nu_pert = mesh.MaterialCF({m: NU0 * (1.0 - 1.0 / mr) for m, mr in mu_r_by_material.items()},
                              default=0.0)
    fes = HCurl(mesh, order=order, dirichlet=dirichlet)
    u, v = fes.TnT()
    a = BilinearForm(fes)
    a += nu_cf * InnerProduct(curl(u), curl(v)) * dx + reg * NU0 * InnerProduct(u, v) * dx
    f = LinearForm(fes)
    f += InnerProduct(nu_pert * B0cf, curl(v)) * dx
    a.Assemble(); f.Assemble()
    gfA = GridFunction(fes)
    gfA.vec.data = _direct_inverse(a.mat, fes.FreeDofs(), "sparsecholesky",
                                   "solve_scattered_uniform_field") * f.vec
    return gfA


def shell_shielding_factor(mu_r, a, b, geometry="sphere"):
    """Exact magnetic shielding factor S = B_applied / B_cavity for a permeable
    shell (inner radius ``a``, outer ``b``, relative permeability ``mu_r``) in a
    uniform applied field -- the COMSOL AC/DC "magnetic shielding" benchmark.

        sphere   : S = [(2 mu_r+1)(mu_r+2) - 2 (a/b)^3 (mu_r-1)^2] / (9 mu_r)
        cylinder : S = [(mu_r+1)^2 - (a/b)^2 (mu_r-1)^2] / (4 mu_r)   (transverse)

    Both -> 1 at mu_r=1 (no shielding) and -> (k mu_r)(1-(a/b)^d) for large mu_r
    (thin-shell limit, d=3/2 sphere, k=2/9 and 1/4 resp.). Larger S = better
    shielding (smaller interior field).
    """
    r = a / b
    if geometry == "sphere":
        return ((2 * mu_r + 1) * (mu_r + 2) - 2 * r**3 * (mu_r - 1)**2) / (9 * mu_r)
    if geometry == "cylinder":
        return ((mu_r + 1)**2 - r**2 * (mu_r - 1)**2) / (4 * mu_r)
    raise ValueError("geometry must be 'sphere' or 'cylinder'")


def solve_magnetostatic_nonlinear(mesh, nu_of_B, source=None, order=2,
                                  dirichlet="outer", reg=1e-2, relax=0.35,
                                  max_iter=80, tol=1e-4, min_iter=4, monitor=None):
    """Picard fixed-point solve of the A-formulation with a FIELD-DEPENDENT
    reluctivity:  curl(nu(|B|) curl A) = J  (nonlinear / saturating B-H iron).

    Cross-validated against an independent reference (Newton) on the
    saturating-iron force case: NGSolve F_z = -4.898 N vs reference -4.930 N
    (0.65 %), B_iron ~1.138 T; see the ``force_validation`` MCP tool
    (cross_validation case 3). The fixed point matches the reference's
    Newton even though we iterate Picard, because both converge the same
    self-consistent nu(|B|).

    Parameters
    ----------
    nu_of_B : callable(B_cf) -> nu_cf
        Given the current flux density CF ``B = curl(A)``, return the reluctivity
        CF for that field, e.g.::

            def nu_of_B(B):
                Bmag = sqrt(InnerProduct(B, B) + 1e-12)
                mur  = 1.0 + (mur0 - 1.0) / (1.0 + (Bmag / Bk)**p)   # saturation
                return iron_ind / (MU0 * mur) + (1.0 - iron_ind) * NU0

    source  : current-density vector CF (RHS int J.v) or None.
    reg     : gradient-nullspace mass regularization (x NU0). Default 1e-2 --
              larger than the linear solve's 1e-6 to keep the under-relaxed fixed
              point well conditioned; B is gauge-invariant so this only fixes A.
    relax   : under-relaxation factor w for A <- (1-w) A + w A_new. The fixed
              point is independent of ``relax`` (it only affects stability /
              convergence rate), BUT a too-large w does not just slow things --
              for a saturating B-H it DIVERGES: overshooting B above the fixed
              point pushes the iron into the low-mu_r / stiff regime and the
              iteration oscillates and blows up (measured: w=0.85 diverged,
              1.14->1.30->...). Approach from below: keep w in ~0.3-0.5 (default
              0.35 is stable and monotone, ~18-20 sweeps to tol 1e-4).
    max_iter, tol : stop when the relative change of ``monitor`` is < ``tol``.
    min_iter : do at least this many Picard sweeps before the convergence test
               can fire -- guards against a premature stop when the probe starts
               near zero (e.g. an iron-region probe before the body magnetizes,
               where ``|cur-prev| < tol*cur`` is trivially true at cur~0).
    monitor : callable(B_cf) -> float, the scalar convergence probe. Default
              None = global field energy int |B|^2 dV. For a SMALL nonlinear body
              in a large air domain pass a probe localized to that body (e.g. the
              mean |B| over its material) -- the global energy is dominated by the
              source field and can mask the body's (force-determining) convergence.

    Returns the converged HCurl GridFunction ``gfu`` (B = curl(gfu)).
    """
    fes = HCurl(mesh, order=order, dirichlet=dirichlet)
    u, v = fes.TnT()
    f = LinearForm(fes)
    if source is not None:
        f += InnerProduct(source, v) * dx
    f.Assemble()
    gfu = GridFunction(fes)
    gfu.vec[:] = 0.0
    probe = monitor if monitor is not None else (
        lambda B: Integrate(InnerProduct(B, B) * dx, mesh))
    prev = None
    for it in range(max_iter):
        a = BilinearForm(fes)
        a += nu_of_B(curl(gfu)) * curl(u) * curl(v) * dx + reg * NU0 * u * v * dx
        a.Assemble()
        gnew = GridFunction(fes)
        gnew.vec.data = _direct_inverse(a.mat, fes.FreeDofs(), "sparsecholesky",
                                        "solve_magnetostatic_nonlinear") * f.vec
        gfu.vec.data = (1.0 - relax) * gfu.vec + relax * gnew.vec
        cur = probe(curl(gfu))
        if (prev is not None and it + 1 >= min_iter
                and abs(cur - prev) < tol * max(abs(cur), 1e-30)):
            break
        prev = cur
    return gfu


def solve_magnetostatic_newton(mesh, source, energy_density, steel_region,
                               order=2, dirichlet="outer", reg=1e-6,
                               max_iter=30, tol=1e-7, load_steps=1,
                               precond="bddc", cg_tol=1e-8, cg_maxiter=4000,
                               verbose=False):
    """Newton energy-minimisation A-form solve with an ITERATIVE CG + BDDC linear
    solver -- the LARGE-3D nonlinear path (the direct factorisation in
    ``solve_magnetostatic_Aform`` / the Picard ``solve_magnetostatic_nonlinear``
    do not scale to 10^5+ DOF saturating problems like TEAM-20).

    Minimises the magnetic co-energy
        E(A) = int_steel phi(|curl A|) dx + int_else 1/2 nu0 |curl A|^2 dx
             + 1/2 reg nu0 int |A|^2 dx  -  int source . A dx
    by Newton with an Armijo line search; each step solves the linearisation with
    CG preconditioned by ``precond`` (default "bddc"). Robust where the Picard
    fixed point needs a tiny relax / hundreds of sweeps (hard saturation).

    Parameters
    ----------
    source         : current-density vector CF (the int J.v excitation).
    energy_density : callable(Bmag_cf) -> co-energy density phi(|B|)=int_0^|B| H db
                     [J/m^3] of the steel, e.g.
                     ``BSpline(2, [0]+B_list, H_list).Integrate()``.
    steel_region   : materials string for the nonlinear steel (e.g. "yoke|pole").
                     Everywhere else gets the linear vacuum energy 1/2 nu0 |B|^2.
    reg            : |A|^2 gauge regularisation (x nu0) fixing the curl nullspace.
    load_steps     : ramp the source 1/n .. 1 over n steps (continuation) so Newton
                     stays in its basin at high excitation (saturated steel).
    precond        : NGSolve Preconditioner type ("bddc" default; "local"/"multigrid").

    Returns the HCurl GridFunction ``gfu`` (A); flux density ``B = curl(gfu)``.
    """
    from ngsolve.krylovspace import CGSolver
    fes = HCurl(mesh, order=order, dirichlet=dirichlet)
    u, v = fes.TnT()
    gfu = GridFunction(fes)
    gfu.vec[:] = 0.0

    cc = InnerProduct(curl(u), curl(u))
    Bmag = sqrt(cc + 1e-12)
    a = BilinearForm(fes, symmetric=True)
    a += Variation(0.5 * NU0 * cc * dx)                                   # linear everywhere
    a += Variation((energy_density(Bmag) - 0.5 * NU0 * cc)
                   * dx(definedon=mesh.Materials(steel_region)))          # steel correction
    a += Variation(0.5 * reg * NU0 * InnerProduct(u, u) * dx)            # gauge fix
    pre = Preconditioner(a, type=precond)

    f0 = LinearForm(fes)
    f0 += InnerProduct(source, v) * dx
    f0.Assemble()

    r = gfu.vec.CreateVector(); w = gfu.vec.CreateVector()
    au = gfu.vec.CreateVector(); xn = gfu.vec.CreateVector(); fv = gfu.vec.CreateVector()

    for step in range(1, load_steps + 1):
        fv.data = (step / load_steps) * f0.vec
        for it in range(max_iter):
            a.AssembleLinearization(gfu.vec)
            pre.Update()
            a.Apply(gfu.vec, au)
            r.data = fv - au
            inv = CGSolver(mat=a.mat, pre=pre.mat, tol=cg_tol,
                           maxiter=cg_maxiter, printrates=False)
            w.data = inv * r
            err = InnerProduct(w, r)
            if abs(err) < tol:
                if verbose:
                    print(f"  [load {step}/{load_steps}] converged it{it+1} err={err:.2e}", flush=True)
                break
            E0 = a.Energy(gfu.vec) - InnerProduct(fv, gfu.vec)
            tau = 1.0
            xn.data = gfu.vec + w
            while (a.Energy(xn) - InnerProduct(fv, xn)) > E0 - 1e-4 * tau * abs(err) and tau > 1e-10:
                tau *= 0.5
                xn.data = gfu.vec + tau * w
            gfu.vec.data = xn
            if verbose:
                print(f"  [load {step}/{load_steps}] it{it+1} err={err:.2e} tau={tau:.3g}", flush=True)
    return gfu


def solve_eddy_current_harmonic_APhi(
    mesh, nu, sigma, omega, add_source,
    order=5, precond="local", dirichlet="", dirichlet_bbbnd="GND",
    periodic=True,
    tol=1e-8, maxsteps=2000, restart=200, reg=1e-10,
):
    """3-D frequency-domain eddy current via A-Φ formulation.

    Governing equations (weak form, e^{jωt} convention):
      ∫ν curl A·curl v dx  +  jωσ(A+∇Φ)·v dx_cond  =  F_source   [A eq.]
      jωσ(A+∇Φ)·∇ψ dx_cond  =  0                                  [Φ eq.]

    Two open-boundary modes:

    ``periodic=True`` (default, Kelvin transform):
      The mesh carries an outer domain whose dome faces are identified with
      ``face.Identify(..., IdentificationType.PERIODIC)``.  ``nu`` must encode
      the Kelvin scaling (ν_outer = (r'/R)² ν₀).  A GND point vertex fixes the
      Φ gauge (``dirichlet_bbbnd``).  TEAM-7 uses this mode.

    ``periodic=False`` (Dirichlet box):
      The outer boundary carries a Dirichlet tag given by ``dirichlet``.  No
      Kelvin transform and no GND vertex are needed; ``dirichlet_bbbnd`` is
      ignored.  TEAM-21 uses this mode.

    ``sigma`` is zero outside conductors — the bilinear form integrates over all
    of dx, which is equivalent to restricting to conductor regions.

    Parameters
    ----------
    mesh     : NGSolve Mesh.
    nu       : reluctivity CF [m/H].
    sigma    : conductivity CF [S/m] (0 outside conductors).
    omega    : angular frequency [rad/s].
    add_source : ``add_source(lf, (v_A, psi))`` appends source terms to LinearForm
                 ``lf``.  ``v_A`` is the HCurl test function, ``psi`` the H1 test
                 function.  Example for a stranded rectangular coil::

                     def add_source(lf, test_fn):
                         v_A, _ = test_fn
                         lf += -N/L * tau_coil * v_A.Trace() * ds("coili")
                         lf += N/(dw*L) * pot_coil * curl(v_A) * dx("coil")

    order    : HCurl / H1 polynomial order (default 5; order=3 is ~4-5x faster).
    precond  : Preconditioner type for GMRes.  "local" (block Jacobi) is required
               for complex=True — "bddc" triggers a Cholesky complex-number failure.
    dirichlet       : Dirichlet boundary tag for HCurl and (when periodic=False)
                      H1 (default "" = none).
    dirichlet_bbbnd : H1 point tag for Φ gauge fixing (default "GND").
                      Only used when ``periodic=True``.
    periodic : If True (default), wrap FE spaces with Periodic() for Kelvin BC.
               If False, use plain Dirichlet box (``dirichlet`` tag on outer faces).
    tol, maxsteps, restart : GMRes convergence parameters.
    reg      : Φ regularisation coefficient (×ν) preventing a singular kernel
               outside the conductor; default 1e-10.

    Returns
    -------
    gfAPhi : compound GridFunction with components ``(gfA, gfPhi)``.
             ``B = curl(gfAPhi.components[0])``
             ``J = -1j*omega*sigma*(gfAPhi.components[0]
                                   + grad(gfAPhi.components[1]))``
    """
    from ngsolve import Periodic, H1, HCurl, BilinearForm, LinearForm, \
        GridFunction, Preconditioner, TaskManager, curl, grad, dx
    from ngsolve.krylovspace import GMRes

    fes_A_base  = HCurl(mesh, order=order, complex=True,
                        dirichlet=dirichlet, nograds=True)
    fes_A       = Periodic(fes_A_base) if periodic else fes_A_base

    if periodic:
        fes_Ph_base = H1(mesh, order=order, complex=True,
                         dirichlet=dirichlet, dirichlet_bbbnd=dirichlet_bbbnd)
    else:
        fes_Ph_base = H1(mesh, order=order, complex=True,
                         dirichlet=dirichlet)
    fes_Phi     = Periodic(fes_Ph_base) if periodic else fes_Ph_base

    fes         = fes_A * fes_Phi
    (A_h, Ph_h) = fes.TrialFunction()
    (v,   psi)  = fes.TestFunction()

    a = BilinearForm(fes)
    a += nu           * curl(A_h)  * curl(v)     * dx
    a += 1j * omega   * sigma      * A_h         * v         * dx
    a += 1j * omega   * sigma      * grad(Ph_h)  * v         * dx
    a += 1j * omega   * sigma      * A_h         * grad(psi) * dx
    a += 1j * omega   * sigma      * grad(Ph_h)  * grad(psi) * dx
    a += reg * nu     * grad(Ph_h) * grad(psi)   * dx

    pre = Preconditioner(a, precond)

    f = LinearForm(fes)
    add_source(f, (v, psi))

    gf = GridFunction(fes)
    gf.vec[:] = 0

    with TaskManager():
        a.Assemble()
        pre.Update()
        f.Assemble()

    GMRes(a.mat, f.vec, pre=pre.mat, x=gf.vec,
          tol=tol, maxsteps=maxsteps, restart=restart, printrates=False)

    return gf
