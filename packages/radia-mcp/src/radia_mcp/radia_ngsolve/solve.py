"""Reusable NGSolve magnetostatic A-formulation solver -- the validated,
high-mu-safe foundation that the force / energy extractors (``force.py``) and the
COMSOL cross-validation tests build on.

#25 safe: put high-permeability bodies DIRECTLY in the surrounding air (no nested
"shell" material) so the interior B is correct; see the ``force_validation`` MCP
tool. Cross-validated against COMSOL (sphere field 0.11 %, solenoid B <0.35 %,
inductance 0.01 %, ...).
"""
import cmath
import math

from ngsolve import (HCurl, H1, NumberSpace, FESpace, BilinearForm, LinearForm,
                     GridFunction, CoefficientFunction, InnerProduct, curl, grad,
                     dx, Integrate, Conj, Variation, Preconditioner, x, y, sqrt)

MU0 = 4.0e-7 * math.pi
NU0 = 1.0 / MU0


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


def solve_planar_magnetostatic(mesh, nu, Jz=None, magnets=None, order=2,
                               dirichlet="outer"):
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
    a.Assemble()
    f.Assemble()
    gfu = GridFunction(fes)
    gfu.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * f.vec
    return gfu


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
                      dirichlet="outer"):
    """2D PLANAR time-harmonic eddy currents (complex A_z) -- FEMM ``prob2big`` analog.

    Solves  -div(nu grad A_z) + j w sigma A_z = J_src .  The eddy current in
    conductors (sigma>0) is J_eddy = -j w sigma A_z (added automatically). Drive
    modes:

    * ``Jz`` : imposed source current density CF [A/m^2] (e.g. a stranded coil).
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
    BEM-Foster (see packages/radia-axifemm/tests/test_disk_eigenvalue.py).
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
        CAVEAT: H1Henrotte's COMPLEX value-eval falls back to a base-class Id
        stub and is NOT reliable for post-processing |A|/B pointwise; use
        ``P_eddy`` for the loss.  (A correct complex field eval would need the
        AxiHenrotte DiffOp Apply implemented in C++.)

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
    gfu.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * f.vec
    return gfu


def solve_magnetostatic_nonlinear(mesh, nu_of_B, source=None, order=2,
                                  dirichlet="outer", reg=1e-2, relax=0.35,
                                  max_iter=80, tol=1e-4, min_iter=4, monitor=None):
    """Picard fixed-point solve of the A-formulation with a FIELD-DEPENDENT
    reluctivity:  curl(nu(|B|) curl A) = J  (nonlinear / saturating B-H iron).

    Cross-validated against COMSOL (Newton) on the saturating-iron force case:
    NGSolve F_z = -4.842 N vs COMSOL -4.930 N (~1.8 %), B_iron ~1.13 T; see the
    ``force_validation`` MCP tool (comsol_xval case 3). The fixed point matches
    COMSOL's FullyCoupled/Newton even though we iterate Picard, because both
    converge the same self-consistent nu(|B|).

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
        gnew.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * f.vec
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
