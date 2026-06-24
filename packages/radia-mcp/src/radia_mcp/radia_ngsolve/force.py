"""Executable electromagnetic force / energy extractors for the radia-ngsolve
A-formulation FEM path.

These are the COMSOL-cross-validated methods, as RUNNABLE code (not just the
theory in ``differential_forms``). Each function takes ``B`` -- the magnetic
flux density CoefficientFunction, ``B = curl(gfu)`` for an HCurl GridFunction
``gfu`` -- plus the NGSolve mesh, and returns an SI quantity (force [N],
energy [J], inductance [H]).

Validated against COMSOL (LiveLink) and analytics; see the ``force_validation``
MCP tool for the agreement table (sphere 0.11 %, coil+iron force ~3 %,
self-inductance 0.01 %, ...). The regression tests in
``validation/force/validate_force_xval.py`` assert these keep matching.

#25 lesson baked in: for a HIGH-permeability body do NOT carve a separate nested
"shell" material around it -- the nested-sphere interface isolates the body and
zeroes its interior B. Put the body directly in the surrounding air;
``eggshell_force`` integrates a radial weight band over that plain air, so no
extra material region is required.
"""
import math

from ngsolve import (CoefficientFunction, InnerProduct, sqrt, dx, ds, Integrate,
                     IfPos, specialcf, Conj, x, y, z)

MU0 = 4.0e-7 * math.pi
EPS0 = 8.8541878128e-12


def _float_vector(values, name):
    vec = [float(value) for value in values]
    if len(vec) not in (2, 3):
        raise ValueError(f"{name} must have length 2 or 3")
    return vec


def _unit_vector(values, name):
    vec = _float_vector(values, name)
    norm = math.sqrt(sum(value * value for value in vec))
    if norm <= 0.0:
        raise ValueError(f"{name} must be nonzero")
    return [value / norm for value in vec]


def maxwell_stress_tensor_air(B, mu=MU0):
    """Pointwise magnetic Maxwell stress tensor in air.

    ``B`` is a 2- or 3-component flux-density vector [T].  The returned nested
    list is

        T_ij = (B_i B_j - 0.5 |B|^2 delta_ij) / mu

    in pascals.  This dependency-free helper mirrors the integrand used by the
    surface and weighted-stress FEM extractors, so examples can teach the local
    traction identity before moving to mesh integrals.
    """

    mu = float(mu)
    if mu <= 0.0:
        raise ValueError("mu must be > 0")
    b = _float_vector(B, "B")
    b2 = sum(value * value for value in b)
    dim = len(b)
    return [
        [
            (b[i] * b[j] - (0.5 * b2 if i == j else 0.0)) / mu
            for j in range(dim)
        ]
        for i in range(dim)
    ]


def maxwell_traction_air(B, normal, mu=MU0):
    """Maxwell traction vector ``T n`` in air for a unit surface normal.

    ``normal`` is normalised internally; it must have the same length as ``B``.
    For a uniform normal field this returns ``p n`` with
    ``p = B^2/(2 mu)``, which is exactly :func:`air_gap_maxwell_pressure`.
    A purely tangential field gives ``-p n`` (magnetic tension).
    """

    b = _float_vector(B, "B")
    n = _unit_vector(normal, "normal")
    if len(b) != len(n):
        raise ValueError("B and normal must have the same length")
    tensor = maxwell_stress_tensor_air(b, mu=mu)
    return [
        sum(tensor[i][j] * n[j] for j in range(len(n)))
        for i in range(len(n))
    ]


def maxwell_traction_summary(B, normal, area_m2=1.0, mu=MU0):
    """JSON-friendly Maxwell traction decomposition for one surface patch.

    The normal component is

        traction . n = (B_n^2 - |B_t|^2) / (2 mu)

    and the tangential component has magnitude ``|B_n B_t| / mu``.  ``area_m2``
    scales the traction to a force vector for simple patch/air-gap examples.
    """

    area = float(area_m2)
    if area < 0.0:
        raise ValueError("area_m2 must be >= 0")
    mu = float(mu)
    if mu <= 0.0:
        raise ValueError("mu must be > 0")
    b = _float_vector(B, "B")
    n = _unit_vector(normal, "normal")
    if len(b) != len(n):
        raise ValueError("B and normal must have the same length")
    traction = maxwell_traction_air(b, n, mu=mu)
    b_normal = sum(bi * ni for bi, ni in zip(b, n))
    b2 = sum(bi * bi for bi in b)
    b_tangent2 = max(0.0, b2 - b_normal * b_normal)
    normal_traction = sum(ti * ni for ti, ni in zip(traction, n))
    tangential_traction = [
        ti - normal_traction * ni
        for ti, ni in zip(traction, n)
    ]
    return {
        "B": b,
        "normal": n,
        "mu": mu,
        "area_m2": area,
        "B_normal_T": b_normal,
        "B_tangent_T": math.sqrt(b_tangent2),
        "traction_Pa": traction,
        "normal_traction_Pa": normal_traction,
        "normal_traction_identity_Pa": (b_normal * b_normal - b_tangent2) / (2.0 * mu),
        "tangential_traction_Pa": tangential_traction,
        "tangential_traction_magnitude_Pa": math.sqrt(
            sum(value * value for value in tangential_traction)
        ),
        "force_N": [area * value for value in traction],
    }


def air_gap_maxwell_pressure(B_T, mu=MU0):
    """Magnetic pressure [Pa] for a normal flux density in an air gap.

    For a locally uniform normal field at an iron/air interface, the Maxwell
    stress gives

        p = B^2 / (2 mu)

    with ``mu=mu0`` for air.  The same value is the magnetic energy density in
    the gap.  This tiny helper is intentionally dependency-free so magnetic
    circuit examples can turn a solved ``B_T`` directly into a holding-force
    estimate before running a full weighted-stress FEM extraction.
    """

    B = float(B_T)
    mu = float(mu)
    if mu <= 0.0:
        raise ValueError("mu must be > 0")
    return B * B / (2.0 * mu)


def air_gap_holding_force(B_T, area_m2, faces=1, mu=MU0):
    """Uniform-gap holding force [N] from flux density and active pole area.

    ``faces`` is the number of active, equal pole faces/gaps contributing the
    same pressure.  Use ``faces=2`` for a symmetric two-pole yoke with two equal
    gaps; keep ``faces=1`` for a single plunger or one pole face.
    """

    area = float(area_m2)
    faces = int(faces)
    if area < 0.0:
        raise ValueError("area_m2 must be >= 0")
    if faces < 1:
        raise ValueError("faces must be >= 1")
    return air_gap_maxwell_pressure(B_T, mu=mu) * area * faces


def air_gap_force_summary(B_T, area_m2, faces=1, mu=MU0):
    """Readable JSON-friendly air-gap force summary."""

    pressure = air_gap_maxwell_pressure(B_T, mu=mu)
    area = float(area_m2)
    faces = int(faces)
    force = air_gap_holding_force(B_T, area, faces=faces, mu=mu)
    return {
        "B_T": float(B_T),
        "mu": float(mu),
        "area_m2": area,
        "faces": faces,
        "pressure_Pa": pressure,
        "energy_density_J_per_m3": pressure,
        "force_N": force,
        "force_per_area_N_per_m2": force / (area * faces) if area > 0.0 else math.inf,
    }


def air_gap_shear_stress(B_radial_T, B_tangential_T, mu=MU0):
    """Tangential Maxwell shear stress [Pa] in a cylindrical air gap.

    For radial and tangential flux-density components ``Br`` and ``Bt`` on a
    cylindrical integration surface, the tangential traction is

        tau = Br Bt / mu

    with sign set by ``Bt``.  This is the local stress used by air-gap motor
    torque estimates and by FE Maxwell-stress post-processing.
    """

    mu = float(mu)
    if mu <= 0.0:
        raise ValueError("mu must be > 0")
    return float(B_radial_T) * float(B_tangential_T) / mu


def air_gap_shear_torque(
    B_radial_T,
    B_tangential_T,
    radius_m,
    axial_length_m=1.0,
    angle_rad=2.0 * math.pi,
    mu=MU0,
):
    """Torque [N.m] from uniform air-gap Maxwell shear stress.

    For a cylindrical surface patch, ``area = radius * angle * axial_length``
    and ``torque = radius * tau * area``.  Use ``angle_rad=2*pi`` for a full
    machine, or a sector angle for a symmetry-sector result before multiplying
    by the sector count.
    """

    radius = float(radius_m)
    length = float(axial_length_m)
    angle = float(angle_rad)
    if radius < 0.0:
        raise ValueError("radius_m must be >= 0")
    if length < 0.0:
        raise ValueError("axial_length_m must be >= 0")
    if angle < 0.0:
        raise ValueError("angle_rad must be >= 0")
    shear = air_gap_shear_stress(B_radial_T, B_tangential_T, mu=mu)
    return shear * radius * radius * angle * length


def air_gap_shear_torque_summary(
    B_radial_T,
    B_tangential_T,
    radius_m,
    axial_length_m=1.0,
    angle_rad=2.0 * math.pi,
    mu=MU0,
):
    """JSON-friendly air-gap shear-stress torque summary."""

    radius = float(radius_m)
    length = float(axial_length_m)
    angle = float(angle_rad)
    if radius < 0.0:
        raise ValueError("radius_m must be >= 0")
    if length < 0.0:
        raise ValueError("axial_length_m must be >= 0")
    if angle < 0.0:
        raise ValueError("angle_rad must be >= 0")
    shear = air_gap_shear_stress(B_radial_T, B_tangential_T, mu=mu)
    area = radius * angle * length
    force = shear * area
    torque = force * radius
    return {
        "B_radial_T": float(B_radial_T),
        "B_tangential_T": float(B_tangential_T),
        "mu": float(mu),
        "radius_m": radius,
        "axial_length_m": length,
        "angle_rad": angle,
        "surface_area_m2": area,
        "shear_stress_Pa": shear,
        "tangential_force_N": force,
        "torque_Nm": torque,
        "torque_per_axial_length_N": torque / length if length > 0.0 else math.inf,
    }


def electrostatic_eggshell_force(E, mesh, gradg, air_region="air"):
    """Weighted Maxwell-stress ("eggshell") ELECTROSTATIC force -- the electric twin
    of :func:`eggshell_force` (ε0 E in place of B/μ0). ``E`` is the electric field
    CoefficientFunction (``E = -grad(gfV)``); ``gradg`` = grad(g) of a smooth weight g
    (=1 on the body side, 0 on the far side of a band that lies in air), a vector CF
    nonzero only inside the band:

        F_k = - int_air [ eps0 E_k (E.gradg) - (eps0/2) |E|^2 d_k g ] dV

    Returns (Fx, Fy, Fz) in newtons. For a compact body use a spherical band
    (gradg = band * (r-center)/|r-center| * -1/(r_outer-r_inner)); for a plate/gap use
    an axis-aligned ramp (e.g. gradg = (0,0,g'(z)) across the gap).
    """
    Edg = InnerProduct(E, gradg)
    E2 = InnerProduct(E, E)
    region = dx(definedon=mesh.Materials(air_region))
    return tuple(-Integrate((EPS0 * E[k] * Edg - 0.5 * EPS0 * E2 * gradg[k]) * region, mesh)
                 for k in range(3))


def electrostatic_eggshell_force_2d(E, mesh, gradg, air_region="air"):
    """2D weighted Maxwell-stress ("eggshell") ELECTROSTATIC force [N/m] -- the 2D
    twin of :func:`electrostatic_eggshell_force` (range 2, per unit out-of-plane
    depth). ``E`` = ``-grad(V)`` (a 2-vector CF); ``gradg`` = grad of a smooth weight
    ``g`` (=1 on the body side, 0 on the far side across a band lying in the
    dielectric), a 2-vector CF nonzero only inside the band:

        F_k = - int [ eps0 E_k (E.gradg) - (eps0/2) |E|^2 d_k g ] dA .

    For a plate/gap use an axis-aligned ramp ``gradg = (0, g'(y))`` (a horizontal
    band in the gap); for a compact body a radial band. Returns (Fx, Fy) in N/m.

    WHY a volume band, not a boundary trace of ``grad(V)``: on a CONDUCTOR face
    ``V`` is constant, so the boundary trace of ``grad(V)`` (in a ``ds`` integral)
    keeps only the TANGENTIAL derivative (=0) and DROPS the normal field -- the
    surface-stress integral then reads ~0. The volume gradient is the true field, so
    the eggshell band (integrating in the dielectric AROUND the conductor) is the
    correct, robust extractor. (DEAD END recorded in ngsolve_usage("electro_mechanical").)
    """
    Edg = InnerProduct(E, gradg)
    E2 = InnerProduct(E, E)
    region = dx(definedon=mesh.Materials(air_region))
    return tuple(-Integrate((EPS0 * E[k] * Edg - 0.5 * EPS0 * E2 * gradg[k]) * region, mesh)
                 for k in range(2))


def magnetic_energy_2d(B, mesh, region=None):
    """2D magnetic field energy PER UNIT LENGTH [J/m]:  W = int |B|^2/(2 mu0) dA  over
    ``region`` (a material name; whole mesh if None).  For a current-driven problem the
    inductance is ``L = 2 W / I^2`` [H/m].  ``B`` = ``CF((grad(A)[1], -grad(A)[0]))``.

    NOTE: for a 1/r field (e.g. a coaxial line) the |B|^2 integrand is sharply peaked at
    small r, so refine the mesh near the inner radius -- the coax energy converged
    1.8 % -> 0.5 % under such refinement (examples/comsol_class/coax_line.py)."""
    dom = dx if region is None else dx(definedon=mesh.Materials(region))
    return Integrate(InnerProduct(B, B) / (2.0 * MU0) * dom, mesh)


def eggshell_force(B, mesh, center, r_inner, r_outer, air_region="air"):
    """Weighted Maxwell-stress ("eggshell") force on the body inside ``r_inner``.

    Robust on unstructured meshes: replaces the sharp surface integral by a
    volume integral over a radial band ``r_inner < |r-center| < r_outer`` with a
    smooth weight g (=1 at r_inner, 0 at r_outer):

        F_k = - int_band [ (1/mu0) B_k (B.grad g) - (1/2mu0) |B|^2 d_k g ] dV

    The band must lie in the air surrounding the body (choose ``r_inner`` >= the
    body radius). Returns (Fx, Fy, Fz) in newtons.
    """
    cx, cy, cz = center
    rho = sqrt((x - cx)**2 + (y - cy)**2 + (z - cz)**2)
    band = IfPos(rho - r_inner, IfPos(r_outer - rho, 1.0, 0.0), 0.0)
    gscale = -1.0 / (r_outer - r_inner)
    gradg = band * CoefficientFunction((x - cx, y - cy, z - cz)) / rho * gscale
    Bdg = InnerProduct(B, gradg)
    B2 = InnerProduct(B, B)
    region = dx(definedon=mesh.Materials(air_region))
    F = []
    for k in range(3):
        integ = (1.0 / MU0) * B[k] * Bdg - (1.0 / (2.0 * MU0)) * B2 * gradg[k]
        F.append(-Integrate(integ * region, mesh))
    return tuple(F)


def eggshell_torque(B, mesh, center, r_inner, r_outer, pivot=(0.0, 0.0, 0.0),
                    air_region="air"):
    """3D weighted Maxwell-stress ("eggshell") TORQUE [N m] about ``pivot`` on
    the body inside ``r_inner``.  3D analogue of :func:`eggshell_torque_2d`:

        tau = - int_band  r' x S  dV,   r' = r - pivot,
        S_k = (1/mu0) B_k (B.grad g) - (1/2mu0) |B|^2 d_k g,

    same radial weight band (g=1 at r_inner, 0 at r_outer) in the air around the
    body. Returns (Tx, Ty, Tz). Use the same band as :func:`eggshell_force`;
    validated on a magnetised cylinder in a uniform field (tau = m x B0,
    examples/comsol_class/motor_torque.py)."""
    cx, cy, cz = center
    px, py, pz = pivot
    rho = sqrt((x - cx)**2 + (y - cy)**2 + (z - cz)**2)
    band = IfPos(rho - r_inner, IfPos(r_outer - rho, 1.0, 0.0), 0.0)
    gscale = -1.0 / (r_outer - r_inner)
    gradg = band * CoefficientFunction((x - cx, y - cy, z - cz)) / rho * gscale
    Bdg = InnerProduct(B, gradg)
    B2 = InnerProduct(B, B)
    S = [(1.0 / MU0) * B[k] * Bdg - (1.0 / (2.0 * MU0)) * B2 * gradg[k] for k in range(3)]
    rp = (x - px, y - py, z - pz)
    cross = ((rp[1] * S[2] - rp[2] * S[1]),     # r' x S
             (rp[2] * S[0] - rp[0] * S[2]),
             (rp[0] * S[1] - rp[1] * S[0]))
    region = dx(definedon=mesh.Materials(air_region))
    return tuple(-Integrate(c * region, mesh) for c in cross)


def eggshell_force_2d(B, mesh, center, r_inner, r_outer, air_region="air"):
    """2D PLANAR weighted Maxwell-stress ("eggshell") force [N/m] on the body
    inside ``r_inner`` (force PER UNIT LENGTH in the out-of-plane direction).

    Same identity as :func:`eggshell_force` but in the (x, y) plane with the
    2-vector ``B = (Bx, By)`` (e.g. from the A_z solver,
    ``B = CF((grad(gfu)[1], -grad(gfu)[0]))``):

        F_k = - int_band [ (1/mu0) B_k (B.grad g) - (1/2mu0) |B|^2 d_k g ] dA

    The radial band ``r_inner < |r-center| < r_outer`` must lie in the air
    surrounding the body and enclose it. Validated on the two-parallel-wire
    benchmark F/L = mu0 I1 I2 / (2 pi d) to ~1 % (see tests/test_planar_force.py).
    Returns (Fx, Fy) in N/m.
    """
    cx, cy = center
    rho = sqrt((x - cx)**2 + (y - cy)**2)
    band = IfPos(rho - r_inner, IfPos(r_outer - rho, 1.0, 0.0), 0.0)
    gscale = -1.0 / (r_outer - r_inner)
    gradg = band * CoefficientFunction((x - cx, y - cy)) / rho * gscale
    Bdg = InnerProduct(B, gradg)
    B2 = InnerProduct(B, B)
    region = dx(definedon=mesh.Materials(air_region))
    F = []
    for k in range(2):
        integ = (1.0 / MU0) * B[k] * Bdg - (1.0 / (2.0 * MU0)) * B2 * gradg[k]
        F.append(-Integrate(integ * region, mesh))
    return tuple(F)


def eggshell_torque_2d(B, mesh, center, r_inner, r_outer, pivot=(0.0, 0.0),
                       air_region="air"):
    """2D PLANAR weighted Maxwell-stress torque [N] (per unit length) about
    ``pivot``, on the body inside ``r_inner``:

        tau_z = - int_band [ x' S_y - y' S_x ] dA,
        S_k = (1/mu0) B_k (B.grad g) - (1/2mu0)|B|^2 d_k g,   r' = r - pivot.

    Same eggshell band convention as :func:`eggshell_force_2d`. The lever-arm
    weighting is validated to reproduce r' x F of the validated force.
    """
    cx, cy = center
    px, py = pivot
    rho = sqrt((x - cx)**2 + (y - cy)**2)
    band = IfPos(rho - r_inner, IfPos(r_outer - rho, 1.0, 0.0), 0.0)
    gscale = -1.0 / (r_outer - r_inner)
    gradg = band * CoefficientFunction((x - cx, y - cy)) / rho * gscale
    Bdg = InnerProduct(B, gradg)
    B2 = InnerProduct(B, B)
    Sx = (1.0 / MU0) * B[0] * Bdg - (1.0 / (2.0 * MU0)) * B2 * gradg[0]
    Sy = (1.0 / MU0) * B[1] * Bdg - (1.0 / (2.0 * MU0)) * B2 * gradg[1]
    integ = (x - px) * Sy - (y - py) * Sx
    return -Integrate(integ * dx(definedon=mesh.Materials(air_region)), mesh)


def eggshell_force_axi(B, mesh, center, r_inner, r_outer, air_region="air"):
    """AXISYMMETRIC net AXIAL force [N] (full 3D torus) on the body whose
    meridional cross-section lies inside ``r_inner`` of ``center=(rc, zc)``,
    via the weighted Maxwell-stress ("eggshell") band in the (r, z) half-plane.

    ``B`` = ``CF((B_r, B_z))`` is the meridional flux density
    (``B_r = -grad(u)[1]``, ``B_z = grad(u)[0] + u/r`` from
    :func:`solve_axi_magnetostatic`).  The band ``r_inner < rho < r_outer``,
    ``rho = |(r, z) - center|``, must lie in the air enclosing the body's
    cross-section.  Only the AXIAL force is returned -- an axisymmetric body has
    no net radial force.  The ``2*pi*r`` toroidal weight makes this a 3D force:

        F_z = -2 pi int_band [ (1/mu0) B_z (B.grad g)
                               - (1/2mu0)|B|^2 d_z g ] r dr dz,

    g = 1 at r_inner, 0 at r_outer.  Validated on two coaxial loops against the
    exact mutual-inductance force I1 I2 dM/dz (M via elliptic integrals)
    (tests/test_axi_force.py)."""
    rc, zc = center
    rho = sqrt((x - rc)**2 + (y - zc)**2)
    band = IfPos(rho - r_inner, IfPos(r_outer - rho, 1.0, 0.0), 0.0)
    gscale = -1.0 / (r_outer - r_inner)
    gradg = band * CoefficientFunction((x - rc, y - zc)) / rho * gscale
    Bdg = InnerProduct(B, gradg)          # B.grad g
    B2 = InnerProduct(B, B)
    integ_z = (1.0 / MU0) * B[1] * Bdg - (1.0 / (2.0 * MU0)) * B2 * gradg[1]
    return -2.0 * math.pi * Integrate(
        integ_z * x * dx(definedon=mesh.Materials(air_region)), mesh)


def maxwell_surface_force(B, mesh, surface):
    """Maxwell-stress force as a surface integral over a named closed boundary
    ``surface`` in air enclosing the body (outward normal):

        F_k = oint (1/mu0) [ B_k (B.n) - 1/2 |B|^2 n_k ] dS

    Use ``eggshell_force`` instead unless you have a clean meshed surface --
    point/surface traces are noisier than the volume-band method. (Fx, Fy, Fz) [N].
    """
    n = specialcf.normal(mesh.dim)
    Bn = InnerProduct(B, n)
    B2 = InnerProduct(B, B)
    bnd = ds(definedon=mesh.Boundaries(surface))
    F = []
    for k in range(3):
        integ = (1.0 / MU0) * B[k] * Bn - (1.0 / (2.0 * MU0)) * B2 * n[k]
        F.append(Integrate(integ * bnd, mesh))
    return tuple(F)


def maxwell_surface_force_harmonic(B, mesh, surface):
    """TIME-AVERAGED Maxwell-stress force [N] over a closed boundary ``surface``
    for a COMPLEX time-harmonic flux density ``B`` (phasor).  The time-average of
    the quadratic Maxwell stress is

        <F_k> = oint_S [ (1/(2 mu0)) Re( B_k conj(B.n) )
                          - (1/(4 mu0)) |B|^2 n_k ] dS ,

    where the extra 1/2 vs the static :func:`maxwell_surface_force` is the
    time-average of cos^2(wt).  For a field oscillating with REAL peak amplitude
    B0 this equals exactly ``0.5 * maxwell_surface_force(B0)``.

    Intended for the time-harmonic eddy force on a body enclosed by an air
    surface -- in particular the ``"sibc"`` hole boundary of
    ``calc_fem_kelvin`` (the workpiece is a HOLE, so the only force handle is the
    Maxwell stress over its surface; ``B = curl(gfu)`` from the SIBC A-solve).
    Returns (Fx, Fy, Fz) in newtons.  Validated by reduction to the
    COMSOL-cross-validated :func:`maxwell_surface_force` (tests/test_maxwell_surface_harmonic.py).
    """
    n = specialcf.normal(mesh.dim)
    Bn = sum(B[k] * n[k] for k in range(3))                  # B . n  (n real)
    B2 = sum((B[k] * Conj(B[k])).real for k in range(3))     # |B|^2 (real)
    bnd = ds(definedon=mesh.Boundaries(surface))
    F = []
    for k in range(3):
        integ = (0.5 / MU0) * (B[k] * Conj(Bn)).real - (0.25 / MU0) * B2 * n[k]
        F.append(Integrate(integ * bnd, mesh))
    return tuple(F)


def ohmic_loss_2d(Ez, mesh, sigma, region=None):
    """Time-averaged ohmic loss PER UNIT LENGTH [W/m] for a 2D planar harmonic
    eddy problem:  P = 1/2 int sigma |E_z|^2 dA  (E_z = -j w A_z (+ Vc)).

    With a current-driven conductor (net current I), the AC resistance per length
    is ``Rac = 2 P / |I|^2``. Validated on the round-wire skin effect (Rac/Rdc vs
    Kelvin functions, 0.07 %; see tests/test_planar_eddy.py)."""
    integrand = 0.5 * sigma * (Ez * Conj(Ez)).real
    dom = dx if region is None else dx(definedon=mesh.Materials(region))
    return Integrate(integrand * dom, mesh)


def lorentz_force_2d(Jz, B, mesh, region):
    """2D PLANAR Lorentz force PER UNIT LENGTH [N/m] on the conductor ``region``
    carrying out-of-plane current density ``Jz`` [A/m^2] in flux density
    ``B = (Bx, By)``:

        F = int J x B dA = int Jz (zhat x B) dA   ->   Fx = -int Jz By,  Fy = int Jz Bx.

    Pass the TOTAL field B (the conductor's own self-field exerts ZERO net force by
    symmetry, so it drops out). The direct current-source twin of the Maxwell-stress
    :func:`eggshell_force_2d`. Returns ``(Fx, Fy)`` [N/m]. Validated on parallel
    busbars: |F| = mu0 I1 I2/(2 pi d) (examples/comsol_class/busbar_force.py)."""
    dom = dx(definedon=mesh.Materials(region))
    Fx = -Integrate(Jz * B[1] * dom, mesh)
    Fy = Integrate(Jz * B[0] * dom, mesh)
    return Fx, Fy


def magnetic_energy(B, mesh, region=None):
    """Field energy  W = 1/2 integral |B|^2 / mu0 dV  [J].
    ``region=None`` integrates the whole domain; else a material name."""
    integrand = 0.5 * InnerProduct(B, B) / MU0
    if region is None:
        return Integrate(integrand * dx, mesh)
    return Integrate(integrand * dx(definedon=mesh.Materials(region)), mesh)


def self_inductance(B, mesh, i_terminal):
    """Self-inductance  L = 2W / I^2  [H] from the field energy (I = terminal
    current = ampere-turns / N)."""
    return 2.0 * magnetic_energy(B, mesh) / (i_terminal * i_terminal)


def inductance_2d(B, mesh, nu, current, region=None):
    """2D PLANAR inductance PER UNIT LENGTH [H/m] via the energy method
    L = 2W/I^2,  W = 1/2 int nu |B|^2 dA  (``nu`` = reluctivity CF, B in-plane).

    ``region=None`` integrates the whole domain (total L); pass a material name
    for a partial energy (e.g. a conductor's internal inductance). Validated:
    round-wire internal inductance L_int = mu0/(8 pi) = 5.0e-8 H/m, radius-
    independent, to 0.06 % (tests/test_planar_inductance.py)."""
    dom = dx if region is None else dx(definedon=mesh.Materials(region))
    W = 0.5 * Integrate(nu * InnerProduct(B, B) * dom, mesh)
    return 2.0 * W / (current * current)


def inductance_axi(B, mesh, nu, current, region=None):
    """3D INDUCTANCE [H] of an axisymmetric coil via the energy method.

    L = 2 W_3D / I^2,  W_3D = 2*pi * int_half-plane  nu/2 |B|^2 r dr dz

    where the factor 2*pi converts the meridional-plane energy to the full
    toroidal (3D) energy.  B = (B_z, B_r) = (grad(u)[0]+u/r, -grad(u)[1])
    from ``solve_axi_magnetostatic``; ``nu`` = reluctivity; ``current`` = total
    current through the coil cross-section [A].

    ``region=None`` uses the whole domain; pass a material name to restrict."""
    dom = dx if region is None else dx(definedon=mesh.Materials(region))
    W_half = 0.5 * Integrate(nu * InnerProduct(B, B) * x * dom, mesh)
    W_3D = 2.0 * math.pi * W_half
    return 2.0 * W_3D / (current * current)


def ohmic_loss_axi(E_phi, mesh, sigma, region=None):
    """Time-averaged ohmic loss [W] (full 3D torus) for an axisymmetric
    time-harmonic eddy problem:

        P = 2*pi * int_half-plane  sigma/2 |E_phi|^2 r dr dz

    E_phi = -j*omega*A_phi (+ Vc/r for a driven conductor with NumberSpace Vc).
    AC resistance of a ring conductor: Rac = 2*P / |I|^2.

    Analogous to ``ohmic_loss_2d`` but accounts for the toroidal (2*pi*r) volume.
    """
    integrand = 0.5 * sigma * (E_phi * Conj(E_phi)).real * x
    dom = dx if region is None else dx(definedon=mesh.Materials(region))
    return 2.0 * math.pi * Integrate(integrand * dom, mesh)
