"""Regression tests for radia_mcp.radia_ngsolve.force -- the executable EM force
/ energy extractors. They ASSERT that the NGSolve magnetostatic force pipeline
keeps reproducing the values cross-validated against an independent reference
solve (see the ``force_validation`` MCP tool for the agreement table).

Skipped automatically if ngsolve / netgen are not installed. These are heavy FEM
solves (sparse-Cholesky on ~40-50k tets), all marked validation: run
the cross-validation suite with ``python validation/.../validate_*.py``, or skip it with
``the default pytest gate``. The two linear solves take ~1-2 min total; the
nonlinear saturating-iron case runs ~20 Picard sweeps (multi-minute).
"""
import math
import pytest

pytest.importorskip("ngsolve")
pytest.importorskip("netgen")

from ngsolve import (Mesh, HCurl, BilinearForm, LinearForm, GridFunction,
                     CoefficientFunction, InnerProduct, Conj, curl, grad,
                     dx, ds, Integrate,
                     L2, BoundaryFromVolumeCF, specialcf, IfPos, BSpline,
                     x, y, z, sqrt)
from netgen.csg import CSGeometry, Sphere, OrthoBrick, Cylinder, Plane, Pnt, Vec
from netgen.geom2d import SplineGeometry

from radia_mcp.radia_ngsolve.force import (eggshell_force, self_inductance, MU0,
                                           electrostatic_eggshell_force)
from radia_mcp.radia_ngsolve.solve import (
    solve_magnetostatic_Aform, solve_magnetostatic_nonlinear,
    azimuthal_coil_current, remanent_source, solve_magnetostatic_newton,
    solve_eddy_current_harmonic_APhi,
    solve_scattered_uniform_field, shell_shielding_factor,
    solve_planar_magnetostatic, halbach_dipole_magnetization, halbach_bore_field)
from radia_mcp.radia_ngsolve.electrostatic3d import (
    solve_electrostatic_3d, capacitance_from_energy, spherical_capacitor_C,
    capacitance_matrix, EPS0)
from radia_mcp.radia_ngsolve.fieldquality import (
    multipole_coefficients, line_current_multipoles, superpose_multipoles,
    field_errors)
from radia_mcp.radia_ngsolve.scalar_fem2d import (
    solve_electrostatic as solve_electrostatic_2d, capacitance as capacitance_2d)

MM = 1e-3
NU0 = 1.0 / MU0


def _ring(R, dR, zc, hw):
    """Rectangular-section azimuthal coil ring (mean radius R, +/-dR, +/-hw in z)."""
    co = Cylinder(Pnt(0, 0, -1), Pnt(0, 0, 1), R + dR)
    ci = Cylinder(Pnt(0, 0, -1), Pnt(0, 0, 1), R - dR)
    slab = Plane(Pnt(0, 0, zc - hw), Vec(0, 0, -1)) * Plane(Pnt(0, 0, zc + hw), Vec(0, 0, 1))
    return (co - ci) * slab


def _solve_Aform(mesh, nu, J0):
    """A-form magnetostatic solve with an azimuthal coil current of density J0."""
    fes = HCurl(mesh, order=2, dirichlet="outer")
    u, v = fes.TnT()
    coil_ind = CoefficientFunction([1.0 if m == "coil" else 0.0 for m in mesh.GetMaterials()])
    r_xy = sqrt(x * x + y * y + 1e-12)
    Jvec = coil_ind * J0 * CoefficientFunction((-y, x, 0)) / r_xy
    a = BilinearForm(fes)
    a += nu * curl(u) * curl(v) * dx + 1e-6 * NU0 * u * v * dx
    f = LinearForm(fes)
    f += InnerProduct(Jvec, v) * dx(definedon=mesh.Materials("coil"))
    a.Assemble(); f.Assemble()
    gfu = GridFunction(fes)
    gfu.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * f.vec
    return gfu


def validate_coil_linear_iron_force():
    """Case 2: coil + linear-iron sphere (mu_r=1000), force via eggshell.
    Recorded: NGSolve F_z = -0.1937 N, ref -0.1995 N (dipole -0.207)."""
    geo = CSGeometry()
    coil = _ring(50 * MM, 2 * MM, 0, 2 * MM).mat("coil").maxh(3 * MM)
    c = Pnt(0, 0, 30 * MM)
    iron = Sphere(c, 5 * MM).mat("iron").maxh(1.2 * MM)
    shell = (Sphere(c, 10 * MM) - Sphere(c, 5 * MM)).mat("shell").maxh(2 * MM)
    box = OrthoBrick(Pnt(-150 * MM, -150 * MM, -150 * MM),
                     Pnt(150 * MM, 150 * MM, 150 * MM)).bc("outer")
    air = (box - coil - iron - shell).mat("air")
    geo.Add(iron); geo.Add(shell); geo.Add(coil); geo.Add(air)
    mesh = Mesh(geo.GenerateMesh(maxh=25 * MM))

    iron_ind = CoefficientFunction([1.0 if m == "iron" else 0.0 for m in mesh.GetMaterials()])
    nu = iron_ind / (MU0 * 1000.0) + (1.0 - iron_ind) * NU0
    J0 = 10000.0 / ((2 * 2 * MM) * (2 * 2 * MM))
    B = curl(_solve_Aform(mesh, nu, J0))

    Fx, Fy, Fz = eggshell_force(B, mesh, center=(0, 0, 30 * MM),
                                r_inner=5 * MM, r_outer=10 * MM, air_region="shell")
    assert Fz < 0.0                                  # attractive toward the coil
    assert abs(Fz - (-0.1937)) < 0.03                # regression band vs recorded
    assert abs(Fx) < 0.15 * abs(Fz)                  # transverse ~0 by axisymmetry
    assert abs(Fy) < 0.15 * abs(Fz)


def validate_solenoid_self_inductance():
    """Cases 4-5: finite solenoid, self-inductance via field energy.
    Recorded: NGSolve L = 269.60 uH, ref 269.62 uH (Nagaoka 279.8)."""
    geo = CSGeometry()
    coil = _ring(30 * MM, 1 * MM, 0, 50 * MM).mat("coil").maxh(4 * MM)
    box = OrthoBrick(Pnt(-150 * MM, -150 * MM, -150 * MM),
                     Pnt(150 * MM, 150 * MM, 150 * MM)).bc("outer")
    air = (box - coil).mat("air")
    geo.Add(coil); geo.Add(air)
    mesh = Mesh(geo.GenerateMesh(maxh=20 * MM))

    J0 = 10000.0 / ((2 * 1 * MM) * (2 * 50 * MM))
    B = curl(_solve_Aform(mesh, NU0, J0))

    L_uH = self_inductance(B, mesh, i_terminal=100.0) * 1e6
    assert abs(L_uH - 269.6) < 8.0                   # regression band vs recorded


def validate_coil_nonlinear_iron_force():
    """Case 3: coil + NONLINEAR saturating-iron sphere, force via the shipped
    Picard solver ``solve_magnetostatic_nonlinear``. Iron permeability follows
    the saturation law mu_r(B) = 1 + 999/(1 + (|B|/1 T)^8) (mur0=1000) -- the SAME
    curve ref used (FullyCoupled Newton, Itot = 50000 A-turns).

    The iron sits DIRECTLY in air (NO nested 'shell' material): a nested shell
    isolates the iron interior (B_iron -> 0 at every sweep) so the saturation law
    never engages -- that #25 artifact is what produced the old, wrong -4.842 N.
    Recorded (iron-in-air, relax 0.35): NGSolve B_iron = 1.138 T, F_z = -4.898 N;
    ref B_iron = 1.13 T, F_z = -4.930 N => agree 0.65 %. ~20 Picard sweeps on
    ~37k tets -> multi-minute, so this stays in validation rather than pytest.
    """
    geo = CSGeometry()
    coil = _ring(50 * MM, 2 * MM, 0, 2 * MM).mat("coil").maxh(3 * MM)
    c = Pnt(0, 0, 30 * MM)
    iron = Sphere(c, 5 * MM).mat("iron").maxh(1.2 * MM)
    box = OrthoBrick(Pnt(-150 * MM, -150 * MM, -150 * MM),
                     Pnt(150 * MM, 150 * MM, 150 * MM)).bc("outer")
    air = (box - coil - iron).mat("air")             # iron DIRECTLY in air (no shell)
    geo.Add(iron); geo.Add(coil); geo.Add(air)
    mesh = Mesh(geo.GenerateMesh(maxh=25 * MM))

    iron_ind = CoefficientFunction([1.0 if m == "iron" else 0.0 for m in mesh.GetMaterials()])
    coil_ind = CoefficientFunction([1.0 if m == "coil" else 0.0 for m in mesh.GetMaterials()])
    r_xy = sqrt(x * x + y * y + 1e-12)
    J0 = 50000.0 / ((2 * 2 * MM) * (2 * 2 * MM))
    Jvec = coil_ind * J0 * CoefficientFunction((-y, x, 0)) / r_xy

    mur0, Bk, pp = 1000.0, 1.0, 8.0

    def nu_of_B(B):
        Bmag = sqrt(InnerProduct(B, B) + 1e-12)
        mur = 1.0 + (mur0 - 1.0) / (1.0 + (Bmag / Bk) ** pp)
        return iron_ind / (MU0 * mur) + (1.0 - iron_ind) * NU0

    vol_iron = Integrate(CoefficientFunction(1) * dx(definedon=mesh.Materials("iron")), mesh)
    iron_probe = (lambda B: Integrate(sqrt(InnerProduct(B, B))
                  * dx(definedon=mesh.Materials("iron")), mesh) / vol_iron)

    # defaults relax=0.35 / tol=1e-4 are the validated (stable) setting; a larger
    # relax overshoots into the low-mu_r regime and DIVERGES (see solve.py).
    gfu = solve_magnetostatic_nonlinear(mesh, nu_of_B, source=Jvec, monitor=iron_probe)
    B = curl(gfu)
    Biron = iron_probe(B)
    Fx, Fy, Fz = eggshell_force(B, mesh, center=(0, 0, 30 * MM),
                                r_inner=5 * MM, r_outer=10 * MM, air_region="air")
    assert 1.05 < Biron < 1.25                        # iron saturated (~1.138 T): law engaged
    assert Fz < 0.0                                   # attractive toward the coil
    assert abs(Fz - (-4.898)) < 0.2                   # regression band vs recorded (ref -4.930)
    assert abs(Fx) < 0.1 * abs(Fz)                    # transverse ~0 by axisymmetry
    assert abs(Fy) < 0.1 * abs(Fz)


def validate_helmholtz_center_field():
    """Case 6: Helmholtz pair (R = 50 mm coils at z = +/-25 mm, separation = R,
    NI = 10000 A-turns each, SAME sense), center field + flat-field uniformity,
    via the shipped ``solve_magnetostatic_Aform`` + ``azimuthal_coil_current``.
    Recorded: NGSolve B_center = 0.17407 T, ref 0.17420 T (agree 0.075 %);
    ideal thin-coil analytic (4/5)^1.5 mu0 NI/R = 0.17984 T (both ~3.2 % below,
    finite 4x4 mm section). The Helmholtz flat-field gives B_z/B_center ~ 0.998
    out to r = 10 mm and z = 10 mm. Linear iron-free single solve (~1-2 min)."""
    geo = CSGeometry()
    ring_p = _ring(50 * MM, 2 * MM, 25 * MM, 2 * MM)
    ring_m = _ring(50 * MM, 2 * MM, -25 * MM, 2 * MM)
    coil = (ring_p + ring_m).mat("coil").maxh(3 * MM)
    box = OrthoBrick(Pnt(-150 * MM, -150 * MM, -150 * MM),
                     Pnt(150 * MM, 150 * MM, 150 * MM)).bc("outer")
    air = (box - coil).mat("air")
    geo.Add(coil); geo.Add(air)
    mesh = Mesh(geo.GenerateMesh(maxh=20 * MM))

    J = azimuthal_coil_current(mesh, total_current=10000.0, section_area=(2 * 2 * MM) * (2 * 2 * MM))
    B = curl(solve_magnetostatic_Aform(mesh, NU0, source=J))
    bz_c = B(mesh(0, 0, 0))[2]
    bz_z = B(mesh(0, 0, 10 * MM))[2]
    bz_r = B(mesh(10 * MM, 0, 0))[2]
    assert abs(bz_c - 0.1741) < 0.012                 # regression band vs recorded (ref 0.17420)
    assert bz_z / bz_c > 0.99                          # flat field along axis (Helmholtz)
    assert bz_r / bz_c > 0.99                          # flat field off axis


def validate_solenoid_axis_field():
    """Case 4: finite solenoid (mean R = 30 mm, 2 mm radial, length 100 mm,
    NI = 10000 A-turns) on-axis B_z(z), via the shipped ``solve_magnetostatic_Aform``
    + ``azimuthal_coil_current``. Recorded NGSolve / ref / analytic [T]:
        z=0 : 0.10701 / 0.10721 / 0.10776
        z=15: 0.10401 / 0.10406 / 0.10475
        z=30: 0.09300 / 0.09274 / 0.09368
        z=45: 0.06947 / 0.06923 / 0.07025
    reference <-> NGSolve agree <0.35 %; both ~0.7 % below the thin-solenoid formula
    (finite coil thickness). Same solenoid as validate_solenoid_self_inductance."""
    geo = CSGeometry()
    coil = _ring(30 * MM, 1 * MM, 0, 50 * MM).mat("coil").maxh(4 * MM)
    box = OrthoBrick(Pnt(-150 * MM, -150 * MM, -150 * MM),
                     Pnt(150 * MM, 150 * MM, 150 * MM)).bc("outer")
    air = (box - coil).mat("air")
    geo.Add(coil); geo.Add(air)
    mesh = Mesh(geo.GenerateMesh(maxh=20 * MM))

    J = azimuthal_coil_current(mesh, total_current=10000.0, section_area=(2 * 1 * MM) * (2 * 50 * MM))
    B = curl(solve_magnetostatic_Aform(mesh, NU0, source=J))
    recorded = {0: 0.10701, 15: 0.10401, 30: 0.09300, 45: 0.06947}
    for zmm, b_rec in recorded.items():
        bz = B(mesh(0, 0, zmm * MM))[2]
        assert abs(bz - b_rec) < 0.006                 # regression band vs recorded (ref agrees <0.35%)


def validate_magnetized_sphere_field():
    """Case 1: uniformly magnetized sphere (Br = 1 T along z, radius 10 mm) in a
    200 mm air cube, via the shipped permanent-magnet source ``remanent_source``
    -> ``solve_magnetostatic_Aform(curl_source=...)``. The interior field is
    uniform and analytic:  B_z = (2/3) Br = 0.66667 T. Recorded: NGSolve volume-
    averaged <B_z> = 0.66391 T, ref 0.66466 T (agree 0.11 %). This is the
    canonical analytic warm-up -- it pins the A-form PM source sign + scaling."""
    geo = CSGeometry()
    magnet = Sphere(Pnt(0, 0, 0), 10 * MM).mat("magnet").maxh(2 * MM)
    box = OrthoBrick(Pnt(-100 * MM, -100 * MM, -100 * MM),
                     Pnt(100 * MM, 100 * MM, 100 * MM)).bc("outer")
    air = (box - magnet).mat("air")
    geo.Add(magnet); geo.Add(air)
    mesh = Mesh(geo.GenerateMesh(maxh=15 * MM))

    curl_src = remanent_source(mesh, {"magnet": (0.0, 0.0, 1.0)})   # Br = 1 T z-hat
    B = curl(solve_magnetostatic_Aform(mesh, NU0, curl_source=curl_src))
    vol = Integrate(CoefficientFunction(1) * dx(definedon=mesh.Materials("magnet")), mesh)
    bz_avg = Integrate(B[2] * dx(definedon=mesh.Materials("magnet")), mesh) / vol
    assert abs(bz_avg - 0.6667) < 0.02                 # interior (2/3)Br; recorded NGSolve 0.6639 / ref 0.6647


def validate_pm_sphere_exterior_dipole_field():
    """Uniformly magnetized PM sphere (Br=1 T, R=10 mm) — exterior dipole field on axis.
    The exact exterior field of a uniformly magnetized sphere is a pure dipole:
        B_z(z) = (μ₀/4π) × 2m / z³  where  m = (Br/μ₀)(4π/3)R³.

    Pass criterion: |B_z_fem − B_z_ref| < 10 % of B_z_ref at z = 20 mm and 35 mm.
    """
    R   = 10 * MM
    Br  = 1.0
    m   = (Br / MU0) * (4.0 * math.pi / 3.0) * R**3
    Lbox = 200 * MM

    geo = CSGeometry()
    magnet = Sphere(Pnt(0, 0, 0), R).mat("magnet").maxh(1.5 * MM)
    box    = OrthoBrick(Pnt(-Lbox, -Lbox, -Lbox), Pnt(Lbox, Lbox, Lbox)).bc("outer")
    air    = (box - magnet).mat("air")
    geo.Add(magnet); geo.Add(air)
    mesh   = Mesh(geo.GenerateMesh(maxh=20 * MM))

    curl_src = remanent_source(mesh, {"magnet": (0.0, 0.0, Br)})
    B = curl(solve_magnetostatic_Aform(mesh, NU0, curl_source=curl_src))

    for z_m in [20 * MM, 35 * MM]:
        bz_fem = B(mesh(1e-4, 0, z_m))[2]
        bz_ref = (MU0 / (4.0 * math.pi)) * 2.0 * m / z_m**3
        tol    = 0.10 * bz_ref
        assert abs(bz_fem - bz_ref) < tol, (
            f"PM dipole exterior Bz at z={z_m*1e3:.0f}mm: "
            f"fem={bz_fem*1e3:.2f} mT  ref={bz_ref*1e3:.2f} mT  "
            f"(tol {tol*1e3:.2f} mT)"
        )


def validate_single_loop_axis_field():
    """Case 7: single circular loop (R = 50 mm, 4x4 mm section, NI = 10000 AT)
    on-axis B_z(z), via the shipped ``solve_magnetostatic_Aform`` +
    ``azimuthal_coil_current``. Thin-filament analytic
    B_z = mu0 NI R^2 / (2 (R^2+z^2)^1.5) = 0.12566 / 0.05983 T at z = 0 / 40 mm;
    the finite 4x4 mm section puts NGSolve ~2-5 % under (measured 0.12278 /
    0.05702, ratio 0.977 / 0.953). ref agrees with NGSolve <0.8 % for z<=40mm."""
    geo = CSGeometry()
    coil = _ring(50 * MM, 2 * MM, 0, 2 * MM).mat("coil").maxh(3 * MM)
    box = OrthoBrick(Pnt(-150 * MM, -150 * MM, -150 * MM),
                     Pnt(150 * MM, 150 * MM, 150 * MM)).bc("outer")
    air = (box - coil).mat("air")
    geo.Add(coil); geo.Add(air)
    mesh = Mesh(geo.GenerateMesh(maxh=20 * MM))

    J = azimuthal_coil_current(mesh, total_current=10000.0, section_area=(2 * 2 * MM) * (2 * 2 * MM))
    B = curl(solve_magnetostatic_Aform(mesh, NU0, source=J))
    R = 50 * MM
    recorded = {0: 0.12278, 40: 0.05702}                 # measured NGSolve (analytic 0.12566 / 0.05983)
    for zmm, b_rec in recorded.items():
        z = zmm * MM
        bz = B(mesh(0, 0, z))[2]
        ana = MU0 * 10000.0 * R ** 2 / (2 * (R ** 2 + z ** 2) ** 1.5)
        assert abs(bz - b_rec) < 0.004                   # regression band vs measured
        assert bz < ana                                  # finite section lowers the center field


def validate_team20_static_force():
    """TEAM Workshop Problem 20 (3-D static force) at NI=1000, via the shipped
    Newton solver ``solve_magnetostatic_newton`` (energy-min + CG/BDDC) on the
    official geometry (yoke + center pole + rectangular coil, 1.5 mm bottom gap),
    NONLINEAR steel B-H (38 pts + Ms=2.16 T extension), force = air-side Maxwell
    stress on the pole surface. Measured Fz = 8.1 N; this pipeline gives ~7.9 N on
    this coarse CI mesh and 8.23 N on the production mesh (+1.7 %), beating the
    lab's own NGSolve A-phi (7.6) and Radia (6.86). Coarse mesh keeps it ~15 s."""
    from netgen.occ import Box, Pnt, Glue, OCCGeometry
    from netgen.meshing import MeshingParameters
    yx, yy, yz = 0.127, 0.025, 0.150
    cav1 = (0.025, 0.102, 0.025, 0.125); cav2 = (0.050, 0.077, 0.125, 0.150)
    pg = (0.051, 0.076, 0.0075, 0.0175, 0.0265, 0.125)
    cgx, cgy, ct = 0.007, 0.0145, 0.018
    ix0, ix1 = pg[0] - cgx, pg[1] + cgx; iy0, iy1 = pg[2] - cgy, pg[3] + cgy
    ox0, ox1 = ix0 - ct, ix1 + ct; oy0, oy1 = iy0 - ct, iy1 + ct
    cz0, cz1 = 0.0267, 0.1233; e = 1e-4; A_coil = ct * (cz1 - cz0)

    def _ys():
        return (Box(Pnt(0, 0, 0), Pnt(yx, yy, yz))
                - Box(Pnt(cav1[0], -e, cav1[2]), Pnt(cav1[1], yy + e, cav1[3]))
                - Box(Pnt(cav2[0], -e, cav2[2]), Pnt(cav2[1], yy + e, cav2[3])))
    def _ps():
        return Box(Pnt(pg[0], pg[2], pg[4]), Pnt(pg[1], pg[3], pg[5]))
    def _cs():
        return (Box(Pnt(ox0, oy0, cz0), Pnt(ox1, oy1, cz1))
                - Box(Pnt(ix0, iy0, cz0 - e), Pnt(ix1, iy1, cz1 + e)))

    yoke = _ys(); yoke.mat('yoke')
    pole = _ps(); pole.mat('pole')
    for f in pole.faces:
        f.name = 'pole_surface'
    coil = _cs(); coil.mat('coil')
    pad = 0.08
    big = Box(Pnt(-pad, oy0 - pad, -pad), Pnt(yx + pad, oy1 + pad, yz + pad))
    for f in big.faces:
        f.name = 'outer'
    air = big - _ys() - _ps() - _cs(); air.mat('air')
    yoke.faces.maxh = 0.012; pole.faces.maxh = 0.012
    coil.faces.maxh = 0.016; air.maxh = 0.050
    mesh = Mesh(OCCGeometry(Glue([yoke, pole, coil, air])).GenerateMesh(
        MeshingParameters(maxh=0.050, grading=0.3))); mesh.Curve(2)
    assert "pole_surface" in set(mesh.GetBoundaries())

    # nonlinear steel B-H -> co-energy density phi(B) = int H dB
    Bd = [0, 0.01, 0.025, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70,
          0.80, 0.90, 1.0, 1.1, 1.2, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55, 1.6, 1.65,
          1.7, 1.75, 1.8, 1.85, 1.9, 1.95, 2.0, 2.05, 2.1, 2.15, 2.2, 2.25, 2.3]
    Hd = [0, 27, 58, 100, 153, 185, 205, 233, 255, 285, 320, 355, 405, 470, 555,
          673, 836, 1065, 1220, 1420, 1720, 2130, 2670, 3480, 4500, 5950, 7650,
          10100, 13000, 15900, 21100, 26300, 32900, 42700, 61700, 84300, 110000,
          135000]
    Ms = 2.16
    for bb in (2.4, 2.6, 2.8, 3.0, 3.5, 4.0):
        Bd.append(bb); Hd.append((bb - Ms) / MU0)
    energy_density = BSpline(2, [0] + list(Bd), list(Hd)).Integrate()

    NI = 1000.0; jj = NI / A_coil
    cxc = (ix0 + ix1) / 2; cyc = (iy0 + iy1) / 2
    rx = x - cxc; ry = y - cyc
    rxa = sqrt(rx * rx + 1e-20); rya = sqrt(ry * ry + 1e-20)
    jx = IfPos(rya - rxa, IfPos(ry, -jj, jj), 0.0)
    jy = IfPos(rxa - rya, IfPos(rx, jj, -jj), 0.0)
    Jv = mesh.MaterialCF({"coil": CoefficientFunction((jx, jy, 0.0))},
                         default=CoefficientFunction((0.0, 0.0, 0.0)))

    gfu = solve_magnetostatic_newton(mesh, Jv, energy_density, "yoke|pole",
                                     load_steps=3, max_iter=25, tol=1e-7)
    B = curl(gfu)
    Vp = Integrate(CoefficientFunction(1) * dx(definedon=mesh.Materials("pole")), mesh)
    Bp = Integrate(sqrt(InnerProduct(B, B)) * dx(definedon=mesh.Materials("pole")), mesh) / Vp

    fl = L2(mesh, order=1, definedon=mesh.Materials("air"))
    gB = [GridFunction(fl) for _ in range(3)]
    for k in range(3):
        gB[k].Set(B[k])
    Ba = CoefficientFunction(tuple(BoundaryFromVolumeCF(gB[k]) for k in range(3)))
    n = specialcf.normal(3)
    Tz = (1.0 / MU0) * (Ba[2] * InnerProduct(Ba, n) - 0.5 * InnerProduct(Ba, Ba) * n[2])
    Fz = -Integrate(Tz * ds(definedon=mesh.Boundaries("pole_surface")), mesh)

    assert 0.45 < Bp < 0.65                            # pole working point ~0.53 T (unsaturated)
    assert Fz > 0.0                                    # pole pulled toward closing the gap
    assert abs(Fz - 8.0) < 1.6                         # measured 8.1; coarse 7.89 / full mesh 8.23


def validate_team13_bflux():
    """TEAM Workshop Problem 13 (3-D nonlinear C-core) at NI=1000 AT, via the shipped
    Newton solver ``solve_magnetostatic_newton`` (energy-min + CG/BDDC) on Yano's
    gap-resolved 37k-element half-model mesh (z >= 0 symmetry plane).

    The critical geometry feature is the 0.5 mm air gap between the vertical iron
    sheet (x = +/-1.6 mm) and the C-shaped return-flux sheet (x = +/-2.1 mm).
    Yano's ``getMeshingParameters`` seeds maxh = 0.8 mm near the gap corners; without
    this refinement the FEM gives |B| roughly 15x too low (13k coarse elements vs
    37k with gap seeds).

    Checks average |Bz| in the vertical sheet cross-section (x in [0, 1.6 mm],
    y in [-25, 25 mm]) at z = 1, 30, 60 mm against TEAM-13 measured values
    1.330, 1.225, 0.655 T.  Solver (bddc, reg=1e-6, 8 Newton iters, ~2 min):
    1.308 / 1.231 / 0.657 T -- all within 2 % of measured.
    Skipped if Yano's Problem13 geometry module is unavailable on this machine.
    """
    import os
    import sys
    import numpy as np
    from ngsolve import TaskManager

    # The external TEAM-13 geometry module is a local, machine-specific asset;
    # its location is taken from an env var (no internal path hardcoded). The
    # test skips cleanly if it is not available on this machine.
    GEO_PATH = os.environ.get("RADIA_TEAM13_GEOMETRY", "")
    if GEO_PATH and GEO_PATH not in sys.path:
        sys.path.insert(0, GEO_PATH)
    try:
        from geometry import generateGeometry, getMeshingParameters
    except ImportError:
        pytest.skip("TEAM-13 geometry module not available "
                    "(set RADIA_TEAM13_GEOMETRY to its directory)")

    geo  = generateGeometry(maxh_iron=100, fullProblem=False)
    mp   = getMeshingParameters(maxh_global=100)
    mesh = Mesh(geo.GenerateMesh(mp=mp))
    mesh.Curve(5)

    # rectangular coil current: 4 straight bricks (constant J) + 4 arc corners (tangential J)
    NI = 1000.0
    J0 = NI / 2500e-6          # A_coil = 50 mm x 50 mm = 2500 mm^2
    mats = mesh.GetMaterials()
    x_rb = x - 0.050; y_rb = y - 0.050; r_rb = (x_rb**2 + y_rb**2)**(0.5)
    x_lb = x + 0.050; y_lb = y - 0.050; r_lb = (x_lb**2 + y_lb**2)**(0.5)
    x_lf = x + 0.050; y_lf = y + 0.050; r_lf = (x_lf**2 + y_lf**2)**(0.5)
    x_rf = x - 0.050; y_rf = y + 0.050; r_rf = (x_rf**2 + y_rf**2)**(0.5)
    _Jm = {
        "corner_right_back":  [y_rb / r_rb, -x_rb / r_rb],
        "corner_left_back":   [y_lb / r_lb, -x_lb / r_lb],
        "corner_left_front":  [y_lf / r_lf, -x_lf / r_lf],
        "corner_right_front": [y_rf / r_rf, -x_rf / r_rf],
        "brick_back": [1, 0], "brick_front": [-1, 0],
        "brick_left": [0, 1], "brick_right": [0, -1],
    }
    J = (J0 * CoefficientFunction([_Jm[m][0] if m in _Jm else 0 for m in mats])
             * CoefficientFunction((1, 0, 0))
       + J0 * CoefficientFunction([_Jm[m][1] if m in _Jm else 0 for m in mats])
             * CoefficientFunction((0, 1, 0)))

    # 77-point B-H curve (Kruger-Lehmann steel, from Yano's notebook)
    H_KL = [
        -4.47197834e-13, 1.60e+01, 3.00e+01, 5.40e+01, 9.30e+01, 1.43e+02, 1.91e+02,
        2.10e+02, 2.22e+02, 2.33e+02, 2.47e+02, 2.58e+02, 2.72e+02, 2.89e+02, 3.13e+02,
        3.42e+02, 3.77e+02, 4.33e+02, 5.09e+02, 6.48e+02, 9.33e+02, 1.228e+03, 1.934e+03,
        2.913e+03, 4.993e+03, 7.189e+03, 9.423e+03, 9.423e+03,
        1.28203768e+04, 1.65447489e+04, 2.07163957e+04, 2.55500961e+04, 3.15206135e+04,
        4.03204637e+04, 7.73038295e+04, 1.29272791e+05, 1.81241752e+05, 2.33210713e+05,
        2.85179674e+05, 3.37148635e+05, 3.89117596e+05, 4.41086557e+05, 4.93055518e+05,
        5.45024479e+05, 5.96993440e+05, 6.48962401e+05, 7.00931362e+05, 7.52900323e+05,
        8.04869284e+05, 8.56838245e+05, 9.08807206e+05, 9.60776167e+05, 1.01274513e+06,
        1.06471409e+06, 1.11668305e+06, 1.16865201e+06, 1.22062097e+06, 1.27258993e+06,
        1.32455889e+06, 1.37652785e+06, 1.42849682e+06, 1.48046578e+06, 1.53243474e+06,
        1.58440370e+06, 1.63637266e+06, 1.68834162e+06, 1.74031058e+06, 1.79227954e+06,
        1.84424850e+06, 1.89621746e+06, 1.94818643e+06, 2.00015539e+06, 2.05212435e+06,
        2.10409331e+06, 2.15606227e+06, 2.20803123e+06, 2.26000019e+06,
    ]
    B_KL = [
        0.0, 2.5e-3, 5.0e-3, 1.25e-2, 2.5e-2, 5.0e-2, 1.0e-1, 2.0e-1, 3.0e-1, 4.0e-1,
        5.0e-1, 6.0e-1, 7.0e-1, 8.0e-1, 9.0e-1, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.55,
        1.6, 1.65, 1.7, 1.75, 1.8, 1.8,
        1.86530612, 1.93061224, 1.99591837, 2.06122449, 2.12653061, 2.19183673,
        2.25714286, 2.32244898, 2.38775510, 2.45306122, 2.51836735, 2.58367347,
        2.64897959, 2.71428571, 2.77959184, 2.84489796, 2.91020408, 2.97551020,
        3.04081633, 3.10612245, 3.17142857, 3.23673469, 3.30204082, 3.36734694,
        3.43265306, 3.49795918, 3.56326531, 3.62857143, 3.69387755, 3.75918367,
        3.82448980, 3.88979592, 3.95510204, 4.02040816, 4.08571429, 4.15102041,
        4.21632653, 4.28163265, 4.34693878, 4.41224490, 4.47755102, 4.54285714,
        4.60816327, 4.67346939, 4.73877551, 4.80408163, 4.86938776, 4.93469388, 5.0,
    ]
    energy_density = BSpline(2, [0] + list(B_KL), list(H_KL)).Integrate()

    with TaskManager():
        gfu = solve_magnetostatic_newton(
            mesh, J, energy_density, "iron",
            reg=1e-6, tol=1e-10, max_iter=15,
        )
    B = curl(gfu)

    # average |Bz| in vertical sheet cross-section at each z plane
    sc = 1e-3
    xi_ = np.linspace(0, 1.6, 10) * sc        # x in [0, 1.6 mm]
    yi_ = np.linspace(-25, 25, 20) * sc        # y in [-25, 25 mm]
    A_sh = 50.0 * 1.6 * sc * sc               # sheet area = 50 mm x 1.6 mm
    for zk, ref in [(0.001, 1.330), (0.030, 1.225), (0.060, 0.655)]:
        vals = np.zeros((len(xi_), len(yi_)))
        for i, xv in enumerate(xi_):
            for j, yv in enumerate(yi_):
                try:
                    vals[i, j] = float(B[2].Norm()(mesh(xv, yv, zk)))
                except Exception:
                    pass
        Bz = np.trapz(np.trapz(vals, yi_), xi_) / A_sh
        assert abs(Bz - ref) < 0.10, (
            f"TEAM-13 z={zk * 1000:.0f}mm: got {Bz:.3f} T, ref {ref:.3f} T"
        )


def validate_team7_eddy_bz():
    """TEAM Workshop Problem 7 (3-D eddy current in a thick aluminium plate,
    50 Hz sinusoidal, 2742-turn rectangular coil) — frequency-domain A-Φ
    formulation with Kelvin open-boundary transform.

    Uses the shipped ``solve_eddy_current_harmonic_APhi`` (complex HCurl×H1,
    Kelvin Periodic identification, GMRes + "local" preconditioner, order=3).
    Solver validated by Yano 2026-03 at order=5 on the same mesh (mesh_maxh=0.05 m);
    order=3 gives a ~4x speed-up with ~20 % accuracy reduction.

    Checks Bz at 17 reference points on measurement line A1-B1 (y=72 mm, z=34 mm)
    against TEAM-7 measured values (Gauss).  Pass criterion: complex RMSE < 12 G
    (reference peak ~78 G, so ~15 % relative).

    Sign convention (Yano plot): Bz_Gauss = -1e4 * B[2].real + 1j * 1e4 * B[2].imag
    (NGSolve's B[2] is in Tesla; the minus sign on real part matches the e^{jωt}
    phase convention used in the reference measurements).
    """
    import numpy as np
    from netgen.occ import (Box, Glue, OCCGeometry,
                            Sphere as OccSphere, Pnt as OccPnt,
                            Vertex, WorkPlane, Axes,
                            X as OccX, Y as OccY, Z as OccZ,
                            IdentificationType)

    mu0      = 4 * np.pi * 1e-7
    nu0      = 1.0 / mu0
    turns    = 2742
    omega    = 2 * np.pi * 50.0
    sigma_Al = 35.26e6        # S/m (TEAM-7 aluminium spec)

    # Kelvin transform geometry (Yano's validated values)
    sphere_cx, sphere_cy, sphere_cz = 0.147, 0.147, 0.050
    kelvin_R   = 0.30
    kelvin_cap = 0.28         # box clip half-width to remove coordinate singularity
    outer_cx   = sphere_cx + 2 * kelvin_R
    big        = 2 * kelvin_R

    # --- Aluminium plate: 294×294 mm outer face, 108×108 mm hole, 19 mm thick ---
    wpplate = WorkPlane(Axes(OccPnt(0, 0, 0), OccZ, OccX))
    f_plate = (wpplate.Rectangle(0.294, 0.294)
                      .MoveTo(0.018, 0.018)
                      .Rectangle(0.108, 0.108)
                      .Reverse()
                      .Face())
    plate = f_plate.Extrude(0.019)
    plate.faces.maxh  = 0.04
    plate.solids.name = "plate"

    # --- Coil: 100×100 mm cross-section (inner+outer), 100 mm height, 2742 turns ---
    wp = WorkPlane(Axes(p=OccPnt(0.294, 0.000, 0.049), n=OccZ, h=OccY))
    coil2di = wp.MoveTo(0.100, 0.100).RectangleC(0.100, 0.100).Offset(0.025).Face()
    coil2di.edges.name = "coili"    # inner boundary surface (source integration site)
    coil2do = wp.MoveTo(0.100, 0.100).RectangleC(0.100, 0.100).Offset(0.050).Face()
    coil2d  = coil2do - coil2di
    cutout  = (wp.MoveTo(0.100, 0.100).RectangleC(0.300, 0.100).Face() +
               wp.MoveTo(0.100, 0.100).RectangleC(0.100, 0.300).Face())
    coil2d  = Glue([coil2d * cutout, coil2d - cutout])
    coil    = coil2d.Extrude(0.100)
    coil.solids.name = "coil"

    # --- Inner Kelvin sphere (2-cap clipped to avoid coord singularities) ---
    inner_sphere = (
        OccSphere(OccPnt(sphere_cx, sphere_cy, sphere_cz), kelvin_R)
        * Box(OccPnt(sphere_cx - kelvin_cap, sphere_cy - kelvin_cap, sphere_cz - big),
              OccPnt(sphere_cx + big,        sphere_cy + big,        sphere_cz + big)))
    air = inner_sphere - plate - coil
    air.mat("air")

    # --- Outer Kelvin sphere (same clipped shape, offset +2R in x) ---
    air_outer = (
        OccSphere(OccPnt(outer_cx, sphere_cy, sphere_cz), kelvin_R)
        * Box(OccPnt(outer_cx - kelvin_cap, sphere_cy - kelvin_cap, sphere_cz - big),
              OccPnt(outer_cx + big,        sphere_cy + big,        sphere_cz + big)))
    air_outer.mat("air_outer")

    # Kelvin periodic identification: largest dome face inner <-> outer
    int_dome = max(air.faces,       key=lambda f: f.mass)
    ext_dome = max(air_outer.faces, key=lambda f: f.mass)
    int_dome.name = "kelvin_int"
    ext_dome.name = "kelvin_ext"
    ext_dome.Identify(int_dome, "periodic", IdentificationType.PERIODIC)

    # GND vertex at outer sphere centre (Φ=0 gauge for the H1 space)
    gnd      = Vertex(OccPnt(outer_cx, sphere_cy, sphere_cz))
    gnd.name = "GND"

    ngmesh = OCCGeometry(Glue([plate, coil, air, air_outer, gnd])).GenerateMesh(maxh=0.05)
    mesh   = Mesh(ngmesh)
    mesh.Curve(3)

    # --- Reluctivity CF (includes Kelvin outer-domain scaling) ---
    r_prime   = sqrt((x - outer_cx)**2 + (y - sphere_cy)**2 + (z - sphere_cz)**2 + 1e-20)
    nu_kelvin = (r_prime / kelvin_R)**2 * nu0
    _nu_map   = {"air": nu0, "plate": nu0, "coil": nu0, "air_outer": nu_kelvin}
    nu    = CoefficientFunction([_nu_map.get(m, nu0) for m in mesh.GetMaterials()])
    sigma = mesh.MaterialCF({"plate": sigma_Al}, default=0.0)

    # --- Coil current source: tau_coil / pot_coil (Yano convention) ---
    # tau_coil: unit tangential winding direction in the rectangular coil
    # pot_coil: vector potential for the curl-source term
    cx_min, cx_max = 0.144, 0.244    # x extent of inner coil square
    cy_min, cy_max = 0.050, 0.150    # y extent of inner coil square
    projx    = IfPos(x - cx_min, IfPos(x - cx_max, cx_max, x), cx_min)
    projy    = IfPos(y - cy_min, IfPos(y - cy_max, cy_max, y), cy_min)
    tau_raw  = CoefficientFunction((projy - y, x - projx, 0.0))
    tau_coil = tau_raw / sqrt(InnerProduct(tau_raw, tau_raw) + 1e-30)
    pot_coil = CoefficientFunction(
        (0.0, 0.0, sqrt((projy - y)**2 + (projx - x)**2) - 0.050))

    def add_source(lf, validate_fn):
        v_A, _ = validate_fn
        lf += (-turns / 0.100 * tau_coil * v_A.Trace()
               * ds("coili", bonus_intorder=4))
        lf += (turns / (0.025 * 0.100) * pot_coil * curl(v_A)
               * dx("coil", bonus_intorder=4))

    # --- Solve A-Phi eddy current ---
    gfAPhi = solve_eddy_current_harmonic_APhi(
        mesh, nu, sigma, omega, add_source, order=3,
    )
    B = curl(gfAPhi.components[0])

    # --- Reference values (TEAM-7 measured, Gauss) on A1-B1 (y=72mm, z=34mm) ---
    xi_msm = np.array([0, 18, 36, 54, 72, 90, 108, 126, 144, 162,
                       180, 198, 216, 234, 252, 270, 288]) * 1e-3
    Bz_A1B1_ref = np.array([
        -4.9  -1.16j, -17.88 +2.48j, -22.13 +4.15j, -20.19 +4.j,   -15.67 +3.07j,
         0.36 +2.31j,  43.64 +1.89j,  78.11 +4.97j,  71.55+12.61j,  60.44+14.15j,
        53.91+13.04j,  52.62+12.4j,   53.81+12.05j,  56.91+12.27j,  59.24+12.66j,
        52.78 +9.96j,  27.61 +2.26j])

    # NGSolve gives B in Tesla; convert to Gauss with Yano's sign convention
    Bz_raw = np.array([complex(B[2](mesh(xi, 0.072, 0.034))) for xi in xi_msm])
    Bz_G   = -1e4 * Bz_raw.real + 1j * 1e4 * Bz_raw.imag

    rmse = float(np.sqrt(np.mean(np.abs(Bz_G - Bz_A1B1_ref) ** 2)))
    assert rmse < 12.0, (
        f"TEAM-7 A1-B1 RMSE = {rmse:.1f} G (limit 12 G); "
        f"peak sim = {-1e4 * Bz_raw[7].real:.1f} G (ref 78.1 G)"
    )


def validate_team21_eddy_loss():
    """TEAM Workshop Problem 21 P21a-0 — frequency-domain eddy current in a
    non-magnetic steel plate (20Mn23Al, σ=1.3889×10⁶ S/m, μr=1) driven by two
    opposite-phase rectangular coils (Case I: flux perpendicular to plate), 50 Hz,
    10 A rms.

    Geometry: coil outer=248 mm, inner=200 mm, height=217 mm×2, gap=24 mm.
    Plate: 458×820×10 mm, placed 12 mm from the coil face.
    Source: volume current ``J_scale × tau_coil`` with opposite signs (Case I).
    Solver: ``solve_eddy_current_harmonic_APhi`` with ``periodic=False``
    (plain Dirichlet box outer boundary, no Kelvin transform).

    Pass criterion: |P_eddy − 9.17 W| < 3.5 W (±38 %).
    Reference: TEAM P21 Family V.2009 Table A2-4
        measured  9.17 W, calculated 9.31 W (Case I, 10 A rms, 50 Hz).
    """
    import numpy as np
    from netgen.occ import Box, Pnt, Glue, OCCGeometry
    from netgen.meshing import MeshingParameters

    mu0         = 4 * np.pi * 1e-7
    nu0         = 1.0 / mu0
    sigma_steel = 1.3889e6       # 20Mn23Al  [S/m]
    omega       = 2 * np.pi * 50.0
    turns       = 300            # per coil
    I_rms       = 10.0           # A rms, Case I reference

    # Geometry (metres)
    half_o    = 0.124            # coil outer half-width  (248 mm / 2)
    half_i    = 0.100            # coil inner half-width  (200 mm / 2)
    coil_h    = 0.217            # each coil height
    coil_gap  = 0.024            # HV/LV inter-coil gap
    plate_gap = 0.012            # coil-face to plate
    plate_x0  = half_o + plate_gap        # 0.136 m
    plate_x1  = plate_x0 + 0.010         # 0.146 m  (10 mm thick)
    plate_y   = 0.410                     # half-width in y  (820 mm / 2)
    plate_z   = coil_h * 2 + coil_gap    # 0.458 m
    coil_w    = half_o - half_i          # 0.024 m  (winding radial width)
    margin    = 0.10

    # Coil 1: z = 0 … coil_h  (hollow annular prism, no fillets)
    coil1 = (Box(Pnt(-half_o, -half_o, 0.0),  Pnt(half_o, half_o, coil_h))
             - Box(Pnt(-half_i, -half_i, 0.0), Pnt(half_i, half_i, coil_h)))
    coil1.mat("coil1")

    # Coil 2: z = coil_h+coil_gap … 2*coil_h+coil_gap
    z2    = coil_h + coil_gap
    coil2 = (Box(Pnt(-half_o, -half_o, z2),  Pnt(half_o, half_o, z2 + coil_h))
             - Box(Pnt(-half_i, -half_i, z2), Pnt(half_i, half_i, z2 + coil_h)))
    coil2.mat("coil2")

    # Non-magnetic steel plate
    plate = Box(Pnt(plate_x0, -plate_y, 0.0), Pnt(plate_x1, plate_y, plate_z))
    plate.faces.maxh = 0.025
    plate.mat("plate")

    # Air box with Dirichlet outer boundary
    air_box = Box(
        Pnt(-half_o - margin, -plate_y - margin, -margin),
        Pnt(plate_x1 + margin,  plate_y + margin, plate_z + margin))
    air_box.faces.name = "outer"
    air = air_box - coil1 - coil2 - plate
    air.mat("air")

    ngmesh = OCCGeometry(Glue([coil1, coil2, plate, air])).GenerateMesh(
        MeshingParameters(maxh=0.07, grading=0.35))
    mesh = Mesh(ngmesh)
    mesh.Curve(2)

    # Material CFs  (μr = 1 everywhere → ν = ν₀ everywhere)
    mats     = mesh.GetMaterials()
    nu       = CoefficientFunction([nu0] * len(mats))
    sigma_cf = mesh.MaterialCF({"plate": sigma_steel}, default=0.0)

    # Coil current: tau_coil winds CCW around the inner 200-mm square
    projx    = IfPos(x - (-half_i), IfPos(x - half_i, half_i, x), -half_i)
    projy    = IfPos(y - (-half_i), IfPos(y - half_i, half_i, y), -half_i)
    tx       = projy - y
    ty       = x - projx
    tau_coil = CoefficientFunction((tx, ty, 0)) / sqrt(tx**2 + ty**2 + 1e-30)

    # Peak NI for 10 A rms sinusoidal; J_scale = peak current density in winding
    NI_peak  = turns * I_rms * math.sqrt(2)
    J_scale  = NI_peak / (coil_w * coil_h)   # [A/m²]

    def add_source(lf, validate_fn):
        v, _ = validate_fn
        lf += J_scale * tau_coil * v * dx("coil1")   # CCW  →  +z flux
        lf += (-J_scale) * tau_coil * v * dx("coil2") # CW   →  −z flux  (Case I)

    gfAPhi = solve_eddy_current_harmonic_APhi(
        mesh, nu, sigma_cf, omega, add_source,
        order=2, precond="local", dirichlet="outer",
        periodic=False, tol=1e-8, maxsteps=2000, restart=200, reg=1e-10,
    )

    # Time-averaged eddy loss: P = ½ σ |E|²  where  E = −jω(A + ∇Φ)
    gfA, gfPhi = gfAPhi.components
    E_cf   = -1j * omega * (gfA + grad(gfPhi))
    P_eddy = (0.5 * Integrate(
        sigma_cf * InnerProduct(E_cf, Conj(E_cf)).real, mesh)).real

    ref = 9.17   # W, measured (P21a-0, Case I, 10 A rms, 50 Hz)
    assert abs(P_eddy - ref) < 3.5, (
        f"TEAM-21 P21a-0 eddy loss = {P_eddy:.2f} W, "
        f"ref = {ref:.2f} W  (tol ±3.5 W)"
    )


def validate_team6_sphere_eddy():
    """TEAM Workshop Problem 6 — hollow conducting sphere (σ=5×10⁷ S/m, μᵣ=1)
    in a uniform sinusoidal magnetic field B₀ = 1 T at 50 Hz.

    Scattered-field A-Φ formulation: A = Ã + A₀ where A₀ = (B₀/2)(−y, x, 0)
    (so curl A₀ = B₀ẑ).  Solve for the scattered Ã with Ã = 0 on a Dirichlet
    box (periodic=False).  Total B = curl(Ã) + B₀ẑ.

    Analytical reference: exact BVP solution using modified spherical Bessel
    functions i₁(κr), k₁(κr) in the shell (κ = sqrt(jωσμ₀)).  The interior
    cavity field is uniform: B_z_inside = 2A where A is the l=1 coefficient.

    Skin depth δ ≈ 10 mm, shell thickness 5 mm (≈ δ/2) → significant but
    partial shielding (|B_in| ≈ 0.5 T, phase ≈ −60°).

    Pass criterion: |B_z_fem(centre) − B_z_ref| < 0.06 T  (6 % of B₀).
    """
    import numpy as np
    from netgen.csg import CSGeometry, Sphere, OrthoBrick, Pnt
    from netgen.meshing import MeshingParameters

    mu0   = 4 * np.pi * 1e-7
    nu0   = 1.0 / mu0
    sigma = 5e7          # copper-like  [S/m]
    omega = 2 * np.pi * 50.0
    B0    = 1.0          # applied uniform field [T]
    R1    = 0.05         # inner radius  [m]
    R2    = 0.055        # outer radius  [m]
    L     = 0.25         # Dirichlet box half-length  [m]

    # ------------------------------------------------------------------ #
    # Analytical reference: modified-SBF BVP for the n=1 (dipole) mode.  #
    # A_φ(r,θ) = f(r) sin θ;  B_z = 2·f(r)/r · const → inside B_z = 2A. #
    # Boundary conditions: f, f' continuous at r = R1 and r = R2.        #
    # ------------------------------------------------------------------ #
    kappa = np.sqrt(1j * omega * sigma * mu0)  # = (1+j)/delta, delta≈10mm

    def _i1(z):  return np.cosh(z)/z - np.sinh(z)/z**2
    def _k1(z):  return np.exp(-z) * (1/z + 1/z**2)
    def _i1p(z): return np.sinh(z)/z - 2*np.cosh(z)/z**2 + 2*np.sinh(z)/z**3
    def _k1p(z): return np.exp(-z) * (-1/z - 2/z**2 - 2/z**3)

    z1, z2 = kappa * R1, kappa * R2
    # Rows: [f cont @ R1, f' cont @ R1, f cont @ R2, f' cont @ R2]
    # Unknowns: [A (inner coeff), C, D (shell), E (outer dipole)]
    M = np.array([
        [R1,  -_i1(z1),          -_k1(z1),          0        ],
        [1,   -kappa * _i1p(z1), -kappa * _k1p(z1), 0        ],
        [0,    _i1(z2),           _k1(z2),           -1/R2**2 ],
        [0,    kappa * _i1p(z2),  kappa * _k1p(z2),  2/R2**3  ],
    ], dtype=complex)
    rhs = np.array([0, 0, B0 * R2 / 2, B0 / 2], dtype=complex)
    B_z_ref = complex(2 * np.linalg.solve(M, rhs)[0])  # uniform B_z inside

    # ------------------------------------------------------------------ #
    # Geometry (CSG)                                                      #
    # ------------------------------------------------------------------ #
    sph_o = Sphere(Pnt(0, 0, 0), R2)
    sph_i = Sphere(Pnt(0, 0, 0), R1)
    box   = OrthoBrick(Pnt(-L, -L, -L), Pnt(L, L, L)).bc("outer")

    shell = (sph_o - sph_i).mat("shell").maxh(0.003)   # 3 mm ≈ 2 layers / 5 mm
    cav   = sph_i.mat("air_in")
    air   = (box - sph_o).mat("air_out")

    geo = CSGeometry()
    geo.Add(shell); geo.Add(cav); geo.Add(air)
    mesh = Mesh(geo.GenerateMesh(maxh=0.05, grading=0.4))
    mesh.Curve(2)

    # ------------------------------------------------------------------ #
    # Material CFs                                                        #
    # ------------------------------------------------------------------ #
    nu       = CoefficientFunction([nu0] * len(mesh.GetMaterials()))
    sigma_cf = mesh.MaterialCF({"shell": sigma}, default=0.0)
    A0       = CoefficientFunction((-B0 * y / 2, B0 * x / 2, 0.0))

    # Scattered-field source: −jωσ A₀ in the shell (σ_cf = 0 elsewhere)
    def add_source(lf, validate_fn):
        v_A, psi = validate_fn
        lf += -1j * omega * sigma_cf * A0 * v_A          * dx
        lf += -1j * omega * sigma_cf * A0 * grad(psi)    * dx

    # ------------------------------------------------------------------ #
    # Solve (Dirichlet box, no Kelvin transform)                          #
    # ------------------------------------------------------------------ #
    gfAPhi = solve_eddy_current_harmonic_APhi(
        mesh, nu, sigma_cf, omega, add_source,
        order=2, precond="local", dirichlet="outer",
        periodic=False, tol=1e-8, maxsteps=2000, restart=200, reg=1e-10,
    )

    # B_total_z = curl(Ã)_z + B₀  (field is uniform inside the cavity)
    gfA    = gfAPhi.components[0]
    B_z_fem = complex(curl(gfA)[2](mesh(0.01, 0.0, 0.0))) + B0

    tol = 0.06   # 6 % of B₀
    assert abs(B_z_fem - B_z_ref) < tol, (
        f"TEAM-6 sphere B_z(centre) = {B_z_fem:.4f} T, "
        f"analytical = {B_z_ref:.4f} T  (tol {tol:.2f} T)"
    )


def validate_team21_p21a1_eddy_loss():
    """TEAM Workshop Problem 21 P21a-1 — same plate as P21a-0 but with one
    central horizontal slit (660 × 10 × 10 mm through-thickness slot centred in
    the plate at z = plate_z/2) that interrupts the main eddy-current loop.

    Geometry: identical coils and plate position to P21a-0.  Slit: 660 mm wide
    (y-direction), 10 mm tall (z-direction), through the full 10 mm thickness.
    Centred at y = 0, z = 229 mm (mid-height of the 458 mm plate).

    Pass criterion: |P_eddy − 3.40 W| < 1.5 W (±44 %).
    Reference: TEAM P21 Family V.2009 Table A2-4
        measured  3.40 W, calculated 3.34 W (Case I, 10 A rms, 50 Hz).
    """
    import numpy as np
    from netgen.occ import Box, Pnt, Glue, OCCGeometry
    from netgen.meshing import MeshingParameters

    mu0         = 4 * np.pi * 1e-7
    nu0         = 1.0 / mu0
    sigma_steel = 1.3889e6
    omega       = 2 * np.pi * 50.0
    turns       = 300
    I_rms       = 10.0

    half_o    = 0.124
    half_i    = 0.100
    coil_h    = 0.217
    coil_gap  = 0.024
    plate_gap = 0.012
    plate_x0  = half_o + plate_gap        # 0.136 m
    plate_x1  = plate_x0 + 0.010         # 0.146 m
    plate_y   = 0.410
    plate_z   = coil_h * 2 + coil_gap    # 0.458 m
    coil_w    = half_o - half_i          # 0.024 m
    margin    = 0.10

    # Slit: 660 mm × 10 mm (y × z), through full x-thickness, centred at z = plate_z/2
    slit_y   = 0.330    # half-width: 660/2 mm
    slit_z   = 0.005    # half-height: 10/2 mm
    ctr_z    = plate_z / 2  # 0.229 m

    coil1 = (Box(Pnt(-half_o, -half_o, 0.0),  Pnt(half_o, half_o, coil_h))
             - Box(Pnt(-half_i, -half_i, 0.0), Pnt(half_i, half_i, coil_h)))
    coil1.mat("coil1")

    z2    = coil_h + coil_gap
    coil2 = (Box(Pnt(-half_o, -half_o, z2),  Pnt(half_o, half_o, z2 + coil_h))
             - Box(Pnt(-half_i, -half_i, z2), Pnt(half_i, half_i, z2 + coil_h)))
    coil2.mat("coil2")

    plate_box = Box(Pnt(plate_x0, -plate_y, 0.0), Pnt(plate_x1, plate_y, plate_z))
    slit_box  = Box(Pnt(plate_x0, -slit_y, ctr_z - slit_z),
                    Pnt(plate_x1,  slit_y, ctr_z + slit_z))
    plate = plate_box - slit_box
    plate.faces.maxh = 0.015
    plate.mat("plate")

    air_box = Box(
        Pnt(-half_o - margin, -plate_y - margin, -margin),
        Pnt(plate_x1 + margin,  plate_y + margin, plate_z + margin))
    air_box.faces.name = "outer"
    air = air_box - coil1 - coil2 - plate
    air.mat("air")

    ngmesh = OCCGeometry(Glue([coil1, coil2, plate, air])).GenerateMesh(
        MeshingParameters(maxh=0.07, grading=0.35))
    mesh = Mesh(ngmesh)
    mesh.Curve(2)

    mats     = mesh.GetMaterials()
    nu       = CoefficientFunction([nu0] * len(mats))
    sigma_cf = mesh.MaterialCF({"plate": sigma_steel}, default=0.0)

    projx    = IfPos(x - (-half_i), IfPos(x - half_i, half_i, x), -half_i)
    projy    = IfPos(y - (-half_i), IfPos(y - half_i, half_i, y), -half_i)
    tx       = projy - y
    ty       = x - projx
    tau_coil = CoefficientFunction((tx, ty, 0)) / sqrt(tx**2 + ty**2 + 1e-30)

    NI_peak  = turns * I_rms * math.sqrt(2)
    J_scale  = NI_peak / (coil_w * coil_h)

    def add_source(lf, validate_fn):
        v, _ = validate_fn
        lf += J_scale * tau_coil * v * dx("coil1")
        lf += (-J_scale) * tau_coil * v * dx("coil2")

    gfAPhi = solve_eddy_current_harmonic_APhi(
        mesh, nu, sigma_cf, omega, add_source,
        order=2, precond="local", dirichlet="outer",
        periodic=False, tol=1e-8, maxsteps=2000, restart=200, reg=1e-10,
    )

    gfA, gfPhi = gfAPhi.components
    E_cf   = -1j * omega * (gfA + grad(gfPhi))
    P_eddy = (0.5 * Integrate(
        sigma_cf * InnerProduct(E_cf, Conj(E_cf)).real, mesh)).real

    ref = 3.40   # W, measured (P21a-1, Case I, 10 A rms, 50 Hz)
    assert abs(P_eddy - ref) < 1.5, (
        f"TEAM-21 P21a-1 eddy loss = {P_eddy:.2f} W, "
        f"ref = {ref:.2f} W  (tol +-1.5 W)"
    )


def validate_cylinder_axial_eddy():
    """Conducting cylinder (R=50 mm, σ=1.6×10⁷ S/m, μᵣ=1) in a uniform axial
    sinusoidal field B₀ = 1 T at 50 Hz.  Scattered-field A-Φ formulation.

    Geometry: finite cylinder height H = 300 mm modelled inside a Dirichlet box.
    The cylinder height is 6× the radius so the middle cross-section (z ≈ 0)
    closely matches the 2-D infinite-cylinder analytical solution.

    Analytical reference (2-D infinite cylinder, axial B₀ field):
        B_z(r) = B₀ · I₀(κr) / I₀(κR)   where κ = √(jωσμ₀) = (1+j)/δ
    Skin depth δ ≈ 17.8 mm at 50 Hz.  κR ≈ 2.81(1+j),  |I₀(κR)| ≈ 3.6.

    This is equivalent to TEAM Workshop Problem 2 (Infinitely Long Cylinder in a
    Sinusoidal Field) at similar parameters.

    Pass criterion: |B_z_fem(centre) − B_z_ref(centre)| < 0.05 T  (5 % of B₀).
    """
    import numpy as np
    from netgen.csg import CSGeometry, Cylinder, OrthoBrick, Pnt, Vec

    mu0   = 4 * np.pi * 1e-7
    nu0   = 1.0 / mu0
    sigma = 1.6e7          # S/m
    omega = 2 * np.pi * 50.0
    B0    = 1.0            # applied uniform axial field [T]
    R     = 0.050          # cylinder radius [m]
    H     = 0.300          # finite model height [m]  (6R)
    L     = 0.25           # Dirichlet box half-length [m]

    # ------------------------------------------------------------------ #
    # Analytical reference: B_z(r) = B₀ I₀(κr) / I₀(κR)               #
    # I₀ via series: I₀(z) = Σ (z/2)^(2k)/(k!)²                        #
    # ------------------------------------------------------------------ #
    kappa = np.sqrt(1j * omega * sigma * mu0)

    def _I0(z, nterms=40):
        z = complex(z)
        result, term, z24 = 1.0+0j, 1.0+0j, (z*z/4)
        for k in range(1, nterms+1):
            term *= z24 / k**2
            result += term
        return result

    B_z_ref = complex(B0 / _I0(kappa * R))   # value at r→0 centre

    # ------------------------------------------------------------------ #
    # Geometry (CSG)                                                      #
    # ------------------------------------------------------------------ #
    cyl_inf = Cylinder(Pnt(0, 0, -1), Pnt(0, 0, 1), R)
    slab    = OrthoBrick(Pnt(-2*R, -2*R, -H/2), Pnt(2*R, 2*R, H/2))
    cond    = (cyl_inf * slab).mat("conductor").maxh(0.010)
    box     = OrthoBrick(Pnt(-L, -L, -L), Pnt(L, L, L)).bc("outer")
    air     = (box - cond).mat("air")

    geo = CSGeometry()
    geo.Add(cond); geo.Add(air)
    mesh = Mesh(geo.GenerateMesh(maxh=0.05, grading=0.4))
    mesh.Curve(2)

    # ------------------------------------------------------------------ #
    # Material CFs                                                        #
    # ------------------------------------------------------------------ #
    nu       = CoefficientFunction([nu0] * len(mesh.GetMaterials()))
    sigma_cf = mesh.MaterialCF({"conductor": sigma}, default=0.0)
    A0       = CoefficientFunction((-B0 * y / 2, B0 * x / 2, 0.0))

    def add_source(lf, validate_fn):
        v_A, psi = validate_fn
        lf += -1j * omega * sigma_cf * A0 * v_A       * dx
        lf += -1j * omega * sigma_cf * A0 * grad(psi) * dx

    # ------------------------------------------------------------------ #
    # Solve                                                               #
    # ------------------------------------------------------------------ #
    gfAPhi = solve_eddy_current_harmonic_APhi(
        mesh, nu, sigma_cf, omega, add_source,
        order=2, precond="local", dirichlet="outer",
        periodic=False, tol=1e-8, maxsteps=2000, restart=200, reg=1e-10,
    )

    # B_z at cylinder centre (r ≈ 0, z = 0)
    gfA    = gfAPhi.components[0]
    B_z_fem = complex(curl(gfA)[2](mesh(0.001, 0.0, 0.0))) + B0

    tol = 0.05   # 5 % of B₀
    assert abs(B_z_fem - B_z_ref) < tol, (
        f"Cylinder axial eddy B_z(centre) = {B_z_fem:.4f} T, "
        f"analytical = {B_z_ref:.4f} T  (tol {tol:.2f} T)"
    )


def validate_sphere_magnetostatic_scattered():
    """Permeable sphere (μr = 10, R = 50 mm) in a uniform applied field B₀ = 1 T.
    Scattered-field A-form magnetostatics: the source arises from the permeability
    jump (not a coil current).

    Analytical reference (demagnetization factor N = 1/3 for a sphere):
        B_z_inside = 3 μr/(μr + 2) × B₀  (uniform inside)
        For μr = 10: B_z = 30/12 × 1 T = 2.5 T

    This validates the ∇×(ν∇×A) = −∇×(νB₀) scattered-field source term for
    non-trivial ν(r) (spatially varying permeability).

    Pass criterion: |B_z_fem − 2.5 T| < 0.15 T  (6 % of B_z_inside).
    """
    import numpy as np
    from ngsolve import (HCurl, BilinearForm, LinearForm, GridFunction,
                         Preconditioner)
    from netgen.csg import CSGeometry, Sphere, OrthoBrick, Pnt

    mu0   = 4 * np.pi * 1e-7
    nu0   = 1.0 / mu0
    mu_r  = 10.0
    nu_sph = nu0 / mu_r
    B0    = 1.0   # [T]
    R     = 0.05  # sphere radius [m]
    L     = 0.25  # Dirichlet box half-length [m]

    # Analytical reference (exact for uniform applied field + sphere)
    B_z_ref = 3 * mu_r / (mu_r + 2) * B0   # = 2.5 T for μr=10

    # -- Geometry --
    sph = Sphere(Pnt(0, 0, 0), R)
    box = OrthoBrick(Pnt(-L, -L, -L), Pnt(L, L, L)).bc("outer")
    sph_reg = sph.mat("sphere").maxh(0.008)
    air_reg = (box - sph).mat("air")
    geo = CSGeometry()
    geo.Add(sph_reg); geo.Add(air_reg)
    mesh = Mesh(geo.GenerateMesh(maxh=0.04, grading=0.4))
    mesh.Curve(2)

    # ν coefficient: ν₀/μr inside sphere, ν₀ outside
    nu_cf   = mesh.MaterialCF({"sphere": nu_sph}, default=nu0)
    # Perturbation ν: positive source comes from (ν₀ − ν_sphere) = ν₀(1 − 1/μr)
    nu_pert = mesh.MaterialCF({"sphere": nu0 * (1.0 - 1.0/mu_r)}, default=0.0)
    B0_vec  = CoefficientFunction((0.0, 0.0, B0))

    # -- Scattered-field A-form: ∇×(ν∇×Ã) = −∇×(νB₀) = +nu_pert × B₀ (source) --
    fes = HCurl(mesh, order=2, dirichlet="outer")
    u, v = fes.TnT()
    a = BilinearForm(fes)
    a += nu_cf * InnerProduct(curl(u), curl(v)) * dx + 1e-8 * nu0 * InnerProduct(u, v) * dx
    f = LinearForm(fes)
    f += InnerProduct(nu_pert * B0_vec, curl(v)) * dx
    a.Assemble(); f.Assemble()
    gfA = GridFunction(fes)
    gfA.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * f.vec

    # B_z_total at sphere centre = scattered B_z + B₀
    B_z_fem = float(curl(gfA)[2](mesh(0.005, 0.0, 0.0))) + B0

    tol = 0.15   # 6 % of 2.5 T
    assert abs(B_z_fem - B_z_ref) < tol, (
        f"Sphere magnetostatic B_z(centre) = {B_z_fem:.4f} T, "
        f"analytical = {B_z_ref:.4f} T  (tol {tol:.2f} T)"
    )


def validate_team21_p21b2n_eddy_loss():
    """TEAM Workshop Problem 21 P21b-2N — two non-magnetic steel plates (each
    458 × 248 × 10 mm, 20Mn23Al) separated by a 200 mm gap centred on the
    coil axis in y.  Identical coils to P21a family; Case I excitation.

    Geometry: plates side-by-side in y with 200 mm gap at y = 0.
      Plate 1: y from −348 mm to −100 mm (outer edge to inner edge)
      Plate 2: y from +100 mm to +348 mm
    The gap allows the coil flux to pass through without inducing eddy currents
    in the gap region; each plate is much narrower than the P21a plate (248 vs 820 mm).

    Pass criterion: |P_eddy − 1.38 W| < 0.6 W (±43 %).
    Reference: TEAM P21 Family V.2009 Table A2-4
        measured  1.38 W, calculated 1.37 W (Case I, 10 A rms, 50 Hz).
    """
    import numpy as np
    from netgen.occ import Box, Pnt, Glue, OCCGeometry
    from netgen.meshing import MeshingParameters

    mu0         = 4 * np.pi * 1e-7
    nu0         = 1.0 / mu0
    sigma_steel = 1.3889e6
    omega       = 2 * np.pi * 50.0
    turns       = 300
    I_rms       = 10.0

    half_o    = 0.124
    half_i    = 0.100
    coil_h    = 0.217
    coil_gap  = 0.024
    plate_gap = 0.012
    plate_x0  = half_o + plate_gap        # 0.136 m
    plate_x1  = plate_x0 + 0.010         # 0.146 m
    plate_z   = coil_h * 2 + coil_gap    # 0.458 m
    coil_w    = half_o - half_i
    margin    = 0.10

    # P21b-2N: two 458×248×10 mm plates, 200 mm gap (centered at y=0)
    half_gap  = 0.100   # = 200 mm / 2
    plate_w   = 0.248   # each plate 248 mm in y

    coil1 = (Box(Pnt(-half_o, -half_o, 0.0),  Pnt(half_o, half_o, coil_h))
             - Box(Pnt(-half_i, -half_i, 0.0), Pnt(half_i, half_i, coil_h)))
    coil1.mat("coil1")

    z2    = coil_h + coil_gap
    coil2 = (Box(Pnt(-half_o, -half_o, z2),  Pnt(half_o, half_o, z2 + coil_h))
             - Box(Pnt(-half_i, -half_i, z2), Pnt(half_i, half_i, z2 + coil_h)))
    coil2.mat("coil2")

    # Plate 1 (y < 0 side)  and  Plate 2 (y > 0 side)
    plate1 = Box(Pnt(plate_x0, -(half_gap + plate_w), 0.0),
                 Pnt(plate_x1, -half_gap,              plate_z))
    plate1.faces.maxh = 0.025
    plate1.mat("plate")

    plate2 = Box(Pnt(plate_x0,  half_gap,              0.0),
                 Pnt(plate_x1,  half_gap + plate_w,    plate_z))
    plate2.faces.maxh = 0.025
    plate2.mat("plate")

    # Air box large enough to contain both plates
    plate_y_max = half_gap + plate_w   # 0.348 m
    air_box = Box(
        Pnt(-half_o - margin, -(plate_y_max + margin), -margin),
        Pnt(plate_x1 + margin,  plate_y_max + margin,  plate_z + margin))
    air_box.faces.name = "outer"
    air = air_box - coil1 - coil2 - plate1 - plate2
    air.mat("air")

    ngmesh = OCCGeometry(Glue([coil1, coil2, plate1, plate2, air])).GenerateMesh(
        MeshingParameters(maxh=0.07, grading=0.35))
    mesh = Mesh(ngmesh)
    mesh.Curve(2)

    mats     = mesh.GetMaterials()
    nu       = CoefficientFunction([nu0] * len(mats))
    sigma_cf = mesh.MaterialCF({"plate": sigma_steel}, default=0.0)

    projx    = IfPos(x - (-half_i), IfPos(x - half_i, half_i, x), -half_i)
    projy    = IfPos(y - (-half_i), IfPos(y - half_i, half_i, y), -half_i)
    tx       = projy - y
    ty       = x - projx
    tau_coil = CoefficientFunction((tx, ty, 0)) / sqrt(tx**2 + ty**2 + 1e-30)

    NI_peak  = turns * I_rms * math.sqrt(2)
    J_scale  = NI_peak / (coil_w * coil_h)

    def add_source(lf, validate_fn):
        v, _ = validate_fn
        lf += J_scale * tau_coil * v * dx("coil1")
        lf += (-J_scale) * tau_coil * v * dx("coil2")

    gfAPhi = solve_eddy_current_harmonic_APhi(
        mesh, nu, sigma_cf, omega, add_source,
        order=2, precond="local", dirichlet="outer",
        periodic=False, tol=1e-8, maxsteps=2000, restart=200, reg=1e-10,
    )

    gfA, gfPhi = gfAPhi.components
    E_cf   = -1j * omega * (gfA + grad(gfPhi))
    P_eddy = (0.5 * Integrate(
        sigma_cf * InnerProduct(E_cf, Conj(E_cf)).real, mesh)).real

    ref = 1.38   # W, measured (P21b-2N, Case I, 10 A rms, 50 Hz)
    assert abs(P_eddy - ref) < 0.6, (
        f"TEAM-21 P21b-2N eddy loss = {P_eddy:.2f} W, "
        f"ref = {ref:.2f} W  (tol +-0.6 W)"
    )


def validate_conducting_sphere_axial_eddy():
    """Conducting sphere (σ=1.6e7 S/m, R=50 mm) in uniform axial AC field B₀=1 T, 50 Hz.
    Scattered-field A-Φ formulation: same source as the cylinder test.

    Analytical reference (modified spherical Bessel, exact for r→0):
        B_z(0) = B₀ · p / sinh(p)
        where p = κR = (1+j) R/δ,  δ = skin depth.
        For σ=1.6e7, f=50 Hz: δ = 17.8 mm → R/δ = 2.81, p = 2.81(1+j).
        B_z_ref ≈ -0.209 - 0.432j T  (amplitude ≈ 0.480 T, phase ≈ -116°).

    Pass criterion: |B_z_fem − B_z_ref| < 0.05 T  (≈ 5 % of B₀).
    """
    import numpy as np

    mu0   = MU0
    nu0   = 1.0 / mu0
    sigma = 1.6e7   # S/m (same material as cylinder test)
    omega = 2 * np.pi * 50.0
    B0    = 1.0
    R     = 0.050   # sphere radius [m]
    L     = 0.25    # Dirichlet box half-length [m]

    # -- Analytical reference: B_z(center) = B0*p/sinh(p) --
    kappa   = np.sqrt(1j * omega * sigma * mu0)
    p       = kappa * R
    B_z_ref = complex(B0 * p / np.sinh(p))

    # -- Geometry (CSG) --
    sph  = Sphere(Pnt(0, 0, 0), R)
    box  = OrthoBrick(Pnt(-L, -L, -L), Pnt(L, L, L)).bc("outer")
    cond = sph.mat("conductor").maxh(0.006)
    air  = (box - sph).mat("air")
    geo  = CSGeometry()
    geo.Add(cond); geo.Add(air)
    mesh = Mesh(geo.GenerateMesh(maxh=0.05, grading=0.4))
    mesh.Curve(2)

    # -- Material CFs --
    nu       = CoefficientFunction([nu0] * len(mesh.GetMaterials()))
    sigma_cf = mesh.MaterialCF({"conductor": sigma}, default=0.0)
    A0       = CoefficientFunction((-B0 * y / 2, B0 * x / 2, 0.0))

    def add_source(lf, validate_fn):
        v_A, psi = validate_fn
        lf += -1j * omega * sigma_cf * A0 * v_A       * dx
        lf += -1j * omega * sigma_cf * A0 * grad(psi) * dx

    # -- Solve --
    from ngsolve import TaskManager
    with TaskManager():
        gfAPhi = solve_eddy_current_harmonic_APhi(
            mesh, nu, sigma_cf, omega, add_source,
            order=2, precond="local", dirichlet="outer",
            periodic=False, tol=1e-8, maxsteps=2000, restart=200, reg=1e-10,
        )

    gfA     = gfAPhi.components[0]
    B_z_fem = complex(curl(gfA)[2](mesh(0.001, 0.0, 0.0))) + B0

    tol = 0.05   # 5 % of B₀
    assert abs(B_z_fem - B_z_ref) < tol, (
        f"Conducting sphere axial eddy B_z(centre) = {B_z_fem:.4f} T, "
        f"analytical = {B_z_ref:.4f} T  (tol {tol:.2f} T)"
    )


def validate_permeable_cylinder_transverse_static():
    """Permeable cylinder (μr=10, R=50 mm, L=500 mm) in transverse uniform field B₀=1 T.
    Scattered-field A-form magnetostatics.

    Analytical reference (infinite-cylinder limit, demagnetization N=1/2):
        B_x_inside = 2μr / (μr+1) × B₀  =  20/11 ≈ 1.818 T

    Pass criterion: |B_x_fem − B_x_ref| < 0.15 T  (accommodates finite-length end effects).
    """
    mu_r = 10.0
    B0   = 1.0
    R    = 0.05    # cylinder radius [m]
    H    = 0.25    # half-height (L=0.5 m, L/R=10)
    Lb   = 0.35    # Dirichlet box half-side [m]

    nu0     = 1.0 / MU0
    B_x_ref = 2.0 * mu_r / (mu_r + 1.0) * B0   # 20/11 ≈ 1.818 T

    # -- Geometry (finite cylinder, axis along z) --
    cyl = Cylinder(Pnt(0, 0, -H), Pnt(0, 0, H), R)
    box = OrthoBrick(Pnt(-Lb, -Lb, -Lb), Pnt(Lb, Lb, Lb)).bc("outer")
    cyl_reg = cyl.mat("cylinder").maxh(0.012)
    air_reg = (box - cyl).mat("air")
    geo = CSGeometry()
    geo.Add(cyl_reg); geo.Add(air_reg)
    mesh = Mesh(geo.GenerateMesh(maxh=0.05, grading=0.4))
    mesh.Curve(2)

    nu_cf   = mesh.MaterialCF({"cylinder": nu0 / mu_r}, default=nu0)
    nu_pert = mesh.MaterialCF({"cylinder": nu0 * (1.0 - 1.0/mu_r)}, default=0.0)
    B0_vec  = CoefficientFunction((B0, 0.0, 0.0))   # B₀ in x-direction

    fes = HCurl(mesh, order=2, dirichlet="outer")
    u, v = fes.TnT()
    a = BilinearForm(fes)
    a += nu_cf * InnerProduct(curl(u), curl(v)) * dx + 1e-8 * NU0 * InnerProduct(u, v) * dx
    f = LinearForm(fes)
    f += InnerProduct(nu_pert * B0_vec, curl(v)) * dx
    a.Assemble(); f.Assemble()

    gfA = GridFunction(fes)
    gfA.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * f.vec

    B_x_fem = float(curl(gfA)[0](mesh(0.001, 0.0, 0.0))) + B0

    tol = 0.15
    assert abs(B_x_fem - B_x_ref) < tol, (
        f"Permeable cylinder transverse B_x(centre) = {B_x_fem:.4f} T, "
        f"analytical = {B_x_ref:.4f} T  (tol {tol:.2f} T)"
    )


def validate_solenoid_center_bz_static():
    """Thick annular solenoid (N=1000 turns, I=1 A, Ri=40 mm, Ro=60 mm, L=200 mm).
    Azimuthal current source J_φ = NI/((Ro-Ri)L) in the winding volume.

    Analytical reference (Biot-Savart, rectangular winding cross-section):
        Bz(centre) = (μ₀J₀/2) × 2h × (arcsinh(Ro/h) − arcsinh(Ri/h))
        where h = L/2.  For these parameters: Bz_ref ≈ 5.62 mT.

    Pass criterion: |Bz_fem − Bz_ref| < 0.30 mT  (≈ 5 % of Bz_ref).
    """
    import numpy as np

    N_turns = 1000
    I_amp   = 1.0
    Ri, Ro  = 0.040, 0.060   # inner / outer radii [m]
    zbot, ztop = 0.0, 0.20   # winding extent [m]
    zc = (zbot + ztop) / 2   # 0.10 m

    h   = (ztop - zbot) / 2
    J0  = N_turns * I_amp / ((Ro - Ri) * (ztop - zbot))  # 250 000 A/m²
    Bz_ref = MU0 * J0 * h * (np.arcsinh(Ro / h) - np.arcsinh(Ri / h))

    Lb = 0.35

    # -- Geometry (annular coil + inner air + outer air) --
    cyl_out = Cylinder(Pnt(0, 0, zbot), Pnt(0, 0, ztop), Ro)
    cyl_in  = Cylinder(Pnt(0, 0, zbot), Pnt(0, 0, ztop), Ri)
    box     = OrthoBrick(Pnt(-Lb, -Lb, -Lb), Pnt(Lb, Lb, Lb)).bc("outer")
    coil_reg = (cyl_out - cyl_in).mat("coil").maxh(0.010)
    core_reg = cyl_in.mat("air")
    air_reg  = (box - cyl_out).mat("air")
    geo = CSGeometry()
    geo.Add(coil_reg); geo.Add(core_reg); geo.Add(air_reg)
    mesh = Mesh(geo.GenerateMesh(maxh=0.05, grading=0.35))
    mesh.Curve(2)

    nu0 = 1.0 / MU0
    fes = HCurl(mesh, order=2, dirichlet="outer")
    u, v = fes.TnT()
    a = BilinearForm(fes)
    a += NU0 * InnerProduct(curl(u), curl(v)) * dx + 1e-8 * NU0 * InnerProduct(u, v) * dx

    r_cf  = sqrt(x**2 + y**2 + 1e-30)
    J_src = J0 * CoefficientFunction((-y / r_cf, x / r_cf, 0.0))
    f = LinearForm(fes)
    f += InnerProduct(J_src, v) * dx(definedon=mesh.Materials("coil"))
    a.Assemble(); f.Assemble()

    gfA = GridFunction(fes)
    gfA.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * f.vec

    Bz_fem = float(curl(gfA)[2](mesh(0.001, 0.001, zc)))

    tol = 0.30e-3   # 0.30 mT ≈ 5 % of Bz_ref
    assert abs(Bz_fem - Bz_ref) < tol, (
        f"Solenoid centre Bz = {Bz_fem*1e3:.4f} mT, "
        f"analytical = {Bz_ref*1e3:.4f} mT  (tol {tol*1e3:.2f} mT)"
    )


def validate_team21_p21a2_eddy_loss():
    """TEAM Workshop Problem 21 P21a-2 — same plate as P21a-0 but with two
    horizontal slits (each 660 × 10 × 10 mm) placed at z = plate_z/3 and
    z = 2*plate_z/3, dividing the plate into 3 equal-height strips.

    Slit positions are estimated by symmetry (the original PDF is password-protected).
    Wider tolerance (±1.2 W) covers geometry uncertainty.

    Pass criterion: |P_eddy − 1.68 W| < 1.2 W.
    Reference: TEAM P21 Family V.2009 Table A2-4
        measured  1.68 W, calculated 1.63 W (Case I, 10 A rms, 50 Hz).
    """
    import numpy as np
    from netgen.occ import Box, Pnt, Glue, OCCGeometry
    from netgen.meshing import MeshingParameters

    nu0         = 1.0 / MU0
    sigma_steel = 1.3889e6
    omega       = 2 * np.pi * 50.0
    turns       = 300
    I_rms       = 10.0

    half_o    = 0.124; half_i = 0.100
    coil_h    = 0.217; coil_gap = 0.024
    plate_gap = 0.012
    plate_x0  = half_o + plate_gap
    plate_x1  = plate_x0 + 0.010
    plate_y   = 0.410
    plate_z   = coil_h * 2 + coil_gap   # 0.458 m
    coil_w    = half_o - half_i
    margin    = 0.10

    slit_y = 0.330; slit_hz = 0.005   # each slit: 660 mm wide, 10 mm tall
    z1 = plate_z / 3                  # 0.1527 m
    z2 = 2.0 * plate_z / 3           # 0.3053 m

    coil1 = (Box(Pnt(-half_o, -half_o, 0.0), Pnt(half_o, half_o, coil_h))
             - Box(Pnt(-half_i, -half_i, 0.0), Pnt(half_i, half_i, coil_h)))
    coil1.mat("coil1")
    z2c   = coil_h + coil_gap
    coil2 = (Box(Pnt(-half_o, -half_o, z2c), Pnt(half_o, half_o, z2c + coil_h))
             - Box(Pnt(-half_i, -half_i, z2c), Pnt(half_i, half_i, z2c + coil_h)))
    coil2.mat("coil2")

    plate_box = Box(Pnt(plate_x0, -plate_y, 0.0), Pnt(plate_x1, plate_y, plate_z))
    slit1 = Box(Pnt(plate_x0, -slit_y, z1 - slit_hz), Pnt(plate_x1, slit_y, z1 + slit_hz))
    slit2 = Box(Pnt(plate_x0, -slit_y, z2 - slit_hz), Pnt(plate_x1, slit_y, z2 + slit_hz))
    plate = plate_box - slit1 - slit2
    plate.faces.maxh = 0.015
    plate.mat("plate")

    air_box = Box(
        Pnt(-half_o - margin, -plate_y - margin, -margin),
        Pnt(plate_x1 + margin,  plate_y + margin,  plate_z + margin))
    air_box.faces.name = "outer"
    air = air_box - coil1 - coil2 - plate
    air.mat("air")

    ngmesh = OCCGeometry(Glue([coil1, coil2, plate, air])).GenerateMesh(
        MeshingParameters(maxh=0.07, grading=0.35))
    mesh = Mesh(ngmesh); mesh.Curve(2)

    mats     = mesh.GetMaterials()
    nu       = CoefficientFunction([nu0] * len(mats))
    sigma_cf = mesh.MaterialCF({"plate": sigma_steel}, default=0.0)

    projx = IfPos(x - (-half_i), IfPos(x - half_i, half_i, x), -half_i)
    projy = IfPos(y - (-half_i), IfPos(y - half_i, half_i, y), -half_i)
    tx = projy - y; ty = x - projx
    tau_coil = CoefficientFunction((tx, ty, 0)) / sqrt(tx**2 + ty**2 + 1e-30)
    NI_peak  = turns * I_rms * math.sqrt(2)
    J_scale  = NI_peak / (coil_w * coil_h)

    def add_source(lf, validate_fn):
        v, _ = validate_fn
        lf += J_scale * tau_coil * v * dx("coil1")
        lf += (-J_scale) * tau_coil * v * dx("coil2")

    from ngsolve import TaskManager
    with TaskManager():
        gfAPhi = solve_eddy_current_harmonic_APhi(
            mesh, nu, sigma_cf, omega, add_source,
            order=2, precond="local", dirichlet="outer",
            periodic=False, tol=1e-8, maxsteps=2000, restart=200, reg=1e-10)

    gfA, gfPhi = gfAPhi.components
    E_cf   = -1j * omega * (gfA + grad(gfPhi))
    P_eddy = (0.5 * Integrate(
        sigma_cf * InnerProduct(E_cf, Conj(E_cf)).real, mesh)).real

    ref = 1.68
    assert abs(P_eddy - ref) < 1.2, (
        f"TEAM-21 P21a-2 eddy loss = {P_eddy:.2f} W, "
        f"ref = {ref:.2f} W  (tol +-1.2 W, slit positions estimated)")


def validate_team21_p21a3_eddy_loss():
    """TEAM Workshop Problem 21 P21a-3 — same plate with four horizontal slits
    (each 660 × 10 × 10 mm) at z = plate_z × k/5 for k=1..4, creating 5 equal strips.

    Slit positions estimated by symmetry (PDF password-protected).

    Pass criterion: |P_eddy − 1.25 W| < 1.0 W.
    Reference: TEAM P21 Family V.2009 Table A2-4
        measured  1.25 W, calculated 1.22 W (Case I, 10 A rms, 50 Hz).
    """
    import numpy as np
    from netgen.occ import Box, Pnt, Glue, OCCGeometry
    from netgen.meshing import MeshingParameters

    nu0         = 1.0 / MU0
    sigma_steel = 1.3889e6
    omega       = 2 * np.pi * 50.0
    turns       = 300
    I_rms       = 10.0

    half_o    = 0.124; half_i = 0.100
    coil_h    = 0.217; coil_gap = 0.024
    plate_gap = 0.012
    plate_x0  = half_o + plate_gap
    plate_x1  = plate_x0 + 0.010
    plate_y   = 0.410
    plate_z   = coil_h * 2 + coil_gap   # 0.458 m
    coil_w    = half_o - half_i
    margin    = 0.10

    slit_y = 0.330; slit_hz = 0.005

    coil1 = (Box(Pnt(-half_o, -half_o, 0.0), Pnt(half_o, half_o, coil_h))
             - Box(Pnt(-half_i, -half_i, 0.0), Pnt(half_i, half_i, coil_h)))
    coil1.mat("coil1")
    z2c   = coil_h + coil_gap
    coil2 = (Box(Pnt(-half_o, -half_o, z2c), Pnt(half_o, half_o, z2c + coil_h))
             - Box(Pnt(-half_i, -half_i, z2c), Pnt(half_i, half_i, z2c + coil_h)))
    coil2.mat("coil2")

    plate = Box(Pnt(plate_x0, -plate_y, 0.0), Pnt(plate_x1, plate_y, plate_z))
    for k in range(1, 5):
        zk = k * plate_z / 5.0
        slit = Box(Pnt(plate_x0, -slit_y, zk - slit_hz), Pnt(plate_x1, slit_y, zk + slit_hz))
        plate = plate - slit
    plate.faces.maxh = 0.015
    plate.mat("plate")

    air_box = Box(
        Pnt(-half_o - margin, -plate_y - margin, -margin),
        Pnt(plate_x1 + margin,  plate_y + margin,  plate_z + margin))
    air_box.faces.name = "outer"
    air = air_box - coil1 - coil2 - plate
    air.mat("air")

    ngmesh = OCCGeometry(Glue([coil1, coil2, plate, air])).GenerateMesh(
        MeshingParameters(maxh=0.07, grading=0.35))
    mesh = Mesh(ngmesh); mesh.Curve(2)

    mats     = mesh.GetMaterials()
    nu       = CoefficientFunction([nu0] * len(mats))
    sigma_cf = mesh.MaterialCF({"plate": sigma_steel}, default=0.0)

    projx = IfPos(x - (-half_i), IfPos(x - half_i, half_i, x), -half_i)
    projy = IfPos(y - (-half_i), IfPos(y - half_i, half_i, y), -half_i)
    tx = projy - y; ty = x - projx
    tau_coil = CoefficientFunction((tx, ty, 0)) / sqrt(tx**2 + ty**2 + 1e-30)
    NI_peak  = turns * I_rms * math.sqrt(2)
    J_scale  = NI_peak / (coil_w * coil_h)

    def add_source(lf, validate_fn):
        v, _ = validate_fn
        lf += J_scale * tau_coil * v * dx("coil1")
        lf += (-J_scale) * tau_coil * v * dx("coil2")

    from ngsolve import TaskManager
    with TaskManager():
        gfAPhi = solve_eddy_current_harmonic_APhi(
            mesh, nu, sigma_cf, omega, add_source,
            order=2, precond="local", dirichlet="outer",
            periodic=False, tol=1e-8, maxsteps=2000, restart=200, reg=1e-10)

    gfA, gfPhi = gfAPhi.components
    E_cf   = -1j * omega * (gfA + grad(gfPhi))
    P_eddy = (0.5 * Integrate(
        sigma_cf * InnerProduct(E_cf, Conj(E_cf)).real, mesh)).real

    ref = 1.25
    assert abs(P_eddy - ref) < 1.0, (
        f"TEAM-21 P21a-3 eddy loss = {P_eddy:.2f} W, "
        f"ref = {ref:.2f} W  (tol +-1.0 W, slit positions estimated)"
    )


def validate_eddy_conducting_slab_tangential():
    """1-D eddy shielding: conducting slab (σ=1.6e7 S/m, 2h=30 mm normal to y)
    in tangential harmonic B₀=1 T (z-dir), f=50 Hz.

    Analytical (1D diffusion in y):
        B_z(0) = B₀ / cosh(p)  where  p = κh,  κ = sqrt(j·ω·σ·μ₀)
    δ ≈ 17.8 mm → h/δ ≈ 0.84 → moderate shielding.
    B_z_ref ≈ 0.6844-0.5280j T, |B|≈0.864 T.

    Pass criterion: |B_z_fem − B_z_ref| < 0.05 T.
    """
    import numpy as np
    from ngsolve import TaskManager

    mu0   = MU0
    sigma = 1.6e7
    omega = 2 * math.pi * 50.0
    B0    = 1.0
    h     = 0.015    # slab half-thickness [m]
    W     = 0.20     # plate half-width in x and z [m]
    L     = 0.35     # Dirichlet box half-side [m]

    kappa = np.sqrt(1j * omega * sigma * mu0)
    p     = kappa * h
    B_z_ref = complex(B0 / np.cosh(p))

    slab_solid = OrthoBrick(Pnt(-W, -h, -W), Pnt(W, h, W))
    box        = OrthoBrick(Pnt(-L, -L, -L), Pnt(L, L, L)).bc("outer")
    geo = CSGeometry()
    geo.Add(slab_solid.mat("slab").maxh(0.005)); geo.Add((box - slab_solid).mat("air"))
    mesh = Mesh(geo.GenerateMesh(maxh=0.06, grading=0.5))
    mesh.Curve(2)

    nu0      = 1.0 / mu0
    nu_cf    = CoefficientFunction([nu0] * len(mesh.GetMaterials()))
    sigma_cf = mesh.MaterialCF({"slab": sigma}, default=0.0)
    A0       = CoefficientFunction((-B0 * y / 2, B0 * x / 2, 0.0))

    def add_source(lf, validate_fn):
        v_A, psi = validate_fn
        lf += -1j * omega * sigma_cf * A0 * v_A       * dx
        lf += -1j * omega * sigma_cf * A0 * grad(psi) * dx

    with TaskManager():
        gfAPhi = solve_eddy_current_harmonic_APhi(
            mesh, nu_cf, sigma_cf, omega, add_source,
            order=2, precond="local", dirichlet="outer",
            periodic=False, tol=1e-8, maxsteps=2000, restart=200, reg=1e-10)

    gfA     = gfAPhi.components[0]
    B_z_fem = complex(curl(gfA)[2](mesh(0.001, 0.001, 0.001))) + B0

    tol = 0.05
    assert abs(B_z_fem - B_z_ref) < tol, (
        f"Slab eddy B_z(0) = {B_z_fem:.5f}  ref = {B_z_ref:.5f}  "
        f"err = {abs(B_z_fem-B_z_ref):.4f} T  (tol {tol} T)"
    )


def validate_helmholtz_pair_bz_center():
    """Helmholtz coil pair — two identical rings (R=50 mm) separated by d=R,
    each N=100 turns at I=100 A — B_z at the geometric centre.

    Analytical (thin-ring Biot-Savart, two loops at z=±R/2):
        B_z(0) = μ₀NI × 8 / (5√5 × R)  ≈ 179.9 mT.

    Pass criterion: |B_z_fem − B_z_ref| < 2.0 mT.
    """
    R = 50 * MM; dR = 1 * MM; hw = 1 * MM
    N = 100; I = 100.0
    B_z_ref = MU0 * N * I * 8.0 / (5.0 * math.sqrt(5.0) * R)

    ring1 = _ring(R, dR,  R / 2, hw).mat("coil").maxh(4 * MM)
    ring2 = _ring(R, dR, -R / 2, hw).mat("coil").maxh(4 * MM)
    box   = OrthoBrick(Pnt(-200 * MM, -200 * MM, -200 * MM),
                       Pnt( 200 * MM,  200 * MM,  200 * MM)).bc("outer")
    air   = (box - ring1 - ring2).mat("air")
    geo   = CSGeometry()
    geo.Add(ring1); geo.Add(ring2); geo.Add(air)
    mesh = Mesh(geo.GenerateMesh(maxh=0.04, grading=0.4))
    mesh.Curve(2)

    J0  = N * I / ((2 * dR) * (2 * hw))
    gfA = _solve_Aform(mesh, NU0, J0)
    Bz_fem = float(curl(gfA)[2](mesh(0.001, 0.001, 0.001)))

    tol = 2e-3
    assert abs(Bz_fem - B_z_ref) < tol, (
        f"Helmholtz B_z = {Bz_fem * 1e3:.3f} mT  ref = {B_z_ref * 1e3:.3f} mT  "
        f"(tol {tol * 1e3:.1f} mT)"
    )


def validate_eddy_solid_cylinder_axial():
    """Solid conducting cylinder (axis ∥ z, σ=1.6e7 S/m, R=25 mm, L=200 mm)
    in uniform harmonic B₀=1 T (z-dir), f=50 Hz.

    Analytical (1D cylindrical diffusion):
        B_z(0) = B₀ / I₀(κR)
    where I₀ is the modified Bessel function of order 0,
    κ = √(j·ω·σ·μ₀), κR ≈ (1+j)·1.405.

    Pass criterion: |B_z_fem − B_z_ref| < 0.08 T (≈ 10 % of |B_z_ref|=0.817 T).
    """
    import numpy as np
    from ngsolve import TaskManager

    mu0   = MU0
    sigma = 1.6e7
    omega = 2 * math.pi * 50.0
    B0    = 1.0
    R     = 0.025    # cylinder radius [m]
    L     = 0.20     # cylinder half-length [m]  (L/R = 8)
    Lb    = 0.35     # box half-side [m]

    kappa = np.sqrt(1j * omega * sigma * mu0)

    def _i0c(z, terms=80):
        """I₀(z) for complex z via power series."""
        result = 1.0 + 0j
        w = (z / 2) ** 2
        power = 1.0 + 0j
        fact_sq = 1
        for k in range(1, terms):
            power *= w
            fact_sq *= k * k
            c = power / fact_sq
            result += c
            if abs(c) < 1e-15 * (abs(result) + 1e-30):
                break
        return result

    B_z_ref = complex(B0 / _i0c(kappa * R))

    cyl = Cylinder(Pnt(0, 0, -L), Pnt(0, 0, L), R)
    box = OrthoBrick(Pnt(-Lb, -Lb, -Lb), Pnt(Lb, Lb, Lb)).bc("outer")
    geo = CSGeometry()
    geo.Add(cyl.mat("cyl").maxh(0.006)); geo.Add((box - cyl).mat("air"))
    mesh = Mesh(geo.GenerateMesh(maxh=0.05, grading=0.5))
    mesh.Curve(2)

    nu0      = 1.0 / mu0
    nu_cf    = CoefficientFunction([nu0] * len(mesh.GetMaterials()))
    sigma_cf = mesh.MaterialCF({"cyl": sigma}, default=0.0)
    A0       = CoefficientFunction((-B0 * y / 2, B0 * x / 2, 0.0))

    def add_source(lf, validate_fn):
        v_A, psi = validate_fn
        lf += -1j * omega * sigma_cf * A0 * v_A       * dx
        lf += -1j * omega * sigma_cf * A0 * grad(psi) * dx

    with TaskManager():
        gfAPhi = solve_eddy_current_harmonic_APhi(
            mesh, nu_cf, sigma_cf, omega, add_source,
            order=2, precond="local", dirichlet="outer",
            periodic=False, tol=1e-8, maxsteps=2000, restart=200, reg=1e-10)

    gfA     = gfAPhi.components[0]
    B_z_fem = complex(curl(gfA)[2](mesh(0.001, 0.001, 0.001))) + B0

    tol = 0.08
    assert abs(B_z_fem - B_z_ref) < tol, (
        f"Cylinder axial eddy B_z(0) = {B_z_fem:.5f}  ref = {B_z_ref:.5f}  "
        f"err = {abs(B_z_fem - B_z_ref):.4f} T  (tol {tol} T)"
    )


def validate_permeable_sphere_highmu_scattered():
    """Permeable sphere at EXTREME contrast (mu_r=1000, R=50 mm) in uniform B0=1 T.
    Exact (demagnetization N=1/3): B_inside = 3 mu_r/(mu_r+2) * B0 = 2.994 T.
    Stress-tests the scattered-field A-form at high permeability (the regime where
    a naive nested-shell mesh zeroed the interior B, issue #25).

    Pass criterion: |B_z_fem - 2.994 T| < 0.12 T  (4 %).
    """
    mu_r = 1000.0; B0 = 1.0; R = 50 * MM; L = 0.25
    B_z_ref = 3.0 * mu_r / (mu_r + 2.0) * B0
    sph = Sphere(Pnt(0, 0, 0), R)
    box = OrthoBrick(Pnt(-L, -L, -L), Pnt(L, L, L)).bc("outer")
    geo = CSGeometry()
    geo.Add(sph.mat("sphere").maxh(6 * MM)); geo.Add((box - sph).mat("air"))
    mesh = Mesh(geo.GenerateMesh(maxh=0.04, grading=0.4)); mesh.Curve(2)

    nu_cf   = mesh.MaterialCF({"sphere": NU0 / mu_r}, default=NU0)
    nu_pert = mesh.MaterialCF({"sphere": NU0 * (1.0 - 1.0 / mu_r)}, default=0.0)
    B0_vec  = CoefficientFunction((0.0, 0.0, B0))
    fes = HCurl(mesh, order=2, dirichlet="outer")
    u, v = fes.TnT()
    a = BilinearForm(fes)
    a += nu_cf * InnerProduct(curl(u), curl(v)) * dx + 1e-8 * NU0 * InnerProduct(u, v) * dx
    f = LinearForm(fes)
    f += InnerProduct(nu_pert * B0_vec, curl(v)) * dx
    a.Assemble(); f.Assemble()
    gfA = GridFunction(fes)
    gfA.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * f.vec

    B_z = float(curl(gfA)[2](mesh(0.003, 0.0, 0.0))) + B0
    assert abs(B_z - B_z_ref) < 0.12, (
        f"high-mu sphere B_z(centre) = {B_z:.4f} T, exact = {B_z_ref:.4f} T")


def validate_coaxial_loops_mutual_inductance():
    """Mutual inductance of two coaxial circular loops (a1=a2=50 mm, gap d=30 mm).
    NGSolve: solve loop-1, M = 2*pi*a2 * A_phi(a2,0,z2) / I1 (flux linkage via A).
    Exact (Maxwell): M = mu0*sqrt(a1 a2)[(2/k - k)K(k) - (2/k)E(k)],
        k^2 = 4 a1 a2 / ((a1+a2)^2 + d^2).   M_ref ~ 46.5 nH.

    Pass criterion: |M_fem - M_ref| / M_ref < 8 %.
    """
    from scipy.special import ellipk, ellipe
    a1 = 50 * MM; a2 = 50 * MM; d = 30 * MM; I1 = 100.0
    dR = 1 * MM; hw = 1 * MM
    k2 = 4 * a1 * a2 / ((a1 + a2)**2 + d**2); k = math.sqrt(k2)
    M_ref = MU0 * math.sqrt(a1 * a2) * ((2 / k - k) * ellipk(k2) - (2 / k) * ellipe(k2))

    coil = _ring(a1, dR, 0.0, hw).mat("coil").maxh(2 * MM)
    box = OrthoBrick(Pnt(-0.2, -0.2, -0.2), Pnt(0.2, 0.2, 0.2)).bc("outer")
    air = (box - coil).mat("air")
    geo = CSGeometry(); geo.Add(coil); geo.Add(air)
    mesh = Mesh(geo.GenerateMesh(maxh=0.04, grading=0.35)); mesh.Curve(2)

    gfu = _solve_Aform(mesh, NU0, I1 / ((2 * dR) * (2 * hw)))
    A_phi = gfu(mesh(a2, 0.0, d))[1]        # phi-component on the +x axis = A_y
    M_fem = 2 * math.pi * a2 * A_phi / I1
    assert abs(M_fem - M_ref) / abs(M_ref) < 0.08, (
        f"coaxial-loop M = {M_fem*1e9:.3f} nH, Maxwell exact = {M_ref*1e9:.3f} nH")


def validate_coaxial_loops_axial_force():
    """Axial force between two coaxial loops (a1=a2=50 mm, gap d=30 mm, I1=I2=100 A).
    NGSolve: solve loop-1, then the Lorentz force on loop-2 is
        F_z = -2*pi*a2*I2 * B_r1(a2,0,d)   (B_r1 = radial field of loop-1 there).
    Exact: F_z = I1 I2 dM/dd with M the Maxwell elliptic mutual inductance
    (parallel currents attract, F_z ~ -16.2 mN).

    Pass criterion: |F_fem - F_ref| / |F_ref| < 8 %.
    """
    from scipy.special import ellipk, ellipe

    def _M(a1, a2, dd):
        k2 = 4 * a1 * a2 / ((a1 + a2)**2 + dd**2); k = math.sqrt(k2)
        return MU0 * math.sqrt(a1 * a2) * ((2 / k - k) * ellipk(k2) - (2 / k) * ellipe(k2))

    a1 = 50 * MM; a2 = 50 * MM; d = 30 * MM; I1 = 100.0; I2 = 100.0
    dR = 1 * MM; hw = 1 * MM; h = 1e-4
    F_ref = I1 * I2 * (_M(a1, a2, d + h) - _M(a1, a2, d - h)) / (2 * h)

    coil = _ring(a1, dR, 0.0, hw).mat("coil").maxh(2 * MM)
    box = OrthoBrick(Pnt(-0.2, -0.2, -0.2), Pnt(0.2, 0.2, 0.2)).bc("outer")
    air = (box - coil).mat("air")
    geo = CSGeometry(); geo.Add(coil); geo.Add(air)
    mesh = Mesh(geo.GenerateMesh(maxh=0.04, grading=0.35)); mesh.Curve(2)

    B = curl(_solve_Aform(mesh, NU0, I1 / ((2 * dR) * (2 * hw))))
    B_r = B(mesh(a2, 0.0, d))[0]            # radial field = B_x on the +x axis
    F_fem = -2 * math.pi * a2 * I2 * B_r
    assert abs(F_fem - F_ref) / abs(F_ref) < 0.08, (
        f"coaxial-loop F_z = {F_fem*1e3:.3f} mN, exact I1 I2 dM/dz = {F_ref*1e3:.3f} mN")


def validate_anti_helmholtz_axial_gradient():
    """Anti-Helmholtz pair: two coaxial coils radius a=50 mm at z=+/-s (s=40 mm)
    with OPPOSING NI=1000 A-turns. On-axis B_z is antisymmetric; the central
    axial gradient is exact:
        G = dB_z/dz|_0 = 3 mu0 NI a^2 s / (a^2+s^2)^(5/2)  ~ 0.350 T/m.
    Measured by central difference of B_z at z=+/-5 mm.

    Pass criterion: |G_fem - G_ref| / G_ref < 8 %.
    """
    a = 50 * MM; s = 40 * MM; NI = 1000.0; dR = 1 * MM; hw = 1 * MM
    G_ref = 3 * MU0 * NI * a**2 * s / (a**2 + s**2)**2.5

    cp = _ring(a, dR,  s, hw).mat("coil_p").maxh(3.5 * MM)
    cm = _ring(a, dR, -s, hw).mat("coil_m").maxh(3.5 * MM)
    box = OrthoBrick(Pnt(-0.15, -0.15, -0.15), Pnt(0.15, 0.15, 0.15)).bc("outer")
    air = (box - cp - cm).mat("air")
    geo = CSGeometry(); geo.Add(cp); geo.Add(cm); geo.Add(air)
    mesh = Mesh(geo.GenerateMesh(maxh=0.06, grading=0.3)); mesh.Curve(2)

    fes = HCurl(mesh, order=2, dirichlet="outer")
    u, v = fes.TnT()
    r_xy = sqrt(x * x + y * y + 1e-12)
    phi = CoefficientFunction((-y, x, 0)) / r_xy
    J0 = NI / ((2 * dR) * (2 * hw))
    Jmag = CoefficientFunction([J0 if m == "coil_p" else (-J0 if m == "coil_m" else 0.0)
                                for m in mesh.GetMaterials()])
    a_bf = BilinearForm(fes)
    a_bf += NU0 * curl(u) * curl(v) * dx + 1e-6 * NU0 * u * v * dx
    f = LinearForm(fes); f += InnerProduct(Jmag * phi, v) * dx
    a_bf.Assemble(); f.Assemble()
    gfu = GridFunction(fes)
    gfu.vec.data = a_bf.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * f.vec

    B = curl(gfu); dl = 5 * MM
    G_fem = (B(mesh(0.001, 0.001, dl))[2] - B(mesh(0.001, 0.001, -dl))[2]) / (2 * dl)
    assert abs(G_fem - G_ref) / abs(G_ref) < 0.08, (
        f"anti-Helmholtz gradient = {G_fem:.4f} T/m, exact = {G_ref:.4f} T/m")


def validate_current_loop_offaxis_field():
    """OFF-axis field of a circular current loop (a=50 mm, I=100 A) at the
    cylindrical point (r=30 mm, z=20 mm) -- both components against the exact
    elliptic-integral field (Simpson / Jackson §5.5):
        B_z = (mu0 I/2pi)/sqrt(Q) [K + (a^2-r^2-z^2)/P * E]
        B_r = (mu0 I/2pi) z/(r sqrt(Q)) [-K + (a^2+r^2+z^2)/P * E]
    with Q=(a+r)^2+z^2, P=(a-r)^2+z^2, k^2=4ar/Q. Validates the full (r,z) field,
    not just the on-axis line.

    Pass criterion: |B_fem - B_exact| / |B_exact| < 5 % for both B_r and B_z.
    """
    from scipy.special import ellipk, ellipe
    a = 50 * MM; I = 100.0; r = 30 * MM; zp = 20 * MM
    dR = 1 * MM; hw = 1 * MM
    Q = (a + r)**2 + zp**2; P = (a - r)**2 + zp**2; k2 = 4 * a * r / Q
    c = MU0 * I / (2 * math.pi)
    Bz_ex = c / math.sqrt(Q) * (ellipk(k2) + (a**2 - r**2 - zp**2) / P * ellipe(k2))
    Br_ex = c * zp / (r * math.sqrt(Q)) * (-ellipk(k2) + (a**2 + r**2 + zp**2) / P * ellipe(k2))

    coil = _ring(a, dR, 0.0, hw).mat("coil").maxh(2 * MM)
    box = OrthoBrick(Pnt(-0.2, -0.2, -0.2), Pnt(0.2, 0.2, 0.2)).bc("outer")
    air = (box - coil).mat("air")
    geo = CSGeometry(); geo.Add(coil); geo.Add(air)
    mesh = Mesh(geo.GenerateMesh(maxh=0.04, grading=0.35)); mesh.Curve(2)

    B = curl(_solve_Aform(mesh, NU0, I / ((2 * dR) * (2 * hw))))
    Bvec = B(mesh(r, 0.0, zp))          # on +x axis: B_x = B_r, B_z = B_z
    Br_fem, Bz_fem = Bvec[0], Bvec[2]
    assert abs(Bz_fem - Bz_ex) / abs(Bz_ex) < 0.05, (
        f"loop off-axis B_z = {Bz_fem*1e3:.4f} mT, exact = {Bz_ex*1e3:.4f} mT")
    assert abs(Br_fem - Br_ex) / abs(Br_ex) < 0.05, (
        f"loop off-axis B_r = {Br_fem*1e3:.4f} mT, exact = {Br_ex*1e3:.4f} mT")


def validate_permeable_sphere_exterior_dipole():
    """EXTERIOR field of a permeable sphere (mu_r=10, R=50 mm) in uniform B0=1 T.
    Outside, the scattered field is an exact dipole; on the polar (z) axis
        B_z(z) = B0 + 2 R^3 (mu_r-1)/(mu_r+2) B0 / z^3   (z > R).
    Complements the INTERIOR check (validate_sphere_magnetostatic_scattered) using the
    same scattered-field A-form solve -- validates the exterior multipole.

    Pass criterion: |B_z_fem - B_z_ref| < 3 % of B_z_ref at z = 100, 150 mm.
    """
    mu_r = 10.0; B0 = 1.0; R = 50 * MM; L = 0.25
    coef = (mu_r - 1.0) / (mu_r + 2.0)
    sph = Sphere(Pnt(0, 0, 0), R)
    box = OrthoBrick(Pnt(-L, -L, -L), Pnt(L, L, L)).bc("outer")
    geo = CSGeometry()
    geo.Add(sph.mat("sphere").maxh(8 * MM)); geo.Add((box - sph).mat("air"))
    mesh = Mesh(geo.GenerateMesh(maxh=0.04, grading=0.4)); mesh.Curve(2)

    nu_cf   = mesh.MaterialCF({"sphere": NU0 / mu_r}, default=NU0)
    nu_pert = mesh.MaterialCF({"sphere": NU0 * (1.0 - 1.0 / mu_r)}, default=0.0)
    B0_vec  = CoefficientFunction((0.0, 0.0, B0))
    fes = HCurl(mesh, order=2, dirichlet="outer")
    u, v = fes.TnT()
    a = BilinearForm(fes)
    a += nu_cf * InnerProduct(curl(u), curl(v)) * dx + 1e-8 * NU0 * InnerProduct(u, v) * dx
    f = LinearForm(fes)
    f += InnerProduct(nu_pert * B0_vec, curl(v)) * dx
    a.Assemble(); f.Assemble()
    gfA = GridFunction(fes)
    gfA.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * f.vec
    B = curl(gfA)

    for zp in (100 * MM, 150 * MM):
        Bz_fem = float(B(mesh(1e-4, 1e-4, zp))[2]) + B0
        Bz_ref = B0 + 2 * R**3 * coef * B0 / zp**3
        assert abs(Bz_fem - Bz_ref) < 0.03 * Bz_ref, (
            f"permeable-sphere exterior B_z at z={zp*1e3:.0f}mm: "
            f"fem={Bz_fem:.4f} ref={Bz_ref:.4f} T")


def validate_coaxial_loops_mutual_parametric():
    """Mutual inductance of coaxial loops (a1=a2=50 mm) at THREE gaps d=20/40/60 mm
    from a SINGLE loop-1 solve -- M(d) = 2*pi*a2 A_phi(a2,0,d)/I1 vs the Maxwell
    elliptic M at each d. Validates the loop-1 field along the whole axis, not
    just one point.

    Pass criterion: |M_fem - M_ref| / M_ref < 8 % at every d.
    """
    from scipy.special import ellipk, ellipe

    def _M(a1, a2, dd):
        k2 = 4 * a1 * a2 / ((a1 + a2)**2 + dd**2); k = math.sqrt(k2)
        return MU0 * math.sqrt(a1 * a2) * ((2 / k - k) * ellipk(k2) - (2 / k) * ellipe(k2))

    a1 = 50 * MM; a2 = 50 * MM; I1 = 100.0; dR = 1 * MM; hw = 1 * MM
    coil = _ring(a1, dR, 0.0, hw).mat("coil").maxh(2 * MM)
    box = OrthoBrick(Pnt(-0.2, -0.2, -0.2), Pnt(0.2, 0.2, 0.2)).bc("outer")
    air = (box - coil).mat("air")
    geo = CSGeometry(); geo.Add(coil); geo.Add(air)
    mesh = Mesh(geo.GenerateMesh(maxh=0.04, grading=0.35)); mesh.Curve(2)
    gfu = _solve_Aform(mesh, NU0, I1 / ((2 * dR) * (2 * hw)))

    for d in (20 * MM, 40 * MM, 60 * MM):
        M_fem = 2 * math.pi * a2 * gfu(mesh(a2, 0.0, d))[1] / I1
        M_ref = _M(a1, a2, d)
        assert abs(M_fem - M_ref) / abs(M_ref) < 0.08, (
            f"M(d={d*1e3:.0f}mm) = {M_fem*1e9:.3f} nH, exact = {M_ref*1e9:.3f} nH")


def validate_magnetic_shielding_spherical_shell():
    """Magnetic shielding factor of a mu-metal SPHERICAL SHELL (a=40, b=60 mm) in a
    uniform applied field B0 -- ref AC/DC "magnetic shielding" benchmark (#5).

    Exact cavity field:
        B_in/B0 = 9 mu_r / [(2 mu_r+1)(mu_r+2) - 2 (a/b)^3 (mu_r-1)^2]
    -> shielding factor S = B0/B_in. For mu_r=100, a/b=2/3 : S ≈ 16.3.

    Scattered-field A-form (solve_scattered_uniform_field, permeability-jump source);
    cavity field read as a VOLUME AVERAGE. Cross-validated to ~1 % on 2026-06-05.
    """
    a, b, L, B0, mu_r = 0.04, 0.06, 0.30, 1.0, 100.0
    sph_a = Sphere(Pnt(0, 0, 0), a); sph_b = Sphere(Pnt(0, 0, 0), b)
    box = OrthoBrick(Pnt(-L, -L, -L), Pnt(L, L, L)).bc("outer")
    geo = CSGeometry()
    geo.Add(sph_a.mat("cavity").maxh(0.012))
    geo.Add((sph_b - sph_a).mat("shell").maxh(0.006))
    geo.Add((box - sph_b).mat("air"))
    mesh = Mesh(geo.GenerateMesh(maxh=0.05, grading=0.3)); mesh.Curve(2)

    gfA = solve_scattered_uniform_field(mesh, {"shell": mu_r}, (0.0, 0.0, B0))
    vol = Integrate(CoefficientFunction(1.0) * dx("cavity"), mesh)
    B_in = Integrate((curl(gfA)[2] + B0) * dx("cavity"), mesh) / vol
    S_fem = B0 / B_in
    S_ref = shell_shielding_factor(mu_r, a, b, "sphere")
    assert abs(S_fem - S_ref) / S_ref < 0.04, (
        f"shielding S_fem={S_fem:.2f}, exact={S_ref:.2f} (mu_r={mu_r})")

    # analytic helper sanity: no shielding at mu_r=1 (sphere AND cylinder)
    assert abs(shell_shielding_factor(1.0, a, b, "sphere") - 1.0) < 1e-12
    assert abs(shell_shielding_factor(1.0, a, b, "cylinder") - 1.0) < 1e-12


def validate_capacitance_3d_spherical():
    """3D capacitance of a SPHERICAL CAPACITOR (inner a, outer b, vacuum gap) --
    ref AC/DC "Computing Capacitance" benchmark (#6).

    Exact:  C = 4 pi eps0 a b / (b - a). Method: 3D electrostatics on H1 with the
    conductors as Dirichlet boundaries, capacitance from the field energy
    C = 2W/V0^2. Cross-validated to <0.2 % on 2026-06-05.
    """
    V0 = 1.0
    for a, b in ((0.05, 0.10), (0.04, 0.12)):
        sph_a = Sphere(Pnt(0, 0, 0), a).bc("inner")
        sph_b = Sphere(Pnt(0, 0, 0), b).bc("outer")
        geo = CSGeometry(); geo.Add((sph_b - sph_a).mat("gap").maxh(0.01))
        mesh = Mesh(geo.GenerateMesh(maxh=0.012, grading=0.3)); mesh.Curve(3)

        gfV = solve_electrostatic_3d(mesh, EPS0, {"inner": V0, "outer": 0.0}, order=2)
        C_fem = capacitance_from_energy(gfV, EPS0, V0, mesh)
        C_ref = spherical_capacitor_C(a, b)
        assert abs(C_fem - C_ref) / C_ref < 0.02, (
            f"C(a={a},b={b}) = {C_fem*1e12:.3f} pF, exact = {C_ref*1e12:.3f} pF")


def validate_capacitance_dielectric_layered():
    """Layered-dielectric spherical capacitor (eps_r1 in a<r<c, eps_r2 in c<r<b) --
    ref "Computing Capacitance" with dielectric materials (#7). Validates the
    piecewise-permittivity (mesh.MaterialCF) path of electrostatic3d. Exact:
    C = 4 pi eps0 / [(1/er1)(1/a-1/c) + (1/er2)(1/c-1/b)]. <0.2 % on 2026-06-05.
    """
    V0 = 1.0
    for a, c, b, er1, er2 in ((0.04, 0.06, 0.10, 2.0, 5.0), (0.03, 0.05, 0.09, 4.0, 2.0)):
        sph_a = Sphere(Pnt(0, 0, 0), a).bc("inner"); sph_c = Sphere(Pnt(0, 0, 0), c)
        sph_b = Sphere(Pnt(0, 0, 0), b).bc("outer")
        geo = CSGeometry()
        geo.Add((sph_c - sph_a).mat("diel1").maxh(0.008))
        geo.Add((sph_b - sph_c).mat("diel2").maxh(0.010))
        mesh = Mesh(geo.GenerateMesh(maxh=0.012, grading=0.3)); mesh.Curve(3)
        eps_cf = mesh.MaterialCF({"diel1": EPS0 * er1, "diel2": EPS0 * er2})
        gfV = solve_electrostatic_3d(mesh, eps_cf, {"inner": V0, "outer": 0.0}, order=2)
        C_fem = capacitance_from_energy(gfV, eps_cf, V0, mesh)
        C_ref = 4 * math.pi * EPS0 / ((1/er1)*(1/a - 1/c) + (1/er2)*(1/c - 1/b))
        assert abs(C_fem - C_ref) / C_ref < 0.02, (
            f"C(layered, er=({er1},{er2})) = {C_fem*1e12:.3f} pF, exact = {C_ref*1e12:.3f} pF")


def validate_halbach_cylinder_bore_field():
    """Ideal dipole HALBACH cylinder (Br=1.2 T, Ri=20, Ro=40 mm) -- accelerator/
    undulator magnet (#8). Rotating remanence (easy axis = 2*phi) gives a UNIFORM
    transverse bore field B = Br ln(Ro/Ri). 2D planar A_z with
    halbach_dipole_magnetization; cross-validated to ~0.3 % + uniformity 2026-06-05.
    """
    Br, Ri, Ro, Rbox = 1.2, 0.02, 0.04, 0.12
    geo = SplineGeometry()
    geo.AddCircle((0, 0), Rbox, leftdomain=3, rightdomain=0, bc="outer")
    geo.AddCircle((0, 0), Ro, leftdomain=2, rightdomain=3)
    geo.AddCircle((0, 0), Ri, leftdomain=1, rightdomain=2)
    geo.SetMaterial(1, "bore"); geo.SetMaterial(2, "magnet"); geo.SetMaterial(3, "air")
    mesh = Mesh(geo.GenerateMesh(maxh=0.004))

    M_cf = halbach_dipole_magnetization(mesh, Br, region="magnet")
    gfu = solve_planar_magnetostatic(mesh, NU0, magnet_cf=M_cf, order=3, dirichlet="outer")
    B = CoefficientFunction((grad(gfu)[1], -grad(gfu)[0]))
    Bmag = sqrt(B[0] * B[0] + B[1] * B[1])
    B_ref = halbach_bore_field(Br, Ri, Ro)

    vals = [float(Bmag(mesh(px, py)))
            for px, py in ((0, 0), (0.01, 0), (0, 0.01), (0.007, 0.007))]
    for v in vals:
        assert abs(v - B_ref) / B_ref < 0.03, f"|B|={v:.4f} T vs exact {B_ref:.4f} T"
    assert (max(vals) - min(vals)) / B_ref < 0.01, f"bore field not uniform: {vals}"


def validate_capacitance_matrix_spherical():
    """Maxwell capacitance MATRIX of a closed 2-conductor spherical capacitor (#9) --
    ref "Computing Capacitance" matrix. Exact C = [[C0,-C0],[-C0,C0]],
    C0 = 4 pi eps0 ab/(b-a). Validates the FEM reaction charge extraction, symmetry,
    and zero row-sum (charge conservation). <0.5 % on 2026-06-05.
    """
    import numpy as np
    a, b = 0.05, 0.10
    sph_a = Sphere(Pnt(0, 0, 0), a).bc("inner"); sph_b = Sphere(Pnt(0, 0, 0), b).bc("outer")
    geo = CSGeometry(); geo.Add((sph_b - sph_a).mat("gap").maxh(0.01))
    mesh = Mesh(geo.GenerateMesh(maxh=0.012, grading=0.3)); mesh.Curve(3)
    C = capacitance_matrix(mesh, EPS0, ["inner", "outer"], order=2)
    C0 = spherical_capacitor_C(a, b)
    C_ref = np.array([[C0, -C0], [-C0, C0]])
    assert np.max(np.abs(C - C_ref)) / C0 < 0.02, f"C matrix off: {C*1e12} pF"
    assert np.max(np.abs(C - C.T)) / C0 < 0.01, "capacitance matrix not symmetric"
    assert np.max(np.abs(C.sum(axis=1))) / C0 < 0.01, "row-sum != 0 (charge not conserved)"


def validate_electrostatic_force_parallel_plate():
    """Electrostatic force on a parallel-plate capacitor (#10, MEMS actuator). A closed
    gap box with NEUMANN side walls gives an exact 1-D field, so |F| = eps0 A V^2/(2d^2)
    exactly (no fringing). Maxwell-stress electrostatic eggshell. <2 % on 2026-06-05.
    """
    from netgen.occ import Box, OCCGeometry, Z
    V0, d, L = 100.0, 0.001, 0.005
    box = Box((-L/2, -L/2, 0), (L/2, L/2, d)); box.mat("air")
    box.faces.Min(Z).name = "bottom"; box.faces.Max(Z).name = "top"
    mesh = Mesh(OCCGeometry(box).GenerateMesh(maxh=0.00025))
    gfV = solve_electrostatic_3d(mesh, EPS0, {"top": V0, "bottom": 0.0}, order=2)
    E = -grad(gfV)
    z0, z1 = 0.3 * d, 0.7 * d
    gz = IfPos(z - z0, IfPos(z1 - z, 1.0 / (z1 - z0), 0.0), 0.0)
    gradg = CoefficientFunction((0, 0, gz))
    Fx, Fy, Fz = electrostatic_eggshell_force(E, mesh, gradg, air_region="air")
    F_ref = EPS0 * (L * L) * V0**2 / (2.0 * d**2)
    assert abs(abs(Fz) - F_ref) / F_ref < 0.03, f"F={Fz:.3e} N vs exact {F_ref:.3e} N"
    assert Fz < 0, "parallel-plate force must be attractive (downward)"


def validate_magnetic_shielding_cylinder():
    """CYLINDER magnetic shielding (2D transverse, #11) -- the infinite-length twin of
    the sphere (#5). mu-metal shell (a=40, b=60 mm) in a uniform transverse B0; exact
    S = [(mu_r+1)^2 - (a/b)^2(mu_r-1)^2]/(4 mu_r). 2D planar A_z, permeability-jump
    magnetization source. ~1.3 % (box at 10*b for the 1/r 2D dipole) on 2026-06-05.
    """
    B0, a, b, Rbox, mu_r = 1.0, 0.04, 0.06, 0.60, 100.0
    geo = SplineGeometry()
    geo.AddCircle((0, 0), Rbox, leftdomain=3, rightdomain=0, bc="outer")
    geo.AddCircle((0, 0), b, leftdomain=2, rightdomain=3)
    geo.AddCircle((0, 0), a, leftdomain=1, rightdomain=2)
    geo.SetMaterial(1, "bore"); geo.SetMaterial(2, "shell"); geo.SetMaterial(3, "air")
    mesh = Mesh(geo.GenerateMesh(maxh=0.01))
    nu_cf = mesh.MaterialCF({"shell": NU0 / mu_r}, default=NU0)
    ind = mesh.MaterialCF({"shell": 1.0}, default=0.0)
    magnet_cf = ind * (NU0 * (1.0 - 1.0 / mu_r)) * CoefficientFunction((B0, 0.0))
    gfu = solve_planar_magnetostatic(mesh, nu_cf, magnet_cf=magnet_cf, order=3)
    vol = Integrate(CoefficientFunction(1.0) * dx("bore"), mesh)
    Bx_in = Integrate((grad(gfu)[1] + B0) * dx("bore"), mesh) / vol
    S_fem = B0 / Bx_in
    S_ref = shell_shielding_factor(mu_r, a, b, "cylinder")
    assert abs(S_fem - S_ref) / S_ref < 0.025, f"S_fem={S_fem:.2f} vs exact {S_ref:.2f}"


def validate_two_wire_line_capacitance():
    """Two-wire transmission-line capacitance per unit length (#12). Exact
    C = pi eps0 / arccosh(D/2a). 2D electrostatics, differential drive (+-V0/2)
    with a far ground ring. <1 % on 2026-06-05.
    """
    a, D, Rbox, V0 = 0.002, 0.01, 0.15, 1.0
    geo = SplineGeometry()
    geo.AddCircle((0, 0), Rbox, leftdomain=1, rightdomain=0, bc="ground")
    geo.AddCircle((-D / 2, 0), a, leftdomain=0, rightdomain=1, bc="wire1", maxh=0.0004)
    geo.AddCircle((+D / 2, 0), a, leftdomain=0, rightdomain=1, bc="wire2", maxh=0.0004)
    geo.SetMaterial(1, "air")
    mesh = Mesh(geo.GenerateMesh(maxh=0.01))
    V = solve_electrostatic_2d(mesh, EPS0, {"wire1": V0 / 2, "wire2": -V0 / 2, "ground": 0.0}, order=2)
    C_fem = capacitance_2d(V, mesh, EPS0, V0)
    C_ref = math.pi * EPS0 / math.acosh(D / (2 * a))
    assert abs(C_fem - C_ref) / C_ref < 0.02, f"C={C_fem*1e12:.3f} pF/m vs exact {C_ref*1e12:.3f}"


def validate_multipole_quadrupole_fem():
    """ref-class #13: an air-cored NORMAL quadrupole (4 line currents
    +I,-I,+I,-I at 0/90/180/270 deg, radius r0) solved as planar A_z; the
    multipoles extracted from the FEM field on R_ref must match the exact
    free-space superposition -- main b2 to <1 %, the allowed 12-pole (n=6) at
    (R_ref/r0)^4, and the forbidden harmonics (n=1,3,4,5) below 1 unit (1e-4)."""
    from netgen.geom2d import SplineGeometry
    NU0 = 1.0 / MU0
    I, r0, rw, R_ref, Rout = 1000.0, 0.05, 0.004, 0.02, 0.40
    signs = [+1, -1, +1, -1]
    angles = [math.radians(90 * k) for k in range(4)]
    geo = SplineGeometry()
    geo.AddCircle((0, 0), Rout, leftdomain=1, rightdomain=0, bc="outer")
    for k, ang in enumerate(angles):
        geo.AddCircle((r0 * math.cos(ang), r0 * math.sin(ang)), rw,
                      leftdomain=2 + k, rightdomain=1)
        geo.SetMaterial(2 + k, f"w{k}")
    geo.SetMaterial(1, "air")
    mesh = Mesh(geo.GenerateMesh(maxh=0.015))
    # area-normalize so int Jz = I exactly on the discrete wire region
    Jz_terms = {}
    for k in range(4):
        area = Integrate(mesh.MaterialCF({f"w{k}": 1.0}, default=0.0) * dx, mesh)
        Jz_terms[f"w{k}"] = signs[k] * I / area
    Jz = mesh.MaterialCF(Jz_terms, default=0.0)
    gfu = solve_planar_magnetostatic(mesh, CoefficientFunction(NU0), Jz=Jz, order=4)
    B = CoefficientFunction((grad(gfu)[1], -grad(gfu)[0]))
    fem = multipole_coefficients(B, R_ref, n_max=10, mesh=mesh)
    refs = [line_current_multipoles(signs[k] * I, r0, angles[k], R_ref, 10)
            for k in range(4)]
    ana = superpose_multipoles(*refs)
    fe = field_errors(fem)
    assert fe["main_n"] == 2                                   # quadrupole
    assert abs(fem["b"][2] - ana["b"][2]) / abs(ana["b"][2]) < 0.01   # main to 1 %
    # allowed 12-pole present at (R_ref/r0)^4 of the main
    assert fem["C"][6] / fem["C"][2] == pytest.approx((R_ref / r0) ** 4, rel=0.03)
    # forbidden harmonics are numerical noise (< 1 unit of 1e-4 of the main)
    for n in (1, 3, 4, 5):
        assert fem["C"][n] / fe["main_C"] * 1e4 < 1.0


def validate_cylinder_magnet_axial_field():
    """ref-class #14: on-axis B_z of an axially-magnetized CYLINDER permanent
    magnet vs the closed form B_z(z)=(Br/2)[(z+L/2)/sqrt(R^2+(z+L/2)^2) -
    (z-L/2)/sqrt(R^2+(z-L/2)^2)]. Axisymmetric A_phi with a rigid magnetization
    M=Br/mu0. The CENTRE field (the headline magnet spec) and the FAR field
    validate tightly; the immediate pole exit (z~L) is corner-singularity
    sensitive and is intentionally not asserted (see cylinder_magnet.py)."""
    from netgen.occ import OCCGeometry, MoveTo, Glue, X, Y
    from ngsolve import L2, GridFunction, x as r_cf
    from radia_mcp.radia_ngsolve.solve import (solve_axi_magnetostatic,
                                               cylinder_magnet_axial_field)
    NU0 = 1.0 / MU0
    Br, R, L = 1.2, 0.010, 0.020
    Rbox, nearR, nearZ = 0.40, 8 * R, 6 * L
    magnet = MoveTo(0, -L / 2).Rectangle(R, L).Face()
    magnet.faces.name = "magnet"; magnet.maxh = 0.001
    magnet.edges.Min(X).name = "axis"; magnet.edges.Min(X).maxh = 0.0003
    magnet.edges.Max(Y).maxh = 0.0004; magnet.edges.Min(Y).maxh = 0.0004
    near = MoveTo(0, -nearZ).Rectangle(nearR, 2 * nearZ).Face() - magnet
    near.faces.name = "air"; near.maxh = 0.001        # external on-axis field needs this
    near.edges.Min(X).name = "axis"; near.edges.Min(X).maxh = 0.0003
    far = MoveTo(0, -Rbox).Rectangle(Rbox, 2 * Rbox).Face() \
        - MoveTo(0, -nearZ).Rectangle(nearR, 2 * nearZ).Face()
    far.faces.name = "air"; far.maxh = 0.05
    far.edges.Min(X).name = "axis"
    far.edges.Max(X).name = "outer"; far.edges.Max(Y).name = "outer"; far.edges.Min(Y).name = "outer"
    mesh = Mesh(OCCGeometry(Glue([far, near, magnet]), dim=2).GenerateMesh(maxh=0.05))
    gfu = solve_axi_magnetostatic(mesh, CoefficientFunction(NU0),
                                  magnets={"magnet": (Br / MU0, 90.0)}, order=2)
    bzgf = GridFunction(L2(mesh, order=2)); bzgf.Set(grad(gfu)[0] + gfu / r_cf)

    def bz0(z):
        return sum(bzgf(mesh(r, z)) for r in (R / 25, R / 22, R / 20, R / 18, R / 16)) / 5

    b0_ex = cylinder_magnet_axial_field(Br, R, L, 0.0)
    assert b0_ex == pytest.approx(Br * L / (2 * math.sqrt(R * R + (L / 2) ** 2)), rel=1e-12)
    # regression guard. A broken magnetization convention would be tens of %.
    assert abs(bz0(0.0) - b0_ex) / b0_ex < 0.025                      # centre field, headline spec
    bf_ex = cylinder_magnet_axial_field(Br, R, L, 2 * L)
    assert abs(bz0(2 * L) - bf_ex) / bf_ex < 0.06                      # far field


def validate_solenoid_field_and_inductance():
    """ref-class #15: finite air-core SOLENOID. On-axis B_z vs the exact closed
    form B_z=(mu0 n I/2)[(L/2-z)/sqrt(..)+(L/2+z)/sqrt(..)] (centre + bore), and the
    self-inductance L=2W/I^2 (energy method) vs the long limit times the Nagaoka
    coefficient (Wheeler 1/(1+0.9 a/L)). Smooth on axis -> clean L2 extraction."""
    from netgen.occ import OCCGeometry, MoveTo, Glue, X, Y
    from ngsolve import L2, GridFunction, x as r_cf
    from radia_mcp.radia_ngsolve.solve import (solve_axi_magnetostatic,
        solenoid_axial_field, solenoid_inductance_long)
    from radia_mcp.radia_ngsolve.force import inductance_axi
    NU0 = 1.0 / MU0
    a, L, t, n, I = 0.05, 0.30, 0.001, 1000.0, 2.0
    Rbox, Zbox, nearR, nearZ = 1.20, 1.20, 0.12, 0.50
    coil = MoveTo(a - t / 2, -L / 2).Rectangle(t, L).Face()
    coil.faces.name = "coil"; coil.maxh = 0.0008
    near = MoveTo(0, -nearZ).Rectangle(nearR, 2 * nearZ).Face() - coil
    near.faces.name = "air"; near.maxh = 0.004
    near.edges.Min(X).name = "axis"; near.edges.Min(X).maxh = 0.001
    far = MoveTo(0, -Zbox).Rectangle(Rbox, 2 * Zbox).Face() \
        - MoveTo(0, -nearZ).Rectangle(nearR, 2 * nearZ).Face()
    far.faces.name = "air"; far.maxh = 0.06
    far.edges.Min(X).name = "axis"
    far.edges.Max(X).name = "outer"; far.edges.Max(Y).name = "outer"; far.edges.Min(Y).name = "outer"
    mesh = Mesh(OCCGeometry(Glue([far, near, coil]), dim=2).GenerateMesh(maxh=0.06))
    Jr = mesh.MaterialCF({"coil": n * I / t}, default=0.0)
    gfu = solve_axi_magnetostatic(mesh, CoefficientFunction(NU0), Jr=Jr, order=2)
    B = CoefficientFunction((grad(gfu)[0] + gfu / r_cf, -grad(gfu)[1]))
    bzgf = GridFunction(L2(mesh, order=2)); bzgf.Set(grad(gfu)[0] + gfu / r_cf)

    def bz0(z):
        return sum(bzgf(mesh(r, z)) for r in (0.001, 0.0015, 0.002, 0.0025, 0.003)) / 5

    for z in (0.0, L / 4, 0.45 * L):                          # centre + bore
        ex = solenoid_axial_field(n, I, a, L, z)
        assert abs(bz0(z) - ex) / abs(ex) < 0.025, f"z={z}: {bz0(z):.6f} vs {ex:.6f}"
    # self-inductance vs the Nagaoka (Wheeler) reference
    kN = inductance_axi(B, mesh, NU0, I) / solenoid_inductance_long(n, a, L)
    assert abs(kN - 1.0 / (1.0 + 0.9 * a / L)) / (1.0 / (1.0 + 0.9 * a / L)) < 0.02, f"k_N={kN:.4f}"


def validate_busbar_lorentz_force():
    """ref-class #16: Lorentz force per length between two parallel busbars vs
    F = mu0 I1 I2/(2 pi d), both attractive (parallel) and repulsive (anti-parallel),
    via the J x B volume integral (force.lorentz_force_2d) on the planar A_z field."""
    from netgen.geom2d import SplineGeometry
    from radia_mcp.radia_ngsolve.force import lorentz_force_2d
    NU0 = 1.0 / MU0
    I1, I2, d, rw, Rbox = 1000.0, 1000.0, 0.04, 0.002, 0.40
    F_exact = MU0 * I1 * I2 / (2.0 * math.pi * d)
    geo = SplineGeometry()
    geo.AddCircle((0, 0), Rbox, leftdomain=1, rightdomain=0, bc="outer")
    geo.AddCircle((-d / 2, 0), rw, leftdomain=2, rightdomain=1); geo.SetMaterial(2, "w1")
    geo.AddCircle((+d / 2, 0), rw, leftdomain=3, rightdomain=1); geo.SetMaterial(3, "w2")
    geo.SetMaterial(1, "air")
    mesh = Mesh(geo.GenerateMesh(maxh=0.02))
    A1 = Integrate(mesh.MaterialCF({"w1": 1.0}, default=0.0) * dx, mesh)
    A2 = Integrate(mesh.MaterialCF({"w2": 1.0}, default=0.0) * dx, mesh)
    for s2, want_sign in ((+1.0, -1.0), (-1.0, +1.0)):       # parallel attract / anti repel
        Jz = mesh.MaterialCF({"w1": I1 / A1, "w2": s2 * I2 / A2}, default=0.0)
        gfu = solve_planar_magnetostatic(mesh, CoefficientFunction(NU0), Jz=Jz, order=3)
        B = CoefficientFunction((grad(gfu)[1], -grad(gfu)[0]))
        Fx, Fy = lorentz_force_2d(Jz, B, mesh, "w2")
        assert abs(abs(Fx) - F_exact) / F_exact < 0.03, f"|Fx|={abs(Fx):.4f} vs {F_exact:.4f}"
        assert Fx * want_sign > 0, f"wrong sign: Fx={Fx:.4f} (s2={s2})"
        assert abs(Fy) < 0.05 * F_exact                       # force is along x


def validate_helmholtz_uniformity():
    """ref-class #17: HELMHOLTZ pair (loops radius a at z=+/-a/2). On-axis B_z vs
    the closed form, centre B0 vs (4/5)^(3/2) mu0 NI/a, and the field UNIFORMITY over
    the central |z|<=a/4 working region (the Helmholtz figure of merit)."""
    from netgen.occ import OCCGeometry, MoveTo, Glue, X, Y
    from ngsolve import L2, GridFunction, x as r_cf
    from radia_mcp.radia_ngsolve.solve import solve_axi_magnetostatic, coil_pair_axial_field
    NU0 = 1.0 / MU0
    a, N, I, w = 0.10, 100.0, 10.0, 0.004
    s, Rbox, Zbox, nearR, nearZ = a, 0.80, 0.80, 0.20, 0.40

    def loop(zc):
        f = MoveTo(a - w / 2, zc - w / 2).Rectangle(w, w).Face(); f.maxh = 0.0015
        return f
    top, bot = loop(+s / 2), loop(-s / 2)
    top.faces.name = "ltop"; bot.faces.name = "lbot"
    near = MoveTo(0, -nearZ).Rectangle(nearR, 2 * nearZ).Face() - top - bot
    near.faces.name = "air"; near.maxh = 0.004
    near.edges.Min(X).name = "axis"; near.edges.Min(X).maxh = 0.0006
    far = MoveTo(0, -Zbox).Rectangle(Rbox, 2 * Zbox).Face() \
        - MoveTo(0, -nearZ).Rectangle(nearR, 2 * nearZ).Face()
    far.faces.name = "air"; far.maxh = 0.06
    far.edges.Min(X).name = "axis"
    far.edges.Max(X).name = "outer"; far.edges.Max(Y).name = "outer"; far.edges.Min(Y).name = "outer"
    mesh = Mesh(OCCGeometry(Glue([far, near, top, bot]), dim=2).GenerateMesh(maxh=0.06))
    Jr = mesh.MaterialCF({"ltop": N * I / (w * w), "lbot": N * I / (w * w)}, default=0.0)
    gfu = solve_axi_magnetostatic(mesh, CoefficientFunction(NU0), Jr=Jr, order=2)
    bzgf = GridFunction(L2(mesh, order=2)); bzgf.Set(grad(gfu)[0] + gfu / r_cf)

    def bz0(z):
        return sum(bzgf(mesh(r, zz)) for r in (0.001, 0.0015, 0.002, 0.0025, 0.003)
                   for zz in (z - 2e-4, z + 2e-4)) / 10

    B0_exact = (4.0 / 5.0) ** 1.5 * MU0 * N * I / a
    assert B0_exact == pytest.approx(coil_pair_axial_field(N, I, a, s, 0.0), rel=1e-12)
    assert abs(bz0(0.0) - B0_exact) / B0_exact < 0.015                # centre field
    assert abs(bz0(a / 4) - coil_pair_axial_field(N, I, a, s, a / 4)) \
        / coil_pair_axial_field(N, I, a, s, a / 4) < 0.015            # off-centre point
    B0 = bz0(0.0)
    unif = max(abs(bz0(k * a / 40) / B0 - 1.0) for k in range(-10, 11))
    assert unif < 0.015, f"uniformity {unif*100:.2f}% over |z|<=a/4"   # Helmholtz flatness


def validate_c_magnet_gap_field():
    """ref-class #18: iron window-frame DIPOLE air-gap field vs the reluctance
    model B_gap = mu0 NI/(g + l_fe/mu_r). 2D planar A_z, high-mu yoke, coil threading
    the left leg; gap field = avg |B_x| over the gap. FEM sits just below the lumped
    model (fringing) and below the mu_r->inf ideal mu0 NI/g."""
    from netgen.occ import OCCGeometry, MoveTo, Glue, X, Y
    from radia_mcp.radia_ngsolve.solve import (solve_planar_magnetostatic,
                                               magnetic_circuit_gap_field)
    NU0 = 1.0 / MU0
    W, ww, b, g, mu_r, NI, Rbox = 0.24, 0.12, 0.06, 0.006, 2000.0, 5000.0, 0.60
    l_fe = 2.0 * ((W + ww) / 2.0) * 2 - g

    def rect(x0, y0, w, h):
        return MoveTo(x0, y0).Rectangle(w, h).Face()
    big, window = rect(-W/2, -W/2, W, W), rect(-ww/2, -ww/2, ww, ww)
    gap = rect(-g/2, ww/2, g, b)
    cin = rect(-ww/2 + 0.005, -0.04, 0.02, 0.08)
    cout = rect(-W/2 - 0.025, -0.04, 0.02, 0.08)
    iron = big - window - gap
    gap.faces.name = "gap"; cin.faces.name = "cin"; cout.faces.name = "cout"; iron.faces.name = "iron"
    air = rect(-Rbox, -Rbox, 2*Rbox, 2*Rbox) - iron - gap - cin - cout
    air.faces.name = "air"
    air.edges.Max(X).name = "outer"; air.edges.Min(X).name = "outer"
    air.edges.Max(Y).name = "outer"; air.edges.Min(Y).name = "outer"
    mesh = Mesh(OCCGeometry(Glue([air, iron, gap, cin, cout]), dim=2).GenerateMesh(maxh=0.04))
    Ain = Integrate(mesh.MaterialCF({"cin": 1.0}, default=0.0) * dx, mesh)
    Aout = Integrate(mesh.MaterialCF({"cout": 1.0}, default=0.0) * dx, mesh)
    Jz = mesh.MaterialCF({"cin": NI / Ain, "cout": -NI / Aout}, default=0.0)
    nu = mesh.MaterialCF({"iron": NU0 / mu_r}, default=NU0)
    gfu = solve_planar_magnetostatic(mesh, nu, Jz=Jz, order=3, dirichlet="outer")
    B = CoefficientFunction((grad(gfu)[1], -grad(gfu)[0]))
    Agap = Integrate(mesh.MaterialCF({"gap": 1.0}, default=0.0) * dx, mesh)
    Bx = abs(Integrate(mesh.MaterialCF({"gap": 1.0}, default=0.0) * B[0] * dx, mesh) / Agap)
    B_ref = magnetic_circuit_gap_field(1.0, NI, g, l_fe, mu_r)
    assert abs(Bx - B_ref) / B_ref < 0.03, f"B_gap {Bx:.4f} vs reluctance {B_ref:.4f}"
    assert Bx < MU0 * NI / g                                  # below the mu_r->inf ideal (fringing)
    assert Bx > 0.85 * MU0 * NI / g                           # but close to it (small gap)


def validate_joule_heating_electrothermal():
    """ref-class #19: ELECTRO-THERMAL Joule heating, chaining the two NGSolve
    solvers solve_current_flow -> joule_heat_source -> solve_heat_steady. A uniform
    bar (voltage V across length L, ends cold, sides insulated) has q=sigma(V/L)^2 and
    the exact parabolic rise dT(x)=(q/2k)x(L-x), peak sigma V^2/(8k)."""
    from netgen.occ import OCCGeometry, MoveTo, X, Y
    from radia_mcp.radia_ngsolve.scalar_fem2d import solve_current_flow
    from radia_mcp.radia_ngsolve.multiphysics import (solve_heat_steady, joule_heat_source,
                                                      joule_bar_temperature_rise)
    L, hh, sigma, k, V0 = 0.10, 0.02, 1.0e6, 50.0, 0.10
    bar = MoveTo(0, 0).Rectangle(L, hh).Face()
    bar.faces.name = "bar"
    bar.edges.Min(X).name = "left"; bar.edges.Max(X).name = "right"
    bar.edges.Min(Y).name = "bottom"; bar.edges.Max(Y).name = "top"
    mesh = Mesh(OCCGeometry(bar, dim=2).GenerateMesh(maxh=0.003))
    V = solve_current_flow(mesh, sigma, {"left": V0, "right": 0.0}, order=2)
    q = joule_heat_source(V, sigma)
    assert abs(Integrate(q * dx, mesh) / (L * hh) - sigma * (V0 / L) ** 2) / (sigma * (V0 / L) ** 2) < 0.02
    dT = solve_heat_steady(mesh, q, k, conductor="bar", dirichlet="left|right", order=2)
    for xl in (0.25, 0.5, 0.75):
        de = joule_bar_temperature_rise(sigma, V0, L, k, xl * L)
        assert abs(dT(mesh(xl * L, hh / 2)) - de) / de < 0.01
    assert abs(dT(mesh(L / 2, hh / 2)) - sigma * V0 ** 2 / (8 * k)) / (sigma * V0 ** 2 / (8 * k)) < 0.01


def validate_electro_thermo_mechanical_chain():
    """ref-class #20: 3-physics chain solve_current_flow -> solve_heat_steady ->
    solve_linear_elasticity. Validates dT_max = sigma V^2/(8k) and the constrained-bar
    axial thermal stress sigma_xx = -E alpha <dT>. Also locks the new elasticity solver
    on a pure uniaxial-tension case (delta = sigma0 L/E exactly)."""
    from netgen.occ import OCCGeometry, MoveTo, X, Y
    from radia_mcp.radia_ngsolve.scalar_fem2d import solve_current_flow
    from radia_mcp.radia_ngsolve.multiphysics import solve_heat_steady, joule_heat_source
    from radia_mcp.radia_ngsolve.elasticity import (solve_linear_elasticity, stress_2d,
                                                    constrained_bar_thermal_stress)
    L, hh, sigma, k, V0 = 0.10, 0.02, 1.0e6, 50.0, 0.10
    E, nu, alpha = 200.0e9, 0.30, 1.2e-5

    def bar():
        f = MoveTo(0, 0).Rectangle(L, hh).Face(); f.faces.name = "bar"
        f.edges.Min(X).name = "left"; f.edges.Max(X).name = "right"
        f.edges.Min(Y).name = "bottom"; f.edges.Max(Y).name = "top"
        return Mesh(OCCGeometry(f, dim=2).GenerateMesh(maxh=0.003))

    # pure tension: locks the elasticity solver (exact linear solution)
    m = bar(); s0 = 1.0e6
    u0 = solve_linear_elasticity(m, E, nu, dirichletx="left", dirichlety="bottom",
                                 traction={"right": (s0, 0.0)}, plane="stress")
    assert abs(u0(m(L, hh/2))[0] - s0 * L / E) / (s0 * L / E) < 1e-3

    # the 3-physics chain
    m = bar()
    V = solve_current_flow(m, sigma, {"left": V0, "right": 0.0}, order=2)
    dT = solve_heat_steady(m, joule_heat_source(V, sigma), k, conductor="bar",
                           dirichlet="left|right", order=2)
    assert abs(dT(m(L/2, hh/2)) - sigma * V0**2 / (8*k)) / (sigma * V0**2 / (8*k)) < 0.01
    dT_mean = Integrate(dT * dx, m) / (L * hh)
    u = solve_linear_elasticity(m, E, nu, dirichletx="left|right", dirichlety="bottom",
                                thermal=(alpha, dT), plane="stress")
    sxx = Integrate(stress_2d(u, E, nu, "stress", thermal=(alpha, dT))[0] * dx, m) / (L * hh)
    sxx_ex = constrained_bar_thermal_stress(E, alpha, dT_mean)
    assert abs(sxx - sxx_ex) / abs(sxx_ex) < 0.01


def validate_magneto_mechanical_beam():
    """ref-class #21: magneto-mechanical -- a current-carrying cantilever in a
    transverse field B0 feels a Lorentz body force f_y=-J_x B0 and deflects; the tip
    deflection matches Euler-Bernoulli w L^4/(8EI), w=I B0, for a slender beam."""
    from netgen.occ import OCCGeometry, MoveTo, X, Y
    from radia_mcp.radia_ngsolve.elasticity import (solve_linear_elasticity,
                                                    cantilever_tip_deflection)
    L, h, E, nu, B0, I = 0.20, 0.008, 70.0e9, 0.33, 0.5, 100.0
    Jx = I / h
    beam = MoveTo(0, 0).Rectangle(L, h).Face(); beam.faces.name = "beam"
    beam.edges.Min(X).name = "left"; beam.edges.Max(X).name = "right"
    beam.edges.Min(Y).name = "bottom"; beam.edges.Max(Y).name = "top"
    mesh = Mesh(OCCGeometry(beam, dim=2).GenerateMesh(maxh=0.002))
    u = solve_linear_elasticity(mesh, E, nu, dirichlet="left",
                                body_force=(0.0, -Jx * B0), plane="stress", order=2)
    tip_fem = u(mesh(L, h / 2))[1]
    tip_eb = -cantilever_tip_deflection(abs(Jx * B0) * h, L, E, h ** 3 / 12.0)
    assert abs(tip_fem - tip_eb) / abs(tip_eb) < 0.03
    assert tip_fem < 0                                     # deflects in the force direction


def validate_mems_electro_mechanical():
    """ref-class #22: MEMS electro-mechanical chain -- solve_electrostatic ->
    electrostatic_eggshell_force_2d (P = 1/2 eps0 (V0/d)^2, EXACT with a whole-gap
    CONSTANT weight ramp gradg=(0,1/d)) -> solve_linear_elasticity (the pull deflects a
    cantilever, tip vs Euler-Bernoulli). The electric twin of #21. A boundary trace of
    grad(V) would read ~0 (V=const on the conductor); the volume band is the fix."""
    from netgen.occ import OCCGeometry, MoveTo, X, Y
    from radia_mcp.radia_ngsolve.scalar_fem2d import solve_electrostatic
    from radia_mcp.radia_ngsolve.force import electrostatic_eggshell_force_2d, EPS0
    from radia_mcp.radia_ngsolve.elasticity import (solve_linear_elasticity,
                                                    cantilever_tip_deflection)
    W, d, V0 = 20e-3, 1e-3, 500.0
    face = MoveTo(0, 0).Rectangle(W, d).Face(); face.faces.name = "air"
    face.edges.Max(Y).name = "top"; face.edges.Min(Y).name = "bottom"
    mesh = Mesh(OCCGeometry(face, dim=2).GenerateMesh(maxh=d / 6))
    V = solve_electrostatic(mesh, EPS0, {"top": V0, "bottom": 0.0}, order=2)
    E = CoefficientFunction((-grad(V)[0], -grad(V)[1]))
    Fx, Fy = electrostatic_eggshell_force_2d(E, mesh, CoefficientFunction((0.0, 1.0 / d)),
                                             air_region="air")
    P_fem = abs(Fy) / W
    P_exact = 0.5 * EPS0 * (V0 / d) ** 2
    assert abs(P_fem - P_exact) / P_exact < 1e-3           # constant ramp -> exact

    Lb, h, Emod, nu = 10e-3, 0.4e-3, 1.0e8, 0.3
    beam = MoveTo(0, 0).Rectangle(Lb, h).Face(); beam.faces.name = "beam"
    beam.edges.Min(X).name = "clamp"; beam.edges.Max(Y).name = "loaded"
    meshB = Mesh(OCCGeometry(beam, dim=2).GenerateMesh(maxh=h / 3))
    u = solve_linear_elasticity(meshB, Emod, nu, dirichlet="clamp",
                                traction={"loaded": (0.0, -P_fem)}, plane="stress", order=2)
    tip_fem = abs(u(meshB(Lb, h / 2))[1])
    tip_eb = cantilever_tip_deflection(P_fem, Lb, Emod, h ** 3 / 12.0)
    assert abs(tip_fem - tip_eb) / tip_eb < 0.03


def validate_electrothermal_sigmaT_2way():
    """ref-class #24: TWO-WAY temperature-dependent sigma(T) electro-thermal.
    sigma(T)=sigma0/(1+alpha T) -> Joule q=(J^2/sigma0)(1+alpha T) feeds back into the heat;
    solve_heat_steady_nonlinear (Picard) matches the exact closed form
    T_max=(1/alpha)(sec(sqrt(b)L/2)-1), b=alpha J^2/(sigma0 k), and shows the feedback raising
    the peak above the constant-sigma parabola."""
    from netgen.occ import OCCGeometry, MoveTo, X
    from radia_mcp.radia_ngsolve.multiphysics import solve_heat_steady_nonlinear
    L, hbar, sigma0, k, alpha, J = 0.10, 0.01, 1.0e6, 50.0, 0.005, 8.0e5
    s = math.sqrt(alpha * J * J / (sigma0 * k))
    Tmax_cf = (1.0/alpha)*(1.0/math.cos(s*L/2.0) - 1.0)
    Tmax_parab = J*J*L*L/(8.0*sigma0*k)
    face = MoveTo(0, 0).Rectangle(L, hbar).Face(); face.faces.name = "bar"
    face.edges.Min(X).name = "ends"; face.edges.Max(X).name = "ends"
    mesh = Mesh(OCCGeometry(face, dim=2).GenerateMesh(maxh=L/60))
    T = solve_heat_steady_nonlinear(mesh, lambda Tcf: (J*J/sigma0)*(1.0 + alpha*Tcf),
                                    k, conductor="bar", dirichlet="ends", order=2)
    Tmax = T(mesh(L/2, hbar/2))
    print(f"[sigma(T) 2-way] T_max FE={Tmax:.4f} K  closed={Tmax_cf:.4f} K  "
          f"(parabola={Tmax_parab:.4f})  err={100*abs(Tmax-Tmax_cf)/Tmax_cf:.3f}%")
    assert abs(Tmax - Tmax_cf) / Tmax_cf < 0.01, "FE vs closed-form sigma(T) mismatch"
    assert Tmax > 1.03 * Tmax_parab, "2-way feedback should raise the peak above the parabola"


def validate_coax_line_inductance():
    """ref-class #25: coaxial-line external inductance from the FE magnetic energy
    matches (mu0/2pi) ln(b/a); with the analytic coax C it gives Z0 = 60 ln(b/a) and v=c."""
    from netgen.occ import OCCGeometry, WorkPlane, Glue
    from radia_mcp.radia_ngsolve.solve import reluctivity
    from radia_mcp.radia_ngsolve.force import magnetic_energy_2d
    EPS0 = 8.8541878128e-12
    a, b, cc, Rout, I = 0.002, 0.008, 0.010, 0.030, 1.0
    da = WorkPlane().Circle(0, 0, a).Face(); da.faces.name = "inner"
    db = WorkPlane().Circle(0, 0, b).Face(); diel = db - da; diel.faces.name = "diel"
    dc = WorkPlane().Circle(0, 0, cc).Face(); shield = dc - db; shield.faces.name = "shield"
    box = WorkPlane().Circle(0, 0, Rout).Face(); box.edges.name = "outer"
    air = box - dc; air.faces.name = "air"
    da.faces.maxh = a/8; diel.faces.maxh = (b-a)/24; shield.faces.maxh = (cc-b)/4; air.faces.maxh = Rout/8
    mesh = Mesh(OCCGeometry(Glue([da, diel, shield, air]), dim=2).GenerateMesh(maxh=Rout/6))
    Jz = CoefficientFunction([{"inner": I/(math.pi*a*a),
                              "shield": -I/(math.pi*(cc*cc-b*b))}.get(m, 0.0)
                             for m in mesh.GetMaterials()])
    A = solve_planar_magnetostatic(mesh, reluctivity(mesh, {}), Jz=Jz, order=3, dirichlet="outer")
    B = CoefficientFunction((grad(A)[1], -grad(A)[0]))
    L_fe = 2.0 * magnetic_energy_2d(B, mesh, "diel") / (I*I)
    L_cf = MU0/(2*math.pi)*math.log(b/a)
    C_cf = 2*math.pi*EPS0/math.log(b/a)
    Z0_fe, Z0_cf = math.sqrt(L_fe/C_cf), 60.0*math.log(b/a)
    print(f"[coax] L_ext FE={L_fe:.4e} closed={L_cf:.4e} ({100*abs(L_fe-L_cf)/L_cf:.2f}%)  "
          f"Z0 {Z0_fe:.2f} vs {Z0_cf:.2f} ohm")
    assert abs(L_fe - L_cf)/L_cf < 0.01, "coax L_ext vs (mu0/2pi)ln(b/a)"
    assert abs(Z0_fe - Z0_cf)/Z0_cf < 0.01, "coax Z0 vs 60 ln(b/a)"


def main():
    for name, fn in sorted(globals().items()):
        if name.startswith("validate_") and callable(fn):
            fn()
            print("ok", name)


if __name__ == "__main__":
    main()
