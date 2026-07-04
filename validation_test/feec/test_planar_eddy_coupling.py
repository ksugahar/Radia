"""Golden lock for the SHARED 2D staggered eddy-current coupling (radia.planar_eddy).

Weak coupling of the analytic MMMM soft-iron demag (radia.mmmm2d) with an NGSolve reduced-potential
complex A_z eddy FEM (Chadebec IEM + Biro reduced potential; the 2D maglev / IM method).  Gates:

 1. STANDALONE EDDY vs analytic conducting-cylinder Bessel  <Bx>/B0 = 2 I1(z)/(z I0(z)),  z=(1+j)a/d.
 2. STAGGERED couple_mmmm reproduces a MONOLITHIC FEM (iron+conductor+air meshed together, direct
    complex solve) on the iron magnetisation M_avg and the conductor-averaged Bx.
 3. sigma->0 limit: the coupled iron M collapses to the pure magnetostatic MMMM demag.

De-risk record: memory mmmm-2d-planar-tri-quad (E1/E2), C:\\temp\\mmmm2d_eddy.
"""
import numpy as np
import pytest

ng = pytest.importorskip("ngsolve")
pytest.importorskip("netgen")
sp = pytest.importorskip("scipy.special")
from netgen.geom2d import SplineGeometry

import radia.mmmm2d as m2
import radia.planar_eddy as pe

MU0 = 4e-7 * np.pi
SIGMA = 3.7e7
A = 0.01
R = 40 * A


def _disk_mesh(c, maxh):
    geo = SplineGeometry(); geo.AddCircle(c, r=A, bc="e")
    return ng.Mesh(geo.GenerateMesh(maxh=maxh))


def _cond_air_mesh(c_cond, maxh_c, maxh_a=R / 12):
    geo = SplineGeometry()
    geo.AddCircle((0, 0), r=R, leftdomain=1, rightdomain=0, bc="outer")
    geo.AddCircle(c_cond, r=A, leftdomain=2, rightdomain=1, bc="cond_ifc")
    geo.SetMaterial(1, "air"); geo.SetMaterial(2, "conductor")
    geo.SetDomainMaxH(1, maxh_a); geo.SetDomainMaxH(2, maxh_c)
    return ng.Mesh(geo.GenerateMesh(maxh=maxh_a))


def _two_body_mesh(c_iron, c_cond, maxh_i, maxh_c, maxh_a=R / 12):
    geo = SplineGeometry()
    geo.AddCircle((0, 0), r=R, leftdomain=1, rightdomain=0, bc="outer")
    geo.AddCircle(c_iron, r=A, leftdomain=2, rightdomain=1, bc="iron_ifc")
    geo.AddCircle(c_cond, r=A, leftdomain=3, rightdomain=1, bc="cond_ifc")
    for i, m in enumerate(("air", "iron", "conductor"), start=1):
        geo.SetMaterial(i, m)
    geo.SetDomainMaxH(1, maxh_a); geo.SetDomainMaxH(2, maxh_i); geo.SetDomainMaxH(3, maxh_c)
    return ng.Mesh(geo.GenerateMesh(maxh=maxh_a))


def _monolithic(mesh, mu_r, sigma, freq, B0=1.0, order=4):
    """Reference: iron mu_r + conductor sigma + air, direct complex A_z solve, A_z=B0 y on outer."""
    w = 2 * np.pi * freq
    mesh.Curve(order)
    fes = ng.H1(mesh, order=order, complex=True, dirichlet="outer")
    u, v = fes.TnT()
    nu = ng.CoefficientFunction([1.0 / (mu_r * MU0) if m == "iron" else 1.0 / MU0
                                 for m in mesh.GetMaterials()])
    sig = ng.CoefficientFunction([sigma if m == "conductor" else 0.0 for m in mesh.GetMaterials()])
    af = ng.BilinearForm(fes, symmetric=True)
    af += nu * ng.grad(u) * ng.grad(v) * ng.dx
    af += 1j * w * sig * u * v * ng.dx
    af.Assemble()
    gfu = ng.GridFunction(fes)
    gfu.Set(B0 * ng.y, definedon=mesh.Boundaries("outer"))
    r = gfu.vec.CreateVector(); r.data = -af.mat * gfu.vec
    gfu.vec.data += af.mat.Inverse(fes.FreeDofs(), inverse="pardiso") * r
    return gfu, mesh


def _avg(gfu, mesh, mat, comp):
    area = ng.Integrate(ng.CF(1), mesh, definedon=mesh.Materials(mat))
    g = ng.grad(gfu)[1] if comp == "x" else -ng.grad(gfu)[0]
    return ng.Integrate(g, mesh, definedon=mesh.Materials(mat)) / area


def test_standalone_eddy_matches_bessel():
    """The reduced-Az eddy FEM (no iron) reproduces the analytic transverse-cylinder Bessel field."""
    ratio = 1.5
    delta = A / ratio
    freq = (2.0 / (MU0 * SIGMA * delta ** 2)) / (2 * np.pi)
    z = (1 + 1j) * A / delta
    bessel = 2.0 * sp.iv(1, z) / (z * sp.iv(0, z))
    with ng.TaskManager():
        fem = _cond_air_mesh((0.0, 0.0), maxh_c=min(A / 6, delta / 3))
        # couple with a mu_r->1 iron (no magnetisation) == standalone eddy
        iron = _disk_mesh((-3 * A, 0.0), maxh=A / 5)
        res = pe.couple_mmmm(iron, mu_r=1.0 + 1e-12, fem_mesh=fem, sigma=SIGMA, freq=freq)
        Bx = 1.0 + ng.Integrate(ng.grad(res["gfu"])[1], fem, definedon=fem.Materials("conductor")) \
            / ng.Integrate(ng.CF(1), fem, definedon=fem.Materials("conductor"))
    assert abs(Bx - bessel) / abs(bessel) < 3e-3, (Bx, bessel)


def test_staggered_matches_monolithic():
    """couple_mmmm == monolithic FEM on iron M_avg and conductor <Bx> (the coupling gate)."""
    mu_r = 100.0
    ratio = 1.5
    delta = A / ratio
    freq = (2.0 / (MU0 * SIGMA * delta ** 2)) / (2 * np.pi)
    c_iron, c_cond = (-1.6 * A, 0.0), (1.6 * A, 0.0)
    with ng.TaskManager():
        # monolithic reference
        mm = _two_body_mesh(c_iron, c_cond, maxh_i=A / 8, maxh_c=min(A / 6, delta / 3))
        gfu_m, mm = _monolithic(mm, mu_r, SIGMA, freq)
        Bx_i, By_i = _avg(gfu_m, mm, "iron", "x"), _avg(gfu_m, mm, "iron", "y")
        M_mono = (mu_r - 1.0) / (mu_r * MU0) * np.array([Bx_i, By_i])
        Bx_c_mono = _avg(gfu_m, mm, "conductor", "x")
        # staggered
        iron = _disk_mesh(c_iron, maxh=A / 8)
        fem = _cond_air_mesh(c_cond, maxh_c=min(A / 6, delta / 3))
        res = pe.couple_mmmm(iron, mu_r=mu_r, fem_mesh=fem, sigma=SIGMA, freq=freq)
    assert res["iters"] <= 12, res["hist"]
    relM = abs(res["M_avg"][0] - M_mono[0]) / abs(M_mono[0])
    assert relM < 3e-3, (res["M_avg"][0], M_mono[0], relM)
    # conductor <Bx>: staggered total = B0 + <iron field>_cond + <eddy reaction>_cond
    from radia.planar_charges import exterior_field
    rng = np.random.default_rng(0)
    r = A * np.sqrt(rng.random(400)); th = 2 * np.pi * rng.random(400)
    P = np.stack([c_cond[0] + r * np.cos(th), c_cond[1] + r * np.sin(th)], axis=1)
    with ng.TaskManager():
        Hi = exterior_field(iron, res["M"].real, P) + 1j * exterior_field(iron, res["M"].imag, P)
        Br = ng.Integrate(ng.grad(res["gfu"])[1], fem, definedon=fem.Materials("conductor")) \
            / ng.Integrate(ng.CF(1), fem, definedon=fem.Materials("conductor"))
    Bx_c_stag = 1.0 + MU0 * Hi[:, 0].mean() + Br
    assert abs(Bx_c_stag - Bx_c_mono) / abs(Bx_c_mono) < 3e-3, (Bx_c_stag, Bx_c_mono)


def test_sigma_to_zero_recovers_magnetostatic_demag():
    """As sigma->0 the coupled iron M collapses to the pure magnetostatic MMMM demag (no eddy)."""
    mu_r = 50.0
    c_iron, c_cond = (-1.6 * A, 0.0), (1.6 * A, 0.0)
    with ng.TaskManager():
        iron = _disk_mesh(c_iron, maxh=A / 8)
        fem = _cond_air_mesh(c_cond, maxh_c=A / 6)
        res = pe.couple_mmmm(iron, mu_r=mu_r, fem_mesh=fem, sigma=1.0, freq=1e-3)  # ~no eddy
        pure = m2.solve_planar_demag(iron, mu_r=mu_r, H_ext=(1.0 / MU0, 0.0))
    # coupled M (real part) matches the standalone magnetostatic demag
    assert abs(res["M_avg"][0].real - pure["M_avg"][0]) / abs(pure["M_avg"][0]) < 1e-3, \
        (res["M_avg"][0], pure["M_avg"][0])
    assert abs(res["M_avg"][0].imag) < 1e-3 * abs(res["M_avg"][0].real)      # eddy phase ~0


# ---- unified PM + soft-iron + eddy rotor (PM-motor / eddy-current brake) -------------------------
MREM = 8.0e5
C_ROT, C_COND = (-2.0 * A, 0.0), (2.0 * A, 0.0)


def _rotor_mesh(maxh):
    """PM core (0..0.5a) inside a soft-iron annulus (0.5a..a) -- a real PM+iron rotor, no air."""
    geo = SplineGeometry()
    geo.AddCircle(C_ROT, r=A, leftdomain=1, rightdomain=0, bc="iron_o")
    geo.AddCircle(C_ROT, r=0.5 * A, leftdomain=2, rightdomain=1, bc="pm_o")
    geo.SetMaterial(1, "iron"); geo.SetMaterial(2, "pm")
    return ng.Mesh(geo.GenerateMesh(maxh=maxh))


def _rotor_cond_air(maxh_body, maxh_c, maxh_a=R / 12):
    geo = SplineGeometry()
    geo.AddCircle((0, 0), r=R, leftdomain=1, rightdomain=0, bc="outer")
    geo.AddCircle(C_ROT, r=A, leftdomain=2, rightdomain=1, bc="iron_o")
    geo.AddCircle(C_ROT, r=0.5 * A, leftdomain=3, rightdomain=2, bc="pm_o")
    geo.AddCircle(C_COND, r=A, leftdomain=4, rightdomain=1, bc="cond_ifc")
    for i, m in enumerate(("air", "iron", "pm", "conductor"), start=1):
        geo.SetMaterial(i, m)
    geo.SetDomainMaxH(2, maxh_body); geo.SetDomainMaxH(3, maxh_body); geo.SetDomainMaxH(4, maxh_c)
    return ng.Mesh(geo.GenerateMesh(maxh=maxh_a))


def _monolithic_pm_ac(mesh, mu_r, sigma, freq, Mpm, B0=1.0, order=4):
    """Monolithic AC+PM: int nu grad(A).grad(v) + j w sigma int_c A v = int_pm (Mx dv/dy - My dv/dx)."""
    w = 2 * np.pi * freq
    mesh.Curve(order)
    nu = ng.CoefficientFunction([1.0 / (mu_r * MU0) if m == "iron" else 1.0 / MU0
                                 for m in mesh.GetMaterials()])
    sig = ng.CoefficientFunction([sigma if m == "conductor" else 0.0 for m in mesh.GetMaterials()])
    fes = ng.H1(mesh, order=order, complex=True, dirichlet="outer")
    u, v = fes.TnT()
    af = ng.BilinearForm(fes, symmetric=True)
    af += nu * ng.grad(u) * ng.grad(v) * ng.dx
    af += 1j * w * sig * u * v * ng.dx
    af.Assemble()
    gfu = ng.GridFunction(fes)
    gfu.Set(B0 * ng.y, definedon=mesh.Boundaries("outer"))
    lf = ng.LinearForm(fes)
    lf += (Mpm[0] * ng.grad(v)[1] - Mpm[1] * ng.grad(v)[0]) * ng.dx("pm")
    lf.Assemble()
    r = lf.vec.CreateVector(); r.data = lf.vec - af.mat * gfu.vec
    gfu.vec.data += af.mat.Inverse(fes.FreeDofs(), inverse="pardiso") * r
    return gfu, mesh


def test_pm_eddy_unified_rotor():
    """A PM+iron rotor coupled to a conductor eddy (couple_mmmm pm=) == monolithic AC+PM FEM."""
    mu_r = 100.0
    delta = A / 1.5
    freq = (2.0 / (MU0 * SIGMA * delta ** 2)) / (2 * np.pi)
    with ng.TaskManager():
        fm = _rotor_cond_air(maxh_body=A / 8, maxh_c=min(A / 6, delta / 3))
        gfu_m, fm = _monolithic_pm_ac(fm, mu_r, SIGMA, freq, [MREM, 0.0])
        M_mono = (mu_r - 1.0) / (mu_r * MU0) * np.array([_avg(gfu_m, fm, "iron", "x"),
                                                         _avg(gfu_m, fm, "iron", "y")])
        rot = _rotor_mesh(maxh=A / 8)
        fem = _cond_air_mesh(C_COND, maxh_c=min(A / 6, delta / 3))
        res = pe.couple_mmmm(rot, mu_r=mu_r, fem_mesh=fem, sigma=SIGMA, freq=freq,
                             pm={"pm": [MREM, 0.0]})
    soft = np.array([i for i, m in enumerate(m2._element_materials(rot)) if m == "iron"], int)
    M_iron = res["M"][soft].mean(axis=0)
    assert res["iters"] <= 12, res["hist"]
    rel = abs(M_iron[0] - M_mono[0]) / abs(M_mono[0])
    assert rel < 1e-2, (M_iron[0], M_mono[0], rel)           # PM magnetises iron + drives eddy
    assert abs(M_iron[0].imag) > 1e-3 * abs(M_iron[0].real)  # eddy induces a real phase lag


# ---- nonlinear soft iron + eddy (effective-chi AC) ----------------------------------------------
BH = [[0.0, 0.0], [200.0, 0.30], [800.0, 1.20], [3000.0, 1.70], [20000.0, 2.00]]   # saturating iron


def test_nonlinear_sigma_to_zero_recovers_dc_demag():
    """couple(bh_table) at sigma->0 collapses to the standalone DC nonlinear MMMM demag (same law)."""
    B0 = 0.6
    with ng.TaskManager():
        iron = _disk_mesh((-1.6 * A, 0.0), maxh=A / 8)
        fem = _cond_air_mesh((1.6 * A, 0.0), maxh_c=A / 6)
        res = pe.couple_mmmm(iron, fem, sigma=1.0, freq=1e-3, bh_table=BH, B0=B0)   # ~no eddy
        pure = m2.solve_planar_demag(iron, bh_table=BH, H_ext=(B0 / MU0, 0.0))
    rel = abs(res["M_avg"][0].real - pure["M_avg"][0]) / abs(pure["M_avg"][0])
    assert rel < 2e-3, (res["M_avg"][0], pure["M_avg"][0], rel)
    assert abs(res["M_avg"][0].imag) < 2e-3 * abs(res["M_avg"][0].real)      # eddy phase ~0


def test_nonlinear_low_drive_recovers_linear():
    """At a low drive the saturating iron stays in its initial-permeability (chi0) linear regime."""
    from radia.mmmm2d import _hm_arrays
    _, _, chi0 = _hm_arrays(BH)
    delta = A / 1.0
    freq = (2.0 / (MU0 * SIGMA * delta ** 2)) / (2 * np.pi)
    B0 = 5e-5                                                 # H_app = B0/mu0 ~ 40 A/m << 200 (knee)
    with ng.TaskManager():
        iron = _disk_mesh((-1.6 * A, 0.0), maxh=A / 7)
        fem = _cond_air_mesh((1.6 * A, 0.0), maxh_c=min(A / 6, delta / 3))
        r_nl = pe.couple_mmmm(iron, fem, sigma=SIGMA, freq=freq, bh_table=BH, B0=B0)
        r_lin = pe.couple_mmmm(iron, fem, sigma=SIGMA, freq=freq, mu_r=1.0 + chi0, B0=B0)
    rel = abs(r_nl["M_avg"][0] - r_lin["M_avg"][0]) / abs(r_lin["M_avg"][0])
    assert rel < 1e-3, (r_nl["M_avg"][0], r_lin["M_avg"][0], rel)            # nonlinear -> linear chi0


def test_nonlinear_picard_fail_loud():
    """The nonlinear effective-chi Picard FAILS LOUD on non-convergence (No-Fallbacks -- never a
    silent unconverged M).  Force it with nl_maxit=1 + a tight tol under a strong saturating drive."""
    with ng.TaskManager():
        iron = _disk_mesh((0.0, 0.0), maxh=A / 6)
        solve = pe._mmmm_iron_solve(iron, bh_table=BH, nl_maxit=1, nl_tol=1e-14)
        H = np.tile([1e6, 0.0], (iron.ne, 1)).astype(complex)   # deep-saturating -> chi0 far from chi*
        with pytest.raises(RuntimeError, match="NOT converged"):
            solve(H)


# ---- per-region (multi-grade) nonlinear iron + eddy ---------------------------------------------
BH_SOFT = [[0.0, 0.0], [150.0, 0.60], [600.0, 1.30], [2500.0, 1.75], [20000.0, 2.05]]   # softer grade


def _two_grade_rotor(maxh):
    """A concentric two-grade iron rotor (at the module A scale, inside the fem air box): inner disk
    'iron_b' (0..0.5a) + outer annulus 'iron_a' (0.5a..a)."""
    geo = SplineGeometry()
    geo.AddCircle((-1.6 * A, 0), r=A, leftdomain=1, rightdomain=0, bc="ao")
    geo.AddCircle((-1.6 * A, 0), r=0.5 * A, leftdomain=2, rightdomain=1, bc="bo")
    geo.SetMaterial(1, "iron_a"); geo.SetMaterial(2, "iron_b")
    return ng.Mesh(geo.GenerateMesh(maxh=maxh))


def test_per_region_nonlinear_sigma_to_zero():
    """couple(bh_table={grade: table}) at sigma->0 == the standalone per-region DC nonlinear demag."""
    B0 = 0.7
    bh = {"iron_a": BH, "iron_b": BH_SOFT}
    with ng.TaskManager():
        rotor = _two_grade_rotor(maxh=A / 7)
        fem = _cond_air_mesh((1.6 * A, 0.0), maxh_c=A / 6)
        res = pe.couple_mmmm(rotor, fem, sigma=1.0, freq=1e-3, bh_table=bh, B0=B0)   # ~no eddy
        pure = m2.solve_planar_demag(rotor, bh_table=bh, H_ext=(B0 / MU0, 0.0))
    rel = abs(res["M_avg"][0].real - pure["M_avg"][0]) / abs(pure["M_avg"][0])
    assert rel < 2e-3, (res["M_avg"][0], pure["M_avg"][0], rel)
    assert abs(res["M_avg"][0].imag) < 2e-3 * abs(res["M_avg"][0].real)      # eddy phase ~0
