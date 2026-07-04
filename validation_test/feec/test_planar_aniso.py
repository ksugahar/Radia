"""Golden lock for the SHARED 2D anisotropic-susceptibility demag solver (radia.planar_aniso).

M = X.H with X a uniaxial tensor (GO steel): the dense demag operator N is assembled on the SHARED
planar_charges kernel and (I - X N) M = X H0 solved directly (well-conditioned for any chi, unlike a
matrix-free Picard).  Gates: isotropic special case == exact Moment2DSolveLinear; anisotropic disk ==
analytic (I + D X)^-1 X H0 (D = 1/2); per-region (multi-grade); fail-loud.
"""
import numpy as np
import pytest

ng = pytest.importorskip("ngsolve")
pytest.importorskip("netgen")
from netgen.geom2d import SplineGeometry

import radia._radia_pybind as _rp
import radia.mmmm2d as m2
import radia.planar_aniso as pa
import radia.planar_materials as pm

MU0 = 4e-7 * np.pi


def _disk(maxh=0.2):
    g = SplineGeometry(); g.AddCircle((0, 0), r=1.0, bc="e")
    return ng.Mesh(g.GenerateMesh(maxh=maxh))


def test_isotropic_matches_exact_moment_solver():
    """X = chi I reproduces the exact scalar Moment2DSolveLinear demag (Gauss-N self-term accuracy)."""
    with ng.TaskManager():
        d = _disk(0.2)
        verts, offsets, centroids, areas = m2._extract_geometry(d)
        w = areas / areas.sum()
        for chi in (10.0, 100.0):
            r = pa.solve_anisotropic_demag(d, chi_par=chi, chi_perp=chi, H0=(1.0, 0.0))
            Mx_ex = w @ _rp.Moment2DSolveLinear(verts, offsets, np.full(len(areas), chi),
                                                np.tile([1.0, 0.0], (len(areas), 1)))[:, 0]
            assert abs(r["M_avg"][0] - Mx_ex) / Mx_ex < 5e-4, (chi, r["M_avg"][0], Mx_ex)


def test_anisotropic_disk_matches_analytic():
    """Anisotropic disk: M_avg == (I + D X)^-1 X H0, D = 1/2 (incl. tilted-easy-axis cross-magnetisn)."""
    D = 0.5
    with ng.TaskManager():
        d = _disk(0.18)
        for cpar, cperp, easy, H0 in [(50, 5, 0, [1, 0]), (50, 5, 0, [0, 1]),
                                      (50, 5, 45, [1, 0]), (30, 10, 30, [1, 0]), (200, 20, 60, [1, 0])]:
            r = pa.solve_anisotropic_demag(d, cpar, cperp, easy, H0=H0)
            X = pm.chi_tensor(cpar, cperp, easy)
            M_ana = np.linalg.solve(np.eye(2) + D * X, X @ np.asarray(H0, float))
            rel = np.linalg.norm(r["M_avg"] - M_ana) / np.linalg.norm(M_ana)
            assert rel < 2e-3, (cpar, cperp, easy, H0, r["M_avg"], M_ana, rel)


def test_easy_axis_alignment():
    """A disk magnetises MORE along the easy axis than across it.  Use LOW chi so the D=1/2 demag
    does not saturate M (at high chi both directions -> M~2, masking the anisotropy in M_avg); at low
    chi  M_easy/M_hard -> chi_par/chi_perp."""
    cpar, cperp = 8.0, 0.5
    with ng.TaskManager():
        d = _disk(0.2)
        r_easy = pa.solve_anisotropic_demag(d, cpar, cperp, easy_deg=0.0, H0=(1.0, 0.0))   # H along easy
        r_hard = pa.solve_anisotropic_demag(d, cpar, cperp, easy_deg=90.0, H0=(1.0, 0.0))  # H across easy
    # analytic disk (D=1/2): M_easy = cpar/(1+cpar/2), M_hard = cperp/(1+cperp/2)
    assert r_easy["M_avg"][0] > 3 * r_hard["M_avg"][0], (r_easy["M_avg"][0], r_hard["M_avg"][0])
    assert np.isclose(r_easy["M_avg"][0], cpar / (1 + cpar / 2), rtol=2e-3)
    assert np.isclose(r_hard["M_avg"][0], cperp / (1 + cperp / 2), rtol=2e-3)


def test_per_region_grades():
    """Two anisotropic grades on one mesh (concentric): the dict path builds per-region X."""
    with ng.TaskManager():
        g = SplineGeometry()
        g.AddCircle((0, 0), r=1.0, leftdomain=1, rightdomain=0, bc="ao")
        g.AddCircle((0, 0), r=0.5, leftdomain=2, rightdomain=1, bc="bo")
        g.SetMaterial(1, "a"); g.SetMaterial(2, "b")
        mesh = ng.Mesh(g.GenerateMesh(maxh=0.18))
        r = pa.solve_anisotropic_demag(mesh, chi_par={"a": 200.0, "b": 50.0},
                                       chi_perp={"a": 20.0, "b": 10.0}, easy_deg=0.0, H0=(1.0, 0.0))
    assert r["M_avg"][0] > 0 and abs(r["M_avg"][1]) < 0.05 * r["M_avg"][0]
    with pytest.raises(ValueError):
        pa.solve_anisotropic_demag(mesh, chi_par={"a": 200.0}, chi_perp=20.0, H0=(1.0, 0.0))  # 'b' missing


def test_bad_chi_fail_loud():
    with ng.TaskManager():
        d = _disk(0.3)
        with pytest.raises(ValueError):
            pa.solve_anisotropic_demag(d, chi_par=-1.0, chi_perp=10.0, H0=(1.0, 0.0))


def _iron_pm_mesh(maxh=0.16):
    """iron disk at x=+1.6 + PM disk at x=-1.6 (two regions, one mesh) -- design-B geometry."""
    g = SplineGeometry()
    g.AddCircle((1.6, 0.0), r=1.0, leftdomain=1, rightdomain=0, bc="ie")
    g.AddCircle((-1.6, 0.0), r=1.0, leftdomain=2, rightdomain=0, bc="pe")
    g.SetMaterial(1, "iron"); g.SetMaterial(2, "pm")
    return ng.Mesh(g.GenerateMesh(maxh=maxh))


def test_design_b_isotropic_matches_mmmm():
    """pm= (embedded PM) in the ISOTROPIC limit (chi_par=chi_perp) == the validated mmmm2d design-B
    (which itself matches a monolithic magnetostatic FEM) -- to the Gauss-N self-term accuracy."""
    MREM, chi = 8.0e5, 199.0
    with ng.TaskManager():
        mesh = _iron_pm_mesh(0.16)
        rA = pa.solve_anisotropic_demag(mesh, chi_par=chi, chi_perp=chi, H0=(0.0, 0.0),
                                        pm={"pm": [MREM, 0.0]})
        rM = m2.solve_planar_demag(mesh, mu_r={"iron": 1.0 + chi}, H_ext=(0.0, 0.0),
                                   pm={"pm": [MREM, 0.0]})
    iron_ids = np.array([i for i, m in enumerate(m2._element_materials(mesh)) if m == "iron"], int)
    MA = rA["M"][iron_ids].mean(axis=0)
    MM = rM["M"][iron_ids].mean(axis=0)
    assert rA["pm"] is True and abs(MA[0] - MM[0]) / abs(MM[0]) < 3e-3, (MA, MM)


def test_design_b_anisotropic_pm():
    """An ANISOTROPIC iron rotor with an embedded PM: the iron magnetises (novel combo neither the
    scalar MMMM design-B nor the HDiv-VIM had); PM stays pinned."""
    MREM = 8.0e5
    with ng.TaskManager():
        mesh = _iron_pm_mesh(0.16)
        r = pa.solve_anisotropic_demag(mesh, chi_par={"iron": 400.0}, chi_perp={"iron": 20.0},
                                       easy_deg={"iron": 0.0}, H0=(0.0, 0.0), pm={"pm": [MREM, 0.0]})
    pm_ids = np.array([i for i, m in enumerate(m2._element_materials(mesh)) if m == "pm"], int)
    iron_ids = np.array([i for i, m in enumerate(m2._element_materials(mesh)) if m == "iron"], int)
    assert np.allclose(r["M"][pm_ids], [MREM, 0.0])          # PM pinned
    assert r["M"][iron_ids].mean(axis=0)[0] > 0              # iron magnetised toward the magnet
