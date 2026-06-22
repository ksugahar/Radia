"""Golden test (productionization M1): radia.vim.hdiv_demag_solve -- the consolidated LINEAR
soft-iron applied-field demag solve (the candidate yano-type replacement for rad.Solve).

Locks:
  (1) PHYSICS: on a uniform-field sphere the volume-average M matches the analytic linear
      magnetization  M = chi/(1 + chi D) H_applied  (D = demag factor ~ 1/3) across mu_r 10..1000;
  (2) per-element M is ~uniform and +z-aligned (uniform applied field on a sphere);
  (3) the scalable path (H-matvec + approximate Jacobi diagonal, no dense N^2) matches the dense
      reference path;
  (4) fail-loud: mu_r <= 1 RAISES (CLAUDE.md No-Fallbacks).
"""
import numpy as np
import pytest

pytest.importorskip("ngsolve")
pytest.importorskip("netgen.csg")
import ngsolve as ng  # noqa: E402
from netgen.csg import CSGeometry, Sphere, Pnt  # noqa: E402

from radia.vim import hdiv_demag_solve  # noqa: E402


def _sphere(h=0.45):
    g = CSGeometry(); g.Add(Sphere(Pnt(0, 0, 0), 1.0))
    with ng.TaskManager():
        return ng.Mesh(g.GenerateMesh(maxh=h))


H0 = 1000.0  # uniform applied field, +z (A/m)
_HEXT = ng.CoefficientFunction((0, 0, H0))


@pytest.mark.parametrize("mu_r", [10.0, 100.0, 1000.0])
def test_sphere_linear_matches_analytic(mu_r):
    """Volume-average M == analytic chi/(1+chi D) H on the uniform-field sphere (dense path)."""
    mesh = _sphere()
    with ng.TaskManager():
        res = hdiv_demag_solve(mesh, mu_r, _HEXT)
    chi = mu_r - 1.0
    D = res["demag"]
    M_analytic = chi / (1.0 + chi * D) * H0
    assert abs(D - 1.0 / 3.0) < 5e-3, f"demag factor off: {D}"
    rel = abs(res["M_avg"][2] - M_analytic) / M_analytic
    assert rel < 2e-3, f"mu_r={mu_r}: Mz_avg {res['M_avg'][2]:.2f} vs analytic {M_analytic:.2f} (rel {rel:.2e})"
    # transverse average components are ~0 (no transverse drive)
    assert abs(res["M_avg"][0]) < 1e-3 * M_analytic and abs(res["M_avg"][1]) < 1e-3 * M_analytic


def test_per_element_M_uniform_and_aligned():
    """Per-element M on the uniform-field sphere is ~uniform and dominated by +z."""
    mesh = _sphere()
    with ng.TaskManager():
        res = hdiv_demag_solve(mesh, 100.0, _HEXT)
    M = res["M"]
    assert M.shape == (res["n_el"], 3)
    mz = M[:, 2]
    # +z dominant, low relative spread (RT0 on a coarse sphere -> a few % element-to-element)
    assert mz.mean() > 0 and np.std(mz) / abs(mz.mean()) < 0.1, \
        f"per-element Mz not uniform: mean {mz.mean():.1f}, std {np.std(mz):.1f}"
    assert np.abs(M[:, :2]).mean() < 0.05 * abs(mz.mean()), "spurious transverse per-element M"


def test_fail_loud_on_nonmagnetic():
    """mu_r <= 1 is not a soft-iron demag problem -> RAISE (no silent fallback)."""
    mesh = _sphere()
    with pytest.raises(ValueError):
        with ng.TaskManager():
            hdiv_demag_solve(mesh, 1.0, _HEXT)


def test_requires_exactly_one_material_spec():
    """Exactly one of mu_r (linear) / bh_table (nonlinear) -> else RAISE."""
    mesh = _sphere()
    with ng.TaskManager():
        with pytest.raises(ValueError):                       # neither
            hdiv_demag_solve(mesh, H_ext=_HEXT)
        with pytest.raises(ValueError):                       # both
            hdiv_demag_solve(mesh, 100.0, _HEXT, bh_table=[[0.0, 0.0], [1e6, 2.0]])


def test_nonlinear_sphere_vs_dense_newton():
    """NONLINEAR hdiv_demag_solve (uniform field + a real BH table) reproduces the dense reference
    Newton (radia.vim._nonlinear.solve_nonlinear_newton, analytic Gram) on the saturating sphere."""
    from radia.vim._nonlinear import solve_nonlinear_newton
    chi0, Msat, H0 = 1000.0, 1.0e6, 2.0e5
    Hs = np.concatenate([[0.0], np.logspace(-1, 7, 60)])
    Ms = chi0 * Hs / (1.0 + chi0 * Hs / Msat)
    Bs = (4e-7 * np.pi) * (Hs + Ms)
    BH = [[float(h), float(b)] for h, b in zip(Hs, Bs)]
    mesh = _sphere()
    with ng.TaskManager():
        res = hdiv_demag_solve(mesh, bh_table=BH, H_ext=ng.CoefficientFunction((0, 0, H0)))
        mz_dense, _nit, _D = solve_nonlinear_newton(mesh, chi0, Msat, H0,
                                                    bh_table=(Hs, Bs), analytic_gram=True, maxit=60)
    assert res["nonlinear"] is True
    rel = abs(res["M_avg"][2] - mz_dense) / abs(mz_dense)
    assert rel < 2e-2, f"nonlinear scalable Mz {res['M_avg'][2]:.1f} vs dense {mz_dense:.1f} rel {rel:.2e}"
    # the entry's nonlinear solver is the Anderson-Hantila fixed point (cheap step), so the outer count
    # is higher than the dense reference Newton's ~5-6 (a fixed point + saturation slowdown) but bounded
    # well inside nl_maxit; the accuracy match above is the real gate.
    assert res["iters"] < 250, f"nonlinear Anderson-Hantila outer iters {res['iters']}"
