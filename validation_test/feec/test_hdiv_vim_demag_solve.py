"""Golden test (productionization M1): radia.vim.Solve -- the consolidated LINEAR
soft-iron applied-field demag solve (the candidate six-face surface-charge replacement for rad.Solve).

Locks:
  (1) PHYSICS: on a uniform-field sphere the volume-average M matches the analytic linear
      magnetization  M = chi/(1 + chi D) H_applied  (D = demag factor ~ 1/3) across mu_r 10..1000;
  (2) per-element M is ~uniform and +z-aligned (uniform applied field on a sphere);
  (3) exactly one of mu_r (linear) / bh_table (nonlinear) must be given -- else RAISE;
  (4) fail-loud: mu_r <= 1 RAISES (CLAUDE.md No-Fallbacks);
  (5) NONLINEAR: a saturating BH table reproduces the analytic uniform-sphere fixed point
      M = Mof(H0 - D M), via the C++ charge-Gram matrix-free Newton.
"""
import numpy as np
import pytest

pytest.importorskip("ngsolve")
pytest.importorskip("netgen.csg")
import ngsolve as ng  # noqa: E402
from netgen.csg import CSGeometry, Sphere, Pnt  # noqa: E402

from radia.vim import Solve  # noqa: E402


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
        res = Solve(mesh, mu_r, _HEXT)
    assert res["linear_solver"] in {"cpp-hlu", "mass-riesz-cg", "mass-riesz-gmres"}  # default 'auto' = symmetric mass-riesz CG
    assert "hmat_stats" in res
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
        res = Solve(mesh, 100.0, _HEXT)
    M = res["M"]
    assert M.shape == (res["n_el"], 3)
    mz = M[:, 2]
    # +z dominant, low relative spread (RT0 on a coarse sphere -> a few % element-to-element)
    assert mz.mean() > 0 and np.std(mz) / abs(mz.mean()) < 0.1, \
        f"per-element Mz not uniform: mean {mz.mean():.1f}, std {np.std(mz):.1f}"
    assert np.abs(M[:, :2]).mean() < 0.05 * abs(mz.mean()), "spurious transverse per-element M"


@pytest.mark.parametrize("mu_r", [1e2, 1e5])
def test_default_symmetric_cg_matches_gmres(mu_r):
    """The default 'auto' is the all-C++ SYMMETRIC mass-Riesz CG (the symmetric-HACApK matvec makes CG
    mathematically valid); 'cpp-cg' is an explicit alias for it, and 'gmres' is the asymmetry-tolerant
    cross-check.  All three converge to the SAME magnetization (the symmetric Gram is a robustness +
    speed change, not an accuracy change)."""
    mesh = _sphere(h=0.5)
    with ng.TaskManager():
        auto = Solve(mesh, mu_r, _HEXT)                          # default -> symmetric C++ CG
        cg = Solve(mesh, mu_r, _HEXT, linear_solver="cpp-cg")    # explicit alias
        gm = Solve(mesh, mu_r, _HEXT, linear_solver="gmres")     # GMRES cross-check
    assert auto["linear_solver"] == "mass-riesz-cg"
    assert cg["linear_solver"] == "mass-riesz-cg"
    assert gm["linear_solver"] == "mass-riesz-gmres"
    rel = abs(auto["M_avg"][2] - gm["M_avg"][2]) / abs(gm["M_avg"][2])
    assert rel < 1e-6, f"symmetric CG vs GMRES M_avg disagree: {rel:.2e}"
    # Two independent TaskManager/PARDISO/CG runs of the same alias can differ by last-bit reduction order.
    assert abs(auto["M_avg"][2] - cg["M_avg"][2]) < 5e-9, "auto must equal the explicit cpp-cg alias"


# (test_explicit_hlu_linear_solver removed 2026-06-29: the 'hlu' system-A H-LU solver was RT0-only and is
#  retired -- HDiv-VIM (RT1) uses the symmetric mass-Riesz CG.  See test_hdiv_vim_rt1_contract.py.)


def test_fail_loud_on_nonmagnetic():
    """mu_r <= 1 is not a soft-iron demag problem -> RAISE (no silent fallback)."""
    mesh = _sphere()
    with pytest.raises(ValueError):
        with ng.TaskManager():
            Solve(mesh, 1.0, _HEXT)


def test_rt1_linear_solve_is_physically_sane():
    """The RT1 LINEAR material solve on the uniform-field sphere gives the sphere demag ~1/3 and a +z M_avg
    (no ~2-4x blow-up -- the per-element change-of-basis fix, 2026-06-28).  RT0 is retired; RT1 is the demag
    spectrum reference (eig in [0,1] locked by test_hdiv_vim_demag_spectrum_psd)."""
    mesh = _sphere(h=0.5)
    with ng.TaskManager():
        r = Solve(mesh, 100.0, _HEXT, order=1)
    assert r["order"] == 1
    assert abs(r["demag"] - 1.0 / 3.0) < 1e-2, f"RT1 demag {r['demag']:.4f} not ~1/3"
    assert r["M_avg"][2] > 0, f"RT1 M_avg_z {r['M_avg'][2]:.0f} should respond +z to the applied field"


# (test_order_gt0_unsupported_combos_fail_loud removed 2026-06-29: the retired-feature fail-loud guards
#  (curved image / hlu / gauss / pm_M / unsupported topology) are now locked by
#  test_hdiv_vim_rt1_contract.py.)


def test_order_gt1_rt2_abolished():
    """RT2+ (HDiv solution order >= 2) is abolished -- order>1 must RAISE ValueError (RT2 gave no per-element
    magnetization gain over RT1 and was slower; No-Fallbacks)."""
    mesh = _sphere(h=0.7)
    with pytest.raises(ValueError, match="RT2"):
        with ng.TaskManager():
            Solve(mesh, 100.0, _HEXT, order=2)


def test_requires_exactly_one_material_spec():
    """Exactly one of mu_r (linear) / bh_table (nonlinear) -> else RAISE."""
    mesh = _sphere()
    with ng.TaskManager():
        with pytest.raises(ValueError):                       # neither
            Solve(mesh, H_ext=_HEXT)
        with pytest.raises(ValueError):                       # both
            Solve(mesh, 100.0, _HEXT, bh_table=[[0.0, 0.0], [1e6, 2.0]])


def test_nonlinear_sphere_vs_analytic_fixed_point():
    """NONLINEAR hdiv_demag_solve (uniform field + a real BH table) reproduces the analytic uniform-sphere
    fixed point M = Mof(H0 - D M) (D = demag factor ~ 1/3) on the saturating sphere."""
    chi0, Msat, H0 = 1000.0, 1.0e6, 2.0e5
    Hs = np.concatenate([[0.0], np.logspace(-1, 7, 60)])
    Ms = chi0 * Hs / (1.0 + chi0 * Hs / Msat)
    Bs = (4e-7 * np.pi) * (Hs + Ms)
    BH = [[float(h), float(b)] for h, b in zip(Hs, Bs)]
    Mof = lambda H: chi0 * H / (1.0 + chi0 * abs(H) / Msat)   # noqa: E731 (the table's M(H))
    mesh = _sphere()
    with ng.TaskManager():
        res = Solve(mesh, bh_table=BH, H_ext=ng.CoefficientFunction((0, 0, H0)))
    assert res["nonlinear"] is True
    # analytic uniform-sphere fixed point with the solved demag factor D
    D = res["demag"]
    lo, hi = -Msat, Msat
    f = lambda M: M - Mof(H0 - D * M)                          # noqa: E731
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if f(lo) * f(mid) <= 0:
            hi = mid
        else:
            lo = mid
    M_fp = 0.5 * (lo + hi)
    rel = abs(res["M_avg"][2] - M_fp) / abs(M_fp)
    assert rel < 2e-2, f"nonlinear Mz {res['M_avg'][2]:.1f} vs analytic fixed point {M_fp:.1f} rel {rel:.2e}"
    # the nonlinear solver is the damped matrix-free Newton (C++ charge Gram); bounded well inside nl_maxit.
    assert res["iters"] < 300, f"nonlinear Newton outer iters {res['iters']}"
