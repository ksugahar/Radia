"""Golden test (productionization M1): radia.vim.hdiv_demag_solve -- the consolidated LINEAR
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
        res = hdiv_demag_solve(mesh, 100.0, _HEXT)
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
        auto = hdiv_demag_solve(mesh, mu_r, _HEXT)                          # default -> symmetric C++ CG
        cg = hdiv_demag_solve(mesh, mu_r, _HEXT, linear_solver="cpp-cg")    # explicit alias
        gm = hdiv_demag_solve(mesh, mu_r, _HEXT, linear_solver="gmres")     # GMRES cross-check
    assert auto["linear_solver"] == "mass-riesz-cg"
    assert cg["linear_solver"] == "mass-riesz-cg"
    assert gm["linear_solver"] == "mass-riesz-gmres"
    rel = abs(auto["M_avg"][2] - gm["M_avg"][2]) / abs(gm["M_avg"][2])
    assert rel < 1e-6, f"symmetric CG vs GMRES M_avg disagree: {rel:.2e}"
    # Two independent TaskManager/PARDISO/CG runs of the same alias can differ by last-bit reduction order.
    assert abs(auto["M_avg"][2] - cg["M_avg"][2]) < 5e-9, "auto must equal the explicit cpp-cg alias"


def test_explicit_hlu_linear_solver():
    """The opt-in system-A H-LU path solves the production linear entry and reports its C++ solver."""
    mesh = _sphere(h=0.6)
    with ng.TaskManager():
        res = hdiv_demag_solve(mesh, 100.0, _HEXT, linear_solver="hlu", gram_eps=1e-8)
    assert res["linear_solver"] == "cpp-hlu"
    assert abs(res["demag"] - 1.0 / 3.0) < 7e-3
    assert res["iters"] == 1


def test_fail_loud_on_nonmagnetic():
    """mu_r <= 1 is not a soft-iron demag problem -> RAISE (no silent fallback)."""
    mesh = _sphere()
    with pytest.raises(ValueError):
        with ng.TaskManager():
            hdiv_demag_solve(mesh, 1.0, _HEXT)


@pytest.mark.parametrize("order", [1, 2])
def test_highorder_linear_solve_matches_rt0(order):
    """High-order (order>0) LINEAR material solve is now production-ready (per-element change-of-basis fix,
    2026-06-28): the order-p demag operator is valid (eig in [0,1]) and the material solve p-converges -- it
    agrees with the order-0 (RT0) volume-average M on the uniform-field sphere (the old ~2-4x blow-up is gone)."""
    mesh = _sphere(h=0.5)
    with ng.TaskManager():
        r0 = hdiv_demag_solve(mesh, 100.0, _HEXT, order=0)
        rp = hdiv_demag_solve(mesh, 100.0, _HEXT, order=order)
    assert rp["order"] == order and rp["ndof"] > r0["ndof"]
    assert abs(rp["demag"] - 1.0 / 3.0) < 1e-2, f"order={order} demag {rp['demag']:.4f} not ~1/3"
    rel = abs(rp["M_avg"][2] - r0["M_avg"][2]) / abs(r0["M_avg"][2])
    assert rel < 0.1, f"order={order} M_avg {rp['M_avg'][2]:.0f} vs RT0 {r0['M_avg'][2]:.0f} (rel {rel:.2f})"
    # higher order resolves more -> M_avg >= RT0 (converges up from below toward the continuum limit)
    assert rp["M_avg"][2] >= r0["M_avg"][2] - 1.0


def test_order_gt0_unsupported_combos_fail_loud():
    """order>0 wires the LINEAR (uniform / per-region) case; the not-yet-validated combos (nonlinear /
    image / HLU) must RAISE, not silently fall back (No-Fallbacks)."""
    mesh = _sphere(h=0.7)
    for kw in (dict(bh_table=[[0.0, 0.0], [1e6, 2.0]]),       # nonlinear at order>0
               dict(mu_r=100.0, image="+x"),                  # image symmetry at order>0
               dict(mu_r=100.0, linear_solver="hlu")):        # HLU is RT0-only
        with pytest.raises(NotImplementedError):
            with ng.TaskManager():
                hdiv_demag_solve(mesh, H_ext=_HEXT, order=1, **kw)


def test_order_gt2_fail_loud():
    """order>=3 needs the Duffy singular quadrature (the analytic-moment potential is exact only to charge
    degree 2 / order<=2) -- it must RAISE, not silently return a wrong M (No-Fallbacks)."""
    mesh = _sphere(h=0.7)
    with pytest.raises(NotImplementedError):
        with ng.TaskManager():
            hdiv_demag_solve(mesh, 100.0, _HEXT, order=3)


def test_requires_exactly_one_material_spec():
    """Exactly one of mu_r (linear) / bh_table (nonlinear) -> else RAISE."""
    mesh = _sphere()
    with ng.TaskManager():
        with pytest.raises(ValueError):                       # neither
            hdiv_demag_solve(mesh, H_ext=_HEXT)
        with pytest.raises(ValueError):                       # both
            hdiv_demag_solve(mesh, 100.0, _HEXT, bh_table=[[0.0, 0.0], [1e6, 2.0]])


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
        res = hdiv_demag_solve(mesh, bh_table=BH, H_ext=ng.CoefficientFunction((0, 0, H0)))
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
