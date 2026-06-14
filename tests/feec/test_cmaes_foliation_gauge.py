"""Golden lock: CMA-ES over the foliation gauge (the F1 remedy).

Runs examples/vim/cmaes_foliation_gauge.py and asserts:

  - the inner field-fit is exact at every gauge (the bilevel inner problem is
    CONVEX -> fit_max ~ 0);
  - the outer coil-complexity objective over foliations is MULTI-BASIN (the F1
    wall: >= 2 local minima);
  - CMA-ES (global, Optuna CmaEsSampler) finds ~the grid-global minimum;
  - a LOCAL optimizer started in a shallow basin stays STUCK at a clearly higher
    complexity -> the demonstration that CMA-ES is needed to escape the wall.

Marked slow: a foliation landscape scan + a CMA-ES run + a local descent.
"""
import os
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
EXDIR = os.path.join(HERE, "..", "..", "examples", "vim")


@pytest.mark.slow
def test_cmaes_foliation_gauge_golden():
    sys.path.insert(0, EXDIR)
    try:
        import cmaes_foliation_gauge as M
    except Exception as e:  # pragma: no cover - env-dependent (optuna/ngsolve)
        pytest.skip(f"CMA-ES foliation example not importable: {e}")

    res, ok = M.main(n_trials=40)
    assert ok, f"example self-check failed: {res}"

    # inner field-fit exact at every gauge (convex inner problem)
    assert res["fit_max"] < 1e-6, f"inner solve not exact: {res['fit_max']:.2e}"

    # outer objective is multi-basin -- the F1 wall
    assert res["n_local_minima"] >= 2, (
        f"foliation landscape not multi-basin: {res['n_local_minima']} minima")

    # CMA-ES finds ~ the global minimum
    assert res["cmaes_vs_grid"] <= 1.05, (
        f"CMA-ES did not reach the grid-global: ratio {res['cmaes_vs_grid']:.3f}")

    # a LOCAL optimizer is stuck clearly above the CMA-ES optimum (the point)
    assert res["local_vs_cmaes"] > 1.03, (
        f"local optimizer not stuck (no F1-wall demonstration): "
        f"g_local/g_cmaes = {res['local_vs_cmaes']:.3f}")
