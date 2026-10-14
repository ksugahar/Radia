"""Golden lock for the SF coil-designer cohomology net-current 'finish'.

Runs validation_test/stream_function/cohomology_net_current.py (torus b1=2) and asserts:

  - the analytic secular generators h_k = n x grad(angle) are div-free
    (a valid surface current) and NON-EXACT (genuine cohomology, not any psi);
  - the surface Euler characteristic (gmsh-free) confirms b1 = 2 (the count);
  - the net-poloidal (TF) generator obeys Ampere (B_tor * R const inside, ~0
    outside);
  - the JOINT design (prescribe net current -> fit psi to the residual)
    reproduces the target B.n to machine precision -- the secular term is now
    INTEGRATED into the design solve.

Marked slow: builds two torus surface meshes + an Euler-characteristic count.
"""
import os
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
EXDIR = os.path.join(HERE, "..", "stream_function")


@pytest.mark.slow
def test_cohomology_net_current_golden():
    sys.path.insert(0, EXDIR)
    try:
        import cohomology_net_current as M
    except Exception as e:  # pragma: no cover - env-dependent (gmsh/ngsolve)
        pytest.skip(f"cohomology net-current example not importable: {e}")

    res, ok = M.main()

    assert ok, f"example self-check failed: {res}"

    # --- generators are valid surface currents (div-free) ---
    for k, dv in res["gen_div_rel"].items():
        assert dv < 5e-2, f"generator {k} not div-free: div_s/|K|={dv:.2e}"

    # --- generators are GENUINE cohomology (non-exact: not any single-valued psi)
    for k, nx in res["gen_nonexactness"].items():
        assert nx > 0.8, f"generator {k} too exact (not cohomology): {nx:.3f}"

    # --- Gmsh confirms the secular-DOF count = b1 = 2 ---
    assert res["b1_euler"] == 2, f"b1 != 2: {res['b1_euler']}"
    assert res["n_secular_dofs"] == 2

    # --- TF generator obeys Ampere (B_tor R const inside, ~0 outside) ---
    assert res["tf_BR_spread_rel"] < 0.05, (
        f"TF field not 1/R: B_tor*R spread={res['tf_BR_spread_rel']:.3f}")
    assert res["tf_outside_over_inside"] < 0.05, (
        f"TF field leaks outside: {res['tf_outside_over_inside']:.2e}")

    # --- the joint design INTEGRATES the secular term: total reproduces target ---
    assert res["joint_total_residual"] < 1e-6, (
        f"joint (psi+secular) design residual too high: "
        f"{res['joint_total_residual']:.2e}")

    # --- prescribed net poloidal current for B_tor=1T (Ampere on the major axis)
    import math
    expect_Ipol = 2 * math.pi * 0.30 * 1.0 / (4.0e-7 * math.pi)
    assert abs(res["net_poloidal_current_A"] - expect_Ipol) / expect_Ipol < 1e-6
