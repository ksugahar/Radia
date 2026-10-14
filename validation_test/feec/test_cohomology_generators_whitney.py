"""Golden lock: cohomology generators computed FROM SCRATCH (matched Whitney complex).

Runs validation_test/stream_function/cohomology_generators_whitney.py (torus) and asserts
the harmonic 1-form generators come out of the EXACT discrete de Rham complex
(H1 order 1 -> HCurl order 0) cleanly:

  - the complex is exact (curl of a gradient ~ 0 to machine precision);
  - exactly b1 = 2 near-zero Hodge-Laplacian eigenvalues with a HUGE gap;
  - each generator is div-free AND curl-free (a genuine harmonic form);
  - the generators reproduce the analytic generators' NET-CURRENT (TF Ampere)
    field -- the SAME cohomology class -- with no analytic ansatz used to build
    them;
  - the surface Euler characteristic (gmsh-free) independently confirms b1 = 2.

Marked slow: a torus surface mesh + a small sparse eigensolve + an Euler-char count.
"""
import os
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
EXDIR = os.path.join(HERE, "..", "stream_function")


@pytest.mark.slow
def test_cohomology_generators_whitney_golden():
    sys.path.insert(0, EXDIR)
    try:
        import cohomology_generators_whitney as M
    except Exception as e:  # pragma: no cover - env-dependent (ngsolve)
        pytest.skip(f"Whitney generators example not importable: {e}")

    res, ok = M.main()
    assert ok, f"example self-check failed: {res}"

    # --- the discrete de Rham complex is EXACT (curl of a gradient ~ 0) ---
    assert res["exactness_resid"] < 1e-8, (
        f"complex not exact: curl(grad) energy = {res['exactness_resid']:.2e} "
        "(mismatched HCurl/H1 order?)")

    # --- exactly b1 = 2 harmonic forms, cleanly isolated ---
    assert res["b1_numeric"] == 2, f"numeric b1 != 2: {res['b1_numeric']}"
    assert res["eig_gap"] > 1e6, (
        f"harmonic forms not isolated: eig gap = {res['eig_gap']:.2e}")

    # --- each generator is a genuine harmonic form (div-free AND curl-free) ---
    for j, (ce, de) in enumerate(zip(res["gen_curl_energy"],
                                     res["gen_div_energy"])):
        assert abs(ce) < 1e-6, f"generator {j} not curl-free: {ce:.2e}"
        assert abs(de) < 1e-6, f"generator {j} not div-free: {de:.2e}"

    # --- the numeric generators carry the right NET CURRENT (same class) ---
    assert res["class_match_resid"] < 0.05, (
        f"numeric generators wrong cohomology class: TF net-current fit "
        f"residual = {res['class_match_resid']:.2e}")

    # --- the Euler characteristic (gmsh-free) confirms the generator count ---
    assert res["b1_euler"] == 2, f"Euler-char b1 != 2: {res['b1_euler']}"
