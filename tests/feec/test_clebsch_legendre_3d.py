"""Golden lock: 3D Clebsch-Legendre (Variant 2) vacuum field-geometry solver.

Runs examples/clebsch_legendre/solve_clebsch_legendre_3d.py and asserts the
manufactured-solution verification:

  - the polynomial vacuum solution (uniform field, z=psi, phi=B0 y) is recovered
    to ~machine precision, and the reconstructed field matches B=(B0,0,0);
  - the non-polynomial vacuum solution (hyperbolic field B=(y,x,B0)) is recovered
    at the FE convergence rate (errz drops monotonically under h-refinement) and
    the reconstructed field matches -- i.e. minimising W gives the VACUUM field.

Marked slow: a few 3-D NGSolve nonlinear (Newton) energy minimisations.
"""
import os
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
EXDIR = os.path.join(HERE, "..", "..", "examples", "clebsch_legendre")


@pytest.mark.slow
def test_clebsch_legendre_3d_golden():
    sys.path.insert(0, EXDIR)
    try:
        import solve_clebsch_legendre_3d as M
    except Exception as e:  # pragma: no cover - env-dependent (ngsolve)
        pytest.skip(f"3D Clebsch-Legendre solver not importable: {e}")

    res, ok = M.main()
    assert ok, f"solver self-check failed: {res}"

    uni, hyp, conv = res["uniform"], res["hyperbolic"], res["convergence"]
    # polynomial vacuum solution -> machine precision, correct reconstructed field
    assert uni["errz"] < 1e-10 and uni["errp"] < 1e-10, f"uniform not exact: {uni}"
    assert uni["Berr"] < 1e-8, f"uniform field reconstruction wrong: {uni['Berr']:.2e}"
    # non-polynomial vacuum solution -> FE-convergent (proves it solves the PDE)
    assert hyp["errz"] < 1e-4, f"hyperbolic z error too large: {hyp['errz']:.2e}"
    assert hyp["Berr"] < 1e-3, f"hyperbolic field reconstruction: {hyp['Berr']:.2e}"
    errs = [c["errz"] for c in conv]
    assert errs[0] > errs[1] > errs[2], f"no h-convergence: {errs}"
