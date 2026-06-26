"""Golden lock: reduced-potential Clebsch hodograph (path A) on the exact sphere.

Runs examples/clebsch_legendre/reduced_clebsch_hodograph.py and asserts:

  - the hodograph recovers the flux psi and potential Phi of the real
    reduced-potential field (sphere = uniform source + induced-dipole reaction);
  - the field reconstructed from the FLUX coordinate and from the POTENTIAL
    coordinate AGREE (bidirectional consistency) and match the EXACT sphere field;
  - the chart Jacobian det[grad psi; grad Phi] = r|grad Phi|^2 is non-degenerate.

Marked slow: two small axisymmetric NGSolve solves on a meridian box.
"""
import os
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
EXDIR = os.path.join(HERE, "..", "..", "examples", "clebsch_legendre")


@pytest.mark.slow
def test_reduced_clebsch_hodograph_golden():
    sys.path.insert(0, EXDIR)
    try:
        import reduced_clebsch_hodograph as M
    except Exception as e:  # pragma: no cover - env-dependent (ngsolve)
        pytest.skip(f"reduced-potential Clebsch hodograph not importable: {e}")

    res, ok = M.main()
    assert ok, f"solver self-check failed: {res}"

    # the hodograph recovers both potential coordinates of the reduced field
    assert res["errPsi"] < 1e-5, f"flux psi not recovered: {res['errPsi']:.2e}"
    assert res["errPhi"] < 1e-5, f"potential Phi not recovered: {res['errPhi']:.2e}"
    # field from flux and from potential AGREE (bidirectional consistency)
    assert res["agree"] < 1e-4, f"flux/potential field disagree: {res['agree']:.2e}"
    # and match the EXACT sphere field
    assert res["acc"] < 1e-4, f"field off the exact sphere solution: {res['acc']:.2e}"
    # valid chart (non-degenerate Jacobian)
    assert res["jac"] < 1e-3, f"chart Jacobian residual too large: {res['jac']:.2e}"
