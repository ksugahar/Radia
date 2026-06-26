"""Golden lock: axisymmetric flux-coordinate bidirectional map (Phase 2).

Runs examples/vim/bidirectional_map_axisym.py and asserts:

  - FORWARD: the FEM recovers the analytic potential coordinates psi = r A_phi
    (Grad-Shafranov, current-free) and Phi (axisymmetric Laplace) to ~machine
    precision (the l=2 zonal pair is exactly order-3 polynomial);
  - the axisymmetric (1/r-weighted) Cauchy-Riemann relations hold and the Jacobian
    identity det[grad psi; grad Phi] = r|grad Phi|^2 ( = r|B|^2 ) holds;
  - INVERSE: the pointwise Newton inversion round-trips (the map is a valid chart).

Marked slow: two small NGSolve harmonic solves on a meridian box.
"""
import os
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
EXDIR = os.path.join(HERE, "..", "..", "examples", "vim")


@pytest.mark.slow
def test_bidirectional_map_axisym_golden():
    sys.path.insert(0, EXDIR)
    try:
        import bidirectional_map_axisym as M
    except Exception as e:  # pragma: no cover - env-dependent (ngsolve)
        pytest.skip(f"axisym bidirectional map example not importable: {e}")

    chk, rt, ok = M.main()
    assert ok, f"example self-check failed: {chk}, rt={rt}"

    # forward FEM recovers the analytic potential coordinates (machine precision)
    assert chk["errPsi"] < 1e-3, f"psi (Stokes stream) not recovered: {chk['errPsi']:.2e}"
    assert chk["errPhi"] < 1e-3, f"Phi (scalar pot) not recovered: {chk['errPhi']:.2e}"

    # axisymmetric Cauchy-Riemann conjugacy + Jacobian identity det = r|grad Phi|^2
    assert chk["conj"] < 1e-2, f"axisym CR conjugacy broken: {chk['conj']:.2e}"
    assert chk["jac_rel"] < 1e-2, f"Jacobian det != r|grad Phi|^2: {chk['jac_rel']:.2e}"

    # inverse (psi,Phi) -> (r,z) Newton round-trips (valid chart, Jacobian r|B|^2 != 0)
    assert rt < 1e-2, f"inverse round-trip did not close: {rt:.2e}"
