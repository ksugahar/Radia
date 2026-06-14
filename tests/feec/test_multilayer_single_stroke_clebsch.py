"""Golden lock: multi-layer single-stroke via the foliated Clebsch (lambda,mu) path.

Runs examples/vim/multilayer_single_stroke_clebsch.py and asserts:

  - the boustrophedon single-stroke (one wire threading M layers) is a VALID
    M-layer solenoid: on-axis Bz ~ M x the single layer (the azimuthal sense is
    the same in every layer, so Bz adds);
  - the ALTERNATING axial direction makes the per-layer axial-lead (mu-step
    secular) currents CANCEL -> the net-axial far-field stray is much smaller
    than naively stacking M same-handed helices -- the rigorous secular
    bookkeeping of the (lambda, mu) path.

Marked slow: a few Radia Biot-Savart evaluations on M-layer helical coils.
"""
import os
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
EXDIR = os.path.join(HERE, "..", "..", "examples", "vim")


@pytest.mark.slow
def test_multilayer_single_stroke_clebsch_golden():
    sys.path.insert(0, EXDIR)
    try:
        import multilayer_single_stroke_clebsch as M
    except Exception as e:  # pragma: no cover - env-dependent (radia)
        pytest.skip(f"multilayer single-stroke example not importable: {e}")

    res, ok = M.main()
    assert ok, f"example self-check failed: {res}"

    # valid M-layer solenoid: Bz ~ M x single layer (azimuthal sense adds)
    assert abs(res["bz_ratio_vs_Mxsingle"] - 1.0) < 0.15, (
        f"boustrophedon not a valid M-layer solenoid: "
        f"Bz/(M*single) = {res['bz_ratio_vs_Mxsingle']:.3f}")

    # alternating mu-steps cancel the per-layer axial-lead secular currents
    assert res["secular_cancellation"] > 5.0, (
        f"axial-lead secular not cancelled vs same-handed stack: "
        f"cancellation = {res['secular_cancellation']:.1f}x")
