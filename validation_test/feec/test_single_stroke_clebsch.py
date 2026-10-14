"""Golden lock: rigorous single-stroke via the Clebsch/cohomology-secular term.

Runs validation_test/feec/vim_legacy/single_stroke_clebsch.py and asserts:

  - the Clebsch single-stroke (the level-set {Psi = psi + secular} = a continuous
    helix) is a VALID coil: its on-axis Bz matches the ideal multi-turn rings;
  - near the connection, the Clebsch level-set (helix, which DISTRIBUTES the
    unavoidable net-axial 'lead' current = the secular term axisymmetrically)
    has clearly LESS stray than the HEURISTIC localised rung -- the point of
    putting the single-stroke on a rigorous cohomology footing.

Marked slow: a few Radia Biot-Savart field evaluations on N-turn coils.
"""
import os
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
EXDIR = os.path.join(HERE, "vim_legacy")


@pytest.mark.slow
def test_single_stroke_clebsch_golden():
    sys.path.insert(0, EXDIR)
    try:
        import single_stroke_clebsch as M
    except Exception as e:  # pragma: no cover - env-dependent (radia)
        pytest.skip(f"single-stroke example not importable: {e}")

    res, ok = M.main(n_turns=20)
    assert ok, f"example self-check failed: {res}"

    # the single-stroke helix is a valid solenoid (on-axis Bz ~ ideal rings)
    assert res["helix_bz_err"] < 0.06, (
        f"single-stroke helix Bz off the ideal design: {res['helix_bz_err']:.2e}")

    # the Clebsch level-set distributes the connection -> less near-join stray
    # than the heuristic rung (the demonstrable advantage of the rigorous route)
    assert res["rungs_over_helix"] > 1.5, (
        f"heuristic rung not worse near the join: "
        f"rungs/helix = {res['rungs_over_helix']:.2f}")

    # the connection current is the secular coefficient (net axial I)
    assert res["net_axial_current_A"] == pytest.approx(100.0)
