"""Golden lock: Stage-B line-topology gate (rotation-number closure defect).

Runs validation_test/feec/vim_legacy/line_topology_gate.py and asserts the gate's discriminator:

  - a RATIONAL-rotation-number torus winding (iota = 1/5, 2/5) CLOSES after q = 5
    toroidal turns (min closure-defect small) -> windable as a q-fold loop;
  - an IRRATIONAL winding (iota = 1/phi) never closes within Q_max turns (closure
    defect bounded away from zero) -> NOT windable as a finite loop family.

The robust closure metric (NOT a chaotic-Poincare 2-D-fill, which is unreliable on
the ABC mixed phase space) -- so the separation margin is wide and seed-stable.
"""
import os
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
EXDIR = os.path.join(HERE, "vim_legacy")


@pytest.mark.slow
def test_line_topology_gate_golden():
    sys.path.insert(0, EXDIR)
    try:
        import line_topology_gate as G
    except Exception as e:  # pragma: no cover - env-dependent (numpy only)
        pytest.skip(f"line topology gate example not importable: {e}")

    res, margin, ok = G.main()
    assert ok, f"example self-check failed: margin={margin}, res={res}"

    rat1 = res["iota = 1/5  (rational)"]
    rat2 = res["iota = 2/5  (rational)"]
    irr = res["iota = 1/phi (irrational)"]

    # rational windings close after q = 5 turns -> windable
    assert rat1["closed"] and rat1["closes_at_turn"] == 5, f"iota=1/5 did not close: {rat1}"
    assert rat2["closed"] and rat2["closes_at_turn"] == 5, f"iota=2/5 did not close: {rat2}"
    assert rat1["min_rel_defect"] < 0.15, f"rational min-defect too large: {rat1['min_rel_defect']:.3f}"

    # irrational winding never closes -> not windable; wide, robust margin
    assert not irr["closed"], f"iota=1/phi wrongly flagged closed: {irr}"
    assert irr["min_rel_defect"] > 0.3, f"irrational min-defect too small: {irr['min_rel_defect']:.3f}"
    assert margin > 3.0, f"closed/open separation margin too small: {margin:.1f}"
