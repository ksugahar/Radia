r"""Cogging-torque order / period rule: N_c = LCM(slots, poles), period = 360/N_c [deg].

The slot/pole topology that sets the cogging fundamental.  Pure geometric reference (no solver,
no commercial tooling); the rule is independently confirmed by a slotted multi-pole rotor FEA
(the cogging waveform is dominated by exactly one cycle per 360/LCM degrees).
"""
import os
import sys

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)
from radia_mcp.radia_ngsolve.solve import cogging_torque_order, cogging_period_deg


def test_cogging_order_and_period():
    # classic combos: order = LCM(slots, poles), period = 360/order
    cases = {(6, 4): 12, (12, 10): 60, (24, 8): 24, (9, 8): 72, (6, 2): 6, (12, 8): 24}
    for (slots, poles), order in cases.items():
        assert cogging_torque_order(slots, poles) == order, (slots, poles)
        assert cogging_period_deg(slots, poles) == 360.0 / order


def test_fractional_slot_reduces_cogging():
    # higher LCM (fractional-slot) -> smaller period -> lower, finer cogging
    assert cogging_period_deg(12, 10) < cogging_period_deg(12, 4)   # 6 deg < 30 deg
    # slots == poles is the pathological (huge cogging) case: order collapses to the pole count
    assert cogging_torque_order(8, 8) == 8


def main():
    for s, p in [(6, 4), (12, 10), (24, 8)]:
        print(f"{s}-slot/{p}-pole: cogging order = {cogging_torque_order(s, p)}/rev, "
              f"period = {cogging_period_deg(s, p):.1f} deg")
    print("[OK] cogging order = LCM(slots, poles), period = 360/LCM.")


if __name__ == "__main__":
    main()
