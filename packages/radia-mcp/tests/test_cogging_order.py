r"""Cogging-torque order / period rule: N_c = LCM(slots, poles), period = 360/N_c [deg].

The slot/pole topology that sets the cogging fundamental.  Pure geometric reference (no solver,
no commercial tooling); the rule is independently confirmed by a slotted multi-pole rotor FEA
(the cogging waveform is dominated by exactly one cycle per 360/LCM degrees).
"""
import os
import sys

import pytest

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)
from radia_mcp.radia_ngsolve.solve import (
    cogging_period_deg,
    cogging_skew_plan,
    cogging_torque_order,
    machine_symmetry_sector,
)


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


def test_machine_symmetry_sector_is_gcd_pair_to_cogging_lcm():
    cases = {
        (36, 4): (4, 90.0, 9, 1, "anti-periodic"),
        (12, 10): (2, 180.0, 6, 5, "anti-periodic"),
        (24, 8): (8, 45.0, 3, 1, "anti-periodic"),
        (18, 6): (6, 60.0, 3, 1, "anti-periodic"),
        (24, 4): (4, 90.0, 6, 1, "anti-periodic"),
        (12, 12): (12, 30.0, 1, 1, "anti-periodic"),
    }
    for (slots, poles), (sectors, angle, qsec, psec, bc) in cases.items():
        out = machine_symmetry_sector(slots, poles)
        assert out["sectors"] == sectors
        assert out["symmetry_factor"] == sectors
        assert out["sector_angle_deg"] == angle
        assert out["slots_per_sector"] == qsec
        assert out["poles_per_sector"] == psec
        assert out["boundary"] == bc
        assert cogging_torque_order(slots, poles) == slots * poles // sectors

    assert machine_symmetry_sector(12, 8)["boundary"] == "periodic"  # two poles per 90-degree sector


def test_machine_symmetry_sector_validation():
    for bad in (
        lambda: machine_symmetry_sector(0, 4),
        lambda: machine_symmetry_sector(12, 0),
    ):
        try:
            bad()
            assert False, "invalid slot/pole count accepted"
        except ValueError:
            pass


def test_cogging_skew_plan_one_slot_cancels_cogging():
    plan = cogging_skew_plan(36, 4)
    assert plan["cogging_order_per_rev"] == 36
    assert plan["cogging_period_mech_deg"] == 10.0
    assert plan["skew_angle_mech_deg"] == 10.0
    assert plan["skew_angle_elec_deg"] == 20.0
    assert abs(plan["cogging_skew_factor"]) < 1.0e-12
    assert plan["fundamental_emf_skew_factor"] > 0.99
    assert plan["symmetry"]["sectors"] == 4


def test_cogging_skew_plan_fractional_slot_tradeoff():
    integer_slot = cogging_skew_plan(36, 4)
    fractional_slot = cogging_skew_plan(12, 10)
    half_slot = cogging_skew_plan(12, 10, skew_slot_pitches=0.5)

    assert fractional_slot["cogging_order_per_rev"] == 60
    assert fractional_slot["cogging_period_mech_deg"] == 6.0
    assert abs(fractional_slot["cogging_skew_factor"]) < 1.0e-12
    assert fractional_slot["fundamental_emf_skew_factor"] < integer_slot["fundamental_emf_skew_factor"]
    assert abs(half_slot["cogging_skew_factor"]) > abs(fractional_slot["cogging_skew_factor"])
    assert half_slot["fundamental_emf_skew_factor"] > fractional_slot["fundamental_emf_skew_factor"]


def test_continuous_loop_skew_tradeoff_gate_for_motor_teaching():
    integer_slot = cogging_skew_plan(36, 4, emf_harmonics=(1, 5, 7, 11, 13))
    fractional_slot = cogging_skew_plan(12, 10, emf_harmonics=(1, 5, 7, 11, 13))
    half_fractional = cogging_skew_plan(12, 10, skew_slot_pitches=0.5, emf_harmonics=(1, 5, 7, 11, 13))

    assert integer_slot["cogging_order_per_rev"] == 36
    assert integer_slot["cogging_period_mech_deg"] == 10.0
    assert integer_slot["skew_angle_mech_deg"] == 10.0
    assert abs(integer_slot["cogging_skew_factor"]) < 1.0e-12
    assert integer_slot["fundamental_emf_skew_factor"] == pytest.approx(0.9949307700452986)

    assert fractional_slot["cogging_order_per_rev"] == 60
    assert fractional_slot["skew_angle_mech_deg"] == pytest.approx(30.0)
    assert abs(fractional_slot["cogging_skew_factor"]) < 1.0e-12
    assert fractional_slot["fundamental_emf_skew_factor"] == pytest.approx(0.7379129755873375)

    assert half_fractional["skew_slot_pitches"] == 0.5
    assert half_fractional["fundamental_emf_skew_factor"] == pytest.approx(0.9301189496680686)
    assert abs(half_fractional["cogging_skew_factor"]) == pytest.approx(0.1273239544735163)


def main():
    for s, p in [(6, 4), (12, 10), (24, 8)]:
        print(f"{s}-slot/{p}-pole: cogging order = {cogging_torque_order(s, p)}/rev, "
              f"period = {cogging_period_deg(s, p):.1f} deg")
    print("[OK] cogging order = LCM(slots, poles), period = 360/LCM.")


if __name__ == "__main__":
    main()
