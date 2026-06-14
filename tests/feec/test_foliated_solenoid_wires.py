"""Golden: wire extraction from a foliated volume current (Stage B / B2).

Locks the Stage-B capstone -- a continuous div-free VOLUME current (the foliated
solenoid's solved lambda) is turned into WINDABLE DISCRETE EQUAL-CURRENT WIRES
whose Biot-Savart reproduces the target field:
  * helicity gate H_rel ~ 0 (extraction is Clebsch-feasible),
  * equal current per wire (I = Delta-lambda * Delta-mu, spread ~ 0),
  * field reproduction vs B0 (discretisation-limited),
  * TWO-CODEBASE agreement: the example's straight-segment Biot-Savart matches
    Radia's rad.ObjFlmCur + rad.Fld (so it cannot pass on a wrong formula).

See examples/vim/foliated_solenoid_wires.py and streamfunction("clebsch_3d").
"""
import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__),
                                "../../examples/vim"))


@pytest.fixture(scope="module")
def result():
    import foliated_solenoid_wires as fw
    res, ok = fw.main()
    return res, ok


def test_extraction_is_helicity_feasible(result):
    res, _ = result
    assert res["h_rel"] < 5e-2, f"gate-1 helicity {res['h_rel']:.2e}"


def test_wires_are_equal_current(result):
    """Every extracted wire carries I = Delta-lambda * Delta-mu (equal)."""
    res, _ = result
    assert res["equal_cur"] < 1e-9, f"equal-current spread {res['equal_cur']:.2e}"
    assert res["n_wires"] >= 9, f"too few wires: {res['n_wires']}"


def test_wires_reproduce_target_field(result):
    """Straight-segment Biot-Savart of the extracted wires ~ target Bz."""
    res, _ = result
    assert res["field_res"] < 8e-2, f"field residual {res['field_res']:.2e}"


def test_two_codebase_biot_savart_agreement(result):
    """The example's Biot-Savart matches Radia's rad.ObjFlmCur + rad.Fld."""
    res, _ = result
    assert res["twocode"] < 1e-3, f"two-codebase disagreement {res['twocode']:.2e}"


def test_overall_pass(result):
    _, ok = result
    assert ok
