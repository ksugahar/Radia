"""Decision contract of the IH fin physical-limit checker.

Guards: ``physical_verdict`` in validation_test/induction_heating/
fin_return_principle.py applies the user-owned 5 percent limit to equivalent
R, allocated partial L, root and tip heat, never adopts while the FEM
reference is unresolved, and fails loudly on malformed input.
Assumes: pure arithmetic only; no solve, no physics accuracy claim and no
surrogate baseline.  The limits themselves are user-owned.
"""

from __future__ import annotations

import importlib.util
import math
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "fin_return_principle", ROOT / "validation_test" / "induction_heating" / "fin_return_principle.py")
fin = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(fin)

REFERENCE = dict(R=2.0, partial_L=-3.0e-9, root_heat=4.0, tip_heat=5.0)


def scaled(factor, name=None):
    return {k: v * (factor if name in (None, k) else 1.0) for k, v in REFERENCE.items()}


def test_limits_are_the_user_owned_five_percent():
    assert fin.PHYSICAL_TARGETS == dict(R=.05, partial_L=.05, root_heat=.05, tip_heat=.05), fin.PHYSICAL_TARGETS


def test_within_limit_passes_only_with_resolved_reference():
    candidate = scaled(1.04)
    resolved = fin.physical_verdict(candidate, REFERENCE, True)
    assert resolved["status"] == "PASS", resolved
    unresolved = fin.physical_verdict(candidate, REFERENCE, False)
    assert unresolved["status"].startswith("NOT_ADOPTED"), unresolved
    assert unresolved["within_raw_limits"] is True, unresolved


@pytest.mark.parametrize("name", sorted(REFERENCE))
def test_any_quantity_above_limit_fails(name):
    result = fin.physical_verdict(scaled(1.06, name), REFERENCE, True)
    assert result["status"] == "FAIL", (name, result)
    assert math.isclose(result["differences"][name]["relative"], 0.06), result["differences"][name]


@pytest.mark.parametrize("flag", ["false", [False], 0, None, 1.0])
def test_reference_flag_must_be_bool(flag):
    with pytest.raises(TypeError, match="reference_resolved"):
        fin.physical_verdict(REFERENCE, REFERENCE, flag)


def test_zero_reference_is_unresolved():
    with pytest.raises(ValueError, match="zero reference"):
        fin.physical_verdict(REFERENCE, dict(REFERENCE, tip_heat=0.0), True)


@pytest.mark.parametrize("bad", [1 + 0.01j, True, "2.0"])
def test_non_real_values_are_rejected(bad):
    with pytest.raises(TypeError, match="real number"):
        fin.physical_verdict(dict(REFERENCE, R=bad), REFERENCE, True)


@pytest.mark.parametrize("bad", [math.nan, math.inf, -math.inf])
@pytest.mark.parametrize("side", ["candidate", "reference"])
def test_nonfinite_values_are_rejected(bad, side):
    values = dict(REFERENCE, root_heat=bad)
    args = (values, REFERENCE) if side == "candidate" else (REFERENCE, values)
    with pytest.raises(ValueError, match="finite"):
        fin.physical_verdict(*args, True)
