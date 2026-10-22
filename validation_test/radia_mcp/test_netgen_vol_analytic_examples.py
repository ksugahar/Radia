"""Nightly analytic acceptance, distinct from the frozen-output regression lock."""
import json
import pytest
from validation_test.radia_mcp.netgen_vol_examples import CASES, example_module

pytestmark = pytest.mark.validation
VALIDATION_TIER = "nightly"

@pytest.mark.parametrize("module_name", CASES, ids=CASES)
def test_analytic_acceptance(module_name, tmp_path):
    with example_module(module_name) as module:
        result = module.validate_summary(module.build_summary())
    (tmp_path / (module_name + ".json")).write_text(json.dumps(result, sort_keys=True), encoding="utf-8")
