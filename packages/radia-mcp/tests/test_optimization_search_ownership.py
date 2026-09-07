"""The ownership migration must not fork algorithms or break historical imports."""

import numpy as np
import pytest

from radia_mcp.optimization import global_optimizers as common
from radia_mcp.topology_optimization import global_optimizers as legacy


@pytest.mark.parametrize("name", [
    "constraint_violation", "best_feasible_record", "differential_evolution",
])
def test_legacy_import_is_canonical(name):
    assert getattr(legacy, name) is getattr(common, name)


def test_seeded_search_preserves_return_contract():
    objective = lambda x: float(np.dot(x, x))
    kwargs = dict(bounds=[(-1, 1), (-1, 1)], popsize=3, maxiter=3, seed=12)
    old = legacy.differential_evolution(objective, **kwargs)
    new = common.differential_evolution(objective, **kwargs)
    assert set(new) == {"x", "fun", "nfev", "nit"}
    np.testing.assert_array_equal(old["x"], new["x"])
    assert old["fun"] == new["fun"]
    assert new["nfev"] == 24


def test_selection_returns_copy_with_explicit_infeasible_fallback():
    records = [{"value": 0, "constraints": [2]}, {"value": 9, "constraints": [1]}]
    chosen = common.best_feasible_record(records)
    assert chosen["value"] == 9
    assert chosen["feasible"] is False
    assert "feasible" not in records[1]


def test_guide_records_completed_ownership_and_regression_gate_limits():
    from radia_mcp.optimization.server import optimization_guide
    guide = optimization_guide()
    assert "radia_mcp.optimization.global_optimizers" in guide["search_ownership"]["helpers"]
    assert "false convergence" in guide["search_ownership"]["regression_gates"]
    assert "global-search helper ownership audit" not in guide["deferred"]
