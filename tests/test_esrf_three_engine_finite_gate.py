"""Execute production gate functions via AST; no native imports or solves."""
from __future__ import annotations

import ast
import itertools
import json
from pathlib import Path

import numpy as np
import pytest

RUNNER = Path(__file__).resolve().parents[1] / "validation_test/esrf_three_engine/run_coil_yoke_three_engine.py"
NAMES = ("hdiv_mmm", "reduced_a", "mixed_total_reduced_omega")


@pytest.fixture
def gate():
    tree = ast.parse(RUNNER.read_text(encoding="utf-8"))
    functions = {"_validated_field", "_validated_tolerance", "_relative_rms", "_pairwise",
                 "_comparison_gate", "_is_converged_result", "_legacy_contract",
                 "_read_checkpoint", "_write_checkpoint"}
    ns = {"np": np, "Path": Path, "json": json,
          "CHECKPOINT_SCHEMA": "radia.validation.esrf-coil-yoke-checkpoint.v3",
          "LEGACY_CHECKPOINT_SCHEMAS": ("radia.validation.esrf-coil-yoke-checkpoint.v2",),
          "LEGACY_CAP_KEY": "nonlinear_maximum_iterations"}
    selected = [n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name in functions]
    exec(compile(ast.Module(body=selected, type_ignores=[]), str(RUNNER), "exec"), ns)  # noqa: S102 - reviewed local AST only
    return ns


@pytest.mark.parametrize("order", list(itertools.permutations(NAMES)))
@pytest.mark.parametrize("bad_engine", NAMES)
@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf])
def test_nonfinite_any_engine_any_order_rejected(gate, order, bad_engine, bad):
    fields = {name: np.ones((2, 3)) for name in order}
    fields[bad_engine][1, 2] = bad
    with pytest.raises(ValueError, match="finite"):
        gate["_comparison_gate"](fields, np.ones(2, dtype=bool), .03)


@pytest.mark.parametrize("bad", [[], [[1, 2]], [[1, 2, 3]], [[True]*3]*2,
                                   [[1j]*3]*2, [["1"]*3]*2])
def test_shape_type_and_empty(gate, bad):
    fields = {name: np.ones((2, 3)) for name in NAMES}
    fields[NAMES[-1]] = bad
    with pytest.raises(ValueError):
        gate["_comparison_gate"](fields, np.ones(2, dtype=bool), .03)


@pytest.mark.parametrize("names", [(), NAMES[:2], (*NAMES, "extra")])
def test_required_engine_set(gate, names):
    with pytest.raises(ValueError, match="three named engines"):
        gate["_comparison_gate"]({n: np.ones((2, 3)) for n in names}, np.ones(2, dtype=bool), .03)


@pytest.mark.parametrize("tol", [0, -1, 1, 2, np.nan, np.inf, -np.inf, True, ".03"])
def test_invalid_tolerance(gate, tol):
    with pytest.raises(ValueError, match="tolerance"):
        gate["_comparison_gate"]({n: np.ones((2, 3)) for n in NAMES}, np.ones(2, dtype=bool), tol)


@pytest.mark.parametrize("mask", [[], [0, 1], [False, False], [True], [[True, True]]])
def test_invalid_selector(gate, mask):
    with pytest.raises(ValueError, match="selector"):
        gate["_comparison_gate"]({n: np.ones((2, 3)) for n in NAMES}, mask, .03)


def test_valid_agreement_and_mismatch(gate):
    fields = {n: np.ones((2, 3)) for n in NAMES}
    assert gate["_comparison_gate"](fields, np.ones(2, dtype=bool), .03)[1:] == (0, True)
    fields[NAMES[-1]] *= 1.1
    pairs, maximum, passed = gate["_comparison_gate"](fields, np.ones(2, dtype=bool), .03)
    assert len(pairs) == 3 and maximum == pytest.approx(.1) and not passed


@pytest.mark.parametrize("value", [0.0, 1e308])
def test_zero_or_overflow_norm(gate, value):
    with np.errstate(over="ignore", invalid="ignore"), pytest.raises(ValueError):
        gate["_comparison_gate"]({n: np.full((2, 3), value) for n in NAMES},
                                  np.ones(2, dtype=bool), .03)


@pytest.mark.parametrize("bad", [np.full((2, 3), np.nan), np.ones((1, 3)), np.ones((2, 2))])
def test_checkpoint_read_write_reject_bad_field_without_overwriting(gate, bad, tmp_path):
    path = tmp_path / "checkpoint.json"
    contract = {"observation_points_m": [[0, 0, 0], [1, 0, 0]]}
    diag = {"nonlinear": True, "nonlinear_stats": {"converged": True}}
    payload = {"schema": gate["CHECKPOINT_SCHEMA"], "contract": contract,
               "diagnostics": diag, "field_T": bad.tolist()}
    path.write_text(json.dumps(payload), encoding="utf-8")
    original = path.read_bytes()
    with pytest.raises(ValueError):
        gate["_read_checkpoint"](path, contract)
    with pytest.raises(ValueError):
        gate["_write_checkpoint"](path, contract, bad, diag, {})
    assert path.read_bytes() == original


@pytest.mark.parametrize("value", ["false", 1, None, False])
def test_truthy_convergence_is_not_true(gate, value):
    assert not gate["_is_converged_result"]({"nonlinear_stats": {"converged": value}})


def test_checkpoint_roundtrip(gate, tmp_path):
    path = tmp_path / "checkpoint.json"
    contract = {"observation_points_m": [[0, 0, 0], [1, 0, 0]]}
    diag = {"nonlinear": True, "nonlinear_stats": {"converged": True}}
    gate["_write_checkpoint"](path, contract, np.ones((2, 3)), diag, {})
    np.testing.assert_array_equal(gate["_read_checkpoint"](path, contract)[0], np.ones((2, 3)))
