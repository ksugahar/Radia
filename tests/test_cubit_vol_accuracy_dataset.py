"""Fast contract checks for the measured complex-shape .vol dataset."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "validation_test" / "cubit" / "vol_accuracy_dataset"
CASES_PATH = DATA_DIR / "cases.json"
DATASET_PATH = DATA_DIR / "volume_accuracy_dataset.json"
JSONL_PATH = DATA_DIR / "volume_accuracy_rows.jsonl"


def _load():
    cases = json.loads(CASES_PATH.read_text(encoding="utf-8-sig"))
    dataset = json.loads(DATASET_PATH.read_text(encoding="utf-8-sig"))
    jsonl_rows = [
        json.loads(line) for line in JSONL_PATH.read_text(
            encoding="utf-8-sig"
        ).splitlines() if line.strip()
    ]
    return cases, dataset, jsonl_rows


def _command_hash(case):
    normalized = "\n".join(value.strip() for value in case["commands"])
    return hashlib.sha256(normalized.encode("utf-8")).hexdigest()


def test_complex_shape_dataset_has_reproducible_measured_rows():
    cases, dataset, jsonl_rows = _load()

    assert cases["schema"] == "radia.cubit-vol-accuracy-cases.v1"
    assert dataset["schema"] == "radia.cubit-vol-accuracy-dataset.v1"
    assert dataset["orders"] == [1, 2, 3]
    assert dataset["rows"] == jsonl_rows

    case_by_id = {case["case_id"]: case for case in cases["cases"]}
    assert len(case_by_id) >= 6
    acceptance = [
        case for case in cases["cases"] if case["role"] == "acceptance"
    ]
    assert len(acceptance) >= 5
    assert all(len(case["features"]) >= 3 for case in cases["cases"])
    assert {
        "multi_tool_boolean_cut",
        "mixed_open_closed_boolean_cut",
        "union_and_multi_tool_boolean_cut",
        "trimmed_loft",
        "trimmed_sweep",
        "boolean_union_boundary_layer",
    }.issubset({case["geometry_family"] for case in acceptance})

    rows = dataset["rows"]
    assert len(rows) == len(case_by_id) * 3
    assert len({row["sample_id"] for row in rows}) == len(rows)
    for case_id, case in case_by_id.items():
        case_rows = [row for row in rows if row["case_id"] == case_id]
        assert {row["curve_order"] for row in case_rows} == {1, 2, 3}
        assert all(row["role"] == case["role"] for row in case_rows)
        assert all(row["case_command_sha256"] == _command_hash(case)
                   for row in case_rows)

    for row in rows:
        assert row["cubit_cad_volume"] > 0.0
        assert row["ngsolve_volume"] > 0.0
        expected = (
            (row["ngsolve_volume"] - row["cubit_cad_volume"])
            / row["cubit_cad_volume"] * 100.0
        )
        assert math.isfinite(expected)
        assert row["signed_error_pct"] == pytest.approx(expected, abs=1e-13)
        assert row["absolute_error_pct"] == pytest.approx(abs(expected), abs=1e-13)
        assert row["cad_metadata_passed"] is True
        assert row["mapping_sample_count"] >= row["volume_element_count"]
        if row["curve_order"] >= 2:
            assert row["mapping_sample_count"] > row["volume_element_count"]


def test_all_acceptance_complex_shapes_are_valid_and_high_order_accurate():
    _, dataset, _ = _load()
    rows = [row for row in dataset["rows"] if row["role"] == "acceptance"]

    assert len(rows) >= 15
    assert dataset["summary"]["all_acceptance_structural_gate_passed"] is True
    assert dataset["summary"][
        "all_acceptance_order2plus_solver_ready_at_1pct"
    ] is True
    assert all(row["structural_gate_passed"] for row in rows)
    assert all(row["quality_passed"] for row in rows)
    assert all(row["invalid_jacobian_sample_count"] == 0 for row in rows)
    assert all(row["orientation_flip_sample_count"] == 0 for row in rows)
    # Cubit tetrahedra use a consistently negative reference orientation;
    # only a sign change within one element is a folded-map failure.
    assert all(
        0 <= row["negative_jacobian_sample_count"]
        <= row["mapping_sample_count"]
        for row in rows
    )
    assert all(row["structural_gate_issues"] == [] for row in rows)
    assert all(
        row["solver_ready_at_accuracy_threshold"]
        for row in rows if row["curve_order"] >= 2
    )

    by_case = {}
    for row in rows:
        by_case.setdefault(row["case_id"], {})[row["curve_order"]] = row
    for order_rows in by_case.values():
        assert order_rows[2]["absolute_error_pct"] < order_rows[1][
            "absolute_error_pct"
        ]
        assert order_rows[3]["absolute_error_pct"] < order_rows[2][
            "absolute_error_pct"
        ]

    challenging_cases = {
        "loft_circle_to_rectangle",
        "gapped_torus_order3_refined",
    }
    strict_rows = [
        row for row in rows if row["case_id"] not in challenging_cases
    ]
    order2 = [
        row["absolute_error_pct"] for row in strict_rows
        if row["curve_order"] == 2
    ]
    order3 = [
        row["absolute_error_pct"] for row in strict_rows
        if row["curve_order"] == 3
    ]
    assert max(order2) < 0.05
    assert max(order3) < 0.003


def test_periodic_seam_loft_is_valid_not_merely_volume_accurate():
    _, dataset, _ = _load()
    loft = [
        row for row in dataset["rows"]
        if row["case_id"] == "loft_circle_to_rectangle"
        and row["curve_order"] >= 2
    ]
    assert {row["curve_order"] for row in loft} == {2, 3}
    assert all(row["role"] == "acceptance" for row in loft)
    assert all(row["structural_gate_passed"] for row in loft)
    assert all(row["invalid_jacobian_sample_count"] == 0 for row in loft)
    assert all(row["orientation_flip_sample_count"] == 0 for row in loft)
    assert all(row["negative_jacobian_sample_count"] == 0 for row in loft)
    assert all(row["minimum_scaled_jacobian"] > 0.5 for row in loft)
    assert all(row["absolute_error_pct"] < 0.05 for row in loft)


def test_gapped_torus_records_the_order3_resolution_boundary():
    _, dataset, _ = _load()
    rows = dataset["rows"]
    coarse = {
        row["curve_order"]: row for row in rows
        if row["case_id"] == "gapped_torus_order3_coarse"
    }
    refined = {
        row["curve_order"]: row for row in rows
        if row["case_id"] == "gapped_torus_order3_refined"
    }
    assert set(coarse) == {1, 2, 3}
    assert set(refined) == {1, 2, 3}

    # At the production order-2 resolution, order 3 looks accurate by volume
    # but is not solver-ready because the cubic map changes orientation.
    assert coarse[3]["absolute_error_pct"] < 0.05
    assert coarse[3]["invalid_jacobian_sample_count"] > 0
    assert coarse[3]["orientation_flip_sample_count"] > 0
    assert coarse[3]["structural_gate_passed"] is False

    # One refinement step removes the fold at both supported high orders.
    for order in (2, 3):
        assert refined[order]["structural_gate_passed"] is True
        assert refined[order]["invalid_jacobian_sample_count"] == 0
        assert refined[order]["orientation_flip_sample_count"] == 0
    assert refined[2]["absolute_error_pct"] < 0.11
    assert refined[3]["absolute_error_pct"] < 0.02
    assert refined[3]["minimum_scaled_jacobian"] > 0.01
    assert dataset["summary"]["diagnostic_structural_failure_count"] >= 1
