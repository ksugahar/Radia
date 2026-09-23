"""Integrity and scope of retained spatial evidence; no numerical solve."""

import hashlib
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent


def read(name):
    return json.loads((ROOT / name).read_bytes())


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_retained_inputs_and_results_match_the_record():
    record = read("compute_record.json")
    assert record["hostname"].lower() == "mdx2"
    for name, expected in record["inputs"].items():
        assert digest(ROOT / "inputs" / name) == expected
    assert [case["frequency_hz"] for case in record["cases"]] == [10000, 50000, 100000]
    for case in record["cases"]:
        assert case["returncode"] == 0
        stem = f"f{case['frequency_hz']}"
        assert digest(ROOT / f"{stem}.json") == case["result_sha256"]
        assert digest(ROOT / f"{stem}_qsurf.sol") == case["qsurf_sol_sha256"]
        result = read(f"{stem}.json")
        assert result["esim_converged"]
        assert result["current_A"] == 100
        assert result["coupling_mode"] == "weak"
        assert result["wp_basis_order"] == 1
        assert len(result["esim_per_panel_H_t"]) == result["wp_ndof"]
        assert len(result["esim_per_panel_Z_s_real"]) == result["wp_ndof"]
        assert min(result["esim_per_panel_Z_s_real"]) > 0


def test_mesh_tolerance_is_scoped_to_the_boundary_problem():
    check = read("vol_check.json")
    assert check["passed"] and check["threshold_pct"] == 10
    for family, name in (("materials", "workpiece"), ("boundaries", "sibc")):
        row = next(row for row in check[family] if row["name"] == name)
        assert abs(row["error_pct"]) < 1


def test_transfer_conservation_does_not_make_local_loss_identical():
    metrics = read("spatial_metrics.json")
    assert len(metrics) == 3
    for result in metrics:
        assert min(result["P_BIE_W"], result["P_local_W"], result["P_transfer_W"]) > 0
        assert result["transfer_relative_gap"] < 1e-5
    middle = next(row for row in metrics if row["frequency_hz"] == 50000)
    assert middle["P_local_W"] != pytest.approx(middle["P_transfer_W"], rel=0.01)
