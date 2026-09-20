"""Replay the independent three-dimensional ECB plate-force evidence."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
REPO = next(path for path in HERE.parents if (path / "src" / "radia").is_dir())
SUMMARY = HERE / "ecb_foster_lorentz_3d_reference_summary.json"
SCRIPT = HERE / "ecb_foster_lorentz_3d_reference.py"
LORENTZ = REPO / "src" / "radia" / "maglev" / "ecb" / "lorentz.py"


def _payload():
    return json.loads(SUMMARY.read_text(encoding="utf-8"))


def _sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def test_3d_reference_is_full_current_and_source_exact():
    payload = _payload()
    assert payload["schema"] == "radia.maglev.ecb-foster-lorentz-3d-reference.v1"
    assert payload["profile"] == "full"
    assert payload["pass"] is True
    assert all(payload["checks"].values()), payload["checks"]
    assert payload["source_sha256"] == {
        "validation_script": _sha256(SCRIPT),
        "lorentz_module": _sha256(LORENTZ),
    }
    for name in (
        "host",
        "python_version",
        "radia_version",
        "ngsolve_version",
        "numpy_version",
        "scipy_version",
    ):
        assert payload["runtime"][name]


def test_3d_rank_and_mesh_refinement_are_converged():
    payload = _payload()
    rows = payload["rows"]
    assert [(row["maxh_m"], row["steps"]) for row in rows] == [
        (0.012, 3),
        (0.012, 6),
        (0.008, 6),
    ]
    convergence = payload["convergence"]
    assert max(convergence["three_d_rank_lift_relative_differences"]) < 0.005
    assert max(convergence["three_d_mesh_lift_relative_differences"]) < 0.005
    assert all(row["passive"] for row in rows)


def test_converged_3d_lift_and_zero_torque_are_locked():
    payload = _payload()
    cases = {row["frequency_hz"]: row for row in payload["rows"][-1]["cases"]}
    expected_lift = {
        50.0: -0.30713,
        500.0: -6.6435,
        5000.0: -8.3740,
    }
    for frequency_hz, expected in expected_lift.items():
        case = cases[frequency_hz]
        assert case["force_N"][2] == pytest.approx(expected, rel=5.0e-3)
        assert case["horizontal_to_vertical_ratio"] < 0.01
        assert case["torque_to_force_length_ratio"] < 0.01
        assert -case["force_N"][2] < payload["problem"][
            "peak_phasor_image_lift_bound_N"
        ]


def test_scalar_model_is_converged_but_rejected_against_3d():
    payload = _payload()
    convergence = payload["convergence"]
    assert max(convergence["scalar_mesh_lift_relative_differences"]) < 0.001
    discrepancies = [
        row["lift_relative_difference"] for row in payload["comparisons"]
    ]
    assert discrepancies == pytest.approx([0.627, 0.866, 0.890], abs=0.01)
    assert min(discrepancies) > 0.5
    conclusion = payload["conclusion"]
    assert conclusion["scalar_model_status"] == (
        "rejected-for-quantitative-3d-force"
    )
    assert conclusion["recommended_backend"] == "hcurl_vim"
    assert conclusion["experiment_status"] == (
        "not-run-no-measurement-dataset-supplied"
    )
