"""Replay the dynamic CoilBuilder/HCurl TEAM 28 Simulink evidence."""

from __future__ import annotations

import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parent.parent
RESULT = HERE / "team28_coilbuilder_dynamic_simulink_results.json"
MODEL = REPO_ROOT / "matlab" / "radia_team28_coilbuilder_dynamic.slx"


def test_dynamic_simulink_result_is_validated_and_replayable():
    payload = json.loads(RESULT.read_text(encoding="utf-8"))

    assert next(iter(payload)) == "radia_version"
    assert payload["schema"] == "cae-ai-lab.solver-run.v1"
    assert payload["pass"] is True
    assert payload["checks"]["validation_passed"] is True
    assert all(payload["checks"]["details"].values())
    assert payload["errors"]["terminal_height_abs_m"] < 0.1e-3
    assert payload["errors"]["terminal_force_balance_abs_N"] < 0.02
    assert payload["model_contract"]["damping_provenance"].startswith(
        "explicit control-oriented"
    )
    for relative_path in payload["result_files"]:
        assert (REPO_ROOT / relative_path).is_file()


def test_dynamic_model_is_nonempty_and_uses_the_validated_family():
    payload = json.loads(RESULT.read_text(encoding="utf-8"))

    assert MODEL.stat().st_size > 20_000
    assert payload["checks"]["details"][
        "common_rank_three_eddy_basis_is_retained"
    ]
    assert payload["observables"]["sample_count"] == 2001
    assert abs(
        payload["observables"]["terminal_upward_lift_N"]
        - payload["model_contract"]["disk_weight_N"]
    ) < 0.02
