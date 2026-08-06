"""Replay the moving CoilBuilder/HCurl TEAM 28 family evidence."""

from __future__ import annotations

import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parent.parent
RESULT = HERE / "team28_coilbuilder_height_family_results.json"
FAMILY = HERE / "team28_coilbuilder_hcurl_eddy_cln_family.json"


def test_coilbuilder_height_family_is_validated_and_replayable():
    payload = json.loads(RESULT.read_text(encoding="utf-8"))

    assert next(iter(payload)) == "radia_version"
    assert payload["schema"] == "cae-ai-lab.solver-run.v1"
    assert payload["pass"] is True
    assert all(payload["checks"].values())
    details = payload["details"]
    assert details["problem"]["evrs_rank"] == 3
    assert len(details["height_family"]) == 25
    assert details["curve_comparison"]["normalized_max_absolute_error"] < 0.02
    assert details["curve_comparison"]["equilibrium_offset_error_mm"] < 0.5
    for relative_path in payload["result_files"]:
        assert (REPO_ROOT / relative_path).is_file()


def test_family_exchange_declares_common_coilbuilder_basis():
    payload = json.loads(FAMILY.read_text(encoding="utf-8"))

    assert payload["schema"] == "radia.hcurl.eddy_cln.family.v1"
    assert payload["shared_state_basis"] is True
    assert payload["snapshot_count"] == 25
    assert payload["state_order"] == 3
    assert payload["metadata"]["coil_source"] == "radia.coil_builder.CoilBuilder"
    heights = [snapshot["height_m"] for snapshot in payload["snapshots"]]
    assert heights == sorted(heights)
