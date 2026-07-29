from __future__ import annotations

import hashlib
import json
from pathlib import Path


ROOT = (
    Path(__file__).resolve().parents[3]
    / "validation_test"
    / "radia_mcp"
    / "artifacts"
    / "field_study_production_v1"
)
FAMILIES = {"P1", "Q1", "P2", "Q2", "P2_curved", "Q2_curved"}
PHYSICS_CASES = {"steady_heat", "current_flow_dc", "current_flow_ac"}


def test_manifest_freezes_all_six_verified_element_families() -> None:
    manifest_bytes = (ROOT / "manifest.json").read_bytes()
    assert b"\r\n" not in manifest_bytes
    manifest = json.loads(manifest_bytes)
    assert manifest["schema"] == "radia.field-study-production-manifest.v1"
    assert manifest["all_passed"] is True
    assert manifest["execution_version"]["radia_mcp"]
    assert manifest["execution_version"]["ngsolve"]
    rows = {row["element_family"]: row for row in manifest["element_families"]}
    assert set(rows) == FAMILIES
    for family, row in rows.items():
        path = ROOT / row["artifact"]
        artifact_bytes = path.read_bytes()
        assert b"\r\n" not in artifact_bytes
        assert hashlib.sha256(artifact_bytes).hexdigest() == row["sha256"]
        artifact = json.loads(path.read_text(encoding="utf-8"))
        assert artifact["schema"] == "cae-ai-lab.solver-run.v1"
        assert artifact["pass"] is True
        assert artifact["checks"]["validation_passed"] is True
        assert artifact["production_contract"]["element_family"] == family
        assert artifact["production_contract"]["generated_vol_git_required"] is False
        assert artifact["tool_versions"]["radia_mcp"]
        assert artifact["tool_versions"]["ngsolve"]
        assert artifact["timing_breakdown_s"]


def test_manifest_freezes_steady_heat_and_dc_ac_current_flow() -> None:
    manifest = json.loads((ROOT / "manifest.json").read_text(encoding="utf-8"))
    rows = {row["physics_case"]: row for row in manifest["physics_cases"]}
    assert set(rows) == PHYSICS_CASES
    for case_id, row in rows.items():
        path = ROOT / row["artifact"]
        assert hashlib.sha256(path.read_bytes()).hexdigest() == row["sha256"]
        artifact = json.loads(path.read_text(encoding="utf-8"))
        assert artifact["schema"] == "cae-ai-lab.solver-run.v1"
        assert artifact["pass"] is True
        assert artifact["checks"]["validation_passed"] is True
        assert artifact["checks"]["replay_gate_accepted"] is True
        assert artifact["production_contract"]["physics_case"] == case_id
        assert artifact["production_contract"]["generated_vol_git_required"] is False
        assert artifact["tool_versions"]["radia_mcp"]
        assert artifact["tool_versions"]["ngsolve"]
        assert artifact["timing_breakdown_s"]
