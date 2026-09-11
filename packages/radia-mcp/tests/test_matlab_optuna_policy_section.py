from radia_mcp.matlab.optuna_oracle import _policy_block
from radia_mcp.matlab import optuna_oracle
from pathlib import Path


def test_oracle_audit_rejects_a_stale_coverage_oracle_digest(monkeypatch):
    root=Path(__file__).resolve().parents[3]
    original=optuna_oracle._json
    def stale(path):
        data=original(path)
        if path.name == "optuna50_api_coverage.json":
            data=dict(data, upstream_oracle_sha256="0"*64)
        return data
    monkeypatch.setattr(optuna_oracle,"_json",stale)
    result=optuna_oracle.matlab_optuna_oracle_audit(root)
    assert result["status"] == "fail"
    assert "MATLAB API coverage was not generated from the checked oracle" in result["errors"]


def test_policy_comparison_excludes_neighboring_sections(tmp_path):
    agents = tmp_path / "AGENTS.md"
    claude = tmp_path / "CLAUDE.md"
    section = "## Optuna\n\nUse pinned upstream fixtures.\n\n### Details\nKeep parity.\n\n"
    agents.write_text("## Before\nAgent context.\n" + section + "## Git\nOne.\n", encoding="utf-8")
    claude.write_text(section + "## Git\nTwo.\n", encoding="utf-8")
    assert _policy_block(agents) == _policy_block(claude) == section.strip()


def test_policy_comparison_detects_real_drift(tmp_path):
    agents = tmp_path / "AGENTS.md"
    claude = tmp_path / "CLAUDE.md"
    agents.write_text("## Optuna\nUse upstream.\n", encoding="utf-8")
    claude.write_text("## Optuna\nUse handwritten values.\n", encoding="utf-8")
    assert _policy_block(agents) != _policy_block(claude)


def test_missing_policy_is_not_matched_by_similar_heading(tmp_path):
    path = tmp_path / "AGENTS.md"
    path.write_text("## Optuna History\nOld policy.\n", encoding="utf-8")
    assert _policy_block(path) == ""
