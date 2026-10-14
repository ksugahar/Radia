from radia_mcp.matlab import optuna_oracle
from pathlib import Path


def test_oracle_audit_rejects_a_stale_coverage_oracle_digest(monkeypatch):
    root=Path(__file__).resolve().parents[2]
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
