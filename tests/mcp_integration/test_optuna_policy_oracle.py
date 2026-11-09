from radia_mcp.matlab import optuna_oracle


def test_oracle_audit_rejects_a_stale_coverage_oracle_digest(monkeypatch):
    original=optuna_oracle._json
    def stale(path):
        data=original(path)
        if path.name == "optuna50_api_coverage.json":
            data=dict(data, upstream_oracle_sha256="0"*64)
        return data
    monkeypatch.setattr(optuna_oracle,"_json",stale)
    result=optuna_oracle.matlab_optuna_oracle_audit()
    assert result["status"] == "fail"
    assert "MATLAB API coverage was not generated from the checked oracle" in result["errors"]
