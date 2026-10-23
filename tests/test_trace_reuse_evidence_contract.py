"""Bind trace reuse evidence to an executable runner and exercised cache."""
import hashlib
import importlib.util
import json
import subprocess
from pathlib import Path


def test_trace_reuse_runner_and_exercised_cache_are_recorded():
    folder = Path(__file__).parents[1] / "validation_test/mixed_omega_trace_reuse"
    record = json.loads((folder / "result_2607_20260930.json").read_text(encoding="utf-8"))
    runner = folder / "run.py"
    if record["runner_sha256"] != hashlib.sha256(runner.read_bytes()).hexdigest():
        root = Path(__file__).parents[1]
        spec = importlib.util.spec_from_file_location("privacy_source", root / "validation_test/radia_mcp/privacy_source_contract.py")
        contract = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(contract)
        relative = "validation_test/mixed_omega_trace_reuse/run.py"
        text = runner.read_text(encoding="utf-8")
        binding = record.get("source_privacy_redactions", {}).get(relative, {})
        if not contract.verifies_redaction(relative, record["runner_sha256"], text, binding):
            status = json.loads((folder / "source-status.json").read_text(encoding="utf-8"))
            assert status["status"] == "stale" and status["stale_reason"]
            assert status["runner"] == relative
            historical = subprocess.check_output(
                ["git", "show", status["historical_source_revision"] + ":" + relative], cwd=root)
            assert hashlib.sha256(historical).hexdigest() == binding["published_text_sha256"]
            text = historical.decode("utf-8")
        assert contract.verifies_redaction(relative, record["runner_sha256"], text, binding)
    assert record["ngsolve"] == "6.2.2607"
    assert len(record["cases"]) == 3
    for case in record["cases"]:
        assert case["cached_primal_factor"] is True
        assert len(case["updated"]) == len(case["baseline"]) == 3
        for row in case["updated"] + case["baseline"]:
            assert 0 <= row["residual"] < 1e-10
