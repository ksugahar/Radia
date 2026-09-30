"""Bind trace reuse evidence to an executable runner and exercised cache."""
import hashlib
import json
from pathlib import Path


def test_trace_reuse_runner_and_exercised_cache_are_recorded():
    folder = Path(__file__).parents[1] / "validation_test/mixed_omega_trace_reuse"
    record = json.loads((folder / "result_2607_20260930.json").read_text(encoding="utf-8"))
    assert record["runner_sha256"] == hashlib.sha256((folder / "run.py").read_bytes()).hexdigest()
    assert record["ngsolve"] == "6.2.2607"
    assert len(record["cases"]) == 3
    for case in record["cases"]:
        assert case["cached_primal_factor"] is True
        assert len(case["updated"]) == len(case["baseline"]) == 3
        for row in case["updated"] + case["baseline"]:
            assert 0 <= row["residual"] < 1e-10
