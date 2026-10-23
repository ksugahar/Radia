"""Package-owned CUDA oracle identity contract.

Guards: Changes to the approved producer or stored CUDA expectations, drift of
its upstream pin, and reintroducing monorepo-relative verification inputs.
Assumes: Stored upstream Optuna 5.0.0 evidence; producer identity uses its
existing LF normalization. This checks source identity, not solver accuracy.
"""

import hashlib
import json
from pathlib import Path
import tomllib


def test_cuda_oracle_preserves_approved_source_identity():
    package = Path(__file__).resolve().parents[1]
    inputs = package / "oracles"
    producer = (inputs / "generate_optuna50_cuda_oracle.py").read_bytes()
    producer_hash = hashlib.sha256(producer.replace(b"\r\n", b"\n")).hexdigest()
    expected_producer = "68b26d2ad46e9c7595f6f9ead22894ab3df3fb9c77ee6f96caec3acca993957b"
    assert producer_hash == expected_producer, (producer_hash, expected_producer)
    fixture = (inputs / "optuna50_cuda_oracle.json").read_bytes()
    fixture_hash = hashlib.sha256(fixture).hexdigest()
    expected_fixture = "e72ccf5cdb84cbd50c61d8b445d1310e02ea1bf9a2eab51baddc2d7eec29e9bf"
    assert fixture_hash == expected_fixture, (fixture_hash, expected_fixture)
    oracle = json.loads(fixture)
    assert oracle["generator_sha256"] == producer_hash, oracle["generator_sha256"]
    assert oracle["optuna_version"] == "5.0.0", oracle["optuna_version"]
    project = tomllib.loads((package / "pyproject.toml").read_text(encoding="utf-8"))
    upstream = project["project"]["optional-dependencies"]["upstream"]
    assert f"optuna=={oracle['optuna_version']}" in upstream, upstream
