"""Package-owned upstream oracle identity contract.

Guards: Changes to the approved producer or stored oracle expectations, drift of
its upstream pin, and reintroducing monorepo-relative verification inputs.
Assumes: Stored upstream Optuna 5.0.0 evidence; producer identity uses its
existing LF normalization. This checks source identity, not solver accuracy.
"""

import hashlib
import json
from pathlib import Path
import tomllib

import pytest


@pytest.mark.parametrize(
    "producer_name,fixture_name,expected_producer,expected_fixture,embedded_hash",
    [
        pytest.param(
            "generate_optuna50_cuda_oracle.py", "optuna50_cuda_oracle.json",
            "68b26d2ad46e9c7595f6f9ead22894ab3df3fb9c77ee6f96caec3acca993957b",
            "e72ccf5cdb84cbd50c61d8b445d1310e02ea1bf9a2eab51baddc2d7eec29e9bf",
            True, id="cuda50",
        ),
        pytest.param(
            "generate_optuna_cuda_density_oracle.py", "optuna_cuda_density_oracle.json",
            "434e0ed2aacb6dec0db5d9e99a99f9a05cf0d75fb02d04e2e9310fc94bfe654f",
            "04d11672c5dfd73f7dc399fc237d2f8e1c8b01c2bff06cb1a7ed5085f807dbfc",
            False, id="density",
        ),
        pytest.param(
            "generate_optuna50_mcp_oracle.py", "optuna50_mcp_oracle.json",
            "7d6a08a91e8fe21a438660e6f359bcfa468578efd9c7579c6e592dbec56dd6cb",
            "890bd0c25e6842427db8c55a6c1cda4135afdacb03ea2d53b5c4785aa15aa8aa",
            False, id="mcp50",
        ),
        pytest.param(
            "generate_optuna50_api_inventory.py", "optuna50_public_api.json",
            "351d63e534dcca09fc9b47493e4942cd728b59e32cea99c6694ca64d946f7e00",
            "81cbd4c33acac95063994cc8de8b5bbef788c40ed23f8eceae690c9f0d14b602",
            False, id="api50",
        ),
    ],
)
def test_cuda_oracle_preserves_approved_source_identity(
    producer_name, fixture_name, expected_producer, expected_fixture, embedded_hash,
):
    package = Path(__file__).resolve().parents[1]
    inputs = package / "oracles"
    producer = (inputs / producer_name).read_bytes()
    producer_hash = hashlib.sha256(producer.replace(b"\r\n", b"\n")).hexdigest()
    assert producer_hash == expected_producer, (producer_hash, expected_producer)
    fixture = (inputs / fixture_name).read_bytes()
    fixture_hash = hashlib.sha256(fixture).hexdigest()
    assert fixture_hash == expected_fixture, (fixture_hash, expected_fixture)
    oracle = json.loads(fixture)
    if embedded_hash:
        assert oracle["generator_sha256"] == producer_hash, oracle["generator_sha256"]
    else:
        # The unchanged legacy schema has no embedded producer fingerprint.
        assert "generator_sha256" not in oracle, oracle.keys()
    assert oracle["optuna_version"] == "5.0.0", oracle["optuna_version"]
    project = tomllib.loads((package / "pyproject.toml").read_text(encoding="utf-8"))
    upstream = project["project"]["optional-dependencies"]["upstream"]
    assert f"optuna=={oracle['optuna_version']}" in upstream, upstream
