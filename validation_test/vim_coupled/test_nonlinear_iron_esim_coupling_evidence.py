"""Evidence gate for the same-region magnetic-conductor local-ESIM lane."""

import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
ARTIFACT = HERE / "results_nonlinear_iron_esim_coupling.json"


def _artifact():
    return json.loads(ARTIFACT.read_text(encoding="utf-8"))


def test_nonlinear_iron_esim_artifact_has_reproducible_identity_and_versions():
    artifact = _artifact()
    assert artifact["schema"] == "radia.validation.nonlinear-iron-local-esim-coupling.v1"
    assert artifact["created_at_utc"]
    assert artifact["tool_versions"]["radia"] != "unknown"
    assert artifact["tool_versions"]["ngsolve"] != "unknown"
    assert len(artifact["source_fingerprints"]) >= 4
    assert artifact["identity"]["magnetic_region"] == "iron"
    assert artifact["identity"]["conductive_region"] == "iron"
    assert artifact["identity"]["mesh"].endswith("no mesh tracked")


def test_nonlinear_iron_esim_artifact_passes_physics_and_replay_gates():
    artifact = _artifact()
    assert artifact["pass"]
    assert all(artifact["checks"].values())
    rows = artifact["amplitude_ladder"]
    assert len(rows) == 3
    assert all(row["local_esim"]["converged"] for row in rows)
    assert all(row["local_surface_impedance"]["passive"] for row in rows)
    assert max(row["mixed"]["residual_relative_norm"] for row in rows) < 1.0e-10
    assert max(
        row["mixed"]["fixed_gram_replay_relative_difference"] for row in rows
    ) < 1.0e-12
    assert min(row["mixed"]["joule_loss_W"] for row in rows) > 0.0
    assert "bulk nonlinear B-H" in artifact["claim_boundary"]["not_established"]
