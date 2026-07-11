from __future__ import annotations

import copy
import json

from radia_mcp.build123d.server import build123d_curved_step_topology_crosscheck_gate


def _summary() -> dict:
    native_volume = 645.0537
    external_volume = native_volume * (1.0 + 1.15e-5)
    native = {
        "volume": native_volume,
        "area": 1497.8833,
        "bbox_size": [18.0, 18.0, 8.44],
        "solids": 1,
        "faces": 70,
        "edges": 171,
        "is_valid": True,
    }
    roundtrip = copy.deepcopy(native)
    roundtrip["volume"] *= 1.0 - 2.2e-8
    roundtrip["area"] *= 1.0 - 2.1e-8
    return {
        "source_kind": "upstream_native_example",
        "upstream_commit": "a" * 40,
        "source_sha256": "b" * 64,
        "build123d_version": "0.10.0",
        "step_sha256": "c" * 64,
        "native": native,
        "roundtrip": roundtrip,
        "external_imports": [
            {
                "mode": mode,
                "volume": external_volume,
                "body_count": 1,
                "volume_count": 1,
                "surface_count": 70,
                "curve_count": 171,
                "step_sha256": "c" * 64,
            }
            for mode in ("noheal", "heal")
        ],
        "external_volume_bias_classification": "cross_kernel_curved_surface_translation",
        "external_volume_bias_tolerance_basis": "curved STEP with exact body, face, edge, and import-mode invariants",
        "timings_s": {
            "source_build": 2.4,
            "step_export": 0.03,
            "step_reimport": 0.1,
            "external_import": 0.2,
        },
    }


def _gate(summary: dict) -> dict:
    return json.loads(build123d_curved_step_topology_crosscheck_gate(json.dumps(summary)))


def test_accepts_classified_curved_step_bias_with_exact_topology():
    result = _gate(_summary())
    assert result["status"] == "ok"
    assert result["diagnosis"] == "portable_with_classified_cross_kernel_bias"
    assert result["metrics"]["maximum_external_volume_relative_error"] > 1.0e-5
    assert result["metrics"]["external_import_mode_volume_spread"] == 0.0


def test_rejects_external_cavity_topology_loss_even_with_volume_in_tolerance():
    summary = _summary()
    summary["external_imports"][1]["surface_count"] = 69
    result = _gate(summary)
    assert result["status"] == "needs_attention"
    assert result["checks"]["external_face_edge_topology_matches"] is False


def test_rejects_unclassified_bias_and_heal_noheal_disagreement():
    summary = _summary()
    summary["external_volume_bias_classification"] = ""
    summary["external_imports"][1]["volume"] *= 1.0 + 5.0e-6
    result = _gate(summary)
    assert result["status"] == "needs_attention"
    assert result["checks"]["external_import_modes_volume_invariant"] is False
    assert result["checks"]["cross_kernel_volume_bias_explicitly_classified"] is False
