from __future__ import annotations

import asyncio
import json

from radia_mcp.matlab_agentic_ml import validate_aicia_catalog


def _catalog():
    required = [
        "seed_or_determinism",
        "units_and_schema",
        "provenance",
        "independent_forward_solver_verification",
    ]
    prior = [
        {
            "id": f"prior-{index:03d}",
            "previously_verified": True,
            "disposition": "previously_verified",
        }
        for index in range(244)
    ]
    pending = []
    for disposition, count in (("candidate", 46), ("review", 14),
                               ("not_promoted", 93)):
        for index in range(count):
            item = {
                "id": f"{disposition}-{index:03d}",
                "previously_verified": False,
                "disposition": disposition,
            }
            if disposition == "candidate":
                item["promotion_requirements"] = required.copy()
            pending.append(item)
    return {
        "channel_id": "UC2lJYodMaAfFeFQrGUwhlaQ",
        "policy": {
            "metadata_only": True,
            "transcripts_downloaded": False,
            "media_downloaded": False,
            "generated_candidate_is_ground_truth": False,
            "promotion_requires_forward_solver_verification": True,
        },
        "counts": {
            "videos": 245,
            "shorts": 2,
            "streams": 150,
            "total": 397,
            "unique_ids": 397,
            "previously_verified": 244,
            "processed_now": 153,
            "candidate": 46,
            "review": 14,
            "not_promoted": 93,
        },
        "items": prior + pending,
    }


def test_all_153_new_channel_items_are_dispositioned():
    result = validate_aicia_catalog(_catalog())
    assert result["status"] == "ok"
    assert result["processed_now"] == 153
    assert result["candidate_count"] == 46


def test_videos_only_snapshot_is_rejected():
    catalog = _catalog()
    catalog["counts"]["streams"] = 0
    catalog["counts"]["total"] = 247
    result = validate_aicia_catalog(catalog)
    assert result["status"] == "needs_attention"
    assert result["checks"]["full_channel_scope_is_397"] is False


def test_candidate_without_forward_solver_gate_is_rejected():
    catalog = _catalog()
    candidate = next(item for item in catalog["items"] if item["disposition"] == "candidate")
    candidate["promotion_requirements"].remove("independent_forward_solver_verification")
    result = validate_aicia_catalog(catalog)
    assert result["status"] == "needs_attention"
    assert result["checks"]["candidate_promotions_keep_solver_gates"] is False


def test_matlab_mcp_exposes_aicia_catalog_gate():
    from radia_mcp.matlab.server import matlab_aicia_catalog_gate, mcp

    names = {tool.name for tool in asyncio.run(mcp.list_tools())}
    assert "matlab_validation_catalog" in names
    catalog = mcp._tool_manager._tools["matlab_validation_catalog"].fn()
    assert "matlab_aicia_catalog_gate" in {
        item["name"] for item in catalog["operations"]
    }
    result = json.loads(matlab_aicia_catalog_gate(json.dumps(_catalog())))
    assert result["status"] == "ok"
