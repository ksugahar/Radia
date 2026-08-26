import asyncio
import json
import os
import sys
from pathlib import Path

from mcp import ClientSession, StdioServerParameters
from mcp.client.stdio import stdio_client

from radia_mcp.cubit.nastran_consumer import evaluate_nastran_consumer_contract
from radia_mcp.cubit.server import (
    cubit_docs,
    cubit_export_decision_guide,
    cubit_nastran_consumer_gate,
    mcp,
)


def _summary():
    return {
        "producer": {
            "export_command": 'export nastran_bdf "mesh.bdf" order 2 dimension 3 overwrite',
            "artifact_scope": "mesh_interchange",
            "independent_parse_ok": True,
            "bdf_sha256": "a" * 64,
            "dimension": 3,
            "order": 2,
            "expected_nodes": 7353,
            "expected_linear_nodes": 1032,
            "expected_primary_elements": 4902,
        },
        "consumer": {
            "import_succeeded": True,
            "imported_dimension": 3,
            "imported_nodes": 7353,
            "imported_primary_elements": 4902,
            "has_second_order_elements": True,
            "remeshed": False,
            "study_has_imported_mesh": True,
            "bbox_max_abs_error": 0.0,
            "material_assignment_available": True,
        },
    }


def test_nastran_consumer_contract_accepts_identity_preserving_handoff():
    result = evaluate_nastran_consumer_contract(_summary())

    assert result["status"] == "pass"
    assert result["passed"] is True
    assert result["producer_valid"] is True
    assert result["consumer_valid"] is True
    assert all(result["checks"].values())


def test_nastran_consumer_contract_rejects_builtin_export_contract():
    summary = _summary()
    summary["producer"]["export_command"] = 'export nastran "mesh.bdf" overwrite everything'

    result = evaluate_nastran_consumer_contract(summary)

    assert result["status"] == "wrong_exporter"
    assert result["checks"]["canonical_export_command"] is False


def test_nastran_consumer_contract_detects_second_order_downgrade():
    summary = _summary()
    summary["consumer"]["imported_nodes"] = 1032
    summary["consumer"]["has_second_order_elements"] = False

    result = evaluate_nastran_consumer_contract(summary)

    assert result["status"] == "consumer_order_downgrade"
    assert result["checks"]["order_preserved"] is False
    assert result["checks"]["primary_element_count_preserved"] is True


def test_nastran_consumer_contract_detects_unretained_3d_mesh():
    summary = _summary()
    summary["consumer"].update(
        study_has_imported_mesh=False,
        imported_nodes=0,
        imported_primary_elements=0,
    )

    result = evaluate_nastran_consumer_contract(summary)

    assert result["status"] == "consumer_mesh_not_retained"
    assert result["checks"]["study_retains_imported_mesh"] is False


def test_nastran_consumer_contract_rejects_post_import_remesh_counts():
    summary = _summary()
    summary["consumer"]["remeshed"] = True

    result = evaluate_nastran_consumer_contract(summary)

    assert result["status"] == "consumer_remeshed"
    assert "before any consumer-side" in result["recommendation"]


def test_nastran_consumer_contract_rejects_complete_deck_claim():
    summary = _summary()
    summary["producer"]["artifact_scope"] = "complete_analysis_deck"

    result = evaluate_nastran_consumer_contract(summary)

    assert result["status"] == "scope_mismatch"
    assert result["checks"]["mesh_interchange_scope"] is False


def test_nastran_consumer_contract_marks_legacy_alias_for_migration():
    summary = _summary()
    summary["producer"]["export_command"] = 'export jmag_nastran "mesh.bdf" order 2 overwrite'

    result = evaluate_nastran_consumer_contract(summary)

    assert result["status"] == "legacy_alias"
    assert result["passed"] is False
    assert "deprecated compatibility alias" in result["warnings"][0]


def test_nastran_consumer_contract_requires_independent_parser_and_digest():
    summary = _summary()
    summary["producer"]["independent_parse_ok"] = False
    summary["producer"]["bdf_sha256"] = "not-a-digest"

    result = evaluate_nastran_consumer_contract(summary)

    assert result["status"] == "producer_invalid"
    assert result["checks"]["independent_parse_ok"] is False
    assert result["checks"]["digest_recorded"] is False


def test_nastran_consumer_contract_optional_set_semantics_is_explicit():
    summary = _summary()

    result = evaluate_nastran_consumer_contract(
        summary, require_set_semantics=True
    )

    assert result["status"] == "consumer_set_semantics_unverified"
    summary["consumer"]["set_semantics_verified"] = True
    assert evaluate_nastran_consumer_contract(
        summary, require_set_semantics=True
    )["status"] == "pass"


def test_nastran_consumer_mcp_wrapper_and_invalid_input_contract():
    payload = json.loads(cubit_nastran_consumer_gate(_summary()))
    invalid = json.loads(cubit_nastran_consumer_gate({}))

    assert payload["status"] == "pass"
    assert invalid["status"] == "invalid_input"
    assert invalid["passed"] is False
    assert "producer must be an object" in invalid["error"]


def test_nastran_mcp_knowledge_uses_tool_neutral_command_and_gate():
    docs = cubit_docs("nastran")
    resource = cubit_export_decision_guide()

    assert "export nastran_bdf" in docs
    assert "compatibility alias" in docs
    assert "mesh-interchange artifact" in docs
    assert "cubit_nastran_consumer_gate" in docs
    assert "export nastran_bdf" in resource
    assert "cubit_nastran_consumer_gate" in resource


def test_nastran_consumer_mcp_tool_has_explicit_contract_metadata():
    tool = mcp._tool_manager._tools["cubit_nastran_consumer_gate"]

    assert tool.title == "Validate Nastran Mesh Consumer Contract"
    assert tool.annotations.readOnlyHint is True
    assert tool.annotations.destructiveHint is False
    assert tool.annotations.idempotentHint is True
    assert tool.annotations.openWorldHint is False


async def _probe_nastran_gate_stdio():
    package_root = Path(__file__).resolve().parents[1]
    env = os.environ.copy()
    env["PYTHONPATH"] = os.pathsep.join(
        [str(package_root / "src"), env.get("PYTHONPATH", "")]
    ).rstrip(os.pathsep)
    params = StdioServerParameters(
        command=sys.executable,
        args=["-m", "radia_mcp.cubit.server"],
        cwd=str(package_root),
        env=env,
    )
    negative = _summary()
    negative["consumer"]["imported_nodes"] = 1032
    negative["consumer"]["has_second_order_elements"] = False
    async with stdio_client(params) as (read_stream, write_stream):
        async with ClientSession(read_stream, write_stream) as session:
            initialized = await session.initialize()
            listed = await session.list_tools()
            tool = next(
                item
                for item in listed.tools
                if item.name == "cubit_nastran_consumer_gate"
            )
            positive = await session.call_tool(
                "cubit_nastran_consumer_gate", {"summary": _summary()}
            )
            rejected = await session.call_tool(
                "cubit_nastran_consumer_gate", {"summary": negative}
            )
            return {
                "server_name": initialized.serverInfo.name,
                "server_version": initialized.serverInfo.version,
                "instructions": initialized.instructions,
                "title": tool.title,
                "annotations": tool.annotations,
                "positive": json.loads(positive.content[0].text),
                "positive_error": bool(positive.isError),
                "negative": json.loads(rejected.content[0].text),
                "negative_error": bool(rejected.isError),
            }


def test_nastran_gate_passes_real_stdio_initialize_list_and_calls():
    result = asyncio.run(asyncio.wait_for(_probe_nastran_gate_stdio(), timeout=45))

    assert "cubit" in result["server_name"].lower()
    assert result["server_version"]
    assert result["instructions"]
    assert result["title"] == "Validate Nastran Mesh Consumer Contract"
    assert result["annotations"].readOnlyHint is True
    assert result["positive_error"] is False
    assert result["positive"]["status"] == "pass"
    assert result["negative_error"] is False
    assert result["negative"]["status"] == "consumer_order_downgrade"
