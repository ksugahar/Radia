import json

from radia_mcp.document_meta import PUBLIC_TOOLS
from radia_mcp.document_meta.tools import document_meta_notebook_result_audit


def _write_notebook(path, *, executed=True, outputs=True, source="print('ok')\n",
                    metadata=None, rich_output=False):
    path.parent.mkdir(parents=True, exist_ok=True)
    cell = {
        "cell_type": "code",
        "execution_count": 1 if executed else None,
        "metadata": {},
        "outputs": (
            [{
                "output_type": "display_data",
                "metadata": {},
                "data": {
                    "application/vnd.jupyter.widget-view+json": {
                        "model_id": "webgui-scene",
                        "version_major": 2,
                        "version_minor": 0,
                    },
                    "text/plain": ["WebGUI scene"],
                },
            }]
            if outputs and rich_output
            else [{"output_type": "stream", "name": "stdout", "text": ["ok\n"]}]
            if outputs
            else []
        ),
        "source": [source],
    }
    path.write_text(json.dumps({
        "cells": [cell],
        "metadata": metadata or {},
        "nbformat": 4,
        "nbformat_minor": 5,
    }), encoding="utf-8")


def test_public_surface_excludes_retired_repository_migration_tools():
    assert PUBLIC_TOOLS == (
        "document_meta_deadline_countdown",
        "document_meta_diff_versions",
        "document_meta_template_loader",
        "document_meta_lint_all",
        "document_meta_notebook_result_audit",
    )


def test_notebook_result_audit_needs_saved_outputs_but_not_json(tmp_path):
    docs = tmp_path / "docs" / "demo"
    nb = docs / "demo.ipynb"
    _write_notebook(nb, executed=True, outputs=True)

    ready = document_meta_notebook_result_audit(str(tmp_path), "docs")
    assert ready["summary"]["notebooks_scanned"] == 1
    assert ready["summary"]["ok_result_saved"] == 1
    assert ready["gaps"] == []
    assert not list(docs.glob("*.json"))

    _write_notebook(nb, executed=False, outputs=False)
    missing = document_meta_notebook_result_audit(str(tmp_path), "docs")
    assert missing["summary"]["needs_saved_outputs"] == 1
    assert missing["gaps"][0]["status"] == "needs_saved_outputs"


def test_default_notebook_audit_does_not_parse_adjacent_json(tmp_path):
    docs = tmp_path / "docs" / "demo"
    nb = docs / "demo.ipynb"
    _write_notebook(nb, executed=True, outputs=True)
    (docs / "demo_results.json").write_text("{broken", encoding="utf-8")

    ready = document_meta_notebook_result_audit(str(tmp_path), "docs")
    assert ready["summary"]["ok_result_saved"] == 1
    assert ready["notebooks"][0]["result_json_count"] == 0

    compatibility = document_meta_notebook_result_audit(
        str(tmp_path), "docs", require_json=True
    )
    assert compatibility["summary"]["result_json_missing_version_or_date"] == 1


def test_notebook_result_audit_requires_saved_webgui_for_examples(tmp_path):
    nb = tmp_path / "docs" / "demo" / "demo.ipynb"
    metadata = {"radia": {"notebook_role": "example", "webgui_required": True}}
    _write_notebook(nb, metadata=metadata)

    missing = document_meta_notebook_result_audit(str(tmp_path), "docs")
    assert missing["summary"]["needs_webgui_draw"] == 1

    _write_notebook(
        nb,
        source="from ngsolve.webgui import Draw\nDraw(mesh)\n",
        metadata=metadata,
        rich_output=True,
    )
    ready = document_meta_notebook_result_audit(str(tmp_path), "docs")
    assert ready["summary"]["webgui_ready"] == 1
    assert ready["gaps"] == []


def test_notebook_result_audit_requires_parameterized_field_scene(tmp_path):
    nb = tmp_path / "docs" / "demo" / "demo.ipynb"
    metadata = {"radia": {
        "notebook_role": "example",
        "webgui_required": True,
        "webgui_field_required": True,
    }}
    _write_notebook(
        nb,
        source="from ngsolve.webgui import Draw\nDraw(field, mesh)\n",
        metadata=metadata,
        rich_output=True,
    )
    missing = document_meta_notebook_result_audit(str(tmp_path), "docs")
    assert missing["gaps"][0]["status"] == "needs_parameterized_webgui_field_draw"

    _write_notebook(
        nb,
        source=(
            "from ngsolve.webgui import Draw\n"
            "Draw(field, mesh, name='B_magnitude', draw_vol=True)\n"
        ),
        metadata=metadata,
        rich_output=True,
    )
    ready = document_meta_notebook_result_audit(str(tmp_path), "docs")
    assert ready["summary"]["webgui_field_ready"] == 1
    assert ready["gaps"] == []
