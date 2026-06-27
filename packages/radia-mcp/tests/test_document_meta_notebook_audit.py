import json

from radia_mcp.document_meta.tools import (
    document_meta_examples_notebook_audit,
    document_meta_notebook_result_audit,
    document_meta_panel_layout_audit,
    document_meta_write_docs_notebook_result_jsons,
    document_meta_write_notebook_result_json,
)


def _write_notebook(path, *, executed=True, outputs=True):
    path.parent.mkdir(parents=True, exist_ok=True)
    cell = {
        "cell_type": "code",
        "execution_count": 1 if executed else None,
        "metadata": {},
        "outputs": (
            [{"output_type": "stream", "name": "stdout", "text": ["ok\n"]}]
            if outputs else []
        ),
        "source": ["print('ok')\n"],
    }
    nb = {
        "cells": [cell],
        "metadata": {},
        "nbformat": 4,
        "nbformat_minor": 5,
    }
    path.write_text(json.dumps(nb), encoding="utf-8")


def test_notebook_result_audit_requires_saved_outputs_and_versioned_json(tmp_path):
    repo = tmp_path
    (repo / "examples").mkdir()
    (repo / "docs" / "demo").mkdir(parents=True)
    nb = repo / "docs" / "demo" / "demo.ipynb"
    _write_notebook(nb, executed=True, outputs=True)
    document_meta_write_notebook_result_json(str(nb))

    result = document_meta_notebook_result_audit(str(repo), "docs")

    assert result["summary"]["notebooks_scanned"] == 1
    assert result["summary"]["ok_result_saved"] == 1
    assert result["gaps"] == []


def test_notebook_result_audit_flags_stale_json_sidecar(tmp_path):
    repo = tmp_path
    (repo / "examples").mkdir()
    docs = repo / "docs" / "demo"
    docs.mkdir(parents=True)
    nb = docs / "demo.ipynb"
    _write_notebook(nb, executed=True, outputs=True)
    document_meta_write_notebook_result_json(str(nb))

    # Simulate editing/re-executing the notebook without refreshing JSON.
    _write_notebook(nb, executed=True, outputs=True)
    data = json.loads(nb.read_text(encoding="utf-8"))
    data["cells"][0]["source"] = ["print('changed')\n"]
    nb.write_text(json.dumps(data), encoding="utf-8")

    result = document_meta_notebook_result_audit(str(repo), "docs")

    assert result["summary"]["result_json_not_synced_to_notebook"] == 1
    assert result["gaps"][0]["status"] == "result_json_not_synced_to_notebook"


def test_write_notebook_result_json_creates_versioned_sidecar(tmp_path):
    repo = tmp_path
    (repo / "examples").mkdir()
    docs = repo / "docs" / "demo"
    docs.mkdir(parents=True)
    nb = docs / "demo.ipynb"
    _write_notebook(nb, executed=True, outputs=True)

    write_result = document_meta_write_notebook_result_json(str(nb))
    assert "error" not in write_result
    sidecar = docs / "demo_result.json"
    assert sidecar.exists()

    data = json.loads(sidecar.read_text(encoding="utf-8"))
    assert data["schema"] == "radia.notebook_result.v1"
    assert data["generated_at_utc"].endswith("Z")
    assert data["versions"]["python_version"]
    assert data["summary"]["result_saved"] is True

    audit = document_meta_notebook_result_audit(str(repo), "docs")
    assert audit["summary"]["ok_result_saved"] == 1


def test_write_docs_notebook_result_jsons_skips_unexecuted(tmp_path):
    repo = tmp_path
    (repo / "examples").mkdir()
    docs = repo / "docs" / "demo"
    docs.mkdir(parents=True)
    _write_notebook(docs / "executed.ipynb", executed=True, outputs=True)
    _write_notebook(docs / "unexecuted.ipynb", executed=False, outputs=False)

    result = document_meta_write_docs_notebook_result_jsons(str(repo))

    assert result["summary"]["written_or_planned"] == 1
    assert result["summary"]["skipped"] == 1
    assert (docs / "executed_result.json").exists()
    assert not (docs / "unexecuted_result.json").exists()


def test_examples_notebook_audit_flags_unexecuted_notebook_and_py_review(tmp_path):
    repo = tmp_path
    topic = repo / "examples" / "demo"
    topic.mkdir(parents=True)
    (topic / "README.md").write_text("# demo\n", encoding="utf-8")
    (topic / "demo.py").write_text("print('demo')\n", encoding="utf-8")
    (repo / "docs" / "demo").mkdir(parents=True)
    _write_notebook(repo / "docs" / "demo" / "demo.ipynb",
                    executed=False, outputs=False)

    result = document_meta_examples_notebook_audit(str(repo), "demo")
    row = result["topics"][0]

    assert row["status"] == "notebook_without_saved_results"
    assert row["python_policy"] == "review_each_example_py_for_docs_helper_or_src_api"
    assert result["notebooks_without_saved_results"]


def test_examples_notebook_audit_reports_validation_test_lane(tmp_path):
    repo = tmp_path
    topic = repo / "examples" / "demo"
    topic.mkdir(parents=True)
    (topic / "validation_demo.py").write_text("print('validation')\n", encoding="utf-8")
    vtest = repo / "validation_test" / "demo"
    vtest.mkdir(parents=True)
    (vtest / "test_demo.py").write_text(
        'SOURCE = "examples/demo/validation_demo.py"\n',
        encoding="utf-8",
    )
    (repo / "docs").mkdir()

    result = document_meta_examples_notebook_audit(str(repo), "demo")
    row = result["topics"][0]

    assert row["validation_class"] is True
    assert row["validation_test_reference_count"] == 1
    assert row["promotion_lane"] == "validation_test_or_keep_validation_corpus"
    assert result["summary"]["validation_class_topics"] == 1
    assert result["validation_reviews"][0]["topic"] == "demo"


def test_examples_notebook_audit_uses_kelvin_docs_alias(tmp_path):
    repo = tmp_path
    topic = repo / "examples" / "kelvin_transformation"
    topic.mkdir(parents=True)
    (topic / "demo.py").write_text("print('kelvin')\n", encoding="utf-8")
    docs = repo / "docs" / "kelvin"
    docs.mkdir(parents=True)
    _write_notebook(docs / "kelvin_examples_migration.ipynb",
                    executed=True, outputs=True)

    result = document_meta_examples_notebook_audit(str(repo), "kelvin_transformation")
    row = result["topics"][0]

    assert row["docs_notebook_count"] == 1
    assert row["result_saved_notebook_count"] == 1
    assert row["status"] == "notebook_result_saved"


def test_panel_layout_audit_reports_old_path_references(tmp_path):
    repo = tmp_path
    (repo / "examples").mkdir()
    (repo / "docs").mkdir()
    panel_dir = repo / "src" / "radia" / "panels"
    panel_dir.mkdir(parents=True)
    (panel_dir / "calc_demo.py").write_text("print('panel')\n", encoding="utf-8")
    (repo / "pyproject.toml").write_text(
        '[tool.demo]\npath = "src/radia/panels/calc_demo.py"\n',
        encoding="utf-8",
    )

    result = document_meta_panel_layout_audit(str(repo))

    assert result["current_counts"]["calc_scripts"] == 1
    assert result["old_path_reference_count_reported"] == 1
    assert result["readiness"] == "needs_staged_migration_or_compat_shim"
