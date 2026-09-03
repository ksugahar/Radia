"""Fast structural checks for the public, result-bearing notebook corpus."""

from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_docs_notebooks_are_parseable_and_free_of_replacement_glyphs():
    notebooks = sorted((ROOT / "docs").rglob("*.ipynb"))
    assert notebooks

    for path in notebooks:
        data = json.loads(path.read_text(encoding="utf-8"))
        assert data.get("nbformat", 0) >= 4, path
        assert data.get("cells"), path
        assert "\ufffd" not in json.dumps(data, ensure_ascii=False), path


def test_public_example_notebooks_keep_their_saved_scene_contract():
    for path in sorted((ROOT / "docs").rglob("*.ipynb")):
        data = json.loads(path.read_text(encoding="utf-8"))
        radia = data.get("metadata", {}).get("radia", {})
        if radia.get("notebook_role") != "example":
            continue

        assert radia.get("webgui_required") is True, path
        if radia.get("webgui_field_required"):
            assert radia.get("webgui_required") is True, path
        assert any(
            cell.get("outputs")
            for cell in data.get("cells", [])
            if cell.get("cell_type") == "code"
        ), path


def test_docs_do_not_own_validation_or_benchmark_json():
    offenders = [
        path.relative_to(ROOT)
        for path in (ROOT / "docs").rglob("*.json")
        if "validation" in path.name.lower() or "benchmark" in path.name.lower()
    ]
    assert not offenders, (
        "Move validation and benchmark JSON records to validation_test/: "
        + ", ".join(map(str, offenders))
    )
