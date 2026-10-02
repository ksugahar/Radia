"""Contract: repository input paths quoted in knowledge and prompts exist."""

import re
from pathlib import Path

import pytest

PACKAGE = Path(__file__).resolve().parents[1]
REPO = PACKAGE.parents[1]
SOURCE = PACKAGE / "src" / "radia_mcp"
PATH = re.compile(
    r"(?:src/radia/panels/samples|validation_test)/[A-Za-z0-9_./-]+"
    r"\.(?:vol|jou|step|stp|txt|json|py|csv|msh)\b")


def _knowledge_files():
    for path in sorted(SOURCE.rglob("*")):
        if not path.is_file():
            continue
        if (path.suffix == ".py" and "knowledge" in path.name) or (
                path.suffix == ".md" and path.parent.name == "prompts"):
            yield path


def test_quoted_repository_paths_exist():
    if not (REPO / "validation_test").is_dir():
        pytest.skip("repository checkout required")
    missing = sorted(
        f"{path.relative_to(SOURCE)}: {ref}"
        for path in _knowledge_files()
        for ref in set(PATH.findall(path.read_text(encoding="utf-8")))
        if not (REPO / ref).exists())
    assert not missing, "\n".join(missing)
