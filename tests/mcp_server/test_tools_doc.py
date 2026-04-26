"""Drift check for packages/radia-mcp/docs/TOOLS.md.

Fails if the committed file no longer matches what
`packages/radia-mcp/scripts/gen_tools_doc.py` would generate now.
Forces the doc to be regenerated whenever a tool is added/renamed/removed.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent.parent
SCRIPTS = REPO / "packages" / "radia-mcp" / "scripts"


@pytest.fixture(scope="module", autouse=True)
def _add_scripts_to_path():
    sys.path.insert(0, str(SCRIPTS))
    yield
    try:
        sys.path.remove(str(SCRIPTS))
    except ValueError:
        pass


def test_radia_mcp_tools_doc_up_to_date():
    import gen_tools_doc

    assert gen_tools_doc.check(), (
        "packages/radia-mcp/docs/TOOLS.md is stale. "
        "Regenerate it: `python packages/radia-mcp/scripts/gen_tools_doc.py`"
    )
