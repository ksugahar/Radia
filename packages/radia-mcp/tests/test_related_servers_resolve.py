"""Status-tool related_servers name existing catalog keys.

Guards: every literal ``related_servers=[...]`` passed in a radia-mcp
``server.py`` names a current catalog key or an external package, so
renamed or merged servers cannot linger in discovery output.
Assumes: related_servers are written as list literals of strings.
"""

from __future__ import annotations

import ast
from pathlib import Path

from radia_mcp.meta.catalog import CATALOG, EXTERNAL_PACKAGES

SRC = Path(__file__).resolve().parents[1] / "src" / "radia_mcp"


def _related_literals():
    for path in sorted(SRC.rglob("server.py")):
        tree = ast.parse(path.read_text(encoding="utf-8"))
        for node in ast.walk(tree):
            if not isinstance(node, ast.keyword) or node.arg != "related_servers":
                continue
            if not isinstance(node.value, ast.List):
                continue
            for elt in node.value.elts:
                if isinstance(elt, ast.Constant) and isinstance(elt.value, str):
                    yield path.relative_to(SRC).as_posix(), node.value.lineno, elt.value


def test_related_servers_are_catalog_keys():
    known = set(CATALOG) | set(EXTERNAL_PACKAGES)
    unknown = [item for item in _related_literals() if item[2] not in known]
    assert unknown == [], unknown
