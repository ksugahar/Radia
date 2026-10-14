"""Knowledge-drift gates for radia_mcp.build123d (gmsh-pattern rollout).

Same idea as the gmsh option-name self-verification: the knowledge a
server ships must match the reality of the installed libraries.

1. The auto-generated build123d API snapshot (build123d_api.py) names
   module-level functions -- every one must still exist in the
   INSTALLED build123d (a vanished name means the snapshot must be
   regenerated via `python -m radia_mcp.build123d._gen_api_reference`).
2. The hand-written knowledge references `modeling.xxx` /
   `archetypes.xxx` lab helpers -- each must exist in the actual
   modules.
"""

import re

import pytest

build123d = pytest.importorskip("build123d")

from radia_mcp.build123d import archetypes, modeling  # noqa: E402
from radia_mcp.build123d.build123d_api import API_REFERENCE  # noqa: E402
from radia_mcp.build123d.build123d_knowledge import (  # noqa: E402
    get_build123d_documentation,
)

_FUNC_RE = re.compile(r"^### ([a-z_][a-z0-9_]*)\((?!self)", re.MULTILINE)
_HELPER_RE = re.compile(r"\b(modeling|archetypes)\.([a-z_][a-z0-9_]*)")


def test_api_snapshot_functions_exist_in_installed_build123d():
    names = sorted(set(_FUNC_RE.findall(API_REFERENCE)))
    assert len(names) > 20, f"extraction looks broken: {names}"

    missing = [n for n in names if not hasattr(build123d, n)]
    assert missing == [], (
        "API snapshot names functions that the installed build123d "
        f"({build123d.__version__}) no longer exports -- regenerate via "
        f"python -m radia_mcp.build123d._gen_api_reference: {missing}")


def test_knowledge_helper_references_exist():
    text = get_build123d_documentation("all")
    refs = set(_HELPER_RE.findall(text))
    assert refs, "no modeling./archetypes. references found in knowledge"

    modules = {"modeling": modeling, "archetypes": archetypes}
    missing = sorted(
        f"{mod}.{name}" for mod, name in refs
        if not hasattr(modules[mod], name))
    assert missing == [], (
        "knowledge references lab helpers that do not exist: "
        f"{missing}")
