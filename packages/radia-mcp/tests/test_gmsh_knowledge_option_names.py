"""Self-verification: every gmsh option named in the knowledge base exists.

Extracts option names from gmsh_knowledge.py / gmsh_reference.py source
(assignment lines in code blocks and backticked table cells) and probes
them against the installed gmsh's option database in one subprocess.
A failure means the knowledge has drifted from gmsh reality -- fix the
text or, for an intentionally-invalid example, extend the allowlist.
"""

import importlib.util
import re
from pathlib import Path

import pytest

from radia_mcp.gmsh.msh_inspect import INVALID_GEO_OPTIONS

_GMSH_AVAILABLE = importlib.util.find_spec("gmsh") is not None

pytestmark = pytest.mark.skipif(not _GMSH_AVAILABLE,
                                reason="gmsh package not installed")

_KNOWLEDGE_DIR = (Path(__file__).parent.parent / "src" / "radia_mcp" / "gmsh")

# Only option roots that gmsh actually has -- keeps prose like
# "e.g. Foo.Bar" or file names out of the extraction.
_ROOTS = ("General", "Mesh", "View", "PostProcessing", "Geometry",
          "Solver", "Print")

_ASSIGN_RE = re.compile(
    rf"^\s*((?:{'|'.join(_ROOTS)})(?:\.[A-Za-z0-9_]+|\[\d+\]\.[A-Za-z0-9_]+)+)"
    rf"\s*=", re.MULTILINE)
_TABLE_RE = re.compile(
    rf"`((?:{'|'.join(_ROOTS)})(?:\.[A-Za-z0-9_]+|\[n\]\.[A-Za-z0-9_]+"
    rf"|\[\d+\]\.[A-Za-z0-9_]+)+)`")

# Names the knowledge mentions ON PURPOSE as removed/renamed examples.
_ALLOWLIST = set(INVALID_GEO_OPTIONS) | {
    "General.ConfirmQuit",  # removed from gmsh; cited as a stale-name example
    "View.ArrowScale",      # pre-4.x name; cited as a rename example
}


def _extract_names() -> set[str]:
    names: set[str] = set()
    for module in ("gmsh_knowledge.py", "gmsh_reference.py"):
        text = (_KNOWLEDGE_DIR / module).read_text(encoding="utf-8",
                                                   errors="replace")
        for match in _ASSIGN_RE.finditer(text):
            names.add(match.group(1))
        for match in _TABLE_RE.finditer(text):
            names.add(match.group(1).replace("[n]", "[0]"))
    return names


def test_every_documented_option_name_exists_in_gmsh():
    from radia_mcp.gmsh.msh_inspect import probe_options

    names = sorted(_extract_names() - _ALLOWLIST)
    assert len(names) > 50, f"extraction looks broken: only {names}"

    result = probe_options(names, timeout_s=300.0)
    assert result["ran"] is True, result.get("error")
    assert result["missing"] == [], (
        "knowledge names gmsh options that do not exist in the installed "
        f"gmsh -- fix the docs or extend the allowlist: {result['missing']}")
