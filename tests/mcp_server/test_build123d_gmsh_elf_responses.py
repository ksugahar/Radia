"""Response regression tests for build123d / gmsh / elf knowledge tools.

Each server exposes a single ``get_*_documentation`` function dispatching
on a topic string.  These tests guard against:

  * accidental deletion of an advertised topic (silent rename, typo)
  * empty or truncated topic body (lost content)
  * unknown-topic fallback crashing / returning empty
  * 'all' concatenation regressions (returns only the first entry)

Run::

    pytest tests/mcp_server/test_build123d_gmsh_elf_responses.py -v
"""
from __future__ import annotations

import pytest

from radia_mcp.build123d.build123d_knowledge import get_build123d_documentation
from radia_mcp.elf.elf_knowledge import get_elf_documentation
from radia_mcp.gmsh.gmsh_knowledge import get_gmsh_documentation


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _assert_substantial(body: str, *, min_chars: int = 200) -> None:
    assert isinstance(body, str), f"expected str, got {type(body).__name__}"
    assert len(body) >= min_chars, (
        f"topic body is too short ({len(body)} chars); "
        f"likely empty/truncated: {body[:100]!r}")


def _assert_contains_any(body: str, markers: list[str], label: str) -> None:
    low = body.lower()
    hits = [m for m in markers if m.lower() in low]
    assert hits, (
        f"{label}: none of {markers} found in the {len(body)}-char body.\n"
        f"First 200 chars:\n{body[:200]!r}")


def _assert_unknown_topic_graceful(fn, *, unknown: str = "definitely_not_a_topic"):
    r = fn(unknown)
    assert isinstance(r, str)
    assert len(r) > 0
    low = r.lower()
    assert "unknown" in low or "available" in low or unknown.lower() in low, (
        f"{fn.__name__}({unknown!r}) does not mention unknown/available: {r[:200]!r}")


# ---------------------------------------------------------------------------
# build123d knowledge
# ---------------------------------------------------------------------------

class TestBuild123dKnowledge:
    TOPICS = {
        "overview":                     ["build123d", "cad"],
        "lab_policy":                   ["lab", "policy", "radia"],
        "cubit_rosetta":                ["cubit", "build123d"],
        "primitives_3d":                ["box", "cylinder", "sphere"],
        "primitives_2d":                ["rectangle", "circle"],
        "operations":                   ["fillet", "chamfer", "extrude", "boolean"],
        "curves":                       ["line", "arc", "spline"],
        "export_import":                ["step", "export", "import"],
        "topology":                     ["vertex", "edge", "face", "solid"],
        "coil_modeling":                ["coil", "helix", "sweep"],
        "plane_axis_location_cookbook": ["plane", "axis", "location"],
        "builder_vs_algebra_rosetta":   ["builder", "algebra"],
        "api_reference":                ["class", "method", "api"],
    }

    def test_all_topic(self):
        body = get_build123d_documentation("all")
        _assert_substantial(body, min_chars=5000)
        _assert_contains_any(body, ["build123d", "cad"], label="build123d all")

    @pytest.mark.parametrize("topic,markers", list(TOPICS.items()))
    def test_topic(self, topic, markers):
        body = get_build123d_documentation(topic)
        _assert_substantial(body, min_chars=100)
        _assert_contains_any(body, markers, label=f"build123d.{topic}")

    def test_unknown_topic(self):
        _assert_unknown_topic_graceful(get_build123d_documentation)


# ---------------------------------------------------------------------------
# gmsh knowledge
# ---------------------------------------------------------------------------

class TestGmshKnowledge:
    TOPICS = {
        "policy":     ["gmsh", "policy", "v4.1"],
        "overview":   ["gmsh", "mesh"],
        "cli":        ["gmsh", "-", "command"],
        "shortcuts":  ["shortcut", "key"],
        "options":    ["option"],
        "opt_file":   [".opt", "options"],
        "msh_format": ["msh", "format", "v4"],
        "geo":        ["geo", "point", "line"],
        "high_order": ["order", "high", "curving", "curve"],
        "workflow":   ["workflow", "radia"],
        "onelab":     ["onelab"],
        "pitfalls":   ["pitfall", "gotcha", "warning", "issue"],
        "animation":  ["animation", "animate", "frame"],
    }

    def test_all_topic(self):
        body = get_gmsh_documentation("all")
        _assert_substantial(body, min_chars=3000)
        _assert_contains_any(body, ["gmsh"], label="gmsh all")

    @pytest.mark.parametrize("topic,markers", list(TOPICS.items()))
    def test_topic(self, topic, markers):
        body = get_gmsh_documentation(topic)
        _assert_substantial(body, min_chars=100)
        _assert_contains_any(body, markers, label=f"gmsh.{topic}")

    def test_unknown_topic(self):
        _assert_unknown_topic_graceful(get_gmsh_documentation)


# ---------------------------------------------------------------------------
# elf knowledge
# ---------------------------------------------------------------------------

class TestElfKnowledge:
    TOPICS = {
        "overview":       ["elf"],
        "mai_format":     ["mai", "format"],
        "mei_format":     ["mei", "format"],
        "meg_format":     ["meg", "format"],
        "magic":          ["magic"],
        "elfin":          ["elfin"],
        "element_types":  ["element", "type"],
        "bh_curves":      ["b-h", "bh curve", "nonlinear"],
        "meg_export":     ["meg", "export"],
        "iemesh":         ["iemesh", "mesh"],
    }

    def test_all_topic(self):
        body = get_elf_documentation("all")
        _assert_substantial(body, min_chars=3000)
        _assert_contains_any(body, ["elf", "magic"], label="elf all")

    @pytest.mark.parametrize("topic,markers", list(TOPICS.items()))
    def test_topic(self, topic, markers):
        body = get_elf_documentation(topic)
        _assert_substantial(body, min_chars=100)
        _assert_contains_any(body, markers, label=f"elf.{topic}")

    def test_unknown_topic(self):
        _assert_unknown_topic_graceful(get_elf_documentation)
