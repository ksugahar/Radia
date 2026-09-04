"""Response regression tests for mcp-server-cubit knowledge tools.

Covers every advertised ``get_*`` knowledge function in the cubit MCP
server package and asserts:

  * 'all' (or default topic) returns non-empty substantial text with
    topic-relevant markers
  * every documented sub-topic returns non-empty content with its own
    specific markers (catches accidental topic deletion / empty body
    / silent rename)
  * unknown topics degrade gracefully (no KeyError / crash)

Run::

    pytest tests/mcp_server/test_cubit_responses.py -v
"""
from __future__ import annotations

import pytest

from radia_mcp.cubit.knowledge.scripting import get_cubit_documentation
from radia_mcp.cubit.knowledge.custom_toolbar import get_toolbar_documentation
from radia_mcp.cubit.knowledge.export import get_export_documentation
from radia_mcp.cubit.knowledge.mesh_diagnostics import (
    get_diagnostics_documentation,
)
from radia_mcp.cubit.knowledge.netgen_workflow import get_netgen_documentation


# ---------------------------------------------------------------------------
# Helpers (duplicated intentionally so each test file is standalone)
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
# cubit_scripting_knowledge — the big one (~40 topics + many aliases)
# ---------------------------------------------------------------------------

class TestCubitScriptingKnowledge:
    # A small curated set of critical topics (one marker must appear in the
    # topic body). Keep this list high-signal; don't try to enumerate every
    # alias.
    TOPICS = {
        "overview":            ["Cubit", "mesh", "geometry"],
        "lab_policy":          ["lab", "policy", "convention"],
        "boolean_policy":      ["boolean", "unite", "imprint", "merge"],
        "trial_error_policy":  ["batch", "nographics", "commit"],
        "trampoline_pitfalls": ["trampoline", "pitfall", "acis", "curving"],
        "blocks":              ["block", "sideset"],
        "hex_workflow":        ["hex", "sweep"],
        "tet_workflow":        ["tet", "netgen"],
        "kelvin_transform":    ["Kelvin", "transform"],
        "ngsolve_panel":       ["Simulink", "check-vol", "DesignSpec"],
        "aprepro_journal":     ["APREPRO", "aprepro", "journal"],
        "batch_processing":    ["batch", "nographics"],
        "helix_coil":          ["helix", "coil"],
        # v4.7.0 operational knowledge (radia-mcp 0.32.0):
        "daemon_persistence":  ["daemon", "pid.lock", "attach", "0.01"],
        "license_warmup":      ["rlm", "warmup", "renewal", "3-day"],
        "utf8_path":           ["utf-8", "cp932", "japanese", "multibytetowidechar"],
        "v4_7_0":              ["4.7.0", "cubit", "release", "0.32.0"],
        "cubit_refresh":       ["warmup", "cubit_refresh", "rlm"],  # alias -> warmup
        "japanese_path":       ["utf-8", "cp932", "japanese"],     # alias -> utf8
    }

    def test_all_topic_is_huge(self):
        # 'all' concatenates every topic body — expect a very large string.
        body = get_cubit_documentation("all")
        _assert_substantial(body, min_chars=10000)
        # Should mention at least a handful of distinct domains.
        for markers in [["boolean"], ["helix", "coil"], ["aprepro"]]:
            _assert_contains_any(body, markers, label="cubit all")

    @pytest.mark.parametrize("topic,markers", list(TOPICS.items()))
    def test_subtopic(self, topic, markers):
        body = get_cubit_documentation(topic)
        _assert_substantial(body, min_chars=100)
        _assert_contains_any(body, markers, label=f"cubit.{topic}")

    def test_unknown_topic(self):
        _assert_unknown_topic_graceful(get_cubit_documentation)

    def test_solver_handoff_does_not_teach_retired_cubit_panels(self):
        body = get_cubit_documentation("ngsolve_panel")
        normalized = " ".join(body.split())

        assert "masked blocks in the Radia Simulink library" in normalized
        assert "Cubit is Radia's CAD, meshing, labeling" in normalized
        assert "QDialog" not in body
        assert "QProcess" not in body
        assert "Panel UI Settings" not in body


# ---------------------------------------------------------------------------
# export_knowledge (format-oriented)
# ---------------------------------------------------------------------------

class TestExportKnowledge:
    FORMATS = {
        "overview":       ["export", "format"],
        "gmsh_v4":        ["msh", "4.1", "gmsh"],
        "curved":         ["curved", "order", "netgen"],
        "netgen":         ["vol", "netgen"],
        "nastran":        ["nastran", "bdf"],
        "exodus":         ["exodus", "genesis", "block", "sideset"],
        "decision_guide": ["format", "choose", "decision", "recommend"],
    }

    def test_all_format(self):
        body = get_export_documentation("all")
        _assert_substantial(body, min_chars=500)
        _assert_contains_any(body, ["export", "gmsh", "netgen"],
                              label="export all")

    @pytest.mark.parametrize("fmt,markers", list(FORMATS.items()))
    def test_format(self, fmt, markers):
        body = get_export_documentation(fmt)
        _assert_substantial(body, min_chars=100)
        _assert_contains_any(body, markers, label=f"export.{fmt}")

    def test_unknown_format(self):
        _assert_unknown_topic_graceful(get_export_documentation)

    def test_gmsh_v4_documents_geo_as_radia_post_launch_artifact(self):
        body = get_export_documentation("gmsh_v4")

        assert ".msh v4.1 + .geo" in body
        assert "standard file to open for Radia post-processing is `.geo`" in body
        assert "`.msh` association is optional raw mesh/data inspection" in body
        assert ".geo` companion" in body
        assert "case.geo.opt" in body
        assert "case.msh.opt" in body
        assert "plain `case.opt` is not" in body


# ---------------------------------------------------------------------------
# mesh_diagnostics_knowledge
# ---------------------------------------------------------------------------

class TestDiagnosticsKnowledge:
    TOPICS = {
        "overview":      ["diagnos", "mesh"],
        "imprint_merge": ["imprint", "merge"],
        "exodus":        ["exodus", "block", "sideset", "nodeset"],
        "small_features":["small", "fillet", "cleanup", "blend"],
        "gaps_overlaps": ["gap", "overlap", "stitch"],
        "quality":       ["quality", "scaled jacobian", "skew"],
    }

    @pytest.mark.parametrize("topic,markers", list(TOPICS.items()))
    def test_topic(self, topic, markers):
        body = get_diagnostics_documentation(topic)
        _assert_substantial(body, min_chars=100)
        _assert_contains_any(body, markers, label=f"diagnostics.{topic}")

    def test_unknown_topic(self):
        _assert_unknown_topic_graceful(get_diagnostics_documentation)


# ---------------------------------------------------------------------------
# netgen_workflow_knowledge
# ---------------------------------------------------------------------------

class TestNetgenWorkflowKnowledge:
    WORKFLOWS = {
        "overview":        ["netgen", "workflow", "mesh"],
        "simple_cylinder": ["cylinder"],
        "simple_sphere":   ["sphere"],
        "simple_torus":    ["torus"],
        "accuracy":        ["order", "accuracy", "curving", "curve"],
        "kelvin_auto":     ["Kelvin", "auto"],
        "troubleshooting": ["troubleshoot", "fix", "error", "problem"],
    }

    @pytest.mark.parametrize("wf,markers", list(WORKFLOWS.items()))
    def test_workflow(self, wf, markers):
        body = get_netgen_documentation(wf)
        _assert_substantial(body, min_chars=100)
        _assert_contains_any(body, markers, label=f"netgen.{wf}")

    def test_unknown_workflow(self):
        _assert_unknown_topic_graceful(get_netgen_documentation)


# ---------------------------------------------------------------------------
# custom_toolbar_knowledge
# ---------------------------------------------------------------------------

class TestToolbarKnowledge:
    TOPICS = {
        "overview":             ["toolbar", "cubit"],
        "python_conventions":   ["python", "script", "shebang", "import"],
        "claro":                ["claro", "parent"],
        "llm_prompt":           ["llm", "prompt", "dialog"],
        "packaging":            ["tar", "distribution", "packaging"],
        "encapsulation":        ["workflow", "pattern", "tire"],
        "quality_convergence":  ["quality", "convergence", "refine"],
        "troubleshooting":      ["debug", "troubleshoot", "error", "problem"],
    }

    @pytest.mark.parametrize("topic,markers", list(TOPICS.items()))
    def test_topic(self, topic, markers):
        body = get_toolbar_documentation(topic)
        _assert_substantial(body, min_chars=100)
        _assert_contains_any(body, markers, label=f"toolbar.{topic}")

    def test_all_topic(self):
        body = get_toolbar_documentation("all")
        _assert_substantial(body, min_chars=500)
