"""Response regression tests for mcp-server-radia-ngsolve knowledge tools.

These tests assert that the canonical knowledge-serving functions
(``get_*_documentation``) return non-empty strings containing topic-
specific markers for every advertised topic, and degrade gracefully
on unknown topics.

They catch regressions that the lint-rule tests miss:
  * accidental empty topic body (e.g. variable rename on the docstring)
  * topic removed but still advertised in the tool's docstring
  * unknown-topic fallback becoming a hard error / KeyError / crash
  * "all" topic returning only the first entry

Run::

    pytest tests/mcp_server/test_radia_ngsolve_responses.py -v

Fast (all imports, no Cubit / NGSolve heavy work).
"""
from __future__ import annotations

import pytest

from radia_mcp.radia_ngsolve.knowledge.kelvin import get_kelvin_documentation
from radia_mcp.radia_ngsolve.knowledge.ngsbem_inductance import (
    get_ngsbem_inductance_documentation,
)
from radia_mcp.radia_ngsolve.knowledge.ngsolve import (
    get_ngsolve_documentation,
)
from radia_mcp.radia_ngsolve.knowledge.panel_gui_pitfalls import (
    get_panel_gui_pitfalls,
)
from radia_mcp.radia_ngsolve.knowledge.peec_inductance import (
    get_peec_inductance_documentation,
)
from radia_mcp.radia_ngsolve.knowledge.radia import get_radia_documentation
from radia_mcp.radia_ngsolve.knowledge.sparsesolv import (
    get_sparsesolv_documentation,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _assert_substantial(body: str, *, min_chars: int = 200) -> None:
    """A topic's body must be at least min_chars and have human-readable text."""
    assert isinstance(body, str), f"expected str, got {type(body).__name__}"
    assert len(body) >= min_chars, (
        f"topic body is too short ({len(body)} chars); "
        f"likely empty/truncated: {body[:100]!r}")


def _assert_contains_any(body: str, markers: list[str], label: str) -> None:
    """Body must contain at least one of the expected markers (case-insensitive)."""
    low = body.lower()
    hits = [m for m in markers if m.lower() in low]
    assert hits, (
        f"{label}: none of {markers} found in the {len(body)}-char body.\n"
        f"First 200 chars:\n{body[:200]!r}")


def _assert_unknown_topic_graceful(fn, *, unknown: str = "definitely_not_a_topic"):
    """Unknown topic must return a helpful string (not raise, not empty)."""
    r = fn(unknown)
    assert isinstance(r, str), f"{fn.__name__}({unknown!r}) returned {type(r).__name__}, want str"
    assert len(r) > 0, f"{fn.__name__}({unknown!r}) returned empty string"
    low = r.lower()
    assert "unknown" in low or "available" in low or unknown.lower() in low, (
        f"{fn.__name__}({unknown!r}) does not mention unknown/available/topic name: {r[:200]!r}")


# ---------------------------------------------------------------------------
# peec_inductance (the new 0.32.0 tool — extra-strict checks)
# ---------------------------------------------------------------------------

class TestPeecInductanceKnowledge:
    TOPICS = {
        "overview":          ["PEEC", "perimeter", "filament"],
        "centerline":        ["TORUS", "loft", "revolution", "spine"],
        "filament_dispatch": ["UV-map", "per-station", "n_peri", "Tier"],
        "japanese_path":     ["UTF-8", "CP_UTF8", "Windows"],
    }

    def test_all_topic_contains_every_subtopic(self):
        body = get_peec_inductance_documentation("all")
        _assert_substantial(body, min_chars=2000)
        # 'all' should concatenate every topic body
        for topic_markers in self.TOPICS.values():
            _assert_contains_any(body, topic_markers, label="peec_inductance all")

    @pytest.mark.parametrize("topic,markers", list(TOPICS.items()))
    def test_subtopic(self, topic, markers):
        body = get_peec_inductance_documentation(topic)
        _assert_substantial(body)
        _assert_contains_any(body, markers, label=f"peec_inductance.{topic}")

    def test_unknown_topic(self):
        _assert_unknown_topic_graceful(get_peec_inductance_documentation)


# ---------------------------------------------------------------------------
# kelvin
# ---------------------------------------------------------------------------

class TestKelvinKnowledge:
    def test_all(self):
        body = get_kelvin_documentation("all")
        _assert_substantial(body, min_chars=1000)
        _assert_contains_any(body, ["Kelvin", "open boundary", "sphere"],
                              label="kelvin all")

    def test_unknown_topic(self):
        _assert_unknown_topic_graceful(get_kelvin_documentation)


# ---------------------------------------------------------------------------
# ngsbem_inductance
# ---------------------------------------------------------------------------

class TestNgsbemInductanceKnowledge:
    def test_all(self):
        body = get_ngsbem_inductance_documentation("all")
        _assert_substantial(body, min_chars=1000)
        _assert_contains_any(body, ["BEM", "inductance", "LaplaceSL", "ngsolve"],
                              label="ngsbem all")

    def test_unknown_topic(self):
        _assert_unknown_topic_graceful(get_ngsbem_inductance_documentation)


# ---------------------------------------------------------------------------
# ngsolve
# ---------------------------------------------------------------------------

class TestNgsolveKnowledge:
    def test_all(self):
        body = get_ngsolve_documentation("all")
        _assert_substantial(body, min_chars=1000)
        _assert_contains_any(body, ["NGSolve", "FESpace", "HCurl", "BilinearForm"],
                              label="ngsolve all")

    def test_unknown_topic(self):
        _assert_unknown_topic_graceful(get_ngsolve_documentation)

    def test_ngsolve_bem_50_vol_visualization_topic(self):
        body = get_ngsolve_documentation("ngsolve_bem_50")
        _assert_substantial(body, min_chars=1000)
        _assert_contains_any(body, ["Netgen", ".vol", "50-case", "NGSolve.BEM"],
                             label="ngsolve ngsolve_bem_50")
        assert "Do not call a cross-code artifact \"analytic\"" in body
        assert "radia-vol-viewer --register" in body
        assert "Ng_LoadMesh" in body
        assert "mesh-free" in body
        assert "GYP-001..010" in body
        assert "struck drum" in body

    def test_vol_double_click_alias(self):
        body = get_ngsolve_documentation("vol_double_click")
        _assert_substantial(body, min_chars=1000)
        assert "raw Windows file association" in body
        assert "solution visual mode" in body

    def test_vibroacoustic_drum_alias(self):
        body = get_ngsolve_documentation("vibroacoustic_drum")
        _assert_substantial(body, min_chars=1000)
        assert "baffled" in body
        assert "normal velocity" in body
        assert "radiation impedance" in body
        assert "MATLAB rather than Gmsh" in body
        assert "animated GIF" in body
        assert "headless-friendly" in body
        assert "high-order impedance" in body
        assert "never as a Kelvin" in body
        assert "Kelvin boundary language is not allowed" in body
        assert "struck top" in body
        assert "axis-equal" in body
        assert "not a hemisphere" in body
        assert "real-drum FEM/BEM" in body
        assert "lower-half radiation" in body
        assert "reduced FEM" in body
        assert "ode45" in body
        assert "retarded boundary" in body
        assert "must not split the observation" in body
        assert "every exterior air observation point" in body
        assert "direction-only painting is a" in body
        assert "color-map only propagating air" in body
        assert "3D axisymmetric" in body
        assert "r-z slice" in body
        assert "not yet a full 3D" in body
        assert "volFemBemIfftResponse" in body
        assert "frequency-domain Helmholtz FEM/BEM" in body
        assert "inverse FFT" in body
        assert "not a periodic sine-wave animation" in body
        assert "not convolution-quadrature TD-BEM" in body
        assert "volTdBemConvolutionQuadrature" in body
        assert "BDF generating function" in body
        assert "Laplace-domain retarded single-layer" in body
        assert "real Lubich CQ" in body
        assert "volFemBemCoupledConvolutionQuadrature" in body
        assert "H1/P1" in body
        assert "interior wave FEM" in body
        assert "(1/2 Mb-K(s))*T" in body
        assert "-S(s)q + D(s)Tu" in body
        assert "Calderon/Johnson-Nedelec coupled CQ" in body
        assert "retarded double-layer `K(s)`" in body
        assert "SingleLayerTeaching" in body
        assert "one-sided baffled top-head" in body
        assert "r-z pressure snapshots" in body

    def test_curved_vol_geometry_alias(self):
        body = get_ngsolve_documentation("curved_vol_geometry")
        _assert_substantial(body, min_chars=1000)
        assert "Curve-only high-order geometry" in body
        assert "CurvedBoundaryView" in body
        assert "EnableCurvedGeometry=true" in body


# ---------------------------------------------------------------------------
# radia
# ---------------------------------------------------------------------------

class TestRadiaKnowledge:
    def test_all(self):
        body = get_radia_documentation("all")
        _assert_substantial(body, min_chars=1000)
        _assert_contains_any(body, ["Radia", "ObjHexahedron", "magnetization", "MMM"],
                              label="radia all")

    def test_unknown_topic(self):
        _assert_unknown_topic_graceful(get_radia_documentation)


# ---------------------------------------------------------------------------
# sparsesolv
# ---------------------------------------------------------------------------

class TestSparsesolvKnowledge:
    def test_all(self):
        body = get_sparsesolv_documentation("all")
        _assert_substantial(body, min_chars=500)
        _assert_contains_any(body, ["sparsesolv", "CompactAMS", "AMS", "HX",
                                     "preconditioner"],
                              label="sparsesolv all")


# ---------------------------------------------------------------------------
# panel_gui_pitfalls
# ---------------------------------------------------------------------------

class TestPanelGuiPitfalls:
    def test_full_document(self):
        # Convention for this tool: empty string returns the full document.
        body = get_panel_gui_pitfalls("")
        _assert_substantial(body, min_chars=500)
        _assert_contains_any(body, ["panel", "combo", "widget", "GUI", "pitfall"],
                              label="panel_gui_pitfalls full")

    def test_unknown_topic(self):
        _assert_unknown_topic_graceful(get_panel_gui_pitfalls)
