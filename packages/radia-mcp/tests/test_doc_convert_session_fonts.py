"""doc_convert_session_font_check / _repair: the session font table, not the file."""
import sys

import pytest

from radia_mcp.doc_convert.plans.T17_session_fonts import (
    _font_files,
    doc_convert_session_font_check,
    doc_convert_session_font_repair,
    probe_face,
)

windows_only = pytest.mark.skipif(sys.platform != "win32", reason="reads the Windows session font table")


@windows_only
def test_probe_reports_a_healthy_shipped_face():
    r = probe_face("Arial")
    assert r["realised_as"].casefold() == "arial"
    assert r["font_data_bytes"] and r["font_data_bytes"] > 100_000
    assert r["glyph_outline_bytes"] and r["glyph_outline_bytes"] > 0
    assert r["ok"] and not r["broken"]


@windows_only
def test_check_lists_broken_and_substituted_faces():
    r = doc_convert_session_font_check("Arial, No Such Face Zz")
    faces = {x["face"]: x for x in r["faces"]}
    assert faces["Arial"]["ok"]
    assert faces["No Such Face Zz"]["substituted"]
    assert "No Such Face Zz" in r["substituted"]
    assert "No Such Face Zz" not in r["broken"]


@windows_only
def test_font_files_are_resolved_for_known_and_registry_faces():
    assert any(p.name.lower() == "times.ttf" for p in _font_files("Times New Roman"))
    assert any(p.name.lower().startswith("georgia") for p in _font_files("Georgia"))


@windows_only
def test_repair_leaves_a_healthy_face_alone_by_default():
    r = doc_convert_session_font_repair("Arial")
    assert r["ok"]
    assert r["results"][0]["action"] == "healthy, left alone"


def test_non_windows_is_refused():
    if sys.platform == "win32":
        pytest.skip("Windows host")
    assert "error" in doc_convert_session_font_check()
    assert "error" in doc_convert_session_font_repair()
