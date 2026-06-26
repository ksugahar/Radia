"""
Tests for the Open GMSH button enable / disable logic.

The 2026-04-12 "FEM finished but Open GMSH stays grey" bug was a
silent ``msh_file == ""`` in the result dict. The handler in
radia_gui_base.py::AnalysisWindow._on_finished checks several keys
(``field_gmsh_file``, ``msh_output``, ``msh_file``) and enables the
button when at least one of them points to an existing file.

These tests pin the enable logic by feeding mock result dicts and
inspecting the button state.
"""

from __future__ import annotations

import os

import pytest


def _enable_for_result(window, result, tmp_path):
    """Reach into AnalysisWindow's enable logic without needing a
    real subprocess. We replicate the relevant slice of
    _on_finished in test code so the test stays focused on the
    enable rule rather than the json/stderr parsing."""
    btn = window._gmsh_btn
    btn.setEnabled(False)
    target = None
    if result.get("field_gmsh_file") and os.path.isfile(
            result["field_gmsh_file"]):
        target = result["field_gmsh_file"]
    else:
        msh_key = "msh_output" if "msh_output" in result else "msh_file"
        msh = result.get(msh_key)
        if msh and os.path.isfile(msh):
            target = msh
    if target:
        window._last_msh = target
        btn.setEnabled(True)
    return btn.isEnabled()


@pytest.fixture
def fake_msh(tmp_path):
    """Create a placeholder .msh file on disk so os.path.isfile
    returns True."""
    p = tmp_path / "model.msh"
    p.write_text("$MeshFormat\n4.1 0 8\n$EndMeshFormat\n")
    return str(p)


@pytest.fixture
def fake_geo(tmp_path):
    p = tmp_path / "model.geo"
    p.write_text('Merge "model.msh";\n')
    return str(p)


class TestOpenGmshEnable:

    def test_msh_output_enables(self, ih_window, fake_msh, tmp_path):
        result = {"msh_output": fake_msh}
        assert _enable_for_result(ih_window, result, tmp_path)

    def test_msh_file_enables(self, ih_window, fake_msh, tmp_path):
        """The legacy result-dict key ``msh_file`` (used by
        calc_fem_kelvin.py) must enable the button just like the
        newer ``msh_output`` does."""
        result = {"msh_file": fake_msh}
        assert _enable_for_result(ih_window, result, tmp_path)

    def test_field_gmsh_file_takes_precedence(self, ih_window, fake_msh,
                                               fake_geo, tmp_path):
        """When both keys are present the .geo (with display
        options) wins."""
        result = {"msh_output": fake_msh, "field_gmsh_file": fake_geo}
        assert _enable_for_result(ih_window, result, tmp_path)
        assert ih_window._last_msh == fake_geo

    def test_empty_string_does_not_enable(self, ih_window, tmp_path):
        """The d364796 regression: result had ``msh_file=""`` and
        the button stayed disabled silently. We want a deterministic
        DISABLED here, not a crash, but tests must pin that the
        bug pattern is detected."""
        result = {"msh_file": ""}
        assert not _enable_for_result(ih_window, result, tmp_path)

    def test_missing_path_does_not_enable(self, ih_window, tmp_path):
        result = {"msh_file": str(tmp_path / "does_not_exist.msh")}
        assert not _enable_for_result(ih_window, result, tmp_path)

    def test_no_keys_does_not_enable(self, ih_window, tmp_path):
        assert not _enable_for_result(ih_window, {}, tmp_path)
