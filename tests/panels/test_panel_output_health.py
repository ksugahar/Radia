"""
test_panel_output_health.py -- panel "health" checks BEYOND widget
wiring: result-output completeness + layout robustness.

Encodes the kubota 2026-05-29 reports that the widget-wiring tests
(test_ih_panel_qt.py etc.) did NOT catch:

  1. OUTPUT showed no element count / DoF / compute-time / heat /
     temperature -- because `_on_finished`'s summary checked wrong key
     names (`n_dofs`/`t_solve`/`P_total_W`) while calc_*.py emit
     `wp_ndof`/`t_bem_solve_s`/`P_wp_W`.  `_append_standard_summary`
     now surfaces them generically (Result Output Policy).
  2. A short window compressed QFormLayout rows -> vertical text
     clipping ("文字潰れ").  AnalysisWindow now wraps the form in a
     QScrollArea (vertical scrollbar, NEVER compress) at 10pt
     (Panel Layout Policy).

LAB caveat: importing PySide6 under pytest crashes with
``0xc0000139 / DLL load failed`` (conftest MKL add_dll_directory).  CI
``--ignore=tests/panels`` so this does not run there either.  On LAB,
validate with the standalone offscreen smoke instead -- see the
``panel-qt-test`` / ``panel-wheel-guard`` skills.  This file is the
clean-env / future-CI gate.
"""

from __future__ import annotations

import os

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PySide6.QtWidgets import QScrollArea


# ---- Panel Layout Policy: 10pt + vertical scrollbar, never compress ----

def test_panel_base_font_is_10pt():
    import radia_gui_base as gb
    assert gb.PANEL_BASE_FONT_POINT_SIZE == 10, (
        "Panel Layout Policy: base panel font must be 10pt")


def test_form_is_in_resizable_scrollarea(ih_window):
    """The parameter form is hosted in a QScrollArea(setWidgetResizable)
    so a short window scrolls instead of compressing rows (which clips
    the field text and makes entered values unconfirmable)."""
    scrolls = ih_window.findChildren(QScrollArea)
    assert scrolls, "AnalysisWindow must wrap the form in a QScrollArea"
    assert any(s.widgetResizable() for s in scrolls), (
        "the panel QScrollArea must be widgetResizable so rows keep "
        "their natural (uncompressed) sizeHint height")


# ---- Result Output Policy: ne / DoF / time + integral quantities ----

def test_output_summary_surfaces_ne_dof_time_heat(ih_window):
    """_append_standard_summary renders element count, DoF, the
    compute-time breakdown and heat for a BEM-A-shaped result --
    keyed on the ACTUAL emitted names, not a fixed cascade."""
    mock = {"P_wp_W": 1.23e-3, "wp_ndof": 5040, "wp_mesh_n_tris": 2520,
            "wp_mesh_nv": 1262, "t_bem_assembly_s": 8.4,
            "t_bem_solve_s": 3.2, "t_coil_solve_s": 1.1}
    ih_window._output.clear()
    ih_window._append_standard_summary(mock)
    txt = ih_window._output.toPlainText()
    assert "Elements = 2520" in txt, txt
    assert "DoF = 5040" in txt, txt
    assert "Heat (P)" in txt, txt
    assert "Compute time" in txt and "t_bem_solve_s" in txt, txt


def test_output_summary_temperature_mean_max_min(ih_window):
    """Thermal: temperature reported as mean (volume-averaged) / max /
    min -- not a single peak value."""
    mock = {"T_mean_C": 85.3, "T_max_C": 142.7, "T_min_C": 23.9,
            "ndof": 31021, "ne": 18044, "t_total_s": 12.4}
    ih_window._output.clear()
    ih_window._append_standard_summary(mock)
    txt = ih_window._output.toPlainText()
    assert "Temperature [C]" in txt, txt
    assert "mean 85.3" in txt and "max 142.7" in txt and "min 23.9" in txt, txt
