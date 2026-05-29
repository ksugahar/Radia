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


# ---- Panel Layout Policy (functional): short window scrolls, never compress ----

def test_short_window_engages_scrollbar_no_compression(ih_window):
    """At a deliberately short window height, the QScrollArea MUST engage
    its vertical scrollbar (``maximum() > 0``) rather than compress the
    form into the viewport.  A compressed form (``maximum() == 0`` with
    the inner widget shrunk below its sizeHint) is the structural cause
    of the vertical text clipping bug (kubota 2026-05-29).

    Functional replacement for a brittle screenshot regression: image
    diff would catch the clipping visually but depends on fonts / DPI /
    OS.  Here we assert the scroll-area invariant directly."""
    from PySide6.QtCore import QCoreApplication

    scrolls = [s for s in ih_window.findChildren(QScrollArea)
               if s.widgetResizable()]
    assert scrolls, "expected at least one widgetResizable QScrollArea"
    panel_scroll = scrolls[0]
    inner = panel_scroll.widget()
    assert inner is not None, "panel scroll area has no inner widget"

    # Resize to a deliberately short height (any IH panel form's natural
    # sizeHint exceeds this).
    ih_window.resize(820, 260)
    ih_window.show()
    QCoreApplication.processEvents()
    QCoreApplication.processEvents()

    inner_hint = inner.sizeHint().height()
    viewport_h = panel_scroll.viewport().height()
    vbar = panel_scroll.verticalScrollBar()

    # Sanity: the test must actually exercise the "form too tall" case.
    assert inner_hint > viewport_h, (
        f"form sizeHint h={inner_hint} does not exceed viewport "
        f"h={viewport_h}; resize did not take effect, test is trivial")

    # Core invariant: vertical scrollbar has scroll range (maximum > 0)
    # <=> the inner widget is taller than the viewport <=> rows kept
    # their natural height <=> NOT compressed.
    assert vbar.maximum() > 0, (
        f"vertical scrollbar maximum == 0: the form was COMPRESSED "
        f"into the viewport (inner.height={inner.height()}, "
        f"sizeHint={inner_hint}, viewport={viewport_h}). Field text "
        f"would clip vertically.")
    # And the inner widget must keep at least its sizeHint height
    # (widgetResizable grows it to sizeHint when viewport is shorter --
    # the "never compress" invariant from CLAUDE.md Panel Layout Policy).
    assert inner.height() >= inner_hint, (
        f"inner widget compressed: height={inner.height()} < "
        f"sizeHint={inner_hint}")
