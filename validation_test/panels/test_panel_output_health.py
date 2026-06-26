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
``--ignore=validation_test/panels`` so this does not run there either.  On LAB,
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


def test_output_summary_shows_radia_version(ih_window):
    """The Result block surfaces the radia version that produced the
    numbers (calc_main stamps result['radia_version']; the panel renders
    it) -- so the persisted .log / screenshot records which radia ran.
    Absent the key (older calc), no version line is shown (no crash)."""
    ih_window._output.clear()
    ih_window._append_standard_summary({"radia_version": "4.89.1",
                                        "wp_ndof": 100})
    txt = ih_window._output.toPlainText()
    assert "Radia version = 4.89.1" in txt, txt
    # Gracefully absent when the calc did not stamp it.
    ih_window._output.clear()
    ih_window._append_standard_summary({"wp_ndof": 100})
    assert "Radia version" not in ih_window._output.toPlainText()


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


# ---- Result Output Persistence Policy (2026-05-30) -------------------
# Every Run leaves a <base>_<suffix>.log next to the --output JSON,
# overwritten on each Run, containing a verbatim copy of the Output
# window.


def test_persist_output_log_writes_log_next_to_json(ih_window, tmp_path):
    """_persist_output_log writes the Output text to a .log file
    derived from _last_output_json by swapping .json -> .log."""
    json_path = str(tmp_path / "model_fem_kelvin.json")
    ih_window._last_output_json = json_path
    ih_window._output.clear()
    ih_window._output.appendPlainText("=== Launching ===")
    ih_window._output.appendPlainText("--- Result ---")
    ih_window._output.appendPlainText("  DoF = 31021")
    ih_window._persist_output_log()

    log_path = str(tmp_path / "model_fem_kelvin.log")
    assert os.path.isfile(log_path), (
        f"persist_output_log did not write {log_path}")
    with open(log_path, encoding="utf-8") as f:
        body = f.read()
    assert "=== Launching ===" in body
    assert "--- Result ---" in body
    assert "DoF = 31021" in body


def test_persist_output_log_overwrites_previous_run(ih_window, tmp_path):
    """Per user-confirmed policy decision (2026-05-30): the .log is
    OVERWRITTEN on each Run, not rotated or appended. Calling
    _persist_output_log twice must leave only the latest content."""
    json_path = str(tmp_path / "model_peec_bem.json")
    ih_window._last_output_json = json_path
    log_path = str(tmp_path / "model_peec_bem.log")

    # First Run
    ih_window._output.clear()
    ih_window._output.appendPlainText("OLD_RUN_MARKER")
    ih_window._persist_output_log()
    with open(log_path, encoding="utf-8") as f:
        assert "OLD_RUN_MARKER" in f.read()

    # Second Run
    ih_window._output.clear()
    ih_window._output.appendPlainText("NEW_RUN_MARKER")
    ih_window._persist_output_log()

    with open(log_path, encoding="utf-8") as f:
        body = f.read()
    assert "NEW_RUN_MARKER" in body, "second Run did not write its content"
    assert "OLD_RUN_MARKER" not in body, (
        "overwrite invariant violated -- old Run content survived")


def test_persist_output_log_silent_when_no_output_json(ih_window, tmp_path):
    """No --output JSON means no log path to derive from. The helper
    must be a no-op (no exception, no spurious file). This protects
    panels that legitimately don't pass --output (e.g. some debug
    paths)."""
    ih_window._last_output_json = None
    ih_window._output.clear()
    ih_window._output.appendPlainText("anything")
    # Must NOT raise.
    ih_window._persist_output_log()
    # And must not have written into tmp_path.
    assert not list(tmp_path.iterdir()), (
        f"persist_output_log wrote something despite no JSON path: "
        f"{list(tmp_path.iterdir())}")


def test_persist_output_log_captures_failed_run(ih_window, tmp_path):
    """The error-return path of _on_finished also persists the log,
    so a failed Run leaves an audit trail (the failure text is in the
    Output box at that point)."""
    json_path = str(tmp_path / "model_omega_reduced.json")
    ih_window._last_output_json = json_path
    ih_window._output.clear()
    ih_window._output.appendPlainText("=== Launching ===")
    ih_window._output.appendPlainText("Traceback (most recent call last):")
    ih_window._output.appendPlainText("  RuntimeError: solver diverged")
    ih_window._output.appendPlainText("*** Process exited with code 1")
    ih_window._persist_output_log()

    log_path = str(tmp_path / "model_omega_reduced.log")
    with open(log_path, encoding="utf-8") as f:
        body = f.read()
    assert "Traceback" in body
    assert "solver diverged" in body
    assert "Process exited with code 1" in body
