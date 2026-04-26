"""Panel quality-assurance checks.

A single entry point ``check_panel_health(window)`` runs a suite of
layout / accessibility / correctness checks on a Radia PySide6 panel
window.  Used by:

  - the ``/deploy`` skill's Step 8 (pre-deploy visual check)
  - the ``pytest tests/panels/test_panel_qa.py`` regression suite

Why a dedicated module instead of inline markdown heredoc:

  - Reusable from pytest + skill
  - Each check has a clear name → failures point at a specific cause
  - Checks can be extended without touching SKILL.md

All checks return a ``CheckResult(ok, detail)``.  ``run_panel_checks``
aggregates them, saves a screenshot, and raises ``PanelQAError`` if
any mandatory check fails.
"""
from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path

# Safe to import at module load — Radia panels all use PySide6.
from PySide6.QtCore import QPoint
from PySide6.QtGui import QPixmap
from PySide6.QtWidgets import (QComboBox, QFormLayout, QLabel, QLineEdit,
                                 QMainWindow, QPushButton, QWidget)


# 2K (2560x1440) minimum display, usable ~1350 px vertical, ~2400 horizontal.
# Thresholds reflect the 11pt baseline font (Segoe UI); the 9pt-era
# numbers (1100 / 900) under-budget by ~20% for the larger glyphs.
MAX_HEIGHT_RED = 1300     # hard fail: bottom buttons off-screen on 2K
MAX_HEIGHT_YELLOW = 1100  # warning: getting tall
MAX_WIDTH_RED = 1300      # hard fail: panel eats too much horizontal real estate
MAX_WIDTH_YELLOW = 1100   # warning: method combo text probably long

# Minimum readable font size for the 2K target display.  Qt's OS default
# of 9pt is unreadable on 2K+ at 100% scaling, so the lab baseline is
# 11pt (set via apply_panel_base_font in radia_gui_base).  10pt leaves
# 1pt of headroom for intentionally smaller status-line text.
MIN_FONT_POINT_SIZE = 10

# Unicode-math characters the console / cp932 cannot render.  If these
# appear in any widget label / tooltip, they fail the ASCII check.
# (Policy: CLAUDE.md "Windows Console Encoding (cp932)".)
FORBIDDEN_UNICODE = ("±", "²", "³", "µ", "·", "×", "÷",
                     "→", "←", "↑", "↓", "≤", "≥", "≠")


class PanelQAError(AssertionError):
    """Raised by run_panel_checks when any mandatory check fails."""


@dataclass
class CheckResult:
    name: str
    ok: bool
    detail: str = ""
    severity: str = "error"  # "error" or "warning"

    def __str__(self):
        flag = "OK  " if self.ok else ("WARN" if self.severity == "warning"
                                         else "FAIL")
        d = f"  -- {self.detail}" if self.detail else ""
        return f"  [{flag}] {self.name}{d}"


# ----------------------------------------------------------------------
# Individual checks
# ----------------------------------------------------------------------

def check_height(window) -> CheckResult:
    h = window.sizeHint().height()
    if h >= MAX_HEIGHT_RED:
        return CheckResult("height", False,
                            f"{h} px >= {MAX_HEIGHT_RED} (2K unusable)")
    if h >= MAX_HEIGHT_YELLOW:
        return CheckResult("height", True,
                            f"{h} px (getting tall; consider trimming)",
                            severity="warning")
    return CheckResult("height", True, f"{h} px")


def check_width(window) -> CheckResult:
    w = window.sizeHint().width()
    if w >= MAX_WIDTH_RED:
        return CheckResult("width", False,
                            f"{w} px >= {MAX_WIDTH_RED} (too wide for 2K)")
    if w >= MAX_WIDTH_YELLOW:
        return CheckResult("width", True,
                            f"{w} px (method combo text probably long)",
                            severity="warning")
    return CheckResult("width", True, f"{w} px")


def check_buttons_reachable(window) -> CheckResult:
    """Run button's bottom must be within the window's sizeHint height.

    Does NOT mutate window size — previous implementation called
    ``window.resize()`` which inflated subsequent sizeHint() calls
    for other checks.  Use ``geometry()`` (which is valid after
    ``adjustSize()`` ran in ``check_panel_health``) to read the
    button's y position.
    """
    btns = [b for b in window.findChildren(QPushButton)
            if "Run" in b.text().strip()]
    if not btns:
        return CheckResult("buttons_reachable", False,
                            "no Run button found in window")
    run_btn = btns[0]
    g = run_btn.geometry()  # local to button's parent
    # Map bottom-left of the button to window coordinates
    btn_bottom = run_btn.mapTo(window, QPoint(0, g.height())).y()
    win_h = window.sizeHint().height()
    if btn_bottom > win_h + 2:  # 2 px tolerance for frame/padding
        return CheckResult("buttons_reachable", False,
                            f"Run bottom y={btn_bottom} > sizeHint h={win_h}")
    return CheckResult("buttons_reachable", True,
                        f"Run y={btn_bottom} of {win_h}")


def _iter_form_layouts(root):
    for w in [root] + root.findChildren(QWidget):
        lay = w.layout() if hasattr(w, "layout") else None
        if isinstance(lay, QFormLayout):
            yield lay


def _is_section_header_label(field_widget) -> bool:
    """A section header is a QLabel whose text starts with ``<b>``."""
    if not isinstance(field_widget, QLabel):
        return False
    return field_widget.text().lstrip().startswith("<b>")


def _row_has_visible_input(lay, row):
    """A row contributes to 'visible content' if either its label column
    contains a non-empty string OR its field column is a visible input
    widget (not another section header label)."""
    label_item = lay.itemAt(row, QFormLayout.LabelRole)
    field_item = lay.itemAt(row, QFormLayout.FieldRole)
    field_w = field_item.widget() if field_item else None
    if field_w is None:
        return False
    if not field_w.isVisibleTo(field_w.window()):
        return False
    if _is_section_header_label(field_w):
        return False
    return True


def check_no_orphan_section_headers(window) -> CheckResult:
    """A section header with no visible content rows below it (before
    the next section or end-of-form) is 'orphan' — cosmetic garbage
    that confuses users."""
    orphans = []
    for lay in _iter_form_layouts(window):
        n = lay.rowCount()
        for i in range(n):
            field_item = lay.itemAt(i, QFormLayout.FieldRole)
            field_w = field_item.widget() if field_item else None
            if not _is_section_header_label(field_w):
                continue
            # Is this section header itself hidden? Then fine.
            if not field_w.isVisibleTo(window):
                continue
            # Scan rows below this one until next section header (or end)
            has_content = False
            for j in range(i + 1, n):
                field_j = lay.itemAt(j, QFormLayout.FieldRole)
                fj = field_j.widget() if field_j else None
                if _is_section_header_label(fj):
                    break
                if _row_has_visible_input(lay, j):
                    has_content = True
                    break
            if not has_content:
                header_text = field_w.text().strip()
                orphans.append(header_text)
    if orphans:
        return CheckResult("no_orphan_section_headers", False,
                            f"orphan headers: {orphans}")
    return CheckResult("no_orphan_section_headers", True)


def check_ascii_only_labels(window) -> CheckResult:
    """Panel UI text must be English / ASCII.  CLAUDE.md policy: no
    Unicode math symbols (cp932 can't render them in the subprocess log)."""
    bad = []
    # QLabel text
    for lbl in window.findChildren(QLabel):
        txt = lbl.text()
        for ch in FORBIDDEN_UNICODE:
            if ch in txt:
                bad.append(("QLabel", txt, ch))
    # QComboBox items
    for combo in window.findChildren(QComboBox):
        for i in range(combo.count()):
            it = combo.itemText(i)
            for ch in FORBIDDEN_UNICODE:
                if ch in it:
                    bad.append(("QComboBox", it, ch))
    # QPushButton text
    for btn in window.findChildren(QPushButton):
        txt = btn.text()
        for ch in FORBIDDEN_UNICODE:
            if ch in txt:
                bad.append(("QPushButton", txt, ch))
    # QLineEdit placeholder
    for le in window.findChildren(QLineEdit):
        ph = le.placeholderText()
        for ch in FORBIDDEN_UNICODE:
            if ch in ph:
                bad.append(("QLineEdit.placeholder", ph, ch))
    if bad:
        return CheckResult("ascii_only_labels", False,
                            f"{len(bad)} Unicode-math char(s) in UI text "
                            f"(e.g. {bad[0]})")
    return CheckResult("ascii_only_labels", True)


def check_method_combo_populated(window) -> CheckResult:
    """If the panel has a Method combo, it must have >=1 items AND
    have something selected (not a blank currentText)."""
    panel = getattr(window, "_panel", None)
    if panel is None:
        return CheckResult("method_combo_populated", True,
                            "(no panel — window-only)")
    combo = getattr(panel, "_method_combo", None)
    if combo is None:
        # Check for a generic "method" widget key
        combo = panel._widgets.get("method") if hasattr(panel, "_widgets") else None
    if combo is None:
        return CheckResult("method_combo_populated", True,
                            "(no method combo)")
    n = combo.count()
    cur = combo.currentText().strip()
    if n == 0:
        return CheckResult("method_combo_populated", False,
                            "Method combo is empty")
    if not cur:
        return CheckResult("method_combo_populated", False,
                            f"Method combo has {n} items but blank current")
    return CheckResult("method_combo_populated", True,
                        f"{n} items, current={cur!r}")


def check_font_size_min(window) -> CheckResult:
    """Every visible label / input / button must render at >= MIN_FONT_POINT_SIZE.

    Qt's OS default font (Windows: 9pt Segoe UI) is unreadable on 2K+
    displays at 100% scaling.  ``apply_panel_base_font`` bumps the
    QApplication baseline to 11pt; this check guards against any
    widget regressing back to the 9pt default via an explicit
    ``setFont`` / ``font-size:`` style override.

    Note: ``QFont.pointSize()`` returns -1 if the font was set in
    pixels (``setPixelSize`` / ``font-size: Npx``).  In that case we
    convert via ``QFontInfo`` which always returns a positive value.
    """
    from PySide6.QtGui import QFontInfo
    too_small = []
    widget_types = (QLabel, QLineEdit, QComboBox, QPushButton)
    # Also include QPlainTextEdit / QSpinBox / QGroupBox if present
    from PySide6.QtWidgets import QPlainTextEdit, QSpinBox, QGroupBox
    widget_types = widget_types + (QPlainTextEdit, QSpinBox, QGroupBox)

    for w in window.findChildren(QWidget):
        if not isinstance(w, widget_types):
            continue
        if not w.isVisibleTo(window):
            continue
        f = w.font()
        pt = f.pointSize()
        if pt <= 0:
            # Set via setPixelSize / font-size: Npx — use QFontInfo
            pt = QFontInfo(f).pointSize()
        if pt > 0 and pt < MIN_FONT_POINT_SIZE:
            label = (w.text()[:30] if hasattr(w, "text") and w.text()
                     else type(w).__name__)
            too_small.append((label, pt))
    if too_small:
        return CheckResult("font_size_min", False,
                            f"{len(too_small)} widget(s) below "
                            f"{MIN_FONT_POINT_SIZE}pt: {too_small[:3]}")
    return CheckResult("font_size_min", True,
                        f"all >= {MIN_FONT_POINT_SIZE}pt")


def check_visible_rows_have_labels(window) -> CheckResult:
    """Every visible input row must have a non-empty label in the
    label column.  Empty labels make the form ambiguous."""
    nameless = []
    for lay in _iter_form_layouts(window):
        n = lay.rowCount()
        for i in range(n):
            field_item = lay.itemAt(i, QFormLayout.FieldRole)
            fw = field_item.widget() if field_item else None
            if fw is None or not fw.isVisibleTo(window):
                continue
            if _is_section_header_label(fw):
                continue
            # Skip plain-text QLabel fields (status lines, footers)
            if isinstance(fw, QLabel):
                continue
            label_item = lay.itemAt(i, QFormLayout.LabelRole)
            lw = label_item.widget() if label_item else None
            label_text = lw.text().strip() if isinstance(lw, QLabel) else ""
            if not label_text:
                nameless.append(type(fw).__name__)
    if nameless:
        return CheckResult("visible_rows_have_labels", False,
                            f"{len(nameless)} unlabeled input row(s): "
                            f"{nameless[:5]}")
    return CheckResult("visible_rows_have_labels", True)


# ----------------------------------------------------------------------
# Runner
# ----------------------------------------------------------------------

DEFAULT_CHECKS = [
    check_height,
    check_width,
    check_buttons_reachable,
    check_no_orphan_section_headers,
    check_ascii_only_labels,
    check_method_combo_populated,
    check_visible_rows_have_labels,
    check_font_size_min,
]


def check_panel_health(window, checks=None):
    """Run all checks against a panel window, return list of CheckResult."""
    checks = checks or DEFAULT_CHECKS
    # Settle layout ONCE so geometry() / sizeHint() return realistic
    # values.  show() is required under the offscreen platform plugin
    # for child widget geometries to be computed; no actual window
    # appears on screen.  Checks must NOT mutate size after this.
    from PySide6.QtWidgets import QApplication
    window.show()
    QApplication.processEvents()
    window.adjustSize()
    QApplication.processEvents()
    results = [c(window) for c in checks]
    window.hide()
    return results


def grab_screenshot(window, path):
    """Save a PNG snapshot of the window at its natural sizeHint."""
    sh = window.sizeHint()
    window.resize(sh.width(), sh.height())
    from PySide6.QtWidgets import QApplication
    QApplication.processEvents()
    pix: QPixmap = window.grab()
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    pix.save(str(path))


def run_panel_checks(tag, window, screenshot_dir, strict=True):
    """Check a single panel and save a screenshot.

    Args:
        tag: short identifier (used in filename and log)
        window: the panel window (QMainWindow or QWidget)
        screenshot_dir: where to save the PNG
        strict: if True (default), raise PanelQAError on any FAIL result.

    Returns:
        (results, png_path)
    """
    results = check_panel_health(window)
    png = Path(screenshot_dir) / f"panel_{tag}.png"
    grab_screenshot(window, png)
    failed = [r for r in results if not r.ok and r.severity == "error"]
    print(f"=== {tag} ===")
    print(f"  screenshot: {png}")
    for r in results:
        print(str(r))
    if failed and strict:
        raise PanelQAError(
            f"{tag}: {len(failed)} mandatory check(s) failed: "
            f"{[r.name for r in failed]}")
    return results, str(png)


# ----------------------------------------------------------------------
# Registry of panels to render
# ----------------------------------------------------------------------

def get_panel_registry():
    """Return a list of (tag, WindowClass, combo_attr_or_key, value).

    Each tuple describes one render of one panel in a specific mode.
    Used by the deploy skill and pytest to iterate every supported
    panel × mode combination.
    """
    import radia_ih, radia_em, radia_pcb
    return [
        # IH — 3 methods
        ("ih_ind", radia_ih.IHWindow,
         "_method_combo", radia_ih.METHOD_PEEC_IND),
        ("ih_bem", radia_ih.IHWindow,
         "_method_combo", radia_ih.METHOD_PEEC_BEM),
        ("ih_fem", radia_ih.IHWindow,
         "_method_combo", radia_ih.METHOD_FEM_FULL),
        # EM — formulations (Omega / A-Phi / MSC / Kelvin Benchmark)
        ("em_omega", radia_em.EMWindow, "formulation", "Omega"),
        ("em_aphi", radia_em.EMWindow, "formulation", "A-Phi"),
        ("em_msc", radia_em.EMWindow, "formulation", "MSC"),
        ("em_kelvin_bench", radia_em.EMWindow,
         "formulation", "Kelvin Benchmark"),
        # PCB — single layout
        ("pcb", radia_pcb.PCBWindow, None, None),
    ]


def _force_mode(window, combo_ref, value):
    """Apply a mode on the panel's combo, whether by attribute or key."""
    if combo_ref is None or value is None:
        return
    panel = window._panel
    combo = getattr(panel, combo_ref, None) \
            or (panel._widgets.get(combo_ref)
                if hasattr(panel, "_widgets") else None)
    if combo is not None:
        combo.setCurrentText(value)


def run_all_panel_checks(screenshot_dir="temp", strict=True):
    """Render every panel × mode in the registry, check each, save PNGs.

    Returns a dict {tag: [CheckResult, ...]}.
    Raises PanelQAError if any panel fails a mandatory check (strict=True).
    """
    from PySide6.QtWidgets import QApplication
    app = QApplication.instance() or QApplication([])
    # Mirror the real panel runtime: apply the lab-standard 11pt
    # baseline before instantiating windows.  Without this, every
    # widget renders at Qt's 9pt OS default and check_font_size_min
    # fails by design.
    from radia_gui_base import apply_panel_base_font
    apply_panel_base_font(app)
    all_results = {}
    fails = []
    for tag, WinCls, combo_ref, value in get_panel_registry():
        window = WinCls("")
        _force_mode(window, combo_ref, value)
        try:
            results, _ = run_panel_checks(tag, window, screenshot_dir,
                                           strict=False)
            all_results[tag] = results
            bad = [r for r in results
                    if not r.ok and r.severity == "error"]
            if bad:
                fails.append((tag, [r.name for r in bad]))
        finally:
            window.close()
            window.deleteLater()
    if fails and strict:
        raise PanelQAError(f"panels failing checks: {fails}")
    return all_results
