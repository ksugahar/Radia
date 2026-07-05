"""Panel quality-assurance checks.

A single entry point ``check_panel_health(window)`` runs a suite of
layout / accessibility / correctness checks on a Radia PySide6 panel
window.  Used by:

  - the ``/deploy`` skill's Step 8 (pre-deploy visual check)
  - the ``pytest validation_test/panels/test_panel_qa.py`` regression suite

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
# Thresholds reflect the 13pt baseline font (Segoe UI).  9pt-era numbers
# (1100 / 900) under-budgeted by ~45% once the user said "ケチるな"
# (don't be stingy) -- a comfortable panel on 2K is allowed to take
# half the screen.  Hard fail still leaves room for output area
# scrolling.
#
# 2026-05-02 update: hard limit relaxed 1350 -> 1700 px after IH BEM /
# Heat 3D / Heat axisym all hit 1450-1610 px in production with the
# v4.18.0+ feature set (back-reaction, axisym mode, etc.). All three
# panels are usable on 2K at the relaxed limit (some scrolling needed
# below 1700 vertical, but Run button reachable).  4K display has
# plenty of room.  Yellow threshold raised in step (1500 px = "design
# warning") so panels don't grow unchecked between releases.
MAX_HEIGHT_RED = 1700     # hard fail: panel cannot fit on 2K even with scrolling
MAX_HEIGHT_YELLOW = 1500  # warning: getting tall, consider tabs or sub-grouping
MAX_WIDTH_RED = 1400      # hard fail: panel eats too much horizontal real estate
MAX_WIDTH_YELLOW = 1200   # warning: method combo text probably long

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
    """Run button's bottom must be within the window's *actual* height.

    Does NOT mutate window size — previous implementation called
    ``window.resize()`` which inflated subsequent sizeHint() calls
    for other checks.  Use ``geometry()`` (which is valid after
    ``adjustSize()`` ran in ``check_panel_health``) to read the
    button's y position.

    Compares against ``max(sizeHint.h, minimumSize.h)`` because Qt
    enforces ``minimumSize`` even when ``sizeHint`` is smaller — for
    panels with few form rows (PCB has 6) the layout sizeHint can
    fall below ``minimumSize`` and Qt forces the window to expand
    to ``minimumSize``, placing the Run button toward the bottom of
    the actual rendered window (NOT the sizeHint).
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
    # Actual rendered height = max(sizeHint, minimumSize) after adjustSize
    win_h = max(window.sizeHint().height(),
                window.minimumSize().height())
    if btn_bottom > win_h + 2:  # 2 px tolerance for frame/padding
        return CheckResult("buttons_reachable", False,
                            f"Run bottom y={btn_bottom} > window h={win_h} "
                            f"(sizeHint={window.sizeHint().height()}, "
                            f"minimumSize={window.minimumSize().height()})")
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


def check_run_button_unique(window) -> CheckResult:
    """Exactly one QPushButton labelled 'Run' (exact match after strip).

    Catches accidental duplicate Run buttons left behind during a
    panel refactor.  Secondary actions ('Run thermal...', 'Run sweep')
    are intentionally excluded by requiring exact 'Run' — they are
    distinct user actions and SHOULD coexist with the main Run.
    """
    runs = [b for b in window.findChildren(QPushButton)
            if b.text().strip() == "Run"]
    if len(runs) == 1:
        return CheckResult("run_button_unique", True,
                            f"text={runs[0].text()!r}")
    return CheckResult("run_button_unique", False,
                        f"found {len(runs)} 'Run' buttons (exact match): "
                        f"{[b.text() for b in runs]}")


def check_widget_to_cli_coverage(window) -> CheckResult:
    """Every visible user-input widget must reach build_command() output.

    Catches the common "widget added to UI but forgotten in
    build_command" silent bug -- the user types a value, hits Run, and
    the calc subprocess never sees the change.

    Algorithm:
      1. For each panel widget key K (from panel._widgets):
         - skip if widget is hidden in current mode
         - skip if K is in panel._cli_orphan_keys (explicit whitelist
           for GUI-only widgets like preset combos that drive *other*
           widgets rather than emitting their own flag)
         - skip if K is in panel._cli_value_only_keys (whitelist for
           widgets whose value appears in cmd but key does not -- e.g.
           filenames passed positionally)
      2. Build the CLI via panel.build_command("").  Skip the check
         entirely if build_command raises (some panels require a real
         file path; that's caught by other tests).
      3. The widget is "covered" if any of these tokens appears in cmd:
            --{key.replace('_', '-')}
            --{key}
            the widget's val() (after kebab-casing for combos)
    """
    panel = getattr(window, "_panel", None)
    if panel is None or not hasattr(panel, "_widgets"):
        return CheckResult("widget_to_cli_coverage", True,
                            "(no panel / no _widgets)")

    if not hasattr(panel, "build_command"):
        return CheckResult("widget_to_cli_coverage", True,
                            "(no build_command)")

    try:
        cmd = panel.build_command("")
    except Exception as e:
        return CheckResult("widget_to_cli_coverage", True,
                            f"(build_command raised: {type(e).__name__}; "
                            f"skipped)")

    cmd_text = " ".join(str(x) for x in cmd)

    orphan_keys = set(getattr(panel, "_cli_orphan_keys", ()) or ())
    value_only_keys = set(
        getattr(panel, "_cli_value_only_keys", ()) or ())

    missing = []
    for key, w in panel._widgets.items():
        if key in orphan_keys:
            continue
        if not w.isVisibleTo(window):
            continue

        kebab = key.replace("_", "-")
        if f"--{kebab}" in cmd_text or f"--{key}" in cmd_text:
            continue

        # Fall back to value match (for value-only widgets and
        # the heat panel's wp_vol that is passed as --wp-vol but
        # we don't want to re-encode the kebab logic for paths).
        try:
            v = panel.val(key)
        except Exception:
            v = ""
        if v and str(v) in cmd_text:
            continue

        if key in value_only_keys:
            continue

        missing.append(key)

    if missing:
        return CheckResult("widget_to_cli_coverage", False,
                            f"visible widgets not in build_command: "
                            f"{missing}",
                            severity="warning")
    return CheckResult("widget_to_cli_coverage", True)


def check_visible_rows_have_labels(window) -> CheckResult:
    """Every visible input row must have a non-empty label in the
    label column.  Empty labels make the form ambiguous.

    Legitimate exemptions:
      - Section header QLabels (handled by _is_section_header_label).
      - Plain-text QLabel fields (status lines, footers).
      - Sub-panel widgets that carry their own QFormLayout — they
        label themselves internally via their own form rows
        (e.g. HeatPanel embedded in IHPanel).
      - QCheckBox with a non-empty text() — the checkbox text IS
        the label by Qt convention; a left-column label would be a
        duplicate (e.g. "Override rho/cp/k").
    """
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
            # Skip sub-panel widgets that have their own QFormLayout —
            # they label themselves internally (e.g. HeatPanel inside
            # IHPanel).  A parent-side label would be redundant with
            # the section header that introduces the sub-panel.
            if isinstance(fw.layout(), QFormLayout):
                continue
            # Skip QScrollArea-wrapped sub-panels: the scroll viewport
            # contains a widget that itself uses QFormLayout for
            # labels (e.g. IHPanel wraps HeatPanel in QScrollArea to
            # bound the window height).
            from PySide6.QtWidgets import QScrollArea
            if isinstance(fw, QScrollArea):
                inner = fw.widget()
                if inner is not None and isinstance(inner.layout(), QFormLayout):
                    continue
            # Skip self-labeled QCheckBoxes — the checkbox text is
            # the label by Qt convention.
            from PySide6.QtWidgets import QCheckBox
            if isinstance(fw, QCheckBox) and fw.text().strip():
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
    check_run_button_unique,
    check_no_orphan_section_headers,
    check_ascii_only_labels,
    check_method_combo_populated,
    check_visible_rows_have_labels,
    check_font_size_min,
    # check_widget_to_cli_coverage is deliberately NOT in the default
    # set: it requires a populated .vol/.step path for build_command()
    # to run, which the deploy-skill render flow does not provide.
    # test_build_command_parses.py invokes it explicitly with fixtures.
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
        # IH — 4 methods (PEEC-IND / PEEC-BEM / FEM-full / Thermal).
        # Thermal is the post-4.59.0 home for heat analysis (formerly
        # the standalone radia_heat HeatWindow).
        ("ih_ind", radia_ih.IHWindow,
         "_method_combo", radia_ih.METHOD_PEEC_IND),
        ("ih_bem", radia_ih.IHWindow,
         "_method_combo", radia_ih.METHOD_PEEC_BEM),
        ("ih_fem", radia_ih.IHWindow,
         "_method_combo", radia_ih.METHOD_FEM_FULL),
        ("ih_thermal_3d_static", radia_ih.IHWindow,
         "_method_combo", radia_ih.METHOD_THERMAL_3D_STATIC),
        ("ih_thermal_3d_rotating", radia_ih.IHWindow,
         "_method_combo", radia_ih.METHOD_THERMAL_3D_ROTATING),
        ("ih_thermal_axisym", radia_ih.IHWindow,
         "_method_combo", radia_ih.METHOD_THERMAL_AXISYM),
        # EM — formulations (Omega / A-Phi / HDiv-VIM / Kelvin Benchmark).
        # _method_combo is the convention-name attribute exposed by
        # every mode-switching panel (EM, IH, ...).
        ("em_omega", radia_em.EMWindow, "_method_combo", "Omega"),
        ("em_aphi", radia_em.EMWindow, "_method_combo", "A-Phi"),
        ("em_hdiv", radia_em.EMWindow, "_method_combo", "HDiv-VIM"),
        ("em_kelvin_bench", radia_em.EMWindow,
         "_method_combo", "Kelvin Benchmark"),
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
