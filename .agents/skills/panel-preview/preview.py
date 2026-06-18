"""panel-preview/preview.py -- render any Radia Layer 3 panel + screenshot.

See SKILL.md for usage.  In short::

    python .claude/skills/panel-preview/preview.py
    python .claude/skills/panel-preview/preview.py --panel radia_ih
    python .claude/skills/panel-preview/preview.py --panel radia_ih --method "PEEC + BEM..."
    python .claude/skills/panel-preview/preview.py --real-qt

Output: ``C:/temp/panel_preview/<panel>_<method-tag>.png`` and ``.txt``.
The .txt file lists the visible widgets per method so AI agents can
verify "did widget X appear" without rendering the PNG.

This script is import-safe -- you can also call ``capture_panel(...)``
directly from a Python session if running ``preview.py`` is awkward
(e.g. inside ``Bash`` that swallows multi-line commands).
"""
from __future__ import annotations

import argparse
import os
import re
import sys
import time

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                          "..", "..", ".."))
SRC_RADIA = os.path.join(REPO_ROOT, "src", "radia")
SRC_RADIA_PANELS = os.path.join(SRC_RADIA, "panels")
DEFAULT_OUTPUT_DIR = "C:/temp/panel_preview"


# Map of panel-module -> WindowClass.  Add new panels here as they
# come online.
PANEL_REGISTRY = {
    "radia_ih":   "IHWindow",
    "radia_em":   "EMWindow",
    "radia_pcb":  "PCBWindow",
    "radia_heat": "HeatWindow",
}


def _slugify(text: str) -> str:
    """Convert method label to a filename-safe tag (lowercase, _-only)."""
    s = re.sub(r"[^A-Za-z0-9]+", "_", text).strip("_").lower()
    return s[:60] if len(s) > 60 else s


def _ensure_paths():
    """Insert src/radia and src/radia/panels into sys.path so panel
    modules can be imported standalone."""
    for p in (SRC_RADIA, SRC_RADIA_PANELS):
        if p not in sys.path:
            sys.path.insert(0, p)


def discover_methods(panel_module: str) -> list[str]:
    """Return the list of METHOD_* string constants exposed by the
    panel module.  Empty list if the panel has no method enum
    (e.g. radia_heat is single-mode)."""
    _ensure_paths()
    mod = __import__(panel_module)
    methods = []
    for name in dir(mod):
        if name.startswith("METHOD_") and name not in (
                "METHOD_TOOLTIP",):
            val = getattr(mod, name)
            if isinstance(val, str):
                methods.append(val)
    return methods


def _walk_visible_widgets(window):
    """Return [(label, class_name, hidden_str), ...] for every QLabel,
    QLineEdit, QSpinBox, QComboBox, QDoubleSpinBox, QPushButton in the
    window.  Hidden widgets are kept (so the report can confirm a
    method-switch hides them) with hidden_str = '[hidden]'."""
    from PySide6.QtWidgets import (QLabel, QLineEdit, QSpinBox,
                                    QComboBox, QDoubleSpinBox,
                                    QPushButton, QCheckBox)
    out = []
    for cls in (QLabel, QLineEdit, QSpinBox, QDoubleSpinBox, QComboBox,
                QPushButton, QCheckBox):
        for w in window.findChildren(cls):
            text = ""
            try:
                if isinstance(w, QLabel):
                    text = w.text()
                elif isinstance(w, (QSpinBox, QDoubleSpinBox)):
                    text = f"value={w.value()}"
                elif isinstance(w, QLineEdit):
                    text = w.text()
                elif isinstance(w, QComboBox):
                    text = f"current={w.currentText()!r} items={[w.itemText(i) for i in range(w.count())]}"
                elif isinstance(w, QPushButton):
                    text = w.text()
                elif isinstance(w, QCheckBox):
                    text = f"{w.text()} checked={w.isChecked()}"
            except Exception:
                text = "<error>"
            visible = "[visible]" if w.isVisible() else "[hidden]"
            out.append((text, cls.__name__, visible))
    return out


def _settle(app, window, n_iter=20, sleep_s=0.05):
    """Pump the Qt event loop a few times so geometry / fonts settle
    before grab().  Empirically 10-20 processEvents passes covers
    layout, dynamic visibility hooks, and tooltip render."""
    for _ in range(n_iter):
        app.processEvents()
        time.sleep(sleep_s)


def capture_panel(panel_module: str,
                   methods: list[str] | None = None,
                   output_dir: str = DEFAULT_OUTPUT_DIR,
                   real_qt: bool = False,
                   vol_path: str = "") -> list[str]:
    """Render ``panel_module`` and save one screenshot per method.

    Args:
        panel_module: e.g. ``"radia_ih"``
        methods: list of METHOD_* string values to capture, or None to
            auto-discover via ``discover_methods``.  Empty list for
            panels without methods (single screenshot).
        output_dir: where to write PNG + TXT pairs.
        real_qt: if True, do NOT force ``QT_QPA_PLATFORM=offscreen``;
            uses the default desktop platform (windows/X11).  Required
            for font-metric-accurate screenshots.

    Returns:
        list of PNG paths written.
    """
    _ensure_paths()
    if not real_qt:
        os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    os.makedirs(output_dir, exist_ok=True)

    from PySide6.QtWidgets import QApplication
    app = QApplication.instance() or QApplication([])

    win_cls_name = PANEL_REGISTRY.get(panel_module)
    if win_cls_name is None:
        raise ValueError(
            f"Unknown panel {panel_module!r}.  Known: "
            f"{sorted(PANEL_REGISTRY)}")
    mod = __import__(panel_module)
    WindowCls = getattr(mod, win_cls_name)

    if methods is None:
        methods = discover_methods(panel_module)
    if not methods:
        methods = [None]   # single capture

    written = []
    for method in methods:
        try:
            w = WindowCls(vol_path)
        except TypeError:
            w = WindowCls()
        # Switch method if applicable.
        panel = getattr(w, "_panel", None)
        if method is not None and panel is not None and \
                hasattr(panel, "_method_combo"):
            panel._method_combo.setCurrentText(method)
            if hasattr(panel, "_on_method_changed"):
                panel._on_method_changed(method)
        w.show()
        _settle(app, w)

        tag = _slugify(method) if method else "default"
        png_path = os.path.join(output_dir, f"{panel_module}_{tag}.png")
        txt_path = os.path.join(output_dir, f"{panel_module}_{tag}.txt")

        pix = w.grab()
        ok = pix.save(png_path)
        if not ok:
            sys.stderr.write(f"WARN: pix.save({png_path}) returned False\n")
            continue
        widgets = _walk_visible_widgets(w)
        with open(txt_path, "w", encoding="utf-8") as f:
            f.write(f"# Panel: {panel_module}\n")
            f.write(f"# Method: {method!r}\n")
            f.write(f"# Window class: {win_cls_name}\n")
            f.write(f"# Real Qt: {real_qt}\n")
            f.write(f"# Widgets ({len(widgets)} total):\n")
            for text, cls, visible in widgets:
                f.write(f"  {visible:10s} {cls:20s} {text}\n")
        written.append(png_path)
        print(f"  written: {png_path}  ({len(widgets)} widgets)")
        w.close()
        app.processEvents()
    return written


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--panel", default=None,
                    help=f"panel module name (one of "
                         f"{sorted(PANEL_REGISTRY)}); default: all")
    p.add_argument("--method", default=None,
                    help="single method label to capture (default: all "
                         "methods of the panel)")
    p.add_argument("--output-dir", default=DEFAULT_OUTPUT_DIR,
                    help=f"output directory (default: {DEFAULT_OUTPUT_DIR})")
    p.add_argument("--real-qt", action="store_true", default=True,
                    help="(default ON)  Use the default desktop Qt "
                         "platform.  This produces PNGs with readable "
                         "fonts; offscreen Qt on Windows renders all "
                         "labels as ☐ box characters because the offscreen "
                         "platform plugin does not load system fonts.  "
                         "Pass --offscreen to override (CI / no-display).")
    p.add_argument("--offscreen", action="store_true",
                    help="force the offscreen Qt platform (CI / SSH "
                         "without display).  Note: PNGs will have ☐ box "
                         "characters instead of readable text on Windows; "
                         "the companion .txt file is the source of truth "
                         "for widget visibility in that case.")
    p.add_argument("--vol-path", default="",
                    help="optional .vol file path passed to the window "
                         "constructor (some panels read it on launch)")
    args = p.parse_args()

    # --offscreen overrides --real-qt (default).
    real_qt = args.real_qt and not args.offscreen

    panels = [args.panel] if args.panel else sorted(PANEL_REGISTRY)
    all_written = []
    for panel in panels:
        print(f"=== {panel} ===")
        methods = [args.method] if args.method else None
        try:
            written = capture_panel(panel, methods=methods,
                                     output_dir=args.output_dir,
                                     real_qt=real_qt,
                                     vol_path=args.vol_path)
            all_written.extend(written)
        except Exception as exc:
            print(f"  FAILED: {exc}")
            import traceback
            traceback.print_exc()

    print(f"\nTotal: {len(all_written)} screenshots in {args.output_dir}")
    return 0 if all_written else 1


if __name__ == "__main__":
    sys.exit(main())
