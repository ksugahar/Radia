"""
test_combo_wheel_guard.py -- detect the "mouse-wheel scroll silently
changes a combo / spin value" pitfall.

A bare QComboBox / QSpinBox / QDoubleSpinBox grabs mouse-wheel events
and changes its value -- even when the user only meant to SCROLL THE
FORM.  In a tall parameter panel this silently mutates a selection the
user already set and believes is fixed, so the run uses an input the
user never intended and cannot easily confirm.

Reported by kubota 2026-05-29: after clicking the "Workpiece impedance
model" combo to inspect it, scrolling the IH form kept changing the
selection unnoticed.

radia_gui_base ships ``_NoWheelComboBox`` / ``_NoWheelSpinBox`` /
``_NoWheelDoubleSpinBox`` which ignore the wheel event (it propagates
to the panel's QScrollArea so the FORM scrolls); the value changes
only via click+dropdown / arrow keys.  ``add_combo`` / ``add_spin``
use them, and the non-ModePanel panels (radia_motor) alias-import them.

This test instantiates every panel window and asserts NO combo /
spinbox changes value on a wheel notch -- focused or not.  It is
behavioural (sends a real QWheelEvent and checks the value), so it is
immune to HOW the guard is implemented (subclass, alias, event filter).

Run::

    QT_QPA_PLATFORM=offscreen python -m pytest validation_test/panels/test_combo_wheel_guard.py -v
"""

from __future__ import annotations

import os

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import importlib

import pytest
from PySide6.QtCore import Qt, QPoint, QPointF
from PySide6.QtGui import QWheelEvent
from PySide6.QtWidgets import QApplication, QComboBox, QAbstractSpinBox

# Panel entry points: each *Window builds the full panel + every tab,
# so findChildren() below sees every combo/spin a user can reach.
PANEL_WINDOWS = [
    ("radia_ih", "IHWindow"),
    ("radia_em", "EMWindow"),
    ("radia_pcb", "PCBWindow"),
    ("radia_motor", "MotorWindow"),
]


def _send_wheel(widget, dy):
    """Send one mouse-wheel notch to *widget* (dy<0 = down, dy>0 = up)."""
    ev = QWheelEvent(
        QPointF(5, 5),
        QPointF(widget.mapToGlobal(QPoint(5, 5))),
        QPoint(0, 0),
        QPoint(0, dy),
        Qt.MouseButton.NoButton,
        Qt.KeyboardModifier.NoModifier,
        Qt.ScrollPhase.NoScrollPhase,
        False,
    )
    QApplication.sendEvent(widget, ev)


def _combo_steals(w):
    """True if the combo changes index on a wheel notch (either dir)."""
    if w.count() < 2:
        return False  # nothing to change into
    for focused in (False, True):
        w.setFocus() if focused else w.clearFocus()
        before = w.currentIndex()
        _send_wheel(w, -120)              # scroll down
        if w.currentIndex() != before:
            return True
        _send_wheel(w, +120)              # scroll up
        if w.currentIndex() != before:
            return True
    return False


def _spin_steals(w):
    """True if the spinbox changes value on a wheel notch (either dir)."""
    lo, hi = w.minimum(), w.maximum()
    if lo == hi:
        return False
    for focused in (False, True):
        w.setFocus() if focused else w.clearFocus()
        before = w.value()
        _send_wheel(w, -120)
        if w.value() != before:
            return True
        _send_wheel(w, +120)
        if w.value() != before:
            return True
    return False


def _name(w):
    return w.objectName() or f"<{type(w).__name__}#{id(w) & 0xffff:x}>"


def _wheel_offenders(root):
    """Return the names of combos/spins that change value on wheel."""
    bad = []
    for w in root.findChildren(QComboBox):
        try:
            if _combo_steals(w):
                bad.append(_name(w))
        except RuntimeError:
            pass  # widget deleted by a mode-switch handler mid-scan
    for w in root.findChildren(QAbstractSpinBox):
        try:
            if _spin_steals(w):
                bad.append(_name(w))
        except RuntimeError:
            pass
    return bad


@pytest.mark.parametrize("module_name,cls_name", PANEL_WINDOWS)
def test_panel_window_no_wheel_steal(qapp, module_name, cls_name):
    """Every combo/spin in every panel window ignores wheel scrolling."""
    mod = importlib.import_module(module_name)
    win = getattr(mod, cls_name)(vol_path="")
    try:
        offenders = _wheel_offenders(win)
        assert not offenders, (
            f"{cls_name}: these combos/spins change value on mouse-wheel "
            f"scroll (scroll-steal pitfall -- use _NoWheel* from "
            f"radia_gui_base): {offenders}"
        )
    finally:
        win.deleteLater()


def test_ih_window_combos_are_wheel_guarded(ih_window):
    """Explicit cover of the reported panel (IH workpiece-impedance combo)."""
    offenders = _wheel_offenders(ih_window)
    assert not offenders, f"IHWindow scroll-steal: {offenders}"


def test_add_combo_and_add_spin_are_guarded(ih_panel):
    """White-box: the base helpers build the guarded subclasses, so any
    panel built via add_combo/add_spin inherits the wheel guard."""
    import radia_gui_base as gb
    combo = ih_panel.add_combo("_wheelguard_probe_c", "probe", ["a", "b"], 0)
    spin = ih_panel.add_spin("_wheelguard_probe_s", "probe", 1, 1, 9)
    assert isinstance(combo, gb._NoWheelComboBox), (
        "add_combo must build _NoWheelComboBox")
    assert isinstance(spin, gb._NoWheelSpinBox), (
        "add_spin must build _NoWheelSpinBox")
