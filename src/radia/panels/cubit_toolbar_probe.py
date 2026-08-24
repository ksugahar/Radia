#!python
"""Probe the Coreform-owned Radia Export toolbar in a real Cubit GUI.

This file is executed inside Cubit's embedded Python by
``radia.cubit_toolbar_smoke``.  Keep it standalone: importing the normal
``radia`` package here would load Python-3.12 extension modules into Cubit's
private Python runtime.
"""

from __future__ import annotations

import json
import os
import time

from PySide6.QtCore import QTimer
from PySide6.QtWidgets import QApplication, QMainWindow, QToolBar

SCHEMA = "radia.cubit-toolbar-probe.v1"
EXPECTED_ACTIONS = [
    "Netgen Vol (.vol)",
    "GMSH (.msh)",
    "Nastran (.bdf)",
    "VTK (.vtk)",
    "FEMEEM",
    "MEG (ELF/MAGIC)",
]
RESULT_PATH = os.environ.get(
    "RADIA_TOOLBAR_PROBE_RESULT",
    r"C:\temp\radia_toolbar_probe_result.json",
)
TIMEOUT_SECONDS = float(os.environ.get("RADIA_TOOLBAR_PROBE_TIMEOUT", "20"))

_started = time.monotonic()
_finished = False


def _plain(text):
    return str(text).replace("&", "")


def _menu_entries(menu):
    entries = []
    for action in menu.actions():
        entries.append(_plain(action.text()))
        submenu = action.menu()
        if submenu is not None:
            entries.extend(_menu_entries(submenu))
    return entries


def _finish(payload):
    global _finished
    if _finished:
        return
    _finished = True
    payload["schema"] = SCHEMA
    payload["elapsed_seconds"] = round(time.monotonic() - _started, 3)
    os.makedirs(os.path.dirname(os.path.abspath(RESULT_PATH)), exist_ok=True)
    with open(RESULT_PATH, "w", encoding="utf-8") as stream:
        json.dump(payload, stream, ensure_ascii=False, indent=2, sort_keys=True)

    # Run the exit command after the result has been flushed.  The external
    # runner verifies that the GUI process actually disappears.
    import cubit
    QTimer.singleShot(0, lambda: cubit.cmd("exit"))


def _snapshot():
    app = QApplication.instance()
    if app is None or not hasattr(app, "topLevelWidgets"):
        return None

    main_windows = [
        widget for widget in app.topLevelWidgets()
        if isinstance(widget, QMainWindow)
    ]
    main = next(
        (widget for widget in main_windows
         if "Coreform Cubit" in widget.windowTitle()),
        None,
    )
    if main is None:
        return None

    candidates = [
        toolbar for toolbar in main.findChildren(QToolBar)
        if _plain(toolbar.objectName()) == "Radia Export"
        or _plain(toolbar.windowTitle()) == "Radia Export"
    ]
    toolbar = candidates[0] if candidates else None
    actions = list(toolbar.actions()) if toolbar is not None else []
    action_names = [_plain(action.text()) for action in actions]

    top_level_menu = [_plain(action.text()) for action in main.menuBar().actions()]
    view_action = next(
        (action for action in main.menuBar().actions()
         if _plain(action.text()) == "View"),
        None,
    )
    view_entries = []
    if view_action is not None and view_action.menu() is not None:
        view_entries = _menu_entries(view_action.menu())

    popup = main.createPopupMenu()
    popup_entries = _menu_entries(popup) if popup is not None else []
    if popup is not None:
        popup.deleteLater()

    toolbar_visible = bool(toolbar is not None and toolbar.isVisible())
    visible_region_nonempty = bool(
        toolbar is not None and not toolbar.visibleRegion().isEmpty()
    )
    size = [0, 0] if toolbar is None else [
        int(toolbar.width()), int(toolbar.height())
    ]
    action_visible = {
        name: bool(action.isVisible())
        for name, action in zip(action_names, actions)
    }
    action_enabled = {
        name: bool(action.isEnabled())
        for name, action in zip(action_names, actions)
    }
    toolbar_menu_has = (
        "Radia Export" in view_entries or "Radia Export" in popup_entries
    )

    payload = {
        "main_window_visible": bool(main.isVisible()),
        "toolbar_count": len(candidates),
        "toolbar_visible": toolbar_visible,
        "toolbar_visible_region_nonempty": visible_region_nonempty,
        "toolbar_size": size,
        "toolbar_actions": action_names,
        "action_visible": action_visible,
        "action_enabled": action_enabled,
        "toolbar_menu_has_radia_export": toolbar_menu_has,
        "view_has_radia_export": "Radia Export" in view_entries,
        "unsupported_top_level_menu_present": "Radia Export" in top_level_menu,
    }
    payload["ok"] = (
        payload["main_window_visible"]
        and payload["toolbar_count"] == 1
        and payload["toolbar_visible"]
        and payload["toolbar_visible_region_nonempty"]
        and all(value > 0 for value in payload["toolbar_size"])
        and payload["toolbar_actions"] == EXPECTED_ACTIONS
        and all(payload["action_visible"].get(name) for name in EXPECTED_ACTIONS)
        and all(payload["action_enabled"].get(name) for name in EXPECTED_ACTIONS)
        and payload["toolbar_menu_has_radia_export"]
        and not payload["unsupported_top_level_menu_present"]
    )
    return payload


def _probe():
    payload = _snapshot()
    if payload is not None and payload.get("ok"):
        _finish(payload)
        return
    if time.monotonic() - _started < TIMEOUT_SECONDS:
        QTimer.singleShot(250, _probe)
        return
    if payload is None:
        payload = {"ok": False, "error": "Cubit main window missing"}
    _finish(payload)


QTimer.singleShot(0, _probe)
