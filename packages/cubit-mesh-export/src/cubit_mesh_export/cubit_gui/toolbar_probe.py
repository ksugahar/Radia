#!python
"""Probe the Coreform-owned Cubit Mesh Export toolbar in a real Cubit GUI.

This file is executed inside Cubit's embedded Python by
``cubit_mesh_export.toolbar_smoke``. Keep it standalone: importing the normal
``radia`` package here would load Python-3.12 extension modules into Cubit's
private Python runtime.
"""

from __future__ import annotations

import json
import os
import time

from PySide6.QtCore import QTimer
from PySide6.QtWidgets import QApplication, QMainWindow, QToolBar

SCHEMA = "cubit-mesh-export.toolbar-probe.v2"
EXPECTED_ACTIONS = [
    "Netgen Vol (.vol)",
    "GMSH (.msh)",
    "Nastran (.bdf)",
    "VTK (.vtk)",
    "FEMEEM",
    "MEG (ELF/MAGIC)",
]
RESULT_PATH = os.environ.get(
    "CUBIT_MESH_EXPORT_TOOLBAR_PROBE_RESULT",
    r"C:\temp\cubit_mesh_export_toolbar_probe_result.json",
)
TIMEOUT_SECONDS = float(os.environ.get("CUBIT_MESH_EXPORT_TOOLBAR_PROBE_TIMEOUT", "20"))

_started = time.monotonic()
_finished = False


def _plain(text):
    return str(text).replace("&", "")


def _finish(payload):
    global _finished
    if _finished:
        return
    _finished = True
    payload["schema"] = SCHEMA
    payload["elapsed_seconds"] = round(time.monotonic() - _started, 3)
    os.makedirs(os.path.dirname(os.path.abspath(RESULT_PATH)), exist_ok=True)
    temporary_result = RESULT_PATH + ".tmp"
    with open(temporary_result, "w", encoding="utf-8") as stream:
        json.dump(payload, stream, ensure_ascii=False, indent=2, sort_keys=True)
    os.replace(temporary_result, RESULT_PATH)

    # Run the exit command after the result has been flushed.  The external
    # runner verifies that the GUI process actually disappears.
    import cubit
    QTimer.singleShot(0, lambda: cubit.cmd("exit 0"))


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
        if _plain(toolbar.objectName()) == "Cubit Mesh Export"
        or _plain(toolbar.windowTitle()) == "Cubit Mesh Export"
    ]
    toolbar = candidates[0] if candidates else None
    actions = list(toolbar.actions()) if toolbar is not None else []
    action_names = [_plain(action.text()) for action in actions]

    # This action belongs to QToolBar itself. Never traverse Claro's menubar:
    # wrapping its C++/SWIG menu actions with shiboken can corrupt GUI teardown.
    toggle = toolbar.toggleViewAction() if toolbar is not None else None

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
    toolbar_menu_has = toggle is not None and _plain(toggle.text()) == "Cubit Mesh Export"

    payload = {
        "main_window_visible": bool(main.isVisible()),
        "toolbar_count": len(candidates),
        "toolbar_visible": toolbar_visible,
        "toolbar_visible_region_nonempty": visible_region_nonempty,
        "toolbar_size": size,
        "toolbar_actions": action_names,
        "action_visible": action_visible,
        "action_enabled": action_enabled,
        "toolbar_menu_has_cubit_mesh_export": toolbar_menu_has,
        "toolbar_owner": toolbar.metaObject().className() if toolbar is not None else None,
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
        and payload["toolbar_menu_has_cubit_mesh_export"]
        and payload["toolbar_owner"] == "WorkflowToolbar"
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


if __name__ in ('__main__', '_coreform_cubit'):
    QTimer.singleShot(0, _probe)
