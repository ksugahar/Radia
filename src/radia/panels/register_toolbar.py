"""
Radia Cubit panel startup hook.

Runs inside Cubit on startup (via ~/.cubit -> startup.py) and performs:
  1. Panel-log initialization (shared across all radia subprocesses)
  2. Export-dialog default-directory seeding for the C++ .ccl component
  3. Cleanup of any legacy Python-side "Solve" / "Radia-NGSolve" /
     "Generate Coil" / "Reload Panels" menus left behind by older installs

All user-facing menus (Export Mesh + Radia-NGSolve analysis launcher) now
live in the C++ .ccl component (src/cubit_plugin/RadiaComp.cpp). This file
no longer creates any menus.
"""

import json
import os
import sys

# ----------------------------------------------------------------------
# Panel debug log: shared writer in radia.panels.panel_log.
# Cubit-side processes truncate the log on session start so one
# Cubit session = one continuous log file across all subprocesses.
# ----------------------------------------------------------------------
_my_dir = os.path.dirname(os.path.abspath(__file__)) if "__file__" in dir() \
    else os.getcwd()
if _my_dir not in sys.path:
    sys.path.insert(0, _my_dir)
from panel_log import (init_panel_log, panel_log as _panel_log,
                       panel_log_exception as _panel_log_exception,
                       PANEL_LOG_PATH as _PANEL_LOG_PATH)

init_panel_log("cubit", truncate=True, banner=True)

try:
    _self_mtime = os.path.getmtime(__file__) if "__file__" in dir() else 0
    _self_path = __file__ if "__file__" in dir() else "(no __file__)"
    _panel_log(f"  this file: {_self_path}")
    _panel_log(f"  this mtime: {int(_self_mtime)}")
except Exception:
    pass
print(f"[Radia] Panel debug log: {_PANEL_LOG_PATH}")

# Determine script location (__file__ is not defined when run via Cubit 'play')
try:
    _this_dir = os.path.dirname(os.path.abspath(__file__))
except NameError:
    import importlib
    try:
        _radia_spec = importlib.util.find_spec("radia")
        if _radia_spec and _radia_spec.origin:
            _this_dir = os.path.join(os.path.dirname(_radia_spec.origin),
                                     "panels")
        else:
            raise ImportError
    except (ImportError, AttributeError):
        import glob as _glob
        for _base in [
            os.path.join(sys.prefix, "Lib", "site-packages"),
            os.path.join(sys.prefix, "lib", "python*", "site-packages"),
        ]:
            for _candidate in _glob.glob(
                    os.path.join(_base, "radia", "panels")):
                if os.path.isdir(_candidate):
                    _this_dir = _candidate
                    break
            else:
                continue
            break
        else:
            _this_dir = os.getcwd()

# Qt bindings: prefer PySide6, fall back to PyQt5 (Cubit ships PyQt5)
try:
    from PySide6.QtWidgets import QApplication, QToolBar
    _QT = "PySide6"
except ImportError:
    from PyQt5.QtWidgets import QApplication, QToolBar  # noqa: F401
    _QT = "PyQt5"


def _find_main_window():
    """Find Cubit's top-level QMainWindow for menu cleanup."""
    app = QApplication.instance()
    if app is None:
        return None
    if _QT == "PySide6":
        from PySide6.QtWidgets import QMainWindow
    else:
        from PyQt5.QtWidgets import QMainWindow
    for w in app.topLevelWidgets():
        if isinstance(w, QMainWindow):
            return w
    return None


def _unc_to_drive(path):
    """Convert UNC path to drive letter if mapped (Windows)."""
    if sys.platform != "win32" or not path.startswith("\\\\"):
        return path
    try:
        import subprocess
        result = subprocess.run(
            ["net", "use"], capture_output=True, text=True, timeout=5)
        for line in result.stdout.splitlines():
            parts = line.split()
            if len(parts) >= 3 and parts[0].endswith(":"):
                drive, unc = parts[0], parts[1]
                if path.lower().startswith(unc.lower()):
                    return drive + path[len(unc):]
    except Exception:
        pass
    return path


def _samples_dir():
    """Return the package samples directory (default working folder)."""
    return _unc_to_drive(os.path.join(_this_dir, "samples"))


def _init_export_default_dir():
    """Write default_dir to C++ export_settings.json.

    The C++ ExportDialog reads default_dir as fallback when no
    journal path and no saved dir exist. Without this, it falls
    back to _getcwd() which may be OneDrive.
    """
    appdata = os.path.join(os.path.expanduser("~"),
                           "AppData", "Roaming", "Radia")
    settings_path = os.path.join(appdata, "export_settings.json")
    samples = _samples_dir().replace("\\", "/")
    try:
        os.makedirs(appdata, exist_ok=True)
        data = {}
        if os.path.isfile(settings_path):
            with open(settings_path, "r", encoding="utf-8") as f:
                data = json.load(f)
        data["default_dir"] = samples
        with open(settings_path, "w", encoding="utf-8") as f:
            json.dump(data, f, indent=2)
    except (OSError, json.JSONDecodeError):
        pass


def _cleanup_legacy_menus():
    """Remove legacy Python-side menus left behind by older installs.

    Menus owned by the C++ .ccl component (Export Mesh) are NOT removed.
    """
    main_window = _find_main_window()
    if main_window is None:
        _panel_log("_cleanup_legacy_menus: main_window is None — skip")
        return
    menu_bar = main_window.menuBar()
    removed = []
    for name in ("Solve", "Radia-NGSolve", "Generate Coil",
                 "Reload Panels"):
        for action in list(menu_bar.actions()):
            if action.text().replace("&", "") == name:
                sub = action.menu()
                menu_bar.removeAction(action)
                if sub:
                    sub.deleteLater()
                removed.append(name)
    for tb in main_window.findChildren(QToolBar, "RadiaToolBar"):
        main_window.removeToolBar(tb)
        tb.deleteLater()
        removed.append("RadiaToolBar")
    app = QApplication.instance()
    if app is not None:
        app.setQuitOnLastWindowClosed(False)
    if removed:
        _panel_log(f"_cleanup_legacy_menus: removed {removed}")


def register_menu():
    """Initialize export defaults + cleanup any legacy Solve menus.

    Despite the name, this function no longer registers any Python menus.
    All user-facing menus are in the C++ .ccl component. The name is
    kept for compatibility with startup.py which calls register_menu().
    """
    _panel_log("register_menu: ENTER (C++-only mode)")
    _cleanup_legacy_menus()
    _init_export_default_dir()
    _panel_log("register_menu: EXIT (export defaults seeded)")


# Auto-register when this script is executed
try:
    register_menu()
except Exception:
    _panel_log_exception("register_menu top-level FAILED")
    raise
