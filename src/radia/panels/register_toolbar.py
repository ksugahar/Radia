"""
Radia Cubit panel startup hook.

Runs inside Cubit on startup (via ~/.cubit -> startup.py) and performs:
  1. Panel-log initialization (shared across all radia subprocesses)
  2. Export-dialog default-directory seeding (export_settings.json)
  3. Cleanup of any legacy menus left behind by older installs
     (including the C++ Qt5 .ccl "Export Mesh" menu, which is now
     removed -- see radia 4.80.0 / src/radia/panels/radia_export_menu.py)
  4. Plugin freshness diagnostics
  5. Installation of the "Export" menu through Cubit's official Claro
     API (radia_export_menu.install_menu() -> emclaro.add_to_menu)

The menu is registered through ``emclaro``, the SWIG binding of Cubit's
Claro GUI framework that ships inside Cubit itself.  That is the same
``Claro::add_to_menu()`` entry point the old C++ ``.ccl`` component used,
so the menu is owned by Cubit and survives cold start.  The Coreform
WorkflowToolbar package remains available as the toolbar surface; the two
are complementary, not alternatives.

The Qt5 C++ Claro component (src/cubit_plugin/RadiaComp.cpp) was
deleted in radia 4.80.0; all GUI work is now PySide6.  The .ccm
(APREPRO commands `export gmsh / netgen / nastran / vtk /
femeem / meg`) is unchanged and remains the backend.
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

# Qt binding: PySide6 only (radia 4.80.0).  Target is Cubit 2025.12 which
# ships PySide6.  Per CLAUDE.md "No Fallbacks — Fail Fast, Fail Loud":
# if PySide6 is not importable, this script raises ImportError loudly so
# the operator sees the real issue (old Cubit, broken PySide6 install,
# wrong Python env) instead of a silent PyQt5 fallback that masks it.
from PySide6.QtWidgets import QApplication, QToolBar, QMainWindow


def _find_main_window():
    """Find Cubit's top-level QMainWindow for menu cleanup.

    Cubit sometimes runs this startup script before it has upgraded its
    Qt application from QCoreApplication (non-GUI) to QApplication
    (GUI). QCoreApplication has no ``topLevelWidgets`` method, so we
    have to bail out gracefully instead of crashing the whole startup.
    """
    app = QApplication.instance()
    if app is None or not hasattr(app, "topLevelWidgets"):
        return None  # no GUI app yet, nothing to clean up
    # MEASURED (Cubit 2025.12): TWO QMainWindows exist -- "Claro" (the
    # real main window) and "QtJournalEditor".  Returning the first match
    # handed back the journal editor, so legacy-menu cleanup silently
    # operated on the wrong menu bar.  Match the main window by name.
    for w in app.topLevelWidgets():
        if isinstance(w, QMainWindow) and w.objectName().lower() == "claro":
            return w
    for w in app.topLevelWidgets():
        if isinstance(w, QMainWindow) and w.objectName() != "QtJournalEditor":
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
    """Write default_dir to export_settings.json.

    Both the (deleted) C++ ExportDialog and the new PySide6
    radia_export_menu read `default_dir` from this file as the
    initial directory for the file dialog when no .jou is loaded.
    Without it, the dialog falls back to _getcwd() which may be
    OneDrive on Windows.
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


def _check_plugin_freshness():
    """Warn if Cubit C++ plugin (.ccm / .pyd) is older than this Python file.

    Silent mismatch between site-packages/radia and the Cubit plugin dir
    has caused hard-to-diagnose bugs (e.g., 2026-04-14: IH BEM failed to
    read source/sink sidesets because an outdated .ccm shipped without
    the FaceDescriptor DomainIn/Out fix).  Fail loudly.

    Note: the Qt5 .ccl is gone (radia 4.80.0).  Only .ccm (APREPRO
    commands) and the .pyd (cubit_mesh_curver) are deployed now.
    """
    # install_panels.py lives in the PARENT of this panels/ dir.  Import it
    # STANDALONE (put src/radia on sys.path) -- NEVER `from radia.install_panels
    # import ...`, which imports the radia PACKAGE and loads its Python-3.12
    # C++ extension into Cubit's bundled Python 3.10, CRASHING Cubit on startup
    # (CLAUDE.md: Layer 2 must NOT import radia / ngsolve -- DLL conflict; the
    # crash is a hard DLL-load failure that try/except CANNOT catch).
    try:
        _src = os.path.dirname(_this_dir)   # src/radia (parent of panels/)
        if _src and _src not in sys.path:
            sys.path.insert(0, _src)
        from install_panels import find_cubit_bin
    except Exception:
        return   # cannot locate install_panels standalone -> skip the check
    try:
        cubit_bin = find_cubit_bin()
    except Exception:
        cubit_bin = None
    if not cubit_bin:
        return
    plugin_dir = os.path.join(cubit_bin, "plugins")
    try:
        self_mtime = os.path.getmtime(__file__)
    except Exception:
        return
    tol_sec = 3600  # 1 h tolerance (one pip install batches within a minute)
    stale = []
    for name in ("cubit_mesh_export.ccm", "cubit_mesh_curver.cp312-win_amd64.pyd"):
        p = os.path.join(plugin_dir, name)
        if not os.path.isfile(p):
            continue
        m = os.path.getmtime(p)
        if m + tol_sec < self_mtime:
            lag_hr = (self_mtime - m) / 3600.0
            stale.append((name, p, m, lag_hr))
    if not stale:
        return
    import time as _t
    lines = [
        "=" * 70,
        "[Radia] WARNING: Cubit C++ plugin is out of date.",
        "-" * 70,
    ]
    for name, p, m, lag_hr in stale:
        lines.append(f"  {name}: {lag_hr:.1f} h older than Python package")
        lines.append(f"    {p}")
        lines.append(f"    plugin mtime: {_t.ctime(m)}")
    lines.append(f"  Python package ({os.path.basename(__file__)}):")
    lines.append(f"    mtime: {_t.ctime(self_mtime)}")
    lines.append("")
    lines.append("  RISK: silent wrong results (e.g., IH BEM source/sink")
    lines.append("        sidesets not read if FaceDescriptor DomainIn/Out")
    lines.append("        is stale).")
    lines.append("  FIX:  rebuild + redeploy the plugin:")
    lines.append("          cd src/cubit_plugin && cmake --build . --config Release")
    lines.append("          cubit-plugin-install   # from cubit-mesh-export pkg")
    lines.append("=" * 70)
    msg = "\n".join(lines)
    # Print to Cubit console (visible in the GUI log tab)
    try:
        print(msg)
    except Exception:
        pass
    # Panel log (timestamped, persistent across subprocess calls)
    for name, p, m, lag_hr in stale:
        _panel_log(
            f"PLUGIN_STALE: {name} is {lag_hr:.1f} h older than Python "
            f"package — rebuild and redeploy required")


def _cleanup_legacy_menus():
    """Remove legacy menus left behind by older installs.

    Targets:
      - Python-side "Solve" / "Radia-NGSolve" / "Generate Coil" /
        "Reload Panels" menus from very old radia (pre-4.x panel
        registry refactor)
      - The Qt5 C++ "Export Mesh" menu from radia <= 4.79.0
        (RadiaComp.cpp was deleted in 4.80.0; the menu may linger
        if an old .ccl is still loaded by Cubit)
      - The transient PySide6 "Export Mesh" / "Radia Export" menus
        previously injected directly into Cubit's QMenuBar
      - "RadiaToolBar" QToolBar from older toolbar prototypes

    The supported persistent surface is now Coreform's WorkflowToolbar,
    installed through Custom Toolbar Editor.  It is owned by Cubit and is not
    one of the QMenu/QToolBar objects removed here.
    """
    main_window = _find_main_window()
    if main_window is None:
        _panel_log("_cleanup_legacy_menus: main_window is None — skip")
        return
    menu_bar = main_window.menuBar()
    removed = []
    # The Claro-owned Export menu is deliberately NOT in this list: it is
    # registered through Cubit's Claro API and owned by Cubit.
    # Enumerating those QMenu/QAction objects from PySide6 deadlocks
    # Cubit, so it is removed via emclaro instead -- see
    # radia_export_menu.install_menu(), which calls remove_menu_items()
    # by component tag before re-registering.  "Export Mesh" below is the
    # historical Qt5 .ccl menu name and is unrelated to the new menu.
    legacy_names = (
        "Solve", "Radia-NGSolve", "Generate Coil", "Reload Panels",
        "Export Mesh",
    )
    for name in legacy_names:
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


def _install_radia_export_menu():
    """Install the Export menu through Cubit's official Claro API.

    Delegates to ``radia_export_menu.install_menu()``, which registers
    PyActions via ``emclaro.add_to_menu()``.  Cubit owns the resulting
    menu, so it survives the cold-start menu-bar rebuild that destroyed
    the previous PySide6 ``QMenuBar`` injection.
    """
    import radia_export_menu
    if radia_export_menu.install_menu():
        _panel_log("_install_radia_export_menu: installed via "
                   "emclaro.add_to_menu (Claro-owned, cold-start safe)")
    else:
        _panel_log("_install_radia_export_menu: Claro GUI unavailable "
                   "(batch mode?) - menu not installed")


def register_menu():
    """Initialize export defaults and install the Export menu.

    ``~/.cubit`` IS the correct hook for this.  MEASURED on Coreform
    Cubit 2025.12: ``emclaro.is_loaded()`` already returns True while
    this hook runs, and a menu registered here through Cubit's Claro API
    is still present after startup finishes.  The historical claim that
    a menu cannot be installed from this hook applied only to injecting
    a ``QMenu`` directly into Cubit's ``QMenuBar`` -- Cubit discarded
    those because it did not own them.
    """
    _panel_log("register_menu: ENTER (PySide6 dialogs + Claro menu API)")
    for step in (_cleanup_legacy_menus, _init_export_default_dir,
                 _check_plugin_freshness, _install_radia_export_menu):
        try:
            step()
        except Exception:
            _panel_log_exception(f"register_menu: {step.__name__} FAILED")
    _panel_log("register_menu: EXIT")


# Auto-register when this script is executed
try:
    register_menu()
except Exception:
    _panel_log_exception("register_menu top-level FAILED")
    raise
