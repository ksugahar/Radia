"""
Register a single "Launch Radia App" button in Coreform Cubit.

This script runs inside Cubit on startup (via ~/.cubit) and adds
a toolbar button that saves the current model as .cub5 and launches
the standalone Radia app (Python 3.12 + tkinter).

Architecture:
  Cubit GUI (Python 3.10 + PySide6):
    1 button -> save .cub5 -> subprocess.run([python3.12, radia_app.py, path.cub5])

  Radia App (Python 3.12 + tkinter):
    All settings, computation, and result display
"""

import sys
import os
import subprocess

# Determine script location (__file__ is not defined when run via Cubit 'play')
try:
    _this_dir = os.path.dirname(os.path.abspath(__file__))
except NameError:
    import importlib
    try:
        _radia_spec = importlib.util.find_spec("radia")
        if _radia_spec and _radia_spec.origin:
            _this_dir = os.path.join(os.path.dirname(_radia_spec.origin), "panels")
        else:
            raise ImportError
    except (ImportError, AttributeError):
        import glob as _glob
        for _base in [
            os.path.join(sys.prefix, "Lib", "site-packages"),
            os.path.join(sys.prefix, "lib", "python*", "site-packages"),
        ]:
            for _candidate in _glob.glob(os.path.join(_base, "radia", "panels")):
                if os.path.isdir(_candidate):
                    _this_dir = _candidate
                    break
            else:
                continue
            break
        else:
            _this_dir = os.getcwd()

import cubit

# Qt bindings: prefer PySide6, fall back to PyQt5 (Cubit ships PyQt5)
try:
    from PySide6.QtWidgets import (QApplication, QMainWindow, QMessageBox,
                                   QToolBar)
    from PySide6.QtGui import QAction
except ImportError:
    from PyQt5.QtWidgets import (QApplication, QMainWindow, QMessageBox,
                                 QToolBar, QAction)


def _find_external_python():
    """Find external Python 3.12 (not Cubit's bundled Python 3.10).

    Search order:
      1. RADIA_PYTHON environment variable
      2. py -3 (Windows Python Launcher)
      3. python (from PATH, skip if Cubit's Python)
    """
    radia_py = os.environ.get("RADIA_PYTHON")
    if radia_py and os.path.isfile(radia_py):
        return radia_py

    if sys.platform == "win32":
        try:
            result = subprocess.run(
                ["py", "-3", "-c", "import sys; print(sys.executable)"],
                capture_output=True, text=True, timeout=5
            )
            if result.returncode == 0:
                py_path = result.stdout.strip()
                if os.path.isfile(py_path):
                    return py_path
        except Exception:
            pass

    try:
        result = subprocess.run(
            ["python", "-c", "import sys; print(sys.executable)"],
            capture_output=True, text=True, timeout=5
        )
        if result.returncode == 0:
            py_path = result.stdout.strip()
            if os.path.isfile(py_path) and "Cubit" not in py_path:
                return py_path
    except Exception:
        pass

    return None


def _find_main_window():
    """Find Cubit's main QMainWindow."""
    app = QApplication.instance()
    if app is None:
        return None
    best = None
    best_count = 0
    for widget in app.topLevelWidgets():
        if isinstance(widget, QMainWindow):
            count = len(widget.menuBar().actions())
            if count > best_count:
                best_count = count
                best = widget
    return best


def _find_radia_app():
    """Find radia_app.py path."""
    # Same directory as this file -> ../radia_app.py
    pkg_root = os.path.dirname(_this_dir)  # src/radia/
    app_path = os.path.join(pkg_root, "radia_app.py")
    if os.path.isfile(app_path):
        return app_path
    # Fallback: try importlib to find installed radia package
    try:
        import importlib.util
        spec = importlib.util.find_spec("radia")
        if spec and spec.origin:
            app_path = os.path.join(os.path.dirname(spec.origin), "radia_app.py")
            if os.path.isfile(app_path):
                return app_path
    except Exception:
        pass
    return None


def _has_model():
    """Check if Cubit has geometry loaded."""
    try:
        return cubit.get_volume_count() > 0 or cubit.get_surface_count() > 0
    except Exception:
        return False


# Remembers the directory of the last .jou file loaded
_last_jou_dir = [None]


def _ensure_model():
    """Ensure Cubit has a model. If empty, prompt for .jou file.

    Returns the .jou directory if a journal was loaded, or cwd if model
    already existed. Returns None if the user cancelled or .jou failed.
    """
    if _has_model():
        return _last_jou_dir[0] or os.getcwd()

    start_dir = _last_jou_dir[0] or os.getcwd()
    jou_path, _ = QFileDialog.getOpenFileName(
        None, "Select Journal File",
        start_dir,
        "Cubit Journal (*.jou);;All Files (*)"
    )
    if not jou_path:
        return None  # cancelled
    jou_path = jou_path.replace("\\", "/")
    jou_dir = os.path.dirname(jou_path)
    _last_jou_dir[0] = jou_dir
    cubit.cmd(f'play "{jou_path}"')

    if not _has_model():
        QMessageBox.warning(
            None, "Warning",
            "Journal file did not create any geometry."
        )
        return None

    return jou_dir


def _launch_radia_app():
    """Save current Cubit model as .cub5 and launch Radia app."""
    ext_python = _find_external_python()
    if not ext_python:
        QMessageBox.critical(
            None, "Error",
            "External Python 3.12 not found.\n\n"
            "Set RADIA_PYTHON environment variable or install Python 3.12."
        )
        return

    work_dir = _ensure_model()
    if work_dir is None:
        return

    # Save .cub5 in the working directory
    cub5_path = os.path.join(work_dir, "radia_cubit_model.cub5")
    cub5_path = cub5_path.replace("\\", "/")
    cubit.cmd(f'save cub5 "{cub5_path}" overwrite')

    # Launch standalone app (non-blocking)
    # Clean environment: remove Cubit's Qt paths to avoid PySide6 conflicts
    env = os.environ.copy()
    for key in list(env.keys()):
        if "QT" in key.upper() or "PYSIDE" in key.upper():
            del env[key]
    # Remove Cubit bin from PATH to avoid DLL conflicts
    cubit_bin = os.path.dirname(ext_python)  # not Cubit's bin
    cubit_install = os.environ.get("CUBIT_DIR", "")
    if cubit_install:
        env["PATH"] = ";".join(
            p for p in env.get("PATH", "").split(";")
            if "Cubit" not in p and "cubit" not in p
        )

    radia_app = _find_radia_app()
    if radia_app:
        cmd = [ext_python, radia_app, cub5_path]
    else:
        cmd = [ext_python, "-m", "radia.radia_app", cub5_path]
    subprocess.Popen(cmd, cwd=work_dir, env=env)



def register_menu():
    """Register three top-level menus in Cubit's menu bar."""
    main_window = _find_main_window()
    if main_window is None:
        return

    menu_bar = main_window.menuBar()

    # Remove existing menus (for reload)
    is_reload = False
    for name in ("Radia-NGSolve", "Export Mesh", "Reload Panels"):
        for action in list(menu_bar.actions()):
            if action.text().replace("&", "") == name:
                sub = action.menu()
                menu_bar.removeAction(action)
                if sub:
                    sub.deleteLater()
                is_reload = True

    # Remove old toolbar (cleanup from previous versions)
    for tb in main_window.findChildren(QToolBar, "RadiaToolBar"):
        main_window.removeToolBar(tb)
        tb.deleteLater()

    app = QApplication.instance()
    if app is not None:
        app.setQuitOnLastWindowClosed(False)

    # Export Mesh menu is provided by the C++ .ccl component (RadiaComp).
    # Do NOT register a Python Export Mesh menu here.

    # --- Menu 1: Radia-NGSolve (direct action, no submenu) ---
    action_launch = QAction("Radia-NGSolve", main_window)
    action_launch.setStatusTip(
        "Save model as .cub5 and launch Radia standalone app (Python 3.12)")
    action_launch.triggered.connect(_launch_radia_app)
    menu_bar.addAction(action_launch)

    # --- Menu 2: Reload Panels (debug, direct action) ---
    def _reload_panels():
        startup = os.path.join(_this_dir, "startup.py").replace("\\", "/")
        cubit.cmd(f'play "{startup}"')
    action_reload = QAction("Reload Panels", main_window)
    action_reload.setStatusTip("Re-read register_toolbar.py from disk (debug)")
    action_reload.triggered.connect(_reload_panels)
    menu_bar.addAction(action_reload)

    if is_reload:
        print("Radia-NGSolve menus re-registered.")
    else:
        print("Radia-NGSolve menus registered.")


# Auto-register when this script is executed
register_menu()
