"""
Register Radia menus in Coreform Cubit.

This script runs inside Cubit on startup (via ~/.cubit) and adds:
  - Radia-NGSolve: save .cub5 -> launch standalone PySide6 app
  - Reload Panels: re-read this script (debug)

Architecture:
  Cubit GUI (Python 3.10 + PySide6):
    1 button -> save .cub5 -> subprocess.run([python3.12, radia_app.py, path.cub5])

  Radia App (Python 3.12 + PySide6):
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
      1. RADIA_PYTHON env var
      2. py -3 launcher (Windows)
      3. python in PATH (if not Cubit's)
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
    """Find radia_app.py path without importing radia package.

    IMPORTANT: Do NOT use importlib.util.find_spec("radia") here.
    It triggers radia/__init__.py which loads _radia_pybind.pyd,
    causing DLL conflicts inside Cubit's process.
    """
    # Same directory as this file -> ../radia_app.py
    pkg_root = os.path.dirname(_this_dir)  # src/radia/ or site-packages/radia/
    app_path = os.path.join(pkg_root, "radia_app.py")
    if os.path.isfile(app_path):
        return app_path
    # Fallback: search site-packages directly (no import)
    for sp in sys.path:
        candidate = os.path.join(sp, "radia", "radia_app.py")
        if os.path.isfile(candidate):
            return candidate
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

    Returns the working directory, or None if user cancelled.
    """
    if _has_model():
        return _last_jou_dir[0] or os.getcwd()

    start_dir = _last_jou_dir[0] or os.getcwd()
    try:
        from PySide6.QtWidgets import QFileDialog
    except ImportError:
        from PyQt5.QtWidgets import QFileDialog
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
        print("WARNING: Journal file did not create any geometry.")
        return None

    return jou_dir


def _launch_radia_app():
    """Export .vol and launch Radia app.

    Pipeline:
      1. No model -> prompt for .jou, play it
      2. Export netgen .vol (C++ plugin, order 2)
      3. Save .cub5 alongside
      4. Launch radia_app.py with .vol path
    """
    try:
        if not _has_model():
            work_dir = _ensure_model()
            if work_dir is None:
                return
        else:
            work_dir = _last_jou_dir[0] or os.getcwd()

        # Set Cubit working directory
        work_dir = work_dir.replace("\\", "/")
        cubit.cmd(f'cd "{work_dir}"')

        # Export .vol via C++ plugin (Path A, includes curvedelements)
        vol_path = os.path.join(work_dir, "radia_model.vol").replace("\\", "/")
        cubit.cmd(f'export netgen "{vol_path}" order 2 overwrite')
        if not os.path.isfile(vol_path):
            print("ERROR: export netgen failed. Check blocks/sidesets.")
            return
        print(f"Exported: {vol_path}")

        # Also save .cub5 (for reference / Path B debugging)
        cub5_path = os.path.join(work_dir, "radia_model.cub5").replace("\\", "/")
        cubit.cmd(f'save cub5 "{cub5_path}" overwrite')

        ext_python = _find_external_python()
        if not ext_python:
            print("ERROR: Python 3.12 not found. Set RADIA_PYTHON env var.")
            return

        radia_app = _find_radia_app()
        if not radia_app:
            print("ERROR: radia_app.py not found. Is radia installed? (pip install radia)")
            return

        # Clean environment to avoid Qt5/Qt6 and MKL DLL conflicts
        env = os.environ.copy()
        for key in list(env.keys()):
            if "QT" in key.upper() or "PYSIDE" in key.upper():
                del env[key]

        cmd = [ext_python, radia_app, vol_path]
        print(f"Launching: {' '.join(cmd)}")
        subprocess.Popen(cmd, cwd=work_dir, env=env)

    except Exception as e:
        print(f"ERROR in _launch_radia_app: {e}")


def register_menu():
    """Register menus in Cubit's menu bar."""
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

    # --- Menu 2: Generate Coil (subprocess -> STEP -> import) ---
    def _generate_coil():
        try:
            ext_python = _find_external_python()
            if not ext_python:
                print("ERROR: Python 3.12 not found.")
                return

            work_dir = _last_jou_dir[0] or os.getcwd()
            work_dir = work_dir.replace("\\", "/")
            cubit.cmd(f'cd "{work_dir}"')

            step_path = os.path.join(work_dir, "coil.step").replace("\\", "/")

            # Find generate_coil.py
            gen_script = os.path.join(_this_dir, "generate_coil.py")
            if not os.path.isfile(gen_script):
                print(f"ERROR: {gen_script} not found")
                return

            # Build command — default racetrack for now
            cmd = [ext_python, gen_script, "--output", step_path]
            print(f"Generating coil: {' '.join(cmd)}")

            env = os.environ.copy()
            for key in list(env.keys()):
                if "QT" in key.upper() or "PYSIDE" in key.upper():
                    del env[key]

            result = subprocess.run(cmd, cwd=work_dir, env=env,
                                     capture_output=True, text=True, timeout=30)
            if result.returncode != 0:
                print(f"ERROR: {result.stderr[:500]}")
                return

            # Parse JSON result
            import json
            info = json.loads(result.stdout.strip().split('\n')[-1])
            print(f"Coil STEP: {info['output']} ({info['n_segments']} segments)")

            # Import STEP into Cubit
            cubit.cmd(f'import step "{step_path}" heal')
            print("Coil imported into Cubit.")

        except Exception as e:
            print(f"ERROR in _generate_coil: {e}")

    action_coil = QAction("Generate Coil", main_window)
    action_coil.setStatusTip(
        "Generate racetrack coil STEP and import into Cubit")
    action_coil.triggered.connect(_generate_coil)
    menu_bar.addAction(action_coil)

    # --- Menu 3: Reload Panels (debug, direct action) ---
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
