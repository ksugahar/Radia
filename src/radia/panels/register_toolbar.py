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

import json
import sys
import os
import subprocess
import time
import traceback


# ============================================================
# Panel debug log: shared writer in radia.panels.panel_log.
# Cubit-side processes truncate the log on session start so one
# Cubit session = one continuous log file across all subprocesses.
# ============================================================
# panel_log.py is sibling to this file inside src/radia/panels/.
# Make sure that directory is on sys.path so we can import it
# even when Cubit's `play` runs us without a parent package context.
_my_dir = os.path.dirname(os.path.abspath(__file__)) if "__file__" in dir() \
    else os.getcwd()
if _my_dir not in sys.path:
    sys.path.insert(0, _my_dir)
from panel_log import (init_panel_log, panel_log as _panel_log,
                       panel_log_exception as _panel_log_exception,
                       PANEL_LOG_PATH as _PANEL_LOG_PATH)

# Truncate (Cubit session start) and tag this process as 'cubit'.
init_panel_log("cubit", truncate=True, banner=True)

# Log the *file mtime* of this script. The user can compare against
# the on-disk mtime via `verify-deploy` skill / `Solve > Verify Deploy`
# menu to confirm Cubit is using the freshly-edited code (vs. a stale
# Python cache or shadow install).
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
    from PySide6.QtWidgets import (QApplication, QMainWindow, QMenu,
                                   QMessageBox, QToolBar, QDialog,
                                   QVBoxLayout, QHBoxLayout, QFormLayout,
                                   QLabel, QLineEdit, QPushButton,
                                   QFileDialog, QDialogButtonBox, QCheckBox,
                                   QPlainTextEdit, QComboBox, QGroupBox,
                                   QWidget)
    from PySide6.QtGui import QAction
    _QT = "PySide6"
except ImportError:
    from PyQt5.QtWidgets import (QApplication, QMainWindow, QMenu,
                                 QMessageBox, QToolBar, QAction, QDialog,
                                 QVBoxLayout, QHBoxLayout, QFormLayout,
                                 QLabel, QLineEdit, QPushButton,
                                 QFileDialog, QDialogButtonBox, QCheckBox,
                                 QPlainTextEdit, QComboBox, QGroupBox,
                                 QWidget)
    _QT = "PyQt5"


def _no_window_kwargs():
    """Return subprocess kwargs to suppress console window on Windows.

    Only CREATE_NO_WINDOW is used (hides the console).
    Do NOT set SW_HIDE — it hides PySide6 GUI windows too.
    """
    if sys.platform == "win32":
        return {"creationflags": subprocess.CREATE_NO_WINDOW}
    return {}


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
                capture_output=True, text=True, timeout=5,
                **_no_window_kwargs()
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
            capture_output=True, text=True, timeout=5,
            **_no_window_kwargs()
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
    # In batch mode (-nographics), Qt uses QCoreApplication which
    # has no topLevelWidgets(). Skip menu registration silently.
    if not hasattr(app, 'topLevelWidgets'):
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


def _find_radia_script(name):
    """Find a radia/*.py script without importing radia package.

    IMPORTANT: Do NOT use importlib.util.find_spec("radia") here.
    It triggers radia/__init__.py which loads _radia_pybind.pyd,
    causing DLL conflicts inside Cubit's process.
    """
    pkg_root = os.path.dirname(_this_dir)  # src/radia/ or site-packages/radia/
    path = os.path.join(pkg_root, name)
    if os.path.isfile(path):
        return path
    for sp in sys.path:
        candidate = os.path.join(sp, "radia", name)
        if os.path.isfile(candidate):
            return candidate
    return None
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

    start_dir = _last_jou_dir[0] or _load_last_dir()
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
    _save_last_dir(jou_dir)
    cubit.cmd(f'play "{jou_path}"')

    if not _has_model():
        print("WARNING: Journal file did not create any geometry.")
        return None

    return jou_dir


def _discover_mode_scripts():
    """Discover radia_*.py analysis windows and their metadata.

    Convention: each radia_*.py in the radia package directory defines:
      TITLE = "Display Name"
      REQUIRED_LABELS = ["source", "sink"]    # must exist as block/sideset
      OPTIONAL_LABELS = ["workpiece", "air"]  # shown but not required
      OPTIONAL_FILES = {"Coil script": "Python (*.py)"}  # optional file inputs

    Returns dict {title: {"file": filename, "required": [...],
                          "optional": [...], "opt_files": {name: filter}}}
    sorted by filename.
    Does NOT import the modules (avoids DLL conflicts inside Cubit).
    """
    import ast
    import glob
    import re
    pkg_root = os.path.dirname(_this_dir)
    modes = {}
    for path in sorted(glob.glob(os.path.join(pkg_root, "radia_*.py"))):
        name = os.path.basename(path)
        if name == "radia_gui_base.py":
            continue
        title = None
        required = []
        optional = []
        opt_files = {}
        try:
            with open(path, "r", encoding="utf-8") as f:
                for line in f:
                    m = re.match(r'^TITLE\s*=\s*["\'](.+?)["\']', line)
                    if m:
                        title = m.group(1)
                    m = re.match(
                        r'^REQUIRED_LABELS\s*=\s*(\[.*\])', line)
                    if m:
                        required = ast.literal_eval(m.group(1))
                    m = re.match(
                        r'^OPTIONAL_LABELS\s*=\s*(\[.*\])', line)
                    if m:
                        optional = ast.literal_eval(m.group(1))
                    m = re.match(
                        r'^OPTIONAL_FILES\s*=\s*(\{.*\})', line)
                    if m:
                        opt_files = ast.literal_eval(m.group(1))
        except (OSError, ValueError, SyntaxError):
            pass
        if title:
            modes[title] = {
                "file": name,
                "required": required,
                "optional": optional,
                "opt_files": opt_files,
            }
    return modes


def _unc_to_drive(path):
    """Convert UNC path to mapped drive letter if possible.

    QFileDialog on Windows may not handle UNC paths correctly,
    defaulting to OneDrive instead. Convert to drive letter.
    """
    if not path or not path.startswith(("//", "\\\\")):
        return path
    norm = path.replace("/", "\\")
    try:
        import ctypes
        buf = ctypes.create_unicode_buffer(512)
        # Enumerate network drives
        for letter in "ABCDEFGHIJKLMNOPQRSTUVWXYZ":
            drive = letter + ":"
            size = ctypes.c_ulong(512)
            ret = ctypes.windll.mpr.WNetGetConnectionW(
                drive, buf, ctypes.byref(size))
            if ret == 0:  # NO_ERROR
                remote = buf.value
                if norm.lower().startswith(remote.lower()):
                    rest = norm[len(remote):]
                    return drive + rest.replace("\\", "/")
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


def _launcher_settings_path():
    """~/.radia/radia_launcher.json"""
    return os.path.join(os.path.expanduser("~"),
                        ".radia", "radia_launcher.json")


def _load_launcher_settings():
    """Load launcher settings from ~/.radia/radia_launcher.json."""
    p = _launcher_settings_path()
    if not os.path.isfile(p):
        return {}
    try:
        with open(p, "r", encoding="utf-8") as f:
            return json.load(f)
    except (OSError, json.JSONDecodeError):
        return {}


def _save_launcher_settings(data):
    """Save launcher settings to ~/.radia/radia_launcher.json."""
    p = _launcher_settings_path()
    os.makedirs(os.path.dirname(p), exist_ok=True)
    try:
        # Merge with existing settings
        existing = _load_launcher_settings()
        existing.update(data)
        with open(p, "w", encoding="utf-8") as f:
            json.dump(existing, f, indent=2)
    except OSError:
        pass


def _load_last_dir():
    """Load last working directory from settings. Fallback to samples dir."""
    default = _samples_dir()
    data = _load_launcher_settings()
    d = data.get("last_jou_dir", "")
    if d and os.path.isdir(d):
        return d
    parent = os.path.dirname(d)
    if parent and os.path.isdir(parent):
        return parent
    return default


def _save_last_dir(d):
    """Save last working directory to settings."""
    _save_launcher_settings({"last_jou_dir": d})


def _launch_radia_ngsolve():
    """Show mode/order/folder dialog, export .vol, launch analysis window.

    Pipeline:
      1. Dialog: mode (IH/EM/PCB) + mesh order + output folder
      2. radia_export netgen .vol (C++ plugin, user-chosen order)
      3. Launch the selected analysis window with .vol path
    """
    _panel_log("_launch_radia_ngsolve: ENTER")
    try:
        if not _has_model():
            work_dir = _ensure_model()
            if work_dir is None:
                return
        else:
            work_dir = _last_jou_dir[0] or _load_last_dir()

        # --- Dialog ---
        dlg = QDialog(_find_main_window())
        dlg.setWindowTitle("Radia-NGSolve")
        dlg.setMinimumWidth(450)
        layout = QVBoxLayout(dlg)

        # Mode (dynamically discovered from radia_*.py TITLE)
        mode_scripts = _discover_mode_scripts()
        if not mode_scripts:
            print("ERROR: No radia_*.py analysis windows found.")
            return

        # Restore previous selections
        prev = _load_launcher_settings()
        last_mode = prev.get("last_mode", "")
        last_order = prev.get("last_order", 2)

        h_mode = QHBoxLayout()
        h_mode.addWidget(QLabel("Analysis:"))
        mode_combo = QComboBox()
        mode_names = list(mode_scripts.keys())
        mode_combo.addItems(mode_names)
        if last_mode in mode_names:
            mode_combo.setCurrentText(last_mode)
        h_mode.addWidget(mode_combo)
        h_mode.addStretch()
        layout.addLayout(h_mode)

        # Order
        h_order = QHBoxLayout()
        h_order.addWidget(QLabel("Mesh order:"))
        order_combo = QComboBox()
        for p in range(1, 6):
            order_combo.addItem(str(p))
        order_combo.setCurrentIndex(max(0, min(4, last_order - 1)))
        h_order.addWidget(order_combo)
        h_order.addStretch()
        layout.addLayout(h_order)

        # --- Label check area ---
        label_group = QGroupBox("Labels (blocks / sidesets)")
        label_layout = QVBoxLayout(label_group)
        label_layout.setContentsMargins(8, 8, 8, 8)
        label_widget = QLabel("")
        label_widget.setWordWrap(True)
        label_layout.addWidget(label_widget)
        layout.addWidget(label_group)

        def _get_model_labels():
            """Get all block and sideset names from current Cubit model.

            ID enumeration: `get_block_id_list()` / `get_sideset_id_list()`.
            Name lookup:    `get_exodus_entity_name(kind, id)` for both.

            We use the single canonical name API
            (`get_exodus_entity_name`) instead of `get_block_name` /
            `get_sideset_name` because the latter pair is inconsistent
            in Cubit 2025.3 — `get_sideset_name` does not exist on this
            build, only `get_block_name` does.

            Heavily traced — every Cubit API call is logged so we can
            diagnose dialog launch failures and label detection issues.
            """
            _panel_log("_get_model_labels: ENTER")
            names = set()
            try:
                bids = cubit.get_block_id_list()
                _panel_log(
                    f"  get_block_id_list -> {type(bids).__name__} {list(bids)}")
            except Exception:
                _panel_log_exception("  get_block_id_list FAILED")
                raise
            for bid in bids:
                try:
                    n = cubit.get_exodus_entity_name("block", bid)
                    _panel_log(
                        f"  get_exodus_entity_name('block',{bid}) -> {n!r}")
                except Exception:
                    _panel_log_exception(
                        f"  get_exodus_entity_name('block',{bid}) FAILED")
                    raise
                if n:
                    names.add(n.lower())
            try:
                sids = cubit.get_sideset_id_list()
                _panel_log(
                    f"  get_sideset_id_list -> {type(sids).__name__} {list(sids)}")
            except Exception:
                _panel_log_exception("  get_sideset_id_list FAILED")
                raise
            for sid in sids:
                try:
                    n = cubit.get_exodus_entity_name("sideset", sid)
                    _panel_log(
                        f"  get_exodus_entity_name('sideset',{sid}) -> {n!r}")
                except Exception:
                    _panel_log_exception(
                        f"  get_exodus_entity_name('sideset',{sid}) FAILED")
                    raise
                if n:
                    names.add(n.lower())
            _panel_log(f"_get_model_labels: EXIT names={sorted(names)}")
            return names

        def _update_labels(_=None):
            _panel_log("_update_labels: ENTER")
            mode = mode_combo.currentText()
            info = mode_scripts.get(mode, {})
            req = info.get("required", [])
            opt = info.get("optional", [])
            _panel_log(f"  mode={mode!r} req={req} opt={opt}")
            try:
                model_labels = _get_model_labels()
            except Exception:
                _panel_log_exception("_update_labels: _get_model_labels failed")
                # Show the error in the dialog rather than crashing it
                label_widget.setText(
                    '<span style="color:red">'
                    f'ERROR reading model labels — see {_PANEL_LOG_PATH}'
                    '</span>')
                ok_btn.setEnabled(False)
                return

            lines = []
            all_ok = True
            for lbl in req:
                found = lbl.lower() in model_labels
                if found:
                    lines.append(
                        f'<span style="color:green">'
                        f'[OK] {lbl} (required)</span>')
                else:
                    lines.append(
                        f'<span style="color:red; font-weight:bold">'
                        f'[MISSING] {lbl} (required)</span>')
                    all_ok = False
            for lbl in opt:
                found = lbl.lower() in model_labels
                if found:
                    lines.append(
                        f'<span style="color:green">'
                        f'[OK] {lbl} (optional)</span>')
                else:
                    lines.append(
                        f'<span style="color:gray">'
                        f'[ - ] {lbl} (optional)</span>')

            if not req and not opt:
                lines.append(
                    '<span style="color:gray">'
                    'No label requirements</span>')

            label_widget.setText("<br>".join(lines))
            ok_btn.setEnabled(all_ok)

        mode_combo.currentTextChanged.connect(_update_labels)

        # --- Optional files (dynamic per mode) ---
        files_group = QGroupBox("Optional files")
        files_layout = QFormLayout(files_group)
        files_layout.setContentsMargins(8, 8, 8, 8)
        layout.addWidget(files_group)
        _file_widgets = {}  # name -> QLineEdit

        def _update_files(_=None):
            # Clear previous rows
            while files_layout.rowCount() > 0:
                files_layout.removeRow(0)
            _file_widgets.clear()

            mode = mode_combo.currentText()
            info = mode_scripts.get(mode, {})
            opt_files = info.get("opt_files", {})

            if not opt_files:
                files_group.setVisible(False)
                return
            files_group.setVisible(True)

            for fname, ffilter in opt_files.items():
                row = QHBoxLayout()
                le = QLineEdit()
                le.setPlaceholderText(ffilter)
                btn = QPushButton("...")
                btn.setFixedWidth(30)
                def _browse(le=le, ff=ffilter):
                    p, _ = QFileDialog.getOpenFileName(
                        dlg, f"Select {fname}", "", ff)
                    if p:
                        le.setText(p.replace("\\", "/"))
                btn.clicked.connect(_browse)
                row.addWidget(le)
                row.addWidget(btn)
                container = QWidget()
                container.setLayout(row)
                row.setContentsMargins(0, 0, 0, 0)
                files_layout.addRow(fname + ":", container)
                _file_widgets[fname] = le

        mode_combo.currentTextChanged.connect(_update_files)
        _update_files()

        # Output folder
        h_dir = QHBoxLayout()
        h_dir.addWidget(QLabel("Output folder:"))
        dir_edit = QLineEdit(work_dir.replace("\\", "/"))
        h_dir.addWidget(dir_edit, 1)
        btn_browse = QPushButton("Browse...")
        def _browse_dir():
            d = QFileDialog.getExistingDirectory(
                dlg, "Select output folder", dir_edit.text())
            if d:
                dir_edit.setText(d.replace("\\", "/"))
        btn_browse.clicked.connect(_browse_dir)
        h_dir.addWidget(btn_browse)
        layout.addLayout(h_dir)

        # Buttons
        buttons = QDialogButtonBox(
            QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        ok_btn = buttons.button(QDialogButtonBox.Ok)
        buttons.accepted.connect(dlg.accept)
        buttons.rejected.connect(dlg.reject)
        layout.addWidget(buttons)

        # Trigger initial label check
        _update_labels()

        if dlg.exec_() != QDialog.Accepted:
            return

        mode_name = mode_combo.currentText()
        script_name = mode_scripts[mode_name]["file"]
        order = order_combo.currentIndex() + 1
        out_dir = dir_edit.text().replace("\\", "/")

        # Save selections for next launch
        _save_launcher_settings({
            "last_mode": mode_name,
            "last_order": order,
            "last_jou_dir": out_dir,
        })

        # --- Export .vol ---
        cubit.cmd(f'cd "{out_dir}"')
        vol_name = "radia_model.vol"
        vol_path = out_dir + "/" + vol_name
        cubit.cmd(
            f'radia_export netgen "{vol_path}" order {order} overwrite')
        if not os.path.isfile(vol_path):
            print("ERROR: radia_export netgen failed. "
                  "Check blocks/sidesets.")
            return
        print(f"Exported: {vol_path} (order {order})")

        # --- Find and launch analysis window ---
        ext_python = _find_external_python()
        if not ext_python:
            print("ERROR: Python 3.12 not found. "
                  "Set RADIA_PYTHON env var.")
            return

        script_path = _find_radia_script(script_name)
        if not script_path:
            print(f"ERROR: {script_name} not found. "
                  "Is radia installed? (pip install radia)")
            return

        # Clean environment to avoid Qt5/Qt6 and MKL DLL conflicts
        env = os.environ.copy()
        for key in list(env.keys()):
            if "QT" in key.upper() or "PYSIDE" in key.upper():
                del env[key]

        cmd = [ext_python, script_path, vol_path]
        # Pass optional files as --key=value arguments
        for fname, le in _file_widgets.items():
            fpath = le.text().strip()
            if fpath:
                # Convert "Coil script" -> "--coil-script"
                arg_key = "--" + fname.lower().replace(
                    " ", "-").replace(".", "")
                cmd += [arg_key, fpath]
        print(f"Launching: {' '.join(cmd)}")
        subprocess.Popen(cmd, cwd=out_dir, env=env,
                         **_no_window_kwargs())

    except Exception as e:
        _panel_log_exception("_launch_radia_ngsolve")
        print(f"ERROR in _launch_radia_ngsolve: {e} "
              f"(full traceback in {_PANEL_LOG_PATH})")


def register_menu():
    """Register menus in Cubit's menu bar."""
    _panel_log("register_menu: ENTER")
    main_window = _find_main_window()
    if main_window is None:
        _panel_log("register_menu: main_window is None — abort")
        return

    menu_bar = main_window.menuBar()

    # Remove existing menus (for reload)
    is_reload = False
    # Export Mesh is owned by C++ .ccl -- do NOT remove it here
    for name in ("Solve", "Radia-NGSolve", "Generate Coil",
                 "Reload Panels"):
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

    # === Top-level "Solve" menu with sub-items ===
    solve_menu = QMenu("Solve", main_window)

    # --- Sub 1: Radia-NGSolve ---
    action_launch = QAction("Radia-NGSolve...", main_window)
    action_launch.setStatusTip(
        "Export .vol and launch analysis window (IH / EM / PCB)")
    action_launch.triggered.connect(_launch_radia_ngsolve)
    solve_menu.addAction(action_launch)

    # --- Sub 2: Generate Coil ---
    def _generate_coil():
        try:
            default_script = os.path.join(_this_dir, "generate_coil.py")
            work_dir = _last_jou_dir[0] or os.getcwd()
            ext_py = _find_external_python() or "python"

            dlg = QDialog(main_window)
            dlg.setWindowTitle("Generate Coil")
            dlg.setMinimumWidth(500)
            layout = QVBoxLayout(dlg)

            # Script path
            layout.addWidget(QLabel("Coil script (.py):"))
            h1 = QHBoxLayout()
            script_edit = QLineEdit(default_script)
            h1.addWidget(script_edit)
            btn_browse = QPushButton("...")
            btn_browse.setFixedWidth(30)
            def _browse_script():
                f, _ = QFileDialog.getOpenFileName(
                    dlg, "Select coil script",
                    os.path.dirname(script_edit.text()),
                    "Python files (*.py)")
                if f:
                    script_edit.setText(f)
            btn_browse.clicked.connect(_browse_script)
            h1.addWidget(btn_browse)
            layout.addLayout(h1)

            # Output STEP path
            layout.addWidget(QLabel("Output STEP:"))
            step_edit = QLineEdit(
                os.path.join(work_dir, "coil.step").replace("\\", "/"))
            layout.addWidget(step_edit)

            # Command (editable + Copy button)
            layout.addWidget(QLabel("Command:"))
            h_cmd = QHBoxLayout()
            cmd_edit = QLineEdit()
            def _update_cmd():
                cmd_edit.setText(
                    ext_py + ' "' + script_edit.text()
                    + '" --output "' + step_edit.text() + '"')
            script_edit.textChanged.connect(_update_cmd)
            step_edit.textChanged.connect(_update_cmd)
            _update_cmd()
            h_cmd.addWidget(cmd_edit)
            btn_copy = QPushButton("Copy")
            btn_copy.setFixedWidth(60)
            def _do_copy():
                QApplication.clipboard().setText(cmd_edit.text())
            btn_copy.clicked.connect(_do_copy)
            h_cmd.addWidget(btn_copy)
            layout.addLayout(h_cmd)

            # Buttons
            buttons = QDialogButtonBox(
                QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
            layout.addWidget(buttons)
            buttons.rejected.connect(dlg.reject)

            def _on_ok():
                dlg.accept()
                try:
                    gen_script = script_edit.text().strip()
                    step_path = step_edit.text().strip().replace("\\", "/")

                    if not os.path.isfile(gen_script):
                        print("ERROR: " + gen_script + " not found")
                        return

                    out_dir = os.path.dirname(step_path) or work_dir
                    cubit.cmd('cd "' + out_dir.replace("\\", "/") + '"')

                    cmd = [ext_py, gen_script, "--output", step_path]
                    print("Generating coil: " + " ".join(cmd))

                    env = os.environ.copy()
                    for key in list(env.keys()):
                        if "QT" in key.upper() or "PYSIDE" in key.upper():
                            del env[key]

                    result = subprocess.run(
                        cmd, cwd=out_dir, env=env,
                        capture_output=True, text=True, timeout=30,
                        **_no_window_kwargs())
                    if result.returncode != 0:
                        print("ERROR: " + result.stderr[:500])
                        return

                    import json
                    info = json.loads(
                        result.stdout.strip().split('\n')[-1])
                    print("Coil STEP: " + info['output']
                          + " (" + str(info['n_segments'])
                          + " segments)")

                    cubit.cmd('import step "' + step_path + '" heal')
                    print("Coil imported into Cubit.")

                except Exception as e:
                    print("ERROR in _generate_coil: " + str(e))

            buttons.accepted.connect(_on_ok)
            dlg.show()

        except Exception as e:
            import traceback
            traceback.print_exc()
            print("ERROR: Generate Coil dialog failed: " + str(e))

    action_coil = QAction("Generate Coil...", main_window)
    action_coil.setStatusTip(
        "Generate racetrack coil STEP and import into Cubit")
    action_coil.triggered.connect(_generate_coil)
    solve_menu.addAction(action_coil)

    # --- Sub 3: Analytic Function ---
    def _analytic_function():
        try:
            if not _has_model():
                work_dir = _ensure_model()
                if work_dir is None:
                    return
            else:
                work_dir = _last_jou_dir[0] or _load_last_dir()

            dlg = QDialog(_find_main_window())
            dlg.setWindowTitle("Analytic Function")
            dlg.setMinimumWidth(480)
            layout = QVBoxLayout(dlg)

            # Curve order
            h_order = QHBoxLayout()
            h_order.addWidget(QLabel("Curve order:"))
            order_combo = QComboBox()
            for p in range(1, 6):
                order_combo.addItem(str(p))
            order_combo.setCurrentIndex(1)  # default 2
            h_order.addWidget(order_combo)
            h_order.addStretch()
            layout.addLayout(h_order)

            # Expression
            layout.addWidget(QLabel("Expression (x, y, z):"))
            expr_edit = QLineEdit("sqrt(x*x + y*y + z*z)")
            expr_edit.setPlaceholderText(
                "e.g. sin(x)*cos(y), x*x+y*y+z*z, exp(-x*x)")
            layout.addWidget(expr_edit)

            # Field name
            h_name = QHBoxLayout()
            h_name.addWidget(QLabel("Field name:"))
            name_edit = QLineEdit("f")
            name_edit.setFixedWidth(150)
            h_name.addWidget(name_edit)
            h_name.addStretch()
            layout.addLayout(h_name)

            # Output folder
            h_dir = QHBoxLayout()
            h_dir.addWidget(QLabel("Output folder:"))
            dir_edit = QLineEdit(work_dir.replace("\\", "/"))
            h_dir.addWidget(dir_edit, 1)
            btn_browse = QPushButton("Browse...")
            def _browse_dir():
                d = QFileDialog.getExistingDirectory(
                    dlg, "Select output folder", dir_edit.text())
                if d:
                    dir_edit.setText(d.replace("\\", "/"))
            btn_browse.clicked.connect(_browse_dir)
            h_dir.addWidget(btn_browse)
            layout.addLayout(h_dir)

            # Buttons
            buttons = QDialogButtonBox(
                QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
            buttons.accepted.connect(dlg.accept)
            buttons.rejected.connect(dlg.reject)
            layout.addWidget(buttons)

            if dlg.exec_() != QDialog.Accepted:
                return

            order = order_combo.currentIndex() + 1
            expression = expr_edit.text().strip()
            field_name = name_edit.text().strip() or "f"
            out_dir = dir_edit.text().replace("\\", "/")

            if not expression:
                print("ERROR: No expression specified.")
                return

            # Export .vol
            cubit.cmd(f'cd "{out_dir}"')
            vol_path = out_dir + "/analytic_model.vol"
            cubit.cmd(
                f'radia_export netgen "{vol_path}" order {order} overwrite')
            if not os.path.isfile(vol_path):
                print("ERROR: radia_export netgen failed.")
                return
            print(f"Exported: {vol_path} (order {order})")
            _save_last_dir(out_dir)

            # Launch calc_analytic.py
            ext_python = _find_external_python()
            if not ext_python:
                print("ERROR: Python 3.12 not found.")
                return

            script = _find_radia_script(
                os.path.join("panels", "calc_analytic.py"))
            if not script:
                # Try relative to _this_dir
                script = os.path.join(_this_dir, "calc_analytic.py")
            if not os.path.isfile(script):
                print("ERROR: calc_analytic.py not found.")
                return

            cmd = [ext_python, script,
                   "--vol", vol_path,
                   "--expr", expression,
                   "--name", field_name]
            print(f"Running: {' '.join(cmd)}")

            import subprocess as _sp
            result = _sp.run(cmd, capture_output=True, text=True,
                             cwd=out_dir, timeout=120)
            if result.stdout:
                print(result.stdout.rstrip())
            if result.stderr:
                print(result.stderr.rstrip())

            # Parse JSON result, open Netgen GUI
            for line in reversed(result.stdout.split("\n")):
                line = line.strip()
                if line.startswith("{"):
                    try:
                        data = json.loads(line)
                        vp = data.get("vol_path", "")
                        if vp and os.path.isfile(vp):
                            print(f"Opening Netgen: {vp}")
                            _sp.Popen(
                                [ext_python, "-m", "netgen", vp],
                                creationflags=0x08000000)
                    except json.JSONDecodeError:
                        pass
                    break

        except Exception as e:
            print(f"ERROR in _analytic_function: {e}")

    action_analytic = QAction("Analytic Function...", main_window)
    action_analytic.setStatusTip(
        "Evaluate mathematical expression on mesh and visualize in GMSH")
    action_analytic.triggered.connect(_analytic_function)
    solve_menu.addAction(action_analytic)

    # --- Sub 4: Kelvin Transform ---
    def _kelvin_transform():
        try:
            dlg = QDialog(main_window)
            dlg.setWindowTitle("Kelvin Transform")
            dlg.setMinimumWidth(520)
            layout = QVBoxLayout(dlg)

            # Radius
            h_r = QHBoxLayout()
            h_r.addWidget(QLabel("Radius [m]:"))
            radius_edit = QLineEdit("0.06")
            radius_edit.setFixedWidth(100)
            h_r.addWidget(radius_edit)
            h_r.addStretch()
            layout.addLayout(h_r)

            # Offset direction + distance
            grp_off = QGroupBox("Exterior sphere offset")
            off_lay = QVBoxLayout(grp_off)
            h_dir = QHBoxLayout()
            h_dir.addWidget(QLabel("Direction:"))
            dir_combo = QComboBox()
            dir_combo.addItems(["+X", "-X", "+Y", "-Y", "+Z", "-Z"])
            h_dir.addWidget(dir_combo)
            h_dir.addWidget(QLabel("  Distance [m]:"))
            dist_edit = QLineEdit("0.15")
            dist_edit.setFixedWidth(80)
            h_dir.addWidget(dist_edit)
            h_dir.addStretch()
            off_lay.addLayout(h_dir)
            layout.addWidget(grp_off)

            # Symmetry planes
            layout.addWidget(QLabel("Symmetry planes (cut sphere):"))
            h_sym = QHBoxLayout()
            chk_x = QCheckBox("X = 0")
            chk_y = QCheckBox("Y = 0")
            chk_z = QCheckBox("Z = 0")
            h_sym.addWidget(chk_x)
            h_sym.addWidget(chk_y)
            h_sym.addWidget(chk_z)
            h_sym.addStretch()
            layout.addLayout(h_sym)

            # Command preview (multi-line) + Copy
            layout.addWidget(QLabel("Commands:"))
            cmd_edit = QPlainTextEdit()
            cmd_edit.setMaximumHeight(160)
            layout.addWidget(cmd_edit)
            h_copy = QHBoxLayout()
            h_copy.addStretch()
            btn_copy = QPushButton("Copy")
            btn_copy.setFixedWidth(60)
            def _do_copy():
                QApplication.clipboard().setText(cmd_edit.toPlainText())
            btn_copy.clicked.connect(_do_copy)
            h_copy.addWidget(btn_copy)
            layout.addLayout(h_copy)

            def _offset_xyz():
                d = dir_combo.currentText()
                v = dist_edit.text().strip()
                sign = 1 if d[0] == "+" else -1
                axis = d[1]
                if axis == "X":
                    return (str(sign) + "*" + v, "0", "0")
                elif axis == "Y":
                    return ("0", str(sign) + "*" + v, "0")
                else:
                    return ("0", "0", str(sign) + "*" + v)

            def _offset_floats():
                d = dir_combo.currentText()
                try:
                    v = float(dist_edit.text().strip())
                except ValueError:
                    v = 0.3
                sign = 1.0 if d[0] == "+" else -1.0
                axis = d[1]
                if axis == "X":
                    return (sign * v, 0.0, 0.0)
                elif axis == "Y":
                    return (0.0, sign * v, 0.0)
                else:
                    return (0.0, 0.0, sign * v)

            def _build_cmds():
                r = radius_edit.text().strip()
                ox, oy, oz = _offset_floats()
                sx, sy, sz = str(ox), str(oy), str(oz)
                sym_x = chk_x.isChecked()
                sym_y = chk_y.isChecked()
                sym_z = chk_z.isChecked()
                lines = []
                lines.append("# kelvin_int (semicircle sweep)")
                lines.append("create vertex 0 0 " + r)
                lines.append("create vertex " + r + " 0 0")
                lines.append("create vertex 0 0 -" + r)
                lines.append("create curve arc three vertex"
                             " {v_top} {v_mid} {v_bot}")
                lines.append("create curve vertex"
                             " {v_top} {v_bot}")
                lines.append("create surface curve"
                             " {arc} {line}")
                lines.append("sweep surface {s}"
                             " axis 0 0 0 0 0 1 angle 360")
                lines.append('volume {id} rename "kelvin_int"')
                lines.append("imprint volume all")
                lines.append("merge volume all")
                lines.append("")
                lines.append("# kelvin_ext (at "
                             + sx + ", " + sy + ", " + sz + ")")
                lines.append("volume {inner} copy move"
                             " x " + sx
                             + " y " + sy
                             + " z " + sz + " nomesh")
                lines.append('volume {id} rename "kelvin_ext"')
                cuts = []
                if sym_x:
                    cuts.append("xplane")
                if sym_y:
                    cuts.append("yplane")
                if sym_z:
                    cuts.append("zplane")
                for pl in cuts:
                    lines.append(
                        "webcut volume {outer} with plane " + pl
                        + " center location "
                        + sx + " " + sy + " " + sz)
                if cuts:
                    lines.append(
                        "# Delete negative-side pieces")
                lines.append("")
                lines.append("# Mesh + block/sideset")
                lines.append(
                    "volume {inner} scheme tetmesh")
                lines.append("mesh volume {inner}")
                lines.append(
                    "copy mesh from volume {inner}"
                    " to volume {outer}")
                lines.append("mesh volume {outer}")
                lines.append(
                    'block {N} add volume {outer}')
                lines.append(
                    'block {N} name "kelvin"')
                lines.append(
                    "# Sidesets: kelvin_int + kelvin_ext")
                return "\n".join(lines)

            def _update_cmd():
                cmd_edit.setPlainText(_build_cmds())

            for w in [radius_edit, dist_edit]:
                w.textChanged.connect(_update_cmd)
            dir_combo.currentIndexChanged.connect(lambda: _update_cmd())
            chk_x.toggled.connect(lambda: _update_cmd())
            chk_y.toggled.connect(lambda: _update_cmd())
            chk_z.toggled.connect(lambda: _update_cmd())
            _update_cmd()

            # State: track generated volume IDs
            _kelvin_state = {"inner": None, "outer": None}
            # Undo stack: list of actions to reverse
            _undo_stack = []

            # --- Buttons: Generate + Mesh + Close ---
            h_btns = QHBoxLayout()
            btn_gen = QPushButton("Generate")
            btn_gen.setFixedHeight(28)
            btn_mesh = QPushButton("Mesh")
            btn_mesh.setFixedHeight(28)
            btn_mesh.setEnabled(False)
            btn_undo = QPushButton("Undo")
            btn_undo.setFixedHeight(28)
            btn_close = QPushButton("Close")
            btn_close.setFixedHeight(28)
            h_btns.addWidget(btn_gen)
            h_btns.addWidget(btn_mesh)
            h_btns.addStretch()
            h_btns.addWidget(btn_undo)
            h_btns.addWidget(btn_close)
            layout.addLayout(h_btns)
            btn_close.clicked.connect(dlg.close)

            # Detect existing kelvin_int / kelvin_ext volumes
            def _detect_kelvin_volumes():
                inner = []
                outer = []
                for vid in cubit.parse_cubit_list("volume", "all"):
                    vname = cubit.get_entity_name("volume", vid)
                    if vname == "kelvin_int":
                        inner.append(vid)
                    elif vname == "kelvin_ext":
                        outer.append(vid)
                return inner, outer

            _det_inner, _det_outer = _detect_kelvin_volumes()
            if _det_inner and _det_outer:
                _kelvin_state["inner"] = _det_inner
                _kelvin_state["outer"] = _det_outer
                btn_mesh.setEnabled(True)
                print("Kelvin detected: inner="
                      + str(_det_inner) + " outer="
                      + str(_det_outer))

            def _on_generate():
                try:
                    r = float(radius_edit.text().strip())
                    ox, oy, oz = _offset_floats()
                    dist = (ox**2 + oy**2 + oz**2)**0.5
                    sym_x = chk_x.isChecked()
                    sym_y = chk_y.isChecked()
                    sym_z = chk_z.isChecked()

                    if dist < 2 * r:
                        print("ERROR: Offset distance ("
                              + str(dist) + ") must be >= 2*radius ("
                              + str(2*r) + ").")
                        return

                    # Sweep angle and seed direction from
                    # symmetry planes.
                    # Arc seed point lies on the equator;
                    # sweep rotates it around the z-axis.
                    #   X=0 + Y=0: seed on x, sweep 90
                    #   Y=0 only:  seed on x, sweep 180
                    #   X=0 only:  seed on y, sweep 180
                    #   neither:   seed on x, sweep 360
                    if sym_x and sym_y:
                        sweep_angle = 90
                        seed_dir = (r, 0, 0)
                    elif sym_y:
                        sweep_angle = 180
                        seed_dir = (r, 0, 0)
                    elif sym_x:
                        sweep_angle = 180
                        seed_dir = (0, r, 0)
                    else:
                        sweep_angle = 360
                        seed_dir = (r, 0, 0)

                    vols_before = set(
                        cubit.parse_cubit_list("volume", "all"))

                    # 1. Create sphere via arc sweep
                    #    (sweep creates seam curves for
                    #     copy mesh mapping)
                    #    Z=0 sym: quarter-arc (90 deg)
                    #    else: semicircle (180 deg)
                    cubit.cmd("create vertex 0 0 " + str(r))
                    v_top = cubit.get_last_id("vertex")
                    sx, sy, sz = seed_dir
                    cubit.cmd("create vertex "
                              + str(sx) + " "
                              + str(sy) + " "
                              + str(sz))
                    v_mid = cubit.get_last_id("vertex")

                    if sym_z:
                        # Quarter-arc: top pole to equator
                        cubit.cmd("create vertex 0 0 0")
                        v_org = cubit.get_last_id("vertex")
                        cubit.cmd(
                            "create curve arc center vertex "
                            + str(v_org) + " "
                            + str(v_top) + " " + str(v_mid))
                        arc_id = cubit.get_last_id("curve")
                        # Close with two lines (z-axis + x-axis)
                        cubit.cmd(
                            "create curve vertex "
                            + str(v_top) + " " + str(v_org))
                        l1 = cubit.get_last_id("curve")
                        cubit.cmd(
                            "create curve vertex "
                            + str(v_org) + " " + str(v_mid))
                        l2 = cubit.get_last_id("curve")
                        cubit.cmd(
                            "create surface curve "
                            + str(arc_id) + " "
                            + str(l1) + " " + str(l2))
                    else:
                        # Semicircle: top pole to bottom pole
                        cubit.cmd(
                            "create vertex 0 0 " + str(-r))
                        v_bot = cubit.get_last_id("vertex")
                        cubit.cmd(
                            "create curve arc three vertex "
                            + str(v_top) + " " + str(v_mid)
                            + " " + str(v_bot))
                        arc_id = cubit.get_last_id("curve")
                        cubit.cmd(
                            "create curve vertex "
                            + str(v_top) + " " + str(v_bot))
                        line_id = cubit.get_last_id("curve")
                        cubit.cmd(
                            "create surface curve "
                            + str(arc_id) + " "
                            + str(line_id))

                    surf_id = cubit.get_last_id("surface")

                    # Sweep around z-axis
                    cubit.cmd(
                        "sweep surface " + str(surf_id)
                        + " axis 0 0 0 0 0 1"
                        + " angle " + str(sweep_angle))

                    vols_after = set(
                        cubit.parse_cubit_list("volume", "all"))
                    inner = sorted(vols_after - vols_before)

                    if not inner:
                        print("ERROR: sphere creation failed")
                        return

                    for vid in inner:
                        cubit.cmd('volume ' + str(vid)
                                  + ' rename "kelvin_int"')

                    cubit.cmd("imprint volume all")
                    cubit.cmd("merge volume all")

                    # 2. Copy to create exterior sphere
                    inner_str = " ".join(
                        str(v) for v in inner)
                    cubit.cmd(
                        "volume " + inner_str
                        + " copy move"
                        + " x " + str(ox)
                        + " y " + str(oy)
                        + " z " + str(oz)
                        + " nomesh")
                    vols_after_copy = set(
                        cubit.parse_cubit_list("volume", "all"))
                    outer = sorted(
                        vols_after_copy - vols_after
                        - vols_before)
                    # vols_after may be stale after zcut/delete
                    # recalculate: outer = new vols not in inner
                    all_now = set(
                        cubit.parse_cubit_list("volume", "all"))
                    outer = sorted(
                        all_now - set(inner) - vols_before)

                    for vid in outer:
                        cubit.cmd('volume ' + str(vid)
                                  + ' rename "kelvin_ext"')

                    # Clipping plane to see interior
                    cubit.cmd(
                        "Graphics Clip on Plane"
                        " location 0 0 0"
                        " direction 0 1 0")
                    cubit.cmd("from 0 " + str(-r) + " 0")
                    cubit.cmd("zoom reset")
                    _kelvin_state["inner"] = inner
                    _kelvin_state["outer"] = outer
                    _kelvin_state["r"] = r
                    _kelvin_state["offset"] = (ox, oy, oz)
                    btn_mesh.setEnabled(True)

                    n_sym = sum([sym_x, sym_y, sym_z])
                    frac = ("1/" + str(2**n_sym)
                            if n_sym > 0 else "full")
                    print("")
                    print("=== Kelvin Generate ===")
                    print("  sweep: " + str(sweep_angle) + " deg")
                    print("  kelvin_int: volumes "
                          + str(inner))
                    print("  kelvin_ext: volumes "
                          + str(outer)
                          + " (" + frac + ")")

                    # Push to undo stack
                    _undo_stack.append({
                        "action": "generate",
                        "volumes": inner + outer,
                    })

                except Exception as e:
                    import traceback
                    traceback.print_exc()
                    print("ERROR: " + str(e))

            def _on_mesh():
                try:
                    inner_halves = _kelvin_state.get("inner")
                    outer_halves = _kelvin_state.get("outer")
                    r = _kelvin_state.get("r", 0.06)
                    if not inner_halves or not outer_halves:
                        print("ERROR: Generate first")
                        return

                    blocks_before = set(
                        cubit.parse_cubit_list("block", "all"))
                    ss_before = set(
                        cubit.parse_cubit_list("sideset", "all"))

                    mesh_size = r / 3.0

                    # Mesh interior halves
                    for hid in inner_halves:
                        cubit.cmd("volume " + str(hid)
                                  + " scheme tetmesh")
                        cubit.cmd("volume " + str(hid)
                                  + " size " + str(mesh_size))
                    cubit.cmd("mesh volume "
                              + " ".join(str(h) for h in inner_halves))

                    # Copy surface mesh from interior to exterior
                    inner_hemi_surfs = []
                    outer_hemi_surfs = []
                    for src_id, dst_id in zip(inner_halves,
                                              outer_halves):
                        # Find hemisphere surface (largest area)
                        src_surfs = list(
                            cubit.get_relatives("volume", src_id,
                                                "surface"))
                        dst_surfs = list(
                            cubit.get_relatives("volume", dst_id,
                                                "surface"))
                        src_hemi = max(src_surfs,
                            key=lambda s: cubit.get_surface_area(s))
                        dst_hemi = max(dst_surfs,
                            key=lambda s: cubit.get_surface_area(s))
                        inner_hemi_surfs.append(src_hemi)
                        outer_hemi_surfs.append(dst_hemi)

                        # Get curve + vertex for mapping
                        src_c = cubit.get_relatives(
                            "surface", src_hemi, "curve")[0]
                        dst_c = cubit.get_relatives(
                            "surface", dst_hemi, "curve")[0]
                        src_v = cubit.get_relatives(
                            "curve", src_c, "vertex")[0]
                        dst_v = cubit.get_relatives(
                            "curve", dst_c, "vertex")[0]

                        cubit.cmd(
                            "copy mesh surface " + str(src_hemi)
                            + " onto surface " + str(dst_hemi)
                            + " source curve " + str(src_c)
                            + " source vertex " + str(src_v)
                            + " target curve " + str(dst_c)
                            + " target vertex " + str(dst_v))

                    # Mesh exterior volumes (surface already copied)
                    for hid in outer_halves:
                        cubit.cmd("volume " + str(hid)
                                  + " scheme tetmesh")
                        cubit.cmd("volume " + str(hid)
                                  + " size " + str(mesh_size))
                    cubit.cmd("mesh volume "
                              + " ".join(str(h) for h in outer_halves))

                    # Create blocks (skip if already exist).
                    # Use get_block_name (NOT get_exodus_entity_name /
                    # get_entity_name) — the latter raises
                    # "Invalid list type for the get_mref_entity function"
                    # on some Cubit versions.
                    existing_block_names = set()
                    for bid in cubit.parse_cubit_list(
                            "block", "all"):
                        try:
                            n = cubit.get_block_name(bid) or ""
                        except Exception:
                            n = ""
                        existing_block_names.add(n.lower())

                    existing_blocks = set(
                        cubit.parse_cubit_list("block", "all"))
                    nb = (max(existing_blocks) + 1
                          if existing_blocks else 1)

                    if "kelvin_int" not in existing_block_names:
                        cubit.cmd(
                            "block " + str(nb) + " add volume "
                            + " ".join(
                                str(h) for h in inner_halves))
                        cubit.cmd(
                            'block ' + str(nb)
                            + ' name "kelvin_int"')
                        nb += 1

                    if "kelvin" not in existing_block_names:
                        cubit.cmd(
                            "block " + str(nb) + " add volume "
                            + " ".join(
                                str(h) for h in outer_halves))
                        cubit.cmd(
                            'block ' + str(nb)
                            + ' name "kelvin"')

                    # Create sidesets for periodic BC.
                    # Use get_sideset_name for the same reason as above.
                    existing_ss_names = set()
                    for sid in cubit.parse_cubit_list(
                            "sideset", "all"):
                        try:
                            n = cubit.get_sideset_name(sid) or ""
                        except Exception:
                            n = ""
                        existing_ss_names.add(n.lower())

                    existing_ss = set(
                        cubit.parse_cubit_list("sideset", "all"))
                    ns = (max(existing_ss) + 1
                          if existing_ss else 1)

                    if "kelvin_int" not in existing_ss_names:
                        cubit.cmd(
                            "sideset " + str(ns)
                            + " add surface "
                            + " ".join(
                                str(s) for s in inner_hemi_surfs))
                        cubit.cmd(
                            'sideset ' + str(ns)
                            + ' name "kelvin_int"')
                        ns += 1

                    if "kelvin_ext" not in existing_ss_names:
                        cubit.cmd(
                            "sideset " + str(ns)
                            + " add surface "
                            + " ".join(
                                str(s) for s in outer_hemi_surfs))
                        cubit.cmd(
                            'sideset ' + str(ns)
                            + ' name "kelvin_ext"')

                    print("")
                    print("=== Kelvin Mesh ===")
                    print("  kelvin_int " + str(inner_halves)
                          + ": meshed (tet)")
                    print("  kelvin_ext " + str(outer_halves)
                          + ": meshed (copy + tet)")
                    print("  blocks: kelvin_int, kelvin")
                    print("  sidesets: kelvin_int, kelvin_ext")

                    # Track created blocks/sidesets
                    new_blocks = sorted(
                        set(cubit.parse_cubit_list(
                            "block", "all")) - blocks_before)
                    new_ss = sorted(
                        set(cubit.parse_cubit_list(
                            "sideset", "all")) - ss_before)

                    _undo_stack.append({
                        "action": "mesh",
                        "volumes": inner_halves + outer_halves,
                        "blocks": new_blocks,
                        "sidesets": new_ss,
                    })

                except Exception as e:
                    import traceback
                    traceback.print_exc()
                    print("ERROR: " + str(e))

            def _on_undo():
                try:
                    if not _undo_stack:
                        print("Nothing to undo")
                        return
                    entry = _undo_stack.pop()
                    action = entry["action"]

                    if action == "mesh":
                        # Delete sidesets, blocks, then mesh
                        for sid in reversed(entry.get(
                                "sidesets", [])):
                            cubit.cmd("delete sideset "
                                      + str(sid))
                        for bid in reversed(entry.get(
                                "blocks", [])):
                            cubit.cmd("delete block "
                                      + str(bid))
                        vol_str = " ".join(
                            str(v) for v in entry["volumes"])
                        cubit.cmd("delete mesh volume "
                                  + vol_str + " propagate")
                        print("Undo: mesh deleted")

                    elif action == "generate":
                        vol_str = " ".join(
                            str(v) for v in entry["volumes"])
                        cubit.cmd("delete volume " + vol_str)
                        _kelvin_state["inner"] = None
                        _kelvin_state["outer"] = None
                        btn_mesh.setEnabled(False)
                        print("Undo: volumes deleted")

                    # Re-detect
                    inner, outer = _detect_kelvin_volumes()
                    if inner and outer:
                        _kelvin_state["inner"] = inner
                        _kelvin_state["outer"] = outer
                        btn_mesh.setEnabled(True)

                except Exception as e:
                    print("ERROR: undo failed: " + str(e))

            btn_gen.clicked.connect(_on_generate)
            btn_mesh.clicked.connect(_on_mesh)
            btn_undo.clicked.connect(_on_undo)

            # Ctrl+Z shortcut within dialog
            try:
                from PySide6.QtGui import QKeySequence, QShortcut
            except ImportError:
                from PyQt5.QtGui import QKeySequence
                from PyQt5.QtWidgets import QShortcut
            shortcut_undo = QShortcut(
                QKeySequence("Ctrl+Z"), dlg)
            shortcut_undo.activated.connect(_on_undo)

            dlg.show()

        except Exception as e:
            import traceback
            traceback.print_exc()
            print("ERROR: Kelvin dialog failed: " + str(e))

    action_kelvin = QAction("Kelvin Transform...", main_window)
    action_kelvin.setStatusTip(
        "Add Periodic Kelvin sphere (open boundary) to the model")
    action_kelvin.triggered.connect(_kelvin_transform)
    solve_menu.addAction(action_kelvin)

    # --- Sub 4: Reload Panels (debug) ---
    solve_menu.addSeparator()

    def _reload_panels():
        # Re-execute startup.py after a QTimer delay so the current
        # Qt event handler returns cleanly before the script re-runs.
        # Do NOT call cubit.cmd("reset") — it destroys the loaded
        # model and can crash Cubit's Python state.
        from PySide6.QtCore import QTimer
        startup = os.path.join(_this_dir, "startup.py").replace("\\", "/")
        _panel_log("Reload Panels: scheduling re-play via QTimer")
        QTimer.singleShot(200, lambda: cubit.cmd('play "' + startup + '"'))
    action_reload = QAction("Reload Panels", main_window)
    action_reload.setStatusTip("Re-read register_toolbar.py from disk (debug)")
    action_reload.triggered.connect(_reload_panels)
    solve_menu.addAction(action_reload)

    # --- Sub 5: Verify Deploy (debug) ---
    def _verify_deploy():
        """Show mtime and short hash of key Radia modules so the user
        can confirm Cubit is reading the freshly-edited source files.

        Output goes to THREE places:
          1. A QMessageBox dialog (so the user actually sees it)
          2. C:/radia_panel_log.txt (permanent record)
          3. Cubit's Python console via print() (sometimes silent)
        """
        import hashlib
        import datetime as _dt

        _panel_log("=" * 70)
        _panel_log("Verify Deploy: checking Radia module mtimes")

        # Files inside the radia package that the panel cares about
        # (from the editable install or wheel install).
        radia_root = os.path.dirname(os.path.dirname(_this_dir))
        radia_pkg = os.path.join(radia_root, "radia")
        if not os.path.isdir(radia_pkg):
            radia_pkg = _this_dir.rsplit(os.sep + "panels", 1)[0]

        candidates = [
            os.path.join(radia_pkg, "radia_ih.py"),
            os.path.join(radia_pkg, "radia_em.py"),
            os.path.join(radia_pkg, "radia_pcb.py"),
            os.path.join(radia_pkg, "radia_gui_base.py"),
            os.path.join(radia_pkg, "bem_coupled_solver.py"),
            os.path.join(radia_pkg, "bem_inductance.py"),
            os.path.join(_this_dir, "register_toolbar.py"),
            os.path.join(_this_dir, "calc_inductance.py"),
            os.path.join(_this_dir, "calc_fem_kelvin.py"),
            os.path.join(_this_dir, "calc_heating_bem.py"),
        ]

        report_lines = ["Radia source-file deployment report", ""]
        for path in candidates:
            if not os.path.isfile(path):
                line = f"  [MISSING] {os.path.basename(path)}"
                report_lines.append(line)
                _panel_log(line)
                continue
            mtime = os.path.getmtime(path)
            ts = _dt.datetime.fromtimestamp(mtime).strftime(
                "%Y-%m-%d %H:%M:%S")
            with open(path, "rb") as fh:
                h = hashlib.sha1(fh.read()).hexdigest()[:8]
            short = os.path.basename(path)
            line = f"  {short:<26} {ts}  {h}"
            report_lines.append(line)
            _panel_log(line)

        # Also dump the editable install pointer if any
        report_lines.append("")
        try:
            import radia
            ver_line = f"  radia.__version__ = {radia.__version__}"
            file_line = f"  radia.__file__    = {radia.__file__}"
            report_lines.append(ver_line)
            report_lines.append(file_line)
            _panel_log(ver_line)
            _panel_log(file_line)
        except Exception as _e:
            err_line = f"  radia import FAILED: {_e}"
            report_lines.append(err_line)
            _panel_log(err_line)

        report_lines.append("")
        report_lines.append(
            "Full log: C:/radia_panel_log.txt")

        _panel_log(
            "Verify Deploy: done. See C:/radia_panel_log.txt for full log")

        # The actual user-visible output: a QMessageBox. Without this
        # the menu action looks like a no-op (the report only goes to
        # the log file, which the user has no reason to be tailing).
        try:
            from PySide6.QtWidgets import QMessageBox
        except ImportError:
            try:
                from PyQt5.QtWidgets import QMessageBox
            except ImportError:
                QMessageBox = None
        if QMessageBox is not None:
            box = QMessageBox(main_window)
            box.setWindowTitle("Radia: Verify Deploy")
            box.setIcon(QMessageBox.Information)
            box.setText("Source-file deployment report")
            box.setInformativeText(
                "Compare the timestamps below against your most "
                "recent edits to confirm Cubit is reading the latest "
                "files. If a row shows an old time, restart Cubit "
                "(or rerun the deploy skill).")
            box.setDetailedText("\n".join(report_lines))
            # Make the box wide enough that the monospace report
            # does not get truncated.
            box.setStyleSheet("QLabel{min-width: 600px;}")
            box.exec()
        else:
            print("\n".join(report_lines))

    action_verify = QAction("Verify Deploy", main_window)
    action_verify.setStatusTip(
        "Print mtime + sha1 of key Radia modules — confirms Cubit is "
        "reading the latest edited files")
    action_verify.triggered.connect(_verify_deploy)
    solve_menu.addAction(action_verify)

    # Insert before Help menu (last built-in menu)
    help_action = None
    for action in menu_bar.actions():
        if action.text().replace("&", "") == "Help":
            help_action = action
            break
    if help_action:
        menu_bar.insertMenu(help_action, solve_menu)
    else:
        menu_bar.addMenu(solve_menu)

    if is_reload:
        print("Solve menu re-registered.")
        _panel_log("register_menu: re-registered Solve menu")
    else:
        print("Solve menu registered.")
        _panel_log("register_menu: registered Solve menu")

    # Write default_dir to C++ export_settings.json so Export Mesh
    # dialogs default to samples/ instead of OneDrive/CWD
    _init_export_default_dir()
    _panel_log("register_menu: EXIT (success)")


# Auto-register when this script is executed
try:
    register_menu()
except Exception:
    _panel_log_exception("register_menu top-level FAILED")
    raise
