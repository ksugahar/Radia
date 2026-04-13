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

        # --- Auto-Kelvin: add if not already present ---
        model_labels = _get_model_labels()
        if "kelvin" not in model_labels:
            try:
                _panel_log("AUTO-KELVIN: 'kelvin' block not found, "
                           "auto-adding Kelvin open boundary")
                panels_dir = os.path.dirname(
                    os.path.abspath(__file__))
                if panels_dir not in sys.path:
                    sys.path.insert(0, panels_dir)
                from add_kelvin import add_kelvin_cubit

                # Detect air sphere radius from Cubit geometry.
                # Find the air block, get its volumes, find the
                # largest unmerged surface -> that's the outer sphere.
                air_bid = None
                for bid in cubit.get_block_id_list():
                    n = cubit.get_exodus_entity_name("block", bid)
                    if n and n.lower() == "air":
                        air_bid = bid
                        break
                if air_bid is None:
                    print("WARNING: No 'air' block found. "
                          "Cannot auto-add Kelvin.")
                else:
                    air_vols = list(cubit.parse_cubit_list(
                        "volume", "in block %d" % air_bid))
                    # Find largest surface area -> outer boundary
                    best_sid, best_area = 0, 0
                    for vid in air_vols:
                        for sid in cubit.get_relatives(
                                "volume", vid, "surface"):
                            a = cubit.surface(sid).area()
                            if a > best_area:
                                best_area = a
                                best_sid = sid
                    if best_sid > 0:
                        import math as _m
                        # R = max vertex distance from origin
                        vids = cubit.get_relatives(
                            "surface", best_sid, "vertex")
                        R = max(
                            _m.sqrt(sum(
                                c**2 for c in cubit.vertex(v).coordinates()))
                            for v in vids)
                        # Detect symmetry from vertex positions
                        # (same logic as calc_common.detect_symmetry)
                        all_verts = set()
                        for vid in air_vols:
                            for v in cubit.get_relatives(
                                    "volume", vid, "vertex"):
                                all_verts.add(v)
                        sym = []
                        for axis, name in enumerate(["x", "y", "z"]):
                            coords = [cubit.vertex(v).coordinates()[axis]
                                       for v in all_verts]
                            if (min(coords) >= -1e-6
                                    and any(abs(c) < 1e-6
                                            for c in coords)):
                                sym.append(name)
                        _panel_log(
                            f"AUTO-KELVIN: R={R:.4f}, symmetry={sym}")
                        info = add_kelvin_cubit(
                            R=R, symmetry=sym)
                        ox, oy, oz = info["center"]
                        print(f"Auto-Kelvin: R={R:.4f}, "
                              f"offset=({ox:.3f}, {oy:.3f}, {oz:.3f}), "
                              f"symmetry={sym}")
            except Exception as e:
                _panel_log_exception("AUTO-KELVIN failed")
                print(f"WARNING: Auto-Kelvin failed: {e}")
                print("Proceeding without Kelvin "
                      "(Dirichlet truncation).")

        # --- Export .vol (detect Kelvin status before export) ---
        model_labels = _get_model_labels()
        has_kelvin = "kelvin" in model_labels
        _panel_log(f"VOL EXPORT: order={order}, "
                   f"has_kelvin={has_kelvin}, "
                   f"labels={sorted(model_labels)}")
        if has_kelvin:
            print(f"Kelvin open boundary detected (order {order})")
        else:
            print(f"No Kelvin domain -- Dirichlet truncation "
                  f"(order {order})")

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
