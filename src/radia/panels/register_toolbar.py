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
    from PySide6.QtWidgets import (QApplication, QMainWindow, QMenu,
                                   QMessageBox, QToolBar, QDialog,
                                   QVBoxLayout, QHBoxLayout, QLabel,
                                   QLineEdit, QPushButton, QFileDialog,
                                   QDialogButtonBox, QCheckBox,
                                   QPlainTextEdit, QComboBox, QGroupBox)
    from PySide6.QtGui import QAction
    _QT = "PySide6"
except ImportError:
    from PyQt5.QtWidgets import (QApplication, QMainWindow, QMenu,
                                 QMessageBox, QToolBar, QAction, QDialog,
                                 QVBoxLayout, QHBoxLayout, QLabel,
                                 QLineEdit, QPushButton, QFileDialog,
                                 QDialogButtonBox, QCheckBox,
                                 QPlainTextEdit, QComboBox, QGroupBox)
    _QT = "PyQt5"


def _no_window_kwargs():
    """Return subprocess kwargs to suppress console window on Windows."""
    if sys.platform == "win32":
        si = subprocess.STARTUPINFO()
        si.dwFlags |= subprocess.STARTF_USESHOWWINDOW
        si.wShowWindow = 0  # SW_HIDE
        return {"startupinfo": si,
                "creationflags": subprocess.CREATE_NO_WINDOW}
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
        subprocess.Popen(cmd, cwd=work_dir, env=env,
                         **_no_window_kwargs())

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
        "Save model as .cub5 and launch Radia standalone app (Python 3.12)")
    action_launch.triggered.connect(_launch_radia_app)
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

    # --- Sub 3: Kelvin Transform ---
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

                    # Sweep angle from symmetry planes
                    # X=0 and Y=0: 90 (quarter)
                    # X=0 or Y=0: 180 (half)
                    # neither: 360 (full)
                    n_xy = int(sym_x) + int(sym_y)
                    if n_xy == 2:
                        sweep_angle = 90
                    elif n_xy == 1:
                        sweep_angle = 180
                    else:
                        sweep_angle = 360

                    vols_before = set(
                        cubit.parse_cubit_list("volume", "all"))

                    # 1. Create sphere via arc sweep
                    #    (sweep creates seam curves for
                    #     copy mesh mapping)
                    #    Z=0 sym: quarter-arc (90 deg)
                    #    else: semicircle (180 deg)
                    cubit.cmd("create vertex 0 0 " + str(r))
                    v_top = cubit.get_last_id("vertex")
                    cubit.cmd("create vertex " + str(r) + " 0 0")
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

                    # Create blocks (skip if already exist)
                    existing_block_names = set()
                    for bid in cubit.parse_cubit_list(
                            "block", "all"):
                        existing_block_names.add(
                            cubit.get_exodus_entity_name(
                                "block", bid).lower())

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

                    # Create sidesets for periodic BC
                    existing_ss_names = set()
                    for sid in cubit.parse_cubit_list(
                            "sideset", "all"):
                        existing_ss_names.add(
                            cubit.get_exodus_entity_name(
                                "sideset", sid).lower())

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
        cubit.cmd("reset")
        startup = os.path.join(_this_dir, "startup.py").replace("\\", "/")
        cubit.cmd('play "' + startup + '"')
    action_reload = QAction("Reload Panels", main_window)
    action_reload.setStatusTip("Re-read register_toolbar.py from disk (debug)")
    action_reload.triggered.connect(_reload_panels)
    solve_menu.addAction(action_reload)

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
    else:
        print("Solve menu registered.")


# Auto-register when this script is executed
register_menu()
