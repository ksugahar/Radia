"""
Shared base for Radia standalone analysis windows (PySide6).

Provides:
  - ModePanel: form layout with add_line/add_combo/add_spin/add_browse helpers
  - AnalysisWindow: QMainWindow with Run/Stop/Save/GMSH buttons + QProcess
"""

import json
import os
import sys

from PySide6.QtCore import Qt, QProcess
from PySide6.QtGui import QFont, QIcon
from PySide6.QtWidgets import (
    QApplication, QMainWindow, QWidget,
    QVBoxLayout, QHBoxLayout, QFormLayout,
    QLabel, QLineEdit, QComboBox, QSpinBox,
    QPushButton, QPlainTextEdit, QFileDialog,
    QMessageBox, QSplitter, QGroupBox, QStyle,
)

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_PANELS_DIR = os.path.join(_THIS_DIR, "panels")
_RESOURCES_DIR = os.path.join(_THIS_DIR, "resources")
_PYTHON = sys.executable
_SETTINGS_DIR = os.path.join(os.path.expanduser("~"), ".radia")

# Shared panel debug log (C:/radia_panel_log.txt on Windows). Append-only
# from this process; the Cubit-side register_toolbar.py truncates it on
# session start so the file holds one continuous Cubit session log.
if _PANELS_DIR not in sys.path:
    sys.path.insert(0, _PANELS_DIR)
try:
    from panel_log import (init_panel_log, panel_log,
                           panel_log_exception, PANEL_LOG_PATH)
    init_panel_log("ih-window", truncate=False, banner=True)
    panel_log(f"radia_gui_base.py file={__file__}")
except Exception:
    # Panel log is best-effort. Define no-op shims so the analysis
    # windows still run if panel_log.py is missing for any reason.
    PANEL_LOG_PATH = ""
    def panel_log(_msg):
        pass
    def panel_log_exception(_prefix=""):
        pass


def _icon_path():
    for ext in (".ico", ".png"):
        p = os.path.join(_RESOURCES_DIR, "radia_icon" + ext)
        if os.path.isfile(p):
            return p
    return ""


def calc_script(name):
    return os.path.join(_PANELS_DIR, name)


def msh_output(vol_path, suffix):
    base = vol_path if vol_path else "output"
    return os.path.splitext(base)[0] + suffix + ".msh"


# ============================================================
# ModePanel: form layout base
# ============================================================

class ModePanel(QWidget):
    """Base class for analysis parameter panels."""

    _RED = "QLineEdit { background-color: #FFD0D0; }"
    _NORMAL = ""

    def __init__(self, parent=None):
        super().__init__(parent)
        self._form = QFormLayout()
        self._form.setLabelAlignment(Qt.AlignRight)
        self.setLayout(self._form)
        self._widgets = {}
        self._row_indices = {}

    def _set_row_visible(self, key, visible):
        row_idx = self._row_indices.get(key)
        if row_idx is None:
            return
        label_item = self._form.itemAt(row_idx, QFormLayout.LabelRole)
        field_item = self._form.itemAt(row_idx, QFormLayout.FieldRole)
        if label_item and label_item.widget():
            label_item.widget().setVisible(visible)
        if field_item and field_item.widget():
            field_item.widget().setVisible(visible)

    def add_line(self, key, label, default="", placeholder=""):
        w = QLineEdit(default)
        if placeholder:
            w.setPlaceholderText(placeholder)
        self._form.addRow(label, w)
        self._widgets[key] = w
        self._row_indices[key] = self._form.rowCount() - 1
        return w

    def add_combo(self, key, label, items, default=0):
        w = QComboBox()
        w.addItems(items)
        w.setCurrentIndex(default)
        self._form.addRow(label, w)
        self._widgets[key] = w
        self._row_indices[key] = self._form.rowCount() - 1
        return w

    def add_spin(self, key, label, value=1, lo=1, hi=999):
        w = QSpinBox()
        w.setRange(lo, hi)
        w.setValue(value)
        self._form.addRow(label, w)
        self._widgets[key] = w
        self._row_indices[key] = self._form.rowCount() - 1
        return w

    def add_browse(self, key, label, default="", filter_str="All files (*.*)"):
        row = QHBoxLayout()
        le = QLineEdit(default)
        btn = QPushButton("...")
        btn.setFixedWidth(30)
        btn.clicked.connect(lambda: self._do_browse(le, filter_str))
        row.addWidget(le)
        row.addWidget(btn)
        container = QWidget()
        container.setLayout(row)
        row.setContentsMargins(0, 0, 0, 0)
        self._form.addRow(label, container)
        self._widgets[key] = le
        self._row_indices[key] = self._form.rowCount() - 1
        return le

    def _do_browse(self, line_edit, filter_str):
        path, _ = QFileDialog.getOpenFileName(self, "Select file", "", filter_str)
        if path:
            line_edit.setText(path)

    def val(self, key):
        w = self._widgets[key]
        if isinstance(w, QLineEdit):
            return w.text().strip()
        elif isinstance(w, QComboBox):
            return w.currentText()
        elif isinstance(w, QSpinBox):
            return str(w.value())
        return ""

    def save_state(self):
        state = {}
        for key, w in self._widgets.items():
            if isinstance(w, QLineEdit):
                state[key] = w.text()
            elif isinstance(w, QComboBox):
                # Save the SELECTED TEXT, not the index. Combo items can
                # be added / removed / reordered between releases (e.g.
                # the IH panel dropping "BEM-SIBC (WP)" 2026-04-12), and
                # a stale index would silently jump to the wrong item or
                # land out-of-range and leave the combo blank.
                state[key] = w.currentText()
            elif isinstance(w, QSpinBox):
                state[key] = w.value()
        return state

    def restore_state(self, state):
        if not state:
            return
        for key, val in state.items():
            w = self._widgets.get(key)
            if w is None:
                continue
            try:
                if isinstance(w, QLineEdit):
                    w.setText(str(val))
                elif isinstance(w, QComboBox):
                    # Two state encodings need to be supported:
                    #  - text   (new format, robust to combo edits)
                    #  - index  (legacy format from before 2026-04-12)
                    # If neither matches a current item we leave the
                    # combo at its panel-default selection rather than
                    # forcing an out-of-range index that blanks the
                    # widget.
                    if isinstance(val, str):
                        idx = w.findText(val)
                        if idx >= 0:
                            w.setCurrentIndex(idx)
                    else:
                        ival = int(val)
                        if 0 <= ival < w.count():
                            w.setCurrentIndex(ival)
                elif isinstance(w, QSpinBox):
                    w.setValue(int(val))
            except (ValueError, TypeError):
                pass

    def is_runnable(self):
        return True

    def build_command(self, vol_path):
        raise NotImplementedError


# ============================================================
# AnalysisWindow: main window with Run/Stop/GMSH
# ============================================================

class AnalysisWindow(QMainWindow):
    """Standalone analysis window with .vol input, parameter panel, Run/Stop."""

    def __init__(self, title, vol_path="", settings_key="default"):
        super().__init__()
        self.setWindowTitle(title)
        self.resize(650, 600)
        self.setMinimumSize(500, 400)
        self._settings_key = settings_key

        icon = _icon_path()
        if icon:
            self.setWindowIcon(QIcon(icon))

        self._process = None
        self._last_msh = None
        self._panel = None  # set by subclass via _set_panel()
        self._vol_path = vol_path

        self._build_ui(vol_path)

    def _set_panel(self, panel):
        """Set the parameter panel (call from subclass __init__)."""
        self._panel = panel
        self._panel_area.addWidget(panel)
        panel.validationChanged = self._update_run_state
        self._update_run_state()

    def _build_ui(self, vol_path):
        central = QWidget()
        self.setCentralWidget(central)
        root = QVBoxLayout(central)
        root.setContentsMargins(10, 10, 10, 5)

        # Model path
        model_row = QHBoxLayout()
        model_row.addWidget(QLabel("Model (.vol):"))
        self._vol_edit = QLineEdit(vol_path)
        self._vol_edit.setPlaceholderText(
            ".vol file (exported from Cubit or Netgen)")
        model_row.addWidget(self._vol_edit, 1)
        browse_btn = QPushButton("Browse...")
        browse_btn.clicked.connect(self._browse_vol)
        model_row.addWidget(browse_btn)
        root.addLayout(model_row)

        # Splitter: panel | output
        splitter = QSplitter(Qt.Vertical)
        root.addWidget(splitter, 1)

        # Panel area (subclass inserts panel here)
        self._panel_container = QWidget()
        self._panel_area = QVBoxLayout(self._panel_container)
        self._panel_area.setContentsMargins(0, 0, 0, 0)
        splitter.addWidget(self._panel_container)

        # Output
        out_group = QGroupBox("Output")
        out_layout = QVBoxLayout(out_group)
        out_layout.setContentsMargins(5, 5, 5, 5)
        self._output = QPlainTextEdit()
        self._output.setReadOnly(True)
        self._output.setFont(QFont("Consolas", 9))
        out_layout.addWidget(self._output)
        splitter.addWidget(out_group)
        splitter.setSizes([300, 250])

        # Buttons
        style = self.style()
        btn_row = QHBoxLayout()

        self._run_btn = QPushButton(
            style.standardIcon(QStyle.SP_MediaPlay), " Run")
        self._run_btn.setFixedHeight(32)
        self._run_btn.clicked.connect(self._on_run)
        btn_row.addWidget(self._run_btn)

        self._stop_btn = QPushButton(
            style.standardIcon(QStyle.SP_MediaStop), " Stop")
        self._stop_btn.setFixedHeight(32)
        self._stop_btn.setEnabled(False)
        self._stop_btn.clicked.connect(self._on_stop)
        btn_row.addWidget(self._stop_btn)

        self._gmsh_btn = QPushButton(
            style.standardIcon(QStyle.SP_ComputerIcon), " Open GMSH")
        self._gmsh_btn.setFixedHeight(32)
        self._gmsh_btn.setEnabled(False)
        self._gmsh_btn.clicked.connect(self._open_gmsh)
        btn_row.addWidget(self._gmsh_btn)

        btn_row.addStretch()
        root.addLayout(btn_row)

        # Status bar
        self._status = self.statusBar()

    def _browse_vol(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Select .vol file",
            os.path.dirname(self._vol_edit.text()),
            "Netgen Vol (*.vol);;All (*)")
        if path:
            self._vol_edit.setText(path)

    def _update_run_state(self):
        if self._process is not None:
            return
        runnable = self._panel.is_runnable() if self._panel else True
        self._run_btn.setEnabled(runnable)

    def _on_run(self):
        if self._process is not None:
            return
        vol = self._vol_edit.text().strip()
        try:
            cmd = self._panel.build_command(vol)
        except ValueError as e:
            panel_log(f"_on_run: build_command FAILED: {e}")
            QMessageBox.warning(self, "Input Error", str(e))
            return

        panel_log(f"_on_run: vol={vol}")
        panel_log(f"  cmd: {' '.join(cmd)}")

        self._save_settings()
        self._output.clear()
        self._output.appendPlainText(f"> {' '.join(cmd)}\n")
        self._last_msh = None
        self._gmsh_btn.setEnabled(False)

        work_dir = os.path.dirname(vol) if vol else os.getcwd()
        if not os.path.isdir(work_dir):
            work_dir = os.getcwd()

        # Cleanup stale GMSH artifacts left by previous Runs from the
        # working directory. Without this the user can accidentally
        # open an old `inductance.geo` (full B-field box from a
        # previous Air=on run) when the current Run is Air=off.
        # We delete files the panel itself produces; the user's own
        # .msh / .geo files (anything not in the known list) are left
        # alone.
        for stale in ("inductance.geo", "inductance_B.msh",
                      "inductance_J.msh", "J_coeffs.npy",
                      "surface_mesh.vol",
                      os.path.basename(msh_output(vol, "_bem"))
                      if vol else ""):
            if not stale:
                continue
            p = os.path.join(work_dir, stale)
            if os.path.isfile(p):
                try:
                    os.remove(p)
                    panel_log(f"  cleaned stale: {stale}")
                except OSError:
                    pass

        self._process = QProcess(self)
        self._process.setWorkingDirectory(work_dir)
        self._process.readyReadStandardOutput.connect(self._read_stdout)
        self._process.readyReadStandardError.connect(self._read_stderr)
        self._process.finished.connect(self._on_finished)
        self._process.start(_PYTHON, cmd[1:])

        self._run_btn.setEnabled(False)
        self._stop_btn.setEnabled(True)
        self._status.showMessage("Running...")

    def _on_stop(self):
        if self._process:
            self._process.kill()
            self._process.waitForFinished(3000)

    def _read_stdout(self):
        data = self._process.readAllStandardOutput().data()
        text = data.decode("utf-8", errors="replace").rstrip("\n")
        if text:
            self._output.appendPlainText(text)

    def _read_stderr(self):
        data = self._process.readAllStandardError().data()
        text = data.decode("utf-8", errors="replace").rstrip("\n")
        if text:
            self._output.appendPlainText(text)

    def _on_finished(self, exit_code, exit_status):
        remaining = self._process.readAllStandardOutput().data()
        if remaining:
            self._output.appendPlainText(
                remaining.decode("utf-8", errors="replace").rstrip("\n"))

        stdout_text = self._output.toPlainText()
        self._process = None
        self._stop_btn.setEnabled(False)
        self._update_run_state()

        panel_log(f"_on_finished: exit_code={exit_code}")

        if exit_code != 0:
            self._status.showMessage(f"Error (exit code {exit_code})")
            self._output.appendPlainText(
                f"\n*** Process exited with code {exit_code}")
            # Log the last ~20 lines so the failure is visible in the
            # panel debug log without the user needing to copy/paste.
            tail = "\n".join(stdout_text.splitlines()[-20:])
            panel_log(f"_on_finished: subprocess FAILED, tail:\n{tail}")
            return

        # Find .msh in output for GMSH button
        result = None
        for line in reversed(stdout_text.split("\n")):
            line = line.strip()
            if line.startswith("{"):
                try:
                    result = json.loads(line)
                except json.JSONDecodeError:
                    pass
                break
        # Prefer the merged .geo (B + J + companion) when air field
        # post-processing ran; otherwise fall back to the BEM .msh.
        gmsh_target = None
        gmsh_reason = ""
        if result and result.get("field_gmsh_file") and os.path.isfile(
                result["field_gmsh_file"]):
            gmsh_target = result["field_gmsh_file"]
            gmsh_reason = "field_gmsh_file (.geo with B + J + air-field)"
        elif result:
            msh_key = "msh_output" if "msh_output" in result else "msh_file"
            msh = result.get(msh_key)
            if msh and os.path.isfile(msh):
                gmsh_target = msh
                gmsh_reason = f"{msh_key}"
            else:
                gmsh_reason = (
                    f"no {msh_key} (value={msh!r}) — Open GMSH disabled")
        if gmsh_target:
            self._last_msh = gmsh_target
            self._gmsh_btn.setEnabled(True)
        panel_log(f"_on_finished: gmsh button = "
                  f"{'ENABLED' if gmsh_target else 'disabled'} "
                  f"({gmsh_reason})")
        if gmsh_target:
            panel_log(f"  -> {gmsh_target}")

        # Log result keys + headline numbers to panel debug log
        if result is not None:
            panel_log(f"_on_finished: result keys = "
                      f"{sorted(result.keys())[:20]}")
            if "error" in result:
                panel_log(f"  ERROR: {result['error']}")
            if "inductance_H" in result:
                panel_log(
                    f"  L = {result['inductance_H']*1e9:.3f} nH")
            if "coupled_dL_H" in result:
                panel_log(
                    f"  delta_L = {result['coupled_dL_H']*1e9:+.3f} nH")
            if "P_total_W" in result:
                panel_log(f"  P_total = {result['P_total_W']:.4e} W")
            if "wp_P_total" in result:
                panel_log(f"  wp_P_total = {result['wp_P_total']:.4e} W")

        # Display result summary in output window
        if result and "error" not in result:
            self._output.appendPlainText("\n--- Result ---")
            if "inductance_H" in result:
                L_nH = result["inductance_H"] * 1e9
                self._output.appendPlainText(f"  L = {L_nH:.2f} nH")
            if "n_dofs" in result:
                self._output.appendPlainText(f"  DOFs = {result['n_dofs']}")
            if "t_solve" in result:
                self._output.appendPlainText(
                    f"  Time: {result['t_solve']:.1f}s total")
            if "t_assembly" in result:
                self._output.appendPlainText(
                    f"    Assembly: {result['t_assembly']:.1f}s")
            if "t_lu" in result:
                self._output.appendPlainText(
                    f"    LU solve: {result['t_lu']:.1f}s")
            # BEM-SIBC workpiece results (calc_heating_bem)
            if "P_total_W" in result:
                self._output.appendPlainText(
                    "  --- Coil ---")
                self._output.appendPlainText(
                    f"  Coil radius:  "
                    f"{result.get('coil_radius_m', 0)*1e3:.1f} mm")
                self._output.appendPlainText(
                    f"  Coil current: "
                    f"{result.get('coil_current_A', 0):.2f} A")
                self._output.appendPlainText(
                    f"  Coil sigma:   "
                    f"{result.get('coil_sigma_Sm', 0):.4e} S/m")
                self._output.appendPlainText(
                    "  --- Workpiece ---")
                self._output.appendPlainText(
                    f"  Material:  {result.get('material', '?')}")
                self._output.appendPlainText(
                    f"  WP sigma:  "
                    f"{result.get('wp_sigma_Sm', 0):.4e} S/m")
                self._output.appendPlainText(
                    f"  WP mu_r:   "
                    f"{result.get('wp_mu_r', 0):.1f}")
                self._output.appendPlainText(
                    "  --- SIBC ---")
                self._output.appendPlainText(
                    f"  Frequency:  "
                    f"{result.get('frequency_Hz', 0):.0f} Hz")
                self._output.appendPlainText(
                    f"  Skin depth: "
                    f"{result.get('skin_depth_mm', 0):.3f} mm")
                self._output.appendPlainText(
                    f"  |Z_s| = {result.get('Z_s_abs_Ohm', 0):.4e} Ohm"
                    f"  (phase = {result.get('Z_s_phase_deg', 0):.1f} deg)")
                self._output.appendPlainText(
                    f"  Re(Z_s) = {result.get('Z_s_real_Ohm', 0):.4e}, "
                    f"Im(Z_s) = {result.get('Z_s_imag_Ohm', 0):.4e}")
                self._output.appendPlainText(
                    "  --- Results ---")
                self._output.appendPlainText(
                    f"  P_total = {result['P_total_W']:.4e} W")
                if "Q_total_var" in result:
                    self._output.appendPlainText(
                        f"  Q_total = {result['Q_total_var']:.4e} var")
                self._output.appendPlainText(
                    f"  H_t_rms = {result.get('H_t_rms_Am', 0):.4f} A/m")
                self._output.appendPlainText(
                    "  --- Solver ---")
                self._output.appendPlainText(
                    f"  Karl iterations = {result.get('n_iter', '?')}")
                self._output.appendPlainText(
                    f"  BEM DOFs = {result.get('ndof', '?')}, "
                    f"elements = {result.get('n_elements', '?')}")
            # FEM-SIBC results (calc_fem_kelvin)
            if "P_total" in result and "P_total_W" not in result:
                self._output.appendPlainText(
                    "  --- FEM-SIBC ---")
                self._output.appendPlainText(
                    f"  P_total = {result['P_total']:.4e} W")
                if "Q_total" in result:
                    self._output.appendPlainText(
                        f"  Q_total = {result['Q_total']:.4e} var")
                if "H_t_rms" in result:
                    self._output.appendPlainText(
                        f"  H_t_rms = {result['H_t_rms']:.4f} A/m")
                if "L" in result:
                    self._output.appendPlainText(
                        f"  L = {result['L']*1e9:.2f} nH")
                if "Z_s" in result:
                    self._output.appendPlainText(
                        f"  Z_s = {result['Z_s']}")
                if "delta" in result:
                    self._output.appendPlainText(
                        f"  Skin depth = {result['delta']*1e3:.3f} mm")
                if "ndof" in result:
                    self._output.appendPlainText(
                        f"  DOFs = {result['ndof']}, "
                        f"Elements = {result.get('ne', '?')}")
                if "iterations" in result:
                    self._output.appendPlainText(
                        f"  Karl iterations = {result['iterations']}")
                if "t_total" in result:
                    self._output.appendPlainText(
                        f"  Time = {result['t_total']:.1f}s")
            # BEM coil + workpiece SIBC/ESIM (calc_inductance.py with
            # --workpiece). Distinct from BEM-SIBC (WP) which uses
            # P_total_W; this path uses wp_P_total + wp_R_effective.
            if "wp_P_total" in result:
                self._output.appendPlainText(
                    "  --- Workpiece SIBC ---")
                self._output.appendPlainText(
                    f"  Model:    {result.get('wp_model', '?')}")
                self._output.appendPlainText(
                    f"  Material: {result.get('wp_material', '?')}")
                self._output.appendPlainText(
                    f"  Frequency: {result.get('wp_frequency', 0):.0f} Hz")
                self._output.appendPlainText(
                    f"  WP sigma:  {result.get('wp_sigma', 0):.4e} S/m")
                if "wp_delta_min" in result:
                    self._output.appendPlainText(
                        f"  Skin depth: "
                        f"{result['wp_delta_min']*1e3:.4f} mm")
                self._output.appendPlainText(
                    f"  Panels: {result.get('wp_n_panels', '?')}")
                if "wp_H_t_max" in result:
                    self._output.appendPlainText(
                        f"  H_t range: {result['wp_H_t_min']:.2f} - "
                        f"{result['wp_H_t_max']:.2f} A/m")
                self._output.appendPlainText(
                    f"  P_total = {result['wp_P_total']:.4e} W")
                if "wp_Q_total" in result:
                    self._output.appendPlainText(
                        f"  Q_total = {result['wp_Q_total']:.4e} var")
                if "wp_R_effective" in result:
                    self._output.appendPlainText(
                        f"  R_eff   = "
                        f"{result['wp_R_effective']:.4e} Ohm")
                if "wp_X_effective" in result:
                    self._output.appendPlainText(
                        f"  X_eff   = "
                        f"{result['wp_X_effective']:.4e} Ohm")
                # Coupled coil-terminal results.
                # Two cases:
                #   1. Linear SIBC (Dowell): full coupled BEM solve via
                #      bem_coupled_solver.CoupledBEMSolver. Reports
                #      L_air, L_total, Delta_L (sign is physically
                #      correct: <0 for non-magnetic Lenz screening,
                #      >0 for ferromagnetic flux concentration).
                #   2. Nonlinear SIBC (ESIM): uncoupled estimator,
                #      reports R only (no Delta L).
                if "coupled_L_total_H" in result:
                    L_air_nH = result.get(
                        "L_air_H", result.get("coupled_L_air_H", 0)) * 1e9
                    L_total_nH = result["coupled_L_total_H"] * 1e9
                    dL_nH = result["coupled_dL_H"] * 1e9
                    self._output.appendPlainText(
                        "  --- Coil terminal (with workpiece) ---")
                    self._output.appendPlainText(
                        f"  L (air)   = {L_air_nH:.3f} nH")
                    self._output.appendPlainText(
                        f"  delta L  = {dL_nH:+.3f} nH")
                    self._output.appendPlainText(
                        f"  L (eff)   = {L_total_nH:.3f} nH")
                    if "coupled_R_effective_Ohm" in result:
                        self._output.appendPlainText(
                            f"  R (added) = "
                            f"{result['coupled_R_effective_Ohm']*1e3:.4f}"
                            f" mOhm")
                    if "coupled_iterations" in result:
                        self._output.appendPlainText(
                            f"  iters     = "
                            f"{result['coupled_iterations']} "
                            f"(Picard, relax=0.5)")
                    if "coupled_delta_skin_m" in result:
                        self._output.appendPlainText(
                            f"  skin depth = "
                            f"{result['coupled_delta_skin_m']*1e3:.4f} mm")
                    # Per-panel curvature SIBC (Phase 5, 2026-04-12):
                    # Surface the local-R range and panel count so the
                    # user can sanity-check that the mesh-driven extractor
                    # picked up the workpiece geometry correctly.
                    if result.get("coupled_use_local_curvature"):
                        rmin = result.get("coupled_R_local_min_m", 0) * 1e3
                        rmax = result.get("coupled_R_local_max_m", 0) * 1e3
                        npan = result.get("coupled_n_wp_panels", 0)
                        self._output.appendPlainText(
                            f"  R_local   = [{rmin:.3f}, {rmax:.3f}] mm "
                            f"(per-panel, {npan} panels)")
                    elif "coupled_use_local_curvature" in result:
                        self._output.appendPlainText(
                            "  R_local   = global half-thickness "
                            "(per-panel curvature: off)")
                    self._output.appendPlainText(
                        "  (Coupled BEM: per-DOF f_back, validated"
                        " 2026-04-12 vs FEM-Kelvin SIBC at 0.3% on"
                        " copper / 1.7% on steel mu_r=100. See"
                        " examples/cubit_panels/inductance/"
                        "compare_bem_coupled_vs_fem_kelvin.py)")
                elif "coupled_R_effective_Ohm" in result:
                    self._output.appendPlainText(
                        "  --- Coil terminal (with workpiece, ESIM) ---")
                    self._output.appendPlainText(
                        f"  L (air-only) = "
                        f"{result.get('inductance_H', 0)*1e9:.3f} nH")
                    self._output.appendPlainText(
                        f"  R (added)   = "
                        f"{result['coupled_R_effective_Ohm']*1e3:.4f}"
                        f" mOhm")
                    self._output.appendPlainText(
                        "  (delta L not reported: ESIM is one-way only;"
                        " coupled BEM only supports linear Dowell SIBC.)")
                if "coupled_error" in result:
                    self._output.appendPlainText(
                        "  --- Coupled BEM ERROR ---")
                    self._output.appendPlainText(
                        f"  {result['coupled_error']}")
            elif "wp_error" in result:
                self._output.appendPlainText(
                    "  --- Workpiece SIBC ---")
                self._output.appendPlainText(
                    f"  ERROR: {result['wp_error']}")

            # Volume/area results (calc_volume, calc_surface)
            if "ng_volume" in result:
                self._output.appendPlainText(
                    f"  Volume = {result['ng_volume']:.6e}")
            if "vol_error_pct" in result:
                self._output.appendPlainText(
                    f"  Volume error = {result['vol_error_pct']:+.4e}%")
            if "ng_area" in result:
                self._output.appendPlainText(
                    f"  Area = {result['ng_area']:.6e}")
            self._output.appendPlainText("")

        self._status.showMessage("Done." if exit_code == 0 else "Failed.")

    def _open_gmsh(self):
        if not self._last_msh:
            return
        import subprocess
        # Build a small Python launcher that:
        #   - merges the .msh / .geo we just produced
        #   - sets sensible defaults so the curved surfaces render
        #     correctly and air-mesh edges do not cover the screen
        #   - shows the first vector view (J) as arrows by default
        launcher = (
            "import gmsh; gmsh.initialize();"
            f" gmsh.merge(r'{self._last_msh}');"
            # 1. Curved-element subdivision (Tri6/10/15 etc.)
            " gmsh.option.setNumber('Mesh.NumSubEdges', 4);"
            # 2. Hide volume mesh edges (kills the black-line storm)
            " gmsh.option.setNumber('Mesh.VolumeEdges', 0);"
            " gmsh.option.setNumber('Mesh.VolumeFaces', 0);"
            # 3. Surface edges off, surface faces on so the coil
            #    looks like a smooth shaded body, not a wireframe
            " gmsh.option.setNumber('Mesh.SurfaceEdges', 0);"
            " gmsh.option.setNumber('Mesh.SurfaceFaces', 1);"
            # 4. If a vector view is present, draw it as small arrows
            #    on the surface (the J view from calc_inductance.py)
            " tags = gmsh.view.getTags();"
            " [gmsh.option.setNumber(f'View[{i}].VectorType', 4)"
            "  for i, _ in enumerate(tags)];"
            " [gmsh.option.setNumber(f'View[{i}].ArrowSizeMax', 60)"
            "  for i, _ in enumerate(tags)];"
            " gmsh.fltk.run(); gmsh.finalize()"
        )
        subprocess.Popen(
            [_PYTHON, "-c", launcher],
            creationflags=(0x08000000 if sys.platform == "win32" else 0))

    # Settings
    def _settings_path(self):
        return os.path.join(_SETTINGS_DIR, f"radia_{self._settings_key}.json")

    def _save_settings(self):
        os.makedirs(_SETTINGS_DIR, exist_ok=True)
        data = {
            "vol": self._vol_edit.text(),
            "panel": self._panel.save_state() if self._panel else {},
        }
        try:
            with open(self._settings_path(), "w", encoding="utf-8") as f:
                json.dump(data, f, indent=2)
        except OSError:
            pass

    def _restore_settings(self):
        p = self._settings_path()
        if not os.path.isfile(p):
            return
        try:
            with open(p, "r", encoding="utf-8") as f:
                data = json.load(f)
        except (OSError, json.JSONDecodeError):
            return
        if not self._vol_edit.text() and "vol" in data:
            self._vol_edit.setText(data["vol"])
        if self._panel:
            self._panel.restore_state(data.get("panel", {}))

    def closeEvent(self, event):
        try:
            self._save_settings()
        except Exception:
            pass
        if self._process:
            try:
                self._process.kill()
                self._process.waitForFinished(3000)
            except Exception:
                pass
        super().closeEvent(event)


def run_app(window_class, vol_path=""):
    """Entry point helper: create QApplication + window, exec."""
    app = QApplication(sys.argv)
    vol = vol_path or (sys.argv[1] if len(sys.argv) > 1 else "")
    window = window_class(vol)
    window.show()
    sys.exit(app.exec())
