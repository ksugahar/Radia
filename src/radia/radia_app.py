"""
Radia standalone analysis app (PySide6).

Launched from Cubit (1-button panel) or command line:
    python -m radia.radia_app model.cub5
    radia-app model.cub5
    radia-app                  # no model, browse later

Modes:
  - Mesh Evaluation    (calc_mesh_eval.py)
  - IH (BEM)           (calc_inductance.py)
  - IH (FEM)           (calc_fem_kelvin.py)
  - Electromagnet (FEM)(calc_accel_magnet.py)
  - Electromagnet (MSC)(calc_accel_msc.py)
  - PCB (PEEC)         (calc_pcb_peec.py)
"""

import json
import os
import sys

from PySide6.QtCore import Qt, QProcess
from PySide6.QtGui import QFont, QIcon
from PySide6.QtWidgets import (
    QApplication, QMainWindow, QWidget,
    QVBoxLayout, QHBoxLayout, QFormLayout,
    QTabWidget, QLabel, QLineEdit, QComboBox, QSpinBox,
    QPushButton, QPlainTextEdit, QFileDialog,
    QStatusBar, QMessageBox, QSplitter, QGroupBox, QStyle,
)

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_PANELS_DIR = os.path.join(_THIS_DIR, "panels")
_RESOURCES_DIR = os.path.join(_THIS_DIR, "resources")
_PYTHON = sys.executable
_SETTINGS_PATH = os.path.join(os.path.expanduser("~"), ".radia", "radia_app.json")


def _icon_path():
    """Return path to radia_icon.png (or .ico on Windows)."""
    for ext in (".ico", ".png"):
        p = os.path.join(_RESOURCES_DIR, "radia_icon" + ext)
        if os.path.isfile(p):
            return p
    return ""


def _calc_script(name):
    return os.path.join(_PANELS_DIR, name)


# ============================================================
# Settings persistence
# ============================================================

def _load_settings():
    if os.path.isfile(_SETTINGS_PATH):
        try:
            with open(_SETTINGS_PATH, "r", encoding="utf-8") as f:
                return json.load(f)
        except Exception:
            pass
    return {}


def _save_settings(data):
    try:
        settings_dir = os.path.dirname(_SETTINGS_PATH)
        if settings_dir and not os.path.isdir(settings_dir):
            os.makedirs(settings_dir, exist_ok=True)
        with open(_SETTINGS_PATH, "w", encoding="utf-8") as f:
            json.dump(data, f, indent=2, ensure_ascii=False)
    except Exception:
        pass


# ============================================================
# Mode Tabs
# ============================================================

class ModeTab(QWidget):
    """Base class for analysis mode tabs."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self._form = QFormLayout()
        self._form.setLabelAlignment(Qt.AlignRight)
        self.setLayout(self._form)
        self._widgets = {}  # key -> widget, for save/restore
        self._build_ui()

    def _build_ui(self):
        """Override to add widgets via _add_*() helpers."""
        raise NotImplementedError

    def build_command(self, cub5):
        """Return [python, script, ...args] list. Override per tab."""
        raise NotImplementedError

    # -- helpers --

    def _add_line(self, key, label, default="", placeholder=""):
        w = QLineEdit(default)
        if placeholder:
            w.setPlaceholderText(placeholder)
        self._form.addRow(label, w)
        self._widgets[key] = w
        return w

    def _add_combo(self, key, label, items, default=0):
        w = QComboBox()
        w.addItems(items)
        w.setCurrentIndex(default)
        self._form.addRow(label, w)
        self._widgets[key] = w
        return w

    def _add_spin(self, key, label, value=1, lo=1, hi=999):
        w = QSpinBox()
        w.setRange(lo, hi)
        w.setValue(value)
        self._form.addRow(label, w)
        self._widgets[key] = w
        return w

    def _add_browse(self, key, label, default="", filter_str="All files (*.*)"):
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
        return le

    def _do_browse(self, line_edit, filter_str):
        path, _ = QFileDialog.getOpenFileName(self, "Select file", "", filter_str)
        if path:
            line_edit.setText(path)

    def _val(self, key):
        """Get string value of a widget by key."""
        w = self._widgets[key]
        if isinstance(w, QLineEdit):
            return w.text().strip()
        elif isinstance(w, QComboBox):
            return w.currentText()
        elif isinstance(w, QSpinBox):
            return str(w.value())
        return ""

    def _msh_output(self, cub5, suffix):
        """Generate .msh output path next to the input file."""
        base = cub5 if cub5 else "output"
        return os.path.splitext(base)[0] + suffix + ".msh"

    def save_state(self):
        state = {}
        for key, w in self._widgets.items():
            if isinstance(w, QLineEdit):
                state[key] = w.text()
            elif isinstance(w, QComboBox):
                state[key] = w.currentIndex()
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
                    w.setCurrentIndex(int(val))
                elif isinstance(w, QSpinBox):
                    w.setValue(int(val))
            except (ValueError, TypeError):
                pass


class IHBEMTab(ModeTab):
    def _build_ui(self):
        self._add_spin("order", "Curve order:", 2, 1, 5)
        self._add_spin("fes_order", "FES order:", 0, 0, 5)
        self._add_line("source", "Source block:", "source")
        self._add_line("sink", "Sink block:", "sink")
        self._add_line("workpiece", "Workpiece block:", "",
                       placeholder="(empty = inductance only)")
        self._add_combo("impedance", "Impedance model:",
                        ["dowell", "esim", "bem-sibc"])
        self._add_line("freq", "Frequency [Hz]:", "50000")
        self._add_line("sigma", "Sigma [S/m]:", "2e6")

        # Material: mu_r / BH (no hysteresis for IH)
        mat_combo = self._add_combo("material", "Material:",
                                    ["mu_r (Linear)", "BH Curve"])
        self._add_line("mu_r", "mu_r:", "100")
        self._add_browse("bh_file", "BH file:",
                         filter_str="Text files (*.txt *.csv);;All (*)")
        mat_combo.currentIndexChanged.connect(self._on_material_changed)

        self._add_line("half_thickness", "Half thickness [m]:", "0.005")
        self._add_combo("esim_geometry", "ESIM geometry:",
                        ["local_curvature", "planar"])
        self._on_material_changed(0)

    def _on_material_changed(self, idx):
        mu_r_w = self._widgets["mu_r"]
        bh_w = self._widgets["bh_file"]
        mu_r_w.parentWidget().setVisible(idx == 0) if mu_r_w.parentWidget() else None
        bh_w.parentWidget().setVisible(idx == 1) if bh_w.parentWidget() else None

    def build_command(self, model_path):
        if not model_path:
            raise ValueError("No model file specified.")
        if model_path.endswith(".vol"):
            cmd = [_PYTHON, _calc_script("calc_inductance.py"),
                   "--vol", model_path,
                   "--source", self._val("source"),
                   "--sink", self._val("sink")]
        else:
            cmd = [_PYTHON, _calc_script("calc_inductance.py"),
                   "--cub5", model_path,
                   "--order", self._val("order"),
                   "--source", self._val("source"),
                   "--sink", self._val("sink")]
        fes = self._val("fes_order")
        if fes and fes != "0":
            cmd += ["--fes-order", fes]
        wp = self._val("workpiece")
        if wp:
            mat_idx = self._widgets["material"].currentIndex()
            mat_name = "custom" if mat_idx == 0 else "steel"
            cmd += ["--workpiece", wp,
                    "--impedance-model", self._val("impedance"),
                    "--frequency", self._val("freq"),
                    "--sigma", self._val("sigma"),
                    "--material", mat_name,
                    "--half-thickness", self._val("half_thickness"),
                    "--esim-geometry", self._val("esim_geometry")]
            if mat_idx == 0:
                cmd += ["--mu-r", self._val("mu_r")]
            else:
                bh = self._val("bh_file")
                if bh:
                    cmd += ["--bh-file", bh]
        cmd += ["--msh-output", self._msh_output(model_path, "_bem")]
        return cmd


class IHFEMTab(ModeTab):
    def _build_ui(self):
        self._add_spin("order", "Curve order:", 2, 1, 5)
        self._add_spin("fes_order", "FES order:", 1, 1, 5)
        self._add_line("freq", "Frequency [Hz]:", "7000")
        self._add_line("sigma", "Sigma [S/m]:", "2e6")
        self._add_combo("impedance", "Impedance:", ["sibc", "esim"])

        # Material: mu_r / BH (no hysteresis for IH)
        mat_combo = self._add_combo("material", "Material:",
                                    ["mu_r (Linear)", "BH Curve"])
        self._add_line("mu_r", "mu_r:", "100")
        self._add_browse("bh_file", "BH file:",
                         filter_str="Text files (*.txt *.csv);;All (*)")
        mat_combo.currentIndexChanged.connect(self._on_material_changed)

        self._add_line("current", "Current [A]:", "1.0")
        self._add_line("a_coil", "Coil radius [m]:", "0.003")
        self._add_line("r_wp", "Workpiece radius [m]:", "0.010")
        self._add_combo("solver", "Solver:",
                        ["pardiso", "bddc", "iccg", "ams"])
        self._add_spin("max_iter", "Max iterations:", 15, 1, 200)
        self._on_material_changed(0)

    def _on_material_changed(self, idx):
        mu_r_w = self._widgets["mu_r"]
        bh_w = self._widgets["bh_file"]
        mu_r_w.parentWidget().setVisible(idx == 0) if mu_r_w.parentWidget() else None
        bh_w.parentWidget().setVisible(idx == 1) if bh_w.parentWidget() else None

    def build_command(self, cub5):
        if not cub5:
            raise ValueError("No .cub5 file specified.")
        mat_idx = self._widgets["material"].currentIndex()
        if mat_idx == 0:
            mat_args = ["--material", "custom", "--mu-r", self._val("mu_r")]
        else:
            mat_args = ["--material", "steel"]
            bh = self._val("bh_file")
            if bh:
                mat_args += ["--bh-file", bh]
        cmd = [_PYTHON, _calc_script("calc_fem_kelvin.py"),
               "--cub5", cub5,
               "--order", self._val("order"),
               "--fes-order", self._val("fes_order"),
               "--frequency", self._val("freq"),
               "--sigma", self._val("sigma"),
               "--impedance", self._val("impedance"),
               "--current", self._val("current"),
               "--a-coil", self._val("a_coil"),
               "--r-wp", self._val("r_wp"),
               "--solver", self._val("solver"),
               "--max-iter", self._val("max_iter"),
               "--msh-output", self._msh_output(cub5, "_fem")]
        cmd += mat_args
        return cmd


class EMFEMTab(ModeTab):
    _OMEGA_SOLVERS = ["Direct (PARDISO)", "AMG (Compact)", "BDDC"]
    _A_SOLVERS = ["Direct (PARDISO)", "AMS (Compact)", "BDDC"]
    _SOLVER_MAP = {
        "Direct (PARDISO)": "pardiso",
        "AMG (Compact)": "amg",
        "AMS (Compact)": "ams",
        "BDDC": "bddc",
    }

    def _build_ui(self):
        form_combo = self._add_combo("formulation", "Formulation:",
                                     ["omega", "a-phi"])
        self._add_spin("order", "Curve order:", 2, 1, 5)
        self._add_spin("fes_order", "FES order:", 1, 1, 5)
        self._add_browse("coil_script", "Coil script:",
                         filter_str="Python (*.py);;All (*)")

        # Material: mu_r / BH / Hysteresis
        mat_combo = self._add_combo("material", "Material:",
                                    ["mu_r (Linear)", "BH Curve",
                                     "Hysteresis (.hys)"])
        self._add_line("mu_r", "mu_r:", "1000")
        self._add_browse("bh_file", "BH file:",
                         filter_str="Text files (*.txt *.csv);;All (*)")
        self._add_browse("hys_file", "Hysteresis file:",
                         filter_str="HYS files (*.hys);;All (*)")
        mat_combo.currentIndexChanged.connect(self._on_material_changed)

        # Solver: depends on formulation
        self._add_combo("solver", "Solver:", self._OMEGA_SOLVERS)
        form_combo.currentTextChanged.connect(self._on_formulation_changed)

        self._add_line("ima", "IMA symmetry:", "",
                       placeholder="e.g. +x-z")
        self._add_line("tol", "Tolerance:", "1e-3")
        self._add_spin("max_iter", "Max iterations:", 100, 1, 9999)
        self._add_spin("n_steps", "N steps:", 1, 1, 100)
        self._add_line("relax", "Relaxation:", "0.0")
        self._add_combo("method", "Method:", ["Picard", "Newton"])

        # Initial visibility
        self._on_material_changed(0)

    def _on_formulation_changed(self, text):
        combo = self._widgets["solver"]
        prev = combo.currentText()
        combo.clear()
        items = self._A_SOLVERS if text == "a-phi" else self._OMEGA_SOLVERS
        combo.addItems(items)
        # Try to keep "Direct" or "BDDC" selected if switching
        for i, item in enumerate(items):
            if prev.split()[0] == item.split()[0]:
                combo.setCurrentIndex(i)
                break

    def _on_material_changed(self, idx):
        # 0=mu_r, 1=BH, 2=Hysteresis
        mu_r_w = self._widgets["mu_r"]
        bh_w = self._widgets["bh_file"]
        hys_w = self._widgets["hys_file"]
        mu_r_w.parentWidget().setVisible(idx == 0) if mu_r_w.parentWidget() else None
        bh_w.parentWidget().setVisible(idx == 1) if bh_w.parentWidget() else None
        hys_w.parentWidget().setVisible(idx == 2) if hys_w.parentWidget() else None

    def build_command(self, cub5):
        coil = self._val("coil_script")
        if not coil:
            raise ValueError("No coil script specified.")
        solver_key = self._SOLVER_MAP.get(self._val("solver"), "pardiso")

        mat_idx = self._widgets["material"].currentIndex()
        if mat_idx == 0:
            mat_args = ["--material", "linear", "--mu-r", self._val("mu_r")]
        elif mat_idx == 1:
            mat_args = ["--material", "steel"]
            bh = self._val("bh_file")
            if bh:
                mat_args += ["--bh-file", bh]
        else:
            mat_args = ["--material", "hysteresis"]
            hys = self._val("hys_file")
            if hys:
                mat_args += ["--hys-file", hys]

        cmd = [_PYTHON, _calc_script("calc_accel_magnet.py"),
               "--formulation", self._val("formulation"),
               "--order", self._val("order"),
               "--fes-order", self._val("fes_order"),
               "--coil-script", coil,
               "--solver", solver_key,
               "--tol", self._val("tol"),
               "--max-iter", self._val("max_iter"),
               "--n-steps", self._val("n_steps"),
               "--relax", self._val("relax")]
        cmd += mat_args
        if cub5:
            cmd += ["--cub5", cub5]
        ima = self._val("ima")
        if ima:
            cmd += ["--ima", ima]
        if self._val("method") == "Newton":
            cmd += ["--newton"]
        cmd += ["--msh-output", self._msh_output(
            cub5 or coil, "_emfem")]
        return cmd


class EMMSCTab(ModeTab):
    def _build_ui(self):
        self._add_browse("coil_script", "Coil script:",
                         filter_str="Python (*.py);;All (*)")

        # Material: mu_r / BH / Hysteresis
        mat_combo = self._add_combo("material", "Material:",
                                    ["mu_r (Linear)", "BH Curve",
                                     "Hysteresis (.hys)"])
        self._add_line("mu_r", "mu_r:", "1000")
        self._add_browse("bh_file", "BH file:",
                         filter_str="Text files (*.txt *.csv);;All (*)")
        self._add_browse("hys_file", "Hysteresis file:",
                         filter_str="HYS files (*.hys);;All (*)")
        mat_combo.currentIndexChanged.connect(self._on_material_changed)

        self._add_combo("solver", "Solver:",
                        ["0 (LU)", "1 (BiCGSTAB)", "2 (HACApK)"])
        self._add_line("ima", "IMA symmetry:", "",
                       placeholder="e.g. +x-z")
        self._add_line("tol", "Tolerance:", "1e-3")
        self._add_spin("max_iter", "Max iterations:", 100, 1, 9999)
        self._add_line("relax", "Relaxation:", "0.0")
        self._add_line("unit_scale", "Unit scale:", "0.001")
        self._on_material_changed(0)

    def _on_material_changed(self, idx):
        mu_r_w = self._widgets["mu_r"]
        bh_w = self._widgets["bh_file"]
        hys_w = self._widgets.get("hys_file")
        mu_r_w.parentWidget().setVisible(idx == 0) if mu_r_w.parentWidget() else None
        bh_w.parentWidget().setVisible(idx == 1) if bh_w.parentWidget() else None
        if hys_w:
            hys_w.parentWidget().setVisible(idx == 2) if hys_w.parentWidget() else None

    def build_command(self, cub5):
        coil = self._val("coil_script")
        if not coil:
            raise ValueError("No coil script specified.")
        solver_id = self._val("solver").split()[0]

        mat_idx = self._widgets["material"].currentIndex()
        if mat_idx == 0:
            mat_args = ["--material", "custom", "--mu-r", self._val("mu_r")]
        elif mat_idx == 1:
            mat_args = ["--material", "steel"]
            bh = self._val("bh_file")
            if bh:
                mat_args += ["--bh-file", bh]
        else:
            mat_args = ["--material", "hysteresis"]
            hys = self._val("hys_file")
            if hys:
                mat_args += ["--hys-file", hys]

        cmd = [_PYTHON, _calc_script("calc_accel_msc.py"),
               "--coil-script", coil,
               "--solver", solver_id,
               "--tol", self._val("tol"),
               "--max-iter", self._val("max_iter"),
               "--relax", self._val("relax"),
               "--unit-scale", self._val("unit_scale")]
        cmd += mat_args
        if cub5:
            cmd += ["--cub5", cub5]
        ima = self._val("ima")
        if ima:
            cmd += ["--ima", ima]
        cmd += ["--msh-output", self._msh_output(
            cub5 or coil, "_msc")]
        return cmd


class PCBPEECTab(ModeTab):
    def _build_ui(self):
        self._add_browse("inp_file", "FastHenry .inp:",
                         filter_str="FastHenry (*.inp);;All (*)")
        self._add_line("freq_min", "Freq min [Hz]:", "1e3")
        self._add_line("freq_max", "Freq max [Hz]:", "1e9")
        self._add_spin("n_freq", "N freq points:", 50, 1, 1000)
        self._add_combo("solver", "Solver:",
                        ["0 (LU)", "1 (BiCGSTAB)", "2 (HACApK)"])
        self._add_line("spice_output", "SPICE output:", "",
                       placeholder="e.g. output.cir")

    def build_command(self, cub5):
        inp = self._val("inp_file")
        if not inp:
            raise ValueError("No FastHenry .inp file specified.")
        solver_id = self._val("solver").split()[0]
        cmd = [_PYTHON, _calc_script("calc_pcb_peec.py"),
               "--inp", inp,
               "--freq-min", self._val("freq_min"),
               "--freq-max", self._val("freq_max"),
               "--n-freq", self._val("n_freq"),
               "--solver-method", solver_id]
        spice = self._val("spice_output")
        if spice:
            cmd += ["--spice-output", spice]
        return cmd


# ============================================================
# Main Window
# ============================================================

class RadiaApp(QMainWindow):

    TAB_CLASSES = [
        ("IH (BEM)",   IHBEMTab),
        ("IH (FEM)",   IHFEMTab),
        ("EM (FEM)",   EMFEMTab),
        ("EM (MSC)",   EMMSCTab),
        ("PCB (PEEC)", PCBPEECTab),
    ]

    def __init__(self, cub5_path=""):
        super().__init__()
        self.setWindowTitle("Radia Analysis Tool")
        self.resize(720, 650)
        self.setMinimumSize(550, 450)

        icon = _icon_path()
        if icon:
            self.setWindowIcon(QIcon(icon))

        self._process = None
        self._last_msh = None

        self._build_ui(cub5_path)
        self._restore_settings()

    def _build_ui(self, cub5_path):
        central = QWidget()
        self.setCentralWidget(central)
        root = QVBoxLayout(central)
        root.setContentsMargins(10, 10, 10, 5)

        # -- Model path --
        model_row = QHBoxLayout()
        model_row.addWidget(QLabel("Model:"))
        self._cub5 = QLineEdit(cub5_path)
        self._cub5.setPlaceholderText(".cub5 (Cubit) or .step (OCC) file")
        model_row.addWidget(self._cub5, 1)
        browse_btn = QPushButton("Browse...")
        browse_btn.clicked.connect(self._browse_model)
        model_row.addWidget(browse_btn)
        root.addLayout(model_row)

        # -- Splitter: tabs | output --
        splitter = QSplitter(Qt.Vertical)
        root.addWidget(splitter, 1)

        # Tabs
        self._tabs = QTabWidget()
        self._tab_instances = []
        for name, cls in self.TAB_CLASSES:
            tab = cls()
            self._tabs.addTab(tab, name)
            self._tab_instances.append(tab)
        splitter.addWidget(self._tabs)

        # Output
        out_group = QGroupBox("Output")
        out_layout = QVBoxLayout(out_group)
        out_layout.setContentsMargins(5, 5, 5, 5)
        self._output = QPlainTextEdit()
        self._output.setReadOnly(True)
        self._output.setFont(QFont("Consolas", 9))
        out_layout.addWidget(self._output)
        splitter.addWidget(out_group)

        splitter.setSizes([350, 250])

        # -- Buttons with icons --
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

        self._save_btn = QPushButton(
            style.standardIcon(QStyle.SP_DialogSaveButton), " Save .py")
        self._save_btn.setFixedHeight(32)
        self._save_btn.clicked.connect(self._on_save_script)
        btn_row.addWidget(self._save_btn)

        btn_row.addStretch()

        self._gmsh_btn = QPushButton(
            style.standardIcon(QStyle.SP_ComputerIcon), " Open GMSH")
        self._gmsh_btn.setFixedHeight(32)
        self._gmsh_btn.setEnabled(False)
        self._gmsh_btn.clicked.connect(self._open_gmsh)
        btn_row.addWidget(self._gmsh_btn)
        root.addLayout(btn_row)

        # Status bar
        self._status = QStatusBar()
        self.setStatusBar(self._status)
        self._status.showMessage("Ready")

    # --------------------------------------------------------
    # Actions
    # --------------------------------------------------------

    def _browse_model(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Select model",
            os.path.dirname(self._cub5.text()),
            "Cubit (*.cub5);;STEP (*.step *.stp);;All (*)")
        if path:
            self._cub5.setText(path)

    def _on_save_script(self):
        tab = self._tab_instances[self._tabs.currentIndex()]
        cub5 = self._cub5.text().strip()
        try:
            cmd = tab.build_command(cub5)
        except ValueError as e:
            QMessageBox.warning(self, "Input Error", str(e))
            return

        default_name = os.path.splitext(os.path.basename(cub5))[0] if cub5 else "solve"
        default_path = os.path.join(
            os.path.dirname(cub5) if cub5 else os.getcwd(),
            f"run_{default_name}.py")
        path, _ = QFileDialog.getSaveFileName(
            self, "Save script", default_path,
            "Python (*.py);;All (*)")
        if not path:
            return

        cmd_str = ", ".join(repr(c) for c in cmd)
        script = (
            '#!/usr/bin/env python\n'
            '"""Generated by Radia App"""\n'
            'import subprocess, sys\n'
            f'cmd = [{cmd_str}]\n'
            'print("Running:", " ".join(cmd))\n'
            'sys.exit(subprocess.run(cmd).returncode)\n'
        )
        with open(path, "w", encoding="utf-8") as f:
            f.write(script)
        self._status.showMessage(f"Saved: {path}")

    def _on_run(self):
        if self._process is not None:
            return
        tab = self._tab_instances[self._tabs.currentIndex()]
        cub5 = self._cub5.text().strip()
        try:
            cmd = tab.build_command(cub5)
        except ValueError as e:
            QMessageBox.warning(self, "Input Error", str(e))
            return

        self._save_settings()
        self._output.clear()
        self._output.appendPlainText(
            f"[{self._tabs.tabText(self._tabs.currentIndex())}]")
        self._output.appendPlainText(f"> {' '.join(cmd)}\n")
        self._last_msh = None
        self._gmsh_btn.setEnabled(False)

        # Determine working directory
        work_dir = os.path.dirname(cub5) if cub5 else os.getcwd()
        if not os.path.isdir(work_dir):
            work_dir = os.getcwd()

        # Launch QProcess
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

    def _read_stdout(self):
        data = self._process.readAllStandardOutput().data()
        text = data.decode("utf-8", errors="replace")
        self._output.appendPlainText(text.rstrip("\n"))

    def _read_stderr(self):
        data = self._process.readAllStandardError().data()
        text = data.decode("utf-8", errors="replace")
        # stderr often carries progress info
        for line in text.strip().split("\n"):
            line = line.strip()
            if line:
                self._status.showMessage(line)

    def _on_finished(self, exit_code, exit_status):
        # Collect any remaining stdout
        remaining = self._process.readAllStandardOutput().data()
        if remaining:
            self._output.appendPlainText(
                remaining.decode("utf-8", errors="replace").rstrip("\n"))

        stdout_text = self._output.toPlainText()
        self._process = None
        self._run_btn.setEnabled(True)
        self._stop_btn.setEnabled(False)

        if exit_code != 0:
            self._status.showMessage(f"Error (exit code {exit_code})")
            self._output.appendPlainText(f"\n*** Process exited with code {exit_code}")
            return

        # Parse JSON result from stdout
        result = None
        for line in reversed(stdout_text.split("\n")):
            line = line.strip()
            if line.startswith("{"):
                try:
                    result = json.loads(line)
                    break
                except json.JSONDecodeError:
                    continue

        if result:
            self._output.appendPlainText("\n--- Results ---")
            for key, val in result.items():
                if isinstance(val, list) and len(val) > 10:
                    self._output.appendPlainText(f"  {key}: [{len(val)} items]")
                else:
                    self._output.appendPlainText(f"  {key}: {val}")

            msh = result.get("msh_file", "")
            if msh and os.path.isfile(msh):
                self._last_msh = msh
                self._gmsh_btn.setEnabled(True)

        self._status.showMessage("Done")

    # --------------------------------------------------------
    # GMSH
    # --------------------------------------------------------

    def _open_gmsh(self):
        if not self._last_msh or not os.path.isfile(self._last_msh):
            return
        pyw = _PYTHON.replace("python.exe", "pythonw.exe")
        if not os.path.exists(pyw):
            pyw = _PYTHON
        code = (
            "import gmsh; gmsh.initialize(); "
            f"gmsh.merge(r'{self._last_msh}'); "
            "gmsh.fltk.run(); gmsh.finalize()")
        try:
            import subprocess
            subprocess.Popen([pyw, "-c", code])
        except Exception:
            try:
                os.startfile(self._last_msh)
            except Exception as e:
                QMessageBox.warning(self, "Error", str(e))

    # --------------------------------------------------------
    # Settings
    # --------------------------------------------------------

    def _save_settings(self):
        data = {"cub5": self._cub5.text(), "active_tab": self._tabs.currentIndex()}
        for i, tab in enumerate(self._tab_instances):
            data[f"tab_{i}"] = tab.save_state()
        _save_settings(data)

    def _restore_settings(self):
        data = _load_settings()
        if not data:
            return
        if not self._cub5.text() and "cub5" in data:
            self._cub5.setText(data["cub5"])
        if "active_tab" in data:
            idx = int(data["active_tab"])
            if 0 <= idx < self._tabs.count():
                self._tabs.setCurrentIndex(idx)
        for i, tab in enumerate(self._tab_instances):
            tab.restore_state(data.get(f"tab_{i}", {}))

    def closeEvent(self, event):
        try:
            self._save_settings()
        except Exception:
            pass  # never crash on close
        if self._process:
            try:
                self._process.kill()
                self._process.waitForFinished(3000)
            except Exception:
                pass
        super().closeEvent(event)


# ============================================================
# Entry point
# ============================================================

def main():
    app = QApplication(sys.argv)
    cub5 = sys.argv[1] if len(sys.argv) > 1 else ""
    window = RadiaApp(cub5)
    window.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
