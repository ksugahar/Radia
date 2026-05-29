"""Radia Motor analysis window (Stage-3 PySide6 wrapper).

Two-tab GUI wrapping the Stage-2 motor CLI scripts:

  Tab 1: Transient (Lange-Henrotte-Hameyer 2009)
         calc_motor_transient.py -- nonlinear FE + circuit ODE
         with PM rotation and Arkkio torque.

  Tab 2: Lamination (Hollaus Effective Material)
         calc_motor_lamination.py -- cell-problem homogenization
         + global FE with effective material.

Both tabs share a file picker for the .vol mesh; each has its own
parameter form and Run button.  The shared output panel shows CLI
stderr (progress) and final JSON.

Per the 4-Layer Panel Architecture (CLAUDE.md), this Layer 3 file
launches Layer 4 calc_*.py via subprocess and does NOT import radia
or ngsolve directly.

Usage:
    python radia_motor.py [mesh.vol]
    python -m radia.radia_motor [mesh.vol]
"""

from __future__ import annotations

import json
import os
import sys

from PySide6.QtCore import QProcess, Qt
from PySide6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QFormLayout, QTabWidget, QPushButton, QLineEdit, QFileDialog,
    QLabel, QTextEdit, QGroupBox,
)
# Wheel-guarded combos/spins: ignore mouse-wheel so scrolling the form
# does not silently change a value (see radia_gui_base, kubota 2026-05-29).
from radia_gui_base import (
    _NoWheelDoubleSpinBox as QDoubleSpinBox,
    _NoWheelSpinBox as QSpinBox,
    _NoWheelComboBox as QComboBox,
)


# ---------------------------------------------------------------- helpers

_PYTHON = sys.executable
_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_PANELS_DIR = os.path.join(_THIS_DIR, "panels")


def _calc_script(name):
    p = os.path.join(_PANELS_DIR, name)
    if not os.path.exists(p):
        raise FileNotFoundError(f"calc script not found: {p}")
    return p


def _double(min_val, max_val, default, step=0.1, decimals=4):
    sb = QDoubleSpinBox()
    sb.setRange(min_val, max_val)
    sb.setValue(default)
    sb.setSingleStep(step)
    sb.setDecimals(decimals)
    return sb


def _spin(min_val, max_val, default):
    sb = QSpinBox()
    sb.setRange(min_val, max_val)
    sb.setValue(default)
    return sb


# ------------------------------------------------------------------ Tab 1

class TransientTab(QWidget):
    """Lange-Henrotte-Hameyer 2009 transient solver tab."""

    def __init__(self, get_vol):
        super().__init__()
        self._get_vol = get_vol
        layout = QVBoxLayout()
        self.setLayout(layout)

        # Geometry group
        geo = QGroupBox("Geometry")
        gf = QFormLayout()
        self.nbr_phases = _spin(1, 6, 3)
        self.n_pole_pairs = _spin(1, 32, 4)
        self.n_turns = _spin(1, 1000, 100)
        self.slot_area = _double(1e-6, 1e-2, 1e-4, step=1e-5, decimals=6)
        self.stack_length = _double(0.001, 1.0, 0.05, step=0.005, decimals=4)
        self.r_airgap_mid = _double(0.001, 1.0, 0.050, step=0.001, decimals=4)
        gf.addRow("Phases", self.nbr_phases)
        gf.addRow("Pole pairs", self.n_pole_pairs)
        gf.addRow("Turns / slot", self.n_turns)
        gf.addRow("Slot area [m²]", self.slot_area)
        gf.addRow("Stack length [m]", self.stack_length)
        gf.addRow("Air-gap mid R [m]", self.r_airgap_mid)
        geo.setLayout(gf)
        layout.addWidget(geo)

        # Circuit group
        circ = QGroupBox("Circuit + Mechanical")
        cf = QFormLayout()
        self.R_phase = _double(0.001, 100, 0.3872, step=0.01, decimals=4)
        self.L_endwinding = _double(0.0, 1.0, 0.0, step=1e-4, decimals=5)
        self.J_inertia = _double(1e-6, 1.0, 1e-3, step=1e-4, decimals=6)
        self.T_load = _double(-100, 100, 0.0, step=0.1, decimals=3)
        self.pm_Br = _double(0.0, 2.0, 1.2, step=0.1, decimals=3)
        cf.addRow("Phase R [Ω]", self.R_phase)
        cf.addRow("End-winding L [H]", self.L_endwinding)
        cf.addRow("Rotor inertia J [kg·m²]", self.J_inertia)
        cf.addRow("Load torque [N·m]", self.T_load)
        cf.addRow("PM Br [T] (0=no PM)", self.pm_Br)
        circ.setLayout(cf)
        layout.addWidget(circ)

        # Excitation group
        src = QGroupBox("Excitation + Time integration")
        sf = QFormLayout()
        self.v_amp = _double(0.0, 10000.0, 50.0, step=10, decimals=2)
        self.v_freq = _double(0.0, 100000.0, 100.0, step=10, decimals=2)
        self.omega_init = _double(-10000, 10000, 62.83, step=10, decimals=2)
        self.t_end = _double(1e-5, 100.0, 0.5e-3, step=1e-4, decimals=6)
        self.dt_FE = _double(1e-6, 1.0, 1e-4, step=1e-5, decimals=6)
        self.dt_circ = _double(1e-7, 1.0, 1e-5, step=1e-6, decimals=7)
        self.n_steps_per_FE = _spin(1, 1000, 10)
        sf.addRow("V amplitude [V]", self.v_amp)
        sf.addRow("V frequency [Hz]", self.v_freq)
        sf.addRow("Initial ω [rad/s]", self.omega_init)
        sf.addRow("t_end [s]", self.t_end)
        sf.addRow("Δt_FE [s]", self.dt_FE)
        sf.addRow("Δt_circ [s]", self.dt_circ)
        sf.addRow("Circuit steps / FE", self.n_steps_per_FE)
        src.setLayout(sf)
        layout.addWidget(src)

        layout.addStretch()

    def build_cmd(self, output_json):
        vol = self._get_vol()
        if not vol:
            return None
        return [
            _PYTHON, _calc_script("calc_motor_transient.py"),
            "--vol", vol,
            "--method", "linearization",
            "--nbr-phases", str(self.nbr_phases.value()),
            "--n-pole-pairs", str(self.n_pole_pairs.value()),
            "--n-turns-per-slot", str(self.n_turns.value()),
            "--slot-area", str(self.slot_area.value()),
            "--stack-length", str(self.stack_length.value()),
            "--r-airgap-mid", str(self.r_airgap_mid.value()),
            "--R-phase", str(self.R_phase.value()),
            "--L-endwinding", str(self.L_endwinding.value()),
            "--J-inertia", str(self.J_inertia.value()),
            "--T-load", str(self.T_load.value()),
            "--pm-Br-value", str(self.pm_Br.value()),
            "--v-amp", str(self.v_amp.value()),
            "--v-freq", str(self.v_freq.value()),
            "--omega-init", str(self.omega_init.value()),
            "--t-end", str(self.t_end.value()),
            "--dt-FE", str(self.dt_FE.value()),
            "--dt-circ", str(self.dt_circ.value()),
            "--n-steps-per-FE", str(self.n_steps_per_FE.value()),
            "--output", output_json,
        ]


# ------------------------------------------------------------------ Tab 2

class LaminationTab(QWidget):
    """Hollaus Effective Material homogenization tab."""

    def __init__(self, get_vol):
        super().__init__()
        self._get_vol = get_vol
        layout = QVBoxLayout()
        self.setLayout(layout)

        # Mode group
        mode = QGroupBox("Mode")
        mf = QFormLayout()
        self.mode = QComboBox()
        self.mode.addItems(["cell (offline only)", "global (need em-table)",
                             "full (cell + global)"])
        self.mode.setCurrentIndex(2)
        mf.addRow("Mode", self.mode)
        mode.setLayout(mf)
        layout.addWidget(mode)

        # Cell-problem group
        cell = QGroupBox("Cell problem (offline)")
        ccf = QFormLayout()
        self.d_iron = _double(1e-5, 5e-3, 0.35e-3, step=1e-5, decimals=6)
        self.d_ins = _double(0.0, 1e-3, 0.05e-3, step=1e-6, decimals=7)
        self.sigma = _double(1e3, 1e8, 4e6, step=1e5, decimals=1)
        self.mu_r_iron = _double(1.0, 100000, 5000, step=100, decimals=1)
        self.B_list = QLineEdit("1.0")
        self.freq_list = QLineEdit("50,500,5000")
        self.cell_n_elements = _spin(10, 1000, 60)
        ccf.addRow("Iron thickness d [m]", self.d_iron)
        ccf.addRow("Insulation thickness [m]", self.d_ins)
        ccf.addRow("σ_iron [S/m]", self.sigma)
        ccf.addRow("μ_r iron (linear)", self.mu_r_iron)
        ccf.addRow("B amplitudes [T] (csv)", self.B_list)
        ccf.addRow("Frequencies [Hz] (csv)", self.freq_list)
        ccf.addRow("Cell mesh elements", self.cell_n_elements)
        cell.setLayout(ccf)
        layout.addWidget(cell)

        # Global FE group
        glob = QGroupBox("Global FE (online)")
        gf = QFormLayout()
        self.H_amp = _double(0.0, 1e6, 1000.0, step=10, decimals=2)
        self.freq = _double(0.1, 1e6, 1000.0, step=10, decimals=2)
        self.fes_order = _spin(1, 5, 2)
        gf.addRow("H amplitude [A/m]", self.H_amp)
        gf.addRow("Operating freq [Hz]", self.freq)
        gf.addRow("FES order", self.fes_order)
        glob.setLayout(gf)
        layout.addWidget(glob)

        layout.addStretch()

    def build_cmd(self, output_json):
        mode_text = self.mode.currentText().split(" ")[0]
        cmd = [
            _PYTHON, _calc_script("calc_motor_lamination.py"),
            "--mode", mode_text,
            "--d-iron", str(self.d_iron.value()),
            "--d-ins", str(self.d_ins.value()),
            "--sigma", str(self.sigma.value()),
            "--mu-r-iron", str(self.mu_r_iron.value()),
            "--B-list", self.B_list.text(),
            "--freq-list", self.freq_list.text(),
            "--cell-n-elements", str(self.cell_n_elements.value()),
            "--H-amplitude", str(self.H_amp.value()),
            "--freq", str(self.freq.value()),
            "--fes-order", str(self.fes_order.value()),
            "--output", output_json,
        ]
        if mode_text in ("global", "full"):
            vol = self._get_vol()
            if not vol:
                return None
            cmd += ["--vol", vol]
        return cmd


# ----------------------------------------------------------------- Window

class MotorWindow(QMainWindow):
    """Top-level window with vol picker + two tabs + output."""

    def __init__(self, vol_path=""):
        super().__init__()
        self.setWindowTitle("Radia -- Motor (Lange-Henrotte-Hameyer + Hollaus)")
        self.resize(720, 880)

        central = QWidget()
        self.setCentralWidget(central)
        layout = QVBoxLayout()
        central.setLayout(layout)

        # File picker bar
        bar = QHBoxLayout()
        bar.addWidget(QLabel(".vol mesh:"))
        self.vol_edit = QLineEdit(vol_path)
        bar.addWidget(self.vol_edit, 1)
        browse = QPushButton("Browse…")
        browse.clicked.connect(self._browse)
        bar.addWidget(browse)
        layout.addLayout(bar)

        # Tabs
        self.tabs = QTabWidget()
        self.transient = TransientTab(get_vol=self.get_vol)
        self.lamination = LaminationTab(get_vol=self.get_vol)
        self.tabs.addTab(self.transient, "Transient (Lange-HH)")
        self.tabs.addTab(self.lamination, "Lamination (Hollaus EM)")
        layout.addWidget(self.tabs)

        # Run button + output
        run_bar = QHBoxLayout()
        self.run_btn = QPushButton("Run")
        self.run_btn.clicked.connect(self._run)
        run_bar.addWidget(self.run_btn)
        self.stop_btn = QPushButton("Stop")
        self.stop_btn.clicked.connect(self._stop)
        self.stop_btn.setEnabled(False)
        run_bar.addWidget(self.stop_btn)
        run_bar.addStretch()
        layout.addLayout(run_bar)

        self.output = QTextEdit()
        self.output.setReadOnly(True)
        self.output.setStyleSheet("QTextEdit { font-family: Consolas, monospace; "
                                  "font-size: 10pt; }")
        layout.addWidget(self.output, 1)

        self._proc = None

    def get_vol(self):
        v = self.vol_edit.text().strip()
        if not v:
            self.output.append(
                "[ERROR] please pick a .vol mesh first (Browse...)"
            )
            return ""
        if not os.path.exists(v):
            self.output.append(f"[ERROR] .vol not found: {v}")
            return ""
        return v

    def _browse(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Choose Netgen .vol mesh", "",
            "Netgen Vol (*.vol);;All files (*)")
        if path:
            self.vol_edit.setText(path)

    def _run(self):
        # Determine which tab is active
        tab = self.tabs.currentWidget()
        # Output JSON in same dir as mesh, otherwise C:/temp
        vol = self.vol_edit.text().strip()
        out_dir = (os.path.dirname(vol) if vol and os.path.exists(vol)
                   else "C:/temp")
        os.makedirs(out_dir, exist_ok=True)
        if tab is self.transient:
            out_json = os.path.join(out_dir, "motor_transient_result.json")
        else:
            out_json = os.path.join(out_dir, "motor_lamination_result.json")

        cmd = tab.build_cmd(out_json)
        if cmd is None:
            return

        self.output.clear()
        self.output.append(f"[CMD] {' '.join(cmd)}")
        self.run_btn.setEnabled(False)
        self.stop_btn.setEnabled(True)

        self._proc = QProcess(self)
        self._proc.setProcessChannelMode(QProcess.MergedChannels)
        self._proc.readyReadStandardOutput.connect(self._on_stdout)
        self._proc.finished.connect(self._on_finished)
        env = os.environ.copy()
        env["PYTHONIOENCODING"] = "utf-8"
        from PySide6.QtCore import QProcessEnvironment
        qenv = QProcessEnvironment.systemEnvironment()
        for k, v in env.items():
            qenv.insert(k, v)
        self._proc.setProcessEnvironment(qenv)
        self._proc.start(cmd[0], cmd[1:])

    def _on_stdout(self):
        if not self._proc:
            return
        data = bytes(self._proc.readAllStandardOutput()).decode(
            "utf-8", errors="replace")
        self.output.moveCursor(self.output.textCursor().End)
        self.output.insertPlainText(data)
        sb = self.output.verticalScrollBar()
        sb.setValue(sb.maximum())

    def _on_finished(self, exit_code, _exit_status):
        self.output.append(f"\n[DONE] exit code = {exit_code}")
        self.run_btn.setEnabled(True)
        self.stop_btn.setEnabled(False)
        self._proc = None

    def _stop(self):
        if self._proc is not None:
            self._proc.terminate()
            self.output.append("[STOPPED]")


def main():
    app = QApplication(sys.argv)
    vol = sys.argv[1] if len(sys.argv) > 1 else ""
    win = MotorWindow(vol_path=vol)
    win.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
