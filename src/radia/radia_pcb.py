"""
Radia PCB (PEEC) analysis window.

Usage:
    python -m radia.radia_pcb
    python radia_pcb.py
"""

import sys
import os

TITLE = "PCB"
REQUIRED_LABELS = []
OPTIONAL_LABELS = []
OPTIONAL_FILES = {"FastHenry .inp": "FastHenry (*.inp)"}

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from radia_gui_base import (
    ModePanel, AnalysisWindow, calc_script, json_output, run_app, _PYTHON,
)


class PCBPanel(ModePanel):

    def __init__(self, parent=None):
        super().__init__(parent)
        self._build_ui()

    def _build_ui(self):
        self.add_browse("inp_file", "FastHenry .inp:",
                        filter_str="FastHenry (*.inp);;All (*)")
        self.add_line("freq_min", "Freq min [Hz]:", "1e3")
        self.add_line("freq_max", "Freq max [Hz]:", "1e9")
        self.add_spin("n_freq", "N freq points:", 50, 1, 1000)
        self.add_combo("solver", "Solver:",
                       ["0 (LU)", "1 (BiCGSTAB)", "2 (HACApK)"])
        self.add_line("spice_output", "SPICE output:", "",
                      placeholder="e.g. output.cir")

    def is_runnable(self):
        return bool(self._widgets["inp_file"].text().strip())

    def build_command(self, vol_path):
        # PCB does not use .vol -- uses FastHenry .inp
        inp = self.val("inp_file")
        if not inp:
            raise ValueError("No FastHenry .inp file specified.")
        solver_id = self.val("solver").split()[0]
        cmd = [_PYTHON, calc_script("calc_pcb_peec.py"),
               "--inp", inp,
               "--freq-min", self.val("freq_min"),
               "--freq-max", self.val("freq_max"),
               "--n-freq", self.val("n_freq"),
               "--solver-method", solver_id,
               "--output", json_output(inp, "_pcb_peec")]
        spice = self.val("spice_output")
        if spice:
            cmd += ["--spice-output", spice]
        return cmd


class PCBWindow(AnalysisWindow):
    def __init__(self, vol_path=""):
        super().__init__("Radia - PCB (PEEC)", vol_path,
                         settings_key="pcb")
        panel = PCBPanel()
        self._set_panel(panel)
        self._restore_settings()
        self._update_run_state()


def main():
    run_app(PCBWindow)


if __name__ == "__main__":
    main()
