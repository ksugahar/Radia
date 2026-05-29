"""
Radia PCB (PEEC) analysis window.

Pilot panel for the argparse-driven generator (POLICY 2026-05-30).  Every
widget below comes from `calc_pcb_peec.build_argparser()` via the
ModePanel.bind_argparser generator -- adding `--new-arg` to the calc
script automatically surfaces it in this panel without touching this
file.

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
    ModePanel, AnalysisWindow, calc_script, json_output, run_app,
)

# Display name <-> CLI value mapping for --solver-method (which argparse
# declares as type=int choices=[0,1,2]).  The generator's choice_map
# shows the readable name in the combo and emits the int on the CLI.
PCB_SOLVERS = [("LU", 0), ("BiCGSTAB", 1), ("HACApK", 2)]


class PCBPanel(ModePanel):

    def __init__(self, parent=None):
        super().__init__(parent)
        from radia.panels.calc_pcb_peec import build_argparser
        self.bind_argparser(
            build_argparser(),
            file_browse={
                "inp": ("FastHenry .inp:", "FastHenry (*.inp);;All (*)"),
            },
            choice_map={
                "solver_method": PCB_SOLVERS,
            },
            labels={
                "freq_min":     "Freq min [Hz]:",
                "freq_max":     "Freq max [Hz]:",
                "n_freq":       "N freq points:",
                "spice_output": "SPICE output:",
            },
        )

    def is_runnable(self):
        return bool(self._widgets["inp"].text().strip())

    def build_command(self, vol_path):
        inp = self.val("inp")
        if not inp:
            raise ValueError("No FastHenry .inp file specified.")
        return self.build_command_from_parser(
            vol_path=vol_path,
            script_path=calc_script("calc_pcb_peec.py"),
            vol_flag=None,          # PCB uses --inp, not --vol
            output_path=json_output(inp, "_pcb_peec"),
        )


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
