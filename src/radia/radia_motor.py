"""radia_motor panel -- Motor analysis (generator-driven, new-panel contract).

Layer-3 PySide6 panel for the two Stage-2 motor CLI scripts.  A composite
top-level panel holds the motor .vol Browse + an Analysis combo + a
QStackedWidget over the two sub-panels (Transient / Lamination), each a
ModePanel that bind_argparser()s its calc script's build_argparser().  Per
the new-panel contract every widget, the subprocess wiring, the .log capture
and the JSON result are derived from those argparsers -- no hand-assembled
argv, no orphan widgets (replaces the old hand-built 2-tab MotorWindow).

  Transient  : calc_motor_transient.py  (Lange-Henrotte-Hameyer nonlinear FE
               + circuit ODE, PM rotation, Arkkio torque)            -> --vol
  Lamination : calc_motor_lamination.py (Hollaus effective-material cell
               problem + global FE)                                  -> --vol

Per the 4-Layer Architecture this Layer-3 file launches Layer-4 calc_*.py via
subprocess and does NOT import radia or ngsolve directly.

Usage:
    python radia_motor.py [mesh.vol]
    python -m radia.radia_motor [mesh.vol]
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from PySide6.QtWidgets import QStackedWidget

from radia_gui_base import (
    ModePanel, AnalysisWindow, calc_script, json_output, run_app,
)

TITLE = "Radia - Motor"

# argparser dest owned by the composite (motor .vol -> --vol via wp_vol) /
# by calc_main (--output via build_command).  Per-mode knobs (transient
# --method, lamination --mode) stay as sub-panel widgets.
_SKIP = ("vol", "output")


def _transient_argparser():
    from radia.panels.calc_motor_transient import build_argparser
    return build_argparser()


def _lamination_argparser():
    from radia.panels.calc_motor_lamination import build_argparser
    return build_argparser()


class _TransientPanel(ModePanel):
    """Nonlinear FE + circuit ODE transient (Lange-Henrotte-Hameyer)."""
    def __init__(self, parent=None):
        super().__init__(parent)
        self.bind_argparser(_transient_argparser(), skip=_SKIP)

    def build_command(self, vol_path):
        return self.build_command_from_parser(
            vol_path=vol_path, vol_flag="--vol",
            script_path=calc_script("calc_motor_transient.py"),
            output_path=json_output(vol_path, "_motor_transient"),
        )


class _LaminationPanel(ModePanel):
    """Hollaus effective-material homogenization (cell / global / full)."""
    def __init__(self, parent=None):
        super().__init__(parent)
        self.bind_argparser(_lamination_argparser(), skip=_SKIP)

    def build_command(self, vol_path):
        return self.build_command_from_parser(
            vol_path=vol_path, vol_flag="--vol",
            script_path=calc_script("calc_motor_lamination.py"),
            output_path=json_output(vol_path, "_motor_lamination"),
        )


_ANALYSIS_FACTORIES = [
    ("Transient", _TransientPanel),
    ("Lamination", _LaminationPanel),
]
# modes whose calc requires a .vol (transient needs the mesh; lamination
# cell-mode does not -- it solves a 1D cell problem with no .vol).
_NEEDS_VOL = {"Transient"}


class MotorPanel(ModePanel):
    """Composite: motor .vol Browse + Analysis combo + QStackedWidget."""
    def __init__(self, parent=None):
        super().__init__(parent)
        self.add_browse("wp_vol", "Motor .vol:",
                        filter_str="Netgen Vol (*.vol);;All (*)")
        self._analysis_combo = self.add_combo(
            "analysis", "Analysis:", [name for name, _ in _ANALYSIS_FACTORIES])
        self._sub_panels = {name: factory()
                            for name, factory in _ANALYSIS_FACTORIES}
        self._stack = QStackedWidget()
        for name, _ in _ANALYSIS_FACTORIES:
            self._stack.addWidget(self._sub_panels[name])
        self._form.addRow(self._stack)
        self._analysis_combo.currentTextChanged.connect(self._on_analysis_changed)
        self._on_analysis_changed(_ANALYSIS_FACTORIES[0][0])

    def _current_name(self):
        return self._analysis_combo.currentText()

    def _current_sub(self):
        return self._sub_panels[self._current_name()]

    def _on_analysis_changed(self, text):
        if text in self._sub_panels:
            self._stack.setCurrentWidget(self._sub_panels[text])
        cb = getattr(self, "validationChanged", None)
        if callable(cb):
            cb()

    def is_runnable(self):
        if self._current_name() in _NEEDS_VOL and not self.val("wp_vol").strip():
            return False
        return self._current_sub().is_runnable()

    def build_command(self, vol_path):
        return self._current_sub().build_command(vol_path)

    def wp_vol_path(self):
        return self.val("wp_vol") if "wp_vol" in self._widgets else ""

    def save_state(self):
        state = super().save_state()
        for name, sub in self._sub_panels.items():
            for k, v in sub.save_state().items():
                state[f"{name}/{k}"] = v
        return state

    def restore_state(self, state):
        if not state:
            return
        super().restore_state({k: v for k, v in state.items() if "/" not in k})
        for name, sub in self._sub_panels.items():
            prefix = f"{name}/"
            sub.restore_state({k[len(prefix):]: v for k, v in state.items()
                               if k.startswith(prefix)})
        self._on_analysis_changed(self._analysis_combo.currentText())


class MotorWindow(AnalysisWindow):
    def __init__(self, vol_path=""):
        super().__init__(TITLE, vol_path, settings_key="motor")
        panel = MotorPanel()
        self._set_panel(panel)
        # Restore last session FIRST, then let an explicit launcher vol_path
        # override the restored wp_vol (the radia_em / radia_ih restore-order
        # bug class -- restore before setText, never the reverse).
        self._restore_settings()
        if vol_path and "wp_vol" in panel._widgets:
            panel._widgets["wp_vol"].setText(self.display_path(vol_path))


def main():
    run_app(MotorWindow, sys.argv)


if __name__ == "__main__":
    main()
