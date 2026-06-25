"""UI-neutral PCB PEEC design state."""

from __future__ import annotations

from dataclasses import dataclass
import sys

from .panel_design_common import append_value, calc_script, json_output


PCB_SOLVERS = (("LU", 0), ("BiCGSTAB", 1), ("HACApK", 2))


@dataclass(slots=True)
class PCBDesignSpec:
    inp: str = ""
    freq_min: str = "1e3"
    freq_max: str = "1e9"
    n_freq: int = 50
    solver_method: int = 0
    spice_output: str = ""

    def visible_fields(self) -> set[str]:
        return {
            "inp", "freq_min", "freq_max", "n_freq",
            "solver_method", "spice_output",
        }

    def missing_required_inputs(self) -> list[str]:
        return [] if self.inp.strip() else ["FastHenry .inp"]

    def is_runnable(self) -> bool:
        return not self.missing_required_inputs()

    def build_command(self, *, python: str | None = None, panels_dir=None) -> list[str]:
        if not self.inp:
            raise ValueError("No FastHenry .inp file specified.")
        cmd = [
            python or sys.executable,
            calc_script("calc_pcb_peec.py", panels_dir),
            "--inp", self.inp,
            "--freq-min", str(self.freq_min),
            "--freq-max", str(self.freq_max),
            "--n-freq", str(self.n_freq),
            "--solver-method", str(self.solver_method),
            "--output", json_output(self.inp, "_pcb_peec"),
        ]
        append_value(cmd, "--spice-output", self.spice_output)
        return cmd
