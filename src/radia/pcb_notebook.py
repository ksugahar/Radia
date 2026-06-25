"""Notebook adapter for the Radia PCB PEEC panel."""

from __future__ import annotations

from .notebook_workbench import CommandWorkbench, NotebookFieldSpec, field_keys
from .pcb_design import PCBDesignSpec, PCB_SOLVERS


PCB_FIELD_SPECS = (
    NotebookFieldSpec("inp", "FastHenry .inp", section="Inputs", width="620px"),
    NotebookFieldSpec("freq_min", "f min", section="Sweep"),
    NotebookFieldSpec("freq_max", "f max", section="Sweep"),
    NotebookFieldSpec("n_freq", "N freq", "int", section="Sweep"),
    NotebookFieldSpec("solver_method", "Solver", "dropdown", PCB_SOLVERS, "Solver"),
    NotebookFieldSpec("spice_output", "SPICE output", section="Outputs", width="620px"),
)
PCB_NOTEBOOK_FIELD_ORDER = field_keys(PCB_FIELD_SPECS)


class PCBWorkbench(CommandWorkbench):
    title = "Radia PCB"
    field_specs = PCB_FIELD_SPECS
    section_order = ("Inputs", "Sweep", "Solver", "Outputs")

    def __init__(self, spec: PCBDesignSpec | None = None):
        super().__init__(spec or PCBDesignSpec(), run_root="runs/radia_pcb")


def display_pcb_workbench(spec: PCBDesignSpec | None = None):
    return PCBWorkbench(spec).display()
