"""Notebook adapter for the Radia motor panel."""

from __future__ import annotations

from .motor_design import (
    LAMINATION_DRIVES,
    LAMINATION_MODES,
    LINEAR_SOLVERS,
    MOTOR_ANALYSES,
    MotorDesignSpec,
    TRANSIENT_METHODS,
)
from .notebook_workbench import CommandWorkbench, NotebookFieldSpec, field_keys


MOTOR_FIELD_SPECS = (
    NotebookFieldSpec("analysis", "Analysis", "dropdown", MOTOR_ANALYSES, "Study"),
    NotebookFieldSpec("vol", "Motor .vol", section="Inputs", width="620px"),
    NotebookFieldSpec("method", "Method", "dropdown", TRANSIENT_METHODS, "Transient"),
    NotebookFieldSpec("fes_order", "FES order", "int", section="Solver"),
    NotebookFieldSpec("linear_solver", "linear solver", "dropdown", LINEAR_SOLVERS, "Solver"),
    NotebookFieldSpec("nbr_phases", "phases", "int", section="Transient"),
    NotebookFieldSpec("n_turns_per_slot", "turns/slot", "int", section="Transient"),
    NotebookFieldSpec("slot_area", "slot area", section="Transient"),
    NotebookFieldSpec("stack_length", "stack length", section="Transient"),
    NotebookFieldSpec("n_pole_pairs", "pole pairs", "int", section="Transient"),
    NotebookFieldSpec("r_airgap_mid", "airgap r", section="Transient"),
    NotebookFieldSpec("r_phase", "R phase", section="Circuit"),
    NotebookFieldSpec("l_endwinding", "L end", section="Circuit"),
    NotebookFieldSpec("pm_br_value", "PM Br", section="Circuit"),
    NotebookFieldSpec("j_inertia", "J inertia", section="Mechanics"),
    NotebookFieldSpec("b_viscous", "B visc", section="Mechanics"),
    NotebookFieldSpec("t_load", "T load", section="Mechanics"),
    NotebookFieldSpec("v_amp", "V amp", section="Circuit"),
    NotebookFieldSpec("v_freq", "V freq", section="Circuit"),
    NotebookFieldSpec("theta_init", "theta0", section="Time"),
    NotebookFieldSpec("omega_init", "omega0", section="Time"),
    NotebookFieldSpec("t_end", "t end", section="Time"),
    NotebookFieldSpec("dt_fe", "dt FE", section="Time"),
    NotebookFieldSpec("dt_circ", "dt circ", section="Time"),
    NotebookFieldSpec("n_steps_per_fe", "steps/FE", "int", section="Time"),
    NotebookFieldSpec("lamination_mode", "Mode", "dropdown", LAMINATION_MODES, "Lamination"),
    NotebookFieldSpec("d_iron", "d iron", section="Lamination"),
    NotebookFieldSpec("d_ins", "d ins", section="Lamination"),
    NotebookFieldSpec("sigma", "sigma", section="Lamination"),
    NotebookFieldSpec("mu_r_iron", "mu_r iron", section="Lamination"),
    NotebookFieldSpec("b_list", "B list", section="Lamination"),
    NotebookFieldSpec("freq_list", "freq list", section="Lamination"),
    NotebookFieldSpec("cell_n_elements", "cell N", "int", section="Lamination"),
    NotebookFieldSpec("drive", "drive", "dropdown", LAMINATION_DRIVES, "Lamination"),
    NotebookFieldSpec("j_s_iron", "J_s iron", section="Lamination"),
    NotebookFieldSpec("em_table", "EM table", section="Lamination", width="620px"),
    NotebookFieldSpec("h_amplitude", "H amp", section="Lamination"),
    NotebookFieldSpec("freq", "freq", section="Lamination"),
)
MOTOR_NOTEBOOK_FIELD_ORDER = field_keys(MOTOR_FIELD_SPECS)


class MotorWorkbench(CommandWorkbench):
    title = "Radia Motor"
    field_specs = MOTOR_FIELD_SPECS
    section_order = (
        "Study",
        "Inputs",
        "Solver",
        "Transient",
        "Circuit",
        "Mechanics",
        "Time",
        "Lamination",
    )

    def __init__(self, spec: MotorDesignSpec | None = None):
        super().__init__(spec or MotorDesignSpec(), run_root="runs/radia_motor")


def display_motor_workbench(spec: MotorDesignSpec | None = None):
    return MotorWorkbench(spec).display()
