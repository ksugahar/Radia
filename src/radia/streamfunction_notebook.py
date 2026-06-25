"""Notebook adapter for the Radia stream-function panel."""

from __future__ import annotations

from .notebook_workbench import CommandWorkbench, NotebookFieldSpec, field_keys
from .streamfunction_design import (
    STREAMFUNCTION_METHODS,
    StreamFunctionDesignSpec,
)


STREAMFUNCTION_FIELD_SPECS = (
    NotebookFieldSpec("method", "Method", "dropdown", STREAMFUNCTION_METHODS, "Study"),
    NotebookFieldSpec("coil_vol", "Coil .vol", section="Surface", width="620px"),
    NotebookFieldSpec("eval_vol", "Eval .vol", section="Surface", width="620px"),
    NotebookFieldSpec("target_cf", "target CF", section="Target", width="420px"),
    NotebookFieldSpec("target_harmonic", "harmonic", section="Target", width="420px"),
    NotebookFieldSpec("harmonic_lmax", "lmax", "int", section="Target"),
    NotebookFieldSpec("order", "psi order", "int", section="Solver"),
    NotebookFieldSpec(
        "regularize",
        "regularize",
        "dropdown",
        ("l2", "h1", "inductance"),
        "Solver",
    ),
    NotebookFieldSpec("confine", "confine", "dropdown", ("abe", "off", "on"), "Solver"),
    NotebookFieldSpec("alpha", "alpha", section="Solver"),
    NotebookFieldSpec("eval_max", "eval max", "int", section="Solver"),
    NotebookFieldSpec(
        "pareto_lever",
        "lever",
        "dropdown",
        ("alpha", "linf", "geometry"),
        "Pareto",
    ),
    NotebookFieldSpec("alpha_min", "alpha min", section="Pareto"),
    NotebookFieldSpec("alpha_max", "alpha max", section="Pareto"),
    NotebookFieldSpec("n_alpha", "N alpha", "int", section="Pareto"),
    NotebookFieldSpec("linf_iter", "linf iter", "int", section="Pareto"),
    NotebookFieldSpec("geom_scale_min", "scale min", section="Pareto"),
    NotebookFieldSpec("geom_scale_max", "scale max", section="Pareto"),
    NotebookFieldSpec("nlevels", "levels", "int", section="Manufacture"),
    NotebookFieldSpec("optimize_levels", "opt levels", "checkbox", section="Manufacture"),
    NotebookFieldSpec("greedy_turns", "greedy turns", section="Manufacture"),
    NotebookFieldSpec(
        "greedy_dict",
        "greedy dict",
        "dropdown",
        ("contour", "pin", "bubble"),
        "Manufacture",
    ),
    NotebookFieldSpec("greedy_connector_weight", "connector w", section="Manufacture"),
    NotebookFieldSpec("greedy_target", "greedy target", section="Manufacture"),
    NotebookFieldSpec("greedy_plot", "greedy plot", section="Manufacture", width="620px"),
    NotebookFieldSpec("pin_tiling", "pin tiling", "checkbox", section="Manufacture"),
    NotebookFieldSpec("pin_tiling_pins", "pins", "int", section="Manufacture"),
    NotebookFieldSpec("pin_tiling_frac", "pin frac", section="Manufacture"),
    NotebookFieldSpec("target_inductance", "target L", section="Manufacture"),
    NotebookFieldSpec("resonance_cap", "C tank", section="Manufacture"),
    NotebookFieldSpec("nlevels_max", "levels max", "int", section="Manufacture"),
    NotebookFieldSpec("contour_sub", "contour sub", "int", section="Manufacture"),
    NotebookFieldSpec("chain", "chain", "dropdown", ("field_aware", "nn"), "Manufacture"),
    NotebookFieldSpec("chain_ncut", "chain cuts", "int", section="Manufacture"),
    NotebookFieldSpec("chain_passes", "passes", "int", section="Manufacture"),
    NotebookFieldSpec("distort", "distort", "checkbox", section="Manufacture"),
    NotebookFieldSpec("distort_grid", "dist grid", "int", section="Manufacture"),
    NotebookFieldSpec("distort_iter", "dist iter", "int", section="Manufacture"),
    NotebookFieldSpec("step_output", "STEP output", section="Outputs", width="620px"),
    NotebookFieldSpec("peec", "PEEC", "checkbox", section="Outputs"),
    NotebookFieldSpec("wire_diam", "wire diam", section="Outputs"),
    NotebookFieldSpec("peec_freq", "PEEC freq", section="Outputs"),
    NotebookFieldSpec("flux_plot", "flux plot", section="Outputs", width="620px"),
    NotebookFieldSpec("flux_plane", "flux plane", "dropdown", ("x", "y", "z"), "Outputs"),
    NotebookFieldSpec("steps_plot", "steps plot", section="Outputs", width="620px"),
    NotebookFieldSpec("iron_vol", "Iron .vol", section="Material-aware", width="620px"),
    NotebookFieldSpec("mu_r", "mu_r", section="Material-aware"),
    NotebookFieldSpec("iron_mat", "iron mat", section="Material-aware"),
    NotebookFieldSpec("iron_exact_source", "exact source", "checkbox", section="Material-aware"),
    NotebookFieldSpec("iron_quad_order", "quad order", section="Material-aware"),
    NotebookFieldSpec("shield_vol", "Shield .vol", section="Shield", width="620px"),
    NotebookFieldSpec("shield_eval_vol", "Shield eval .vol", section="Shield", width="620px"),
    NotebookFieldSpec("shield_weight", "shield w", section="Shield"),
    NotebookFieldSpec("volume_vol", "Conductor .vol", section="Volume 3D", width="620px"),
    NotebookFieldSpec("target_bz", "target Bz", section="Volume 3D"),
    NotebookFieldSpec("n_leaves", "leaves", "int", section="Volume 3D"),
    NotebookFieldSpec("fes_order", "lambda order", "int", section="Volume 3D"),
    NotebookFieldSpec("aca_eps", "ACA eps", section="Volume 3D"),
    NotebookFieldSpec("rings_per_span", "rings/span", "int", section="Volume 3D"),
    NotebookFieldSpec("n_targets", "targets", "int", section="Volume 3D"),
    NotebookFieldSpec("nphi", "nphi", "int", section="Volume 3D"),
    NotebookFieldSpec("nz", "nz", "int", section="Volume 3D"),
    NotebookFieldSpec("n_threads", "threads", "int", section="Volume 3D"),
)
STREAMFUNCTION_NOTEBOOK_FIELD_ORDER = field_keys(STREAMFUNCTION_FIELD_SPECS)


class StreamFunctionWorkbench(CommandWorkbench):
    title = "Radia Stream Function"
    field_specs = STREAMFUNCTION_FIELD_SPECS
    section_order = (
        "Study",
        "Surface",
        "Target",
        "Solver",
        "Pareto",
        "Manufacture",
        "Material-aware",
        "Shield",
        "Volume 3D",
        "Outputs",
    )

    def __init__(self, spec: StreamFunctionDesignSpec | None = None):
        super().__init__(
            spec or StreamFunctionDesignSpec(),
            run_root="runs/radia_streamfunction",
        )


def display_streamfunction_workbench(spec: StreamFunctionDesignSpec | None = None):
    return StreamFunctionWorkbench(spec).display()
