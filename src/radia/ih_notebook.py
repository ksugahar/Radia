"""Jupyter widgets for the Radia IH design workbench."""

from __future__ import annotations

from pathlib import Path

from .ih_design import (
    FEM_SOLVERS,
    HEAT_SRC_SPATIAL,
    HEAT_SRC_UNIFORM,
    IHDesignSpec,
    IH_METHODS,
    MESH_TYPE_3D,
    MESH_TYPE_AXISYM,
    METHOD_PEEC_IND,
    PEEC_SOLVERS,
    THERMAL_PRESET_TO_CLI,
    TIME_SCHEME_TO_CLI,
)
from .notebook_workbench import (
    CommandWorkbench,
    NotebookFieldSpec,
    field_keys,
    shell_join,
)


IH_FIELD_SPECS = (
    NotebookFieldSpec("method", "Method", "dropdown", IH_METHODS, "Study", "520px"),
    NotebookFieldSpec(
        "solver",
        "Solver",
        "dropdown",
        tuple(PEEC_SOLVERS) + tuple(FEM_SOLVERS),
        "Study",
        "340px",
    ),
    NotebookFieldSpec("frequency", "f [Hz]", section="Study", width="180px"),
    NotebookFieldSpec("current", "I [A]", section="Study", width="180px"),
    NotebookFieldSpec("fes_order", "order", "int", section="Study", width="180px"),
    NotebookFieldSpec(
        "coil_material",
        "Coil material",
        "dropdown",
        ("Copper", "Aluminum", "Custom"),
        "Coil",
    ),
    NotebookFieldSpec("coil_sigma", "Coil sigma", section="Coil"),
    NotebookFieldSpec("peec_step", "Coil STEP", section="Coil", width="620px"),
    NotebookFieldSpec("peec_n_peri", "PEEC n peri", "int", section="Coil"),
    NotebookFieldSpec("peec_nwinc", "PEEC nwinc", "int", section="Coil"),
    NotebookFieldSpec("peec_nhinc", "PEEC nhinc", "int", section="Coil"),
    NotebookFieldSpec("coil_vol", "Coil .vol", section="Coil", width="620px"),
    NotebookFieldSpec("coil_source_name", "Source BC", section="Coil"),
    NotebookFieldSpec("coil_sink_name", "Sink BC", section="Coil"),
    NotebookFieldSpec("wp_vol", "Workpiece .vol", section="Workpiece", width="620px"),
    NotebookFieldSpec(
        "wp_material",
        "WP material",
        "dropdown",
        ("Steel (mu_r=100)", "Copper", "Aluminum", "Custom"),
        "Workpiece",
    ),
    NotebookFieldSpec("wp_sigma", "WP sigma", section="Workpiece"),
    NotebookFieldSpec("mu_r", "mu_r", section="Workpiece"),
    NotebookFieldSpec("half_thickness", "SIBC half t", section="Workpiece"),
    NotebookFieldSpec(
        "impedance_model",
        "Impedance",
        "dropdown",
        ("Linear SIBC", "Nonlinear ESIM (Karl iteration)"),
        "ESIM",
        "360px",
    ),
    NotebookFieldSpec("bh_file", "BH file", section="ESIM", width="620px"),
    NotebookFieldSpec("esim_max_iter", "max iter", "int", section="ESIM"),
    NotebookFieldSpec("esim_tol", "tol", section="ESIM"),
    NotebookFieldSpec("esim_per_panel", "per panel", "checkbox", section="ESIM"),
    NotebookFieldSpec("esim_anderson_m", "Anderson m", "int", section="ESIM"),
    NotebookFieldSpec("esim_relax", "relax", section="ESIM"),
    NotebookFieldSpec(
        "thermal_mesh_type",
        "Mesh type",
        "dropdown",
        (MESH_TYPE_3D, MESH_TYPE_AXISYM),
        "Thermal",
        "360px",
    ),
    NotebookFieldSpec(
        "heat_source",
        "Heat source",
        "dropdown",
        (HEAT_SRC_UNIFORM, HEAT_SRC_SPATIAL),
        "Thermal",
        "360px",
    ),
    NotebookFieldSpec("q_uniform", "q uniform", section="Thermal"),
    NotebookFieldSpec("qsurf_sol", "qsurf .sol", section="Thermal", width="620px"),
    NotebookFieldSpec("em_vol", "EM .vol", section="Thermal", width="620px"),
    NotebookFieldSpec("qsurf_order", "qsurf order", "int", section="Thermal"),
    NotebookFieldSpec("q_phi_average", "phi average", "checkbox", section="Thermal"),
    NotebookFieldSpec("n_phi_samples", "n phi", "int", section="Thermal"),
    NotebookFieldSpec("surface_label", "surface label", section="Thermal"),
    NotebookFieldSpec(
        "thermal_material",
        "Material",
        "dropdown",
        tuple(THERMAL_PRESET_TO_CLI),
        "Thermal",
    ),
    NotebookFieldSpec("override_kcprho", "custom rho/cp/k", "checkbox", section="Thermal"),
    NotebookFieldSpec("rho", "rho", section="Thermal"),
    NotebookFieldSpec("cp", "cp", section="Thermal"),
    NotebookFieldSpec("k", "k", section="Thermal"),
    NotebookFieldSpec("rotation_rpm", "rpm", "float", section="Motion"),
    NotebookFieldSpec("rotation_axis", "axis", "dropdown", ("x", "y", "z"), "Motion"),
    NotebookFieldSpec("h_conv", "h conv", section="Thermal BC"),
    NotebookFieldSpec("t_ext", "T ext", section="Thermal BC"),
    NotebookFieldSpec("emissivity", "emissivity", section="Thermal BC"),
    NotebookFieldSpec("t_init", "T init", section="Thermal BC"),
    NotebookFieldSpec(
        "time_scheme",
        "time scheme",
        "dropdown",
        tuple(TIME_SCHEME_TO_CLI),
        "Time",
        "260px",
    ),
    NotebookFieldSpec("dt", "dt", section="Time"),
    NotebookFieldSpec("t_end", "t end", section="Time"),
    NotebookFieldSpec(
        "linear_solver",
        "linear solver",
        "dropdown",
        ("sparsecholesky", "pardiso", "umfpack"),
        "Time",
    ),
    NotebookFieldSpec("probe_point", "probe point", section="Outputs", width="420px"),
    NotebookFieldSpec("csv_output", "CSV output", section="Outputs", width="620px"),
    NotebookFieldSpec("vtu_prefix", "VTU prefix", section="Outputs", width="620px"),
)
IH_NOTEBOOK_FIELD_ORDER = field_keys(IH_FIELD_SPECS)

__all__ = [
    "IHWorkbench",
    "IH_FIELD_SPECS",
    "IH_NOTEBOOK_FIELD_ORDER",
    "display_ih_workbench",
    "load_spec_from_paths",
    "shell_join",
]


class IHWorkbench(CommandWorkbench):
    """Browser-stable IH design panel for JupyterHub."""

    title = "IH"
    field_specs = IH_FIELD_SPECS
    section_order = (
        "Study",
        "Coil",
        "Workpiece",
        "ESIM",
        "Thermal",
        "Motion",
        "Thermal BC",
        "Time",
        "Outputs",
    )

    def __init__(self, spec: IHDesignSpec | None = None):
        super().__init__(spec or IHDesignSpec(), run_root="runs/radia_ih")


def display_ih_workbench(spec: IHDesignSpec | None = None):
    """Convenience entry point for notebooks."""

    return IHWorkbench(spec).display()


def load_spec_from_paths(
    *,
    method: str = METHOD_PEEC_IND,
    peec_step: str | Path = "",
    coil_vol: str | Path = "",
    wp_vol: str | Path = "",
) -> IHDesignSpec:
    """Small helper for notebook templates and examples."""

    return IHDesignSpec(
        method=method,
        peec_step=str(peec_step),
        coil_vol=str(coil_vol),
        wp_vol=str(wp_vol),
    )
