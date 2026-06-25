"""Notebook adapter for the Radia electromagnet panel."""

from __future__ import annotations

from .em_design import (
    EMDesignSpec,
    EM_METHODS,
    FEM_SOLVERS,
    MATERIALS,
    MSC_SOLVERS,
)
from .notebook_workbench import CommandWorkbench, NotebookFieldSpec, field_keys


EM_FIELD_SPECS = (
    NotebookFieldSpec("method", "Method", "dropdown", EM_METHODS, "Study"),
    NotebookFieldSpec("vol", "Mesh .vol", section="Inputs", width="620px"),
    NotebookFieldSpec("coil_script", "Coil script", section="Inputs", width="620px"),
    NotebookFieldSpec("material", "Material", "dropdown", MATERIALS, "Material"),
    NotebookFieldSpec("sigma", "sigma", section="Material"),
    NotebookFieldSpec("mu_r", "mu_r", section="Material"),
    NotebookFieldSpec("bh_file", "BH file", section="Material", width="620px"),
    NotebookFieldSpec("hys_file", "HYS file", section="Material", width="620px"),
    NotebookFieldSpec("fes_order", "FES order", "int", section="Solver"),
    NotebookFieldSpec("n_steps", "n steps", "int", section="Solver"),
    NotebookFieldSpec("max_iter", "max iter", "int", section="Solver"),
    NotebookFieldSpec("tol", "tol", section="Solver"),
    NotebookFieldSpec("relax", "relax", section="Solver"),
    NotebookFieldSpec("newton", "Newton", "checkbox", section="Solver"),
    NotebookFieldSpec("solver", "FEM solver", "dropdown", FEM_SOLVERS, "Solver"),
    NotebookFieldSpec("ima", "IMA", section="MSC"),
    NotebookFieldSpec("msc_solver", "MSC solver", "dropdown", MSC_SOLVERS, "MSC"),
    NotebookFieldSpec("demag_backend", "demag", "dropdown", ("hdiv",), "MSC"),
    NotebookFieldSpec("kelvin_mu_r", "mu_r", section="Kelvin"),
    NotebookFieldSpec("h0", "H0", section="Kelvin"),
    NotebookFieldSpec("field_axis", "axis", "dropdown", ("x", "y", "z"), "Kelvin"),
    NotebookFieldSpec("r_kelvin", "R Kelvin", section="Kelvin"),
    NotebookFieldSpec(
        "clebsch_geometry",
        "geometry",
        "dropdown",
        ("cylinder", "sphere"),
        "Clebsch",
    ),
    NotebookFieldSpec("clebsch_mu_r", "mu_r", section="Clebsch"),
    NotebookFieldSpec("clebsch_maxh", "maxh", section="Clebsch"),
    NotebookFieldSpec("clebsch_fes_order", "FES order", "int", section="Clebsch"),
)
EM_NOTEBOOK_FIELD_ORDER = field_keys(EM_FIELD_SPECS)


class EMWorkbench(CommandWorkbench):
    title = "Radia EM"
    field_specs = EM_FIELD_SPECS
    section_order = ("Study", "Inputs", "Material", "Solver", "MSC", "Kelvin", "Clebsch")

    def __init__(self, spec: EMDesignSpec | None = None):
        super().__init__(spec or EMDesignSpec(), run_root="runs/radia_em")


def display_em_workbench(spec: EMDesignSpec | None = None):
    return EMWorkbench(spec).display()
