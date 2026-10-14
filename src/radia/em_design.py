"""UI-neutral electromagnet panel design state."""

from __future__ import annotations

from dataclasses import dataclass
import sys

from .panel_design_common import (
    append_switch,
    append_value,
    calc_script,
    json_output,
    msh_output,
)


METHOD_OMEGA = "Omega"
METHOD_APHI = "A-Phi"
METHOD_HDIV = "HDiv-VIM"
METHOD_KELVIN_BENCH = "Kelvin Benchmark"
METHOD_CLEBSCH = "Clebsch hodograph"
EM_METHODS = (
    METHOD_OMEGA,
    METHOD_APHI,
    METHOD_HDIV,
    METHOD_KELVIN_BENCH,
    METHOD_CLEBSCH,
)

FEM_SOLVERS = ("auto", "sparsecholesky", "bddc_ams", "bddc", "iccg")
MATERIALS = ("steel", "copper", "aluminum", "elf_steel", "linear", "hysteresis")
# Options of the Omega/A-Phi FEM iteration that the HDiv-VIM route cannot apply.
_FEM_ONLY_OPTIONS = ("n_steps", "max_iter", "tol", "relax", "newton", "solver")


@dataclass(slots=True)
class EMDesignSpec:
    method: str = METHOD_OMEGA
    vol: str = ""
    coil_script: str = ""
    material: str = "steel"
    sigma: str = ""
    mu_r: str = ""
    bh_file: str = ""
    hys_file: str = ""
    fes_order: int = 1
    n_steps: int = 1
    max_iter: int = 30
    tol: str = "1e-3"
    relax: str = "0.3"
    newton: bool = False
    solver: str = "auto"
    ima: str = ""
    demag_backend: str = "hdiv"
    kelvin_mu_r: str = "100"
    h0: str = "1.0"
    field_axis: str = "z"
    r_kelvin: str = "0.20"
    clebsch_geometry: str = "cylinder"
    clebsch_mu_r: str = "1000"
    clebsch_maxh: str = "0.08"
    clebsch_fes_order: int = 3

    def visible_fields(self) -> set[str]:
        fields = {"method"}
        if self.method in (METHOD_OMEGA, METHOD_APHI):
            fields.update({
                "vol", "coil_script", "material", "sigma", "mu_r",
                "bh_file", "hys_file", "fes_order", "n_steps",
                "max_iter", "tol", "relax", "newton", "solver",
            })
        elif self.method == METHOD_HDIV:
            # HDiv-VIM owns its solver and convergence contract
            # (radia.vim.Solve); the FEM iteration options do not apply.
            fields.update({
                "vol", "coil_script", "material", "sigma", "mu_r",
                "bh_file", "hys_file", "fes_order", "ima",
                "demag_backend",
            })
        elif self.method == METHOD_KELVIN_BENCH:
            fields.update({
                "vol", "fes_order", "kelvin_mu_r", "h0",
                "field_axis", "r_kelvin",
            })
        elif self.method == METHOD_CLEBSCH:
            fields.update({
                "clebsch_geometry", "clebsch_mu_r",
                "clebsch_maxh", "clebsch_fes_order",
            })
        return fields

    def missing_required_inputs(self) -> list[str]:
        if self.method in (METHOD_OMEGA, METHOD_APHI, METHOD_HDIV):
            return [] if self.coil_script.strip() else ["Coil script"]
        if self.method == METHOD_KELVIN_BENCH:
            return [] if self.vol.strip() else ["Mesh .vol"]
        return []

    def is_runnable(self) -> bool:
        return not self.missing_required_inputs()

    def build_command(self, *, python: str | None = None, panels_dir=None) -> list[str]:
        py = python or sys.executable
        if self.method in (METHOD_OMEGA, METHOD_APHI):
            return self._build_fem_command(py, panels_dir)
        if self.method == METHOD_HDIV:
            return self._build_hdiv_command(py, panels_dir)
        if self.method == METHOD_KELVIN_BENCH:
            return self._build_kelvin_command(py, panels_dir)
        if self.method == METHOD_CLEBSCH:
            return self._build_clebsch_command(py, panels_dir)
        raise ValueError(f"Unknown EM method: {self.method}")

    def _append_material(self, cmd: list[str]) -> None:
        cmd += ["--material", self.material]
        append_value(cmd, "--sigma", self.sigma)
        append_value(cmd, "--mu-r", self.mu_r)
        append_value(cmd, "--bh-file", self.bh_file)
        append_value(cmd, "--hys-file", self.hys_file)

    def _build_fem_command(self, py: str, panels_dir) -> list[str]:
        if self.solver not in FEM_SOLVERS:
            raise ValueError(f"Unsupported FE solver: {self.solver}")
        if not self.coil_script:
            raise ValueError("No coil script specified.")
        stem = self.vol or self.coil_script
        cmd = [
            py,
            calc_script("calc_accel_magnet.py", panels_dir),
            "--coil-script", self.coil_script,
            "--formulation", "a" if self.method == METHOD_APHI else "omega",
            "--fes-order", str(self.fes_order),
            "--n-steps", str(self.n_steps),
            "--max-iter", str(self.max_iter),
            "--tol", str(self.tol),
            "--relax", str(self.relax),
            "--solver", self.solver,
            "--msh-output", msh_output(stem, "_emfem"),
            "--output", json_output(stem, "_emfem"),
        ]
        append_value(cmd, "--vol", self.vol)
        self._append_material(cmd)
        append_switch(cmd, "--newton", self.newton)
        return cmd

    def _build_hdiv_command(self, py: str, panels_dir) -> list[str]:
        if not self.coil_script:
            raise ValueError("No coil script specified.")
        # HDiv-VIM owns its solver and convergence contract; an explicitly
        # changed FEM iteration option would otherwise be dropped silently.
        changed = [
            name for name in _FEM_ONLY_OPTIONS
            if getattr(self, name) != EMDesignSpec.__dataclass_fields__[name].default
        ]
        if changed:
            raise ValueError(
                "HDiv-VIM does not apply the FEM options "
                + ", ".join(changed)
                + "; its convergence contract is radia.vim.Solve's. Reset them "
                "or choose the Omega or A-Phi method.")
        stem = self.vol or self.coil_script
        cmd = [
            py,
            calc_script("calc_accel_hdiv.py", panels_dir),
            "--coil-script", self.coil_script,
            "--hdiv-order", str(self.fes_order),
            "--demag-backend", self.demag_backend,
            "--msh-output", msh_output(stem, "_hdiv"),
            "--output", json_output(stem, "_hdiv"),
        ]
        append_value(cmd, "--vol", self.vol)
        append_value(cmd, "--ima", self.ima)
        self._append_material(cmd)
        return cmd

    def _build_kelvin_command(self, py: str, panels_dir) -> list[str]:
        if not self.vol:
            raise ValueError("Kelvin Benchmark requires a .vol mesh.")
        return [
            py,
            calc_script("calc_kelvin_benchmark.py", panels_dir),
            "--vol", self.vol,
            "--fes-order", str(self.fes_order),
            "--mu-r", str(self.kelvin_mu_r),
            "--H0", str(self.h0),
            "--field-axis", self.field_axis,
            "--R-kelvin", str(self.r_kelvin),
            "--msh-output", msh_output(self.vol, "_kelvin_bench"),
            "--output", json_output(self.vol, "_kelvin_bench"),
        ]

    def _build_clebsch_command(self, py: str, panels_dir) -> list[str]:
        stem = self.clebsch_geometry or "clebsch"
        return [
            py,
            calc_script("calc_clebsch_hodograph.py", panels_dir),
            "--geometry", self.clebsch_geometry,
            "--mu-r", str(self.clebsch_mu_r),
            "--maxh", str(self.clebsch_maxh),
            "--fes-order", str(self.clebsch_fes_order),
            "--msh-output", msh_output(stem, "_clebsch"),
            "--output", json_output(stem, "_clebsch"),
        ]
