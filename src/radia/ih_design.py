"""UI-neutral induction-heating design specification.

This module is the first step toward a browser/Jupyter based IH design
workbench.  It intentionally imports no Qt/PySide modules: desktop GUI
panels, notebooks, Voila apps, and tests should all be able to share the
same design state and command-building rules.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import sys


METHOD_PEEC_IND = "PEEC inductance (coil only, STEP)"
METHOD_BEMA_IND = "BEM-A inductance (coil only, .vol)"
METHOD_THERMAL_3D_STATIC = "Thermal: 3D static (no rotation)"
METHOD_THERMAL_3D_ROTATING = "Thermal: 3D + rotation (q_surf re-sampled per step)"
METHOD_THERMAL_AXISYM = "Thermal: 2D axisymmetric (rotation implicit)"
METHOD_PEEC_BEM = "PEEC + BEM weak coupling (workpiece)"
METHOD_BEMA_BEM = "BEM-A + BEM weak coupling (workpiece)"
METHOD_BEMA_BEM_STRONG = "BEM-A + BEM strong coupling (workpiece, CoupledBEMSolver)"
METHOD_PEEC_FEM_KELVIN = "PEEC coil + FEM wp (SIBC) + Kelvin"
METHOD_FEM_FULL = "Full simulation (FEM A-V + wp SIBC + Kelvin)"

IH_METHODS = (
    METHOD_PEEC_IND,
    METHOD_BEMA_IND,
    METHOD_PEEC_BEM,
    METHOD_BEMA_BEM,
    METHOD_BEMA_BEM_STRONG,
    METHOD_PEEC_FEM_KELVIN,
    METHOD_FEM_FULL,
    METHOD_THERMAL_3D_STATIC,
    METHOD_THERMAL_3D_ROTATING,
    METHOD_THERMAL_AXISYM,
)
THERMAL_METHODS = frozenset({
    METHOD_THERMAL_3D_STATIC,
    METHOD_THERMAL_3D_ROTATING,
    METHOD_THERMAL_AXISYM,
})
WORKPIECE_METHODS = frozenset({
    METHOD_PEEC_BEM,
    METHOD_BEMA_BEM,
    METHOD_BEMA_BEM_STRONG,
    METHOD_PEEC_FEM_KELVIN,
    METHOD_FEM_FULL,
    *THERMAL_METHODS,
})
PEEC_STEP_METHODS = frozenset({
    METHOD_PEEC_IND,
    METHOD_PEEC_BEM,
    METHOD_PEEC_FEM_KELVIN,
})
BEMA_COIL_VOL_METHODS = frozenset({
    METHOD_BEMA_IND,
    METHOD_BEMA_BEM,
    METHOD_BEMA_BEM_STRONG,
})

HEAT_SRC_UNIFORM = "Uniform q_surf [W/m^2]"
HEAT_SRC_SPATIAL = "Spatial q_surf .sol (from IH)"
MESH_TYPE_3D = "3D volume"
MESH_TYPE_AXISYM = "2D axisymmetric (r, z)"

THERMAL_PRESET_TO_CLI = {
    "Steel (Fe)": "steel",
    "Aluminum": "aluminum",
    "Copper": "copper",
    "Stainless 304": "stainless",
    "Brass": "brass",
    "Custom": "custom",
}
TIME_SCHEME_TO_CLI = {
    "Backward Euler": "backward-euler",
    "Crank-Nicolson": "crank-nicolson",
}

PEEC_SOLVERS = {
    "Dense LU (small)": {
        "coil_bem_solver": "dense-lu",
        "coil_saddle_solver": "auto",
        "wp_bem_backend": "intree-dense",
    },
    "HACApK (large)": {
        "coil_bem_solver": "hacapk-gmres",
        "coil_saddle_solver": "hacapk_cocr",
        "wp_bem_backend": "hacapk",
    },
}
FEM_SOLVERS = {
    "pardiso (direct)": "pardiso",
    "AMS (iterative, p=1)": "ams",
    "BDDC (iterative, p>=2)": "bddc",
    "iccg (fallback)": "iccg",
}


def _panels_dir() -> Path:
    return Path(__file__).resolve().parent / "panels"


def calc_script(name: str, panels_dir: str | Path | None = None) -> str:
    base = Path(panels_dir) if panels_dir is not None else _panels_dir()
    return str(base / name)


def msh_output(base_path: str | Path, suffix: str) -> str:
    base = Path(base_path) if base_path else Path("output")
    return str(base.with_suffix("")) + suffix + ".msh"


def json_output(base_path: str | Path, suffix: str) -> str:
    base = Path(base_path) if base_path else Path("output")
    return str(base.with_suffix("")) + suffix + ".json"


@dataclass(slots=True)
class IHDesignSpec:
    """All user-facing IH design inputs, independent of any UI toolkit."""

    method: str = METHOD_PEEC_IND
    frequency: str = "7000"
    current: str = "1.0"

    coil_material: str = "Copper"
    coil_sigma: str = "5.8e7"
    peec_step: str = ""
    peec_n_peri: int = 16
    peec_nwinc: int = 2
    peec_nhinc: int = 2

    coil_vol: str = ""
    coil_source_name: str = "source"
    coil_sink_name: str = "sink"

    wp_vol: str = ""
    wp_material: str = "Steel (mu_r=100)"
    wp_sigma: str = "5.0e6"
    mu_r: str = "100"
    half_thickness: str = "0.0125"
    # calc_inductance weak-path genus-1 loop DOF. "auto" (default)
    # applies it whenever it can: with the
    # "Dense LU (small)" solver preset (= intree-dense backend) + linear
    # SIBC + order 1, a genus-1 workpiece gets the solved shorted-turn
    # current automatically; inapplicable cases skip with a recorded
    # wp_loop_dof_skip_reason + caveat.  "on" makes unmet prerequisites
    # fail loud.  There is deliberately no "off" because the unextended
    # genus-1 solve omits the required cohomology mode, and the incident potential
    # is basis-determined in calc_inductance (P1 -> surface-Poisson),
    # not a panel knob.
    wp_loop_dof: str | bool = "auto"

    # calc_inductance weak-path genus-1 loop DOF (Takahashi-validated
    # 2026-07-17).  "auto" (default) applies it whenever it can: with the
    # "Dense LU (small)" solver preset (= intree-dense backend) + linear
    # SIBC + order 1, a genus-1 workpiece gets the solved shorted-turn
    # current automatically; inapplicable cases skip with a recorded
    # wp_loop_dof_skip_reason + caveat.  "on" makes unmet prerequisites
    # fail loud.  There is deliberately no "off" (the un-extended genus-1
    # solve is a known +25-30% over-estimate), and the incident potential
    # is basis-determined in calc_inductance (P1 -> surface-Poisson),
    # not a panel knob.
    wp_loop_dof: str = "auto"

    impedance_model: str = "Linear SIBC"
    bh_file: str = ""
    esim_max_iter: int = 30
    esim_tol: str = "1e-3"
    esim_per_panel: bool = False
    esim_anderson_m: int = 5
    esim_relax: str = "0.5"

    solver: str = "Dense LU (small)"
    fes_order: int = 1

    thermal_mesh_type: str = MESH_TYPE_3D
    heat_source: str = HEAT_SRC_UNIFORM
    q_uniform: str = "1.0e6"
    qsurf_sol: str = ""
    em_vol: str = ""
    qsurf_order: int = 1
    q_phi_average: bool = False
    n_phi_samples: int = 8
    surface_label: str = ""
    thermal_material: str = "Steel (Fe)"
    override_kcprho: bool = False
    rho: str = "7800"
    cp: str = "467"
    k: str = "46.6"
    rotation_rpm: float = 0.0
    rotation_axis: str = "z"
    h_conv: str = "10"
    t_ext: str = "20"
    emissivity: str = "0"
    t_init: str = "20"
    time_scheme: str = "Backward Euler"
    dt: str = "0.5"
    t_end: str = "5.0"
    linear_solver: str = "sparsecholesky"
    probe_point: str = ""
    csv_output: str = ""
    vtu_prefix: str = ""

    def __post_init__(self) -> None:
        # Preserve the pre-dropdown constructor contract without reviving the
        # removed unextended genus-1 route.
        if isinstance(self.wp_loop_dof, bool):
            self.wp_loop_dof = "on" if self.wp_loop_dof else "auto"
        if self.wp_loop_dof not in {"auto", "on"}:
            raise ValueError("wp_loop_dof must be auto or on")

    def impedance_model_cli(self) -> str:
        return "esim" if self.impedance_model.startswith("Nonlinear ESIM") else "sibc"

    def coil_solver_cli(self) -> str:
        if self.method in (METHOD_BEMA_IND, METHOD_BEMA_BEM,
                            METHOD_BEMA_BEM_STRONG):
            return "bem-a"
        return "peec"

    def visible_fields(self) -> set[str]:
        """Return spec fields that should be exposed for the current method."""

        fields = {"method"}
        thermal = self.method in THERMAL_METHODS
        if not thermal:
            fields.update({"frequency", "current", "solver"})
        if self.method in PEEC_STEP_METHODS:
            fields.update({
                "peec_step", "coil_material", "coil_sigma",
                "peec_n_peri",
            })
        if self.method == METHOD_PEEC_FEM_KELVIN:
            fields.update({"peec_nwinc", "peec_nhinc"})
        if self.method in BEMA_COIL_VOL_METHODS:
            fields.update({
                "coil_vol", "coil_source_name", "coil_sink_name",
                "coil_material", "coil_sigma",
            })
        if self.method in WORKPIECE_METHODS:
            fields.update({"wp_vol"})
        if self.method in (
            METHOD_PEEC_BEM, METHOD_BEMA_BEM,
            METHOD_PEEC_FEM_KELVIN, METHOD_FEM_FULL,
        ):
            fields.update({
                "wp_material", "wp_sigma", "mu_r", "half_thickness",
                "impedance_model", "fes_order",
            })
            if (
                self.method in (METHOD_PEEC_BEM, METHOD_BEMA_BEM)
                and self.impedance_model_cli() == "sibc"
                and self.fes_order == 1
            ):
                # calc_inductance weak path only (strong / FEM-Kelvin /
                # FEM-full do not take this flag; ESIM and P2 cannot apply it).
                fields.add("wp_loop_dof")
            if self.impedance_model_cli() == "esim":
                fields.update({
                    "bh_file", "esim_max_iter", "esim_per_panel",
                })
                if self.method != METHOD_PEEC_FEM_KELVIN:
                    fields.update({"esim_tol", "esim_relax"})
                if self.method in (METHOD_PEEC_BEM, METHOD_BEMA_BEM):
                    fields.add("esim_anderson_m")
        if self.method == METHOD_BEMA_BEM_STRONG:
            # Strong coupling (CoupledBEMSolver) is linear-SIBC only:
            # workpiece material + Leontovich Z_s knobs, no impedance-model /
            # fes-order / ESIM (the coupled solver fixes wp_order=1 and passes
            # a single global scalar Z_s).
            fields.update({
                "wp_material", "wp_sigma", "mu_r", "half_thickness",
            })
        if thermal:
            fields.update({
                "thermal_mesh_type", "heat_source", "surface_label",
                "thermal_material", "override_kcprho", "rho", "cp", "k",
                "h_conv", "t_ext", "emissivity", "t_init", "time_scheme",
                "dt", "t_end", "linear_solver", "fes_order",
                "probe_point", "csv_output", "vtu_prefix",
            })
            if self.method != METHOD_THERMAL_3D_STATIC:
                fields.update({"rotation_rpm", "rotation_axis"})
            if self.heat_source == HEAT_SRC_UNIFORM:
                fields.add("q_uniform")
            else:
                fields.update({"qsurf_sol", "em_vol", "qsurf_order"})
                if self.method == METHOD_THERMAL_AXISYM:
                    fields.add("n_phi_samples")
                else:
                    fields.add("q_phi_average")
        return fields

    def missing_required_inputs(
        self,
        *,
        check_exists: bool = False,
    ) -> list[str]:
        """Return human-readable missing input labels for the current method."""

        missing: list[str] = []

        def need(value: str | Path, label: str) -> None:
            text = str(value).strip()
            if not text:
                missing.append(label)
                return
            if check_exists and not Path(text).is_file():
                missing.append(f"{label} (not found)")

        if self.method in PEEC_STEP_METHODS:
            need(self.peec_step, "Coil STEP")
        if self.method in BEMA_COIL_VOL_METHODS:
            need(self.coil_vol, "Coil .vol")
        if self.method in WORKPIECE_METHODS:
            need(self.wp_vol, "Workpiece .vol")
        if self.impedance_model_cli() == "esim" and self.method not in THERMAL_METHODS:
            need(self.bh_file, "BH file")
        if self.method in THERMAL_METHODS and self.heat_source == HEAT_SRC_SPATIAL:
            need(self.qsurf_sol, "qsurf .sol")
            need(self.em_vol, "EM .vol")
        return missing

    def is_runnable(self, *, check_exists: bool = False) -> bool:
        """Return whether the current spec has the required user inputs."""

        return not self.missing_required_inputs(check_exists=check_exists)

    def build_command(
        self,
        vol_path: str | None = None,
        *,
        python: str | None = None,
        panels_dir: str | Path | None = None,
    ) -> list[str]:
        """Build the headless calc command for the current design.

        All IH methods, including thermal modes, are represented here so
        the browser/Jupyter workbench and the legacy PySide shell can share
        a single command contract.
        """

        py = python or sys.executable
        if self.method in THERMAL_METHODS:
            return self._build_thermal_command(py, panels_dir)
        if self.method in (METHOD_PEEC_IND, METHOD_BEMA_IND):
            return self._build_inductance_command(py, panels_dir)

        wp_vol = vol_path or self.wp_vol
        if not wp_vol:
            raise ValueError("No workpiece .vol file specified.")
        if self.method in (METHOD_PEEC_BEM, METHOD_BEMA_BEM):
            return self._build_weak_bem_command(wp_vol, py, panels_dir)
        if self.method == METHOD_BEMA_BEM_STRONG:
            return self._build_strong_bem_command(wp_vol, py, panels_dir)
        if self.method == METHOD_PEEC_FEM_KELVIN:
            return self._build_fem_kelvin_command(wp_vol, py, panels_dir)
        if self.method == METHOD_FEM_FULL:
            return self._build_fem_full_command(wp_vol, py, panels_dir)
        raise ValueError(f"Unknown IH method: {self.method}")

    def _coil_input_args(self) -> tuple[list[str], str]:
        coil_solver = self.coil_solver_cli()
        if coil_solver == "peec":
            if not self.peec_step:
                raise ValueError("Coil STEP file is required for PEEC mode.")
            coil_in = self.peec_step
            return ["--coil-step", coil_in], coil_in
        if not self.coil_vol:
            raise ValueError("Coil .vol file is required for BEM-A mode.")
        coil_in = self.coil_vol
        return [
            "--coil-vol", coil_in,
            "--coil-source-name", self.coil_source_name,
            "--coil-sink-name", self.coil_sink_name,
        ], coil_in

    def _bem_size(self) -> dict[str, str]:
        return PEEC_SOLVERS.get(
            self.solver,
            {"coil_bem_solver": "auto",
             "coil_saddle_solver": "auto",
             "wp_bem_backend": "hacapk"},
        )

    def _fem_solver(self) -> str:
        return FEM_SOLVERS.get(self.solver, "pardiso")

    def _build_inductance_command(
        self,
        py: str,
        panels_dir: str | Path | None,
    ) -> list[str]:
        coil_arg, coil_in = self._coil_input_args()
        coil_solver = self.coil_solver_cli()
        bem_size = self._bem_size()
        cmd = [
            py,
            calc_script("calc_inductance.py", panels_dir),
            *coil_arg,
            "--coil-solver", coil_solver,
            "--coil-bem-solver", bem_size["coil_bem_solver"],
            "--frequency", self.frequency,
            "--current", self.current,
            "--coil-sigma", self.coil_sigma,
            "--msh-output", msh_output(coil_in, "_peec_ind"),
            "--output", json_output(coil_in, "_peec_ind"),
        ]
        if coil_solver == "peec":
            cmd += ["--peec-n-peri", str(self.peec_n_peri)]
        return cmd

    def _build_weak_bem_command(
        self,
        vol_path: str,
        py: str,
        panels_dir: str | Path | None,
    ) -> list[str]:
        coil_arg, coil_in = self._coil_input_args()
        coil_solver = self.coil_solver_cli()
        bem_size = self._bem_size()
        cmd = [
            py,
            calc_script("calc_inductance.py", panels_dir),
            *coil_arg,
            "--coil-solver", coil_solver,
            "--coil-bem-solver", bem_size["coil_bem_solver"],
            "--wp-bem-backend", bem_size["wp_bem_backend"],
            "--frequency", self.frequency,
            "--current", self.current,
            "--coil-sigma", self.coil_sigma,
            "--vol", vol_path,
            "--wp-label", "sibc",
            "--sigma", self.wp_sigma,
            "--half-thickness", self.half_thickness,
            "--mu-r", self.mu_r,
            "--impedance-model", self.impedance_model_cli(),
            "--h1-order", str(self.fes_order),
            "--msh-output", msh_output(vol_path, "_peec_bem"),
            "--output", json_output(vol_path, "_peec_bem"),
        ]
        if coil_solver == "peec":
            cmd += ["--peec-n-peri", str(self.peec_n_peri)]
        if self.wp_loop_dof != "auto":
            cmd += ["--wp-loop-dof", self.wp_loop_dof]
        self._append_esim_args(cmd, include_anderson=True, kelvin=False)
        return cmd

    def _build_strong_bem_command(
        self,
        vol_path: str,
        py: str,
        panels_dir: str | Path | None,
    ) -> list[str]:
        """BEM-A coil + iterative CoupledBEMSolver strong coupling.

        Linear-SIBC only (a single global Leontovich Z_s); the coupled
        solver recomputes the coil surface current against the workpiece
        reaction field each Picard iteration, so ΔL includes the workpiece
        magnetic-energy term and P_wp is self-consistent.  Requires a BEM-A
        coil .vol (source/sink labels) -- no PEEC, no ESIM.
        """
        coil_arg, _coil_in = self._coil_input_args()
        bem_size = self._bem_size()
        return [
            py,
            calc_script("calc_inductance.py", panels_dir),
            *coil_arg,
            "--coil-solver", "bem-a",
            "--coil-bem-solver", bem_size["coil_bem_solver"],
            "--coil-saddle-solver", bem_size["coil_saddle_solver"],
            "--wp-bem-backend", bem_size["wp_bem_backend"],
            "--frequency", self.frequency,
            "--current", self.current,
            "--coil-sigma", self.coil_sigma,
            "--vol", vol_path,
            "--wp-label", "sibc",
            "--sigma", self.wp_sigma,
            "--half-thickness", self.half_thickness,
            "--mu-r", self.mu_r,
            "--coupling-mode", "strong",
            "--msh-output", msh_output(vol_path, "_bema_bem_strong"),
            "--output", json_output(vol_path, "_bema_bem_strong"),
        ]

    def _build_fem_kelvin_command(
        self,
        vol_path: str,
        py: str,
        panels_dir: str | Path | None,
    ) -> list[str]:
        if not self.peec_step:
            raise ValueError("PEEC STEP file is required for PEEC+FEM+Kelvin.")
        cmd = [
            py,
            calc_script("calc_fem_kelvin.py", panels_dir),
            "--vol", vol_path,
            "--fes-order", str(self.fes_order),
            "--frequency", self.frequency,
            "--material", "custom",
            "--sigma", self.wp_sigma,
            "--mu-r", self.mu_r,
            "--impedance", self.impedance_model_cli(),
            "--formulation", "total",
            "--current", self.current,
            "--half-thickness", self.half_thickness,
            "--solver", self._fem_solver(),
            "--peec-step", self.peec_step,
            "--peec-sigma", self.coil_sigma,
            "--peec-nwinc", str(self.peec_nwinc),
            "--peec-nhinc", str(self.peec_nhinc),
            "--peec-n-peri", str(self.peec_n_peri),
            "--require-kelvin",
            "--msh-output", msh_output(vol_path, "_fem_kelvin"),
            "--output", json_output(vol_path, "_fem_kelvin"),
        ]
        self._append_esim_args(cmd, include_anderson=False, kelvin=True)
        return cmd

    def _build_fem_full_command(
        self,
        vol_path: str,
        py: str,
        panels_dir: str | Path | None,
    ) -> list[str]:
        cmd = [
            py,
            calc_script("calc_fem_coilmesh.py", panels_dir),
            "--vol", vol_path,
            "--frequency", self.frequency,
            "--current", self.current,
            "--coil-sigma", self.coil_sigma,
            "--sigma", self.wp_sigma,
            "--mu-r", self.mu_r,
            "--half-thickness", self.half_thickness,
            "--fes-order", str(self.fes_order),
            "--solver", self._fem_solver(),
            "--sibc-bnd", "sibc",
            "--source-bnd", "source",
            "--sink-bnd", "sink",
            "--coil-mat", "coil",
            "--impedance-model", self.impedance_model_cli(),
            "--require-kelvin",
            "--msh-output", msh_output(vol_path, "_fem_full"),
            "--output", json_output(vol_path, "_fem_full"),
        ]
        self._append_esim_args(cmd, include_anderson=False, kelvin=False)
        return cmd

    def _append_esim_args(
        self,
        cmd: list[str],
        *,
        include_anderson: bool,
        kelvin: bool,
    ) -> None:
        if self.impedance_model_cli() != "esim":
            return
        if not self.bh_file:
            raise ValueError(
                "Nonlinear ESIM impedance model requires a BH file."
            )
        if kelvin:
            cmd += ["--bh-file", self.bh_file, "--max-iter", str(self.esim_max_iter)]
        else:
            cmd += [
                "--bh-file", self.bh_file,
                "--esim-max-iter", str(self.esim_max_iter),
                "--esim-tol", self.esim_tol,
                "--esim-relax", self.esim_relax,
            ]
            if include_anderson:
                cmd += ["--esim-anderson-m", str(self.esim_anderson_m)]
        if self.esim_per_panel:
            cmd.append("--esim-per-panel")

    def _thermal_mesh_type_for_method(self) -> str:
        if self.method == METHOD_THERMAL_AXISYM:
            return MESH_TYPE_AXISYM
        if self.method in (METHOD_THERMAL_3D_STATIC, METHOD_THERMAL_3D_ROTATING):
            return MESH_TYPE_3D
        return self.thermal_mesh_type

    def _build_thermal_command(
        self,
        py: str,
        panels_dir: str | Path | None,
    ) -> list[str]:
        wp = self.wp_vol
        if not wp:
            raise ValueError("Workpiece .vol is required for Heat.")

        mesh_type = self._thermal_mesh_type_for_method()
        is_axisym = mesh_type == MESH_TYPE_AXISYM
        calc = "calc_heat_axisym.py" if is_axisym else "calc_heat.py"
        material_cli = THERMAL_PRESET_TO_CLI.get(self.thermal_material, "custom")
        scheme_cli = TIME_SCHEME_TO_CLI.get(self.time_scheme, "backward-euler")
        rotation_rpm = 0.0 if self.method == METHOD_THERMAL_3D_STATIC else self.rotation_rpm

        cmd = [
            py,
            calc_script(calc, panels_dir),
            "--wp-vol", wp,
            "--surface-label", self.surface_label,
            "--material", material_cli,
            "--h-conv", self.h_conv,
            "--t-ext", self.t_ext,
            "--emissivity", self.emissivity,
            "--t-initial", self.t_init,
            "--dt", self.dt,
            "--t-end", self.t_end,
            "--time-scheme", scheme_cli,
            "--linear-solver", self.linear_solver,
            "--fes-order", str(self.fes_order),
            "--rotation-rpm", str(rotation_rpm),
            "--rotation-axis", str(self.rotation_axis),
            "--msh-output", msh_output(wp, "_heat"),
            "--output", json_output(wp, "_heat"),
        ]

        if material_cli == "custom" or self.override_kcprho:
            cmd += ["--rho", self.rho, "--cp", self.cp, "--k", self.k]

        if self.heat_source == HEAT_SRC_UNIFORM:
            cmd += ["--q-uniform", self.q_uniform]
        else:
            if not (self.qsurf_sol and self.em_vol):
                raise ValueError(
                    "Spatial qsurf mode requires a qsurf .sol and its companion EM .vol."
                )
            cmd += [
                "--qsurf-sol", self.qsurf_sol,
                "--em-vol", self.em_vol,
                "--qsurf-order", str(self.qsurf_order),
            ]
            if is_axisym:
                cmd += ["--n-phi-samples", str(self.n_phi_samples)]
            elif self.q_phi_average:
                cmd += ["--q-phi-average"]

        probe = self.probe_point.strip()
        if probe:
            cmd += ["--probe-point", probe]
        if self.csv_output:
            cmd += ["--csv-output", self.csv_output]
        if self.vtu_prefix:
            cmd += ["--vtu-prefix", self.vtu_prefix]
        return cmd
