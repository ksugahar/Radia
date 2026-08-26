"""UI-neutral stream-function coil design state."""

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


METHOD_DESIGN = "Design"
METHOD_PARETO = "Pareto"
METHOD_MANUFACTURE = "Manufacture"
METHOD_VOLUME_3D = "Volume 3D"
STREAMFUNCTION_METHODS = (
    METHOD_DESIGN,
    METHOD_PARETO,
    METHOD_MANUFACTURE,
    METHOD_VOLUME_3D,
)
_SURFACE_CLI = {
    METHOD_DESIGN: "design",
    METHOD_PARETO: "pareto",
    METHOD_MANUFACTURE: "manufacture",
}


@dataclass(slots=True)
class StreamFunctionDesignSpec:
    method: str = METHOD_DESIGN
    coil_vol: str = ""
    eval_vol: str = ""
    target_cf: str = "1"
    target_harmonic: str = ""
    harmonic_lmax: int = 3
    order: int = 3
    regularize: str = "h1"
    confine: str = "abe"
    alpha: str = "0.0"
    eval_max: int = 400
    pareto_lever: str = "alpha"
    alpha_min: str = "1e-4"
    alpha_max: str = "3e1"
    n_alpha: int = 12
    linf_iter: int = 8
    geom_scale_min: str = "0.6"
    geom_scale_max: str = "1.6"
    nlevels: int = 12
    optimize_levels: bool = False
    greedy_turns: str = ""
    greedy_dict: str = "contour"
    greedy_connector_weight: str = "0.0"
    greedy_target: str = ""
    greedy_plot: str = ""
    pin_tiling: bool = False
    pin_tiling_pins: int = 80
    pin_tiling_frac: str = "0.7"
    target_inductance: str = ""
    resonance_cap: str = ""
    nlevels_max: int = 60
    contour_sub: int = 1
    chain: str = "field_aware"
    chain_ncut: int = 0
    chain_passes: int = 0
    distort: bool = False
    distort_grid: int = 3
    distort_iter: int = 5
    step_output: str = ""
    peec: bool = False
    wire_diam: str = "1e-3"
    peec_freq: str = "1e5"
    flux_plot: str = ""
    flux_plane: str = "y"
    steps_plot: str = ""
    iron_vol: str = ""
    mu_r: str = "1000.0"
    iron_mat: str = "iron"
    iron_exact_source: bool = False
    iron_quad_order: str = ""
    shield_vol: str = ""
    shield_eval_vol: str = ""
    shield_weight: str = "1.0"
    volume_vol: str = ""
    target_bz: str = "1.0e-3"
    n_leaves: int = 3
    fes_order: int = 2
    aca_eps: str = "1.0e-10"
    rings_per_span: int = 14
    n_targets: int = 9
    nphi: int = 49
    nz: int = 41
    n_threads: int = 4

    def visible_fields(self) -> set[str]:
        fields = {"method"}
        if self.method == METHOD_VOLUME_3D:
            fields.update({
                "volume_vol", "target_bz", "n_leaves", "fes_order",
                "aca_eps", "rings_per_span", "n_targets", "nphi",
                "nz", "n_threads",
            })
            return fields
        fields.update({
            "coil_vol", "eval_vol", "target_cf", "target_harmonic",
            "harmonic_lmax", "order", "regularize", "confine",
            "alpha", "aca_eps", "eval_max", "iron_vol", "mu_r", "iron_mat",
            "iron_exact_source", "iron_quad_order", "shield_vol",
            "shield_eval_vol", "shield_weight",
        })
        if self.method == METHOD_PARETO:
            fields.update({
                "pareto_lever", "alpha_min", "alpha_max", "n_alpha",
                "linf_iter", "geom_scale_min", "geom_scale_max",
            })
        if self.method == METHOD_MANUFACTURE:
            fields.update({
                "nlevels", "optimize_levels", "greedy_turns",
                "greedy_dict", "greedy_connector_weight", "greedy_target",
                "greedy_plot", "pin_tiling", "pin_tiling_pins",
                "pin_tiling_frac", "target_inductance", "resonance_cap",
                "nlevels_max", "contour_sub", "chain", "chain_ncut",
                "chain_passes", "distort", "distort_grid", "distort_iter",
                "step_output", "peec", "wire_diam", "peec_freq",
                "flux_plot", "flux_plane", "steps_plot",
            })
        return fields

    def missing_required_inputs(self) -> list[str]:
        if self.method == METHOD_VOLUME_3D:
            return [] if self.volume_vol.strip() else ["Conductor .vol"]
        missing = []
        if not self.coil_vol.strip():
            missing.append("Coil surface .vol")
        if not self.eval_vol.strip():
            missing.append("Eval region .vol")
        return missing

    def is_runnable(self) -> bool:
        return not self.missing_required_inputs()

    def build_command(self, *, python: str | None = None, panels_dir=None) -> list[str]:
        if self.method == METHOD_VOLUME_3D:
            return self._build_volume(python or sys.executable, panels_dir)
        return self._build_surface(python or sys.executable, panels_dir)

    def _build_surface(self, py: str, panels_dir) -> list[str]:
        if not self.coil_vol:
            raise ValueError("Surface stream-function design requires --coil-vol.")
        if not self.eval_vol:
            raise ValueError("Surface stream-function design requires --eval-vol.")
        cli_method = _SURFACE_CLI[self.method]
        cmd = [
            py,
            calc_script("calc_streamfunction.py", panels_dir),
            "--coil-vol", self.coil_vol,
            "--eval-vol", self.eval_vol,
            "--method", cli_method,
            "--harmonic-lmax", str(self.harmonic_lmax),
            "--order", str(self.order),
            "--regularize", self.regularize,
            "--confine", self.confine,
            "--alpha", str(self.alpha),
            "--aca-eps", str(self.aca_eps),
            "--eval-max", str(self.eval_max),
            "--msh-output", msh_output(self.coil_vol, f"_sf_{cli_method}"),
            "--output", json_output(self.coil_vol, f"_sf_{cli_method}"),
        ]
        if self.target_harmonic:
            append_value(cmd, "--target-harmonic", self.target_harmonic)
        else:
            append_value(cmd, "--target-cf", self.target_cf)
        self._append_pareto(cmd)
        self._append_manufacture(cmd)
        self._append_material_aware(cmd)
        return cmd

    def _append_pareto(self, cmd: list[str]) -> None:
        if self.method != METHOD_PARETO:
            return
        cmd += [
            "--pareto-lever", self.pareto_lever,
            "--alpha-min", str(self.alpha_min),
            "--alpha-max", str(self.alpha_max),
            "--n-alpha", str(self.n_alpha),
            "--linf-iter", str(self.linf_iter),
            "--geom-scale-min", str(self.geom_scale_min),
            "--geom-scale-max", str(self.geom_scale_max),
        ]

    def _append_manufacture(self, cmd: list[str]) -> None:
        if self.method != METHOD_MANUFACTURE:
            return
        cmd += [
            "--nlevels", str(self.nlevels),
            "--greedy-dict", self.greedy_dict,
            "--greedy-connector-weight", str(self.greedy_connector_weight),
            "--pin-tiling-pins", str(self.pin_tiling_pins),
            "--pin-tiling-frac", str(self.pin_tiling_frac),
            "--nlevels-max", str(self.nlevels_max),
            "--contour-sub", str(self.contour_sub),
            "--chain", self.chain,
            "--chain-ncut", str(self.chain_ncut),
            "--chain-passes", str(self.chain_passes),
            "--distort-grid", str(self.distort_grid),
            "--distort-iter", str(self.distort_iter),
            "--wire-diam", str(self.wire_diam),
            "--peec-freq", str(self.peec_freq),
            "--flux-plane", self.flux_plane,
        ]
        append_switch(cmd, "--optimize-levels", self.optimize_levels)
        append_switch(cmd, "--pin-tiling", self.pin_tiling)
        append_switch(cmd, "--distort", self.distort)
        append_switch(cmd, "--peec", self.peec)
        append_value(cmd, "--greedy-turns", self.greedy_turns)
        append_value(cmd, "--greedy-target", self.greedy_target)
        append_value(cmd, "--greedy-plot", self.greedy_plot)
        append_value(cmd, "--target-inductance", self.target_inductance)
        append_value(cmd, "--resonance-cap", self.resonance_cap)
        append_value(cmd, "--step-output", self.step_output)
        append_value(cmd, "--flux-plot", self.flux_plot)
        append_value(cmd, "--steps-plot", self.steps_plot)

    def _append_material_aware(self, cmd: list[str]) -> None:
        append_value(cmd, "--iron-vol", self.iron_vol)
        if self.iron_vol:
            cmd += ["--mu-r", str(self.mu_r), "--iron-mat", self.iron_mat]
            append_switch(cmd, "--iron-exact-source", self.iron_exact_source)
            append_value(cmd, "--iron-quad-order", self.iron_quad_order)
        append_value(cmd, "--shield-vol", self.shield_vol)
        append_value(cmd, "--shield-eval-vol", self.shield_eval_vol)
        if self.shield_vol:
            cmd += ["--shield-weight", str(self.shield_weight)]

    def _build_volume(self, py: str, panels_dir) -> list[str]:
        if not self.volume_vol:
            raise ValueError("Volume 3D stream-function design requires --vol.")
        return [
            py,
            calc_script("calc_streamfunction_volume.py", panels_dir),
            "--vol", self.volume_vol,
            "--target-bz", str(self.target_bz),
            "--n-leaves", str(self.n_leaves),
            "--fes-order", str(self.fes_order),
            "--aca-eps", str(self.aca_eps),
            "--rings-per-span", str(self.rings_per_span),
            "--n-targets", str(self.n_targets),
            "--nphi", str(self.nphi),
            "--nz", str(self.nz),
            "--n-threads", str(self.n_threads),
            "--msh-output", msh_output(self.volume_vol, "_sf_volume"),
            "--output", json_output(self.volume_vol, "_sf_volume"),
        ]
