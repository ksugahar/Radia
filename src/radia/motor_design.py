"""UI-neutral motor panel design state."""

from __future__ import annotations

from dataclasses import dataclass
import sys

from .panel_design_common import append_value, calc_script, json_output


ANALYSIS_TRANSIENT = "Transient"
ANALYSIS_LAMINATION = "Lamination"
ANALYSIS_HDIV_REDUCED = "HDiv Reduced"
MOTOR_ANALYSES = (ANALYSIS_TRANSIENT, ANALYSIS_LAMINATION, ANALYSIS_HDIV_REDUCED)
TRANSIENT_METHODS = ("linearization", "coupled")
LINEAR_SOLVERS = ("pardiso", "sparsecholesky", "umfpack")
LAMINATION_MODES = ("cell", "global", "full")
LAMINATION_DRIVES = ("meanB", "current", "voltage")


@dataclass(slots=True)
class MotorDesignSpec:
    analysis: str = ANALYSIS_TRANSIENT
    vol: str = ""
    rotor_vol: str = ""
    method: str = "linearization"
    fes_order: int = 1
    linear_solver: str = "pardiso"
    nbr_phases: int = 3
    n_turns_per_slot: int = 100
    slot_area: str = "1e-4"
    stack_length: str = "0.05"
    n_pole_pairs: int = 4
    r_airgap_mid: str = "0.05"
    r_phase: str = "0.3872"
    l_endwinding: str = "0.0"
    pm_br_value: str = "0.0"
    j_inertia: str = "1e-3"
    b_viscous: str = "0.0"
    t_load: str = "0.0"
    v_amp: str = "0.0"
    v_freq: str = "0.0"
    theta_init: str = "0.0"
    omega_init: str = "0.0"
    t_end: str = "0.05"
    dt_fe: str = "1e-4"
    dt_circ: str = "1e-5"
    n_steps_per_fe: int = 10
    lamination_mode: str = "cell"
    d_iron: str = "0.35e-3"
    d_ins: str = "0.05e-3"
    sigma: str = "4e6"
    mu_r_iron: str = "5000.0"
    b_list: str = "0.5,1.0,1.5"
    freq_list: str = "50,500,5000"
    cell_n_elements: int = 100
    drive: str = "meanB"
    j_s_iron: str = "0.0"
    em_table: str = ""
    h_amplitude: str = "1000.0"
    freq: str = "1000.0"
    field_angle_deg: str = "0.0"
    rotor_angle_start_deg: str = "-45.0"
    rotor_angle_stop_deg: str = "45.0"
    rotor_angle_steps: int = 7
    energy_delta_deg: str = "0.25"
    circle_points: int = 1440
    center_x: str = "0.0"
    center_y: str = "0.0"
    hdiv_mu_r: str = "1000.0"
    hdiv_h_amplitude: str = "80000.0"
    hdiv_eta: str = "2.0"

    def visible_fields(self) -> set[str]:
        fields = {"analysis"}
        if self.analysis == ANALYSIS_TRANSIENT:
            fields.update({
                "vol", "method", "fes_order", "linear_solver",
                "nbr_phases", "n_turns_per_slot", "slot_area",
                "stack_length", "n_pole_pairs", "r_airgap_mid",
                "r_phase", "l_endwinding", "pm_br_value", "j_inertia",
                "b_viscous", "t_load", "v_amp", "v_freq", "theta_init",
                "omega_init", "t_end", "dt_fe", "dt_circ",
                "n_steps_per_fe",
            })
        elif self.analysis == ANALYSIS_LAMINATION:
            fields.update({
                "lamination_mode", "d_iron", "d_ins", "sigma",
                "mu_r_iron", "b_list", "freq_list", "cell_n_elements",
                "drive", "j_s_iron", "linear_solver",
            })
            if self.lamination_mode in ("global", "full"):
                fields.update({"vol", "h_amplitude", "freq", "fes_order"})
            if self.lamination_mode == "global":
                fields.add("em_table")
        elif self.analysis == ANALYSIS_HDIV_REDUCED:
            fields.update({
                "rotor_vol", "hdiv_mu_r", "hdiv_h_amplitude", "field_angle_deg",
                "rotor_angle_start_deg", "rotor_angle_stop_deg",
                "rotor_angle_steps", "r_airgap_mid", "stack_length",
                "energy_delta_deg", "circle_points", "center_x", "center_y",
                "hdiv_eta",
            })
        else:
            raise ValueError(f"Unknown motor analysis: {self.analysis}")
        return fields

    def missing_required_inputs(self) -> list[str]:
        if self.analysis == ANALYSIS_TRANSIENT:
            return [] if self.vol.strip() else ["Motor .vol"]
        if self.analysis == ANALYSIS_HDIV_REDUCED:
            return [] if self.rotor_vol.strip() else ["Rotor-only 2D .vol"]
        if self.analysis == ANALYSIS_LAMINATION:
            if self.lamination_mode in ("global", "full") and not self.vol.strip():
                return ["Motor .vol"]
            if self.lamination_mode == "global" and not self.em_table.strip():
                return ["EM table JSON"]
            return []
        raise ValueError(f"Unknown motor analysis: {self.analysis}")

    def is_runnable(self) -> bool:
        return not self.missing_required_inputs()

    def build_command(self, *, python: str | None = None, panels_dir=None) -> list[str]:
        if self.analysis == ANALYSIS_TRANSIENT:
            return self._build_transient(python or sys.executable, panels_dir)
        if self.analysis == ANALYSIS_LAMINATION:
            return self._build_lamination(python or sys.executable, panels_dir)
        if self.analysis == ANALYSIS_HDIV_REDUCED:
            return self._build_hdiv_reduced(python or sys.executable, panels_dir)
        raise ValueError(f"Unknown motor analysis: {self.analysis}")

    def _build_transient(self, py: str, panels_dir) -> list[str]:
        if not self.vol:
            raise ValueError("Transient motor analysis requires a .vol mesh.")
        return [
            py,
            calc_script("calc_motor_transient.py", panels_dir),
            "--vol", self.vol,
            "--method", self.method,
            "--fes-order", str(self.fes_order),
            "--linear-solver", self.linear_solver,
            "--nbr-phases", str(self.nbr_phases),
            "--n-turns-per-slot", str(self.n_turns_per_slot),
            "--slot-area", str(self.slot_area),
            "--stack-length", str(self.stack_length),
            "--n-pole-pairs", str(self.n_pole_pairs),
            "--r-airgap-mid", str(self.r_airgap_mid),
            "--R-phase", str(self.r_phase),
            "--L-endwinding", str(self.l_endwinding),
            "--pm-Br-value", str(self.pm_br_value),
            "--J-inertia", str(self.j_inertia),
            "--B-viscous", str(self.b_viscous),
            "--T-load", str(self.t_load),
            "--v-amp", str(self.v_amp),
            "--v-freq", str(self.v_freq),
            "--theta-init", str(self.theta_init),
            "--omega-init", str(self.omega_init),
            "--t-end", str(self.t_end),
            "--dt-FE", str(self.dt_fe),
            "--dt-circ", str(self.dt_circ),
            "--n-steps-per-FE", str(self.n_steps_per_fe),
            "--output", json_output(self.vol, "_motor_transient"),
        ]

    def _build_lamination(self, py: str, panels_dir) -> list[str]:
        if self.lamination_mode in ("global", "full") and not self.vol:
            raise ValueError("Global lamination analysis requires a .vol mesh.")
        if self.lamination_mode == "global" and not self.em_table:
            raise ValueError("Global lamination analysis requires --em-table.")
        stem = self.vol or "motor_lamination"
        cmd = [
            py,
            calc_script("calc_motor_lamination.py", panels_dir),
            "--mode", self.lamination_mode,
            "--d-iron", str(self.d_iron),
            "--d-ins", str(self.d_ins),
            "--sigma", str(self.sigma),
            "--mu-r-iron", str(self.mu_r_iron),
            "--B-list", self.b_list,
            "--freq-list", self.freq_list,
            "--cell-n-elements", str(self.cell_n_elements),
            "--drive", self.drive,
            "--J-s-iron", str(self.j_s_iron),
            "--linear-solver", self.linear_solver,
            "--output", json_output(stem, "_motor_lamination"),
        ]
        append_value(cmd, "--vol", self.vol)
        append_value(cmd, "--em-table", self.em_table)
        if self.lamination_mode in ("global", "full"):
            cmd += [
                "--H-amplitude", str(self.h_amplitude),
                "--freq", str(self.freq),
                "--fes-order", str(self.fes_order),
            ]
        return cmd

    def _build_hdiv_reduced(self, py: str, panels_dir) -> list[str]:
        if not self.rotor_vol:
            raise ValueError("HDiv Reduced motor analysis requires a rotor-only 2D .vol mesh.")
        return [
            py,
            calc_script("calc_motor_hdiv_reduced.py", panels_dir),
            "--vol", self.rotor_vol,
            "--mu-r", str(self.hdiv_mu_r),
            "--H-amplitude", str(self.hdiv_h_amplitude),
            "--field-angle-deg", str(self.field_angle_deg),
            "--rotor-angle-start-deg", str(self.rotor_angle_start_deg),
            "--rotor-angle-stop-deg", str(self.rotor_angle_stop_deg),
            "--rotor-angle-steps", str(self.rotor_angle_steps),
            "--maxwell-radius", str(self.r_airgap_mid),
            "--stack-length", str(self.stack_length),
            "--energy-delta-deg", str(self.energy_delta_deg),
            "--circle-points", str(self.circle_points),
            "--center-x", str(self.center_x),
            "--center-y", str(self.center_y),
            "--eta", str(self.hdiv_eta),
            "--output", json_output(self.rotor_vol, "_motor_hdiv_reduced"),
        ]
