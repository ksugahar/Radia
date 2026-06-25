"""Small 2D magnetic-circuit / MMM-like motor quick checks.

The routines here are intentionally lightweight. They are not a replacement
for NGSolve AGE, JMAG, or external motor-deck solvers. They provide public-safe
first-order anchors that help an MCP client decide whether a deck or AGE solve is
physically plausible before spending solver time.
"""

from __future__ import annotations

from dataclasses import dataclass, asdict
import math
from typing import Any


MU0 = 4.0 * math.pi * 1e-7


@dataclass(frozen=True)
class MmmQuickInput:
    motor_type: str = "spm"
    pole_pairs: int = 4
    airgap_radius_m: float = 0.05
    stack_length_m: float = 0.05
    airgap_m: float = 1.0e-3
    turns_per_phase: float = 50.0
    phase_current_a: float = 10.0
    electrical_angle_deg: float = 0.0
    winding_factor: float = 0.933
    magnet_br_t: float = 1.2
    magnet_thickness_m: float = 3.0e-3
    magnet_arc_fraction: float = 0.75
    recoil_mu_r: float = 1.05
    carter_factor: float = 1.1
    saliency_ratio_lq_over_ld: float = 1.5
    base_inductance_h: float | None = None
    slip_hz: float = 5.0
    rotor_time_constant_s: float = 0.05
    rotor_resistance_ohm: float = 0.2


def _positive(value: float, fallback: float) -> float:
    return value if value > 0 else fallback


def _infer_family(motor_type: str) -> str:
    t = motor_type.strip().lower()
    if any(key in t for key in ("ipm", "interior")):
        return "ipm"
    if any(key in t for key in ("im", "induction", "cage")):
        return "induction"
    if any(key in t for key in ("srm", "sr motor", "switched")):
        return "srm"
    if any(key in t for key in ("synrm", "reluctance")):
        return "synrm"
    if "hysteresis" in t:
        return "hysteresis"
    return "spm"


def evaluate_mmm_quick_check(inp: MmmQuickInput) -> dict[str, Any]:
    """Evaluate a first-order 2D motor quick check."""
    family = _infer_family(inp.motor_type)
    p = max(1, int(inp.pole_pairs))
    radius = _positive(inp.airgap_radius_m, 0.05)
    stack = _positive(inp.stack_length_m, 0.05)
    gap = _positive(inp.airgap_m, 1.0e-3)
    turns = _positive(inp.turns_per_phase, 1.0)
    current = inp.phase_current_a
    theta = math.radians(inp.electrical_angle_deg)
    kw = _positive(inp.winding_factor, 0.933)
    carter_gap = _positive(inp.carter_factor, 1.0) * gap

    pole_pitch = math.pi * radius / p
    pole_area = max(1e-12, inp.magnet_arc_fraction * pole_pitch * stack)
    magnet_thickness = _positive(inp.magnet_thickness_m, 1.0e-4)
    br = inp.magnet_br_t
    bgap = br * magnet_thickness / (
        magnet_thickness + _positive(inp.recoil_mu_r, 1.0) * carter_gap
    )
    phi_pole_wb = bgap * pole_area
    lambda_pm_peak_wb = turns * kw * phi_pole_wb
    lambda_phase_wb = lambda_pm_peak_wb * math.cos(theta)
    back_emf_v_per_rad_s_mech = p * lambda_pm_peak_wb
    torque_constant_nm_per_a_peak = 1.5 * p * lambda_pm_peak_wb

    if inp.base_inductance_h and inp.base_inductance_h > 0:
        ld = inp.base_inductance_h
    else:
        ld = MU0 * turns * turns * pole_area / carter_gap
    lq = ld * _positive(inp.saliency_ratio_lq_over_ld, 1.0)
    iq = current * math.cos(theta)
    id_ = -current * math.sin(theta)
    pm_torque_nm = 1.5 * p * lambda_pm_peak_wb * iq
    reluctance_torque_nm = 1.5 * p * (ld - lq) * id_ * iq

    slip_rad_s = 2.0 * math.pi * abs(inp.slip_hz)
    x = slip_rad_s * _positive(inp.rotor_time_constant_s, 1e-9)
    rotor_loss_proxy_w = 3.0 * current * current * max(inp.rotor_resistance_ohm, 0.0) * (
        x * x / (1.0 + x * x)
    )
    induction_torque_proxy_nm = rotor_loss_proxy_w / max(slip_rad_s, 1e-9)

    if family == "induction":
        primary_quantity = "slip_loss_proxy"
        recommended_age = ("induction_machine", "airgap_eddy_machine", "deep_bar")
        applicability = "Use as a slip-frequency trend check only; AGE is required for field validation."
    elif family in ("srm", "synrm"):
        primary_quantity = "reluctance_torque_proxy"
        recommended_age = ("reluctance_torque", "saturating_inductance", "cross_saturation")
        applicability = "Use for saliency sign and current-angle trend checks; AGE is required for torque maps."
    elif family == "ipm":
        primary_quantity = "pm_plus_reluctance_torque_proxy"
        recommended_age = ("ld_lq", "mtpa", "field_weakening", "demag_margin")
        applicability = "Use for PM flux and saliency sanity checks; AGE is required for dq maps."
    elif family == "hysteresis":
        primary_quantity = "pm_flux_proxy"
        recommended_age = ("hysteresis_motor_loss", "hysteresis_play", "core_loss")
        applicability = "Use only for geometry flux sanity; hysteresis needs a stateful material model."
    else:
        primary_quantity = "pm_flux_linkage_proxy"
        recommended_age = ("back_emf", "cogging_torque", "ld_lq", "mtpa")
        applicability = "Use for PM flux-linkage and EMF constants; AGE is required for slotting/cogging."

    warnings = [
        "first-order 2D quick check, not a production solver",
        "no saturation iteration, no motion band, no end effects, no skew",
    ]
    if family == "induction":
        warnings.append("cage/end-ring coupling is represented only by a single-pole proxy")
    if family == "hysteresis":
        warnings.append("hysteresis loop area and vector history are not modeled")

    return {
        "schema_version": "radia-motor-mmm-quick/v1",
        "input": asdict(inp),
        "family": family,
        "primary_quantity": primary_quantity,
        "outputs": {
            "carter_gap_m": carter_gap,
            "pole_pitch_m": pole_pitch,
            "pole_area_m2": pole_area,
            "airgap_flux_density_t": bgap,
            "flux_per_pole_wb": phi_pole_wb,
            "lambda_pm_peak_wb": lambda_pm_peak_wb,
            "lambda_phase_wb": lambda_phase_wb,
            "back_emf_v_per_rad_s_mech": back_emf_v_per_rad_s_mech,
            "torque_constant_nm_per_a_peak": torque_constant_nm_per_a_peak,
            "ld_h": ld,
            "lq_h": lq,
            "id_a": id_,
            "iq_a": iq,
            "pm_torque_nm": pm_torque_nm,
            "reluctance_torque_nm": reluctance_torque_nm,
            "total_dq_torque_proxy_nm": pm_torque_nm + reluctance_torque_nm,
            "rotor_loss_proxy_w": rotor_loss_proxy_w,
            "induction_torque_proxy_nm": induction_torque_proxy_nm,
        },
        "recommended_age_targets": list(recommended_age),
        "applicability": applicability,
        "warnings": warnings,
        "public_boundary": (
            "This evaluator is an open, approximate magnetic-circuit quick check. "
            "It does not call commercial or external deck solvers, and it "
            "does not encode private benchmark numbers."
        ),
    }


def format_mmm_quick_check(result: dict[str, Any]) -> str:
    """Format a quick-check result as Markdown."""
    out = result["outputs"]
    lines = [
        "# 2D MMM/BEM-like motor quick check",
        "",
        f"- schema: `{result['schema_version']}`",
        f"- inferred family: `{result['family']}`",
        f"- primary quantity: `{result['primary_quantity']}`",
        f"- applicability: {result['applicability']}",
        f"- public boundary: {result['public_boundary']}",
        "",
        "## Key Estimates",
        f"- effective Carter gap: `{out['carter_gap_m']:.6g}` m",
        f"- air-gap flux density: `{out['airgap_flux_density_t']:.6g}` T",
        f"- flux per pole: `{out['flux_per_pole_wb']:.6g}` Wb",
        f"- PM flux-linkage peak: `{out['lambda_pm_peak_wb']:.6g}` Wb-turn",
        f"- phase flux-linkage at angle: `{out['lambda_phase_wb']:.6g}` Wb-turn",
        f"- back-EMF constant: `{out['back_emf_v_per_rad_s_mech']:.6g}` V/(rad/s mechanical)",
        f"- torque constant: `{out['torque_constant_nm_per_a_peak']:.6g}` N.m/A peak",
        f"- Ld/Lq proxy: `{out['ld_h']:.6g}` H / `{out['lq_h']:.6g}` H",
        f"- dq torque proxy: `{out['total_dq_torque_proxy_nm']:.6g}` N.m",
        f"- induction rotor-loss proxy: `{out['rotor_loss_proxy_w']:.6g}` W",
        "",
        "## Route To AGE Validation",
    ]
    lines.extend(f"- `ngsolve_usage(\"{target}\")`" for target in result["recommended_age_targets"])
    lines.extend(["", "## Warnings"])
    lines.extend(f"- {warning}" for warning in result["warnings"])
    return "\n".join(lines).rstrip()


def route_motor_validation(goal: str) -> dict[str, Any]:
    """Route a motor prompt to a public deck, MMM quick check, and AGE validation."""
    g = goal.lower()
    if any(term in g for term in ("induction", " cage", " im ", "slip", "deep bar")):
        family = "induction"
        deck_hint = "application/motor/emdlab_induction_bar_10"
        mmm = "motor_mmm_quick_check(motor_type='induction', slip_hz=...)"
        age = ("induction_machine", "airgap_eddy_machine", "deep_bar")
    elif any(term in g for term in ("srm", "switched reluctance", "sr motor")):
        family = "srm"
        deck_hint = "application/motor/emdlab_srm_pole_variants_10"
        mmm = "motor_mmm_quick_check(motor_type='srm', electrical_angle_deg=...)"
        age = ("reluctance_torque", "saturating_inductance")
    elif any(term in g for term in ("synrm", "reluctance motor")):
        family = "synrm"
        deck_hint = "application/motor/emdlab_synrm_flux_barrier_10"
        mmm = "motor_mmm_quick_check(motor_type='synrm', saliency_ratio_lq_over_ld=...)"
        age = ("synchronous_power_angle", "mtpa", "cross_saturation")
    elif any(term in g for term in ("ipm", "interior", "hairpin")):
        family = "ipm"
        deck_hint = "application/motor/emdlab_ipm_hairpin_10"
        mmm = "motor_mmm_quick_check(motor_type='ipm', electrical_angle_deg=...)"
        age = ("ld_lq", "mtpa", "field_weakening", "demag_margin")
    elif "hysteresis" in g:
        family = "hysteresis"
        deck_hint = "application/motor/hysteresis_motor_10"
        mmm = "motor_mmm_quick_check(motor_type='hysteresis')"
        age = ("hysteresis_motor_loss", "hysteresis_play")
    else:
        family = "spm"
        deck_hint = "application/motor/spm_surface_pm_10"
        mmm = "motor_mmm_quick_check(motor_type='spm', electrical_angle_deg=...)"
        age = ("back_emf", "cogging_torque", "ld_lq", "mtpa")

    return {
        "schema_version": "radia-motor-validation-router/v1",
        "goal": goal,
        "family": family,
        "deck_hint": deck_hint,
        "mmm_quick_check": mmm,
        "age_validation_targets": list(age),
        "workflow": [
            "Select and inspect a public motor input deck.",
            "Run motor_mmm_quick_check for a first-order sign/scale sanity check.",
            "Call motor_age_validation_plan(goal) to select the public AGE quality gates.",
            "Use NGSolve AGE / radia-ngsolve for the independent validation anchor.",
            "Only after the reduced quantities agree, move to local product runs.",
        ],
    }


def format_motor_validation_route(route: dict[str, Any]) -> str:
    """Format the hybrid motor validation route."""
    lines = [
        "# Motor validation route",
        "",
        f"- schema: `{route['schema_version']}`",
        f"- goal: {route['goal']}",
        f"- inferred family: `{route['family']}`",
        f"- public deck hint: `{route['deck_hint']}`",
        f"- MMM quick check: `{route['mmm_quick_check']}`",
        "",
        "## AGE Validation Targets",
    ]
    lines.extend(f"- `ngsolve_usage(\"{target}\")`" for target in route["age_validation_targets"])
    lines.extend(["", "## Workflow"])
    for i, step in enumerate(route["workflow"], 1):
        lines.append(f"{i}. {step}")
    return "\n".join(lines).rstrip()
