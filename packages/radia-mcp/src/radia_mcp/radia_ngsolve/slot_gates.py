# -*- coding: utf-8 -*-
"""Small solver-independent gates used by CAE loop slots.

These helpers intentionally avoid any commercial-tool provenance.  They encode
the physics checks that a loop slot can reuse before the source-tool/private
lane records where the reference data came from.
"""
from __future__ import annotations

import cmath
import math
from datetime import datetime, timezone

from .air_gap import (
    carter_coefficient,
    effective_air_gap,
    slotted_air_gap_permeance_factor,
)


MU0 = 4.0e-7 * math.pi
EPS0 = 8.8541878128e-12


def _parse_utc_like_datetime(value):
    if value in (None, ""):
        return None
    text = str(value).strip()
    if not text:
        return None
    try:
        parsed = datetime.fromisoformat(text.replace("Z", "+00:00"))
    except ValueError:
        return None
    if parsed.tzinfo is None:
        return parsed.replace(tzinfo=timezone.utc)
    return parsed.astimezone(timezone.utc)


def parallel_wire_force_per_length(i1, i2, separation_m, mu0=MU0):
    """Return the signed force per length between two long parallel wires.

    Positive means attraction for equal current directions.  The identity is
    ``F/L = mu0 * I1 * I2 / (2*pi*d)``.
    """

    d = float(separation_m)
    if d <= 0.0:
        raise ValueError("separation_m must be > 0")
    return float(mu0) * float(i1) * float(i2) / (2.0 * math.pi * d)


def coaxial_rc_duality_gate(
    inner_radius,
    outer_radius,
    eps_r=1.0,
    sigma=1.0,
    length=1.0,
    measured_capacitance=None,
    measured_resistance=None,
    rtol=1.0e-6,
    atol=0.0,
    eps0=EPS0,
):
    """Check the coaxial radial EC/ES duality.

    For a cylindrical annulus of length ``L`` with inner radius ``a`` and outer
    radius ``b``,

    ``C = 2*pi*eps0*eps_r*L/log(b/a)`` and
    ``R = log(b/a)/(2*pi*sigma*L)``.

    Their product removes the geometry: ``R*C = eps0*eps_r/sigma``.  The gate is
    solver-independent and is useful when a source tool reports either terminal
    charge capacitance, radial current resistance, or both.
    """

    a = float(inner_radius)
    b = float(outer_radius)
    er = float(eps_r)
    sig = float(sigma)
    ell = float(length)
    rel_tol = float(rtol)
    abs_tol = float(atol)
    e0 = float(eps0)
    if a <= 0.0:
        raise ValueError("inner_radius must be > 0")
    if b <= a:
        raise ValueError("outer_radius must be greater than inner_radius")
    if er <= 0.0:
        raise ValueError("eps_r must be > 0")
    if sig <= 0.0:
        raise ValueError("sigma must be > 0")
    if ell <= 0.0:
        raise ValueError("length must be > 0")
    if rel_tol < 0.0 or abs_tol < 0.0:
        raise ValueError("rtol and atol must be non-negative")
    if e0 <= 0.0:
        raise ValueError("eps0 must be > 0")

    log_ratio = math.log(b / a)
    capacitance = 2.0 * math.pi * e0 * er * ell / log_ratio
    resistance = log_ratio / (2.0 * math.pi * sig * ell)
    rc_product = resistance * capacitance
    rc_reference = e0 * er / sig
    rc_abs_error = abs(rc_product - rc_reference)
    rc_rel_error = rc_abs_error / max(abs(rc_product), abs(rc_reference), 1.0e-300)
    checks = {
        "capacitance_positive": capacitance > 0.0,
        "resistance_positive": resistance > 0.0,
        "geometry_free_rc_duality_ok": rc_abs_error <= abs_tol or rc_rel_error <= rel_tol,
    }

    cap_abs_error = None
    cap_rel_error = None
    if measured_capacitance is not None:
        cap_meas = float(measured_capacitance)
        cap_abs_error = abs(cap_meas - capacitance)
        cap_rel_error = cap_abs_error / max(abs(cap_meas), abs(capacitance), 1.0e-300)
        checks["measured_capacitance_ok"] = cap_abs_error <= abs_tol or cap_rel_error <= rel_tol

    resistance_abs_error = None
    resistance_rel_error = None
    if measured_resistance is not None:
        resistance_meas = float(measured_resistance)
        resistance_abs_error = abs(resistance_meas - resistance)
        resistance_rel_error = resistance_abs_error / max(abs(resistance_meas), abs(resistance), 1.0e-300)
        checks["measured_resistance_ok"] = resistance_abs_error <= abs_tol or resistance_rel_error <= rel_tol

    measured_rc_abs_error = None
    measured_rc_rel_error = None
    if measured_capacitance is not None and measured_resistance is not None:
        measured_rc = float(measured_capacitance) * float(measured_resistance)
        measured_rc_abs_error = abs(measured_rc - rc_reference)
        measured_rc_rel_error = measured_rc_abs_error / max(abs(measured_rc), abs(rc_reference), 1.0e-300)
        checks["measured_rc_duality_ok"] = (
            measured_rc_abs_error <= abs_tol or measured_rc_rel_error <= rel_tol
        )

    return {
        "policy": "coaxial_rc_duality_gate",
        "inner_radius_m": a,
        "outer_radius_m": b,
        "eps_r": er,
        "sigma_S_per_m": sig,
        "length_m": ell,
        "log_radius_ratio": log_ratio,
        "capacitance_F": capacitance,
        "resistance_ohm": resistance,
        "rc_product_s": rc_product,
        "rc_reference_s": rc_reference,
        "rc_abs_error_s": rc_abs_error,
        "rc_rel_error": rc_rel_error,
        "measured_capacitance_abs_error_F": cap_abs_error,
        "measured_capacitance_rel_error": cap_rel_error,
        "measured_resistance_abs_error_ohm": resistance_abs_error,
        "measured_resistance_rel_error": resistance_rel_error,
        "measured_rc_abs_error_s": measured_rc_abs_error,
        "measured_rc_rel_error": measured_rc_rel_error,
        "checks": checks,
        "status": "ok" if all(checks.values()) else "needs_attention",
        "rtol": rel_tol,
        "atol": abs_tol,
    }


def central_difference_periodic(values, spacing):
    """Return periodic central differences for equally spaced samples."""

    values = [float(value) for value in values]
    h = float(spacing)
    if len(values) < 3:
        raise ValueError("at least three samples are required")
    if h <= 0.0:
        raise ValueError("spacing must be > 0")
    n = len(values)
    return [
        (values[(i + 1) % n] - values[(i - 1) % n]) / (2.0 * h)
        for i in range(n)
    ]


def coenergy_torque_periodic_summary(
    theta_rad,
    coenergy_j,
    torque_nm,
    rtol=1.0e-6,
    atol=1.0e-9,
    near_zero_threshold=1.0e-12,
    near_zero_abs_tolerance_schema_id="coenergy_torque_near_zero_abs_tolerance_v1",
):
    """Check torque samples against ``T = dWprime/dtheta``.

    The samples must be equally spaced over one periodic cycle.  Near torque
    zero crossings the absolute tolerance is decisive; away from zero crossings
    the relative tolerance is reported as well.  The near-zero schema id is
    returned so source-tool table readers can store the same zero-crossing
    policy next to their private export metadata without publishing provenance.
    """

    theta = [float(value) for value in theta_rad]
    w = [float(value) for value in coenergy_j]
    torque = [float(value) for value in torque_nm]
    rel_tol = float(rtol)
    abs_tol = float(atol)
    near_zero_tol = float(near_zero_threshold)
    schema_id = str(near_zero_abs_tolerance_schema_id).strip()
    if not (len(theta) == len(w) == len(torque)):
        raise ValueError("theta_rad, coenergy_j, and torque_nm must have the same length")
    if len(theta) < 3:
        raise ValueError("at least three samples are required")
    if rel_tol < 0.0 or abs_tol < 0.0:
        raise ValueError("rtol and atol must be non-negative")
    if near_zero_tol < 0.0:
        raise ValueError("near_zero_threshold must be non-negative")
    if not schema_id:
        raise ValueError("near_zero_abs_tolerance_schema_id must not be empty")
    steps = [theta[i + 1] - theta[i] for i in range(len(theta) - 1)]
    h = sum(steps) / len(steps)
    if any(abs(step - h) > max(1.0e-12, 1.0e-9 * abs(h)) for step in steps):
        raise ValueError("theta samples must be equally spaced")
    estimated = central_difference_periodic(w, h)
    rows = []
    for angle, ref, est in zip(theta, torque, estimated):
        abs_error = abs(est - ref)
        rel_error = abs_error / max(abs(est), abs(ref), 1.0e-300)
        passed = abs_error <= abs_tol or rel_error <= rel_tol
        rows.append({
            "theta_rad": angle,
            "reference_torque_nm": ref,
            "estimated_torque_nm": est,
            "abs_error": abs_error,
            "rel_error": rel_error,
            "passed": passed,
        })
    near_zero_rows = [
        row for row in rows
        if abs(row["reference_torque_nm"]) <= near_zero_tol
    ]
    near_zero_rows_pass_absolute = all(row["abs_error"] <= abs_tol for row in near_zero_rows)
    checks = {
        "periodic_derivative_rows_pass": all(row["passed"] for row in rows),
        "near_zero_rows_use_absolute_tolerance": near_zero_rows_pass_absolute,
    }
    return {
        "policy": "coenergy_torque_periodic_derivative_gate",
        "n_samples": len(rows),
        "rtol": rel_tol,
        "atol": abs_tol,
        "near_zero_threshold_Nm": near_zero_tol,
        "near_zero_abs_tolerance_schema_id": schema_id,
        "near_zero_row_count": len(near_zero_rows),
        "near_zero_rows_pass_absolute_tolerance": near_zero_rows_pass_absolute,
        "max_abs_error": max(row["abs_error"] for row in rows),
        "max_rel_error": max(row["rel_error"] for row in rows),
        "n_passed": sum(1 for row in rows if row["passed"]),
        "checks": checks,
        "status": "ok" if all(checks.values()) else "needs_attention",
        "rows": rows,
    }


def carter_slot_opening_sweep_gate(slot_pitch, gap, slot_openings, tol=1.0e-12):
    """Summarize the Carter slot-opening correction before motor FEA.

    FEMM and radia-ngsolve both need the same preflight idea for slotted
    machines: slot openings reduce average air-gap permeance, so the physical
    gap ``g`` should be accompanied by ``g_eff = k_C g`` and
    ``P_slot/P_smooth = 1/k_C`` before interpreting magnetising inductance,
    no-load flux, or cogging/carter slotting studies.
    """

    tau = float(slot_pitch)
    g = float(gap)
    openings = [float(value) for value in slot_openings]
    tolerance = float(tol)
    if tau <= 0.0:
        raise ValueError("slot_pitch must be > 0")
    if g <= 0.0:
        raise ValueError("gap must be > 0")
    if not openings:
        raise ValueError("slot_openings must not be empty")
    if any(value < 0.0 for value in openings):
        raise ValueError("slot openings must be >= 0")
    if any(value >= tau for value in openings):
        raise ValueError("slot openings must be smaller than slot pitch")

    rows = []
    for opening in openings:
        kc = carter_coefficient(tau, g, opening)
        geff = effective_air_gap(tau, g, opening)
        permeance = slotted_air_gap_permeance_factor(tau, g, opening)
        rows.append({
            "slot_opening_m": opening,
            "opening_to_gap": opening / g,
            "carter_coefficient": kc,
            "effective_gap_m": geff,
            "permeance_factor": permeance,
            "inverse_identity_error": abs(kc * permeance - 1.0),
        })

    openings_sorted = all(rows[i]["slot_opening_m"] <= rows[i + 1]["slot_opening_m"] + tolerance for i in range(len(rows) - 1))
    kc_nondecreasing = all(rows[i]["carter_coefficient"] <= rows[i + 1]["carter_coefficient"] + tolerance for i in range(len(rows) - 1))
    permeance_nonincreasing = all(rows[i]["permeance_factor"] + tolerance >= rows[i + 1]["permeance_factor"] for i in range(len(rows) - 1))
    checks = {
        "zero_opening_identity": abs(carter_coefficient(tau, g, 0.0) - 1.0) <= tolerance,
        "all_kc_ge_one": all(row["carter_coefficient"] >= 1.0 - tolerance for row in rows),
        "all_permeance_le_one": all(row["permeance_factor"] <= 1.0 + tolerance for row in rows),
        "openings_sorted": openings_sorted,
        "kc_nondecreasing_with_opening": kc_nondecreasing,
        "permeance_nonincreasing_with_opening": permeance_nonincreasing,
        "inverse_identity": max(row["inverse_identity_error"] for row in rows) <= tolerance,
    }
    return {
        "policy": "carter_slot_opening_pre_fea_gate",
        "slot_pitch_m": tau,
        "gap_m": g,
        "tol": tolerance,
        "rows": rows,
        "max_carter_coefficient": max(row["carter_coefficient"] for row in rows),
        "min_permeance_factor": min(row["permeance_factor"] for row in rows),
        "checks": checks,
        "status": "ok" if all(checks.values()) else "needs_attention",
    }


def double_layer_winding_pitch_harmonic_gate(
    q,
    coil_pitch,
    pole_pitch,
    phases=3,
    harmonics=(1, 3, 5, 7),
    expected_kw_signs=None,
    expected_kp_signs=None,
    tol=1.0e-12,
):
    """Check double-layer distributed-winding harmonic factors.

    ``q`` is slots per pole per phase, ``coil_pitch/pole_pitch`` is the
    short-pitch fraction, and ``phases`` defaults to three.  The gate keeps the
    distribution factor, pitch factor, and their product visible before a FEMM
    or radia-ngsolve rotor sweep consumes winding data.
    """

    q_value = float(q)
    y = float(coil_pitch)
    tau = float(pole_pitch)
    m = int(phases)
    tolerance = float(tol)
    if q_value <= 0.0:
        raise ValueError("q must be > 0")
    if y <= 0.0:
        raise ValueError("coil_pitch must be > 0")
    if tau <= 0.0:
        raise ValueError("pole_pitch must be > 0")
    if y > tau:
        raise ValueError("coil_pitch must not exceed pole_pitch")
    if m <= 0:
        raise ValueError("phases must be positive")
    if tolerance < 0.0:
        raise ValueError("tol must be non-negative")
    harmonic_list = [int(value) for value in harmonics]
    if not harmonic_list:
        raise ValueError("harmonics must not be empty")
    if any(value <= 0 for value in harmonic_list):
        raise ValueError("harmonic numbers must be positive")

    gamma = math.pi / (m * q_value)
    beta = y / tau

    def sign(value):
        if abs(value) <= tolerance:
            return 0
        return 1 if value > 0.0 else -1

    rows = []
    for harmonic in sorted(set(harmonic_list)):
        denom = q_value * math.sin(harmonic * gamma / 2.0)
        numerator = math.sin(harmonic * q_value * gamma / 2.0)
        if abs(denom) <= tolerance:
            kd = 1.0 if abs(numerator) <= tolerance else math.copysign(float("inf"), numerator)
        else:
            kd = numerator / denom
        kp = math.sin(harmonic * beta * math.pi / 2.0)
        kw = kd * kp
        rows.append({
            "harmonic": harmonic,
            "distribution_factor_kd": kd,
            "pitch_factor_kp": kp,
            "winding_factor_kw": kw,
            "kd_sign": sign(kd),
            "kp_sign": sign(kp),
            "kw_sign": sign(kw),
        })

    row_by_harmonic = {row["harmonic"]: row for row in rows}
    checks = {
        "has_fundamental": 1 in row_by_harmonic,
        "fundamental_kw_positive": row_by_harmonic.get(1, {}).get("kw_sign") == 1,
    }
    expected_kw = {int(k): int(v) for k, v in (expected_kw_signs or {}).items()}
    for harmonic, expected in expected_kw.items():
        checks[f"kw_sign_h{harmonic}_ok"] = (
            harmonic in row_by_harmonic and row_by_harmonic[harmonic]["kw_sign"] == expected
        )
    expected_kp = {int(k): int(v) for k, v in (expected_kp_signs or {}).items()}
    for harmonic, expected in expected_kp.items():
        checks[f"kp_sign_h{harmonic}_ok"] = (
            harmonic in row_by_harmonic and row_by_harmonic[harmonic]["kp_sign"] == expected
        )

    return {
        "policy": "double_layer_winding_pitch_harmonic_gate",
        "q_slots_per_pole_per_phase": q_value,
        "coil_pitch_slots": y,
        "pole_pitch_slots": tau,
        "short_pitch_fraction": beta,
        "phases": m,
        "slot_angle_electrical_deg": math.degrees(gamma),
        "rows": rows,
        "checks": checks,
        "status": "ok" if all(checks.values()) else "needs_attention",
        "tol": tolerance,
    }


def acoustic_interface_result_package_gate(
    artifacts,
    expected_case_id=None,
    expected_run_id=None,
    expected_export_id=None,
    expected_frequency_hz=None,
    required_kinds=("material_impedance", "interface_continuity", "power_split"),
    frequency_rtol=1.0e-12,
    residual_tol=1.0e-9,
):
    """Check a reusable acoustic FEM/BEM interface result package.

    The algebraic interface identities are handled by
    ``acoustic_normal_incidence_interface_gate``. This manifest-level gate keeps
    material impedance, interface continuity, and passive power split rows tied
    to one case/run/export/frequency identity before reuse.
    """

    rows_in = list(artifacts)
    if not rows_in:
        raise ValueError("artifacts must not be empty")

    def _norm(value):
        return str(value or "").strip().lower().replace("-", "_").replace(" ", "_")

    def _first(row, names):
        for name in names:
            if name in row and row[name] is not None:
                return row[name]
        return None

    def _float_or_none(value):
        if value is None:
            return None
        return float(value)

    def _consistent(values, rtol):
        if len(values) <= 1:
            return True
        scale = max([1.0] + [abs(value) for value in values])
        return (max(values) - min(values)) <= rtol * scale

    required = tuple(_norm(kind) for kind in required_kinds)
    if not required:
        raise ValueError("required_kinds must not be empty")
    frequency_tolerance = float(frequency_rtol)
    residual_tolerance = float(residual_tol)
    if frequency_tolerance < 0.0:
        raise ValueError("frequency_rtol must be non-negative")
    if residual_tolerance < 0.0:
        raise ValueError("residual_tol must be non-negative")

    expected_policies = {
        "material_impedance": {
            "acoustic_material_impedance_contract",
            "acoustic_normal_incidence_interface_gate",
        },
        "interface_continuity": {
            "acoustic_interface_continuity_contract",
            "acoustic_normal_incidence_interface_gate",
        },
        "power_split": {
            "acoustic_power_split_contract",
            "acoustic_normal_incidence_interface_gate",
        },
    }

    details = []
    kind_counts = {}
    case_ids = []
    run_ids = []
    export_ids = []
    frequencies = []
    missing_case_id = []
    missing_run_id = []
    missing_export_id = []
    missing_frequency = []
    missing_source_tool = []
    missing_paths = []
    unknown_kinds = []
    bad_upstream_status = []
    bad_upstream_policy = []
    bad_material_impedance = []
    bad_pressure_residual = []
    bad_velocity_residual = []
    bad_power_split = []
    pressure_residuals = []
    velocity_residuals = []
    power_balance_errors = []
    reflected_ratios = []
    transmitted_ratios = []

    for index, row in enumerate(rows_in, start=1):
        if not isinstance(row, dict):
            raise ValueError("each artifact must be a dictionary")
        kind = _norm(_first(row, ("kind", "artifact_kind", "type")))
        case_id = _first(row, ("case_id", "model_id", "problem_id"))
        run_id = _first(row, ("run_id", "solver_run_id", "simulation_id"))
        export_id = _first(row, ("export_id", "result_id", "dataset_id"))
        source_tool = _first(row, ("source_tool", "tool", "source"))
        path = _first(row, ("path", "file", "artifact_path", "table_path"))
        gate_policy = _first(row, ("gate_policy", "policy", "validator"))
        gate_policy_norm = _norm(gate_policy)
        status = _first(row, ("status", "gate_status", "validation_status"))
        status_norm = _norm(status)
        pass_flag = bool(row.get("pass", False))
        frequency = _float_or_none(_first(row, ("frequency_hz", "design_frequency_hz", "freq_hz")))

        if kind:
            kind_counts[kind] = kind_counts.get(kind, 0) + 1
        else:
            unknown_kinds.append(index)
        if kind and kind not in set(required) | set(expected_policies):
            unknown_kinds.append(kind)
        if not case_id:
            missing_case_id.append(index)
        else:
            case_ids.append(str(case_id))
        if not run_id:
            missing_run_id.append(index)
        else:
            run_ids.append(str(run_id))
        if not export_id:
            missing_export_id.append(index)
        else:
            export_ids.append(str(export_id))
        if not source_tool:
            missing_source_tool.append(index)
        if not path:
            missing_paths.append(index)
        if frequency is None:
            missing_frequency.append(index)
        else:
            frequencies.append(frequency)
        if not (pass_flag or status_norm in {"ok", "pass", "passed", "verified"}):
            bad_upstream_status.append({"index": index, "kind": kind, "status": status})
        if kind in expected_policies and gate_policy_norm not in expected_policies[kind]:
            bad_upstream_policy.append({
                "index": index,
                "kind": kind,
                "gate_policy": gate_policy,
                "expected": sorted(expected_policies[kind]),
            })

        if kind == "material_impedance":
            z_left = _float_or_none(_first(row, ("z_left_Pa_s_per_m", "z_left", "impedance_left")))
            z_right = _float_or_none(_first(row, ("z_right_Pa_s_per_m", "z_right", "impedance_right")))
            rho_left = _float_or_none(_first(row, ("rho_left_kg_per_m3", "rho_left")))
            c_left = _float_or_none(_first(row, ("c_left_m_per_s", "c_left")))
            rho_right = _float_or_none(_first(row, ("rho_right_kg_per_m3", "rho_right")))
            c_right = _float_or_none(_first(row, ("c_right_m_per_s", "c_right")))
            if z_left is None and rho_left is not None and c_left is not None:
                z_left = rho_left * c_left
            if z_right is None and rho_right is not None and c_right is not None:
                z_right = rho_right * c_right
            if z_left is None or z_left <= 0.0 or z_right is None or z_right <= 0.0:
                bad_material_impedance.append({"index": index, "z_left": z_left, "z_right": z_right})

        if kind == "interface_continuity":
            pressure = abs(_float_or_none(_first(row, (
                "pressure_jump_Pa",
                "pressure_continuity_residual_Pa",
                "pressure_residual_Pa",
            ))) or 0.0)
            velocity = abs(_float_or_none(_first(row, (
                "velocity_jump_m_per_s",
                "normal_velocity_continuity_residual_m_per_s",
                "velocity_residual_m_per_s",
            ))) or 0.0)
            pressure_residuals.append(pressure)
            velocity_residuals.append(velocity)
            if pressure > residual_tolerance:
                bad_pressure_residual.append({"index": index, "pressure_residual_Pa": pressure})
            if velocity > residual_tolerance:
                bad_velocity_residual.append({"index": index, "velocity_residual_m_per_s": velocity})

        if kind == "power_split":
            reflected = _float_or_none(_first(row, ("reflected_power_ratio", "reflection_power_ratio")))
            transmitted = _float_or_none(_first(row, ("transmitted_power_ratio", "transmission_power_ratio")))
            balance = _float_or_none(_first(row, ("power_balance_error", "power_split_residual")))
            if reflected is None or transmitted is None:
                bad_power_split.append({"index": index, "reason": "missing_power_ratios"})
            else:
                reflected_ratios.append(reflected)
                transmitted_ratios.append(transmitted)
                if balance is None:
                    balance = reflected + transmitted - 1.0
                power_balance_errors.append(balance)
                if reflected < -residual_tolerance or transmitted < -residual_tolerance:
                    bad_power_split.append({"index": index, "reason": "negative_power_ratio"})
                if abs(balance) > residual_tolerance:
                    bad_power_split.append({"index": index, "reason": "power_balance", "power_balance_error": balance})

        details.append({
            "index": index,
            "kind": kind,
            "case_id": None if case_id is None else str(case_id),
            "run_id": None if run_id is None else str(run_id),
            "export_id": None if export_id is None else str(export_id),
            "frequency_hz": frequency,
            "source_tool": source_tool,
            "path": path,
            "gate_policy": gate_policy,
            "status": status,
        })

    unique_case_ids = sorted(set(case_ids))
    unique_run_ids = sorted(set(run_ids))
    unique_export_ids = sorted(set(export_ids))
    expected_frequency_ok = True
    if expected_frequency_hz is not None:
        expected_frequency_ok = all(
            abs(value - float(expected_frequency_hz))
            <= frequency_tolerance * max(1.0, abs(float(expected_frequency_hz)))
            for value in frequencies
        )
    checks = {
        "required_kinds_present": all(kind_counts.get(kind, 0) > 0 for kind in required),
        "single_case_id": len(unique_case_ids) == 1 and not missing_case_id,
        "single_run_id": len(unique_run_ids) == 1 and not missing_run_id,
        "single_export_id": len(unique_export_ids) == 1 and not missing_export_id,
        "expected_case_id_ok": expected_case_id is None or unique_case_ids == [str(expected_case_id)],
        "expected_run_id_ok": expected_run_id is None or unique_run_ids == [str(expected_run_id)],
        "expected_export_id_ok": expected_export_id is None or unique_export_ids == [str(expected_export_id)],
        "frequencies_recorded": not missing_frequency and bool(frequencies),
        "frequencies_consistent": _consistent(frequencies, frequency_tolerance),
        "expected_frequency_ok": expected_frequency_ok,
        "source_tools_recorded": not missing_source_tool,
        "paths_recorded": not missing_paths,
        "all_rows_verified": not bad_upstream_status,
        "upstream_policies_match": not bad_upstream_policy,
        "material_impedances_positive": not bad_material_impedance,
        "pressure_continuity_residual_ok": not bad_pressure_residual,
        "normal_velocity_continuity_residual_ok": not bad_velocity_residual,
        "passive_power_split_ok": not bad_power_split,
        "no_unknown_kinds": not unknown_kinds,
    }
    issues = [name for name, ok in checks.items() if not ok]

    return {
        "policy": "acoustic_interface_result_package_gate",
        "status": "ok" if not issues else "needs_attention",
        "required_kinds": list(required),
        "present_kinds": dict(sorted(kind_counts.items())),
        "case_ids": unique_case_ids,
        "run_ids": unique_run_ids,
        "export_ids": unique_export_ids,
        "expected_case_id": None if expected_case_id is None else str(expected_case_id),
        "expected_run_id": None if expected_run_id is None else str(expected_run_id),
        "expected_export_id": None if expected_export_id is None else str(expected_export_id),
        "expected_frequency_hz": None if expected_frequency_hz is None else float(expected_frequency_hz),
        "frequency_hz_values": frequencies,
        "max_pressure_residual_Pa": max(pressure_residuals) if pressure_residuals else None,
        "max_normal_velocity_residual_m_per_s": max(velocity_residuals) if velocity_residuals else None,
        "max_abs_power_balance_error": max((abs(value) for value in power_balance_errors), default=None),
        "min_reflected_power_ratio": min(reflected_ratios) if reflected_ratios else None,
        "min_transmitted_power_ratio": min(transmitted_ratios) if transmitted_ratios else None,
        "missing_case_id_rows": missing_case_id,
        "missing_run_id_rows": missing_run_id,
        "missing_export_id_rows": missing_export_id,
        "missing_frequency_rows": missing_frequency,
        "missing_source_tool_rows": missing_source_tool,
        "missing_path_rows": missing_paths,
        "unknown_kinds": unknown_kinds,
        "bad_upstream_status_rows": bad_upstream_status,
        "bad_upstream_policy_rows": bad_upstream_policy,
        "bad_material_impedance_rows": bad_material_impedance,
        "bad_pressure_residual_rows": bad_pressure_residual,
        "bad_velocity_residual_rows": bad_velocity_residual,
        "bad_power_split_rows": bad_power_split,
        "artifacts": details,
        "checks": checks,
        "issues": issues,
        "frequency_rtol": frequency_tolerance,
        "residual_tol": residual_tolerance,
        "version_note": (
            "Run after material impedance, interface continuity, and power "
            "split sub-gates so FEM/BEM acoustic data cannot mix case identity, "
            "frequency, or passive power accounting."
        ),
    }


def acoustic_impedance_power_result_package_gate(
    artifacts,
    expected_case_id=None,
    expected_run_id=None,
    expected_export_id=None,
    expected_frequency_hz=None,
    required_kinds=("impedance_power",),
    frequency_rtol=1.0e-12,
    residual_tol=1.0e-9,
):
    """Check an acoustic one-port impedance/power result package.

    This gate is the one-port companion to the two-medium acoustic interface
    package. It keeps a radiation/impedance boundary row tied to one
    case/run/export/frequency identity and verifies the passive power identity

    ``A = 1-|Gamma|^2`` and
    ``I_abs = I_inc - I_ref = Re(0.5*p*conj(v))``.

    It is intentionally solver-independent; source-tool provenance belongs in
    the private lane.
    """

    rows_in = list(artifacts)
    if not rows_in:
        raise ValueError("artifacts must not be empty")

    def _norm(value):
        return str(value or "").strip().lower().replace("-", "_").replace(" ", "_")

    def _first(row, names):
        for name in names:
            if name in row and row[name] is not None:
                return row[name]
        return None

    def _float_or_none(value):
        if value is None:
            return None
        return float(value)

    def _complex_or_none(row, base_names):
        for base in base_names:
            value = row.get(base)
            if value is not None:
                if isinstance(value, dict):
                    return complex(float(value.get("real", 0.0)), float(value.get("imag", 0.0)))
                return complex(value)
            real = row.get(f"{base}_real")
            imag = row.get(f"{base}_imag")
            if real is not None or imag is not None:
                return complex(float(real or 0.0), float(imag or 0.0))
        return None

    def _consistent(values, rtol):
        if len(values) <= 1:
            return True
        scale = max([1.0] + [abs(value) for value in values])
        return (max(values) - min(values)) <= rtol * scale

    required = tuple(_norm(kind) for kind in required_kinds)
    if not required:
        raise ValueError("required_kinds must not be empty")
    frequency_tolerance = float(frequency_rtol)
    residual_tolerance = float(residual_tol)
    if frequency_tolerance < 0.0:
        raise ValueError("frequency_rtol must be non-negative")
    if residual_tolerance < 0.0:
        raise ValueError("residual_tol must be non-negative")

    expected_policies = {
        "impedance_power": {
            "acoustic_impedance_power_contract",
            "velocity_positive_into_load_pressure_reflection_gamma_zload_minus_zn_over_zload_plus_zn",
        },
        "radiation_impedance": {
            "acoustic_impedance_power_contract",
            "velocity_positive_into_load_pressure_reflection_gamma_zload_minus_zn_over_zload_plus_zn",
        },
        "boundary_power": {
            "acoustic_boundary_power_contract",
            "outward_active_power_positive_for_pressure_times_conjugate_normal_velocity",
        },
    }

    details = []
    kind_counts = {}
    case_ids = []
    run_ids = []
    export_ids = []
    frequencies = []
    missing_case_id = []
    missing_run_id = []
    missing_export_id = []
    missing_frequency = []
    missing_source_tool = []
    missing_paths = []
    unknown_kinds = []
    bad_upstream_status = []
    bad_upstream_policy = []
    bad_impedance = []
    bad_absorption = []
    bad_power_balance = []
    bad_active_power = []
    absorption_values = []
    power_balance_errors = []
    active_powers = []
    impedance_real_values = []

    for index, row in enumerate(rows_in, start=1):
        if not isinstance(row, dict):
            raise ValueError("each artifact must be a dictionary")
        kind = _norm(_first(row, ("kind", "artifact_kind", "type")))
        case_id = _first(row, ("case_id", "model_id", "problem_id"))
        run_id = _first(row, ("run_id", "solver_run_id", "simulation_id"))
        export_id = _first(row, ("export_id", "result_id", "dataset_id"))
        source_tool = _first(row, ("source_tool", "tool", "source"))
        path = _first(row, ("path", "file", "artifact_path", "table_path"))
        gate_policy = _first(row, ("gate_policy", "policy", "validator"))
        gate_policy_norm = _norm(gate_policy)
        status = _first(row, ("status", "gate_status", "validation_status"))
        status_norm = _norm(status)
        pass_flag = bool(row.get("pass", False))
        frequency = _float_or_none(_first(row, ("frequency_hz", "design_frequency_hz", "freq_hz")))

        if kind:
            kind_counts[kind] = kind_counts.get(kind, 0) + 1
        else:
            unknown_kinds.append(index)
        if kind and kind not in set(required) | set(expected_policies):
            unknown_kinds.append(kind)
        if not case_id:
            missing_case_id.append(index)
        else:
            case_ids.append(str(case_id))
        if not run_id:
            missing_run_id.append(index)
        else:
            run_ids.append(str(run_id))
        if not export_id:
            missing_export_id.append(index)
        else:
            export_ids.append(str(export_id))
        if not source_tool:
            missing_source_tool.append(index)
        if not path:
            missing_paths.append(index)
        if frequency is None:
            missing_frequency.append(index)
        else:
            frequencies.append(frequency)
        if not (pass_flag or status_norm in {"ok", "pass", "passed", "verified"}):
            bad_upstream_status.append({"index": index, "kind": kind, "status": status})
        if kind in expected_policies and gate_policy_norm not in expected_policies[kind]:
            bad_upstream_policy.append({
                "index": index,
                "kind": kind,
                "gate_policy": gate_policy,
                "expected": sorted(expected_policies[kind]),
            })

        if kind in {"impedance_power", "radiation_impedance", "boundary_power"}:
            z_load = _complex_or_none(row, ("specific_impedance", "radiation_impedance", "z_load"))
            gamma = _complex_or_none(row, ("pressure_reflection_coefficient", "reflection_coefficient", "gamma"))
            z_normal = _float_or_none(_first(row, (
                "characteristic_normal_impedance",
                "z_normal",
                "z0",
            )))
            if gamma is None and z_load is not None and z_normal is not None:
                denom = z_load + z_normal
                if denom != 0.0:
                    gamma = (z_load - z_normal) / denom
            absorption = _float_or_none(_first(row, ("absorption_coefficient", "absorption")))
            incident = _float_or_none(_first(row, ("incident_intensity", "incident_power_density")))
            reflected = _float_or_none(_first(row, ("reflected_intensity", "reflected_power_density")))
            active = _float_or_none(_first(row, (
                "boundary_active_intensity_into_load",
                "active_intensity_into_load",
                "active_power_density",
            )))
            balance = _float_or_none(_first(row, ("power_balance_residual", "power_residual")))

            if z_load is None or not math.isfinite(z_load.real) or not math.isfinite(z_load.imag):
                bad_impedance.append({"index": index, "reason": "missing_or_nonfinite_impedance"})
            else:
                impedance_real_values.append(z_load.real)
                if z_load.real < -residual_tolerance:
                    bad_impedance.append({"index": index, "reason": "negative_real_impedance", "real": z_load.real})

            if gamma is None:
                bad_absorption.append({"index": index, "reason": "missing_reflection_coefficient"})
            else:
                expected_absorption = 1.0 - abs(gamma) ** 2
                if absorption is None:
                    absorption = expected_absorption
                absorption_values.append(absorption)
                if abs(absorption - expected_absorption) > residual_tolerance:
                    bad_absorption.append({
                        "index": index,
                        "reason": "absorption_mismatch",
                        "absorption": absorption,
                        "expected": expected_absorption,
                    })
                if absorption < -residual_tolerance or absorption > 1.0 + residual_tolerance:
                    bad_absorption.append({"index": index, "reason": "absorption_outside_passive_range"})

            if active is not None:
                active_powers.append(active)
                if active < -residual_tolerance:
                    bad_active_power.append({"index": index, "active_power_density": active})
            if None not in (incident, reflected, active):
                expected_balance = incident - reflected - active
                if balance is None:
                    balance = expected_balance
                power_balance_errors.append(balance)
                if abs(balance) > residual_tolerance or abs(expected_balance) > residual_tolerance:
                    bad_power_balance.append({
                        "index": index,
                        "power_balance_error": balance,
                        "expected_balance": expected_balance,
                    })
            else:
                bad_power_balance.append({"index": index, "reason": "missing_incident_reflected_or_active_power"})

        details.append({
            "index": index,
            "kind": kind,
            "case_id": None if case_id is None else str(case_id),
            "run_id": None if run_id is None else str(run_id),
            "export_id": None if export_id is None else str(export_id),
            "frequency_hz": frequency,
            "source_tool": source_tool,
            "path": path,
            "gate_policy": gate_policy,
            "status": status,
        })

    unique_case_ids = sorted(set(case_ids))
    unique_run_ids = sorted(set(run_ids))
    unique_export_ids = sorted(set(export_ids))
    expected_frequency_ok = True
    if expected_frequency_hz is not None:
        expected_frequency_ok = all(
            abs(value - float(expected_frequency_hz))
            <= frequency_tolerance * max(1.0, abs(float(expected_frequency_hz)))
            for value in frequencies
        )
    checks = {
        "required_kinds_present": all(kind_counts.get(kind, 0) > 0 for kind in required),
        "single_case_id": len(unique_case_ids) == 1 and not missing_case_id,
        "single_run_id": len(unique_run_ids) == 1 and not missing_run_id,
        "single_export_id": len(unique_export_ids) == 1 and not missing_export_id,
        "expected_case_id_ok": expected_case_id is None or unique_case_ids == [str(expected_case_id)],
        "expected_run_id_ok": expected_run_id is None or unique_run_ids == [str(expected_run_id)],
        "expected_export_id_ok": expected_export_id is None or unique_export_ids == [str(expected_export_id)],
        "frequencies_recorded": not missing_frequency and bool(frequencies),
        "frequencies_consistent": _consistent(frequencies, frequency_tolerance),
        "expected_frequency_ok": expected_frequency_ok,
        "source_tools_recorded": not missing_source_tool,
        "paths_recorded": not missing_paths,
        "all_rows_verified": not bad_upstream_status,
        "upstream_policies_match": not bad_upstream_policy,
        "impedance_passive_or_reactive": not bad_impedance,
        "absorption_matches_reflection": not bad_absorption,
        "boundary_active_power_nonnegative": not bad_active_power,
        "power_balance_ok": not bad_power_balance,
        "no_unknown_kinds": not unknown_kinds,
    }
    issues = [name for name, ok in checks.items() if not ok]

    return {
        "policy": "acoustic_impedance_power_result_package_gate",
        "status": "ok" if not issues else "needs_attention",
        "required_kinds": list(required),
        "present_kinds": dict(sorted(kind_counts.items())),
        "case_ids": unique_case_ids,
        "run_ids": unique_run_ids,
        "export_ids": unique_export_ids,
        "expected_case_id": None if expected_case_id is None else str(expected_case_id),
        "expected_run_id": None if expected_run_id is None else str(expected_run_id),
        "expected_export_id": None if expected_export_id is None else str(expected_export_id),
        "expected_frequency_hz": None if expected_frequency_hz is None else float(expected_frequency_hz),
        "frequency_hz_values": frequencies,
        "min_impedance_real": min(impedance_real_values) if impedance_real_values else None,
        "min_absorption_coefficient": min(absorption_values) if absorption_values else None,
        "max_abs_power_balance_error": max((abs(value) for value in power_balance_errors), default=None),
        "min_boundary_active_power": min(active_powers) if active_powers else None,
        "bad_impedance_rows": bad_impedance,
        "bad_absorption_rows": bad_absorption,
        "bad_power_balance_rows": bad_power_balance,
        "bad_active_power_rows": bad_active_power,
        "bad_upstream_status_rows": bad_upstream_status,
        "bad_upstream_policy_rows": bad_upstream_policy,
        "missing_case_id_rows": missing_case_id,
        "missing_run_id_rows": missing_run_id,
        "missing_export_id_rows": missing_export_id,
        "missing_frequency_rows": missing_frequency,
        "missing_source_tool_rows": missing_source_tool,
        "missing_path_rows": missing_paths,
        "unknown_kinds": unknown_kinds,
        "artifacts": details,
        "checks": checks,
        "issues": issues,
        "frequency_rtol": frequency_tolerance,
        "residual_tol": residual_tolerance,
        "version_note": (
            "Run after an acoustic impedance/reflection sub-gate so one-port "
            "FEM/BEM or radiation-boundary results cannot mix case identity, "
            "frequency, sign convention, or passive power accounting."
        ),
    }


def _complex(value):
    if isinstance(value, dict):
        return complex(float(value.get("real", 0.0)), float(value.get("imag", 0.0)))
    return complex(value)


def _complex_row(value):
    return {"real": float(value.real), "imag": float(value.imag), "abs": abs(value)}


def touchstone_sparameter_to_complex(value, fmt="MA"):
    """Convert one explicit Touchstone-style S-parameter value to complex.

    ``fmt`` must be explicit to avoid mixing real/imag, magnitude/angle, and
    dB/angle exports.  Supported formats are ``RI`` (real, imag), ``MA``
    (magnitude, phase degrees), and ``DB`` (20log10 magnitude, phase degrees).
    """

    key = str(fmt or "MA").strip().upper()
    if isinstance(value, dict):
        if key == "RI":
            pair = (value.get("real", value.get("re", 0.0)), value.get("imag", value.get("im", 0.0)))
        else:
            pair = (value.get("magnitude", value.get("mag", value.get("db", 0.0))), value.get("phase_deg", value.get("phase", 0.0)))
    else:
        pair = value
    try:
        a, b = pair
    except (TypeError, ValueError) as exc:
        raise ValueError("Touchstone value must contain two numbers") from exc

    x = float(a)
    y = float(b)
    if key == "RI":
        return complex(x, y)
    if key == "MA":
        magnitude = x
    elif key == "DB":
        magnitude = 10.0 ** (x / 20.0)
    else:
        raise ValueError(f"unsupported Touchstone format: {fmt!r}")
    phase = math.radians(y)
    return complex(magnitude * math.cos(phase), magnitude * math.sin(phase))


def two_port_sparameter_health(s11, s21, s12=None, s22=None, tol=1.0e-9):
    """Return reciprocity/passivity metrics for a 2x2 S-parameter matrix."""

    a = _complex(s11)
    b = _complex(s21)
    c = b if s12 is None else _complex(s12)
    d = a if s22 is None else _complex(s22)
    # S^H S for a 2x2 matrix.
    m00 = abs(a) ** 2 + abs(c) ** 2
    m11 = abs(b) ** 2 + abs(d) ** 2
    m01 = a.conjugate() * b + c.conjugate() * d
    trace = m00 + m11
    det = m00 * m11 - abs(m01) ** 2
    disc = max(trace * trace - 4.0 * det, 0.0)
    lambda_max = 0.5 * (trace + math.sqrt(disc))
    reciprocity_error = abs(b - c)
    passive_margin = 1.0 - lambda_max
    return {
        "policy": "two_port_sparameter_reciprocity_passivity_gate",
        "s11": _complex_row(a),
        "s21": _complex_row(b),
        "s12": _complex_row(c),
        "s22": _complex_row(d),
        "reciprocity_error": reciprocity_error,
        "max_singular_value_squared": lambda_max,
        "passive_margin": passive_margin,
        "reciprocal": reciprocity_error <= float(tol),
        "passive": lambda_max <= 1.0 + float(tol),
        "status": "ok" if reciprocity_error <= float(tol) and lambda_max <= 1.0 + float(tol) else "needs_attention",
        "tol": float(tol),
    }


def touchstone_power_wave_balance_gate(
    s11,
    s21,
    s12=None,
    s22=None,
    data_format="MA",
    power_balance_basis="power_waves_unit_incident_port",
    expected_power_balance_basis="power_waves_unit_incident_port",
    export_artifact_id=None,
    expected_export_artifact_id=None,
    result_set_id=None,
    expected_result_set_id=None,
    tol=1.0e-9,
):
    """Check two-port Touchstone power-wave balance per unit incident port.

    This gate keeps the S-parameter basis explicit before a row is reused for
    power accounting.  For power waves with one unit incident port, each column
    of S maps incident power to reflected/transmitted output power; the
    remaining fraction is absorbed or dissipated by a passive network.
    """

    def norm(value):
        return str(value or "").strip().lower().replace("-", "_").replace(" ", "_")

    tolerance = float(tol)
    if tolerance < 0.0:
        raise ValueError("tol must be non-negative")

    fmt = str(data_format or "MA").strip().upper()
    a = touchstone_sparameter_to_complex(s11, fmt=fmt)
    b = touchstone_sparameter_to_complex(s21, fmt=fmt)
    c = b if s12 is None else touchstone_sparameter_to_complex(s12, fmt=fmt)
    d = a if s22 is None else touchstone_sparameter_to_complex(s22, fmt=fmt)
    health = two_port_sparameter_health(a, b, s12=c, s22=d, tol=tolerance)

    basis = norm(power_balance_basis)
    expected_basis = None if expected_power_balance_basis is None else norm(expected_power_balance_basis)
    export_artifact = str(export_artifact_id or "").strip()
    expected_export_artifact = (
        None if expected_export_artifact_id is None else str(expected_export_artifact_id).strip()
    )
    result_set = str(result_set_id or "").strip()
    expected_result_set = (
        None if expected_result_set_id is None else str(expected_result_set_id).strip()
    )

    def power_row(port, reflected, transmitted):
        reflected_power = abs(reflected) ** 2
        transmitted_power = abs(transmitted) ** 2
        output_power = reflected_power + transmitted_power
        absorbed_power = 1.0 - output_power
        return {
            "excitation_port": port,
            "incident_power_fraction": 1.0,
            "reflected_power_fraction": reflected_power,
            "transmitted_power_fraction": transmitted_power,
            "output_power_sum": output_power,
            "absorbed_power_fraction": absorbed_power,
            "power_balance_residual": abs(1.0 - output_power - absorbed_power),
        }

    rows = [
        power_row("port1", a, b),
        power_row("port2", d, c),
    ]
    max_residual = max(row["power_balance_residual"] for row in rows)
    min_absorbed = min(row["absorbed_power_fraction"] for row in rows)
    checks = {
        "touchstone_format_recorded": fmt in {"RI", "MA", "DB"},
        "power_balance_basis_recorded": bool(basis),
        "power_balance_basis_matches_expected": expected_basis is None or basis == expected_basis,
        "column_power_not_active": min_absorbed >= -tolerance,
        "power_balance_rows_close": max_residual <= tolerance,
        "sparameter_passivity_ok": bool(health["passive"]),
        "sparameter_reciprocity_ok": bool(health["reciprocal"]),
        "export_artifact_id_recorded": expected_export_artifact is None or bool(export_artifact),
        "expected_export_artifact_id_matches": expected_export_artifact is None
        or export_artifact == expected_export_artifact,
        "result_set_id_recorded": expected_result_set is None or bool(result_set),
        "expected_result_set_id_matches": expected_result_set is None
        or result_set == expected_result_set,
    }
    return {
        "policy": "touchstone_power_wave_balance_gate",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "data_format": fmt,
        "power_balance_basis": basis,
        "expected_power_balance_basis": expected_basis,
        "export_artifact_id": export_artifact or None,
        "expected_export_artifact_id": expected_export_artifact,
        "result_set_id": result_set or None,
        "expected_result_set_id": expected_result_set,
        "rows": rows,
        "max_power_balance_residual": max_residual,
        "min_absorbed_power_fraction": min_absorbed,
        "health": health,
        "checks": checks,
        "tol": tolerance,
        "version_note": (
            "Run after explicit RI/MA/DB normalization.  The compact power "
            "ledger assumes power-wave S-parameters with one unit incident port."
        ),
    }


def one_port_match_quality_gate(s11=None, reflection_coefficient=None, vswr=None, return_loss_db=None, tol=1.0e-12):
    """Convert one one-port match descriptor into VSWR, return loss, and mismatch loss.

    Exactly one of ``s11``, ``reflection_coefficient``, ``vswr``, or
    ``return_loss_db`` must be supplied.  This gate is useful before a CST,
    Touchstone, VNA, or open-solver port row is promoted into a circuit or
    motor-drive protection validation record.
    """

    provided = [
        s11 is not None,
        reflection_coefficient is not None,
        vswr is not None,
        return_loss_db is not None,
    ]
    if sum(provided) != 1:
        raise ValueError("provide exactly one of s11, reflection_coefficient, vswr, or return_loss_db")
    tolerance = float(tol)
    if tolerance < 0.0:
        raise ValueError("tol must be non-negative")

    source = None
    if s11 is not None:
        gamma = abs(_complex(s11))
        source = "s11"
    elif reflection_coefficient is not None:
        gamma = float(reflection_coefficient)
        source = "reflection_coefficient"
    elif vswr is not None:
        standing_wave_ratio = float(vswr)
        if standing_wave_ratio < 1.0:
            raise ValueError("vswr must be >= 1")
        gamma = (standing_wave_ratio - 1.0) / (standing_wave_ratio + 1.0)
        source = "vswr"
    else:
        rl = float(return_loss_db)
        if rl < 0.0:
            raise ValueError("return_loss_db must be >= 0")
        gamma = 10.0 ** (-rl / 20.0)
        source = "return_loss_db"

    if gamma < -tolerance:
        raise ValueError("reflection coefficient magnitude must be non-negative")
    gamma = max(gamma, 0.0)
    passive_reflection_ok = gamma <= 1.0 + tolerance
    gamma_for_ratios = min(gamma, 1.0)
    vswr_out = math.inf if gamma_for_ratios >= 1.0 else (1.0 + gamma_for_ratios) / (1.0 - gamma_for_ratios)
    return_loss_out = math.inf if gamma == 0.0 else -20.0 * math.log10(gamma)
    transmitted_power_fraction = 1.0 - gamma * gamma
    mismatch_loss = (
        -10.0 * math.log10(transmitted_power_fraction)
        if transmitted_power_fraction > 0.0
        else math.inf
    )
    gamma_from_vswr = 1.0 if math.isinf(vswr_out) else (vswr_out - 1.0) / (vswr_out + 1.0)
    gamma_from_return_loss = 0.0 if math.isinf(return_loss_out) else 10.0 ** (-return_loss_out / 20.0)
    checks = {
        "passive_reflection_ok": passive_reflection_ok,
        "vswr_round_trip_ok": abs(gamma_from_vswr - gamma_for_ratios) <= tolerance,
        "return_loss_round_trip_ok": abs(gamma_from_return_loss - gamma) <= tolerance,
        "transmitted_power_fraction_nonnegative": transmitted_power_fraction >= -tolerance,
    }
    return {
        "policy": "one_port_match_quality_gate",
        "source": source,
        "reflection_coefficient": gamma,
        "s11_abs": gamma,
        "vswr": vswr_out,
        "return_loss_db": return_loss_out,
        "mismatch_loss_db": mismatch_loss,
        "transmitted_power_fraction": transmitted_power_fraction,
        "reflected_power_fraction": gamma * gamma,
        "checks": checks,
        "status": "ok" if all(checks.values()) else "needs_attention",
        "tol": tolerance,
    }


def touchstone_port_metadata_gate(
    metadata,
    required_ports=("P1", "P2"),
    expected_network="S",
    data_format=None,
    frequency_unit=None,
    reference_impedance_ohm=None,
    port_order=None,
    reference_plane=None,
    port_mode_basis=None,
    tol=1.0e-12,
):
    """Check Touchstone option-line and port metadata before parsing rows.

    CST, VNA, and notebook exports can all produce plausible S-parameter
    values while silently swapping ports, losing the reference impedance, or
    mixing RI/MA/DB numeric columns.  This gate keeps the table contract ahead
    of row-level passivity, reciprocity, and match checks.
    """

    if not isinstance(metadata, dict):
        raise ValueError("metadata must be a mapping")
    tolerance = float(tol)
    if tolerance < 0.0:
        raise ValueError("tol must be non-negative")

    def pick(*keys, default=None):
        for key in keys:
            if key in metadata and metadata[key] is not None:
                return metadata[key]
        return default

    def normalize_ports(value):
        if value is None:
            return []
        if isinstance(value, str):
            return [part.strip() for part in value.replace(";", ",").split(",") if part.strip()]
        ports = []
        for item in value:
            if isinstance(item, dict):
                raw = item.get("name", item.get("label", item.get("id", "")))
            else:
                raw = item
            text = str(raw).strip()
            if text:
                ports.append(text)
        return ports

    def normalize_unit(value):
        if value is None:
            return None
        key = str(value).strip().lower()
        aliases = {
            "hz": "Hz",
            "khz": "kHz",
            "mhz": "MHz",
            "ghz": "GHz",
        }
        return aliases.get(key)

    ports = normalize_ports(pick("ports", "port_names", default=[]))
    expected_required = [str(item).strip() for item in required_ports if str(item).strip()]
    expected_order = normalize_ports(port_order)
    if not expected_order:
        expected_order = normalize_ports(pick("port_order", default=[]))

    port_count_raw = pick("port_count", "n_ports", default=None)
    port_count = None if port_count_raw is None else int(port_count_raw)
    network_present = any(key in metadata and metadata[key] is not None for key in ("network_parameter", "network", "parameter"))
    network = pick("network_parameter", "network", "parameter", default=None)
    network_key = None if network is None else str(network).strip().upper()
    expected_network_key = None if expected_network is None else str(expected_network).strip().upper()
    format_present = any(key in metadata and metadata[key] is not None for key in ("data_format", "format"))
    fmt = pick("data_format", "format", default=data_format)
    fmt_key = None if fmt is None else str(fmt).strip().upper()
    frequency_present = any(key in metadata and metadata[key] is not None for key in ("frequency_unit", "freq_unit"))
    freq = pick("frequency_unit", "freq_unit", default=frequency_unit)
    freq_unit = normalize_unit(freq)
    expected_freq_unit = normalize_unit(frequency_unit)
    z0_present = any(key in metadata and metadata[key] is not None for key in ("reference_impedance_ohm", "reference_impedance", "z0_ohm", "z0"))
    z_raw = pick("reference_impedance_ohm", "reference_impedance", "z0_ohm", "z0", default=reference_impedance_ohm)
    z0 = None if z_raw is None else float(z_raw)
    expected_z0 = None if reference_impedance_ohm is None else float(reference_impedance_ohm)
    ref_plane_raw = pick("reference_plane", "deembed_reference_plane", "port_reference_plane", default=None)
    ref_plane = None if ref_plane_raw is None else str(ref_plane_raw).strip()
    expected_ref_plane = None if reference_plane is None else str(reference_plane).strip()
    mode_basis_raw = pick("port_mode_basis", "mode_basis", "port_modes", default=None)
    mode_basis = None if mode_basis_raw is None else str(mode_basis_raw).strip()
    expected_mode_basis = None if port_mode_basis is None else str(port_mode_basis).strip()

    port_counts_match = port_count is None or port_count == len(ports)
    required_present = all(item in ports for item in expected_required)
    order_ok = True if not expected_order else ports[:len(expected_order)] == expected_order
    network_ok = network_key in {"S", "Y", "Z", "H", "G"}
    if expected_network_key is not None:
        network_ok = network_ok and network_key == expected_network_key
    format_ok = fmt_key in {"RI", "MA", "DB"}
    if data_format is not None:
        format_ok = format_ok and fmt_key == str(data_format).strip().upper()
    frequency_ok = freq_unit in {"Hz", "kHz", "MHz", "GHz"}
    if expected_freq_unit is not None:
        frequency_ok = frequency_ok and freq_unit == expected_freq_unit
    z0_ok = z0 is not None and z0 > 0.0
    if expected_z0 is not None:
        z0_ok = z0_ok and abs(z0 - expected_z0) <= tolerance * max(1.0, abs(expected_z0))
    ref_plane_recorded = expected_ref_plane is None or bool(ref_plane)
    ref_plane_ok = expected_ref_plane is None or (
        bool(ref_plane) and ref_plane.lower() == expected_ref_plane.lower()
    )
    mode_basis_recorded = expected_mode_basis is None or bool(mode_basis)
    mode_basis_ok = expected_mode_basis is None or (
        bool(mode_basis) and mode_basis.lower() == expected_mode_basis.lower()
    )

    checks = {
        "port_names_recorded": bool(ports),
        "ports_unique": len(ports) == len(set(ports)),
        "port_count_matches_names": port_counts_match,
        "required_ports_present": required_present,
        "port_order_matches_expected": order_ok,
        "network_parameter_recorded": network_present,
        "network_parameter_matches_expected": network_ok,
        "touchstone_format_recorded": format_present,
        "touchstone_format_matches_expected": format_ok,
        "frequency_unit_recorded": frequency_present,
        "frequency_unit_matches_expected": frequency_ok,
        "reference_impedance_recorded": z0_present,
        "reference_impedance_matches_expected": z0_ok,
        "reference_plane_recorded": ref_plane_recorded,
        "reference_plane_matches_expected": ref_plane_ok,
        "port_mode_basis_recorded": mode_basis_recorded,
        "port_mode_basis_matches_expected": mode_basis_ok,
    }
    return {
        "policy": "touchstone_port_metadata_gate",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "ports": ports,
        "port_count": port_count,
        "required_ports": expected_required,
        "expected_port_order": expected_order,
        "network_parameter": network_key,
        "expected_network": expected_network_key,
        "data_format": fmt_key,
        "frequency_unit": freq_unit,
        "reference_impedance_ohm": z0,
        "reference_plane": ref_plane or None,
        "expected_reference_plane": expected_ref_plane,
        "port_mode_basis": mode_basis or None,
        "expected_port_mode_basis": expected_mode_basis,
        "checks": checks,
        "tol": tolerance,
        "notes": [
            "Run this before converting RI/MA/DB numeric pairs into complex S-parameters.",
            "A passive row can still be unusable if ports are swapped, z0 is missing, the reference plane moved, or the port-mode basis is ambiguous.",
        ],
    }


def farfield_pattern_metadata_gate(
    metadata,
    expected_quantity="gain",
    expected_quantity_unit="dBi",
    expected_angle_unit="deg",
    expected_coordinate_system="spherical",
    expected_polarization_basis="theta_phi",
    expected_normalization="accepted_power",
    required_components=("Etheta", "Ephi"),
    required_phi_values_deg=None,
    min_theta_span_deg=180.0,
    frequency_hz=None,
    tol=1.0e-9,
):
    """Check far-field/radiation-pattern table metadata before row values.

    CST, ngsolve.bem, antenna notebooks, and measurement exports can all report
    plausible gain/directivity numbers while mixing angle units, theta/phi cuts,
    polarization bases, or power normalizations.  This public-safe gate keeps
    those table contracts visible before comparing far-field rows.
    """

    if not isinstance(metadata, dict):
        raise ValueError("metadata must be a mapping")
    tolerance = float(tol)
    if tolerance < 0.0:
        raise ValueError("tol must be non-negative")

    def pick(*keys, default=None):
        for key in keys:
            if key in metadata and metadata[key] is not None:
                return metadata[key]
        return default

    def as_list(value):
        if value is None:
            return []
        if isinstance(value, str):
            parts = value.replace(";", ",").split(",")
            return [part.strip() for part in parts if part.strip()]
        if isinstance(value, (list, tuple)):
            return list(value)
        return [value]

    def number_list(value):
        return [float(item) for item in as_list(value)]

    def angle_unit_label(value):
        if value is None:
            return None
        key = str(value).strip().lower()
        aliases = {
            "deg": "deg",
            "degree": "deg",
            "degrees": "deg",
            "rad": "rad",
            "radian": "rad",
            "radians": "rad",
        }
        return aliases.get(key)

    def values_to_deg(values, unit):
        if unit == "rad":
            return [math.degrees(value) for value in values]
        return list(values)

    def normalize_token(value):
        if value is None:
            return None
        return str(value).strip().replace("-", "_").replace("/", "_").lower()

    def normalize_basis(value):
        token = normalize_token(value)
        aliases = {
            "theta_phi": "theta_phi",
            "etheta_ephi": "theta_phi",
            "spherical_theta_phi": "theta_phi",
            "lhcp_rhcp": "circular",
            "rhcp_lhcp": "circular",
            "linear": "linear",
        }
        return aliases.get(token, token)

    def contains_angle(values, target):
        return any(abs(value - target) <= max(tolerance, 1.0e-9) for value in values)

    angle_unit_present = any(key in metadata and metadata[key] is not None for key in ("angle_unit", "angle_units"))
    angle_unit_raw = pick("angle_unit", "angle_units", default=None)
    angle_unit = angle_unit_label(angle_unit_raw)
    expected_angle = angle_unit_label(expected_angle_unit)
    theta_values = number_list(pick("theta_values_deg", "theta_deg", "theta_values", "theta_grid", default=[]))
    phi_values = number_list(pick("phi_values_deg", "phi_deg", "phi_values", "phi_grid", default=[]))
    theta_deg = values_to_deg(theta_values, angle_unit)
    phi_deg = values_to_deg(phi_values, angle_unit)
    components = [str(item).strip() for item in as_list(pick("field_components", "components", default=[])) if str(item).strip()]
    required = [str(item).strip() for item in required_components if str(item).strip()]
    quantity_present = any(key in metadata and metadata[key] is not None for key in ("quantity", "field_quantity", "pattern_quantity"))
    quantity = normalize_token(pick("quantity", "field_quantity", "pattern_quantity", default=None))
    expected_quantity_key = normalize_token(expected_quantity)
    quantity_unit_present = any(key in metadata and metadata[key] is not None for key in ("quantity_unit", "gain_unit", "directivity_unit", "field_unit"))
    quantity_unit = str(pick("quantity_unit", "gain_unit", "directivity_unit", "field_unit", default=None) or "").strip()
    coordinate_present = any(key in metadata and metadata[key] is not None for key in ("coordinate_system", "coordinates"))
    coordinate_system = normalize_token(pick("coordinate_system", "coordinates", default=None))
    expected_coordinates = normalize_token(expected_coordinate_system)
    basis_present = any(key in metadata and metadata[key] is not None for key in ("polarization_basis", "basis", "polarization"))
    basis = normalize_basis(pick("polarization_basis", "basis", "polarization", default=None))
    expected_basis = normalize_basis(expected_polarization_basis)
    normalization_present = any(key in metadata and metadata[key] is not None for key in ("normalization", "power_normalization"))
    normalization = normalize_token(pick("normalization", "power_normalization", default=None))
    expected_norm = normalize_token(expected_normalization)

    freq_present = any(key in metadata and metadata[key] is not None for key in ("frequency_hz", "frequency_Hz", "freq_hz"))
    freq_raw = pick("frequency_hz", "frequency_Hz", "freq_hz", default=frequency_hz)
    freq = None if freq_raw is None else float(freq_raw)
    expected_freq = None if frequency_hz is None else float(frequency_hz)
    row_count_raw = pick("row_count", "n_rows", default=None)
    row_count = None if row_count_raw is None else int(row_count_raw)
    expected_grid_rows = len(theta_deg) * max(1, len(phi_deg)) if theta_deg else 0
    required_phi = [] if required_phi_values_deg is None else [float(value) for value in required_phi_values_deg]

    theta_finite = all(math.isfinite(value) for value in theta_deg)
    phi_finite = all(math.isfinite(value) for value in phi_deg)
    theta_strict = all(a < b for a, b in zip(theta_deg, theta_deg[1:]))
    phi_unique = len(phi_deg) == len(set(round(value, 12) for value in phi_deg))
    theta_span = (max(theta_deg) - min(theta_deg)) if theta_deg else 0.0
    theta_range_ok = bool(theta_deg) and theta_finite and min(theta_deg) >= -tolerance and max(theta_deg) <= 180.0 + tolerance
    phi_range_ok = bool(phi_deg) and phi_finite and min(phi_deg) >= -360.0 - tolerance and max(phi_deg) <= 360.0 + tolerance
    required_phi_present = all(contains_angle(phi_deg, value) for value in required_phi)

    checks = {
        "frequency_recorded": freq_present,
        "frequency_positive": freq is not None and freq > 0.0,
        "frequency_matches_expected": expected_freq is None or (
            freq is not None and abs(freq - expected_freq) <= tolerance * max(1.0, abs(expected_freq))
        ),
        "angle_unit_recorded": angle_unit_present and angle_unit is not None,
        "angle_unit_matches_expected": angle_unit == expected_angle,
        "coordinate_system_recorded": coordinate_present and coordinate_system is not None,
        "coordinate_system_matches_expected": coordinate_system == expected_coordinates,
        "polarization_basis_recorded": basis_present and basis is not None,
        "polarization_basis_matches_expected": basis == expected_basis,
        "quantity_recorded": quantity_present and quantity is not None,
        "quantity_matches_expected": quantity == expected_quantity_key,
        "quantity_unit_recorded": quantity_unit_present and bool(quantity_unit),
        "quantity_unit_matches_expected": quantity_unit == str(expected_quantity_unit),
        "normalization_recorded": normalization_present and normalization is not None,
        "normalization_matches_expected": normalization == expected_norm,
        "field_components_recorded": bool(components),
        "required_components_present": all(item in components for item in required),
        "theta_grid_recorded": bool(theta_deg),
        "theta_grid_strictly_increasing": theta_strict,
        "theta_range_degrees_ok": theta_range_ok,
        "theta_span_ok": theta_span + tolerance >= float(min_theta_span_deg),
        "phi_grid_recorded": bool(phi_deg),
        "phi_values_unique": phi_unique,
        "phi_range_degrees_ok": phi_range_ok,
        "required_phi_values_present": required_phi_present,
        "row_count_recorded": row_count is not None,
        "row_count_covers_grid": row_count is not None and row_count >= expected_grid_rows,
    }
    return {
        "policy": "farfield_pattern_metadata_gate",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "frequency_hz": freq,
        "expected_frequency_hz": expected_freq,
        "angle_unit": angle_unit,
        "theta_values_deg": theta_deg,
        "phi_values_deg": phi_deg,
        "theta_span_deg": theta_span,
        "coordinate_system": coordinate_system,
        "polarization_basis": basis,
        "quantity": quantity,
        "quantity_unit": quantity_unit,
        "normalization": normalization,
        "field_components": components,
        "required_components": required,
        "required_phi_values_deg": required_phi,
        "row_count": row_count,
        "expected_grid_rows": expected_grid_rows,
        "checks": checks,
        "tol": tolerance,
        "notes": [
            "Run this before comparing gain, directivity, RCS, or Etheta/Ephi far-field rows.",
            "A radiation pattern can look plausible while angle units, cuts, polarization basis, or normalization are wrong.",
        ],
    }


def farfield_lobe_notebook_handoff_gate(
    metadata,
    row,
    lobe_id_key="lobe_id",
    tol=1.0e-9,
):
    """Bundle far-field metadata and one lobe row before interface use.

    An application block or result notebook should not receive a naked gain scalar. The lobe row must
    carry a stable row identity, frequency, theta/phi location, polarization
    basis, accepted-power normalization, gain/directivity units, and the
    radiation-efficiency identity ``G = eta_rad * D``.
    """

    if not isinstance(row, dict):
        raise ValueError("row must be a mapping")
    tolerance = float(tol)
    if tolerance < 0.0:
        raise ValueError("tol must be non-negative")

    meta_gate = farfield_pattern_metadata_gate(metadata, tol=tolerance)

    def pick(*keys, default=None):
        for key in keys:
            if key in row and row[key] is not None:
                return row[key]
        return default

    def normalize_token(value):
        if value is None:
            return None
        return str(value).strip().replace("-", "_").replace("/", "_").lower()

    def normalize_basis(value):
        token = normalize_token(value)
        aliases = {
            "theta_phi": "theta_phi",
            "etheta_ephi": "theta_phi",
            "spherical_theta_phi": "theta_phi",
            "lhcp_rhcp": "circular",
            "rhcp_lhcp": "circular",
            "linear": "linear",
        }
        return aliases.get(token, token)

    def contains_angle(values, target):
        return any(abs(value - target) <= max(tolerance, 1.0e-9) for value in values)

    lobe_id = pick(lobe_id_key, "row_id", "case_id", "label", default=None)
    row_frequency = pick("frequency_hz", "frequency_Hz", "freq_hz", default=None)
    row_frequency = None if row_frequency is None else float(row_frequency)
    theta = pick("theta_deg", "theta", default=None)
    phi = pick("phi_deg", "phi", default=None)
    theta = None if theta is None else float(theta)
    phi = None if phi is None else float(phi)
    row_basis = normalize_basis(pick("polarization_basis", "basis", "polarization", default=None))
    row_norm = normalize_token(pick("normalization", "power_normalization", default=None))
    gain_unit = str(pick("gain_unit", "gain_quantity_unit", default="") or "").strip()
    directivity_unit = str(pick("directivity_unit", default="") or "").strip()

    gain_raw = pick("gain_dbi", "gain_dBi", "gain_db", "gain_DB", default=None)
    directivity_raw = pick("directivity_dbi", "directivity_dBi", "directivity_db", "directivity_DB", default=None)
    radiated_raw = pick("radiated_power_w", "radiated_power_W", "prad_w", "P_rad_W", default=None)
    accepted_raw = pick("accepted_power_w", "accepted_power_W", "pacc_w", "P_acc_W", default=None)

    checks = {
        "metadata_gate_ok": meta_gate["status"] == "ok",
        "lobe_id_recorded": lobe_id is not None and str(lobe_id).strip() != "",
        "row_frequency_recorded": row_frequency is not None,
        "row_frequency_matches_metadata": (
            row_frequency is not None
            and meta_gate["frequency_hz"] is not None
            and abs(row_frequency - meta_gate["frequency_hz"]) <= tolerance * max(1.0, abs(meta_gate["frequency_hz"]))
        ),
        "theta_recorded": theta is not None,
        "theta_on_export_grid": theta is not None and contains_angle(meta_gate["theta_values_deg"], theta),
        "phi_recorded": phi is not None,
        "phi_on_export_grid": phi is not None and contains_angle(meta_gate["phi_values_deg"], phi),
        "polarization_basis_recorded": row_basis is not None,
        "polarization_basis_matches_metadata": row_basis is not None and row_basis == meta_gate["polarization_basis"],
        "normalization_recorded": row_norm is not None,
        "normalization_is_accepted_power": row_norm == "accepted_power",
        "normalization_matches_metadata": row_norm is not None and row_norm == meta_gate["normalization"],
        "gain_unit_recorded": bool(gain_unit),
        "gain_unit_is_dbi": gain_unit == "dBi",
        "directivity_unit_recorded": bool(directivity_unit),
        "directivity_unit_is_dbi": directivity_unit == "dBi",
        "gain_recorded": gain_raw is not None,
        "directivity_recorded": directivity_raw is not None,
        "radiated_power_recorded": radiated_raw is not None,
        "accepted_power_recorded": accepted_raw is not None,
        "radiated_power_nonnegative": False,
        "accepted_power_positive": False,
        "radiation_efficiency_in_0_1": False,
        "gain_not_above_directivity": False,
        "gain_matches_directivity_times_efficiency": False,
    }
    values = {
        "gain_dbi": None,
        "directivity_dbi": None,
        "gain_linear": None,
        "directivity_linear": None,
        "radiated_power_w": None,
        "accepted_power_w": None,
        "radiation_efficiency": None,
        "expected_gain_linear": None,
        "efficiency_from_gain_over_directivity": None,
        "gain_relative_error": None,
    }
    if None not in (gain_raw, directivity_raw, radiated_raw, accepted_raw):
        gain_dbi = float(gain_raw)
        directivity_dbi = float(directivity_raw)
        radiated = float(radiated_raw)
        accepted = float(accepted_raw)
        gain_linear = 10.0 ** (gain_dbi / 10.0)
        directivity_linear = 10.0 ** (directivity_dbi / 10.0)
        checks["radiated_power_nonnegative"] = radiated >= 0.0
        checks["accepted_power_positive"] = accepted > 0.0
        if accepted > 0.0 and directivity_linear > 0.0:
            eta = radiated / accepted
            expected_gain = directivity_linear * eta
            rel_error = abs(gain_linear - expected_gain) / max(abs(expected_gain), abs(gain_linear), 1.0)
            checks["radiation_efficiency_in_0_1"] = 0.0 <= eta <= 1.0 + tolerance
            checks["gain_not_above_directivity"] = gain_linear <= directivity_linear + tolerance
            checks["gain_matches_directivity_times_efficiency"] = rel_error <= tolerance
            values.update(
                {
                    "gain_dbi": gain_dbi,
                    "directivity_dbi": directivity_dbi,
                    "gain_linear": gain_linear,
                    "directivity_linear": directivity_linear,
                    "radiated_power_w": radiated,
                    "accepted_power_w": accepted,
                    "radiation_efficiency": eta,
                    "expected_gain_linear": expected_gain,
                    "efficiency_from_gain_over_directivity": gain_linear / directivity_linear,
                    "gain_relative_error": rel_error,
                }
            )

    return {
        "policy": "farfield_lobe_notebook_handoff_gate",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "metadata_gate": meta_gate,
        "lobe_id_key": lobe_id_key,
        "lobe_id": lobe_id,
        "frequency_hz": row_frequency,
        "theta_deg": theta,
        "phi_deg": phi,
        "polarization_basis": row_basis,
        "normalization": row_norm,
        "gain_unit": gain_unit,
        "directivity_unit": directivity_unit,
        "checks": checks,
        "tol": tolerance,
        "notes": [
            "Run this after the far-field export metadata gate and before an application block or result notebook plots or ranks lobes.",
            "The row identity, angular location, accepted-power normalization, and G=eta_rad*D identity must travel together.",
        ],
        **values,
    }

def cst_result_export_package_gate(
    artifacts,
    expected_project_id=None,
    expected_run_id=None,
    expected_export_id=None,
    expected_frequency_Hz=None,
    required_kinds=("touchstone_metadata", "touchstone_row", "farfield_metadata", "farfield_lobe"),
    frequency_rtol=1.0e-12,
):
    """Check that CST RF result artifacts belong to one export package."""

    rows_in = list(artifacts)
    if not rows_in:
        raise ValueError("artifacts must not be empty")

    def _norm(value):
        return str(value or "").strip().lower().replace("-", "_").replace(" ", "_")

    def _first(row, names):
        for name in names:
            if name in row and row[name] is not None:
                return row[name]
        return None

    required = tuple(_norm(kind) for kind in required_kinds)
    expected_policies = {
        "touchstone_metadata": {"touchstone_port_metadata_gate"},
        "touchstone_row": {"touchstone_row_solver_ready_preflight_gate", "two_port_sparameter_health"},
        "farfield_metadata": {"farfield_pattern_metadata_gate"},
        "farfield_lobe": {"farfield_lobe_notebook_handoff_gate"},
    }
    tolerance = float(frequency_rtol)
    if not required:
        raise ValueError("required_kinds must not be empty")
    if tolerance < 0.0:
        raise ValueError("frequency_rtol must be non-negative")

    details = []
    kind_counts = {}
    project_ids = []
    run_ids = []
    export_ids = []
    frequencies = []
    missing_project_id = []
    missing_run_id = []
    missing_export_id = []
    missing_frequency = []
    bad_source_tool = []
    missing_paths = []
    unknown_kinds = []
    bad_upstream_status = []
    bad_upstream_policy = []

    for index, row in enumerate(rows_in, start=1):
        if not isinstance(row, dict):
            raise ValueError("each artifact must be a dictionary")
        kind = _norm(_first(row, ("kind", "artifact_kind", "type")))
        project_id = _first(row, ("project_id", "cst_project_id", "model_id"))
        run_id = _first(row, ("run_id", "solver_run_id", "study_id"))
        export_id = _first(row, ("export_id", "result_set_id", "dataset_id"))
        source_tool = _first(row, ("source_tool", "tool", "source"))
        source_tool_norm = _norm(source_tool)
        path = _first(row, ("path", "file", "artifact_path", "table_path"))
        frequency = _first(row, ("frequency_Hz", "frequency_hz", "freq_Hz", "freq_hz"))
        gate_policy = _first(row, ("gate_policy", "policy", "validator"))
        gate_policy_norm = _norm(gate_policy)
        status = _first(row, ("status", "gate_status", "validation_status"))
        status_norm = _norm(status)
        pass_flag = bool(row.get("pass", False))
        if kind:
            kind_counts[kind] = kind_counts.get(kind, 0) + 1
        else:
            unknown_kinds.append(index)
        if kind and kind not in set(required) | set(expected_policies):
            unknown_kinds.append(kind)
        if not project_id:
            missing_project_id.append(index)
        else:
            project_ids.append(str(project_id))
        if not run_id:
            missing_run_id.append(index)
        else:
            run_ids.append(str(run_id))
        if not export_id:
            missing_export_id.append(index)
        else:
            export_ids.append(str(export_id))
        if frequency is None:
            missing_frequency.append(index)
        else:
            frequencies.append(float(frequency))
        if source_tool_norm not in {"cst", "cst_studio", "cst_studio_suite"}:
            bad_source_tool.append({"index": index, "source_tool": source_tool})
        if not path:
            missing_paths.append(index)
        if not (pass_flag or status_norm in {"ok", "pass", "passed", "verified"}):
            bad_upstream_status.append({"index": index, "kind": kind, "status": status})
        if kind in expected_policies and gate_policy_norm not in expected_policies[kind]:
            bad_upstream_policy.append({
                "index": index,
                "kind": kind,
                "gate_policy": gate_policy,
                "expected": sorted(expected_policies[kind]),
            })
        details.append({
            "index": index,
            "kind": kind,
            "project_id": None if project_id is None else str(project_id),
            "run_id": None if run_id is None else str(run_id),
            "export_id": None if export_id is None else str(export_id),
            "frequency_Hz": None if frequency is None else float(frequency),
            "source_tool": source_tool,
            "path": path,
            "gate_policy": gate_policy,
            "status": status,
            "pass": pass_flag,
        })

    unique_project_ids = sorted(set(project_ids))
    unique_run_ids = sorted(set(run_ids))
    unique_export_ids = sorted(set(export_ids))
    max_frequency_rel_span = 0.0
    if frequencies:
        max_frequency_rel_span = abs(max(frequencies) - min(frequencies)) / max(max(abs(f) for f in frequencies), 1.0)
    checks = {
        "required_kinds_present": set(required).issubset(set(kind_counts)),
        "no_unknown_kinds": not unknown_kinds,
        "project_ids_present": not missing_project_id,
        "project_ids_unique": len(unique_project_ids) == 1,
        "run_ids_present": not missing_run_id,
        "run_ids_unique": len(unique_run_ids) == 1,
        "export_ids_present": not missing_export_id,
        "export_ids_unique": len(unique_export_ids) == 1,
        "frequencies_present": not missing_frequency,
        "frequencies_match": bool(frequencies) and max_frequency_rel_span <= tolerance,
        "source_tool_is_cst": not bad_source_tool,
        "paths_present": not missing_paths,
        "upstream_gate_status_ok": not bad_upstream_status,
        "upstream_gate_policy_known": not bad_upstream_policy,
    }
    if expected_project_id is not None:
        checks["expected_project_id_matches"] = unique_project_ids == [str(expected_project_id)]
    if expected_run_id is not None:
        checks["expected_run_id_matches"] = unique_run_ids == [str(expected_run_id)]
    if expected_export_id is not None:
        checks["expected_export_id_matches"] = unique_export_ids == [str(expected_export_id)]
    if expected_frequency_Hz is not None:
        expected_f = float(expected_frequency_Hz)
        checks["expected_frequency_matches"] = (
            bool(frequencies)
            and max(abs(freq - expected_f) / max(abs(freq), abs(expected_f), 1.0) for freq in frequencies) <= tolerance
        )

    return {
        "policy": "cst_result_export_package_gate",
        "required_kinds": list(required),
        "present_kinds": dict(sorted(kind_counts.items())),
        "project_ids": unique_project_ids,
        "run_ids": unique_run_ids,
        "export_ids": unique_export_ids,
        "frequencies_Hz": sorted(set(frequencies)),
        "max_frequency_rel_span": max_frequency_rel_span,
        "expected_project_id": None if expected_project_id is None else str(expected_project_id),
        "expected_run_id": None if expected_run_id is None else str(expected_run_id),
        "expected_export_id": None if expected_export_id is None else str(expected_export_id),
        "expected_frequency_Hz": None if expected_frequency_Hz is None else float(expected_frequency_Hz),
        "missing_project_id_rows": missing_project_id,
        "missing_run_id_rows": missing_run_id,
        "missing_export_id_rows": missing_export_id,
        "missing_frequency_rows": missing_frequency,
        "bad_source_tool_rows": bad_source_tool,
        "missing_path_rows": missing_paths,
        "unknown_kinds": unknown_kinds,
        "bad_upstream_status_rows": bad_upstream_status,
        "bad_upstream_policy_rows": bad_upstream_policy,
        "artifacts": details,
        "checks": checks,
        "status": "ok" if all(checks.values()) else "needs_attention",
        "version_note": (
            "Run after Touchstone and far-field sub-gates so CST notebook rows "
            "cannot mix S-parameter and radiation evidence from different "
            "project/run/export/frequency packages."
        ),
    }


def cst_export_manifest_solver_ready_gate(
    manifest,
    expected_project_id=None,
    expected_run_id=None,
    expected_export_id=None,
    expected_solver_kind=None,
    expected_design_frequency_Hz=None,
    required_file_kinds=("touchstone", "farfield"),
    frequency_rtol=1.0e-12,
):
    """Check a CST export manifest before row-level result gates consume it.

    A manifest is the handoff object between CST automation/export and public
    result validation.  It should freeze identity, solver kind, frequency grid,
    design frequency, file kinds, and the minimal normalization metadata needed
    by Touchstone and far-field gates.
    """

    if not isinstance(manifest, dict):
        raise ValueError("manifest must be a dictionary")

    def _norm(value):
        return str(value or "").strip().lower().replace("-", "_").replace(" ", "_")

    def _first(row, names):
        for name in names:
            if name in row and row[name] is not None:
                return row[name]
        return None

    tolerance = float(frequency_rtol)
    if tolerance < 0.0:
        raise ValueError("frequency_rtol must be non-negative")

    project_id = _first(manifest, ("project_id", "cst_project_id", "model_id"))
    run_id = _first(manifest, ("run_id", "solver_run_id", "study_id"))
    export_id = _first(manifest, ("export_id", "result_set_id", "dataset_id"))
    source_tool = _first(manifest, ("source_tool", "tool", "source"))
    source_tool_norm = _norm(source_tool)
    solver_kind = _norm(_first(manifest, ("solver_kind", "solver", "solver_type")))
    frequency_unit = _norm(_first(manifest, ("frequency_unit", "freq_unit", "grid_unit")))
    design_frequency = _first(manifest, ("design_frequency_Hz", "design_frequency_hz", "frequency_Hz", "frequency_hz"))
    grid = [float(value) for value in manifest.get("frequency_grid_Hz", manifest.get("frequency_grid_hz", ()))]
    files = list(manifest.get("files", manifest.get("artifacts", ())))
    required = tuple(_norm(kind) for kind in required_file_kinds)
    if not required:
        raise ValueError("required_file_kinds must not be empty")

    allowed_solver_kinds = {
        "frequency_domain",
        "hf_frequency_domain",
        "time_domain",
        "integral_equation",
        "eigenmode",
        "asymptotic",
    }
    file_details = []
    file_kind_counts = {}
    missing_file_paths = []
    missing_file_status = []
    bad_file_source = []
    touchstone_missing_metadata = []
    farfield_missing_metadata = []

    for index, item in enumerate(files, start=1):
        if not isinstance(item, dict):
            raise ValueError("each manifest file row must be a dictionary")
        kind = _norm(_first(item, ("kind", "artifact_kind", "type")))
        path = _first(item, ("path", "file", "artifact_path", "table_path"))
        status = _first(item, ("status", "export_status", "validation_status"))
        status_norm = _norm(status)
        file_source = _first(item, ("source_tool", "tool", "source"))
        file_source_norm = _norm(file_source or source_tool)
        if kind:
            file_kind_counts[kind] = file_kind_counts.get(kind, 0) + 1
        if not path:
            missing_file_paths.append(index)
        if status_norm not in {"ok", "pass", "passed", "verified", "exported"}:
            missing_file_status.append({"index": index, "kind": kind, "status": status})
        if file_source_norm not in {"cst", "cst_studio", "cst_studio_suite"}:
            bad_file_source.append({"index": index, "kind": kind, "source_tool": file_source})
        if kind == "touchstone":
            if not item.get("data_format") or item.get("z0_ohm") is None or not item.get("port_order"):
                touchstone_missing_metadata.append(index)
        if kind == "farfield":
            if not item.get("normalization") or not item.get("angle_unit") or not item.get("polarization_basis"):
                farfield_missing_metadata.append(index)
        file_details.append({
            "index": index,
            "kind": kind,
            "path": path,
            "status": status,
            "source_tool": file_source or source_tool,
        })

    grid_strict = bool(grid) and all(a < b for a, b in zip(grid, grid[1:]))
    design_f = None if design_frequency is None else float(design_frequency)
    design_bracketed = False
    design_on_grid = False
    if grid and design_f is not None:
        design_bracketed = min(grid) <= design_f <= max(grid)
        design_on_grid = any(
            abs(freq - design_f) / max(abs(freq), abs(design_f), 1.0) <= tolerance
            for freq in grid
        )

    checks = {
        "project_id_present": bool(project_id),
        "run_id_present": bool(run_id),
        "export_id_present": bool(export_id),
        "source_tool_is_cst": source_tool_norm in {"cst", "cst_studio", "cst_studio_suite"},
        "solver_kind_known": solver_kind in allowed_solver_kinds,
        "frequency_unit_normalized_to_hz": frequency_unit in {"hz", "frequency_hz"},
        "frequency_grid_present": bool(grid),
        "frequency_grid_strictly_increasing": grid_strict,
        "design_frequency_present": design_f is not None,
        "design_frequency_bracketed": design_bracketed,
        "design_frequency_on_grid": design_on_grid,
        "required_file_kinds_present": set(required).issubset(set(file_kind_counts)),
        "file_paths_present": not missing_file_paths,
        "file_status_ok": not missing_file_status,
        "file_source_tool_is_cst": not bad_file_source,
        "touchstone_metadata_present": not touchstone_missing_metadata,
        "farfield_metadata_present": not farfield_missing_metadata,
    }
    if expected_project_id is not None:
        checks["expected_project_id_matches"] = str(project_id) == str(expected_project_id)
    if expected_run_id is not None:
        checks["expected_run_id_matches"] = str(run_id) == str(expected_run_id)
    if expected_export_id is not None:
        checks["expected_export_id_matches"] = str(export_id) == str(expected_export_id)
    if expected_solver_kind is not None:
        checks["expected_solver_kind_matches"] = solver_kind == _norm(expected_solver_kind)
    if expected_design_frequency_Hz is not None:
        expected_f = float(expected_design_frequency_Hz)
        checks["expected_design_frequency_matches"] = (
            design_f is not None
            and abs(design_f - expected_f) / max(abs(design_f), abs(expected_f), 1.0) <= tolerance
        )

    return {
        "policy": "cst_export_manifest_solver_ready_gate",
        "project_id": None if project_id is None else str(project_id),
        "run_id": None if run_id is None else str(run_id),
        "export_id": None if export_id is None else str(export_id),
        "source_tool": source_tool,
        "solver_kind": solver_kind or None,
        "frequency_unit": frequency_unit or None,
        "frequency_grid_Hz": grid,
        "design_frequency_Hz": design_f,
        "required_file_kinds": list(required),
        "present_file_kinds": dict(sorted(file_kind_counts.items())),
        "missing_file_path_rows": missing_file_paths,
        "bad_file_status_rows": missing_file_status,
        "bad_file_source_tool_rows": bad_file_source,
        "touchstone_missing_metadata_rows": touchstone_missing_metadata,
        "farfield_missing_metadata_rows": farfield_missing_metadata,
        "files": file_details,
        "checks": checks,
        "status": "ok" if all(checks.values()) else "needs_attention",
        "version_note": (
            "Run before CST result row gates so solver-ready notebooks receive "
            "one normalized export manifest instead of inferring project, run, "
            "frequency, port, or far-field metadata from filenames."
        ),
    }


def touchstone_row_solver_ready_preflight_gate(
    row,
    data_format="MA",
    z0=50.0,
    return_loss_min_db=None,
    vswr_max=None,
    tol=1.0e-9,
):
    """Bundle public-safe Touchstone row checks before solver reuse.

    This is the compact row-level preflight used for CST, VNA, and open RF
    examples before a Touchstone row is promoted into a solver-ready validation
    record.  It intentionally keeps format, reference impedance, S-matrix
    passivity/reciprocity, and S11 match quality in one payload.
    """

    if not isinstance(row, dict):
        raise ValueError("row must be a mapping with s11/s21 and optional s12/s22")
    fmt = str(data_format or "MA").strip().upper()
    z_ref = float(z0)
    tolerance = float(tol)
    if z_ref <= 0.0:
        raise ValueError("z0 must be > 0")
    if tolerance < 0.0:
        raise ValueError("tol must be non-negative")

    s11 = touchstone_sparameter_to_complex(row.get("s11", [0.0, 0.0]), fmt)
    s21 = touchstone_sparameter_to_complex(row["s21"], fmt)
    s12 = touchstone_sparameter_to_complex(row.get("s12", row["s21"]), fmt)
    s22 = touchstone_sparameter_to_complex(row.get("s22", row.get("s11", [0.0, 0.0])), fmt)
    health = two_port_sparameter_health(s11, s21, s12=s12, s22=s22, tol=tolerance)
    match = one_port_match_quality_gate(s11=s11, tol=tolerance)

    return_loss_limit_ok = True
    if return_loss_min_db is not None:
        return_loss_limit_ok = match["return_loss_db"] >= float(return_loss_min_db)
    vswr_limit_ok = True
    if vswr_max is not None:
        vswr_limit_ok = match["vswr"] <= float(vswr_max)

    checks = {
        "touchstone_format_recorded": fmt in {"RI", "MA", "DB"},
        "reference_impedance_recorded": z_ref > 0.0,
        "sparameter_passivity_ok": health["passive"],
        "sparameter_reciprocity_ok": health["reciprocal"],
        "port_match_passive_reflection_ok": match["checks"]["passive_reflection_ok"],
        "return_loss_limit_ok": return_loss_limit_ok,
        "vswr_limit_ok": vswr_limit_ok,
    }
    return {
        "policy": "touchstone_row_solver_ready_preflight_gate",
        "frequency": row.get("frequency"),
        "data_format": fmt,
        "z0_ohm": z_ref,
        "s11": _complex_row(s11),
        "s21": _complex_row(s21),
        "s12": _complex_row(s12),
        "s22": _complex_row(s22),
        "sparameter_health": health,
        "port_match": match,
        "return_loss_min_db": None if return_loss_min_db is None else float(return_loss_min_db),
        "vswr_max": None if vswr_max is None else float(vswr_max),
        "checks": checks,
        "status": "ok" if all(checks.values()) else "needs_attention",
        "tol": tolerance,
    }


def touchstone_frequency_grid_interpolation_gate(
    frequencies_hz,
    design_frequency_hz,
    max_relative_spacing=0.05,
    require_bracket=True,
):
    """Check that a Touchstone sweep brackets the intended design frequency.

    A good S-parameter row is not enough when the sweep grid is too sparse near
    the frequency reused by a solver, equivalent circuit, or notebook example.
    This gate records the nearest row and the lower/upper bracket used for
    interpolation before a sweep is promoted to validation evidence.
    """

    freqs = [float(value) for value in frequencies_hz]
    design = float(design_frequency_hz)
    spacing_limit = float(max_relative_spacing)
    require = bool(require_bracket)
    if not freqs:
        raise ValueError("frequencies_hz must contain at least one row")
    if design <= 0.0:
        raise ValueError("design_frequency_hz must be > 0")
    if spacing_limit < 0.0:
        raise ValueError("max_relative_spacing must be non-negative")
    if any(freq < 0.0 for freq in freqs):
        raise ValueError("frequencies_hz must be non-negative")
    if any(b <= a for a, b in zip(freqs, freqs[1:])):
        raise ValueError("frequencies_hz must be strictly increasing")

    nearest_index = min(range(len(freqs)), key=lambda idx: abs(freqs[idx] - design))
    nearest_abs_error = abs(freqs[nearest_index] - design)
    nearest_rel_error = nearest_abs_error / design

    exact_tol = max(1.0e-12, abs(design) * 1.0e-12)
    exact_indices = [idx for idx, freq in enumerate(freqs) if abs(freq - design) <= exact_tol]
    if exact_indices:
        lower_index = upper_index = exact_indices[0]
    else:
        lower_candidates = [idx for idx, freq in enumerate(freqs) if freq < design]
        upper_candidates = [idx for idx, freq in enumerate(freqs) if freq > design]
        lower_index = lower_candidates[-1] if lower_candidates else None
        upper_index = upper_candidates[0] if upper_candidates else None

    bracketed = lower_index is not None and upper_index is not None
    if bracketed:
        bracket_gap_hz = freqs[upper_index] - freqs[lower_index]
        bracket_gap_rel = bracket_gap_hz / design
    else:
        bracket_gap_hz = None
        bracket_gap_rel = None

    spacing_ok = bracketed and bracket_gap_rel <= spacing_limit
    if not require and not bracketed:
        spacing_ok = True

    checks = {
        "frequency_grid_strictly_increasing": True,
        "design_frequency_bracketed": bracketed if require else True,
        "design_spacing_ok": spacing_ok,
        "nearest_row_recorded": nearest_index is not None,
    }
    return {
        "policy": "touchstone_frequency_grid_interpolation_gate",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "design_frequency_hz": design,
        "frequency_min_hz": freqs[0],
        "frequency_max_hz": freqs[-1],
        "n_rows": len(freqs),
        "nearest_index": nearest_index,
        "nearest_frequency_hz": freqs[nearest_index],
        "nearest_abs_error_hz": nearest_abs_error,
        "nearest_rel_error": nearest_rel_error,
        "lower_index": lower_index,
        "upper_index": upper_index,
        "lower_frequency_hz": None if lower_index is None else freqs[lower_index],
        "upper_frequency_hz": None if upper_index is None else freqs[upper_index],
        "bracket_gap_hz": bracket_gap_hz,
        "bracket_gap_rel": bracket_gap_rel,
        "max_relative_spacing": spacing_limit,
        "require_bracket": require,
        "checks": checks,
        "notes": [
            "Use the lower/upper bracket, not only the nearest row, before interpolating CST/Touchstone evidence.",
            "A solver-ready row can still fail when the sweep is too coarse around the design frequency.",
        ],
    }


def _touchstone_frequency_unit_scale(unit):
    key = str(unit or "Hz").strip().lower()
    aliases = {
        "hz": 1.0,
        "khz": 1.0e3,
        "mhz": 1.0e6,
        "ghz": 1.0e9,
    }
    if key not in aliases:
        raise ValueError(f"unsupported Touchstone frequency unit: {unit!r}")
    return key, aliases[key]


def touchstone_frequency_unit_normalization_gate(
    frequencies,
    frequency_unit,
    design_frequency,
    design_frequency_unit=None,
    expected_frequency_unit=None,
    expected_design_frequency_hz=None,
    selected_row_index=None,
    max_relative_spacing=0.05,
    require_bracket=True,
    frequency_rtol=1.0e-12,
):
    """Normalize raw Touchstone sweep frequencies to Hz before grid checks.

    CST and Touchstone exports often store compact raw frequency columns such as
    ``1.0`` with an option-line unit like ``GHz``.  This gate keeps the raw
    column, declared unit, normalized Hz grid, design frequency, and optional
    selected design-row index together before interpolation or row reuse.
    """

    raw = [float(value) for value in frequencies]
    unit_key, scale = _touchstone_frequency_unit_scale(frequency_unit)
    design_unit_key, design_scale = _touchstone_frequency_unit_scale(design_frequency_unit or frequency_unit)
    expected_unit_key = None
    if expected_frequency_unit is not None:
        expected_unit_key, _ = _touchstone_frequency_unit_scale(expected_frequency_unit)
    normalized = [value * scale for value in raw]
    design_hz = float(design_frequency) * design_scale
    grid = touchstone_frequency_grid_interpolation_gate(
        normalized,
        design_hz,
        max_relative_spacing=max_relative_spacing,
        require_bracket=require_bracket,
    )
    tolerance = float(frequency_rtol)
    if tolerance < 0.0:
        raise ValueError("frequency_rtol must be non-negative")
    selected_valid = selected_row_index is None
    selected_matches_design = selected_row_index is None
    selected_frequency_hz = None
    if selected_row_index is not None:
        index = int(selected_row_index)
        selected_valid = 0 <= index < len(normalized)
        if selected_valid:
            selected_frequency_hz = normalized[index]
            selected_matches_design = (
                abs(selected_frequency_hz - design_hz)
                <= tolerance * max(1.0, abs(design_hz))
            )
    expected_design_matches = True
    if expected_design_frequency_hz is not None:
        expected = float(expected_design_frequency_hz)
        expected_design_matches = (
            abs(design_hz - expected)
            <= tolerance * max(1.0, abs(expected))
        )
    checks = {
        "frequency_unit_matches_expected": expected_unit_key is None or unit_key == expected_unit_key,
        "design_frequency_matches_expected_hz": expected_design_matches,
        "frequency_grid_contract_ok": grid["status"] == "ok",
        "selected_row_index_valid": selected_valid,
        "selected_row_matches_design_frequency": selected_matches_design,
    }
    return {
        "policy": "touchstone_frequency_unit_normalization_gate",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "frequency_unit": unit_key,
        "design_frequency_unit": design_unit_key,
        "expected_frequency_unit": expected_unit_key,
        "raw_frequencies": raw,
        "frequency_hz": normalized,
        "design_frequency_raw": float(design_frequency),
        "design_frequency_hz": design_hz,
        "expected_design_frequency_hz": None if expected_design_frequency_hz is None else float(expected_design_frequency_hz),
        "selected_row_index": None if selected_row_index is None else int(selected_row_index),
        "selected_frequency_hz": selected_frequency_hz,
        "grid_contract": grid,
        "checks": checks,
        "frequency_rtol": tolerance,
        "notes": [
            "Normalize the raw Touchstone frequency column exactly once before selecting a design row.",
            "A row labelled 1.0 is not solver-ready until the option-line unit, normalized Hz value, and selected row index travel together.",
        ],
    }


def cst_touchstone_solver_ready_manifest_gate(
    artifacts,
    expected_project_id=None,
    expected_run_id=None,
    expected_export_id=None,
    expected_model_input_artifact_id=None,
    expected_model_input_digest=None,
    expected_model_input_path=None,
    expected_parameter_set_artifact_id=None,
    expected_parameter_set_digest=None,
    expected_parameter_set_path=None,
    expected_objective_observable_id=None,
    expected_objective_observable_family=None,
    expected_design_frequency_hz=None,
    expected_network_kind=None,
    expected_port_order=None,
    expected_data_format=None,
    expected_reference_impedance_ohm=None,
    expected_touchstone_option_line_artifact_id=None,
    expected_touchstone_option_line_digest=None,
    expected_reference_plane=None,
    expected_reference_plane_geometry_digest=None,
    expected_port_face_centers_xyz_m=None,
    expected_port_mode_basis=None,
    expected_touchstone_port_mode_basis_schema_id=None,
    expected_incident_wave_convention=None,
    expected_power_balance_basis=None,
    expected_touchstone_export_method=None,
    expected_export_recipe_artifact_id=None,
    expected_export_recipe_digest=None,
    expected_export_recipe_path=None,
    expected_result_tree_path=None,
    expected_result_item_id=None,
    expected_solver_setup_artifact_id=None,
    expected_mesh_setup_artifact_id=None,
    expected_port_definition_artifact_id=None,
    expected_excitation_setup_artifact_id=None,
    expected_frequency_grid_id=None,
    expected_frequency_grid_digest=None,
    expected_interpolation_policy=None,
    expected_touchstone_file_id=None,
    expected_touchstone_observable_id=None,
    expected_touchstone_observable_family=None,
    expected_touchstone_output_artifact_id=None,
    expected_touchstone_output_digest=None,
    expected_touchstone_output_schema_id=None,
    expected_touchstone_output_columns=None,
    expected_touchstone_output_units=None,
    expected_touchstone_convention_schema_id=None,
    expected_touchstone_postprocess_row_convention_schema_id=None,
    expected_created_at_utc=None,
    expected_run_timestamp_utc=None,
    max_created_run_skew_s=None,
    require_touchstone_output_artifact=False,
    require_touchstone_output_schema=False,
    require_touchstone_convention_schema=False,
    require_touchstone_port_mode_basis_schema=False,
    require_touchstone_postprocess_row_convention_schema=False,
    require_export_recipe_artifact=False,
    require_execution_metadata=False,
    expected_renormalized_reference_impedance_ohm=None,
    expected_renormalization_method=None,
    expected_renormalization_artifact_id=None,
    expected_deembedding_method=None,
    expected_deembedding_artifact_id=None,
    expected_deembedding_length_m=None,
    require_model_input_artifact=False,
    require_parameter_set_artifact=False,
    required_kinds=("port_metadata", "frequency_grid", "design_row"),
):
    """Check that CST Touchstone metadata, sweep grid, and design row align.

    This is a Touchstone-only manifest for cases where no far-field export is
    involved.  It keeps port metadata, frequency-grid bracketing, and the
    selected S-parameter design row tied to one CST project/run/export before a
    row is reused by an equivalent-circuit, BEM, or notebook validation.
    """

    rows_in = list(artifacts)
    if not rows_in:
        raise ValueError("artifacts must not be empty")
    if max_created_run_skew_s is not None and float(max_created_run_skew_s) < 0.0:
        raise ValueError("max_created_run_skew_s must be non-negative")

    def _norm(value):
        return str(value or "").strip().lower().replace("-", "_").replace(" ", "_")

    def _first(row, names):
        for name in names:
            if name in row and row[name] is not None:
                return row[name]
        return None

    def _network_kind(value):
        text = str(value or "").strip().upper().replace("-", "_").replace(" ", "_")
        aliases = {
            "SPARAM": "S",
            "SPARAMS": "S",
            "S_PARAM": "S",
            "S_PARAMS": "S",
            "S_PARAMETER": "S",
            "S_PARAMETERS": "S",
            "SCATTERING": "S",
            "SCATTERING_PARAMETER": "S",
            "SCATTERING_PARAMETERS": "S",
        }
        return aliases.get(text, text)

    def _port_order(value):
        if value is None:
            return []
        if isinstance(value, str):
            text = value.strip()
            if not text:
                return []
            text = text.strip("[]()")
            parts = text.replace(";", ",").split(",")
        elif isinstance(value, (list, tuple)):
            parts = value
        else:
            parts = [value]
        return [str(part).strip() for part in parts if str(part).strip()]

    def _data_format(value):
        return str(value or "").strip().upper()

    def _string_list(value):
        if value in (None, ""):
            return []
        if isinstance(value, str):
            text = value.strip()
            if not text:
                return []
            for sep in (";", "\n"):
                text = text.replace(sep, ",")
            return [item.strip() for item in text.split(",") if item.strip()]
        if isinstance(value, dict):
            return [str(item).strip() for item in value.keys() if str(item).strip()]
        return [str(item).strip() for item in value if str(item).strip()]

    def _unit_mapping(value):
        if value in (None, ""):
            return {}
        if isinstance(value, dict):
            return {
                str(key).strip(): str(unit).strip()
                for key, unit in value.items()
                if str(key).strip()
            }
        pairs = {}
        for item in _string_list(value):
            if ":" in item:
                key, unit = item.split(":", 1)
            elif "=" in item:
                key, unit = item.split("=", 1)
            else:
                continue
            key = key.strip()
            if key:
                pairs[key] = unit.strip()
        return pairs

    def _coordinate_tuple(value):
        if value in (None, ""):
            return None
        if isinstance(value, dict):
            lower_keys = ("x", "y", "z")
            upper_keys = ("X", "Y", "Z")
            if all(key in value for key in lower_keys[:2]):
                coords = [value[key] for key in lower_keys if key in value]
            elif all(key in value for key in upper_keys[:2]):
                coords = [value[key] for key in upper_keys if key in value]
            else:
                raise ValueError("coordinate dictionaries must include x/y or X/Y")
        elif isinstance(value, str):
            coords = [item for item in re.split(r"[,;\s]+", value.strip()) if item]
        else:
            coords = list(value)
        if len(coords) not in (2, 3):
            raise ValueError("coordinates must contain two or three values")
        return tuple(float(coord) for coord in coords)

    def _coordinate_sequence(value):
        if value in (None, ""):
            return None
        if isinstance(value, dict):
            items = [value[key] for key in sorted(value)]
        else:
            items = list(value)
        if not items:
            return None
        return tuple(_coordinate_tuple(item) for item in items)

    required = tuple(_norm(kind) for kind in required_kinds)
    if not required:
        raise ValueError("required_kinds must not be empty")

    expected_policies = {
        "port_metadata": {"touchstone_port_metadata_gate"},
        "frequency_grid": {
            "touchstone_frequency_grid_interpolation_gate",
            "cst_touchstone_design_frequency_grid_contract",
        },
        "design_row": {
            "touchstone_row_solver_ready_preflight_gate",
            "two_port_sparameter_health",
        },
        "touchstone_row": {
            "touchstone_row_solver_ready_preflight_gate",
            "two_port_sparameter_health",
        },
    }

    details = []
    kind_counts = {}
    project_ids = []
    run_ids = []
    export_ids = []
    model_input_artifact_ids = []
    model_input_digests = []
    model_input_paths = []
    parameter_set_artifact_ids = []
    parameter_set_digests = []
    parameter_set_paths = []
    objective_observable_ids = []
    objective_observable_families = []
    design_frequencies = []
    missing_project_id = []
    missing_run_id = []
    missing_export_id = []
    missing_model_input_artifact_id = []
    missing_model_input_digest = []
    missing_model_input_path = []
    missing_parameter_set_artifact_id = []
    missing_parameter_set_digest = []
    missing_parameter_set_path = []
    missing_objective_observable_id = []
    missing_objective_observable_family = []
    missing_frequency = []
    missing_paths = []
    bad_source_tool = []
    unknown_kinds = []
    bad_upstream_status = []
    bad_upstream_policy = []
    missing_port_metadata = []
    unbracketed_grid = []
    nonpassive_rows = []
    nonreciprocal_rows = []
    nonpassive_touchstone_rows = []
    nonreciprocal_touchstone_rows = []
    missing_frequency_grid_row_count = []
    missing_design_row_index = []
    design_row_index_out_of_range = []
    missing_touchstone_row_index = []
    missing_touchstone_row_frequency = []
    touchstone_row_index_mismatch = []
    touchstone_row_frequency_mismatch = []
    frequency_grid_row_counts = []
    design_row_indices = []
    touchstone_row_indices = []
    touchstone_row_frequencies = []
    touchstone_row_index_records = []
    touchstone_row_frequency_records = []
    missing_network_kind = []
    missing_port_order = []
    missing_data_format = []
    missing_reference_impedance = []
    missing_touchstone_option_line_artifact_id = []
    missing_touchstone_option_line_digest = []
    missing_reference_plane = []
    missing_reference_plane_geometry_digest = []
    missing_port_face_centers_xyz_m = []
    missing_port_mode_basis = []
    missing_touchstone_port_mode_basis_schema_id = []
    missing_incident_wave_convention = []
    missing_power_balance_basis = []
    missing_touchstone_export_method = []
    missing_export_recipe_artifact_id = []
    missing_export_recipe_digest = []
    missing_export_recipe_path = []
    missing_result_tree_path = []
    missing_result_item_id = []
    missing_solver_setup_artifact_id = []
    missing_mesh_setup_artifact_id = []
    missing_port_definition_artifact_id = []
    missing_excitation_setup_artifact_id = []
    missing_frequency_grid_id = []
    missing_frequency_grid_digest = []
    missing_interpolation_policy = []
    missing_selected_frequency = []
    selected_frequency_mismatch = []
    network_kinds = []
    port_orders = []
    data_formats = []
    reference_impedances = []
    touchstone_option_line_artifact_ids = []
    touchstone_option_line_digests = []
    reference_planes = []
    reference_plane_geometry_digests = []
    port_face_centers_xyz_m = []
    port_mode_bases = []
    touchstone_port_mode_basis_schema_ids = []
    incident_wave_conventions = []
    power_balance_bases = []
    touchstone_export_methods = []
    export_recipe_artifact_ids = []
    export_recipe_digests = []
    export_recipe_paths = []
    result_tree_paths = []
    result_item_ids = []
    solver_setup_artifact_ids = []
    mesh_setup_artifact_ids = []
    port_definition_artifact_ids = []
    excitation_setup_artifact_ids = []
    frequency_grid_ids = []
    frequency_grid_digests = []
    interpolation_policies = []
    selected_frequencies = []
    touchstone_file_ids = []
    missing_touchstone_file_id = []
    touchstone_observable_ids = []
    touchstone_observable_families = []
    missing_touchstone_observable_id = []
    missing_touchstone_observable_family = []
    touchstone_output_artifact_ids = []
    touchstone_output_digests = []
    touchstone_output_paths = []
    touchstone_output_schema_ids = []
    touchstone_convention_schema_ids = []
    touchstone_postprocess_row_convention_schema_ids = []
    touchstone_output_column_sets = []
    touchstone_output_unit_maps = []
    missing_touchstone_output_artifact_id = []
    missing_touchstone_output_digest = []
    missing_touchstone_output_path = []
    missing_touchstone_output_schema_id = []
    missing_touchstone_convention_schema_id = []
    missing_touchstone_postprocess_row_convention_schema_id = []
    missing_touchstone_output_columns = []
    missing_touchstone_output_units = []
    created_at_utc_values = []
    run_timestamp_utc_values = []
    bad_created_at_utc = []
    bad_run_timestamp_utc = []
    missing_created_at_utc = []
    missing_run_timestamp_utc = []
    created_run_timestamp_skews_s = []
    renormalized_reference_impedances = []
    renormalization_methods = []
    renormalization_artifact_ids = []
    missing_renormalization_method = []
    missing_renormalization_artifact_id = []
    missing_renormalized_reference_impedance = []
    deembedding_methods = []
    deembedding_artifact_ids = []
    deembedding_lengths_m = []
    missing_deembedding_method = []
    missing_deembedding_artifact_id = []
    missing_deembedding_length = []

    for index, row in enumerate(rows_in, start=1):
        if not isinstance(row, dict):
            raise ValueError("each artifact must be a dictionary")
        kind = _norm(_first(row, ("kind", "artifact_kind", "type")))
        project_id = _first(row, ("project_id", "cst_project_id", "model_id"))
        run_id = _first(row, ("run_id", "solver_run_id", "simulation_id"))
        export_id = _first(row, ("export_id", "result_id", "dataset_id"))
        model_input_artifact_id = _first(row, (
            "model_input_artifact_id",
            "cst_project_artifact_id",
            "project_artifact_id",
            "input_project_artifact_id",
            "model_artifact_id",
            "cst_model_artifact_id",
        ))
        model_input_digest = _first(row, (
            "model_input_digest",
            "cst_project_digest",
            "project_digest",
            "input_project_digest",
            "model_digest",
            "cst_model_digest",
        ))
        model_input_path = _first(row, (
            "model_input_path",
            "cst_project_path",
            "project_path",
            "input_project_path",
            "model_path",
            "cst_model_path",
        ))
        parameter_set_artifact_id = _first(row, (
            "parameter_set_artifact_id",
            "design_parameter_set_artifact_id",
            "optimization_parameter_set_artifact_id",
            "tuning_parameter_set_artifact_id",
            "parameter_artifact_id",
        ))
        parameter_set_digest = _first(row, (
            "parameter_set_digest",
            "parameter_set_sha256",
            "design_parameter_set_digest",
            "optimization_parameter_set_digest",
            "tuning_parameter_set_digest",
            "parameter_digest",
        ))
        parameter_set_path = _first(row, (
            "parameter_set_path",
            "parameter_set_file",
            "design_parameter_set_path",
            "optimization_parameter_set_path",
            "tuning_parameter_set_path",
            "parameter_path",
        ))
        objective_observable_id = _first(row, (
            "objective_observable_id",
            "optimization_objective_observable_id",
            "design_objective_observable_id",
            "target_observable_id",
            "observable_objective_id",
        ))
        objective_observable_family = _first(row, (
            "objective_observable_family",
            "optimization_objective_observable_family",
            "design_objective_observable_family",
            "target_observable_family",
            "observable_objective_family",
        ))
        touchstone_export_method = _first(row, (
            "touchstone_export_method",
            "sparameter_export_method",
            "cst_export_method",
            "result_export_method",
            "export_method",
        ))
        export_recipe_artifact_id = _first(row, (
            "export_recipe_artifact_id",
            "touchstone_export_recipe_artifact_id",
            "export_script_artifact_id",
            "export_macro_artifact_id",
            "result_export_recipe_artifact_id",
            "postprocess_recipe_artifact_id",
            "vba_macro_artifact_id",
        ))
        export_recipe_digest = _first(row, (
            "export_recipe_digest",
            "export_recipe_sha256",
            "touchstone_export_recipe_digest",
            "touchstone_export_recipe_sha256",
            "export_script_digest",
            "export_script_sha256",
            "export_macro_digest",
            "export_macro_sha256",
            "result_export_recipe_digest",
            "postprocess_recipe_digest",
            "vba_macro_digest",
        ))
        export_recipe_path = _first(row, (
            "export_recipe_path",
            "export_recipe_file",
            "touchstone_export_recipe_path",
            "export_script_path",
            "export_macro_path",
            "result_export_recipe_path",
            "postprocess_recipe_path",
            "vba_macro_path",
        ))
        result_tree_path = _first(row, (
            "result_tree_path",
            "touchstone_result_tree_path",
            "cst_result_tree_path",
            "result_item_path",
            "monitor_path",
            "result_path",
        ))
        result_item_id = _first(row, (
            "result_item_id",
            "touchstone_result_item_id",
            "cst_result_item_id",
            "monitor_id",
            "result_tree_item_id",
            "result_id",
        ))
        solver_setup_artifact_id = _first(row, (
            "solver_setup_artifact_id",
            "solver_configuration_artifact_id",
            "simulation_setup_artifact_id",
            "study_solver_setup_artifact_id",
            "cst_solver_setup_artifact_id",
        ))
        mesh_setup_artifact_id = _first(row, (
            "mesh_setup_artifact_id",
            "mesh_artifact_id",
            "mesh_revision_id",
            "adaptive_mesh_artifact_id",
            "mesh_configuration_artifact_id",
            "cst_mesh_setup_artifact_id",
        ))
        port_definition_artifact_id = _first(row, (
            "port_definition_artifact_id",
            "port_setup_artifact_id",
            "port_metadata_artifact_id",
            "discrete_port_artifact_id",
            "waveguide_port_artifact_id",
            "cst_port_definition_artifact_id",
        ))
        excitation_setup_artifact_id = _first(row, (
            "excitation_setup_artifact_id",
            "port_excitation_artifact_id",
            "source_setup_artifact_id",
            "stimulus_artifact_id",
            "cst_excitation_setup_artifact_id",
        ))
        source_tool = _first(row, ("source_tool", "tool", "source"))
        source_tool_norm = _norm(source_tool)
        path = _first(row, ("path", "file", "artifact_path", "table_path"))
        gate_policy = _first(row, ("gate_policy", "policy", "validator"))
        gate_policy_norm = _norm(gate_policy)
        status = _first(row, ("status", "gate_status", "validation_status"))
        status_norm = _norm(status)
        pass_flag = bool(row.get("pass", False))
        frequency = _first(row, ("design_frequency_hz", "frequency_hz", "freq_hz"))
        row_count = _first(row, ("row_count", "n_rows", "sweep_row_count"))
        selected_row_index = _first(row, ("selected_row_index", "design_row_index", "row_index"))
        network_kind = _network_kind(_first(row, ("network_kind", "network_parameter", "parameter_kind")))
        port_order = _port_order(_first(row, ("port_order", "ports", "port_names")))
        data_format = _data_format(_first(row, ("data_format", "format", "touchstone_format")))
        reference_impedance = _first(row, ("reference_impedance_ohm", "z0_ohm", "z0"))
        touchstone_option_line_artifact_id = _first(row, (
            "touchstone_option_line_artifact_id",
            "option_line_artifact_id",
            "touchstone_header_artifact_id",
            "format_option_line_artifact_id",
        ))
        touchstone_option_line_digest = _first(row, (
            "touchstone_option_line_digest",
            "option_line_digest",
            "touchstone_header_digest",
            "format_option_line_digest",
            "option_line_sha256",
        ))
        reference_plane = _first(row, ("reference_plane", "port_reference_plane", "deembed_reference_plane"))
        reference_plane_geometry_digest = _first(row, (
            "reference_plane_geometry_digest",
            "port_reference_plane_geometry_digest",
            "port_face_geometry_digest",
            "port_geometry_digest",
        ))
        port_face_centers = _first(row, (
            "port_face_centers_xyz_m",
            "port_centers_xyz_m",
            "port_reference_points_xyz_m",
            "reference_plane_points_xyz_m",
        ))
        port_mode_basis = _first(row, ("port_mode_basis", "mode_basis", "wave_mode_basis"))
        touchstone_port_mode_basis_schema_id = _first(row, (
            "touchstone_port_mode_basis_schema_id",
            "port_mode_basis_schema_id",
            "mode_basis_schema_id",
            "wave_mode_basis_schema_id",
            "sparameter_port_mode_basis_schema_id",
        ))
        incident_wave_convention = _first(row, (
            "incident_wave_convention",
            "incident_power_convention",
            "incident_wave_basis",
            "sparameter_incident_convention",
            "excitation_convention",
        ))
        power_balance_basis = _first(row, (
            "power_balance_basis",
            "power_wave_basis",
            "power_normalization",
            "power_ledger_basis",
        ))
        frequency_grid_id = _first(row, ("frequency_grid_id", "grid_id", "sweep_grid_id", "frequency_axis_id"))
        frequency_grid_digest = _first(row, (
            "frequency_grid_digest",
            "frequency_grid_sha256",
            "frequency_axis_digest",
            "frequency_axis_sha256",
            "sweep_grid_digest",
            "sweep_grid_sha256",
        ))
        interpolation_policy = _first(row, (
            "interpolation_policy",
            "frequency_interpolation_policy",
            "row_selection_policy",
        ))
        selected_frequency = _first(row, (
            "selected_frequency_hz",
            "selected_frequency_Hz",
            "selected_row_frequency_hz",
            "design_row_frequency_hz",
        ))
        touchstone_row_index = _first(row, (
            "touchstone_row_index",
            "row_index",
            "selected_row_index",
            "design_row_index",
        ))
        touchstone_row_frequency = _first(row, (
            "touchstone_row_frequency_hz",
            "row_frequency_hz",
            "frequency_hz",
            "selected_frequency_hz",
            "design_row_frequency_hz",
        ))
        touchstone_file_id = _first(row, (
            "touchstone_file_id",
            "raw_touchstone_file_id",
            "snp_file_id",
            "touchstone_sha256",
            "raw_file_sha256",
        ))
        touchstone_observable_id = _first(row, (
            "touchstone_observable_id",
            "sparameter_observable_id",
            "network_observable_id",
            "observable_id",
        ))
        touchstone_observable_family = _first(row, (
            "touchstone_observable_family",
            "sparameter_observable_family",
            "network_observable_family",
            "observable_family",
        ))
        touchstone_output_artifact_id = _first(row, (
            "touchstone_output_artifact_id",
            "output_artifact_id",
            "selected_row_output_artifact_id",
            "postprocess_output_artifact_id",
            "table_artifact_id",
            "export_output_artifact_id",
        ))
        touchstone_output_digest = _first(row, (
            "touchstone_output_digest",
            "touchstone_output_sha256",
            "output_digest",
            "output_sha256",
            "selected_row_output_digest",
            "postprocess_output_digest",
            "table_sha256",
            "export_output_digest",
        ))
        touchstone_output_path = _first(row, (
            "touchstone_output_path",
            "output_path",
            "selected_row_output_path",
            "postprocess_output_path",
            "table_path",
            "export_output_path",
        ))
        touchstone_output_schema_id = _first(row, (
            "touchstone_output_schema_id",
            "touchstone_table_schema_id",
            "sparameter_output_schema_id",
            "postprocess_output_schema_id",
            "output_schema_id",
            "table_schema_id",
        ))
        touchstone_convention_schema_id = _first(row, (
            "touchstone_convention_schema_id",
            "touchstone_physics_convention_schema_id",
            "network_convention_schema_id",
            "rf_convention_schema_id",
            "physics_convention_schema_id",
        ))
        touchstone_postprocess_row_convention_schema_id = _first(row, (
            "touchstone_postprocess_row_convention_schema_id",
            "postprocess_row_convention_schema_id",
            "touchstone_row_convention_schema_id",
            "sparameter_row_convention_schema_id",
            "postprocess_convention_schema_id",
        ))
        touchstone_output_columns = _string_list(_first(row, (
            "touchstone_output_columns",
            "touchstone_table_columns",
            "sparameter_output_columns",
            "postprocess_output_columns",
            "output_columns",
            "table_columns",
            "columns",
        )))
        touchstone_output_units = _unit_mapping(_first(row, (
            "touchstone_output_units",
            "touchstone_table_units",
            "sparameter_output_units",
            "postprocess_output_units",
            "output_units",
            "table_units",
            "column_units",
            "units",
        )))
        created_at_utc = _first(row, (
            "created_at_utc",
            "artifact_created_at_utc",
            "created_at",
        ))
        run_timestamp_utc = _first(row, (
            "run_timestamp_utc",
            "executed_at_utc",
            "run_date_utc",
            "date_utc",
            "run_date",
        ))
        renormalized_reference_impedance = _first(row, (
            "renormalized_reference_impedance_ohm",
            "renormalized_z0_ohm",
            "analysis_reference_impedance_ohm",
            "postprocess_reference_impedance_ohm",
        ))
        renormalization_method = _first(row, (
            "renormalization_method",
            "renormalization_policy",
            "sparameter_renormalization",
            "z0_renormalization_method",
        ))
        renormalization_artifact_id = _first(row, (
            "renormalization_artifact_id",
            "renormalized_touchstone_file_id",
            "renormalization_source_artifact_id",
            "renormalized_snp_file_id",
        ))
        deembedding_method = _first(row, (
            "deembedding_method",
            "deembed_method",
            "port_extension_method",
            "reference_plane_shift_method",
        ))
        deembedding_artifact_id = _first(row, (
            "deembedding_artifact_id",
            "deembed_artifact_id",
            "port_extension_artifact_id",
            "reference_plane_shift_artifact_id",
            "calibration_artifact_id",
        ))
        deembedding_length = _first(row, (
            "deembedding_length_m",
            "deembed_length_m",
            "port_extension_length_m",
            "reference_plane_shift_m",
            "reference_plane_offset_m",
        ))
        deembedding_length_value = None

        if kind:
            kind_counts[kind] = kind_counts.get(kind, 0) + 1
        else:
            unknown_kinds.append(index)
        if kind and kind not in set(required) | set(expected_policies):
            unknown_kinds.append(kind)
        if not project_id:
            missing_project_id.append(index)
        else:
            project_ids.append(str(project_id))
        if not run_id:
            missing_run_id.append(index)
        else:
            run_ids.append(str(run_id))
        if not export_id:
            missing_export_id.append(index)
        else:
            export_ids.append(str(export_id))
        model_input_expected = (
            require_model_input_artifact
            or expected_model_input_artifact_id is not None
            or expected_model_input_digest is not None
            or expected_model_input_path is not None
        )
        if model_input_artifact_id is None or not str(model_input_artifact_id).strip():
            if model_input_expected:
                missing_model_input_artifact_id.append(index)
        else:
            model_input_artifact_ids.append(str(model_input_artifact_id).strip())
        if model_input_digest is None or not str(model_input_digest).strip():
            if model_input_expected:
                missing_model_input_digest.append(index)
        else:
            model_input_digests.append(str(model_input_digest).strip())
        if model_input_path is None or not str(model_input_path).strip():
            if model_input_expected:
                missing_model_input_path.append(index)
        else:
            model_input_paths.append(str(model_input_path).strip())
        parameter_set_expected = (
            require_parameter_set_artifact
            or expected_parameter_set_artifact_id is not None
            or expected_parameter_set_digest is not None
            or expected_parameter_set_path is not None
        )
        if parameter_set_artifact_id is None or not str(parameter_set_artifact_id).strip():
            if parameter_set_expected:
                missing_parameter_set_artifact_id.append(index)
        else:
            parameter_set_artifact_ids.append(str(parameter_set_artifact_id).strip())
        if parameter_set_digest is None or not str(parameter_set_digest).strip():
            if parameter_set_expected:
                missing_parameter_set_digest.append(index)
        else:
            parameter_set_digests.append(str(parameter_set_digest).strip())
        if parameter_set_path is None or not str(parameter_set_path).strip():
            if parameter_set_expected:
                missing_parameter_set_path.append(index)
        else:
            parameter_set_paths.append(str(parameter_set_path).strip())
        if objective_observable_id is None or not str(objective_observable_id).strip():
            if expected_objective_observable_id is not None:
                missing_objective_observable_id.append(index)
        else:
            objective_observable_ids.append(str(objective_observable_id).strip())
        if objective_observable_family is None or not str(objective_observable_family).strip():
            if expected_objective_observable_family is not None:
                missing_objective_observable_family.append(index)
        else:
            objective_observable_families.append(_norm(objective_observable_family))
        if source_tool_norm not in {"cst", "cst_studio", "cst_studio_suite"}:
            bad_source_tool.append({"index": index, "kind": kind, "source_tool": source_tool})
        if not path:
            missing_paths.append(index)
        if frequency is None:
            missing_frequency.append(index)
        else:
            design_frequencies.append(float(frequency))
        if not network_kind:
            missing_network_kind.append(index)
        else:
            network_kinds.append(network_kind)
        if not port_order:
            missing_port_order.append(index)
        else:
            port_orders.append(tuple(port_order))
        if not data_format:
            missing_data_format.append(index)
        else:
            data_formats.append(data_format)
        if reference_impedance is None or float(reference_impedance) <= 0.0:
            missing_reference_impedance.append(index)
        else:
            reference_impedances.append(float(reference_impedance))
        if touchstone_option_line_artifact_id is None or not str(touchstone_option_line_artifact_id).strip():
            if expected_touchstone_option_line_artifact_id is not None:
                missing_touchstone_option_line_artifact_id.append(index)
        else:
            touchstone_option_line_artifact_ids.append(str(touchstone_option_line_artifact_id).strip())
        if touchstone_option_line_digest is None or not str(touchstone_option_line_digest).strip():
            if expected_touchstone_option_line_digest is not None:
                missing_touchstone_option_line_digest.append(index)
        else:
            touchstone_option_line_digests.append(str(touchstone_option_line_digest).strip())
        if reference_plane is None or not str(reference_plane).strip():
            if expected_reference_plane is not None:
                missing_reference_plane.append(index)
        else:
            reference_planes.append(str(reference_plane).strip())
        if reference_plane_geometry_digest is None or not str(reference_plane_geometry_digest).strip():
            if expected_reference_plane_geometry_digest is not None:
                missing_reference_plane_geometry_digest.append(index)
        else:
            reference_plane_geometry_digests.append(str(reference_plane_geometry_digest).strip())
        if port_face_centers is None:
            port_face_centers_value = None
            if expected_port_face_centers_xyz_m is not None:
                missing_port_face_centers_xyz_m.append(index)
        else:
            port_face_centers_value = _coordinate_sequence(port_face_centers)
            if port_face_centers_value is None:
                if expected_port_face_centers_xyz_m is not None:
                    missing_port_face_centers_xyz_m.append(index)
            else:
                port_face_centers_xyz_m.append(port_face_centers_value)
        if port_mode_basis is None or not str(port_mode_basis).strip():
            if expected_port_mode_basis is not None:
                missing_port_mode_basis.append(index)
        else:
            port_mode_bases.append(str(port_mode_basis).strip())
        if incident_wave_convention is None or not str(incident_wave_convention).strip():
            if expected_incident_wave_convention is not None:
                missing_incident_wave_convention.append(index)
        else:
            incident_wave_conventions.append(str(incident_wave_convention).strip())
        if power_balance_basis is None or not str(power_balance_basis).strip():
            if expected_power_balance_basis is not None:
                missing_power_balance_basis.append(index)
        else:
            power_balance_bases.append(str(power_balance_basis).strip())
        if touchstone_export_method is None or not str(touchstone_export_method).strip():
            if expected_touchstone_export_method is not None:
                missing_touchstone_export_method.append(index)
        else:
            touchstone_export_methods.append(str(touchstone_export_method).strip())
        export_recipe_expected = (
            require_export_recipe_artifact
            or expected_export_recipe_artifact_id is not None
            or expected_export_recipe_digest is not None
            or expected_export_recipe_path is not None
        )
        if export_recipe_artifact_id is None or not str(export_recipe_artifact_id).strip():
            if export_recipe_expected:
                missing_export_recipe_artifact_id.append(index)
        else:
            export_recipe_artifact_ids.append(str(export_recipe_artifact_id).strip())
        if export_recipe_digest is None or not str(export_recipe_digest).strip():
            if export_recipe_expected:
                missing_export_recipe_digest.append(index)
        else:
            export_recipe_digests.append(str(export_recipe_digest).strip())
        if export_recipe_path is None or not str(export_recipe_path).strip():
            if export_recipe_expected:
                missing_export_recipe_path.append(index)
        else:
            export_recipe_paths.append(str(export_recipe_path).strip())
        if result_tree_path is None or not str(result_tree_path).strip():
            if expected_result_tree_path is not None:
                missing_result_tree_path.append(index)
        else:
            result_tree_paths.append(str(result_tree_path).strip())
        if result_item_id is None or not str(result_item_id).strip():
            if expected_result_item_id is not None:
                missing_result_item_id.append(index)
        else:
            result_item_ids.append(str(result_item_id).strip())
        if solver_setup_artifact_id is None or not str(solver_setup_artifact_id).strip():
            if expected_solver_setup_artifact_id is not None:
                missing_solver_setup_artifact_id.append(index)
        else:
            solver_setup_artifact_ids.append(str(solver_setup_artifact_id).strip())
        if mesh_setup_artifact_id is None or not str(mesh_setup_artifact_id).strip():
            if expected_mesh_setup_artifact_id is not None:
                missing_mesh_setup_artifact_id.append(index)
        else:
            mesh_setup_artifact_ids.append(str(mesh_setup_artifact_id).strip())
        if port_definition_artifact_id is None or not str(port_definition_artifact_id).strip():
            if expected_port_definition_artifact_id is not None:
                missing_port_definition_artifact_id.append(index)
        else:
            port_definition_artifact_ids.append(str(port_definition_artifact_id).strip())
        if excitation_setup_artifact_id is None or not str(excitation_setup_artifact_id).strip():
            if expected_excitation_setup_artifact_id is not None:
                missing_excitation_setup_artifact_id.append(index)
        else:
            excitation_setup_artifact_ids.append(str(excitation_setup_artifact_id).strip())
        if frequency_grid_id is None or not str(frequency_grid_id).strip():
            if expected_frequency_grid_id is not None:
                missing_frequency_grid_id.append(index)
        else:
            frequency_grid_ids.append(str(frequency_grid_id).strip())
        if frequency_grid_digest is None or not str(frequency_grid_digest).strip():
            if expected_frequency_grid_digest is not None:
                missing_frequency_grid_digest.append(index)
        else:
            frequency_grid_digests.append(str(frequency_grid_digest).strip())
        if interpolation_policy is None or not str(interpolation_policy).strip():
            if expected_interpolation_policy is not None:
                missing_interpolation_policy.append(index)
        else:
            interpolation_policies.append(str(interpolation_policy).strip())
        if touchstone_file_id is None or not str(touchstone_file_id).strip():
            if expected_touchstone_file_id is not None:
                missing_touchstone_file_id.append(index)
        else:
            touchstone_file_ids.append(str(touchstone_file_id).strip())
        if touchstone_observable_id is None or not str(touchstone_observable_id).strip():
            if expected_touchstone_observable_id is not None:
                missing_touchstone_observable_id.append(index)
        else:
            touchstone_observable_ids.append(str(touchstone_observable_id).strip())
        if touchstone_observable_family is None or not str(touchstone_observable_family).strip():
            if expected_touchstone_observable_family is not None:
                missing_touchstone_observable_family.append(index)
        else:
            touchstone_observable_families.append(_norm(touchstone_observable_family))
        touchstone_output_expected = (
            require_touchstone_output_artifact
            or expected_touchstone_output_artifact_id is not None
            or expected_touchstone_output_digest is not None
        )
        touchstone_output_schema_expected = (
            require_touchstone_output_schema
            or expected_touchstone_output_schema_id is not None
            or bool(_string_list(expected_touchstone_output_columns))
            or bool(_unit_mapping(expected_touchstone_output_units))
        )
        touchstone_convention_schema_expected = (
            require_touchstone_convention_schema
            or expected_touchstone_convention_schema_id is not None
        )
        touchstone_port_mode_basis_schema_expected = (
            require_touchstone_port_mode_basis_schema
            or expected_touchstone_port_mode_basis_schema_id is not None
        )
        touchstone_postprocess_row_convention_schema_expected = (
            require_touchstone_postprocess_row_convention_schema
            or expected_touchstone_postprocess_row_convention_schema_id is not None
        )
        if touchstone_output_artifact_id is None or not str(touchstone_output_artifact_id).strip():
            if touchstone_output_expected:
                missing_touchstone_output_artifact_id.append(index)
        else:
            touchstone_output_artifact_ids.append(str(touchstone_output_artifact_id).strip())
        if touchstone_output_digest is None or not str(touchstone_output_digest).strip():
            if touchstone_output_expected:
                missing_touchstone_output_digest.append(index)
        else:
            touchstone_output_digests.append(str(touchstone_output_digest).strip())
        if touchstone_output_path is None or not str(touchstone_output_path).strip():
            if touchstone_output_expected:
                missing_touchstone_output_path.append(index)
        else:
            touchstone_output_paths.append(str(touchstone_output_path).strip())
        if touchstone_output_schema_id is None or not str(touchstone_output_schema_id).strip():
            if touchstone_output_schema_expected:
                missing_touchstone_output_schema_id.append(index)
        else:
            touchstone_output_schema_ids.append(str(touchstone_output_schema_id).strip())
        if touchstone_convention_schema_id is None or not str(touchstone_convention_schema_id).strip():
            if touchstone_convention_schema_expected:
                missing_touchstone_convention_schema_id.append(index)
        else:
            touchstone_convention_schema_ids.append(str(touchstone_convention_schema_id).strip())
        if (
            touchstone_port_mode_basis_schema_id is None
            or not str(touchstone_port_mode_basis_schema_id).strip()
        ):
            if touchstone_port_mode_basis_schema_expected:
                missing_touchstone_port_mode_basis_schema_id.append(index)
        else:
            touchstone_port_mode_basis_schema_ids.append(
                str(touchstone_port_mode_basis_schema_id).strip()
            )
        if (
            touchstone_postprocess_row_convention_schema_id is None
            or not str(touchstone_postprocess_row_convention_schema_id).strip()
        ):
            if touchstone_postprocess_row_convention_schema_expected:
                missing_touchstone_postprocess_row_convention_schema_id.append(index)
        else:
            touchstone_postprocess_row_convention_schema_ids.append(
                str(touchstone_postprocess_row_convention_schema_id).strip()
            )
        if not touchstone_output_columns:
            if touchstone_output_schema_expected:
                missing_touchstone_output_columns.append(index)
        else:
            touchstone_output_column_sets.append(tuple(touchstone_output_columns))
        if not touchstone_output_units:
            if touchstone_output_schema_expected:
                missing_touchstone_output_units.append(index)
        else:
            touchstone_output_unit_maps.append(tuple(sorted(touchstone_output_units.items())))
        execution_expected = (
            require_execution_metadata
            or expected_created_at_utc is not None
            or expected_run_timestamp_utc is not None
            or max_created_run_skew_s is not None
        )
        created_dt = None
        run_dt = None
        if created_at_utc is None or not str(created_at_utc).strip():
            if execution_expected:
                missing_created_at_utc.append(index)
        else:
            created_text = str(created_at_utc).strip()
            created_at_utc_values.append(created_text)
            created_dt = _parse_utc_like_datetime(created_text)
            if created_dt is None:
                bad_created_at_utc.append(index)
        if run_timestamp_utc is None or not str(run_timestamp_utc).strip():
            if execution_expected:
                missing_run_timestamp_utc.append(index)
        else:
            run_text = str(run_timestamp_utc).strip()
            run_timestamp_utc_values.append(run_text)
            run_dt = _parse_utc_like_datetime(run_text)
            if run_dt is None:
                bad_run_timestamp_utc.append(index)
        if created_dt is not None and run_dt is not None:
            created_run_timestamp_skews_s.append(abs((created_dt - run_dt).total_seconds()))
        if renormalized_reference_impedance is None or not str(renormalized_reference_impedance).strip():
            if expected_renormalized_reference_impedance_ohm is not None:
                missing_renormalized_reference_impedance.append(index)
        else:
            renormalized_reference_impedances.append(float(renormalized_reference_impedance))
            if renormalization_method is None or not str(renormalization_method).strip():
                missing_renormalization_method.append(index)
            else:
                renormalization_methods.append(str(renormalization_method).strip())
            if renormalization_artifact_id is None or not str(renormalization_artifact_id).strip():
                missing_renormalization_artifact_id.append(index)
            else:
                renormalization_artifact_ids.append(str(renormalization_artifact_id).strip())
        if expected_renormalization_method is not None and (
            renormalization_method is None or not str(renormalization_method).strip()
        ):
            missing_renormalization_method.append(index)
        if expected_renormalization_artifact_id is not None and (
            renormalization_artifact_id is None or not str(renormalization_artifact_id).strip()
        ):
            missing_renormalization_artifact_id.append(index)
        has_deembedding_metadata = any(
            item is not None and str(item).strip()
            for item in (deembedding_method, deembedding_artifact_id, deembedding_length)
        )
        deembedding_expected = any(
            item is not None
            for item in (
                expected_deembedding_method,
                expected_deembedding_artifact_id,
                expected_deembedding_length_m,
            )
        )
        if deembedding_method is None or not str(deembedding_method).strip():
            if has_deembedding_metadata or deembedding_expected:
                missing_deembedding_method.append(index)
        else:
            deembedding_methods.append(str(deembedding_method).strip())
        if deembedding_artifact_id is None or not str(deembedding_artifact_id).strip():
            if has_deembedding_metadata or deembedding_expected:
                missing_deembedding_artifact_id.append(index)
        else:
            deembedding_artifact_ids.append(str(deembedding_artifact_id).strip())
        if deembedding_length is None or not str(deembedding_length).strip():
            if has_deembedding_metadata or deembedding_expected:
                missing_deembedding_length.append(index)
        else:
            try:
                deembedding_length_value = round(float(deembedding_length), 15)
                deembedding_lengths_m.append(deembedding_length_value)
            except (TypeError, ValueError):
                missing_deembedding_length.append(index)
        if not (pass_flag or status_norm in {"ok", "pass", "passed", "verified"}):
            bad_upstream_status.append({"index": index, "kind": kind, "status": status})
        if kind in expected_policies and gate_policy_norm not in expected_policies[kind]:
            bad_upstream_policy.append({
                "index": index,
                "kind": kind,
                "gate_policy": gate_policy,
                "expected": sorted(expected_policies[kind]),
            })
        if kind == "port_metadata" and not (
            row.get("data_format")
            and _first(row, ("reference_impedance_ohm", "z0_ohm", "z0")) is not None
            and row.get("port_order")
        ):
            missing_port_metadata.append(index)
        if kind == "frequency_grid" and row.get("design_frequency_bracketed") is not True:
            unbracketed_grid.append(index)
        if kind == "frequency_grid":
            if row_count is None or int(row_count) <= 0:
                missing_frequency_grid_row_count.append(index)
            else:
                frequency_grid_row_counts.append(int(row_count))
        if kind == "design_row" and row.get("sparameter_passivity_ok") is not True:
            nonpassive_rows.append(index)
        if kind == "design_row" and row.get("sparameter_reciprocity_ok") is not True:
            nonreciprocal_rows.append(index)
        if kind == "touchstone_row" and row.get("sparameter_passivity_ok") is not True:
            nonpassive_touchstone_rows.append(index)
        if kind == "touchstone_row" and row.get("sparameter_reciprocity_ok") is not True:
            nonreciprocal_touchstone_rows.append(index)
        if kind == "design_row":
            if selected_row_index is None:
                missing_design_row_index.append(index)
            else:
                design_row_indices.append(int(selected_row_index))
                if frequency_grid_row_counts and not (0 <= int(selected_row_index) < max(frequency_grid_row_counts)):
                    design_row_index_out_of_range.append(index)
            if selected_frequency is None:
                missing_selected_frequency.append(index)
            else:
                selected_hz = float(selected_frequency)
                selected_frequencies.append(selected_hz)
                if frequency is not None and abs(selected_hz - float(frequency)) > max(
                    1.0e-9,
                    abs(float(frequency)) * 1.0e-12,
                ):
                    selected_frequency_mismatch.append(index)
        if kind == "touchstone_row":
            if touchstone_row_index is None:
                missing_touchstone_row_index.append(index)
            else:
                touchstone_index = int(touchstone_row_index)
                touchstone_row_indices.append(touchstone_index)
                touchstone_row_index_records.append((index, touchstone_index))
            if touchstone_row_frequency is None:
                missing_touchstone_row_frequency.append(index)
            else:
                touchstone_frequency_hz = float(touchstone_row_frequency)
                touchstone_row_frequencies.append(touchstone_frequency_hz)
                touchstone_row_frequency_records.append((index, touchstone_frequency_hz))

        details.append({
            "index": index,
            "kind": kind,
            "project_id": None if project_id is None else str(project_id),
            "run_id": None if run_id is None else str(run_id),
            "export_id": None if export_id is None else str(export_id),
            "model_input_artifact_id": (
                None if model_input_artifact_id is None else str(model_input_artifact_id).strip()
            ),
            "model_input_digest": None if model_input_digest is None else str(model_input_digest).strip(),
            "model_input_path": None if model_input_path is None else str(model_input_path).strip(),
            "parameter_set_artifact_id": (
                None
                if parameter_set_artifact_id is None
                else str(parameter_set_artifact_id).strip()
            ),
            "parameter_set_digest": (
                None if parameter_set_digest is None else str(parameter_set_digest).strip()
            ),
            "parameter_set_path": None if parameter_set_path is None else str(parameter_set_path).strip(),
            "objective_observable_id": (
                None
                if objective_observable_id is None
                else str(objective_observable_id).strip()
            ),
            "objective_observable_family": (
                None
                if objective_observable_family is None
                else _norm(objective_observable_family)
            ),
            "source_tool": source_tool,
            "path": path,
            "gate_policy": gate_policy,
            "status": status,
            "pass": pass_flag,
            "design_frequency_hz": None if frequency is None else float(frequency),
            "row_count": None if row_count is None else int(row_count),
            "selected_row_index": None if selected_row_index is None else int(selected_row_index),
            "network_kind": None if not network_kind else network_kind,
            "port_order": port_order,
            "data_format": None if not data_format else data_format,
            "reference_impedance_ohm": None if reference_impedance is None else float(reference_impedance),
            "touchstone_option_line_artifact_id": (
                None
                if touchstone_option_line_artifact_id is None
                else str(touchstone_option_line_artifact_id).strip()
            ),
            "touchstone_option_line_digest": (
                None
                if touchstone_option_line_digest is None
                else str(touchstone_option_line_digest).strip()
            ),
            "reference_plane": None if reference_plane is None else str(reference_plane).strip(),
            "reference_plane_geometry_digest": (
                None
                if reference_plane_geometry_digest is None
                else str(reference_plane_geometry_digest).strip()
            ),
            "port_face_centers_xyz_m": (
                None
                if port_face_centers_value is None
                else [list(point) for point in port_face_centers_value]
            ),
            "port_mode_basis": None if port_mode_basis is None else str(port_mode_basis).strip(),
            "incident_wave_convention": None if incident_wave_convention is None else str(incident_wave_convention).strip(),
            "power_balance_basis": None if power_balance_basis is None else str(power_balance_basis).strip(),
            "touchstone_export_method": None if touchstone_export_method is None else str(touchstone_export_method).strip(),
            "export_recipe_artifact_id": (
                None if export_recipe_artifact_id is None else str(export_recipe_artifact_id).strip()
            ),
            "export_recipe_digest": None if export_recipe_digest is None else str(export_recipe_digest).strip(),
            "export_recipe_path": None if export_recipe_path is None else str(export_recipe_path).strip(),
            "result_tree_path": None if result_tree_path is None else str(result_tree_path).strip(),
            "result_item_id": None if result_item_id is None else str(result_item_id).strip(),
            "solver_setup_artifact_id": None if solver_setup_artifact_id is None else str(solver_setup_artifact_id).strip(),
            "mesh_setup_artifact_id": None if mesh_setup_artifact_id is None else str(mesh_setup_artifact_id).strip(),
            "port_definition_artifact_id": None if port_definition_artifact_id is None else str(port_definition_artifact_id).strip(),
            "excitation_setup_artifact_id": None if excitation_setup_artifact_id is None else str(excitation_setup_artifact_id).strip(),
            "frequency_grid_id": None if frequency_grid_id is None else str(frequency_grid_id).strip(),
            "frequency_grid_digest": None if frequency_grid_digest is None else str(frequency_grid_digest).strip(),
            "interpolation_policy": None if interpolation_policy is None else str(interpolation_policy).strip(),
            "selected_frequency_hz": None if selected_frequency is None else float(selected_frequency),
            "touchstone_row_index": None if touchstone_row_index is None else int(touchstone_row_index),
            "touchstone_row_frequency_hz": None if touchstone_row_frequency is None else float(touchstone_row_frequency),
            "touchstone_file_id": None if touchstone_file_id is None else str(touchstone_file_id).strip(),
            "touchstone_observable_id": None if touchstone_observable_id is None else str(touchstone_observable_id).strip(),
            "touchstone_observable_family": None if touchstone_observable_family is None else _norm(touchstone_observable_family),
            "touchstone_output_artifact_id": None if touchstone_output_artifact_id is None else str(touchstone_output_artifact_id).strip(),
            "touchstone_output_digest": None if touchstone_output_digest is None else str(touchstone_output_digest).strip(),
            "touchstone_output_path": None if touchstone_output_path is None else str(touchstone_output_path).strip(),
            "touchstone_output_schema_id": None if touchstone_output_schema_id is None else str(touchstone_output_schema_id).strip(),
            "touchstone_output_columns": touchstone_output_columns,
            "touchstone_output_units": touchstone_output_units,
            "created_at_utc": None if created_at_utc is None else str(created_at_utc).strip(),
            "run_timestamp_utc": None if run_timestamp_utc is None else str(run_timestamp_utc).strip(),
            "renormalized_reference_impedance_ohm": None if renormalized_reference_impedance is None else float(renormalized_reference_impedance),
            "renormalization_method": None if renormalization_method is None else str(renormalization_method).strip(),
            "renormalization_artifact_id": None if renormalization_artifact_id is None else str(renormalization_artifact_id).strip(),
            "deembedding_method": None if deembedding_method is None else str(deembedding_method).strip(),
            "deembedding_artifact_id": None if deembedding_artifact_id is None else str(deembedding_artifact_id).strip(),
            "deembedding_length_m": deembedding_length_value,
        })

    unique_project_ids = sorted(set(project_ids))
    unique_run_ids = sorted(set(run_ids))
    unique_export_ids = sorted(set(export_ids))
    unique_model_input_artifact_ids = sorted(set(model_input_artifact_ids))
    unique_model_input_digests = sorted(set(model_input_digests))
    unique_model_input_paths = sorted(set(model_input_paths))
    unique_parameter_set_artifact_ids = sorted(set(parameter_set_artifact_ids))
    unique_parameter_set_digests = sorted(set(parameter_set_digests))
    unique_parameter_set_paths = sorted(set(parameter_set_paths))
    unique_objective_observable_ids = sorted(set(objective_observable_ids))
    unique_objective_observable_families = sorted(set(objective_observable_families))
    unique_design_frequencies = sorted(set(design_frequencies))
    unique_network_kinds = sorted(set(network_kinds))
    unique_port_orders = sorted(set(port_orders))
    unique_data_formats = sorted(set(data_formats))
    unique_reference_impedances = sorted(set(reference_impedances))
    unique_touchstone_option_line_artifact_ids = sorted(set(touchstone_option_line_artifact_ids))
    unique_touchstone_option_line_digests = sorted(set(touchstone_option_line_digests))
    unique_reference_planes = sorted(set(reference_planes))
    unique_reference_plane_geometry_digests = sorted(set(reference_plane_geometry_digests))
    unique_port_face_centers_xyz_m = sorted(set(port_face_centers_xyz_m))
    unique_port_mode_bases = sorted(set(port_mode_bases))
    unique_touchstone_port_mode_basis_schema_ids = sorted(
        set(touchstone_port_mode_basis_schema_ids)
    )
    unique_incident_wave_conventions = sorted(set(incident_wave_conventions))
    unique_power_balance_bases = sorted(set(power_balance_bases))
    unique_touchstone_export_methods = sorted(set(touchstone_export_methods))
    unique_export_recipe_artifact_ids = sorted(set(export_recipe_artifact_ids))
    unique_export_recipe_digests = sorted(set(export_recipe_digests))
    unique_export_recipe_paths = sorted(set(export_recipe_paths))
    unique_result_tree_paths = sorted(set(result_tree_paths))
    unique_result_item_ids = sorted(set(result_item_ids))
    unique_solver_setup_artifact_ids = sorted(set(solver_setup_artifact_ids))
    unique_mesh_setup_artifact_ids = sorted(set(mesh_setup_artifact_ids))
    unique_port_definition_artifact_ids = sorted(set(port_definition_artifact_ids))
    unique_excitation_setup_artifact_ids = sorted(set(excitation_setup_artifact_ids))
    unique_frequency_grid_ids = sorted(set(frequency_grid_ids))
    unique_frequency_grid_digests = sorted(set(frequency_grid_digests))
    unique_interpolation_policies = sorted(set(interpolation_policies))
    unique_selected_frequencies = sorted(set(selected_frequencies))
    unique_touchstone_row_indices = sorted(set(touchstone_row_indices))
    unique_touchstone_row_frequencies = sorted(set(touchstone_row_frequencies))
    unique_touchstone_file_ids = sorted(set(touchstone_file_ids))
    unique_touchstone_observable_ids = sorted(set(touchstone_observable_ids))
    unique_touchstone_observable_families = sorted(set(touchstone_observable_families))
    unique_touchstone_output_artifact_ids = sorted(set(touchstone_output_artifact_ids))
    unique_touchstone_output_digests = sorted(set(touchstone_output_digests))
    unique_touchstone_output_paths = sorted(set(touchstone_output_paths))
    unique_touchstone_output_schema_ids = sorted(set(touchstone_output_schema_ids))
    unique_touchstone_convention_schema_ids = sorted(set(touchstone_convention_schema_ids))
    unique_touchstone_postprocess_row_convention_schema_ids = sorted(
        set(touchstone_postprocess_row_convention_schema_ids)
    )
    unique_touchstone_output_column_sets = sorted(set(touchstone_output_column_sets))
    unique_touchstone_output_unit_maps = sorted(set(touchstone_output_unit_maps))
    unique_created_at_utc_values = sorted(set(created_at_utc_values))
    unique_run_timestamp_utc_values = sorted(set(run_timestamp_utc_values))
    unique_renormalized_reference_impedances = sorted(set(renormalized_reference_impedances))
    unique_renormalization_methods = sorted(set(renormalization_methods))
    unique_renormalization_artifact_ids = sorted(set(renormalization_artifact_ids))
    unique_deembedding_methods = sorted(set(deembedding_methods))
    unique_deembedding_artifact_ids = sorted(set(deembedding_artifact_ids))
    unique_deembedding_lengths_m = sorted(set(deembedding_lengths_m))
    if touchstone_row_index_records and design_row_indices:
        expected_indices = set(design_row_indices)
        touchstone_row_index_mismatch = [
            index
            for index, row_index in touchstone_row_index_records
            if row_index not in expected_indices
        ]
    if touchstone_row_frequency_records and selected_frequencies:
        touchstone_row_frequency_mismatch = [
            index
            for index, row_frequency in touchstone_row_frequency_records
            if not any(
                abs(row_frequency - selected_frequency) <= max(1.0e-9, abs(selected_frequency) * 1.0e-12)
                for selected_frequency in selected_frequencies
            )
        ]
    required_set = set(required)
    present_set = set(kind_counts)
    checks = {
        "required_kinds_present": required_set.issubset(present_set),
        "no_unknown_kinds": not unknown_kinds,
        "project_ids_present": not missing_project_id,
        "project_ids_unique": len(unique_project_ids) == 1,
        "run_ids_present": not missing_run_id,
        "run_ids_unique": len(unique_run_ids) == 1,
        "export_ids_present": not missing_export_id,
        "export_ids_unique": len(unique_export_ids) == 1,
        "model_input_artifact_id_consistent_when_present": len(unique_model_input_artifact_ids) <= 1,
        "model_input_digest_consistent_when_present": len(unique_model_input_digests) <= 1,
        "model_input_path_consistent_when_present": len(unique_model_input_paths) <= 1,
        "parameter_set_artifact_id_consistent_when_present": len(unique_parameter_set_artifact_ids) <= 1,
        "parameter_set_digest_consistent_when_present": len(unique_parameter_set_digests) <= 1,
        "parameter_set_path_consistent_when_present": len(unique_parameter_set_paths) <= 1,
        "objective_observable_id_consistent_when_present": len(unique_objective_observable_ids) <= 1,
        "objective_observable_family_consistent_when_present": len(unique_objective_observable_families) <= 1,
        "source_tool_is_cst": not bad_source_tool,
        "paths_present": not missing_paths,
        "design_frequencies_present": not missing_frequency,
        "design_frequencies_unique": len(unique_design_frequencies) == 1,
        "upstream_gate_status_ok": not bad_upstream_status,
        "upstream_gate_policy_known": not bad_upstream_policy,
        "port_metadata_complete": not missing_port_metadata,
        "frequency_grid_brackets_design": not unbracketed_grid,
        "design_row_passive": not nonpassive_rows,
        "design_row_reciprocal": not nonreciprocal_rows,
        "touchstone_row_passive_when_present": not nonpassive_touchstone_rows,
        "touchstone_row_reciprocal_when_present": not nonreciprocal_touchstone_rows,
        "frequency_grid_row_count_recorded": not missing_frequency_grid_row_count,
        "design_row_index_recorded": not missing_design_row_index,
        "design_row_index_within_grid": not design_row_index_out_of_range,
        "touchstone_row_index_recorded_when_present": not missing_touchstone_row_index,
        "touchstone_row_frequency_recorded_when_present": not missing_touchstone_row_frequency,
        "touchstone_row_index_matches_design_row": not touchstone_row_index_mismatch,
        "touchstone_row_frequency_matches_selected_frequency": not touchstone_row_frequency_mismatch,
        "network_kind_recorded": not missing_network_kind,
        "network_kind_unique": len(unique_network_kinds) == 1,
        "port_order_recorded": not missing_port_order,
        "port_order_unique": len(unique_port_orders) == 1,
        "data_format_recorded": not missing_data_format,
        "data_format_unique": len(unique_data_formats) == 1,
        "reference_impedance_recorded": not missing_reference_impedance,
        "reference_impedance_unique": len(unique_reference_impedances) == 1,
        "touchstone_option_line_artifact_id_consistent_when_present": (
            len(unique_touchstone_option_line_artifact_ids) <= 1
        ),
        "touchstone_option_line_digest_consistent_when_present": (
            len(unique_touchstone_option_line_digests) <= 1
        ),
        "reference_plane_consistent_when_present": len(unique_reference_planes) <= 1,
        "reference_plane_geometry_digest_consistent_when_present": (
            len(unique_reference_plane_geometry_digests) <= 1
        ),
        "port_face_centers_xyz_consistent_when_present": len(unique_port_face_centers_xyz_m) <= 1,
        "port_mode_basis_consistent_when_present": len(unique_port_mode_bases) <= 1,
        "touchstone_port_mode_basis_schema_id_consistent_when_present": (
            len(unique_touchstone_port_mode_basis_schema_ids) <= 1
        ),
        "incident_wave_convention_consistent_when_present": len(unique_incident_wave_conventions) <= 1,
        "power_balance_basis_consistent_when_present": len(unique_power_balance_bases) <= 1,
        "touchstone_export_method_consistent_when_present": len(unique_touchstone_export_methods) <= 1,
        "export_recipe_artifact_id_consistent_when_present": len(unique_export_recipe_artifact_ids) <= 1,
        "export_recipe_digest_consistent_when_present": len(unique_export_recipe_digests) <= 1,
        "export_recipe_path_consistent_when_present": len(unique_export_recipe_paths) <= 1,
        "result_tree_path_consistent_when_present": len(unique_result_tree_paths) <= 1,
        "result_item_id_consistent_when_present": len(unique_result_item_ids) <= 1,
        "solver_setup_artifact_id_consistent_when_present": len(unique_solver_setup_artifact_ids) <= 1,
        "mesh_setup_artifact_id_consistent_when_present": len(unique_mesh_setup_artifact_ids) <= 1,
        "port_definition_artifact_id_consistent_when_present": len(unique_port_definition_artifact_ids) <= 1,
        "excitation_setup_artifact_id_consistent_when_present": len(unique_excitation_setup_artifact_ids) <= 1,
        "frequency_grid_id_consistent_when_present": len(unique_frequency_grid_ids) <= 1,
        "frequency_grid_digest_consistent_when_present": len(unique_frequency_grid_digests) <= 1,
        "interpolation_policy_consistent_when_present": len(unique_interpolation_policies) <= 1,
        "touchstone_file_id_consistent_when_present": len(unique_touchstone_file_ids) <= 1,
        "touchstone_observable_id_consistent_when_present": len(unique_touchstone_observable_ids) <= 1,
        "touchstone_observable_family_consistent_when_present": len(unique_touchstone_observable_families) <= 1,
        "touchstone_output_artifact_id_consistent_when_present": len(unique_touchstone_output_artifact_ids) <= 1,
        "touchstone_output_digest_consistent_when_present": len(unique_touchstone_output_digests) <= 1,
        "touchstone_output_path_consistent_when_present": len(unique_touchstone_output_paths) <= 1,
        "touchstone_output_schema_id_consistent_when_present": len(unique_touchstone_output_schema_ids) <= 1,
        "touchstone_convention_schema_id_consistent_when_present": len(unique_touchstone_convention_schema_ids) <= 1,
        "touchstone_postprocess_row_convention_schema_id_consistent_when_present": (
            len(unique_touchstone_postprocess_row_convention_schema_ids) <= 1
        ),
        "touchstone_output_columns_consistent_when_present": len(unique_touchstone_output_column_sets) <= 1,
        "touchstone_output_units_consistent_when_present": len(unique_touchstone_output_unit_maps) <= 1,
        "created_at_utc_consistent_when_present": len(unique_created_at_utc_values) <= 1,
        "run_timestamp_utc_consistent_when_present": len(unique_run_timestamp_utc_values) <= 1,
        "created_at_utc_parseable_when_present": not bad_created_at_utc,
        "run_timestamp_utc_parseable_when_present": not bad_run_timestamp_utc,
        "renormalized_reference_impedance_consistent_when_present": len(unique_renormalized_reference_impedances) <= 1,
        "renormalization_method_consistent_when_present": len(unique_renormalization_methods) <= 1,
        "renormalization_artifact_id_consistent_when_present": len(unique_renormalization_artifact_ids) <= 1,
        "renormalization_method_recorded_when_renormalized": not missing_renormalization_method,
        "renormalization_artifact_id_recorded_when_renormalized": not missing_renormalization_artifact_id,
        "deembedding_method_consistent_when_present": len(unique_deembedding_methods) <= 1,
        "deembedding_artifact_id_consistent_when_present": len(unique_deembedding_artifact_ids) <= 1,
        "deembedding_length_consistent_when_present": len(unique_deembedding_lengths_m) <= 1,
        "deembedding_method_recorded_when_deembedded": not missing_deembedding_method,
        "deembedding_artifact_id_recorded_when_deembedded": not missing_deembedding_artifact_id,
        "deembedding_length_recorded_when_deembedded": not missing_deembedding_length,
        "selected_frequency_recorded": not missing_selected_frequency,
        "selected_frequency_matches_design": not selected_frequency_mismatch,
    }
    if expected_project_id is not None:
        checks["expected_project_id_matches"] = unique_project_ids == [str(expected_project_id)]
    if expected_run_id is not None:
        checks["expected_run_id_matches"] = unique_run_ids == [str(expected_run_id)]
    if expected_export_id is not None:
        checks["expected_export_id_matches"] = unique_export_ids == [str(expected_export_id)]
    if require_model_input_artifact:
        checks["model_input_artifact_id_recorded_when_required"] = (
            not missing_model_input_artifact_id and bool(unique_model_input_artifact_ids)
        )
        checks["model_input_digest_recorded_when_required"] = (
            not missing_model_input_digest and bool(unique_model_input_digests)
        )
        checks["model_input_path_recorded_when_required"] = (
            not missing_model_input_path and bool(unique_model_input_paths)
        )
    if expected_model_input_artifact_id is not None:
        expected = str(expected_model_input_artifact_id).strip()
        checks["model_input_artifact_id_recorded_when_expected"] = (
            not missing_model_input_artifact_id and bool(unique_model_input_artifact_ids)
        )
        checks["expected_model_input_artifact_id_matches"] = (
            unique_model_input_artifact_ids == [expected]
            and not missing_model_input_artifact_id
        )
    if expected_model_input_digest is not None:
        expected = str(expected_model_input_digest).strip()
        checks["model_input_digest_recorded_when_expected"] = (
            not missing_model_input_digest and bool(unique_model_input_digests)
        )
        checks["expected_model_input_digest_matches"] = (
            unique_model_input_digests == [expected]
            and not missing_model_input_digest
        )
    if expected_model_input_path is not None:
        expected = str(expected_model_input_path).strip()
        checks["model_input_path_recorded_when_expected"] = (
            not missing_model_input_path and bool(unique_model_input_paths)
        )
        checks["expected_model_input_path_matches"] = (
            unique_model_input_paths == [expected]
            and not missing_model_input_path
        )
    if require_parameter_set_artifact:
        checks["parameter_set_artifact_id_recorded_when_required"] = (
            not missing_parameter_set_artifact_id
            and bool(unique_parameter_set_artifact_ids)
        )
        checks["parameter_set_digest_recorded_when_required"] = (
            not missing_parameter_set_digest
            and bool(unique_parameter_set_digests)
        )
        checks["parameter_set_path_recorded_when_required"] = (
            not missing_parameter_set_path
            and bool(unique_parameter_set_paths)
        )
    if expected_parameter_set_artifact_id is not None:
        expected = str(expected_parameter_set_artifact_id).strip()
        checks["parameter_set_artifact_id_recorded_when_expected"] = (
            not missing_parameter_set_artifact_id
            and bool(unique_parameter_set_artifact_ids)
        )
        checks["expected_parameter_set_artifact_id_matches"] = (
            unique_parameter_set_artifact_ids == [expected]
            and not missing_parameter_set_artifact_id
        )
    if expected_parameter_set_digest is not None:
        expected = str(expected_parameter_set_digest).strip()
        checks["parameter_set_digest_recorded_when_expected"] = (
            not missing_parameter_set_digest
            and bool(unique_parameter_set_digests)
        )
        checks["expected_parameter_set_digest_matches"] = (
            unique_parameter_set_digests == [expected]
            and not missing_parameter_set_digest
        )
    if expected_parameter_set_path is not None:
        expected = str(expected_parameter_set_path).strip()
        checks["parameter_set_path_recorded_when_expected"] = (
            not missing_parameter_set_path
            and bool(unique_parameter_set_paths)
        )
        checks["expected_parameter_set_path_matches"] = (
            unique_parameter_set_paths == [expected]
            and not missing_parameter_set_path
        )
    if expected_objective_observable_id is not None:
        expected = str(expected_objective_observable_id).strip()
        checks["objective_observable_id_recorded_when_expected"] = (
            not missing_objective_observable_id
            and bool(unique_objective_observable_ids)
        )
        checks["expected_objective_observable_id_matches"] = (
            unique_objective_observable_ids == [expected]
            and not missing_objective_observable_id
        )
    if expected_objective_observable_family is not None:
        expected = _norm(expected_objective_observable_family)
        checks["objective_observable_family_recorded_when_expected"] = (
            not missing_objective_observable_family
            and bool(unique_objective_observable_families)
        )
        checks["expected_objective_observable_family_matches"] = (
            unique_objective_observable_families == [expected]
            and not missing_objective_observable_family
        )
    if expected_design_frequency_hz is not None:
        expected = float(expected_design_frequency_hz)
        checks["expected_design_frequency_matches"] = (
            len(unique_design_frequencies) == 1
            and abs(unique_design_frequencies[0] - expected) <= max(1.0e-9, abs(expected) * 1.0e-12)
        )
        checks["selected_frequency_matches_expected_design"] = (
            len(unique_selected_frequencies) == 1
            and abs(unique_selected_frequencies[0] - expected) <= max(1.0e-9, abs(expected) * 1.0e-12)
        )
    if expected_network_kind is not None:
        expected = _network_kind(expected_network_kind)
        checks["expected_network_kind_matches"] = unique_network_kinds == [expected]
    if expected_port_order is not None:
        expected = tuple(_port_order(expected_port_order))
        checks["expected_port_order_matches"] = unique_port_orders == [expected]
    if expected_data_format is not None:
        expected = _data_format(expected_data_format)
        checks["expected_data_format_matches"] = unique_data_formats == [expected]
    if expected_reference_impedance_ohm is not None:
        expected = float(expected_reference_impedance_ohm)
        checks["expected_reference_impedance_matches"] = (
            len(unique_reference_impedances) == 1
            and abs(unique_reference_impedances[0] - expected) <= max(1.0e-12, abs(expected) * 1.0e-12)
        )
    if expected_touchstone_option_line_artifact_id is not None:
        expected = str(expected_touchstone_option_line_artifact_id).strip()
        checks["touchstone_option_line_artifact_id_recorded_when_expected"] = (
            not missing_touchstone_option_line_artifact_id
            and bool(unique_touchstone_option_line_artifact_ids)
        )
        checks["expected_touchstone_option_line_artifact_id_matches"] = (
            unique_touchstone_option_line_artifact_ids == [expected]
            and not missing_touchstone_option_line_artifact_id
        )
    if expected_touchstone_option_line_digest is not None:
        expected = str(expected_touchstone_option_line_digest).strip()
        checks["touchstone_option_line_digest_recorded_when_expected"] = (
            not missing_touchstone_option_line_digest
            and bool(unique_touchstone_option_line_digests)
        )
        checks["expected_touchstone_option_line_digest_matches"] = (
            unique_touchstone_option_line_digests == [expected]
            and not missing_touchstone_option_line_digest
        )
    if expected_reference_plane is not None:
        expected = str(expected_reference_plane).strip()
        checks["reference_plane_recorded_when_expected"] = not missing_reference_plane and bool(unique_reference_planes)
        checks["expected_reference_plane_matches"] = unique_reference_planes == [expected] and not missing_reference_plane
    if expected_reference_plane_geometry_digest is not None:
        expected = str(expected_reference_plane_geometry_digest).strip()
        checks["reference_plane_geometry_digest_recorded_when_expected"] = (
            not missing_reference_plane_geometry_digest
            and bool(unique_reference_plane_geometry_digests)
        )
        checks["expected_reference_plane_geometry_digest_matches"] = (
            unique_reference_plane_geometry_digests == [expected]
            and not missing_reference_plane_geometry_digest
        )
    if expected_port_face_centers_xyz_m is not None:
        expected = _coordinate_sequence(expected_port_face_centers_xyz_m)
        tolerance_m = 1.0e-9
        matches = (
            len(unique_port_face_centers_xyz_m) == 1
            and expected is not None
            and len(unique_port_face_centers_xyz_m[0]) == len(expected)
            and all(
                len(actual_point) == len(expected_point)
                and all(
                    abs(actual - expected_value) <= tolerance_m
                    for actual, expected_value in zip(actual_point, expected_point)
                )
                for actual_point, expected_point in zip(unique_port_face_centers_xyz_m[0], expected)
            )
            and not missing_port_face_centers_xyz_m
        )
        checks["port_face_centers_xyz_recorded_when_expected"] = (
            not missing_port_face_centers_xyz_m and bool(unique_port_face_centers_xyz_m)
        )
        checks["expected_port_face_centers_xyz_matches"] = matches
    if expected_port_mode_basis is not None:
        expected = str(expected_port_mode_basis).strip()
        checks["port_mode_basis_recorded_when_expected"] = not missing_port_mode_basis and bool(unique_port_mode_bases)
        checks["expected_port_mode_basis_matches"] = unique_port_mode_bases == [expected] and not missing_port_mode_basis
    if require_touchstone_port_mode_basis_schema:
        checks["touchstone_port_mode_basis_schema_id_recorded_when_required"] = (
            not missing_touchstone_port_mode_basis_schema_id
            and bool(unique_touchstone_port_mode_basis_schema_ids)
        )
    if expected_touchstone_port_mode_basis_schema_id is not None:
        expected = str(expected_touchstone_port_mode_basis_schema_id).strip()
        checks["touchstone_port_mode_basis_schema_id_recorded_when_expected"] = (
            not missing_touchstone_port_mode_basis_schema_id
            and bool(unique_touchstone_port_mode_basis_schema_ids)
        )
        checks["expected_touchstone_port_mode_basis_schema_id_matches"] = (
            unique_touchstone_port_mode_basis_schema_ids == [expected]
            and not missing_touchstone_port_mode_basis_schema_id
        )
    if expected_incident_wave_convention is not None:
        expected = str(expected_incident_wave_convention).strip()
        checks["incident_wave_convention_recorded_when_expected"] = (
            not missing_incident_wave_convention and bool(unique_incident_wave_conventions)
        )
        checks["expected_incident_wave_convention_matches"] = (
            unique_incident_wave_conventions == [expected] and not missing_incident_wave_convention
        )
    if expected_power_balance_basis is not None:
        expected = str(expected_power_balance_basis).strip()
        checks["power_balance_basis_recorded_when_expected"] = (
            not missing_power_balance_basis and bool(unique_power_balance_bases)
        )
        checks["expected_power_balance_basis_matches"] = (
            unique_power_balance_bases == [expected] and not missing_power_balance_basis
        )
    if expected_touchstone_export_method is not None:
        expected = str(expected_touchstone_export_method).strip()
        checks["touchstone_export_method_recorded_when_expected"] = (
            not missing_touchstone_export_method and bool(unique_touchstone_export_methods)
        )
        checks["expected_touchstone_export_method_matches"] = (
            unique_touchstone_export_methods == [expected] and not missing_touchstone_export_method
        )
    if require_export_recipe_artifact:
        checks["export_recipe_artifact_id_recorded_when_required"] = (
            not missing_export_recipe_artifact_id and bool(unique_export_recipe_artifact_ids)
        )
        checks["export_recipe_digest_recorded_when_required"] = (
            not missing_export_recipe_digest and bool(unique_export_recipe_digests)
        )
        checks["export_recipe_path_recorded_when_required"] = (
            not missing_export_recipe_path and bool(unique_export_recipe_paths)
        )
    if expected_export_recipe_artifact_id is not None:
        expected = str(expected_export_recipe_artifact_id).strip()
        checks["export_recipe_artifact_id_recorded_when_expected"] = (
            not missing_export_recipe_artifact_id and bool(unique_export_recipe_artifact_ids)
        )
        checks["expected_export_recipe_artifact_id_matches"] = (
            unique_export_recipe_artifact_ids == [expected]
            and not missing_export_recipe_artifact_id
        )
    if expected_export_recipe_digest is not None:
        expected = str(expected_export_recipe_digest).strip()
        checks["export_recipe_digest_recorded_when_expected"] = (
            not missing_export_recipe_digest and bool(unique_export_recipe_digests)
        )
        checks["expected_export_recipe_digest_matches"] = (
            unique_export_recipe_digests == [expected]
            and not missing_export_recipe_digest
        )
    if expected_export_recipe_path is not None:
        expected = str(expected_export_recipe_path).strip()
        checks["export_recipe_path_recorded_when_expected"] = (
            not missing_export_recipe_path and bool(unique_export_recipe_paths)
        )
        checks["expected_export_recipe_path_matches"] = (
            unique_export_recipe_paths == [expected]
            and not missing_export_recipe_path
        )
    if expected_result_tree_path is not None:
        expected = str(expected_result_tree_path).strip()
        checks["result_tree_path_recorded_when_expected"] = (
            not missing_result_tree_path and bool(unique_result_tree_paths)
        )
        checks["expected_result_tree_path_matches"] = (
            unique_result_tree_paths == [expected] and not missing_result_tree_path
        )
    if expected_result_item_id is not None:
        expected = str(expected_result_item_id).strip()
        checks["result_item_id_recorded_when_expected"] = (
            not missing_result_item_id and bool(unique_result_item_ids)
        )
        checks["expected_result_item_id_matches"] = (
            unique_result_item_ids == [expected] and not missing_result_item_id
        )
    if expected_solver_setup_artifact_id is not None:
        expected = str(expected_solver_setup_artifact_id).strip()
        checks["solver_setup_artifact_id_recorded_when_expected"] = (
            not missing_solver_setup_artifact_id and bool(unique_solver_setup_artifact_ids)
        )
        checks["expected_solver_setup_artifact_id_matches"] = (
            unique_solver_setup_artifact_ids == [expected] and not missing_solver_setup_artifact_id
        )
    if expected_mesh_setup_artifact_id is not None:
        expected = str(expected_mesh_setup_artifact_id).strip()
        checks["mesh_setup_artifact_id_recorded_when_expected"] = (
            not missing_mesh_setup_artifact_id and bool(unique_mesh_setup_artifact_ids)
        )
        checks["expected_mesh_setup_artifact_id_matches"] = (
            unique_mesh_setup_artifact_ids == [expected] and not missing_mesh_setup_artifact_id
        )
    if expected_port_definition_artifact_id is not None:
        expected = str(expected_port_definition_artifact_id).strip()
        checks["port_definition_artifact_id_recorded_when_expected"] = (
            not missing_port_definition_artifact_id and bool(unique_port_definition_artifact_ids)
        )
        checks["expected_port_definition_artifact_id_matches"] = (
            unique_port_definition_artifact_ids == [expected] and not missing_port_definition_artifact_id
        )
    if expected_excitation_setup_artifact_id is not None:
        expected = str(expected_excitation_setup_artifact_id).strip()
        checks["excitation_setup_artifact_id_recorded_when_expected"] = (
            not missing_excitation_setup_artifact_id and bool(unique_excitation_setup_artifact_ids)
        )
        checks["expected_excitation_setup_artifact_id_matches"] = (
            unique_excitation_setup_artifact_ids == [expected] and not missing_excitation_setup_artifact_id
        )
    if expected_frequency_grid_id is not None:
        expected = str(expected_frequency_grid_id).strip()
        checks["frequency_grid_id_recorded_when_expected"] = (
            not missing_frequency_grid_id and bool(unique_frequency_grid_ids)
        )
        checks["expected_frequency_grid_id_matches"] = (
            unique_frequency_grid_ids == [expected] and not missing_frequency_grid_id
        )
    if expected_frequency_grid_digest is not None:
        expected = str(expected_frequency_grid_digest).strip()
        checks["frequency_grid_digest_recorded_when_expected"] = (
            not missing_frequency_grid_digest and bool(unique_frequency_grid_digests)
        )
        checks["expected_frequency_grid_digest_matches"] = (
            unique_frequency_grid_digests == [expected] and not missing_frequency_grid_digest
        )
    if expected_interpolation_policy is not None:
        expected = str(expected_interpolation_policy).strip()
        checks["interpolation_policy_recorded_when_expected"] = (
            not missing_interpolation_policy and bool(unique_interpolation_policies)
        )
        checks["expected_interpolation_policy_matches"] = (
            unique_interpolation_policies == [expected] and not missing_interpolation_policy
        )
    if expected_touchstone_file_id is not None:
        expected = str(expected_touchstone_file_id).strip()
        checks["touchstone_file_id_recorded_when_expected"] = (
            not missing_touchstone_file_id and bool(unique_touchstone_file_ids)
        )
        checks["expected_touchstone_file_id_matches"] = (
            unique_touchstone_file_ids == [expected] and not missing_touchstone_file_id
        )
    if expected_touchstone_observable_id is not None:
        expected = str(expected_touchstone_observable_id).strip()
        checks["touchstone_observable_id_recorded_when_expected"] = (
            not missing_touchstone_observable_id and bool(unique_touchstone_observable_ids)
        )
        checks["expected_touchstone_observable_id_matches"] = (
            unique_touchstone_observable_ids == [expected] and not missing_touchstone_observable_id
        )
    if expected_touchstone_observable_family is not None:
        expected = _norm(expected_touchstone_observable_family)
        checks["touchstone_observable_family_recorded_when_expected"] = (
            not missing_touchstone_observable_family and bool(unique_touchstone_observable_families)
        )
        checks["expected_touchstone_observable_family_matches"] = (
            unique_touchstone_observable_families == [expected] and not missing_touchstone_observable_family
        )
    if require_touchstone_output_artifact:
        checks["touchstone_output_artifact_id_recorded_when_required"] = (
            not missing_touchstone_output_artifact_id and bool(unique_touchstone_output_artifact_ids)
        )
        checks["touchstone_output_digest_recorded_when_required"] = (
            not missing_touchstone_output_digest and bool(unique_touchstone_output_digests)
        )
        checks["touchstone_output_path_recorded_when_required"] = (
            not missing_touchstone_output_path and bool(unique_touchstone_output_paths)
        )
    if expected_touchstone_output_artifact_id is not None:
        expected = str(expected_touchstone_output_artifact_id).strip()
        checks["touchstone_output_artifact_id_recorded_when_expected"] = (
            not missing_touchstone_output_artifact_id and bool(unique_touchstone_output_artifact_ids)
        )
        checks["expected_touchstone_output_artifact_id_matches"] = (
            unique_touchstone_output_artifact_ids == [expected] and not missing_touchstone_output_artifact_id
        )
        checks["touchstone_output_path_recorded_when_expected"] = (
            not missing_touchstone_output_path and bool(unique_touchstone_output_paths)
        )
    if expected_touchstone_output_digest is not None:
        expected = str(expected_touchstone_output_digest).strip()
        checks["touchstone_output_digest_recorded_when_expected"] = (
            not missing_touchstone_output_digest and bool(unique_touchstone_output_digests)
        )
        checks["expected_touchstone_output_digest_matches"] = (
            unique_touchstone_output_digests == [expected] and not missing_touchstone_output_digest
        )
        checks["touchstone_output_path_recorded_when_expected"] = (
            not missing_touchstone_output_path and bool(unique_touchstone_output_paths)
        )
    if require_touchstone_output_schema:
        checks["touchstone_output_schema_id_recorded_when_required"] = (
            not missing_touchstone_output_schema_id
            and bool(unique_touchstone_output_schema_ids)
        )
        checks["touchstone_output_columns_recorded_when_required"] = (
            not missing_touchstone_output_columns
            and bool(unique_touchstone_output_column_sets)
        )
        checks["touchstone_output_units_recorded_when_required"] = (
            not missing_touchstone_output_units
            and bool(unique_touchstone_output_unit_maps)
        )
    if expected_touchstone_output_schema_id is not None:
        expected = str(expected_touchstone_output_schema_id).strip()
        checks["touchstone_output_schema_id_recorded_when_expected"] = (
            not missing_touchstone_output_schema_id
            and bool(unique_touchstone_output_schema_ids)
        )
        checks["expected_touchstone_output_schema_id_matches"] = (
            unique_touchstone_output_schema_ids == [expected]
            and not missing_touchstone_output_schema_id
        )
    if require_touchstone_convention_schema:
        checks["touchstone_convention_schema_id_recorded_when_required"] = (
            not missing_touchstone_convention_schema_id
            and bool(unique_touchstone_convention_schema_ids)
        )
    if expected_touchstone_convention_schema_id is not None:
        expected = str(expected_touchstone_convention_schema_id).strip()
        checks["touchstone_convention_schema_id_recorded_when_expected"] = (
            not missing_touchstone_convention_schema_id
            and bool(unique_touchstone_convention_schema_ids)
        )
        checks["expected_touchstone_convention_schema_id_matches"] = (
            unique_touchstone_convention_schema_ids == [expected]
            and not missing_touchstone_convention_schema_id
        )
    if require_touchstone_postprocess_row_convention_schema:
        checks["touchstone_postprocess_row_convention_schema_id_recorded_when_required"] = (
            not missing_touchstone_postprocess_row_convention_schema_id
            and bool(unique_touchstone_postprocess_row_convention_schema_ids)
        )
    if expected_touchstone_postprocess_row_convention_schema_id is not None:
        expected = str(expected_touchstone_postprocess_row_convention_schema_id).strip()
        checks["touchstone_postprocess_row_convention_schema_id_recorded_when_expected"] = (
            not missing_touchstone_postprocess_row_convention_schema_id
            and bool(unique_touchstone_postprocess_row_convention_schema_ids)
        )
        checks["expected_touchstone_postprocess_row_convention_schema_id_matches"] = (
            unique_touchstone_postprocess_row_convention_schema_ids == [expected]
            and not missing_touchstone_postprocess_row_convention_schema_id
        )
    expected_output_columns = tuple(_string_list(expected_touchstone_output_columns))
    if expected_output_columns:
        checks["touchstone_output_columns_recorded_when_expected"] = (
            not missing_touchstone_output_columns
            and bool(unique_touchstone_output_column_sets)
        )
        checks["expected_touchstone_output_columns_match"] = (
            unique_touchstone_output_column_sets == [expected_output_columns]
            and not missing_touchstone_output_columns
        )
    expected_output_units = tuple(sorted(_unit_mapping(expected_touchstone_output_units).items()))
    if expected_output_units:
        checks["touchstone_output_units_recorded_when_expected"] = (
            not missing_touchstone_output_units
            and bool(unique_touchstone_output_unit_maps)
        )
        checks["expected_touchstone_output_units_match"] = (
            unique_touchstone_output_unit_maps == [expected_output_units]
            and not missing_touchstone_output_units
        )
    if require_execution_metadata:
        checks["created_at_utc_recorded_when_required"] = (
            not missing_created_at_utc and bool(unique_created_at_utc_values)
        )
        checks["run_timestamp_utc_recorded_when_required"] = (
            not missing_run_timestamp_utc and bool(unique_run_timestamp_utc_values)
        )
    if expected_created_at_utc is not None:
        expected = str(expected_created_at_utc).strip()
        expected_dt = _parse_utc_like_datetime(expected)
        checks["created_at_utc_recorded_when_expected"] = (
            not missing_created_at_utc and bool(unique_created_at_utc_values)
        )
        checks["expected_created_at_utc_matches"] = (
            not missing_created_at_utc
            and not bad_created_at_utc
            and bool(unique_created_at_utc_values)
            and (
                unique_created_at_utc_values == [expected]
                if expected_dt is None
                else all(_parse_utc_like_datetime(value) == expected_dt for value in unique_created_at_utc_values)
            )
        )
    if expected_run_timestamp_utc is not None:
        expected = str(expected_run_timestamp_utc).strip()
        expected_dt = _parse_utc_like_datetime(expected)
        checks["run_timestamp_utc_recorded_when_expected"] = (
            not missing_run_timestamp_utc and bool(unique_run_timestamp_utc_values)
        )
        checks["expected_run_timestamp_utc_matches"] = (
            not missing_run_timestamp_utc
            and not bad_run_timestamp_utc
            and bool(unique_run_timestamp_utc_values)
            and (
                unique_run_timestamp_utc_values == [expected]
                if expected_dt is None
                else all(_parse_utc_like_datetime(value) == expected_dt for value in unique_run_timestamp_utc_values)
            )
        )
    if max_created_run_skew_s is not None:
        max_skew_s = float(max_created_run_skew_s)
        checks["created_run_timestamp_skew_recorded"] = (
            not missing_created_at_utc
            and not missing_run_timestamp_utc
            and not bad_created_at_utc
            and not bad_run_timestamp_utc
            and len(created_run_timestamp_skews_s) == len(rows_in)
        )
        checks["created_run_timestamp_skew_within_limit"] = (
            checks["created_run_timestamp_skew_recorded"]
            and all(skew <= max_skew_s for skew in created_run_timestamp_skews_s)
        )
    if expected_renormalized_reference_impedance_ohm is not None:
        expected = float(expected_renormalized_reference_impedance_ohm)
        checks["renormalized_reference_impedance_recorded_when_expected"] = (
            not missing_renormalized_reference_impedance and bool(unique_renormalized_reference_impedances)
        )
        checks["expected_renormalized_reference_impedance_matches"] = (
            len(unique_renormalized_reference_impedances) == 1
            and abs(unique_renormalized_reference_impedances[0] - expected)
            <= max(1.0e-12, abs(expected) * 1.0e-12)
            and not missing_renormalized_reference_impedance
        )
    if expected_renormalization_method is not None:
        expected = str(expected_renormalization_method).strip()
        checks["renormalization_method_recorded_when_expected"] = (
            not missing_renormalization_method and bool(unique_renormalization_methods)
        )
        checks["expected_renormalization_method_matches"] = (
            unique_renormalization_methods == [expected] and not missing_renormalization_method
        )
    if expected_renormalization_artifact_id is not None:
        expected = str(expected_renormalization_artifact_id).strip()
        checks["renormalization_artifact_id_recorded_when_expected"] = (
            not missing_renormalization_artifact_id and bool(unique_renormalization_artifact_ids)
        )
        checks["expected_renormalization_artifact_id_matches"] = (
            unique_renormalization_artifact_ids == [expected] and not missing_renormalization_artifact_id
        )
    if expected_deembedding_method is not None:
        expected = str(expected_deembedding_method).strip()
        checks["deembedding_method_recorded_when_expected"] = (
            not missing_deembedding_method and bool(unique_deembedding_methods)
        )
        checks["expected_deembedding_method_matches"] = (
            unique_deembedding_methods == [expected] and not missing_deembedding_method
        )
    if expected_deembedding_artifact_id is not None:
        expected = str(expected_deembedding_artifact_id).strip()
        checks["deembedding_artifact_id_recorded_when_expected"] = (
            not missing_deembedding_artifact_id and bool(unique_deembedding_artifact_ids)
        )
        checks["expected_deembedding_artifact_id_matches"] = (
            unique_deembedding_artifact_ids == [expected] and not missing_deembedding_artifact_id
        )
    if expected_deembedding_length_m is not None:
        expected = round(float(expected_deembedding_length_m), 15)
        checks["deembedding_length_recorded_when_expected"] = (
            not missing_deembedding_length and bool(unique_deembedding_lengths_m)
        )
        checks["expected_deembedding_length_matches"] = (
            len(unique_deembedding_lengths_m) == 1
            and abs(unique_deembedding_lengths_m[0] - expected)
            <= max(1.0e-15, abs(expected) * 1.0e-12)
            and not missing_deembedding_length
        )

    return {
        "policy": "cst_touchstone_solver_ready_manifest_gate",
        "required_kinds": list(required),
        "present_kinds": dict(sorted(kind_counts.items())),
        "project_ids": unique_project_ids,
        "run_ids": unique_run_ids,
        "export_ids": unique_export_ids,
        "model_input_artifact_ids": unique_model_input_artifact_ids,
        "model_input_digests": unique_model_input_digests,
        "model_input_paths": unique_model_input_paths,
        "parameter_set_artifact_ids": unique_parameter_set_artifact_ids,
        "parameter_set_digests": unique_parameter_set_digests,
        "parameter_set_paths": unique_parameter_set_paths,
        "objective_observable_ids": unique_objective_observable_ids,
        "objective_observable_families": unique_objective_observable_families,
        "design_frequencies_hz": unique_design_frequencies,
        "network_kinds": unique_network_kinds,
        "port_orders": [list(order) for order in unique_port_orders],
        "data_formats": unique_data_formats,
        "reference_impedances_ohm": unique_reference_impedances,
        "touchstone_option_line_artifact_ids": unique_touchstone_option_line_artifact_ids,
        "touchstone_option_line_digests": unique_touchstone_option_line_digests,
        "reference_planes": unique_reference_planes,
        "reference_plane_geometry_digests": unique_reference_plane_geometry_digests,
        "port_face_centers_xyz_m": [
            [list(point) for point in centers] for centers in unique_port_face_centers_xyz_m
        ],
        "port_mode_bases": unique_port_mode_bases,
        "touchstone_port_mode_basis_schema_ids": unique_touchstone_port_mode_basis_schema_ids,
        "incident_wave_conventions": unique_incident_wave_conventions,
        "power_balance_bases": unique_power_balance_bases,
        "touchstone_export_methods": unique_touchstone_export_methods,
        "export_recipe_artifact_ids": unique_export_recipe_artifact_ids,
        "export_recipe_digests": unique_export_recipe_digests,
        "export_recipe_paths": unique_export_recipe_paths,
        "result_tree_paths": unique_result_tree_paths,
        "result_item_ids": unique_result_item_ids,
        "solver_setup_artifact_ids": unique_solver_setup_artifact_ids,
        "mesh_setup_artifact_ids": unique_mesh_setup_artifact_ids,
        "port_definition_artifact_ids": unique_port_definition_artifact_ids,
        "excitation_setup_artifact_ids": unique_excitation_setup_artifact_ids,
        "frequency_grid_ids": unique_frequency_grid_ids,
        "frequency_grid_digests": unique_frequency_grid_digests,
        "interpolation_policies": unique_interpolation_policies,
        "selected_frequencies_hz": unique_selected_frequencies,
        "touchstone_row_indices": unique_touchstone_row_indices,
        "touchstone_row_frequencies_hz": unique_touchstone_row_frequencies,
        "touchstone_file_ids": unique_touchstone_file_ids,
        "touchstone_observable_ids": unique_touchstone_observable_ids,
        "touchstone_observable_families": unique_touchstone_observable_families,
        "touchstone_output_artifact_ids": unique_touchstone_output_artifact_ids,
        "touchstone_output_digests": unique_touchstone_output_digests,
        "touchstone_output_paths": unique_touchstone_output_paths,
        "touchstone_output_schema_ids": unique_touchstone_output_schema_ids,
        "touchstone_convention_schema_ids": unique_touchstone_convention_schema_ids,
        "touchstone_postprocess_row_convention_schema_ids": (
            unique_touchstone_postprocess_row_convention_schema_ids
        ),
        "touchstone_output_columns": [list(columns) for columns in unique_touchstone_output_column_sets],
        "touchstone_output_units": [dict(unit_map) for unit_map in unique_touchstone_output_unit_maps],
        "created_at_utc_values": unique_created_at_utc_values,
        "run_timestamp_utc_values": unique_run_timestamp_utc_values,
        "created_run_timestamp_skews_s": created_run_timestamp_skews_s,
        "max_created_run_skew_s": None if max_created_run_skew_s is None else float(max_created_run_skew_s),
        "renormalized_reference_impedances_ohm": unique_renormalized_reference_impedances,
        "renormalization_methods": unique_renormalization_methods,
        "renormalization_artifact_ids": unique_renormalization_artifact_ids,
        "deembedding_methods": unique_deembedding_methods,
        "deembedding_artifact_ids": unique_deembedding_artifact_ids,
        "deembedding_lengths_m": unique_deembedding_lengths_m,
        "expected_project_id": None if expected_project_id is None else str(expected_project_id),
        "expected_run_id": None if expected_run_id is None else str(expected_run_id),
        "expected_export_id": None if expected_export_id is None else str(expected_export_id),
        "expected_model_input_artifact_id": (
            None if expected_model_input_artifact_id is None else str(expected_model_input_artifact_id).strip()
        ),
        "expected_model_input_digest": (
            None if expected_model_input_digest is None else str(expected_model_input_digest).strip()
        ),
        "expected_model_input_path": (
            None if expected_model_input_path is None else str(expected_model_input_path).strip()
        ),
        "require_model_input_artifact": bool(require_model_input_artifact),
        "expected_parameter_set_artifact_id": (
            None
            if expected_parameter_set_artifact_id is None
            else str(expected_parameter_set_artifact_id).strip()
        ),
        "expected_parameter_set_digest": (
            None
            if expected_parameter_set_digest is None
            else str(expected_parameter_set_digest).strip()
        ),
        "expected_parameter_set_path": (
            None
            if expected_parameter_set_path is None
            else str(expected_parameter_set_path).strip()
        ),
        "expected_objective_observable_id": (
            None
            if expected_objective_observable_id is None
            else str(expected_objective_observable_id).strip()
        ),
        "expected_objective_observable_family": (
            None
            if expected_objective_observable_family is None
            else _norm(expected_objective_observable_family)
        ),
        "require_parameter_set_artifact": bool(require_parameter_set_artifact),
        "expected_design_frequency_hz": None if expected_design_frequency_hz is None else float(expected_design_frequency_hz),
        "expected_network_kind": None if expected_network_kind is None else _network_kind(expected_network_kind),
        "expected_port_order": None if expected_port_order is None else _port_order(expected_port_order),
        "expected_data_format": None if expected_data_format is None else _data_format(expected_data_format),
        "expected_reference_impedance_ohm": None if expected_reference_impedance_ohm is None else float(expected_reference_impedance_ohm),
        "expected_touchstone_option_line_artifact_id": (
            None
            if expected_touchstone_option_line_artifact_id is None
            else str(expected_touchstone_option_line_artifact_id).strip()
        ),
        "expected_touchstone_option_line_digest": (
            None
            if expected_touchstone_option_line_digest is None
            else str(expected_touchstone_option_line_digest).strip()
        ),
        "expected_reference_plane": None if expected_reference_plane is None else str(expected_reference_plane).strip(),
        "expected_reference_plane_geometry_digest": (
            None
            if expected_reference_plane_geometry_digest is None
            else str(expected_reference_plane_geometry_digest).strip()
        ),
        "expected_port_face_centers_xyz_m": (
            None
            if expected_port_face_centers_xyz_m is None
            else [list(point) for point in _coordinate_sequence(expected_port_face_centers_xyz_m)]
        ),
        "expected_port_mode_basis": None if expected_port_mode_basis is None else str(expected_port_mode_basis).strip(),
        "expected_touchstone_port_mode_basis_schema_id": (
            None
            if expected_touchstone_port_mode_basis_schema_id is None
            else str(expected_touchstone_port_mode_basis_schema_id).strip()
        ),
        "require_touchstone_port_mode_basis_schema": bool(
            require_touchstone_port_mode_basis_schema
        ),
        "expected_incident_wave_convention": None if expected_incident_wave_convention is None else str(expected_incident_wave_convention).strip(),
        "expected_power_balance_basis": None if expected_power_balance_basis is None else str(expected_power_balance_basis).strip(),
        "expected_touchstone_export_method": None if expected_touchstone_export_method is None else str(expected_touchstone_export_method).strip(),
        "expected_export_recipe_artifact_id": (
            None
            if expected_export_recipe_artifact_id is None
            else str(expected_export_recipe_artifact_id).strip()
        ),
        "expected_export_recipe_digest": (
            None
            if expected_export_recipe_digest is None
            else str(expected_export_recipe_digest).strip()
        ),
        "expected_export_recipe_path": (
            None
            if expected_export_recipe_path is None
            else str(expected_export_recipe_path).strip()
        ),
        "require_export_recipe_artifact": bool(require_export_recipe_artifact),
        "expected_result_tree_path": None if expected_result_tree_path is None else str(expected_result_tree_path).strip(),
        "expected_result_item_id": None if expected_result_item_id is None else str(expected_result_item_id).strip(),
        "expected_solver_setup_artifact_id": None if expected_solver_setup_artifact_id is None else str(expected_solver_setup_artifact_id).strip(),
        "expected_mesh_setup_artifact_id": None if expected_mesh_setup_artifact_id is None else str(expected_mesh_setup_artifact_id).strip(),
        "expected_port_definition_artifact_id": None if expected_port_definition_artifact_id is None else str(expected_port_definition_artifact_id).strip(),
        "expected_excitation_setup_artifact_id": None if expected_excitation_setup_artifact_id is None else str(expected_excitation_setup_artifact_id).strip(),
        "expected_frequency_grid_id": None if expected_frequency_grid_id is None else str(expected_frequency_grid_id).strip(),
        "expected_frequency_grid_digest": None if expected_frequency_grid_digest is None else str(expected_frequency_grid_digest).strip(),
        "expected_interpolation_policy": None if expected_interpolation_policy is None else str(expected_interpolation_policy).strip(),
        "expected_touchstone_file_id": None if expected_touchstone_file_id is None else str(expected_touchstone_file_id).strip(),
        "expected_touchstone_observable_id": None if expected_touchstone_observable_id is None else str(expected_touchstone_observable_id).strip(),
        "expected_touchstone_observable_family": None if expected_touchstone_observable_family is None else _norm(expected_touchstone_observable_family),
        "expected_touchstone_output_artifact_id": None if expected_touchstone_output_artifact_id is None else str(expected_touchstone_output_artifact_id).strip(),
        "expected_touchstone_output_digest": None if expected_touchstone_output_digest is None else str(expected_touchstone_output_digest).strip(),
        "require_touchstone_output_artifact": bool(require_touchstone_output_artifact),
        "expected_touchstone_output_schema_id": None if expected_touchstone_output_schema_id is None else str(expected_touchstone_output_schema_id).strip(),
        "expected_touchstone_output_columns": list(expected_output_columns),
        "expected_touchstone_output_units": dict(expected_output_units),
        "require_touchstone_output_schema": bool(require_touchstone_output_schema),
        "expected_touchstone_convention_schema_id": None if expected_touchstone_convention_schema_id is None else str(expected_touchstone_convention_schema_id).strip(),
        "require_touchstone_convention_schema": bool(require_touchstone_convention_schema),
        "expected_touchstone_postprocess_row_convention_schema_id": (
            None
            if expected_touchstone_postprocess_row_convention_schema_id is None
            else str(expected_touchstone_postprocess_row_convention_schema_id).strip()
        ),
        "require_touchstone_postprocess_row_convention_schema": bool(
            require_touchstone_postprocess_row_convention_schema
        ),
        "expected_created_at_utc": None if expected_created_at_utc is None else str(expected_created_at_utc).strip(),
        "expected_run_timestamp_utc": None if expected_run_timestamp_utc is None else str(expected_run_timestamp_utc).strip(),
        "require_execution_metadata": bool(require_execution_metadata),
        "expected_renormalized_reference_impedance_ohm": None if expected_renormalized_reference_impedance_ohm is None else float(expected_renormalized_reference_impedance_ohm),
        "expected_renormalization_method": None if expected_renormalization_method is None else str(expected_renormalization_method).strip(),
        "expected_renormalization_artifact_id": None if expected_renormalization_artifact_id is None else str(expected_renormalization_artifact_id).strip(),
        "expected_deembedding_method": None if expected_deembedding_method is None else str(expected_deembedding_method).strip(),
        "expected_deembedding_artifact_id": None if expected_deembedding_artifact_id is None else str(expected_deembedding_artifact_id).strip(),
        "expected_deembedding_length_m": None if expected_deembedding_length_m is None else float(expected_deembedding_length_m),
        "missing_project_id_rows": missing_project_id,
        "missing_run_id_rows": missing_run_id,
        "missing_export_id_rows": missing_export_id,
        "missing_model_input_artifact_id_rows": missing_model_input_artifact_id,
        "missing_model_input_digest_rows": missing_model_input_digest,
        "missing_model_input_path_rows": missing_model_input_path,
        "missing_parameter_set_artifact_id_rows": missing_parameter_set_artifact_id,
        "missing_parameter_set_digest_rows": missing_parameter_set_digest,
        "missing_parameter_set_path_rows": missing_parameter_set_path,
        "missing_objective_observable_id_rows": missing_objective_observable_id,
        "missing_objective_observable_family_rows": missing_objective_observable_family,
        "missing_frequency_rows": missing_frequency,
        "missing_network_kind_rows": missing_network_kind,
        "missing_port_order_rows": missing_port_order,
        "missing_data_format_rows": missing_data_format,
        "missing_reference_impedance_rows": missing_reference_impedance,
        "missing_touchstone_option_line_artifact_id_rows": missing_touchstone_option_line_artifact_id,
        "missing_touchstone_option_line_digest_rows": missing_touchstone_option_line_digest,
        "missing_reference_plane_rows": missing_reference_plane,
        "missing_reference_plane_geometry_digest_rows": missing_reference_plane_geometry_digest,
        "missing_port_face_centers_xyz_rows": missing_port_face_centers_xyz_m,
        "missing_port_mode_basis_rows": missing_port_mode_basis,
        "missing_touchstone_port_mode_basis_schema_id_rows": (
            missing_touchstone_port_mode_basis_schema_id
        ),
        "missing_incident_wave_convention_rows": missing_incident_wave_convention,
        "missing_power_balance_basis_rows": missing_power_balance_basis,
        "missing_touchstone_export_method_rows": missing_touchstone_export_method,
        "missing_export_recipe_artifact_id_rows": missing_export_recipe_artifact_id,
        "missing_export_recipe_digest_rows": missing_export_recipe_digest,
        "missing_export_recipe_path_rows": missing_export_recipe_path,
        "missing_result_tree_path_rows": missing_result_tree_path,
        "missing_result_item_id_rows": missing_result_item_id,
        "missing_solver_setup_artifact_id_rows": missing_solver_setup_artifact_id,
        "missing_mesh_setup_artifact_id_rows": missing_mesh_setup_artifact_id,
        "missing_port_definition_artifact_id_rows": missing_port_definition_artifact_id,
        "missing_excitation_setup_artifact_id_rows": missing_excitation_setup_artifact_id,
        "missing_frequency_grid_id_rows": missing_frequency_grid_id,
        "missing_frequency_grid_digest_rows": missing_frequency_grid_digest,
        "missing_interpolation_policy_rows": missing_interpolation_policy,
        "missing_touchstone_file_id_rows": missing_touchstone_file_id,
        "missing_touchstone_observable_id_rows": missing_touchstone_observable_id,
        "missing_touchstone_observable_family_rows": missing_touchstone_observable_family,
        "missing_touchstone_output_artifact_id_rows": missing_touchstone_output_artifact_id,
        "missing_touchstone_output_digest_rows": missing_touchstone_output_digest,
        "missing_touchstone_output_path_rows": missing_touchstone_output_path,
        "missing_touchstone_output_schema_id_rows": missing_touchstone_output_schema_id,
        "missing_touchstone_convention_schema_id_rows": missing_touchstone_convention_schema_id,
        "missing_touchstone_postprocess_row_convention_schema_id_rows": (
            missing_touchstone_postprocess_row_convention_schema_id
        ),
        "missing_touchstone_output_columns_rows": missing_touchstone_output_columns,
        "missing_touchstone_output_units_rows": missing_touchstone_output_units,
        "missing_created_at_utc_rows": missing_created_at_utc,
        "missing_run_timestamp_utc_rows": missing_run_timestamp_utc,
        "bad_created_at_utc_rows": bad_created_at_utc,
        "bad_run_timestamp_utc_rows": bad_run_timestamp_utc,
        "missing_renormalized_reference_impedance_rows": sorted(set(missing_renormalized_reference_impedance)),
        "missing_renormalization_method_rows": sorted(set(missing_renormalization_method)),
        "missing_renormalization_artifact_id_rows": sorted(set(missing_renormalization_artifact_id)),
        "missing_deembedding_method_rows": sorted(set(missing_deembedding_method)),
        "missing_deembedding_artifact_id_rows": sorted(set(missing_deembedding_artifact_id)),
        "missing_deembedding_length_rows": sorted(set(missing_deembedding_length)),
        "missing_selected_frequency_rows": missing_selected_frequency,
        "selected_frequency_mismatch_rows": selected_frequency_mismatch,
        "missing_touchstone_row_index_rows": missing_touchstone_row_index,
        "missing_touchstone_row_frequency_rows": missing_touchstone_row_frequency,
        "touchstone_row_index_mismatch_rows": touchstone_row_index_mismatch,
        "touchstone_row_frequency_mismatch_rows": touchstone_row_frequency_mismatch,
        "bad_source_tool_rows": bad_source_tool,
        "missing_path_rows": missing_paths,
        "unknown_kinds": unknown_kinds,
        "bad_upstream_status_rows": bad_upstream_status,
        "bad_upstream_policy_rows": bad_upstream_policy,
        "missing_port_metadata_rows": missing_port_metadata,
        "unbracketed_grid_rows": unbracketed_grid,
        "nonpassive_design_rows": nonpassive_rows,
        "nonreciprocal_design_rows": nonreciprocal_rows,
        "nonpassive_touchstone_rows": nonpassive_touchstone_rows,
        "nonreciprocal_touchstone_rows": nonreciprocal_touchstone_rows,
        "missing_frequency_grid_row_count_rows": missing_frequency_grid_row_count,
        "missing_design_row_index_rows": missing_design_row_index,
        "design_row_index_out_of_range_rows": design_row_index_out_of_range,
        "frequency_grid_row_counts": frequency_grid_row_counts,
        "design_row_indices": design_row_indices,
        "artifacts": details,
        "checks": checks,
        "status": "ok" if all(checks.values()) else "needs_attention",
        "version_note": (
            "Run after Touchstone port metadata, frequency-grid, and design-row "
            "preflights so S-parameter evidence cannot mix project/run/export "
            "identity, model-input artifact id/digest/path, parameter-set "
            "artifact id/digest/path, objective observable id/family, network kind, port order, data format, reference impedance, "
            "Touchstone option-line artifact/digest, "
            "reference plane, reference-plane geometry/port-face centers, "
            "port-mode basis, port-mode basis schema, incident-wave convention, "
            "power-balance basis, Touchstone export method, export recipe/macro/script "
            "artifact id/digest/path, solver setup artifact, "
            "CST result-tree path/item identity, mesh setup artifact, port definition "
            "artifact, excitation setup artifact, frequency-grid id, interpolation "
            "policy, raw Touchstone file identity, Touchstone observable "
            "id/family, optional Touchstone row index/frequency identity, "
            "execution created/run timestamp identity, "
            "optional renormalized Z0 identity, optional de-embedding identity, "
            "output artifact identity, output schema/column/unit identity, "
            "Touchstone physics convention schema, postprocess-row convention schema, "
            "or reuse an unbracketed, stale-frequency, "
            "active, or non-reciprocal row."
        ),
    }


def balanced_mcp_learning_profile_gate(profile):
    """Validate the ten-stage equal public/source MCP learning profile."""

    if not isinstance(profile, dict):
        raise ValueError("profile must be a mapping")
    expected_ids = [
        "baseline_gap",
        "source_controls",
        "structured_output",
        "input_validation",
        "security_boundary",
        "timeout_cancel_progress",
        "source_provenance",
        "artifact_feedback",
        "protocol_smoke",
        "balance_audit",
    ]
    stages = profile.get("stages")
    if not isinstance(stages, list):
        stages = []
    ids = [str(row.get("capability_id") or "") for row in stages if isinstance(row, dict)]
    rounds = [row.get("round") for row in stages if isinstance(row, dict)]
    controls_complete = all(
        str(row.get("positive_control") or "").strip()
        and str(row.get("negative_control") or "").strip()
        for row in stages
        if isinstance(row, dict)
    ) and len(stages) == 10
    roles = profile.get("workflow_roles")
    role_names = set(roles) if isinstance(roles, dict) else set()
    protocol = profile.get("protocol_policy")
    checks = {
        "schema_matches": profile.get("schema") == "cae-ai-lab.balanced-mcp-learning-profile.v1",
        "policy_matches": profile.get("policy") == "equal_capability_gain_v1",
        "stage_count_is_ten": profile.get("stage_count") == 10 and len(stages) == 10,
        "capability_ids_match": ids == expected_ids,
        "capability_ids_unique": len(set(ids)) == 10,
        "rounds_ordered": rounds == list(range(1, 11)),
        "controls_complete": controls_complete,
        "workflow_roles_complete": role_names == {"detect", "check", "run", "test"},
        "protocol_policy_complete": isinstance(protocol, dict)
        and all(str(protocol.get(key) or "").strip() for key in ("inspector_cli", "conformance", "fallback")),
        "owners_recorded": bool(str(profile.get("public_owner") or "").strip())
        and bool(str(profile.get("source_owner") or "").strip()),
        "completion_rule_recorded": bool(str(profile.get("completion_rule") or "").strip()),
    }
    return {
        "policy": "balanced_mcp_learning_profile_gate",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "server": profile.get("server"),
        "capability_ids": ids,
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
    }


def licensed_solver_agentic_profile_gate(
    agentic_profile_schema,
    recommended_first_calls,
    *,
    required_first_call_kinds=("status", "license", "live"),
    preflight_consumes_license=None,
    default_gui=None,
    owned_instance_cleanup_only=None,
):
    """Validate a seat-aware agentic profile for a licensed solver.

    The gate is vendor-neutral. It requires status/license/live discovery before
    execution, a non-consuming preflight, headless-by-default automation, and
    cleanup restricted to an instance owned by the helper.
    """

    schema = str(agentic_profile_schema or "").strip()
    calls = [str(item).strip() for item in (recommended_first_calls or []) if str(item).strip()]
    required = [str(item).strip().lower() for item in required_first_call_kinds if str(item).strip()]
    checks = {
        "agentic_profile_schema_recorded": bool(schema),
        "recommended_first_calls_recorded": bool(calls),
        "required_first_call_kinds_recorded": bool(required),
        "recommended_first_calls_cover_required_kinds": bool(required) and all(
            any(kind in call.lower() for call in calls) for kind in required
        ),
        "preflight_license_policy_recorded": preflight_consumes_license is not None,
        "preflight_is_seat_non_consuming": preflight_consumes_license is False,
        "default_gui_policy_recorded": default_gui is not None,
        "headless_by_default": default_gui is False,
        "cleanup_ownership_policy_recorded": owned_instance_cleanup_only is not None,
        "cleanup_is_owned_instance_only": owned_instance_cleanup_only is True,
    }
    return {
        "policy": "licensed_solver_agentic_profile_gate",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "agentic_profile_schema": schema,
        "recommended_first_calls": calls,
        "required_first_call_kinds": required,
        "preflight_consumes_license": preflight_consumes_license,
        "default_gui": default_gui,
        "owned_instance_cleanup_only": owned_instance_cleanup_only,
        "checks": checks,
        "notes": [
            "Discover status and license availability before a licensed live operation.",
            "Never quit a user-owned solver instance from an MCP smoke helper.",
        ],
    }


def mcp_artifact_access_behavior_gate(
    *,
    path_registered,
    path_under_allowed_root,
    access_status,
    sensitive_payload_visible=False,
):
    """Check that an MCP path reader enforces its artifact boundary.

    A path request is allowed only when the artifact is registered or resides
    below an explicit owned root. Unregistered out-of-root requests must return
    a denial without echoing data from the requested file.
    """

    registered = bool(path_registered)
    under_root = bool(path_under_allowed_root)
    observed = str(access_status or "").strip().lower()
    expected = "ok" if registered or under_root else "denied"
    checks = {
        "access_status_recorded": observed in {"ok", "denied"},
        "decision_matches_registration_or_root": observed == expected,
        "unregistered_out_of_root_path_denied": registered or under_root or observed == "denied",
        "denied_response_hides_sensitive_payload": observed != "denied" or not bool(sensitive_payload_visible),
    }
    return {
        "policy": "mcp_artifact_access_behavior_gate",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "path_registered": registered,
        "path_under_allowed_root": under_root,
        "access_status": observed,
        "expected_access_status": expected,
        "sensitive_payload_visible": bool(sensitive_payload_visible),
        "checks": checks,
        "notes": [
            "Prefer text/data inputs for public MCP tools; path readers need owned roots or manifest registration.",
            "A denial response may report policy and allowed roots but must not echo requested-file payload values.",
        ],
    }


def shared_solver_session_health_gate(
    connected,
    api_visible,
    discovered_engines=None,
    matlab_engine_find_matlab=None,
    shared_engine_name="",
    livelink_matlab_pid=None,
    status="",
    started_new_process=False,
    killed_process=False,
    direct_discovery_status="",
    passive_diagnostic_verdict="",
    version="",
    model_tags=None,
    passive_server_pid=None,
    passive_matlab_pid=None,
    passive_worker_pid=None,
    passive_matlab_parent_pid=None,
    passive_worker_parent_pid=None,
    target_port=None,
    established_connection_count=None,
    shared_engine_eval="",
    livelink_out_fields=None,
    matlab_version_source="",
    solver_version_source="",
    passive_diagnostic_timestamp="",
    passive_machine_policy="",
    passive_port_owner_pid=None,
    shared_engine_eval_status="",
    previous_shared_engine_eval_status="",
    shared_engine_eval_timeout_s=None,
    shared_engine_eval_timeout_mode="",
    matlab_process_count=None,
    matlab_mcp_server_count=None,
    livelink_candidate_count=None,
    session_selection_basis="",
    agentic_profile_schema="",
    recommended_first_calls=None,
    required_first_call_kinds=None,
    session_policy_mode="",
    profile_allows_new_solver_process=None,
    profile_allows_kill_solver_process=None,
):
    """Check solver-session health without treating it as physics validation.

    This public-safe gate records whether an external solver session was reused
    cleanly before a numerical validation row is trusted.  It is deliberately
    generic: COMSOL LiveLink, MATLAB Engine, Jupyter kernels, and similar
    long-lived solver sessions can all use the same separation between session
    health and physics residuals.
    """

    engines = [] if discovered_engines is None else [str(item) for item in discovered_engines]
    find_matlab_engines = (
        []
        if matlab_engine_find_matlab is None
        else [str(item) for item in matlab_engine_find_matlab]
    )
    find_matlab_requested = matlab_engine_find_matlab is not None
    shared_engine_name_text = str(shared_engine_name or "").strip()
    state = str(status or "").strip().lower()
    direct_status = str(direct_discovery_status or "").strip()
    direct_status_norm = direct_status.lower()
    passive_verdict = str(passive_diagnostic_verdict or "").strip()
    version_text = str(version or "").strip()
    tags = None if model_tags is None else [str(item) for item in model_tags]
    out_fields = [] if livelink_out_fields is None else [str(item) for item in livelink_out_fields]
    out_field_keys = {item.strip().lower() for item in out_fields}
    matlab_version_source_text = str(matlab_version_source or "").strip()
    solver_version_source_text = str(solver_version_source or "").strip()
    diagnostic_timestamp_text = str(passive_diagnostic_timestamp or "").strip()
    machine_policy_text = str(passive_machine_policy or "").strip()
    shared_eval_status_text = str(shared_engine_eval_status or "").strip().lower().replace("-", "_")
    previous_shared_eval_status_text = (
        str(previous_shared_engine_eval_status or "").strip().lower().replace("-", "_")
    )
    shared_eval_timeout_mode_text = str(shared_engine_eval_timeout_mode or "").strip()
    selection_basis_text = str(session_selection_basis or "").strip()
    profile_schema_text = str(agentic_profile_schema or "").strip()
    first_calls = (
        []
        if recommended_first_calls is None
        else [str(item).strip() for item in recommended_first_calls if str(item).strip()]
    )
    required_call_kinds = (
        []
        if required_first_call_kinds is None
        else [str(item).strip().lower() for item in required_first_call_kinds if str(item).strip()]
    )
    session_policy_text = str(session_policy_mode or "").strip()

    def positive_int_or_none(value):
        if value is None:
            return None
        try:
            converted = int(value)
        except (TypeError, ValueError):
            return None
        return converted if converted > 0 else None

    def positive_float_or_none(value):
        if value is None:
            return None
        try:
            converted = float(value)
        except (TypeError, ValueError):
            return None
        return converted if converted > 0.0 else None

    passive_server_pid_int = positive_int_or_none(passive_server_pid)
    passive_matlab_pid_int = positive_int_or_none(passive_matlab_pid)
    passive_worker_pid_int = positive_int_or_none(passive_worker_pid)
    passive_matlab_parent_pid_int = positive_int_or_none(passive_matlab_parent_pid)
    passive_worker_parent_pid_int = positive_int_or_none(passive_worker_parent_pid)
    livelink_matlab_pid_int = positive_int_or_none(livelink_matlab_pid)
    target_port_int = positive_int_or_none(target_port)
    established_count_int = positive_int_or_none(established_connection_count)
    passive_port_owner_pid_int = positive_int_or_none(passive_port_owner_pid)
    matlab_process_count_int = positive_int_or_none(matlab_process_count)
    matlab_mcp_server_count_int = positive_int_or_none(matlab_mcp_server_count)
    livelink_candidate_count_int = positive_int_or_none(livelink_candidate_count)
    shared_eval_timeout = positive_float_or_none(shared_engine_eval_timeout_s)
    shared_eval_path = str(shared_engine_eval or "").strip()

    checks = {
        "session_connected": bool(connected),
        "api_visible": bool(api_visible),
        "engine_discovered": len(engines) > 0,
        "status_allows_reuse": state in {"already-connected", "connected", "live", "ok"},
        "started_no_new_process": not bool(started_new_process),
        "killed_no_process": not bool(killed_process),
    }
    if find_matlab_requested:
        checks["matlab_engine_find_matlab_recorded"] = len(find_matlab_engines) > 0
        checks["selected_shared_engine_visible_in_find_matlab"] = (
            bool(shared_engine_name_text)
            and shared_engine_name_text in find_matlab_engines
        )
    if direct_status:
        checks["direct_discovery_status_recorded"] = True
    if shared_engine_name_text or livelink_matlab_pid is not None:
        checks["shared_engine_name_recorded"] = bool(shared_engine_name_text)
        checks["shared_engine_name_discovered"] = shared_engine_name_text in engines
        checks["livelink_matlab_pid_recorded"] = livelink_matlab_pid_int is not None
        checks["livelink_matlab_pid_matches_worker_pid"] = (
            livelink_matlab_pid_int is not None
            and (
                passive_worker_pid_int is None
                or livelink_matlab_pid_int == passive_worker_pid_int
            )
        )
        if shared_engine_name_text.lower().startswith("matlab_"):
            try:
                engine_pid = int(shared_engine_name_text.split("_", 1)[1])
            except (IndexError, ValueError):
                engine_pid = None
            checks["shared_engine_name_matches_pid"] = (
                engine_pid is not None
                and livelink_matlab_pid_int is not None
                and engine_pid == livelink_matlab_pid_int
            )
    if passive_verdict:
        checks["passive_diagnostic_recorded"] = True
    if version_text:
        checks["version_recorded"] = True
    if out_fields:
        checks["livelink_core_fields_present"] = {"connected", "status", "reason"}.issubset(out_field_keys)
        checks["livelink_version_field_not_required"] = (
            "version" not in out_field_keys
            and bool(matlab_version_source_text)
            and bool(solver_version_source_text)
        )
    if matlab_version_source_text or solver_version_source_text:
        checks["matlab_version_source_recorded"] = bool(matlab_version_source_text)
        checks["solver_version_source_recorded"] = bool(solver_version_source_text)
    if diagnostic_timestamp_text:
        checks["passive_diagnostic_timestamp_recorded"] = True
    if machine_policy_text:
        policy_norm = machine_policy_text.lower()
        checks["passive_no_tcp_probe_policy_recorded"] = (
            "passive" in policy_norm and "no" in policy_norm and "tcp" in policy_norm
        )
    if tags is not None:
        checks["model_tags_recorded"] = True
        checks["model_tags_are_introspection_only"] = True
    passive_evidence_requested = any(
        value is not None
        for value in (
            passive_server_pid,
            passive_matlab_pid,
            passive_worker_pid,
            passive_matlab_parent_pid,
            passive_worker_parent_pid,
            target_port,
            established_connection_count,
        )
    ) or bool(shared_eval_path)
    if passive_evidence_requested:
        checks["passive_server_pid_recorded"] = passive_server_pid_int is not None
        checks["passive_matlab_pid_recorded"] = passive_matlab_pid_int is not None
        checks["passive_worker_pid_recorded"] = passive_worker_pid_int is not None
        if passive_matlab_parent_pid is not None:
            checks["passive_matlab_parent_pid_recorded"] = passive_matlab_parent_pid_int is not None
            checks["livelink_matlab_parent_is_server"] = (
                passive_matlab_parent_pid_int is not None
                and passive_server_pid_int is not None
                and passive_matlab_parent_pid_int == passive_server_pid_int
            )
        if passive_worker_parent_pid is not None:
            checks["passive_worker_parent_pid_recorded"] = passive_worker_parent_pid_int is not None
            checks["worker_parent_is_livelink_matlab"] = (
                passive_worker_parent_pid_int is not None
                and passive_matlab_pid_int is not None
                and passive_worker_parent_pid_int == passive_matlab_pid_int
            )
        checks["target_port_recorded"] = target_port_int is not None
        checks["established_connection_recorded"] = established_count_int is not None
        checks["shared_engine_eval_recorded"] = bool(shared_eval_path)
        checks["passive_session_evidence_complete"] = all(
            checks[key]
            for key in (
                "passive_server_pid_recorded",
                "passive_matlab_pid_recorded",
                "passive_worker_pid_recorded",
                "target_port_recorded",
                "established_connection_recorded",
                "shared_engine_eval_recorded",
            )
        )
    if passive_port_owner_pid is not None:
        checks["passive_port_owner_pid_recorded"] = passive_port_owner_pid_int is not None
        checks["target_port_owned_by_server_pid"] = (
            passive_port_owner_pid_int is not None
            and passive_server_pid_int is not None
            and passive_port_owner_pid_int == passive_server_pid_int
        )
    if shared_eval_status_text:
        checks["shared_engine_eval_status_recorded"] = shared_eval_status_text in {
            "ok",
            "timeout",
            "failed",
            "not_run",
        }
    if previous_shared_eval_status_text:
        checks["previous_shared_engine_eval_status_recorded"] = previous_shared_eval_status_text in {
            "ok",
            "timeout",
            "failed",
            "not_run",
            "error",
            "needs_attention",
        }
        checks["shared_engine_eval_recovered_after_timeout"] = (
            previous_shared_eval_status_text != "timeout"
            or (
                shared_eval_status_text == "ok"
                and bool(connected)
                and bool(api_visible)
            )
        )
    if shared_eval_status_text == "timeout":
        checks["shared_engine_eval_timeout_recorded"] = shared_eval_timeout is not None
        checks["shared_engine_eval_timeout_mode_recorded"] = bool(shared_eval_timeout_mode_text)
        checks["shared_engine_eval_timeout_is_diagnostic"] = (
            not bool(started_new_process) and not bool(killed_process)
        )
    if direct_status_norm == "no matlab session discovered":
        find_matlab_evidence_ok = (
            not find_matlab_requested
            or (
                bool(shared_engine_name_text)
                and shared_engine_name_text in find_matlab_engines
            )
        )
        checks["direct_discovery_false_negative_reconciled"] = (
            len(engines) > 0
            and bool(passive_verdict)
            and bool(connected)
            and find_matlab_evidence_ok
        )
        checks["direct_discovery_false_negative_has_selected_engine"] = bool(shared_engine_name_text)
        checks["direct_discovery_false_negative_has_ok_shared_eval"] = (
            bool(shared_eval_path) and shared_eval_status_text == "ok"
        )
        if find_matlab_requested:
            checks["direct_discovery_false_negative_has_find_matlab_engine"] = (
                shared_engine_name_text in find_matlab_engines
            )
    multi_process_context_requested = any(
        value is not None
        for value in (
            matlab_process_count,
            matlab_mcp_server_count,
            livelink_candidate_count,
        )
    ) or bool(selection_basis_text)
    if multi_process_context_requested:
        basis_norm = selection_basis_text.lower()
        checks["matlab_process_count_recorded"] = matlab_process_count_int is not None
        checks["matlab_mcp_server_count_recorded"] = matlab_mcp_server_count_int is not None
        checks["livelink_candidate_count_recorded"] = livelink_candidate_count_int is not None
        checks["session_selection_basis_recorded"] = bool(selection_basis_text)
        checks["selection_uses_parent_or_port_chain"] = (
            "parent" in basis_norm or "port" in basis_norm or "connection" in basis_norm
        )
        if matlab_process_count_int is not None and livelink_candidate_count_int is not None:
            checks["multiple_matlab_processes_do_not_create_ambiguity"] = (
                livelink_candidate_count_int >= 1
                and matlab_process_count_int >= livelink_candidate_count_int
                and checks["selection_uses_parent_or_port_chain"]
            )
    profile_contract_requested = any(
        (
            profile_schema_text,
            first_calls,
            required_call_kinds,
            session_policy_text,
            profile_allows_new_solver_process is not None,
            profile_allows_kill_solver_process is not None,
        )
    )
    if profile_contract_requested:
        policy_norm = session_policy_text.lower()
        checks["agentic_profile_schema_recorded"] = bool(profile_schema_text)
        checks["recommended_first_calls_recorded"] = bool(first_calls)
        checks["required_first_call_kinds_recorded"] = bool(required_call_kinds)
        checks["recommended_first_calls_cover_required_kinds"] = bool(required_call_kinds) and all(
            any(kind in call.lower() for call in first_calls)
            for kind in required_call_kinds
        )
        checks["session_policy_mode_recorded"] = bool(session_policy_text)
        checks["session_policy_prefers_reuse"] = any(
            token in policy_norm for token in ("attach", "reuse", "existing", "shared")
        )
        checks["profile_new_process_policy_recorded"] = (
            profile_allows_new_solver_process is not None
        )
        checks["profile_prohibits_new_solver_process"] = (
            profile_allows_new_solver_process is False
        )
        checks["profile_kill_policy_recorded"] = (
            profile_allows_kill_solver_process is not None
        )
        checks["profile_prohibits_killing_solver_process"] = (
            profile_allows_kill_solver_process is False
        )
    return {
        "policy": "shared_solver_session_health_gate",
        "connected": bool(connected),
        "api_visible": bool(api_visible),
        "discovered_engines": engines,
        "matlab_engine_find_matlab": find_matlab_engines,
        "shared_engine_name": shared_engine_name_text,
        "livelink_matlab_pid": livelink_matlab_pid_int,
        "status": state,
        "started_new_process": bool(started_new_process),
        "killed_process": bool(killed_process),
        "direct_discovery_status": direct_status,
        "passive_diagnostic_verdict": passive_verdict,
        "version": version_text,
        "model_tags": tags,
        "passive_server_pid": passive_server_pid_int,
        "passive_matlab_pid": passive_matlab_pid_int,
        "passive_worker_pid": passive_worker_pid_int,
        "passive_matlab_parent_pid": passive_matlab_parent_pid_int,
        "passive_worker_parent_pid": passive_worker_parent_pid_int,
        "target_port": target_port_int,
        "established_connection_count": established_count_int,
        "shared_engine_eval": shared_eval_path,
        "livelink_out_fields": out_fields,
        "matlab_version_source": matlab_version_source_text,
        "solver_version_source": solver_version_source_text,
        "passive_diagnostic_timestamp": diagnostic_timestamp_text,
        "passive_machine_policy": machine_policy_text,
        "passive_port_owner_pid": passive_port_owner_pid_int,
        "shared_engine_eval_status": shared_eval_status_text,
        "previous_shared_engine_eval_status": previous_shared_eval_status_text,
        "shared_engine_eval_timeout_s": shared_eval_timeout,
        "shared_engine_eval_timeout_mode": shared_eval_timeout_mode_text,
        "matlab_process_count": matlab_process_count_int,
        "matlab_mcp_server_count": matlab_mcp_server_count_int,
        "livelink_candidate_count": livelink_candidate_count_int,
        "session_selection_basis": selection_basis_text,
        "agentic_profile_schema": profile_schema_text,
        "recommended_first_calls": first_calls,
        "required_first_call_kinds": required_call_kinds,
        "session_policy_mode": session_policy_text,
        "profile_allows_new_solver_process": profile_allows_new_solver_process,
        "profile_allows_kill_solver_process": profile_allows_kill_solver_process,
        "checks": checks,
        "status_label": "ok" if all(checks.values()) else "needs_attention",
        "notes": [
            "session health is a preflight, not a physical residual",
            "already-connected is valid reuse evidence when no new process was started or killed",
            "model tags are useful session identity evidence but are not the connection pass/fail gate",
            "a direct discovery false negative can be healthy when passive diagnostics and shared-engine evaluation agree",
            "close a direct discovery false negative only after selecting the shared engine and running an ok shared-engine eval",
            "record passive process and port evidence when a direct discovery path misses a reusable session",
            "do not assume cc_livelink returns a version field; record MATLAB and solver version sources separately",
            "a shared-engine timeout is diagnostic evidence, not a reason to kill solver-owned processes",
            "a previous shared-engine timeout is closed only by a later ok eval with visible API evidence",
            "when many MATLAB processes exist, select the LiveLink target by parent/port evidence, not by the first MATLAB executable",
            "record the selected shared engine name and verify it matches the LiveLink MATLAB worker PID when available",
            "record matlab.engine.find_matlab() output separately when direct MCP discovery misses a shared MATLAB session",
            "record the COMSOL server -> LiveLink MATLAB -> shared worker MATLAB parent chain when passive diagnostics expose it",
            "an external solver MCP profile should expose health-first calls and prohibit unrequested process start or kill",
        ],
    }


def solver_submodel_boundary_handoff_gate(
    parent_model_id,
    submodel_region_id,
    zoom_boundary_id,
    boundary_transfer_quantity,
    boundary_transfer_error_estimate,
    local_refinement_rule,
    target_observable_id,
    *,
    parent_mesh_id="",
    local_mesh_id="",
    boundary_trace_id="",
    boundary_condition_source="",
    boundary_transfer_error_unit="",
    expected_boundary_transfer_quantity=None,
    max_boundary_transfer_error=None,
):
    """Check that a local submodel solve records its boundary handoff.

    A refined local/zoom solve is not self-validating just because the local
    mesh is fine.  Its artifact must also identify the global/parent model and
    the boundary data inherited by the local solve, including an error estimate
    or budget for that transfer.
    """

    parent_model_text = str(parent_model_id or "").strip()
    submodel_region_text = str(submodel_region_id or "").strip()
    zoom_boundary_text = str(zoom_boundary_id or "").strip()
    transfer_quantity_text = str(boundary_transfer_quantity or "").strip()
    refinement_rule_text = str(local_refinement_rule or "").strip()
    observable_text = str(target_observable_id or "").strip()
    parent_mesh_text = str(parent_mesh_id or "").strip()
    local_mesh_text = str(local_mesh_id or "").strip()
    trace_text = str(boundary_trace_id or "").strip()
    source_text = str(boundary_condition_source or "").strip()
    error_unit_text = str(boundary_transfer_error_unit or "").strip()
    expected_quantity_text = (
        None
        if expected_boundary_transfer_quantity is None
        else str(expected_boundary_transfer_quantity or "").strip()
    )

    def finite_float_or_none(value):
        if value is None:
            return None
        try:
            converted = float(value)
        except (TypeError, ValueError):
            return None
        return converted if math.isfinite(converted) else None

    transfer_error = finite_float_or_none(boundary_transfer_error_estimate)
    max_transfer_error = finite_float_or_none(max_boundary_transfer_error)

    checks = {
        "parent_model_id_recorded": bool(parent_model_text),
        "submodel_region_id_recorded": bool(submodel_region_text),
        "zoom_boundary_id_recorded": bool(zoom_boundary_text),
        "boundary_transfer_quantity_recorded": bool(transfer_quantity_text),
        "boundary_transfer_error_estimate_recorded": transfer_error is not None,
        "boundary_transfer_error_nonnegative": (
            transfer_error is not None and transfer_error >= 0.0
        ),
        "local_refinement_rule_recorded": bool(refinement_rule_text),
        "target_observable_id_recorded": bool(observable_text),
        "boundary_handoff_not_value_only": all(
            bool(item)
            for item in (
                parent_model_text,
                submodel_region_text,
                zoom_boundary_text,
                transfer_quantity_text,
                refinement_rule_text,
                observable_text,
            )
        ),
    }
    if error_unit_text or transfer_error is not None:
        checks["boundary_transfer_error_unit_recorded"] = bool(error_unit_text)
    if expected_quantity_text is not None:
        checks["boundary_transfer_quantity_matches_expected"] = (
            transfer_quantity_text.lower() == expected_quantity_text.lower()
        )
    if max_boundary_transfer_error is not None:
        checks["max_boundary_transfer_error_recorded"] = (
            max_transfer_error is not None and max_transfer_error >= 0.0
        )
        checks["boundary_transfer_error_within_limit"] = (
            transfer_error is not None
            and max_transfer_error is not None
            and transfer_error <= max_transfer_error
        )
    if parent_mesh_text or local_mesh_text:
        checks["parent_mesh_id_recorded"] = bool(parent_mesh_text)
        checks["local_mesh_id_recorded"] = bool(local_mesh_text)
        checks["parent_local_mesh_identity_separated"] = (
            bool(parent_mesh_text)
            and bool(local_mesh_text)
            and parent_mesh_text != local_mesh_text
        )
    if trace_text:
        checks["boundary_trace_id_recorded"] = True
    if source_text:
        checks["boundary_condition_source_recorded"] = True

    return {
        "policy": "solver_submodel_boundary_handoff_gate",
        "parent_model_id": parent_model_text,
        "submodel_region_id": submodel_region_text,
        "zoom_boundary_id": zoom_boundary_text,
        "boundary_transfer_quantity": transfer_quantity_text,
        "boundary_transfer_error_estimate": transfer_error,
        "boundary_transfer_error_unit": error_unit_text,
        "max_boundary_transfer_error": max_transfer_error,
        "local_refinement_rule": refinement_rule_text,
        "target_observable_id": observable_text,
        "parent_mesh_id": parent_mesh_text,
        "local_mesh_id": local_mesh_text,
        "boundary_trace_id": trace_text,
        "boundary_condition_source": source_text,
        "expected_boundary_transfer_quantity": expected_quantity_text,
        "checks": checks,
        "status": "ok" if all(checks.values()) else "needs_attention",
        "notes": [
            "local refinement and parent-to-local boundary transfer are separate error sources",
            "record the boundary handoff before interpreting the local target observable",
            "a submodel validation row should not be only a refined local value",
        ],
    }


def owned_solver_model_tag_preflight_gate(
    existing_tags,
    requested_tag,
    owned_prefix,
    owner,
    max_tag_length=63,
):
    """Reject unsafe or colliding temporary solver tags before mutation."""

    tags = [str(tag).strip() for tag in existing_tags]
    tag = str(requested_tag or "").strip()
    prefix = str(owned_prefix or "").strip()
    owner_text = str(owner or "").strip()
    limit = int(max_tag_length)
    if limit <= 0:
        raise ValueError("max_tag_length must be positive")
    safe_chars = bool(tag) and all(
        char.isascii() and (char.isalnum() or char == "_") for char in tag
    )
    lower_tags = {item.lower() for item in tags}
    checks = {
        "requested_tag_recorded": bool(tag),
        "owned_prefix_recorded": bool(prefix),
        "owner_recorded": bool(owner_text),
        "requested_tag_uses_owned_prefix": bool(prefix) and tag.startswith(prefix),
        "requested_tag_safe_for_solver": safe_chars and tag[:1].isalpha(),
        "requested_tag_within_length_limit": bool(tag) and len(tag) <= limit,
        "requested_tag_collision_free": bool(tag) and tag.lower() not in lower_tags,
        "existing_tags_unique_case_insensitive": len(lower_tags) == len(tags),
    }
    return {
        "policy": "owned_solver_model_tag_preflight_gate",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "may_create": all(checks.values()),
        "requested_tag": tag or None,
        "owned_prefix": prefix or None,
        "owner": owner_text or None,
        "existing_tag_count": len(tags),
        "collision": bool(tag) and tag.lower() in lower_tags,
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "cleanup_contract": {
            "remove_only_requested_tag": tag or None,
            "preserve_preexisting_tags": True,
            "verify_owned_prefix_absent_after_cleanup": True,
        },
    }


def owned_solver_model_tag_lifecycle_gate(
    artifact,
    expected_artifact_id=None,
    expected_model_tag_prefix=None,
    expected_model_tag_owner=None,
):
    """Check that a temporary solver model tag was owned and cleaned up.

    Long-lived commercial or notebook-backed solver sessions often retain model
    tags between runs.  A validation slot may create a temporary tag for
    introspection or a compact solve, but the result should prove that the tag
    belonged to that slot, was visible while needed, and was removed without
    disturbing pre-existing session tags.  The gate is solver-independent; the
    private source-tool lane records how the tag was created.
    """

    if not isinstance(artifact, dict):
        raise ValueError("artifact must be a mapping")

    def _first(row, names, default=None):
        for name in names:
            if name in row and row[name] is not None:
                return row[name]
        return default

    def _norm(value):
        return str(value or "").strip().lower().replace("-", "_").replace(" ", "_")

    def _string_list(value):
        if value is None:
            return []
        if isinstance(value, str):
            return [value]
        if isinstance(value, (list, tuple)):
            return [str(item) for item in value]
        return [str(value)]

    def _nonnegative_int_or_none(value):
        if value is None:
            return None
        try:
            converted = int(value)
        except (TypeError, ValueError):
            return None
        return converted if converted >= 0 else None

    status_text = _norm(_first(artifact, ("status", "session_status", "livelink_status"), ""))
    artifact_id = str(_first(artifact, ("artifact_id", "model_artifact_id", "result_artifact_id"), "")).strip()
    model_tag = str(_first(artifact, ("model_tag", "tag", "temporary_model_tag"), "")).strip()
    owner = str(_first(artifact, ("model_tag_owner", "tag_owner", "owner"), "")).strip()
    owned_prefix = str(
        expected_model_tag_prefix
        if expected_model_tag_prefix is not None
        else _first(artifact, ("owned_model_tag_prefix", "model_tag_prefix", "tag_prefix"), "")
    ).strip()
    tags_before = _string_list(_first(artifact, ("tags_before", "model_tags_before"), []))
    tags_after = _string_list(_first(artifact, ("tags_after", "model_tags_after"), []))
    before_count = _nonnegative_int_or_none(
        _first(artifact, ("tags_before_count", "model_tags_before_count"))
    )
    created_count = _nonnegative_int_or_none(
        _first(artifact, ("tags_created_count", "model_tags_created_count"))
    )
    after_count = _nonnegative_int_or_none(
        _first(artifact, ("tags_after_count", "model_tags_after_count"))
    )
    pre_existing = bool(_first(artifact, ("pre_existing_owned_tag", "tag_pre_existing"), False))
    created_present = bool(_first(artifact, ("created_present", "tag_created_present"), False))
    removed_after = bool(_first(artifact, ("removed_after_probe", "tag_removed_after_probe"), False))
    preserved = bool(
        _first(artifact, ("preexisting_tags_preserved", "existing_tags_preserved"), False)
    )
    connected = bool(_first(artifact, ("connected", "session_connected"), False))
    owned_prefix_tags_after = [
        tag for tag in tags_after if owned_prefix and tag.startswith(owned_prefix)
    ]

    checks = {
        "session_connected": connected,
        "status_allows_reuse": status_text in {"already_connected", "connected", "live", "ok"},
        "artifact_id_recorded": bool(artifact_id),
        "model_tag_recorded": bool(model_tag),
        "model_tag_owner_recorded": bool(owner),
        "started_no_new_process": not bool(artifact.get("started_new_process", False))
        and not bool(artifact.get("started_new_matlab", False))
        and not bool(artifact.get("started_new_comsol", False)),
        "killed_no_process": not bool(artifact.get("killed_process", False)),
        "tag_was_created_and_visible": created_present,
        "owned_tag_removed_after_probe": removed_after,
        "preexisting_tags_preserved": preserved,
        "tag_absent_after_cleanup": bool(model_tag) and model_tag not in set(tags_after),
        "tag_counts_recorded": before_count is not None and created_count is not None and after_count is not None,
    }
    if checks["tag_counts_recorded"]:
        checks["created_count_increased_by_one"] = (
            created_count == before_count + 1 if not pre_existing else created_count >= before_count
        )
        checks["after_count_restored"] = after_count == before_count
    if expected_artifact_id is not None:
        checks["expected_artifact_id_matches"] = artifact_id == str(expected_artifact_id)
    if owned_prefix:
        checks["owned_model_tag_prefix_recorded"] = True
        checks["expected_model_tag_prefix_matches"] = model_tag.startswith(owned_prefix)
        checks["owned_prefix_absent_after_cleanup"] = not owned_prefix_tags_after
    if expected_model_tag_owner is not None:
        checks["expected_model_tag_owner_matches"] = owner == str(expected_model_tag_owner)

    return {
        "policy": "owned_solver_model_tag_lifecycle_gate",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "artifact_id": artifact_id or None,
        "model_tag": model_tag or None,
        "owned_model_tag_prefix": owned_prefix or None,
        "model_tag_owner": owner or None,
        "connected": connected,
        "session_status": status_text,
        "pre_existing_owned_tag": pre_existing,
        "created_present": created_present,
        "removed_after_probe": removed_after,
        "preexisting_tags_preserved": preserved,
        "tags_before_count": before_count,
        "tags_created_count": created_count,
        "tags_after_count": after_count,
        "tags_before": tags_before,
        "tags_after": tags_after,
        "owned_prefix_tags_after": owned_prefix_tags_after,
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "notes": [
            "temporary solver-owned tags are provenance, not physics residuals",
            "create/remove only the owned tag and preserve pre-existing session tags",
            "use a unique owned prefix for slot-created tags and verify that prefix leaves no residue",
            "bind model_tag_owner and artifact_id before trusting a model-tag probe",
        ],
    }


def solver_result_artifact_provenance_timing_gate(
    artifact,
    required_versions=("solver",),
    required_timing_stages=(),
    min_timing_stages=1,
    max_timing_stages=4,
    max_missing_timing_stages=0,
    expected_created_at_utc=None,
    expected_run_date_utc=None,
    max_created_run_skew_s=None,
    require_run_date=True,
    expected_execution_session_id=None,
    require_execution_session_id=False,
    expected_result_output_schema_id=None,
    expected_result_output_columns=None,
    expected_result_output_units=None,
    require_result_output_schema=False,
):
    """Check that a solver result JSON records provenance and timing.

    Result numbers are much easier to reuse when the artifact records when it
    was run, which solver/tool versions produced it, and where the wall time
    went.  This gate is intentionally solver-independent; private COMSOL,
    MATLAB, CST, or other lanes can keep their provenance private while public
    helpers replay the structural contract.
    """

    if not isinstance(artifact, dict):
        raise ValueError("artifact must be a mapping")

    def parse_datetime(value):
        text = str(value or "").strip()
        if not text:
            return None
        try:
            return datetime.fromisoformat(text.replace("Z", "+00:00"))
        except ValueError:
            return None

    def parse_time(value):
        return parse_datetime(value) is not None

    def string_list(value):
        if value is None:
            return ()
        if isinstance(value, str):
            text = value.strip()
            if not text:
                return ()
            if "," in text or ";" in text:
                return tuple(
                    part.strip()
                    for part in text.replace(";", ",").split(",")
                    if part.strip()
                )
            return (text,)
        if isinstance(value, (list, tuple)):
            return tuple(str(item).strip() for item in value if str(item).strip())
        return (str(value).strip(),) if str(value).strip() else ()

    def unit_mapping(value):
        if value is None:
            return ()
        if isinstance(value, str):
            text = value.strip()
            if not text:
                return ()
            try:
                import json

                parsed = json.loads(text)
            except (TypeError, ValueError):
                parsed = None
            if isinstance(parsed, dict):
                value = parsed
            else:
                pairs = []
                for part in text.replace(";", ",").split(","):
                    if not part.strip():
                        continue
                    if ":" in part:
                        key, unit = part.split(":", 1)
                    elif "=" in part:
                        key, unit = part.split("=", 1)
                    else:
                        continue
                    pairs.append((key.strip(), unit.strip()))
                return tuple(sorted((key, unit) for key, unit in pairs if key))
        if isinstance(value, dict):
            return tuple(
                sorted(
                    (str(key).strip(), str(unit).strip())
                    for key, unit in value.items()
                    if str(key).strip()
                )
            )
        return ()

    def first_from(records, names):
        for record in records:
            if not isinstance(record, dict):
                continue
            for name in names:
                if name in record and record[name] is not None:
                    return record[name]
        return None

    versions = artifact.get("versions", {})
    if versions is None:
        versions = {}
    if not isinstance(versions, dict):
        raise ValueError("artifact['versions'] must be a mapping when present")
    required_version_keys = [str(key).strip() for key in required_versions if str(key).strip()]

    timing = artifact.get("timing_breakdown_s", {})
    if timing is None:
        timing = {}
    if not isinstance(timing, dict):
        raise ValueError("artifact['timing_breakdown_s'] must be a mapping when present")
    required_stages = [str(key).strip() for key in required_timing_stages if str(key).strip()]
    min_stages = int(min_timing_stages)
    max_stages = int(max_timing_stages) if max_timing_stages is not None else 0
    missing_allowed = int(max_missing_timing_stages)
    max_skew = None
    if max_created_run_skew_s is not None:
        max_skew = float(max_created_run_skew_s)
    if min_stages < 0:
        raise ValueError("min_timing_stages must be non-negative")
    if max_stages < 0:
        raise ValueError("max_timing_stages must be non-negative")
    if missing_allowed < 0:
        raise ValueError("max_missing_timing_stages must be non-negative")
    if max_skew is not None and max_skew < 0.0:
        raise ValueError("max_created_run_skew_s must be non-negative")

    timing_values = {}
    timing_values_ok = True
    for key, value in timing.items():
        try:
            numeric = float(value)
        except (TypeError, ValueError):
            timing_values_ok = False
            continue
        if not math.isfinite(numeric) or numeric < 0.0:
            timing_values_ok = False
        timing_values[str(key)] = numeric

    missing_versions = [
        key for key in required_version_keys
        if not str(versions.get(key, "")).strip()
    ]
    missing_stages = [
        key for key in required_stages
        if key not in timing_values
    ]
    dominant_stages = [
        {"stage": key, "seconds": value}
        for key, value in sorted(timing_values.items(), key=lambda item: item[1], reverse=True)[:4]
    ]
    created_at = artifact.get("created_at_utc", artifact.get("created_at", ""))
    execution = artifact.get("execution", {})
    if execution is None:
        execution = {}
    if not isinstance(execution, dict):
        raise ValueError("artifact['execution'] must be a mapping when present")
    result_block = artifact.get("result", {})
    if not isinstance(result_block, dict):
        result_block = {}
    output_block = artifact.get("output", artifact.get("outputs", {}))
    if not isinstance(output_block, dict):
        output_block = {}
    run_date = execution.get("run_date_utc", artifact.get("run_date_utc", created_at))
    execution_session_id = (
        execution.get(
            "execution_session_id",
            execution.get(
                "session_id",
                execution.get(
                    "shared_engine",
                    execution.get(
                        "matlab_engine_session",
                        execution.get("engine_session", artifact.get("execution_session_id", "")),
                    ),
                ),
            ),
        )
    )
    created_text = str(created_at or "").strip()
    run_date_text = str(run_date or "").strip()
    execution_session_text = str(execution_session_id or "").strip()
    run_date_recorded = bool(
        str(execution.get("run_date_utc", "")).strip()
        or str(artifact.get("run_date_utc", "")).strip()
    )
    created_dt = parse_datetime(created_at)
    run_dt = parse_datetime(run_date)
    created_run_skew_s = None
    if created_dt is not None and run_dt is not None:
        created_run_skew_s = abs((run_dt - created_dt).total_seconds())
    expected_created = (
        None
        if expected_created_at_utc is None
        else str(expected_created_at_utc).strip()
    )
    expected_run_date = (
        None
        if expected_run_date_utc is None
        else str(expected_run_date_utc).strip()
    )
    expected_execution_session = (
        None
        if expected_execution_session_id is None
        else str(expected_execution_session_id).strip()
    )
    result_records = (artifact, execution, result_block, output_block)
    result_output_schema_id = first_from(
        result_records,
        (
            "result_output_schema_id",
            "resultOutputSchemaId",
            "output_schema_id",
            "outputSchemaId",
            "table_schema_id",
            "tableSchemaId",
        ),
    )
    result_output_columns = string_list(
        first_from(
            result_records,
            (
                "result_output_columns",
                "resultOutputColumns",
                "output_columns",
                "outputColumns",
                "table_columns",
                "tableColumns",
                "columns",
            ),
        )
    )
    result_output_units = unit_mapping(
        first_from(
            result_records,
            (
                "result_output_units",
                "resultOutputUnits",
                "output_units",
                "outputUnits",
                "table_units",
                "tableUnits",
                "column_units",
                "columnUnits",
                "units",
            ),
        )
    )
    result_output_schema_text = str(result_output_schema_id or "").strip()
    expected_result_output_schema = (
        None
        if expected_result_output_schema_id is None
        else str(expected_result_output_schema_id).strip()
    )
    expected_result_output_columns_tuple = (
        None
        if expected_result_output_columns is None
        else string_list(expected_result_output_columns)
    )
    expected_result_output_units_tuple = (
        None
        if expected_result_output_units is None
        else unit_mapping(expected_result_output_units)
    )
    result_output_schema_required = (
        bool(require_result_output_schema)
        or expected_result_output_schema is not None
        or expected_result_output_columns_tuple is not None
        or expected_result_output_units_tuple is not None
    )

    checks = {
        "schema_recorded": bool(str(artifact.get("schema", "")).strip()),
        "created_at_utc_recorded": bool(created_text),
        "created_at_utc_parseable": parse_time(created_at),
        "run_date_utc_recorded_when_required": (
            not bool(require_run_date) or run_date_recorded
        ),
        "run_date_utc_parseable": parse_time(run_date),
        "expected_created_at_utc_matches": (
            expected_created is None or created_text == expected_created
        ),
        "expected_run_date_utc_matches": (
            expected_run_date is None or run_date_text == expected_run_date
        ),
        "execution_session_id_recorded_when_required": (
            not bool(require_execution_session_id) or bool(execution_session_text)
        ),
        "expected_execution_session_id_matches": (
            expected_execution_session is None
            or execution_session_text == expected_execution_session
        ),
        "created_run_timestamp_skew_within_limit": (
            max_skew is None
            or (
                created_run_skew_s is not None
                and created_run_skew_s <= max_skew
            )
        ),
        "versions_mapping_recorded": bool(versions),
        "required_versions_recorded": not missing_versions,
        "timing_breakdown_recorded": len(timing_values) >= min_stages,
        "timing_stage_count_reasonable": max_stages == 0 or len(timing_values) <= max_stages,
        "timing_values_nonnegative": timing_values_ok,
        "required_timing_stages_recorded": len(missing_stages) <= missing_allowed,
        "dominant_timing_stages_available": len(dominant_stages) > 0,
        "result_output_schema_id_recorded_when_required": (
            not result_output_schema_required or bool(result_output_schema_text)
        ),
        "result_output_columns_recorded_when_required": (
            not result_output_schema_required or bool(result_output_columns)
        ),
        "result_output_units_recorded_when_required": (
            not result_output_schema_required or bool(result_output_units)
        ),
        "expected_result_output_schema_id_matches": (
            expected_result_output_schema is None
            or result_output_schema_text == expected_result_output_schema
        ),
        "expected_result_output_columns_match": (
            expected_result_output_columns_tuple is None
            or result_output_columns == expected_result_output_columns_tuple
        ),
        "expected_result_output_units_match": (
            expected_result_output_units_tuple is None
            or result_output_units == expected_result_output_units_tuple
        ),
    }

    return {
        "policy": "solver_result_artifact_provenance_timing_gate",
        "schema": artifact.get("schema"),
        "created_at_utc": created_at,
        "run_date_utc": run_date,
        "expected_created_at_utc": expected_created,
        "expected_run_date_utc": expected_run_date,
        "require_run_date": bool(require_run_date),
        "execution_session_id": execution_session_id,
        "expected_execution_session_id": expected_execution_session,
        "require_execution_session_id": bool(require_execution_session_id),
        "result_output_schema_id": result_output_schema_text,
        "result_output_columns": list(result_output_columns),
        "result_output_units": dict(result_output_units),
        "expected_result_output_schema_id": expected_result_output_schema,
        "expected_result_output_columns": (
            None
            if expected_result_output_columns_tuple is None
            else list(expected_result_output_columns_tuple)
        ),
        "expected_result_output_units": (
            None
            if expected_result_output_units_tuple is None
            else dict(expected_result_output_units_tuple)
        ),
        "require_result_output_schema": bool(require_result_output_schema),
        "created_run_skew_s": created_run_skew_s,
        "max_created_run_skew_s": max_skew,
        "versions": versions,
        "required_versions": required_version_keys,
        "missing_versions": missing_versions,
        "timing_breakdown_s": timing_values,
        "required_timing_stages": required_stages,
        "missing_timing_stages": missing_stages,
        "max_timing_stages": max_stages,
        "timing_stage_count": len(timing_values),
        "dominant_timing_stages": dominant_stages,
        "total_recorded_timing_s": sum(timing_values.values()),
        "checks": checks,
        "status": "ok" if all(checks.values()) else "needs_attention",
        "notes": [
            "record solver/tool versions before result values are reused",
            "record run date and created_at in parseable ISO form",
            "record the execution session id when reused results depend on an existing live solver or MATLAB session",
            "record result output schema id, columns, and units before importing saved JSON/table rows into notebooks",
            "keep created_at and run_date close enough to describe the same executed result artifact",
            "keep the heaviest timing stages visible so slow slots can be compared later",
        ],
    }


def source_native_seed_queue_gate(
    queue_artifact,
    expected_tools=(),
    expected_rounds=None,
    expected_total_slots=None,
    require_all_local_present=True,
    require_public_safe_sources=False,
    allow_verified_slots=False,
):
    """Gate source-native loop seed queues before calling them learned.

    A source-native seed queue is useful learning material, but it is not a
    solver result by itself.  This gate keeps that distinction explicit: queued
    examples can be accepted as replay seeds while solver/cross-validation
    claims must wait for a promoted result artifact.
    """

    if not isinstance(queue_artifact, dict):
        raise ValueError("queue_artifact must be a mapping")

    def as_int(value):
        try:
            return int(value)
        except (TypeError, ValueError):
            return None

    def clean_text(value):
        return str(value or "").strip()

    def list_tools(value):
        if value is None:
            return []
        if isinstance(value, str):
            return [item.strip() for item in value.split(",") if item.strip()]
        if isinstance(value, (list, tuple, set)):
            return [clean_text(item) for item in value if clean_text(item)]
        return [clean_text(value)] if clean_text(value) else []

    def looks_private_source(value):
        text = clean_text(value)
        if not text:
            return False
        normalized = text.replace("\\", "/").lower()
        if len(text) >= 3 and text[1] == ":" and text[2] in ("\\", "/"):
            return True
        return any(
            marker in normalized
            for marker in (
                "_cross" + "val",
                "internal://",
                "private://",
                "lab_private",
                "unpublished",
            )
        )

    valid_source_types = {
        "local_path",
        "local_project",
        "public_url",
        "public_doc",
        "upstream_example",
        "training_project",
        "source_native",
        "repository",
        "manual_example",
    }
    required_slot_fields = (
        "tool",
        "source_native_example",
        "source_type",
        "lesson_axis",
        "intended_validation",
        "status",
    )
    solver_claim_tokens = (
        "verified",
        "learned",
        "solver_pass",
        "solver_verified",
        "crossval_passed",
        "validated",
    )

    slots = queue_artifact.get("slots")
    if not isinstance(slots, list):
        slots = []

    expected_tool_names = list_tools(expected_tools)
    actual_tools = []
    tool_counts = {}
    missing_fields = []
    invalid_source_types = []
    missing_local_sources = []
    solver_claim_slots = []
    private_source_slots = []
    slot_ids = []
    laps = []

    for index, slot in enumerate(slots):
        if not isinstance(slot, dict):
            missing_fields.append(
                {
                    "index": index,
                    "slot_id": "",
                    "tool": "",
                    "missing": list(required_slot_fields),
                }
            )
            continue
        tool = clean_text(slot.get("tool"))
        actual_tools.append(tool)
        if tool:
            tool_counts[tool] = tool_counts.get(tool, 0) + 1
        slot_id = clean_text(slot.get("slot_id") or slot.get("id") or index)
        slot_ids.append(slot_id)
        lap = as_int(slot.get("lap") or slot.get("round"))
        if lap is not None:
            laps.append(lap)

        missing = [field for field in required_slot_fields if not clean_text(slot.get(field))]
        if missing:
            missing_fields.append(
                {
                    "index": index,
                    "slot_id": slot_id,
                    "tool": tool,
                    "missing": missing,
                }
            )

        source_type = clean_text(slot.get("source_type"))
        if source_type and source_type not in valid_source_types:
            invalid_source_types.append(
                {
                    "index": index,
                    "slot_id": slot_id,
                    "tool": tool,
                    "source_type": source_type,
                }
            )

        if (
            bool(require_all_local_present)
            and source_type in {"local_path", "local_project"}
            and slot.get("local_exists") is False
        ):
            missing_local_sources.append(
                {
                    "index": index,
                    "slot_id": slot_id,
                    "tool": tool,
                    "source_native_example": clean_text(slot.get("source_native_example")),
                }
            )

        status_text = clean_text(slot.get("status")).lower()
        if (
            not bool(allow_verified_slots)
            and any(token in status_text for token in solver_claim_tokens)
        ):
            solver_claim_slots.append(
                {
                    "index": index,
                    "slot_id": slot_id,
                    "tool": tool,
                    "status": status_text,
                }
            )

        if bool(require_public_safe_sources) and looks_private_source(
            slot.get("source_native_example")
        ):
            private_source_slots.append(
                {
                    "index": index,
                    "slot_id": slot_id,
                    "tool": tool,
                    "source_type": source_type,
                }
            )

    expected_round_count = as_int(expected_rounds)
    expected_slot_count = as_int(expected_total_slots)
    declared_rounds = as_int(queue_artifact.get("rounds"))
    declared_total_slots = as_int(queue_artifact.get("total_slots"))
    actual_rounds = max(laps) if laps else 0
    missing_expected_tools = [
        tool for tool in expected_tool_names if tool not in set(actual_tools)
    ]
    unexpected_tools = [
        tool for tool in sorted(set(actual_tools)) if expected_tool_names and tool not in expected_tool_names
    ]
    duplicate_slot_ids = sorted(
        slot_id
        for slot_id in set(slot_ids)
        if slot_id and slot_ids.count(slot_id) > 1
    )

    checks = {
        "created_at_recorded": bool(
            clean_text(queue_artifact.get("created_at_utc"))
            or clean_text(queue_artifact.get("created_at"))
        ),
        "slots_list_recorded": isinstance(queue_artifact.get("slots"), list),
        "slots_nonempty": len(slots) > 0,
        "declared_total_matches_slot_count": (
            declared_total_slots is None or declared_total_slots == len(slots)
        ),
        "expected_total_matches_slot_count": (
            expected_slot_count is None or expected_slot_count == len(slots)
        ),
        "declared_rounds_match_slots": (
            declared_rounds is None or actual_rounds == 0 or declared_rounds == actual_rounds
        ),
        "expected_rounds_match_slots": (
            expected_round_count is None
            or actual_rounds == 0
            or expected_round_count == actual_rounds
        ),
        "expected_tools_present": not missing_expected_tools,
        "no_unexpected_tools": not unexpected_tools,
        "required_slot_fields_present": not missing_fields,
        "source_types_allowed": not invalid_source_types,
        "local_sources_exist_when_required": not missing_local_sources,
        "no_solver_or_learning_overclaim": not solver_claim_slots,
        "public_safe_sources_when_required": not private_source_slots,
        "slot_ids_unique_when_recorded": not duplicate_slot_ids,
    }

    ok = all(checks.values())
    return {
        "policy": "source_native_seed_queue_gate",
        "status": "ok" if ok else "needs_attention",
        "learning_stage": "queued_not_learned" if ok else "queue_needs_attention",
        "slot_count": len(slots),
        "declared_total_slots": declared_total_slots,
        "expected_total_slots": expected_slot_count,
        "declared_rounds": declared_rounds,
        "actual_rounds": actual_rounds,
        "expected_rounds": expected_round_count,
        "expected_tools": expected_tool_names,
        "actual_tools": sorted(set(actual_tools)),
        "tool_counts": tool_counts,
        "missing_expected_tools": missing_expected_tools,
        "unexpected_tools": unexpected_tools,
        "missing_fields": missing_fields,
        "invalid_source_types": invalid_source_types,
        "missing_local_sources": missing_local_sources,
        "solver_claim_slots": solver_claim_slots,
        "private_source_slots": private_source_slots,
        "duplicate_slot_ids": duplicate_slot_ids,
        "checks": checks,
        "notes": [
            "This gate verifies replay seeds, not solver results.",
            "Keep seed queues in candidate/queued state until a result artifact passes the MCP feedback gate.",
            "Use require_public_safe_sources=True only for scrubbed public artifacts; private loop queues may contain internal source paths.",
        ],
    }


def computed_reference_rows_gate(
    artifact,
    rows=None,
    require_pass=True,
    max_global_rel_error=None,
    rel_error_atol=1.0e-12,
):
    """Check actual computed/reference cross-validation rows.

    The feedback gate verifies provenance and MCP closure.  This row gate checks
    the numerical table itself: every row must expose a computed value, a
    reference value, a tolerance, and a relative error that is within tolerance.
    """

    if rows is None:
        if not isinstance(artifact, dict):
            raise ValueError("artifact must be a mapping when rows is omitted")
        candidates = []
        for key in ("rows", "table", "checks"):
            value = artifact.get(key)
            if isinstance(value, list):
                candidates.extend(value)
        results = artifact.get("results")
        if isinstance(results, list):
            for result in results:
                if not isinstance(result, dict):
                    continue
                case_name = str(result.get("name") or result.get("case") or "").strip()
                for check in result.get("checks", []):
                    if isinstance(check, dict):
                        row = dict(check)
                        if case_name and not row.get("case"):
                            row["case"] = case_name
                        candidates.append(row)
        rows = candidates

    if rows is None:
        rows = []
    if not isinstance(rows, list):
        raise ValueError("rows must be a list")

    def as_float(value):
        try:
            return float(value)
        except (TypeError, ValueError):
            return None

    def calc_rel(value, reference):
        return abs(value - reference) / max(abs(reference), 1.0e-300)

    required_fields = ("quantity", "computed", "reference", "tolerance")
    normalized_rows = []
    missing_fields = []
    nonnumeric_rows = []
    rel_error_mismatch_rows = []
    tolerance_fail_rows = []
    pass_flag_fail_rows = []

    for index, raw in enumerate(rows):
        if not isinstance(raw, dict):
            missing_fields.append(
                {
                    "index": index,
                    "case": "",
                    "quantity": "",
                    "missing": list(required_fields),
                }
            )
            continue
        case_name = str(raw.get("case") or "").strip()
        quantity = str(raw.get("quantity") or raw.get("name") or "").strip()
        missing = [field for field in required_fields if field not in raw or raw.get(field) in (None, "")]
        if missing:
            missing_fields.append(
                {
                    "index": index,
                    "case": case_name,
                    "quantity": quantity,
                    "missing": missing,
                }
            )
            continue

        computed = as_float(raw.get("computed"))
        reference = as_float(raw.get("reference"))
        tolerance = as_float(raw.get("tolerance"))
        recorded_rel = as_float(raw.get("rel_error"))
        if computed is None or reference is None or tolerance is None:
            nonnumeric_rows.append(
                {
                    "index": index,
                    "case": case_name,
                    "quantity": quantity,
                    "computed": raw.get("computed"),
                    "reference": raw.get("reference"),
                    "tolerance": raw.get("tolerance"),
                }
            )
            continue

        computed_rel = calc_rel(computed, reference)
        effective_rel = computed_rel if recorded_rel is None else recorded_rel
        if recorded_rel is not None:
            mismatch = abs(recorded_rel - computed_rel)
            if mismatch > float(rel_error_atol) * max(1.0, abs(computed_rel)):
                rel_error_mismatch_rows.append(
                    {
                        "index": index,
                        "case": case_name,
                        "quantity": quantity,
                        "recorded_rel_error": recorded_rel,
                        "computed_rel_error": computed_rel,
                    }
                )
        row_pass = bool(raw.get("pass")) if "pass" in raw else effective_rel <= tolerance
        if effective_rel > tolerance:
            tolerance_fail_rows.append(
                {
                    "index": index,
                    "case": case_name,
                    "quantity": quantity,
                    "rel_error": effective_rel,
                    "tolerance": tolerance,
                }
            )
        if bool(require_pass) and not row_pass:
            pass_flag_fail_rows.append(
                {
                    "index": index,
                    "case": case_name,
                    "quantity": quantity,
                    "pass": raw.get("pass"),
                }
            )
        normalized_rows.append(
            {
                "index": index,
                "case": case_name,
                "quantity": quantity,
                "computed": computed,
                "reference": reference,
                "rel_error": effective_rel,
                "computed_rel_error": computed_rel,
                "tolerance": tolerance,
                "pass": row_pass,
            }
        )

    max_rel = max((row["rel_error"] for row in normalized_rows), default=None)
    global_limit = as_float(max_global_rel_error)
    checks = {
        "rows_recorded": len(rows) > 0,
        "required_fields_present": not missing_fields,
        "numeric_values_recorded": not nonnumeric_rows,
        "rel_error_matches_computed_reference": not rel_error_mismatch_rows,
        "row_errors_within_tolerance": not tolerance_fail_rows,
        "pass_flags_true_when_required": not pass_flag_fail_rows,
        "max_global_rel_error_within_limit": (
            global_limit is None
            or max_rel is not None
            and max_rel <= global_limit
        ),
    }
    return {
        "policy": "computed_reference_rows_gate",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "row_count": len(rows),
        "valid_row_count": len(normalized_rows),
        "max_rel_error": max_rel,
        "max_global_rel_error": global_limit,
        "missing_fields": missing_fields,
        "nonnumeric_rows": nonnumeric_rows,
        "rel_error_mismatch_rows": rel_error_mismatch_rows,
        "tolerance_fail_rows": tolerance_fail_rows,
        "pass_flag_fail_rows": pass_flag_fail_rows,
        "rows": normalized_rows,
        "checks": checks,
        "notes": [
            "Use this after a real solver run has produced computed/reference rows.",
            "Pair this row gate with the artifact feedback gate so numerical validity and MCP closure are both checked.",
        ],
    }


def cross_validation_artifact_to_mcp_feedback_gate(
    artifact,
    expected_public_status="verified",
    require_pass=True,
    require_public_lesson=True,
    require_learning_target=True,
    require_verification=True,
    require_result_provenance=True,
    require_result_artifact=True,
    require_notebook_source=False,
    require_result_output_schema=True,
    required_versions=("solver", "radia_mcp"),
    min_timing_stages=1,
    max_timing_stages=4,
    require_replayable_verification_commands=False,
):
    """Check that a cross-validation or notebook result was fed back to MCP.

    This gate sits above the result provenance/table gates.  It answers a
    workflow question: did the artifact become durable MCP knowledge, or was it
    merely collected as a JSON/notebook output?
    """

    if not isinstance(artifact, dict):
        raise ValueError("artifact must be a mapping")

    def first_from(records, names, default=""):
        for record in records:
            if not isinstance(record, dict):
                continue
            for name in names:
                value = record.get(name)
                if value is not None and str(value).strip():
                    return value
        return default

    def string_list(value):
        if value is None:
            return []
        if isinstance(value, str):
            text = value.strip()
            if not text:
                return []
            if "," in text or ";" in text:
                return [
                    item.strip()
                    for item in text.replace(";", ",").split(",")
                    if item.strip()
                ]
            return [text]
        if isinstance(value, (list, tuple, set)):
            return [str(item).strip() for item in value if str(item).strip()]
        return [str(value).strip()] if str(value).strip() else []

    def rank_status(value):
        text = str(value or "").strip().lower()
        ranks = {
            "": -1,
            "none": 0,
            "candidate": 1,
            "collected": 1,
            "distilled": 1,
            "encoded": 2,
            "verified": 3,
            "learned": 3,
        }
        return ranks.get(text, -1), text

    def replay_command_candidates(record):
        if not isinstance(record, dict):
            return []
        commands = []
        nested = record.get("commands")
        if isinstance(nested, (list, tuple)):
            for item in nested:
                if isinstance(item, dict):
                    value = item.get("command")
                    if value is not None and str(value).strip():
                        commands.append(str(value).strip())
                elif str(item).strip():
                    commands.append(str(item).strip())
        command = record.get("command")
        if command is not None and str(command).strip():
            commands.append(str(command).strip())
        return commands

    def is_normalized_replay_command(command):
        text = str(command or "").strip()
        if not text:
            return False
        if "->" in text:
            return False
        return any(
            text.startswith(prefix)
            for prefix in (
                "python ",
                "python -m ",
                "pytest ",
                "pwsh ",
                "powershell ",
                "git ",
                "matlab -batch",
            )
        )

    learning_lanes = artifact.get("learning_lanes", {})
    if not isinstance(learning_lanes, dict):
        learning_lanes = {}
    verification = artifact.get("verification", {})
    if not isinstance(verification, dict):
        verification = {}
    mcp_feedback = artifact.get("mcp_feedback", {})
    if not isinstance(mcp_feedback, dict):
        mcp_feedback = {}

    public_status = learning_lanes.get(
        "public",
        learning_lanes.get("open", learning_lanes.get("radia_mcp", "")),
    )
    source_tool_status = learning_lanes.get(
        "source_tool",
        learning_lanes.get("source", learning_lanes.get("private", "")),
    )
    public_rank, public_status_text = rank_status(public_status)
    expected_public_rank, expected_public_status_text = rank_status(expected_public_status)

    learning_targets = [
        *string_list(artifact.get("learning_targets")),
        *string_list(mcp_feedback.get("encoded_targets")),
        *string_list(mcp_feedback.get("knowledge_topics")),
    ]
    public_target_recorded = any(
        "radia-mcp" in target.lower()
        or "radia_ngsolve" in target.lower()
        or "radia-ngsolve" in target.lower()
        for target in learning_targets
    )
    public_verification = first_from(
        (verification, mcp_feedback),
        (
            "public",
            "open",
            "radia_mcp",
            "radia-mcp",
            "radia_ngsolve",
            "radia-ngsolve",
            "test",
            "tests",
        ),
    )
    replay_commands = [
        *replay_command_candidates(verification),
        *replay_command_candidates(mcp_feedback),
    ]
    normalized_replay_commands = [
        command for command in replay_commands if is_normalized_replay_command(command)
    ]
    public_lesson = first_from(
        (artifact, mcp_feedback),
        (
            "public_lesson",
            "public_summary",
            "teaching_note",
            "knowledge_note",
            "lesson",
        ),
    )
    if not public_lesson:
        notes = artifact.get("notes", [])
        note_list = string_list(notes)
        public_lesson = note_list[0] if note_list else ""

    execution = artifact.get("execution", {})
    if not isinstance(execution, dict):
        execution = {}
    result_block = artifact.get("result", {})
    if not isinstance(result_block, dict):
        result_block = {}
    notebook_block = artifact.get("notebook", {})
    if not isinstance(notebook_block, dict):
        notebook_block = {}
    result_records = (artifact, execution, result_block, notebook_block)
    result_artifact_id = first_from(
        result_records,
        (
            "result_artifact_id",
            "execution_result_artifact_id",
            "notebook_result_artifact_id",
            "run_result_artifact_id",
        ),
    )
    notebook_source_artifact_id = first_from(
        result_records,
        (
            "notebook_source_artifact_id",
            "notebook_artifact_id",
            "notebook_id",
            "source_notebook_artifact_id",
        ),
    )
    notebook_source_digest = first_from(
        result_records,
        (
            "notebook_source_digest",
            "notebook_digest",
            "source_notebook_digest",
        ),
    )
    notebook_source_path = first_from(
        result_records,
        (
            "notebook_source_path",
            "notebook_path",
            "source_notebook_path",
        ),
    )

    provenance_gate = solver_result_artifact_provenance_timing_gate(
        artifact,
        required_versions=required_versions,
        min_timing_stages=min_timing_stages,
        max_timing_stages=max_timing_stages,
        require_run_date=True,
        require_result_output_schema=require_result_output_schema,
    )

    checks = {
        "pass_recorded_true_when_required": (
            not bool(require_pass) or artifact.get("pass") is True
        ),
        "public_lane_status_at_expected_level": (
            expected_public_rank < 0 or public_rank >= expected_public_rank
        ),
        "public_lesson_recorded_when_required": (
            not bool(require_public_lesson) or bool(str(public_lesson).strip())
        ),
        "public_learning_target_recorded_when_required": (
            not bool(require_learning_target) or public_target_recorded
        ),
        "public_verification_recorded_when_required": (
            not bool(require_verification) or bool(str(public_verification).strip())
        ),
        "result_provenance_gate_ok_when_required": (
            not bool(require_result_provenance) or provenance_gate["status"] == "ok"
        ),
        "result_artifact_id_recorded_when_required": (
            not bool(require_result_artifact) or bool(str(result_artifact_id).strip())
        ),
        "notebook_source_artifact_id_recorded_when_required": (
            not bool(require_notebook_source) or bool(str(notebook_source_artifact_id).strip())
        ),
        "notebook_source_digest_recorded_when_required": (
            not bool(require_notebook_source) or bool(str(notebook_source_digest).strip())
        ),
        "notebook_source_path_recorded_when_required": (
            not bool(require_notebook_source) or bool(str(notebook_source_path).strip())
        ),
        "replay_command_recorded_when_required": (
            not bool(require_replayable_verification_commands) or bool(replay_commands)
        ),
        "replay_commands_normalized_when_required": (
            not bool(require_replayable_verification_commands)
            or bool(replay_commands)
            and len(normalized_replay_commands) == len(replay_commands)
        ),
    }

    learned = all(checks.values())
    return {
        "policy": "cross_validation_artifact_to_mcp_feedback_gate",
        "status": "ok" if learned else "needs_attention",
        "learning_stage": "learned" if learned else "collected_or_encoded",
        "expected_public_status": expected_public_status_text,
        "public_lane_status": public_status_text,
        "source_tool_lane_status": str(source_tool_status or "").strip().lower(),
        "public_lesson": str(public_lesson or ""),
        "learning_targets": learning_targets,
        "public_verification": str(public_verification or ""),
        "result_artifact_id": str(result_artifact_id or ""),
        "notebook_source_artifact_id": str(notebook_source_artifact_id or ""),
        "notebook_source_digest": str(notebook_source_digest or ""),
        "notebook_source_path": str(notebook_source_path or ""),
        "replay_commands": replay_commands,
        "normalized_replay_commands": normalized_replay_commands,
        "provenance_gate_status": provenance_gate["status"],
        "provenance_gate_policy": provenance_gate["policy"],
        "dominant_timing_stages": provenance_gate.get("dominant_timing_stages", []),
        "checks": checks,
        "notes": [
            "Use this after a JSON, notebook, or cross-validation result exists.",
            "The artifact is learned only after a public-safe lesson, MCP target, and focused verification are recorded.",
            "Store replayable commands separately from human result notes; keep annotations like '-> passed' in result fields.",
            "Pair this feedback gate with the solver/result provenance and table metadata gates before notebook reuse.",
            "Keep private source-tool provenance in the owning private lane; public radia-mcp records only the scrubbed lesson and gate.",
        ],
    }


def solver_result_table_metadata_gate(
    metadata,
    required_columns=(),
    required_units=None,
    independent_axis=None,
    expected_source=None,
    expected_dataset_id=None,
    expected_solution_tag=None,
    expected_solution_artifact_id=None,
    expected_solution_digest=None,
    expected_sweep_axis_id=None,
    expected_sweep_axis_digest=None,
    expected_sweep_axis_row_count=None,
    expected_parameter_set_artifact_id=None,
    expected_parameter_set_digest=None,
    expected_parameter_set_path=None,
    expected_objective_observable_id=None,
    expected_objective_observable_family=None,
    expected_solver_configuration_artifact_id=None,
    expected_solver_configuration_digest=None,
    expected_solver_sequence_tag=None,
    expected_linear_solver=None,
    expected_relative_tolerance=None,
    expected_study_tag=None,
    expected_study_step_tag=None,
    expected_table_id=None,
    expected_selection_tags=None,
    expected_entity_dimensions=None,
    expected_expressions=None,
    expected_operator_tags=None,
    expected_result_table_schema_id=None,
    expected_result_output_artifact_id=None,
    expected_result_output_digest=None,
    expected_result_observable_id=None,
    expected_result_observable_family=None,
    expected_result_row_convention=None,
    expected_result_normalization_basis=None,
    expected_result_evaluation_method=None,
    expected_physics_convention_schema_id=None,
    expected_result_postprocess_row_convention_schema_id=None,
    expected_result_component_basis_schema_id=None,
    expected_result_artifact_id=None,
    expected_comsol_version=None,
    require_solver_configuration=False,
    require_parameter_set_artifact=False,
    require_result_provenance=False,
    require_result_table_schema=False,
    require_result_output_artifact=False,
    require_physics_convention_schema=False,
    require_result_postprocess_row_convention_schema=False,
    require_result_component_basis_schema=False,
    min_rows=1,
):
    """Check solver result-table metadata before interpreting numeric rows.

    This is a solver-independent companion to COMSOL ``mphglobal``/``mphtable``,
    MATLAB tables, CSV exports, and measurement logs.  It verifies that columns,
    units, an independent sweep axis, row count, optional source label, and
    optional dataset/solution/sweep/parameter/objective/solver-configuration/study/table identity,
    optional selection/entity scope, optional expression/operator identity, and
    optional table-schema identity are explicit before a table is used as
    validation evidence.  When a table is
    promoted to a replayable result package, the run artifact id, solver version,
    timestamp, and coarse timing breakdown should travel with it too.
    """

    if not isinstance(metadata, dict):
        raise ValueError("metadata must be a mapping")

    def pick(*keys, default=None):
        for key in keys:
            if key in metadata and metadata[key] is not None:
                return metadata[key]
        return default

    def normalize_columns(value):
        if value is None:
            return []
        if isinstance(value, str):
            return [part.strip() for part in value.replace(";", ",").split(",") if part.strip()]
        out = []
        for item in value:
            if isinstance(item, dict):
                raw = item.get("name", item.get("label", item.get("id", "")))
            else:
                raw = item
            text = str(raw).strip()
            if text:
                out.append(text)
        return out

    def normalize_string_list(value):
        if value is None:
            return []
        if isinstance(value, str):
            return [part.strip() for part in value.replace(";", ",").split(",") if part.strip()]
        out = []
        try:
            iterator = iter(value)
        except TypeError:
            iterator = iter((value,))
        for item in iterator:
            if isinstance(item, dict):
                raw = item.get("tag", item.get("name", item.get("label", item.get("id", ""))))
            else:
                raw = item
            text = str(raw).strip()
            if text:
                out.append(text)
        return out

    def collect_string_values(*keys):
        values = []
        for key in keys:
            if key not in metadata or metadata[key] is None:
                continue
            raw_values = normalize_string_list(metadata[key])
            values.extend(raw_values)
        return values

    def unique_strings(values):
        return sorted({str(value).strip() for value in values if str(value).strip()})

    def unique_ints(values):
        out = set()
        for value in values:
            text = str(value).strip()
            if not text:
                continue
            out.add(int(text))
        return sorted(out)

    def unique_floats(values):
        out = set()
        for value in values:
            text = str(value).strip()
            if not text:
                continue
            out.add(float(text))
        return sorted(out)

    def normalize_label(value):
        return str(value).strip().lower().replace("-", "_").replace(" ", "_")

    def parse_datetime(value):
        text = str(value or "").strip()
        if not text:
            return None
        try:
            return datetime.fromisoformat(text.replace("Z", "+00:00"))
        except ValueError:
            return None

    def parse_time(value):
        return parse_datetime(value) is not None

    def pick_from(mapping, *keys, default=None):
        if not isinstance(mapping, dict):
            return default
        for key in keys:
            if key in mapping and mapping[key] is not None:
                return mapping[key]
        return default

    def collect_string_values_from(mapping, *keys):
        values = []
        if not isinstance(mapping, dict):
            return values
        for key in keys:
            if key not in mapping or mapping[key] is None:
                continue
            values.extend(normalize_string_list(mapping[key]))
        return values

    def pick_mapping(*values):
        for value in values:
            if value is None:
                continue
            if isinstance(value, dict):
                return value
            raise ValueError("timing breakdown must be a mapping when present")
        return {}

    def normalize_entity_dimensions(value):
        aliases = {
            "3": "domain",
            "domain": "domain",
            "domains": "domain",
            "vol": "domain",
            "volume": "domain",
            "volumes": "domain",
            "2": "boundary",
            "boundary": "boundary",
            "boundaries": "boundary",
            "bnd": "boundary",
            "face": "boundary",
            "faces": "boundary",
            "surface": "boundary",
            "surfaces": "boundary",
            "1": "edge",
            "edge": "edge",
            "edges": "edge",
            "0": "point",
            "point": "point",
            "points": "point",
            "vertex": "point",
            "vertices": "point",
        }
        return [aliases.get(item.strip().lower(), item.strip().lower()) for item in normalize_string_list(value)]

    columns = normalize_columns(pick("columns", "column_names", default=[]))
    required = [str(item).strip() for item in required_columns if str(item).strip()]
    units = pick("units", "column_units", default={})
    if units is None:
        units = {}
    if not isinstance(units, dict):
        raise ValueError("units must be a mapping from column name to unit")
    expected_units = {} if required_units is None else {str(k): str(v) for k, v in dict(required_units).items()}
    axis = pick("independent_axis", "sweep_axis", "x_axis", default=independent_axis)
    axis_text = None if axis is None else str(axis).strip()
    source = pick("source", "solver", "tool", default=None)
    source_text = None if source is None else str(source).strip()
    dataset_id = pick("dataset_id", "dataset", "dataset_tag", "dataset_name", default=None)
    dataset_text = None if dataset_id is None else str(dataset_id).strip()
    solution_tag = pick("solution_tag", "solution_id", "solution", "sol_tag", default=None)
    solution_text = None if solution_tag is None else str(solution_tag).strip()
    execution = pick("execution", "run", "result_provenance", default={}) or {}
    if not isinstance(execution, dict):
        raise ValueError("execution/result_provenance must be a mapping when present")
    solution_artifact_ids = unique_strings(
        collect_string_values(
            "solution_artifact_id",
            "solutionArtifactId",
            "solution_data_artifact_id",
            "solutionDataArtifactId",
            "solver_solution_artifact_id",
            "solverSolutionArtifactId",
            "solution_output_artifact_id",
            "solutionOutputArtifactId",
        )
        + collect_string_values_from(
            execution,
            "solution_artifact_id",
            "solutionArtifactId",
            "solution_data_artifact_id",
            "solutionDataArtifactId",
            "solver_solution_artifact_id",
            "solverSolutionArtifactId",
            "solution_output_artifact_id",
            "solutionOutputArtifactId",
        )
    )
    solution_digests = unique_strings(
        collect_string_values(
            "solution_digest",
            "solutionDigest",
            "solution_sha256",
            "solutionSha256",
            "solution_data_digest",
            "solutionDataDigest",
            "solution_data_sha256",
            "solutionDataSha256",
            "solver_solution_digest",
            "solverSolutionDigest",
        )
        + collect_string_values_from(
            execution,
            "solution_digest",
            "solutionDigest",
            "solution_sha256",
            "solutionSha256",
            "solution_data_digest",
            "solutionDataDigest",
            "solution_data_sha256",
            "solutionDataSha256",
            "solver_solution_digest",
            "solverSolutionDigest",
        )
    )
    solution_artifact_id = solution_artifact_ids[0] if solution_artifact_ids else None
    solution_digest = solution_digests[0] if solution_digests else None
    sweep_axis_ids = unique_strings(
        collect_string_values(
            "sweep_axis_id",
            "sweepAxisId",
            "sweep_grid_id",
            "sweepGridId",
            "parameter_axis_id",
            "parameterAxisId",
            "parameter_grid_id",
            "parameterGridId",
            "frequency_grid_id",
            "frequencyGridId",
            "independent_axis_id",
            "independentAxisId",
        )
        + collect_string_values_from(
            execution,
            "sweep_axis_id",
            "sweepAxisId",
            "sweep_grid_id",
            "sweepGridId",
            "parameter_axis_id",
            "parameterAxisId",
            "parameter_grid_id",
            "parameterGridId",
            "frequency_grid_id",
            "frequencyGridId",
            "independent_axis_id",
            "independentAxisId",
        )
    )
    sweep_axis_digests = unique_strings(
        collect_string_values(
            "sweep_axis_digest",
            "sweepAxisDigest",
            "sweep_axis_sha256",
            "sweepAxisSha256",
            "sweep_grid_digest",
            "sweepGridDigest",
            "parameter_axis_digest",
            "parameterAxisDigest",
            "parameter_grid_digest",
            "parameterGridDigest",
            "frequency_grid_digest",
            "frequencyGridDigest",
            "independent_axis_digest",
            "independentAxisDigest",
        )
        + collect_string_values_from(
            execution,
            "sweep_axis_digest",
            "sweepAxisDigest",
            "sweep_axis_sha256",
            "sweepAxisSha256",
            "sweep_grid_digest",
            "sweepGridDigest",
            "parameter_axis_digest",
            "parameterAxisDigest",
            "parameter_grid_digest",
            "parameterGridDigest",
            "frequency_grid_digest",
            "frequencyGridDigest",
            "independent_axis_digest",
            "independentAxisDigest",
        )
    )
    sweep_axis_row_counts = unique_ints(
        collect_string_values(
            "sweep_axis_row_count",
            "sweepAxisRowCount",
            "sweep_axis_count",
            "sweepAxisCount",
            "sweep_grid_row_count",
            "sweepGridRowCount",
            "parameter_axis_row_count",
            "parameterAxisRowCount",
            "parameter_grid_row_count",
            "parameterGridRowCount",
            "frequency_grid_row_count",
            "frequencyGridRowCount",
            "independent_axis_row_count",
            "independentAxisRowCount",
        )
        + collect_string_values_from(
            execution,
            "sweep_axis_row_count",
            "sweepAxisRowCount",
            "sweep_axis_count",
            "sweepAxisCount",
            "sweep_grid_row_count",
            "sweepGridRowCount",
            "parameter_axis_row_count",
            "parameterAxisRowCount",
            "parameter_grid_row_count",
            "parameterGridRowCount",
            "frequency_grid_row_count",
            "frequencyGridRowCount",
            "independent_axis_row_count",
            "independentAxisRowCount",
        )
    )
    sweep_axis_id = sweep_axis_ids[0] if sweep_axis_ids else None
    sweep_axis_digest = sweep_axis_digests[0] if sweep_axis_digests else None
    sweep_axis_row_count = sweep_axis_row_counts[0] if sweep_axis_row_counts else None
    optimization = pick(
        "optimization",
        "optimizer",
        "parameter_study",
        "parameterStudy",
        "parameter_sweep",
        "parameterSweep",
        default={},
    ) or {}
    if not isinstance(optimization, dict):
        optimization = {}
    objective_metadata = pick(
        "objective",
        "objective_metadata",
        "objectiveMetadata",
        "optimization_objective",
        "optimizationObjective",
        default={},
    ) or {}
    if not isinstance(objective_metadata, dict):
        objective_metadata = {}
    parameter_set_artifact_ids = unique_strings(
        collect_string_values(
            "parameter_set_artifact_id",
            "parameterSetArtifactId",
            "parameter_artifact_id",
            "parameterArtifactId",
            "parameter_table_artifact_id",
            "parameterTableArtifactId",
            "optimization_parameter_set_artifact_id",
            "optimizationParameterSetArtifactId",
        )
        + collect_string_values_from(
            execution,
            "parameter_set_artifact_id",
            "parameterSetArtifactId",
            "parameter_artifact_id",
            "parameterArtifactId",
            "parameter_table_artifact_id",
            "parameterTableArtifactId",
            "optimization_parameter_set_artifact_id",
            "optimizationParameterSetArtifactId",
        )
        + collect_string_values_from(
            optimization,
            "parameter_set_artifact_id",
            "parameterSetArtifactId",
            "parameter_artifact_id",
            "parameterArtifactId",
            "parameter_table_artifact_id",
            "parameterTableArtifactId",
            "optimization_parameter_set_artifact_id",
            "optimizationParameterSetArtifactId",
        )
    )
    parameter_set_digests = unique_strings(
        collect_string_values(
            "parameter_set_digest",
            "parameterSetDigest",
            "parameter_set_sha256",
            "parameterSetSha256",
            "parameter_digest",
            "parameterDigest",
            "parameter_table_digest",
            "parameterTableDigest",
            "optimization_parameter_set_digest",
            "optimizationParameterSetDigest",
        )
        + collect_string_values_from(
            execution,
            "parameter_set_digest",
            "parameterSetDigest",
            "parameter_set_sha256",
            "parameterSetSha256",
            "parameter_digest",
            "parameterDigest",
            "parameter_table_digest",
            "parameterTableDigest",
            "optimization_parameter_set_digest",
            "optimizationParameterSetDigest",
        )
        + collect_string_values_from(
            optimization,
            "parameter_set_digest",
            "parameterSetDigest",
            "parameter_set_sha256",
            "parameterSetSha256",
            "parameter_digest",
            "parameterDigest",
            "parameter_table_digest",
            "parameterTableDigest",
            "optimization_parameter_set_digest",
            "optimizationParameterSetDigest",
        )
    )
    parameter_set_paths = unique_strings(
        collect_string_values(
            "parameter_set_path",
            "parameterSetPath",
            "parameter_path",
            "parameterPath",
            "parameter_table_path",
            "parameterTablePath",
            "optimization_parameter_set_path",
            "optimizationParameterSetPath",
        )
        + collect_string_values_from(
            execution,
            "parameter_set_path",
            "parameterSetPath",
            "parameter_path",
            "parameterPath",
            "parameter_table_path",
            "parameterTablePath",
            "optimization_parameter_set_path",
            "optimizationParameterSetPath",
        )
        + collect_string_values_from(
            optimization,
            "parameter_set_path",
            "parameterSetPath",
            "parameter_path",
            "parameterPath",
            "parameter_table_path",
            "parameterTablePath",
            "optimization_parameter_set_path",
            "optimizationParameterSetPath",
        )
    )
    objective_observable_ids = unique_strings(
        collect_string_values(
            "objective_observable_id",
            "objectiveObservableId",
            "objective_id",
            "objectiveId",
            "optimization_objective_id",
            "optimizationObjectiveId",
            "target_observable_id",
            "targetObservableId",
        )
        + collect_string_values_from(
            execution,
            "objective_observable_id",
            "objectiveObservableId",
            "objective_id",
            "objectiveId",
            "optimization_objective_id",
            "optimizationObjectiveId",
            "target_observable_id",
            "targetObservableId",
        )
        + collect_string_values_from(
            optimization,
            "objective_observable_id",
            "objectiveObservableId",
            "objective_id",
            "objectiveId",
            "optimization_objective_id",
            "optimizationObjectiveId",
            "target_observable_id",
            "targetObservableId",
        )
        + collect_string_values_from(
            objective_metadata,
            "objective_observable_id",
            "objectiveObservableId",
            "objective_id",
            "objectiveId",
            "optimization_objective_id",
            "optimizationObjectiveId",
            "target_observable_id",
            "targetObservableId",
        )
    )
    objective_observable_families = sorted({
        normalize_label(item)
        for item in (
            collect_string_values(
                "objective_observable_family",
                "objectiveObservableFamily",
                "objective_family",
                "objectiveFamily",
                "optimization_objective_family",
                "optimizationObjectiveFamily",
                "target_observable_family",
                "targetObservableFamily",
            )
            + collect_string_values_from(
                execution,
                "objective_observable_family",
                "objectiveObservableFamily",
                "objective_family",
                "objectiveFamily",
                "optimization_objective_family",
                "optimizationObjectiveFamily",
                "target_observable_family",
                "targetObservableFamily",
            )
            + collect_string_values_from(
                optimization,
                "objective_observable_family",
                "objectiveObservableFamily",
                "objective_family",
                "objectiveFamily",
                "optimization_objective_family",
                "optimizationObjectiveFamily",
                "target_observable_family",
                "targetObservableFamily",
            )
            + collect_string_values_from(
                objective_metadata,
                "objective_observable_family",
                "objectiveObservableFamily",
                "objective_family",
                "objectiveFamily",
                "optimization_objective_family",
                "optimizationObjectiveFamily",
                "target_observable_family",
                "targetObservableFamily",
            )
        )
        if str(item).strip()
    })
    parameter_set_artifact_id = (
        parameter_set_artifact_ids[0] if parameter_set_artifact_ids else None
    )
    parameter_set_digest = parameter_set_digests[0] if parameter_set_digests else None
    parameter_set_path = parameter_set_paths[0] if parameter_set_paths else None
    objective_observable_id = objective_observable_ids[0] if objective_observable_ids else None
    objective_observable_family = (
        objective_observable_families[0] if objective_observable_families else None
    )
    solver_configuration = pick("solver_configuration", "solver_config", "solver_settings", default={}) or {}
    if not isinstance(solver_configuration, dict):
        raise ValueError("solver_configuration/solver_config/solver_settings must be a mapping when present")
    execution_solver_configuration = pick_from(
        execution,
        "solver_configuration",
        "solverConfiguration",
        "solver_config",
        "solverConfig",
        "solver_settings",
        "solverSettings",
        default={},
    ) or {}
    if not isinstance(execution_solver_configuration, dict):
        raise ValueError("execution solver configuration must be a mapping when present")
    solver_configuration_artifact_ids = unique_strings(
        collect_string_values(
            "solver_configuration_artifact_id",
            "solverConfigurationArtifactId",
            "solver_config_artifact_id",
            "solverConfigArtifactId",
            "solver_settings_artifact_id",
            "solverSettingsArtifactId",
        )
        + collect_string_values_from(
            execution,
            "solver_configuration_artifact_id",
            "solverConfigurationArtifactId",
            "solver_config_artifact_id",
            "solverConfigArtifactId",
            "solver_settings_artifact_id",
            "solverSettingsArtifactId",
        )
        + collect_string_values_from(
            solver_configuration,
            "artifact_id",
            "artifactId",
            "solver_configuration_artifact_id",
            "solverConfigurationArtifactId",
            "solver_config_artifact_id",
            "solverConfigArtifactId",
        )
        + collect_string_values_from(
            execution_solver_configuration,
            "artifact_id",
            "artifactId",
            "solver_configuration_artifact_id",
            "solverConfigurationArtifactId",
            "solver_config_artifact_id",
            "solverConfigArtifactId",
        )
    )
    solver_configuration_digests = unique_strings(
        collect_string_values(
            "solver_configuration_digest",
            "solverConfigurationDigest",
            "solver_configuration_sha256",
            "solverConfigurationSha256",
            "solver_config_digest",
            "solverConfigDigest",
            "solver_settings_digest",
            "solverSettingsDigest",
        )
        + collect_string_values_from(
            execution,
            "solver_configuration_digest",
            "solverConfigurationDigest",
            "solver_configuration_sha256",
            "solverConfigurationSha256",
            "solver_config_digest",
            "solverConfigDigest",
            "solver_settings_digest",
            "solverSettingsDigest",
        )
        + collect_string_values_from(
            solver_configuration,
            "digest",
            "sha256",
            "solver_configuration_digest",
            "solverConfigurationDigest",
            "solver_config_digest",
            "solverConfigDigest",
        )
        + collect_string_values_from(
            execution_solver_configuration,
            "digest",
            "sha256",
            "solver_configuration_digest",
            "solverConfigurationDigest",
            "solver_config_digest",
            "solverConfigDigest",
        )
    )
    solver_sequence_tags = unique_strings(
        collect_string_values(
            "solver_sequence_tag",
            "solverSequenceTag",
            "solver_sequence",
            "solverSequence",
            "solver_tag",
            "solverTag",
            "solution_sequence_tag",
            "solutionSequenceTag",
        )
        + collect_string_values_from(
            execution,
            "solver_sequence_tag",
            "solverSequenceTag",
            "solver_sequence",
            "solverSequence",
            "solver_tag",
            "solverTag",
            "solution_sequence_tag",
            "solutionSequenceTag",
        )
        + collect_string_values_from(
            solver_configuration,
            "solver_sequence_tag",
            "solverSequenceTag",
            "solver_sequence",
            "solverSequence",
            "solver_tag",
            "solverTag",
        )
        + collect_string_values_from(
            execution_solver_configuration,
            "solver_sequence_tag",
            "solverSequenceTag",
            "solver_sequence",
            "solverSequence",
            "solver_tag",
            "solverTag",
        )
    )
    linear_solvers = unique_strings(
        collect_string_values(
            "linear_solver",
            "linearSolver",
            "linear_solver_name",
            "linearSolverName",
            "solver_linear_method",
            "solverLinearMethod",
        )
        + collect_string_values_from(
            execution,
            "linear_solver",
            "linearSolver",
            "linear_solver_name",
            "linearSolverName",
            "solver_linear_method",
            "solverLinearMethod",
        )
        + collect_string_values_from(
            solver_configuration,
            "linear_solver",
            "linearSolver",
            "linear_solver_name",
            "linearSolverName",
            "solver_linear_method",
            "solverLinearMethod",
        )
        + collect_string_values_from(
            execution_solver_configuration,
            "linear_solver",
            "linearSolver",
            "linear_solver_name",
            "linearSolverName",
            "solver_linear_method",
            "solverLinearMethod",
        )
    )
    relative_tolerances = unique_floats(
        collect_string_values(
            "relative_tolerance",
            "relativeTolerance",
            "solver_relative_tolerance",
            "solverRelativeTolerance",
            "rtol",
        )
        + collect_string_values_from(
            execution,
            "relative_tolerance",
            "relativeTolerance",
            "solver_relative_tolerance",
            "solverRelativeTolerance",
            "rtol",
        )
        + collect_string_values_from(
            solver_configuration,
            "relative_tolerance",
            "relativeTolerance",
            "solver_relative_tolerance",
            "solverRelativeTolerance",
            "rtol",
        )
        + collect_string_values_from(
            execution_solver_configuration,
            "relative_tolerance",
            "relativeTolerance",
            "solver_relative_tolerance",
            "solverRelativeTolerance",
            "rtol",
        )
    )
    solver_configuration_artifact_id = (
        solver_configuration_artifact_ids[0] if solver_configuration_artifact_ids else None
    )
    solver_configuration_digest = solver_configuration_digests[0] if solver_configuration_digests else None
    solver_sequence_tag = solver_sequence_tags[0] if solver_sequence_tags else None
    linear_solver = linear_solvers[0] if linear_solvers else None
    relative_tolerance = relative_tolerances[0] if relative_tolerances else None
    study_tag = pick("study_tag", "study_id", "study", "study_name", default=None)
    study_text = None if study_tag is None else str(study_tag).strip()
    study_step_tag = pick(
        "study_step_tag",
        "study_feature_tag",
        "study_step",
        "study_feature",
        "study_step_id",
        default=None,
    )
    study_step_text = None if study_step_tag is None else str(study_step_tag).strip()
    table_id = pick("table_id", "table_tag", "result_table_id", "export_id", default=None)
    table_text = None if table_id is None else str(table_id).strip()
    selection_tags = normalize_string_list(
        pick("selection_tags", "selection_tag", "selection_ids", "selection_id", "selection", default=[])
    )
    entity_dimensions = normalize_entity_dimensions(
        pick("entity_dimensions", "entity_dimension", "entity_dims", "entity_dim", "edim", default=[])
    )
    expressions = normalize_string_list(
        pick("expressions", "expression_names", "result_expressions", "exprs", "expr", default=[])
    )
    operator_tags = normalize_string_list(
        pick(
            "operator_tags",
            "postprocess_operator_tags",
            "integration_operator_tags",
            "probe_tags",
            "operator_tag",
            "probe_tag",
            default=[],
        )
    )
    result_table_schema_ids = unique_strings(
        collect_string_values(
            "result_table_schema_id",
            "resultTableSchemaId",
            "table_schema_id",
            "tableSchemaId",
            "output_table_schema_id",
            "outputTableSchemaId",
            "derived_value_table_schema_id",
            "derivedValueTableSchemaId",
        )
        + collect_string_values_from(
            execution,
            "result_table_schema_id",
            "resultTableSchemaId",
            "table_schema_id",
            "tableSchemaId",
            "output_table_schema_id",
            "outputTableSchemaId",
            "derived_value_table_schema_id",
            "derivedValueTableSchemaId",
        )
    )
    result_table_schema_id = (
        result_table_schema_ids[0] if result_table_schema_ids else None
    )
    physics_convention_schema_ids = unique_strings(
        collect_string_values(
            "physics_convention_schema_id",
            "physicsConventionSchemaId",
            "result_physics_convention_schema_id",
            "resultPhysicsConventionSchemaId",
            "derived_value_convention_schema_id",
            "derivedValueConventionSchemaId",
            "comsol_derived_value_convention_schema_id",
            "comsolDerivedValueConventionSchemaId",
            "field_convention_schema_id",
            "fieldConventionSchemaId",
        )
        + collect_string_values_from(
            execution,
            "physics_convention_schema_id",
            "physicsConventionSchemaId",
            "result_physics_convention_schema_id",
            "resultPhysicsConventionSchemaId",
            "derived_value_convention_schema_id",
            "derivedValueConventionSchemaId",
            "comsol_derived_value_convention_schema_id",
            "comsolDerivedValueConventionSchemaId",
            "field_convention_schema_id",
            "fieldConventionSchemaId",
        )
    )
    physics_convention_schema_id = (
        physics_convention_schema_ids[0] if physics_convention_schema_ids else None
    )
    result_postprocess_row_convention_schema_ids = unique_strings(
        collect_string_values(
            "result_postprocess_row_convention_schema_id",
            "resultPostprocessRowConventionSchemaId",
            "postprocess_row_convention_schema_id",
            "postprocessRowConventionSchemaId",
            "derived_value_postprocess_row_convention_schema_id",
            "derivedValuePostprocessRowConventionSchemaId",
            "comsol_postprocess_row_convention_schema_id",
            "comsolPostprocessRowConventionSchemaId",
        )
        + collect_string_values_from(
            execution,
            "result_postprocess_row_convention_schema_id",
            "resultPostprocessRowConventionSchemaId",
            "postprocess_row_convention_schema_id",
            "postprocessRowConventionSchemaId",
            "derived_value_postprocess_row_convention_schema_id",
            "derivedValuePostprocessRowConventionSchemaId",
            "comsol_postprocess_row_convention_schema_id",
            "comsolPostprocessRowConventionSchemaId",
        )
    )
    result_postprocess_row_convention_schema_id = (
        result_postprocess_row_convention_schema_ids[0]
        if result_postprocess_row_convention_schema_ids
        else None
    )
    result_component_basis_schema_ids = unique_strings(
        collect_string_values(
            "result_component_basis_schema_id",
            "resultComponentBasisSchemaId",
            "component_basis_schema_id",
            "componentBasisSchemaId",
            "derived_value_component_basis_schema_id",
            "derivedValueComponentBasisSchemaId",
            "table_component_basis_schema_id",
            "tableComponentBasisSchemaId",
            "field_component_basis_schema_id",
            "fieldComponentBasisSchemaId",
        )
        + collect_string_values_from(
            execution,
            "result_component_basis_schema_id",
            "resultComponentBasisSchemaId",
            "component_basis_schema_id",
            "componentBasisSchemaId",
            "derived_value_component_basis_schema_id",
            "derivedValueComponentBasisSchemaId",
            "table_component_basis_schema_id",
            "tableComponentBasisSchemaId",
            "field_component_basis_schema_id",
            "fieldComponentBasisSchemaId",
        )
    )
    result_component_basis_schema_id = (
        result_component_basis_schema_ids[0]
        if result_component_basis_schema_ids
        else None
    )
    result_output_artifact_ids = unique_strings(
        collect_string_values(
            "result_output_artifact_id",
            "output_artifact_id",
            "table_output_artifact_id",
            "postprocess_output_artifact_id",
            "export_output_artifact_id",
        )
    )
    result_output_digests = unique_strings(
        collect_string_values(
            "result_output_digest",
            "output_digest",
            "result_output_sha256",
            "output_sha256",
            "table_output_digest",
            "postprocess_output_digest",
            "export_output_digest",
        )
    )
    result_output_paths = unique_strings(
        collect_string_values(
            "result_output_path",
            "output_path",
            "table_output_path",
            "postprocess_output_path",
            "export_output_path",
        )
    )
    result_output_artifact_id = result_output_artifact_ids[0] if result_output_artifact_ids else None
    result_output_digest = result_output_digests[0] if result_output_digests else None
    result_output_path = result_output_paths[0] if result_output_paths else None
    result_observable_ids = unique_strings(
        collect_string_values(
            "result_observable_id",
            "table_observable_id",
            "postprocess_observable_id",
            "observable_id",
        )
    )
    result_observable_families = sorted({
        normalize_label(item)
        for item in collect_string_values(
            "result_observable_family",
            "table_observable_family",
            "postprocess_observable_family",
            "observable_family",
        )
        if item.strip()
    })
    result_observable_id = result_observable_ids[0] if result_observable_ids else None
    result_observable_family = result_observable_families[0] if result_observable_families else None
    result_row_conventions = sorted({
        normalize_label(item)
        for item in collect_string_values(
            "result_row_convention",
            "table_row_convention",
            "result_table_row_convention",
            "row_convention",
            "row_identity_convention",
            "row_basis",
            "row_semantics",
        )
        if item.strip()
    })
    result_normalization_bases = sorted({
        normalize_label(item)
        for item in collect_string_values(
            "result_normalization_basis",
            "result_table_normalization_basis",
            "normalization_basis",
            "quantity_normalization",
            "power_normalization",
            "power_balance_basis",
            "normalization",
        )
        if item.strip()
    })
    result_evaluation_methods = sorted({
        normalize_label(item)
        for item in collect_string_values(
            "result_evaluation_method",
            "table_evaluation_method",
            "postprocess_method",
            "evaluation_method",
            "readout_method",
            "result_api",
            "table_api",
        )
        if item.strip()
    })
    result_row_convention = result_row_conventions[0] if result_row_conventions else None
    result_normalization_basis = result_normalization_bases[0] if result_normalization_bases else None
    result_evaluation_method = result_evaluation_methods[0] if result_evaluation_methods else None
    result_artifact_ids = unique_strings(
        collect_string_values(
            "result_artifact_id",
            "resultArtifactId",
            "execution_result_artifact_id",
            "notebook_result_artifact_id",
            "run_result_artifact_id",
        )
        + collect_string_values_from(
            execution,
            "result_artifact_id",
            "resultArtifactId",
            "execution_result_artifact_id",
            "notebook_result_artifact_id",
            "run_result_artifact_id",
        )
    )
    result_artifact_id = result_artifact_ids[0] if result_artifact_ids else None
    run_started = (
        pick("run_started_at", "runStartedAt", "run_date", "runDate", "run_date_utc", "created_at", default=None)
        or pick_from(
            execution,
            "run_started_at",
            "runStartedAt",
            "run_date",
            "runDate",
            "run_date_utc",
            "created_at",
            default=None,
        )
    )
    run_started_text = None if run_started is None else str(run_started).strip()
    comsol_versions = unique_strings(
        collect_string_values(
            "comsol_version",
            "comsolVersion",
            "solver_version",
            "solverVersion",
        )
        + collect_string_values_from(
            execution,
            "comsol_version",
            "comsolVersion",
            "solver_version",
            "solverVersion",
        )
    )
    comsol_version = comsol_versions[0] if comsol_versions else None
    timing_breakdown = pick_mapping(
        pick("timing_breakdown_s", "timing_breakdown", "timingBreakdown", "timings", default=None),
        pick_from(
            execution,
            "timing_breakdown_s",
            "timing_breakdown",
            "timingBreakdown",
            "timings",
            default=None,
        ),
    )
    timing_breakdown_seconds = {}
    for name, value in timing_breakdown.items():
        try:
            timing_breakdown_seconds[str(name)] = float(value)
        except (TypeError, ValueError):
            timing_breakdown_seconds[str(name)] = math.nan
    timing_breakdown_names = sorted(timing_breakdown_seconds)
    row_count_raw = pick("row_count", "n_rows", default=None)
    row_count = None if row_count_raw is None else int(row_count_raw)
    min_row_count = int(min_rows)
    if min_row_count < 0:
        raise ValueError("min_rows must be non-negative")

    required_present = all(item in columns for item in required)
    units_for_required = all(item in units and str(units[item]).strip() for item in required)
    expected_units_ok = all(str(units.get(column, "")).strip() == unit for column, unit in expected_units.items())
    row_count_ok = row_count is not None and row_count >= min_row_count
    source_ok = True
    if expected_source is not None:
        source_ok = source_text == str(expected_source)
    expected_dataset = None if expected_dataset_id is None else str(expected_dataset_id).strip()
    expected_solution = None if expected_solution_tag is None else str(expected_solution_tag).strip()
    expected_solution_artifact = (
        None
        if expected_solution_artifact_id is None
        else str(expected_solution_artifact_id).strip()
    )
    expected_solution_data_digest = (
        None if expected_solution_digest is None else str(expected_solution_digest).strip()
    )
    expected_sweep_axis = None if expected_sweep_axis_id is None else str(expected_sweep_axis_id).strip()
    expected_sweep_axis_data_digest = (
        None if expected_sweep_axis_digest is None else str(expected_sweep_axis_digest).strip()
    )
    expected_sweep_axis_count = (
        None if expected_sweep_axis_row_count is None else int(expected_sweep_axis_row_count)
    )
    expected_parameter_set_artifact = (
        None
        if expected_parameter_set_artifact_id is None
        else str(expected_parameter_set_artifact_id).strip()
    )
    expected_parameter_set_data_digest = (
        None
        if expected_parameter_set_digest is None
        else str(expected_parameter_set_digest).strip()
    )
    expected_parameter_set_file_path = (
        None
        if expected_parameter_set_path is None
        else str(expected_parameter_set_path).strip()
    )
    expected_objective_observable = (
        None
        if expected_objective_observable_id is None
        else str(expected_objective_observable_id).strip()
    )
    expected_objective_observable_group = (
        None
        if expected_objective_observable_family is None
        else normalize_label(expected_objective_observable_family)
    )
    expected_solver_configuration_artifact = (
        None
        if expected_solver_configuration_artifact_id is None
        else str(expected_solver_configuration_artifact_id).strip()
    )
    expected_solver_configuration_data_digest = (
        None
        if expected_solver_configuration_digest is None
        else str(expected_solver_configuration_digest).strip()
    )
    expected_solver_sequence = (
        None if expected_solver_sequence_tag is None else str(expected_solver_sequence_tag).strip()
    )
    expected_linear_solver_name = (
        None if expected_linear_solver is None else str(expected_linear_solver).strip()
    )
    expected_solver_relative_tolerance = (
        None if expected_relative_tolerance is None else float(expected_relative_tolerance)
    )
    expected_study = None if expected_study_tag is None else str(expected_study_tag).strip()
    expected_study_step = (
        None if expected_study_step_tag is None else str(expected_study_step_tag).strip()
    )
    expected_table = None if expected_table_id is None else str(expected_table_id).strip()
    expected_selections = normalize_string_list(expected_selection_tags)
    expected_edims = normalize_entity_dimensions(expected_entity_dimensions)
    expected_exprs = normalize_string_list(expected_expressions)
    expected_ops = normalize_string_list(expected_operator_tags)
    expected_result_table_schema = (
        None
        if expected_result_table_schema_id is None
        else str(expected_result_table_schema_id).strip()
    )
    expected_result_output = (
        None
        if expected_result_output_artifact_id is None
        else str(expected_result_output_artifact_id).strip()
    )
    expected_result_digest = (
        None if expected_result_output_digest is None else str(expected_result_output_digest).strip()
    )
    expected_result_observable = (
        None if expected_result_observable_id is None else str(expected_result_observable_id).strip()
    )
    expected_result_observable_group = (
        None
        if expected_result_observable_family is None
        else normalize_label(expected_result_observable_family)
    )
    expected_row_convention = (
        None
        if expected_result_row_convention is None
        else normalize_label(expected_result_row_convention)
    )
    expected_normalization_basis = (
        None
        if expected_result_normalization_basis is None
        else normalize_label(expected_result_normalization_basis)
    )
    expected_evaluation_method = (
        None
        if expected_result_evaluation_method is None
        else normalize_label(expected_result_evaluation_method)
    )
    expected_physics_convention_schema = (
        None
        if expected_physics_convention_schema_id is None
        else str(expected_physics_convention_schema_id).strip()
    )
    expected_result_postprocess_row_convention_schema = (
        None
        if expected_result_postprocess_row_convention_schema_id is None
        else str(expected_result_postprocess_row_convention_schema_id).strip()
    )
    expected_result_component_basis_schema = (
        None
        if expected_result_component_basis_schema_id is None
        else str(expected_result_component_basis_schema_id).strip()
    )
    expected_result_artifact = (
        None
        if expected_result_artifact_id is None
        else str(expected_result_artifact_id).strip()
    )
    expected_solver_version = (
        None if expected_comsol_version is None else str(expected_comsol_version).strip()
    )
    result_table_schema_required = bool(
        require_result_table_schema or expected_result_table_schema is not None
    )
    physics_convention_schema_required = bool(
        require_physics_convention_schema or expected_physics_convention_schema is not None
    )
    result_postprocess_row_convention_schema_required = bool(
        require_result_postprocess_row_convention_schema
        or expected_result_postprocess_row_convention_schema is not None
    )
    result_component_basis_schema_required = bool(
        require_result_component_basis_schema
        or expected_result_component_basis_schema is not None
    )
    result_output_required = bool(require_result_output_artifact)
    result_provenance_required = bool(
        require_result_provenance
        or expected_result_artifact is not None
        or expected_solver_version is not None
    )
    solver_configuration_required = bool(
        require_solver_configuration
        or expected_solver_configuration_artifact is not None
        or expected_solver_configuration_data_digest is not None
        or expected_solver_sequence is not None
        or expected_solver_relative_tolerance is not None
    )
    parameter_set_required = bool(
        require_parameter_set_artifact
        or expected_parameter_set_artifact is not None
        or expected_parameter_set_data_digest is not None
        or expected_parameter_set_file_path is not None
    )
    selection_scope_consistent = (
        not selection_tags
        or not entity_dimensions
        or len(entity_dimensions) == 1
        or len(entity_dimensions) == len(selection_tags)
    )
    checks = {
        "columns_recorded": bool(columns),
        "columns_unique": len(columns) == len(set(columns)),
        "required_columns_present": required_present,
        "units_recorded_for_required_columns": units_for_required,
        "expected_units_match": expected_units_ok,
        "independent_axis_recorded": bool(axis_text),
        "independent_axis_is_column": axis_text in columns if axis_text else False,
        "row_count_recorded": row_count is not None,
        "row_count_at_least_minimum": row_count_ok,
        "source_matches_expected": source_ok,
        "dataset_id_recorded": expected_dataset is None or bool(dataset_text),
        "expected_dataset_id_matches": expected_dataset is None or dataset_text == expected_dataset,
        "solution_tag_recorded": expected_solution is None or bool(solution_text),
        "expected_solution_tag_matches": expected_solution is None or solution_text == expected_solution,
        "solution_artifact_id_consistent_when_present": len(solution_artifact_ids) <= 1,
        "solution_digest_consistent_when_present": len(solution_digests) <= 1,
        "solution_artifact_id_recorded_when_expected": (
            expected_solution_artifact is None or bool(solution_artifact_id)
        ),
        "expected_solution_artifact_id_matches": (
            expected_solution_artifact is None or solution_artifact_id == expected_solution_artifact
        ),
        "solution_digest_recorded_when_expected": (
            expected_solution_data_digest is None or bool(solution_digest)
        ),
        "expected_solution_digest_matches": (
            expected_solution_data_digest is None or solution_digest == expected_solution_data_digest
        ),
        "sweep_axis_id_consistent_when_present": len(sweep_axis_ids) <= 1,
        "sweep_axis_digest_consistent_when_present": len(sweep_axis_digests) <= 1,
        "sweep_axis_row_count_consistent_when_present": len(sweep_axis_row_counts) <= 1,
        "sweep_axis_id_recorded_when_expected": (
            expected_sweep_axis is None or bool(sweep_axis_id)
        ),
        "expected_sweep_axis_id_matches": (
            expected_sweep_axis is None or sweep_axis_id == expected_sweep_axis
        ),
        "sweep_axis_digest_recorded_when_expected": (
            expected_sweep_axis_data_digest is None or bool(sweep_axis_digest)
        ),
        "expected_sweep_axis_digest_matches": (
            expected_sweep_axis_data_digest is None or sweep_axis_digest == expected_sweep_axis_data_digest
        ),
        "sweep_axis_row_count_recorded_when_expected": (
            expected_sweep_axis_count is None or sweep_axis_row_count is not None
        ),
        "expected_sweep_axis_row_count_matches": (
            expected_sweep_axis_count is None or sweep_axis_row_count == expected_sweep_axis_count
        ),
        "sweep_axis_row_count_matches_table_rows_when_present": (
            sweep_axis_row_count is None or row_count is None or sweep_axis_row_count == row_count
        ),
        "parameter_set_artifact_id_consistent_when_present": (
            len(parameter_set_artifact_ids) <= 1
        ),
        "parameter_set_digest_consistent_when_present": len(parameter_set_digests) <= 1,
        "parameter_set_path_consistent_when_present": len(parameter_set_paths) <= 1,
        "objective_observable_id_consistent_when_present": (
            len(objective_observable_ids) <= 1
        ),
        "objective_observable_family_consistent_when_present": (
            len(objective_observable_families) <= 1
        ),
        "parameter_set_artifact_id_recorded_when_required": (
            not parameter_set_required or bool(parameter_set_artifact_id)
        ),
        "parameter_set_digest_recorded_when_required": (
            not parameter_set_required or bool(parameter_set_digest)
        ),
        "parameter_set_path_recorded_when_required": (
            not parameter_set_required or bool(parameter_set_path)
        ),
        "parameter_set_artifact_id_recorded_when_expected": (
            expected_parameter_set_artifact is None or bool(parameter_set_artifact_id)
        ),
        "expected_parameter_set_artifact_id_matches": (
            expected_parameter_set_artifact is None
            or parameter_set_artifact_id == expected_parameter_set_artifact
        ),
        "parameter_set_digest_recorded_when_expected": (
            expected_parameter_set_data_digest is None or bool(parameter_set_digest)
        ),
        "expected_parameter_set_digest_matches": (
            expected_parameter_set_data_digest is None
            or parameter_set_digest == expected_parameter_set_data_digest
        ),
        "parameter_set_path_recorded_when_expected": (
            expected_parameter_set_file_path is None or bool(parameter_set_path)
        ),
        "expected_parameter_set_path_matches": (
            expected_parameter_set_file_path is None
            or parameter_set_path == expected_parameter_set_file_path
        ),
        "objective_observable_id_recorded_when_expected": (
            expected_objective_observable is None or bool(objective_observable_id)
        ),
        "expected_objective_observable_id_matches": (
            expected_objective_observable is None
            or objective_observable_id == expected_objective_observable
        ),
        "objective_observable_family_recorded_when_expected": (
            expected_objective_observable_group is None or bool(objective_observable_family)
        ),
        "expected_objective_observable_family_matches": (
            expected_objective_observable_group is None
            or objective_observable_family == expected_objective_observable_group
        ),
        "solver_configuration_artifact_id_consistent_when_present": (
            len(solver_configuration_artifact_ids) <= 1
        ),
        "solver_configuration_digest_consistent_when_present": (
            len(solver_configuration_digests) <= 1
        ),
        "solver_sequence_tag_consistent_when_present": len(solver_sequence_tags) <= 1,
        "linear_solver_consistent_when_present": len(linear_solvers) <= 1,
        "relative_tolerance_consistent_when_present": len(relative_tolerances) <= 1,
        "solver_configuration_artifact_id_recorded_when_required": (
            not solver_configuration_required or bool(solver_configuration_artifact_id)
        ),
        "solver_configuration_digest_recorded_when_required": (
            not solver_configuration_required or bool(solver_configuration_digest)
        ),
        "solver_sequence_tag_recorded_when_required": (
            not solver_configuration_required or bool(solver_sequence_tag)
        ),
        "relative_tolerance_recorded_when_required": (
            not solver_configuration_required or relative_tolerance is not None
        ),
        "expected_solver_configuration_artifact_id_matches": (
            expected_solver_configuration_artifact is None
            or solver_configuration_artifact_id == expected_solver_configuration_artifact
        ),
        "expected_solver_configuration_digest_matches": (
            expected_solver_configuration_data_digest is None
            or solver_configuration_digest == expected_solver_configuration_data_digest
        ),
        "expected_solver_sequence_tag_matches": (
            expected_solver_sequence is None or solver_sequence_tag == expected_solver_sequence
        ),
        "expected_linear_solver_matches": (
            expected_linear_solver_name is None or linear_solvers == [expected_linear_solver_name]
        ),
        "expected_relative_tolerance_matches": (
            expected_solver_relative_tolerance is None
            or (
                relative_tolerance is not None
                and math.isclose(
                    relative_tolerance,
                    expected_solver_relative_tolerance,
                    rel_tol=1.0e-12,
                    abs_tol=1.0e-15,
                )
            )
        ),
        "relative_tolerance_finite_positive_when_present": all(
            math.isfinite(value) and value > 0.0 for value in relative_tolerances
        ),
        "study_tag_recorded": expected_study is None or bool(study_text),
        "expected_study_tag_matches": expected_study is None or study_text == expected_study,
        "study_step_tag_recorded": expected_study_step is None or bool(study_step_text),
        "expected_study_step_tag_matches": (
            expected_study_step is None or study_step_text == expected_study_step
        ),
        "table_id_recorded": expected_table is None or bool(table_text),
        "expected_table_id_matches": expected_table is None or table_text == expected_table,
        "selection_tags_recorded": not expected_selections or bool(selection_tags),
        "expected_selection_tags_match": not expected_selections
        or sorted(selection_tags) == sorted(expected_selections),
        "entity_dimensions_recorded": not expected_edims or bool(entity_dimensions),
        "expected_entity_dimensions_match": not expected_edims
        or sorted(entity_dimensions) == sorted(expected_edims),
        "selection_entity_scope_consistent": selection_scope_consistent,
        "expressions_recorded_when_expected": not expected_exprs or bool(expressions),
        "expected_expressions_present": not expected_exprs
        or set(expected_exprs).issubset(set(expressions)),
        "operator_tags_recorded_when_expected": not expected_ops or bool(operator_tags),
        "expected_operator_tags_match": not expected_ops or sorted(operator_tags) == sorted(expected_ops),
        "result_table_schema_id_consistent_when_present": len(result_table_schema_ids) <= 1,
        "result_table_schema_id_recorded_when_required": (
            not result_table_schema_required or bool(result_table_schema_id)
        ),
        "result_table_schema_id_recorded_when_expected": (
            expected_result_table_schema is None or bool(result_table_schema_id)
        ),
        "expected_result_table_schema_id_matches": (
            expected_result_table_schema is None
            or result_table_schema_id == expected_result_table_schema
        ),
        "physics_convention_schema_id_consistent_when_present": (
            len(physics_convention_schema_ids) <= 1
        ),
        "physics_convention_schema_id_recorded_when_required": (
            not physics_convention_schema_required or bool(physics_convention_schema_id)
        ),
        "physics_convention_schema_id_recorded_when_expected": (
            expected_physics_convention_schema is None
            or bool(physics_convention_schema_id)
        ),
        "expected_physics_convention_schema_id_matches": (
            expected_physics_convention_schema is None
            or physics_convention_schema_id == expected_physics_convention_schema
        ),
        "result_postprocess_row_convention_schema_id_consistent_when_present": (
            len(result_postprocess_row_convention_schema_ids) <= 1
        ),
        "result_postprocess_row_convention_schema_id_recorded_when_required": (
            not result_postprocess_row_convention_schema_required
            or bool(result_postprocess_row_convention_schema_id)
        ),
        "result_postprocess_row_convention_schema_id_recorded_when_expected": (
            expected_result_postprocess_row_convention_schema is None
            or bool(result_postprocess_row_convention_schema_id)
        ),
        "expected_result_postprocess_row_convention_schema_id_matches": (
            expected_result_postprocess_row_convention_schema is None
            or result_postprocess_row_convention_schema_id
            == expected_result_postprocess_row_convention_schema
        ),
        "result_component_basis_schema_id_consistent_when_present": (
            len(result_component_basis_schema_ids) <= 1
        ),
        "result_component_basis_schema_id_recorded_when_required": (
            not result_component_basis_schema_required
            or bool(result_component_basis_schema_id)
        ),
        "result_component_basis_schema_id_recorded_when_expected": (
            expected_result_component_basis_schema is None
            or bool(result_component_basis_schema_id)
        ),
        "expected_result_component_basis_schema_id_matches": (
            expected_result_component_basis_schema is None
            or result_component_basis_schema_id
            == expected_result_component_basis_schema
        ),
        "result_output_artifact_id_consistent_when_present": len(result_output_artifact_ids) <= 1,
        "result_output_digest_consistent_when_present": len(result_output_digests) <= 1,
        "result_output_path_consistent_when_present": len(result_output_paths) <= 1,
        "result_output_artifact_id_recorded_when_required": (
            not result_output_required or bool(result_output_artifact_id)
        ),
        "result_output_digest_recorded_when_required": (
            not result_output_required or bool(result_output_digest)
        ),
        "result_output_path_recorded_when_required": (
            not result_output_required or bool(result_output_path)
        ),
        "result_output_artifact_id_recorded_when_expected": (
            expected_result_output is None or bool(result_output_artifact_id)
        ),
        "expected_result_output_artifact_id_matches": (
            expected_result_output is None or result_output_artifact_id == expected_result_output
        ),
        "result_output_digest_recorded_when_expected": (
            expected_result_digest is None or bool(result_output_digest)
        ),
        "expected_result_output_digest_matches": (
            expected_result_digest is None or result_output_digest == expected_result_digest
        ),
        "result_observable_id_consistent_when_present": len(result_observable_ids) <= 1,
        "result_observable_family_consistent_when_present": len(result_observable_families) <= 1,
        "result_row_convention_consistent_when_present": len(result_row_conventions) <= 1,
        "result_normalization_basis_consistent_when_present": len(result_normalization_bases) <= 1,
        "result_evaluation_method_consistent_when_present": len(result_evaluation_methods) <= 1,
        "result_observable_id_recorded_when_expected": (
            expected_result_observable is None or bool(result_observable_id)
        ),
        "expected_result_observable_id_matches": (
            expected_result_observable is None or result_observable_ids == [expected_result_observable]
        ),
        "result_observable_family_recorded_when_expected": (
            expected_result_observable_group is None or bool(result_observable_family)
        ),
        "expected_result_observable_family_matches": (
            expected_result_observable_group is None
            or result_observable_families == [expected_result_observable_group]
        ),
        "result_row_convention_recorded_when_expected": (
            expected_row_convention is None or bool(result_row_convention)
        ),
        "expected_result_row_convention_matches": (
            expected_row_convention is None or result_row_conventions == [expected_row_convention]
        ),
        "result_normalization_basis_recorded_when_expected": (
            expected_normalization_basis is None or bool(result_normalization_basis)
        ),
        "expected_result_normalization_basis_matches": (
            expected_normalization_basis is None
            or result_normalization_bases == [expected_normalization_basis]
        ),
        "result_evaluation_method_recorded_when_expected": (
            expected_evaluation_method is None or bool(result_evaluation_method)
        ),
        "expected_result_evaluation_method_matches": (
            expected_evaluation_method is None
            or result_evaluation_methods == [expected_evaluation_method]
        ),
        "result_artifact_id_consistent_when_present": len(result_artifact_ids) <= 1,
        "comsol_version_consistent_when_present": len(comsol_versions) <= 1,
        "result_artifact_id_recorded_when_required": (
            not result_provenance_required or bool(result_artifact_id)
        ),
        "expected_result_artifact_id_matches": (
            expected_result_artifact is None or result_artifact_id == expected_result_artifact
        ),
        "run_started_at_recorded_when_required": (
            not result_provenance_required or bool(run_started_text)
        ),
        "run_started_at_parseable_when_present": (
            not run_started_text or parse_time(run_started_text)
        ),
        "comsol_version_recorded_when_required": (
            not result_provenance_required or bool(comsol_version)
        ),
        "expected_comsol_version_matches": (
            expected_solver_version is None or comsol_versions == [expected_solver_version]
        ),
        "timing_breakdown_recorded_when_required": (
            not result_provenance_required or bool(timing_breakdown_seconds)
        ),
        "timing_breakdown_has_at_least_four_items": (
            not result_provenance_required or len(timing_breakdown_seconds) >= 4
        ),
        "timing_breakdown_values_finite_nonnegative": all(
            math.isfinite(value) and value >= 0.0
            for value in timing_breakdown_seconds.values()
        ),
    }
    return {
        "policy": "solver_result_table_metadata_gate",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "columns": columns,
        "required_columns": required,
        "units": {str(k): str(v) for k, v in units.items()},
        "expected_units": expected_units,
        "independent_axis": axis_text,
        "row_count": row_count,
        "min_rows": min_row_count,
        "source": source_text,
        "expected_source": None if expected_source is None else str(expected_source),
        "dataset_id": dataset_text,
        "expected_dataset_id": expected_dataset,
        "solution_tag": solution_text,
        "expected_solution_tag": expected_solution,
        "solution_artifact_id": solution_artifact_id,
        "solution_digest": solution_digest,
        "solution_artifact_ids": solution_artifact_ids,
        "solution_digests": solution_digests,
        "expected_solution_artifact_id": expected_solution_artifact,
        "expected_solution_digest": expected_solution_data_digest,
        "sweep_axis_id": sweep_axis_id,
        "sweep_axis_digest": sweep_axis_digest,
        "sweep_axis_row_count": sweep_axis_row_count,
        "sweep_axis_ids": sweep_axis_ids,
        "sweep_axis_digests": sweep_axis_digests,
        "sweep_axis_row_counts": sweep_axis_row_counts,
        "expected_sweep_axis_id": expected_sweep_axis,
        "expected_sweep_axis_digest": expected_sweep_axis_data_digest,
        "expected_sweep_axis_row_count": expected_sweep_axis_count,
        "parameter_set_artifact_id": parameter_set_artifact_id,
        "parameter_set_digest": parameter_set_digest,
        "parameter_set_path": parameter_set_path,
        "parameter_set_artifact_ids": parameter_set_artifact_ids,
        "parameter_set_digests": parameter_set_digests,
        "parameter_set_paths": parameter_set_paths,
        "expected_parameter_set_artifact_id": expected_parameter_set_artifact,
        "expected_parameter_set_digest": expected_parameter_set_data_digest,
        "expected_parameter_set_path": expected_parameter_set_file_path,
        "objective_observable_id": objective_observable_id,
        "objective_observable_family": objective_observable_family,
        "objective_observable_ids": objective_observable_ids,
        "objective_observable_families": objective_observable_families,
        "expected_objective_observable_id": expected_objective_observable,
        "expected_objective_observable_family": expected_objective_observable_group,
        "solver_configuration_artifact_id": solver_configuration_artifact_id,
        "solver_configuration_digest": solver_configuration_digest,
        "solver_sequence_tag": solver_sequence_tag,
        "linear_solver": linear_solver,
        "relative_tolerance": relative_tolerance,
        "solver_configuration_artifact_ids": solver_configuration_artifact_ids,
        "solver_configuration_digests": solver_configuration_digests,
        "solver_sequence_tags": solver_sequence_tags,
        "linear_solvers": linear_solvers,
        "relative_tolerances": relative_tolerances,
        "expected_solver_configuration_artifact_id": expected_solver_configuration_artifact,
        "expected_solver_configuration_digest": expected_solver_configuration_data_digest,
        "expected_solver_sequence_tag": expected_solver_sequence,
        "expected_linear_solver": expected_linear_solver_name,
        "expected_relative_tolerance": expected_solver_relative_tolerance,
        "study_tag": study_text,
        "expected_study_tag": expected_study,
        "study_step_tag": study_step_text,
        "expected_study_step_tag": expected_study_step,
        "table_id": table_text,
        "expected_table_id": expected_table,
        "selection_tags": selection_tags,
        "expected_selection_tags": expected_selections,
        "entity_dimensions": entity_dimensions,
        "expected_entity_dimensions": expected_edims,
        "expressions": expressions,
        "expected_expressions": expected_exprs,
        "operator_tags": operator_tags,
        "expected_operator_tags": expected_ops,
        "result_table_schema_id": result_table_schema_id,
        "result_table_schema_ids": result_table_schema_ids,
        "expected_result_table_schema_id": expected_result_table_schema,
        "physics_convention_schema_id": physics_convention_schema_id,
        "physics_convention_schema_ids": physics_convention_schema_ids,
        "expected_physics_convention_schema_id": expected_physics_convention_schema,
        "result_postprocess_row_convention_schema_id": result_postprocess_row_convention_schema_id,
        "result_postprocess_row_convention_schema_ids": (
            result_postprocess_row_convention_schema_ids
        ),
        "expected_result_postprocess_row_convention_schema_id": (
            expected_result_postprocess_row_convention_schema
        ),
        "result_component_basis_schema_id": result_component_basis_schema_id,
        "result_component_basis_schema_ids": result_component_basis_schema_ids,
        "expected_result_component_basis_schema_id": (
            expected_result_component_basis_schema
        ),
        "result_output_artifact_id": result_output_artifact_id,
        "result_output_digest": result_output_digest,
        "result_output_path": result_output_path,
        "result_output_artifact_ids": result_output_artifact_ids,
        "result_output_digests": result_output_digests,
        "result_output_paths": result_output_paths,
        "expected_result_output_artifact_id": expected_result_output,
        "expected_result_output_digest": expected_result_digest,
        "result_observable_id": result_observable_id,
        "result_observable_family": result_observable_family,
        "result_observable_ids": result_observable_ids,
        "result_observable_families": result_observable_families,
        "expected_result_observable_id": expected_result_observable,
        "expected_result_observable_family": expected_result_observable_group,
        "result_row_convention": result_row_convention,
        "result_row_conventions": result_row_conventions,
        "expected_result_row_convention": expected_row_convention,
        "result_normalization_basis": result_normalization_basis,
        "result_normalization_bases": result_normalization_bases,
        "expected_result_normalization_basis": expected_normalization_basis,
        "result_evaluation_method": result_evaluation_method,
        "result_evaluation_methods": result_evaluation_methods,
        "expected_result_evaluation_method": expected_evaluation_method,
        "result_artifact_id": result_artifact_id,
        "result_artifact_ids": result_artifact_ids,
        "expected_result_artifact_id": expected_result_artifact,
        "run_started_at": run_started_text,
        "comsol_version": comsol_version,
        "comsol_versions": comsol_versions,
        "expected_comsol_version": expected_solver_version,
        "timing_breakdown_seconds": timing_breakdown_seconds,
        "timing_breakdown_names": timing_breakdown_names,
        "require_result_provenance": result_provenance_required,
        "require_solver_configuration": solver_configuration_required,
        "require_parameter_set_artifact": parameter_set_required,
        "require_result_table_schema": result_table_schema_required,
        "require_result_output_artifact": result_output_required,
        "require_physics_convention_schema": physics_convention_schema_required,
        "require_result_postprocess_row_convention_schema": (
            result_postprocess_row_convention_schema_required
        ),
        "require_result_component_basis_schema": (
            result_component_basis_schema_required
        ),
        "checks": checks,
        "notes": [
            "Run this before numeric residuals so column position does not become hidden solver knowledge.",
            "A plausible result row can still be unusable if units, sweep axis, source table identity, dataset tag, solution tag, study tag, study step tag, selection tag, or entity dimension were lost.",
            "Dataset tags and solution tags are not enough for replay; record the solution data artifact id/digest separately from the exported table artifact.",
            "For sweep, frequency, and parameter tables, record the independent sweep-axis artifact id/digest and row count separately from both the solution data and the exported table.",
            "When a table row is produced by a parameter study or optimization objective, record the parameter-set artifact id/digest/path and objective observable id/family before reusing it.",
            "Solver configuration is also evidence: record the solver settings artifact id/digest, solver sequence tag, linear solver, and relative tolerance before reusing solver-derived result rows.",
            "The result-table schema id distinguishes the full table layout from scalar or cut-down exports even when columns, units, and output artifact ids look plausible.",
            "The physics convention schema id distinguishes the operator/source/material/coupling/evaluation meaning from the table layout.",
            "The result postprocess-row convention schema id distinguishes selected-row ownership, scalar reduction, and notebook/objective row semantics from both the table schema and the physics convention.",
            "The result component-basis schema id distinguishes component columns, complex representation, and basis/normalization ordering from the table layout and row-reduction convention.",
            "When a table is exported or read by a notebook, the concrete output artifact id, digest, and path should travel with the table metadata.",
            "The result observable id/family says what the table row represents; it is distinct from the output file and from the expression/operator provenance.",
            "The result row convention and normalization basis say how rows should be interpreted before values are reused by a notebook, panel, or cross-validation replay.",
            "The result evaluation method records which API or postprocess path produced the table before values are compared.",
            "Replayable result-table packages should also carry result_artifact_id, run_started_at, solver version, and a compact timing breakdown.",
        ],
    }


def two_port_s_to_yz_equivalent_gate(s11, s21, s12=None, s22=None, z0=50.0, tol=1.0e-9):
    """Convert a two-port S row to Y/Z matrices and equivalent pi/T values.

    This is a public-safe Touchstone post-processing gate for CST and
    ngsolve.bem style port data.  The reference impedance is part of the
    physical contract: the same S row gives different admittance/impedance
    elements when ``z0`` changes.
    """

    z = float(z0)
    tolerance = float(tol)
    if z <= 0.0:
        raise ValueError("z0 must be > 0")
    if tolerance < 0.0:
        raise ValueError("tol must be non-negative")

    S11 = _complex(s11)
    S21 = _complex(s21)
    S12 = S21 if s12 is None else _complex(s12)
    S22 = S11 if s22 is None else _complex(s22)
    inf = complex(float("inf"), float("inf"))

    d_y = (1.0 + S11) * (1.0 + S22) - S12 * S21
    if abs(d_y) <= tolerance:
        y11 = y12 = y21 = y22 = inf
        y_shunt1 = y_shunt2 = y_series = inf
        y_defined = False
    else:
        y11 = ((1.0 - S11) * (1.0 + S22) + S12 * S21) / (z * d_y)
        y12 = (-2.0 * S12) / (z * d_y)
        y21 = (-2.0 * S21) / (z * d_y)
        y22 = ((1.0 + S11) * (1.0 - S22) + S12 * S21) / (z * d_y)
        y_shunt1 = y11 + y12
        y_shunt2 = y22 + y12
        y_series = -y12
        y_defined = True

    d_z = (1.0 - S11) * (1.0 - S22) - S12 * S21
    if abs(d_z) <= tolerance:
        z11 = z12 = z21 = z22 = inf
        z_series1 = z_series2 = z_shunt = inf
        z_defined = False
    else:
        z11 = z * ((1.0 + S11) * (1.0 - S22) + S12 * S21) / d_z
        z12 = z * (2.0 * S12) / d_z
        z21 = z * (2.0 * S21) / d_z
        z22 = z * ((1.0 - S11) * (1.0 + S22) + S12 * S21) / d_z
        z_series1 = z11 - z12
        z_series2 = z22 - z12
        z_shunt = z12
        z_defined = True

    health = two_port_sparameter_health(S11, S21, s12=S12, s22=S22, tol=tolerance)
    checks = {
        "reference_impedance_positive": z > 0.0,
        "sparameter_reciprocity_ok": health["reciprocal"],
        "sparameter_passivity_ok": health["passive"],
        "y_matrix_defined": y_defined,
        "z_matrix_defined": z_defined,
        "y_reciprocal": (abs(y12 - y21) <= tolerance) if y_defined else False,
        "z_reciprocal": (abs(z12 - z21) <= tolerance) if z_defined else False,
    }
    return {
        "policy": "two_port_s_to_yz_equivalent_gate",
        "z0_ohm": z,
        "sparameter_health": health,
        "y11": _complex_row(y11),
        "y12": _complex_row(y12),
        "y21": _complex_row(y21),
        "y22": _complex_row(y22),
        "z11": _complex_row(z11),
        "z12": _complex_row(z12),
        "z21": _complex_row(z21),
        "z22": _complex_row(z22),
        "y_shunt1": _complex_row(y_shunt1),
        "y_shunt2": _complex_row(y_shunt2),
        "y_series": _complex_row(y_series),
        "z_series1": _complex_row(z_series1),
        "z_series2": _complex_row(z_series2),
        "z_shunt": _complex_row(z_shunt),
        "checks": checks,
        "status": "ok" if all(checks.values()) else "needs_attention",
        "tol": tolerance,
    }


def two_port_abcd_cascade_gate(abcd_list, z0=50.0, expect_lossless=False, tol=1.0e-9):
    """Cascade two-port ABCD matrices and check the resulting S-parameters.

    This is a compact public gate for CST/Touchstone post-processing and open
    RF examples: if each section is reciprocal, the determinant of the product
    stays one and the converted S-matrix should remain reciprocal.  Lossless
    cascades may additionally require ``|S11|^2 + |S21|^2 = 1``.
    """

    z = float(z0)
    tolerance = float(tol)
    if z <= 0.0:
        raise ValueError("z0 must be > 0")
    if tolerance < 0.0:
        raise ValueError("tol must be non-negative")
    sections = list(abcd_list)
    if not sections:
        raise ValueError("abcd_list must not be empty")

    A = complex(1.0)
    B = complex(0.0)
    C = complex(0.0)
    D = complex(1.0)
    for index, section in enumerate(sections):
        if isinstance(section, dict):
            try:
                a = _complex(section["A"])
                b = _complex(section["B"])
                c = _complex(section["C"])
                d = _complex(section["D"])
            except KeyError as exc:
                raise ValueError(f"abcd_list[{index}] is missing {exc.args[0]!r}") from exc
        else:
            try:
                a, b, c, d = (_complex(value) for value in section)
            except (TypeError, ValueError) as exc:
                raise ValueError(f"abcd_list[{index}] must contain A,B,C,D") from exc
        A, B, C, D = A * a + B * c, A * b + B * d, C * a + D * c, C * b + D * d

    denominator = A + B / z + C * z + D
    if denominator == 0:
        raise ValueError("cascaded network is singular at the ports")
    s11 = (A + B / z - C * z - D) / denominator
    s21 = 2.0 / denominator
    det = A * D - B * C
    s12 = 2.0 * det / denominator
    s22 = (-A + B / z - C * z + D) / denominator
    health = two_port_sparameter_health(s11, s21, s12=s12, s22=s22, tol=tolerance)
    lossless_power_sum = abs(s11) ** 2 + abs(s21) ** 2
    lossless_error = abs(lossless_power_sum - 1.0)
    checks = {
        "reciprocal_abcd_determinant_ok": abs(det - 1.0) <= tolerance,
        "sparameter_reciprocity_ok": health["reciprocal"],
        "sparameter_passivity_ok": health["passive"],
        "lossless_power_sum_ok": lossless_error <= tolerance if expect_lossless else True,
    }
    return {
        "policy": "two_port_abcd_cascade_gate",
        "n_sections": len(sections),
        "z0_ohm": z,
        "expect_lossless": bool(expect_lossless),
        "A": _complex_row(A),
        "B": _complex_row(B),
        "C": _complex_row(C),
        "D": _complex_row(D),
        "abcd_determinant": _complex_row(det),
        "s11": _complex_row(s11),
        "s21": _complex_row(s21),
        "s12": _complex_row(s12),
        "s22": _complex_row(s22),
        "lossless_power_sum": lossless_power_sum,
        "lossless_power_sum_abs_error": lossless_error,
        "sparameter_health": health,
        "checks": checks,
        "status": "ok" if all(checks.values()) else "needs_attention",
        "tol": tolerance,
    }


def cst_abcd_cascade_solver_ready_manifest_gate(
    artifacts,
    expected_project_id=None,
    expected_run_id=None,
    expected_export_id=None,
    expected_design_frequency_hz=None,
    required_kinds=("port_metadata", "abcd_cascade"),
    frequency_rtol=1.0e-12,
):
    """Check that CST ABCD cascade evidence belongs to one export package.

    This is the manifest-level companion to ``two_port_abcd_cascade_gate``.
    It keeps port metadata and the ordered ABCD cascade tied to one
    project/run/export/design-frequency identity before a de-embedded chain is
    reused by an equivalent-circuit, BEM, or notebook validation.
    """

    rows_in = list(artifacts)
    if not rows_in:
        raise ValueError("artifacts must not be empty")

    def _norm(value):
        return str(value or "").strip().lower().replace("-", "_").replace(" ", "_")

    def _first(row, names):
        for name in names:
            if name in row and row[name] is not None:
                return row[name]
        return None

    def _row_check(row, name):
        checks = row.get("checks", {})
        if isinstance(checks, dict) and name in checks:
            return checks[name]
        return row.get(name)

    required = tuple(_norm(kind) for kind in required_kinds)
    if not required:
        raise ValueError("required_kinds must not be empty")
    tolerance = float(frequency_rtol)
    if tolerance < 0.0:
        raise ValueError("frequency_rtol must be non-negative")

    expected_policies = {
        "port_metadata": {
            "touchstone_port_metadata_gate",
            "cst_touchstone_port_metadata_contract",
        },
        "abcd_cascade": {
            "two_port_abcd_cascade_gate",
            "cst_abcd_cascade_touchstone_contract",
        },
        "touchstone_row": {
            "touchstone_row_solver_ready_preflight_gate",
            "two_port_sparameter_health",
        },
    }

    details = []
    kind_counts = {}
    project_ids = []
    run_ids = []
    export_ids = []
    design_frequencies = []
    z0_values = []
    missing_project_id = []
    missing_run_id = []
    missing_export_id = []
    missing_frequency = []
    missing_paths = []
    bad_source_tool = []
    unknown_kinds = []
    bad_upstream_status = []
    bad_upstream_policy = []
    missing_port_metadata = []
    missing_abcd_metadata = []
    bad_abcd_determinant = []
    nonreciprocal_abcd = []
    nonpassive_abcd = []
    nonlossless_abcd = []
    nonpassive_touchstone_rows = []
    nonreciprocal_touchstone_rows = []

    for index, row in enumerate(rows_in, start=1):
        if not isinstance(row, dict):
            raise ValueError("each artifact must be a dictionary")
        kind = _norm(_first(row, ("kind", "artifact_kind", "type")))
        project_id = _first(row, ("project_id", "cst_project_id", "model_id"))
        run_id = _first(row, ("run_id", "solver_run_id", "simulation_id"))
        export_id = _first(row, ("export_id", "result_id", "dataset_id"))
        source_tool = _first(row, ("source_tool", "tool", "source"))
        source_tool_norm = _norm(source_tool)
        path = _first(row, ("path", "file", "artifact_path", "table_path"))
        gate_policy = _first(row, ("gate_policy", "policy", "validator"))
        gate_policy_norm = _norm(gate_policy)
        status = _first(row, ("status", "gate_status", "validation_status"))
        status_norm = _norm(status)
        pass_flag = bool(row.get("pass", False))
        frequency = _first(row, ("design_frequency_hz", "frequency_hz", "freq_hz"))

        if kind:
            kind_counts[kind] = kind_counts.get(kind, 0) + 1
        else:
            unknown_kinds.append(index)
        if kind and kind not in expected_policies:
            unknown_kinds.append(kind)
        if not project_id:
            missing_project_id.append(index)
        else:
            project_ids.append(str(project_id))
        if not run_id:
            missing_run_id.append(index)
        else:
            run_ids.append(str(run_id))
        if not export_id:
            missing_export_id.append(index)
        else:
            export_ids.append(str(export_id))
        if source_tool_norm not in {"cst", "cst_studio", "cst_studio_suite"}:
            bad_source_tool.append({"index": index, "kind": kind, "source_tool": source_tool})
        if not path:
            missing_paths.append(index)
        if frequency is None:
            missing_frequency.append(index)
        else:
            design_frequencies.append(float(frequency))
        if not (pass_flag or status_norm in {"ok", "pass", "passed", "verified"}):
            bad_upstream_status.append({"index": index, "kind": kind, "status": status})
        if kind in expected_policies and gate_policy_norm not in expected_policies[kind]:
            bad_upstream_policy.append({
                "index": index,
                "kind": kind,
                "gate_policy": gate_policy,
                "expected": sorted(expected_policies[kind]),
            })

        if kind == "port_metadata":
            z0 = _first(row, ("reference_impedance_ohm", "z0_ohm", "z0"))
            if z0 is not None:
                z0_values.append(float(z0))
            if not (row.get("data_format") and z0 is not None and row.get("port_order")):
                missing_port_metadata.append(index)
        if kind == "abcd_cascade":
            z0 = _first(row, ("z0_ohm", "reference_impedance_ohm", "z0"))
            n_sections = _first(row, ("n_sections", "section_count"))
            expect_lossless = bool(row.get("expect_lossless", False))
            if z0 is not None:
                z0_values.append(float(z0))
            if z0 is None or not n_sections or int(n_sections) < 1:
                missing_abcd_metadata.append(index)
            if _row_check(row, "reciprocal_abcd_determinant_ok") is not True:
                bad_abcd_determinant.append(index)
            if _row_check(row, "sparameter_reciprocity_ok") is not True:
                nonreciprocal_abcd.append(index)
            if _row_check(row, "sparameter_passivity_ok") is not True:
                nonpassive_abcd.append(index)
            if expect_lossless and _row_check(row, "lossless_power_sum_ok") is not True:
                nonlossless_abcd.append(index)
        if kind == "touchstone_row":
            if _row_check(row, "sparameter_passivity_ok") is not True:
                nonpassive_touchstone_rows.append(index)
            if _row_check(row, "sparameter_reciprocity_ok") is not True:
                nonreciprocal_touchstone_rows.append(index)

        details.append({
            "index": index,
            "kind": kind,
            "project_id": None if project_id is None else str(project_id),
            "run_id": None if run_id is None else str(run_id),
            "export_id": None if export_id is None else str(export_id),
            "source_tool": source_tool,
            "path": path,
            "gate_policy": gate_policy,
            "status": status,
            "pass": pass_flag,
            "design_frequency_hz": None if frequency is None else float(frequency),
        })

    unique_project_ids = sorted(set(project_ids))
    unique_run_ids = sorted(set(run_ids))
    unique_export_ids = sorted(set(export_ids))
    unique_design_frequencies = sorted(set(design_frequencies))
    unique_z0_values = sorted(set(z0_values))
    required_set = set(required)
    present_set = set(kind_counts)
    checks = {
        "required_kinds_present": required_set.issubset(present_set),
        "no_unknown_kinds": not unknown_kinds,
        "project_ids_present": not missing_project_id,
        "project_ids_unique": len(unique_project_ids) == 1,
        "run_ids_present": not missing_run_id,
        "run_ids_unique": len(unique_run_ids) == 1,
        "export_ids_present": not missing_export_id,
        "export_ids_unique": len(unique_export_ids) == 1,
        "source_tool_is_cst": not bad_source_tool,
        "paths_present": not missing_paths,
        "design_frequencies_present": not missing_frequency,
        "design_frequencies_unique": len(unique_design_frequencies) == 1,
        "z0_recorded": bool(unique_z0_values),
        "z0_consistent": len(unique_z0_values) == 1,
        "upstream_gate_status_ok": not bad_upstream_status,
        "upstream_gate_policy_known": not bad_upstream_policy,
        "port_metadata_complete": not missing_port_metadata,
        "abcd_metadata_complete": not missing_abcd_metadata,
        "abcd_determinant_reciprocal": not bad_abcd_determinant,
        "abcd_sparameters_reciprocal": not nonreciprocal_abcd,
        "abcd_sparameters_passive": not nonpassive_abcd,
        "abcd_lossless_when_requested": not nonlossless_abcd,
        "touchstone_rows_passive": not nonpassive_touchstone_rows,
        "touchstone_rows_reciprocal": not nonreciprocal_touchstone_rows,
    }
    if expected_project_id is not None:
        checks["expected_project_id_matches"] = unique_project_ids == [str(expected_project_id)]
    if expected_run_id is not None:
        checks["expected_run_id_matches"] = unique_run_ids == [str(expected_run_id)]
    if expected_export_id is not None:
        checks["expected_export_id_matches"] = unique_export_ids == [str(expected_export_id)]
    if expected_design_frequency_hz is not None:
        expected = float(expected_design_frequency_hz)
        checks["expected_design_frequency_matches"] = (
            len(unique_design_frequencies) == 1
            and abs(unique_design_frequencies[0] - expected) <= max(1.0e-9, abs(expected) * tolerance)
        )

    return {
        "policy": "cst_abcd_cascade_solver_ready_manifest_gate",
        "required_kinds": list(required),
        "present_kinds": dict(sorted(kind_counts.items())),
        "project_ids": unique_project_ids,
        "run_ids": unique_run_ids,
        "export_ids": unique_export_ids,
        "design_frequencies_hz": unique_design_frequencies,
        "z0_values_ohm": unique_z0_values,
        "expected_project_id": None if expected_project_id is None else str(expected_project_id),
        "expected_run_id": None if expected_run_id is None else str(expected_run_id),
        "expected_export_id": None if expected_export_id is None else str(expected_export_id),
        "expected_design_frequency_hz": None if expected_design_frequency_hz is None else float(expected_design_frequency_hz),
        "missing_project_id_rows": missing_project_id,
        "missing_run_id_rows": missing_run_id,
        "missing_export_id_rows": missing_export_id,
        "missing_frequency_rows": missing_frequency,
        "bad_source_tool_rows": bad_source_tool,
        "missing_path_rows": missing_paths,
        "unknown_kinds": unknown_kinds,
        "bad_upstream_status_rows": bad_upstream_status,
        "bad_upstream_policy_rows": bad_upstream_policy,
        "missing_port_metadata_rows": missing_port_metadata,
        "missing_abcd_metadata_rows": missing_abcd_metadata,
        "bad_abcd_determinant_rows": bad_abcd_determinant,
        "nonreciprocal_abcd_rows": nonreciprocal_abcd,
        "nonpassive_abcd_rows": nonpassive_abcd,
        "nonlossless_abcd_rows": nonlossless_abcd,
        "nonpassive_touchstone_rows": nonpassive_touchstone_rows,
        "nonreciprocal_touchstone_rows": nonreciprocal_touchstone_rows,
        "artifacts": details,
        "checks": checks,
        "status": "ok" if all(checks.values()) else "needs_attention",
        "version_note": (
            "Run after ABCD cascade and port-metadata sub-gates so de-embedded "
            "two-port chains cannot mix project/run/export identity, reference "
            "impedance, or active/passivity-failed rows."
        ),
    }


def mesh_import_quality_manifest_gate(
    surface_element_types=("tri",),
    volume_element_types=("tet",),
    order=1,
    min_scaled_jacobian_before=None,
    min_scaled_jacobian_after=1.0,
    min_scaled_jacobian_threshold=0.1,
    negative_jacobian_count_before=0,
    negative_jacobian_count_after=0,
    cad_connectivity_recorded=True,
    cad_compliance_recorded=True,
    boundary_conformity_tolerance=1.0e-8,
    max_boundary_distance=0.0,
    source="third_party_mesh",
    mesh_format=".vol",
):
    """Check a mesh-import quality manifest before solver or notebook reuse.

    The public lesson is intentionally conservative: record CAD
    connectivity/conformity/compliance and element validity, but do not silently
    convert unsupported high-order or non-tri/tet topology into a first-order
    teaching lane.
    """

    def _norm_kind(value):
        key = str(value or "").strip().lower()
        aliases = {
            "triangle": "tri",
            "tri3": "tri",
            "tetra": "tet",
            "tetrahedron": "tet",
            "tet4": "tet",
        }
        return aliases.get(key, key)

    surface = [_norm_kind(value) for value in surface_element_types]
    volume = [_norm_kind(value) for value in volume_element_types]
    polynomial_order = int(order) if float(order).is_integer() else float(order)
    min_j_after = float(min_scaled_jacobian_after)
    min_j_threshold = float(min_scaled_jacobian_threshold)
    neg_before = int(negative_jacobian_count_before)
    neg_after = int(negative_jacobian_count_after)
    tolerance = float(boundary_conformity_tolerance)
    max_distance = float(max_boundary_distance)
    if min_j_threshold < -1.0:
        raise ValueError("min_scaled_jacobian_threshold is unexpectedly low")
    if neg_before < 0 or neg_after < 0:
        raise ValueError("negative Jacobian counts must be non-negative")
    if tolerance < 0.0 or max_distance < 0.0:
        raise ValueError("boundary distances/tolerances must be non-negative")

    if min_scaled_jacobian_before is None:
        min_j_before = None
        improvement = None
    else:
        min_j_before = float(min_scaled_jacobian_before)
        improvement = min_j_after - min_j_before

    checks = {
        "source_recorded": bool(str(source).strip()),
        "format_recorded": bool(str(mesh_format).strip()),
        "surface_triangles_only": bool(surface) and all(kind == "tri" for kind in surface),
        "volume_tetrahedra_only": bool(volume) and all(kind == "tet" for kind in volume),
        "first_order_only": polynomial_order == 1,
        "no_negative_jacobian_after": neg_after == 0,
        "negative_jacobian_count_improved": neg_after <= neg_before,
        "min_scaled_jacobian_recorded": math.isfinite(min_j_after),
        "min_scaled_jacobian_above_threshold": min_j_after >= min_j_threshold,
        "cad_connectivity_recorded": bool(cad_connectivity_recorded),
        "cad_compliance_recorded": bool(cad_compliance_recorded),
        "boundary_conformity_within_tolerance": max_distance <= tolerance,
    }
    checks["first_order_tri_tet_policy_honored"] = (
        checks["surface_triangles_only"]
        and checks["volume_tetrahedra_only"]
        and checks["first_order_only"]
    )
    issues = [name for name, ok in checks.items() if not ok]

    return {
        "policy": "mesh_import_quality_manifest_gate",
        "status": "ok" if not issues else "needs_attention",
        "source": str(source),
        "mesh_format": str(mesh_format),
        "surface_element_types": surface,
        "volume_element_types": volume,
        "order": polynomial_order,
        "min_scaled_jacobian_before": min_j_before,
        "min_scaled_jacobian_after": min_j_after,
        "min_scaled_jacobian_threshold": min_j_threshold,
        "min_scaled_jacobian_improvement": improvement,
        "negative_jacobian_count_before": neg_before,
        "negative_jacobian_count_after": neg_after,
        "boundary_conformity_tolerance": tolerance,
        "max_boundary_distance": max_distance,
        "checks": checks,
        "issues": issues,
        "notes": [
            "Use this before physics residuals: invalid or unsupported mesh topology is not a solver failure.",
            "High-order/hex/prism evidence can be useful engineering context while still being rejected by a first-order tri/tet lane.",
        ],
    }


def netgen_vol_first_order_fem_bem_trace_package_handoff(
    package,
    expected_mesh_id=None,
    expected_export_id=None,
    expected_trace_artifact_id=None,
    expected_surface_mesh_id=None,
    expected_source_file_id=None,
    expected_source_format=".vol",
    expected_coupling_kind=None,
    expected_formulation_id=None,
    expected_bem_kernel_family=None,
    expected_coupling_convention_schema_id=None,
    expected_fem_bem_postprocess_row_convention_schema_id=None,
    expected_trace_basis_schema_id=None,
    expected_assembly_rule_id=None,
    expected_quadrature_rule_id=None,
    expected_volume_space=None,
    expected_surface_space=None,
    expected_boundary_numbers=None,
    expected_boundary_names=None,
    expected_boundary_row_identity=None,
    expected_trace_operator_artifact_id=None,
    expected_trace_operator_policy=None,
    expected_trace_output_artifact_id=None,
    expected_trace_output_digest=None,
    expected_trace_observable_id=None,
    expected_trace_observable_family=None,
    expected_coupled_system_artifact_id=None,
    expected_coupled_system_digest=None,
    expected_result_artifact_id=None,
    expected_matlab_version=None,
    expected_linear_solver_report_artifact_id=None,
    expected_linear_solver_report_digest=None,
    expected_linear_solver_name=None,
    expected_linear_solver_tolerance=None,
    expected_linear_solver_residual_norm_max=None,
    expected_parameter_set_artifact_id=None,
    expected_parameter_set_digest=None,
    expected_parameter_set_path=None,
    expected_objective_observable_id=None,
    expected_objective_observable_family=None,
    require_result_provenance=False,
    require_trace_output_artifact=False,
    require_linear_solver_report=False,
    require_parameter_set_artifact=False,
    require_coupling_convention_schema=False,
    require_fem_bem_postprocess_row_convention_schema=False,
    require_trace_basis_schema=False,
    tol=1.0e-12,
):
    """Check a readable Netgen .vol FEM/BEM trace handoff package.

    This public-safe gate encodes the MATLAB/Gypsilab/Lukas mesh policy: a
    first-order Netgen ``.vol`` package should keep volume tetrahedra,
    boundary triangles, and the H1-to-boundary trace gather together with the
    same one-based node ids.  The gate does not run MATLAB; it replays the
    package contract from archived JSON/notebook results.
    """

    if not isinstance(package, dict):
        raise ValueError("package must be a dictionary")

    def _norm(value):
        return str(value or "").strip().lower().replace("-", "_").replace(" ", "_")

    def _as_list(value):
        if value is None:
            return []
        if isinstance(value, (list, tuple)):
            return list(value)
        return [value]

    def _int_list(value):
        return [int(item) for item in _as_list(value)]

    def _string_list(value):
        return [str(item).strip() for item in _as_list(value) if str(item).strip()]

    def _maybe_int(value):
        if value is None:
            return None
        try:
            return int(value)
        except (TypeError, ValueError):
            return None

    def _maybe_float(value):
        if value is None:
            return None
        try:
            numeric = float(value)
        except (TypeError, ValueError):
            return None
        if not math.isfinite(numeric):
            return None
        return numeric

    def _record_value(record, *names):
        if not isinstance(record, dict):
            return None
        for name in names:
            if name in record:
                return record[name]
        return None

    def _record_values(record, *names):
        if not isinstance(record, dict):
            return []
        return [record[name] for name in names if name in record and record[name] is not None]

    def _first_mapping(value):
        if isinstance(value, dict):
            return value
        for item in _as_list(value):
            if isinstance(item, dict):
                return item
        return {}

    def _unique_strings(values):
        return sorted({str(value) for value in values if str(value)})

    def _trace_row_identity_records(value):
        if value is None:
            return []
        if isinstance(value, dict):
            row_indices = _as_list(
                _record_value(value, "trace_row_index", "row_index", "row")
            )
            fem_ids = _as_list(_record_value(value, "fem_node_id", "fem_node"))
            bem_ids = _as_list(_record_value(value, "bem_node_id", "bem_node"))
            surface_ids = _as_list(
                _record_value(value, "surface_node_index", "surface_node_id")
            )
            count = max(
                len(row_indices),
                len(fem_ids),
                len(bem_ids),
                len(surface_ids),
                0,
            )
            records = []
            for index in range(count):
                records.append(
                    {
                        "trace_row_index": row_indices[index]
                        if index < len(row_indices)
                        else None,
                        "fem_node_id": fem_ids[index] if index < len(fem_ids) else None,
                        "bem_node_id": bem_ids[index] if index < len(bem_ids) else None,
                        "surface_node_index": surface_ids[index]
                        if index < len(surface_ids)
                        else None,
                    }
                )
            return records
        return [record for record in _as_list(value) if isinstance(record, dict)]

    def _boundary_row_identity_records(value):
        if value is None:
            return []
        if isinstance(value, dict):
            row_indices = _as_list(
                _record_value(value, "surface_triangle_index", "triangle_index", "row_index", "row")
            )
            triangles = _as_list(
                _record_value(value, "surface_triangle_nodes", "triangle_nodes", "nodes")
            )
            numbers = _as_list(_record_value(value, "boundary_number", "boundary_id", "bc_id"))
            names = _as_list(_record_value(value, "boundary_name", "bc_name", "name"))
            adjacent = _as_list(_record_value(value, "adjacent_tet_index", "adjacent_tet", "tet_index"))
            count = max(len(row_indices), len(triangles), len(numbers), len(names), len(adjacent), 0)
            records = []
            for index in range(count):
                records.append(
                    {
                        "surface_triangle_index": row_indices[index] if index < len(row_indices) else None,
                        "surface_triangle_nodes": triangles[index] if index < len(triangles) else None,
                        "boundary_number": numbers[index] if index < len(numbers) else None,
                        "boundary_name": names[index] if index < len(names) else None,
                        "adjacent_tet_index": adjacent[index] if index < len(adjacent) else None,
                    }
                )
            return records
        return [record for record in _as_list(value) if isinstance(record, dict)]

    def _normalize_boundary_row_identity(records):
        normalized = []
        for record in records:
            triangle = _record_value(record, "surface_triangle_nodes", "triangle_nodes", "nodes")
            normalized.append(
                {
                    "surface_triangle_index": _maybe_int(
                        _record_value(record, "surface_triangle_index", "triangle_index", "row_index", "row")
                    ),
                    "surface_triangle_nodes": _int_list(triangle),
                    "boundary_number": _maybe_int(_record_value(record, "boundary_number", "boundary_id", "bc_id")),
                    "boundary_name": str(_record_value(record, "boundary_name", "bc_name", "name") or ""),
                    "adjacent_tet_index": _maybe_int(
                        _record_value(record, "adjacent_tet_index", "adjacent_tet", "tet_index")
                    ),
                }
            )
        return normalized

    mesh_id = package.get("mesh_id") or package.get("vol_id") or package.get("grid_id")
    export_id = package.get("export_id") or package.get("handoff_id") or package.get("dataset_id")
    source_path = package.get("source_path") or package.get("vol_path") or package.get("path")
    source_format = package.get("source_format") or package.get("mesh_format") or package.get("format")
    source_file_id = (
        package.get("source_file_id")
        or package.get("sourceFileId")
        or package.get("vol_file_id")
        or package.get("vol_sha256")
        or package.get("source_file_sha256")
        or package.get("source_sha256")
    )
    policy = package.get("policy")
    trace_artifact_id = (
        package.get("trace_artifact_id")
        or package.get("trace_matrix_artifact_id")
        or package.get("trace_id")
    )
    surface_mesh_id = (
        package.get("surface_mesh_id")
        or package.get("boundary_mesh_id")
        or package.get("bem_surface_mesh_id")
    )
    polynomial_order = int(package.get("polynomial_order") or package.get("order") or 1)
    curved_element_count = int(
        package.get("curved_element_count")
        or package.get("curvedelements_count")
        or package.get("high_order_element_count")
        or 0
    )
    geo = _first_mapping(package.get("geo") or package.get("lukas"))
    gypsilab = _first_mapping(package.get("gypsilab"))
    trace = _first_mapping(package.get("trace"))
    operators = _first_mapping(package.get("operators"))
    operator_trace = _first_mapping(operators.get("trace") if isinstance(operators, dict) else None)
    operator_bem = _first_mapping(operators.get("bem") if isinstance(operators, dict) else None)
    execution = _first_mapping(
        package.get("execution")
        or package.get("run")
        or package.get("result_provenance")
        or package.get("provenance")
    )
    result_artifact_id = (
        package.get("result_artifact_id")
        or package.get("resultArtifactId")
        or package.get("result_id")
        or execution.get("result_artifact_id")
        or execution.get("resultArtifactId")
        or execution.get("result_id")
    )
    run_started_at = (
        package.get("run_started_at")
        or package.get("runStartedAt")
        or package.get("run_date")
        or package.get("runDate")
        or execution.get("run_started_at")
        or execution.get("runStartedAt")
        or execution.get("run_date")
        or execution.get("runDate")
    )
    matlab_version = (
        package.get("matlab_version")
        or package.get("matlabVersion")
        or execution.get("matlab_version")
        or execution.get("matlabVersion")
    )
    timing_breakdown = _first_mapping(
        package.get("timing_breakdown")
        or package.get("timingBreakdown")
        or package.get("timings")
        or execution.get("timing_breakdown")
        or execution.get("timingBreakdown")
        or execution.get("timings")
    )
    timing_breakdown_seconds = {}
    for name, value in timing_breakdown.items():
        try:
            timing_breakdown_seconds[str(name)] = float(value)
        except (TypeError, ValueError):
            timing_breakdown_seconds[str(name)] = math.nan
    timing_breakdown_names = sorted(timing_breakdown_seconds)
    linear_solver_report = _first_mapping(
        package.get("linear_solver_report")
        or package.get("linearSolverReport")
        or package.get("solver_report")
        or package.get("solverReport")
        or execution.get("linear_solver_report")
        or execution.get("linearSolverReport")
        or execution.get("solver_report")
        or execution.get("solverReport")
        or operators.get("linear_solver_report")
        or operators.get("linearSolverReport")
        or operators.get("solver_report")
        or operators.get("solverReport")
    )
    linear_solver_report_artifact_id_names = (
        "linear_solver_report_artifact_id",
        "linearSolverReportArtifactId",
        "solver_report_artifact_id",
        "solverReportArtifactId",
        "solve_report_artifact_id",
        "solveReportArtifactId",
    )
    linear_solver_report_digest_names = (
        "linear_solver_report_digest",
        "linearSolverReportDigest",
        "linear_solver_report_sha256",
        "linearSolverReportSha256",
        "solver_report_digest",
        "solverReportDigest",
        "solver_report_sha256",
        "solverReportSha256",
        "solve_report_digest",
        "solveReportDigest",
    )
    linear_solver_name_names = (
        "linear_solver_name",
        "linearSolverName",
        "solver_name",
        "solverName",
        "solve_method",
        "solveMethod",
    )
    linear_solver_tolerance_names = (
        "linear_solver_tolerance",
        "linearSolverTolerance",
        "solver_tolerance",
        "solverTolerance",
        "residual_tolerance",
        "residualTolerance",
    )
    linear_solver_residual_norm_names = (
        "linear_solver_residual_norm",
        "linearSolverResidualNorm",
        "solver_residual_norm",
        "solverResidualNorm",
        "relative_residual_norm",
        "relativeResidualNorm",
        "residual_norm",
        "residualNorm",
    )
    linear_solver_iteration_count_names = (
        "linear_solver_iteration_count",
        "linearSolverIterationCount",
        "solver_iteration_count",
        "solverIterationCount",
        "iteration_count",
        "iterationCount",
    )
    linear_solver_report_artifact_id_values = (
        _record_values(package, *linear_solver_report_artifact_id_names)
        + _record_values(execution, *linear_solver_report_artifact_id_names)
        + _record_values(operators, *linear_solver_report_artifact_id_names)
        + _record_values(linear_solver_report, *linear_solver_report_artifact_id_names)
    )
    linear_solver_report_digest_values = (
        _record_values(package, *linear_solver_report_digest_names)
        + _record_values(execution, *linear_solver_report_digest_names)
        + _record_values(operators, *linear_solver_report_digest_names)
        + _record_values(linear_solver_report, *linear_solver_report_digest_names)
    )
    linear_solver_name_values = (
        _record_values(package, *linear_solver_name_names)
        + _record_values(execution, *linear_solver_name_names)
        + _record_values(operators, *linear_solver_name_names)
        + _record_values(linear_solver_report, *linear_solver_name_names)
    )
    linear_solver_tolerance_values = [
        value
        for value in (
            _record_values(package, *linear_solver_tolerance_names)
            + _record_values(execution, *linear_solver_tolerance_names)
            + _record_values(operators, *linear_solver_tolerance_names)
            + _record_values(linear_solver_report, *linear_solver_tolerance_names)
        )
        if _maybe_float(value) is not None
    ]
    linear_solver_residual_norm_values = [
        value
        for value in (
            _record_values(package, *linear_solver_residual_norm_names)
            + _record_values(execution, *linear_solver_residual_norm_names)
            + _record_values(operators, *linear_solver_residual_norm_names)
            + _record_values(linear_solver_report, *linear_solver_residual_norm_names)
        )
        if _maybe_float(value) is not None
    ]
    linear_solver_iteration_count_values = [
        value
        for value in (
            _record_values(package, *linear_solver_iteration_count_names)
            + _record_values(execution, *linear_solver_iteration_count_names)
            + _record_values(operators, *linear_solver_iteration_count_names)
            + _record_values(linear_solver_report, *linear_solver_iteration_count_names)
        )
        if _maybe_int(value) is not None
    ]
    linear_solver_report_artifact_ids = _unique_strings(linear_solver_report_artifact_id_values)
    linear_solver_report_digests = _unique_strings(linear_solver_report_digest_values)
    linear_solver_names = sorted({_norm(value) for value in linear_solver_name_values if str(value)})
    linear_solver_tolerances = sorted({_maybe_float(value) for value in linear_solver_tolerance_values})
    linear_solver_residual_norms = sorted({_maybe_float(value) for value in linear_solver_residual_norm_values})
    linear_solver_iteration_counts = sorted({_maybe_int(value) for value in linear_solver_iteration_count_values})
    linear_solver_report_artifact_id = (
        linear_solver_report_artifact_ids[0] if linear_solver_report_artifact_ids else None
    )
    linear_solver_report_digest = linear_solver_report_digests[0] if linear_solver_report_digests else None
    linear_solver_name = linear_solver_names[0] if linear_solver_names else None
    linear_solver_tolerance = linear_solver_tolerances[0] if linear_solver_tolerances else None
    linear_solver_residual_norm = linear_solver_residual_norms[0] if linear_solver_residual_norms else None
    linear_solver_iteration_count = linear_solver_iteration_counts[0] if linear_solver_iteration_counts else None
    optimization = _first_mapping(
        package.get("optimization")
        or package.get("optimizer")
        or package.get("objective")
        or package.get("design")
        or execution.get("optimization")
        or execution.get("optimizer")
        or execution.get("objective")
        or execution.get("design")
    )
    parameter_set_artifact_id_names = (
        "parameter_set_artifact_id",
        "parameterSetArtifactId",
        "initial_value_artifact_id",
        "initialValueArtifactId",
        "design_parameter_artifact_id",
        "designParameterArtifactId",
        "optimization_parameter_set_artifact_id",
        "optimizationParameterSetArtifactId",
    )
    parameter_set_digest_names = (
        "parameter_set_digest",
        "parameterSetDigest",
        "parameter_set_sha256",
        "parameterSetSha256",
        "initial_value_digest",
        "initialValueDigest",
        "design_parameter_digest",
        "designParameterDigest",
        "optimization_parameter_set_digest",
        "optimizationParameterSetDigest",
    )
    parameter_set_path_names = (
        "parameter_set_path",
        "parameterSetPath",
        "initial_value_path",
        "initialValuePath",
        "design_parameter_path",
        "designParameterPath",
        "optimization_parameter_set_path",
        "optimizationParameterSetPath",
    )
    objective_observable_id_names = (
        "objective_observable_id",
        "objectiveObservableId",
        "objective_artifact_id",
        "objectiveArtifactId",
        "optimization_objective_id",
        "optimizationObjectiveId",
    )
    objective_observable_family_names = (
        "objective_observable_family",
        "objectiveObservableFamily",
        "objective_family",
        "objectiveFamily",
        "optimization_objective_family",
        "optimizationObjectiveFamily",
    )
    parameter_set_artifact_id_values = (
        _record_values(package, *parameter_set_artifact_id_names)
        + _record_values(execution, *parameter_set_artifact_id_names)
        + _record_values(optimization, *parameter_set_artifact_id_names)
    )
    parameter_set_digest_values = (
        _record_values(package, *parameter_set_digest_names)
        + _record_values(execution, *parameter_set_digest_names)
        + _record_values(optimization, *parameter_set_digest_names)
    )
    parameter_set_path_values = (
        _record_values(package, *parameter_set_path_names)
        + _record_values(execution, *parameter_set_path_names)
        + _record_values(optimization, *parameter_set_path_names)
    )
    objective_observable_id_values = (
        _record_values(package, *objective_observable_id_names)
        + _record_values(execution, *objective_observable_id_names)
        + _record_values(optimization, *objective_observable_id_names)
    )
    objective_observable_family_values = (
        _record_values(package, *objective_observable_family_names)
        + _record_values(execution, *objective_observable_family_names)
        + _record_values(optimization, *objective_observable_family_names)
    )
    parameter_set_artifact_ids = _unique_strings(parameter_set_artifact_id_values)
    parameter_set_digests = _unique_strings(parameter_set_digest_values)
    parameter_set_paths = _unique_strings(parameter_set_path_values)
    objective_observable_ids = _unique_strings(objective_observable_id_values)
    objective_observable_families = sorted(
        {_norm(value) for value in objective_observable_family_values if str(value)}
    )
    parameter_set_artifact_id = (
        parameter_set_artifact_ids[0] if parameter_set_artifact_ids else None
    )
    parameter_set_digest = parameter_set_digests[0] if parameter_set_digests else None
    parameter_set_path = parameter_set_paths[0] if parameter_set_paths else None
    objective_observable_id = objective_observable_ids[0] if objective_observable_ids else None
    objective_observable_family = (
        objective_observable_families[0] if objective_observable_families else None
    )
    trace_source_file_id = (
        trace.get("source_file_id")
        or trace.get("sourceFileId")
        or trace.get("vol_file_id")
        or trace.get("vol_sha256")
        or trace.get("source_file_sha256")
        or trace.get("source_sha256")
    )
    trace_artifact_id = (
        trace_artifact_id
        or trace.get("trace_artifact_id")
        or trace.get("trace_matrix_artifact_id")
        or trace.get("artifact_id")
    )
    surface_mesh_id = (
        surface_mesh_id
        or trace.get("surface_mesh_id")
        or trace.get("boundary_mesh_id")
        or trace.get("bem_surface_mesh_id")
        or gypsilab.get("surface_mesh_id")
        or gypsilab.get("boundary_mesh_id")
    )
    trace_operator_artifact_id = (
        package.get("trace_operator_artifact_id")
        or package.get("operator_trace_artifact_id")
        or package.get("trace_assembly_artifact_id")
        or trace.get("trace_operator_artifact_id")
        or trace.get("operator_trace_artifact_id")
        or trace.get("trace_assembly_artifact_id")
        or operator_trace.get("trace_operator_artifact_id")
        or operator_trace.get("traceOperatorArtifactId")
        or operator_trace.get("operator_artifact_id")
    )
    trace_operator_policy = (
        package.get("trace_operator_policy")
        or package.get("operator_trace_policy")
        or package.get("trace_assembly_policy")
        or trace.get("trace_operator_policy")
        or trace.get("operator_trace_policy")
        or trace.get("trace_assembly_policy")
        or operator_trace.get("trace_operator_policy")
        or operator_trace.get("traceOperatorPolicy")
        or operator_trace.get("operator_policy")
    )
    trace_output_artifact_id_names = (
        "trace_output_artifact_id",
        "traceOutputArtifactId",
        "output_artifact_id",
        "outputArtifactId",
        "trace_matrix_output_artifact_id",
        "traceMatrixOutputArtifactId",
        "postprocess_output_artifact_id",
    )
    trace_output_digest_names = (
        "trace_output_digest",
        "traceOutputDigest",
        "trace_output_sha256",
        "output_digest",
        "outputDigest",
        "output_sha256",
        "trace_matrix_output_digest",
        "traceMatrixOutputDigest",
        "postprocess_output_digest",
    )
    trace_observable_id_names = (
        "trace_observable_id",
        "traceObservableId",
        "fem_bem_trace_observable_id",
        "observable_id",
        "observableId",
    )
    trace_observable_family_names = (
        "trace_observable_family",
        "traceObservableFamily",
        "fem_bem_trace_observable_family",
        "observable_family",
        "observableFamily",
    )
    coupled_system_artifact_id_names = (
        "coupled_system_artifact_id",
        "coupledSystemArtifactId",
        "fem_bem_system_artifact_id",
        "femBemSystemArtifactId",
        "fem_bem_coupled_system_artifact_id",
        "femBemCoupledSystemArtifactId",
        "schur_system_artifact_id",
        "schurSystemArtifactId",
        "linear_system_artifact_id",
        "linearSystemArtifactId",
    )
    coupled_system_digest_names = (
        "coupled_system_digest",
        "coupledSystemDigest",
        "coupled_system_sha256",
        "coupledSystemSha256",
        "fem_bem_system_digest",
        "femBemSystemDigest",
        "fem_bem_system_sha256",
        "femBemSystemSha256",
        "schur_system_digest",
        "schurSystemDigest",
        "schur_system_sha256",
        "schurSystemSha256",
        "linear_system_digest",
        "linearSystemDigest",
    )
    trace_output_path_names = (
        "trace_output_path",
        "traceOutputPath",
        "output_path",
        "outputPath",
        "trace_matrix_output_path",
        "traceMatrixOutputPath",
        "postprocess_output_path",
    )
    trace_output_artifact_id_values = (
        _record_values(package, *trace_output_artifact_id_names)
        + _record_values(trace, *trace_output_artifact_id_names)
        + _record_values(operator_trace, *trace_output_artifact_id_names)
    )
    trace_output_digest_values = (
        _record_values(package, *trace_output_digest_names)
        + _record_values(trace, *trace_output_digest_names)
        + _record_values(operator_trace, *trace_output_digest_names)
    )
    trace_observable_id_values = (
        _record_values(package, *trace_observable_id_names)
        + _record_values(trace, *trace_observable_id_names)
        + _record_values(operator_trace, *trace_observable_id_names)
    )
    trace_observable_family_values = (
        _record_values(package, *trace_observable_family_names)
        + _record_values(trace, *trace_observable_family_names)
        + _record_values(operator_trace, *trace_observable_family_names)
    )
    trace_output_path_values = (
        _record_values(package, *trace_output_path_names)
        + _record_values(trace, *trace_output_path_names)
        + _record_values(operator_trace, *trace_output_path_names)
    )
    coupled_system_artifact_id_values = (
        _record_values(package, *coupled_system_artifact_id_names)
        + _record_values(trace, *coupled_system_artifact_id_names)
        + _record_values(operators, *coupled_system_artifact_id_names)
        + _record_values(operator_bem, *coupled_system_artifact_id_names)
    )
    coupled_system_digest_values = (
        _record_values(package, *coupled_system_digest_names)
        + _record_values(trace, *coupled_system_digest_names)
        + _record_values(operators, *coupled_system_digest_names)
        + _record_values(operator_bem, *coupled_system_digest_names)
    )
    trace_output_artifact_ids = _unique_strings(trace_output_artifact_id_values)
    trace_output_digests = _unique_strings(trace_output_digest_values)
    trace_observable_ids = _unique_strings(trace_observable_id_values)
    trace_observable_families = sorted({_norm(value) for value in trace_observable_family_values if str(value)})
    trace_output_paths = _unique_strings(trace_output_path_values)
    coupled_system_artifact_ids = _unique_strings(coupled_system_artifact_id_values)
    coupled_system_digests = _unique_strings(coupled_system_digest_values)
    trace_output_artifact_id = trace_output_artifact_ids[0] if trace_output_artifact_ids else None
    trace_output_digest = trace_output_digests[0] if trace_output_digests else None
    trace_observable_id = trace_observable_ids[0] if trace_observable_ids else None
    trace_observable_family = trace_observable_families[0] if trace_observable_families else None
    trace_output_path = trace_output_paths[0] if trace_output_paths else None
    coupled_system_artifact_id = coupled_system_artifact_ids[0] if coupled_system_artifact_ids else None
    coupled_system_digest = coupled_system_digests[0] if coupled_system_digests else None
    coupling_kind = package.get("coupling_kind") or trace.get("coupling_kind") or package.get("coupling_id")
    formulation_id = package.get("formulation_id") or trace.get("formulation_id") or package.get("formulation")
    bem_kernel_family = (
        package.get("bem_kernel_family")
        or trace.get("bem_kernel_family")
        or package.get("kernel_family")
    )
    coupling_convention_schema_id_names = (
        "coupling_convention_schema_id",
        "couplingConventionSchemaId",
        "fem_bem_coupling_convention_schema_id",
        "femBemCouplingConventionSchemaId",
        "physics_convention_schema_id",
        "physicsConventionSchemaId",
        "fem_bem_physics_convention_schema_id",
        "femBemPhysicsConventionSchemaId",
    )
    coupling_convention_schema_id_values = (
        _record_values(package, *coupling_convention_schema_id_names)
        + _record_values(trace, *coupling_convention_schema_id_names)
        + _record_values(operators, *coupling_convention_schema_id_names)
        + _record_values(operator_bem, *coupling_convention_schema_id_names)
    )
    coupling_convention_schema_ids = _unique_strings(coupling_convention_schema_id_values)
    coupling_convention_schema_id = (
        coupling_convention_schema_ids[0] if coupling_convention_schema_ids else None
    )
    postprocess_row_convention_schema_id_names = (
        "fem_bem_postprocess_row_convention_schema_id",
        "femBemPostprocessRowConventionSchemaId",
        "postprocess_row_convention_schema_id",
        "postprocessRowConventionSchemaId",
        "trace_postprocess_row_convention_schema_id",
        "tracePostprocessRowConventionSchemaId",
    )
    postprocess_row_convention_schema_id_values = (
        _record_values(package, *postprocess_row_convention_schema_id_names)
        + _record_values(trace, *postprocess_row_convention_schema_id_names)
        + _record_values(operators, *postprocess_row_convention_schema_id_names)
        + _record_values(operator_trace, *postprocess_row_convention_schema_id_names)
    )
    postprocess_row_convention_schema_ids = _unique_strings(
        postprocess_row_convention_schema_id_values
    )
    postprocess_row_convention_schema_id = (
        postprocess_row_convention_schema_ids[0]
        if postprocess_row_convention_schema_ids
        else None
    )
    trace_basis_schema_id_names = (
        "trace_basis_schema_id",
        "traceBasisSchemaId",
        "fem_bem_trace_basis_schema_id",
        "femBemTraceBasisSchemaId",
        "h1_bem_trace_basis_schema_id",
        "h1BemTraceBasisSchemaId",
        "surface_trace_basis_schema_id",
        "surfaceTraceBasisSchemaId",
        "basis_schema_id",
        "basisSchemaId",
    )
    trace_basis_schema_id_values = (
        _record_values(package, *trace_basis_schema_id_names)
        + _record_values(trace, *trace_basis_schema_id_names)
        + _record_values(operators, *trace_basis_schema_id_names)
        + _record_values(operator_trace, *trace_basis_schema_id_names)
        + _record_values(gypsilab, *trace_basis_schema_id_names)
    )
    trace_basis_schema_ids = _unique_strings(trace_basis_schema_id_values)
    trace_basis_schema_id = trace_basis_schema_ids[0] if trace_basis_schema_ids else None
    assembly_rule_id_names = (
        "assembly_rule_id",
        "assemblyRuleId",
        "operator_assembly_rule_id",
        "operatorAssemblyRuleId",
        "bem_operator_assembly_rule_id",
        "bemOperatorAssemblyRuleId",
        "fem_bem_assembly_rule_id",
        "femBemAssemblyRuleId",
    )
    quadrature_rule_id_names = (
        "quadrature_rule_id",
        "quadratureRuleId",
        "surface_quadrature_rule_id",
        "surfaceQuadratureRuleId",
        "bem_quadrature_rule_id",
        "bemQuadratureRuleId",
        "singular_quadrature_rule_id",
        "singularQuadratureRuleId",
    )
    assembly_rule_id_values = (
        _record_values(package, *assembly_rule_id_names)
        + _record_values(trace, *assembly_rule_id_names)
        + _record_values(operator_trace, *assembly_rule_id_names)
        + _record_values(operator_bem, *assembly_rule_id_names)
        + _record_values(gypsilab, *assembly_rule_id_names)
    )
    quadrature_rule_id_values = (
        _record_values(package, *quadrature_rule_id_names)
        + _record_values(trace, *quadrature_rule_id_names)
        + _record_values(operator_trace, *quadrature_rule_id_names)
        + _record_values(operator_bem, *quadrature_rule_id_names)
        + _record_values(gypsilab, *quadrature_rule_id_names)
    )
    assembly_rule_ids = _unique_strings(assembly_rule_id_values)
    quadrature_rule_ids = _unique_strings(quadrature_rule_id_values)
    assembly_rule_id = assembly_rule_ids[0] if assembly_rule_ids else None
    quadrature_rule_id = quadrature_rule_ids[0] if quadrature_rule_ids else None
    volume_space = package.get("volume_space") or trace.get("volume_space") or package.get("fem_space")
    surface_space = package.get("surface_space") or trace.get("surface_space") or package.get("bem_space")
    boundary_numbers = _int_list(
        trace.get("boundary_numbers")
        or trace.get("boundaryNumbers")
        or trace.get("boundary_ids")
        or gypsilab.get("boundary_numbers")
        or gypsilab.get("boundaryNumbers")
        or gypsilab.get("col")
        or package.get("boundary_numbers")
    )
    boundary_names = _string_list(
        trace.get("boundary_names")
        or trace.get("boundaryNames")
        or gypsilab.get("boundary_names")
        or gypsilab.get("boundaryNames")
        or package.get("boundary_names")
    )
    expected_boundary_numbers_list = _int_list(expected_boundary_numbers)
    expected_boundary_names_list = _string_list(expected_boundary_names)

    conn_matrix = [list(row) for row in _as_list(geo.get("conn_matrix") or geo.get("tets"))]
    gypsilab_elt = [list(row) for row in _as_list(gypsilab.get("elt") or gypsilab.get("triangles"))]
    surface_triangles = [
        list(row) for row in _as_list(trace.get("surface_triangles") or gypsilab_elt)
    ]
    if boundary_numbers and len(boundary_numbers) == 1 and len(surface_triangles) > 1:
        boundary_numbers = boundary_numbers * len(surface_triangles)
    if boundary_names and len(boundary_names) == 1 and len(surface_triangles) > 1:
        boundary_names = boundary_names * len(surface_triangles)
    fem_node_ids = _int_list(trace.get("fem_node_ids") or trace.get("boundary_node_ids"))
    bem_node_ids = _int_list(trace.get("bem_node_ids") or trace.get("boundary_node_ids"))
    trace_row_identity = _trace_row_identity_records(
        trace.get("trace_row_identity")
        or trace.get("trace_rows")
        or package.get("trace_row_identity")
    )
    operator_trace_row_identity = _trace_row_identity_records(
        operator_trace.get("trace_row_identity")
        or operator_trace.get("traceRowIdentity")
        or trace.get("operator_trace_row_identity")
        or package.get("operator_trace_row_identity")
    )
    boundary_row_identity = _boundary_row_identity_records(
        trace.get("boundary_row_identity")
        or trace.get("boundaryRowIdentity")
        or gypsilab.get("boundary_row_identity")
        or gypsilab.get("boundaryRowIdentity")
        or package.get("boundary_row_identity")
        or package.get("boundaryRowIdentity")
    )
    operator_boundary_row_identity = _boundary_row_identity_records(
        operator_trace.get("boundary_row_identity")
        or operator_trace.get("boundaryRowIdentity")
        or trace.get("operator_boundary_row_identity")
        or package.get("operator_boundary_row_identity")
    )
    expected_boundary_row_identity_normalized = _normalize_boundary_row_identity(
        _boundary_row_identity_records(expected_boundary_row_identity)
    )
    trace_matrix = trace.get("trace_matrix") or trace.get("matrix")
    if trace_matrix is None:
        trace_rows = []
    else:
        trace_rows = [[float(value) for value in row] for row in trace_matrix]

    n_fem = int(geo.get("N") or geo.get("n_nodes") or len(geo.get("nodes") or []))
    if n_fem <= 0:
        max_conn_node = max((int(node) for row in conn_matrix for node in row), default=0)
        max_trace_node = max(fem_node_ids, default=0)
        n_fem = max(max_conn_node, max_trace_node)

    tolerance = float(tol)
    if tolerance < 0.0:
        raise ValueError("tol must be non-negative")

    boundary_nodes_from_triangles = sorted({int(node) for tri in surface_triangles for node in tri})
    trace_nodes_sorted = sorted(fem_node_ids)
    one_based_ids_in_range = all(1 <= int(node) <= n_fem for node in fem_node_ids + bem_node_ids)
    trace_matrix_shape_ok = (
        len(trace_rows) == len(fem_node_ids)
        and all(len(row) == n_fem for row in trace_rows)
    )
    trace_matrix_one_hot = False
    trace_matrix_matches_ids = False
    trace_matrix_nonzero_node_ids = []
    if trace_matrix_shape_ok:
        row_one_hot = []
        row_matches = []
        for row, node_id in zip(trace_rows, fem_node_ids):
            nonzero = [i for i, value in enumerate(row, start=1) if abs(value) > tolerance]
            trace_matrix_nonzero_node_ids.append(nonzero[0] if len(nonzero) == 1 else None)
            row_one_hot.append(len(nonzero) == 1 and abs(row[nonzero[0] - 1] - 1.0) <= tolerance)
            row_matches.append(nonzero == [int(node_id)])
        trace_matrix_one_hot = all(row_one_hot)
        trace_matrix_matches_ids = all(row_matches)

    trace_row_identity_present = bool(trace_row_identity)
    trace_row_identity_normalized = [
        {
            "trace_row_index": _maybe_int(
                _record_value(record, "trace_row_index", "row_index", "row")
            ),
            "fem_node_id": _maybe_int(_record_value(record, "fem_node_id", "fem_node")),
            "bem_node_id": _maybe_int(_record_value(record, "bem_node_id", "bem_node")),
            "surface_node_index": _maybe_int(
                _record_value(record, "surface_node_index", "surface_node_id")
            ),
        }
        for record in trace_row_identity
    ]
    trace_row_identity_row_indices = [
        record["trace_row_index"] for record in trace_row_identity_normalized
    ]
    trace_row_identity_fem_node_ids = [
        record["fem_node_id"] for record in trace_row_identity_normalized
    ]
    trace_row_identity_bem_node_ids = [
        record["bem_node_id"] for record in trace_row_identity_normalized
    ]
    trace_row_identity_mismatch_rows = []
    if trace_row_identity_present:
        expected_row_indices = list(range(1, len(fem_node_ids) + 1))
        for index, record in enumerate(trace_row_identity_normalized, start=1):
            row_ok = (
                record["trace_row_index"] == index
                and index <= len(fem_node_ids)
                and record["fem_node_id"] == fem_node_ids[index - 1]
                and index <= len(bem_node_ids)
                and record["bem_node_id"] == bem_node_ids[index - 1]
            )
            if trace_matrix_shape_ok and index <= len(trace_matrix_nonzero_node_ids):
                row_ok = row_ok and record["fem_node_id"] == trace_matrix_nonzero_node_ids[index - 1]
            if not row_ok:
                trace_row_identity_mismatch_rows.append(index)
        for index in range(len(trace_row_identity_normalized) + 1, len(expected_row_indices) + 1):
            trace_row_identity_mismatch_rows.append(index)
    operator_trace_row_identity_present = bool(operator_trace_row_identity)
    operator_trace_row_identity_normalized = [
        {
            "trace_row_index": _maybe_int(
                _record_value(record, "trace_row_index", "row_index", "row")
            ),
            "fem_node_id": _maybe_int(_record_value(record, "fem_node_id", "fem_node")),
            "bem_node_id": _maybe_int(_record_value(record, "bem_node_id", "bem_node")),
            "surface_node_index": _maybe_int(
                _record_value(record, "surface_node_index", "surface_node_id")
            ),
        }
        for record in operator_trace_row_identity
    ]
    operator_trace_row_identity_row_indices = [
        record["trace_row_index"] for record in operator_trace_row_identity_normalized
    ]
    operator_trace_row_identity_fem_node_ids = [
        record["fem_node_id"] for record in operator_trace_row_identity_normalized
    ]
    operator_trace_row_identity_bem_node_ids = [
        record["bem_node_id"] for record in operator_trace_row_identity_normalized
    ]
    operator_trace_row_identity_mismatch_rows = []
    if operator_trace_row_identity_present:
        expected_row_indices = list(range(1, len(fem_node_ids) + 1))
        for index, record in enumerate(operator_trace_row_identity_normalized, start=1):
            row_ok = (
                record["trace_row_index"] == index
                and index <= len(fem_node_ids)
                and record["fem_node_id"] == fem_node_ids[index - 1]
                and index <= len(bem_node_ids)
                and record["bem_node_id"] == bem_node_ids[index - 1]
            )
            if trace_matrix_shape_ok and index <= len(trace_matrix_nonzero_node_ids):
                row_ok = row_ok and record["fem_node_id"] == trace_matrix_nonzero_node_ids[index - 1]
            if trace_row_identity_present and index <= len(trace_row_identity_normalized):
                row_ok = row_ok and record == trace_row_identity_normalized[index - 1]
            if not row_ok:
                operator_trace_row_identity_mismatch_rows.append(index)
        for index in range(len(operator_trace_row_identity_normalized) + 1, len(expected_row_indices) + 1):
            operator_trace_row_identity_mismatch_rows.append(index)

    boundary_row_identity_present = bool(boundary_row_identity)
    boundary_row_identity_normalized = _normalize_boundary_row_identity(boundary_row_identity)
    operator_boundary_row_identity_present = bool(operator_boundary_row_identity)
    operator_boundary_row_identity_normalized = _normalize_boundary_row_identity(operator_boundary_row_identity)
    boundary_row_identity_mismatch_rows = []
    if boundary_row_identity_present:
        for index, record in enumerate(boundary_row_identity_normalized, start=1):
            row_ok = record["surface_triangle_index"] == index
            if index <= len(surface_triangles):
                row_ok = row_ok and record["surface_triangle_nodes"] == [int(node) for node in surface_triangles[index - 1]]
            if index <= len(boundary_numbers):
                row_ok = row_ok and record["boundary_number"] == int(boundary_numbers[index - 1])
            if index <= len(boundary_names):
                row_ok = row_ok and record["boundary_name"] == str(boundary_names[index - 1])
            if not row_ok:
                boundary_row_identity_mismatch_rows.append(index)
        for index in range(len(boundary_row_identity_normalized) + 1, len(surface_triangles) + 1):
            boundary_row_identity_mismatch_rows.append(index)
    operator_boundary_row_identity_mismatch_rows = []
    if operator_boundary_row_identity_present:
        for index, record in enumerate(operator_boundary_row_identity_normalized, start=1):
            row_ok = record["surface_triangle_index"] == index
            if index <= len(surface_triangles):
                row_ok = row_ok and record["surface_triangle_nodes"] == [int(node) for node in surface_triangles[index - 1]]
            if index <= len(boundary_numbers):
                row_ok = row_ok and record["boundary_number"] == int(boundary_numbers[index - 1])
            if index <= len(boundary_names):
                row_ok = row_ok and record["boundary_name"] == str(boundary_names[index - 1])
            if boundary_row_identity_present and index <= len(boundary_row_identity_normalized):
                row_ok = row_ok and record == boundary_row_identity_normalized[index - 1]
            if not row_ok:
                operator_boundary_row_identity_mismatch_rows.append(index)
        for index in range(len(operator_boundary_row_identity_normalized) + 1, len(surface_triangles) + 1):
            operator_boundary_row_identity_mismatch_rows.append(index)

    checks = {
        "mesh_id_recorded": bool(mesh_id),
        "export_id_recorded": bool(export_id),
        "trace_artifact_id_recorded": bool(trace_artifact_id),
        "surface_mesh_id_recorded": bool(surface_mesh_id),
        "trace_operator_artifact_id_recorded": bool(trace_operator_artifact_id),
        "trace_operator_policy_recorded": bool(trace_operator_policy),
        "trace_output_artifact_id_consistent_when_present": len(trace_output_artifact_ids) <= 1,
        "trace_output_digest_consistent_when_present": len(trace_output_digests) <= 1,
        "trace_observable_id_consistent_when_present": len(trace_observable_ids) <= 1,
        "trace_observable_family_consistent_when_present": len(trace_observable_families) <= 1,
        "trace_output_path_consistent_when_present": len(trace_output_paths) <= 1,
        "coupled_system_artifact_id_consistent_when_present": len(coupled_system_artifact_ids) <= 1,
        "coupled_system_digest_consistent_when_present": len(coupled_system_digests) <= 1,
        "linear_solver_report_artifact_id_consistent_when_present": (
            len(linear_solver_report_artifact_ids) <= 1
        ),
        "linear_solver_report_digest_consistent_when_present": (
            len(linear_solver_report_digests) <= 1
        ),
        "linear_solver_name_consistent_when_present": len(linear_solver_names) <= 1,
        "linear_solver_tolerance_consistent_when_present": len(linear_solver_tolerances) <= 1,
        "linear_solver_residual_norm_consistent_when_present": (
            len(linear_solver_residual_norms) <= 1
        ),
        "linear_solver_iteration_count_consistent_when_present": (
            len(linear_solver_iteration_counts) <= 1
        ),
        "parameter_set_artifact_id_consistent_when_present": (
            len(parameter_set_artifact_ids) <= 1
        ),
        "parameter_set_digest_consistent_when_present": len(parameter_set_digests) <= 1,
        "parameter_set_path_consistent_when_present": len(parameter_set_paths) <= 1,
        "objective_observable_id_consistent_when_present": (
            len(objective_observable_ids) <= 1
        ),
        "objective_observable_family_consistent_when_present": (
            len(objective_observable_families) <= 1
        ),
        "linear_solver_tolerance_finite_positive_when_present": (
            linear_solver_tolerance is None or linear_solver_tolerance > 0.0
        ),
        "linear_solver_residual_norm_finite_nonnegative_when_present": (
            linear_solver_residual_norm is None or linear_solver_residual_norm >= 0.0
        ),
        "linear_solver_iteration_count_nonnegative_when_present": (
            linear_solver_iteration_count is None or linear_solver_iteration_count >= 0
        ),
        "assembly_rule_id_consistent_when_present": len(assembly_rule_ids) <= 1,
        "quadrature_rule_id_consistent_when_present": len(quadrature_rule_ids) <= 1,
        "coupling_convention_schema_id_consistent_when_present": (
            len(coupling_convention_schema_ids) <= 1
        ),
        "coupling_convention_schema_id_recorded_when_required": (
            not require_coupling_convention_schema or bool(coupling_convention_schema_id)
        ),
        "fem_bem_postprocess_row_convention_schema_id_consistent_when_present": (
            len(postprocess_row_convention_schema_ids) <= 1
        ),
        "fem_bem_postprocess_row_convention_schema_id_recorded_when_required": (
            not require_fem_bem_postprocess_row_convention_schema
            or bool(postprocess_row_convention_schema_id)
        ),
        "trace_basis_schema_id_consistent_when_present": (
            len(trace_basis_schema_ids) <= 1
        ),
        "trace_basis_schema_id_recorded_when_required": (
            not require_trace_basis_schema or bool(trace_basis_schema_id)
        ),
        "trace_output_artifact_id_recorded_when_required": (
            not require_trace_output_artifact or bool(trace_output_artifact_id)
        ),
        "trace_output_digest_recorded_when_required": (
            not require_trace_output_artifact or bool(trace_output_digest)
        ),
        "trace_output_path_recorded_when_required": (
            not require_trace_output_artifact or bool(trace_output_path)
        ),
        "linear_solver_report_artifact_id_recorded_when_required": (
            not require_linear_solver_report or bool(linear_solver_report_artifact_id)
        ),
        "linear_solver_report_digest_recorded_when_required": (
            not require_linear_solver_report or bool(linear_solver_report_digest)
        ),
        "linear_solver_name_recorded_when_required": (
            not require_linear_solver_report or bool(linear_solver_name)
        ),
        "linear_solver_tolerance_recorded_when_required": (
            not require_linear_solver_report or linear_solver_tolerance is not None
        ),
        "linear_solver_residual_norm_recorded_when_required": (
            not require_linear_solver_report or linear_solver_residual_norm is not None
        ),
        "parameter_set_artifact_id_recorded_when_required": (
            not require_parameter_set_artifact or bool(parameter_set_artifact_id)
        ),
        "parameter_set_digest_recorded_when_required": (
            not require_parameter_set_artifact or bool(parameter_set_digest)
        ),
        "parameter_set_path_recorded_when_required": (
            not require_parameter_set_artifact or bool(parameter_set_path)
        ),
        "source_path_recorded": bool(source_path),
        "source_file_id_consistent_when_present": (
            not source_file_id
            or not trace_source_file_id
            or str(source_file_id) == str(trace_source_file_id)
        ),
        "source_format_recorded": bool(source_format),
        "source_format_is_vol": _norm(source_format) in {"vol", ".vol", "netgen_vol"},
        "policy_recorded": _norm(policy) == "netgen_vol_tri_tet_only_shared_one_based_nodes",
        "polynomial_order_first_order": polynomial_order == 1,
        "curvedelements_absent": curved_element_count == 0,
        "volume_tetrahedra_first_order": bool(conn_matrix) and all(len(row) == 4 for row in conn_matrix),
        "boundary_triangles_first_order": bool(surface_triangles) and all(len(row) == 3 for row in surface_triangles),
        "gypsilab_triangles_match_trace_triangles": gypsilab_elt == surface_triangles,
        "trace_node_ids_recorded": bool(fem_node_ids) and bool(bem_node_ids),
        "fem_bem_node_ids_identical": fem_node_ids == bem_node_ids,
        "trace_node_ids_sorted_unique": fem_node_ids == sorted(set(fem_node_ids)),
        "trace_covers_boundary_triangle_nodes": trace_nodes_sorted == boundary_nodes_from_triangles,
        "one_based_node_ids_in_range": one_based_ids_in_range,
        "trace_matrix_shape_ok": trace_matrix_shape_ok,
        "trace_matrix_one_hot": trace_matrix_one_hot,
        "trace_matrix_matches_fem_node_ids": trace_matrix_matches_ids,
    }
    if trace_row_identity_present:
        checks["trace_row_identity_rows_match_trace_rows"] = len(trace_row_identity_normalized) == len(fem_node_ids)
        checks["trace_row_identity_row_indices_match"] = trace_row_identity_row_indices == list(range(1, len(fem_node_ids) + 1))
        checks["trace_row_identity_fem_nodes_match"] = trace_row_identity_fem_node_ids == fem_node_ids
        checks["trace_row_identity_bem_nodes_match"] = trace_row_identity_bem_node_ids == bem_node_ids
        checks["trace_row_identity_matches_trace_matrix"] = (
            trace_matrix_shape_ok
            and trace_row_identity_fem_node_ids == trace_matrix_nonzero_node_ids
        )
    if operator_trace_row_identity_present:
        checks["operator_trace_row_identity_rows_match_trace_rows"] = (
            len(operator_trace_row_identity_normalized) == len(fem_node_ids)
        )
        checks["operator_trace_row_identity_row_indices_match"] = (
            operator_trace_row_identity_row_indices == list(range(1, len(fem_node_ids) + 1))
        )
        checks["operator_trace_row_identity_fem_nodes_match"] = operator_trace_row_identity_fem_node_ids == fem_node_ids
        checks["operator_trace_row_identity_bem_nodes_match"] = operator_trace_row_identity_bem_node_ids == bem_node_ids
        checks["operator_trace_row_identity_matches_trace_matrix"] = (
            trace_matrix_shape_ok
            and operator_trace_row_identity_fem_node_ids == trace_matrix_nonzero_node_ids
        )
        checks["operator_trace_row_identity_matches_trace_identity"] = (
            not trace_row_identity_present
            or operator_trace_row_identity_normalized == trace_row_identity_normalized
        )
    if boundary_numbers:
        checks["boundary_numbers_match_surface_triangles"] = len(boundary_numbers) == len(surface_triangles)
        checks["boundary_numbers_positive"] = all(number > 0 for number in boundary_numbers)
    if boundary_names:
        checks["boundary_names_match_surface_triangles"] = len(boundary_names) == len(surface_triangles)
        checks["boundary_names_recorded"] = all(bool(name) for name in boundary_names)
    if boundary_numbers and boundary_names and len(boundary_numbers) == len(boundary_names):
        number_to_names = {}
        for number, name in zip(boundary_numbers, boundary_names):
            number_to_names.setdefault(number, set()).add(name)
        checks["boundary_number_name_pairs_consistent"] = all(
            len(names) == 1 for names in number_to_names.values()
        )
    if boundary_row_identity_present:
        checks["boundary_row_identity_rows_match_surface_triangles"] = (
            len(boundary_row_identity_normalized) == len(surface_triangles)
        )
        checks["boundary_row_identity_matches_surface_triangles"] = not boundary_row_identity_mismatch_rows
        checks["boundary_row_identity_boundary_numbers_match"] = (
            not boundary_numbers
            or [record["boundary_number"] for record in boundary_row_identity_normalized] == boundary_numbers
        )
        checks["boundary_row_identity_boundary_names_match"] = (
            not boundary_names
            or [record["boundary_name"] for record in boundary_row_identity_normalized] == boundary_names
        )
    if operator_boundary_row_identity_present:
        checks["operator_boundary_row_identity_rows_match_surface_triangles"] = (
            len(operator_boundary_row_identity_normalized) == len(surface_triangles)
        )
        checks["operator_boundary_row_identity_matches_surface_triangles"] = (
            not operator_boundary_row_identity_mismatch_rows
        )
        checks["operator_boundary_row_identity_matches_trace_identity"] = (
            not boundary_row_identity_present
            or operator_boundary_row_identity_normalized == boundary_row_identity_normalized
        )
    if expected_mesh_id is not None:
        checks["expected_mesh_id_matches"] = str(mesh_id) == str(expected_mesh_id)
    if expected_export_id is not None:
        checks["expected_export_id_matches"] = str(export_id) == str(expected_export_id)
    if expected_trace_artifact_id is not None:
        checks["expected_trace_artifact_id_matches"] = str(trace_artifact_id) == str(expected_trace_artifact_id)
    if expected_surface_mesh_id is not None:
        checks["expected_surface_mesh_id_matches"] = str(surface_mesh_id) == str(expected_surface_mesh_id)
    if expected_source_file_id is not None:
        expected = str(expected_source_file_id)
        checks["source_file_id_recorded_when_expected"] = bool(source_file_id)
        checks["expected_source_file_id_matches"] = str(source_file_id) == expected
        if trace_source_file_id:
            checks["expected_trace_source_file_id_matches"] = str(trace_source_file_id) == expected
    if expected_source_format is not None:
        checks["expected_source_format_matches"] = _norm(source_format) == _norm(expected_source_format)
    if expected_coupling_kind is not None:
        checks["coupling_kind_recorded_when_expected"] = bool(coupling_kind)
        checks["expected_coupling_kind_matches"] = _norm(coupling_kind) == _norm(expected_coupling_kind)
    if expected_formulation_id is not None:
        checks["formulation_id_recorded_when_expected"] = bool(formulation_id)
        checks["expected_formulation_id_matches"] = _norm(formulation_id) == _norm(expected_formulation_id)
    if expected_bem_kernel_family is not None:
        checks["bem_kernel_family_recorded_when_expected"] = bool(bem_kernel_family)
        checks["expected_bem_kernel_family_matches"] = _norm(bem_kernel_family) == _norm(expected_bem_kernel_family)
    if expected_coupling_convention_schema_id is not None:
        expected = str(expected_coupling_convention_schema_id)
        checks["coupling_convention_schema_id_recorded_when_expected"] = bool(
            coupling_convention_schema_id
        )
        checks["expected_coupling_convention_schema_id_matches"] = (
            coupling_convention_schema_ids == [expected]
        )
    if expected_fem_bem_postprocess_row_convention_schema_id is not None:
        expected = str(expected_fem_bem_postprocess_row_convention_schema_id)
        checks["fem_bem_postprocess_row_convention_schema_id_recorded_when_expected"] = bool(
            postprocess_row_convention_schema_id
        )
        checks["expected_fem_bem_postprocess_row_convention_schema_id_matches"] = (
            postprocess_row_convention_schema_ids == [expected]
        )
    if expected_trace_basis_schema_id is not None:
        expected = str(expected_trace_basis_schema_id)
        checks["trace_basis_schema_id_recorded_when_expected"] = bool(
            trace_basis_schema_id
        )
        checks["expected_trace_basis_schema_id_matches"] = (
            trace_basis_schema_ids == [expected]
        )
    if expected_assembly_rule_id is not None:
        checks["assembly_rule_id_recorded_when_expected"] = bool(assembly_rule_id)
        checks["expected_assembly_rule_id_matches"] = str(assembly_rule_id) == str(expected_assembly_rule_id)
    if expected_quadrature_rule_id is not None:
        checks["quadrature_rule_id_recorded_when_expected"] = bool(quadrature_rule_id)
        checks["expected_quadrature_rule_id_matches"] = str(quadrature_rule_id) == str(expected_quadrature_rule_id)
    if expected_volume_space is not None:
        checks["volume_space_recorded_when_expected"] = bool(volume_space)
        checks["expected_volume_space_matches"] = _norm(volume_space) == _norm(expected_volume_space)
    if expected_surface_space is not None:
        checks["surface_space_recorded_when_expected"] = bool(surface_space)
        checks["expected_surface_space_matches"] = _norm(surface_space) == _norm(expected_surface_space)
    if expected_boundary_numbers is not None:
        checks["boundary_numbers_recorded_when_expected"] = bool(boundary_numbers)
        checks["expected_boundary_numbers_match"] = sorted(set(boundary_numbers)) == sorted(set(expected_boundary_numbers_list))
    if expected_boundary_names is not None:
        checks["boundary_names_recorded_when_expected"] = bool(boundary_names)
        checks["expected_boundary_names_match"] = sorted(set(boundary_names)) == sorted(set(expected_boundary_names_list))
    if expected_boundary_row_identity is not None:
        checks["boundary_row_identity_recorded_when_expected"] = bool(boundary_row_identity_normalized)
        checks["expected_boundary_row_identity_matches"] = (
            boundary_row_identity_normalized == expected_boundary_row_identity_normalized
        )
    if expected_trace_operator_artifact_id is not None:
        checks["expected_trace_operator_artifact_id_matches"] = (
            str(trace_operator_artifact_id) == str(expected_trace_operator_artifact_id)
        )
    if expected_trace_operator_policy is not None:
        checks["expected_trace_operator_policy_matches"] = (
            _norm(trace_operator_policy) == _norm(expected_trace_operator_policy)
        )
    if expected_trace_output_artifact_id is not None:
        checks["trace_output_artifact_id_recorded_when_expected"] = bool(trace_output_artifact_id)
        checks["expected_trace_output_artifact_id_matches"] = (
            str(trace_output_artifact_id) == str(expected_trace_output_artifact_id)
        )
    if expected_trace_output_digest is not None:
        checks["trace_output_digest_recorded_when_expected"] = bool(trace_output_digest)
        checks["expected_trace_output_digest_matches"] = (
            str(trace_output_digest) == str(expected_trace_output_digest)
        )
    if expected_trace_observable_id is not None:
        checks["trace_observable_id_recorded_when_expected"] = bool(trace_observable_id)
        checks["expected_trace_observable_id_matches"] = (
            trace_observable_ids == [str(expected_trace_observable_id)]
        )
    if expected_trace_observable_family is not None:
        checks["trace_observable_family_recorded_when_expected"] = bool(trace_observable_family)
        checks["expected_trace_observable_family_matches"] = (
            trace_observable_families == [_norm(expected_trace_observable_family)]
        )
    if expected_coupled_system_artifact_id is not None:
        checks["coupled_system_artifact_id_recorded_when_expected"] = bool(coupled_system_artifact_id)
        checks["expected_coupled_system_artifact_id_matches"] = (
            str(coupled_system_artifact_id) == str(expected_coupled_system_artifact_id)
        )
    if expected_coupled_system_digest is not None:
        checks["coupled_system_digest_recorded_when_expected"] = bool(coupled_system_digest)
        checks["expected_coupled_system_digest_matches"] = (
            str(coupled_system_digest) == str(expected_coupled_system_digest)
        )
    if expected_linear_solver_report_artifact_id is not None:
        checks["linear_solver_report_artifact_id_recorded_when_expected"] = bool(
            linear_solver_report_artifact_id
        )
        checks["expected_linear_solver_report_artifact_id_matches"] = (
            str(linear_solver_report_artifact_id)
            == str(expected_linear_solver_report_artifact_id)
        )
    if expected_linear_solver_report_digest is not None:
        checks["linear_solver_report_digest_recorded_when_expected"] = bool(
            linear_solver_report_digest
        )
        checks["expected_linear_solver_report_digest_matches"] = (
            str(linear_solver_report_digest) == str(expected_linear_solver_report_digest)
        )
    if expected_linear_solver_name is not None:
        checks["linear_solver_name_recorded_when_expected"] = bool(linear_solver_name)
        checks["expected_linear_solver_name_matches"] = (
            linear_solver_names == [_norm(expected_linear_solver_name)]
        )
    if expected_linear_solver_tolerance is not None:
        expected_tolerance = float(expected_linear_solver_tolerance)
        checks["linear_solver_tolerance_recorded_when_expected"] = (
            linear_solver_tolerance is not None
        )
        checks["expected_linear_solver_tolerance_matches"] = (
            linear_solver_tolerance is not None
            and abs(linear_solver_tolerance - expected_tolerance)
            <= max(tolerance, abs(expected_tolerance) * 1.0e-12)
        )
    if expected_linear_solver_residual_norm_max is not None:
        max_residual_norm = float(expected_linear_solver_residual_norm_max)
        checks["linear_solver_residual_norm_recorded_when_expected"] = (
            linear_solver_residual_norm is not None
        )
        checks["linear_solver_residual_norm_below_expected_max"] = (
            linear_solver_residual_norm is not None
            and linear_solver_residual_norm <= max_residual_norm
        )
    if expected_parameter_set_artifact_id is not None:
        checks["parameter_set_artifact_id_recorded_when_expected"] = bool(
            parameter_set_artifact_id
        )
        checks["expected_parameter_set_artifact_id_matches"] = (
            str(parameter_set_artifact_id) == str(expected_parameter_set_artifact_id)
        )
    if expected_parameter_set_digest is not None:
        checks["parameter_set_digest_recorded_when_expected"] = bool(parameter_set_digest)
        checks["expected_parameter_set_digest_matches"] = (
            str(parameter_set_digest) == str(expected_parameter_set_digest)
        )
    if expected_parameter_set_path is not None:
        checks["parameter_set_path_recorded_when_expected"] = bool(parameter_set_path)
        checks["expected_parameter_set_path_matches"] = (
            str(parameter_set_path) == str(expected_parameter_set_path)
        )
    if expected_objective_observable_id is not None:
        checks["objective_observable_id_recorded_when_expected"] = bool(
            objective_observable_id
        )
        checks["expected_objective_observable_id_matches"] = (
            objective_observable_ids == [str(expected_objective_observable_id)]
        )
    if expected_objective_observable_family is not None:
        checks["objective_observable_family_recorded_when_expected"] = bool(
            objective_observable_family
        )
        checks["expected_objective_observable_family_matches"] = (
            objective_observable_families == [_norm(expected_objective_observable_family)]
        )
    if require_result_provenance or expected_result_artifact_id is not None:
        checks["result_artifact_id_recorded_when_required"] = bool(result_artifact_id)
    if expected_result_artifact_id is not None:
        checks["expected_result_artifact_id_matches"] = (
            str(result_artifact_id) == str(expected_result_artifact_id)
        )
    if require_result_provenance:
        checks["run_started_at_recorded_when_required"] = bool(run_started_at)
        checks["matlab_version_recorded_when_required"] = bool(matlab_version)
        checks["timing_breakdown_recorded_when_required"] = bool(timing_breakdown_seconds)
        checks["timing_breakdown_has_at_least_four_items"] = len(timing_breakdown_seconds) >= 4
        checks["timing_breakdown_values_finite_nonnegative"] = all(
            math.isfinite(value) and value >= 0.0
            for value in timing_breakdown_seconds.values()
        )
    if expected_matlab_version is not None:
        checks["expected_matlab_version_matches"] = str(matlab_version) == str(expected_matlab_version)

    return {
        "policy": "netgen_vol_first_order_fem_bem_trace_package_handoff",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "mesh_id": None if mesh_id is None else str(mesh_id),
        "export_id": None if export_id is None else str(export_id),
        "trace_artifact_id": None if trace_artifact_id is None else str(trace_artifact_id),
        "trace_operator_artifact_id": None if trace_operator_artifact_id is None else str(trace_operator_artifact_id),
        "trace_operator_policy": None if trace_operator_policy is None else str(trace_operator_policy),
        "trace_output_artifact_id": None if trace_output_artifact_id is None else str(trace_output_artifact_id),
        "trace_output_digest": None if trace_output_digest is None else str(trace_output_digest),
        "trace_observable_id": None if trace_observable_id is None else str(trace_observable_id),
        "trace_observable_family": None if trace_observable_family is None else str(trace_observable_family),
        "trace_output_path": None if trace_output_path is None else str(trace_output_path),
        "coupled_system_artifact_id": None if coupled_system_artifact_id is None else str(coupled_system_artifact_id),
        "coupled_system_digest": None if coupled_system_digest is None else str(coupled_system_digest),
        "trace_output_artifact_ids": trace_output_artifact_ids,
        "trace_output_digests": trace_output_digests,
        "trace_observable_ids": trace_observable_ids,
        "trace_observable_families": trace_observable_families,
        "trace_output_paths": trace_output_paths,
        "coupled_system_artifact_ids": coupled_system_artifact_ids,
        "coupled_system_digests": coupled_system_digests,
        "linear_solver_report_artifact_id": None
        if linear_solver_report_artifact_id is None
        else str(linear_solver_report_artifact_id),
        "linear_solver_report_digest": None
        if linear_solver_report_digest is None
        else str(linear_solver_report_digest),
        "linear_solver_name": None if linear_solver_name is None else str(linear_solver_name),
        "linear_solver_tolerance": linear_solver_tolerance,
        "linear_solver_residual_norm": linear_solver_residual_norm,
        "linear_solver_iteration_count": linear_solver_iteration_count,
        "linear_solver_report_artifact_ids": linear_solver_report_artifact_ids,
        "linear_solver_report_digests": linear_solver_report_digests,
        "linear_solver_names": linear_solver_names,
        "linear_solver_tolerances": linear_solver_tolerances,
        "linear_solver_residual_norms": linear_solver_residual_norms,
        "linear_solver_iteration_counts": linear_solver_iteration_counts,
        "parameter_set_artifact_id": None
        if parameter_set_artifact_id is None
        else str(parameter_set_artifact_id),
        "parameter_set_digest": None if parameter_set_digest is None else str(parameter_set_digest),
        "parameter_set_path": None if parameter_set_path is None else str(parameter_set_path),
        "objective_observable_id": None
        if objective_observable_id is None
        else str(objective_observable_id),
        "objective_observable_family": None
        if objective_observable_family is None
        else str(objective_observable_family),
        "parameter_set_artifact_ids": parameter_set_artifact_ids,
        "parameter_set_digests": parameter_set_digests,
        "parameter_set_paths": parameter_set_paths,
        "objective_observable_ids": objective_observable_ids,
        "objective_observable_families": objective_observable_families,
        "result_artifact_id": None if result_artifact_id is None else str(result_artifact_id),
        "run_started_at": None if run_started_at is None else str(run_started_at),
        "matlab_version": None if matlab_version is None else str(matlab_version),
        "timing_breakdown_seconds": timing_breakdown_seconds,
        "timing_breakdown_names": timing_breakdown_names,
        "surface_mesh_id": None if surface_mesh_id is None else str(surface_mesh_id),
        "source_file_id": None if source_file_id is None else str(source_file_id),
        "trace_source_file_id": None if trace_source_file_id is None else str(trace_source_file_id),
        "coupling_kind": None if coupling_kind is None else str(coupling_kind),
        "formulation_id": None if formulation_id is None else str(formulation_id),
        "bem_kernel_family": None if bem_kernel_family is None else str(bem_kernel_family),
        "coupling_convention_schema_id": None
        if coupling_convention_schema_id is None
        else str(coupling_convention_schema_id),
        "coupling_convention_schema_ids": coupling_convention_schema_ids,
        "fem_bem_postprocess_row_convention_schema_id": None
        if postprocess_row_convention_schema_id is None
        else str(postprocess_row_convention_schema_id),
        "fem_bem_postprocess_row_convention_schema_ids": (
            postprocess_row_convention_schema_ids
        ),
        "trace_basis_schema_id": None
        if trace_basis_schema_id is None
        else str(trace_basis_schema_id),
        "trace_basis_schema_ids": trace_basis_schema_ids,
        "assembly_rule_id": None if assembly_rule_id is None else str(assembly_rule_id),
        "quadrature_rule_id": None if quadrature_rule_id is None else str(quadrature_rule_id),
        "assembly_rule_ids": assembly_rule_ids,
        "quadrature_rule_ids": quadrature_rule_ids,
        "volume_space": None if volume_space is None else str(volume_space),
        "surface_space": None if surface_space is None else str(surface_space),
        "expected_coupling_kind": None if expected_coupling_kind is None else str(expected_coupling_kind),
        "expected_formulation_id": None if expected_formulation_id is None else str(expected_formulation_id),
        "expected_bem_kernel_family": None if expected_bem_kernel_family is None else str(expected_bem_kernel_family),
        "expected_coupling_convention_schema_id": None
        if expected_coupling_convention_schema_id is None
        else str(expected_coupling_convention_schema_id),
        "expected_fem_bem_postprocess_row_convention_schema_id": None
        if expected_fem_bem_postprocess_row_convention_schema_id is None
        else str(expected_fem_bem_postprocess_row_convention_schema_id),
        "expected_trace_basis_schema_id": None
        if expected_trace_basis_schema_id is None
        else str(expected_trace_basis_schema_id),
        "expected_assembly_rule_id": None if expected_assembly_rule_id is None else str(expected_assembly_rule_id),
        "expected_quadrature_rule_id": None if expected_quadrature_rule_id is None else str(expected_quadrature_rule_id),
        "expected_volume_space": None if expected_volume_space is None else str(expected_volume_space),
        "expected_surface_space": None if expected_surface_space is None else str(expected_surface_space),
        "expected_source_file_id": None if expected_source_file_id is None else str(expected_source_file_id),
        "expected_boundary_numbers": expected_boundary_numbers_list,
        "expected_boundary_names": expected_boundary_names_list,
        "expected_boundary_row_identity": expected_boundary_row_identity_normalized,
        "expected_trace_operator_artifact_id": None
        if expected_trace_operator_artifact_id is None
        else str(expected_trace_operator_artifact_id),
        "expected_trace_operator_policy": None
        if expected_trace_operator_policy is None
        else str(expected_trace_operator_policy),
        "expected_trace_output_artifact_id": None
        if expected_trace_output_artifact_id is None
        else str(expected_trace_output_artifact_id),
        "expected_trace_output_digest": None
        if expected_trace_output_digest is None
        else str(expected_trace_output_digest),
        "expected_trace_observable_id": None
        if expected_trace_observable_id is None
        else str(expected_trace_observable_id),
        "expected_trace_observable_family": None
        if expected_trace_observable_family is None
        else _norm(expected_trace_observable_family),
        "expected_coupled_system_artifact_id": None
        if expected_coupled_system_artifact_id is None
        else str(expected_coupled_system_artifact_id),
        "expected_coupled_system_digest": None
        if expected_coupled_system_digest is None
        else str(expected_coupled_system_digest),
        "expected_linear_solver_report_artifact_id": None
        if expected_linear_solver_report_artifact_id is None
        else str(expected_linear_solver_report_artifact_id),
        "expected_linear_solver_report_digest": None
        if expected_linear_solver_report_digest is None
        else str(expected_linear_solver_report_digest),
        "expected_linear_solver_name": None
        if expected_linear_solver_name is None
        else _norm(expected_linear_solver_name),
        "expected_linear_solver_tolerance": None
        if expected_linear_solver_tolerance is None
        else float(expected_linear_solver_tolerance),
        "expected_linear_solver_residual_norm_max": None
        if expected_linear_solver_residual_norm_max is None
        else float(expected_linear_solver_residual_norm_max),
        "expected_parameter_set_artifact_id": None
        if expected_parameter_set_artifact_id is None
        else str(expected_parameter_set_artifact_id),
        "expected_parameter_set_digest": None
        if expected_parameter_set_digest is None
        else str(expected_parameter_set_digest),
        "expected_parameter_set_path": None
        if expected_parameter_set_path is None
        else str(expected_parameter_set_path),
        "expected_objective_observable_id": None
        if expected_objective_observable_id is None
        else str(expected_objective_observable_id),
        "expected_objective_observable_family": None
        if expected_objective_observable_family is None
        else _norm(expected_objective_observable_family),
        "expected_result_artifact_id": None
        if expected_result_artifact_id is None
        else str(expected_result_artifact_id),
        "expected_matlab_version": None
        if expected_matlab_version is None
        else str(expected_matlab_version),
        "require_result_provenance": bool(require_result_provenance),
        "require_trace_output_artifact": bool(require_trace_output_artifact),
        "require_linear_solver_report": bool(require_linear_solver_report),
        "require_parameter_set_artifact": bool(require_parameter_set_artifact),
        "require_coupling_convention_schema": bool(require_coupling_convention_schema),
        "require_fem_bem_postprocess_row_convention_schema": bool(
            require_fem_bem_postprocess_row_convention_schema
        ),
        "require_trace_basis_schema": bool(require_trace_basis_schema),
        "source_path": None if source_path is None else str(source_path),
        "source_format": None if source_format is None else str(source_format),
        "polynomial_order": polynomial_order,
        "curved_element_count": curved_element_count,
        "n_fem_nodes": n_fem,
        "volume_tet_count": len(conn_matrix),
        "boundary_triangle_count": len(surface_triangles),
        "fem_node_ids": fem_node_ids,
        "bem_node_ids": bem_node_ids,
        "boundary_nodes_from_triangles": boundary_nodes_from_triangles,
        "trace_shape": [len(trace_rows), len(trace_rows[0]) if trace_rows else 0],
        "trace_row_identity_present": trace_row_identity_present,
        "trace_row_identity": trace_row_identity_normalized,
        "trace_row_identity_mismatch_rows": trace_row_identity_mismatch_rows,
        "operator_trace_row_identity_present": operator_trace_row_identity_present,
        "operator_trace_row_identity": operator_trace_row_identity_normalized,
        "operator_trace_row_identity_mismatch_rows": operator_trace_row_identity_mismatch_rows,
        "boundary_row_identity_present": boundary_row_identity_present,
        "boundary_row_identity": boundary_row_identity_normalized,
        "boundary_row_identity_mismatch_rows": boundary_row_identity_mismatch_rows,
        "operator_boundary_row_identity_present": operator_boundary_row_identity_present,
        "operator_boundary_row_identity": operator_boundary_row_identity_normalized,
        "operator_boundary_row_identity_mismatch_rows": operator_boundary_row_identity_mismatch_rows,
        "boundary_numbers": boundary_numbers,
        "boundary_names": boundary_names,
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "tol": tolerance,
        "notes": [
            "Keep Lukas volume tetrahedra, Gypsilab boundary triangles, and trace gather as one package.",
            "A valid tri/tet mesh is still source-incomplete if mesh_id/export_id/source .vol path are missing.",
            "When recorded, source_file_id or source hash must remain consistent between the .vol package and trace artifact.",
            "The trace matrix and compact surface mesh need their own identities; FEM/BEM coupling reuses both later.",
            "The operator that assembled the trace matrix is a separate artifact from the trace matrix values and from any remote observation field map.",
            "The trace output artifact id/digest/path identify the concrete matrix or table consumed by a notebook or downstream coupling step.",
            "The trace observable id/family state that this output is a boundary trace, not a remote field map, residual curve, or optimization observable.",
            "Once FEM and BEM rows are assembled together, record the coupled system artifact id and digest so a trace package is not mistaken for the actual Schur/coupled solve.",
            "Once the coupled system is solved, keep the linear solver report artifact id, digest, solver name, tolerance, residual norm, and iteration count with the result.",
            "If trace_row_identity is recorded, every row must bind row index, FEM node id, BEM node id, and the trace matrix column.",
            "If operator_trace_row_identity is recorded, it must match the trace row identity and the assembled trace matrix.",
            "If boundary_numbers or boundary_names are recorded, they must stay aligned with the boundary triangles before boundary-condition or BEM-kernel rows use the trace.",
            "If boundary_row_identity is recorded, every boundary triangle row must bind row index, triangle nodes, boundary number, and boundary name.",
            "When coupling leaves the trace-only stage, record coupling kind, formulation id, BEM kernel family, and volume/surface spaces before comparing values.",
            "Record a FEM/BEM coupling convention schema id so stale value-only packages are not reused when the operator/source/normal/kernel convention changed.",
            "Record a trace basis schema id so row identity is not reused with a stale H1-to-surface basis ordering or surface-basis convention.",
            "Readable FEM/BEM notebooks should record the assembly rule id and quadrature rule id before reusing operator rows.",
            "Executed MATLAB notebook/result packages should carry result artifact id, run timestamp, MATLAB version, and a compact timing breakdown before reuse.",
            "Before a FEM/BEM result is reused as an optimization or notebook-default row, record the parameter-set artifact id/digest/path and the objective observable id/family.",
            "The MATLAB/Gypsilab teaching lane is first-order; Netgen curvedelements/high-order blocks require a separate mesh-quality manifest.",
        ],
    }


def netgen_vol_boundary_orientation_trace_package_gate(
    package,
    expected_boundary_orientation=None,
):
    """Check that a .vol FEM/BEM handoff records boundary-triangle orientation.

    The MATLAB/Gypsilab lane can accept Netgen boundary triangles whose stored
    order is inward or outward, but it must say which.  BEM kernels, RWG signs,
    pressure/flux integrals, and exterior normal derivatives should not infer a
    normal convention from filenames or plotting order.
    """

    if not isinstance(package, dict):
        raise ValueError("package must be a dictionary")

    def _as_list(value):
        if value is None:
            return []
        if isinstance(value, (list, tuple)):
            return list(value)
        return [value]

    def _norm(value):
        return str(value or "").strip().lower().replace("-", "_").replace(" ", "_")

    trace = package.get("trace") or {}
    orientation = package.get("boundary_orientation") or package.get("orientation") or {}
    gypsilab = package.get("gypsilab") or {}
    surface_triangles = [
        list(row)
        for row in _as_list(trace.get("surface_triangles") or gypsilab.get("elt") or gypsilab.get("triangles"))
    ]
    boundary_orientation = (
        trace.get("boundary_orientation")
        or orientation.get("boundary_orientation")
        or orientation.get("boundaryOrientation")
    )
    signs = [
        int(value)
        for value in _as_list(
            trace.get("triangle_orientation_signs_to_outward")
            or trace.get("orientation_signs_to_outward")
            or orientation.get("triangle_orientation_signs_to_outward")
            or orientation.get("triangleOrientationSignsToOutward")
        )
    ]
    adjacent = [
        int(value)
        for value in _as_list(
            trace.get("adjacent_tet_indices")
            or orientation.get("adjacent_tet_indices")
            or orientation.get("adjacentTetIndices")
        )
    ]
    rows = _as_list(orientation.get("rows") or trace.get("orientation_rows"))
    triangle_count = len(surface_triangles) or len(rows) or len(signs)
    allowed = {"outward", "inward", "mixed", "unknown_or_open"}
    checks = {
        "surface_triangles_recorded": triangle_count > 0,
        "boundary_orientation_recorded": _norm(boundary_orientation) in allowed,
        "orientation_signs_recorded": len(signs) == triangle_count and bool(signs),
        "orientation_signs_are_unit": bool(signs) and all(abs(sign) == 1 for sign in signs),
        "adjacent_tet_indices_recorded": len(adjacent) == triangle_count and bool(adjacent),
        "no_orphan_boundary_triangles": bool(adjacent) and all(index > 0 for index in adjacent),
    }
    if _norm(boundary_orientation) in {"outward", "inward"} and signs:
        expected_sign = 1 if _norm(boundary_orientation) == "outward" else -1
        checks["uniform_orientation_matches_signs"] = all(sign == expected_sign for sign in signs)
    if expected_boundary_orientation is not None:
        checks["expected_boundary_orientation_matches"] = (
            _norm(boundary_orientation) == _norm(expected_boundary_orientation)
        )

    return {
        "policy": "netgen_vol_boundary_orientation_trace_package_gate",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "boundary_orientation": None if boundary_orientation is None else str(boundary_orientation),
        "expected_boundary_orientation": (
            None if expected_boundary_orientation is None else str(expected_boundary_orientation)
        ),
        "triangle_count": triangle_count,
        "triangle_orientation_signs_to_outward": signs,
        "adjacent_tet_indices": adjacent,
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "notes": [
            "Boundary triangle order is allowed to be inward or outward only when the handoff says so explicitly.",
            "Use orientation_sign_to_outward before exterior BEM normals, RWG signs, flux integrals, or pressure resultants.",
        ],
    }


def netgen_vol_fem_bem_normal_flux_sign_package_gate(
    package,
    expected_normal_convention=None,
    tol=1.0e-12,
):
    """Check that FEM/BEM normal-flux rows consume .vol orientation signs.

    The boundary-orientation gate records whether stored Netgen triangles are
    inward or outward.  This companion gate verifies that a scalar/vector flux
    handoff actually uses ``orientation_sign_to_outward`` before comparing
    exterior normal derivatives, pressure resultants, or charge/current rows.
    """

    if not isinstance(package, dict):
        raise ValueError("package must be a dictionary")

    def _as_list(value):
        if value is None:
            return []
        if isinstance(value, (list, tuple)):
            return list(value)
        return [value]

    def _float_list(value):
        return [float(item) for item in _as_list(value)]

    trace = package.get("trace") or {}
    orientation = package.get("boundary_orientation") or package.get("orientation") or {}
    flux = package.get("normal_flux") or package.get("flux") or {}
    normal_convention = (
        flux.get("normal_convention")
        or flux.get("normalConvention")
        or package.get("normal_convention")
        or package.get("normalConvention")
        or trace.get("normal_convention")
        or orientation.get("normal_convention")
    )
    normal_key = None if normal_convention is None else str(normal_convention).strip().lower().replace("-", "_").replace(" ", "_")
    expected_normal_key = (
        None
        if expected_normal_convention is None
        else str(expected_normal_convention).strip().lower().replace("-", "_").replace(" ", "_")
    )
    signs = [
        int(value)
        for value in _as_list(
            flux.get("triangle_orientation_signs_to_outward")
            or trace.get("triangle_orientation_signs_to_outward")
            or trace.get("orientation_signs_to_outward")
            or orientation.get("triangle_orientation_signs_to_outward")
            or orientation.get("triangleOrientationSignsToOutward")
        )
    ]
    stored_flux = _float_list(
        flux.get("stored_normal_flux")
        or flux.get("storedNormalFlux")
        or package.get("stored_normal_flux")
    )
    corrected_flux = _float_list(
        flux.get("orientation_corrected_normal_flux")
        or flux.get("outward_normal_flux")
        or flux.get("correctedNormalFlux")
        or package.get("orientation_corrected_normal_flux")
    )
    reference_flux = _float_list(
        flux.get("outward_normal_flux_reference")
        or flux.get("outwardReferenceFlux")
        or package.get("outward_normal_flux_reference")
    )
    triangle_count = max(len(signs), len(stored_flux), len(corrected_flux), len(reference_flux))
    tolerance = float(tol)
    if tolerance < 0.0:
        raise ValueError("tol must be non-negative")

    sign_corrected = [
        sign * value
        for sign, value in zip(signs, stored_flux)
    ]
    corrected_residuals = [
        abs(a - b)
        for a, b in zip(sign_corrected, corrected_flux)
    ]
    reference_residuals = [
        abs(a - b)
        for a, b in zip(corrected_flux, reference_flux)
    ]
    closed_surface_sum = float(
        flux.get(
            "closed_surface_flux_sum",
            package.get("closed_surface_flux_sum", sum(corrected_flux) if corrected_flux else 0.0),
        )
    )
    expected_closed_surface_sum = float(
        flux.get(
            "expected_closed_surface_flux_sum",
            package.get("expected_closed_surface_flux_sum", 0.0),
        )
    )
    closed_surface_residual = abs(closed_surface_sum - expected_closed_surface_sum)
    checks = {
        "orientation_signs_recorded": len(signs) == triangle_count and bool(signs),
        "orientation_signs_are_unit": bool(signs) and all(abs(sign) == 1 for sign in signs),
        "stored_flux_rows_recorded": len(stored_flux) == triangle_count and bool(stored_flux),
        "corrected_flux_rows_recorded": len(corrected_flux) == triangle_count and bool(corrected_flux),
        "outward_reference_rows_recorded": len(reference_flux) == triangle_count and bool(reference_flux),
        "sign_corrected_flux_matches_rows": (
            len(corrected_residuals) == triangle_count
            and max(corrected_residuals, default=float("inf")) <= tolerance
        ),
        "corrected_flux_matches_outward_reference": (
            len(reference_residuals) == triangle_count
            and max(reference_residuals, default=float("inf")) <= tolerance
        ),
        "closed_surface_flux_balance_ok": closed_surface_residual <= tolerance,
    }
    if expected_normal_key is not None:
        checks["normal_convention_recorded"] = bool(normal_key)
        checks["normal_convention_matches_expected"] = normal_key == expected_normal_key
    return {
        "policy": "netgen_vol_fem_bem_normal_flux_sign_package_gate",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "triangle_count": triangle_count,
        "triangle_orientation_signs_to_outward": signs,
        "stored_normal_flux": stored_flux,
        "sign_corrected_normal_flux": sign_corrected,
        "orientation_corrected_normal_flux": corrected_flux,
        "outward_normal_flux_reference": reference_flux,
        "max_sign_correction_abs_error": max(corrected_residuals, default=None),
        "max_outward_reference_abs_error": max(reference_residuals, default=None),
        "closed_surface_flux_sum": closed_surface_sum,
        "expected_closed_surface_flux_sum": expected_closed_surface_sum,
        "closed_surface_flux_residual": closed_surface_residual,
        "normal_convention": normal_key,
        "expected_normal_convention": expected_normal_key,
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "tol": tolerance,
        "notes": [
            "Do not feed stored .vol triangle normals directly into exterior BEM quantities.",
            "Use orientation_sign_to_outward row-by-row before normal derivatives, flux integrals, or surface source balances.",
        ],
    }


def quarter_wave_directional_coupler_gate(coupling, z0=50.0):
    """Return algebraic checks for an ideal coupled-line directional coupler.

    ``coupling`` is the linear voltage coupling coefficient.  At the design
    electrical length, the ideal lossless coupler has matched and isolated
    ports, ``|S21|^2 + |S31|^2 = 1``, and even/odd impedances satisfying
    ``Z0e * Z0o = z0^2``.
    """

    c = float(coupling)
    z = float(z0)
    if not (0.0 < c < 1.0):
        raise ValueError("coupling must be in (0, 1)")
    if z <= 0.0:
        raise ValueError("z0 must be > 0")
    z_even = z * math.sqrt((1.0 + c) / (1.0 - c))
    z_odd = z * math.sqrt((1.0 - c) / (1.0 + c))
    s11 = 0.0 + 0.0j
    s41 = 0.0 + 0.0j
    s21 = 1j * math.sqrt(1.0 - c * c)
    s31 = complex(-c, 0.0)
    product_error = z_even * z_odd - z * z
    power_sum = abs(s21) ** 2 + abs(s31) ** 2
    return {
        "policy": "quarter_wave_directional_coupler_algebra_gate",
        "coupling": c,
        "z0": z,
        "z0_even": z_even,
        "z0_odd": z_odd,
        "impedance_product": z_even * z_odd,
        "impedance_product_error": product_error,
        "coupling_db": -20.0 * math.log10(c),
        "s11": _complex_row(s11),
        "s21": _complex_row(s21),
        "s31": _complex_row(s31),
        "s41": _complex_row(s41),
        "through_phase_deg": math.degrees(cmath.phase(s21)),
        "coupled_phase_deg": math.degrees(cmath.phase(s31)),
        "power_sum": power_sum,
        "matched": abs(s11) <= 1.0e-15,
        "isolated": abs(s41) <= 1.0e-15,
        "lossless": abs(power_sum - 1.0) <= 1.0e-15,
        "status": "ok" if abs(product_error) <= 1.0e-12 and abs(power_sum - 1.0) <= 1.0e-15 else "needs_attention",
    }


def branch_line_hybrid_gate(z0=50.0, tol=1.0e-12):
    """Return algebraic checks for an ideal 90-degree branch-line hybrid.

    Port convention follows the common quadrature hybrid row: port 1 input,
    port 2 through, port 3 coupled, port 4 isolated.  At the design frequency,
    ``|S21| = |S31| = 1/sqrt(2)``, ``S11 = S41 = 0``, and the through/coupled
    phase difference is -90 degrees.
    """

    z = float(z0)
    if z <= 0.0:
        raise ValueError("z0 must be > 0")
    s = 1.0 / math.sqrt(2.0)
    smatrix = [
        [0.0 + 0.0j, -1j * s, -s + 0.0j, 0.0 + 0.0j],
        [-1j * s, 0.0 + 0.0j, 0.0 + 0.0j, -s + 0.0j],
        [-s + 0.0j, 0.0 + 0.0j, 0.0 + 0.0j, -1j * s],
        [0.0 + 0.0j, -s + 0.0j, -1j * s, 0.0 + 0.0j],
    ]
    column_powers = [
        sum(abs(smatrix[row][col]) ** 2 for row in range(4))
        for col in range(4)
    ]
    max_column_power_error = max(abs(power - 1.0) for power in column_powers)
    max_orthogonality_error = 0.0
    for col_a in range(4):
        for col_b in range(col_a + 1, 4):
            dot = sum(smatrix[row][col_a].conjugate() * smatrix[row][col_b] for row in range(4))
            max_orthogonality_error = max(max_orthogonality_error, abs(dot))
    s11 = smatrix[0][0]
    s21 = smatrix[1][0]
    s31 = smatrix[2][0]
    s41 = smatrix[3][0]
    split_power_sum = abs(s21) ** 2 + abs(s31) ** 2
    phase_difference = math.degrees(cmath.phase(s31) - cmath.phase(s21))
    while phase_difference <= -180.0:
        phase_difference += 360.0
    while phase_difference > 180.0:
        phase_difference -= 360.0
    checks = {
        "matched": abs(s11) <= float(tol),
        "isolated": abs(s41) <= float(tol),
        "equal_split": abs(abs(s21) - s) <= float(tol) and abs(abs(s31) - s) <= float(tol),
        "lossless_split": abs(split_power_sum - 1.0) <= float(tol),
        "quadrature_phase": abs(phase_difference + 90.0) <= 1.0e-9,
        "unitary": max_column_power_error <= float(tol) and max_orthogonality_error <= float(tol),
    }
    return {
        "policy": "branch_line_hybrid_quadrature_gate",
        "z0": z,
        "through_branch_impedance": z,
        "shunt_branch_impedance": z / math.sqrt(2.0),
        "s11": _complex_row(s11),
        "s21": _complex_row(s21),
        "s31": _complex_row(s31),
        "s41": _complex_row(s41),
        "split_power_sum": split_power_sum,
        "through_phase_deg": math.degrees(cmath.phase(s21)),
        "coupled_phase_deg": math.degrees(cmath.phase(s31)),
        "phase_difference_deg": phase_difference,
        "column_powers": column_powers,
        "max_column_power_error": max_column_power_error,
        "max_orthogonality_error": max_orthogonality_error,
        "checks": checks,
        "status": "ok" if all(checks.values()) else "needs_attention",
    }


def dq_to_three_phase_currents(id_current, iq_current, theta_e_rad, phases=("U", "V", "W")):
    """Return balanced three-phase currents from amplitude-invariant dq values.

    The convention is ``i_phase = id*cos(theta_phase) - iq*sin(theta_phase)``
    with phase angles ``U=theta``, ``V=theta-120deg``, and ``W=theta+120deg``.
    It is intentionally a small handoff contract for motor FEM and drive
    notebooks before any solver-specific circuit-current API is called.
    """

    phase_names = tuple(phases)
    if len(phase_names) != 3:
        raise ValueError("exactly three phase names are required")
    theta = float(theta_e_rad)
    angles = (theta, theta - 2.0 * math.pi / 3.0, theta + 2.0 * math.pi / 3.0)
    d = float(id_current)
    q = float(iq_current)
    return {
        phase: d * math.cos(angle) - q * math.sin(angle)
        for phase, angle in zip(phase_names, angles)
    }


def three_phase_currents_to_dq_summary(
    currents,
    theta_e_rad,
    expected_id=None,
    expected_iq=None,
    phases=("U", "V", "W"),
    tol=1.0e-12,
):
    """Recover dq currents and balance checks from a three-phase current row."""

    phase_names = tuple(phases)
    if len(phase_names) != 3:
        raise ValueError("exactly three phase names are required")
    theta = float(theta_e_rad)
    angles = (theta, theta - 2.0 * math.pi / 3.0, theta + 2.0 * math.pi / 3.0)
    values = [float(currents[phase]) for phase in phase_names]
    d = (2.0 / 3.0) * sum(value * math.cos(angle) for value, angle in zip(values, angles))
    q = -(2.0 / 3.0) * sum(value * math.sin(angle) for value, angle in zip(values, angles))
    i0 = sum(values) / 3.0
    abc_square_sum = sum(value * value for value in values)
    dq_square_sum = 1.5 * (d * d + q * q)
    payload = {
        "policy": "three_phase_dq_current_handoff_gate",
        "phase_order": list(phase_names),
        "theta_e_rad": theta,
        "id": d,
        "iq": q,
        "i0": i0,
        "abc_square_sum": abc_square_sum,
        "dq_square_sum_scaled": dq_square_sum,
        "abc_vs_dq_square_sum_error": abs(abc_square_sum - dq_square_sum),
        "zero_sequence_abs": abs(i0),
        "tol": float(tol),
    }
    checks = {
        "zero_sequence_ok": abs(i0) <= float(tol),
        "abc_square_sum_ok": abs(abc_square_sum - dq_square_sum) <= float(tol),
    }
    if expected_id is not None:
        payload["expected_id"] = float(expected_id)
        payload["id_abs_error"] = abs(d - float(expected_id))
        checks["id_ok"] = payload["id_abs_error"] <= float(tol)
    if expected_iq is not None:
        payload["expected_iq"] = float(expected_iq)
        payload["iq_abs_error"] = abs(q - float(expected_iq))
        checks["iq_ok"] = payload["iq_abs_error"] <= float(tol)
    payload["checks"] = checks
    payload["status"] = "ok" if all(checks.values()) else "needs_attention"
    return payload


def femm_static_current_circuit_rows_gate(
    currents,
    theta_e_rad,
    circuit_rows,
    expected_id=None,
    expected_iq=None,
    phases=("U", "V", "W"),
    expected_current_kind="instantaneous",
    tol=1.0e-12,
):
    """Check solver-ready static circuit-current rows before a FEM snapshot.

    The physics part is solver-independent: a balanced U/V/W row must round-trip
    through dq.  The handoff part records the small FEM-style contract that is
    easy to lose in converters: every phase row needs an explicit circuit name,
    an instantaneous snapshot current (not an RMS label), and a non-zero turns
    entry before a solver-specific API such as ``mi_modifycircprop`` is called.
    """

    phase_names = tuple(phases)
    if len(phase_names) != 3:
        raise ValueError("exactly three phase names are required")
    tolerance = float(tol)
    if tolerance < 0.0:
        raise ValueError("tol must be non-negative")

    def _normalize_kind(kind):
        key = str(kind).strip().lower().replace("-", "_").replace(" ", "_")
        aliases = {
            "instantaneous": "instantaneous",
            "instant": "instantaneous",
            "snapshot": "instantaneous",
            "sample": "instantaneous",
            "sampled": "instantaneous",
            "peak": "instantaneous",
            "peak_sample": "instantaneous",
            "instantaneous_sample": "instantaneous",
            "rms": "rms",
            "root_mean_square": "rms",
        }
        return aliases.get(key, key)

    def _first_present(row, names):
        for name in names:
            if name in row and row[name] is not None:
                return row[name]
        return None

    dq = three_phase_currents_to_dq_summary(
        currents,
        theta_e_rad,
        expected_id=expected_id,
        expected_iq=expected_iq,
        phases=phase_names,
        tol=tolerance,
    )
    expected_kind = _normalize_kind(expected_current_kind)

    row_by_phase = {}
    duplicate_phases = []
    unknown_phases = []
    if hasattr(circuit_rows, "items"):
        iterator = circuit_rows.items()
    else:
        iterator = enumerate(circuit_rows)
    for key, row in iterator:
        if not isinstance(row, dict):
            row = {"current_A": row}
        phase = row.get("phase", key)
        if phase not in phase_names:
            unknown_phases.append(phase)
            continue
        if phase in row_by_phase:
            duplicate_phases.append(phase)
        row_by_phase[phase] = dict(row)

    details = []
    current_errors = []
    for phase in phase_names:
        row = row_by_phase.get(phase)
        if row is None:
            details.append({"phase": phase, "present": False})
            continue
        circuit_name = _first_present(row, ("circuit_name", "circuit", "name"))
        current_value = _first_present(row, ("current_A", "current", "amps", "current_amps"))
        turns = _first_present(row, ("turns", "series_turns", "nturns"))
        multiplier = _first_present(row, ("current_multiplier", "multiplier"))
        multiplier = 1.0 if multiplier is None else float(multiplier)
        kind = row.get("current_kind")
        normalized_kind = None if kind is None else _normalize_kind(kind)
        current_error = None
        if current_value is not None:
            current_error = abs(float(current_value) - float(currents[phase]) * multiplier)
            current_errors.append(current_error)
        details.append({
            "phase": phase,
            "present": True,
            "circuit_name": circuit_name,
            "current_A": None if current_value is None else float(current_value),
            "expected_current_A": float(currents[phase]) * multiplier,
            "current_multiplier": multiplier,
            "current_abs_error_A": current_error,
            "turns": None if turns is None else float(turns),
            "current_kind": kind,
            "normalized_current_kind": normalized_kind,
        })

    checks = {
        "phase_set_ok": set(row_by_phase) == set(phase_names),
        "no_unknown_phase_rows": not unknown_phases,
        "no_duplicate_phase_rows": not duplicate_phases,
        "circuit_names_present": all(bool(row.get("circuit_name")) or bool(row.get("circuit")) or bool(row.get("name")) for row in row_by_phase.values()),
        "current_values_present": all(_first_present(row, ("current_A", "current", "amps", "current_amps")) is not None for row in row_by_phase.values()) and set(row_by_phase) == set(phase_names),
        "current_kind_matches": all(
            _normalize_kind(row.get("current_kind")) == expected_kind
            for row in row_by_phase.values()
        ) and set(row_by_phase) == set(phase_names),
        "turns_present": all(_first_present(row, ("turns", "series_turns", "nturns")) is not None for row in row_by_phase.values()) and set(row_by_phase) == set(phase_names),
        "turns_nonzero": all(
            _first_present(row, ("turns", "series_turns", "nturns")) is not None
            and abs(float(_first_present(row, ("turns", "series_turns", "nturns")))) > 0.0
            for row in row_by_phase.values()
        ) and set(row_by_phase) == set(phase_names),
        "circuit_currents_match": bool(current_errors) and len(current_errors) == len(phase_names) and max(current_errors) <= tolerance,
        "dq_roundtrip_ok": dq["status"] == "ok",
    }
    return {
        "policy": "femm_static_current_circuit_rows_gate",
        "phase_order": list(phase_names),
        "theta_e_rad": float(theta_e_rad),
        "expected_current_kind": expected_kind,
        "currents": {phase: float(currents[phase]) for phase in phase_names},
        "dq": dq,
        "circuit_rows": details,
        "unknown_phases": unknown_phases,
        "duplicate_phases": duplicate_phases,
        "max_circuit_current_abs_error_A": None if not current_errors else max(current_errors),
        "tol": tolerance,
        "checks": checks,
        "status": "ok" if all(checks.values()) else "needs_attention",
    }


def femm_block_label_source_contract_gate(
    block_rows,
    required_regions=None,
    allowed_source_kinds=("air", "passive", "coil", "pm"),
):
    """Check block-label rows before a finite-element motor model handoff.

    The row contract is intentionally solver-independent: every region needs a
    stable name, material, and group id; coil rows need circuit and non-zero
    turns; PM rows need an explicit magnetization direction.  This prevents
    converters from inferring sources from material names after the fact.
    """

    rows_in = list(block_rows)
    if not rows_in:
        raise ValueError("block_rows must not be empty")
    allowed = {str(kind).strip().lower().replace("-", "_").replace(" ", "_") for kind in allowed_source_kinds}
    if not allowed:
        raise ValueError("allowed_source_kinds must not be empty")
    required = [] if required_regions is None else [str(region) for region in required_regions]

    normalized_rows = []
    region_counts = {}
    unknown_source_kinds = []
    for index, row in enumerate(rows_in, start=1):
        if not isinstance(row, dict):
            raise ValueError("each block row must be a dictionary")
        region = str(row.get("region") or row.get("name") or row.get("label") or "").strip()
        material = str(row.get("material") or row.get("material_name") or "").strip()
        source_kind = str(row.get("source_kind") or row.get("kind") or "passive").strip().lower().replace("-", "_").replace(" ", "_")
        group = row.get("group", row.get("group_id"))
        circuit_name = row.get("circuit_name", row.get("circuit"))
        turns = row.get("turns", row.get("nturns", row.get("series_turns")))
        magnetization = row.get("magnetization_deg", row.get("magnetization_direction_deg", row.get("magdir_deg")))
        if source_kind not in allowed:
            unknown_source_kinds.append(source_kind)
        if region:
            region_counts[region] = region_counts.get(region, 0) + 1
        normalized_rows.append({
            "index": index,
            "region": region,
            "material": material,
            "group": None if group is None else int(group),
            "source_kind": source_kind,
            "circuit_name": None if circuit_name is None else str(circuit_name),
            "turns": None if turns is None else float(turns),
            "magnetization_deg": None if magnetization is None else float(magnetization),
        })

    regions = {row["region"] for row in normalized_rows if row["region"]}
    duplicate_regions = sorted(name for name, count in region_counts.items() if count > 1)
    missing_required_regions = sorted(region for region in required if region not in regions)
    coil_rows = [row for row in normalized_rows if row["source_kind"] == "coil"]
    pm_rows = [row for row in normalized_rows if row["source_kind"] == "pm"]
    passive_source_rows = [
        row for row in normalized_rows
        if row["source_kind"] in {"air", "passive"}
        and (row["circuit_name"] or row["turns"] not in (None, 0.0) or row["magnetization_deg"] is not None)
    ]
    source_counts = {}
    for row in normalized_rows:
        source_counts[row["source_kind"]] = source_counts.get(row["source_kind"], 0) + 1
    checks = {
        "region_names_present": all(row["region"] for row in normalized_rows),
        "region_names_unique": not duplicate_regions,
        "required_regions_present": not missing_required_regions,
        "materials_present": all(row["material"] for row in normalized_rows),
        "groups_present": all(row["group"] is not None for row in normalized_rows),
        "source_kinds_allowed": not unknown_source_kinds,
        "coil_rows_have_circuit": all(bool(row["circuit_name"]) for row in coil_rows),
        "coil_rows_have_nonzero_turns": all(row["turns"] is not None and abs(row["turns"]) > 0.0 for row in coil_rows),
        "pm_rows_have_magnetization_direction": all(row["magnetization_deg"] is not None for row in pm_rows),
        "air_passive_rows_have_no_source_metadata": not passive_source_rows,
    }
    return {
        "policy": "femm_block_label_source_contract_gate",
        "n_rows": len(normalized_rows),
        "required_regions": required,
        "missing_required_regions": missing_required_regions,
        "duplicate_regions": duplicate_regions,
        "unknown_source_kinds": unknown_source_kinds,
        "source_counts": dict(sorted(source_counts.items())),
        "rows": normalized_rows,
        "passive_source_rows": passive_source_rows,
        "checks": checks,
        "status": "ok" if all(checks.values()) else "needs_attention",
        "version_note": (
            "Use this after FEMM .fem block-label parsing and before ModelIR or "
            "radia-ngsolve emission so coils and PMs stay explicit sources."
        ),
    }


def femm_group_motion_selection_gate(
    rows,
    expected_group_id,
    required_entity_kinds=("block_label", "segment", "arc_segment"),
    required_motion_command="mi_moverotate",
):
    """Check FEMM group metadata before rotating or translating a moving part.

    FEMM motion commands operate on the current selection. For rotor sweeps and
    cogging examples, every moving entity must carry the intended nonzero group
    id before selecting that group; otherwise labels can move while geometry or
    sources remain stale.
    """

    entity_rows = [dict(row) for row in rows]
    if not entity_rows:
        raise ValueError("rows must not be empty")

    def _norm(value):
        return str(value or "").strip().lower().replace("-", "_").replace(" ", "_")

    expected_group = int(expected_group_id)
    if expected_group <= 0:
        raise ValueError("expected_group_id must be positive")
    required = {_norm(kind) for kind in required_entity_kinds if str(kind).strip()}
    motion_command = _norm(required_motion_command)
    kinds = []
    bad_group_rows = []
    group_zero_rows = []
    unselected_rows = []
    missing_identity_rows = []
    bad_motion_rows = []
    for index, row in enumerate(entity_rows, start=1):
        kind = _norm(row.get("entity_kind", row.get("kind")))
        kinds.append(kind)
        group_id = row.get("group_id", row.get("group"))
        try:
            group_value = int(group_id)
        except (TypeError, ValueError):
            group_value = None
        if group_value != expected_group:
            bad_group_rows.append({"index": index, "kind": kind, "group_id": group_id})
        if group_value == 0:
            group_zero_rows.append({"index": index, "kind": kind})
        if row.get("selected_for_motion") is not True:
            unselected_rows.append({"index": index, "kind": kind})
        if not (row.get("entity_id") or row.get("name") or row.get("region_name")):
            missing_identity_rows.append(index)
        row_motion = _norm(row.get("motion_command", required_motion_command))
        if row_motion != motion_command:
            bad_motion_rows.append({"index": index, "kind": kind, "motion_command": row.get("motion_command")})

    kind_set = set(kinds)
    checks = {
        "rows_present": bool(entity_rows),
        "required_entity_kinds_present": required.issubset(kind_set),
        "expected_group_is_positive": expected_group > 0,
        "all_rows_use_expected_group": not bad_group_rows,
        "no_group_zero_rows": not group_zero_rows,
        "all_rows_selected_for_motion": not unselected_rows,
        "entity_identity_recorded": not missing_identity_rows,
        "motion_command_matches": not bad_motion_rows,
    }
    return {
        "policy": "femm_group_motion_selection_gate",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "expected_group_id": expected_group,
        "required_entity_kinds": sorted(required),
        "present_entity_kinds": sorted(kind_set),
        "required_motion_command": required_motion_command,
        "bad_group_rows": bad_group_rows,
        "group_zero_rows": group_zero_rows,
        "unselected_rows": unselected_rows,
        "missing_identity_rows": missing_identity_rows,
        "bad_motion_rows": bad_motion_rows,
        "checks": checks,
        "notes": [
            "Run before mi_selectgroup/mi_moverotate rotor sweeps so geometry, block labels, and magnetization sources move together.",
            "A row with group 0 is not solver-ready motion metadata even if the label text looks correct.",
        ],
    }


def femm_pm_magnetization_convention_gate(
    pm_rows,
    required_regions=None,
    allowed_frames=("global_xy", "rotor_xy", "local_radial", "local_tangential"),
):
    """Check permanent-magnet direction rows before emission.

    FEMM stores a PM source through ``mi_setblockprop(..., magdir, ...)``.
    The numeric angle is not solver-ready by itself: it must carry degree
    units, a coordinate frame, and a Br/Hc strength column before conversion to
    a radia-ngsolve magnetization vector.
    """

    rows_in = list(pm_rows)
    if not rows_in:
        raise ValueError("pm_rows must not be empty")
    required = [] if required_regions is None else [str(region) for region in required_regions]
    allowed = {str(frame).strip().lower().replace("-", "_").replace(" ", "_") for frame in allowed_frames}
    if not allowed:
        raise ValueError("allowed_frames must not be empty")

    normalized_rows = []
    region_counts = {}
    bad_angle_regions = []
    bad_frame_regions = []
    missing_strength_regions = []
    for index, row in enumerate(rows_in, start=1):
        if not isinstance(row, dict):
            raise ValueError("each PM row must be a dictionary")
        region = str(row.get("region") or row.get("name") or row.get("label") or "").strip()
        unit = str(row.get("angle_unit") or row.get("magdir_unit") or "deg").strip().lower()
        unit = unit.replace("degrees", "deg").replace("degree", "deg")
        frame = str(row.get("frame") or row.get("angle_frame") or row.get("coordinate_frame") or "").strip().lower()
        frame = frame.replace("-", "_").replace(" ", "_")
        magdir = row.get("magdir_deg", row.get("magnetization_deg", row.get("magnetization_direction_deg")))
        br = row.get("br_T", row.get("Br_T", row.get("remanence_T")))
        hc = row.get("hc_A_per_m", row.get("Hc_A_per_m", row.get("coercivity_A_per_m")))
        angle = None
        if magdir is not None:
            try:
                angle = float(magdir)
            except (TypeError, ValueError):
                angle = None
        if region:
            region_counts[region] = region_counts.get(region, 0) + 1
        if unit != "deg" or angle is None or not math.isfinite(angle):
            bad_angle_regions.append(region or f"row_{index}")
        if frame not in allowed:
            bad_frame_regions.append(region or f"row_{index}")
        strength_present = br is not None or hc is not None
        if not strength_present:
            missing_strength_regions.append(region or f"row_{index}")
        theta = math.radians(angle or 0.0)
        normalized_rows.append({
            "index": index,
            "region": region,
            "angle_unit": unit,
            "frame": frame or None,
            "magdir_deg": angle,
            "unit_vector_xy": [math.cos(theta), math.sin(theta)] if angle is not None else None,
            "br_T": None if br is None else float(br),
            "hc_A_per_m": None if hc is None else float(hc),
            "strength_present": strength_present,
        })

    regions = {row["region"] for row in normalized_rows if row["region"]}
    duplicate_regions = sorted(region for region, count in region_counts.items() if count > 1)
    missing_required_regions = sorted(region for region in required if region not in regions)
    vector_errors = []
    for row in normalized_rows:
        vec = row["unit_vector_xy"]
        if vec is None:
            continue
        vector_errors.append(abs(math.hypot(vec[0], vec[1]) - 1.0))
    checks = {
        "regions_present": all(row["region"] for row in normalized_rows),
        "regions_unique": not duplicate_regions,
        "required_regions_present": not missing_required_regions,
        "angles_are_degrees_and_finite": not bad_angle_regions,
        "frames_allowed": not bad_frame_regions,
        "strength_present": not missing_strength_regions,
        "unit_vectors_normalized": (max(vector_errors) if vector_errors else 0.0) < 1.0e-12,
    }
    return {
        "policy": "femm_pm_magnetization_convention_gate",
        "n_rows": len(normalized_rows),
        "allowed_frames": sorted(allowed),
        "required_regions": required,
        "missing_required_regions": missing_required_regions,
        "duplicate_regions": duplicate_regions,
        "bad_angle_regions": sorted(bad_angle_regions),
        "bad_frame_regions": sorted(bad_frame_regions),
        "missing_strength_regions": sorted(missing_strength_regions),
        "max_unit_vector_norm_error": max(vector_errors) if vector_errors else 0.0,
        "rows": normalized_rows,
        "checks": checks,
        "status": "ok" if all(checks.values()) else "needs_attention",
        "version_note": (
            "Use this after femm_block_label_source_contract_gate for PM rows: "
            "FEMM magdir must be degrees plus an explicit coordinate frame and "
            "Br/Hc strength before radia-ngsolve vector emission."
        ),
    }


def jmag_motor_table_column_metadata_gate(
    metadata,
    required_columns=(),
    *,
    angle_column=None,
    torque_column=None,
    angle_unit=None,
    angle_basis=None,
    pole_pairs=None,
    symmetry_factor=1.0,
    current_basis=None,
    torque_sign_convention=None,
):
    """Check motor table column metadata before reading values."""

    if not isinstance(metadata, dict):
        raise ValueError("metadata must be a dictionary")
    columns_raw = metadata.get("columns", metadata.get("column_names", ()))
    columns = [str(column).strip() for column in columns_raw]
    if not columns:
        raise ValueError("metadata must include columns")
    column_set = set(columns)
    required = [str(column).strip() for column in required_columns]
    missing_required = [column for column in required if column not in column_set]
    angle_col = str(angle_column or metadata.get("angle_column") or "").strip()
    torque_col = str(torque_column or metadata.get("torque_column") or "").strip()
    unit = str(angle_unit or metadata.get("angle_unit") or "").strip().lower().replace("degrees", "deg").replace("degree", "deg")
    basis = str(angle_basis or metadata.get("angle_basis") or "").strip().lower().replace("-", "_").replace(" ", "_")
    pp_raw = pole_pairs if pole_pairs is not None else metadata.get("pole_pairs")
    sf = float(metadata.get("symmetry_factor", symmetry_factor))
    current = current_basis if current_basis is not None else metadata.get("current_basis")
    current_key = None if current is None else str(current).strip().lower().replace("-", "_").replace(" ", "_")
    sign = torque_sign_convention if torque_sign_convention is not None else metadata.get("torque_sign_convention")
    sign_key = None if sign is None else str(sign).strip().lower().replace("-", "_").replace(" ", "_")
    angle_units = {"deg", "rad"}
    angle_bases = {"mechanical", "electrical"}
    current_bases = {"instantaneous", "peak", "rms"}
    sign_conventions = {"positive_motoring", "positive_generator", "as_exported"}
    pole_pair_value = None if pp_raw is None else float(pp_raw)
    checks = {
        "required_columns_present": not missing_required,
        "angle_column_present": bool(angle_col) and angle_col in column_set,
        "torque_column_present": (not torque_col) or torque_col in column_set,
        "angle_unit_valid": unit in angle_units,
        "angle_basis_valid": basis in angle_bases,
        "pole_pairs_positive": pole_pair_value is not None and pole_pair_value > 0.0,
        "symmetry_factor_positive": sf > 0.0,
        "current_basis_valid": current_key is None or current_key in current_bases,
        "torque_sign_convention_valid": sign_key is None or sign_key in sign_conventions,
    }
    return {
        "policy": "jmag_motor_table_column_metadata_gate",
        "columns": columns,
        "required_columns": required,
        "missing_required_columns": missing_required,
        "angle_column": angle_col,
        "torque_column": torque_col or None,
        "angle_unit": unit or None,
        "angle_basis": basis or None,
        "pole_pairs": pole_pair_value,
        "symmetry_factor": sf,
        "current_basis": current_key,
        "torque_sign_convention": sign_key,
        "checks": checks,
        "status": "ok" if all(checks.values()) else "needs_attention",
        "version_note": (
            "Use this before JMAG torque/current/efficiency table parsing so "
            "column units and basis are explicit before numerical gates run."
        ),
    }


def jmag_force_table_metadata_gate(
    metadata,
    required_columns=(),
    *,
    force_columns=None,
    position_columns=None,
    force_unit=None,
    position_unit=None,
    component_frame=None,
    projection_axis=None,
    source_tool=None,
    symmetry_factor=1.0,
    force_sign_convention=None,
    force_kind=None,
    quantity_dimension=None,
    expected_case_id=None,
    expected_study_id=None,
    expected_operating_point_id=None,
    expected_analysis_type=None,
    expected_frequency_hz=None,
    expected_export_artifact_id=None,
    expected_result_set_id=None,
    expected_target_region_id=None,
    expected_target_region_name=None,
    expected_target_material=None,
    expected_target_region_artifact_id=None,
    expected_target_region_geometry_digest=None,
    expected_target_region_centroid_xyz_m=None,
    expected_mesh_id=None,
    expected_solver_run_id=None,
    expected_result_revision_id=None,
    expected_solver_setup_artifact_id=None,
    expected_material_state_artifact_id=None,
    expected_material_state_digest=None,
    expected_excitation_source_artifact_id=None,
    expected_current_definition_method=None,
    expected_export_trace_id=None,
    expected_export_command_digest=None,
    expected_export_output_artifact_id=None,
    expected_export_output_digest=None,
    expected_force_observable_id=None,
    expected_force_observable_family=None,
    expected_component_frame=None,
    expected_projection_axis=None,
    expected_force_sign_convention=None,
    force_report_method=None,
    expected_force_report_method=None,
    require_export_command_trace=False,
    require_export_output_artifact=False,
    table_rows=None,
    min_table_rows=None,
    require_unique_row_identity=None,
):
    """Check force table metadata before reading force values."""

    if not isinstance(metadata, dict):
        raise ValueError("metadata must be a dictionary")
    columns_raw = metadata.get("columns", metadata.get("column_names", ()))
    columns = [str(column).strip() for column in columns_raw]
    if not columns:
        raise ValueError("metadata must include columns")
    column_set = set(columns)
    required = [str(column).strip() for column in required_columns]
    missing_required = [column for column in required if column not in column_set]
    force_cols_raw = force_columns if force_columns is not None else metadata.get("force_columns", ())
    if isinstance(force_cols_raw, str):
        force_cols = [force_cols_raw.strip()]
    else:
        force_cols = [str(column).strip() for column in force_cols_raw if str(column).strip()]
    position_cols_raw = position_columns if position_columns is not None else metadata.get("position_columns", ())
    if isinstance(position_cols_raw, str):
        position_cols = [position_cols_raw.strip()]
    else:
        position_cols = [str(column).strip() for column in position_cols_raw if str(column).strip()]
    missing_force = [column for column in force_cols if column not in column_set]
    missing_position = [column for column in position_cols if column not in column_set]
    f_unit = str(force_unit or metadata.get("force_unit") or "").strip().lower().replace(" ", "")
    p_unit = str(position_unit or metadata.get("position_unit") or "").strip().lower().replace(" ", "")
    frame = str(component_frame or metadata.get("component_frame") or "").strip().lower().replace("-", "_").replace(" ", "_")
    projection_axis_text = str(
        projection_axis
        or metadata.get("projection_axis")
        or metadata.get("force_projection_axis")
        or metadata.get("normal_convention")
        or ""
    ).strip()
    projection_axis_key = projection_axis_text.lower().replace("-", "_").replace(" ", "_")
    tool = str(source_tool or metadata.get("source_tool") or "").strip().lower().replace("-", "_").replace(" ", "_")
    sf = float(metadata.get("symmetry_factor", symmetry_factor))
    sign = force_sign_convention if force_sign_convention is not None else metadata.get("force_sign_convention")
    sign_key = None if sign is None else str(sign).strip().lower().replace("-", "_").replace(" ", "_")
    kind = force_kind if force_kind is not None else metadata.get("force_kind")
    kind_key = None if kind is None else str(kind).strip().lower().replace("-", "_").replace(" ", "_")
    qdim_raw = (
        quantity_dimension
        if quantity_dimension is not None
        else metadata.get("quantity_dimension", metadata.get("force_quantity_dimension"))
    )
    qdim_key = (
        None
        if qdim_raw is None
        else str(qdim_raw).strip().lower().replace("-", "_").replace(" ", "_")
    )
    qdim_aliases = {
        "2d": "2d_per_length",
        "2d_per_m": "2d_per_length",
        "2d_per_meter": "2d_per_length",
        "2d_per_metre": "2d_per_length",
        "2d_per_length": "2d_per_length",
        "per_length": "2d_per_length",
        "per_meter": "2d_per_length",
        "per_metre": "2d_per_length",
        "3d": "3d_total",
        "3d_total": "3d_total",
        "axisymmetric_total": "3d_total",
        "total": "3d_total",
        "total_force": "3d_total",
    }
    qdim = qdim_aliases.get(qdim_key, qdim_key)

    def _coordinate_tuple(value):
        if value in (None, ""):
            return None
        if isinstance(value, dict):
            keys = ("x", "y", "z")
            upper_keys = ("X", "Y", "Z")
            if all(key in value for key in keys[:2]):
                coords = [value[key] for key in keys if key in value]
            elif all(key in value for key in upper_keys[:2]):
                coords = [value[key] for key in upper_keys if key in value]
            else:
                raise ValueError("coordinate dictionaries must include x/y or X/Y")
        elif isinstance(value, str):
            coords = [item for item in re.split(r"[,;\s]+", value.strip()) if item]
        else:
            coords = list(value)
        if len(coords) not in (2, 3):
            raise ValueError("coordinates must contain two or three values")
        return tuple(float(coord) for coord in coords)

    def _normalize_analysis_type(value):
        text = str(value).strip().lower().replace("-", "_").replace(" ", "_")
        aliases = {
            "static": "magnetostatic",
            "magnetic_static": "magnetostatic",
            "magneto_static": "magnetostatic",
            "freq": "frequency_response",
            "frequency": "frequency_response",
            "frequency_domain": "frequency_response",
            "ac": "frequency_response",
            "time_transient": "transient",
            "transient_magnetic": "transient",
        }
        return aliases.get(text, text)

    case_id = str(
        metadata.get("case_id", metadata.get("analysis_case_id", metadata.get("design_case_id", "")))
    ).strip()
    study_id = str(metadata.get("study_id", metadata.get("study_name", ""))).strip()
    analysis_type_raw = str(
        metadata.get(
            "analysis_type",
            metadata.get("solver_type", metadata.get("physics_type", metadata.get("study_type", ""))),
        )
    ).strip()
    analysis_type_key = _normalize_analysis_type(analysis_type_raw) if analysis_type_raw else ""
    frequency_source = metadata.get(
        "frequency_hz",
        metadata.get("excitation_frequency_hz", metadata.get("problem_frequency_hz")),
    )
    frequency_hz = None
    if frequency_source not in (None, ""):
        frequency_hz = float(frequency_source)
    operating_point_id = str(
        metadata.get(
            "operating_point_id",
            metadata.get("op_id", metadata.get("operating_condition_id", "")),
        )
    ).strip()
    export_artifact_id = str(
        metadata.get(
            "export_artifact_id",
            metadata.get(
                "artifact_id",
                metadata.get("case_artifact_id", metadata.get("report_artifact_id", "")),
            ),
        )
    ).strip()
    result_set_id = str(metadata.get("result_set_id", "")).strip()
    mesh_id = str(
        metadata.get(
            "mesh_id",
            metadata.get("mesh_revision_id", metadata.get("mesh_artifact_id", "")),
        )
    ).strip()
    solver_run_id = str(
        metadata.get(
            "solver_run_id",
            metadata.get("run_id", metadata.get("jmag_solver_run_id", metadata.get("simulation_id", ""))),
        )
    ).strip()
    result_revision_id = str(
        metadata.get(
            "result_revision_id",
            metadata.get("result_version_id", metadata.get("result_id", "")),
        )
    ).strip()
    solver_setup_artifact_id = str(
        metadata.get(
            "solver_setup_artifact_id",
            metadata.get(
                "solver_settings_artifact_id",
                metadata.get("study_solver_setup_artifact_id", metadata.get("solver_configuration_artifact_id", "")),
            ),
        )
    ).strip()
    material_state_artifact_id = str(
        metadata.get(
            "material_state_artifact_id",
            metadata.get(
                "material_model_artifact_id",
                metadata.get("material_law_artifact_id", metadata.get("bh_curve_artifact_id", "")),
            ),
        )
    ).strip()
    material_state_digest = str(
        metadata.get(
            "material_state_digest",
            metadata.get(
                "material_model_digest",
                metadata.get(
                    "material_law_digest",
                    metadata.get("bh_curve_digest", metadata.get("hysteresis_model_digest", "")),
                ),
            ),
        )
    ).strip()
    excitation_source_artifact_id = str(
        metadata.get(
            "excitation_source_artifact_id",
            metadata.get(
                "current_source_artifact_id",
                metadata.get(
                    "source_current_artifact_id",
                    metadata.get(
                        "current_snapshot_artifact_id",
                        metadata.get("excitation_table_artifact_id", ""),
                    ),
                ),
            ),
        )
    ).strip()
    current_definition_raw = metadata.get(
        "current_definition_method",
        metadata.get(
            "current_basis",
            metadata.get(
                "current_kind",
                metadata.get("current_convention", metadata.get("excitation_current_convention")),
            ),
        ),
    )
    current_definition_method_key = (
        None
        if current_definition_raw in (None, "")
        else str(current_definition_raw).strip().lower().replace("-", "_").replace(" ", "_")
    )
    export_trace_id = str(
        metadata.get(
            "export_trace_id",
            metadata.get(
                "export_macro_trace_id",
                metadata.get("table_export_trace_id", metadata.get("macro_trace_id", "")),
            ),
        )
    ).strip()
    export_command_digest = str(
        metadata.get(
            "export_command_digest",
            metadata.get(
                "export_macro_digest",
                metadata.get("macro_command_digest", metadata.get("export_script_sha256", "")),
            ),
        )
    ).strip()
    export_commands_raw = metadata.get(
        "export_commands",
        metadata.get(
            "export_macro_commands",
            metadata.get("export_command_sequence", metadata.get("macro_command_sequence", ())),
        ),
    )
    if isinstance(export_commands_raw, str):
        command_chunks = export_commands_raw.replace(";", "\n").splitlines()
        export_commands = [command.strip() for command in command_chunks if command.strip()]
    else:
        export_commands = [
            str(command).strip()
            for command in export_commands_raw
            if str(command).strip()
        ] if export_commands_raw else []
    export_output_artifact_id = str(
        metadata.get(
            "export_output_artifact_id",
            metadata.get(
                "table_artifact_id",
                metadata.get(
                    "output_artifact_id",
                    metadata.get("value_table_artifact_id", metadata.get("export_table_artifact_id", "")),
                ),
            ),
        )
    ).strip()
    export_output_digest = str(
        metadata.get(
            "export_output_digest",
            metadata.get(
                "export_output_sha256",
                metadata.get(
                    "table_sha256",
                    metadata.get("output_digest", metadata.get("value_table_digest", "")),
                ),
            ),
        )
    ).strip()
    export_output_path = str(
        metadata.get(
            "export_output_path",
            metadata.get("table_path", metadata.get("output_path", metadata.get("value_table_path", ""))),
        )
    ).strip()
    force_observable_id = str(
        metadata.get(
            "force_observable_id",
            metadata.get(
                "force_report_observable_id",
                metadata.get(
                    "jmag_force_observable_id",
                    metadata.get("observable_id", metadata.get("report_observable_id", "")),
                ),
            ),
        )
    ).strip()
    force_observable_family = str(
        metadata.get(
            "force_observable_family",
            metadata.get(
                "force_report_family",
                metadata.get(
                    "jmag_force_observable_family",
                    metadata.get("observable_family", metadata.get("report_family", "")),
                ),
            ),
        )
    ).strip().lower().replace("-", "_").replace(" ", "_")
    method_raw = (
        force_report_method
        if force_report_method is not None
        else metadata.get(
            "force_report_method",
            metadata.get(
                "force_extraction_method",
                metadata.get("report_method", metadata.get("extraction_method", metadata.get("postprocess_method"))),
            ),
        )
    )
    force_report_method_key = (
        None
        if method_raw is None
        else str(method_raw).strip().lower().replace("-", "_").replace(" ", "_")
    )
    target_region_id = str(
        metadata.get(
            "target_region_id",
            metadata.get("force_region_id", metadata.get("selected_body_id", "")),
        )
    ).strip()
    target_region_name = str(
        metadata.get(
            "target_region_name",
            metadata.get(
                "target_body_name",
                metadata.get("force_region_name", metadata.get("selected_body_name", "")),
            ),
        )
    ).strip()
    target_material = str(
        metadata.get(
            "target_material",
            metadata.get(
                "target_material_name",
                metadata.get("material_name", metadata.get("selected_body_material", "")),
            ),
        )
    ).strip()
    target_region_artifact_id = str(
        metadata.get(
            "target_region_artifact_id",
            metadata.get(
                "target_body_artifact_id",
                metadata.get("region_label_artifact_id", metadata.get("body_label_artifact_id", "")),
            ),
        )
    ).strip()
    target_region_geometry_digest = str(
        metadata.get(
            "target_region_geometry_digest",
            metadata.get(
                "target_body_geometry_digest",
                metadata.get("target_geometry_digest", metadata.get("geometry_digest", "")),
            ),
        )
    ).strip()
    target_region_centroid_xyz_m = _coordinate_tuple(
        metadata.get(
            "target_region_centroid_xyz_m",
            metadata.get(
                "target_body_centroid_xyz_m",
                metadata.get(
                    "target_centroid_xyz_m",
                    metadata.get("target_region_centroid_m", metadata.get("target_centroid_m")),
                ),
            ),
        )
    )
    identity_cols_raw = metadata.get("identity_columns", metadata.get("region_columns", ()))
    if isinstance(identity_cols_raw, str):
        identity_cols = [identity_cols_raw.strip()]
    else:
        identity_cols = [str(column).strip() for column in identity_cols_raw if str(column).strip()]
    missing_identity = [column for column in identity_cols if column not in column_set]
    rows_raw = table_rows
    if rows_raw is None:
        rows_raw = metadata.get("table_rows", metadata.get("sample_rows", ()))
    rows = [dict(row) for row in rows_raw] if rows_raw else []
    table_row_count = len(rows)
    min_rows_raw = metadata.get("min_table_rows", 0) if min_table_rows is None else min_table_rows
    min_rows = int(min_rows_raw or 0)
    if min_rows < 0:
        raise ValueError("min_table_rows must be non-negative")

    def _truthy(value):
        if isinstance(value, str):
            return value.strip().lower() in {"1", "true", "yes", "on", "required"}
        return bool(value)

    unique_required = _truthy(
        metadata.get("require_unique_row_identity", False)
        if require_unique_row_identity is None
        else require_unique_row_identity
    )
    expected_result_set = (
        None if expected_result_set_id is None else str(expected_result_set_id).strip()
    )
    expected_export_artifact = (
        None if expected_export_artifact_id is None else str(expected_export_artifact_id).strip()
    )
    expected_case = None if expected_case_id is None else str(expected_case_id).strip()
    expected_study = None if expected_study_id is None else str(expected_study_id).strip()
    expected_operating_point = (
        None if expected_operating_point_id is None else str(expected_operating_point_id).strip()
    )
    expected_analysis = (
        None if expected_analysis_type is None else _normalize_analysis_type(expected_analysis_type)
    )
    expected_frequency = (
        None if expected_frequency_hz is None else float(expected_frequency_hz)
    )
    expected_target_region = (
        None if expected_target_region_id is None else str(expected_target_region_id).strip()
    )
    expected_target_name = (
        None if expected_target_region_name is None else str(expected_target_region_name).strip()
    )
    expected_target_mat = (
        None if expected_target_material is None else str(expected_target_material).strip()
    )
    expected_target_region_artifact = (
        None
        if expected_target_region_artifact_id is None
        else str(expected_target_region_artifact_id).strip()
    )
    expected_target_geometry_digest = (
        None
        if expected_target_region_geometry_digest is None
        else str(expected_target_region_geometry_digest).strip()
    )
    expected_target_centroid = _coordinate_tuple(expected_target_region_centroid_xyz_m)
    expected_mesh = None if expected_mesh_id is None else str(expected_mesh_id).strip()
    expected_solver_run = (
        None if expected_solver_run_id is None else str(expected_solver_run_id).strip()
    )
    expected_result_revision = (
        None if expected_result_revision_id is None else str(expected_result_revision_id).strip()
    )
    expected_solver_setup_artifact = (
        None
        if expected_solver_setup_artifact_id is None
        else str(expected_solver_setup_artifact_id).strip()
    )
    expected_material_state_artifact = (
        None
        if expected_material_state_artifact_id is None
        else str(expected_material_state_artifact_id).strip()
    )
    expected_material_state_hash = (
        None
        if expected_material_state_digest is None
        else str(expected_material_state_digest).strip()
    )
    expected_excitation_source_artifact = (
        None
        if expected_excitation_source_artifact_id is None
        else str(expected_excitation_source_artifact_id).strip()
    )
    expected_current_definition_key = (
        None
        if expected_current_definition_method is None
        else str(expected_current_definition_method).strip().lower().replace("-", "_").replace(" ", "_")
    )
    expected_export_trace = (
        None if expected_export_trace_id is None else str(expected_export_trace_id).strip()
    )
    expected_export_digest = (
        None
        if expected_export_command_digest is None
        else str(expected_export_command_digest).strip()
    )
    expected_export_output_artifact = (
        None
        if expected_export_output_artifact_id is None
        else str(expected_export_output_artifact_id).strip()
    )
    expected_export_output_hash = (
        None
        if expected_export_output_digest is None
        else str(expected_export_output_digest).strip()
    )
    expected_force_observable = (
        None if expected_force_observable_id is None else str(expected_force_observable_id).strip()
    )
    expected_force_observable_group = (
        None
        if expected_force_observable_family is None
        else str(expected_force_observable_family).strip().lower().replace("-", "_").replace(" ", "_")
    )
    expected_component_frame_key = (
        None
        if expected_component_frame is None
        else str(expected_component_frame).strip().lower().replace("-", "_").replace(" ", "_")
    )
    expected_projection_axis_key = (
        None
        if expected_projection_axis is None
        else str(expected_projection_axis).strip().lower().replace("-", "_").replace(" ", "_")
    )
    expected_sign_key = (
        None
        if expected_force_sign_convention is None
        else str(expected_force_sign_convention).strip().lower().replace("-", "_").replace(" ", "_")
    )
    expected_force_report_method_key = (
        None
        if expected_force_report_method is None
        else str(expected_force_report_method).strip().lower().replace("-", "_").replace(" ", "_")
    )
    case_required = expected_case is not None
    study_required = expected_study is not None
    operating_point_required = expected_operating_point is not None
    analysis_required = expected_analysis is not None
    frequency_required = expected_frequency is not None
    export_artifact_required = expected_export_artifact is not None
    result_set_required = expected_result_set is not None
    target_region_required = expected_target_region is not None
    target_region_name_required = expected_target_name is not None
    target_material_required = expected_target_mat is not None
    target_region_artifact_required = expected_target_region_artifact is not None
    target_region_geometry_digest_required = expected_target_geometry_digest is not None
    target_region_centroid_required = expected_target_centroid is not None
    mesh_required = expected_mesh is not None
    solver_run_required = expected_solver_run is not None
    result_revision_required = expected_result_revision is not None
    solver_setup_artifact_required = expected_solver_setup_artifact is not None
    material_state_artifact_required = expected_material_state_artifact is not None
    material_state_digest_required = expected_material_state_hash is not None
    excitation_source_artifact_required = expected_excitation_source_artifact is not None
    current_definition_method_required = expected_current_definition_key is not None
    export_trace_required = expected_export_trace is not None or _truthy(require_export_command_trace)
    export_digest_required = expected_export_digest is not None or _truthy(require_export_command_trace)
    export_commands_required = _truthy(require_export_command_trace)
    export_output_required = (
        expected_export_output_artifact is not None
        or expected_export_output_hash is not None
        or _truthy(require_export_output_artifact)
    )
    export_output_digest_required = expected_export_output_hash is not None or _truthy(
        require_export_output_artifact
    )
    export_output_path_required = _truthy(require_export_output_artifact)
    force_observable_required = expected_force_observable is not None
    force_observable_family_required = expected_force_observable_group is not None
    force_report_method_required = expected_force_report_method_key is not None
    target_region_for_rows = expected_target_region or target_region_id
    target_region_name_for_rows = expected_target_name or target_region_name
    target_material_for_rows = expected_target_mat or target_material
    target_region_artifact_for_rows = expected_target_region_artifact or target_region_artifact_id
    case_for_rows = expected_case or case_id
    operating_point_for_rows = expected_operating_point or operating_point_id
    analysis_for_rows = expected_analysis or analysis_type_key
    frequency_for_rows = expected_frequency if expected_frequency is not None else frequency_hz
    mesh_for_rows = expected_mesh or mesh_id
    solver_run_for_rows = expected_solver_run or solver_run_id
    result_revision_for_rows = expected_result_revision or result_revision_id
    solver_setup_artifact_for_rows = expected_solver_setup_artifact or solver_setup_artifact_id
    material_state_artifact_for_rows = expected_material_state_artifact or material_state_artifact_id
    material_state_digest_for_rows = expected_material_state_hash or material_state_digest
    excitation_source_artifact_for_rows = (
        expected_excitation_source_artifact or excitation_source_artifact_id
    )
    current_definition_method_for_rows = (
        expected_current_definition_key or current_definition_method_key
    )
    centroid_tolerance_m = 1.0e-9
    target_region_centroid_matches_expected = True
    if expected_target_centroid is not None:
        target_region_centroid_matches_expected = (
            target_region_centroid_xyz_m is not None
            and len(target_region_centroid_xyz_m) == len(expected_target_centroid)
            and all(
                abs(actual - expected) <= centroid_tolerance_m
                for actual, expected in zip(target_region_centroid_xyz_m, expected_target_centroid)
            )
        )
    case_row_columns = ("CaseId", "case_id", "case", "CaseID")
    operating_point_row_columns = (
        "OperatingPointId",
        "operating_point_id",
        "OpId",
        "op_id",
        "OperatingConditionId",
    )
    analysis_row_columns = ("AnalysisType", "analysis_type", "SolverType", "solver_type", "StudyType")
    frequency_row_columns = ("FrequencyHz", "frequency_hz", "Freq_Hz", "freq_hz")
    mesh_row_columns = ("MeshId", "mesh_id", "MeshRevisionId", "mesh_revision_id", "MeshArtifactId")
    solver_run_row_columns = ("SolverRunId", "solver_run_id", "RunId", "run_id", "SimulationId")
    result_revision_row_columns = (
        "ResultRevisionId",
        "result_revision_id",
        "ResultVersionId",
        "result_version_id",
        "ResultId",
        "result_id",
    )
    solver_setup_artifact_row_columns = (
        "SolverSetupArtifactId",
        "solver_setup_artifact_id",
        "SolverSettingsArtifactId",
        "solver_settings_artifact_id",
        "SolverConfigurationArtifactId",
        "solver_configuration_artifact_id",
    )
    material_state_artifact_row_columns = (
        "MaterialStateArtifactId",
        "material_state_artifact_id",
        "MaterialModelArtifactId",
        "material_model_artifact_id",
        "MaterialLawArtifactId",
        "material_law_artifact_id",
        "BhCurveArtifactId",
        "bh_curve_artifact_id",
    )
    material_state_digest_row_columns = (
        "MaterialStateDigest",
        "material_state_digest",
        "MaterialModelDigest",
        "material_model_digest",
        "MaterialLawDigest",
        "material_law_digest",
        "BhCurveDigest",
        "bh_curve_digest",
        "HysteresisModelDigest",
        "hysteresis_model_digest",
    )
    excitation_source_artifact_row_columns = (
        "ExcitationSourceArtifactId",
        "excitation_source_artifact_id",
        "CurrentSourceArtifactId",
        "current_source_artifact_id",
        "SourceCurrentArtifactId",
        "source_current_artifact_id",
        "CurrentSnapshotArtifactId",
        "current_snapshot_artifact_id",
        "ExcitationTableArtifactId",
        "excitation_table_artifact_id",
    )
    current_definition_method_row_columns = (
        "CurrentDefinitionMethod",
        "current_definition_method",
        "CurrentBasis",
        "current_basis",
        "CurrentKind",
        "current_kind",
        "CurrentConvention",
        "current_convention",
        "ExcitationCurrentConvention",
        "excitation_current_convention",
    )
    target_region_name_row_columns = (
        "TargetRegionName",
        "target_region_name",
        "TargetBodyName",
        "target_body_name",
        "BodyName",
        "body_name",
        "RegionName",
        "region_name",
    )
    target_material_row_columns = (
        "TargetMaterial",
        "target_material",
        "TargetMaterialName",
        "target_material_name",
        "MaterialName",
        "material_name",
        "BodyMaterial",
        "body_material",
    )
    target_region_artifact_row_columns = (
        "TargetRegionArtifactId",
        "target_region_artifact_id",
        "TargetBodyArtifactId",
        "target_body_artifact_id",
        "RegionArtifactId",
        "region_artifact_id",
        "BodyArtifactId",
        "body_artifact_id",
    )
    case_column_present = any(column in column_set for column in case_row_columns)
    operating_point_column_present = any(
        column in column_set for column in operating_point_row_columns
    )
    analysis_column_present = any(column in column_set for column in analysis_row_columns)
    frequency_column_present = any(column in column_set for column in frequency_row_columns)
    mesh_column_present = any(column in column_set for column in mesh_row_columns)
    solver_run_column_present = any(column in column_set for column in solver_run_row_columns)
    result_revision_column_present = any(
        column in column_set for column in result_revision_row_columns
    )
    solver_setup_artifact_column_present = any(
        column in column_set for column in solver_setup_artifact_row_columns
    )
    material_state_artifact_column_present = any(
        column in column_set for column in material_state_artifact_row_columns
    )
    material_state_digest_column_present = any(
        column in column_set for column in material_state_digest_row_columns
    )
    excitation_source_artifact_column_present = any(
        column in column_set for column in excitation_source_artifact_row_columns
    )
    current_definition_method_column_present = any(
        column in column_set for column in current_definition_method_row_columns
    )
    target_region_name_column_present = any(
        column in column_set for column in target_region_name_row_columns
    )
    target_material_column_present = any(
        column in column_set for column in target_material_row_columns
    )
    target_region_artifact_column_present = any(
        column in column_set for column in target_region_artifact_row_columns
    )
    missing_row_identity = []
    row_identity_mismatches = []
    row_case_id_mismatches = []
    missing_row_case_id = []
    row_operating_point_id_mismatches = []
    missing_row_operating_point_id = []
    row_analysis_type_mismatches = []
    missing_row_analysis_type = []
    row_frequency_hz_mismatches = []
    missing_row_frequency_hz = []
    row_mesh_id_mismatches = []
    missing_row_mesh_id = []
    row_solver_run_id_mismatches = []
    missing_row_solver_run_id = []
    row_result_revision_id_mismatches = []
    missing_row_result_revision_id = []
    row_solver_setup_artifact_id_mismatches = []
    missing_row_solver_setup_artifact_id = []
    row_material_state_artifact_id_mismatches = []
    missing_row_material_state_artifact_id = []
    row_material_state_digest_mismatches = []
    missing_row_material_state_digest = []
    row_excitation_source_artifact_id_mismatches = []
    missing_row_excitation_source_artifact_id = []
    row_current_definition_method_mismatches = []
    missing_row_current_definition_method = []
    row_target_region_name_mismatches = []
    missing_row_target_region_name = []
    row_target_material_mismatches = []
    missing_row_target_material = []
    row_target_region_artifact_id_mismatches = []
    missing_row_target_region_artifact_id = []
    row_identity_values = []
    duplicate_row_identity_values = []
    seen_identity_rows = {}
    frequency_tol = 1.0e-12

    def _first_row_value(row, names):
        for name in names:
            value = row.get(name)
            if value not in (None, ""):
                return str(value).strip()
        return ""

    for index, row in enumerate(rows):
        values = [str(row.get(column, "")).strip() for column in identity_cols]
        row_identity_values.append({"row": index, "values": values})
        if identity_cols and not all(values):
            missing_row_identity.append(index)
        if target_region_for_rows and identity_cols and target_region_for_rows not in values:
            row_identity_mismatches.append(index)
        if unique_required and identity_cols and all(values):
            identity_key = tuple(values)
            if identity_key in seen_identity_rows:
                duplicate_row_identity_values.append(
                    {
                        "first_row": seen_identity_rows[identity_key],
                        "row": index,
                        "values": values,
                    }
                )
            else:
                seen_identity_rows[identity_key] = index
        if case_for_rows and case_column_present:
            row_case = _first_row_value(row, case_row_columns)
            if not row_case:
                missing_row_case_id.append(index)
            elif row_case != case_for_rows:
                row_case_id_mismatches.append(index)
        if operating_point_for_rows and operating_point_column_present:
            row_op = _first_row_value(row, operating_point_row_columns)
            if not row_op:
                missing_row_operating_point_id.append(index)
            elif row_op != operating_point_for_rows:
                row_operating_point_id_mismatches.append(index)
        if analysis_for_rows and analysis_column_present:
            row_analysis = _first_row_value(row, analysis_row_columns)
            if not row_analysis:
                missing_row_analysis_type.append(index)
            elif _normalize_analysis_type(row_analysis) != analysis_for_rows:
                row_analysis_type_mismatches.append(index)
        if frequency_for_rows is not None and frequency_column_present:
            row_frequency = _first_row_value(row, frequency_row_columns)
            if not row_frequency:
                missing_row_frequency_hz.append(index)
            else:
                try:
                    row_frequency_value = float(row_frequency)
                except (TypeError, ValueError):
                    row_frequency_value = math.nan
                if not math.isfinite(row_frequency_value) or abs(row_frequency_value - frequency_for_rows) > frequency_tol:
                    row_frequency_hz_mismatches.append(index)
        if mesh_for_rows and mesh_column_present:
            row_mesh = _first_row_value(row, mesh_row_columns)
            if not row_mesh:
                missing_row_mesh_id.append(index)
            elif row_mesh != mesh_for_rows:
                row_mesh_id_mismatches.append(index)
        if solver_run_for_rows and solver_run_column_present:
            row_solver_run = _first_row_value(row, solver_run_row_columns)
            if not row_solver_run:
                missing_row_solver_run_id.append(index)
            elif row_solver_run != solver_run_for_rows:
                row_solver_run_id_mismatches.append(index)
        if result_revision_for_rows and result_revision_column_present:
            row_result_revision = _first_row_value(row, result_revision_row_columns)
            if not row_result_revision:
                missing_row_result_revision_id.append(index)
            elif row_result_revision != result_revision_for_rows:
                row_result_revision_id_mismatches.append(index)
        if solver_setup_artifact_for_rows and solver_setup_artifact_column_present:
            row_solver_setup_artifact = _first_row_value(row, solver_setup_artifact_row_columns)
            if not row_solver_setup_artifact:
                missing_row_solver_setup_artifact_id.append(index)
            elif row_solver_setup_artifact != solver_setup_artifact_for_rows:
                row_solver_setup_artifact_id_mismatches.append(index)
        if material_state_artifact_for_rows and material_state_artifact_column_present:
            row_material_state_artifact = _first_row_value(row, material_state_artifact_row_columns)
            if not row_material_state_artifact:
                missing_row_material_state_artifact_id.append(index)
            elif row_material_state_artifact != material_state_artifact_for_rows:
                row_material_state_artifact_id_mismatches.append(index)
        if material_state_digest_for_rows and material_state_digest_column_present:
            row_material_state_hash = _first_row_value(row, material_state_digest_row_columns)
            if not row_material_state_hash:
                missing_row_material_state_digest.append(index)
            elif row_material_state_hash != material_state_digest_for_rows:
                row_material_state_digest_mismatches.append(index)
        if excitation_source_artifact_for_rows and excitation_source_artifact_column_present:
            row_excitation_source_artifact = _first_row_value(row, excitation_source_artifact_row_columns)
            if not row_excitation_source_artifact:
                missing_row_excitation_source_artifact_id.append(index)
            elif row_excitation_source_artifact != excitation_source_artifact_for_rows:
                row_excitation_source_artifact_id_mismatches.append(index)
        if current_definition_method_for_rows and current_definition_method_column_present:
            row_current_definition = _first_row_value(row, current_definition_method_row_columns)
            if not row_current_definition:
                missing_row_current_definition_method.append(index)
            elif row_current_definition.strip().lower().replace("-", "_").replace(" ", "_") != current_definition_method_for_rows:
                row_current_definition_method_mismatches.append(index)
        if target_region_name_for_rows and target_region_name_column_present:
            row_target_name = _first_row_value(row, target_region_name_row_columns)
            if not row_target_name:
                missing_row_target_region_name.append(index)
            elif row_target_name != target_region_name_for_rows:
                row_target_region_name_mismatches.append(index)
        if target_material_for_rows and target_material_column_present:
            row_target_mat = _first_row_value(row, target_material_row_columns)
            if not row_target_mat:
                missing_row_target_material.append(index)
            elif row_target_mat != target_material_for_rows:
                row_target_material_mismatches.append(index)
        if target_region_artifact_for_rows and target_region_artifact_column_present:
            row_target_artifact = _first_row_value(row, target_region_artifact_row_columns)
            if not row_target_artifact:
                missing_row_target_region_artifact_id.append(index)
            elif row_target_artifact != target_region_artifact_for_rows:
                row_target_region_artifact_id_mismatches.append(index)
    force_units = {"n", "n/m", "n_per_m", "nperm", "n_per_meter", "n_per_metre"}
    force_units_by_dimension = {
        "2d_per_length": {"n/m", "n_per_m", "nperm", "n_per_meter", "n_per_metre"},
        "3d_total": {"n"},
    }
    position_units = {"m", "mm", "deg", "degree", "rad", "s", "sec"}
    frames = {"global_xy", "global_xyz", "cylindrical_rt", "local_rt", "as_exported"}
    signs = {
        "positive_attraction",
        "positive_repulsion",
        "positive_motoring",
        "positive_generator",
        "positive_x",
        "positive_y",
        "positive_z",
        "as_exported",
    }
    kinds = {"block_integral", "maxwell_stress", "virtual_work", "nodal_force", "line_probe", "as_exported"}
    def _column_axis_flags(force_columns):
        flags = {"x": False, "y": False, "z": False, "radial": False, "tangential": False}
        for column in force_columns:
            token = str(column).strip().lower().replace("-", "_").replace(" ", "_")
            compact = token.replace("_", "")
            flags["x"] = flags["x"] or token.startswith("fx") or compact.startswith("forcex")
            flags["y"] = flags["y"] or token.startswith("fy") or compact.startswith("forcey")
            flags["z"] = flags["z"] or token.startswith("fz") or compact.startswith("forcez")
            flags["radial"] = flags["radial"] or token.startswith("fr") or "radial" in token or compact.startswith("forcer")
            flags["tangential"] = (
                flags["tangential"]
                or token.startswith("ft")
                or "tangential" in token
                or token.startswith("f_theta")
                or token.startswith("fth")
                or compact.startswith("forcet")
            )
        return flags

    axis_flags = _column_axis_flags(force_cols)
    expected_axes_by_frame = {
        "global_xy": ("x", "y"),
        "global_xyz": ("x", "y", "z"),
        "cylindrical_rt": ("radial", "tangential"),
        "local_rt": ("radial", "tangential"),
    }
    expected_axes = expected_axes_by_frame.get(frame, ())
    force_columns_match_frame = not expected_axes or all(axis_flags[axis] for axis in expected_axes)
    projection_axis_required = sign_key not in (None, "as_exported")
    projection_axis_has_physics = (
        not projection_axis_required
        or any(
            token in projection_axis_key
            for token in ("gap", "normal", "radial", "target", "source", "rotor", "stator", "x", "y", "z")
        )
    )
    export_command_blob = "\n".join(export_commands).lower().replace("_", "").replace(" ", "")
    export_commands_include_table_export = (
        not export_commands_required
        or "writeallcasetable" in export_command_blob
        or ("export" in export_command_blob and "table" in export_command_blob)
    )
    export_commands_reference_force_report = (
        not export_commands_required
        or "force" in export_command_blob
        or "maxwell" in export_command_blob
        or "stress" in export_command_blob
    )
    checks = {
        "required_columns_present": not missing_required,
        "force_columns_recorded": bool(force_cols),
        "force_columns_present": bool(force_cols) and not missing_force,
        "force_columns_match_component_frame": force_columns_match_frame,
        "position_columns_present": not missing_position,
        "force_unit_valid": f_unit in force_units,
        "position_unit_valid": not p_unit or p_unit in position_units,
        "component_frame_valid": not frame or frame in frames,
        "component_frame_recorded_when_expected": expected_component_frame_key is None
        or bool(frame),
        "expected_component_frame_matches": expected_component_frame_key is None
        or frame == expected_component_frame_key,
        "force_projection_axis_recorded": not projection_axis_required or bool(projection_axis_text),
        "force_projection_axis_descriptive": projection_axis_has_physics,
        "projection_axis_recorded_when_expected": expected_projection_axis_key is None
        or bool(projection_axis_text),
        "expected_projection_axis_matches": expected_projection_axis_key is None
        or projection_axis_key == expected_projection_axis_key,
        "source_tool_is_jmag": tool in {"jmag", "jmag_designer", "jmagdesigner"},
        "symmetry_factor_positive": sf > 0.0,
        "force_sign_convention_valid": sign_key is None or sign_key in signs,
        "force_sign_convention_recorded_when_expected": expected_sign_key is None
        or bool(sign_key),
        "expected_force_sign_convention_matches": expected_sign_key is None
        or sign_key == expected_sign_key,
        "force_report_method_recorded_when_expected": not force_report_method_required
        or bool(force_report_method_key),
        "expected_force_report_method_matches": expected_force_report_method_key is None
        or force_report_method_key == expected_force_report_method_key,
        "force_kind_valid": kind_key is None or kind_key in kinds,
        "quantity_dimension_valid": qdim is None or qdim in force_units_by_dimension,
        "force_unit_matches_quantity_dimension_when_present": qdim is None
        or f_unit in force_units_by_dimension.get(qdim, set()),
        "case_id_recorded": not case_required or bool(case_id),
        "expected_case_id_matches": expected_case is None or case_id == expected_case,
        "study_id_recorded": not study_required or bool(study_id),
        "expected_study_id_matches": expected_study is None or study_id == expected_study,
        "operating_point_id_recorded": not operating_point_required or bool(operating_point_id),
        "expected_operating_point_id_matches": expected_operating_point is None
        or operating_point_id == expected_operating_point,
        "analysis_type_recorded": not analysis_required or bool(analysis_type_key),
        "expected_analysis_type_matches": expected_analysis is None
        or analysis_type_key == expected_analysis,
        "frequency_hz_recorded": not frequency_required or frequency_hz is not None,
        "expected_frequency_hz_matches": expected_frequency is None
        or (frequency_hz is not None and abs(frequency_hz - expected_frequency) <= frequency_tol),
        "export_artifact_id_recorded": not export_artifact_required or bool(export_artifact_id),
        "expected_export_artifact_id_matches": expected_export_artifact is None
        or export_artifact_id == expected_export_artifact,
        "result_set_id_recorded": not result_set_required or bool(result_set_id),
        "expected_result_set_id_matches": expected_result_set is None
        or result_set_id == expected_result_set,
        "target_region_id_recorded": not target_region_required or bool(target_region_id),
        "expected_target_region_id_matches": expected_target_region is None
        or target_region_id == expected_target_region,
        "target_region_name_recorded": not target_region_name_required or bool(target_region_name),
        "expected_target_region_name_matches": expected_target_name is None
        or target_region_name == expected_target_name,
        "target_material_recorded": not target_material_required or bool(target_material),
        "expected_target_material_matches": expected_target_mat is None
        or target_material == expected_target_mat,
        "target_region_artifact_id_recorded": not target_region_artifact_required
        or bool(target_region_artifact_id),
        "expected_target_region_artifact_id_matches": expected_target_region_artifact is None
        or target_region_artifact_id == expected_target_region_artifact,
        "target_region_geometry_digest_recorded": not target_region_geometry_digest_required
        or bool(target_region_geometry_digest),
        "expected_target_region_geometry_digest_matches": expected_target_geometry_digest is None
        or target_region_geometry_digest == expected_target_geometry_digest,
        "target_region_centroid_xyz_recorded_when_expected": not target_region_centroid_required
        or target_region_centroid_xyz_m is not None,
        "expected_target_region_centroid_xyz_matches": target_region_centroid_matches_expected,
        "mesh_id_recorded": not mesh_required or bool(mesh_id),
        "expected_mesh_id_matches": expected_mesh is None or mesh_id == expected_mesh,
        "solver_run_id_recorded": not solver_run_required or bool(solver_run_id),
        "expected_solver_run_id_matches": expected_solver_run is None
        or solver_run_id == expected_solver_run,
        "result_revision_id_recorded": not result_revision_required or bool(result_revision_id),
        "expected_result_revision_id_matches": expected_result_revision is None
        or result_revision_id == expected_result_revision,
        "solver_setup_artifact_id_recorded": not solver_setup_artifact_required
        or bool(solver_setup_artifact_id),
        "expected_solver_setup_artifact_id_matches": expected_solver_setup_artifact is None
        or solver_setup_artifact_id == expected_solver_setup_artifact,
        "material_state_artifact_id_recorded": not material_state_artifact_required
        or bool(material_state_artifact_id),
        "expected_material_state_artifact_id_matches": expected_material_state_artifact is None
        or material_state_artifact_id == expected_material_state_artifact,
        "material_state_digest_recorded": not material_state_digest_required
        or bool(material_state_digest),
        "expected_material_state_digest_matches": expected_material_state_hash is None
        or material_state_digest == expected_material_state_hash,
        "excitation_source_artifact_id_recorded": not excitation_source_artifact_required
        or bool(excitation_source_artifact_id),
        "expected_excitation_source_artifact_id_matches": expected_excitation_source_artifact is None
        or excitation_source_artifact_id == expected_excitation_source_artifact,
        "current_definition_method_recorded": not current_definition_method_required
        or bool(current_definition_method_key),
        "expected_current_definition_method_matches": expected_current_definition_key is None
        or current_definition_method_key == expected_current_definition_key,
        "export_trace_id_recorded": not export_trace_required or bool(export_trace_id),
        "expected_export_trace_id_matches": expected_export_trace is None
        or export_trace_id == expected_export_trace,
        "export_command_digest_recorded": not export_digest_required or bool(export_command_digest),
        "expected_export_command_digest_matches": expected_export_digest is None
        or export_command_digest == expected_export_digest,
        "export_commands_recorded": not export_commands_required or bool(export_commands),
        "export_commands_include_table_export": export_commands_include_table_export,
        "export_commands_reference_force_report": export_commands_reference_force_report,
        "export_output_artifact_id_recorded": not export_output_required
        or bool(export_output_artifact_id),
        "expected_export_output_artifact_id_matches": expected_export_output_artifact is None
        or export_output_artifact_id == expected_export_output_artifact,
        "export_output_digest_recorded": not export_output_digest_required
        or bool(export_output_digest),
        "expected_export_output_digest_matches": expected_export_output_hash is None
        or export_output_digest == expected_export_output_hash,
        "export_output_path_recorded": not export_output_path_required
        or bool(export_output_path),
        "force_observable_id_recorded": not force_observable_required or bool(force_observable_id),
        "expected_force_observable_id_matches": expected_force_observable is None
        or force_observable_id == expected_force_observable,
        "force_observable_family_recorded": not force_observable_family_required
        or bool(force_observable_family),
        "expected_force_observable_family_matches": expected_force_observable_group is None
        or force_observable_family == expected_force_observable_group,
        "identity_columns_present": not missing_identity,
        "row_identity_columns_recorded_for_uniqueness": not unique_required or bool(identity_cols),
        "table_rows_meet_minimum": table_row_count >= min_rows,
        "row_identity_columns_populated": not missing_row_identity,
        "row_identity_unique_when_required": not duplicate_row_identity_values,
        "row_identity_matches_target_region": not row_identity_mismatches,
        "row_case_id_matches_package": not missing_row_case_id and not row_case_id_mismatches,
        "row_operating_point_id_matches_package": (
            not missing_row_operating_point_id and not row_operating_point_id_mismatches
        ),
        "row_analysis_type_matches_package": (
            not missing_row_analysis_type and not row_analysis_type_mismatches
        ),
        "row_frequency_hz_matches_package": (
            not missing_row_frequency_hz and not row_frequency_hz_mismatches
        ),
        "row_mesh_id_matches_package": not missing_row_mesh_id and not row_mesh_id_mismatches,
        "row_solver_run_id_matches_package": (
            not missing_row_solver_run_id and not row_solver_run_id_mismatches
        ),
        "row_result_revision_id_matches_package": (
            not missing_row_result_revision_id and not row_result_revision_id_mismatches
        ),
        "row_solver_setup_artifact_id_matches_package": (
            not missing_row_solver_setup_artifact_id and not row_solver_setup_artifact_id_mismatches
        ),
        "row_material_state_artifact_id_matches_package": (
            not missing_row_material_state_artifact_id and not row_material_state_artifact_id_mismatches
        ),
        "row_material_state_digest_matches_package": (
            not missing_row_material_state_digest and not row_material_state_digest_mismatches
        ),
        "row_excitation_source_artifact_id_matches_package": (
            not missing_row_excitation_source_artifact_id
            and not row_excitation_source_artifact_id_mismatches
        ),
        "row_current_definition_method_matches_package": (
            not missing_row_current_definition_method
            and not row_current_definition_method_mismatches
        ),
        "row_target_region_name_matches_package": (
            not missing_row_target_region_name and not row_target_region_name_mismatches
        ),
        "row_target_material_matches_package": (
            not missing_row_target_material and not row_target_material_mismatches
        ),
        "row_target_region_artifact_id_matches_package": (
            not missing_row_target_region_artifact_id
            and not row_target_region_artifact_id_mismatches
        ),
    }
    return {
        "policy": "jmag_force_table_metadata_gate",
        "columns": columns,
        "required_columns": required,
        "missing_required_columns": missing_required,
        "force_columns": force_cols,
        "position_columns": position_cols,
        "missing_force_columns": missing_force,
        "missing_position_columns": missing_position,
        "force_unit": f_unit or None,
        "position_unit": p_unit or None,
        "component_frame": frame or None,
        "projection_axis": projection_axis_text or None,
        "force_column_axes_present": axis_flags,
        "expected_force_axes_for_frame": list(expected_axes),
        "source_tool": tool or None,
        "symmetry_factor": sf,
        "force_sign_convention": sign_key,
        "force_kind": kind_key,
        "quantity_dimension": qdim,
        "case_id": case_id or None,
        "study_id": study_id or None,
        "analysis_type": analysis_type_key or None,
        "raw_analysis_type": analysis_type_raw or None,
        "frequency_hz": frequency_hz,
        "operating_point_id": operating_point_id or None,
        "export_artifact_id": export_artifact_id or None,
        "result_set_id": result_set_id or None,
        "mesh_id": mesh_id or None,
        "solver_run_id": solver_run_id or None,
        "result_revision_id": result_revision_id or None,
        "solver_setup_artifact_id": solver_setup_artifact_id or None,
        "material_state_artifact_id": material_state_artifact_id or None,
        "material_state_digest": material_state_digest or None,
        "excitation_source_artifact_id": excitation_source_artifact_id or None,
        "current_definition_method": current_definition_method_key,
        "export_trace_id": export_trace_id or None,
        "export_command_digest": export_command_digest or None,
        "export_commands": export_commands,
        "export_output_artifact_id": export_output_artifact_id or None,
        "export_output_digest": export_output_digest or None,
        "export_output_path": export_output_path or None,
        "force_observable_id": force_observable_id or None,
        "force_observable_family": force_observable_family or None,
        "force_report_method": force_report_method_key,
        "target_region_id": target_region_id or None,
        "target_region_name": target_region_name or None,
        "target_material": target_material or None,
        "target_region_artifact_id": target_region_artifact_id or None,
        "target_region_geometry_digest": target_region_geometry_digest or None,
        "target_region_centroid_xyz_m": (
            list(target_region_centroid_xyz_m)
            if target_region_centroid_xyz_m is not None
            else None
        ),
        "expected_case_id": expected_case,
        "expected_study_id": expected_study,
        "expected_operating_point_id": expected_operating_point,
        "expected_analysis_type": expected_analysis,
        "expected_frequency_hz": expected_frequency,
        "expected_export_artifact_id": expected_export_artifact,
        "expected_result_set_id": expected_result_set,
        "expected_mesh_id": expected_mesh,
        "expected_solver_run_id": expected_solver_run,
        "expected_result_revision_id": expected_result_revision,
        "expected_solver_setup_artifact_id": expected_solver_setup_artifact,
        "expected_material_state_artifact_id": expected_material_state_artifact,
        "expected_material_state_digest": expected_material_state_hash,
        "expected_excitation_source_artifact_id": expected_excitation_source_artifact,
        "expected_current_definition_method": expected_current_definition_key,
        "expected_export_trace_id": expected_export_trace,
        "expected_export_command_digest": expected_export_digest,
        "expected_export_output_artifact_id": expected_export_output_artifact,
        "expected_export_output_digest": expected_export_output_hash,
        "expected_force_observable_id": expected_force_observable,
        "expected_force_observable_family": expected_force_observable_group,
        "expected_component_frame": expected_component_frame_key,
        "expected_projection_axis": expected_projection_axis_key,
        "expected_force_sign_convention": expected_sign_key,
        "expected_force_report_method": expected_force_report_method_key,
        "require_export_command_trace": _truthy(require_export_command_trace),
        "require_export_output_artifact": _truthy(require_export_output_artifact),
        "expected_target_region_id": expected_target_region,
        "expected_target_region_name": expected_target_name,
        "expected_target_material": expected_target_mat,
        "expected_target_region_artifact_id": expected_target_region_artifact,
        "expected_target_region_geometry_digest": expected_target_geometry_digest,
        "expected_target_region_centroid_xyz_m": (
            list(expected_target_centroid) if expected_target_centroid is not None else None
        ),
        "identity_columns": identity_cols,
        "missing_identity_columns": missing_identity,
        "table_row_count": table_row_count,
        "min_table_rows": min_rows,
        "require_unique_row_identity": unique_required,
        "row_identity_values": row_identity_values,
        "missing_row_identity_rows": missing_row_identity,
        "duplicate_row_identity_values": duplicate_row_identity_values,
        "row_identity_mismatch_rows": row_identity_mismatches,
        "missing_row_case_id_rows": missing_row_case_id,
        "row_case_id_mismatch_rows": row_case_id_mismatches,
        "missing_row_operating_point_id_rows": missing_row_operating_point_id,
        "row_operating_point_id_mismatch_rows": row_operating_point_id_mismatches,
        "missing_row_analysis_type_rows": missing_row_analysis_type,
        "row_analysis_type_mismatch_rows": row_analysis_type_mismatches,
        "missing_row_frequency_hz_rows": missing_row_frequency_hz,
        "row_frequency_hz_mismatch_rows": row_frequency_hz_mismatches,
        "missing_row_mesh_id_rows": missing_row_mesh_id,
        "row_mesh_id_mismatch_rows": row_mesh_id_mismatches,
        "missing_row_solver_run_id_rows": missing_row_solver_run_id,
        "row_solver_run_id_mismatch_rows": row_solver_run_id_mismatches,
        "missing_row_result_revision_id_rows": missing_row_result_revision_id,
        "row_result_revision_id_mismatch_rows": row_result_revision_id_mismatches,
        "missing_row_solver_setup_artifact_id_rows": missing_row_solver_setup_artifact_id,
        "row_solver_setup_artifact_id_mismatch_rows": row_solver_setup_artifact_id_mismatches,
        "missing_row_material_state_artifact_id_rows": missing_row_material_state_artifact_id,
        "row_material_state_artifact_id_mismatch_rows": row_material_state_artifact_id_mismatches,
        "missing_row_material_state_digest_rows": missing_row_material_state_digest,
        "row_material_state_digest_mismatch_rows": row_material_state_digest_mismatches,
        "missing_row_excitation_source_artifact_id_rows": missing_row_excitation_source_artifact_id,
        "row_excitation_source_artifact_id_mismatch_rows": row_excitation_source_artifact_id_mismatches,
        "missing_row_current_definition_method_rows": missing_row_current_definition_method,
        "row_current_definition_method_mismatch_rows": row_current_definition_method_mismatches,
        "missing_row_target_region_name_rows": missing_row_target_region_name,
        "row_target_region_name_mismatch_rows": row_target_region_name_mismatches,
        "missing_row_target_material_rows": missing_row_target_material,
        "row_target_material_mismatch_rows": row_target_material_mismatches,
        "missing_row_target_region_artifact_id_rows": missing_row_target_region_artifact_id,
        "row_target_region_artifact_id_mismatch_rows": row_target_region_artifact_id_mismatches,
        "checks": checks,
        "status": "ok" if all(checks.values()) else "needs_attention",
        "version_note": (
            "Use this before JMAG force/NVH/contact-force table parsing so "
            "columns, units, force dimensionality, frame, symmetry scaling, "
            "projection axis, sign convention, mesh/run/result revision identity, "
            "solver setup, material-state, and excitation-source identity, "
            "target body/material label identity, target geometry digest/centroid identity, export output artifact identity, "
            "force observable identity, force report method, expected frame/projection/sign convention, analysis type, and excitation "
            "frequency are explicit before comparing Maxwell-stress, virtual-work, or radia-ngsolve force results."
        ),
    }


def jmag_airgap_flux_sample_metadata_gate(
    rows,
    *,
    expected_result_set_id=None,
    expected_export_artifact_id=None,
    expected_field_probe_id=None,
    expected_field_probe_method=None,
    expected_field_probe_output_artifact_id=None,
    expected_field_probe_output_digest=None,
    expected_sample_grid_id=None,
    expected_sample_grid_digest=None,
    expected_sample_count=None,
    expected_angle_unit="deg",
    expected_angle_basis=None,
    expected_component_frame="cylindrical_rt",
    expected_torque_sign_convention=None,
    require_field_probe_output_artifact=False,
    min_rows=3,
):
    """Check air-gap flux sample rows before torque use.

    The gate is intentionally solver-independent: it does not trust Br/Bt
    values until the row package also carries the result/export/probe/output
    identity that produced the sampled field table.
    """

    data = list(rows)
    if not data:
        raise ValueError("at least one air-gap flux sample row is required")

    def _norm(value):
        return str(value).strip().lower().replace("-", "_").replace(" ", "_")

    def _first(row, names):
        for name in names:
            if name in row and row[name] not in (None, ""):
                return row[name]
        return None

    def _float_or_none(value):
        if value in (None, ""):
            return None
        return float(value)

    def _finite(value):
        return value is not None and math.isfinite(float(value))

    def _positive(value):
        return value is not None and math.isfinite(float(value)) and float(value) > 0.0

    def _truthy(value):
        if isinstance(value, bool):
            return value
        return _norm(value) in {"1", "true", "yes", "y", "required", "require"}

    expected_result = None if expected_result_set_id is None else str(expected_result_set_id).strip()
    expected_export = None if expected_export_artifact_id is None else str(expected_export_artifact_id).strip()
    expected_probe = None if expected_field_probe_id is None else str(expected_field_probe_id).strip()
    expected_method = (
        None if expected_field_probe_method is None else _norm(expected_field_probe_method)
    )
    expected_output_artifact = (
        None
        if expected_field_probe_output_artifact_id is None
        else str(expected_field_probe_output_artifact_id).strip()
    )
    expected_output_digest = (
        None
        if expected_field_probe_output_digest is None
        else str(expected_field_probe_output_digest).strip()
    )
    expected_grid = None if expected_sample_grid_id is None else str(expected_sample_grid_id).strip()
    expected_grid_digest = (
        None if expected_sample_grid_digest is None else str(expected_sample_grid_digest).strip()
    )
    expected_count = None if expected_sample_count is None else int(expected_sample_count)
    expected_angle = None if expected_angle_unit is None else _norm(expected_angle_unit)
    expected_basis = None if expected_angle_basis is None else _norm(expected_angle_basis)
    expected_frame = None if expected_component_frame is None else _norm(expected_component_frame)
    expected_sign = (
        None if expected_torque_sign_convention is None else _norm(expected_torque_sign_convention)
    )
    output_required = (
        _truthy(require_field_probe_output_artifact)
        or expected_output_artifact is not None
        or expected_output_digest is not None
    )

    normalized_rows = []
    source_tools = []
    result_set_ids = []
    export_artifact_ids = []
    probe_ids = []
    probe_methods = []
    output_artifact_ids = []
    output_digests = []
    output_paths = []
    sample_grid_ids = []
    sample_grid_digests = []
    sample_counts = []
    angle_units = []
    angle_bases = []
    component_frames = []
    torque_sign_conventions = []
    angles = []
    br_values = []
    bt_values = []
    missing_angles = []
    missing_br = []
    missing_bt = []
    missing_source_tool = []
    missing_result_set = []
    missing_export = []
    missing_probe = []
    missing_probe_method = []
    missing_output_artifact = []
    missing_output_digest = []
    missing_output_path = []
    missing_sample_grid = []
    missing_sample_grid_digest = []
    missing_sample_count = []
    missing_radius = []
    missing_axial_length = []
    bad_radius = []
    bad_axial_length = []
    bad_symmetry_factor = []

    for index, raw in enumerate(data):
        row = dict(raw)
        source_tool = str(_first(row, ("source_tool", "tool", "source")) or "").strip()
        source_tool_key = _norm(source_tool) if source_tool else ""
        result_set_id = str(_first(row, ("result_set_id", "resultset_id", "result_id")) or "").strip()
        export_artifact_id = str(
            _first(row, ("export_artifact_id", "field_export_artifact_id", "graph_export_artifact_id")) or ""
        ).strip()
        probe_id = str(_first(row, ("field_probe_id", "probe_id", "line_probe_id")) or "").strip()
        probe_method = str(
            _first(row, ("field_probe_method", "probe_method", "export_method")) or ""
        ).strip()
        output_artifact_id = str(
            _first(
                row,
                (
                    "field_probe_output_artifact_id",
                    "probe_output_artifact_id",
                    "field_table_artifact_id",
                    "export_output_artifact_id",
                ),
            )
            or ""
        ).strip()
        output_digest = str(
            _first(
                row,
                (
                    "field_probe_output_digest",
                    "probe_output_digest",
                    "field_table_digest",
                    "export_output_digest",
                ),
            )
            or ""
        ).strip()
        output_path = str(
            _first(
                row,
                (
                    "field_probe_output_path",
                    "probe_output_path",
                    "field_table_path",
                    "export_output_path",
                ),
            )
            or ""
        ).strip()
        sample_grid_id = str(
            _first(row, ("sample_grid_id", "angle_grid_id", "sampling_grid_id")) or ""
        ).strip()
        sample_grid_digest = str(
            _first(row, ("sample_grid_digest", "angle_grid_digest", "sampling_grid_digest")) or ""
        ).strip()
        sample_count = _float_or_none(
            _first(row, ("sample_count", "sample_row_count", "angle_sample_count"))
        )
        angle_unit = _norm(_first(row, ("angle_unit", "theta_unit", "sample_angle_unit")) or "")
        angle_basis = _norm(_first(row, ("angle_basis", "theta_basis", "sample_angle_basis")) or "")
        component_frame = _norm(
            _first(row, ("component_frame", "field_component_frame", "B_component_frame")) or ""
        )
        torque_sign = _norm(
            _first(row, ("torque_sign_convention", "sign_convention", "shear_sign_convention")) or ""
        )
        angle_value = _float_or_none(
            _first(
                row,
                (
                    "RotorAngle_deg",
                    "theta_mech_deg",
                    "angle_deg",
                    "theta_deg",
                    "RotorAngle_rad",
                    "theta_rad",
                    "angle_rad",
                ),
            )
        )
        br = _float_or_none(_first(row, ("Br_T", "br_T", "B_radial_T", "radial_flux_density_T")))
        bt = _float_or_none(
            _first(row, ("Bt_T", "bt_T", "B_tangential_T", "tangential_flux_density_T"))
        )
        radius_m = _float_or_none(_first(row, ("radius_m", "sample_radius_m", "airgap_radius_m")))
        axial_length_m = _float_or_none(_first(row, ("axial_length_m", "stack_length_m", "length_m")))
        symmetry_factor = _float_or_none(_first(row, ("symmetry_factor", "sector_symmetry_factor")))

        if source_tool:
            source_tools.append(source_tool_key)
        else:
            missing_source_tool.append(index)
        if result_set_id:
            result_set_ids.append(result_set_id)
        elif expected_result is not None:
            missing_result_set.append(index)
        if export_artifact_id:
            export_artifact_ids.append(export_artifact_id)
        elif expected_export is not None:
            missing_export.append(index)
        if probe_id:
            probe_ids.append(probe_id)
        elif expected_probe is not None:
            missing_probe.append(index)
        if probe_method:
            probe_methods.append(_norm(probe_method))
        elif expected_method is not None:
            missing_probe_method.append(index)
        if output_artifact_id:
            output_artifact_ids.append(output_artifact_id)
        elif output_required:
            missing_output_artifact.append(index)
        if output_digest:
            output_digests.append(output_digest)
        elif output_required:
            missing_output_digest.append(index)
        if output_path:
            output_paths.append(output_path)
        elif output_required:
            missing_output_path.append(index)
        if sample_grid_id:
            sample_grid_ids.append(sample_grid_id)
        elif expected_grid is not None:
            missing_sample_grid.append(index)
        if sample_grid_digest:
            sample_grid_digests.append(sample_grid_digest)
        elif expected_grid_digest is not None:
            missing_sample_grid_digest.append(index)
        if sample_count is not None:
            sample_counts.append(int(sample_count))
        elif expected_count is not None:
            missing_sample_count.append(index)
        if angle_unit:
            angle_units.append(angle_unit)
        if angle_basis:
            angle_bases.append(angle_basis)
        if component_frame:
            component_frames.append(component_frame)
        if torque_sign:
            torque_sign_conventions.append(torque_sign)
        if _finite(angle_value):
            angles.append(float(angle_value))
        else:
            missing_angles.append(index)
        if _finite(br):
            br_values.append(float(br))
        else:
            missing_br.append(index)
        if _finite(bt):
            bt_values.append(float(bt))
        else:
            missing_bt.append(index)
        if radius_m is None:
            missing_radius.append(index)
        elif not _positive(radius_m):
            bad_radius.append(index)
        if axial_length_m is None:
            missing_axial_length.append(index)
        elif not _positive(axial_length_m):
            bad_axial_length.append(index)
        if symmetry_factor is None or not _positive(symmetry_factor):
            bad_symmetry_factor.append(index)

        normalized_rows.append({
            "index": index,
            "angle": angle_value,
            "Br_T": br,
            "Bt_T": bt,
            "source_tool": source_tool or None,
            "result_set_id": result_set_id or None,
            "export_artifact_id": export_artifact_id or None,
            "field_probe_id": probe_id or None,
            "field_probe_method": _norm(probe_method) if probe_method else None,
            "field_probe_output_artifact_id": output_artifact_id or None,
            "field_probe_output_digest": output_digest or None,
            "field_probe_output_path": output_path or None,
            "sample_grid_id": sample_grid_id or None,
            "sample_grid_digest": sample_grid_digest or None,
            "sample_count": int(sample_count) if sample_count is not None else None,
        })

    unique_source_tools = sorted(set(source_tools))
    unique_result_sets = sorted(set(result_set_ids))
    unique_exports = sorted(set(export_artifact_ids))
    unique_probe_ids = sorted(set(probe_ids))
    unique_probe_methods = sorted(set(probe_methods))
    unique_output_artifacts = sorted(set(output_artifact_ids))
    unique_output_digests = sorted(set(output_digests))
    unique_output_paths = sorted(set(output_paths))
    unique_sample_grid_ids = sorted(set(sample_grid_ids))
    unique_sample_grid_digests = sorted(set(sample_grid_digests))
    unique_sample_counts = sorted(set(sample_counts))
    unique_angle_units = sorted(set(angle_units))
    unique_angle_bases = sorted(set(angle_bases))
    unique_component_frames = sorted(set(component_frames))
    unique_signs = sorted(set(torque_sign_conventions))
    duplicate_angles = sorted(value for value in set(angles) if angles.count(value) > 1)

    checks = {
        "sample_rows_present": len(data) >= int(min_rows),
        "source_tool_recorded": not missing_source_tool,
        "source_tool_is_jmag": (
            not missing_source_tool
            and set(unique_source_tools).issubset({"jmag", "jmag_designer", "jmagdesigner"})
        ),
        "angles_recorded_and_finite": not missing_angles,
        "sample_angles_unique": not duplicate_angles,
        "br_t_recorded_and_finite": not missing_br,
        "bt_t_recorded_and_finite": not missing_bt,
        "angle_unit_recorded": bool(unique_angle_units),
        "expected_angle_unit_matches": expected_angle is None or unique_angle_units == [expected_angle],
        "component_frame_recorded": bool(unique_component_frames),
        "expected_component_frame_matches": (
            expected_frame is None or unique_component_frames == [expected_frame]
        ),
        "radius_m_recorded": not missing_radius,
        "radius_m_positive": not missing_radius and not bad_radius,
        "axial_length_m_recorded": not missing_axial_length,
        "axial_length_m_positive": not missing_axial_length and not bad_axial_length,
        "symmetry_factor_positive": not bad_symmetry_factor,
        "result_set_id_recorded_when_expected": expected_result is None or not missing_result_set,
        "expected_result_set_id_matches": expected_result is None or unique_result_sets == [expected_result],
        "export_artifact_id_recorded_when_expected": expected_export is None or not missing_export,
        "expected_export_artifact_id_matches": expected_export is None or unique_exports == [expected_export],
        "field_probe_id_recorded_when_expected": expected_probe is None or not missing_probe,
        "expected_field_probe_id_matches": expected_probe is None or unique_probe_ids == [expected_probe],
        "field_probe_method_recorded_when_expected": expected_method is None or not missing_probe_method,
        "expected_field_probe_method_matches": (
            expected_method is None or unique_probe_methods == [expected_method]
        ),
        "field_probe_output_artifact_id_recorded": (
            not output_required or not missing_output_artifact
        ),
        "expected_field_probe_output_artifact_id_matches": (
            expected_output_artifact is None or unique_output_artifacts == [expected_output_artifact]
        ),
        "field_probe_output_digest_recorded": not output_required or not missing_output_digest,
        "expected_field_probe_output_digest_matches": (
            expected_output_digest is None or unique_output_digests == [expected_output_digest]
        ),
        "field_probe_output_path_recorded": not output_required or not missing_output_path,
        "sample_grid_id_recorded_when_expected": expected_grid is None or not missing_sample_grid,
        "expected_sample_grid_id_matches": expected_grid is None or unique_sample_grid_ids == [expected_grid],
        "sample_grid_digest_recorded_when_expected": (
            expected_grid_digest is None or not missing_sample_grid_digest
        ),
        "expected_sample_grid_digest_matches": (
            expected_grid_digest is None or unique_sample_grid_digests == [expected_grid_digest]
        ),
        "sample_count_recorded_when_expected": expected_count is None or not missing_sample_count,
        "expected_sample_count_matches": (
            expected_count is None
            or unique_sample_counts == [expected_count]
            and expected_count == len(data)
        ),
    }
    if expected_basis is not None:
        checks["expected_angle_basis_matches"] = unique_angle_bases == [expected_basis]
    else:
        checks["angle_basis_recorded"] = bool(unique_angle_bases)
    if expected_sign is not None:
        checks["expected_torque_sign_convention_matches"] = unique_signs == [expected_sign]
    else:
        checks["torque_sign_convention_recorded"] = bool(unique_signs)

    return {
        "policy": "jmag_airgap_flux_sample_metadata_gate",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "n_rows": len(data),
        "rows": normalized_rows,
        "source_tools": unique_source_tools,
        "result_set_ids": unique_result_sets,
        "export_artifact_ids": unique_exports,
        "field_probe_ids": unique_probe_ids,
        "field_probe_methods": unique_probe_methods,
        "field_probe_output_artifact_ids": unique_output_artifacts,
        "field_probe_output_digests": unique_output_digests,
        "field_probe_output_paths": unique_output_paths,
        "sample_grid_ids": unique_sample_grid_ids,
        "sample_grid_digests": unique_sample_grid_digests,
        "sample_counts": unique_sample_counts,
        "angle_units": unique_angle_units,
        "angle_bases": unique_angle_bases,
        "component_frames": unique_component_frames,
        "torque_sign_conventions": unique_signs,
        "duplicate_angles": duplicate_angles,
        "missing_angle_rows": missing_angles,
        "missing_br_t_rows": missing_br,
        "missing_bt_t_rows": missing_bt,
        "missing_output_artifact_rows": missing_output_artifact,
        "missing_output_digest_rows": missing_output_digest,
        "missing_output_path_rows": missing_output_path,
        "checks": checks,
        "version_note": (
            "Use this after JMAG column/symmetry metadata and before Arkkio, "
            "Maxwell-shear, or radia-ngsolve air-gap torque checks so Br/Bt "
            "sample rows remain bound to result/export/probe/output identity."
        ),
    }


def jmag_airgap_torque_integration_package_gate(
    sample_metadata_gate,
    torque_package,
    *,
    expected_input_field_table_artifact_id=None,
    expected_input_field_table_digest=None,
    expected_model_input_artifact_id=None,
    expected_model_input_digest=None,
    expected_model_input_path=None,
    expected_export_recipe_artifact_id=None,
    expected_export_recipe_digest=None,
    expected_export_recipe_path=None,
    expected_parameter_set_artifact_id=None,
    expected_parameter_set_digest=None,
    expected_parameter_set_path=None,
    expected_objective_observable_id=None,
    expected_objective_observable_family=None,
    expected_sample_grid_id=None,
    expected_sample_grid_digest=None,
    expected_integration_method="maxwell_shear_from_br_bt_samples",
    expected_torque_output_artifact_id=None,
    expected_torque_output_digest=None,
    expected_torque_output_schema_id=None,
    expected_torque_output_columns=None,
    expected_torque_output_units=None,
    expected_torque_convention_schema_id=None,
    expected_torque_component_basis_schema_id=None,
    expected_torque_postprocess_row_convention_schema_id=None,
    expected_torque_nm=None,
    torque_abs_tol=1.0e-12,
    expected_created_at_utc=None,
    expected_run_timestamp_utc=None,
    expected_solver_version=None,
    expected_radia_mcp_version=None,
    max_created_run_skew_s=None,
    min_timing_sections=4,
    require_torque_output_artifact=False,
    require_torque_output_schema=False,
    require_torque_convention_schema=False,
    require_torque_component_basis_schema=False,
    require_torque_postprocess_row_convention_schema=False,
    require_model_input_artifact=False,
    require_export_recipe_artifact=False,
    require_parameter_set_artifact=False,
    require_execution_metadata=False,
    require_timing_breakdown=False,
):
    """Bind a air-gap torque result to its Br/Bt input package.

    This gate sits after :func:`jmag_airgap_flux_sample_metadata_gate`.  It
    keeps the computed torque result from becoming a free-floating scalar by
    requiring the field table artifact, sample-grid artifact, and integration
    method to be repeated in the result package.
    """

    sample_gate = dict(sample_metadata_gate or {})
    package = dict(torque_package or {})

    def _norm(value):
        return str(value).strip().lower().replace("-", "_").replace(" ", "_")

    def _first(mapping, names):
        for name in names:
            if name in mapping and mapping[name] not in (None, ""):
                return mapping[name]
        return None

    def _string_or_none(value):
        if value in (None, ""):
            return None
        return str(value).strip()

    def _string_list(value):
        if value in (None, ""):
            return []
        if isinstance(value, str):
            text = value.strip()
            if not text:
                return []
            for sep in (";", "\n"):
                text = text.replace(sep, ",")
            return [item.strip() for item in text.split(",") if item.strip()]
        if isinstance(value, dict):
            return [str(item).strip() for item in value.keys() if str(item).strip()]
        return [str(item).strip() for item in value if str(item).strip()]

    def _unit_mapping(value):
        if value in (None, ""):
            return {}
        if isinstance(value, dict):
            return {
                str(key).strip(): str(unit).strip()
                for key, unit in value.items()
                if str(key).strip()
            }
        pairs = {}
        for item in _string_list(value):
            if ":" in item:
                key, unit = item.split(":", 1)
            elif "=" in item:
                key, unit = item.split("=", 1)
            else:
                continue
            key = key.strip()
            if key:
                pairs[key] = unit.strip()
        return pairs

    def _float_or_none(value):
        if value in (None, ""):
            return None
        return float(value)

    def _timing_duration(value):
        if isinstance(value, dict):
            for key in ("duration_s", "elapsed_s", "runtime_s", "seconds", "s"):
                if value.get(key) not in (None, ""):
                    return float(value[key])
            return None
        if value in (None, ""):
            return None
        return float(value)

    def _timing_rows(value):
        if value in (None, ""):
            return []
        rows = []
        if isinstance(value, dict):
            iterable = value.items()
        else:
            iterable = enumerate(value)
        for key, entry in iterable:
            if isinstance(entry, dict):
                name = str(
                    entry.get(
                        "name",
                        entry.get("stage", entry.get("phase", entry.get("label", key))),
                    )
                ).strip()
                duration = _timing_duration(entry)
            elif isinstance(entry, (list, tuple)) and len(entry) >= 2:
                name = str(entry[0]).strip()
                duration = _timing_duration(entry[1])
            else:
                name = str(key).strip()
                duration = _timing_duration(entry)
            if name and duration is not None:
                rows.append({"name": name, "duration_s": duration})
        return rows

    def _truthy(value):
        if isinstance(value, bool):
            return value
        return _norm(value) in {"1", "true", "yes", "y", "required", "require"}

    sample_status = sample_gate.get("status")
    sample_field_artifacts = list(sample_gate.get("field_probe_output_artifact_ids") or [])
    sample_field_digests = list(sample_gate.get("field_probe_output_digests") or [])
    sample_field_paths = list(sample_gate.get("field_probe_output_paths") or [])
    sample_grid_ids = list(sample_gate.get("sample_grid_ids") or [])
    sample_grid_digests = list(sample_gate.get("sample_grid_digests") or [])
    sample_counts = list(sample_gate.get("sample_counts") or [])
    sample_frames = list(sample_gate.get("component_frames") or [])
    sample_signs = list(sample_gate.get("torque_sign_conventions") or [])

    expected_input_artifact = (
        None
        if expected_input_field_table_artifact_id is None
        else str(expected_input_field_table_artifact_id).strip()
    )
    expected_input_digest = (
        None
        if expected_input_field_table_digest is None
        else str(expected_input_field_table_digest).strip()
    )
    expected_model_input_artifact = (
        None
        if expected_model_input_artifact_id is None
        else str(expected_model_input_artifact_id).strip()
    )
    expected_model_input_digest_text = (
        None
        if expected_model_input_digest is None
        else str(expected_model_input_digest).strip()
    )
    expected_model_input_path_text = (
        None
        if expected_model_input_path is None
        else str(expected_model_input_path).strip()
    )
    expected_export_recipe_artifact = (
        None
        if expected_export_recipe_artifact_id is None
        else str(expected_export_recipe_artifact_id).strip()
    )
    expected_export_recipe_digest_text = (
        None
        if expected_export_recipe_digest is None
        else str(expected_export_recipe_digest).strip()
    )
    expected_export_recipe_path_text = (
        None
        if expected_export_recipe_path is None
        else str(expected_export_recipe_path).strip()
    )
    expected_parameter_set_artifact = (
        None
        if expected_parameter_set_artifact_id is None
        else str(expected_parameter_set_artifact_id).strip()
    )
    expected_parameter_set_digest_text = (
        None
        if expected_parameter_set_digest is None
        else str(expected_parameter_set_digest).strip()
    )
    expected_parameter_set_path_text = (
        None
        if expected_parameter_set_path is None
        else str(expected_parameter_set_path).strip()
    )
    expected_objective_observable = (
        None
        if expected_objective_observable_id is None
        else str(expected_objective_observable_id).strip()
    )
    expected_objective_family = (
        None
        if expected_objective_observable_family is None
        else _norm(expected_objective_observable_family)
    )
    expected_grid = None if expected_sample_grid_id is None else str(expected_sample_grid_id).strip()
    expected_grid_digest = (
        None if expected_sample_grid_digest is None else str(expected_sample_grid_digest).strip()
    )
    expected_method = (
        None if expected_integration_method is None else _norm(expected_integration_method)
    )
    expected_output_artifact = (
        None
        if expected_torque_output_artifact_id is None
        else str(expected_torque_output_artifact_id).strip()
    )
    expected_output_digest = (
        None
        if expected_torque_output_digest is None
        else str(expected_torque_output_digest).strip()
    )
    expected_output_schema = (
        None
        if expected_torque_output_schema_id is None
        else str(expected_torque_output_schema_id).strip()
    )
    expected_output_columns = _string_list(expected_torque_output_columns)
    expected_output_units = _unit_mapping(expected_torque_output_units)
    expected_torque_convention_schema = (
        None
        if expected_torque_convention_schema_id is None
        else str(expected_torque_convention_schema_id).strip()
    )
    expected_torque_component_basis_schema = (
        None
        if expected_torque_component_basis_schema_id is None
        else str(expected_torque_component_basis_schema_id).strip()
    )
    expected_torque_postprocess_row_convention_schema = (
        None
        if expected_torque_postprocess_row_convention_schema_id is None
        else str(expected_torque_postprocess_row_convention_schema_id).strip()
    )
    expected_created_at = (
        None
        if expected_created_at_utc is None
        else str(expected_created_at_utc).strip()
    )
    expected_run_timestamp = (
        None
        if expected_run_timestamp_utc is None
        else str(expected_run_timestamp_utc).strip()
    )
    expected_solver_version_text = (
        None if expected_solver_version is None else str(expected_solver_version).strip()
    )
    expected_radia_mcp_version_text = (
        None if expected_radia_mcp_version is None else str(expected_radia_mcp_version).strip()
    )
    output_required = (
        _truthy(require_torque_output_artifact)
        or expected_output_artifact is not None
        or expected_output_digest is not None
    )
    output_schema_required = (
        _truthy(require_torque_output_schema)
        or expected_output_schema is not None
        or bool(expected_output_columns)
        or bool(expected_output_units)
    )
    torque_convention_schema_required = (
        _truthy(require_torque_convention_schema)
        or expected_torque_convention_schema is not None
    )
    torque_component_basis_schema_required = (
        _truthy(require_torque_component_basis_schema)
        or expected_torque_component_basis_schema is not None
    )
    torque_postprocess_row_convention_schema_required = (
        _truthy(require_torque_postprocess_row_convention_schema)
        or expected_torque_postprocess_row_convention_schema is not None
    )
    model_input_required = (
        _truthy(require_model_input_artifact)
        or expected_model_input_artifact is not None
        or expected_model_input_digest_text is not None
        or expected_model_input_path_text is not None
    )
    model_input_digest_required = (
        _truthy(require_model_input_artifact)
        or expected_model_input_digest_text is not None
    )
    model_input_path_required = (
        _truthy(require_model_input_artifact)
        or expected_model_input_path_text is not None
    )
    export_recipe_required = (
        _truthy(require_export_recipe_artifact)
        or expected_export_recipe_artifact is not None
        or expected_export_recipe_digest_text is not None
        or expected_export_recipe_path_text is not None
    )
    export_recipe_digest_required = (
        _truthy(require_export_recipe_artifact)
        or expected_export_recipe_digest_text is not None
    )
    export_recipe_path_required = (
        _truthy(require_export_recipe_artifact)
        or expected_export_recipe_path_text is not None
    )
    parameter_set_required = (
        _truthy(require_parameter_set_artifact)
        or expected_parameter_set_artifact is not None
        or expected_parameter_set_digest_text is not None
        or expected_parameter_set_path_text is not None
    )
    parameter_set_digest_required = (
        _truthy(require_parameter_set_artifact)
        or expected_parameter_set_digest_text is not None
    )
    parameter_set_path_required = (
        _truthy(require_parameter_set_artifact)
        or expected_parameter_set_path_text is not None
    )
    execution_metadata_required = (
        _truthy(require_execution_metadata)
        or expected_created_at is not None
        or expected_run_timestamp is not None
        or expected_solver_version_text is not None
        or expected_radia_mcp_version_text is not None
        or max_created_run_skew_s is not None
    )
    created_at_required = expected_created_at is not None or max_created_run_skew_s is not None
    timing_breakdown_required = _truthy(require_timing_breakdown)

    input_artifact = _string_or_none(
        _first(
            package,
            (
                "input_field_table_artifact_id",
                "field_probe_output_artifact_id",
                "field_table_artifact_id",
                "source_field_table_artifact_id",
            ),
        )
    )
    input_digest = _string_or_none(
        _first(
            package,
            (
                "input_field_table_digest",
                "field_probe_output_digest",
                "field_table_digest",
                "source_field_table_digest",
            ),
        )
    )
    input_path = _string_or_none(
        _first(
            package,
            (
                "input_field_table_path",
                "field_probe_output_path",
                "field_table_path",
                "source_field_table_path",
            ),
        )
    )
    model_input_artifact = _string_or_none(
        _first(
            package,
            (
                "model_input_artifact_id",
                "project_artifact_id",
                "jproj_artifact_id",
                "input_project_artifact_id",
                "source_project_artifact_id",
            ),
        )
    )
    model_input_digest = _string_or_none(
        _first(
            package,
            (
                "model_input_digest",
                "project_digest",
                "jproj_digest",
                "input_project_digest",
                "source_project_digest",
            ),
        )
    )
    model_input_path = _string_or_none(
        _first(
            package,
            (
                "model_input_path",
                "project_path",
                "jproj_path",
                "input_project_path",
                "source_project_path",
            ),
        )
    )
    export_recipe_artifact = _string_or_none(
        _first(
            package,
            (
                "export_recipe_artifact_id",
                "export_script_artifact_id",
                "export_macro_artifact_id",
                "table_export_recipe_artifact_id",
                "postprocess_recipe_artifact_id",
            ),
        )
    )
    export_recipe_digest = _string_or_none(
        _first(
            package,
            (
                "export_recipe_digest",
                "export_recipe_sha256",
                "export_script_digest",
                "export_macro_digest",
                "table_export_recipe_digest",
                "postprocess_recipe_digest",
            ),
        )
    )
    export_recipe_path = _string_or_none(
        _first(
            package,
            (
                "export_recipe_path",
                "export_recipe_file",
                "export_script_path",
                "export_macro_path",
                "table_export_recipe_path",
                "postprocess_recipe_path",
            ),
        )
    )
    parameter_set_artifact = _string_or_none(
        _first(
            package,
            (
                "parameter_set_artifact_id",
                "design_parameter_set_artifact_id",
                "torque_parameter_set_artifact_id",
                "optimization_parameter_set_artifact_id",
                "operating_point_parameter_set_artifact_id",
            ),
        )
    )
    parameter_set_digest = _string_or_none(
        _first(
            package,
            (
                "parameter_set_digest",
                "parameter_set_sha256",
                "design_parameter_set_digest",
                "torque_parameter_set_digest",
                "optimization_parameter_set_digest",
            ),
        )
    )
    parameter_set_path = _string_or_none(
        _first(
            package,
            (
                "parameter_set_path",
                "parameter_set_file",
                "design_parameter_set_path",
                "torque_parameter_set_path",
                "optimization_parameter_set_path",
            ),
        )
    )
    objective_observable = _string_or_none(
        _first(
            package,
            (
                "objective_observable_id",
                "torque_objective_observable_id",
                "optimization_objective_observable_id",
                "objective_id",
                "objective_function_id",
            ),
        )
    )
    objective_family_raw = _string_or_none(
        _first(
            package,
            (
                "objective_observable_family",
                "torque_objective_observable_family",
                "optimization_objective_observable_family",
                "objective_family",
                "objective_kind",
            ),
        )
    )
    objective_family = _norm(objective_family_raw) if objective_family_raw else None
    input_grid = _string_or_none(
        _first(package, ("sample_grid_id", "input_sample_grid_id", "angle_grid_id"))
    )
    input_grid_digest = _string_or_none(
        _first(package, ("sample_grid_digest", "input_sample_grid_digest", "angle_grid_digest"))
    )
    input_sample_count = _float_or_none(
        _first(package, ("sample_count", "input_sample_count", "angle_sample_count"))
    )
    integration_method = _norm(
        _first(package, ("integration_method", "torque_integration_method", "method")) or ""
    )
    integration_policy = _norm(
        _first(package, ("integration_policy", "torque_integration_policy", "policy")) or ""
    )
    component_frame = _norm(
        _first(package, ("component_frame", "field_component_frame", "B_component_frame")) or ""
    )
    torque_sign = _norm(
        _first(package, ("torque_sign_convention", "sign_convention", "shear_sign_convention")) or ""
    )
    output_artifact = _string_or_none(
        _first(
            package,
            (
                "torque_output_artifact_id",
                "output_artifact_id",
                "integration_output_artifact_id",
            ),
        )
    )
    output_digest = _string_or_none(
        _first(package, ("torque_output_digest", "output_digest", "integration_output_digest"))
    )
    output_path = _string_or_none(
        _first(package, ("torque_output_path", "output_path", "integration_output_path"))
    )
    output_schema = _string_or_none(
        _first(
            package,
            (
                "torque_output_schema_id",
                "torque_table_schema_id",
                "integration_output_schema_id",
                "output_schema_id",
            ),
        )
    )
    output_columns = _string_list(
        _first(
            package,
            (
                "torque_output_columns",
                "torque_table_columns",
                "integration_output_columns",
                "output_columns",
                "columns",
            ),
        )
    )
    output_units = _unit_mapping(
        _first(
            package,
            (
                "torque_output_units",
                "torque_table_units",
                "integration_output_units",
                "output_units",
                "column_units",
                "units",
            ),
        )
    )
    torque_convention_schema = _string_or_none(
        _first(
            package,
            (
                "torque_convention_schema_id",
                "torque_physics_convention_schema_id",
                "airgap_torque_convention_schema_id",
                "physics_convention_schema_id",
            ),
        )
    )
    torque_component_basis_schema = _string_or_none(
        _first(
            package,
            (
                "torque_component_basis_schema_id",
                "airgap_torque_component_basis_schema_id",
                "br_bt_component_basis_schema_id",
                "field_component_basis_schema_id",
                "component_basis_schema_id",
            ),
        )
    )
    torque_postprocess_row_convention_schema = _string_or_none(
        _first(
            package,
            (
                "torque_postprocess_row_convention_schema_id",
                "postprocess_row_convention_schema_id",
                "torque_row_convention_schema_id",
                "airgap_torque_row_convention_schema_id",
                "postprocess_convention_schema_id",
            ),
        )
    )
    torque_nm = _float_or_none(_first(package, ("torque_Nm", "torque_nm", "T_Nm")))
    created_at = _string_or_none(
        _first(package, ("created_at_utc", "artifact_created_at_utc", "created_at"))
    )
    run_timestamp = _string_or_none(
        _first(
            package,
            ("run_timestamp_utc", "executed_at_utc", "run_date_utc", "date_utc", "run_date"),
        )
    )
    solver_version = _string_or_none(
        _first(package, ("solver_version", "jmag_version", "source_tool_version"))
    )
    radia_mcp_version = _string_or_none(
        _first(package, ("radia_mcp_version", "radia_ngsolve_version", "mcp_server_version"))
    )
    run_duration_s = _float_or_none(
        _first(package, ("run_duration_s", "elapsed_s", "runtime_s", "wall_time_s"))
    )
    timing_breakdown_rows = _timing_rows(
        _first(package, ("timing_breakdown_s", "timing_breakdown", "timing_breakdown_rows", "timings"))
    )
    timing_durations = [row["duration_s"] for row in timing_breakdown_rows]
    timing_total_s = sum(timing_durations) if timing_durations else None
    timing_sections_min = int(min_timing_sections or 0)
    top_durations = timing_durations[: max(timing_sections_min, 1)]
    timing_top_sections_descending = all(
        top_durations[index] >= top_durations[index + 1]
        for index in range(len(top_durations) - 1)
    )
    created_at_dt = _parse_utc_like_datetime(created_at)
    run_timestamp_dt = _parse_utc_like_datetime(run_timestamp)
    max_created_run_skew = (
        None if max_created_run_skew_s is None else float(max_created_run_skew_s)
    )
    created_run_skew_s = None
    if created_at_dt is not None and run_timestamp_dt is not None:
        created_run_skew_s = abs((created_at_dt - run_timestamp_dt).total_seconds())
    tol = float(torque_abs_tol)
    if tol < 0.0:
        raise ValueError("torque_abs_tol must be non-negative")

    sample_field_artifact = sample_field_artifacts[0] if len(sample_field_artifacts) == 1 else None
    sample_field_digest = sample_field_digests[0] if len(sample_field_digests) == 1 else None
    sample_grid_id = sample_grid_ids[0] if len(sample_grid_ids) == 1 else None
    sample_grid_digest = sample_grid_digests[0] if len(sample_grid_digests) == 1 else None
    sample_count = sample_counts[0] if len(sample_counts) == 1 else None
    sample_frame = sample_frames[0] if len(sample_frames) == 1 else None
    sample_sign = sample_signs[0] if len(sample_signs) == 1 else None

    torque_error = None
    if expected_torque_nm is not None and torque_nm is not None:
        torque_error = abs(float(torque_nm) - float(expected_torque_nm))

    checks = {
        "sample_metadata_gate_ok": sample_status == "ok",
        "sample_field_table_artifact_unique": len(sample_field_artifacts) == 1,
        "input_field_table_artifact_recorded": input_artifact is not None,
        "input_field_table_artifact_matches_sample_gate": (
            input_artifact is not None
            and sample_field_artifact is not None
            and input_artifact == sample_field_artifact
        ),
        "expected_input_field_table_artifact_matches": (
            expected_input_artifact is None or input_artifact == expected_input_artifact
        ),
        "input_field_table_digest_recorded": input_digest is not None,
        "input_field_table_digest_matches_sample_gate": (
            input_digest is not None
            and sample_field_digest is not None
            and input_digest == sample_field_digest
        ),
        "expected_input_field_table_digest_matches": (
            expected_input_digest is None or input_digest == expected_input_digest
        ),
        "input_field_table_path_recorded": input_path is not None or not sample_field_paths,
        "model_input_artifact_id_recorded": not model_input_required
        or model_input_artifact is not None,
        "model_input_digest_recorded": not model_input_digest_required
        or model_input_digest is not None,
        "model_input_path_recorded": not model_input_path_required
        or model_input_path is not None,
        "expected_model_input_artifact_id_matches": (
            expected_model_input_artifact is None
            or model_input_artifact == expected_model_input_artifact
        ),
        "expected_model_input_digest_matches": (
            expected_model_input_digest_text is None
            or model_input_digest == expected_model_input_digest_text
        ),
        "expected_model_input_path_matches": (
            expected_model_input_path_text is None
            or model_input_path == expected_model_input_path_text
        ),
        "export_recipe_artifact_id_recorded": not export_recipe_required
        or export_recipe_artifact is not None,
        "export_recipe_digest_recorded": not export_recipe_digest_required
        or export_recipe_digest is not None,
        "export_recipe_path_recorded": not export_recipe_path_required
        or export_recipe_path is not None,
        "expected_export_recipe_artifact_id_matches": (
            expected_export_recipe_artifact is None
            or export_recipe_artifact == expected_export_recipe_artifact
        ),
        "expected_export_recipe_digest_matches": (
            expected_export_recipe_digest_text is None
            or export_recipe_digest == expected_export_recipe_digest_text
        ),
        "expected_export_recipe_path_matches": (
            expected_export_recipe_path_text is None
            or export_recipe_path == expected_export_recipe_path_text
        ),
        "parameter_set_artifact_id_recorded": not parameter_set_required
        or parameter_set_artifact is not None,
        "parameter_set_digest_recorded": not parameter_set_digest_required
        or parameter_set_digest is not None,
        "parameter_set_path_recorded": not parameter_set_path_required
        or parameter_set_path is not None,
        "expected_parameter_set_artifact_id_matches": (
            expected_parameter_set_artifact is None
            or parameter_set_artifact == expected_parameter_set_artifact
        ),
        "expected_parameter_set_digest_matches": (
            expected_parameter_set_digest_text is None
            or parameter_set_digest == expected_parameter_set_digest_text
        ),
        "expected_parameter_set_path_matches": (
            expected_parameter_set_path_text is None
            or parameter_set_path == expected_parameter_set_path_text
        ),
        "objective_observable_id_recorded": expected_objective_observable is None
        or objective_observable is not None,
        "expected_objective_observable_id_matches": (
            expected_objective_observable is None
            or objective_observable == expected_objective_observable
        ),
        "objective_observable_family_recorded": expected_objective_family is None
        or objective_family is not None,
        "expected_objective_observable_family_matches": (
            expected_objective_family is None
            or objective_family == expected_objective_family
        ),
        "sample_grid_artifact_unique": len(sample_grid_ids) == 1,
        "sample_grid_id_recorded": input_grid is not None,
        "sample_grid_id_matches_sample_gate": (
            input_grid is not None and sample_grid_id is not None and input_grid == sample_grid_id
        ),
        "expected_sample_grid_id_matches": expected_grid is None or input_grid == expected_grid,
        "sample_grid_digest_recorded": input_grid_digest is not None,
        "sample_grid_digest_matches_sample_gate": (
            input_grid_digest is not None
            and sample_grid_digest is not None
            and input_grid_digest == sample_grid_digest
        ),
        "expected_sample_grid_digest_matches": (
            expected_grid_digest is None or input_grid_digest == expected_grid_digest
        ),
        "sample_count_matches_sample_gate": (
            input_sample_count is not None
            and sample_count is not None
            and int(input_sample_count) == int(sample_count)
        ),
        "integration_method_recorded": bool(integration_method),
        "expected_integration_method_matches": (
            expected_method is None or integration_method == expected_method
        ),
        "integration_policy_recorded": bool(integration_policy),
        "component_frame_matches_sample_gate": (
            not sample_frame or component_frame == sample_frame
        ),
        "torque_sign_convention_matches_sample_gate": (
            not sample_sign or torque_sign == sample_sign
        ),
        "torque_output_artifact_recorded": not output_required or output_artifact is not None,
        "expected_torque_output_artifact_matches": (
            expected_output_artifact is None or output_artifact == expected_output_artifact
        ),
        "torque_output_digest_recorded": not output_required or output_digest is not None,
        "expected_torque_output_digest_matches": (
            expected_output_digest is None or output_digest == expected_output_digest
        ),
        "torque_output_path_recorded": not output_required or output_path is not None,
        "torque_output_schema_id_recorded": not output_schema_required
        or bool(output_schema),
        "expected_torque_output_schema_id_matches": (
            expected_output_schema is None or output_schema == expected_output_schema
        ),
        "torque_output_columns_recorded": not output_schema_required
        or bool(output_columns),
        "expected_torque_output_columns_match": not expected_output_columns
        or output_columns == expected_output_columns,
        "torque_output_units_recorded": not output_schema_required
        or bool(output_units),
        "expected_torque_output_units_match": not expected_output_units
        or output_units == expected_output_units,
        "torque_convention_schema_id_recorded": not torque_convention_schema_required
        or bool(torque_convention_schema),
        "expected_torque_convention_schema_id_matches": (
            expected_torque_convention_schema is None
            or torque_convention_schema == expected_torque_convention_schema
        ),
        "torque_component_basis_schema_id_recorded": (
            not torque_component_basis_schema_required
            or bool(torque_component_basis_schema)
        ),
        "expected_torque_component_basis_schema_id_matches": (
            expected_torque_component_basis_schema is None
            or torque_component_basis_schema == expected_torque_component_basis_schema
        ),
        "torque_postprocess_row_convention_schema_id_recorded": (
            not torque_postprocess_row_convention_schema_required
            or bool(torque_postprocess_row_convention_schema)
        ),
        "expected_torque_postprocess_row_convention_schema_id_matches": (
            expected_torque_postprocess_row_convention_schema is None
            or torque_postprocess_row_convention_schema
            == expected_torque_postprocess_row_convention_schema
        ),
        "torque_nm_recorded_and_finite": torque_nm is not None and math.isfinite(torque_nm),
        "expected_torque_nm_matches": expected_torque_nm is None or (
            torque_error is not None and torque_error <= tol
        ),
        "created_at_utc_recorded": not created_at_required or bool(created_at),
        "created_at_utc_parseable": not created_at or created_at_dt is not None,
        "expected_created_at_utc_matches": expected_created_at is None
        or created_at == expected_created_at,
        "run_timestamp_utc_recorded": not execution_metadata_required
        or bool(run_timestamp),
        "run_timestamp_utc_parseable": not run_timestamp
        or run_timestamp_dt is not None,
        "expected_run_timestamp_utc_matches": expected_run_timestamp is None
        or run_timestamp == expected_run_timestamp,
        "created_run_timestamp_skew_within_limit": max_created_run_skew is None
        or (
            created_run_skew_s is not None
            and created_run_skew_s <= max_created_run_skew
        ),
        "solver_version_recorded": not execution_metadata_required
        or bool(solver_version),
        "expected_solver_version_matches": expected_solver_version_text is None
        or solver_version == expected_solver_version_text,
        "radia_mcp_version_recorded": not execution_metadata_required
        or bool(radia_mcp_version),
        "expected_radia_mcp_version_matches": expected_radia_mcp_version_text is None
        or radia_mcp_version == expected_radia_mcp_version_text,
        "run_duration_s_recorded": not (
            execution_metadata_required or timing_breakdown_required
        )
        or run_duration_s is not None,
        "run_duration_s_positive": run_duration_s is None or run_duration_s > 0.0,
        "timing_breakdown_recorded": not timing_breakdown_required
        or bool(timing_breakdown_rows),
        "timing_breakdown_has_required_sections": not timing_breakdown_required
        or len(timing_breakdown_rows) >= timing_sections_min,
        "timing_breakdown_sections_named": not timing_breakdown_required
        or all(bool(row["name"]) for row in timing_breakdown_rows),
        "timing_breakdown_values_nonnegative": not timing_breakdown_required
        or all(value >= 0.0 for value in timing_durations),
        "timing_breakdown_top_sections_descending": not timing_breakdown_required
        or timing_top_sections_descending,
        "timing_breakdown_total_within_run_duration": (
            run_duration_s is None
            or timing_total_s is None
            or timing_total_s <= run_duration_s * (1.0 + tol) + 1.0e-12
        ),
    }

    return {
        "policy": "jmag_airgap_torque_integration_package_gate",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "input_field_table_artifact_id": input_artifact,
        "input_field_table_digest": input_digest,
        "input_field_table_path": input_path,
        "model_input_artifact_id": model_input_artifact,
        "model_input_digest": model_input_digest,
        "model_input_path": model_input_path,
        "expected_model_input_artifact_id": expected_model_input_artifact,
        "expected_model_input_digest": expected_model_input_digest_text,
        "expected_model_input_path": expected_model_input_path_text,
        "export_recipe_artifact_id": export_recipe_artifact,
        "export_recipe_digest": export_recipe_digest,
        "export_recipe_path": export_recipe_path,
        "expected_export_recipe_artifact_id": expected_export_recipe_artifact,
        "expected_export_recipe_digest": expected_export_recipe_digest_text,
        "expected_export_recipe_path": expected_export_recipe_path_text,
        "parameter_set_artifact_id": parameter_set_artifact,
        "parameter_set_digest": parameter_set_digest,
        "parameter_set_path": parameter_set_path,
        "expected_parameter_set_artifact_id": expected_parameter_set_artifact,
        "expected_parameter_set_digest": expected_parameter_set_digest_text,
        "expected_parameter_set_path": expected_parameter_set_path_text,
        "objective_observable_id": objective_observable,
        "objective_observable_family": objective_family,
        "expected_objective_observable_id": expected_objective_observable,
        "expected_objective_observable_family": expected_objective_family,
        "sample_grid_id": input_grid,
        "sample_grid_digest": input_grid_digest,
        "sample_count": int(input_sample_count) if input_sample_count is not None else None,
        "integration_method": integration_method or None,
        "integration_policy": integration_policy or None,
        "component_frame": component_frame or None,
        "torque_sign_convention": torque_sign or None,
        "torque_output_artifact_id": output_artifact,
        "torque_output_digest": output_digest,
        "torque_output_path": output_path,
        "torque_output_schema_id": output_schema,
        "torque_output_columns": output_columns,
        "torque_output_units": output_units,
        "expected_torque_output_schema_id": expected_output_schema,
        "expected_torque_output_columns": expected_output_columns,
        "expected_torque_output_units": expected_output_units,
        "torque_convention_schema_id": torque_convention_schema,
        "expected_torque_convention_schema_id": expected_torque_convention_schema,
        "torque_component_basis_schema_id": torque_component_basis_schema,
        "expected_torque_component_basis_schema_id": expected_torque_component_basis_schema,
        "torque_postprocess_row_convention_schema_id": torque_postprocess_row_convention_schema,
        "expected_torque_postprocess_row_convention_schema_id": expected_torque_postprocess_row_convention_schema,
        "torque_Nm": torque_nm,
        "expected_torque_Nm": (
            None if expected_torque_nm is None else float(expected_torque_nm)
        ),
        "torque_abs_error": torque_error,
        "torque_abs_tol": tol,
        "created_at_utc": created_at,
        "expected_created_at_utc": expected_created_at,
        "run_timestamp_utc": run_timestamp,
        "expected_run_timestamp_utc": expected_run_timestamp,
        "created_run_timestamp_skew_s": created_run_skew_s,
        "max_created_run_skew_s": max_created_run_skew,
        "solver_version": solver_version,
        "expected_solver_version": expected_solver_version_text,
        "radia_mcp_version": radia_mcp_version,
        "expected_radia_mcp_version": expected_radia_mcp_version_text,
        "run_duration_s": run_duration_s,
        "timing_breakdown_s": timing_breakdown_rows,
        "timing_total_s": timing_total_s,
        "min_timing_sections": timing_sections_min,
        "model_input_artifact_required": model_input_required,
        "export_recipe_artifact_required": export_recipe_required,
        "parameter_set_artifact_required": parameter_set_required,
        "torque_output_schema_required": output_schema_required,
        "torque_convention_schema_required": torque_convention_schema_required,
        "torque_component_basis_schema_required": torque_component_basis_schema_required,
        "torque_postprocess_row_convention_schema_required": torque_postprocess_row_convention_schema_required,
        "execution_metadata_required": execution_metadata_required,
        "timing_breakdown_required": timing_breakdown_required,
        "checks": checks,
        "version_note": (
            "Use this after Br/Bt sample metadata and before promoting an "
            "air-gap torque scalar to a notebook, optimizer, or cross-solver "
            "comparison. The torque result must repeat field-table identity, "
            "sample-grid identity, integration method/policy, output artifact "
            "identity, output schema/column/unit identity, physics-convention "
            "schema identity, component-basis schema identity, model-input "
            "project identity, parameter-set/objective identity, sign/frame "
            "conventions, "
             "and execution metadata/timing when the row is promoted to "
             "reusable cross-validation JSON.  Keep created_at_utc close to "
             "run_timestamp_utc so copied torque-result artifacts are not "
             "mistaken for fresh notebook or cross-validation runs."
        ),
    }


def jmag_symmetry_sweep_coverage_gate(
    rows,
    pole_pairs,
    symmetry_factor,
    angle_column="RotorAngle_deg",
    angle_unit="deg",
    angle_basis="mechanical",
    endpoint_policy="included",
    step_rtol=1.0e-9,
):
    """Check that a JMAG sector sweep covers the expected symmetry span.

    JMAG motor exports often come from a symmetry-sector model.  Before torque
    harmonics or current snapshots consume the rows, verify that the angle
    column covers exactly ``360/symmetry_factor`` mechanical degrees, with the
    electrical span reported as ``pole_pairs`` times that value.
    """

    data = list(rows)
    if len(data) < 3:
        raise ValueError("at least three angle rows are required")
    pp = float(pole_pairs)
    sf = float(symmetry_factor)
    if pp <= 0.0:
        raise ValueError("pole_pairs must be > 0")
    if sf <= 0.0:
        raise ValueError("symmetry_factor must be > 0")
    unit = str(angle_unit).strip().lower().replace("degrees", "deg").replace("degree", "deg")
    basis = str(angle_basis).strip().lower().replace("-", "_").replace(" ", "_")
    endpoint = str(endpoint_policy).strip().lower().replace("-", "_").replace(" ", "_")
    column = str(angle_column)
    angles_raw = []
    missing_angle_rows = []
    for index, row in enumerate(data, start=1):
        if column not in row:
            missing_angle_rows.append(index)
            continue
        angles_raw.append(float(row[column]))
    angles_rad = []
    if unit == "deg":
        angles_rad = [math.radians(value) for value in angles_raw]
    elif unit == "rad":
        angles_rad = list(angles_raw)
    else:
        angles_rad = list(angles_raw)
    if basis == "mechanical":
        theta_mech = angles_rad
    elif basis == "electrical":
        theta_mech = [value / pp for value in angles_rad]
    else:
        theta_mech = angles_rad

    strictly_increasing = all(theta_mech[i + 1] > theta_mech[i] for i in range(len(theta_mech) - 1))
    steps = [theta_mech[i + 1] - theta_mech[i] for i in range(len(theta_mech) - 1)]
    mean_step = sum(steps) / len(steps) if steps else math.nan
    max_step_error = max((abs(step - mean_step) for step in steps), default=math.inf)
    expected_mech_span = 2.0 * math.pi / sf
    if endpoint == "included":
        covered_mech_span = theta_mech[-1] - theta_mech[0] if theta_mech else math.nan
    elif endpoint == "excluded":
        covered_mech_span = theta_mech[-1] - theta_mech[0] + mean_step if theta_mech else math.nan
    else:
        covered_mech_span = math.nan
    span_error = abs(covered_mech_span - expected_mech_span) if math.isfinite(covered_mech_span) else math.inf
    electrical_span = covered_mech_span * pp if math.isfinite(covered_mech_span) else math.nan
    expected_electrical_span = expected_mech_span * pp
    step_tol = float(step_rtol) * max(abs(mean_step), 1.0)
    span_tol = max(1.0e-12, float(step_rtol) * max(abs(expected_mech_span), 1.0))
    checks = {
        "angle_column_present": not missing_angle_rows and len(angles_raw) == len(data),
        "angle_unit_valid": unit in {"deg", "rad"},
        "angle_basis_valid": basis in {"mechanical", "electrical"},
        "endpoint_policy_valid": endpoint in {"included", "excluded"},
        "angles_strictly_increasing": strictly_increasing,
        "angle_step_uniform": max_step_error <= step_tol,
        "mechanical_sector_span_matches_symmetry": span_error <= span_tol,
        "electrical_span_matches_pole_pairs": (
            math.isfinite(electrical_span)
            and abs(electrical_span - expected_electrical_span) <= pp * span_tol
        ),
    }
    return {
        "policy": "jmag_symmetry_sweep_coverage_gate",
        "n_rows": len(data),
        "angle_column": column,
        "angle_unit": unit,
        "angle_basis": basis,
        "endpoint_policy": endpoint,
        "pole_pairs": pp,
        "symmetry_factor": sf,
        "expected_mechanical_span_deg": math.degrees(expected_mech_span),
        "covered_mechanical_span_deg": math.degrees(covered_mech_span) if math.isfinite(covered_mech_span) else None,
        "mechanical_span_abs_error_deg": math.degrees(span_error) if math.isfinite(span_error) else None,
        "expected_electrical_span_deg": math.degrees(expected_electrical_span),
        "covered_electrical_span_deg": math.degrees(electrical_span) if math.isfinite(electrical_span) else None,
        "mean_step_deg_mechanical": math.degrees(mean_step) if math.isfinite(mean_step) else None,
        "max_step_error_deg_mechanical": math.degrees(max_step_error) if math.isfinite(max_step_error) else None,
        "missing_angle_rows": missing_angle_rows,
        "checks": checks,
        "status": "ok" if all(checks.values()) else "needs_attention",
        "version_note": (
            "Use this after jmag_motor_table_column_metadata_gate and before "
            "torque/current/harmonic parsing so a symmetry-sector export covers "
            "the intended mechanical/electrical angle span."
        ),
    }


def jmag_angle_alignment_contract_gate(
    rows,
    pole_pairs,
    *,
    mechanical_angle_column="theta_mech_deg",
    electrical_angle_column="theta_e_deg",
    gamma_export_column="gamma_jmag_deg",
    gamma_reference_column="gamma_reference_deg",
    expected_gamma_offset_deg=0.0,
    rotor_electrical_offset_deg=0.0,
    expected_symmetry_factor=None,
    angle_tol_deg=1.0e-9,
):
    """Check JMAG/open-FE motor angle alignment before waveform comparison."""

    data = list(rows)
    if len(data) < 2:
        raise ValueError("at least two angle rows are required")
    pp = float(pole_pairs)
    if pp <= 0.0:
        raise ValueError("pole_pairs must be > 0")
    tol = float(angle_tol_deg)
    if tol < 0.0:
        raise ValueError("angle_tol_deg must be >= 0")

    mech_candidates = (
        mechanical_angle_column,
        "theta_mech_deg",
        "rotor_angle_mech_deg",
        "mechanical_angle_deg",
        "RotorAngle_deg",
    )
    elec_candidates = (
        electrical_angle_column,
        "theta_e_deg",
        "theta_electrical_deg",
        "electrical_angle_deg",
        "initial_electrical_angle_deg",
    )
    gamma_export_candidates = (
        gamma_export_column,
        "gamma_jmag_deg",
        "gamma_export_deg",
        "current_angle_jmag_deg",
        "current_angle_export_deg",
    )
    gamma_reference_candidates = (
        gamma_reference_column,
        "gamma_reference_deg",
        "gamma_open_fe_deg",
        "gamma_ref_deg",
        "current_angle_reference_deg",
    )

    def _first(row, names):
        for name in names:
            if name in row and row[name] is not None:
                return row[name]
        return None

    def _delta_deg(actual, expected):
        return ((float(actual) - float(expected) + 180.0) % 360.0) - 180.0

    normalized = []
    missing_rows = []
    theta_e_errors = []
    gamma_offset_errors = []
    symmetry_errors = []
    symmetry_values = []
    expected_sf = None if expected_symmetry_factor is None else float(expected_symmetry_factor)

    for index, row in enumerate(data, start=1):
        theta_m_raw = _first(row, mech_candidates)
        theta_e_raw = _first(row, elec_candidates)
        gamma_export_raw = _first(row, gamma_export_candidates)
        gamma_ref_raw = _first(row, gamma_reference_candidates)
        missing = []
        if theta_m_raw is None:
            missing.append("mechanical_angle")
        if theta_e_raw is None:
            missing.append("electrical_angle")
        if gamma_export_raw is None:
            missing.append("gamma_export")
        if gamma_ref_raw is None:
            missing.append("gamma_reference")
        if missing:
            missing_rows.append({"row": index, "missing": missing})
            continue

        theta_m = float(theta_m_raw)
        theta_e = float(theta_e_raw)
        gamma_export = float(gamma_export_raw)
        gamma_ref = float(gamma_ref_raw)
        expected_theta_e = pp * theta_m + float(rotor_electrical_offset_deg)
        theta_e_error = abs(_delta_deg(theta_e, expected_theta_e))
        gamma_offset = _delta_deg(gamma_export, gamma_ref)
        gamma_offset_error = abs(_delta_deg(gamma_offset, expected_gamma_offset_deg))
        sf_raw = row.get("symmetry_factor")
        sf_value = None if sf_raw is None else float(sf_raw)
        if sf_value is not None:
            symmetry_values.append(sf_value)
            if expected_sf is not None:
                symmetry_errors.append(abs(sf_value - expected_sf))

        theta_e_errors.append(theta_e_error)
        gamma_offset_errors.append(gamma_offset_error)
        normalized.append(
            {
                "row": index,
                "theta_mech_deg": theta_m,
                "theta_e_deg": theta_e,
                "expected_theta_e_deg": expected_theta_e,
                "theta_e_error_deg": theta_e_error,
                "gamma_export_deg": gamma_export,
                "gamma_reference_deg": gamma_ref,
                "gamma_offset_deg": gamma_offset,
                "gamma_offset_error_deg": gamma_offset_error,
                "symmetry_factor": sf_value,
            }
        )

    mech_angles = [row["theta_mech_deg"] for row in sorted(normalized, key=lambda item: item["theta_mech_deg"])]
    steps = [mech_angles[i + 1] - mech_angles[i] for i in range(len(mech_angles) - 1)]
    mean_step = sum(steps) / len(steps) if steps else math.nan
    max_step_error = max((abs(step - mean_step) for step in steps), default=0.0)
    max_theta_e_error = max(theta_e_errors, default=math.inf)
    max_gamma_error = max(gamma_offset_errors, default=math.inf)
    max_symmetry_error = max(symmetry_errors, default=0.0)
    symmetry_positive = all(value > 0.0 for value in symmetry_values)
    checks = {
        "required_angle_columns_present": not missing_rows and len(normalized) == len(data),
        "theta_e_matches_pole_pairs_and_offset": max_theta_e_error <= tol,
        "gamma_offset_matches_expected": max_gamma_error <= tol,
        "mechanical_angle_step_uniform": max_step_error <= tol,
        "symmetry_factor_positive_when_present": symmetry_positive,
        "symmetry_factor_matches_expected": expected_sf is None
        or (len(symmetry_values) == len(data) and max_symmetry_error <= tol),
    }
    return {
        "policy": "jmag_angle_alignment_contract_gate",
        "n_rows": len(data),
        "pole_pairs": pp,
        "rotor_electrical_offset_deg": float(rotor_electrical_offset_deg),
        "expected_gamma_offset_deg": float(expected_gamma_offset_deg),
        "expected_symmetry_factor": expected_sf,
        "mean_mechanical_step_deg": mean_step if math.isfinite(mean_step) else None,
        "max_mechanical_step_error_deg": max_step_error,
        "max_theta_e_error_deg": max_theta_e_error if math.isfinite(max_theta_e_error) else None,
        "max_gamma_offset_error_deg": max_gamma_error if math.isfinite(max_gamma_error) else None,
        "max_symmetry_factor_error": max_symmetry_error,
        "missing_rows": missing_rows,
        "rows": normalized,
        "checks": checks,
        "status": "ok" if all(checks.values()) else "needs_attention",
        "version_note": (
            "Use this after JMAG column/symmetry metadata and before torque or "
            "flux waveform comparison: mechanical angle, electrical angle, "
            "current-angle gamma offset, and symmetry factor must be explicit."
        ),
    }


def jmag_export_case_package_gate(
    artifacts,
    expected_case_id=None,
    expected_study_id=None,
    expected_result_set_id=None,
    required_kinds=("column_metadata", "symmetry_coverage", "value_table", "notebook_row"),
):
    """Check that JMAG-derived export artifacts belong to one case package.

    JMAG postprocessing often exports column metadata, sector/symmetry coverage,
    value tables, and selected notebook rows separately.  This gate keeps the
    case/study/result-set identity explicit so a downstream notebook cannot
    join a stale table from another sweep or design point.
    """

    rows_in = list(artifacts)
    if not rows_in:
        raise ValueError("artifacts must not be empty")

    def _norm(value):
        return str(value or "").strip().lower().replace("-", "_").replace(" ", "_")

    def _first(row, names):
        for name in names:
            if name in row and row[name] is not None:
                return row[name]
        return None

    required = tuple(_norm(kind) for kind in required_kinds)
    if not required:
        raise ValueError("required_kinds must not be empty")

    expected_policies = {
        "column_metadata": {"jmag_motor_table_column_metadata_gate"},
        "symmetry_coverage": {"jmag_symmetry_sweep_coverage_gate"},
        "value_table": {
            "pm_drive_terminal_table_health_gate",
            "pm_drive_loss_bucket_efficiency_gate",
            "pm_drive_efficiency_map_health_gate",
            "dq_torque_table_health",
            "torque_angle_table_export_health",
        },
        "notebook_row": {
            "pm_drive_operating_point_notebook_handoff_gate",
            "pm_drive_terminal_table_health_gate",
            "pm_drive_loss_bucket_efficiency_gate",
        },
    }

    details = []
    kind_counts = {}
    case_ids = []
    study_ids = []
    result_set_ids = []
    notebook_operating_points = []
    table_operating_points = []
    missing_case_id = []
    missing_study_id = []
    missing_result_set_id = []
    missing_notebook_operating_point = []
    bad_source_tool = []
    missing_paths = []
    unknown_kinds = []
    bad_upstream_status = []
    bad_upstream_policy = []

    for index, row in enumerate(rows_in, start=1):
        if not isinstance(row, dict):
            raise ValueError("each artifact must be a dictionary")
        kind = _norm(_first(row, ("kind", "artifact_kind", "type")))
        case_id = _first(row, ("case_id", "jmag_case_id", "design_case_id"))
        study_id = _first(row, ("study_id", "study_name", "study"))
        result_set_id = _first(row, ("result_set_id", "result_id", "dataset_id", "export_id"))
        source_tool = _first(row, ("source_tool", "tool", "source"))
        source_tool_norm = _norm(source_tool)
        path = _first(row, ("path", "file", "table_path", "artifact_path"))
        gate_policy = _first(row, ("gate_policy", "policy", "validator"))
        gate_policy_norm = _norm(gate_policy)
        status = _first(row, ("status", "gate_status", "validation_status"))
        status_norm = _norm(status)
        pass_flag = bool(row.get("pass", False))

        if kind:
            kind_counts[kind] = kind_counts.get(kind, 0) + 1
        else:
            unknown_kinds.append(index)
        if kind and kind not in set(required) | set(expected_policies):
            unknown_kinds.append(kind)
        if not case_id:
            missing_case_id.append(index)
        else:
            case_ids.append(str(case_id))
        if not study_id:
            missing_study_id.append(index)
        else:
            study_ids.append(str(study_id))
        if not result_set_id:
            missing_result_set_id.append(index)
        else:
            result_set_ids.append(str(result_set_id))
        if source_tool_norm not in {"jmag", "jmag_designer", "jmagdesigner"}:
            bad_source_tool.append({"index": index, "source_tool": source_tool})
        if not path:
            missing_paths.append(index)
        if not (pass_flag or status_norm in {"ok", "pass", "passed", "verified"}):
            bad_upstream_status.append({"index": index, "kind": kind, "status": status})
        if kind in expected_policies and gate_policy_norm not in expected_policies[kind]:
            bad_upstream_policy.append({
                "index": index,
                "kind": kind,
                "gate_policy": gate_policy,
                "expected": sorted(expected_policies[kind]),
            })

        if kind == "value_table":
            op_ids = row.get("operating_point_ids")
            if op_ids is not None:
                table_operating_points.extend(str(item) for item in op_ids)
            op_single = _first(row, ("operating_point_id", "op_id"))
            if op_single is not None:
                table_operating_points.append(str(op_single))
        if kind == "notebook_row":
            op_id = _first(row, ("operating_point_id", "op_id"))
            if not op_id:
                missing_notebook_operating_point.append(index)
            else:
                notebook_operating_points.append(str(op_id))

        details.append({
            "index": index,
            "kind": kind,
            "case_id": None if case_id is None else str(case_id),
            "study_id": None if study_id is None else str(study_id),
            "result_set_id": None if result_set_id is None else str(result_set_id),
            "source_tool": source_tool,
            "path": path,
            "gate_policy": gate_policy,
            "status": status,
            "pass": pass_flag,
        })

    required_set = set(required)
    present_set = set(kind_counts)
    unique_case_ids = sorted(set(case_ids))
    unique_study_ids = sorted(set(study_ids))
    unique_result_set_ids = sorted(set(result_set_ids))
    table_op_set = set(table_operating_points)
    checks = {
        "required_kinds_present": required_set.issubset(present_set),
        "no_unknown_kinds": not unknown_kinds,
        "case_ids_present": not missing_case_id,
        "case_ids_unique": len(unique_case_ids) == 1,
        "study_ids_present": not missing_study_id,
        "study_ids_unique": len(unique_study_ids) == 1,
        "result_set_ids_present": not missing_result_set_id,
        "result_set_ids_unique": len(unique_result_set_ids) == 1,
        "source_tool_is_jmag": not bad_source_tool,
        "paths_present": not missing_paths,
        "upstream_gate_status_ok": not bad_upstream_status,
        "upstream_gate_policy_known": not bad_upstream_policy,
        "notebook_operating_point_present": not missing_notebook_operating_point,
        "notebook_operating_point_in_value_table": (
            bool(notebook_operating_points)
            and bool(table_op_set)
            and all(op_id in table_op_set for op_id in notebook_operating_points)
        ),
    }
    if expected_case_id is not None:
        checks["expected_case_id_matches"] = unique_case_ids == [str(expected_case_id)]
    if expected_study_id is not None:
        checks["expected_study_id_matches"] = unique_study_ids == [str(expected_study_id)]
    if expected_result_set_id is not None:
        checks["expected_result_set_id_matches"] = unique_result_set_ids == [str(expected_result_set_id)]

    return {
        "policy": "jmag_export_case_package_gate",
        "required_kinds": list(required),
        "present_kinds": dict(sorted(kind_counts.items())),
        "case_ids": unique_case_ids,
        "study_ids": unique_study_ids,
        "result_set_ids": unique_result_set_ids,
        "table_operating_point_ids": sorted(table_op_set),
        "notebook_operating_point_ids": sorted(set(notebook_operating_points)),
        "expected_case_id": None if expected_case_id is None else str(expected_case_id),
        "expected_study_id": None if expected_study_id is None else str(expected_study_id),
        "expected_result_set_id": None if expected_result_set_id is None else str(expected_result_set_id),
        "missing_case_id_rows": missing_case_id,
        "missing_study_id_rows": missing_study_id,
        "missing_result_set_id_rows": missing_result_set_id,
        "missing_notebook_operating_point_rows": missing_notebook_operating_point,
        "bad_source_tool_rows": bad_source_tool,
        "missing_path_rows": missing_paths,
        "unknown_kinds": unknown_kinds,
        "bad_upstream_status_rows": bad_upstream_status,
        "bad_upstream_policy_rows": bad_upstream_policy,
        "artifacts": details,
        "checks": checks,
        "status": "ok" if all(checks.values()) else "needs_attention",
        "version_note": (
            "Run after JMAG column metadata, symmetry coverage, value-table, "
            "and notebook-row gates so a panel cannot join rows from different "
            "case/study/result-set exports."
        ),
    }


def jmag_current_torque_solver_ready_manifest_gate(
    artifacts,
    expected_case_id=None,
    expected_result_set_id=None,
    expected_operating_point_id=None,
    expected_phases=("U", "V", "W"),
    required_kinds=("column_metadata", "symmetry_coverage", "current_snapshot", "torque_table"),
):
    """Check JMAG current/torque export metadata before value parsing.

    This is the motor-table counterpart of the FEMM pre-solve manifest.  It
    keeps column metadata, symmetry-sector coverage, instantaneous current
    snapshots, and torque-angle tables tied to one JMAG case/result set and one
    operating point before a notebook or validation script joins those rows.
    """

    rows_in = list(artifacts)
    if not rows_in:
        raise ValueError("artifacts must not be empty")

    def _norm(value):
        return str(value or "").strip().lower().replace("-", "_").replace(" ", "_")

    def _first(row, names):
        for name in names:
            if name in row and row[name] is not None:
                return row[name]
        return None

    def _phase_list(value):
        if value is None:
            return []
        if isinstance(value, str):
            return [part.strip() for part in value.replace(";", ",").split(",") if part.strip()]
        try:
            return [str(part).strip() for part in value if str(part).strip()]
        except TypeError:
            return [str(value).strip()]

    required = tuple(_norm(kind) for kind in required_kinds)
    expected_phase_list = [str(phase).strip() for phase in expected_phases]
    if not required:
        raise ValueError("required_kinds must not be empty")
    if not expected_phase_list:
        raise ValueError("expected_phases must not be empty")

    expected_policies = {
        "column_metadata": {"jmag_motor_table_column_metadata_gate"},
        "symmetry_coverage": {"jmag_symmetry_sweep_coverage_gate"},
        "current_snapshot": {
            "motor_current_snapshot_table_contract_gate",
            "spwm_snapshot_current_handoff_gate",
        },
        "torque_table": {
            "torque_angle_table_export_health",
            "torque_angle_sweep_health_summary",
            "dq_torque_table_health",
            "dq_torque_table_schema_health_gate",
            "coenergy_torque_periodic_derivative_gate",
        },
    }

    details = []
    kind_counts = {}
    case_ids = []
    result_set_ids = []
    operating_point_ids = []
    missing_case_id = []
    missing_result_set_id = []
    missing_operating_point_id = []
    bad_source_tool = []
    missing_paths = []
    unknown_kinds = []
    bad_upstream_status = []
    bad_upstream_policy = []
    bad_current_kind = []
    phase_mismatches = []
    bad_torque_lock = []

    for index, row in enumerate(rows_in, start=1):
        if not isinstance(row, dict):
            raise ValueError("each artifact must be a dictionary")
        kind = _norm(_first(row, ("kind", "artifact_kind", "type")))
        case_id = _first(row, ("case_id", "jmag_case_id", "design_case_id"))
        result_set_id = _first(row, ("result_set_id", "result_id", "dataset_id", "export_id"))
        operating_point_id = _first(row, ("operating_point_id", "op_id", "snapshot_id"))
        source_tool = _first(row, ("source_tool", "tool", "source"))
        source_tool_norm = _norm(source_tool)
        path = _first(row, ("path", "file", "table_path", "artifact_path"))
        gate_policy = _first(row, ("gate_policy", "policy", "validator"))
        gate_policy_norm = _norm(gate_policy)
        status = _first(row, ("status", "gate_status", "validation_status"))
        status_norm = _norm(status)
        pass_flag = bool(row.get("pass", False))
        phases = _phase_list(_first(row, ("phase_set", "phases", "phase_order", "circuit_names")))

        if kind:
            kind_counts[kind] = kind_counts.get(kind, 0) + 1
        else:
            unknown_kinds.append(index)
        if kind and kind not in set(required) | set(expected_policies):
            unknown_kinds.append(kind)
        if not case_id:
            missing_case_id.append(index)
        else:
            case_ids.append(str(case_id))
        if not result_set_id:
            missing_result_set_id.append(index)
        else:
            result_set_ids.append(str(result_set_id))
        if kind in {"current_snapshot", "torque_table"}:
            if not operating_point_id:
                missing_operating_point_id.append(index)
            else:
                operating_point_ids.append(str(operating_point_id))
        if source_tool_norm not in {"jmag", "jmag_designer", "jmagdesigner"}:
            bad_source_tool.append({"index": index, "kind": kind, "source_tool": source_tool})
        if not path:
            missing_paths.append(index)
        if not (pass_flag or status_norm in {"ok", "pass", "passed", "verified"}):
            bad_upstream_status.append({"index": index, "kind": kind, "status": status})
        if kind in expected_policies and gate_policy_norm not in expected_policies[kind]:
            bad_upstream_policy.append({
                "index": index,
                "kind": kind,
                "gate_policy": gate_policy,
                "expected": sorted(expected_policies[kind]),
            })
        if kind == "current_snapshot" and _norm(row.get("current_kind")) not in {
            "instantaneous",
            "instant",
            "snapshot",
            "sample",
            "peak",
        }:
            bad_current_kind.append(index)
        if phases and phases != expected_phase_list:
            phase_mismatches.append({"index": index, "kind": kind, "phases": phases})
        if kind == "torque_table" and row.get("rotor_current_phase_locked") is not True:
            bad_torque_lock.append(index)

        details.append({
            "index": index,
            "kind": kind,
            "case_id": None if case_id is None else str(case_id),
            "result_set_id": None if result_set_id is None else str(result_set_id),
            "operating_point_id": None if operating_point_id is None else str(operating_point_id),
            "source_tool": source_tool,
            "path": path,
            "gate_policy": gate_policy,
            "status": status,
            "pass": pass_flag,
            "phases": phases,
            "rotor_current_phase_locked": row.get("rotor_current_phase_locked"),
        })

    required_set = set(required)
    present_set = set(kind_counts)
    unique_case_ids = sorted(set(case_ids))
    unique_result_set_ids = sorted(set(result_set_ids))
    unique_operating_point_ids = sorted(set(operating_point_ids))
    checks = {
        "required_kinds_present": required_set.issubset(present_set),
        "no_unknown_kinds": not unknown_kinds,
        "case_ids_present": not missing_case_id,
        "case_ids_unique": len(unique_case_ids) == 1,
        "result_set_ids_present": not missing_result_set_id,
        "result_set_ids_unique": len(unique_result_set_ids) == 1,
        "source_tool_is_jmag": not bad_source_tool,
        "paths_present": not missing_paths,
        "operating_point_ids_present_for_current_and_torque": not missing_operating_point_id,
        "operating_point_ids_unique": len(unique_operating_point_ids) == 1,
        "upstream_gate_status_ok": not bad_upstream_status,
        "upstream_gate_policy_known": not bad_upstream_policy,
        "current_snapshot_is_instantaneous": not bad_current_kind,
        "phase_sets_match_expected": not phase_mismatches,
        "torque_table_locked_to_current_phase": not bad_torque_lock,
    }
    if expected_case_id is not None:
        checks["expected_case_id_matches"] = unique_case_ids == [str(expected_case_id)]
    if expected_result_set_id is not None:
        checks["expected_result_set_id_matches"] = unique_result_set_ids == [str(expected_result_set_id)]
    if expected_operating_point_id is not None:
        checks["expected_operating_point_id_matches"] = unique_operating_point_ids == [str(expected_operating_point_id)]

    return {
        "policy": "jmag_current_torque_solver_ready_manifest_gate",
        "required_kinds": list(required),
        "present_kinds": dict(sorted(kind_counts.items())),
        "case_ids": unique_case_ids,
        "result_set_ids": unique_result_set_ids,
        "operating_point_ids": unique_operating_point_ids,
        "expected_case_id": None if expected_case_id is None else str(expected_case_id),
        "expected_result_set_id": None if expected_result_set_id is None else str(expected_result_set_id),
        "expected_operating_point_id": None if expected_operating_point_id is None else str(expected_operating_point_id),
        "expected_phases": expected_phase_list,
        "missing_case_id_rows": missing_case_id,
        "missing_result_set_id_rows": missing_result_set_id,
        "missing_operating_point_id_rows": missing_operating_point_id,
        "bad_source_tool_rows": bad_source_tool,
        "missing_path_rows": missing_paths,
        "unknown_kinds": unknown_kinds,
        "bad_upstream_status_rows": bad_upstream_status,
        "bad_upstream_policy_rows": bad_upstream_policy,
        "bad_current_kind_rows": bad_current_kind,
        "phase_mismatch_rows": phase_mismatches,
        "bad_torque_lock_rows": bad_torque_lock,
        "artifacts": details,
        "checks": checks,
        "status": "ok" if all(checks.values()) else "needs_attention",
        "version_note": (
            "Run after JMAG column metadata, symmetry coverage, current-snapshot, "
            "and torque-table gates, before table values are joined, so a "
            "torque row cannot drift away from the current phase or result set."
        ),
    }


def jmag_efficiency_operating_point_package_gate(
    artifacts,
    expected_case_id=None,
    expected_result_set_id=None,
    required_kinds=("terminal_table", "loss_bucket_table", "drive_cycle_weights", "notebook_row"),
):
    """Check that JMAG efficiency-map artifacts share operating-point ids.

    Terminal dq rows, loss-bucket tables, drive-cycle weights, and selected
    notebook rows are often exported or curated separately.  This gate keeps
    their operating-point sets aligned before an efficiency map or panel
    compares maxima, losses, or drive-cycle scores.
    """

    rows_in = list(artifacts)
    if not rows_in:
        raise ValueError("artifacts must not be empty")

    def _norm(value):
        return str(value or "").strip().lower().replace("-", "_").replace(" ", "_")

    def _first(row, names):
        for name in names:
            if name in row and row[name] is not None:
                return row[name]
        return None

    def _op_ids(row):
        values = row.get("operating_point_ids", row.get("point_ids"))
        if values is None:
            single = _first(row, ("operating_point_id", "point_id", "op_id"))
            return [] if single is None else [str(single)]
        if isinstance(values, str):
            return [part.strip() for part in values.replace(";", ",").split(",") if part.strip()]
        try:
            return [str(value).strip() for value in values if str(value).strip()]
        except TypeError:
            return [str(values).strip()]

    required = tuple(_norm(kind) for kind in required_kinds)
    if not required:
        raise ValueError("required_kinds must not be empty")

    expected_policies = {
        "terminal_table": {"pm_drive_terminal_table_health", "pm_drive_terminal_table_health_gate"},
        "loss_bucket_table": {"pm_drive_loss_bucket_efficiency_gate"},
        "drive_cycle_weights": {"drive_cycle_weighted_efficiency_gate"},
        "notebook_row": {"pm_drive_operating_point_notebook_handoff_gate"},
    }

    details = []
    kind_counts = {}
    case_ids = []
    result_set_ids = []
    op_sets = {}
    selected_notebook_ops = []
    missing_case_id = []
    missing_result_set_id = []
    missing_operating_points = []
    bad_source_tool = []
    missing_paths = []
    unknown_kinds = []
    bad_upstream_status = []
    bad_upstream_policy = []

    for index, row in enumerate(rows_in, start=1):
        if not isinstance(row, dict):
            raise ValueError("each artifact must be a dictionary")
        kind = _norm(_first(row, ("kind", "artifact_kind", "type")))
        case_id = _first(row, ("case_id", "jmag_case_id", "design_case_id"))
        result_set_id = _first(row, ("result_set_id", "result_id", "dataset_id", "export_id"))
        source_tool = _first(row, ("source_tool", "tool", "source"))
        source_tool_norm = _norm(source_tool)
        path = _first(row, ("path", "file", "table_path", "artifact_path"))
        gate_policy = _first(row, ("gate_policy", "policy", "validator"))
        gate_policy_norm = _norm(gate_policy)
        status = _first(row, ("status", "gate_status", "validation_status"))
        status_norm = _norm(status)
        pass_flag = bool(row.get("pass", False))
        ops = _op_ids(row)

        if kind:
            kind_counts[kind] = kind_counts.get(kind, 0) + 1
        else:
            unknown_kinds.append(index)
        if kind and kind not in set(required) | set(expected_policies):
            unknown_kinds.append(kind)
        if not case_id:
            missing_case_id.append(index)
        else:
            case_ids.append(str(case_id))
        if not result_set_id:
            missing_result_set_id.append(index)
        else:
            result_set_ids.append(str(result_set_id))
        if source_tool_norm not in {"jmag", "jmag_designer", "jmagdesigner"}:
            bad_source_tool.append({"index": index, "source_tool": source_tool})
        if not path:
            missing_paths.append(index)
        if not ops:
            missing_operating_points.append(index)
        else:
            op_sets[kind] = set(ops)
            if kind == "notebook_row":
                selected_notebook_ops.extend(ops)
        if not (pass_flag or status_norm in {"ok", "pass", "passed", "verified"}):
            bad_upstream_status.append({"index": index, "kind": kind, "status": status})
        if kind in expected_policies and gate_policy_norm not in expected_policies[kind]:
            bad_upstream_policy.append({
                "index": index,
                "kind": kind,
                "gate_policy": gate_policy,
                "expected": sorted(expected_policies[kind]),
            })
        details.append({
            "index": index,
            "kind": kind,
            "case_id": None if case_id is None else str(case_id),
            "result_set_id": None if result_set_id is None else str(result_set_id),
            "source_tool": source_tool,
            "path": path,
            "gate_policy": gate_policy,
            "status": status,
            "pass": pass_flag,
            "operating_point_ids": ops,
        })

    required_set = set(required)
    present_set = set(kind_counts)
    unique_case_ids = sorted(set(case_ids))
    unique_result_set_ids = sorted(set(result_set_ids))
    table_kinds = [kind for kind in ("terminal_table", "loss_bucket_table", "drive_cycle_weights") if kind in op_sets]
    reference_ops = op_sets.get("terminal_table") or (op_sets[table_kinds[0]] if table_kinds else set())
    mismatched_op_sets = [
        {"kind": kind, "operating_point_ids": sorted(op_sets[kind])}
        for kind in table_kinds
        if op_sets[kind] != reference_ops
    ]
    notebook_op_set = set(selected_notebook_ops)
    checks = {
        "required_kinds_present": required_set.issubset(present_set),
        "no_unknown_kinds": not unknown_kinds,
        "case_ids_present": not missing_case_id,
        "case_ids_unique": len(unique_case_ids) == 1,
        "result_set_ids_present": not missing_result_set_id,
        "result_set_ids_unique": len(unique_result_set_ids) == 1,
        "source_tool_is_jmag": not bad_source_tool,
        "paths_present": not missing_paths,
        "operating_point_ids_present": not missing_operating_points,
        "table_operating_point_sets_match": bool(reference_ops) and not mismatched_op_sets,
        "notebook_operating_points_in_tables": bool(notebook_op_set) and notebook_op_set.issubset(reference_ops),
        "upstream_gate_status_ok": not bad_upstream_status,
        "upstream_gate_policy_known": not bad_upstream_policy,
    }
    if expected_case_id is not None:
        checks["expected_case_id_matches"] = unique_case_ids == [str(expected_case_id)]
    if expected_result_set_id is not None:
        checks["expected_result_set_id_matches"] = unique_result_set_ids == [str(expected_result_set_id)]

    return {
        "policy": "jmag_efficiency_operating_point_package_gate",
        "required_kinds": list(required),
        "present_kinds": dict(sorted(kind_counts.items())),
        "case_ids": unique_case_ids,
        "result_set_ids": unique_result_set_ids,
        "expected_case_id": None if expected_case_id is None else str(expected_case_id),
        "expected_result_set_id": None if expected_result_set_id is None else str(expected_result_set_id),
        "reference_operating_point_ids": sorted(reference_ops),
        "notebook_operating_point_ids": sorted(notebook_op_set),
        "mismatched_operating_point_sets": mismatched_op_sets,
        "missing_case_id_rows": missing_case_id,
        "missing_result_set_id_rows": missing_result_set_id,
        "missing_operating_point_rows": missing_operating_points,
        "bad_source_tool_rows": bad_source_tool,
        "missing_path_rows": missing_paths,
        "unknown_kinds": unknown_kinds,
        "bad_upstream_status_rows": bad_upstream_status,
        "bad_upstream_policy_rows": bad_upstream_policy,
        "artifacts": details,
        "checks": checks,
        "status": "ok" if all(checks.values()) else "needs_attention",
        "version_note": (
            "Run after terminal, loss-bucket, drive-cycle, and notebook-row "
            "gates so a JMAG efficiency map cannot mix operating-point sets "
            "from different exports or curated tables."
        ),
    }


def spwm_snapshot_current_handoff_summary(
    id_current,
    iq_current,
    sample_count,
    theta0_rad=0.0,
    sample_offset_fraction=0.5,
    sampling_mode=None,
    timer_alignment=None,
    carrier_ratio=None,
    phases=("U", "V", "W"),
    tol=1.0e-12,
):
    """Build and check a balanced current snapshot table for FEMM/static solves.

    SPWM and timer details live in the drive/control layer, but each FEM
    snapshot still receives one balanced three-phase current row.  This helper
    samples one electrical period, converts dq -> U/V/W, immediately recovers
    dq, and checks zero-sequence, RMS, and phase-order-sensitive round-trip
    invariants before a solver-specific circuit-current API is called.
    """

    n = int(sample_count)
    if n < 3:
        raise ValueError("sample_count must be at least 3")
    d = float(id_current)
    q = float(iq_current)
    offset = float(sample_offset_fraction)
    theta0 = float(theta0_rad)
    tolerance = float(tol)
    if tolerance < 0.0:
        raise ValueError("tol must be non-negative")
    sampling_aliases = {
        "symmetrical": "symmetrical",
        "symmetric": "symmetrical",
        "center_aligned": "symmetrical",
        "centered": "symmetrical",
        "asymmetrical": "asymmetrical",
        "asymmetric": "asymmetrical",
        "edge_aligned": "asymmetrical",
        "edge": "asymmetrical",
    }
    sampling = None
    if sampling_mode is not None:
        key = str(sampling_mode).strip().lower().replace("-", "_").replace(" ", "_")
        if key not in sampling_aliases:
            raise ValueError("sampling_mode must be symmetrical/asymmetrical or a center/edge alias")
        sampling = sampling_aliases[key]
    timer = None if timer_alignment is None else str(timer_alignment).strip().lower().replace("-", "_").replace(" ", "_")
    carrier = None if carrier_ratio is None else float(carrier_ratio)
    if carrier is not None and carrier <= 0.0:
        raise ValueError("carrier_ratio must be positive when provided")

    rows = []
    id_errors = []
    iq_errors = []
    zero_sequence = []
    square_errors = []
    phase_squares = {phase: [] for phase in phases}
    for index in range(n):
        theta = theta0 + 2.0 * math.pi * (index + offset) / n
        currents = dq_to_three_phase_currents(d, q, theta, phases=phases)
        dq = three_phase_currents_to_dq_summary(
            currents,
            theta,
            expected_id=d,
            expected_iq=q,
            phases=phases,
            tol=tolerance,
        )
        id_errors.append(dq["id_abs_error"])
        iq_errors.append(dq["iq_abs_error"])
        zero_sequence.append(dq["zero_sequence_abs"])
        square_errors.append(dq["abc_vs_dq_square_sum_error"])
        for phase in phases:
            phase_squares[phase].append(currents[phase] * currents[phase])
        rows.append({
            "sample": index,
            "theta_e_rad": theta,
            "theta_e_deg": math.degrees(theta),
            "currents": currents,
            "dq": {key: dq[key] for key in ("id", "iq", "i0", "status")},
        })

    amplitude = math.hypot(d, q)
    expected_phase_rms = amplitude / math.sqrt(2.0)
    phase_rms = {
        phase: math.sqrt(sum(values) / n)
        for phase, values in phase_squares.items()
    }
    phase_rms_errors = {
        phase: abs(value - expected_phase_rms)
        for phase, value in phase_rms.items()
    }
    checks = {
        "dq_roundtrip_ok": max(id_errors + iq_errors) <= tolerance,
        "zero_sequence_ok": max(zero_sequence) <= tolerance,
        "abc_square_sum_ok": max(square_errors) <= tolerance,
        "phase_rms_ok": max(phase_rms_errors.values()) <= max(tolerance, 1.0e-12 * max(1.0, expected_phase_rms)),
    }
    expected_offset_for_sampling = None
    if sampling == "symmetrical":
        expected_offset_for_sampling = 0.5
    elif sampling == "asymmetrical":
        expected_offset_for_sampling = 0.0
    if expected_offset_for_sampling is not None:
        checks["sampling_offset_matches_mode"] = (
            abs(offset - expected_offset_for_sampling) <= max(tolerance, 1.0e-12)
        )
    return {
        "policy": "spwm_snapshot_current_handoff_gate",
        "id": d,
        "iq": q,
        "current_amplitude": amplitude,
        "sample_count": n,
        "theta0_rad": theta0,
        "sample_offset_fraction": offset,
        "sampling_mode": sampling,
        "timer_alignment": timer,
        "carrier_ratio": carrier,
        "phase_order": list(phases),
        "expected_phase_rms": expected_phase_rms,
        "phase_rms": phase_rms,
        "max_phase_rms_abs_error": max(phase_rms_errors.values()),
        "max_id_abs_error": max(id_errors),
        "max_iq_abs_error": max(iq_errors),
        "max_zero_sequence_abs": max(zero_sequence),
        "max_abc_square_sum_error": max(square_errors),
        "checks": checks,
        "rows": rows,
        "status": "ok" if all(checks.values()) else "needs_attention",
    }


def motor_current_snapshot_table_contract_gate(
    rows,
    pole_pairs,
    phases=("U", "V", "W"),
    expected_current_kind="instantaneous",
    require_sampling_metadata=True,
    tol=1.0e-12,
):
    """Check solver-ready motor current snapshot table rows.

    This is the public-safe table contract that sits between drive/control
    sampling rows and any FEM motor map export.  It keeps electrical angle,
    mechanical angle, pole pairs, U/V/W currents, dq recovery, current-angle
    gamma, current magnitude, and sampling metadata in the same row so that a
    table is not promoted merely because the three-phase row is balanced.
    """

    table = list(rows)
    if not table:
        raise ValueError("rows must not be empty")
    phase_names = tuple(phases)
    if len(phase_names) != 3:
        raise ValueError("exactly three phase names are required")
    pp = float(pole_pairs)
    tolerance = float(tol)
    if pp <= 0.0:
        raise ValueError("pole_pairs must be > 0")
    if tolerance < 0.0:
        raise ValueError("tol must be non-negative")

    def _normalize_kind(kind):
        key = str(kind).strip().lower().replace("-", "_").replace(" ", "_")
        aliases = {
            "instantaneous": "instantaneous",
            "instant": "instantaneous",
            "snapshot": "instantaneous",
            "sample": "instantaneous",
            "sampled": "instantaneous",
            "peak": "instantaneous",
            "peak_sample": "instantaneous",
            "rms": "rms",
            "root_mean_square": "rms",
        }
        return aliases.get(key, key)

    def _angle_rad(row, rad_key, deg_key):
        if rad_key in row:
            return float(row[rad_key])
        if deg_key in row:
            return math.radians(float(row[deg_key]))
        return None

    expected_kind = _normalize_kind(expected_current_kind)
    required_sampling = ("sampling_mode", "timer_alignment", "carrier_ratio", "sample_offset_fraction")
    row_summaries = []
    sample_indices = []
    max_angle_error = 0.0
    max_id_error = 0.0
    max_iq_error = 0.0
    max_gamma_error = 0.0
    max_current_error = 0.0
    max_zero_sequence = 0.0
    checks = {
        "sample_index_present": True,
        "angle_columns_present": True,
        "mechanical_electrical_angle_ok": True,
        "current_columns_present": True,
        "dq_columns_present": True,
        "dq_recovery_ok": True,
        "gamma_column_ok": True,
        "current_magnitude_ok": True,
        "current_kind_matches": True,
        "sampling_metadata_present": True,
    }
    for index, row in enumerate(table):
        sample_index = row.get("sample_index", row.get("sample"))
        if sample_index is None:
            checks["sample_index_present"] = False
        else:
            sample_indices.append(int(sample_index))
        theta_e = _angle_rad(row, "theta_e_rad", "theta_e_deg")
        theta_m = _angle_rad(row, "theta_mech_rad", "theta_mech_deg")
        if theta_e is None or theta_m is None:
            checks["angle_columns_present"] = False
            theta_e = 0.0 if theta_e is None else theta_e
            theta_m = 0.0 if theta_m is None else theta_m
        angle_error = abs(theta_e - pp * theta_m)
        max_angle_error = max(max_angle_error, angle_error)
        if angle_error > tolerance:
            checks["mechanical_electrical_angle_ok"] = False

        current_keys = [f"current_{phase}_A" for phase in phase_names]
        if not all(key in row for key in current_keys):
            checks["current_columns_present"] = False
            currents = {phase: 0.0 for phase in phase_names}
        else:
            currents = {
                phase: float(row[f"current_{phase}_A"])
                for phase in phase_names
            }
        if "id_A" not in row or "iq_A" not in row:
            checks["dq_columns_present"] = False
            id_row = 0.0
            iq_row = 0.0
        else:
            id_row = float(row["id_A"])
            iq_row = float(row["iq_A"])
        dq = three_phase_currents_to_dq_summary(
            currents,
            theta_e,
            expected_id=id_row,
            expected_iq=iq_row,
            phases=phase_names,
            tol=tolerance,
        )
        max_id_error = max(max_id_error, dq.get("id_abs_error", 0.0))
        max_iq_error = max(max_iq_error, dq.get("iq_abs_error", 0.0))
        max_zero_sequence = max(max_zero_sequence, dq["zero_sequence_abs"])
        if dq["status"] != "ok":
            checks["dq_recovery_ok"] = False

        gamma_expected = math.degrees(math.atan2(-id_row, iq_row))
        gamma_error = None
        if "gamma_deg" in row:
            gamma_error = abs(float(row["gamma_deg"]) - gamma_expected)
            max_gamma_error = max(max_gamma_error, gamma_error)
            if gamma_error > max(tolerance, 1.0e-9):
                checks["gamma_column_ok"] = False
        else:
            checks["gamma_column_ok"] = False

        current_expected = math.hypot(id_row, iq_row)
        current_error = None
        if "current_A" in row:
            current_error = abs(float(row["current_A"]) - current_expected)
            max_current_error = max(max_current_error, current_error)
            if current_error > max(tolerance, 1.0e-12 * max(1.0, current_expected)):
                checks["current_magnitude_ok"] = False
        else:
            checks["current_magnitude_ok"] = False

        kind = row.get("current_kind")
        normalized_kind = None if kind is None else _normalize_kind(kind)
        if normalized_kind != expected_kind:
            checks["current_kind_matches"] = False
        missing_sampling = [key for key in required_sampling if key not in row]
        if require_sampling_metadata and missing_sampling:
            checks["sampling_metadata_present"] = False
        row_summaries.append({
            "row": index,
            "sample_index": sample_index,
            "theta_e_rad": theta_e,
            "theta_mech_rad": theta_m,
            "angle_abs_error_rad": angle_error,
            "dq": {key: dq[key] for key in ("id", "iq", "i0", "status")},
            "id_abs_error_A": dq.get("id_abs_error"),
            "iq_abs_error_A": dq.get("iq_abs_error"),
            "gamma_expected_deg": gamma_expected,
            "gamma_abs_error_deg": gamma_error,
            "current_expected_A": current_expected,
            "current_abs_error_A": current_error,
            "current_kind": kind,
            "normalized_current_kind": normalized_kind,
            "missing_sampling_metadata": missing_sampling,
        })

    checks["sample_index_monotone"] = all(a < b for a, b in zip(sample_indices, sample_indices[1:]))
    return {
        "policy": "motor_current_snapshot_table_contract_gate",
        "pole_pairs": pp,
        "phase_order": list(phase_names),
        "expected_current_kind": expected_kind,
        "require_sampling_metadata": bool(require_sampling_metadata),
        "n_rows": len(row_summaries),
        "max_angle_abs_error_rad": max_angle_error,
        "max_id_abs_error_A": max_id_error,
        "max_iq_abs_error_A": max_iq_error,
        "max_zero_sequence_abs_A": max_zero_sequence,
        "max_gamma_abs_error_deg": max_gamma_error,
        "max_current_abs_error_A": max_current_error,
        "checks": checks,
        "rows": row_summaries,
        "status": "ok" if all(checks.values()) else "needs_attention",
    }


def femm_motor_model_artifact_package_gate(
    artifacts,
    expected_model_id=None,
    expected_operating_point_id=None,
    required_kinds=("block_labels", "current_snapshot", "torque_table"),
):
    """Check that FEMM-derived motor artifacts belong to one model snapshot.

    Block labels, current snapshots, and torque-angle tables are often exported
    by different scripts.  This package gate keeps the model identity and the
    operating-point identity explicit before the rows are promoted to
    radia-ngsolve or notebook examples.
    """

    rows_in = list(artifacts)
    if not rows_in:
        raise ValueError("artifacts must not be empty")

    def _norm(value):
        return str(value or "").strip().lower().replace("-", "_").replace(" ", "_")

    def _first(row, names):
        for name in names:
            if name in row and row[name] is not None:
                return row[name]
        return None

    required = tuple(_norm(kind) for kind in required_kinds)
    if not required:
        raise ValueError("required_kinds must not be empty")

    expected_policies = {
        "block_labels": {"femm_block_label_source_contract_gate"},
        "current_snapshot": {
            "motor_current_snapshot_table_contract_gate",
            "femm_static_current_circuit_rows_gate",
            "spwm_snapshot_current_handoff_gate",
        },
        "torque_table": {
            "torque_angle_table_export_health",
            "torque_angle_sweep_health_summary",
        },
    }

    details = []
    kind_counts = {}
    model_ids = []
    operating_point_ids = []
    missing_model_id = []
    missing_operating_point_id = []
    bad_source_tool = []
    missing_paths = []
    unknown_kinds = []
    bad_upstream_status = []
    bad_upstream_policy = []
    bad_current_kind = []
    bad_torque_metadata = []

    for index, row in enumerate(rows_in, start=1):
        if not isinstance(row, dict):
            raise ValueError("each artifact must be a dictionary")
        kind = _norm(_first(row, ("kind", "artifact_kind", "type")))
        model_id = _first(row, ("model_id", "motor_model_id", "geometry_id"))
        operating_point_id = _first(row, ("operating_point_id", "op_id", "snapshot_id"))
        source_tool = _first(row, ("source_tool", "tool", "source"))
        source_tool_norm = _norm(source_tool)
        path = _first(row, ("path", "file", "table_path", "artifact_path"))
        gate_policy = _first(row, ("gate_policy", "policy", "validator"))
        gate_policy_norm = _norm(gate_policy)
        status = _first(row, ("status", "gate_status", "validation_status"))
        status_norm = _norm(status)
        pass_flag = bool(row.get("pass", False))

        if kind:
            kind_counts[kind] = kind_counts.get(kind, 0) + 1
        else:
            unknown_kinds.append(index)
        if kind and kind not in set(required) | set(expected_policies):
            unknown_kinds.append(kind)
        if not model_id:
            missing_model_id.append(index)
        else:
            model_ids.append(str(model_id))
        if kind in {"current_snapshot", "torque_table"}:
            if not operating_point_id:
                missing_operating_point_id.append(index)
            else:
                operating_point_ids.append(str(operating_point_id))
        if source_tool_norm not in {"femm", "pyfemm"}:
            bad_source_tool.append({"index": index, "source_tool": source_tool})
        if not path:
            missing_paths.append(index)
        if not (pass_flag or status_norm in {"ok", "pass", "passed", "verified"}):
            bad_upstream_status.append({"index": index, "kind": kind, "status": status})
        if kind in expected_policies and gate_policy_norm not in expected_policies[kind]:
            bad_upstream_policy.append({
                "index": index,
                "kind": kind,
                "gate_policy": gate_policy,
                "expected": sorted(expected_policies[kind]),
            })
        if kind == "current_snapshot" and _norm(row.get("current_kind")) not in {"instantaneous", "instant", "snapshot", "sample", "peak"}:
            bad_current_kind.append(index)
        if kind == "torque_table":
            angle_basis_ok = _norm(row.get("angle_basis")) == "mechanical"
            source_function_ok = str(row.get("source_function") or "").strip() == "mo_blockintegral(22)"
            locked = row.get("rotor_current_phase_locked")
            locked_ok = locked is True or _norm(locked) == "true"
            if not (angle_basis_ok and source_function_ok and locked_ok):
                bad_torque_metadata.append({
                    "index": index,
                    "angle_basis": row.get("angle_basis"),
                    "source_function": row.get("source_function"),
                    "rotor_current_phase_locked": row.get("rotor_current_phase_locked"),
                })
        details.append({
            "index": index,
            "kind": kind,
            "model_id": None if model_id is None else str(model_id),
            "operating_point_id": None if operating_point_id is None else str(operating_point_id),
            "source_tool": source_tool,
            "path": path,
            "gate_policy": gate_policy,
            "status": status,
            "pass": pass_flag,
        })

    required_set = set(required)
    present_set = set(kind_counts)
    unique_model_ids = sorted(set(model_ids))
    unique_operating_point_ids = sorted(set(operating_point_ids))
    checks = {
        "required_kinds_present": required_set.issubset(present_set),
        "no_unknown_kinds": not unknown_kinds,
        "model_ids_present": not missing_model_id,
        "model_ids_unique": len(unique_model_ids) == 1,
        "source_tool_is_femm": not bad_source_tool,
        "paths_present": not missing_paths,
        "operating_point_ids_present_for_current_and_torque": not missing_operating_point_id,
        "operating_point_ids_unique": len(unique_operating_point_ids) == 1,
        "upstream_gate_status_ok": not bad_upstream_status,
        "upstream_gate_policy_known": not bad_upstream_policy,
        "current_snapshot_is_instantaneous": not bad_current_kind,
        "torque_table_metadata_solver_ready": not bad_torque_metadata,
    }
    if expected_model_id is not None:
        checks["expected_model_id_matches"] = unique_model_ids == [str(expected_model_id)]
    if expected_operating_point_id is not None:
        checks["expected_operating_point_id_matches"] = unique_operating_point_ids == [str(expected_operating_point_id)]

    return {
        "policy": "femm_motor_model_artifact_package_gate",
        "required_kinds": list(required),
        "present_kinds": dict(sorted(kind_counts.items())),
        "model_ids": unique_model_ids,
        "operating_point_ids": unique_operating_point_ids,
        "expected_model_id": None if expected_model_id is None else str(expected_model_id),
        "expected_operating_point_id": None if expected_operating_point_id is None else str(expected_operating_point_id),
        "missing_model_id_rows": missing_model_id,
        "missing_operating_point_id_rows": missing_operating_point_id,
        "bad_source_tool_rows": bad_source_tool,
        "missing_path_rows": missing_paths,
        "unknown_kinds": unknown_kinds,
        "bad_upstream_status_rows": bad_upstream_status,
        "bad_upstream_policy_rows": bad_upstream_policy,
        "bad_current_kind_rows": bad_current_kind,
        "bad_torque_metadata_rows": bad_torque_metadata,
        "artifacts": details,
        "checks": checks,
        "status": "ok" if all(checks.values()) else "needs_attention",
        "version_note": (
            "Run after FEMM block-label, current-snapshot, and torque-table "
            "gates so a notebook or radia-ngsolve handoff cannot mix rows from "
            "different motor models or operating points."
        ),
    }


def femm_winding_current_package_gate(
    artifacts,
    expected_model_id=None,
    expected_phases=("U", "V", "W"),
    required_kinds=("winding_table", "block_labels", "current_snapshot"),
):
    """Check that winding, FEMM source labels, and currents describe one motor.

    A winding-factor table is usually computed before FEMM solves, while block
    labels and circuit currents are exported by FEMM/pyFEMM scripts.  This gate
    keeps those rows tied to one ``model_id`` and one phase set before the
    package is promoted to torque, flux-linkage, or back-EMF notebooks.
    """

    rows_in = list(artifacts)
    if not rows_in:
        raise ValueError("artifacts must not be empty")

    def _norm(value):
        return str(value or "").strip().lower().replace("-", "_").replace(" ", "_")

    def _first(row, names):
        for name in names:
            if name in row and row[name] is not None:
                return row[name]
        return None

    def _phase_list(value):
        if value is None:
            return []
        if isinstance(value, str):
            return [part.strip() for part in value.replace(";", ",").split(",") if part.strip()]
        try:
            return [str(part).strip() for part in value if str(part).strip()]
        except TypeError:
            return [str(value).strip()]

    required = tuple(_norm(kind) for kind in required_kinds)
    expected_phase_list = [str(phase).strip() for phase in expected_phases]
    if not required:
        raise ValueError("required_kinds must not be empty")
    if len(expected_phase_list) < 1:
        raise ValueError("expected_phases must not be empty")

    expected_policies = {
        "winding_table": {
            "double_layer_winding_pitch_harmonic_gate",
            "integral_slot_winding_factor",
            "winding_factor_table_gate",
        },
        "block_labels": {"femm_block_label_source_contract_gate"},
        "current_snapshot": {
            "motor_current_snapshot_table_contract_gate",
            "femm_static_current_circuit_rows_gate",
            "spwm_snapshot_current_handoff_gate",
        },
    }
    source_tools = {
        "winding_table": {"analytic", "radia_ngsolve", "radia-ngsolve", "femm", "pyfemm"},
        "block_labels": {"femm", "pyfemm"},
        "current_snapshot": {"femm", "pyfemm"},
    }

    details = []
    kind_counts = {}
    model_ids = []
    missing_model_id = []
    bad_source_tool = []
    missing_paths = []
    unknown_kinds = []
    bad_upstream_status = []
    bad_upstream_policy = []
    phase_mismatches = []
    bad_winding_metadata = []
    bad_current_kind = []

    for index, row in enumerate(rows_in, start=1):
        if not isinstance(row, dict):
            raise ValueError("each artifact must be a dictionary")
        kind = _norm(_first(row, ("kind", "artifact_kind", "type")))
        model_id = _first(row, ("model_id", "motor_model_id", "geometry_id"))
        source_tool = _first(row, ("source_tool", "tool", "source"))
        source_tool_norm = _norm(source_tool)
        path = _first(row, ("path", "file", "table_path", "artifact_path"))
        gate_policy = _first(row, ("gate_policy", "policy", "validator"))
        gate_policy_norm = _norm(gate_policy)
        status = _first(row, ("status", "gate_status", "validation_status"))
        status_norm = _norm(status)
        pass_flag = bool(row.get("pass", False))
        phases = _phase_list(_first(row, ("phase_set", "phases", "phase_order", "circuit_names")))

        if kind:
            kind_counts[kind] = kind_counts.get(kind, 0) + 1
        else:
            unknown_kinds.append(index)
        if kind and kind not in set(required) | set(expected_policies):
            unknown_kinds.append(kind)
        if not model_id:
            missing_model_id.append(index)
        else:
            model_ids.append(str(model_id))
        if source_tool_norm not in source_tools.get(kind, set()):
            bad_source_tool.append({"index": index, "kind": kind, "source_tool": source_tool})
        if not path:
            missing_paths.append(index)
        if not (pass_flag or status_norm in {"ok", "pass", "passed", "verified"}):
            bad_upstream_status.append({"index": index, "kind": kind, "status": status})
        if kind in expected_policies and gate_policy_norm not in expected_policies[kind]:
            bad_upstream_policy.append({
                "index": index,
                "kind": kind,
                "gate_policy": gate_policy,
                "expected": sorted(expected_policies[kind]),
            })
        if phases and phases != expected_phase_list:
            phase_mismatches.append({"index": index, "kind": kind, "phases": phases})
        if kind == "winding_table":
            slots = row.get("slots", row.get("slot_count"))
            poles = row.get("poles", row.get("pole_count"))
            phase_count = row.get("phase_count", row.get("n_phases", len(phases) if phases else None))
            try:
                metadata_ok = int(slots) > 0 and int(poles) > 0 and int(phase_count) == len(expected_phase_list)
            except (TypeError, ValueError):
                metadata_ok = False
            if not metadata_ok:
                bad_winding_metadata.append(index)
        if kind == "current_snapshot" and _norm(row.get("current_kind")) not in {
            "instantaneous",
            "instant",
            "snapshot",
            "sample",
            "peak",
        }:
            bad_current_kind.append(index)

        details.append({
            "index": index,
            "kind": kind,
            "model_id": None if model_id is None else str(model_id),
            "source_tool": source_tool,
            "path": path,
            "gate_policy": gate_policy,
            "status": status,
            "pass": pass_flag,
            "phases": phases,
        })

    required_set = set(required)
    present_set = set(kind_counts)
    unique_model_ids = sorted(set(model_ids))
    checks = {
        "required_kinds_present": required_set.issubset(present_set),
        "no_unknown_kinds": not unknown_kinds,
        "model_ids_present": not missing_model_id,
        "model_ids_unique": len(unique_model_ids) == 1,
        "source_tools_match_kind": not bad_source_tool,
        "paths_present": not missing_paths,
        "upstream_gate_status_ok": not bad_upstream_status,
        "upstream_gate_policy_known": not bad_upstream_policy,
        "phase_sets_match_expected": not phase_mismatches,
        "winding_geometry_metadata_present": not bad_winding_metadata,
        "current_snapshot_is_instantaneous": not bad_current_kind,
    }
    if expected_model_id is not None:
        checks["expected_model_id_matches"] = unique_model_ids == [str(expected_model_id)]

    return {
        "policy": "femm_winding_current_package_gate",
        "required_kinds": list(required),
        "present_kinds": dict(sorted(kind_counts.items())),
        "model_ids": unique_model_ids,
        "expected_model_id": None if expected_model_id is None else str(expected_model_id),
        "expected_phases": expected_phase_list,
        "missing_model_id_rows": missing_model_id,
        "bad_source_tool_rows": bad_source_tool,
        "missing_path_rows": missing_paths,
        "unknown_kinds": unknown_kinds,
        "bad_upstream_status_rows": bad_upstream_status,
        "bad_upstream_policy_rows": bad_upstream_policy,
        "phase_mismatch_rows": phase_mismatches,
        "bad_winding_metadata_rows": bad_winding_metadata,
        "bad_current_kind_rows": bad_current_kind,
        "artifacts": details,
        "checks": checks,
        "status": "ok" if all(checks.values()) else "needs_attention",
        "version_note": (
            "Run before FEMM torque, flux-linkage, or back-EMF tables so a "
            "winding-factor table cannot be mixed with block labels or currents "
            "from another motor or phase convention."
        ),
    }


def femm_source_current_solver_ready_manifest_gate(
    artifacts,
    expected_model_id=None,
    expected_operating_point_id=None,
    expected_phases=("U", "V", "W"),
    required_kinds=("block_labels", "pm_magnetization", "current_snapshot"),
):
    """Check FEMM source/current metadata before a static solve is trusted.

    This is the pre-solve companion to
    :func:`femm_motor_model_artifact_package_gate`.  It bundles the upstream
    source contracts that must already be true before FEMM/pyFEMM starts
    solving: block-label source rows, PM magnetization convention rows, and
    instantaneous circuit-current rows.
    """

    rows_in = list(artifacts)
    if not rows_in:
        raise ValueError("artifacts must not be empty")

    def _norm(value):
        return str(value or "").strip().lower().replace("-", "_").replace(" ", "_")

    def _first(row, names):
        for name in names:
            if name in row and row[name] is not None:
                return row[name]
        return None

    def _phase_list(value):
        if value is None:
            return []
        if isinstance(value, str):
            return [part.strip() for part in value.replace(";", ",").split(",") if part.strip()]
        try:
            return [str(part).strip() for part in value if str(part).strip()]
        except TypeError:
            return [str(value).strip()]

    required = tuple(_norm(kind) for kind in required_kinds)
    expected_phase_list = [str(phase).strip() for phase in expected_phases]
    if not required:
        raise ValueError("required_kinds must not be empty")
    if not expected_phase_list:
        raise ValueError("expected_phases must not be empty")

    expected_policies = {
        "block_labels": {"femm_block_label_source_contract_gate"},
        "pm_magnetization": {"femm_pm_magnetization_convention_gate"},
        "current_snapshot": {
            "femm_static_current_circuit_rows_gate",
            "motor_current_snapshot_table_contract_gate",
            "spwm_snapshot_current_handoff_gate",
        },
    }
    allowed_source_tools = {
        "block_labels": {"femm", "pyfemm"},
        "pm_magnetization": {"femm", "pyfemm"},
        "current_snapshot": {"femm", "pyfemm"},
    }

    details = []
    kind_counts = {}
    model_ids = []
    operating_point_ids = []
    missing_model_id = []
    missing_operating_point_id = []
    missing_paths = []
    bad_source_tool = []
    unknown_kinds = []
    bad_upstream_status = []
    bad_upstream_policy = []
    bad_current_kind = []
    phase_mismatches = []

    for index, row in enumerate(rows_in, start=1):
        if not isinstance(row, dict):
            raise ValueError("each artifact must be a dictionary")
        kind = _norm(_first(row, ("kind", "artifact_kind", "type")))
        model_id = _first(row, ("model_id", "motor_model_id", "geometry_id"))
        operating_point_id = _first(row, ("operating_point_id", "op_id", "snapshot_id"))
        source_tool = _first(row, ("source_tool", "tool", "source"))
        source_tool_norm = _norm(source_tool)
        path = _first(row, ("path", "file", "table_path", "artifact_path"))
        gate_policy = _first(row, ("gate_policy", "policy", "validator"))
        gate_policy_norm = _norm(gate_policy)
        status = _first(row, ("status", "gate_status", "validation_status"))
        status_norm = _norm(status)
        pass_flag = bool(row.get("pass", False))
        phases = _phase_list(_first(row, ("phase_set", "phases", "phase_order", "circuit_names")))

        if kind:
            kind_counts[kind] = kind_counts.get(kind, 0) + 1
        else:
            unknown_kinds.append(index)
        if kind and kind not in set(required) | set(expected_policies):
            unknown_kinds.append(kind)
        if not model_id:
            missing_model_id.append(index)
        else:
            model_ids.append(str(model_id))
        if kind == "current_snapshot":
            if not operating_point_id:
                missing_operating_point_id.append(index)
            else:
                operating_point_ids.append(str(operating_point_id))
        if source_tool_norm not in allowed_source_tools.get(kind, set()):
            bad_source_tool.append({"index": index, "kind": kind, "source_tool": source_tool})
        if not path:
            missing_paths.append(index)
        if not (pass_flag or status_norm in {"ok", "pass", "passed", "verified"}):
            bad_upstream_status.append({"index": index, "kind": kind, "status": status})
        if kind in expected_policies and gate_policy_norm not in expected_policies[kind]:
            bad_upstream_policy.append({
                "index": index,
                "kind": kind,
                "gate_policy": gate_policy,
                "expected": sorted(expected_policies[kind]),
            })
        if kind == "current_snapshot" and _norm(row.get("current_kind")) not in {
            "instantaneous",
            "instant",
            "snapshot",
            "sample",
            "peak",
        }:
            bad_current_kind.append(index)
        if phases and phases != expected_phase_list:
            phase_mismatches.append({"index": index, "kind": kind, "phases": phases})

        details.append({
            "index": index,
            "kind": kind,
            "model_id": None if model_id is None else str(model_id),
            "operating_point_id": None if operating_point_id is None else str(operating_point_id),
            "source_tool": source_tool,
            "path": path,
            "gate_policy": gate_policy,
            "status": status,
            "pass": pass_flag,
            "phases": phases,
        })

    required_set = set(required)
    present_set = set(kind_counts)
    unique_model_ids = sorted(set(model_ids))
    unique_operating_point_ids = sorted(set(operating_point_ids))
    checks = {
        "required_kinds_present": required_set.issubset(present_set),
        "no_unknown_kinds": not unknown_kinds,
        "model_ids_present": not missing_model_id,
        "model_ids_unique": len(unique_model_ids) == 1,
        "source_tools_match_kind": not bad_source_tool,
        "paths_present": not missing_paths,
        "current_operating_point_id_present": not missing_operating_point_id,
        "current_operating_point_id_unique": len(unique_operating_point_ids) == 1,
        "upstream_gate_status_ok": not bad_upstream_status,
        "upstream_gate_policy_known": not bad_upstream_policy,
        "current_snapshot_is_instantaneous": not bad_current_kind,
        "phase_sets_match_expected": not phase_mismatches,
    }
    if expected_model_id is not None:
        checks["expected_model_id_matches"] = unique_model_ids == [str(expected_model_id)]
    if expected_operating_point_id is not None:
        checks["expected_operating_point_id_matches"] = unique_operating_point_ids == [str(expected_operating_point_id)]

    return {
        "policy": "femm_source_current_solver_ready_manifest_gate",
        "required_kinds": list(required),
        "present_kinds": dict(sorted(kind_counts.items())),
        "model_ids": unique_model_ids,
        "operating_point_ids": unique_operating_point_ids,
        "expected_model_id": None if expected_model_id is None else str(expected_model_id),
        "expected_operating_point_id": None if expected_operating_point_id is None else str(expected_operating_point_id),
        "expected_phases": expected_phase_list,
        "missing_model_id_rows": missing_model_id,
        "missing_operating_point_id_rows": missing_operating_point_id,
        "bad_source_tool_rows": bad_source_tool,
        "missing_path_rows": missing_paths,
        "unknown_kinds": unknown_kinds,
        "bad_upstream_status_rows": bad_upstream_status,
        "bad_upstream_policy_rows": bad_upstream_policy,
        "bad_current_kind_rows": bad_current_kind,
        "phase_mismatch_rows": phase_mismatches,
        "artifacts": details,
        "checks": checks,
        "status": "ok" if all(checks.values()) else "needs_attention",
        "version_note": (
            "Run after FEMM block-label, PM magnetization, and static current "
            "gates, but before FEMM solve/result tables, so source metadata and "
            "instantaneous circuit currents cannot drift apart."
        ),
    }


def femm_air_gap_sample_solver_ready_manifest_gate(
    artifacts,
    expected_model_id=None,
    expected_operating_point_id=None,
    required_kinds=("source_current_manifest", "air_gap_sample_table", "torque_summary"),
):
    """Check FEMM air-gap Br/Bt samples before torque/ripple promotion.

    FEMM/pyFEMM air-gap workflows often produce source/current manifests,
    degree-sampled ``mo_getgapb`` Br/Bt rows, and radia-ngsolve Maxwell-shear
    torque summaries as separate files.  This gate keeps them locked to one
    model and one operating point before a notebook or radia validation consumes
    the package.
    """

    rows_in = list(artifacts)
    if not rows_in:
        raise ValueError("artifacts must not be empty")

    def _norm(value):
        return str(value or "").strip().lower().replace("-", "_").replace(" ", "_")

    def _first(row, names):
        for name in names:
            if name in row and row[name] is not None:
                return row[name]
        return None

    required = tuple(_norm(kind) for kind in required_kinds)
    if not required:
        raise ValueError("required_kinds must not be empty")

    expected_policies = {
        "source_current_manifest": {"femm_source_current_solver_ready_manifest_gate"},
        "air_gap_sample_table": {
            "femm_air_gap_sample_metadata_contract",
            "femm_air_gap_sample_metadata_contract_gate",
        },
        "torque_summary": {
            "air_gap_shear_torque_from_angle_samples",
            "air_gap_shear_torque_summary",
        },
    }
    allowed_source_tools = {
        "source_current_manifest": {"femm", "pyfemm"},
        "air_gap_sample_table": {"femm", "pyfemm"},
        "torque_summary": {"radia_ngsolve", "radia_mcp", "radia"},
    }

    details = []
    kind_counts = {}
    model_ids = []
    operating_point_ids = []
    missing_model_id = []
    missing_operating_point_id = []
    missing_paths = []
    unknown_kinds = []
    bad_source_tool = []
    bad_upstream_status = []
    bad_upstream_policy = []
    bad_sample_angle_unit = []
    bad_sample_frame = []
    bad_sample_radius_length = []
    bad_torque_sign = []
    torque_signs = []

    for index, row in enumerate(rows_in, start=1):
        if not isinstance(row, dict):
            raise ValueError("each artifact must be a dictionary")
        kind = _norm(_first(row, ("kind", "artifact_kind", "type")))
        model_id = _first(row, ("model_id", "motor_model_id", "geometry_id"))
        operating_point_id = _first(row, ("operating_point_id", "op_id", "snapshot_id"))
        source_tool = _first(row, ("source_tool", "tool", "source"))
        source_tool_norm = _norm(source_tool)
        path = _first(row, ("path", "file", "table_path", "artifact_path"))
        gate_policy = _first(row, ("gate_policy", "policy", "validator"))
        gate_policy_norm = _norm(gate_policy)
        status = _first(row, ("status", "gate_status", "validation_status"))
        status_norm = _norm(status)
        pass_flag = bool(row.get("pass", False))
        sign = _norm(_first(row, ("torque_sign_convention", "sign_convention")))

        if kind:
            kind_counts[kind] = kind_counts.get(kind, 0) + 1
        else:
            unknown_kinds.append(index)
        if kind and kind not in set(required) | set(expected_policies):
            unknown_kinds.append(kind)
        if not model_id:
            missing_model_id.append(index)
        else:
            model_ids.append(str(model_id))
        if not operating_point_id:
            missing_operating_point_id.append(index)
        else:
            operating_point_ids.append(str(operating_point_id))
        if source_tool_norm not in allowed_source_tools.get(kind, set()):
            bad_source_tool.append({"index": index, "kind": kind, "source_tool": source_tool})
        if not path:
            missing_paths.append(index)
        if not (pass_flag or status_norm in {"ok", "pass", "passed", "verified"}):
            bad_upstream_status.append({"index": index, "kind": kind, "status": status})
        if kind in expected_policies and gate_policy_norm not in expected_policies[kind]:
            bad_upstream_policy.append({
                "index": index,
                "kind": kind,
                "gate_policy": gate_policy,
                "expected": sorted(expected_policies[kind]),
            })
        if sign:
            torque_signs.append(sign)

        if kind == "air_gap_sample_table":
            angle_unit = _norm(_first(row, ("angle_unit", "angle_units")))
            frame = _norm(_first(row, ("component_frame", "frame")))
            radius = _first(row, ("radius_m", "air_gap_radius_m", "r_m"))
            axial_length = _first(row, ("axial_length_m", "length_m", "stack_length_m"))
            if angle_unit != "deg":
                bad_sample_angle_unit.append({"index": index, "angle_unit": angle_unit or None})
            if frame != "cylindrical_rt":
                bad_sample_frame.append({"index": index, "component_frame": frame or None})
            try:
                radius_ok = float(radius) > 0.0
                length_ok = float(axial_length) > 0.0
            except (TypeError, ValueError):
                radius_ok = False
                length_ok = False
            if not (radius_ok and length_ok):
                bad_sample_radius_length.append({
                    "index": index,
                    "radius_m": radius,
                    "axial_length_m": axial_length,
                })

        details.append({
            "index": index,
            "kind": kind,
            "model_id": None if model_id is None else str(model_id),
            "operating_point_id": None if operating_point_id is None else str(operating_point_id),
            "source_tool": source_tool,
            "path": path,
            "gate_policy": gate_policy,
            "status": status,
            "pass": pass_flag,
            "torque_sign_convention": sign or None,
        })

    unique_model_ids = sorted(set(model_ids))
    unique_operating_point_ids = sorted(set(operating_point_ids))
    unique_torque_signs = sorted(set(torque_signs))
    if len(unique_torque_signs) > 1:
        bad_torque_sign = unique_torque_signs

    checks = {
        "required_kinds_present": set(required).issubset(set(kind_counts)),
        "no_unknown_kinds": not unknown_kinds,
        "model_ids_present": not missing_model_id,
        "model_ids_unique": len(unique_model_ids) == 1,
        "operating_point_ids_present": not missing_operating_point_id,
        "operating_point_ids_unique": len(unique_operating_point_ids) == 1,
        "source_tools_match_kind": not bad_source_tool,
        "paths_present": not missing_paths,
        "upstream_gate_status_ok": not bad_upstream_status,
        "upstream_gate_policy_known": not bad_upstream_policy,
        "sample_angle_unit_is_deg": not bad_sample_angle_unit,
        "sample_component_frame_is_cylindrical_rt": not bad_sample_frame,
        "sample_radius_and_length_positive": not bad_sample_radius_length,
        "torque_sign_convention_consistent": not bad_torque_sign,
    }
    if expected_model_id is not None:
        checks["expected_model_id_matches"] = unique_model_ids == [str(expected_model_id)]
    if expected_operating_point_id is not None:
        checks["expected_operating_point_id_matches"] = unique_operating_point_ids == [str(expected_operating_point_id)]

    return {
        "policy": "femm_air_gap_sample_solver_ready_manifest_gate",
        "required_kinds": list(required),
        "present_kinds": dict(sorted(kind_counts.items())),
        "model_ids": unique_model_ids,
        "operating_point_ids": unique_operating_point_ids,
        "expected_model_id": None if expected_model_id is None else str(expected_model_id),
        "expected_operating_point_id": None if expected_operating_point_id is None else str(expected_operating_point_id),
        "torque_sign_conventions": unique_torque_signs,
        "missing_model_id_rows": missing_model_id,
        "missing_operating_point_id_rows": missing_operating_point_id,
        "missing_path_rows": missing_paths,
        "unknown_kinds": unknown_kinds,
        "bad_source_tool_rows": bad_source_tool,
        "bad_upstream_status_rows": bad_upstream_status,
        "bad_upstream_policy_rows": bad_upstream_policy,
        "bad_sample_angle_unit_rows": bad_sample_angle_unit,
        "bad_sample_frame_rows": bad_sample_frame,
        "bad_sample_radius_length_rows": bad_sample_radius_length,
        "bad_torque_sign_conventions": bad_torque_sign,
        "artifacts": details,
        "checks": checks,
        "status": "ok" if all(checks.values()) else "needs_attention",
        "version_note": (
            "Run after FEMM source/current pre-solve gates and FEMM air-gap "
            "sample metadata checks, before promoting Br/Bt rows to Maxwell "
            "shear torque, cogging, or torque-ripple notebooks."
        ),
    }


def balanced_back_emf_line_voltage_handoff_gate(
    phase_harmonic_peaks,
    measured_phase_rms=None,
    measured_line_line_rms=None,
    tol=1.0e-12,
):
    """Check balanced three-phase phase/back-EMF harmonics before line-line use.

    FEMM flux-linkage/back-EMF postprocessing often starts from phase-to-neutral
    harmonic amplitudes.  For a balanced three-phase set, the line-line peak
    factor of harmonic ``n`` is ``2*abs(sin(n*pi/3))``.  Triplen harmonics can
    be present in each phase, but cancel from line-line voltage and balanced
    instantaneous power.
    """

    if not phase_harmonic_peaks:
        raise ValueError("phase_harmonic_peaks must not be empty")
    tolerance = float(tol)
    if tolerance < 0.0:
        raise ValueError("tol must be non-negative")

    rows = []
    phase_rms_sq = 0.0
    line_rms_sq = 0.0
    triplen_line_peaks = []
    for harmonic, peak in sorted((int(k), float(v)) for k, v in phase_harmonic_peaks.items()):
        if harmonic <= 0:
            raise ValueError("harmonic numbers must be positive")
        if peak < 0.0:
            raise ValueError("phase harmonic peaks must be non-negative")
        line_factor = 2.0 * abs(math.sin(harmonic * math.pi / 3.0))
        line_peak = line_factor * peak
        phase_rms = peak / math.sqrt(2.0)
        line_rms = line_peak / math.sqrt(2.0)
        if harmonic % 3 == 0:
            triplen_line_peaks.append(abs(line_peak))
        phase_rms_sq += phase_rms * phase_rms
        line_rms_sq += line_rms * line_rms
        rows.append({
            "harmonic": harmonic,
            "phase_peak": peak,
            "line_line_peak_factor": line_factor,
            "line_line_peak": line_peak,
            "phase_rms_contribution": phase_rms,
            "line_line_rms_contribution": line_rms,
            "triplen": harmonic % 3 == 0,
        })

    phase_rms_total = math.sqrt(phase_rms_sq)
    line_rms_total = math.sqrt(line_rms_sq)
    checks = {
        "has_fundamental": any(row["harmonic"] == 1 and row["phase_peak"] > 0.0 for row in rows),
        "triplen_cancel_from_line_line": max(triplen_line_peaks or [0.0]) <= tolerance,
        "line_line_rms_not_from_triplen": line_rms_total <= math.sqrt(3.0) * phase_rms_total + tolerance,
    }
    out = {
        "policy": "balanced_back_emf_line_voltage_handoff_gate",
        "rows": rows,
        "phase_rms_total": phase_rms_total,
        "line_line_rms_total": line_rms_total,
        "fundamental_line_line_factor": math.sqrt(3.0),
        "max_triplen_line_line_peak": max(triplen_line_peaks or [0.0]),
        "tol": tolerance,
        "checks": checks,
    }
    if measured_phase_rms is not None:
        measured = float(measured_phase_rms)
        out["measured_phase_rms"] = measured
        out["measured_phase_rms_abs_error"] = abs(measured - phase_rms_total)
        checks["measured_phase_rms_ok"] = out["measured_phase_rms_abs_error"] <= tolerance
    if measured_line_line_rms is not None:
        measured = float(measured_line_line_rms)
        out["measured_line_line_rms"] = measured
        out["measured_line_line_rms_abs_error"] = abs(measured - line_rms_total)
        checks["measured_line_line_rms_ok"] = out["measured_line_line_rms_abs_error"] <= tolerance

    out["status"] = "ok" if all(checks.values()) else "needs_attention"
    return out


def flux_linkage_back_emf_derivative_gate(
    theta_rad,
    flux_linkage_wb,
    back_emf_v,
    omega_rad_per_s,
    sign_convention="faraday",
    rtol=2.0e-3,
    atol=1.0e-9,
    period_rad=2.0 * math.pi,
):
    """Check phase back-EMF rows against the derivative of flux linkage.

    For a periodic phase flux linkage table sampled over one electrical period,
    Faraday's convention is ``e = -omega * d(lambda)/d(theta)``.  The gate keeps
    angle spacing, endpoint duplication, sign convention, and speed scaling
    explicit before a motor workflow promotes flux-linkage rows to a Ke table.
    """

    theta = [float(value) for value in theta_rad]
    flux = [float(value) for value in flux_linkage_wb]
    emf = [float(value) for value in back_emf_v]
    omega = float(omega_rad_per_s)
    rel_tol = float(rtol)
    abs_tol = float(atol)
    period = float(period_rad)
    if not (len(theta) == len(flux) == len(emf)):
        raise ValueError("theta_rad, flux_linkage_wb, and back_emf_v must have the same length")
    if len(theta) < 3:
        raise ValueError("at least three samples are required")
    if omega <= 0.0:
        raise ValueError("omega_rad_per_s must be > 0")
    if rel_tol < 0.0 or abs_tol < 0.0:
        raise ValueError("rtol and atol must be non-negative")
    if period <= 0.0:
        raise ValueError("period_rad must be > 0")

    steps = [theta[i + 1] - theta[i] for i in range(len(theta) - 1)]
    if any(step <= 0.0 for step in steps):
        raise ValueError("theta samples must be strictly increasing")
    h = sum(steps) / len(steps)
    theta_spacing_abs_error = max(abs(step - h) for step in steps)
    if theta_spacing_abs_error > max(1.0e-12, 1.0e-9 * abs(h)):
        raise ValueError("theta samples must be equally spaced")
    covered_period = theta[-1] - theta[0] + h
    period_abs_error = abs(covered_period - period)
    period_rel_error = period_abs_error / max(abs(period), 1.0e-300)
    if period_abs_error > max(1.0e-12, 1.0e-9 * abs(period)):
        raise ValueError("theta samples must cover one period without a duplicate endpoint")

    convention = str(sign_convention).strip().lower().replace("-", "_")
    negative_aliases = {"faraday", "negative_derivative", "generator", "passive"}
    positive_aliases = {"positive_derivative", "motor"}
    if convention in negative_aliases:
        derivative_sign = -1.0
        normalized_convention = "negative_derivative"
    elif convention in positive_aliases:
        derivative_sign = 1.0
        normalized_convention = "positive_derivative"
    else:
        raise ValueError("sign_convention must be 'faraday'/'negative_derivative' or 'positive_derivative'")

    dflux_dtheta = central_difference_periodic(flux, h)
    expected = [derivative_sign * omega * value for value in dflux_dtheta]
    rows = []
    max_abs_error = 0.0
    max_rel_error = 0.0
    rms_error_sq = 0.0
    for index, (angle, lam, measured, reference, derivative) in enumerate(
        zip(theta, flux, emf, expected, dflux_dtheta)
    ):
        abs_error = abs(measured - reference)
        rel_error = abs_error / max(abs(measured), abs(reference), 1.0e-300)
        passed = abs_error <= abs_tol or rel_error <= rel_tol
        max_abs_error = max(max_abs_error, abs_error)
        max_rel_error = max(max_rel_error, rel_error)
        rms_error_sq += abs_error * abs_error
        rows.append({
            "sample": index,
            "theta_rad": angle,
            "flux_linkage_wb": lam,
            "dflux_dtheta_wb_per_rad": derivative,
            "expected_back_emf_v": reference,
            "back_emf_v": measured,
            "abs_error_v": abs_error,
            "rel_error": rel_error,
            "passed": passed,
        })

    checks = {
        "omega_positive": omega > 0.0,
        "theta_equally_spaced": theta_spacing_abs_error <= max(1.0e-12, 1.0e-9 * abs(h)),
        "one_period_without_duplicate_endpoint": period_abs_error <= max(1.0e-12, 1.0e-9 * abs(period)),
        "sign_convention_recorded": normalized_convention in {"negative_derivative", "positive_derivative"},
        "all_rows_within_tolerance": all(row["passed"] for row in rows),
    }
    return {
        "policy": "flux_linkage_back_emf_derivative_gate",
        "n_samples": len(rows),
        "omega_rad_per_s": omega,
        "sign_convention": normalized_convention,
        "theta_spacing_rad": h,
        "period_rad": period,
        "covered_period_rad": covered_period,
        "period_abs_error_rad": period_abs_error,
        "period_rel_error": period_rel_error,
        "max_abs_error_v": max_abs_error,
        "max_rel_error": max_rel_error,
        "rms_abs_error_v": math.sqrt(rms_error_sq / len(rows)),
        "rtol": rel_tol,
        "atol": abs_tol,
        "checks": checks,
        "rows": rows,
        "status": "ok" if all(checks.values()) else "needs_attention",
    }


def inverter_dc_bus_voltage_limit_gate(
    dc_bus_v,
    modulation_index,
    method="spwm",
    measured_line_line_rms=None,
    tol=1.0e-12,
):
    """Check fundamental line-line RMS voltage against a DC-bus limit.

    This is a small public-safe bridge from FEMM/radia phase back-EMF tables to
    drive-control notebooks.  In the linear region, sine-triangle SPWM can
    deliver ``sqrt(3)/(2*sqrt(2)) * Vdc`` line-line RMS at modulation index 1,
    while SVPWM can deliver ``Vdc/sqrt(2)``.  The two methods are separated so
    a table cannot silently use the more generous SVPWM limit for an SPWM run.
    """

    vdc = float(dc_bus_v)
    m = float(modulation_index)
    tolerance = float(tol)
    if vdc <= 0.0:
        raise ValueError("dc_bus_v must be > 0")
    if not (0.0 <= m <= 1.0):
        raise ValueError("modulation_index must be in the linear range [0, 1]")
    if tolerance < 0.0:
        raise ValueError("tol must be non-negative")

    normalized = str(method).strip().lower().replace("-", "_")
    if normalized in {"spwm", "sine_pwm", "sinusoidal_pwm"}:
        method_name = "spwm"
        factor = math.sqrt(3.0) / (2.0 * math.sqrt(2.0))
    elif normalized in {"svpwm", "space_vector_pwm", "space_vector"}:
        method_name = "svpwm"
        factor = 1.0 / math.sqrt(2.0)
    else:
        raise ValueError("method must be 'spwm' or 'svpwm'")

    line_line_rms_limit = m * vdc * factor
    phase_rms_limit = line_line_rms_limit / math.sqrt(3.0)
    line_line_peak_limit = math.sqrt(2.0) * line_line_rms_limit
    phase_peak_limit = math.sqrt(2.0) * phase_rms_limit
    payload = {
        "policy": "inverter_dc_bus_voltage_limit_gate",
        "method": method_name,
        "dc_bus_v": vdc,
        "modulation_index": m,
        "line_line_rms_factor_at_m1": factor,
        "line_line_rms_limit": line_line_rms_limit,
        "line_line_peak_limit": line_line_peak_limit,
        "phase_rms_limit": phase_rms_limit,
        "phase_peak_limit": phase_peak_limit,
        "tol": tolerance,
    }
    checks = {
        "linear_modulation_index_ok": True,
    }
    if measured_line_line_rms is not None:
        measured = float(measured_line_line_rms)
        margin = line_line_rms_limit - measured
        payload["measured_line_line_rms"] = measured
        payload["line_line_rms_margin"] = margin
        payload["line_line_rms_utilization"] = measured / line_line_rms_limit if line_line_rms_limit > 0.0 else math.inf
        checks["measured_line_line_rms_within_limit"] = margin >= -tolerance
    payload["checks"] = checks
    payload["status"] = "ok" if all(checks.values()) else "needs_attention"
    return payload


def dq_current_from_gamma_deg(current, gamma_deg):
    """Return ``(id, iq)`` for an amplitude-invariant motor current angle.

    Positive ``gamma`` means field weakening for the common IPM convention:
    ``id = -I sin(gamma)``, ``iq = I cos(gamma)``.
    """

    i = float(current)
    g = math.radians(float(gamma_deg))
    return -i * math.sin(g), i * math.cos(g)


def lumped_pm_dq_torque(lambda_m, Ld, Lq, id_current, iq_current, pole_pairs):
    """Return the closed-form PM/saliency dq torque.

    ``T = 3/2 p [lambda_m iq + (Ld - Lq) id iq]``.  This is deliberately
    duplicated here instead of importing the full NGSolve-backed solve module so
    loop table checks remain lightweight.
    """

    return 1.5 * float(pole_pairs) * (
        float(lambda_m) * float(iq_current)
        + (float(Ld) - float(Lq)) * float(id_current) * float(iq_current)
    )


def ipm_saliency_torque_component_gate(
    lambda_m,
    Ld,
    Lq,
    id_current,
    iq_current,
    pole_pairs,
    require_ipm_saliency=True,
    tol=1.0e-12,
):
    """Decompose IPM dq torque into magnet and reluctance components.

    For the common IPM convention ``Ld < Lq`` and field-weakening current
    ``id < 0``, the reluctance term
    ``1.5*p*(Ld-Lq)*id*iq`` is positive and adds to the magnet torque.  This
    is the small pre-FE gate to run before asking FEMM/JMAG/radia-ngsolve for
    an Ld/Lq map or V-type IPM torque table.
    """

    lm = float(lambda_m)
    ld = float(Ld)
    lq = float(Lq)
    id_a = float(id_current)
    iq_a = float(iq_current)
    p = float(pole_pairs)
    tolerance = float(tol)
    if ld <= 0.0 or lq <= 0.0:
        raise ValueError("Ld and Lq must be positive")
    if p <= 0.0:
        raise ValueError("pole_pairs must be positive")
    if tolerance < 0.0:
        raise ValueError("tol must be non-negative")
    magnet_torque = 1.5 * p * lm * iq_a
    reluctance_torque = 1.5 * p * (ld - lq) * id_a * iq_a
    total = magnet_torque + reluctance_torque
    direct_total = lumped_pm_dq_torque(lm, ld, lq, id_a, iq_a, p)
    saliency_ratio = lq / ld
    current = math.hypot(id_a, iq_a)
    reluctance_fraction = (
        reluctance_torque / total if abs(total) > 1.0e-300 else None
    )
    checks = {
        "torque_sum_identity": abs(total - direct_total) <= tolerance,
        "positive_current": current > 0.0,
        "ipm_saliency_ratio_gt_one": saliency_ratio > 1.0 if require_ipm_saliency else True,
        "field_weakening_negative_id": id_a < 0.0 if require_ipm_saliency else True,
        "positive_iq_for_motoring": iq_a > 0.0,
        "reluctance_torque_adds_to_magnet": reluctance_torque >= -tolerance if require_ipm_saliency else True,
    }
    return {
        "policy": "ipm_saliency_torque_component_gate",
        "lambda_m_Wb": lm,
        "Ld_H": ld,
        "Lq_H": lq,
        "saliency_ratio_Lq_over_Ld": saliency_ratio,
        "id_A": id_a,
        "iq_A": iq_a,
        "current_A": current,
        "pole_pairs": p,
        "magnet_torque_Nm": magnet_torque,
        "reluctance_torque_Nm": reluctance_torque,
        "total_torque_Nm": total,
        "direct_total_torque_Nm": direct_total,
        "reluctance_fraction_of_total": reluctance_fraction,
        "checks": checks,
        "status": "ok" if all(checks.values()) else "needs_attention",
        "tol": tolerance,
    }


def jmag_pm_short_circuit_fault_table_gate(
    rows,
    R,
    Ld,
    Lq,
    lambda_m,
    pole_pairs,
    *,
    residual_tol=1.0e-9,
    value_tol=1.0e-9,
    require_high_speed_demag_fraction=0.95,
):
    """Check a PM-machine short-circuit fault table before JMAG comparison.

    The table should describe a shorted-terminal dq sweep, usually exported
    from a JMAG fault/protection study.  This public-safe gate verifies the
    closed-form short-circuit currents, residual ``vd=vq=0`` equations,
    characteristic-current ratio, d-axis demagnetizing fraction, and braking
    torque trend before private FE values are compared.
    """

    table = list(rows)
    if len(table) < 3:
        raise ValueError("rows must contain at least three speed points")
    r = float(R)
    ld = float(Ld)
    lq = float(Lq)
    lm = float(lambda_m)
    pp = float(pole_pairs)
    if r < 0.0:
        raise ValueError("R must be non-negative")
    if ld <= 0.0 or lq <= 0.0:
        raise ValueError("Ld and Lq must be positive")
    if pp <= 0.0:
        raise ValueError("pole_pairs must be positive")
    residual_tol = float(residual_tol)
    value_tol = float(value_tol)
    ich = lm / ld

    def _first(row, *names):
        for name in names:
            if name in row and row[name] is not None:
                return row[name]
        return None

    summaries = []
    omega_values = []
    residual_errors = []
    id_errors = []
    iq_errors = []
    torque_errors = []
    ratio_errors = []
    demag_fraction_errors = []
    mech_speed_errors = []
    missing_required = []
    negative_id_ok = []
    negative_iq_ok = []

    for index, row in enumerate(table, start=1):
        omega_raw = _first(row, "omega_e", "omega_e_rad_per_s", "omega_e_rad_s")
        id_raw = _first(row, "id_A", "id")
        iq_raw = _first(row, "iq_A", "iq")
        torque_raw = _first(row, "torque_Nm", "torque")
        ratio_raw = _first(row, "current_ratio_to_characteristic", "current_ratio")
        demag_raw = _first(row, "d_axis_demag_fraction", "demag_fraction")
        vd_raw = _first(row, "vd_residual", "vd_residual_V", "vd_V")
        vq_raw = _first(row, "vq_residual", "vq_residual_V", "vq_V")
        if any(value is None for value in (omega_raw, id_raw, iq_raw, torque_raw, ratio_raw, demag_raw, vd_raw, vq_raw)):
            missing_required.append(index)
            continue
        omega = float(omega_raw)
        id_got = float(id_raw)
        iq_got = float(iq_raw)
        torque_got = float(torque_raw)
        ratio_got = float(ratio_raw)
        demag_got = float(demag_raw)
        vd = float(vd_raw)
        vq = float(vq_raw)
        den = r * r + omega * omega * ld * lq
        id_ref = -omega * omega * lq * lm / den
        iq_ref = -omega * r * lm / den
        imag_ref = math.hypot(id_ref, iq_ref)
        ratio_ref = math.inf if ich == 0.0 else imag_ref / abs(ich)
        demag_ref = math.inf if ich == 0.0 else -id_ref / ich
        torque_ref = lumped_pm_dq_torque(lm, ld, lq, id_ref, iq_ref, pp)
        vd_ref = r * id_got - omega * lq * iq_got
        vq_ref = r * iq_got + omega * (ld * id_got + lm)
        residual = max(abs(vd), abs(vq), abs(vd_ref), abs(vq_ref))
        omega_mech = _first(row, "omega_mech", "omega_mech_rad_per_s", "omega_mech_rad_s")
        mech_error = 0.0
        if omega_mech is not None:
            mech_error = abs(float(omega_mech) - omega / pp)
            mech_speed_errors.append(mech_error)

        omega_values.append(omega)
        residual_errors.append(residual)
        id_errors.append(abs(id_got - id_ref))
        iq_errors.append(abs(iq_got - iq_ref))
        torque_errors.append(abs(torque_got - torque_ref))
        ratio_errors.append(abs(ratio_got - ratio_ref))
        demag_fraction_errors.append(abs(demag_got - demag_ref))
        negative_id_ok.append(id_got <= value_tol)
        negative_iq_ok.append(iq_got <= value_tol)
        summaries.append({
            "row": index,
            "omega_e": omega,
            "omega_mech": None if omega_mech is None else float(omega_mech),
            "id_A": id_got,
            "iq_A": iq_got,
            "expected_id_A": id_ref,
            "expected_iq_A": iq_ref,
            "torque_Nm": torque_got,
            "expected_torque_Nm": torque_ref,
            "current_ratio_to_characteristic": ratio_got,
            "expected_current_ratio_to_characteristic": ratio_ref,
            "d_axis_demag_fraction": demag_got,
            "expected_d_axis_demag_fraction": demag_ref,
            "max_short_terminal_residual_V": residual,
            "omega_mech_abs_error": mech_error,
        })

    torque_magnitudes = [abs(row["torque_Nm"]) for row in summaries]
    peak_index = max(range(len(torque_magnitudes)), key=lambda index: torque_magnitudes[index]) if torque_magnitudes else -1
    high_speed_row = max(summaries, key=lambda row: row["omega_e"]) if summaries else {}
    omega_strict = all(a < b for a, b in zip(omega_values, omega_values[1:]))
    checks = {
        "required_columns_present": not missing_required,
        "omega_e_strictly_increasing": omega_strict,
        "omega_mech_matches_pole_pairs": max(mech_speed_errors, default=0.0) <= value_tol,
        "short_terminal_residuals_ok": max(residual_errors, default=math.inf) <= residual_tol,
        "id_column_matches_closed_form": max(id_errors, default=math.inf) <= value_tol,
        "iq_column_matches_closed_form": max(iq_errors, default=math.inf) <= value_tol,
        "torque_column_matches_closed_form": max(torque_errors, default=math.inf) <= value_tol,
        "current_ratio_matches_characteristic_current": max(ratio_errors, default=math.inf) <= value_tol,
        "d_axis_demag_fraction_matches_characteristic_current": max(demag_fraction_errors, default=math.inf) <= value_tol,
        "short_circuit_currents_are_demagnetizing": all(negative_id_ok),
        "q_axis_current_is_braking_direction": all(negative_iq_ok),
        "braking_torque_peaks_at_intermediate_speed": 0 < peak_index < len(summaries) - 1,
        "high_speed_demag_fraction_near_one": (
            float(high_speed_row.get("d_axis_demag_fraction", 0.0)) >= float(require_high_speed_demag_fraction)
        ),
    }
    return {
        "policy": "jmag_pm_short_circuit_fault_table_gate",
        "row_count": len(table),
        "R_ohm": r,
        "Ld_H": ld,
        "Lq_H": lq,
        "lambda_m_Wb": lm,
        "pole_pairs": pp,
        "characteristic_current_A": ich,
        "residual_tol": residual_tol,
        "value_tol": value_tol,
        "require_high_speed_demag_fraction": float(require_high_speed_demag_fraction),
        "missing_required_rows": missing_required,
        "max_short_terminal_residual_V": max(residual_errors, default=math.inf),
        "max_id_abs_error_A": max(id_errors, default=math.inf),
        "max_iq_abs_error_A": max(iq_errors, default=math.inf),
        "max_torque_abs_error_Nm": max(torque_errors, default=math.inf),
        "max_current_ratio_abs_error": max(ratio_errors, default=math.inf),
        "max_demag_fraction_abs_error": max(demag_fraction_errors, default=math.inf),
        "peak_braking_row": None if peak_index < 0 else summaries[peak_index],
        "high_speed_row": high_speed_row,
        "checks": checks,
        "status": "ok" if all(checks.values()) else "needs_attention",
        "rows": summaries,
    }


def coaxial_pm_force_gap_sweep_gate(
    rows,
    gap_key="gap_m",
    force_key="force_N",
    expected_sign="attractive_negative",
    rtol_invariant=1.0e-12,
):
    """Check a PM force-gap sweep against the dipole-limit fourth-power law.

    This is a public-safe preflight for ELF/MAGIC, radia-ngsolve, or other
    magnetostatic force tables before comparing detailed BEM/FEM runs.  In the
    far-field dipole limit the axial force magnitude scales as ``1 / gap**4``.
    """

    table = sorted(list(rows), key=lambda row: float(row[gap_key]))
    if len(table) < 2:
        raise ValueError("rows must contain at least two force-gap samples")
    gaps = [float(row[gap_key]) for row in table]
    forces = [float(row[force_key]) for row in table]
    if any(gap <= 0.0 for gap in gaps):
        raise ValueError("all gaps must be positive")
    if any(force == 0.0 for force in forces):
        raise ValueError("force samples must be nonzero")

    magnitudes = [abs(force) for force in forces]
    invariants = [mag * gap**4 for mag, gap in zip(magnitudes, gaps)]
    invariant_mean = sum(invariants) / len(invariants)
    max_invariant_rel_error = max(
        abs(value - invariant_mean) / max(abs(invariant_mean), 1.0e-300)
        for value in invariants
    )
    first_last_ratio = magnitudes[0] / magnitudes[-1]
    expected_ratio = (gaps[-1] / gaps[0]) ** 4
    ratio_rel_error = abs(first_last_ratio - expected_ratio) / max(abs(expected_ratio), 1.0e-300)
    monotone_decreasing = all(a > b for a, b in zip(magnitudes, magnitudes[1:]))

    if expected_sign == "attractive_negative":
        sign_ok = all(force < 0.0 for force in forces)
    elif expected_sign == "repulsive_positive":
        sign_ok = all(force > 0.0 for force in forces)
    elif expected_sign in (None, "any"):
        sign_ok = all(force > 0.0 for force in forces) or all(force < 0.0 for force in forces)
    else:
        raise ValueError("expected_sign must be attractive_negative, repulsive_positive, any, or None")

    checks = {
        "gap_axis_positive": True,
        "force_sign_ok": sign_ok,
        "force_magnitude_monotone_decreasing": monotone_decreasing,
        "fourth_power_invariant_ok": max_invariant_rel_error <= float(rtol_invariant),
        "first_last_ratio_ok": ratio_rel_error <= float(rtol_invariant),
    }
    return {
        "policy": "coaxial_pm_force_gap_sweep_gate",
        "row_count": len(table),
        "gap_key": gap_key,
        "force_key": force_key,
        "expected_sign": expected_sign,
        "gaps_m": gaps,
        "forces_N": forces,
        "force_magnitudes_N": magnitudes,
        "force_gap4_invariants": invariants,
        "mean_force_gap4_invariant": invariant_mean,
        "max_force_gap4_invariant_rel_error": max_invariant_rel_error,
        "force_ratio_first_last": first_last_ratio,
        "expected_force_ratio_first_last": expected_ratio,
        "force_ratio_rel_error": ratio_rel_error,
        "checks": checks,
        "status": "ok" if all(checks.values()) else "needs_attention",
        "rtol_invariant": float(rtol_invariant),
    }


def maxwell_stress_surface_package_gate(
    artifacts,
    *,
    required_kinds=("stress_surface", "solve_command", "observable_table"),
    expected_case_id=None,
    expected_surface_id=None,
    expected_result_set_id=None,
    expected_observable_id=None,
    required_axis=None,
    expected_normal_orientation=None,
    expected_formulation_id=None,
    expected_kernel_family=None,
    expected_run_artifact_id=None,
    expected_result_revision_id=None,
):
    """Check a Maxwell-stress force/torque package before reading values.

    This is a public-safe pre-value gate for ELF/MAGIC, NGSolve, or other
    magnetostatic workflows.  It verifies that the stress surface, solve
    command, and observable rows are one package with a closed surface and an
    explicit Maxwell-stress method before numerical force/torque comparisons.
    """

    rows = list(artifacts)
    if not rows:
        raise ValueError("artifacts must not be empty")

    def _norm(value):
        return str(value or "").strip().lower().replace("-", "_").replace(" ", "_")

    def _first(row, names):
        for name in names:
            if name in row and row[name] is not None:
                return row[name]
        return None

    def _status_ok(row):
        status = _norm(row.get("status", "ok"))
        return status in {"ok", "pass", "passed", "verified"}

    def _force_unit_norm(value):
        unit = _norm(value)
        unit = unit.replace("newtons", "n").replace("newton", "n")
        unit = unit.replace("per_meter", "per_m")
        return unit

    required = tuple(_norm(kind) for kind in required_kinds)
    kinds = [_norm(row.get("kind", row.get("artifact_kind"))) for row in rows]
    kind_set = set(kinds)
    missing_kinds = [kind for kind in required if kind not in kind_set]
    case_ids = [str(value) for value in (_first(row, ("case_id", "package_id")) for row in rows) if value]
    surface_ids = [
        str(value)
        for value in (
            _first(row, ("stress_surface_id", "surface_id", "selection_id", "mcm_id", "ecm_id"))
            for row in rows
        )
        if value
    ]
    result_set_ids = [
        str(value)
        for value in (
            _first(row, ("result_set_id", "solve_result_id", "run_result_id")) for row in rows
        )
        if value
    ]
    run_artifact_ids = [
        str(value)
        for value in (
            _first(row, ("run_artifact_id", "run_manifest_id", "solver_run_id", "execution_id"))
            for row in rows
        )
        if value
    ]
    result_revision_ids = [
        str(value)
        for value in (
            _first(row, ("result_revision_id", "result_version_id", "postprocess_revision_id"))
            for row in rows
        )
        if value
    ]
    observable_ids = [
        str(value)
        for value in (
            _first(row, ("observable_id", "observable_table_id", "result_table_id", "force_observable_id"))
            for row in rows
        )
        if value
    ]
    axes = [
        _norm(value)
        for value in (
            _first(row, ("axis", "force_axis", "torque_axis", "component_axis", "observable_axis"))
            for row in rows
        )
        if value
    ]
    methods = [
        _norm(value)
        for value in (
            _first(row, ("method", "force_method", "solve_method", "postprocess_method", "sol_command"))
            for row in rows
        )
        if value
    ]
    symmetry_values = [
        float(value)
        for value in (_first(row, ("symmetry_factor", "scale_factor")) for row in rows)
        if value is not None
    ]
    closed_flags = [
        bool(value)
        for row, kind in zip(rows, kinds)
        for value in (row.get("closed_surface", row.get("surface_closed")),)
        if kind == "stress_surface" and value is not None
    ]
    normal_orientations = [
        _norm(value)
        for value in (
            _first(row, ("normal_orientation", "surface_normal_orientation")) for row in rows
        )
        if value
    ]
    sign_conventions = [
        _norm(value)
        for value in (_first(row, ("sign_convention", "torque_sign_convention")) for row in rows)
        if value
    ]
    quantity_dimensions = [
        _norm(value)
        for value in (_first(row, ("quantity_dimension", "force_quantity_dimension")) for row in rows)
        if value
    ]
    force_units = [
        _force_unit_norm(value)
        for value in (_first(row, ("force_unit", "force_units", "unit")) for row in rows)
        if value
    ]
    formulation_ids = [
        _norm(value)
        for value in (
            _first(row, ("formulation_id", "bem_formulation", "solver_formulation"))
            for row in rows
        )
        if value
    ]
    kernel_families = [
        _norm(value)
        for value in (
            _first(row, ("kernel_family", "bem_kernel_family", "integral_kernel", "kernel"))
            for row in rows
        )
        if value
    ]
    singular_treatments = [
        _norm(value)
        for value in (
            _first(row, ("singular_treatment", "singular_quadrature", "regularization_id"))
            for row in rows
        )
        if value
    ]
    allowed_force_units_by_dimension = {
        "2d_per_length": {"n/m", "n_per_m"},
        "3d_total": {"n"},
    }
    force_quantity_metadata_rows = []
    bad_force_quantity_dimensions = []
    missing_force_units = []
    missing_quantity_dimensions = []
    bad_force_unit_dimension_pairs = []
    for index, row in enumerate(rows, start=1):
        quantity_dimension = _norm(_first(row, ("quantity_dimension", "force_quantity_dimension")))
        force_unit = _force_unit_norm(_first(row, ("force_unit", "force_units", "unit")))
        if not quantity_dimension and not force_unit:
            continue
        force_quantity_metadata_rows.append(index)
        if quantity_dimension not in allowed_force_units_by_dimension:
            if quantity_dimension:
                bad_force_quantity_dimensions.append({"index": index, "quantity_dimension": quantity_dimension})
            else:
                missing_quantity_dimensions.append(index)
            continue
        if not force_unit:
            missing_force_units.append(index)
            continue
        if force_unit not in allowed_force_units_by_dimension[quantity_dimension]:
            bad_force_unit_dimension_pairs.append({
                "index": index,
                "quantity_dimension": quantity_dimension,
                "force_unit": force_unit,
                "expected_force_units": sorted(allowed_force_units_by_dimension[quantity_dimension]),
            })
    maxwell_methods = {"fort", "sol_fort", "maxwell_stress", "maxwell_stress_fort"}
    method_ok = bool(methods) and all(
        method in maxwell_methods or method.endswith("_fort") or "maxwell" in method
        for method in methods
    )
    case_consistent = bool(case_ids) and len(set(case_ids)) == 1
    surface_consistent = bool(surface_ids) and len(set(surface_ids)) == 1
    result_set_consistent = not result_set_ids or len(set(result_set_ids)) == 1
    run_artifact_consistent = not run_artifact_ids or len(set(run_artifact_ids)) == 1
    result_revision_consistent = not result_revision_ids or len(set(result_revision_ids)) == 1
    observable_consistent = not observable_ids or len(set(observable_ids)) == 1
    axis_consistent = bool(axes) and len(set(axes)) == 1
    expected_case_ok = expected_case_id is None or (case_consistent and case_ids[0] == str(expected_case_id))
    expected_surface_ok = expected_surface_id is None or (
        surface_consistent and surface_ids[0] == str(expected_surface_id)
    )
    expected_result_set_ok = expected_result_set_id is None or (
        result_set_consistent and bool(result_set_ids) and result_set_ids[0] == str(expected_result_set_id)
    )
    expected_run_artifact_ok = expected_run_artifact_id is None or (
        run_artifact_consistent
        and bool(run_artifact_ids)
        and run_artifact_ids[0] == str(expected_run_artifact_id)
    )
    expected_result_revision_ok = expected_result_revision_id is None or (
        result_revision_consistent
        and bool(result_revision_ids)
        and result_revision_ids[0] == str(expected_result_revision_id)
    )
    expected_observable_ok = expected_observable_id is None or (
        observable_consistent and bool(observable_ids) and observable_ids[0] == str(expected_observable_id)
    )
    required_axis_key = None if required_axis is None else _norm(required_axis)
    required_axis_ok = required_axis_key is None or (axis_consistent and axes[0] == required_axis_key)
    expected_normal_key = None if expected_normal_orientation is None else _norm(expected_normal_orientation)
    normal_orientation_consistent = not normal_orientations or len(set(normal_orientations)) == 1
    expected_normal_ok = expected_normal_key is None or (
        normal_orientation_consistent
        and bool(normal_orientations)
        and normal_orientations[0] == expected_normal_key
    )
    expected_formulation_key = None if expected_formulation_id is None else _norm(expected_formulation_id)
    expected_kernel_key = None if expected_kernel_family is None else _norm(expected_kernel_family)
    formulation_consistent = not formulation_ids or len(set(formulation_ids)) == 1
    kernel_family_consistent = not kernel_families or len(set(kernel_families)) == 1
    singular_treatment_consistent = not singular_treatments or len(set(singular_treatments)) == 1
    expected_formulation_ok = expected_formulation_key is None or (
        formulation_consistent and bool(formulation_ids) and formulation_ids[0] == expected_formulation_key
    )
    expected_kernel_ok = expected_kernel_key is None or (
        kernel_family_consistent and bool(kernel_families) and kernel_families[0] == expected_kernel_key
    )
    observable_axis_mismatches = []
    seen_observable_axis_mismatches = set()
    if axis_consistent and axes and axes[0] in {"x", "y", "z"}:
        axis_key = axes[0]
        for observable_id in observable_ids:
            tokens = set(_norm(observable_id).split("_"))
            named_axes = tokens.intersection({"x", "y", "z"})
            if named_axes and axis_key not in named_axes:
                mismatch_key = (observable_id, axis_key, tuple(sorted(named_axes)))
                if mismatch_key in seen_observable_axis_mismatches:
                    continue
                seen_observable_axis_mismatches.add(mismatch_key)
                observable_axis_mismatches.append({
                    "observable_id": observable_id,
                    "axis": axis_key,
                    "observable_axis_tokens": sorted(named_axes),
                })
    checks = {
        "required_kinds_present": not missing_kinds,
        "case_id_consistent": case_consistent,
        "case_id_matches_expected": expected_case_ok,
        "stress_surface_id_consistent": surface_consistent,
        "stress_surface_id_matches_expected": expected_surface_ok,
        "result_set_id_consistent_when_present": result_set_consistent,
        "result_set_id_matches_expected": expected_result_set_ok,
        "run_artifact_id_consistent_when_present": run_artifact_consistent,
        "run_artifact_id_matches_expected": expected_run_artifact_ok,
        "result_revision_id_consistent_when_present": result_revision_consistent,
        "result_revision_id_matches_expected": expected_result_revision_ok,
        "observable_id_consistent_when_present": observable_consistent,
        "observable_id_matches_expected": expected_observable_ok,
        "observable_id_axis_matches_axis_when_named": not observable_axis_mismatches,
        "closed_surface_confirmed": bool(closed_flags) and all(closed_flags),
        "maxwell_stress_method_declared": method_ok,
        "axis_consistent": axis_consistent,
        "axis_matches_required": required_axis_ok,
        "normal_orientation_recorded": bool(normal_orientations),
        "normal_orientation_consistent_when_present": normal_orientation_consistent,
        "normal_orientation_matches_expected": expected_normal_ok,
        "formulation_id_consistent_when_present": formulation_consistent,
        "formulation_id_matches_expected": expected_formulation_ok,
        "kernel_family_consistent_when_present": kernel_family_consistent,
        "kernel_family_matches_expected": expected_kernel_ok,
        "singular_treatment_consistent_when_present": singular_treatment_consistent,
        "sign_convention_recorded": bool(sign_conventions),
        "symmetry_factor_positive_when_present": all(value > 0.0 for value in symmetry_values),
        "force_quantity_dimension_allowed_when_present": not bad_force_quantity_dimensions,
        "force_unit_recorded_when_quantity_dimension_present": not missing_force_units,
        "force_quantity_dimension_recorded_when_unit_present": not missing_quantity_dimensions,
        "force_unit_matches_quantity_dimension_when_present": not bad_force_unit_dimension_pairs,
        "all_artifacts_status_ok": all(_status_ok(row) for row in rows),
    }
    return {
        "policy": "maxwell_stress_surface_package_gate",
        "row_count": len(rows),
        "required_kinds": list(required),
        "present_kinds": sorted(kind_set),
        "missing_kinds": missing_kinds,
        "case_ids": sorted(set(case_ids)),
        "stress_surface_ids": sorted(set(surface_ids)),
        "result_set_ids": sorted(set(result_set_ids)),
        "run_artifact_ids": sorted(set(run_artifact_ids)),
        "result_revision_ids": sorted(set(result_revision_ids)),
        "observable_ids": sorted(set(observable_ids)),
        "axes": sorted(set(axes)),
        "methods": sorted(set(methods)),
        "symmetry_factors": symmetry_values,
        "normal_orientations": sorted(set(normal_orientations)),
        "expected_normal_orientation": expected_normal_key,
        "formulation_ids": sorted(set(formulation_ids)),
        "kernel_families": sorted(set(kernel_families)),
        "singular_treatments": sorted(set(singular_treatments)),
        "expected_formulation_id": expected_formulation_key,
        "expected_kernel_family": expected_kernel_key,
        "sign_conventions": sorted(set(sign_conventions)),
        "quantity_dimensions": sorted(set(quantity_dimensions)),
        "force_units": sorted(set(force_units)),
        "force_quantity_metadata_rows": force_quantity_metadata_rows,
        "bad_force_quantity_dimensions": bad_force_quantity_dimensions,
        "missing_force_unit_rows": missing_force_units,
        "missing_quantity_dimension_rows": missing_quantity_dimensions,
        "bad_force_unit_dimension_pairs": bad_force_unit_dimension_pairs,
        "observable_axis_mismatches": observable_axis_mismatches,
        "checks": checks,
        "status": "ok" if all(checks.values()) else "needs_attention",
        "version_note": (
            "Use this before Maxwell-stress force/torque values: a closed "
            "stress surface, method, axis, sign convention, symmetry factor, "
            "package identity, run/result revision identity, formulation/kernel "
            "metadata, and any supplied force quantity dimension/unit pair must be explicit."
        ),
    }


def dq_torque_table_health(
    rows,
    lambda_m,
    Ld,
    Lq,
    current,
    pole_pairs,
    tol=1.0e-9,
    gamma_tol_deg=1.0e-9,
):
    """Check a dq torque table against the closed-form current-angle contract.

    The expected public-safe schema is one row per current angle with
    ``gamma_deg``, ``id_A``, ``iq_A``, and ``torque_Nm``.  This catches the
    common export mistakes before a private FE map is compared: degree/radian
    confusion, peak/RMS current mismatch, phase/sign convention drift, and
    choosing the wrong optimum row.
    """

    required = ("gamma_deg", "id_A", "iq_A", "torque_Nm")
    table = list(rows)
    if not table:
        raise ValueError("rows must not be empty")
    for index, row in enumerate(table):
        missing = [name for name in required if name not in row]
        if missing:
            raise ValueError(f"row {index} is missing required columns: {missing}")

    i_ref = float(current)
    tol = float(tol)
    row_summaries = []
    max_current_abs_error = 0.0
    max_id_abs_error = 0.0
    max_iq_abs_error = 0.0
    max_torque_abs_error = 0.0
    max_torque_rel_error = 0.0
    for row in table:
        gamma = float(row["gamma_deg"])
        id_ref, iq_ref = dq_current_from_gamma_deg(i_ref, gamma)
        torque_ref = lumped_pm_dq_torque(lambda_m, Ld, Lq, id_ref, iq_ref, pole_pairs)
        id_got = float(row["id_A"])
        iq_got = float(row["iq_A"])
        torque_got = float(row["torque_Nm"])
        current_abs_error = abs(math.hypot(id_got, iq_got) - i_ref)
        id_abs_error = abs(id_got - id_ref)
        iq_abs_error = abs(iq_got - iq_ref)
        torque_abs_error = abs(torque_got - torque_ref)
        torque_rel_error = torque_abs_error / max(abs(torque_got), abs(torque_ref), 1.0e-300)
        max_current_abs_error = max(max_current_abs_error, current_abs_error)
        max_id_abs_error = max(max_id_abs_error, id_abs_error)
        max_iq_abs_error = max(max_iq_abs_error, iq_abs_error)
        max_torque_abs_error = max(max_torque_abs_error, torque_abs_error)
        max_torque_rel_error = max(max_torque_rel_error, torque_rel_error)
        row_summaries.append({
            "gamma_deg": gamma,
            "id_A": id_got,
            "iq_A": iq_got,
            "torque_Nm": torque_got,
            "expected_id_A": id_ref,
            "expected_iq_A": iq_ref,
            "expected_torque_Nm": torque_ref,
            "current_abs_error_A": current_abs_error,
            "id_abs_error_A": id_abs_error,
            "iq_abs_error_A": iq_abs_error,
            "torque_abs_error_Nm": torque_abs_error,
            "torque_rel_error": torque_rel_error,
        })

    pure_q = min(row_summaries, key=lambda row: abs(row["gamma_deg"]))
    peak = max(row_summaries, key=lambda row: row["torque_Nm"])
    expected_peak = max(row_summaries, key=lambda row: row["expected_torque_Nm"])
    positive_gamma_rows = [row for row in row_summaries if row["gamma_deg"] > gamma_tol_deg]
    if float(Lq) > float(Ld) and positive_gamma_rows:
        positive_gamma_sign_ok = all(row["id_A"] <= tol for row in positive_gamma_rows)
    else:
        positive_gamma_sign_ok = True

    checks = {
        "required_columns_ok": True,
        "current_magnitude_ok": max_current_abs_error <= tol,
        "id_column_ok": max_id_abs_error <= tol,
        "iq_column_ok": max_iq_abs_error <= tol,
        "torque_column_ok": max_torque_abs_error <= tol,
        "pure_q_row_present": abs(pure_q["gamma_deg"]) <= gamma_tol_deg,
        "peak_row_matches_closed_form": abs(peak["gamma_deg"] - expected_peak["gamma_deg"]) <= gamma_tol_deg,
        "positive_gamma_negative_id_ok": positive_gamma_sign_ok,
    }
    return {
        "policy": "dq_torque_table_schema_health_gate",
        "row_count": len(row_summaries),
        "required_columns": list(required),
        "lambda_m_Wb": float(lambda_m),
        "Ld_H": float(Ld),
        "Lq_H": float(Lq),
        "current_A": i_ref,
        "pole_pairs": float(pole_pairs),
        "tol": tol,
        "gamma_tol_deg": float(gamma_tol_deg),
        "max_current_abs_error_A": max_current_abs_error,
        "max_id_abs_error_A": max_id_abs_error,
        "max_iq_abs_error_A": max_iq_abs_error,
        "max_torque_abs_error_Nm": max_torque_abs_error,
        "max_torque_rel_error": max_torque_rel_error,
        "pure_q_row": pure_q,
        "peak_row": peak,
        "expected_peak_row": expected_peak,
        "checks": checks,
        "status": "ok" if all(checks.values()) else "needs_attention",
        "rows": row_summaries,
    }


def pm_drive_terminal_table_health(rows, pole_pairs, tol=1.0e-9):
    """Check PM drive terminal-quantity rows before reading an FE map.

    The expected public-safe schema is one row per operating point with
    electrical/mechanical speed, dq current, dq terminal voltage, power split,
    power factor, and efficiency.  This catches the table mistakes that often
    happen between motor FEA and drive/control notebooks: electrical vs
    mechanical speed, voltage magnitude, current magnitude, ``P_in=P_em+P_cu``,
    and ``P_em=T*omega_mech``.

    If ``voltage_utilization_lossless`` is present, it is treated as the
    selector/feasibility voltage.  The R-included terminal
    ``voltage_utilization`` is a report value and may be slightly above one
    without making the row fail.
    """

    required = (
        "omega_e",
        "omega_mech",
        "id_A",
        "iq_A",
        "vd_V",
        "vq_V",
        "Vmag_V",
        "Imag_A",
        "torque_Nm",
        "P_in_W",
        "P_em_W",
        "P_cu_W",
        "power_factor",
        "efficiency",
    )
    table = list(rows)
    if not table:
        raise ValueError("rows must not be empty")
    for index, row in enumerate(table):
        missing = [name for name in required if name not in row]
        if missing:
            raise ValueError(f"row {index} is missing required columns: {missing}")

    pp = float(pole_pairs)
    if pp <= 0.0:
        raise ValueError("pole_pairs must be > 0")
    tolerance = float(tol)

    row_summaries = []
    speed_errors = []
    voltage_errors = []
    current_errors = []
    power_balance_errors = []
    torque_speed_errors = []
    power_factor_errors = []
    efficiency_errors = []
    voltage_utilization_errors = []
    lossless_voltage_utils = []
    terminal_voltage_utils = []
    speed_order = []
    pf_values = []

    for row in table:
        omega_e = float(row["omega_e"])
        omega_mech = float(row["omega_mech"])
        id_current = float(row["id_A"])
        iq_current = float(row["iq_A"])
        vd = float(row["vd_V"])
        vq = float(row["vq_V"])
        vmag = float(row["Vmag_V"])
        imag = float(row["Imag_A"])
        torque = float(row["torque_Nm"])
        p_in = float(row["P_in_W"])
        p_em = float(row["P_em_W"])
        p_cu = float(row["P_cu_W"])
        power_factor = float(row["power_factor"])
        efficiency = float(row["efficiency"])

        vmag_ref = math.hypot(vd, vq)
        imag_ref = math.hypot(id_current, iq_current)
        speed_error = abs(omega_e - pp * omega_mech) / max(abs(omega_e), 1.0)
        voltage_error = abs(vmag - vmag_ref)
        current_error = abs(imag - imag_ref)
        power_scale = max(abs(p_in), abs(p_em), abs(p_cu), 1.0)
        power_balance_error = abs(p_in - p_em - p_cu) / power_scale
        torque_speed_error = abs(p_em - torque * omega_mech) / max(abs(p_em), abs(torque * omega_mech), 1.0)
        denom = 1.5 * vmag * imag
        pf_ref = p_in / denom if denom > 0.0 else 0.0
        pf_error = abs(power_factor - pf_ref)
        eta_ref = p_em / p_in if p_in > 0.0 else math.nan
        eta_error = abs(efficiency - eta_ref) if math.isfinite(eta_ref) else math.nan

        speed_errors.append(speed_error)
        voltage_errors.append(voltage_error)
        current_errors.append(current_error)
        power_balance_errors.append(power_balance_error)
        torque_speed_errors.append(torque_speed_error)
        power_factor_errors.append(pf_error)
        if math.isfinite(eta_error):
            efficiency_errors.append(eta_error)
        speed_order.append(omega_e)
        pf_values.append(power_factor)

        vmax = row.get("Vmax_V")
        util_report = row.get("voltage_utilization")
        util_lossless = row.get("voltage_utilization_lossless")
        if vmax is not None and util_report is not None:
            util_ref = vmag / float(vmax)
            voltage_utilization_errors.append(abs(float(util_report) - util_ref))
        if util_report is not None:
            terminal_voltage_utils.append(float(util_report))
        if util_lossless is not None:
            lossless_voltage_utils.append(float(util_lossless))

        row_summaries.append({
            "omega_e": omega_e,
            "omega_mech": omega_mech,
            "region": str(row.get("region", "")),
            "speed_multiple": row.get("speed_multiple"),
            "id_A": id_current,
            "iq_A": iq_current,
            "Vmag_V": vmag,
            "Imag_A": imag,
            "torque_Nm": torque,
            "P_in_W": p_in,
            "P_em_W": p_em,
            "P_cu_W": p_cu,
            "power_factor": power_factor,
            "efficiency": efficiency,
            "voltage_utilization": util_report,
            "voltage_utilization_lossless": util_lossless,
            "speed_contract_rel_error": speed_error,
            "voltage_magnitude_abs_error_V": voltage_error,
            "current_magnitude_abs_error_A": current_error,
            "power_balance_rel_error": power_balance_error,
            "torque_speed_rel_error": torque_speed_error,
            "power_factor_abs_error": pf_error,
            "efficiency_abs_error": eta_error,
        })

    max_lossless_util = max(lossless_voltage_utils) if lossless_voltage_utils else math.nan
    max_terminal_util = max(terminal_voltage_utils) if terminal_voltage_utils else math.nan
    terminal_over_limit_rows = [
        row for row in row_summaries
        if row["voltage_utilization"] is not None and float(row["voltage_utilization"]) > 1.0 + tolerance
    ]
    checks = {
        "required_columns_ok": True,
        "speed_contract_ok": max(speed_errors) <= tolerance,
        "voltage_magnitude_ok": max(voltage_errors) <= tolerance,
        "current_magnitude_ok": max(current_errors) <= tolerance,
        "power_balance_ok": max(power_balance_errors) <= tolerance,
        "torque_speed_ok": max(torque_speed_errors) <= tolerance,
        "power_factor_formula_ok": max(power_factor_errors) <= tolerance,
        "power_factor_bounded": min(pf_values) >= -tolerance and max(pf_values) <= 1.0 + tolerance,
        "efficiency_formula_ok": max(efficiency_errors) <= tolerance if efficiency_errors else True,
        "speed_non_decreasing": all(
            speed_order[i] <= speed_order[i + 1] + tolerance
            for i in range(len(speed_order) - 1)
        ),
        "voltage_utilization_report_ok": max(voltage_utilization_errors) <= tolerance
        if voltage_utilization_errors else True,
        "lossless_voltage_constraint_ok": max_lossless_util <= 1.0 + tolerance
        if lossless_voltage_utils else True,
    }
    return {
        "policy": "pm_drive_terminal_table_health_gate",
        "row_count": len(row_summaries),
        "required_columns": list(required),
        "pole_pairs": pp,
        "tol": tolerance,
        "max_speed_contract_rel_error": max(speed_errors),
        "max_voltage_magnitude_abs_error_V": max(voltage_errors),
        "max_current_magnitude_abs_error_A": max(current_errors),
        "max_power_balance_rel_error": max(power_balance_errors),
        "max_torque_speed_rel_error": max(torque_speed_errors),
        "max_power_factor_abs_error": max(power_factor_errors),
        "max_efficiency_abs_error": max(efficiency_errors) if efficiency_errors else math.nan,
        "max_voltage_utilization_abs_error": max(voltage_utilization_errors) if voltage_utilization_errors else math.nan,
        "max_lossless_voltage_utilization": max_lossless_util,
        "max_terminal_voltage_utilization": max_terminal_util,
        "terminal_voltage_over_limit_row_count": len(terminal_over_limit_rows),
        "terminal_voltage_over_limit_policy": "report_only_use_lossless_voltage_utilization_for_selector_feasibility",
        "checks": checks,
        "status": "ok" if all(checks.values()) else "needs_attention",
        "rows": row_summaries,
    }


def pm_drive_loss_bucket_efficiency_gate(rows, loss_columns=None, tol=1.0e-9):
    """Check PM drive efficiency-map rows with explicit loss buckets.

    This is the public-safe contract for FE efficiency maps before they are
    handed to control notebooks: normalize tool-specific export columns into a
    small table, then verify ``P_in = P_out + sum(losses)`` and
    ``efficiency = P_out / P_in``.  If ``omega_mech`` and ``torque_Nm`` are
    present, ``P_out = torque*omega_mech`` is checked as an additional contract.
    """

    table = list(rows)
    if not table:
        raise ValueError("rows must not be empty")
    if loss_columns is None:
        jmag_style = ("P_cu_W", "P_iron_W", "P_magnet_W", "P_mechanical_loss_W")
        elf_style = ("copper_loss_w", "iron_loss_w", "magnet_loss_w", "mechanical_loss_w")
        if all(name in table[0] for name in jmag_style):
            loss_columns = jmag_style
        elif all(name in table[0] for name in elf_style):
            loss_columns = elf_style
        else:
            loss_columns = jmag_style
    loss_columns = tuple(str(name) for name in loss_columns)
    required = ("P_out_W", "P_in_W", "efficiency", *loss_columns)
    for index, row in enumerate(table):
        missing = [name for name in required if name not in row]
        if missing:
            raise ValueError(f"row {index} is missing required columns: {missing}")

    tolerance = float(tol)
    row_summaries = []
    power_balance_errors = []
    efficiency_errors = []
    torque_speed_errors = []
    loss_fractions = []
    efficiencies = []
    negative_loss_rows = []
    pin_positive = []
    pout_nonnegative = []

    for index, row in enumerate(table):
        p_out = float(row["P_out_W"])
        p_in = float(row["P_in_W"])
        efficiency = float(row["efficiency"])
        loss_values = {name: float(row[name]) for name in loss_columns}
        total_loss = sum(loss_values.values())
        scale = max(abs(p_in), abs(p_out), abs(total_loss), 1.0)
        power_balance_error = abs(p_in - p_out - total_loss) / scale
        efficiency_ref = p_out / p_in if p_in > 0.0 else math.nan
        efficiency_error = abs(efficiency - efficiency_ref) if math.isfinite(efficiency_ref) else math.nan
        loss_fraction = total_loss / p_in if p_in > 0.0 else math.nan
        dominant_loss_bucket = max(loss_values, key=lambda name: loss_values[name])
        negative_losses = [
            name for name, value in loss_values.items()
            if value < -tolerance
        ]
        if negative_losses:
            negative_loss_rows.append(index)

        torque_speed_error = math.nan
        if "omega_mech" in row and "torque_Nm" in row:
            p_out_ref = float(row["omega_mech"]) * float(row["torque_Nm"])
            torque_speed_error = abs(p_out - p_out_ref) / max(abs(p_out), abs(p_out_ref), 1.0)
            torque_speed_errors.append(torque_speed_error)

        power_balance_errors.append(power_balance_error)
        if math.isfinite(efficiency_error):
            efficiency_errors.append(efficiency_error)
        if math.isfinite(loss_fraction):
            loss_fractions.append(loss_fraction)
        efficiencies.append(efficiency)
        pin_positive.append(p_in > 0.0)
        pout_nonnegative.append(p_out >= -tolerance)

        row_summaries.append({
            "index": index,
            "operating_point": row.get("operating_point", row.get("region", "")),
            "speed_rpm": row.get("speed_rpm"),
            "omega_mech": row.get("omega_mech"),
            "torque_Nm": row.get("torque_Nm"),
            "P_out_W": p_out,
            "P_in_W": p_in,
            "losses_W": loss_values,
            "total_loss_W": total_loss,
            "loss_fraction": loss_fraction,
            "dominant_loss_bucket": dominant_loss_bucket,
            "efficiency": efficiency,
            "efficiency_ref": efficiency_ref,
            "power_balance_rel_error": power_balance_error,
            "efficiency_abs_error": efficiency_error,
            "torque_speed_rel_error": torque_speed_error,
            "negative_loss_buckets": negative_losses,
        })

    max_efficiency_row = max(row_summaries, key=lambda item: item["efficiency"])
    max_loss_fraction_row = max(row_summaries, key=lambda item: item["loss_fraction"])
    checks = {
        "required_columns_ok": True,
        "pin_positive": all(pin_positive),
        "pout_nonnegative": all(pout_nonnegative),
        "loss_buckets_nonnegative": not negative_loss_rows,
        "power_balance_ok": max(power_balance_errors) <= tolerance,
        "efficiency_formula_ok": max(efficiency_errors) <= tolerance if efficiency_errors else False,
        "efficiency_bounded": min(efficiencies) >= -tolerance and max(efficiencies) <= 1.0 + tolerance,
        "torque_speed_ok": max(torque_speed_errors) <= tolerance if torque_speed_errors else True,
    }
    return {
        "policy": "pm_drive_loss_bucket_efficiency_gate",
        "row_count": len(row_summaries),
        "required_columns": list(required),
        "loss_columns": list(loss_columns),
        "tol": tolerance,
        "max_power_balance_rel_error": max(power_balance_errors),
        "max_efficiency_abs_error": max(efficiency_errors) if efficiency_errors else math.nan,
        "max_torque_speed_rel_error": max(torque_speed_errors) if torque_speed_errors else math.nan,
        "max_efficiency": max_efficiency_row["efficiency"],
        "max_efficiency_row": max_efficiency_row,
        "max_loss_fraction": max(loss_fractions) if loss_fractions else math.nan,
        "max_loss_fraction_row": max_loss_fraction_row,
        "negative_loss_row_count": len(negative_loss_rows),
        "checks": checks,
        "status": "ok" if all(checks.values()) else "needs_attention",
        "rows": row_summaries,
    }


def drive_cycle_weighted_efficiency_gate(rows, tol=1.0e-9):
    """Compute weighted drive-cycle efficiency from operating-point rows.

    Rows may provide ``P_out_W``/``P_in_W`` directly, or ``P_out_W`` plus loss
    buckets such as ``P_cu_W`` and ``P_iron_W``.  Weights are normalized before
    scoring so an ELF/JMAG/MATLAB optimization table can use percentages or
    fractions without changing the result.
    """

    table = list(rows)
    tolerance = float(tol)
    if not table:
        raise ValueError("rows must not be empty")
    if tolerance < 0.0:
        raise ValueError("tol must be non-negative")

    loss_aliases = {
        "P_cu_W": ("P_cu_W", "copper_loss_w", "copper_loss_W"),
        "P_iron_W": ("P_iron_W", "iron_loss_w", "iron_loss_W"),
        "P_magnet_W": ("P_magnet_W", "magnet_loss_w", "magnet_loss_W"),
        "P_mechanical_loss_W": ("P_mechanical_loss_W", "mechanical_loss_w", "mechanical_loss_W"),
    }

    missing_default = object()

    def value(row, names, default=missing_default):
        for name in names:
            if name in row and row[name] is not None:
                return float(row[name])
        if default is missing_default:
            raise ValueError(f"row is missing one of {names!r}")
        return default

    weights = [value(row, ("weight", "cycle_weight", "duty_weight"), 1.0) for row in table]
    if any(weight < 0.0 for weight in weights):
        raise ValueError("weights must be non-negative")
    weight_sum = sum(weights)
    if weight_sum <= 0.0:
        raise ValueError("at least one weight must be positive")

    normalized_weights = [weight / weight_sum for weight in weights]
    rows_out = []
    weighted_output = 0.0
    weighted_input = 0.0
    weighted_losses = {name: 0.0 for name in loss_aliases}
    max_balance_error = 0.0
    max_efficiency_error = 0.0
    for index, (row, weight) in enumerate(zip(table, normalized_weights)):
        p_out = value(row, ("P_out_W", "p_out_w", "output_power_w", "mechanical_power_w"))
        loss_terms = {
            canonical: value(row, aliases, 0.0)
            for canonical, aliases in loss_aliases.items()
        }
        explicit_total_loss = value(row, ("total_loss_W", "total_loss_w", "loss_w"), None)
        total_loss = sum(loss_terms.values())
        if explicit_total_loss is not None and total_loss == 0.0:
            total_loss = explicit_total_loss
        p_in = value(row, ("P_in_W", "p_in_w", "input_power_w"), p_out + total_loss)
        efficiency = p_out / p_in if p_in > 0.0 else math.inf
        reported_efficiency = value(row, ("efficiency", "eta"), efficiency)
        balance_error = abs(p_in - (p_out + total_loss))
        balance_rel_error = balance_error / max(abs(p_in), abs(p_out + total_loss), 1.0e-300)
        efficiency_error = abs(reported_efficiency - efficiency)
        max_balance_error = max(max_balance_error, balance_rel_error)
        max_efficiency_error = max(max_efficiency_error, efficiency_error)
        weighted_output += weight * p_out
        weighted_input += weight * p_in
        for key, loss in loss_terms.items():
            weighted_losses[key] += weight * loss
        rows_out.append({
            "index": index,
            "point_id": str(row.get("point_id", index)),
            "normalized_weight": weight,
            "P_out_W": p_out,
            "P_in_W": p_in,
            "total_loss_W": total_loss,
            "efficiency": efficiency,
            "reported_efficiency": reported_efficiency,
            "power_balance_rel_error": balance_rel_error,
            "efficiency_abs_error": efficiency_error,
        })

    cycle_efficiency = weighted_output / weighted_input if weighted_input > 0.0 else math.inf
    weighted_total_loss = weighted_input - weighted_output
    dominant_loss_bucket = max(weighted_losses, key=lambda key: weighted_losses[key])
    worst_efficiency_row = min(rows_out, key=lambda row: row["efficiency"])
    checks = {
        "weights_sum_positive": weight_sum > 0.0,
        "weighted_input_positive": weighted_input > 0.0,
        "all_efficiency_in_0_1": all(0.0 <= row["efficiency"] <= 1.0 + tolerance for row in rows_out),
        "power_balance_ok": max_balance_error <= tolerance,
        "reported_efficiency_ok": max_efficiency_error <= tolerance,
    }
    return {
        "policy": "drive_cycle_weighted_efficiency_gate",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "row_count": len(rows_out),
        "input_weight_sum": weight_sum,
        "weighted_output_W": weighted_output,
        "weighted_input_W": weighted_input,
        "weighted_total_loss_W": weighted_total_loss,
        "cycle_efficiency": cycle_efficiency,
        "weighted_losses_W": weighted_losses,
        "dominant_weighted_loss_bucket": dominant_loss_bucket,
        "worst_efficiency_row": worst_efficiency_row,
        "max_power_balance_rel_error": max_balance_error,
        "max_efficiency_abs_error": max_efficiency_error,
        "checks": checks,
        "rows": rows_out,
        "tol": tolerance,
    }


def pm_recoil_demag_step_summary(step_fields_A_per_m, H_knee_A_per_m):
    """Summarize a three-step PM recoil-demag teaching contract.

    The expected step order is the irreversible-demag workflow:
    step 0 nominal field, step 1 opposing field applied, step 2 opposing field
    removed and returning on a recoil line.  The gate is solver-independent:
    if step 1 crosses below the knee, step 2 must carry a reduced remanence
    ratio instead of claiming full recovery.
    """

    required = (0, 1, 2)
    fields = {int(key): float(value) for key, value in dict(step_fields_A_per_m).items()}
    missing = [step for step in required if step not in fields]
    if missing:
        raise ValueError(f"step_fields_A_per_m is missing steps: {missing}")
    knee = float(H_knee_A_per_m)
    margins = {step: fields[step] - knee for step in required}
    crossed = fields[1] < knee
    recovered_to_nominal_side = fields[2] >= fields[0]
    if crossed:
        # A compact proxy for teaching and linting: full remanence at or above
        # the knee, linearly reduced in proportion to the step-1 overshoot.
        denom = max(abs(knee), 1.0e-300)
        remanence_ratio = max(0.0, 1.0 - (knee - fields[1]) / denom)
    else:
        remanence_ratio = 1.0
    checks = {
        "three_steps_present": True,
        "step1_is_worst_field": fields[1] <= fields[0] and fields[1] <= fields[2],
        "irreversible_flag_matches_knee": crossed == (margins[1] < 0.0),
        "step2_not_stronger_than_nominal_after_crossing": (not crossed) or not recovered_to_nominal_side,
        "remanence_ratio_in_unit_interval": 0.0 <= remanence_ratio <= 1.0,
    }
    return {
        "policy": "pm_recoil_demag_three_step_gate",
        "step_fields_A_per_m": {str(step): fields[step] for step in required},
        "H_knee_A_per_m": knee,
        "margins_A_per_m": {str(step): margins[step] for step in required},
        "irreversible_demag": crossed,
        "recoil_remanence_ratio_proxy": remanence_ratio,
        "checks": checks,
        "status": "ok" if all(checks.values()) else "needs_attention",
    }


def pm_loadline_metadata_gate(
    metadata,
    required_columns=(
        "temperature_C",
        "B_gap_T",
        "H_pm_A_per_m",
        "H_knee_A_per_m",
    ),
    *,
    h_field_unit=None,
    b_flux_density_unit=None,
    temperature_unit=None,
    field_sign_convention=None,
    magnetization_axis=None,
    knee_reference=None,
    recoil_mu_r=None,
):
    """Check PM load-line table metadata before reading demag values.

    This is a solver-independent preflight for demagnetization
    workflows and open motor examples.  It keeps the meaning of ``H`` and the
    knee field explicit before any table values are compared against
    ``pm_recoil_demag_step_summary`` or a temperature/load-line sweep.
    """

    if not isinstance(metadata, dict):
        raise ValueError("metadata must be a dictionary")
    columns_raw = metadata.get("columns", metadata.get("column_names", ()))
    columns = [str(column).strip() for column in columns_raw]
    if not columns:
        raise ValueError("metadata must include columns")
    column_set = set(columns)
    required = [str(column).strip() for column in required_columns]
    missing_required = [column for column in required if column not in column_set]

    h_unit = str(h_field_unit or metadata.get("h_field_unit") or "").strip().lower()
    h_unit = h_unit.replace("a/m", "a_per_m").replace("ka/m", "ka_per_m").replace(" ", "_")
    b_unit = str(b_flux_density_unit or metadata.get("b_flux_density_unit") or "").strip().lower()
    b_unit = b_unit.replace("tesla", "t").replace(" ", "_")
    temp_unit = str(temperature_unit or metadata.get("temperature_unit") or "").strip().lower()
    temp_unit = temp_unit.replace("degc", "c").replace("celsius", "c").replace("kelvin", "k")
    sign = str(field_sign_convention or metadata.get("field_sign_convention") or "").strip().lower()
    sign = sign.replace("-", "_").replace(" ", "_")
    axis = str(magnetization_axis or metadata.get("magnetization_axis") or "").strip().lower()
    axis = axis.replace("-", "_").replace(" ", "_")
    knee = str(knee_reference or metadata.get("knee_reference") or "").strip().lower()
    knee = knee.replace("-", "_").replace(" ", "_")
    mu_raw = recoil_mu_r if recoil_mu_r is not None else metadata.get("recoil_mu_r")
    mu = None if mu_raw is None else float(mu_raw)

    allowed_h_units = {"a_per_m", "ka_per_m"}
    allowed_b_units = {"t", "mt"}
    allowed_temperature_units = {"c", "k"}
    allowed_signs = {"negative_is_demag", "positive_is_demag"}
    allowed_axes = {"x", "y", "z", "radial", "tangential", "axial"}
    allowed_knees = {"intrinsic_bh_knee", "loadline_knee", "h_knee_field"}
    checks = {
        "required_columns_present": not missing_required,
        "h_field_unit_valid": h_unit in allowed_h_units,
        "b_flux_density_unit_valid": b_unit in allowed_b_units,
        "temperature_unit_valid": temp_unit in allowed_temperature_units,
        "field_sign_convention_valid": sign in allowed_signs,
        "magnetization_axis_valid": axis in allowed_axes,
        "knee_reference_valid": knee in allowed_knees,
        "recoil_mu_r_positive": mu is not None and mu > 0.0,
    }
    return {
        "policy": "pm_loadline_metadata_gate",
        "columns": columns,
        "required_columns": required,
        "missing_required_columns": missing_required,
        "h_field_unit": h_unit or None,
        "b_flux_density_unit": b_unit or None,
        "temperature_unit": temp_unit or None,
        "field_sign_convention": sign or None,
        "magnetization_axis": axis or None,
        "knee_reference": knee or None,
        "recoil_mu_r": mu,
        "checks": checks,
        "status": "ok" if all(checks.values()) else "needs_attention",
        "version_note": (
            "Use this before PM demag/load-line rows so H direction, units, "
            "knee reference, and recoil permeability are explicit."
        ),
    }


def pm_bem_surface_normal_metadata_gate(
    surface_rows,
    expected_normal_convention="outward_from_magnet",
    expected_surface_source_artifact_id=None,
    expected_magnetization_source_id=None,
    expected_material_id=None,
    expected_material_name=None,
    require_closed_charge_balance=True,
    tol=1.0e-12,
):
    """Check PM BEM surface-normal metadata before demag/source assembly.

    Magnetic-charge BEM uses ``sigma_m = M dot n``.  A surface row is not
    solver-ready unless it records outward-normal convention, positive area,
    unit normal, and magnetization direction.  For a closed PM with uniform
    magnetization, the signed surface-charge proxy must sum to zero.
    If expected identity fields are supplied, every row must also bind the
    equivalent-source artifact, magnetization source, and material identity.
    """

    rows_in = list(surface_rows)
    if not rows_in:
        raise ValueError("surface_rows must not be empty")
    tolerance = float(tol)
    if tolerance < 0.0:
        raise ValueError("tol must be non-negative")
    expected = str(expected_normal_convention).strip().lower().replace("-", "_").replace(" ", "_")
    normalized_rows = []
    name_counts = {}
    missing_convention = []
    bad_normals = []
    missing_magnetization = []
    nonpositive_area = []
    missing_identity = {
        "surface_source_artifact_id": [],
        "magnetization_source_id": [],
        "material_id": [],
        "material_name": [],
    }
    mismatched_identity = {
        "surface_source_artifact_id": [],
        "magnetization_source_id": [],
        "material_id": [],
        "material_name": [],
    }
    observed_identity = {
        "surface_source_artifact_id": set(),
        "magnetization_source_id": set(),
        "material_id": set(),
        "material_name": set(),
    }
    charge_sum = 0.0
    total_area = 0.0

    expected_identity = {
        "surface_source_artifact_id": expected_surface_source_artifact_id,
        "magnetization_source_id": expected_magnetization_source_id,
        "material_id": expected_material_id,
        "material_name": expected_material_name,
    }
    expected_identity = {
        key: str(value).strip() if value is not None else None
        for key, value in expected_identity.items()
    }

    def vector3(value):
        if value is None or len(value) != 3:
            return None
        return [float(value[0]), float(value[1]), float(value[2])]

    def first_text(row, names):
        for field in names:
            if field in row and row[field] is not None:
                value = str(row[field]).strip()
                if value:
                    return value
        return None

    for index, row in enumerate(rows_in, start=1):
        if not isinstance(row, dict):
            raise ValueError("each surface row must be a dictionary")
        name = str(row.get("surface") or row.get("name") or row.get("label") or "").strip()
        row_label = name or f"row_{index}"
        area = float(row.get("area_m2", row.get("area", 0.0)))
        convention = str(row.get("normal_convention") or row.get("normal_orientation") or "").strip().lower()
        convention = convention.replace("-", "_").replace(" ", "_")
        normal = vector3(row.get("normal") or row.get("normal_vector"))
        magnetization = vector3(row.get("magnetization") or row.get("magnetization_vector") or row.get("m_hat"))
        if name:
            name_counts[name] = name_counts.get(name, 0) + 1
        if convention != expected:
            missing_convention.append(name or f"row_{index}")
        if area <= 0.0:
            nonpositive_area.append(name or f"row_{index}")
        normal_norm = math.inf
        if normal is None:
            bad_normals.append(name or f"row_{index}")
        else:
            normal_norm = math.sqrt(sum(value * value for value in normal))
            if abs(normal_norm - 1.0) > tolerance:
                bad_normals.append(name or f"row_{index}")
        magnetization_norm = math.inf
        if magnetization is None:
            missing_magnetization.append(name or f"row_{index}")
            m_hat = None
        else:
            magnetization_norm = math.sqrt(sum(value * value for value in magnetization))
            if magnetization_norm <= 0.0:
                missing_magnetization.append(name or f"row_{index}")
                m_hat = None
            else:
                m_hat = [value / magnetization_norm for value in magnetization]
        dot = None
        signed_proxy = None
        if normal is not None and m_hat is not None and area > 0.0 and normal_norm > 0.0:
            n_hat = [value / normal_norm for value in normal]
            dot = sum(a * b for a, b in zip(m_hat, n_hat))
            signed_proxy = area * dot
            charge_sum += signed_proxy
            total_area += area
        row_identity = {
            "surface_source_artifact_id": first_text(
                row,
                (
                    "surface_source_artifact_id",
                    "source_balance_artifact_id",
                    "bem_source_artifact_id",
                ),
            ),
            "magnetization_source_id": first_text(
                row,
                (
                    "magnetization_source_id",
                    "magnetization_artifact_id",
                    "magnetization_map_id",
                ),
            ),
            "material_id": first_text(row, ("material_id", "magnet_material_id", "pm_material_id")),
            "material_name": first_text(
                row,
                ("material_name", "magnet_material_name", "pm_material_name"),
            ),
        }
        for key, value in row_identity.items():
            if value:
                observed_identity[key].add(value)
            expected_value = expected_identity[key]
            if expected_value is None:
                continue
            if not value:
                missing_identity[key].append(row_label)
            elif value != expected_value:
                mismatched_identity[key].append(row_label)
        normalized_rows.append({
            "index": index,
            "surface": name,
            "area_m2": area,
            "normal_convention": convention or None,
            "normal": normal,
            "normal_norm": normal_norm,
            "magnetization_unit": m_hat,
            "magnetization_norm": magnetization_norm,
            "m_dot_n": dot,
            "signed_charge_proxy": signed_proxy,
            **row_identity,
        })

    duplicate_surfaces = sorted(name for name, count in name_counts.items() if count > 1)
    balance_scale = max(total_area, 1.0)
    checks = {
        "surface_names_present": all(row["surface"] for row in normalized_rows),
        "surface_names_unique": not duplicate_surfaces,
        "normal_convention_matches": not missing_convention,
        "areas_positive": not nonpositive_area,
        "normals_unit": not bad_normals,
        "magnetization_vectors_present": not missing_magnetization,
        "surface_source_artifact_id_matches_expected": not missing_identity["surface_source_artifact_id"]
        and not mismatched_identity["surface_source_artifact_id"],
        "magnetization_source_id_matches_expected": not missing_identity["magnetization_source_id"]
        and not mismatched_identity["magnetization_source_id"],
        "material_id_matches_expected": not missing_identity["material_id"]
        and not mismatched_identity["material_id"],
        "material_name_matches_expected": not missing_identity["material_name"]
        and not mismatched_identity["material_name"],
        "closed_surface_charge_balances": (abs(charge_sum) <= tolerance * balance_scale)
        if require_closed_charge_balance
        else True,
    }
    return {
        "policy": "pm_bem_surface_normal_metadata_gate",
        "n_rows": len(normalized_rows),
        "expected_normal_convention": expected,
        "expected_surface_source_artifact_id": expected_identity["surface_source_artifact_id"],
        "expected_magnetization_source_id": expected_identity["magnetization_source_id"],
        "expected_material_id": expected_identity["material_id"],
        "expected_material_name": expected_identity["material_name"],
        "require_closed_charge_balance": bool(require_closed_charge_balance),
        "duplicate_surfaces": duplicate_surfaces,
        "missing_or_wrong_convention_surfaces": sorted(missing_convention),
        "nonpositive_area_surfaces": sorted(nonpositive_area),
        "bad_normal_surfaces": sorted(bad_normals),
        "missing_magnetization_surfaces": sorted(missing_magnetization),
        "observed_surface_source_artifact_ids": sorted(observed_identity["surface_source_artifact_id"]),
        "observed_magnetization_source_ids": sorted(observed_identity["magnetization_source_id"]),
        "observed_material_ids": sorted(observed_identity["material_id"]),
        "observed_material_names": sorted(observed_identity["material_name"]),
        "missing_identity_surfaces": {
            key: sorted(value) for key, value in missing_identity.items()
        },
        "mismatched_identity_surfaces": {
            key: sorted(value) for key, value in mismatched_identity.items()
        },
        "signed_charge_proxy_sum": charge_sum,
        "signed_charge_proxy_abs": abs(charge_sum),
        "total_area_m2": total_area,
        "tol": tolerance,
        "rows": normalized_rows,
        "checks": checks,
        "status": "ok" if all(checks.values()) else "needs_attention",
        "version_note": (
            "Use this before PM magnetic-charge BEM assembly so surface normals, "
            "magnetization direction, outward convention, source artifact, "
            "magnetization map, and material identity are explicit."
        ),
    }


def pm_demag_package_identity_gate(
    artifacts,
    expected_case_id=None,
    expected_magnet_id=None,
    required_kinds=("run_result", "loadline_metadata", "bem_surface", "recoil_steps"),
):
    """Check that PM demag artifacts describe one magnet/case package."""

    rows_in = list(artifacts)
    if not rows_in:
        raise ValueError("artifacts must not be empty")

    def _norm(value):
        return str(value or "").strip().lower().replace("-", "_").replace(" ", "_")

    def _first(row, names):
        for name in names:
            if name in row and row[name] is not None:
                return row[name]
        return None

    required = tuple(_norm(kind) for kind in required_kinds)
    if not required:
        raise ValueError("required_kinds must not be empty")

    expected_policies = {
        "run_result": {
            "elf_python_run_result_parse",
            "elf_python_run_result_parse_path",
            "run_result_loadline_handoff",
        },
        "loadline_metadata": {"pm_loadline_metadata_gate"},
        "bem_surface": {"pm_bem_surface_normal_metadata_gate"},
        "recoil_steps": {"pm_recoil_demag_three_step_gate", "pm_recoil_demag_step_summary"},
    }
    required_run_columns = {"H_pm_A_per_m", "H_knee_A_per_m", "safe_against_knee"}

    details = []
    kind_counts = {}
    case_ids = []
    magnet_ids = []
    missing_case_id = []
    missing_magnet_id = []
    missing_paths = []
    unknown_kinds = []
    bad_upstream_status = []
    bad_upstream_policy = []
    missing_run_columns = []
    bad_recoil_steps = []

    for index, row in enumerate(rows_in, start=1):
        if not isinstance(row, dict):
            raise ValueError("each artifact must be a dictionary")
        kind = _norm(_first(row, ("kind", "artifact_kind", "type")))
        case_id = _first(row, ("case_id", "run_case_id", "demag_case_id"))
        magnet_id = _first(row, ("magnet_id", "material_id", "pm_id", "mid"))
        path = _first(row, ("path", "file", "artifact_path", "table_path"))
        gate_policy = _first(row, ("gate_policy", "policy", "validator"))
        gate_policy_norm = _norm(gate_policy)
        status = _first(row, ("status", "gate_status", "validation_status"))
        status_norm = _norm(status)
        pass_flag = bool(row.get("pass", False))

        if kind:
            kind_counts[kind] = kind_counts.get(kind, 0) + 1
        else:
            unknown_kinds.append(index)
        if kind and kind not in set(required) | set(expected_policies):
            unknown_kinds.append(kind)
        if not case_id:
            missing_case_id.append(index)
        else:
            case_ids.append(str(case_id))
        if not magnet_id:
            missing_magnet_id.append(index)
        else:
            magnet_ids.append(str(magnet_id))
        if not path:
            missing_paths.append(index)
        if not (pass_flag or status_norm in {"ok", "pass", "passed", "verified"}):
            bad_upstream_status.append({"index": index, "kind": kind, "status": status})
        if kind in expected_policies and gate_policy_norm not in expected_policies[kind]:
            bad_upstream_policy.append({
                "index": index,
                "kind": kind,
                "gate_policy": gate_policy,
                "expected": sorted(expected_policies[kind]),
            })
        if kind == "run_result":
            columns = set(str(name) for name in row.get("normalized_columns", row.get("columns", ())))
            missing = sorted(required_run_columns - columns)
            if missing:
                missing_run_columns.append({"index": index, "missing": missing})
        if kind == "recoil_steps":
            steps = [int(step) for step in row.get("steps", row.get("required_steps", ()))]
            if set(steps) != {0, 1, 2}:
                bad_recoil_steps.append({"index": index, "steps": steps})

        details.append({
            "index": index,
            "kind": kind,
            "case_id": None if case_id is None else str(case_id),
            "magnet_id": None if magnet_id is None else str(magnet_id),
            "path": path,
            "gate_policy": gate_policy,
            "status": status,
            "pass": pass_flag,
        })

    required_set = set(required)
    present_set = set(kind_counts)
    unique_case_ids = sorted(set(case_ids))
    unique_magnet_ids = sorted(set(magnet_ids))
    checks = {
        "required_kinds_present": required_set.issubset(present_set),
        "no_unknown_kinds": not unknown_kinds,
        "case_ids_present": not missing_case_id,
        "case_ids_unique": len(unique_case_ids) == 1,
        "magnet_ids_present": not missing_magnet_id,
        "magnet_ids_unique": len(unique_magnet_ids) == 1,
        "paths_present": not missing_paths,
        "upstream_gate_status_ok": not bad_upstream_status,
        "upstream_gate_policy_known": not bad_upstream_policy,
        "run_result_has_loadline_columns": not missing_run_columns,
        "recoil_steps_are_three_step": not bad_recoil_steps,
    }
    if expected_case_id is not None:
        checks["expected_case_id_matches"] = unique_case_ids == [str(expected_case_id)]
    if expected_magnet_id is not None:
        checks["expected_magnet_id_matches"] = unique_magnet_ids == [str(expected_magnet_id)]

    return {
        "policy": "pm_demag_package_identity_gate",
        "required_kinds": list(required),
        "present_kinds": dict(sorted(kind_counts.items())),
        "case_ids": unique_case_ids,
        "magnet_ids": unique_magnet_ids,
        "expected_case_id": None if expected_case_id is None else str(expected_case_id),
        "expected_magnet_id": None if expected_magnet_id is None else str(expected_magnet_id),
        "missing_case_id_rows": missing_case_id,
        "missing_magnet_id_rows": missing_magnet_id,
        "missing_path_rows": missing_paths,
        "unknown_kinds": unknown_kinds,
        "bad_upstream_status_rows": bad_upstream_status,
        "bad_upstream_policy_rows": bad_upstream_policy,
        "missing_run_result_loadline_columns": missing_run_columns,
        "bad_recoil_step_rows": bad_recoil_steps,
        "artifacts": details,
        "checks": checks,
        "status": "ok" if all(checks.values()) else "needs_attention",
        "version_note": (
            "Run after RunResult parsing, PM load-line metadata, BEM surface "
            "normal metadata, and recoil-step gates so demag notebooks cannot "
            "mix rows from different magnets or cases."
        ),
    }


def pm_demag_margin_screening_package_gate(
    artifacts,
    expected_case_id=None,
    expected_magnet_id=None,
    expected_temperature_c=None,
    expected_bem_normal_convention="outward_from_magnet",
    expected_bem_surface_mesh_id=None,
    expected_bem_surface_mesh_digest=None,
    expected_bem_surface_row_count=None,
    expected_bem_source_balance_artifact_id=None,
    expected_bem_source_balance_digest=None,
    expected_bem_source_convention=None,
    expected_field_probe_id=None,
    expected_field_probe_family=None,
    expected_observation_region_id=None,
    expected_observation_component=None,
    expected_field_axis_convention=None,
    expected_field_sign_convention=None,
    expected_field_probe_method=None,
    expected_averaging_rule=None,
    expected_field_probe_geometry_digest=None,
    expected_field_probe_point_xyz_m=None,
    expected_field_probe_output_artifact_id=None,
    expected_field_probe_output_digest=None,
    expected_material_state_artifact_id=None,
    expected_material_state_digest=None,
    expected_load_step_id=None,
    expected_fault_step_id=None,
    expected_demag_step_id=None,
    require_field_probe_identity=False,
    require_field_probe_output_artifact=False,
    required_kinds=(
        "loadline_sweep",
        "fault_current_screening",
        "bem_source_balance",
        "demag_package",
    ),
):
    """Check that PM demag-margin screening artifacts form one package.

    This gate sits one level above ``pm_demag_package_identity_gate``.  It keeps
    load-line sweep summaries, negative-Id/fault-current screening rows, BEM
    magnetic-charge balance metadata, and the final demag package tied to the
    same case, magnet, and temperature before a notebook or MCP tool promotes a
    demag-margin result.
    """

    rows_in = list(artifacts)
    if not rows_in:
        raise ValueError("artifacts must not be empty")

    def _norm(value):
        return str(value or "").strip().lower().replace("-", "_").replace(" ", "_")

    def _first(row, names):
        for name in names:
            if name in row and row[name] is not None:
                return row[name]
        return None

    def _coordinate_tuple(value):
        if value in (None, ""):
            return None
        if isinstance(value, dict):
            lower_keys = ("x", "y", "z")
            upper_keys = ("X", "Y", "Z")
            if all(key in value for key in lower_keys[:2]):
                coords = [value[key] for key in lower_keys if key in value]
            elif all(key in value for key in upper_keys[:2]):
                coords = [value[key] for key in upper_keys if key in value]
            else:
                raise ValueError("coordinate dictionaries must include x/y or X/Y")
        elif isinstance(value, str):
            coords = [item for item in re.split(r"[,;\s]+", value.strip()) if item]
        else:
            coords = list(value)
        if len(coords) not in (2, 3):
            raise ValueError("coordinates must contain two or three values")
        return tuple(float(coord) for coord in coords)

    required = tuple(_norm(kind) for kind in required_kinds)
    if not required:
        raise ValueError("required_kinds must not be empty")

    expected_policies = {
        "loadline_sweep": {
            "pm_temperature_demag_sweep_summary",
            "pm_loadline_metadata_gate",
            "pm_loadline_risk_summary",
        },
        "fault_current_screening": {
            "fault_current_demag_screening",
            "motor_fault_current_demag_screening",
            "short_circuit_operating_point",
            "field_weakening_speed_capability",
        },
        "bem_source_balance": {
            "pm_bem_surface_source_balance_gate",
            "pm_bem_surface_normal_metadata_gate",
        },
        "demag_package": {"pm_demag_package_identity_gate"},
    }
    required_demag_artifacts = {"run_result", "loadline_metadata", "bem_surface", "recoil_steps"}
    required_fault_observables = {"field_probe", "demag_margin_A_per_m", "recoil_remanence_ratio_proxy"}
    expected_bem_normal = _norm(expected_bem_normal_convention)
    expected_bem_surface_mesh = (
        None if expected_bem_surface_mesh_id is None else str(expected_bem_surface_mesh_id).strip()
    )
    expected_bem_surface_mesh_digest_value = (
        None
        if expected_bem_surface_mesh_digest is None
        else str(expected_bem_surface_mesh_digest).strip()
    )
    expected_bem_surface_count = (
        None if expected_bem_surface_row_count is None else int(expected_bem_surface_row_count)
    )
    expected_bem_source_balance_artifact = (
        None
        if expected_bem_source_balance_artifact_id is None
        else str(expected_bem_source_balance_artifact_id).strip()
    )
    expected_bem_source_balance_digest_value = (
        None
        if expected_bem_source_balance_digest is None
        else str(expected_bem_source_balance_digest).strip()
    )
    expected_bem_source_convention_value = (
        None if expected_bem_source_convention is None else _norm(expected_bem_source_convention)
    )
    expected_field_probe = (
        None if expected_field_probe_id is None else str(expected_field_probe_id).strip()
    )
    expected_probe_family = (
        None if expected_field_probe_family is None else _norm(expected_field_probe_family)
    )
    expected_region = (
        None
        if expected_observation_region_id is None
        else str(expected_observation_region_id).strip()
    )
    expected_component = (
        None if expected_observation_component is None else _norm(expected_observation_component)
    )
    expected_axis_convention = (
        None if expected_field_axis_convention is None else _norm(expected_field_axis_convention)
    )
    expected_sign_convention = (
        None if expected_field_sign_convention is None else _norm(expected_field_sign_convention)
    )
    expected_probe_method = (
        None if expected_field_probe_method is None else _norm(expected_field_probe_method)
    )
    expected_average = None if expected_averaging_rule is None else _norm(expected_averaging_rule)
    expected_probe_geometry_digest = (
        None
        if expected_field_probe_geometry_digest is None
        else str(expected_field_probe_geometry_digest).strip()
    )
    expected_probe_point = _coordinate_tuple(expected_field_probe_point_xyz_m)
    expected_probe_output_artifact = (
        None
        if expected_field_probe_output_artifact_id is None
        else str(expected_field_probe_output_artifact_id).strip()
    )
    expected_probe_output_digest = (
        None
        if expected_field_probe_output_digest is None
        else str(expected_field_probe_output_digest).strip()
    )
    expected_material_state_artifact = (
        None
        if expected_material_state_artifact_id is None
        else str(expected_material_state_artifact_id).strip()
    )
    expected_material_state_digest_value = (
        None
        if expected_material_state_digest is None
        else str(expected_material_state_digest).strip()
    )
    expected_load_step = (
        None if expected_load_step_id is None else str(expected_load_step_id).strip()
    )
    expected_fault_step = (
        None if expected_fault_step_id is None else str(expected_fault_step_id).strip()
    )
    expected_demag_step = (
        None if expected_demag_step_id is None else str(expected_demag_step_id).strip()
    )
    field_probe_output_required = bool(
        require_field_probe_output_artifact
        or expected_probe_output_artifact is not None
        or expected_probe_output_digest is not None
    )
    field_probe_identity_required = bool(
        require_field_probe_identity
        or field_probe_output_required
        or expected_field_probe is not None
        or expected_probe_family is not None
        or expected_region is not None
        or expected_component is not None
        or expected_axis_convention is not None
        or expected_sign_convention is not None
        or expected_probe_method is not None
        or expected_average is not None
        or expected_probe_geometry_digest is not None
        or expected_probe_point is not None
    )
    field_probe_method_required = expected_probe_method is not None
    field_probe_geometry_digest_required = expected_probe_geometry_digest is not None
    field_probe_point_required = expected_probe_point is not None
    material_state_artifact_required = expected_material_state_artifact is not None
    material_state_digest_required = expected_material_state_digest_value is not None
    bem_surface_mesh_digest_required = expected_bem_surface_mesh_digest_value is not None
    bem_surface_row_count_required = expected_bem_surface_count is not None
    bem_source_balance_digest_required = expected_bem_source_balance_digest_value is not None
    bem_source_convention_required = expected_bem_source_convention_value is not None
    load_step_required = expected_load_step is not None
    fault_step_required = expected_fault_step is not None
    demag_step_required = expected_demag_step is not None

    details = []
    kind_counts = {}
    case_ids = []
    magnet_ids = []
    temperatures = []
    missing_case_id = []
    missing_magnet_id = []
    missing_temperature = []
    missing_paths = []
    unknown_kinds = []
    bad_upstream_status = []
    bad_upstream_policy = []
    missing_loadline_margin = []
    missing_fault_direction = []
    missing_fault_observables = []
    missing_bem_tolerance = []
    missing_bem_surface_area = []
    missing_bem_normal_convention = []
    wrong_bem_normal_convention = []
    missing_bem_source_unit = []
    missing_bem_balance_value = []
    missing_bem_surface_mesh_id = []
    bem_surface_mesh_ids = []
    bem_surface_mesh_digests = []
    bem_surface_row_counts = []
    missing_bem_surface_mesh_digest = []
    missing_bem_surface_row_count = []
    missing_bem_source_balance_artifact_id = []
    bem_source_balance_artifact_ids = []
    bem_source_balance_digests = []
    bem_source_conventions = []
    missing_bem_source_balance_digest = []
    missing_bem_source_convention = []
    missing_demag_artifacts = []
    material_state_rows = []
    missing_material_state = []
    material_state_artifact_ids = []
    material_state_digests = []
    missing_material_state_artifact_id = []
    missing_material_state_digest = []
    load_step_ids = []
    fault_step_ids = []
    demag_step_ids = []
    missing_load_step_id = []
    missing_fault_step_id = []
    missing_demag_step_id = []
    field_probe_ids = []
    field_probe_families = []
    observation_region_ids = []
    observation_components = []
    field_axis_conventions = []
    field_sign_conventions = []
    field_probe_methods = []
    averaging_rules = []
    field_probe_geometry_digests = []
    field_probe_points_xyz_m = []
    field_probe_output_artifact_ids = []
    field_probe_output_digests = []
    field_probe_output_paths = []
    missing_field_probe_id = []
    missing_field_probe_family = []
    missing_observation_region_id = []
    missing_observation_component = []
    missing_field_axis_convention = []
    missing_field_sign_convention = []
    missing_field_probe_method = []
    missing_averaging_rule = []
    missing_field_probe_geometry_digest = []
    missing_field_probe_point_xyz_m = []
    missing_field_probe_output_artifact_id = []
    missing_field_probe_output_digest = []
    missing_field_probe_output_path = []

    def _material_state(row):
        state = row.get("material_state")
        source = state if isinstance(state, dict) else row
        br = _first(source, ("Br_T", "br_T", "br_hot_T", "remanence_T"))
        h_knee = _first(source, ("H_knee_A_per_m", "h_knee_hot_A_per_m"))
        recoil = _first(source, ("recoil_mu_r", "mu_rec"))
        if br is None and h_knee is None and recoil is None:
            return None
        if br is None or h_knee is None or recoil is None:
            return "incomplete"
        return (
            round(float(br), 12),
            round(float(h_knee), 6),
            round(float(recoil), 12),
        )

    def _material_state_identity(row):
        state = row.get("material_state")
        source = state if isinstance(state, dict) else {}
        artifact_id = _first(
            row,
            (
                "material_state_artifact_id",
                "material_artifact_id",
                "bh_curve_artifact_id",
                "hbrm_hbcn_artifact_id",
                "hbrm_hbcn_contract_artifact_id",
            ),
        )
        if artifact_id is None:
            artifact_id = _first(
                source,
                (
                    "material_state_artifact_id",
                    "material_artifact_id",
                    "bh_curve_artifact_id",
                    "hbrm_hbcn_artifact_id",
                    "hbrm_hbcn_contract_artifact_id",
                ),
            )
        digest = _first(
            row,
            (
                "material_state_digest",
                "material_state_sha256",
                "material_digest",
                "material_sha256",
                "bh_curve_digest",
                "bh_curve_sha256",
                "hbrm_hbcn_digest",
                "hbrm_hbcn_sha256",
                "hbrm_hbcn_contract_digest",
            ),
        )
        if digest is None:
            digest = _first(
                source,
                (
                    "material_state_digest",
                    "material_state_sha256",
                    "material_digest",
                    "material_sha256",
                    "bh_curve_digest",
                    "bh_curve_sha256",
                    "hbrm_hbcn_digest",
                    "hbrm_hbcn_sha256",
                    "hbrm_hbcn_contract_digest",
                ),
            )
        return {
            "material_state_artifact_id": (
                None if artifact_id is None else str(artifact_id).strip()
            ),
            "material_state_digest": None if digest is None else str(digest).strip(),
        }

    def _step_id(row, names):
        value = _first(row, names)
        if value is None:
            for container in ("step_identity", "step", "solver_step", "workflow_step"):
                source = row.get(container)
                if isinstance(source, dict):
                    value = _first(source, names)
                    if value is not None:
                        break
        if value is None:
            return None
        value = str(value).strip()
        return value or None

    def _fault_probe_identity(row):
        probe = row.get("field_probe")
        probe_source = probe if isinstance(probe, dict) else {}

        probe_id = _first(row, ("field_probe_id", "field_probe_artifact_id", "observation_id", "probe_id"))
        if probe_id is None:
            probe_id = _first(probe_source, ("field_probe_id", "field_probe_artifact_id", "id", "probe_id"))
        if probe_id is None and probe is not None and not isinstance(probe, dict):
            probe_id = probe

        probe_family = _first(
            row,
            (
                "field_probe_family",
                "field_probe_observable_family",
                "probe_family",
                "observation_family",
                "observable_family",
            ),
        )
        if probe_family is None:
            probe_family = _first(
                probe_source,
                (
                    "field_probe_family",
                    "field_probe_observable_family",
                    "probe_family",
                    "observation_family",
                    "observable_family",
                ),
            )

        region_id = _first(
            row,
            (
                "observation_region_id",
                "field_probe_region_id",
                "probe_region_id",
                "region_id",
            ),
        )
        if region_id is None:
            region_id = _first(
                probe_source,
                (
                    "observation_region_id",
                    "field_probe_region_id",
                    "probe_region_id",
                    "region_id",
                ),
            )

        component = _first(
            row,
            ("field_component", "observation_component", "probe_component", "component"),
        )
        if component is None:
            component = _first(
                probe_source,
                ("field_component", "observation_component", "probe_component", "component"),
            )

        axis_convention = _first(
            row,
            (
                "field_axis_convention",
                "observation_axis_convention",
                "probe_axis_convention",
                "axis_convention",
            ),
        )
        if axis_convention is None:
            axis_convention = _first(
                probe_source,
                (
                    "field_axis_convention",
                    "observation_axis_convention",
                    "probe_axis_convention",
                    "axis_convention",
                ),
            )

        sign_convention = _first(
            row,
            (
                "field_sign_convention",
                "observation_sign_convention",
                "probe_sign_convention",
                "sign_convention",
            ),
        )
        if sign_convention is None:
            sign_convention = _first(
                probe_source,
                (
                    "field_sign_convention",
                    "observation_sign_convention",
                    "probe_sign_convention",
                    "sign_convention",
                ),
            )

        probe_method = _first(
            row,
            (
                "field_probe_method",
                "probe_method",
                "field_sampling_method",
                "sampling_method",
                "extraction_method",
            ),
        )
        if probe_method is None:
            probe_method = _first(
                probe_source,
                (
                    "field_probe_method",
                    "probe_method",
                    "field_sampling_method",
                    "sampling_method",
                    "extraction_method",
                ),
            )

        averaging = _first(
            row,
            (
                "averaging_rule",
                "field_averaging_rule",
                "probe_averaging_rule",
                "spatial_averaging",
            ),
        )
        if averaging is None:
            averaging = _first(
                probe_source,
                (
                    "averaging_rule",
                    "field_averaging_rule",
                    "probe_averaging_rule",
                    "spatial_averaging",
                ),
            )

        geometry_digest = _first(
            row,
            (
                "field_probe_geometry_digest",
                "field_probe_region_geometry_digest",
                "probe_geometry_digest",
                "observation_geometry_digest",
                "geometry_digest",
            ),
        )
        if geometry_digest is None:
            geometry_digest = _first(
                probe_source,
                (
                    "field_probe_geometry_digest",
                    "field_probe_region_geometry_digest",
                    "probe_geometry_digest",
                    "observation_geometry_digest",
                    "geometry_digest",
                ),
            )

        probe_point = _first(
            row,
            (
                "field_probe_point_xyz_m",
                "field_probe_location_xyz_m",
                "probe_point_xyz_m",
                "probe_location_xyz_m",
                "observation_point_xyz_m",
                "representative_point_xyz_m",
                "centroid_xyz_m",
            ),
        )
        if probe_point is None:
            probe_point = _first(
                probe_source,
                (
                    "field_probe_point_xyz_m",
                    "field_probe_location_xyz_m",
                    "probe_point_xyz_m",
                    "probe_location_xyz_m",
                    "observation_point_xyz_m",
                    "representative_point_xyz_m",
                    "centroid_xyz_m",
                ),
            )

        output_artifact_id = _first(
            row,
            (
                "field_probe_output_artifact_id",
                "probe_output_artifact_id",
                "field_probe_table_artifact_id",
                "table_artifact_id",
                "output_artifact_id",
                "observation_artifact_id",
            ),
        )
        if output_artifact_id is None:
            output_artifact_id = _first(
                probe_source,
                (
                    "field_probe_output_artifact_id",
                    "probe_output_artifact_id",
                    "field_probe_table_artifact_id",
                    "table_artifact_id",
                    "output_artifact_id",
                    "observation_artifact_id",
                ),
            )

        output_digest = _first(
            row,
            (
                "field_probe_output_digest",
                "field_probe_output_sha256",
                "probe_output_digest",
                "probe_output_sha256",
                "field_probe_table_digest",
                "field_probe_table_sha256",
                "table_digest",
                "table_sha256",
                "output_digest",
                "output_sha256",
                "observation_digest",
            ),
        )
        if output_digest is None:
            output_digest = _first(
                probe_source,
                (
                    "field_probe_output_digest",
                    "field_probe_output_sha256",
                    "probe_output_digest",
                    "probe_output_sha256",
                    "field_probe_table_digest",
                    "field_probe_table_sha256",
                    "table_digest",
                    "table_sha256",
                    "output_digest",
                    "output_sha256",
                    "observation_digest",
                ),
            )

        output_path = _first(
            row,
            (
                "field_probe_output_path",
                "probe_output_path",
                "field_probe_table_path",
                "table_path",
                "output_path",
                "observation_path",
            ),
        )
        if output_path is None:
            output_path = _first(
                probe_source,
                (
                    "field_probe_output_path",
                    "probe_output_path",
                    "field_probe_table_path",
                    "table_path",
                    "output_path",
                    "observation_path",
                ),
            )

        return {
            "field_probe_id": None if probe_id is None else str(probe_id).strip(),
            "field_probe_family": None if probe_family is None else _norm(probe_family),
            "observation_region_id": None if region_id is None else str(region_id).strip(),
            "observation_component": None if component is None else _norm(component),
            "field_axis_convention": None if axis_convention is None else _norm(axis_convention),
            "field_sign_convention": None if sign_convention is None else _norm(sign_convention),
            "field_probe_method": None if probe_method is None else _norm(probe_method),
            "averaging_rule": None if averaging is None else _norm(averaging),
            "field_probe_geometry_digest": (
                None if geometry_digest is None else str(geometry_digest).strip()
            ),
            "field_probe_point_xyz_m": _coordinate_tuple(probe_point),
            "field_probe_output_artifact_id": (
                None if output_artifact_id is None else str(output_artifact_id).strip()
            ),
            "field_probe_output_digest": (
                None if output_digest is None else str(output_digest).strip()
            ),
            "field_probe_output_path": None if output_path is None else str(output_path).strip(),
        }

    for index, row in enumerate(rows_in, start=1):
        if not isinstance(row, dict):
            raise ValueError("each artifact must be a dictionary")
        kind = _norm(_first(row, ("kind", "artifact_kind", "type")))
        case_id = _first(row, ("case_id", "run_case_id", "demag_case_id"))
        magnet_id = _first(row, ("magnet_id", "material_id", "pm_id", "mid"))
        temperature = _first(row, ("temperature_C", "temperature_c", "temp_C", "temp_c"))
        path = _first(row, ("path", "file", "artifact_path", "table_path"))
        gate_policy = _first(row, ("gate_policy", "policy", "validator"))
        gate_policy_norm = _norm(gate_policy)
        status = _first(row, ("status", "gate_status", "validation_status"))
        status_norm = _norm(status)
        pass_flag = bool(row.get("pass", False))
        material_state = _material_state(row)
        material_identity = _material_state_identity(row)
        load_step_id = None
        fault_step_id = None
        demag_step_id = None

        if kind:
            kind_counts[kind] = kind_counts.get(kind, 0) + 1
        else:
            unknown_kinds.append(index)
        if kind and kind not in set(required) | set(expected_policies):
            unknown_kinds.append(kind)
        if not case_id:
            missing_case_id.append(index)
        else:
            case_ids.append(str(case_id))
        if not magnet_id:
            missing_magnet_id.append(index)
        else:
            magnet_ids.append(str(magnet_id))
        if temperature is None:
            missing_temperature.append(index)
        else:
            temperatures.append(float(temperature))
        if not path:
            missing_paths.append(index)
        if not (pass_flag or status_norm in {"ok", "pass", "passed", "verified"}):
            bad_upstream_status.append({"index": index, "kind": kind, "status": status})
        if kind in expected_policies and gate_policy_norm not in expected_policies[kind]:
            bad_upstream_policy.append({
                "index": index,
                "kind": kind,
                "gate_policy": gate_policy,
                "expected": sorted(expected_policies[kind]),
            })

        if kind == "loadline_sweep" and not any(
            row.get(name) is not None
            for name in ("first_unsafe_gap_m", "minimum_demag_margin_A_per_m", "risk_label")
        ):
            missing_loadline_margin.append(index)
        if kind == "loadline_sweep":
            load_step_id = _step_id(
                row,
                (
                    "load_step_id",
                    "loadline_step_id",
                    "loadline_sweep_step_id",
                    "load_case_step_id",
                    "load_step_artifact_id",
                    "step_id",
                ),
            )
            if load_step_id:
                load_step_ids.append(load_step_id)
            elif load_step_required:
                missing_load_step_id.append(index)
        if kind == "fault_current_screening":
            fault_step_id = _step_id(
                row,
                (
                    "fault_step_id",
                    "fault_current_step_id",
                    "fault_current_screening_step_id",
                    "fault_case_step_id",
                    "operating_step_id",
                    "step_id",
                ),
            )
            if fault_step_id:
                fault_step_ids.append(fault_step_id)
            elif fault_step_required:
                missing_fault_step_id.append(index)
            if row.get("negative_id_is_demag_direction") is not True:
                missing_fault_direction.append(index)
            observables = set(str(item) for item in row.get("recommended_observable_keys", ()))
            missing = sorted(required_fault_observables - observables)
            if missing:
                missing_fault_observables.append({"index": index, "missing": missing})
            probe_identity = _fault_probe_identity(row)
            field_probe_id = probe_identity["field_probe_id"]
            field_probe_family = probe_identity["field_probe_family"]
            observation_region_id = probe_identity["observation_region_id"]
            observation_component = probe_identity["observation_component"]
            field_axis_convention = probe_identity["field_axis_convention"]
            field_sign_convention = probe_identity["field_sign_convention"]
            field_probe_method = probe_identity["field_probe_method"]
            averaging_rule = probe_identity["averaging_rule"]
            geometry_digest = probe_identity["field_probe_geometry_digest"]
            probe_point = probe_identity["field_probe_point_xyz_m"]
            output_artifact_id = probe_identity["field_probe_output_artifact_id"]
            output_digest = probe_identity["field_probe_output_digest"]
            output_path = probe_identity["field_probe_output_path"]
            if field_probe_id:
                field_probe_ids.append(field_probe_id)
            elif field_probe_identity_required:
                missing_field_probe_id.append(index)
            if field_probe_family:
                field_probe_families.append(field_probe_family)
            elif field_probe_identity_required:
                missing_field_probe_family.append(index)
            if observation_region_id:
                observation_region_ids.append(observation_region_id)
            elif field_probe_identity_required:
                missing_observation_region_id.append(index)
            if observation_component:
                observation_components.append(observation_component)
            elif field_probe_identity_required:
                missing_observation_component.append(index)
            if field_axis_convention:
                field_axis_conventions.append(field_axis_convention)
            elif field_probe_identity_required:
                missing_field_axis_convention.append(index)
            if field_sign_convention:
                field_sign_conventions.append(field_sign_convention)
            elif field_probe_identity_required:
                missing_field_sign_convention.append(index)
            if field_probe_method:
                field_probe_methods.append(field_probe_method)
            elif field_probe_method_required:
                missing_field_probe_method.append(index)
            if averaging_rule:
                averaging_rules.append(averaging_rule)
            elif field_probe_identity_required:
                missing_averaging_rule.append(index)
            if geometry_digest:
                field_probe_geometry_digests.append(geometry_digest)
            elif field_probe_geometry_digest_required:
                missing_field_probe_geometry_digest.append(index)
            if probe_point is not None:
                field_probe_points_xyz_m.append(probe_point)
            elif field_probe_point_required:
                missing_field_probe_point_xyz_m.append(index)
            if output_artifact_id:
                field_probe_output_artifact_ids.append(output_artifact_id)
            elif field_probe_output_required:
                missing_field_probe_output_artifact_id.append(index)
            if output_digest:
                field_probe_output_digests.append(output_digest)
            elif field_probe_output_required:
                missing_field_probe_output_digest.append(index)
            if output_path:
                field_probe_output_paths.append(output_path)
            elif field_probe_output_required:
                missing_field_probe_output_path.append(index)
        if kind == "bem_source_balance":
            surface_mesh_id = _first(
                row,
                ("surface_mesh_id", "boundary_mesh_id", "bem_surface_mesh_id", "mesh_id"),
            )
            if not surface_mesh_id:
                missing_bem_surface_mesh_id.append(index)
            else:
                bem_surface_mesh_ids.append(str(surface_mesh_id).strip())
            surface_mesh_digest = _first(
                row,
                (
                    "surface_mesh_digest",
                    "surface_mesh_sha256",
                    "bem_surface_mesh_digest",
                    "bem_surface_mesh_sha256",
                    "surface_mesh_artifact_digest",
                    "surface_mesh_artifact_sha256",
                ),
            )
            if surface_mesh_digest:
                bem_surface_mesh_digests.append(str(surface_mesh_digest).strip())
            elif bem_surface_mesh_digest_required:
                missing_bem_surface_mesh_digest.append(index)
            surface_row_count = _first(
                row,
                (
                    "surface_row_count",
                    "surface_mesh_row_count",
                    "bem_surface_row_count",
                    "surface_count",
                    "n_surface_rows",
                    "n_rows",
                ),
            )
            try:
                surface_count_ok = surface_row_count is not None and int(surface_row_count) > 0
            except (TypeError, ValueError):
                surface_count_ok = False
            if surface_count_ok:
                bem_surface_row_counts.append(int(surface_row_count))
            elif bem_surface_row_count_required:
                missing_bem_surface_row_count.append(index)
            source_balance_artifact_id = _first(
                row,
                (
                    "source_balance_artifact_id",
                    "bem_source_balance_artifact_id",
                    "artifact_id",
                    "export_artifact_id",
                ),
            )
            if not source_balance_artifact_id:
                missing_bem_source_balance_artifact_id.append(index)
            else:
                bem_source_balance_artifact_ids.append(str(source_balance_artifact_id).strip())
            source_balance_digest = _first(
                row,
                (
                    "source_balance_digest",
                    "source_balance_sha256",
                    "bem_source_balance_digest",
                    "bem_source_balance_sha256",
                    "source_balance_artifact_digest",
                    "source_balance_artifact_sha256",
                ),
            )
            if source_balance_digest:
                bem_source_balance_digests.append(str(source_balance_digest).strip())
            elif bem_source_balance_digest_required:
                missing_bem_source_balance_digest.append(index)
            source_convention = _first(
                row,
                (
                    "source_convention",
                    "bem_source_convention",
                    "source_density_convention",
                    "equivalent_source_convention",
                    "operator_source_convention",
                    "surface_source_density",
                ),
            )
            if source_convention:
                bem_source_conventions.append(_norm(source_convention))
            elif bem_source_convention_required:
                missing_bem_source_convention.append(index)
            if row.get("signed_charge_balance_rel_tol") is None:
                missing_bem_tolerance.append(index)
            surface_area = _first(row, ("total_area_m2", "surface_area_m2", "bem_surface_area_m2"))
            try:
                surface_area_ok = surface_area is not None and float(surface_area) > 0.0
            except (TypeError, ValueError):
                surface_area_ok = False
            if not surface_area_ok:
                missing_bem_surface_area.append(index)
            normal_convention = _first(
                row,
                (
                    "normal_convention",
                    "normal_orientation",
                    "surface_normal_convention",
                    "expected_normal_convention",
                ),
            )
            if not normal_convention:
                missing_bem_normal_convention.append(index)
            elif expected_bem_normal and _norm(normal_convention) != expected_bem_normal:
                wrong_bem_normal_convention.append({
                    "index": index,
                    "normal_convention": str(normal_convention),
                    "expected": expected_bem_normal,
                })
            source_unit = _first(
                row,
                ("source_balance_unit", "source_unit", "signed_charge_unit", "source_quantity_unit"),
            )
            if not source_unit:
                missing_bem_source_unit.append(index)
            balance_value = _first(
                row,
                (
                    "signed_charge_balance_abs",
                    "signed_charge_proxy_abs",
                    "source_balance_abs",
                    "source_balance_residual_abs",
                ),
            )
            try:
                balance_value_ok = balance_value is not None and float(balance_value) >= 0.0
            except (TypeError, ValueError):
                balance_value_ok = False
            if not balance_value_ok:
                missing_bem_balance_value.append(index)
        if kind == "demag_package":
            demag_step_id = _step_id(
                row,
                (
                    "demag_step_id",
                    "demag_package_step_id",
                    "hbrm_hbcn_step_id",
                    "hbcn_step_id",
                    "recoil_demag_step_id",
                    "step_id",
                ),
            )
            if demag_step_id:
                demag_step_ids.append(demag_step_id)
            elif demag_step_required:
                missing_demag_step_id.append(index)
            required_artifacts = set(str(item) for item in row.get("required_artifacts", ()))
            missing = sorted(required_demag_artifacts - required_artifacts)
            if missing:
                missing_demag_artifacts.append({"index": index, "missing": missing})
        if material_state is not None:
            if material_state == "incomplete":
                missing_material_state.append(index)
            else:
                material_state_rows.append(material_state)
        material_state_artifact_id = material_identity["material_state_artifact_id"]
        material_state_digest = material_identity["material_state_digest"]
        if material_state_artifact_id:
            material_state_artifact_ids.append(material_state_artifact_id)
        elif material_state_artifact_required:
            missing_material_state_artifact_id.append(index)
        if material_state_digest:
            material_state_digests.append(material_state_digest)
        elif material_state_digest_required:
            missing_material_state_digest.append(index)

        detail = {
            "index": index,
            "kind": kind,
            "case_id": None if case_id is None else str(case_id),
            "magnet_id": None if magnet_id is None else str(magnet_id),
            "temperature_C": None if temperature is None else float(temperature),
            "material_state": None if material_state in (None, "incomplete") else {
                "Br_T": material_state[0],
                "H_knee_A_per_m": material_state[1],
                "recoil_mu_r": material_state[2],
            },
            "material_state_artifact_id": material_state_artifact_id,
            "material_state_digest": material_state_digest,
            "load_step_id": load_step_id,
            "fault_step_id": fault_step_id,
            "demag_step_id": demag_step_id,
            "path": path,
            "gate_policy": gate_policy,
            "status": status,
            "pass": pass_flag,
        }
        if kind == "fault_current_screening":
            detail.update(probe_identity)
        details.append(detail)

    required_set = set(required)
    present_set = set(kind_counts)
    unique_case_ids = sorted(set(case_ids))
    unique_magnet_ids = sorted(set(magnet_ids))
    unique_temperatures = sorted(set(temperatures))
    material_state_seen = bool(material_state_rows or missing_material_state)
    unique_material_states = sorted(set(material_state_rows))
    unique_material_state_artifact_ids = sorted(set(material_state_artifact_ids))
    unique_material_state_digests = sorted(set(material_state_digests))
    unique_load_step_ids = sorted(set(load_step_ids))
    unique_fault_step_ids = sorted(set(fault_step_ids))
    unique_demag_step_ids = sorted(set(demag_step_ids))
    unique_field_probe_ids = sorted(set(field_probe_ids))
    unique_field_probe_families = sorted(set(field_probe_families))
    unique_observation_region_ids = sorted(set(observation_region_ids))
    unique_observation_components = sorted(set(observation_components))
    unique_field_axis_conventions = sorted(set(field_axis_conventions))
    unique_field_sign_conventions = sorted(set(field_sign_conventions))
    unique_field_probe_methods = sorted(set(field_probe_methods))
    unique_averaging_rules = sorted(set(averaging_rules))
    unique_field_probe_geometry_digests = sorted(set(field_probe_geometry_digests))
    unique_field_probe_points_xyz_m = sorted(set(field_probe_points_xyz_m))
    unique_field_probe_output_artifact_ids = sorted(set(field_probe_output_artifact_ids))
    unique_field_probe_output_digests = sorted(set(field_probe_output_digests))
    unique_field_probe_output_paths = sorted(set(field_probe_output_paths))
    unique_bem_surface_mesh_ids = sorted(set(bem_surface_mesh_ids))
    unique_bem_surface_mesh_digests = sorted(set(bem_surface_mesh_digests))
    unique_bem_surface_row_counts = sorted(set(bem_surface_row_counts))
    unique_bem_source_balance_artifact_ids = sorted(set(bem_source_balance_artifact_ids))
    unique_bem_source_balance_digests = sorted(set(bem_source_balance_digests))
    unique_bem_source_conventions = sorted(set(bem_source_conventions))
    point_tolerance_m = 1.0e-9
    field_probe_point_matches_expected = True
    if expected_probe_point is not None:
        field_probe_point_matches_expected = (
            len(unique_field_probe_points_xyz_m) == 1
            and len(unique_field_probe_points_xyz_m[0]) == len(expected_probe_point)
            and all(
                abs(actual - expected) <= point_tolerance_m
                for actual, expected in zip(unique_field_probe_points_xyz_m[0], expected_probe_point)
            )
        )
    if material_state_seen:
        missing_material_state.extend(
            index
            for index, row in enumerate(rows_in, start=1)
            if _material_state(row) is None
        )
        missing_material_state = sorted(set(missing_material_state))
    checks = {
        "required_kinds_present": required_set.issubset(present_set),
        "no_unknown_kinds": not unknown_kinds,
        "case_ids_present": not missing_case_id,
        "case_ids_unique": len(unique_case_ids) == 1,
        "magnet_ids_present": not missing_magnet_id,
        "magnet_ids_unique": len(unique_magnet_ids) == 1,
        "temperatures_present": not missing_temperature,
        "temperatures_unique": len(unique_temperatures) == 1,
        "paths_present": not missing_paths,
        "upstream_gate_status_ok": not bad_upstream_status,
        "upstream_gate_policy_known": not bad_upstream_policy,
        "loadline_margin_summary_recorded": not missing_loadline_margin,
        "fault_current_demag_direction_recorded": not missing_fault_direction,
        "fault_current_observables_complete": not missing_fault_observables,
        "bem_balance_tolerance_recorded": not missing_bem_tolerance,
        "bem_surface_area_recorded": not missing_bem_surface_area,
        "bem_normal_convention_recorded": not missing_bem_normal_convention,
        "bem_normal_convention_matches_expected": not wrong_bem_normal_convention,
        "bem_source_unit_recorded": not missing_bem_source_unit,
        "bem_balance_value_recorded": not missing_bem_balance_value,
        "bem_surface_mesh_id_recorded": not missing_bem_surface_mesh_id,
        "bem_surface_mesh_id_unique_when_present": (
            not unique_bem_surface_mesh_ids or len(unique_bem_surface_mesh_ids) == 1
        ),
        "bem_surface_mesh_digest_recorded": (
            not bem_surface_mesh_digest_required or not missing_bem_surface_mesh_digest
        ),
        "bem_surface_mesh_digest_unique_when_present": (
            not unique_bem_surface_mesh_digests or len(unique_bem_surface_mesh_digests) == 1
        ),
        "bem_surface_row_count_recorded": (
            not bem_surface_row_count_required or not missing_bem_surface_row_count
        ),
        "bem_surface_row_count_unique_when_present": (
            not unique_bem_surface_row_counts or len(unique_bem_surface_row_counts) == 1
        ),
        "bem_source_balance_artifact_id_recorded": not missing_bem_source_balance_artifact_id,
        "bem_source_balance_artifact_id_unique_when_present": (
            not unique_bem_source_balance_artifact_ids
            or len(unique_bem_source_balance_artifact_ids) == 1
        ),
        "bem_source_balance_digest_recorded": (
            not bem_source_balance_digest_required or not missing_bem_source_balance_digest
        ),
        "bem_source_balance_digest_unique_when_present": (
            not unique_bem_source_balance_digests
            or len(unique_bem_source_balance_digests) == 1
        ),
        "bem_source_convention_recorded": (
            not bem_source_convention_required or not missing_bem_source_convention
        ),
        "bem_source_convention_unique_when_present": (
            not unique_bem_source_conventions or len(unique_bem_source_conventions) == 1
        ),
        "demag_package_artifacts_complete": not missing_demag_artifacts,
        "material_state_complete_when_present": not material_state_seen or not missing_material_state,
        "material_state_unique_when_present": not material_state_seen or len(unique_material_states) == 1,
        "material_state_artifact_id_recorded": (
            not material_state_artifact_required or not missing_material_state_artifact_id
        ),
        "material_state_artifact_id_unique_when_present": (
            not unique_material_state_artifact_ids
            or len(unique_material_state_artifact_ids) == 1
        ),
        "material_state_digest_recorded": (
            not material_state_digest_required or not missing_material_state_digest
        ),
        "material_state_digest_unique_when_present": (
            not unique_material_state_digests or len(unique_material_state_digests) == 1
        ),
        "load_step_id_recorded": not load_step_required or not missing_load_step_id,
        "load_step_id_unique_when_present": (
            not unique_load_step_ids or len(unique_load_step_ids) == 1
        ),
        "fault_step_id_recorded": not fault_step_required or not missing_fault_step_id,
        "fault_step_id_unique_when_present": (
            not unique_fault_step_ids or len(unique_fault_step_ids) == 1
        ),
        "demag_step_id_recorded": not demag_step_required or not missing_demag_step_id,
        "demag_step_id_unique_when_present": (
            not unique_demag_step_ids or len(unique_demag_step_ids) == 1
        ),
        "field_probe_id_recorded": not field_probe_identity_required or not missing_field_probe_id,
        "field_probe_family_recorded": (
            not field_probe_identity_required or not missing_field_probe_family
        ),
        "observation_region_id_recorded": (
            not field_probe_identity_required or not missing_observation_region_id
        ),
        "observation_component_recorded": (
            not field_probe_identity_required or not missing_observation_component
        ),
        "field_axis_convention_recorded": (
            not field_probe_identity_required or not missing_field_axis_convention
        ),
        "field_sign_convention_recorded": (
            not field_probe_identity_required or not missing_field_sign_convention
        ),
        "field_probe_method_recorded_when_expected": (
            not field_probe_method_required or not missing_field_probe_method
        ),
        "averaging_rule_recorded": not field_probe_identity_required or not missing_averaging_rule,
        "field_probe_geometry_digest_recorded": (
            not field_probe_geometry_digest_required or not missing_field_probe_geometry_digest
        ),
        "field_probe_point_xyz_recorded_when_expected": (
            not field_probe_point_required or not missing_field_probe_point_xyz_m
        ),
        "field_probe_output_artifact_id_recorded": (
            not field_probe_output_required or not missing_field_probe_output_artifact_id
        ),
        "field_probe_output_digest_recorded": (
            not field_probe_output_required or not missing_field_probe_output_digest
        ),
        "field_probe_output_path_recorded": (
            not field_probe_output_required or not missing_field_probe_output_path
        ),
    }
    if expected_case_id is not None:
        checks["expected_case_id_matches"] = unique_case_ids == [str(expected_case_id)]
    if expected_magnet_id is not None:
        checks["expected_magnet_id_matches"] = unique_magnet_ids == [str(expected_magnet_id)]
    if expected_temperature_c is not None:
        checks["expected_temperature_c_matches"] = unique_temperatures == [float(expected_temperature_c)]
    if expected_bem_surface_mesh is not None:
        checks["expected_bem_surface_mesh_id_matches"] = (
            unique_bem_surface_mesh_ids == [expected_bem_surface_mesh]
        )
    if expected_bem_surface_mesh_digest_value is not None:
        checks["expected_bem_surface_mesh_digest_matches"] = (
            unique_bem_surface_mesh_digests == [expected_bem_surface_mesh_digest_value]
        )
    if expected_bem_surface_count is not None:
        checks["expected_bem_surface_row_count_matches"] = (
            unique_bem_surface_row_counts == [expected_bem_surface_count]
        )
    if expected_bem_source_balance_artifact is not None:
        checks["expected_bem_source_balance_artifact_id_matches"] = (
            unique_bem_source_balance_artifact_ids == [expected_bem_source_balance_artifact]
        )
    if expected_bem_source_balance_digest_value is not None:
        checks["expected_bem_source_balance_digest_matches"] = (
            unique_bem_source_balance_digests == [expected_bem_source_balance_digest_value]
        )
    if expected_bem_source_convention_value is not None:
        checks["expected_bem_source_convention_matches"] = (
            unique_bem_source_conventions == [expected_bem_source_convention_value]
        )
    if expected_field_probe is not None:
        checks["expected_field_probe_id_matches"] = unique_field_probe_ids == [expected_field_probe]
    if expected_probe_family is not None:
        checks["expected_field_probe_family_matches"] = (
            unique_field_probe_families == [expected_probe_family]
        )
    if expected_region is not None:
        checks["expected_observation_region_id_matches"] = unique_observation_region_ids == [expected_region]
    if expected_component is not None:
        checks["expected_observation_component_matches"] = unique_observation_components == [expected_component]
    if expected_axis_convention is not None:
        checks["expected_field_axis_convention_matches"] = (
            unique_field_axis_conventions == [expected_axis_convention]
        )
    if expected_sign_convention is not None:
        checks["expected_field_sign_convention_matches"] = (
            unique_field_sign_conventions == [expected_sign_convention]
        )
    if expected_probe_method is not None:
        checks["expected_field_probe_method_matches"] = (
            unique_field_probe_methods == [expected_probe_method]
        )
    if expected_average is not None:
        checks["expected_averaging_rule_matches"] = unique_averaging_rules == [expected_average]
    if expected_probe_geometry_digest is not None:
        checks["expected_field_probe_geometry_digest_matches"] = (
            unique_field_probe_geometry_digests == [expected_probe_geometry_digest]
        )
    if expected_probe_point is not None:
        checks["expected_field_probe_point_xyz_matches"] = field_probe_point_matches_expected
    if expected_probe_output_artifact is not None:
        checks["expected_field_probe_output_artifact_id_matches"] = (
            unique_field_probe_output_artifact_ids == [expected_probe_output_artifact]
        )
    if expected_probe_output_digest is not None:
        checks["expected_field_probe_output_digest_matches"] = (
            unique_field_probe_output_digests == [expected_probe_output_digest]
        )
    if expected_material_state_artifact is not None:
        checks["expected_material_state_artifact_id_matches"] = (
            unique_material_state_artifact_ids == [expected_material_state_artifact]
        )
    if expected_material_state_digest_value is not None:
        checks["expected_material_state_digest_matches"] = (
            unique_material_state_digests == [expected_material_state_digest_value]
        )
    if expected_load_step is not None:
        checks["expected_load_step_id_matches"] = unique_load_step_ids == [expected_load_step]
    if expected_fault_step is not None:
        checks["expected_fault_step_id_matches"] = unique_fault_step_ids == [expected_fault_step]
    if expected_demag_step is not None:
        checks["expected_demag_step_id_matches"] = unique_demag_step_ids == [expected_demag_step]

    return {
        "policy": "pm_demag_margin_screening_package_gate",
        "required_kinds": list(required),
        "present_kinds": dict(sorted(kind_counts.items())),
        "case_ids": unique_case_ids,
        "magnet_ids": unique_magnet_ids,
        "temperatures_C": unique_temperatures,
        "field_probe_identity_required": field_probe_identity_required,
        "field_probe_method_required": field_probe_method_required,
        "field_probe_output_required": field_probe_output_required,
        "field_probe_ids": unique_field_probe_ids,
        "field_probe_families": unique_field_probe_families,
        "observation_region_ids": unique_observation_region_ids,
        "observation_components": unique_observation_components,
        "field_axis_conventions": unique_field_axis_conventions,
        "field_sign_conventions": unique_field_sign_conventions,
        "field_probe_methods": unique_field_probe_methods,
        "averaging_rules": unique_averaging_rules,
        "field_probe_geometry_digests": unique_field_probe_geometry_digests,
        "field_probe_points_xyz_m": [list(point) for point in unique_field_probe_points_xyz_m],
        "field_probe_output_artifact_ids": unique_field_probe_output_artifact_ids,
        "field_probe_output_digests": unique_field_probe_output_digests,
        "field_probe_output_paths": unique_field_probe_output_paths,
        "material_states": [
            {
                "Br_T": state[0],
                "H_knee_A_per_m": state[1],
                "recoil_mu_r": state[2],
            }
            for state in unique_material_states
        ],
        "material_state_artifact_ids": unique_material_state_artifact_ids,
        "material_state_digests": unique_material_state_digests,
        "load_step_ids": unique_load_step_ids,
        "fault_step_ids": unique_fault_step_ids,
        "demag_step_ids": unique_demag_step_ids,
        "expected_case_id": None if expected_case_id is None else str(expected_case_id),
        "expected_magnet_id": None if expected_magnet_id is None else str(expected_magnet_id),
        "expected_temperature_C": (
            None if expected_temperature_c is None else float(expected_temperature_c)
        ),
        "expected_bem_normal_convention": expected_bem_normal or None,
        "expected_bem_surface_mesh_id": expected_bem_surface_mesh,
        "expected_bem_surface_mesh_digest": expected_bem_surface_mesh_digest_value,
        "expected_bem_surface_row_count": expected_bem_surface_count,
        "bem_surface_mesh_ids": unique_bem_surface_mesh_ids,
        "bem_surface_mesh_digests": unique_bem_surface_mesh_digests,
        "bem_surface_row_counts": unique_bem_surface_row_counts,
        "expected_bem_source_balance_artifact_id": expected_bem_source_balance_artifact,
        "expected_bem_source_balance_digest": expected_bem_source_balance_digest_value,
        "expected_bem_source_convention": expected_bem_source_convention_value,
        "bem_source_balance_artifact_ids": unique_bem_source_balance_artifact_ids,
        "bem_source_balance_digests": unique_bem_source_balance_digests,
        "bem_source_conventions": unique_bem_source_conventions,
        "expected_field_probe_id": expected_field_probe,
        "expected_field_probe_family": expected_probe_family,
        "expected_observation_region_id": expected_region,
        "expected_observation_component": expected_component,
        "expected_field_axis_convention": expected_axis_convention,
        "expected_field_sign_convention": expected_sign_convention,
        "expected_field_probe_method": expected_probe_method,
        "expected_averaging_rule": expected_average,
        "expected_field_probe_geometry_digest": expected_probe_geometry_digest,
        "expected_field_probe_point_xyz_m": (
            list(expected_probe_point) if expected_probe_point is not None else None
        ),
        "expected_field_probe_output_artifact_id": expected_probe_output_artifact,
        "expected_field_probe_output_digest": expected_probe_output_digest,
        "expected_material_state_artifact_id": expected_material_state_artifact,
        "expected_material_state_digest": expected_material_state_digest_value,
        "expected_load_step_id": expected_load_step,
        "expected_fault_step_id": expected_fault_step,
        "expected_demag_step_id": expected_demag_step,
        "missing_case_id_rows": missing_case_id,
        "missing_magnet_id_rows": missing_magnet_id,
        "missing_temperature_rows": missing_temperature,
        "missing_path_rows": missing_paths,
        "unknown_kinds": unknown_kinds,
        "bad_upstream_status_rows": bad_upstream_status,
        "bad_upstream_policy_rows": bad_upstream_policy,
        "missing_loadline_margin_rows": missing_loadline_margin,
        "missing_fault_direction_rows": missing_fault_direction,
        "missing_fault_observable_rows": missing_fault_observables,
        "missing_bem_tolerance_rows": missing_bem_tolerance,
        "missing_bem_surface_area_rows": missing_bem_surface_area,
        "missing_bem_normal_convention_rows": missing_bem_normal_convention,
        "wrong_bem_normal_convention_rows": wrong_bem_normal_convention,
        "missing_bem_source_unit_rows": missing_bem_source_unit,
        "missing_bem_balance_value_rows": missing_bem_balance_value,
        "missing_bem_surface_mesh_id_rows": missing_bem_surface_mesh_id,
        "missing_bem_surface_mesh_digest_rows": missing_bem_surface_mesh_digest,
        "missing_bem_surface_row_count_rows": missing_bem_surface_row_count,
        "missing_bem_source_balance_artifact_id_rows": missing_bem_source_balance_artifact_id,
        "missing_bem_source_balance_digest_rows": missing_bem_source_balance_digest,
        "missing_bem_source_convention_rows": missing_bem_source_convention,
        "missing_demag_artifact_rows": missing_demag_artifacts,
        "missing_material_state_rows": missing_material_state,
        "missing_material_state_artifact_id_rows": missing_material_state_artifact_id,
        "missing_material_state_digest_rows": missing_material_state_digest,
        "missing_load_step_id_rows": missing_load_step_id,
        "missing_fault_step_id_rows": missing_fault_step_id,
        "missing_demag_step_id_rows": missing_demag_step_id,
        "missing_field_probe_id_rows": missing_field_probe_id,
        "missing_field_probe_family_rows": missing_field_probe_family,
        "missing_observation_region_id_rows": missing_observation_region_id,
        "missing_observation_component_rows": missing_observation_component,
        "missing_field_axis_convention_rows": missing_field_axis_convention,
        "missing_field_sign_convention_rows": missing_field_sign_convention,
        "missing_field_probe_method_rows": missing_field_probe_method,
        "missing_averaging_rule_rows": missing_averaging_rule,
        "missing_field_probe_geometry_digest_rows": missing_field_probe_geometry_digest,
        "missing_field_probe_point_xyz_rows": missing_field_probe_point_xyz_m,
        "missing_field_probe_output_artifact_id_rows": missing_field_probe_output_artifact_id,
        "missing_field_probe_output_digest_rows": missing_field_probe_output_digest,
        "missing_field_probe_output_path_rows": missing_field_probe_output_path,
        "artifacts": details,
        "checks": checks,
        "status": "ok" if all(checks.values()) else "needs_attention",
        "version_note": (
            "Run after PM load-line sweep, negative-Id/fault-current screening, "
            "BEM source balance with surface area, expected normal convention, source unit, "
            "source convention, source-balance artifact/digest, and balance residual "
            "metadata, and demag package gates so demag-margin "
            "panels cannot mix temperatures, magnets, material states, cases, or "
            "load/fault/demag step identities, field-probe family, observation "
            "axis/sign convention, probe method, and probe geometry/output-artifact "
            "identity.  For ELF/MAGIC demag screening, the material-state artifact "
            "id and digest should also travel with HBRM/HBCN-derived B-H/recoil "
            "state."
        ),
    }


def _float_vector(name, values):
    vector = [float(value) for value in values]
    if not vector:
        raise ValueError(f"{name} must not be empty")
    return vector


def _dense_matvec(matrix, vector):
    return [sum(float(a) * x for a, x in zip(row, vector)) for row in matrix]


def _dense_transpose_matvec(matrix, vector):
    n_cols = len(matrix[0])
    return [
        sum(float(row[col]) * value for row, value in zip(matrix, vector))
        for col in range(n_cols)
    ]


def trace_surface_mass_energy_gate(trace_matrix, surface_mass, fem_values, tol=1.0e-12):
    """Check that a FEM-to-boundary trace preserves surface-mass energy.

    ``trace_matrix`` maps volume/FEM nodal values to compact boundary/BEM
    values.  ``surface_mass`` is the boundary P1 mass matrix.  The identity is
    purely algebraic and public-safe:

    ``(T u)^T M_b (T u) == u^T T^T M_b T u``.
    """

    trace = [list(row) for row in trace_matrix]
    mass = [list(row) for row in surface_mass]
    u = _float_vector("fem_values", fem_values)
    tolerance = float(tol)
    if tolerance < 0.0:
        raise ValueError("tol must be non-negative")
    if not trace or not trace[0]:
        raise ValueError("trace_matrix must be a non-empty dense row list")
    n_bem = len(trace)
    n_fem = len(trace[0])
    if any(len(row) != n_fem for row in trace):
        raise ValueError("trace_matrix rows must have a consistent column count")
    if len(u) != n_fem:
        raise ValueError("fem_values length must match trace_matrix column count")
    if len(mass) != n_bem or any(len(row) != n_bem for row in mass):
        raise ValueError("surface_mass must be square with trace_matrix row count")

    boundary_values = _dense_matvec(trace, u)
    mass_boundary_values = _dense_matvec(mass, boundary_values)
    boundary_energy = sum(a * b for a, b in zip(boundary_values, mass_boundary_values))
    embedded_action = _dense_transpose_matvec(trace, mass_boundary_values)
    embedded_energy = sum(a * b for a, b in zip(u, embedded_action))
    ones_fem = [1.0] * n_fem
    traced_ones = _dense_matvec(trace, ones_fem)
    ones_bem = [1.0] * n_bem
    mass_ones = _dense_matvec(mass, ones_bem)
    surface_area_from_mass = sum(mass_ones)
    symmetry_errors = [
        abs(float(mass[i][j]) - float(mass[j][i]))
        for i in range(n_bem)
        for j in range(n_bem)
    ]
    traced_constant_errors = [abs(value - 1.0) for value in traced_ones]
    interior_columns = [
        index + 1
        for index in range(n_fem)
        if all(abs(float(row[index])) <= tolerance for row in trace)
    ]
    embedded_interior_action = [
        abs(embedded_action[index - 1])
        for index in interior_columns
    ]
    energy_abs_error = abs(boundary_energy - embedded_energy)
    energy_rel_error = energy_abs_error / max(abs(boundary_energy), abs(embedded_energy), 1.0e-300)
    checks = {
        "surface_mass_symmetric": max(symmetry_errors) <= tolerance,
        "constant_trace_ok": max(traced_constant_errors) <= tolerance,
        "energy_identity_ok": energy_abs_error <= tolerance or energy_rel_error <= tolerance,
        "interior_columns_have_zero_boundary_action": (
            max(embedded_interior_action) <= tolerance if embedded_interior_action else True
        ),
        "surface_area_positive": surface_area_from_mass > 0.0,
    }
    return {
        "policy": "trace_surface_mass_energy_gate",
        "trace_shape": [n_bem, n_fem],
        "surface_mass_shape": [n_bem, n_bem],
        "interior_fem_node_ids": interior_columns,
        "boundary_values": boundary_values,
        "embedded_boundary_action": embedded_action,
        "boundary_energy": boundary_energy,
        "embedded_energy": embedded_energy,
        "energy_abs_error": energy_abs_error,
        "energy_rel_error": energy_rel_error,
        "surface_area_from_mass": surface_area_from_mass,
        "max_surface_mass_symmetry_error": max(symmetry_errors),
        "max_traced_constant_error": max(traced_constant_errors),
        "max_interior_boundary_action": max(embedded_interior_action) if embedded_interior_action else 0.0,
        "checks": checks,
        "status": "ok" if all(checks.values()) else "needs_attention",
        "tol": tolerance,
    }


def mqs_coulomb_gauge_efield_postprocess_gate(
    artifacts,
    expected_case_id=None,
    expected_mesh_id=None,
    expected_frequency_Hz=None,
    required_kinds=(
        "mqs_solution",
        "coulomb_gauge",
        "spatial_potential",
        "electric_field",
        "validity_envelope",
    ),
    max_frequency_ratio_to_fullwave=0.1,
    frequency_rtol=1.0e-12,
):
    """Check an MQS A-phi to Coulomb-gauge E-field postprocess package.

    The gate encodes a readable teaching rule from low-frequency EM practice:
    an electric-field row derived after a magneto-quasistatic solve is not a
    standalone result.  It must carry the A-phi solve identity, Coulomb-gauge
    postprocess, surface-potential boundary source, E-field units, and a
    validity envelope against Darwin/full-wave regimes.
    """

    rows_in = list(artifacts)
    if not rows_in:
        raise ValueError("artifacts must not be empty")

    def _norm(value):
        return str(value or "").strip().lower().replace("-", "_").replace(" ", "_")

    def _first(row, names):
        for name in names:
            if name in row and row[name] is not None:
                return row[name]
        return None

    required = tuple(_norm(kind) for kind in required_kinds)
    if not required:
        raise ValueError("required_kinds must not be empty")
    tolerance = float(frequency_rtol)
    max_ratio = float(max_frequency_ratio_to_fullwave)
    if tolerance < 0.0:
        raise ValueError("frequency_rtol must be non-negative")
    if max_ratio < 0.0:
        raise ValueError("max_frequency_ratio_to_fullwave must be non-negative")

    expected_policies = {
        "mqs_solution": {"mqs_a_phi_solution", "a_phi_mqs_solution", "mqs_solution"},
        "coulomb_gauge": {"coulomb_gauge_postprocess", "coulomb_gauge_condition"},
        "spatial_potential": {"electrostatic_potential_postprocess", "spatial_potential_solve"},
        "electric_field": {"efield_gradient_recovery", "electric_field_postprocess"},
        "validity_envelope": {"mqs_darwin_fullwave_validity_envelope", "quasistatic_validity_envelope"},
    }

    details = []
    kind_counts = {}
    case_ids = []
    mesh_ids = []
    frequencies = []
    missing_case_id = []
    missing_mesh_id = []
    missing_frequency = []
    missing_paths = []
    unknown_kinds = []
    bad_upstream_status = []
    bad_upstream_policy = []
    missing_mqs_formulation = []
    bad_gauge_condition = []
    missing_surface_potential_bc = []
    bad_efield_unit = []
    bad_validity_envelope = []

    for index, row in enumerate(rows_in, start=1):
        if not isinstance(row, dict):
            raise ValueError("each artifact must be a dictionary")
        kind = _norm(_first(row, ("kind", "artifact_kind", "type")))
        case_id = _first(row, ("case_id", "run_case_id", "model_id"))
        mesh_id = _first(row, ("mesh_id", "vol_id", "grid_id"))
        frequency = _first(row, ("frequency_Hz", "frequency_hz", "freq_Hz", "freq_hz"))
        path = _first(row, ("path", "file", "artifact_path", "table_path"))
        gate_policy = _first(row, ("gate_policy", "policy", "validator"))
        gate_policy_norm = _norm(gate_policy)
        status = _first(row, ("status", "gate_status", "validation_status"))
        status_norm = _norm(status)
        pass_flag = bool(row.get("pass", False))

        if kind:
            kind_counts[kind] = kind_counts.get(kind, 0) + 1
        else:
            unknown_kinds.append(index)
        if kind and kind not in set(required) | set(expected_policies):
            unknown_kinds.append(kind)
        if not case_id:
            missing_case_id.append(index)
        else:
            case_ids.append(str(case_id))
        if not mesh_id:
            missing_mesh_id.append(index)
        else:
            mesh_ids.append(str(mesh_id))
        if frequency is None:
            missing_frequency.append(index)
        else:
            frequencies.append(float(frequency))
        if not path:
            missing_paths.append(index)
        if not (pass_flag or status_norm in {"ok", "pass", "passed", "verified"}):
            bad_upstream_status.append({"index": index, "kind": kind, "status": status})
        if kind in expected_policies and gate_policy_norm not in expected_policies[kind]:
            bad_upstream_policy.append({
                "index": index,
                "kind": kind,
                "gate_policy": gate_policy,
                "expected": sorted(expected_policies[kind]),
            })
        if kind == "mqs_solution" and _norm(row.get("formulation")) not in {"a_phi", "a_phi_mqs", "mqs_a_phi"}:
            missing_mqs_formulation.append(index)
        if kind == "coulomb_gauge" and _norm(row.get("gauge_condition")) not in {"coulomb", "div_a_zero", "divergence_a_zero"}:
            bad_gauge_condition.append(index)
        if kind == "spatial_potential" and not row.get("boundary_condition_source"):
            missing_surface_potential_bc.append(index)
        if kind == "electric_field" and _norm(row.get("e_unit", row.get("E_unit"))) not in {"v_per_m", "v/m"}:
            bad_efield_unit.append(index)
        if kind == "validity_envelope":
            ratio = row.get("frequency_ratio_to_fullwave_limit")
            dominant_inductive = row.get("dominant_inductive")
            reference = _norm(row.get("comparison_reference"))
            valid_reference = reference in {"darwin", "full_wave", "darwin_full_wave", "fullwave"}
            if ratio is None or float(ratio) > max_ratio or dominant_inductive is not True or not valid_reference:
                bad_validity_envelope.append({
                    "index": index,
                    "frequency_ratio_to_fullwave_limit": ratio,
                    "dominant_inductive": dominant_inductive,
                    "comparison_reference": row.get("comparison_reference"),
                })

        details.append({
            "index": index,
            "kind": kind,
            "case_id": None if case_id is None else str(case_id),
            "mesh_id": None if mesh_id is None else str(mesh_id),
            "frequency_Hz": None if frequency is None else float(frequency),
            "path": path,
            "gate_policy": gate_policy,
            "status": status,
            "pass": pass_flag,
        })

    unique_case_ids = sorted(set(case_ids))
    unique_mesh_ids = sorted(set(mesh_ids))
    max_frequency_rel_span = 0.0
    if frequencies:
        max_frequency_rel_span = abs(max(frequencies) - min(frequencies)) / max(max(abs(f) for f in frequencies), 1.0)
    checks = {
        "required_kinds_present": set(required).issubset(set(kind_counts)),
        "no_unknown_kinds": not unknown_kinds,
        "case_ids_present": not missing_case_id,
        "case_ids_unique": len(unique_case_ids) == 1,
        "mesh_ids_present": not missing_mesh_id,
        "mesh_ids_unique": len(unique_mesh_ids) == 1,
        "frequencies_present": not missing_frequency,
        "frequencies_match": bool(frequencies) and max_frequency_rel_span <= tolerance,
        "paths_present": not missing_paths,
        "upstream_gate_status_ok": not bad_upstream_status,
        "upstream_gate_policy_known": not bad_upstream_policy,
        "mqs_a_phi_formulation_recorded": not missing_mqs_formulation,
        "coulomb_gauge_condition_recorded": not bad_gauge_condition,
        "spatial_potential_bc_recorded": not missing_surface_potential_bc,
        "electric_field_unit_recorded": not bad_efield_unit,
        "validity_envelope_ok": not bad_validity_envelope,
    }
    if expected_case_id is not None:
        checks["expected_case_id_matches"] = unique_case_ids == [str(expected_case_id)]
    if expected_mesh_id is not None:
        checks["expected_mesh_id_matches"] = unique_mesh_ids == [str(expected_mesh_id)]
    if expected_frequency_Hz is not None:
        expected_f = float(expected_frequency_Hz)
        checks["expected_frequency_matches"] = (
            bool(frequencies)
            and max(abs(freq - expected_f) / max(abs(freq), abs(expected_f), 1.0) for freq in frequencies) <= tolerance
        )

    return {
        "policy": "mqs_coulomb_gauge_efield_postprocess_gate",
        "required_kinds": list(required),
        "present_kinds": dict(sorted(kind_counts.items())),
        "case_ids": unique_case_ids,
        "mesh_ids": unique_mesh_ids,
        "frequencies_Hz": sorted(set(frequencies)),
        "max_frequency_rel_span": max_frequency_rel_span,
        "expected_case_id": None if expected_case_id is None else str(expected_case_id),
        "expected_mesh_id": None if expected_mesh_id is None else str(expected_mesh_id),
        "expected_frequency_Hz": None if expected_frequency_Hz is None else float(expected_frequency_Hz),
        "missing_case_id_rows": missing_case_id,
        "missing_mesh_id_rows": missing_mesh_id,
        "missing_frequency_rows": missing_frequency,
        "missing_path_rows": missing_paths,
        "unknown_kinds": unknown_kinds,
        "bad_upstream_status_rows": bad_upstream_status,
        "bad_upstream_policy_rows": bad_upstream_policy,
        "missing_mqs_formulation_rows": missing_mqs_formulation,
        "bad_gauge_condition_rows": bad_gauge_condition,
        "missing_surface_potential_bc_rows": missing_surface_potential_bc,
        "bad_efield_unit_rows": bad_efield_unit,
        "bad_validity_envelope_rows": bad_validity_envelope,
        "artifacts": details,
        "checks": checks,
        "status": "ok" if all(checks.values()) else "needs_attention",
        "version_note": (
            "Use this when an MQS A-phi solve is postprocessed into spatial E "
            "fields.  Coulomb gauge, conductor-surface potential boundary data, "
            "E-field units, and the Darwin/full-wave validity envelope must "
            "travel with the result."
        ),
    }


def thermal_annulus_conductance_gate(
    inner_radius,
    outer_radius,
    conductivity,
    delta_temperature=1.0,
    length=1.0,
    probe_radius=None,
    measured_temperature=None,
    measured_conductance=None,
    rtol=1.0e-6,
    atol=0.0,
):
    """Check the closed-form radial heat-conduction annulus identity.

    Inner boundary temperature is ``delta_temperature`` and outer boundary
    temperature is zero.  The exact field is
    ``T(r)=DeltaT*log(b/r)/log(b/a)`` and the thermal conductance is
    ``G=2*pi*k*L/log(b/a)``.  The geometric-mean radius ``sqrt(a*b)`` is the
    useful teaching probe because the exact temperature there is ``DeltaT/2``.
    """

    a = float(inner_radius)
    b = float(outer_radius)
    k = float(conductivity)
    dT = float(delta_temperature)
    ell = float(length)
    tolerance_r = float(rtol)
    tolerance_a = float(atol)
    if a <= 0.0:
        raise ValueError("inner_radius must be > 0")
    if b <= a:
        raise ValueError("outer_radius must be greater than inner_radius")
    if k <= 0.0:
        raise ValueError("conductivity must be > 0")
    if ell <= 0.0:
        raise ValueError("length must be > 0")
    if tolerance_r < 0.0 or tolerance_a < 0.0:
        raise ValueError("rtol and atol must be non-negative")

    radius = math.sqrt(a * b) if probe_radius is None else float(probe_radius)
    if radius < a or radius > b:
        raise ValueError("probe_radius must lie inside the annulus")
    log_ratio = math.log(b / a)
    conductance = 2.0 * math.pi * k * ell / log_ratio
    temperature = dT * math.log(b / radius) / log_ratio
    heat_rate = conductance * dT
    radial_flux = heat_rate / (2.0 * math.pi * radius * ell)
    geometric_mean_temperature = dT / 2.0
    temp_abs_error = None
    temp_rel_error = None
    conductance_abs_error = None
    conductance_rel_error = None
    checks = {
        "positive_conductance": conductance > 0.0,
        "probe_temperature_between_boundaries": min(0.0, dT) - tolerance_a <= temperature <= max(0.0, dT) + tolerance_a,
        "geometric_mean_half_temperature": (
            abs(temperature - geometric_mean_temperature) <= max(tolerance_a, tolerance_r * max(1.0, abs(geometric_mean_temperature)))
            if probe_radius is None else True
        ),
    }
    if measured_temperature is not None:
        temp_abs_error = abs(float(measured_temperature) - temperature)
        temp_rel_error = temp_abs_error / max(abs(temperature), abs(float(measured_temperature)), 1.0e-300)
        checks["measured_temperature_ok"] = temp_abs_error <= tolerance_a or temp_rel_error <= tolerance_r
    if measured_conductance is not None:
        conductance_abs_error = abs(float(measured_conductance) - conductance)
        conductance_rel_error = conductance_abs_error / max(abs(conductance), abs(float(measured_conductance)), 1.0e-300)
        checks["measured_conductance_ok"] = conductance_abs_error <= tolerance_a or conductance_rel_error <= tolerance_r
    return {
        "policy": "thermal_annulus_radial_conductance_gate",
        "inner_radius_m": a,
        "outer_radius_m": b,
        "conductivity_W_per_m_K": k,
        "length_m": ell,
        "delta_temperature_K": dT,
        "probe_radius_m": radius,
        "log_radius_ratio": log_ratio,
        "temperature_at_probe_K": temperature,
        "geometric_mean_temperature_K": geometric_mean_temperature,
        "thermal_conductance_W_per_K": conductance,
        "heat_rate_W": heat_rate,
        "radial_heat_flux_W_per_m2": radial_flux,
        "temperature_abs_error_K": temp_abs_error,
        "temperature_rel_error": temp_rel_error,
        "conductance_abs_error_W_per_K": conductance_abs_error,
        "conductance_rel_error": conductance_rel_error,
        "checks": checks,
        "status": "ok" if all(checks.values()) else "needs_attention",
        "rtol": tolerance_r,
        "atol": tolerance_a,
    }


def thermal_layer_stack_conductance_gate(
    area,
    thicknesses,
    conductivities,
    delta_temperature,
    measured_conductance=None,
    measured_interface_temperatures=None,
    measured_heat_rate=None,
    rtol=1.0e-9,
    atol=1.0e-12,
):
    """Check a 1-D layered heat-conduction stack by series resistances.

    The layers are stacked normal to the heat-flow direction.  The heat rate is
    constant through the stack and the resistance of layer ``i`` is
    ``R_i = d_i / (k_i A)``.  This is the thermal counterpart of the layered
    dielectric COMSOL/radia gates: keep the interface temperatures visible
    before trusting a full heat-transfer model.
    """

    surface_area = float(area)
    ds = [float(value) for value in thicknesses]
    ks = [float(value) for value in conductivities]
    dT = float(delta_temperature)
    rel_tol = float(rtol)
    abs_tol = float(atol)
    if surface_area <= 0.0:
        raise ValueError("area must be > 0")
    if not ds:
        raise ValueError("thicknesses must not be empty")
    if len(ds) != len(ks):
        raise ValueError("thicknesses and conductivities must have the same length")
    if any(value <= 0.0 for value in ds):
        raise ValueError("all layer thicknesses must be > 0")
    if any(value <= 0.0 for value in ks):
        raise ValueError("all layer conductivities must be > 0")
    if rel_tol < 0.0 or abs_tol < 0.0:
        raise ValueError("rtol and atol must be non-negative")

    layer_resistances = [d / (k * surface_area) for d, k in zip(ds, ks)]
    total_resistance = sum(layer_resistances)
    conductance = 1.0 / total_resistance
    heat_rate = conductance * dT
    heat_flux = heat_rate / surface_area
    drops = [heat_rate * resistance for resistance in layer_resistances]
    interface_temperatures = []
    running_drop = 0.0
    for drop in drops[:-1]:
        running_drop += drop
        interface_temperatures.append(dT - running_drop)
    total_thickness = sum(ds)
    effective_conductivity = total_thickness / sum(d / k for d, k in zip(ds, ks))
    checks = {
        "positive_conductance": conductance > 0.0,
        "temperature_drops_sum_to_drive": abs(sum(drops) - dT) <= max(abs_tol, rel_tol * max(abs(dT), 1.0)),
        "interface_temperatures_between_boundaries": all(
            -abs_tol <= temperature <= dT + abs_tol
            for temperature in interface_temperatures
        ),
    }
    out = {
        "policy": "thermal_layer_stack_series_resistance_gate",
        "area_m2": surface_area,
        "thicknesses_m": ds,
        "conductivities_W_per_m_K": ks,
        "delta_temperature_K": dT,
        "layer_resistances_K_per_W": layer_resistances,
        "total_resistance_K_per_W": total_resistance,
        "thermal_conductance_W_per_K": conductance,
        "effective_conductivity_W_per_m_K": effective_conductivity,
        "heat_rate_W": heat_rate,
        "heat_flux_W_per_m2": heat_flux,
        "layer_temperature_drops_K": drops,
        "interface_temperatures_K": interface_temperatures,
        "checks": checks,
        "rtol": rel_tol,
        "atol": abs_tol,
    }

    if measured_conductance is not None:
        measured_g = float(measured_conductance)
        error = abs(measured_g - conductance)
        out["measured_conductance_W_per_K"] = measured_g
        out["measured_conductance_abs_error_W_per_K"] = error
        out["measured_conductance_rel_error"] = error / max(abs(conductance), 1.0e-300)
        checks["measured_conductance_ok"] = (
            error <= abs_tol or out["measured_conductance_rel_error"] <= rel_tol
        )
    if measured_heat_rate is not None:
        measured_q = float(measured_heat_rate)
        error = abs(measured_q - heat_rate)
        out["measured_heat_rate_W"] = measured_q
        out["measured_heat_rate_abs_error_W"] = error
        out["measured_heat_rate_rel_error"] = error / max(abs(heat_rate), 1.0e-300)
        checks["measured_heat_rate_ok"] = (
            error <= abs_tol or out["measured_heat_rate_rel_error"] <= rel_tol
        )
    if measured_interface_temperatures is not None:
        measured_t = [float(value) for value in measured_interface_temperatures]
        if len(measured_t) != len(interface_temperatures):
            raise ValueError("measured_interface_temperatures must have one value per internal interface")
        errors = [
            abs(measured - expected)
            for measured, expected in zip(measured_t, interface_temperatures)
        ]
        out["measured_interface_temperatures_K"] = measured_t
        out["measured_interface_temperature_abs_errors_K"] = errors
        checks["measured_interface_temperatures_ok"] = all(
            error <= max(abs_tol, rel_tol * max(abs(expected), 1.0))
            for error, expected in zip(errors, interface_temperatures)
        )

    out["status"] = "ok" if all(checks.values()) else "needs_attention"
    return out


def thermal_conduction_convection_robin_gate(
    area,
    thicknesses,
    conductivities,
    heat_transfer_coefficient,
    hot_temperature,
    ambient_temperature=0.0,
    measured_conductance=None,
    measured_heat_rate=None,
    measured_heat_flux=None,
    measured_surface_temperature=None,
    measured_heat_transfer_coefficient=None,
    rtol=1.0e-9,
    atol=1.0e-12,
):
    """Check a 1-D conduction stack terminated by a convective Robin boundary.

    The total resistance is the sum of the layer conduction resistances and the
    film resistance: ``R = sum(d_i/(k_i*A)) + 1/(h*A)``.  This keeps Robin
    boundary checks teachable before moving to a full heat-transfer model.
    """

    surface_area = float(area)
    ds = [float(value) for value in thicknesses]
    ks = [float(value) for value in conductivities]
    h = float(heat_transfer_coefficient)
    hot = float(hot_temperature)
    ambient = float(ambient_temperature)
    rel_tol = float(rtol)
    abs_tol = float(atol)
    if surface_area <= 0.0:
        raise ValueError("area must be > 0")
    if not ds:
        raise ValueError("thicknesses must not be empty")
    if len(ds) != len(ks):
        raise ValueError("thicknesses and conductivities must have the same length")
    if any(value <= 0.0 for value in ds):
        raise ValueError("all layer thicknesses must be > 0")
    if any(value <= 0.0 for value in ks):
        raise ValueError("all layer conductivities must be > 0")
    if h <= 0.0:
        raise ValueError("heat_transfer_coefficient must be > 0")
    if rel_tol < 0.0 or abs_tol < 0.0:
        raise ValueError("rtol and atol must be non-negative")

    delta_temperature = hot - ambient
    conduction_resistances = [d / (k * surface_area) for d, k in zip(ds, ks)]
    convection_resistance = 1.0 / (h * surface_area)
    total_resistance = sum(conduction_resistances) + convection_resistance
    conductance = 1.0 / total_resistance
    heat_rate = conductance * delta_temperature
    heat_flux = heat_rate / surface_area
    conduction_drops = [heat_rate * resistance for resistance in conduction_resistances]
    convection_drop = heat_rate * convection_resistance
    interface_temperatures = []
    running_drop = 0.0
    for drop in conduction_drops[:-1]:
        running_drop += drop
        interface_temperatures.append(hot - running_drop)
    surface_temperature = ambient + convection_drop
    all_drops = conduction_drops + [convection_drop]
    low = min(hot, ambient) - abs_tol
    high = max(hot, ambient) + abs_tol
    checks = {
        "positive_conductance": conductance > 0.0,
        "temperature_drops_sum_to_drive": abs(sum(all_drops) - delta_temperature)
        <= max(abs_tol, rel_tol * max(abs(delta_temperature), 1.0)),
        "surface_temperature_between_hot_and_ambient": low <= surface_temperature <= high,
        "interface_temperatures_between_boundaries": all(
            low <= temperature <= high for temperature in interface_temperatures
        ),
    }
    out = {
        "policy": "thermal_conduction_convection_robin_gate",
        "area_m2": surface_area,
        "thicknesses_m": ds,
        "conductivities_W_per_m_K": ks,
        "heat_transfer_coefficient_W_per_m2_K": h,
        "hot_temperature_K": hot,
        "ambient_temperature_K": ambient,
        "delta_temperature_K": delta_temperature,
        "conduction_resistances_K_per_W": conduction_resistances,
        "convection_resistance_K_per_W": convection_resistance,
        "total_resistance_K_per_W": total_resistance,
        "thermal_conductance_W_per_K": conductance,
        "heat_rate_W": heat_rate,
        "heat_flux_W_per_m2": heat_flux,
        "conduction_temperature_drops_K": conduction_drops,
        "convection_temperature_drop_K": convection_drop,
        "solid_interface_temperatures_K": interface_temperatures,
        "robin_surface_temperature_K": surface_temperature,
        "checks": checks,
        "rtol": rel_tol,
        "atol": abs_tol,
    }

    if measured_conductance is not None:
        measured_g = float(measured_conductance)
        error = abs(measured_g - conductance)
        out["measured_conductance_W_per_K"] = measured_g
        out["measured_conductance_abs_error_W_per_K"] = error
        out["measured_conductance_rel_error"] = error / max(abs(conductance), 1.0e-300)
        checks["measured_conductance_ok"] = (
            error <= abs_tol or out["measured_conductance_rel_error"] <= rel_tol
        )
    if measured_heat_rate is not None:
        measured_q = float(measured_heat_rate)
        error = abs(measured_q - heat_rate)
        out["measured_heat_rate_W"] = measured_q
        out["measured_heat_rate_abs_error_W"] = error
        out["measured_heat_rate_rel_error"] = error / max(abs(heat_rate), 1.0e-300)
        checks["measured_heat_rate_ok"] = (
            error <= abs_tol or out["measured_heat_rate_rel_error"] <= rel_tol
        )
    if measured_heat_flux is not None:
        measured_flux = float(measured_heat_flux)
        error = abs(measured_flux - heat_flux)
        out["measured_heat_flux_W_per_m2"] = measured_flux
        out["measured_heat_flux_abs_error_W_per_m2"] = error
        out["measured_heat_flux_rel_error"] = error / max(abs(heat_flux), 1.0e-300)
        checks["measured_heat_flux_ok"] = (
            error <= abs_tol or out["measured_heat_flux_rel_error"] <= rel_tol
        )
    if measured_surface_temperature is not None:
        measured_surface = float(measured_surface_temperature)
        error = abs(measured_surface - surface_temperature)
        out["measured_robin_surface_temperature_K"] = measured_surface
        out["measured_robin_surface_temperature_abs_error_K"] = error
        checks["measured_robin_surface_temperature_ok"] = (
            error <= max(abs_tol, rel_tol * max(abs(surface_temperature), 1.0))
        )
    if measured_heat_transfer_coefficient is not None:
        measured_h = float(measured_heat_transfer_coefficient)
        error = abs(measured_h - h)
        out["measured_heat_transfer_coefficient_W_per_m2_K"] = measured_h
        out["measured_heat_transfer_coefficient_abs_error_W_per_m2_K"] = error
        out["measured_heat_transfer_coefficient_rel_error"] = error / max(abs(h), 1.0e-300)
        checks["measured_heat_transfer_coefficient_ok"] = (
            error <= abs_tol or out["measured_heat_transfer_coefficient_rel_error"] <= rel_tol
        )

    out["status"] = "ok" if all(checks.values()) else "needs_attention"
    return out


def spherical_dirichlet_laplacian_eigen_gate(
    radius,
    eigenvalues,
    rtol_first=2.0e-3,
    rtol_l1=2.0e-3,
    atol=0.0,
):
    """Check scalar Dirichlet Laplace eigenvalues on a ball.

    The first radial mode of ``-Delta u = lambda u`` with ``u=0`` on a ball of
    radius ``R`` is ``lambda_1 = (pi/R)^2``.  The next family is the triply
    degenerate ``l=1`` mode at ``(4.493409457909064/R)^2``.  Coarse hex meshes
    may split the degeneracy, so this gate checks the first eigenvalue and the
    nearest reported higher eigenvalue rather than assuming a perfectly
    repeated triplet.
    """

    r = float(radius)
    if r <= 0.0:
        raise ValueError("radius must be > 0")
    vals = sorted(float(value) for value in eigenvalues)
    if len(vals) < 2:
        raise ValueError("at least two eigenvalues are required")
    if any(value <= 0.0 for value in vals):
        raise ValueError("eigenvalues must be positive")
    tol_first = float(rtol_first)
    tol_l1 = float(rtol_l1)
    tolerance_a = float(atol)
    if tol_first < 0.0 or tol_l1 < 0.0 or tolerance_a < 0.0:
        raise ValueError("tolerances must be non-negative")

    l0_root = math.pi
    l1_root = 4.493409457909064
    lambda_first = (l0_root / r) ** 2
    lambda_l1 = (l1_root / r) ** 2
    measured_first = vals[0]
    l1_candidates = vals[1:]
    measured_l1 = min(l1_candidates, key=lambda value: abs(value - lambda_l1))
    first_abs_error = abs(measured_first - lambda_first)
    l1_abs_error = abs(measured_l1 - lambda_l1)
    first_rel_error = first_abs_error / max(abs(lambda_first), abs(measured_first), 1.0e-300)
    l1_rel_error = l1_abs_error / max(abs(lambda_l1), abs(measured_l1), 1.0e-300)
    checks = {
        "positive_radius": r > 0.0,
        "eigenvalues_positive": all(value > 0.0 for value in vals),
        "first_mode_matches_ball": first_abs_error <= tolerance_a or first_rel_error <= tol_first,
        "l1_mode_matches_ball": l1_abs_error <= tolerance_a or l1_rel_error <= tol_l1,
        "l1_mode_above_first_mode": measured_l1 > measured_first,
    }
    return {
        "policy": "spherical_dirichlet_laplacian_eigen_gate",
        "radius_m": r,
        "eigenvalues": vals,
        "reference": {
            "first_radial_lambda": lambda_first,
            "first_l1_lambda": lambda_l1,
            "first_l1_root": l1_root,
        },
        "measured_first_lambda": measured_first,
        "measured_l1_lambda": measured_l1,
        "first_abs_error": first_abs_error,
        "first_rel_error": first_rel_error,
        "l1_abs_error": l1_abs_error,
        "l1_rel_error": l1_rel_error,
        "checks": checks,
        "status": "ok" if all(checks.values()) else "needs_attention",
        "rtol_first": tol_first,
        "rtol_l1": tol_l1,
        "atol": tolerance_a,
    }


def _project_box(point, lower, upper):
    return [min(max(x, lo), hi) for x, lo, hi in zip(point, lower, upper)]


def _least_squares_objective(matrix, rhs, point):
    residual = [value - target for value, target in zip(_dense_matvec(matrix, point), rhs)]
    return 0.5 * sum(value * value for value in residual)


def _least_squares_gradient(matrix, rhs, point):
    residual = [value - target for value, target in zip(_dense_matvec(matrix, point), rhs)]
    return _dense_transpose_matvec(matrix, residual)


def _box_kkt_residual(point, gradient, lower, upper, active_tol):
    residuals = []
    for x, g, lo, hi in zip(point, gradient, lower, upper):
        if x <= lo + active_tol:
            residuals.append(max(0.0, -g))
        elif x >= hi - active_tol:
            residuals.append(max(0.0, g))
        else:
            residuals.append(abs(g))
    return residuals


def box_projected_gradient_least_squares_gate(
    matrix,
    rhs,
    lower,
    upper,
    initial=None,
    step_size=None,
    max_iterations=200,
    tol=1.0e-12,
    active_tol=1.0e-10,
):
    """Solve a small box-constrained least-squares teaching gate.

    The objective is ``0.5*||A*x-b||^2`` with componentwise bounds.  This is a
    lightweight public check for MATLAB/Gypsilab optimization lessons: it keeps
    the projection step and KKT sign convention visible without depending on
    Optuna or a commercial solver.
    """

    rows = [list(row) for row in matrix]
    if not rows or not rows[0]:
        raise ValueError("matrix must be a non-empty dense row list")
    n_cols = len(rows[0])
    if any(len(row) != n_cols for row in rows):
        raise ValueError("matrix rows must have a consistent column count")
    b = _float_vector("rhs", rhs)
    lo = _float_vector("lower", lower)
    hi = _float_vector("upper", upper)
    if len(rows) != len(b):
        raise ValueError("matrix row count must match rhs length")
    if not (len(lo) == len(hi) == n_cols):
        raise ValueError("lower and upper must have one entry per matrix column")
    if any(lo_i > hi_i for lo_i, hi_i in zip(lo, hi)):
        raise ValueError("each lower bound must be <= upper bound")
    if initial is None:
        raw_initial = [0.0] * n_cols
    else:
        raw_initial = _float_vector("initial", initial)
        if len(raw_initial) != n_cols:
            raise ValueError("initial must have one entry per matrix column")
    if max_iterations < 0:
        raise ValueError("max_iterations must be >= 0")

    if step_size is None:
        frob_squared = sum(float(value) ** 2 for row in rows for value in row)
        step = 1.0 / frob_squared if frob_squared > 0.0 else 1.0
    else:
        step = float(step_size)
    if step <= 0.0:
        raise ValueError("step_size must be > 0")

    x = _project_box(raw_initial, lo, hi)
    objective_history = [_least_squares_objective(rows, b, x)]
    iterations = 0
    last_step_norm = 0.0
    for iterations in range(1, int(max_iterations) + 1):
        gradient = _least_squares_gradient(rows, b, x)
        raw_step = [value - step * g for value, g in zip(x, gradient)]
        next_x = _project_box(raw_step, lo, hi)
        last_step_norm = math.sqrt(sum((new - old) ** 2 for new, old in zip(next_x, x)))
        x = next_x
        objective_history.append(_least_squares_objective(rows, b, x))
        if last_step_norm <= float(tol):
            break

    gradient = _least_squares_gradient(rows, b, x)
    projected_step = _project_box([value - step * g for value, g in zip(x, gradient)], lo, hi)
    projected_gradient_residual = math.sqrt(
        sum((value - projected) ** 2 for value, projected in zip(x, projected_step))
    ) / step
    kkt_residuals = _box_kkt_residual(x, gradient, lo, hi, float(active_tol))
    active_lower = [value <= bound + float(active_tol) for value, bound in zip(x, lo)]
    active_upper = [value >= bound - float(active_tol) for value, bound in zip(x, hi)]
    decreases = [
        objective_history[i + 1] - objective_history[i]
        for i in range(len(objective_history) - 1)
    ]
    checks = {
        "bounds_ok": all(lo_i - active_tol <= x_i <= hi_i + active_tol for x_i, lo_i, hi_i in zip(x, lo, hi)),
        "kkt_residual_ok": max(kkt_residuals) <= float(tol),
        "projected_gradient_residual_ok": projected_gradient_residual <= float(tol),
        "objective_monotone": all(delta <= max(float(tol), 1.0e-14) for delta in decreases),
    }
    return {
        "policy": "box_projected_gradient_least_squares_gate",
        "matrix_shape": [len(rows), n_cols],
        "rhs": b,
        "lower": lo,
        "upper": hi,
        "initial_raw": raw_initial,
        "initial_projected": _project_box(raw_initial, lo, hi),
        "x": x,
        "gradient": gradient,
        "objective": objective_history[-1],
        "objective_history": objective_history,
        "step_size": step,
        "iterations": iterations,
        "last_step_norm": last_step_norm,
        "projected_gradient_residual": projected_gradient_residual,
        "kkt_residuals": kkt_residuals,
        "max_kkt_residual": max(kkt_residuals),
        "active_lower": active_lower,
        "active_upper": active_upper,
        "checks": checks,
        "status": "ok" if all(checks.values()) else "needs_attention",
        "tol": float(tol),
        "active_tol": float(active_tol),
    }


def geometric_integrator_energy_drift_gate(
    method_rows,
    required_geometric_methods=("symplectic_euler", "implicit_midpoint"),
    explicit_method="explicit_euler",
    max_geometric_rel_drift=5.0e-2,
    min_explicit_to_geometric_drift_ratio=10.0,
    tol=1.0e-12,
):
    """Check a small geometric time-integration energy-drift teaching table.

    The gate is intentionally solver-independent.  MATLAB/Gypsilab notebooks
    can run a harmonic-oscillator time integration and export plain method rows
    with energy drift metrics; this helper replays the contract without needing
    MATLAB.  The expected lesson is that geometric methods keep Hamiltonian
    energy bounded on a fixed time grid, while explicit Euler is a negative
    control with much larger drift.
    """

    def canonical_method(value):
        return str(value).strip().lower().replace("-", "_").replace(" ", "_")

    rows = []
    by_method = {}
    tolerance = float(tol)
    max_allowed = float(max_geometric_rel_drift)
    min_ratio = float(min_explicit_to_geometric_drift_ratio)
    required = tuple(canonical_method(method) for method in required_geometric_methods)
    explicit = canonical_method(explicit_method)

    if not required:
        raise ValueError("required_geometric_methods must not be empty")
    if max_allowed < 0.0:
        raise ValueError("max_geometric_rel_drift must be non-negative")
    if min_ratio < 0.0:
        raise ValueError("min_explicit_to_geometric_drift_ratio must be non-negative")
    if tolerance < 0.0:
        raise ValueError("tol must be non-negative")

    for raw in method_rows:
        if not isinstance(raw, dict):
            raise ValueError("each method row must be a dictionary")
        method = canonical_method(raw.get("method", ""))
        if not method:
            raise ValueError("each method row needs a method name")
        if method in by_method:
            raise ValueError(f"duplicate method row: {method}")
        if "energy_initial" in raw:
            energy_initial = float(raw["energy_initial"])
        elif "energy_initial_J" in raw:
            energy_initial = float(raw["energy_initial_J"])
        else:
            raise ValueError(f"{method} row needs energy_initial")
        max_rel_drift = float(raw["max_rel_energy_drift"])
        if not math.isfinite(energy_initial) or energy_initial <= 0.0:
            raise ValueError(f"{method} energy_initial must be finite and > 0")
        if not math.isfinite(max_rel_drift) or max_rel_drift < 0.0:
            raise ValueError(f"{method} max_rel_energy_drift must be finite and >= 0")

        final_rel_drift = None
        if "final_rel_energy_drift" in raw:
            final_rel_drift = float(raw["final_rel_energy_drift"])
        elif "energy_final" in raw:
            final_rel_drift = abs(float(raw["energy_final"]) - energy_initial) / energy_initial
        elif "energy_final_J" in raw:
            final_rel_drift = abs(float(raw["energy_final_J"]) - energy_initial) / energy_initial
        if final_rel_drift is not None:
            if not math.isfinite(final_rel_drift) or final_rel_drift < 0.0:
                raise ValueError(f"{method} final energy drift must be finite and >= 0")

        row = {
            "method": method,
            "energy_initial": energy_initial,
            "max_rel_energy_drift": max_rel_drift,
            "final_rel_energy_drift": final_rel_drift,
        }
        for key in ("steps", "step_size_s", "omega_rad_per_s"):
            if key in raw:
                row[key] = int(raw[key]) if key == "steps" else float(raw[key])
        rows.append(row)
        by_method[method] = row

    if not rows:
        raise ValueError("method_rows must not be empty")

    def same_numeric(key):
        values = [row.get(key) for row in rows if key in row]
        if len(values) != len(rows):
            return False, None
        if key == "steps":
            first = int(values[0])
            return all(int(value) == first for value in values), first
        first = float(values[0])
        scale = max(1.0, *(abs(float(value)) for value in values))
        return all(abs(float(value) - first) <= tolerance * scale for value in values), first

    steps_consistent, common_steps = same_numeric("steps")
    step_size_consistent, common_step_size = same_numeric("step_size_s")
    omega_consistent, common_omega = same_numeric("omega_rad_per_s")

    missing_required = [method for method in required if method not in by_method]
    geometric_rows = [by_method[method] for method in required if method in by_method]
    geometric_drifts = [row["max_rel_energy_drift"] for row in geometric_rows]
    worst_geometric = max(geometric_drifts) if geometric_drifts else None
    best_geometric = min(geometric_drifts) if geometric_drifts else None
    best_geometric_method = None
    if geometric_rows:
        best_geometric_method = min(geometric_rows, key=lambda row: row["max_rel_energy_drift"])["method"]

    explicit_row = by_method.get(explicit)
    explicit_to_worst_ratio = None
    if explicit_row is not None and worst_geometric is not None:
        explicit_to_worst_ratio = explicit_row["max_rel_energy_drift"] / max(worst_geometric, tolerance, 1.0e-300)

    checks = {
        "required_geometric_methods_present": not missing_required,
        "geometric_energy_drift_bounded": bool(geometric_rows)
        and all(row["max_rel_energy_drift"] <= max_allowed + tolerance for row in geometric_rows),
        "explicit_negative_control_present": explicit_row is not None,
        "explicit_drift_larger_than_geometric": explicit_to_worst_ratio is not None
        and explicit_to_worst_ratio >= min_ratio - tolerance,
        "steps_consistent": steps_consistent,
        "step_size_consistent": step_size_consistent,
        "omega_consistent": omega_consistent,
    }

    return {
        "policy": "geometric_integrator_energy_drift_gate",
        "rows": rows,
        "required_geometric_methods": list(required),
        "explicit_method": explicit,
        "missing_required_methods": missing_required,
        "common_steps": common_steps,
        "common_step_size_s": common_step_size,
        "common_omega_rad_per_s": common_omega,
        "max_geometric_rel_energy_drift": worst_geometric,
        "best_geometric_rel_energy_drift": best_geometric,
        "best_geometric_method": best_geometric_method,
        "explicit_rel_energy_drift": None if explicit_row is None else explicit_row["max_rel_energy_drift"],
        "explicit_to_worst_geometric_drift_ratio": explicit_to_worst_ratio,
        "checks": checks,
        "status": "ok" if all(checks.values()) else "needs_attention",
        "max_allowed_geometric_rel_drift": max_allowed,
        "min_explicit_to_geometric_drift_ratio": min_ratio,
        "tol": tolerance,
    }


def lcurve_corner_choice(alphas, residual_norms, solution_norms, tol=1.0e-12, eps=1.0e-300):
    """Choose an L-curve corner from a recorded Tikhonov path.

    MATLAB/Gypsilab teaching notebooks already expose the path as alpha values,
    residual norms, and solution norms.  This public gate consumes those plain
    arrays directly, computes the same log-log three-point curvature rule, and
    records enough checks for JSON/notebook replay without depending on MATLAB.
    """

    alpha_values = [float(value) for value in alphas]
    residual_values = [float(value) for value in residual_norms]
    solution_values = [float(value) for value in solution_norms]
    tolerance = float(tol)
    floor = float(eps)
    if len(alpha_values) < 3:
        raise ValueError("at least three alpha values are required")
    if not (len(alpha_values) == len(residual_values) == len(solution_values)):
        raise ValueError("alphas, residual_norms, and solution_norms must have the same length")
    if tolerance < 0.0:
        raise ValueError("tol must be non-negative")
    if floor <= 0.0:
        raise ValueError("eps must be > 0")
    if any(value <= 0.0 for value in alpha_values):
        raise ValueError("alphas must be positive for log-log L-curve curvature")
    if any(value < 0.0 for value in residual_values):
        raise ValueError("residual_norms must be non-negative")
    if any(value < 0.0 for value in solution_values):
        raise ValueError("solution_norms must be non-negative")

    curvature = [0.0 for _ in alpha_values]
    log_points = [
        (math.log(max(residual, floor)), math.log(max(solution, floor)))
        for residual, solution in zip(residual_values, solution_values)
    ]
    for index in range(1, len(alpha_values) - 1):
        left = (
            log_points[index][0] - log_points[index - 1][0],
            log_points[index][1] - log_points[index - 1][1],
        )
        right = (
            log_points[index + 1][0] - log_points[index][0],
            log_points[index + 1][1] - log_points[index][1],
        )
        chord = (
            log_points[index + 1][0] - log_points[index - 1][0],
            log_points[index + 1][1] - log_points[index - 1][1],
        )
        denom = (
            math.hypot(*left)
            * math.hypot(*right)
            * math.hypot(*chord)
        )
        if denom > 0.0:
            curvature[index] = 2.0 * abs(left[0] * right[1] - left[1] * right[0]) / denom

    selected_index0 = max(range(len(curvature)), key=lambda index: curvature[index])
    alphas_strictly_increasing = all(
        b > a for a, b in zip(alpha_values, alpha_values[1:])
    )
    residuals_nondecreasing = all(
        b + tolerance >= a for a, b in zip(residual_values, residual_values[1:])
    )
    solution_norms_nonincreasing = all(
        a + tolerance >= b for a, b in zip(solution_values, solution_values[1:])
    )
    checks = {
        "alphas_strictly_increasing": alphas_strictly_increasing,
        "residuals_nondecreasing": residuals_nondecreasing,
        "solution_norms_nonincreasing": solution_norms_nonincreasing,
        "endpoint_curvatures_zero": curvature[0] == 0.0 and curvature[-1] == 0.0,
        "interior_selected": 0 < selected_index0 < len(alpha_values) - 1,
        "positive_corner_curvature": curvature[selected_index0] > 0.0,
    }
    return {
        "policy": "lcurve_corner_regularization_choice",
        "alphas": alpha_values,
        "residual_norms": residual_values,
        "solution_norms": solution_values,
        "curvature": curvature,
        "selected_index": selected_index0 + 1,
        "selected_index0": selected_index0,
        "selected_alpha": alpha_values[selected_index0],
        "selected_residual_norm": residual_values[selected_index0],
        "selected_solution_norm": solution_values[selected_index0],
        "max_curvature": curvature[selected_index0],
        "checks": checks,
        "status": "ok" if all(checks.values()) else "needs_attention",
        "tol": tolerance,
        "eps": floor,
    }


def morozov_discrepancy_choice(alphas, residual_norms, noise_norm, tol=1.0e-12):
    """Choose a regularization weight by the Morozov discrepancy principle.

    This is intentionally solver-independent: the caller supplies a monotone
    residual path, and the helper chooses the row whose residual is closest to
    the estimated noise norm.  It is a teaching gate for readable FEM/BEM
    Tikhonov paths, not a replacement for a full inverse-problem workflow.
    """

    alpha_values = [float(value) for value in alphas]
    residual_values = [float(value) for value in residual_norms]
    noise = float(noise_norm)
    tolerance = float(tol)
    if not alpha_values:
        raise ValueError("alphas must not be empty")
    if len(alpha_values) != len(residual_values):
        raise ValueError("alphas and residual_norms must have the same length")
    if noise <= 0.0:
        raise ValueError("noise_norm must be > 0")
    if tolerance < 0.0:
        raise ValueError("tol must be non-negative")
    if any(value < 0.0 for value in alpha_values):
        raise ValueError("alphas must be non-negative")
    if any(value < 0.0 for value in residual_values):
        raise ValueError("residual_norms must be non-negative")

    alphas_sorted = all(b >= a for a, b in zip(alpha_values, alpha_values[1:]))
    residuals_nondecreasing = all(
        b + tolerance >= a for a, b in zip(residual_values, residual_values[1:])
    )
    errors = [abs(value - noise) for value in residual_values]
    selected_index0 = min(range(len(errors)), key=lambda index: errors[index])
    bracketed = min(residual_values) <= noise <= max(residual_values)
    lower_index0 = None
    upper_index0 = None
    for index, value in enumerate(residual_values):
        if value <= noise:
            lower_index0 = index
        if upper_index0 is None and value >= noise:
            upper_index0 = index
    checks = {
        "alphas_sorted": alphas_sorted,
        "residuals_nondecreasing": residuals_nondecreasing,
        "noise_bracketed": bracketed,
        "selected_minimizes_discrepancy": errors[selected_index0] == min(errors),
    }
    return {
        "policy": "morozov_discrepancy_regularization_choice",
        "alphas": alpha_values,
        "residual_norms": residual_values,
        "noise_norm": noise,
        "selected_index": selected_index0 + 1,
        "selected_index0": selected_index0,
        "selected_alpha": alpha_values[selected_index0],
        "selected_residual_norm": residual_values[selected_index0],
        "selected_abs_error": errors[selected_index0],
        "lower_bracket_index": None if lower_index0 is None else lower_index0 + 1,
        "upper_bracket_index": None if upper_index0 is None else upper_index0 + 1,
        "lower_bracket_alpha": None if lower_index0 is None else alpha_values[lower_index0],
        "upper_bracket_alpha": None if upper_index0 is None else alpha_values[upper_index0],
        "checks": checks,
        "status": "ok" if all(checks.values()) else "needs_attention",
        "tol": tolerance,
    }


def acoustic_plane_wave_intensity_convention_gate(
    pressure_peak,
    rho=1.2,
    c=343.0,
    area=1.0,
    tol=1.0e-12,
):
    """Check peak/RMS acoustic plane-wave power conventions.

    For a normally incident plane wave, ``Z0=rho*c``, ``v_peak=p_peak/Z0``,
    and the time-averaged intensity is ``0.5*p_peak^2/Z0``.  The same value
    must be obtained from RMS pressure as ``p_rms^2/Z0``.  This keeps COMSOL
    acoustic LiveLink slots and open radia-ngsolve notebooks aligned on
    amplitude convention before comparing any solver field.
    """

    p_peak = float(pressure_peak)
    density = float(rho)
    sound_speed = float(c)
    boundary_area = float(area)
    if density <= 0.0:
        raise ValueError("rho must be > 0")
    if sound_speed <= 0.0:
        raise ValueError("c must be > 0")
    if boundary_area <= 0.0:
        raise ValueError("area must be > 0")
    z0 = density * sound_speed
    v_peak = p_peak / z0
    p_rms = p_peak / math.sqrt(2.0)
    v_rms = v_peak / math.sqrt(2.0)
    intensity_from_peak = 0.5 * p_peak * v_peak
    intensity_from_rms = p_rms * v_rms
    power_from_peak = boundary_area * intensity_from_peak
    power_from_rms = boundary_area * intensity_from_rms
    residual = abs(power_from_peak - power_from_rms)
    checks = {
        "specific_impedance_positive": z0 > 0.0,
        "peak_rms_power_match": residual <= float(tol),
        "intensity_nonnegative": intensity_from_peak >= -float(tol),
    }
    return {
        "policy": "acoustic_plane_wave_peak_rms_intensity_gate",
        "pressure_peak_Pa": p_peak,
        "pressure_rms_Pa": p_rms,
        "velocity_peak_m_per_s": v_peak,
        "velocity_rms_m_per_s": v_rms,
        "rho_kg_per_m3": density,
        "c_m_per_s": sound_speed,
        "specific_impedance_Pa_s_per_m": z0,
        "area_m2": boundary_area,
        "intensity_from_peak_W_per_m2": intensity_from_peak,
        "intensity_from_rms_W_per_m2": intensity_from_rms,
        "power_from_peak_W": power_from_peak,
        "power_from_rms_W": power_from_rms,
        "peak_rms_power_residual_W": residual,
        "checks": checks,
        "status": "ok" if all(checks.values()) else "needs_attention",
        "tol": float(tol),
    }


def acoustic_normal_incidence_interface_gate(
    rho_left,
    c_left,
    rho_right,
    c_right,
    pressure_incident_peak=1.0,
    area=1.0,
    tol=1.0e-12,
):
    """Check normal-incidence acoustic pressure/velocity interface identities."""

    rho_l = float(rho_left)
    rho_r = float(rho_right)
    c_l = float(c_left)
    c_r = float(c_right)
    p_inc = float(pressure_incident_peak)
    surface_area = float(area)
    tolerance = float(tol)
    if rho_l <= 0.0 or rho_r <= 0.0:
        raise ValueError("densities must be > 0")
    if c_l <= 0.0 or c_r <= 0.0:
        raise ValueError("sound speeds must be > 0")
    if surface_area <= 0.0:
        raise ValueError("area must be > 0")
    if tolerance < 0.0:
        raise ValueError("tol must be non-negative")

    z_left = rho_l * c_l
    z_right = rho_r * c_r
    reflection = (z_right - z_left) / (z_right + z_left)
    transmission_pressure = 2.0 * z_right / (z_right + z_left)
    p_reflected = reflection * p_inc
    p_transmitted = transmission_pressure * p_inc
    v_left_interface = (p_inc - p_reflected) / z_left
    v_right_interface = p_transmitted / z_right
    incident_intensity = 0.5 * p_inc * p_inc / z_left
    reflected_intensity = reflection * reflection * incident_intensity
    transmitted_intensity = 0.5 * p_transmitted * p_transmitted / z_right
    reflected_power_ratio = reflected_intensity / incident_intensity
    transmitted_power_ratio = transmitted_intensity / incident_intensity
    power_balance_error = reflected_power_ratio + transmitted_power_ratio - 1.0
    pressure_jump = p_inc + p_reflected - p_transmitted
    velocity_jump = v_left_interface - v_right_interface
    checks = {
        "pressure_continuity": abs(pressure_jump) <= tolerance * max(1.0, abs(p_inc), abs(p_transmitted)),
        "velocity_continuity": abs(velocity_jump) <= tolerance * max(1.0, abs(v_left_interface), abs(v_right_interface)),
        "energy_flux_balance": abs(power_balance_error) <= tolerance,
        "passive_power_split": reflected_power_ratio >= -tolerance and transmitted_power_ratio >= -tolerance,
    }
    return {
        "policy": "acoustic_normal_incidence_interface_gate",
        "rho_left_kg_per_m3": rho_l,
        "c_left_m_per_s": c_l,
        "rho_right_kg_per_m3": rho_r,
        "c_right_m_per_s": c_r,
        "z_left_Pa_s_per_m": z_left,
        "z_right_Pa_s_per_m": z_right,
        "pressure_incident_peak_Pa": p_inc,
        "pressure_reflection_coefficient": reflection,
        "pressure_transmission_coefficient": transmission_pressure,
        "pressure_reflected_peak_Pa": p_reflected,
        "pressure_transmitted_peak_Pa": p_transmitted,
        "velocity_left_interface_peak_m_per_s": v_left_interface,
        "velocity_right_interface_peak_m_per_s": v_right_interface,
        "incident_intensity_W_per_m2": incident_intensity,
        "reflected_intensity_W_per_m2": reflected_intensity,
        "transmitted_intensity_W_per_m2": transmitted_intensity,
        "incident_power_W": surface_area * incident_intensity,
        "reflected_power_W": surface_area * reflected_intensity,
        "transmitted_power_W": surface_area * transmitted_intensity,
        "reflected_power_ratio": reflected_power_ratio,
        "transmitted_power_ratio": transmitted_power_ratio,
        "power_balance_error": power_balance_error,
        "pressure_jump_Pa": pressure_jump,
        "velocity_jump_m_per_s": velocity_jump,
        "checks": checks,
        "status": "ok" if all(checks.values()) else "needs_attention",
        "tol": tolerance,
    }
