# -*- coding: utf-8 -*-
"""Small solver-independent gates used by CAE loop slots.

These helpers intentionally avoid any commercial-tool provenance.  They encode
the physics checks that a loop slot can reuse before the source-tool/private
lane records where the reference data came from.
"""
from __future__ import annotations

import cmath
import math

from radia_mcp.radia_ngsolve.solve import (
    carter_coefficient,
    effective_air_gap,
    slotted_air_gap_permeance_factor,
)


MU0 = 4.0e-7 * math.pi
EPS0 = 8.8541878128e-12


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


def coenergy_torque_periodic_summary(theta_rad, coenergy_j, torque_nm, rtol=1.0e-6, atol=1.0e-9):
    """Check torque samples against ``T = dWprime/dtheta``.

    The samples must be equally spaced over one periodic cycle.  Near torque
    zero crossings the absolute tolerance is decisive; away from zero crossings
    the relative tolerance is reported as well.
    """

    theta = [float(value) for value in theta_rad]
    w = [float(value) for value in coenergy_j]
    torque = [float(value) for value in torque_nm]
    if not (len(theta) == len(w) == len(torque)):
        raise ValueError("theta_rad, coenergy_j, and torque_nm must have the same length")
    if len(theta) < 3:
        raise ValueError("at least three samples are required")
    steps = [theta[i + 1] - theta[i] for i in range(len(theta) - 1)]
    h = sum(steps) / len(steps)
    if any(abs(step - h) > max(1.0e-12, 1.0e-9 * abs(h)) for step in steps):
        raise ValueError("theta samples must be equally spaced")
    estimated = central_difference_periodic(w, h)
    rows = []
    for angle, ref, est in zip(theta, torque, estimated):
        abs_error = abs(est - ref)
        rel_error = abs_error / max(abs(est), abs(ref), 1.0e-300)
        passed = abs_error <= float(atol) or rel_error <= float(rtol)
        rows.append({
            "theta_rad": angle,
            "reference_torque_nm": ref,
            "estimated_torque_nm": est,
            "abs_error": abs_error,
            "rel_error": rel_error,
            "passed": passed,
        })
    return {
        "policy": "coenergy_torque_periodic_derivative_gate",
        "n_samples": len(rows),
        "rtol": float(rtol),
        "atol": float(atol),
        "max_abs_error": max(row["abs_error"] for row in rows),
        "max_rel_error": max(row["rel_error"] for row in rows),
        "n_passed": sum(1 for row in rows if row["passed"]),
        "status": "ok" if all(row["passed"] for row in rows) else "needs_attention",
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
        "checks": checks,
        "tol": tolerance,
        "notes": [
            "Run this before converting RI/MA/DB numeric pairs into complex S-parameters.",
            "A passive row can still be unusable if ports are swapped, z0 is missing, or the frequency unit is ambiguous.",
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
    """Bundle far-field metadata and one lobe row before notebook use.

    A notebook panel should not receive a naked gain scalar.  The lobe row must
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
            "Run this after the far-field export metadata gate and before a notebook panel plots or ranks lobes.",
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


def shared_solver_session_health_gate(
    connected,
    api_visible,
    discovered_engines=None,
    status="",
    started_new_process=False,
    killed_process=False,
):
    """Check solver-session health without treating it as physics validation.

    This public-safe gate records whether an external solver session was reused
    cleanly before a numerical validation row is trusted.  It is deliberately
    generic: COMSOL LiveLink, MATLAB Engine, Jupyter kernels, and similar
    long-lived solver sessions can all use the same separation between session
    health and physics residuals.
    """

    engines = [] if discovered_engines is None else [str(item) for item in discovered_engines]
    state = str(status or "").strip().lower()
    checks = {
        "session_connected": bool(connected),
        "api_visible": bool(api_visible),
        "engine_discovered": len(engines) > 0,
        "status_allows_reuse": state in {"already-connected", "connected", "live", "ok"},
        "started_no_new_process": not bool(started_new_process),
        "killed_no_process": not bool(killed_process),
    }
    return {
        "policy": "shared_solver_session_health_gate",
        "connected": bool(connected),
        "api_visible": bool(api_visible),
        "discovered_engines": engines,
        "status": state,
        "started_new_process": bool(started_new_process),
        "killed_process": bool(killed_process),
        "checks": checks,
        "status_label": "ok" if all(checks.values()) else "needs_attention",
        "notes": [
            "session health is a preflight, not a physical residual",
            "already-connected is valid reuse evidence when no new process was started or killed",
        ],
    }


def solver_result_table_metadata_gate(
    metadata,
    required_columns=(),
    required_units=None,
    independent_axis=None,
    expected_source=None,
    min_rows=1,
):
    """Check solver result-table metadata before interpreting numeric rows.

    This is a solver-independent companion to COMSOL ``mphglobal``/``mphtable``,
    MATLAB tables, CSV exports, and measurement logs.  It verifies that columns,
    units, an independent sweep axis, row count, and optional source label are
    explicit before a table is used as validation evidence.
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
        "checks": checks,
        "notes": [
            "Run this before numeric residuals so column position does not become hidden solver knowledge.",
            "A plausible result row can still be unusable if units, sweep axis, or source table identity were lost.",
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
    """Check block-label rows before a FEMM-style motor model handoff.

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


def femm_pm_magnetization_convention_gate(
    pm_rows,
    required_regions=None,
    allowed_frames=("global_xy", "rotor_xy", "local_radial", "local_tangential"),
):
    """Check FEMM-style permanent-magnet direction rows before emission.

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

    The expected step order is the ELF/MAGIC-style irreversible-demag workflow:
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

    This is a solver-independent preflight for ELF/MAGIC-style demagnetization
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
    require_closed_charge_balance=True,
    tol=1.0e-12,
):
    """Check PM BEM surface-normal metadata before demag/source assembly.

    Magnetic-charge BEM uses ``sigma_m = M dot n``.  A surface row is not
    solver-ready unless it records outward-normal convention, positive area,
    unit normal, and magnetization direction.  For a closed PM with uniform
    magnetization, the signed surface-charge proxy must sum to zero.
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
    charge_sum = 0.0
    total_area = 0.0

    def vector3(value):
        if value is None or len(value) != 3:
            return None
        return [float(value[0]), float(value[1]), float(value[2])]

    for index, row in enumerate(rows_in, start=1):
        if not isinstance(row, dict):
            raise ValueError("each surface row must be a dictionary")
        name = str(row.get("surface") or row.get("name") or row.get("label") or "").strip()
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
        "closed_surface_charge_balances": (abs(charge_sum) <= tolerance * balance_scale)
        if require_closed_charge_balance
        else True,
    }
    return {
        "policy": "pm_bem_surface_normal_metadata_gate",
        "n_rows": len(normalized_rows),
        "expected_normal_convention": expected,
        "require_closed_charge_balance": bool(require_closed_charge_balance),
        "duplicate_surfaces": duplicate_surfaces,
        "missing_or_wrong_convention_surfaces": sorted(missing_convention),
        "nonpositive_area_surfaces": sorted(nonpositive_area),
        "bad_normal_surfaces": sorted(bad_normals),
        "missing_magnetization_surfaces": sorted(missing_magnetization),
        "signed_charge_proxy_sum": charge_sum,
        "signed_charge_proxy_abs": abs(charge_sum),
        "total_area_m2": total_area,
        "tol": tolerance,
        "rows": normalized_rows,
        "checks": checks,
        "status": "ok" if all(checks.values()) else "needs_attention",
        "version_note": (
            "Use this before PM magnetic-charge BEM assembly so surface normals, "
            "magnetization direction, and outward convention are explicit."
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
