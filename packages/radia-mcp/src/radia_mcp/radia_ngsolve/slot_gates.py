# -*- coding: utf-8 -*-
"""Small solver-independent gates used by CAE loop slots.

These helpers intentionally avoid any commercial-tool provenance.  They encode
the physics checks that a loop slot can reuse before the source-tool/private
lane records where the reference data came from.
"""
from __future__ import annotations

import cmath
import math


MU0 = 4.0e-7 * math.pi


def parallel_wire_force_per_length(i1, i2, separation_m, mu0=MU0):
    """Return the signed force per length between two long parallel wires.

    Positive means attraction for equal current directions.  The identity is
    ``F/L = mu0 * I1 * I2 / (2*pi*d)``.
    """

    d = float(separation_m)
    if d <= 0.0:
        raise ValueError("separation_m must be > 0")
    return float(mu0) * float(i1) * float(i2) / (2.0 * math.pi * d)


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


def _complex(value):
    if isinstance(value, dict):
        return complex(float(value.get("real", 0.0)), float(value.get("imag", 0.0)))
    return complex(value)


def _complex_row(value):
    return {"real": float(value.real), "imag": float(value.imag), "abs": abs(value)}


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


def spwm_snapshot_current_handoff_summary(
    id_current,
    iq_current,
    sample_count,
    theta0_rad=0.0,
    sample_offset_fraction=0.5,
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
    return {
        "policy": "spwm_snapshot_current_handoff_gate",
        "id": d,
        "iq": q,
        "current_amplitude": amplitude,
        "sample_count": n,
        "theta0_rad": theta0,
        "sample_offset_fraction": offset,
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
