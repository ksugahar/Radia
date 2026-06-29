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
