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
