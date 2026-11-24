"""Compare sampled no-current cogging torque on a common mechanical grid."""

from math import lcm

import numpy as np


_METHODS = {"maxwell_stress", "virtual_work", "arkkio"}


def _read_lane(lane):
    if not isinstance(lane, dict):
        raise ValueError("each lane must be an object")
    slots, poles = lane["slots"], lane["poles"]
    if any(type(count) is not int or count < 1 for count in (slots, poles)):
        raise ValueError("slots and poles must be positive integers")
    symmetry = lane["slot_pole_symmetry"]
    drop_endpoint = lane.get("drop_duplicate_endpoint", False)
    if type(symmetry) is not bool or type(drop_endpoint) is not bool:
        raise ValueError("symmetry and endpoint declarations must be boolean")
    method, skew = lane["extraction_method"], lane["skew"]
    if not isinstance(method, str) or method not in _METHODS:
        raise ValueError("unsupported extraction_method")
    if not isinstance(skew, str) or not skew.strip():
        raise ValueError("declare skew as 'none' or describe its geometry and averaging")
    skew = skew.strip()
    angle = np.asarray(lane["mechanical_angles_deg"], dtype=float)
    torque = np.asarray(lane["torque_nm"], dtype=float)
    if angle.ndim != 1 or angle.shape != torque.shape or angle.size < 3:
        raise ValueError("provide equal one-dimensional angle and torque arrays")
    if not np.all(np.isfinite(angle)) or not np.all(np.isfinite(torque)):
        raise ValueError("angle and torque must be finite")
    return slots, poles, symmetry, method, skew, angle, torque, drop_endpoint


def _periodic_samples(angle, torque, drop_endpoint, period):
    # The absolute tolerance accommodates rounding without resampling either lane.
    grid_atol = 1e-7
    step = np.diff(angle)
    if np.any(step <= 0) or not np.allclose(step, step[0], atol=grid_atol, rtol=0):
        raise ValueError("mechanical angles must be increasing and uniformly spaced")
    if abs(angle[-1] - angle[0] - period) <= grid_atol:
        if not drop_endpoint:
            raise ValueError("duplicate periodic endpoint; declare drop_duplicate_endpoint to discard it")
        angle, torque = angle[:-1], torque[:-1]
    if angle.size < 3 or abs(angle.size * step[0] - period) > grid_atol:
        raise ValueError("grid must cover one complete comparison period")
    return angle, torque


def cogging_torque_cross_lane_gate(first, second, tolerances):
    """Compare signed waveforms, sampled peak-to-peak and resolved slot/pole orders."""
    try:
        a, b = _read_lane(first), _read_lane(second)
        if a[:2] != b[:2]:
            raise ValueError("slot and pole counts must match")
        if not isinstance(tolerances, dict):
            raise ValueError("tolerances must be an object")
        limits = {}
        for name in ("peak_to_peak_nm", "waveform_rms_nm", "harmonic_nm"):
            value = tolerances[name]
            if type(value) not in (float, int) or not np.isfinite(float(value)) or value < 0:
                raise ValueError(f"{name} tolerance must be finite and nonnegative")
            limits[name] = float(value)

        slot_pole_order = lcm(*a[:2])
        base_order = slot_pole_order if a[2] and b[2] else 1
        period = 360.0 / base_order
        angle_a, torque_a = _periodic_samples(*a[5:], period)
        angle_b, torque_b = _periodic_samples(*b[5:], period)
        if angle_a.shape != angle_b.shape or not np.allclose(angle_a, angle_b, atol=1e-7, rtol=0):
            raise ValueError("lanes must use a common mechanical angle grid and angular datum")

        # Exclude Nyquist, where a signed phase cannot be resolved. FFT phase is
        # referenced to the common first sample, so complex differences keep sign.
        bins = np.arange(1, (angle_a.size + 1) // 2)
        orders = bins * base_order
        selected = orders % slot_pole_order == 0
        if not np.any(selected):
            raise ValueError("grid cannot resolve the first slot/pole harmonic below Nyquist")
        bins, orders = bins[selected], orders[selected]
        with np.errstate(over="raise", invalid="raise"):
            first_fft = np.fft.rfft(torque_a / angle_a.size)[bins] * 2
            second_fft = np.fft.rfft(torque_b / angle_a.size)[bins] * 2
            errors = np.abs(first_fft - second_fft)
            peak_a, peak_b = float(np.ptp(torque_a)), float(np.ptp(torque_b))
            waveform_rms = float(np.sqrt(np.mean(np.square(torque_a - torque_b))))
            if not np.all(np.isfinite(errors)) or not np.all(np.isfinite(np.abs(first_fft))) or not np.all(np.isfinite(np.abs(second_fft))):
                raise ValueError("harmonic calculation exceeds numeric range")
    except (KeyError, ValueError, TypeError, OverflowError, FloatingPointError) as exc:
        return {"status": "invalid_input", "checks": {"inputs": False}, "error": str(exc)}

    checks = {
        "same_skew": a[4] == b[4],
        "peak_to_peak": abs(peak_a - peak_b) <= limits["peak_to_peak_nm"],
        "signed_waveform_rms": waveform_rms <= limits["waveform_rms_nm"],
        "slot_pole_harmonics": bool(np.all(errors <= limits["harmonic_nm"])),
    }
    findings = []
    if a[3] != b[3]:
        findings.append("extraction_methods_differ")
    if not checks["same_skew"]:
        findings.append("skew_declarations_differ")
    return {
        "status": "pass" if all(checks.values()) else "mismatch",
        "checks": checks,
        "findings": findings,
        "extraction_methods": [a[3], b[3]],
        "skew": [a[4], b[4]],
        "slot_pole_symmetry_applied": bool(a[2] and b[2]),
        "slot_pole_order": slot_pole_order,
        "comparison_period_deg_mechanical": period,
        "angle_origin_deg_mechanical": float(angle_a[0]),
        "peak_to_peak_nm": [peak_a, peak_b],
        "signed_waveform_rms_difference_nm": waveform_rms,
        "harmonics": [
            {
                "order_per_mechanical_revolution": int(order),
                "amplitudes_nm": [float(abs(left)), float(abs(right))],
                "complex_difference_nm": float(error),
            }
            for order, left, right, error in zip(orders, first_fft, second_fft, errors)
        ],
    }
