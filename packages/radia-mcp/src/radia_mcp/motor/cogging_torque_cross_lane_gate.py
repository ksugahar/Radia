"""Compare sampled no-current cogging torque on a common mechanical grid."""

from dataclasses import dataclass
from math import isfinite, lcm

import numpy as np


_METHODS = {"maxwell_stress", "virtual_work", "arkkio"}


@dataclass(frozen=True)
class _Lane:
    slots: int
    pole_count: int
    symmetry: bool
    method: str
    skew: str
    angle: np.ndarray
    torque: np.ndarray
    drop_endpoint: bool
    torque_basis: str
    stack_length: float
    model_dimension: str
    model_fraction: float

    @property
    def torque_factor(self):
        length_factor = self.stack_length if self.torque_basis == "per_metre" else 1.0
        return length_factor / self.model_fraction


def _positive_number(value, name):
    if type(value) not in (int, float) or not isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be finite and positive")
    return float(value)


def _read_lane(lane):
    if not isinstance(lane, dict):
        raise ValueError("each lane must be an object")
    slots, pole_count = lane["slots"], lane["pole_count"]
    if any(type(count) is not int or count < 1 for count in (slots, pole_count)):
        raise ValueError("slots and pole_count must be positive integers")
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
    if skew.casefold() == "none":
        skew = "none"
    basis, dimension = lane["torque_basis"], lane["model_dimension"]
    if basis not in ("per_metre", "total"):
        raise ValueError("torque_basis must be per_metre or total")
    if dimension not in ("2d", "3d"):
        raise ValueError("model_dimension must be 2d or 3d")
    length = _positive_number(lane["stack_length_m"], "stack_length_m")
    fraction = _positive_number(lane["model_fraction"], "model_fraction")
    if fraction > 1:
        raise ValueError("model_fraction must be at most 1")
    if "excitation" in lane:
        excitation = lane["excitation"]
        if (
            not isinstance(excitation, dict)
            or type(excitation.get("phase_current_A")) not in (int, float)
            or excitation["phase_current_A"] != 0
        ):
            raise ValueError("no-current cogging requires excitation.phase_current_A = 0")
    angle = np.asarray(lane["mechanical_angles_deg"], dtype=float)
    torque = np.asarray(lane["torque_nm"], dtype=float)
    if angle.ndim != 1 or angle.shape != torque.shape or angle.size < 3:
        raise ValueError("provide equal one-dimensional angle and torque arrays")
    if not np.all(np.isfinite(angle)) or not np.all(np.isfinite(torque)):
        raise ValueError("angle and torque must be finite")
    return _Lane(slots, pole_count, symmetry, method, skew, angle, torque, drop_endpoint, basis, length, dimension, fraction)


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


def _stack_scaling_suspected(first, second, lengths, tolerance):
    # A candidate scale must align the signed waveform, including its phase.
    for length in lengths:
        factor = min(length, 1 / length)
        for left, right in ((first, second), (second, first)):
            if np.max(np.abs(left - right * factor)) <= tolerance:
                return True
    return False


def cogging_torque_cross_lane_gate(first, second, tolerances):
    """Compare signed waveforms, sampled peak-to-peak and resolved slot/pole orders."""
    findings = []
    try:
        a, b = _read_lane(first), _read_lane(second)
        if a.pole_count != b.pole_count:
            if max(a.pole_count, b.pole_count) == 2 * min(a.pole_count, b.pole_count):
                findings.append("pole_count_vs_pole_pairs_ambiguity")
            raise ValueError("pole_count must match")
        if a.slots != b.slots:
            raise ValueError("slot counts must match")
        if not isinstance(tolerances, dict):
            raise ValueError("tolerances must be an object")
        limits = {}
        for name in ("peak_to_peak_nm", "waveform_rms_nm", "harmonic_nm"):
            value = tolerances[name]
            if type(value) not in (float, int) or not np.isfinite(float(value)) or value < 0:
                raise ValueError(f"{name} tolerance must be finite and nonnegative")
            limits[name] = float(value)

        slot_pole_order = lcm(a.slots, a.pole_count)
        base_order = slot_pole_order if a.symmetry and b.symmetry else 1
        period = 360.0 / base_order
        angle_a, torque_a = _periodic_samples(a.angle, a.torque, a.drop_endpoint, period)
        angle_b, torque_b = _periodic_samples(b.angle, b.torque, b.drop_endpoint, period)
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
            factors = [a.torque_factor, b.torque_factor]
            torque_a = torque_a * factors[0]
            torque_b = torque_b * factors[1]
            means = [float(np.mean(torque_a)), float(np.mean(torque_b))]
            first_fft = np.fft.rfft(torque_a / angle_a.size)[bins] * 2
            second_fft = np.fft.rfft(torque_b / angle_a.size)[bins] * 2
            errors = np.abs(first_fft - second_fft)
            peak_a, peak_b = float(np.ptp(torque_a)), float(np.ptp(torque_b))
            waveform_rms = float(np.sqrt(np.mean(np.square(torque_a - torque_b))))
            if not np.all(np.isfinite(errors)) or not np.all(np.isfinite(np.abs(first_fft))) or not np.all(np.isfinite(np.abs(second_fft))):
                raise ValueError("harmonic calculation exceeds numeric range")
    except (KeyError, ValueError, TypeError, OverflowError, FloatingPointError) as exc:
        return {"status": "invalid_input", "checks": {"inputs": False}, "error": str(exc), "findings": findings}

    checks = {
        "same_skew": a.skew == b.skew,
        "zero_mean_torque": all(abs(mean) <= limits["harmonic_nm"] for mean in means),
        "peak_to_peak": abs(peak_a - peak_b) <= limits["peak_to_peak_nm"],
        "signed_waveform_rms": waveform_rms <= limits["waveform_rms_nm"],
        "slot_pole_harmonics": bool(np.all(errors <= limits["harmonic_nm"])),
    }
    if a.method != b.method:
        findings.append("extraction_methods_differ")
    if not checks["same_skew"]:
        findings.append("skew_declarations_differ")
    if not checks["zero_mean_torque"]:
        findings.append("nonzero_mean_torque")
    if a.model_dimension != b.model_dimension:
        findings.append("model_difference_declared_required")
    if not checks["signed_waveform_rms"] and _stack_scaling_suspected(
        torque_a, torque_b, (a.stack_length, b.stack_length), limits["waveform_rms_nm"]
    ):
        findings.append("stack_scaling_suspected")
    return {
        "status": "pass" if all(checks.values()) else "mismatch",
        "checks": checks,
        "findings": findings,
        "extraction_methods": [a.method, b.method],
        "skew": [a.skew, b.skew],
        "slot_pole_symmetry_applied": a.symmetry and b.symmetry,
        "slot_pole_order": slot_pole_order,
        "comparison_period_deg_mechanical": period,
        "angle_origin_deg_mechanical": float(angle_a[0]),
        "model_dimensions": [a.model_dimension, b.model_dimension],
        "torque_to_full_machine_factors": factors,
        "mean_torque_nm": means,
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
