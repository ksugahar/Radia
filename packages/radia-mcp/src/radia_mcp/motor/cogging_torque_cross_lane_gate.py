"""Compare sampled no-current cogging torque on a common mechanical grid."""

from dataclasses import dataclass
from cmath import rect
from math import fsum, isclose, isfinite, lcm, sqrt, tau


_METHODS = {"maxwell_stress", "virtual_work", "arkkio"}


@dataclass(frozen=True)
class _Lane:
    slots: int
    pole_count: int
    symmetry: bool
    method: str
    skew: str
    angle: list[float]
    torque: list[float]
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
    angles, samples = lane["mechanical_angles_deg"], lane["torque_nm"]
    if not isinstance(angles, (list, tuple)) or not isinstance(samples, (list, tuple)):
        raise ValueError("provide one-dimensional angle and torque arrays")
    if len(angles) != len(samples) or len(angles) < 3:
        raise ValueError("provide equal one-dimensional angle and torque arrays")
    angle = [float(value) for value in angles]
    torque = [float(value) for value in samples]
    if not all(map(isfinite, angle)) or not all(map(isfinite, torque)):
        raise ValueError("angle and torque must be finite")
    return _Lane(slots, pole_count, symmetry, method, skew, angle, torque, drop_endpoint, basis, length, dimension, fraction)


def _periodic_samples(angle, torque, drop_endpoint, period):
    # The absolute tolerance accommodates rounding without resampling either lane.
    grid_atol = 1e-7
    steps = [right - left for left, right in zip(angle, angle[1:])]
    step = steps[0]
    if any(delta <= 0 or not isclose(delta, step, abs_tol=grid_atol, rel_tol=0) for delta in steps):
        raise ValueError("mechanical angles must be increasing and uniformly spaced")
    if abs(angle[-1] - angle[0] - period) <= grid_atol:
        if not drop_endpoint:
            raise ValueError("duplicate periodic endpoint; declare drop_duplicate_endpoint to discard it")
        angle, torque = angle[:-1], torque[:-1]
    if len(angle) < 3 or abs(len(angle) * step - period) > grid_atol:
        raise ValueError("grid must cover one complete comparison period")
    return angle, torque


def _stack_scaling_suspected(first, second, lengths, tolerance):
    # A candidate scale must align the signed waveform, including its phase.
    for length in lengths:
        factor = min(length, 1 / length)
        for left, right in ((first, second), (second, first)):
            if all(abs(a - b * factor) <= tolerance for a, b in zip(left, right)):
                return True
    return False


def _harmonic_pair(first, second, phase):
    count = len(first)
    weights = [rect(2 / count, phase * index) for index in range(count)]
    return tuple(
        complex(
            fsum(value * weight.real for value, weight in zip(torque, weights)),
            fsum(value * weight.imag for value, weight in zip(torque, weights)),
        )
        for torque in (first, second)
    )


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
            if type(value) not in (float, int) or not isfinite(value) or value < 0:
                raise ValueError(f"{name} tolerance must be finite and nonnegative")
            limits[name] = float(value)

        slot_pole_order = lcm(a.slots, a.pole_count)
        base_order = slot_pole_order if a.symmetry and b.symmetry else 1
        period = 360.0 / base_order
        angle_a, torque_a = _periodic_samples(a.angle, a.torque, a.drop_endpoint, period)
        angle_b, torque_b = _periodic_samples(b.angle, b.torque, b.drop_endpoint, period)
        if len(angle_a) != len(angle_b) or any(
            not isclose(left, right, abs_tol=1e-7, rel_tol=0)
            for left, right in zip(angle_a, angle_b)
        ):
            raise ValueError("lanes must use a common mechanical angle grid and angular datum")

        # Skip Nyquist, whose phase is unresolved. Only project the bins
        # selected by the slot/pole order, using the first sample as phase zero.
        count = len(angle_a)
        bin_stride = slot_pole_order // base_order
        bins = range(bin_stride, (count + 1) // 2, bin_stride)
        if not bins:
            raise ValueError("grid cannot resolve the first slot/pole harmonic below Nyquist")
        factors = [a.torque_factor, b.torque_factor]
        torque_a = [value * factors[0] for value in torque_a]
        torque_b = [value * factors[1] for value in torque_b]
        means = [fsum(value / count for value in torque) for torque in (torque_a, torque_b)]
        peaks = [max(torque) - min(torque) for torque in (torque_a, torque_b)]
        waveform_rms = sqrt(fsum(((left - right) / sqrt(count)) ** 2 for left, right in zip(torque_a, torque_b)))
        harmonics = []
        for bin_index in bins:
            phase = -tau * bin_index / count
            left, right = _harmonic_pair(torque_a, torque_b, phase)
            harmonics.append({
                "order_per_mechanical_revolution": bin_index * base_order,
                "amplitudes_nm": [abs(left), abs(right)],
                "complex_difference_nm": abs(left - right),
            })
        metrics = [waveform_rms, *factors, *means, *peaks]
        for harmonic in harmonics:
            metrics.extend((*harmonic["amplitudes_nm"], harmonic["complex_difference_nm"]))
        if not all(map(isfinite, metrics)):
            raise ValueError("torque calculation exceeds numeric range")
    except (KeyError, ValueError, TypeError, OverflowError, FloatingPointError) as exc:
        return {"status": "invalid_input", "checks": {"inputs": False}, "error": str(exc), "findings": findings}

    checks = {
        "same_skew": a.skew == b.skew,
        "zero_mean_torque": all(abs(mean) <= limits["harmonic_nm"] for mean in means),
        "peak_to_peak": abs(peaks[0] - peaks[1]) <= limits["peak_to_peak_nm"],
        "signed_waveform_rms": waveform_rms <= limits["waveform_rms_nm"],
        "slot_pole_harmonics": all(h["complex_difference_nm"] <= limits["harmonic_nm"] for h in harmonics),
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
        "peak_to_peak_nm": peaks,
        "signed_waveform_rms_difference_nm": waveform_rms,
        "harmonics": harmonics,
    }
