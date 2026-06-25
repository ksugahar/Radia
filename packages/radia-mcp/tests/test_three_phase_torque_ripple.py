r"""Three-phase back-EMF harmonic mixing -> 6k torque ripple.

Pure power-balance reference: balanced sinusoidal phase currents multiply the
phase back-EMF waveform, and the three-phase sum leaves only n=6k+/-1 harmonic
pairs as 6k instantaneous-power ripples.  At constant speed the torque ripple
has the same normalized amplitudes.
"""
import math
import os
import sys

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.radia_ngsolve.solve import (
    three_phase_torque_ripple_harmonics,
    three_phase_torque_ripple_pair_table,
    torque_angle_sweep_comparison_summary,
    torque_angle_sweep_health_summary,
    torque_angle_sweep_summary,
)


def _sample_power(emf_harmonics, current_peak, samples=4096):
    vals = []
    shifts = (0.0, -2.0 * math.pi / 3.0, 2.0 * math.pi / 3.0)
    for i in range(samples):
        th = 2.0 * math.pi * i / samples
        p = 0.0
        for s in shifts:
            e = sum(complex(E * complex(math.cos(n * (th + s)), math.sin(n * (th + s)))).real
                    for n, E in emf_harmonics.items())
            cur = current_peak * math.cos(th + s)
            p += e * cur
        vals.append(p)
    return vals


def _fourier_amplitude(vals, order):
    n = len(vals)
    c = sum(vals[i] * complex(math.cos(-2.0 * math.pi * order * i / n),
                              math.sin(-2.0 * math.pi * order * i / n))
            for i in range(n)) / n
    return 2.0 * abs(c)


def test_three_phase_harmonic_pairs_map_to_6k_ripple():
    harmonics = {1: 1.0, 5: 0.10, 7: 0.07, 11: 0.02, 13: 0.03, 3: 0.50}
    out = three_phase_torque_ripple_harmonics(harmonics, current_peak=2.0, mechanical_speed=5.0)

    assert math.isclose(out["mean_power"], 3.0, abs_tol=1e-14)
    assert math.isclose(out["power_ripple"][6], 0.51, abs_tol=1e-14)
    assert math.isclose(out["power_ripple"][12], 0.15, abs_tol=1e-14)
    assert math.isclose(out["normalized_ripple"][6], 0.17, abs_tol=1e-14)
    assert math.isclose(out["normalized_ripple"][12], 0.05, abs_tol=1e-14)
    assert 3 not in out["power_ripple"]                    # triplen harmonics cancel in 3-phase power
    assert math.isclose(out["torque_ripple"][6], 0.51 / 5.0, abs_tol=1e-14)


def test_pair_table_exposes_harmonic_budget():
    harmonics = {1: 1.0, 3: 0.5, 5: 0.10, 7: 0.07, 11: 0.02, 13: 0.03}
    out = three_phase_torque_ripple_harmonics(harmonics, current_peak=2.0, mechanical_speed=5.0)
    table = three_phase_torque_ripple_pair_table(harmonics, current_peak=2.0, mechanical_speed=5.0)
    by_order = {row["ripple_order"]: row for row in table}

    assert sorted(by_order) == [6, 12]
    assert by_order[6]["contributing_harmonics"] == [5, 7]
    assert by_order[12]["contributing_harmonics"] == [11, 13]
    assert math.isclose(by_order[6]["emf_phasor_abs"], 0.17, abs_tol=1e-14)
    assert math.isclose(by_order[6]["power_ripple"], out["power_ripple"][6], abs_tol=1e-14)
    assert math.isclose(by_order[6]["torque_ripple"], out["torque_ripple"][6], abs_tol=1e-14)
    assert math.isclose(by_order[12]["normalized_ripple"], out["normalized_ripple"][12], abs_tol=1e-14)


def test_pair_table_keeps_phasor_cancellation_visible():
    harmonics = {1: 1.0, 5: 0.10, 7: -0.10}
    table = three_phase_torque_ripple_pair_table(harmonics)
    assert table[0]["ripple_order"] == 6
    assert table[0]["contributing_harmonics"] == [5, 7]
    assert math.isclose(table[0]["emf_phasor_abs"], 0.0, abs_tol=1e-14)
    assert math.isclose(table[0]["normalized_ripple"], 0.0, abs_tol=1e-14)


def test_closed_form_matches_direct_time_waveform_fourier():
    harmonics = {1: 0.8, 5: 0.12, 7: 0.04, 11: -0.03, 13: 0.01}
    current_peak = 1.7
    out = three_phase_torque_ripple_harmonics(harmonics, current_peak=current_peak)
    p = _sample_power(harmonics, current_peak)
    mean = sum(p) / len(p)

    assert math.isclose(mean, out["mean_power"], abs_tol=1e-12)
    assert math.isclose(_fourier_amplitude(p, 6), out["power_ripple"][6], abs_tol=1e-12)
    assert math.isclose(_fourier_amplitude(p, 12), out["power_ripple"][12], abs_tol=1e-12)


def test_torque_angle_sweep_summary_extracts_ripple_harmonics():
    samples = 720
    mean = 10.0
    ripple6 = 0.4
    ripple12 = 0.1
    torque = [
        mean
        + ripple6 * math.cos(6 * 2.0 * math.pi * idx / samples)
        + ripple12 * math.sin(12 * 2.0 * math.pi * idx / samples)
        for idx in range(samples)
    ]

    summary = torque_angle_sweep_summary(torque, max_harmonic=18)
    by_order = {row["order"]: row for row in summary["harmonic_rows"]}

    assert math.isclose(summary["mean_torque_Nm"], mean, abs_tol=1.0e-14)
    assert math.isclose(summary["ac_rms_torque_Nm"], math.sqrt((ripple6 * ripple6 + ripple12 * ripple12) / 2.0), rel_tol=1.0e-12)
    assert summary["dominant_harmonic"] == 6
    assert math.isclose(by_order[6]["cos_coefficient_Nm"], ripple6, abs_tol=1.0e-14)
    assert math.isclose(by_order[6]["sin_coefficient_Nm"], 0.0, abs_tol=1.0e-14)
    assert math.isclose(by_order[6]["amplitude_Nm"], ripple6, abs_tol=1.0e-14)
    assert math.isclose(by_order[12]["cos_coefficient_Nm"], 0.0, abs_tol=1.0e-14)
    assert math.isclose(by_order[12]["sin_coefficient_Nm"], ripple12, abs_tol=1.0e-14)
    assert math.isclose(by_order[12]["amplitude_Nm"], ripple12, abs_tol=1.0e-14)


def test_torque_angle_sweep_health_summary_flags_ripple_limits():
    samples = 720
    mean = 10.0
    ripple6 = 0.4
    ripple12 = 0.1
    torque = [
        mean
        + ripple6 * math.cos(6 * 2.0 * math.pi * idx / samples)
        + ripple12 * math.sin(12 * 2.0 * math.pi * idx / samples)
        for idx in range(samples)
    ]

    health = torque_angle_sweep_health_summary(
        torque,
        max_harmonic=18,
        max_ac_rms_over_mean=0.04,
        allowed_dominant_harmonics=[6],
        min_mean_abs_torque_Nm=9.0,
        top_harmonics=2,
    )

    assert health["status"] == "ok"
    assert health["dominant_harmonic"] == 6
    assert health["top_harmonic_rows"][0]["order"] == 6
    assert health["top_harmonic_rows"][1]["order"] == 12
    assert math.isclose(health["ac_rms_over_mean"], math.sqrt((ripple6 * ripple6 + ripple12 * ripple12) / 2.0) / mean, rel_tol=1.0e-12)
    assert math.isclose(health["top_harmonic_rows"][0]["ac_variance_fraction"], (ripple6 * ripple6) / (ripple6 * ripple6 + ripple12 * ripple12), rel_tol=1.0e-12)

    bad = torque_angle_sweep_health_summary(
        torque,
        max_harmonic=18,
        max_ac_rms_over_mean=0.02,
        allowed_dominant_harmonics=[12],
    )
    assert bad["status"] == "needs_attention"
    assert bad["checks"]["ac_rms_over_mean"] is False
    assert bad["checks"]["dominant_harmonic"] is False


def test_torque_angle_sweep_comparison_summary_tracks_mean_and_ripple_deltas():
    samples = 720
    reference = []
    candidate = []
    for idx in range(samples):
        theta = 2.0 * math.pi * idx / samples
        reference.append(
            10.0
            + 0.4 * math.cos(6 * theta)
            + 0.1 * math.sin(12 * theta)
        )
        candidate.append(
            10.2
            + 0.3 * math.cos(6 * theta)
            + 0.12 * math.sin(12 * theta)
        )

    comparison = torque_angle_sweep_comparison_summary(
        reference,
        candidate,
        max_harmonic=18,
    )
    rows = {row["order"]: row for row in comparison["harmonic_delta_rows"]}

    expected_delta_rms = math.sqrt(0.2 * 0.2 + (0.1 * 0.1 + 0.02 * 0.02) / 2.0)
    expected_delta_ac_rms = math.sqrt((0.1 * 0.1 + 0.02 * 0.02) / 2.0)

    assert math.isclose(comparison["mean_delta_Nm"], 0.2, abs_tol=1.0e-13)
    assert math.isclose(comparison["sample_delta_rms_Nm"], expected_delta_rms, abs_tol=1.0e-13)
    assert math.isclose(comparison["difference_summary"]["ac_rms_torque_Nm"], expected_delta_ac_rms, abs_tol=1.0e-13)
    assert comparison["dominant_harmonic_changed"] is False
    assert comparison["worst_harmonic_order"] == 6
    assert math.isclose(rows[6]["amplitude_delta_Nm"], -0.1, abs_tol=1.0e-13)
    assert math.isclose(rows[6]["delta_waveform_amplitude_Nm"], 0.1, abs_tol=1.0e-13)
    assert math.isclose(rows[12]["amplitude_delta_Nm"], 0.02, abs_tol=1.0e-13)
    assert math.isclose(rows[12]["sin_coefficient_delta_Nm"], 0.02, abs_tol=1.0e-13)


def test_invalid_inputs():
    for bad in ({0: 1.0}, {-1: 1.0}):
        try:
            three_phase_torque_ripple_harmonics(bad)
        except ValueError:
            pass
        else:
            raise AssertionError("invalid harmonic order accepted")
    for kwargs in ({"current_peak": -1.0}, {"mechanical_speed": 0.0}):
        try:
            three_phase_torque_ripple_harmonics({1: 1.0}, **kwargs)
        except ValueError:
            pass
        else:
            raise AssertionError("invalid operating input accepted")
    for kwargs in ({"torque_Nm": [1.0, 2.0]}, {"torque_Nm": [1.0, 2.0, 3.0], "max_harmonic": 0}):
        try:
            torque_angle_sweep_summary(**kwargs)
        except ValueError:
            pass
        else:
            raise AssertionError("invalid torque sweep input accepted")
    try:
        torque_angle_sweep_health_summary([1.0, 2.0, 3.0], top_harmonics=-1)
    except ValueError:
        pass
    else:
        raise AssertionError("invalid torque health input accepted")
    try:
        torque_angle_sweep_comparison_summary([1.0, 2.0, 3.0], [1.0, 2.0])
    except ValueError:
        pass
    else:
        raise AssertionError("mismatched torque sweep lengths accepted")
