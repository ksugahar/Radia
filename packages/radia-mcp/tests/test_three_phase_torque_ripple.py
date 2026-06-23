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

from radia_mcp.radia_ngsolve.solve import three_phase_torque_ripple_harmonics


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


def test_closed_form_matches_direct_time_waveform_fourier():
    harmonics = {1: 0.8, 5: 0.12, 7: 0.04, 11: -0.03, 13: 0.01}
    current_peak = 1.7
    out = three_phase_torque_ripple_harmonics(harmonics, current_peak=current_peak)
    p = _sample_power(harmonics, current_peak)
    mean = sum(p) / len(p)

    assert math.isclose(mean, out["mean_power"], abs_tol=1e-12)
    assert math.isclose(_fourier_amplitude(p, 6), out["power_ripple"][6], abs_tol=1e-12)
    assert math.isclose(_fourier_amplitude(p, 12), out["power_ripple"][12], abs_tol=1e-12)


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
