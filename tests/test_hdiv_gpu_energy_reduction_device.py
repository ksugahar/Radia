"""Fixed-order compensated material energy reduction on the device.

Guards: the material energy and maximum amplitude of the nonlinear HDiv-MMM Newton
loop are reduced by a parallel partition that depends on the point count only.
The device result equals, bitwise, an independent host replay of that fixed order,
repeats bitwise, stays within 2 ulp of the exactly rounded sum (math.fsum) on
cancelling and positive data, and keeps the serial loop's
maximum semantics (start at 0, fmax). An Armijo-style threshold placed 2 ulp on
either side of the exact energy is decided as the exact sum decides it.
Assumes: math.fsum is the exactly rounded reference; IEEE-754
binary64 round-to-nearest on host and device without contraction in the merge.
Optional GPU lane: skips visibly unless RADIA_HDIV_NATIVE_TEST_DIR is set.
"""

import math
import os

import numpy as np
import pytest

THREADS, MAX_BLOCKS = 256, 1024


@pytest.fixture(scope="module")
def cuda():
    if not os.environ.get("RADIA_HDIV_NATIVE_TEST_DIR"):
        pytest.skip("Optional GPU lane requires RADIA_HDIV_NATIVE_TEST_DIR")
    from radia.gpu import _native
    return _native("_hdiv_cuda")


def _merge(s, c, s2, c2):
    t = s + s2
    e = (s - t) + s2 if abs(s) >= abs(s2) else (s2 - t) + s
    return t, (c + c2) + e


def _fmax(a, b):
    return b if math.isnan(a) else a if math.isnan(b) else max(a, b)


def _tree(values):
    """Shared-memory tree of one 256-thread block, in the device's order."""
    ss, sc, sm = [list(x) for x in zip(*values)]
    k = THREADS // 2
    while k:
        for t in range(k):
            ss[t], sc[t] = _merge(ss[t], sc[t], ss[t + k], sc[t + k])
            sm[t] = _fmax(sm[t], sm[t + k])
        k //= 2
    return ss[0], sc[0], sm[0]


def _replay(energy, amplitude):
    """Host replay of the fixed device order: an oracle of order, not of accuracy."""
    n = len(energy)
    blocks = min(MAX_BLOCKS, -(-n // THREADS))
    stride = blocks * THREADS
    energy, amplitude = energy.tolist(), amplitude.tolist()
    partials = []
    for b in range(blocks):
        lanes = []
        for t in range(THREADS):
            s = c = m = 0.
            for q in range(b * THREADS + t, n, stride):
                s, c = _merge(s, c, energy[q], 0.)
                m = _fmax(m, amplitude[q])
            lanes.append((s, c, m))
        partials.append(_tree(lanes))
    lanes = []
    for t in range(THREADS):
        s = c = m = 0.
        for i in range(t, blocks, THREADS):
            s, c = _merge(s, c, partials[i][0], partials[i][1])
            m = _fmax(m, partials[i][2])
        lanes.append((s, c, m))
    s, c, m = _tree(lanes)
    return s + c, m, blocks


def _reduce(cuda, energy, amplitude):
    return cuda._diagnostic_material_energy_reduction(
        np.ascontiguousarray(energy, dtype=np.float64), np.ascontiguousarray(amplitude, dtype=np.float64), backend="gpu")


def _cases():
    rng = np.random.default_rng(20261112)
    for n in (1, 255, 256, 257, 4097, MAX_BLOCKS * THREADS + 3):
        positive = rng.random(n) * 10.0 ** rng.uniform(-6, 6, n)
        signs = np.where(rng.random(n) < .5, -1., 1.)
        cancelling = positive * signs
        if n > 4:
            cancelling[1], cancelling[3] = 1e16, -1e16
        yield n, "positive", positive, rng.random(n) * 1e6
        yield n, "cancelling", cancelling, rng.random(n) * 1e6


@pytest.mark.parametrize("n,kind,energy,amplitude", list(_cases()),
                         ids=lambda v: str(v) if isinstance(v, (int, str)) else "")
def test_device_reduction_matches_fixed_order_and_exact_sum(cuda, n, kind, energy, amplitude):
    report = _reduce(cuda, energy, amplitude)
    again = _reduce(cuda, energy, amplitude)
    replay, maximum, blocks = _replay(energy, amplitude)
    exact = math.fsum(energy)
    assert report["blocks"] == blocks, (report["blocks"], blocks)
    assert report["energy"] == again["energy"] and report["maximum"] == again["maximum"], (report, again)
    assert report["energy"] == replay, ("fixed order", n, kind, report["energy"], replay)
    assert report["maximum"] == maximum == max(0., float(amplitude.max())), (report["maximum"], maximum)
    ulp = math.ulp(exact)
    assert abs(report["energy"] - exact) <= 2 * ulp, ("exact sum", n, kind, report["energy"], exact, ulp)


def test_maximum_keeps_serial_fmax_semantics(cuda):
    energy = np.ones(300)
    amplitude = -np.ones(300)
    assert _reduce(cuda, energy, amplitude)["maximum"] == 0.
    amplitude[7] = np.nan
    amplitude[290] = 3.5
    assert _reduce(cuda, energy, amplitude)["maximum"] == 3.5


def test_armijo_threshold_two_ulp_either_side_follows_the_exact_energy(cuda):
    rng = np.random.default_rng(7)
    energy = rng.random(MAX_BLOCKS * THREADS + 11) * 10.0 ** rng.uniform(-3, 3, MAX_BLOCKS * THREADS + 11)
    exact = math.fsum(energy)
    value = _reduce(cuda, energy, np.zeros_like(energy))["energy"]
    above, below = exact + 2 * math.ulp(exact), exact - 2 * math.ulp(exact)
    assert value <= above, ("accept at exact + 2 ulp", value, above)
    assert not value <= below, ("reject at exact - 2 ulp", value, below)


@pytest.mark.parametrize("energy,amplitude", [(np.zeros(0), np.zeros(0)), (np.zeros(3), np.zeros(2)),
                                              (np.zeros((2, 2)), np.zeros((2, 2)))])
def test_diagnostic_refuses_bad_shapes(cuda, energy, amplitude):
    with pytest.raises(ValueError, match="matching nonempty"):
        _reduce(cuda, energy, amplitude)
