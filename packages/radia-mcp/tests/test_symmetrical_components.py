r"""Fortescue symmetrical components for three-phase phasors.

This is the phasor diagnostic layer before dq/Park extraction: a balanced abc set is positive
sequence, the reversed set is negative sequence, and common-mode content is zero sequence.
"""
import math
import os
import sys

import pytest

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.radia_ngsolve.solve import (
    inverse_symmetrical_components,
    symmetrical_components,
)

ALPHA = complex(-0.5, math.sqrt(3.0) / 2.0)


def _close(z, w=0.0, tol=1e-12):
    assert abs(z - w) < tol


def test_balanced_positive_sequence():
    v = 1.2 - 0.4j
    seq = symmetrical_components(v, ALPHA * ALPHA * v, ALPHA * v)
    _close(seq["positive"], v)
    _close(seq["negative"])
    _close(seq["zero"])
    assert seq["negative_unbalance"] < 1e-12
    assert seq["zero_unbalance"] < 1e-12


def test_balanced_negative_and_zero_sequence():
    vneg = -0.3 + 0.7j
    seq_neg = symmetrical_components(vneg, ALPHA * vneg, ALPHA * ALPHA * vneg)
    _close(seq_neg["negative"], vneg)
    _close(seq_neg["positive"])
    _close(seq_neg["zero"])
    assert math.isinf(seq_neg["negative_unbalance"])

    v0 = 0.15 - 0.2j
    seq_zero = symmetrical_components(v0, v0, v0)
    _close(seq_zero["zero"], v0)
    _close(seq_zero["positive"])
    _close(seq_zero["negative"])
    assert math.isinf(seq_zero["zero_unbalance"])


def test_roundtrip_mixed_components_and_unbalance_factors():
    v0 = 0.08 - 0.03j
    v1 = 1.0 + 0.2j
    v2 = -0.12 + 0.18j
    phases = inverse_symmetrical_components(v0, v1, v2)
    seq = symmetrical_components(*phases)

    _close(seq["zero"], v0)
    _close(seq["positive"], v1)
    _close(seq["negative"], v2)
    assert math.isclose(seq["negative_unbalance"], abs(v2) / abs(v1), rel_tol=1e-12)
    assert math.isclose(seq["zero_unbalance"], abs(v0) / abs(v1), rel_tol=1e-12)
    for got, want in zip(inverse_symmetrical_components(seq["zero"], seq["positive"], seq["negative"]), phases):
        _close(got, want)


def test_invalid_base_unbalance_behavior():
    seq = symmetrical_components(0.0, 0.0, 0.0)
    assert seq["negative_unbalance"] == 0.0
    assert seq["zero_unbalance"] == 0.0

    # No positive sequence but nonzero negative/zero sequence: report infinite relative unbalance.
    assert math.isinf(symmetrical_components(1.0, ALPHA, ALPHA * ALPHA)["negative_unbalance"])
    assert math.isinf(symmetrical_components(2.0, 2.0, 2.0)["zero_unbalance"])


if __name__ == "__main__":
    pytest.main([__file__])
