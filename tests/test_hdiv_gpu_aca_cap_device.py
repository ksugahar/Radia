"""A device ACA leaf that settles at its rank cap still answers its probe rows.

Guards: the stop of `acaStep` at the rank cap: a settled last term must not mark the leaf converged before the unchecked rows
agree, since a row the factors miss then stays out of the stored block.
Assumes: leaves described as diag(d) + w w^T, evaluated on the device. Status 1 is converged and 2 the rank cap, which stores
the leaf dense. Probe rows are samples, so this covers the stopping path and is no bound on the error of a leaf.
Optional GPU lane: skips visibly unless RADIA_HDIV_NATIVE_TEST_DIR is set.
"""

import os

import numpy as np
import pytest

CONVERGED, RANK_CAP = 1, 2
EPS = 1e-5
DELTA = 1e-6


@pytest.fixture(scope="module")
def cuda():
    if not os.environ.get("RADIA_HDIV_NATIVE_TEST_DIR"):
        pytest.skip("Optional GPU lane requires RADIA_HDIV_NATIVE_TEST_DIR")
    from radia.gpu import _native
    return _native("_hdiv_cuda")


def _leaf(cuda, diagonal, outer, cap):
    return cuda._diagnostic_aca_leaf(np.array(diagonal, dtype=float), np.array(outer, dtype=float), EPS, cap, "gpu")


def test_a_row_missed_at_the_cap_sends_the_leaf_dense(cuda):
    report = _leaf(cuda, [1.0, DELTA, 1.0], [0.0, 0.0, 0.0], 2)
    assert (report["rank"], report["status"]) == (2, RANK_CAP), f"pivots 1, delta, 1 at cap 2 gave {report}, expected rank 2 and status {RANK_CAP}"


def test_room_for_the_missed_row_reaches_it(cuda):
    report = _leaf(cuda, [1.0, DELTA, 1.0], [0.0, 0.0, 0.0], 3)
    assert (report["rank"], report["status"]) == (3, CONVERGED), f"pivots 1, delta, 1 at cap 3 gave {report}, expected rank 3 and status {CONVERGED}"


def test_rows_that_agree_at_the_cap_keep_the_factors(cuda):
    outer = np.cos(np.linspace(0.3, 2.5, 32))
    diagonal = np.zeros(32)
    diagonal[5] = DELTA
    report = _leaf(cuda, diagonal, outer, 2)
    assert (report["rank"], report["status"]) == (2, CONVERGED), f"rank-2 leaf at cap 2 gave {report}, expected rank 2 and status {CONVERGED}"


def test_the_diagnostic_refuses_mismatched_vectors(cuda):
    with pytest.raises(ValueError, match="matching vectors"):
        cuda._diagnostic_aca_leaf(np.zeros(3), np.zeros(2), EPS, 2, "gpu")
    with pytest.raises(ValueError, match="backend"):
        cuda._diagnostic_aca_leaf(np.zeros(3), np.zeros(3), EPS, 2, "cpu")
