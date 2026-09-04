"""Fast contracts for the coil-driven ESRF three-engine cases."""

from __future__ import annotations

import importlib.util
from pathlib import Path
import sys

import numpy as np
import pytest


MODULE_PATH = (
    Path(__file__).resolve().parents[1]
    / "validation_test"
    / "esrf_three_engine"
    / "esrf_coil_yoke.py"
)


@pytest.fixture(scope="module")
def contracts():
    spec = importlib.util.spec_from_file_location("_esrf_coil_yoke", MODULE_PATH)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


@pytest.mark.parametrize("number,axis", ((6, 0), (7, 2)))
def test_esrf_coil_yoke_observations_are_symmetric_gap_stencils(
    contracts, number: int, axis: int
):
    case = contracts.get_case(number)
    points = contracts.observation_points(number)
    assert case.beam_axis == axis
    assert points.shape == (45, 3)
    assert np.any(np.isclose(points[:, axis], 0.0))
    transverse = tuple(index for index in range(3) if index != axis)
    for coordinate in transverse:
        assert set(np.round(points[:, coordinate], 12)) == {
            -case.transverse_offsets_m[-1],
            0.0,
            case.transverse_offsets_m[-1],
        }
    core = contracts.core_selector(number, points)
    assert core.dtype == np.dtype(bool)
    assert 0 < int(np.count_nonzero(core)) < len(points)


def test_esrf_coil_yoke_rejects_a_non_coil_case(contracts):
    with pytest.raises(ValueError, match="6 and 7"):
        contracts.get_case(5)
