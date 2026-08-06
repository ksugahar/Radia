"""Fast source-level gates for the TEAM 28 CoilBuilder example."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
MAGLEV_DIR = REPO_ROOT / "validation_test" / "maglev"
if str(MAGLEV_DIR) not in sys.path:
    sys.path.insert(0, str(MAGLEV_DIR))

from team28_coilbuilder_eddy_bubble import (  # noqa: E402
    build_team28_coils,
    compare_coil_fields,
)


def test_team28_coilbuilder_paths_are_closed_and_counter_wound():
    coils = build_team28_coils(20.0)

    assert len(coils) == 2
    assert all(coil.is_closed for coil in coils)
    assert max(coil.gap for coil in coils) < 1.0e-12
    assert coils[0].current == 960.0 * 20.0
    assert coils[1].current == -576.0 * 20.0


def test_team28_coilbuilder_fields_match_independent_winding_quadrature():
    points = np.asarray(
        [
            [0.012, 0.003, 0.0120],
            [0.030, -0.010, 0.0130],
            [0.050, 0.015, 0.0140],
            [0.063, -0.004, 0.0125],
        ]
    )

    _, _, report = compare_coil_fields(
        points,
        coil_current_A=20.0,
        arc_max_segment_length_m=0.002,
    )
    cross_check = report["field_cross_check"]

    assert report["all_paths_closed"]
    assert cross_check["vector_potential_relative_l2"] < 5.0e-4
    assert cross_check["flux_density_relative_l2"] < 1.0e-3
    assert cross_check["flux_density_max_pointwise_over_reference_max"] < 2.0e-3
