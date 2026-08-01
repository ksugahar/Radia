"""Regression gates for CoilBuilder field fidelity and yoke clearance."""

import radia as rad
from netgen.occ import Box, Pnt

from radia.coil_builder import (
    CoilBuilder,
    audit_coil_field_consistency,
    audit_coil_yoke_clearance,
)


def _racetrack():
    radius = 0.02
    return (
        CoilBuilder(current=1000.0)
        .set_start([0.08, -0.10, 0.05])
        .set_cross_section(0.004, 0.004)
        .add_straight(0.20)
        .add_arc(radius, 90)
        .add_straight(0.12)
        .add_arc(radius, 90)
        .add_straight(0.20)
        .add_arc(radius, 90)
        .add_straight(0.12)
        .add_arc(radius, 90)
    )


def test_adaptive_radia_arc_matches_finite_filament_field():
    rad.UtiDelAll()
    try:
        report = audit_coil_field_consistency(
            _racetrack(),
            [[0.0, 0.0, 0.0]],
            n_arc=400,
            relative_tolerance=1.0e-3,
        )
        assert report["passed"]
        assert report["closed"]
        assert report["field_consistent"]
        assert report["max_relative_error"] < 1.0e-3
    finally:
        rad.UtiDelAll()


def test_coil_yoke_clearance_passes_separation_and_rejects_overlap():
    coil = _racetrack()
    separated_yoke = Box(Pnt(-0.01, -0.01, -0.01), Pnt(0.01, 0.01, 0.01))
    separated = audit_coil_yoke_clearance(
        coil, separated_yoke, minimum_clearance=0.02
    )
    assert separated["passed"]
    assert separated["no_overlap"]
    assert separated["measured_clearance"] > 0.07

    intersecting_yoke = Box(
        Pnt(0.075, -0.02, 0.045), Pnt(0.085, 0.02, 0.055)
    )
    intersecting = audit_coil_yoke_clearance(coil, intersecting_yoke)
    assert not intersecting["passed"]
    assert not intersecting["no_overlap"]
    assert intersecting["intersection_volume"] > 0.0
