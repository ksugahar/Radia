"""Regression gates for CoilBuilder field fidelity and yoke clearance."""

import radia as rad
from netgen.occ import Box, Pnt

from radia.coil_builder import (
    CoilBuilder,
    audit_coil_field_consistency,
    audit_coil_yoke_clearance,
)


def _racetrack(current=1000.0, width=0.004):
    radius = 0.02
    return (
        CoilBuilder(current=current)
        .set_start([0.08, -0.10, 0.05])
        .set_cross_section(width, 0.004)
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


def test_rounded_rectangle_racetrack_has_closed_digest_bound_si_identity():
    coil = CoilBuilder.rounded_rectangle_racetrack(
        250.0,
        corner_centerline_radius=0.012,
        straight_x=0.08,
        straight_y=0.12,
        width=0.004,
        height=0.006,
        centre=(0.01, -0.02, 0.03),
    )
    manifest = coil.to_identity_manifest(turns=8, region_label="drive_coil")

    assert manifest["closed"] is True
    assert manifest["closure_error_m"] < 1.0e-12
    assert manifest["ampere_turns_A"] == 2000.0
    assert [item["kind"] for item in manifest["segments"]] == [
        "straight", "arc", "straight", "arc",
        "straight", "arc", "straight", "arc",
    ]
    assert len(manifest["identity_sha256"]) == 64


def test_coil_identity_changes_with_current_orientation_and_cross_section():
    base = _racetrack().to_identity_manifest(region_label="coil")
    changed_current = _racetrack(current=1001.0)
    assert changed_current.to_identity_manifest()["identity_sha256"] != base["identity_sha256"]

    changed_section = _racetrack(width=0.00404)
    assert changed_section.to_identity_manifest()["identity_sha256"] != base["identity_sha256"]
