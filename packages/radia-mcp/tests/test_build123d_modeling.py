# -*- coding: utf-8 -*-
r"""Parametric modeling ops for build123d (radia_mcp.build123d.modeling) -- geometry-gated.

Each op is checked against a closed-form volume / count where possible, for OCCT validity, for the
region label, and -- the point of "CAE-safe" -- that the result MESHES cleanly in Netgen (the
build123d -> Netgen -> Radia/NGSolve tet pipeline).
"""
import math
import os
import sys
import json

import pytest

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)
from radia_mcp.build123d.modeling import (annular_segment, tube, racetrack_coil, polar_array,
                                          linear_array, mirrored, assembly,
                                          shape_envelope_row, enclosing_box,
                                          enclosure_clearance_row, enclosure_difference_region,
                                          shape_measurement_row, shape_measurement_rows,
                                          box_through_cylinder_reference_row,
                                          mounting_plate_boss_reference_row,
                                          keyed_terminal_plate_reference_row,
                                          flanged_sleeve_reference_row,
                                          coax_annular_sleeve_reference_row,
                                          ribbed_busbar_heat_sink_reference_row,
                                          three_phase_busbar_snubber_plate_reference_row,
                                          rcd_snubber_heat_spreader_reference_row,
                                          rcd_snubber_capacitance_sweep_rows,
                                          thermal_robin_cooling_plate_reference_row,
                                          v_type_ipm_rotor_coupon_reference_row,
                                          box_face_vector_area_rows,
                                          box_face_pressure_force_rows,
                                          box_face_pressure_moment_rows,
                                          box_face_pressure_resultant_summary,
                                          box_face_traction_moment_rows,
                                          compare_boundary_vector_area_rows,
                                          compare_shape_measurement_rows,
                                          compare_shape_volume_rows,
                                          shape_name_identity_gate,
                                          shape_role_metadata_gate,
                                          shape_transition_role_metadata_gate,
                                          shape_volume_crosscheck_summary,
                                          shape_mass_property_crosscheck_summary,
                                          shape_measurement_comparison_summary,
                                          shape_measurement_inventory_summary,
                                          worst_shape_measurement_comparison_rows,
                                          shape_measurement_health_summary,
                                          shape_bbox_pair_clearance_summary,
                                          shape_parameter_sweep_summary)
from radia_mcp.build123d.build123d_knowledge import get_build123d_documentation
from build123d import Box, Compound, Pos


def _vol(obj):
    return sum(s.volume for s in obj.solids()) if isinstance(obj, Compound) else obj.volume


def _valid(obj):
    return all(s.is_valid for s in obj.solids()) if isinstance(obj, Compound) else obj.is_valid


def test_build123d_lab_policy_routes_tet_to_netgen_and_mixed_to_cubit():
    doc = get_build123d_documentation("lab_policy")

    assert "tet-only path" in doc
    assert "build123d → Netgen (tet) → Radia" in doc
    assert "hex path + mixed-mesh CAD fallback" in doc
    assert "pyramid transition elements" in doc
    assert "cubit_vol_inventory" in doc
    assert "box_through_cylinder_reference_row" in doc
    assert "volume as the common currency" in doc
    assert "build123d_volume_crosscheck" in doc
    assert "shape_mass_property_crosscheck_summary" in doc
    assert "build123d_mass_property_crosscheck" in doc
    assert "coreform_cubit.com -nographics -batch" in doc
    assert '"box_hole", "volume": 22.994690350851265' in doc
    assert '"l_bracket_two_holes", "volume": 2.8982123980236905' in doc
    assert "boolean union/overlap accounting" in doc
    assert "mounting_plate_boss_reference_row" in doc
    assert "three_phase_busbar_snubber_plate_reference_row" in doc
    assert '"mounting_plate_boss_five_holes", "volume": 12.786811880091562' in doc
    assert "min_edge_over_characteristic" in doc
    assert "max_volume_rel_error = 0.0" in doc
    assert "stepped cylindrical spacer" in doc
    assert "rel_error = 5.854366888206992e-6" in doc
    assert "shape_volume_crosscheck_summary(..., rtol=1e-5)" in doc
    assert "shape_name_identity_gate" in doc
    assert "shape_role_metadata_gate" in doc
    assert "shape_transition_role_metadata_gate" in doc
    assert "hex_region" in doc
    assert "transition_kind` set to `\"pyramid\"" in doc
    assert "missing, extra, duplicate, or unnamed imported solids" in doc
    assert "keyed terminal plate" in doc
    assert '"keyed_terminal_plate_two_bosses", "volume": 9.364087557965556' in doc
    assert "flanged sleeve" in doc
    assert '"flanged_sleeve_four_bolt_holes", "volume": 5.085955762823052' in doc
    assert "coax_annular_sleeve_reference_row" in doc
    assert '"coax_annular_sleeve", "volume": 38.453094079939064' in doc
    assert "max_volume_rel_error = 3.693299710979559e-5" in doc
    assert "`1e-4` external-CAD volume gate" in doc
    assert "Ribbed busbar / heat-sink" in doc
    assert "ribbed_busbar_heat_sink_reference_row" in doc
    assert '"ribbed_busbar_heat_sink_four_holes", "volume": 13.617497357233166' in doc
    assert "RCD snubber heat-spreader" in doc
    assert "rcd_snubber_heat_spreader_reference_row" in doc
    assert '"rcd_snubber_heat_spreader", "volume": 15.52205629192717' in doc
    assert "RCD capacitance sweep design-table gate" in doc
    assert "rcd_snubber_capacitance_sweep_rows" in doc
    assert "parameter_key=\"capacitance_uF\"" in doc
    assert "a `0.114` volume mismatch" in doc
    assert "Thermal Robin cooling plate" in doc
    assert "thermal_robin_cooling_plate_reference_row" in doc
    assert '"thermal_robin_cooling_plate", "volume": 15.11290974078757' in doc
    assert "multi-line dict literals" in doc
    assert "V-type IPM rotor-coupon" in doc
    assert "v_type_ipm_rotor_coupon_reference_row" in doc
    assert '"v_type_ipm_rotor_coupon", "volume": 13.238339620676824' in doc


def test_annular_segment_volume_and_validity():
    seg = annular_segment(40, 55, 20, 0, 30, label="seg0")
    ana = (30/360)*math.pi*(55**2-40**2)*20
    assert abs(_vol(seg)-ana)/ana < 1e-6, f"wedge volume {_vol(seg)} vs {ana}"
    assert seg.is_valid and seg.label == "seg0"


def test_tube_volume():
    t = tube(40, 55, 20)
    assert abs(_vol(t) - math.pi*(55**2-40**2)*20)/_vol(t) < 1e-9 and t.is_valid


def test_racetrack_coil_builds():
    c = racetrack_coil(120, 80, 8, 15, 20, label="rc")
    assert c.is_valid and c.volume > 0 and c.label == "rc"


def test_polar_array_full_ring_sums_to_annulus():
    seg = annular_segment(40, 55, 20, 0, 30, label="seg")
    ring = polar_array(seg, 12, 360.0, label="halbach")
    assert len(ring.solids()) == 12 and _valid(ring)
    assert abs(_vol(ring) - math.pi*(55**2-40**2)*20)/_vol(ring) < 1e-6, "12 x 30deg = full ring"
    # build123d keeps region labels on the Compound's CHILDREN (.solids() flattens & drops them)
    assert [c.label for c in ring.children] == [f"halbach_{k:02d}" for k in range(12)]


def test_polar_array_partial_fan():
    fan = polar_array(Box(10, 4, 4).translate((30, 0, 0)), 4, total_angle=90.0)
    assert len(fan.solids()) == 4 and _valid(fan)


def test_linear_array_count_and_spacing():
    arr = linear_array(Box(10, 10, 10), 5, 15.0, (1, 0, 0), label="row")
    assert len(arr.solids()) == 5 and _valid(arr)
    assert abs(_vol(arr) - 5*1000.0) < 1e-6
    xs = sorted(s.center().X for s in arr.solids())
    assert abs((xs[-1]-xs[0]) - 4*15.0) < 1e-6, "5 copies span 4 spacings"


def test_mirrored_doubles_volume():
    half = Box(20, 10, 10).translate((0, 20, 0))
    both = mirrored(half, label="sym")
    assert len(both.solids()) == 2 and abs(_vol(both) - 2*2000.0) < 1e-6


def test_assembly_keeps_regions_separate():
    a = annular_segment(40, 55, 20, 0, 90, label="iron")
    b = tube(10, 20, 20, label="coil")
    asm = assembly(a, b, label="machine")
    assert len(asm.solids()) == 2 and asm.label == "machine"
    labels = {c.label for c in asm.children}
    assert labels == {"iron", "coil"}, "regions stay separate and labelled on the children (not fused)"


def test_assembly_accepts_raw_primitives():
    """A raw build123d Part (Box / Cylinder / extrude result) is a Compound with its solid in .solids()
    but EMPTY .children -- assembly must still keep it (fall back to .solids()), not silently drop it,
    and a label set on the raw Part must survive onto its child solid."""
    from build123d import Box, Cylinder, Pos
    core = Cylinder(0.5, 2); core.label = "core"
    asm = assembly(Box(1, 1, 1), Pos(3, 0, 0) * core, label="mixed")
    assert len(asm.solids()) == 2 and len(asm.children) == 2, "raw primitives are kept, not dropped"
    assert abs(asm.volume - (1.0 + math.pi * 0.5 ** 2 * 2)) < 1e-6
    assert "core" in {c.label for c in asm.children}, "a label on a raw Part carries onto its solid"


def test_shape_measurement_row_matches_box_geometry():
    box = Box(2, 3, 4).solid()
    box.label = "box"
    row = shape_measurement_row(box)

    assert row["name"] == "box"
    assert row["is_valid"]
    assert row["volume"] == pytest.approx(24.0)
    assert row["area"] == pytest.approx(52.0)
    assert row["faces"] == 6
    assert row["edges"] == 12
    assert row["vertices"] == 8
    assert row["solids"] == 1
    assert row["bounding_box"]["min"] == pytest.approx([-1.0, -1.5, -2.0])
    assert row["bounding_box"]["max"] == pytest.approx([1.0, 1.5, 2.0])
    assert row["bounding_box"]["center"] == pytest.approx([0.0, 0.0, 0.0])
    assert row["bounding_box"]["size"] == pytest.approx([2.0, 3.0, 4.0])
    assert row["bounding_box"]["diagonal"] == pytest.approx(math.sqrt(29.0))
    assert row["characteristic_length"] == pytest.approx(4.0)


def test_box_through_cylinder_reference_matches_build123d_mass_properties():
    from build123d import Align, Cylinder

    x, y, z, radius = 4.0, 3.0, 2.0, 0.4
    part = Box(x, y, z) - Cylinder(
        radius=radius,
        height=1.5 * z,
        align=(Align.CENTER, Align.CENTER, Align.CENTER),
    )
    part = part.solid()
    part.label = "box_hole"

    reference = [box_through_cylinder_reference_row(x, y, z, radius, axis="z", label="box_hole")]
    measured = [shape_measurement_row(part)]
    health = shape_measurement_health_summary(reference, measured, rtol=1.0e-12, bbox_atol=1.0e-12)

    assert health["status"] == "ok"
    assert health["comparison_summary"]["max_volume_rel_error"] < 1.0e-14
    assert health["comparison_summary"]["max_area_rel_error"] < 1.0e-14
    assert health["comparison_summary"]["max_bbox_abs_error"] == pytest.approx(0.0)
    assert reference[0]["policy"] == "analytic_box_through_cylinder_mass_property_reference"


def test_shape_measurement_rows_follow_assembly_children():
    left = Box(1, 2, 3).solid()
    left.label = "left"
    right = Box(2, 2, 1).solid()
    right.label = "right"
    asm = assembly(left, right, label="two_region")
    rows = shape_measurement_rows(asm)

    assert [row["name"] for row in rows] == ["left", "right"]
    assert [row["volume"] for row in rows] == pytest.approx([6.0, 4.0])
    assert [row["area"] for row in rows] == pytest.approx([22.0, 16.0])


def test_shape_measurement_inventory_summary_reports_assembly_fractions():
    left = (Pos(-2, 0, 0) * Box(2, 2, 2)).solid()
    left.label = "left"
    right = (Pos(2, 0, 0) * Box(1, 2, 2)).solid()
    right.label = "right"
    rows = shape_measurement_rows(assembly(left, right, label="two_region"))

    summary = shape_measurement_inventory_summary(rows)

    assert summary["n_shapes"] == 2
    assert summary["n_valid"] == 2
    assert summary["total_volume"] == pytest.approx(12.0)
    assert summary["total_area"] == pytest.approx(40.0)
    assert summary["bounding_box"]["min"] == pytest.approx([-3.0, -1.0, -1.0])
    assert summary["bounding_box"]["max"] == pytest.approx([2.5, 1.0, 1.0])
    assert summary["bbox_volume"] == pytest.approx(22.0)
    assert summary["bbox_fill_fraction"] == pytest.approx(12.0 / 22.0)
    assert summary["largest_volume_name"] == "left"
    assert summary["smallest_volume_name"] == "right"
    fractions = {row["name"]: row for row in summary["volume_fraction_rows"]}
    assert fractions["left"]["volume_fraction"] == pytest.approx(8.0 / 12.0)
    assert fractions["right"]["volume_fraction"] == pytest.approx(4.0 / 12.0)


def test_shape_bbox_pair_clearance_summary_flags_overlap_for_precise_check():
    left = Box(1, 1, 1).solid()
    left.label = "left"
    right = (Pos(2.0, 0, 0) * Box(1, 1, 1)).solid()
    right.label = "right"
    overlap = (Pos(0.4, 0, 0) * Box(1, 1, 1)).solid()
    overlap.label = "overlap"
    rows = shape_measurement_rows(assembly(left, right, overlap, label="three_region"))

    summary = shape_bbox_pair_clearance_summary(rows)
    pairs = {row["pair"]: row for row in summary["pair_rows"]}

    assert summary["status"] == "needs_attention"
    assert summary["n_pairs"] == 3
    assert summary["separated_pair_count"] == 2
    assert summary["bbox_overlap_pair_count"] == 1
    assert summary["touching_pair_count"] == 0
    assert summary["min_positive_gap"] == pytest.approx(0.6)

    assert pairs["left::right"]["status"] == "separated"
    assert pairs["left::right"]["axis_gaps"]["x"] == pytest.approx(1.0)
    assert pairs["left::overlap"]["status"] == "bbox_overlap_needs_precise_check"
    assert pairs["left::overlap"]["bbox_intersection_volume"] == pytest.approx(0.6)
    assert pairs["left::overlap"]["axis_overlaps"]["x"] == pytest.approx(0.6)


def test_shape_parameter_sweep_summary_tracks_monotonic_metrics_and_limits():
    rows = []
    for height in [1.0, 2.0, 3.0, 4.0]:
        box = Box(2.0, 3.0, height).solid()
        row = shape_measurement_row(box, name=f"h_{height:g}")
        row["height"] = height
        rows.append(row)

    summary = shape_parameter_sweep_summary(
        rows,
        "height",
        metric_keys=("volume", "area"),
        limits_by_metric={"volume": {"min": 12.0, "max": 24.0}},
    )
    metrics = {row["metric"]: row for row in summary["metric_rows"]}

    assert summary["parameter_values"] == [1.0, 2.0, 3.0, 4.0]
    assert summary["parameter_strictly_increasing"] is True
    assert summary["status"] == "needs_attention"
    assert summary["constraint_violation_count"] == 1
    assert summary["constraint_violations"][0]["kind"] == "below_min"
    assert metrics["volume"]["min"] == pytest.approx(6.0)
    assert metrics["volume"]["max"] == pytest.approx(24.0)
    assert metrics["volume"]["monotonic_non_decreasing"] is True
    assert metrics["volume"]["min_step_delta"] == pytest.approx(6.0)
    assert metrics["area"]["first"] == pytest.approx(22.0)
    assert metrics["area"]["last"] == pytest.approx(52.0)
    assert metrics["area"]["monotonic_non_decreasing"] is True

    clean = shape_parameter_sweep_summary(rows, "height", metric_keys=("volume",))
    assert clean["status"] == "ok"
    assert clean["ok_for_design_table"] is True

    duplicate = shape_parameter_sweep_summary(rows + [dict(rows[-1])], "height", metric_keys=("volume",))
    assert duplicate["duplicate_parameter_values"] is True
    assert duplicate["status"] == "needs_attention"


def test_box_face_vector_area_rows_match_centered_box():
    box = Box(2, 3, 5).solid()
    measurement = shape_measurement_row(box)
    rows = box_face_vector_area_rows(measurement["bounding_box"]["size"])
    by_name = {row["name"]: row for row in rows}

    assert sum(row["surface_area"] for row in rows) == pytest.approx(measurement["area"])
    assert by_name["xmin"]["vector_area"] == pytest.approx((-15.0, 0.0, 0.0))
    assert by_name["xmax"]["vector_area"] == pytest.approx((15.0, 0.0, 0.0))
    assert by_name["ymin"]["vector_area"] == pytest.approx((0.0, -10.0, 0.0))
    assert by_name["ymax"]["vector_area"] == pytest.approx((0.0, 10.0, 0.0))
    assert by_name["zmin"]["vector_area"] == pytest.approx((0.0, 0.0, -6.0))
    assert by_name["zmax"]["vector_area"] == pytest.approx((0.0, 0.0, 6.0))
    assert all(row["vector_area_norm_over_area"] == pytest.approx(1.0) for row in rows)


def test_compare_boundary_vector_area_rows_marks_direction_mismatch():
    reference = box_face_vector_area_rows((2, 3, 5))
    measured = [dict(row) for row in reference]
    measured[0]["vector_area"] = (15.0, 0.0, 0.0)
    measured[0]["unit_normal"] = (1.0, 0.0, 0.0)

    rows = compare_boundary_vector_area_rows(reference, measured, vector_atol=1.0e-12)
    by_name = {row["name"]: row for row in rows}

    assert not by_name["xmin"]["passed"]
    assert by_name["xmin"]["vector_abs_error"] == pytest.approx(30.0)
    assert by_name["xmin"]["unit_normal_abs_error"] == pytest.approx(2.0)
    assert by_name["xmax"]["passed"]


def test_box_face_pressure_force_rows_integrate_pressure_loads():
    uniform = box_face_pressure_force_rows(
        (2, 3, 5),
        {"xmin": 2.0, "xmax": 2.0, "ymin": 2.0, "ymax": 2.0, "zmin": 2.0, "zmax": 2.0},
    )
    total = [sum(row["force_N"][axis] for row in uniform) for axis in range(3)]
    assert total == pytest.approx([0.0, 0.0, 0.0])

    zmax = box_face_pressure_force_rows((2, 3, 5), {"zmax": 2.0}, default_pressure=0.0)
    by_name = {row["name"]: row for row in zmax}
    assert by_name["zmax"]["force_N"] == pytest.approx((0.0, 0.0, 12.0))
    assert by_name["zmax"]["force_magnitude_N"] == pytest.approx(12.0)
    assert by_name["zmax"]["pressure_source"] == "name"
    assert by_name["xmin"]["pressure_source"] == "default"

    by_index = box_face_pressure_force_rows((2, 3, 5), {6: 3.0}, default_pressure=0.0)
    assert {row["name"]: row for row in by_index}["zmax"]["force_N"] == pytest.approx((0.0, 0.0, 18.0))

    with pytest.raises(KeyError):
        box_face_pressure_force_rows((2, 3, 5), {"zmax": 2.0}, default_pressure=None)


def test_box_face_pressure_moment_rows_integrate_pivot_moments():
    uniform = box_face_pressure_moment_rows(
        (2, 3, 5),
        {"xmin": 2.0, "xmax": 2.0, "ymin": 2.0, "ymax": 2.0, "zmin": 2.0, "zmax": 2.0},
        center=(1.0, 1.5, 2.5),
    )
    total_force = [sum(row["force_N"][axis] for row in uniform) for axis in range(3)]
    total_moment = [sum(row["moment_about_pivot_Nm"][axis] for row in uniform) for axis in range(3)]
    assert total_force == pytest.approx((0.0, 0.0, 0.0))
    assert total_moment == pytest.approx((0.0, 0.0, 0.0))

    zmax = box_face_pressure_moment_rows(
        (2, 3, 5),
        {"zmax": 2.0},
        center=(1.0, 1.5, 2.5),
        default_pressure=0.0,
    )
    by_name = {row["name"]: row for row in zmax}
    assert by_name["zmax"]["face_center"] == pytest.approx((1.0, 1.5, 5.0))
    assert by_name["zmax"]["force_N"] == pytest.approx((0.0, 0.0, 12.0))
    assert by_name["zmax"]["moment_about_pivot_Nm"] == pytest.approx((18.0, -12.0, 0.0))

    shifted = box_face_pressure_moment_rows(
        (2, 3, 5),
        {"zmax": 2.0},
        center=(1.0, 1.5, 2.5),
        default_pressure=0.0,
        pivot_m=(1.0, 1.5, 0.0),
    )
    assert {row["name"]: row for row in shifted}["zmax"]["moment_about_pivot_Nm"] == pytest.approx((0.0, 0.0, 0.0))

    with pytest.raises(KeyError):
        box_face_pressure_moment_rows((2, 3, 5), {"zmax": 2.0}, default_pressure=None)


def test_box_face_pressure_resultant_summary_matches_closed_box_balance():
    box = (Pos(1.0, 1.5, 2.5) * Box(2, 3, 5)).solid()
    measurement = shape_measurement_row(box)
    size = measurement["bounding_box"]["size"]
    center = measurement["bounding_box"]["center"]

    uniform = box_face_pressure_resultant_summary(size, {}, center=center, default_pressure=2.0)
    assert uniform["boundary_count"] == 6
    assert uniform["box_size"] == pytest.approx((2.0, 3.0, 5.0))
    assert uniform["box_center"] == pytest.approx((1.0, 1.5, 2.5))
    assert uniform["total_force_N"] == pytest.approx((0.0, 0.0, 0.0))
    assert uniform["total_moment_about_pivot_Nm"] == pytest.approx((0.0, 0.0, 0.0))
    assert uniform["absolute_force_sum_N"] == pytest.approx(124.0)
    assert uniform["force_balance_ratio"] == pytest.approx(0.0)
    assert uniform["surface_vector_area"] == pytest.approx((0.0, 0.0, 0.0))

    zmax = box_face_pressure_resultant_summary(
        size,
        {"zmax": 2.0},
        center=center,
        default_pressure=0.0,
    )
    assert zmax["total_force_N"] == pytest.approx((0.0, 0.0, 12.0))
    assert zmax["total_moment_about_pivot_Nm"] == pytest.approx((18.0, -12.0, 0.0))
    assert zmax["force_balance_ratio"] == pytest.approx(1.0)

    shifted = box_face_pressure_resultant_summary(
        size,
        {"zmax": 2.0},
        center=center,
        default_pressure=0.0,
        pivot_m=(1.0, 1.5, 0.0),
    )
    assert shifted["total_moment_about_pivot_Nm"] == pytest.approx((0.0, 0.0, 0.0))

    with pytest.raises(KeyError):
        box_face_pressure_resultant_summary(size, {"zmax": 2.0}, default_pressure=None)


def test_box_face_traction_moment_rows_integrate_vector_tractions():
    rows = box_face_traction_moment_rows(
        (2, 3, 5),
        {"zmax": (1.0, -2.0, 3.0)},
        center=(1.0, 1.5, 2.5),
        default_traction=(0.0, 0.0, 0.0),
    )
    by_name = {row["name"]: row for row in rows}
    assert by_name["zmax"]["face_center"] == pytest.approx((1.0, 1.5, 5.0))
    assert by_name["zmax"]["traction_N_per_m2"] == pytest.approx((1.0, -2.0, 3.0))
    assert by_name["zmax"]["force_N"] == pytest.approx((6.0, -12.0, 18.0))
    assert by_name["zmax"]["moment_about_pivot_Nm"] == pytest.approx((87.0, 12.0, -21.0))
    assert by_name["zmax"]["traction_source"] == "name"
    assert by_name["xmin"]["traction_source"] == "default"

    shifted = box_face_traction_moment_rows(
        (2, 3, 5),
        {"zmax": (1.0, -2.0, 3.0)},
        center=(1.0, 1.5, 2.5),
        default_traction=(0.0, 0.0, 0.0),
        pivot_m=(1.0, 1.5, 5.0),
    )
    assert {row["name"]: row for row in shifted}["zmax"]["moment_about_pivot_Nm"] == pytest.approx((0.0, 0.0, 0.0))

    by_index = box_face_traction_moment_rows((2, 3, 5), {6: (0.0, 0.0, 2.0)}, default_traction=(0.0, 0.0, 0.0))
    assert {row["name"]: row for row in by_index}["zmax"]["force_N"] == pytest.approx((0.0, 0.0, 12.0))

    with pytest.raises(KeyError):
        box_face_traction_moment_rows((2, 3, 5), {"zmax": (1.0, -2.0, 3.0)}, default_traction=None)
    with pytest.raises(ValueError):
        box_face_traction_moment_rows((2, 3, 5), {"zmax": (1.0, -2.0)}, default_traction=(0.0, 0.0, 0.0))


def test_shape_envelope_row_and_enclosing_box_use_union_bbox_margin():
    left = (Pos(-2, 0, 0) * Box(2, 4, 6)).solid()
    right = (Pos(3, 0, 1) * Box(2, 2, 2)).solid()

    row = shape_envelope_row([left, right], margin=(1, 2, 3), name="outer")
    assert row["n_shapes"] == 2
    assert row["min"] == pytest.approx([-4.0, -4.0, -6.0])
    assert row["max"] == pytest.approx([5.0, 4.0, 6.0])
    assert row["size"] == pytest.approx([9.0, 8.0, 12.0])
    assert row["center"] == pytest.approx([0.5, 0.0, 0.0])
    assert row["volume"] == pytest.approx(864.0)

    outer = enclosing_box([left, right], margin=(1, 2, 3), label="outer")
    assert outer.label == "outer"
    assert outer.volume == pytest.approx(864.0)
    size = outer.bounding_box().size
    assert [size.X, size.Y, size.Z] == pytest.approx([9.0, 8.0, 12.0])


def test_enclosure_difference_region_volume_area_and_clearance():
    inner = Box(2, 2, 2).solid()
    inner.label = "solid"
    outer = enclosing_box([inner], margin=1.0, label="outer")

    clear = enclosure_clearance_row(outer, [inner])
    assert clear["contained_by_bbox"]
    assert clear["min_clearance"] == pytest.approx(1.0)
    assert clear["clearances"] == pytest.approx({
        "xmin": 1.0,
        "xmax": 1.0,
        "ymin": 1.0,
        "ymax": 1.0,
        "zmin": 1.0,
        "zmax": 1.0,
    })
    assert clear["enclosure_volume"] == pytest.approx(64.0)
    assert clear["inner_volume_sum"] == pytest.approx(8.0)
    assert clear["nominal_void_volume"] == pytest.approx(56.0)
    assert clear["inner_volume_fraction"] == pytest.approx(0.125)

    void = enclosure_difference_region(outer, [inner], label="void")
    assert void.label == "void"
    assert void.is_valid
    assert void.volume == pytest.approx(56.0)
    assert void.area == pytest.approx(96.0 + 24.0)


def test_compare_shape_measurement_rows_marks_pass_fail_and_missing():
    reference = [
        {"name": "ok", "volume": 10.0, "area": 20.0},
        {"name": "bad", "volume": 10.0, "area": 20.0},
        {"name": "missing", "volume": 1.0, "area": 2.0},
    ]
    measured = [
        {"name": "ok", "volume": 10.00001, "area": 20.00001},
        {"name": "bad", "volume": 10.2, "area": 20.0},
    ]
    rows = compare_shape_measurement_rows(reference, measured, rtol=1.0e-4, measured_label="external")
    by_name = {row["name"]: row for row in rows}

    assert by_name["ok"]["passed"]
    assert by_name["ok"]["measured_label"] == "external"
    assert by_name["ok"]["volume_rel_error"] < 1.0e-4
    assert not by_name["bad"]["passed"]
    assert by_name["bad"]["reason"] == "outside tolerance"
    assert not by_name["missing"]["passed"]
    assert by_name["missing"]["reason"] == "missing measured row"


def test_shape_volume_crosscheck_summary_accepts_cubit_and_external_cad_rows():
    reference = [
        {"name": "iron", "volume": 100.0, "area": 220.0},
        {"name": "coil", "volume": 25.0, "area": 60.0},
    ]
    measured = {
        "cubit": [
            {"name": "iron", "volume": 100.0000001},
            {"name": "coil", "volume": 24.9999999},
        ],
        "external_cad": {
            "rows": [
                {"name": "iron", "volume": 100.0},
                {"name": "coil", "volume": 25.0},
            ],
        },
    }

    summary = shape_volume_crosscheck_summary(reference, measured, rtol=1.0e-8)

    assert summary["policy"] == "build123d_external_cad_volume_crosscheck"
    assert summary["status"] == "ok"
    assert summary["ok_for_cad_roundtrip_volume"] is True
    assert summary["sources"] == ["cubit", "external_cad"]
    assert summary["max_volume_rel_error"] < 1.0e-8

    rows = compare_shape_volume_rows(reference, [{"name": "iron", "volume": 100.0}], measured_label="cubit")
    by_name = {row["name"]: row for row in rows}
    assert by_name["iron"]["passed"]
    assert not by_name["coil"]["passed"]
    assert by_name["coil"]["reason"] == "missing measured row"


def test_shape_name_identity_gate_rejects_extra_missing_and_duplicate_shapes():
    reference = [
        {"name": "stator", "volume": 100.0, "area": 220.0},
        {"name": "coil", "volume": 25.0, "area": 60.0},
    ]
    measured = [
        {"name": "stator", "volume": 100.0, "area": 220.0},
        {"name": "coil", "volume": 25.0, "area": 60.0},
    ]

    ok = shape_name_identity_gate(reference, measured, measured_label="cubit")
    assert ok["status"] == "ok"
    assert ok["policy"] == "build123d_cad_roundtrip_named_shape_identity_gate"
    assert ok["checks"]["same_name_multiset"] is True
    assert ok["version_note"].startswith("Use this before volume/area/bbox gates")

    bad = shape_name_identity_gate(
        reference,
        [
            {"name": "stator", "volume": 100.0, "area": 220.0},
            {"name": "stator", "volume": 100.0, "area": 220.0},
            {"name": "bolt", "volume": 1.0, "area": 2.0},
        ],
        measured_label="cst_import",
    )
    assert bad["status"] == "needs_attention"
    assert bad["measured_label"] == "cst_import"
    assert bad["missing_names"] == ["coil"]
    assert bad["extra_names"] == ["bolt", "stator"]
    assert bad["duplicate_measured_names"] == ["stator"]
    assert bad["checks"]["measured_names_unique"] is False

    summary = shape_mass_property_crosscheck_summary(
        reference,
        {"cst_import": [
            {"name": "stator", "volume": 100.0, "area": 220.0},
            {"name": "coil", "volume": 25.0, "area": 60.0},
            {"name": "bolt", "volume": 1.0, "area": 2.0},
        ]},
        rtol=1.0e-12,
    )
    assert summary["status"] == "needs_attention"
    assert summary["checks"]["all_sources_present_and_within_tolerance"] is True
    assert summary["checks"]["all_sources_preserve_named_shape_identity"] is False
    assert "missing, extra, duplicate, or unnamed shapes" in summary["issues"][0]
    assert summary["comparison_sets"][0]["name_identity_gate"]["extra_names"] == ["bolt"]


def test_shape_role_metadata_gate_requires_solver_handoff_semantics():
    rows = [
        {"name": "iron_core", "role": "magnetic_core", "material": "electrical_steel"},
        {"name": "phase_a_coil", "role": "conductor", "material_name": "copper"},
        {"name": "air_box", "solver_role": "air_region", "mat": "air"},
    ]

    ok = shape_role_metadata_gate(
        rows,
        required_names=["iron_core", "phase_a_coil"],
        required_roles=["magnetic_core", "conductor"],
        required_materials=["electrical_steel", "copper"],
        source_label="slot123_build123d",
    )

    assert ok["policy"] == "build123d_solver_handoff_role_material_metadata_gate"
    assert ok["status"] == "ok"
    assert ok["source_label"] == "slot123_build123d"
    assert ok["checks"]["required_roles_present"] is True
    assert ok["checks"]["required_materials_present"] is True
    assert ok["materials"] == ["air", "copper", "electrical_steel"]
    assert ok["version_note"].startswith("Run this after shape_name_identity_gate")

    bad = shape_role_metadata_gate(
        [
            {"name": "iron_core", "role": "magnetic_core", "material": "electrical_steel"},
            {"name": "iron_core", "role": "fixture", "material": "steel"},
            {"name": "phase_a_coil", "role": "conductor"},
            {"name": "air_box", "material": "air"},
        ],
        required_names=["iron_core", "phase_a_coil", "air_box", "shaft"],
        required_roles=["magnetic_core", "conductor", "air_region"],
        required_materials=["electrical_steel", "copper", "air"],
    )

    assert bad["status"] == "needs_attention"
    assert bad["duplicate_names"] == ["iron_core"]
    assert bad["rows_missing_material"] == ["phase_a_coil"]
    assert bad["rows_missing_role"] == ["air_box"]
    assert bad["missing_required_names"] == ["shaft"]
    assert bad["missing_required_roles"] == ["air_region"]
    assert bad["missing_required_materials"] == ["copper"]


def test_shape_transition_role_metadata_gate_preserves_hex_tet_handoff_intent():
    hex_body = (Pos(-1.25, 0, 0) * Box(1.0, 1.0, 1.0)).solid()
    hex_body.label = "hex_core"
    transition = Box(0.5, 1.0, 1.0).solid()
    transition.label = "pyramid_transition_envelope"
    tet_body = (Pos(1.25, 0, 0) * Box(1.0, 1.0, 1.0)).solid()
    tet_body.label = "tet_region"

    rows = shape_measurement_rows(assembly(hex_body, transition, tet_body, label="handoff"))
    by_name = {row["name"]: row for row in rows}
    by_name["hex_core"].update({"role": "hex_region", "material": "core_steel"})
    by_name["pyramid_transition_envelope"].update({
        "role": "mesh_transition",
        "material": "transition_air",
        "transition_kind": "pyramid",
        "connects_roles": ["hex_region", "tet_region"],
    })
    by_name["tet_region"].update({"role": "tet_region", "material": "air"})

    ok = shape_transition_role_metadata_gate(rows, source_label="slot131_build123d")

    assert ok["policy"] == "build123d_hex_tet_transition_role_metadata_gate"
    assert ok["status"] == "ok"
    assert ok["source_label"] == "slot131_build123d"
    assert ok["checks"]["required_roles_present"] is True
    assert ok["checks"]["transition_kind_matches"] is True
    assert ok["checks"]["transition_connects_required_roles"] is True
    assert ok["roles"] == ["hex_region", "mesh_transition", "tet_region"]
    assert ok["transition_kinds"] == ["pyramid"]
    assert ok["connected_roles"] == ["hex_region", "tet_region"]

    wrong_kind = [dict(row) for row in rows]
    for row in wrong_kind:
        if row["name"] == "pyramid_transition_envelope":
            row["transition_kind"] = "unknown"
    bad_kind = shape_transition_role_metadata_gate(wrong_kind)
    assert bad_kind["status"] == "needs_attention"
    assert bad_kind["checks"]["transition_kind_matches"] is False

    missing_connection = [dict(row) for row in rows]
    for row in missing_connection:
        if row["name"] == "pyramid_transition_envelope":
            row["connects_roles"] = ["hex_region"]
    bad_connection = shape_transition_role_metadata_gate(missing_connection)
    assert bad_connection["status"] == "needs_attention"
    assert bad_connection["checks"]["transition_connects_required_roles"] is False


def test_build123d_l_bracket_slot_volume_crosscheck_accepts_cubit_roundtrip():
    reference = [{"name": "l_bracket_two_holes", "volume": 2.898212398023691}]
    measured = {"cubit": [{"name": "l_bracket_two_holes", "volume": 2.8982123980236905}]}

    summary = shape_volume_crosscheck_summary(reference, measured, rtol=1.0e-12)

    assert summary["status"] == "ok"
    assert summary["ok_for_cad_roundtrip_volume"] is True
    assert summary["max_volume_rel_error"] == pytest.approx(1.5322866265870984e-16)


def test_mounting_plate_boss_reference_matches_slot19_volume_gate():
    reference = [mounting_plate_boss_reference_row(6.0, 4.0, 0.5, 0.8, 0.6, 0.25, 0.18, 2.2, 1.2)]
    measured = {
        "build123d": [{"name": "mounting_plate_boss_five_holes", "volume": 12.78681188009156}],
        "cubit": [{"name": "mounting_plate_boss_five_holes", "volume": 12.786811880091562}],
    }

    row = reference[0]
    summary = shape_volume_crosscheck_summary(reference, measured, rtol=1.0e-12)

    assert row["volume"] == pytest.approx(12.786811880091562)
    assert row["terms"]["base"] == pytest.approx(12.0)
    assert row["terms"]["central_hole"] < 0.0
    assert row["policy"] == "analytic_mounting_plate_boss_volume_reference"
    assert row["bounding_box"]["size"] == pytest.approx([6.0, 4.0, 1.1])
    assert summary["status"] == "ok"
    assert summary["sources"] == ["build123d", "cubit"]
    assert summary["max_volume_rel_error"] < 2.0e-16


def test_stepped_spacer_slot27_uses_external_kernel_roundtrip_tolerance():
    reference = [{"name": "stepped_spacer_four_bolt_holes", "volume": 4.8536349860900865}]
    measured = {
        "build123d": [{"name": "stepped_spacer_four_bolt_holes", "volume": 4.8536349860900865}],
        "cubit": [{"name": "stepped_spacer_four_bolt_holes", "volume": 4.853663401216389}],
    }

    strict_summary = shape_volume_crosscheck_summary(reference, measured, rtol=1.0e-6)
    roundtrip_summary = shape_volume_crosscheck_summary(reference, measured, rtol=1.0e-5)

    assert strict_summary["status"] == "needs_attention"
    assert strict_summary["comparison_sets"][0]["status"] == "ok"
    assert strict_summary["comparison_sets"][1]["status"] == "needs_attention"
    assert roundtrip_summary["status"] == "ok"
    assert roundtrip_summary["ok_for_cad_roundtrip_volume"] is True
    assert roundtrip_summary["max_volume_rel_error"] == pytest.approx(5.854366888206992e-6)


def test_keyed_terminal_plate_slot35_matches_volume_roundtrip_gate():
    reference = [keyed_terminal_plate_reference_row(
        8.0, 3.0, 0.4,
        0.55, 0.5, 2.5, 0.22,
        1.8, 0.8,
        0.7, 1.0,
        0.15, 3.4, -0.9,
    )]
    measured = {
        "build123d": [{"name": "keyed_terminal_plate_two_bosses", "volume": 9.364087557965554}],
        "cubit": [{"name": "keyed_terminal_plate_two_bosses", "volume": 9.364087557965554}],
    }

    row = reference[0]
    summary = shape_volume_crosscheck_summary(reference, measured, rtol=1.0e-12)

    assert row["volume"] == pytest.approx(9.364087557965556)
    assert row["terms"]["base"] == pytest.approx(9.6)
    assert row["terms"]["rectangular_window"] < 0.0
    assert row["terms"]["edge_key_slot"] < 0.0
    assert row["policy"] == "analytic_keyed_terminal_plate_volume_reference"
    assert row["bounding_box"]["size"] == pytest.approx([8.0, 3.0, 0.9])
    assert summary["status"] == "ok"
    assert summary["sources"] == ["build123d", "cubit"]
    assert summary["max_volume_rel_error"] < 2.0e-16


def test_flanged_sleeve_slot43_matches_volume_roundtrip_gate():
    reference = [flanged_sleeve_reference_row(
        1.8, 0.35,
        0.75, 1.1,
        0.28,
        1.25, 0.12,
        bolt_count=4,
    )]
    measured = {
        "build123d": [{"name": "flanged_sleeve_four_bolt_holes", "volume": 5.085955762823052}],
        "cubit": [{"name": "flanged_sleeve_four_bolt_holes", "volume": 5.085958320184421}],
    }

    row = reference[0]
    strict_summary = shape_volume_crosscheck_summary(reference, measured, rtol=1.0e-7)
    roundtrip_summary = shape_volume_crosscheck_summary(reference, measured, rtol=1.0e-6)

    assert row["volume"] == pytest.approx(5.085955762823052)
    assert row["terms"]["flange_annulus"] == pytest.approx(3.476360766756322)
    assert row["terms"]["hub_annulus"] == pytest.approx(1.6729295039631007)
    assert row["terms"]["bolt_holes"] < 0.0
    assert row["bolt_count"] == 4
    assert row["policy"] == "analytic_flanged_sleeve_volume_reference"
    assert row["bounding_box"]["size"] == pytest.approx([3.6, 3.6, 1.45])
    assert strict_summary["status"] == "needs_attention"
    assert roundtrip_summary["status"] == "ok"
    assert roundtrip_summary["sources"] == ["build123d", "cubit"]
    assert roundtrip_summary["max_volume_rel_error"] == pytest.approx(5.028278267134254e-7)


def test_coax_annular_sleeve_slot91_uses_physical_rc_geometry_volume_gate():
    reference = [coax_annular_sleeve_reference_row(0.9, 2.1, 3.4)]
    sleeve = tube(0.9, 2.1, 3.4, label="coax_annular_sleeve")
    measured = {
        "build123d": [{"name": "coax_annular_sleeve", "volume": _vol(sleeve)}],
        "cubit": [{"name": "coax_annular_sleeve", "volume": 38.451673891926546}],
    }

    row = reference[0]
    strict_summary = shape_volume_crosscheck_summary(reference, measured, rtol=1.0e-5)
    roundtrip_summary = shape_volume_crosscheck_summary(reference, measured, rtol=1.0e-4)

    assert row["volume"] == pytest.approx(38.453094079939064)
    assert _vol(sleeve) == pytest.approx(row["volume"])
    assert row["area"] == pytest.approx(86.70795723907828)
    assert row["terms"]["inner_void"] < 0.0
    assert row["parameters"] == {"inner_radius": 0.9, "outer_radius": 2.1, "height": 3.4}
    assert row["policy"] == "analytic_coax_annular_sleeve_volume_reference"
    assert row["bounding_box"]["size"] == pytest.approx([4.2, 4.2, 3.4])
    assert strict_summary["status"] == "needs_attention"
    assert roundtrip_summary["status"] == "ok"
    assert roundtrip_summary["sources"] == ["build123d", "cubit"]
    assert roundtrip_summary["max_volume_rel_error"] == pytest.approx(3.693299710979559e-5)


def test_ribbed_busbar_heat_sink_reference_matches_volume_roundtrip_gate():
    reference = [ribbed_busbar_heat_sink_reference_row(
        8.0, 3.2, 0.35,
        4, 0.2, 0.75, 0.35,
        0.18, 3.2, 1.25,
    )]
    measured = {
        "build123d": [{"name": "ribbed_busbar_heat_sink_four_holes", "volume": reference[0]["volume"]}],
        "cubit": [{"name": "ribbed_busbar_heat_sink_four_holes", "volume": reference[0]["volume"]}],
    }

    row = reference[0]
    summary = shape_volume_crosscheck_summary(reference, measured, rtol=1.0e-12)

    assert row["volume"] == pytest.approx(13.617497357233166)
    assert row["terms"]["base"] == pytest.approx(8.96)
    assert row["terms"]["straight_ribs"] == pytest.approx(4.8)
    assert row["terms"]["four_base_holes"] < 0.0
    assert row["fin_count"] == 4
    assert row["clearances"]["hole_to_fin_band"] > 0.0
    assert row["policy"] == "analytic_ribbed_busbar_heat_sink_volume_reference"
    assert row["bounding_box"]["size"] == pytest.approx([8.0, 3.2, 1.1])
    assert summary["status"] == "ok"
    assert summary["sources"] == ["build123d", "cubit"]


def test_three_phase_busbar_snubber_plate_slot59_matches_volume_roundtrip_gate():
    reference = [three_phase_busbar_snubber_plate_reference_row(
        9.0, 3.6, 0.35,
        3, 1.0, 0.6, 0.45, 2.4,
        2, 1.2, 0.45, 0.25, 2.0,
        0.16, 3.9, 1.35,
        phase_tab_y0=0.15,
        snubber_pad_y0=-1.1,
    )]
    measured = {
        "build123d": [{"name": "three_phase_busbar_snubber_plate", "volume": 12.307405319295338}],
        "cubit": [{"name": "three_phase_busbar_snubber_plate", "volume": 12.307405319295343}],
    }

    row = reference[0]
    summary = shape_volume_crosscheck_summary(reference, measured, rtol=1.0e-12)

    assert row["volume"] == pytest.approx(12.30740531929534)
    assert row["terms"]["base"] == pytest.approx(11.34)
    assert row["terms"]["three_phase_tabs"] == pytest.approx(0.81)
    assert row["terms"]["two_snubber_pads"] == pytest.approx(0.27)
    assert row["terms"]["four_mount_holes"] < 0.0
    assert row["counts"] == {"phase_tabs": 3, "snubber_pads": 2, "mount_holes": 4}
    assert row["clearances"]["mount_hole_to_y_edge"] == pytest.approx(0.29)
    assert row["clearances"]["phase_tab_gap"] == pytest.approx(1.4)
    assert row["policy"] == "analytic_three_phase_busbar_snubber_plate_volume_reference"
    assert row["bounding_box"]["size"] == pytest.approx([9.0, 3.6, 0.8])
    assert summary["status"] == "ok"
    assert summary["sources"] == ["build123d", "cubit"]
    assert summary["max_volume_rel_error"] < 5.0e-16


def test_rcd_snubber_heat_spreader_slot75_matches_volume_roundtrip_gate():
    reference = [rcd_snubber_heat_spreader_reference_row(
        10.0, 4.0, 0.32,
        5, 8.0, 0.12, 0.45, 0.45,
        2, 1.25, 0.70, 0.38, 2.4,
        0.16, 4.2, 1.55,
        snubber_pad_y0=-1.55,
    )]
    measured = {
        "build123d": [{"name": "rcd_snubber_heat_spreader", "volume": 15.522056291927165}],
        "cubit": [{"name": "rcd_snubber_heat_spreader", "volume": 15.522056291927173}],
    }

    row = reference[0]
    summary = shape_volume_crosscheck_summary(reference, measured, rtol=1.0e-12)

    assert row["volume"] == pytest.approx(15.52205629192717)
    assert row["terms"]["base"] == pytest.approx(12.8)
    assert row["terms"]["straight_ribs"] == pytest.approx(2.16)
    assert row["terms"]["snubber_pads"] == pytest.approx(0.665)
    assert row["terms"]["four_mount_holes"] < 0.0
    assert row["counts"] == {"ribs": 5, "snubber_pads": 2, "mount_holes": 4}
    assert row["clearances"]["snubber_pad_to_rib_band"] == pytest.approx(0.24)
    assert row["clearances"]["mount_hole_to_rib_band"] == pytest.approx(0.43)
    assert row["parameters"]["rib_pitch"] == pytest.approx(0.45)
    assert row["parameters"]["snubber_pad_y0"] == pytest.approx(-1.55)
    assert row["policy"] == "analytic_rcd_snubber_heat_spreader_volume_reference"
    assert row["bounding_box"]["size"] == pytest.approx([10.0, 4.0, 0.77])
    assert summary["status"] == "ok"
    assert summary["sources"] == ["build123d", "cubit"]
    assert summary["max_volume_rel_error"] < 6.0e-16


def test_rcd_snubber_capacitance_sweep_keeps_cad_variant_provenance():
    rows = rcd_snubber_capacitance_sweep_rows(
        [0.047, 0.10, 0.22],
        [1.05, 1.25, 1.55],
    )

    assert [row["capacitance_uF"] for row in rows] == pytest.approx([0.047, 0.10, 0.22])
    assert [row["snubber_pad_x"] for row in rows] == pytest.approx([1.05, 1.25, 1.55])
    assert rows[0]["name"] == "rcd_snubber_heat_spreader_0.047uF"
    assert rows[-1]["terms"]["snubber_pads"] > rows[0]["terms"]["snubber_pads"]
    assert rows[-1]["volume"] > rows[0]["volume"]
    assert all(row["design_table_role"].startswith("RCD snubber capacitance") for row in rows)

    summary = shape_parameter_sweep_summary(
        rows,
        "capacitance_uF",
        metric_keys=("volume", "snubber_pad_volume"),
    )
    metrics = {row["metric"]: row for row in summary["metric_rows"]}

    assert summary["status"] == "ok"
    assert summary["parameter_strictly_increasing"] is True
    assert metrics["volume"]["monotonic_non_decreasing"] is True
    assert metrics["snubber_pad_volume"]["monotonic_non_decreasing"] is True
    assert metrics["snubber_pad_volume"]["delta_first_to_last"] == pytest.approx(
        rows[-1]["snubber_pad_volume"] - rows[0]["snubber_pad_volume"]
    )


def test_thermal_robin_cooling_plate_slot83_matches_volume_roundtrip_gate():
    reference = [thermal_robin_cooling_plate_reference_row(
        9.0, 4.5, 0.30,
        4, 7.2, 0.14, 0.55, 0.42,
        2, 1.6, 0.80, 0.32, 2.6,
        0.14, 3.8, 1.8,
        fin_y0=0.45,
        device_pad_y0=-1.45,
    )]
    measured = {
        "build123d": [{"name": "thermal_robin_cooling_plate", "volume": 15.112909740787572}],
        "cubit": [{"name": "thermal_robin_cooling_plate", "volume": 15.112909740787561}],
    }

    row = reference[0]
    summary = shape_volume_crosscheck_summary(reference, measured, rtol=1.0e-12)

    assert row["volume"] == pytest.approx(15.11290974078757)
    assert row["terms"]["base"] == pytest.approx(12.15)
    assert row["terms"]["straight_cooling_fins"] == pytest.approx(2.2176)
    assert row["terms"]["device_pads"] == pytest.approx(0.8192)
    assert row["terms"]["four_mount_holes"] < 0.0
    assert row["counts"] == {"fins": 4, "device_pads": 2, "mount_holes": 4}
    assert row["clearances"]["device_pad_to_fin_band"] == pytest.approx(0.8)
    assert row["clearances"]["mount_hole_to_y_edge"] == pytest.approx(0.31)
    assert row["parameters"]["fin_y0"] == pytest.approx(0.45)
    assert row["parameters"]["device_pad_y0"] == pytest.approx(-1.45)
    assert row["policy"] == "analytic_thermal_robin_cooling_plate_volume_reference"
    assert row["bounding_box"]["size"] == pytest.approx([9.0, 4.5, 0.85])
    assert summary["status"] == "ok"
    assert summary["sources"] == ["build123d", "cubit"]
    assert summary["max_volume_rel_error"] < 8.0e-16


def test_v_type_ipm_rotor_coupon_reference_matches_volume_roundtrip_gate():
    reference = [v_type_ipm_rotor_coupon_reference_row(
        8.0, 5.0, 0.35,
        2.2, 0.35, 28.0,
        1.65, 0.55,
        0.45,
    )]
    measured = {
        "build123d": [{"name": "v_type_ipm_rotor_coupon", "volume": reference[0]["volume"]}],
        "cubit": [{"name": "v_type_ipm_rotor_coupon", "volume": reference[0]["volume"]}],
    }

    row = reference[0]
    summary = shape_volume_crosscheck_summary(reference, measured, rtol=1.0e-12)

    assert row["volume"] == pytest.approx(13.238339620676824)
    assert row["terms"]["coupon"] == pytest.approx(14.0)
    assert row["terms"]["two_v_magnet_pockets"] == pytest.approx(-0.539)
    assert row["terms"]["central_bore"] == pytest.approx(-0.2226603793231766)
    assert row["counts"] == {"magnet_pockets": 2, "bore": 1}
    assert row["clearances"]["pocket_to_bore_x"] == pytest.approx(0.1466001243676493)
    assert row["clearances"]["mirrored_pocket_gap"] == pytest.approx(1.1932002487352986)
    assert row["parameters"]["magnet_slot_angle_deg"] == pytest.approx(28.0)
    assert row["policy"] == "analytic_v_type_ipm_rotor_coupon_volume_reference"
    assert row["bounding_box"]["size"] == pytest.approx([8.0, 5.0, 0.35])
    assert summary["status"] == "ok"
    assert summary["sources"] == ["build123d", "cubit"]


def test_build123d_volume_crosscheck_mcp_tool_dispatches_json():
    from radia_mcp.build123d.server import build123d_volume_crosscheck

    reference = [{"name": "box", "volume": 24.0}]
    measured = {"cubit": [{"name": "box", "volume": 24.0}]}
    payload = json.loads(build123d_volume_crosscheck(json.dumps(reference), json.dumps(measured)))

    assert payload["status"] == "ok"
    assert payload["n_sources"] == 1
    assert payload["comparison_sets"][0]["source"] == "cubit"
    assert payload["comparison_sets"][0]["rows"][0]["passed"] is True


def test_build123d_mass_property_crosscheck_mcp_tool_dispatches_json():
    from radia_mcp.build123d.server import build123d_mass_property_crosscheck

    reference = [{
        "name": "box",
        "volume": 24.0,
        "area": 52.0,
        "is_valid": True,
        "bounding_box": {
            "min": [-1.0, -1.5, -2.0],
            "max": [1.0, 1.5, 2.0],
            "center": [0.0, 0.0, 0.0],
            "size": [2.0, 3.0, 4.0],
        },
    }]
    measured = {"cubit": [dict(reference[0])], "cst_import": [dict(reference[0])]}
    payload = json.loads(build123d_mass_property_crosscheck(json.dumps(reference), json.dumps(measured)))

    assert payload["status"] == "ok"
    assert payload["n_sources"] == 2
    assert payload["sources"] == ["cubit", "cst_import"]
    assert payload["comparison_sets"][0]["rows"][0]["passed"] is True


def test_compare_shape_measurement_rows_compares_bbox_when_present():
    reference = [{
        "name": "box",
        "volume": 1.0,
        "area": 6.0,
        "bounding_box": {
            "min": [0.0, 0.0, 0.0],
            "max": [1.0, 1.0, 1.0],
            "center": [0.5, 0.5, 0.5],
            "size": [1.0, 1.0, 1.0],
        },
    }]
    measured = [{
        "name": "box",
        "volume": 1.0,
        "area": 6.0,
        "bounding_box": {
            "min": [0.0, 0.0, 0.0],
            "max": [1.0, 1.0, 1.0000002],
            "center": [0.5, 0.5, 0.5000001],
            "size": [1.0, 1.0, 1.0000002],
        },
    }]

    row = compare_shape_measurement_rows(reference, measured, bbox_atol=1.0e-5)[0]
    assert row["passed"]
    assert row["bbox_compared"]
    assert row["bbox_abs_error"] == pytest.approx(2.0e-7)

    measured[0]["bounding_box"]["max"] = [1.0, 1.0, 1.01]
    measured[0]["bounding_box"]["size"] = [1.0, 1.0, 1.01]
    bad = compare_shape_measurement_rows(reference, measured, bbox_atol=1.0e-5)[0]
    assert not bad["passed"]
    assert bad["reason"] == "bbox outside tolerance"


def test_shape_measurement_comparison_summary_compacts_errors():
    reference = [{"name": "box", "volume": 24.0, "area": 52.0}]
    measured = [{"name": "box", "volume": 24.0, "area": 52.0}]
    summary = shape_measurement_comparison_summary(reference, measured, measured_label="cubit")

    assert summary["measured_label"] == "cubit"
    assert summary["n_cases"] == 1
    assert summary["n_passed"] == 1
    assert summary["n_bbox_compared"] == 0
    assert summary["max_volume_rel_error"] == pytest.approx(0.0)
    assert summary["max_area_rel_error"] == pytest.approx(0.0)
    assert summary["max_bbox_abs_error"] == pytest.approx(0.0)


def test_shape_measurement_health_summary_reports_worst_mismatches():
    reference = [
        {"name": "ok", "volume": 10.0, "area": 20.0, "is_valid": True},
        {"name": "bad", "volume": 10.0, "area": 20.0, "is_valid": True},
        {"name": "missing", "volume": 1.0, "area": 2.0, "is_valid": True},
    ]
    measured = [
        {"name": "ok", "volume": 10.000001, "area": 20.000001},
        {"name": "bad", "volume": 10.2, "area": 20.0},
    ]

    health = shape_measurement_health_summary(
        reference,
        measured,
        rtol=1.0e-4,
        measured_label="external",
        worst_limit=2,
    )

    assert health["status"] == "needs_attention"
    assert health["ok_for_geometry_roundtrip"] is False
    assert health["checks"]["all_reference_shapes_valid"] is True
    assert health["checks"]["all_measurements_present_and_within_tolerance"] is False
    assert health["comparison_summary"]["n_passed"] == 1
    assert [row["name"] for row in health["worst_comparisons"]] == ["missing", "bad"]
    assert worst_shape_measurement_comparison_rows([], limit=0) == []
    with pytest.raises(ValueError, match="limit must be non-negative"):
        worst_shape_measurement_comparison_rows([], limit=-1)


def test_shape_mass_property_crosscheck_summary_handles_multiple_cad_sources():
    box = Box(2, 3, 5).solid()
    reference = [shape_measurement_row(box, name="cad_block")]
    measured = {
        "cubit": [dict(reference[0])],
        "cst_import": [dict(reference[0])],
    }

    summary = shape_mass_property_crosscheck_summary(
        reference,
        measured,
        rtol=1.0e-12,
        bbox_atol=1.0e-12,
    )

    assert summary["policy"] == "build123d_external_cad_volume_area_bbox_crosscheck"
    assert summary["status"] == "ok"
    assert summary["sources"] == ["cubit", "cst_import"]
    assert summary["n_sources"] == 2
    assert summary["n_failed_rows"] == 0
    assert summary["max_area_rel_error"] == pytest.approx(0.0)
    assert all(row["n_bbox_compared"] == 1 for row in summary["comparison_sets"])

    bad_cst = dict(reference[0])
    bad_cst["area"] *= 1.01
    bad_summary = shape_mass_property_crosscheck_summary(
        reference,
        {"cubit": [dict(reference[0])], "cst_import": [bad_cst]},
        rtol=1.0e-4,
        bbox_atol=1.0e-12,
    )
    by_source = {row["source"]: row for row in bad_summary["comparison_sets"]}

    assert bad_summary["status"] == "needs_attention"
    assert bad_summary["ok_for_cad_roundtrip_mass_properties"] is False
    assert by_source["cubit"]["status"] == "ok"
    assert by_source["cst_import"]["status"] == "needs_attention"
    assert by_source["cst_import"]["worst_comparisons"][0]["reason"] == "outside tolerance"


def test_segment_meshes_in_netgen():
    """CAE gate: the wedge tet-meshes cleanly through the build123d -> STEP -> Netgen pipeline
    (build123d's OCP kernel and Netgen's OCC binding are different, so STEP is the bridge)."""
    import tempfile
    from build123d import export_step
    from netgen.occ import OCCGeometry
    from ngsolve import Mesh
    seg = annular_segment(40, 55, 20, 0, 45)
    with tempfile.TemporaryDirectory() as d:
        f = os.path.join(d, "seg.step")
        export_step(seg, f)
        mesh = Mesh(OCCGeometry(f).GenerateMesh(maxh=5.0))
    assert mesh.ne > 50, f"wedge should tet-mesh (got {mesh.ne} elements)"


def main():
    test_annular_segment_volume_and_validity()
    test_tube_volume()
    test_racetrack_coil_builds()
    test_polar_array_full_ring_sums_to_annulus()
    test_polar_array_partial_fan()
    test_linear_array_count_and_spacing()
    test_mirrored_doubles_volume()
    test_assembly_keeps_regions_separate()
    test_segment_meshes_in_netgen()
    print("[OK] build123d modeling ops: annular wedge, tube, racetrack coil, polar/linear arrays, "
          "mirror, labelled assembly -- analytic volumes, OCCT-valid, region labels, Netgen-meshable.")


if __name__ == "__main__":
    main()
