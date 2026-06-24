# -*- coding: utf-8 -*-
r"""Parametric modeling ops for build123d (radia_mcp.build123d.modeling) -- geometry-gated.

Each op is checked against a closed-form volume / count where possible, for OCCT validity, for the
region label, and -- the point of "CAE-safe" -- that the result MESHES cleanly in Netgen (the
build123d -> Netgen -> Radia/NGSolve tet pipeline).
"""
import math
import os
import sys

import pytest

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)
from radia_mcp.build123d.modeling import (annular_segment, tube, racetrack_coil, polar_array,
                                          linear_array, mirrored, assembly,
                                          shape_envelope_row, enclosing_box,
                                          enclosure_clearance_row, enclosure_difference_region,
                                          shape_measurement_row, shape_measurement_rows,
                                          box_face_vector_area_rows,
                                          box_face_pressure_force_rows,
                                          box_face_pressure_moment_rows,
                                          box_face_traction_moment_rows,
                                          compare_boundary_vector_area_rows,
                                          compare_shape_measurement_rows,
                                          shape_measurement_comparison_summary)
from build123d import Box, Compound, Pos


def _vol(obj):
    return sum(s.volume for s in obj.solids()) if isinstance(obj, Compound) else obj.volume


def _valid(obj):
    return all(s.is_valid for s in obj.solids()) if isinstance(obj, Compound) else obj.is_valid


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
