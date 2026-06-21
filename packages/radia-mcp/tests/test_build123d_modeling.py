# -*- coding: utf-8 -*-
r"""Parametric modeling ops for build123d (radia_mcp.build123d.modeling) -- geometry-gated.

Each op is checked against a closed-form volume / count where possible, for OCCT validity, for the
region label, and -- the point of "CAE-safe" -- that the result MESHES cleanly in Netgen (the
build123d -> Netgen -> Radia/NGSolve tet pipeline).
"""
import math
import os
import sys

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)
from radia_mcp.build123d.modeling import (annular_segment, tube, racetrack_coil, polar_array,
                                          linear_array, mirrored, assembly)
from build123d import Box, Compound


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
