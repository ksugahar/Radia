# -*- coding: utf-8 -*-
r"""EM-device archetypes for build123d (radia_mcp.build123d.archetypes) -- geometry + physics gated.

The Halbach ring's per-segment easy-axis angles must follow the Mallinson (pole_pairs+1)*theta law and
round-trip through the magnetization-label convention; the C-core must actually have its pole gap; the
PM/coil primitives carry their magnetization/region labels; and a representative archetype tet-meshes
through the build123d -> STEP -> Netgen pipeline.
"""
import math
import os
import sys

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)
from radia_mcp.build123d.archetypes import (magnetization_tag, parse_magnetization, cylindrical_magnet,
                                            block_magnet, halbach_ring, c_core, solenoid)
from build123d import Box


def test_magnetization_label_roundtrip():
    for ang in (0.0, 37.5, 359.9, 720.0 + 12.0):
        lab = magnetization_tag("seg", 3, ang)
        assert abs(parse_magnetization(lab) - (ang % 360.0)) < 1e-6, lab
    assert parse_magnetization("plain_iron") is None


def test_pm_primitives_carry_magnetization():
    cyl = cylindrical_magnet(10, 20, m_angle_deg=45.0, name="pm", index=2)
    blk = block_magnet(10, 8, 20, m_angle_deg=90.0, name="pm", index=5)
    assert cyl.is_valid and abs(cyl.volume - math.pi*10**2*20)/cyl.volume < 1e-9
    assert abs(parse_magnetization(cyl.label) - 45.0) < 1e-6
    assert blk.is_valid and abs(blk.volume - 10*8*20) < 1e-9
    assert abs(parse_magnetization(blk.label) - 90.0) < 1e-6


def test_halbach_ring_mallinson_easy_axes():
    n, p = 12, 1                                          # dipole Halbach
    hb = halbach_ring(40, 55, 20, n, pole_pairs=p, name="hb")
    assert len(hb.children) == n
    vol = sum(s.volume for s in hb.solids())
    assert abs(vol - math.pi*(55**2-40**2)*20)/vol < 1e-6, "segments tile the full ring"
    angles = [parse_magnetization(c.label) for c in hb.children]
    for k, a in enumerate(angles):
        theta_c = k*360.0/n + 0.5*360.0/n
        assert abs(a - ((p+1)*theta_c) % 360.0) < 1e-6, f"seg {k}: easy axis {a} != Mallinson"
    # dipole signature: the easy axis advances at TWICE the mechanical step
    step = (angles[1]-angles[0]) % 360.0
    assert abs(step - 2*(360.0/n)) < 1e-6, "dipole Halbach: easy axis advances 2x mechanical angle"


def test_halbach_quadrupole_advances_3x():
    n = 16
    hb = halbach_ring(40, 55, 20, n, pole_pairs=2, name="q")   # quadrupole: (2+1)=3x
    angles = [parse_magnetization(c.label) for c in hb.children]
    assert abs(((angles[1]-angles[0]) % 360.0) - 3*(360.0/n)) < 1e-6


def test_c_core_has_pole_gap():
    W, H, D, leg, gap = 100, 80, 60, 18, 24
    core = c_core(W, H, D, leg, gap, name="yoke")
    assert core.is_valid and core.label == "yoke"
    # the gap removes material: volume < the closed O-frame (window-frame) volume
    frame_vol = W*H*D - (W-2*leg)*(H-2*leg)*D
    assert core.volume < frame_vol - 0.5*gap*leg*D, "the pole gap must actually open the frame"


def test_solenoid_is_tube():
    s = solenoid(20, 30, 40, name="winding")
    assert s.is_valid and s.label == "winding"
    assert abs(s.volume - math.pi*(30**2-20**2)*40)/s.volume < 1e-9


def test_halbach_meshes_in_netgen():
    """CAE gate: a Halbach segment tet-meshes through build123d -> STEP -> Netgen."""
    import tempfile
    from build123d import export_step
    from netgen.occ import OCCGeometry
    from ngsolve import Mesh
    hb = halbach_ring(40, 55, 20, 12, name="hb")
    seg0 = hb.children[0]
    with tempfile.TemporaryDirectory() as d:
        f = os.path.join(d, "seg.step")
        export_step(seg0, f)
        mesh = Mesh(OCCGeometry(f).GenerateMesh(maxh=6.0))
    assert mesh.ne > 50, f"Halbach segment should tet-mesh (got {mesh.ne})"


def main():
    test_magnetization_label_roundtrip()
    test_pm_primitives_carry_magnetization()
    test_halbach_ring_mallinson_easy_axes()
    test_halbach_quadrupole_advances_3x()
    test_c_core_has_pole_gap()
    test_solenoid_is_tube()
    test_halbach_meshes_in_netgen()
    print("[OK] build123d EM archetypes: PM cylinder/block, Halbach ring (Mallinson easy axes via "
          "magnetization labels), C-core yoke with pole gap, solenoid bundle -- meshable, label-driven.")


if __name__ == "__main__":
    main()
