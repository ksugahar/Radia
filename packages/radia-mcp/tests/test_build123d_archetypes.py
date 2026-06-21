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
from radia_mcp.build123d.archetypes import (magnetization_tag, parse_magnetization, magnetization_map,
                                            cylindrical_magnet, block_magnet, halbach_ring, c_core,
                                            solenoid, pole_tip, multipole_yoke, h_dipole, helmholtz_pair,
                                            cos_theta_dipole, e_core, slotted_stator, spm_rotor)
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


def test_cos_theta_dipole_arcsin_spacing():
    ct = cos_theta_dipole(45, 4, 8, 120, n_per_half=8, name="ct")
    assert len(ct.children) == 16 and all(s.is_valid for s in ct.solids())
    gos = [c for c in ct.children if "_go_" in c.label]
    assert len(gos) == 8
    sins = sorted(math.sin(math.atan2(c.center().Y, c.center().X)) for c in gos)
    diffs = [sins[i+1]-sins[i] for i in range(len(sins)-1)]
    # cos-theta layout: sin(theta) is UNIFORMLY spaced (density ~ cos theta)
    assert max(diffs)/min(diffs) < 1.05, "bars must be arcsin-spaced (uniform sin theta)"


def test_e_core_two_windows():
    W, H, D, lw, bt = 120, 80, 50, 20, 25
    core = e_core(W, H, D, lw, bt, name="iron_core")
    assert core.is_valid and core.label == "iron_core"
    ww = (W - 3*lw)/2.0
    assert abs(core.volume - (W*H*D - 2*ww*(H-bt)*D)) < 1e-6, "block minus two windows"


def test_slotted_stator_removes_slots():
    s = slotted_stator(30, 55, 12, 10, 12, 40, name="stator")
    assert s.is_valid and s.label == "stator"
    full = math.pi*(55**2-30**2)*40
    assert 0 < s.volume < full, "the slots must remove material from the full ring"


def test_spm_rotor_alternating_radial_magnets():
    n = 8
    rotor = spm_rotor(10, 30, n, 6, 35, 40, name="pm")
    assert len(rotor.children) == n + 1                      # hub + n magnets
    assert any(c.label == "rotor_iron" for c in rotor.children)
    mags = [c for c in rotor.children if parse_magnetization(c.label) is not None]
    assert len(mags) == n and all(s.is_valid for s in rotor.solids())
    # radial + alternating: magnet k easy axis = k*(360/n) + (180 if k odd) (relative to its position)
    for k, c in enumerate(mags):
        expect = (k*360.0/n + (180.0 if k % 2 else 0.0)) % 360.0
        d = abs((parse_magnetization(c.label) - expect + 180.0) % 360.0 - 180.0)
        assert d < 1e-2, f"magnet {k} easy axis not radial-alternating"


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


def test_slotted_stator_meshes_in_netgen():
    """CAE gate: the slotted (toothed) stator -- a non-trivial multiply-cut solid -- tet-meshes."""
    import tempfile
    from build123d import export_step
    from netgen.occ import OCCGeometry
    from ngsolve import Mesh
    s = slotted_stator(30, 55, 12, 10, 12, 40, name="stator")
    with tempfile.TemporaryDirectory() as d:
        f = os.path.join(d, "stator.step")
        export_step(s, f)
        mesh = Mesh(OCCGeometry(f).GenerateMesh(maxh=6.0))
    assert mesh.ne > 200, f"slotted stator should tet-mesh (got {mesh.ne})"


def test_magnetization_map_follows_halbach():
    hb = halbach_ring(40, 55, 20, 12, pole_pairs=1, name="hb")
    Mmap = magnetization_map(hb, Br=1.2)
    assert len(Mmap) == 12
    for c in hb.children:
        mx, my = Mmap[c.label]
        assert abs(math.hypot(mx, my) - 1.2) < 1e-9, "|M| = Br"
        ang = math.degrees(math.atan2(my, mx))
        d = abs((ang - parse_magnetization(c.label) + 180.0) % 360.0 - 180.0)   # circular distance
        assert d < 1e-3, "M direction matches the label"
    # non-magnet regions are skipped
    assert magnetization_map([c_core(60, 50, 40, 10, 14)]) == {}


def test_pole_tip_trapezoid_volume():
    p = pole_tip(30, 16, 40, 50, name="pole")
    assert p.is_valid and p.label == "pole"
    assert abs(p.volume - 0.5*(30+16)*40*50) < 1e-6, "trapezoid area x depth"


def test_multipole_yoke_regions():
    for n in (2, 4, 6):
        y = multipole_yoke(n, 25, 20, 14, 12, 60, name="m")
        assert len(y.children) == n + 1, f"{n} poles + 1 return ring"
        assert all(s.is_valid for s in y.solids())
        labels = [c.label for c in y.children]
        assert "yoke" in labels and sum(l.startswith("pole") for l in labels) == n


def test_h_dipole_has_gap_and_poles():
    y = h_dipole(120, 100, 60, 20, 30, 24, name="dip")
    assert y.is_valid and y.label == "dip" and y.volume > 0
    # poles + frame fused into one solid that spans more than the bare frame
    frame_vol = 120*100*60 - (120-40)*(100-40)*60
    assert y.volume > frame_vol, "the protruding poles add iron beyond the window frame"


def test_helmholtz_pair_two_coils():
    hh = helmholtz_pair(40, 50, 15, 45, name="hh")
    assert len(hh.children) == 2 and all(s.is_valid for s in hh.solids())
    zc = sorted(c.center().Z for c in hh.children)
    assert abs((zc[1]-zc[0]) - 45) < 1e-6, "coils separated along z"


def main():
    test_magnetization_label_roundtrip()
    test_pm_primitives_carry_magnetization()
    test_halbach_ring_mallinson_easy_axes()
    test_halbach_quadrupole_advances_3x()
    test_c_core_has_pole_gap()
    test_solenoid_is_tube()
    test_magnetization_map_follows_halbach()
    test_pole_tip_trapezoid_volume()
    test_multipole_yoke_regions()
    test_h_dipole_has_gap_and_poles()
    test_helmholtz_pair_two_coils()
    test_cos_theta_dipole_arcsin_spacing()
    test_e_core_two_windows()
    test_slotted_stator_removes_slots()
    test_spm_rotor_alternating_radial_magnets()
    test_halbach_meshes_in_netgen()
    test_slotted_stator_meshes_in_netgen()
    print("[OK] build123d EM archetypes: PM cylinder/block, Halbach ring (Mallinson easy axes via "
          "magnetization labels) + magnetization_map, C-core / multipole / H-dipole yokes, pole tip, "
          "solenoid & Helmholtz pair -- meshable, label-driven, region-separated.")


if __name__ == "__main__":
    main()
