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

import pytest

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)
from radia_mcp.build123d.archetypes import (magnetization_tag, parse_magnetization, magnetization_map,
                                            cylindrical_magnet, block_magnet, halbach_ring, c_core,
                                            solenoid, pole_tip, multipole_yoke, h_dipole, helmholtz_pair,
                                            cos_theta_dipole, e_core, slotted_stator, spm_rotor,
                                            litz_wire, litz_packing_radius, litz_fill_factor,
                                            litz_single_layer_metrics, hierarchical_litz, rectangular_litz,
                                            rectangular_litz_fill_factor, litz_serving,
                                            cos_theta_dipole_layout)
from radia_mcp.build123d.archetypes import _carried_centerline, _superposed_centerline
from radia_mcp.build123d.archetypes import involute_gear, threaded_rod, airfoil, blade
from radia_mcp.build123d.archetypes import gear_rack, bevel_gear, worm, chain_sprocket, vbelt_pulley
from radia_mcp.build123d.archetypes import ipm_rotor, squirrel_cage_rotor, claw_pole_rotor
from build123d import Box, Rectangle, RegularPolygon, extrude
import math as _math


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


def test_cos_theta_dipole_layout_matches_geometry():
    layout = cos_theta_dipole_layout(8)
    assert len(layout) == 16
    assert [row["current_sign"] for row in layout[:2]] == [1, -1]
    go = [row for row in layout if row["group"] == "go"]
    assert [row["index"] for row in go] == list(range(8))
    diffs = [go[i + 1]["sin_theta"] - go[i]["sin_theta"] for i in range(len(go) - 1)]
    assert max(diffs) == pytest.approx(min(diffs))

    ct = cos_theta_dipole(45, 4, 8, 120, n_per_half=8, name="ct")
    got = {c.label: math.degrees(math.atan2(c.center().Y, c.center().X)) % 360.0
           for c in ct.children}
    for row in layout:
        label = f"ct_{row['group']}_{row['index']:02d}"
        assert got[label] == pytest.approx(row["angle_deg"])

    with pytest.raises(ValueError):
        cos_theta_dipole_layout(1)


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


def test_litz_wire_twisted_strands():
    N, rs, Rb, L, pitch = 7, 0.5, 1.6, 30.0, 12.0
    litz = litz_wire(N, rs, Rb, L, pitch, name="litz")
    assert len(litz.children) == N and all(s.is_valid for s in litz.solids())
    turns = L / pitch
    strand_len = math.sqrt(L ** 2 + (turns * 2 * math.pi * Rb) ** 2)
    vol = sum(s.volume for s in litz.solids())
    assert abs(vol - N * math.pi * rs ** 2 * strand_len) / vol < 1e-3, "N strands x circle x helix length"
    assert [c.label for c in litz.children] == [f"litz_{k:02d}" for k in range(N)], "per-strand regions"


def test_litz_packing_radius_strands_touch():
    rs = 0.5
    for n in (6, 12, 19):
        R = litz_packing_radius(n, rs)
        sep = 2 * R * math.sin(math.pi / n)              # neighbour strand-centre distance
        assert abs(sep - 2 * rs) < 1e-9, "single-layer packing: neighbours just touch"


def test_litz_strand_meshes_in_netgen():
    """CAE gate: a Litz strand (a swept helix) tet-meshes through build123d -> STEP -> Netgen."""
    import tempfile
    from build123d import export_step
    from netgen.occ import OCCGeometry
    from ngsolve import Mesh
    litz = litz_wire(3, 0.6, 1.8, 12.0, 12.0, name="litz")
    with tempfile.TemporaryDirectory() as d:
        f = os.path.join(d, "strand.step")
        export_step(litz.children[0], f)
        mesh = Mesh(OCCGeometry(f).GenerateMesh(maxh=1.5))
    assert mesh.ne > 50, f"a Litz strand should tet-mesh (got {mesh.ne})"


def test_litz_wire_noncircular_section():
    """Non-circular strand: a rectangular (flat / edge) wire swept on the helix has volume
    rect_area x helix_length, the same twist machinery with a non-round section."""
    N, w, h, Rb, L, pitch = 6, 0.8, 0.4, 1.6, 24.0, 12.0
    litz = litz_wire(N, 0.5, Rb, L, pitch, name="flat", strand_section=Rectangle(w, h))
    assert len(litz.children) == N and all(s.is_valid for s in litz.solids())
    turns = L / pitch
    strand_len = math.sqrt(L ** 2 + (turns * 2 * math.pi * Rb) ** 2)
    vol = sum(s.volume for s in litz.solids())
    assert abs(vol - N * (w * h) * strand_len) / vol < 1e-3, "N strands x rect area x helix length"


def test_litz_fill_factor():
    """Copper fill = n*(rs/Rb)^2; single-layer PHYSICAL envelope (ring radius + rs) gives the closed form
    n sin^2(pi/n)/(1+sin(pi/n))^2, a fraction < 1 (the centre-ring radius alone would overcount > 1)."""
    assert abs(litz_fill_factor(7, 0.5, 1.5) - 7 * (0.5 / 1.5) ** 2) < 1e-12
    for n in (6, 12, 19):
        rs = 0.5
        env = litz_packing_radius(n, rs) + rs                 # physical outer envelope
        ff = litz_fill_factor(n, rs, env)
        s = math.sin(math.pi / n)
        assert abs(ff - n * s ** 2 / (1 + s) ** 2) < 1e-12, "single-layer fill via physical envelope"
        assert ff < 1.0, "physical fill factor is a fraction"


def test_litz_single_layer_metrics_round_envelope():
    """Single-layer round Litz metrics bridge the strand centre ring, bare envelope, serving envelope,
    and fill-factor report used by the CAD + solver handoff."""
    rs = 0.5
    for n in (6, 12, 19):
        m = litz_single_layer_metrics(n, rs)
        assert abs(m["center_radius"] - litz_packing_radius(n, rs)) < 1e-12
        assert abs(m["envelope_radius"] - (litz_packing_radius(n, rs) + rs)) < 1e-12
        assert abs(m["center_spacing"] - 2.0 * rs) < 1e-12
        assert abs(m["fill_factor"] - litz_fill_factor(n, rs, m["envelope_radius"])) < 1e-12

    gapped = litz_single_layer_metrics(12, rs, strand_gap=0.2, serving_thickness=0.3)
    touch = litz_single_layer_metrics(12, rs)
    assert gapped["center_spacing"] == 2.0 * rs + 0.2
    assert gapped["envelope_radius"] > touch["envelope_radius"]
    assert gapped["served_radius"] == gapped["envelope_radius"] + 0.3
    assert gapped["fill_factor"] < touch["fill_factor"]
    assert gapped["served_fill_factor"] < gapped["fill_factor"]

    for bad in (
        lambda: litz_single_layer_metrics(2, rs),
        lambda: litz_single_layer_metrics(6, 0.0),
        lambda: litz_single_layer_metrics(6, rs, strand_gap=-1e-3),
        lambda: litz_single_layer_metrics(6, rs, serving_thickness=-1e-3),
    ):
        with pytest.raises(ValueError):
            bad()


def test_litz_single_layer_metrics_match_build123d_bbox():
    """The analytic envelope radius matches the actual build123d strand bundle XY bounding box."""
    n, rs, L, pitch = 6, 0.35, 14.0, 14.0
    m = litz_single_layer_metrics(n, rs, strand_gap=0.1)
    litz = litz_wire(n, rs, m["center_radius"], L, pitch, name="round")
    bb = litz.bounding_box()
    assert abs(bb.size.X / 2.0 - m["envelope_radius"]) < 1e-6
    assert abs(bb.size.Y / 2.0 - m["envelope_radius"]) < 1e-6


def test_rectangular_litz_fill_factor():
    """Rectangular Litz fill is copper area over the tight rectangular envelope; for a touching square
    grid it tends to pi/4, and stretched pitch lowers the fill."""
    rs = 0.4
    assert abs(rectangular_litz_fill_factor(1, 1, rs, 1.0) - math.pi / 4.0) < 1e-12

    ff_touch = rectangular_litz_fill_factor(5, 4, rs, 2 * rs)
    assert abs(ff_touch - math.pi / 4.0) < 1e-12

    ff_stretched = rectangular_litz_fill_factor(5, 4, rs, 1.1, 1.4)
    expected = 5 * 4 * math.pi * rs ** 2 / (((5 - 1) * 1.1 + 2 * rs) * ((4 - 1) * 1.4 + 2 * rs))
    assert abs(ff_stretched - expected) < 1e-12
    assert ff_stretched < ff_touch

    with pytest.raises(ValueError):
        rectangular_litz_fill_factor(2, 1, rs, 0.5)


def _superposed_len(levels, indices, length, n):                 # independent re-impl for the test
    pts = []
    for s in range(n + 1):
        z = length * s / n
        x = y = 0.0
        for (count, radius, pitch), idx in zip(levels, indices):
            ph = 2 * math.pi * z / pitch + 2 * math.pi * idx / count
            x += radius * math.cos(ph); y += radius * math.sin(ph)
        pts.append((x, y, z))
    return sum(math.dist(a, b) for a, b in zip(pts, pts[1:]))


def test_hierarchical_litz_coiled_coil():
    """A 3x2 coiled-coil carries prod(count)=6 strands; each strand volume = pi rs^2 x its superposed
    centreline length (Pappus tube), and the strands are separate labelled regions."""
    levels = [(3, 2.0, 30.0), (2, 0.7, 8.0)]
    rs, L, n = 0.25, 24.0, 160
    cab = hierarchical_litz(levels, rs, L, name="hl", n_axial=n)
    import itertools
    combos = list(itertools.product(range(3), range(2)))
    assert len(cab.children) == len(combos) == 6 and all(s.is_valid for s in cab.solids())
    assert [c.label for c in cab.children] == [f"hl_{i:02d}_{j:02d}" for (i, j) in combos]
    for c, (i, j) in zip(cab.children, combos):
        exp = math.pi * rs ** 2 * _superposed_len(levels, (i, j), L, n)
        assert abs(c.volume - exp) / exp < 0.02, "strand vol = pi rs^2 x centreline length"


def test_rectangular_litz_grid_and_twist():
    """Rectangular bundle: nx*ny strands on a grid. Straight -> each volume = pi rs^2 L exactly and the
    bundle bbox spans the grid; twisted -> strands stay valid following per-strand helices."""
    nx, ny, rs, pitch, L = 2, 3, 0.4, 1.0, 20.0
    rl = rectangular_litz(nx, ny, rs, pitch, L, twist_pitch=None, name="rl")
    assert len(rl.children) == nx * ny and all(s.is_valid for s in rl.solids())
    for s in rl.solids():
        assert abs(s.volume - math.pi * rs ** 2 * L) / (math.pi * rs ** 2 * L) < 1e-6, "straight extrude"
    bb = rl.bounding_box()
    assert abs(bb.size.X - ((nx - 1) * pitch + 2 * rs)) < 1e-6, "x spans the grid"
    assert abs(bb.size.Y - ((ny - 1) * pitch + 2 * rs)) < 1e-6, "y spans the grid"
    tw = rectangular_litz(nx, ny, rs, pitch, L, twist_pitch=40.0, name="rl", n_axial=120)
    assert len(tw.children) == nx * ny and all(s.is_valid for s in tw.solids())


def test_hierarchical_litz_carried():
    """carried=True carries each inner orbit in the parent's rotation-minimizing frame, so the orbit is
    perpendicular to the local parent tangent everywhere (additive keeps it in the lab plane)."""
    import numpy as np
    levels = [(3, 2.0, 30.0), (2, 0.7, 8.0)]
    L, n, combo = 24.0, 300, (1, 1)

    def perp_cos(fn):
        strand = np.array(fn(levels, combo, L, n))
        parent = np.array(fn(levels[:-1], combo[:-1], L, n))      # sub-bundle (outer level only)
        T = np.zeros_like(parent)
        T[1:-1] = parent[2:] - parent[:-2]; T[0] = parent[1] - parent[0]; T[-1] = parent[-1] - parent[-2]
        T /= np.linalg.norm(T, axis=1, keepdims=True)
        off = strand - parent; off /= np.linalg.norm(off, axis=1, keepdims=True)
        return np.abs(np.einsum("ij,ij->i", off, T))

    assert perp_cos(_carried_centerline).max() < 1e-6, "carried orbit perpendicular to parent tangent"
    assert perp_cos(_superposed_centerline).max() > 0.1, "additive orbit is NOT perpendicular (lab frame)"
    cab = hierarchical_litz(levels, 0.25, L, name="hc", n_axial=160, carried=True)
    assert len(cab.children) == 6 and all(s.is_valid for s in cab.solids())
    assert [c.label for c in cab.children] == [f"hc_{i:02d}_{j:02d}" for i in range(3) for j in range(2)]


def test_litz_wire_insulation():
    """insulation>0 enamels each strand: copper core {name}_kk + coaxial shell {name}_kk_ins, a flat
    2N-region compound; core and shell volumes match circle / annulus x helix length."""
    N, rs, Rb, L, pitch, t = 5, 0.5, 1.6, 24.0, 12.0, 0.12
    litz = litz_wire(N, rs, Rb, L, pitch, name="enam", insulation=t)
    assert len(litz.children) == 2 * N and all(s.is_valid for s in litz.solids())
    cores = [c for c in litz.children if not c.label.endswith("_ins")]
    shells = [c for c in litz.children if c.label.endswith("_ins")]
    assert [c.label for c in cores] == [f"enam_{k:02d}" for k in range(N)]
    assert [c.label for c in shells] == [f"enam_{k:02d}_ins" for k in range(N)]
    turns = L / pitch
    hlen = math.sqrt(L ** 2 + (turns * 2 * math.pi * Rb) ** 2)
    cv = sum(c.volume for c in cores); sv = sum(s.volume for s in shells)
    assert abs(cv - N * math.pi * rs ** 2 * hlen) / cv < 1e-3, "cores = N circle x helix length"
    assert abs(sv - N * math.pi * ((rs + t) ** 2 - rs ** 2) * hlen) / sv < 1e-3, "shells = N annulus x len"


def test_litz_serving_tube():
    """Bundle serving = a concentric tube of given radial thickness around the bundle envelope."""
    env, t, L = 2.0, 0.2, 24.0
    serv = litz_serving(env, t, L, name="serve")
    assert serv.is_valid and serv.label == "serve"
    assert abs(serv.volume - math.pi * ((env + t) ** 2 - env ** 2) * L) / serv.volume < 1e-9, "annulus tube"


def test_hierarchical_litz_meshes_in_netgen():
    """CAE gate: a coiled-coil strand (superposed-helix spline sweep) tet-meshes via STEP -> Netgen."""
    import tempfile
    from build123d import export_step
    from netgen.occ import OCCGeometry
    from ngsolve import Mesh
    cab = hierarchical_litz([(3, 1.6, 24.0), (2, 0.6, 8.0)], 0.25, 16.0, name="hl", n_axial=120)
    with tempfile.TemporaryDirectory() as d:
        f = os.path.join(d, "coil.step")
        export_step(cab.children[0], f)
        mesh = Mesh(OCCGeometry(f).GenerateMesh(maxh=0.8))
    assert mesh.ne > 50, f"a coiled-coil strand should tet-mesh (got {mesh.ne})"


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


def test_involute_gear_spur_and_helical():
    spur = involute_gear(18, 4.0, 1.0, 1.2, name="spur")
    assert spur.is_valid and spur.label == "spur"
    assert spur.volume > math.pi * (4.0 * 1.06) ** 2 * 1.2, "teeth add material beyond the root cylinder"
    helical = involute_gear(16, 3.5, 0.9, 2.0, twist_deg=20.0, name="hel")
    assert helical.is_valid and helical.volume > 0


def test_threaded_rod_adds_thread():
    rod = threaded_rod(1.5, 1.0, 6.0, name="bolt")
    assert rod.is_valid and rod.label == "bolt"
    assert rod.volume > math.pi * 1.5 ** 2 * 6.0, "the V-thread adds material to the core"


def test_airfoil_and_blade():
    sec = extrude(airfoil(4.0, 0.12), 1.0)
    assert sec.is_valid and sec.volume > 0
    bl = blade([(4.0, 0.14, 0, 0), (3.0, 0.12, 3, 12), (2.0, 0.10, 6, 25)], name="vane")
    assert bl.is_valid and bl.label == "vane" and bl.volume > 0


def test_transmission_family():
    rack = gear_rack(8, 1.4, 1.2, 1.0, name="rack")
    assert rack.is_valid and rack.volume > 1.2 * 1.0 * (9 * 1.4) * 0.99, "bar + teeth"
    bevel = bevel_gear(14, 3.0, 0.8, 2.5, name="bevel")
    assert bevel.is_valid and bevel.volume > 0
    wm = worm(1.0, 2.0, 8.0, name="worm")
    assert wm.is_valid and wm.volume > _math.pi * 1.0 ** 2 * 8.0, "thread adds material"
    spr = chain_sprocket(16, 4.0, 0.8, 0.5, name="spr")
    assert spr.is_valid and spr.volume < _math.pi * 4.0 ** 2 * 0.8, "roller seats remove material"
    pul = vbelt_pulley(3.5, 2.0, 0.7, 0.7, name="pul")
    assert pul.is_valid and pul.volume < _math.pi * 3.5 ** 2 * 2.0, "groove + bore remove material"


def test_ipm_rotor_regions_and_magnetization():
    n = 4
    rotor = ipm_rotor(0.8, 4.0, n, 1.2, 0.4, 40.0, 2.0, name="pm")
    assert len(rotor.children) == 1 + 2 * n and all(s.is_valid for s in rotor.solids())
    iron = [c for c in rotor.children if c.label == "rotor_iron"]
    mags = [c for c in rotor.children if c.label != "rotor_iron"]
    assert len(iron) == 1 and len(mags) == 2 * n
    angles = sorted({round(parse_magnetization(m.label)) % 360 for m in mags})
    assert len(angles) >= 2, "magnets alternate N/S (more than one easy-axis angle)"


def test_squirrel_cage_regions():
    nb = 16
    rotor = squirrel_cage_rotor(0.8, 3.0, nb, 0.3, 0.4, 4.0, name="bar")
    assert len(rotor.children) == 1 + nb + 2 and all(s.is_valid for s in rotor.solids())
    labels = {c.label for c in rotor.children}
    assert "rotor_iron" in labels and "end_ring_hi" in labels and "end_ring_lo" in labels
    assert sum(1 for c in rotor.children if c.label.startswith("bar_")) == nb


def test_claw_pole_two_halves():
    rotor = claw_pole_rotor(0.6, 3.0, 6, 1.0, 5.0)
    assert {c.label for c in rotor.children} == {"claw_north", "claw_south"}
    assert all(s.is_valid for s in rotor.solids())


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
    test_litz_wire_twisted_strands()
    test_litz_packing_radius_strands_touch()
    test_litz_strand_meshes_in_netgen()
    test_litz_wire_noncircular_section()
    test_litz_fill_factor()
    test_litz_single_layer_metrics_round_envelope()
    test_litz_single_layer_metrics_match_build123d_bbox()
    test_hierarchical_litz_coiled_coil()
    test_rectangular_litz_grid_and_twist()
    test_hierarchical_litz_carried()
    test_litz_wire_insulation()
    test_litz_serving_tube()
    test_involute_gear_spur_and_helical()
    test_threaded_rod_adds_thread()
    test_airfoil_and_blade()
    test_transmission_family()
    test_ipm_rotor_regions_and_magnetization()
    test_squirrel_cage_regions()
    test_claw_pole_two_halves()
    test_hierarchical_litz_meshes_in_netgen()
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
