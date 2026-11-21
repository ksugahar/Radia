"""Physical lane correspondence of the cap-march station pipeline (no CAD).

Sections are generated in a twist-free reference frame built independently
(SciPy minimal rotations between successive tangents), sampled from random
wire seams, and pushed through the production pipeline: frames, outline
origin, registration and lattice construction.  Each lane point is then read
back in the reference section coordinates; a lane that follows one material
line keeps the same coordinates along the sweep.
"""

from __future__ import annotations

import math

import numpy as np
import pytest
from scipy.spatial.transform import Rotation

from radia.fin_cap_march import MarchError, outline_origin, section_frames
from radia.fin_section import resample_outline
from radia.fin_sweep import (
    StationOutline,
    align_station_origins,
    build_fin_graph,
    cap_march_station_outlines,
)

N_OUT = 1024


def _rect(a=3.0e-3, b=1.0e-3, n=N_OUT, phase=0.0):
    per = 4 * (a + b)
    q = (np.arange(n) + phase) / n * per
    c = np.array([[a, -b], [a, b], [-a, b], [-a, -b], [a, -b]])
    s = np.r_[0.0, np.cumsum(np.linalg.norm(np.diff(c, axis=0), axis=1))]
    return np.column_stack([np.interp(q, s, c[:, k]) for k in range(2)])


def _rounded_rect(a=3.0e-3, b=2.0e-3, r=0.4e-3, n=N_OUT):
    t = np.linspace(0.0, 2 * np.pi, 8192, endpoint=False)
    c, s = np.cos(t), np.sin(t)
    p = np.column_stack([np.sign(c) * (a - r), np.sign(s) * (b - r)]) + r * np.column_stack([c, s])
    return resample_outline(p, n)


def _beak(n=N_OUT):
    """Rectangle with a round-tipped fin along +a (the auto-grading path)."""
    tl = (3.8, -(0.25**2 - 0.05**2) ** 0.5)
    tu = (3.8, +(0.25**2 - 0.05**2) ** 0.5)
    a0 = math.atan2(tl[1], tl[0] - 3.75)
    a1 = math.atan2(tu[1], tu[0] - 3.75)
    pts = [(-4, -2), (1.6, -2), (1.6, -0.7), tl]
    pts += [(3.75 + 0.25 * math.cos(x), 0.25 * math.sin(x))
            for x in np.linspace(a0, a1, 48)[1:-1]]
    pts += [tu, (1.6, 0.7), (1.6, 2), (-4, 2)]
    return resample_outline(np.array(pts) * 1e-3, n)


def _path(kind, n):
    if kind == "straight":
        c = np.c_[np.zeros(n), np.zeros(n), np.linspace(0, 0.03, n)]
        return c, np.tile([0.0, 0.0, 1.0], (n, 1))
    if kind == "planar_arc":
        th = np.linspace(0, 1.5 * np.pi, n)
        return (np.c_[0.02 * np.cos(th), 0.02 * np.sin(th), np.zeros(n)],
                np.c_[-np.sin(th), np.cos(th), np.zeros(n)])
    if kind == "helix":
        th = np.linspace(0, 2.5 * np.pi, n)
        c = np.c_[0.02 * np.cos(th), 0.02 * np.sin(th), 0.004 * th]
        t = np.c_[-0.02 * np.sin(th), 0.02 * np.cos(th), 0.004 * np.ones(n)]
        return c, t / np.linalg.norm(t, axis=1)[:, None]
    if kind == "arc_then_vertical_lead":
        # 270-degree ring in z = 0, then a vertical lead: the centroids'
        # least-variance direction is z, parallel to the lead tangent.
        m = n - 4
        th = np.linspace(0, 1.5 * np.pi, m)
        c = np.c_[0.02 * np.cos(th), 0.02 * np.sin(th), np.zeros(m)]
        t = np.c_[-np.sin(th), np.cos(th), np.zeros(m)]
        end = c[-1]
        lead = np.array([end + [0.0, 0.0, 0.002 * (k + 1)] for k in range(4)])
        return np.vstack([c, lead]), np.vstack([t, np.tile([0.0, 0.0, 1.0], (4, 1))])
    raise ValueError(kind)


def _reference_frames(t):
    """Twist-free (U, V) along tangents, by SciPy minimal rotations."""
    seed = np.array([0.3, -0.4, 0.86])
    U = [seed - (seed @ t[0]) * t[0]]
    U[0] /= np.linalg.norm(U[0])
    for i in range(1, len(t)):
        rot, _ = Rotation.align_vectors([t[i]], [t[i - 1]])
        ui = rot.apply(U[-1])
        ui -= (ui @ t[i]) * t[i]
        U.append(ui / np.linalg.norm(ui))
    U = np.array(U)
    return U, np.cross(t, U)


def _pipeline(kind, section, n_st=12, n_lanes=24, grading="uniform", seed=0):
    rng = np.random.default_rng(seed)
    c, t = _path(kind, n_st)
    U, V = _reference_frames(t)
    xyz = []
    for i in range(n_st):
        ab = np.roll(section, int(rng.integers(0, N_OUT)), axis=0)
        if rng.random() < 0.5:
            ab = ab[::-1]  # clockwise wire
        xyz.append(c[i] + ab[:, :1] * U[i] + ab[:, 1:2] * V[i])
    area = 1.0
    per = float(np.sum(np.linalg.norm(np.roll(section, -1, axis=0) - section, axis=1)))
    stations = cap_march_station_outlines(xyz, c, t, [area] * n_st, [per] * n_st)
    stations = align_station_origins(stations)
    g = build_fin_graph(stations, n_lanes=n_lanes, lane_grading=grading,
                        lane_phase="anchored")
    d = g.rings - c[:, None, :]
    ab_lanes = np.stack([np.einsum("ikj,ij->ik", d, U), np.einsum("ikj,ij->ik", d, V)], axis=-1)
    return g, ab_lanes, per / n_lanes


@pytest.mark.parametrize("kind", ["straight", "planar_arc", "helix", "arc_then_vertical_lead"])
@pytest.mark.parametrize("shape", ["rect", "rounded"])
def test_uniform_lanes_follow_material_lines(kind, shape):
    """Guards: lane k reads back the same reference section coordinates at
    every station (frames, origin, registration and lattice together), on
    planar, straight, helical and lead-parallel-to-seed paths, for random
    wire seams and orientations."""
    section = _rect() if shape == "rect" else _rounded_rect()
    _g, ab, h = _pipeline(kind, section)
    drift = np.linalg.norm(ab - ab[:1], axis=-1).max()
    assert drift <= 0.05 * h, (drift / h, kind, shape)


@pytest.mark.parametrize("kind", ["straight", "helix", "arc_then_vertical_lead"])
def test_graded_fin_lanes_follow_material_lines(kind):
    """Guards: the automatic fin-graded lanes keep physical correspondence."""
    g, ab, h = _pipeline(kind, _beak(), grading="auto", n_lanes=32)
    assert all(f.graded for f in g.features)
    drift = np.linalg.norm(ab - ab[:1], axis=-1).max()
    assert drift <= 0.05 * h, drift / h


def test_registration_ignores_independently_rotated_local_frames():
    """Guards: the R4 review witness (geometry-probe): stations whose local
    in-plane frames are rotated independently still register to the same
    physical sample; comparing raw (u, v) arrays shifted a 3:1 rectangle by a
    quarter outline at 90 degrees."""
    base = _rect()
    o = outline_origin(base, np.array([1.0, 0.0]))
    base = np.roll(base, -o, axis=0)
    st = []
    for i, th in enumerate((0.0, 0.57, np.pi / 2, np.pi, 4.0)):
        c, s = math.cos(th), math.sin(th)
        u = np.array([c, s, 0.0])
        v = np.array([-s, c, 0.0])
        local = np.roll(base @ np.array([[c, -s], [s, c]]), 37 * i, axis=0)  # (a,b) -> (u,v)
        st.append(StationOutline(i, np.array([0, 0, 0.01 * i]), u, v,
                                 np.array([0, 0, 1.0]), local, 1.0, 1.0))
    out = align_station_origins(st)
    first = np.array([s.to_xyz(s.uv[:1])[0] - s.centroid for s in out])
    assert np.abs(first - first[0]).max() < 1e-12, first


@pytest.mark.parametrize("anchor", [(0.0, 1.0), (1.0, 0.0), (0.6, 0.8)])
def test_corner_ties_are_seam_and_frame_invariant(anchor):
    """Guards: exact mirror-symmetric corner ties (anchor normal to a flat
    side) pick one physical corner for every cyclic roll of the wire seam,
    every sub-sample seam offset and every frame rotation (the 1e-9 tie with
    sample-order argmax picked seam-dependent corners)."""
    ref = np.array(anchor)
    picks = []
    for phase in (0.0, 0.31, 0.77):
        uv = _rect(phase=phase)
        for roll in (0, 128, 256, 384, 77):
            for th in (0.0, np.pi / 2, 1.1):
                c, s = math.cos(th), math.sin(th)
                rot = np.array([[c, -s], [s, c]])
                x = np.roll(uv, roll, axis=0) @ rot.T
                k = outline_origin(x, ref @ rot.T)
                picks.append(x[k] @ rot)  # back to the physical frame
    picks = np.array(picks)
    spacing = 8 * 4e-3 / N_OUT
    assert np.abs(picks - picks[0]).max() <= 1.01 * spacing, picks[:3]
    assert np.allclose(np.abs(picks[0]), [3e-3, 1e-3], atol=spacing)


def test_square_corner_tie_left_open_by_both_directions_raises():
    """Guards: a tie the two reference directions cannot split fails loudly
    (a square read along a diagonal: two corners tie in reach and reach2 is
    symmetric about the diagonal only if they coincide; a degenerate equal
    pair is constructed directly)."""
    sq = _rect(a=1e-3, b=1e-3)
    ok = outline_origin(sq, np.array([1.0, 1.0]) / math.sqrt(2))
    assert np.allclose(sq[ok], [1e-3, 1e-3], atol=8e-3 / N_OUT)
    doubled = np.vstack([sq, sq])  # two coincident copies of each corner
    with pytest.raises(MarchError, match="ambiguous"):
        outline_origin(doubled, np.array([1.0, 0.0]))


@pytest.mark.parametrize("kind", ["straight", "planar_arc", "helix", "arc_then_vertical_lead"])
def test_section_frames_stay_transverse_and_twist_free(kind):
    """Guards: frames are unit, transverse, and match an independent
    minimal-rotation transport up to one constant in-plane rotation (no
    planarity threshold, no zero projection at a lead parallel to the
    centroids' least-variance direction)."""
    c, t = _path(kind, 14)
    u, v = section_frames(c, t)
    assert np.allclose(np.linalg.norm(u, axis=1), 1.0)
    assert np.abs(np.einsum("ij,ij->i", u, t)).max() < 1e-12
    U, V = _reference_frames(t)
    ang = np.arctan2(np.einsum("ij,ij->i", u, V), np.einsum("ij,ij->i", u, U))
    assert np.ptp(np.unwrap(ang)) < 1e-9, np.unwrap(ang)


def test_sweep_reversal_between_stations_raises():
    """Guards: an undefined transport (tangent reverses) fails loudly."""
    c = np.c_[np.zeros(3), np.zeros(3), [0.0, 0.01, 0.02]]
    t = np.array([[0, 0, 1.0], [0, 0, 1.0], [0, 0, -1.0]])
    with pytest.raises(MarchError, match="reverses"):
        section_frames(c, t)
