"""Contracts of the terminal-to-terminal section march (radia.fin_cap_march)."""
from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest


def _straight_line():
    from radia.fin_cap_march import Centerline

    pts = np.c_[np.linspace(0.0, 10.0, 11), np.zeros(11), np.zeros(11)]
    cap = SimpleNamespace(face="cap", centroid=pts[0], inward=np.array([1.0, 0, 0]),
                          area=1.0)
    return Centerline(pts, np.tile([1.0, 0.0, 0.0], (11, 1)), np.ones(11),
                      cap, cap, 0)


def _stub(monkeypatch, good_when_tilted):
    import radia.fin_cap_march as fcm

    def components(solid, point, normal):
        tilted = abs(normal[0] - 1.0) > 1e-12
        if tilted and good_when_tilted:
            return [("good", np.asarray(point, float), 1.05)]
        # One merged face across the conductor: 5x the march area, far centroid.
        return [("merged", np.asarray(point, float) + [0.0, 6.0, 0.0], 5.0)]

    monkeypatch.setattr(fcm, "section_components", components)
    monkeypatch.setattr(fcm, "_containing", lambda comps, p, caps=(): comps[0])


def test_station_section_rejects_a_merged_face_and_retries_tilted(monkeypatch):
    """Guards: a station section inconsistent with the dense march (area
    outside its band, centroid far from the centreline point) is not used;
    a slightly tilted cut that matches the march is.
    Assumes: tilts up to 1e-3 rad are geometrically negligible."""
    from radia.fin_cap_march import station_sections

    _stub(monkeypatch, good_when_tilted=True)
    out = station_sections(None, _straight_line(), 5)
    assert [o[0] for o in out[1:-1]] == ["good"] * 3, out
    for _, _, n, _ in out[1:-1]:
        assert 0 < np.degrees(np.arccos(min(1.0, n[0]))) < 0.1


def test_station_section_raises_when_no_cut_matches_the_march(monkeypatch):
    """Guards: no silent substitution -- a station whose every cut is the
    merged face raises MarchError naming the station.
    Assumes: the stub returns only the merged face."""
    from radia.fin_cap_march import MarchError, station_sections

    _stub(monkeypatch, good_when_tilted=False)
    with pytest.raises(MarchError, match="station 1: no section"):
        station_sections(None, _straight_line(), 5)


def _fin_line(start, end, n=401):
    """Straight 0..100 sweep, area 1 except 1.1 on [start, end] (a fin)."""
    from radia.fin_cap_march import Centerline

    x = np.linspace(0.0, 100.0, n)
    pts = np.c_[x, np.zeros(n), np.zeros(n)]
    cap = SimpleNamespace(face="cap", centroid=pts[0], inward=np.array([1.0, 0, 0]),
                          area=1.0)
    areas = np.where((x >= start) & (x <= end), 1.1, 1.0)
    return Centerline(pts, np.tile([1.0, 0.0, 0.0], (n, 1)), areas, cap, cap, 0)


@pytest.mark.parametrize("n", [16, 20, 23, 25, 27, 29])
def test_station_targets_bracket_each_fin_end(n):
    """Guards: every fin end gets a station pair +/- h/16 within one dense step,
    so the captured fin length no longer depends on the station count.
    Assumes: feature ends are located on the dense centerline areas."""
    from radia.fin_cap_march import station_targets

    t = station_targets(_fin_line(31.4, 68.9), n)
    h = 100.0 / (n - 1)
    assert t[0] == 0.0 and t[-1] == 100.0 and np.all(np.diff(t) > 0), t
    for e in (31.4, 68.9):
        lo, hi = t[t < e].max(), t[t > e].min()
        assert hi - lo <= h / 8 + 0.25 + 1e-9, (n, e, lo, hi)
    assert np.min(np.diff(t)) >= h / 8 - 0.25 - 1e-9, np.min(np.diff(t))


def test_station_targets_without_features_are_equal_arc():
    from radia.fin_cap_march import station_targets

    np.testing.assert_array_equal(station_targets(_fin_line(200, 300), 7),
                                  np.linspace(0.0, 100.0, 7))


@pytest.mark.parametrize("n", [16, 23, 29])
def test_station_targets_have_no_half_pair_or_near_cap_station(n):
    """Guards: a feature end within h/3 of a cap gets no pair, so no one-sided
    pair and no station nearly coincident with a cap; spacing >= 2h/16."""
    from radia.fin_cap_march import station_targets

    h = 100.0 / (n - 1)
    t = station_targets(_fin_line(h / 6, 60.3), n)
    assert np.min(np.diff(t)) >= h / 8 * (1 - 1e-9), (n, np.min(np.diff(t)))
    assert t[1] > h / 4 and t[-2] < 100.0 - h / 4, t[:3]
    lo, hi = t[t < 60.3].max(), t[t > 60.3].min()
    assert hi - lo <= h / 8 + 0.25 + 1e-9, (lo, hi)


def test_station_targets_reject_features_closer_than_two_thirds_spacing():
    from radia.fin_cap_march import MarchError, station_targets

    with pytest.raises(MarchError, match="closer than 2h/3"):
        station_targets(_fin_line(40.0, 41.0), 16)


def _rounded_rect(n=1024, a=3.0e-3, b=2.0e-3, r=0.4e-3, rot=0.0):
    t = np.linspace(0.0, 2 * np.pi, 4096, endpoint=False)
    c, s = np.cos(t), np.sin(t)
    p = np.column_stack([np.sign(c) * (a - r), np.sign(s) * (b - r)]) + r * np.column_stack([c, s])
    seg = np.linalg.norm(np.diff(np.vstack([p, p[:1]]), axis=0), axis=1)
    cum = np.concatenate([[0.0], np.cumsum(seg)])
    q = np.linspace(0.0, cum[-1], n, endpoint=False)
    uv = np.column_stack([np.interp(q, cum, np.append(p[:, k], p[0, k])) for k in range(2)])
    k0 = np.cos(rot), np.sin(rot)
    return uv @ np.array([[k0[0], k0[1]], [-k0[1], k0[0]]])


@pytest.mark.parametrize("rot", [0.0, 0.3, 1.1, 2.5])
def test_outline_origin_is_a_corner_in_any_frame(rot):
    """Guards: the lane origin is a section corner whatever in-plane frame the
    station got (a frame-tied origin let uniform lanes skip every corner, R -7 %)."""
    from radia.fin_cap_march import outline_corner_start

    uv = _rounded_rect(rot=rot)
    p = uv[outline_corner_start(uv)]
    radius = np.linalg.norm(uv, axis=1).max()
    assert np.linalg.norm(p) >= 0.995 * radius  # arc ends sit at 0.985, (np.linalg.norm(p), radius)
