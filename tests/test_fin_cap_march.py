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
