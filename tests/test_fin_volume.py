"""Contracts of the experimental unreduced volume-lattice fin PEEC."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

FIXTURES = Path(__file__).resolve().parent / "coil_from_cad" / "fixtures"


def _shoelace(uv):
    x, y = uv[:, 0], uv[:, 1]
    return 0.5 * float(np.dot(x, np.roll(y, -1)) - np.dot(np.roll(x, -1), y))


@pytest.fixture(scope="module")
def short_fin():
    from radia.fin_volume import solve_fin_volume

    return solve_fin_volume(FIXTURES / "beak_fin_short.step", frequency=222e3,
                            sigma=5.8e7, current_A=2.0, n_layers=2, n_lanes=16,
                            n_stations=7, lane_grading="uniform")


def test_volume_lattice_topology_and_tiling(short_fin):
    """Guards: the lattice is one connected graph whose longitudinal bars
    tile each section (dual areas sum to the lane polygon area, computed
    independently by the shoelace formula) and whose branch count follows
    lattice_branch_count.
    Assumes: beak_fin_short is star-shaped about each section centroid."""
    from radia.fin_volume import LONGITUDINAL, lattice_branch_count

    lat, sweep = short_fin["lattice"], short_fin["sweep"]
    gates = lat.gates
    assert gates["n_branches"] == lattice_branch_count(7, 16, 2) == len(lat.branches)
    assert gates["connected_components"] == 1
    assert gates["incidence_rank"] == gates["n_nodes"] - 1
    assert gates["min_signed_cell_area_m2"] > 0

    # Interior longitudinal bar areas between stations 2 and 3 average the
    # two section tilings; compare with the two lane-polygon areas.
    lanes = [np.asarray(f.lane_uv) for f in sweep.features]
    expected = 0.5 * (_shoelace(lanes[2]) + _shoelace(lanes[3]))
    sel = (lat.kind == LONGITUDINAL) & (lat.station == 2)
    got = float(np.sum(lat.width[sel] * lat.height[sel]))
    assert abs(got - expected) <= 1e-12 * expected, (got, expected)


def test_volume_currents_balance_and_residuals(short_fin):
    """Guards: the full port current crosses every station gap, and the
    solve meets its own KCL/KVL true-residual limit.
    Assumes: one port between the first and last station."""
    from radia.fin_volume import LONGITUDINAL

    lat, current = short_fin["lattice"], short_fin["current"]
    for s in range(6):
        sel = (lat.kind == LONGITUDINAL) & (lat.station == s)
        assert abs(np.sum(current[sel]) - 2.0) < 1e-9, (s, np.sum(current[sel]))
    res = short_fin["coil_true_residual"]
    assert res["kcl_relative"] < 1e-9 and res["kvl_relative"] < 1e-9, res
    assert short_fin["R_ohm"] > 0 and short_fin["L_H"] > 0


def test_volume_lattice_rejects_points_outside_the_section():
    """Guards: a lattice point outside its section raises instead of
    meshing air (the radial map assumes a star-shaped section).
    Assumes: shrinking a station outline to half scale exposes the
    outer lattice ring."""
    import dataclasses

    from radia.fin_sweep import fin_graph_from_step
    from radia.fin_volume import build_volume_lattice

    sweep = fin_graph_from_step(FIXTURES / "beak_fin_short.step", n_lanes=16,
                                n_stations=7, lane_grading="uniform")
    stations = list(sweep.stations)
    stations[3] = dataclasses.replace(stations[3], uv=0.5 * stations[3].uv)
    bad = dataclasses.replace(sweep, stations=tuple(stations))
    with pytest.raises(ValueError, match="outside their section"):
        build_volume_lattice(bad, 3)


def test_branch_budget_refuses_actual_stations_before_the_lattice(monkeypatch):
    """Guards: the dense-kernel budget uses the stations actually placed
    (feature-aligned placement can add some), before any lattice is built.
    Assumes: lattice_branch_count(29, 24, 6) = 11836 <= 12000 < 12702 (31)."""
    from types import SimpleNamespace
    import radia.fin_sweep as fsw
    import radia.fin_volume as fv

    def graph(*a, n_stations, **k):
        return SimpleNamespace(meta={"n_stations": graph.actual})

    def lattice(*a, **k):
        raise AssertionError("lattice built")

    monkeypatch.setattr(fsw, "fin_graph_from_step", graph)
    monkeypatch.setattr(fv, "build_volume_lattice", lattice)
    kw = dict(sigma=5.8e7, n_layers=6, n_lanes=24, n_stations=29)
    graph.actual = 31
    with pytest.raises(ValueError, match="31 stations after feature-aligned"):
        fv.assemble_fin_volume("x.step", **kw)
    graph.actual = 29
    with pytest.raises(AssertionError, match="lattice built"):
        fv.assemble_fin_volume("x.step", **kw)
