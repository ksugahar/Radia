"""Conservative transition tests for the experimental fin PEEC graph."""

import sys
import types
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from radia.peec_fin_topology import (
    assemble_experimental_fin_peec,
    build_hybrid_surface_topology,
    build_hybrid_surface_topology_from_step,
    rings_from_filament_paths,
)


def _rings(n_stations=5, n_lanes=4):
    theta = 2 * np.pi * np.arange(n_lanes) / n_lanes
    section = np.c_[np.cos(theta), np.sin(theta)]
    return np.array([
        np.c_[section, np.full(n_lanes, float(s))]
        for s in range(n_stations)
    ])


def _beak_rings():
    # Off-centre beak point is intentionally retained as its own CAD lane.
    section = np.array([
        [-1.0, -1.0], [0.8, -1.0], [1.1, -0.25],
        [2.4, 0.0], [1.1, 0.25], [0.8, 1.0], [-1.0, 1.0],
    ]) * 1e-3
    return np.array([
        np.c_[section, np.full(len(section), s * 1e-3)]
        for s in range(5)
    ])


def test_no_mesh_reduces_to_independent_series_lanes():
    graph = build_hybrid_surface_topology(_rings(), ())
    assert graph.branches.shape == (16, 2)
    currents, residual = graph.solve_reference(
        0.0, np.ones(16), np.zeros((16, 16)), 4.0)
    np.testing.assert_allclose(currents, 1.0, atol=1e-12)
    np.testing.assert_allclose(residual, 0.0, atol=1e-12)


def test_mesh_transition_is_not_an_equipotential_short():
    graph = build_hybrid_surface_topology(_rings(), (2,))
    assert graph.branches.shape == (20, 2)
    # Four distinct surface nodes at the transition, not one tied node.
    assert len(set(graph.branches[4:8, 1])) == 4
    # With an asymmetric induced EMF, transverse branches carry current,
    # but every transition node still obeys KCL.
    emf = np.zeros(20)
    emf[4] = 1.0
    currents, residual = graph.solve_reference(
        0.0, np.ones(20), np.zeros((20, 20)), 0.0, emf)
    assert np.max(np.abs(currents[16:])) > 1e-3
    np.testing.assert_allclose(residual, 0.0, atol=1e-12)
    assert abs(currents[:4].sum()) < 1e-12


def test_selected_mesh_region_only_adds_local_transverse_branches():
    graph = build_hybrid_surface_topology(_rings(), (1, 2, 3))
    assert graph.branch_kind.count("longitudinal") == 16
    assert graph.branch_kind.count("transverse") == 12
    assert graph.mesh_stations == (1, 2, 3)
    a = graph.incidence()
    np.testing.assert_allclose(np.asarray(a.sum(axis=0)), 0.0)


def test_cad_paths_recover_stations_and_mesh_unknown_section(monkeypatch):
    rings = _rings()
    paths = [[(rings[s, k], rings[s + 1, k])
              for s in range(len(rings) - 1)]
             for k in range(rings.shape[1])]
    np.testing.assert_allclose(rings_from_filament_paths(paths), rings)
    fake = types.ModuleType("radia.coil_from_cad")
    fake.filaments_from_step = lambda *args, **kwargs: {
        "source": "step_per_station", "cross_section_kind": "unknown",
        "filament_paths": paths,
    }
    monkeypatch.setitem(sys.modules, "radia.coil_from_cad", fake)
    graph, info = build_hybrid_surface_topology_from_step(
        "fin.step", sigma=5.8e7, n_peri=4)
    assert graph.mesh_stations == (1, 2, 3)
    assert info["cad_source"] == "step_per_station"


def test_cad_equivalent_circle_fallback_is_rejected(monkeypatch):
    fake = types.ModuleType("radia.coil_from_cad")
    fake.filaments_from_step = lambda *args, **kwargs: {
        "source": "step_longest_edge", "cross_section_kind": "equivalent-circle",
        "filament_paths": [],
    }
    monkeypatch.setitem(sys.modules, "radia.coil_from_cad", fake)
    with pytest.raises(ValueError, match="direct CAD perimeter"):
        build_hybrid_surface_topology_from_step(
            "fin.step", sigma=5.8e7, n_peri=4)


def test_experimental_cpp_peec_keeps_physical_geometry_and_kcl():
    rings = _rings() * 1e-3
    graph = build_hybrid_surface_topology(rings, (2,))
    # The port is one logical node, but the four physical starting points
    # remain distinct in the partial-inductance geometry.
    assert len({tuple(p) for p in graph.branch_xyz[:4, 0]}) == 4
    solver, widths = assemble_experimental_fin_peec(
        graph, rings, sigma=5.8e7, sheet_depth=0.1e-3)
    assert widths.shape == (20,)
    currents = np.asarray(solver.compute_branch_currents(1e3, [1.0]))
    assert currents.shape == (20,)
    rhs = np.zeros(len(graph.nodes), dtype=complex)
    rhs[graph.port_plus] = -1.0
    rhs[graph.port_minus] = 1.0
    np.testing.assert_allclose(graph.incidence() @ currents, rhs,
                               atol=1e-7)


def test_beak_polygon_keeps_tip_lane_and_transverse_redistribution():
    rings = _beak_rings()
    graph = build_hybrid_surface_topology(rings, (1, 2, 3))
    assert graph.n_lanes == 7
    np.testing.assert_allclose(graph.branch_xyz[3, 0], rings[0, 3])
    n = len(graph.branches)
    emf = np.zeros(n)
    emf[graph.n_lanes + 3] = 0.2
    currents, residual = graph.solve_reference(
        0.0, np.ones(n), np.zeros((n, n)), 1.0, emf)
    assert np.linalg.norm(currents[graph.branch_kind.count("longitudinal"):]) > 1e-3
    np.testing.assert_allclose(residual, 0.0, atol=1e-12)


@pytest.mark.parametrize("stations", [(-1,), (0,), (4,), (5,)])
def test_no_mesh_at_terminals_or_outside(stations):
    with pytest.raises(ValueError, match="interior"):
        build_hybrid_surface_topology(_rings(), stations)
