"""Contract for the equal-flux seed adapter on a probed Gmsh B view.

Assumes: P1 node data reproduce the linear field B = (0, 0, 1 + x/2 + y/4)
exactly, so probed seeds must equal the seeds of the analytic field.
"""

import importlib.util

import numpy as np
import pytest
from msh42_fixtures import field_section
from test_gmsh_flux_streamlines import _CUBE, _P

pytestmark = pytest.mark.skipif(importlib.util.find_spec("gmsh") is None,
                                reason="gmsh package not installed")
flux_seeding = pytest.importorskip("radia.flux_seeding")

from radia_mcp.gmsh.equal_flux import equal_flux_seeds_on_plane  # noqa: E402

ORIGIN, E_XI, E_ETA = [-0.4, -0.4, 0.0], [0.8, 0.0, 0.0], [0.0, 0.8, 0.0]


def bz(p):
    return 1 + p[..., 0] / 2 + p[..., 1] / 4


def _msh(tmp_path):
    rows = [(i, [0.0, 0.0, float(bz(np.asarray(p, float)))]) for i, p in _P.items()]
    path = tmp_path / "linear.msh"
    path.write_text(_CUBE + field_section("NodeData", "B", 3, rows), encoding="utf-8")
    return path


def test_probed_seeds_match_the_analytic_field_and_feed_the_tracer(tmp_path):
    """Guards: probe sampling, seed IDs, SI fluxes and the explicit-seed tracer."""
    result = equal_flux_seeds_on_plane(
        _msh(tmp_path), ORIGIN, E_XI, E_ETA, n_u=3, n_v=2, table_nodes=5,
        trace=True, trace_options={"step_size": 0.05, "max_steps": 40,
                                   "return_points": True})
    assert result["ok"] is True, result.get("error")

    patch = flux_seeding.planar_patch(ORIGIN, E_XI, E_ETA)
    field = lambda p: np.stack(np.broadcast_arrays(0.0, 0.0, bz(np.asarray(p))), axis=-1)
    nodes = np.linspace(0, 1, 5)
    reference = flux_seeding.equal_flux_seeds(patch.flux_density(field), xi_nodes=nodes,
                                              eta_nodes=nodes, n_u=3, n_v=2,
                                              chart=patch.chart)
    points = np.asarray([s["point_m"] for s in result["seeds"]])
    np.testing.assert_allclose(points, reference.points_m, atol=1e-12)
    assert [s["id"] for s in result["seeds"]] == list(range(6))
    assert result["orientation"] == 1
    assert result["total_flux_Wb"] == pytest.approx(0.64)
    assert result["equal_flux_quadrature"]["max_relative_cell_deviation"] < 1e-12

    trace = result["trace"]
    assert trace["ok"] is True, trace.get("error")
    assert len(trace["polylines"]) == 6
    for seed, line in zip(points, trace["polylines"]):
        xy = np.asarray(line["points"])[:, :2]
        np.testing.assert_allclose(xy, np.broadcast_to(seed[:2], xy.shape), atol=1e-9)


def test_patch_outside_the_data_fails_loudly(tmp_path):
    """Guards: samples outside the mesh are not silently filled."""
    result = equal_flux_seeds_on_plane(
        _msh(tmp_path), [0.2, 0.2, 0.0], E_XI, E_ETA, n_u=2, n_v=2, table_nodes=3)
    assert result["ok"] is False
    assert "outside the field data" in result["error"]


def _two_view_msh(tmp_path):
    """B is vertical at step 0 and tilted at step 1, beside a scalar view."""
    vertical = []
    tilted = []
    for tag, point in _P.items():
        base_bz = float(bz(np.asarray(point, float)))
        vertical.append((tag, [0.0, 0.0, base_bz]))
        tilted.append((tag, [base_bz, 0.0, 2 * base_bz]))
    other = [(i, [1.0]) for i in _P]
    path = tmp_path / "two_views.msh"
    path.write_text(_CUBE + field_section("NodeData", "B", 3, vertical, step=0)
                    + field_section("NodeData", "B", 3, tilted, time=1.0, step=1)
                    + field_section("NodeData", "other", 1, other), encoding="utf-8")
    return path


@pytest.mark.parametrize("view", [None, 0, "B"])
def test_one_view_and_step_serve_probe_and_trace(tmp_path, view):
    """Guards: implicit and named views resolve once; step 1 of B is used throughout."""
    result = equal_flux_seeds_on_plane(
        _two_view_msh(tmp_path), ORIGIN, E_XI, E_ETA, n_u=2, n_v=2, table_nodes=3,
        view=view, time_step=1, quadrature_order=0, trace=True,
        trace_options={"step_size": 0.05, "max_steps": 40, "return_points": True})
    assert result["ok"] is True, result.get("error")
    assert result["total_flux_Wb"] == pytest.approx(2 * 0.64)
    assert result["provenance"]["view"] == (0 if view is None else view)
    assert result["trace"]["ok"] is True, result["trace"].get("error")
    assert len(result["trace"]["polylines"]) == 4
    for line in result["trace"]["polylines"]:
        points = np.asarray(line["points"])
        displacement = points - points[0]
        assert np.ptp(points[:, 2]) > 0.05
        np.testing.assert_allclose(displacement[:, 0], 0.5 * displacement[:, 2],
                                   atol=1e-9)


@pytest.mark.parametrize("two_steps", [False, True])
def test_negative_time_step_is_rejected_before_probing(tmp_path, monkeypatch, two_steps):
    """Guards: the all-steps sentinel never reaches the single-step probe."""
    path = _two_view_msh(tmp_path) if two_steps else _msh(tmp_path)

    def unexpected_probe(*args, **kwargs):
        pytest.fail("negative time step reached the probe")

    monkeypatch.setattr("radia_mcp.gmsh.equal_flux.probe_field", unexpected_probe)
    result = equal_flux_seeds_on_plane(
        path, ORIGIN, E_XI, E_ETA, n_u=2, n_v=2, time_step=-1, trace=True)
    assert result == {"ok": False, "error": "time_step must be a non-negative integer"}


@pytest.mark.parametrize("key", ["view", "time_step", "seeds", "seed_start", "unknown"])
def test_reserved_and_unknown_trace_options_are_rejected_first(tmp_path, key):
    """Guards: reserved or unknown tracer keys fail before any probing."""
    missing = tmp_path / "absent_field.msh"
    missing.write_text(_CUBE, encoding="utf-8")  # no view: probing would fail differently
    result = equal_flux_seeds_on_plane(
        missing, ORIGIN, E_XI, E_ETA, n_u=2, n_v=2, trace=True, trace_options={key: 0})
    assert result["ok"] is False
    assert "unsupported trace_options" in result["error"] and key in result["error"]
