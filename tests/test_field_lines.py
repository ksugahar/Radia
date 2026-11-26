"""Surrogate contracts for field-agnostic field-line tracing (radia.field_lines).

Assumes analytic fields with exact answers:
- B = (-y, x, 0): circles about z through every seed;
- the straight wire B = mu0 I / (2 pi r) phi_hat: circles of constant r, and
  the plane y = 0, x in [r1, r2], z in [-h/2, h/2] carries
  mu0 I h ln(r2 / r1) / (2 pi);
- B = z_hat inside the ball |p| < R and no data (NaN) outside it.
"""

import numpy as np
import pytest

from radia.field_lines import (flux_account, surface_crossings, trace_field_lines,
                               write_field_lines_msh)
from radia.flux_seeding import balanced_equal_flux_seeds, planar_patch

MU0 = 4e-7 * np.pi


def _flat(field):
    """Accept points of any leading shape (..., 3), as seeding and tracing pass."""
    def wrapped(p):
        p = np.asarray(p, dtype=float)
        return field(p.reshape(-1, 3)).reshape(p.shape)
    return wrapped


@_flat
def rotation(p):
    return np.column_stack((-p[:, 1], p[:, 0], np.zeros(len(p))))


def wire(current=10.0):
    @_flat
    def field(p):
        r2 = p[:, 0] ** 2 + p[:, 1] ** 2
        scale = MU0 * current / (2 * np.pi * r2)
        return np.column_stack((-p[:, 1] * scale, p[:, 0] * scale, np.zeros(len(p))))
    return field


def ball(radius):
    @_flat
    def field(p):
        out = np.tile([0.0, 0.0, 1.0], (len(p), 1))
        out[np.linalg.norm(p, axis=1) >= radius] = np.nan
        return out
    return field


def test_equal_flux_wire_lines_are_closed_circles_traced_in_batches():
    """Guards: batched RK4, closure, +B orientation and the call count per step."""
    r1, r2, h, current = 0.01, 0.05, 0.04, 10.0
    patch = planar_patch([r1, 0, -h / 2], [r2 - r1, 0, 0], [0, 0, h])  # normal -y
    calls = []

    def counted(p):
        calls.append(len(p))
        return wire(current)(p)

    g = patch.flux_density(wire(current))
    nodes = np.linspace(0, 1, 33)
    seeds = balanced_equal_flux_seeds(g, xi_nodes=nodes, eta_nodes=nodes, n_lines=6,
                                      chart=patch.chart, metric=patch.metric)
    assert seeds.total_flux_Wb == pytest.approx(-MU0 * current * h * np.log(r2 / r1) / (2 * np.pi),
                                                rel=2e-3)
    lines, record = trace_field_lines(counted, seeds.points_m, step_m=0.002, max_steps=1000)
    assert len(lines) == 6 and all(line.closed for line in lines)
    longest = max(len(line.points_m) for line in lines)
    # One batched call per RK stage for all lines, not one per line.
    assert record["field_calls"] <= 4 * 2 * longest + 2, (record, longest)
    for line in lines:
        radius = np.hypot(line.points_m[:, 0], line.points_m[:, 1])
        np.testing.assert_allclose(radius, radius[0], rtol=1e-5)  # RK4 drift, Euler ~1e-2
        direction = np.diff(line.points_m[:-1], axis=0)
        b = wire(current)(line.points_m[:-2])
        assert np.all(np.sum(direction * b, axis=1) > 0)  # points run along +B


def test_rotation_field_returns_once_with_the_opposite_sign():
    """Guards: crossing signs, seed-crossing exclusion and the opposite-flux share."""
    patch = planar_patch([-0.4, 0, -0.2], [0.8, 0, 0], [0, 0, 0.4])
    g = patch.flux_density(rotation)
    nodes = np.linspace(0, 1, 5)
    seeds = balanced_equal_flux_seeds(g, xi_nodes=nodes, eta_nodes=nodes, n_lines=4,
                                      chart=patch.chart, metric=patch.metric, sign="positive")
    lines, _ = trace_field_lines(rotation, seeds.points_m, step_m=0.02, max_steps=400)
    account, warnings = flux_account(lines, patch, n_seeds=4,
                                     flux_per_line_Wb=seeds.target_cell_flux_Wb,
                                     orientation=seeds.orientation,
                                     opposite_flux_Wb=seeds.negative_flux_Wb)
    assert account["lines_by_end"]["closed"] == 4
    assert (account["repeat_crossings"], account["opposite_crossings"]) == (0, 4)
    assert account["opposite_crossing_flux_ratio"] == pytest.approx(1.0, abs=1e-12)
    assert warnings == []
    for line in lines:
        assert surface_crossings(line, patch) == [-1]


def test_lines_end_within_the_minimum_step_of_the_data_and_the_extent():
    """Guards: step halving at the data boundary and at the extent box."""
    radius, step = 0.1, 0.02
    seeds = [[0.0, 0.0, 0.0], [0.03, 0.02, 0.0]]
    lines, record = trace_field_lines(ball(radius), seeds, step_m=step, max_steps=100)
    for line in lines:
        assert line.reasons == {"forward": "left_data", "backward": "left_data"}
        ends = np.linalg.norm(line.points_m[[0, -1]], axis=1)
        assert np.all(ends < radius) and np.all(radius - ends <= record["min_step_m"] * 1.01)
    box = {"min": [-1, -1, -0.05], "max": [1, 1, 0.07]}
    lines, record = trace_field_lines(ball(radius), seeds, step_m=step, max_steps=100, extent_m=box)
    for line in lines:
        assert line.reasons == {"forward": "left_extent", "backward": "left_extent"}
        z = line.points_m[:, 2]
        tol = record["min_step_m"] * 1.01
        assert 0.07 - z.max() <= tol and z.min() + 0.05 <= tol


def test_truncation_one_sided_tracing_and_untraced_seeds_are_warned():
    patch = planar_patch([-0.4, 0, -0.2], [0.8, 0, 0], [0, 0, 0.4])
    seeds = [[-0.2, 0, 0], [5.0, 0, 0]]  # the second seed lies outside the data
    field = lambda p: np.where(np.linalg.norm(p, axis=1)[:, None] < 1, rotation(p), np.nan)
    lines, _ = trace_field_lines(field, seeds, step_m=0.02, max_steps=5, both_directions=False)
    account, warnings = flux_account(lines, patch, n_seeds=2, flux_per_line_Wb=0.01,
                                     orientation=1, opposite_flux_Wb=-0.02)
    assert account["lines_by_end"]["truncated"] == 1 and account["not_traced_lines"] == 1
    joined = " ".join(warnings)
    assert "max_steps" in joined and "no line" in joined and "forward only" in joined


def test_closed_lines_are_written_without_duplicate_nodes(tmp_path):
    lines, _ = trace_field_lines(rotation, [[-0.2, 0, 0]], step_m=0.02, max_steps=400)
    assert lines[0].closed
    out = write_field_lines_msh(tmp_path / "lines.msh", lines, length_scale=1e-3)
    n = len(lines[0].points_m) - 1
    assert (out["n_nodes"], out["n_segments"]) == (n, n)
    text = (tmp_path / "lines.msh").read_text(encoding="utf-8")
    assert '"|B| (T)"' in text and '"B (T)"' in text and "$Elements" in text
    # Millimetres: the seed at x = -0.2 m is written as -200.
    assert "-2.000000000000000e+02" in text


@pytest.mark.parametrize("kwargs, message", [
    ({"seeds_m": [[0, 0]]}, "seeds_m"),
    ({"step_m": 0}, "step_m"),
    ({"max_steps": 0}, "max_steps"),
    ({"max_turn_deg": 0}, "max_turn_deg"),
    ({"extent_m": {"min": [0, 0, 0], "max": [0, 1, 1]}}, "extent_m"),
])
def test_invalid_tracer_inputs_fail_loudly(kwargs, message):
    args = {"seeds_m": [[0.1, 0, 0]], "step_m": 0.01, **kwargs}
    with pytest.raises(ValueError, match=message):
        trace_field_lines(rotation, **args)
    with pytest.raises(ValueError, match="shape"):
        trace_field_lines(lambda p: np.zeros((len(p), 2)), [[0.1, 0, 0]], step_m=0.01)


def _layered(lower, upper):
    """B = lower for z < 0 and upper for z >= 0; equal Bz keeps div B = 0 at z = 0."""
    @_flat
    def field(p):
        out = np.tile(lower, (len(p), 1)).astype(float)
        out[p[:, 2] >= 0] = upper
        return out
    return field


def test_lines_step_across_refraction_but_stop_where_the_direction_reverses():
    """Guards: a kink below 90 degrees at the minimum step is crossed, a reversal stops."""
    box = {"min": [-1, -1, -0.1], "max": [1, 1, 0.1]}
    field = _layered([0.0, 0.0, 1.0], [1.0, 0.0, 1.0])
    lines, record = trace_field_lines(field, [[0, 0, -0.05]], step_m=0.01, max_steps=100,
                                      extent_m=box, both_directions=False)
    assert lines[0].reasons == {"forward": "left_extent"} and record["kinks"] >= 1, record
    assert lines[0].points_m[-1, 2] > 0.09 and lines[0].points_m[-1, 0] > 0.09
    lines, record = trace_field_lines(_layered([1.0, 0.0, 0.1], [-1.0, 0.0, 0.1]),
                                      [[0, 0, -0.005]], step_m=0.01,
                                      max_steps=100, extent_m=box, both_directions=False)
    assert lines[0].reasons == {"forward": "stagnation"}, (lines[0].reasons, record)
    assert abs(lines[0].points_m[-1, 2]) <= 2 * record["min_step_m"]
