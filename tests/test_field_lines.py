"""Surrogate contracts for field-agnostic field-line tracing (radia.field_lines).

Assumes analytic fields with exact answers:
- B = (-y, x, 0): circles about z through every seed;
- the straight wire B = mu0 I / (2 pi r) phi_hat: circles of constant r, and
  the plane y = 0, x in [r1, r2], z in [-h/2, h/2] carries
  mu0 I h ln(r2 / r1) / (2 pi);
- B = z_hat inside the ball |p| < R and no data (NaN) outside it;
- a uniformly magnetised ball: 2/3 b_r z_hat inside, a dipole outside.
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
    n = len(lines[0].points_m)
    assert (out["n_nodes"], out["n_segments"]) == (n, n)
    assert 0 <= lines[0].closure_gap_m <= 0.02 / 1000
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
    """Guards: refraction is crossed at the minimum step, a head-on reversal stops.

    Across z = 0 with equal Bz, B may turn by 90 degrees or more (tangential B
    flips, Bz = 0.1 carries on): the line crosses and continues along the upper
    B. Where B = +x meets B = -x at x = 0 the far direction leads back: the
    line stops there as a reversal.
    """
    box = {"min": [-1, -1, -0.1], "max": [1, 1, 0.1]}
    field = _layered([0.0, 0.0, 1.0], [1.0, 0.0, 1.0])
    lines, record = trace_field_lines(field, [[0, 0, -0.05]], step_m=0.01, max_steps=100,
                                      extent_m=box, both_directions=False)
    assert lines[0].reasons == {"forward": "left_extent"}, record
    assert record["interface_crossings"] >= 1 and record["unresolved_turns"] == 0, record
    assert lines[0].points_m[-1, 2] > 0.09 and lines[0].points_m[-1, 0] > 0.09
    lines, record = trace_field_lines(_layered([1.0, 0.0, 0.1], [-1.0, 0.0, 0.1]),
                                      [[0, 0, -0.005]], step_m=0.01,
                                      max_steps=400, extent_m=box, both_directions=False)
    assert lines[0].reasons == {"forward": "left_extent"}, (lines[0].reasons, record)
    assert record["interface_crossings"] == 1, record
    upper = lines[0].points_m[lines[0].points_m[:, 2] > 2 * record["min_step_m"]]
    assert len(upper) > 5 and np.all(np.diff(upper[:, 0]) < 0)  # along (-1, 0, 0.1)
    np.testing.assert_allclose(np.diff(upper[:, 2]) / np.diff(upper[:, 0]), -0.1, rtol=1e-9)
    head_on = _flat(lambda p: np.where((p[:, 0] < 0)[:, None], [1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]))
    lines, record = trace_field_lines(head_on, [[-0.005, 0, 0]], step_m=0.01, max_steps=100,
                                      extent_m=box, both_directions=False)
    assert lines[0].reasons == {"forward": "reversal"}, (lines[0].reasons, record)
    assert lines[0].end == "reversal"
    assert abs(lines[0].points_m[-1, 0]) <= 2 * record["min_step_m"]
    assert record["interface_crossings"] == 0, record


@pytest.mark.parametrize("eps", [0.005, 1e-6])
def test_a_near_tangential_interface_is_crossed_without_a_threshold(eps):
    """Guards: any shared normal B crosses; only the far side's direction decides.

    B = (1, 0, eps) for z < 0 and (-1, 0, eps) above: the exact line reaches
    z = 0 and goes on along (-1, 0, eps), however small eps > 0 is.
    """
    box = {"min": [-1, -1, -0.1], "max": [1, 1, 0.1]}
    lines, record = trace_field_lines(_layered([1.0, 0.0, eps], [-1.0, 0.0, eps]),
                                      [[0, 0, -1e-3 * eps / 0.005]], step_m=0.04, max_steps=30,
                                      extent_m=box, both_directions=False)
    pts = lines[0].points_m
    assert lines[0].end != "reversal" and record["interface_crossings"] == 1, (lines[0].reasons, record)
    upper = pts[pts[:, 2] > 0]
    assert len(upper) >= 5 and np.all(np.diff(upper[:, 0]) < 0)
    np.testing.assert_allclose(np.diff(upper[:, 2]) / np.diff(upper[:, 0]), -eps, rtol=1e-6)


def _curved_interface(p):
    """Interface z = q^2 (q = min(x, 0)); B = (1, 0, 0) below and (-1, 0, -4q) above.

    Each side is divergence-free and both share the normal component -2q
    against (-2q, 0, 1), so a line at z = a^2 meets the interface at x = -a
    and goes on above it along z = 2 x^2 - a^2 (Codex review of dbd26e866).
    """
    q = np.minimum(p[:, 0], 0.0)
    far = p[:, 2] >= q * q
    return np.column_stack((np.where(far, -1.0, 1.0), np.zeros(len(p)), np.where(far, -4 * q, 0.0)))


@pytest.mark.parametrize("a, x0", [(1e-4, -0.0031), (1e-4, -0.0015), (4e-7, -0.0015),
                                   (4e-7, -0.0031), (1e-8, -0.0008)])
@pytest.mark.parametrize("step", [0.064, 0.0064])
def test_a_curved_interface_is_crossed_or_left_unresolved_never_reversed(a, x0, step):
    """Guards: a traversable jump is never a reversal, and every point after it is on the far side.

    The interface bends back across a straight far piece longer than about 2a:
    the probe and the far piece shrink to what stays across, and a second jump
    within one minimum step stops the line as unresolved, not reversed.
    """
    lines, record = trace_field_lines(_flat(_curved_interface), [[x0, 0, a * a]], step_m=step,
                                      max_steps=30, both_directions=False, closure=False)
    assert lines, record
    line = lines[0]
    assert line.end != "reversal", (line.reasons, record)
    assert line.reasons["forward"] in ("unresolved_turn", "max_steps", "left_extent"), line.reasons
    pts = line.points_m
    above = pts[:, 2] >= np.minimum(pts[:, 0], 0) ** 2
    assert record["interface_crossings"] == 1, record
    first = int(np.argmax(above))
    assert first > 0 and np.all(above[first:]), (pts[first:][~above[first:]], record)


@pytest.mark.parametrize("step", [1.92, 0.96, 0.192])
def test_a_coarse_step_never_crosses_a_smooth_separatrix(step):
    """Guards: a smooth turn is no interface; the line stays on its side of x = 0.

    B = (-x, y, 0) is divergence-free with lines x y = const; from (0.02, 0.001)
    x > 0 for ever. A step far coarser than the turn must not step across it.
    """
    saddle = _flat(lambda p: p * np.array([-1.0, 1.0, 0.0]))
    lines, record = trace_field_lines(saddle, [[0.02, 0.001, 0]], step_m=step, max_steps=12,
                                      both_directions=False)
    assert record["interface_crossings"] == 0, record
    for line in lines:
        assert np.all(line.points_m[:, 0] > 0), (line.points_m[:, 0].min(), line.reasons, record)
        assert line.reasons["forward"] in ("unresolved_turn", "max_steps"), line.reasons


def test_a_smooth_head_on_null_stops_as_an_unresolved_turn():
    """Guards: a smooth reversal of direction (B = -x x_hat through x = 0) is not an interface."""
    null = _flat(lambda p: np.column_stack((-p[:, 0], 0 * p[:, 0], 0 * p[:, 0])))
    lines, record = trace_field_lines(null, [[-0.3, 0, 0]], step_m=0.05, max_steps=200,
                                      both_directions=False)
    assert lines[0].end in ("unresolved", "stagnation"), (lines[0].reasons, record)
    assert record["interface_crossings"] == 0
    assert np.all(lines[0].points_m[:, 0] < 0)


def _magnetised_ball(radius, b_r):
    """Uniformly magnetised ball: B = 2/3 b_r z_hat inside, the dipole of m = 4/3 pi a^3 b_r / mu0 outside."""
    @_flat
    def field(p):
        r = np.linalg.norm(p, axis=1)
        out = np.tile([0.0, 0.0, 2.0 * b_r / 3.0], (len(p), 1))
        o = r >= radius
        rhat = p[o] / r[o, None]
        c = b_r * radius ** 3 / (3.0 * r[o] ** 3)
        out[o] = c[:, None] * (3.0 * rhat * rhat[:, 2:3] - np.array([0.0, 0.0, 1.0]))
        return out
    return field


def test_lines_cross_a_magnetised_ball_surface_where_tangential_b_flips():
    """Guards: a line seeded near the equator crosses the surface, turning by 90 degrees or more.

    At polar angle theta the line turns by atan(tan theta) + atan(tan theta / 2),
    90 degrees or more beyond sin theta = sqrt(2/3); normal B (2/3 b_r cos theta)
    is shared, so the line goes on outside and returns, closed, through the seed.
    """
    a = 0.05
    seeds = [[0.95 * a, 0, 0], [0.999 * a, 0, 0]]
    lines, record = trace_field_lines(_magnetised_ball(a, 1.0), seeds, step_m=0.002, max_steps=3000)
    assert [line.end for line in lines] == ["closed", "closed"], ([l.reasons for l in lines], record)
    assert record["interface_crossings"] >= 4 and record["unresolved_turns"] == 0, record
    for line, rho in zip(lines, (0.95 * a, 0.999 * a)):
        inside = line.points_m[np.linalg.norm(line.points_m, axis=1) < a]
        # Every pass keeps the radius: the kink is located, not smeared over a
        # step (an RK step across it errs by about the minimum step, 3e-5 m).
        np.testing.assert_allclose(np.hypot(inside[:, 0], inside[:, 1]), rho,
                                   atol=1e-3 * record["min_step_m"])
        top = np.sqrt(a ** 2 - rho ** 2)
        assert 0 <= top - np.max(np.abs(inside[:, 2])) <= 2 * record["min_step_m"]
        assert np.linalg.norm(line.points_m, axis=1).max() > a  # it leaves the ball


def test_a_low_pitch_helix_is_not_reported_closed():
    """Guards: closure needs a near return within closure_tol_m, and nothing is snapped.

    B = (-y, x, 0.001) has helical lines of pitch 2 pi 0.001 m; at step 0.05 m
    each turn passes 6 mm from the seed, so the line stays open and z keeps
    increasing.
    """
    helix = _flat(lambda p: np.column_stack((-p[:, 1], p[:, 0], np.full(len(p), 1e-3))))
    lines, record = trace_field_lines(helix, [[1.0, 0, 0]], step_m=0.05, max_steps=300,
                                      both_directions=False)
    assert not lines[0].closed and lines[0].reasons == {"forward": "max_steps"}, record
    assert np.all(np.diff(lines[0].points_m[:, 2]) > 0)
    loose, _ = trace_field_lines(helix, [[1.0, 0, 0]], step_m=0.05, max_steps=300,
                                 both_directions=False, closure_tol_m=0.01)
    assert loose[0].closed and loose[0].closure_gap_m > 0.006  # a tolerance above the pitch
    assert np.all(np.diff(loose[0].points_m[:, 2]) > 0)  # still not snapped to the seed


def test_a_null_met_only_by_a_trial_stage_does_not_stop_the_line():
    """Guards: stage failures are retried with a shorter step before a line stops.

    B = (-y, x, 0) |p - (1, 0.05, 0)|^2 vanishes off the unit circle through the
    seed (1, 0, 0); a 0.1 m trial step samples that null at its first midpoint.
    """
    def field(p):
        scale = (p[:, 0] - 1) ** 2 + (p[:, 1] - 0.05) ** 2
        return np.column_stack((-p[:, 1] * scale, p[:, 0] * scale, 0 * scale))

    lines, record = trace_field_lines(field, [[1.0, 0, 0]], step_m=0.1, max_steps=400,
                                      both_directions=False)
    assert len(lines) == 1 and lines[0].closed, (lines, record)
    np.testing.assert_allclose(np.hypot(lines[0].points_m[:, 0], lines[0].points_m[:, 1]), 1,
                               atol=1e-3)
