"""Field lines of any vector field, seeded by equal flux and checked for it.

The field is a callable field_T(points) -> values with points of shape
(k, 3) in metres and values of shape (k, 3); a non-finite row means the
point lies outside the field's data. Any source fits: an HDiv-MMM solution
evaluated by its integral operator in unmeshed air, radia.Fld of a Radia
object, Biot-Savart of a coil, or an analytic field. Nothing is resampled
onto a mesh.

trace_field_lines integrates all seeds together: each RK4 stage evaluates
every active line in one call, so an expensive integral field is called
four times per step, not four times per step and line. The arc-length step
halves where the direction turns faster than max_turn_deg, and also before
a line leaves the data or the extent, so lines end within the minimum step
of the boundary. A trial step whose RK stages meet a field null or leave
the data is retried with half the step; only at the minimum step does the
line stop there. The turn of a step is the largest turn of its RK stages
and its end, since a stage beyond a kink skews the step even when the end
stops short of it. At the minimum step a turn the step cannot resolve is
accepted as an unresolved turn (a material interface, or field detail finer
than the minimum step) and counted in the record; the step then runs
straight to where the direction turns (located by bisection) and on along
the new direction, so a kink costs an error of order step^2 curvature, not
of order the step. A turn of 90 degrees or more is accepted only as an interface
crossing: B before and after must share a normal component, the one normal
to B_before - B_after, of at least interface_normal_fraction of |B| (normal
B is continuous at an interface, tangential B may flip, as at the equator
of a magnetised ball). Otherwise B runs head-on into its reverse and the
line stops as a reversal, distinct from a stagnation where |B| vanishes.

A traced line cannot prove that the exact line closes. A line counts as
closed when, after leaving its seed, a step passes within closure_tol_m of
the seed in the seed's direction (measured on the cubic Hermite arc of the
step, not its chord); the line ends at the point of that step nearest the
seed, not at the seed, and the gap is recorded. A tolerance at
or above the drift of a non-closing line per turn (a helix of small pitch)
marks it closed.

flux_account reports how the lines represent the seeded flux: how each
ended, how often it crosses the seed surface again in the seeded sense
(drawn twice) and a flux-weighted count of its observed returns through
the opposite-sign part (a count, not a coverage). Line density stands for
flux only along lines that are complete and drawn once.

write_field_lines_msh writes the lines as MSH 4.2 line elements with |B|
node data for display in Gmsh.

Seeding lines so that each carries equal flux is known: While and Forbes
(J. Eng. Math. 76, 81-100, 2012, doi:10.1007/s10665-011-9530-0) place seeds
in the volume from an updated density map so that the flux per line is
nearly constant; Zhong, Hu and Zou (Data Sci. J. 14, 6, 2015,
doi:10.5334/dsj-2015-006) apportion lines to the cells of a seed-surface
grid in proportion to their flux. radia.flux_seeding instead divides the
seed surface in flux coordinates so that every cell carries exactly equal
flux.

Seeding one sign of the surface flux draws each flux tube once only if
every line crosses the surface once in each sense. Magnetostatic lines need
not close: they may cover a surface or wander without returning, and then
cross again or never return. flux_account reports what the traced lines
actually do; it does not assume closure.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .flux_seeding import surface_crossings as _surface_crossings

LINE_ENDS = ("closed", "open", "truncated", "stagnation", "reversal")
REASONS = ("left_data", "left_extent", "closed", "stagnation", "reversal", "max_steps")
_OK, _LEFT_DATA, _LEFT_EXTENT, _STAGNANT = 0, 1, 2, 3
_STATUS_REASON = {_LEFT_DATA: "left_data", _LEFT_EXTENT: "left_extent",
                  _STAGNANT: "stagnation"}


@dataclass(frozen=True)
class FieldLine:
    """One traced line; points run along +B and seed_vertex is the seed's row.

    reasons holds the termination reason of each traced direction
    ("forward", and "backward" unless the line closed or was traced one way).
    field_T holds B at the points (None when the line came without values,
    e.g. read back from a file; such a line cannot be written).
    closure_gap_m is the distance by which a closed line misses its seed.
    """

    seed_index: int
    points_m: NDArray[np.float64]
    field_T: NDArray[np.float64] | None
    seed_vertex: int
    closed: bool
    reasons: dict = field(default_factory=dict)
    closure_gap_m: float | None = None

    @property
    def end(self) -> str:
        if self.closed:
            return "closed"
        reasons = set(self.reasons.values())
        if "max_steps" in reasons:
            return "truncated"
        if "stagnation" in reasons:
            return "stagnation"
        if "reversal" in reasons:
            return "reversal"
        return "open"


class _Field:
    """Batched evaluation with outside, extent and stagnation status."""

    def __init__(self, field_T: Callable, extent_m, stagnation_T: float):
        self.field_T = field_T
        self.lo = self.hi = None
        if extent_m is not None:
            self.lo, self.hi = (np.asarray(extent_m[k], dtype=float) for k in ("min", "max"))
        self.stagnation_T = stagnation_T
        self.calls = 0
        self.points = 0

    def __call__(self, points: NDArray):
        values = np.full(points.shape, np.nan)
        status = np.zeros(len(points), dtype=np.int64)
        inside = np.ones(len(points), dtype=bool)
        if self.lo is not None:
            inside = np.all((points >= self.lo) & (points <= self.hi), axis=1)
            status[~inside] = _LEFT_EXTENT
        if np.any(inside):
            got = np.asarray(self.field_T(points[inside]), dtype=float)
            if got.shape != (int(inside.sum()), 3):
                raise ValueError(f"field_T must return shape (k, 3), got {got.shape}")
            self.calls += 1
            self.points += int(inside.sum())
            values[inside] = got
        finite = np.all(np.isfinite(values), axis=1)
        status[inside & ~finite] = _LEFT_DATA
        norm = np.linalg.norm(np.where(finite[:, None], values, 0.0), axis=1)
        status[inside & finite & (norm <= self.stagnation_T)] = _STAGNANT
        return values, norm, status


def _nearest_on_step(p0, t0, p1, t1, target):
    """Point of a step nearest target, on the cubic Hermite arc of its end tangents.

    The chord would miss a point of a curved line by its sagitta (step^2
    curvature / 8); the Hermite arc follows the line to fourth order.
    Returns (s in [0, 1], point, distance).
    """
    length = float(np.linalg.norm(p1 - p0))
    m0, m1 = length * t0, length * t1

    def arc(s):
        s = np.atleast_1d(s)[:, None]
        s2, s3 = s * s, s * s * s
        return ((2 * s3 - 3 * s2 + 1) * p0 + (s3 - 2 * s2 + s) * m0
                + (3 * s2 - 2 * s3) * p1 + (s3 - s2) * m1)

    def dist(s):
        return np.linalg.norm(arc(s) - target, axis=1)

    grid = np.linspace(0.0, 1.0, 17)
    k = int(np.argmin(dist(grid)))
    lo, hi = grid[max(k - 1, 0)], grid[min(k + 1, len(grid) - 1)]
    ratio = (np.sqrt(5.0) - 1) / 2
    for _ in range(60):  # golden section on the bracketing interval
        a, b = hi - ratio * (hi - lo), lo + ratio * (hi - lo)
        if dist(a)[0] <= dist(b)[0]:
            hi = b
        else:
            lo = a
    s = 0.5 * (lo + hi)
    point = arc(s)[0]
    return s, point, float(np.linalg.norm(point - target))


def _shared_normal_fraction(b0, b1):
    """|B . n| / max(|B0|, |B1|) for the normal n to B0 - B1, shared by B0 and B1.

    The only plane across which B0 and B1 have equal normal components (div B
    = 0) is normal to B0 - B1 projected out; zero for B1 = -B0 (a head-on
    reversal), cos(theta) at a magnetised ball's surface.
    """
    b0, b1 = np.asarray(b0, dtype=float), np.asarray(b1, dtype=float)
    d = b0 - b1
    scale = max(float(np.linalg.norm(b0)), float(np.linalg.norm(b1)))
    if scale == 0 or not np.all(np.isfinite(d)):
        return 0.0
    dn = float(np.linalg.norm(d))
    normal = b0 - (b0 @ d) * d / dn ** 2 if dn > 0 else b0
    return float(np.linalg.norm(normal)) / scale


def _step_through_kink(probe, p0, t0, ds, sign, iterations=40):
    """One step of length ds that turns where B turns on the straight path p0 + s ds t0.

    The turn is located by bisection (where the direction has turned by half
    the end-to-end turn), so a kink at an interface is crossed with an error
    of order ds^2 curvature instead of the RK order-ds error. Returns
    (point, B, unit direction) or None (then the RK step is kept).
    """
    def direction(points):
        v, nv, st = probe(np.atleast_2d(points))
        if np.any(st != _OK):
            return None, None
        return v, _unit(v, nv, np.full(len(v), sign))

    _, u_end = direction(p0 + ds * t0)
    if u_end is None:
        return None
    half = 0.5 * np.arccos(np.clip(u_end[0] @ t0, -1.0, 1.0))
    if half <= 0:
        return None
    lo, hi = 0.0, 1.0
    for _ in range(iterations):
        mid = 0.5 * (lo + hi)
        _, u = direction(p0 + mid * ds * t0)
        if u is None:
            return None
        if np.arccos(np.clip(u[0] @ t0, -1.0, 1.0)) > half:
            hi = mid
        else:
            lo = mid
    _, u_after = direction(p0 + hi * ds * t0)
    if u_after is None:
        return None
    end = p0 + lo * ds * t0 + (1.0 - lo) * ds * u_after[0]
    v, u = direction(end)
    if v is None:
        return None
    return end, v[0], u[0]


def _stages_reverse(k, j, fraction):
    """True if a stage direction turns by 90 degrees or more with no shared normal."""
    first = k[0][j]
    for other in k[1:]:
        if first @ other[j] <= 0 and _shared_normal_fraction(first, other[j]) < fraction:
            return True
    return False


def _unit(values, norm, sign):
    return sign[:, None] * values / np.where(norm > 0, norm, 1.0)[:, None]


def trace_field_lines(
    field_T: Callable,
    seeds_m: ArrayLike,
    *,
    step_m: float,
    max_steps: int = 2000,
    both_directions: bool = True,
    adaptive: bool = True,
    max_turn_deg: float = 10.0,
    closure: bool = True,
    closure_tol_m: float | None = None,
    extent_m: dict | None = None,
    stagnation_T: float | None = None,
    interface_normal_fraction: float = 1e-2,
) -> tuple[list[FieldLine], dict]:
    """Trace field lines from seeds through any callable field.

    step_m is the maximum arc-length step (metres); it halves down to
    step_m / 64. extent_m = {"min": [x, y, z], "max": [x, y, z]} stops
    lines at a box. stagnation_T is the |B| at or below which a line stops
    (default 1e-12 of the median seed |B|). closure_tol_m is the closure
    tolerance (default step_m / 1000). interface_normal_fraction is the
    shared normal |B| fraction above which a turn of 90 degrees or more at
    the minimum step is crossed as an interface. Returns the lines in seed order
    (seeds that produced no step are left out) and an evaluation record.
    """
    seeds = np.asarray(seeds_m, dtype=float)
    if seeds.ndim != 2 or seeds.shape[1] != 3 or not np.all(np.isfinite(seeds)):
        raise ValueError("seeds_m must be finite points of shape (n, 3)")
    step = float(step_m)
    if not (np.isfinite(step) and step > 0):
        raise ValueError("step_m must be positive")
    if int(max_steps) < 1:
        raise ValueError("max_steps must be positive")
    if not 0 < float(max_turn_deg) <= 180:
        raise ValueError("max_turn_deg must be in (0, 180]")
    if not 0 < float(interface_normal_fraction) < 1:
        raise ValueError("interface_normal_fraction must be in (0, 1)")
    closure_tol = 1e-3 * step if closure_tol_m is None else float(closure_tol_m)
    if not (np.isfinite(closure_tol) and closure_tol > 0):
        raise ValueError("closure_tol_m must be positive")
    if extent_m is not None:
        if not isinstance(extent_m, dict) or set(extent_m) != {"min", "max"}:
            raise ValueError('extent_m must be {"min": [x, y, z], "max": [x, y, z]}')
        lo, hi = (np.asarray(extent_m[k], dtype=float) for k in ("min", "max"))
        if lo.shape != (3,) or hi.shape != (3,) or not np.all(lo < hi):
            raise ValueError("extent_m min must be below max in every coordinate")

    probe = _Field(field_T, extent_m, 0.0)
    b0, n0, s0 = probe(seeds)
    if stagnation_T is None:
        good = n0[s0 == _OK]
        stagnation_T = 1e-12 * float(np.median(good)) if len(good) else 0.0
    probe.stagnation_T = float(stagnation_T)
    s0[(s0 == _OK) & (n0 <= probe.stagnation_T)] = _STAGNANT

    directions = (1.0, -1.0) if both_directions else (1.0,)
    owner = np.repeat(np.arange(len(seeds)), len(directions))
    sign = np.tile(np.asarray(directions), len(seeds))
    n = len(owner)
    p = seeds[owner].copy()
    t = np.zeros((n, 3))
    ok0 = s0[owner] == _OK
    t[ok0] = _unit(b0[owner][ok0], n0[owner][ok0], sign[ok0])
    t_start = t.copy()
    ds = np.full(n, step)
    ds_min = step / 64.0
    arc = np.zeros(n)
    steps = np.zeros(n, dtype=np.int64)
    paths = [[seeds[owner[i]].copy()] for i in range(n)]
    values = [[b0[owner[i]].copy()] for i in range(n)]
    reason = np.array([("" if ok0[i] else _STATUS_REASON[int(s0[owner[i]])]) for i in range(n)],
                      dtype=object)
    closed = np.zeros(n, dtype=bool)
    gap = np.full(n, np.nan)
    active = ok0.copy()
    angle_max = np.radians(float(max_turn_deg))
    turns = np.zeros(n, dtype=np.int64)
    sharp = np.zeros(n, dtype=np.int64)

    while np.any(active):
        idx = np.nonzero(active)[0]
        # RK4 on the unit direction field, one batched evaluation per stage.
        k = [t[idx]]
        failed = np.zeros(len(idx), dtype=np.int64)
        for stage, frac in enumerate((0.5, 0.5, 1.0)):
            live = failed == _OK
            q = p[idx] + frac * ds[idx, None] * k[-1]
            v, nv, st = probe(q[live])
            ks = np.zeros_like(q)
            sub = np.nonzero(live)[0]
            okv = st == _OK
            ks[sub[okv]] = _unit(v[okv], nv[okv], sign[idx][sub[okv]])
            failed[sub[~okv]] = st[~okv]
            k.append(ks)
        live = failed == _OK
        p_new = p[idx] + ds[idx, None] / 6.0 * (k[0] + 2 * k[1] + 2 * k[2] + k[3])
        v_new = np.full((len(idx), 3), np.nan)
        t_new = np.zeros((len(idx), 3))
        sub = np.nonzero(live)[0]
        v, nv, st = probe(p_new[sub])
        okv = st == _OK
        v_new[sub[okv]] = v[okv]
        t_new[sub[okv]] = _unit(v[okv], nv[okv], sign[idx][sub[okv]])
        failed[sub[~okv]] = st[~okv]

        angle = np.arccos(np.clip(np.sum(t[idx] * t_new, axis=1), -1.0, 1.0))
        # A stage beyond a kink skews the RK step even when its end stops short of it.
        for stage in k[1:]:
            angle = np.maximum(angle, np.arccos(np.clip(np.sum(k[0] * stage, axis=1), -1.0, 1.0)))
        chord = np.linalg.norm(p_new - p[idx], axis=1)
        at_min = ds[idx] <= ds_min * 1.0001
        for j, i in enumerate(idx):
            if failed[j] != _OK:
                if not at_min[j]:
                    # A stage met a boundary or a null the exact line may miss: retry shorter.
                    ds[i] = max(ds_min, 0.5 * ds[i])
                    continue
                reason[i] = _STATUS_REASON[int(failed[j])]
                active[i] = False
                continue
            if adaptive and (angle[j] > angle_max or chord[j] < 0.5 * ds[i]):
                if not at_min[j]:
                    ds[i] = max(ds_min, 0.5 * ds[i])
                    continue
                if _stages_reverse(k, j, interface_normal_fraction):
                    reason[i] = "reversal"  # B runs head-on into its reverse: the step stalls
                    active[i] = False
                    continue
                # RK smears a kink over the step; step straight to the kink, then on.
                kinked = _step_through_kink(probe, p[i], t[i], ds[i], sign[i])
                if kinked is not None:
                    p_new[j], v_new[j], t_new[j] = kinked
                if t[i] @ t_new[j] <= 0:
                    if _shared_normal_fraction(values[i][-1], v_new[j]) < interface_normal_fraction:
                        reason[i] = "reversal"
                        active[i] = False
                        continue
                    sharp[i] += 1  # an interface where tangential B flips
                turns[i] += 1  # unresolved at the minimum step: accepted, less accurate
            p_old, t_old, v_old = p[i].copy(), t[i].copy(), values[i][-1]
            p[i] = p_new[j]
            t[i] = t_new[j]
            arc[i] += ds[i]
            steps[i] += 1
            paths[i].append(p_new[j].copy())
            values[i].append(v_new[j].copy())
            if adaptive and angle[j] < 0.5 * angle_max and ds[i] < step:
                ds[i] = min(step, 2.0 * ds[i])
            seed = seeds[owner[i]]
            if (closure and arc[i] > 4.0 * step and t[i] @ t_start[i] > 0.8
                    and np.linalg.norm(seed - p_old) < 2.0 * np.linalg.norm(p[i] - p_old)):
                s, nearest, miss = _nearest_on_step(p_old, t_old, p[i], t[i], seed)
            else:
                miss = np.inf
            if miss <= closure_tol:
                paths[i].pop()
                values[i].pop()
                if s > 0:  # end at the point of the step nearest the seed
                    paths[i].append(nearest)
                    values[i].append((1 - s) * v_old + s * v_new[j])
                gap[i] = miss
                closed[i] = True
                reason[i] = "closed"
                active[i] = False
                if len(directions) == 2:  # the backward partner is no longer needed
                    partner = i + 1 if sign[i] > 0 else i - 1
                    if active[partner]:
                        active[partner] = False
                        reason[partner] = "superseded"
                continue
            if steps[i] >= int(max_steps):
                reason[i] = "max_steps"
                active[i] = False

    lines = []
    width = len(directions)
    for s in range(len(seeds)):
        fwd = s * width
        if closed[fwd] or width == 1:
            if len(paths[fwd]) > 1:
                lines.append(FieldLine(s, np.asarray(paths[fwd]), np.asarray(values[fwd]), 0,
                                       bool(closed[fwd]), {"forward": str(reason[fwd])},
                                       float(gap[fwd]) if closed[fwd] else None))
            continue
        bwd = fwd + 1
        if closed[bwd]:
            # A loop closed backward: reverse it so the points run along +B.
            pts, vals = np.asarray(paths[bwd])[::-1], np.asarray(values[bwd])[::-1]
            lines.append(FieldLine(s, pts, vals, len(pts) - 1, True,
                                   {"forward": str(reason[bwd])}, float(gap[bwd])))
            continue
        merged = np.asarray(paths[bwd][::-1] + paths[fwd][1:])
        merged_values = np.asarray(values[bwd][::-1] + values[fwd][1:])
        if len(merged) > 1:
            lines.append(FieldLine(s, merged, merged_values, len(paths[bwd]) - 1, False,
                                   {"forward": str(reason[fwd]), "backward": str(reason[bwd])}))
    record = {"field_calls": probe.calls, "field_points": probe.points,
              "step_m": step, "min_step_m": ds_min, "stagnation_T": probe.stagnation_T,
              "n_seeds": len(seeds), "n_lines": len(lines),
              "closure_tol_m": closure_tol, "unresolved_turns": int(turns.sum()),
              "interface_turns_ge_90deg": int(sharp.sum()),
              "interface_normal_fraction": float(interface_normal_fraction)}
    return lines, record


def surface_crossings(line: FieldLine, chart) -> list[int]:
    """Signs (+1 along the chart's area normal) of a line's crossings of a seed surface.

    See radia.flux_seeding.surface_crossings; the crossing at the seed
    itself is not counted.
    """
    return _surface_crossings(line.points_m, chart, seed_vertex=line.seed_vertex,
                              closed=line.closed)


def flux_account(lines: list[FieldLine], chart, *, n_seeds: int, flux_per_line_Wb: float,
                 orientation: int, opposite_flux_Wb: float) -> tuple[dict, list[str]]:
    """How traced lines represent the seeded and the opposite flux.

    Returns the account (flux by line end, repeat and opposite crossings of
    the seed surface, the flux-weighted count of observed opposite crossings,
    which may be partial or repeated and is not a coverage, per-line detail) and
    warnings where line density stops standing for flux.
    """
    q = float(flux_per_line_Wb)
    ends = dict.fromkeys(LINE_ENDS, 0)
    per_line, repeats, opposite, one_sided = [], 0, 0, 0
    for line in lines:
        ends[line.end] += 1
        signs = surface_crossings(line, chart)
        same = sum(1 for c in signs if c == orientation)
        repeats += same
        opposite += len(signs) - same
        if not line.closed and "backward" not in line.reasons:
            one_sided += 1
        per_line.append({"seed_index": line.seed_index, "end": line.end,
                         "repeat_crossings": same, "opposite_crossings": len(signs) - same})
    not_traced = n_seeds - len(lines)
    observed = -q * opposite
    account = {
        "flux_per_line_Wb": q,
        "flux_by_end_Wb": {end: q * count for end, count in ends.items()},
        "lines_by_end": ends,
        "not_traced_lines": not_traced,
        "repeat_crossings": repeats,
        "opposite_crossings": opposite,
        "opposite_flux_Wb": opposite_flux_Wb,
        "opposite_crossing_flux_Wb": observed,
        "opposite_crossing_flux_ratio": (observed / opposite_flux_Wb if opposite_flux_Wb
                                         else None),
        "lines": per_line,
    }
    warnings = []
    if ends["truncated"]:
        warnings.append(f"{ends['truncated']} line(s) stopped at max_steps "
                        f"({abs(q) * ends['truncated']:.4g} Wb): beyond the cut the line "
                        "density under-represents flux; raise max_steps")
    if not_traced:
        warnings.append(f"{not_traced} seed(s) produced no line (outside the data or the "
                        "extent)")
    if one_sided:
        warnings.append(f"{one_sided} open line(s) were traced forward only "
                        "(both_directions=False): behind the seed surface they are missing")
    if repeats:
        warnings.append(f"lines cross the seed surface {repeats} more time(s) in the seeded "
                        "direction: there each tube is drawn again and line density "
                        "over-represents flux")
    return account, warnings


def write_field_lines_msh(path: str | Path, lines: list[FieldLine], *,
                          length_scale: float = 1.0, name: str = "field_lines") -> dict:
    """Write lines as MSH 4.2 two-node line elements with |B| and B node data.

    Each line is its own curve entity in one physical group; a closed line
    gets one more element from its last point back to its first, spanning
    the recorded closure gap. Coordinates are written
    as points_m / length_scale (1e-3 writes millimetres). Returns counts.
    """
    from radia.gmsh_post_export import (DataField, ElementBlock, Entity, NodeBlock,
                                        PayloadV42, PhysicalGroup, _emit_v42_msh)

    if not lines:
        raise ValueError("no field lines to write")
    if any(line.field_T is None for line in lines):
        raise ValueError("every line needs field_T values to be written")
    payload = PayloadV42(groups=[PhysicalGroup(1, 1, name)])
    magnitude, vectors = [], []
    next_id = 1
    for tag, line in enumerate(lines, start=1):
        pts = line.points_m / float(length_scale)
        vals = line.field_T
        ids = list(range(next_id, next_id + len(pts)))
        connectivity = [[a, b] for a, b in zip(ids[:-1], ids[1:])]
        if line.closed:
            connectivity.append([ids[-1], ids[0]])
        lo, hi = pts.min(axis=0), pts.max(axis=0)
        payload.entities.append(Entity(1, tag, (*map(float, lo), *map(float, hi)), [1]))
        payload.node_blocks.append(NodeBlock(1, tag, next_id, [tuple(map(float, q)) for q in pts]))
        payload.element_blocks.append(ElementBlock(1, tag, 1, connectivity))
        for node, b in zip(ids, vals):
            magnitude.append((node, float(np.linalg.norm(b))))
            vectors.append((node, tuple(map(float, b))))
        next_id += len(pts)
    payload.data_fields += [DataField("$NodeData", "|B| (T)", 1, magnitude),
                            DataField("$NodeData", "B (T)", 3, vectors)]
    with open(path, "w", encoding="utf-8", newline="\n") as stream:
        _emit_v42_msh(stream, payload)
    return {"path": str(path), "n_lines": len(lines), "n_nodes": next_id - 1,
            "n_segments": sum(len(b.connectivity) for b in payload.element_blocks)}
