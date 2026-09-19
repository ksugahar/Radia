"""Automatic fin ("beak") detection on a swept-conductor cross-section.

The IH fin-PEEC comparison needs three things that were previously hard
coded for one fixture (``x >= 1.6 mm`` beak, ``x >= 3.5 mm`` tip, probe at
``x = 5 mm``): which perimeter samples belong to a fin, where its tip is, and
where the workpiece-side field probe sits.  This module derives all of them
from the closed section outline alone, so any STEP whose section carries a
thin protrusion is handled without per-geometry constants.

Method (pure geometry, no field solve):

1. Local thickness ``t(s)``: from every outline sample cast a ray along the
   inward normal and take the distance to the first opposite edge.
2. Fin candidates: contiguous runs where ``t < fin_ratio * median(t)``.  The
   body dominates the perimeter, so the median is the body thickness.
3. The rounded tip has a *large* thickness along its normal (the ray runs
   back through the whole fin and body), so it splits a fin into two thin
   runs.  Short gaps between thin runs are merged when the gap arc length is
   below ``gap_factor`` times the neighbouring fin thickness.
4. Per fin: root points (run ends), axis (root midpoint -> farthest sample),
   tip radius (circumradius around the extremity), tip arc (samples whose
   local circumradius stays within ``1.5 * r_tip``), and a probe segment
   perpendicular to the axis beyond the extremity.

Everything is expressed in the outline's own units (metres in Radia), and
the analysis is invariant to rigid motion of the section.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

# ----------------------------------------------------------------------
# outline helpers
# ----------------------------------------------------------------------

def _as_closed_ccw(xy):
    p = np.asarray(xy, dtype=float)
    if p.ndim != 2 or p.shape[1] != 2 or p.shape[0] < 3:
        raise ValueError("outline must be an (n>=3, 2) array")
    if not np.all(np.isfinite(p)):
        raise ValueError("outline must be finite")
    if np.linalg.norm(p[0] - p[-1]) <= 1e-12 * np.ptp(p, axis=0).max():
        p = p[:-1]
    x, y = p[:, 0], p[:, 1]
    area2 = np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y)
    if abs(area2) <= 0:
        raise ValueError("outline has zero area")
    if area2 < 0:
        p = p[::-1].copy()
    return p


def outline_arclength(xy):
    """Cumulative arc length at each sample and the closed perimeter."""
    p = np.asarray(xy, dtype=float)
    seg = np.linalg.norm(np.roll(p, -1, axis=0) - p, axis=1)
    return np.r_[0.0, np.cumsum(seg)[:-1]], float(np.sum(seg))


def resample_outline(xy, n):
    """Equal-arc-length resampling of a closed CCW polygon (cell centred)."""
    p = _as_closed_ccw(xy)
    _s, per = outline_arclength(p)
    return outline_points_at(p, (np.arange(int(n)) + 0.5) * per / int(n))


def outline_points_at(xy, s_query):
    """Points on the closed polygon ``xy`` at arc lengths ``s_query``."""
    p = np.asarray(xy, dtype=float)
    seg = np.roll(p, -1, axis=0) - p
    seg_len = np.linalg.norm(seg, axis=1)
    cum = np.r_[0.0, np.cumsum(seg_len)]
    per = cum[-1]
    s = np.mod(np.asarray(s_query, dtype=float), per)
    idx = np.clip(np.searchsorted(cum, s, side="right") - 1, 0, len(p) - 1)
    t = (s - cum[idx]) / np.maximum(seg_len[idx], 1e-300)
    return p[idx] + t[:, None] * seg[idx]


def _inward_normals(p):
    tangent = np.roll(p, -1, axis=0) - np.roll(p, 1, axis=0)
    tangent /= np.maximum(np.linalg.norm(tangent, axis=1), 1e-300)[:, None]
    # CCW polygon: interior is to the left of the tangent.
    return np.column_stack((-tangent[:, 1], tangent[:, 0]))


def local_thickness(xy):
    """Distance from each sample along its inward normal to the outline."""
    p = _as_closed_ccw(xy)
    n = len(p)
    normals = _inward_normals(p)
    a = p
    b = np.roll(p, -1, axis=0)
    e = b - a                                            # (n, 2) edges
    scale = np.ptp(p, axis=0).max()
    eps = 1e-9 * scale
    # ray: p_i + t * n_i ; edge j: a_j + u * e_j
    # solve [n_i, -e_j] [t, u]^T = a_j - p_i
    d = a[None, :, :] - p[:, None, :]                    # (n, n, 2)
    nx, ny = normals[:, 0][:, None], normals[:, 1][:, None]
    ex, ey = e[:, 0][None, :], e[:, 1][None, :]
    det = nx * (-ey) - (-ex) * ny                        # (n, n)
    with np.errstate(divide="ignore", invalid="ignore"):
        t = (d[:, :, 0] * (-ey) - (-ex) * d[:, :, 1]) / det
        u = (nx * d[:, :, 1] - ny * d[:, :, 0]) / det
    ok = np.isfinite(t) & np.isfinite(u) & (t > eps) & (u >= -1e-9) & (u <= 1 + 1e-9)
    idx = np.arange(n)
    for shift in (0, -1):                                # own two edges
        ok[idx, (idx + shift) % n] = False
    t = np.where(ok, t, np.inf)
    thickness = t.min(axis=1)
    if not np.all(np.isfinite(thickness)):
        raise ValueError("inward ray missed the outline; is the polygon closed and simple?")
    return thickness


def _circumradius(p0, p1, p2):
    a = np.linalg.norm(p1 - p0, axis=-1)
    b = np.linalg.norm(p2 - p1, axis=-1)
    c = np.linalg.norm(p2 - p0, axis=-1)
    cross = np.abs((p1 - p0)[..., 0] * (p2 - p0)[..., 1]
                   - (p1 - p0)[..., 1] * (p2 - p0)[..., 0])
    with np.errstate(divide="ignore", invalid="ignore"):
        r = a * b * c / (2.0 * cross)
    return np.where(cross > 0, r, np.inf)


def _fit_circle(points):
    """Kasa algebraic circle fit; returns (radius, centre)."""
    q = np.asarray(points, dtype=float)
    if len(q) < 3:
        return np.inf, q.mean(axis=0)
    x, y = q[:, 0], q[:, 1]
    a_mat = np.column_stack((2 * x, 2 * y, np.ones_like(x)))
    rhs = x * x + y * y
    sol, *_ = np.linalg.lstsq(a_mat, rhs, rcond=None)
    cx, cy, c0 = sol
    r2 = c0 + cx * cx + cy * cy
    if r2 <= 0:
        return np.inf, np.array([cx, cy])
    return float(np.sqrt(r2)), np.array([cx, cy])


def _circular_runs(mask):
    """Start index and length of True runs on a circular boolean array."""
    m = np.asarray(mask, dtype=bool)
    n = len(m)
    if m.all():
        return [(0, n)]
    if not m.any():
        return []
    start = int(np.flatnonzero(~m)[0])           # rotate so we start on False
    rolled = np.roll(m, -start)
    runs = []
    i = 0
    while i < n:
        if rolled[i]:
            j = i
            while j < n and rolled[j]:
                j += 1
            runs.append(((i + start) % n, j - i))
            i = j
        else:
            i += 1
    return runs


# ----------------------------------------------------------------------
# result containers
# ----------------------------------------------------------------------

@dataclass(frozen=True)
class FinFeature:
    """One detected fin on a section outline (all lengths in outline units)."""
    root_a: np.ndarray
    root_b: np.ndarray
    root_mid: np.ndarray
    extremity: np.ndarray
    axis: np.ndarray                 # unit vector root_mid -> extremity
    normal: np.ndarray               # unit vector perpendicular to axis
    length: float                    # root_mid -> extremity
    root_thickness: float
    min_thickness: float
    tip_radius: float
    tip_arc_length: float
    tip_span: float                  # tip region depth along the axis (2 r_tip)
    arc_start: float                 # arc length (on the analysed outline)
    arc_end: float
    tip_arc_start: float
    tip_arc_end: float
    body_half_width: float           # half extent of the section normal to the axis
    root_tolerance: float = 0.0      # root-plane inclusion tolerance (~2 outline steps)

    # --- classification of arbitrary points (lanes, triangle centroids) ---
    def project(self, xy):
        q = np.asarray(xy, dtype=float)
        return (q - self.root_mid) @ self.axis

    def fin_mask(self, xy):
        """Points on or beyond the root plane (perpendicular to the axis).

        Samples lying on the root plane itself (e.g. the shoulder faces of a
        stepped root) count as fin side; ``root_tolerance`` (two steps of the
        analysed outline) absorbs the sampling error of the root corner.
        This is a half-space classification for one selected feature. If two
        fins extend into the same half-space, use the feature's connected
        outline interval (for example :func:`feature_panel_weights`) rather
        than combining multiple ``fin_mask`` results.
        """
        return self.project(xy) >= -self.root_tolerance

    def tip_mask(self, xy):
        return self.project(xy) >= self.length - self.tip_span - 1e-9 * self.length

    def probe_points(self, offset=None, n=41, half_width=None):
        """Segment perpendicular to the axis, ``offset`` beyond the tip.

        ``offset`` defaults to four tip radii (>= a few skin depths for the
        IH fixture); callers with a known skin depth should pass
        ``max(5*delta, 4*r_tip)`` so the probe stays outside the SIBC
        breakdown zone of the tip.
        """
        if offset is None:
            offset = 4.0 * self.tip_radius
        if half_width is None:
            half_width = self.body_half_width
        centre = self.extremity + float(offset) * self.axis
        if int(n) <= 1:
            return centre[None, :]
        eta = np.linspace(-float(half_width), float(half_width), int(n))
        return centre[None, :] + eta[:, None] * self.normal[None, :]

    def on_tip_arc(self, arclength, perimeter):
        """True where an outline arc length falls on the tip arc."""
        u = np.mod(np.asarray(arclength, dtype=float) - self.tip_arc_start,
                   perimeter)
        return u <= self.tip_arc_length + 1e-12 * perimeter

    def as_dict(self):
        return {
            "root_a_m": self.root_a.tolist(), "root_b_m": self.root_b.tolist(),
            "root_mid_m": self.root_mid.tolist(),
            "extremity_m": self.extremity.tolist(),
            "axis": self.axis.tolist(), "length_m": self.length,
            "root_thickness_m": self.root_thickness,
            "min_thickness_m": self.min_thickness,
            "tip_radius_m": self.tip_radius,
            "tip_arc_length_m": self.tip_arc_length,
            "tip_span_m": self.tip_span,
            "body_half_width_m": self.body_half_width,
            "root_tolerance_m": self.root_tolerance,
        }


@dataclass(frozen=True)
class SectionAnalysis:
    outline: np.ndarray
    arclength: np.ndarray
    perimeter: float
    thickness: np.ndarray
    body_thickness: float
    fins: tuple[FinFeature, ...] = field(default_factory=tuple)

    def arclength_of(self, xy):
        """Continuous arc coordinate from nearest outline-segment projection."""
        q = np.asarray(xy, dtype=float)
        p0 = self.outline
        edge = np.roll(p0, -1, axis=0) - p0
        edge2 = np.sum(edge * edge, axis=1)
        if np.any(edge2 <= 0):
            raise ValueError("outline contains a zero-length edge")
        rel = q[:, None, :] - p0[None, :, :]
        fraction = np.clip(np.sum(rel * edge[None, :, :], axis=2) /
                           edge2[None, :], 0.0, 1.0)
        projected = p0[None, :, :] + fraction[:, :, None] * edge[None, :, :]
        distance2 = np.sum((q[:, None, :] - projected) ** 2, axis=2)
        nearest = np.argmin(distance2, axis=1)
        edge_length = np.sqrt(edge2)
        return np.mod(self.arclength[nearest]
                      + fraction[np.arange(len(q)), nearest]
                      * edge_length[nearest], self.perimeter)

    @property
    def primary(self):
        """The longest fin, or None."""
        if not self.fins:
            return None
        return max(self.fins, key=lambda f: f.length)

    def as_dict(self):
        return {
            "perimeter_m": self.perimeter,
            "body_thickness_m": self.body_thickness,
            "n_fins": len(self.fins),
            "fins": [f.as_dict() for f in self.fins],
        }


# ----------------------------------------------------------------------
# detection
# ----------------------------------------------------------------------

def analyze_section(xy, *, fin_ratio=0.5, gap_factor=3.0, min_samples=4,
                    tip_window=3):
    """Detect fins on a closed section outline.

    Args:
        xy: (n, 2) closed outline, reasonably dense (>= ~256 samples for a
            0.25 mm tip on a 23 mm perimeter; use :func:`resample_outline`).
        fin_ratio: a sample is "thin" when its local thickness is below
            ``fin_ratio * median(thickness)``.
        gap_factor: thin runs separated by less than ``gap_factor`` times the
            neighbouring fin thickness are one fin (this absorbs the rounded
            tip, whose own normal thickness is large).
        min_samples: discard thin runs shorter than this.
        tip_window: retained for API compatibility (unused; the tip radius
            is a least-squares circle fit).
    """
    p = _as_closed_ccw(xy)
    n = len(p)
    if n < 64:
        raise ValueError("analyze_section needs a dense outline (>= 64 samples); "
                         "use resample_outline first")
    s, per = outline_arclength(p)
    ds = per / n
    t = local_thickness(p)
    t_body = float(np.median(t))
    thin = t < fin_ratio * t_body
    runs = _circular_runs(thin)
    runs = [r for r in runs if r[1] >= min_samples]
    if not runs:
        return SectionAnalysis(p, s, per, t, t_body, ())

    # merge runs across short gaps (the tip arc)
    runs.sort()
    merged = []
    for start, length in runs:
        if merged:
            ps, pl = merged[-1]
            gap = start - (ps + pl)
            neighbour_t = max(t[(ps + pl - 1) % n], t[start % n])
            if gap >= 0 and gap * ds < gap_factor * neighbour_t:
                merged[-1] = (ps, start + length - ps)
                continue
        merged.append((start, length))
    if len(merged) > 1:                                   # wrap-around merge
        (fs, fl), (ls, ll) = merged[0], merged[-1]
        gap = fs + n - (ls + ll)
        neighbour_t = max(t[(ls + ll - 1) % n], t[fs % n])
        if gap * ds < gap_factor * neighbour_t:
            merged[-1] = (ls, ll + gap + fl)
            merged.pop(0)

    # Extend each run back to the geometric root: the thin criterion stops
    # where the inward ray first reaches the body, which is inside the fin.
    # Walk outward while the thickness is still clearly below the body
    # thickness (a concave root corner or a body face ends the walk).
    extended = []
    for start, length in merged:
        limit = 0.9 * t_body
        max_ext = int(np.ceil(2.0 * max(t[(start + np.arange(length)) % n]) / ds))
        k = 0
        while k < max_ext and length < n and t[(start - 1) % n] < limit:
            start -= 1
            length += 1
            k += 1
        k = 0
        while k < max_ext and length < n and t[(start + length) % n] < limit:
            length += 1
            k += 1
        extended.append((start % n, min(length, n)))
    merged = extended

    centroid = p.mean(axis=0)
    fins = []
    for start, length in merged:
        idx = (start + np.arange(length)) % n
        q = p[idx]
        root_a, root_b = p[idx[0]], p[idx[-1]]
        root_mid = 0.5 * (root_a + root_b)
        far = int(np.argmax(np.linalg.norm(q - root_mid, axis=1)))
        extremity = q[far]
        axis = extremity - root_mid
        fin_len = float(np.linalg.norm(axis))
        if fin_len <= 0:
            continue
        axis /= fin_len
        # orient so the fin points away from the section centroid
        if (extremity - centroid) @ axis < 0:
            axis = -axis
        normal = np.array([-axis[1], axis[0]])
        # tip radius: algebraic (Kasa) circle fit on the samples around the
        # extremity, within +-0.6 * (minimum fin thickness) of arc length.
        # This is insensitive to chord polygonisation of the CAD arc, unlike
        # a three-point circumradius.
        t_min = float(t[idx].min())
        half_win = max(round(0.6 * t_min / ds), 3)
        win = idx[max(far - half_win, 0):min(far + half_win + 1, length)]
        r_tip, c_tip = _fit_circle(p[win])
        if not (np.isfinite(r_tip) and 0 < r_tip < 2.0 * t_min):
            r_tip = 0.5 * t_min
            c_tip = extremity - axis * r_tip
        # tip arc: contiguous samples around the extremity staying on the
        # fitted circle (2 % radial tolerance, at least half an outline step)
        tol = max(0.02 * r_tip, 0.5 * ds)
        on_circle = np.abs(np.linalg.norm(p - c_tip, axis=1) - r_tip) <= tol
        lo = far
        while lo - 1 >= 0 and on_circle[idx[lo - 1]]:
            lo -= 1
        hi = far
        while hi + 1 < length and on_circle[idx[hi + 1]]:
            hi += 1
        tip_arc_len = float((hi - lo + 1) * ds)
        proj_all = (p - root_mid) @ normal
        body_half = 0.5 * float(np.ptp(proj_all))
        fins.append(FinFeature(
            root_a=root_a, root_b=root_b, root_mid=root_mid,
            extremity=extremity, axis=axis, normal=normal, length=fin_len,
            root_thickness=float(np.linalg.norm(root_b - root_a)),
            min_thickness=float(t[idx].min()), tip_radius=r_tip,
            tip_arc_length=tip_arc_len, tip_span=2.0 * r_tip,
            arc_start=float(s[idx[0]]), arc_end=float(s[idx[-1]] + ds),
            tip_arc_start=float(s[idx[lo]]), tip_arc_end=float(s[idx[hi]] + ds),
            body_half_width=body_half, root_tolerance=2.0 * ds))
    return SectionAnalysis(p, s, per, t, t_body, tuple(fins))


# ----------------------------------------------------------------------
# lane placement
# ----------------------------------------------------------------------

def graded_lane_arclengths(analysis, n_lanes, *, tip_lanes=8, taper=None):
    """Arc-length positions of ``n_lanes`` lanes graded toward each fin tip.

    Body: uniform spacing ``h_body``.  Each fin face: spacing decreases
    geometrically from ``h_body`` at the root to ``h_tip`` at the tip arc.
    Tip arc: ``tip_lanes`` equal cells.  ``h_body`` is solved so that the
    total lane count is exactly ``n_lanes`` (cell-centred positions).
    Without fins this is plain equal-arc-length sampling.
    """
    n_lanes = int(n_lanes)
    if n_lanes < 8:
        raise ValueError("n_lanes must be >= 8")
    per = analysis.perimeter
    s = analysis.arclength
    n = len(s)
    if not analysis.fins:
        return (np.arange(n_lanes) + 0.5) * per / n_lanes

    # spacing target h(s) on the dense samples, up to the unknown body scale
    def spacing(h_body):
        h = np.full(n, h_body)
        for fin in analysis.fins:
            a, b = fin.arc_start, fin.arc_end
            arc = np.mod(s - a, per)
            on_fin = arc < (b - a) % per if b != a else np.zeros(n, bool)
            fin_len = (b - a) % per
            h_tip = fin.tip_arc_length / tip_lanes
            if h_tip >= h_body:
                continue
            u = arc[on_fin]
            hh = np.empty_like(u)
            t0 = (fin.tip_arc_start - a) % per
            t1 = t0 + fin.tip_arc_length
            ratio = h_tip / h_body if taper is None else float(taper)
            # face 1: root -> tip arc
            m1 = u < t0
            hh[m1] = h_body * ratio ** (u[m1] / max(t0, 1e-300))
            # tip arc
            m2 = (u >= t0) & (u <= t1)
            hh[m2] = h_tip
            # face 2: tip arc -> root
            m3 = u > t1
            hh[m3] = h_body * ratio ** ((fin_len - u[m3]) / max(fin_len - t1, 1e-300))
            h[on_fin] = hh
        return h

    seg = np.roll(s, -1) - s
    seg[-1] += per

    def count(h_body):
        return float(np.sum(seg / spacing(h_body)))

    lo, hi = per / (50.0 * n_lanes), per
    if count(hi) > n_lanes:
        raise ValueError("n_lanes too small for the requested tip resolution")
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if count(mid) > n_lanes:
            lo = mid
        else:
            hi = mid
    h = spacing(0.5 * (lo + hi))
    density = np.r_[0.0, np.cumsum(seg / h)]
    density *= n_lanes / density[-1]
    s_nodes = np.r_[s, per]
    return np.interp(np.arange(n_lanes) + 0.5, density, s_nodes)


def graded_lanes(outline, n_lanes, **kwargs):
    """Convenience: analyse ``outline`` and return (lanes_xy, analysis)."""
    analysis = analyze_section(outline)
    s_l = graded_lane_arclengths(analysis, n_lanes, **kwargs)
    return outline_points_at(analysis.outline, s_l), analysis


def feature_panel_weights(analysis, fin, lane_xy, *, tip=False):
    """Fraction of each perimeter lane's dual cell inside a fin feature.

    Centre-point masks jump when an irregular/graded cell crosses the root.
    These overlap weights make integrated current and loss independent of
    that arbitrary lane phase. Values lie in [0, 1].
    """
    xy = np.asarray(lane_xy, dtype=float)
    if xy.ndim != 2 or xy.shape[1] != 2 or len(xy) < 3:
        raise ValueError("lane_xy must have shape (n>=3, 2)")
    outline = analysis.outline
    s = analysis.arclength_of(xy)
    per = float(analysis.perimeter)
    s = np.unwrap(2.0 * np.pi * s / per) * per / (2.0 * np.pi)
    if np.any(np.diff(s) <= 1e-12 * per):
        raise ValueError("perimeter lanes must be ordered and distinct")
    prev = np.r_[s[-1] - per, s[:-1]]
    following = np.r_[s[1:], s[0] + per]
    left = 0.5 * (prev + s)
    right = 0.5 * (s + following)
    region = fin.tip_mask(outline) if tip else fin.fin_mask(outline)
    runs = _circular_runs(region)
    if not runs:
        return np.zeros(len(s))
    start_idx, run_length = max(runs, key=lambda item: item[1])
    dense_ds = per / len(outline)
    start = float(analysis.arclength[start_idx] - 0.5 * dense_ds)
    span = float(run_length * dense_ds)
    weight = np.zeros(len(s))
    for shift in range(-2, 4):
        a = start + shift * per
        b = a + span
        weight += np.maximum(0.0, np.minimum(right, b) - np.maximum(left, a))
    return np.clip(weight / (right - left), 0.0, 1.0)
