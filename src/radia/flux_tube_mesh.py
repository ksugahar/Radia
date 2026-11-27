"""Experimental: conforming hexahedral flux-tube meshes for Gmsh post-processing.

Status: experimental post-processing. The mesh shows where equal-flux tubes
run; it is not an analysis-ready volume mesh and no solver consumes it.

Input. A field callable field_T(points_m) -> (..., 3) in tesla, a seed chart
from radia.flux_seeding and equal-flux seed cells on it (ring_equal_flux_cells
for a disk-like chart, balanced_equal_flux_seeds otherwise). Every cell is the
cross-section of one tube.

Construction.
1. Lattice. A cell side carries every chart node of the cell sides it touches,
   so neighbouring tubes share their wall lines (no T-junctions). A cell with
   its xi = 0 side collapsed to a point (a disk centre) is split into three
   quads. A transition cell (RingFluxCells.transition_cells) keeps its own eta
   nodes; its outer side meets exactly twice as many rim nodes, joined by a
   conforming 2:1 layer of six quads per two coarse segments.
2. Walls. Every lattice point is traced as a field line (radia.field_lines).
   Sections are global: on every wall the section-k point is where the
   caller's angle angle_fn about one nesting curve equals level k, so a wall
   shared by several tubes carries the same points in all of them. Closed
   tubes use all n_sections levels. Optional pieces, all model-specific and
   given by the caller: an interface sphere whose crossings are mapped to two
   fixed levels, a radial family for open tubes (levels are spheres about
   radial_centre), and a target radius for the outermost open level.
3. Assembly. Wall nodes are numbered once per (wall, level); section
   interiors are Coons grids per tube; consecutive sections give 8-node hexes
   in Gmsh order.

Checks reported (not asserted): per-tube section flux by 3 x 3 Gauss on the
bilinear section quads over the seeded cell flux, scaled corner Jacobian,
edge ratio, a Bernstein lower bound of the trilinear Jacobian inside every
hex, and a face census (faces used more than twice, shared faces with equal
orientation, boundary faces whose nodes all belong to another tube).

Observable failure modes (raised): a lattice that does not settle, a
transition cell whose outer side is not an exact 2:1 refinement, untraced
walls (reported per tube). Collapsed centre cells are brought to one common
even eta count by halving their longest segments. Reported:
closed walls that do not make one turn of angle_fn, angle decreasing along a
wall, hexes with a non-positive Jacobian bound.

Supported (tested): a uniformly magnetized sphere seeded on its equatorial
disk or a spherical cap with ring cells, optionally with a 2:1 rim; a
straight wire seeded on a planar patch with balanced cells. Other fields need
their own angle_fn and nesting curve; nothing here detects a wrong choice
beyond the reports above.
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field
from typing import Callable

import numpy as np
from numpy.typing import NDArray

from .field_lines import trace_field_lines

STATUS = "experimental"
EPS = 1e-9
# Gmsh hexahedron: corner -> (edge neighbours in right-handed order)
_CORNERS = ((0, 1, 3, 4), (1, 2, 0, 5), (2, 3, 1, 6), (3, 0, 2, 7),
            (4, 7, 5, 0), (5, 4, 6, 1), (6, 5, 7, 2), (7, 6, 4, 3))
_EDGES = ((0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4),
          (0, 4), (1, 5), (2, 6), (3, 7))
# outward faces of a positively oriented Gmsh hexahedron
_FACES = ((0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7))
_REF = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
                 [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]], float)


# ---------------------------------------------------------------- seed lattice
class _Canon:
    """One spelling per chart coordinate: the first registered value within tol."""

    def __init__(self, tol=1e-8):
        self.tol, self.values = tol, []

    def __call__(self, v):
        v = float(v)
        for c in self.values:
            if abs(c - v) < self.tol:
                return c
        self.values.append(v)
        return v


def _graded_x(c, grading):
    """xi nodes of a rim cell, geometric in the distance to xi_s: grading = (xi_s, n) or
    (xi_s, n, to_p, from_p) to grade in p = to_p(xi) (e.g. a projected radius)."""
    x_s, n = grading[0], int(grading[1])
    if len(grading) == 4:
        to_p, from_p = grading[2:]
        d0, d1 = x_s - to_p(c[0]), x_s - to_p(c[1])
        x = from_p(x_s - d0 * (d1 / d0) ** (np.arange(n + 1) / n))
        x[0], x[-1] = c[0], c[1]
        return x
    d0, d1 = x_s - c[0], x_s - c[1]
    return x_s - d0 * (d1 / d0) ** (np.arange(n + 1) / n)


def _conforming_lattice(cells, m, periodic_eta, collapsed_xi0, b_split, grading, transition, key):
    """Per cell sorted xi nodes X_c and eta nodes Y_c closed under neighbour sides, and for
    each transition cell the fine eta nodes YF_c of its outer side."""
    transition = set(transition)
    X = [set(map(key, _graded_x(c, grading) if (grading and abs(c[1] - 1.0) < EPS)
                 else np.linspace(c[0], c[1], m + 1))) for c in cells]
    Y = [set(map(key, np.linspace(c[2], c[3], m + 1))) for c in cells]
    shifts = (0.0, 1.0, -1.0) if periodic_eta else (0.0,)

    def eta_touch(c, d):
        if min(c[1], d[1]) - max(c[0], d[0]) <= EPS:
            return False
        return any(abs(c[3] - (d[2] + s)) < EPS or abs(c[2] - (d[3] + s)) < EPS for s in shifts)

    def xi_touch(c, d):
        if not (abs(c[1] - d[0]) < EPS or abs(c[0] - d[1]) < EPS):
            return None
        for s in shifts:
            if min(c[3], d[3] + s) - max(c[2], d[2] + s) > EPS:
                return s
        return None

    n = len(cells)
    centre = [i for i in range(n) if collapsed_xi0 and abs(cells[i][0]) < EPS]
    if len({cells[i][1] for i in centre}) > 1:
        raise ValueError("collapsed centre cells must share one outer xi; got "
                         f"{sorted({float(cells[i][1]) for i in centre})}")
    for _ in range(100):
        changed = False
        for i in range(n):
            ci = cells[i]
            for j in range(n):
                if i == j:
                    continue
                cj = cells[j]
                if eta_touch(ci, cj):
                    add = {x for x in X[j] if ci[0] - EPS <= x <= ci[1] + EPS} - X[i]
                    if add:
                        X[i] |= add
                        changed = True
                s = xi_touch(ci, cj)
                if s is not None and i in transition and abs(ci[1] - cj[0]) < EPS:
                    s = None
                if s is not None:
                    add = {key(y + s) for y in Y[j] if ci[2] - EPS <= y + s <= ci[3] + EPS} - Y[i]
                    if add:
                        Y[i] |= add
                        changed = True
        if centre:
            # Collapsed cells share their xi nodes, so they need one even eta count: fill each
            # cell's whole deficit to the largest even count by halving its longest segment.
            target = max(len(Y[i]) - 1 for i in centre)
            target += target % 2
            for i in centre:
                while len(Y[i]) - 1 < target:
                    ys = np.array(sorted(Y[i]))
                    k = int(np.argmax(np.diff(ys)))
                    before = len(Y[i])
                    Y[i].add(key(0.5 * (ys[k] + ys[k + 1])))
                    if len(Y[i]) == before:
                        raise RuntimeError(f"collapsed cell {i}: eta segment [{ys[k]}, {ys[k + 1]}] "
                                           "is too short to split")
                    changed = True
        for i in centre:
            ny = len(Y[i]) - 1
            if ny % 2 or ny != max(len(Y[j]) - 1 for j in centre):
                continue
            nx = ny // 2 + b_split
            if len(X[i]) - 1 != nx:
                X[i] = set(map(key, np.linspace(cells[i][0], cells[i][1], nx + 1)))
                changed = True
        if not changed:
            YF = {}
            for i in transition:
                ci, f = cells[i], set()
                for j in range(n):
                    if j != i and abs(ci[1] - cells[j][0]) < EPS:
                        s = xi_touch(ci, cells[j])
                        if s is not None:
                            f |= {key(y + s) for y in Y[j] if ci[2] - EPS <= y + s <= ci[3] + EPS}
                f = np.array(sorted(f))
                nc = len(Y[i]) - 1
                if len(f) - 1 != 2 * nc or nc % 2:
                    raise RuntimeError(f"transition cell {i}: {nc} coarse vs {len(f) - 1} fine "
                                       "segments; a 2:1 layer needs an even coarse count")
                YF[i] = f
            return [np.array(sorted(x)) for x in X], [np.array(sorted(y)) for y in Y], YF
    raise RuntimeError("the seed lattice did not settle in 100 passes; eta counts "
                       f"{[len(y) - 1 for y in Y]}")


def _merge_points(points, tol):
    """Cluster points closer than tol; returns representatives and the owner of each point."""
    from scipy.spatial import cKDTree
    parent = np.arange(len(points))

    def root(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i
    for i, j in cKDTree(points).query_pairs(tol):
        ri, rj = root(i), root(j)
        if ri != rj:
            parent[max(ri, rj)] = min(ri, rj)
    roots = np.array([root(i) for i in range(len(points))])
    rep, owner = np.unique(roots, return_inverse=True)
    return points[rep], owner.ravel()


# ---------------------------------------------------------------- section grids
def _coons(bottom, right, top, left, u, v):
    """Grid (nu+1, nv+1, 3): bottom/top along u, left/right along v."""
    U, V = np.meshgrid(u, v, indexing="ij")
    U, V = U[..., None], V[..., None]
    p00, p10, p11, p01 = bottom[0], bottom[-1], top[-1], top[0]
    return ((1 - V) * bottom[:, None] + V * top[:, None] + (1 - U) * left[None, :] + U * right[None, :]
            - ((1 - U) * (1 - V) * p00 + U * (1 - V) * p10 + U * V * p11 + (1 - U) * V * p01))


def _seg(p, q, h):
    return p + np.linspace(0, 1, h + 1)[:, None] * (q - p)


def _par(xs):
    xs = np.asarray(xs, float)
    return (xs - xs[0]) / (xs[-1] - xs[0])


# one 2:1 unit: coarse nodes 0, 1, 2; fine f0..f4; interior a, m, b. Quads CCW in (row, layer).
_TEMPLATE = ((0, 1, "m", "a"), (1, 2, "b", "m"), (0, "a", "f1", "f0"), ("a", "m", "f2", "f1"),
             ("m", "b", "f3", "f2"), (2, "f4", "f3", "b"))


def _template_quads(lookup):
    """The six quads of one 2:1 unit, each a 2 x 2 grid (layer, row): p0, p3, p2, p1."""
    out = []
    for q in _TEMPLATE:
        p = [lookup(v) for v in q]
        g = (np.empty((2, 2) + p[0].shape) if isinstance(p[0], np.ndarray)
             else np.full((2, 2), None, dtype=object))
        g[0, 0], g[1, 0], g[1, 1], g[0, 1] = p[0], p[3], p[2], p[1]
        out.append(g)
    return out


def _transition_units(c, f):
    """2:1 layer between a coarse row c (2u + 1 nodes) and a fine row f (4u + 1 nodes)."""
    quads = []
    for u in range((len(c) - 1) // 2):
        cc, ff = c[2 * u: 2 * u + 3], f[4 * u: 4 * u + 5]
        pts = {"a": 0.5 * (0.5 * (cc[0] + cc[1]) + ff[1]), "m": 0.5 * (cc[1] + ff[2]),
               "b": 0.5 * (0.5 * (cc[1] + cc[2]) + ff[3])}
        quads += _template_quads(lambda v: cc[v] if isinstance(v, int) else
                                 (ff[int(v[1])] if v[0] == "f" else pts[v]))
    return quads


def _cell_blocks(sides, X, Y, collapsed, trans):
    """Structured blocks of one section from the cell's boundary polylines bottom (+xi),
    right (+eta), top (+xi), left (+eta)."""
    bottom, right, top, left = sides
    if trans:
        G = _coons(bottom, right[::2], top, left, _par(X), _par(Y))
        return [G[:-1]] + _transition_units(G[-2], right)
    if not collapsed:
        return [_coons(bottom, right, top, left, _par(X), _par(Y))]
    nx, ny = len(X) - 1, len(Y) - 1
    a, b = ny // 2, nx - ny // 2
    p1, p2, p3 = bottom[a], right[a], top[a]
    g = (p1 + p2 + p3) / 3.0
    qa = _coons(bottom[a:], right[:a + 1], _seg(g, p2, b), _seg(p1, g, a), _par(X[a:]), _par(Y[:a + 1]))
    qb = _coons(right[a:], top[a:][::-1], _seg(g, p3, a), _seg(p2, g, b), _par(Y[a:]), _par(X[a:][::-1]))
    qo = _coons(top[:a + 1][::-1], bottom[:a + 1], _seg(g, p1, a), _seg(p3, g, a),
                _par(X[:a + 1][::-1]), _par(X[:a + 1]))
    return [qa, qb, qo]


def _block_labels(side_ids, nx, ny, collapsed, trans):
    """Per block an object grid of node labels: ("w", wall) on the cell boundary, ("s", ...)
    on a seam inside the cell, None inside a block."""
    bot, rig, top, lef = side_ids

    def grid(nu, nv):
        return np.full((nu + 1, nv + 1), None, dtype=object)
    if trans:
        L = grid(nx - 1, ny)
        for i in range(nx):
            L[i, 0], L[i, ny] = ("w", bot[i]), ("w", top[i])
        for j in range(ny + 1):
            L[0, j] = ("w", lef[j])
            if 0 < j < ny:
                L[nx - 1, j] = ("s", 10, j)
        c = [L[nx - 1, j] for j in range(ny + 1)]
        f = [("w", r) for r in rig]
        out = [L]
        for u in range(ny // 2):
            cc, ff = c[2 * u: 2 * u + 3], f[4 * u: 4 * u + 5]
            pts = {k: ("s", 11, u, k) for k in "amb"}
            out += _template_quads(lambda v: cc[v] if isinstance(v, int) else
                                   (ff[int(v[1])] if v[0] == "f" else pts[v]))
        return out
    if not collapsed:
        L = grid(nx, ny)
        for i in range(nx + 1):
            L[i, 0], L[i, ny] = ("w", bot[i]), ("w", top[i])
        for j in range(ny + 1):
            L[0, j], L[nx, j] = ("w", lef[j]), ("w", rig[j])
        return [L]
    a, b = ny // 2, nx - ny // 2

    def s(k, j):
        return ("s", k, j)
    qa = grid(b, a)
    for i in range(b + 1):
        qa[i, 0], qa[i, a] = ("w", bot[a + i]), (s(2, i) if i < b else ("w", rig[a]))
    for j in range(a + 1):
        qa[b, j] = ("w", rig[j])
        qa[0, j] = ("w", bot[a]) if j == 0 else s(1, a - j)
    qb = grid(a, b)
    for i in range(a + 1):
        qb[i, 0] = ("w", rig[a + i])
        qb[i, b] = s(3, i) if i < a else ("w", top[a])
    for j in range(b + 1):
        qb[a, j] = ("w", top[nx - j])
        qb[0, j] = ("w", rig[a]) if j == 0 else s(2, b - j)
    qo = grid(a, a)
    for i in range(a + 1):
        qo[i, 0] = ("w", top[a - i])
        qo[i, a] = s(1, i) if i < a else ("w", bot[a])
    for j in range(a + 1):
        qo[a, j] = ("w", bot[j])
        qo[0, j] = ("w", top[a]) if j == 0 else s(3, a - j)
    qa[0, a] = qb[0, b] = qo[0, a] = ("s", 0, 0)
    return [qa, qb, qo]


def _grid_flux(field_T, grid, order=3):
    """Flux through the bilinear quads of a structured grid (order x order Gauss)."""
    t, w = np.polynomial.legendre.leggauss(order)
    t, w = 0.5 * (t + 1), 0.5 * w
    uu, vv = np.meshgrid(t, t, indexing="ij")
    ww = np.outer(w, w)
    a, b, c, d = grid[:-1, :-1], grid[1:, :-1], grid[1:, 1:], grid[:-1, 1:]
    U, V = uu[None, None, ..., None], vv[None, None, ..., None]
    A, B, C, D = (x[:, :, None, None, :] for x in (a, b, c, d))
    p = (1 - U) * (1 - V) * A + U * (1 - V) * B + U * V * C + (1 - U) * V * D
    du = (1 - V) * (B - A) + V * (C - D)
    dv = (1 - U) * (D - A) + U * (C - B)
    bf = np.asarray(field_T(p.reshape(-1, 3)), dtype=float).reshape(p.shape)
    return float(np.sum(ww[None, None] * np.einsum("...k,...k->...", bf, np.cross(du, dv))))


# ---------------------------------------------------------------- walls
class _Wall:
    """A traced wall line with its angle phi (increasing along +B)."""

    def __init__(self, line, angle_fn, sense, phi_ref):
        self.closed = line.closed
        if line.closed:
            pts = np.roll(line.points_m, -line.seed_vertex, axis=0)
            pts = np.vstack([pts, pts[:1]])
            k0 = 0
            if phi_ref is not None:
                pts = self._reroot(pts, angle_fn, sense, phi_ref)
        else:
            pts, k0 = line.points_m, line.seed_vertex
        th = np.unwrap(angle_fn(pts))
        if phi_ref is None or line.closed:
            th0 = float(th[k0])
        else:
            u = sense * (th[k0] - phi_ref)
            th0 = float(phi_ref + sense * 2 * np.pi * np.round(u / (2 * np.pi)))
        phi = sense * (th - th0)
        self.pts, self.phi, self.th_ref = pts, phi, th0
        self.k_ref = k0 if (phi_ref is None or line.closed) else int(
            np.clip(np.searchsorted(np.maximum.accumulate(phi), 0.0), 0, len(phi) - 1))
        d = np.diff(phi)
        self.monotone_violation = float(max(0.0, -d.min())) if len(d) else 0.0
        self.turn = float(phi[-1]) if self.closed else None
        self.lo, self.hi = float(phi.min()), float(phi.max())
        self.kmap = None

    @staticmethod
    def _reroot(cyc, angle_fn, sense, phi_ref):
        """Roll a closed cycle (first == last) to start where the angle is phi_ref."""
        th = np.unwrap(angle_fn(cyc))
        u = sense * (th - phi_ref)
        T = 2 * np.pi * np.ceil(u[0] / (2 * np.pi) - 1e-12)
        if abs(T - u[0]) < 1e-12:
            return cyc
        j = int(np.nonzero((u[:-1] <= T) & (u[1:] > T))[0][0])
        a, b = cyc[j], cyc[j + 1]
        lo, hi = 0.0, 1.0
        for _ in range(55):
            mid = 0.5 * (lo + hi)
            v = sense * (angle_fn((a + mid * (b - a))[None])[0] - phi_ref)
            v += 2 * np.pi * np.round((u[j] - v) / (2 * np.pi))
            lo, hi = (mid, hi) if v < T else (lo, mid)
        x = a + 0.5 * (lo + hi) * (b - a)
        body = cyc[:-1]
        return np.vstack([x[None], body[j + 1:], body[:j + 1], x[None]])

    def point_at(self, targets, angle_fn, sense):
        """Points where phi equals the targets (bisection of the angle along a chord)."""
        targets = np.atleast_1d(np.asarray(targets, float))
        phi, pts = self.phi, self.pts
        if self.closed:
            targets = np.mod(targets, 2 * np.pi)
            if self.kmap is not None:
                targets = np.interp(targets, self.kmap[1], self.kmap[0])
            targets[targets > phi[-1]] = phi[-1]
        run = np.maximum.accumulate(phi)
        k = np.clip(np.searchsorted(run, targets, side="right") - 1, 0, len(phi) - 2)
        a, b = pts[k], pts[k + 1]
        lo, hi = np.zeros(len(k)), np.ones(len(k))
        for _ in range(55):
            mid = 0.5 * (lo + hi)
            ph = sense * (angle_fn(a + mid[:, None] * (b - a)) - self.th_ref)
            ph = ph + 2 * np.pi * np.round((phi[k] - ph) / (2 * np.pi))
            below = ph < targets
            lo, hi = np.where(below, mid, lo), np.where(below, hi, mid)
        w = 0.5 * (lo + hi)
        self.last_frac_index = k + w
        return a + w[:, None] * (b - a)

    def interface_angles(self, radius, centre, angle_fn, sense):
        """phi where a closed wall crosses |x - centre| = radius, in order."""
        r = np.linalg.norm(self.pts - centre, axis=1) - radius
        sgn = np.where(np.abs(r) <= 1e-7 * radius, 0, np.sign(r))
        nz = np.nonzero(sgn)[0]
        out = []
        for p_, q_ in zip(nz[:-1], nz[1:]):
            if sgn[p_] == sgn[q_]:
                continue
            on = [j for j in range(p_ + 1, q_) if sgn[j] == 0]
            if on:
                out.append(float(self.phi[on[len(on) // 2]]))
                continue
            a, e = self.pts[p_] - centre, self.pts[p_ + 1] - self.pts[p_]
            A, B, C = e @ e, 2 * a @ e, a @ a - radius ** 2
            disc = np.sqrt(max(B * B - 4 * A * C, 0.0))
            ws = [w for w in ((-B - disc) / (2 * A), (-B + disc) / (2 * A)) if 0 <= w <= 1]
            x = self.pts[p_] + (ws[0] if ws else 0.5) * e
            ph = sense * (angle_fn(x[None])[0] - self.th_ref)
            ph += 2 * np.pi * np.round((self.phi[p_] - ph) / (2 * np.pi))
            out.append(float(ph))
        return out

    def radial_point(self, R, sign, centre):
        """First point at |x - centre| = R going forward (sign > 0) or backward from the
        phi = 0 point (a closed wall: up to its farthest point); (point, fractional index)
        or None."""
        r = np.linalg.norm(self.pts - centre, axis=1)
        if self.closed:
            j = int(np.argmax(r))
            idx = np.arange(0, j + 1) if sign > 0 else np.arange(len(r) - 1, j - 1, -1)
        else:
            idx = np.arange(self.k_ref, len(r)) if sign > 0 else np.arange(self.k_ref, -1, -1)
        hit = np.nonzero(r[idx][1:] >= R)[0]
        if len(hit) == 0:
            return None
        j0, j1 = idx[hit[0]], idx[hit[0] + 1]
        a, e = self.pts[j0] - centre, self.pts[j1] - self.pts[j0]
        A, B, C = e @ e, 2 * a @ e, a @ a - R * R
        w = (-B + np.sqrt(max(B * B - 4 * A * C, 0.0))) / (2 * A) if A > 0 else 0.0
        w = min(max(w, 0.0), 1.0)
        return self.pts[j0] + w * e, j0 + (w if j1 > j0 else -w)


# ---------------------------------------------------------------- hex checks
def hex_corner_jacobians(nodes_m: NDArray, hexes: NDArray) -> NDArray:
    """Corner Jacobian determinants (n, 8) of Gmsh-ordered hexes (0-based node indices)."""
    h = np.asarray(nodes_m, float)[np.asarray(hexes)]
    return np.stack([np.einsum("ij,ij->i", np.cross(h[:, i] - h[:, c], h[:, j] - h[:, c]),
                               h[:, k] - h[:, c]) for c, i, j, k in _CORNERS], axis=1)


def hex_quality(nodes_m: NDArray, hexes: NDArray) -> tuple[NDArray, NDArray]:
    """Minimum scaled corner Jacobian and max/min edge ratio of every hex."""
    h = np.asarray(nodes_m, float)[np.asarray(hexes)]
    scj = []
    for c, i, j, k in _CORNERS:
        e = [h[:, i] - h[:, c], h[:, j] - h[:, c], h[:, k] - h[:, c]]
        det = np.einsum("ij,ij->i", np.cross(e[0], e[1]), e[2])
        scj.append(det / np.prod([np.linalg.norm(x, axis=1) for x in e], axis=0))
    el = np.stack([np.linalg.norm(h[:, i] - h[:, j], axis=1) for i, j in _EDGES], 1)
    return np.min(scj, axis=0), el.max(1) / el.min(1)


def _bernstein_matrix():
    """Values at {0, 1/2, 1}^3 of the trilinear Jacobian -> its degree-(2, 2, 2) Bernstein
    coefficients (per axis b0 = f0, b2 = f1, b1 = 2 f(1/2) - (f0 + f1) / 2)."""
    m1 = np.array([[1.0, 0.0, 0.0], [-0.5, 2.0, -0.5], [0.0, 0.0, 1.0]])
    return np.einsum("ai,bj,ck->abcijk", m1, m1, m1).reshape(27, 27)


def hex_jacobian_lower_bound(nodes_m: NDArray, hexes: NDArray, chunk: int = 20000) -> NDArray:
    """Lower bound of det(dx/du) over the whole reference cube of every Gmsh-ordered
    trilinear hex: the smallest Bernstein coefficient of the degree-(2, 2, 2) polynomial,
    less a floating-point margin. Positive means the Jacobian is positive everywhere inside,
    not only at corners. Each element is moved to its first node before differentiating, so
    far-from-origin coordinates do not cancel, and the margin 2048 eps S^3 (S the sum of the
    node distances from that node, which bounds every Jacobian column) covers the rounding
    of the derivatives, determinants and Bernstein combination."""
    t = np.array([0.0, 0.5, 1.0])
    U, V, W = (x.ravel() for x in np.meshgrid(t, t, t, indexing="ij"))
    P = np.stack([U, V, W], 1)                                        # (27, 3)
    f = np.where(_REF[None] > 0, P[:, None], 1 - P[:, None])          # (27, 8, 3)
    s = 2 * _REF - 1                                                  # d/dt of the factor
    D = np.stack([s[None, :, a] * np.prod(np.delete(f, a, axis=2), axis=2) for a in range(3)], -1)
    M = _bernstein_matrix()
    X = np.asarray(nodes_m, float)
    hexes = np.asarray(hexes)
    out = np.empty(len(hexes))
    for s0 in range(0, len(hexes), chunk):
        h = X[hexes[s0:s0 + chunk]]                                   # (n, 8, 3)
        h = h - h[:, :1]                                              # translation-free
        J = np.einsum("pna,hnk->hpka", D, h)                          # (n, 27, 3, 3)
        margin = 2048 * np.finfo(float).eps * np.linalg.norm(h, axis=2).sum(axis=1) ** 3
        out[s0:s0 + chunk] = (np.linalg.det(J) @ M.T).min(axis=1) - margin
    return out


def face_census(hexes: NDArray, tube_of: NDArray) -> dict:
    """Faces of positively oriented hexes: a conforming mesh uses every face at most twice and
    a shared face in opposite orientations. boundary_faces_on_another_tube counts boundary
    faces whose four nodes all belong to one other tube (a non-matching joint)."""
    hexes, tube_of = np.asarray(hexes), np.asarray(tube_of)
    F = hexes[:, np.array(_FACES)].reshape(-1, 4)                     # (6n, 4)
    owner = np.repeat(np.arange(len(hexes)), 6)
    key = np.sort(F, axis=1)
    _, inv, cnt = np.unique(key, axis=0, return_inverse=True, return_counts=True)
    inv = inv.ravel()
    used = cnt[inv]
    two = np.nonzero(used == 2)[0]
    order = two[np.argsort(inv[two], kind="stable")]
    a, b = order[0::2], order[1::2]

    def rot_min(f):          # cyclic sequence starting at the smallest node
        k = np.argmin(f, axis=1)
        return np.stack([f[np.arange(len(f)), (k + i) % 4] for i in range(4)], 1)
    same = np.all(rot_min(F[a]) == rot_min(F[b]), axis=1)
    cross = tube_of[owner[a]] != tube_of[owner[b]]
    node_tubes: dict = {}
    for t_, h in zip(tube_of, hexes):
        for n in h:
            node_tubes.setdefault(int(n), set()).add(int(t_))
    suspicious = 0
    for i in np.nonzero(used == 1)[0]:
        common = set.intersection(*[node_tubes[int(n)] for n in F[i]]) - {int(tube_of[owner[i]])}
        suspicious += bool(common)
    return {"faces_used_more_than_twice": int(np.sum(cnt > 2)),
            "interior_faces": int(len(a)),
            "interior_faces_same_orientation": int(np.sum(same)),
            "tube_to_tube_shared_faces": int(np.sum(cross)),
            "boundary_faces": int(np.sum(used == 1)),
            "boundary_faces_on_another_tube": int(suspicious),
            "shared_face_pairs": (owner[a], owner[b])}


# ---------------------------------------------------------------- build
@dataclass
class FluxTubeMesh:
    """Conforming hex flux tubes. hexes are 0-based node indices in Gmsh order, oriented so
    that the corner Jacobians are positive; tube_of is the seed cell of every hex."""

    nodes_m: NDArray[np.float64]
    hexes: NDArray[np.int64]
    tube_of: NDArray[np.int64]
    tubes: list
    info: dict = field(default_factory=dict)
    status: str = STATUS


def build_flux_tube_mesh(
    field_T: Callable, patch, seeds, *, angle_fn: Callable, n_sections: int, step_m: float,
    max_steps: int, extent_m: dict | None = None, cells_per_side: int = 3,
    periodic_eta: bool = False, collapsed_xi0: bool = False, centre_split: int | None = None,
    rim_grading: tuple | None = None, phi_ref: float | None = None,
    interface: tuple | None = None, radial_centre=None, radial_max: float | None = None,
    level_target_radius: float | None = None,
) -> FluxTubeMesh:
    """Build conforming hex flux tubes from equal-flux seed cells (experimental).

    seeds: radia.flux_seeding cells (RingFluxCells or BalancedFluxSeeds) on patch.
    angle_fn(points) -> angle about the nesting curve; n_sections levels per turn.
    cells_per_side: chart segments per cell side before neighbour closure.
    periodic_eta: eta wraps (disk, cap). collapsed_xi0: the xi = 0 side is one point.
    centre_split: b of the three-quad split (nx = ny / 2 + b), default max(1, m // 2).
    rim_grading: (xi_s, n[, to_p, from_p]) xi nodes of cells at xi = 1, geometric toward xi_s.
    phi_ref: angle where phi = 0 when the seeds do not lie on one angle level.
    interface: (centre, radius) sphere whose two crossings per closed wall map to two levels.
    radial_centre: open tubes use spheres about it as levels; radial_max caps their radius.
    level_target_radius: move the outermost open level to this radius when it fits.
    Model-specific choices are the caller's; see the module docstring for what is checked.
    """
    t_start = time.perf_counter()
    m = int(cells_per_side)
    if int(n_sections) < 4:
        raise ValueError("n_sections must be at least 4")
    n_sec = int(n_sections)
    key = _Canon()
    cells = np.array([[key(v) for v in c] for c in np.asarray(seeds.cells, float)])
    transition = tuple(getattr(seeds, "transition_cells", ()))
    X, Y, YF = _conforming_lattice(cells, m, periodic_eta, collapsed_xi0,
                                   centre_split or max(1, m // 2), rim_grading, transition, key)
    side_chart = []
    for ci, (c, xs, ys) in enumerate(zip(cells, X, Y)):
        yr = YF.get(ci, ys)
        side_chart.append([(xs, np.full_like(xs, c[2])), (np.full_like(yr, c[1]), yr),
                           (xs, np.full_like(xs, c[3])), (np.full_like(ys, c[0]), ys)])
    allxe = np.vstack([np.column_stack(s) for sc in side_chart for s in sc])
    allpts = np.asarray(patch.chart(allxe[:, 0], allxe[:, 1]), float)
    scale = float(np.max(np.ptp(allpts, axis=0)))
    unique, owner = _merge_points(allpts, 1e-9 * scale)
    lines, record = trace_field_lines(field_T, unique, step_m=step_m, max_steps=max_steps,
                                      extent_m=extent_m)
    by_seed = {ln.seed_index: ln for ln in lines}
    votes = []
    for i in range(len(unique)):
        ln = by_seed.get(i)
        if ln is not None and ln.closed and len(votes) < 50:
            cyc = np.roll(ln.points_m, -ln.seed_vertex, axis=0)
            cyc = np.vstack([cyc, cyc[:1]])
            votes.append(np.sign(np.unwrap(angle_fn(cyc))[-1] - angle_fn(cyc[:1])[0]))
    sense = 1.0 if np.sum(votes) >= 0 else -1.0
    walls = {i: _Wall(by_seed[i], angle_fn, sense, phi_ref) for i in range(len(unique)) if i in by_seed}
    closed_walls = [w for w in walls.values() if w.closed]

    cell_walls, pos = [], 0
    for sc in side_chart:
        ids = []
        for s_ in sc:
            ids.append([int(owner[pos + t]) for t in range(len(s_[0]))])
            pos += len(s_[0])
        cell_walls.append(ids)

    kink = None
    if interface is not None:
        ic, ir = np.asarray(interface[0], float), float(interface[1])
        cr = {i: w.interface_angles(ir, ic, angle_fn, sense) for i, w in walls.items() if w.closed}
        two = {i: v for i, v in cr.items() if len(v) == 2}
        if not two:
            raise RuntimeError("no closed wall crosses the interface twice")
        step = 2 * np.pi / n_sec
        k_in = max(1, int(round(np.median([v[0] for v in two.values()]) / step)))
        Lu, Ld = k_in * step, 2 * np.pi - k_in * step
        for i, (pu, pd) in two.items():
            walls[i].kmap = (np.array([0.0, pu, pd, walls[i].phi[-1]]), np.array([0.0, Lu, Ld, 2 * np.pi]))
        kink = {"k_in": k_in, "walls_aligned": len(two), "closed_walls": len(cr),
                "crossings_count_hist": {str(n): sum(1 for v in cr.values() if len(v) == n)
                                         for n in sorted(set(len(v) for v in cr.values()))}}
    levels = 2 * np.pi * np.arange(n_sec) / n_sec

    level_shift = None
    if level_target_radius is not None and radial_centre is not None:
        ctr0 = np.asarray(radial_centre, float)
        fam0, closed0 = set(), set()
        for ids in cell_walls:
            ws_ = set(sum(ids, []))
            if all(i in walls for i in ws_):
                (closed0 if all(walls[i].closed for i in ws_) else fam0).update(ws_)
        K = n_sec // 2 - 1
        new_levels, level_shift = levels.copy(), {}
        for sgn, kidx in ((1, K), (-1, n_sec - K)):
            ph = []
            for i in sorted(fam0 & closed0):
                hit = walls[i].radial_point(level_target_radius, sgn, ctr0)
                if hit is None:
                    continue
                j = min(max(int(np.floor(hit[1])), 0), len(walls[i].phi) - 2)
                w_ = hit[1] - j
                p_raw = (1 - w_) * walls[i].phi[j] + w_ * walls[i].phi[j + 1]
                km = walls[i].kmap
                ph.append(np.interp(p_raw, km[0], km[1]) if km is not None else p_raw)
            if ph:
                cand = float(np.mean(ph))
                lo_ = levels[kidx - 1]
                hi_ = levels[(kidx + 1) % n_sec] + (2 * np.pi if kidx + 1 >= n_sec else 0)
                if lo_ < cand < hi_:
                    new_levels[kidx] = cand
                level_shift[str(kidx)] = [float(np.degrees(levels[kidx])),
                                          float(np.degrees(new_levels[kidx])), bool(lo_ < cand < hi_)]
        levels = new_levels

    def level(k):
        q, r = divmod(k, n_sec)
        return levels[r] + 2 * np.pi * q

    tubes = []
    for c, ids in enumerate(cell_walls):
        tw = sorted(set(sum(ids, [])))
        t = {"cell": c, "side_walls": ids, "walls": tw, "ok": all(i in walls for i in tw),
             "collapsed": len(set(ids[3])) == 1}
        if t["ok"]:
            if all(walls[i].closed for i in tw):
                t.update(closed=True, ks=list(range(n_sec)))
            else:
                hi = min(walls[i].hi if not walls[i].closed else np.inf for i in tw)
                lo = max(walls[i].lo if not walls[i].closed else -np.inf for i in tw)
                kf = max(k for k in range(0, 10 * n_sec) if level(k) <= hi + 1e-12) if hi >= 0 else 0
                kb = min(k for k in range(0, -10 * n_sec, -1) if level(k) >= lo - 1e-12) if lo <= 0 else 0
                t.update(closed=False, ks=list(range(kb, kf + 1)))
        tubes.append(t)

    radial = {"K_fwd": 0, "K_bwd": 0, "R_fwd": [], "R_bwd": []}
    fam = set()
    if radial_centre is not None:
        ctr = np.asarray(radial_centre, float)
        open_t = [t for t in tubes if t["ok"] and not t["closed"]]
        fam = set(i for t in open_t for i in t["walls"])
        bw = sorted(fam & set(i for t in tubes if t["ok"] and t["closed"] for i in t["walls"]))
        for sgn, name in ((1, "fwd"), (-1, "bwd")):
            Ks, Rs = 0, []
            for kk in range(1, n_sec // 2) if bw else ():
                k = sgn * kk
                Rk = float(np.mean([np.linalg.norm(walls[i].point_at([level(k)], angle_fn, sense)[0] - ctr)
                                    for i in bw]))
                if (Rs and Rk <= Rs[-1]) or (radial_max is not None and Rk > radial_max):
                    break
                ok = True
                for i in fam:
                    hit = walls[i].radial_point(Rk, sgn, ctr)
                    if hit is None:
                        ok = False
                        break
                    if walls[i].closed:
                        walls[i].point_at([level(k + sgn)], angle_fn, sense)
                        if (float(walls[i].last_frac_index[0]) - hit[1]) * sgn <= 0:
                            ok = False
                            break
                if not ok:
                    break
                Ks, Rs = kk, Rs + [Rk]
            radial["K_" + name], radial["R_" + name] = Ks, Rs
        for t in open_t:
            t["ks"] = list(range(-radial["K_bwd"], radial["K_fwd"] + 1))

    def signed(i, k):
        if not walls[i].closed:
            return k
        r = k % n_sec
        return r if r <= n_sec // 2 else r - n_sec

    def is_radial(i, k):
        sk = signed(i, k)
        return i in fam and sk != 0 and -radial["K_bwd"] <= sk <= radial["K_fwd"]

    def wkey(i, k):
        return (i, k % n_sec) if walls[i].closed else (i, k)

    need: dict = {}
    for t in tubes:
        if t["ok"]:
            for i in t["walls"]:
                need.setdefault(i, set()).update(t["ks"])
    wall_pt = {}
    for i, ks in need.items():
        ks = sorted({wkey(i, k)[1] for k in ks})
        ka = [k for k in ks if not is_radial(i, k)]
        if ka:
            for k, p in zip(ka, walls[i].point_at([level(k) for k in ka], angle_fn, sense)):
                wall_pt[(i, k)] = p
        for k in ks:
            if is_radial(i, k):
                sk = signed(i, k)
                R = radial["R_fwd"][sk - 1] if sk > 0 else radial["R_bwd"][-sk - 1]
                wall_pt[(i, k)] = walls[i].radial_point(R, 1 if sk > 0 else -1,
                                                        np.asarray(radial_centre, float))[0]

    node_id, coords = {}, []

    def nid(k_, p):
        if k_ not in node_id:
            node_id[k_] = len(coords)
            coords.append(p)
        return node_id[k_]

    spheres = []
    if interface is not None:
        spheres.append((np.asarray(interface[0], float), float(interface[1])))
    if radial_centre is not None:
        spheres += [(np.asarray(radial_centre, float), float(R)) for R in radial["R_fwd"] + radial["R_bwd"]]
    n_projected = 0
    elements, elem_tube, tube_stats = [], [], []
    q = abs(seeds.target_cell_flux_Wb)
    for t in tubes:
        c = t["cell"]
        if not t["ok"]:
            tube_stats.append({"cell": c, "why": "untraced wall"})
            continue
        trans = c in YF
        grids = []
        for k in t["ks"]:
            sides = [np.array([wall_pt[wkey(i, k)] for i in side]) for side in t["side_walls"]]
            blocks = _cell_blocks(sides, X[c], Y[c], t["collapsed"], trans)
            allw = np.vstack(sides)
            for ctr_, rad_ in spheres:    # a section on one sphere keeps its interior on it
                if np.all(np.abs(np.linalg.norm(allw - ctr_, axis=1) - rad_) <= 1e-9 * rad_):
                    blocks = [ctr_ + (g - ctr_) * (rad_ / np.linalg.norm(g - ctr_, axis=-1))[..., None]
                              for g in blocks]
                    n_projected += 1
                    break
            grids.append(blocks)
        labels = _block_labels(t["side_walls"], len(X[c]) - 1, len(Y[c]) - 1, t["collapsed"], trans)
        sec_ids = []
        for k, blocks in zip(t["ks"], grids):
            ids_k = []
            for bi, (grd, lab) in enumerate(zip(blocks, labels)):
                if grd.shape[:2] != lab.shape:
                    raise RuntimeError(f"tube {c}: block {bi} grid {grd.shape[:2]} vs labels {lab.shape}")
                ids_b = np.zeros(lab.shape, dtype=np.int64)
                for a_ in range(lab.shape[0]):
                    for b_ in range(lab.shape[1]):
                        lb = lab[a_, b_]
                        if lb is None:
                            ids_b[a_, b_] = nid(("t", c, k, bi, a_, b_), grd[a_, b_])
                        elif lb[0] == "w":
                            wk = wkey(lb[1], k)
                            ids_b[a_, b_] = nid(("w",) + wk, wall_pt[wk])
                        else:
                            ids_b[a_, b_] = nid(("s", c, k) + tuple(lb[1:]), grd[a_, b_])
                ids_k.append(ids_b)
            sec_ids.append(ids_k)
        n_layers = len(t["ks"]) if t["closed"] else len(t["ks"]) - 1
        hexes = []
        for L in range(n_layers):
            L1 = (L + 1) % len(t["ks"])
            for i0, i1 in zip(sec_ids[L], sec_ids[L1]):
                for a in range(i0.shape[0] - 1):
                    for bb in range(i0.shape[1] - 1):
                        qd = ((a, bb), (a + 1, bb), (a + 1, bb + 1), (a, bb + 1))
                        hexes.append([i0[x] for x in qd] + [i1[x] for x in qd])
        hexes = np.array(hexes, dtype=np.int64).reshape(-1, 8)
        C = np.asarray(coords)
        det = hex_corner_jacobians(C, hexes)
        if np.sum(det < 0) > np.sum(det > 0):
            hexes = np.hstack([hexes[:, 4:], hexes[:, :4]])
            det = -det
        scj, er = hex_quality(C, hexes)
        flux = np.array([sum(_grid_flux(field_T, g) for g in blocks) for blocks in grids])
        ratio = flux * np.sign(flux[0]) / q
        stat = {"cell": c, "closed": t["closed"], "k_range": [t["ks"][0], t["ks"][-1]],
                "n_hex": int(len(hexes)), "grid_nx_ny": [len(X[c]) - 1, len(Y[c]) - 1],
                "collapsed": t["collapsed"], "transition": trans,
                "invalid_hex": int(np.sum(det.min(axis=1) <= 0)),
                "min_scaled_corner_J": float(scj.min()), "max_edge_ratio": float(er.max()),
                "section_flux_Wb": flux.tolist(),
                "flux_ratio_min_max": [float(ratio.min()), float(ratio.max())]}
        if trans:
            stat["fine_outer_segments"] = int(len(YF[c]) - 1)
        tube_stats.append(stat)
        elements.append(hexes)
        elem_tube += [c] * len(hexes)
    if not elements:
        raise RuntimeError("no tube could be built: every tube has an untraced wall")
    E = np.vstack(elements)
    tube_of = np.array(elem_tube, dtype=np.int64)
    order = np.argsort(tube_of, kind="stable")
    info = {"status": STATUS, "n_walls_traced": int(len(unique)), "sense": sense,
            "closed_walls": len(closed_walls),
            "closed_walls_not_one_turn": sum(1 for w in closed_walls if abs(w.turn - 2 * np.pi) > 1e-3),
            "max_phi_decrease_along_wall_rad": max(w.monotone_violation for w in walls.values()),
            "levels_deg": np.degrees(levels).tolist(), "n_sections": n_sec,
            "radial_open_family": {"n_walls": len(fam), **radial},
            "interface_alignment": kink, "level_shift_deg": level_shift,
            "transition_cells": sorted(YF), "sections_projected_to_sphere": n_projected,
            "target_cell_flux_Wb": float(seeds.target_cell_flux_Wb),
            "grid_nx_ny_per_cell": [[len(x) - 1, len(y) - 1] for x, y in zip(X, Y)],
            "build_s": time.perf_counter() - t_start, "trace_record": {"trace_s": record.get("trace_s")}}
    return FluxTubeMesh(np.asarray(coords), E[order], tube_of[order], tube_stats, info)


def flux_tube_report(mesh: FluxTubeMesh) -> dict:
    """Quality, Jacobian bound and face census of a built mesh (nothing is asserted)."""
    scj, er = hex_quality(mesh.nodes_m, mesh.hexes)
    det = hex_corner_jacobians(mesh.nodes_m, mesh.hexes)
    lb = hex_jacobian_lower_bound(mesh.nodes_m, mesh.hexes)
    census = face_census(mesh.hexes, mesh.tube_of)
    census.pop("shared_face_pairs")
    built = [s for s in mesh.tubes if "n_hex" in s]
    return {"status": mesh.status, "n_hex": int(len(mesh.hexes)), "n_nodes": int(len(mesh.nodes_m)),
            "invalid_corner_hex": int(np.sum(det.min(axis=1) <= 0)),
            "nonpositive_jacobian_bound_hex": int(np.sum(lb <= 0)),
            "min_jacobian_bound_m3": float(lb.min()),
            "min_scaled_corner_J": float(scj.min()), "max_edge_ratio": float(er.max()),
            "flux_ratio_min_max": [min(s["flux_ratio_min_max"][0] for s in built),
                                   max(s["flux_ratio_min_max"][1] for s in built)],
            "untraced_tubes": len(mesh.tubes) - len(built), "faces": census}


def write_flux_tube_msh(path, mesh: FluxTubeMesh, length_scale: float = 1.0) -> None:
    """Write the tubes as a Gmsh MSH 4.2 file: one volume entity per tube, physical group
    "flux_tubes", element data "tube" (the seed cell). Coordinates are divided by
    length_scale (1e-3 writes millimetres)."""
    from .gmsh_post_export import (DataField, ElementBlock, Entity, NodeBlock, PayloadV42,
                                   PhysicalGroup, _emit_v42_msh)
    payload = PayloadV42(groups=[PhysicalGroup(3, 1, "flux_tubes")])
    pts = np.asarray(mesh.nodes_m, float) / float(length_scale)
    lo, hi = pts.min(axis=0), pts.max(axis=0)
    tags = sorted(set(int(t) for t in mesh.tube_of))
    for tg in tags:
        payload.entities.append(Entity(3, tg + 1, (*map(float, lo), *map(float, hi)), [1]))
    payload.node_blocks.append(NodeBlock(3, tags[0] + 1, 1, [tuple(map(float, p)) for p in pts]))
    tube_sorted = []
    for tg in tags:
        sel = mesh.hexes[mesh.tube_of == tg] + 1
        payload.element_blocks.append(ElementBlock(3, tg + 1, 5, [[int(x) for x in h] for h in sel]))
        tube_sorted += [tg] * len(sel)
    payload.data_fields.append(DataField("$ElementData", "tube", 1,
                                         [(i + 1, float(t)) for i, t in enumerate(tube_sorted)]))
    with open(path, "w", encoding="utf-8", newline="\n") as stream:
        _emit_v42_msh(stream, payload)
