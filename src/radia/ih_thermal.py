"""Field artifacts and surface heat-source transfer for induction heating.

The electromagnetic step writes the surface loss density ``q_surf`` [W/m^2]
as an NGSolve ``GridFunction.Save`` file.  That format is a bare coefficient
vector: it carries no mesh, no order and no units, and ``GridFunction.Load``
does not check its length.  Pairing it with the wrong mesh, the wrong order,
or a mesh that was modified after loading (for example by a failed
``Mesh.Curve`` call) silently produces nonsense.  This module owns the
contract that prevents that:

* :func:`write_field_sidecar` / :func:`load_field` -- a JSON sidecar
  ``<file>.sol.json`` records the mesh file hash, the H1 order, the number of
  DOFs, the physical quantity and, for heat sources, the electromagnetic
  power.  :func:`load_field` rebuilds the GridFunction exactly as it was
  written and never curves the mesh.
* :class:`SurfaceP1Field` -- a P1 boundary field that can be evaluated at
  arbitrary points by projection onto the nearest source triangle.  Every
  target point must lie within an explicit distance of the source surface;
  there is no zero-flux fallback.
* :class:`AxisymmetricSurfaceProfile` -- the exact circumferential average of
  a surface field on a body of revolution.  The average is integrated over
  the source triangles (not sampled at a few angles) and binned along the
  meridian arc length, so it conserves the source power by construction and
  the value at a target point depends only on its meridian position.

Only NumPy, SciPy and NGSolve are used, so the same code serves the 3D and
axisymmetric thermal solvers, the Simulink operator assembly and user
post-processing.
"""
from __future__ import annotations

import hashlib
import json
import math
import os
from dataclasses import dataclass, field
from typing import Iterable, Sequence

import numpy as np

SIDECAR_SCHEMA = "radia.field-artifact/1"
SIDECAR_SUFFIX = ".json"

# Irrational section angle: CAD-derived meshes often place a seam and
# vertices exactly on the phi = 0 half-plane, where edge/plane tests become
# degenerate.
_SECTION_ANGLE_RAD = 0.3717028531


# ---------------------------------------------------------------------------
# Artifacts
# ---------------------------------------------------------------------------

def file_sha256(path: str, chunk: int = 1 << 20) -> str:
    """Return the hex SHA-256 digest of ``path``."""
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        while True:
            block = handle.read(chunk)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def _h1(mesh, order, definedon=None):
    from ngsolve import H1
    if definedon:
        return H1(mesh, order=int(order), definedon=definedon)
    return H1(mesh, order=int(order))


def sidecar_path(sol_path: str) -> str:
    """Return the sidecar path that belongs to ``sol_path``."""
    return os.fspath(sol_path) + SIDECAR_SUFFIX


def write_field_sidecar(sol_path: str, *, mesh_path: str, mesh, fes_order: int,
                        quantity: str, unit: str, space: str = "H1",
                        boundaries: Sequence[str] | None = None,
                        extra: dict | None = None,
                        definedon: str | None = None) -> str:
    """Write ``<sol_path>.json`` describing a saved H1 GridFunction.

    ``mesh`` is the NGSolve mesh the field was saved on and ``mesh_path`` the
    ``.vol`` file it was loaded from (or written to).  The digest ties the
    coefficient vector to that exact file.
    """
    from ngsolve import H1

    if space != "H1":
        raise ValueError(f"only H1 field artifacts are supported, got {space!r}")
    sol_path = os.path.abspath(sol_path)
    mesh_path = os.path.abspath(mesh_path)
    if not os.path.isfile(sol_path):
        raise FileNotFoundError(f"field file not found: {sol_path}")
    if not os.path.isfile(mesh_path):
        raise FileNotFoundError(f"mesh file not found: {mesh_path}")
    ndof = int(_h1(mesh, fes_order, definedon).ndof)
    size = os.path.getsize(sol_path)
    if size != 8 * ndof:
        raise ValueError(
            f"{os.path.basename(sol_path)} holds {size} bytes but H1 order "
            f"{fes_order} on this mesh has {ndof} DOFs ({8 * ndof} bytes)")
    record = {
        "schema": SIDECAR_SCHEMA,
        "quantity": str(quantity),
        "unit": str(unit),
        "space": space,
        "fes_order": int(fes_order),
        "ndof": ndof,
        "sol_file": os.path.basename(sol_path),
        "sol_sha256": file_sha256(sol_path),
        "mesh_file": mesh_path.replace("\\", "/"),
        "mesh_sha256": file_sha256(mesh_path),
        "mesh_nv": int(mesh.nv),
        "mesh_ne": int(mesh.ne),
        "mesh_dim": int(mesh.dim),
        "mesh_curve_order": int(mesh.GetCurveOrder()),
        "boundaries": list(boundaries or []),
        "definedon": definedon,
    }
    if extra:
        overlap = set(extra) & set(record)
        if overlap:
            raise ValueError(f"sidecar extra keys collide with {sorted(overlap)}")
        record.update(extra)
    out = sidecar_path(sol_path)
    with open(out, "w", encoding="utf-8") as handle:
        json.dump(record, handle, indent=1, ensure_ascii=False)
    return out


def read_field_sidecar(sol_path: str) -> dict | None:
    """Return the sidecar record for ``sol_path`` or ``None`` when absent."""
    path = sidecar_path(os.path.abspath(sol_path))
    if not os.path.isfile(path):
        return None
    with open(path, encoding="utf-8") as handle:
        record = json.load(handle)
    if record.get("schema") != SIDECAR_SCHEMA:
        raise ValueError(
            f"{os.path.basename(path)} has schema {record.get('schema')!r}; "
            f"expected {SIDECAR_SCHEMA!r}")
    return record


def verify_field_pair(sol_path: str, mesh_path: str, mesh, fes_order: int,
                      *, quantity: str) -> dict:
    """Check that ``sol_path`` really belongs to ``mesh``/``fes_order``.

    The sidecar is required: the mesh digest, order, DOF count and quantity
    must all match.  A legacy ``.sol`` without one gets a sidecar through an
    explicit ``python -m radia.ih_thermal sidecar ...`` in which the user
    states the pairing; nothing is inferred here.
    """
    sol_path = os.path.abspath(sol_path)
    record = read_field_sidecar(sol_path)
    if record is None:
        raise ValueError(
            f"{os.path.basename(sol_path)} has no {SIDECAR_SUFFIX} sidecar, so "
            "its mesh, order and quantity cannot be verified. Write one with "
            "'python -m radia.ih_thermal sidecar --sol ... --mesh ... "
            "--order ... --quantity ...' after confirming the pairing.")
    definedon = record.get("definedon")
    ndof = int(_h1(mesh, fes_order, definedon).ndof)
    size = os.path.getsize(sol_path)
    if size != 8 * ndof:
        raise ValueError(
            f"{os.path.basename(sol_path)} holds {size} bytes but H1 order "
            f"{fes_order} on {os.path.basename(mesh_path)} has {ndof} DOFs "
            f"({8 * ndof} bytes); the field and mesh do not belong together")
    problems = []
    if int(record["fes_order"]) != int(fes_order):
        problems.append(f"order {record['fes_order']} != {fes_order}")
    if int(record["ndof"]) != ndof:
        problems.append(f"ndof {record['ndof']} != {ndof}")
    mesh_digest = file_sha256(mesh_path)
    if record["mesh_sha256"] != mesh_digest:
        problems.append("mesh file digest differs from the one the field "
                        f"was saved on ({record['mesh_file']})")
    if record["sol_sha256"] != file_sha256(sol_path):
        problems.append("field file changed after its sidecar was written")
    if record["quantity"] != quantity:
        problems.append(f"quantity {record['quantity']!r} != {quantity!r}")
    if problems:
        raise ValueError(
            f"{os.path.basename(sol_path)} does not match its sidecar: "
            + "; ".join(problems))
    return {"provenance": "sidecar-verified", "ndof": ndof,
            "sidecar": sidecar_path(sol_path),
            "quantity": record["quantity"]}


def load_field(sol_path: str, *, quantity: str, mesh_path: str | None = None):
    """Load a saved H1 field exactly as it was written.

    The sidecar names the mesh, the order and the ``definedon`` region.
    ``mesh_path`` relocates the mesh file only; its digest must still match.
    The mesh is used as loaded: its stored geometry order is preserved and
    ``Mesh.Curve`` is never called (a failed ``Curve`` on a CAD-derived
    ``.vol`` corrupts the element maps while leaving the mesh usable).

    Returns ``(mesh, gridfunction, audit)``; ``audit["sidecar_record"]`` is
    the sidecar content (for ``definedon``, units, producer).
    """
    from ngsolve import GridFunction, Mesh

    record = read_field_sidecar(sol_path)
    if record is None:
        raise ValueError(
            f"{os.path.basename(sol_path)} has no {SIDECAR_SUFFIX} sidecar; "
            "write one with 'python -m radia.ih_thermal sidecar ...' after "
            "confirming the pairing")
    if mesh_path is None:
        mesh_path = record["mesh_file"]
    fes_order = int(record["fes_order"])
    mesh = Mesh(mesh_path)
    audit = verify_field_pair(sol_path, mesh_path, mesh, fes_order,
                              quantity=quantity)
    gf = GridFunction(_h1(mesh, fes_order, record.get("definedon")))
    gf.Load(os.path.abspath(sol_path))
    audit.update({"mesh_file": os.path.abspath(mesh_path),
                  "fes_order": int(fes_order),
                  "mesh_curve_order": int(mesh.GetCurveOrder()),
                  "sidecar_record": record})
    return mesh, gf, audit


# ---------------------------------------------------------------------------
# Boundary geometry
# ---------------------------------------------------------------------------

def mesh_vertices(mesh) -> np.ndarray:
    """Return vertex coordinates as an ``(nv, 3)`` array (2D meshes get z=0)."""
    pts = np.array([tuple(v.point) for v in mesh.vertices], dtype=float)
    if pts.shape[1] == 2:
        pts = np.c_[pts, np.zeros(len(pts))]
    return pts


def boundary_triangles(mesh, boundary_names: Iterable[str]) -> np.ndarray:
    """Vertex numbers of the triangles on the named boundaries, ``(m, 3)``."""
    from ngsolve import BND

    names = set(boundary_names)
    tris = []
    for el in mesh.Elements(BND):
        if el.mat not in names:
            continue
        nrs = [v.nr for v in el.vertices]
        if len(nrs) == 3:
            tris.append(nrs)
        elif len(nrs) == 4:                     # quad face -> two triangles
            tris.append([nrs[0], nrs[1], nrs[2]])
            tris.append([nrs[0], nrs[2], nrs[3]])
        else:
            raise ValueError(f"unsupported boundary element with {len(nrs)} "
                             "vertices")
    if not tris:
        raise ValueError(f"no boundary triangles on {sorted(names)}")
    return np.asarray(tris, dtype=int)


def boundary_vertex_numbers(mesh, boundary_names: Iterable[str]) -> np.ndarray:
    """Sorted vertex numbers touching the named boundaries."""
    from ngsolve import BND

    names = set(boundary_names)
    out = set()
    for el in mesh.Elements(BND):
        if el.mat in names:
            out.update(v.nr for v in el.vertices)
    if not out:
        raise ValueError(f"no boundary vertices on {sorted(names)}")
    return np.asarray(sorted(out), dtype=int)


def boundary_edges_2d(mesh, boundary_names: Iterable[str]) -> np.ndarray:
    """Vertex numbers of the segments on named boundaries of a 2D mesh."""
    from ngsolve import BND

    names = set(boundary_names)
    segs = [[v.nr for v in el.vertices] for el in mesh.Elements(BND)
            if el.mat in names]
    if not segs:
        raise ValueError(f"no boundary segments on {sorted(names)}")
    return np.asarray(segs, dtype=int)


def median_edge_length(points: np.ndarray, tris: np.ndarray) -> float:
    """Median triangle edge length."""
    e = np.concatenate([
        np.linalg.norm(points[tris[:, 1]] - points[tris[:, 0]], axis=1),
        np.linalg.norm(points[tris[:, 2]] - points[tris[:, 1]], axis=1),
        np.linalg.norm(points[tris[:, 0]] - points[tris[:, 2]], axis=1)])
    return float(np.median(e))


def triangle_areas(points: np.ndarray, tris: np.ndarray) -> np.ndarray:
    a, b, c = points[tris[:, 0]], points[tris[:, 1]], points[tris[:, 2]]
    return 0.5 * np.linalg.norm(np.cross(b - a, c - a), axis=1)


def p1_surface_power(points: np.ndarray, tris: np.ndarray,
                     q_nodal: np.ndarray) -> float:
    """Exact integral of a P1 field over flat triangles."""
    area = triangle_areas(points, tris)
    return float(np.sum(area * q_nodal[tris].mean(axis=1)))


# ---------------------------------------------------------------------------
# Rotation about a coordinate axis
# ---------------------------------------------------------------------------

_AXES = {"x": (1, 2, 0), "y": (2, 0, 1), "z": (0, 1, 2)}


def _axis_frame(axis: str):
    """Return index triple (u, w, a): (u, w) span the plane normal to axis a,
    ordered so rotation by +theta is right-handed about +a."""
    key = str(axis).lower().strip()
    if key not in _AXES:
        raise ValueError(f"rotation axis must be x, y or z, got {axis!r}")
    return _AXES[key]


def rotate_about_axis(points: np.ndarray, theta: float, axis: str = "z"):
    """Rotate ``points`` by ``theta`` [rad] about the coordinate ``axis``."""
    iu, iw, _ = _axis_frame(axis)
    c, s = math.cos(theta), math.sin(theta)
    out = np.array(points, dtype=float, copy=True)
    u, w = points[:, iu], points[:, iw]
    out[:, iu] = c * u - s * w
    out[:, iw] = s * u + c * w
    return out


def to_meridian(points: np.ndarray, axis: str = "z") -> np.ndarray:
    """Map 3D points to meridian coordinates ``(r, a)``."""
    iu, iw, ia = _axis_frame(axis)
    return np.c_[np.hypot(points[:, iu], points[:, iw]), points[:, ia]]


# ---------------------------------------------------------------------------
# Point-to-triangle projection
# ---------------------------------------------------------------------------

def _closest_point_on_triangles(p, a, b, c):
    """Closest points on triangles (a, b, c) to points p (all (n, 3)).

    Returns (closest, barycentric (n, 3)).  Ericson, *Real-Time Collision
    Detection*, section 5.1.5, vectorised.
    """
    ab, ac, ap = b - a, c - a, p - a
    d1 = np.einsum("ij,ij->i", ab, ap)
    d2 = np.einsum("ij,ij->i", ac, ap)
    bp = p - b
    d3 = np.einsum("ij,ij->i", ab, bp)
    d4 = np.einsum("ij,ij->i", ac, bp)
    cp = p - c
    d5 = np.einsum("ij,ij->i", ab, cp)
    d6 = np.einsum("ij,ij->i", ac, cp)
    va = d3 * d6 - d5 * d4
    vb = d5 * d2 - d1 * d6
    vc = d1 * d4 - d3 * d2
    n = len(p)
    bary = np.zeros((n, 3))
    done = np.zeros(n, dtype=bool)

    def _set(mask, l0, l1, l2):
        m = mask & ~done
        bary[m, 0], bary[m, 1], bary[m, 2] = l0[m], l1[m], l2[m]
        done[m] = True

    one, zero = np.ones(n), np.zeros(n)
    with np.errstate(divide="ignore", invalid="ignore"):
        _set((d1 <= 0) & (d2 <= 0), one, zero, zero)                 # A
        _set((d3 >= 0) & (d4 <= d3), zero, one, zero)                # B
        v = d1 / (d1 - d3)
        _set((vc <= 0) & (d1 >= 0) & (d3 <= 0), 1 - v, v, zero)     # AB
        _set((d6 >= 0) & (d5 <= d6), zero, zero, one)                # C
        w = d2 / (d2 - d6)
        _set((vb <= 0) & (d2 >= 0) & (d6 <= 0), 1 - w, zero, w)     # AC
        w = (d4 - d3) / ((d4 - d3) + (d5 - d6))
        _set((va <= 0) & ((d4 - d3) >= 0) & ((d5 - d6) >= 0),
             zero, 1 - w, w)                                          # BC
        denom = 1.0 / (va + vb + vc)
        v = vb * denom
        w = vc * denom
        _set(np.ones(n, dtype=bool), 1 - v - w, v, w)                # face
    closest = (bary[:, :1] * a + bary[:, 1:2] * b + bary[:, 2:] * c)
    return closest, bary


class SurfaceP1Field:
    """A P1 field on a triangulated surface, evaluable anywhere nearby.

    ``evaluate`` projects each query point onto the closest source triangle
    and interpolates there.  Query points farther than ``max_distance`` from
    the source surface raise ``ValueError``: the thermal and electromagnetic
    meshes must describe the same physical surface.
    """

    def __init__(self, points: np.ndarray, tris: np.ndarray,
                 values: np.ndarray, *, candidates: int = 12):
        from scipy.spatial import cKDTree

        self.points = np.asarray(points, float)
        self.tris = np.asarray(tris, int)
        self.values = np.asarray(values, float)
        a = self.points[self.tris[:, 0]]
        b = self.points[self.tris[:, 1]]
        c = self.points[self.tris[:, 2]]
        self._abc = (a, b, c)
        self._tree = cKDTree((a + b + c) / 3.0)
        self._k = int(min(candidates, len(self.tris)))
        self.h_median = median_edge_length(self.points, self.tris)

    @classmethod
    def from_gridfunction(cls, mesh, gf, boundary_names):
        """Build from an H1 order-1 GridFunction on the named boundaries."""
        if gf.space.globalorder != 1:
            raise ValueError("SurfaceP1Field needs an H1 order-1 field")
        vals = np.asarray(gf.vec.FV().NumPy(), float)[:mesh.nv]
        return cls(mesh_vertices(mesh), boundary_triangles(mesh, boundary_names),
                   vals)

    @property
    def power(self) -> float:
        """Integral of the field over the (flat) source triangles."""
        return p1_surface_power(self.points, self.tris, self.values)

    @property
    def magnitude_power(self) -> float:
        """Integral of the P1 interpolant of ``|q|`` (normalisation scale)."""
        return p1_surface_power(self.points, self.tris, np.abs(self.values))

    def locate(self, query: np.ndarray):
        """Return (triangle index, barycentric, distance) for each query."""
        q = np.asarray(query, float)
        _, cand = self._tree.query(q, k=self._k)
        cand = np.atleast_2d(cand).reshape(len(q), -1)
        a, b, c = self._abc
        best_d = np.full(len(q), np.inf)
        best_t = np.zeros(len(q), dtype=int)
        best_l = np.zeros((len(q), 3))
        for j in range(cand.shape[1]):
            t = cand[:, j]
            cp, lam = _closest_point_on_triangles(q, a[t], b[t], c[t])
            d = np.linalg.norm(q - cp, axis=1)
            better = d < best_d
            best_d[better] = d[better]
            best_t[better] = t[better]
            best_l[better] = lam[better]
        return best_t, best_l, best_d

    def evaluate(self, query: np.ndarray, *, max_distance: float,
                 what: str = "target points"):
        """Interpolate at ``query``; returns (values, distances)."""
        t, lam, d = self.locate(query)
        bad = d > max_distance
        if np.any(bad):
            worst = np.argsort(d)[::-1][:5]
            raise ValueError(
                f"{int(bad.sum())}/{len(d)} {what} lie farther than "
                f"{max_distance:.3e} m from the source surface (worst "
                f"{d[worst[0]]:.3e} m at {np.round(query[worst[0]], 6).tolist()}). "
                "The thermal and electromagnetic meshes must describe the "
                "same physical surface; no zero-flux fallback is applied.")
        v = np.einsum("ij,ij->i", lam, self.values[self.tris[t]])
        return v, d


# ---------------------------------------------------------------------------
# Exact circumferential average on a body of revolution
# ---------------------------------------------------------------------------

@dataclass
class _Chain:
    rz: np.ndarray          # (n, 2) meridian polyline
    s: np.ndarray           # (n,) arc length
    closed: bool
    centres: np.ndarray = field(default=None)   # bin centres (arc length)
    q: np.ndarray = field(default=None)         # bin-averaged value
    power: float = 0.0


def _cumulative_below(level, s, q, area):
    """Integrals over the part of linear triangles where s < level.

    ``s`` and ``q`` are ``(n, 3)`` vertex values sorted by ``s`` and
    ``area`` the triangle areas; ``level`` is ``(n,)``.  Both ``s`` and
    ``q`` are linear on each triangle, so the cut regions are triangles
    whose integrals are exact.  Returns ``(area_below, q_below)``.
    """
    sa, sb, sc = s[:, 0], s[:, 1], s[:, 2]
    qa, qb, qc = q[:, 0], q[:, 1], q[:, 2]
    q_total = area * (qa + qb + qc) / 3.0
    a_out = np.zeros_like(level)
    q_out = np.zeros_like(level)
    with np.errstate(divide="ignore", invalid="ignore"):
        low = (level > sa) & (level <= sb) & (sb > sa)
        t_ab = np.where(low, (level - sa) / (sb - sa), 0.0)
        t_ac = np.where(low, (level - sa) / (sc - sa), 0.0)
        a_low = area * t_ab * t_ac
        q_low = a_low * (3 * qa + t_ab * (qb - qa) + t_ac * (qc - qa)) / 3.0
        high = (level > sb) & (level < sc)
        u_cb = np.where(high, (sc - level) / (sc - sb), 0.0)
        u_ca = np.where(high, (sc - level) / (sc - sa), 0.0)
        a_up = area * u_cb * u_ca
        q_up = a_up * (3 * qc + u_cb * (qb - qc) + u_ca * (qa - qc)) / 3.0
    a_out = np.where(low, a_low, a_out)
    q_out = np.where(low, q_low, q_out)
    a_out = np.where(high, area - a_up, a_out)
    q_out = np.where(high, q_total - q_up, q_out)
    full = level >= sc
    a_out = np.where(full, area, a_out)
    q_out = np.where(full, q_total, q_out)
    return a_out, q_out


def _split_at_axis(pts, tris, vals, axis):
    """Split every triangle pierced by ``axis`` into three at the piercing
    point, with the P1 value there.  Returns new ``(pts, tris, vals)``."""
    iu, iw, _ia = _axis_frame(axis)
    uw = pts[:, [iu, iw]][tris]                          # (n, 3, 2)
    a, b, c = uw[:, 0], uw[:, 1], uw[:, 2]
    det = (b[:, 0] - a[:, 0]) * (c[:, 1] - a[:, 1]) - \
        (c[:, 0] - a[:, 0]) * (b[:, 1] - a[:, 1])
    scale = np.max(np.abs(uw).reshape(len(tris), -1), axis=1)
    live = np.abs(det) > 1e-12 * np.maximum(scale, 1e-300) ** 2
    d = np.where(live, det, 1.0)          # tangent-to-axis triangles: no hit
    l1 = (a[:, 0] * (a[:, 1] - c[:, 1]) - a[:, 1] * (a[:, 0] - c[:, 0])) / d
    l2 = (a[:, 1] * (a[:, 0] - b[:, 0]) - a[:, 0] * (a[:, 1] - b[:, 1])) / d
    lam = np.stack([1.0 - l1 - l2, l1, l2], axis=1)
    # interior piercing only: at a vertex the kink is already a vertex, and on
    # an edge the neighbouring triangle is split instead
    eps = 1e-9
    hit = live & np.all(lam > eps, axis=1)
    if not np.any(hit):
        return pts, tris, vals
    idx = np.flatnonzero(hit)
    lam_h = lam[idx]
    new_pts = np.einsum("ij,ijk->ik", lam_h, pts[tris[idx]])
    new_vals = np.einsum("ij,ij->i", lam_h, vals[tris[idx]])
    base = len(pts)
    p_ids = base + np.arange(len(idx))
    t = tris[idx]
    sub = np.concatenate([np.c_[p_ids, t[:, 0], t[:, 1]],
                          np.c_[p_ids, t[:, 1], t[:, 2]],
                          np.c_[p_ids, t[:, 2], t[:, 0]]])
    keep = np.delete(tris, idx, axis=0)
    return (np.vstack([pts, new_pts]), np.vstack([keep, sub]),
            np.concatenate([vals, new_vals]))


class AxisymmetricSurfaceProfile:
    """Circumferential average of a P1 surface field about a coordinate axis.

    Every source vertex is given its meridian arc length ``s``; on each
    triangle ``s`` is taken as linear, and the triangle is cut exactly at
    the bin edges ``s = const``.  The bin averages are therefore exact
    integrals of the P1 source, with no sampling in azimuth, and the
    profile power equals the source power to rounding.  The meridian itself
    is the section of the source surface by a half-plane at a fixed azimuth.

    ``evaluate(rz)`` returns the profile at meridian points and raises when a
    point is farther than ``max_distance`` from the meridian.  Source
    vertices are checked the same way, which rejects surfaces that are not
    bodies of revolution about ``axis``.
    """

    def __init__(self, points: np.ndarray, tris: np.ndarray,
                 values: np.ndarray, *, axis: str = "z",
                 bin_width: float | None = None,
                 axisymmetry_tolerance: float | None = None):
        from scipy.spatial import cKDTree

        self.axis = str(axis).lower()
        _axis_frame(self.axis)
        pts = np.asarray(points, float)
        tris = np.asarray(tris, int)
        vals = np.asarray(values, float)
        h = median_edge_length(pts, tris)
        self.h_median = h
        requested = float(bin_width) if bin_width else h / 4.0
        if requested <= 0:
            raise ValueError("bin_width must be positive")
        self.axisymmetry_tolerance = (float(axisymmetry_tolerance)
                                      if axisymmetry_tolerance else 0.25 * h)

        self.chains = self._section_chains(pts, tris)
        dense_rz, dense_chain, dense_s = [], [], []
        step = min(requested, h) / 20.0
        for k, ch in enumerate(self.chains):
            n = max(2, int(math.ceil(ch.s[-1] / step)) + 1)
            s = np.linspace(0.0, ch.s[-1], n)
            dense_rz.append(np.c_[np.interp(s, ch.s, ch.rz[:, 0]),
                                  np.interp(s, ch.s, ch.rz[:, 1])])
            dense_chain.append(np.full(n, k))
            dense_s.append(s)
        self._dense_rz = np.concatenate(dense_rz)
        self._dense_chain = np.concatenate(dense_chain)
        self._dense_s = np.concatenate(dense_s)
        self._tree = cKDTree(self._dense_rz)

        # The meridian coordinate is not linear on a triangle the axis
        # pierces (it has a kink at r = 0).  Such triangles are split at the
        # piercing point; q is P1, so the split is exact for q.
        pts, tris, vals = _split_at_axis(pts, tris, vals, self.axis)

        # Meridian position of every source vertex.
        used = np.unique(tris)
        rz_v = to_meridian(pts[used], self.axis)
        d_v, idx_v = self._tree.query(rz_v)
        self.max_sample_distance = float(d_v.max())
        if self.max_sample_distance > self.axisymmetry_tolerance:
            worst = int(np.argmax(d_v))
            raise ValueError(
                "the source surface is not a body of revolution about the "
                f"{self.axis} axis: source vertex (r, a) = "
                f"{np.round(rz_v[worst], 6).tolist()} lies "
                f"{self.max_sample_distance:.3e} m from the meridian section "
                f"(tolerance {self.axisymmetry_tolerance:.3e} m)")
        chain_of = np.full(len(pts), -1)
        s_of = np.zeros(len(pts))
        chain_of[used] = self._dense_chain[idx_v]
        s_of[used] = self._dense_s[idx_v]

        tri_chain = chain_of[tris]
        if np.any(tri_chain != tri_chain[:, :1]):
            raise ValueError("a source triangle spans two meridian sections; "
                             "the heated surfaces touch or overlap")
        tri_chain = tri_chain[:, 0]
        area = triangle_areas(pts, tris)
        self.bin_width = requested
        self._widths = []
        for k, ch in enumerate(self.chains):
            m = tri_chain == k
            length = float(ch.s[-1])
            nb = max(1, int(round(length / requested)))
            w = length / nb
            self._widths.append(w)
            num = np.zeros(nb)
            den = np.zeros(nb)
            if np.any(m):
                s_t = s_of[tris[m]]
                q_t = vals[tris[m]]
                if ch.closed:          # unwrap across the seam
                    ref = s_t[:, :1]
                    s_t = ref + np.mod(s_t - ref + 0.5 * length, length) \
                        - 0.5 * length
                order = np.argsort(s_t, axis=1)
                s_t = np.take_along_axis(s_t, order, axis=1)
                q_t = np.take_along_axis(q_t, order, axis=1)
                a_t = area[m]
                k0 = np.floor(s_t[:, 0] / w).astype(int)
                k1 = np.floor(s_t[:, 2] / w).astype(int)
                k1 = np.maximum(k1, k0)
                count = k1 - k0 + 1
                tri_idx = np.repeat(np.arange(len(a_t)), count)
                offset = np.arange(count.sum()) - np.repeat(
                    np.cumsum(count) - count, count)
                bin_k = k0[tri_idx] + offset
                upper = (bin_k + 1) * w
                lower = bin_k * w
                a_hi, q_hi = _cumulative_below(
                    upper, s_t[tri_idx], q_t[tri_idx], a_t[tri_idx])
                a_lo, q_lo = _cumulative_below(
                    lower, s_t[tri_idx], q_t[tri_idx], a_t[tri_idx])
                # Degenerate triangles with constant s fall in one bin.
                flat = s_t[tri_idx, 2] <= s_t[tri_idx, 0]
                q_tot = a_t[tri_idx] * q_t[tri_idx].mean(axis=1)
                a_bin = np.where(flat, a_t[tri_idx], a_hi - a_lo)
                q_bin = np.where(flat, q_tot, q_hi - q_lo)
                if ch.closed:
                    bin_k = np.mod(bin_k, nb)
                else:
                    bin_k = np.clip(bin_k, 0, nb - 1)
                num = np.bincount(bin_k, weights=q_bin, minlength=nb)
                den = np.bincount(bin_k, weights=a_bin, minlength=nb)
            empty = np.flatnonzero(den <= 0)
            if empty.size:
                raise ValueError(
                    f"meridian section {k} has {empty.size}/{nb} bins without "
                    "source area (bin width "
                    f"{w:.3e} m); the section and the source triangles "
                    "disagree")
            ch.centres = (np.arange(nb) + 0.5) * w
            ch.q = num / den
            ch.power = float(num.sum())
        self.source_power = p1_surface_power(pts, tris, vals)
        self.power = float(sum(ch.power for ch in self.chains))

    @classmethod
    def from_gridfunction(cls, mesh, gf, boundary_names, **kwargs):
        if gf.space.globalorder != 1:
            raise ValueError("AxisymmetricSurfaceProfile needs an H1 order-1 "
                             "field")
        vals = np.asarray(gf.vec.FV().NumPy(), float)[:mesh.nv]
        return cls(mesh_vertices(mesh), boundary_triangles(mesh, boundary_names),
                   vals, **kwargs)

    # -- meridian section ---------------------------------------------------
    def _section_chains(self, pts, tris):
        iu, iw, ia = _axis_frame(self.axis)
        c0, s0 = math.cos(_SECTION_ANGLE_RAD), math.sin(_SECTION_ANGLE_RAD)
        side = -s0 * pts[:, iu] + c0 * pts[:, iw]      # signed distance
        along = c0 * pts[:, iu] + s0 * pts[:, iw]      # radial direction
        side = np.where(side == 0.0, 1e-300, side)
        points = {}
        adjacency: dict = {}
        for ti, tri in enumerate(tris):
            hits = []
            for e0, e1 in ((tri[0], tri[1]), (tri[1], tri[2]),
                           (tri[2], tri[0])):
                s_a, s_b = side[e0], side[e1]
                if (s_a > 0) == (s_b > 0):
                    continue
                t = s_a / (s_a - s_b)
                r = (1 - t) * along[e0] + t * along[e1]
                key = (min(e0, e1), max(e0, e1))
                if key not in points:
                    points[key] = (r, (1 - t) * pts[e0, ia] + t * pts[e1, ia])
                hits.append(key)
            if len(hits) != 2:
                continue
            (ra, za), (rb, zb) = points[hits[0]], points[hits[1]]
            if ra <= 0 and rb <= 0:
                continue                        # the opposite half-plane
            a, b = hits
            if ra <= 0 or rb <= 0:
                # The section crosses the axis inside this triangle: it ends
                # on the axis at the crossing point.
                if ra <= 0:
                    (ra, za), (rb, zb), a, b = (rb, zb), (ra, za), b, a
                f = ra / (ra - rb)
                b = ("axis", ti)
                points[b] = (0.0, za + f * (zb - za))
            adjacency.setdefault(a, []).append(b)
            adjacency.setdefault(b, []).append(a)
        if not adjacency:
            raise ValueError("the source surface does not cross the meridian "
                             f"half-plane about the {self.axis} axis")
        bad = [k for k, v in adjacency.items() if len(v) > 2]
        if bad:
            raise ValueError(
                "the meridian section branches (non-manifold source surface "
                f"at {len(bad)} points); select a single closed or open "
                "surface of revolution")
        chains = []
        seen = set()
        starts = [k for k, v in adjacency.items() if len(v) == 1]
        for start in starts + list(adjacency):
            if start in seen:
                continue
            order = [start]
            seen.add(start)
            prev, cur = None, start
            closed = False
            while True:
                nxt = [n for n in adjacency[cur] if n != prev]
                if not nxt:
                    break
                n = nxt[0]
                if n == start:
                    closed = True
                    break
                if n in seen:
                    break
                order.append(n)
                seen.add(n)
                prev, cur = cur, n
            rz = np.asarray([points[k] for k in order], float)
            if closed:
                rz = np.vstack([rz, rz[:1]])
            if len(rz) < 2:
                continue
            s = np.r_[0.0, np.cumsum(np.hypot(*np.diff(rz, axis=0).T))]
            chains.append(_Chain(rz=rz, s=s, closed=closed))
        return chains

    # -- evaluation ---------------------------------------------------------
    def locate(self, rz: np.ndarray):
        d, idx = self._tree.query(np.asarray(rz, float))
        return self._dense_chain[idx], self._dense_s[idx], d

    def evaluate(self, rz: np.ndarray, *, max_distance: float | None = None,
                 what: str = "target points"):
        """Profile value at meridian points ``rz``; returns (values, dist)."""
        rz = np.asarray(rz, float)
        tol = self.axisymmetry_tolerance if max_distance is None else max_distance
        chain, s, d = self.locate(rz)
        bad = d > tol
        if np.any(bad):
            worst = int(np.argmax(d))
            raise ValueError(
                f"{int(bad.sum())}/{len(d)} {what} lie farther than "
                f"{tol:.3e} m from the source meridian (worst {d[worst]:.3e} m "
                f"at (r, a) = {np.round(rz[worst], 6).tolist()}). The target "
                "surface is not the same body of revolution as the source.")
        out = np.empty(len(rz))
        for k in np.unique(chain):
            m = chain == k
            ch = self.chains[k]
            if ch.closed:
                period = ch.s[-1]
                out[m] = np.interp(s[m], ch.centres, ch.q, period=period)
            else:
                out[m] = np.interp(s[m], ch.centres, ch.q)
        return out, d

    def evaluate_xyz(self, xyz: np.ndarray, **kwargs):
        return self.evaluate(to_meridian(np.asarray(xyz, float), self.axis),
                             **kwargs)

    def audit(self) -> dict:
        return {
            "method": "exact-meridian-bin-integration",
            "axis": self.axis,
            "bin_width_m": self.bin_width,
            "chain_bin_widths_m": [float(w) for w in self._widths],
            "source_h_median_m": self.h_median,
            "n_chains": len(self.chains),
            "chain_lengths_m": [float(ch.s[-1]) for ch in self.chains],
            "source_power_W": self.source_power,
            "profile_power_W": self.power,
            "max_source_to_meridian_m": self.max_sample_distance,
            "axisymmetry_tolerance_m": self.axisymmetry_tolerance,
        }


# ---------------------------------------------------------------------------
# Power bookkeeping
# ---------------------------------------------------------------------------

def power_balance(target_power: float, reference_power: float,
                  tolerance: float, *, what: str,
                  scale: float | None = None, hint: str = "") -> dict:
    """Compare two powers; raise when they differ by more than ``tolerance``.

    The difference is normalised by ``scale`` (default ``|reference|``).
    Pass the integral of ``|q|`` for signed test fields whose net power is
    near zero.
    """
    if not math.isfinite(target_power) or not math.isfinite(reference_power):
        raise ValueError(f"{what}: non-finite power ({target_power!r}, "
                         f"{reference_power!r})")
    norm = abs(reference_power) if scale is None else max(abs(scale),
                                                          abs(reference_power))
    if norm == 0.0:
        rel = 0.0 if target_power == 0.0 else math.inf
    else:
        rel = (target_power - reference_power) / norm
    record = {"target_W": float(target_power),
              "reference_W": float(reference_power),
              "normalisation_W": float(norm),
              "relative_error": float(rel), "tolerance": float(tolerance)}
    if abs(rel) > tolerance:
        raise ValueError(
            f"{what}: {target_power:.6e} W against {reference_power:.6e} W "
            f"(relative error {rel:+.3%}, tolerance {tolerance:.3%}). The heat "
            "source was not transferred faithfully; refusing to solve."
            + (f" {hint}" if hint else ""))
    return record


# ---------------------------------------------------------------------------
# Electromagnetic heat source -> thermal mesh
# ---------------------------------------------------------------------------

QSURF_QUANTITY = "surface_loss_density"
TEMPERATURE_QUANTITY = "temperature"
TEMPERATURE_UNIT = "degC"
QSURF_UNIT = "W/m^2"
DEFAULT_POWER_TOLERANCE = 0.02
TRANSFER_POWER_HINT = (
    "If the thermal heat-flux boundaries cover only part of the EM heated "
    "surface, select the matching EM boundaries with --em-heat-boundaries.")
SIDECAR_POWER_TOLERANCE = 1.0e-3


def check_tolerance(value, name, upper=0.5):
    """A relative tolerance must satisfy 0 < value <= upper."""
    v = float(value)
    if not (math.isfinite(v) and 0.0 < v <= upper):
        raise ValueError(f"{name} must lie in (0, {upper}], got {value!r}")
    return v


def boundaries_with_nonzero_field(mesh, gf) -> list[str]:
    """Boundary names on which a P1 field is not identically zero."""
    from ngsolve import BND

    vals = np.asarray(gf.vec.FV().NumPy(), float)[:mesh.nv]
    out = set()
    for el in mesh.Elements(BND):
        if el.mat in out:
            continue
        if any(vals[v.nr] != 0.0 for v in el.vertices):
            out.add(el.mat)
    return sorted(out)


def resolve_em_heat_boundaries(em_mesh, gf_em, *, requested=None,
                               sidecar=None) -> tuple:
    """The EM boundaries that carry the heat source.

    They are the explicit ``requested`` selector, or else the list the EM
    solver recorded in the sidecar.  Nothing is borrowed from the thermal
    mesh's names.  Returns ``(names, how)``.
    """
    import re

    available = sorted(set(em_mesh.GetBoundaries()))
    if requested:
        pattern = re.compile(str(requested))
        names = [n for n in available if pattern.fullmatch(n)]
        if not names:
            raise ValueError(f"--em-heat-boundaries={requested!r} matches no "
                             f"EM boundary; available: {available}")
        return names, "explicit"
    if sidecar and sidecar.get("boundaries"):
        names = list(sidecar["boundaries"])
        missing = [n for n in names if n not in available]
        if missing:
            raise ValueError(f"sidecar boundaries {missing} are not on the EM "
                             f"mesh; available: {available}")
        return names, "sidecar"
    raise ValueError(
        "the q_surf sidecar records no heated EM boundaries and "
        "--em-heat-boundaries was not given. The field is non-zero on "
        f"{boundaries_with_nonzero_field(em_mesh, gf_em)}; state the heated "
        "boundaries explicitly.")


class EMHeatSource:
    """A verified EM surface heat source ready to be transferred.

    ``mode`` is ``"direct"`` (pointwise, optionally at a rotated body
    angle) or ``"phi-average"`` (exact circumferential average).
    """

    def __init__(self, qsurf_sol: str, em_vol: str, *,
                 em_boundaries=None, qsurf_order: int = 1,
                 mode: str = "direct", axis: str = "z",
                 q_scale: float = 1.0, transfer_tolerance=None,
                 bin_width=None):
        from ngsolve import GridFunction, H1, Mesh

        if int(qsurf_order) != 1:
            raise ValueError("the q_surf handoff is H1 order 1 only")
        if mode not in ("direct", "phi-average"):
            raise ValueError(f"unknown transfer mode {mode!r}")
        if not (math.isfinite(q_scale) and q_scale > 0):
            raise ValueError(f"q_scale must be positive, got {q_scale!r}")
        self.qsurf_sol = os.path.abspath(qsurf_sol)
        self.em_vol = os.path.abspath(em_vol)
        for p, what in ((self.qsurf_sol, "--qsurf-sol"),
                        (self.em_vol, "--em-vol")):
            if not os.path.isfile(p):
                raise FileNotFoundError(f"{what} not found: {p}")
        self.mode = mode
        self.axis = axis
        self.q_scale = float(q_scale)
        self.em_mesh = Mesh(self.em_vol)
        self.pair_audit = verify_field_pair(
            self.qsurf_sol, self.em_vol, self.em_mesh, 1,
            quantity=QSURF_QUANTITY)
        self.sidecar = read_field_sidecar(self.qsurf_sol)
        self.gf = GridFunction(H1(self.em_mesh, order=1))
        self.gf.Load(self.qsurf_sol)
        if self.q_scale != 1.0:
            self.gf.vec.data = self.q_scale * self.gf.vec
        self.em_boundaries, self.boundary_rule = resolve_em_heat_boundaries(
            self.em_mesh, self.gf, requested=em_boundaries,
            sidecar=self.sidecar)
        self.surface = SurfaceP1Field.from_gridfunction(
            self.em_mesh, self.gf, self.em_boundaries)
        # Reference power on the mesh geometry as loaded (curved elements
        # included); the flat-triangle value only feeds the transfer.
        from ngsolve import BND, Integrate
        self.power = float(Integrate(
            self.gf, self.em_mesh, BND,
            definedon=self.em_mesh.Boundaries("|".join(self.em_boundaries))
        ).real)
        self.power_flat = self.surface.power
        if self.sidecar.get("P_wp_W") is None:
            raise ValueError(
                f"{os.path.basename(self.qsurf_sol)}: the sidecar records no "
                "EM power P_wp_W, so the heat source cannot be checked; "
                "write the sidecar with the EM run's power")
        if self.boundary_rule == "sidecar":
            self.sidecar_power = power_balance(
                self.power, float(self.sidecar["P_wp_W"]) * self.q_scale,
                SIDECAR_POWER_TOLERANCE,
                what="q_surf integral against the EM power in its sidecar")
        else:
            # explicit boundaries may be a subset of the recorded heated
            # surface; the thermal power gate still compares like with like
            self.sidecar_power = {"skipped": "explicit --em-heat-boundaries",
                                  "sidecar_P_wp_W": float(
                                      self.sidecar["P_wp_W"])}
        h = self.surface.h_median
        if transfer_tolerance is None:
            self.transfer_tolerance = 0.25 * h
        else:
            self.transfer_tolerance = float(transfer_tolerance)
            if not (0.0 < self.transfer_tolerance <= h):
                raise ValueError(
                    f"transfer tolerance must lie in (0, {h:.3e}] m (the "
                    "median EM surface edge); a larger value would accept "
                    f"points off the surface, got {transfer_tolerance!r}")
        self.profile = None
        if mode == "phi-average":
            self.profile = AxisymmetricSurfaceProfile(
                self.surface.points, self.surface.tris, self.surface.values,
                axis=axis, bin_width=bin_width,
                axisymmetry_tolerance=self.transfer_tolerance)

    def transfer(self, target_mesh, heat_names, gf_target, *,
                 power_tolerance, theta: float = 0.0, weight=None,
                 meridian: bool = False) -> dict:
        """Write the source onto the thermal heat-flux vertices and gate it.

        ``gf_target`` is an H1 order-1 GridFunction on ``target_mesh``.  For
        an axisymmetric (r, z) thermal mesh pass ``meridian=True`` and the
        revolution weight ``2 pi r``.  The transferred power must match the
        EM power within ``power_tolerance`` or ``ValueError`` is raised.
        """
        from ngsolve import BND, CF, Integrate

        tol = check_tolerance(power_tolerance, "power tolerance")
        vnrs = boundary_vertex_numbers(target_mesh, heat_names)
        pts = mesh_vertices(target_mesh)[vnrs]
        if meridian:
            values, dist = self.values_at_meridian(pts[:, :2])
        else:
            values, dist = self.values_at_xyz(pts, theta)
        vec = gf_target.vec.FV().NumPy()
        vec[:] = 0.0
        vec[vnrs] = values
        w = CF(1.0) if weight is None else weight
        power = float(Integrate(gf_target * w, target_mesh, BND,
                                definedon=target_mesh.Boundaries(
                                    "|".join(heat_names))).real)
        balance = power_balance(
            power, self.power, tol,
            what="thermal heat input against the EM source power",
            scale=self.surface.magnitude_power, hint=TRANSFER_POWER_HINT)
        return {"theta_rad": float(theta), "target_vertices": int(len(vnrs)),
                "max_transfer_distance_m": float(dist.max()),
                "power_balance": balance, "vertex_numbers": vnrs}

    def values_at_xyz(self, xyz: np.ndarray, theta: float = 0.0):
        """Source values at thermal points ``xyz`` (body angle ``theta``).

        Ungated: use :meth:`transfer`, which also checks the power."""
        xyz = np.asarray(xyz, float)
        if self.profile is not None:
            return self.profile.evaluate_xyz(xyz, what="thermal heat-flux "
                                             "vertices")
        world = rotate_about_axis(xyz, theta, self.axis) if theta else xyz
        return self.surface.evaluate(world,
                                     max_distance=self.transfer_tolerance,
                                     what="thermal heat-flux vertices")

    def values_at_meridian(self, rz: np.ndarray):
        """Circumferential average at meridian points (2D thermal mesh)."""
        if self.profile is None:
            raise ValueError("values_at_meridian needs mode='phi-average'")
        return self.profile.evaluate(np.asarray(rz, float),
                                     what="axisymmetric heat-flux vertices")

    def audit(self) -> dict:
        out = {
            "mode": self.mode,
            "qsurf_sol": self.qsurf_sol.replace("\\", "/"),
            "em_vol": self.em_vol.replace("\\", "/"),
            "provenance": self.pair_audit["provenance"],
            "em_heat_boundaries": list(self.em_boundaries),
            "em_heat_boundary_rule": self.boundary_rule,
            "q_scale": self.q_scale,
            "source_power_W": self.power,
            "source_power_flat_W": self.power_flat,
            "sidecar_power_check": self.sidecar_power,
            "transfer_tolerance_m": self.transfer_tolerance,
            "source_h_median_m": self.surface.h_median,
        }
        if self.profile is not None:
            out["phi_average"] = self.profile.audit()
        return out


# ---------------------------------------------------------------------------
# Rotating workpiece: angle-resolved source and its exact step average
# ---------------------------------------------------------------------------

ROTOR_STATES_SCHEMA = "radia.ih-rotor-states/1"
ROTOR_FRAMES = ("world", "body")
DEFAULT_ANGLE_STEP_TOLERANCE = 0.25


def load_rotor_states(path: str, *, em_boundaries=None, q_scale=1.0,
                      transfer_tolerance=None):
    """EM solutions of a rotating workpiece, one per rotor angle.

    The manifest (schema ``radia.ih-rotor-states/1``) lists ``states`` of
    ``{"angle_rad", "qsurf_sol", "em_vol"}`` (paths relative to the manifest),
    the rotation ``axis``, the ``period_rad`` after which the arrangement
    repeats (2 pi, or 2 pi / n for an n-fold part) and the ``frame``:

    ``world``  each EM mesh holds the part turned by ``angle_rad``;
    ``body``   each EM mesh holds the part unturned and the coil turned by
               ``-angle_rad``.

    Returns ``(states, meta)`` with ``states`` a list of
    ``(angle_rad, EMHeatSource)`` sorted by angle.
    """
    with open(path, encoding="utf-8") as fh:
        man = json.load(fh)
    if man.get("schema") != ROTOR_STATES_SCHEMA:
        raise ValueError(f"{path}: schema must be {ROTOR_STATES_SCHEMA!r}, "
                         f"got {man.get('schema')!r}")
    for key in ("axis", "period_rad", "frame", "states"):
        if key not in man:
            raise ValueError(f"{path}: missing {key!r}")
    frame = man["frame"]
    if frame not in ROTOR_FRAMES:
        raise ValueError(f"{path}: frame must be one of {ROTOR_FRAMES}")
    period = float(man["period_rad"])
    n_fold = 2 * math.pi / period if period > 0 else 0.0
    if not (period > 0 and abs(n_fold - round(n_fold)) < 1e-9):
        raise ValueError(f"{path}: period_rad must be 2 pi / n, got {period!r}")
    base = os.path.dirname(os.path.abspath(path))
    states = []
    for i, st in enumerate(man["states"]):
        for key in ("angle_rad", "qsurf_sol", "em_vol"):
            if key not in st:
                raise ValueError(f"{path}: state {i} lacks {key!r}")
        ang = float(st["angle_rad"])
        if not (0.0 <= ang < period):
            raise ValueError(f"{path}: state {i} angle {ang!r} is outside "
                             f"[0, period_rad = {period!r})")
        src = EMHeatSource(os.path.join(base, st["qsurf_sol"]),
                           os.path.join(base, st["em_vol"]),
                           em_boundaries=em_boundaries, q_scale=q_scale,
                           transfer_tolerance=transfer_tolerance)
        recorded = src.sidecar.get("rotor_angle_rad")
        if recorded is not None and abs(float(recorded) - ang) > 1e-9:
            raise ValueError(f"{st['qsurf_sol']}: its sidecar records rotor "
                             f"angle {recorded!r}, the manifest {ang!r}")
        states.append((ang, src))
    if len(states) < 2:
        raise ValueError(f"{path}: a rotor-states manifest needs at least two "
                         "states; a body of revolution needs none (one EM "
                         "solution is turned with the part)")
    states.sort(key=lambda s: s[0])
    angles = [a for a, _ in states]
    if np.any(np.diff(angles) <= 0):
        raise ValueError(f"{path}: rotor angles must be distinct")
    return states, {"manifest": os.path.abspath(path).replace("\\", "/"),
                    "axis": str(man["axis"]).lower(), "period_rad": period,
                    "frame": frame, "n_states": len(states)}


class RotatingSurfaceSource:
    """Heat flux on a rotating part, resolved in its angle.

    The source at body angle ``theta`` is tabulated at knots over one
    revolution and taken piecewise linear in between.  Knots come either
    from one EM solution of a body of revolution, turned with the part at a
    spacing that resolves the EM surface mesh, or from one EM solution per
    rotor angle over the symmetry ``period`` (:func:`load_rotor_states`).
    An n-fold part repeats its EM solution every ``period`` = 2 pi / n, but a
    point of the part returns only after a revolution: the knot at
    ``a_k + m period`` reads state ``k`` at the point turned by
    ``a_k + m period`` (world frame) or ``m period`` (body frame), which lies
    on the part only if it really is n-fold -- otherwise the transfer fails.
    Every knot is transferred with the power gate.

    :meth:`average` writes the exact time average over the angles a time
    step sweeps -- a step longer than a revolution averages whole
    revolutions -- so a coarse time step does not sample the rotation
    stroboscopically.  :meth:`at` writes the value at one angle.
    """

    def __init__(self, target_mesh, heat_names, gf_target, *, states,
                 power_tolerance, axis="z", period=2 * math.pi,
                 frame="world", angle_step_tolerance=None):
        import time

        t0 = time.perf_counter()
        self.axis = str(axis).lower()
        _axis_frame(self.axis)
        if frame not in ROTOR_FRAMES:
            raise ValueError(f"frame must be one of {ROTOR_FRAMES}")
        self.state_period = float(period)
        n_fold = 2 * math.pi / self.state_period
        if abs(n_fold - round(n_fold)) > 1e-9:
            raise ValueError("period must be 2 pi / n")
        n_fold = int(round(n_fold))
        self.period = 2 * math.pi          # of the body-frame source
        self.gf = gf_target
        self.tol_step = check_tolerance(
            DEFAULT_ANGLE_STEP_TOLERANCE if angle_step_tolerance is None
            else angle_step_tolerance, "angle step tolerance", upper=1.0)
        self.vnrs = boundary_vertex_numbers(target_mesh, heat_names)
        pts = mesh_vertices(target_mesh)[self.vnrs]
        r_max = float(np.max(to_meridian(pts, self.axis)[:, 0]))
        if len(states) == 1:
            _, src = states[0]
            if src.mode != "direct":
                raise ValueError("a rotating source is transferred pointwise; "
                                 "the phi-average is already rotation-free")
            # resolve the EM surface mesh along the path of the outermost
            # point: half an EM edge per knot
            if n_fold != 1:
                raise ValueError("one turned EM solution needs period 2 pi")
            n = max(8, int(math.ceil(self.period * r_max
                                     / (0.5 * src.surface.h_median))))
            knots = [(self.period * k / n, src, self.period * k / n)
                     for k in range(n)]
            self.kind = "turned-single-state"
        else:
            knots = [(a + m * self.state_period, src,
                      a + m * self.state_period if frame == "world"
                      else m * self.state_period)
                     for m in range(n_fold) for a, src in states]
            self.kind = "rotor-states"
        self.n_fold = n_fold
        self.angles = np.asarray([a for a, _, _ in knots], float)
        vals, records = [], []
        for ang, src, theta in knots:
            try:
                rec = src.transfer(target_mesh, heat_names, gf_target,
                                   power_tolerance=power_tolerance,
                                   theta=theta)
            except ValueError as exc:
                if len(states) == 1:
                    raise ValueError(
                        f"at body angle {ang:.4f} rad: {exc}  The turned part "
                        "leaves the EM surface, so it is not a body of "
                        "revolution about the rotation axis; give one EM "
                        "solution per rotor angle (--rotor-states).") from exc
                raise
            if not np.array_equal(rec.pop("vertex_numbers"), self.vnrs):
                raise RuntimeError("heat-flux vertex numbering changed")
            vals.append(gf_target.vec.FV().NumPy()[self.vnrs].copy())
            records.append(rec)
        self.values = np.asarray(vals)
        self.knot_records = records
        # angular resolution: change between neighbouring knots, relative
        # to the largest knot field
        ref = float(np.max(np.linalg.norm(self.values, axis=1)))
        nxt = np.roll(self.values, -1, axis=0)
        change = np.linalg.norm(nxt - self.values, axis=1) / max(ref, 1e-300)
        self.max_step_change = float(change.max())
        worst = int(np.argmax(change))
        if self.max_step_change > self.tol_step:
            raise ValueError(
                f"the rotor angles are too coarse: the source changes by "
                f"{self.max_step_change:.1%} between {self.angles[worst]:.4f} "
                f"and the next knot (limit {self.tol_step:.0%}); add EM "
                "solutions between them")
        # segment integrals for the exact average
        a_ext = np.r_[self.angles, self.angles[0] + self.period]
        self._L = np.diff(a_ext)
        v_ext = np.vstack([self.values, self.values[:1]])
        self._seg = 0.5 * self._L[:, None] * (v_ext[:-1] + v_ext[1:])
        self._cum = np.vstack([np.zeros(len(self.vnrs)),
                               np.cumsum(self._seg, axis=0)])
        self._v_ext = v_ext
        self.build_s = time.perf_counter() - t0
        self.at(self.angles[0])

    # -- evaluation ---------------------------------------------------------
    def _locate(self, theta):
        u = (float(theta) - self.angles[0]) % self.period
        j = int(np.searchsorted(self.angles - self.angles[0], u,
                                side="right") - 1)
        j = min(max(j, 0), len(self.angles) - 1)
        s = u - (self.angles[j] - self.angles[0])
        return j, s, u

    def _value(self, theta):
        j, s, _ = self._locate(theta)
        v0, v1 = self._v_ext[j], self._v_ext[j + 1]
        return v0 + (v1 - v0) * (s / self._L[j])

    def _integral(self, theta):
        """Integral of the source from the first knot to ``theta``."""
        m = math.floor((float(theta) - self.angles[0]) / self.period)
        j, s, _ = self._locate(theta)
        v0, v1 = self._v_ext[j], self._v_ext[j + 1]
        part = s * v0 + (v1 - v0) * s * s / (2.0 * self._L[j])
        return m * self._cum[-1] + self._cum[j] + part

    def _write(self, values):
        vec = self.gf.vec.FV().NumPy()
        vec[:] = 0.0
        vec[self.vnrs] = values

    def at(self, theta):
        self._write(self._value(theta))

    def average(self, theta0, theta1):
        """Write the time average over body angles ``theta0 .. theta1``."""
        if theta1 < theta0:
            raise ValueError("theta1 must not precede theta0")
        if theta1 - theta0 <= 1e-15 * self.period:
            self.at(theta0)
            return
        self._write((self._integral(theta1) - self._integral(theta0))
                    / (theta1 - theta0))

    def audit(self) -> dict:
        return {"kind": self.kind, "axis": self.axis,
                "state_period_rad": self.state_period, "n_fold": self.n_fold,
                "knots": int(len(self.angles)),
                "knot_spacing_max_rad": float(self._L.max()),
                "max_neighbour_change": self.max_step_change,
                "angle_step_tolerance": self.tol_step,
                "knot_power_relative_error_max": float(max(
                    abs(r["power_balance"]["relative_error"])
                    for r in self.knot_records)),
                "build_s": round(self.build_s, 3)}


def lumped_surface_p1(mesh, q_cf, region, *, bonus_intorder=4):
    """P1 handoff of a surface loss density by the lumped L2 projection.

    ``q_i = int(q phi_i) / int(phi_i)`` over the ``region`` boundary
    (a ``mesh.Boundaries`` region).  It conserves the integral of ``q``
    exactly and keeps ``q >= 0``, also where ``q`` is singular at convex
    edges, where an interpolating ``Set`` overshoots.  Off the region the
    field is zero.  Returns ``(gf, on_region_mask)``.
    """
    from ngsolve import GridFunction, H1, LinearForm, ds
    fes = H1(mesh, order=1)
    gf = GridFunction(fes)
    _, v = fes.TnT()
    num_lf = LinearForm(fes)
    num_lf += q_cf * v * ds(definedon=region, bonus_intorder=bonus_intorder)
    num_lf.Assemble()
    den_lf = LinearForm(fes)
    den_lf += v * ds(definedon=region)
    den_lf.Assemble()
    num = np.asarray(num_lf.vec.FV().NumPy()).real
    den = np.asarray(den_lf.vec.FV().NumPy()).real
    on = den > 0
    vals = np.zeros_like(den)
    vals[on] = num[on] / den[on]
    gf.vec.FV().NumPy()[:] = vals
    return gf, on


def scale_field_artifact(sol_path: str, out_path: str, factor: float) -> str:
    """Write ``factor * field`` with a sidecar recording the scaling.

    Used for current sweeps (q scales with I^2) so a rescaled source stays a
    verifiable artifact instead of an anonymous coefficient vector.
    """
    record = read_field_sidecar(sol_path)
    if record is None:
        raise ValueError(f"{os.path.basename(sol_path)} has no sidecar; "
                         "write one before scaling")
    data = np.fromfile(sol_path, dtype="<f8")
    if data.size != int(record["ndof"]):
        raise ValueError("field size disagrees with its sidecar")
    (float(factor) * data).astype("<f8").tofile(out_path)
    extra = {k: v for k, v in record.items() if k not in {
        "schema", "quantity", "unit", "space", "fes_order", "ndof",
        "sol_file", "sol_sha256", "mesh_file", "mesh_sha256", "mesh_nv",
        "mesh_ne", "mesh_dim", "mesh_curve_order", "boundaries",
        "definedon"}}
    if extra.get("P_wp_W") is not None:
        extra["P_wp_W"] = float(extra["P_wp_W"]) * float(factor)
    extra["scaled_from"] = {"sol_file": record["sol_file"],
                            "sol_sha256": record["sol_sha256"],
                            "factor": float(factor)}
    from ngsolve import Mesh
    mesh = Mesh(record["mesh_file"])
    return write_field_sidecar(
        out_path, mesh_path=record["mesh_file"], mesh=mesh,
        fes_order=record["fes_order"], quantity=record["quantity"],
        unit=record["unit"], boundaries=record["boundaries"], extra=extra,
        definedon=record.get("definedon"))

# ---------------------------------------------------------------------------
# Temperature-dependent surface source (local surface-impedance model)
# ---------------------------------------------------------------------------

HT_QUANTITY = "surface_tangential_H_amplitude"
HT_UNIT = "A/m"
FROZEN_HT_REFUSAL = (
    "the --em-table source freezes the EM run's |H_t|; for a ferromagnetic "
    "workpiece it overstates heating below the Curie point and predicts a "
    "false self-limit above it (validation_test/induction_heating/results/"
    "coupled_curie_cylinder_frozen_ht.json). Pass --allow-frozen-ht to use "
    "it knowingly (non-magnetic parts, or between EM re-solves), or use the "
    "axisymmetric coupled EM-thermal solver.")


def load_impedance_table(path: str):
    """Load an (|H_t|, T) surface-impedance table written by calc_em_table.

    Returns a dict with ``H`` (A/m, peak), ``T`` (degC), ``q`` (W/m^2 at
    that peak |H_t|, ``q = Re(Z_s) |H_t|^2 / 2``) and ``meta``.
    """
    data = np.load(path, allow_pickle=False)
    for key in ("H_grid", "T_grid", "Zs_re", "q_surf", "meta"):
        if key not in data.files:
            raise ValueError(f"{path} is not an impedance table (no {key})")
    H = np.asarray(data["H_grid"], float)
    T = np.asarray(data["T_grid"], float)
    q = np.asarray(data["q_surf"], float)
    if not (np.all(np.diff(H) > 0) and np.all(np.diff(T) > 0)):
        raise ValueError(f"{path}: table axes must increase")
    if q.shape != (len(H), len(T)) or np.any(q <= 0):
        raise ValueError(f"{path}: q table must be positive, shape (nH, nT)")
    if np.any(np.diff(q, axis=0) <= 0):
        raise ValueError(f"{path}: q must increase with |H_t| at every T")
    return {"path": os.path.abspath(path), "H": H, "T": T, "q": q,
            "meta": json.loads(str(data["meta"]))}


class ImpedanceTable:
    """Bilinear lookup of q(|H_t|, T) in (log H, T); strict range checks."""

    def __init__(self, table: dict, *, allow_extrapolation: bool = False):
        self.table = table
        self.logH = np.log(table["H"])
        self.T = table["T"]
        self.logq = np.log(table["q"])
        self.allow = bool(allow_extrapolation)
        self.max_T_excursion = 0.0
        self.max_H_excursion = 0.0

    def _check_T(self, T):
        lo, hi = self.T[0], self.T[-1]
        over = max(0.0, float(np.max(T) - hi), float(lo - np.min(T)))
        if over > 0:
            if not self.allow:
                raise ValueError(
                    f"surface temperature {np.min(T):.1f}..{np.max(T):.1f} C "
                    f"leaves the impedance table {lo:.1f}..{hi:.1f} C; "
                    "extend the table or allow extrapolation explicitly")
            self.max_T_excursion = max(self.max_T_excursion, over)
        return np.clip(T, lo, hi)

    def q(self, H, T):
        """q [W/m^2] at peak |H_t| ``H`` and temperature ``T`` (broadcast).

        ``H = 0`` gives ``q = 0``.  Any other |H_t| outside the table's H
        grid is an error unless extrapolation is allowed, in which case q
        continues with the local ``q ~ |H_t|^2`` trend and the excursion is
        recorded.
        """
        H = np.asarray(H, float)
        if np.any(H < 0) or not np.all(np.isfinite(H)):
            raise ValueError("|H_t| must be finite and non-negative")
        T = self._check_T(np.asarray(T, float))
        Hmin, Hmax = np.exp(self.logH[0]), np.exp(self.logH[-1])
        live = H > 0
        Hs = np.where(live, H, Hmin)
        lh = np.log(Hs)
        out_lo = live & (lh < self.logH[0] - 1e-12)
        out_hi = lh > self.logH[-1] + 1e-12
        if np.any(out_lo) or np.any(out_hi):
            if not self.allow:
                raise ValueError(
                    f"|H_t| {float(H[live].min()):.3e}..{float(H.max()):.3e} "
                    f"A/m leaves the impedance table {Hmin:.3e}..{Hmax:.3e} "
                    "A/m; extend the table or allow extrapolation explicitly")
            if np.any(out_hi):
                self.max_H_excursion = max(self.max_H_excursion,
                                           float(H.max()) / Hmax)
            if np.any(out_lo):
                self.max_H_excursion = max(self.max_H_excursion,
                                           Hmin / float(H[out_lo].min()))
        lhc = np.clip(lh, self.logH[0], self.logH[-1])
        i = np.clip(np.searchsorted(self.logH, lhc) - 1, 0, len(self.logH) - 2)
        j = np.clip(np.searchsorted(self.T, T) - 1, 0, len(self.T) - 2)
        u = (lhc - self.logH[i]) / (self.logH[i + 1] - self.logH[i])
        w = (T - self.T[j]) / (self.T[j + 1] - self.T[j])
        L = self.logq
        lq = ((1 - u) * (1 - w) * L[i, j] + u * (1 - w) * L[i + 1, j]
              + (1 - u) * w * L[i, j + 1] + u * w * L[i + 1, j + 1])
        q = np.exp(lq)
        # (reached only with extrapolation allowed) q ~ |H_t|^2 outside
        q = np.where(out_lo, q * np.exp(2.0 * (lh - self.logH[0])), q)
        q = np.where(out_hi, q * np.exp(2.0 * (lh - self.logH[-1])), q)
        return q * live


class TemperatureDependentSource:
    """q_i(T) = q_ref,i * <q_tab(s H_ik, T_i)>_k / <q_tab(H_ik, T_ref)>_k.

    ``q_ref`` is the transferred EM source at the reference temperature
    (exact circumferential average or pointwise), ``H_ik`` the |H_t| samples
    that feed vertex ``i`` (one, or one per azimuth for a rotating part),
    ``s`` the current scale.  At ``T = T_ref`` and ``s = 1`` the source
    equals the EM solution exactly; the table only supplies the local
    temperature and current dependence.
    """

    def __init__(self, table: ImpedanceTable, q_ref, H_samples, vertex_dofs,
                 gf_target, *, T_ref: float, acknowledge_frozen_ht: bool,
                 H_scale: float = 1.0, gf_damping=None, dT_fd: float = 1.0):
        if acknowledge_frozen_ht is not True:
            raise ValueError(FROZEN_HT_REFUSAL)
        if not (math.isfinite(H_scale) and H_scale > 0):
            raise ValueError(f"H_scale must be positive, got {H_scale!r}")
        self.table = table
        self.q_ref = np.asarray(q_ref, float)
        self.H = np.atleast_2d(np.asarray(H_samples, float))
        if self.H.shape[0] != len(self.q_ref):
            self.H = self.H.T
        self.dofs = np.asarray(vertex_dofs, int)
        self.gf = gf_target
        self.T_ref = float(T_ref)
        self.s = float(H_scale)
        denom = table.q(self.H, np.full(self.H.shape, self.T_ref)).mean(axis=1)
        dead = (denom <= 0) & (self.q_ref != 0)
        if np.any(dead):
            raise ValueError(
                f"{int(dead.sum())} heat-flux vertices carry q_surf but zero "
                "|H_t|; the |H_t| field does not belong to this source")
        self.ratio_ref = np.where(denom > 0, 1.0 / np.maximum(denom, 1e-300),
                                  0.0)
        self.T_surface_range = [np.inf, -np.inf]
        # max(-dq/dT, 0) on the same vertices: the stabilising part of the
        # source derivative, added to the Newton matrix by the stepper.
        self.gf_damping = gf_damping
        self.dT_fd = float(dT_fd)

    def values(self, T_vertices, *, record_temperature=True):
        T = np.asarray(T_vertices, float)
        if record_temperature:
            self.T_surface_range = [min(self.T_surface_range[0], float(T.min())),
                                    max(self.T_surface_range[1], float(T.max()))]
        num = self.table.q(self.s * self.H,
                           np.broadcast_to(T[:, None], self.H.shape))
        return self.q_ref * num.mean(axis=1) * self.ratio_ref

    def update(self, gfT, temperature_dofs):
        T = np.asarray(gfT.vec.FV().NumPy())[temperature_dofs]
        vec = self.gf.vec.FV().NumPy()
        q = self.values(T)
        vec[self.dofs] = q
        if self.gf_damping is not None:
            lo, hi = self.table.T[0], self.table.T[-1]
            Tp = np.minimum(T + self.dT_fd, hi)
            Tm = np.maximum(T - self.dT_fd, lo)
            span = np.maximum(Tp - Tm, 1e-12)
            dq = (self.values(Tp, record_temperature=False)
                  - self.values(Tm, record_temperature=False)) / span
            self.gf_damping.vec.FV().NumPy()[self.dofs] = np.maximum(-dq, 0.0)

    def audit(self) -> dict:
        return {"model": "ratio-of-local-surface-impedance",
                "table": self.table.table["path"].replace("\\", "/"),
                "table_meta": self.table.table["meta"],
                "T_ref_C": self.T_ref, "H_scale": self.s,
                "azimuth_samples": int(self.H.shape[1]),
                "surface_T_range_C": self.T_surface_range,
                "table_T_excursion_C": self.table.max_T_excursion,
                "table_H_excursion_ratio": self.table.max_H_excursion}

def build_temperature_dependent_source(source: "EMHeatSource", *, table_path,
                                       target_xyz, q_ref, gf_target, dofs,
                                       T_ref, H_scale=1.0, ht_sol="",
                                       azimuths=0, allow_extrapolation=False,
                                       consistency_tolerance=0.05,
                                       acknowledge_frozen_ht=False):
    """Couple a transferred EM source to the local surface temperature.

    ``source`` is the verified reference source, ``target_xyz`` the thermal
    heat-flux points (for an axisymmetric thermal mesh, ``(r, 0, z)``),
    ``q_ref`` the transferred reference values there.  ``|H_t|`` comes from
    ``ht_sol``, the ``_Ht.sol`` artifact the EM run wrote next to its
    ``q_surf``.  With
    ``azimuths > 0`` each point takes |H_t| on that many azimuths of its
    ring (a rotating part); otherwise at the point itself.

    The table is checked against the EM run: its power at ``T_ref`` must
    match the EM power within ``consistency_tolerance``, otherwise the table
    does not describe the material that was solved.

    **Validity.**  The model freezes the spatial |H_t| of the EM run.  For
    a current-driven coil around a ferromagnetic workpiece that is wrong:
    |H_t| falls as sigma(T) falls and rises when the surface passes the
    Curie band.  In the validation case the true power stays within +5 %
    until the surface reaches the Curie band and peaks at +34 % as the band
    passes, while this model rises to +53 % below the Curie point and then
    stalls the surface at the Curie temperature, a false self-limit
    (validation_test/induction_heating/results/
    coupled_curie_cylinder_frozen_ht.json).  It is therefore refused unless
    ``acknowledge_frozen_ht`` is set; use it for non-magnetic workpieces,
    or between EM re-solves of a staggered EM-thermal run.
    """
    from ngsolve import GridFunction, H1

    if acknowledge_frozen_ht is not True:
        raise ValueError(FROZEN_HT_REFUSAL)
    consistency_tolerance = check_tolerance(consistency_tolerance,
                                            "EM table tolerance")
    if int(azimuths) < 0:
        raise ValueError("azimuths must be >= 0")
    if source.q_scale != 1.0:
        raise ValueError("with a temperature-dependent source, scale the "
                         "current with --ht-scale (|H_t| ~ I), not --q-scale")
    table = ImpedanceTable(load_impedance_table(table_path),
                           allow_extrapolation=allow_extrapolation)
    freq_table = float(table.table["meta"].get("frequency", 0.0))
    audit = {"table_frequency_Hz": freq_table}
    em = source.surface
    if not ht_sol:
        raise ValueError("a temperature-dependent source needs the EM run's "
                         "|H_t| field (--ht-sol <stem>_Ht.sol)")
    pair = verify_field_pair(ht_sol, source.em_vol, source.em_mesh, 1,
                             quantity=HT_QUANTITY)
    gf_h = GridFunction(H1(source.em_mesh, order=1))
    gf_h.Load(os.path.abspath(ht_sol))
    H_v = np.asarray(gf_h.vec.FV().NumPy(), float)[:source.em_mesh.nv]
    rec = read_field_sidecar(ht_sol)
    f_ht = rec.get("frequency_Hz")
    if f_ht is None or not freq_table:
        raise ValueError("both the impedance table and the |H_t| sidecar must "
                         "record their frequency")
    if abs(f_ht - freq_table) > 1e-6 * f_ht:
        raise ValueError(f"impedance table is for {freq_table:g} Hz but "
                         f"the |H_t| field was solved at {f_ht:g} Hz")
    q_tab = table.q(H_v, np.full(H_v.shape, float(T_ref)))
    P_tab = p1_surface_power(em.points, em.tris, q_tab)
    audit["consistency"] = power_balance(
        P_tab, source.power_flat, consistency_tolerance,
        what="impedance table at the reference temperature against the "
             "EM power",
        hint="The table must be built for the material, frequency and "
             "|H_t| range of the EM run (same sigma, BH curve).")
    audit["H_t_source"] = {"file": os.path.abspath(ht_sol)
                           .replace("\\", "/"), **pair}
    h_field = SurfaceP1Field(em.points, em.tris, H_v)
    xyz = np.asarray(target_xyz, float)
    if azimuths and azimuths > 0:
        cols = []
        for k in range(int(azimuths)):
            pts = rotate_about_axis(xyz, 2.0 * math.pi * k / azimuths,
                                    source.axis)
            vals, _ = h_field.evaluate(pts, max_distance=source.transfer_tolerance,
                                       what="ring samples of |H_t|")
            cols.append(vals)
        H_samples = np.stack(cols, axis=1)
    else:
        H_samples, _ = h_field.evaluate(xyz, max_distance=source.transfer_tolerance,
                                        what="thermal heat-flux vertices")
        H_samples = H_samples[:, None]
    gf_damp = GridFunction(gf_target.space)
    gf_damp.vec[:] = 0.0
    coupled = TemperatureDependentSource(
        table, q_ref, H_samples, dofs, gf_target, T_ref=T_ref,
        H_scale=H_scale, gf_damping=gf_damp, acknowledge_frozen_ht=True)
    audit.update({"H_t_range_A_m": [float(H_samples.min()),
                                    float(H_samples.max())]})
    return coupled, audit


# ---------------------------------------------------------------------------
# Command line: explicit sidecar for a legacy field file
# ---------------------------------------------------------------------------

_QUANTITIES = {QSURF_QUANTITY: QSURF_UNIT, HT_QUANTITY: HT_UNIT,
               TEMPERATURE_QUANTITY: TEMPERATURE_UNIT}


def main(argv=None) -> int:
    """``python -m radia.ih_thermal sidecar ...``

    Writes the sidecar of an existing ``.sol`` whose pairing the user
    states: mesh file, H1 order, quantity and, for a heat source, the
    heated EM boundaries and the EM power from the EM run's result.
    """
    import argparse
    from ngsolve import Mesh

    ap = argparse.ArgumentParser(prog="python -m radia.ih_thermal",
                                 description=main.__doc__)
    sub = ap.add_subparsers(dest="command", required=True)
    sc = sub.add_parser("sidecar", help="write <sol>.json for a legacy .sol")
    sc.add_argument("--sol", required=True)
    sc.add_argument("--mesh", required=True, help="the .vol it was saved on")
    sc.add_argument("--order", type=int, required=True)
    sc.add_argument("--quantity", required=True, choices=sorted(_QUANTITIES))
    sc.add_argument("--boundaries", default="",
                    help="heated EM boundaries, '|'-separated (heat source "
                         "and |H_t|)")
    sc.add_argument("--p-wp", type=float, default=None,
                    help="EM power [W] from the EM run (heat source)")
    sc.add_argument("--frequency", type=float, default=None,
                    help="EM frequency [Hz] (|H_t|)")
    sc.add_argument("--definedon", default=None,
                    help="region of a definedon temperature space")
    sc.add_argument("--geometry", default=None,
                    choices=("3d", "axisymmetric-rz"),
                    help="temperature fields: 3d or axisymmetric-rz")
    a = ap.parse_args(argv)
    names = [b for b in a.boundaries.split("|") if b.strip()]
    extra = {"producer": "radia.ih_thermal sidecar (user-stated pairing)"}
    if a.quantity in (QSURF_QUANTITY, HT_QUANTITY) and not names:
        ap.error("--boundaries is required for a surface field")
    if a.quantity == QSURF_QUANTITY:
        if a.p_wp is None:
            ap.error("--p-wp (the EM run's workpiece power) is required")
        extra["P_wp_W"] = float(a.p_wp)
    if a.quantity == TEMPERATURE_QUANTITY:
        if a.geometry is None:
            ap.error("--geometry is required for a temperature field")
        extra["geometry"] = a.geometry
    if a.quantity == HT_QUANTITY:
        if a.frequency is None:
            ap.error("--frequency is required for |H_t|")
        extra["frequency_Hz"] = float(a.frequency)
    path = write_field_sidecar(
        a.sol, mesh_path=a.mesh, mesh=Mesh(a.mesh), fes_order=a.order,
        quantity=a.quantity, unit=_QUANTITIES[a.quantity], boundaries=names,
        extra=extra, definedon=a.definedon)
    print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
