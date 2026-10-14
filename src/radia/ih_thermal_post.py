"""Post-processing of induction-heating temperature fields.

* :func:`thermal_exposure` -- how much material exceeds given temperatures,
  where it is, and how uniform the hottest ring is around the axis.  The
  volumes are quadrature sums of the indicator ``T > T_limit`` over the
  element integration points, so a thin overheated skin is counted even
  when no vertex lies inside it.
* :func:`case_depth` -- depth below the surface at which the temperature
  falls through a threshold (for example the austenitising temperature
  that sets the hardened case), marched along the inward normal from
  stations on the surface.

Both work on 3D meshes and on axisymmetric (r, z) meshes.  Temperatures are
evaluated with the mesh exactly as loaded (see
:func:`radia.ih_thermal.load_field`); nothing here calls ``Mesh.Curve``.
"""
from __future__ import annotations

import math
from typing import Iterable, Sequence

import numpy as np

from . import ih_thermal as _it


def _ref_rule(element_type, order):
    from ngsolve import IntegrationRule
    rule = IntegrationRule(element_type, int(order))
    w = np.asarray(list(rule.weights), float)
    return rule, w / w.sum()


def _region_mask(mesh, region):
    """Boolean mask over volume elements belonging to ``region`` (or all)."""
    if region is None:
        return np.ones(mesh.ne, dtype=bool)
    import re
    pat = re.compile(str(region))
    mask = np.asarray([bool(pat.fullmatch(el.mat)) for el in mesh.Elements()])
    if not mask.any():
        raise ValueError(f"region {region!r} matches no volume element")
    return mask


def field_samples(mesh, gf, *, order: int | None = None,
                  axisymmetric: bool = False, region: str | None = None
                  ) -> dict:
    """Values, volume weights and coordinates at element integration points.

    Weights are the element measure (volume, or revolved volume ``2 pi r``
    for an axisymmetric mesh) distributed by the reference quadrature
    weights, so they sum to the domain measure for straight elements.
    """
    from ngsolve import CF, Integrate, VOL, x, y, z

    if order is None:
        order = max(6, 2 * int(gf.space.globalorder) + 2)
    measure = np.asarray(Integrate(CF(1.0), mesh, VOL, element_wise=True),
                         float)
    xyz_cf = CF((x, y, z)) if mesh.dim == 3 else CF((x, y))
    types = {el.type for el in mesh.Elements(VOL)}
    if len(types) == 1:
        rule, wref = _ref_rule(types.pop(), order)
        mips = mesh.MapToAllElements(rule, VOL)
        nq = len(wref)
        v = np.real(np.asarray(gf(mips))).reshape(-1)
        X = np.asarray(xyz_cf(mips), float).reshape(-1, mesh.dim)
        elem = np.repeat(np.arange(mesh.ne), nq)
        if v.size != mesh.ne * nq:
            raise RuntimeError("unexpected integration point count")
        w = np.tile(wref, mesh.ne) * measure[elem]
    else:
        vs, Xs, es, ws = [], [], [], []
        rules = {}
        for el in mesh.Elements(VOL):
            if el.type not in rules:
                rules[el.type] = _ref_rule(el.type, order)
            rule, wref = rules[el.type]
            mir = mesh.GetTrafo(el)(rule)
            vs.append(np.real(np.asarray(gf(mir))).reshape(-1))
            Xs.append(np.asarray(xyz_cf(mir), float).reshape(-1, mesh.dim))
            es.append(np.full(len(wref), el.nr))
            ws.append(wref * measure[el.nr])
        v, X = np.concatenate(vs), np.concatenate(Xs)
        elem, w = np.concatenate(es), np.concatenate(ws)
    if axisymmetric:
        # measure[] is the (r, z) area; revolve each point with its radius.
        w = w * 2.0 * math.pi * X[:, 0]
    if region is not None:
        keep = _region_mask(mesh, region)[elem]
        v, w, X, elem = v[keep], w[keep], X[keep], elem[keep]
    return {"values": v, "weights": w, "coords": X, "element": elem,
            "order": int(order)}


def _boundary_peak(mesh, gf, order, vertices=None):
    """Largest value and its location over boundary integration points
    (of the boundary elements whose vertices all lie in ``vertices``)."""
    from ngsolve import BND, CF, x, y, z

    xyz_cf = CF((x, y, z)) if mesh.dim == 3 else CF((x, y))
    best, where = -np.inf, None
    rules = {}
    for el in mesh.Elements(BND):
        if vertices is not None and not all(vertices[v.nr]
                                            for v in el.vertices):
            continue
        if el.type not in rules:
            rules[el.type] = _ref_rule(el.type, order)[0]
        mir = mesh.GetTrafo(el)(rules[el.type])
        v = np.real(np.asarray(gf(mir))).reshape(-1)
        k = int(np.argmax(v))
        if v[k] > best:
            best = float(v[k])
            where = np.asarray(xyz_cf(mir), float).reshape(-1, mesh.dim)[k]
    return best, where


def _vertex_values(mesh, gf):
    """Nodal values of an H1 field at every mesh vertex."""
    from ngsolve import NodeId, VERTEX
    dofs = np.asarray([gf.space.GetDofNrs(NodeId(VERTEX, v.nr))[0]
                       for v in mesh.vertices])
    return np.asarray(gf.vec.FV().NumPy(), float)[dofs]


def _require_region_for_definedon(mesh, gf, region):
    """A field on a subdomain has meaningless zeros elsewhere: its region
    must be named."""
    from ngsolve import VOL
    defined = gf.space.GetDefinedOnRegion(VOL)
    # NGSolve returns an all-zero sentinel for an unrestricted space.
    mask = list(defined.Mask()) if defined is not None else []
    if region is None and any(mask) and not all(mask):
        raise ValueError(
            "the field is defined on a subdomain (definedon space); pass its "
            "region so the rest of the mesh is not counted as material")


def thermal_exposure(mesh, gf, thresholds_C: Iterable[float] = (), *,
                     axisymmetric: bool = False, axis: str = "z",
                     order: int | None = None, ring_samples: int = 72,
                     samples: dict | None = None,
                     region: str | None = None) -> dict:
    """Summarise where a temperature field exceeds given limits.

    Returns a JSON-serialisable dictionary with the peak temperature and
    its location, and for each threshold the volume above it, its fraction
    of the domain and its bounding box (xyz, and meridian (r, a) for 3D
    meshes about ``axis``).  For a 3D mesh the ring through the hottest
    point is sampled on ``ring_samples`` azimuths and its spread reported:
    a large spread under an axisymmetric heat source means the source was
    not axisymmetric.
    """
    _require_region_for_definedon(mesh, gf, region)
    # The threshold indicator is discontinuous even for a polynomial
    # temperature. Its sampling needs more points than integration of T.
    # This remains a quadrature estimate; use explicit order refinement
    # when a thin layer's volume is an acceptance criterion.
    if order is None and samples is None:
        order = max(12, 2 * int(gf.space.globalorder) + 2)
    s = samples if samples is not None else field_samples(
        mesh, gf, order=order, axisymmetric=axisymmetric, region=region)
    vals, w, X = s["values"], s["weights"], s["coords"]
    if not np.all(np.isfinite(vals)):
        raise RuntimeError("the temperature field contains non-finite values")
    pts = _it.mesh_vertices(mesh)[:, :mesh.dim]
    vv = _vertex_values(mesh, gf)
    in_region = None
    if region is not None:
        mask = _region_mask(mesh, region)
        in_region = np.zeros(mesh.nv, dtype=bool)
        for el in mesh.Elements():
            if mask[el.nr]:
                for vtx in el.vertices:
                    in_region[vtx.nr] = True
        pts, vv = pts[in_region], vv[in_region]
    i_s, i_v = int(np.argmax(vals)), int(np.argmax(vv))
    candidates = [(float(vals[i_s]), X[i_s]), (float(vv[i_v]), pts[i_v])]
    if gf.space.globalorder > 1:
        candidates.append(_boundary_peak(mesh, gf, s["order"], in_region))
    t_max, loc = max(candidates, key=lambda c: c[0])
    total = float(w.sum())
    out = {
        "method": "integration-point-indicator",
        "region": region,
        "integration_order": s["order"],
        "axisymmetric": bool(axisymmetric),
        "domain_measure_m3": total,
        "T_max_C": t_max,
        "T_max_location_m": [float(c) for c in loc],
        "T_min_C": float(min(vals.min(), vv.min())),
        "T_mean_C": float(np.sum(vals * w) / total),
        "thresholds": [],
    }
    loc3 = np.r_[loc, np.zeros(3 - len(loc))]
    if not axisymmetric and mesh.dim == 3:
        out["T_max_meridian_m"] = [float(c) for c in
                                   _it.to_meridian(loc3[None, :], axis)[0]]
    for thr in sorted(float(t) for t in thresholds_C):
        hot = vals > thr
        entry = {"T_C": thr,
                 "volume_m3": float(w[hot].sum()),
                 "volume_fraction": float(w[hot].sum() / total),
                 "n_elements": int(np.unique(s["element"][hot]).size)}
        if np.any(hot):
            entry["bbox_min_m"] = [float(c) for c in X[hot].min(axis=0)]
            entry["bbox_max_m"] = [float(c) for c in X[hot].max(axis=0)]
            if not axisymmetric and mesh.dim == 3:
                rz = _it.to_meridian(X[hot], axis)
                entry["meridian_min_m"] = [float(c) for c in rz.min(axis=0)]
                entry["meridian_max_m"] = [float(c) for c in rz.max(axis=0)]
        out["thresholds"].append(entry)
    if not axisymmetric and mesh.dim == 3 and ring_samples > 0:
        out["hottest_ring"] = _ring_spread(mesh, gf, loc3, axis, ring_samples)
    return out


def _ring_spread(mesh, gf, point, axis, n):
    """Temperature on the ring through ``point`` about ``axis``.

    The ring is taken a quarter of the median boundary edge radially inside
    ``point``: rotated copies of a surface point lie off a faceted surface,
    so a ring exactly on it would leave the mesh.  Points of the ring that
    are still outside the mesh are counted, not replaced.
    """
    from ngsolve import BND

    iu, iw, ia = _it._axis_frame(axis)
    r = math.hypot(point[iu], point[iw])
    tris = np.asarray([[v.nr for v in el.vertices][:3]
                       for el in mesh.Elements(BND)], int)
    inset = min(0.25 * _it.median_edge_length(_it.mesh_vertices(mesh), tris),
                0.5 * r)
    q0 = np.array(point, float)
    if r > 0:
        q0[iu] *= (r - inset) / r
        q0[iw] *= (r - inset) / r
    ring = np.concatenate([_it.rotate_about_axis(q0[None, :],
                                                 2 * math.pi * k / n, axis)
                           for k in range(n)])
    vals = _evaluate_points(mesh, gf, ring)
    ok = ~np.isnan(vals)
    rec = {"r_m": r, "ring_r_m": r - inset, "a_m": float(point[ia]),
           "n_requested": n, "n_inside": int(ok.sum())}
    if np.any(ok):
        v = vals[ok]
        rec.update({"T_min_C": float(v.min()), "T_max_C": float(v.max()),
                    "T_mean_C": float(v.mean()),
                    "spread_C": float(v.max() - v.min())})
    return rec


def limit_check(exposure: dict, limit_C: float) -> dict:
    """Constraint record for an optimiser: violation volume above a limit."""
    for e in exposure["thresholds"]:
        if e["T_C"] == float(limit_C):
            return {"limit_C": float(limit_C),
                    "exceeded": e["volume_m3"] > 0.0,
                    "excess_C": max(0.0, exposure["T_max_C"] - float(limit_C)),
                    "volume_above_m3": e["volume_m3"]}
    raise ValueError(f"limit {limit_C} C was not among the thresholds")


# ---------------------------------------------------------------------------
# Case depth
# ---------------------------------------------------------------------------

DEPTH_STATUS = ("ok", "not_reached", "through", "beyond_span", "edge")


def _evaluate_points(mesh, gf, pts, region_mask=None):
    """Field values at ``pts`` (n, dim); NaN where a point is outside the
    mesh or outside the region."""
    cols = [np.ascontiguousarray(pts[:, k]) for k in range(pts.shape[1])]
    mips = mesh(*cols)
    inside = mips["nr"] >= 0
    if region_mask is not None:
        nr = mips["nr"]
        inside &= np.where(nr >= 0, region_mask[np.maximum(nr, 0)], False)
    out = np.full(len(pts), np.nan)
    if np.any(inside):
        out[inside] = np.real(np.asarray(gf(mips[inside]))).reshape(-1)
    return out


def _region_boundary(mesh, rmask):
    """Boundary facets of the region: triangles (3D, on the mesh boundary)
    or segments (2D, including interfaces with other regions)."""
    from collections import Counter
    from ngsolve import BND
    pts = _it.mesh_vertices(mesh)[:, :mesh.dim]
    region_vertices = set()
    for el in mesh.Elements():
        if rmask is None or rmask[el.nr]:
            region_vertices.update(v.nr for v in el.vertices)
    if mesh.dim == 3:
        facets = [[v.nr for v in el.vertices] for el in mesh.Elements(BND)
                  if all(v.nr in region_vertices for v in el.vertices)]
        tris = []
        for f in facets:
            tris.append(f[:3])
            if len(f) == 4:
                tris.append([f[0], f[2], f[3]])
        return pts, np.asarray(tris, int)
    edges = Counter()
    for el in mesh.Elements():
        if rmask is not None and not rmask[el.nr]:
            continue
        vv = [v.nr for v in el.vertices]
        for i in range(len(vv)):
            edges[tuple(sorted((vv[i], vv[(i + 1) % len(vv)])))] += 1
    return pts, np.asarray([e for e, n in edges.items() if n == 1], int)


def _project_to_boundary(query, pts, facets):
    """Closest points on facets (triangles in 3D, segments in 2D)."""
    q = np.asarray(query, float)
    if facets.shape[1] == 3:
        field = _it.SurfaceP1Field(np.c_[pts, np.zeros((len(pts), 3 - pts.shape[1]))]
                                   if pts.shape[1] < 3 else pts,
                                   facets, np.zeros(len(pts)))
        tri, lam, d = field.locate(q)
        a, b, c = field._abc
        cp = lam[:, :1] * a[tri] + lam[:, 1:2] * b[tri] + lam[:, 2:] * c[tri]
        return cp, d
    a, b = pts[facets[:, 0]], pts[facets[:, 1]]
    best_d = np.full(len(q), np.inf)
    best = np.zeros_like(q)
    for lo in range(0, len(q), 256):
        qq = q[lo:lo + 256]
        ab = b - a
        t = np.einsum("ijk,jk->ij", qq[:, None, :] - a[None], ab) / \
            np.maximum(np.einsum("ij,ij->i", ab, ab), 1e-300)[None]
        t = np.clip(t, 0.0, 1.0)
        cp = a[None] + t[..., None] * ab[None]
        d = np.linalg.norm(qq[:, None, :] - cp, axis=2)
        j = np.argmin(d, axis=1)
        best_d[lo:lo + 256] = d[np.arange(len(qq)), j]
        best[lo:lo + 256] = cp[np.arange(len(qq)), j]
    return best, best_d


def _march(mesh, gf, origins, normals, threshold, span, step,
           region_mask=None):
    """Depth of the first fall through ``threshold`` along each ray.

    Samples sit at ``(k + 1/2) step`` below the (projected) surface point.
    A ray whose first sample is not in the region points out of the
    material; that is an error, not a skipped sample.
    """
    ts = (np.arange(int(math.ceil(span / step))) + 0.5) * step
    n_st = len(origins)
    P = (origins[:, None, :] + ts[None, :, None] * normals[:, None, :])
    T = _evaluate_points(mesh, gf, P.reshape(-1, origins.shape[1]),
                         region_mask).reshape(n_st, len(ts))
    bad = np.isnan(T[:, 0])
    if np.any(bad):
        i = int(np.flatnonzero(bad)[0])
        raise ValueError(
            f"{int(bad.sum())}/{n_st} depth rays leave the material at their "
            f"first sample (station {np.round(origins[i], 6).tolist()}, normal "
            f"{np.round(normals[i], 3).tolist()}); the normal must point into "
            "the region")
    depth = np.full(n_st, np.nan)
    status = np.empty(n_st, dtype=object)
    surface = T[:, 0].copy()
    reentrant = np.zeros(n_st, dtype=bool)
    for i in range(n_st):
        row = T[i]
        gap = np.flatnonzero(np.isnan(row))
        last = gap[0] - 1 if gap.size else len(row) - 1
        seg = row[:last + 1]
        if seg[0] < threshold:
            depth[i], status[i] = 0.0, "not_reached"
            continue
        below = np.flatnonzero(seg < threshold)
        if below.size == 0:
            depth[i] = ts[last]
            status[i] = "through" if gap.size else "beyond_span"
            continue
        j = below[0]
        t0, t1 = ts[j - 1], ts[j]
        v0, v1 = seg[j - 1], seg[j]
        depth[i] = t0 + (t1 - t0) * (v0 - threshold) / (v0 - v1)
        status[i] = "ok"
        reentrant[i] = bool(np.any(seg[j:] >= threshold))
    # A ray that leaves the material (the far surface, or a hole) still hot
    # is reported at the exit point, found by bisection between the last
    # sample inside and the first outside.
    out = np.flatnonzero(status == "through")
    if out.size:
        lo = depth[out].copy()
        hi = lo + step
        for _ in range(40):
            mid = 0.5 * (lo + hi)
            v = _evaluate_points(mesh, gf, origins[out] + mid[:, None]
                                 * normals[out], region_mask)
            inside = ~np.isnan(v)
            lo = np.where(inside, mid, lo)
            hi = np.where(inside, hi, mid)
        depth[out] = 0.5 * (lo + hi)
    return depth, status, surface, reentrant


def surface_stations(mesh, boundary_names: Sequence[str],
                     region_mask=None, edge_angle_deg: float = 30.0):
    """Stations on the named boundaries with unit inward normals.

    3D: every boundary vertex, normal = area-weighted mean of the adjacent
    triangle normals; a vertex whose adjacent face normals differ from that
    mean by more than ``edge_angle_deg`` sits on an edge or corner, where a
    single ray does not measure the depth -- it is returned flagged.
    2D (r, z): every boundary segment midpoint.  The sign is chosen so that
    a short step along the normal stays in the region.
    Returns ``(origins, normals, h, edge_flags)``.
    """
    pts = _it.mesh_vertices(mesh)[:, :mesh.dim]
    if mesh.dim == 3:
        tris = _it.boundary_triangles(mesh, boundary_names)
        a, b, c = pts[tris[:, 0]], pts[tris[:, 1]], pts[tris[:, 2]]
        fn = np.cross(b - a, c - a)
        acc = np.zeros_like(pts)
        for k in range(3):
            np.add.at(acc, tris[:, k], fn)
        idx = np.unique(tris)
        mean = acc / np.maximum(np.linalg.norm(acc, axis=1), 1e-300)[:, None]
        unit_fn = fn / np.linalg.norm(fn, axis=1)[:, None]
        worst = np.ones(len(pts))
        for k in range(3):
            cosang = np.einsum("ij,ij->i", unit_fn, mean[tris[:, k]])
            np.minimum.at(worst, tris[:, k], cosang)
        edge = worst[idx] < math.cos(math.radians(edge_angle_deg))
        origins, normals = pts[idx], acc[idx]
        h = _it.median_edge_length(pts, tris)
    else:
        segs = _it.boundary_edges_2d(mesh, boundary_names)
        a, b = pts[segs[:, 0]], pts[segs[:, 1]]
        origins = 0.5 * (a + b)
        d = b - a
        normals = np.c_[-d[:, 1], d[:, 0]]
        h = float(np.median(np.linalg.norm(d, axis=1)))
        edge = np.zeros(len(origins), dtype=bool)
    normals = normals / np.linalg.norm(normals, axis=1)[:, None]
    probe = origins + 0.25 * h * normals
    cols = [np.ascontiguousarray(probe[:, k]) for k in range(mesh.dim)]
    nr = mesh(*cols)["nr"]
    inside = nr >= 0
    if region_mask is not None:
        inside &= np.where(nr >= 0, region_mask[np.maximum(nr, 0)], False)
    normals[~inside] *= -1.0
    return origins, normals, h, edge


def case_depth(mesh, gf, threshold_C: float, *, span: float,
               origins=None, normals=None,
               boundary_names: Sequence[str] | None = None,
               step: float | None = None, axis: str = "z", n_phi: int = 0,
               region: str | None = None, edge_angle_deg: float = 30.0,
               projection_tolerance: float | None = None) -> dict:
    """Depth at which the temperature falls below ``threshold_C``.

    Stations are either given (``origins``, ``normals``; unit inward
    normals in mesh coordinates) or taken from ``boundary_names`` with
    :func:`surface_stations`.  For a 3D mesh ``origins`` may instead be
    meridian points ``(r, a)`` with meridian normals ``(n_r, n_a)`` about
    ``axis``, marched on ``n_phi`` azimuths each.  Every origin is projected
    onto the region boundary; a projection farther than
    ``projection_tolerance`` (default a quarter of the median boundary edge)
    is an error, and the projection distances are reported.

    Rays are marched inward for ``span`` in steps of ``step``.  Status:

    ``ok``          the temperature falls through the threshold at ``depth``;
    ``not_reached`` the first sample is below the threshold (depth 0);
    ``through``     the ray leaves the material still above the threshold
                    (depth = distance to the far surface or hole wall);
    ``beyond_span`` still above the threshold after ``span`` inside the
                    material (depth is only a lower bound);
    ``edge``        an automatic station on an edge or corner (no depth).
    """
    if (origins is None) != (normals is None):
        raise ValueError("pass both origins and normals, or neither")
    if not (span and span > 0):
        raise ValueError("span (probe length) must be positive; choose it "
                         "longer than the expected depth ('beyond_span' marks "
                         "rays it cut short)")
    _require_region_for_definedon(mesh, gf, region)
    rmask = _region_mask(mesh, region) if region is not None else None
    meridian = False
    edge = None
    if origins is None:
        if not boundary_names:
            raise ValueError("give stations or boundary_names")
        origins, normals, h, edge = surface_stations(
            mesh, boundary_names, rmask, edge_angle_deg)
    else:
        origins = np.asarray(origins, float)
        normals = np.asarray(normals, float)
        nrm = np.linalg.norm(normals, axis=1)
        if np.any(nrm == 0):
            raise ValueError("station normals must be non-zero")
        normals = normals / nrm[:, None]
        if origins.shape != normals.shape or origins.shape[1] not in (
                (2,) if mesh.dim == 2 else (2, 3)):
            raise ValueError(
                f"stations of shape {origins.shape}/{normals.shape} do not fit "
                f"a {mesh.dim}D mesh (a 2D (r, z) mesh takes (r, z) stations)")
        meridian = mesh.dim == 3 and origins.shape[1] == 2
    bpts, facets = _region_boundary(mesh, rmask)
    fe = bpts[facets[:, 1]] - bpts[facets[:, 0]]
    h = float(np.median(np.linalg.norm(fe, axis=1)))
    tol = 0.25 * h if projection_tolerance is None else float(
        projection_tolerance)
    step = float(step) if step else min(h / 8.0, span / 50.0)
    n_st = len(origins)
    if meridian:
        if n_phi < 1:
            raise ValueError("meridian stations need n_phi >= 1")
        iu, iw, ia = _it._axis_frame(axis)
        o3 = np.zeros((n_st, 3))
        n3 = np.zeros((n_st, 3))
        o3[:, iu], o3[:, ia] = origins[:, 0], origins[:, 1]
        n3[:, iu], n3[:, ia] = normals[:, 0], normals[:, 1]
        O = np.concatenate([_it.rotate_about_axis(o3, 2 * math.pi * k / n_phi,
                                                  axis) for k in range(n_phi)])
        N = np.concatenate([_it.rotate_about_axis(n3, 2 * math.pi * k / n_phi,
                                                  axis) for k in range(n_phi)])
    else:
        O, N = origins, normals
    Op, dproj = _project_to_boundary(O, bpts, facets)
    if np.any(dproj > tol):
        i = int(np.argmax(dproj))
        raise ValueError(
            f"{int((dproj > tol).sum())}/{len(O)} depth stations lie farther "
            f"than {tol:.3e} m from the material surface (worst {dproj[i]:.3e}"
            f" m at {np.round(O[i], 6).tolist()})")
    keep = np.ones(len(Op), dtype=bool) if edge is None else ~edge
    depth = np.full(len(Op), np.nan)
    status = np.full(len(Op), "edge", dtype=object)
    surf = np.full(len(Op), np.nan)
    reent = np.zeros(len(Op), dtype=bool)
    if np.any(keep):
        d_, s_, t_, r_ = _march(mesh, gf, Op[keep], N[keep],
                                float(threshold_C), float(span), step, rmask)
        depth[keep], status[keep], surf[keep], reent[keep] = d_, s_, t_, r_
    out = {"threshold_C": float(threshold_C), "span_m": float(span),
           "step_m": step, "n_stations": n_st, "region": region,
           "projection_tolerance_m": tol,
           "max_projection_distance_m": float(dproj.max()),
           "origins": origins.tolist(), "normals": normals.tolist()}
    if meridian:
        D = depth.reshape(n_phi, n_st).T
        S = status.reshape(n_phi, n_st).T
        out.update({
            "n_phi": int(n_phi),
            "depth_min_m": np.min(D, axis=1).tolist(),
            "depth_mean_m": np.mean(D, axis=1).tolist(),
            "depth_max_m": np.max(D, axis=1).tolist(),
            "depth_m": D.tolist(),
            "status": [[str(x) for x in row] for row in S],
            "status_counts": {k: int(np.sum(S == k)) for k in DEPTH_STATUS},
            "surface_T_C": surf.reshape(n_phi, n_st).T.tolist(),
        })
    else:
        out.update({
            "depth_m": depth.tolist(),
            "status": [str(x) for x in status],
            "status_counts": {k: int(np.sum(status == k))
                              for k in DEPTH_STATUS},
            "surface_T_C": surf.tolist(),
            "reentrant": reent.tolist(),
        })
    return out

# ---------------------------------------------------------------------------
# Command line
# ---------------------------------------------------------------------------

_STATION_LAYOUTS = {
    ("x", "y", "z", "nx", "ny", "nz"): "xyz",
    ("r", "z", "n_r", "n_z"): "meridian",
}


def _read_stations(path: str):
    """Stations CSV.  The header states the layout and the length unit:
    ``x_m,y_m,z_m,nx,ny,nz`` or ``r_m,z_m,n_r,n_z`` (``_mm`` for
    millimetres), optionally preceded by ``s_m`` / ``s_mm``."""
    import csv

    with open(path, encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.reader(handle))
    header = [h.strip() for h in rows[0]]
    data = np.asarray([[float(v) for v in r] for r in rows[1:] if r])
    if data.size == 0:
        raise ValueError(f"{path} has no stations")
    s_col = None
    if header[0] in ("s_m", "s_mm"):
        s_col = data[:, 0] * (1e-3 if header[0] == "s_mm" else 1.0)
        header, data = header[1:], data[:, 1:]
    units = {h.rsplit("_", 1)[1] for h in header
             if h.endswith(("_m", "_mm"))}
    if len(units) != 1:
        raise ValueError(f"{path}: coordinate columns must all be _m or all "
                         f"_mm, got {header}")
    scale = 1e-3 if units == {"mm"} else 1.0
    bare = tuple(h[:-3] if h.endswith("_mm") else h[:-2]
                 if h.endswith("_m") else h for h in header)
    if bare not in _STATION_LAYOUTS:
        raise ValueError(f"{path}: header must be x_m,y_m,z_m,nx,ny,nz or "
                         f"r_m,z_m,n_r,n_z (or _mm), got {header}")
    nc = 3 if _STATION_LAYOUTS[bare] == "xyz" else 2
    return data[:, :nc] * scale, data[:, nc:], s_col


def main(argv=None) -> int:
    import argparse
    import json

    ap = argparse.ArgumentParser(
        prog="python -m radia.ih_thermal_post",
        description="Exposure and case depth of a saved IH temperature "
                    "field (<T>.sol with its .sol.json sidecar).")
    ap.add_argument("command", choices=("exposure", "depth"))
    ap.add_argument("--temperature", required=True, help="temperature .sol")
    ap.add_argument("--mesh", default=None,
                    help="relocated mesh .vol (its digest must match the "
                         "sidecar)")
    ap.add_argument("--axis", default="z", choices=("x", "y", "z"))
    ap.add_argument("--thresholds", default="",
                    help="exposure: comma-separated temperatures [degC]")
    ap.add_argument("--threshold", type=float, default=None,
                    help="depth: threshold temperature [degC]")
    ap.add_argument("--stations", default="",
                    help="depth: stations CSV with header x_m,y_m,z_m,nx,ny,nz"
                         " or r_m,z_m,n_r,n_z (_mm for millimetres); meridian "
                         "stations are marched on --n-phi azimuths")
    ap.add_argument("--boundaries", default="",
                    help="depth: '|'-separated boundary names for automatic "
                         "stations")
    ap.add_argument("--span", type=float, default=None,
                    help="depth: probe length [m]")
    ap.add_argument("--step", type=float, default=None,
                    help="depth: probe step [m]")
    ap.add_argument("--n-phi", type=int, default=0,
                    help="depth: azimuths per meridian station")
    ap.add_argument("--output", default="", help="JSON output file")
    ap.add_argument("--csv", default="", help="depth: CSV output file")
    a = ap.parse_args(argv)

    mesh, gf, audit = _it.load_field(a.temperature,
                                     quantity=_it.TEMPERATURE_QUANTITY,
                                     mesh_path=a.mesh)
    record = audit.pop("sidecar_record")
    region = record.get("definedon")
    geometry = record.get("geometry")
    if geometry not in ("3d", "axisymmetric-rz"):
        raise SystemExit(
            f"{a.temperature}: the sidecar records no geometry ('3d' or "
            "'axisymmetric-rz'); rewrite it with 'python -m radia.ih_thermal "
            "sidecar ... --geometry'")
    axisymmetric = geometry == "axisymmetric-rz"
    if mesh.dim != (2 if axisymmetric else 3):
        raise SystemExit(f"{a.temperature}: the sidecar records geometry "
                         f"{geometry!r} but the mesh is {mesh.dim}D")
    if a.command == "exposure":
        thr = [float(t) for t in a.thresholds.split(",") if t.strip()]
        result = thermal_exposure(mesh, gf, thr, axisymmetric=axisymmetric,
                                  axis=a.axis, region=region)
    else:
        if a.threshold is None:
            ap.error("depth needs --threshold")
        if a.span is None:
            ap.error("depth needs --span")
        s_col = None
        if a.stations:
            o, n, s_col = _read_stations(a.stations)
            result = case_depth(mesh, gf, a.threshold, origins=o, normals=n,
                                span=a.span, step=a.step, axis=a.axis,
                                n_phi=a.n_phi, region=region)
        else:
            names = [b for b in a.boundaries.split("|") if b.strip()]
            result = case_depth(mesh, gf, a.threshold, boundary_names=names,
                                span=a.span, step=a.step, region=region)
        if a.csv:
            _write_depth_csv(a.csv, result, s_col)
    result["field"] = audit
    text = json.dumps(result, indent=1, ensure_ascii=False)
    if a.output:
        with open(a.output, "w", encoding="utf-8") as handle:
            handle.write(text)
    else:
        print(text)
    return 0


def _write_depth_csv(path, result, s_col):
    import csv

    with open(path, "w", newline="", encoding="utf-8") as handle:
        w = csv.writer(handle)
        if "n_phi" in result:
            w.writerow(["s", "r_m", "a_m", "depth_min_m", "depth_mean_m",
                        "depth_max_m", "statuses"])
            for i, o in enumerate(result["origins"]):
                sts = sorted(set(result["status"][i]))
                w.writerow([s_col[i] if s_col is not None else "", o[0], o[1],
                            result["depth_min_m"][i], result["depth_mean_m"][i],
                            result["depth_max_m"][i], "|".join(sts)])
        else:
            w.writerow(["s", *[f"o{k}" for k in range(len(result["origins"][0]))],
                        "depth_m", "status"])
            for i, o in enumerate(result["origins"]):
                w.writerow([s_col[i] if s_col is not None else "", *o,
                            result["depth_m"][i], result["status"][i]])


if __name__ == "__main__":
    raise SystemExit(main())
