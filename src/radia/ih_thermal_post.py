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

try:                                     # package or panels-style import
    from . import ih_thermal as _it
except ImportError:                      # pragma: no cover
    import ih_thermal as _it


def _ref_rule(element_type, order):
    from ngsolve import IntegrationRule
    rule = IntegrationRule(element_type, int(order))
    w = np.asarray(list(rule.weights), float)
    return rule, w / w.sum()


def field_samples(mesh, gf, *, order: int | None = None,
                  axisymmetric: bool = False) -> dict:
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
    return {"values": v, "weights": w, "coords": X, "element": elem,
            "order": int(order)}


def _boundary_peak(mesh, gf, order):
    """Largest value and its location over boundary integration points."""
    from ngsolve import BND, CF, x, y, z

    xyz_cf = CF((x, y, z)) if mesh.dim == 3 else CF((x, y))
    best, where = -np.inf, None
    rules = {}
    for el in mesh.Elements(BND):
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
    from ngsolve import NodeId, VERTEX
    dofs = [gf.space.GetDofNrs(NodeId(VERTEX, v.nr))[0] for v in mesh.vertices]
    return np.asarray(gf.vec.FV().NumPy(), float)[dofs]


def thermal_exposure(mesh, gf, thresholds_C: Iterable[float] = (), *,
                     axisymmetric: bool = False, axis: str = "z",
                     order: int | None = None, ring_samples: int = 72,
                     samples: dict | None = None) -> dict:
    """Summarise where a temperature field exceeds given limits.

    Returns a JSON-serialisable dictionary with the peak temperature and
    its location, and for each threshold the volume above it, its fraction
    of the domain and its bounding box (xyz, and meridian (r, a) for 3D
    meshes about ``axis``).  For a 3D mesh the ring through the hottest
    point is sampled on ``ring_samples`` azimuths and its spread reported:
    a large spread under an axisymmetric heat source means the source was
    not axisymmetric.
    """
    s = samples if samples is not None else field_samples(
        mesh, gf, order=order, axisymmetric=axisymmetric)
    vals, w, X = s["values"], s["weights"], s["coords"]
    if not np.all(np.isfinite(vals)):
        raise RuntimeError("the temperature field contains non-finite values")
    pts = _it.mesh_vertices(mesh)[:, :mesh.dim]
    vv = _vertex_values(mesh, gf)
    i_s, i_v = int(np.argmax(vals)), int(np.argmax(vv))
    candidates = [(float(vals[i_s]), X[i_s]), (float(vv[i_v]), pts[i_v])]
    if gf.space.globalorder > 1:
        candidates.append(_boundary_peak(mesh, gf, s["order"]))
    t_max, loc = max(candidates, key=lambda c: c[0])
    total = float(w.sum())
    out = {
        "method": "integration-point-indicator",
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
    """Sample the ring through ``point`` about ``axis``.

    On a faceted surface the rotated copies of a surface point lie just
    outside the chords; those are pulled radially inward by up to a few
    chord sags (``h^2 / 8r``), and the largest inset used is reported.
    """
    from ngsolve import BND

    iu, iw, ia = _it._axis_frame(axis)
    r = math.hypot(point[iu], point[iw])
    tris = np.asarray([[v.nr for v in el.vertices][:3]
                       for el in mesh.Elements(BND)], int)
    h = _it.median_edge_length(_it.mesh_vertices(mesh), tris)
    sag = h * h / (8.0 * r) if r > 0 else 0.0
    insets = [0.0] + [min(f * sag, 0.5 * r) for f in (1.0, 2.0, 4.0)]
    vals, used = [], 0.0
    for k in range(n):
        p = _it.rotate_about_axis(point[None, :], 2 * math.pi * k / n, axis)[0]
        for inset in insets:
            q = p.copy()
            if inset and r > 0:
                scale = (r - inset) / r
                q[iu] *= scale
                q[iw] *= scale
            mip = mesh(*q)
            if mip.nr >= 0:
                vals.append(float(np.real(gf(mip))))
                used = max(used, inset)
                break
    rec = {"r_m": r, "a_m": float(point[ia]), "n_requested": n,
           "n_inside": len(vals), "max_radial_inset_m": used}
    if vals:
        rec.update({"T_min_C": min(vals), "T_max_C": max(vals),
                    "T_mean_C": float(np.mean(vals)),
                    "spread_C": max(vals) - min(vals)})
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

DEPTH_STATUS = ("ok", "not_reached", "through", "beyond_span")


def _evaluate_points(mesh, gf, pts):
    """Field values at ``pts`` (n, dim); NaN where a point is outside."""
    cols = [np.ascontiguousarray(pts[:, k]) for k in range(pts.shape[1])]
    mips = mesh(*cols)
    inside = mips["nr"] >= 0
    out = np.full(len(pts), np.nan)
    if np.any(inside):
        out[inside] = np.real(np.asarray(gf(mips[inside]))).reshape(-1)
    return out


def _march(mesh, gf, origins, normals, threshold, span, step, entry):
    """Depth of the first fall through ``threshold`` along each ray."""
    ts = np.arange(0.0, span + 0.5 * step, step)
    n_st = len(origins)
    P = (origins[:, None, :] + ts[None, :, None] * normals[:, None, :])
    T = _evaluate_points(mesh, gf, P.reshape(-1, origins.shape[1])) \
        .reshape(n_st, len(ts))
    depth = np.full(n_st, np.nan)
    status = np.empty(n_st, dtype=object)
    surface = np.full(n_st, np.nan)
    exit_at = np.full(n_st, np.nan)
    reentrant = np.zeros(n_st, dtype=bool)
    entry_steps = int(math.ceil(entry / step))
    for i in range(n_st):
        row = T[i]
        valid = np.flatnonzero(~np.isnan(row))
        if valid.size == 0 or valid[0] > entry_steps:
            status[i] = "off_mesh"
            continue
        first = valid[0]
        # material ends at the first invalid sample after entering it
        gap = np.flatnonzero(np.isnan(row[first:]))
        last = first + (gap[0] - 1 if gap.size else len(row) - first - 1)
        seg = row[first:last + 1]
        surface[i] = seg[0]
        if seg[0] < threshold:
            depth[i], status[i] = 0.0, "not_reached"
            continue
        below = np.flatnonzero(seg < threshold)
        if below.size == 0:
            if gap.size:
                depth[i], status[i] = ts[last], "through"
                exit_at[i] = ts[last]
            else:
                depth[i], status[i] = ts[last], "beyond_span"
            continue
        j = below[0]
        t0, t1 = ts[first + j - 1], ts[first + j]
        v0, v1 = seg[j - 1], seg[j]
        depth[i] = t0 + (t1 - t0) * (v0 - threshold) / (v0 - v1)
        status[i] = "ok"
        reentrant[i] = bool(np.any(seg[j:] >= threshold))
    return depth, status, surface, exit_at, reentrant


def surface_stations(mesh, boundary_names: Sequence[str]):
    """Stations on the named boundaries with unit inward normals.

    3D: every boundary vertex, normal = area-weighted mean of the adjacent
    triangle normals.  2D (r, z): every boundary segment midpoint.  The
    sign is chosen so that a short step along the normal stays in the mesh.
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
        origins, normals = pts[idx], acc[idx]
        h = _it.median_edge_length(pts, tris)
    else:
        segs = _it.boundary_edges_2d(mesh, boundary_names)
        a, b = pts[segs[:, 0]], pts[segs[:, 1]]
        origins = 0.5 * (a + b)
        d = b - a
        normals = np.c_[-d[:, 1], d[:, 0]]
        h = float(np.median(np.linalg.norm(d, axis=1)))
    normals = normals / np.linalg.norm(normals, axis=1)[:, None]
    probe = origins + 0.25 * h * normals
    cols = [np.ascontiguousarray(probe[:, k]) for k in range(mesh.dim)]
    inside = mesh(*cols)["nr"] >= 0
    normals[~inside] *= -1.0
    return origins, normals, h


def case_depth(mesh, gf, threshold_C: float, *, origins=None, normals=None,
               boundary_names: Sequence[str] | None = None,
               span: float | None = None, step: float | None = None,
               axis: str = "z", n_phi: int = 0) -> dict:
    """Depth at which the temperature falls below ``threshold_C``.

    Stations are either given (``origins``, ``normals``; unit inward
    normals in mesh coordinates) or taken from ``boundary_names`` with
    :func:`surface_stations`.  For a 3D mesh ``origins`` may instead be
    meridian points ``(r, a)`` with meridian normals ``(n_r, n_a)`` about
    ``axis``, marched on ``n_phi`` azimuths each.

    Each ray is marched from the surface in steps of ``step`` for at most
    ``span``.  Status per ray:

    ``ok``          the temperature falls through the threshold at ``depth``;
    ``not_reached`` the surface itself is below the threshold (depth 0);
    ``through``     the ray leaves the material still above the threshold
                    (depth = distance to the far surface);
    ``beyond_span`` still above the threshold after ``span`` inside the
                    material (depth is only a lower bound);
    ``off_mesh``    the station could not be located on the mesh.
    """
    if (origins is None) != (normals is None):
        raise ValueError("pass both origins and normals, or neither")
    meridian = False
    if origins is None:
        if not boundary_names:
            raise ValueError("give stations or boundary_names")
        origins, normals, h = surface_stations(mesh, boundary_names)
    else:
        origins = np.asarray(origins, float)
        normals = np.asarray(normals, float)
        nrm = np.linalg.norm(normals, axis=1)
        if np.any(nrm == 0):
            raise ValueError("station normals must be non-zero")
        normals = normals / nrm[:, None]
        meridian = mesh.dim == 3 and origins.shape[1] == 2
        from ngsolve import BND
        tris = [[v.nr for v in el.vertices][:mesh.dim]
                for el in mesh.Elements(BND)]
        pts = _it.mesh_vertices(mesh)[:, :mesh.dim]
        tris = np.asarray(tris, int)
        e = pts[tris[:, 1]] - pts[tris[:, 0]]
        h = float(np.median(np.linalg.norm(e, axis=1)))
    if span is None:
        raise ValueError("span (probe length) is required: choose it longer "
                         "than the expected depth; 'beyond_span' marks rays "
                         "that it cut short")
    step = float(step) if step else min(h / 8.0, span / 50.0)
    entry = 0.5 * h
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
    depth, status, surf, exit_at, reent = _march(
        mesh, gf, O, N, float(threshold_C), float(span), step, entry)
    off = np.flatnonzero(status == "off_mesh")
    if off.size:
        raise ValueError(
            f"{off.size}/{len(status)} depth stations could not be located "
            f"within {entry:.3e} m of the mesh surface (first at "
            f"{np.round(O[off[0]], 6).tolist()})")
    out = {"threshold_C": float(threshold_C), "span_m": float(span),
           "step_m": step, "n_stations": n_st,
           "origins": origins.tolist(), "normals": normals.tolist()}
    if meridian:
        D = depth.reshape(n_phi, n_st).T
        S = status.reshape(n_phi, n_st).T
        out.update({
            "n_phi": int(n_phi),
            "depth_min_m": np.nanmin(D, axis=1).tolist(),
            "depth_mean_m": np.nanmean(D, axis=1).tolist(),
            "depth_max_m": np.nanmax(D, axis=1).tolist(),
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

def _read_stations(path: str, unit_scale: float):
    """Stations CSV with a header: x,y,z,nx,ny,nz or r,z,n_r,n_z (optional s)."""
    import csv

    with open(path, encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError(f"{path} has no stations")
    keys = {k.strip().lower().split("_mm")[0].split("[")[0]: k
            for k in rows[0]}

    def col(*names):
        for n in names:
            if n in keys:
                return np.asarray([float(r[keys[n]]) for r in rows])
        return None

    s = col("s")
    if col("x") is not None:
        o = np.c_[col("x"), col("y"), col("z")] * unit_scale
        n = np.c_[col("nx"), col("ny"), col("nz")]
    elif col("r") is not None:
        a = col("z", "a")
        o = np.c_[col("r"), a] * unit_scale
        n = np.c_[col("n_r", "nr"), col("n_z", "n_a", "nz")]
    else:
        raise ValueError(f"{path}: need columns x,y,z,nx,ny,nz or "
                         "r,z,n_r,n_z")
    if np.any(np.isnan(n)):
        raise ValueError(f"{path}: missing normal columns")
    return o, n, s


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
                    help="mesh .vol (default: from the sidecar)")
    ap.add_argument("--order", type=int, default=None,
                    help="H1 order (default: from the sidecar)")
    ap.add_argument("--axisymmetric", action="store_true",
                    help="the mesh is an (r, z) axisymmetric section")
    ap.add_argument("--axis", default="z", choices=("x", "y", "z"))
    ap.add_argument("--thresholds", default="",
                    help="exposure: comma-separated temperatures [degC]")
    ap.add_argument("--threshold", type=float, default=None,
                    help="depth: threshold temperature [degC]")
    ap.add_argument("--stations", default="",
                    help="depth: stations CSV (x,y,z,nx,ny,nz or "
                         "r,z,n_r,n_z; meridian stations are marched on "
                         "--n-phi azimuths)")
    ap.add_argument("--stations-unit", default="m", choices=("m", "mm"))
    ap.add_argument("--boundaries", default="",
                    help="depth: boundary names for automatic stations")
    ap.add_argument("--span", type=float, default=None,
                    help="depth: probe length [m]")
    ap.add_argument("--step", type=float, default=None,
                    help="depth: probe step [m]")
    ap.add_argument("--n-phi", type=int, default=0,
                    help="depth: azimuths per meridian station")
    ap.add_argument("--output", default="", help="JSON output file")
    ap.add_argument("--csv", default="", help="depth: CSV output file")
    a = ap.parse_args(argv)

    mesh, gf, audit = _it.load_field(
        a.temperature, a.mesh, fes_order=a.order,
        quantity=_it.TEMPERATURE_QUANTITY if _it.read_field_sidecar(
            a.temperature) else None)
    if a.command == "exposure":
        thr = [float(t) for t in a.thresholds.split(",") if t.strip()]
        result = thermal_exposure(mesh, gf, thr, axisymmetric=a.axisymmetric,
                                  axis=a.axis)
    else:
        if a.threshold is None:
            ap.error("depth needs --threshold")
        s_col = None
        if a.stations:
            o, n, s_col = _read_stations(
                a.stations, 1e-3 if a.stations_unit == "mm" else 1.0)
            result = case_depth(mesh, gf, a.threshold, origins=o, normals=n,
                                span=a.span, step=a.step, axis=a.axis,
                                n_phi=a.n_phi)
        else:
            names = [b for b in a.boundaries.split(",") if b.strip()]
            result = case_depth(mesh, gf, a.threshold, boundary_names=names,
                                span=a.span, step=a.step)
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