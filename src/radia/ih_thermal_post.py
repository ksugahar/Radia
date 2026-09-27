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
