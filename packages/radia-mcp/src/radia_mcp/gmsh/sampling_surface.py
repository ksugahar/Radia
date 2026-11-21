"""Resample a planar field-animation record with native Gmsh 5.

A ``radia.field-animation.v1`` record holds a signed vector field at fixed
sample points per solver step.  When those samples are the nodes of a
rectilinear grid on one plane (the Radia rotating-magnet demo samples the
plate mid-plane this way), the grid is a declared 2-D *sampling surface*:
its Quad4 connectivity is implied by the sample ids ``id = i * n2 + k``.
This module

1. checks that association exactly (planarity, grid order, ids, finite data,
   increasing time) and refuses anything else;
2. writes the surface and one NodeData time step per solver step to ASCII
   MSH 4.2 with Radia's ``gmsh_post_export`` writer;
3. runs the qualified Gmsh 5 runtime (``qualify_runtime``) to probe signed
   values at explicit query points, step by step, with ``distanceMax=0`` and
   ``dim=2`` so only exact in-element matches count;
4. writes a separate field-animation record with the original time, solver
   steps, magnet poses and region, plus the provenance of every input.

Limits (declared in the output): the result is the bilinear Quad4
interpolant of the *sampled* data on that plane within one solver step.  It is
not the solver's finite-element field, nothing is interpolated in time, and
points off the plane or outside the grid are rejected rather than projected
or extrapolated.
"""

from __future__ import annotations

import hashlib
import json
import os
import tempfile
from pathlib import Path

import numpy as np

SCHEMA = "radia.field-animation.v1"
FIELD_KEY = "J_A_per_m2"
ASSOCIATION = ("NodeData on the declared planar rectilinear sampling surface "
               "(nodes = source samples, id = i*n2 + k, Quad4 cells); Gmsh "
               "Quad4 interpolation inside one solver step; no time interpolation")
LIMITS = ("2-D sampling surface only: values are the bilinear interpolant of "
          "the sampled solver output on this plane, not the solver FE field; "
          "off-plane, outside-grid and boundary points are rejected")


def _sha256(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def grid_association(points, ids, rel_tol=1e-9):
    """Return (normal_axis, axes, coords, shape) of a planar sample grid.

    ``points[j]`` must equal ``(u_i, w_k)`` on the two in-plane axes with
    ``ids[j] = i * n2 + k`` and a constant coordinate on the normal axis.
    Raises ValueError for anything else (scattered points, permuted ids,
    non-planar or non-rectilinear samples).
    """
    pts = np.asarray(points, float)
    ids = np.asarray(ids)
    if pts.ndim != 2 or pts.shape[1] != 3 or len(ids) != len(pts):
        raise ValueError(f"samples must be (m, 3) with m ids, got {pts.shape} and {ids.shape}")
    if not np.all(np.isfinite(pts)):
        raise ValueError("sample points contain non-finite values")
    if not np.array_equal(np.sort(ids), np.arange(len(ids))):
        raise ValueError("sample ids must be 0..m-1 for the declared grid association")
    span = pts.max(0) - pts.min(0)
    scale = float(span.max())
    if scale <= 0:
        raise ValueError("sample points are coincident")
    flat = np.flatnonzero(span <= rel_tol * scale)
    if len(flat) != 1:
        raise ValueError(f"samples are not on one coordinate plane (extent {span.tolist()})")
    normal = int(flat[0])
    axes = [a for a in range(3) if a != normal]
    ordered = pts[np.argsort(ids)]
    c1 = np.unique(np.round(ordered[:, axes[0]] / (rel_tol * scale)))
    c2 = np.unique(np.round(ordered[:, axes[1]] / (rel_tol * scale)))
    n1, n2 = len(c1), len(c2)
    if n1 < 2 or n2 < 2 or n1 * n2 != len(ids):
        raise ValueError(f"samples are not a rectilinear grid ({n1} x {n2} != {len(ids)})")
    u = ordered[::n2, axes[0]]
    w = ordered[:n2, axes[1]]
    expected = np.empty_like(ordered)
    expected[:, normal] = ordered[0, normal]
    expected[:, axes[0]] = np.repeat(u, n2)
    expected[:, axes[1]] = np.tile(w, n1)
    if np.abs(ordered - expected).max() > rel_tol * scale:
        raise ValueError("sample ids do not follow id = i * n2 + k on the grid")
    if not (np.all(np.diff(u) > 0) and np.all(np.diff(w) > 0)):
        raise ValueError("grid coordinates must increase with i and k")
    return normal, axes, (u, w, float(ordered[0, normal])), (n1, n2)


def write_surface_msh(path, points, ids, time_s, field, steps):
    """Write the sampling surface and the selected steps as MSH 4.2 NodeData."""
    from radia.gmsh_post_export import (DataField, ElementBlock, Entity, NodeBlock,
                                        PayloadV42, PhysicalGroup, _emit_v42_msh)
    normal, axes, (u, w, c0), (n1, n2) = grid_association(points, ids)
    pts = np.asarray(points, float)[np.argsort(ids)]
    lo, hi = pts.min(0), pts.max(0)
    node = lambda i, k: i * n2 + k + 1          # gmsh node id = sample id + 1
    quads = [[node(i, k), node(i + 1, k), node(i + 1, k + 1), node(i, k + 1)]
             for i in range(n1 - 1) for k in range(n2 - 1)]
    payload = PayloadV42(
        groups=[PhysicalGroup(2, 1, "sampling_surface")],
        entities=[Entity(2, 1, (*lo, *hi), physical_tags=[1])],
        node_blocks=[NodeBlock(2, 1, 1, [tuple(p) for p in pts])],
        element_blocks=[ElementBlock(2, 1, 3, quads)])
    field = np.asarray(field, float)
    order = np.argsort(ids)
    for s, n in enumerate(steps):      # Gmsh step index = position in ``steps``
        payload.data_fields.append(DataField(
            "$NodeData", FIELD_KEY, 3,
            [(j + 1, field[n, order[j]]) for j in range(len(pts))],
            time=float(time_s[n]), timestep=s))
    with open(path, "w", encoding="ascii", newline="\n") as f:
        _emit_v42_msh(f, payload)
    return {"normal_axis": normal, "in_plane_axes": axes, "plane_coordinate_m": c0,
            "grid_shape": [n1, n2], "n_quad4": len(quads)}


_PROBE_SCRIPT = r'''
import json, sys
import gmsh
msh, req_path, out = sys.argv[1], sys.argv[2], sys.argv[-1]
req = json.load(open(req_path, encoding="utf-8"))
gmsh.initialize(["-noconfig"], readConfigFiles=False)
try:
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.open(msh)
    views = gmsh.view.getTags()
    if len(views) != 1:
        raise RuntimeError(f"expected one J view, found {len(views)}")
    view = views[0]
    times, values, distances = [], [], []
    for s, expected_step in enumerate(req["steps"]):
        kind, tags, data, time, ncomp = gmsh.view.getModelData(view, s)
        if kind != "NodeData" or ncomp != 3:
            raise RuntimeError(f"step {s}: {kind} with {ncomp} components")
        times.append(time)
        row, drow = [], []
        for x, y, z in req["points"]:
            v, d = gmsh.view.probe(view, x, y, z, step=s, numComp=3,
                                   distanceMax=0.0, dim=2)
            row.append(list(v)); drow.append(d)
        values.append(row); distances.append(drow)
    json.dump({"ok": True, "ran": True, "times": times, "values": values,
               "distances": distances}, open(out, "w", encoding="utf-8"))
except Exception as exc:
    json.dump({"ok": False, "ran": True, "error": f"{type(exc).__name__}: {exc}"},
              open(out, "w", encoding="utf-8"))
finally:
    gmsh.finalize()
'''


def _check_output_target(source_h5, out_h5):
    """Refuse an output that is, or could replace, the source or another file."""
    src, out = Path(source_h5), Path(out_h5)
    if not src.is_file():
        raise FileNotFoundError(f"source record {src} does not exist")
    if out.exists() or out.is_symlink():
        if out.exists() and os.path.samefile(src, out):
            raise ValueError(f"output {out} is the source record (same file); refusing to overwrite it")
        raise FileExistsError(f"output {out} already exists; choose a new path")
    if out.resolve() == src.resolve():
        raise ValueError(f"output {out} resolves to the source record {src}")
    if not out.parent.is_dir():
        raise FileNotFoundError(f"output directory {out.parent} does not exist")


def _publish_exclusive(tmp, out):
    """Make ``tmp`` visible as ``out`` only if ``out`` does not exist yet."""
    try:
        os.link(tmp, out)          # fails with FileExistsError if out appeared meanwhile
    finally:
        os.unlink(tmp)


def check_query_points(points, ids, query, rel_tol=1e-9):
    """Raise unless every query point is strictly inside the planar grid."""
    normal, axes, (u, w, c0), _ = grid_association(points, ids)
    q = np.asarray(query, float)
    if q.size == 0:
        raise ValueError("no query points selected")
    if q.ndim != 2 or q.shape[1] != 3 or not np.all(np.isfinite(q)):
        raise ValueError(f"query points must be finite (q, 3), got {q.shape}")
    scale = max(u[-1] - u[0], w[-1] - w[0])
    off = np.abs(q[:, normal] - c0) > rel_tol * scale
    if off.any():
        raise ValueError(f"{int(off.sum())} query points are off the sampling plane "
                         f"(first index {int(np.flatnonzero(off)[0])})")
    inside = ((q[:, axes[0]] > u[0]) & (q[:, axes[0]] < u[-1])
              & (q[:, axes[1]] > w[0]) & (q[:, axes[1]] < w[-1]))
    if not inside.all():
        raise ValueError(f"{int((~inside).sum())} query points are outside or on the "
                         f"grid boundary (first index {int(np.flatnonzero(~inside)[0])})")


def resample(source_h5, out_h5, query_points, query_ids, *, steps=None,
             workdir=None, timeout_s=600.0):
    """Resample /field/J_A_per_m2 of ``source_h5`` at query points via Gmsh 5.

    ``source_h5`` is only read; its SHA-256 is taken before and after reading
    and again before publishing.  ``out_h5`` must be a new path: the record is
    written to a temporary file in the same directory and published only if
    ``out_h5`` still does not exist.  ``workdir``, if given, must be an
    existing directory; scratch files go to a new unique subdirectory of it
    that is removed afterwards.  Returns a summary dict; raises on any
    target, association, time, field, runtime or probe failure.
    """
    import h5py
    from radia_mcp.gmsh._gmsh_subprocess import run_gmsh_json_subprocess
    _check_output_target(source_h5, out_h5)
    source_sha = _sha256(source_h5)
    with h5py.File(source_h5, "r") as f:
        schema = f.attrs.get("schema")
        schema = schema.decode() if isinstance(schema, bytes) else schema
        if schema != SCHEMA:
            raise ValueError(f"source schema {schema!r} is not {SCHEMA!r}")
        if f"field/{FIELD_KEY}" not in f:
            raise ValueError(f"source has no /field/{FIELD_KEY}")
        t = f["time_s"][()]
        J = f[f"field/{FIELD_KEY}"][()]
        ids = f["samples/id"][()]
        pts = f["samples/point_m"][()]
        if J.shape != (len(t), len(ids), 3):
            raise ValueError(f"J shape {J.shape} != ({len(t)}, {len(ids)}, 3)")
        if not (np.all(np.isfinite(t)) and np.all(np.diff(t) > 0)):
            raise ValueError("time_s must be finite and strictly increasing")
        if not np.all(np.isfinite(J)):
            raise ValueError("J contains non-finite values")
        copy = {k: f[k][()] for k in ("solver_step", "magnet/center_m", "magnet/rotation",
                                      "magnet/size_m", "magnet/magnetization_A_per_m",
                                      "conductor/box_min_m", "conductor/box_max_m") if k in f}
        conductor_attrs = dict(f["conductor"].attrs) if "conductor" in f else {}
        src_attrs = {k: (v.decode() if isinstance(v, bytes) else v) for k, v in f.attrs.items()}
    if _sha256(source_h5) != source_sha:
        raise RuntimeError("source record changed while it was being read")
    steps = list(range(len(t))) if steps is None else [int(s) for s in steps]
    if not steps:
        raise ValueError("no solver steps selected")
    if sorted(set(steps)) != steps or steps[0] < 0 or steps[-1] >= len(t) or len(steps) < 2:
        raise ValueError(f"steps must be >= 2 increasing indices into time_s, got {steps}")
    q = np.asarray(query_points, float)
    qid = np.asarray(query_ids, dtype=np.int64)
    if q.size == 0 or qid.size == 0:
        raise ValueError("no query points selected")
    if len(qid) != len(q) or len(np.unique(qid)) != len(qid):
        raise ValueError("query ids must be unique, one per query point")
    check_query_points(pts, ids, q)
    if workdir is not None and not Path(workdir).is_dir():
        raise NotADirectoryError(f"workdir {workdir!s} is not an existing directory")
    # Scratch always lives in a fresh, job-owned directory (under workdir when
    # given), so existing files there, including the source, are never opened.
    with tempfile.TemporaryDirectory(prefix="radia_surface_", dir=workdir) as tmp:
        work = Path(tmp)
        msh = work / "sampling_surface.msh"
        grid = write_surface_msh(msh, pts, ids, t, J, steps)
        req = work / "probe_request.json"
        req.write_text(json.dumps({"steps": steps, "points": q.tolist()}), encoding="utf-8")
        res = run_gmsh_json_subprocess(_PROBE_SCRIPT, [str(msh), str(req)],
                                       timeout_s=timeout_s, prefix="radia_surface_probe_")
        if not res.get("ok"):
            raise RuntimeError(f"Gmsh probe failed: {res.get('error')}")
        msh_sha = _sha256(msh)
    V = _native_values(res, len(steps), len(q))
    times = np.asarray(res["times"], float)
    if not np.allclose(times, t[steps], rtol=1e-14, atol=0):
        raise RuntimeError(f"Gmsh step times {times.tolist()} differ from source {t[steps].tolist()}")
    dist = np.asarray(res["distances"], float)
    if np.any(dist != 0.0):
        raise RuntimeError(f"probe returned non-exact matches (max distance {dist.max():.3e})")
    runtime = res["runtime"]
    out_h5 = Path(out_h5)
    fd, tmp = tempfile.mkstemp(prefix=f".{out_h5.name}.", suffix=".tmp", dir=out_h5.parent)
    os.close(fd)
    try:
        _write_output(h5py, tmp, src_attrs, source_h5, source_sha, msh_sha, grid, runtime,
                      t, steps, copy, conductor_attrs, qid, q, V)
        if _sha256(source_h5) != source_sha:
            raise RuntimeError("source record changed during resampling")
        _publish_exclusive(tmp, out_h5)
    except BaseException:
        if os.path.exists(tmp):
            os.unlink(tmp)
        raise
    return {"out_h5": str(out_h5), "out_sha256": _sha256(out_h5), "steps": steps,
            "n_query": len(q), "grid": grid, "surface_msh_sha256": msh_sha,
            "source_h5_sha256": source_sha, "gmsh_runtime": runtime}


def _native_values(res, n_steps, n_points):
    """Return the probe values as a finite (steps, points, 3) array or raise."""
    values, times, dist = res.get("values"), res.get("times"), res.get("distances")
    expected = (n_steps, n_points, 3)
    if not (isinstance(values, list) and isinstance(times, list) and isinstance(dist, list)
            and len(values) == len(times) == len(dist) == n_steps
            and all(isinstance(r, list) and len(r) == n_points for r in values)
            and all(isinstance(r, list) and len(r) == n_points for r in dist)):
        raise RuntimeError(f"Gmsh values do not have shape {expected}")
    bad = [(s, j) for s, row in enumerate(values) for j, v in enumerate(row)
           if not (isinstance(v, list) and len(v) == 3)]
    if bad:
        raise RuntimeError(f"{len(bad)} probes found no element or not 3 components "
                           f"(first step/point {bad[0]}); expected shape {expected}")
    V = np.asarray(values, float)
    if not np.all(np.isfinite(V)):
        raise RuntimeError("Gmsh values contain non-finite entries")
    return V


def _write_output(h5py, path, src_attrs, source_h5, source_sha, msh_sha, grid, runtime,
                  t, steps, copy, conductor_attrs, qid, q, V):
    with h5py.File(path, "w") as g:
        for k, v in src_attrs.items():
            g.attrs[f"source_{k}"] = v
        g.attrs.update({
            "schema": SCHEMA, "units": "SI", "index_base": 0,
            "case": src_attrs.get("case", ""),
            "reference_status": src_attrs.get("reference_status", ""),
            "derivation": "Gmsh 5 resampling of the source record on its declared sampling surface",
            "interpolation_association": ASSOCIATION, "limits": LIMITS,
            "source_h5": Path(source_h5).name, "source_h5_sha256": source_sha,
            "surface_msh_sha256": msh_sha, "surface_msh_format": "ASCII MSH 4.2",
            "surface_grid": json.dumps(grid),
            "probe": "gmsh.view.probe(step=s, numComp=3, distanceMax=0, dim=2)",
            "gmsh_runtime": json.dumps(runtime, sort_keys=True),
            "gmsh_wheel_sha256": runtime.get("wheel_sha256", ""),
            "gmsh_version": runtime.get("gmsh_version", ""),
        })
        g["time_s"] = t[steps]
        g["source_step_index"] = np.asarray(steps, np.int64)
        for k, v in copy.items():
            if k == "solver_step":
                g[k] = v[steps]
            elif k.startswith("magnet/") and k.split("/")[1] in ("center_m", "rotation"):
                g[k] = v[steps]
            else:
                g[k] = v
        for k, v in conductor_attrs.items():
            g["conductor"].attrs[k] = v
        g["samples/id"] = qid
        g["samples/point_m"] = q
        g["field/J_A_per_m2"] = V
