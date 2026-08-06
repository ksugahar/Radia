"""Field post-processing tools: probe, integrate, transform, extract.

Wraps gmsh's ACTUAL post-processing machinery (gmsh.view.probe and the
gmsh.plugin.* family) as one-shot subprocess tools -- the missing
"post-processing verbs" next to the inspection/rendering lanes.

Empirically verified semantics (gmsh 4.15.2, locked by tests):

- ``gmsh.view.probe`` returns ``(values, distance)``; with ``step=-1``
  the values of ALL time steps are concatenated.  Outside the mesh the
  values are empty and ``distance`` is the gap to the nearest element.
- Plugin(Integrate) sums over elements of EVERY dimension by default
  (volume integral + surface integral + ...); pass ``dimension`` to
  select one.  Its result is a list-based SP view read via
  ``getListData`` (x, y, z, value_step0, value_step1, ...).
- Plugin(MathEval) applies the expression to NODAL values; the result
  view interpolates f(node values) -- interp(T^2) at a tet center is
  the average of squared vertex values, not (interp T)^2.
"""

from __future__ import annotations

import json
import tempfile
from pathlib import Path
from typing import Any

from ._gmsh_subprocess import run_gmsh_json_subprocess

_POST_SCRIPT = r"""
import json
import sys

cfg_path, out_path = sys.argv[1], sys.argv[2]
with open(cfg_path, encoding="utf-8") as f:
    cfg = json.load(f)
result = {"ok": False, "ran": False}


def _view_tag_by_selector(tags, selector):
    # selector: None -> index 0; int -> view index; str -> view name.
    import gmsh
    if selector is None:
        return tags[0]
    if isinstance(selector, int):
        if not 0 <= selector < len(tags):
            raise RuntimeError(
                f"view index {selector} out of range (0..{len(tags)-1})")
        return tags[selector]
    names = {}
    for tag in tags:
        idx = gmsh.view.getIndex(tag)
        names[gmsh.option.getString(f"View[{idx}].Name")] = tag
    if selector not in names:
        raise RuntimeError(
            f"view {selector!r} not found; available: {sorted(names)}")
    return names[selector]


def _n_steps(tag):
    import gmsh
    idx = gmsh.view.getIndex(tag)
    return int(gmsh.option.getNumber(f"View[{idx}].NbTimeStep"))


def _view_name(tag):
    import gmsh
    idx = gmsh.view.getIndex(tag)
    return gmsh.option.getString(f"View[{idx}].Name")


def _view_ncomp(tag):
    # Storage-kind dispatch: list-based views encode the component count
    # in the data-type letter (S/V/T); model-based views report it via
    # getHomogeneousModelData.
    import gmsh
    try:
        dtypes, _ne, _data = gmsh.view.getListData(tag)
    except Exception:
        dtypes = ()
    if dtypes:
        return {"S": 1, "V": 3, "T": 9}[str(dtypes[0])[0]]
    _dt, _tags, _data, _time, ncomp = gmsh.view.getHomogeneousModelData(
        tag, 0)
    return int(ncomp)


def _materialize(tag, exprs=None, name=None):
    # Identity (or sign-flipped) MathEval -> NEW list-based copy, leaving
    # the input view untouched.  Needed because Transform / Warp / Smooth
    # / ModulusPhase modify their input IN PLACE (measured gmsh 4.15.2).
    import gmsh
    if exprs is None:
        exprs = [f"v{i}" for i in range(_view_ncomp(tag))]
    gmsh.plugin.setNumber("MathEval", "View", gmsh.view.getIndex(tag))
    gmsh.plugin.setNumber("MathEval", "OtherView", -1)
    gmsh.plugin.setNumber("MathEval", "TimeStep", -1)
    for i in range(9):
        gmsh.plugin.setString("MathEval", f"Expression{i}",
                              exprs[i] if i < len(exprs) else "")
    out = gmsh.plugin.run("MathEval")
    if name:
        gmsh.view.option.setString(out, "Name", name)
    return out


def _parse_point_rows(tag):
    # Parse an SP/VP list view into points + per-step samples.
    import gmsh
    dtypes, nels, data = gmsh.view.getListData(tag)
    n_steps = max(1, _n_steps(tag))
    points, values = [], []
    for dtype, ne, block in zip(dtypes, nels, data):
        kind = str(dtype)[0]
        ncomp = {"S": 1, "V": 3, "T": 9}[kind]
        width = 3 + ncomp * n_steps
        ne = int(ne)
        for e in range(ne):
            row = [float(v) for v in block[e * width:(e + 1) * width]]
            points.append(row[:3])
            vals = row[3:]
            values.append([vals[s * ncomp:(s + 1) * ncomp]
                           for s in range(n_steps)])
    return points, values


try:
    import gmsh
    gmsh.initialize(["-noconfig"])
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.open(cfg["path"])
        tags = list(gmsh.view.getTags())
        op = cfg["op"]

        if op == "probe":
            sel = cfg.get("view")
            targets = tags if sel is None else [_view_tag_by_selector(tags, sel)]
            if not targets:
                raise RuntimeError("no post-processing views in the file")
            step = int(cfg.get("step", -1))
            dmax = float(cfg.get("distance_max", 0.0))
            views_out = []
            for tag in targets:
                n_steps = _n_steps(tag)
                per_point = []
                for pt in cfg["points"]:
                    values, distance = gmsh.view.probe(
                        tag, pt[0], pt[1], pt[2], step=step,
                        distanceMax=dmax)
                    values = [float(v) for v in values]
                    # List-based VECTOR views return the NEAREST value at
                    # any distance (measured; scalar views return empty
                    # outside) -- gate "found" on the distance, not just
                    # on values being present.
                    entry = {"point": pt,
                             "found": bool(values)
                             and float(distance) <= max(dmax, 0.0),
                             "distance": float(distance)}
                    if values:
                        if step == -1 and n_steps > 1:
                            ncomp = len(values) // n_steps
                            entry["steps"] = [
                                values[s * ncomp:(s + 1) * ncomp]
                                for s in range(n_steps)]
                        else:
                            entry["values"] = values
                    per_point.append(entry)
                views_out.append({"name": _view_name(tag),
                                  "n_steps": n_steps,
                                  "points": per_point})
            result.update({"ok": True, "ran": True, "views": views_out})

        elif op == "integrate":
            tag = _view_tag_by_selector(tags, cfg.get("view"))
            gmsh.plugin.setNumber("Integrate", "View",
                                  gmsh.view.getIndex(tag))
            gmsh.plugin.setNumber("Integrate", "Dimension",
                                  int(cfg.get("dimension", -1)))
            out_tag = gmsh.plugin.run("Integrate")
            _dt, _ne, data = gmsh.view.getListData(out_tag)
            values = [float(v) for v in data[0][3:]] if data else []
            result.update({"ok": True, "ran": True,
                           "view": _view_name(tag),
                           "dimension": int(cfg.get("dimension", -1)),
                           "integral_per_step": values})

        elif op == "math_eval":
            tag = _view_tag_by_selector(tags, cfg.get("view"))
            gmsh.plugin.setNumber("MathEval", "View",
                                  gmsh.view.getIndex(tag))
            other = cfg.get("other_view")
            if other is not None:
                other_tag = _view_tag_by_selector(tags, other)
                gmsh.plugin.setNumber("MathEval", "OtherView",
                                      gmsh.view.getIndex(other_tag))
            gmsh.plugin.setNumber("MathEval", "TimeStep",
                                  int(cfg.get("time_step", -1)))
            expressions = cfg["expressions"]
            for i in range(9):
                expr = expressions[i] if i < len(expressions) else ""
                gmsh.plugin.setString("MathEval", f"Expression{i}", expr)
            out_tag = gmsh.plugin.run("MathEval")
            gmsh.view.option.setString(out_tag, "Name", cfg.get(
                "result_name", "math_eval"))
            gmsh.view.write(out_tag, cfg["out_file"])
            result.update({"ok": True, "ran": True,
                           "out_file": cfg["out_file"],
                           "n_steps": _n_steps(out_tag)})

        elif op == "isosurface":
            tag = _view_tag_by_selector(tags, cfg.get("view"))
            gmsh.plugin.setNumber("Isosurface", "View",
                                  gmsh.view.getIndex(tag))
            gmsh.plugin.setNumber("Isosurface", "Value",
                                  float(cfg["value"]))
            out_tag = gmsh.plugin.run("Isosurface")
            dtypes, nels, _data = gmsh.view.getListData(out_tag)
            gmsh.view.write(out_tag, cfg["out_file"])
            result.update({"ok": True, "ran": True,
                           "out_file": cfg["out_file"],
                           "pieces": {str(t): int(n)
                                      for t, n in zip(dtypes, nels)}})

        elif op == "cut_plane":
            tag = _view_tag_by_selector(tags, cfg.get("view"))
            a, b, c = cfg["normal"]
            gmsh.plugin.setNumber("CutPlane", "View",
                                  gmsh.view.getIndex(tag))
            gmsh.plugin.setNumber("CutPlane", "A", float(a))
            gmsh.plugin.setNumber("CutPlane", "B", float(b))
            gmsh.plugin.setNumber("CutPlane", "C", float(c))
            gmsh.plugin.setNumber("CutPlane", "D", float(cfg["offset"]))
            out_tag = gmsh.plugin.run("CutPlane")
            dtypes, nels, _ = gmsh.view.getListData(out_tag)
            gmsh.view.write(out_tag, cfg["out_file"])
            result.update({"ok": True, "ran": True,
                           "out_file": cfg["out_file"],
                           "pieces": {str(t): int(n)
                                      for t, n in zip(dtypes, nels)}})

        elif op == "harmonic_to_time":
            tag = _view_tag_by_selector(tags, cfg.get("view"))
            gmsh.plugin.setNumber("HarmonicToTime", "View",
                                  gmsh.view.getIndex(tag))
            gmsh.plugin.setNumber("HarmonicToTime", "RealPart",
                                  int(cfg.get("real_step", 0)))
            gmsh.plugin.setNumber("HarmonicToTime", "ImaginaryPart",
                                  int(cfg.get("imag_step", 1)))
            gmsh.plugin.setNumber("HarmonicToTime", "NumSteps",
                                  int(cfg.get("n_steps", 20)))
            out_tag = gmsh.plugin.run("HarmonicToTime")
            gmsh.view.write(out_tag, cfg["out_file"])
            result.update({"ok": True, "ran": True,
                           "out_file": cfg["out_file"],
                           "n_steps": _n_steps(out_tag)})

        elif op == "streamlines":
            # Probe-driven arc-length RK4 (Plugin(StreamLines) on this
            # gmsh build only re-emits the seed points): unit-tangent
            # integration, both directions, |v| carried as line color.
            tag = _view_tag_by_selector(tags, cfg.get("view"))
            step = int(cfg.get("time_step", 0))
            ds = cfg.get("step_size")
            if ds is None:
                bb = gmsh.model.getBoundingBox(-1, -1)
                diag = sum((bb[i + 3] - bb[i]) ** 2 for i in range(3)) ** 0.5
                ds = diag / 200.0
            ds = float(ds)
            max_steps = int(cfg.get("max_steps", 400))
            n_seeds = max(1, int(cfg.get("n_seeds", 10)))
            s0, s1 = cfg["seed_start"], cfg["seed_end"]
            seeds = [[s0[k] + (s1[k] - s0[k]) *
                      (i / (n_seeds - 1) if n_seeds > 1 else 0.0)
                      for k in range(3)] for i in range(n_seeds)]
            directions = ([1, -1] if cfg.get("both_directions", True)
                          else [1])

            def _sample(q):
                values, dist = gmsh.view.probe(tag, q[0], q[1], q[2],
                                               step=step)
                # dist > 0 = outside the data (list-based vector views
                # report the nearest value at ANY distance -- without
                # this gate the tracer would never stop).
                if len(values) < 3 or float(dist) > 0.0:
                    return None, 0.0
                vx, vy, vz = (float(values[0]), float(values[1]),
                              float(values[2]))
                norm = (vx * vx + vy * vy + vz * vz) ** 0.5
                if norm < 1e-30:
                    return None, 0.0
                return (vx / norm, vy / norm, vz / norm), norm

            polylines = []
            for seed in seeds:
                for sgn in directions:
                    path = [list(seed)]
                    norms = [_sample(seed)[1]]
                    p = list(seed)
                    for _ in range(max_steps):
                        def _f(q):
                            t, _n = _sample(q)
                            if t is None:
                                return None
                            return [sgn * c for c in t]
                        k1 = _f(p)
                        if k1 is None:
                            break
                        q2 = [p[i] + 0.5 * ds * k1[i] for i in range(3)]
                        k2 = _f(q2)
                        if k2 is None:
                            break
                        q3 = [p[i] + 0.5 * ds * k2[i] for i in range(3)]
                        k3 = _f(q3)
                        if k3 is None:
                            break
                        q4 = [p[i] + ds * k3[i] for i in range(3)]
                        k4 = _f(q4)
                        if k4 is None:
                            break
                        p = [p[i] + ds / 6.0 * (k1[i] + 2 * k2[i]
                                                + 2 * k3[i] + k4[i])
                             for i in range(3)]
                        _t, norm = _sample(p)
                        if _t is None:
                            break
                        path.append(list(p))
                        norms.append(norm)
                    if len(path) > 1:
                        if sgn < 0:
                            path.reverse()
                            norms.reverse()
                        polylines.append((path, norms))

            sl_data = []
            n_segments = 0
            for path, norms in polylines:
                for i in range(len(path) - 1):
                    a, b = path[i], path[i + 1]
                    sl_data += [a[0], b[0], a[1], b[1], a[2], b[2],
                                norms[i], norms[i + 1]]
                    n_segments += 1
            out_tag = gmsh.view.add("streamlines")
            if n_segments:
                gmsh.view.addListData(out_tag, "SL", n_segments, sl_data)
            gmsh.view.write(out_tag, cfg["out_file"])
            result.update({
                "ok": True, "ran": True,
                "out_file": cfg["out_file"],
                "n_polylines": len(polylines),
                "n_segments": n_segments,
                "step_size": ds,
                "polylines": [
                    {"points": path, "field_norms": norms}
                    for path, norms in polylines
                ] if cfg.get("return_points") else None,
            })

        elif op == "derived":
            tag = _view_tag_by_selector(tags, cfg.get("view"))
            plugin = cfg["plugin"]
            gmsh.plugin.setNumber(plugin, "View", gmsh.view.getIndex(tag))
            before = set(gmsh.view.getTags())
            gmsh.plugin.run(plugin)
            new = [t for t in gmsh.view.getTags() if t not in before]
            if not new:
                raise RuntimeError(f"Plugin({plugin}) produced no view")
            names, pieces = [], {}
            for k, t in enumerate(new):
                gmsh.view.write(t, cfg["out_file"], append=(k > 0))
                names.append(_view_name(t))
                dt, ne, _d = gmsh.view.getListData(t)
                for typ, n in zip(dt, ne):
                    pieces[str(typ)] = pieces.get(str(typ), 0) + int(n)
            result.update({"ok": True, "ran": True,
                           "out_file": cfg["out_file"],
                           "views": names, "pieces": pieces,
                           "probe": None})
            pt = cfg.get("check_point")
            if pt:
                vals, _d = gmsh.view.probe(new[0], pt[0], pt[1], pt[2])
                result["probe"] = [float(v) for v in vals]

        elif op == "threshold":
            tag = _view_tag_by_selector(tags, cfg.get("view"))
            gmsh.plugin.setNumber("ExtractElements", "View",
                                  gmsh.view.getIndex(tag))
            gmsh.plugin.setNumber("ExtractElements", "MinVal",
                                  float(cfg["min_val"]))
            gmsh.plugin.setNumber("ExtractElements", "MaxVal",
                                  float(cfg["max_val"]))
            gmsh.plugin.setNumber("ExtractElements", "TimeStep",
                                  int(cfg.get("time_step", 0)))
            gmsh.plugin.setNumber("ExtractElements", "Dimension",
                                  int(cfg.get("dimension", -1)))
            out_tag = gmsh.plugin.run("ExtractElements")
            dt, ne, _d = gmsh.view.getListData(out_tag)
            gmsh.view.write(out_tag, cfg["out_file"])
            result.update({"ok": True, "ran": True,
                           "out_file": cfg["out_file"],
                           "pieces": {str(t): int(n)
                                      for t, n in zip(dt, ne)},
                           "n_kept": int(sum(int(n) for n in ne))})

        elif op == "skin":
            tag = _view_tag_by_selector(tags, cfg.get("view"))
            gmsh.plugin.setNumber("Skin", "View", gmsh.view.getIndex(tag))
            gmsh.plugin.setNumber("Skin", "FromMesh",
                                  int(cfg.get("from_mesh", 0)))
            out_tag = gmsh.plugin.run("Skin")
            dt, ne, _d = gmsh.view.getListData(out_tag)
            gmsh.view.write(out_tag, cfg["out_file"])
            result.update({"ok": True, "ran": True,
                           "out_file": cfg["out_file"],
                           "pieces": {str(t): int(n)
                                      for t, n in zip(dt, ne)}})

        elif op == "mirror_expand":
            tag = _view_tag_by_selector(tags, cfg.get("view"))
            ncomp = _view_ncomp(tag)
            parity = cfg.get("parity", "scalar")
            if ncomp == 1:
                parity = "scalar"
            elif ncomp != 3:
                raise RuntimeError(
                    "mirror_expand supports scalar and 3-component vector "
                    f"views (got {ncomp} components); transform tensors "
                    "with transform_affine + explicit value_expressions")
            planes = cfg["planes"]
            origin = cfg.get("origin", [0.0, 0.0, 0.0])
            name = cfg.get("result_name", "mirrored")
            axis_idx = {"x": 0, "y": 1, "z": 2}
            made = [_materialize(tag, name=name)]
            subsets = []
            for mask in range(1, 1 << len(planes)):
                subsets.append([planes[i] for i in range(len(planes))
                                if mask >> i & 1])
            for subset in subsets:
                m = [1.0, 1.0, 1.0]
                for ax in subset:
                    m[axis_idx[ax]] = -1.0
                det = m[0] * m[1] * m[2]
                if parity == "scalar":
                    signs = [1.0]
                elif parity == "vector":
                    signs = list(m)
                else:  # pseudovector (B, H): v' = det(M) * M v
                    signs = [det * c for c in m]
                exprs = [("-" if signs[i] < 0 else "") + f"v{i}"
                         for i in range(ncomp)]
                copy = _materialize(tag, exprs=exprs, name=name)
                gmsh.plugin.setNumber("Transform", "View",
                                      gmsh.view.getIndex(copy))
                for r in range(3):
                    for c in range(3):
                        gmsh.plugin.setNumber(
                            "Transform", f"A{r + 1}{c + 1}",
                            m[r] if r == c else 0.0)
                for ax_name, k in (("Tx", 0), ("Ty", 1), ("Tz", 2)):
                    off = 2.0 * origin[k] if m[k] < 0 else 0.0
                    gmsh.plugin.setNumber("Transform", ax_name, off)
                gmsh.plugin.setNumber("Transform", "SwapOrientation",
                                      1 if det < 0 else 0)
                gmsh.plugin.run("Transform")
                made.append(copy)
            gmsh.view.combine("elements", "name", True, True)
            # combine renames the merged view to "<name>_Combine"
            # (measured gmsh 4.15.2)
            merged = [t for t in gmsh.view.getTags()
                      if _view_name(t) in (name, name + "_Combine")]
            if len(merged) != 1:
                raise RuntimeError(
                    f"view combine produced {len(merged)} views named "
                    f"{name!r} (expected 1)")
            gmsh.view.option.setString(merged[0], "Name", name)
            gmsh.view.write(merged[0], cfg["out_file"])
            result.update({"ok": True, "ran": True,
                           "out_file": cfg["out_file"],
                           "n_copies": len(subsets),
                           "parity": parity,
                           "n_steps": _n_steps(merged[0])})

        elif op == "transform_affine":
            tag = _view_tag_by_selector(tags, cfg.get("view"))
            exprs = cfg.get("value_expressions")
            copy = _materialize(tag, exprs=exprs,
                                name=cfg.get("result_name", "transformed"))
            A = cfg["matrix"]
            T = cfg.get("translation", [0.0, 0.0, 0.0])
            gmsh.plugin.setNumber("Transform", "View",
                                  gmsh.view.getIndex(copy))
            for r in range(3):
                for c in range(3):
                    gmsh.plugin.setNumber("Transform", f"A{r + 1}{c + 1}",
                                          float(A[3 * r + c]))
            gmsh.plugin.setNumber("Transform", "Tx", float(T[0]))
            gmsh.plugin.setNumber("Transform", "Ty", float(T[1]))
            gmsh.plugin.setNumber("Transform", "Tz", float(T[2]))
            gmsh.plugin.setNumber("Transform", "SwapOrientation",
                                  int(cfg.get("swap_orientation", 0)))
            gmsh.plugin.run("Transform")
            gmsh.view.write(copy, cfg["out_file"])
            result.update({"ok": True, "ran": True,
                           "out_file": cfg["out_file"],
                           "n_steps": _n_steps(copy)})

        elif op == "warp":
            tag = _view_tag_by_selector(tags, cfg.get("view"))
            gmsh.plugin.setNumber("Warp", "View", gmsh.view.getIndex(tag))
            gmsh.plugin.setNumber("Warp", "Factor", float(cfg["factor"]))
            gmsh.plugin.setNumber("Warp", "TimeStep",
                                  int(cfg.get("time_step", 0)))
            gmsh.plugin.run("Warp")  # in place: MOVES THE MODEL NODES
            # view.write converts the (displaced) model view to list
            # format directly -- a MathEval copy here would re-read
            # geometry through a different path (measured).
            gmsh.view.write(tag, cfg["out_file"])
            result.update({"ok": True, "ran": True,
                           "out_file": cfg["out_file"],
                           "n_steps": _n_steps(tag)})

        elif op == "smooth":
            tag = _view_tag_by_selector(tags, cfg.get("view"))
            gmsh.plugin.setNumber("Smooth", "View", gmsh.view.getIndex(tag))
            gmsh.plugin.run("Smooth")  # in place: element -> nodal
            gmsh.view.write(tag, cfg["out_file"])
            result.update({"ok": True, "ran": True,
                           "out_file": cfg["out_file"],
                           "n_steps": _n_steps(tag)})

        elif op == "modulus_phase":
            tag = _view_tag_by_selector(tags, cfg.get("view"))
            gmsh.plugin.setNumber("ModulusPhase", "View",
                                  gmsh.view.getIndex(tag))
            gmsh.plugin.setNumber("ModulusPhase", "RealPart",
                                  int(cfg.get("real_step", 0)))
            gmsh.plugin.setNumber("ModulusPhase", "ImaginaryPart",
                                  int(cfg.get("imag_step", 1)))
            gmsh.plugin.run("ModulusPhase")  # in place: steps -> mod/phase
            gmsh.view.write(tag, cfg["out_file"])
            result.update({"ok": True, "ran": True,
                           "out_file": cfg["out_file"],
                           "n_steps": _n_steps(tag)})

        elif op == "min_max":
            tag = _view_tag_by_selector(tags, cfg.get("view"))
            gmsh.plugin.setNumber("MinMax", "View", gmsh.view.getIndex(tag))
            gmsh.plugin.setNumber("MinMax", "OverTime",
                                  int(cfg.get("over_time", 0)))
            gmsh.plugin.setNumber("MinMax", "Argument", 1)
            before = set(gmsh.view.getTags())
            gmsh.plugin.run("MinMax")
            new = [t for t in gmsh.view.getTags() if t not in before]
            found = {}
            for t in new:
                nm = _view_name(t)
                kind = ("min" if nm.endswith("_Min") else
                        "max" if nm.endswith("_Max") else nm)
                _dt, _ne, data = gmsh.view.getListData(t)
                rows = [[float(v) for v in d] for d in data]
                entry = {"raw": rows}
                if rows and len(rows[0]) >= 4:
                    entry["point"] = rows[0][:3]
                    entry["values"] = rows[0][3:]
                found[kind] = entry
            result.update({"ok": True, "ran": True, **found})

        elif op == "curve_profile":
            tag = _view_tag_by_selector(tags, cfg.get("view"))
            n = int(cfg["n"])
            u0, u1 = float(cfg["u_min"]), float(cfg["u_max"])
            for key, val in (("View", gmsh.view.getIndex(tag)),
                             ("MinU", u0), ("MaxU", u1),
                             ("NumPointsU", n), ("ConnectPoints", 0),
                             ("MinV", 0.0), ("MaxV", 0.0),
                             ("NumPointsV", 1)):
                gmsh.plugin.setNumber("CutParametric", key, float(val))
            for key in ("X", "Y", "Z"):
                gmsh.plugin.setString("CutParametric", key,
                                      str(cfg[key.lower() + "_expr"]))
            out_tag = gmsh.plugin.run("CutParametric")
            points, values = _parse_point_rows(out_tag)
            if cfg.get("out_file"):
                gmsh.plugin.setNumber("CutParametric", "ConnectPoints", 1)
                line_tag = gmsh.plugin.run("CutParametric")
                gmsh.view.write(line_tag, cfg["out_file"])
            result.update({"ok": True, "ran": True,
                           "out_file": cfg.get("out_file"),
                           "points": points, "values": values,
                           "n_steps": max(1, _n_steps(out_tag))})

        elif op == "resample_grid":
            tag = _view_tag_by_selector(tags, cfg.get("view"))
            o = cfg["origin"]
            corners = {"1": cfg["u_point"], "2": cfg["v_point"],
                       "3": cfg["w_point"]}
            gmsh.plugin.setNumber("CutBox", "View",
                                  gmsh.view.getIndex(tag))
            for ax, val in zip(("X0", "Y0", "Z0"), o):
                gmsh.plugin.setNumber("CutBox", ax, float(val))
            for suffix, pt in corners.items():
                for ax, val in zip(("X", "Y", "Z"), pt):
                    gmsh.plugin.setNumber("CutBox", ax + suffix, float(val))
            gmsh.plugin.setNumber("CutBox", "NumPointsU", int(cfg["nu"]))
            gmsh.plugin.setNumber("CutBox", "NumPointsV", int(cfg["nv"]))
            gmsh.plugin.setNumber("CutBox", "NumPointsW", int(cfg["nw"]))
            gmsh.plugin.setNumber("CutBox", "ConnectPoints", 0)
            gmsh.plugin.setNumber("CutBox", "Boundary", 0)
            out_tag = gmsh.plugin.run("CutBox")
            points, values = _parse_point_rows(out_tag)
            if cfg.get("out_file"):
                gmsh.view.write(out_tag, cfg["out_file"])
            result.update({"ok": True, "ran": True,
                           "out_file": cfg.get("out_file"),
                           "points": points, "values": values,
                           "order": "w fastest, then v, then u",
                           "n_steps": max(1, _n_steps(out_tag))})

        else:
            raise RuntimeError(f"unknown op {op!r}")

        if cfg.get("plot_png") and op == "probe":
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            distances = cfg["plot_distances"]
            fig, ax = plt.subplots(figsize=(7, 4.5))
            view0 = result["views"][0]
            n_steps = view0["n_steps"]
            series = {}
            for d, entry in zip(distances, view0["points"]):
                rows = entry.get("steps") or (
                    [entry["values"]] if entry.get("values") else [])
                for s, row in enumerate(rows):
                    if len(row) > 1:
                        mag = sum(v * v for v in row) ** 0.5
                        series.setdefault(f"step {s} |v|", []).append((d, mag))
                    elif row:
                        series.setdefault(f"step {s}", []).append((d, row[0]))
            for label, pts in series.items():
                xs = [p[0] for p in pts]
                ys = [p[1] for p in pts]
                ax.plot(xs, ys, marker="." if len(xs) < 60 else None,
                        label=label)
            ax.set_xlabel("distance along line")
            ax.set_ylabel(view0["name"])
            ax.grid(True, alpha=0.4)
            if len(series) > 1 and n_steps <= 8:
                ax.legend(fontsize=8)
            fig.tight_layout()
            fig.savefig(cfg["plot_png"], dpi=130)
            result["plot_png"] = cfg["plot_png"]
    finally:
        gmsh.finalize()
except Exception as exc:
    result["error"] = f"{type(exc).__name__}: {exc}"
with open(out_path, "w", encoding="utf-8") as f:
    json.dump(result, f)
"""


def _run_post(cfg: dict[str, Any], timeout_s: float) -> dict[str, Any]:
    with tempfile.TemporaryDirectory(prefix="radia_mcp_gmsh_post_") as work:
        cfg_path = Path(work) / "post.json"
        cfg_path.write_text(json.dumps(cfg), encoding="utf-8")
        result = run_gmsh_json_subprocess(
            _POST_SCRIPT, [str(cfg_path)],
            timeout_s=timeout_s, prefix="radia_mcp_gmsh_post_")
    result.setdefault("input", cfg.get("path"))
    return result


def _check_input(path: str | Path) -> dict[str, Any] | None:
    p = Path(path)
    if not p.is_file():
        return {"ok": False, "error": f"file not found: {p}"}
    return None


def _default_out(path: str | Path, out_file: str | Path | None,
                 suffix: str) -> str:
    if out_file is not None:
        out = Path(out_file)
    else:
        src = Path(path)
        out = src.with_name(f"{src.stem}_{suffix}.pos")
    out.parent.mkdir(parents=True, exist_ok=True)
    return str(out)


def probe_field(path: str | Path, points: list[list[float]], *,
                view: str | int | None = None, step: int = -1,
                distance_max: float = 0.0,
                timeout_s: float = 300.0) -> dict[str, Any]:
    """Probe post views at arbitrary points (interpolated values)."""
    err = _check_input(path)
    if err:
        return err
    pts = [[float(c) for c in p] for p in points]
    if any(len(p) != 3 for p in pts):
        return {"ok": False, "error": "every point needs exactly 3 coords"}
    return _run_post({"op": "probe", "path": str(path), "points": pts,
                      "view": view, "step": step,
                      "distance_max": float(distance_max)}, timeout_s)


def line_profile(path: str | Path, start: list[float], end: list[float],
                 n: int = 100, *, view: str | int | None = None,
                 step: int = -1, plot_png: str | Path | None = None,
                 timeout_s: float = 300.0) -> dict[str, Any]:
    """Sample a view along a straight line; optionally plot a PNG graph."""
    err = _check_input(path)
    if err:
        return err
    n = max(2, int(n))
    p0 = [float(c) for c in start]
    p1 = [float(c) for c in end]
    if len(p0) != 3 or len(p1) != 3:
        return {"ok": False, "error": "start/end need exactly 3 coords"}
    length = sum((a - b) ** 2 for a, b in zip(p0, p1)) ** 0.5
    points = [[p0[k] + (p1[k] - p0[k]) * i / (n - 1) for k in range(3)]
              for i in range(n)]
    distances = [length * i / (n - 1) for i in range(n)]
    cfg: dict[str, Any] = {"op": "probe", "path": str(path),
                           "points": points, "view": view, "step": step,
                           "distance_max": 0.0,
                           "plot_distances": distances}
    if plot_png is not None:
        out = Path(plot_png)
        out.parent.mkdir(parents=True, exist_ok=True)
        cfg["plot_png"] = str(out)
    result = _run_post(cfg, timeout_s)
    if result.get("ok"):
        result["distances"] = distances
        result["length"] = length
    return result


def integrate_view(path: str | Path, *, view: str | int | None = None,
                   dimension: int = -1,
                   timeout_s: float = 300.0) -> dict[str, Any]:
    """Integrate a view over its elements (per time step).

    dimension=-1 SUMS integrals over every element dimension present
    (volume + surface + line); pass 3 or 2 to select one -- almost
    always what a physical quantity needs.  Accuracy (measured): the
    plugin integrates at piecewise-linear accuracy even on high-order
    elements (nonlinear integrands carry O(h^2) error); exact FE
    integrals belong to NGSolve Integrate on the solver side.
    """
    err = _check_input(path)
    if err:
        return err
    return _run_post({"op": "integrate", "path": str(path), "view": view,
                      "dimension": int(dimension)}, timeout_s)


def math_eval(path: str | Path, expressions: list[str] | str, *,
              view: str | int | None = None,
              other_view: str | int | None = None,
              time_step: int = -1,
              result_name: str = "math_eval",
              out_file: str | Path | None = None,
              timeout_s: float = 300.0) -> dict[str, Any]:
    """Create a derived view with Plugin(MathEval) and save it.

    Expressions use v0..v8 (view components) and w0..w8 (other view);
    they are applied to NODAL values (the result interpolates
    f(node values), the standard FEM post semantics).
    """
    err = _check_input(path)
    if err:
        return err
    exprs = [expressions] if isinstance(expressions, str) else list(expressions)
    return _run_post({"op": "math_eval", "path": str(path), "view": view,
                      "other_view": other_view, "time_step": int(time_step),
                      "expressions": [str(e) for e in exprs],
                      "result_name": str(result_name),
                      "out_file": _default_out(path, out_file, "math")},
                     timeout_s)


def isosurface(path: str | Path, value: float, *,
               view: str | int | None = None,
               out_file: str | Path | None = None,
               timeout_s: float = 300.0) -> dict[str, Any]:
    """Extract the isosurface of a scalar view and save it."""
    err = _check_input(path)
    if err:
        return err
    return _run_post({"op": "isosurface", "path": str(path), "view": view,
                      "value": float(value),
                      "out_file": _default_out(path, out_file, "iso")},
                     timeout_s)


def cut_plane_extract(path: str | Path, normal: list[float], offset: float,
                      *, view: str | int | None = None,
                      out_file: str | Path | None = None,
                      timeout_s: float = 300.0) -> dict[str, Any]:
    """Cut a view with the plane A*x+B*y+C*z+D=0 and save the section."""
    err = _check_input(path)
    if err:
        return err
    n = [float(c) for c in normal]
    if len(n) != 3 or not any(n):
        return {"ok": False, "error": "normal must be 3 coords, nonzero"}
    return _run_post({"op": "cut_plane", "path": str(path), "view": view,
                      "normal": n, "offset": float(offset),
                      "out_file": _default_out(path, out_file, "cut")},
                     timeout_s)


def harmonic_to_time(path: str | Path, *, view: str | int | None = None,
                     real_step: int = 0, imag_step: int = 1,
                     n_steps: int = 20,
                     out_file: str | Path | None = None,
                     timeout_s: float = 300.0) -> dict[str, Any]:
    """Expand a complex (re/im two-step) view into a time animation.

    v(t) = re*cos(2*pi*k/n) - im*sin(2*pi*k/n): the standard AC-field
    time expansion. Feed the output to gmsh_export_animation.
    """
    err = _check_input(path)
    if err:
        return err
    return _run_post({"op": "harmonic_to_time", "path": str(path),
                      "view": view, "real_step": int(real_step),
                      "imag_step": int(imag_step),
                      "n_steps": int(n_steps),
                      "out_file": _default_out(path, out_file, "time")},
                     timeout_s)


def streamlines(path: str | Path, seed_start: list[float],
                seed_end: list[float], *, n_seeds: int = 10,
                view: str | int | None = None,
                step_size: float | None = None, max_steps: int = 400,
                both_directions: bool = True, time_step: int = 0,
                return_points: bool = False,
                out_file: str | Path | None = None,
                timeout_s: float = 600.0) -> dict[str, Any]:
    """Trace field lines of a vector view from seeds on a line segment.

    Probe-driven arc-length RK4 (unit-tangent integration): step_size
    defaults to bbox diagonal / 200, lines march in both directions
    until they leave the data or hit max_steps, and the local |v| is
    carried as the line color. The traced polylines are saved as an SL
    view; return_points=True additionally returns the coordinates.
    """
    err = _check_input(path)
    if err:
        return err
    s0 = [float(c) for c in seed_start]
    s1 = [float(c) for c in seed_end]
    if len(s0) != 3 or len(s1) != 3:
        return {"ok": False, "error": "seed_start/seed_end need 3 coords"}
    return _run_post({"op": "streamlines", "path": str(path), "view": view,
                      "seed_start": s0, "seed_end": s1,
                      "n_seeds": int(n_seeds),
                      "step_size": step_size,
                      "max_steps": int(max_steps),
                      "both_directions": bool(both_directions),
                      "time_step": int(time_step),
                      "return_points": bool(return_points),
                      "out_file": _default_out(path, out_file, "stream")},
                     timeout_s)


# ======================================================================
# ParaView-parity verbs (plugin semantics measured on gmsh 4.15.2)
# ======================================================================

_DERIVED_PLUGINS = {"gradient": "Gradient", "curl": "Curl",
                    "divergence": "Divergence",
                    "eigenvalues": "Eigenvalues"}


def derived_field(path: str | Path, operation: str, *,
                  view: str | int | None = None,
                  check_point: list[float] | None = None,
                  out_file: str | Path | None = None,
                  timeout_s: float = 300.0) -> dict[str, Any]:
    """Create a derived view: gradient, curl, divergence, or eigenvalues.

    The ParaView "Gradient of Unstructured DataSet" analog.  The result
    is the exact derivative of the piecewise-linear interpolant --
    element-wise constant per simplex for P1 data, discontinuous across
    elements.  ``eigenvalues`` (tensor views) writes THREE views
    (min/mid/max) into one file.  ``check_point`` additionally probes
    the first result view there (quick sanity check).
    """
    err = _check_input(path)
    if err:
        return err
    op_key = str(operation).lower()
    if op_key not in _DERIVED_PLUGINS:
        return {"ok": False,
                "error": f"unknown operation {operation!r}; available: "
                         f"{sorted(_DERIVED_PLUGINS)}"}
    return _run_post({"op": "derived", "path": str(path), "view": view,
                      "plugin": _DERIVED_PLUGINS[op_key],
                      "check_point": check_point,
                      "out_file": _default_out(path, out_file, op_key)},
                     timeout_s)


def threshold(path: str | Path, min_val: float, max_val: float, *,
              view: str | int | None = None, time_step: int = 0,
              dimension: int = -1,
              out_file: str | Path | None = None,
              timeout_s: float = 300.0) -> dict[str, Any]:
    """Keep only elements whose MEAN view value lies in [min_val, max_val].

    The ParaView Threshold analog (Plugin ExtractElements).  Selection
    uses the ELEMENT MEAN of the scalar at ``time_step`` (measured
    semantics) -- a tet with nodal values 1..4 is selected by 2.5.  For
    vector views threshold |v| via math_eval first.
    """
    err = _check_input(path)
    if err:
        return err
    return _run_post({"op": "threshold", "path": str(path), "view": view,
                      "min_val": float(min_val), "max_val": float(max_val),
                      "time_step": int(time_step),
                      "dimension": int(dimension),
                      "out_file": _default_out(path, out_file, "thresh")},
                     timeout_s)


def extract_skin(path: str | Path, *, view: str | int | None = None,
                 from_mesh: bool = False,
                 out_file: str | Path | None = None,
                 timeout_s: float = 300.0) -> dict[str, Any]:
    """Extract the boundary skin of a volume view (ParaView ExtractSurface).

    Produces the boundary triangles/quads of the view's volume elements
    with the field interpolated on them; ``from_mesh=True`` skins the
    model mesh instead of the view data.
    """
    err = _check_input(path)
    if err:
        return err
    return _run_post({"op": "skin", "path": str(path), "view": view,
                      "from_mesh": int(bool(from_mesh)),
                      "out_file": _default_out(path, out_file, "skin")},
                     timeout_s)


def mirror_expand(path: str | Path, planes: list[str], *,
                  parity: str = "scalar",
                  view: str | int | None = None,
                  origin: list[float] | None = None,
                  result_name: str = "mirrored",
                  out_file: str | Path | None = None,
                  timeout_s: float = 600.0) -> dict[str, Any]:
    """Expand a half/quarter/eighth model view by mirror symmetry.

    The ParaView Reflect analog with the field physics done right:
    every subset of ``planes`` ("x" mirrors across the plane x=origin_x,
    ...) produces one mirrored copy, and the data components transform
    by ``parity``:

    - "scalar": values copied unchanged (phi, |B|, T)
    - "vector": polar vector, v' = M v (A, J, force density)
    - "pseudovector": axial vector, v' = det(M) M v  (B and H!)

    All copies are merged into ONE view named ``result_name``.  Mirrored
    element orientation is repaired (SwapOrientation) so the expanded
    view renders and probes correctly on both sides.
    """
    err = _check_input(path)
    if err:
        return err
    planes = [str(p).lower() for p in planes]
    bad = [p for p in planes if p not in ("x", "y", "z")]
    if bad or not planes or len(planes) != len(set(planes)):
        return {"ok": False,
                "error": f"planes must be a non-repeating subset of "
                         f"['x','y','z'], got {planes}"}
    if parity not in ("scalar", "vector", "pseudovector"):
        return {"ok": False,
                "error": f"parity must be scalar|vector|pseudovector, "
                         f"got {parity!r}"}
    org = [float(c) for c in (origin or [0.0, 0.0, 0.0])]
    if len(org) != 3:
        return {"ok": False, "error": "origin needs exactly 3 coords"}
    return _run_post({"op": "mirror_expand", "path": str(path),
                      "view": view, "planes": planes, "parity": parity,
                      "origin": org, "result_name": str(result_name),
                      "out_file": _default_out(path, out_file, "full")},
                     timeout_s)


def transform_view(path: str | Path, matrix: list[float], *,
                   translation: list[float] | None = None,
                   view: str | int | None = None,
                   value_expressions: list[str] | None = None,
                   swap_orientation: bool = False,
                   result_name: str = "transformed",
                   out_file: str | Path | None = None,
                   timeout_s: float = 300.0) -> dict[str, Any]:
    """Apply an affine transform x' = A x + t to a COPY of a view.

    The ParaView Transform analog.  The input view is first
    materialized (Plugin Transform works IN PLACE, measured), so the
    source file stays valid.  ``value_expressions`` optionally rewrites
    the data during the copy (v0..v8 syntax) -- e.g. rotate vector
    components to match a rotated geometry, which Plugin(Transform)
    itself does NOT do.  Set ``swap_orientation=True`` when det(A) < 0.
    """
    err = _check_input(path)
    if err:
        return err
    A = [float(v) for v in matrix]
    if len(A) != 9:
        return {"ok": False,
                "error": "matrix must have 9 entries (row-major 3x3)"}
    t = [float(v) for v in (translation or [0.0, 0.0, 0.0])]
    if len(t) != 3:
        return {"ok": False, "error": "translation needs exactly 3 coords"}
    return _run_post({"op": "transform_affine", "path": str(path),
                      "view": view, "matrix": A, "translation": t,
                      "value_expressions": value_expressions,
                      "swap_orientation": int(bool(swap_orientation)),
                      "result_name": str(result_name),
                      "out_file": _default_out(path, out_file, "xform")},
                     timeout_s)


def warp_view(path: str | Path, factor: float = 1.0, *,
              view: str | int | None = None, time_step: int = 0,
              out_file: str | Path | None = None,
              timeout_s: float = 300.0) -> dict[str, Any]:
    """Displace a vector view's geometry by factor * its own vectors.

    The ParaView WarpByVector analog (deformed-shape display for
    displacement / force-density fields).  The warped view is
    materialized to a self-contained .pos for rendering.
    """
    err = _check_input(path)
    if err:
        return err
    return _run_post({"op": "warp", "path": str(path), "view": view,
                      "factor": float(factor), "time_step": int(time_step),
                      "out_file": _default_out(path, out_file, "warp")},
                     timeout_s)


def smooth_to_nodes(path: str | Path, *, view: str | int | None = None,
                    out_file: str | Path | None = None,
                    timeout_s: float = 300.0) -> dict[str, Any]:
    """Average element-wise data to nodes (CellDataToPointData analog).

    Plugin(Smooth): each node receives the mean of the adjacent
    elements' values (measured: elements 10/20 -> shared node 15).
    Use on ElementData views (per-element loss density, |J| per cell)
    before probing or isosurfacing, which need nodal continuity.
    """
    err = _check_input(path)
    if err:
        return err
    return _run_post({"op": "smooth", "path": str(path), "view": view,
                      "out_file": _default_out(path, out_file, "nodal")},
                     timeout_s)


def view_min_max(path: str | Path, *, view: str | int | None = None,
                 over_time: bool = False,
                 timeout_s: float = 300.0) -> dict[str, Any]:
    """Locate the min and max of a scalar view WITH their coordinates.

    The ParaView "find data / spreadsheet max" analog: returns
    {"min": {point, values}, "max": {point, values}} (values per time
    step).  For vector views build |v| with math_eval first.
    """
    err = _check_input(path)
    if err:
        return err
    return _run_post({"op": "min_max", "path": str(path), "view": view,
                      "over_time": int(bool(over_time))}, timeout_s)


def modulus_phase(path: str | Path, *, view: str | int | None = None,
                  real_step: int = 0, imag_step: int = 1,
                  out_file: str | Path | None = None,
                  timeout_s: float = 300.0) -> dict[str, Any]:
    """Convert a complex re/im two-step view to modulus and phase steps.

    AC post companion to harmonic_to_time: step 0 becomes |Z| =
    sqrt(re^2+im^2), step 1 becomes atan2(im, re) (measured: 3/4 ->
    5, 0.9273).  Amplitude maps of eddy-current solutions in one call.
    """
    err = _check_input(path)
    if err:
        return err
    return _run_post({"op": "modulus_phase", "path": str(path),
                      "view": view, "real_step": int(real_step),
                      "imag_step": int(imag_step),
                      "out_file": _default_out(path, out_file, "modphase")},
                     timeout_s)


def curve_profile(path: str | Path, x_expr: str, y_expr: str, z_expr: str,
                  u_min: float, u_max: float, n: int = 100, *,
                  view: str | int | None = None,
                  plot_png: str | Path | None = None,
                  csv_out: str | Path | None = None,
                  out_file: str | Path | None = None,
                  timeout_s: float = 300.0) -> dict[str, Any]:
    """Sample a view along a parametric curve x(u), y(u), z(u).

    The ParaView PlotOverLine generalized to curves (Plugin
    CutParametric, MathEx syntax: ``0.05*Cos(u)``).  THE tool for
    air-gap profiles: B(theta) on a circle needs one call with
    x=r*Cos(u), y=r*Sin(u), u in [0, 2*Pi].  Optionally writes a PNG
    graph (value / |v| vs u), a CSV, and an SL line view for rendering.
    """
    err = _check_input(path)
    if err:
        return err
    n = max(2, int(n))
    cfg: dict[str, Any] = {
        "op": "curve_profile", "path": str(path), "view": view,
        "x_expr": str(x_expr), "y_expr": str(y_expr),
        "z_expr": str(z_expr), "u_min": float(u_min),
        "u_max": float(u_max), "n": n,
        "out_file": (str(Path(out_file)) if out_file is not None else
                     _default_out(path, None, "curve")),
    }
    result = _run_post(cfg, timeout_s)
    if not result.get("ok"):
        return result
    pts = result.get("points", [])
    if len(pts) == n:
        us = [u_min + (u_max - u_min) * i / (n - 1) for i in range(n)]
    else:
        us = list(range(len(pts)))
        result["note"] = (
            f"{len(pts)} of {n} requested samples returned (curve leaves "
            f"the data region); the parameter column is the sample index")
    result["u"] = us
    if csv_out is not None:
        result["csv"] = _write_samples_csv(
            csv_out, us, pts, result["values"], param_name="u")
    if plot_png is not None:
        series = _series_from_samples(us, result["values"])
        plot = _plot_series(series, str(plot_png), xlabel="u",
                            ylabel="value", title="curve profile")
        result["plot_png"] = plot.get("png")
        if not plot.get("ok"):
            result["plot_error"] = plot.get("error")
    return result


def resample_grid(path: str | Path, origin: list[float],
                  u_point: list[float], v_point: list[float],
                  w_point: list[float], nu: int, nv: int, nw: int, *,
                  view: str | int | None = None,
                  csv_out: str | Path | None = None,
                  out_file: str | Path | None = None,
                  timeout_s: float = 600.0) -> dict[str, Any]:
    """Resample a view on a regular grid spanned by three box edges.

    The ParaView ResampleToImage analog (Plugin CutBox): origin +
    endpoints of the U/V/W edges, nu x nv x nw samples (W varies
    fastest in the returned order).  Feeds uniform-grid exports (CSV
    for numpy/MATLAB) from any unstructured result.
    """
    err = _check_input(path)
    if err:
        return err
    pts = {"origin": origin, "u_point": u_point, "v_point": v_point,
           "w_point": w_point}
    clean: dict[str, list[float]] = {}
    for key, val in pts.items():
        v = [float(c) for c in val]
        if len(v) != 3:
            return {"ok": False, "error": f"{key} needs exactly 3 coords"}
        clean[key] = v
    result = _run_post({"op": "resample_grid", "path": str(path),
                        "view": view, **clean,
                        "nu": max(1, int(nu)), "nv": max(1, int(nv)),
                        "nw": max(1, int(nw)),
                        "out_file": (str(Path(out_file))
                                     if out_file is not None else None)},
                       timeout_s)
    if not result.get("ok"):
        return result
    if csv_out is not None:
        result["csv"] = _write_samples_csv(
            csv_out, list(range(len(result["points"]))), result["points"],
            result["values"], param_name="index")
    return result


# ======================================================================
# Pure-Python post verbs (no gmsh needed: our own MSH v4.1 parser)
# ======================================================================

def _write_samples_csv(csv_out: str | Path, params: list[Any],
                       points: list[list[float]],
                       values: list[list[list[float]]],
                       param_name: str = "u") -> str:
    import csv as _csv
    out = Path(csv_out)
    out.parent.mkdir(parents=True, exist_ok=True)
    n_steps = len(values[0]) if values else 0
    ncomp = len(values[0][0]) if n_steps else 0
    header = [param_name, "x", "y", "z"] + [
        f"v_s{s}_c{c}" for s in range(n_steps) for c in range(ncomp)]
    with open(out, "w", newline="", encoding="utf-8") as f:
        w = _csv.writer(f)
        w.writerow(header)
        for p, pt, vals in zip(params, points, values):
            row = [p, *pt]
            for step_vals in vals:
                row.extend(step_vals)
            w.writerow(row)
    return str(out)


def _series_from_samples(xs: list[float],
                         values: list[list[list[float]]]
                         ) -> list[dict[str, Any]]:
    """Per-step series; scalar -> value, multi-component -> magnitude."""
    series: list[dict[str, Any]] = []
    if not values:
        return series
    n_steps = len(values[0])
    ncomp = len(values[0][0]) if n_steps else 0
    for s in range(n_steps):
        ys = []
        for vals in values:
            row = vals[s]
            ys.append(row[0] if ncomp == 1
                      else sum(v * v for v in row) ** 0.5)
        label = f"step {s}" + ("" if ncomp == 1 else " |v|")
        series.append({"x": xs, "y": ys, "label": label})
    return series


_PLOT_SCRIPT = r"""
import json
import sys

cfg_path, out_path = sys.argv[1], sys.argv[2]
with open(cfg_path, encoding="utf-8") as f:
    cfg = json.load(f)
result = {"ok": False}
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7, 4.5))
    if cfg["kind"] == "bar":
        centers, counts, width = cfg["centers"], cfg["counts"], cfg["width"]
        ax.bar(centers, counts, width=width * 0.95, align="center")
    else:
        for s in cfg["series"]:
            xs, ys = s["x"], s["y"]
            ax.plot(xs, ys, marker="." if len(xs) < 60 else None,
                    label=s.get("label"))
        if 1 < len(cfg["series"]) <= 8:
            ax.legend(fontsize=8)
    ax.set_xlabel(cfg.get("xlabel", ""))
    ax.set_ylabel(cfg.get("ylabel", ""))
    if cfg.get("title"):
        ax.set_title(cfg["title"], fontsize=10)
    ax.grid(True, alpha=0.4)
    fig.tight_layout()
    fig.savefig(cfg["png"], dpi=130)
    result = {"ok": True, "png": cfg["png"]}
except Exception as exc:
    result["error"] = f"{type(exc).__name__}: {exc}"
with open(out_path, "w", encoding="utf-8") as f:
    json.dump(result, f)
"""


def _plot_series(series: list[dict[str, Any]], png: str, *,
                 xlabel: str = "", ylabel: str = "", title: str = "",
                 timeout_s: float = 120.0) -> dict[str, Any]:
    out = Path(png)
    out.parent.mkdir(parents=True, exist_ok=True)
    cfg = {"kind": "line", "series": series, "xlabel": xlabel,
           "ylabel": ylabel, "title": title, "png": str(out)}
    with tempfile.TemporaryDirectory(prefix="radia_mcp_gmsh_plot_") as work:
        cfg_path = Path(work) / "plot.json"
        cfg_path.write_text(json.dumps(cfg), encoding="utf-8")
        return run_gmsh_json_subprocess(
            _PLOT_SCRIPT, [str(cfg_path)],
            timeout_s=timeout_s, prefix="radia_mcp_gmsh_plot_")


def export_view_csv(path: str | Path, csv_out: str | Path, *,
                    view: str | None = None,
                    kind: str = "auto") -> dict[str, Any]:
    """Dump view data to CSV (the ParaView spreadsheet / SaveData analog).

    Pure Python on our own MSH v4.1 parser -- no gmsh needed.
    ``kind="nodes"`` writes tag,x,y,z + one column per view step and
    component from $NodeData; ``kind="elements"`` writes element tag +
    CENTROID coordinates + $ElementData columns.  ``auto`` picks nodes
    when NodeData exists, else elements.  List-based .pos files have no
    node table -- resample_grid/curve_profile with csv_out cover those.
    """
    from .msh_inspect import read_msh_data
    src = Path(path)
    if not src.is_file():
        return {"ok": False, "error": f"file not found: {src}"}
    try:
        data = read_msh_data(src, include_elements=True)
    except ValueError as exc:
        return {"ok": False, "error": str(exc)}
    node_views = [v for v in data["views"] if v["section"] == "NodeData"]
    elem_views = [v for v in data["views"] if v["section"] == "ElementData"]
    if view is not None:
        node_views = [v for v in node_views if v["name"] == view]
        elem_views = [v for v in elem_views if v["name"] == view]
    if kind == "auto":
        kind = "nodes" if node_views else "elements"
    if kind == "nodes":
        chosen, coords = node_views, data["nodes"]
    elif kind == "elements":
        chosen = elem_views
        coords = {}
        for tag, el in data["elements"].items():
            pts = [data["nodes"][r] for r in el["nodes"]
                   if r in data["nodes"]]
            if pts:
                coords[tag] = [sum(p[k] for p in pts) / len(pts)
                               for k in range(3)]
    else:
        return {"ok": False, "error": f"kind must be auto|nodes|elements, "
                                      f"got {kind!r}"}
    if not chosen:
        avail = sorted({v["name"] for v in data["views"]})
        return {"ok": False,
                "error": f"no matching {kind} data sections"
                         + (f" for view {view!r}" if view else "")
                         + f"; available views: {avail}"}

    import csv as _csv
    out = Path(csv_out)
    out.parent.mkdir(parents=True, exist_ok=True)
    columns: list[tuple[dict[str, Any], int]] = []
    header = ["tag", "x", "y", "z"]
    for v in chosen:
        ncomp = int(v["components"] or 1)
        for c in range(ncomp):
            header.append(f"{v['name']}_s{v['step']}_c{c}")
            columns.append((v, c))
    tags = sorted(coords)
    n_rows = 0
    with open(out, "w", newline="", encoding="utf-8") as f:
        w = _csv.writer(f)
        w.writerow(header)
        for tag in tags:
            row: list[Any] = [tag, *coords[tag]]
            for v, c in columns:
                vals = v["rows"].get(tag)
                row.append(vals[c] if vals and c < len(vals) else "")
            w.writerow(row)
            n_rows += 1
    return {"ok": True, "csv": str(out), "kind": kind, "n_rows": n_rows,
            "columns": header,
            "views": [{"name": v["name"], "step": v["step"],
                       "components": v["components"]} for v in chosen]}


def field_histogram(path: str | Path, *, view: str | None = None,
                    step: int | None = None, component: int | None = None,
                    bins: int = 32,
                    value_range: list[float] | None = None,
                    plot_png: str | Path | None = None) -> dict[str, Any]:
    """Histogram of a view's values (the ParaView Histogram analog).

    Pure Python on the MSH parser.  Scalars bin the value, multi-
    component views bin the euclidean magnitude unless ``component``
    selects one.  ``step=None`` pools every time step; an int selects
    one.  Optional PNG bar chart.
    """
    from .msh_inspect import read_msh_data
    src = Path(path)
    if not src.is_file():
        return {"ok": False, "error": f"file not found: {src}"}
    try:
        data = read_msh_data(src)
    except ValueError as exc:
        return {"ok": False, "error": str(exc)}
    views = data["views"]
    if view is not None:
        views = [v for v in views if v["name"] == view]
    if step is not None:
        views = [v for v in views if v["step"] == step]
    if not views:
        avail = sorted({v["name"] for v in data["views"]})
        return {"ok": False,
                "error": "no matching data sections; available views: "
                         f"{avail}"}
    name = views[0]["name"]
    views = [v for v in views if v["name"] == name]

    samples: list[float] = []
    for v in views:
        ncomp = int(v["components"] or 1)
        for vals in v["rows"].values():
            for i in range(0, len(vals) - ncomp + 1, ncomp):
                chunk = vals[i:i + ncomp]
                if component is not None:
                    if component >= ncomp:
                        return {"ok": False,
                                "error": f"component {component} out of "
                                         f"range (view has {ncomp})"}
                    samples.append(chunk[component])
                elif ncomp == 1:
                    samples.append(chunk[0])
                else:
                    samples.append(sum(c * c for c in chunk) ** 0.5)
    finite = [s for s in samples if s == s and abs(s) != float("inf")]
    if not finite:
        return {"ok": False, "error": "view contains no finite samples"}
    lo, hi = ((float(value_range[0]), float(value_range[1]))
              if value_range else (min(finite), max(finite)))
    if hi <= lo:
        hi = lo + 1.0
    bins = max(1, int(bins))
    width = (hi - lo) / bins
    counts = [0] * bins
    for s in finite:
        k = int((s - lo) / width)
        if 0 <= k < bins:
            counts[k] += 1
        elif k == bins:  # right edge inclusive
            counts[-1] += 1
    mean = sum(finite) / len(finite)
    var = sum((s - mean) ** 2 for s in finite) / len(finite)
    result: dict[str, Any] = {
        "ok": True, "view": name, "n_samples": len(finite),
        "n_nonfinite": len(samples) - len(finite),
        "bin_edges": [lo + k * width for k in range(bins + 1)],
        "counts": counts,
        "stats": {"min": min(finite), "max": max(finite), "mean": mean,
                  "std": var ** 0.5},
    }
    if plot_png is not None:
        out = Path(plot_png)
        out.parent.mkdir(parents=True, exist_ok=True)
        cfg = {"kind": "bar",
               "centers": [lo + (k + 0.5) * width for k in range(bins)],
               "counts": counts, "width": width,
               "xlabel": name, "ylabel": "count",
               "title": f"histogram of {name}", "png": str(out)}
        with tempfile.TemporaryDirectory(
                prefix="radia_mcp_gmsh_plot_") as work:
            cfg_path = Path(work) / "plot.json"
            cfg_path.write_text(json.dumps(cfg), encoding="utf-8")
            plot = run_gmsh_json_subprocess(
                _PLOT_SCRIPT, [str(cfg_path)],
                timeout_s=120.0, prefix="radia_mcp_gmsh_plot_")
        result["plot_png"] = plot.get("png")
        if not plot.get("ok"):
            result["plot_error"] = plot.get("error")
    return result


def point_history(path: str | Path, point: list[float], *,
                  view: str | int | None = None,
                  plot_png: str | Path | None = None,
                  timeout_s: float = 300.0) -> dict[str, Any]:
    """Value of a view at one point across ALL time steps.

    The ParaView PlotDataOverTime analog for a probe point: returns the
    per-step values plus the step TIMES when the file records them
    ($NodeData headers), and optionally plots value vs time/step.
    """
    result = probe_field(path, [point], view=view, step=-1,
                         timeout_s=timeout_s)
    if not result.get("ok"):
        return result
    entry = result["views"][0]["points"][0]
    if not entry.get("found"):
        return {"ok": False,
                "error": f"point {point} lies outside the data (nearest "
                         f"element {entry.get('distance'):.3g} away)"}
    steps = entry.get("steps") or [entry["values"]]
    times: list[float] | None = None
    src = Path(path)
    if src.suffix.lower() == ".msh":
        from .msh_inspect import read_msh_data
        try:
            data = read_msh_data(src)
        except ValueError:
            data = None
        if data is not None:
            name = result["views"][0]["name"]
            per_step = {v["step"]: v["time"] for v in data["views"]
                        if v["name"] == name and v["time"] is not None}
            if len(per_step) == len(steps):
                times = [per_step[s] for s in sorted(per_step)]
    out: dict[str, Any] = {"ok": True, "point": point,
                           "view": result["views"][0]["name"],
                           "n_steps": len(steps), "steps": steps,
                           "times": times}
    if plot_png is not None:
        xs = times if times is not None else list(range(len(steps)))
        ncomp = len(steps[0])
        series = []
        for c in range(ncomp):
            series.append({"x": xs, "y": [s[c] for s in steps],
                           "label": f"c{c}"})
        if ncomp > 1:
            series.append({"x": xs,
                           "y": [sum(v * v for v in s) ** 0.5
                                 for s in steps],
                           "label": "|v|"})
        plot = _plot_series(series, str(plot_png),
                            xlabel="time" if times is not None else "step",
                            ylabel=out["view"], title="point history")
        out["plot_png"] = plot.get("png")
        if not plot.get("ok"):
            out["plot_error"] = plot.get("error")
    return out
