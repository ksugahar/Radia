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

import ast
import json
import keyword
import math
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


def _set_adaptive_extraction(tag, cfg):
    # Plugin(Isosurface) honors RecurLevel ONLY when the source view
    # carries adaptive visualization data (measured gmsh 4.15.2: on a
    # TET10 quadratic field the radial error drops 0.21 -> 0.008 once
    # View.AdaptVisualizationGrid is enabled before the plugin runs;
    # without it RecurLevel is silently ignored).
    import gmsh
    recur = int(cfg.get("recur_level", 0))
    if recur <= 0:
        return
    idx = gmsh.view.getIndex(tag)
    gmsh.option.setNumber(f"View[{idx}].AdaptVisualizationGrid", 1)
    gmsh.option.setNumber(f"View[{idx}].MaxRecursionLevel", recur)
    gmsh.option.setNumber(f"View[{idx}].TargetError",
                          float(cfg.get("target_error", 1e-4)))


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
            _set_adaptive_extraction(tag, cfg)
            gmsh.plugin.setNumber("Isosurface", "View",
                                  gmsh.view.getIndex(tag))
            gmsh.plugin.setNumber("Isosurface", "Value",
                                  float(cfg["value"]))
            gmsh.plugin.setNumber("Isosurface", "RecurLevel",
                                  int(cfg.get("recur_level", 0)))
            gmsh.plugin.setNumber("Isosurface", "TargetError",
                                  float(cfg.get("target_error", 1e-4)))
            out_tag = gmsh.plugin.run("Isosurface")
            dtypes, nels, _data = gmsh.view.getListData(out_tag)
            gmsh.view.write(out_tag, cfg["out_file"])
            # Detect the common open-shell case: the level set reaches the
            # model's outer bounding box.  Gmsh's Isosurface plugin emits
            # element-local polygons without shared topology, so this is
            # deliberately reported as an outer-boundary contact check;
            # internal openings cannot be classified from this view alone.
            bb = gmsh.model.getBoundingBox(-1, -1)
            span = max(bb[3] - bb[0], bb[4] - bb[1], bb[5] - bb[2], 1e-30)
            tol = 1e-6 * span
            on_outer_boundary = 0
            n_vertices = 0
            for kind, nel, arr in zip(dtypes, nels, _data):
                nv = {"SP": 1, "SL": 2, "ST": 3, "SQ": 4,
                      "VP": 1, "VL": 2, "VT": 3, "VQ": 4}.get(str(kind))
                if not nv or not nel:
                    continue
                stride = len(arr) // int(nel)
                for e in range(int(nel)):
                    base = e * stride
                    for v in range(nv):
                        px = arr[base + v]
                        py = arr[base + nv + v]
                        pz = arr[base + 2 * nv + v]
                        n_vertices += 1
                        if (abs(px - bb[0]) < tol or abs(px - bb[3]) < tol
                                or abs(py - bb[1]) < tol
                                or abs(py - bb[4]) < tol
                                or abs(pz - bb[2]) < tol
                                or abs(pz - bb[5]) < tol):
                            on_outer_boundary += 1
            result.update({"ok": True, "ran": True,
                           "out_file": cfg["out_file"],
                           "recur_level": int(cfg.get("recur_level", 0)),
                           "pieces": {str(t): int(n)
                                      for t, n in zip(dtypes, nels)},
                           "n_vertices": n_vertices,
                           "boundary_vertices": on_outer_boundary,
                           "outer_boundary_vertices": on_outer_boundary,
                           "touches_outer_boundary": bool(on_outer_boundary),
                           "open_surface_check": "outer_bbox_contact_only",
                           "open_surface": bool(on_outer_boundary)})
            if on_outer_boundary:
                result["note"] = (
                    "the isosurface is CUT OPEN at the OUTER MODEL BOUNDARY "
                    f"({on_outer_boundary} of {n_vertices} polygon vertices): "
                    "rendered semi-transparent you will see straight "
                    "through the opening. That is the geometry, not a "
                    "rendering artefact -- pick a level whose shell "
                    "closes inside the domain, or clip the view.")

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
            # Probe-driven arc-length tracer with curvature-ADAPTIVE RK4
            # (step_size is the MAXIMUM step; it halves while the turn
            # per step exceeds max_turn_deg), CLOSED-LOOP detection (the
            # magnetic field-line case), one merged polyline per seed,
            # and per-line termination reasons.  Plugin(StreamLines) on
            # this gmsh build only re-emits the seed points.
            import math as _math
            tag = _view_tag_by_selector(tags, cfg.get("view"))
            step = int(cfg.get("time_step", 0))
            ds0 = cfg.get("step_size")
            if ds0 is None:
                bb = gmsh.model.getBoundingBox(-1, -1)
                diag = sum((bb[i + 3] - bb[i]) ** 2 for i in range(3)) ** 0.5
                ds0 = diag / 200.0
            ds0 = float(ds0)
            ds_min = ds0 / 64.0
            adaptive = bool(cfg.get("adaptive", True))
            detect_closure = bool(cfg.get("closure", True))
            angle_max = _math.radians(float(cfg.get("max_turn_deg", 10.0)))
            max_steps = int(cfg.get("max_steps", 400))
            n_seeds = max(1, int(cfg.get("n_seeds", 10)))
            s0, s1 = cfg["seed_start"], cfg["seed_end"]
            seeds = [[s0[k] + (s1[k] - s0[k]) *
                      (i / (n_seeds - 1) if n_seeds > 1 else 0.0)
                      for k in range(3)] for i in range(n_seeds)]
            both = bool(cfg.get("both_directions", True))

            def _sample(q):
                values, dist = gmsh.view.probe(tag, q[0], q[1], q[2],
                                               step=step)
                # dist > 0 = outside the data (list-based vector views
                # report the nearest value at ANY distance -- without
                # this gate the tracer would never stop).
                if len(values) < 3 or float(dist) > 0.0:
                    return None, 0.0, "left_data"
                vx, vy, vz = (float(values[0]), float(values[1]),
                              float(values[2]))
                norm = (vx * vx + vy * vy + vz * vz) ** 0.5
                if norm < 1e-30:
                    return None, 0.0, "stagnation"
                return (vx / norm, vy / norm, vz / norm), norm, None

            def _trace(seed, sgn):
                t0, n0, why = _sample(seed)
                if t0 is None:
                    return [list(seed)], [n0], why, False

                def _f(q):
                    # returns (direction, None) or (None, reason) so a
                    # mid-step field zero reports "stagnation", not
                    # "left_data"
                    t, _n, w = _sample(q)
                    if t is None:
                        return None, w
                    return [sgn * c for c in t], None

                path = [list(seed)]
                norms = [n0]
                p = list(seed)
                start_dir = [sgn * c for c in t0]
                t_prev = start_dir[:]
                ds = ds0
                arc = 0.0
                reason = "max_steps"
                closed = False
                for _ in range(max_steps):
                    outcome = None
                    while True:
                        k1, w1 = _f(p)
                        if k1 is None:
                            outcome = w1
                            break
                        q2 = [p[i] + 0.5 * ds * k1[i] for i in range(3)]
                        k2, w2 = _f(q2)
                        if k2 is None:
                            outcome = w2
                            break
                        q3 = [p[i] + 0.5 * ds * k2[i] for i in range(3)]
                        k3, w3 = _f(q3)
                        if k3 is None:
                            outcome = w3
                            break
                        q4 = [p[i] + ds * k3[i] for i in range(3)]
                        k4, w4 = _f(q4)
                        if k4 is None:
                            outcome = w4
                            break
                        p_new = [p[i] + ds / 6.0 * (k1[i] + 2 * k2[i]
                                                    + 2 * k3[i] + k4[i])
                                 for i in range(3)]
                        t_new, n_new, why = _sample(p_new)
                        if t_new is None:
                            outcome = why
                            break
                        t_dir = [sgn * c for c in t_new]
                        d = sum(a * b for a, b in zip(t_prev, t_dir))
                        angle = _math.acos(max(-1.0, min(1.0, d)))
                        if (not adaptive or angle <= angle_max
                                or ds <= ds_min * 1.0001):
                            outcome = "ok" if (angle <= angle_max
                                               or not adaptive) \
                                else "incoherent"
                            break
                        ds *= 0.5
                    if outcome == "left_data":
                        reason = "left_data"
                        break
                    if outcome == "stagnation":
                        reason = "stagnation"
                        break
                    if outcome == "incoherent":
                        # ds_min reached and the direction still flips:
                        # the field direction is no longer coherent
                        # (approaching a zero / reversal point)
                        reason = "stagnation"
                        break
                    p = p_new
                    arc += ds
                    path.append(list(p))
                    norms.append(n_new)
                    t_prev = t_dir
                    if adaptive and angle < 0.5 * angle_max and ds < ds0:
                        ds = min(ds0, ds * 2.0)
                    if detect_closure and arc > 4.0 * ds0:
                        dd = sum((a - b) ** 2
                                 for a, b in zip(p, seed)) ** 0.5
                        if (dd < max(ds, 0.75 * ds0)
                                and sum(a * b for a, b in
                                        zip(t_dir, start_dir)) > 0.8):
                            path.append(list(seed))
                            norms.append(norms[0])
                            closed = True
                            reason = "closed"
                            break
                return path, norms, reason, closed

            polylines = []
            skipped_seeds = 0
            for seed in seeds:
                fwd_path, fwd_norms, fwd_reason, closed = _trace(seed, 1)
                if closed or not both:
                    if len(fwd_path) > 1:
                        polylines.append((fwd_path, fwd_norms, closed,
                                          {"forward": fwd_reason}))
                    else:
                        skipped_seeds += 1
                    continue
                bwd_path, bwd_norms, bwd_reason, _bc = _trace(seed, -1)
                merged = list(reversed(bwd_path)) + fwd_path[1:]
                merged_norms = (list(reversed(bwd_norms))
                                + fwd_norms[1:])
                if len(merged) > 1:
                    polylines.append((merged, merged_norms, False,
                                      {"forward": fwd_reason,
                                       "backward": bwd_reason}))
                else:
                    skipped_seeds += 1

            sl_data = []
            n_segments = 0
            for path, norms, _cl, _rs in polylines:
                for i in range(len(path) - 1):
                    a, b = path[i], path[i + 1]
                    sl_data += [a[0], b[0], a[1], b[1], a[2], b[2],
                                norms[i], norms[i + 1]]
                    n_segments += 1
            out_tag = gmsh.view.add("streamlines")
            if n_segments:
                gmsh.view.addListData(out_tag, "SL", n_segments, sl_data)
            gmsh.view.write(out_tag, cfg["out_file"])

            arrows_every = int(cfg.get("arrows_every", 0))
            n_arrows = 0
            if arrows_every > 0:
                vp_data = []
                for path, norms, _cl, _rs in polylines:
                    for i in range(1, len(path) - 1, arrows_every):
                        t = [path[i + 1][k] - path[i - 1][k]
                             for k in range(3)]
                        tn = sum(c * c for c in t) ** 0.5
                        if tn < 1e-30:
                            continue
                        vp_data += list(path[i]) + [
                            c / tn * norms[i] for c in t]
                        n_arrows += 1
                arrow_tag = gmsh.view.add("streamline_arrows")
                if n_arrows:
                    gmsh.view.addListData(arrow_tag, "VP", n_arrows,
                                          vp_data)
                gmsh.view.write(arrow_tag, cfg["out_file"], append=True)

            reason_counts = {}
            for _path, _norms, _cl, rs in polylines:
                for r in rs.values():
                    reason_counts[r] = reason_counts.get(r, 0) + 1
            result.update({
                "ok": True, "ran": True,
                "out_file": cfg["out_file"],
                "n_polylines": len(polylines),
                "n_segments": n_segments,
                "n_closed": sum(1 for _p, _n, cl, _r in polylines if cl),
                "n_arrows": n_arrows,
                "skipped_seeds": skipped_seeds,
                "reasons": reason_counts,
                "step_size": ds0,
                "lines": [
                    {"n_points": len(path), "closed": cl, "reasons": rs,
                     "arc_length": sum(
                         sum((path[i + 1][k] - path[i][k]) ** 2
                             for k in range(3)) ** 0.5
                         for i in range(len(path) - 1))}
                    for path, _norms, cl, rs in polylines],
                "polylines": [
                    {"points": path, "field_norms": norms,
                     "closed": cl, "reasons": rs}
                    for path, norms, cl, rs in polylines
                ] if cfg.get("return_points") else None,
            })

        elif op == "particle_trace":
            # Relativistic Boris pusher over the probed field views:
            # the charged-particle ORBIT verb (dp/dt = q(E + v x B)),
            # distinct from "streamlines" (massless tangent curves of
            # the field itself).  The magnetic rotation is an exact
            # isometry of u = gamma*v, so in a pure-B trace the
            # per-track speed_change_rel is an integrator health
            # metric (~1e-15); kinetic energy changes only through E.
            import math as _math
            C_LIGHT = 299792458.0
            J_PER_EV = 1.602176634e-19
            b_tag = _view_tag_by_selector(tags, cfg.get("view"))
            e_tag = (None if cfg.get("e_view") is None else
                     _view_tag_by_selector(tags, cfg["e_view"]))
            step = int(cfg.get("time_step", 0))
            q = float(cfg["charge_c"])
            m = float(cfg["mass_kg"])
            ke0_j = float(cfg["kinetic_energy_ev"]) * J_PER_EV
            gamma0 = 1.0 + ke0_j / (m * C_LIGHT * C_LIGHT)
            u0_abs = C_LIGHT * _math.sqrt(gamma0 * gamma0 - 1.0)
            d = cfg["direction"]
            dn = _math.sqrt(sum(c * c for c in d))
            dirn = [c / dn for c in d]
            seeds = cfg["seeds"]
            spg = int(cfg.get("steps_per_gyration", 64))
            max_steps = int(cfg.get("max_steps", 20000))
            max_time = cfg.get("max_time_s")
            color_by = cfg.get("color_by", "time")
            dt_fixed = cfg.get("dt_s")
            arrows_every = int(cfg.get("arrows_every", 0))

            def _probe_vec(tag_, p):
                vals, dist = gmsh.view.probe(tag_, p[0], p[1], p[2],
                                             step=step)
                # dist > 0 = outside the data (list-based vector views
                # report the nearest value at ANY distance).
                if len(vals) < 3 or float(dist) > 0.0:
                    return None
                return [float(vals[0]), float(vals[1]), float(vals[2])]

            def _cross(a, b):
                return [a[1] * b[2] - a[2] * b[1],
                        a[2] * b[0] - a[0] * b[2],
                        a[0] * b[1] - a[1] * b[0]]

            def _gamma_of(u):
                return _math.sqrt(
                    1.0 + (u[0] * u[0] + u[1] * u[1] + u[2] * u[2])
                    / (C_LIGHT * C_LIGHT))

            def _color(t, u):
                g = _gamma_of(u)
                if color_by == "speed":
                    return _math.sqrt(sum(c * c for c in u)) / g
                if color_by == "energy":
                    return (g - 1.0) * m * C_LIGHT * C_LIGHT / J_PER_EV
                return t

            tracks_out = []
            polylines = []
            skipped_seeds = 0
            for seed in seeds:
                B0 = _probe_vec(b_tag, seed)
                if B0 is None:
                    skipped_seeds += 1
                    tracks_out.append({"seed": list(seed), "n_steps": 0,
                                       "reason": "seed_outside_data"})
                    continue
                b0_mag = _math.sqrt(sum(c * c for c in B0))
                if dt_fixed is not None:
                    dt = float(dt_fixed)
                elif b0_mag > 0.0:
                    dt = (2.0 * _math.pi * gamma0 * m
                          / (abs(q) * b0_mag) / spg)
                else:
                    raise RuntimeError(
                        f"|B| = 0 at seed {list(seed)}: no gyration "
                        "period to derive the time step from -- pass "
                        "dt_s explicitly")
                # gyroradius from the momentum component perp to B(seed)
                u_vec = [u0_abs * c for c in dirn]
                if b0_mag > 0.0:
                    bh = [c / b0_mag for c in B0]
                    upar = sum(a * b for a, b in zip(u_vec, bh))
                    uperp = _math.sqrt(max(u0_abs**2 - upar**2, 0.0))
                    r_gyro = m * uperp / (abs(q) * b0_mag)
                else:
                    r_gyro = None
                p = list(seed)
                t_now = 0.0
                path_len = 0.0
                pts = [list(p)]
                cols = [_color(t_now, u_vec)]
                vels = [[c / gamma0 for c in u_vec]]
                reason = "max_steps"
                for _ in range(max_steps):
                    if e_tag is not None:
                        E = _probe_vec(e_tag, p)
                        if E is None:
                            reason = "left_data"
                            break
                    else:
                        E = [0.0, 0.0, 0.0]
                    B = _probe_vec(b_tag, p)
                    if B is None:
                        reason = "left_data"
                        break
                    qmdt2 = q * dt / (2.0 * m)
                    um = [u_vec[k] + qmdt2 * E[k] for k in range(3)]
                    g_minus = _gamma_of(um)
                    f = q * dt / (2.0 * m * g_minus)
                    tv = [f * B[k] for k in range(3)]
                    t2 = sum(c * c for c in tv)
                    upr = [um[k] + _cross(um, tv)[k] for k in range(3)]
                    sv = [2.0 * c / (1.0 + t2) for c in tv]
                    upl = [um[k] + _cross(upr, sv)[k] for k in range(3)]
                    u_vec = [upl[k] + qmdt2 * E[k] for k in range(3)]
                    g_new = _gamma_of(u_vec)
                    v = [c / g_new for c in u_vec]
                    p = [p[k] + v[k] * dt for k in range(3)]
                    t_now += dt
                    v_abs = _math.sqrt(sum(c * c for c in v))
                    path_len += v_abs * dt
                    pts.append(list(p))
                    cols.append(_color(t_now, u_vec))
                    vels.append(v)
                    if max_time is not None and t_now >= float(max_time):
                        reason = "max_time"
                        break
                g_end = _gamma_of(u_vec)
                v0_abs = u0_abs / gamma0
                v_end_abs = _math.sqrt(sum(c * c for c in u_vec)) / g_end
                tracks_out.append({
                    "seed": list(seed),
                    "dt_s": dt,
                    "n_steps": len(pts) - 1,
                    "time_s": t_now,
                    "path_length_m": path_len,
                    "b_seed_t": b0_mag,
                    "gyroradius_seed_m": r_gyro,
                    "ke_start_ev": float(cfg["kinetic_energy_ev"]),
                    "ke_end_ev": (g_end - 1.0) * m * C_LIGHT * C_LIGHT
                    / J_PER_EV,
                    "speed_change_rel": abs(v_end_abs - v0_abs)
                    / max(v0_abs, 1e-300),
                    "reason": reason,
                })
                if len(pts) > 1:
                    polylines.append((pts, cols, vels, dt))

            sl_data = []
            n_segments = 0
            for pts, cols, _vels, _dt in polylines:
                for i in range(len(pts) - 1):
                    a, b = pts[i], pts[i + 1]
                    sl_data += [a[0], b[0], a[1], b[1], a[2], b[2],
                                cols[i], cols[i + 1]]
                    n_segments += 1
            out_tag = gmsh.view.add("particle_tracks")
            if n_segments:
                gmsh.view.addListData(out_tag, "SL", n_segments, sl_data)
            gmsh.view.write(out_tag, cfg["out_file"])

            # Beam animation: a MULTI-STEP copy whose step k shows only
            # the part of each orbit the particle has already flown.
            # Rendered with a FIXED colour range and SaturateValues=0,
            # gmsh does not draw an element whose value is outside the
            # range, so a far out-of-range SENTINEL hides the future
            # while the colour scale stays put (sweeping the range
            # instead would rescale the colours every frame).  The
            # sentinel must be FINITE: gmsh writes a NaN as the literal
            # "nan" and its own parser then rejects the file.
            n_anim = int(cfg.get("animation_frames", 0))
            animation = None
            if n_anim and n_segments:
                # Gate on ABSOLUTE flight time, not on the fraction of
                # each track: seeds sitting in different |B| get
                # different time steps, and only a shared clock makes
                # the frames a physical snapshot of the whole beam.
                flat = []          # (p_a, p_b, colour, t_mid)
                for pts, cols, _v, dt_track in polylines:
                    for i in range(len(pts) - 1):
                        flat.append((pts[i], pts[i + 1],
                                     0.5 * (cols[i] + cols[i + 1]),
                                     (i + 0.5) * dt_track))
                c_vals = [row[2] for row in flat]
                lo, hi = min(c_vals), max(c_vals)
                # Far outside on the scale of the VALUES, not of their
                # spread: a monoenergetic beam coloured by energy has a
                # spread of ~1e-13, and a sentinel one spread above the
                # maximum would carry the beam's own colour.
                scale = max(hi - lo, abs(lo), abs(hi), 1.0)
                if hi - lo < 1e-12 * scale:
                    hi = lo + 1e-3 * scale
                sentinel = hi + 10.0 * scale
                t_end = max(row[3] for row in flat)
                mode = cfg.get("animation_mode", "trail")
                window = float(cfg.get("comet_window", 0.15)) * t_end
                anim_data = []
                for a, b, cval, t_mid in flat:
                    row = [a[0], b[0], a[1], b[1], a[2], b[2]]
                    for k in range(n_anim):
                        head = t_end * (k + 1) / n_anim
                        vis = t_mid <= head
                        if mode == "comet":
                            vis = vis and t_mid >= head - window
                        v = cval if vis else sentinel
                        row += [v, v]
                    anim_data += row
                anim_name = f"beam ({color_by})"
                anim_tag = gmsh.view.add(anim_name)
                gmsh.view.addListData(anim_tag, "SL", len(flat), anim_data)
                gmsh.view.write(anim_tag, cfg["out_file"], append=True)
                animation = {
                    "view": anim_name, "n_steps": n_anim, "mode": mode,
                    "flight_time_s": t_end,
                    "color_range": [lo, hi], "sentinel": sentinel,
                    "render_hint": {"color": {"range": [lo, hi],
                                              "saturate": False}},
                }

            n_arrows = 0
            if arrows_every > 0:
                vp_data = []
                for pts, _cols, vels, _dt in polylines:
                    for i in range(0, len(pts), arrows_every):
                        vn = _math.sqrt(sum(c * c for c in vels[i]))
                        if vn < 1e-30:
                            continue
                        vp_data += list(pts[i]) + list(vels[i])
                        n_arrows += 1
                arrow_tag = gmsh.view.add("particle_arrows")
                if n_arrows:
                    gmsh.view.addListData(arrow_tag, "VP", n_arrows,
                                          vp_data)
                gmsh.view.write(arrow_tag, cfg["out_file"], append=True)

            reason_counts = {}
            for tr in tracks_out:
                reason_counts[tr["reason"]] = (
                    reason_counts.get(tr["reason"], 0) + 1)
            result.update({
                "ok": True, "ran": True,
                "out_file": cfg["out_file"],
                "n_tracks": len(polylines),
                "n_segments": n_segments,
                "n_arrows": n_arrows,
                "skipped_seeds": skipped_seeds,
                "reasons": reason_counts,
                "color_by": color_by,
                "animation": animation,
                "tracks": tracks_out,
                "polylines": [
                    {"points": pts, "colors": cols}
                    for pts, cols, _v, _dt in polylines
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

        elif op == "flux_lines":
            # N-level isoline/isosurface extraction merged into ONE
            # view.  On a 2D scalar (A_z, the axisymmetric flux psi)
            # the equally spaced levels ARE the exact equal-flux field
            # lines (the FEMM-style motor plot) -- integration-free.
            tag = _view_tag_by_selector(tags, cfg.get("view"))
            levels = cfg.get("levels")
            if levels is None:
                gmsh.plugin.setNumber("MinMax", "View",
                                      gmsh.view.getIndex(tag))
                gmsh.plugin.setNumber("MinMax", "OverTime", 0)
                gmsh.plugin.setNumber("MinMax", "Argument", 1)
                before = set(gmsh.view.getTags())
                gmsh.plugin.run("MinMax")
                lo = hi = None
                for t in list(gmsh.view.getTags()):
                    if t in before:
                        continue
                    nm = _view_name(t)
                    _dt, _ne, data = gmsh.view.getListData(t)
                    val = (float(data[0][3])
                           if data and len(data[0]) >= 4 else None)
                    if nm.endswith("_Min"):
                        lo = val
                    elif nm.endswith("_Max"):
                        hi = val
                    gmsh.view.remove(t)
                if lo is None or hi is None or hi <= lo:
                    raise RuntimeError(
                        f"cannot derive levels from view range "
                        f"[{lo}, {hi}]")
                n = max(1, int(cfg.get("n_levels", 20)))
                levels = [lo + (hi - lo) * (k + 1) / (n + 1)
                          for k in range(n)]
            levels = [float(lv) for lv in levels]
            _set_adaptive_extraction(tag, cfg)
            acc = {}
            pieces_per_level = []
            for level in levels:
                gmsh.plugin.setNumber("Isosurface", "View",
                                      gmsh.view.getIndex(tag))
                gmsh.plugin.setNumber("Isosurface", "Value", level)
                gmsh.plugin.setNumber("Isosurface", "RecurLevel",
                                      int(cfg.get("recur_level", 0)))
                gmsh.plugin.setNumber("Isosurface", "TargetError",
                                      float(cfg.get("target_error", 1e-4)))
                iso_tag = gmsh.plugin.run("Isosurface")
                dtypes, nels, data = gmsh.view.getListData(iso_tag)
                per = {}
                for dt_, ne_, block in zip(dtypes, nels, data):
                    key = str(dt_)
                    per[key] = per.get(key, 0) + int(ne_)
                    old_n, old_data = acc.get(key, (0, []))
                    acc[key] = (old_n + int(ne_),
                                old_data + [float(v) for v in block])
                pieces_per_level.append(per)
                gmsh.view.remove(iso_tag)
            out_tag = gmsh.view.add(cfg.get("result_name", "flux_lines"))
            for key, (n_el, flat) in acc.items():
                if n_el:
                    gmsh.view.addListData(out_tag, key, n_el, flat)
            gmsh.view.write(out_tag, cfg["out_file"])
            result.update({"ok": True, "ran": True,
                           "out_file": cfg["out_file"],
                           "levels": levels,
                           "pieces_per_level": pieces_per_level,
                           "pieces": {k: v[0] for k, v in acc.items()}})

        elif op == "streamlines_2d":
            # Evenly-spaced streamlines on a plane patch (Jobard-Lefer
            # 1997): trace the IN-PLANE projection of the vector view;
            # new seeds spawn d_sep away from accepted lines and lines
            # stop d_test (= d_sep/2) from other lines -- the uniform-
            # density picture ParaView only offers for native-2D data.
            # Exact field lines where the plane is a symmetry plane
            # (B.n = 0); elsewhere it is the projected-field portrait.
            import collections as _collections
            tag = _view_tag_by_selector(tags, cfg.get("view"))
            step = int(cfg.get("time_step", 0))
            o = cfg["origin"]
            u_vec = [cfg["u_point"][k] - o[k] for k in range(3)]
            v_raw = [cfg["v_point"][k] - o[k] for k in range(3)]
            Lu = sum(c * c for c in u_vec) ** 0.5
            if Lu < 1e-30:
                raise RuntimeError("u_point coincides with origin")
            uh = [c / Lu for c in u_vec]
            dot_uv = sum(a * b for a, b in zip(v_raw, uh))
            v_ortho = [v_raw[k] - dot_uv * uh[k] for k in range(3)]
            Lv = sum(c * c for c in v_ortho) ** 0.5
            if Lv < 1e-30:
                raise RuntimeError("v_point is collinear with u_point")
            vh = [c / Lv for c in v_ortho]
            diag = (Lu * Lu + Lv * Lv) ** 0.5
            d_sep = float(cfg.get("d_sep") or diag / 30.0)
            d_test = float(cfg.get("d_test") or 0.5 * d_sep)
            ds = float(cfg.get("step_size")
                       or max(d_sep / 4.0, diag / 800.0))
            max_steps = int(cfg.get("max_steps", 1500))
            max_lines = int(cfg.get("max_lines", 200))
            budget = int(cfg.get("max_total_steps", 150000))
            min_arc = 1.5 * d_sep

            def _p3(a, b):
                return [o[k] + a * uh[k] + b * vh[k] for k in range(3)]

            def _sample2(a, b):
                if a < 0.0 or a > Lu or b < 0.0 or b > Lv:
                    return None, 0.0
                q = _p3(a, b)
                values, dist = gmsh.view.probe(tag, q[0], q[1], q[2],
                                               step=step)
                if len(values) < 3 or float(dist) > 0.0:
                    return None, 0.0
                vx, vy, vz = (float(values[0]), float(values[1]),
                              float(values[2]))
                fa = vx * uh[0] + vy * uh[1] + vz * uh[2]
                fb = vx * vh[0] + vy * vh[1] + vz * vh[2]
                nip = (fa * fa + fb * fb) ** 0.5
                if nip < 1e-30:
                    return None, 0.0
                n3 = (vx * vx + vy * vy + vz * vz) ** 0.5
                return (fa / nip, fb / nip), n3

            grid = {}

            def _add_to_grid(pts):
                for (a, b, _n) in pts:
                    grid.setdefault((int(a // d_sep), int(b // d_sep)),
                                    []).append((a, b))

            def _too_close(a, b, radius):
                ci, cj = int(a // d_sep), int(b // d_sep)
                r2 = radius * radius
                for i in range(ci - 1, ci + 2):
                    for j in range(cj - 1, cj + 2):
                        for (pa, pb) in grid.get((i, j), ()):
                            if (pa - a) ** 2 + (pb - b) ** 2 < r2:
                                return True
                return False

            steps_used = [0]  # list: the op body runs at module level
            budget_reached = [False]

            def _trace2(seed):
                t0, n0 = _sample2(seed[0], seed[1])
                if t0 is None:
                    return None, False
                closed = False
                halves = []
                for sgn in (1, -1):
                    pth = [(seed[0], seed[1], n0)]
                    a, b = seed
                    arc = 0.0

                    def _f(qa, qb):
                        t, _n = _sample2(qa, qb)
                        if t is None:
                            return None
                        return (sgn * t[0], sgn * t[1])

                    for _ in range(max_steps):
                        if steps_used[0] >= budget:
                            budget_reached[0] = True
                            break
                        steps_used[0] += 1
                        k1 = _f(a, b)
                        if k1 is None:
                            break
                        k2 = _f(a + 0.5 * ds * k1[0],
                                b + 0.5 * ds * k1[1])
                        if k2 is None:
                            break
                        k3 = _f(a + 0.5 * ds * k2[0],
                                b + 0.5 * ds * k2[1])
                        if k3 is None:
                            break
                        k4 = _f(a + ds * k3[0], b + ds * k3[1])
                        if k4 is None:
                            break
                        a2 = a + ds / 6.0 * (k1[0] + 2 * k2[0]
                                             + 2 * k3[0] + k4[0])
                        b2 = b + ds / 6.0 * (k1[1] + 2 * k2[1]
                                             + 2 * k3[1] + k4[1])
                        t2, n2 = _sample2(a2, b2)
                        if t2 is None:
                            break
                        if _too_close(a2, b2, d_test):
                            break
                        a, b = a2, b2
                        arc += ds
                        pth.append((a, b, n2))
                        if arc > 4.0 * d_test:
                            dd = ((a - seed[0]) ** 2
                                  + (b - seed[1]) ** 2) ** 0.5
                            # march-direction alignment: dot(sgn*t2,
                            # sgn*t0) = t2.t0 (sgn cancels)
                            if (dd < max(ds, 0.5 * d_test)
                                    and (t2[0] * t0[0]
                                         + t2[1] * t0[1]) > 0.8):
                                pth.append((seed[0], seed[1], n0))
                                closed = True
                                break
                    halves.append(pth)
                    if budget_reached[0]:
                        if len(halves) == 1:
                            return halves[0], False
                        merged = list(reversed(halves[1]))[:-1] + halves[0]
                        return merged, False
                    if closed:
                        return halves[0], True
                merged = list(reversed(halves[1]))[:-1] + halves[0]
                return merged, False

            # first seed: the strongest in-plane point of a coarse scan
            seed0 = cfg.get("first_seed")
            if seed0 is None:
                best, best_n = None, -1.0
                for i in range(1, 8):
                    for j in range(1, 8):
                        a, b = Lu * i / 8.0, Lv * j / 8.0
                        t, n3 = _sample2(a, b)
                        if t is not None and n3 > best_n:
                            best, best_n = (a, b), n3
                if best is None:
                    raise RuntimeError(
                        "no in-plane field found on the patch (all scan "
                        "points outside the data or |v_inplane| = 0)")
                seed0 = best
            queue = _collections.deque([tuple(seed0)])
            lines = []
            budget_hit = False
            while queue and len(lines) < max_lines:
                if steps_used[0] >= budget:
                    budget_hit = True
                    break
                cand = queue.popleft()
                # 0.95: candidates sit EXACTLY d_sep from their parent
                # line, so testing at the full radius rejects them all
                # on floating-point rounding
                if _too_close(cand[0], cand[1], 0.95 * d_sep):
                    continue
                pts, closed = _trace2(cand)
                if budget_reached[0]:
                    budget_hit = True
                if pts is None:
                    continue
                arc = ds * (len(pts) - 1)
                if arc < min_arc and not closed:
                    continue
                _add_to_grid(pts)
                lines.append((pts, closed))
                stride = max(1, int(d_sep / ds))
                for i in range(0, len(pts), stride):
                    a, b, _n = pts[i]
                    inext = min(i + 1, len(pts) - 1)
                    ta = pts[inext][0] - pts[i - 1][0] if i > 0 else \
                        pts[inext][0] - pts[i][0]
                    tb = pts[inext][1] - pts[i - 1][1] if i > 0 else \
                        pts[inext][1] - pts[i][1]
                    tn = (ta * ta + tb * tb) ** 0.5
                    if tn < 1e-30:
                        continue
                    na, nb = -tb / tn, ta / tn
                    for s in (1.0, -1.0):
                        queue.append((a + s * na * d_sep,
                                      b + s * nb * d_sep))
                if budget_hit:
                    break

            sl_data = []
            n_segments = 0
            for pts, _cl in lines:
                for i in range(len(pts) - 1):
                    pa = _p3(pts[i][0], pts[i][1])
                    pb = _p3(pts[i + 1][0], pts[i + 1][1])
                    sl_data += [pa[0], pb[0], pa[1], pb[1], pa[2], pb[2],
                                pts[i][2], pts[i + 1][2]]
                    n_segments += 1
            out_tag = gmsh.view.add(cfg.get("result_name",
                                            "streamlines_2d"))
            if n_segments:
                gmsh.view.addListData(out_tag, "SL", n_segments, sl_data)
            gmsh.view.write(out_tag, cfg["out_file"])

            arrows_every = int(cfg.get("arrows_every", 0))
            n_arrows = 0
            if arrows_every > 0:
                vp_data = []
                for pts, _cl in lines:
                    for i in range(1, len(pts) - 1, arrows_every):
                        ta = pts[i + 1][0] - pts[i - 1][0]
                        tb = pts[i + 1][1] - pts[i - 1][1]
                        tn = (ta * ta + tb * tb) ** 0.5
                        if tn < 1e-30:
                            continue
                        t3 = [(ta * uh[k] + tb * vh[k]) / tn * pts[i][2]
                              for k in range(3)]
                        vp_data += _p3(pts[i][0], pts[i][1]) + t3
                        n_arrows += 1
                arrow_tag = gmsh.view.add("streamlines_2d_arrows")
                if n_arrows:
                    gmsh.view.addListData(arrow_tag, "VP", n_arrows,
                                          vp_data)
                gmsh.view.write(arrow_tag, cfg["out_file"], append=True)

            result.update({
                "ok": True, "ran": True,
                "out_file": cfg["out_file"],
                "n_lines": len(lines),
                "n_segments": n_segments,
                "n_closed": sum(1 for _p, cl in lines if cl),
                "n_arrows": n_arrows,
                "d_sep": d_sep, "d_test": d_test, "step_size": ds,
                "steps_used": steps_used[0],
                "budget_exceeded": budget_hit,
                "plane": {"origin": o, "u_axis": uh, "v_axis": vh,
                          "extent": [Lu, Lv]},
                "lines": [
                    {"n_points": len(pts), "closed": cl,
                     "points": ([[round(c, 12) for c in _p3(a, b)]
                                 for a, b, _n in pts]
                                if cfg.get("return_points") else None)}
                    for pts, cl in lines],
            })

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


def _finite_float(value: Any, name: str) -> float:
    try:
        clean = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be a finite number") from exc
    if not math.isfinite(clean):
        raise ValueError(f"{name} must be a finite number")
    return clean


def _finite_vector(value: Any, name: str, size: int) -> list[float]:
    try:
        clean = [_finite_float(component, name) for component in value]
    except TypeError as exc:
        raise ValueError(f"{name} needs exactly {size} coordinates") from exc
    if len(clean) != size:
        raise ValueError(f"{name} needs exactly {size} coordinates")
    return clean


def _adaptive_extraction_controls(recur_level: Any,
                                  target_error: Any) -> tuple[int, float]:
    level_value = _finite_float(recur_level, "recur_level")
    if not level_value.is_integer():
        raise ValueError("recur_level must be an integer")
    level = int(level_value)
    if not 0 <= level <= 6:
        raise ValueError("recur_level must be in [0, 6]")
    error = _finite_float(target_error, "target_error")
    if error <= 0.0:
        raise ValueError("target_error must be positive")
    return level, error


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
               recur_level: int = 0, target_error: float = 1e-4,
               out_file: str | Path | None = None,
               timeout_s: float = 300.0) -> dict[str, Any]:
    """Extract the isosurface of a scalar view and save it.

    ``0 < recur_level <= 6`` enables ADAPTIVE extraction on high-order data
    (order-2 NodeData from GmshPostExport): each element is recursively
    subdivided using the actual high-order interpolant, so the surface
    follows the curved field instead of the P1 chord (measured on a
    TET10 quadratic: radial error 0.21 -> 0.008 at level 4).  Requires
    setting the view adaptive BEFORE the plugin runs -- handled here.
    """
    err = _check_input(path)
    if err:
        return err
    try:
        clean_value = _finite_float(value, "value")
        clean_level, clean_error = _adaptive_extraction_controls(
            recur_level, target_error)
    except ValueError as exc:
        return {"ok": False, "error": str(exc)}
    return _run_post({"op": "isosurface", "path": str(path), "view": view,
                      "value": clean_value,
                      "recur_level": clean_level,
                      "target_error": clean_error,
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
                adaptive: bool = True, closure: bool = True,
                max_turn_deg: float = 10.0, arrows_every: int = 0,
                return_points: bool = False,
                out_file: str | Path | None = None,
                timeout_s: float = 600.0) -> dict[str, Any]:
    """Trace field lines of a vector view from seeds on a line segment.

    Probe-driven arc-length RK4 with curvature ADAPTIVITY (step_size =
    bbox diagonal / 200 by default and is the MAXIMUM step; it halves
    wherever the turn per step exceeds max_turn_deg) and CLOSED-LOOP
    detection -- a magnetic field line that returns to its seed is
    closed exactly and reported as such instead of overdrawing or
    stopping mid-loop.  One merged polyline per seed (backward +
    forward), local |v| as the line color, per-line termination
    reasons (left_data | closed | stagnation | max_steps), and
    optional direction arrows every arrows_every-th point as a
    companion VP view.  return_points=True returns the coordinates.
    """
    err = _check_input(path)
    if err:
        return err
    try:
        s0 = _finite_vector(seed_start, "seed_start", 3)
        s1 = _finite_vector(seed_end, "seed_end", 3)
        clean_step = (None if step_size is None else
                      _finite_float(step_size, "step_size"))
        clean_turn = _finite_float(max_turn_deg, "max_turn_deg")
    except ValueError as exc:
        return {"ok": False, "error": str(exc)}
    if clean_step is not None and clean_step <= 0.0:
        return {"ok": False, "error": "step_size must be positive"}
    if int(n_seeds) < 1:
        return {"ok": False, "error": "n_seeds must be positive"}
    if int(max_steps) < 1:
        return {"ok": False, "error": "max_steps must be positive"}
    if not 0.0 < clean_turn <= 180.0:
        return {"ok": False,
                "error": "max_turn_deg must be in (0, 180]"}
    if int(arrows_every) < 0:
        return {"ok": False, "error": "arrows_every must be non-negative"}
    if int(time_step) < 0:
        return {"ok": False, "error": "time_step must be non-negative"}
    return _run_post({"op": "streamlines", "path": str(path), "view": view,
                      "seed_start": s0, "seed_end": s1,
                      "n_seeds": int(n_seeds),
                      "step_size": clean_step,
                      "max_steps": int(max_steps),
                      "both_directions": bool(both_directions),
                      "time_step": int(time_step),
                      "adaptive": bool(adaptive),
                      "closure": bool(closure),
                      "max_turn_deg": clean_turn,
                      "arrows_every": int(arrows_every),
                      "return_points": bool(return_points),
                      "out_file": _default_out(path, out_file, "stream")},
                     timeout_s)


_PARTICLE_SPECIES = {
    # name: (charge in elementary charges, mass in kg)
    "electron": (-1.0, 9.1093837015e-31),
    "positron": (1.0, 9.1093837015e-31),
    "proton": (1.0, 1.67262192369e-27),
    "antiproton": (-1.0, 1.67262192369e-27),
    "alpha": (2.0, 6.6446573357e-27),
}
_ELEMENTARY_CHARGE_C = 1.602176634e-19
_AMU_KG = 1.66053906660e-27


def particle_trace(path: str | Path, seeds: list[list[float]],
                   direction: list[float], kinetic_energy_ev: float, *,
                   species: str = "electron",
                   charge_e: float | None = None,
                   mass_amu: float | None = None,
                   view: str | int | None = None,
                   e_view: str | int | None = None,
                   time_step: int = 0,
                   dt_s: float | None = None,
                   steps_per_gyration: int = 64,
                   max_steps: int = 20000,
                   max_time_s: float | None = None,
                   color_by: str = "time",
                   arrows_every: int = 0,
                   animation_frames: int = 0,
                   animation_mode: str = "trail",
                   comet_window: float = 0.15,
                   return_points: bool = False,
                   out_file: str | Path | None = None,
                   timeout_s: float = 600.0) -> dict[str, Any]:
    """Trace charged-particle ORBITS through the probed B (and E) views.

    Relativistic Boris pusher on ``dp/dt = q(E + v x B)`` -- the
    particle-dynamics companion to ``streamlines`` (which draws the
    massless tangent curves of the field itself; a particle GYRATES
    around those lines instead of following them).  The B view is in
    Tesla on a mesh in meters; the optional ``e_view`` adds an electric
    field in V/m from the same file.

    Each seed launches one particle of the given ``species`` (or a
    custom ``charge_e``/``mass_amu`` pair) with the given kinetic
    energy along ``direction``.  The time step defaults to 1 /
    ``steps_per_gyration`` of the LOCAL gyration period at the seed;
    where B vanishes at a seed an explicit ``dt_s`` is required.  A
    trace ends on ``left_data`` (the orbit exits the field data),
    ``max_steps``, or ``max_time``.

    Output: a ``particle_tracks`` SL view colored by ``color_by``
    (``time`` | ``speed`` | ``energy``), optional velocity arrows every
    ``arrows_every``-th sample, and per-track diagnostics including the
    seed gyroradius and ``speed_change_rel`` -- in a pure magnetic
    field the Boris rotation conserves speed exactly, so a nonzero
    value there measures integrator health, not physics.

    ``animation_frames > 0`` adds a multi-step ``beam (...)`` view: the
    particles FLY along their orbits on a SHARED clock (absolute flight
    time, so seeds with different gyration periods stay in step), ready
    for ``gmsh_export_animation``.  ``animation_mode="trail"`` leaves
    the flown path behind; ``"comet"`` lights only a moving window of
    ``comet_window`` x the flight time.  The beam view renders
    correctly ONLY with the fixed colour range returned in
    ``animation["render_hint"]`` -- future segments are hidden by an
    out-of-range sentinel value, which an autoscaled colour bar would
    fold back into the picture.
    """
    err = _check_input(path)
    if err:
        return err
    try:
        clean_seeds = [_finite_vector(s, "seeds", 3) for s in seeds]
        clean_dir = _finite_vector(direction, "direction", 3)
        clean_ke = _finite_float(kinetic_energy_ev, "kinetic_energy_ev")
        clean_dt = None if dt_s is None else _finite_float(dt_s, "dt_s")
        clean_tmax = (None if max_time_s is None
                      else _finite_float(max_time_s, "max_time_s"))
    except (TypeError, ValueError) as exc:
        return {"ok": False, "error": str(exc)}
    if not clean_seeds:
        return {"ok": False, "error": "seeds must contain at least one point"}
    if not any(clean_dir):
        return {"ok": False, "error": "direction must be nonzero"}
    if clean_ke <= 0.0:
        return {"ok": False, "error": "kinetic_energy_ev must be positive"}
    if (charge_e is None) != (mass_amu is None):
        return {"ok": False, "error":
                "custom particles need BOTH charge_e and mass_amu"}
    if charge_e is not None:
        try:
            q_c = _finite_float(charge_e, "charge_e") * _ELEMENTARY_CHARGE_C
            m_kg = _finite_float(mass_amu, "mass_amu") * _AMU_KG
        except ValueError as exc:
            return {"ok": False, "error": str(exc)}
        if q_c == 0.0:
            return {"ok": False, "error": "charge_e must be nonzero"}
        if m_kg <= 0.0:
            return {"ok": False, "error": "mass_amu must be positive"}
    else:
        if species not in _PARTICLE_SPECIES:
            return {"ok": False, "error":
                    f"unknown species {species!r}; available: "
                    f"{sorted(_PARTICLE_SPECIES)} (or give charge_e + "
                    "mass_amu)"}
        q_e, m_kg = _PARTICLE_SPECIES[species]
        q_c = q_e * _ELEMENTARY_CHARGE_C
    if clean_dt is not None and clean_dt <= 0.0:
        return {"ok": False, "error": "dt_s must be positive"}
    if clean_tmax is not None and clean_tmax <= 0.0:
        return {"ok": False, "error": "max_time_s must be positive"}
    if int(steps_per_gyration) < 4:
        return {"ok": False, "error": "steps_per_gyration must be >= 4"}
    if int(max_steps) < 1:
        return {"ok": False, "error": "max_steps must be positive"}
    if color_by not in ("time", "speed", "energy"):
        return {"ok": False, "error":
                "color_by must be one of: time, speed, energy"}
    if int(arrows_every) < 0:
        return {"ok": False, "error": "arrows_every must be non-negative"}
    if int(time_step) < 0:
        return {"ok": False, "error": "time_step must be non-negative"}
    if int(animation_frames) < 0:
        return {"ok": False,
                "error": "animation_frames must be non-negative"}
    if animation_mode not in ("trail", "comet"):
        return {"ok": False,
                "error": "animation_mode must be 'trail' or 'comet'"}
    try:
        clean_window = _finite_float(comet_window, "comet_window")
    except ValueError as exc:
        return {"ok": False, "error": str(exc)}
    if not 0.0 < clean_window <= 1.0:
        return {"ok": False, "error": "comet_window must be in (0, 1]"}
    return _run_post({"op": "particle_trace", "path": str(path),
                      "view": view, "e_view": e_view,
                      "seeds": clean_seeds, "direction": clean_dir,
                      "kinetic_energy_ev": clean_ke,
                      "charge_c": q_c, "mass_kg": m_kg,
                      "time_step": int(time_step),
                      "dt_s": clean_dt,
                      "steps_per_gyration": int(steps_per_gyration),
                      "max_steps": int(max_steps),
                      "max_time_s": clean_tmax,
                      "color_by": color_by,
                      "arrows_every": int(arrows_every),
                      "animation_frames": int(animation_frames),
                      "animation_mode": animation_mode,
                      "comet_window": clean_window,
                      "return_points": bool(return_points),
                      "out_file": _default_out(path, out_file, "tracks")},
                     timeout_s)


def flux_lines(path: str | Path, *, n_levels: int = 20,
               levels: list[float] | None = None,
               view: str | int | None = None,
               recur_level: int = 0, target_error: float = 1e-4,
               result_name: str = "flux_lines",
               out_file: str | Path | None = None,
               timeout_s: float = 600.0) -> dict[str, Any]:
    """Equally spaced isolines of a scalar view, merged into ONE view.

    THE 2D-magnetics field-line tool: on a planar A_z view (or the
    axisymmetric flux function psi = r A_theta) the equally spaced
    levels ARE the exact field lines with EQUAL FLUX between adjacent
    lines -- the FEMM-style motor plot, no integration, no seeds, no
    density tuning.  Levels default to n_levels interior values
    between the view min and max; on a 3D scalar the same call yields
    a multi-level isosurface stack.  ``recur_level > 0`` extracts
    adaptively on high-order data (see ``isosurface``).
    """
    err = _check_input(path)
    if err:
        return err
    try:
        clean_level, clean_error = _adaptive_extraction_controls(
            recur_level, target_error)
    except ValueError as exc:
        return {"ok": False, "error": str(exc)}
    if levels is not None:
        try:
            levels = [_finite_float(lv, "levels") for lv in levels]
        except (TypeError, ValueError) as exc:
            return {"ok": False, "error": str(exc)}
        if not levels:
            return {"ok": False, "error": "levels must not be empty"}
        if len(set(levels)) != len(levels):
            return {"ok": False, "error": "levels must be unique"}
    elif int(n_levels) < 1:
        return {"ok": False, "error": "n_levels must be positive"}
    if not str(result_name).strip():
        return {"ok": False, "error": "result_name must not be empty"}
    return _run_post({"op": "flux_lines", "path": str(path), "view": view,
                      "n_levels": int(n_levels), "levels": levels,
                      "recur_level": clean_level,
                      "target_error": clean_error,
                      "result_name": str(result_name),
                      "out_file": _default_out(path, out_file, "flux")},
                     timeout_s)


def streamlines_2d(path: str | Path, origin: list[float],
                   u_point: list[float], v_point: list[float], *,
                   d_sep: float | None = None,
                   d_test: float | None = None,
                   step_size: float | None = None,
                   view: str | int | None = None,
                   first_seed: list[float] | None = None,
                   max_lines: int = 200, max_steps: int = 1500,
                   max_total_steps: int = 150000,
                   arrows_every: int = 0, time_step: int = 0,
                   result_name: str = "streamlines_2d",
                   return_points: bool = False,
                   out_file: str | Path | None = None,
                   timeout_s: float = 900.0) -> dict[str, Any]:
    """Evenly spaced streamlines on a plane patch (Jobard-Lefer).

    The uniform-density streamline picture ParaView only offers for
    native-2D datasets, here on ANY plane slice of a 3D field: the
    IN-PLANE projection of the vector view is traced, new seeds spawn
    d_sep away from accepted lines, and lines stop d_test (= d_sep/2)
    from other lines -- no manual seed placement, no bunching, closed
    loops handled.  The patch follows the resample_grid convention
    (origin + u/v edge endpoints; the v axis is orthogonalized).
    Exact field lines where the plane is a symmetry plane (B.n = 0);
    elsewhere it is the standard projected-field portrait.
    """
    err = _check_input(path)
    if err:
        return err
    pts = {"origin": origin, "u_point": u_point, "v_point": v_point}
    clean: dict[str, list[float]] = {}
    try:
        for key, val in pts.items():
            clean[key] = _finite_vector(val, key, 3)
        clean_sep = (None if d_sep is None else
                     _finite_float(d_sep, "d_sep"))
        clean_test = (None if d_test is None else
                      _finite_float(d_test, "d_test"))
        clean_step = (None if step_size is None else
                      _finite_float(step_size, "step_size"))
        if first_seed is not None:
            first_seed = _finite_vector(first_seed, "first_seed", 2)
    except ValueError as exc:
        return {"ok": False, "error": str(exc)}

    u_vec = [clean["u_point"][k] - clean["origin"][k] for k in range(3)]
    v_raw = [clean["v_point"][k] - clean["origin"][k] for k in range(3)]
    lu = math.sqrt(sum(component * component for component in u_vec))
    if lu < 1e-30:
        return {"ok": False, "error": "u_point must differ from origin"}
    u_hat = [component / lu for component in u_vec]
    dot_uv = sum(a * b for a, b in zip(v_raw, u_hat))
    v_ortho = [v_raw[k] - dot_uv * u_hat[k] for k in range(3)]
    lv = math.sqrt(sum(component * component for component in v_ortho))
    if lv < 1e-30:
        return {"ok": False,
                "error": "v_point must not be collinear with the U edge"}
    diag = math.hypot(lu, lv)
    effective_sep = clean_sep if clean_sep is not None else diag / 30.0
    effective_test = (clean_test if clean_test is not None else
                      0.5 * effective_sep)
    effective_step = (clean_step if clean_step is not None else
                      max(effective_sep / 4.0, diag / 800.0))
    for name, value in (("d_sep", effective_sep),
                        ("d_test", effective_test),
                        ("step_size", effective_step)):
        if value <= 0.0:
            return {"ok": False, "error": f"{name} must be positive"}
    if effective_test > effective_sep:
        return {"ok": False, "error": "d_test must not exceed d_sep"}
    if effective_step > effective_test:
        return {"ok": False,
                "error": "step_size must not exceed d_test"}
    if int(max_lines) < 1:
        return {"ok": False, "error": "max_lines must be positive"}
    if int(max_steps) < 1:
        return {"ok": False, "error": "max_steps must be positive"}
    if int(max_total_steps) < 1:
        return {"ok": False, "error": "max_total_steps must be positive"}
    if int(arrows_every) < 0:
        return {"ok": False, "error": "arrows_every must be non-negative"}
    if int(time_step) < 0:
        return {"ok": False, "error": "time_step must be non-negative"}
    if not str(result_name).strip():
        return {"ok": False, "error": "result_name must not be empty"}
    if first_seed is not None and not (
            0.0 <= first_seed[0] <= lu and 0.0 <= first_seed[1] <= lv):
        return {"ok": False,
                "error": "first_seed must lie inside the [0, Lu] x "
                         "[0, Lv] plane patch"}
    return _run_post({"op": "streamlines_2d", "path": str(path),
                      "view": view, **clean,
                      "d_sep": effective_sep, "d_test": effective_test,
                      "step_size": effective_step,
                      "first_seed": first_seed,
                      "max_lines": int(max_lines),
                      "max_steps": int(max_steps),
                      "max_total_steps": int(max_total_steps),
                      "arrows_every": int(arrows_every),
                      "time_step": int(time_step),
                      "result_name": str(result_name),
                      "return_points": bool(return_points),
                      "out_file": _default_out(path, out_file, "stream2d")},
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
                            ylabel="value")
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
                 xlabel: str = "", ylabel: str = "",
                 timeout_s: float = 120.0) -> dict[str, Any]:
    out = Path(png)
    out.parent.mkdir(parents=True, exist_ok=True)
    cfg = {"kind": "line", "series": series, "xlabel": xlabel,
           "ylabel": ylabel, "png": str(out)}
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
    finite = [s for s in samples if math.isfinite(s)]
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
               "xlabel": name, "ylabel": "count", "png": str(out)}
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
                            ylabel=out["view"])
        out["plot_png"] = plot.get("png")
        if not plot.get("ok"):
            out["plot_error"] = plot.get("error")
    return out


# =====================================================================
# Cross-file colour range, compound selection, flow texture
# =====================================================================

def _row_samples(view: dict[str, Any], vals: list[float],
                 component: int | None) -> list[float]:
    """Return the scalar samples represented by one data row."""
    ncomp = int(view["components"])
    if component is not None:
        if isinstance(component, bool) or component < 0 or component >= ncomp:
            raise IndexError(
                f"view {view['name']!r} has {ncomp} components, "
                f"component {component} requested")
        samples = [float(v) for v in vals[component::ncomp]]
    elif ncomp == 1:
        samples = [float(v) for v in vals]
    else:
        samples = [
            math.sqrt(sum(float(v) * float(v)
                          for v in vals[i:i + ncomp]))
            for i in range(0, len(vals), ncomp)]
    if any(not math.isfinite(value) for value in samples):
        raise ValueError(f"view {view['name']!r} contains NaN or infinity")
    return samples


def _view_values(view: dict[str, Any],
                 component: int | None) -> dict[int, float]:
    """Reduce a view to one representative value per node/element tag.

    ``ElementNodeData`` has one sample per local element node.  Selection
    and file-series statistics operate per element, so those local samples
    are averaged instead of being mistaken for one large vector.
    """
    out: dict[int, float] = {}
    for tag, vals in view["rows"].items():
        samples = _row_samples(view, vals, component)
        if samples:
            out[tag] = sum(samples) / len(samples)
    return out


def field_range(paths: list[str | Path], *,
                view: str | int | None = None,
                component: int | None = None,
                time_step: int | None = None) -> dict[str, Any]:
    """Union colour range across several .msh files (and their views).

    gmsh autoscales EVERY view to its own extrema, so two renders of the
    same quantity are not comparable until they are pinned to one range.
    This computes that range without launching gmsh (pure-Python
    reader), so the result feeds straight into
    ``gmsh_render(color={"range": [...]})``.

    Args:
        paths: .msh files to combine.
        view: view name or index to restrict to (default: every view).
        component: component index, or None for the magnitude.
        time_step: restrict to one step (default: every step).

    Returns:
        ``{"range": [lo, hi], "per_file": {...}, "views": [...]}``.
    """
    from .msh_inspect import read_msh_data

    if not paths:
        raise ValueError("field_range needs at least one path")
    lo, hi = math.inf, -math.inf
    per_file: dict[str, Any] = {}
    seen_views: list[str] = []
    for p in paths:
        src = Path(p)
        if not src.is_file():
            return {"ok": False, "error": f"file not found: {src}"}
        data = read_msh_data(src)
        f_lo, f_hi, used = math.inf, -math.inf, []
        for idx, v in enumerate(data["views"]):
            if isinstance(view, bool):
                raise TypeError("view must be a name or an index")
            if isinstance(view, int) and idx != view:
                continue
            if isinstance(view, str) and v["name"] != view:
                continue
            if time_step is not None and v.get("step") != time_step:
                continue
            vals = [sample for row in v["rows"].values()
                    for sample in _row_samples(v, row, component)]
            if not vals:
                continue
            f_lo = min(f_lo, min(vals))
            f_hi = max(f_hi, max(vals))
            used.append(v["name"])
            if v["name"] not in seen_views:
                seen_views.append(v["name"])
        if not used:
            return {"ok": False,
                    "error": f"no matching view in {src.name} "
                             f"(views: {[v['name'] for v in data['views']]})"}
        per_file[str(src)] = {"range": [f_lo, f_hi], "views": used}
        lo, hi = min(lo, f_lo), max(hi, f_hi)
    return {"ok": True, "range": [lo, hi], "per_file": per_file,
            "views": seen_views, "component": component,
            "n_files": len(paths)}


_SELECT_ALLOWED = {
    "abs": abs, "min": min, "max": max, "sqrt": math.sqrt,
    "log": math.log, "log10": math.log10, "exp": math.exp,
    "sin": math.sin, "cos": math.cos, "atan2": math.atan2,
    "hypot": math.hypot, "pi": math.pi, "e": math.e,
}


class _SelectExpressionValidator(ast.NodeVisitor):
    """Allow arithmetic/boolean expressions without Python object access."""

    _binary = (ast.Add, ast.Sub, ast.Mult, ast.Div, ast.FloorDiv,
               ast.Mod, ast.Pow)
    _unary = (ast.UAdd, ast.USub, ast.Not)
    _boolean = (ast.And, ast.Or)
    _compare = (ast.Eq, ast.NotEq, ast.Lt, ast.LtE, ast.Gt, ast.GtE)

    def __init__(self, names: set[str]):
        self.names = names

    def generic_visit(self, node):
        raise ValueError(
            f"unsupported expression syntax: {type(node).__name__}")

    def visit_Expression(self, node):
        self.visit(node.body)

    def visit_Constant(self, node):
        if not isinstance(node.value, (bool, int, float)):
            raise ValueError("only numeric and boolean constants are allowed")

    def visit_Name(self, node):
        if node.id not in self.names:
            raise NameError(f"name {node.id!r} is not defined")

    def visit_Call(self, node):
        if not isinstance(node.func, ast.Name) \
                or node.func.id not in _SELECT_ALLOWED \
                or not callable(_SELECT_ALLOWED[node.func.id]):
            raise ValueError("only documented math functions may be called")
        if node.keywords:
            raise ValueError("keyword arguments are not allowed")
        for arg in node.args:
            self.visit(arg)

    def visit_BinOp(self, node):
        if not isinstance(node.op, self._binary):
            raise ValueError(
                f"unsupported arithmetic operator: {type(node.op).__name__}")
        self.visit(node.left)
        self.visit(node.right)

    def visit_UnaryOp(self, node):
        if not isinstance(node.op, self._unary):
            raise ValueError(
                f"unsupported unary operator: {type(node.op).__name__}")
        self.visit(node.operand)

    def visit_BoolOp(self, node):
        if not isinstance(node.op, self._boolean):
            raise ValueError(
                f"unsupported boolean operator: {type(node.op).__name__}")
        for value in node.values:
            self.visit(value)

    def visit_Compare(self, node):
        if any(not isinstance(op, self._compare) for op in node.ops):
            raise ValueError("unsupported comparison operator")
        self.visit(node.left)
        for comparator in node.comparators:
            self.visit(comparator)


def _parse_select_expression(expression: str,
                             available: list[str]) -> ast.expr:
    try:
        tree = ast.parse(expression, mode="eval")
        _SelectExpressionValidator(
            set(available) | set(_SELECT_ALLOWED) | {"True", "False"}
        ).visit(tree)
    except (SyntaxError, ValueError, NameError) as exc:
        raise ValueError(str(exc)) from exc
    return tree.body


def _eval_select_expression(node: ast.expr,
                            env: dict[str, Any]) -> Any:
    """Evaluate an already validated selection AST without Python eval."""
    if isinstance(node, ast.Constant):
        return node.value if isinstance(node.value, bool) else float(node.value)
    if isinstance(node, ast.Name):
        return env[node.id]
    if isinstance(node, ast.Call):
        func = env[node.func.id]
        return func(*[_eval_select_expression(arg, env) for arg in node.args])
    if isinstance(node, ast.BinOp):
        left = _eval_select_expression(node.left, env)
        right = _eval_select_expression(node.right, env)
        if isinstance(node.op, ast.Add):
            return left + right
        if isinstance(node.op, ast.Sub):
            return left - right
        if isinstance(node.op, ast.Mult):
            return left * right
        if isinstance(node.op, ast.Div):
            return left / right
        if isinstance(node.op, ast.FloorDiv):
            return left // right
        if isinstance(node.op, ast.Mod):
            return left % right
        if isinstance(node.op, ast.Pow):
            return left ** right
    if isinstance(node, ast.UnaryOp):
        value = _eval_select_expression(node.operand, env)
        if isinstance(node.op, ast.UAdd):
            return +value
        if isinstance(node.op, ast.USub):
            return -value
        if isinstance(node.op, ast.Not):
            return not value
    if isinstance(node, ast.BoolOp):
        if isinstance(node.op, ast.And):
            for value in node.values:
                result = _eval_select_expression(value, env)
                if not result:
                    return result
            return result
        for value in node.values:
            result = _eval_select_expression(value, env)
            if result:
                return result
        return result
    if isinstance(node, ast.Compare):
        left = _eval_select_expression(node.left, env)
        for op, comparator in zip(node.ops, node.comparators):
            right = _eval_select_expression(comparator, env)
            if isinstance(op, ast.Eq):
                ok = left == right
            elif isinstance(op, ast.NotEq):
                ok = left != right
            elif isinstance(op, ast.Lt):
                ok = left < right
            elif isinstance(op, ast.LtE):
                ok = left <= right
            elif isinstance(op, ast.Gt):
                ok = left > right
            else:
                ok = left >= right
            if not ok:
                return False
            left = right
        return True
    raise ValueError(f"unsupported expression syntax: {type(node).__name__}")


def _select_ident(name: str) -> str:
    out = "".join(c if c.isalnum() else "_" for c in name).strip("_")
    if not out:
        return "view"
    out = out.lower()
    return ("v_" + out if out[:1].isdigit() or not out.isidentifier()
            or keyword.iskeyword(out) else out)


def _select_idents(names: list[str]) -> list[str]:
    reserved = set(_SELECT_ALLOWED) | {"x", "y", "z"}
    reserved.update(f"v{i}" for i in range(len(names)))
    used = set(reserved)
    result = []
    for name in names:
        base = _select_ident(name)
        candidate = base
        suffix = 1
        while candidate in used:
            candidate = (f"{base}_view" if suffix == 1
                         else f"{base}_view{suffix}")
            suffix += 1
        used.add(candidate)
        result.append(candidate)
    return result


def select(path: str | Path, expression: str, *,
           out_file: str | Path | None = None,
           result_name: str = "selection",
           extract: bool = True,
           carry: str | int | None = 0,
           timeout_s: float = 300.0) -> dict[str, Any]:
    """Select elements with a compound condition (ParaView "Find Data").

    ``threshold`` filters ONE view; a real query mixes fields with each
    other and with position -- "where |B| exceeds 1.5 T on the upper
    half".  The expression is evaluated per element in Python, so it can
    reference every view at once, and the surviving elements are written
    as a 1/0 mask view.  ``extract=True`` then runs the threshold
    (Plugin(ExtractElements), MinVal 0.5) on that mask to materialise
    the selection as its own view.

    Names available in the expression:
        ``x``, ``y``, ``z``  element centroid coordinates
        ``v0``, ``v1``, ...  per-view value (magnitude for vectors)
        the view's own name, lowercased with non-word characters turned
        into ``_`` (``B`` -> ``b``, ``|J| [A/m^2]`` -> ``j_a_m_2``)
        plus abs/min/max/sqrt/log/log10/exp/sin/cos/atan2/hypot/pi/e.

    Python boolean operators apply (``and``, ``or``, ``not``).  An
    unknown name raises with the available list rather than evaluating
    to something silently wrong.
    """
    from .msh_inspect import read_msh_data

    src = Path(path)
    if not src.is_file():
        return {"ok": False, "error": f"file not found: {src}"}
    data = read_msh_data(src, include_elements=True)
    if not data["elements"]:
        return {"ok": False, "error": f"{src.name} holds no elements"}
    nodes = data["nodes"]

    names: list[str] = []
    per_view: list[dict[int, float]] = []
    for v in data["views"]:
        vals = _view_values(v, None)
        if v["section"] == "NodeData":
            elem_vals: dict[int, float] = {}
            for tag, el in data["elements"].items():
                got = [vals[n] for n in el["nodes"] if n in vals]
                if got:
                    elem_vals[tag] = sum(got) / len(got)
            per_view.append(elem_vals)
        else:
            per_view.append(vals)
        names.append(v["name"])

    idents = _select_idents(names)
    available = (["x", "y", "z"] + [f"v{i}" for i in range(len(names))]
                 + idents)
    try:
        expression_tree = _parse_select_expression(expression, available)
    except ValueError as exc:
        return {"ok": False, "error": str(exc),
                "available_names": available, "view_names": names}
    selected: list[int] = []
    for tag, el in data["elements"].items():
        pts = [nodes[n] for n in el["nodes"] if n in nodes]
        if not pts:
            continue
        env = dict(_SELECT_ALLOWED)
        env["x"] = sum(p[0] for p in pts) / len(pts)
        env["y"] = sum(p[1] for p in pts) / len(pts)
        env["z"] = sum(p[2] for p in pts) / len(pts)
        for k, (ident, vals) in enumerate(zip(idents, per_view)):
            val = vals.get(tag, 0.0)
            env[f"v{k}"] = val
            env[ident] = val
        try:
            keep = bool(_eval_select_expression(expression_tree, env))
        except Exception as exc:                        # noqa: BLE001
            return {"ok": False,
                    "error": f"expression failed on element {tag}: {exc}",
                    "available_names": available,
                    "view_names": names}
        if keep:
            selected.append(tag)

    out = (Path(out_file) if out_file is not None
           else src.with_name(f"{src.stem}_{result_name}.msh"))
    out.parent.mkdir(parents=True, exist_ok=True)
    chosen = set(selected)
    mask_rows = "".join(f"{tag} {1.0 if tag in chosen else 0.0:.1f}\n"
                        for tag in data["elements"])
    mask = ("$ElementData\n"
            f'1\n"{result_name}"\n1\n0\n'
            f"3\n0\n1\n{len(data['elements'])}\n"
            f"{mask_rows}$EndElementData\n")
    # A second view carries the CHOSEN FIELD on the selected elements and
    # a far-below sentinel elsewhere, so extracting it by value yields a
    # selection that still holds the physics (extracting the 1/0 mask
    # alone gives a flat blob whose colour bar reads "1").
    carry_idx = None
    if carry is not None and names:
        if isinstance(carry, int):
            carry_idx = carry if 0 <= carry < len(names) else None
        else:
            carry_idx = names.index(carry) if carry in names else None
        if carry_idx is None:
            return {"ok": False,
                    "error": f"unknown carry view {carry!r} "
                             f"(views: {names})",
                    "view_names": names}
    carried = ""
    carry_lo = carry_hi = None
    if carry_idx is not None and selected:
        vals = per_view[carry_idx]
        picked = [vals.get(t, 0.0) for t in selected]
        carry_lo, carry_hi = min(picked), max(picked)
        sentinel = carry_lo - 1.0 - abs(carry_lo)
        rows = "".join(
            f"{tag} {(vals.get(tag, 0.0) if tag in chosen else sentinel):.9e}\n"
            for tag in data["elements"])
        carried = ("$ElementData\n"
                   f'1\n"{result_name}_{_select_ident(names[carry_idx])}"\n'
                   f"1\n0\n3\n0\n1\n{len(data['elements'])}\n"
                   f"{rows}$EndElementData\n")
    out.write_text(src.read_text(encoding="utf-8", errors="replace")
                   + mask + carried, encoding="utf-8")

    result: dict[str, Any] = {
        "ok": True, "msh": str(out), "expression": expression,
        "n_elements": len(data["elements"]),
        "n_selected": len(selected),
        "fraction": (len(selected) / len(data["elements"])
                     if data["elements"] else 0.0),
        "mask_view": result_name,
        "available_names": available, "view_names": names}
    if carried:
        result["carried_view"] = names[carry_idx]
        result["carried_range"] = [carry_lo, carry_hi]
    if extract and selected:
        if carried:
            span = max(abs(carry_hi - carry_lo), 1e-30)
            result["extract"] = threshold(
                out, carry_lo - 1e-6 * span, carry_hi + 1e-6 * span,
                view=f"{result_name}_{_select_ident(names[carry_idx])}",
                out_file=str(out.with_name(out.stem + "_extract.pos")),
                timeout_s=timeout_s)
        else:
            result["extract"] = threshold(
                out, 0.5, 1.5, view=result_name,
                out_file=str(out.with_name(out.stem + "_extract.pos")),
                timeout_s=timeout_s)
    return result


def flow_texture(path: str | Path, *,
                 view: str | int | None = None,
                 plane: str = "xy", offset: float = 0.0,
                 density: float = 60.0,
                 out_file: str | Path | None = None,
                 timeout_s: float = 900.0) -> dict[str, Any]:
    """Dense evenly-spaced streamline texture -- the LIC alternative.

    gmsh has no line integral convolution: it cannot smear a noise
    texture along a vector field.  The classical substitute is this --
    Jobard-Lefer evenly spaced streamlines packed densely enough that
    the eye reads them as a flow texture rather than as countable
    curves.  ``density`` is how many line spacings fit across the
    domain diagonal (60 reads as a texture; 15-20 stays countable).

    This is NOT LIC and does not claim to be.  It is strictly better in
    one respect and worse in another: every curve here is a real
    trajectory (quantitative, probe-able), whereas LIC is a purely
    visual convolution that fills every pixel.
    """
    from .msh_inspect import read_msh_data

    src = Path(path)
    if not src.is_file():
        return {"ok": False, "error": f"file not found: {src}"}
    if src.suffix.lower() != ".msh":
        return {"ok": False,
                "error": "flow_texture requires an ASCII MSH v4.x file"}
    density = float(density)
    if not math.isfinite(density) or density <= 0:
        raise ValueError(f"density must be > 0, got {density}")
    axes = {"xy": (0, 1, 2), "yz": (1, 2, 0), "xz": (0, 2, 1)}
    if plane not in axes:
        raise ValueError(
            f"unknown plane {plane!r} (available: {', '.join(sorted(axes))})")
    data = read_msh_data(src)
    if not data["views"]:
        return {"ok": False, "error": f"{src.name} carries no view"}
    if view is None:
        chosen = data["views"][0]
    elif isinstance(view, bool):
        return {"ok": False, "error": "view must be a name or an index"}
    elif isinstance(view, int):
        if view < 0 or view >= len(data["views"]):
            return {"ok": False, "error": f"view index {view} out of range"}
        chosen = data["views"][view]
    else:
        matches = [entry for entry in data["views"]
                   if entry["name"] == view]
        if len(matches) != 1:
            return {"ok": False,
                    "error": f"expected one view {view!r}, got {len(matches)}"}
        chosen = matches[0]
    if chosen["components"] != 3:
        return {"ok": False,
                "error": f"flow_texture needs a 3-component vector view; "
                         f"{chosen['name']!r} has {chosen['components']}"}
    pts = list(data["nodes"].values())
    if not pts:
        return {"ok": False, "error": f"{src.name} holds no nodes"}
    lo = [min(p[i] for p in pts) for i in range(3)]
    hi = [max(p[i] for p in pts) for i in range(3)]
    a0, a1, an = axes[plane]
    span = [hi[i] - lo[i] for i in range(3)]
    diag = math.sqrt(span[a0] ** 2 + span[a1] ** 2)
    if diag <= 0.0:
        return {"ok": False,
                "error": f"{src.name} is degenerate in the {plane} plane"}
    clean_offset = float(offset)
    if not math.isfinite(clean_offset):
        raise ValueError("offset must be finite")
    mid = 0.5 * (lo[an] + hi[an]) + clean_offset

    def _pt(u, v):
        p = [0.0, 0.0, 0.0]
        p[a0], p[a1], p[an] = u, v, mid
        return p

    origin = _pt(lo[a0], lo[a1])
    u_point = _pt(hi[a0], lo[a1])
    v_point = _pt(lo[a0], hi[a1])
    d_sep = diag / float(density)
    res = streamlines_2d(src, origin, u_point, v_point, view=view,
                         d_sep=d_sep, max_lines=max(1, int(4 * density)),
                         out_file=out_file, timeout_s=timeout_s)
    if isinstance(res, dict):
        res["d_sep"] = d_sep
        res["density"] = float(density)
        res["method"] = ("evenly spaced streamlines (Jobard-Lefer); "
                         "NOT line integral convolution")
    return res


_TIME_STATS = ("min", "max", "mean", "std", "rms", "ptp", "argmax_time",
               "argmin_time")


def time_series(paths: list[str | Path], *,
                view: str | int | None = None,
                component: int | None = None,
                times: list[float] | None = None,
                stats: tuple[str, ...] | list[str] = _TIME_STATS,
                out_file: str | Path | None = None,
                points: list[list[float]] | None = None,
                plot_png: str | Path | None = None,
                timeout_s: float = 600.0) -> dict[str, Any]:
    """Temporal statistics over a FILE SERIES (one .msh per step).

    A transient solver writes one mesh per step, which gmsh has no verb
    for: its own time steps live INSIDE a single view.  This treats an
    ordered list of files as the time axis and reduces it two ways:

    * per-tag statistics (``min``, ``max``, ``mean``, ``std``, ``rms``,
      ``ptp``, ``argmin_time``, ``argmax_time``) written as views into
      one output .msh, so "where does the peak occur, and when" is a
      picture rather than a table;
    * per-step global aggregates (min/max/mean/rms over the domain),
      which is the "plot data over time" series -- returned as arrays
      and, with ``plot_png``, drawn.

    The files must share one node/element numbering: a series where the
    mesh changed is not a time series of the same quantity, and mixing
    them silently would average unrelated tags.  That is checked, not
    assumed.

    Args:
        paths: ordered .msh files, one per step.
        view: view name or index (default: the first view of file 0,
            matched by NAME in the others).
        component: component index, or None for the magnitude.
        times: the time value of each file (default: 0, 1, 2, ...).
        stats: which per-tag statistics to write.
        out_file: output .msh (default: ``<first stem>_timestats.msh``).
        points: optional [x, y, z] list probed in EVERY file, giving an
            interpolated history per point (uses the gmsh probe, so it
            is a real field evaluation, not a nearest-node lookup).
        plot_png: draw the aggregate (and point) histories.

    Returns:
        ``{"msh":, "times":, "aggregate": {...}, "point_history": [...],
        "stats_written": [...]}``.
    """
    from .msh_inspect import read_msh_data

    srcs = [Path(p) for p in paths]
    if len(srcs) < 2:
        raise ValueError(
            f"a time series needs at least 2 files, got {len(srcs)}")
    missing = [str(p) for p in srcs if not p.is_file()]
    if missing:
        return {"ok": False, "error": f"file(s) not found: {missing}"}
    bad = [s for s in stats if s not in _TIME_STATS]
    if bad:
        raise ValueError(
            f"unknown stat(s) {bad} (available: {', '.join(_TIME_STATS)})")
    if times is not None and len(times) != len(srcs):
        raise ValueError(
            f"times has {len(times)} entries for {len(srcs)} files")
    t = [float(v) for v in (times if times is not None
                            else range(len(srcs)))]
    if any(not math.isfinite(value) for value in t):
        raise ValueError("times must contain only finite values")

    series: list[dict[int, float]] = []
    view_name = None
    section = None
    components = None
    for k, src in enumerate(srcs):
        data = read_msh_data(src, include_elements=True)
        if not data["views"]:
            return {"ok": False, "error": f"{src.name} carries no view"}
        if k == 0:
            if view is None:
                chosen = data["views"][0]
            elif isinstance(view, bool):
                return {"ok": False,
                        "error": "view must be a name or a non-negative index"}
            elif isinstance(view, int):
                if view < 0 or view >= len(data["views"]):
                    return {"ok": False,
                            "error": f"view index {view} out of range "
                                     f"({len(data['views'])} views)"}
                chosen = data["views"][view]
            else:
                match = [v for v in data["views"] if v["name"] == view]
                if not match:
                    return {"ok": False,
                            "error": f"no view {view!r} in {src.name} "
                                     f"(views: "
                                     f"{[v['name'] for v in data['views']]})"}
                if len(match) > 1:
                    return {"ok": False,
                            "error": f"view {view!r} is ambiguous in "
                                     f"{src.name} ({len(match)} matches)"}
                chosen = match[0]
            view_name = chosen["name"]
            section = chosen["section"]
            components = chosen["components"]
        else:
            match = [v for v in data["views"] if v["name"] == view_name]
            if not match:
                return {"ok": False,
                        "error": f"{src.name} has no view {view_name!r} "
                                 f"(views: "
                                 f"{[v['name'] for v in data['views']]})"}
            if len(match) > 1:
                return {"ok": False,
                        "error": f"view {view_name!r} is ambiguous in "
                                 f"{src.name} ({len(match)} matches)"}
            chosen = match[0]
        if chosen["section"] == "ElementNodeData":
            return {"ok": False,
                    "error": "time_series does not support ElementNodeData; "
                             "convert it to NodeData or ElementData first"}
        if k and chosen["section"] != section:
            return {"ok": False,
                    "error": f"{src.name} stores {view_name!r} as "
                             f"{chosen['section']}, expected {section}"}
        if k and chosen["components"] != components:
            return {"ok": False,
                    "error": f"{src.name} stores {view_name!r} with "
                             f"{chosen['components']} components, expected "
                             f"{components}"}
        vals = _view_values(chosen, component)
        if not vals:
            return {"ok": False,
                    "error": f"view {view_name!r} in {src.name} has no values"}
        if k == 0:
            mesh_nodes = data["nodes"]
            mesh_elements = data["elements"]
            if mesh_nodes:
                mesh_span = max(
                    max(point[axis] for point in mesh_nodes.values())
                    - min(point[axis] for point in mesh_nodes.values())
                    for axis in range(3))
                coordinate_scale = max(
                    abs(value) for point in mesh_nodes.values()
                    for value in point)
                mesh_tolerance = max(
                    mesh_span * 1e-12,
                    math.ulp(max(coordinate_scale, mesh_span, 1e-300)) * 16)
            else:
                mesh_tolerance = 0.0
        elif set(data["nodes"]) != set(mesh_nodes) \
                or data["elements"] != mesh_elements:
            return {"ok": False,
                    "error": f"{src.name} does not share the node/element "
                             f"tag numbering and connectivity of "
                             f"{srcs[0].name} -- a series whose mesh changed "
                             "is not one time series"}
        elif any(abs(coord - mesh_nodes[tag][axis]) > mesh_tolerance
                 for tag, point in data["nodes"].items()
                 for axis, coord in enumerate(point)):
            return {"ok": False,
                    "error": f"{src.name} does not share the node geometry "
                             f"of {srcs[0].name} -- a series whose mesh "
                             "changed is not one time series"}
        if series and set(vals) != set(series[0]):
            return {"ok": False,
                    "error": f"{src.name} does not share the tag numbering "
                             f"of {srcs[0].name} ({len(vals)} vs "
                             f"{len(series[0])} entries) -- a series whose "
                             f"mesh changed is not one time series"}
        series.append(vals)

    tags = sorted(series[0])
    n = len(series)
    per_tag: dict[str, dict[int, float]] = {s: {} for s in stats}
    for tag in tags:
        v = [step[tag] for step in series]
        mean = sum(v) / n
        var = sum((x - mean) ** 2 for x in v) / n
        hi_i = max(range(n), key=lambda i: v[i])
        lo_i = min(range(n), key=lambda i: v[i])
        computed = {"min": v[lo_i], "max": v[hi_i], "mean": mean,
                    "std": math.sqrt(var),
                    "rms": math.sqrt(sum(x * x for x in v) / n),
                    "ptp": v[hi_i] - v[lo_i],
                    "argmax_time": t[hi_i], "argmin_time": t[lo_i]}
        for s in stats:
            per_tag[s][tag] = computed[s]

    aggregate = {"time": t, "min": [], "max": [], "mean": [], "rms": []}
    for step in series:
        v = list(step.values())
        aggregate["min"].append(min(v))
        aggregate["max"].append(max(v))
        aggregate["mean"].append(sum(v) / len(v))
        aggregate["rms"].append(math.sqrt(sum(x * x for x in v) / len(v)))

    out = (Path(out_file) if out_file is not None
           else srcs[0].with_name(srcs[0].stem + "_timestats.msh"))
    out.parent.mkdir(parents=True, exist_ok=True)
    text = srcs[0].read_text(encoding="utf-8", errors="replace")
    starts = [pos for marker in ("$NodeData", "$ElementData",
                                 "$ElementNodeData")
              if (pos := text.find(marker)) >= 0]
    head = text[:min(starts)] if starts else text
    blocks = []
    for s in stats:
        rows = "".join(f"{tag} {per_tag[s][tag]:.9e}\n" for tag in tags)
        blocks.append(f"${section}\n"
                      f'1\n"{view_name}_{s}"\n1\n0\n'
                      f"3\n0\n1\n{len(tags)}\n{rows}$End{section}\n")
    out.write_text(head + "".join(blocks), encoding="utf-8")

    result: dict[str, Any] = {
        "ok": True, "msh": str(out), "view": view_name,
        "section": section, "n_steps": n, "times": t,
        "stats_written": [f"{view_name}_{s}" for s in stats],
        "aggregate": aggregate, "n_tags": len(tags),
        "inputs": [str(p) for p in srcs]}

    if points:
        histories = []
        for p in points:
            histories.append({"point": [float(v) for v in p], "values": []})
        for src in srcs:
            pr = probe_field(src, [list(h["point"]) for h in histories],
                             view=view_name, timeout_s=timeout_s)
            if not pr.get("ok") or not pr.get("views"):
                return {"ok": False,
                        "error": f"probe failed on {src.name}: {pr}"}
            for h, entry in zip(histories, pr["views"][0]["points"]):
                vals = entry.get("values") or []
                if not entry.get("found") or not vals:
                    h["values"].append(None)
                elif len(vals) == 1:
                    h["values"].append(float(vals[0]))
                else:
                    h["values"].append(
                        math.sqrt(sum(float(v) ** 2 for v in vals)))
        result["point_history"] = histories

    if plot_png:
        result["plot"] = _plot_time_series(result, Path(plot_png))
    return result


def _plot_time_series(result: dict[str, Any], png: Path) -> dict[str, Any]:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return {"ok": False, "error": "matplotlib not installed"}
    agg = result["aggregate"]
    hist = result.get("point_history") or []
    fig, axes = plt.subplots(1, 2 if hist else 1,
                             figsize=(9 if hist else 5, 3.4), dpi=150)
    ax = axes[0] if hist else axes
    for key, style in (("max", "-"), ("rms", "--"), ("mean", "-."),
                       ("min", ":")):
        ax.plot(agg["time"], agg[key], style, label=key)
    ax.set_xlabel("time")
    ax.set_ylabel(result["view"])
    ax.legend(fontsize=7)
    ax.grid(alpha=0.3)
    if hist:
        ax2 = axes[1]
        for h in hist:
            lbl = "({:.3g}, {:.3g}, {:.3g})".format(*h["point"])
            xs = [x for x, y in zip(agg["time"], h["values"]) if y is not None]
            ys = [y for y in h["values"] if y is not None]
            ax2.plot(xs, ys, marker="o", ms=2.5, label=lbl)
        ax2.set_xlabel("time")
        ax2.set_ylabel(result["view"])
        ax2.legend(fontsize=7)
        ax2.grid(alpha=0.3)
    png.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(png)
    plt.close(fig)
    return {"ok": True, "png": str(png)}
