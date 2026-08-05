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
                    entry = {"point": pt, "found": bool(values),
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
                values, _dist = gmsh.view.probe(tag, q[0], q[1], q[2],
                                                step=step)
                if len(values) < 3:
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
