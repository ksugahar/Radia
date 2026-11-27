"""Equal-flux field-line seeds and field-line tracing on a Gmsh vector view.

Thin adapter: B is probed from the MSH view, and seeds come from the
validated radia.flux_seeding layouts on a seed surface (a parallelogram,
a planar annulus or a cylinder sector):

- n_lines: recursive equal-flux bisection (any count, near-square cells,
  one sign of a sign-changing surface);
- n_u, n_v: the two-stage inverse grid (one-signed surfaces).

trace_field_lines is the field-line entry point: equal-flux seeding is the
default choice, so line density follows flux; uniform seeds along a segment
stay available only on explicit request and are labelled as not
flux-weighted. After tracing, a flux account reports how each line ended
(closed, open at the data or extent boundary, truncated by max_steps,
stagnated), how often lines cross the seed surface again and how much of
the opposite-sign flux the returning lines draw. Line density stands for
flux only along lines that are neither truncated nor drawn twice.

MSH carries no length unit. The picture and the seed positions do not
depend on it; length_unit ("m" or "mm") only scales the reported fluxes
to Wb.
"""

from __future__ import annotations

import hashlib
from pathlib import Path
from typing import Any

import numpy as np

from .post_process import _check_input, probe_field, streamlines

# Tracer keywords a caller may set. The adapter owns path, seeds, view and
# time_step so that placement and tracing read the same view and step.
TRACE_OPTIONS = frozenset({"step_size", "max_steps", "both_directions", "adaptive",
                           "closure", "max_turn_deg", "arrows_every", "return_points",
                           "extent", "out_file", "timeout_s"})
LENGTH_UNITS = {"m": 1.0, "mm": 1e-3}
PATCH_SHAPES = {
    "parallelogram": (("origin", "e_xi", "e_eta"), ()),
    "annulus": (("center", "normal", "r_outer"), ("r_inner", "theta_range", "reference")),
    "cylinder": (("axis_origin", "axis", "radius"), ("theta_range", "reference")),
}
LENGTH_KEYS = frozenset({"origin", "e_xi", "e_eta", "center", "r_outer", "r_inner",
                         "axis_origin", "axis", "radius"})
PATCH_COMMON = frozenset({"shape", "n_lines", "n_u", "n_v", "sign", "table_nodes",
                          "quadrature_order"})
SEED_LINE_KEYS = frozenset({"start", "end", "n"})


def _radia_version() -> str:
    try:
        import importlib.metadata
        return importlib.metadata.version("radia")
    except Exception:  # noqa: BLE001 - diagnostic only
        return "unknown"


def _probe_b(path: Path, view, time_step: int, timeout_s: float, scale: float = 1.0):
    """B at points given in metres; the file's coordinates are metres / scale."""
    def field(points):
        pts = np.asarray(points, dtype=float)
        flat = pts.reshape(-1, 3) / scale
        result = probe_field(path, flat.tolist(), view=view, step=time_step,
                             timeout_s=timeout_s)
        if not result.get("ok"):
            raise RuntimeError(f"probe failed: {result.get('error')}")
        rows = result["views"][0]["points"]
        missing = [r["point"] for r in rows if not r.get("found")]
        if missing:
            raise ValueError(f"{len(missing)} patch samples lie outside the field data, "
                             f"first {missing[0]}")
        values = np.asarray([r["values"] for r in rows], dtype=float)
        if values.shape != (len(flat), 3):
            raise ValueError("the selected view is not a 3-component vector field")
        return values.reshape(pts.shape)

    return field


def _check_patch(spec: Any) -> str | None:
    """Return an error for a malformed seed_patch, before any probing."""
    if not isinstance(spec, dict):
        return "seed_patch must be a dict"
    shape = spec.get("shape", "parallelogram")
    if not isinstance(shape, str) or shape not in PATCH_SHAPES:
        return f"seed_patch shape must be one of {sorted(PATCH_SHAPES)}, got {shape!r}"
    required, optional = PATCH_SHAPES[shape]
    missing = sorted(set(required) - set(spec))
    extra = sorted(set(spec) - set(required) - set(optional) - PATCH_COMMON)
    if missing or extra:
        return (f"seed_patch ({shape}) needs {sorted(required)}, optional "
                f"{sorted(set(optional) | PATCH_COMMON)}; missing {missing}, unknown {extra}")
    has_lines, has_grid = "n_lines" in spec, ("n_u" in spec or "n_v" in spec)
    if has_lines == has_grid or (has_grid and not {"n_u", "n_v"} <= set(spec)):
        return ("seed_patch needs either n_lines (balanced bisection) or both n_u and n_v "
                "(two-stage grid)")
    if has_grid and "sign" in spec:
        return ("sign applies to n_lines only; the n_u x n_v grid needs a one-signed surface")
    for key, value in spec.items():
        problem = _check_value(key, value)
        if problem:
            return f"seed_patch {key} {problem}"
    if has_grid and spec["n_u"] * spec["n_v"] > MAX_LINES:
        return f"seed_patch n_u x n_v = {spec['n_u'] * spec['n_v']} exceeds {MAX_LINES} lines"
    # The table and the independent total quadrature are known before anything is probed;
    # the rest is counted as it is probed (radia.flux_seeding.SampleBudget).
    nodes, order = spec.get("table_nodes", 33), spec.get("quadrature_order", 3)
    need = nodes ** 2 + ((nodes - 1) * order) ** 2
    if need > FIELD_SAMPLE_LIMIT:
        return (f"seed_patch table_nodes={nodes} with quadrature_order={order} needs about "
                f"{need:,} field samples, above the {FIELD_SAMPLE_LIMIT:,} limit; reduce either")
    return None


VECTOR_KEYS = frozenset({"origin", "e_xi", "e_eta", "center", "normal", "axis_origin", "axis",
                         "reference", "start", "end"})
SCALAR_KEYS = frozenset({"r_outer", "r_inner", "radius"})
COUNT_MINIMUM = {"n_lines": 1, "n_u": 1, "n_v": 1, "n": 1, "table_nodes": 2,
                 "quadrature_order": 0}
# Resource guards against runaway allocation, not accuracy limits; radia.flux_seeding
# also bounds the field samples of every grid it evaluates (MAX_FIELD_SAMPLES).
MAX_LINES = 1_000  # user-set 2026-11-26: up to about 1000 field lines per seeding
FIELD_SAMPLE_LIMIT = 100_000  # user-set 2026-11-26: field samples per seeding, all probes
COUNT_MAXIMUM = {"n_lines": MAX_LINES, "n_u": MAX_LINES, "n_v": MAX_LINES, "n": MAX_LINES,
                 "table_nodes": 316, "quadrature_order": 64}


def _is_number(value: Any) -> bool:
    """A finite int or float; an int too large for a float is not one."""
    if not isinstance(value, (int, float)) or isinstance(value, bool):
        return False
    try:
        return bool(np.isfinite(float(value)))
    except OverflowError:
        return False


def _check_value(key: str, value: Any) -> str | None:
    """Return why one seed_patch or seed_line value is malformed, or None.

    JSON null, booleans, strings and non-finite numbers are rejected here so
    that a schema-accepted payload never reaches int() or float() unchecked.
    """
    if key in VECTOR_KEYS:
        if (not isinstance(value, (list, tuple)) or len(value) != 3
                or not all(_is_number(v) for v in value)):
            return f"must be a list of 3 finite numbers, got {value!r}"
    elif key in SCALAR_KEYS:
        if not _is_number(value):
            return f"must be a finite number, got {value!r}"
    elif key == "theta_range":
        if (not isinstance(value, (list, tuple)) or len(value) != 2
                or not all(_is_number(v) for v in value)):
            return f"must be [theta0, theta1] in radians, got {value!r}"
    elif key in COUNT_MINIMUM:
        if (not isinstance(value, int) or isinstance(value, bool)
                or not COUNT_MINIMUM[key] <= value <= COUNT_MAXIMUM[key]):
            shown = value if not isinstance(value, int) or abs(value) < 10**30 else "a huge integer"
            return (f"must be an integer in [{COUNT_MINIMUM[key]}, {COUNT_MAXIMUM[key]}], "
                    f"got {shown!r}")
    elif key == "sign":
        if value not in ("auto", "positive", "negative"):
            return f"must be 'auto', 'positive' or 'negative', got {value!r}"
    return None


def _chart(spec: dict, scale: float):
    from radia.flux_seeding import annular_patch, cylindrical_patch, planar_patch

    v = {k: (scale * np.asarray(val, dtype=float) if k in LENGTH_KEYS else val)
         for k, val in spec.items()}
    shape = spec.get("shape", "parallelogram")
    if shape == "parallelogram":
        return planar_patch(v["origin"], v["e_xi"], v["e_eta"])
    extra = {k: v[k] for k in ("theta_range", "reference") if k in v}
    if shape == "annulus":
        return annular_patch(v["center"], v["normal"], float(v["r_outer"]),
                             r_inner_m=float(v.get("r_inner", 0.0)), **extra)
    return cylindrical_patch(v["axis_origin"], v["axis"], float(v["radius"]), **extra)


def _place(path, spec, *, view, time_step, length_unit, timeout_s):
    """Place seeds; returns (result, chart, scale). chart is None on failure."""
    err = _check_input(path)
    if err:
        return err, None, None
    problem = _check_patch(spec)
    if problem:
        return {"ok": False, "error": problem}, None, None
    try:
        time_step = int(time_step)
        if time_step < 0:
            raise ValueError
    except (TypeError, ValueError, OverflowError):
        return {"ok": False, "error": "time_step must be a non-negative integer"}, None, None
    if length_unit not in LENGTH_UNITS:
        return {"ok": False, "error": f"length_unit must be one of {sorted(LENGTH_UNITS)}; "
                                      "MSH carries no unit, so state the file's coordinate "
                                      "unit"}, None, None
    scale = LENGTH_UNITS[length_unit]
    view = 0 if view is None else view
    try:
        import radia.field_lines  # noqa: F401  (the flux account needs it)
        from radia.flux_seeding import (SampleBudget, balanced_cell_flux_quadrature,
                                        balanced_equal_flux_seeds, cell_flux_quadrature,
                                        equal_flux_seeds)
    except ImportError as exc:
        return {"ok": False, "error": "radia.flux_seeding or radia.field_lines is unavailable "
                                      f"or too old in radia {_radia_version()}: {exc}; "
                                      "update radia"}, None, None
    path = Path(path)
    table_nodes = spec.get("table_nodes", 33)
    order = int(spec.get("quadrature_order", 3))
    balanced = "n_lines" in spec
    try:
        if int(table_nodes) < 2:
            raise ValueError("table_nodes must be at least 2")
        chart = _chart(spec, scale)
        # One budget for the whole seeding: placement and cell check share it.
        g = SampleBudget(chart.flux_density(_probe_b(path, view, time_step, timeout_s, scale)),
                         FIELD_SAMPLE_LIMIT)
        nodes = np.linspace(0.0, 1.0, int(table_nodes))
        if balanced:
            seeds = balanced_equal_flux_seeds(g, xi_nodes=nodes, eta_nodes=nodes,
                                              n_lines=spec["n_lines"], chart=chart.chart,
                                              metric=chart.metric,
                                              sign=spec.get("sign", "auto"))
            check_fn = balanced_cell_flux_quadrature
        else:
            seeds = equal_flux_seeds(g, xi_nodes=nodes, eta_nodes=nodes, n_u=spec["n_u"],
                                     n_v=spec["n_v"], chart=chart.chart, metric=chart.metric)
            check_fn = cell_flux_quadrature
        check = None
        if order > 0:
            q = check_fn(g, seeds, xi_nodes=nodes, eta_nodes=nodes, order=order)
            check = {k: q[k] for k in ("independent_total_flux_Wb",
                                       "max_relative_cell_deviation", "coverage", "order")}
    except (ValueError, RuntimeError) as exc:
        return {"ok": False, "error": str(exc)}, None, None

    entries = []
    for i, p, c in zip(seeds.ids, seeds.points_m, seeds.chart_coords):
        entries.append({"id": int(i), "point": (p / scale).tolist(), "point_m": p.tolist(),
                        "chart": c.tolist()})
    if balanced:
        counts = {"layout": "balanced", "n_lines": seeds.n_lines}
        parts = {"positive_flux_Wb": seeds.positive_flux_Wb,
                 "negative_flux_Wb": seeds.negative_flux_Wb}
    else:
        for entry, f in zip(entries, seeds.flux_coords):
            entry["flux"] = f.tolist()
        counts = {"layout": "grid", "n_lines": seeds.n_u * seeds.n_v,
                  "n_u": seeds.n_u, "n_v": seeds.n_v}
        one = seeds.total_flux_Wb
        parts = {"positive_flux_Wb": max(one, 0.0), "negative_flux_Wb": min(one, 0.0)}
    result: dict[str, Any] = {
        "ok": True,
        "seeds": entries,
        **counts,
        "orientation": seeds.orientation,
        "total_flux_Wb": seeds.total_flux_Wb,
        "flux_per_line_Wb": seeds.target_cell_flux_Wb,
        **parts,
        "placement_diagnostics": seeds.diagnostics,
        "equal_flux_quadrature": check,
        "provenance": {**seeds.provenance, "field_samples": "gmsh view probe",
                       "msh_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                       "view": view, "time_step": time_step,
                       "table_nodes": int(table_nodes), "length_unit": length_unit,
                       "seed_patch": {k: (np.asarray(v, float).tolist()
                                          if k in LENGTH_KEYS or k in ("normal", "reference",
                                                                       "theta_range")
                                          else v) for k, v in spec.items()}},
        "limitations": "field samples are the view's interpolation; seeds cover one flux "
                       "sign; no field-line closure or CAD splitting is implied",
    }
    return result, chart, scale


def place_equal_flux_seeds(path: str | Path, seed_patch: dict, *,
                           view: str | int | None = None, time_step: int = 0,
                           length_unit: str = "m", timeout_s: float = 300.0) -> dict[str, Any]:
    """Place equal-flux seeds on seed_patch (no tracing); see trace_field_lines."""
    return _place(path, seed_patch, view=view, time_step=time_step,
                  length_unit=length_unit, timeout_s=timeout_s)[0]


def _flux_account(traced: dict, placed: dict, chart, scale: float) -> tuple[dict, list[str]]:
    """How the traced lines represent the seeded and the opposite flux (radia.field_lines)."""
    from radia.field_lines import FieldLine, flux_account

    lines = [FieldLine(line["seed_index"], np.asarray(poly["points"], dtype=float) * scale, None,
                       poly["seed_vertex"], poly["closed"], dict(line["reasons"]))
             for line, poly in zip(traced["lines"], traced["polylines"])]
    sign = placed["orientation"]
    opposite = placed["negative_flux_Wb"] if sign > 0 else placed["positive_flux_Wb"]
    return flux_account(lines, chart, n_seeds=len(placed["seeds"]),
                        flux_per_line_Wb=placed["flux_per_line_Wb"], orientation=sign,
                        opposite_flux_Wb=opposite)


def trace_field_lines(path: str | Path, *, seed_patch: dict | None = None,
                      seed_line: dict | None = None, view: str | int | None = None,
                      time_step: int = 0, length_unit: str = "m",
                      **trace_options) -> dict[str, Any]:
    """Trace field lines of a vector view; equal-flux seeding is the default choice.

    seed_patch places seeds that each carry the same flux through a seed
    surface (radia.flux_seeding), so line density follows the field:
    shape "parallelogram" {origin, e_xi, e_eta} (default), "annulus"
    {center, normal, r_outer[, r_inner, theta_range, reference]} or
    "cylinder" {axis_origin, axis, radius[, theta_range, reference]}, with
    either n_lines (balanced bisection; sign "auto" | "positive" |
    "negative") or n_u and n_v (two-stage grid), and optional table_nodes
    and quadrature_order. Lengths are in the file's coordinates; angles in
    radians. seed_line = {start, end, n} places n uniform seeds on a
    segment; it is not flux-weighted and is reported as such. Exactly one
    must be given. trace_options are the tracer keywords in TRACE_OPTIONS
    (extent = {"min", "max"} stops lines at a box). The result is the
    tracer's result plus "seeding" and, for equal flux, "flux_account" and
    "warnings". flux_per_line_Wb is the signed nominal flux of each seed's
    cell; tracer lengths and returned points are in the file's coordinates,
    point_m in metres, fluxes in Wb. Malformed values (JSON null, booleans,
    non-finite numbers, wrong lengths) are rejected before probing.
    """
    if (seed_patch is None) == (seed_line is None):
        return {"ok": False, "error": "give exactly one of seed_patch (equal-flux seeds on a "
                                      "seed surface; the default choice) or seed_line (uniform "
                                      "seeds on a segment, not flux-weighted)"}
    bad = sorted(set(trace_options) - TRACE_OPTIONS)
    if bad:
        return {"ok": False, "error": f"unsupported tracer options {bad}; "
                                      f"allowed: {sorted(TRACE_OPTIONS)}"}
    if seed_patch is not None:
        options = dict(trace_options)
        timeout_s = options.get("timeout_s", 300.0)
        placed, chart, scale = _place(path, seed_patch, view=view, time_step=time_step,
                                      length_unit=length_unit, timeout_s=timeout_s)
        if chart is None:
            return placed
        keep_points = bool(options.pop("return_points", False))
        traced = streamlines(path, seeds=[s["point"] for s in placed["seeds"]],
                             view=0 if view is None else view, time_step=time_step,
                             return_points=True, **options)
        if not traced.get("ok"):
            return traced
        account, warnings = _flux_account(traced, placed, chart, scale)
        result = dict(traced)
        if not keep_points:
            result["polylines"] = None
        result["seeding"] = {
            "method": "equal_flux", "flux_weighted": True,
            **{k: placed[k] for k in ("layout", "n_lines", "orientation", "total_flux_Wb",
                                      "flux_per_line_Wb", "positive_flux_Wb",
                                      "negative_flux_Wb", "equal_flux_quadrature",
                                      "placement_diagnostics", "seeds", "provenance",
                                      "limitations")},
            **{k: placed[k] for k in ("n_u", "n_v") if k in placed}}
        result["flux_account"] = account
        result["warnings"] = warnings
        return result
    if not isinstance(seed_line, dict) or set(seed_line) != SEED_LINE_KEYS:
        return {"ok": False, "error": f"seed_line needs exactly {sorted(SEED_LINE_KEYS)}"}
    line = dict(seed_line)
    for key, value in line.items():
        problem = _check_value(key, value)
        if problem:
            return {"ok": False, "error": f"seed_line {key} {problem}"}
    if length_unit not in LENGTH_UNITS:
        return {"ok": False, "error": f"length_unit must be one of {sorted(LENGTH_UNITS)}"}
    result = streamlines(path, line["start"], line["end"], n_seeds=line["n"], view=view,
                         time_step=time_step, **trace_options)
    if result.get("ok"):
        result["seeding"] = {"method": "uniform_line", "flux_weighted": False,
                             "length_unit": length_unit,
                             "note": "line density does not encode flux; use seed_patch "
                                     "for equal-flux seeds"}
    return result
