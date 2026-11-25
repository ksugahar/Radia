"""Equal-flux field-line seeds and field-line tracing on a Gmsh vector view.

Thin adapter: B is probed from the MSH view, and seeds come from the
validated radia.flux_seeding two-stage inverse. Optionally the seeds feed
the existing 3D streamline tracer. Seed placement, its equal-flux
quadrature check and trace termination are reported separately.

trace_field_lines is the field-line entry point: equal-flux seeding on a
patch is the default choice, so line density follows flux; uniform seeds
along a segment stay available only on explicit request and are labelled
as not flux-weighted.

MSH carries no length unit. length_unit states the unit of the file's
coordinates ("m" or "mm"); fluxes are always reported in Wb.
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
                           "out_file", "timeout_s"})
LENGTH_UNITS = {"m": 1.0, "mm": 1e-3}
SEED_PATCH_KEYS = frozenset({"origin", "e_xi", "e_eta", "n_u", "n_v", "table_nodes",
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


def equal_flux_seeds_on_plane(path: str | Path, origin: list[float], e_xi: list[float],
                              e_eta: list[float], *, n_u: int, n_v: int,
                              table_nodes: int = 33, view: str | int | None = None,
                              time_step: int = 0, quadrature_order: int = 6,
                              trace: bool = False, trace_options: dict | None = None,
                              length_unit: str = "m",
                              timeout_s: float = 300.0) -> dict[str, Any]:
    """Seed the patch origin + xi e_xi + eta e_eta, xi, eta in [0, 1], with
    n_u * n_v equal-flux seeds of the probed B view (B in T, fluxes in Wb).

    origin, e_xi and e_eta are in the file's coordinates, whose unit is
    length_unit ("m" or "mm"). Each seed reports "point" in file
    coordinates and "point_m" in metres; tracing uses file coordinates.
    view None selects the first view (index 0) for both probing and tracing.
    trace_options may hold only TRACE_OPTIONS keys; others are rejected
    before any probing.
    """
    err = _check_input(path)
    if err:
        return err
    options = dict(trace_options or {})
    bad = sorted(set(options) - TRACE_OPTIONS)
    if bad:
        return {"ok": False, "error": f"unsupported trace_options {bad}; "
                                      f"allowed: {sorted(TRACE_OPTIONS)}"}
    try:
        time_step = int(time_step)
        if time_step < 0:
            raise ValueError
    except (TypeError, ValueError, OverflowError):
        return {"ok": False, "error": "time_step must be a non-negative integer"}
    if length_unit not in LENGTH_UNITS:
        return {"ok": False, "error": f"length_unit must be one of {sorted(LENGTH_UNITS)}; "
                                      "MSH carries no unit, so state the file's coordinate unit"}
    scale = LENGTH_UNITS[length_unit]
    view = 0 if view is None else view
    try:
        from radia.flux_seeding import cell_flux_quadrature, equal_flux_seeds, planar_patch
    except ImportError as exc:
        return {"ok": False, "error": f"radia.flux_seeding is unavailable in radia "
                                      f"{_radia_version()}: {exc}; update radia"}
    path = Path(path)
    try:
        if int(table_nodes) < 2:
            raise ValueError("table_nodes must be at least 2")
        patch = planar_patch(*(scale * np.asarray(v, dtype=float) for v in (origin, e_xi, e_eta)))
        g = patch.flux_density(_probe_b(path, view, time_step, timeout_s, scale))
        nodes = np.linspace(0.0, 1.0, int(table_nodes))
        seeds = equal_flux_seeds(g, xi_nodes=nodes, eta_nodes=nodes, n_u=n_u, n_v=n_v,
                                 chart=patch.chart)
        check = None
        if int(quadrature_order) > 0:
            q = cell_flux_quadrature(g, seeds, xi_nodes=nodes, eta_nodes=nodes,
                                     order=int(quadrature_order))
            check = {k: q[k] for k in ("independent_total_flux_Wb",
                                       "max_relative_cell_deviation", "coverage", "order")}
    except (ValueError, RuntimeError) as exc:
        return {"ok": False, "error": str(exc)}

    result: dict[str, Any] = {
        "ok": True,
        "seeds": [{"id": int(i), "point": (p / scale).tolist(), "point_m": p.tolist(),
                   "chart": c.tolist(), "flux": f.tolist()}
                  for i, p, c, f in zip(seeds.ids, seeds.points_m, seeds.chart_coords,
                                        seeds.flux_coords)],
        "n_u": seeds.n_u, "n_v": seeds.n_v,
        "orientation": seeds.orientation,
        "total_flux_Wb": seeds.total_flux_Wb,
        "target_cell_flux_Wb": seeds.target_cell_flux_Wb,
        "placement_diagnostics": seeds.diagnostics,
        "equal_flux_quadrature": check,
        "provenance": {**seeds.provenance, "field_samples": "gmsh view probe",
                       "msh_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                       "view": view, "time_step": time_step,
                       "table_nodes": int(table_nodes), "length_unit": length_unit,
                       "patch": {"origin": list(map(float, origin)),
                                 "e_xi": list(map(float, e_xi)),
                                 "e_eta": list(map(float, e_eta))}},
        "limitations": "planar patch only; field samples are the view's interpolation; "
                       "no field-line closure or CAD splitting is implied",
    }
    if trace:
        traced = streamlines(path, seeds=(seeds.points_m / scale).tolist(), view=view,
                             time_step=time_step, **options)
        result["trace"] = traced
    return result


def trace_field_lines(path: str | Path, *, seed_patch: dict | None = None,
                      seed_line: dict | None = None, view: str | int | None = None,
                      time_step: int = 0, length_unit: str = "m",
                      **trace_options) -> dict[str, Any]:
    """Trace field lines of a vector view; equal-flux seeding is the default choice.

    seed_patch = {origin, e_xi, e_eta, n_u, n_v[, table_nodes, quadrature_order]}
    places n_u * n_v seeds that each carry the same flux through the patch
    (radia.flux_seeding), so line density follows the field. seed_line =
    {start, end, n} places n uniform seeds on a segment; it is not
    flux-weighted and is reported as such. Exactly one must be given.
    trace_options are the tracer keywords in TRACE_OPTIONS. The result is
    the tracer's result plus "seeding", which records the method and, for
    equal flux, total and per-line flux in Wb and the placement provenance.
    """
    if (seed_patch is None) == (seed_line is None):
        return {"ok": False, "error": "give exactly one of seed_patch (equal-flux seeds on a "
                                      "planar patch; the default choice) or seed_line (uniform "
                                      "seeds on a segment, not flux-weighted)"}
    bad = sorted(set(trace_options) - TRACE_OPTIONS)
    if bad:
        return {"ok": False, "error": f"unsupported tracer options {bad}; "
                                      f"allowed: {sorted(TRACE_OPTIONS)}"}
    if seed_patch is not None:
        patch = dict(seed_patch)
        missing = sorted({"origin", "e_xi", "e_eta", "n_u", "n_v"} - set(patch))
        extra = sorted(set(patch) - SEED_PATCH_KEYS)
        if missing or extra:
            return {"ok": False, "error": f"seed_patch needs {sorted(SEED_PATCH_KEYS)}; "
                                          f"missing {missing}, unknown {extra}"}
        placed = equal_flux_seeds_on_plane(
            path, patch["origin"], patch["e_xi"], patch["e_eta"], n_u=patch["n_u"],
            n_v=patch["n_v"], table_nodes=patch.get("table_nodes", 33),
            quadrature_order=patch.get("quadrature_order", 6), view=view,
            time_step=time_step, trace=True, trace_options=trace_options,
            length_unit=length_unit)
        if not placed.get("ok"):
            return placed
        result = dict(placed["trace"])
        result["seeding"] = {
            "method": "equal_flux", "flux_weighted": True,
            "n_u": placed["n_u"], "n_v": placed["n_v"], "orientation": placed["orientation"],
            "total_flux_Wb": placed["total_flux_Wb"],
            "flux_per_line_Wb": placed["target_cell_flux_Wb"],
            "equal_flux_quadrature": placed["equal_flux_quadrature"],
            "placement_diagnostics": placed["placement_diagnostics"],
            "seeds": placed["seeds"], "provenance": placed["provenance"],
            "limitations": placed["limitations"]}
        return result
    line = dict(seed_line)
    if set(line) != SEED_LINE_KEYS:
        return {"ok": False, "error": f"seed_line needs exactly {sorted(SEED_LINE_KEYS)}"}
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
