"""Equal-flux field-line seeds on a planar patch of a Gmsh vector view.

Thin adapter: B is probed from the MSH view, and seeds come from the
validated radia.flux_seeding two-stage inverse. Optionally the seeds feed
the existing 3D streamline tracer. Seed placement, its equal-flux
quadrature check and trace termination are reported separately.
"""

from __future__ import annotations

import hashlib
from pathlib import Path
from typing import Any

import numpy as np

from .post_process import _check_input, probe_field, streamlines


def _probe_b(path: Path, view, time_step: int, timeout_s: float):
    def field(points):
        pts = np.asarray(points, dtype=float)
        flat = pts.reshape(-1, 3)
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
                              timeout_s: float = 300.0) -> dict[str, Any]:
    """Seed the patch origin + xi e_xi + eta e_eta, xi, eta in [0, 1], with
    n_u * n_v equal-flux seeds of the probed B view (SI: m, T, Wb)."""
    err = _check_input(path)
    if err:
        return err
    try:
        from radia.flux_seeding import cell_flux_quadrature, equal_flux_seeds, planar_patch
    except ImportError as exc:
        return {"ok": False, "error": f"radia.flux_seeding is unavailable: {exc}"}
    path = Path(path)
    try:
        if int(table_nodes) < 2:
            raise ValueError("table_nodes must be at least 2")
        patch = planar_patch(origin, e_xi, e_eta)
        g = patch.flux_density(_probe_b(path, view, int(time_step), timeout_s))
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
        "seeds": [{"id": int(i), "point_m": p.tolist(), "chart": c.tolist(),
                   "flux": f.tolist()}
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
                       "view": view, "time_step": int(time_step),
                       "table_nodes": int(table_nodes),
                       "patch": {"origin_m": list(map(float, origin)),
                                 "e_xi_m": list(map(float, e_xi)),
                                 "e_eta_m": list(map(float, e_eta))}},
        "limitations": "planar patch only; field samples are the view's interpolation; "
                       "no field-line closure or CAD splitting is implied",
    }
    if trace:
        traced = streamlines(path, seeds=seeds.points_m.tolist(), view=view,
                             time_step=int(time_step), **(trace_options or {}))
        result["trace"] = traced
    return result
