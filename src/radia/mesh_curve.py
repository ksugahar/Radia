"""Guarded geometry-order handling for NGSolve meshes loaded from ``.vol``.

``ngsolve.Mesh("x.vol")`` already applies the curved-element data stored in
the file, so the saved geometry order is in effect without calling
``Mesh.Curve``.  Calling ``Curve`` on a loaded mesh rebuilds the curved nodes
from whatever geometry Netgen currently associates with it, which is only
safe when the file embeds its own CAD in the same units as the mesh:

* Netgen/OCC ``.vol`` saved with its geometry: ``Curve(p)`` is correct.
* ``.vol`` meshed in millimetres and rescaled with ``ngmesh.Scale(1e-3)``:
  the embedded CAD stays in millimetres and ``Curve(p >= 2)`` projects the
  nodes onto it (surface area grows by ~1e6..1e8).
* ``.vol`` without CAD (e.g. Cubit ``export netgen``): in a fresh process
  ``Curve(p)`` flattens the stored curving while ``GetCurveOrder()`` reports
  ``p``; after another geometry-bearing mesh was loaded in the same process,
  Netgen reuses that process-global geometry instead.

``ensure_curve_order`` keeps a sufficient stored order untouched and curves
only after checking that an embedded CAD with matching extent exists, then
verifies boundary area and volume.  See
``docs/ngsolve_integration/curve_order.md``.
"""
from __future__ import annotations

import math

__all__ = ["CurveOrderError", "ensure_curve_order", "mesh_measures"]


class CurveOrderError(RuntimeError):
    """The requested geometry order cannot be provided safely."""


def mesh_measures(mesh):
    """Return ``(boundary_area, volume)`` integrated on the current geometry map."""
    import ngsolve as ng

    one = ng.CoefficientFunction(1.0)
    with ng.TaskManager():
        area = float(ng.Integrate(one, mesh, ng.BND))
        volume = float(ng.Integrate(one, mesh, ng.VOL))
    return area, volume


def _vertex_extent(mesh):
    lo = [math.inf] * 3
    hi = [-math.inf] * 3
    for v in mesh.vertices:
        p = v.point
        for k in range(len(p)):
            lo[k] = min(lo[k], p[k])
            hi[k] = max(hi[k], p[k])
    return [h - l for l, h in zip(lo, hi) if math.isfinite(l)]


def _shape_extent(geometry):
    shape = getattr(geometry, "shape", None)
    if shape is None:
        return None
    lo, hi = shape.bounding_box
    return [float(hi[k] - lo[k]) for k in range(3)]


def ensure_curve_order(mesh, order, *, extent_rtol=0.05, measure_rtol=0.5, what="mesh"):
    """Make the geometry order of ``mesh`` at least ``order`` without silent damage.

    Args:
        mesh: ``ngsolve.Mesh`` (typically loaded from ``.vol``).
        order: required geometry order (>= 1).
        extent_rtol: allowed relative mismatch between the embedded CAD
            bounding box and the mesh vertex bounding box (unit check).
        measure_rtol: allowed relative change of boundary area and volume
            caused by curving. Curving changes them by ``O((h/R)^2)`` (about
            10% for two elements per radius); the wrong geometry changes them
            by orders of magnitude.
        what: label used in error messages.

    Returns:
        dict with ``action`` (``"kept"`` or ``"curved"``), ``curve_order`` and,
        when curved, the relative area and volume changes.

    Raises:
        CurveOrderError: no embedded CAD, CAD/mesh extent mismatch, or an
            implausible change after curving. In the last case the mesh has
            already been modified and must be reloaded.
    """
    import ngsolve as ng

    order = int(order)
    if order < 1:
        raise ValueError("order must be >= 1")
    current = int(mesh.GetCurveOrder())
    if current >= order:
        return {"action": "kept", "curve_order": current}

    shape_extent = _shape_extent(mesh.ngmesh.GetGeometry())
    if shape_extent is None:
        raise CurveOrderError(
            f"{what}: stored geometry order {current} < required {order}, and the "
            f"mesh carries no CAD to curve against. Re-export the .vol at order "
            f"{order} (Cubit: export netgen ... order {order}) or regenerate it from "
            f"the CAD; Mesh.Curve would flatten it or borrow another geometry.")
    mesh_extent = _vertex_extent(mesh)
    for s, m in zip(shape_extent, mesh_extent):
        scale = max(abs(s), abs(m), 1e-300)
        if abs(s - m) > extent_rtol * scale:
            raise CurveOrderError(
                f"{what}: embedded CAD extent {shape_extent} does not match the mesh "
                f"extent {mesh_extent}. The mesh was probably rescaled after meshing "
                f"(ngmesh.Scale) or the geometry belongs to another mesh; scale the "
                f"CAD shape before GenerateMesh instead.")

    area0, vol0 = mesh_measures(mesh)
    with ng.TaskManager():
        mesh.Curve(order)
    area1, vol1 = mesh_measures(mesh)
    d_area = area1 / area0 - 1.0 if area0 else math.inf
    d_vol = vol1 / vol0 - 1.0 if vol0 else 0.0
    if (not math.isfinite(d_area) or not math.isfinite(d_vol)
            or abs(d_area) > measure_rtol or abs(d_vol) > measure_rtol):
        raise CurveOrderError(
            f"{what}: Curve({order}) changed boundary area by {d_area:+.3e} and "
            f"volume by {d_vol:+.3e} (limit {measure_rtol}). The geometry map is "
            f"invalid; reload the mesh and fix the CAD/mesh pairing.")
    return {"action": "curved", "curve_order": int(mesh.GetCurveOrder()),
            "from_order": current, "area_rel_change": d_area, "volume_rel_change": d_vol}
