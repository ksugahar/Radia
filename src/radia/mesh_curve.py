"""Guarded geometry-order handling for NGSolve meshes loaded from ``.vol``.

``ngsolve.Mesh("x.vol")`` already applies the curved-element data stored in
the file, so the saved geometry order is in effect without calling
``Mesh.Curve``.  Calling ``Curve`` on a loaded mesh rebuilds the curved nodes
from the geometry Netgen associates with the mesh, which is safe only when
that file embeds its own CAD in the same units as the mesh:

* Netgen/OCC ``.vol`` saved with its geometry: the CAD stays attached to that
  mesh even after other files are loaded, and ``Curve(p)`` is correct.
* ``.vol`` meshed in millimetres and rescaled with ``ngmesh.Scale(1e-3)``:
  the embedded CAD stays in millimetres and ``Curve(p >= 2)`` projects the
  nodes onto it (surface area grows by ~1e6..1e8).
* ``.vol`` without CAD (e.g. Cubit ``export netgen``): in a fresh process
  ``Curve(p)`` flattens the stored curving while ``GetCurveOrder()`` reports
  ``p``; after another geometry-bearing mesh was loaded in the same process,
  ``GetGeometry()`` returns that other geometry and ``Curve`` uses it.

Because ``GetGeometry()`` cannot tell an embedded CAD from a borrowed one,
``ensure_curve_order`` takes the provenance from the ``.vol`` file itself.
Only OCC geometry is supported for curving.  The caller owns ``TaskManager``.
See ``docs/ngsolve_integration/curve_order.md``.
"""
from __future__ import annotations

import math

__all__ = ["CurveOrderError", "ensure_curve_order", "mesh_measures", "vol_embeds_geometry"]


class CurveOrderError(RuntimeError):
    """The requested geometry order cannot be provided safely.

    ``mesh_modified`` is True when ``Curve`` already ran and the result was
    rejected; that mesh must be discarded and reloaded.
    """

    def __init__(self, message, *, mesh_modified=False):
        super().__init__(message)
        self.mesh_modified = mesh_modified


def mesh_measures(mesh):
    """Return ``(boundary_measure, domain_measure)`` on the current geometry map.

    For a 3D mesh these are boundary area and volume; for a 2D mesh boundary
    length and area.  Runs in the caller's ``TaskManager`` if any.
    """
    import ngsolve as ng

    one = ng.CoefficientFunction(1.0)
    return float(ng.Integrate(one, mesh, ng.BND)), float(ng.Integrate(one, mesh, ng.VOL))


def vol_embeds_geometry(vol_path):
    """True when a Netgen ``.vol`` file stores a geometry archive after ``endmesh``."""
    seen_end = False
    with open(vol_path, "r", errors="ignore") as handle:
        for line in handle:
            if seen_end:
                if line.strip():
                    return True
            elif line.strip() == "endmesh":
                seen_end = True
    return False


def _vertex_box(mesh):
    lo = [math.inf] * 3
    hi = [-math.inf] * 3
    dim = 0
    for v in mesh.vertices:
        p = v.point
        dim = len(p)
        for k in range(dim):
            lo[k] = min(lo[k], p[k])
            hi[k] = max(hi[k], p[k])
    return lo[:dim], hi[:dim]


def ensure_curve_order(mesh, order, *, vol_path, extent_ratio=2.0, measure_ratio=3.0,
                       what="mesh"):
    """Make the geometry order of a mesh loaded from ``vol_path`` at least ``order``.

    Args:
        mesh: ``ngsolve.Mesh`` loaded from ``vol_path`` (possibly re-wrapped
            with ``ngsolve.Mesh(mesh.ngmesh)``).
        order: required geometry order (>= 1).
        vol_path: the ``.vol`` the mesh was loaded from; its embedded geometry
            archive is the only accepted CAD provenance.
        extent_ratio: allowed factor between the CAD and mesh bounding-box
            extents per axis (unit check; coarse meshes stay well inside 2).
        measure_ratio: allowed factor by which curving may change the boundary
            and domain measures.  Correct curving changes them by
            ``O((h/R)^2)``; unit errors or a wrong CAD change them by orders of
            magnitude.  This is a coarse sanity check, not a proof of validity.
        what: label used in error messages.

    Returns:
        dict with ``action`` (``"kept"``: the stored order suffices and the
        geometry was not re-examined, or ``"curved"``) and ``curve_order``.

    Raises:
        CurveOrderError: no embedded CAD, unsupported geometry type, CAD/mesh
            extent mismatch, or an implausible change after curving
            (``mesh_modified=True``; reload the mesh).
    """
    import ngsolve as ng

    order = int(order)
    if order < 1:
        raise ValueError("order must be >= 1")
    current = int(mesh.GetCurveOrder())
    if current >= order:
        return {"action": "kept", "curve_order": current}

    if not vol_embeds_geometry(vol_path):
        raise CurveOrderError(
            f"{what}: stored geometry order {current} < required {order}, and the "
            f".vol stores no CAD to curve against. Re-export it at order {order} "
            f"(Cubit: export netgen ... order {order}) or regenerate it from the CAD; "
            f"Mesh.Curve would flatten it or use another mesh's geometry.")
    geometry = mesh.ngmesh.GetGeometry()
    shape = getattr(geometry, "shape", None)
    if shape is None:
        raise CurveOrderError(
            f"{what}: embedded geometry type {type(geometry).__name__} is not supported "
            f"for guarded curving (OCC only). Curve it in the originating mesher.")

    lo_c, hi_c = shape.bounding_box
    lo_m, hi_m = _vertex_box(mesh)
    for k in range(len(lo_m)):
        ec, em = float(hi_c[k] - lo_c[k]), hi_m[k] - lo_m[k]
        if em <= 0.0 or ec <= 0.0:
            continue
        shift = abs(0.5 * (hi_c[k] + lo_c[k]) - 0.5 * (hi_m[k] + lo_m[k]))
        if not (1.0 / extent_ratio <= ec / em <= extent_ratio) or shift > 0.1 * max(ec, em):
            raise CurveOrderError(
                f"{what}: embedded CAD box {tuple(lo_c)}..{tuple(hi_c)} does not match "
                f"the mesh box {tuple(lo_m)}..{tuple(hi_m)}. The mesh was probably "
                f"rescaled after meshing (ngmesh.Scale); scale the CAD shape before "
                f"GenerateMesh instead.")

    b0, d0 = mesh_measures(mesh)
    mesh.Curve(order)
    b1, d1 = mesh_measures(mesh)
    rb = b1 / b0 if b0 else math.inf
    rd = d1 / d0 if d0 else math.inf
    ok = all(math.isfinite(r) and 1.0 / measure_ratio <= r <= measure_ratio for r in (rb, rd))
    if not ok:
        raise CurveOrderError(
            f"{what}: Curve({order}) scaled the boundary measure by {rb:.3e} and the "
            f"domain measure by {rd:.3e} (allowed factor {measure_ratio}). The geometry "
            f"map is invalid; reload the mesh and fix the CAD/mesh pairing.",
            mesh_modified=True)
    return {"action": "curved", "curve_order": int(mesh.GetCurveOrder()),
            "from_order": current, "boundary_ratio": rb, "domain_ratio": rd}
