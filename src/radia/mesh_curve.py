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

__all__ = ["CurveOrderError", "ensure_curve_order", "mesh_measures", "vol_geometry_kind"]


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


def _read_vol(vol_path, *, with_points=True):
    """Return ``(points, dim, tail_lines)`` of a text Netgen ``.vol``.

    ``points`` is a flat ``array('d')`` of the stored vertex coordinates
    (``dim`` values per point; empty when ``with_points`` is False);
    ``tail_lines`` are the non-empty lines after ``endmesh`` (the geometry
    archive, if any).
    """
    from array import array

    points, dim, tail = array("d"), 0, []
    with open(vol_path, "r", errors="ignore") as handle:
        lines = iter(handle)
        for line in lines:
            s = line.strip()
            if s == "points":
                count = int(next(lines).split()[0])
                for _ in range(count):
                    values = next(lines).split()
                    if with_points:
                        dim = len(values)
                        points.extend(float(v) for v in values)
            elif s == "endmesh":
                tail = [x.strip() for x in lines if x.strip()]
                break
    return points, dim, tail


def vol_geometry_kind(vol_path):
    """Name the geometry archive stored after ``endmesh`` in a Netgen ``.vol``.

    Returns ``"occ"`` for an OCC archive (``TextOutArchive`` header, the
    ``netgen::OCCGeometry`` class tag and a ``CASCADE Topology`` block),
    ``"csg"`` for a CSG description, ``"other"`` for unrecognised trailing
    text and ``None`` when nothing follows ``endmesh``.
    """
    return _geometry_kind(_read_vol(vol_path, with_points=False)[2])


def _geometry_kind(tail):
    if not tail:
        return None
    if (tail[0] == "TextOutArchive" and "netgen::OCCGeometry" in tail[:12]
            and any(x.startswith("CASCADE Topology") for x in tail[:40])):
        return "occ"
    if tail[0].startswith("csgsurfaces"):
        return "csg"
    return "other"


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
        vol_path: the ``.vol`` the mesh was loaded from.  Its embedded OCC
            archive is the only accepted CAD provenance, and the current mesh
            vertices must still equal the points stored in it (a mesh moved or
            rescaled after loading no longer matches its CAD).
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

    stored, dim, tail = _read_vol(vol_path)
    kind = _geometry_kind(tail)
    if kind is None:
        raise CurveOrderError(
            f"{what}: stored geometry order {current} < required {order}, and the "
            f".vol stores no CAD to curve against. Re-export it at order {order} "
            f"(Cubit: export netgen ... order {order}) or regenerate it from the CAD; "
            f"Mesh.Curve would flatten it or use another mesh's geometry.")
    if kind != "occ":
        raise CurveOrderError(
            f"{what}: embedded geometry ({kind}) is not supported for guarded curving "
            f"(OCC archives only). Curve it in the originating mesher.")
    count = len(stored) // dim if dim else 0
    if count != mesh.nv:
        raise CurveOrderError(
            f"{what}: mesh has {mesh.nv} vertices but {vol_path} stores "
            f"{count}; it is not the mesh loaded from that file.")
    scale = max((abs(x) for x in stored), default=0.0)
    worst = 0.0
    for i, v in enumerate(mesh.vertices):
        p = v.point
        base = i * dim
        for k in range(min(dim, len(p))):
            worst = max(worst, abs(p[k] - stored[base + k]))
    if worst > 1e-9 * max(scale, 1e-300):
        raise CurveOrderError(
            f"{what}: mesh vertices differ from the points stored in {vol_path} by up "
            f"to {worst:.3e}; the mesh was moved or rescaled after loading and no "
            f"longer matches its embedded CAD.")
    geometry = mesh.ngmesh.GetGeometry()
    shape = getattr(geometry, "shape", None)
    if shape is None:
        raise CurveOrderError(
            f"{what}: the .vol stores an OCC archive but the loaded mesh exposes "
            f"{type(geometry).__name__}; reload the mesh from {vol_path}.")

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
