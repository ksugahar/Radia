"""Density-design -> CAD bridge: nodal level sets, Exodus hand-off, iso STL.

This module is the DESIGN-side half of the topology-optimization shape
regeneration loop (design record:
``docs/hdiv_vim/TOPOLOGY_OPTIMIZATION_ISOCHRONOUS.md``; the MESH-side half
lives on the radia-mcp cubit server, which turns the artifacts written here
into all-hex / tet solver meshes with closure and quality gates):

    per-element density  --(``nodal_from_element_density``)-->  P1 level set
        --> ``write_levelset_exodus``  (Cubit ``create tri iso`` route,
                                        mesh-native, no resampling)
        --> ``iso_stl_from_grid``      (marching-cubes route, pure Python)

Verified facts this module is built on (2026-08-08, LAB, Cubit 2025.12):

* ``import mesh ... nodal_var`` + ``create tri iso tet all`` work headlessly
  on Cubit 2025.12 (beta, needs ``set dev on``); on a radius-0.6 sphere level
  set the extracted surface area was -4.6 % of exact at h=0.25 (linear cut).
* ``import stl`` + ``sculpt volume all`` gives all-hex with 0.28 % volume
  closure; ``scheme tetmesh`` on the same STL gives tets at ~1 % closure.
* netgen's ``STLGeometry`` REJECTS marching-cubes STL ("STL Triangle
  degenerated", 0 tets) even after vertex weld + Taubin smoothing -- recorded
  negative; do not re-route the tet leg through netgen.
* ``trimesh`` Taubin smoothing at 10 iterations shrank a marching-cubes body
  by 12.8 % -- smoothing iterations are capped and the volume drift is
  reported instead of silently accepted.

Optional dependencies (each raises with install guidance when missing, per
the fail-fast policy): ``netCDF4`` (Exodus), ``scikit-image`` (marching
cubes), ``trimesh`` (STL hygiene), ``fast-simplification`` (quadric
decimation via ``target_faces``).
"""
from __future__ import annotations

import numpy as np

__all__ = [
    "nodal_from_element_density",
    "write_levelset_exodus",
    "iso_stl_from_grid",
]


def _require(module_name: str, pip_name: str):
    try:
        return __import__(module_name)
    except ImportError as exc:
        raise ImportError(
            f"radia.topopt_cad requires the optional dependency "
            f"'{module_name}' for this function; install it with "
            f"`pip install {pip_name}`") from exc


# ----------------------------------------------------------------------
# P0 -> P1: volume-weighted nodal average of the element density
# ----------------------------------------------------------------------

def nodal_from_element_density(mesh, density):
    """Volume-weighted vertex average of a per-element density.

    Parameters
    ----------
    mesh : ngsolve.Mesh
        Straight (non-curved) 3D mesh the design lives on.
    density : array-like, shape (mesh.ne,)
        Per-element design density in NGSolve VOL element order.

    Returns
    -------
    numpy.ndarray, shape (mesh.nv,)
        Nodal density: ``sum_e rho_e V_e / sum_e V_e`` over the elements
        incident to each vertex.  This is the P1 field the iso-surface
        routes cut at the threshold level.
    """
    import ngsolve as ng

    rho = np.asarray(density, dtype=float).ravel()
    if mesh.dim != 3:
        raise ValueError("nodal_from_element_density: 3D meshes only")
    if mesh.GetCurveOrder() >= 2:
        raise ValueError(
            "nodal_from_element_density: curved meshes are not supported; "
            "pass the straight design mesh")
    if rho.size != mesh.ne:
        raise ValueError(
            "nodal_from_element_density: density has %d entries, mesh has "
            "%d volume elements" % (rho.size, mesh.ne))
    if not np.all(np.isfinite(rho)):
        raise ValueError("nodal_from_element_density: density has non-finite "
                         "entries")

    volumes = np.asarray(
        ng.Integrate(ng.CoefficientFunction(1.0), mesh, ng.VOL,
                     element_wise=True), dtype=float)
    if volumes.size != mesh.ne or not np.all(volumes > 0.0):
        raise RuntimeError(
            "nodal_from_element_density: element volume integration "
            "returned %d entries / min %g -- mesh is inconsistent"
            % (volumes.size, volumes.min() if volumes.size else np.nan))

    num = np.zeros(mesh.nv)
    den = np.zeros(mesh.nv)
    for el, rho_e, vol_e in zip(mesh.Elements(ng.VOL), rho, volumes):
        for v in el.vertices:
            num[v.nr] += rho_e * vol_e
            den[v.nr] += vol_e
    if not np.all(den > 0.0):
        raise RuntimeError(
            "nodal_from_element_density: %d vertices have no incident "
            "volume element" % int(np.sum(den <= 0.0)))
    return num / den


# ----------------------------------------------------------------------
# Exodus II hand-off for the Cubit `create tri iso` (ATO) route
# ----------------------------------------------------------------------

def write_levelset_exodus(mesh, nodal, path, *, level=0.5, varname="LSD"):
    """Write a minimal Exodus II file: the TET mesh + one nodal variable.

    The variable holds ``nodal - level`` so the Cubit command family

        set dev on
        import mesh "<path>" nodal_var "LSD" no_geom
        create tri iso tet all nodal_var "LSD"

    cuts the iso-surface at its native level 0 (the official ATO recipe,
    ``docs/help/appendix/ato_to_mesh.htm`` of Cubit 2025.8+).  Two triangle
    blocks appear after the command: per the official documentation the
    LARGER new block id carries the optimized free surface, the smaller one
    the triangles on fixed (design-boundary) portions.

    Returns a dict with the written counts.
    """
    _require("netCDF4", "netCDF4")
    from netCDF4 import Dataset
    import ngsolve as ng

    nod = np.asarray(nodal, dtype=float).ravel()
    if nod.size != mesh.nv:
        raise ValueError("write_levelset_exodus: nodal field has %d entries, "
                         "mesh has %d vertices" % (nod.size, mesh.nv))
    if mesh.GetCurveOrder() >= 2:
        raise ValueError("write_levelset_exodus: curved meshes are not "
                         "supported")

    pts = np.array([list(v.point) for v in mesh.vertices], dtype=float)
    conn = []
    for el in mesh.Elements(ng.VOL):
        vs = [v.nr for v in el.vertices]
        if len(vs) != 4:
            raise NotImplementedError(
                "write_levelset_exodus: TET meshes only (the Cubit "
                "`create tri iso tet all` command cuts tets); got a %d-vertex "
                "element" % len(vs))
        conn.append([v + 1 for v in vs])          # Exodus is 1-based
    conn = np.asarray(conn, dtype=np.int32)

    d = Dataset(str(path), "w", format="NETCDF3_CLASSIC")
    try:
        d.api_version = np.float32(4.98)
        d.version = np.float32(4.98)
        d.floating_point_word_size = np.int32(8)
        d.file_size = np.int32(1)
        d.title = "radia.topopt_cad level set"
        d.createDimension("len_string", 33)
        d.createDimension("len_line", 81)
        d.createDimension("four", 4)
        d.createDimension("time_step", None)
        d.createDimension("num_dim", 3)
        d.createDimension("num_nodes", pts.shape[0])
        d.createDimension("num_elem", conn.shape[0])
        d.createDimension("num_el_blk", 1)
        d.createDimension("num_el_in_blk1", conn.shape[0])
        d.createDimension("num_nod_per_el1", 4)
        d.createDimension("num_nod_var", 1)
        for i, name in enumerate("xyz"):
            v = d.createVariable(f"coord{name}", "f8", ("num_nodes",))
            v[:] = pts[:, i]
        c = d.createVariable("connect1", "i4",
                             ("num_el_in_blk1", "num_nod_per_el1"))
        c.elem_type = "TETRA"
        c[:, :] = conn
        st = d.createVariable("eb_status", "i4", ("num_el_blk",))
        st[:] = [1]
        pr = d.createVariable("eb_prop1", "i4", ("num_el_blk",))
        pr.setncattr("name", "ID")
        pr[:] = [1]
        tw = d.createVariable("time_whole", "f8", ("time_step",))
        tw[0] = 0.0
        nv = d.createVariable("name_nod_var", "S1",
                              ("num_nod_var", "len_string"))
        nv[0, :] = np.array(list(varname.ljust(33, "\x00")), dtype="S1")
        vv = d.createVariable("vals_nod_var1", "f8",
                              ("time_step", "num_nodes"))
        vv[0, :] = nod - float(level)
    finally:
        d.close()
    return {"path": str(path), "n_nodes": int(pts.shape[0]),
            "n_tets": int(conn.shape[0]), "varname": varname,
            "level_shift": float(level)}


# ----------------------------------------------------------------------
# Pure-Python iso-surface route: resample -> marching cubes -> clean STL
# ----------------------------------------------------------------------

def iso_stl_from_grid(mesh, nodal, out_stl, *, level=0.5, resolution=64,
                      smooth_iterations=0, cutoff_factor=1.2,
                      target_faces=0):
    """Marching-cubes STL of the nodal density at ``level``.

    The nodal field is resampled onto a regular grid by nearest-VERTEX
    lookup; grid points farther than ``cutoff_factor`` x (median edge
    length) from every mesh vertex read as void, which closes the surface
    at (approximately) the design-domain boundary.  This route therefore
    blurs the domain boundary by O(grid spacing + h) -- the Cubit
    ``create tri iso`` route via :func:`write_levelset_exodus` is the
    mesh-native alternative without resampling.

    The raw marching-cubes surface is welded, degenerate faces are dropped,
    normals are made outward-consistent, and an optional Taubin smoothing
    pass (``smooth_iterations`` <= 5; volume drift is measured and
    reported) is applied.  Raises if the result is not watertight.

    ``target_faces > 0`` additionally quadric-decimates the surface to
    about that face count (requires ``fast-simplification``).  This is
    the lever that controls the DOWNSTREAM mesh size: Cubit meshes
    faceted (STL) geometry to the facets, so a marching-cubes surface
    with tens of thousands of facets forces tens of thousands of solver
    elements regardless of the requested element size (measured: a
    2996-facet body meshed to 27991 tets at size 0.06; an undecimated
    resolution-96 design surface forced 47k tets).  Decimation happens
    BEFORE the watertightness gate and its volume drift is reported.

    Returns a dict with counts, volumes, and the smoothing/decimation
    drifts.
    """
    _require("skimage", "scikit-image")
    trimesh = _require("trimesh", "trimesh")
    from skimage import measure
    from scipy.spatial import cKDTree
    import ngsolve as ng
    from trimesh import smoothing as tri_smoothing

    if smooth_iterations < 0 or smooth_iterations > 5:
        raise ValueError(
            "iso_stl_from_grid: smooth_iterations must be in [0, 5] -- "
            "10 Taubin iterations shrank a test body by 12.8 % (measured "
            "2026-08-08); keep smoothing light and report the drift")

    nod = np.asarray(nodal, dtype=float).ravel()
    if nod.size != mesh.nv:
        raise ValueError("iso_stl_from_grid: nodal field has %d entries, "
                         "mesh has %d vertices" % (nod.size, mesh.nv))

    pts = np.array([list(v.point) for v in mesh.vertices], dtype=float)
    tree = cKDTree(pts)

    # median edge length from a sample of elements (straight tets)
    edge_lengths = []
    for el in mesh.Elements(ng.VOL):
        vs = [pts[v.nr] for v in el.vertices]
        for a in range(len(vs)):
            for b in range(a + 1, len(vs)):
                edge_lengths.append(np.linalg.norm(vs[a] - vs[b]))
        if len(edge_lengths) > 6000:
            break
    h_med = float(np.median(edge_lengths))
    cutoff = cutoff_factor * h_med

    lo = pts.min(axis=0)
    hi = pts.max(axis=0)
    span = hi - lo
    pad = 2.0 * span.max() / (resolution - 1)
    axes = [np.linspace(lo[i] - pad, hi[i] + pad, resolution)
            for i in range(3)]
    X, Y, Z = np.meshgrid(*axes, indexing="ij")
    grid_pts = np.stack([X.ravel(), Y.ravel(), Z.ravel()], axis=1)

    # Linear (Delaunay) interpolation of the nodal field.  Nearest-vertex
    # lookup is NOT usable here: it produces a piecewise-constant Voronoi
    # field whose cell walls explode into degenerate marching-cubes
    # triangles (measured: 43982-face non-watertight surface on a 194-tet
    # design).  Outside the point cloud's hull LinearND returns NaN ->
    # void; the distance cutoff additionally voids hull regions that are
    # not covered by the (possibly non-convex) mesh.
    from scipy.interpolate import LinearNDInterpolator
    from scipy import ndimage
    interp = LinearNDInterpolator(pts, nod)
    vals = interp(grid_pts)
    dist, _ = tree.query(grid_pts)
    vals[~np.isfinite(vals)] = 0.0
    vals[dist > cutoff] = 0.0
    field = vals.reshape(X.shape)
    # one-cell Gaussian blur turns the void cutoff jump into a smooth
    # transition marching cubes can triangulate cleanly; this is part of
    # the documented O(grid spacing) boundary blur of this route
    field = ndimage.gaussian_filter(field, sigma=1.0)

    spacing = tuple(float(a[1] - a[0]) for a in axes)
    if not (field.min() < level < field.max()):
        raise ValueError(
            "iso_stl_from_grid: level %g is outside the sampled field range "
            "[%g, %g] -- nothing to extract" % (level, field.min(),
                                                field.max()))
    verts, faces, _, _ = measure.marching_cubes(field, level=level,
                                                spacing=spacing)
    verts += np.array([axes[0][0], axes[1][0], axes[2][0]])

    m = trimesh.Trimesh(vertices=verts, faces=faces, process=True)
    m.merge_vertices()
    m.update_faces(m.nondegenerate_faces(height=1e-9 * span.max()))
    m.remove_unreferenced_vertices()
    m.fix_normals()
    vol_raw = float(abs(m.volume))
    if smooth_iterations:
        tri_smoothing.filter_taubin(m, lamb=0.5, nu=-0.53,
                                    iterations=smooth_iterations)
        m.update_faces(m.nondegenerate_faces(height=1e-9 * span.max()))
        m.remove_unreferenced_vertices()
        m.fix_normals()
    vol_smooth = float(abs(m.volume))
    if target_faces:
        if target_faces < 100:
            raise ValueError("iso_stl_from_grid: target_faces must be 0 "
                             "(off) or >= 100")
        _require("fast_simplification", "fast-simplification")
        m = m.simplify_quadric_decimation(face_count=int(target_faces))
        m.merge_vertices()
        m.remove_unreferenced_vertices()
        m.fix_normals()
    if not m.is_watertight:
        raise RuntimeError(
            "iso_stl_from_grid: extracted surface is not watertight "
            "(%d faces); raise `resolution`, lower `target_faces` less "
            "aggressively, or inspect the density" % len(m.faces))
    vol = float(abs(m.volume))
    m.export(str(out_stl))
    return {"path": str(out_stl), "n_faces": int(len(m.faces)),
            "n_vertices": int(len(m.vertices)),
            "volume": vol, "volume_before_smoothing": vol_raw,
            "smoothing_volume_drift": (vol_smooth - vol_raw) / vol_raw
                                      if vol_raw else 0.0,
            "decimation_volume_drift": (vol - vol_smooth) / vol_smooth
                                       if (target_faces and vol_smooth)
                                       else 0.0,
            "grid_resolution": int(resolution),
            "median_edge_length": h_med, "cutoff": cutoff,
            "watertight": True}
