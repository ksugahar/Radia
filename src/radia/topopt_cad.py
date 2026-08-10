"""Density-design -> CAD bridge: nodal level sets, Exodus hand-off, iso STL.

This module is the DESIGN-side half of the topology-optimization shape
regeneration loop (design record:
``docs/hdiv_vim/TOPOLOGY_OPTIMIZATION_ISOCHRONOUS.md``; the MESH-side half
lives on the radia-mcp cubit server, which turns the artifacts written here
into all-hex / tet solver meshes with closure and quality gates):

    per-element density  --(``nodal_from_element_density``)-->  P1 level set
        --> ``write_levelset_exodus``  (Cubit ``create tri iso`` route,
                                        mesh-native, no resampling)
        --> ``iso_stl_from_grid``      (standalone marching-cubes route)
        --> ``write_vfrac_exodus``     (Sculpt-native volume-fraction route:
                                        no STL at all -- Sculpt's interface
                                        reconstruction places the boundary
                                        sub-cell from the cell fractions;
                                        measured volume fidelity ~3e-4 vs
                                        the STL route's percent-level
                                        closure)

Both routes are deliberately gated: generated surfaces must remain watertight,
smoothing is capped and reports volume drift, and the downstream Cubit bridge
checks closure, inversion, and exported boundary faces before solver use.

Optional dependencies (each raises with install guidance when missing, per
the fail-fast policy): ``netCDF4`` (Exodus), ``scikit-image`` (marching
cubes), ``trimesh`` (STL hygiene), ``fast-simplification`` (quadric
decimation via ``target_faces``).
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

__all__ = [
    "iso_stl_from_grid",
    "nodal_from_element_density",
    "write_levelset_exodus",
    "write_vfrac_exodus",
]


def _require(module_name: str, pip_name: str):
    try:
        return __import__(module_name)
    except ImportError as exc:
        raise ImportError(
            f"radia.topopt_cad requires the optional dependency "
            f"'{module_name}' for this function; install it with "
            f"`pip install {pip_name}`") from exc


def _bounded_int(name, value, minimum, maximum):
    import operator

    try:
        parsed = operator.index(value)
    except TypeError as exc:
        raise ValueError(f"{name} must be an integer, got {value!r}") from exc
    if not minimum <= parsed <= maximum:
        raise ValueError(
            f"{name} must be in [{minimum}, {maximum}], got {parsed}")
    return parsed


def _finite_float(name, value, *, positive=False):
    try:
        parsed = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be a finite number") from exc
    if not np.isfinite(parsed) or (positive and parsed <= 0.0):
        qualifier = "finite and positive" if positive else "finite"
        raise ValueError(f"{name} must be {qualifier}, got {value!r}")
    return parsed


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
            f"nodal_from_element_density: density has {rho.size} entries, "
            f"mesh has {mesh.ne} volume elements")
    if not np.all(np.isfinite(rho)):
        raise ValueError("nodal_from_element_density: density has non-finite "
                         "entries")

    volumes = np.asarray(
        ng.Integrate(ng.CoefficientFunction(1.0), mesh, ng.VOL,
                     element_wise=True), dtype=float)
    if volumes.size != mesh.ne or not np.all(volumes > 0.0):
        minimum = volumes.min() if volumes.size else np.nan
        raise RuntimeError(
            "nodal_from_element_density: element volume integration "
            f"returned {volumes.size} entries / min {minimum:g} -- mesh is "
            "inconsistent")

    num = np.zeros(mesh.nv)
    den = np.zeros(mesh.nv)
    for el, rho_e, vol_e in zip(mesh.Elements(ng.VOL), rho, volumes):
        for v in el.vertices:
            num[v.nr] += rho_e * vol_e
            den[v.nr] += vol_e
    if not np.all(den > 0.0):
        isolated = int(np.sum(den <= 0.0))
        raise RuntimeError(
            f"nodal_from_element_density: {isolated} vertices have no "
            "incident volume element")
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
    nod = np.asarray(nodal, dtype=float).ravel()
    if mesh.dim != 3:
        raise ValueError("write_levelset_exodus: 3D meshes only")
    if nod.size != mesh.nv:
        raise ValueError(
            f"write_levelset_exodus: nodal field has {nod.size} entries, "
            f"mesh has {mesh.nv} vertices")
    if mesh.GetCurveOrder() >= 2:
        raise ValueError("write_levelset_exodus: curved meshes are not "
                         "supported")
    if not np.isfinite(nod).all():
        raise ValueError("write_levelset_exodus: nodal field has non-finite "
                         "entries")
    level = _finite_float("write_levelset_exodus: level", level)
    if (not isinstance(varname, str) or not varname
            or len(varname.encode("ascii", errors="ignore")) != len(varname)
            or len(varname) > 32
            or any(not 32 <= ord(ch) <= 126 for ch in varname)):
        raise ValueError(
            "write_levelset_exodus: varname must be 1-32 printable ASCII "
            "characters")

    _require("netCDF4", "netCDF4")
    import ngsolve as ng
    from netCDF4 import Dataset

    pts = np.array([list(v.point) for v in mesh.vertices], dtype=float)
    if pts.shape != (mesh.nv, 3) or not np.isfinite(pts).all():
        raise ValueError("write_levelset_exodus: mesh vertices must be finite "
                         "3D coordinates")
    conn = []
    for el in mesh.Elements(ng.VOL):
        vs = [v.nr for v in el.vertices]
        if len(vs) != 4:
            raise NotImplementedError(
                "write_levelset_exodus: TET meshes only (the Cubit "
                "`create tri iso tet all` command cuts tets); got a "
                f"{len(vs)}-vertex element")
        conn.append([v + 1 for v in vs])          # Exodus is 1-based
    if not conn:
        raise ValueError("write_levelset_exodus: mesh has no TET volume "
                         "elements")
    conn = np.asarray(conn, dtype=np.int32)

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

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
        vv[0, :] = nod - level
    finally:
        d.close()
    return {"path": str(path), "n_nodes": int(pts.shape[0]),
            "n_tets": int(conn.shape[0]), "varname": varname,
            "level_shift": level}


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
    with tens of thousands of facets forces a correspondingly large solver
    mesh regardless of the requested element size.  Decimation happens
    BEFORE the watertightness gate and its volume drift is reported.

    Returns a dict with counts, volumes, and the smoothing/decimation
    drifts.
    """
    if mesh.dim != 3:
        raise ValueError("iso_stl_from_grid: 3D meshes only")
    if mesh.GetCurveOrder() >= 2:
        raise ValueError("iso_stl_from_grid: curved meshes are not supported")
    resolution = _bounded_int(
        "iso_stl_from_grid: resolution", resolution, 16, 128)
    smooth_iterations = _bounded_int(
        "iso_stl_from_grid: smooth_iterations", smooth_iterations, 0, 5)
    target_faces = _bounded_int(
        "iso_stl_from_grid: target_faces", target_faces, 0, 1_000_000)
    if 0 < target_faces < 100:
        raise ValueError("iso_stl_from_grid: target_faces must be 0 "
                         "(off) or >= 100")
    cutoff_factor = _finite_float(
        "iso_stl_from_grid: cutoff_factor", cutoff_factor, positive=True)
    level = _finite_float("iso_stl_from_grid: level", level)
    nod = np.asarray(nodal, dtype=float).ravel()
    if nod.size != mesh.nv:
        raise ValueError(
            f"iso_stl_from_grid: nodal field has {nod.size} entries, "
            f"mesh has {mesh.nv} vertices")
    if not np.isfinite(nod).all():
        raise ValueError("iso_stl_from_grid: nodal field has non-finite "
                         "entries")

    _require("skimage", "scikit-image")
    trimesh = _require("trimesh", "trimesh")
    import ngsolve as ng
    from scipy.spatial import cKDTree
    from skimage import measure
    from trimesh import smoothing as tri_smoothing

    pts = np.array([list(v.point) for v in mesh.vertices], dtype=float)
    if pts.shape != (mesh.nv, 3) or not np.isfinite(pts).all():
        raise ValueError("iso_stl_from_grid: mesh vertices must be finite "
                         "3D coordinates")
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
    if not edge_lengths:
        raise ValueError("iso_stl_from_grid: mesh has no volume-element edges")
    h_med = float(np.median(edge_lengths))
    if not np.isfinite(h_med) or h_med <= 0.0:
        raise ValueError("iso_stl_from_grid: mesh has no finite positive "
                         "edge length")
    cutoff = cutoff_factor * h_med

    lo = pts.min(axis=0)
    hi = pts.max(axis=0)
    span = hi - lo
    if not np.all(np.isfinite(span)) or not np.all(span > 0.0):
        raise ValueError("iso_stl_from_grid: mesh bounding box is degenerate")
    pad = 2.0 * span.max() / (resolution - 1)
    axes = [np.linspace(lo[i] - pad, hi[i] + pad, resolution)
            for i in range(3)]
    X, Y, Z = np.meshgrid(*axes, indexing="ij")
    grid_pts = np.stack([X.ravel(), Y.ravel(), Z.ravel()], axis=1)

    # Linear (Delaunay) interpolation of the nodal field.  Nearest-vertex
    # lookup is NOT usable here: its piecewise-constant Voronoi field can
    # create degenerate marching-cubes triangles.  Outside the point cloud's
    # hull LinearND returns NaN ->
    # void; the distance cutoff additionally voids hull regions that are
    # not covered by the (possibly non-convex) mesh.
    from scipy import ndimage
    from scipy.interpolate import LinearNDInterpolator
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
            f"iso_stl_from_grid: level {level:g} is outside the sampled "
            f"field range [{field.min():g}, {field.max():g}] -- nothing to "
            "extract")
    verts, faces, _, _ = measure.marching_cubes(field, level=level,
                                                spacing=spacing)
    verts += np.array([axes[0][0], axes[1][0], axes[2][0]])

    m = trimesh.Trimesh(vertices=verts, faces=faces, process=True)
    m.merge_vertices()
    m.update_faces(m.nondegenerate_faces(height=1e-9 * span.max()))
    m.remove_unreferenced_vertices()
    m.fix_normals()
    vol_raw = float(abs(m.volume))
    if not np.isfinite(vol_raw) or vol_raw <= 0.0:
        raise RuntimeError("iso_stl_from_grid: extracted surface has zero or "
                           "non-finite volume")
    if smooth_iterations:
        tri_smoothing.filter_taubin(m, lamb=0.5, nu=-0.53,
                                    iterations=smooth_iterations)
        m.update_faces(m.nondegenerate_faces(height=1e-9 * span.max()))
        m.remove_unreferenced_vertices()
        m.fix_normals()
    vol_smooth = float(abs(m.volume))
    if not np.isfinite(vol_smooth) or vol_smooth <= 0.0:
        raise RuntimeError("iso_stl_from_grid: smoothing produced zero or "
                           "non-finite volume")
    if target_faces and target_faces < len(m.faces):
        _require("fast_simplification", "fast-simplification")
        m = m.simplify_quadric_decimation(face_count=target_faces)
        m.merge_vertices()
        m.remove_unreferenced_vertices()
        m.fix_normals()
    if not m.is_watertight:
        raise RuntimeError(
            "iso_stl_from_grid: extracted surface is not watertight "
            f"({len(m.faces)} faces); raise `resolution`, lower "
            "`target_faces` less aggressively, or inspect the density")
    vol = float(abs(m.volume))
    if not np.isfinite(vol) or vol <= 0.0:
        raise RuntimeError("iso_stl_from_grid: final surface has zero or "
                           "non-finite volume")
    out_stl = Path(out_stl)
    out_stl.parent.mkdir(parents=True, exist_ok=True)
    m.export(str(out_stl))
    return {"path": str(out_stl), "n_faces": len(m.faces),
            "n_vertices": len(m.vertices),
            "volume": vol, "volume_before_smoothing": vol_raw,
            "smoothing_volume_drift": (vol_smooth - vol_raw) / vol_raw
                                      if vol_raw else 0.0,
            "decimation_volume_drift": (vol - vol_smooth) / vol_smooth
                                       if (target_faces and vol_smooth)
                                       else 0.0,
            "grid_resolution": int(resolution),
            "median_edge_length": h_med, "cutoff": cutoff,
            "watertight": True}


# ----------------------------------------------------------------------
# Sculpt-native volume-fraction route: no surface extraction at all
# ----------------------------------------------------------------------

def write_vfrac_exodus(mesh, nodal, path, *, level=0.5, cells=48,
                       supersample=4, pad_cells=2, cutoff_factor=2.0):
    """Cartesian volume-fraction Exodus of the body ``{nodal >= level}``
    for Sculpt's ``--input_vfrac`` reader.

    Sculpt is natively a volume-fraction mesher ("all-hex parallel mesh
    generation from volume fractions"); the STL route is one of its front
    ends.  Feeding the design occupancy DIRECTLY skips marching cubes,
    STL welding, and watertightness entirely, and Sculpt's interface
    reconstruction places the boundary sub-cell: measured on a voxel
    sphere (24^3 grid), the meshed material volume tracks the input
    fractions to ~1e-4 relative, against -8 % shrink for the equivalent
    0/1 SPN input and the percent-level closure of the marching-cubes STL
    route at comparable resolution.

    The occupancy uses the SAME body definition as ``iso_stl_from_grid``
    (linear interpolation of the P1 nodal field over the design mesh
    vertices, void beyond ``cutoff_factor`` x median edge from any
    vertex, body = field >= ``level``), sampled at ``supersample``^3
    sub-cell centers per Cartesian cell, so the two routes regenerate the
    same design and stay A/B comparable.

    File layout (verified against Sculpt's own ``--volfrac_file`` dump,
    Coreform Cubit 2025.12): NETCDF3 Exodus, one HEX8 block on the
    Cartesian lattice in PHYSICAL coordinates, element fields ``VOID`` +
    ``MAT_1`` summing to 1, and the 22 global variables (grid bounds,
    intervals, material table, serial parallel bookkeeping) its reader
    requires.  Sculpt resolves the file from ``<base>.e.1.0``; ``path``
    may be either the base name or the full ``.e.1.0`` name.

    Args:
        mesh: straight (non-curved) 3D NGSolve design mesh.
        nodal: P1 nodal field (``nodal_from_element_density`` output).
        path: output base name or full ``.e.1.0`` path.
        level: body threshold (the iso contract's 0.5 default).
        cells: Cartesian cells along the LONGEST bbox axis (cells are
            cubic; other axes scale by extent).
        supersample: sub-samples per axis per cell for the occupancy.
        pad_cells: empty guard cells added around the bbox on each side.
        cutoff_factor: void beyond this multiple of the median mesh edge
            from the nearest design vertex (matches ``iso_stl_from_grid``).

    Returns:
        Dict with the resolved ``path``, per-axis cell counts, cell size,
        grid bounds, and the volume-fraction integral ``v_vfrac`` (the
        closure reference for the downstream mesh gate).
    """
    nc = _require("netCDF4", "netCDF4")
    import ngsolve as ng
    from scipy.interpolate import LinearNDInterpolator
    from scipy.spatial import cKDTree

    if mesh.dim != 3:
        raise ValueError("write_vfrac_exodus: 3D meshes only")
    if mesh.GetCurveOrder() > 1:
        raise ValueError("write_vfrac_exodus: curved meshes are not "
                         "supported")
    cells = _bounded_int("write_vfrac_exodus: cells", cells, 8, 512)
    supersample = _bounded_int("write_vfrac_exodus: supersample",
                               supersample, 1, 8)
    pad_cells = _bounded_int("write_vfrac_exodus: pad_cells", pad_cells,
                             0, 16)
    cutoff_factor = _finite_float("write_vfrac_exodus: cutoff_factor",
                                  cutoff_factor, positive=True)
    level = _finite_float("write_vfrac_exodus: level", level)
    nod = np.asarray(nodal, dtype=float).ravel()
    if nod.size != mesh.nv:
        raise ValueError(
            f"write_vfrac_exodus: nodal field has {nod.size} entries, "
            f"mesh has {mesh.nv} vertices")
    if not np.isfinite(nod).all():
        raise ValueError("write_vfrac_exodus: nodal field has non-finite "
                         "entries")

    pts = np.array([list(v.point) for v in mesh.vertices], dtype=float)
    if pts.shape != (mesh.nv, 3) or not np.isfinite(pts).all():
        raise ValueError("write_vfrac_exodus: mesh vertices must be "
                         "finite 3D coordinates")
    tree = cKDTree(pts)
    edge_lengths = []
    for el in mesh.Elements(ng.VOL):
        vs = [pts[v.nr] for v in el.vertices]
        for a in range(len(vs)):
            for b in range(a + 1, len(vs)):
                edge_lengths.append(np.linalg.norm(vs[a] - vs[b]))
        if len(edge_lengths) > 6000:
            break
    if not edge_lengths:
        raise ValueError("write_vfrac_exodus: mesh has no volume-element "
                         "edges")
    h_med = float(np.median(edge_lengths))
    if not np.isfinite(h_med) or h_med <= 0.0:
        raise ValueError("write_vfrac_exodus: mesh has no finite positive "
                         "edge length")
    cutoff = cutoff_factor * h_med

    lo = pts.min(axis=0)
    hi = pts.max(axis=0)
    span = hi - lo
    if not np.all(np.isfinite(span)) or not np.all(span > 0.0):
        raise ValueError("write_vfrac_exodus: mesh bounding box is "
                         "degenerate")
    h_c = float(span.max()) / cells
    nel = [int(np.ceil(span[i] / h_c)) + 2 * pad_cells for i in range(3)]
    gmin = [float(lo[i]) - pad_cells * h_c for i in range(3)]
    gmax = [gmin[i] + nel[i] * h_c for i in range(3)]

    # occupancy at supersampled sub-cell centers (i fastest, matching the
    # exodus lattice ordering below)
    ss = supersample
    frac = None
    axes_fine = [gmin[i] + (np.arange(nel[i] * ss) + 0.5) * (h_c / ss)
                 for i in range(3)]
    Xf, Yf, Zf = np.meshgrid(*axes_fine, indexing="ij")
    fine_pts = np.stack([Xf.ravel(), Yf.ravel(), Zf.ravel()], axis=1)
    interp = LinearNDInterpolator(pts, nod)
    vals = interp(fine_pts)
    dist, _ = tree.query(fine_pts)
    vals[~np.isfinite(vals)] = 0.0
    vals[dist > cutoff] = 0.0
    inside = (vals >= level).astype(float).reshape(
        nel[0], ss, nel[1], ss, nel[2], ss)
    frac = inside.mean(axis=(1, 3, 5))
    v_vfrac = float(frac.sum()) * h_c ** 3
    if v_vfrac <= 0.0:
        raise ValueError(
            "write_vfrac_exodus: the thresholded body has zero volume on "
            "the sampling grid -- inspect the density/level")

    out = Path(path)
    if not out.name.endswith(".e.1.0"):
        out = out.with_name(out.name + ".e.1.0")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.unlink(missing_ok=True)

    nx, ny, nz = nel
    nn = (nx + 1) * (ny + 1) * (nz + 1)
    ne = nx * ny * nz
    ds = nc.Dataset(str(out), "w", format="NETCDF3_64BIT_OFFSET")
    try:
        ds.setncatts({"api_version": np.float32(8.25),
                      "version": np.float32(8.25),
                      "floating_point_word_size": np.int32(8),
                      "file_size": np.int32(1),
                      "title": "radia.topopt_cad volume fractions"})
        for name, size in (("len_name", 256), ("time_step", None),
                           ("num_dim", 3), ("num_nodes", nn),
                           ("num_elem", ne), ("num_el_blk", 1),
                           ("num_el_in_blk1", ne), ("num_nod_per_el1", 8),
                           ("num_elem_var", 2), ("num_glo_var", 22)):
            ds.createDimension(name, size)

        def chararr(var, strings, width=256):
            arr = np.zeros((len(strings), width), dtype="S1")
            for row, text in enumerate(strings):
                for col, ch in enumerate(text):
                    arr[row, col] = ch.encode()
            var[:] = arr

        axis = [gmin[i] + np.arange(nel[i] + 1) * h_c for i in range(3)]
        Xn, Yn, Zn = np.meshgrid(*axis, indexing="ij")
        flat = lambda A: A.ravel(order="F")   # noqa: E731  (i fastest)
        for nm, arr in (("coordx", Xn), ("coordy", Yn), ("coordz", Zn)):
            var = ds.createVariable(nm, "f8", ("num_nodes",))
            var[:] = flat(arr)
        node_id = np.arange(nn).reshape(nx + 1, ny + 1, nz + 1, order="F")
        conn = np.empty((ne, 8), dtype=np.int32)
        e = 0
        for k in range(nz):
            for j in range(ny):
                for i in range(nx):
                    conn[e] = (node_id[i, j, k], node_id[i + 1, j, k],
                               node_id[i + 1, j + 1, k], node_id[i, j + 1, k],
                               node_id[i, j, k + 1], node_id[i + 1, j, k + 1],
                               node_id[i + 1, j + 1, k + 1],
                               node_id[i, j + 1, k + 1])
                    e += 1
        var = ds.createVariable("connect1", "i4",
                                ("num_el_in_blk1", "num_nod_per_el1"))
        var.elem_type = "HEX8"
        var[:] = conn + 1
        var = ds.createVariable("eb_prop1", "i4", ("num_el_blk",))
        var.setncattr("name", "ID")
        var[:] = [1]
        var = ds.createVariable("eb_status", "i4", ("num_el_blk",))
        var[:] = [1]
        chararr(ds.createVariable("eb_names", "S1",
                                  ("num_el_blk", "len_name")), [""])
        chararr(ds.createVariable("coor_names", "S1",
                                  ("num_dim", "len_name")), ["x", "y", "z"])
        var = ds.createVariable("time_whole", "f8", ("time_step",))
        var[:] = [0.0]
        chararr(ds.createVariable("name_glo_var", "S1",
                                  ("num_glo_var", "len_name")),
                ["num_procs", "rank",
                 "neighbor_plus_i", "neighbor_negative_i",
                 "neighbor_plus_j", "neighbor_negative_j",
                 "neighbor_plus_k", "neighbor_negative_k",
                 "xmin", "ymin", "zmin", "xmax", "ymax", "zmax",
                 "xint", "yint", "zint", "gxint", "gyint", "gzint",
                 "num_mats", "mat[00]"])
        var = ds.createVariable("vals_glo_var", "f8",
                                ("time_step", "num_glo_var"))
        var[:] = [[1, 0, -1, -1, -1, -1, -1, -1,
                   gmin[0], gmin[1], gmin[2], gmax[0], gmax[1], gmax[2],
                   nx, ny, nz, nx, ny, nz, 1, 1]]
        chararr(ds.createVariable("name_elem_var", "S1",
                                  ("num_elem_var", "len_name")),
                ["VOID", "MAT_1"])
        var = ds.createVariable("elem_var_tab", "i4",
                                ("num_el_blk", "num_elem_var"))
        var[:] = [[1, 1]]
        mat = flat(frac)
        var = ds.createVariable("vals_elem_var1eb1", "f8",
                                ("time_step", "num_el_in_blk1"))
        var[:] = (1.0 - mat)[None, :]
        var = ds.createVariable("vals_elem_var2eb1", "f8",
                                ("time_step", "num_el_in_blk1"))
        var[:] = mat[None, :]
    finally:
        ds.close()
    return {"path": str(out), "nel": [nx, ny, nz], "cell_size": h_c,
            "bounds_min": gmin, "bounds_max": gmax,
            "v_vfrac": v_vfrac, "supersample": ss,
            "median_edge_length": h_med, "cutoff": cutoff,
            "occupied_cells": int((frac > 0).sum())}
