"""NGSolve integration utilities for Radia.

The primary CoefficientFunction is ``radia.RadiaField`` (C++ class in
_radia_pybind.pyd).  This module re-exports it for backward compatibility
and adds a standalone ``create_voxel_cf()`` helper for trajectory
calculations where O(1) evaluation per step is needed.

Usage (primary -- CoefficientFunction):
    import radia as rad
    from ngsolve import *

    B_cf = rad.RadiaField(magnet, 'b')      # C++ CoefficientFunction
    gf = GridFunction(HDiv(mesh, order=2))
    gf.Set(B_cf)                             # direct integration-point eval

Usage (option -- VoxelCoefficient for trajectory):
    from radia.radia_ngsolve import create_voxel_cf

    B_voxel = create_voxel_cf(magnet, 'b', mesh=mesh, resolution=61)
    # O(1) trilinear interpolation per evaluation
"""

import numpy as np


def create_voxel_cf(radia_obj, field_type='b', mesh=None, bbox=None,
                    resolution=41):
    """Create a VoxelCoefficient CoefficientFunction from a Radia object.

    Pre-computes field values on a regular 3-D grid, then wraps them in
    ``ngsolve.VoxelCoefficient`` with trilinear interpolation.  Evaluation
    cost is O(1) per point -- ideal for particle trajectory integration.

    Parameters
    ----------
    radia_obj : int
        Radia object handle.
    field_type : str
        Field type: 'b', 'h', 'a', 'm', or 'phi'.
    mesh : ngsolve.Mesh, optional
        Used to determine bounding box automatically.
    bbox : list of [min, max] pairs, optional
        Explicit bounding box ``[[xmin,xmax], [ymin,ymax], [zmin,zmax]]``.
        Overrides ``mesh`` if given.
    resolution : int
        Number of voxels per dimension (default 41).

    Returns
    -------
    ngsolve.CoefficientFunction
        VoxelCoefficient (scalar) or CF(tuple) of VoxelCoefficients (vector).
    """
    import radia as rad
    from ngsolve import VoxelCoefficient, CF

    if bbox is None:
        if mesh is None:
            raise ValueError("Provide mesh or bbox")
        pmin, pmax = mesh.ngmesh.bounding_box
        pmin = [pmin[i] for i in range(3)]
        pmax = [pmax[i] for i in range(3)]
        margin = 0.01 * max(pmax[i] - pmin[i] for i in range(3))
        bbox = [[pmin[i] - margin, pmax[i] + margin] for i in range(3)]

    nx = ny = nz = resolution
    x = np.linspace(bbox[0][0], bbox[0][1], nx)
    y = np.linspace(bbox[1][0], bbox[1][1], ny)
    z = np.linspace(bbox[2][0], bbox[2][1], nz)

    xx, yy, zz = np.meshgrid(x, y, z, indexing='ij')
    points = np.column_stack([xx.ravel(), yy.ravel(), zz.ravel()])

    field = np.asarray(rad.Fld(radia_obj, field_type, points))

    start = (bbox[0][0], bbox[1][0], bbox[2][0])
    end = (bbox[0][1], bbox[1][1], bbox[2][1])

    if field_type == 'phi':
        data = field.reshape(nx, ny, nz).transpose(2, 1, 0).copy()
        return VoxelCoefficient(start, end, data, linear=True)
    else:
        cfs = []
        for comp in range(3):
            data = field[:, comp].reshape(nx, ny, nz)
            data = data.transpose(2, 1, 0).copy()
            cfs.append(VoxelCoefficient(start, end, data, linear=True))
        return CF(tuple(cfs))


def prepare_cache_hmatrix(cf, points, eps=1e-6):
    """Pre-cache a RadiaField CoefficientFunction at `points` via the O(N log N) HACApK
    ``_FieldEvalHMatrix`` (instead of the direct O(N_pts*N_src) rad.Fld inside ``cf.PrepareCache``), so a
    subsequent ``gf.Set(cf)`` hits the cache.  This is the path-B (magnet-container field -> GridFunction)
    acceleration for the RadiaField CF workflow; pass the FES integration points as `points`.

    cf     : a ``radia.RadiaField(obj, 'b'|'a')`` CoefficientFunction (flat container, no per-object
             transform; 'h'/'phi' are not yet in the H-matrix -> use cf.PrepareCache for those).
    points : (N,3) array-like of global evaluation points (the FES integration points gf.Set will use).
    eps    : ACA tolerance.

    CALLER wraps in ``with ngsolve.TaskManager():`` -- the H-matrix build/matvec parallelise under it.
    """
    import numpy as np
    import radia._radia_pybind as _rp

    ft = cf.field_type
    if ft not in ("b", "a"):
        raise ValueError("prepare_cache_hmatrix: field_type must be 'b' or 'a' (got %r); use "
                         "cf.PrepareCache for 'h'/'phi'." % ft)
    if getattr(cf, "use_transform", False):
        raise ValueError("prepare_cache_hmatrix: the H-matrix supports a FLAT global-coordinate container "
                         "only (cf.use_transform is True); use cf.PrepareCache.")
    pts = np.asarray(points, float).reshape(-1, 3)
    G = _rp._FieldEvalHMatrix(cf.radia_obj, pts.reshape(-1).tolist(), eps=eps, field_type=ft)
    x = [0.0] * (3 * G.n_obs()) + list(G.src_magnetization())
    vals = np.asarray(G.matvec(x), float)[:3 * G.n_obs()].reshape(-1, 3)
    cf.PrepareCacheFromValues(pts, vals)


# Backward compatibility: ``from radia_ngsolve import RadiaField``
# now returns the C++ CoefficientFunction from _radia_pybind.pyd.
try:
    from radia import RadiaField
except ImportError:
    pass
