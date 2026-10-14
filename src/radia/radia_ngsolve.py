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


# Backward compatibility: ``from radia_ngsolve import RadiaField``
# now returns the C++ CoefficientFunction from _radia_pybind.pyd.
try:
    from radia import RadiaField
except ImportError:
    pass
