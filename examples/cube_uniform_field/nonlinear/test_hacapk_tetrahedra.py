#!/usr/bin/env python
"""
Test HACApK solver for tetrahedra (3 DOF per element).

Quick test to verify HACApK H-matrix works with tetrahedral meshes.
"""

import sys
import os
import time

# Path setup
_src_path = os.path.join(os.path.dirname(__file__), '../../../src/radia')
sys.path.insert(0, _src_path)

import numpy as np
import radia as rad

# Physical constants
MU_0 = 4 * np.pi * 1e-7  # T/(A/m)

# Problem parameters
CUBE_SIZE = 1.0      # 1.0 m cube
CUBE_HALF = 0.5      # half size
H_EXT = 50000.0      # External field (A/m)
B_EXT = MU_0 * H_EXT  # External B field (T)

# B-H curve data
BH_DATA = [
    [0.0, 0.0],
    [100.0, 0.1],
    [200.0, 0.3],
    [500.0, 0.8],
    [1000.0, 1.2],
    [2000.0, 1.5],
    [5000.0, 1.7],
    [10000.0, 1.8],
    [50000.0, 2.0],
    [100000.0, 2.1],
]


def test_hacapk_tetrahedra():
    """Test HACApK solver with tetrahedral mesh."""
    try:
        from netgen.occ import Box, Pnt, OCCGeometry
        from ngsolve import Mesh
        from netgen_mesh_import import netgen_mesh_to_radia
    except ImportError as e:
        print('SKIP: %s' % e)
        return

    rad.FldUnits('m')
    rad.UtiDelAll()

    # Create mesh with Netgen (small for quick test)
    maxh = 0.5
    print('Creating tetrahedral mesh (maxh=%.2f)...' % maxh)

    t_mesh_start = time.time()
    cube_solid = Box(Pnt(-CUBE_HALF, -CUBE_HALF, -CUBE_HALF),
                      Pnt(CUBE_HALF, CUBE_HALF, CUBE_HALF))
    cube_solid.mat('magnetic')
    geo = OCCGeometry(cube_solid)

    ngmesh = geo.GenerateMesh(maxh=maxh)
    mesh = Mesh(ngmesh)
    n_elements = mesh.ne

    cube = netgen_mesh_to_radia(mesh,
                                 material={'magnetization': [0, 0, 0]},
                                 units='m',
                                 material_filter='magnetic',
                                 verbose=False)
    t_mesh = time.time() - t_mesh_start

    print('Generated %d tetrahedral elements' % n_elements)
    print('DOF: %d (3 per element)' % (n_elements * 3))

    # Apply nonlinear material
    mat = rad.MatSatIsoTab(BH_DATA)
    rad.MatApl(cube, mat)

    # External field
    ext = rad.ObjBckg([0, 0, B_EXT])
    grp = rad.ObjCnt([cube, ext])

    # Test LU first (for reference)
    print('\n--- Testing LU solver ---')
    rad.UtiDelAll()

    cube_solid = Box(Pnt(-CUBE_HALF, -CUBE_HALF, -CUBE_HALF),
                      Pnt(CUBE_HALF, CUBE_HALF, CUBE_HALF))
    cube_solid.mat('magnetic')
    geo = OCCGeometry(cube_solid)
    ngmesh = geo.GenerateMesh(maxh=maxh)
    mesh = Mesh(ngmesh)

    cube = netgen_mesh_to_radia(mesh,
                                 material={'magnetization': [0, 0, 0]},
                                 units='m',
                                 material_filter='magnetic',
                                 verbose=False)
    mat = rad.MatSatIsoTab(BH_DATA)
    rad.MatApl(cube, mat)
    ext = rad.ObjBckg([0, 0, B_EXT])
    grp = rad.ObjCnt([cube, ext])

    t_lu_start = time.time()
    result = rad.Solve(grp, 0.001, 100, 0)  # Method 0 = LU
    t_lu = time.time() - t_lu_start

    all_M = rad.ObjM(cube)
    M_list = [m[1] for m in all_M]
    M_avg_z_lu = np.mean([m[2] for m in M_list])

    print('LU time: %.3f s' % t_lu)
    print('LU M_avg_z: %.0f A/m' % M_avg_z_lu)

    # Test HACApK
    print('\n--- Testing HACApK solver (Method 2) ---')
    rad.UtiDelAll()

    cube_solid = Box(Pnt(-CUBE_HALF, -CUBE_HALF, -CUBE_HALF),
                      Pnt(CUBE_HALF, CUBE_HALF, CUBE_HALF))
    cube_solid.mat('magnetic')
    geo = OCCGeometry(cube_solid)
    ngmesh = geo.GenerateMesh(maxh=maxh)
    mesh = Mesh(ngmesh)

    cube = netgen_mesh_to_radia(mesh,
                                 material={'magnetization': [0, 0, 0]},
                                 units='m',
                                 material_filter='magnetic',
                                 verbose=False)
    mat = rad.MatSatIsoTab(BH_DATA)
    rad.MatApl(cube, mat)
    ext = rad.ObjBckg([0, 0, B_EXT])
    grp = rad.ObjCnt([cube, ext])

    try:
        # Set HACApK parameters
        rad.SetHACApKParams(1e-4, 10, 2.0)
        print('HACApK parameters set: eps=1e-4, leaf_size=10, eta=2.0')
    except AttributeError:
        print('WARNING: SetHACApKParams not available')

    t_hacapk_start = time.time()
    result = rad.Solve(grp, 0.001, 100, 2)  # Method 2 = HACApK
    t_hacapk = time.time() - t_hacapk_start

    all_M = rad.ObjM(cube)
    M_list = [m[1] for m in all_M]
    M_avg_z_hacapk = np.mean([m[2] for m in M_list])

    print('HACApK time: %.3f s' % t_hacapk)
    print('HACApK M_avg_z: %.0f A/m' % M_avg_z_hacapk)

    # Get H-matrix stats
    try:
        stats = rad.GetHACApKStats()
        if stats:
            print('H-matrix stats:')
            print('  n_lowrank: %d' % stats.get('n_lowrank', 0))
            print('  n_dense: %d' % stats.get('n_dense', 0))
            print('  max_rank: %d' % stats.get('max_rank', 0))
            print('  compression: %.4f' % stats.get('compression', 1.0))
    except AttributeError:
        print('GetHACApKStats not available')

    # Compare results
    print('\n--- Comparison ---')
    diff = abs(M_avg_z_lu - M_avg_z_hacapk) / max(abs(M_avg_z_lu), 1)
    print('M_avg_z difference: %.2f%%' % (diff * 100))

    if diff < 0.01:
        print('OK: HACApK result matches LU (within 1%%)')
    else:
        print('WARNING: Results differ significantly!')

    print('\nTest complete.')


if __name__ == '__main__':
    test_hacapk_tetrahedra()
