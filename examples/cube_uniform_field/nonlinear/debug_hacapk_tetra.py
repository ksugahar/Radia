#!/usr/bin/env python
"""Debug HACApK for tetrahedra - check if H-matrix builds correctly."""

import sys
import os

_src_path = os.path.join(os.path.dirname(__file__), '../../../src/radia')
sys.path.insert(0, _src_path)

import numpy as np
import radia as rad

MU_0 = 4 * np.pi * 1e-7
CUBE_HALF = 0.5
H_EXT = 50000.0
B_EXT = MU_0 * H_EXT

BH_DATA = [
    [0.0, 0.0],
    [100.0, 0.1],
    [500.0, 0.8],
    [1000.0, 1.2],
    [5000.0, 1.7],
    [50000.0, 2.0],
]

def main():
    try:
        from netgen.occ import Box, Pnt, OCCGeometry
        from ngsolve import Mesh
        from netgen_mesh_import import netgen_mesh_to_radia
    except ImportError as e:
        print('SKIP: %s' % e)
        return

    rad.FldUnits('m')
    rad.UtiDelAll()

    # Create a small mesh
    maxh = 0.5
    print('Creating tetrahedral mesh (maxh=%.2f)...' % maxh)

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

    print('Generated %d tetrahedral elements' % n_elements)
    print('Expected DOF: %d (3 per element)' % (n_elements * 3))

    # Apply material
    mat = rad.MatSatIsoTab(BH_DATA)
    rad.MatApl(cube, mat)

    # External field
    ext = rad.ObjBckg([0, 0, B_EXT])
    grp = rad.ObjCnt([cube, ext])

    # Check Radia object info
    try:
        info = rad.UtiDmp(cube)
        print('Object info: %s' % (info[:200] if len(info) > 200 else info))
    except Exception as e:
        print('UtiDmp failed: %s' % e)

    # Check interaction info
    print('Setting HACApK params...')
    try:
        rad.SetHACApKParams(1e-4, 10, 2.0)
        print('  SetHACApKParams: OK')
    except Exception as e:
        print('  SetHACApKParams failed: %s' % e)

    # Run Solve with Method 2
    print('Running Solve with Method 2 (HACApK)...')
    try:
        # Run 1 iteration at a time to see progress
        for i in range(5):
            result = rad.Solve(grp, 0.001, 1, 2)
            all_M = rad.ObjM(cube)
            M_list = [m[1] for m in all_M]
            M_avg_z = np.mean([m[2] for m in M_list])
            print('  Iter %d: result=%s, M_avg_z=%.1f' % (i+1, str(result), M_avg_z))
    except Exception as e:
        import traceback
        print('  Solve failed: %s' % e)
        traceback.print_exc()

    # Check H-matrix stats
    print('Getting HACApK stats...')
    try:
        stats = rad.GetHACApKStats()
        print('  stats: %s' % stats)
    except Exception as e:
        print('  GetHACApKStats failed: %s' % e)

    # Check magnetization
    print('Checking magnetization...')
    all_M = rad.ObjM(cube)
    print('  Number of M entries: %d' % len(all_M))
    if all_M:
        M_list = [m[1] for m in all_M]
        M_avg = [np.mean([m[0] for m in M_list]),
                 np.mean([m[1] for m in M_list]),
                 np.mean([m[2] for m in M_list])]
        print('  M_avg: [%.1f, %.1f, %.1f] A/m' % tuple(M_avg))

if __name__ == '__main__':
    main()
