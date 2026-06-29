#!/usr/bin/env python
"""
Compare H-matrix statistics between ELF and Radia.
This script runs both solvers and compares leaf-by-leaf rank information.
"""

import sys
import os
import json

# Add both paths
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/radia'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../build-msvc/Release'))

import radia as rad
from netgen.occ import Box, Pnt, OCCGeometry
from ngsolve import Mesh, VOL
from ngsolve import TaskManager

# Parameters
CUBE_SIZE = 1.0
MAXH = 0.20
H_EXT = 200000.0

# B-H curve data
BH_DATA = [
    [0.0, 0.0],
    [50.0, 0.1],
    [100.0, 0.35],
    [200.0, 0.65],
    [500.0, 1.0],
    [1000.0, 1.2],
    [2000.0, 1.35],
    [5000.0, 1.5],
    [10000.0, 1.6],
    [50000.0, 1.8],
    [100000.0, 1.9],
    [200000.0, 2.0],
]

def create_tetra_mesh():
    """Create tetrahedral mesh using Netgen"""

    half = CUBE_SIZE / 2
    box = Box(Pnt(-half, -half, -half), Pnt(half, half, half))
    box.mat('magnetic')
    geo = OCCGeometry(box)
    with TaskManager():
        mesh = Mesh(geo.GenerateMesh(maxh=MAXH))

        # Extract vertices
        vertices = []
        for v in mesh.vertices:
            pt = v.point
            vertices.append([pt[0], pt[1], pt[2]])

        # Create tetrahedra
        container = rad.ObjCnt([])
        for el in mesh.Elements(VOL):
            v_indices = [v.nr for v in el.vertices]
            tet_verts = [vertices[i] for i in v_indices]
            tet = rad.ObjTetrahedron(tet_verts, [0, 0, 0])
            rad.ObjAddToCnt(container, [tet])

        return container, len(list(mesh.Elements(VOL)))

def run_radia_hacapk(eps):
    """Run Radia HACApK solver"""
    rad.UtiDelAll()

    container, n_elem = create_tetra_mesh()

    # Apply material
    mat = rad.MatSatIsoTab(BH_DATA)
    rad.MatApl(container, mat)

    # Apply external field
    rad.ObjAddToCnt(container, [rad.ObjBckg(lambda p: [0, 0, H_EXT])])

    # Set HACApK parameters
    rad.SolverConfig(hacapk_eps=eps, hacapk_leaf=10, hacapk_eta=2.0)

    # Solve with HACApK
    result = rad.Solve(container, 0.001, 100, 2)  # Method 2 = HACApK

    # Get statistics
    config = rad.GetSolverConfig()
    stats = config.get('hacapk_stats')

    return {
        'n_elem': n_elem,
        'eps': eps,
        'nonl_iter': result[3] if len(result) > 3 else 0,
        'M_avg_z': rad.FldIntM(container, 'mxmymz')[2] / (CUBE_SIZE ** 3),
        'hmatrix': stats
    }

def main():
    import sys
    sys.stdout.flush()
    print("=" * 70, flush=True)
    print("H-matrix Statistics Comparison", flush=True)
    print("=" * 70, flush=True)

    # Test with different eps values
    eps_values = [1e-4, 1e-5, 1e-6]

    for eps in eps_values:
        print(f"\n--- eps = {eps} ---", flush=True)
        try:
            result = run_radia_hacapk(eps)
        except Exception as e:
            print(f"Error: {e}", flush=True)
            import traceback
            traceback.print_exc()
            continue

        print(f"Elements: {result['n_elem']}")
        print(f"Nonlinear iterations: {result['nonl_iter']}")
        print(f"M_avg_z: {result['M_avg_z']:.0f} A/m")

        if result['hmatrix']:
            hm = result['hmatrix']
            print(f"H-matrix stats:")
            print(f"  n_lowrank: {hm.get('n_lowrank', 'N/A')}")
            print(f"  n_dense: {hm.get('n_dense', 'N/A')}")
            print(f"  max_rank: {hm.get('max_rank', 'N/A')}")
            print(f"  nlf: {hm.get('nlf', 'N/A')}")

if __name__ == '__main__':
    main()
