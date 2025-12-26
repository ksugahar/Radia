#!/usr/bin/env python
"""
Test Mixed Hexahedra + Tetrahedra Elements

This test creates a simple case with both hexahedral (6DOF) and tetrahedral (3DOF)
elements to verify the VariableDOF solver correctly handles mixed element types.

Setup:
- 1 hexahedral cube (center)
- Tetrahedral mesh surrounding it (or vice versa)

For simplicity, we create separate hex and tetra objects and combine them.
"""

import sys
import os
import time

# Add paths for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../src/radia'))

import radia as rad

# Material: nonlinear isotropic (steel B-H curve)
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

MU_0 = 4 * 3.14159265358979 * 1e-7
H_EXT = 200000.0  # A/m (same as benchmarks)
B_ext = MU_0 * H_EXT

# Face indices (1-indexed for Radia)
HEX_FACES = [
    [1, 4, 3, 2],  # bottom (-z)
    [5, 6, 7, 8],  # top (+z)
    [1, 2, 6, 5],  # front (-y)
    [3, 4, 8, 7],  # back (+y)
    [1, 5, 8, 4],  # left (-x)
    [2, 3, 7, 6],  # right (+x)
]

TETRA_FACES = [
    [1, 3, 2],  # base (bottom)
    [1, 2, 4],  # front
    [2, 3, 4],  # right
    [3, 1, 4],  # left
]


def create_hex_element(center, size):
    """Create a hexahedral element (6DOF MSC)."""
    cx, cy, cz = center
    dx, dy, dz = size

    vertices = [
        [cx - dx/2, cy - dy/2, cz - dz/2],  # 1
        [cx + dx/2, cy - dy/2, cz - dz/2],  # 2
        [cx + dx/2, cy + dy/2, cz - dz/2],  # 3
        [cx - dx/2, cy + dy/2, cz - dz/2],  # 4
        [cx - dx/2, cy - dy/2, cz + dz/2],  # 5
        [cx + dx/2, cy - dy/2, cz + dz/2],  # 6
        [cx + dx/2, cy + dy/2, cz + dz/2],  # 7
        [cx - dx/2, cy + dy/2, cz + dz/2],  # 8
    ]

    return rad.ObjPolyhdr(vertices, HEX_FACES, [0, 0, 0])


def create_tetra_element(vertices_list):
    """Create a tetrahedral element (3DOF MMM)."""
    return rad.ObjPolyhdr(vertices_list, TETRA_FACES, [0, 0, 0])


def test_mixed_elements(solver_type='lu'):
    """Test mixed hex + tetra elements."""
    print('=' * 70)
    print('MIXED ELEMENT TEST: %s solver' % solver_type.upper())
    print('=' * 70)

    rad.FldUnits('m')
    rad.UtiDelAll()

    # Create hexahedral elements (2x2x2 = 8 hexahedra in center)
    hex_objs = []
    hex_size = 0.2  # Each hex is 0.2m on a side
    for iz in range(2):
        for iy in range(2):
            for ix in range(2):
                cx = -0.1 + ix * hex_size
                cy = -0.1 + iy * hex_size
                cz = -0.1 + iz * hex_size
                obj = create_hex_element([cx, cy, cz], [hex_size, hex_size, hex_size])
                hex_objs.append(obj)

    # Create tetrahedral elements (4 tetrahedra around the hex core)
    # Place tetrahedra at the corners extending outward
    tetra_objs = []

    # Tetrahedra at +x side
    tetra_objs.append(create_tetra_element([
        [0.2, -0.2, -0.2],  # 1
        [0.4, -0.2, -0.2],  # 2
        [0.3, 0.0, -0.1],   # 3
        [0.3, -0.1, 0.1],   # 4
    ]))

    # Tetrahedra at -x side
    tetra_objs.append(create_tetra_element([
        [-0.2, -0.2, -0.2],  # 1
        [-0.4, -0.2, -0.2],  # 2
        [-0.3, 0.0, -0.1],   # 3
        [-0.3, -0.1, 0.1],   # 4
    ]))

    # Tetrahedra at +y side
    tetra_objs.append(create_tetra_element([
        [-0.2, 0.2, -0.2],   # 1
        [0.2, 0.2, -0.2],    # 2
        [0.0, 0.4, -0.1],    # 3
        [0.0, 0.3, 0.1],     # 4
    ]))

    # Tetrahedra at -y side
    tetra_objs.append(create_tetra_element([
        [-0.2, -0.2, -0.2],  # 1
        [0.2, -0.2, -0.2],   # 2
        [0.0, -0.4, -0.1],   # 3
        [0.0, -0.3, 0.1],    # 4
    ]))

    n_hex = len(hex_objs)
    n_tetra = len(tetra_objs)
    print('Created %d hexahedral elements (6DOF each)' % n_hex)
    print('Created %d tetrahedral elements (3DOF each)' % n_tetra)
    print('Total elements: %d' % (n_hex + n_tetra))
    print('Total DOF: %d' % (n_hex * 6 + n_tetra * 3))

    # Combine all elements
    all_objs = hex_objs + tetra_objs
    container = rad.ObjCnt(all_objs)

    # Apply nonlinear material
    mat = rad.MatSatIsoTab(BH_DATA)
    rad.MatApl(container, mat)

    # Add background field
    ext = rad.ObjBckg([0, 0, B_ext])
    grp = rad.ObjCnt([container, ext])

    # Set solver method
    solver_methods = {'lu': 0, 'bicgstab': 1, 'hacapk': 2}
    method = solver_methods.get(solver_type, 0)

    if solver_type == 'hacapk':
        rad.SetHACApKParams(1e-4, 10, 2.0)
        print('H-matrix: Enabled (eps=1e-4, leaf_size=10, eta=2.0)')

    # Solve
    print('Solving with method %d (%s)...' % (method, solver_type.upper()))
    t_start = time.time()
    result = rad.Solve(grp, 0.001, 100, method)
    t_solve = time.time() - t_start

    # Get statistics
    stats = rad.GetSolveStats()

    print('')
    print('Results:')
    print('  Solve time:      %.3f s' % t_solve)
    print('  Nonl iterations: %d' % result[3])
    print('  Converged:       %s' % ('Yes' if result[1] < 0.001 else 'No'))
    print('  Residual:        %.6f' % result[1])

    # Compute average magnetization in center region
    # Sample at center of hex region
    B_center = rad.Fld(container, 'b', [0, 0, 0])
    print('  B at center:     [%.4f, %.4f, %.4f] T' % tuple(B_center))

    # Timing breakdown
    print('')
    print('Timing breakdown:')
    print('  t_matrix_build:  %.4f s' % stats.get('t_matrix_build', 0))
    print('  t_lu_decomp:     %.4f s' % stats.get('t_lu_decomp', 0))
    print('  t_linear_solve:  %.4f s' % stats.get('t_linear_solve', 0))

    if solver_type == 'hacapk':
        print('  t_hmatrix_build: %.4f s' % stats.get('t_hmatrix_build', 0))

    return {
        'converged': result[1] < 0.001,
        'residual': result[1],
        'iterations': result[3],
        't_solve': t_solve,
        'B_center': B_center,
        'n_hex': n_hex,
        'n_tetra': n_tetra,
    }


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Test mixed hex+tetra elements')
    parser.add_argument('--lu', action='store_true', help='Use LU solver')
    parser.add_argument('--bicgstab', action='store_true', help='Use BiCGSTAB solver')
    parser.add_argument('--hacapk', action='store_true', help='Use HACApK solver')
    parser.add_argument('--all', action='store_true', help='Test all solvers')
    args = parser.parse_args()

    solvers = []
    if args.lu:
        solvers.append('lu')
    if args.bicgstab:
        solvers.append('bicgstab')
    if args.hacapk:
        solvers.append('hacapk')
    if args.all or not solvers:
        solvers = ['lu', 'bicgstab', 'hacapk']

    results = {}
    for solver in solvers:
        results[solver] = test_mixed_elements(solver)
        print('')

    # Summary
    print('=' * 70)
    print('SUMMARY')
    print('=' * 70)
    print('')
    print('%-10s %8s %8s %8s %s' % ('Solver', 'Iter', 'Time(s)', 'Conv', 'B_center_z'))
    print('-' * 60)
    for solver, res in results.items():
        if res:
            print('%-10s %8d %8.3f %8s %.4f T' % (
                solver.upper(),
                res['iterations'],
                res['t_solve'],
                'Yes' if res['converged'] else 'No',
                res['B_center'][2]
            ))
    print('')
