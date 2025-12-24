#!/usr/bin/env python
"""
Test single matrix element comparison between LU and HACApK.

This tests if the matrix elements computed by both solvers are identical.
"""

import sys
import os
import numpy as np

# Add Radia to path
_src_path = os.path.join(os.path.dirname(__file__), '../../../src/radia')
sys.path.insert(0, _src_path)
import radia as rad

MU_0 = 4 * np.pi * 1e-7
H_EXT = 200000.0

# B-H curve data
BH_DATA = [
    [0.0, 0.0], [100.0, 0.1], [200.0, 0.3], [500.0, 0.8], [1000.0, 1.2],
    [2000.0, 1.5], [5000.0, 1.7], [10000.0, 1.8], [50000.0, 2.0], [100000.0, 2.1],
]

# Hexahedral face topology (1-indexed for Radia)
HEX_FACES = [
    [1, 4, 3, 2], [5, 6, 7, 8], [1, 2, 6, 5], [3, 4, 8, 7], [1, 5, 8, 4], [2, 3, 7, 6]
]


def create_single_hexahedron():
    """Create a single 1m cube centered at origin."""
    vertices = [
        [-0.5, -0.5, -0.5], [0.5, -0.5, -0.5], [0.5, 0.5, -0.5], [-0.5, 0.5, -0.5],
        [-0.5, -0.5, 0.5], [0.5, -0.5, 0.5], [0.5, 0.5, 0.5], [-0.5, 0.5, 0.5]
    ]
    return rad.ObjPolyhdr(vertices, HEX_FACES, [0, 0, 0])


def test_single_element():
    """Test single element with LU."""
    rad.FldUnits('m')
    rad.UtiDelAll()

    # Create single hexahedron
    hex_obj = create_single_hexahedron()
    container = rad.ObjCnt([hex_obj])

    # Apply linear material (fixed chi for testing)
    mu_r = 1001.0  # chi = 1000
    mat = rad.MatLin(mu_r)
    rad.MatApl(container, mat)

    # External field
    B_ext = MU_0 * H_EXT
    ext = rad.ObjBckg([0, 0, B_ext])
    grp = rad.ObjCnt([container, ext])

    # Solve with LU (single iteration for linear problem)
    result = rad.Solve(grp, 0.001, 1, 0)

    # Get magnetization
    M_data = rad.ObjM(container)
    M = M_data[0][1] if M_data else [0, 0, 0]

    print('Single element LU test:')
    print(f'  M = [{M[0]:.4f}, {M[1]:.4f}, {M[2]:.4f}]')
    print(f'  Expected M_z ~ chi * H_ext = 1000 * 200000 = 2e8 A/m (approx)')

    return M


def test_single_element_hacapk():
    """Test single element with HACApK."""
    rad.FldUnits('m')
    rad.UtiDelAll()

    # Create single hexahedron
    hex_obj = create_single_hexahedron()
    container = rad.ObjCnt([hex_obj])

    # Apply linear material (fixed chi for testing)
    mu_r = 1001.0  # chi = 1000
    mat = rad.MatLin(mu_r)
    rad.MatApl(container, mat)

    # External field
    B_ext = MU_0 * H_EXT
    ext = rad.ObjBckg([0, 0, B_ext])
    grp = rad.ObjCnt([container, ext])

    # Configure HACApK
    rad.SetHACApKParams(1e-8, 10, 2.0)

    # Solve with HACApK
    result = rad.Solve(grp, 0.001, 1, 2)

    # Get magnetization
    M_data = rad.ObjM(container)
    M = M_data[0][1] if M_data else [0, 0, 0]

    print('Single element HACApK test:')
    print(f'  M = [{M[0]:.4f}, {M[1]:.4f}, {M[2]:.4f}]')

    return M


def main():
    print('=' * 70)
    print('TEST: Single Element Matrix Comparison')
    print('=' * 70)

    M_lu = test_single_element()
    print()
    M_hacapk = test_single_element_hacapk()

    print('\n' + '=' * 70)
    print('COMPARISON')
    print('=' * 70)
    diff = np.array(M_hacapk) - np.array(M_lu)
    print(f'  LU M:      [{M_lu[0]:.4f}, {M_lu[1]:.4f}, {M_lu[2]:.4f}]')
    print(f'  HACApK M:  [{M_hacapk[0]:.4f}, {M_hacapk[1]:.4f}, {M_hacapk[2]:.4f}]')
    print(f'  Diff:      [{diff[0]:.4e}, {diff[1]:.4e}, {diff[2]:.4e}]')

    if np.linalg.norm(M_lu) > 0:
        rel_err = np.linalg.norm(diff) / np.linalg.norm(M_lu)
        print(f'  Relative error: {rel_err:.6e}')


if __name__ == '__main__':
    main()
