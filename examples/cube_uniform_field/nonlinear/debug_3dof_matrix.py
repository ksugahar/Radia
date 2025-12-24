"""
Debug script to compare matrix elements between LU and HACApK for 3DOF tetrahedra.

This script creates a simple tetrahedral mesh and compares the matrix elements
used by LU solver vs HACApK.
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../src/radia'))

import radia as rad
import math

# Simple tetrahedra test
def create_simple_tetra():
    """Create a single tetrahedron for matrix element debugging."""
    rad.UtiDelAll()
    rad.FldUnits('m')

    # Single tetrahedron
    vertices = [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.5, 1.0, 0.0],
        [0.5, 0.5, 1.0]
    ]
    faces = [[1, 3, 2], [1, 2, 4], [2, 3, 4], [3, 1, 4]]  # 1-indexed

    elem = rad.ObjPolyhdr(vertices, faces, [0.0, 0.0, 0.0])
    return elem


def create_two_tetra():
    """Create two tetrahedra for interaction matrix debugging."""
    rad.UtiDelAll()
    rad.FldUnits('m')

    # Tetrahedron 1 at origin
    v1 = [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.5, 1.0, 0.0],
        [0.5, 0.5, 1.0]
    ]

    # Tetrahedron 2 offset in x direction
    v2 = [
        [2.0, 0.0, 0.0],
        [3.0, 0.0, 0.0],
        [2.5, 1.0, 0.0],
        [2.5, 0.5, 1.0]
    ]

    faces = [[1, 3, 2], [1, 2, 4], [2, 3, 4], [3, 1, 4]]  # 1-indexed

    elem1 = rad.ObjPolyhdr(v1, faces, [0.0, 0.0, 0.0])
    elem2 = rad.ObjPolyhdr(v2, faces, [0.0, 0.0, 0.0])

    grp = rad.ObjCnt([elem1, elem2])
    return grp, 2


def main():
    print("=" * 60)
    print("3DOF Tetrahedra Matrix Element Debug")
    print("=" * 60)

    # Create two tetrahedra
    container, n_elem = create_two_tetra()
    print(f"Created {n_elem} tetrahedra ({n_elem * 3} DOF)")

    # Apply linear material (no BH curve complications)
    mat = rad.MatLin(99)  # mu_r = 100
    rad.MatApl(container, mat)

    # Apply external field
    H_EXT = 10000.0  # 10 kA/m
    MU_0 = 4 * math.pi * 1e-7
    B_ext = MU_0 * H_EXT
    ext = rad.ObjBckg([0, 0, B_ext])
    grp = rad.ObjCnt([container, ext])

    # Test LU solver
    print("\n--- LU Solver (Method 0) ---")
    result_lu = rad.Solve(grp, 0.0001, 100, 0)
    M_lu = rad.ObjM(container)
    print(f"LU iterations: {result_lu[3]:.0f}")
    print(f"LU residual: {result_lu[0]:.6e}")
    for i, m in enumerate(M_lu):
        print(f"  Element {i}: M = ({m[1][0]:.1f}, {m[1][1]:.1f}, {m[1][2]:.1f}) A/m")

    # Reset and test HACApK solver
    rad.UtiDelAll()
    container, n_elem = create_two_tetra()
    mat = rad.MatLin(99)
    rad.MatApl(container, mat)
    ext = rad.ObjBckg([0, 0, B_ext])
    grp = rad.ObjCnt([container, ext])

    print("\n--- HACApK Solver (Method 2) ---")
    result_hacapk = rad.Solve(grp, 0.0001, 100, 2)
    M_hacapk = rad.ObjM(container)
    print(f"HACApK iterations: {result_hacapk[3]:.0f}")
    print(f"HACApK residual: {result_hacapk[0]:.6e}")
    for i, m in enumerate(M_hacapk):
        print(f"  Element {i}: M = ({m[1][0]:.1f}, {m[1][1]:.1f}, {m[1][2]:.1f}) A/m")

    # Compare
    print("\n--- Comparison ---")
    for i in range(len(M_lu)):
        lu_M = M_lu[i][1]
        hacapk_M = M_hacapk[i][1]
        diff = [hacapk_M[j] - lu_M[j] for j in range(3)]
        print(f"Element {i} difference: ({diff[0]:.1f}, {diff[1]:.1f}, {diff[2]:.1f}) A/m")


if __name__ == '__main__':
    main()
