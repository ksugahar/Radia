"""
Test HACApK (Method 2) with 3DOF tetrahedra MMM.

This script tests the H-matrix accelerated BiCGSTAB solver (Method 2)
and compares results with LU solver (Method 0) for tetrahedral meshes.
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../src/radia'))

import radia as rad
import time
import math

# Test parameters
N_TETRA = 28  # Small mesh for quick testing
H_EXT = 200000.0  # A/m

# BH curve (soft iron)
BH_CURVE = [
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

# Tetrahedron faces (1-indexed for Radia)
TETRA_FACES = [[1, 3, 2], [1, 2, 4], [2, 3, 4], [3, 1, 4]]


def create_tetrahedral_mesh():
    """Create simple tetrahedral mesh from a cube."""
    rad.UtiDelAll()
    rad.FldUnits('m')

    size = 0.1  # 10cm cube
    half = size / 2.0
    n = 2  # 2x2x2 division

    dx = size / n

    elements = []

    # Create tetrahedra from hexahedral sub-cubes
    # Each cube is divided into 6 tetrahedra
    for iz in range(n):
        for iy in range(n):
            for ix in range(n):
                x0 = -half + ix * dx
                y0 = -half + iy * dx
                z0 = -half + iz * dx
                x1 = x0 + dx
                y1 = y0 + dx
                z1 = z0 + dx

                # 8 vertices of cube
                v = [
                    [x0, y0, z0],  # 0
                    [x1, y0, z0],  # 1
                    [x1, y1, z0],  # 2
                    [x0, y1, z0],  # 3
                    [x0, y0, z1],  # 4
                    [x1, y0, z1],  # 5
                    [x1, y1, z1],  # 6
                    [x0, y1, z1],  # 7
                ]

                # Divide cube into 6 tetrahedra (standard BCC decomposition)
                tetra_indices = [
                    [0, 1, 3, 4],  # tetra 1
                    [1, 2, 3, 6],  # tetra 2
                    [1, 5, 4, 6],  # tetra 3
                    [3, 4, 6, 7],  # tetra 4
                    [1, 4, 3, 6],  # tetra 5 (center tetra connecting all)
                    [4, 5, 1, 6],  # tetra 6 (redundant with 3, need different decomp)
                ]

                # Actually use a simpler 5-tetra decomposition
                tetra_indices = [
                    [0, 1, 2, 5],
                    [0, 2, 3, 7],
                    [0, 5, 7, 4],
                    [2, 5, 6, 7],
                    [0, 2, 5, 7],  # center tetra
                ]

                for ti in tetra_indices:
                    verts = [v[ti[0]], v[ti[1]], v[ti[2]], v[ti[3]]]
                    elem = rad.ObjPolyhdr(verts, TETRA_FACES, [0.0, 0.0, 0.0])
                    elements.append(elem)

    grp = rad.ObjCnt(elements)
    return grp, len(elements)


def run_test(solver_method, solver_name):
    """Run test with specified solver method."""
    print(f"\n--- Testing {solver_name} solver ---")

    # Create mesh
    container, n_elem = create_tetrahedral_mesh()
    print(f"Created {n_elem} tetrahedra ({n_elem * 3} DOF)")

    # Apply nonlinear material
    mat = rad.MatSatIsoTab(BH_CURVE)
    rad.MatApl(container, mat)

    # Apply external field
    MU_0 = 4 * math.pi * 1e-7
    B_ext = MU_0 * H_EXT
    ext = rad.ObjBckg([0, 0, B_ext])
    grp = rad.ObjCnt([container, ext])

    # Solve
    t0 = time.time()
    try:
        res = rad.Solve(grp, 0.0001, 100, solver_method)
        t_solve = time.time() - t0

        # Get average magnetization
        M_all = rad.ObjM(container)
        M_z_sum = sum(m[1][2] for m in M_all)
        M_avg_z = M_z_sum / len(M_all)

        print(f"{solver_name} time: {t_solve:.3f} s")
        print(f"{solver_name} M_avg_z: {M_avg_z:.0f} A/m")
        print(f"{solver_name} residual: {res[0]:.6e}")
        print(f"{solver_name} iterations: {res[3]:.0f}")

        return {
            'n_elem': n_elem,
            'ndof': n_elem * 3,
            't_solve': t_solve,
            'M_avg_z': M_avg_z,
            'residual': res[0],
            'iterations': int(res[3]),
        }

    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()
        return None


def main():
    print("=" * 60)
    print("3DOF Tetrahedra HACApK Test")
    print("=" * 60)
    print(f"H_ext = {H_EXT} A/m")

    # Test LU solver (Method 0) as reference
    result_lu = run_test(0, 'LU')

    # Test HACApK solver (Method 2)
    result_hacapk = run_test(2, 'HACApK')

    # Summary
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)

    if result_lu and result_hacapk:
        diff = abs(result_hacapk['M_avg_z'] - result_lu['M_avg_z'])
        diff_pct = diff / result_lu['M_avg_z'] * 100
        print(f"LU M_avg_z:     {result_lu['M_avg_z']:.0f} A/m")
        print(f"HACApK M_avg_z: {result_hacapk['M_avg_z']:.0f} A/m")
        print(f"Difference:     {diff:.0f} A/m ({diff_pct:.2f}%)")

        if diff_pct < 1.0:
            print("\nSUCCESS: Results match within 1%")
        else:
            print("\nFAILURE: Results differ by more than 1%")


if __name__ == '__main__':
    main()
