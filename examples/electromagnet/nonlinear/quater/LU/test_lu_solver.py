#!/usr/bin/env python
"""
Test LU solver (Method 0) for nonlinear electromagnet.

LU solver: Dense direct solver using LAPACK dgesv.
- Best for small problems (N < 500 elements)
- Guaranteed convergence for well-posed problems
- Memory: O(N^2)
"""

import sys
import os
import time
import numpy as np

# Add Radia path
script_dir = os.path.dirname(os.path.abspath(__file__))
repo_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(script_dir)))))
sys.path.insert(0, os.path.join(repo_root, 'src', 'radia'))
sys.path.insert(0, os.path.join(repo_root, 'src'))
import radia as rad
from meg_to_vol import parse_meg_file


def load_bh_curve(filepath):
    """Load B-H curve from text file."""
    data = np.loadtxt(filepath, comments='#')
    return data.tolist()


def main():
    print("=" * 70)
    print("Nonlinear Electromagnet - LU Solver (Method 0)")
    print("=" * 70)

    # Initialize Radia
    rad.UtiDelAll()
    rad.FldUnits('m')

    # File paths
    script_dir = os.path.dirname(os.path.abspath(__file__))
    nonlinear_dir = os.path.dirname(os.path.dirname(script_dir))  # nonlinear folder
    # Use quarter model (13 elements instead of 52)
    meg_file = r"S:/ELF_MAGIC/2020_03_07_CEFC_2020/model_C-Type/nonlinear_20000AT/ELF_MMB8T_EIEM2_1x1x1/quater/ELF_magic.meg"
    bh_file = os.path.join(nonlinear_dir, "BH.txt")

    # Load B-H curve
    print("\nLoading B-H curve...")
    bh_data = load_bh_curve(bh_file)
    print(f"  {len(bh_data)} data points")

    # Parse MEG file
    print("\nParsing MEG file...")
    nodes, elements, scale = parse_meg_file(meg_file)

    # Extract yoke elements
    yoke_elements = [e for e in elements if e['type'] == 'MMB8T']
    print(f"  Yoke elements: {len(yoke_elements)}")
    print(f"  DOF: {len(yoke_elements) * 6}")

    # Create nonlinear material
    mat = rad.MatSatIsoTab(bh_data)

    # Create yoke geometry
    print("\nCreating yoke geometry...")
    yoke_objs = []
    for elem in yoke_elements:
        vertices = [list(nodes[nid]) for nid in elem['nodes']]
        hex_obj = rad.ObjHexahedron(vertices, [0, 0, 0])
        rad.MatApl(hex_obj, mat)
        yoke_objs.append(hex_obj)

    yoke = rad.ObjCnt(yoke_objs)
    print(f"  Created {len(yoke_objs)} elements")

    # Create coil
    print("\nCreating coil...")
    coil_center = [0, 0.13125, 0.02625]
    r_inner, r_outer = 0.060, 0.0725
    straight_half, height = 0.0525, 0.070
    current = 10000.0
    nseg = 8

    coil = rad.ObjRaceTrk(
        coil_center,
        [r_inner, r_outer],
        [straight_half * 2, height],
        height, nseg, 'man', 'z',
        current / (height * (r_outer - r_inner))
    )

    model = rad.ObjCnt([yoke, coil])

    # Solve with LU solver
    print("\n" + "-" * 50)
    print("Solving with LU solver (Method 0)...")
    print("-" * 50)
    print("  Precision: 0.001")
    print("  Max iterations: 100")
    print("  Image: '+x-z'")

    t_start = time.time()
    try:
        result = rad.Solve(model, 0.001, 100, 0, image='+x-z')
        t_solve = time.time() - t_start
        print(f"\nResult:")
        print(f"  Max |M|: {result[0]:.2f} A/m")
        print(f"  Max |dM|/|M|: {result[1]:.6f}")
        print(f"  Solve time: {t_solve:.3f} s")
        converged = True
    except Exception as e:
        t_solve = time.time() - t_start
        print(f"\n  Solve failed: {e}")
        print(f"  Solve time: {t_solve:.3f} s")
        converged = False

    # Compute field at gap center
    print("\nField at gap center (0, 0, 0):")
    B = rad.Fld(model, 'b', [0, 0, 0])
    B_mag = np.sqrt(B[0]**2 + B[1]**2 + B[2]**2)
    print(f"  Bx = {B[0]*1000:.2f} mT")
    print(f"  By = {B[1]*1000:.2f} mT")
    print(f"  Bz = {B[2]*1000:.2f} mT")
    print(f"  |B| = {B_mag*1000:.2f} mT")

    print("\n" + "=" * 70)
    print("LU Solver Test Complete")
    print("=" * 70)

    return {
        'solver': 'LU',
        'method': 0,
        'converged': converged,
        't_solve': t_solve,
        'B_gap': B_mag
    }


if __name__ == "__main__":
    main()
