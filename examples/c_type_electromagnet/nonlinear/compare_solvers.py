#!/usr/bin/env python
"""
Compare all three solvers for nonlinear electromagnet.

Solver Methods:
- Method 0: LU (Dense direct solver)
- Method 1: BiCGSTAB (Iterative solver)
- Method 2: HACApK (H-matrix accelerated)
"""

import sys
import os
import time
import numpy as np

# Add Radia path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../src'))
import radia as rad
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../src/radia'))
from meg_to_vol import parse_meg_file


def load_bh_curve(filepath):
    """Load B-H curve from text file."""
    data = np.loadtxt(filepath, comments='#')
    return data.tolist()


def create_model():
    """Create the electromagnet model."""
    rad.UtiDelAll()
    rad.FldUnits('m')

    # File paths
    script_dir = os.path.dirname(os.path.abspath(__file__))
    meg_file = r"S:/ELF_MAGIC/2020_03_07_CEFC_2020/model_C-Type/nonlinear_20000AT/ELF_MMB8T_EIEM2_1x1x1/ELF_magic.meg"
    bh_file = os.path.join(script_dir, "BH.txt")

    # Load B-H curve
    bh_data = load_bh_curve(bh_file)

    # Parse MEG file
    nodes, elements, scale = parse_meg_file(meg_file)

    # Extract yoke elements
    yoke_elements = [e for e in elements if e['type'] == 'MMB8T']

    # Create nonlinear material
    mat = rad.MatSatIsoTab(bh_data)

    # Create yoke geometry
    yoke_objs = []
    for elem in yoke_elements:
        vertices = [list(nodes[nid]) for nid in elem['nodes']]
        hex_obj = rad.ObjHexahedron(vertices, [0, 0, 0])
        rad.MatApl(hex_obj, mat)
        yoke_objs.append(hex_obj)

    yoke = rad.ObjCnt(yoke_objs)

    # Create coil
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

    return model, len(yoke_elements)


def test_solver(method, name, prec=0.001, max_iter=100):
    """Test a specific solver method."""
    model, n_elem = create_model()

    # Set HACApK parameters if needed
    if method == 2:
        rad.SetHACApKParams(1e-4, 10, 2.0)

    t_start = time.time()
    try:
        result = rad.Solve(model, prec, max_iter, method, image='+x-z')
        t_solve = time.time() - t_start
        converged = True
        max_M = result[0]
        rel_err = result[1]
    except Exception as e:
        t_solve = time.time() - t_start
        converged = False
        max_M = 0
        rel_err = 0

    # Compute field at gap center
    B = rad.Fld(model, 'b', [0, 0, 0])
    B_mag = np.sqrt(B[0]**2 + B[1]**2 + B[2]**2)

    return {
        'name': name,
        'method': method,
        'converged': converged,
        't_solve': t_solve,
        'max_M': max_M,
        'rel_err': rel_err,
        'B_mag': B_mag,
        'Bx': B[0],
        'By': B[1],
        'Bz': B[2]
    }


def main():
    print("=" * 70)
    print("Nonlinear Electromagnet - Solver Comparison")
    print("=" * 70)

    # Test parameters
    prec = 0.001
    max_iter = 100

    print(f"\nTest parameters:")
    print(f"  Precision: {prec}")
    print(f"  Max iterations: {max_iter}")
    print(f"  Image symmetry: '+x-z'")

    results = []

    # Test LU solver
    print("\n" + "-" * 50)
    print("Testing LU Solver (Method 0)...")
    print("-" * 50)
    r = test_solver(0, 'LU', prec, max_iter)
    results.append(r)
    print(f"  Converged: {r['converged']}")
    print(f"  Time: {r['t_solve']:.3f} s")
    print(f"  |B| at gap: {r['B_mag']*1000:.2f} mT")

    # Test BiCGSTAB solver
    print("\n" + "-" * 50)
    print("Testing BiCGSTAB Solver (Method 1)...")
    print("-" * 50)
    r = test_solver(1, 'BiCGSTAB', prec, max_iter * 5)  # More iterations for BiCGSTAB
    results.append(r)
    print(f"  Converged: {r['converged']}")
    print(f"  Time: {r['t_solve']:.3f} s")
    print(f"  |B| at gap: {r['B_mag']*1000:.2f} mT")

    # Test HACApK solver
    print("\n" + "-" * 50)
    print("Testing HACApK Solver (Method 2)...")
    print("-" * 50)
    r = test_solver(2, 'HACApK', prec, max_iter * 5)  # More iterations for HACApK
    results.append(r)
    print(f"  Converged: {r['converged']}")
    print(f"  Time: {r['t_solve']:.3f} s")
    print(f"  |B| at gap: {r['B_mag']*1000:.2f} mT")

    # Summary table
    print("\n" + "=" * 70)
    print("Summary")
    print("=" * 70)
    print(f"\n{'Solver':<12} {'Method':<8} {'Converged':<10} {'Time (s)':<10} {'|B| (mT)':<10}")
    print("-" * 60)
    for r in results:
        conv_str = 'Yes' if r['converged'] else 'No'
        print(f"{r['name']:<12} {r['method']:<8} {conv_str:<10} {r['t_solve']:<10.3f} {r['B_mag']*1000:<10.2f}")
    print("-" * 60)

    print("\nField components at gap center (mT):")
    print(f"{'Solver':<12} {'Bx':<10} {'By':<10} {'Bz':<10}")
    print("-" * 42)
    for r in results:
        print(f"{r['name']:<12} {r['Bx']*1000:<10.2f} {r['By']*1000:<10.2f} {r['Bz']*1000:<10.2f}")
    print("-" * 42)

    print("\n" + "=" * 70)
    print("Solver Comparison Complete")
    print("=" * 70)

    return results


if __name__ == "__main__":
    main()
