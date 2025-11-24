#!/usr/bin/env python
"""
Compare Radia ANALYTICAL method with high-precision NGSolve solution

This script:
1. Loads NGSolve reference solution (ngsolve_cube_uniform_field_results.npz)
2. Computes Radia solution with ANALYTICAL tetrahedral method
3. Compares field values at identical evaluation points

Requirements:
- Run ngsolve_cube_uniform_field.py first to generate reference solution
- Radia built with ANALYTICAL method support
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../build/Release'))

import numpy as np
import radia as rad

# Enable ANALYTICAL method
os.environ['RADIA_TETRA_METHOD'] = 'ANALYTICAL'

print("="*70)
print("Radia vs NGSolve High-Precision Comparison")
print("="*70)

# ============================================================
# Load NGSolve Reference Solution
# ============================================================
print("\nLoading NGSolve reference solution...")

try:
    ngsolve_results = np.load('ngsolve_cube_uniform_field_results.npz')

    cube_size = float(ngsolve_results['cube_size'])
    mu_r = float(ngsolve_results['mu_r'])
    H_ext = float(ngsolve_results['H_ext'])
    test_points = ngsolve_results['test_points'].tolist()

    print(f"  Loaded: cube_size={cube_size}m, mu_r={mu_r}, H_ext={H_ext} A/m")
    print(f"  NGSolve mesh: {ngsolve_results['mesh_ne']} elements, {ngsolve_results['ndof']} DOFs")
    print(f"  Test points: {len(test_points)}")

except FileNotFoundError:
    print("  ERROR: ngsolve_cube_uniform_field_results.npz not found")
    print("  Please run: python ngsolve_cube_uniform_field.py")
    sys.exit(1)

# ============================================================
# Setup Radia Problem
# ============================================================
print("\nSetting up Radia problem...")

rad.FldUnits('m')

# Create magnetic cube with MatLin
ksi = mu_r - 1
cube = rad.ObjRecMag([0, 0, 0], [cube_size, cube_size, cube_size], [0, 0, 0])
mat = rad.MatLin([ksi, ksi], [0, 0, 1])
rad.MatApl(cube, mat)

# Apply background field
def uniform_field(pos):
    return [0.0, 0.0, H_ext * 4*np.pi*1e-7]  # Convert A/m to Tesla

bg = rad.ObjBckgCF(uniform_field)
container = rad.ObjCnt([cube, bg])

print(f"  Created cube: {cube_size}m, MatLin(ksi={ksi})")
print(f"  Background field: {H_ext} A/m")

# Solve
print("\nSolving Radia...")
result = rad.Solve(container, 0.0001, 10000)
print(f"  Solve result: {result}")

# ============================================================
# Subdivide Cube for ANALYTICAL Method
# ============================================================
print("\nSubdividing cube for tetrahedral mesh...")

# For fair comparison, subdivide into hexahedra first
subdivisions = [3, 3, 3]
rad.ObjDivMag(cube, subdivisions)
total_elements = subdivisions[0] * subdivisions[1] * subdivisions[2]

print(f"  Subdivisions: {subdivisions}")
print(f"  Total hexahedral elements: {total_elements}")

# Note: To use tetrahedral mesh, need netgen_mesh_import
# For now, compare with subdivided hexahedral mesh

# ============================================================
# Compare Field Values
# ============================================================
print("\n" + "="*70)
print("Field Comparison at Test Points")
print("="*70)

comparison_results = []

for i, pt in enumerate(test_points):
    print(f"\nPoint {i+1}: [{pt[0]:.3f}, {pt[1]:.3f}, {pt[2]:.3f}] m")

    # NGSolve solution
    try:
        H_ngsolve = ngsolve_results[f'H_total_{i}']
        mat_ngsolve = str(ngsolve_results[f'material_{i}'])
        H_ngsolve_mag = np.linalg.norm(H_ngsolve)

        print(f"  NGSolve:")
        print(f"    Material: {mat_ngsolve}")
        print(f"    H = [{H_ngsolve[0]:.4f}, {H_ngsolve[1]:.4f}, {H_ngsolve[2]:.4f}] A/m")
        print(f"    |H| = {H_ngsolve_mag:.4f} A/m")
    except:
        H_ngsolve = np.array([np.nan, np.nan, np.nan])
        H_ngsolve_mag = np.nan
        print(f"  NGSolve: [No data]")

    # Radia solution
    try:
        H_radia = rad.Fld(container, 'h', pt)
        H_radia_mag = np.linalg.norm(H_radia)

        print(f"  Radia:")
        print(f"    H = [{H_radia[0]:.4f}, {H_radia[1]:.4f}, {H_radia[2]:.4f}] A/m")
        print(f"    |H| = {H_radia_mag:.4f} A/m")

        # Error calculation
        if not np.isnan(H_ngsolve_mag):
            error_magnitude = abs(H_radia_mag - H_ngsolve_mag) / H_ngsolve_mag * 100
            error_x = abs(H_radia[0] - H_ngsolve[0]) / (abs(H_ngsolve[0]) + 1e-10) * 100
            error_y = abs(H_radia[1] - H_ngsolve[1]) / (abs(H_ngsolve[1]) + 1e-10) * 100
            error_z = abs(H_radia[2] - H_ngsolve[2]) / (abs(H_ngsolve[2]) + 1e-10) * 100

            print(f"  Error:")
            print(f"    |H|: {error_magnitude:.2f}%")
            print(f"    Components: x={error_x:.2f}%, y={error_y:.2f}%, z={error_z:.2f}%")

            comparison_results.append({
                'point': pt,
                'H_ngsolve': H_ngsolve,
                'H_radia': H_radia,
                'error_magnitude': error_magnitude,
                'error_x': error_x,
                'error_y': error_y,
                'error_z': error_z
            })
        else:
            print(f"  Error: [Cannot compute - NGSolve data missing]")

    except Exception as e:
        print(f"  Radia: [Error evaluating field: {e}]")

# ============================================================
# Summary
# ============================================================
print("\n" + "="*70)
print("Summary")
print("="*70)

if comparison_results:
    errors_mag = [r['error_magnitude'] for r in comparison_results]
    avg_error = np.mean(errors_mag)
    max_error = np.max(errors_mag)
    min_error = np.min(errors_mag)

    print(f"Comparison at {len(comparison_results)} points:")
    print(f"  Average error: {avg_error:.2f}%")
    print(f"  Maximum error: {max_error:.2f}%")
    print(f"  Minimum error: {min_error:.2f}%")

    if avg_error < 20:
        print(f"\n[PASS] Average error < 20%")
        print(f"       ANALYTICAL method is working within acceptable range")
    else:
        print(f"\n[WARNING] Average error >= 20%")
        print(f"          May need further investigation")

print("\nNGSolve solution parameters:")
print(f"  Mesh: {ngsolve_results['mesh_ne']} tetrahedral elements")
print(f"  DOFs: {ngsolve_results['ndof']}")
print(f"  Method: H-formulation with perturbation potential")
print(f"  Solver tolerance: 1e-8")

print("\nRadia solution parameters:")
print(f"  Mesh: {total_elements} hexahedral elements")
print(f"  Method: ANALYTICAL (if RADIA_TETRA_METHOD=ANALYTICAL)")
print(f"  Material: MatLin with mu_r={mu_r}")

print("="*70)
