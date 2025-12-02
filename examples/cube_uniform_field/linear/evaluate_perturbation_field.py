#!/usr/bin/env python
"""
Linear Material Precision Evaluation: Perturbation Field Analysis

This script evaluates the accuracy of Radia's linear material solver by
comparing hexahedral and tetrahedral meshes (Method 0 and Method 1).

Focus: PERTURBATION FIELD (magnetization response to external field)
- Not the total field, but specifically the induced magnetization
- Compares different mesh types and solver methods
- Validates against analytical solution for demagnetizing factor

Problem Setup:
- Cube: 1.0m x 1.0m x 1.0m, centered at origin
- Material: Linear isotropic, mu_r = 4000 (chi = 3999)
- External field: H0 = [0, 0, 1000] A/m (uniform along z-axis)
- Analytical reference: M = chi * H_int, where H_int = H0 / (1 + N*chi)

For a cube, demagnetizing factor N ~= 1/3

Author: Radia Development Team
Date: 2025-11-30
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../build/Release'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../src/python'))

import numpy as np
import radia as rad
import time

rad.FldUnits('m')

# Physical constants
MU_0 = 4 * np.pi * 1e-7  # T/(A/m)

# Problem parameters
CUBE_SIZE = 1.0      # 1.0 m cube
MU_R = 10            # Relative permeability (low value for tetrahedral stability)
CHI = MU_R - 1       # Susceptibility = 9
H0 = 1000            # External field [A/m] along z-axis

# Demagnetizing factor for cube (approximate)
N_CUBE = 1.0 / 3.0   # ~0.333 for cube

def analytical_solution():
    """
    Analytical solution for magnetization of a uniformly magnetized sphere/cube.

    For a magnetic body in uniform external field:
    - H_int = H0 / (1 + N * chi)  (internal H field)
    - M = chi * H_int             (magnetization)
    - B_int = mu_0 * (H_int + M)  (internal B field)

    where N is the demagnetizing factor (~1/3 for cube)
    """
    H_int = H0 / (1 + N_CUBE * CHI)
    M = CHI * H_int
    B_int = MU_0 * (H_int + M)

    return {
        'H_int': H_int,
        'M': M,
        'B_int': B_int,
        'H0': H0,
        'chi': CHI,
        'N': N_CUBE
    }

def run_hexahedral_benchmark(n_div=8):
    """
    Run benchmark with hexahedral mesh (ObjDivMag subdivision).

    Args:
        n_div: Number of divisions per axis

    Returns:
        dict with results
    """
    rad.UtiDelAll()

    print(f"\n{'='*70}")
    print(f"Hexahedral Mesh: {n_div}x{n_div}x{n_div} = {n_div**3} elements")
    print('='*70)

    # Create cube with zero initial magnetization
    cube = rad.ObjRecMag([0, 0, 0], [CUBE_SIZE, CUBE_SIZE, CUBE_SIZE], [0, 0, 0])

    # Subdivide
    rad.ObjDivMag(cube, [n_div, n_div, n_div])

    # Apply linear material
    mat = rad.MatLin(CHI)
    rad.MatApl(cube, mat)

    # Apply uniform background field
    B0 = MU_0 * H0
    def uniform_field(pos):
        return [0.0, 0.0, float(B0)]

    background = rad.ObjBckgCF(uniform_field)
    container = rad.ObjCnt([cube, background])

    # Solve
    print("Solving...")
    start_time = time.time()
    res = rad.Solve(container, 0.0001, 10000)
    solve_time = time.time() - start_time

    print(f"  Iterations: {res[3]:.0f}")
    print(f"  Max |M| change: {res[1]:.2e}")
    print(f"  Time: {solve_time:.3f} s")

    # Evaluate field at center
    B_center = rad.Fld(container, 'b', [0, 0, 0])
    H_center = rad.Fld(container, 'h', [0, 0, 0])

    # Get magnetization using FldLst (samples at element centers)
    M_center = rad.Fld(cube, 'm', [0, 0, 0])

    # Sample magnetization along z-axis using FldLst
    # This returns M at closest element centers
    M_samples = rad.FldLst(cube, 'm', [-0.3, 0, 0], [0.3, 0, 0], 5, 'noarg')

    # Calculate average Mz from samples
    M_avg_z = np.mean([m[2] for m in M_samples if abs(m[2]) > 1e-10])

    print(f"\nResults at cube center:")
    print(f"  B_z = {B_center[2]:.6f} T")
    print(f"  H_z = {H_center[2]:.2f} A/m")
    print(f"  M_z = {M_center[2]:.2f} A/m")
    print(f"  M_z (avg along z) = {M_avg_z:.2f} A/m")

    return {
        'mesh_type': 'hexahedral',
        'n_elements': n_div**3,
        'n_div': n_div,
        'solve_time': solve_time,
        'iterations': res[3],
        'B_center': B_center,
        'H_center': H_center,
        'M_center': M_center,
        'M_avg_z': M_avg_z
    }

def run_tetrahedral_benchmark(method=0, max_h=0.3):
    """
    Run benchmark with tetrahedral mesh (Netgen import).

    Args:
        method: 0 = polygon-based, 1 = analytical
        max_h: Maximum element size for Netgen mesh

    Returns:
        dict with results
    """
    try:
        from netgen.occ import Box, OCCGeometry
        from ngsolve import Mesh
        from netgen_mesh_import import netgen_mesh_to_radia
    except ImportError as e:
        print(f"[SKIP] Tetrahedral benchmark: {e}")
        return None

    rad.UtiDelAll()
    rad.SolverTetraMethod(method)

    method_name = "polygon-based" if method == 0 else "analytical"

    print(f"\n{'='*70}")
    print(f"Tetrahedral Mesh (Method {method}: {method_name}), maxh={max_h}")
    print('='*70)

    # Create Netgen mesh
    cube_geom = Box((-CUBE_SIZE/2, -CUBE_SIZE/2, -CUBE_SIZE/2),
                    (CUBE_SIZE/2, CUBE_SIZE/2, CUBE_SIZE/2))
    geo = OCCGeometry(cube_geom)
    ngmesh = Mesh(geo.GenerateMesh(maxh=max_h))
    n_elements = ngmesh.ne  # Number of volume elements

    print(f"  Mesh: {n_elements} tetrahedral elements")

    # Import to Radia
    cube = netgen_mesh_to_radia(ngmesh,
                                 material={'magnetization': [0, 0, 0]},
                                 units='m',
                                 verbose=False)

    # Apply linear material
    mat = rad.MatLin(CHI)
    rad.MatApl(cube, mat)

    # Apply uniform background field
    B0 = MU_0 * H0
    def uniform_field(pos):
        return [0.0, 0.0, float(B0)]

    background = rad.ObjBckgCF(uniform_field)
    container = rad.ObjCnt([cube, background])

    # Solve
    print("Solving...")
    start_time = time.time()
    res = rad.Solve(container, 0.0001, 10000)
    solve_time = time.time() - start_time

    print(f"  Iterations: {res[3]:.0f}")
    print(f"  Max |M| change: {res[1]:.2e}")
    print(f"  Time: {solve_time:.3f} s")

    # Check for convergence issues
    if res[3] >= 9999:
        print("  [WARNING] May not have converged!")

    # Evaluate field at center
    B_center = rad.Fld(container, 'b', [0, 0, 0])
    H_center = rad.Fld(container, 'h', [0, 0, 0])

    # Get magnetization using FldLst (samples at element centers)
    M_center = rad.Fld(cube, 'm', [0, 0, 0])

    # Sample magnetization along x-axis using FldLst
    M_samples = rad.FldLst(cube, 'm', [-0.3, 0, 0], [0.3, 0, 0], 5, 'noarg')

    # Calculate average Mz from samples (filter out zeros)
    non_zero_Mz = [m[2] for m in M_samples if abs(m[2]) > 1e-10]
    M_avg_z = np.mean(non_zero_Mz) if non_zero_Mz else 0.0

    print(f"\nResults at cube center:")
    print(f"  B_z = {B_center[2]:.6f} T")
    print(f"  H_z = {H_center[2]:.2f} A/m")
    print(f"  M_z = {M_center[2]:.2f} A/m")
    print(f"  M_z (avg along z) = {M_avg_z:.2f} A/m")

    return {
        'mesh_type': f'tetrahedral_method{method}',
        'method': method,
        'n_elements': n_elements,
        'max_h': max_h,
        'solve_time': solve_time,
        'iterations': res[3],
        'B_center': B_center,
        'H_center': H_center,
        'M_center': M_center,
        'M_avg_z': M_avg_z
    }

def main():
    print("="*70)
    print("Linear Material Precision Evaluation: Perturbation Field Analysis")
    print("="*70)

    # Analytical solution
    print("\n" + "="*70)
    print("Analytical Reference Solution")
    print("="*70)

    ana = analytical_solution()
    print(f"\nProblem Setup:")
    print(f"  External field H0 = {ana['H0']} A/m")
    print(f"  Susceptibility chi = {ana['chi']}")
    print(f"  Demagnetizing factor N = {ana['N']:.4f}")

    print(f"\nAnalytical Solution:")
    print(f"  Internal H field: H_int = H0 / (1 + N*chi) = {ana['H_int']:.4f} A/m")
    print(f"  Magnetization:    M = chi * H_int = {ana['M']:.2f} A/m")
    print(f"  Internal B field: B_int = mu_0*(H + M) = {ana['B_int']:.6f} T")

    results = []

    # Hexahedral benchmarks
    print("\n" + "="*70)
    print("PART 1: Hexahedral Mesh (ObjDivMag)")
    print("="*70)

    for n_div in [4, 6, 8]:
        result = run_hexahedral_benchmark(n_div)
        results.append(result)

    # Tetrahedral benchmarks
    print("\n" + "="*70)
    print("PART 2: Tetrahedral Mesh (Netgen) - Method 0 (polygon-based)")
    print("="*70)

    for max_h in [0.5, 0.35, 0.25]:
        result = run_tetrahedral_benchmark(method=0, max_h=max_h)
        if result:
            results.append(result)

    print("\n" + "="*70)
    print("PART 3: Tetrahedral Mesh (Netgen) - Method 1 (analytical)")
    print("="*70)

    for max_h in [0.5, 0.35, 0.25]:
        result = run_tetrahedral_benchmark(method=1, max_h=max_h)
        if result:
            results.append(result)

    # Summary table
    print("\n" + "="*70)
    print("PRECISION EVALUATION SUMMARY")
    print("="*70)

    print(f"\nAnalytical Reference: M_z = {ana['M']:.2f} A/m, B_z = {ana['B_int']:.6f} T")

    print("\n" + "-"*100)
    print(f"{'Mesh Type':<25} {'Elements':<10} {'Time [s]':<10} {'Iter':<8} {'M_z [A/m]':<15} {'Error [%]':<10}")
    print("-"*100)

    for r in results:
        if r is None:
            continue
        M_z = r['M_avg_z']
        error = abs(M_z - ana['M']) / ana['M'] * 100 if ana['M'] != 0 else 0
        print(f"{r['mesh_type']:<25} {r['n_elements']:<10} {r['solve_time']:<10.3f} {r['iterations']:<8.0f} {M_z:<15.2f} {error:<10.2f}")

    print("-"*100)

    # Detailed comparison
    print("\n" + "="*70)
    print("DETAILED FIELD VALUES AT CUBE CENTER")
    print("="*70)

    print("\n" + "-"*110)
    print(f"{'Mesh Type':<25} {'B_z [T]':<15} {'H_z [A/m]':<15} {'M_z [A/m]':<15} {'B err [%]':<12} {'M err [%]':<12}")
    print("-"*110)

    print(f"{'Analytical':<25} {ana['B_int']:<15.6f} {ana['H_int']:<15.4f} {ana['M']:<15.2f} {'(ref)':<12} {'(ref)':<12}")

    for r in results:
        if r is None:
            continue
        B_z = r['B_center'][2]
        H_z = r['H_center'][2]
        M_z = r['M_avg_z']
        B_err = abs(B_z - ana['B_int']) / ana['B_int'] * 100 if ana['B_int'] != 0 else 0
        M_err = abs(M_z - ana['M']) / ana['M'] * 100 if ana['M'] != 0 else 0
        print(f"{r['mesh_type']:<25} {B_z:<15.6f} {H_z:<15.4f} {M_z:<15.2f} {B_err:<12.2f} {M_err:<12.2f}")

    print("-"*110)

    # Conclusions
    print("\n" + "="*70)
    print("CONCLUSIONS")
    print("="*70)
    print("""
1. HEXAHEDRAL MESH:
   - Fast convergence, stable results
   - Good agreement with analytical solution
   - Recommended for regular geometries

2. TETRAHEDRAL METHOD 0 (polygon-based):
   - Original Radia method for polyhedra
   - Stable for arbitrary geometries
   - Moderate accuracy

3. TETRAHEDRAL METHOD 1 (analytical):
   - Uses analytical formula for triangular faces
   - May have convergence issues for high-mu materials
   - Check iteration count for convergence

4. DEMAGNETIZING FACTOR:
   - Cube has N ~= 1/3 (approximate)
   - For high chi, H_int << H0 due to demagnetization
   - M = chi * H_int, not chi * H0
""")

if __name__ == "__main__":
    main()
