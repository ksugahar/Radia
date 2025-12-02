#!/usr/bin/env python
"""
Linear Material: External Field Comparison with NGSolve

This script evaluates the perturbation field OUTSIDE the magnetic cube,
which is the proper way to compare with NGSolve FEM reference solutions.

Problem Setup:
- Cube: 1.0m x 1.0m x 1.0m, centered at origin (boundaries at +/-0.5m)
- Material: Linear isotropic, mu_r = 100
- External field: H0 = [0, 0, 1] A/m (uniform along z-axis)
- Evaluation: Points OUTSIDE the cube (r > 0.5m)

Key Points:
- Radia computes perturbation field from magnetization
- Total field: H_total = H0 + H_pert
- Compare H_pert with NGSolve at points outside cube

Reference: ELF/MAGIC compare_with_ngsolve.py

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

# Problem parameters (match NGSolve reference)
CUBE_SIZE = 1.0      # 1.0 m cube
CUBE_HALF = 0.5      # half size
MU_R = 100           # Relative permeability
CHI = MU_R - 1       # Susceptibility = 99
H0 = np.array([0.0, 0.0, 1.0])  # External field [A/m] along z-axis
B0 = MU_0 * H0       # Background B field


def run_hexahedral_benchmark(n_div=8):
    """
    Run benchmark with hexahedral mesh and evaluate external field.

    Args:
        n_div: Number of divisions per axis

    Returns:
        dict with results including external field values
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
    def uniform_field(pos):
        return [float(B0[0]), float(B0[1]), float(B0[2])]

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

    # Get average magnetization
    M_samples = rad.FldLst(cube, 'm', [-0.3, 0, 0], [0.3, 0, 0], 5, 'noarg')
    M_avg = np.mean([np.array(m) for m in M_samples], axis=0)
    print(f"\nAverage magnetization: M = [{M_avg[0]:.4f}, {M_avg[1]:.4f}, {M_avg[2]:.4f}] A/m")
    print(f"|M| = {np.linalg.norm(M_avg):.4f} A/m")

    # Evaluation points OUTSIDE the cube
    mesh_step = CUBE_SIZE / n_div
    eval_points = [
        # On z-axis (above cube)
        {'point': [0.0, 0.0, CUBE_HALF + mesh_step], 'desc': 'z-axis, 1 step'},
        {'point': [0.0, 0.0, CUBE_HALF + 2*mesh_step], 'desc': 'z-axis, 2 steps'},
        {'point': [0.0, 0.0, 0.7], 'desc': 'z-axis, r=0.7m'},
        {'point': [0.0, 0.0, 0.8], 'desc': 'z-axis, r=0.8m'},
        {'point': [0.0, 0.0, 1.0], 'desc': 'z-axis, r=1.0m'},
        # On x-axis (outside cube)
        {'point': [CUBE_HALF + mesh_step, 0.0, 0.0], 'desc': 'x-axis, 1 step'},
        {'point': [CUBE_HALF + 2*mesh_step, 0.0, 0.0], 'desc': 'x-axis, 2 steps'},
        {'point': [0.7, 0.0, 0.0], 'desc': 'x-axis, r=0.7m'},
        {'point': [0.8, 0.0, 0.0], 'desc': 'x-axis, r=0.8m'},
        {'point': [1.0, 0.0, 0.0], 'desc': 'x-axis, r=1.0m'},
    ]

    print(f"\nExternal field evaluation (mesh_step = {mesh_step:.3f} m):")
    print("-" * 80)
    print(f"{'Point':<25} {'Description':<15} {'Hz_pert [A/m]':>15} {'Hz_total [A/m]':>15}")
    print("-" * 80)

    field_results = []
    for ep in eval_points:
        pt = ep['point']

        # Get field from Radia (this is total B field)
        B = rad.Fld(container, 'b', pt)
        H = rad.Fld(container, 'h', pt)

        # H from Radia is the total H including magnetization contribution
        # H_pert = H_total - H0
        H_pert = np.array(H) - H0

        # Also compute B_pert = B - B0
        B_pert = np.array(B) - B0

        pt_str = f"({pt[0]:.2f}, {pt[1]:.2f}, {pt[2]:.2f})"
        print(f"{pt_str:<25} {ep['desc']:<15} {H_pert[2]:>15.6f} {H[2]:>15.6f}")

        field_results.append({
            'point': pt,
            'description': ep['desc'],
            'H_total': np.array(H),
            'H_pert': H_pert,
            'B_total': np.array(B),
            'B_pert': B_pert
        })

    print("-" * 80)

    return {
        'mesh_type': 'hexahedral',
        'n_elements': n_div**3,
        'n_div': n_div,
        'solve_time': solve_time,
        'iterations': res[3],
        'M_avg': M_avg,
        'field_results': field_results
    }


def run_tetrahedral_benchmark(method=0, max_h=0.3):
    """
    Run benchmark with tetrahedral mesh and evaluate external field.

    Args:
        method: 0 = polygon-based, 1 = analytical
        max_h: Maximum element size for Netgen mesh

    Returns:
        dict with results including external field values, or None if failed
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
    cube_geom = Box((-CUBE_HALF, -CUBE_HALF, -CUBE_HALF),
                    (CUBE_HALF, CUBE_HALF, CUBE_HALF))
    geo = OCCGeometry(cube_geom)
    ngmesh = Mesh(geo.GenerateMesh(maxh=max_h))
    n_elements = ngmesh.ne

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
    def uniform_field(pos):
        return [float(B0[0]), float(B0[1]), float(B0[2])]

    background = rad.ObjBckgCF(uniform_field)
    container = rad.ObjCnt([cube, background])

    # Solve
    print("Solving...")
    start_time = time.time()
    res = rad.Solve(container, 0.01, 10000)  # Looser tolerance for tetrahedra
    solve_time = time.time() - start_time

    print(f"  Iterations: {res[3]:.0f}")
    print(f"  Max |M| change: {res[1]:.2e}")
    print(f"  Time: {solve_time:.3f} s")

    # Check for convergence issues
    if res[3] >= 9999 or np.isnan(res[1]) or np.isinf(res[1]):
        print("  [WARNING] Solver did not converge properly!")
        return None

    # Get average magnetization
    M_samples = rad.FldLst(cube, 'm', [-0.3, 0, 0], [0.3, 0, 0], 5, 'noarg')
    non_zero = [np.array(m) for m in M_samples if np.linalg.norm(m) > 1e-10]
    if non_zero:
        M_avg = np.mean(non_zero, axis=0)
    else:
        M_avg = np.array([0.0, 0.0, 0.0])
    print(f"\nAverage magnetization: M = [{M_avg[0]:.4f}, {M_avg[1]:.4f}, {M_avg[2]:.4f}] A/m")

    # Evaluation points OUTSIDE the cube
    eval_points = [
        {'point': [0.0, 0.0, 0.7], 'desc': 'z-axis, r=0.7m'},
        {'point': [0.0, 0.0, 0.8], 'desc': 'z-axis, r=0.8m'},
        {'point': [0.0, 0.0, 1.0], 'desc': 'z-axis, r=1.0m'},
        {'point': [0.7, 0.0, 0.0], 'desc': 'x-axis, r=0.7m'},
        {'point': [0.8, 0.0, 0.0], 'desc': 'x-axis, r=0.8m'},
        {'point': [1.0, 0.0, 0.0], 'desc': 'x-axis, r=1.0m'},
    ]

    print(f"\nExternal field evaluation:")
    print("-" * 80)
    print(f"{'Point':<25} {'Description':<15} {'Hz_pert [A/m]':>15} {'Hz_total [A/m]':>15}")
    print("-" * 80)

    field_results = []
    for ep in eval_points:
        pt = ep['point']

        # Get field from Radia
        H = rad.Fld(container, 'h', pt)
        H_pert = np.array(H) - H0

        pt_str = f"({pt[0]:.2f}, {pt[1]:.2f}, {pt[2]:.2f})"
        if np.any(np.isnan(H)):
            print(f"{pt_str:<25} {ep['desc']:<15} {'NaN':>15} {'NaN':>15}")
        else:
            print(f"{pt_str:<25} {ep['desc']:<15} {H_pert[2]:>15.6f} {H[2]:>15.6f}")

        field_results.append({
            'point': pt,
            'description': ep['desc'],
            'H_total': np.array(H),
            'H_pert': H_pert,
        })

    print("-" * 80)

    return {
        'mesh_type': f'tetrahedral_method{method}',
        'method': method,
        'n_elements': n_elements,
        'max_h': max_h,
        'solve_time': solve_time,
        'iterations': res[3],
        'M_avg': M_avg,
        'field_results': field_results
    }


def main():
    print("="*70)
    print("Linear Material: External Field Comparison")
    print("="*70)

    print(f"\nProblem Setup:")
    print(f"  Cube size:  {CUBE_SIZE} m (centered at origin)")
    print(f"  mu_r:       {MU_R}")
    print(f"  H0:         {H0} A/m")
    print(f"  B0:         {B0} T")

    results = []

    # Hexahedral benchmarks
    print("\n" + "="*70)
    print("PART 1: Hexahedral Mesh (ObjDivMag)")
    print("="*70)

    for n_div in [5, 8, 10]:
        result = run_hexahedral_benchmark(n_div)
        results.append(result)

    # Tetrahedral benchmarks
    print("\n" + "="*70)
    print("PART 2: Tetrahedral Mesh (Netgen)")
    print("="*70)

    for method in [0, 1]:
        for max_h in [0.4, 0.3]:
            result = run_tetrahedral_benchmark(method=method, max_h=max_h)
            if result:
                results.append(result)

    # Summary table
    print("\n" + "="*70)
    print("EXTERNAL FIELD COMPARISON SUMMARY")
    print("="*70)

    # Compare at fixed points
    ref_points = [
        [0.0, 0.0, 0.7],
        [0.0, 0.0, 0.8],
        [0.7, 0.0, 0.0],
        [0.8, 0.0, 0.0],
    ]

    for ref_pt in ref_points:
        print(f"\n--- Point ({ref_pt[0]:.1f}, {ref_pt[1]:.1f}, {ref_pt[2]:.1f}) ---")
        print(f"{'Mesh Type':<30} {'Elements':<10} {'Hz_pert [A/m]':>15}")
        print("-" * 60)

        for r in results:
            if r is None:
                continue
            for fr in r['field_results']:
                if np.allclose(fr['point'], ref_pt, atol=0.01):
                    Hz_pert = fr['H_pert'][2]
                    if np.isnan(Hz_pert):
                        print(f"{r['mesh_type']:<30} {r['n_elements']:<10} {'NaN':>15}")
                    else:
                        print(f"{r['mesh_type']:<30} {r['n_elements']:<10} {Hz_pert:>15.6f}")

    # NGSolve reference values (placeholder - need to run NGSolve to get actual values)
    print("\n" + "="*70)
    print("NGSolve REFERENCE (from ELF/MAGIC)")
    print("="*70)
    print("""
To compare with NGSolve:
1. Run cube_ngsolve_reference_with_Kelvin.py to get reference values
2. Compare Hz_pert at same evaluation points
3. Agreement within ~1-5% validates the solution

Expected behavior:
- Hz_pert should be NEGATIVE on z-axis (above cube)
  because the magnetized cube creates a demagnetizing field
- Hz_pert should be POSITIVE on x-axis (beside cube)
  because the dipole field points in +z at the equator
""")

    print("\n" + "="*70)
    print("PRECISION REPORT")
    print("="*70)

    # Get hex 10x10x10 as reference
    hex_ref = None
    for r in results:
        if r and r['mesh_type'] == 'hexahedral' and r.get('n_div') == 10:
            hex_ref = r
            break

    if hex_ref:
        print(f"\nReference: Hexahedral 10x10x10 mesh")
        print("-" * 70)

        for ref_pt in ref_points:
            ref_Hz = None
            for fr in hex_ref['field_results']:
                if np.allclose(fr['point'], ref_pt, atol=0.01):
                    ref_Hz = fr['H_pert'][2]
                    break

            if ref_Hz is not None:
                print(f"\nPoint ({ref_pt[0]:.1f}, {ref_pt[1]:.1f}, {ref_pt[2]:.1f}): Hz_pert = {ref_Hz:.6f} A/m")
                print(f"{'Mesh Type':<30} {'Hz_pert [A/m]':>15} {'Error [%]':>12}")
                print("-" * 60)

                for r in results:
                    if r is None:
                        continue
                    for fr in r['field_results']:
                        if np.allclose(fr['point'], ref_pt, atol=0.01):
                            Hz = fr['H_pert'][2]
                            if np.isnan(Hz):
                                print(f"{r['mesh_type']:<30} {'NaN':>15} {'---':>12}")
                            else:
                                err = abs(Hz - ref_Hz) / abs(ref_Hz) * 100 if ref_Hz != 0 else 0
                                print(f"{r['mesh_type']:<30} {Hz:>15.6f} {err:>12.2f}")


if __name__ == "__main__":
    main()
