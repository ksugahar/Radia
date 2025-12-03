#!/usr/bin/env python
"""
Linear Material Benchmark: Tetrahedral vs Hexahedral Meshes

Compares Radia tetrahedral and hexahedral mesh results for a linear
magnetic material in uniform external field.

Evaluation points are OUTSIDE the cube (air region) - proper MMM validation.

Problem Setup:
- Cube: 0.1m x 0.1m x 0.1m, centered at origin
- Material: Linear isotropic, mu_r = 100
- External field: B0 = 1 T (uniform along z-axis)
- Evaluation: Points OUTSIDE the cube at various distances

Author: Radia Development Team
Date: 2025-12-03
"""
import sys
import os
# Add build/Release first for radia.pyd, then src/python for Python modules only
# (netgen_mesh_import.py, etc.)
_build_path = os.path.join(os.path.dirname(__file__), '../../../build/Release')
_src_path = os.path.join(os.path.dirname(__file__), '../../../src/python')
sys.path.insert(0, _build_path)
# Only add src/python AFTER build/Release so radia.pyd is loaded from build/Release
sys.path.append(_src_path)  # Use append, not insert!

import numpy as np
import time
import radia as rad

# Physical constants
MU_0 = 4 * np.pi * 1e-7  # T/(A/m)

# Problem parameters
CUBE_SIZE = 0.1      # 0.1 m cube (100 mm)
CUBE_HALF = 0.05     # half size
MU_R = 100           # Relative permeability
CHI = MU_R - 1       # Susceptibility = 99
B0_z = 1.0           # 1 T external field
B0 = np.array([0.0, 0.0, B0_z])


def create_evaluation_points():
    """Create evaluation points outside the cube."""
    points = []
    descriptions = []

    # Points on z-axis (above cube)
    for z in [0.06, 0.08, 0.10, 0.15, 0.20, 0.30]:
        points.append([0.0, 0.0, z])
        descriptions.append(f'z={z*1000:.0f}mm')

    # Points on x-axis (beside cube)
    for x in [0.06, 0.08, 0.10, 0.15, 0.20]:
        points.append([x, 0.0, 0.0])
        descriptions.append(f'x={x*1000:.0f}mm')

    return np.array(points), descriptions


def run_and_evaluate_hex(n_div, points):
    """Run hexahedral mesh and evaluate field."""
    rad.FldUnits('m')
    rad.UtiDelAll()

    # Create cube
    cube = rad.ObjRecMag([0, 0, 0], [CUBE_SIZE, CUBE_SIZE, CUBE_SIZE], [0, 0, 0])
    rad.ObjDivMag(cube, [n_div, n_div, n_div])
    n_elem = n_div ** 3

    # Apply material
    mat = rad.MatLin(CHI)
    rad.MatApl(cube, mat)

    # Background field (convert to Python float for Radia API)
    background = rad.ObjBckg([float(B0[0]), float(B0[1]), float(B0[2])])
    container = rad.ObjCnt([cube, background])

    # Solve
    start = time.time()
    res = rad.Solve(container, 0.0001, 10000, 9)  # Method 9
    solve_time = time.time() - start

    # Evaluate field at points
    B_values = []
    for pt in points:
        B = rad.Fld(container, 'b', list(pt))
        B_pert = np.array(B) - B0  # Subtract background
        B_values.append(B_pert)

    return {
        'n_elem': n_elem,
        'time': solve_time,
        'iter': res[3],
        'B_pert': np.array(B_values)
    }


def run_and_evaluate_tet(max_h, points):
    """Run tetrahedral mesh and evaluate field."""
    try:
        from netgen.occ import Box, OCCGeometry
        from ngsolve import Mesh
        from netgen_mesh_import import netgen_mesh_to_radia
    except ImportError as e:
        print(f"[SKIP] Tetrahedral: {e}")
        return None

    rad.FldUnits('m')
    rad.UtiDelAll()
    rad.SolverTetraMethod(0)

    # Create mesh
    cube_geom = Box((-CUBE_HALF, -CUBE_HALF, -CUBE_HALF),
                    (CUBE_HALF, CUBE_HALF, CUBE_HALF))
    geo = OCCGeometry(cube_geom)
    ngmesh = Mesh(geo.GenerateMesh(maxh=max_h))
    n_elem = ngmesh.ne

    # Import to Radia
    cube = netgen_mesh_to_radia(ngmesh,
                                 material={'magnetization': [0, 0, 0]},
                                 units='m',
                                 verbose=False)

    # Apply material
    mat = rad.MatLin(CHI)
    rad.MatApl(cube, mat)

    # Background field (convert to Python float for Radia API)
    background = rad.ObjBckg([float(B0[0]), float(B0[1]), float(B0[2])])
    container = rad.ObjCnt([cube, background])

    # Solve
    start = time.time()
    res = rad.Solve(container, 0.001, 10000, 9)  # Method 9
    solve_time = time.time() - start

    # Evaluate field at points
    B_values = []
    for pt in points:
        B = rad.Fld(container, 'b', list(pt))
        B_pert = np.array(B) - B0  # Subtract background
        B_values.append(B_pert)

    return {
        'n_elem': n_elem,
        'time': solve_time,
        'iter': res[3],
        'B_pert': np.array(B_values)
    }


def main():
    print("=" * 80)
    print("Linear Material Benchmark: Tetrahedral vs Hexahedral Meshes")
    print("=" * 80)

    print(f"\nProblem Setup:")
    print(f"  Cube size:   {CUBE_SIZE*1000:.0f} mm (centered at origin)")
    print(f"  mu_r:        {MU_R}")
    print(f"  chi:         {CHI}")
    print(f"  B0:          {B0_z:.1f} T (along z-axis)")
    print(f"  Solver:      Method 9 (Direct LU)")

    # Evaluation points
    points, descriptions = create_evaluation_points()

    print(f"\nEvaluation points ({len(points)} points outside cube):")
    for i, (pt, desc) in enumerate(zip(points[:5], descriptions[:5])):
        print(f"  {i+1}. {desc}: ({pt[0]*1000:.0f}, {pt[1]*1000:.0f}, {pt[2]*1000:.0f}) mm")
    print(f"  ... and {len(points)-5} more")

    results = {}

    # 1. Hexahedral meshes
    print("\n" + "=" * 80)
    print("HEXAHEDRAL MESHES (ObjDivMag)")
    print("=" * 80)

    for n_div in [4, 6, 8, 10]:
        key = f'hex_{n_div}x{n_div}x{n_div}'
        print(f"\n{key}:", end='')

        r = run_and_evaluate_hex(n_div, points)
        results[key] = r

        print(f" {r['n_elem']} elements, {r['time']:.3f}s")

    # 2. Tetrahedral meshes
    print("\n" + "=" * 80)
    print("TETRAHEDRAL MESHES (Netgen import)")
    print("=" * 80)

    for max_h in [0.08, 0.06, 0.04, 0.03]:
        key = f'tet_maxh{int(max_h*1000)}mm'
        print(f"\n{key}:", end='')

        r = run_and_evaluate_tet(max_h, points)
        if r is None:
            continue
        results[key] = r

        print(f" {r['n_elem']} elements, {r['time']:.3f}s")

    # 3. Use hex_10x10x10 as reference
    ref_key = 'hex_10x10x10'
    if ref_key not in results:
        ref_key = 'hex_8x8x8'

    B_ref = results[ref_key]['B_pert']

    # 4. Results table
    print("\n" + "=" * 80)
    print(f"RESULTS: Bz_pert at External Points (Reference: {ref_key})")
    print("=" * 80)

    # Header
    print(f"\n{'Point':<12}", end='')
    for key in results:
        short_name = key.replace('hex_', 'H').replace('tet_maxh', 'T').replace('x', '')
        print(f"{short_name:>12}", end='')
    print()
    print("-" * (12 + 12 * len(results)))

    # Data rows - Bz values in mT
    for i, desc in enumerate(descriptions):
        print(f"{desc:<12}", end='')
        for key in results:
            Bz = results[key]['B_pert'][i, 2]
            print(f"{Bz*1000:>12.4f}", end='')  # Convert to mT
        print()

    print("-" * (12 + 12 * len(results)))
    print("(Bz_pert in mT)")

    # 5. Error summary
    print("\n" + "=" * 80)
    print(f"ERROR SUMMARY (vs {ref_key})")
    print("=" * 80)

    print(f"\n{'Mesh Type':<20} {'Elements':>10} {'Avg Err':>10} {'Max Err':>10} {'Time':>10}")
    print("-" * 65)

    for key in results:
        if key == ref_key:
            continue

        B = results[key]['B_pert']
        errors = []
        for i in range(len(points)):
            Bz = B[i, 2]
            Bz_ref = B_ref[i, 2]
            if abs(Bz_ref) > 1e-10:
                err = abs(Bz - Bz_ref) / abs(Bz_ref) * 100
                errors.append(err)

        if errors:
            avg_err = np.mean(errors)
            max_err = np.max(errors)
        else:
            avg_err = max_err = 0

        n_elem = results[key]['n_elem']
        t = results[key]['time']
        print(f"{key:<20} {n_elem:>10} {avg_err:>9.2f}% {max_err:>9.2f}% {t:>9.3f}s")

    # Return for README update
    return results, points, descriptions


if __name__ == "__main__":
    results, points, descriptions = main()
