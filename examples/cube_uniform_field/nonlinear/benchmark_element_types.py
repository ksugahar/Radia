#!/usr/bin/env python
"""
Benchmark: Element Type Comparison for Nonlinear Material

Compares three element types in Radia for nonlinear material solving:
1. Tetrahedra (ObjPolyhdr) - using Netgen mesh
2. Hexahedra (ObjRecMag + ObjDivMag)
3. ThickPgn (ObjThckPgn) - extruded polygons

Problem Setup:
- Cube: 1.0m x 1.0m x 1.0m centered at origin
- Material: Nonlinear B-H curve (soft iron)
- External field: Hz = 50,000 A/m

Author: Radia Development Team
Date: 2025-12-05
"""

import sys
import os
import time
import json
import argparse

# Path setup
_build_path = os.path.join(os.path.dirname(__file__), '../../../build/Release')
_src_path = os.path.join(os.path.dirname(__file__), '../../../src/python')
sys.path.insert(0, _build_path)
sys.path.append(_src_path)

import numpy as np
import radia as rad

# Physical constants
MU_0 = 4 * np.pi * 1e-7  # T/(A/m)

# Problem parameters
CUBE_SIZE = 1.0      # 1.0 m cube
CUBE_HALF = 0.5      # half size
H_EXT = 50000.0      # External field (A/m)
B_EXT = MU_0 * H_EXT  # External B field (T)

# B-H curve data: [H (A/m), B (T)]
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

# Convert B-H to H-M format: M = B/mu_0 - H
HM_DATA = [[h, b/MU_0 - h] for h, b in BH_DATA]


def benchmark_hexahedra(n_div, solver_method=1):
    """Benchmark hexahedral mesh (ObjRecMag + ObjDivMag)."""
    rad.FldUnits('m')
    rad.UtiDelAll()

    element_type = 'hexa'
    n_elements = n_div ** 3

    print(f'\n{"="*70}')
    print(f'HEXAHEDRAL MESH: {n_div}x{n_div}x{n_div} = {n_elements} elements')
    print(f'{"="*70}')

    # Create mesh
    t_mesh_start = time.time()
    cube = rad.ObjRecMag([0, 0, 0], [CUBE_SIZE, CUBE_SIZE, CUBE_SIZE], [0, 0, 0])
    rad.ObjDivMag(cube, [n_div, n_div, n_div])
    t_mesh = time.time() - t_mesh_start

    # Apply nonlinear material
    mat = rad.MatSatIsoTab(HM_DATA)
    rad.MatApl(cube, mat)

    # External field
    ext = rad.ObjBckg([0, 0, B_EXT])
    grp = rad.ObjCnt([cube, ext])

    # Solve
    t_solve_start = time.time()
    result = rad.Solve(grp, 0.001, 1000, solver_method)
    t_solve = time.time() - t_solve_start

    # Get magnetization
    all_M = rad.ObjM(cube)
    M_list = [m[1] for m in all_M]
    M_avg_z = np.mean([m[2] for m in M_list])

    n_iter = int(result[3]) if result[3] else 0
    converged = n_iter < 1000
    residual = result[0] if result[0] else 0.0

    print(f'Mesh time:    {t_mesh:.4f} s')
    print(f'Solve time:   {t_solve:.3f} s')
    print(f'Iterations:   {n_iter}')
    print(f'Converged:    {converged}')
    print(f'M_avg_z:      {M_avg_z:.0f} A/m')

    return {
        'element_type': element_type,
        'mesh_description': f'{n_div}x{n_div}x{n_div}',
        'n_elements': n_elements,
        'ndof': n_elements * 3,
        'H_ext': H_EXT,
        't_mesh': t_mesh,
        't_solve': t_solve,
        'solver_method': solver_method,
        'converged': converged,
        'residual': residual,
        'nonl_iterations': n_iter,
        'M_avg_z': M_avg_z,
    }


def benchmark_tetrahedra(maxh=0.3, solver_method=1):
    """Benchmark tetrahedral mesh (Netgen + ObjPolyhdr)."""
    try:
        from netgen.occ import Box, Pnt, OCCGeometry
        from ngsolve import Mesh
        from netgen_mesh_import import netgen_mesh_to_radia
    except ImportError as e:
        print(f'[SKIP] Tetrahedra: {e}')
        return None

    rad.FldUnits('m')
    rad.UtiDelAll()
    rad.SolverTetraMethod(0)  # Original method

    element_type = 'tetra'

    print(f'\n{"="*70}')
    print(f'TETRAHEDRAL MESH: maxh={maxh}m')
    print(f'{"="*70}')

    # Create mesh with Netgen
    t_mesh_start = time.time()
    cube_solid = Box(Pnt(-CUBE_HALF, -CUBE_HALF, -CUBE_HALF),
                      Pnt(CUBE_HALF, CUBE_HALF, CUBE_HALF))
    cube_solid.mat('magnetic')
    geo = OCCGeometry(cube_solid)

    ngmesh = geo.GenerateMesh(maxh=maxh)
    mesh = Mesh(ngmesh)
    n_elements = mesh.ne

    cube = netgen_mesh_to_radia(mesh,
                                 material={'magnetization': [0, 0, 0]},
                                 units='m',
                                 material_filter='magnetic',
                                 verbose=False)
    t_mesh = time.time() - t_mesh_start

    print(f'Generated {n_elements} tetrahedral elements')

    # Apply nonlinear material
    mat = rad.MatSatIsoTab(HM_DATA)
    rad.MatApl(cube, mat)

    # External field
    ext = rad.ObjBckg([0, 0, B_EXT])
    grp = rad.ObjCnt([cube, ext])

    # Solve
    t_solve_start = time.time()
    try:
        result = rad.Solve(grp, 0.001, 1000, solver_method)
        t_solve = time.time() - t_solve_start

        # Get magnetization
        all_M = rad.ObjM(cube)
        M_list = [m[1] for m in all_M]
        M_avg_z = np.mean([m[2] for m in M_list])

        n_iter = int(result[3]) if result[3] else 0
        converged = n_iter < 1000 and not np.isnan(M_avg_z)
        residual = result[0] if result[0] else 0.0
    except Exception as e:
        print(f'Solve failed: {e}')
        return None

    print(f'Mesh time:    {t_mesh:.4f} s')
    print(f'Solve time:   {t_solve:.3f} s')
    print(f'Iterations:   {n_iter}')
    print(f'Converged:    {converged}')
    print(f'M_avg_z:      {M_avg_z:.0f} A/m')

    return {
        'element_type': element_type,
        'mesh_description': f'maxh={maxh}m',
        'n_elements': n_elements,
        'ndof': n_elements * 3,
        'H_ext': H_EXT,
        't_mesh': t_mesh,
        't_solve': t_solve,
        'solver_method': solver_method,
        'converged': converged,
        'residual': residual,
        'nonl_iterations': n_iter,
        'M_avg_z': M_avg_z,
    }


def benchmark_thickpgn(n_xy, n_z, solver_method=1):
    """
    Benchmark thick polygon mesh (ObjThckPgn).

    Creates a cube by stacking square prisms along z-axis,
    each prism created using ObjThckPgn.
    """
    rad.FldUnits('m')
    rad.UtiDelAll()

    element_type = 'thickpgn'

    print(f'\n{"="*70}')
    print(f'THICK POLYGON MESH: {n_xy}x{n_xy}x{n_z} extruded squares')
    print(f'{"="*70}')

    # Create mesh
    t_mesh_start = time.time()

    cell_xy = CUBE_SIZE / n_xy  # xy dimension of each cell
    cell_z = CUBE_SIZE / n_z    # z thickness of each cell

    all_objs = []

    for iz in range(n_z):
        z_center = -CUBE_HALF + cell_z/2 + iz * cell_z

        for ix in range(n_xy):
            for iy in range(n_xy):
                # Center of this cell in xy plane
                x_center = -CUBE_HALF + cell_xy/2 + ix * cell_xy
                y_center = -CUBE_HALF + cell_xy/2 + iy * cell_xy

                # Define square polygon (2D points in yz plane for 'x' extrusion)
                # Or xy plane for 'z' extrusion
                # Using 'z' extrusion: points are in (x, y) plane
                half = cell_xy / 2
                points = [
                    [x_center - half, y_center - half],
                    [x_center + half, y_center - half],
                    [x_center + half, y_center + half],
                    [x_center - half, y_center + half],
                ]

                # Create thick polygon extruded along z
                # ObjThckPgn(z_center, thickness, points, axis='z', magnetization)
                obj = rad.ObjThckPgn(z_center, cell_z, points, 'z', [0, 0, 0])
                all_objs.append(obj)

    n_elements = len(all_objs)
    cube = rad.ObjCnt(all_objs)
    t_mesh = time.time() - t_mesh_start

    print(f'Generated {n_elements} thick polygon elements')

    # Apply nonlinear material
    mat = rad.MatSatIsoTab(HM_DATA)
    rad.MatApl(cube, mat)

    # External field
    ext = rad.ObjBckg([0, 0, B_EXT])
    grp = rad.ObjCnt([cube, ext])

    # Solve
    t_solve_start = time.time()
    result = rad.Solve(grp, 0.001, 1000, solver_method)
    t_solve = time.time() - t_solve_start

    # Get magnetization
    all_M = rad.ObjM(cube)
    M_list = [m[1] for m in all_M]
    M_avg_z = np.mean([m[2] for m in M_list])

    n_iter = int(result[3]) if result[3] else 0
    converged = n_iter < 1000
    residual = result[0] if result[0] else 0.0

    print(f'Mesh time:    {t_mesh:.4f} s')
    print(f'Solve time:   {t_solve:.3f} s')
    print(f'Iterations:   {n_iter}')
    print(f'Converged:    {converged}')
    print(f'M_avg_z:      {M_avg_z:.0f} A/m')

    return {
        'element_type': element_type,
        'mesh_description': f'{n_xy}x{n_xy}x{n_z}',
        'n_elements': n_elements,
        'ndof': n_elements * 3,
        'H_ext': H_EXT,
        't_mesh': t_mesh,
        't_solve': t_solve,
        'solver_method': solver_method,
        'converged': converged,
        'residual': residual,
        'nonl_iterations': n_iter,
        'M_avg_z': M_avg_z,
    }


def save_result(result, output_dir):
    """Save result to JSON file in appropriate folder."""
    if result is None:
        return

    elem_type = result['element_type']
    mesh_desc = result['mesh_description'].replace('=', '').replace('.', '_')
    filename = f'{elem_type}_{mesh_desc}_results.json'

    filepath = os.path.join(output_dir, elem_type, filename)
    os.makedirs(os.path.dirname(filepath), exist_ok=True)

    with open(filepath, 'w') as f:
        json.dump(result, f, indent=2)

    print(f'Saved: {filepath}')
    return filepath


def main():
    parser = argparse.ArgumentParser(description='Benchmark element types for nonlinear material')
    parser.add_argument('--hexa', action='store_true', help='Run hexahedral benchmarks')
    parser.add_argument('--tetra', action='store_true', help='Run tetrahedral benchmarks')
    parser.add_argument('--thickpgn', action='store_true', help='Run thick polygon benchmarks')
    parser.add_argument('--all', action='store_true', help='Run all benchmarks')
    parser.add_argument('--method', type=int, default=1, help='Solver method (0=LU, 1=BiCGSTAB)')
    args = parser.parse_args()

    if args.all:
        args.hexa = args.tetra = args.thickpgn = True

    if not (args.hexa or args.tetra or args.thickpgn):
        args.all = args.hexa = args.tetra = args.thickpgn = True

    output_dir = os.path.join(os.path.dirname(__file__), 'results')

    print('='*70)
    print('ELEMENT TYPE BENCHMARK: Nonlinear Material')
    print('='*70)
    print(f'Cube size:    {CUBE_SIZE} m')
    print(f'H_ext:        {H_EXT} A/m')
    print(f'Solver:       Method {args.method} ({"LU" if args.method == 0 else "BiCGSTAB"})')

    results = []

    # Hexahedral benchmarks
    if args.hexa:
        for n_div in [6, 8, 10]:
            r = benchmark_hexahedra(n_div, args.method)
            save_result(r, output_dir)
            results.append(r)

    # Tetrahedral benchmarks
    if args.tetra:
        for maxh in [0.4, 0.3, 0.25]:
            r = benchmark_tetrahedra(maxh, args.method)
            save_result(r, output_dir)
            if r:
                results.append(r)

    # Thick polygon benchmarks
    if args.thickpgn:
        for n_xy, n_z in [(6, 6), (8, 8), (10, 10)]:
            r = benchmark_thickpgn(n_xy, n_z, args.method)
            save_result(r, output_dir)
            results.append(r)

    # Summary
    print('\n' + '='*70)
    print('SUMMARY')
    print('='*70)
    print(f'{"Type":<10} {"Mesh":<15} {"Elements":>10} {"Time (s)":>10} {"M_avg_z":>12}')
    print('-'*70)

    for r in results:
        if r:
            print(f'{r["element_type"]:<10} {r["mesh_description"]:<15} '
                  f'{r["n_elements"]:>10} {r["t_solve"]:>10.3f} {r["M_avg_z"]:>12.0f}')

    print('='*70)
    return results


if __name__ == '__main__':
    main()
