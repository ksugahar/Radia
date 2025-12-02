#!/usr/bin/env python3
"""
High Permeability Benchmark: Tetrahedra vs Hexahedra

Evaluate MMM accuracy for high permeability materials (mu_r = 100 to 4000).

Test case: Linear magnetic cube in uniform external field
- Cube: 100mm x 100mm x 100mm
- B0 = 1.0 T (z-direction)
- Validation: External field at points outside cube
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../build/Release'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../src/python'))

import numpy as np
import radia as rad
import time

MU_0 = 4 * np.pi * 1e-7
CUBE_SIZE = 0.1
B0 = 1.0


def uniform_B(pos):
    return [0, 0, B0]


OBS_POINTS = [
    ([0.10, 0, 0], '+x, r=0.10m'),
    ([0.12, 0, 0], '+x, r=0.12m'),
    ([0.15, 0, 0], '+x, r=0.15m'),
    ([0, 0, 0.10], '+z, r=0.10m'),
    ([0, 0, 0.12], '+z, r=0.12m'),
    ([0, 0, 0.15], '+z, r=0.15m'),
]


def test_hexahedral(mu_r, n_div=8, tol=0.001, max_iter=20000):
    chi = mu_r - 1.0
    rad.UtiDelAll()
    rad.FldUnits('m')

    cube = rad.ObjRecMag([0, 0, 0], [CUBE_SIZE]*3, [0, 0, 0])
    rad.ObjDivMag(cube, [n_div, n_div, n_div])
    mat = rad.MatLin(chi)
    rad.MatApl(cube, mat)

    bg = rad.ObjBckgCF(uniform_B)
    system = rad.ObjCnt([cube, bg])

    t0 = time.time()
    result = rad.Solve(system, tol, max_iter)
    solve_time = time.time() - t0

    fields = {}
    for pt, desc in OBS_POINTS:
        H = np.array(rad.Fld(system, 'h', pt))
        fields[tuple(pt)] = H

    return {
        'n_div': n_div,
        'n_elements': n_div**3,
        'solve_time': solve_time,
        'iterations': result[3],
        'max_change': result[1],
        'converged': result[3] < max_iter - 1,
        'fields': fields
    }


def test_tetra_manual(mu_r, n_subdivs=2, tol=0.001, max_iter=20000):
    chi = mu_r - 1.0
    rad.UtiDelAll()
    rad.FldUnits('m')
    rad.SolverTetraMethod(0)

    TETRA_FACES = [[1, 3, 2], [1, 2, 4], [2, 3, 4], [3, 1, 4]]
    hs = CUBE_SIZE / 2

    def create_cell_tetras(cx, cy, cz, size):
        h = size / 2
        v = [
            [cx - h, cy - h, cz - h], [cx + h, cy - h, cz - h],
            [cx - h, cy + h, cz - h], [cx + h, cy + h, cz - h],
            [cx - h, cy - h, cz + h], [cx + h, cy - h, cz + h],
            [cx - h, cy + h, cz + h], [cx + h, cy + h, cz + h],
        ]
        indices = [[0,1,3,5], [0,3,2,6], [0,5,4,6], [3,5,6,7], [0,3,5,6]]
        return [[v[i] for i in idx] for idx in indices]

    cell_size = CUBE_SIZE / n_subdivs
    all_objs = []

    for ix in range(n_subdivs):
        for iy in range(n_subdivs):
            for iz in range(n_subdivs):
                cx = -hs + cell_size/2 + ix * cell_size
                cy = -hs + cell_size/2 + iy * cell_size
                cz = -hs + cell_size/2 + iz * cell_size
                for verts in create_cell_tetras(cx, cy, cz, cell_size):
                    obj = rad.ObjPolyhdr(verts, TETRA_FACES, [0, 0, 0])
                    all_objs.append(obj)

    n_elements = len(all_objs)
    cube = rad.ObjCnt(all_objs)
    mat = rad.MatLin(chi)
    rad.MatApl(cube, mat)

    bg = rad.ObjBckgCF(uniform_B)
    system = rad.ObjCnt([cube, bg])

    t0 = time.time()
    result = rad.Solve(system, tol, max_iter)
    solve_time = time.time() - t0

    fields = {}
    valid = True
    for pt, desc in OBS_POINTS:
        H = np.array(rad.Fld(system, 'h', pt))
        if np.any(np.isnan(H)):
            valid = False
        fields[tuple(pt)] = H

    return {
        'n_subdivs': n_subdivs,
        'n_elements': n_elements,
        'solve_time': solve_time,
        'iterations': result[3],
        'max_change': result[1],
        'converged': result[3] < max_iter - 1 and valid,
        'fields': fields,
        'valid': valid
    }


def test_tetra_netgen(mu_r, maxh=0.05, tol=0.001, max_iter=20000):
    try:
        from netgen.occ import Box, Pnt, OCCGeometry
        from ngsolve import Mesh
        from netgen_mesh_import import netgen_mesh_to_radia
    except ImportError:
        return None

    chi = mu_r - 1.0
    rad.UtiDelAll()
    rad.FldUnits('m')
    rad.SolverTetraMethod(0)

    hs = CUBE_SIZE / 2
    cube_solid = Box(Pnt(-hs, -hs, -hs), Pnt(hs, hs, hs))
    cube_solid.mat('magnetic')
    geo = OCCGeometry(cube_solid)

    ngmesh = geo.GenerateMesh(maxh=maxh)
    mesh = Mesh(ngmesh)
    n_elements = mesh.ne

    cube = netgen_mesh_to_radia(mesh, material={'magnetization': [0, 0, 0]},
                                 units='m', material_filter='magnetic', verbose=False)

    mat = rad.MatLin(chi)
    rad.MatApl(cube, mat)

    bg = rad.ObjBckgCF(uniform_B)
    system = rad.ObjCnt([cube, bg])

    t0 = time.time()
    result = rad.Solve(system, tol, max_iter)
    solve_time = time.time() - t0

    fields = {}
    valid = True
    for pt, desc in OBS_POINTS:
        H = np.array(rad.Fld(system, 'h', pt))
        if np.any(np.isnan(H)):
            valid = False
        fields[tuple(pt)] = H

    return {
        'maxh': maxh,
        'n_elements': n_elements,
        'solve_time': solve_time,
        'iterations': result[3],
        'max_change': result[1],
        'converged': result[3] < max_iter - 1 and valid,
        'fields': fields,
        'valid': valid
    }


def compute_error(test_fields, ref_fields):
    """Compute average and max error vs reference."""
    errors = []
    for pt, desc in OBS_POINTS:
        H_ref = ref_fields[tuple(pt)]
        H_test = test_fields[tuple(pt)]
        err = np.linalg.norm(H_test - H_ref) / np.linalg.norm(H_ref) * 100
        errors.append(err)
    return np.mean(errors), max(errors)


def main():
    print("=" * 80)
    print("HIGH PERMEABILITY BENCHMARK (tol=0.001, max_iter=20000)")
    print("=" * 80)
    print(f"Cube: {CUBE_SIZE*1000:.0f} mm, B0 = {B0} T")
    print()

    results_table = []
    mu_r_list = [100, 500, 1000, 2000, 4000]

    for mu_r in mu_r_list:
        chi = mu_r - 1
        print(f"\nmu_r = {mu_r} (chi = {chi})")
        print("-" * 70)

        row = {'mu_r': mu_r}

        # Hexahedral reference (fine mesh)
        hex_ref = test_hexahedral(mu_r, 8)
        status = 'OK' if hex_ref['converged'] else f"NOT CONV ({hex_ref['max_change']:.1e})"
        print(f"  Hex 8x8x8 (512): {hex_ref['iterations']:.0f} iter, {hex_ref['solve_time']:.2f}s, {status}")
        row['hex_8'] = hex_ref

        # Coarser hex meshes
        for n_div in [4, 6]:
            hex_test = test_hexahedral(mu_r, n_div)
            status = 'OK' if hex_test['converged'] else 'NOT CONV'

            if hex_ref['converged'] and hex_test['converged']:
                avg_err, max_err = compute_error(hex_test['fields'], hex_ref['fields'])
                print(f"  Hex {n_div}x{n_div}x{n_div} ({n_div**3:3d}): {hex_test['iterations']:.0f} iter, "
                      f"{hex_test['solve_time']:.2f}s, {status}, err={avg_err:.2f}%")
            else:
                print(f"  Hex {n_div}x{n_div}x{n_div} ({n_div**3:3d}): {hex_test['iterations']:.0f} iter, "
                      f"{hex_test['solve_time']:.2f}s, {status}")
            row[f'hex_{n_div}'] = hex_test

        # Tetra manual
        for n_sub in [1, 2]:
            tet = test_tetra_manual(mu_r, n_sub)
            if tet['valid']:
                status = 'OK' if tet['converged'] else f"NOT CONV ({tet['max_change']:.1e})"
                if hex_ref['converged']:
                    avg_err, max_err = compute_error(tet['fields'], hex_ref['fields'])
                    print(f"  Tet {n_sub}^3x5 ({tet['n_elements']:3d}): {tet['iterations']:.0f} iter, "
                          f"{tet['solve_time']:.2f}s, {status}, err={avg_err:.2f}%")
                else:
                    print(f"  Tet {n_sub}^3x5 ({tet['n_elements']:3d}): {tet['iterations']:.0f} iter, "
                          f"{tet['solve_time']:.2f}s, {status}")
            else:
                print(f"  Tet {n_sub}^3x5 ({tet['n_elements']:3d}): {tet['iterations']:.0f} iter, NaN")
            row[f'tet_manual_{n_sub}'] = tet

        # Tetra Netgen
        for maxh in [0.08, 0.06]:
            tet = test_tetra_netgen(mu_r, maxh)
            if tet and tet['valid']:
                status = 'OK' if tet['converged'] else f"NOT CONV ({tet['max_change']:.1e})"
                if hex_ref['converged']:
                    avg_err, max_err = compute_error(tet['fields'], hex_ref['fields'])
                    print(f"  Tet NG h={maxh} ({tet['n_elements']:3d}): {tet['iterations']:.0f} iter, "
                          f"{tet['solve_time']:.2f}s, {status}, err={avg_err:.2f}%")
                else:
                    print(f"  Tet NG h={maxh} ({tet['n_elements']:3d}): {tet['iterations']:.0f} iter, "
                          f"{tet['solve_time']:.2f}s, {status}")
            elif tet:
                print(f"  Tet NG h={maxh} ({tet['n_elements']:3d}): {tet['iterations']:.0f} iter, NaN")
            row[f'tet_ng_{maxh}'] = tet

        results_table.append(row)

    # Summary table
    print()
    print("=" * 80)
    print("SUMMARY: External Field Error vs Hex 8x8x8 Reference")
    print("=" * 80)
    header = f"{'mu_r':<8} {'Hex 4x4x4':<12} {'Hex 6x6x6':<12} {'Tet 1^3x5':<12} {'Tet 2^3x5':<12} {'Tet NG 28':<12}"
    print(header)
    print("-" * 70)

    for row in results_table:
        mu_r = row['mu_r']
        hex_ref = row['hex_8']

        cols = [f"{mu_r}"]

        for key in ['hex_4', 'hex_6', 'tet_manual_1', 'tet_manual_2', 'tet_ng_0.06']:
            if key not in row or row[key] is None:
                cols.append('N/A')
                continue

            test = row[key]
            if not test.get('valid', True):
                cols.append('NaN')
            elif not hex_ref['converged'] or not test.get('converged', False):
                cols.append('N/C')
            else:
                avg_err, _ = compute_error(test['fields'], hex_ref['fields'])
                cols.append(f"{avg_err:.2f}%")

        print(f"{cols[0]:<8} {cols[1]:<12} {cols[2]:<12} {cols[3]:<12} {cols[4]:<12} {cols[5]:<12}")

    print("=" * 80)
    print()
    print("Legend: N/A = Not available, N/C = Not converged, NaN = Numerical failure")
    print()
    print("Notes:")
    print("  - High permeability (mu_r > 100) requires more iterations")
    print("  - Tet 1^3x5 (5 elements) converges fast but has coarse approximation")
    print("  - Finer meshes may not converge within iteration limit")


if __name__ == '__main__':
    main()
