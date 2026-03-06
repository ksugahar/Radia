#!/usr/bin/env python
"""
Debug script for Tetra maxh=0.15 HACApK comparison with ELF_MAGIC.
"""

import sys
import os
import time
import math

_src_path = os.path.join(os.path.dirname(__file__), '../../../src/radia')
sys.path.insert(0, _src_path)

import radia as rad

try:
    import numpy as np
    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False

try:
    import psutil
    HAS_PSUTIL = True
except ImportError:
    HAS_PSUTIL = False


def get_peak_memory_mb():
    if not HAS_PSUTIL:
        return None
    process = psutil.Process(os.getpid())
    mem_info = process.memory_info()
    if hasattr(mem_info, 'peak_wset'):
        return mem_info.peak_wset / (1024 * 1024)
    else:
        return mem_info.rss / (1024 * 1024)


MU_0 = 4 * math.pi * 1e-7
H_EXT = 200000.0
MAXH = 0.15
HACAPK_EPS = 1e-6

# B-H curve data (same as ELF_MAGIC)
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


def get_M_avg_z(container):
    """Get average Mz."""
    all_M = rad.ObjM(container)
    M_list = [m[1] for m in all_M]
    if HAS_NUMPY:
        return float(np.mean([m[2] for m in M_list]))
    else:
        return sum(m[2] for m in M_list) / len(M_list) if M_list else 0.0


def run_solver_batch(solver_method, solver_name, max_iter=1000):
    """Run solver in batch mode (all iterations at once)."""
    try:
        from netgen.occ import Box, Pnt, OCCGeometry
        from ngsolve import Mesh
        from netgen_mesh_import import netgen_mesh_to_radia
    except ImportError as e:
        print(f'[SKIP] Tetrahedra: {e}')
        return None

    rad.FldUnits('m')
    rad.UtiDelAll()

    print(f'\n{"="*70}')
    print(f'{solver_name} Solver (BATCH) - maxh={MAXH}')
    print(f'{"="*70}')

    # Create mesh with Netgen
    t0 = time.time()
    cube_solid = Box(Pnt(-0.5, -0.5, -0.5), Pnt(0.5, 0.5, 0.5))
    cube_solid.mat('magnetic')
    geo = OCCGeometry(cube_solid)

    ngmesh = geo.GenerateMesh(maxh=MAXH)
    mesh = Mesh(ngmesh)
    n_elements = mesh.ne

    cube = netgen_mesh_to_radia(mesh,
                                 material={'magnetization': [0, 0, 0]},
                                 units='m',
                                 material_filter='magnetic',
                                 verbose=False)
    t_mesh = time.time() - t0
    print(f'Mesh creation: {t_mesh:.4f} s ({n_elements} elements, {n_elements*3} DOF)')

    # Apply nonlinear material
    mat = rad.MatSatIsoTab(BH_DATA)
    rad.MatApl(cube, mat)

    # External field
    B_ext = MU_0 * H_EXT
    ext = rad.ObjBckg([0, 0, B_ext])
    grp = rad.ObjCnt([cube, ext])

    # Configure HACApK if Method 2
    if solver_method == 2:
        rad.SetHACApKParams(HACAPK_EPS, 10, 2.0)
        print(f'HACApK: eps={HACAPK_EPS:.0e}, leaf_size=10, eta=2.0')

    # Run solver in batch
    t_start = time.time()
    res = rad.Solve(grp, 0.001, max_iter, solver_method)
    t_solve = time.time() - t_start

    M_avg_z = get_M_avg_z(cube)
    residual = res[0] if res[0] else 0.0
    n_iter = int(res[3]) if res[3] else 0
    peak_mem = get_peak_memory_mb()

    print(f'Solve time: {t_solve:.3f} s')
    print(f'Iterations: {n_iter}')
    print(f'Residual: {residual:.4e}')
    print(f'M_avg_z: {M_avg_z:.0f} A/m')
    if peak_mem:
        print(f'Peak memory: {peak_mem:.1f} MB')

    # Get H-matrix stats if available
    if solver_method == 2:
        try:
            stats = rad.GetHACApKStats()
            if stats:
                print(f'\nH-matrix stats:')
                print(f'  Leaves: {stats["n_leaves"]} (low-rank: {stats["n_lowrank"]}, dense: {stats["n_dense"]})')
                print(f'  Max rank: {stats["max_rank"]}')
                print(f'  Compression: {stats["compression"]:.4f}')
                print(f'  Build time: {stats["build_time"]:.4f} s')
        except:
            pass

    return {
        'n_elements': n_elements,
        't_solve': t_solve,
        'n_iter': n_iter,
        'residual': residual,
        'M_avg_z': M_avg_z,
        'peak_mem': peak_mem
    }


def main():
    print('=' * 70)
    print('TETRA maxh=0.15 DEBUG COMPARISON')
    print('=' * 70)
    print(f'H_ext = {H_EXT} A/m')
    print(f'maxh = {MAXH}')
    print(f'HACApK eps = {HACAPK_EPS:.0e}')

    # Run in BATCH mode (all iterations at once)
    result_lu = run_solver_batch(0, 'LU')
    result_hacapk = run_solver_batch(2, 'HACApK')

    # Summary
    print('\n' + '=' * 70)
    print('SUMMARY (BATCH MODE)')
    print('=' * 70)

    if result_lu and result_hacapk:
        print(f'{"Solver":<12} {"Elements":>8} {"Time (s)":>10} {"Iter":>6} {"M_avg_z":>12} {"Residual":>12}')
        print('-' * 70)
        print(f'{"LU":<12} {result_lu["n_elements"]:>8} {result_lu["t_solve"]:>10.3f} {result_lu["n_iter"]:>6} {result_lu["M_avg_z"]:>12.0f} {result_lu["residual"]:>12.4e}')
        print(f'{"HACApK":<12} {result_hacapk["n_elements"]:>8} {result_hacapk["t_solve"]:>10.3f} {result_hacapk["n_iter"]:>6} {result_hacapk["M_avg_z"]:>12.0f} {result_hacapk["residual"]:>12.4e}')

        # Differences
        diff = abs(result_hacapk['M_avg_z'] - result_lu['M_avg_z'])
        print(f'\nDifference from LU: {diff:.0f} A/m ({diff/result_lu["M_avg_z"]*100:.3f}%)')


if __name__ == '__main__':
    main()
