#!/usr/bin/env python
"""Benchmark: Cube model with Block Jacobi + Newton.

Re-runs the cube uniform field benchmark (Table 1 in SA-26-010)
with the new Block Jacobi preconditioner and Newton-Raphson.

Usage: python -u benchmark_cube_block_jacobi.py [5|10|15|20|25|all]
"""
import sys, os, time
import numpy as np

repo_root = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')
sys.path.insert(0, os.path.join(repo_root, 'src'))

import radia as rad

MU_0 = 4 * np.pi * 1e-7
H_EXT = 200000.0  # A/m

# Same BH curve as benchmark_common.py
BH_DATA = [
    [0.0, 0.0], [100.0, 0.1], [200.0, 0.3], [500.0, 0.8], [1000.0, 1.2],
    [2000.0, 1.5], [5000.0, 1.7], [10000.0, 1.8], [50000.0, 2.0], [100000.0, 2.1],
]

def generate_hex_mesh(n_div, size=1.0):
    vertices_list = []
    dx = size / n_div
    offset = size / 2
    for iz in range(n_div):
        for iy in range(n_div):
            for ix in range(n_div):
                x0 = ix * dx - offset
                y0 = iy * dx - offset
                z0 = iz * dx - offset
                verts = [
                    [x0, y0, z0],
                    [x0 + dx, y0, z0],
                    [x0 + dx, y0 + dx, z0],
                    [x0, y0 + dx, z0],
                    [x0, y0, z0 + dx],
                    [x0 + dx, y0, z0 + dx],
                    [x0 + dx, y0 + dx, z0 + dx],
                    [x0, y0 + dx, z0 + dx]
                ]
                vertices_list.append(verts)
    return vertices_list

def run_cube_benchmark(n_div, solver_method, solver_name, use_newton=False, max_iter=100):
    rad.UtiDelAll()

    # Create mesh
    t0 = time.time()
    mesh = generate_hex_mesh(n_div)
    n_elem = len(mesh)
    n_dof = n_elem * 6

    mat = rad.MatSatIsoTab(BH_DATA)
    all_objects = []
    for verts in mesh:
        obj = rad.ObjHexahedron(verts, [0, 0, 0])
        rad.MatApl(obj, mat)
        all_objects.append(obj)
    cube = rad.ObjCnt(all_objects)

    # External field
    B_ext = MU_0 * H_EXT
    ext = rad.ObjBckg(lambda p: [0, 0, B_ext])
    grp = rad.ObjCnt([cube, ext])
    t_mesh = time.time() - t0

    if solver_method == 2:
        rad.SolverConfig(hacapk_eps=1e-4, hacapk_leaf=10, hacapk_eta=2.0)

    rad.SolverConfig(newton_method=use_newton)

    t0 = time.time()
    try:
        result = rad.Solve(grp, 0.001, max_iter, solver_method)
        converged = True
    except RuntimeError as e:
        converged = False
        print(f"  WARNING: {e}", flush=True)
    t_solve = time.time() - t0

    # Get magnetization
    all_M = rad.ObjM(cube)
    M_total_z = sum(m[1][2] for m in all_M)
    M_avg_z = M_total_z / n_elem

    stats = rad.GetSolveStats()

    hacapk_stats = None
    if solver_method == 2:
        try:
            config = rad.GetSolverConfig()
            hacapk_stats = config.get('hacapk_stats')
        except:
            pass

    mem_mb = 'N/A'
    compression = 'N/A'
    if hacapk_stats:
        mem_mb = f"{hacapk_stats.get('memory_mb', 0):.0f}"
        comp = hacapk_stats.get('compression', 0.0)
        compression = f"{comp*100:.0f}%"
    elif solver_method in [0, 1]:
        mem_mb = f"{n_dof * n_dof * 8 / 1024 / 1024:.0f}"

    rad.SolverConfig(newton_method=False)

    return {
        'name': solver_name,
        'n_div': n_div,
        'n_elem': n_elem,
        'n_dof': n_dof,
        'M_avg_z': M_avg_z,
        't_solve': t_solve,
        't_mesh': t_mesh,
        'nonl_iter': stats.get('nonl_iterations', 'N/A') if stats else 'N/A',
        'linear_iter': stats.get('linear_iterations', 'N/A') if stats else 'N/A',
        'converged': converged,
        'mem_mb': mem_mb,
        'compression': compression,
    }

def print_result(r):
    status = "CONVERGED" if r['converged'] else "NOT CONVERGED"
    print(f"  {r['name']}:", flush=True)
    print(f"    M_avg_z = {r['M_avg_z']:.0f} A/m", flush=True)
    print(f"    Time = {r['t_solve']:.1f} s (mesh: {r['t_mesh']:.1f} s)", flush=True)
    print(f"    NL iter = {r['nonl_iter']}", flush=True)
    print(f"    Linear iter = {r['linear_iter']}", flush=True)
    print(f"    Memory = {r['mem_mb']} MB", flush=True)
    print(f"    Compression = {r['compression']}", flush=True)
    print(f"    Status: {status}", flush=True)
    print(flush=True)


# --- Main ---
mesh_arg = sys.argv[1] if len(sys.argv) > 1 else 'all'

if mesh_arg == 'all':
    sizes = [5, 10, 15, 20, 25]
else:
    sizes = [int(mesh_arg)]

results_all = {}

for n in sizes:
    n_elem = n**3
    n_dof = n_elem * 6
    print(f"{'='*60}", flush=True)
    print(f"N={n} (Elements: {n_elem}, DOF: {n_dof})", flush=True)
    print(f"{'='*60}", flush=True)

    results = []

    # Determine which solvers to run based on problem size
    configs = []

    if n <= 15:
        configs.append((0, 'LU', False, 100))
    if n <= 20:
        configs.append((1, 'BiCGSTAB (dense)', False, 100))

    configs.append((2, 'HACApK + Block Jacobi', False, 100))
    configs.append((2, 'HACApK + Block Jacobi + Newton', True, 200))

    for method, name, newton, max_iter in configs:
        print(f"Running {name}...", flush=True)
        r = run_cube_benchmark(n, method, name, use_newton=newton, max_iter=max_iter)
        print_result(r)
        results.append(r)

    results_all[f'N={n}'] = results

# Summary
print(f"\n{'='*60}", flush=True)
print("SUMMARY", flush=True)
print(f"{'='*60}", flush=True)
for mesh, results in results_all.items():
    print(f"\n--- {mesh} ---", flush=True)
    for r in results:
        status = "OK" if r['converged'] else "FAIL"
        print(f"  {r['name']:40s} Mz={r['M_avg_z']:10.0f}  "
              f"Time={r['t_solve']:8.1f}s  NL={str(r['nonl_iter']):>4s}  "
              f"Lin={str(r['linear_iter']):>6s}  Mem={r['mem_mb']:>8s} MB  "
              f"Comp={r['compression']:>6s}  [{status}]", flush=True)

print("\nDone.", flush=True)
