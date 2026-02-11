#!/usr/bin/env python
"""Benchmark: Block Jacobi + Newton for C-type electromagnet.

Compares:
  - LU (method 0)
  - Dense BiCGSTAB (method 1)
  - HACApK + Block Jacobi (method 2, Picard)
  - HACApK + Block Jacobi + Newton hybrid (method 2, Newton)

Usage: python -u benchmark_block_jacobi_newton.py [6|10|20|all]
"""
import sys, os, time
import numpy as np

work_dir = r"S:\Radia\01_GitHub\examples\electromagnet\nonlinear"
repo_root = r"S:\Radia\01_GitHub"
sys.path.insert(0, os.path.join(repo_root, 'src'))
sys.path.insert(0, work_dir)

import radia as rad
from coil_model import create_racetrack_coil

scale = 0.001

def load_elf_geometry(path):
    nodes = {}
    hex_elements = []
    with open(os.path.join(path, "ELF_magic.meg"), 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('MGR1'):
                parts = line.split()
                node_id = int(parts[1])
                x, y, z = float(parts[3]), float(parts[4]), float(parts[5])
                nodes[node_id] = np.array([x, y, z])
            elif line.startswith('MMB8T'):
                parts = line.split()
                elem_id = int(parts[1])
                node_ids = [int(parts[i]) for i in range(4, 12)]
                hex_elements.append((elem_id, node_ids))
    return nodes, hex_elements

def load_bh_curve(filepath):
    data = np.loadtxt(filepath, comments='#')
    return data.tolist()

def run_benchmark(nodes, hex_elements, bh_data, solver_method, solver_name,
                  use_newton=False, max_iter=100):
    rad.UtiDelAll()
    rad.FldUnits('m')
    mat = rad.MatSatIsoTab(bh_data)
    all_objects = []
    for elem_id, node_ids in hex_elements:
        verts = [[nodes[nid][0]*scale, nodes[nid][1]*scale, nodes[nid][2]*scale]
                 for nid in node_ids]
        obj = rad.ObjHexahedron(verts, [0, 0, 0])
        rad.MatApl(obj, mat)
        all_objects.append(obj)
    yoke = rad.ObjCnt(all_objects)
    coil = create_racetrack_coil(20000.0)
    model = rad.ObjCnt([yoke, coil])

    if solver_method == 2:
        rad.SetHACApKParams(1e-4, 10, 2.0)

    rad.SetNewtonMethod(use_newton)

    t0 = time.time()
    try:
        result = rad.Solve(model, 0.001, max_iter, solver_method, image='+x-z')
        converged = True
    except RuntimeError as e:
        converged = False
        print(f"  WARNING: {e}", flush=True)
    t_solve = time.time() - t0

    B = rad.Fld(model, 'b', [0, 0, 0])
    stats = rad.GetSolveStats()

    hacapk_stats = None
    if solver_method == 2:
        try:
            hacapk_stats = rad.GetHACApKStats()
        except:
            pass

    mem_mb = 'N/A'
    if hacapk_stats:
        mem_mb = f"{hacapk_stats.get('memory_mb', 'N/A'):.1f}"
    elif solver_method in [0, 1]:
        # Dense matrix memory estimate: N*N*8 bytes
        n = len(hex_elements) * 6
        mem_mb = f"{n*n*8 / 1024 / 1024:.0f}"

    compression = 'N/A'
    if hacapk_stats and 'compression_ratio' in hacapk_stats:
        compression = f"{hacapk_stats['compression_ratio']*100:.1f}%"

    rad.SetNewtonMethod(False)

    return {
        'name': solver_name,
        'Bz_mT': B[2] * 1000,
        't_solve': t_solve,
        'nonl_iter': stats.get('nonl_iterations', 'N/A') if stats else 'N/A',
        'linear_iter': stats.get('linear_iterations', 'N/A') if stats else 'N/A',
        'converged': converged,
        'mem_mb': mem_mb,
        'compression': compression,
    }

def print_result(r):
    status = "CONVERGED" if r['converged'] else "NOT CONVERGED"
    print(f"  {r['name']}:", flush=True)
    print(f"    Bz = {r['Bz_mT']:.2f} mT", flush=True)
    print(f"    Time = {r['t_solve']:.1f} s", flush=True)
    print(f"    NL iter = {r['nonl_iter']}", flush=True)
    print(f"    Linear iter = {r['linear_iter']}", flush=True)
    print(f"    Memory = {r['mem_mb']} MB", flush=True)
    print(f"    Compression = {r['compression']}", flush=True)
    print(f"    Status: {status}", flush=True)
    print(flush=True)


# --- Main ---
bh_file = os.path.join(work_dir, "BH.txt")
bh_data = load_bh_curve(bh_file)

mesh_arg = sys.argv[1] if len(sys.argv) > 1 else 'all'

ELF_base = r"S:\ELF_MAGIC\2020_03_07_CEFC_2020\model_C-Type\nonlinear_20000AT"

results_all = {}

# ============================================================
# 6x6x6
# ============================================================
if mesh_arg in ['6', 'all']:
    ELF_6 = os.path.join(ELF_base, "ELF_MMB8T_EIEM2_6x6x6")
    nodes, hex_elements = load_elf_geometry(ELF_6)
    n_elem = len(hex_elements)
    n_dof = n_elem * 6
    print(f"{'='*60}", flush=True)
    print(f"6x6x6 (Elements: {n_elem}, DOF: {n_dof})", flush=True)
    print(f"{'='*60}", flush=True)

    results_6 = []

    configs = [
        (0, 'LU', False, 100),
        (1, 'BiCGSTAB (dense)', False, 100),
        (2, 'HACApK + Block Jacobi', False, 100),
        (2, 'HACApK + Block Jacobi + Newton', True, 200),
    ]

    for method, name, newton, max_iter in configs:
        print(f"Running {name}...", flush=True)
        r = run_benchmark(nodes, hex_elements, bh_data, method, name,
                          use_newton=newton, max_iter=max_iter)
        print_result(r)
        results_6.append(r)

    results_all['6x6x6'] = results_6

# ============================================================
# 10x10x10
# ============================================================
if mesh_arg in ['10', 'all']:
    ELF_10 = os.path.join(ELF_base, "ELF_MMB8T_EIEM2_10x10x10")
    nodes, hex_elements = load_elf_geometry(ELF_10)
    n_elem = len(hex_elements)
    n_dof = n_elem * 6
    print(f"{'='*60}", flush=True)
    print(f"10x10x10 (Elements: {n_elem}, DOF: {n_dof})", flush=True)
    print(f"{'='*60}", flush=True)

    results_10 = []

    configs = [
        (0, 'LU', False, 100),
        (1, 'BiCGSTAB (dense)', False, 100),
        (2, 'HACApK + Block Jacobi', False, 100),
        (2, 'HACApK + Block Jacobi + Newton', True, 200),
    ]

    for method, name, newton, max_iter in configs:
        print(f"Running {name}...", flush=True)
        r = run_benchmark(nodes, hex_elements, bh_data, method, name,
                          use_newton=newton, max_iter=max_iter)
        print_result(r)
        results_10.append(r)

    results_all['10x10x10'] = results_10

# ============================================================
# 20x20x20
# ============================================================
if mesh_arg in ['20', 'all']:
    ELF_20 = os.path.join(ELF_base, "ELF_MMB8T_EIEM2_20x20x20")
    nodes, hex_elements = load_elf_geometry(ELF_20)
    n_elem = len(hex_elements)
    n_dof = n_elem * 6
    print(f"{'='*60}", flush=True)
    print(f"20x20x20 (Elements: {n_elem}, DOF: {n_dof})", flush=True)
    print(f"{'='*60}", flush=True)

    results_20 = []

    configs = [
        (2, 'HACApK + Block Jacobi', False, 100),
        (2, 'HACApK + Block Jacobi + Newton', True, 200),
    ]

    for method, name, newton, max_iter in configs:
        print(f"Running {name}...", flush=True)
        r = run_benchmark(nodes, hex_elements, bh_data, method, name,
                          use_newton=newton, max_iter=max_iter)
        print_result(r)
        results_20.append(r)

    results_all['20x20x20'] = results_20

# ============================================================
# Summary
# ============================================================
print(f"\n{'='*60}", flush=True)
print("SUMMARY", flush=True)
print(f"{'='*60}", flush=True)
for mesh, results in results_all.items():
    print(f"\n--- {mesh} ---", flush=True)
    for r in results:
        status = "OK" if r['converged'] else "FAIL"
        print(f"  {r['name']:40s} Bz={r['Bz_mT']:8.2f} mT  "
              f"Time={r['t_solve']:8.1f}s  NL={str(r['nonl_iter']):>4s}  "
              f"Lin={str(r['linear_iter']):>6s}  Mem={r['mem_mb']:>8s} MB  "
              f"Comp={r['compression']:>6s}  [{status}]", flush=True)

print("\nDone.", flush=True)
