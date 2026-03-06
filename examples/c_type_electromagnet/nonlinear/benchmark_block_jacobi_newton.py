#!/usr/bin/env python
"""Benchmark: C-type electromagnet with LU, BiCGSTAB, and HACApK solvers.

Results are saved per solver per mesh as JSON files:

  quarter/{solver}/{n_dof}DOF.json   (with image='+x-z')
  full/{solver}/{n_dof}DOF.json      (without image symmetry)

Solver directories: lu/, bicgstab/, hacapk/

Timing and H-matrix statistics follow the convention used in
cube_uniform_field/benchmark_common.py.

IMPORTANT: Run sequentially (one job at a time) for reliable timing benchmarks.

Usage: python -u benchmark_block_jacobi_newton.py [6|10|20|all] [quarter|full]
"""
import sys
import os
import time
import json
from datetime import datetime
import numpy as np

work_dir = os.path.dirname(os.path.abspath(__file__))
repo_root = os.path.join(work_dir, '..', '..', '..')
sys.path.insert(0, os.path.join(repo_root, 'src'))
sys.path.insert(0, work_dir)

import radia as rad
from coil_model import create_racetrack_coil

try:
    import psutil
    HAS_PSUTIL = True
except ImportError:
    HAS_PSUTIL = False

scale = 0.001  # mm to m

# HACApK parameters
HMAT_EPS = 1e-4
HMAT_LEAF_SIZE = 10
HMAT_ETA = 2.0

# Nonlinear solver parameters
NONL_TOL = 0.001
CURRENT_AT = 20000.0

# ELF reference values: Bz [mT] at origin (with MIMA, 20000 AT)
ELF_REFERENCE = {
    '6x6x6':    -960.78,
    '10x10x10': None,
    '20x20x20': None,
}


# =============================================================================
# Memory measurement (same approach as cube_uniform_field/benchmark_common.py)
# =============================================================================

def get_current_memory_mb():
    if not HAS_PSUTIL:
        return None
    return psutil.Process(os.getpid()).memory_info().rss / (1024 * 1024)


def get_peak_memory_mb():
    if not HAS_PSUTIL:
        return None
    mem_info = psutil.Process(os.getpid()).memory_info()
    if hasattr(mem_info, 'peak_wset'):
        return mem_info.peak_wset / (1024 * 1024)
    return mem_info.rss / (1024 * 1024)


# =============================================================================
# Geometry and material loading
# =============================================================================

def load_elf_geometry(path):
    """Load ELF yoke geometry from .meg file (MMB8T elements only)."""
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
    """Load B-H curve from text file."""
    data = np.loadtxt(filepath, comments='#')
    return data.tolist()


# =============================================================================
# Result saving
# =============================================================================

SOLVER_DIRS = {0: 'lu', 1: 'bicgstab', 2: 'hacapk'}


def save_result(base_dir, solver_method, n_dof, result_data):
    """Save result to {base_dir}/{solver}/{n_dof}DOF.json"""
    solver_dir = SOLVER_DIRS.get(solver_method, 'unknown')
    out_dir = os.path.join(base_dir, solver_dir)
    os.makedirs(out_dir, exist_ok=True)
    filepath = os.path.join(out_dir, f"{n_dof}DOF.json")
    with open(filepath, 'w', encoding='utf-8') as f:
        json.dump(result_data, f, indent=2, ensure_ascii=False, default=str)
    print(f"  -> Saved: {filepath}", flush=True)
    return filepath


# =============================================================================
# Benchmark runner
# =============================================================================

def run_benchmark(nodes, hex_elements, bh_data, solver_method, solver_name,
                  mesh_name, use_ima=True, max_iter=100):
    """Run a single solver benchmark and return result dict.

    Parameters
    ----------
    solver_method : int
        0=LU, 1=BiCGSTAB, 2=HACApK
    use_ima : bool
        True for quarter model with image='+x-z'
    """
    rad.UtiDelAll()
    mat = rad.MatSatIsoTab(bh_data)
    all_objects = []
    for elem_id, node_ids in hex_elements:
        verts = [[nodes[nid][0]*scale, nodes[nid][1]*scale, nodes[nid][2]*scale]
                 for nid in node_ids]
        obj = rad.ObjHexahedron(verts, [0, 0, 0])
        rad.MatApl(obj, mat)
        all_objects.append(obj)
    yoke = rad.ObjCnt(all_objects)
    coil = create_racetrack_coil(CURRENT_AT)
    model = rad.ObjCnt([yoke, coil])

    n_elem = len(hex_elements)
    n_dof = n_elem * 6
    dense_mem_mb = n_dof * n_dof * 8 / 1024 / 1024

    # Configure HACApK
    hmatrix_enabled = (solver_method == 2)
    if hmatrix_enabled:
        rad.SolverConfig(hacapk_eps=HMAT_EPS, hacapk_leaf=HMAT_LEAF_SIZE, hacapk_eta=HMAT_ETA)

    # HACApK uses Newton (hybrid Picard+Newton)
    use_newton = hmatrix_enabled
    rad.SolverConfig(newton_method=use_newton)

    # Measure memory before solve
    mem_before = get_current_memory_mb()

    t0 = time.time()
    converged = True
    try:
        if use_ima:
            result = rad.Solve(model, NONL_TOL, max_iter, solver_method,
                               image='+x-z')
        else:
            result = rad.Solve(model, NONL_TOL, max_iter, solver_method)
    except RuntimeError as e:
        converged = False
        result = [0, 0]
        print(f"  WARNING: {e}", flush=True)
    t_solve = time.time() - t0

    # Measure memory after solve
    mem_after = get_current_memory_mb()
    peak_mem = get_peak_memory_mb()
    solver_memory_mb = None
    if mem_before is not None and mem_after is not None:
        solver_memory_mb = mem_after - mem_before

    # Field at origin
    B = rad.Fld(model, 'b', [0, 0, 0])

    # Solver statistics (iterations + timing breakdown)
    stats = rad.GetSolveStats() or {}
    n_iter = stats.get('nonl_iterations', None)
    n_linear_iter = stats.get('linear_iterations', None)
    residual = result[1] if len(result) > 1 else None

    # H-matrix statistics
    hmat_info = None
    if hmatrix_enabled:
        try:
            config = rad.GetSolverConfig()
            hmat_info = config.get('hacapk_stats')
        except (AttributeError, Exception):
            pass

    # Reset Newton method
    rad.SolverConfig(newton_method=False)

    # ELF reference
    elf_ref = ELF_REFERENCE.get(mesh_name)

    # Build result dict (matching cube_uniform_field convention)
    result_data = {
        'timestamp': datetime.now().isoformat(),
        'mesh_description': mesh_name,
        'n_elements': n_elem,
        'n_dof': n_dof,
        'current_at': CURRENT_AT,
        'material_type': 'nonlinear',
        'image': '+x-z' if use_ima else None,
        # Solver parameters
        'solver_method': solver_method,
        'solver_name': solver_name,
        'use_newton': use_newton,
        'nonl_tol': NONL_TOL,
        'max_iter': max_iter,
        # HACApK parameters (if applicable)
        'hmat_eps': HMAT_EPS if hmatrix_enabled else None,
        'hmat_leaf_size': HMAT_LEAF_SIZE if hmatrix_enabled else None,
        'hmat_eta': HMAT_ETA if hmatrix_enabled else None,
        # Convergence
        'converged': converged,
        'residual': residual,
        'nonl_iterations': n_iter,
        'linear_iterations': n_linear_iter,
        # Results
        'Bz_mT': B[2] * 1000,
        'elf_reference_Bz_mT': elf_ref,
        # Timing (wall clock)
        't_solve': t_solve,
    }

    # Timing breakdown (from GetSolveStats)
    result_data['timing'] = {
        't_matrix_build': stats.get('t_matrix_build', 0.0),
        't_lu_decomp': stats.get('t_lu_decomp', 0.0),
        't_hmatrix_build': stats.get('t_hmatrix_build', 0.0),
        't_hmatrix_cluster': stats.get('t_hmatrix_cluster', 0.0),
        't_hmatrix_frame': stats.get('t_hmatrix_frame', 0.0),
        't_hmatrix_fill': stats.get('t_hmatrix_fill', 0.0),
        't_linear_solve': stats.get('t_linear_solve', 0.0),
        't_total': t_solve,
        'num_threads': stats.get('num_threads', 1),
        'taskmanager_enabled': stats.get('taskmanager_enabled', False),
    }

    # Memory
    result_data['dense_memory_mb'] = dense_mem_mb
    if peak_mem is not None:
        result_data['peak_memory_mb'] = peak_mem
    if solver_memory_mb is not None:
        result_data['solver_memory_mb'] = solver_memory_mb

    # H-matrix detailed statistics
    if hmat_info:
        result_data['hmatrix'] = {
            'n_lowrank': hmat_info.get('n_lowrank', 0),
            'n_dense': hmat_info.get('n_dense', 0),
            'max_rank': hmat_info.get('max_rank', 0),
            'compression_ratio': hmat_info.get('compression', 0.0),
            'build_time': hmat_info.get('build_time', 0.0),
            'memory_mb': hmat_info.get('memory_mb', 0.0),
            'dense_memory_mb': hmat_info.get('dense_memory_mb', 0.0),
            'nlf': hmat_info.get('n_leaves', 0),
            'hmat_eps': HMAT_EPS,
            'leaf_size': HMAT_LEAF_SIZE,
            'eta': HMAT_ETA,
        }

    return result_data


def print_result(r):
    """Print formatted result to stdout."""
    status = "CONVERGED" if r['converged'] else "NOT CONVERGED"
    print(f"  {r['solver_name']}:", flush=True)
    print(f"    Bz         = {r['Bz_mT']:.2f} mT", flush=True)
    if r.get('elf_reference_Bz_mT') is not None:
        err = (r['Bz_mT'] - r['elf_reference_Bz_mT']) / abs(r['elf_reference_Bz_mT']) * 100
        print(f"    ELF ref    = {r['elf_reference_Bz_mT']:.2f} mT (err: {err:+.2f}%)", flush=True)
    print(f"    Time       = {r['t_solve']:.1f} s", flush=True)
    print(f"    NL iter    = {r['nonl_iterations']}", flush=True)
    print(f"    Lin iter   = {r['linear_iterations']}", flush=True)
    print(f"    Dense mem  = {r['dense_memory_mb']:.1f} MB", flush=True)
    # Timing breakdown
    t = r.get('timing', {})
    if t.get('t_matrix_build', 0) > 0:
        print(f"    t_mat_build= {t['t_matrix_build']:.3f} s", flush=True)
    if t.get('t_lu_decomp', 0) > 0:
        print(f"    t_lu_decomp= {t['t_lu_decomp']:.3f} s", flush=True)
    if t.get('t_hmatrix_build', 0) > 0:
        print(f"    t_hmat_bld = {t['t_hmatrix_build']:.3f} s", flush=True)
    if t.get('t_linear_solve', 0) > 0:
        print(f"    t_lin_solve= {t['t_linear_solve']:.3f} s", flush=True)
    # H-matrix stats
    hm = r.get('hmatrix')
    if hm:
        print(f"    H-mat mem  = {hm['memory_mb']:.1f} MB", flush=True)
        print(f"    Compress   = {hm['compression_ratio']*100:.1f}%", flush=True)
        print(f"    Max rank   = {hm['max_rank']}", flush=True)
    # Memory
    if r.get('peak_memory_mb') is not None:
        print(f"    Peak mem   = {r['peak_memory_mb']:.0f} MB", flush=True)
    if r.get('solver_memory_mb') is not None:
        print(f"    Solver mem = {r['solver_memory_mb']:.1f} MB", flush=True)
    print(f"    Threads    = {t.get('num_threads', '?')} (TM: {t.get('taskmanager_enabled', '?')})", flush=True)
    print(f"    Status     = {status}", flush=True)
    print(flush=True)


# =============================================================================
# Main
# =============================================================================

if __name__ == '__main__':
    bh_file = os.path.join(work_dir, "BH.txt")
    bh_data = load_bh_curve(bh_file)

    mesh_arg = sys.argv[1] if len(sys.argv) > 1 else 'all'
    mode_arg = sys.argv[2] if len(sys.argv) > 2 else 'quarter'
    use_ima = (mode_arg == 'quarter')

    base_dir = os.path.join(work_dir, mode_arg)

    ELF_base = r"S:\ELF_MAGIC\2020_03_07_CEFC_2020\model_C-Type\nonlinear_20000AT"

    # Solver configurations: (method, name, max_iter)
    CONFIGS_ALL = [
        (0, 'LU',       100),
        (1, 'BiCGSTAB', 100),
        (2, 'HACApK',   200),
    ]

    # For 20x20x20: dense solvers will fail (204 GB > 128 GB RAM)
    CONFIGS_HACAPK_ONLY = [
        (2, 'HACApK', 1000),
    ]

    MESHES = {
        '6':  ('ELF_MMB8T_EIEM2_6x6x6',   '6x6x6',   CONFIGS_ALL),
        '10': ('ELF_MMB8T_EIEM2_10x10x10', '10x10x10', CONFIGS_ALL),
        '20': ('ELF_MMB8T_EIEM2_20x20x20', '20x20x20', CONFIGS_HACAPK_ONLY),
    }

    print(f"Mode: {mode_arg} (IMA: {'+x-z' if use_ima else 'none'})", flush=True)
    print(f"Output: {base_dir}", flush=True)
    print(flush=True)

    results_all = {}

    for mesh_key, (elf_dir, mesh_name, configs) in MESHES.items():
        if mesh_arg not in [mesh_key, 'all']:
            continue

        elf_path = os.path.join(ELF_base, elf_dir)
        if not os.path.isdir(elf_path):
            print(f"Skipping {mesh_name}: {elf_path} not found", flush=True)
            continue

        nodes, hex_elements = load_elf_geometry(elf_path)
        n_elem = len(hex_elements)
        n_dof = n_elem * 6

        print(f"{'='*60}", flush=True)
        print(f"{mesh_name} (Elements: {n_elem}, DOF: {n_dof})", flush=True)
        print(f"{'='*60}", flush=True)

        results = []
        for method, name, max_iter in configs:
            print(f"Running {name}...", flush=True)
            r = run_benchmark(nodes, hex_elements, bh_data, method, name,
                              mesh_name, use_ima=use_ima, max_iter=max_iter)
            print_result(r)
            save_result(base_dir, method, n_dof, r)
            results.append(r)

        results_all[mesh_name] = results

    # Summary
    print(f"\n{'='*60}", flush=True)
    print(f"SUMMARY ({mode_arg})", flush=True)
    print(f"{'='*60}", flush=True)

    summary = []
    for mesh, results in results_all.items():
        print(f"\n--- {mesh} ---", flush=True)
        for r in results:
            status = "OK" if r['converged'] else "FAIL"
            hm = r.get('hmatrix')
            comp_str = f"{hm['compression_ratio']*100:.1f}%" if hm else "---"
            nl = str(r['nonl_iterations']) if r['nonl_iterations'] is not None else "N/A"
            lin = str(r['linear_iterations']) if r['linear_iterations'] is not None else "N/A"
            print(f"  {r['solver_name']:12s} Bz={r['Bz_mT']:8.2f} mT  "
                  f"Time={r['t_solve']:8.1f}s  NL={nl:>4s}  Lin={lin:>6s}  "
                  f"Dense={r['dense_memory_mb']:8.1f} MB  Comp={comp_str:>6s}  [{status}]",
                  flush=True)
            summary.append(r)

    # Save combined summary
    os.makedirs(base_dir, exist_ok=True)
    summary_path = os.path.join(base_dir, "summary.json")
    with open(summary_path, 'w', encoding='utf-8') as f:
        json.dump(summary, f, indent=2, ensure_ascii=False, default=str)
    print(f"\nSummary saved: {summary_path}", flush=True)

    print("\nDone.", flush=True)
