#!/usr/bin/env python
"""
Common benchmark functions for Radia nonlinear solver testing.

This module provides shared functionality for hexahedral and tetrahedral
benchmarks, matching the ELF output format for fair comparison.

Usage:
    from benchmark_common import run_nonlinear_benchmark, BH_DATA, H_EXT
"""

import os
import sys
import time
import json
from typing import List, Dict, Any, Optional

# Add Radia to path
_src_path = os.path.join(os.path.dirname(__file__), '../../../src/radia')
sys.path.insert(0, _src_path)

import numpy as np
import radia as rad

try:
    import psutil
    HAS_PSUTIL = True
except ImportError:
    HAS_PSUTIL = False


# Problem parameters - shared between hex and tetra benchmarks
CUBE_SIZE = 1.0      # 1.0 m cube
CUBE_HALF = 0.5      # half size
H_EXT = 200000.0     # External field (A/m) - enough for nonlinear behavior

# Physical constants
MU_0 = 4 * np.pi * 1e-7  # T/(A/m)

# B-H curve data - soft iron saturation curve (matches ELF exactly)
BH_DATA = [
    [0.0, 0.0], [100.0, 0.1], [200.0, 0.3], [500.0, 0.8], [1000.0, 1.2],
    [2000.0, 1.5], [5000.0, 1.7], [10000.0, 1.8], [50000.0, 2.0], [100000.0, 2.1],
]

# Hexahedral face topology (1-indexed for Radia)
HEX_FACES = [
    [1, 4, 3, 2],  # Bottom face (z=0)
    [5, 6, 7, 8],  # Top face (z=1)
    [1, 2, 6, 5],  # Front face (y=0)
    [3, 4, 8, 7],  # Back face (y=1)
    [1, 5, 8, 4],  # Left face (x=0)
    [2, 3, 7, 6]   # Right face (x=1)
]


def get_current_memory_mb() -> Optional[float]:
    """Get current memory usage in MB."""
    if not HAS_PSUTIL:
        return None
    process = psutil.Process(os.getpid())
    mem_info = process.memory_info()
    return mem_info.rss / (1024 * 1024)  # MB


def get_peak_memory_mb() -> Optional[float]:
    """Get peak memory usage in MB (Windows: peak_wset)"""
    if not HAS_PSUTIL:
        return None
    process = psutil.Process(os.getpid())
    mem_info = process.memory_info()
    if hasattr(mem_info, 'peak_wset'):
        return mem_info.peak_wset / (1024 * 1024)  # MB
    else:
        return mem_info.rss / (1024 * 1024)  # MB (fallback)


def run_nonlinear_benchmark(
    radia_obj,
    n_elements: int,
    solver_type: str,
    output_dir: str,
    element_type: str,
    mesh_description: str,
    t_mesh: float,
    nonl_tol: float = 0.001,
    bicg_tol: float = 1e-4,
    hmat_eps: float = 1e-4,
    hmat_leaf_size: int = 10,
    hmat_eta: float = 2.0,
    extra_data: Optional[Dict[str, Any]] = None
) -> Optional[Dict[str, Any]]:
    """Run nonlinear benchmark with specified solver.

    Args:
        radia_obj: Radia object containing elements
        n_elements: Number of elements
        solver_type: 'lu', 'bicgstab', or 'hacapk'
        output_dir: Directory to save results
        element_type: 'hex' or 'tetra'
        mesh_description: Human-readable mesh description (e.g., 'N=10' or 'maxh=0.35m')
        t_mesh: Mesh generation time in seconds
        nonl_tol: Nonlinear iteration convergence tolerance (ELF default: 0.001)
        bicg_tol: BiCGSTAB convergence tolerance (ELF default: 1e-4)
        hmat_eps: ACA tolerance for H-matrix (ELF default: 1e-4)
        hmat_leaf_size: H-matrix leaf size (ELF default: 10)
        hmat_eta: H-matrix admissibility parameter (ELF default: 2.0)
        extra_data: Additional data to include in result JSON

    Returns:
        Result dictionary or None if failed
    """
    if n_elements == 0:
        print('ERROR: No elements provided!')
        return None

    # Map solver type to method number
    solver_method_map = {'lu': 0, 'bicgstab': 1, 'hacapk': 2}
    solver_method = solver_method_map.get(solver_type, 0)

    # Setup material (BH curve)
    mat = rad.MatSatIsoTab(BH_DATA)
    rad.MatApl(radia_obj, mat)

    # External field H_z
    B_ext = MU_0 * H_EXT
    ext = rad.ObjBckg([0, 0, B_ext])
    grp = rad.ObjCnt([radia_obj, ext])

    # Configure H-matrix if using HACApK
    hmatrix_enabled = False
    if solver_method == 2:
        try:
            rad.SetHACApKParams(hmat_eps, hmat_leaf_size, hmat_eta)
            hmatrix_enabled = True
            print('H-matrix: Enabled (eps=%.0e, leaf_size=%d, eta=%.1f)' % (hmat_eps, hmat_leaf_size, hmat_eta))
        except AttributeError:
            print('H-matrix: Not available (API not found)')

    # Set solver tolerances
    MAX_ITER = 100
    rad.SetBiCGSTABTol(bicg_tol)
    rad.SetRelaxParam(0.0)

    # Measure memory before solve
    mem_before = get_current_memory_mb()

    # Solve
    print('Solving...')
    t_solve_start = time.time()
    result = rad.Solve(grp, nonl_tol, MAX_ITER, solver_method)
    t_solve = time.time() - t_solve_start

    # Measure memory after solve
    mem_after = get_current_memory_mb()
    if mem_before is not None and mem_after is not None:
        solver_memory_mb = mem_after
    else:
        solver_memory_mb = get_peak_memory_mb()

    # Get solve statistics
    stats = rad.GetSolveStats()
    n_iter = stats.get('nonl_iterations', 0)
    n_linear_iter = stats.get('linear_iterations', 0)
    converged = n_iter < MAX_ITER
    residual = result[0] if result[0] else 0.0

    # Get average magnetization
    all_M = rad.ObjM(radia_obj)
    M_total_z = sum(m[1][2] for m in all_M)
    M_avg_z = M_total_z / n_elements

    # Get H-matrix info if applicable
    hmat_info = None
    if hmatrix_enabled:
        try:
            hmat_info = rad.GetHACApKStats()
        except AttributeError:
            pass

    # Print results
    print('Mesh time:       %.4f s' % t_mesh)
    print('Solve time:      %.3f s' % t_solve)
    print('Nonl iter:       %d' % n_iter)
    print('Linear iter:     %d' % n_linear_iter)
    print('Converged:       %s' % ('Yes' if converged else 'No'))
    print('M_avg_z:         %.0f A/m' % M_avg_z)
    if solver_memory_mb is not None:
        print('Memory (RSS):    %.1f MB' % solver_memory_mb)
    print()

    # Build result dictionary (matching ELF format exactly)
    # Determine DOF based on element type
    if element_type == 'hex':
        ndof = n_elements * 6  # 6 DOF per hexahedron
    else:
        ndof = n_elements * 3  # 3 DOF per tetrahedron

    result_data = {
        'element_type': element_type,
        'mesh_description': mesh_description,
        'n_elements': n_elements,
        'ndof': ndof,
        'H_ext': H_EXT,
        # Solver parameters (ELF-compatible)
        'nonl_tol': nonl_tol,
        'bicg_tol': bicg_tol if solver_type in ['bicgstab', 'hacapk'] else None,
        'hmat_eps': hmat_eps if solver_type == 'hacapk' else None,
        'hmat_leaf_size': hmat_leaf_size if solver_type == 'hacapk' else None,
        'hmat_eta': hmat_eta if solver_type == 'hacapk' else None,
        # Timing
        't_mesh': t_mesh,
        't_solve': t_solve,
        # Solver identification
        'solver_method': solver_method,
        'solver_name': solver_type,
        # Convergence
        'converged': converged,
        'residual': residual,
        'nonl_iterations': n_iter,
        'linear_iterations': n_linear_iter,
        # Results
        'M_avg_z': M_avg_z,
    }

    if solver_memory_mb is not None:
        result_data['peak_memory_mb'] = solver_memory_mb

    # Add detailed timing (matching ELF format)
    result_data['timing'] = {
        't_matrix_build': stats.get('t_matrix_build', 0.0),
        't_lu_decomp': stats.get('t_lu_decomp', 0.0),
        't_hmatrix_build': stats.get('t_hmatrix_build', 0.0),
        't_hmatrix_cluster': stats.get('t_hmatrix_cluster', 0.0),
        't_hmatrix_frame': stats.get('t_hmatrix_frame', 0.0),
        't_hmatrix_fill': stats.get('t_hmatrix_fill', 0.0),
        't_linear_solve': stats.get('t_linear_solve', 0.0),
        't_total': t_solve
    }

    # Add H-matrix stats if available (matching ELF format)
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
            # Solver parameters used
            'hmat_eps': hmat_eps,
            'leaf_size': hmat_leaf_size,
            'eta': hmat_eta,
        }

    # Add extra data if provided
    if extra_data:
        result_data.update(extra_data)

    # Save result
    os.makedirs(output_dir, exist_ok=True)

    # Generate filename matching ELF format
    if element_type == 'hex':
        # Extract N from mesh_description like 'N=10'
        n_div = extra_data.get('n_div', 0) if extra_data else 0
        filename = 'hex_N%d_results.json' % n_div
    else:
        # Extract maxh from mesh_description like 'maxh=0.35m'
        maxh = extra_data.get('maxh', 0.0) if extra_data else 0.0
        # Convert to format like '0_35m'
        maxh_str = ('%.2fm' % maxh).replace('.', '_')
        filename = 'tetra_maxh%s_results.json' % maxh_str

    filepath = os.path.join(output_dir, filename)

    with open(filepath, 'w') as f:
        json.dump(result_data, f, indent=2)
    print('Saved: %s' % filepath)

    return result_data


def print_summary(results: List[Dict[str, Any]], solver_name: str, output_dir: str):
    """Print benchmark summary table.

    Args:
        results: List of result dictionaries
        solver_name: Human-readable solver name
        output_dir: Output directory name for display
    """
    if not results:
        return

    print('\n%s Solver (%s/):\n' % (solver_name, output_dir))

    # Determine column width based on mesh_description
    desc_width = max(len(r['mesh_description']) for r in results)
    desc_width = max(desc_width, 10)

    print('%-*s %10s %10s %10s %12s %12s %10s' % (
        desc_width, 'Mesh', 'Elements', 'Time (s)', 'Nonl Iter', 'Linear Iter', 'M_avg_z', 'Conv'))
    print('-' * (desc_width + 75))

    for r in results:
        print('%-*s %10d %10.3f %10d %12d %12.0f %10s' % (
            desc_width,
            r['mesh_description'],
            r['n_elements'],
            r['t_solve'],
            r['nonl_iterations'],
            r.get('linear_iterations', 0),
            r['M_avg_z'],
            'Yes' if r['converged'] else 'No'
        ))
