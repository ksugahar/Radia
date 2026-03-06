#!/usr/bin/env python
"""
Common benchmark functions for Radia solver testing.

This module provides shared functionality for hexahedral and tetrahedral
benchmarks, supporting both linear and nonlinear materials.

Usage:
    from benchmark_common import run_benchmark, BH_DATA, H_EXT, MU_R
"""

import os
import sys
import time
import json
from typing import List, Dict, Any, Optional

# Add Radia to path
_src_path = os.path.join(os.path.dirname(__file__), '../../src')
sys.path.insert(0, _src_path)

import numpy as np
import radia as rad


try:
    import psutil
    HAS_PSUTIL = True
except ImportError:
    HAS_PSUTIL = False


# =============================================================================
# Physical Constants and Problem Parameters
# =============================================================================

MU_0 = 4 * np.pi * 1e-7  # Vacuum permeability [T/(A/m)]

# Geometry parameters
CUBE_SIZE = 1.0      # 1.0 m cube
CUBE_HALF = 0.5      # half size

# Field parameters
H_EXT = 200000.0     # External field (A/m) - sufficient for nonlinear saturation

# Linear material parameters
MU_R = 1000          # Relative permeability for linear material
CHI = MU_R - 1       # Magnetic susceptibility

# Analytical solution for linear material (demagnetizing factor N = 1/3 for cube)
N_DEMAG = 1.0 / 3.0
M_ANALYTICAL_Z = CHI * H_EXT / (1 + CHI * N_DEMAG)

# B-H curve data for nonlinear material - soft iron saturation curve
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


# =============================================================================
# Memory Measurement Functions
# =============================================================================

def get_current_memory_mb() -> Optional[float]:
    """Get current memory usage in MB (RSS)."""
    if not HAS_PSUTIL:
        return None
    process = psutil.Process(os.getpid())
    mem_info = process.memory_info()
    return mem_info.rss / (1024 * 1024)


def get_peak_memory_mb() -> Optional[float]:
    """Get peak memory usage in MB (Windows: peak_wset, Linux: max_rss)."""
    if not HAS_PSUTIL:
        return None
    process = psutil.Process(os.getpid())
    mem_info = process.memory_info()
    if hasattr(mem_info, 'peak_wset'):
        return mem_info.peak_wset / (1024 * 1024)
    else:
        try:
            import resource
            return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
        except (ImportError, AttributeError):
            return mem_info.rss / (1024 * 1024)


# Baseline memory at module load
_BASELINE_MEMORY_MB = get_current_memory_mb()


def get_solver_memory_mb() -> Optional[float]:
    """Get memory used by solver operations (current - baseline)."""
    current = get_current_memory_mb()
    if current is None or _BASELINE_MEMORY_MB is None:
        return None
    return max(0.0, current - _BASELINE_MEMORY_MB)


# =============================================================================
# Mesh Generation Functions
# =============================================================================

def generate_hex_mesh(n_div: int, size: float = 1.0) -> List[List[List[float]]]:
    """Generate hexahedral mesh vertices for a cube.

    Args:
        n_div: Number of divisions per edge
        size: Cube edge length

    Returns:
        List of vertex lists, each containing 8 vertices for one hexahedron
    """
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


# =============================================================================
# Main Benchmark Function
# =============================================================================

def run_benchmark(
    radia_obj,
    n_elements: int,
    solver_type: str,
    output_dir: str,
    element_type: str,
    mesh_description: str,
    t_mesh: float,
    is_linear: bool = False,
    nonl_tol: float = 0.001,
    bicg_tol: float = 1e-4,
    hmat_eps: float = 1e-4,
    hmat_leaf_size: int = 10,
    hmat_eta: float = 2.0,
    extra_data: Optional[Dict[str, Any]] = None
) -> Optional[Dict[str, Any]]:
    """Run benchmark with specified solver.

    Args:
        radia_obj: Radia object containing elements
        n_elements: Number of elements
        solver_type: 'lu', 'bicgstab', or 'hacapk'
        output_dir: Directory to save results
        element_type: 'hex' or 'tetra'
        mesh_description: Human-readable mesh description (e.g., 'N=10' or 'maxh=0.35m')
        t_mesh: Mesh generation time in seconds
        is_linear: True for linear material, False for nonlinear (BH curve)
        nonl_tol: Nonlinear iteration convergence tolerance
        bicg_tol: BiCGSTAB convergence tolerance
        hmat_eps: ACA tolerance for H-matrix
        hmat_leaf_size: H-matrix leaf size
        hmat_eta: H-matrix admissibility parameter
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

    # Apply material
    if is_linear:
        mat = rad.MatLin(MU_R)
    else:
        mat = rad.MatSatIsoTab(BH_DATA)
    rad.MatApl(radia_obj, mat)

    # External field (ObjBckg callback returns B in Tesla)
    B_ext = MU_0 * H_EXT
    ext = rad.ObjBckg(lambda p: [0, 0, B_ext])
    grp = rad.ObjCnt([radia_obj, ext])

    # Configure H-matrix if using HACApK
    hmatrix_enabled = False
    if solver_method == 2:
        try:
            rad.SolverConfig(hacapk_eps=hmat_eps, hacapk_leaf=hmat_leaf_size, hacapk_eta=hmat_eta)
            hmatrix_enabled = True
            print('H-matrix: Enabled (eps=%.0e, leaf_size=%d, eta=%.1f)' % (
                hmat_eps, hmat_leaf_size, hmat_eta))
        except AttributeError:
            print('H-matrix: Not available (API not found)')

    # Set solver tolerances
    MAX_ITER = 100
    rad.SolverConfig(bicgstab_tol=bicg_tol, relax_param=0.0)

    # Solve
    print('Solving...')
    t_solve_start = time.time()
    result = rad.Solve(grp, nonl_tol, MAX_ITER, solver_method)
    t_solve = time.time() - t_solve_start

    # Memory measurement
    solver_memory_mb = get_solver_memory_mb()
    peak_memory_mb = get_peak_memory_mb()

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

    # Calculate error for linear material
    error_percent = None
    if is_linear:
        error_percent = abs(M_avg_z - M_ANALYTICAL_Z) / M_ANALYTICAL_Z * 100

    # Get H-matrix info if applicable
    hmat_info = None
    if hmatrix_enabled:
        try:
            config = rad.GetSolverConfig()
            hmat_info = config.get('hacapk_stats')
        except AttributeError:
            pass

    # Print results
    print('Mesh time:       %.4f s' % t_mesh)
    print('Solve time:      %.3f s' % t_solve)
    print('Nonl iter:       %d' % n_iter)
    print('Linear iter:     %d' % n_linear_iter)
    print('Converged:       %s' % ('Yes' if converged else 'No'))
    print('M_avg_z:         %.0f A/m' % M_avg_z)
    if is_linear:
        print('Analytical:      %.0f A/m' % M_ANALYTICAL_Z)
        print('Error:           %.2f%%' % error_percent)
    if peak_memory_mb is not None:
        print('Peak memory:     %.1f MB' % peak_memory_mb)
    if solver_memory_mb is not None:
        print('Solver memory:   %.1f MB (excluding baseline)' % solver_memory_mb)

    # Print H-matrix statistics if available
    if hmat_info:
        print('--- H-matrix ---')
        print('  lowrank:       %d' % hmat_info.get('n_lowrank', 0))
        print('  dense:         %d' % hmat_info.get('n_dense', 0))
        print('  max_rank:      %d' % hmat_info.get('max_rank', 0))
        hmat_mem = hmat_info.get('memory_mb', 0.0)
        dense_mem = hmat_info.get('dense_memory_mb', 0.0)
        compression = hmat_info.get('compression', 0.0)
        print('  H-mat memory:  %.2f MB' % hmat_mem)
        print('  Dense memory:  %.2f MB' % dense_mem)
        print('  Compression:   %.1f%%' % (compression * 100))
    print()

    # Build result dictionary
    ndof = n_elements * (6 if element_type == 'hex' else 3)

    result_data = {
        'element_type': element_type,
        'mesh_description': mesh_description,
        'n_elements': n_elements,
        'ndof': ndof,
        'H_ext': H_EXT,
        'material_type': 'linear' if is_linear else 'nonlinear',
        # Linear material parameters
        'mu_r': MU_R if is_linear else None,
        'chi': CHI if is_linear else None,
        # Solver parameters
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

    if is_linear:
        result_data['M_analytical_z'] = M_ANALYTICAL_Z
        result_data['error_percent'] = error_percent

    if peak_memory_mb is not None:
        result_data['peak_memory_mb'] = peak_memory_mb
    if solver_memory_mb is not None:
        result_data['solver_memory_mb'] = solver_memory_mb

    # Add detailed timing
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

    # Add H-matrix stats if available
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
            'hmat_eps': hmat_eps,
            'leaf_size': hmat_leaf_size,
            'eta': hmat_eta,
        }

    # Add extra data if provided
    if extra_data:
        result_data.update(extra_data)

    # Save result
    os.makedirs(output_dir, exist_ok=True)

    # Generate filename
    if element_type == 'hex':
        n_div = extra_data.get('n_div', 0) if extra_data else 0
        filename = 'hex_N%d_results.json' % n_div
    else:
        maxh = extra_data.get('maxh', 0.0) if extra_data else 0.0
        maxh_str = ('%.2fm' % maxh).replace('.', '_')
        filename = 'tetra_maxh%s_results.json' % maxh_str

    filepath = os.path.join(output_dir, filename)

    with open(filepath, 'w') as f:
        json.dump(result_data, f, indent=2)
    print('Saved: %s' % filepath)

    return result_data


def print_summary(results: List[Dict[str, Any]], solver_name: str, output_dir: str):
    """Print benchmark summary table."""
    if not results:
        return

    print('\n%s Solver (%s/):\n' % (solver_name, output_dir))

    desc_width = max(len(r['mesh_description']) for r in results)
    desc_width = max(desc_width, 10)

    # Check if linear (has error_percent)
    is_linear = 'error_percent' in results[0]

    if is_linear:
        print('%-*s %10s %10s %10s %12s %10s' % (
            desc_width, 'Mesh', 'Elements', 'Time (s)', 'Iterations', 'M_avg_z', 'Error %'))
        print('-' * (desc_width + 65))
        for r in results:
            print('%-*s %10d %10.3f %10d %12.0f %10.2f' % (
                desc_width,
                r['mesh_description'],
                r['n_elements'],
                r['t_solve'],
                r['nonl_iterations'],
                r['M_avg_z'],
                r.get('error_percent', 0)
            ))
    else:
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
