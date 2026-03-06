#!/usr/bin/env python
"""
Aggregate hexahedron benchmark JSON results and plot solver scaling laws.

Fits power-law scaling t = a * N^alpha using the last FIT_LAST_N data points,
and plots time and memory vs DOF.

Usage:
    python plot_scaling.py
"""

import os
import json
import glob
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Global font size multiplier
FONT_SCALE = 1.5
plt.rcParams.update({
    'font.family': 'serif',
    'font.serif': ['Times New Roman'],
    'mathtext.fontset': 'stix',
    'font.size': 10 * FONT_SCALE,
    'axes.labelsize': 11 * FONT_SCALE,
    'xtick.labelsize': 9 * FONT_SCALE,
    'ytick.labelsize': 9 * FONT_SCALE,
    'legend.fontsize': 9 * FONT_SCALE,
})

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# Number of trailing data points used for power-law fitting
FIT_LAST_N = 9

# Solver display settings
SOLVER_STYLES = {
    'lu':       {'color': '#1f77b4', 'marker': 'o', 'label': 'LU'},
    'bicgstab': {'color': '#ff7f0e', 'marker': 's', 'label': 'BiCGSTAB'},
    'hacapk':   {'color': '#2ca02c', 'marker': '^', 'label': 'HACApK'},
}


def load_results():
    """Load all hexahedron JSON result files organized by solver."""
    data = {}
    base = os.path.join(SCRIPT_DIR, 'nonlinear')
    if not os.path.isdir(base):
        return data

    for solver in ['lu', 'bicgstab', 'hacapk']:
        solver_dir = os.path.join(base, solver)
        if not os.path.isdir(solver_dir):
            continue

        results = []
        for jf in sorted(glob.glob(os.path.join(solver_dir, '*_results.json'))):
            with open(jf, 'r') as f:
                results.append(json.load(f))

        if results:
            results.sort(key=lambda r: r['ndof'])
            data[solver] = results

    return data


def fit_power_law(ndofs, times):
    """Fit t = a * N^alpha using least-squares in log space.

    Returns (alpha, a, r_squared).
    """
    log_n = np.log(ndofs)
    log_t = np.log(times)
    coeffs = np.polyfit(log_n, log_t, 1)
    alpha = coeffs[0]
    a = np.exp(coeffs[1])

    log_t_pred = np.polyval(coeffs, log_n)
    ss_res = np.sum((log_t - log_t_pred) ** 2)
    ss_tot = np.sum((log_t - np.mean(log_t)) ** 2)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0

    return alpha, a, r2


def plot_scaling(data, quantity='t_solve', ylabel='Solve time [s]'):
    """Plot scaling for a given quantity."""
    fig, ax = plt.subplots(figsize=(8, 6))

    # First pass: plot data points
    fit_info = []
    for solver, results in data.items():
        style = SOLVER_STYLES.get(solver, {'color': 'gray', 'marker': 'x', 'label': solver})

        ndofs = np.array([r['ndof'] for r in results], dtype=float)
        values = np.array([r[quantity] for r in results if r.get(quantity) is not None], dtype=float)

        if len(values) < 2:
            continue

        ndofs = ndofs[:len(values)]

        ax.loglog(ndofs, values, style['marker'],
                  color=style['color'], markersize=7, alpha=0.7)

        n_fit_pts = min(FIT_LAST_N, len(ndofs))
        alpha, a, r2 = fit_power_law(ndofs[-n_fit_pts:], values[-n_fit_pts:])
        fit_info.append((ndofs[-n_fit_pts], alpha, a, style))

    # Get auto-scaled x range, then plot fit lines to the right edge
    ax.set_xlim(left=1e3)
    x_right = ax.get_xlim()[1]
    for n_start, alpha, a, style in fit_info:
        n_fit = np.logspace(np.log10(n_start), np.log10(x_right), 50)
        t_fit = a * n_fit ** alpha
        ax.loglog(n_fit, t_fit, '--', color=style['color'], linewidth=2,
                  label='%s: $O(N^{%.2f})$' % (style['label'], alpha))

    ax.set_xlabel('DOF')
    ax.set_ylabel(ylabel)
    ax.legend()
    ax.grid(True, which='both', alpha=0.3)
    fig.tight_layout()
    return fig


def plot_memory_scaling(data):
    """Plot memory scaling."""
    fig, ax = plt.subplots(figsize=(8, 6))

    # First pass: plot data points
    fit_info = []
    for solver, results in data.items():
        style = SOLVER_STYLES.get(solver, {'color': 'gray', 'marker': 'x', 'label': solver})

        ndofs = []
        mems = []
        for r in results:
            mem = r.get('peak_memory_mb')
            if mem is not None and mem > 0:
                ndofs.append(r['ndof'])
                mems.append(mem)

        if len(ndofs) < 2:
            continue

        ndofs = np.array(ndofs, dtype=float)
        mems = np.array(mems, dtype=float)

        ax.loglog(ndofs, mems, style['marker'],
                  color=style['color'], markersize=7, alpha=0.7)

        n_fit_pts = min(FIT_LAST_N, len(ndofs))
        alpha, a, r2 = fit_power_law(ndofs[-n_fit_pts:], mems[-n_fit_pts:])
        fit_info.append((ndofs[-n_fit_pts], alpha, a, style))

    # Get auto-scaled x range, then plot fit lines to the right edge
    ax.set_xlim(left=1e3)
    x_right = ax.get_xlim()[1]
    for n_start, alpha, a, style in fit_info:
        n_fit = np.logspace(np.log10(n_start), np.log10(x_right), 50)
        m_fit = a * n_fit ** alpha
        ax.loglog(n_fit, m_fit, '--', color=style['color'], linewidth=2,
                  label='%s: $O(N^{%.2f})$' % (style['label'], alpha))

    # 128 GB memory limit line
    mem_limit_mb = 128 * 1024
    ax.axhline(y=mem_limit_mb, color='red', linestyle=':', linewidth=1.5,
               label='128 GB')

    ax.set_xlabel('DOF')
    ax.set_ylabel('Peak memory [MB]')
    ax.legend()
    ax.grid(True, which='both', alpha=0.3)
    fig.tight_layout()
    return fig


def print_summary_table(data):
    """Print summary table of all results."""
    print('=' * 80)
    print('HEXAHEDRON BENCHMARK SUMMARY')
    print('=' * 80)

    print('\n%-8s %8s %8s %10s %8s %10s %10s' % (
        'Solver', 'N_elem', 'NDOF', 't_solve', 'Nonl_it', 'Lin_it', 'Mem_MB'))
    print('-' * 70)

    for solver, results in sorted(data.items()):
        for r in results:
            print('%-8s %8d %8d %10.3f %8d %10d %10.1f' % (
                solver,
                r['n_elements'],
                r['ndof'],
                r['t_solve'],
                r['nonl_iterations'],
                r.get('linear_iterations', 0),
                r.get('solver_memory_mb', 0.0),
            ))

    print('\n--- SCALING EXPONENTS (last %d points, t = a * N^alpha) ---\n' % FIT_LAST_N)
    print('%-10s %10s %10s %10s %10s' % ('Solver', 'Time_exp', 'Time_R2', 'Mem_exp', 'Mem_R2'))
    print('-' * 55)

    for solver, results in sorted(data.items()):
        ndofs = np.array([r['ndof'] for r in results], dtype=float)
        times = np.array([r['t_solve'] for r in results], dtype=float)
        mems = []
        mem_dofs = []
        for r in results:
            mem = r.get('peak_memory_mb')
            if mem is not None and mem > 0:
                mems.append(mem)
                mem_dofs.append(r['ndof'])

        n_fit_pts = min(FIT_LAST_N, len(times))
        if len(times) >= 2:
            alpha_t, _, r2_t = fit_power_law(ndofs[-n_fit_pts:], times[-n_fit_pts:])
        else:
            alpha_t, r2_t = 0, 0

        mem_dofs = np.array(mem_dofs, dtype=float)
        mems = np.array(mems, dtype=float)
        n_fit_mem = min(FIT_LAST_N, len(mems))
        if len(mems) >= 2:
            alpha_m, _, r2_m = fit_power_law(mem_dofs[-n_fit_mem:], mems[-n_fit_mem:])
        else:
            alpha_m, r2_m = 0, 0

        print('%-10s %10.4f %10.5f %10.4f %10.5f' % (solver, alpha_t, r2_t, alpha_m, r2_m))


def main():
    print('Loading hexahedron benchmark results...')
    data = load_results()

    total = sum(len(results) for results in data.values())
    print('Loaded %d result files' % total)

    if total == 0:
        print('No results found. Run benchmarks first.')
        return

    print_summary_table(data)

    fig_time = plot_scaling(data, 't_solve', 'Solve time [s]')
    time_path = os.path.join(SCRIPT_DIR, 'scaling_time.png')
    fig_time.savefig(time_path, dpi=150, bbox_inches='tight')
    print('\nSaved: %s' % time_path)

    fig_mem = plot_memory_scaling(data)
    mem_path = os.path.join(SCRIPT_DIR, 'scaling_memory.png')
    fig_mem.savefig(mem_path, dpi=150, bbox_inches='tight')
    print('Saved: %s' % mem_path)

    plt.close('all')


if __name__ == '__main__':
    main()
