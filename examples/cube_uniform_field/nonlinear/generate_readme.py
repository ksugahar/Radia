#!/usr/bin/env python
"""
README Generator for Benchmark Results

Generates README.md from JSON benchmark results.
Works for both ELF and Radia results directories.

Usage:
    python generate_readme.py                    # Generate README.md for Radia results
    python generate_readme.py --elf              # Generate README.md for ELF results (run from ELF directory)
    python generate_readme.py --compare          # Generate comparison between ELF and Radia
"""

import os
import sys
import json
import argparse
from glob import glob
from datetime import datetime


def load_results(base_dir, element_type, solver_type):
    """Load all JSON results for a given element type and solver."""
    if element_type == 'tetra':
        pattern = os.path.join(base_dir, 'tetrahedron', solver_type, 'tetra_*.json')
    else:
        pattern = os.path.join(base_dir, 'hexahedron', solver_type, 'hex_*.json')

    results = []
    for filepath in sorted(glob(pattern)):
        try:
            with open(filepath, 'r') as f:
                data = json.load(f)
                results.append(data)
        except (json.JSONDecodeError, IOError) as e:
            print('Warning: Failed to load %s: %s' % (filepath, e))

    return results


def format_table_row(r, include_hmatrix=False):
    """Format a result as a table row."""
    if r['element_type'] == 'hex':
        mesh_col = 'N=%d' % r.get('n_div', 0)
    else:
        mesh_col = 'maxh=%.2fm' % r.get('maxh', 0.0)

    row = '| %s | %d | %d | %.3f | %d | %d | %.0f | %s |' % (
        mesh_col,
        r['n_elements'],
        r['ndof'],
        r['t_solve'],
        r['nonl_iterations'],
        r.get('linear_iterations', 0),
        r['M_avg_z'],
        'Yes' if r['converged'] else 'No'
    )

    if include_hmatrix:
        hm = r.get('hmatrix', {})
        row = row[:-1] + ' %.4f | %d |' % (
            hm.get('compression_ratio', 0.0),
            hm.get('nlf', 0)
        )

    return row


def generate_solver_section(results, solver_name, include_hmatrix=False):
    """Generate markdown section for a solver."""
    if not results:
        return ''

    lines = []
    lines.append('### %s Solver\n' % solver_name)

    # Header
    header = '| Mesh | Elements | DOF | Time (s) | Nonl Iter | Linear Iter | M_avg_z (A/m) | Converged |'
    separator = '|------|----------|-----|----------|-----------|-------------|---------------|-----------|'

    if include_hmatrix:
        header += ' Compression | Leaves |'
        separator += '-------------|--------|'

    lines.append(header)
    lines.append(separator)

    for r in results:
        lines.append(format_table_row(r, include_hmatrix))

    lines.append('')
    return '\n'.join(lines)


def generate_readme(base_dir, title='Benchmark Results'):
    """Generate README.md from benchmark results."""
    lines = []
    lines.append('# %s\n' % title)
    lines.append('Generated: %s\n' % datetime.now().strftime('%Y-%m-%d %H:%M:%S'))

    # Problem description
    lines.append('## Problem Description\n')
    lines.append('- **Geometry**: 1.0 m x 1.0 m x 1.0 m cube')
    lines.append('- **External field**: H_ext = 200,000 A/m (z-direction)')
    lines.append('- **Material**: Nonlinear B-H curve (soft iron saturation)')
    lines.append('- **Convergence**: nonl_tol = 0.001, bicg_tol = 1e-4')
    lines.append('')

    # Hexahedron results
    lines.append('## Hexahedron Results\n')

    hex_lu = load_results(base_dir, 'hex', 'lu')
    hex_bicgstab = load_results(base_dir, 'hex', 'bicgstab')
    hex_hacapk = load_results(base_dir, 'hex', 'hacapk')

    if hex_lu:
        lines.append(generate_solver_section(hex_lu, 'LU'))
    if hex_bicgstab:
        lines.append(generate_solver_section(hex_bicgstab, 'BiCGSTAB'))
    if hex_hacapk:
        lines.append(generate_solver_section(hex_hacapk, 'HACApK (H-matrix)', include_hmatrix=True))

    if not (hex_lu or hex_bicgstab or hex_hacapk):
        lines.append('No hexahedron results found.\n')

    # Tetrahedron results
    lines.append('## Tetrahedron Results\n')

    tetra_lu = load_results(base_dir, 'tetra', 'lu')
    tetra_bicgstab = load_results(base_dir, 'tetra', 'bicgstab')
    tetra_hacapk = load_results(base_dir, 'tetra', 'hacapk')

    if tetra_lu:
        lines.append(generate_solver_section(tetra_lu, 'LU'))
    if tetra_bicgstab:
        lines.append(generate_solver_section(tetra_bicgstab, 'BiCGSTAB'))
    if tetra_hacapk:
        lines.append(generate_solver_section(tetra_hacapk, 'HACApK (H-matrix)', include_hmatrix=True))

    if not (tetra_lu or tetra_bicgstab or tetra_hacapk):
        lines.append('No tetrahedron results found.\n')

    return '\n'.join(lines)


def generate_comparison(radia_dir, elf_dir):
    """Generate comparison README between Radia and ELF results."""
    lines = []
    lines.append('# Radia vs ELF Benchmark Comparison\n')
    lines.append('Generated: %s\n' % datetime.now().strftime('%Y-%m-%d %H:%M:%S'))

    # Problem description
    lines.append('## Problem Description\n')
    lines.append('- **Geometry**: 1.0 m x 1.0 m x 1.0 m cube')
    lines.append('- **External field**: H_ext = 200,000 A/m (z-direction)')
    lines.append('- **Material**: Nonlinear B-H curve (soft iron saturation)')
    lines.append('- **Convergence**: nonl_tol = 0.001, bicg_tol = 1e-4')
    lines.append('')

    # Load all results
    for element_type, element_name in [('hex', 'Hexahedron'), ('tetra', 'Tetrahedron')]:
        lines.append('## %s Comparison\n' % element_name)

        for solver_type, solver_name in [('lu', 'LU'), ('bicgstab', 'BiCGSTAB'), ('hacapk', 'HACApK')]:
            radia_results = load_results(radia_dir, element_type, solver_type)
            elf_results = load_results(elf_dir, element_type, solver_type)

            if not (radia_results or elf_results):
                continue

            lines.append('### %s Solver\n' % solver_name)

            # Create comparison table
            header = '| Mesh | Library | Elements | DOF | Time (s) | Nonl Iter | Linear Iter | M_avg_z (A/m) |'
            separator = '|------|---------|----------|-----|----------|-----------|-------------|---------------|'

            lines.append(header)
            lines.append(separator)

            # Match results by mesh description
            radia_by_mesh = {r['mesh_description']: r for r in radia_results}
            elf_by_mesh = {r['mesh_description']: r for r in elf_results}

            all_meshes = sorted(set(list(radia_by_mesh.keys()) + list(elf_by_mesh.keys())))

            for mesh in all_meshes:
                if mesh in radia_by_mesh:
                    r = radia_by_mesh[mesh]
                    lines.append('| %s | Radia | %d | %d | %.3f | %d | %d | %.0f |' % (
                        mesh, r['n_elements'], r['ndof'], r['t_solve'],
                        r['nonl_iterations'], r.get('linear_iterations', 0), r['M_avg_z']
                    ))

                if mesh in elf_by_mesh:
                    r = elf_by_mesh[mesh]
                    lines.append('| %s | ELF | %d | %d | %.3f | %d | %d | %.0f |' % (
                        mesh, r['n_elements'], r['ndof'], r['t_solve'],
                        r['nonl_iterations'], r.get('linear_iterations', 0), r['M_avg_z']
                    ))

                # Add difference if both exist
                if mesh in radia_by_mesh and mesh in elf_by_mesh:
                    r_radia = radia_by_mesh[mesh]
                    r_elf = elf_by_mesh[mesh]
                    m_diff_pct = 100.0 * (r_radia['M_avg_z'] - r_elf['M_avg_z']) / r_elf['M_avg_z'] if r_elf['M_avg_z'] != 0 else 0
                    t_ratio = r_radia['t_solve'] / r_elf['t_solve'] if r_elf['t_solve'] > 0 else 0
                    lines.append('| %s | **Diff** | - | - | %.2fx | - | - | %+.4f%% |' % (
                        mesh, t_ratio, m_diff_pct
                    ))

            lines.append('')

    return '\n'.join(lines)


def main():
    parser = argparse.ArgumentParser(description='Generate README from benchmark results')
    parser.add_argument('--elf', action='store_true', help='Generate for ELF results')
    parser.add_argument('--compare', action='store_true', help='Generate comparison between ELF and Radia')
    parser.add_argument('--radia-dir', type=str, default=None, help='Radia results directory')
    parser.add_argument('--elf-dir', type=str, default=None, help='ELF results directory')
    parser.add_argument('--output', type=str, default='README.md', help='Output filename')
    args = parser.parse_args()

    script_dir = os.path.dirname(os.path.abspath(__file__))

    if args.compare:
        radia_dir = args.radia_dir or script_dir
        elf_dir = args.elf_dir or os.path.join(script_dir, '..', '..', '..', '..', 'ELF_MAGIC', '01_GitHub', 'examples', 'cube_uniform_field', 'nonlinear')

        if not os.path.exists(elf_dir):
            # Try alternate path
            elf_dir = r's:\ELF_MAGIC\01_GitHub\examples\cube_uniform_field\nonlinear'

        content = generate_comparison(radia_dir, elf_dir)
        title = 'Radia vs ELF Comparison'
    elif args.elf:
        elf_dir = args.elf_dir or script_dir
        content = generate_readme(elf_dir, 'ELF Benchmark Results')
        title = 'ELF Results'
    else:
        radia_dir = args.radia_dir or script_dir
        content = generate_readme(radia_dir, 'Radia Benchmark Results')
        title = 'Radia Results'

    output_path = os.path.join(script_dir, args.output)
    with open(output_path, 'w') as f:
        f.write(content)

    print('Generated: %s' % output_path)


if __name__ == '__main__':
    main()
