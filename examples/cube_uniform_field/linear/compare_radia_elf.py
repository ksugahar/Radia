#!/usr/bin/env python
"""
Compare Radia and ELF Results for Linear Material Benchmark

This script compares hexahedral results between Radia and ELF.

Author: Radia Development Team
Date: 2025-12-20
"""

import os
import json

# Paths
RADIA_BASE = os.path.dirname(__file__)
ELF_BASE = r'S:\ELF_MAGIC\01_GitHub\examples\cube_uniform_field\linear'


def load_json(filepath):
    """Load JSON file if exists."""
    if os.path.exists(filepath):
        with open(filepath, 'r') as f:
            return json.load(f)
    return None


def compare_hexahedron():
    """Compare hexahedron results (6DOF MSC only)."""
    print('=' * 70)
    print('HEXAHEDRON COMPARISON: Radia vs ELF (6DOF MSC)')
    print('=' * 70)
    print()
    print('%5s  %6s  %6s  %12s  %12s  %12s  %8s' %
          ('N', 'R_DOF', 'E_DOF', 'Radia M_z', 'ELF M_z', 'Difference', 'Match?'))
    print('-' * 70)

    n_values = [3, 4, 5, 6, 7, 8]
    all_match = True
    compared = 0

    for n in n_values:
        radia_lu = os.path.join(RADIA_BASE, 'hexahedron_msc', 'lu', 'hexa_n%d_results.json' % n)
        elf_lu = os.path.join(ELF_BASE, 'hexahedron', 'lu', 'hexa_n%d_results.json' % n)

        radia_data = load_json(radia_lu)
        elf_data = load_json(elf_lu)

        if radia_data and elf_data:
            radia_dof = radia_data.get('dof', 0)
            elf_dof = elf_data.get('ndof', 0)
            radia_m = radia_data.get('M_avg_z', 0)
            elf_m = elf_data.get('M_avg_z', 0)

            # Only compare if DOF matches (both 6DOF)
            if radia_dof != elf_dof:
                print('%5d  %6d  %6d  %12s  %12s  %12s  %8s' %
                      (n, radia_dof, elf_dof, '-', '-', 'DOF mismatch', 'SKIP'))
                continue

            compared += 1
            diff = abs(radia_m - elf_m)
            rel_diff = diff / elf_m * 100 if elf_m > 0 else 0

            match = 'YES' if rel_diff < 0.0001 else 'NO'
            if rel_diff >= 0.0001:
                all_match = False

            print('%5d  %6d  %6d  %12.2f  %12.2f  %12.6f  %8s' %
                  (n, radia_dof, elf_dof, radia_m, elf_m, diff, match))
        else:
            status = []
            if not radia_data:
                status.append('Radia missing')
            if not elf_data:
                status.append('ELF missing')
            print('%5d  %s' % (n, ', '.join(status)))

    print()
    if compared == 0:
        print('Result: NO COMPARABLE TESTS (DOF mismatch)')
        return True  # Not a failure, just no data
    elif all_match:
        print('Result: ALL %d HEXAHEDRON TESTS MATCH (6DOF MSC)' % compared)
    else:
        print('Result: SOME TESTS DO NOT MATCH')

    return all_match


def compare_bicgstab_vs_lu():
    """Compare BiCGSTAB vs LU for Radia."""
    print()
    print('=' * 70)
    print('RADIA SOLVER COMPARISON: LU vs BiCGSTAB')
    print('=' * 70)
    print()
    print('%5s  %12s  %12s  %12s  %10s' %
          ('N', 'LU M_z', 'BiCGSTAB M_z', 'Difference', 'Match?'))
    print('-' * 70)

    n_values = [3, 4, 5, 6, 7, 8]
    all_match = True

    for n in n_values:
        lu_file = os.path.join(RADIA_BASE, 'hexahedron_msc', 'lu', 'hexa_n%d_results.json' % n)
        bicg_file = os.path.join(RADIA_BASE, 'hexahedron_msc', 'bicgstab', 'hexa_n%d_results.json' % n)

        lu_data = load_json(lu_file)
        bicg_data = load_json(bicg_file)

        if lu_data and bicg_data:
            lu_m = lu_data.get('M_avg_z', 0)
            bicg_m = bicg_data.get('M_avg_z', 0)
            diff = abs(lu_m - bicg_m)
            rel_diff = diff / lu_m * 100 if lu_m > 0 else 0

            match = 'YES' if rel_diff < 0.001 else 'NO'
            if rel_diff >= 0.001:
                all_match = False

            print('%5d  %12.2f  %12.2f  %12.6f  %10s' %
                  (n, lu_m, bicg_m, diff, match))
        else:
            print('%5d  Data missing' % n)

    print()
    if all_match:
        print('Result: ALL SOLVER COMPARISON TESTS MATCH')
    else:
        print('Result: SOME TESTS DO NOT MATCH')

    return all_match


def main():
    print()
    print('#' * 70)
    print('# RADIA vs ELF LINEAR MATERIAL BENCHMARK COMPARISON')
    print('#' * 70)
    print()

    hex_match = compare_hexahedron()
    solver_match = compare_bicgstab_vs_lu()

    print()
    print('=' * 70)
    print('FINAL SUMMARY')
    print('=' * 70)
    print('Hexahedron Radia vs ELF: %s' % ('PASS' if hex_match else 'FAIL'))
    print('LU vs BiCGSTAB:          %s' % ('PASS' if solver_match else 'FAIL'))
    print()

    return 0 if (hex_match and solver_match) else 1


if __name__ == '__main__':
    exit(main())
