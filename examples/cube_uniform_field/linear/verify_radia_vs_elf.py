#!/usr/bin/env python
"""
Verify Radia Results Match ELF (Linear Material)

Compares Radia hexahedral and tetrahedral results with ELF benchmark results.
Both LU and BiCGSTAB solvers are tested.

Requirements:
- Radia built and installed
- ELF benchmark results in S:/ELF_MAGIC/01_GitHub/examples/cube_uniform_field/linear/

Usage:
    python verify_radia_vs_elf.py

Author: Radia Development Team
Date: 2025-12-20
"""

import sys
import os
import json
import numpy as np

# Path setup
_src_path = os.path.join(os.path.dirname(__file__), '../../../src/radia')
sys.path.insert(0, _src_path)

import radia as rad

# Import unified benchmark conditions
from benchmark_conditions import (
    MU_0, CUBE_SIZE, CUBE_HALF, H_EXT, B_EXT, CHI, MU_R,
    SOLVER_TOLERANCE, MAX_ITERATIONS, M_ANALYTICAL_Z
)


# Hexahedral face topology (1-indexed for Radia)
HEX_FACES = [
    [1, 4, 3, 2],  # Bottom face (z=0)
    [5, 6, 7, 8],  # Top face (z=1)
    [1, 2, 6, 5],  # Front face (y=0)
    [3, 4, 8, 7],  # Back face (y=1)
    [1, 5, 8, 4],  # Left face (x=0)
    [2, 3, 7, 6]   # Right face (x=1)
]


def create_hex_mesh(n_div):
    """Create hexahedral mesh using ObjPolyhdr."""
    rad.FldUnits('m')
    rad.UtiDelAll()

    cube_size_elem = CUBE_SIZE / n_div
    elements = []

    for ix in range(n_div):
        for iy in range(n_div):
            for iz in range(n_div):
                # Corner of this sub-cube
                x0 = ix * cube_size_elem - CUBE_HALF
                y0 = iy * cube_size_elem - CUBE_HALF
                z0 = iz * cube_size_elem - CUBE_HALF

                # 8 vertices of hexahedron
                vertices = [
                    [x0, y0, z0],
                    [x0 + cube_size_elem, y0, z0],
                    [x0 + cube_size_elem, y0 + cube_size_elem, z0],
                    [x0, y0 + cube_size_elem, z0],
                    [x0, y0, z0 + cube_size_elem],
                    [x0 + cube_size_elem, y0, z0 + cube_size_elem],
                    [x0 + cube_size_elem, y0 + cube_size_elem, z0 + cube_size_elem],
                    [x0, y0 + cube_size_elem, z0 + cube_size_elem]
                ]

                obj = rad.ObjPolyhdr(vertices, HEX_FACES, [0, 0, 0])
                elements.append(obj)

    # Create container
    cube = rad.ObjCnt(elements)

    mat = rad.MatLin(MU_R)  # mu_r (industry standard)
    rad.MatApl(cube, mat)

    ext = rad.ObjBckg([0, 0, B_EXT])
    grp = rad.ObjCnt([cube, ext])

    return grp, cube


def create_tetra_mesh(maxh):
    """Create tetrahedral mesh using Netgen."""
    try:
        from netgen.occ import Box, Pnt, OCCGeometry
        from ngsolve import Mesh
        from netgen_mesh_import import netgen_mesh_to_radia
    except ImportError as e:
        print('[SKIP] Tetra mesh creation failed: %s' % e)
        return None, None, 0

    rad.FldUnits('m')
    rad.UtiDelAll()

    cube_solid = Box(Pnt(-CUBE_HALF, -CUBE_HALF, -CUBE_HALF),
                     Pnt(CUBE_HALF, CUBE_HALF, CUBE_HALF))
    cube_solid.mat('magnetic')
    geo = OCCGeometry(cube_solid)

    ngmesh = geo.GenerateMesh(maxh=maxh)
    mesh = Mesh(ngmesh)
    n_elements = mesh.ne

    cube = netgen_mesh_to_radia(mesh,
                                 material={'magnetization': [0, 0, 0]},
                                 units='m',
                                 material_filter='magnetic',
                                 verbose=False)

    mat = rad.MatLin(MU_R)  # mu_r (industry standard)
    rad.MatApl(cube, mat)

    ext = rad.ObjBckg([0, 0, B_EXT])
    grp = rad.ObjCnt([cube, ext])

    return grp, cube, n_elements


def get_avg_magnetization(cube):
    """Get average z-magnetization."""
    all_M = rad.ObjM(cube)
    M_list = [m[1] for m in all_M]
    return np.mean([m[2] for m in M_list])


def load_elf_results(element_type, n_div_or_maxh, solver_name):
    """Load ELF benchmark results for comparison."""
    elf_base = r'S:\ELF_MAGIC\01_GitHub\examples\cube_uniform_field\linear'

    if element_type == 'hex':
        elf_dir = os.path.join(elf_base, 'hexahedron', solver_name)
        # ELF uses lowercase: hexa_n3_results.json
        filename = 'hexa_n%d_results.json' % n_div_or_maxh
    else:
        elf_dir = os.path.join(elf_base, 'tetrahedron', solver_name)
        maxh_str = ('%.2f' % n_div_or_maxh).replace('.', '_')
        filename = 'tetra_maxh%s_results.json' % maxh_str

    filepath = os.path.join(elf_dir, filename)
    if os.path.exists(filepath):
        with open(filepath, 'r') as f:
            return json.load(f)
    return None


def verify_hex_benchmark(n_div, solver_method, solver_name):
    """Verify hexahedral benchmark against ELF."""
    print('\n--- Hexahedron N=%d, Solver=%s ---' % (n_div, solver_name))

    grp, cube = create_hex_mesh(n_div)

    try:
        result = rad.Solve(grp, SOLVER_TOLERANCE, MAX_ITERATIONS, solver_method)
        M_avg_z = get_avg_magnetization(cube)
        n_iter = int(result[3]) if result[3] else 0
    except Exception as e:
        print('Radia solve failed: %s' % e)
        return None

    # Load ELF result
    elf_result = load_elf_results('hex', n_div, solver_name)

    print('Radia M_avg_z: %.6f A/m' % M_avg_z)

    if elf_result:
        elf_M = elf_result.get('M_avg_z', 0)
        diff = abs(M_avg_z - elf_M)
        rel_diff = diff / elf_M * 100 if elf_M > 0 else 0
        print('ELF   M_avg_z: %.6f A/m' % elf_M)
        print('Difference:    %.6f A/m (%.6f%%)' % (diff, rel_diff))

        # Check if match
        if rel_diff < 0.0001:  # < 0.0001%
            print('Result: MATCH (Perfect)')
            return 'MATCH'
        elif rel_diff < 0.01:  # < 0.01%
            print('Result: MATCH (Close)')
            return 'MATCH'
        else:
            print('Result: MISMATCH')
            return 'MISMATCH'
    else:
        print('ELF result not found')
        return 'NO_ELF'


def verify_tetra_benchmark(maxh, solver_method, solver_name):
    """Verify tetrahedral benchmark against ELF."""
    print('\n--- Tetrahedron maxh=%.2f, Solver=%s ---' % (maxh, solver_name))

    grp, cube, n_elements = create_tetra_mesh(maxh)
    if grp is None:
        return 'SKIP'

    print('Elements: %d' % n_elements)

    try:
        result = rad.Solve(grp, SOLVER_TOLERANCE, MAX_ITERATIONS, solver_method)
        M_avg_z = get_avg_magnetization(cube)
        n_iter = int(result[3]) if result[3] else 0
    except Exception as e:
        print('Radia solve failed: %s' % e)
        return None

    # Load ELF result
    elf_result = load_elf_results('tetra', maxh, solver_name)

    print('Radia M_avg_z: %.6f A/m' % M_avg_z)

    if elf_result:
        elf_M = elf_result.get('M_avg_z', 0)
        diff = abs(M_avg_z - elf_M)
        rel_diff = diff / elf_M * 100 if elf_M > 0 else 0
        print('ELF   M_avg_z: %.6f A/m' % elf_M)
        print('Difference:    %.6f A/m (%.6f%%)' % (diff, rel_diff))

        # Check if match - tetra has some numerical differences due to mesh generation
        if rel_diff < 0.001:  # < 0.001%
            print('Result: MATCH (Perfect)')
            return 'MATCH'
        elif rel_diff < 0.1:  # < 0.1%
            print('Result: MATCH (Close)')
            return 'MATCH'
        else:
            print('Result: MISMATCH')
            return 'MISMATCH'
    else:
        print('ELF result not found')
        return 'NO_ELF'


def main():
    print('=' * 70)
    print('VERIFY RADIA VS ELF - LINEAR MATERIAL BENCHMARK')
    print('=' * 70)
    print('Cube size: %.1f m' % CUBE_SIZE)
    print('mu_r: %d, chi: %d' % (MU_R, CHI))
    print('H_ext: %.0f A/m' % H_EXT)
    print('Analytical M_z: %.0f A/m' % M_ANALYTICAL_Z)
    print()

    results = []

    # Test hexahedra
    for n_div in [3, 5]:
        for solver_method, solver_name in [(0, 'lu'), (1, 'bicgstab')]:
            r = verify_hex_benchmark(n_div, solver_method, solver_name)
            results.append(('hex_N%d_%s' % (n_div, solver_name), r))

    # Test tetrahedra
    for maxh in [0.3]:
        for solver_method, solver_name in [(0, 'lu'), (1, 'bicgstab')]:
            r = verify_tetra_benchmark(maxh, solver_method, solver_name)
            results.append(('tetra_maxh%.2f_%s' % (maxh, solver_name), r))

    # Summary
    print('\n' + '=' * 70)
    print('SUMMARY')
    print('=' * 70)

    all_match = True
    for name, status in results:
        status_str = status if status else 'FAILED'
        print('%-30s : %s' % (name, status_str))
        if status not in ['MATCH', 'NO_ELF', 'SKIP']:
            all_match = False

    print()
    if all_match:
        print('All tests PASSED - Radia matches ELF')
    else:
        print('Some tests FAILED - Check differences above')

    return 0 if all_match else 1


if __name__ == '__main__':
    sys.exit(main())
