#!/usr/bin/env python
"""
Test MSC Element Solvers: LU vs BiCGSTAB
Tests both hexahedral MSC and tetrahedral MSC elements.
"""

import sys
import os
sys.path.insert(0, os.path.join('s:/Radia/01_GitHub', 'build/Release'))
sys.path.append(os.path.join('s:/Radia/01_GitHub', 'src/radia'))
import numpy as np
import radia as rad
import time

print('='*70)
print('MSC Element Solver Comparison: LU vs BiCGSTAB')
print('='*70)

rad.FldUnits('m')

chi = 999.0
H_ext = 50000.0
MU_0 = 4 * np.pi * 1e-7
B_ext = MU_0 * H_ext

N_cube = 1.0/3.0
M_analytical = chi * H_ext / (1 + chi * N_cube)
print(f'Analytical M_z (linear): {M_analytical:,.0f} A/m')
print()

# =====================================================================
# 1. Hexahedral MSC (ObjPolyhdr with 6 faces)
# =====================================================================
print('='*70)
print('1. HEXAHEDRAL MSC (ObjPolyhdr - 6 quadrilateral faces)')
print('='*70)

# Hexahedron face definition (6 faces, each with 4 vertices - 1-based indexing!)
# Radia ObjPolyhdr uses 1-based vertex indices
HEX_FACES = [
    [1, 4, 3, 2],  # Bottom (z-)
    [5, 6, 7, 8],  # Top (z+)
    [1, 2, 6, 5],  # Front (y-)
    [3, 4, 8, 7],  # Back (y+)
    [1, 5, 8, 4],  # Left (x-)
    [2, 3, 7, 6],  # Right (x+)
]

def create_hex_mesh_msc(n_div):
    '''Create hexahedral mesh using ObjPolyhdr (MSC method)'''
    elements = []
    step = 1.0 / n_div
    half = 0.5

    for iz in range(n_div):
        for iy in range(n_div):
            for ix in range(n_div):
                x0 = -half + ix * step
                y0 = -half + iy * step
                z0 = -half + iz * step
                x1 = x0 + step
                y1 = y0 + step
                z1 = z0 + step

                vertices = [
                    [x0, y0, z0], [x1, y0, z0], [x1, y1, z0], [x0, y1, z0],
                    [x0, y0, z1], [x1, y0, z1], [x1, y1, z1], [x0, y1, z1],
                ]

                elem = rad.ObjPolyhdr(vertices, HEX_FACES, [0, 0, 0])
                elements.append(elem)

    return rad.ObjCnt(elements)

for n_div in [3]:
    print(f'\n--- n_div={n_div} ({n_div**3} hex elements) ---')
    results = []

    for method, name in [(0, 'LU'), (1, 'BiCGSTAB')]:
        rad.UtiDelAll()

        hex_mesh = create_hex_mesh_msc(n_div)
        mat = rad.MatLin(chi)
        rad.MatApl(hex_mesh, mat)

        ext = rad.ObjBckg([0, 0, B_ext])
        grp = rad.ObjCnt([hex_mesh, ext])

        t0 = time.time()
        res = rad.Solve(grp, 0.0001, 1000, method)
        dt = time.time() - t0

        M_raw = rad.ObjM(hex_mesh)
        M_z_avg = np.mean([elem[1][2] for elem in M_raw])
        results.append(M_z_avg)

        print(f'  {name:10s}: M_z = {M_z_avg:,.0f} A/m, iter={int(res[3]):3d}, time={dt:.4f}s')

    if len(results) == 2:
        diff = abs(results[0] - results[1]) / max(abs(results[0]), 1) * 100
        print(f'  Difference: {diff:.4f}%')

# =====================================================================
# 2. Tetrahedral MSC (ObjPolyhdr with 4 faces)
# =====================================================================
print()
print('='*70)
print('2. TETRAHEDRAL MSC (ObjPolyhdr - 4 triangular faces)')
print('='*70)

try:
    from netgen_mesh_import import netgen_mesh_to_radia
    from netgen.occ import Box, Pnt, OCCGeometry
    from ngsolve import Mesh

    for maxh in [0.4]:
        print(f'\n--- maxh={maxh} ---')

        # Generate tetrahedral mesh
        cube_solid = Box(Pnt(-0.5, -0.5, -0.5), Pnt(0.5, 0.5, 0.5))
        cube_solid.mat('magnetic')
        geo = OCCGeometry(cube_solid)
        ngmesh = geo.GenerateMesh(maxh=maxh)
        mesh = Mesh(ngmesh)
        n_elem = mesh.ne
        print(f'  Elements: {n_elem}')

        results = []

        for method, name in [(0, 'LU'), (1, 'BiCGSTAB')]:
            rad.UtiDelAll()

            mag_obj = netgen_mesh_to_radia(mesh,
                                           material={'magnetization': [0, 0, 0]},
                                           units='m',
                                           material_filter='magnetic',
                                           verbose=False)

            mat = rad.MatLin(chi)
            rad.MatApl(mag_obj, mat)

            ext = rad.ObjBckg([0, 0, B_ext])
            grp = rad.ObjCnt([mag_obj, ext])

            t0 = time.time()
            res = rad.Solve(grp, 0.0001, 1000, method)
            dt = time.time() - t0

            M_raw = rad.ObjM(mag_obj)
            M_z_avg = np.mean([elem[1][2] for elem in M_raw])
            results.append(M_z_avg)

            print(f'  {name:10s}: M_z = {M_z_avg:,.0f} A/m, iter={int(res[3]):3d}, time={dt:.4f}s')

        if len(results) == 2:
            diff = abs(results[0] - results[1]) / max(abs(results[0]), 1) * 100
            print(f'  Difference: {diff:.4f}%')

except ImportError as e:
    print(f'  NGSolve/Netgen not available: {e}')

print()
print('='*70)
print('SUMMARY')
print('='*70)
