#!/usr/bin/env python
"""
Analytical Reference Solution for Magnetized Cube Field

This script computes Radia's B field from a cube with known magnetization
and compares tetrahedral mesh results with hexahedral (manual mesh) results.

Test Case: Uniformly magnetized cube
- Cube size: 1.0m x 1.0m x 1.0m
- Uniform magnetization: M = [0, 0, 1.0e6] A/m

This evaluates whether rad.Fld from tetrahedral mesh elements gives
the same result as rad.Fld from hexahedral (analytical formula) elements.

Author: Radia Development Team
Date: 2025-12-13
"""
import os
import sys
import json
import time
from pathlib import Path
import numpy as np

# Path setup
_repo_root = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_repo_root / 'build' / 'Release'))
sys.path.append(str(_repo_root / 'src' / 'radia'))

import radia as rad

print('=' * 70)
print('Analytical Reference: Tetra vs Hexa Field Comparison')
print('=' * 70)

# Check for NGSolve
try:
    from netgen.occ import Box, Pnt, OCCGeometry
    from ngsolve import Mesh
    from ngsolve import TaskManager
    NGSOLVE_AVAILABLE = True
except ImportError as e:
    print('NGSolve not available: %s' % e)
    NGSOLVE_AVAILABLE = False
    sys.exit(1)

from netgen_mesh_import import netgen_mesh_to_radia

# =============================================================================
# Parameters
# =============================================================================
MU_0 = 4 * np.pi * 1e-7  # Vacuum permeability [T*m/A]

# Geometry
CUBE_SIZE = 1.0      # Cube edge length [m]
CUBE_HALF = 0.5      # Half size

# Magnetization (permanent magnet, uniform)
M_Z = 1.0e6          # Magnetization magnitude [A/m]
MAGNETIZATION = [0.0, 0.0, M_Z]

# Mesh parameters
TETRA_MAXH = 0.25    # Tetrahedral mesh size
HEXA_NDIV = 4        # Hexahedral subdivisions


def generate_test_points():
    """Generate test points for field evaluation - all outside the cube."""
    points = []
    labels = []

    # On-axis points (z-axis, outside cube: z > 0.5)
    z_values = [0.6, 0.8, 1.0, 1.5, 2.0, 3.0]
    for z in z_values:
        points.append([0.0, 0.0, z])
        labels.append('z-axis (z=%.1f)' % z)

    # Off-axis points (all outside cube)
    off_axis_points = [
        [0.6, 0.0, 0.6],
        [0.0, 0.6, 0.6],
        [0.6, 0.6, 0.6],
        [1.0, 0.0, 1.0],
        [0.0, 1.0, 1.0],
        [1.0, 1.0, 1.0],
        [2.0, 0.0, 2.0],
    ]
    for pt in off_axis_points:
        points.append(pt)
        labels.append('off-axis (%.1f,%.1f,%.1f)' % tuple(pt))

    return points, labels


def hex_vertices(cx, cy, cz, dx, dy, dz):
    """Generate hexahedron vertices from center and dimensions."""
    hx, hy, hz = dx/2, dy/2, dz/2
    return [
        [cx-hx, cy-hy, cz-hz], [cx+hx, cy-hy, cz-hz],
        [cx+hx, cy+hy, cz-hz], [cx-hx, cy+hy, cz-hz],
        [cx-hx, cy-hy, cz+hz], [cx+hx, cy-hy, cz+hz],
        [cx+hx, cy+hy, cz+hz], [cx-hx, cy+hy, cz+hz]
    ]


def compute_hexa_solution(test_points):
    """
    Compute B field using Radia hexahedral (manual mesh) method.
    This uses analytical formulas for rectangular blocks.
    """
    rad.UtiDelAll()

    # Create uniformly magnetized cube: 1m side, centered at origin
    half = CUBE_SIZE / 2
    elem_size = CUBE_SIZE / HEXA_NDIV
    n_elements = HEXA_NDIV ** 3

    # Create cube mesh manually
    elements = []
    for i in range(HEXA_NDIV):
        for j in range(HEXA_NDIV):
            for k in range(HEXA_NDIV):
                cx = -half + (i + 0.5) * elem_size
                cy = -half + (j + 0.5) * elem_size
                cz = -half + (k + 0.5) * elem_size
                vertices = hex_vertices(cx, cy, cz, elem_size, elem_size, elem_size)
                elem = rad.ObjHexahedron(vertices, MAGNETIZATION)
                elements.append(elem)

    cube = rad.ObjCnt(elements)

    # Evaluate field at test points
    B_values = []
    for pt in test_points:
        try:
            B = rad.Fld(cube, 'b', pt)
            B_values.append(list(B))
        except Exception as e:
            print('  Error at point %s: %s' % (pt, e))
            B_values.append([np.nan, np.nan, np.nan])

    return {
        'method': 'hexahedral (manual mesh)',
        'n_elements': n_elements,
        'B_values': B_values
    }


def compute_tetra_solution(test_points):
    """
    Compute B field using Radia tetrahedral (ObjTetrahedron) method.
    Uses Netgen for mesh generation with uniform magnetization.
    """
    rad.UtiDelAll()

    # Generate tetrahedral mesh with Netgen
    cube_solid = Box(Pnt(-CUBE_HALF, -CUBE_HALF, -CUBE_HALF),
                     Pnt(CUBE_HALF, CUBE_HALF, CUBE_HALF))
    cube_solid.mat('magnetic')
    geo = OCCGeometry(cube_solid)
    with TaskManager():
        ngmesh = geo.GenerateMesh(maxh=TETRA_MAXH)
        mesh = Mesh(ngmesh)

        n_elements = mesh.ne

        # Create Radia tetrahedra with UNIFORM magnetization
        mag_obj = netgen_mesh_to_radia(
            mesh,
            material={'magnetization': MAGNETIZATION},
            units='m',
            material_filter='magnetic',
            verbose=False
        )

        # Evaluate field at test points
        B_values = []
        for pt in test_points:
            try:
                B = rad.Fld(mag_obj, 'b', pt)
                B_values.append(list(B))
            except Exception as e:
                print('  Error at point %s: %s' % (pt, e))
                B_values.append([np.nan, np.nan, np.nan])

        return {
            'method': 'tetrahedral (ObjTetrahedron MSC)',
            'n_elements': n_elements,
            'maxh': TETRA_MAXH,
            'B_values': B_values
        }


def main():
    """Run the comparison."""
    print()
    print('Parameters:')
    print('  Cube size:      %.1f m' % CUBE_SIZE)
    print('  Magnetization:  [0, 0, %.0e] A/m (uniform)' % M_Z)
    print()

    # Generate test points
    test_points, labels = generate_test_points()
    print('Test points: %d (all outside the magnetic cube)' % len(test_points))

    # Compute hexahedral reference
    print()
    print('Computing hexahedral reference (n_div=%d)...' % HEXA_NDIV)
    t0 = time.time()
    hexa_result = compute_hexa_solution(test_points)
    t_hexa = time.time() - t0
    print('  Elements: %d, Time: %.3f s' % (hexa_result['n_elements'], t_hexa))

    # Compute tetrahedral solution
    print()
    print('Computing tetrahedral solution (maxh=%.2f)...' % TETRA_MAXH)
    t0 = time.time()
    tetra_result = compute_tetra_solution(test_points)
    t_tetra = time.time() - t0
    print('  Elements: %d, Time: %.3f s' % (tetra_result['n_elements'], t_tetra))

    # ==========================================================================
    # Compare Results
    # ==========================================================================
    print()
    print('=' * 70)
    print('Field Comparison: Tetrahedral vs Hexahedral')
    print('=' * 70)

    print()
    print('Reference: Hexahedral (manual mesh with analytical formula)')
    print()

    print('%s  %s  %s  %s  %s' % (
        'Point'.center(25),
        '|B| Hexa'.center(14),
        '|B| Tetra'.center(14),
        'Error (%)'.center(10),
        'Status'.center(8)
    ))
    print('-' * 75)

    errors = []
    for i, (pt, label) in enumerate(zip(test_points, labels)):
        B_hexa = np.array(hexa_result['B_values'][i])
        B_tetra = np.array(tetra_result['B_values'][i])

        B_hexa_mag = np.linalg.norm(B_hexa)
        B_tetra_mag = np.linalg.norm(B_tetra)

        if B_hexa_mag > 1e-15 and not np.isnan(B_tetra_mag):
            error = abs(B_tetra_mag - B_hexa_mag) / B_hexa_mag * 100
            status = 'OK' if error < 5 else ('WARN' if error < 10 else 'CHECK')
            errors.append(error)
        else:
            error = np.nan
            status = 'NaN'

        print('%s  %14.6e  %14.6e  %10.2f  %8s' % (
            label[:25].ljust(25),
            B_hexa_mag,
            B_tetra_mag,
            error if not np.isnan(error) else 0.0,
            status
        ))

    # ==========================================================================
    # Summary
    # ==========================================================================
    print()
    print('=' * 70)
    print('Summary')
    print('=' * 70)

    if errors:
        avg_error = np.mean(errors)
        max_error = np.max(errors)
        min_error = np.min(errors)

        print()
        print('Field Comparison Statistics:')
        print('  Valid points: %d / %d' % (len(errors), len(test_points)))
        print('  Average error: %.4f%%' % avg_error)
        print('  Maximum error: %.4f%%' % max_error)
        print('  Minimum error: %.4f%%' % min_error)

        if avg_error < 0.01:
            print()
            print('[PASS] Tetrahedral MSC field MATCHES hexahedral (< 0.01%%)')
            print('       rad.Fld from tetrahedral mesh is ACCURATE')
        elif avg_error < 1.0:
            print()
            print('[PASS] Tetrahedral MSC field is GOOD (< 1%%)')
        elif avg_error < 5.0:
            print()
            print('[PASS] Tetrahedral MSC field is ACCEPTABLE (< 5%%)')
        else:
            print()
            print('[CHECK] Tetrahedral MSC field needs investigation')

    # ==========================================================================
    # Save Results
    # ==========================================================================
    output_file = os.path.join(os.path.dirname(__file__), 'analytical_reference_results.json')

    output_data = {
        'parameters': {
            'cube_size': CUBE_SIZE,
            'magnetization': MAGNETIZATION,
            'mu_0': MU_0
        },
        'hexa': {
            'n_div': HEXA_NDIV,
            'n_elements': hexa_result['n_elements'],
            'time': t_hexa,
            'B_values': hexa_result['B_values']
        },
        'tetra': {
            'maxh': TETRA_MAXH,
            'n_elements': tetra_result['n_elements'],
            'time': t_tetra,
            'B_values': tetra_result['B_values']
        },
        'test_points': test_points,
        'labels': labels,
        'errors': errors,
        'avg_error': float(np.mean(errors)) if errors else None
    }

    with open(output_file, 'w') as f:
        json.dump(output_data, f, indent=2)

    print()
    print('Results saved to: %s' % output_file)
    print()
    print('=' * 70)

    return output_data


if __name__ == '__main__':
    results = main()
