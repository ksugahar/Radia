#!/usr/bin/env python
"""
Tetrahedral vs Hexahedral Field Accuracy Evaluation

This script evaluates the accuracy of magnetic field computation from
tetrahedral (Netgen) and hexahedral (manual mesh) meshes in Radia.

Both mesh types are solved using Radia's MMM solver to obtain magnetization,
then the B field is compared at external test points.

Test Case: Linear magnetic cube in uniform external field
- Cube size: 1.0m x 1.0m x 1.0m
- Linear material: mu_r = 1000 (chi = 999)
- External field: H_ext = 50,000 A/m (z-direction)

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
print('Tetrahedral vs Hexahedral Field Accuracy Evaluation')
print('=' * 70)

# Check for NGSolve (for Netgen mesh generation)
try:
    from netgen.occ import Box, Pnt, OCCGeometry
    from ngsolve import Mesh
    from ngsolve import TaskManager
    NETGEN_AVAILABLE = True
except ImportError as e:
    print('NGSolve/Netgen not available: %s' % e)
    NETGEN_AVAILABLE = False

from netgen_mesh_import import netgen_mesh_to_radia

# =============================================================================
# Parameters (matching cube_uniform_field/linear benchmarks)
# =============================================================================
MU_0 = 4 * np.pi * 1e-7  # Vacuum permeability [T*m/A]

# Geometry
CUBE_SIZE = 1.0      # Cube edge length [m]
CUBE_HALF = 0.5      # Half size

# Material (Linear)
MU_R = 1000          # Relative permeability (industry standard input)
CHI = MU_R - 1       # For internal calculations

# External field
H_EXT = 50000.0      # External H field [A/m]
B_EXT = MU_0 * H_EXT  # External B field [T]

# Solver parameters
SOLVER_TOLERANCE = 1e-6
MAX_ITERATIONS = 1000
SOLVER_METHOD = 1  # BiCGSTAB

# Mesh parameters
TETRA_MAXH = 0.3     # Tetrahedral mesh size
HEXA_NDIV = 5        # Hexahedral subdivisions

# Analytical solution for comparison
# M_z = chi * H_int, where H_int ~ H_ext (for high mu_r)
# Approximate: M_z ~ chi * H_ext * 3 / (mu_r + 2)  (demagnetization)
# For mu_r >> 1: M_z ~ chi * H_ext * 3 / mu_r = 3 * H_ext * (mu_r - 1) / mu_r
M_ANALYTICAL_Z = CHI * H_EXT * 3.0 / (MU_R + 2.0)  # ~149,850 A/m


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


def create_hexa_solution(test_points):
    """
    Create hexahedral mesh manually and solve with Radia.
    """
    rad.UtiDelAll()

    print()
    print('Creating hexahedral mesh (n_div=%d)...' % HEXA_NDIV)

    # Create cube mesh manually: HEXA_NDIV x HEXA_NDIV x HEXA_NDIV elements
    half = CUBE_SIZE / 2
    elem_size = CUBE_SIZE / HEXA_NDIV
    n_elements = HEXA_NDIV ** 3

    elements = []
    for i in range(HEXA_NDIV):
        for j in range(HEXA_NDIV):
            for k in range(HEXA_NDIV):
                cx = -half + (i + 0.5) * elem_size
                cy = -half + (j + 0.5) * elem_size
                cz = -half + (k + 0.5) * elem_size
                vertices = hex_vertices(cx, cy, cz, elem_size, elem_size, elem_size)
                elem = rad.ObjHexahedron(vertices, [0, 0, 0])
                elements.append(elem)

    cube = rad.ObjCnt(elements)

    # Apply linear material
    mat = rad.MatLin(MU_R)  # relative permeability
    rad.MatApl(cube, mat)

    # External field
    ext = rad.ObjBckg(lambda p: [0, 0, B_EXT])
    grp = rad.ObjCnt([cube, ext])

    # Solve
    print('  Solving with BiCGSTAB...')
    t0 = time.time()
    result = rad.Solve(grp, SOLVER_TOLERANCE, MAX_ITERATIONS, SOLVER_METHOD)
    t_solve = time.time() - t0

    n_iter = int(result[3]) if result[3] else 0

    # Get magnetization
    all_M = rad.ObjM(cube)
    M_list = [m[1] for m in all_M]
    M_avg_z = np.mean([m[2] for m in M_list])

    print('  Elements: %d' % n_elements)
    print('  Solve time: %.3f s' % t_solve)
    print('  Iterations: %d' % n_iter)
    print('  M_avg_z: %.0f A/m' % M_avg_z)
    print('  Analytical M_z: %.0f A/m' % M_ANALYTICAL_Z)
    print('  M error: %.2f%%' % (abs(M_avg_z - M_ANALYTICAL_Z) / M_ANALYTICAL_Z * 100))

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
        'n_div': HEXA_NDIV,
        't_solve': t_solve,
        'iterations': n_iter,
        'M_avg_z': M_avg_z,
        'B_values': B_values
    }


def create_tetra_solution(test_points):
    """
    Create tetrahedral mesh using Netgen and solve with Radia.
    """
    if not NETGEN_AVAILABLE:
        print('Netgen not available - skipping tetrahedral test')
        return None

    rad.UtiDelAll()

    print()
    print('Creating tetrahedral mesh (maxh=%.2f)...' % TETRA_MAXH)

    # Generate tetrahedral mesh with Netgen
    cube_solid = Box(Pnt(-CUBE_HALF, -CUBE_HALF, -CUBE_HALF),
                     Pnt(CUBE_HALF, CUBE_HALF, CUBE_HALF))
    cube_solid.mat('magnetic')
    geo = OCCGeometry(cube_solid)
    with TaskManager():
        ngmesh = geo.GenerateMesh(maxh=TETRA_MAXH)
        mesh = Mesh(ngmesh)

        n_elements = mesh.ne

        # Import to Radia with initial zero magnetization
        cube = netgen_mesh_to_radia(mesh,
                                     material={'magnetization': [0, 0, 0]},
                                     units='m',
                                     material_filter='magnetic',
                                     verbose=False)

        # Apply linear material
        mat = rad.MatLin(MU_R)  # relative permeability
        rad.MatApl(cube, mat)

        # External field
        ext = rad.ObjBckg(lambda p: [0, 0, B_EXT])
        grp = rad.ObjCnt([cube, ext])

        # Solve
        print('  Solving with BiCGSTAB...')
        t0 = time.time()
        result = rad.Solve(grp, SOLVER_TOLERANCE, MAX_ITERATIONS, SOLVER_METHOD)
        t_solve = time.time() - t0

        n_iter = int(result[3]) if result[3] else 0

        # Get magnetization
        all_M = rad.ObjM(cube)
        M_list = [m[1] for m in all_M]
        M_avg_z = np.mean([m[2] for m in M_list])

        print('  Elements: %d' % n_elements)
        print('  Solve time: %.3f s' % t_solve)
        print('  Iterations: %d' % n_iter)
        print('  M_avg_z: %.0f A/m' % M_avg_z)
        print('  Analytical M_z: %.0f A/m' % M_ANALYTICAL_Z)
        print('  M error: %.2f%%' % (abs(M_avg_z - M_ANALYTICAL_Z) / M_ANALYTICAL_Z * 100))

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
            'method': 'tetrahedral (Netgen + ObjTetrahedron)',
            'n_elements': n_elements,
            'maxh': TETRA_MAXH,
            't_solve': t_solve,
            'iterations': n_iter,
            'M_avg_z': M_avg_z,
            'B_values': B_values
        }


def compare_results(hexa_result, tetra_result, test_points, labels):
    """
    Compare B field from hexahedral and tetrahedral meshes.
    """
    print()
    print('=' * 70)
    print('Field Comparison: Tetrahedral vs Hexahedral')
    print('=' * 70)

    print()
    print('Reference: Hexahedral (manual mesh)')
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

    return errors


def main():
    """Run the complete evaluation."""
    print()
    print('Parameters:')
    print('  Cube size:      %.1f m' % CUBE_SIZE)
    print('  mu_r:           %d (chi = %d)' % (MU_R, CHI))
    print('  H_ext:          %.0f A/m (z-direction)' % H_EXT)
    print('  B_ext:          %.6f T' % B_EXT)
    print('  Analytical M_z: %.0f A/m' % M_ANALYTICAL_Z)
    print()

    # Generate test points
    test_points, labels = generate_test_points()
    print('Test points: %d (all outside the magnetic cube)' % len(test_points))

    # Hexahedral solution
    hexa_result = create_hexa_solution(test_points)

    # Tetrahedral solution
    tetra_result = create_tetra_solution(test_points)

    if tetra_result is None:
        print('Tetrahedral test skipped - Netgen not available')
        return

    # Compare results
    errors = compare_results(hexa_result, tetra_result, test_points, labels)

    # Summary
    print()
    print('=' * 70)
    print('Summary')
    print('=' * 70)

    if errors:
        avg_error = np.mean(errors)
        max_error = np.max(errors)
        min_error = np.min(errors)

        print()
        print('Mesh Comparison:')
        print('  Hexahedral: %d elements (n_div=%d)' % (
            hexa_result['n_elements'], hexa_result['n_div']))
        print('  Tetrahedral: %d elements (maxh=%.2f)' % (
            tetra_result['n_elements'], tetra_result['maxh']))

        print()
        print('Magnetization Comparison:')
        print('  Hexa M_avg_z:  %.0f A/m' % hexa_result['M_avg_z'])
        print('  Tetra M_avg_z: %.0f A/m' % tetra_result['M_avg_z'])
        print('  Analytical:    %.0f A/m' % M_ANALYTICAL_Z)

        print()
        print('Field Comparison (Tetra vs Hexa):')
        print('  Valid points:  %d / %d' % (len(errors), len(test_points)))
        print('  Average error: %.4f%%' % avg_error)
        print('  Maximum error: %.4f%%' % max_error)
        print('  Minimum error: %.4f%%' % min_error)

        if avg_error < 1.0:
            print()
            print('[PASS] Tetrahedral MSC field matches Hexahedral (< 1%%)')
            print('       Both mesh types produce accurate results')
        elif avg_error < 5.0:
            print()
            print('[PASS] Tetrahedral MSC field is GOOD (< 5%%)')
        else:
            print()
            print('[CHECK] Field comparison shows differences')

    # Save results
    output_file = os.path.join(os.path.dirname(__file__), 'solver_comparison_results.json')

    output_data = {
        'parameters': {
            'cube_size': CUBE_SIZE,
            'mu_r': MU_R,
            'chi': CHI,
            'H_ext': H_EXT,
            'M_analytical_z': M_ANALYTICAL_Z
        },
        'hexa': {
            'n_elements': hexa_result['n_elements'],
            'n_div': hexa_result['n_div'],
            't_solve': hexa_result['t_solve'],
            'iterations': hexa_result['iterations'],
            'M_avg_z': hexa_result['M_avg_z'],
            'B_values': hexa_result['B_values']
        },
        'tetra': {
            'n_elements': tetra_result['n_elements'],
            'maxh': tetra_result['maxh'],
            't_solve': tetra_result['t_solve'],
            'iterations': tetra_result['iterations'],
            'M_avg_z': tetra_result['M_avg_z'],
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
