#!/usr/bin/env python
"""
NGSolve Reference Solution for Uniformly Magnetized Cube

This script computes the magnetic field from a uniformly magnetized cube
using NGSolve's A-formulation (vector potential).

The problem:
- Cube size: 1.0m x 1.0m x 1.0m centered at origin
- Magnetization: M = [0, 0, 1.0e6] A/m (uniform, permanent magnet)
- Far field boundary: A x n = 0 (zero tangential A)

Uses the A-formulation with HCurl space:
  curl(curl(A)) = mu_0 * curl(M)
  B = curl(A)

Author: Radia Development Team
Date: 2025-12-13
"""
import os
import sys
import json
import time
import numpy as np

print('=' * 70)
print('NGSolve Reference Solution: Uniformly Magnetized Cube')
print('=' * 70)

try:
    from ngsolve import *
    from ngsolve import TaskManager
    from netgen.occ import Box, Pnt, OCCGeometry, Glue
    NGSOLVE_AVAILABLE = True
except ImportError as e:
    print('NGSolve not available: %s' % e)
    NGSOLVE_AVAILABLE = False

# =============================================================================
# Parameters
# =============================================================================
MU_0 = 4 * np.pi * 1e-7  # Vacuum permeability [T*m/A]

# Geometry
CUBE_SIZE = 1.0      # Cube edge length [m]
CUBE_HALF = 0.5      # Half size
AIR_SIZE = 12.0      # Air domain size [m] - larger to include all test points
AIR_HALF = 6.0       # Air domain half size

# Magnetization
M_Z = 1.0e6          # Magnetization magnitude [A/m]

# Mesh
MAXH_CUBE = 0.12     # Fine mesh size in cube
MAXH_AIR = 0.8       # Coarse mesh size in air


def generate_test_points():
    """Generate test points for field evaluation - all outside the cube."""
    points = []
    labels = []

    # On-axis points (z-axis, outside cube: z > 0.5)
    z_values = [0.6, 0.8, 1.0, 1.5, 2.0, 3.0, 5.0]
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
        [0.0, 2.0, 2.0],
    ]
    for pt in off_axis_points:
        points.append(pt)
        labels.append('off-axis (%.1f,%.1f,%.1f)' % tuple(pt))

    return points, labels


def compute_ngsolve_solution(test_points):
    """
    Compute magnetic field using NGSolve A-formulation.

    For a permanent magnet with uniform M:
      curl(1/mu_0 * curl(A)) = curl(M)

    Weak form:
      integral(curl(A) . curl(v)) = mu_0 * integral(M . curl(v))

    with B = curl(A).
    """
    if not NGSOLVE_AVAILABLE:
        print('NGSolve not available')
        return None

    print()
    print('Creating geometry...')

    # Create magnetic cube
    mag_cube = Box(Pnt(-CUBE_HALF, -CUBE_HALF, -CUBE_HALF),
                   Pnt(CUBE_HALF, CUBE_HALF, CUBE_HALF))
    mag_cube.mat('magnetic')
    mag_cube.maxh = MAXH_CUBE

    # Create air domain (larger than all test points)
    air_box = Box(Pnt(-AIR_HALF, -AIR_HALF, -AIR_HALF),
                  Pnt(AIR_HALF, AIR_HALF, AIR_HALF))

    # Name outer boundary before boolean operation
    for face in air_box.faces:
        face.name = 'outer'
    air_box.maxh = MAXH_AIR

    # Create air region (air - cube)
    air_region = air_box - mag_cube
    air_region.mat('air')

    # Combine
    geo = Glue([mag_cube, air_region])

    print('Generating mesh...')
    ngmesh = OCCGeometry(geo).GenerateMesh(maxh=MAXH_AIR, grading=0.4)
    mesh = Mesh(ngmesh)

    print('  Elements: %d' % mesh.ne)
    print('  Vertices: %d' % mesh.nv)
    print('  Materials: %s' % list(mesh.GetMaterials()))

    # ==========================================================================
    # A-Formulation (Vector Potential) with HCurl space
    # ==========================================================================
    print()
    print('Setting up A-formulation (HCurl)...')

    # HCurl space for vector potential A, with Dirichlet BC on outer boundary
    fes = HCurl(mesh, order=2, dirichlet='outer')
    print('  DOFs: %d' % fes.ndof)

    u = fes.TrialFunction()
    v = fes.TestFunction()

    # Magnetization: M = [0, 0, M_z] in magnetic region, 0 in air
    M_cf = mesh.MaterialCF({'magnetic': (0, 0, M_Z)}, default=(0, 0, 0))

    # Bilinear form: a(u,v) = integral(curl(u) . curl(v))
    a = BilinearForm(fes)
    a += curl(u) * curl(v) * dx

    # Linear form: f(v) = mu_0 * integral(M . curl(v))
    f = LinearForm(fes)
    f += MU_0 * M_cf * curl(v) * dx

    print('Assembling...')
    with TaskManager():
        a.Assemble()
        f.Assemble()

        print('Solving...')
        gfA = GridFunction(fes)

        # Use direct solver for robustness
        gfA.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse='sparsecholesky') * f.vec

        print('Computing B = curl(A)...')
        # B = curl(A) as a CoefficientFunction
        B_cf = curl(gfA)

        # Project B onto HDiv space for better evaluation
        print('Projecting B onto HDiv space...')
        fes_hdiv = HDiv(mesh, order=2)
        gfB = GridFunction(fes_hdiv)
        gfB.Set(B_cf)

        # ==========================================================================
        # Evaluate at Test Points
        # ==========================================================================
        print()
        print('Evaluating field at test points...')

        B_values = []
        for pt in test_points:
            try:
                mip = mesh(pt[0], pt[1], pt[2])
                # Use HDiv GridFunction for evaluation
                Bx = gfB[0](mip)
                By = gfB[1](mip)
                Bz = gfB[2](mip)
                B_values.append([float(Bx), float(By), float(Bz)])
                B_mag = np.sqrt(Bx**2 + By**2 + Bz**2)
                print('  [%.2f, %.2f, %.2f]: |B| = %.6e T' %
                      (pt[0], pt[1], pt[2], B_mag))
            except Exception as e:
                print('  Error at [%.2f, %.2f, %.2f]: %s' % (pt[0], pt[1], pt[2], e))
                B_values.append([np.nan, np.nan, np.nan])

        return {
            'method': 'NGSolve A-formulation (HCurl->HDiv)',
            'mesh_ne': mesh.ne,
            'mesh_nv': mesh.nv,
            'ndof': fes.ndof,
            'maxh_cube': MAXH_CUBE,
            'maxh_air': MAXH_AIR,
            'B_values': B_values
        }


def main():
    """Run NGSolve reference computation."""
    print()
    print('Parameters:')
    print('  Cube size:      %.1f m' % CUBE_SIZE)
    print('  Magnetization:  [0, 0, %.0e] A/m' % M_Z)
    print('  Air domain:     %.1f m' % AIR_SIZE)
    print()

    # Generate test points
    test_points, labels = generate_test_points()

    # Compute NGSolve solution
    t0 = time.time()
    result = compute_ngsolve_solution(test_points)
    t_elapsed = time.time() - t0

    if result is None:
        print('NGSolve computation failed')
        return None

    result['time'] = t_elapsed

    # ==========================================================================
    # Display Results
    # ==========================================================================
    print()
    print('=' * 70)
    print('NGSolve Reference Results')
    print('=' * 70)

    print()
    print('Computation time: %.3f s' % t_elapsed)
    print()

    print('%s  %s  %s  %s  %s' % (
        'Point'.center(25),
        'Bx (T)'.center(12),
        'By (T)'.center(12),
        'Bz (T)'.center(12),
        '|B| (T)'.center(12)
    ))
    print('-' * 75)

    for i, (pt, label) in enumerate(zip(test_points, labels)):
        B = np.array(result['B_values'][i])
        B_mag = np.linalg.norm(B)
        print('%s  %12.4e  %12.4e  %12.4e  %12.4e' % (
            label[:25].ljust(25),
            B[0], B[1], B[2], B_mag
        ))

    # ==========================================================================
    # Save Results
    # ==========================================================================
    output_file = os.path.join(os.path.dirname(__file__), 'ngsolve_reference_results.json')

    output_data = {
        'parameters': {
            'cube_size': CUBE_SIZE,
            'magnetization': [0, 0, M_Z],
            'air_size': AIR_SIZE,
            'mu_0': MU_0
        },
        'test_points': test_points,
        'labels': labels,
        'ngsolve': {
            'method': result['method'],
            'mesh_ne': result['mesh_ne'],
            'mesh_nv': result['mesh_nv'],
            'ndof': result['ndof'],
            'maxh_cube': result['maxh_cube'],
            'maxh_air': result['maxh_air'],
            'time': result['time'],
            'B_values': result['B_values']
        }
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
