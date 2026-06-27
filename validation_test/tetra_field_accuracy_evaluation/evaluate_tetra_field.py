#!/usr/bin/env python
"""
Tetrahedral Field Accuracy Evaluation

This script evaluates the accuracy of magnetic field computation from a
tetrahedral mesh in Radia by:

1. Running NGSolve H-formulation (perturbation potential) for a linear
   magnetic cube in uniform external field
2. Extracting the magnetization M = chi * H from the NGSolve solution
3. Transferring the magnetization to Radia tetrahedral mesh
4. Evaluating B field using rad.Fld at points OUTSIDE the magnet
5. Comparing with NGSolve B field (B = mu_0 * mu_r * H inside, mu_0 * H outside)

Test Case: Linear magnetic cube in uniform external field
- Cube size: 0.1m x 0.1m x 0.1m (matching ngsolve_cube_uniform_field.py)
- Linear material: mu_r = 100 (chi = 99)
- External field: H_ext = 1000 A/m (z-direction)

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
print('Tetrahedral Field Accuracy Evaluation')
print('=' * 70)

# Check for NGSolve
try:
    from ngsolve import *
    from ngsolve import TaskManager
    from netgen.occ import Box, Pnt, OCCGeometry, Glue
    NGSOLVE_AVAILABLE = True
except ImportError as e:
    print('NGSolve not available: %s' % e)
    NGSOLVE_AVAILABLE = False
    sys.exit(1)

from netgen_mesh_import import netgen_mesh_to_radia

# =============================================================================
# Parameters (matching ngsolve_cube_uniform_field.py)
# =============================================================================
MU_0 = 4 * np.pi * 1e-7  # Vacuum permeability [T*m/A]

# Geometry
CUBE_SIZE = 0.1      # Cube edge length [m] (matching existing script)
CUBE_HALF = 0.05     # Half size

# Material (Linear)
MU_R = 100           # Relative permeability (matching existing script)
CHI = MU_R - 1       # Magnetic susceptibility = 99

# External field
H_EXT = 1000.0       # External H field [A/m] (matching existing script)
B_EXT = MU_0 * H_EXT  # External B field [T]

# Mesh parameters (finer for better accuracy)
MAXH_CUBE = 0.015    # Mesh size in cube
MAXH_INNER = 0.03    # Mesh size in inner air
MAXH_OUTER = 0.08    # Mesh size in outer air
AIR_INNER = 0.25     # Inner air half-size
AIR_OUTER = 0.5      # Outer air half-size


def generate_test_points():
    """Generate test points for field evaluation - all outside the cube."""
    points = []
    labels = []

    # On-axis points (z-axis, outside cube: z > 0.05)
    z_values = [0.06, 0.08, 0.10, 0.15, 0.20, 0.30]
    for z in z_values:
        points.append([0.0, 0.0, z])
        labels.append('z-axis (z=%.2f)' % z)

    # Off-axis points (all outside cube)
    off_axis_points = [
        [0.06, 0.0, 0.06],
        [0.0, 0.06, 0.06],
        [0.06, 0.06, 0.06],
        [0.10, 0.0, 0.10],
        [0.0, 0.10, 0.10],
        [0.15, 0.0, 0.15],
    ]
    for pt in off_axis_points:
        points.append(pt)
        labels.append('off-axis (%.2f,%.2f,%.2f)' % tuple(pt))

    return points, labels


def solve_ngsolve_h_formulation():
    """
    Solve for H field using NGSolve H-formulation with perturbation potential.

    H-formulation:
      div(mu * grad(phi)) = 0 in all regions
      H_total = H_ext - grad(phi)
      M = chi * H (in magnetic region)

    Returns mesh, H_total, M, and B fields.
    """
    print()
    print('Setting up NGSolve H-formulation...')

    # Create geometry (matching ngsolve_cube_uniform_field.py)
    mag_cube = Box(Pnt(-CUBE_HALF, -CUBE_HALF, -CUBE_HALF),
                   Pnt(CUBE_HALF, CUBE_HALF, CUBE_HALF))
    mag_cube.mat('magnetic')
    mag_cube.maxh = MAXH_CUBE

    # Inner air (fine mesh)
    air_inner_box = Box(Pnt(-AIR_INNER, -AIR_INNER, -AIR_INNER),
                        Pnt(AIR_INNER, AIR_INNER, AIR_INNER))
    air_inner_box.maxh = MAXH_INNER

    # Outer air (coarse mesh)
    air_outer_box = Box(Pnt(-AIR_OUTER, -AIR_OUTER, -AIR_OUTER),
                        Pnt(AIR_OUTER, AIR_OUTER, AIR_OUTER))
    for face in air_outer_box.faces:
        face.name = 'outer'
    air_outer_box.maxh = MAXH_OUTER

    # Boolean operations
    air_inner = air_inner_box - mag_cube
    air_inner.mat('air_inner')

    air_outer = air_outer_box - air_inner_box
    air_outer.mat('air_outer')

    geo = Glue([air_outer, air_inner, mag_cube])

    print('Generating mesh...')
    ngmesh = OCCGeometry(geo).GenerateMesh(maxh=MAXH_OUTER, grading=0.7)
    mesh = Mesh(ngmesh)

    print('  Elements: %d' % mesh.ne)
    print('  Vertices: %d' % mesh.nv)
    print('  Materials: %s' % list(mesh.GetMaterials()))

    # H-formulation setup
    print()
    print('Assembling H-formulation...')

    n = specialcf.normal(mesh.dim)
    fes = H1(mesh, order=2)
    print('  DOFs: %d' % fes.ndof)

    u = fes.TrialFunction()
    v = fes.TestFunction()

    # Material permeability
    mu_dict = {'air_inner': MU_0, 'air_outer': MU_0, 'magnetic': MU_R * MU_0}
    mu = CoefficientFunction([mu_dict[mat] for mat in mesh.GetMaterials()])

    # Background field H_s = [0, 0, H_ext]
    Hs = CoefficientFunction((0, 0, H_EXT))
    Hsb = BoundaryFromVolumeCF(Hs)

    # Bilinear form
    a = BilinearForm(fes)
    a += mu * grad(u) * grad(v) * dx

    # Linear form (perturbation formulation)
    f = LinearForm(fes)
    f += mu * InnerProduct(grad(v), Hs) * dx
    f += -mu * v * InnerProduct(n, Hsb) * ds('outer')

    with TaskManager():
        a.Assemble()
        f.Assemble()

        # Solve
        print('Solving...')
        gfu = GridFunction(fes)
        c = Preconditioner(a, type='local')
        solvers.CG(sol=gfu.vec, rhs=f.vec, mat=a.mat, pre=c.mat,
                   tol=1e-10, printrates=False, maxsteps=10000)

        # Compute fields
        print('Computing fields...')

        # H_pert = -grad(phi)
        H_pert = -grad(gfu)

        # H_total = H_s + H_pert
        H_total = Hs + H_pert

        # Magnetization: M = chi * H (in magnetic region), 0 elsewhere
        chi_dict = {'air_inner': 0, 'air_outer': 0, 'magnetic': CHI}
        chi_cf = CoefficientFunction([chi_dict[mat] for mat in mesh.GetMaterials()])
        M_cf = chi_cf * H_total

        # B = mu_0 * (H + M) = mu_0 * H outside, mu_0 * mu_r * H inside
        B_cf = MU_0 * (H_total + M_cf)

        # Project B onto HDiv for evaluation
        fes_hdiv = HDiv(mesh, order=2)
        gfB = GridFunction(fes_hdiv)
        gfB.Set(B_cf)

        return mesh, H_total, M_cf, gfB


def extract_magnetic_elements(mesh, M_cf):
    """
    Extract tetrahedral elements from magnetic region with their magnetization.

    Uses netgen_mesh_import.extract_elements for correct mesh access.
    """
    from netgen_mesh_import import extract_elements, compute_element_centroid

    print()
    print('Extracting magnetic elements...')

    # Use centralized mesh extraction from netgen_mesh_import
    raw_elements, _ = extract_elements(mesh, material_filter='magnetic')

    elements = []
    for el_data in raw_elements:
        vertices = el_data['vertices']
        centroid = compute_element_centroid(vertices)

        # Evaluate M at centroid
        try:
            mip = mesh(centroid[0], centroid[1], centroid[2])
            Mx = float(M_cf[0](mip))
            My = float(M_cf[1](mip))
            Mz = float(M_cf[2](mip))
            magnetization = [Mx, My, Mz]
        except:
            magnetization = [0.0, 0.0, 0.0]

        elements.append({
            'vertices': vertices,
            'magnetization': magnetization
        })

    print('  Extracted %d magnetic elements' % len(elements))

    # Statistics
    M_z_values = [e['magnetization'][2] for e in elements]
    print('  M_z range: %.1f to %.1f A/m' % (min(M_z_values), max(M_z_values)))
    print('  M_z mean:  %.1f A/m' % np.mean(M_z_values))

    return elements


def create_radia_tetra_mesh(elements):
    """
    Create Radia tetrahedral mesh with element-wise magnetization.
    """
    print()
    print('Creating Radia tetrahedral mesh...')

    rad.UtiDelAll()

    polyhedra = []
    for el in elements:
        try:
            obj_id = rad.ObjTetrahedron(el['vertices'], el['magnetization'])
            polyhedra.append(obj_id)
        except Exception as e:
            pass  # Skip failed elements

    print('  Created %d Radia polyhedra' % len(polyhedra))

    container = rad.ObjCnt(polyhedra)
    return container


def evaluate_and_compare(mesh, gfB, radia_obj, test_points, labels):
    """
    Evaluate B field at test points using both NGSolve and Radia.
    """
    print()
    print('=' * 70)
    print('Field Comparison at Test Points (outside magnet)')
    print('=' * 70)

    results = []

    print()
    print('%-30s  %-12s  %-12s  %-10s' % (
        'Point', '|B| NGSolve', '|B| Radia', 'Error %'))
    print('-' * 70)

    for i, (pt, label) in enumerate(zip(test_points, labels)):
        # NGSolve B from HDiv GridFunction
        try:
            mip = mesh(pt[0], pt[1], pt[2])
            Bx_ng = float(gfB[0](mip))
            By_ng = float(gfB[1](mip))
            Bz_ng = float(gfB[2](mip))
            B_ng_mag = np.sqrt(Bx_ng**2 + By_ng**2 + Bz_ng**2)
        except:
            B_ng_mag = np.nan

        # Radia B from tetrahedral mesh
        try:
            B_radia = rad.Fld(radia_obj, 'b', pt)
            B_rad_mag = np.linalg.norm(B_radia)
        except:
            B_rad_mag = np.nan

        # Error
        if not np.isnan(B_ng_mag) and not np.isnan(B_rad_mag) and B_ng_mag > 1e-15:
            error = abs(B_rad_mag - B_ng_mag) / B_ng_mag * 100
        else:
            error = np.nan

        results.append({
            'point': pt,
            'label': label,
            'B_ngsolve_mag': B_ng_mag,
            'B_radia_mag': B_rad_mag,
            'error_percent': error
        })

        print('%-30s  %12.4e  %12.4e  %10.2f' % (
            label, B_ng_mag, B_rad_mag, error if not np.isnan(error) else 0.0))

    return results


def main():
    """Run the complete evaluation."""
    print()
    print('Parameters (matching ngsolve_cube_uniform_field.py):')
    print('  Cube size:      %.2f m' % CUBE_SIZE)
    print('  mu_r:           %d (chi = %d)' % (MU_R, CHI))
    print('  H_ext:          %.0f A/m (z-direction)' % H_EXT)
    print()

    # Generate test points
    test_points, labels = generate_test_points()
    print('Test points: %d (all outside the magnetic cube)' % len(test_points))

    # Step 1: Solve with NGSolve H-formulation
    t0 = time.time()
    mesh, H_total, M_cf, gfB = solve_ngsolve_h_formulation()
    t_ngsolve = time.time() - t0
    print('NGSolve time: %.3f s' % t_ngsolve)

    # Step 2: Extract magnetic elements with magnetization
    elements = extract_magnetic_elements(mesh, M_cf)

    # Step 3: Create Radia tetrahedral mesh
    t0 = time.time()
    radia_obj = create_radia_tetra_mesh(elements)
    t_radia = time.time() - t0
    print('Radia creation time: %.3f s' % t_radia)

    # Step 4: Compare fields
    results = evaluate_and_compare(mesh, gfB, radia_obj, test_points, labels)

    # Summary
    print()
    print('=' * 70)
    print('Summary')
    print('=' * 70)

    valid_errors = [r['error_percent'] for r in results if not np.isnan(r['error_percent'])]

    if valid_errors:
        avg_err = np.mean(valid_errors)
        max_err = np.max(valid_errors)

        print()
        print('Valid comparison points: %d / %d' % (len(valid_errors), len(results)))
        print('Average error: %.2f%%' % avg_err)
        print('Maximum error: %.2f%%' % max_err)

        if avg_err < 10.0:
            print()
            print('[PASS] Tetrahedral MSC field matches NGSolve (avg error < 10%%)')
        else:
            print()
            print('[CHECK] Field comparison shows larger errors')

    # Save results
    output_file = os.path.join(os.path.dirname(__file__), 'evaluation_results.json')
    output_data = {
        'parameters': {
            'cube_size': CUBE_SIZE,
            'mu_r': MU_R,
            'H_ext': H_EXT
        },
        'n_elements': len(elements),
        'results': results
    }
    with open(output_file, 'w') as f:
        json.dump(output_data, f, indent=2, default=float)

    print()
    print('Results saved to: %s' % output_file)
    print('=' * 70)

    return results


if __name__ == '__main__':
    results = main()
