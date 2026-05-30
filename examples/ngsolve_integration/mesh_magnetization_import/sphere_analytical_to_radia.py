#!/usr/bin/env python
"""
Sphere Analytical Magnetization to Radia Field Comparison

This script tests if Radia's MSC (Magnetic Surface Charge) method correctly
computes the B field from a uniformly magnetized sphere by:

1. Creating a tetrahedral mesh of a sphere using Netgen
2. Assigning analytical uniform magnetization M = [0, 0, M_z] to each element
3. Computing B field at external points using Radia
4. Comparing with analytical solution (dipole field)

The analytical solution for B field outside a uniformly magnetized sphere:
  B = (mu_0/4pi) * (3(m*r)r/r^5 - m/r^3)
where m = (4/3)*pi*a^3*M is the magnetic dipole moment.

On the z-axis (r = [0, 0, z], z > a):
  B_z = (mu_0/4pi) * (2*m)/(z^3) = (mu_0/4pi) * (8/3)*pi*a^3*M/z^3
      = (2/3) * mu_0 * (a/z)^3 * M

IMPORTANT: Uses netgen_mesh_import module for correct mesh transfer.

Author: Radia Development Team
Date: 2025-12-13
"""
import os
import sys
import json
import numpy as np

# Path setup
_script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_script_dir, '..', '..', '..', 'build', 'Release'))
sys.path.insert(0, os.path.join(_script_dir, '..', '..', '..', 'src', 'radia'))

os.chdir(_script_dir)

import radia as rad

# Netgen imports
from netgen.occ import Sphere, Pnt, OCCGeometry
from ngsolve import Mesh
from ngsolve import TaskManager

# Radia mesh import - use the module for correct mesh transfer!
from netgen_mesh_import import netgen_mesh_to_radia

print("=" * 70)
print("Sphere: Analytical Magnetization to Radia Field Comparison")
print("=" * 70)

# =============================================================================
# Parameters
# =============================================================================
MU_0 = 4 * np.pi * 1e-7  # Vacuum permeability

sphere_radius = 0.5      # Sphere radius [m]
maxh = 0.08              # Mesh size [m]
M_z = 1000.0             # Uniform magnetization [A/m]

# Dipole moment
V_sphere = (4.0/3.0) * np.pi * sphere_radius**3
m_z = V_sphere * M_z  # A*m^2

print()
print("Parameters:")
print("  Sphere radius: %.2f m" % sphere_radius)
print("  Magnetization: M = [0, 0, %.1f] A/m" % M_z)
print("  Volume:        %.6f m^3" % V_sphere)
print("  Dipole moment: m = %.6f A*m^2" % m_z)

# =============================================================================
# Analytical B field outside uniformly magnetized sphere
# =============================================================================
def B_analytical(r, M_z, a):
    """
    Analytical B field outside a uniformly magnetized sphere.

    For a sphere of radius a with uniform magnetization M = [0, 0, M_z]:
    The field outside is a dipole field.

    B(r) = (mu_0 / 4*pi) * (3*(m*r_hat)*r_hat - m) / r^3

    where m = (4/3)*pi*a^3 * M is the magnetic moment.
    """
    x, y, z = r
    r_mag = np.sqrt(x**2 + y**2 + z**2)

    if r_mag < a:
        # Inside sphere: B = (2/3)*mu_0*M
        return [0.0, 0.0, (2.0/3.0) * MU_0 * M_z]

    # Outside: dipole field
    # m = (4/3)*pi*a^3*M (pointing in z)
    m = (4.0/3.0) * np.pi * a**3 * M_z

    # B = (mu_0/4pi) * (3(m*r)r/r^5 - m/r^3)
    # With m along z: m*r = m*z
    m_dot_r = m * z
    factor = MU_0 / (4.0 * np.pi)

    Bx = factor * 3 * m_dot_r * x / r_mag**5
    By = factor * 3 * m_dot_r * y / r_mag**5
    Bz = factor * (3 * m_dot_r * z / r_mag**5 - m / r_mag**3)

    return [Bx, By, Bz]

# =============================================================================
# Create tetrahedral mesh with Netgen
# =============================================================================
print()
print("[Step 1] Creating tetrahedral mesh")
print("-" * 70)

sphere = Sphere(Pnt(0, 0, 0), sphere_radius)
sphere.mat("magnetic")
geo = OCCGeometry(sphere)
with TaskManager():
    ngmesh = geo.GenerateMesh(maxh=maxh)
    mesh = Mesh(ngmesh)

    print("  Elements: %d" % mesh.ne)
    print("  Vertices: %d" % mesh.nv)

    # =============================================================================
    # Create Radia tetrahedral mesh using netgen_mesh_import module
    # =============================================================================
    print()
    print("[Step 2] Creating Radia tetrahedral mesh")
    print("-" * 70)

    rad.UtiDelAll()

    # Use netgen_mesh_import for correct mesh transfer
    radia_container = netgen_mesh_to_radia(
        mesh,
        material={'magnetization': [0, 0, M_z]},
        units='m',
        material_filter='magnetic',
        verbose=True
    )

    # Verify magnetization
    all_M_radia = rad.ObjM(radia_container)
    M_radia_list = [m[1] for m in all_M_radia]
    M_radia_avg_z = np.mean([m[2] for m in M_radia_list])
    print("  Radia M_avg_z: %.2f A/m (expected: %.2f A/m)" % (M_radia_avg_z, M_z))

    # =============================================================================
    # Compare B field at external points
    # =============================================================================
    print()
    print("[Step 3] Comparing B field: Analytical vs Radia")
    print("-" * 70)

    test_points = [
        # On z-axis (outside sphere)
        [0.0, 0.0, 0.6],
        [0.0, 0.0, 0.7],
        [0.0, 0.0, 0.8],
        [0.0, 0.0, 1.0],
        [0.0, 0.0, 1.5],
        [0.0, 0.0, 2.0],
        # On x-axis
        [0.6, 0.0, 0.0],
        [0.8, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        # Off-axis
        [0.5, 0.5, 0.5],
        [0.6, 0.6, 0.6],
        [1.0, 1.0, 1.0],
    ]

    print()
    print("  %-20s  %-15s  %-15s  %-10s" % ("Point (m)", "|B| Analytical", "|B| Radia", "Error %"))
    print("  " + "-" * 65)

    errors = []
    results = []

    for pt in test_points:
        # Analytical B field
        B_anal = B_analytical(pt, M_z, sphere_radius)
        B_anal_mag = np.linalg.norm(B_anal)

        # Radia B field
        try:
            B_radia = rad.Fld(radia_container, 'b', pt)
            B_radia_mag = np.linalg.norm(B_radia)
        except:
            B_radia = [np.nan, np.nan, np.nan]
            B_radia_mag = np.nan

        # Error
        if B_anal_mag > 1e-15 and not np.isnan(B_radia_mag):
            error = abs(B_radia_mag - B_anal_mag) / B_anal_mag * 100
            errors.append(error)
        else:
            error = np.nan

        print("  [%.1f, %.1f, %.1f]  %15.6e  %15.6e  %10.2f" % (
            pt[0], pt[1], pt[2], B_anal_mag, B_radia_mag,
            error if not np.isnan(error) else 0.0))

        results.append({
            'point': pt,
            'B_analytical': list(B_anal),
            'B_radia': list(B_radia) if not np.isnan(B_radia_mag) else None,
            'error_percent': error if not np.isnan(error) else None
        })

    # =============================================================================
    # Summary
    # =============================================================================
    print()
    print("=" * 70)
    print("Summary")
    print("=" * 70)

    if errors:
        avg_error = np.mean(errors)
        max_error = np.max(errors)
        min_error = np.min(errors)

        print()
        print("Mesh:")
        print("  Tetrahedral elements: %d" % mesh.ne)
        print("  maxh:                 %.2f m" % maxh)
        print()
        print("Magnetization:")
        print("  Uniform M_z: %.1f A/m" % M_z)
        print("  Radia M_avg_z: %.1f A/m" % M_radia_avg_z)
        print()
        print("Field Comparison (Analytical dipole vs Radia MSC):")
        print("  Test points:   %d" % len(test_points))
        print("  Valid points:  %d" % len(errors))
        print("  Average error: %.4f%%" % avg_error)
        print("  Maximum error: %.4f%%" % max_error)
        print("  Minimum error: %.4f%%" % min_error)

        if avg_error < 5.0:
            print()
            print("[PASS] Radia MSC matches analytical dipole field (< 5%% error)")
        elif avg_error < 10.0:
            print()
            print("[GOOD] Radia MSC is acceptable (< 10%% error)")
        else:
            print()
            print("[CHECK] Radia MSC shows differences from analytical")

    # Save results
    output_data = {
        'parameters': {
            'sphere_radius': sphere_radius,
            'maxh': maxh,
            'M_z': M_z,
            'V_sphere': V_sphere,
            'm_z': m_z
        },
        'mesh': {
            'n_elements': mesh.ne,
            'n_vertices': mesh.nv
        },
        'field_comparison': {
            'test_points': len(test_points),
            'valid_points': len(errors),
            'avg_error_percent': float(avg_error) if errors else None,
            'max_error_percent': float(max_error) if errors else None,
            'min_error_percent': float(min_error) if errors else None
        },
        'results': results
    }

    output_file = 'sphere_analytical_to_radia_results.json'
    with open(output_file, 'w') as f:
        json.dump(output_data, f, indent=2)

    print()
    print("Results saved to: %s" % output_file)
    print()
    print("=" * 70)

    rad.UtiDelAll()
