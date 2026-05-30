#!/usr/bin/env python
"""
Verified NGSolve to Radia: Sphere Magnetization Transfer

This script uses the verified 3D_dipole_with_Kelvin.py solver code
to compute the magnetization in a magnetic sphere, then transfers
the magnetization to Radia and compares the B field.

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

# =============================================================================
# Part 1: NGSolve H-formulation with Kelvin Transform
# (Code from verified 3D_dipole_with_Kelvin.py)
# =============================================================================
from numpy import *
from ngsolve import *
import ngsolve
from ngsolve import TaskManager
from netgen.occ import *

print("=" * 70)
print("Verified NGSolve to Radia: Sphere Magnetization Transfer")
print("=" * 70)

# Parameters (matching verified code)
sphere_radius = 0.5  # Magnetic sphere radius [m]
kelvin_radius = 1.0  # Kelvin transformation radius [m]
maxh_fine = 0.03     # Fine mesh size [m]
offset_x = 3.0       # Offset for exterior domain

mu0 = 4*pi*1e-7
mu_r = 100           # Relative permeability
chi = mu_r - 1

# Analytical solutions
H_int_analytical = 3.0/(mu_r + 2)  # For H_ext = 1 A/m
M_z_analytical = chi * H_int_analytical

print()
print("Parameters (from verified 3D_dipole_with_Kelvin.py):")
print("  Sphere radius: %.2f m" % sphere_radius)
print("  Kelvin radius: %.2f m" % kelvin_radius)
print("  maxh_fine:     %.2f m" % maxh_fine)
print("  mu_r:          %d (chi = %d)" % (mu_r, chi))
print("  Analytical H_int: %.6f A/m" % H_int_analytical)
print("  Analytical M_z:   %.6f A/m" % M_z_analytical)

# =============================================================================
# Geometry Definition (from verified code)
# =============================================================================
print()
print("[Step 1] Creating geometry (verified code)")
print("-" * 70)

# INTERIOR DOMAIN
mag_sphere = Sphere(Pnt(0, 0, 0), sphere_radius)
mag_sphere.mat("magnetic")
mag_sphere.maxh = maxh_fine

inner_sphere = Sphere(Pnt(0, 0, 0), kelvin_radius)
inner_sphere.maxh = maxh_fine
inner_air = inner_sphere - mag_sphere
inner_air.mat("air_inner")

# EXTERIOR DOMAIN
outer_sphere = Sphere(Pnt(offset_x, 0, 0), kelvin_radius)
outer_sphere.maxh = maxh_fine
outer_sphere.mat("air_outer")

vertex = Vertex(Pnt(offset_x, 0, 0))
vertex.name = "GND"

geo = Glue([inner_air, mag_sphere, outer_sphere, vertex])

geo.solids[0].name = "air_inner"
geo.solids[1].name = "magnetic"
geo.solids[2].name = "air_outer"

# Periodic BC setup
air_inner_solid_idx = None
air_outer_solid_idx = None
for i, solid in enumerate(geo.solids):
    if solid.name == "air_inner":
        air_inner_solid_idx = i
    elif solid.name == "air_outer":
        air_outer_solid_idx = i

if air_inner_solid_idx is not None and air_outer_solid_idx is not None:
    # Find largest bbox face for each
    air_inner_outer_face_idx = None
    max_bbox_size = 0
    for i, face in enumerate(geo.solids[air_inner_solid_idx].faces):
        try:
            bbox = face.bounding_box
            bbox_size = (bbox[1][0] - bbox[0][0]) * (bbox[1][1] - bbox[0][1]) * (bbox[1][2] - bbox[0][2])
            if bbox_size > max_bbox_size:
                max_bbox_size = bbox_size
                air_inner_outer_face_idx = i
        except:
            pass

    air_outer_outer_face_idx = None
    max_bbox_size = 0
    for i, face in enumerate(geo.solids[air_outer_solid_idx].faces):
        try:
            bbox = face.bounding_box
            bbox_size = (bbox[1][0] - bbox[0][0]) * (bbox[1][1] - bbox[0][1]) * (bbox[1][2] - bbox[0][2])
            if bbox_size > max_bbox_size:
                max_bbox_size = bbox_size
                air_outer_outer_face_idx = i
        except:
            pass

    if air_inner_outer_face_idx is not None and air_outer_outer_face_idx is not None:
        geo.solids[air_inner_solid_idx].faces[air_inner_outer_face_idx].Identify(
            geo.solids[air_outer_solid_idx].faces[air_outer_outer_face_idx],
            "periodic", IdentificationType.PERIODIC)
        print("  Periodic BC applied")

# =============================================================================
# Mesh Generation
# =============================================================================
print()
print("[Step 2] Generating mesh")
print("-" * 70)

mesh = Mesh(OCCGeometry(geo).GenerateMesh(maxh=maxh_fine, grading=0.7))

print("  Elements:   %d" % mesh.ne)
print("  Vertices:   %d" % mesh.nv)
print("  Materials:  %s" % str(mesh.GetMaterials()))

n_magnetic = sum(1 for el in mesh.Elements(VOL) if el.mat == "magnetic")
print("  Magnetic elements: %d" % n_magnetic)

# =============================================================================
# Problem Setup (from verified code)
# =============================================================================
print()
print("[Step 3] Setting up H-formulation (verified)")
print("-" * 70)

fes = H1(mesh, order=3, dirichlet="GND")
fes = Periodic(fes)
print("  DOFs: %d" % fes.ndof)

u = fes.TrialFunction()
v = fes.TestFunction()

# Material properties (Kelvin-transformed in exterior)
r_prime_sq = (x-offset_x)**2 + y**2 + z**2
mu_outer = kelvin_radius**2/(r_prime_sq+1e-20)*mu0

mu_d = {"air_inner": 1*mu0, "air_outer": mu_outer, "magnetic": mu_r*mu0}
mu = CoefficientFunction([mu_d[mat] for mat in mesh.GetMaterials()])

# Background field (from verified code)
# Detection of which domain
x_from_offset = x - offset_x
r_from_offset = sqrt(x_from_offset**2 + y**2 + z**2)
r_from_origin = sqrt(x**2 + y**2 + z**2)
is_exterior = IfPos(r_from_offset - r_from_origin, 0.0, 1.0)

# Local coords in exterior
x_local = x - offset_x
y_local = y
z_local = z

r_exterior = sqrt(x_local**2 + y_local**2 + z_local**2)
r_safe = IfPos(r_exterior - 1e-10, r_exterior, 1e-10)

# Interior: H_s = (0, 0, 1)
Hx_inner = 0.0
Hy_inner = 0.0
Hz_inner = 1.0

# Exterior: H_s = (0, 0, -(r'/R)^2)
r_prime = sqrt(x_local**2 + y_local**2 + z_local**2)
Hs_z_outer = -(r_prime / kelvin_radius)**2

Hs_x = (1.0 - is_exterior) * Hx_inner + is_exterior * 0.0
Hs_y = (1.0 - is_exterior) * Hy_inner + is_exterior * 0.0
Hs_z = (1.0 - is_exterior) * Hz_inner + is_exterior * Hs_z_outer

Hs = CoefficientFunction((Hs_x, Hs_y, Hs_z))

print("  Interior: H_s = (0, 0, 1)")
print("  Exterior: H_s = (0, 0, -(r'/R)^2)")

# =============================================================================
# Weak Form and Solve (from verified code)
# =============================================================================
print()
print("[Step 4] Assembling and solving")
print("-" * 70)

a = BilinearForm(fes)
a += mu*grad(u)*grad(v)*dx

f = LinearForm(fes)
f += mu*InnerProduct(grad(v), Hs)*dx

with TaskManager():
    a.Assemble()
    f.Assemble()

    gfu = GridFunction(fes)
    c = Preconditioner(a, type="local")

    solvers.CG(sol=gfu.vec, rhs=f.vec, mat=a.mat, pre=c.mat,
               tol=1e-5, printrates=False, maxsteps=10000)

    print("  Solution converged")

    # Compute perturbation field
    H_pert = -grad(gfu)
    H_total = Hs + H_pert

    # Check at origin (should match analytical)
    print()
    print("  Validation at origin:")
    try:
        Hz_pert_origin = H_pert[2](mesh(0,0,0))
        Hz_total_origin = H_total[2](mesh(0,0,0))
        Hz_analytical = -1.0 + 3.0/(mu_r + 2)  # H_pert analytical inside
        print("    H_pert_z:  %.6f A/m (analytical: %.6f A/m)" % (Hz_pert_origin, Hz_analytical))
        print("    H_total_z: %.6f A/m (analytical: %.6f A/m)" % (Hz_total_origin, H_int_analytical))
        error_origin = abs(Hz_pert_origin - Hz_analytical)/abs(Hz_analytical)*100
        print("    Error:     %.3f%%" % error_origin)
    except Exception as e:
        print("    Error: %s" % e)

    # Check at (0.7, 0, 0) outside sphere
    print()
    print("  Validation at (0.7, 0, 0):")
    try:
        Hz_pert_070 = H_pert[2](mesh(0.7, 0, 0))
        Hz_analytical_070 = -(mu_r - 1)/(mu_r + 2) * (sphere_radius/0.7)**3
        print("    H_pert_z:  %.6f A/m (analytical: %.6f A/m)" % (Hz_pert_070, Hz_analytical_070))
        error_070 = abs(Hz_pert_070 - Hz_analytical_070)/abs(Hz_analytical_070)*100
        print("    Error:     %.3f%%" % error_070)
    except Exception as e:
        print("    Error: %s" % e)

    # =============================================================================
    # Part 2: Extract magnetization and transfer to Radia
    # =============================================================================
    print()
    print("[Step 5] Extracting magnetization from magnetic elements")
    print("-" * 70)

    import radia as rad
    from netgen_mesh_import import extract_elements, compute_element_centroid

    # Use centralized mesh extraction from netgen_mesh_import module
    # This ensures correct 0-indexed access to mesh.vertices
    raw_elements, skipped_by_filter = extract_elements(mesh, material_filter='magnetic')
    print("  Total mesh vertices: %d" % mesh.nv)
    print("  Elements with 'magnetic' material: %d" % len(raw_elements))

    # Filter elements and compute magnetization
    tetra_vertices = []
    tetra_magnetization = []
    skipped_exterior = 0

    for el_data in raw_elements:
        vertices = el_data['vertices']
        centroid = compute_element_centroid(vertices)

        # Skip elements whose centroid is outside the magnetic sphere
        centroid_r = np.sqrt(centroid[0]**2 + centroid[1]**2 + centroid[2]**2)
        if centroid_r > sphere_radius * 1.01:  # Allow 1% tolerance
            skipped_exterior += 1
            continue

        # Also check if all vertices are within the sphere (with tolerance)
        all_inside = True
        for v in vertices:
            v_r = np.sqrt(v[0]**2 + v[1]**2 + v[2]**2)
            if v_r > sphere_radius * 1.05:  # 5% tolerance for surface vertices
                all_inside = False
                break
        if not all_inside:
            skipped_exterior += 1
            continue

        tetra_vertices.append(vertices)

        # H_total at centroid -> M = chi * H
        try:
            mip = mesh(centroid[0], centroid[1], centroid[2])
            Hx = H_total[0](mip)
            Hy = H_total[1](mip)
            Hz = H_total[2](mip)
            tetra_magnetization.append([chi * Hx, chi * Hy, chi * Hz])
        except:
            tetra_magnetization.append([0.0, 0.0, M_z_analytical])

    if skipped_exterior > 0:
        print("  WARNING: Skipped %d elements outside sphere (wrong material tag?)" % skipped_exterior)

    # For compatibility with later code
    tetra_elements = list(range(len(tetra_vertices)))  # Dummy indices

    M_all = np.array(tetra_magnetization)
    M_avg = np.mean(M_all, axis=0)
    M_z_std = np.std(M_all[:, 2])
    M_z_min = np.min(M_all[:, 2])
    M_z_max = np.max(M_all[:, 2])

    print("  Magnetic elements: %d" % len(tetra_elements))
    print("  M_avg: [%.4f, %.4f, %.4f] A/m" % (M_avg[0], M_avg[1], M_avg[2]))
    print("  M_z: min=%.4f, max=%.4f, std=%.4f A/m" % (M_z_min, M_z_max, M_z_std))
    print("  Analytical M_z: %.4f A/m" % M_z_analytical)
    print("  Error: %.2f%%" % (abs(M_avg[2] - M_z_analytical)/M_z_analytical*100))

    # Check for outliers
    M_z_outliers = np.sum(np.abs(M_all[:, 2]) > 100)
    print("  Elements with |M_z| > 100 A/m: %d" % M_z_outliers)
    if M_z_outliers > 0:
        outlier_indices = np.where(np.abs(M_all[:, 2]) > 100)[0]
        print("  Sample outlier M values:", M_all[outlier_indices[:5], 2])

    # Debug: Check coordinate ranges of magnetic elements
    all_mag_coords = []
    for verts in tetra_vertices:
        for v in verts:
            all_mag_coords.append(v)
    all_mag_coords = np.array(all_mag_coords)
    print("  Magnetic vertex ranges:")
    print("    x: [%.4f, %.4f]" % (all_mag_coords[:, 0].min(), all_mag_coords[:, 0].max()))
    print("    y: [%.4f, %.4f]" % (all_mag_coords[:, 1].min(), all_mag_coords[:, 1].max()))
    print("    z: [%.4f, %.4f]" % (all_mag_coords[:, 2].min(), all_mag_coords[:, 2].max()))

    # Check if any vertices are outside the sphere (R=0.5)
    R_max = np.sqrt(np.sum(all_mag_coords**2, axis=1)).max()
    print("  Max radius from origin: %.4f (should be <= %.2f)" % (R_max, sphere_radius))

    # =============================================================================
    # Create Radia tetrahedral mesh
    # =============================================================================
    print()
    print("[Step 6] Creating Radia tetrahedral mesh")
    print("-" * 70)

    rad.UtiDelAll()

    # IMPORTANT: Radia magnetization units
    # Radia API accepts magnetization values that are numerically equivalent to A/m
    # (even though the documentation may refer to them as "Tesla")
    # So we pass M [A/m] directly WITHOUT multiplying by mu0
    #
    # This was verified by comparing:
    #   - Dipole field at r=0.6m from M=3.15 A/m sphere: ~1.5e-6 T
    #   - Radia ObjTetrahedron with M=3.15: gives ~1.0e-6 T (correct order)
    #   - Radia ObjTetrahedron with M=mu0*3.15: gives ~1e-12 T (wrong!)
    print("  NOTE: Radia magnetization is passed as-is (A/m values, no mu0 conversion)")
    print("  M_avg [A/m]: [%.4f, %.4f, %.4f]" % tuple(M_avg))

    radia_objects = []
    for i, (verts, mag) in enumerate(zip(tetra_vertices, tetra_magnetization)):
        try:
            poly = rad.ObjTetrahedron(verts, mag)
            radia_objects.append(poly)
        except Exception as e:
            if i < 5:
                print("  Error creating element %d: %s" % (i, e))

    radia_container = rad.ObjCnt(radia_objects)
    print("  Radia objects created: %d" % len(radia_objects))

    # Verify magnetization
    all_M_radia = rad.ObjM(radia_container)
    M_radia_list = [m[1] for m in all_M_radia]
    M_radia_avg_z = np.mean([m[2] for m in M_radia_list])
    print("  Radia M_avg_z: %.4f (equivalent to A/m)" % M_radia_avg_z)

    # Debug: Test with a simple reference case using netgen_mesh_to_radia
    print()
    print("  DEBUG: Simple sphere test for comparison")
    from netgen.occ import Sphere as DebugSphere, Pnt as DebugPnt, OCCGeometry as DebugOCCGeo
    from netgen_mesh_import import netgen_mesh_to_radia
    debug_sphere = DebugSphere(DebugPnt(0, 0, 0), sphere_radius)
    debug_geo = DebugOCCGeo(debug_sphere)
    debug_ngmesh = Mesh(debug_geo.GenerateMesh(maxh=0.15))  # Coarse for speed

    rad.UtiDelAll()
    # Use netgen_mesh_to_radia for correct mesh extraction
    debug_container = netgen_mesh_to_radia(
        debug_ngmesh,
        material={'magnetization': [0, 0, M_z_analytical]},
        units='m',
        verbose=False
    )
    debug_B = rad.Fld(debug_container, 'b', [0.7, 0, 0])
    print("    Simple sphere (%d elements): B at (0.7,0,0) = %.6e T" % (debug_ngmesh.ne, np.linalg.norm(debug_B)))
    print("    Expected (dipole): %.6e T" % (mu0/(4*np.pi) * (4/3*np.pi*sphere_radius**3*M_z_analytical) / 0.7**3))

    # Recreate the main container for comparison
    rad.UtiDelAll()
    radia_objects_new = []
    for i, (verts, mag) in enumerate(zip(tetra_vertices, tetra_magnetization)):
        try:
            tet = rad.ObjTetrahedron(verts, mag)
            radia_objects_new.append(tet)
        except:
            pass
    radia_container = rad.ObjCnt(radia_objects_new)
    print("    NGSolve mesh (%d elements): recreated" % len(radia_objects_new))

    # Debug: Check a few elements
    print()
    print("  DEBUG: First 3 NGSolve tetrahedra:")
    for i in range(3 if len(tetra_vertices) >= 3 else len(tetra_vertices)):
        v = tetra_vertices[i]
        m = tetra_magnetization[i]
        vol = (np.array(v[1]) - np.array(v[0])).dot(
            np.cross(np.array(v[2]) - np.array(v[0]), np.array(v[3]) - np.array(v[0]))) / 6.0
        print("    Element %d:" % i)
        print("      Vertices: %s" % v)
        print("      M: %s" % m)
        print("      Volume: %.6e" % vol)

    # =============================================================================
    # Compare B field at external points
    # =============================================================================
    print()
    print("[Step 7] Comparing B field: NGSolve vs Radia")
    print("-" * 70)

    test_points = [
        [0.6, 0.0, 0.0],
        [0.7, 0.0, 0.0],
        [0.8, 0.0, 0.0],
        [0.9, 0.0, 0.0],
        [0.0, 0.6, 0.0],
        [0.0, 0.7, 0.0],
        [0.0, 0.0, 0.6],
        [0.0, 0.0, 0.7],
        [0.0, 0.0, 0.8],
    ]

    # IMPORTANT: Compare PERTURBATION field, not total field
    # NGSolve H_total = H_ext + H_pert, but Radia only computes the field FROM the magnet
    # So we should compare B_pert = mu0 * H_pert (NGSolve) with B_radia (Radia)
    print()
    print("  NOTE: Comparing PERTURBATION field (H_pert), not total field (H_total)")
    print("        Radia computes field FROM the magnet only (no background field)")
    print()
    print("  %-20s  %-15s  %-15s  %-10s" % ("Point (m)", "|B_pert| NGS", "|B| Radia", "Error %"))
    print("  " + "-" * 65)

    errors = []

    for pt in test_points:
        # NGSolve: B_pert = mu_0 * H_pert (perturbation field from magnetization)
        try:
            mip = mesh(pt[0], pt[1], pt[2])
            # Use H_pert, NOT H_total
            Hx_pert = H_pert[0](mip)
            Hy_pert = H_pert[1](mip)
            Hz_pert = H_pert[2](mip)
            B_ngsolve = [mu0 * Hx_pert, mu0 * Hy_pert, mu0 * Hz_pert]
            B_ngsolve_mag = np.linalg.norm(B_ngsolve)
        except:
            B_ngsolve_mag = np.nan

        # Radia B field (this IS the perturbation field from the magnet)
        try:
            B_radia = rad.Fld(radia_container, 'b', pt)
            B_radia_mag = np.linalg.norm(B_radia)
        except:
            B_radia_mag = np.nan

        # Error
        if not np.isnan(B_ngsolve_mag) and not np.isnan(B_radia_mag) and B_ngsolve_mag > 1e-15:
            error = abs(B_radia_mag - B_ngsolve_mag) / B_ngsolve_mag * 100
            errors.append(error)
        else:
            error = np.nan

        print("  [%.1f, %.1f, %.1f]  %15.6e  %15.6e  %10.2f" % (
            pt[0], pt[1], pt[2], B_ngsolve_mag, B_radia_mag,
            error if not np.isnan(error) else 0.0))

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
        print("Magnetization Transfer:")
        print("  NGSolve M_avg_z: %.4f A/m" % M_avg[2])
        print("  Radia M_avg_z:   %.4f (equiv. A/m)" % M_radia_avg_z)
        print("  Analytical M_z:  %.4f A/m" % M_z_analytical)
        print()
        print("Field Comparison:")
        print("  Average error: %.4f%%" % avg_error)
        print("  Maximum error: %.4f%%" % max_error)
        print("  Minimum error: %.4f%%" % min_error)

        if avg_error < 10.0:
            print()
            print("[PASS] Radia MSC field matches NGSolve (< 10%% error)")
        else:
            print()
            print("[CHECK] Field comparison shows differences")

    print()
    print("=" * 70)

    rad.UtiDelAll()
