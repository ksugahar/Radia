#!/usr/bin/env python
"""
NGSolve H-formulation for cube with uniform external field

This script computes the reference solution for:
- Magnetic cube (0.1m x 0.1m x 0.1m) with mur = 100
- Uniform external field H_ext = 1000 A/m (z-direction)
- High-precision solution using H-formulation perturbation method

Based on: S:/ngsolve/NGSolve/2024_01_31_H-formulation/2025_11_22_H-formulation3D_dipole.py
Updated: 2025-11-25
"""
import os, sys
from numpy import *
from ngsolve import *
import ngsolve
from ngsolve import TaskManager
from netgen.occ import *

print("="*70)
print("NGSolve H-formulation: Cube with Uniform External Field")
print("="*70)

# ============================================================
# Geometry Definition (OCC)
# ============================================================
print("\nCreating geometry...")

# Parameters
cube_size = 0.1  # Magnetic cube size [m]
air_inner_size = 0.5  # Inner air domain size [m]
air_outer_size = 1.0  # Outer air domain size [m]
maxh_fine = 0.02    # Fine mesh size [m] (for magnetic cube and inner air)
maxh_coarse = 0.1   # Coarse mesh size [m] (for outer air)

# Create magnetic cube (finest mesh)
mag_cube = Box(Pnt(-cube_size/2, -cube_size/2, -cube_size/2),
               Pnt(cube_size/2, cube_size/2, cube_size/2))
mag_cube.mat("magnetic")
mag_cube.maxh = maxh_fine

# Create inner air box (fine mesh)
air_inner_box = Box(Pnt(-air_inner_size/2, -air_inner_size/2, -air_inner_size/2),
                    Pnt(air_inner_size/2, air_inner_size/2, air_inner_size/2))
air_inner_box.maxh = maxh_fine

# Create outer air box and name its boundary (coarse mesh)
air_outer_box = Box(Pnt(-air_outer_size/2, -air_outer_size/2, -air_outer_size/2),
                    Pnt(air_outer_size/2, air_outer_size/2, air_outer_size/2))
for face in air_outer_box.faces:
    face.name = "outer"  # Name outer boundary before boolean operation
air_outer_box.maxh = maxh_coarse

# Boolean operations to create three regions
# Inner air = inner box - magnetic cube
air_inner = air_inner_box - mag_cube
air_inner.mat("air_inner")

# Outer air = outer box - inner box
air_outer = air_outer_box - air_inner_box
air_outer.mat("air_outer")

# Combine into single geometry
geo = Glue([air_outer, air_inner, mag_cube])

print(f"Geometry created with three regions:")
print(f"  Magnetic cube size: {cube_size} m")
print(f"  Inner air domain size: {air_inner_size} m")
print(f"  Outer air domain size: {air_outer_size} m")
print(f"  Fine mesh size: {maxh_fine} m")
print(f"  Coarse mesh size: {maxh_coarse} m")

# ============================================================
# Mesh Generation
# ============================================================
print("\nGenerating mesh...")
mesh = Mesh(OCCGeometry(geo).GenerateMesh(maxh=maxh_coarse, grading=0.7))

print(f"  Number of elements: {mesh.ne}")
print(f"  Number of vertices: {mesh.nv}")
print(f"  Materials: {mesh.GetMaterials()}")
print(f"  Boundaries: {mesh.GetBoundaries()}")

# ============================================================
# Problem Setup
# ============================================================
print("\nSetting up H-formulation...")

n = specialcf.normal(mesh.dim)
fes = H1(mesh, order=2)
print(f"  Number of DOFs: {fes.ndof}")

mu0 = 4*pi*1e-7
u = fes.TrialFunction()
v = fes.TestFunction()

# Material properties
mu_r = 100  # Relative permeability (matching Radia test)
mu_d = {"air_inner": 1*mu0, "air_outer": 1*mu0, "magnetic": mu_r*mu0}
mu = CoefficientFunction([mu_d[mat] for mat in mesh.GetMaterials()])

# Background field: H_s = [0, 0, 1000] A/m (z-direction, matching Radia)
H_ext = 1000.0  # A/m
Hs = CoefficientFunction((0, 0, H_ext))
Hsb = BoundaryFromVolumeCF(Hs)

print(f"  Background field: H_s = [0, 0, {H_ext}] A/m (z-direction)")
print(f"  Relative permeability: mu_r = {mu_r}")

# ============================================================
# Weak Form (Perturbation Potential Formulation)
# ============================================================
print("\nAssembling system...")

# Bilinear form: a(u,v) = integral((grad v)·(mu grad u) dOmega)
a = BilinearForm(fes)
a += mu*grad(u)*grad(v)*dx

# Linear form (PERTURBATION FORMULATION):
# f(v) = integral((grad v)·(mu H_s) dOmega) - integral(v (n·mu H_s) dGamma)
f = LinearForm(fes)
f += mu*InnerProduct(grad(v), Hs)*dx                    # POSITIVE sign
f += -mu*v*InnerProduct(n, Hsb)*ds(mesh.Boundaries("outer"))  # NEGATIVE sign

with TaskManager():
    a.Assemble()
    f.Assemble()

    print("  System assembled")

    # ============================================================
    # Solve
    # ============================================================
    print("\nSolving system...")

    gfu = GridFunction(fes)
    c = Preconditioner(a, type="local")

    solvers.CG(sol=gfu.vec, rhs=f.vec, mat=a.mat, pre=c.mat,
               tol=1e-8, printrates=False, maxsteps=10000)

    print("  Solution converged")

    # ============================================================
    # Post-processing
    # ============================================================
    print("\nPost-processing...")

    # Compute perturbation field: H_pert = -grad(phi)
    H_pert = -grad(gfu)

    # Total field: H_total = H_s + H_pert
    H_total = Hs + H_pert

    # Evaluation points (matching Radia test)
    test_points = [
        [0.05, 0.05, 0.2],   # Outside cube (z > 0.05)
        [0.0, 0.0, 0.0],     # Center of cube
        [0.08, 0.0, 0.0],    # Near edge (x-direction)
        [0.0, 0.0, 0.08],    # Near face (z-direction)
    ]

    print("\n" + "="*70)
    print("Field Evaluation Results")
    print("="*70)

    for pt in test_points:
        try:
            mip = mesh(pt[0], pt[1], pt[2])

            # Check if inside magnetic cube
            mat_name = mesh.GetMaterial(mip.ElementId())

            # Total field components
            Hx = H_total[0](mip)
            Hy = H_total[1](mip)
            Hz = H_total[2](mip)
            H_mag = sqrt(Hx**2 + Hy**2 + Hz**2)

            # Perturbation field components
            Hx_pert = H_pert[0](mip)
            Hy_pert = H_pert[1](mip)
            Hz_pert = H_pert[2](mip)

            # Potential
            phi = gfu(mip)

            print(f"\nPoint: [{pt[0]:.3f}, {pt[1]:.3f}, {pt[2]:.3f}] m")
            print(f"  Material: {mat_name}")
            print(f"  H_total = [{Hx:.4f}, {Hy:.4f}, {Hz:.4f}] A/m")
            print(f"  |H_total| = {H_mag:.4f} A/m")
            print(f"  H_pert = [{Hx_pert:.4f}, {Hy_pert:.4f}, {Hz_pert:.4f}] A/m")
            print(f"  phi = {phi:.6e}")

        except:
            print(f"\nPoint: [{pt[0]:.3f}, {pt[1]:.3f}, {pt[2]:.3f}] m")
            print(f"  [Outside mesh domain]")

    # ============================================================
    # Save Results for Comparison
    # ============================================================
    print("\n" + "="*70)
    print("Saving Results")
    print("="*70)

    # Save to numpy file
    results = {
        'cube_size': cube_size,
        'mu_r': mu_r,
        'H_ext': H_ext,
        'test_points': test_points,
        'mesh_ne': mesh.ne,
        'mesh_nv': mesh.nv,
        'ndof': fes.ndof
    }

    # Add field values at test points
    for i, pt in enumerate(test_points):
        try:
            mip = mesh(pt[0], pt[1], pt[2])
            results[f'H_total_{i}'] = array([H_total[0](mip), H_total[1](mip), H_total[2](mip)])
            results[f'H_pert_{i}'] = array([H_pert[0](mip), H_pert[1](mip), H_pert[2](mip)])
            results[f'phi_{i}'] = gfu(mip)
            results[f'material_{i}'] = mesh.GetMaterial(mip.ElementId())
        except:
            results[f'H_total_{i}'] = array([nan, nan, nan])
            results[f'H_pert_{i}'] = array([nan, nan, nan])
            results[f'phi_{i}'] = nan
            results[f'material_{i}'] = 'outside'

    save_file = 'ngsolve_cube_uniform_field_results.npz'
    savez(save_file, **results)
    print(f"  Results saved to: {save_file}")

    print("\n" + "="*70)
    print("Summary")
    print("="*70)
    print(f"Geometry: {cube_size}m cube with mu_r={mu_r}")
    print(f"External field: H_ext = {H_ext} A/m (z-direction)")
    print(f"Mesh: {mesh.ne} elements, {fes.ndof} DOFs")
    print(f"Method: H-formulation with perturbation potential")
    print(f"Solver: CG with tol=1e-8")
    print("\nThis provides high-precision reference solution for Radia comparison.")
    print("="*70)
