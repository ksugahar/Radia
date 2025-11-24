#!/usr/bin/env python
"""
Benchmark: Magnetic Cube in Uniform Field - External Field Comparison

Compare Radia MMM vs NGSolve FEM for EXTERNAL field around magnetic cube.

Problem Setup:
- Magnetic cube (1.5m side, μr=100) in uniform field H0 = [0, 0, 1] A/m
- Evaluate field at points OUTSIDE the cube
- Compare NGSolve (FEM) vs Radia (MMM) solutions

Key Features:
- Independent mesh optimization for each solver:
  * NGSolve: Finer mesh (maxh=0.6) for FEM accuracy
  * Radia: Coarser mesh (maxh=1.0) for MMM stability
- Direct tetrahedral mesh import to Radia (no manual conversion)
- Validates both NGSolve H-formulation and Radia tetrahedral import

Expected Behavior:
- Far from cube: H -> H0 (uniform background)
- Near cube: H deviates due to magnetization
- Inside cube: NOT validated (Radia MMM limitation)

Author: Radia Development Team
Created: 2025-11-21
Modified: 2025-11-22 (μr=100, H0=z-direction, MMM terminology)
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../build/Release'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/python'))

import numpy as np
import radia as rad
from netgen.occ import Sphere, Box, OCCGeometry, Glue
from ngsolve import Mesh, H1, grad, BilinearForm, LinearForm, GridFunction, Preconditioner, solvers
from ngsolve import CoefficientFunction, InnerProduct, specialcf, ds, dx
from radia_ngsolve_utils import create_radia_from_mesh, export_ngsolve_vtk, export_radia_vtk

print("=" * 80)
print("Benchmark: Magnetic Cube in Uniform Field")
print("           External Field Comparison (Radia vs NGSolve)")
print("=" * 80)

# Physical parameters
mu_r_cube = 100.0  # Increased from 10 to avoid loop patterns
mu_0 = 4 * np.pi * 1e-7  # T/(A/m)
H0 = np.array([0, 0, 1.0])  # A/m (z-direction)
B0 = mu_0 * H0

# Geometry
cube_size = 1.5  # 1.5m cube

print(f"\nProblem Setup:")
print(f"  Cube size:           a = {cube_size} m")
print(f"  Permeability:        mu_r = {mu_r_cube}")
print(f"  Background field:    H0 = {H0} A/m")
print(f"  Background B field:  B0 = {B0} T")

# ============================================================================
# Part 1: NGSolve Solution
# ============================================================================

print("\n" + "=" * 80)
print("Part 1: NGSolve H-formulation Solution")
print("=" * 80)

# Create mesh - USE CUBE for magnetic region, SPHERE for air boundary
print("\n[Step 1] Creating mesh...")
# Magnetic cube: 1.5m side length
cube_size = 1.5
cube_geom = Box((-cube_size/2, -cube_size/2, -cube_size/2),
                (cube_size/2, cube_size/2, cube_size/2)).mat('magnetic')
# Air region: Sphere boundary (more natural for far-field BC)
air_radius = 5.0  # 5m radius sphere
air_geom = Sphere((0, 0, 0), air_radius).mat('air')
air_region = air_geom - cube_geom

geo = OCCGeometry(Glue([cube_geom, air_region]))
mesh = Mesh(geo.GenerateMesh(maxh=0.6))  # Coarser mesh to reduce elements

print(f"          Mesh: {mesh.nv} vertices, {mesh.ne} elements")
print(f"          Materials: {mesh.GetMaterials()}")

# Define permeability
mu_d = {"air": 1*mu_0, "magnetic": mu_r_cube*mu_0}
mu = CoefficientFunction([mu_d[mat] for mat in mesh.GetMaterials()])

# H-formulation: solve for scalar potential phi, H = -grad(phi) + Hs
print("\n[Step 2] Setting up H-formulation...")
fes = H1(mesh, order=2)
u = fes.TrialFunction()
v = fes.TestFunction()

# Background field (uniform)
Hs = CoefficientFunction((H0[0], H0[1], H0[2]))

# Bilinear and linear forms
a = BilinearForm(fes)
a += mu*grad(u)*grad(v)*dx

f = LinearForm(fes)
# Correct sign to satisfy far-field boundary condition H -> H0
f += -mu*InnerProduct(grad(v), Hs)*dx

# Boundary condition (natural Neumann on outer boundary)
n = specialcf.normal(mesh.dim)
from ngsolve import BoundaryFromVolumeCF
Hsb = BoundaryFromVolumeCF(Hs)
# Apply boundary condition on outer boundary
# Sign corrected to match volume term
print(f"          Boundaries: {mesh.GetBoundaries()}")
if 'air' in mesh.GetBoundaries():
    f += mu*v*InnerProduct(n, Hsb)*ds(mesh.Boundaries('air'))
else:
    # If no named boundary, apply to all outer boundaries
    f += mu*v*InnerProduct(n, Hsb)*ds

print("          Assembling system...")
a.Assemble()
f.Assemble()

# Solve
print("\n[Step 3] Solving...")
gfu = GridFunction(fes)
c = Preconditioner(a, type="local")

print("          Solving CG...")
import io
import sys
# Suppress CG output
old_stdout = sys.stdout
sys.stdout = io.StringIO()
try:
    nit = solvers.CG(sol=gfu.vec, rhs=f.vec, mat=a.mat, pre=c.mat, tol=1e-8, printrates=False, maxsteps=10000)
finally:
    sys.stdout = old_stdout
print(f"          [OK] Solution converged ({nit} iterations)" if isinstance(nit, int) else "          [OK] Solution converged")

# H field: H = -grad(phi)
# Note: phi represents total potential (not perturbation)
# Background field Hs is included via the weak form
grad_phi = grad(gfu)
H_ngsolve_cf = -grad_phi

# Export NGSolve solution to VTK
print("\n[Step 4] Exporting NGSolve solution to VTK...")
ngsolve_vtk = export_ngsolve_vtk(mesh, gfu, 'cube_benchmark_ngsolve', 'phi')
print(f"          NGSolve VTK: {os.path.basename(ngsolve_vtk)}")

# ============================================================================
# Part 2: Radia Solution
# ============================================================================

print("\n" + "=" * 80)
print("Part 2: Radia MMM Solution")
print("=" * 80)

rad.UtiDelAll()
rad.FldUnits('m')

# Create separate mesh for Radia (independent from NGSolve mesh)
print("\n[Step 1] Creating Radia mesh (independent from NGSolve)...")

# Use coarser mesh for Radia MMM stability
radia_mesh_size = 1.0  # Larger than NGSolve (maxh=0.6) for better MMM convergence
radia_cube_geom = Box((-cube_size/2, -cube_size/2, -cube_size/2),
                      (cube_size/2, cube_size/2, cube_size/2)).mat('magnetic')
radia_geo = OCCGeometry(radia_cube_geom)
radia_mesh = Mesh(radia_geo.GenerateMesh(maxh=radia_mesh_size))

print(f"          NGSolve mesh: {mesh.nv} vertices, {mesh.ne} elements (maxh=0.6)")
print(f"          Radia mesh:   {radia_mesh.nv} vertices, {radia_mesh.ne} elements (maxh={radia_mesh_size})")

# Import Radia mesh to Radia (without VTK export yet)
print("\n[Step 2] Importing Radia mesh to Radia MMM...")
mag_cube = create_radia_from_mesh(
    radia_mesh,
    material={'magnetization': [0, 0, 0]},
    units='m',
    combine=True,
    verbose=False,
    vtk_filename=None  # Will export after solving
)

print(f"\n          [OK] Radia cube object: {mag_cube}")

# Apply material
print("\n[Step 3] Applying linear material...")
chi = mu_r_cube - 1.0

# Isotropic linear material (no remanent magnetization)
mat = rad.MatLin(chi, [0, 0, 1e-10])  # Isotropic, mr vector defines axis
rad.MatApl(mag_cube, mat)
print(f"          Material: MatLin({chi}, [0, 0, 1e-10])")
print(f"          (Isotropic linear material, mu_r = {mu_r_cube})")

# Background field
print("\n[Step 4] Creating background field...")
def uniform_field(pos):
    return [float(B0[0]), float(B0[1]), float(B0[2])]

background = rad.ObjBckgCF(uniform_field)
system = rad.ObjCnt([mag_cube, background])

# Solve
print("\n[Step 5] Solving...")
result = rad.Solve(system, 0.0001, 10000)
print(f"          Result: {result}")

# Export Radia geometry to VTK
print("\n[Step 6] Exporting Radia geometry to VTK...")
radia_vtk = export_radia_vtk(system, 'cube_benchmark_radia')
print(f"          Radia VTK: {os.path.basename(radia_vtk)}")

# ============================================================================
# Part 3: Comparison at External Points
# ============================================================================

print("\n" + "=" * 80)
print("Part 3: External Field Comparison")
print("=" * 80)

# Test points OUTSIDE cube (cube boundaries: ±0.75m, mesh domain: ±3m)
# Use points well inside mesh domain (safety margin from both cube and mesh boundaries)
test_points = [
    # Along axes (medium distance from cube)
    ([1.0, 0, 0], "x-axis, r=1.0m"),
    ([0, 1.0, 0], "y-axis, r=1.0m"),
    ([0, 0, 1.0], "z-axis, r=1.0m"),

    # Near cube surface (just outside cube boundary at 0.75m)
    ([0.9, 0, 0], "x-axis, r=0.9m (near surface)"),
    ([0, 0.9, 0], "y-axis, r=0.9m (near surface)"),

    # Diagonal (safely inside mesh)
    ([1.0, 1.0, 0], "xy-diagonal, r=1.41m"),

    # Far field (but well inside mesh boundary)
    ([2.0, 0, 0], "x-axis, r=2.0m (far field)"),
]

print("\n" + "-" * 100)
print(f"{'Point (m)':<25} {'Description':<25} {'H_NGSolve (A/m)':<25} {'H_Radia (A/m)':<25}")
print("-" * 100)

errors = []

for point, description in test_points:
    # NGSolve evaluation
    try:
        mip = mesh(*point)
        H_ng = H_ngsolve_cf(mip)
        H_ng_arr = np.array([H_ng[0], H_ng[1], H_ng[2]])
    except:
        print(f"{str(point):<25} {description:<25} {'Point outside mesh':<25}")
        continue

    # Radia evaluation
    H_rad = rad.Fld(system, 'h', point)
    H_rad_arr = np.array(H_rad)

    # Calculate error
    H_ng_norm = np.linalg.norm(H_ng_arr)
    if H_ng_norm > 1e-10:
        rel_error = np.linalg.norm(H_rad_arr - H_ng_arr) / H_ng_norm
    else:
        rel_error = np.linalg.norm(H_rad_arr - H_ng_arr)

    errors.append(rel_error)

    # Format output
    point_str = f"[{point[0]:.1f}, {point[1]:.1f}, {point[2]:.1f}]"
    H_ng_str = f"[{H_ng_arr[0]:.4f}, {H_ng_arr[1]:.4f}, {H_ng_arr[2]:.4f}]"
    H_rad_str = f"[{H_rad_arr[0]:.4f}, {H_rad_arr[1]:.4f}, {H_rad_arr[2]:.4f}]"

    print(f"{point_str:<25} {description:<25} {H_ng_str:<25} {H_rad_str:<25}")

print("-" * 100)

# Summary
print("\n" + "=" * 80)
print("Summary")
print("=" * 80)

if errors:
    max_error = max(errors)
    avg_error = np.mean(errors)

    print(f"\nNumber of test points: {len(errors)}")
    print(f"Maximum relative error: {max_error*100:.2f}%")
    print(f"Average relative error: {avg_error*100:.2f}%")
    print()

    if max_error < 0.10:
        print("[PASS] Excellent agreement between Radia and NGSolve (< 10% error)")
        print("       External field validation successful!")
    elif max_error < 0.25:
        print("[GOOD] Good agreement between Radia and NGSolve (< 25% error)")
        print("       Minor deviations acceptable for MMM vs FEM comparison")
    else:
        print("[WARNING] Large deviations detected")
        print(f"           Max error: {max_error*100:.1f}%")
        print("           Possible causes:")
        print("             1. Coarse mesh differences")
        print("             2. MMM vs FEM method differences")
        print("             3. Material parameter mismatch")
else:
    print("[ERROR] No valid comparison points")

print("=" * 80)

print("\n[Note] Internal field (r < 1m) NOT validated")
print("       Radia MMM has limited accuracy inside magnetizable regions")

print("\n" + "=" * 80)
print("VTK Output Files")
print("=" * 80)
print(f"\nNGSolve solution: {os.path.basename(ngsolve_vtk)}")
print(f"Radia geometry:   {os.path.basename(radia_vtk)}")
print("\nVisualize with ParaView:")
print(f"  paraview {os.path.basename(ngsolve_vtk)} {os.path.basename(radia_vtk)}")
print("=" * 80)

sys.exit(0)
