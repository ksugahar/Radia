#!/usr/bin/env python
"""Test: Raw STEP (Cubit coords) + full model (no symmetry) + Dirichlet.

No symmetry, no Kelvin - simplest possible FEM test.
Coil coordinates match raw Cubit system.
"""
import sys, os, math, time
sys.path.insert(0, os.path.join('..','..','..','src'))
sys.path.insert(0, os.path.join('..','..','..','src','radia'))
import numpy as np
import radia as rad
from step_mesh_builder import build_mesh_from_step
from ngsolve import *

MU_0 = 4e-7 * math.pi

# Build mesh: FULL model, no symmetry, no Kelvin
print("Building mesh (full, no symmetry, no Kelvin)...")
mesh, info = build_mesh_from_step(
    'yoke.step',
    symmetry='full',
    curve_order=1,
    add_kelvin=False,
    mesh_size=0.025,
    mesh_size_air=0.05)
print(f'Materials: {set(mesh.GetMaterials())}')
print(f'Boundaries: {set(mesh.GetBoundaries())}')
print(f'ne: {mesh.ne}, nv: {mesh.nv}')

# Coil source (Cubit coordinate system)
print("Creating coil...")
rad.UtiDelAll()
from coil_builder import CoilBuilder
mm = 1e-3
coil = (CoilBuilder(current=20000.0)
    .set_start([47.5*mm, 100*mm, 0])
    .set_cross_section(width=35*mm, height=105*mm)
    .add_straight(62.5*mm).add_arc(radius=22.5*mm, arc_angle=90)
    .add_straight(50*mm).add_arc(radius=22.5*mm, arc_angle=90)
    .add_straight(62.5*mm).add_arc(radius=22.5*mm, arc_angle=90)
    .add_straight(50*mm).add_arc(radius=22.5*mm, arc_angle=90))
container = rad.ObjCnt(coil.to_radia())
H_s = rad.RadiaField(container, 'h')

# Check source field at various points
print("\nCoil source field (no yoke):")
test_pts = [
    [0, 0, 0],
    [0, -0.04, 0],     # below pole piece (gap region)
    [0.131, 0, 0],     # near main leg
    [0.0475, 0.131, 0], # coil center
]
for p in test_pts:
    B = rad.Fld(container, 'b', p)
    print(f'  ({p[0]*1e3:.0f},{p[1]*1e3:.0f},{p[2]*1e3:.0f})mm: '
          f'B=[{B[0]*1e3:.1f}, {B[1]*1e3:.1f}, {B[2]*1e3:.1f}] mT')

# Material
mu_dict = {}
for m in set(mesh.GetMaterials()):
    if 'yoke' in m.lower():
        mu_dict[m] = CF(MU_0 * 1000)
    else:
        mu_dict[m] = CF(MU_0)
mu_cf = mesh.MaterialCF(mu_dict, default=CF(MU_0))

# BC: only outer Dirichlet (no symmetry)
boundaries = set(mesh.GetBoundaries())
dirichlet_bnd = 'outer'
print(f'\nDirichlet: {dirichlet_bnd}')

fes = H1(mesh, order=1, dirichlet=dirichlet_bnd)
u, v = fes.TnT()
print(f'DOFs: {fes.ndof}')

# Assemble + solve
print("Assembling...")
t0 = time.perf_counter()
a_bf = BilinearForm(fes)
a_bf += mu_cf * grad(u) * grad(v) * dx(bonus_intorder=4)
f_lf = LinearForm(fes)
f_lf += mu_cf * H_s * grad(v) * dx(bonus_intorder=4)

with TaskManager():
    a_bf.Assemble()
    f_lf.Assemble()
print(f'Assembly: {time.perf_counter()-t0:.1f}s')

print("Solving...")
t0 = time.perf_counter()
gfu = GridFunction(fes)
with TaskManager():
    if fes.ndof < 300000:
        gfu.vec.data = a_bf.mat.Inverse(fes.FreeDofs(), inverse='pardiso') * f_lf.vec
    else:
        from ngsolve import solvers, Preconditioner
        pre = Preconditioner(a_bf, "bddc")
        a_bf.Assemble()
        solvers.BVP(bf=a_bf, lf=f_lf, gf=gfu, pre=pre, maxsteps=500, tol=1e-8)
print(f'Solve: {time.perf_counter()-t0:.1f}s')

# Post-process
H_total = H_s - grad(gfu)
B_field = mu_cf * H_total

print(f'\nMax Omega: {max(gfu.vec):.4e}')
print(f'Min Omega: {min(gfu.vec):.4e}')

# Field evaluation
# In Cubit coords: gap is at Y < -52.5mm (below the pole piece)
# Pole piece: X in [-31.25, 31.25], Y in [-52.5, 52.5], Z in [-12.5, 12.5]
# Gap center: (X=0, Y=-60, Z=0) approximately
# MSC evaluates after transform at origin, which maps to (0,0,0) in transformed coords
# The origin in Cubit coords = MSC gap center evaluation point
print('\nB field (FEM with yoke):')
eval_pts = [
    (0, 0, 0, 'origin'),
    (0, -0.055, 0, 'gap Y=-55mm'),
    (0, -0.060, 0, 'gap Y=-60mm'),
    (0, -0.065, 0, 'gap Y=-65mm'),
    (0, -0.070, 0, 'gap Y=-70mm'),
    (0, -0.080, 0, 'gap Y=-80mm'),
    (0, -0.100, 0, 'gap Y=-100mm'),
]
for px, py, pz, label in eval_pts:
    try:
        B = [float(b) for b in B_field(mesh(px, py, pz))]
        print(f'  ({px*1e3:.0f},{py*1e3:.0f},{pz*1e3:.0f})mm [{label}]: '
              f'By={B[1]*1e3:.1f} mT, |B|={np.linalg.norm(B)*1e3:.1f} mT')
    except Exception as e:
        print(f'  ({px*1e3:.0f},{py*1e3:.0f},{pz*1e3:.0f})mm [{label}]: FAILED ({e})')
