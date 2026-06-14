#!/usr/bin/env python
"""Quick test: STEP (full model) -> quarter FEM without Kelvin."""
import sys, os, math, time
sys.path.insert(0, os.path.join('..','..','..','src'))
sys.path.insert(0, os.path.join('..','..','..','src','radia'))
import numpy as np
import radia as rad
from step_mesh_builder import build_mesh_from_step
from ngsolve import *

MU_0 = 4e-7 * math.pi

# Build mesh: quarter x-z from full STEP
print("Building mesh...")
mesh, info = build_mesh_from_step('yoke.step', symmetry='quarter_xz',
                                   curve_order=1, add_kelvin=False,
                                   mesh_size=0.015)
print('Materials:', set(mesh.GetMaterials()))
print('Boundaries:', set(mesh.GetBoundaries()))
print('ne:', mesh.ne)

# Coil source
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

B_src = rad.Fld(container, 'b', [0, 0, 0])
print(f'Source B at origin: [{B_src[0]*1e3:.1f}, {B_src[1]*1e3:.1f}, {B_src[2]*1e3:.1f}] mT')

# Material
mu_dict = {}
for m in set(mesh.GetMaterials()):
    if 'yoke' in m.lower():
        mu_dict[m] = CF(MU_0 * 1000)
    else:
        mu_dict[m] = CF(MU_0)
mu_cf = mesh.MaterialCF(mu_dict, default=CF(MU_0))

# Boundary conditions
boundaries = set(mesh.GetBoundaries())
dir_parts = []
if 'outer' in boundaries: dir_parts.append('outer')
if 'sym_normal' in boundaries: dir_parts.append('sym_normal')
dirichlet_bnd = '|'.join(dir_parts)
print(f'Dirichlet: {dirichlet_bnd}')

fes = H1(mesh, order=1, dirichlet=dirichlet_bnd)
u, v = fes.TnT()
print(f'DOFs: {fes.ndof}')

# Assemble
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

# Solve
print("Solving...")
t0 = time.perf_counter()
gfu = GridFunction(fes)
with TaskManager():
    if fes.ndof < 200000:
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
print('\nField at key points:')
pts = [
    (0.001, 0, 0.001, 'near origin (on yoke surface)'),
    (0.001, -0.005, 0.001, 'gap y=-5mm'),
    (0.001, -0.010, 0.001, 'gap y=-10mm'),
    (0.005, 0, 0.01, 'x=5 z=10 y=0'),
    (0.005, -0.005, 0.01, 'x=5 z=10 y=-5'),
    (0.005, -0.010, 0.01, 'x=5 z=10 y=-10'),
    (0.010, 0, 0.025, 'x=10 z=25 y=0'),
    (0.010, -0.005, 0.025, 'x=10 z=25 y=-5'),
    (0.010, -0.010, 0.025, 'x=10 z=25 y=-10'),
]
for px, py, pz, label in pts:
    try:
        B = [float(b) for b in B_field(mesh(px, py, pz))]
        print(f'  ({px*1e3:.0f},{py*1e3:.0f},{pz*1e3:.0f})mm [{label}]: '
              f'By={B[1]*1e3:.1f} mT, |B|={np.linalg.norm(B)*1e3:.1f} mT')
    except Exception as e:
        print(f'  ({px*1e3:.0f},{py*1e3:.0f},{pz*1e3:.0f})mm [{label}]: FAILED ({e})')

# MSC reference is -976 mT in the gap (y-direction in this coord system)
print(f'\nMSC ref: By = -976 mT at gap center')
