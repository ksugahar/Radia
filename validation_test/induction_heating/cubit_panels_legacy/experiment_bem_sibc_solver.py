"""Quick test of ScalarBIESIBCSolver on sphere (analytical validation)."""
import math
import os
import sys
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../src/radia'))

MU_0 = 4e-7 * np.pi

from ngsolve import Mesh, H1, ds, BND, TaskManager, z
from netgen.occ import Sphere, Pnt, OCCGeometry, Glue
from radia.bem_sibc_solver import ScalarBIESIBCSolver

R = 0.01
B0 = 0.001
freq = 1000
omega = 2 * math.pi * freq
H0 = B0 / MU_0
beta = 1j * omega * MU_0 * R

print(f"Sphere R={R*1e3}mm, B0={B0*1e3}mT, f={freq}Hz, |beta|={abs(beta):.4e}")

sph = Sphere(Pnt(0, 0, 0), R)
for f in sph.faces:
    f.name = "sphere"
geo = OCCGeometry(Glue(sph.faces))
mesh = Mesh(geo.GenerateMesh(maxh=R / 5))
with TaskManager():
    mesh.Curve(3)
    print(f"Mesh: {mesh.GetNE(BND)} elements")

    solver = ScalarBIESIBCSolver(mesh, order=1)
    print(f"Assembled: {solver.ndof} DOFs in {solver.t_assembly:.1f}s")

    # phi_inc = -H0 * z
    phi_inc_cf = -H0 * z

    print(f"\n{'Zs/|beta|':>10s} {'Ana':>10s} {'BIE':>10s} {'error%':>8s}")
    print("-" * 42)

    for ratio in [0.0, 0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0]:
        Zs = ratio * abs(beta) * (1 + 1j) / math.sqrt(2) if ratio > 0 else 0

        J_max_ana = abs(1.5 * 1j * omega * B0 * R / (beta + Zs))
        H_ana = J_max_ana * math.sqrt(2.0 / 3.0)

        result = solver.solve(phi_inc_cf, Z_s=Zs, omega=omega)
        H_bie = result['H_t_rms']

        err = (H_bie / H_ana - 1) * 100 if H_ana > 0 else 0
        print(f"{ratio:10.3f} {H_ana:10.2f} {H_bie:10.2f} {err:+8.2f}%")

    print("\nAll errors should be < 0.1% (mesh discretization).")
