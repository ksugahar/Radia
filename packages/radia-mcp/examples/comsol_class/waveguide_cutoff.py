"""COMSOL-class #53 -- rectangular waveguide cutoff frequencies (2D Helmholtz eigenmodes).

The wave-physics entry in the series: a hollow metallic waveguide's modes solve the transverse
Helmholtz eigenproblem  -nabla_t^2 psi = k_c^2 psi  on the cross-section, with TE = NEUMANN
(dH_z/dn=0) and TM = DIRICHLET (E_z=0) walls. The squared cutoff wavenumbers k_c^2 are the
Laplacian eigenvalues; f_c = c k_c/(2 pi). For a rectangle a x b the spectrum is EXACT:

    f_c,mn = (c/2) sqrt((m/a)^2 + (n/b)^2),   TE: m,n>=0 (not 0,0); TM: m,n>=1.

WR-90 X-band guide (22.86 x 10.16 mm): dominant TE10 at c/(2a) = 6.557 GHz, next mode TE20 at
13.11 GHz -> the single-mode band. Self-contained, headless; SI units.
"""
import os
import sys
import math

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src"))
from ngsolve import Mesh, TaskManager
from netgen.occ import OCCGeometry, MoveTo
from radia_mcp.radia_ngsolve.waveguide import (rectangular_waveguide_cutoff, cutoff_frequency,
                                               helmholtz_cutoff_wavenumbers_2d, C0)

A, B = 0.02286, 0.01016          # WR-90 inner dimensions [m]


def build(maxh_frac=18):
    face = MoveTo(0, 0).Rectangle(A, B).Face()
    face.edges.name = "wall"
    return Mesh(OCCGeometry(face, dim=2).GenerateMesh(maxh=min(A, B) / maxh_frac))


with TaskManager():
    mesh = build()
    te_kc = helmholtz_cutoff_wavenumbers_2d(mesh, 3, bc="neumann")     # TE10, TE20, TE01
    tm_kc = helmholtz_cutoff_wavenumbers_2d(mesh, 2, bc="dirichlet")   # TM11, TM21

print(f"WR-90 guide  a={A*1e3:.2f} mm  b={B*1e3:.2f} mm   single-mode band = TE10..TE20")
print("mode   FE f_c [GHz]   analytic [GHz]    err")
for kc, (m, n, nm) in zip(te_kc, [(1, 0, "TE10"), (2, 0, "TE20"), (0, 1, "TE01")]):
    fc, fa = cutoff_frequency(kc), rectangular_waveguide_cutoff(A, B, m, n)
    print(f"{nm}   {fc/1e9:9.4f}     {fa/1e9:9.4f}     {(fc/fa-1)*100:+.3f}%")
for kc, (m, n, nm) in zip(tm_kc, [(1, 1, "TM11"), (2, 1, "TM21")]):
    fc, fa = cutoff_frequency(kc), rectangular_waveguide_cutoff(A, B, m, n)
    print(f"{nm}   {fc/1e9:9.4f}     {fa/1e9:9.4f}     {(fc/fa-1)*100:+.3f}%")

f_te10 = cutoff_frequency(te_kc[0])
print(f"\ndominant TE10 cutoff = {f_te10/1e9:.4f} GHz  (closed form c/2a = {C0/(2*A)/1e9:.4f} GHz)")
print("OK -- TE=Neumann / TM=Dirichlet Laplacian eigenvalues == exact rectangular-guide cutoffs.")
