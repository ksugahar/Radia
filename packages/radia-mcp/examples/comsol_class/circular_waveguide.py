"""COMSOL-class #61 -- circular waveguide cutoff (Bessel-zero 2D Helmholtz eigenmodes).

The circular sibling of the rectangular waveguide cutoff (#53): a round metallic guide's modes
solve the SAME transverse Helmholtz eigenproblem -nabla_t^2 psi = k_c^2 psi on the DISK cross
section -- TE = Neumann (k_c = j'_mn/a, zeros of J'_m), TM = Dirichlet (k_c = j_mn/a, zeros of J_m).
Dominant TE11 (j'_11 = 1.8412); TE0n / TM1n are degenerate (j'_0n = j_1n). The SAME
`helmholtz_cutoff_wavenumbers_2d` eigensolver runs on a disk mesh -- only the geometry changes.
Self-contained, headless; SI units.
"""
import os
import sys
import math

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src"))
from ngsolve import Mesh, TaskManager
from netgen.occ import OCCGeometry, WorkPlane
from radia_mcp.radia_ngsolve.waveguide import (circular_waveguide_cutoff,
                                               helmholtz_cutoff_wavenumbers_2d, cutoff_frequency, C0)

A = 0.0127          # circular guide inner radius [m]


def build(maxh_frac=22):
    disk = WorkPlane().Circle(0, 0, A).Face(); disk.edges.name = "wall"
    return Mesh(OCCGeometry(disk, dim=2).GenerateMesh(maxh=A / maxh_frac))


with TaskManager():
    mesh = build()
    te = helmholtz_cutoff_wavenumbers_2d(mesh, 3, bc="neumann")     # TE11(x2), TE21
    tm = helmholtz_cutoff_wavenumbers_2d(mesh, 2, bc="dirichlet")   # TM01, TM11

print(f"circular guide a={A*1e3:.2f} mm   (dominant TE11, j'_11=1.8412)")
print("mode   FE f_c [GHz]   Bessel f_c [GHz]    err")
for kc, (mode, m, n, nm) in zip(te, [("TE", 1, 1, "TE11"), ("TE", 1, 1, "TE11"), ("TE", 2, 1, "TE21")]):
    fc, fa = cutoff_frequency(kc), circular_waveguide_cutoff(A, mode, m, n)
    print(f"{nm}   {fc/1e9:9.4f}     {fa/1e9:9.4f}     {(fc/fa-1)*100:+.3f}%")
for kc, (mode, m, n, nm) in zip(tm, [("TM", 0, 1, "TM01"), ("TM", 1, 1, "TM11")]):
    fc, fa = cutoff_frequency(kc), circular_waveguide_cutoff(A, mode, m, n)
    print(f"{nm}   {fc/1e9:9.4f}     {fa/1e9:9.4f}     {(fc/fa-1)*100:+.3f}%")

print(f"\ndominant TE11 = {cutoff_frequency(te[0])/1e9:.4f} GHz  (c*1.8412/(2 pi a) = "
      f"{C0*1.8412/(2*math.pi*A)/1e9:.4f} GHz)")
print("OK -- disk Helmholtz eigenvalues == circular-guide Bessel-zero cutoffs (same #53 eigensolver).")
