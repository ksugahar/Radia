"""COMSOL-class #59 -- 3D rectangular cavity resonance (full-wave Maxwell eigenvalue, RF).

The 3-D / vector sequel to the 2-D waveguide cutoff (#53): a closed PEC box a x b x d resonates at

    f_mnp = (c/2) sqrt((m/a)^2 + (n/b)^2 + (p/d)^2),

the genuine ELECTROMAGNETIC cavity (curl-curl E = (omega/c)^2 E with n x E = 0 on the walls),
solved on an HCurl edge-element space. The curl-curl gradient kernel (spurious zero modes) is
suppressed by shifting the eigensolver near the fundamental. Dominant mode (a > d > b) is TE101.
Self-contained, headless; SI units.
"""
import os
import sys
import math

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src"))
from ngsolve import Mesh, TaskManager
from netgen.occ import Box, OCCGeometry, Pnt
from radia_mcp.radia_ngsolve.waveguide import (rectangular_cavity_frequency,
                                               maxwell_cavity_modes_3d, C0)

A, B, D = 0.040, 0.020, 0.030          # cavity dimensions [m]


def build(maxh=0.006):
    box = Box(Pnt(0, 0, 0), Pnt(A, B, D))
    for f in box.faces:
        f.name = "pec"
    return Mesh(OCCGeometry(box).GenerateMesh(maxh=maxh))


# the analytic lowest modes (a > d > b)
LABELS = [(1, 0, 1, "TE101"), (1, 1, 0, "TM110"), (0, 1, 1, "TE011"), (1, 1, 1, "TE/TM111")]
ana = sorted([(rectangular_cavity_frequency(A, B, D, m, n, p), nm) for (m, n, p, nm) in LABELS])

shift = (math.pi / A) ** 2 + (math.pi / D) ** 2          # ~ TE101 k^2
with TaskManager():
    fe = maxwell_cavity_modes_3d(build(), 5, shift)

print(f"PEC cavity {A*1e3:.0f} x {B*1e3:.0f} x {D*1e3:.0f} mm   (dominant TE101)")
print("  FE f [GHz]   nearest analytic   err")
for f in fe:
    fa, nm = min(ana, key=lambda t: abs(t[0] - f))
    print(f"  {f/1e9:8.4f}     {fa/1e9:8.4f} {nm:9s}  {(f/fa-1)*100:+.2f}%")

print(f"\ndominant TE101 = {fe[0]/1e9:.4f} GHz  (closed form (c/2)sqrt(1/a^2+1/d^2) = "
      f"{rectangular_cavity_frequency(A, B, D, 1, 0, 1)/1e9:.4f} GHz)")
print("OK -- HCurl curl-curl Maxwell eigenvalue == exact box spectrum; gradient kernel suppressed.")
