"""COMSOL-class #9 -- CAPACITANCE MATRIX (AC/DC "Computing Capacitance" matrix).

The Maxwell capacitance matrix C (Q_i = sum_j C_ij V_j) of a multi-conductor system,
computed the COMSOL way: energise each conductor to 1 V (others grounded) and read every
conductor's charge -- here via the FEM reaction Q_i = <K V, chi_i>.

Validation: a spherical capacitor as a CLOSED 2-conductor system (inner sphere, outer
shell, vacuum gap, nothing outside) has the exact matrix

    C = [[ C0, -C0 ],
         [-C0,  C0 ]] ,   C0 = 4 pi eps0 a b / (b - a)

which exercises the diagonal (energy), the off-diagonal (induced charge), symmetry, and
the zero row-sum (charge conservation). Generalises to N conductors.
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src"))
import numpy as np
from ngsolve import Mesh
from netgen.csg import CSGeometry, Sphere, Pnt
from radia_mcp.radia_ngsolve.electrostatic3d import (
    capacitance_matrix, spherical_capacitor_C, EPS0)

a, b = 0.05, 0.10
sph_a = Sphere(Pnt(0, 0, 0), a).bc("inner")
sph_b = Sphere(Pnt(0, 0, 0), b).bc("outer")
geo = CSGeometry(); geo.Add((sph_b - sph_a).mat("gap").maxh(0.01))
mesh = Mesh(geo.GenerateMesh(maxh=0.012, grading=0.3)); mesh.Curve(3)
print(f"mesh: {mesh.ne} elems")

C = capacitance_matrix(mesh, EPS0, ["inner", "outer"], order=2)
C0 = spherical_capacitor_C(a, b)
C_ref = np.array([[C0, -C0], [-C0, C0]])

np.set_printoptions(precision=4, suppress=True)
print("C_fem  [pF]:\n", C * 1e12)
print("C_exact[pF]:\n", C_ref * 1e12)
err = np.max(np.abs(C - C_ref)) / C0 * 100
sym = np.max(np.abs(C - C.T)) / C0 * 100
rowsum = np.max(np.abs(C.sum(axis=1))) / C0 * 100
print(f"max entry err {err:.2f}%   asymmetry {sym:.3f}%   |row-sum|/C0 {rowsum:.3f}%")
print("OK" if err < 3 and sym < 1 and rowsum < 1 else "OFF")
