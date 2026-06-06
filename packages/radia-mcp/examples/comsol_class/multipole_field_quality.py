"""COMSOL-class #13 -- MULTIPOLE field-quality analysis of a 2D magnet.

The post-processing an accelerator-magnet designer runs after every solve: take
the FEM flux density on a reference circle and decompose it into normal/skew
multipoles (b_n, a_n), then read the main harmonic and the field errors. Here we
build an air-cored NORMAL QUADRUPOLE from four line currents (+I,-I,+I,-I at
45/135/225/315 deg, radius r0), solve the planar A_z field, and extract the
harmonics about the origin.

Reference (exact, free space): outside a cylindrically-symmetric wire the field
equals a line current at its centre, so for R_ref < r0 the bore harmonics are
EXACTLY the superposition of line_current_multipoles -- no FEM, no convention
ambiguity. The 4-fold symmetry leaves only the ALLOWED quadrupole harmonics
n = 2, 6, 10, ...  (the 12-pole n=6 sits at (R_ref/r0)^4 of the main); the
"forbidden" n = 1,3,4,5,7,8,9 must vanish. Convention (CERN/European):

    B_y + i B_x = sum_{n>=1} (b_n + i a_n) (z/R_ref)^{n-1},  n=1 dipole, n=2 quad.
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src"))
import math
from ngsolve import Mesh, CoefficientFunction, grad, Integrate, dx
from netgen.geom2d import SplineGeometry
from radia_mcp.radia_ngsolve.solve import solve_planar_magnetostatic, NU0
from radia_mcp.radia_ngsolve.fieldquality import (multipole_coefficients,
    line_current_multipoles, superpose_multipoles, field_errors)

I, r0, rw, R_ref, Rout = 1000.0, 0.05, 0.004, 0.02, 0.40   # A, m
signs = [+1, -1, +1, -1]
angles = [math.radians(90 * k) for k in range(4)]           # 0,90,180,270 deg -> NORMAL quad
n_max = 10

# --- geometry: 4 wire disks in air, far Dirichlet boundary -------------------
geo = SplineGeometry()
geo.AddCircle((0, 0), Rout, leftdomain=1, rightdomain=0, bc="outer")
for k, ang in enumerate(angles):
    c = (r0 * math.cos(ang), r0 * math.sin(ang))
    geo.AddCircle(c, rw, leftdomain=2 + k, rightdomain=1)
    geo.SetMaterial(2 + k, f"w{k}")
geo.SetMaterial(1, "air")
mesh = Mesh(geo.GenerateMesh(maxh=0.015))
print(f"mesh: {mesh.ne} elems, materials={mesh.GetMaterials()}")

# --- source: normalize J on the ACTUAL meshed wire area so int Jz = I exactly
# (the field at R_ref depends only on each wire's TOTAL current, not its shape;
#  using pi*rw^2 instead of the discrete material area is the dominant FEM error).
Jz_terms = {}
for k in range(4):
    area_k = Integrate(mesh.MaterialCF({f"w{k}": 1.0}, default=0.0) * dx, mesh)
    Jz_terms[f"w{k}"] = signs[k] * I / area_k
Jz = mesh.MaterialCF(Jz_terms, default=0.0)
gfu = solve_planar_magnetostatic(mesh, CoefficientFunction(NU0), Jz=Jz, order=4)
B = CoefficientFunction((grad(gfu)[1], -grad(gfu)[0]))      # (Bx, By)

# --- extract multipoles from the FEM field -----------------------------------
fem = multipole_coefficients(B, R_ref, n_max=n_max, mesh=mesh)

# --- exact free-space reference (superpose the four line currents) -----------
refs = [line_current_multipoles(signs[k] * I, r0, angles[k], R_ref, n_max)
        for k in range(4)]
ana = superpose_multipoles(*refs)

fe = field_errors(fem)
print(f"\nmain harmonic: n={fe['main_n']} (expect 2 = quadrupole), "
      f"|C2|={fe['main_C']:.4e} T at R_ref={R_ref*1e3:.0f} mm")
print(" n  type      C_fem[T]    C_ana[T]   rel_fem   rel_ana   (units 1e-4 of main)")
ALLOWED = {2, 6, 10}
worst_allowed = 0.0
worst_forbidden = 0.0
for n in range(1, n_max + 1):
    rel_fem = fem["C"][n] / fe["main_C"] * 1e4
    rel_ana = ana["C"][n] / abs(complex(ana["b"][2], ana["a"][2])) * 1e4
    tag = "ALLOWED" if n in ALLOWED else "forbid "
    print(f"{n:2d}  {tag}  {fem['C'][n]:.4e}  {ana['C'][n]:.4e}  "
          f"{rel_fem:8.2f}  {rel_ana:8.2f}")
    if n in ALLOWED:
        if ana["C"][n] > 0:
            worst_allowed = max(worst_allowed,
                                abs(fem["C"][n] - ana["C"][n]) / ana["C"][n])
    else:
        worst_forbidden = max(worst_forbidden, abs(rel_fem))

b2_err = abs(fem["b"][2] - ana["b"][2]) / abs(ana["b"][2]) * 100
print(f"\nb2 (gradient*R_ref): FEM {fem['b'][2]:+.4e}  analytic {ana['b'][2]:+.4e}  "
      f"err {b2_err:.2f}%")
print(f"allowed-harmonic worst err: {worst_allowed*100:.2f}%   "
      f"forbidden leakage worst: {worst_forbidden:.2f} units (1e-4)")
ok = b2_err < 1.0 and worst_allowed < 0.03 and worst_forbidden < 5.0
print("OK" if ok else "OFF")
