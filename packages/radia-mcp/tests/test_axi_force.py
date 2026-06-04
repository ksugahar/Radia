"""Axisymmetric eggshell force (force.eggshell_force_axi) validated on two
coaxial current loops against the EXACT mutual-inductance force

    F_z = I1 I2 dM/dz,

with M(z) the Maxwell mutual inductance of two coaxial loops (radii a0,
separation z) via complete elliptic integrals K, E:

    k^2 = 4 a0^2 / (4 a0^2 + z^2),
    M   = mu0 a0 [ (2/k - k) K(k) - (2/k) E(k) ].

dM/dz is taken by high-precision finite difference of the analytic M, so the
reference is exact (filamentary). The FE coils have a small cross-section to
approximate filaments. Validates the 2*pi*r-weighted meridional Maxwell stress.
"""
import math
import os
import sys

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

import scipy.special as sp
from ngsolve import Mesh, CoefficientFunction, grad, TaskManager
from ngsolve import x as r_cf
from netgen.occ import OCCGeometry, MoveTo, WorkPlane, Glue, X, Y
from radia_mcp.radia_ngsolve.solve import solve_axi_magnetostatic
from radia_mcp.radia_ngsolve.force import eggshell_force_axi

MU0 = 4e-7 * math.pi
A0 = 1.0        # loop radius
D = 1.5         # axial separation
W = 0.08        # coil cross-section side
I1 = I2 = 1.0   # loop currents
R_OUT = 8.0


def M_coax(a0, z):
    m = 4.0 * a0 * a0 / (4.0 * a0 * a0 + z * z)   # = k^2 (scipy uses parameter m)
    k = math.sqrt(m)
    return MU0 * a0 * ((2.0 / k - k) * sp.ellipk(m) - (2.0 / k) * sp.ellipe(m))


def Fz_exact(a0, z, i1, i2):
    de = 1e-6 * z
    return i1 * i2 * (M_coax(a0, z + de) - M_coax(a0, z - de)) / (2.0 * de)


def build_mesh():
    half = MoveTo(0, -R_OUT).Rectangle(R_OUT, 2 * R_OUT).Face()
    c1 = MoveTo(A0 - W / 2, -W / 2).Rectangle(W, W).Face()
    c1.faces.name = "coil1"; c1.maxh = 0.02
    c2 = MoveTo(A0 - W / 2, D - W / 2).Rectangle(W, W).Face()
    c2.faces.name = "coil2"; c2.maxh = 0.02
    near2 = (WorkPlane().Circle(A0, D, 0.5).Face() * half) - c2
    near2.faces.name = "air"; near2.maxh = 0.04
    outer = half - c1 - c2 - near2
    outer.faces.name = "air"; outer.maxh = 0.5
    outer.edges.Min(X).name = "axis"
    outer.edges.Max(X).name = "outer"
    outer.edges.Max(Y).name = "outer"
    outer.edges.Min(Y).name = "outer"
    return Mesh(OCCGeometry(Glue([outer, near2, c1, c2]), dim=2).GenerateMesh(maxh=0.5))


def main():
    mesh = build_mesh()
    area = W * W
    Jr = mesh.MaterialCF({"coil1": I1 / area, "coil2": I2 / area}, default=0.0)
    nu = CoefficientFunction(1.0 / MU0)

    with TaskManager():
        gfu = solve_axi_magnetostatic(mesh, nu, Jr=Jr, order=2)
    B = CoefficientFunction((-grad(gfu)[1], grad(gfu)[0] + gfu / r_cf))  # (B_r, B_z)

    Fz = eggshell_force_axi(B, mesh, center=(A0, D), r_inner=0.15, r_outer=0.35)
    Fz_ref = Fz_exact(A0, D, I1, I2)
    err = (Fz - Fz_ref) / abs(Fz_ref)
    print(f"  coaxial loops a0={A0} d={D}: F_z FEM={Fz:.6e} N  "
          f"exact(I1 I2 dM/dz)={Fz_ref:.6e} N  err={100*err:+.2f}%")
    assert abs(err) < 0.05, f"axi eggshell force off {100*err:.2f}% (>5%)"
    print(f"\n[OK] eggshell_force_axi validated vs coaxial-loop dM/dz "
          f"({100*abs(err):.2f}%).")


if __name__ == "__main__":
    main()
