"""2D PLANAR A_z magnetostatics with PERIODIC / ANTI-PERIODIC boundaries -- the
FEMM "periodic / anti-periodic" BC used to solve ONE pole/slot pitch of a motor
or machine instead of the full cross-section. Exercises
radia_mcp.radia_ngsolve.solve.periodic_h1 and validates against a closed form.

Strip [0,W]x[0,H]: left (x=0) identified onto right (x=W), top/bottom A_z=0, nu=1:
    -div(grad A) = Jz
  PERIODIC       Jz = sin(k x), k = 2pi/W  (W-periodic):
        A = sin(kx)/k^2 * [1 - cosh(k(y-H/2))/cosh(kH/2)]
  ANTI-PERIODIC  Jz = cos(q x), q =  pi/W  (A(x+W) = -A(x)):
        A = cos(qx)/q^2 * [1 - cosh(q(y-H/2))/cosh(qH/2)]
Both recovered to ~1e-7 relative L2 -- machine precision for the discretization.
"""
import math
import os
import sys

# Test the DEV-repo source (radia_mcp is installed non-editable in site-packages).
_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from ngsolve import (Mesh, GridFunction, BilinearForm, LinearForm, grad, dx,
                     sin, cos, exp, x, y, H1, Integrate, TaskManager)
from netgen.geom2d import SplineGeometry
from radia_mcp.radia_ngsolve.solve import periodic_h1

W, H = 1.0, 1.0


def _cosh(z):
    return (exp(z) + exp(-z)) / 2


def periodic_strip_mesh(maxh=0.025):
    """Rectangle [0,W]x[0,H] with the left edge periodically identified onto the
    right edge. The boundary loop is built CCW (p0->p1->p2->p3) so the periodic
    ``copy=`` identification meshes cleanly (a reversed loop stalls the mesher)."""
    geo = SplineGeometry()
    p = [geo.AppendPoint(*q) for q in [(0, 0), (W, 0), (W, H), (0, H)]]
    geo.Append(["line", p[0], p[1]], bc="bottom")
    rt = geo.Append(["line", p[1], p[2]], bc="right")        # master x=W
    geo.Append(["line", p[2], p[3]], bc="top")
    geo.Append(["line", p[3], p[0]], bc="left", copy=rt)     # slave x=0 -> master
    return Mesh(geo.GenerateMesh(maxh=maxh))


def _solve(fes, source):
    u, v = fes.TnT()
    a = BilinearForm(fes, symmetric=True)
    a += grad(u) * grad(v) * dx
    f = LinearForm(fes)
    f += source * v * dx
    a.Assemble()
    f.Assemble()
    gfu = GridFunction(fes)
    gfu.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * f.vec
    return gfu


def _relerr(gf, exact, mesh):
    return math.sqrt(Integrate((gf - exact) ** 2, mesh) / Integrate(exact ** 2, mesh))


def main():
    mesh = periodic_strip_mesh()
    with TaskManager():
        # --- PERIODIC ---
        k = 2 * math.pi / W
        A_per = sin(k * x) / k**2 * (1 - _cosh(k * (y - H / 2)) / math.cosh(k * H / 2))
        fes = periodic_h1(mesh, order=3, dirichlet="bottom|top")
        gp = _solve(fes, sin(k * x))
        e_per = _relerr(gp, A_per, mesh)
        seam = abs(gp(mesh(1e-9, 0.5)) - gp(mesh(W - 1e-9, 0.5)))
        print(f"  PERIODIC      rel L2 = {e_per:.2e}   |A(0)-A(W)| = {seam:.1e}")

        # --- ANTI-PERIODIC (phase=-1, complex space) ---
        q = math.pi / W
        A_anti = cos(q * x) / q**2 * (1 - _cosh(q * (y - H / 2)) / math.cosh(q * H / 2))
        fesc = periodic_h1(mesh, order=3, dirichlet="bottom|top", antiperiodic=True)
        ga = _solve(fesc, cos(q * x))
        gre = GridFunction(H1(mesh, order=3))
        gre.Set(ga.real)
        e_anti = _relerr(gre, A_anti, mesh)
        aseam = abs(ga(mesh(1e-9, 0.5)).real + ga(mesh(W - 1e-9, 0.5)).real)
        print(f"  ANTI-PERIODIC rel L2 = {e_anti:.2e}   |A(0)+A(W)| = {aseam:.1e}")

    assert e_per < 1e-3, f"periodic off by rel L2 {e_per:.2e}"
    assert e_anti < 1e-3, f"anti-periodic off by rel L2 {e_anti:.2e}"
    print("\n[OK] planar periodic + anti-periodic BC (FEMM motor-sector analog) "
          f"validated on NGSolve Periodic (periodic {e_per:.1e}, anti {e_anti:.1e}).")


if __name__ == "__main__":
    main()
