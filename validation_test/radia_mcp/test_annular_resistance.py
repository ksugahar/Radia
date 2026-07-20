r"""Annular conductor radial DC resistance: NGSolve conduction solve vs the closed form.

An annulus a<r<b, conductivity sigma, depth L, V0 on the inner ring and ground on the outer
carries a radial current; -div(sigma grad V)=0 gives R = ln(b/a)/(2 pi sigma L).  Gate: the
order-4 H1 conduction solve (R = V0^2 / power) reproduces annular_radial_resistance to ~1e-8.
Pure analytic reference (no commercial tooling).
"""
import math
import os
import sys

from ngsolve import (Mesh, H1, GridFunction, BilinearForm, LinearForm, grad, dx,
                     BND, Integrate, TaskManager)
from netgen.geom2d import SplineGeometry

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)
from radia_mcp.radia_ngsolve.solve import annular_radial_resistance

SIGMA, A, B, L, V0 = 1.0, 0.01, 0.05, 1.0, 1.0


def _solve_R_ngsolve():
    geo = SplineGeometry()
    geo.AddCircle((0, 0), B, leftdomain=1, rightdomain=0, bc="outer")
    geo.AddCircle((0, 0), A, leftdomain=0, rightdomain=1, bc="inner")
    mesh = Mesh(geo.GenerateMesh(maxh=0.002)); mesh.Curve(4)
    fes = H1(mesh, order=4, dirichlet="inner|outer")
    u, v = fes.TnT()
    a = BilinearForm(fes); a += SIGMA * grad(u) * grad(v) * dx
    f = LinearForm(fes); f.Assemble()
    gfu = GridFunction(fes)
    gfu.Set(mesh.BoundaryCF({"inner": V0, "outer": 0.0}), BND)
    with TaskManager():
        a.Assemble()
        r = f.vec.CreateVector(); r.data = f.vec - a.mat * gfu.vec
        gfu.vec.data += a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * r
    power = Integrate(SIGMA * grad(gfu) * grad(gfu) * dx, mesh)   # per unit depth
    return V0 * V0 / (power * L)


def test_annular_resistance_ngsolve_matches_closed_form():
    R_closed = annular_radial_resistance(SIGMA, A, B, L)
    assert R_closed == math.log(B / A) / (2 * math.pi * SIGMA * L)
    R_ng = _solve_R_ngsolve()
    assert abs(R_ng - R_closed) / R_closed < 1e-6


def main():
    R_closed = annular_radial_resistance(SIGMA, A, B, L)
    R_ng = _solve_R_ngsolve()
    print(f"closed  R = ln(b/a)/(2 pi sigma L) = {R_closed:.8f} ohm")
    print(f"NGSolve R = V0^2/(L*power)         = {R_ng:.8f} ohm")
    print(f"rel err = {abs(R_ng-R_closed)/R_closed:.2e}")
    print("[OK] NGSolve radial conduction == ln(b/a)/(2 pi sigma L).")


if __name__ == "__main__":
    main()
