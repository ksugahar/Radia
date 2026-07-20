r"""AC skin effect in a conducting slab -- the complex NGSolve eddy solve vs the closed form.

A slab x in [-a, a], conductivity sigma, mu_r=1, driven by the surface phasor A_z=A0 on
x=+/-a (natural BC on the transverse faces) obeys the frequency-domain eddy diffusion
    lap(A_z) = j*omega*mu0*sigma*A_z .
The 1-D mid-plane attenuation is closed-form:
    A_z(0)/A0 = 1/cosh(gamma a),  gamma = (1+j)/delta,  delta = sqrt(2/(omega mu0 sigma)).

Gate: the complex H1 (order 5) NGSolve solve reproduces skin_effect_slab_attenuation at the
mid-plane to ~1e-9 (and the surface phasor is held).  Pure analytic reference; this is the
NGSolve leg of the slab-skin-effect cross-check (the field-diffusion / lamination-loss case).
"""
import cmath
import math
import os
import sys

from ngsolve import (Mesh, H1, GridFunction, BilinearForm, LinearForm, grad, dx,
                     BND, TaskManager)
from netgen.geom2d import SplineGeometry

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)
from radia_mcp.radia_ngsolve.solve import skin_effect_slab_attenuation

MU0 = 4e-7 * math.pi
SIGMA, F, A_HALF, B_H, A0 = 1.0e7, 60.0, 0.03, 0.02, 1.0


def _solve_slab_ngsolve():
    omega = 2 * math.pi * F
    geo = SplineGeometry()
    p = [geo.AppendPoint(*q) for q in
         [(-A_HALF, 0), (A_HALF, 0), (A_HALF, B_H), (-A_HALF, B_H)]]
    geo.Append(["line", p[0], p[1]], bc="yface")   # natural BC -> 1-D in x
    geo.Append(["line", p[1], p[2]], bc="surf")
    geo.Append(["line", p[2], p[3]], bc="yface")
    geo.Append(["line", p[3], p[0]], bc="surf")
    mesh = Mesh(geo.GenerateMesh(maxh=0.003))
    fes = H1(mesh, order=5, complex=True, dirichlet="surf")
    u, v = fes.TnT()
    a = BilinearForm(fes)
    a += (1.0 / MU0) * grad(u) * grad(v) * dx + 1j * omega * SIGMA * u * v * dx
    f = LinearForm(fes); f.Assemble()
    gfu = GridFunction(fes)
    gfu.Set(A0, BND, definedon=mesh.Boundaries("surf"))
    with TaskManager():
        a.Assemble()
        r = f.vec.CreateVector(); r.data = f.vec - a.mat * gfu.vec
        gfu.vec.data += a.mat.Inverse(fes.FreeDofs(), inverse="umfpack") * r
    return mesh, gfu


def test_skin_effect_slab_ngsolve_matches_closed_form():
    mesh, gfu = _solve_slab_ngsolve()
    atten_closed = skin_effect_slab_attenuation(SIGMA, F, A_HALF)
    atten_ng = gfu(mesh(0.0, B_H / 2)) / A0
    rel = abs(atten_ng - atten_closed) / abs(atten_closed)
    assert abs(atten_closed) < 0.6, "expected genuine attenuation (a/delta ~ 1.46)"
    assert rel < 1e-7, f"NGSolve mid-plane off closed form by {rel:.2e}"
    surf = gfu(mesh(A_HALF - 1e-6, B_H / 2))            # surface phasor held ~ A0
    assert abs(surf - A0) / A0 < 5e-3


def main():
    mesh, gfu = _solve_slab_ngsolve()
    atten_closed = skin_effect_slab_attenuation(SIGMA, F, A_HALF)
    atten_ng = gfu(mesh(0.0, B_H / 2)) / A0
    delta = math.sqrt(2.0 / (2 * math.pi * F * MU0 * SIGMA))
    print(f"delta={delta*1e3:.3f} mm  a/delta={A_HALF/delta:.3f}")
    print(f"closed  A(0)/A0 = {atten_closed:.8f}")
    print(f"NGSolve A(0)/A0 = {atten_ng:.8f}")
    print(f"rel err = {abs(atten_ng-atten_closed)/abs(atten_closed):.2e}")
    print("[OK] complex NGSolve eddy slab == 1/cosh(gamma a) (skin-effect attenuation).")


if __name__ == "__main__":
    main()
