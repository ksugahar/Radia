"""2D bidirectional potential-coordinate map (Tampere -> Radia, Stage A / Phase 1).

Reproduces the exterior-calculus "potential coordinates" idea of Dervisha,
Marjamaki, Rasilo, Tarhasaari (CEFC 2026 WA-O1 #2) in NGSolve:

  FORWARD  (geometry -> potentials):  given (x, y), solve the flux function
           A(x, y) and the scalar potential phi(x, y).  A-isolines = field
           lines (B = grad(A) x z_hat), phi-isolines = equipotentials.
  INVERSE  (potentials -> geometry):  treat (A, phi) as the independent
           coordinates and solve x(A, phi), y(A, phi) -- "redraw the geometry
           in potential coordinates".

KEY (verified sympy 2026-06-11):
  * With mu = 1 the Tampere inverse-map PDE  d/dA(mu dx/dA) + d/dphi((1/mu)
    dx/dphi) = 0  reduces to  x, y HARMONIC in (A, phi)  -- a plain Laplace
    solve on the (A, phi) image domain.
  * Jacobian  d(A,phi)/d(x,y) = +-|B||H|  (= |grad A| |grad phi| under the
    conjugate-harmonic orthogonality grad(A) . grad(phi) = 0); the map is a
    valid chart exactly where |B||H| != 0.

Validation uses the analytic conjugate-harmonic pair from the complex potential
w(zeta) = zeta^2,  zeta = x + i y:
      A   = Re(w) = x^2 - y^2          (flux function)
      phi = Im(w) = 2 x y              (conjugate scalar potential, mu = 1)
on the square [1,2] x [1,2] (away from the |B|=0 degeneracy at the origin).
The exact inverse is  x + i y = sqrt(A + i phi).

Run:  python validation_test/feec/vim_legacy/bidirectional_map_2d.py
"""
import os
import sys
import math

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../../src/radia"))

from ngsolve import (
    Mesh, H1, GridFunction, CF, grad, InnerProduct, dx, sqrt, atan2,
    x, y, BilinearForm, LinearForm, TaskManager, Integrate, cos, sin,
)
from netgen.geom2d import SplineGeometry


# ---------------------------------------------------------------------------
# Analytic conjugate-harmonic pair  w = zeta^2
# ---------------------------------------------------------------------------
X0, X1, Y0, Y1 = 1.0, 2.0, 1.0, 2.0   # physical square, away from origin

A_exact = x**2 - y**2          # flux function (NGSolve coords x,y = physical)
PHI_exact = 2 * x * y          # conjugate scalar potential


def _inverse_cf():
    """Analytic inverse x(A,phi), y(A,phi) as CFs on the (A,phi) mesh.

    On the image mesh NGSolve's coordinates (x, y) ARE (A, phi).  The principal
    branch of sqrt(A + i phi):  r = |w|, x = sqrt(r) cos(theta/2), etc.
    """
    A, phi = x, y                       # image-mesh coords are (A, phi)
    r = sqrt(A * A + phi * phi)
    rt = sqrt(r)
    th = atan2(phi, A)
    return rt * cos(th / 2), rt * sin(th / 2)


# ---------------------------------------------------------------------------
# FORWARD: geometry -> potentials
# ---------------------------------------------------------------------------
def solve_forward(maxh=0.04, order=3):
    geo = SplineGeometry()
    geo.AddRectangle((X0, Y0), (X1, Y1), bc="bnd")
    mesh = Mesh(geo.GenerateMesh(maxh=maxh))

    fes = H1(mesh, order=order, dirichlet="bnd")
    u, v = fes.TnT()
    a = BilinearForm(fes, symmetric=True)
    a += InnerProduct(grad(u), grad(v)) * dx
    f = LinearForm(fes)
    a.Assemble()
    f.Assemble()
    inv = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky")

    gfA = GridFunction(fes, name="A")
    gfP = GridFunction(fes, name="phi")
    gfA.Set(A_exact, definedon=mesh.Boundaries("bnd"))
    gfP.Set(PHI_exact, definedon=mesh.Boundaries("bnd"))
    # harmonic extension: u = -K^-1 (a * u_bnd) on free dofs
    for gf, ex in ((gfA, A_exact), (gfP, PHI_exact)):
        r = gf.vec.CreateVector()
        r.data = f.vec - a.mat * gf.vec
        gf.vec.data += inv * r
    return mesh, gfA, gfP


def forward_checks(mesh, gfA, gfP):
    """Return dict of forward-map diagnostics."""
    # L2 accuracy vs analytic
    errA = sqrt(Integrate((gfA - A_exact) ** 2, mesh) /
                Integrate(A_exact ** 2, mesh))
    errP = sqrt(Integrate((gfP - PHI_exact) ** 2, mesh) /
                Integrate(PHI_exact ** 2, mesh))

    gA, gP = grad(gfA), grad(gfP)
    gA2 = InnerProduct(gA, gA)
    gP2 = InnerProduct(gP, gP)
    # orthogonality  grad(A) . grad(phi) ~ 0  (relative)
    ortho = sqrt(Integrate(InnerProduct(gA, gP) ** 2, mesh) /
                 Integrate(gA2 * gP2, mesh))
    # Jacobian det[grad A; grad phi] vs |grad A||grad phi| (= |B||H|, mu=1)
    detJ = gA[0] * gP[1] - gA[1] * gP[0]
    jac_rel = sqrt(Integrate((detJ ** 2 - gA2 * gP2) ** 2, mesh) /
                   Integrate((gA2 * gP2) ** 2, mesh))
    return {
        "errA": float(errA), "errP": float(errP),
        "ortho": float(ortho), "jac_rel": float(jac_rel),
    }


# ---------------------------------------------------------------------------
# INVERSE: potentials -> geometry  (x, y harmonic in (A, phi), mu = 1)
# ---------------------------------------------------------------------------
def _image_mesh(n_per_edge=60, maxh=0.4):
    """Mesh the (A, phi) image of the physical square (dense polygon boundary)."""
    def fwd(px, py):
        return (px * px - py * py, 2 * px * py)

    pts = []
    # walk the square boundary CCW: bottom, right, top, left
    edges = [
        [(X0 + t * (X1 - X0), Y0) for t in np.linspace(0, 1, n_per_edge, endpoint=False)],
        [(X1, Y0 + t * (Y1 - Y0)) for t in np.linspace(0, 1, n_per_edge, endpoint=False)],
        [(X1 - t * (X1 - X0), Y1) for t in np.linspace(0, 1, n_per_edge, endpoint=False)],
        [(X0, Y1 - t * (Y1 - Y0)) for t in np.linspace(0, 1, n_per_edge, endpoint=False)],
    ]
    for e in edges:
        for (px, py) in e:
            pts.append(fwd(px, py))

    geo = SplineGeometry()
    ids = [geo.AppendPoint(ax, ay) for (ax, ay) in pts]
    n = len(ids)
    for i in range(n):
        geo.Append(["line", ids[i], ids[(i + 1) % n]], bc="bnd")
    return Mesh(geo.GenerateMesh(maxh=maxh))


def solve_inverse(order=3, n_per_edge=60, maxh=0.4):
    mesh = _image_mesh(n_per_edge=n_per_edge, maxh=maxh)
    x_ex, y_ex = _inverse_cf()

    fes = H1(mesh, order=order, dirichlet="bnd")
    u, v = fes.TnT()
    a = BilinearForm(fes, symmetric=True)
    a += InnerProduct(grad(u), grad(v)) * dx
    f = LinearForm(fes)
    a.Assemble()
    f.Assemble()
    inv = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky")

    gfx = GridFunction(fes, name="x")
    gfy = GridFunction(fes, name="y")
    for gf, ex in ((gfx, x_ex), (gfy, y_ex)):
        gf.Set(ex, definedon=mesh.Boundaries("bnd"))
        r = gf.vec.CreateVector()
        r.data = f.vec - a.mat * gf.vec
        gf.vec.data += inv * r

    errx = sqrt(Integrate((gfx - x_ex) ** 2, mesh) / Integrate(x_ex ** 2, mesh))
    erry = sqrt(Integrate((gfy - y_ex) ** 2, mesh) / Integrate(y_ex ** 2, mesh))
    return mesh, gfx, gfy, {"errx": float(errx), "erry": float(erry)}


def roundtrip_error(mesh_inv, gfx, gfy, samples=None):
    """Physical point -> (A,phi) -> inverse map -> physical point; max error."""
    if samples is None:
        samples = [(1.3, 1.2), (1.5, 1.5), (1.7, 1.4), (1.2, 1.8), (1.9, 1.1)]
    worst = 0.0
    for (px, py) in samples:
        A = px * px - py * py
        phi = 2 * px * py
        mp = mesh_inv(A, phi)
        xr = float(gfx(mp))
        yr = float(gfy(mp))
        worst = max(worst, math.hypot(xr - px, yr - py))
    return worst


# ---------------------------------------------------------------------------
def main():
    print("=== 2D bidirectional potential-coordinate map (Tampere) ===")
    with TaskManager():
        mesh_f, gfA, gfP = solve_forward()
        fchk = forward_checks(mesh_f, gfA, gfP)
        mesh_i, gfx, gfy, ichk = solve_inverse()
        rt = roundtrip_error(mesh_i, gfx, gfy)

    print("\n[forward]  geometry -> (A, phi)")
    print(f"   L2 err A           = {fchk['errA']:.2e}")
    print(f"   L2 err phi         = {fchk['errP']:.2e}")
    print(f"   orthogonality      = {fchk['ortho']:.2e}   (grad A . grad phi ~ 0)")
    print(f"   Jacobian rel-err   = {fchk['jac_rel']:.2e}   (det = |B||H|)")
    print("\n[inverse]  (A, phi) -> geometry  (x,y harmonic in (A,phi))")
    print(f"   L2 err x(A,phi)    = {ichk['errx']:.2e}")
    print(f"   L2 err y(A,phi)    = {ichk['erry']:.2e}")
    print(f"\n[round-trip]  p -> (A,phi) -> p   max error = {rt:.2e}")

    ok = (fchk["ortho"] < 1e-3 and fchk["jac_rel"] < 1e-2 and
          fchk["errA"] < 1e-3 and ichk["errx"] < 2e-2 and rt < 2e-2)
    print(f"\n=== PASS: {ok} ===")
    return ok


if __name__ == "__main__":
    raise SystemExit(0 if main() else 1)
