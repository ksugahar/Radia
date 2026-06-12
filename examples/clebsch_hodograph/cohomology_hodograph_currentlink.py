"""Does the hodograph need cohomology?  -- a current-linking demonstration.

RESEARCH example.  Answers the question "is cohomology needed for the
hodograph?".

The hodograph uses the magnetic SCALAR potential V (current-free region,
H = grad V) as one of its two conjugate coordinates (the other is the flux
psi / A_z).  V is single-valued ONLY IF the field has zero PERIOD -- zero
linked current -- around every hole:

    oint_C H . dl = (current linked by C).

So the rule is NOT "multiply connected => cohomology".  It is:

    simply connected (b1 = 0)                 -> V single-valued  (no cohomology)
    multiply connected, ZERO period           -> V single-valued  (no cut needed)
    multiply connected, NON-zero period       -> V MULTI-valued -> COHOMOLOGY
                                                 (the scalar coordinate IS a
                                                  1st-cohomology class)

This is why the simply-/zero-period demos (cylinder & sphere in a uniform
field, no current) need no cohomology, while a real magnet with a coil
window (a current threads the air) does.

We show it on a washer (b1 = 1) using the lab's gmsh-free engine
``radia.cohomology``: the unit-circulation harmonic generator ``h`` is
curl-free but NOT a gradient (oint h.dl = 1), so the current-linking field
``H = I h`` has no single-valued scalar coordinate.  Quantitatively, the
relative residual of the best single-valued grad-fit ``||F - grad V|| /
||F||`` is ~1 for the cohomology generator (NOT representable as a
gradient) and ~0 for a zero-period field (representable).  The solid
cylinder (b1 = 0) has no generator at all.

``radia.cohomology`` is the interim home for this (NGSolve/Netgen have the
de Rham complex but no turnkey homology generator; gmsh has Pellikka 2013).
"""
import math

from netgen.occ import Cylinder, Pnt, Dir, OCCGeometry
from ngsolve import (Mesh, H1, GridFunction, grad, InnerProduct, dx, CF,
                     BilinearForm, LinearForm, TaskManager, Integrate, curl)
from radia.cohomology import cohomology_basis, betti_numbers, circulation

_AX = Dir(0, 0, 1)


def _washer(maxh=0.02):
    shape = (Cylinder(Pnt(0, 0, -0.015), _AX, r=0.05, h=0.03)
             - Cylinder(Pnt(0, 0, -0.025), _AX, r=0.02, h=0.05))
    return Mesh(OCCGeometry(shape).GenerateMesh(maxh=maxh))


def _solid(maxh=0.02):
    return Mesh(OCCGeometry(Cylinder(Pnt(0, 0, -0.015), _AX, r=0.05, h=0.03))
                .GenerateMesh(maxh=maxh))


def _grad_fit_residual(F, mesh, order=2):
    """Best single-valued V minimising ||grad V - F|| (Galerkin projection of
    F onto gradients), returns ||F - grad V|| / ||F||.  ~0 if F IS a
    gradient (exact), ~1 if F is a pure cohomology generator (harmonic,
    L2-orthogonal to gradients)."""
    fes = H1(mesh, order=order)
    u, w = fes.TnT()
    # +eps mass removes the constant nullspace (RHS is already orthogonal to it).
    a = BilinearForm(InnerProduct(grad(u), grad(w)) * dx + 1e-8 * u * w * dx,
                     symmetric=True)
    f = LinearForm(InnerProduct(F, grad(w)) * dx)
    a.Assemble()
    f.Assemble()
    gV = GridFunction(fes)
    gV.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * f.vec
    num = Integrate(InnerProduct(F - grad(gV), F - grad(gV)), mesh)
    den = Integrate(InnerProduct(F, F), mesh)
    return math.sqrt(num / den)


def solve(maxh=0.02):
    """Return the verification dict (b1, generator checks, grad-fit residuals)."""
    with TaskManager():
        b0s, b1s = betti_numbers(_solid(maxh))

        mw = _washer(maxh)
        basis, b1, fes, ctx, loops = cohomology_basis(mw)
        h = basis[0]
        curl_rel = math.sqrt(Integrate(InnerProduct(curl(h), curl(h)), mw)
                             / Integrate(InnerProduct(h, h), mw))
        oint_hole = circulation(h, mw, 0, 0, 0.035)
        oint_contr = circulation(h, mw, 0.035, 0, 0.006)

        # current-linking field (period 1) -- NOT a gradient
        rho_h = _grad_fit_residual(h, mw)
        # zero-period curl-free field (uniform x_hat = grad x) -- IS a gradient
        F0 = CF((1.0, 0.0, 0.0))
        oint_F0 = circulation(F0, mw, 0, 0, 0.035)
        rho_g = _grad_fit_residual(F0, mw)

    return {
        "b1_solid": int(b1s), "b1_washer": int(b1),
        "curl_rel": float(curl_rel),
        "oint_hole": float(oint_hole), "oint_contractible": float(oint_contr),
        "residual_cohomology_field": float(rho_h),    # ~1: NOT a gradient
        "oint_zeroperiod": float(oint_F0),
        "residual_gradient_field": float(rho_g),       # ~0: IS a gradient
    }


def main():
    r = solve()
    print("Does the hodograph need cohomology?  -- current-linking demo\n")
    print(f"  solid cylinder b1 = {r['b1_solid']}  "
          f"(simply connected: scalar V always single-valued)")
    print(f"  washer         b1 = {r['b1_washer']}  (one hole)\n")
    print(f"  cohomology generator h (radia.cohomology, gmsh-free):")
    print(f"    ||curl h|| / ||h|| = {r['curl_rel']:.1e}    (curl-free)")
    print(f"    oint h.dl  hole = {r['oint_hole']:+.4f}   "
          f"contractible = {r['oint_contractible']:+.1e}\n")
    print(f"  Single-valued scalar coordinate V exists?  "
          f"(grad-fit residual ||F - grad V|| / ||F||)")
    print(f"    current-linking field (period {r['oint_hole']:+.2f}): "
          f"residual = {r['residual_cohomology_field']:.3f}  -> NO  "
          f"(cohomology REQUIRED)")
    print(f"    zero-period field     (period {r['oint_zeroperiod']:+.2f}): "
          f"residual = {r['residual_gradient_field']:.1e}  -> yes "
          f"(single-valued V)")
    print("\n  => the hodograph's SCALAR coordinate is a 1st-cohomology class")
    print("     exactly when a current threads a hole (period != 0); else it")
    print("     stays single-valued (as in the cylinder / sphere demos).")


if __name__ == "__main__":
    main()
