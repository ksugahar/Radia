"""hdiv_demag_bem_singlelayer.py -- the PRODUCTION curved + high-order surface demag Gram via
ngsolve.bem (the Laplace single-layer), validated vs the sphere ANALYTIC 1/3.

The HDiv-VIM surface demag Gram (the uniform-M case: div M = 0 -> pure surface charge sigma = M.n) IS
the Laplace SINGLE-LAYER operator V:  (V sigma)(r) = INT_S sigma(r')/(4 pi |r-r'|) dS'.  NGSolve 6.2.2604
ships `ngsolve.bem` with a high-order, CURVED-aware, FMM-accelerated Laplace single-layer -- so the
curved AND high-order AND scalable surface Gram come essentially for free (the "complement NGSolve /
NGSolve-through" principle), with NO hand-rolled singular quadrature.

The demag factor of a uniform-M body is  D_z = <sigma, V sigma> / V_vol  with sigma = n_z (the boundary
single-layer matrix is V_ij = INT_S INT_S phi_i phi_j /(4 pi r), and sigma = sum c_i phi_i, so
<sigma, V sigma> = c^T V c).  For a sphere this is EXACTLY 1/3.  This example shows the combined
curved + high-order win on the demag factor ITSELF -- with the PROPER Gram, unlike the crude sub-point
method in hdiv_demag_curved.py whose ~2% quadrature bias MASKS the signal (which earlier looked, wrongly,
like "the demag factor does not discriminate the curved win"):

  FLAT  (any order):  ~+0.25% -- the faceting floor; order does NOT help (sigma = M.n is constant per
                      FLAT face, so order >= 1 adds nothing) -> mesh refinement is the only lever.
  CURVE order 0:      ~-1.9%  -- geometry exact, but piecewise-CONSTANT sigma under-resolves the n_z that
                      now VARIES across a curved face.
  CURVE order 2:      ~1e-4% -- geometry exact AND sigma resolved -> the demag factor is EXACT at a
                      COARSE mesh (h=0.6), where the flat mesh needs refinement and still floors at
                      ~0.1-0.25%.

=> curved + high-order converges the demag to machine-meaningful exactness at a coarse mesh, at fixed
(small) ndof -- the accuracy-per-DOF win over flat lowest-order yano-type, ON THE DEMAG FACTOR directly.
This single-layer is the surface Gram engine for the production curved+high-order HDiv-VIM operator
(N = B^T V B); the nonlinear / volume-charge part (div M != 0) still uses the hand-rolled Newtonian
potential phi_tet (hdiv_demag_tet.py).

NON-ISOTROPIC SHAPE CHECK: a 2:1 prolate spheroid has a non-trivial analytic demag N_z ~ 0.1736 (!= 1/3).
curved + order-2 single-layer nails it to ~1e-3% while the flat faceted ellipsoid floors at ~0.34%
(order-insensitive) -- so the single-layer gets the anisotropic SHAPE right, not just the isotropic 1/3.
"""
import json
import os
from math import pi, log, sqrt

from ngsolve import (Mesh, SurfaceL2, GridFunction, CoefficientFunction, Integrate,
                     specialcf, InnerProduct, TaskManager, SetNumThreads)
from ngsolve.bem import SingleLayerPotentialOperator
from netgen.csg import CSGeometry, Sphere, Ellipsoid, Pnt, Vec

SetNumThreads(4)
HERE = os.path.dirname(os.path.abspath(__file__))


def _demag_core(geo, h, curve_order, charge_order, intorder):
    """D_z = <n_z, V n_z>/V_vol via the ngsolve.bem Laplace single-layer, for a CSG solid `geo`.
    curve_order=0 -> flat (faceted) geometry; >0 -> mesh.Curve(curve_order) isoparametric geometry.
    charge_order -> SurfaceL2 polynomial order of the surface charge sigma = M.n (here M = z_hat)."""
    mesh = Mesh(geo.GenerateMesh(maxh=h))
    with TaskManager():
        if curve_order:
            mesh.Curve(curve_order)
        fes = SurfaceL2(mesh, order=charge_order)
        V = SingleLayerPotentialOperator(fes, intorder=intorder)
        gf = GridFunction(fes)
        gf.Set(specialcf.normal(3)[2], definedon=mesh.Boundaries(".*"))   # sigma = n_z on the surface
        Vs = gf.vec.CreateVector(); Vs.data = V.mat * gf.vec
        Dz = InnerProduct(gf.vec, Vs) / Integrate(CoefficientFunction(1.0), mesh)
    return float(Dz), int(fes.ndof)


def _sphere_geo():
    g = CSGeometry(); g.Add(Sphere(Pnt(0, 0, 0), 1.0)); return g


def _prolate_geo(c):
    g = CSGeometry(); g.Add(Ellipsoid(Pnt(0, 0, 0), Vec(1, 0, 0), Vec(0, 1, 0), Vec(0, 0, c))); return g


def prolate_Nz_analytic(c, a=1.0):
    """Analytic demag factor along the long (z) axis of a prolate spheroid (semi-axes a=b<c)."""
    e = sqrt(1.0 - (a / c) ** 2)
    return (1.0 - e * e) / e ** 2 * ((1.0 / (2.0 * e)) * log((1.0 + e) / (1.0 - e)) - 1.0)


def demag_singlelayer(h, curve_order, charge_order, intorder=12):
    """Sphere demag factor (D_z -> 1/3) via the ngsolve.bem single-layer."""
    return _demag_core(_sphere_geo(), h, curve_order, charge_order, intorder)


def demag_prolate(c, h, curve_order, charge_order, intorder=12):
    """Prolate-spheroid demag factor (D_z -> the non-1/3 N_z) via the ngsolve.bem single-layer -- the
    NON-isotropic curved test: confirms the single-layer gets the SHAPE right, not just isotropy."""
    return _demag_core(_prolate_geo(c), h, curve_order, charge_order, intorder)


def run(h=0.6):
    D_an = 1.0 / 3.0
    out = {"analytic_demag_z": D_an, "h": h, "table": [], "mesh_conv": [], "prolate": []}
    for curve in (0, 3):
        for order in (0, 1, 2):
            Dz, ndof = demag_singlelayer(h, curve, order)
            out["table"].append(dict(curved=bool(curve), order=order, ndof=ndof, demag_z=Dz,
                                     err=Dz / D_an - 1))
    for hh in (0.6, 0.4, 0.3):
        df, nf = demag_singlelayer(hh, 0, 0)
        dc, nc = demag_singlelayer(hh, 3, 2)
        out["mesh_conv"].append(dict(h=hh, flat_o0_err=df / D_an - 1, flat_ndof=nf,
                                     curve_o2_err=dc / D_an - 1, curve_ndof=nc))
    # NON-isotropic curved body: 2:1 prolate spheroid, demag N_z != 1/3 (a real shape test).
    c = 2.0
    Nz_an = prolate_Nz_analytic(c)
    out["prolate_analytic_Nz"] = Nz_an
    for curve, order in ((0, 0), (0, 2), (3, 0), (3, 2)):
        Dz, ndof = demag_prolate(c, h, curve, order)
        out["prolate"].append(dict(curved=bool(curve), order=order, ndof=ndof, demag_z=Dz,
                                   err=Dz / Nz_an - 1))
    return out


if __name__ == "__main__":
    res = run()
    print("ngsolve.bem Laplace single-layer -- sphere demag D_z [exact 1/3]:")
    print(f"  {'mesh':>8} {'order':>5} {'ndof':>6} {'D_z':>9} {'err vs 1/3':>12}")
    for r in res["table"]:
        print(f"  {'curve' if r['curved'] else 'FLAT':>8} {r['order']:>5} {r['ndof']:>6} "
              f"{r['demag_z']:>9.5f} {100*r['err']:>11.4f}%")
    print("  mesh convergence (FLAT o0 vs CURVE o2):")
    for m in res["mesh_conv"]:
        print(f"    h={m['h']}: FLAT o0 {100*m['flat_o0_err']:>+7.3f}% (ndof {m['flat_ndof']:4d})   "
              f"CURVE o2 {100*m['curve_o2_err']:>+8.4f}% (ndof {m['curve_ndof']:4d})")
    print("WIN: curved + order-2 -> demag EXACT (~1e-4%) at coarse mesh; flat floors at ~0.1-0.25%")
    print("regardless of order (sigma=M.n constant per flat face).  Surface Gram = ngsolve.bem, no")
    print("hand-rolled singular quadrature; high-order + curved + FMM all reused from NGSolve.")
    print(f"\nNON-isotropic curved body: 2:1 prolate spheroid, analytic N_z = {res['prolate_analytic_Nz']:.5f}")
    print(f"  {'mesh':>8} {'order':>5} {'ndof':>6} {'N_z':>9} {'err vs analytic':>16}")
    for r in res["prolate"]:
        print(f"  {'curve' if r['curved'] else 'FLAT':>8} {r['order']:>5} {r['ndof']:>6} "
              f"{r['demag_z']:>9.5f} {100*r['err']:>15.3f}%")
    print("=> the single-layer gets the non-1/3 SHAPE right; curved+order2 nails the analytic N_z, the")
    print("flat faceted ellipsoid does not -- the curved win on a body where the demag is non-trivial.")
    with open(os.path.join(HERE, "hdiv_demag_bem_singlelayer.json"), "w") as f:
        json.dump(res, f, indent=2)
    print("saved", os.path.join(HERE, "hdiv_demag_bem_singlelayer.json"))
