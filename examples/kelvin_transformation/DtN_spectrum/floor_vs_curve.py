# Prove the 5-6 digit floor (at p>=n) is GEOMETRY: hold FE order p>=n and the
# coarse mesh fixed, raise ONLY the isoparametric (Curve) order k of the sphere.
# If the floor is curved-geometry error, it must DROP as k rises. If it were a
# multipole/method limit it would not move.
import ngsolve as ng
from netgen.occ import Sphere, Pnt, OCCGeometry
from radia_mcp.radia_ngsolve.fem_bem_coupling import _solid_harmonic
R = 1.0

def lam_err(n, p, maxh, curve_k):
    mesh = ng.Mesh(OCCGeometry(Sphere(Pnt(0, 0, 0), R)).GenerateMesh(maxh=maxh))
    mesh.Curve(curve_k)                       # isoparametric GEOMETRY order = k (the only knob)
    fes = ng.H1(mesh, order=p, dirichlet=".*"); u, v = fes.TnT()
    a = ng.BilinearForm(ng.grad(u) * ng.grad(v) * ng.dx); a.Assemble()
    gf = ng.GridFunction(fes); gf.Set(_solid_harmonic(n), ng.BND)
    r = gf.vec.CreateVector(); r.data = -(a.mat * gf.vec)
    gf.vec.data += a.mat.Inverse(fes.FreeDofs()) * r
    bi = 2 * (p + curve_k) + 6                 # integrate richly so the floor is geometry, not quadrature
    energy = float(ng.Integrate(ng.grad(gf) * ng.grad(gf) * ng.dx(bonus_intorder=bi), mesh))
    bmass = float(ng.Integrate(gf * gf * ng.ds(bonus_intorder=bi), mesh))
    lam = -1.0 / R - energy / bmass; lam_ex = -(n + 1) / R
    return abs(lam - lam_ex) / abs(lam_ex)

results = {}
for n, p in [(2, 3), (3, 4)]:
    print(f"n={n}, p={p} (>=n), fixed coarse maxh=0.5; sweep isoparametric Curve order k:")
    ks, es = [], []
    for k in [1, 2, 3, 4, 5]:
        e = lam_err(n, p, 0.5, k)
        ks.append(k); es.append(e)
        print(f"  Curve k={k}: rel_err={e:.3e}")
    results[(n, p)] = (ks, es)
    print()

import json
print("JSON " + json.dumps({f"n{n}_p{p}": {"k": ks, "err": es} for (n, p), (ks, es) in results.items()}))
