# Cleanest test: FIX the Gamma surface (per-face maxh fixed), refine ONLY interior volume.
import ngsolve as ng
from netgen.occ import Sphere, Pnt, OCCGeometry
from radia_mcp.radia_ngsolve.fem_bem_coupling import _solid_harmonic
R = 1.0

def solve_lam(mesh, n, p):
    mesh.Curve(min(p + 1, 3))
    fes = ng.H1(mesh, order=p, dirichlet=".*"); u, v = fes.TnT()
    a = ng.BilinearForm(ng.grad(u) * ng.grad(v) * ng.dx); a.Assemble()
    gf = ng.GridFunction(fes); gf.Set(_solid_harmonic(n), ng.BND)
    r = gf.vec.CreateVector(); r.data = -(a.mat * gf.vec)
    gf.vec.data += a.mat.Inverse(fes.FreeDofs()) * r
    energy = float(ng.Integrate(ng.grad(gf) * ng.grad(gf) * ng.dx(bonus_intorder=10), mesh))
    bmass = float(ng.Integrate(gf * gf * ng.ds(bonus_intorder=10), mesh))
    lam = -1.0 / R - energy / bmass; lam_ex = -(n + 1) / R
    nse = sum(1 for _ in mesh.Elements(ng.BND))
    return fes.ndof, mesh.ne, nse, lam, abs(lam - lam_ex) / abs(lam_ex)

n, p = 2, 3
print(f"mode n={n}, order p={p}>=n.  Gamma surface FIXED (face maxh=0.5), refine INTERIOR volume only:\n")
for hv in [0.5, 0.3, 0.2, 0.13, 0.09]:
    sph = Sphere(Pnt(0, 0, 0), R); sph.faces.maxh = 0.5      # FIX the surface element size
    mesh = ng.Mesh(OCCGeometry(sph).GenerateMesh(maxh=hv))   # interior gets hv
    ndof, ne, nse, lam, rel = solve_lam(mesh, n, p)
    print(f"  interior hv={hv:.2f}: surf_tris={nse:5d}  vol_tets={ne:6d}  ndof={ndof:7d}  lambda={lam:+.8f}  rel_err={rel:.3e}")
