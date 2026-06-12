# Does refining the EXTERIOR (Kelvin-ball) mesh help, if Gamma is fixed?
# DtN-spectral prediction: the open-BC accuracy is a property of Gamma (trace+geometry)
# and the order p, NOT of the exterior VOLUME density. Once p>=n the volume solve is
# Galerkin-EXACT (solution is the degree-n polynomial, in V_h on ANY mesh).
import ngsolve as ng
from netgen.occ import OCCGeometry, Sphere, Pnt
from radia_mcp.radia_ngsolve.fem_bem_coupling import _solid_harmonic

R = 1.0

def solve_lam(mesh, n, p, curve):
    if curve:
        mesh.Curve(min(p + 1, 3))
    fes = ng.H1(mesh, order=p, dirichlet=".*")
    u, v = fes.TnT()
    a = ng.BilinearForm(ng.grad(u) * ng.grad(v) * ng.dx); a.Assemble()
    gf = ng.GridFunction(fes); gf.Set(_solid_harmonic(n), ng.BND)
    r = gf.vec.CreateVector(); r.data = -(a.mat * gf.vec)
    gf.vec.data += a.mat.Inverse(fes.FreeDofs()) * r
    energy = float(ng.Integrate(ng.grad(gf) * ng.grad(gf) * ng.dx(bonus_intorder=10), mesh))
    bmass = float(ng.Integrate(gf * gf * ng.ds(bonus_intorder=10), mesh))
    lam = -1.0 / R - energy / bmass
    lam_ex = -(n + 1) / R
    # Galerkin exactness probe: ||u_h - exact polynomial||_L2 over the volume
    ef = float(ng.Integrate((gf - _solid_harmonic(n)) ** 2 * ng.dx(bonus_intorder=10), mesh)) ** 0.5
    return fes.ndof, mesh.ne, lam, abs(lam - lam_ex) / abs(lam_ex), ef

n, p = 2, 3  # quadrupole, order p=3 >= n  (so the polynomial image is in V_h)
print(f"mode n={n} (image = degree-{n} poly), order p={p} (>= n)\n")

print("=== AFFINE ball: Gamma GEOMETRY FROZEN, refine exterior mesh uniformly ===")
mesh = ng.Mesh(OCCGeometry(Sphere(Pnt(0, 0, 0), R)).GenerateMesh(maxh=0.6))
for k in range(4):
    if k > 0:
        mesh.Refine()
    ndof, ne, lam, rel, ef = solve_lam(mesh, n, p, curve=False)
    print(f"  refine {k}: ne={ne:6d} ndof={ndof:7d}  lam={lam:+.8f}  rel_err={rel:.3e}  ||u_h-poly||={ef:.2e}")

print("\n=== CURVED ball: Gamma GEOMETRY IMPROVES as the mesh refines ===")
mesh2 = ng.Mesh(OCCGeometry(Sphere(Pnt(0, 0, 0), R)).GenerateMesh(maxh=0.6))
for k in range(4):
    if k > 0:
        mesh2.Refine()
    ndof, ne, lam, rel, ef = solve_lam(mesh2, n, p, curve=True)
    print(f"  refine {k}: ne={ne:6d} ndof={ndof:7d}  lam={lam:+.8f}  rel_err={rel:.3e}  ||u_h-poly||={ef:.2e}")
