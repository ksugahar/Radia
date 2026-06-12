"""3D Clebsch-Legendre (Variant 2) vacuum field-geometry solver, in NGSolve.

Solves the FIELD-region inverse problem of the partial Legendre transform (Sec. 6
of README.md): independent variables (x, y, psi), unknown functions

    z(x, y, psi)    -- the height/shape of the flux surface psi = const
    phi(x, y, psi)  -- the in-surface potential

for a CURRENT-FREE (vacuum, mu0) field B = grad(phi) x grad(psi), by minimizing
the magnetic energy expressed in the (x, y, psi) coordinates:

    W[z, phi] = (1/2 mu0) Int ( phi_x^2 + phi_y^2 + D^2 ) / z_psi  dx dy dpsi,
                D = z_x phi_y - z_y phi_x.

Reconstructed field:  B = (1/z_psi) ( phi_y, -phi_x, D ).  (phi_psi does NOT enter.)

Because Clebsch fields are solenoidal by construction, the energy minimiser among
them (with the boundary data fixed) is the VACUUM field (curl-free) -- NOT merely a
force-free field; this is confirmed numerically below (the exact vacuum solutions
are recovered to machine precision / at the FE convergence rate).

This is the 3-D field-side dual of the stream-function current design; it is the
distinctive piece of this formulation (the 2-D case is Tampere CEFC 2026, the 3-D
extension is their announced-but-unpublished program -- see LITERATURE.md).

Verification: NGSolve coordinate `z` is reused as the independent `psi`. Box
x in [0,1], y in [2,3], psi in [1,2] (y-x>0, psi>0 -> non-degenerate z_psi>0).
Manufactured vacuum solutions (sympy-checked in verify_clebsch_partial_legendre_transform.py):
  uniform    B=(B0,0,0):  z = psi,                       phi = B0 y          (polynomial)
  hyperbolic B=(y,x,B0):  z = B0 log(2 psi/(B0(y-x))),   phi = B0(y^2-x^2)/(2 psi)

Run:  python solve_clebsch_legendre_3d.py
"""
import math

from ngsolve import (Mesh, H1, GridFunction, grad, dx, Integrate, sqrt, log,
                     x, y, z, CF, BilinearForm, LinearForm, Variation,
                     InnerProduct, TaskManager)
from ngsolve.solvers import Newton
from netgen.csg import CSGeometry, OrthoBrick, Pnt

B0 = 1.0
X0, X1, Y0, Y1, P0, P1 = 0.0, 1.0, 2.0, 3.0, 1.0, 2.0
psi = z                                   # the third box coordinate IS psi


def exact(case):
    """Manufactured vacuum solution (z_field, phi) and its target field B."""
    if case == "uniform":
        return psi, B0 * y, CF((B0, 0, 0))                 # z=psi, phi=B0 y
    ze = B0 * log(2 * psi / (B0 * (y - x)))
    pe = B0 * (y * y - x * x) / (2 * psi)
    return ze, pe, CF((y, x, B0))


def _make_mesh(maxh):
    geo = CSGeometry()
    geo.Add(OrthoBrick(Pnt(X0, Y0, P0), Pnt(X1, Y1, P1)))
    return Mesh(geo.GenerateMesh(maxh=maxh))


def _harmonic_init(mesh, fes, bc_cf):
    """Init = harmonic extension of the exact boundary data (uses ONLY the
    boundary; the interior is NOT the exact solution, so Newton has real work)."""
    gf = GridFunction(fes)
    gf.Set(bc_cf, definedon=mesh.Boundaries(".*"))
    u, v = fes.TnT()
    a = BilinearForm(InnerProduct(grad(u), grad(v)) * dx, symmetric=True)
    a.Assemble()
    f = LinearForm(fes); f.Assemble()
    r = gf.vec.CreateVector()
    r.data = f.vec - a.mat * gf.vec
    gf.vec.data += a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * r
    return gf


def solve(case, maxh=0.18, order=3):
    """Minimise W over (z, phi) with Dirichlet data from the exact solution;
    return L2 errors and the reconstructed-field error."""
    mesh = _make_mesh(maxh)
    ze, pe, Btar = exact(case)

    fes = H1(mesh, order=order, dirichlet=".*") * H1(mesh, order=order, dirichlet=".*")
    gfu = GridFunction(fes)
    gfu.components[0].vec.data = _harmonic_init(mesh, fes.components[0], ze).vec
    gfu.components[1].vec.data = _harmonic_init(mesh, fes.components[1], pe).vec

    u = fes.TrialFunction()
    zc, ph = u[0], u[1]
    zx, zy, zp = grad(zc)[0], grad(zc)[1], grad(zc)[2]
    px, py = grad(ph)[0], grad(ph)[1]
    D = zx * py - zy * px
    a = BilinearForm(fes, symmetric=False)
    a += Variation(((px * px + py * py + D * D) / zp) * dx)
    Newton(a, gfu, maxit=40, dampfactor=0.4, printing=False)

    gz, gp = gfu.components
    errz = sqrt(Integrate((gz - ze) ** 2, mesh) / Integrate(ze ** 2, mesh))
    errp = sqrt(Integrate((gp - pe) ** 2, mesh) / Integrate(pe ** 2, mesh))
    gzx, gzy, gzp = grad(gz)[0], grad(gz)[1], grad(gz)[2]
    gpx, gpy = grad(gp)[0], grad(gp)[1]
    Bre = CF((gpy / gzp, -gpx / gzp, (gzx * gpy - gzy * gpx) / gzp))
    Berr = sqrt(Integrate(InnerProduct(Bre - Btar, Bre - Btar), mesh) /
                Integrate(InnerProduct(Btar, Btar), mesh))
    return {"case": case, "maxh": maxh, "ndof": fes.ndof,
            "errz": float(errz), "errp": float(errp), "Berr": float(Berr)}


def main():
    print("=== 3D Clebsch-Legendre (Variant 2) vacuum solver ===")
    with TaskManager():
        uni = solve("uniform")
        hyp = solve("hyperbolic")
        conv = [solve("hyperbolic", maxh=h) for h in (0.30, 0.20, 0.13)]

    print(f"uniform    : errz={uni['errz']:.2e} errp={uni['errp']:.2e} "
          f"Berr={uni['Berr']:.2e}  (polynomial -> machine precision)")
    print(f"hyperbolic : errz={hyp['errz']:.2e} errp={hyp['errp']:.2e} "
          f"Berr={hyp['Berr']:.2e}")
    print("h-refinement (hyperbolic): errz should drop monotonically")
    for c in conv:
        print(f"   maxh={c['maxh']:.2f} ndof={c['ndof']:6d} errz={c['errz']:.2e} "
              f"Berr={c['Berr']:.2e}")
    errs = [c["errz"] for c in conv]
    rate = math.log(errs[0] / errs[-1]) / math.log(0.30 / 0.13)
    print(f"observed convergence order (errz) ~ {rate:.1f}")
    print("  => minimising W recovers the exact VACUUM solutions (machine precision")
    print("     for the polynomial case, FE-convergent for the non-polynomial case)")
    print("     -> energy-min on Clebsch fields IS the vacuum field, not force-free.")

    ok = (uni["errz"] < 1e-10 and uni["errp"] < 1e-10 and uni["Berr"] < 1e-8
          and hyp["errz"] < 1e-4 and hyp["Berr"] < 1e-3
          and errs[0] > errs[1] > errs[2])      # monotone h-convergence
    result = {"uniform": uni, "hyperbolic": hyp, "convergence": conv,
              "conv_order": float(rate)}
    print(f"\n=== PASS: {ok} ===")
    return result, ok


if __name__ == "__main__":
    _, ok = main()
    raise SystemExit(0 if ok else 1)
