"""Hodograph + Kelvin: axisymmetric sphere with an EXACT open boundary.

RESEARCH example (NOT a panel mode).  The verified forward "Clebsch
hodograph" panel mode (radia-electromagnet) self-meshes a sphere with a
FAR-TRUNCATED Dirichlet boundary (r/a = 6), which leaves ~3e-3 field
error.  This script replaces that truncation with the lab's verified
axisymmetric Kelvin transformation (open boundary), so the field is exact
and the hodograph is computed on an exact field.

It reuses the verified Omega-Reduced-Omega + Kelvin solve from
``docs/kelvin/legacy_assets/kelvin_transformation/Omega_ReducedOmega/Sphere/`` (interior
physical half-disk + exterior Kelvin half-disk offset in z, periodic
kelvin_int <-> kelvin_ext, mu_ext = mu0 (R/rho')^2) and ADDS the hodograph
post-processing on the physical region (air_inner):

  scalar potential Phi  : ALGEBRAIC,  Phi = H0 z + Omega_r  (H = grad Phi);
  flux function psi      : axisym Stokes,  grad(psi) = r (B_z, -B_r),
                           psi = 0 on the axis (exact: zero flux on axis);
  B(from psi) = (-psi_z/r, psi_r/r),   B(from Phi) = mu0 grad(Phi).

Verified (mu_r = 100, sphere R = 0.5, Kelvin R = 1.0, order 3, maxh 0.03):
  interior Hz = 3/(mu_r+2) H0 = 0.029412,  field_error ~1e-7
      (vs ~3e-3 for the truncated panel-mode sphere -- the Kelvin open
       boundary is the win);
  hodograph consistency B(psi) vs B(Phi) ~2e-5  (air_inner, off-axis).

The full 3-D inverse pole-design (geometry-as-unknown hodograph PDE) is a
separate research line and is NOT here.
"""
import math
import os

from numpy import pi, sqrt as npsqrt
from ngsolve import (Mesh, H1, Periodic, GridFunction, grad, InnerProduct,
                     dx, ds, x, y, CF, IfPos, BilinearForm, LinearForm,
                     TaskManager, Integrate, VOL, BND, specialcf)
from netgen.occ import (Circle, MoveTo, Glue, Vertex, Pnt, OCCGeometry,
                        IdentificationType)
from radia.kelvin_source import kelvin_mu_factor_axisym_cf, build_material_cf


def _geometry(sphere_r, kelvin_r, offset_z):
    cut_int = MoveTo(-kelvin_r - 1, -kelvin_r - 1).Rectangle(
        kelvin_r + 1, 2 * kelvin_r + 2).Face()
    outer_int = Circle((0, 0), kelvin_r).Face() - cut_int
    inner_int = Circle((0, 0), sphere_r).Face() - cut_int
    air_inner = outer_int - inner_int
    for e in air_inner.edges:
        if e.center.x < 1e-6:
            e.name = "axis_int"
        elif abs(npsqrt(e.center.x**2 + e.center.y**2) - kelvin_r) < kelvin_r * 0.2:
            e.name = "kelvin_int"
        else:
            e.name = "sphere"
    air_inner.faces.name = "air_inner"
    for e in inner_int.edges:
        e.name = "axis_int" if e.center.x < 1e-6 else "sphere"
    inner_int.faces.name = "magnetic"
    cut_ext = MoveTo(-kelvin_r - 1, offset_z - kelvin_r - 1).Rectangle(
        kelvin_r + 1, 2 * kelvin_r + 2).Face()
    outer_ext = Circle((0, offset_z), kelvin_r).Face() - cut_ext
    for e in outer_ext.edges:
        e.name = "axis_ext" if e.center.x < 1e-6 else "kelvin_ext"
    outer_ext.faces.name = "air_outer"
    gnd = Vertex(Pnt(0, offset_z, 0)); gnd.name = "GND"
    shape = Glue([air_inner, inner_int, outer_ext, gnd])
    ki = [e for e in shape.edges if e.name == "kelvin_int"]
    ke = [e for e in shape.edges if e.name == "kelvin_ext"]
    for ie in ki:
        for ee in ke:
            if (ie.center.y > 0) == (ee.center.y - offset_z > 0):
                ie.Identify(ee, "periodic", IdentificationType.PERIODIC)
                break
    return OCCGeometry(shape, dim=2)


def solve(mu_r=100.0, order=3, maxh=0.03, sphere_r=0.5, kelvin_r=1.0,
          offset_z=3.0, H0=1.0, x_eval=0.12, plot=False):
    """Solve the Kelvin sphere + hodograph and return the verification dict."""
    mu0 = 4 * pi * 1e-7
    geo = _geometry(sphere_r, kelvin_r, offset_z)
    with TaskManager():
        mesh = Mesh(geo.GenerateMesh(maxh=maxh, grading=0.7))
        mesh.Curve(order)
        mu_fac = kelvin_mu_factor_axisym_cf(z_offset=offset_z, R=kelvin_r,
                                            r_coord=x, z_coord=y)
        Mu = build_material_cf(mesh, mu0, mu_fac, outer_keyword="air_outer",
                               overrides={"magnetic": mu_r * mu0})

        fes = Periodic(H1(mesh, order=order, dirichlet_bbnd="GND"))
        Om, ps = fes.TnT()
        rw = IfPos(x - 1e-10, x, 1e-10)
        Omega_s = H0 * y
        Hs = CF((0.0, H0))
        Bs = CF((0.0, mu0 * H0))

        a = BilinearForm(fes)
        a += Mu * InnerProduct(grad(Om), grad(ps)) * rw * dx
        a.Assemble()
        gfOm = GridFunction(fes, name="Omega")
        gfOm.Set(Omega_s, BND, mesh.Boundaries("sphere"))
        f = LinearForm(fes)
        f += Mu * InnerProduct(grad(gfOm), grad(ps)) * rw * dx("air_inner")
        nrm = specialcf.normal(mesh.dim)
        f += InnerProduct(nrm, Bs) * ps * rw * ds("sphere")
        f.Assemble()
        gfOm.vec.data = a.mat.Inverse(fes.FreeDofs(),
                                      inverse="sparsecholesky") * f.vec

        # air field (reduced region): H = grad(Phi), Phi = H0 z + Omega_r
        fesA = H1(mesh, order=order, definedon="air_inner|air_outer")
        Orr = GridFunction(fesA); Orr.Set(gfOm, VOL, definedon="air_inner|air_outer")
        Oxr = GridFunction(fesA); Oxr.Set(Omega_s, BND, mesh.Boundaries("sphere"))
        H_air = Hs + grad(Orr) - grad(Oxr)
        B_air = mu0 * H_air

        ironmask = mesh.MaterialCF({"magnetic": 1.0, "air_inner": 0.0,
                                    "air_outer": 0.0}, default=0.0)
        Hz_in = Integrate(ironmask * grad(gfOm)[1], mesh) / Integrate(ironmask, mesh)
        Hz_analytic = 3.0 / (mu_r + 2.0) * H0
        field_error = abs(Hz_in - Hz_analytic) / abs(Hz_analytic)

        # hodograph: solve the Stokes flux psi on air_inner (psi=0 on axis)
        fesP = H1(mesh, order=order, definedon="air_inner", dirichlet="axis_int")
        uP, vP = fesP.TnT()
        aP = BilinearForm(InnerProduct(grad(uP), grad(vP)) * dx("air_inner"),
                          symmetric=True)
        RB = CF((x * B_air[1], -x * B_air[0]))     # r (B_z, -B_r) = grad(psi)
        fP = LinearForm(InnerProduct(RB, grad(vP)) * dx("air_inner"))
        aP.Assemble(); fP.Assemble()
        gP = GridFunction(fesP, name="psi")
        gP.vec.data = aP.mat.Inverse(fesP.FreeDofs(),
                                     inverse="sparsecholesky") * fP.vec
        dP = grad(gP)
        B_from_psi = CF((-dP[1] / x, dP[0] / x))

        airmask = mesh.MaterialCF({"magnetic": 0.0, "air_inner": 1.0,
                                   "air_outer": 0.0}, default=0.0)
        offax = airmask * IfPos(x - x_eval, 1.0, 0.0)
        den = Integrate(offax * InnerProduct(B_air, B_air), mesh)
        consistency = math.sqrt(
            Integrate(offax * InnerProduct(B_from_psi - B_air,
                                           B_from_psi - B_air), mesh) / den)

        if plot:
            _plot_flux(mesh, B_air, gfOm, sphere_r, kelvin_r)

    return {
        "mu_r": mu_r, "order": order, "ne": mesh.ne,
        "Hz_in": float(Hz_in), "Hz_analytic": float(Hz_analytic),
        "field_error": float(field_error),
        "consistency": float(consistency),
    }


def _plot_flux(mesh, B_air, gfOm, sphere_r, kelvin_r):
    """Save flux streamlines of the exact (Kelvin) field in the physical
    meridional half-plane next to this script."""
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    n = 60
    rs = np.linspace(0.02, kelvin_r - 0.02, n)
    zs = np.linspace(-kelvin_r + 0.02, kelvin_r - 0.02, n)
    RR, ZZ = np.meshgrid(rs, zs)
    Br = np.full(RR.shape, np.nan)
    Bz = np.full(RR.shape, np.nan)
    for i in range(n):
        for j in range(n):
            if math.hypot(rs[j], zs[i]) < kelvin_r - 0.02:
                mip = mesh(rs[j], zs[i])
                if math.hypot(rs[j], zs[i]) < sphere_r:
                    g = grad(gfOm)
                    Br[i, j] = g[0](mip); Bz[i, j] = g[1](mip)
                else:
                    b = B_air(mip)
                    Br[i, j] = b[0]; Bz[i, j] = b[1]
    fig, ax = plt.subplots(figsize=(4.6, 5.2), dpi=150)
    ax.streamplot(RR, ZZ, Br, Bz, color="0.2", linewidth=0.7, density=1.3,
                  arrowsize=0.7)
    th = np.linspace(0, math.pi, 100)
    ax.fill_betweenx(sphere_r * np.cos(th), 0, sphere_r * np.sin(th),
                     color="lightblue", alpha=0.5)
    ax.plot(kelvin_r * np.sin(th), kelvin_r * np.cos(th), "g--", lw=1.2)
    ax.set_aspect("equal")
    ax.set_xlabel("r")
    ax.set_ylabel("z")
    ax.set_title("Kelvin (exact open bdry) flux lines")
    png = os.path.splitext(os.path.abspath(__file__))[0] + ".png"
    fig.tight_layout()
    fig.savefig(png, bbox_inches="tight")
    plt.close(fig)
    print(f"  figure saved: {png}")


def main():
    res = solve(plot=True)
    print("Hodograph + Kelvin (axisymmetric sphere, exact open boundary)")
    print(f"  mu_r={res['mu_r']:.0f}  ne={res['ne']}")
    print(f"  interior Hz = {res['Hz_in']:.6f}  analytic = "
          f"{res['Hz_analytic']:.6f}  field_error = {res['field_error']:.2e}")
    print(f"  hodograph consistency B(psi) vs B(Phi) = {res['consistency']:.2e}")
    print("  (field_error ~1e-7 with Kelvin vs ~3e-3 truncated -- exact "
          "open boundary)")


if __name__ == "__main__":
    main()
