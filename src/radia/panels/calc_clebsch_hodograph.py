"""calc_clebsch_hodograph.py -- Clebsch hodograph (reduced-potential) forward mode.

CLI entry point for the radia-electromagnet panel's "Clebsch hodograph" mode.

Forward, verified mode
----------------------
Solve the reduced scalar potential Omega for a magnetizable body in a
uniform applied field (the same physics the Kelvin Benchmark mode uses),
then post-process the field into its CONJUGATE Clebsch potentials:

  - the flux function   psi  (the Stokes / 2-D stream function),
  - the scalar potential Phi (= H0 . x - Omega, since H = H_s - grad Omega
    is curl-free in the current-free region, so H = grad Phi).

These are de Rham / Hodge conjugates: flux-lines are level sets of psi,
pole equipotentials are level sets of Phi.  The headline output is the

  BIDIRECTIONAL CONSISTENCY = || B(from psi) - B(from Phi) || / || B ||

evaluated on the SAME FEM field -- a self-consistency (a-posteriori field
quality) metric that does NOT need the exact solution and tends to zero as
the FEM field becomes exact.  A separate field_error compares the FEM field
against the known analytical solution.

Scope (honest)
--------------
2-D / axisymmetric, where the flux function psi is well defined.  Verified
against exact solutions:

  cylinder (2-D)    mu_r disk in uniform field, 2-D dipole reaction
                    -> consistency ~ 7.8e-4   (path B)
  sphere   (axisym) mu_r sphere in uniform field, induced-dipole reaction
                    -> consistency ~ 1e-3     (path A, axis-anchored psi)

The full 3-D INVERSE pole-design (geometry-as-unknown hodograph PDE,
"redraw the geometry in potential coordinates") is a separate research
line (branch clebsch_legendre); it is NOT this forward mode and is NOT
shipped in the panel -- per the repo-first "Don't Publish the Unfinished"
policy it stays on the branch until verified on a real device geometry.

Layer 4 contract (CLAUDE.md): NO `import cubit`, NO PySide6.  Heavy imports
(ngsolve / netgen / radia) live inside run() / the _solve_* helpers, never
at module top, so the panel can introspect build_argparser() cheaply.

References
----------
- Clebsch / Euler potentials: B = grad(phi) x grad(psi); A = phi grad(psi).
- Dervisha, Marjamaki, Rasilo, Tarhasaari, "Bidirectional Coordinate
  Transformation and Its Application to 2-D Magnetic Field Problems",
  CEFC 2026 (the 2-D "potential coordinate" precedent; this forward mode
  is the verified NGSolve realisation, not a novelty claim).
"""

import argparse
import math
import os
import sys

_this_dir = os.path.dirname(os.path.abspath(__file__))
if _this_dir not in sys.path:
    sys.path.insert(0, _this_dir)

from calc_common import progress, calc_main


def _log(msg):
    progress("CLEBSCH", msg)


# ============================================================
# argparse surface (single source of truth for the panel widgets)
# ============================================================

def build_argparser() -> argparse.ArgumentParser:
    """Return the argparse.  Importable WITHOUT side effects (no NGSolve /
    Radia / Cubit imports) so the panel can introspect the widgets."""
    p = argparse.ArgumentParser(
        description="Clebsch hodograph forward mode: reduced-Omega solve, "
                    "then conjugate flux psi / scalar Phi potentials + "
                    "bidirectional consistency (a-posteriori field quality) "
                    "+ flux-line visualization.  Scope: 2-D / axisymmetric, "
                    "verified on cylinder / sphere.")
    p.add_argument("--geometry", choices=["cylinder", "sphere"],
                   default="cylinder",
                   help="Self-meshed verified geometry: 'cylinder' (2-D disk "
                        "in uniform field) or 'sphere' (axisymmetric).")
    p.add_argument("--mu-r", type=float, default=1000.0,
                   help="Relative permeability of the magnetizable body.")
    p.add_argument("--maxh", type=float, default=0.08,
                   help="Mesh size (body radius a = 1 in normalized units).")
    p.add_argument("--fes-order", type=int, default=3,
                   help="H1 polynomial order (>= 2 recommended).")
    p.add_argument("--msh-output", default=None,
                   help="GMSH .msh path for flux-line / equipotential viz.")
    p.add_argument("--output", default=None,
                   help="Result JSON path (calc_main writes the return dict).")
    return p


# ============================================================
# run() -- heavy imports happen here, not at module top
# ============================================================

def run(args) -> dict:
    import time
    from ngsolve import TaskManager   # first ngsolve import of the function

    t0 = time.perf_counter()
    with TaskManager():
        if args.geometry == "cylinder":
            result = _solve_cylinder(args)
        elif args.geometry == "sphere":
            result = _solve_sphere(args)
        else:
            raise ValueError(f"unknown geometry {args.geometry!r}")
    result["t_total_s"] = round(time.perf_counter() - t0, 3)
    return result


# ============================================================
# cylinder (2-D) -- genuine reduced-Omega FEM solve + hodograph (path B)
# ============================================================

def _solve_cylinder(args) -> dict:
    import time
    from ngsolve import (Mesh, H1, GridFunction, grad, InnerProduct, dx,
                         x, y, CF, Norm, BilinearForm, LinearForm, Integrate)
    from netgen.geom2d import CSG2d, Circle

    mu_r = float(args.mu_r)
    order = int(args.fes_order)
    maxh = float(args.maxh)
    a_cyl, r_out, H0 = 1.0, 6.0, 1.0
    beta2 = (mu_r - 1.0) / (mu_r + 1.0)        # 2-D induced-dipole strength

    r2 = x * x + y * y
    phi_red_out = beta2 * a_cyl ** 2 * x / r2  # exact reduced potential (air)
    phi_red_in = beta2 * x                     # exact reduced potential (iron)
    psi_air_exact = y * (1.0 + beta2 * a_cyl ** 2 / r2)   # exact 2-D flux (air)

    # --- mesh: iron disk in an air annulus -----------------------------
    t_mesh0 = time.perf_counter()
    geo = CSG2d()
    iron = Circle(center=(0, 0), radius=a_cyl, bc="iface").Mat("iron")
    air = (Circle(center=(0, 0), radius=r_out, bc="outer").Mat("air") - iron)
    geo.Add(iron)
    geo.Add(air)
    mesh = Mesh(geo.GenerateMesh(maxh=maxh))
    t_mesh_s = time.perf_counter() - t_mesh0
    _log(f"cylinder mesh: {mesh.ne} elements, mu_r={mu_r:.0f}")

    mu = mesh.MaterialCF({"iron": mu_r, "air": 1.0}, default=1.0)
    Hs = CF((H0, 0.0))                          # uniform applied source

    # --- reduced-Omega FEM solve:  int mu grad(w).grad(v) = int mu Hs.grad(v)
    t_solve0 = time.perf_counter()
    fes = H1(mesh, order=order, dirichlet="outer")
    u, v = fes.TnT()
    a = BilinearForm(mu * InnerProduct(grad(u), grad(v)) * dx, symmetric=True)
    f = LinearForm(mu * InnerProduct(Hs, grad(v)) * dx)
    a.Assemble()
    f.Assemble()
    gf = GridFunction(fes)
    gf.Set(CF(phi_red_out), definedon=mesh.Boundaries("outer"))
    rr = gf.vec.CreateVector()
    rr.data = f.vec - a.mat * gf.vec
    gf.vec.data += a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * rr
    t_solve_s = time.perf_counter() - t_solve0

    airmask = mesh.MaterialCF({"iron": 0.0, "air": 1.0}, default=0.0)
    H_fem = Hs - grad(gf)                       # B = H here (mu0 = 1)

    # --- accuracy vs the exact reduced potential / field ---------------
    phi_ex = mesh.MaterialCF({"iron": phi_red_in, "air": phi_red_out},
                             default=phi_red_out)
    err_phi = math.sqrt(Integrate((gf - phi_ex) ** 2, mesh)
                        / Integrate(phi_ex ** 2, mesh))
    H_ex = Hs - CF((phi_red_out.Diff(x), phi_red_out.Diff(y)))
    den = Integrate(airmask * InnerProduct(H_ex, H_ex), mesh)
    field_error = math.sqrt(
        Integrate(airmask * InnerProduct(H_fem - H_ex, H_fem - H_ex), mesh)
        / den)

    # --- hodograph: conjugate Clebsch potentials -----------------------
    t_post0 = time.perf_counter()
    Phi_scalar = H0 * x - gf                    # B = grad(Phi) in air
    # flux psi: grad(psi) = R90(B); DUAL BC -- Dirichlet on the clean far
    # boundary ONLY, the iron interface is NATURAL (it is ~an equipotential
    # where flux enters, so psi varies there -- dual to Phi's Dirichlet).
    fes2 = H1(mesh, order=order, definedon=mesh.Materials("air"),
              dirichlet="outer")
    u2, v2 = fes2.TnT()
    a2 = BilinearForm(InnerProduct(grad(u2), grad(v2)) * dx, symmetric=True)
    RB = CF((-H_fem[1], H_fem[0]))              # R90 of B = grad(psi)
    f2 = LinearForm(InnerProduct(RB, grad(v2)) * dx)
    a2.Assemble()
    f2.Assemble()
    gpsi = GridFunction(fes2)
    gpsi.Set(CF(psi_air_exact), definedon=mesh.Boundaries("outer"))
    r2v = gpsi.vec.CreateVector()
    r2v.data = f2.vec - a2.mat * gpsi.vec
    gpsi.vec.data += a2.mat.Inverse(fes2.FreeDofs(),
                                    inverse="sparsecholesky") * r2v
    dpsi = grad(gpsi)
    B_from_psi = CF((dpsi[1], -dpsi[0]))        # (psi_y, -psi_x)
    consistency = math.sqrt(
        Integrate(airmask * InnerProduct(B_from_psi - H_fem,
                                         B_from_psi - H_fem), mesh) / den)
    t_post_s = time.perf_counter() - t_post0
    _log(f"cylinder: consistency={consistency:.2e}, "
         f"field_error={field_error:.2e}, err_phi={err_phi:.2e}")

    gmsh_file = ""
    if args.msh_output:
        gmsh_file = _emit_msh(mesh, gpsi, Phi_scalar, H_fem, args.msh_output)

    return {
        "geometry": "cylinder", "mu_r": mu_r, "fes_order": order,
        "ne": mesh.ne, "ndof": int(fes.ndof + fes2.ndof),
        "t_mesh_s": round(t_mesh_s, 3),
        "t_solve_s": round(t_solve_s, 3),
        "t_post_s": round(t_post_s, 3),
        "consistency": float(consistency),
        "field_error": float(field_error),
        "err_phi": float(err_phi),
        "gmsh_file": gmsh_file,
    }


# ============================================================
# sphere (axisymmetric) -- genuine reduced-Omega FEM solve + hodograph (path A)
# ============================================================

def _solve_sphere(args) -> dict:
    import time
    from ngsolve import (Mesh, H1, GridFunction, grad, InnerProduct, dx,
                         x, y, CF, IfPos, BilinearForm, LinearForm, Integrate)
    from netgen.geom2d import CSG2d, Circle, Rectangle

    mu_r = float(args.mu_r)
    order = int(args.fes_order)
    maxh = float(args.maxh)
    a_sph, r_out, H0 = 1.0, 6.0, 1.0
    # axisym coords: r = x (>= 0), z = y.  H_s = H0 z_hat.
    x_eval = 0.30        # evaluate consistency off-axis (1/r reconstruction)

    # --- mesh: right half-disk of iron in a half-annulus of air --------
    t_mesh0 = time.perf_counter()
    geo = CSG2d()
    box = Rectangle(pmin=(0, -r_out), pmax=(r_out, r_out),
                    left="axis", right="outer", top="outer", bottom="outer")
    iron = (Circle(center=(0, 0), radius=a_sph, bc="iface").Mat("iron") * box)
    air = (Rectangle(pmin=(0, -r_out), pmax=(r_out, r_out),
                     left="axis", right="outer", top="outer", bottom="outer")
           .Mat("air") - iron)
    geo.Add(iron)
    geo.Add(air)
    mesh = Mesh(geo.GenerateMesh(maxh=maxh))
    t_mesh_s = time.perf_counter() - t_mesh0
    _log(f"sphere (axisym) mesh: {mesh.ne} elements, mu_r={mu_r:.0f}")

    mu = mesh.MaterialCF({"iron": mu_r, "air": 1.0}, default=1.0)
    Hs = CF((0.0, H0))                          # uniform axial applied source

    # --- reduced-Omega axisym FEM solve (2 pi r weight; far Dirichlet=0) -
    t_solve0 = time.perf_counter()
    fes = H1(mesh, order=order, dirichlet="outer")
    u, v = fes.TnT()
    a = BilinearForm(mu * InnerProduct(grad(u), grad(v)) * x * dx,
                     symmetric=True)
    f = LinearForm(mu * InnerProduct(Hs, grad(v)) * x * dx)
    a.Assemble()
    f.Assemble()
    gf = GridFunction(fes)                       # far-field Dirichlet phi = 0
    gf.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * f.vec
    t_solve_s = time.perf_counter() - t_solve0

    airmask = mesh.MaterialCF({"iron": 0.0, "air": 1.0}, default=0.0)
    ironmask = mesh.MaterialCF({"iron": 1.0, "air": 0.0}, default=0.0)
    H_fem = Hs - grad(gf)

    # --- accuracy: interior axial field vs analytical 3/(mu_r+2) H0 ----
    Hz_in_analytic = 3.0 / (mu_r + 2.0) * H0
    vol_iron = Integrate(ironmask, mesh)
    Hz_in = Integrate(ironmask * H_fem[1], mesh) / vol_iron
    field_error = abs(Hz_in - Hz_in_analytic) / abs(Hz_in_analytic)

    # --- hodograph: conjugate Clebsch potentials (axisym) --------------
    t_post0 = time.perf_counter()
    Phi_scalar = H0 * y - gf                    # H = grad(Phi) in air
    # axisym Stokes flux psi:  grad(psi) = r (B_z, -B_r);  B_r=-(1/r)psi_z,
    # B_z=(1/r)psi_r.  Anchor psi = 0 on the axis (exact: flux through the
    # axis is zero), natural elsewhere -- the DUAL BC in axisym.
    fes2 = H1(mesh, order=order, definedon=mesh.Materials("air"),
              dirichlet="axis")
    u2, v2 = fes2.TnT()
    a2 = BilinearForm(InnerProduct(grad(u2), grad(v2)) * dx, symmetric=True)
    RB = CF((x * H_fem[1], -x * H_fem[0]))      # r (B_z, -B_r) = grad(psi)
    f2 = LinearForm(InnerProduct(RB, grad(v2)) * dx)
    a2.Assemble()
    f2.Assemble()
    gpsi = GridFunction(fes2)                    # psi = 0 on axis
    gpsi.vec.data = a2.mat.Inverse(fes2.FreeDofs(),
                                   inverse="sparsecholesky") * f2.vec
    dpsi = grad(gpsi)
    # reconstruct B from the flux (off-axis): B_r=-(1/r)psi_z, B_z=(1/r)psi_r
    B_from_psi = CF((-dpsi[1] / x, dpsi[0] / x))
    offaxis = airmask * IfPos(x - x_eval, 1.0, 0.0)
    den = Integrate(offaxis * InnerProduct(H_fem, H_fem), mesh)
    consistency = math.sqrt(
        Integrate(offaxis * InnerProduct(B_from_psi - H_fem,
                                         B_from_psi - H_fem), mesh) / den)
    t_post_s = time.perf_counter() - t_post0
    _log(f"sphere: consistency={consistency:.2e} (x>{x_eval}), "
         f"Hz_in={Hz_in:.4f} vs {Hz_in_analytic:.4f}, "
         f"field_error={field_error:.2e}")

    gmsh_file = ""
    if args.msh_output:
        gmsh_file = _emit_msh(mesh, gpsi, Phi_scalar, H_fem, args.msh_output)

    return {
        "geometry": "sphere", "mu_r": mu_r, "fes_order": order,
        "ne": mesh.ne, "ndof": int(fes.ndof + fes2.ndof),
        "t_mesh_s": round(t_mesh_s, 3),
        "t_solve_s": round(t_solve_s, 3),
        "t_post_s": round(t_post_s, 3),
        "consistency": float(consistency),
        "field_error": float(field_error),
        "Hz_in": float(Hz_in),
        "Hz_in_analytic": float(Hz_in_analytic),
        "gmsh_file": gmsh_file,
    }


# ============================================================
# visualization
# ============================================================

def _emit_msh(mesh, gpsi, Phi_scalar, H_fem, path) -> str:
    """Write a GMSH .msh v4.1 with the flux psi, scalar Phi, and |B| fields
    so GMSH renders flux-lines (psi iso-contours) and pole equipotentials
    (Phi iso-contours).  Returns the absolute path written."""
    from ngsolve import Norm
    from radia.gmsh_post_export import GmshPostExport
    post = GmshPostExport(mesh)
    post.add_scalar_field("psi_flux", gpsi)
    post.add_scalar_field("Phi_scalar", Phi_scalar)
    post.add_scalar_field("absB", Norm(H_fem))
    post.write(path)
    return os.path.abspath(path)


def main():
    calc_main(run, build_argparser())


if __name__ == "__main__":
    main()
