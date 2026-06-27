"""Solve Omega-Reduced Omega + Kelvin on Cubit-meshed 1/4 sector
sphere_1_4_p<N>.vol files for N = 1, 2, 3 and compare Hz at origin
to the analytical solution.

Analytical (uniform external Hz, infinite domain, sphere of mu_r at origin):
  Hz_inside = 3 / (mu_r + 2) * H0    (uniform in interior)

Formulation: Omega-Reduced Omega (cf. examples/kelvin_transformation/
Omega_ReducedOmega/Sphere/3D_sphere_with_Kelvin.py).

  Total region (magnetic):  mu * grad(Omega_t) . grad(psi) dx
  Reduced region (air):     mu * grad(Omega_r) . grad(psi) dx
                          - mu * grad(Omega_s) . grad(psi) dx
  Kelvin region (image):   (R/rho')^2 * mu * grad(Omega_t) . grad(psi) dx

Source: Omega_s = H0 * z so grad Omega_s = (0, 0, H0) => H_s = H0 z-hat.

1/4 sector BCs (x>=0, y>=0, full z):
  GND        Dirichlet (Omega = 0, kelvin centre = "infinity")
  sym_bn=0_x Natural   (Hx = 0 on x=0 plane; rotational symmetry around z)
  sym_bn=0_y Natural   (Hy = 0 on y=0 plane)
  kelvin_int <-> kelvin_ext   Periodic (couples interior r=R with image r'=R)
"""

import argparse
import json
import math
import os
import sys
import time

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(THIS_DIR, "..", "..", "..", "src", "radia",
                                 "panels"))


def solve_at_order(vol_path, order, mu_r=100.0, H0=1.0, R_kelvin=0.20):
    """Solve Omega-reduced + Kelvin and return Hz at origin."""
    if "calc_common" in sys.modules:
        del sys.modules["calc_common"]
    from calc_common import detect_kelvin_offset, add_periodic_kelvin
    from ngsolve import (Mesh, H1, Periodic, BilinearForm, LinearForm,
                         GridFunction, CoefficientFunction, BND, ds,
                         specialcf, x, y, z,
                         grad, dx, sqrt, Integrate, TaskManager)

    mu0 = 4 * math.pi * 1e-7

    print(f"\n=== order={order}: {vol_path} ===")
    t0 = time.perf_counter()
    mesh = Mesh(vol_path)
    print(f"  ne={mesh.ne}, materials={mesh.GetMaterials()}")
    print(f"  boundaries={mesh.GetBoundaries()}")

    # Add periodic identification (handles unmatched pair via fresh IDPB)
    offset = detect_kelvin_offset(mesh)
    print(f"  detect_kelvin_offset = {offset}")
    if add_periodic_kelvin(mesh, offset):
        import ngsolve
        mesh = ngsolve.Mesh(mesh.ngmesh)

    # Curve the mesh for high-order geometry
    mesh.Curve(order)

    # Build Periodic FES with combined Dirichlet:
    #   GND       -- always (Kelvin centre = "infinity")
    #   sym_ht=0_*-- antimirror plane (Omega = 0 anti-symmetric)
    #   kelvin_far-- 1/8 only: kelvin "infinity plane"
    # sym_bn=0_*  -- Natural BC for Omega (no Dirichlet entry needed)
    bnd = list(mesh.GetBoundaries())
    bbb = list(mesh.GetBBBoundaries()) if hasattr(mesh, "GetBBBoundaries") else []
    dir_parts = []
    for label in ("GND", "sym_ht=0_x", "sym_ht=0_y", "sym_ht=0_z",
                  "kelvin_far"):
        if label in bnd or label in bbb:
            dir_parts.append(label)
    dirichlet = "|".join(dir_parts) if dir_parts else ""
    print(f"  dirichlet = {dirichlet!r}")

    fes = Periodic(H1(mesh, order=order, dirichlet=dirichlet))
    free_total = sum(1 for d in fes.FreeDofs() if d)
    print(f"  ndof={fes.ndof}, free={free_total}")

    # Material coefficient functions -- use the canonical helpers from
    # radia.kelvin_source so we exactly match the reference
    # 3D_sphere_with_Kelvin.py setup.
    sys.path.insert(0, os.path.join(os.path.dirname(__file__),
                                     "..", "..", "..", "src", "radia"))
    from kelvin_source import (kelvin_mu_factor_3d_cf, build_material_cf)
    mu_kelvin_factor = kelvin_mu_factor_3d_cf(center=tuple(offset),
                                               R=R_kelvin)
    mu_cf = build_material_cf(
        mesh, mu0, mu_kelvin_factor,
        outer_keyword="kelvin",
        overrides={"magnetic": mu_r * mu0},
    )

    # === Omega-Reduced Omega weak form (matches reference
    # examples/kelvin_transformation/Omega_ReducedOmega/Sphere/
    # 3D_sphere_with_Kelvin.py).
    #
    # The "sphere" sideset (mag-air interface) is NOT a Dirichlet BC.
    # gfOmega.Set() on sphere is used as a source LIFT only -- the
    # value is captured in the LinearForm via grad(gfOmega), then
    # gfOmega is overwritten by the solve.
    #
    #   a(u, v) = int_dom  mu * grad(u) . grad(v) dx           (all materials)
    #   f(v)    = int_air  mu * grad(gfOmega_lift) . grad(v) dx
    #           + int_sphere  (-n . Bs) * v ds                 (Neumann fix)
    Omega_s = H0 * z
    Bs = CoefficientFunction((0.0, 0.0, mu0 * H0))

    u = fes.TrialFunction()
    v = fes.TestFunction()
    a = BilinearForm(fes, symmetric=True)
    a += mu_cf * grad(u) * grad(v) * dx

    # gfOmega: source lift on sphere (NOT a Dirichlet enforce; sphere
    # is not in fes.dirichlet so its DOFs are free).
    gfOmega = GridFunction(fes)
    gfOmega.vec[:] = 0
    gfOmega.Set(Omega_s, BND, mesh.Boundaries("sphere"))

    f = LinearForm(fes)
    f += mu_cf * grad(gfOmega) * grad(v) * dx("air")
    f.Assemble()

    # Remove components of the Dirichlet boundary from f.vec.
    # (Identical pattern to base class Omega_ReducedOmega.py:104-105.)
    import numpy as _np
    fcut = _np.array(f.vec.FV())[fes.FreeDofs()]
    _np.array(f.vec.FV(), copy=False)[fes.FreeDofs()] = fcut

    normal = -specialcf.normal(mesh.dim)
    f += (normal * Bs) * v * ds("sphere")
    f.Assemble()

    print(f"  assembling...")
    with TaskManager():
        a.Assemble()
        print(f"  solving (pardiso)...")
        gfOmega.vec.data = a.mat.Inverse(fes.FreeDofs(),
                                         inverse="pardiso") * f.vec
        sol = gfOmega

    # H = grad(Omega) in total region (magnetic), H = grad(Omega_r) + Hs in air
    grad_sol = grad(sol)
    # At origin, magnetic region: Hz = grad(sol)_z
    pt = mesh(1e-12, 1e-12, 1e-12)
    Hz_origin = float(grad_sol[2](pt))
    Hz_analytical = 3 / (mu_r + 2) * H0
    error_pct = (Hz_origin - Hz_analytical) / Hz_analytical * 100.0

    elapsed = time.perf_counter() - t0
    print(f"  Hz_origin = {Hz_origin:.6e}  (analytical = "
          f"{Hz_analytical:.6e}, error = {error_pct:+.3f}%)")
    print(f"  elapsed = {elapsed:.2f}s")

    return {
        "order": order,
        "ne": mesh.ne,
        "ndof": fes.ndof,
        "free": free_total,
        "Hz_origin": Hz_origin,
        "Hz_analytical": Hz_analytical,
        "error_pct": error_pct,
        "elapsed_sec": elapsed,
    }


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--in-dir", default=THIS_DIR)
    p.add_argument("--orders", default="1,2,3")
    p.add_argument("--mu-r", type=float, default=100.0)
    p.add_argument("--H0", type=float, default=1.0)
    p.add_argument("--R-kelvin", type=float, default=0.20)
    p.add_argument("--out-json",
                   default=os.path.join(THIS_DIR, "p_convergence.json"))
    args = p.parse_args()

    orders = tuple(int(s) for s in args.orders.split(","))
    results = []
    for order in orders:
        vol = os.path.join(args.in_dir, f"sphere_1_4_p{order}.vol")
        if not os.path.exists(vol):
            print(f"\nWARNING: skipping order={order}: {vol} not found")
            continue
        results.append(solve_at_order(vol, order, mu_r=args.mu_r,
                                       H0=args.H0, R_kelvin=args.R_kelvin))

    print("\n=== p-convergence summary ===")
    print(f"{'order':>6} {'ne':>8} {'ndof':>10} {'Hz_origin':>14} "
          f"{'error %':>10} {'time s':>8}")
    for r in results:
        print(f"{r['order']:>6d} {r['ne']:>8d} {r['ndof']:>10d} "
              f"{r['Hz_origin']:>14.6e} {r['error_pct']:>10.3f} "
              f"{r['elapsed_sec']:>8.1f}")

    with open(args.out_json, "w") as fh:
        json.dump({"mu_r": args.mu_r, "H0": args.H0,
                   "results": results}, fh, indent=2)
    print(f"\nSaved: {args.out_json}")
