"""Kameari-CLN v22e: axisymmetric Cu disk WITHOUT Kelvin (simple outer Dirichlet).

Diagnostic: if v22b's τ_dot vs τ_B2 mismatch (15-43% from n=1) is caused by
Kelvin pullback, then this version (no Kelvin) should give <1% match.
If still mismatched, Kelvin is innocent and J leakage / iteration is suspect.

Setup identical to v22b except:
- No outer Kelvin half-circle, no GND vertex
- Single half-circle of radius A_outer = 5*R_disk = 50mm
- Outer arc has Dirichlet u=0 (truncation BC)
"""
import sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace", line_buffering=True)
sys.path.insert(0, r"S:/Radia/01_GitHub/src/radia")
print("Starting v22e (no Kelvin, outer Dirichlet)...", flush=True)

from netgen.occ import (WorkPlane, Pnt, Axes, X, Z, MoveTo,
                         Glue, OCCGeometry)
from ngsolve import (Mesh, H1, BilinearForm, LinearForm,
                     GridFunction, CoefficientFunction, grad, dx,
                     x, y, IfPos, sqrt, Integrate, TaskManager, ngsglobals)
from math import pi
import json, time
from pathlib import Path

# === Parameters (matching v22b for direct comparison) ===
R_disk = 10e-3
half_thick = 1e-3
sigma_Cu = 5.8e7
mu0 = 4 * pi * 1e-7
nu0 = 1.0 / mu0

A_outer = 50e-3  # outer truncation radius (= a_kelvin in v22b)
maxh = 2e-3
maxh_disk = 0.1e-3
ORDER = 3
N_STAGES = 6
B0 = 1.0

TAU_INF_CYL_US = (mu0 * sigma_Cu * R_disk**2) / 3.83171**2 * 1e6


def build_geo():
    """Single half-circle: disk + air, with outer Dirichlet truncation."""
    # Outer half-circle (replaces the Kelvin two-sphere setup)
    wp_outer = WorkPlane()
    outer_full = wp_outer.Circle(A_outer).Face()
    outer_full.name = "air"

    # Disk
    wp_disk = WorkPlane(Axes((0, -half_thick, 0), n=Z, h=X))
    disk_full = wp_disk.Rectangle(R_disk, 2 * half_thick).Face()
    disk_full.name = "conductor"
    disk_full.maxh = maxh_disk

    # Cut to half-circle (r >= 0)
    cutter = MoveTo(-A_outer - 0.1, -A_outer - 0.1).Rectangle(
        A_outer + 0.1, 2 * A_outer + 0.2).Face()
    outer_half = outer_full - cutter
    disk_half = disk_full - cutter

    # Air region
    air_only = outer_half - disk_half
    air_only.name = "air"

    # Edge naming
    for edge in air_only.edges:
        cx = edge.center.x
        try:
            v0, v1 = edge.vertices
            d0 = sqrt(v0.p.x ** 2 + v0.p.y ** 2)
            d1 = sqrt(v1.p.x ** 2 + v1.p.y ** 2)
            is_outer_arc = (abs(d0 - A_outer) < 0.01 * A_outer
                             and abs(d1 - A_outer) < 0.01 * A_outer
                             and cx > 0.001 * A_outer)
        except Exception:
            is_outer_arc = False

        if cx < 0.001 * A_outer:
            edge.name = "axis"
        elif is_outer_arc:
            edge.name = "outer_arc"
        else:
            edge.name = "interface"

    for edge in disk_half.edges:
        cx = edge.center.x
        if cx < 0.001 * A_outer:
            edge.name = "axis"
        else:
            edge.name = "disk_bnd"

    shape = Glue([air_only, disk_half])
    return OCCGeometry(shape, dim=2)


def main():
    ngsglobals.msg_level = 0
    print(f"=== v22e: axisymmetric Cu disk WITHOUT Kelvin ===", flush=True)
    print(f"  Disk: R={R_disk*1000:.1f}mm, t={2*half_thick*1000:.1f}mm",
          flush=True)
    print(f"  Outer truncation radius: {A_outer*1000:.0f}mm",
          flush=True)
    print(f"  N_STAGES = {N_STAGES}\n", flush=True)

    geo = build_geo()
    print("Generating mesh...", flush=True)
    t0 = time.time()
    with TaskManager():
        ngmesh = geo.GenerateMesh(maxh=maxh, grading=0.5)
    mesh = Mesh(ngmesh)
    print(f"  ne = {mesh.ne}  ({time.time()-t0:.1f}s)", flush=True)
    print(f"  materials = {mesh.GetMaterials()}", flush=True)
    print(f"  boundaries = {mesh.GetBoundaries()}\n", flush=True)

    r_weight = IfPos(x - 1e-12, x, 1e-12)
    nu_cf = CoefficientFunction(nu0)  # uniform vacuum (no Kelvin)
    sigma_cf = mesh.MaterialCF({"conductor": sigma_Cu}, default=0.0)
    sigma_inv_cf = mesh.MaterialCF({"conductor": 1.0 / sigma_Cu}, default=0.0)

    # FE spaces — no Periodic (no Kelvin); Dirichlet on axis + outer_arc
    fes_u = H1(mesh, order=ORDER, dirichlet="axis|outer_arc")
    fes_J = H1(mesh, order=ORDER, definedon=mesh.Materials("conductor"),
                dirichlet="axis")
    print(f"  fes_u (H1, no Periodic) ndof = {fes_u.ndof}", flush=True)
    print(f"  fes_J (H1 on conductor) ndof = {fes_J.ndof}", flush=True)

    u, v = fes_u.TnT()

    a_form = BilinearForm(fes_u)
    a_form += nu_cf / r_weight * grad(u) * grad(v) * dx
    print("  Assembling+factor...", flush=True)
    t0 = time.time()
    with TaskManager():
        a_form.Assemble()
        try:
            inv_u = a_form.mat.Inverse(fes_u.FreeDofs(), inverse="pardiso")
        except Exception:
            inv_u = a_form.mat.Inverse(fes_u.FreeDofs(),
                                        inverse="sparsecholesky")
    print(f"    {time.time()-t0:.1f}s\n", flush=True)

    A_phi_imposed_cf = B0 * x / 2
    gfJ = GridFunction(fes_J)
    gfJ.Set(sigma_Cu * A_phi_imposed_cf,
             definedon=mesh.Materials("conductor"))

    proj_err = float(Integrate(
        (gfJ - sigma_Cu * A_phi_imposed_cf) ** 2 * dx("conductor"), mesh))
    proj_norm = float(Integrate(
        (sigma_Cu * A_phi_imposed_cf) ** 2 * dx("conductor"), mesh))
    print(f"  ||gfJ - σA||/||σA|| = {(proj_err/proj_norm)**0.5:.4e}",
          flush=True)

    twopi = 2 * pi
    R_inv_check = twopi * float(Integrate(
        gfJ * gfJ * sigma_inv_cf * r_weight * dx("conductor"), mesh))
    R0_test = 1.0 / R_inv_check
    R_inv_anal = pi * sigma_Cu * B0**2 * half_thick * R_disk**4 / 4
    print(f"  R_0: {R0_test:.6e} vs analytical {1/R_inv_anal:.6e}, "
          f"ratio {R0_test * R_inv_anal:.6f}\n", flush=True)

    diag = []
    J_history = []
    gfu_pot = GridFunction(fes_u)
    gfu_pot.vec[:] = 0.0
    sign_flip_first_n = None

    for n in range(N_STAGES):
        t_stage = time.time()

        f = LinearForm(fes_u)
        f += gfJ * v * dx("conductor")
        with TaskManager():
            f.Assemble()
        gfu_n = GridFunction(fes_u)
        with TaskManager():
            gfu_n.vec.data = inv_u * f.vec

        R_inv = twopi * float(Integrate(
            gfJ * gfJ * sigma_inv_cf * r_weight * dx("conductor"), mesh))
        if R_inv < 1e-30:
            print(f"  [{n}] R_inv collapsed", flush=True)
            break
        R_n = 1.0 / R_inv

        cross = []
        for gfJ_m in J_history:
            xnm = twopi * float(Integrate(
                gfJ * gfJ_m * sigma_inv_cf * r_weight * dx("conductor"), mesh))
            cross.append(xnm)
        max_cross = max((abs(c) for c in cross), default=0.0)
        rel_drift = max_cross / abs(R_inv) if R_inv else 0.0

        gfu_pot.vec.data += R_n * gfu_n.vec

        Ln_int = twopi * float(Integrate(
            gfJ * gfu_pot * dx("conductor"), mesh))
        L_n = R_n * Ln_int
        tau_dot_us = Ln_int * 1e6

        B2_int = twopi * float(Integrate(
            nu_cf / r_weight * grad(gfu_n) * grad(gfu_n) * dx, mesh))
        L_n_B2 = R_n * R_n * B2_int
        tau_B2_us = L_n_B2 / R_n * 1e6 if R_n > 0 else None

        signL = "+" if L_n > 0 else "-"
        elapsed = time.time() - t_stage
        td = f"{tau_dot_us:.4f}" if tau_dot_us > 0 else "N/A"
        tb = f"{tau_B2_us:.4f}" if tau_B2_us and tau_B2_us > 0 else "N/A"
        ratio = (tau_dot_us / tau_B2_us) if (tau_B2_us and tau_B2_us > 0) else 0
        print(f"  [{n}] drift={rel_drift:.3e}  En={R_inv:.3e}  "
              f"R={R_n:.3e}  L={L_n:.3e}({signL})  "
              f"τ_dot={td}μs  τ_B2={tb}μs  ratio={ratio:.4f}  "
              f"({elapsed:.1f}s)", flush=True)
        diag.append({"n": n, "R_n": R_n, "L_n": L_n, "L_n_B2": L_n_B2,
                     "tau_dot_us": tau_dot_us, "tau_B2_us": tau_B2_us,
                     "ratio_dot_over_B2": ratio,
                     "energy_norm": R_inv, "rel_drift": rel_drift,
                     "sign_L": signL})

        if L_n <= 0 and sign_flip_first_n is None:
            sign_flip_first_n = n

        gfJ_save = GridFunction(fes_J)
        gfJ_save.vec.data = gfJ.vec.data
        J_history.append(gfJ_save)

        if abs(L_n) > 1e-30:
            J_new_cf = gfJ - sigma_cf * (gfu_pot / r_weight) * (1.0 / L_n)
            gfJ_new = GridFunction(fes_J)
            gfJ_new.Set(J_new_cf, definedon=mesh.Materials("conductor"))
            gfJ = gfJ_new
        else:
            break

    print("\n" + "=" * 78, flush=True)
    print("v22e (no Kelvin) RESULTS — diagnostic for Kelvin innocence", flush=True)
    print("=" * 78, flush=True)
    print(f"{'n':>3} {'drift':>12} {'τ_dot μs':>12} {'τ_B2 μs':>12} "
          f"{'ratio':>10}", flush=True)
    for d in diag:
        td = f"{d['tau_dot_us']:.4f}" if d['tau_dot_us'] > 0 else "N/A"
        tb = f"{d['tau_B2_us']:.4f}" if d['tau_B2_us'] and d['tau_B2_us'] > 0 else "N/A"
        print(f"{d['n']:>3} {d['rel_drift']:>12.3e} {td:>12} {tb:>12} "
              f"{d['ratio_dot_over_B2']:>10.4f}", flush=True)

    if diag:
        worst = max(abs(d['ratio_dot_over_B2'] - 1.0)
                    for d in diag if d['tau_B2_us'] and d['tau_B2_us'] > 0)
        print(f"\nWorst |τ_dot/τ_B2 - 1| = {worst:.4f} ({worst*100:.2f}%)",
              flush=True)
        if worst < 0.01:
            print("→ <1% match: KELVIN was the culprit!", flush=True)
        else:
            print("→ Still mismatched: Kelvin is INNOCENT, look elsewhere "
                  "(J leakage, FE basis, etc.)", flush=True)

    out = {"method": "v22e: axisymmetric disk, NO Kelvin (outer Dirichlet)",
           "purpose": "diagnostic: isolate Kelvin contribution to τ mismatch",
           "params": {"R_disk_mm": R_disk*1000,
                       "thickness_mm": 2*half_thick*1000,
                       "A_outer_mm": A_outer*1000,
                       "order": ORDER, "ne": mesh.ne,
                       "fes_u_ndof": fes_u.ndof, "fes_J_ndof": fes_J.ndof},
           "stages": diag,
           "sign_flip_first_n": sign_flip_first_n}
    out_path = Path(__file__).parent / "cuboid_521_kameari_kelvin_v22e_no_kelvin.json"
    out_path.write_text(json.dumps(out, indent=2))
    print(f"\nResults: {out_path}", flush=True)


if __name__ == "__main__":
    main()
