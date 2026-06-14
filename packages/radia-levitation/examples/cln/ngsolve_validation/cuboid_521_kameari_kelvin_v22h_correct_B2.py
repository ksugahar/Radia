"""Kameari + Kelvin v22: Axisymmetric Cu disk Kameari-CLN proof of concept (2026-05-05).

Sugahara guidance: 3D HCurl Kameari-CLN (v15-v21) struggled to reach
canonical answers due to J·n leakage at conductor surfaces. In axisymmetric
the J = J_φ ê_φ has only azimuthal component, so J·n = 0 is satisfied
by geometry — no leakage possible. v22 demonstrates Kameari-CLN works
cleanly in axisymmetric.

Setup:
- Cu disk (R=10mm, full thickness 2mm) at origin
- Inner air half-circle (radius a_kelvin = 50mm) on r >= 0 side of (r,z) plane
- Outer Kelvin half-circle (same radius) at z = z_offset (Sugahara two-sphere)
- Periodic identification on Kelvin arcs
- GND vertex at exterior center (image of infinity)

Discretization:
- u = r × A_φ in H1 (Periodic, full domain, dirichlet="axis|axis_ext|GND")
- J_φ in H1 with definedon="conductor" — DISCONTINUOUS at conductor surface
  (J = 0 outside, naturally enforced — no current leak into air)

Weak form (axisymmetric A-formulation):
- a(u, v) = ∫_2D (ν/r) ∇u · ∇v dr dz
- f(v) = ∫_2D J_φ × v dr dz (in conductor)

Inner products (3D = 2π × 2D-with-r-weight):
- R_inv = 2π × ∫_2D (J_φ)²/σ × r dr dz (in conductor)
- L_n via dot: τ = 2π × ∫_2D J_φ × u_pot dr dz (no r — A_φ=u/r cancels with r)
- L_n via B²: 2π × ∫_2D (ν/r) × |∇u|² dr dz (per stage)

Schmidt update (J side):
- J_{n+1} = J_n - σ × (u_pot_accum / r) / L_n, projected onto fes_J
"""
import sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace", line_buffering=True)
sys.path.insert(0, r"S:/Radia/01_GitHub/src/radia")
print("Starting v22h (B² cross form + material breakdown) (refined mesh + Order=3)...", flush=True)

from netgen.occ import (WorkPlane, Vertex, Pnt, Axes, X, Y, Z, MoveTo,
                         Glue, OCCGeometry, IdentificationType, Dir)
from ngsolve import (Mesh, H1, Periodic, BilinearForm, LinearForm,
                     GridFunction, CoefficientFunction, grad, dx,
                     x, y, IfPos, sqrt, Integrate, TaskManager, ngsglobals)
from kelvin_source import kelvin_nu_factor_axisym_cf, build_material_cf
from math import pi
import json, time
from pathlib import Path
import numpy as np

# === Parameters ===
R_disk = 10e-3       # disk radius [m]
half_thick = 1e-3    # disk half-thickness [m] (full 2*half_thick = 2mm)
sigma_Cu = 5.8e7
mu0 = 4 * pi * 1e-7
nu0 = 1.0 / mu0

a_kelvin = 50e-3     # Kelvin boundary radius
z_offset = 5 * a_kelvin
maxh = 2e-3       # 5mm -> 2mm (refined)
maxh_disk = 0.1e-3  # 0.5mm -> 0.1mm (5x finer in disk)
ORDER = 3         # 2 -> 3 (Kameari/Sugahara recommendation)
N_STAGES = 6
B0 = 1.0

# Analytical reference (infinite cylinder eddy mode):
# τ_n = μ_0 σ R² / j_{1,n}²,  j_{1,1} ≈ 3.83171
TAU_INF_CYL_US = (mu0 * sigma_Cu * R_disk**2) / 3.83171**2 * 1e6


def build_geo():
    """Two-half-circle axisymmetric geometry per sphere_with_Kelvin pattern."""
    # Inner half-circle
    wp1 = WorkPlane()
    inner_full = wp1.Circle(a_kelvin).Face()
    inner_full.name = "air"

    # Disk (rectangle in r-z plane, r in [0, R_disk], z in [-half_thick, half_thick])
    wp_disk = WorkPlane(Axes((0, -half_thick, 0), n=Z, h=X))
    disk_full = wp_disk.Rectangle(R_disk, 2 * half_thick).Face()
    disk_full.name = "conductor"
    disk_full.maxh = maxh_disk

    # Outer Kelvin half-circle
    wp3 = WorkPlane(Axes((0, z_offset, 0), n=Z, h=X))
    outer_full = wp3.Circle(a_kelvin).Face()
    outer_full.name = "kelvin"

    # GND vertex
    gnd_point = Vertex(Pnt(0, z_offset, 0))
    gnd_point.name = "GND"

    # Cut to half (keep r >= 0)
    cutter = MoveTo(-a_kelvin - 0.1, -a_kelvin - 0.1).Rectangle(
        a_kelvin + 0.1, 2 * a_kelvin + 0.2).Face()
    cutter_ext = MoveTo(-a_kelvin - 0.1, z_offset - a_kelvin - 0.1).Rectangle(
        a_kelvin + 0.1, 2 * a_kelvin + 0.2).Face()
    inner_half = inner_full - cutter
    disk_half = disk_full - cutter
    outer_half = outer_full - cutter_ext

    # Inner air = inner_half - disk_half
    inner_air = inner_half - disk_half
    inner_air.name = "air"

    # Edge naming (axis, kelvin_int, kelvin_ext, GND adjacent edges)
    kelvin_inner_edges = []
    for edge in inner_air.edges:
        cx = edge.center.x
        try:
            v0, v1 = edge.vertices
            d0 = sqrt(v0.p.x ** 2 + v0.p.y ** 2)
            d1 = sqrt(v1.p.x ** 2 + v1.p.y ** 2)
            is_kelvin_arc = (abs(d0 - a_kelvin) < 0.01 * a_kelvin
                              and abs(d1 - a_kelvin) < 0.01 * a_kelvin
                              and cx > 0.001 * a_kelvin)
        except Exception:
            is_kelvin_arc = False

        if cx < 0.001 * a_kelvin:
            edge.name = "axis"
        elif is_kelvin_arc:
            edge.name = "kelvin_int"
            kelvin_inner_edges.append(edge)
        else:
            edge.name = "interface"

    for edge in disk_half.edges:
        cx = edge.center.x
        if cx < 0.001 * a_kelvin:
            edge.name = "axis"
        else:
            edge.name = "disk_bnd"

    kelvin_outer_edges = []
    for edge in outer_half.edges:
        cx = edge.center.x
        try:
            v0, v1 = edge.vertices
            d0 = sqrt(v0.p.x ** 2 + (v0.p.y - z_offset) ** 2)
            d1 = sqrt(v1.p.x ** 2 + (v1.p.y - z_offset) ** 2)
            is_kelvin_arc = (abs(d0 - a_kelvin) < 0.01 * a_kelvin
                              and abs(d1 - a_kelvin) < 0.01 * a_kelvin
                              and cx > 0.001 * a_kelvin)
        except Exception:
            is_kelvin_arc = False

        if cx < 0.001 * a_kelvin:
            edge.name = "axis_ext"
        elif is_kelvin_arc:
            edge.name = "kelvin_ext"
            kelvin_outer_edges.append(edge)
        else:
            edge.name = "default"

    # Periodic identification of Kelvin arcs (match by y sign)
    matched = 0
    for ie in kelvin_inner_edges:
        iy = ie.center.y
        for oe in kelvin_outer_edges:
            oy = oe.center.y - z_offset
            if (iy > 0 and oy > 0) or (iy < 0 and oy < 0):
                ie.Identify(oe, "kelvin", IdentificationType.PERIODIC)
                matched += 1
                break
    print(f"  Periodic edge pairs matched: {matched}", flush=True)

    shape = Glue([inner_air, disk_half, outer_half, gnd_point])
    return OCCGeometry(shape, dim=2)


def main():
    ngsglobals.msg_level = 0
    print(f"=== v22: axisymmetric Cu disk Kameari-CLN ===", flush=True)
    print(f"  Disk: R={R_disk*1000:.1f}mm, thickness={2*half_thick*1000:.1f}mm",
          flush=True)
    print(f"  Kelvin: a={a_kelvin*1000:.0f}mm, z_offset={z_offset*1000:.0f}mm",
          flush=True)
    print(f"  Analytical (∞-cylinder leading) τ_1 = {TAU_INF_CYL_US:.2f} μs",
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

    # === Material with Kelvin pullback ===
    r_weight = IfPos(x - 1e-12, x, 1e-12)
    nu_kelvin_factor = kelvin_nu_factor_axisym_cf(z_offset=z_offset, R=a_kelvin)
    nu_cf = build_material_cf(mesh, nu0, nu_kelvin_factor)
    sigma_cf = mesh.MaterialCF({"conductor": sigma_Cu}, default=0.0)
    sigma_inv_cf = mesh.MaterialCF({"conductor": 1.0 / sigma_Cu}, default=0.0)

    # === FE spaces ===
    # u = r * A_phi (continuous, full domain)
    fes_u = Periodic(H1(mesh, order=ORDER, dirichlet="axis|axis_ext|GND"))
    # J_phi (discontinuous at conductor surface: P2 on conductor, 0 outside)
    fes_J = H1(mesh, order=ORDER, definedon=mesh.Materials("conductor"),
                dirichlet="axis")
    print(f"  fes_u (Periodic H1) ndof = {fes_u.ndof}", flush=True)
    print(f"  fes_J (H1 on conductor) ndof = {fes_J.ndof}", flush=True)

    u, v = fes_u.TnT()

    # === Bilinear form (axisymmetric A-formulation) ===
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

    # === Initial J_φ = σ × A_φ_imposed = σ × B_0 × r/2 (in conductor) ===
    A_phi_imposed_cf = B0 * x / 2  # = B_0 × r / 2
    gfJ = GridFunction(fes_J, name="J_imp")
    gfJ.Set(sigma_Cu * A_phi_imposed_cf,
             definedon=mesh.Materials("conductor"))

    # Diagnostic: how well does the H1 P2 projection match analytical?
    proj_err_int = float(Integrate(
        (gfJ - sigma_Cu * A_phi_imposed_cf) ** 2 * dx("conductor"), mesh))
    proj_norm_int = float(Integrate(
        (sigma_Cu * A_phi_imposed_cf) ** 2 * dx("conductor"), mesh))
    rel_proj_err = (proj_err_int / proj_norm_int) ** 0.5 if proj_norm_int else 0
    print(f"  ||gfJ - σA_phys||/||σA_phys|| in conductor = "
          f"{rel_proj_err:.4e}", flush=True)

    # R_0 verification (Sugahara CLN convention with 3D = 2π × 2D-with-r):
    twopi = 2 * pi
    R_inv_check = twopi * float(Integrate(
        gfJ * gfJ * sigma_inv_cf * r_weight * dx("conductor"), mesh))
    R0_test = 1.0 / R_inv_check if R_inv_check > 1e-30 else float('inf')
    # Analytical R_0 for cylinder (sigma A_phys with A_phys = B_0 r/2):
    # R_inv_anal = 2π × σ × (B_0/2)² × ∫ r³ × 2 half_thick dr
    #            = 2π × σ × B_0²/4 × (R^4/4) × (2 half_thick)
    #            = π × σ × B_0² × half_thick × R^4 / 4
    R_inv_anal = pi * sigma_Cu * B0**2 * half_thick * R_disk**4 / 4
    R0_anal = 1.0 / R_inv_anal
    print(f"  R_0 verification: {R0_test:.6e} vs analytical {R0_anal:.6e}, "
          f"ratio {R0_test/R0_anal:.6f}\n", flush=True)

    # === Kameari iteration ===
    diag = []
    J_history = []
    gfu_pot = GridFunction(fes_u, name="u_pot_accum")
    gfu_pot.vec[:] = 0.0
    sign_flip_first_n = None

    for n in range(N_STAGES):
        t_stage = time.time()

        # Linear form: f(v) = ∫_2D J_φ × v dr dz (axisymmetric reduction
        # of ∫_3D J_φ × v_φ dV after v_φ = v/r reduction)
        f = LinearForm(fes_u)
        f += gfJ * v * dx("conductor")
        with TaskManager():
            f.Assemble()
        gfu_n = GridFunction(fes_u)
        with TaskManager():
            gfu_n.vec.data = inv_u * f.vec

        # R_n
        R_inv = twopi * float(Integrate(
            gfJ * gfJ * sigma_inv_cf * r_weight * dx("conductor"), mesh))
        if R_inv < 1e-30:
            print(f"  [{n}] R_inv collapsed ({R_inv:.2e})", flush=True)
            break
        R_n = 1.0 / R_inv

        # Schmidt drift: <J_n, J_m/sigma>_C for m<n  (3D inner product)
        cross = []
        for gfJ_m in J_history:
            xnm = twopi * float(Integrate(
                gfJ * gfJ_m * sigma_inv_cf * r_weight * dx("conductor"), mesh))
            cross.append(xnm)
        max_cross = max((abs(c) for c in cross), default=0.0)
        rel_drift = max_cross / abs(R_inv) if R_inv else 0.0

        gfu_pot.vec.data += R_n * gfu_n.vec  # accumulator

        # τ_n via dot form: 3D inner <J, A>_C with A_φ = u/r
        # = 2π × ∫ J × (u/r) × r dr dz = 2π × ∫ J × u dr dz (no r weight)
        Ln_int = twopi * float(Integrate(
            gfJ * gfu_pot * dx("conductor"), mesh))
        L_n = R_n * Ln_int
        tau_dot_us = Ln_int * 1e6  # = L_n / R_n × 1e6

        # === v22h diagnostic: B² in CROSS form (FE identity test) ===
        # By FE identity:
        #   ∫_full (ν/r) ∇u_n · ∇u_pot dr dz × 2π
        #   = ∫_C J_n × u_pot dr dz × 2π = Ln_int  (machine precision!)
        # Material breakdown to verify Kelvin contribution.
        B2_cross_full = twopi * float(Integrate(
            nu_cf / r_weight * grad(gfu_n) * grad(gfu_pot) * dx, mesh))
        B2_cross_C = twopi * float(Integrate(
            nu_cf / r_weight * grad(gfu_n) * grad(gfu_pot) * dx("conductor"), mesh))
        B2_cross_air = twopi * float(Integrate(
            nu_cf / r_weight * grad(gfu_n) * grad(gfu_pot) * dx("air"), mesh))
        B2_cross_kelvin = twopi * float(Integrate(
            nu_cf / r_weight * grad(gfu_n) * grad(gfu_pot) * dx("kelvin"), mesh))
        # FE identity check
        FE_id_check = abs(B2_cross_full - Ln_int) / max(abs(Ln_int), 1e-30)

        # original per-stage diagonal (legacy)
        B2_int = twopi * float(Integrate(
            nu_cf / r_weight * grad(gfu_n) * grad(gfu_n) * dx, mesh))
        L_n_B2 = R_n * R_n * B2_int
        tau_B2_us = L_n_B2 / R_n * 1e6 if R_n > 0 else None
        tau_B2_cross_us = R_n * B2_cross_full * 1e6  # should = τ_dot

        signL = "+" if L_n > 0 else ("-" if L_n < 0 else "0")
        elapsed = time.time() - t_stage
        tau_dot_str = f"{tau_dot_us:.4f}" if tau_dot_us is not None and tau_dot_us > 0 else "N/A"
        tau_B2_str = f"{tau_B2_us:.4f}" if tau_B2_us is not None and tau_B2_us > 0 else "N/A"
        print(f"  [{n}] drift={rel_drift:.3e}  R={R_n:.3e}  "
              f"τ_dot={tau_dot_str}μs  τ_B2_diag={tau_B2_str}μs  "
              f"τ_B2_cross={tau_B2_cross_us:.4f}μs  FE_id_err={FE_id_check:.3e}",
              flush=True)
        print(f"      B²_cross breakdown: C={B2_cross_C:.4e}, "
              f"air={B2_cross_air:.4e}, kelvin={B2_cross_kelvin:.4e}, "
              f"sum={B2_cross_C+B2_cross_air+B2_cross_kelvin:.4e}",
              flush=True)
        diag.append({"n": n, "R_n": R_n, "L_n": L_n, "L_n_B2": L_n_B2,
                     "B2_int": B2_int, "B2_cross_full": B2_cross_full,
                     "B2_cross_C": B2_cross_C,
                     "B2_cross_air": B2_cross_air,
                     "B2_cross_kelvin": B2_cross_kelvin,
                     "tau_dot_us": tau_dot_us, "tau_B2_us": tau_B2_us,
                     "tau_B2_cross_us": tau_B2_cross_us,
                     "FE_id_err": FE_id_check,
                     "energy_norm": R_inv, "rel_drift": rel_drift,
                     "sign_L": signL})

        if L_n <= 0 and sign_flip_first_n is None:
            sign_flip_first_n = n
            print(f"      ** FIRST SIGN FLIP at n={n} **", flush=True)

        # Save J_n
        gfJ_save = GridFunction(fes_J)
        gfJ_save.vec.data = gfJ.vec.data
        J_history.append(gfJ_save)

        # Schmidt update on J side: J_{n+1} = J_n - σ × A_φ_pot_accum / L_n
        # A_φ_pot_accum = gfu_pot / r
        if abs(L_n) > 1e-30:
            J_new_cf = gfJ - sigma_cf * (gfu_pot / r_weight) * (1.0 / L_n)
            gfJ_new = GridFunction(fes_J)
            gfJ_new.Set(J_new_cf, definedon=mesh.Materials("conductor"))
            gfJ = gfJ_new
        else:
            print(f"  [{n}] |L_n| collapsed, stopping", flush=True)
            break

    # === Summary ===
    print("\n" + "=" * 78, flush=True)
    print("AXISYMMETRIC Cu DISK KAMEARI-CLN RESULTS", flush=True)
    print("=" * 78, flush=True)
    print(f"{'n':>3} {'rel_drift':>12} {'En':>12} {'R_n':>12} "
          f"{'τ_dot μs':>12} {'τ_B2 μs':>12}", flush=True)
    for d in diag:
        td = f"{d['tau_dot_us']:.4f}" if d['tau_dot_us'] and d['tau_dot_us'] > 0 else "N/A"
        tb = f"{d['tau_B2_us']:.4f}" if d['tau_B2_us'] and d['tau_B2_us'] > 0 else "N/A"
        print(f"{d['n']:>3} {d['rel_drift']:>12.3e} {d['energy_norm']:>12.3e} "
              f"{d['R_n']:>12.3e} {td:>12} {tb:>12}", flush=True)

    print(f"\n  Reference τ_1 (∞ cylinder) = {TAU_INF_CYL_US:.4f} μs",
          flush=True)
    if diag:
        print(f"  Stage 0 τ_dot / reference = "
              f"{diag[0]['tau_dot_us']/TAU_INF_CYL_US:.4f}", flush=True)

    out = {"method": "Axisymmetric Cu disk Kameari-CLN (v22)",
           "geometry": {"R_disk_mm": R_disk*1000,
                         "thickness_mm": 2*half_thick*1000,
                         "a_kelvin_mm": a_kelvin*1000,
                         "z_offset_mm": z_offset*1000},
           "params": {"sigma_S_per_m": sigma_Cu, "B0_T": B0,
                       "order": ORDER, "ne": mesh.ne,
                       "fes_u_ndof": fes_u.ndof, "fes_J_ndof": fes_J.ndof,
                       "rel_proj_err_J": rel_proj_err},
           "stages": diag,
           "R0_anal": R0_anal,
           "tau_inf_cylinder_us": TAU_INF_CYL_US,
           "sign_flip_first_n": sign_flip_first_n}
    out_path = Path(__file__).parent / "cuboid_521_kameari_kelvin_v22h_correct_B2.json"
    out_path.write_text(json.dumps(out, indent=2))
    print(f"\nResults: {out_path}", flush=True)


if __name__ == "__main__":
    main()
