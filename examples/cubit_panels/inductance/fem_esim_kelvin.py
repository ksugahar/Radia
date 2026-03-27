"""
FEM prototype: A-formulation + Kelvin + ESIM Robin BC (2D axisymmetric).

Cross-check against BEM + ESIM (impedance_esim.py).

Two modes:
  A) FEM-full:  workpiece volume meshed with sigma (reference)
  B) FEM-ESIM:  workpiece replaced by ESIM Robin BC on its surface

Both use Kelvin transformation for open boundary (no PML).

2D axisymmetric A-formulation:
  Variable: u = r * A_theta
  Equation: -div(nu/r * grad(u)) + j*omega*sigma*u/r = J_theta  (in workpiece)
  Kelvin: nu_K = nu0 * (rho'/a)^2 in exterior domain

Usage:
    python fem_esim_kelvin.py                  # Both modes, compare
    python fem_esim_kelvin.py --mode full      # FEM-full only
    python fem_esim_kelvin.py --mode esim      # FEM-ESIM only
"""

import math
import os
import sys
import time
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../src/radia'))

from ngsolve import *
from netgen.occ import *

MU_0 = 4e-7 * np.pi
NU_0 = 1.0 / MU_0


# ============================================================
# Geometry: air + coil + workpiece + Kelvin exterior
# ============================================================
def build_geometry(R_coil, a_coil, R_wp, H_wp, a_kelvin, z_offset,
                   include_workpiece_volume=True, maxh=0.003):
    """Build 2D axisymmetric geometry.

    Args:
        R_coil: Coil center radius [m]
        a_coil: Coil wire radius [m] (modeled as square cross-section)
        R_wp: Workpiece cylinder radius [m]
        H_wp: Workpiece cylinder height [m]
        a_kelvin: Kelvin boundary radius [m]
        z_offset: Z-offset for exterior domain [m]
        include_workpiece_volume: If True, mesh workpiece interior (FEM-full)
        maxh: Mesh size [m]

    Returns:
        NGSolve Mesh
    """
    # Interior domain: half-circle
    wp_inner = WorkPlane()
    inner_full = wp_inner.Circle(a_kelvin).Face()

    # Cut to half (x >= 0)
    cutter = MoveTo(-a_kelvin - 0.1, -a_kelvin - 0.1).Rectangle(
        a_kelvin + 0.1, 2 * a_kelvin + 0.2).Face()
    inner_half = inner_full - cutter

    # Coil cross-section
    w = a_coil * 2
    coil = MoveTo(R_coil - w / 2, -w / 2).Rectangle(w, w).Face()
    coil.name = "coil"
    coil.maxh = maxh / 2

    # Workpiece cross-section (rectangle in r-z plane)
    wp_rect = MoveTo(0, -H_wp / 2).Rectangle(R_wp, H_wp).Face()
    # Name workpiece edges BEFORE boolean (names survive subtraction)
    for edge in wp_rect.edges:
        cx, cy = edge.center.x, edge.center.y
        if cx < 0.001:
            edge.name = "axis"  # left edge on axis
        else:
            edge.name = "workpiece_bnd"

    if include_workpiece_volume:
        # FEM-full: workpiece is a meshed region
        wp_rect.name = "workpiece"
        wp_rect.maxh = maxh / 3  # finer mesh for skin effect

        # Air = inner_half - coil - workpiece
        air = inner_half - coil - wp_rect
        air.name = "air_inner"
        shapes = [air, coil, wp_rect]
    else:
        # FEM-ESIM: workpiece is a hole (not meshed)
        air = inner_half - coil - wp_rect
        air.name = "air_inner"
        shapes = [air, coil]

    # Exterior domain (Kelvin)
    wp_ext = WorkPlane(Axes((0, z_offset, 0), n=Z, h=X))
    outer_full = wp_ext.Circle(a_kelvin).Face()
    cutter_ext = MoveTo(-a_kelvin - 0.1, z_offset - a_kelvin - 0.1).Rectangle(
        a_kelvin + 0.1, 2 * a_kelvin + 0.2).Face()
    outer_half = outer_full - cutter_ext
    outer_half.name = "air_outer"

    # GND point
    gnd = Vertex(Pnt(0, z_offset, 0))
    gnd.name = "GND"

    # Name edges - detect workpiece boundary by vertex positions
    for edge in air.edges:
        cx = edge.center.x
        try:
            v0, v1 = edge.vertices
            d0 = math.sqrt(v0.p.x**2 + v0.p.y**2)
            d1 = math.sqrt(v1.p.x**2 + v1.p.y**2)
            is_arc = abs(d0 - a_kelvin) < 0.01 and abs(d1 - a_kelvin) < 0.01 and cx > 0.01

            # Workpiece boundary: edge on the workpiece rectangle with at least
            # one vertex at x > 0 (exclude pure axis edges)
            on_wp_rect = False
            if not include_workpiece_volume and max(v0.p.x, v1.p.x) > 1e-6:
                # Check if vertices lie on workpiece rectangle boundary
                def on_rect(px, py):
                    on_right = abs(px - R_wp) < 1e-4 and -H_wp/2 - 1e-4 <= py <= H_wp/2 + 1e-4
                    on_top = abs(py - H_wp/2) < 1e-4 and -1e-4 <= px <= R_wp + 1e-4
                    on_bot = abs(py + H_wp/2) < 1e-4 and -1e-4 <= px <= R_wp + 1e-4
                    return on_right or on_top or on_bot
                on_wp_rect = on_rect(v0.p.x, v0.p.y) and on_rect(v1.p.x, v1.p.y)
        except:
            is_arc = False
            on_wp_rect = False

        if on_wp_rect:
            edge.name = "workpiece_bnd"
        elif is_arc:
            edge.name = "kelvin_int"
        elif cx < 1e-4:
            edge.name = "axis"
        else:
            edge.name = "default"

    for edge in coil.edges:
        edge.name = "coil_bnd" if edge.center.x > 0.01 else "axis"

    if include_workpiece_volume:
        for edge in wp_rect.edges:
            edge.name = "wp_bnd" if edge.center.x > 0.001 else "axis"

    # Kelvin exterior edges
    kelvin_int_edges = []
    kelvin_ext_edges = []

    for edge in air.edges:
        if hasattr(edge, 'name') and edge.name == "kelvin_int":
            kelvin_int_edges.append(edge)

    for edge in outer_half.edges:
        cx = edge.center.x
        cy = edge.center.y - z_offset
        try:
            v0, v1 = edge.vertices
            d0 = math.sqrt(v0.p.x**2 + (v0.p.y - z_offset)**2)
            d1 = math.sqrt(v1.p.x**2 + (v1.p.y - z_offset)**2)
            is_arc = abs(d0 - a_kelvin) < 0.01 and abs(d1 - a_kelvin) < 0.01 and cx > 0.01
        except:
            is_arc = False

        if cx < 0.01:
            edge.name = "axis_ext"
        elif is_arc:
            edge.name = "kelvin_ext"
            kelvin_ext_edges.append(edge)
        else:
            edge.name = "default"

    # Periodic identification
    for int_e in kelvin_int_edges:
        for ext_e in kelvin_ext_edges:
            iy = int_e.center.y
            ey = ext_e.center.y - z_offset
            if (iy > 0 and ey > 0) or (iy < 0 and ey < 0):
                int_e.Identify(ext_e, "kelvin", IdentificationType.PERIODIC)
                break

    shape = Glue(shapes + [outer_half, gnd])
    geo = OCCGeometry(shape, dim=2)
    mesh = Mesh(geo.GenerateMesh(maxh=maxh, grading=0.4))
    return mesh


# ============================================================
# FEM-full: workpiece with sigma (direct eddy current solve)
# ============================================================
def solve_fem_full(R_coil, a_coil, R_wp, H_wp, sigma, frequency,
                   I_total=1.0, a_kelvin=0.15, z_offset=0.4,
                   maxh=0.002, order=3):
    """Solve frequency-domain eddy current with meshed workpiece.

    Returns dict with P_total, Q_total, L, mesh info.
    """
    omega = 2 * np.pi * frequency
    J0 = I_total / (2 * a_coil)**2  # Current density

    print("  Building geometry (workpiece meshed)...")
    t0 = time.perf_counter()
    mesh = build_geometry(R_coil, a_coil, R_wp, H_wp, a_kelvin, z_offset,
                          include_workpiece_volume=True, maxh=maxh)
    t_mesh = time.perf_counter() - t0

    materials = mesh.GetMaterials()
    boundaries = mesh.GetBoundaries()
    print(f"  Mesh: {mesh.ne} elements, materials={materials}")
    print(f"  Boundaries: {boundaries}")

    # FE space
    fes_base = H1(mesh, order=order, complex=True,
                  dirichlet="axis|axis_ext", dirichlet_bbnd="GND")
    fes = Periodic(fes_base)
    u, v = fes.TnT()
    print(f"  DOFs: {fes.ndof}")

    # Coefficient functions
    r_cf = IfPos(x - 1e-10, x, 1e-10)

    # Kelvin factor for exterior
    y_local = y - z_offset
    rho_prime = sqrt(x**2 + y_local**2)
    rho_safe = IfPos(rho_prime - 1e-10, rho_prime, 1e-10)
    kelvin_fac = (rho_safe / a_kelvin)**2

    # Reluctivity per material
    nu_dict = {}
    for mat in materials:
        if "outer" in mat.lower():
            nu_dict[mat] = NU_0 * kelvin_fac
        else:
            nu_dict[mat] = NU_0
    nu_cf = mesh.MaterialCF(nu_dict, default=NU_0)

    # Conductivity (only in workpiece)
    sigma_dict = {}
    for mat in materials:
        if "workpiece" in mat.lower():
            sigma_dict[mat] = sigma
        else:
            sigma_dict[mat] = 0.0
    sigma_cf = mesh.MaterialCF(sigma_dict, default=0.0)

    # Bilinear form
    a_bf = BilinearForm(fes)
    a_bf += nu_cf / r_cf * grad(u) * grad(v) * dx
    # Eddy current term: j*omega*sigma * (u/r)*(v/r) * r = j*omega*sigma * u*v/r
    a_bf += 1j * omega * sigma_cf * u * v / r_cf * dx

    # Source
    f_lf = LinearForm(fes)
    f_lf += J0 * v * dx("coil")

    print("  Assembling...")
    t0 = time.perf_counter()
    a_bf.Assemble()
    f_lf.Assemble()
    t_asm = time.perf_counter() - t0

    print("  Solving...")
    t0 = time.perf_counter()
    gfu = GridFunction(fes)
    gfu.vec.data = a_bf.mat.Inverse(fes.FreeDofs(), inverse="pardiso") * f_lf.vec
    t_solve = time.perf_counter() - t0

    # Post-process: power loss in workpiece
    # P = 0.5 * sigma * omega^2 * |A_theta|^2 * volume
    # In 2D axi: P = 2*pi * integral( 0.5 * sigma * omega^2 * |u/r|^2 * r dr dz )
    #           = 2*pi * integral( 0.5 * sigma * omega^2 * |u|^2 / r dr dz )
    P_density = 0.5 * sigma_cf * omega**2 * Conj(gfu) * gfu / r_cf
    P_per_radian = Integrate(P_density, mesh, definedon=mesh.Materials("workpiece")).real
    P_total = 2 * np.pi * P_per_radian

    # Reactive power (stored magnetic energy rate)
    # Q = 2*pi * integral( 0.5 * omega * nu * |grad(u)|^2 / r )
    # Only meaningful for the workpiece region
    Q_density = 0.5 * omega * sigma_cf * Conj(gfu) * gfu / r_cf
    Q_per_radian = Integrate(Q_density, mesh, definedon=mesh.Materials("workpiece")).real
    Q_total = 2 * np.pi * Q_per_radian

    # Inductance from stored energy
    # W_mag = 0.5 * L * I^2 = 2*pi * integral( 0.5 * nu * |grad(u)|^2 / r )
    W_per_radian = Integrate(
        0.5 * nu_cf / r_cf * (grad(gfu) * Conj(grad(gfu))),
        mesh).real
    W_mag = 2 * np.pi * W_per_radian
    L = 2 * W_mag / (I_total**2)

    return {
        'mode': 'FEM-full',
        'P_total': float(P_total),
        'Q_total': float(Q_total),
        'L': float(L),
        'ndof': fes.ndof,
        'ne': mesh.ne,
        't_mesh': t_mesh,
        't_asm': t_asm,
        't_solve': t_solve,
    }


# ============================================================
# FEM-ESIM: workpiece as Robin BC
# ============================================================
def solve_fem_sibc(R_coil, a_coil, R_wp, H_wp, sigma, frequency,
                   bh_curve=None, mu_r=1.0,
                   I_total=1.0, a_kelvin=0.15, z_offset=0.4,
                   maxh=0.003, order=3, max_iter=15, tol=1e-4):
    """Solve with SIBC Robin BC + Karl Hollaus iteration.

    Karl iteration (Hollaus, IEEE Trans. Mag., 2025):
      1. Solve FEM with current Z_s as Robin BC on workpiece surface
      2. Compute H_t from solution (gradient sampling near boundary)
      3. Update Z_s from ESIM cell problem at average H_t
      4. Repeat until Z_s converges

    Robin BC weak form (2D axi, u = r*A_theta):
      int (nu/r) grad(u).grad(v) dx - (jw/Zs) int u*v/r ds = int J*v dx

    Power (Poynting flux, correct factor = pi not 2*pi):
      P = pi * w^2 * Re(Zs)/|Zs|^2 * int |u|^2/r ds

    Returns dict with P_total, L, iteration info.
    """
    from esim_cell_problem import ESIMFiniteSlabSolver

    omega = 2 * np.pi * frequency
    J0 = I_total / (2 * a_coil)**2

    print("  Building geometry (workpiece = hole)...")
    t0 = time.perf_counter()
    mesh = build_geometry(R_coil, a_coil, R_wp, H_wp, a_kelvin, z_offset,
                          include_workpiece_volume=False, maxh=maxh)
    t_mesh = time.perf_counter() - t0

    materials = mesh.GetMaterials()
    print(f"  Mesh: {mesh.ne} elements")

    r_cf = IfPos(x - 1e-10, x, 1e-10)

    # Kelvin
    y_local = y - z_offset
    rho_prime = sqrt(x**2 + y_local**2)
    rho_safe = IfPos(rho_prime - 1e-10, rho_prime, 1e-10)
    kelvin_fac = (rho_safe / a_kelvin)**2
    nu_dict = {m: NU_0 * kelvin_fac if "outer" in m.lower() else NU_0
               for m in materials}
    nu_cf = mesh.MaterialCF(nu_dict, default=NU_0)

    # ESIM solver
    esim_solver = ESIMFiniteSlabSolver(
        half_thickness=R_wp, bh_curve=bh_curve, sigma=sigma,
        frequency=frequency,
        mu_r=mu_r if bh_curve is None else None, n_nodes=200)

    # Initial Z_s
    sol0 = esim_solver.solve(5.0)
    Z_s = sol0['Z']
    print(f"  Initial Z_s = {Z_s:.4e}, delta = {sol0['delta']*1e3:.3f} mm")

    # Sample points for H_t evaluation (just outside workpiece boundary)
    eps = maxh * 0.2
    pts_lat = [(R_wp + eps, z) for z in
               np.linspace(-H_wp / 2 + 0.001, H_wp / 2 - 0.001, 20)]
    pts_cap = [(r, H_wp / 2 + eps) for r in np.linspace(0.002, R_wp - 0.001, 8)]
    pts_cap += [(r, -H_wp / 2 - eps) for r in np.linspace(0.002, R_wp - 0.001, 8)]
    sample_pts = pts_lat + pts_cap

    # --- Karl iteration ---
    print(f"  {'Iter':>4s} {'|Z_s|':>12s} {'P [W]':>12s} {'H_t_avg':>10s} {'dZ/Z':>10s}")

    t0_solve = time.perf_counter()
    for iteration in range(max_iter):
        # Robin BC: -(jw/Zs) * int u*v/r ds
        robin = -1j * omega / Z_s

        fesc = Periodic(H1(mesh, order=order, complex=True,
                           dirichlet="axis|axis_ext", dirichlet_bbnd="GND"))
        uc, vc = fesc.TnT()
        bf = BilinearForm(fesc)
        bf += nu_cf / r_cf * grad(uc) * grad(vc) * dx
        bf += robin * uc * vc / r_cf * ds("workpiece_bnd")
        lfc = LinearForm(fesc)
        lfc += J0 * vc * dx("coil")
        bf.Assemble()
        lfc.Assemble()
        gfu = GridFunction(fesc)
        gfu.vec.data = bf.mat.Inverse(fesc.FreeDofs(), inverse="pardiso") * lfc.vec

        # Power: P = pi * w^2 * Re(Zs)/|Zs|^2 * int|u|^2/r ds
        int_u2r = Integrate(gfu * Conj(gfu) / r_cf, mesh,
                            definedon=mesh.Boundaries("workpiece_bnd")).real
        P = np.pi * omega**2 * Z_s.real / abs(Z_s)**2 * int_u2r

        # Sample H_t
        grad_u = grad(gfu)
        H_t_vals = []
        for (r_pt, z_pt) in sample_pts:
            try:
                mip = mesh(r_pt, z_pt)
                H_z = abs(NU_0 * grad_u[0](mip) / max(r_pt, 1e-10))
                H_r = abs(NU_0 * grad_u[1](mip) / max(r_pt, 1e-10))
                H_t_vals.append(max(H_z, H_r))
            except:
                H_t_vals.append(1e-3)
        H_t_avg = np.mean(H_t_vals)

        # Update Z_s (under-relaxation)
        Z_s_old = Z_s
        sol_new = esim_solver.solve(max(float(H_t_avg), 1e-3), relaxation=0.3)
        Z_s = 0.5 * Z_s_old + 0.5 * sol_new['Z']

        dZ = abs(Z_s - Z_s_old) / abs(Z_s_old)
        print(f"  {iteration:4d} {abs(Z_s):12.4e} {P:12.4e} {H_t_avg:10.2f} {dZ:10.4e}")

        if dZ < tol and iteration > 0:
            break

    t_solve = time.perf_counter() - t0_solve

    # Inductance
    W = 2 * np.pi * Integrate(
        0.5 * nu_cf / r_cf * grad(gfu) * Conj(grad(gfu)), mesh).real
    L = 2 * W / I_total**2

    return {
        'mode': 'FEM-SIBC (Karl)',
        'P_total': float(P),
        'L': float(L),
        'Z_s': Z_s,
        'H_t_avg': float(H_t_avg),
        'iterations': iteration + 1,
        'ndof': fesc.ndof,
        'ne': mesh.ne,
        't_mesh': t_mesh,
        't_solve': t_solve,
    }


# ============================================================
# Main
# ============================================================
def main(mode='both'):
    # Copper at 1kHz: delta = 2.1 mm, resolvable with maxh=0.5mm
    R_coil = 0.030      # Coil radius [m]
    a_coil = 0.003      # Coil wire radius [m]
    R_wp = 0.010         # Workpiece radius [m]
    H_wp = 0.020         # Workpiece height [m]
    sigma = 5.8e7        # Copper [S/m]
    frequency = 1000     # 1 kHz -> delta ~ 2.1 mm (resolvable)
    mu_r = 1.0
    I_total = 1.0

    rho = 1.0 / sigma
    delta = math.sqrt(2 * rho / (2 * np.pi * frequency * MU_0 * mu_r))

    print("=" * 65)
    print("FEM-full vs FEM-ESIM (2D axisymmetric + Kelvin)")
    print("=" * 65)
    print(f"Coil:      R={R_coil*1e3:.0f} mm, a={a_coil*1e3:.1f} mm")
    print(f"Workpiece: R={R_wp*1e3:.0f} mm, H={H_wp*1e3:.0f} mm")
    print(f"Material:  Copper (sigma={sigma:.1e} S/m, mu_r={mu_r})")
    print(f"Frequency: {frequency} Hz")
    print(f"Skin depth: {delta*1e3:.2f} mm (xi = R_wp/delta = {R_wp/delta:.2f})")
    print()

    results = {}

    if mode in ('both', 'full'):
        print("[FEM-full] Workpiece volume meshed (direct eddy current)")
        print("-" * 65)
        # maxh_wp = delta/3 for adequate skin layer resolution
        r_full = solve_fem_full(
            R_coil, a_coil, R_wp, H_wp, sigma, frequency,
            I_total=I_total, maxh=0.0005, order=3,
            a_kelvin=0.15, z_offset=0.4)
        results['full'] = r_full
        print(f"  P = {r_full['P_total']:.6e} W")
        print(f"  L = {r_full['L']*1e9:.2f} nH")
        print(f"  Time: mesh={r_full['t_mesh']:.1f}s, "
              f"asm={r_full['t_asm']:.1f}s, solve={r_full['t_solve']:.1f}s")
        print()

    if mode in ('both', 'esim'):
        print("[FEM-SIBC] External field -> ESIM on workpiece surface")
        print("-" * 65)
        r_esim = solve_fem_sibc(
            R_coil, a_coil, R_wp, H_wp, sigma, frequency,
            mu_r=mu_r,
            I_total=I_total, maxh=0.003, order=3,
            a_kelvin=0.15, z_offset=0.4)
        results['esim'] = r_esim
        print(f"  P = {r_esim['P_total']:.6e} W")
        print(f"  L = {r_esim['L']*1e9:.2f} nH")
        print(f"  Time: mesh={r_esim['t_mesh']:.1f}s, "
              f"solve={r_esim['t_solve']:.1f}s")
        print()

    # Comparison
    if 'full' in results and 'esim' in results:
        print("=" * 65)
        print("Comparison: FEM-full vs FEM-ESIM")
        print("=" * 65)
        P_f = results['full']['P_total']
        P_e = results['esim']['P_total']
        L_f = results['full']['L']
        L_e = results['esim']['L']
        print(f"  {'':>20s}  {'FEM-full':>12s}  {'FEM-SIBC':>12s}  {'diff':>8s}")
        print(f"  {'P [W]':>20s}  {P_f:12.4e}  {P_e:12.4e}  "
              f"{(P_e-P_f)/max(abs(P_f),1e-30)*100:+8.2f}%")
        print(f"  {'L [nH]':>20s}  {L_f*1e9:12.4f}  {L_e*1e9:12.4f}  "
              f"{(L_e-L_f)/max(abs(L_f),1e-30)*100:+8.2f}%")
        print(f"  {'DOFs':>20s}  {results['full']['ndof']:>12d}  "
              f"{results['esim']['ndof']:>12d}")
        print(f"  {'Elements':>20s}  {results['full']['ne']:>12d}  "
              f"{results['esim']['ne']:>12d}")
        print()
        print(f"  FEM-full: resolves skin layer (maxh << delta)")
        print(f"  FEM-ESIM: static solve + ESIM 1D cell problem")
        print(f"  delta = {delta*1e3:.2f} mm, xi = {R_wp/delta:.2f}")

    return results


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description="FEM A-formulation + Kelvin + ESIM prototype")
    parser.add_argument("--mode", default="both",
                        choices=["both", "full", "esim"])
    args = parser.parse_args()
    main(mode=args.mode)
