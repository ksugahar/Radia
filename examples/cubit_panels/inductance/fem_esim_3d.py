"""
3D FEM-SIBC: gapped torus coil + cylindrical workpiece with ESIM surface impedance.

SIBC = Robin BC on conductor surface.  Conductor interior NOT solved.
Validated 2026-04-14 (2D axisym Kelvin): hole + Robin gives L<1%, P<2%.
For 3D HCurl, Robin = (jw/Z_s) * A_t . v_t on hole boundary (ds("sibc")).

Formulation (thin conducting shell on internal interface):
  int nu0 * curl(A) . curl(v) dx + (jw/Z_s) * int A_t . v_t ds(wp) = int J . v dx

  dx covers ALL domains (air + workpiece_interior + coil)
  ds("wp_surface") is the internal interface between air and workpiece

Physics:
  [n x H] = J_s = -(jw/Z_s) * A_t  (jump of tangential H = surface current)
  E_t = Z_s * J_s  (Ohm's law in the thin shell)

Karl iteration for nonlinear Z_s(H_t):
  1. Solve FEM with SIBC penalty
  2. H_t = |n x H| = |jw/Z_s| * |A_t| on wp_surface
  3. Update Z_s from ESIM cell problem
  4. Repeat until Z_s converges

Usage:
    python fem_esim_3d.py
    python fem_esim_3d.py --freq 7000 --material steel
    python fem_esim_3d.py --freq 1000 --material copper
"""

import math
import os
import sys
import time
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../src/radia'))

MU_0 = 4e-7 * np.pi


def run(R_coil=0.030, a_coil=0.003, gap_deg=5,
        R_wp=0.010, H_wp=0.020,
        sigma=2e6, frequency=7000, material="steel",
        R_air=0.06, maxh_air=0.012, maxh_coil=0.005,
        maxh_wp_bnd=None,
        R_kelvin=0.12,
        order=1, esim_geometry='cylinder',
        approach='hole'):
    """3D FEM-SIBC with Kelvin transform for open boundary.

    Args:
        R_coil, a_coil: Torus major/minor radius [m]
        gap_deg: Gap angle [degrees] (same as BEM)
        R_wp, H_wp: Workpiece cylinder radius/height [m]
        sigma: Workpiece conductivity [S/m]
        frequency: [Hz]
        material: steel/copper/aluminum
        R_air: Air sphere / Kelvin boundary radius [m]
        maxh_air: Air mesh size [m]
        maxh_coil: Coil mesh size [m]
        maxh_wp_bnd: WP surface mesh size [m]. If None, uses min(R_wp/4, maxh_coil).
        R_kelvin: Kelvin exterior domain outer radius [m] (0 = no Kelvin, Dirichlet on R_air)
        order: HCurl polynomial order
        esim_geometry: 'cylinder' (Bessel I0/I1) or 'slab' (cosh/sinh)
    """
    from ngsolve import (Mesh, HCurl, BilinearForm, LinearForm, GridFunction,
                         Integrate, curl, dx, ds, CF, BND, VOL, TaskManager,
                         InnerProduct)
    from netgen.occ import (WorkPlane, Axes, Axis, Pnt, Dir, Sphere,
                             Cylinder, OCCGeometry, Glue)
    from esim_cell_problem import ESIMFiniteSlabSolver

    STEEL_BH = [
        [0, 0], [50, 0.1], [100, 0.25], [200, 0.55],
        [500, 0.95], [1000, 1.2], [2000, 1.4], [5000, 1.55],
        [10000, 1.65], [20000, 1.75], [50000, 1.9], [100000, 2.0],
    ]
    bh_curve = STEEL_BH if material == "steel" else None
    mu_r = 1.0 if material in ("copper", "aluminum") else None

    I_total = 1.0
    J0 = I_total / (math.pi * a_coil**2)

    omega = 2 * np.pi * frequency
    nu0 = 1.0 / MU_0

    # Skin depth estimate (linear mu)
    mu_eff = MU_0 * (mu_r if mu_r else 100.0)  # rough estimate for steel
    delta = math.sqrt(2.0 / (omega * mu_eff * sigma))

    print("=" * 65)
    print("3D FEM-SIBC: Gapped Torus + Workpiece (internal interface)")
    print("=" * 65)
    print(f"Coil:      R={R_coil*1e3:.0f}mm, a={a_coil*1e3:.1f}mm, gap={gap_deg}deg")
    print(f"Workpiece: R={R_wp*1e3:.0f}mm, H={H_wp*1e3:.0f}mm")
    print(f"Material:  {material}, sigma={sigma:.0e}, f={frequency}Hz")
    print(f"Skin depth: delta={delta*1e3:.2f}mm, delta/R={delta/R_wp:.2f}")
    print(f"Air sphere: R={R_air*1e3:.0f}mm")
    print()

    # === Geometry ===
    print(f"[1/3] Building 3D geometry (approach={approach})...")
    t0 = time.perf_counter()

    # Gapped torus (coil volume)
    wp_plane = WorkPlane(Axes(p=Pnt(R_coil, 0, 0),
                              n=Dir(0, 1, 0), h=Dir(0, 0, 1)))
    circle = wp_plane.Circle(a_coil).Face()
    torus = circle.Revolve(Axis(Pnt(0, 0, 0), Dir(0, 0, 1)), 360 - gap_deg)
    torus.name = "coil"
    torus.maxh = maxh_coil

    # Air sphere (inner physical domain)
    air_sphere = Sphere(Pnt(0, 0, 0), R_air)
    air_sphere.name = "air"

    # Kelvin exterior domain (exterior mapped domain)
    use_kelvin = R_kelvin > R_air
    if use_kelvin:
        kelvin_outer = Sphere(Pnt(0, 0, 0), R_kelvin)
        kelvin_ext = kelvin_outer - Sphere(Pnt(0, 0, 0), R_air)
        kelvin_ext.name = "kelvin"
        kelvin_ext.maxh = maxh_air * 2
        for f in kelvin_outer.faces:
            f.name = "outer"
    else:
        for f in air_sphere.faces:
            f.name = "outer"

    # Workpiece cylinder
    wp_cyl = Cylinder(Pnt(0, 0, -H_wp / 2), Dir(0, 0, 1), R_wp, H_wp)
    wp_bnd_h = maxh_wp_bnd if maxh_wp_bnd is not None else min(R_wp / 4, maxh_coil)
    for f in wp_cyl.faces:
        f.name = "wp_surface"
        f.maxh = wp_bnd_h

    if approach == 'hole':
        air = air_sphere - torus - wp_cyl
        air.name = "air"
        parts = [air, torus]
    else:
        wp_cyl.name = "workpiece"
        air = air_sphere - torus - wp_cyl
        air.name = "air"
        parts = [air, wp_cyl, torus]

    if use_kelvin:
        parts.append(kelvin_ext)

    shape = Glue(parts)

    geo = OCCGeometry(shape)
    with TaskManager():
        ngmesh = geo.GenerateMesh(maxh=maxh_air, grading=0.3)
    mesh = Mesh(ngmesh)
    # Curve the mesh geometry. Practical default Curve(>=3) for 3D IH:
    # the torus+sphere+cylinder combination needs order 3 to match
    # 2D axisym reference to within 1% P (measured 2026-04-14, Cu 7kHz).
    # Curve(2) under-predicts L by ~10%; Curve(3) brings it to ~-5%.
    # Raising FES order to 2 instead (at Curve(2)) is far more expensive
    # (60x solve time) and does NOT improve P. Prefer Curve(3) over
    # HCurl order=2 for this geometry class. See
    # `src/radia/mcp_server/ih/ih_knowledge.py` -> "Curve Order vs FES Order"
    # for the measured tradeoff table.
    mesh.Curve(max(order + 1, 3))
    t_mesh = time.perf_counter() - t0

    ne = mesh.ne
    mats = mesh.GetMaterials()
    bnds = mesh.GetBoundaries()
    print(f"  Mesh: {ne} tets, {mesh.nv} vertices ({t_mesh:.1f}s)")
    print(f"  Materials: {mats}")
    print(f"  Boundaries: {bnds}")

    if approach == 'interface' and "workpiece" not in mats:
        print("  ERROR: 'workpiece' material not found!")
        print("  Workpiece was probably not meshed as a separate volume.")
        return None

    # Check wp_surface boundary
    has_wp = "wp_surface" in bnds
    if not has_wp:
        print("  WARNING: 'wp_surface' boundary not found.")
        print("  The interface may have a different name.")
        # List all boundaries and their types
        for i, b in enumerate(bnds):
            if b and b != "default":
                print(f"    [{i}] {b}")

    # === FEM-SIBC solve with Karl iteration ===
    print("[2/3] FEM-SIBC solve...")
    t0 = time.perf_counter()

    fes = HCurl(mesh, order=order, dirichlet="outer", complex=True)
    u, v = fes.TnT()
    print(f"  DOFs: {fes.ndof}")

    # Source: J = J0 * e_phi in coil
    from ngsolve import x, y, z, sqrt, IfPos
    r_xy = sqrt(x * x + y * y)
    r_safe = IfPos(r_xy - 1e-10, r_xy, 1e-10)
    J_source = J0 * CF((-y / r_safe, x / r_safe, 0))

    # RHS (assembled once, constant)
    f_lf = LinearForm(fes)
    f_lf += J_source * v * dx("coil")
    f_lf.Assemble()

    # Pre-compute wp_surface area
    if has_wp:
        wp_region = mesh.Boundaries("wp_surface")
        A_wp = Integrate(CF(1), mesh, BND, definedon=wp_region).real
        print(f"  A_wp = {A_wp:.6e} m^2")
    else:
        A_wp = 2 * math.pi * R_wp * H_wp  # analytical cylinder lateral area
        print(f"  A_wp (analytical) = {A_wp:.6e} m^2")

    # ESIM solver
    esim_solver = ESIMFiniteSlabSolver(
        half_thickness=R_wp, bh_curve=bh_curve, sigma=sigma,
        frequency=frequency,
        mu_r=mu_r if bh_curve is None else None, n_nodes=200,
        geometry=esim_geometry)
    print(f"  ESIM geometry: {esim_geometry}")

    # Karl iteration
    Z_s = esim_solver.solve(5.0)['Z']
    gfu = GridFunction(fes)
    max_iter = 15
    tol = 1e-3
    relax = 0.5

    print(f"  |Z_s| init = {abs(Z_s):.4e}")
    print(f"  Z_s/(jw*mu0*a) = {abs(Z_s)/(omega*MU_0*R_wp):.2f}")
    print(f"  {'Iter':>4s} {'|Z_s|':>12s} {'P [W]':>12s} "
          f"{'H_t_rms':>10s} {'dZ/Z':>10s}")

    P_total = 0.0
    Q_total = 0.0
    H_t_rms = 5.0

    # Determine SIBC boundary name
    sibc_bnd = "wp_surface" if has_wp else None

    # Reluctivity: nu0 everywhere, nu0*(r/R_air)^4 in Kelvin exterior domain
    has_kelvin = "kelvin" in mats
    if has_kelvin:
        r_sq = x * x + y * y + z * z
        kelvin_fac = R_air**2 / (r_sq + 1e-20)  # (a/r)^2
        nu_dict = {}
        for m in mats:
            if m == "kelvin":
                nu_dict[m] = nu0 * kelvin_fac
            elif m == "workpiece" and approach == 'interface':
                mu_r_eff = 1.0
                if bh_curve:
                    for i in range(len(bh_curve) - 1):
                        if bh_curve[i][0] <= 500 <= bh_curve[i + 1][0]:
                            B_mid = (bh_curve[i][1] +
                                     (bh_curve[i + 1][1] - bh_curve[i][1]) *
                                     (500 - bh_curve[i][0]) /
                                     (bh_curve[i + 1][0] - bh_curve[i][0]))
                            mu_r_eff = B_mid / (MU_0 * 500)
                            break
                    if mu_r_eff < 1:
                        mu_r_eff = 1.0
                elif mu_r is not None:
                    mu_r_eff = mu_r
                nu_dict[m] = nu0 / mu_r_eff
                print(f"  Workpiece mu_r_eff = {mu_r_eff:.0f}")
            else:
                nu_dict[m] = nu0
        nu_cf = mesh.MaterialCF(nu_dict, default=nu0)
        print(f"  Kelvin transform: R_air={R_air*1e3:.0f}mm, "
              f"R_kelvin={R_kelvin*1e3:.0f}mm")
    elif approach == 'interface' and "workpiece" in mats:
        mu_r_eff = 1.0
        if bh_curve:
            for i in range(len(bh_curve) - 1):
                if bh_curve[i][0] <= 500 <= bh_curve[i + 1][0]:
                    B_mid = (bh_curve[i][1] +
                             (bh_curve[i + 1][1] - bh_curve[i][1]) *
                             (500 - bh_curve[i][0]) /
                             (bh_curve[i + 1][0] - bh_curve[i][0]))
                    mu_r_eff = B_mid / (MU_0 * 500)
                    break
            if mu_r_eff < 1:
                mu_r_eff = 1.0
        elif mu_r is not None:
            mu_r_eff = mu_r
        nu_cf = mesh.MaterialCF({"workpiece": nu0 / mu_r_eff}, default=nu0)
        print(f"  Workpiece mu_r_eff = {mu_r_eff:.0f}")
    else:
        nu_cf = nu0

    for iteration in range(max_iter):
        # Assemble: curl-curl (with mu_r in workpiece) + SIBC penalty
        a_bf = BilinearForm(fes)
        a_bf += nu_cf * curl(u) * curl(v) * dx
        if sibc_bnd:
            sibc_coeff = 1j * omega / Z_s
            a_bf += sibc_coeff * u.Trace() * v.Trace() * ds(sibc_bnd)
        a_bf.Assemble()

        with TaskManager():
            gfu.vec.data = a_bf.mat.Inverse(
                fes.FreeDofs(), inverse="pardiso") * f_lf.vec

        # Compute H_t on workpiece surface using SIBC relation:
        #   J_s = -(jw/Z_s) * A_t  =>  H_t = |J_s| = |jw/Z_s| * |A_t|
        #
        # NOTE: Do NOT use curl(A) for H_t. On an internal interface,
        # curl(A)/mu0 gives the FULL tangential H (incident + scattered),
        # dominated by the incident coil field (~18 A/m). The SIBC-based
        # H_t gives the physical surface current (~0.7 A/m), which is the
        # correct input for the ESIM cell problem.
        if sibc_bnd:
            # |A_t|^2 = |A|^2 - |A.n|^2 (project out normal component).
            # HCurl basis has NON-zero normal trace on a surface; summing
            # gfu[i]^2 over i=0..2 gives |A|^2 (full vector), not the
            # tangential magnitude. The physical SIBC relation
            # J_s = -(jw/Z_s) A_t needs the TANGENTIAL A only.
            from ngsolve import specialcf
            nrm = specialcf.normal(3)
            A_dot_n_re = sum(gfu[i].real * nrm[i] for i in range(3))
            A_dot_n_im = sum(gfu[i].imag * nrm[i] for i in range(3))
            A_n_sq = A_dot_n_re * A_dot_n_re + A_dot_n_im * A_dot_n_im
            A_full_sq = sum(gfu[i].real * gfu[i].real +
                            gfu[i].imag * gfu[i].imag for i in range(3))
            At_sq = A_full_sq - A_n_sq
            At_int = Integrate(At_sq, mesh, BND,
                               definedon=wp_region).real
            At_rms = math.sqrt(max(At_int, 0) / max(A_wp, 1e-30))
            H_t_rms = abs(1j * omega / Z_s) * At_rms

            if iteration == 0:
                print(f"  DEBUG: |A_t| rms = {At_rms:.4e}, "
                      f"|jw/Z_s| = {abs(1j*omega/Z_s):.4e}")
                print(f"  DEBUG: H_t (SIBC) = {H_t_rms:.2f} A/m")
                # Also show curl-based H for reference (includes incident)
                H_cf = nu0 * curl(gfu)
                H_mag_sq = sum(H_cf[i].real * H_cf[i].real +
                               H_cf[i].imag * H_cf[i].imag
                               for i in range(3))
                Hmag_int = Integrate(H_mag_sq, mesh, BND,
                                     definedon=wp_region).real
                H_curl = math.sqrt(max(Hmag_int, 0) / max(A_wp, 1e-30))
                print(f"  DEBUG: H (curl, includes incident) = {H_curl:.2f}")
        else:
            H_t_rms = 5.0

        # ESIM update
        Z_s_old = Z_s
        sol = esim_solver.solve(max(H_t_rms, 1e-3))
        Z_s_new = sol['Z']
        Z_s = relax * Z_s_new + (1 - relax) * Z_s_old

        P_total = sol['P_prime'] * A_wp
        Q_total = sol['Q_prime'] * A_wp

        dZ = abs(Z_s - Z_s_old) / max(abs(Z_s_old), 1e-30)
        print(f"  {iteration:4d} {abs(Z_s):12.4e} {P_total:12.4e} "
              f"{H_t_rms:10.2f} {dZ:10.4e}")

        if dZ < tol and iteration > 0:
            print(f"  Converged at iteration {iteration}")
            break

    t_solve = time.perf_counter() - t0

    # Inductance: volume magnetic energy + skin-layer reactive energy.
    # SIBC excludes the conductor interior from the mesh, so the skin
    # layer's stored magnetic energy is missing from the volume integral.
    # Add surface term from Im(Z_s) (same formula as 2D axisym, see
    # fem_esim_kelvin.solve_fem_sibc). Without this, L is under-predicted
    # by ~5% at R/delta ~ 30.
    curl_A = curl(gfu)
    B_sq = sum(curl_A[i].real * curl_A[i].real +
               curl_A[i].imag * curl_A[i].imag for i in range(3))
    W = 0.5 * Integrate(nu0 * B_sq, mesh).real
    L_vol = 2 * W / I_total**2
    if sibc_bnd:
        # L_skin = omega * Im(Z_s) / |Z_s|^2 * int |A_t|^2 dA  (3D surface)
        # At_int was computed above for H_t_rms (phasor |A_t|^2 integral).
        L_skin = (omega * Z_s.imag / abs(Z_s)**2 * At_int
                  / I_total**2)
    else:
        L_skin = 0.0
    L = L_vol + L_skin

    # === Results ===
    print()
    print("[3/3] Results")
    print("-" * 65)
    print(f"  P (FEM-SIBC 3D)   = {P_total:.6e} W")
    print(f"  Q                  = {Q_total:.6e} var")
    print(f"  H_t_rms            = {H_t_rms:.2f} A/m")
    print(f"  |Z_s|              = {abs(Z_s):.4e} Ohm")
    print(f"  L (volume only)    = {L_vol*1e9:.3f} nH")
    print(f"  L (+skin surface)  = {L*1e9:.3f} nH  (skin adds {L_skin*1e9:+.3f} nH)")
    print(f"  DOFs: {fes.ndof}, Elements: {ne}")
    print(f"  Time: mesh={t_mesh:.1f}s, solve={t_solve:.1f}s")
    print("-" * 65)

    return {
        'P_total': float(P_total),
        'Q_total': float(Q_total),
        'L': float(L),
        'H_t_rms': float(H_t_rms),
        'Z_s': complex(Z_s),
        'ndof': fes.ndof,
        'ne': ne,
        't_mesh': t_mesh,
        't_solve': t_solve,
    }


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="3D FEM-SIBC")
    parser.add_argument("--freq", type=float, default=7000)
    parser.add_argument("--material", default="steel")
    parser.add_argument("--maxh", type=float, default=0.012)
    parser.add_argument("--order", type=int, default=1,
                        help="HCurl polynomial order (1 or 2)")
    parser.add_argument("--geometry", default="local_curvature",
                        choices=["local_curvature", "none"],
                        help="ESIM curvature: local_curvature (Bessel) or none (flat slab)")
    parser.add_argument("--approach", default="hole",
                        choices=["hole", "interface"],
                        help="hole: cavity (PMC baseline, Karl's FEM-ESIM), "
                             "interface: workpiece as air (transparent baseline)")
    args = parser.parse_args()

    # Material-dependent sigma
    sigma_map = {'steel': 2e6, 'copper': 5.8e7, 'aluminum': 3.5e7}
    sigma = sigma_map.get(args.material, 2e6)
    _geom_map = {"local_curvature": "cylinder", "none": "slab"}
    esim_geom = _geom_map.get(args.geometry, "cylinder")
    r = run(frequency=args.freq, material=args.material, sigma=sigma,
            maxh_air=args.maxh, esim_geometry=esim_geom,
            approach=args.approach, order=args.order)

    if r is None:
        sys.exit(1)

    # Compare with EFIE-SIBC (BEM two-way)
    print()
    print("Running EFIE-SIBC (BEM) for comparison...")
    sys.path.insert(0, os.path.dirname(__file__))
    from pmchwt_sibc import run as efie_run
    r_bem = efie_run(material=args.material, frequency=args.freq, verbose=False)

    print()
    print("=" * 65)
    print(f"Comparison: FEM-SIBC vs EFIE-SIBC ({args.material}, {args.freq}Hz)")
    print("=" * 65)
    print(f"{'':>20s} {'FEM-SIBC':>14s} {'EFIE-SIBC':>14s} {'diff':>8s}")
    P_f, P_b = r['P_total'], r_bem['P_total']
    if abs(P_b) > 1e-30:
        print(f"{'P [W]':>20s} {P_f:14.4e} {P_b:14.4e} "
              f"{(P_f-P_b)/P_b*100:+8.1f}%")
    else:
        print(f"{'P [W]':>20s} {P_f:14.4e} {P_b:14.4e}")
    print(f"{'H_t_rms [A/m]':>20s} {r['H_t_rms']:14.2f} "
          f"{r_bem['H_t_rms']:14.2f}")
