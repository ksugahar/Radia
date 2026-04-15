"""
FEM scattered-field SIBC: gapped torus coil + cylindrical workpiece.

Scattered-field decomposition A = A_inc + A_scat:
  - A_inc: Biot-Savart from coil (interpolated into H1 GridFunctions)
  - H_inc: Biot-Savart H-field (interpolated, for n x H_inc boundary term)
  - A_scat: solved with Kelvin open boundary

RHS:
  f(v) = -(jw/Z_s) * <A_inc, v>_sibc  -  <n_cond x H_inc, v>_sibc

where n_cond = -specialcf.normal(3)  (outward from conductor).

Usage:
    python fem_scattered_coil.py
    python fem_scattered_coil.py --freq 7000 --material steel
"""

import math
import os
import sys
import time
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../src/radia'))

MU_0 = 4e-7 * np.pi


def generate_torus_segments(R_major, a_minor, gap_deg=5, n_per_turn=200):
    """Generate wire segments approximating a gapped torus coil.

    Returns:
        segments: list of ((x1,y1,z1), (x2,y2,z2)) in meters
    """
    gap_rad = math.radians(gap_deg)
    theta_max = 2 * math.pi - gap_rad
    segments = []
    for i in range(n_per_turn):
        t1 = theta_max * i / n_per_turn
        t2 = theta_max * (i + 1) / n_per_turn
        p1 = (R_major * math.cos(t1), R_major * math.sin(t1), 0)
        p2 = (R_major * math.cos(t2), R_major * math.sin(t2), 0)
        segments.append((p1, p2))
    return segments


def interpolate_biot_savart(mesh, segments, current):
    """Compute A_inc and H_inc from Biot-Savart, interpolate into H1 GFs.

    Returns:
        A_inc_cf: CoefficientFunction (dim=3) for vector potential
        H_inc_cf: CoefficientFunction (dim=3) for H-field
    """
    from ngsolve import H1, GridFunction, CF
    from biot_savart import h_segments_batch, a_segments_batch

    # Get all mesh vertex coordinates
    nv = mesh.nv
    pts = np.zeros((nv, 3))
    for i in range(nv):
        p = mesh.vertices[i].point
        pts[i] = [p[0], p[1], p[2]]

    # Compute fields at vertices
    seg_arr = np.array(segments)
    t0 = time.perf_counter()
    H_vals = h_segments_batch(seg_arr, pts, current=current)
    A_vals = a_segments_batch(seg_arr, pts, current=current)
    t_bs = time.perf_counter() - t0
    print(f"  Biot-Savart: {nv} points, {len(segments)} segments ({t_bs:.2f}s)")

    # Interpolate into H1 GridFunctions (order=1, DOFs at vertices)
    fes_h1 = H1(mesh, order=1)
    assert fes_h1.ndof >= nv, f"H1 ndof={fes_h1.ndof} < nv={nv}"

    gf_Hx = GridFunction(fes_h1)
    gf_Hy = GridFunction(fes_h1)
    gf_Hz = GridFunction(fes_h1)
    gf_Ax = GridFunction(fes_h1)
    gf_Ay = GridFunction(fes_h1)
    gf_Az = GridFunction(fes_h1)

    # H1 order=1: first nv DOFs are vertex DOFs
    gf_Hx.vec.FV().NumPy()[:nv] = H_vals[:, 0]
    gf_Hy.vec.FV().NumPy()[:nv] = H_vals[:, 1]
    gf_Hz.vec.FV().NumPy()[:nv] = H_vals[:, 2]
    gf_Ax.vec.FV().NumPy()[:nv] = A_vals[:, 0]
    gf_Ay.vec.FV().NumPy()[:nv] = A_vals[:, 1]
    gf_Az.vec.FV().NumPy()[:nv] = A_vals[:, 2]

    H_inc_cf = CF((gf_Hx, gf_Hy, gf_Hz))
    A_inc_cf = CF((gf_Ax, gf_Ay, gf_Az))

    # Report field magnitudes at origin and at workpiece surface
    H_center = h_segments_batch(seg_arr, np.array([[0, 0, 0]]), current)[0]
    print(f"  |H_inc| at center: {np.linalg.norm(H_center):.2f} A/m")

    return A_inc_cf, H_inc_cf


def run(R_coil=0.030, a_coil=0.003, gap_deg=5,
        R_wp=0.010, H_wp=0.020,
        sigma=2e6, frequency=7000, material="steel",
        R_air=0.10, maxh_air=0.012,
        R_kelvin=0.0, order=1,
        rhs_mode='boundary', maxh_sibc=None):
    """Scattered-field FEM-SIBC with Kelvin transform.

    Workpiece removed as a hole. SIBC on hole boundary.
    """
    from ngsolve import (Mesh, HCurl, H1, BilinearForm, LinearForm, GridFunction,
                         Integrate, curl, dx, ds, CF, BND, TaskManager,
                         specialcf, InnerProduct)
    from netgen.occ import (WorkPlane, Axes, Axis, Pnt, Dir, Sphere,
                             Cylinder, OCCGeometry, Glue, Vertex)
    from netgen.meshing import IdentificationType
    from esim_cell_problem import ESIMFiniteSlabSolver

    I_total = 1.0
    omega = 2 * np.pi * frequency
    nu0 = 1.0 / MU_0

    # Skin depth estimate
    mu_eff = MU_0 * (100.0 if material == "steel" else 1.0)
    delta = math.sqrt(2.0 / (omega * mu_eff * sigma))

    print("=" * 65)
    print("FEM Scattered-Field SIBC: Coil + Workpiece (hole approach)")
    print("=" * 65)
    print(f"Coil:      R={R_coil*1e3:.0f}mm, a={a_coil*1e3:.1f}mm, gap={gap_deg}deg")
    print(f"Workpiece: R={R_wp*1e3:.0f}mm, H={H_wp*1e3:.0f}mm")
    print(f"Material:  {material}, sigma={sigma:.0e}, f={frequency}Hz")
    print(f"Skin depth: delta={delta*1e3:.2f}mm")
    print()

    # === Generate coil segments ===
    print("[1/4] Coil segments...")
    segments = generate_torus_segments(R_coil, a_coil, gap_deg, n_per_turn=200)
    print(f"  {len(segments)} segments")

    # === Build geometry ===
    print("[2/4] Building geometry (hole approach)...")
    t0 = time.perf_counter()

    # Inner air sphere
    air_sphere = Sphere(Pnt(0, 0, 0), R_air)
    air_sphere.name = "air"

    # Workpiece cylinder (will be subtracted as hole)
    wp_cyl = Cylinder(Pnt(0, 0, -H_wp / 2), Dir(0, 0, 1), R_wp, H_wp)
    _maxh_sibc = maxh_sibc if maxh_sibc else min(R_wp / 4, maxh_air / 2)
    for f in wp_cyl.faces:
        f.name = "sibc"
        f.maxh = _maxh_sibc

    # Coil torus (for reference, but we don't mesh it - it's a filament)
    # Coil region would be meshed for J-source in total-field, but here
    # the coil is replaced by Biot-Savart. No coil volume needed.
    air = air_sphere - wp_cyl
    air.name = "air"

    # Kelvin domain
    use_kelvin = R_kelvin > R_air
    parts = [air]

    if use_kelvin:
        kelvin_outer = Sphere(Pnt(0, 0, 0), R_kelvin)
        kelvin_ext = kelvin_outer - Sphere(Pnt(0, 0, 0), R_air)
        kelvin_ext.name = "kelvin"
        kelvin_ext.maxh = maxh_air * 2

        for f_air in air.faces:
            if f_air.name != "sibc":
                f_air.name = "kelvin_int"
        for f_kel in kelvin_outer.faces:
            f_kel.name = "outer"

        for f_int in [f for f in air.faces if f.name == "kelvin_int"]:
            for f_ext in kelvin_ext.faces:
                if f_ext.name != "outer":
                    try:
                        f_int.Identify(f_ext, "kelvin",
                                       IdentificationType.CLOSESURFACES)
                        break
                    except Exception:
                        pass

        gnd = Vertex(Pnt(0, 0, 0))
        gnd.name = "GND"
        parts.extend([kelvin_ext, gnd])
        dirichlet_bc = "GND"
    else:
        for f in air_sphere.faces:
            f.name = "outer"
        dirichlet_bc = "outer"

    shape = Glue(parts)
    geo = OCCGeometry(shape)
    with TaskManager():
        ngmesh = geo.GenerateMesh(maxh=maxh_air, grading=0.3)
    mesh = Mesh(ngmesh)
    mesh.Curve(order + 1)
    t_mesh = time.perf_counter() - t0

    ne = mesh.ne
    mats = mesh.GetMaterials()
    bnds = mesh.GetBoundaries()
    print(f"  Mesh: {ne} tets, {mesh.nv} verts ({t_mesh:.1f}s)")
    print(f"  Materials: {mats}")
    print(f"  Boundaries: {bnds}")

    # Check sibc boundary
    has_sibc = "sibc" in bnds
    if not has_sibc:
        print("  ERROR: 'sibc' boundary not found!")
        return None

    # === Biot-Savart interpolation ===
    print("[3/4] Biot-Savart interpolation...")
    A_inc_cf, H_inc_cf = interpolate_biot_savart(mesh, segments, I_total)

    # === Scattered-field solve ===
    print("[4/4] Scattered-field solve...")

    # ESIM for initial Z_s
    STEEL_BH = [
        [0, 0], [50, 0.1], [100, 0.25], [200, 0.55],
        [500, 0.95], [1000, 1.2], [2000, 1.4], [5000, 1.55],
        [10000, 1.65], [20000, 1.75], [50000, 1.9], [100000, 2.0],
    ]
    bh_curve = STEEL_BH if material == "steel" else None
    mu_r = 1.0 if material in ("copper", "aluminum") else None
    esim_solver = ESIMFiniteSlabSolver(
        half_thickness=R_wp, bh_curve=bh_curve, sigma=sigma,
        frequency=frequency,
        mu_r=mu_r if bh_curve is None else None, n_nodes=200,
        geometry="cylinder")

    Z_s = esim_solver.solve(5.0)['Z']
    sibc_coeff = 1j * omega / Z_s

    from ngsolve import x, y, z, sqrt, IfPos

    # Kelvin reluctivity
    r_cf = sqrt(x * x + y * y + z * z)
    r_safe = IfPos(r_cf - 1e-12, r_cf, 1e-12)
    if use_kelvin:
        nu_kelvin = nu0 * (r_safe / R_air)**4
        nu_cf = mesh.MaterialCF({"air": nu0, "kelvin": nu_kelvin}, default=nu0)
    else:
        nu_cf = nu0

    fes = HCurl(mesh, order=order, dirichlet=dirichlet_bc, complex=True)
    u, v = fes.TnT()
    print(f"  DOFs: {fes.ndof}")
    print(f"  |Z_s| = {abs(Z_s):.4e}, |sibc_coeff| = {abs(sibc_coeff):.4e}")

    # Bilinear form (same for all iterations)
    # Updated with current Z_s each Karl iteration

    # n_cond = -specialcf.normal (outward from conductor = inward into mesh)
    n_mesh = specialcf.normal(3)
    n_cond = -n_mesh

    # n_cond x H_inc
    nxH = CF((n_cond[1] * H_inc_cf[2] - n_cond[2] * H_inc_cf[1],
              n_cond[2] * H_inc_cf[0] - n_cond[0] * H_inc_cf[2],
              n_cond[0] * H_inc_cf[1] - n_cond[1] * H_inc_cf[0]))

    # Precompute workpiece area
    sibc_region = mesh.Boundaries("sibc")
    A_wp = Integrate(CF(1), mesh, BND, definedon=sibc_region).real
    print(f"  A_wp = {A_wp:.6e} m^2")

    # Karl iteration
    gfu = GridFunction(fes)
    max_iter = 15
    tol = 1e-3
    relax = 0.5
    H_t_rms = 5.0

    print(f"  {'Iter':>4s} {'|Z_s|':>12s} {'P [W]':>12s} "
          f"{'H_t_rms':>10s} {'dZ/Z':>10s}")

    for iteration in range(max_iter):
        sibc_coeff = 1j * omega / Z_s

        # Bilinear form
        a_bf = BilinearForm(fes)
        a_bf += nu_cf * curl(u) * curl(v) * dx
        a_bf += sibc_coeff * u.Trace() * v.Trace() * ds("sibc")
        a_bf.Assemble()

        # RHS: scattered-field formulation
        # f(v) = -a(A_inc, v) = -∫ nu curl(A_inc)·curl(v) dx - (jw/Z_s)∫ A_inc·v ds
        # For boundary form: -∫ nu curl(A_inc)·curl(v) dx = -(n_cond x H_inc)·v ds
        #                     (exact when v=0 on outer boundary, i.e., Dirichlet)
        # For volume form: directly use -H_inc·curl(v) dx (no normal needed)
        f_lf = LinearForm(fes)
        f_lf += (-sibc_coeff) * A_inc_cf * v.Trace() * ds("sibc")
        if rhs_mode == 'volume':
            # Volume form: -∫ H_inc · curl(v) dx  (H_inc = nu0*B_inc)
            B_inc = MU_0 * H_inc_cf
            f_lf += (-nu_cf) * B_inc * curl(v) * dx
        else:
            # Boundary form: -(n_cond x H_inc) · v ds("sibc")
            f_lf += (-1.0) * nxH * v.Trace() * ds("sibc")
        f_lf.Assemble()

        t0 = time.perf_counter()
        with TaskManager():
            gfu.vec.data = a_bf.mat.Inverse(
                fes.FreeDofs(), inverse="pardiso") * f_lf.vec
        t_solve = time.perf_counter() - t0

        # H_t from total field: A_total = A_inc + A_scat (gfu)
        # CRITICAL: use TANGENTIAL component only.
        # |A_t|^2 = |A|^2 - (A . n)^2  where n = specialcf.normal
        A_total = A_inc_cf + gfu
        A_sq = sum(A_total[i].real * A_total[i].real +
                   A_total[i].imag * A_total[i].imag for i in range(3))
        # Normal projection: A.n = sum(A_i * n_i) (complex)
        A_dot_n_re = sum(A_total[i].real * n_mesh[i] for i in range(3))
        A_dot_n_im = sum(A_total[i].imag * n_mesh[i] for i in range(3))
        An_sq = A_dot_n_re * A_dot_n_re + A_dot_n_im * A_dot_n_im
        At_sq = A_sq - An_sq  # tangential |A_t|^2
        At_int = Integrate(At_sq, mesh, BND, definedon=sibc_region).real
        At_rms = math.sqrt(max(At_int, 0) / max(A_wp, 1e-30))
        H_t_rms = abs(1j * omega / Z_s) * At_rms

        # ESIM update
        Z_s_old = Z_s
        sol = esim_solver.solve(max(H_t_rms, 1e-3))
        Z_s_new = sol['Z']
        Z_s = relax * Z_s_new + (1 - relax) * Z_s_old

        P_total = sol['P_prime'] * A_wp
        Q_total = sol['Q_prime'] * A_wp

        dZ = abs(Z_s - Z_s_old) / max(abs(Z_s_old), 1e-30)
        print(f"  {iteration:4d} {abs(Z_s):12.4e} {P_total:12.4e} "
              f"{H_t_rms:10.4f} {dZ:10.4e}")

        if dZ < tol and iteration > 0:
            print(f"  Converged at iteration {iteration}")
            break

    # Inductance from scattered + incident field
    A_total_cf = A_inc_cf + gfu
    curl_A = curl(gfu)  # curl of scattered only (HCurl GF)
    # For total B: need curl(A_inc) + curl(A_scat)
    # curl(A_inc) = B_inc (from Biot-Savart) = mu0 * H_inc
    B_inc = MU_0 * H_inc_cf
    B_total_sq = sum((B_inc[i] + curl_A[i]).real * (B_inc[i] + curl_A[i]).real +
                     (B_inc[i] + curl_A[i]).imag * (B_inc[i] + curl_A[i]).imag
                     for i in range(3))
    W = 0.5 * Integrate(nu0 * B_total_sq, mesh).real
    L = 2 * W / I_total**2

    print()
    print("Results:")
    print("-" * 65)
    print(f"  P (scattered FEM)  = {P_total:.6e} W")
    print(f"  H_t_rms            = {H_t_rms:.4f} A/m")
    print(f"  |Z_s|              = {abs(Z_s):.4e} Ohm")
    print(f"  L (FEM)            = {L*1e9:.2f} nH")
    print(f"  DOFs: {fes.ndof}, Elements: {ne}")
    print(f"  Time: mesh={t_mesh:.1f}s, solve/iter={t_solve:.1f}s")
    print("-" * 65)

    return {
        'P_total': float(P_total),
        'H_t_rms': float(H_t_rms),
        'L': float(L),
        'Z_s': complex(Z_s),
        'ndof': fes.ndof,
        'ne': ne,
    }


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--freq", type=float, default=7000)
    parser.add_argument("--material", default="steel")
    parser.add_argument("--maxh", type=float, default=0.012)
    parser.add_argument("--order", type=int, default=1)
    parser.add_argument("--rhs", default="boundary", choices=["boundary", "volume"])
    parser.add_argument("--maxh-sibc", type=float, default=None,
                        help="Mesh size on workpiece surface (default: R_wp/4)")
    args = parser.parse_args()

    sigma_map = {'steel': 2e6, 'copper': 5.8e7, 'aluminum': 3.5e7}
    sigma = sigma_map.get(args.material, 2e6)
    r = run(frequency=args.freq, material=args.material, sigma=sigma,
            maxh_air=args.maxh, order=args.order, rhs_mode=args.rhs,
            maxh_sibc=args.maxh_sibc)

    if r is None:
        sys.exit(1)

    # Compare with BEM (Scalar BIE SIBC)
    print()
    print("Running Scalar BIE + SIBC (BEM) for comparison...")
    try:
        sys.path.insert(0, os.path.dirname(__file__))
        from pmchwt_sibc import run as bem_run
        r_bem = bem_run(material=args.material, frequency=args.freq, verbose=False)
        print()
        print("=" * 65)
        print(f"Comparison ({args.material}, {args.freq}Hz)")
        print("=" * 65)
        print(f"{'':>20s} {'Scattered FEM':>14s} {'BEM':>14s} {'diff':>8s}")
        P_f, P_b = r['P_total'], r_bem['P_total']
        if abs(P_b) > 1e-30:
            print(f"{'P [W]':>20s} {P_f:14.4e} {P_b:14.4e} "
                  f"{(P_f-P_b)/P_b*100:+8.1f}%")
        print(f"{'H_t_rms [A/m]':>20s} {r['H_t_rms']:14.4f} "
              f"{r_bem['H_t_rms']:14.4f}")
    except Exception as e:
        print(f"BEM comparison failed: {e}")
