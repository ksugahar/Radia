"""
FEM total-field SIBC: gapped torus coil + cylindrical workpiece.

Total-field formulation (no scattered field decomposition):
  curl(nu curl A) + (jw/Z_s) * A_t = J_coil  on wp_surface

Workpiece is a hole (removed from mesh). SIBC on hole boundary.
H_t extraction via pointwise evaluation at air-side surface points
(avoids BND integration trace issue).

Usage:
    python fem_total_field.py
    python fem_total_field.py --freq 7000 --material copper
"""

import math
import os
import sys
import time
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../src/radia'))

MU_0 = 4e-7 * np.pi


def generate_torus_segments(R_major, gap_deg=5, n_per_turn=200):
    """Generate wire segments for a circular coil (same plane z=0)."""
    gap_rad = math.radians(gap_deg)
    theta_max = 2 * math.pi - gap_rad
    segments = []
    for i in range(n_per_turn):
        t1 = theta_max * i / n_per_turn
        t2 = theta_max * (i + 1) / n_per_turn
        segments.append(((R_major * math.cos(t1), R_major * math.sin(t1), 0),
                         (R_major * math.cos(t2), R_major * math.sin(t2), 0)))
    return segments


def compute_Ht_pointwise(mesh, gfu, R_wp, H_wp, Z_s, omega, n_theta=60, n_z=30):
    """Extract H_t from pointwise evaluation on air side of workpiece surface.

    Evaluates A_total at points slightly outside the workpiece,
    computes |A_t| (tangential component), then H_t = |jw/Z_s| * |A_t|.

    Returns:
        H_t_rms: RMS of tangential H over workpiece surface
        detail: dict with per-point data
    """
    eps = R_wp * 0.05  # 5% offset into air

    # Sample points on cylinder lateral surface + caps
    points = []
    normals = []
    areas = []

    # Lateral surface: (R_wp+eps)*cos(theta), sin(theta), z
    dz = H_wp / n_z
    dtheta = 2 * math.pi / n_theta
    for i in range(n_theta):
        theta = 2 * math.pi * (i + 0.5) / n_theta
        ct, st = math.cos(theta), math.sin(theta)
        for j in range(n_z):
            z = -H_wp / 2 + dz * (j + 0.5)
            r = R_wp + eps
            points.append([r * ct, r * st, z])
            normals.append([ct, st, 0])  # outward from cylinder
            areas.append(R_wp * dtheta * dz)  # area element

    # Top cap: r, theta, H_wp/2 + eps
    n_r = 10
    dr = R_wp / n_r
    for i in range(n_theta):
        theta = 2 * math.pi * (i + 0.5) / n_theta
        ct, st = math.cos(theta), math.sin(theta)
        for j in range(n_r):
            r = dr * (j + 0.5)
            points.append([r * ct, r * st, H_wp / 2 + eps])
            normals.append([0, 0, 1])
            areas.append(r * dtheta * dr)

    # Bottom cap
    for i in range(n_theta):
        theta = 2 * math.pi * (i + 0.5) / n_theta
        ct, st = math.cos(theta), math.sin(theta)
        for j in range(n_r):
            r = dr * (j + 0.5)
            points.append([r * ct, r * st, -H_wp / 2 - eps])
            normals.append([0, 0, -1])
            areas.append(r * dtheta * dr)

    points = np.array(points)
    normals = np.array(normals)
    areas = np.array(areas)

    # Evaluate A at each point
    At_sq_vals = []
    for k in range(len(points)):
        p = points[k]
        n_vec = normals[k]
        try:
            mip = mesh(*p)
            A = np.array([complex(gfu[i](mip)) for i in range(3)])
            # Tangential component: A_t = A - n*(A.n)
            A_dot_n = np.dot(A, n_vec)
            A_t = A - A_dot_n * n_vec
            At_sq = float(np.sum(np.abs(A_t)**2))
            At_sq_vals.append(At_sq)
        except Exception:
            At_sq_vals.append(0.0)
            areas[k] = 0.0

    At_sq_vals = np.array(At_sq_vals)
    total_area = np.sum(areas)

    # Debug: check for outliers
    valid = areas > 0
    if np.any(valid):
        vals = At_sq_vals[valid]
        med = np.median(vals)
        # Reject outliers > 100x median
        good = vals < 100 * max(med, 1e-30)
        if not np.all(good):
            n_bad = int(np.sum(~good))
            print(f"  Pointwise: rejecting {n_bad} outliers "
                  f"(median={med:.4e}, max={np.max(vals):.4e})")
            areas_good = areas.copy()
            areas_good[valid] = np.where(good, areas[valid], 0)
            At_sq_vals_good = At_sq_vals.copy()
            At_sq_vals_good[valid] = np.where(good, At_sq_vals[valid], 0)
            total_area = np.sum(areas_good)
            At_sq_avg = np.sum(At_sq_vals_good * areas_good) / max(total_area, 1e-30)
        else:
            At_sq_avg = np.sum(At_sq_vals * areas) / max(total_area, 1e-30)
    else:
        At_sq_avg = 0

    At_rms = math.sqrt(max(At_sq_avg, 0))
    H_t_rms = abs(1j * omega / Z_s) * At_rms

    return H_t_rms, {
        'At_rms': At_rms,
        'total_area': total_area,
        'n_points': len(points),
        'n_valid': int(np.sum(areas > 0)),
    }


def run(R_coil=0.030, a_coil=0.003, gap_deg=5,
        R_wp=0.010, H_wp=0.020,
        sigma=2e6, frequency=7000, material="steel",
        R_air=0.10, maxh_air=0.012, order=1):
    """Total-field FEM-SIBC with Dirichlet outer BC.

    Workpiece removed as hole. J_source from coil volume.
    """
    from ngsolve import (Mesh, HCurl, BilinearForm, LinearForm, GridFunction,
                         Integrate, curl, dx, ds, CF, BND, TaskManager,
                         specialcf)
    from netgen.occ import (WorkPlane, Axes, Axis, Pnt, Dir, Sphere,
                             Cylinder, OCCGeometry, Glue)
    from esim_cell_problem import ESIMFiniteSlabSolver

    I_total = 1.0
    J0 = I_total / (math.pi * a_coil**2)
    omega = 2 * np.pi * frequency
    nu0 = 1.0 / MU_0

    mu_eff = MU_0 * (100.0 if material == "steel" else 1.0)
    delta = math.sqrt(2.0 / (omega * mu_eff * sigma))

    print("=" * 65)
    print("FEM Total-Field SIBC: Coil + Workpiece (hole, Dirichlet)")
    print("=" * 65)
    print(f"Coil:      R={R_coil*1e3:.0f}mm, a={a_coil*1e3:.1f}mm, gap={gap_deg}deg")
    print(f"Workpiece: R={R_wp*1e3:.0f}mm, H={H_wp*1e3:.0f}mm")
    print(f"Material:  {material}, sigma={sigma:.0e}, f={frequency}Hz")
    print(f"Skin depth: delta={delta*1e3:.2f}mm")
    print()

    # === Geometry ===
    print("[1/3] Building geometry...")
    t0 = time.perf_counter()

    # Coil torus (volume with J source)
    wp_plane = WorkPlane(Axes(p=Pnt(R_coil, 0, 0),
                              n=Dir(0, 1, 0), h=Dir(0, 0, 1)))
    circle = wp_plane.Circle(a_coil).Face()
    torus = circle.Revolve(Axis(Pnt(0, 0, 0), Dir(0, 0, 1)), 360 - gap_deg)
    torus.name = "coil"
    torus.maxh = min(a_coil, maxh_air)

    # Air sphere
    air_sphere = Sphere(Pnt(0, 0, 0), R_air)
    for f in air_sphere.faces:
        f.name = "outer"

    # Workpiece cylinder (hole)
    wp_cyl = Cylinder(Pnt(0, 0, -H_wp / 2), Dir(0, 0, 1), R_wp, H_wp)
    for f in wp_cyl.faces:
        f.name = "sibc"
        f.maxh = min(R_wp / 4, maxh_air / 2)

    air = air_sphere - torus - wp_cyl
    air.name = "air"

    shape = Glue([air, torus])
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

    has_sibc = "sibc" in bnds
    if not has_sibc:
        print("  ERROR: 'sibc' boundary not found!")
        return None

    # === FEM-SIBC solve ===
    print("[2/3] FEM-SIBC solve...")

    from ngsolve import x, y, z, sqrt, IfPos
    r_xy = sqrt(x * x + y * y)
    r_safe = IfPos(r_xy - 1e-10, r_xy, 1e-10)
    J_source = J0 * CF((-y / r_safe, x / r_safe, 0))

    # ESIM solver
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

    fes = HCurl(mesh, order=order, dirichlet="outer", complex=True)
    u, v = fes.TnT()
    print(f"  DOFs: {fes.ndof}")
    print(f"  |Z_s| = {abs(Z_s):.4e}")

    # RHS: coil current source
    f_lf = LinearForm(fes)
    f_lf += J_source * v * dx("coil")
    f_lf.Assemble()

    # Workpiece area
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

        a_bf = BilinearForm(fes)
        a_bf += nu0 * curl(u) * curl(v) * dx
        a_bf += sibc_coeff * u.Trace() * v.Trace() * ds("sibc")
        a_bf.Assemble()

        t0 = time.perf_counter()
        with TaskManager():
            gfu.vec.data = a_bf.mat.Inverse(
                fes.FreeDofs(), inverse="pardiso") * f_lf.vec
        t_solve = time.perf_counter() - t0

        # H_t via pointwise evaluation (air side)
        H_t_rms, detail = compute_Ht_pointwise(
            mesh, gfu, R_wp, H_wp, Z_s, omega)

        # Also compute BND integral for comparison
        n_mesh = specialcf.normal(3)
        A_sq = sum(gfu[i].real * gfu[i].real +
                   gfu[i].imag * gfu[i].imag for i in range(3))
        A_dot_n_re = sum(gfu[i].real * n_mesh[i] for i in range(3))
        A_dot_n_im = sum(gfu[i].imag * n_mesh[i] for i in range(3))
        An_sq = A_dot_n_re * A_dot_n_re + A_dot_n_im * A_dot_n_im
        At_sq_bnd = A_sq - An_sq
        At_int_bnd = Integrate(At_sq_bnd, mesh, BND, definedon=sibc_region).real
        At_rms_bnd = math.sqrt(max(At_int_bnd, 0) / max(A_wp, 1e-30))
        H_t_bnd = abs(1j * omega / Z_s) * At_rms_bnd

        if iteration == 0:
            print(f"  DEBUG: H_t (pointwise) = {H_t_rms:.4f}, "
                  f"H_t (BND integral) = {H_t_bnd:.4f}")

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

    # Inductance
    curl_A = curl(gfu)
    B_sq = sum(curl_A[i].real * curl_A[i].real +
               curl_A[i].imag * curl_A[i].imag for i in range(3))
    W = 0.5 * Integrate(nu0 * B_sq, mesh).real
    L = 2 * W / I_total**2

    print()
    print("[3/3] Results")
    print("-" * 65)
    print(f"  P (total-field)    = {P_total:.6e} W")
    print(f"  H_t_rms (ptwise)   = {H_t_rms:.4f} A/m")
    print(f"  H_t_rms (BND int)  = {H_t_bnd:.4f} A/m")
    print(f"  |Z_s|              = {abs(Z_s):.4e} Ohm")
    print(f"  L (FEM)            = {L*1e9:.2f} nH")
    print(f"  DOFs: {fes.ndof}, Elements: {ne}")
    print(f"  Pointwise: {detail['n_valid']}/{detail['n_points']} valid")
    print("-" * 65)

    return {
        'P_total': float(P_total),
        'H_t_rms': float(H_t_rms),
        'H_t_bnd': float(H_t_bnd),
        'L': float(L),
        'Z_s': complex(Z_s),
        'ndof': fes.ndof,
        'ne': ne,
    }


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--freq", type=float, default=7000)
    parser.add_argument("--material", default="copper")
    parser.add_argument("--maxh", type=float, default=0.012)
    parser.add_argument("--order", type=int, default=1)
    args = parser.parse_args()

    sigma_map = {'steel': 2e6, 'copper': 5.8e7, 'aluminum': 3.5e7}
    sigma = sigma_map.get(args.material, 2e6)
    r = run(frequency=args.freq, material=args.material, sigma=sigma,
            maxh_air=args.maxh, order=args.order)

    if r is None:
        sys.exit(1)

    # Compare with BEM
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
        print(f"{'':>20s} {'Total FEM':>14s} {'BEM':>14s} {'diff':>8s}")
        P_f, P_b = r['P_total'], r_bem['P_total']
        if abs(P_b) > 1e-30:
            print(f"{'P [W]':>20s} {P_f:14.4e} {P_b:14.4e} "
                  f"{(P_f-P_b)/P_b*100:+8.1f}%")
        print(f"{'H_t_rms [A/m]':>20s} {r['H_t_rms']:14.4f} "
              f"{r_bem['H_t_rms']:14.4f}")
        if abs(r_bem['H_t_rms']) > 0:
            diff = (r['H_t_rms'] - r_bem['H_t_rms']) / r_bem['H_t_rms'] * 100
            print(f"{'H_t diff':>20s} {'':14s} {'':14s} {diff:+8.1f}%")
    except Exception as e:
        print(f"BEM comparison failed: {e}")
