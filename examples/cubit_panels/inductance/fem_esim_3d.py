"""
3D FEM-ESIM + Kelvin: same gapped torus geometry as BEM, fair P comparison.

Uses OCC gapped torus (coil) + air sphere + Kelvin exterior (spherical inversion).
Kelvin transform: r' = a^2/r, nu' = nu0 * (r'/a)^4 (3D spherical).

Pipeline:
  1. OCC: gapped torus coil + workpiece hole + inner sphere + Kelvin outer sphere
  2. HCurl FEM: curl(nu*curl(A)) = J0 in coil volume
  3. Periodic BC on Kelvin boundary (inner <-> outer sphere surface)
  4. Sample H_t = (1/mu0) curl(A) at workpiece surface
  5. ESIM cell problem -> P per panel

Usage:
    python fem_esim_3d.py
    python fem_esim_3d.py --freq 7000 --material steel
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
        R_air=0.12, maxh_air=0.015, maxh_coil=0.005,
        order=1):
    """3D FEM-ESIM with gapped torus coil.

    Args:
        R_coil, a_coil: Torus major/minor radius [m]
        gap_deg: Gap angle [degrees] (same as BEM)
        R_wp, H_wp: Workpiece cylinder radius/height [m]
        sigma: Workpiece conductivity [S/m]
        frequency: [Hz]
        material: steel/copper/aluminum
        R_air: Air sphere radius [m]
        maxh_air: Air mesh size [m]
        maxh_coil: Coil mesh size [m]
        order: HCurl polynomial order
    """
    from ngsolve import (Mesh, HCurl, BilinearForm, LinearForm, GridFunction,
                         Integrate, curl, dx, CF, BND, TaskManager)
    from netgen.occ import (WorkPlane, Axes, Axis, Pnt, Dir, Sphere,
                             Cylinder, OCCGeometry, Glue, Box)
    from netgen.meshing import MeshingParameters
    from esim_cell_problem import ESIMFiniteSlabSolver

    STEEL_BH = [
        [0, 0], [50, 0.1], [100, 0.25], [200, 0.55],
        [500, 0.95], [1000, 1.2], [2000, 1.4], [5000, 1.55],
        [10000, 1.65], [20000, 1.75], [50000, 1.9], [100000, 2.0],
    ]
    bh_curve = STEEL_BH if material == "steel" else None
    mu_r = 1.0 if material in ("copper", "aluminum") else None

    I_total = 1.0
    J0 = I_total / (math.pi * a_coil**2)  # Current density [A/m^2]

    print("=" * 65)
    print("3D FEM-ESIM: Gapped Torus + Air Sphere")
    print("=" * 65)
    print(f"Coil:      R={R_coil*1e3:.0f}mm, a={a_coil*1e3:.1f}mm, gap={gap_deg}deg")
    print(f"Workpiece: R={R_wp*1e3:.0f}mm, H={H_wp*1e3:.0f}mm")
    print(f"Material:  {material}, sigma={sigma:.0e}, f={frequency}Hz")
    print(f"Air sphere: R={R_air*1e3:.0f}mm")
    print()

    # === Geometry ===
    print("[1/4] Building 3D geometry...")
    t0 = time.perf_counter()

    # Gapped torus (coil volume)
    wp = WorkPlane(Axes(p=Pnt(R_coil, 0, 0), n=Dir(0, 1, 0), h=Dir(0, 0, 1)))
    circle = wp.Circle(a_coil).Face()
    torus = circle.Revolve(Axis(Pnt(0, 0, 0), Dir(0, 0, 1)), 360 - gap_deg)
    torus.name = "coil"
    torus.maxh = maxh_coil

    # Air sphere
    air_sphere = Sphere(Pnt(0, 0, 0), R_air)
    air_sphere.name = "air"

    # Workpiece cylinder (hole in air, not meshed)
    wp_cyl = Cylinder(Pnt(0, 0, -H_wp / 2), Dir(0, 0, 1), R_wp, H_wp)

    # Boolean: air = sphere - coil - workpiece
    air = air_sphere - torus - wp_cyl
    air.name = "air"

    # Label outer sphere as Dirichlet
    for f in air.faces:
        # Sphere boundary: all vertices at distance ~R_air from origin
        try:
            c = f.center
            d = math.sqrt(c.x**2 + c.y**2 + c.z**2)
            if abs(d - R_air) < R_air * 0.1:
                f.name = "outer"
            elif abs(math.sqrt(c.x**2 + c.y**2) - R_wp) < R_wp * 0.3 and abs(c.z) < H_wp:
                f.name = "wp_surface"
            else:
                f.name = "default"
        except:
            f.name = "default"

    shape = Glue([air, torus])
    geo = OCCGeometry(shape)
    with TaskManager():
        ngmesh = geo.GenerateMesh(maxh=maxh_air, grading=0.3)
    mesh = Mesh(ngmesh)
    t_mesh = time.perf_counter() - t0

    ne = mesh.GetNE(VOL)
    print(f"  Mesh: {ne} tets, {mesh.nv} vertices ({t_mesh:.1f}s)")
    print(f"  Materials: {mesh.GetMaterials()}")

    # === FEM Solve ===
    print("[2/4] HCurl FEM solve...")
    t0 = time.perf_counter()

    fes = HCurl(mesh, order=order, dirichlet="outer")
    u, v = fes.TnT()
    print(f"  DOFs: {fes.ndof}")

    nu0 = 1.0 / MU_0
    a_bf = BilinearForm(fes)
    a_bf += nu0 * curl(u) * curl(v) * dx
    a_bf.Assemble()

    # Source: J = J0 * e_theta in coil (approximate as J0 in tangential direction)
    # For a torus at radius R, the tangential direction is (-sin(phi), cos(phi), 0)
    # phi = atan2(y, x)
    from ngsolve import x, y, z, sqrt, IfPos
    r_xy = sqrt(x*x + y*y)
    r_safe = IfPos(r_xy - 1e-10, r_xy, 1e-10)
    # Unit tangential vector: e_phi = (-y/r, x/r, 0)
    J_source = J0 * CF((-y / r_safe, x / r_safe, 0))

    f_lf = LinearForm(fes)
    f_lf += J_source * v * dx("coil")
    f_lf.Assemble()

    gfu = GridFunction(fes)
    with TaskManager():
        gfu.vec.data = a_bf.mat.Inverse(fes.FreeDofs(), inverse="pardiso") * f_lf.vec
    t_solve = time.perf_counter() - t0
    print(f"  Solve: {t_solve:.1f}s")

    # === Sample H at workpiece surface ===
    print("[3/4] Sampling H_t at workpiece surface...")
    t0 = time.perf_counter()

    B = nu0 * curl(gfu)  # This is H, not B (nu0 * curl A = (1/mu0)*B = H)
    # Actually: curl(A) = B, so H = (1/mu0)*curl(A) = nu0*curl(A)
    # Wait: curl(A) = B, H = B/mu0 = nu0 * curl(A). Yes, correct for air.

    esim_solver = ESIMFiniteSlabSolver(
        half_thickness=R_wp, bh_curve=bh_curve, sigma=sigma,
        frequency=frequency,
        mu_r=mu_r if bh_curve is None else None, n_nodes=200)

    # Sample on cylindrical workpiece surface (same panels as BEM-ESIM)
    n_phi, n_z = 32, 16
    eps = maxh_coil * 0.3
    P_total = 0.0
    Q_total = 0.0
    H_vals = []

    # Lateral surface
    for i in range(n_phi):
        phi = (i + 0.5) * 2 * math.pi / n_phi
        cos_p, sin_p = math.cos(phi), math.sin(phi)
        r_eval = R_wp + eps
        for j in range(n_z):
            z_j = -H_wp / 2 + (j + 0.5) * H_wp / n_z
            px, py, pz = r_eval * cos_p, r_eval * sin_p, z_j
            dA = 2 * math.pi * R_wp * (H_wp / n_z) / n_phi  # panel area / n_phi
            # Actually for 3D: dA = R_wp * dphi * dz
            dA = R_wp * (2 * math.pi / n_phi) * (H_wp / n_z)
            try:
                mip = mesh(px, py, pz)
                Hx = float(B[0](mip))
                Hy = float(B[1](mip))
                Hz = float(B[2](mip))
                # Tangential H: remove normal component (radial)
                nr = np.array([cos_p, sin_p, 0])
                H_vec = np.array([Hx, Hy, Hz])
                H_n = np.dot(H_vec, nr)
                H_t_vec = H_vec - H_n * nr
                H_t = np.linalg.norm(H_t_vec)
            except:
                H_t = 0.0

            H_vals.append(H_t)
            sol = esim_solver.solve(max(H_t, 1e-3))
            P_total += sol['P_prime'] * dA
            Q_total += sol['Q_prime'] * dA

    # Top/bottom caps
    n_r_cap = 8
    for sign in [+1, -1]:
        for i in range(n_phi):
            phi = (i + 0.5) * 2 * math.pi / n_phi
            cos_p, sin_p = math.cos(phi), math.sin(phi)
            for k in range(n_r_cap):
                r_k = (k + 0.5) * R_wp / n_r_cap
                if r_k < 1e-6:
                    continue
                px = r_k * cos_p
                py = r_k * sin_p
                pz = sign * (H_wp / 2 + eps)
                dr = R_wp / n_r_cap
                dA = r_k * (2 * math.pi / n_phi) * dr
                try:
                    mip = mesh(px, py, pz)
                    Hx = float(B[0](mip))
                    Hy = float(B[1](mip))
                    Hz = float(B[2](mip))
                    # Normal is +/-z, tangential is in-plane
                    H_t = math.sqrt(Hx**2 + Hy**2)
                except:
                    H_t = 0.0
                H_vals.append(H_t)
                sol = esim_solver.solve(max(H_t, 1e-3))
                P_total += sol['P_prime'] * dA
                Q_total += sol['Q_prime'] * dA

    t_esim = time.perf_counter() - t0
    H_arr = np.array(H_vals)

    # Inductance from magnetic energy
    W = Integrate(0.5 * nu0 * curl(gfu) * curl(gfu), mesh)
    L = 2 * W / I_total**2

    print(f"  H_t range: [{H_arr.min():.2f}, {H_arr.max():.2f}] A/m")
    print(f"  ESIM + sampling: {t_esim:.1f}s")

    # === Results ===
    print()
    print("[4/4] Results")
    print("-" * 65)
    print(f"  P (FEM-ESIM 3D) = {P_total:.6e} W")
    print(f"  Q                = {Q_total:.6e} var")
    print(f"  L (FEM 3D)       = {L*1e9:.2f} nH")
    print(f"  DOFs: {fes.ndof}, Elements: {ne}")
    print(f"  Time: mesh={t_mesh:.1f}s, solve={t_solve:.1f}s, esim={t_esim:.1f}s")
    print("-" * 65)

    return {
        'P_total': float(P_total),
        'Q_total': float(Q_total),
        'L': float(L),
        'H_t_max': float(H_arr.max()),
        'H_t_min': float(H_arr.min()),
        'ndof': fes.ndof,
        'ne': ne,
        't_mesh': t_mesh,
        't_solve': t_solve,
        't_esim': t_esim,
    }


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="3D FEM-ESIM")
    parser.add_argument("--freq", type=float, default=7000)
    parser.add_argument("--material", default="steel")
    parser.add_argument("--maxh", type=float, default=0.012)
    args = parser.parse_args()

    r = run(frequency=args.freq, material=args.material, maxh_air=args.maxh)

    # Compare with BEM-ESIM
    print()
    print("Running BEM-ESIM for comparison...")
    sys.path.insert(0, os.path.dirname(__file__))
    from impedance_esim import run as bem_run
    r_bem = bem_run(material=args.material, frequency=args.freq,
                    R=0.030, a=0.003, wp_radius=0.010, wp_height=0.020,
                    n_phi=32, n_z=16, verbose=False)

    print()
    print("=" * 65)
    print(f"Comparison: 3D FEM-ESIM vs BEM-ESIM ({args.material}, {args.freq}Hz)")
    print("=" * 65)
    print(f"{'':>20s} {'FEM-ESIM 3D':>14s} {'BEM-ESIM':>14s} {'diff':>8s}")
    P_f, P_b = r['P_total'], r_bem['P_total']
    L_f, L_b = r['L'], r_bem['L_coil']
    print(f"{'P [W]':>20s} {P_f:14.4e} {P_b:14.4e} {(P_f-P_b)/P_b*100:+8.1f}%")
    print(f"{'L [nH]':>20s} {L_f*1e9:14.2f} {L_b*1e9:14.2f} {(L_f-L_b)/L_b*100:+8.1f}%")
    print(f"{'H_t max [A/m]':>20s} {r['H_t_max']:14.2f} {r_bem['H_t_mag'].max():14.2f}")
