"""
ESIM verification against NGSolve 1D FEM (independent method).

Test 1: Linear material - ESIM vs analytical Z_s = rho*gamma*tanh(gamma*a)
Test 2: Linear material - ESIM vs NGSolve H1 FEM (method independence)
Test 3: Nonlinear (steel BH) - ESIM vs NGSolve H1 FEM with Picard iteration
Test 4: H(z) profile comparison (linear + nonlinear)

Key distinction:
  ESIM = finite difference on uniform mesh (tridiagonal solver)
  NGSolve = H1 finite elements, order p=4, Galerkin method
  These are fundamentally different discretizations.

Usage:
    python verify_esim.py
    python verify_esim.py --test 1
"""

import math
import os
import sys
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../src/radia'))
from esim_cell_problem import ESIMFiniteSlabSolver, BHCurveInterpolator

MU_0 = 4e-7 * np.pi

STEEL_BH = [
    [0.0, 0.0], [50.0, 0.10], [100.0, 0.25], [200.0, 0.55],
    [500.0, 0.95], [1000.0, 1.20], [2000.0, 1.40], [5000.0, 1.55],
    [10000.0, 1.65], [20000.0, 1.75], [50000.0, 1.90], [100000.0, 2.00],
]


# ============================================================
# Analytical reference
# ============================================================
def analytical_Zs(xi, rho, a):
    """Z_s = rho * gamma * tanh(gamma * a) for finite slab."""
    ga = complex(1, 1) * xi
    try:
        return (rho / a) * ga * np.tanh(ga)
    except OverflowError:
        return (rho / a) * ga


def analytical_H(z_pts, H0, gamma, a):
    """H(z) = H0 * cosh(gamma*(a-z)) / cosh(gamma*a)."""
    try:
        return H0 * np.cosh(gamma * (a - z_pts)) / np.cosh(gamma * a)
    except (OverflowError, FloatingPointError):
        return H0 * np.exp(-z_pts / abs(1.0 / gamma.real))


# ============================================================
# NGSolve 1D FEM solver
# ============================================================
def build_1d_mesh(a, n_elem):
    """Build 1D NGSolve mesh on [0, a] with named boundaries."""
    from netgen.meshing import (Mesh as NMesh, MeshPoint, Pnt,
                                Element1D, Element0D, FaceDescriptor)
    from ngsolve import Mesh

    ngm = NMesh(dim=1)
    pids = [ngm.Add(MeshPoint(Pnt(a * i / n_elem, 0, 0)))
            for i in range(n_elem + 1)]
    ngm.Add(FaceDescriptor(surfnr=1, domin=1, bc=1))
    ngm.Add(FaceDescriptor(surfnr=2, domin=1, bc=2))
    for i in range(n_elem):
        ngm.Add(Element1D([pids[i], pids[i + 1]], index=1))
    ngm.Add(Element0D(pids[0], index=1))
    ngm.Add(Element0D(pids[-1], index=2))
    ngm.SetBCName(0, 'surface')
    ngm.SetBCName(1, 'center')
    return Mesh(ngm)


def solve_fem_linear(a, sigma, frequency, mu_r, H0, order=4, n_elem=200):
    """NGSolve H1 FEM for linear 1D slab.

    rho*d^2H/dz^2 + j*omega*mu*H = 0, H(0)=H0, dH/dz(a)=0
    """
    from ngsolve import (H1, GridFunction, BilinearForm,
                         grad, dx, CF, Integrate, Conj)

    rho = 1.0 / sigma
    omega = 2 * np.pi * frequency
    mu = MU_0 * mu_r

    mesh = build_1d_mesh(a, n_elem)
    fes = H1(mesh, order=order, complex=True, dirichlet='surface')
    u, v = fes.TnT()

    bf = BilinearForm(fes)
    bf += rho * grad(u) * grad(v) * dx
    bf += (-1j * omega * mu) * u * v * dx
    with TaskManager():
        bf.Assemble()

        gf = GridFunction(fes)
        gf.Set(CF(H0), definedon=mesh.Boundaries('surface'))
        r = -bf.mat * gf.vec
        gf.vec.data += bf.mat.Inverse(fes.FreeDofs()) * r

        P = Integrate(0.5 * rho * grad(gf) * Conj(grad(gf)), mesh).real
        Q = Integrate(0.5 * omega * mu * gf * Conj(gf), mesh).real
        Z = 2 * (P + 1j * Q) / (H0 ** 2)

        # Sample H(z)
        z_pts = np.linspace(0, a, 201)
        H_z = np.array([gf(mesh(z)) for z in z_pts])

        return {'Z': Z, 'P_prime': P, 'Q_prime': Q, 'H_z': H_z, 'z_pts': z_pts}


def solve_fem_nonlinear(a, sigma, frequency, bh_curve, H0,
                        order=4, n_elem=200, max_iter=80, tol=1e-6):
    """NGSolve H1 FEM with Picard iteration for nonlinear BH.

    Iterates: solve with current mu(z) -> update mu from |H(z)|.
    Uses piecewise-constant mu on elements (updated each iteration).
    """
    from ngsolve import (H1, L2, GridFunction, BilinearForm,
                         grad, dx, CF, Integrate, Conj, Parameter)
    from ngsolve import TaskManager

    rho = 1.0 / sigma
    omega = 2 * np.pi * frequency
    bh_interp = BHCurveInterpolator(bh_curve)

    mesh = build_1d_mesh(a, n_elem)
    fes = H1(mesh, order=order, complex=True, dirichlet='surface')
    u, v = fes.TnT()

    # Piecewise-constant mu on L2
    fes_mu = L2(mesh, order=0)
    gf_mu = GridFunction(fes_mu)
    gf_mu.Set(CF(bh_interp.mu_initial))

    gf_H = GridFunction(fes)
    converged = False

    for iteration in range(max_iter):
        bf = BilinearForm(fes)
        bf += rho * grad(u) * grad(v) * dx
        bf += (-1j * omega * gf_mu) * u * v * dx
        bf.Assemble()

        gf_H.Set(CF(H0), definedon=mesh.Boundaries('surface'))
        r = -bf.mat * gf_H.vec
        gf_H.vec.data += bf.mat.Inverse(fes.FreeDofs()) * r

        # Update mu from H(z) at element centers
        mu_old = gf_mu.vec.FV().NumPy().copy()
        mu_new = np.zeros(n_elem)
        for i in range(n_elem):
            z_mid = a * (i + 0.5) / n_elem
            H_val = abs(gf_H(mesh(z_mid)))
            mu_new[i] = bh_interp.mu(max(H_val, 1e-3))

        rel_change = np.max(np.abs(mu_new - mu_old)) / np.max(np.abs(mu_old))

        # Under-relaxation
        gf_mu.vec.FV().NumPy()[:] = 0.5 * mu_old + 0.5 * mu_new

        if rel_change < tol:
            converged = True
            break

    # Power losses
    P = Integrate(0.5 * rho * grad(gf_H) * Conj(grad(gf_H)), mesh).real
    Q = Integrate(0.5 * omega * gf_mu * gf_H * Conj(gf_H), mesh).real
    Z = 2 * (P + 1j * Q) / (H0 ** 2)

    z_pts = np.linspace(0, a, 201)
    H_z = np.array([gf_H(mesh(z)) for z in z_pts])

    return {'Z': Z, 'P_prime': P, 'Q_prime': Q, 'H_z': H_z, 'z_pts': z_pts,
            'converged': converged, 'iterations': iteration + 1}


# ============================================================
# Test 1: Linear - ESIM vs analytical
# ============================================================
def test1_analytical():
    """ESIM Z_s vs analytical Z_s = rho*gamma*tanh(gamma*a)."""
    print("=" * 72)
    print("Test 1: Linear - ESIM vs Analytical Z_s")
    print("=" * 72)

    materials = [("Copper", 5.8e7, 1.0), ("Steel (mu_r=500)", 2.0e6, 500.0)]
    xi_values = [0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0]
    a = 0.005
    H0 = 1000.0

    for name, sigma, mu_r in materials:
        rho = 1.0 / sigma
        print(f"\n  {name}: sigma={sigma:.1e}, mu_r={mu_r}")
        print(f"  {'xi':>5s}  {'|Z| anal':>12s}  {'|Z| ESIM':>12s}  {'err[%]':>8s}")
        print("  " + "-" * 45)
        max_err = 0

        for xi in xi_values:
            delta = a / xi
            f = rho / (np.pi * MU_0 * mu_r * delta**2)
            Z_anal = analytical_Zs(xi, rho, a)

            solver = ESIMFiniteSlabSolver(
                half_thickness=a, sigma=sigma, frequency=f,
                mu_r=mu_r, n_nodes=500)
            sol = solver.solve(H0)
            Z_esim = sol['Z']

            err = (abs(Z_esim) - abs(Z_anal)) / abs(Z_anal) * 100
            max_err = max(max_err, abs(err))
            print(f"  {xi:5.1f}  {abs(Z_anal):12.6e}  {abs(Z_esim):12.6e}  {err:+8.4f}")

        print(f"  Max error: {max_err:.4f}% [{'PASS' if max_err < 1 else 'FAIL'}]")


# ============================================================
# Test 2: Linear - ESIM vs NGSolve FEM
# ============================================================
def test2_fem_linear():
    """ESIM (FD) vs NGSolve H1 FEM for linear material."""
    print("\n" + "=" * 72)
    print("Test 2: Linear - ESIM (FD) vs NGSolve H1 FEM (p=4)")
    print("  Different discretization methods on same equation")
    print("=" * 72)

    a = 0.005
    H0 = 1000.0
    sigma = 5.8e7
    mu_r = 1.0
    rho = 1.0 / sigma

    xi_values = [0.5, 1.0, 2.0, 5.0, 10.0, 20.0]

    print(f"\n  Copper: sigma={sigma:.1e}, mu_r={mu_r}")
    print(f"  {'xi':>5s}  {'|Z| ESIM':>12s}  {'|Z| FEM':>12s}  "
          f"{'err[%]':>8s}  {'P ESIM':>12s}  {'P FEM':>12s}  {'P err[%]':>8s}")
    print("  " + "-" * 80)

    max_err = 0
    for xi in xi_values:
        delta = a / xi
        f = rho / (np.pi * MU_0 * mu_r * delta**2)

        # ESIM
        solver = ESIMFiniteSlabSolver(
            half_thickness=a, sigma=sigma, frequency=f,
            mu_r=mu_r, n_nodes=500)
        sol_esim = solver.solve(H0)

        # NGSolve FEM
        sol_fem = solve_fem_linear(a, sigma, f, mu_r, H0,
                                   order=4, n_elem=200)

        err_Z = (abs(sol_esim['Z']) - abs(sol_fem['Z'])) / abs(sol_fem['Z']) * 100
        err_P = (sol_esim['P_prime'] - sol_fem['P_prime']) / sol_fem['P_prime'] * 100
        max_err = max(max_err, abs(err_Z))

        print(f"  {xi:5.1f}  {abs(sol_esim['Z']):12.6e}  {abs(sol_fem['Z']):12.6e}  "
              f"{err_Z:+8.4f}  {sol_esim['P_prime']:12.6e}  {sol_fem['P_prime']:12.6e}  "
              f"{err_P:+8.4f}")

    print(f"  Max |Z| error: {max_err:.4f}% [{'PASS' if max_err < 1 else 'FAIL'}]")

    # Also test mu_r=500
    sigma_s = 2.0e6
    mu_r_s = 500.0
    rho_s = 1.0 / sigma_s
    print(f"\n  Steel (linear): sigma={sigma_s:.1e}, mu_r={mu_r_s}")
    print(f"  {'xi':>5s}  {'|Z| ESIM':>12s}  {'|Z| FEM':>12s}  "
          f"{'err[%]':>8s}  {'P ESIM':>12s}  {'P FEM':>12s}  {'P err[%]':>8s}")
    print("  " + "-" * 80)

    max_err2 = 0
    for xi in xi_values:
        delta = a / xi
        f = rho_s / (np.pi * MU_0 * mu_r_s * delta**2)

        solver = ESIMFiniteSlabSolver(
            half_thickness=a, sigma=sigma_s, frequency=f,
            mu_r=mu_r_s, n_nodes=500)
        sol_esim = solver.solve(H0)

        sol_fem = solve_fem_linear(a, sigma_s, f, mu_r_s, H0,
                                   order=4, n_elem=200)

        err_Z = (abs(sol_esim['Z']) - abs(sol_fem['Z'])) / abs(sol_fem['Z']) * 100
        err_P = (sol_esim['P_prime'] - sol_fem['P_prime']) / sol_fem['P_prime'] * 100
        max_err2 = max(max_err2, abs(err_Z))

        print(f"  {xi:5.1f}  {abs(sol_esim['Z']):12.6e}  {abs(sol_fem['Z']):12.6e}  "
              f"{err_Z:+8.4f}  {sol_esim['P_prime']:12.6e}  {sol_fem['P_prime']:12.6e}  "
              f"{err_P:+8.4f}")

    print(f"  Max |Z| error: {max_err2:.4f}% [{'PASS' if max_err2 < 1 else 'FAIL'}]")


# ============================================================
# Test 3: Nonlinear - ESIM vs NGSolve FEM + Picard
# ============================================================
def test3_fem_nonlinear():
    """ESIM (FD + Picard) vs NGSolve FEM (H1 p=4 + Picard) for steel BH."""
    print("\n" + "=" * 72)
    print("Test 3: Nonlinear (Steel BH) - ESIM vs NGSolve H1 FEM + Picard")
    print("  Both iterate mu(|H|), but discretization differs (FD vs FEM)")
    print("=" * 72)

    sigma = 2.0e6
    a = 0.002  # 2 mm

    test_cases = [
        (1000,   100.0,  "f=1kHz, H0=100 (linear)"),
        (1000,   5000.0, "f=1kHz, H0=5000 (saturated)"),
        (10000,  100.0,  "f=10kHz, H0=100"),
        (10000,  5000.0, "f=10kHz, H0=5000 (saturated)"),
        (50000,  100.0,  "f=50kHz, H0=100"),
        (50000,  5000.0, "f=50kHz, H0=5000"),
        (100000, 1000.0, "f=100kHz, H0=1000"),
    ]

    print(f"\n  sigma={sigma:.0e}, a={a*1e3:.0f} mm")
    print(f"  {'Case':>32s}  {'|Z| ESIM':>12s}  {'|Z| FEM':>12s}  "
          f"{'err[%]':>8s}  {'P ESIM':>12s}  {'P FEM':>12s}  {'P err[%]':>8s}")
    print("  " + "-" * 100)

    max_err = 0
    for f, H0, desc in test_cases:
        # ESIM
        solver = ESIMFiniteSlabSolver(
            half_thickness=a, bh_curve=STEEL_BH, sigma=sigma,
            frequency=f, n_nodes=500)
        sol_e = solver.solve(H0, tol=1e-8, max_iter=100, relaxation=0.3)

        # NGSolve FEM
        sol_f = solve_fem_nonlinear(a, sigma, f, STEEL_BH, H0,
                                     order=4, n_elem=200, max_iter=80)

        err_Z = (abs(sol_e['Z']) - abs(sol_f['Z'])) / max(abs(sol_f['Z']), 1e-30) * 100
        err_P = (sol_e['P_prime'] - sol_f['P_prime']) / max(sol_f['P_prime'], 1e-30) * 100
        max_err = max(max_err, abs(err_Z))

        cvg = "ok" if sol_f.get('converged', True) else f"nc({sol_f.get('iterations',0)})"
        print(f"  {desc:>32s}  {abs(sol_e['Z']):12.4e}  {abs(sol_f['Z']):12.4e}  "
              f"{err_Z:+8.3f}  {sol_e['P_prime']:12.4e}  {sol_f['P_prime']:12.4e}  "
              f"{err_P:+8.3f}  {cvg}")

    print(f"\n  Max |Z| error: {max_err:.3f}% [{'PASS' if max_err < 5 else 'CHECK'}]")


# ============================================================
# Test 4: H(z) profile comparison
# ============================================================
def test4_profiles():
    """H(z) depth profiles: ESIM vs analytical (linear), ESIM vs FEM (nonlinear)."""
    print("\n" + "=" * 72)
    print("Test 4: H(z) Depth Profile")
    print("=" * 72)

    # Linear: ESIM vs analytical
    sigma = 5.8e7; rho = 1.0/sigma; a = 0.005; H0 = 1000.0
    for xi, label in [(1.0, "xi=1"), (5.0, "xi=5")]:
        delta = a / xi
        f = rho / (np.pi * MU_0 * delta**2)
        gamma = complex(1, 1) / delta

        z_pts = np.linspace(0, a, 501)
        H_anal = analytical_H(z_pts, H0, gamma, a)

        solver = ESIMFiniteSlabSolver(
            half_thickness=a, sigma=sigma, frequency=f,
            mu_r=1.0, n_nodes=501)
        H_esim = solver.solve(H0)['H_solution']

        print(f"\n  Linear copper, {label}:")
        print(f"  {'z/a':>5s}  {'|H| anal':>10s}  {'|H| ESIM':>10s}  {'err[%]':>8s}")
        for zf in [0.0, 0.25, 0.5, 0.75, 1.0]:
            i = int(zf * 500)
            ha, he = abs(H_anal[i]), abs(H_esim[i])
            err = (he - ha) / max(ha, 1e-10) * 100 if ha > 1e-10 else 0
            print(f"  {zf:5.2f}  {ha:10.4f}  {he:10.4f}  {err:+8.4f}")

    # Nonlinear: ESIM vs FEM
    sigma_s = 2e6; a_s = 0.002; f_s = 10000; H0_s = 1000.0
    solver_nl = ESIMFiniteSlabSolver(
        half_thickness=a_s, bh_curve=STEEL_BH, sigma=sigma_s,
        frequency=f_s, n_nodes=501)
    sol_e = solver_nl.solve(H0_s, tol=1e-8, max_iter=100, relaxation=0.3)
    sol_f = solve_fem_nonlinear(a_s, sigma_s, f_s, STEEL_BH, H0_s,
                                 order=4, n_elem=200, max_iter=80)

    print(f"\n  Nonlinear steel BH, f=10kHz, H0=1000:")
    print(f"  {'z/a':>5s}  {'|H| ESIM':>10s}  {'|H| FEM':>10s}  {'err[%]':>8s}")
    z_e = sol_e['mesh_points']
    z_f = sol_f['z_pts']
    for zf in [0.0, 0.25, 0.5, 0.75, 1.0]:
        ie = int(zf * 500)
        iff = int(zf * 200)
        he = abs(sol_e['H_solution'][ie])
        hf = abs(sol_f['H_z'][iff])
        err = (he - hf) / max(hf, 1e-10) * 100 if hf > 1e-10 else 0
        print(f"  {zf:5.2f}  {he:10.4f}  {hf:10.4f}  {err:+8.4f}")


# ============================================================
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="ESIM verification")
    parser.add_argument("--test", type=int, default=0,
                        help="0=all, 1=analytical, 2=FEM linear, 3=FEM nonlinear, 4=profiles")
    args = parser.parse_args()

    tests = {1: test1_analytical, 2: test2_fem_linear,
             3: test3_fem_nonlinear, 4: test4_profiles}

    if args.test == 0:
        for t in tests.values():
            t()
    else:
        tests[args.test]()

    print("\n" + "=" * 72)
    print("Verification complete.")
    print("=" * 72)
