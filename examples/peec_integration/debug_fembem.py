"""Debug: FEM-BEM coupling using SEPARATE spaces + BlockMatrix.

Follows the ngsolve BEM tutorial approach exactly:
- Separate H1 and SurfaceL2 spaces (no product space)
- SingleLayerPotentialOperator / DoubleLayerPotentialOperator / HypersingularOperator
- BlockMatrix assembly
- GMRes solver with block preconditioner

Adapted for diffusion interior (eddy current) + Laplace exterior.
"""

import sys
import os
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/radia'))

MU_0 = 4.0 * np.pi * 1e-7


def main():
    from ngbem_eddy import create_conductor_mesh
    from ngsolve import (H1, SurfaceL2, BilinearForm, LinearForm,
                          GridFunction, TaskManager, ds, dx, grad,
                          specialcf, CF, Integrate, BlockMatrix, BlockVector)
    from ngsolve.bem import (SingleLayerPotentialOperator,
                              DoubleLayerPotentialOperator,
                              HypersingularOperator)
    from ngsolve.solvers import GMRes

    mu_r = 1.0
    sigma = 3.5e7
    depth = 0.002
    width = 0.01
    height = 0.01
    freq = 5000.0
    omega = 2.0 * np.pi * freq
    k_sq = 1j * omega * mu_r * MU_0 * sigma
    Hz_inc = 1.0 / MU_0

    mesh = create_conductor_mesh(width, height, depth, maxh=0.004,
                                  conductor_label="conductor",
                                  surface_label="surface")
    volume = Integrate(CF(1), mesh).real
    print(f"Volume: {volume:.6e} m^3")

    # === SEPARATE function spaces (like ngsolve tutorial) ===
    order = 2
    fesH1 = H1(mesh, order=order, complex=True)
    fesL2 = SurfaceL2(mesh, order=order - 1, complex=True,
                        dual_mapping=True,
                        definedon=mesh.Boundaries("surface"))

    print(f"H1 DOFs: {fesH1.ndof}, L2 DOFs: {fesL2.ndof}")

    u, v = fesH1.TnT()
    uL2, vL2 = fesL2.TnT()

    # === FEM bilinear form (diffusion: stiffness + k^2 mass) ===
    a = BilinearForm(grad(u) * grad(v) * dx + k_sq * u * v * dx)
    a.Assemble()

    # === BEM operators (EXACTLY like tutorial) ===
    print("Assembling BEM operators...")
    with TaskManager():
        V = SingleLayerPotentialOperator(fesL2, intorder=12)
        K = DoubleLayerPotentialOperator(
            fesH1, fesL2,
            trial_definedon=mesh.Boundaries("surface"),
            test_definedon=mesh.Boundaries("surface"),
            intorder=12)
        D = HypersingularOperator(
            fesH1,
            definedon=mesh.Boundaries("surface"),
            intorder=12)
        M = BilinearForm(
            fesH1.TrialFunction() * fesL2.TestFunction().Trace()
            * ds("surface")).Assemble()
    print("  Done.")

    # === Scattered field formulation ===
    # FEM RHS: -k^2 * Hz_inc * (1, v)
    f_H1 = LinearForm(fesH1)
    f_H1 += CF(-k_sq * Hz_inc) * v * dx
    f_H1.Assemble()

    # BEM RHS: 0
    f_L2 = LinearForm(fesL2)
    f_L2.Assemble()

    # === Block system (Costabel symmetric coupling) ===
    # | A + D           (-0.5*M + K)^T |   | u_scat |   | f_H1 |
    # | (-0.5*M + K)    -V             | * | lambda | = | 0    |
    print("Assembling block system...")
    lhs = BlockMatrix([
        [a.mat + D.mat, (-0.5 * M.mat + K.mat).T],
        [(-0.5 * M.mat + K.mat), (-1) * V.mat]
    ])
    rhs = BlockVector([f_H1.vec, f_L2.vec])
    print("  Done.")

    # === Preconditioner (like tutorial) ===
    l2mass = BilinearForm(
        uL2 * vL2.Trace() * ds("surface")).Assemble()
    astab = BilinearForm(
        (grad(u) * grad(v) + 1e-6 * u * v) * dx).Assemble()

    pre = BlockMatrix([
        [astab.mat.Inverse(freedofs=fesH1.FreeDofs(), inverse="pardiso"),
         None],
        [None,
         l2mass.mat.Inverse(freedofs=fesL2.FreeDofs())]
    ])

    # === Solve with GMRes ===
    print("Solving with GMRes...")
    with TaskManager():
        sol = GMRes(A=lhs, b=rhs, pre=pre,
                    tol=1e-8, maxsteps=500, printrates=True)

    gfu = GridFunction(fesH1)
    gfu.vec[:] = sol[0]  # H1 part = Hz_scat

    gfv = GridFunction(fesL2)
    gfv.vec[:] = sol[1]  # L2 part = lambda = dHz_scat/dn

    # === Compute mu_eff ===
    Hz_total_cf = gfu + CF(Hz_inc)
    Hz_avg = Integrate(Hz_total_cf, mesh) / volume
    mu_eff_bem = mu_r * Hz_avg / Hz_inc

    print(f"\nFEM-BEM (Costabel, separate spaces):")
    print(f"  mu_eff = {mu_eff_bem.real:.6f}{mu_eff_bem.imag:+.6f}j")

    # === FEM reference ===
    from ngbem_eddy import EddyCurrentFEMBEM
    eddy_fem = EddyCurrentFEMBEM(mesh, sigma=sigma, mu_r=mu_r, order=order,
                                   conductor_label="conductor",
                                   surface_label="surface")
    eddy_fem.assemble_fem(freq)
    eddy_fem.solve(B_ext=[0, 0, 1.0], mode='fem')
    Hz_avg_fem = Integrate(eddy_fem.get_total_field(), mesh) / volume
    mu_eff_fem = mu_r * Hz_avg_fem / Hz_inc

    mu = mu_r * MU_0
    gamma = np.sqrt(1j * omega * mu * sigma)
    gd2 = gamma * depth / 2.0
    mu_eff_ana = mu_r * np.tanh(gd2) / gd2

    print(f"\nReferences:")
    print(f"  mu_eff_FEM = {mu_eff_fem.real:.6f}{mu_eff_fem.imag:+.6f}j")
    print(f"  mu_eff_ana = {mu_eff_ana.real:.6f}{mu_eff_ana.imag:+.6f}j")
    print(f"\nErrors:")
    print(f"  |BEM - FEM|: {abs(mu_eff_bem - mu_eff_fem):.6f}")
    print(f"  |BEM - ana|: {abs(mu_eff_bem - mu_eff_ana):.6f}")

    # === Also try dense direct solve ===
    print("\n=== Dense direct solve (for comparison) ===")
    n1 = fesH1.ndof
    n2 = fesL2.ndof
    ntot = n1 + n2

    A_dense = np.zeros((ntot, ntot), dtype=complex)

    # Extract blocks manually
    print("Extracting dense blocks...")
    # [H1, H1] block: a + D
    e = GridFunction(fesH1)
    r = GridFunction(fesH1)
    for j in range(n1):
        e.vec[:] = 0
        e.vec[j] = 1.0
        r.vec.data = (a.mat + D.mat) * e.vec
        A_dense[:n1, j] = r.vec.FV().NumPy()

    # [H1, L2] block: (-0.5*M + K)^T
    eL = GridFunction(fesL2)
    for j in range(n2):
        eL.vec[:] = 0
        eL.vec[j] = 1.0
        r.vec.data = (-0.5 * M.mat + K.mat).T * eL.vec
        A_dense[:n1, n1 + j] = r.vec.FV().NumPy()

    # [L2, H1] block: (-0.5*M + K)
    rL = GridFunction(fesL2)
    for j in range(n1):
        e.vec[:] = 0
        e.vec[j] = 1.0
        rL.vec.data = (-0.5 * M.mat + K.mat) * e.vec
        A_dense[n1:, j] = rL.vec.FV().NumPy()

    # [L2, L2] block: -V
    for j in range(n2):
        eL.vec[:] = 0
        eL.vec[j] = 1.0
        rL.vec.data = (-1) * V.mat * eL.vec
        A_dense[n1:, n1 + j] = rL.vec.FV().NumPy()

    # RHS
    rhs_dense = np.zeros(ntot, dtype=complex)
    rhs_dense[:n1] = f_H1.vec.FV().NumPy()
    rhs_dense[n1:] = 0.0

    print("Solving dense system...")
    sol_dense = np.linalg.solve(A_dense, rhs_dense)

    gfu2 = GridFunction(fesH1)
    gfu2.vec.FV().NumPy()[:] = sol_dense[:n1]

    Hz_total_cf2 = gfu2 + CF(Hz_inc)
    Hz_avg2 = Integrate(Hz_total_cf2, mesh) / volume
    mu_eff_dense = mu_r * Hz_avg2 / Hz_inc

    print(f"  mu_eff (dense) = {mu_eff_dense.real:.6f}{mu_eff_dense.imag:+.6f}j")
    print(f"  |dense - FEM|: {abs(mu_eff_dense - mu_eff_fem):.6f}")

    # Check block norms
    print(f"\n  Block norms:")
    print(f"    [H1,H1] = {np.linalg.norm(A_dense[:n1,:n1]):.4e}")
    print(f"    [H1,L2] = {np.linalg.norm(A_dense[:n1,n1:]):.4e}")
    print(f"    [L2,H1] = {np.linalg.norm(A_dense[n1:,:n1]):.4e}")
    print(f"    [L2,L2] = {np.linalg.norm(A_dense[n1:,n1:]):.4e}")
    print(f"    Condition = {np.linalg.cond(A_dense):.2e}")


if __name__ == '__main__':
    main()
