"""
Test: Analytical SIBC loss from VectorFEMBEM vs ShieldBEMSIBC.

Three loss computation methods:
1. compute_loss() - volume integral P = 0.5*w^2*s*|A_total|^2 (needs fine mesh)
2. compute_loss_sibc() - analytical SIBC using H_inc on boundary faces
3. ShieldBEMSIBC - BEM EFIE with SIBC (includes mutual coupling)

Expected: (2) gives correct order of magnitude but misses edge/corner effects.
ShieldBEMSIBC (3) is the most accurate for thin-skin problems.
"""

import sys
import os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/radia'))

MU_0 = 4e-7 * np.pi


def uniform_B_to_A_inc(B_ext):
    """Create A_inc_func for uniform B field: A = 0.5 * (B x r)."""
    Bx, By, Bz = B_ext

    def A_inc_func(points):
        points = np.atleast_2d(points)
        x, y, z = points[:, 0], points[:, 1], points[:, 2]
        Ax = 0.5 * (By * z - Bz * y)
        Ay = 0.5 * (Bz * x - Bx * z)
        Az = 0.5 * (Bx * y - By * x)
        return np.column_stack([Ax, Ay, Az])

    return A_inc_func


def analytical_loss_half_space(freq, sigma, mu_r, B0, area):
    """Analytical SIBC loss for infinite half-space.

    P = 0.5 * Re(Zs) * |H_t|^2 * area
    For B_ext perpendicular to face normal: H_t = B0 / (mu_r * mu_0)
    """
    omega = 2 * np.pi * freq
    mu = mu_r * MU_0
    delta = np.sqrt(2.0 / (omega * mu * sigma))
    Zs = (1 + 1j) / (sigma * delta)
    H0 = B0 / mu
    return 0.5 * Zs.real * H0**2 * area


def main():
    try:
        from netgen.occ import Box, Pnt, OCCGeometry
        from ngsolve import Mesh, BND
    except ImportError:
        print("SKIP: NGSolve/Netgen not available")
        return

    try:
        from ngsolve.bem import HelmholtzSL
    except ImportError:
        print("SKIP: ngsolve.bem not available")
        return

    from ngbem_eddy import VectorEddyCurrentFEMBEM, ShieldBEMSIBC

    print("=" * 70)
    print("Test: Loss computation methods comparison")
    print("=" * 70)

    # Same geometry as cross-validation
    width, height, depth = 0.02, 0.02, 0.01
    sigma = 3.7e7
    mu_r = 1.0
    maxh = 0.012

    block = Box(Pnt(-width/2, -height/2, -depth/2),
                Pnt(width/2, height/2, depth/2))
    block.solids.name = "conductor"
    block.faces.name = "surface"
    geo = OCCGeometry(block)
    mesh = Mesh(geo.GenerateMesh(maxh=maxh))

    n_vol = sum(1 for _ in mesh.Elements())
    n_bnd = sum(1 for _ in mesh.Elements(BND))
    print(f"\nMesh: {n_vol} vol, {n_bnd} bnd, maxh={maxh*1e3:.0f}mm")

    # Surface area for analytical reference
    # For B_ext = [0,0,1], H_t is nonzero on 4 side faces (not top/bottom)
    side_area = 4 * width * depth  # 4 side faces
    print(f"Side face area (tangential H): {side_area*1e6:.0f} mm^2")

    # --- Build solvers ---
    solver = VectorEddyCurrentFEMBEM(mesh, sigma, mu_r=mu_r, order=1)
    solver.assemble(kappa=0.01, intorder=4)

    sibc = ShieldBEMSIBC(mesh, sigma, mu_r=mu_r)
    sibc.assemble(intorder=4)

    B_ext = [0, 0, 1.0]
    A_inc_func = uniform_B_to_A_inc(B_ext)
    freqs = [1e3, 10e3, 100e3, 1e6]

    print(f"\n{'freq':>10s}  {'delta(mm)':>10s}  {'P_vol(W)':>12s}  "
          f"{'P_analSIBC':>12s}  {'P_shield(W)':>12s}  {'P_halfspace':>12s}  "
          f"{'analSIBC/shld':>13s}")
    print("-" * 100)

    for freq in freqs:
        omega = 2 * np.pi * freq
        delta = np.sqrt(2.0 / (omega * MU_0 * sigma * mu_r))

        # Method 1: Volume integral (requires fine mesh)
        solver.solve(freq, B_ext=B_ext)
        P_vol = solver.compute_loss()

        # Method 2: Analytical SIBC on boundary faces
        P_anal_sibc = solver.compute_loss_sibc()

        # Method 3: ShieldBEMSIBC (BEM + SIBC, includes mutual coupling)
        sibc.solve(freq, A_inc_func=A_inc_func)
        P_shield = sibc.compute_loss()

        # Reference: infinite half-space formula
        P_halfspace = analytical_loss_half_space(
            freq, sigma, mu_r, 1.0, side_area)

        ratio = P_anal_sibc / P_shield if P_shield > 0 else float('inf')

        print(f"{freq:10.0f}  {delta*1e3:10.4f}  {P_vol:12.4e}  "
              f"{P_anal_sibc:12.4e}  {P_shield:12.4e}  {P_halfspace:12.4e}  "
              f"{ratio:13.4f}")

    print(f"\n--- Summary ---")
    print(f"P_vol:       Volume integral. -> 0 on coarse mesh (PEC limit)")
    print(f"P_analSIBC:  Analytical SIBC on boundary faces (half-space approx)")
    print(f"P_shield:    ShieldBEMSIBC (BEM + SIBC with mutual coupling)")
    print(f"P_halfspace: Infinite half-space formula (side faces only)")
    print(f"\nExpected: P_analSIBC ~ P_halfspace (same formula)")
    print(f"Expected: P_shield > P_analSIBC (edge/corner enhancement)")


if __name__ == '__main__':
    main()
