"""
Diagnostic: Why VectorEddyCurrentFEMBEM loss is frequency-independent?
=====================================================================

Dumps internal variables at each frequency to identify the root cause.

Hypothesis:
  A) Missing curl-curl RHS term: f_1 should be -a_FEM*A_inc, not just -jws*M*A_inc
  B) Coarse mesh can't resolve skin layer, A_total = A_scat + A_inc ~ 0
  C) Both effects combined

Part of Radia project
"""

import sys
import os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/radia'))

MU_0 = 4e-7 * np.pi


def main():
    try:
        from netgen.occ import Box, Pnt, OCCGeometry
        from ngsolve import Mesh, BND
    except ImportError:
        print("SKIP: NGSolve/Netgen not available")
        return

    try:
        from ngsolve.bem import HelmholtzSL
        from ngsolve import TaskManager
    except ImportError:
        print("SKIP: ngsolve.bem not available")
        return

    from ngbem_eddy import VectorEddyCurrentFEMBEM

    print("=" * 70)
    print("Diagnostic: VectorFEMBEM frequency-independent loss")
    print("=" * 70)

    # Same mesh as cross-validation
    width, height, depth = 0.02, 0.02, 0.01
    sigma = 3.7e7
    mu_r = 1.0
    maxh = 0.012

    block = Box(Pnt(-width/2, -height/2, -depth/2),
                Pnt(width/2, height/2, depth/2))
    block.solids.name = "conductor"
    block.faces.name = "surface"
    geo = OCCGeometry(block)
    with TaskManager():
        mesh = Mesh(geo.GenerateMesh(maxh=maxh))

        n_vol = sum(1 for _ in mesh.Elements())
        n_bnd = sum(1 for _ in mesh.Elements(BND))
        print(f"\nMesh: {n_vol} vol, {n_bnd} bnd, maxh={maxh*1e3:.0f}mm")

        solver = VectorEddyCurrentFEMBEM(mesh, sigma, mu_r=mu_r, order=1)
        solver.assemble(kappa=0.01, intorder=4)

        B_ext = [0, 0, 1.0]
        freqs = [1e3, 10e3, 100e3, 1e6]

        print(f"\n--- Internal variable dump ---")
        print(f"\n  {'freq':>10s}  {'omega':>12s}  {'delta(mm)':>10s}  "
              f"{'|A_inc|':>12s}  {'|A_scat|':>12s}  {'|A_total|':>12s}  "
              f"{'|j_surf|':>12s}  {'|rho|':>12s}  "
              f"{'w^2*|At|^2':>12s}  {'P_loss(W)':>12s}")
        print("  " + "-" * 140)

        for freq in freqs:
            omega = 2 * np.pi * freq
            delta = np.sqrt(2.0 / (omega * MU_0 * sigma * mu_r))

            result = solver.solve(freq, B_ext=B_ext)

            A_inc = solver._A_inc_coeffs
            A_scat = solver._A_coeffs
            A_total = A_scat + A_inc
            j_surf = solver._j_coeffs
            rho = solver._rho_coeffs

            norm_inc = np.linalg.norm(A_inc)
            norm_scat = np.linalg.norm(A_scat)
            norm_total = np.linalg.norm(A_total)
            norm_j = np.linalg.norm(j_surf)
            norm_rho = np.linalg.norm(rho)

            # Mass-weighted |A_total|^2
            mass_At2 = np.real(A_total.conj() @ solver._a_mass_np @ A_total)
            w2_At2 = omega**2 * mass_At2

            P = solver.compute_loss()

            print(f"  {freq:10.0f}  {omega:12.2f}  {delta*1e3:10.4f}  "
                  f"{norm_inc:12.4e}  {norm_scat:12.4e}  {norm_total:12.4e}  "
                  f"{norm_j:12.4e}  {norm_rho:12.4e}  "
                  f"{w2_At2:12.4e}  {P:12.4e}")

        # --- Diagnosis A: Missing curl-curl RHS term ---
        print(f"\n\n--- Diagnosis A: Missing curl-curl RHS contribution ---")
        print(f"The current code has:  f_1 = -j*w*s * M_mass @ A_inc")
        print(f"Correct should be:    f_1 = -(a_curl + j*w*s*M_mass) @ A_inc")
        print()

        # Check magnitude of the missing term
        A_inc_coeffs = solver._A_inc_coeffs

        for freq in [10e3, 100e3, 1e6]:
            omega = 2 * np.pi * freq
            rhs_mass = np.linalg.norm(omega * sigma * solver._a_mass_np @ A_inc_coeffs)
            rhs_curl = np.linalg.norm(solver._a_curl_np @ A_inc_coeffs)

            print(f"  f={freq/1e3:8.0f} kHz: "
                  f"|w*s*M*A_inc| = {rhs_mass:.4e}, "
                  f"|a_curl*A_inc| = {rhs_curl:.4e}, "
                  f"ratio = {rhs_curl/rhs_mass:.4e}")

        # --- Diagnosis B: What if we add the missing curl-curl term? ---
        print(f"\n\n--- Diagnosis B: Solve WITH curl-curl RHS correction ---")
        print(f"\n  {'freq':>10s}  {'P_orig(W)':>12s}  {'P_fixed(W)':>12s}  "
              f"{'ratio':>8s}")
        print("  " + "-" * 50)

        n1 = solver._n_hcurl
        n2 = solver._n_hdiv_free
        n3 = solver._n_l2
        N = n1 + n2 + n3

        for freq in freqs:
            omega = 2 * np.pi * freq

            # Build system (same as original solve)
            a_FEM = solver._a_curl_np + 1j * omega * sigma * solver._a_mass_np

            A_sys = np.zeros((N, N), dtype=complex)
            A_sys[:n1, :n1] = a_FEM
            A_sys[:n1, n1:n1+n2] = solver._B_np.T
            A_sys[n1:n1+n2, :n1] = solver._B_np
            A_sys[n1:, n1:] = -solver._S_bem_np

            # ORIGINAL RHS (missing curl-curl)
            rhs_orig = np.zeros(N, dtype=complex)
            rhs_orig[:n1] = -1j * omega * sigma * solver._a_mass_np @ A_inc_coeffs
            rhs_orig[n1:n1+n2] = -solver._B_np @ A_inc_coeffs

            # FIXED RHS (with curl-curl)
            rhs_fixed = np.zeros(N, dtype=complex)
            rhs_fixed[:n1] = -(solver._a_curl_np
                               + 1j * omega * sigma * solver._a_mass_np
                               ) @ A_inc_coeffs
            rhs_fixed[n1:n1+n2] = -solver._B_np @ A_inc_coeffs

            # Diagonal scaling (same as original)
            diag_abs = np.abs(np.diag(A_sys))
            diag_abs[diag_abs < 1e-30] = 1.0
            D_scale = np.sqrt(diag_abs)
            D_inv = 1.0 / D_scale

            A_s = A_sys * D_inv[np.newaxis, :] * D_inv[:, np.newaxis]

            # Solve original
            sol_o = np.linalg.solve(A_s, rhs_orig * D_inv) * D_inv
            A_total_o = sol_o[:n1] + A_inc_coeffs
            P_orig = 0.5 * omega**2 * sigma * np.real(
                A_total_o.conj() @ solver._a_mass_np @ A_total_o)

            # Solve fixed
            sol_f = np.linalg.solve(A_s, rhs_fixed * D_inv) * D_inv
            A_total_f = sol_f[:n1] + A_inc_coeffs
            P_fixed = 0.5 * omega**2 * sigma * np.real(
                A_total_f.conj() @ solver._a_mass_np @ A_total_f)

            ratio = P_fixed / P_orig if P_orig > 0 else float('inf')
            print(f"  {freq:10.0f}  {P_orig:12.4e}  {P_fixed:12.4e}  "
                  f"{ratio:8.4f}")

        # --- Diagnosis C: Shielding ratio ---
        print(f"\n\n--- Diagnosis C: Shielding (A_scat canceling A_inc) ---")

        for freq in freqs:
            omega = 2 * np.pi * freq
            solver.solve(freq, B_ext=B_ext)

            A_inc = solver._A_inc_coeffs
            A_scat = solver._A_coeffs
            A_total = A_scat + A_inc

            cancel = np.linalg.norm(A_total) / np.linalg.norm(A_inc)
            cos_angle = np.real(A_scat.conj() @ A_inc) / (
                np.linalg.norm(A_scat) * np.linalg.norm(A_inc) + 1e-30)

            print(f"  f={freq/1e3:8.0f} kHz: "
                  f"|A_total|/|A_inc| = {cancel:.6e}, "
                  f"cos(A_scat, A_inc) = {cos_angle:.6f}")

        # --- Diagnosis D: Compare loss from j (surface current) vs A (volume) ---
        print(f"\n\n--- Diagnosis D: Loss from surface current j ---")
        print(f"  P_Zs = 0.5 * Re(Zs) * j^H * M_hdiv * j")
        print(f"  Compare with P_vol = 0.5 * w^2 * s * A_total^H * M * A_total")

        # We need the HDivSurface mass matrix
        # Approximate: P_j ~ 0.5 * Re(Zs) * |j|^2 * (surface area / n_dofs)
        A_surface = 2 * (width*height + width*depth + height*depth)

        for freq in freqs:
            omega = 2 * np.pi * freq
            delta = np.sqrt(2.0 / (omega * MU_0 * sigma * mu_r))
            Zs = (1 + 1j) / (sigma * delta)

            solver.solve(freq, B_ext=B_ext)
            P_vol = solver.compute_loss()
            j_surf = solver._j_coeffs

            # Quick estimate: |j|^2 * Re(Zs) * scale
            # This is rough but shows the trend
            j_norm2 = np.real(j_surf.conj() @ j_surf)
            P_j_est = 0.5 * Zs.real * j_norm2

            print(f"  f={freq/1e3:8.0f} kHz: "
                  f"P_vol = {P_vol:.4e} W, "
                  f"|j|^2 = {j_norm2:.4e}, "
                  f"Re(Zs) = {Zs.real:.4e}, "
                  f"0.5*Re(Zs)*|j|^2 = {P_j_est:.4e}")


if __name__ == '__main__':
    main()
