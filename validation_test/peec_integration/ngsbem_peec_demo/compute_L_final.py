"""
Compute loop inductance of rectangular frame using ngsbem EFIE.

Uses the EXACT ngsbem patterns from the official demos:
  V1 = HelmholtzSL(u.Trace()*ds, kappa) * v.Trace() * ds   [vector SL = V_A]
  V2 = HelmholtzSL(div(u.Trace())*ds, kappa) * div(v.Trace()) * ds  [scalar SL = V_Phi]
  V_EFIE = kappa * V1 - (1/kappa) * V2

For div-free harmonic c_h: V_Phi(c_h) = 0, so V_A(c_h) = V1(c_h).
L = mu_0 * c_h^T * V1 * c_h / I^2

Reference: FastHenry L ~ 21.6 nH (3D cross-section), Grover ~24 nH (single filament)
"""
import numpy as np
from scipy.linalg import null_space

MU_0 = 4.0 * np.pi * 1e-7
WIDTH = 0.01
TRACE_W = 1e-3
hw = TRACE_W / 2.0
SIGMA = 5.8e7
THICKNESS = 35e-6


def extract_dense(mat, n):
    ei = mat.CreateColVector()
    col = mat.CreateColVector()
    M = np.zeros((n, n), dtype=complex)
    for i in range(n):
        ei[:] = 0; ei[i] = 1.0
        mat.Mult(ei, col)
        for j in range(n):
            M[j, i] = col[j]
    return M


def extract_rect(mat, nrow, ncol):
    ei = mat.CreateRowVector()
    col = mat.CreateColVector()
    M = np.zeros((nrow, ncol))
    for i in range(ncol):
        ei[:] = 0; ei[i] = 1.0
        mat.Mult(ei, col)
        for j in range(nrow):
            M[j, i] = col[j]
    return M


def main():
    from netgen.occ import WorkPlane, OCCGeometry, Axes, Pnt, Dir
    from ngsolve import (Mesh, HDivSurface, SurfaceL2, BilinearForm,
                         ds, BND, TaskManager, InnerProduct)
    from ngsolve import div as ng_div
    from ngsolve.bem import HelmholtzSL

    # --- Frame geometry ---
    wp_o = WorkPlane(Axes(Pnt(-hw, -hw, 0), Dir(0, 0, 1), Dir(1, 0, 0)))
    outer = wp_o.Rectangle(WIDTH + TRACE_W, WIDTH + TRACE_W).Face()
    wp_i = WorkPlane(Axes(Pnt(hw, hw, 0), Dir(0, 0, 1), Dir(1, 0, 0)))
    inner = wp_i.Rectangle(WIDTH - TRACE_W, WIDTH - TRACE_W).Face()
    frame = outer - inner
    frame.faces.name = "conductor"
    geo = OCCGeometry(frame)

    for maxh in [0.001]:
        print(f"\n{'='*60}")
        print(f"maxh = {maxh*1000:.1f} mm")
        print(f"{'='*60}")

        mesh = Mesh(geo.GenerateMesh(maxh=maxh))
        n_el = sum(1 for _ in mesh.Elements(BND))
        fes_J = HDivSurface(mesh, order=0)
        n_J = fes_J.ndof
        n_v = mesh.nv
        fes_L2 = SurfaceL2(mesh, order=0)
        n_f = fes_L2.ndof
        print(f"  {n_el} el, {n_v} vert, {n_J} edges, {n_f} faces")

        # === D matrix (divergence) ===
        u_J = fes_J.TrialFunction()
        q_L2 = fes_L2.TestFunction()
        bf_D = BilinearForm(trialspace=fes_J, testspace=fes_L2)
        bf_D += ng_div(u_J.Trace()) * q_L2 * ds
        bf_D.Assemble()
        D = extract_rect(bf_D.mat, n_f, n_J)

        # === C matrix (surface curl = signed incidence) ===
        C = np.zeros((n_J, n_v))
        for e_idx, edge in enumerate(mesh.edges):
            verts = list(edge.vertices)
            C[e_idx, verts[0].nr] = -1
            C[e_idx, verts[1].nr] = +1

        # === M_J (mass matrix) ===
        u2, v2 = fes_J.TnT()
        bf_M = BilinearForm(fes_J)
        bf_M += InnerProduct(u2.Trace(), v2.Trace()) * ds
        bf_M.Assemble()
        M_J = np.real(extract_dense(bf_M.mat, n_J))

        # === Harmonic (M_J-orthogonal Hodge decomposition) ===
        constraint = np.vstack([D, C.T @ M_J])
        null_h = null_space(constraint, rcond=1e-10)
        c_h = null_h[:, 0]
        energy = c_h @ M_J @ c_h
        print(f"  ||D*c_h|| = {np.linalg.norm(D@c_h):.1e}")
        print(f"  ||C^T*M*c_h|| = {np.linalg.norm(C.T@M_J@c_h):.1e}")
        print(f"  energy = {energy:.6e}")

        # === V_A via HelmholtzSL (demo pattern, with bonus_intorder) ===
        # V1 = HelmholtzSL(u.Trace()*ds, kappa) * v.Trace() * ds
        # This is the PURE vector SL (no V_Phi)
        print(f"\n  --- V_A via HelmholtzSL (ngsbem demo pattern) ---")
        kappa_test = 0.01  # Small kappa -> Laplace limit
        u3, v3 = fes_J.TnT()
        with TaskManager():
            V1_op = HelmholtzSL(
                u3.Trace() * ds("conductor", bonus_intorder=6),
                kappa_test
            ) * v3.Trace() * ds("conductor", bonus_intorder=6)

        V1_mat = extract_dense(V1_op.mat, n_J)
        V_A_helmholtz = np.real(c_h @ V1_mat @ c_h)
        print(f"  V_A (HelmholtzSL, k={kappa_test}) = {V_A_helmholtz:.6e}")

        # Verify: V_Phi should be ~0 for div-free c_h
        with TaskManager():
            V2_op = HelmholtzSL(
                ng_div(u3.Trace()) * ds("conductor", bonus_intorder=6),
                kappa_test
            ) * ng_div(v3.Trace()) * ds("conductor", bonus_intorder=6)
        V2_mat = extract_dense(V2_op.mat, n_J)
        V_Phi_check = np.real(c_h @ V2_mat @ c_h)
        print(f"  V_Phi check = {V_Phi_check:.6e} (should be ~0)")

        # === Full EFIE check: V = kappa * V1 - (1/kappa) * V2 ===
        V_EFIE_proj = kappa_test * V_A_helmholtz - (1/kappa_test) * V_Phi_check
        print(f"  V_EFIE proj = {V_EFIE_proj:.6e}")
        print(f"  V_EFIE/kappa = {V_EFIE_proj/kappa_test:.6e} (should = V_A)")

        # === Inductance ===
        R_sheet = 1.0 / (SIGMA * THICKNESS)
        perimeter = 4 * WIDTH
        R_loop = R_sheet * perimeter / TRACE_W

        # L/R ratio (normalization-independent)
        LR_ratio = MU_0 * V_A_helmholtz / (R_sheet * energy)
        L_val = LR_ratio * R_loop

        print(f"\n  --- Results ---")
        print(f"  R_sheet = {R_sheet:.4f} Ohm/sq")
        print(f"  R_loop = {R_loop*1e3:.2f} mOhm")
        print(f"  L/R = {LR_ratio*1e6:.4f} us")
        print(f"  L = {L_val*1e9:.2f} nH")
        print(f"  FastHenry: L=21.6 nH, R=19.2 mOhm (3D cross-section)")
        print(f"  Grover: L~24 nH (single filament, same GMD)")

    print("\nDone.")


if __name__ == '__main__':
    main()
