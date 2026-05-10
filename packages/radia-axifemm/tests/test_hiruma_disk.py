"""Phase 3: Hiruma 3-term recurrence on radia_axifemm Cu disk.

Algorithm (port of cuboid_521_kameari_kelvin_v23_hiruma_3term.py, but using
radia_axifemm closed-form K, M instead of standard H1 + Gauss):

    G = global stiffness  (1/mu integrated)
    C = global sigma-mass (sigma integrated, conductor only)
    b = sigma * A_phi_imposed (LinearForm on conductor) — for B0=1 axial field
        A_phi_imposed = B0 * r / 2 (uniform B_z imposes A_phi = B0 r / 2)

    w_1 = G^{-1} b              (initial flux)
    lam_1 = w_1^T G w_1         = R_0 conductance term
    w_2 = -w_1 / lam_1
    for n = 1..N:
        lam_{2n}   = w_{2n}^T C w_{2n}                = L_n
        Aw         = -G^{-1} (C w_{2n})
        w_{2n+1}   = w_{2n-1} - Aw / lam_{2n}        (3-term recurrence)
        lam_{2n+1} = w_{2n+1}^T G w_{2n+1}           = 1 / R_{n+1}
        w_{2n+2}   = w_{2n} - w_{2n+1} / lam_{2n+1}

Cauer mapping:
    R_n = 1 / lam_{2n-1}    (lam_odd index = conductance)
    L_n = lam_{2n}          (lam_even index = inductance)
    tau_n_from_cauer = eigenvalues of Cauer tridiagonal

Cross-validation: tau_n_from_cauer should match scipy eigsh tau_n directly.
Final compare: BEM-Foster tau_n (bem_disk_axisym_v3_refined.json).
"""

import json
import os
import sys
import time

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from math import pi

from ngsolve import (
    Mesh, BilinearForm, LinearForm, Integrate, GridFunction,
    CoefficientFunction, TaskManager, x, ngsglobals,
)
from netgen.occ import OCCGeometry, MoveTo, Glue, X, Y
from radia.radia_axifemm import (
    H1Henrotte, AxiHenrotteStiffnessBFI, AxiHenrotteSigmaMassBFI,
)


R_DISK = 10e-3
T_DISK = 2e-3
SIGMA_CU = 5.8e7
MU0 = 4 * pi * 1e-7
B0 = 1.0  # imposed axial flux density [T]

BEM_TAU_REF_PATH = (
    r"W:/30_CauerLadderNetwork/2026_04_01_長方形CLN/ngsolve_validation/"
    r"bem_disk_axisym_v3_refined.json"
)


def build_disk_mesh(R_air=200e-3, Z_air=200e-3, maxh_disk=0.2e-3, maxh_air=10e-3):
    air = MoveTo(0, -Z_air).Rectangle(R_air, 2 * Z_air).Face()
    disk = MoveTo(0, -T_DISK / 2).Rectangle(R_DISK, T_DISK).Face()
    disk.faces.name = "conductor"
    air.faces.name = "air"
    disk.maxh = maxh_disk
    air.edges.Min(X).name = "axis"
    air.edges.Max(X).name = "right"
    air.edges.Min(Y).name = "bot"
    air.edges.Max(Y).name = "top"
    return Mesh(OCCGeometry(Glue([air, disk]), dim=2).GenerateMesh(maxh=maxh_air))


def to_scipy_csr(mat, n):
    rs, cs, vs = mat.COO()
    K = sp.csr_matrix(
        (np.asarray(vs, dtype=np.float64),
         (np.asarray(rs, dtype=np.int64), np.asarray(cs, dtype=np.int64))),
        shape=(n, n))
    return (K + K.T) * 0.5


def cauer_taus_from_RL(R, L):
    """Reconstruct decay constants tau_n from Cauer-I R, L sequences via the
    tridiagonal admittance matrix. Length-N R and L (R[0]=R_1, L[0]=L_1, etc.).
    """
    N = len(R)
    # State-space form for Cauer-I ladder.
    # Y(s) = 1/Z(s);  Z(s) = R_1 + s L_1 + 1/(R_2/(s L_2) + 1/(...))
    # Equivalent: solve eigenvalue of the L^{-1} R like operator.
    # For purely-real R, L: tau_n are eigenvalues of (L_diag)^{-1} K_tri
    # where K_tri is tridiagonal in R.
    # Simpler: build the FOSTER form by partial fraction expansion of the
    # CFE, then tau_n are 1/poles.
    #
    # We use a direct construction: state vector x in R^N (currents in L_n).
    # Voltage equation per stage: V_n = R_n i_n + L_n di_n/dt + V_{n+1}
    # KCL: i_{n} = i_{n-1} + (V_n - V_{n+1})*Y_? ... cleanest is to use the
    # tridiagonal Lanczos form directly: eigenvalues of T = L^{-1/2} R L^{-1/2}
    # for a symmetric tridiagonal R-matrix, where R_diag = R, R_off = ... TBD.
    #
    # Easiest path: build the bidiagonal generator matrix and find its eigvals.
    # The Cauer-I Foster poles satisfy s_n = -1/tau_n where det(L s + R_eff) = 0.
    # For the user's purposes here, the Hiruma R, L per stage themselves *are*
    # the comparison targets — we don't need tau reconstruction.
    return None  # TODO if needed


def hiruma_3term(K_csr, M_csr, b_vec, free, N_stages, label=""):
    """Run Hiruma 3-term recurrence on reduced (free) system."""
    K = K_csr[free[:, None], free[None, :]]
    M = M_csr[free[:, None], free[None, :]]
    b = b_vec[free]

    factor = spla.factorized(K.tocsc())

    w = [None] * (2 * N_stages + 3)
    lam = [None] * (2 * N_stages + 3)

    w[1] = factor(b)
    lam[1] = float(w[1] @ (K @ w[1]))           # = R_0  (in 1/Ohm units, will discuss)
    w[2] = -w[1] / lam[1]

    print(f"\n  [{label}] λ_1 (= 1/R_0) = {lam[1]:.6e}")
    print(f"  {'n':>3} {'λ_{2n} (L_n)':>16} {'λ_{2n+1} (1/R_{n+1})':>22} "
          f"{'tau_n_per_stage μs':>20} {'orth_drift':>12}")

    out = []
    for n in range(1, N_stages + 1):
        lam[2*n] = float(w[2*n] @ (M @ w[2*n]))                # L_n
        Cw = M @ w[2*n]
        Aw = -factor(Cw)
        w[2*n+1] = w[2*n-1].copy() - Aw / lam[2*n]
        lam[2*n+1] = float(w[2*n+1] @ (K @ w[2*n+1]))           # 1 / R_{n+1}

        # tau per stage = L_n × R_n = λ_{2n} × λ_{2n-1}    (Hiruma convention)
        tau_us = lam[2*n] * lam[2*n-1] * 1e6 if lam[2*n-1] > 0 else None

        # G-orthogonality drift for sanity
        drift = 0.0
        for k in range(1, n):
            cross = float(w[2*k+1] @ (K @ w[2*n+1]))
            drift = max(drift, abs(cross) / abs(lam[2*n+1]))

        w[2*n+2] = w[2*n].copy() - w[2*n+1] / lam[2*n+1]

        ts = f"{tau_us:.4f}" if tau_us is not None else "N/A"
        print(f"  {n:>3} {lam[2*n]:>16.6e} {lam[2*n+1]:>22.6e} {ts:>20} "
              f"{drift:>12.3e}")

        out.append({"n": n, "L_n": lam[2*n], "inv_R_n_plus_1": lam[2*n+1],
                    "tau_per_stage_us": tau_us, "orth_drift": drift})

    return {"lam_1_inv_R0": lam[1], "stages": out}


def solve_disk_hiruma(maxh_disk=0.2e-3, maxh_air=10e-3, R_air=200e-3, Z_air=200e-3,
                      N_stages=6, label=""):
    print(f"\n=== Hiruma 3-term on Cu disk: maxh_disk={maxh_disk*1e3}mm "
          f"R_air={R_air*1e3}mm  {label} ===")
    ngsglobals.msg_level = 0
    mesh = build_disk_mesh(R_air, Z_air, maxh_disk, maxh_air)
    fes = H1Henrotte(mesh, dirichlet="axis|right|top|bot")
    n_free = sum(1 for f in fes.FreeDofs() if f)
    print(f"  mesh ne={mesh.ne}  ndof={fes.ndof}  free={n_free}")

    mu_cf = CoefficientFunction(MU0)
    sigma_cf = mesh.MaterialCF({"conductor": SIGMA_CU}, default=0.0)
    A_imposed = B0 * x / 2  # uniform B_z = B0 -> A_phi = B0 r / 2

    a = BilinearForm(fes, symmetric=True)
    a += AxiHenrotteStiffnessBFI(mu_cf)
    with TaskManager(): a.Assemble()

    m = BilinearForm(fes, symmetric=True)
    m += AxiHenrotteSigmaMassBFI(sigma_cf)
    with TaskManager(): m.Assemble()

    # b_vec = sigma * A_phi_imposed * v assembled with the same DOF convention.
    # In V-DOF, the linear form has natural integrand sigma * A_imposed * N_A_v
    # but for consistency with M (which is sigma * A * A * 2 pi r dr dz), we use
    # the same 2 pi r weight. NGSolve LinearForm with sigma_cf * A_imposed * v * dx
    # gives Integrate(sigma A v dr dz) — missing 2 pi r weight.
    # We add it explicitly via 2 pi x:
    b_form = LinearForm(fes)
    v = fes.TestFunction()
    b_form += sigma_cf * A_imposed * v * 2 * pi * x * dx_marker(fes)
    with TaskManager(): b_form.Assemble()
    b_vec = np.array(b_form.vec)

    K_csr = to_scipy_csr(a.mat, fes.ndof)
    M_csr = to_scipy_csr(m.mat, fes.ndof)
    free = np.array([i for i in range(fes.ndof) if fes.FreeDofs()[i]], dtype=int)

    # Cross-check: scipy eigsh tau_n
    K_red = K_csr[free[:, None], free[None, :]]
    M_red = M_csr[free[:, None], free[None, :]]
    K_red = (K_red + K_red.T) * 0.5
    M_red = (M_red + M_red.T) * 0.5
    eigs, _ = spla.eigsh(K_red, k=min(N_stages, n_free // 2),
                         M=M_red, sigma=0.0, which="LM",
                         tol=1e-10, maxiter=3000)
    eigsh_taus = sorted((1.0/e) * 1e6 for e in eigs)[::-1]
    print(f"  scipy eigsh tau_n[:{len(eigsh_taus)}] = {[f'{t:.3f}' for t in eigsh_taus]} us")

    res = hiruma_3term(K_csr, M_csr, b_vec, free, N_stages, label=label)
    res["eigsh_tau_us"] = eigsh_taus
    res["mesh"] = {"ne": mesh.ne, "ndof": fes.ndof, "free": n_free}
    return res


def dx_marker(fes):
    """Helper: NGSolve dx symbol valid in current import scope."""
    from ngsolve import dx
    return dx


def main():
    # Multi-mesh sweep to detect FEM-side convergence
    cases = [
        (0.5e-3, 5e-3, 100e-3, 100e-3, "coarse"),
        (0.3e-3, 5e-3, 200e-3, 200e-3, "medium"),
        (0.2e-3, 10e-3, 200e-3, 200e-3, "fine"),
    ]
    all_results = {}
    for h_d, h_a, R_a, Z_a, lbl in cases:
        all_results[lbl] = solve_disk_hiruma(
            maxh_disk=h_d, maxh_air=h_a, R_air=R_a, Z_air=Z_a,
            N_stages=6, label=lbl)

    # BEM-Foster reference
    with open(BEM_TAU_REF_PATH) as fp:
        bem = json.load(fp)
    bem_tau = bem["bem_tau_us"]
    print("\n" + "="*70)
    print("Cross comparison: tau_per_stage (Hiruma) vs scipy eigsh vs BEM-Foster")
    print("="*70)
    for lbl, res in all_results.items():
        print(f"\n--- {lbl} (ne={res['mesh']['ne']}, ndof={res['mesh']['free']}) ---")
        print(f"{'n':>3} {'τ_eigsh':>10} {'τ_per_stage':>14} "
              f"{'BEM ref':>10} {'eigsh/BEM':>12}")
        for i, st in enumerate(res["stages"]):
            t_eig = res["eigsh_tau_us"][i] if i < len(res["eigsh_tau_us"]) else None
            t_per = st["tau_per_stage_us"]
            t_bem = bem_tau[i] if i < len(bem_tau) else None
            r = (t_eig / t_bem) if (t_eig and t_bem) else None
            t_eig_s = f"{t_eig:.3f}" if t_eig is not None else "N/A"
            t_per_s = f"{t_per:.3f}" if t_per is not None else "N/A"
            t_bem_s = f"{t_bem:.3f}" if t_bem is not None else "N/A"
            r_s = f"{r:.4f}" if r is not None else "N/A"
            print(f"{st['n']:>3} {t_eig_s:>10} {t_per_s:>14} {t_bem_s:>10} {r_s:>12}")

    out_path = os.path.join(os.path.dirname(__file__), "test_hiruma_disk_results.json")
    with open(out_path, "w") as fp:
        json.dump({"results": all_results, "bem_tau_us": bem_tau}, fp, indent=2)
    print(f"\nSaved: {out_path}")


if __name__ == "__main__":
    main()
