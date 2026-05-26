"""Implement Kuriyama Multi-K + Schur block on the 2D Cu cylinder and
compare against the canonical 1-DOF Schur-augmented CLN3, the pure
rank-10 L-terminated CLN, and the exact Bessel reference.

Multi-K (Kuriyama 2019) construction in modal space:
  Foster eigenmodes phi_n with eigenvalues 1/tau_n; source coefficient g_n;
  set K_0 = I, M = diag(tau_n), b_n = g_n so that
      Y_exact(s) = b^T (K_0 + s M)^{-1} b = sum_n g_n^2 / (1 + s tau_n).
  - Krylov vector at expansion point s_i:
        q_i[n] = (K_0 + s_i M)^{-1} b |_n = g_n / (1 + s_i tau_n).
  - K_0-orthogonalised basis Q (rows = basis vectors in mode space) built
    from {q_0, q_1, ..., q_M} via standard Gram-Schmidt in the K_0 = I
    inner product <u, v> = sum u_n v_n.
  - Reduced K_0_r = Q Q^T = I (M+1)x(M+1)
  - Reduced M_r   = Q diag(tau) Q^T
  - Reduced b_r   = Q b
  - Y_CLN(s) = b_r^T (K_0_r + s M_r)^{-1} b_r.

Schur block (Theorem 1, Sugahara 2026 Paper I):
  Y_R(s) = Y_CLN(s) + K_SIBC * sqrt(s) / (s + d)

Output:
  circle_multik_schur_comparison.pdf  — |Y(jw)| log-log comparison.

References:
  - Kameari et al. 2018 (canonical CLN)
  - Kuriyama et al. 2019 (multi-expansion-point Krylov)
  - Koster, Konig, Biro 2021 (CLN as PGD)
  - Sugahara 2026 Paper I, Theorem 1 (sqrt(s) Schur block)
"""
from __future__ import annotations
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
import matplotlib.font_manager as fm
for _p in (r"C:\Windows\Fonts\times.ttf", r"C:\Windows\Fonts\timesbd.ttf",
           r"C:\Windows\Fonts\timesi.ttf", r"C:\Windows\Fonts\timesbi.ttf"):
    try:
        fm.fontManager.addfont(_p)
    except Exception:
        pass
matplotlib.rcParams["font.family"] = "serif"
matplotlib.rcParams["font.serif"] = ["Times New Roman", "DejaVu Serif"]
matplotlib.rcParams["mathtext.fontset"] = "stix"
matplotlib.rcParams["axes.unicode_minus"] = False
import matplotlib.pyplot as plt
import numpy as np
from scipy.special import iv, jn_zeros

# Cu cylinder parameters
SIGMA = 5.8e7
MU = 4 * math.pi * 1e-7
A = 5e-3
Y_DC = math.pi * A**2 * SIGMA
K_SIBC = 2 * math.pi * A * math.sqrt(SIGMA / MU)
T_char = MU * SIGMA * A**2

# Bessel reference
def Y_exact(s):
    gamma = np.sqrt(s * MU * SIGMA)
    gA = gamma * A
    result = np.empty_like(gA, dtype=complex)
    small = np.abs(gA) < 50
    if np.any(small):
        result[small] = (Y_DC * 2.0 * iv(1, gA[small])
                         / (gA[small] * iv(0, gA[small])))
    large = ~small
    if np.any(large):
        result[large] = Y_DC * 2.0 / gA[large]
    return result


def foster_spectrum(M_modes):
    """Return (x_k, tau_k, g_k^2) for first M modes of the Bessel eigenvalue
    problem on a 2D disk with Dirichlet BC."""
    x = jn_zeros(0, M_modes)             # Bessel J_0 zeros
    tau = MU * SIGMA * A**2 / x**2
    g_sq = Y_DC * 4.0 / x**2
    return x, tau, g_sq


# Pure rank-N Foster sum (L-terminated Cauer-I)
def Y_Foster(s, N):
    _, tau, g_sq = foster_spectrum(N)
    return np.sum(g_sq[None, :] / (1.0 + s[:, None] * tau[None, :]), axis=1)


# --------------------------------------------------------------------
# Kuriyama Multi-K modal-space implementation
# --------------------------------------------------------------------
def multik_basis_modal(M0_order, extra_points, M_modes=400):
    """Build a Kuriyama Multi-K K_0-orthonormal basis in mode space.

    M0_order: number of single-point Krylov vectors at s=0
              {b, K_0^{-1} M b, (K_0^{-1} M)^2 b, ..., (K_0^{-1} M)^{M0-1} b}
              In mode basis with K_0=I, M=diag(tau): {g, tau*g, tau^2*g, ...}.
    extra_points: list of additional s_i (complex) to add as first-order
                  Krylov: q_i[n] = g_n / (1 + s_i tau_n).
    Total basis size = M0_order + len(extra_points).

    Inner product for orthogonalization: K_0 = I  ->  <u, v> = sum u_n v_n.
    """
    _, tau, g_sq = foster_spectrum(M_modes)
    g = np.sqrt(g_sq)
    n_total = M0_order + len(extra_points)
    V = np.zeros((n_total, M_modes), dtype=complex)
    # M0_order single-point Krylov vectors at s=0:
    #   K_0^{-1} = I, so (K_0^{-1} M)^k b = (M^k) b = (tau^k) * g
    for k in range(M0_order):
        V[k, :] = (tau ** k) * g
    # Extra-point first-order Krylov vectors:
    for i, s_i in enumerate(extra_points):
        V[M0_order + i, :] = g / (1.0 + s_i * tau)
    # MGS in K_0 = I inner product
    Q = np.zeros_like(V)
    for i in range(n_total):
        u = V[i].copy()
        for j in range(i):
            u = u - np.vdot(Q[j], u) * Q[j]
        nrm = np.sqrt(np.vdot(u, u).real)
        if nrm < 1e-14:
            raise RuntimeError(
                f"Multi-K basis #{i} numerically dependent on prior basis "
                f"(norm = {nrm:.2e})."
            )
        Q[i] = u / nrm
    return Q


def Y_MultiK(s, M0_order, extra_points, M_modes=400):
    """Y_CLN(s) = b_r^T (K_0_r + s M_r)^{-1} b_r in modal space.

    K_0_r = Q Q^T = I, M_r = Q diag(tau) Q^T, b_r = Q g.
    Basis size = M0_order + len(extra_points).
    """
    _, tau, g_sq = foster_spectrum(M_modes)
    g = np.sqrt(g_sq)
    Q = multik_basis_modal(M0_order, extra_points, M_modes=M_modes)
    n_basis = Q.shape[0]
    K0_r = np.eye(n_basis, dtype=complex)               # = Q Q^T (orthonormal)
    M_r  = (Q * tau[None, :]) @ Q.conj().T              # Q diag(tau) Q^T
    b_r  = Q @ g
    Y = np.empty(s.shape, dtype=complex)
    for k, sv in enumerate(s):
        Asys = K0_r + sv * M_r
        x = np.linalg.solve(Asys, b_r)
        Y[k] = b_r @ x
    return Y


def Y_MultiK_Schur(s, M0_order, extra_points, d, M_modes=400):
    """Multi-K base + sqrt(s) Schur block (Theorem 1)."""
    Y_cln = Y_MultiK(s, M0_order, extra_points, M_modes=M_modes)
    return Y_cln + K_SIBC * np.sqrt(s) / (s + d)


# --------------------------------------------------------------------
# Optimize d for each Multi-K configuration
# --------------------------------------------------------------------
def optimize_d(s, Y_ref, mag_ref, M0_order, extra_points, M_modes=400):
    """Sweep d over log-spaced range, return (best_err, best_d)."""
    best = (1e9, None)
    for d_decade in np.linspace(np.log10(2*math.pi*1e2),
                                 np.log10(2*math.pi*1e10), 81):
        d_try = 10**d_decade
        Y_try = Y_MultiK_Schur(s, M0_order, extra_points, d_try, M_modes=M_modes)
        err = (np.abs(Y_try - Y_ref) / mag_ref).max() * 100
        if err < best[0]:
            best = (err, d_try)
    return best


def main():
    f = np.logspace(0, 8, 240)
    omega = 2 * math.pi * f
    s = 1j * omega

    Y_ex = Y_exact(s)
    mag_ex = np.abs(Y_ex)

    # Configurations
    print("=" * 70)
    print("Configuration comparison on 2D Cu cylinder a = 5 mm")
    print("=" * 70)

    # 1. Pure rank-10 Foster (L-terminated CLN baseline)
    Y_L10 = Y_Foster(s, 10)
    err_L10 = (np.abs(Y_L10 - Y_ex) / mag_ex).max() * 100
    print(f"Pure CLN L-term basis 10  : wall-band err = {err_L10:7.3f}%")

    # 2. Canonical Schur-augmented CLN3 (basis 4 = single-point Pade order 3 + Schur)
    # Single-point Padé[N-1/N] at s=0: M0_order=3 Krylov vectors, no extra points.
    M0_canon = 3
    extra_canon = []
    best_canon = optimize_d(s, Y_ex, mag_ex, M0_canon, extra_canon)
    Y_canon = Y_MultiK_Schur(s, M0_canon, extra_canon, best_canon[1])
    print(f"Schur-aug. CLN3 (basis 4) : wall-band err = {best_canon[0]:7.3f}%  "
          f"d/(2pi) = {best_canon[1]/(2*math.pi):.2e} Hz")

    # 3. Kuriyama Multi-K M_0+4pt + Schur (basis 8)
    # Paper 1 §V.G hybrid: 3 Krylov at s=0 + 4 extra points + 1 Schur DOF
    M0_multik = 3
    extra_multik = [1j * 2 * math.pi * 10**k for k in (3, 4, 5, 6)]
    best_multik = optimize_d(s, Y_ex, mag_ex, M0_multik, extra_multik)
    Y_multik = Y_MultiK_Schur(s, M0_multik, extra_multik, best_multik[1])
    print(f"Multi-K M0=3+4pt + Schur (basis 8): wall-band err = "
          f"{best_multik[0]:7.3f}%  d/(2pi) = {best_multik[1]/(2*math.pi):.2e} Hz")

    # 4. Kuriyama Multi-K M_0+7pt + Schur (basis 11)
    M0_dense = 3
    extra_dense = [1j * 2 * math.pi * 10**k
                   for k in (2.5, 3, 3.5, 4, 4.5, 5, 5.5)]
    best_dense = optimize_d(s, Y_ex, mag_ex, M0_dense, extra_dense)
    Y_dense = Y_MultiK_Schur(s, M0_dense, extra_dense, best_dense[1])
    print(f"Multi-K M0=3+7pt + Schur (basis 11): wall-band err = "
          f"{best_dense[0]:7.3f}%  d/(2pi) = {best_dense[1]/(2*math.pi):.2e} Hz")

    # Asymptote ratios at high f
    print()
    print("C-axis check at f = 10^12 Hz:")
    s_hi = np.array([1j * 2 * math.pi * 1e12])
    K_target_hi = K_SIBC / np.sqrt(2 * math.pi * 1e12)
    for label, Y_func in [
        ("Bessel exact      ", lambda s_: Y_exact(s_)),
        ("CLN L-term  10    ", lambda s_: Y_Foster(s_, 10)),
        ("Schur-aug CLN3 (4)", lambda s_: Y_MultiK_Schur(s_, M0_canon, extra_canon, best_canon[1])),
        ("Multi-K 3+4 (b=8) ", lambda s_: Y_MultiK_Schur(s_, M0_multik, extra_multik, best_multik[1])),
        ("Multi-K 3+7 (b=11)", lambda s_: Y_MultiK_Schur(s_, M0_dense, extra_dense, best_dense[1])),
    ]:
        Y_hi = Y_func(s_hi)[0]
        r = abs(Y_hi) / K_target_hi
        print(f"  {label} : |Y| = {abs(Y_hi):.4e}, r(10^12) = {r:.4f}")

    # ----------------------------------------------------------------
    # PLOT
    # ----------------------------------------------------------------
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.0, 2.8))

    # Panel (a): |Y(jw)|
    ax1.loglog(f, np.abs(Y_ex), color="black", linestyle="-", lw=1.8,
               label=r"$Y_{\mathrm{exact}}$ (Bessel)", zorder=5)
    ax1.loglog(f, np.abs(Y_L10), color="royalblue", linestyle=":", lw=1.3,
               label="Pure CLN L-term basis 10", zorder=3)
    ax1.loglog(f, np.abs(Y_canon), color="darkorange", linestyle="-.",
               lw=1.4, label="Schur-aug. CLN3 (basis 4)", zorder=4)
    ax1.loglog(f, np.abs(Y_multik), color="red", linestyle=(0, (6, 4)),
               lw=1.8, label=r"Multi-K $M_0\!=\!3\!+\!4$pt + Schur (basis 8)",
               zorder=6)

    # SIBC asymptote
    ax1.loglog(f, K_SIBC / np.sqrt(omega), color="gray", linestyle=(0, (1, 3)),
               lw=1.0, label=r"$K_{\mathrm{SIBC}}/\sqrt{\omega}$ asymptote",
               zorder=2)

    ax1.set_xlabel(r"$f$ (Hz)")
    ax1.set_ylabel(r"$|Y(j\omega)|$ (S/m)")
    ax1.set_xlim(1e0, 1e8)
    ax1.grid(True, which="both", linestyle=":", alpha=0.35)
    ax1.legend(loc="lower left", fontsize=6.2, frameon=False,
               handlelength=2.0, labelspacing=0.25)
    ax1.set_title("(a) Admittance comparison", fontsize=9)

    # Panel (b): relative error
    err_canon = np.abs(Y_canon - Y_ex) / mag_ex * 100
    err_multik = np.abs(Y_multik - Y_ex) / mag_ex * 100
    err_dense = np.abs(Y_dense - Y_ex) / mag_ex * 100
    err_L10_arr = np.abs(Y_L10 - Y_ex) / mag_ex * 100

    ax2.loglog(f, err_L10_arr, color="royalblue", linestyle=":",
               lw=1.3, label="Pure CLN L-term basis 10")
    ax2.loglog(f, err_canon, color="darkorange", linestyle="-.",
               lw=1.4, label="Schur-aug. CLN3 (basis 4)")
    ax2.loglog(f, err_multik, color="red", linestyle=(0, (6, 4)),
               lw=1.8, label=r"Multi-K $M_0\!=\!3\!+\!4$pt + Schur (basis 8)")
    ax2.loglog(f, err_dense, color="darkgreen", linestyle="--",
               lw=1.3, label=r"Multi-K $M_0\!=\!3\!+\!7$pt + Schur (basis 11)")
    ax2.axhline(1.0, color="gray", linestyle="-", lw=0.6, alpha=0.5)

    ax2.set_xlabel(r"$f$ (Hz)")
    ax2.set_ylabel(r"Relative error $|Y_{R}-Y_{\rm exact}|/|Y_{\rm exact}|$ (%)")
    ax2.set_xlim(1e0, 1e8)
    ax2.set_ylim(1e-4, 1e3)
    ax2.grid(True, which="both", linestyle=":", alpha=0.35)
    ax2.legend(loc="lower right", fontsize=6.2, frameon=False,
               handlelength=2.0, labelspacing=0.25)
    ax2.set_title("(b) Relative error vs Bessel", fontsize=9)

    fig.tight_layout(pad=0.4)
    out = Path(__file__).parent / "circle_multik_schur_comparison"
    fig.savefig(str(out) + ".pdf")
    fig.savefig(str(out) + ".png", dpi=220)
    print()
    print(f"Wrote {out}.pdf, {out}.png")


if __name__ == "__main__":
    main()
