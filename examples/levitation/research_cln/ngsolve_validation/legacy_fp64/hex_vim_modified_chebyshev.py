"""hex_vim_modified_chebyshev.py — Gautschi's Modified Chebyshev algorithm
for VIM Cauer extraction.

Strategy: replace POWER moments alpha_n = sum g^2 tau^n (which suffer
exponential dynamic range and cancellation) with MODIFIED moments

  m_k = sum_j g_j^2 pi_k(tau_j)

where {pi_k} are shifted-Chebyshev polynomials on a safe interval [c, d]
containing the Foster support. Since |pi_k(tau)| is bounded on [c, d],
modified moments do not blow up exponentially — much better conditioning.

Wheeler's algorithm then converts modified moments + the reference
recurrence (alpha_k, beta_k of shifted Chebyshev) to the Stieltjes
recurrence (a_n, b_n) of orthogonal polynomials w.r.t. the discrete
measure mu = sum g^2 delta(tau - tau_j).

Once we have (a_n, b_n), we reconstruct power moments alpha_n in a
NUMERICALLY STABLE way (no direct exponentials of tau), then apply
QD-Padé to get Cauer-I rungs.

Reference: W. Gautschi, "Orthogonal Polynomials: Computation and
Approximation" (Oxford 2004), Algorithm 2.45 (Modified Chebyshev).
Also J.C. Wheeler, "Modified moments and Gaussian quadrature", Rocky
Mt. J. Math. 4 (1974) 287-296.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")

import mpmath as mp
mp.mp.dps = 80


def shifted_chebyshev_recurrence(c, d, n_stages):
    """Monic shifted Chebyshev T*_n on [c, d]:
    T*_{n+1}(tau) = (tau - alpha_n) T*_n(tau) - beta_n T*_{n-1}(tau)

    Standard monic Chebyshev on [-1, 1]: alpha = 0, beta_0 = 2, beta_n = 1/4 for n >= 1.
    (Note: beta_0 = <T_0, T_0> = pi for orthogonal-but-not-monic; for monic
    beta_0 = 1, but the inner product weight differs. We just need the
    recurrence coefficients here.)

    Shifted to [c, d] via x = (2 tau - c - d) / (d - c) so tau = (c+d)/2 + (d-c)/2 x:
    alpha_n* = (c + d) / 2
    beta_0* = arbitrary (only used for normalization)
    beta_1* = (d - c)^2 / 8
    beta_n* = (d - c)^2 / 16 for n >= 2

    Returns: alpha (list), beta (list)
    """
    midpoint = (c + d) / 2
    half_width = (d - c) / 2

    alpha = [midpoint] * n_stages
    beta = [mp.mpf(1)]   # beta_0 conventional 1
    if n_stages >= 2:
        beta.append(half_width ** 2 / 2)   # beta_1 = (d-c)^2 / 8 = half_width^2 / 2
    for _ in range(2, n_stages):
        beta.append(half_width ** 2 / 4)   # beta_n = (d-c)^2 / 16 = half_width^2 / 4
    return alpha, beta


def chebyshev_value(tau, c, d, k_max):
    """Evaluate monic shifted Chebyshev polynomials T*_k(tau) for k=0..k_max
    via the recurrence. tau can be scalar or list of mp.mpf.

    Returns list of mp.mpf with k_max+1 entries: T*_0, T*_1, ..., T*_{k_max}.
    """
    midpoint = (c + d) / 2
    half_width = (d - c) / 2

    p_prev = mp.mpf(0)
    p_curr = mp.mpf(1)   # T*_0 = 1
    out = [p_curr]
    for k in range(k_max):
        # T*_{k+1}(tau) = (tau - alpha_k) T*_k - beta_k T*_{k-1}
        alpha_k = midpoint
        if k == 0:
            beta_k = mp.mpf(0)  # T*_1 = tau - midpoint (no beta term)
        elif k == 1:
            beta_k = half_width ** 2 / 2
        else:
            beta_k = half_width ** 2 / 4
        p_next = (tau - alpha_k) * p_curr - beta_k * p_prev
        out.append(p_next)
        p_prev = p_curr
        p_curr = p_next
    return out


def modified_moments(tau, g2, c, d, n_moments):
    """Compute m_k = sum_j g_j^2 T*_k(tau_j) for k = 0..n_moments-1."""
    m = [mp.mpf(0)] * n_moments
    for j in range(len(tau)):
        cheb_vals = chebyshev_value(tau[j], c, d, n_moments - 1)
        for k in range(n_moments):
            m[k] = m[k] + g2[j] * cheb_vals[k]
    return m


def wheeler_algorithm(modified_m, alpha_ref, beta_ref, n_stages):
    """Convert modified moments m_k = sum g^2 pi_k(tau) + reference
    recurrence (alpha_ref, beta_ref) of pi_k to the desired Stieltjes
    recurrence (a_n, b_n) of orthogonal polynomials w.r.t. dmu.

    Wheeler's algorithm (Gautschi Algorithm 2.45 / Wheeler 1974):

      sigma[-1, k] = 0  (k = 0, 1, ..., 2n-1)
      sigma[ 0, k] = m_k  (k = 0, 1, ..., 2n-1)
      a_0 = alpha_ref[0] + m_1 / m_0
      b_0 = m_0

      For n = 1, 2, ..., n_stages-1:
        sigma[n, k] = sigma[n-1, k+1]
                     - (a[n-1] - alpha_ref[k]) * sigma[n-1, k]
                     - b[n-1] * sigma[n-2, k]
                     + beta_ref[k] * sigma[n-1, k-1]
                       (for k = n, n+1, ..., 2n-n-1)

        a[n] = alpha_ref[n] + sigma[n, n+1] / sigma[n, n]
                            - sigma[n-1, n] / sigma[n-1, n-1]
        b[n] = sigma[n, n] / sigma[n-1, n-1]

    Returns (a, b) as lists of length n_stages.
    """
    N = len(modified_m)
    assert N >= 2 * n_stages, f"need 2*n_stages={2*n_stages} moments, have {N}"

    # sigma[n, k] indexed; store as 2D list
    sigma = [[mp.mpf(0)] * N for _ in range(n_stages + 1)]
    # Row 0 (sigma_{-1}) already zeros
    # Row 1 corresponds to sigma_0; sigma_0[k] = m_k
    for k in range(N):
        sigma[1][k] = modified_m[k]

    a = []
    b = []
    a.append(alpha_ref[0] + sigma[1][1] / sigma[1][0])
    b.append(sigma[1][0])   # b_0 = m_0 = <p_0, p_0>

    for n in range(1, n_stages):
        # Compute sigma[n+1, k] for k = n, n+1, ..., 2*n_stages - n - 1
        # Using formula:
        # sigma[n+1, k] = sigma[n, k+1] - (a[n-1] - alpha_ref[k]) sigma[n, k]
        #                              - b[n-1] sigma[n-1, k]
        #                              + beta_ref[k] sigma[n, k-1]
        for k in range(n, 2 * n_stages - n):
            term1 = sigma[n][k + 1] if k + 1 < N else mp.mpf(0)
            term2 = (a[n - 1] - alpha_ref[k]) * sigma[n][k]
            term3 = b[n - 1] * sigma[n - 1][k]
            term4 = beta_ref[k] * sigma[n][k - 1]
            sigma[n + 1][k] = term1 - term2 - term3 + term4

        a_n = (alpha_ref[n]
               + sigma[n + 1][n + 1] / sigma[n + 1][n]
               - sigma[n][n] / sigma[n][n - 1])
        b_n = sigma[n + 1][n] / sigma[n][n - 1]
        a.append(a_n)
        b.append(b_n)

    return a, b


def alpha_from_jacobi(a, b, n_moments):
    """Reconstruct power moments alpha_n = sum g^2 tau^n from Stieltjes
    (a, b) of monic OG polynomials. Uses recurrence:
       p_{n+1}(0) = -a_n p_n(0) - b_n p_{n-1}(0)
    NOT what we want. Instead reconstruct from Gauss quadrature:

    Better: use the Jacobi matrix J = tridiag(a, sqrt(b[1:])), compute
    eigenvalues lambda_k and weights w_k = b_0 * (first component of
    eigenvector_k)^2. Then alpha_n = sum_k w_k lambda_k^n.

    For our purposes (extracting Cauer rungs via QD-Padé), we want
    high-precision alpha_n that match the original Foster spectrum but
    come through the stable (a, b) chain.
    """
    import numpy as np

    N = len(a)
    # Build Jacobi matrix (mpmath floats converted)
    diag = [float(ai) for ai in a]
    offdiag = [float(mp.sqrt(b[k])) for k in range(1, N)]

    # Convert to numpy for eigh_tridiagonal
    from scipy.linalg import eigh_tridiagonal
    eigvals, eigvecs = eigh_tridiagonal(np.array(diag), np.array(offdiag))

    # Weights: w_k = b_0 * (first component of normalized eigenvector_k)^2
    b0_f = float(b[0])
    weights = b0_f * eigvecs[0, :] ** 2

    # alpha_n = sum_k w_k lambda_k^n
    # Use mpmath for high precision
    eigvals_mp = [mp.mpf(ev) for ev in eigvals]
    weights_mp = [mp.mpf(w) for w in weights]

    alphas = []
    for n in range(n_moments):
        a_val = mp.mpf(0)
        for w, l in zip(weights_mp, eigvals_mp):
            a_val = a_val + w * (l ** n)
        alphas.append((mp.mpf(-1)) ** n * a_val)
    return alphas


def qd_pade_cauer(alphas, n_stages):
    c = [mp.mpf(a) for a in alphas]
    p = []
    for _ in range(n_stages):
        if len(c) < 2 or abs(c[0]) < mp.mpf(10) ** (-mp.mp.dps + 10):
            break
        n_c = len(c)
        e = [mp.mpf(0)] * n_c
        e[0] = 1 / c[0]
        for k in range(1, n_c):
            e[k] = -sum(c[j + 1] * e[k - j - 1] for j in range(k)
                        if k - j - 1 >= 0) / c[0]
        p.append(e[0])
        c = e[1:]
    return p


def main():
    if len(sys.argv) < 2:
        json_path = (Path("S:/Radia/01_GitHub/src/ext/radia_vim/scripts")
                     / "cuboid521_order4_results.json")
    else:
        json_path = Path(sys.argv[1])
    print(f"=== Modified Chebyshev (Gautschi/Wheeler) Cauer extraction ===")
    print(f"  Foster JSON: {json_path}")

    data = json.loads(json_path.read_text(encoding="utf-8"))
    modes = data["top_modes_by_g2"]
    n_use = len([m for m in modes if abs(m["g2"]) > 1e-32])
    print(f"  Using {n_use} Foster modes")

    tau = [mp.mpf(m["tau_us"]) * mp.mpf("1e-6") for m in modes[:n_use]]
    g2 = [mp.mpf(m["g2"]) for m in modes[:n_use]]

    tau_min = min(tau)
    tau_max = max(tau)
    print(f"  tau range: [{float(tau_min)*1e6:.4e}, {float(tau_max)*1e6:.4e}] us")

    # Choose Chebyshev support [c, d] with small margin
    c_supp = tau_min * mp.mpf("0.5")
    d_supp = tau_max * mp.mpf("1.5")
    print(f"  Chebyshev support: [c={float(c_supp)*1e6:.4e}, "
          f"d={float(d_supp)*1e6:.4e}] us")

    n_stages = 10
    n_moments = 2 * n_stages + 4

    # Modified moments
    print(f"  Computing {n_moments} modified Chebyshev moments...")
    m_mod = modified_moments(tau, g2, c_supp, d_supp, n_moments)
    print(f"  m_0 = {float(m_mod[0]):.4e}, m_5 = {float(m_mod[5]):.4e}, "
          f"m_{n_moments-1} = {float(m_mod[-1]):.4e}")
    print(f"  (Note: modified moments are BOUNDED, no exponential growth)")

    # Reference recurrence (shifted Chebyshev)
    alpha_ref, beta_ref = shifted_chebyshev_recurrence(c_supp, d_supp,
                                                       n_moments)

    # Wheeler's algorithm
    print(f"\n  Wheeler's modified Chebyshev algorithm -> Stieltjes (a, b)...")
    a_wheeler, b_wheeler = wheeler_algorithm(m_mod, alpha_ref, beta_ref,
                                              n_stages)
    print(f"  {'n':>3} {'a_n (us)':>15} {'b_n':>20}")
    for n in range(n_stages):
        a_us = float(a_wheeler[n]) * 1e6
        b_val = float(b_wheeler[n])
        print(f"  {n:>3} {a_us:>15.6e} {b_val:>20.6e}")

    # Reconstruct alpha_n via Gauss quadrature on Jacobi matrix
    print(f"\n  Reconstructing alpha_n from (a, b) via Gauss quadrature...")
    alphas_reconstr = alpha_from_jacobi(a_wheeler, b_wheeler, n_moments)
    print(f"  alpha_0 reconstructed = {float(alphas_reconstr[0]):.4e}")
    # Direct comparison
    alpha0_direct = sum(g2[k] for k in range(len(g2)))
    print(f"  alpha_0 direct        = {float(alpha0_direct):.4e}")
    rel = abs(alphas_reconstr[0] - alpha0_direct) / abs(alpha0_direct)
    print(f"  relative diff         = {float(rel):.2e}")

    # QD-Padé on reconstructed alpha
    print(f"\n  QD-Padé on reconstructed alpha (via Wheeler)...")
    p_wheeler = qd_pade_cauer(alphas_reconstr, n_stages)
    print(f"  {'k':>3} {'R_2k':>15} {'L_2k+1':>15} {'tau_pair us':>15}")
    for k in range(min(7, len(p_wheeler) // 2)):
        Rv = p_wheeler[2 * k]
        Linv = p_wheeler[2 * k + 1]
        if abs(Linv) < mp.mpf(10) ** (-50):
            break
        Lv = 1 / Linv
        tau_p = float(Lv / Rv) * 1e6
        print(f"  {k:>3} {float(Rv):>15.4e} {float(Lv):>15.4e} {tau_p:>15.4f}")

    # Direct QD-Padé reference
    print(f"\n  Direct QD-Padé on alpha = sum g^2 tau^n (reference)...")
    alphas_direct = []
    for n in range(n_moments):
        a_val = mp.mpf(0)
        for gv, tv in zip(g2, tau):
            a_val = a_val + gv * (tv ** n)
        alphas_direct.append((mp.mpf(-1)) ** n * a_val)
    p_direct = qd_pade_cauer(alphas_direct, n_stages)
    print(f"  {'k':>3} {'R_2k':>15} {'L_2k+1':>15} {'tau_pair us':>15}")
    for k in range(min(7, len(p_direct) // 2)):
        Rv = p_direct[2 * k]
        Linv = p_direct[2 * k + 1]
        if abs(Linv) < mp.mpf(10) ** (-50):
            break
        Lv = 1 / Linv
        tau_p = float(Lv / Rv) * 1e6
        print(f"  {k:>3} {float(Rv):>15.4e} {float(Lv):>15.4e} {tau_p:>15.4f}")


if __name__ == "__main__":
    main()
