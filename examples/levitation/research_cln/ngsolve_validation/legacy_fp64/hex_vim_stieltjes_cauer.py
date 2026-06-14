"""hex_vim_stieltjes_cauer.py — VIM-native Cauer extraction via Stieltjes
3-term recurrence on the Foster spectrum.

Mathematically equivalent to Lanczos on (K, M, b), but operates directly on
the Foster spectrum {tau_k, g_k^2} from the VIM eigenproblem. AVOIDS the
moment alpha_n = (-1)^n sum g_k^2 tau_k^n computation, which suffers from
alternating-sign cancellation and dynamic range issues at high n.

Algorithm (Stieltjes / Gragg / modified Chebyshev):
  Given {tau_k, g_k^2}, k = 1..N (Foster spectrum, all positive weights)
  Initialize p_0 = 1, p_{-1} = 0 (as functions on the discrete grid)
  For n = 0, 1, 2, ...:
    a_n = sum_k g_k^2 tau_k p_n(tau_k)^2 / sum_k g_k^2 p_n(tau_k)^2
    b_n = sum_k g_k^2 p_n(tau_k)^2 / sum_k g_k^2 p_{n-1}(tau_k)^2  (b_0 = 0)
    p_{n+1}(tau) = (tau - a_n) p_n(tau) - b_n p_{n-1}(tau)

  The (a_n, b_n) coefficients give the Jacobi matrix
    J = [[a_0,   sqrt(b_1), 0,     ...],
         [sqrt(b_1), a_1, sqrt(b_2), ...],
         ...]
  whose eigenvalues are tau_k (Gauss quadrature on the discrete measure).

  Cauer-I rungs in KAMEARI convention:
    R_2k     ←  ?
    L_{2k+1} ←  ?
  The mapping depends on which generating function we Padé-fit. For our
  problem, χ(s) = sum g_k^2 / (1 + s tau_k) is the generating function;
  the Cauer-I continued fraction of χ(s) gives:
    1/χ(s) = R_0 + 1/(s L_1 + 1/(R_2 + 1/(s L_3 + ...)))
  whose coefficients are related to (a_n, b_n) by the Stieltjes-Cauer
  mapping. See Akhiezer "The Classical Moment Problem", or
  Kim-Heinz "Recursive computation of Cauer ladder networks".

Reference test: extract from cuboid 5x2x1 / A1 Foster JSONs (radia-vim
GPU output). Compare against QD-Padé reference.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")

import mpmath as mp
mp.mp.dps = 80


def stieltjes_recurrence(tau, g2, n_stages):
    """Compute 3-term recurrence coefficients (a_n, b_n) for orthogonal
    polynomials on the discrete measure dmu(tau) = sum g_k^2 delta(tau-tau_k).

    Args:
      tau: list of mp.mpf, Foster time constants (all > 0)
      g2:  list of mp.mpf, Foster amplitudes squared (all > 0)
      n_stages: number of recurrence steps

    Returns:
      a: list of length n_stages, diagonal of Jacobi matrix
      b: list of length n_stages, off-diagonal (b[0] = 0, b[n] = ||p_n||^2 / ||p_{n-1}||^2)
      p_norm2: list of <p_n, p_n> values (for diagnostics)
    """
    N = len(tau)
    assert len(g2) == N

    # Initialize: p_0(tau_k) = 1 for all k, p_{-1}(tau_k) = 0
    p_prev = [mp.mpf(0)] * N         # p_{-1}
    p_curr = [mp.mpf(1)] * N         # p_0

    a = []
    b = []
    p_norm2 = []

    # <p_{-1}, p_{-1}> set to 1 by convention (won't be used since b_0 = 0)
    norm2_prev = mp.mpf(1)
    norm2_curr = sum(g2[k] * p_curr[k] ** 2 for k in range(N))  # = sum g_k^2
    p_norm2.append(norm2_curr)

    for n in range(n_stages):
        # a_n = <tau p_n, p_n> / <p_n, p_n>
        num = sum(g2[k] * tau[k] * p_curr[k] ** 2 for k in range(N))
        a_n = num / norm2_curr
        a.append(a_n)

        # b_n = <p_n, p_n> / <p_{n-1}, p_{n-1}>  (b_0 = 0 by convention)
        if n == 0:
            b_n = mp.mpf(0)
        else:
            b_n = norm2_curr / norm2_prev
        b.append(b_n)

        # p_{n+1}(tau) = (tau - a_n) p_n - b_n p_{n-1}
        p_next = [(tau[k] - a_n) * p_curr[k] - b_n * p_prev[k]
                  for k in range(N)]

        norm2_next = sum(g2[k] * p_next[k] ** 2 for k in range(N))
        p_norm2.append(norm2_next)

        # Shift
        p_prev = p_curr
        p_curr = p_next
        norm2_prev = norm2_curr
        norm2_curr = norm2_next

    return a, b, p_norm2


def stieltjes_to_cauer_kameari(a, b, p_norm2):
    """Convert Stieltjes recurrence (a_n, b_n) to Kameari Cauer-I rungs.

    For chi(s) = sum_k g_k^2 / (1 + s tau_k), the Cauer-I expansion is
        chi(s) = c_0 / (1 + s a_0 - b_1 * c_2(s) / (s - ...))
    Eqivalently, via the Hankel-Padé approach:
        chi(s) -> 1/chi(s) is approximated by truncated continued fraction
        R_0 + s L_1 + ...  ish

    The mapping from (a_n, b_n) to (R_{2n}, L_{2n+1}):
      R_0      = 1 / chi(0)              = 1 / sum g_k^2
                                          = 1 / p_norm2[0]
      Then the continued fraction coefficients of 1/chi(s) involve
      the orthogonal polynomial recurrence.

    For DEMONSTRATIVE purposes, we report (a, b) directly and observe
    convergence behavior; rigorous Stieltjes-Cauer mapping requires
    additional theoretical setup (see Akhiezer Ch. 3).

    Quick check: extract via Padé continued fraction of f(s) = chi(s)
    using the (a, b) coefficients to construct the polynomial recurrence
    in the variable z = 1/s, which gives Cauer ladder directly.
    """
    # Use Jacobi matrix eigenvalue ratios as proxy for stability
    # Build tridiagonal Jacobi matrix J = diag(a) + offdiag(sqrt(b[1:]))
    # Its largest eigenvalue ≈ tau_max
    # For verification, just report a, b
    n_stages = len(a)
    out = []
    for n in range(n_stages):
        out.append({
            "stage": n,
            "a_n": float(a[n]),
            "b_n": float(b[n]) if n > 0 else 0.0,
            "p_norm2": float(p_norm2[n + 1]),
        })
    return out


def hankel_pade_qd_on_alpha(tau, g2, n_moments, n_stages):
    """Reference QD-Padé: compute alpha_n moments then QD continued fraction.
    For comparison with Stieltjes approach.
    """
    n_use = len(tau)
    alphas = []
    for n in range(n_moments):
        a = mp.mpf(0)
        for gv, tv in zip(g2, tau):
            a = a + gv * (tv ** n)
        alphas.append((mp.mpf(-1)) ** n * a)

    # QD on alpha
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
        json_path = (Path("S:/Radia/01_GitHub/packages/radia-vim/scripts")
                     / "cuboid521_order4_results.json")
    else:
        json_path = Path(sys.argv[1])
    print(f"=== Stieltjes vs QD-Padé Cauer comparison ===")
    print(f"  Foster JSON: {json_path}")

    data = json.loads(json_path.read_text(encoding="utf-8"))
    modes = data["top_modes_by_g2"]
    n_use = len([m for m in modes if abs(m["g2"]) > 1e-32])
    print(f"  Using {n_use} Foster modes")

    tau = [mp.mpf(m["tau_us"]) * mp.mpf("1e-6") for m in modes[:n_use]]
    g2 = [mp.mpf(m["g2"]) for m in modes[:n_use]]

    n_stages = 12

    # Stieltjes 3-term recurrence (VIM-native)
    print(f"\n--- Stieltjes 3-term on Foster spectrum ---")
    a, b, p_norm2 = stieltjes_recurrence(tau, g2, n_stages)
    print(f"  {'n':>3} {'a_n (us)':>15} {'b_n':>15} {'||p_n||^2':>15}")
    for n in range(n_stages):
        a_us = float(a[n]) * 1e6
        b_val = float(b[n]) if n > 0 else 0.0
        norm2 = float(p_norm2[n + 1])
        print(f"  {n:>3} {a_us:>15.6e} {b_val:>15.6e} {norm2:>15.4e}")

    # QD-Padé reference (alpha-based)
    print(f"\n--- QD-Padé on alpha_n (reference) ---")
    p_cauer = hankel_pade_qd_on_alpha(tau, g2, n_moments=40, n_stages=n_stages)
    print(f"  {'k':>3} {'R_2k':>15} {'L_2k+1':>15} {'tau_pair us':>15}")
    for k in range(min(6, len(p_cauer) // 2)):
        Rv = p_cauer[2 * k]
        Linv = p_cauer[2 * k + 1]
        if abs(Linv) < mp.mpf(10) ** (-50):
            break
        Lv = 1 / Linv
        tau_p = float(Lv / Rv) * 1e6
        print(f"  {k:>3} {float(Rv):>15.4e} {float(Lv):>15.4e} {tau_p:>15.4f}")

    # Observation: a_n should be in time-constant units, around the
    # "typical" tau. b_n should be of order tau^2.
    print(f"\n--- Observations ---")
    print(f"  Foster top tau = {float(tau[0]) * 1e6:.4f} us")
    print(f"  Stieltjes a_0  = {float(a[0]) * 1e6:.4f} us "
          f"(weighted-mean tau)")
    print(f"  ||p_0||^2 = sum g^2 = {float(p_norm2[0]):.4e}")
    print(f"  R_0_canonical = 1 / sum g^2 = {1.0 / float(p_norm2[0]):.4e}")


if __name__ == "__main__":
    main()
