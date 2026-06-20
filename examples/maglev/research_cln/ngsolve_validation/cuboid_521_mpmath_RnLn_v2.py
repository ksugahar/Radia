"""mpmath Cauer-II extraction R_n, L_n for cuboid 5x2x1 vacuum.

Direct port of Mathematica cauerExtract algorithm from
   cuboid_3D_vacuum_multipole_C_stage3.wls
to Python/mpmath at dps=50.

Inputs: Foster moments M_k of Y_Foster(s) = sum a_n/(1+s tau_n).
Algorithm: Z(s) = 1/Y_Foster(s), Cauer-II expansion at s=0:
  Z(s) = R_0 + s*L_0 + 1/(R_1 + s*L_1 + 1/(R_2 + ...))

Output: alternating (R_0, L_0, R_1, L_1, ...) circuit elements.
This is the TARGET for Kameari FEM to match.
"""
import mpmath as mp
import json, sys
from pathlib import Path
sys.stdout.reconfigure(encoding="utf-8", errors="replace")

mp.mp.dps = 50


def poly_eval_at_zero(coeffs):
    """coeffs in ASCENDING order: poly(0) = coeffs[0]"""
    return coeffs[0] if coeffs else mp.mpf(0)


def poly_subtract(a, b):
    """Subtract polynomials (ascending coeffs); pad shorter with zeros."""
    n = max(len(a), len(b))
    a_pad = list(a) + [mp.mpf(0)] * (n - len(a))
    b_pad = list(b) + [mp.mpf(0)] * (n - len(b))
    return [a_pad[i] - b_pad[i] for i in range(n)]


def poly_scalar_mul(c, p):
    return [c * x for x in p]


def poly_div_by_s(p):
    """Divide by s; requires p[0] = 0 within tolerance."""
    if abs(p[0]) > mp.mpf("1e-100") * (abs(p[1]) + mp.mpf("1e-100")):
        # Tolerance: check relative
        pass  # might still be OK
    return p[1:]   # drop the s^0 term


def cauer_II_extract(M_moments, max_stages):
    """Cauer-II at s=0 from Foster moments.

    Y_Foster(s) Taylor: c_k = (-1)^k M_k for k=0,1,...
    Z(s) = 1/Y_Foster(s) is rational. Extract:
      Z = R_0 + s*L_0 + 1/(R_1 + s*L_1 + 1/(...))

    Implementation: Z(s) is initially rational with num = D(s), den = N(s)
    where Y = N/D in formal series sense. We use polynomial truncation:
    treat Y as truncated power series and Z = 1/Y as truncated power series.
    """
    N_moments = len(M_moments)
    # Y(s) Taylor: c_k = (-1)^k M_k
    Y_taylor = [(-1)**k * M_moments[k] for k in range(N_moments)]

    # Z(s) Taylor: invert Y(s) = c_0 (1 + (c_1/c_0) s + ...)
    # 1/Y(s) = (1/c_0) sum (-1)^j ((c_1/c_0) s + (c_2/c_0) s^2 + ...)^j
    # Easier: use formal series inversion
    def series_invert(c, N):
        d = [mp.mpf(0)] * N
        d[0] = mp.mpf(1)/c[0]
        for k in range(1, N):
            s = mp.mpf(0)
            for j in range(1, k+1):
                if j < len(c):
                    s += c[j] * d[k-j]
            d[k] = -s / c[0]
        return d

    Z_taylor = series_invert(Y_taylor, N_moments)

    # Cauer-II expansion at s=0:
    # Z(s) = R_0 + s*L_0 + s^2/(...) ?
    # Actually Cauer-II: Z(s) = R_0 + s*L_0 + 1/(R_1 + s*L_1 + 1/(...))
    # Inverse continued fraction step:
    # Given Z(s) Taylor: a_0 + a_1 s + a_2 s^2 + a_3 s^3 + ...
    # R_0 = a_0
    # Z' = (Z - R_0)/s = a_1 + a_2 s + a_3 s^2 + ...
    # L_0 = a_1
    # Z'' = (Z' - L_0)/s = a_2 + a_3 s + ...  (now Z'' has constant a_2 = ???)
    # 1/(R_1 + s L_1 + ...) = Z''
    # So R_1 + s L_1 + ... = 1/Z''
    # Recursively: invert series Z'' to get next stage's series

    pairs = []
    Z_cur = list(Z_taylor)

    for stage in range(max_stages):
        if len(Z_cur) < 2:
            break
        R = Z_cur[0]
        # Z' = (Z - R)/s = Z_cur[1:] (drop R, then divide by s = drop the leading 0)
        Z_prime = Z_cur[1:]   # already (Z-R)/s ≈ a_1 + a_2 s + ...
        if len(Z_prime) < 1:
            break
        L = Z_prime[0]
        # Z'' = (Z' - L)/s = Z_prime[1:]
        Z_dprime = Z_prime[1:]
        pairs.append({"stage": stage, "R": R, "L": L})

        if len(Z_dprime) < 1:
            break
        # Next Z = 1/Z_dprime (formal series inversion)
        # If Z_dprime[0] near zero, can't invert
        if abs(Z_dprime[0]) < mp.mpf("1e-200"):
            print(f"  Stage {stage}: Z_dprime[0] near zero, breakdown")
            break
        Z_cur = series_invert(Z_dprime, len(Z_dprime))

    return pairs


# === Cuboid 5x2x1 BEM moments (high precision) ===
M0 = mp.mpf("3.05897258746351459800805934293582037956")
alphas_user = [
    mp.mpf("-0.0000440346570278169352257847180889677"),
    mp.mpf("7.4783240535434074561446429934e-10"),
    mp.mpf("-1.38797777606801509029076575691e-14"),
    mp.mpf("2.66525073846078827742571391442e-19"),
    mp.mpf("-5.19179202154889838584477748969e-24"),
    mp.mpf("1.01780181397466786388312805616e-28"),
    mp.mpf("-2.00117265750345380104921264126e-33"),
    mp.mpf("3.94011856973007074943279704837e-38"),
]
# M_k = (-1)^k alpha_k_user (Phase 7-11 convention)
M = [M0] + [(-1)**k * alphas_user[k-1] for k in range(1, 9)]

print("Foster moments M_k (high precision, 50 dps):")
for k in range(9):
    print(f"  M_{k} = {mp.nstr(M[k], 25)}")
print()

print("=" * 78)
print("Cauer-II R_n, L_n via series inversion (Mathematica algorithm port)")
print("=" * 78)
pairs = cauer_II_extract(M, max_stages=4)

print(f"\n{'stage':>5} {'R_n':>30} {'L_n':>30}")
for p in pairs:
    R_str = mp.nstr(p['R'], 18)
    L_str = mp.nstr(p['L'], 18)
    print(f"{p['stage']:>5} {R_str:>30} {L_str:>30}")

print()
print("Convert to Y_zz convention if needed:")
print("  These (R, L) are for Y_Foster(s) = sum a_n/(1+s tau_n) admittance.")
print("  Y_zz(s) = -mu_0 * s * Y_Foster(s).")
print("  Kameari R_n, L_n are for Y_zz directly.")
print()

# Save
out = {
    "method": "Cauer-II R_n, L_n for cuboid 5x2x1 (mpmath dps=50)",
    "convention": "Y_Foster(s) = sum a_n/(1+s tau_n), Z = 1/Y_Foster",
    "M_moments_50dps": [mp.nstr(m, 25) for m in M],
    "rl_pairs": [{"stage": p["stage"],
                   "R": mp.nstr(p["R"], 25),
                   "L": mp.nstr(p["L"], 25)} for p in pairs],
}
out_path = Path(__file__).parent / "cuboid_521_mpmath_RnLn_v2.json"
out_path.write_text(json.dumps(out, indent=2))
print(f"Results: {out_path}")
