"""2D rect Case B (uniform H_0 applied) via DUALITY from Case A.

Skin function g(s) shared between Case A and Case B:
  Y_A(s) = M_0^A * g(s)            (current-driven port, M_0^A = sigma a b)
  Y_B(s) = (1 - g(s)) / s            (magnetic admittance port)

Moment relation (simple shift):
  M_k^B = M_{k+1}^A / M_0^A

Procedure:
  1. Read NGSolve Case A Cauer ladder (R_n^A, L_n^A) from JSON.
  2. Reconstruct Z_A(s) as rational poly num/den by walking the ladder
     from innermost outward. Then Y_A(s) = 1/Z_A(s) Taylor -> M_k^A.
  3. Apply shift M_k^B = M_{k+1}^A / M_0^A.
  4. Build Y_B(s) = sum (-s)^k M_k^B and extract Cauer-I via polynomial
     division.

Validation: compare to rect_2D_caseB_moments.wls (Mathematica reference).
"""
import json
from pathlib import Path
import sys

from mpmath import mp, mpf

mp.dps = 50


def poly_add(p, q):
    if len(p) < len(q):
        p, q = q, p
    out = list(p)
    for i, c in enumerate(q):
        out[i] = out[i] + c
    return out


def poly_mul(p, q):
    if not p or not q:
        return [mpf(0)]
    out = [mpf(0)] * (len(p) + len(q) - 1)
    for i, ci in enumerate(p):
        for j, cj in enumerate(q):
            out[i+j] = out[i+j] + ci * cj
    return out


def poly_scale(p, s):
    return [s * c for c in p]


def cauer_to_Y_taylor(R_list, L_list, K_max):
    """Reconstruct Y(s) Taylor coefficients (sum y_k s^k) from Cauer-I.

    Cauer-I structure (R-first, RL ladder):
      Z_0 = R_0 + (sL_0 || Z_1)
      Z_k = R_k + (sL_k || Z_{k+1})
      Z_{N-1} = R_{N-1} + sL_{N-1}    (terminate with parallel Z_N=infinity)

    Returns list y_0, y_1, ..., y_{K_max} where Y(s) = sum y_k s^k.
    """
    N = min(len(R_list), len(L_list))
    if N == 0:
        return None

    # Innermost: Z_{N-1} = R_{N-1} + s L_{N-1}
    num = [mpf(R_list[N-1]), mpf(L_list[N-1])]
    den = [mpf(1)]

    # Walk outward: k = N-2 down to 0
    for k in range(N-2, -1, -1):
        R_k = mpf(R_list[k])
        L_k = mpf(L_list[k])
        # Step 1: parallel Z_in || sL_k
        # Z_par_num/Z_par_den = (num × [0, L_k]) / (num + den × [0, L_k])
        sLk = [mpf(0), L_k]
        par_num = poly_mul(num, sLk)
        par_den = poly_add(num, poly_mul(den, sLk))
        # Step 2: series + R_k
        # Z_new_num/Z_new_den = (R_k × par_den + par_num) / par_den
        new_num = poly_add(poly_scale(par_den, R_k), par_num)
        new_den = par_den
        num, den = new_num, new_den

    # Y(s) = 1/Z(s) = den/num
    Y_num = den
    Y_den = num

    # Compute Taylor: Y(s) = Y_num/Y_den; y_0 = Y_num[0]/Y_den[0]
    # y_k = (Y_num[k] - sum_{j=1..k} y_{k-j} Y_den[j]) / Y_den[0]
    Yn = list(Y_num) + [mpf(0)] * (K_max + 1)
    Yd = list(Y_den) + [mpf(0)] * (K_max + 1)
    y = [mpf(0)] * (K_max + 1)
    y[0] = Yn[0] / Yd[0]
    for k in range(1, K_max + 1):
        s = Yn[k] if k < len(Y_num) else mpf(0)
        for j in range(1, k+1):
            if j < len(Yd):
                s = s - y[k-j] * Yd[j]
        y[k] = s / Yd[0]
    return y


def Y_taylor_to_moments(y_list):
    """Y(s) = sum y_k s^k = sum (-s)^k M_k.  So M_k = (-1)^k y_k."""
    return [(-1)**k * y_list[k] for k in range(len(y_list))]


def moments_to_cauer(M_list, n_stages):
    """Extract Cauer-I from moment series Y(s) = sum (-s)^k M_k.

    Same algorithm as Mathematica wls cauerExtract.
    """
    K = len(M_list)
    # Y(s) polynomial in s: y_k = (-1)^k M_k
    Ypoly = [(-1)**k * M_list[k] for k in range(K)]
    # Z(s) = 1/Y(s). For Cauer-I extraction, treat as rational num/den.
    # Initially: znum = [1], zden = Ypoly
    znum = [mpf(1)]
    zden = list(Ypoly)

    R_list = []
    L_list = []
    EPS = mpf(10)**(-mp.dps + 10)

    for kk in range(n_stages):
        while len(znum) > 1 and abs(znum[-1]) < EPS:
            znum.pop()
        while len(zden) > 1 and abs(zden[-1]) < EPS:
            zden.pop()
        if not zden or abs(zden[0]) < EPS:
            break
        R_k = znum[0] / zden[0]
        R_list.append(float(R_k))
        # presid = znum - R_k zden
        Lpresid = max(len(znum), len(zden))
        presid = []
        for i in range(Lpresid):
            zn = znum[i] if i < len(znum) else mpf(0)
            zd = zden[i] if i < len(zden) else mpf(0)
            presid.append(zn - R_k * zd)
        # pas = presid / s (drop first)
        if abs(presid[0]) > EPS:
            pass  # should be ~0
        pas = presid[1:]
        if not pas or abs(pas[0]) < EPS:
            break
        L_k = pas[0] / zden[0]
        L_list.append(float(L_k))
        # nresid = zden L_k - pas
        Lnresid = max(len(zden), len(pas))
        nresid = []
        for i in range(Lnresid):
            zd = zden[i] if i < len(zden) else mpf(0)
            ps = pas[i] if i < len(pas) else mpf(0)
            nresid.append(zd * L_k - ps)
        nresid_shifted = nresid[1:]
        # New: znum = L_k pas, zden = nresid_shifted
        znum_new = [L_k * c for c in pas]
        zden_new = nresid_shifted
        znum, zden = znum_new, zden_new

    return R_list, L_list


def main():
    here = Path(__file__).parent
    ngs_A = json.loads((here / "rect_2D_kameari_results.json").read_text())
    R_A = ngs_A["R"]
    L_A = ngs_A["L"]   # in H/m

    print(f"Source: rect_2D_kameari.py output, h={ngs_A['h_mm']}mm, order={ngs_A['order']}")
    print(f"  Case A R_0 = {R_A[0]:.10e} Ohm/m, L_0 = {L_A[0]*1e9:.10e} nH/m")
    print(f"  Expected: M_0^A = 1/R_0 = {1/R_A[0]:.10e} = sigma*a*b")

    # Reconstruct Y_A(s) Taylor and moments
    K_recon = 14
    y_A = cauer_to_Y_taylor(R_A, L_A, K_recon)
    if y_A is None:
        print("Failed")
        return
    M_A = Y_taylor_to_moments(y_A)
    print(f"\nReconstructed M_k^A:")
    for k in range(min(6, len(M_A))):
        print(f"  M_{k}^A = {float(M_A[k]):.10e}")

    # Sanity: M_0^A should equal 1/R_0^A
    expected_M0 = 1.0 / R_A[0]
    actual_M0 = float(M_A[0])
    print(f"  expected M_0^A = 1/R_0 = {expected_M0:.10e}")
    print(f"  actual = {actual_M0:.10e}")
    if abs(actual_M0 - expected_M0)/expected_M0 > 1e-6:
        print(f"  WARNING: discrepancy {abs(actual_M0 - expected_M0)/expected_M0:.2e}")
    else:
        print(f"  OK (within 1e-6 relative)")

    # Apply duality: M_k^B = M_{k+1}^A / M_0^A
    M0A = M_A[0]
    M_B = [M_A[k+1] / M0A for k in range(min(len(M_A) - 1, 12))]
    print(f"\nCase B moments via shift:")
    for k in range(min(6, len(M_B))):
        print(f"  M_{k}^B = {float(M_B[k]):.10e}")

    # Cauer-I from M_B
    R_B, L_B = moments_to_cauer(M_B, 12)
    print(f"\nCase B Cauer-I from NGSolve duality ({len(R_B)} stages):")
    for n in range(min(len(R_B), len(L_B))):
        print(f"  n={n}: R_B={R_B[n]:.10e}, L_B={L_B[n]:.10e}")

    out = {
        "method": "Case B Cauer-I via duality from Case A NGSolve, M_k^B = M_{k+1}^A / M_0^A",
        "source": "rect_2D_kameari_results.json",
        "h_mm": ngs_A["h_mm"],
        "order": ngs_A["order"],
        "ne": ngs_A["ne"],
        "M_A_reconstructed": [float(m) for m in M_A[:8]],
        "M_B_shifted": [float(m) for m in M_B[:8]],
        "R_B": R_B,
        "L_B": L_B,
    }
    out_path = here / "rect_2D_caseB_via_duality_results.json"
    out_path.write_text(json.dumps(out, indent=2))
    print(f"\nResults: {out_path}")


if __name__ == "__main__":
    main()
