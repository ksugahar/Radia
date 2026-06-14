"""Python pure Foster-CLN library (no Mathematica dependency).

Full pipeline at ~17 digit precision via mpmath:
  1. Cell-cell Newton K matrix via mp.quad (default Tanh-Sinh)
  2. Legendre basis transformation
  3. Symmetric eigenvalue (mp.eighe)
  4. Foster amplitudes
  5. alpha_n moments
  6. CLN-I (s=0) continued fraction
  7. CLN ladder

Same precision as MATLAB Symbolic Toolbox + vpaintegral.
For 30 digit, use Mathematica path (bem_foster_spectral_v5.wls).

Date: 2026-05-02
"""

import time
import json

import mpmath as mp
import sympy as sp

mp.mp.dps = 50
SDPS = 50

sigma = mp.mpf("5.8e7")
mu0 = 4 * mp.pi * mp.mpf("1e-7")
a = mp.mpf("5e-3")
b = mp.mpf("2e-3")
c = mp.mpf("1e-3")
V = a * b * c
mu0sigma = mu0 * sigma

# References
alpha1_ref = -mu0sigma * (a*a + b*b) / 48
alpha2_ref_new = mp.mpf("7.4783240535434074691212813799e-10")
tau1_ref = alpha2_ref_new / (-alpha1_ref)

P_MAX = 4
N_MOMENTS = 8

print("=" * 64)
print(f"Python pure Foster-CLN library (P={P_MAX} dps={SDPS})")
print("=" * 64)
print()

# === gAxis polynomial coefficients (SymPy) ===
_eta_s, _zeta_s, _xi_s, _u_s = sp.symbols('eta zeta xi u', positive=True)
_gAxis_cache = {}

def gAxis_xi_coeffs(p1, p2, h_mp):
    if (p1, p2) not in _gAxis_cache:
        expr = sp.integrate(
            (_u_s + _xi_s/2)**p1 * (_u_s - _xi_s/2)**p2,
            (_u_s, _xi_s/2 - sp.Symbol('h')/2, sp.Symbol('h')/2 - _xi_s/2))
        expr = sp.expand(expr)
        poly_xi = sp.Poly(expr, _xi_s)
        _gAxis_cache[(p1, p2)] = poly_xi.all_coeffs()[::-1]
    sym_coeffs = _gAxis_cache[(p1, p2)]
    h_subs = sp.Float(mp.nstr(h_mp, SDPS), SDPS)
    h_sym = sp.Symbol('h')
    return [mp.mpf(str(c.subs(h_sym, h_subs).evalf(SDPS))) for c in sym_coeffs]


def I_xi(k, h, R):
    if k == 0: return mp.asinh(h / R)
    if k == 1: return mp.sqrt(h*h + R*R) - R
    Ikm2 = I_xi(k - 2, h, R)
    return (h**(k-1) * mp.sqrt(h*h + R*R) - (k-1) * R*R * Ikm2) / k


def IxiClosed(coeffs_xi, h, R):
    return mp.fsum(c * I_xi(k, h, R) for k, c in enumerate(coeffs_xi) if c != 0)


def poly_eval(coeffs, x):
    result = mp.mpf(0)
    for c in reversed(coeffs):
        result = result * x + c
    return result


_K_cache = {}
def K_mon(p1q1r1, p2q2r2):
    key = (p1q1r1, p2q2r2)
    if key in _K_cache:
        return _K_cache[key]
    p1, q1, r1 = p1q1r1
    p2, q2, r2 = p2q2r2
    cx = gAxis_xi_coeffs(p1, p2, a)
    cy = gAxis_xi_coeffs(q1, q2, b)
    cz = gAxis_xi_coeffs(r1, r2, c)

    def integrand(eta, zeta):
        R = mp.sqrt(eta*eta + zeta*zeta)
        ixi = IxiClosed(cx, a, R)
        return ixi * poly_eval(cy, eta) * poly_eval(cz, zeta)

    val = 8 * mp.quad(integrand, [0, b], [0, c])
    _K_cache[key] = val
    return val


def build_basis(P_max, parity):
    basis = []
    px, py, pz = parity
    for p in range(P_max + 1):
        for q in range(P_max - p + 1):
            for r in range(P_max - p - q + 1):
                if p % 2 == px and q % 2 == py and r % 2 == pz:
                    basis.append((p, q, r))
    return basis


def legendre_T(P_max, h):
    x = sp.Symbol('x')
    T = mp.matrix(P_max + 1, P_max + 1)
    h_sub = sp.Float(mp.nstr(h, SDPS), SDPS)
    for n in range(P_max + 1):
        Ln = sp.legendre(n, x)
        cofs = sp.Poly(Ln, x).all_coeffs()[::-1]
        for k, ck in enumerate(cofs):
            v = ck * (sp.Float(2, SDPS) / h_sub)**k
            T[n, k] = mp.mpf(str(v.evalf(SDPS)))
    return T


def solve_block(basis, axis_label):
    n = len(basis)
    print(f"  Block {axis_label}: {n} basis functions")
    print(f"    Building K_mon ({n*(n+1)//2} unique elements)...")
    t0 = time.time()
    K = mp.matrix(n, n)
    for i in range(n):
        for j in range(i, n):
            val = K_mon(basis[i], basis[j])
            K[i, j] = val
            K[j, i] = val
    print(f"      K_mon time: {time.time()-t0:.1f} s")

    Tx = legendre_T(P_MAX, a)
    Ty = legendre_T(P_MAX, b)
    Tz = legendre_T(P_MAX, c)
    Tblock = mp.matrix(n, n)
    for i in range(n):
        P, Q, R = basis[i]
        for j in range(n):
            p, q, r = basis[j]
            Tblock[i, j] = Tx[P, p] * Ty[Q, q] * Tz[R, r]

    K_leg = Tblock * K * Tblock.T
    mass = [a / (2*P + 1) * b / (2*Q + 1) * c / (2*R + 1) for (P, Q, R) in basis]
    K_std = mp.matrix(n, n)
    for i in range(n):
        for j in range(n):
            K_std[i, j] = K_leg[i, j] / mp.sqrt(mass[i] * mass[j])

    eigvals, eigvecs = mp.eighe(K_std)
    return basis, eigvals, eigvecs, mass


def project_aext(basis, mass, axis):
    n = len(basis)
    vec = mp.matrix(n, 1)
    if axis == "x":
        target = (0, 1, 0); coef = -b/4
    else:
        target = (1, 0, 0); coef = a/4
    if target in basis:
        idx = basis.index(target)
        vec[idx, 0] = coef * mp.sqrt(mass[idx])
    return vec


def invert_taylor(c):
    n = len(c)
    e = [mp.mpf(0)] * n
    e[0] = 1 / c[0]
    for k in range(1, n):
        s = mp.fsum(c[j] * e[k - j] for j in range(1, k + 1))
        e[k] = -s / c[0]
    return e


def cauer_extract(c_taylor, max_stages):
    c = list(c_taylor)
    p_list = []
    for _ in range(max_stages):
        if len(c) < 2 or abs(c[0]) < mp.mpf("1e-300"):
            break
        e = invert_taylor(c)
        p_list.append(e[0])
        c = e[1:]
    return p_list


def cauer_to_poles(p_list, K):
    if K < 1: return None, None
    N = [mp.mpf(1)]
    D = [p_list[K - 1]]
    for k in range(K - 2, -1, -1):
        sN = [mp.mpf(0)] + N
        Dnew_pk = [p_list[k] * d for d in D]
        max_len = max(len(Dnew_pk), len(sN))
        Dnew = [mp.mpf(0)] * max_len
        for i in range(len(Dnew_pk)): Dnew[i] += Dnew_pk[i]
        for i in range(len(sN)):       Dnew[i] += sN[i]
        N, D = D, Dnew
    return N, D


# === Main pipeline ===
print("=== Solving Foster eigenvalue problems ===")
bX, eX, vX, mX = solve_block(build_basis(P_MAX, (0, 1, 0)), "x-coupling (e,o,e)")
bY, eY, vY, mY = solve_block(build_basis(P_MAX, (1, 0, 0)), "y-coupling (o,e,e)")

aExtX = project_aext(bX, mX, "x")
aExtY = project_aext(bY, mY, "y")

# Foster amplitudes (eigvecs columns are eigenvectors)
gX = [mp.fsum(vX[i, k] * aExtX[i, 0] for i in range(len(bX))) for k in range(len(bX))]
gY = [mp.fsum(vY[i, k] * aExtY[i, 0] for i in range(len(bY))) for k in range(len(bY))]

lambdaX = [eX[k] / (4 * mp.pi) for k in range(len(bX))]
lambdaY = [eY[k] / (4 * mp.pi) for k in range(len(bY))]

# alpha_n
alphas = []
for nM in range(1, N_MOMENTS + 1):
    sX = mp.fsum(gX[k]**2 * lambdaX[k]**(nM - 1) for k in range(len(bX)))
    sY = mp.fsum(gY[k]**2 * lambdaY[k]**(nM - 1) for k in range(len(bY)))
    alphas.append((-1)**nM * mu0sigma**nM / V * (sX + sY))

# alpha_0
alpha_0 = mp.fsum(gX[k]**2 / (V * lambdaX[k]) for k in range(len(bX))) + \
          mp.fsum(gY[k]**2 / (V * lambdaY[k]) for k in range(len(bY)))

print()
print("=== alpha_n vs reference ===")
print(f"  alpha_1 = {mp.nstr(alphas[0], 18)}")
print(f"    err = {mp.nstr(abs((alphas[0] - alpha1_ref)/alpha1_ref), 6)}")
print(f"  alpha_2 = {mp.nstr(alphas[1], 18)}")
print(f"    err = {mp.nstr(abs((alphas[1] - alpha2_ref_new)/alpha2_ref_new), 6)}")
for n in range(3, N_MOMENTS + 1):
    print(f"  alpha_{n} = {mp.nstr(alphas[n-1], 18)}")
tau1_BEM = alphas[1] / (-alphas[0])
print(f"  tau_1   = {mp.nstr(tau1_BEM*1e6, 18)} us")
print(f"  tau_1 err = {mp.nstr(abs((tau1_BEM - tau1_ref)/tau1_ref), 6)}")

# === Cauer extraction ===
print()
print("=== CLN-I (s=0) continued fraction ===")
c_taylor = [alpha_0] + alphas
p_cauer = cauer_extract(c_taylor, len(c_taylor) - 1)
print(f"  Extracted {len(p_cauer)} stages")
for k, p in enumerate(p_cauer):
    print(f"    p_{k+1} = {mp.nstr(p, 14)}")

# Cauer truncation -> poles
print()
print("=== Cauer truncation -> tau_k [us] ===")
print(f"  K  poles  tau_k [us]")
for K in range(2, len(p_cauer) + 1):
    Nc, Dc = cauer_to_poles(p_cauer, K)
    coefs_hi_to_lo = list(reversed(Dc))
    try:
        roots = mp.polyroots(coefs_hi_to_lo, maxsteps=200, extraprec=20)
        real_neg = [r.real if hasattr(r, 'real') else r
                    for r in roots
                    if (hasattr(r, 'real') and r.real < 0 and abs(r.imag) < mp.mpf("1e-15") * abs(r.real))
                    or (not hasattr(r, 'imag') and r < 0)]
        if real_neg:
            taus_us = sorted([-1/p * 1e6 for p in real_neg], reverse=True)
            tau_str = ", ".join(mp.nstr(t, 10) for t in taus_us)
            print(f"  {K}  {len(real_neg)}  {tau_str}")
    except Exception as e:
        print(f"  {K}  failed: {e}")

# === Save ===
results = {
    "P_max": P_MAX,
    "dps": SDPS,
    "alpha_n": [mp.nstr(a, 30) for a in alphas],
    "alpha_0": mp.nstr(alpha_0, 30),
    "tau_1_us": mp.nstr(tau1_BEM*1e6, 30),
    "tau_1_err": mp.nstr(abs((tau1_BEM - tau1_ref)/tau1_ref), 6),
    "p_cauer": [mp.nstr(p, 18) for p in p_cauer],
}
with open("foster_cln_python_pure_results.json", "w") as f:
    json.dump(results, f, indent=2)
print("\nResults saved to foster_cln_python_pure_results.json")
