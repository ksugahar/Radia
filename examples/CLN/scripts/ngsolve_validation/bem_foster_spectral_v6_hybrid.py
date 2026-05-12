"""Spectral BEM v6 — Hybrid Mathematica + Python (mpmath dps=50).

PRODUCTION pipeline:
  1. Mathematica computes K matrix at WP=50 (30+ digit precision via GlobalAdaptive)
     -> bem_foster_K_export_P{P}.json
  2. Python loads K matrix (this script), runs:
     - Legendre basis transformation
     - Generalized eigenvalue problem (mpmath.eighe)
     - Foster amplitude projection
     - alpha_n moment calculation (n=1..8 at full precision)
     - tau_1 = alpha_2 / (-alpha_1)

Why hybrid:
  - Mathematica's NIntegrate handles corner log-singularity (arcsinh(h/R))
  - mpmath.quad limited to ~15 digits at this kind of integrand (documented)
  - Python is faster for matrix algebra and final analysis
  - Best of both worlds: 18+ digit alpha_n, fast pipeline

Date: 2026-05-02
"""

import time
import json

import mpmath as mp
import sympy as sp

mp.mp.dps = 50
SDPS = 50

# === Constants ===
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

import re

def parse_mma_number(s):
    """Convert Mathematica InputForm number to Python-readable scientific notation.

    Mathematica: '5.24918`30.*^-21' (backtick = precision marker, *^ = power)
                 '5.24918e-21 + 0`*i' (complex with zero imag part)
    Python:      '5.24918e-21' (just the real part)
    """
    # Strip imaginary part if zero. Mathematica formats:
    #   '+ 0`*I'  or  '+ 0``62.64*I'  (zero with arbitrary precision marker)
    # Match: ' + 0' followed by any number of `...digits, then optional *, then I/i
    s = re.sub(r'\s*[+-]\s*0`+[0-9.]*\s*\*?\s*[iI]\s*', '', s)
    s = re.sub(r'\s*[+-]\s*0\s*\*?\s*[iI]\s*', '', s)
    # Remove backtick + digits (precision marker, possibly with double backticks)
    s = re.sub(r'`+[0-9.]+\.?', '', s)
    # Replace *^ with e for exponent
    s = s.replace('*^', 'e')
    return s.strip()


P_MAX = 4
JSON_FILE = f"bem_foster_K_export_P{P_MAX}.json"

print("=" * 64)
print(f"Spectral BEM v6 (HYBRID Math+Python, P={P_MAX} dps={SDPS})")
print("=" * 64)
print()

# === Load K matrices from Mathematica JSON ===
print(f"Loading K matrices from {JSON_FILE}...")
t0 = time.time()
with open(JSON_FILE, "r") as f:
    data = json.load(f)
print(f"  WP from Mathematica: {data['WP']}")
print(f"  Loaded in {time.time()-t0:.2f} s")
print()


def parse_kmat(block_data):
    basis = [tuple(b) for b in block_data["basis"]]
    n = len(basis)
    K = mp.matrix(n, n)
    for entry in block_data["K_elements"]:
        i, j = entry["i"], entry["j"]
        K[i, j] = mp.mpf(parse_mma_number(entry["value"]))
    return basis, K


basisX, KMonX = parse_kmat(data["x_block"])
basisY, KMonY = parse_kmat(data["y_block"])
nX = len(basisX)
nY = len(basisY)
print(f"  x-block (e,o,e): {nX} basis functions")
print(f"  y-block (o,e,e): {nY} basis functions")
print()


# === Legendre transformation ===
def legendreT(P_max, h):
    """T[P, p] = coefficient of x^p in L_P(2x/h) for x in [-h/2, h/2]."""
    x = sp.Symbol('x')
    T = mp.matrix(P_max + 1, P_max + 1)
    h_sub = sp.Float(mp.nstr(h, SDPS), SDPS)
    for n in range(P_max + 1):
        Ln = sp.legendre(n, x)
        coeffs = sp.Poly(Ln, x).all_coeffs()[::-1]
        for k, ck in enumerate(coeffs):
            v = ck * (sp.Float(2, SDPS) / h_sub)**k
            T[n, k] = mp.mpf(str(v.evalf(SDPS)))
    return T


Tx = legendreT(P_MAX, a)
Ty = legendreT(P_MAX, b)
Tz = legendreT(P_MAX, c)


def legendre_block(basis, K_mon):
    n = len(basis)
    Tblock = mp.matrix(n, n)
    for i in range(n):
        P, Q, R = basis[i]
        for j in range(n):
            p, q, r = basis[j]
            Tblock[i, j] = Tx[P, p] * Ty[Q, q] * Tz[R, r]

    K_leg = Tblock * K_mon * Tblock.T
    massDiag = [a / (2*P + 1) * b / (2*Q + 1) * c / (2*R + 1)
                for (P, Q, R) in basis]

    K_leg_std = mp.matrix(n, n)
    for i in range(n):
        for j in range(n):
            K_leg_std[i, j] = K_leg[i, j] / mp.sqrt(massDiag[i] * massDiag[j])
    return K_leg_std, massDiag


# === Solve eigenvalue and compute alpha_n ===
def solve_block(basis, K_mon, axis):
    n = len(basis)
    K_leg_std, massDiag = legendre_block(basis, K_mon)

    print(f"  Solving symmetric eigenvalue ({axis}-block, n={n})...")
    t0 = time.time()
    eigvals, eigvecs = mp.eighe(K_leg_std)
    print(f"    eigvalue solve: {time.time()-t0:.4f} s")

    # Project A_ext onto M-orthonormal Legendre basis
    aExt = mp.matrix(n, 1)
    if axis == "x":
        target = (0, 1, 0)
        coef = -b / 4
    elif axis == "y":
        target = (1, 0, 0)
        coef = a / 4
    if target in basis:
        idx = basis.index(target)
        aExt[idx, 0] = coef * mp.sqrt(massDiag[idx])

    g = [mp.fsum(eigvecs[i, k] * aExt[i, 0] for i in range(n))
         for k in range(n)]
    lambdaK = [eigvals[k] / (4 * mp.pi) for k in range(n)]
    return g, lambdaK


print("=== Solving Foster eigenvalue problems ===")
gX, lambdaX = solve_block(basisX, KMonX, "x")
gY, lambdaY = solve_block(basisY, KMonY, "y")
print()

# === Compute alpha_n ===
NMoments = 8
alphas = []
for nM in range(1, NMoments + 1):
    sX = mp.fsum(gX[k]**2 * lambdaX[k]**(nM - 1) for k in range(nX))
    sY = mp.fsum(gY[k]**2 * lambdaY[k]**(nM - 1) for k in range(nY))
    a_n = (-1)**nM * mu0sigma**nM / V * (sX + sY)
    alphas.append(a_n)

print("=== alpha_n vs reference (NEW closed form) ===")
print()
print(f"  alpha_1 = {mp.nstr(alphas[0], 30)}")
print(f"    err = {mp.nstr(abs((alphas[0] - alpha1_ref)/alpha1_ref), 6)}")
print(f"  alpha_2 = {mp.nstr(alphas[1], 30)}")
print(f"    err = {mp.nstr(abs((alphas[1] - alpha2_ref_new)/alpha2_ref_new), 6)}")
for n in range(3, NMoments + 1):
    print(f"  alpha_{n} = {mp.nstr(alphas[n-1], 30)}")
print()

tau1_BEM = alphas[1] / (-alphas[0])
print(f"  tau_1   = {mp.nstr(tau1_BEM * 1e6, 30)} us")
print(f"  tau_1 err = {mp.nstr(abs((tau1_BEM - tau1_ref)/tau1_ref), 6)}")
print()

# === Save results ===
out = {
    "method": "Math K + Python eigvalue/Cauer (hybrid v6)",
    "P_max": P_MAX,
    "dps": SDPS,
    "alpha_n": [mp.nstr(a, 30) for a in alphas],
    "tau_1_us": mp.nstr(tau1_BEM * 1e6, 30),
    "alpha_2_ref_new": mp.nstr(alpha2_ref_new, 30),
    "tau_1_ref_us": mp.nstr(tau1_ref * 1e6, 30),
    "alpha_2_err": mp.nstr(abs((alphas[1] - alpha2_ref_new)/alpha2_ref_new), 6),
    "tau_1_err": mp.nstr(abs((tau1_BEM - tau1_ref)/tau1_ref), 6),
}
with open("bem_foster_spectral_v6_hybrid_results.json", "w") as f:
    json.dump(out, f, indent=2)
print("Results saved to bem_foster_spectral_v6_hybrid_results.json")
