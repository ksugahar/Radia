"""Spectral BEM v5 — Python production version (mpmath dps=50).

This is the PRIMARY production implementation.
Mathematica equivalent (bem_foster_spectral_v5.wls) serves as VERIFICATION ONLY.

Pipeline:
  1. gAxis polynomial (SymPy at startup, cached)
  2. IxiClosed inner-xi integration (recursive closed form via mpmath)
  3. Kmon[(p,q,r),(p',q',r')] via 2D mpmath.quad
  4. Legendre basis transformation (mass matrix diagonal)
  5. Standard eigenvalue problem via mpmath.eig
  6. Project A_ext onto Foster modes
  7. Compute alpha_n (n=1..8) via Foster sum

Validates against Mathematica v5 result:
  alpha_2 = 7.4783240535434074561e-10 (relative error ~1.7e-18 vs new closed)
  tau_1 = 16.98281435193036473 us

Date: 2026-05-02
"""

import time
import json

import mpmath as mp
import sympy as sp

# === Set precision (50 decimal digits, restore default) ===
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

# References (Mathematica WP=60 closed form, NEW value)
alpha1_ref = -mu0sigma * (a*a + b*b) / 48
alpha2_ref_new = mp.mpf("7.4783240535434074691212813799e-10")
tau1_ref = alpha2_ref_new / (-alpha1_ref)

print("=" * 64)
print("Spectral BEM v5 (Python mpmath, dps={})".format(SDPS))
print("=" * 64)
print()
print("References:")
print("  alpha_1 = {}".format(mp.nstr(alpha1_ref, 18)))
print("  alpha_2 = {}".format(mp.nstr(alpha2_ref_new, 18)))
print("  tau_1   = {} us".format(mp.nstr(tau1_ref * 1e6, 18)))
print()

# === gAxis polynomial coefficients (SymPy derivation, cached) ===
_h_sym, _xi_sym, _u_sym = sp.symbols('h xi u', positive=True)
_gAxis_cache_h = {}  # symbolic: (p1, p2) -> list of coeffs in h
_gAxis_cache_h_eval = {}  # numerical: (p1, p2, h_mpf) -> list of mpf


def gAxis_xi_coeffs_symbolic(p1, p2):
    """Symbolic polynomial coefficients of xi^k (ascending order)."""
    if (p1, p2) not in _gAxis_cache_h:
        expr = sp.integrate(
            (_u_sym + _xi_sym/2)**p1 * (_u_sym - _xi_sym/2)**p2,
            (_u_sym, _xi_sym/2 - _h_sym/2, _h_sym/2 - _xi_sym/2))
        expr = sp.expand(expr)
        poly_xi = sp.Poly(expr, _xi_sym)
        coeffs = poly_xi.all_coeffs()[::-1]  # ascending
        _gAxis_cache_h[(p1, p2)] = coeffs
    return _gAxis_cache_h[(p1, p2)]


def gAxis_xi_coeffs(p1, p2, h_mp):
    """mpmath coefficients of xi^k after substituting numerical h."""
    key = (p1, p2, str(h_mp))
    if key not in _gAxis_cache_h_eval:
        sym_coeffs = gAxis_xi_coeffs_symbolic(p1, p2)
        h_subs = sp.Float(mp.nstr(h_mp, SDPS), SDPS)
        result = []
        for c in sym_coeffs:
            v = c.subs(_h_sym, h_subs).evalf(SDPS)
            result.append(mp.mpf(str(v)))
        _gAxis_cache_h_eval[key] = result
    return _gAxis_cache_h_eval[key]


# === Inner xi closed-form integration ===
_I_xi_cache = {}


def I_xi(k, h, R):
    """Closed form: int_0^h xi^k / sqrt(xi^2 + R^2) dxi.

    Recursion: I_k = (h^{k-1} sqrt(h^2+R^2) - (k-1) R^2 I_{k-2}) / k for k >= 2.
    Base: I_0 = arcsinh(h/R), I_1 = sqrt(h^2+R^2) - R.
    """
    key = (k, str(h), str(R))
    if key in _I_xi_cache:
        return _I_xi_cache[key]
    if k == 0:
        v = mp.asinh(h / R)
    elif k == 1:
        v = mp.sqrt(h*h + R*R) - R
    else:
        Ikm2 = I_xi(k - 2, h, R)
        v = (h**(k - 1) * mp.sqrt(h*h + R*R) - (k - 1) * R*R * Ikm2) / k
    _I_xi_cache[key] = v
    return v


def IxiClosed_poly(coeffs_xi, h, R):
    """Sum_k c_k * I_xi(k, h, R)."""
    result = mp.mpf(0)
    for k, c in enumerate(coeffs_xi):
        if c != 0:
            result += c * I_xi(k, h, R)
    return result


def poly_eval(coeffs, x):
    """Horner-style polynomial evaluation."""
    result = mp.mpf(0)
    for c in reversed(coeffs):
        result = result * x + c
    return result


# === K element (cached) ===
_Kmon_cache = {}


def Kmon(p1q1r1, p2q2r2):
    """8 * int_0^a int_0^b int_0^c gAxis_x gAxis_y gAxis_z / r dV (octant)."""
    key = (p1q1r1, p2q2r2)
    if key in _Kmon_cache:
        return _Kmon_cache[key]

    p1, q1, r1 = p1q1r1
    p2, q2, r2 = p2q2r2
    cx = gAxis_xi_coeffs(p1, p2, a)
    cy = gAxis_xi_coeffs(q1, q2, b)
    cz = gAxis_xi_coeffs(r1, r2, c)

    def integrand(eta, zeta):
        R = mp.sqrt(eta*eta + zeta*zeta)
        ixi = IxiClosed_poly(cx, a, R)
        gy = poly_eval(cy, eta)
        gz = poly_eval(cz, zeta)
        return ixi * gy * gz

    # mpmath default Tanh-Sinh (best for our case among tried)
    val = 8 * mp.quad(integrand, [0, b], [0, c])
    _Kmon_cache[key] = val
    return val


# === Build basis (e,o,e) parity for x-coupling, (o,e,e) for y-coupling ===
def build_basis(P_max, parity):
    px, py, pz = parity
    basis = []
    for p in range(P_max + 1):
        for q in range(P_max - p + 1):
            for r in range(P_max - p - q + 1):
                if p % 2 == px and q % 2 == py and r % 2 == pz:
                    basis.append((p, q, r))
    return basis


# === Legendre transformation matrix ===
_legendreT_cache = {}


def legendreT(P_max, h):
    """T[P, p] = coefficient of x^p in L_P(2x/h) for x in [-h/2, h/2]."""
    key = (P_max, str(h))
    if key in _legendreT_cache:
        return _legendreT_cache[key]
    x = sp.Symbol('x')
    T = mp.matrix(P_max + 1, P_max + 1)
    h_sub = sp.Float(mp.nstr(h, SDPS), SDPS)
    for n in range(P_max + 1):
        Ln = sp.legendre(n, x)
        coeffs = sp.Poly(Ln, x).all_coeffs()[::-1]
        for k, ck in enumerate(coeffs):
            v = ck * (sp.Float(2, SDPS) / h_sub)**k
            T[n, k] = mp.mpf(str(v.evalf(SDPS)))
    _legendreT_cache[key] = T
    return T


# === Solve one parity block ===
def solve_block(P_max, parity, label):
    basis = build_basis(P_max, parity)
    n = len(basis)
    print(f"  Block {label} at P={P_max}: {n} functions")

    print(f"    Building K_mon matrix ({n*(n+1)//2} unique)...")
    t0 = time.time()
    K_mon = mp.matrix(n, n)
    for i in range(n):
        for j in range(i, n):
            val = Kmon(basis[i], basis[j])
            K_mon[i, j] = val
            K_mon[j, i] = val
    print(f"      K_mon time: {time.time()-t0:.2f} s")

    Tx = legendreT(P_max, a)
    Ty = legendreT(P_max, b)
    Tz = legendreT(P_max, c)
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

    print(f"    Solving symmetric eigenvalue (mpmath.eighe)...")
    t0 = time.time()
    # K_leg_std is symmetric (forced) - use Hermitian eigvalue solver
    # which returns orthonormal eigenvectors with sorted real eigenvalues
    try:
        eigvals, eigvecs = mp.eighe(K_leg_std)
    except (AttributeError, NameError):
        # Fallback to mp.eig + manual orthogonalization
        eigvals, eigvecs = mp.eig(K_leg_std)
        # Verify and normalize each eigenvector
        for k in range(n):
            norm = mp.sqrt(mp.fsum(eigvecs[i, k]**2 for i in range(n)))
            for i in range(n):
                eigvecs[i, k] = eigvecs[i, k] / norm
        # Note: mp.eig should return mutually orthogonal vectors for symmetric matrix
        # but normalization isn't guaranteed; we just normalized.
    print(f"      eig solve: {time.time()-t0:.2f} s")

    # Sanity: verify orthonormality
    max_orth_err = mp.mpf(0)
    for i in range(min(n, 4)):
        for j in range(min(n, 4)):
            dot = mp.fsum(eigvecs[k, i] * eigvecs[k, j] for k in range(n))
            target = mp.mpf(1) if i == j else mp.mpf(0)
            err = abs(dot - target)
            if err > max_orth_err:
                max_orth_err = err
    print(f"      max orthonormality error: {mp.nstr(max_orth_err, 4)}")

    return basis, eigvals, eigvecs, massDiag


def project_aext(basis, massDiag, axis):
    """Project A_ext on M-orthonormal Legendre basis.

    A_ext^x = -B_0 y/2 = (-b/4) L_1(2y/b) so coef at (0,1,0) Legendre is -b/4.
    M-orthonormal: c_norm = c * sqrt(mass[(0,1,0)]).
    """
    n = len(basis)
    vec = mp.matrix(n, 1)
    if axis == "x":
        target = (0, 1, 0)
        coef = -b / 4
    elif axis == "y":
        target = (1, 0, 0)
        coef = a / 4
    if target in basis:
        idx = basis.index(target)
        vec[idx, 0] = coef * mp.sqrt(massDiag[idx])
    return vec


def runForPmax(P_max):
    print(f"=== P_max = {P_max} ===")

    bX, eX, vX, mX = solve_block(P_max, (0, 1, 0), "x-coupling (e,o,e)")
    bY, eY, vY, mY = solve_block(P_max, (1, 0, 0), "y-coupling (o,e,e)")

    aExt_X = project_aext(bX, mX, "x")
    aExt_Y = project_aext(bY, mY, "y")

    # Foster amplitudes: g_k = sum_I phi_kI * aExt_I
    # vX, vY have eigenvectors as COLUMNS (mpmath convention).
    # We want g_k = phi_k . aExt for each k. So g = vX.T * aExt_X (as matrix)
    nX = len(bX)
    nY = len(bY)
    gX = [mp.fsum(vX[i, k] * aExt_X[i, 0] for i in range(nX)) for k in range(nX)]
    gY = [mp.fsum(vY[i, k] * aExt_Y[i, 0] for i in range(nY)) for k in range(nY)]

    lambdaX = [eX[k] / (4 * mp.pi) for k in range(nX)]
    lambdaY = [eY[k] / (4 * mp.pi) for k in range(nY)]

    NMoments = 8
    alphas = []
    for nM in range(1, NMoments + 1):
        sX = mp.fsum(gX[k]**2 * lambdaX[k]**(nM - 1) for k in range(nX))
        sY = mp.fsum(gY[k]**2 * lambdaY[k]**(nM - 1) for k in range(nY))
        a_n = (-1)**nM * mu0sigma**nM / V * (sX + sY)
        alphas.append(a_n)

    print()
    print("  Eigenvalues (x-block, sorted descending):")
    eX_sorted = sorted([e.real if hasattr(e, 'real') else e for e in eX], reverse=True)
    for k, ev in enumerate(eX_sorted[:6]):
        print(f"    eX[{k}] = {mp.nstr(ev, 25)}")
    print("  Eigenvalues (y-block, sorted descending):")
    eY_sorted = sorted([e.real if hasattr(e, 'real') else e for e in eY], reverse=True)
    for k, ev in enumerate(eY_sorted[:6]):
        print(f"    eY[{k}] = {mp.nstr(ev, 25)}")
    print()
    print("  alpha_n vs reference (NEW closed):")
    print(f"    alpha_1 = {mp.nstr(alphas[0], 18)}")
    print(f"      err = {mp.nstr(abs((alphas[0] - alpha1_ref)/alpha1_ref), 6)}")
    print(f"    alpha_2 = {mp.nstr(alphas[1], 18)}")
    print(f"      err = {mp.nstr(abs((alphas[1] - alpha2_ref_new)/alpha2_ref_new), 6)}")
    for n in range(3, NMoments + 1):
        print(f"    alpha_{n} = {mp.nstr(alphas[n-1], 18)}")

    tau1_BEM = alphas[1] / (-alphas[0])
    print(f"    tau_1   = {mp.nstr(tau1_BEM * 1e6, 18)} us")
    print(f"    tau_1 err = {mp.nstr(abs((tau1_BEM - tau1_ref)/tau1_ref), 6)}")
    print()
    return alphas, (bX, eX, vX, mX), (bY, eY, vY, mY)


# === Run P=4 ===
total_t0 = time.time()
result4 = runForPmax(4)
total_dt = time.time() - total_t0
print(f"Total runtime: {total_dt:.1f} s")

# Save results
alphas4, (bX, eX, vX, mX), (bY, eY, vY, mY) = result4
with open("bem_foster_spectral_v5_python_results.json", "w") as f:
    json.dump({
        "P_max": 4,
        "alpha_n": [mp.nstr(a, 30) for a in alphas4],
        "tau_1_us": mp.nstr(alphas4[1] / (-alphas4[0]) * 1e6, 18),
        "alpha_2_ref": mp.nstr(alpha2_ref_new, 30),
        "tau_1_ref_us": mp.nstr(tau1_ref * 1e6, 18),
        "runtime_s": total_dt,
        "dps": SDPS,
    }, f, indent=2)
print("Results saved to bem_foster_spectral_v5_python_results.json")
