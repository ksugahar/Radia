"""Spectral BEM v4 — Python (mpmath) version.

Demonstrates that the Mathematica v4 approach (1D-closed + 2D-NIntegrate)
runs equivalently in Python at arbitrary precision via mpmath.

Validates against Mathematica result:
  D_yy (single K element) at WP=50 should give 17+ digit agreement
  with v1 closed form reference.

Date: 2026-05-02
"""

import mpmath as mp
import time

# Set precision (50 decimal digits)
mp.mp.dps = 50

# Constants
sigma = mp.mpf("5.8e7")
mu0 = 4 * mp.pi * mp.mpf("1e-7")
a = mp.mpf("5e-3")
b = mp.mpf("2e-3")
c = mp.mpf("1e-3")
V = a * b * c
mu0sigma = mu0 * sigma

# References (Mathematica v1 closed form at WP=60)
alpha2_ref_new = mp.mpf("7.4783240535434074691212813799e-10")
Dyy_ref = mp.mpf("5.24918270035279310099903723219621782347721846299e-21")
Dxx_ref = mp.mpf("6.55125437660237267362579845589291716611473968610625e-20")

print("=" * 64)
print("Spectral BEM v4 (Python mpmath, dps=50)")
print("=" * 64)
print()
print(f"Reference D_yy (Mathematica WP=60) = {mp.nstr(Dyy_ref, 16)}")
print()

# === Inner u-integration for (p1, p2) = (0, 0) and (1, 1) ===
# gAxis[0, 0, h, xi] = h - xi  (for xi > 0)
# gAxis[1, 1, h, xi] = (h^3 - 3 h^2 xi + 2 xi^3) / 12  (for xi > 0)

def gAxis_00(h, xi):
    """g_x for p1=p2=0 (constant weight, after octant inner u integration)."""
    return h - xi

def gAxis_11(h, xi):
    """g_x for p1=p2=1."""
    return (h**3 - 3 * h**2 * xi + 2 * xi**3) / 12

# === 1D-closed inner xi integration ===
# I_k(h, R) = int_0^h xi^k / sqrt(xi^2 + R^2) dxi
# Closed forms:
#   I_0 = arcsinh(h/R)
#   I_1 = sqrt(h^2+R^2) - R
#   I_2 = (h*sqrt(h^2+R^2) - R^2 * arcsinh(h/R))/2
#   I_3 = ((h^2 - 2*R^2) * sqrt(h^2+R^2) + 2*R^3) / 3
#   I_k = (h^{k-1}*sqrt(h^2+R^2) - (k-1)*R^2*I_{k-2}) / k  for k >= 2

def I_xi(k, h, R):
    """Closed-form 1D integral int_0^h xi^k / sqrt(xi^2+R^2) dxi."""
    if k == 0:
        return mp.asinh(h / R)
    if k == 1:
        return mp.sqrt(h*h + R*R) - R
    # Recursion for k >= 2
    Ikm2 = I_xi(k - 2, h, R)
    return (h**(k - 1) * mp.sqrt(h*h + R*R) - (k - 1) * R*R * Ikm2) / k

def IxiClosed_00(h, R):
    """int_0^h gAxis_00(xi) / sqrt(xi^2 + R^2) dxi
       = int_0^h (h - xi)/sqrt(xi^2+R^2) dxi
       = h * I_0 - I_1
    """
    return h * I_xi(0, h, R) - I_xi(1, h, R)

def IxiClosed_11(h, R):
    """int_0^h gAxis_11(xi) / sqrt(xi^2 + R^2) dxi
       = (h^3 * I_0 - 3 h^2 * I_1 + 2 * I_3) / 12
    """
    return (h**3 * I_xi(0, h, R) - 3 * h**2 * I_xi(1, h, R) + 2 * I_xi(3, h, R)) / 12

# === Test 1: K[(0,1,0),(0,1,0)] = D_yy ===
# A_ext^x = -B_0 y/2; basis function (0,1,0) corresponds to y.
# Weights: g_x is gAxis_00 (h=a), g_y is gAxis_11 (h=b), g_z is gAxis_00 (h=c).
# K = 8 * int_0^b int_0^c IxiClosed_00(a, sqrt(eta^2+zeta^2)) * gAxis_11(b, eta) * gAxis_00(c, zeta) d(eta) d(zeta)

print("Test 1: K[(0,1,0),(0,1,0)] = D_yy")
def integrand_Dyy(eta, zeta):
    R = mp.sqrt(eta*eta + zeta*zeta)
    ixi = IxiClosed_00(a, R)
    gy = gAxis_11(b, eta)
    gz = gAxis_00(c, zeta)
    return ixi * gy * gz

t0 = time.time()
Dyy_test = 8 * mp.quad(integrand_Dyy, [0, b], [0, c])
dt = time.time() - t0
print(f"  D_yy (Python mpmath dps=50) = {mp.nstr(Dyy_test, 30)}")
print(f"  Time: {dt:.2f} s")
print(f"  Reference D_yy             = {mp.nstr(Dyy_ref, 30)}")
print(f"  Relative diff = {mp.nstr(abs((Dyy_test - Dyy_ref) / Dyy_ref), 6)}")
print()

# === Test 2: K[(1,0,0),(1,0,0)] = D_xx ===
print("Test 2: K[(1,0,0),(1,0,0)] = D_xx")
def integrand_Dxx(eta, zeta):
    R = mp.sqrt(eta*eta + zeta*zeta)
    ixi = IxiClosed_11(a, R)  # gAxis_11 for x
    gy = gAxis_00(b, eta)
    gz = gAxis_00(c, zeta)
    return ixi * gy * gz

t0 = time.time()
Dxx_test = 8 * mp.quad(integrand_Dxx, [0, b], [0, c])
dt = time.time() - t0
print(f"  D_xx (Python mpmath dps=50) = {mp.nstr(Dxx_test, 30)}")
print(f"  Time: {dt:.2f} s")
print(f"  Reference D_xx             = {mp.nstr(Dxx_ref, 30)}")
print(f"  Relative diff = {mp.nstr(abs((Dxx_test - Dxx_ref) / Dxx_ref), 6)}")
print()

# === Compute alpha_2 ===
alpha2_test = (mu0sigma)**2 / (16 * mp.pi * V) * (Dxx_test + Dyy_test)
print(f"alpha_2 (Python single-element) = {mp.nstr(alpha2_test, 30)}")
print(f"  vs alpha_2 closed (NEW)      = {mp.nstr(alpha2_ref_new, 30)}")
print(f"  Relative diff = {mp.nstr(abs((alpha2_test - alpha2_ref_new) / alpha2_ref_new), 6)}")
