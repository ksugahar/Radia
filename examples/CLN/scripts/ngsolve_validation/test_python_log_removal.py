"""Python analytic log singularity removal test.

Strategy (same as MATLAB version):
  ixi(R) = a*arcsinh(a/R) - sqrt(a^2+R^2) + R
         = -a*log(R) + ixi_smooth(R)
  where ixi_smooth(R) = a*log(a + sqrt(a^2+R^2)) - sqrt(a^2+R^2) + R (smooth at R=0)

  K = -8a * INT_log + 8 * INT_smooth
  INT_log = int int log(R) g_y g_z deta dzeta -> closed form via SymPy
  INT_smooth = int int ixi_smooth(R) g_y g_z deta dzeta -> mp.quad (smooth, fast)

Goal: D_yy = K[(0,1,0),(0,1,0)] should match Mathematica's
  5.24918270035279297821091928151251136119220629893627828114e-21 to 30+ digits.

Date: 2026-05-02
"""

import time
import mpmath as mp
import sympy as sp

mp.mp.dps = 50

a = mp.mpf("5e-3")
b = mp.mpf("2e-3")
c = mp.mpf("1e-3")

print("=" * 64)
print("Python (sympy + mpmath): analytic log removal D_yy test")
print("=" * 64)
print()

# === Part 1: log integral via SymPy ===
print("Part 1: int int log(R) * g_y * g_z deta dzeta (SymPy closed form)...")
t0 = time.time()
eta_s, zeta_s = sp.symbols('eta zeta', positive=True)
a_s = sp.Rational(5, 1000)
b_s = sp.Rational(2, 1000)
c_s = sp.Rational(1, 1000)

g_y_s = (b_s**3 - 3*b_s**2*eta_s + 2*eta_s**3) / 12  # for (q1,q2)=(1,1)
g_z_s = c_s - zeta_s                                  # for (r1,r2)=(0,0)

log_R = sp.log(sp.sqrt(eta_s**2 + zeta_s**2))

# Inner integral over zeta
inner_log = sp.integrate(log_R * g_y_s * g_z_s, (zeta_s, 0, c_s))
print(f"  Inner integral done: {time.time()-t0:.2f} s")

# Outer integral over eta
outer_log = sp.integrate(inner_log, (eta_s, 0, b_s))
print(f"  Outer integral done: {time.time()-t0:.2f} s")

INT_log_sp = outer_log
INT_log = mp.mpf(str(INT_log_sp.evalf(50)))
print(f"  INT_log = {mp.nstr(INT_log, 30)}")
print()

# === Part 2: smooth integral (mp.quad) ===
print("Part 2: int int ixi_smooth(R) * g_y * g_z deta dzeta (mp.quad)...")
t0 = time.time()


def ixi_smooth(R):
    return a * mp.log(a + mp.sqrt(a*a + R*R)) - mp.sqrt(a*a + R*R) + R


def g_y_num(eta):
    return (b**3 - 3*b**2*eta + 2*eta**3) / 12


def g_z_num(zeta):
    return c - zeta


def integrand_smooth(eta, zeta):
    R = mp.sqrt(eta*eta + zeta*zeta)
    return ixi_smooth(R) * g_y_num(eta) * g_z_num(zeta)


INT_smooth = mp.quad(integrand_smooth, [0, b], [0, c])
print(f"  INT_smooth = {mp.nstr(INT_smooth, 30)}")
print(f"  Time: {time.time()-t0:.2f} s")
print()

# === Combine ===
K_log = -8 * a * INT_log
K_smooth = 8 * INT_smooth
D_yy_total = K_log + K_smooth

print("=" * 64)
print("Final result")
print("=" * 64)
print(f"  D_yy (Python log removal)  = {mp.nstr(D_yy_total, 35)}")

ref = mp.mpf("5.24918270035279297821091928151251136119220629893627828114e-21")
print(f"  Reference (Math WP=50)     = {mp.nstr(ref, 35)}")
print(f"  Relative diff = {mp.nstr(abs((D_yy_total - ref)/ref), 6)}")

# Compare with new closed reference (alpha_2 closed)
alpha2_ref_new = mp.mpf("7.4783240535434074691212813799e-10")
mu0sigma = 4 * mp.pi * mp.mpf("1e-7") * mp.mpf("5.8e7")
V = a * b * c

# Need both D_xx and D_yy for alpha_2; here only D_yy. Just check D_yy alone.
# Or proceed to log_removal_K_full computation later.

print()
print("=" * 64)
print("Diagnostic: did log removal achieve 30+ digits?")
print("=" * 64)
err = abs((D_yy_total - ref) / ref)
if err < mp.mpf("1e-30"):
    print(f"SUCCESS: Python achieves {-mp.log10(err):.1f} digit precision (>= 30)")
elif err < mp.mpf("1e-20"):
    print(f"PARTIAL: Python achieves {-mp.log10(err):.1f} digit precision (better than 18 but < 30)")
elif err < mp.mpf("1e-15"):
    print(f"BARELY: Python achieves {-mp.log10(err):.1f} digit precision (similar to mp.quad alone)")
else:
    print(f"FAIL: Python achieves only {-mp.log10(err):.1f} digit precision")
