"""Python log removal v2: use mpmath polar integration for INT_log.

Strategy:
  Convert ∫∫ log(R) f(η,ζ) dη dζ over rectangle to polar coordinates,
  splitting into two triangular regions (no domain issues).

  Region 1: phi in [0, atan(c/b)], rho in [0, b/cos(phi)]
  Region 2: phi in [atan(c/b), π/2], rho in [0, c/sin(phi)]

  In polar, log(R) integration becomes:
    ∫ rho * log(rho) drho * ∫ f(rho cos phi, rho sin phi) dphi
  The 1D rho * log(rho) integral has closed form (eliminates corner singularity).

Goal: 30+ digit D_yy without SymPy bug.

Date: 2026-05-02
"""

import time
import mpmath as mp

mp.mp.dps = 50

a = mp.mpf("5e-3")
b = mp.mpf("2e-3")
c = mp.mpf("1e-3")

print("=" * 64)
print("Python (mpmath polar): D_yy via log removal v2")
print("=" * 64)
print()


def ixi_smooth(R):
    """Smooth part of ixi after log removal."""
    return a * mp.log(a + mp.sqrt(a*a + R*R)) - mp.sqrt(a*a + R*R) + R


def g_y(eta):
    return (b**3 - 3*b**2*eta + 2*eta**3) / 12


def g_z(zeta):
    return c - zeta


# === Compute INT_log via polar coordinates ===
# INT_log = ∫_0^b ∫_0^c log(sqrt(η²+ζ²)) g_y(η) g_z(ζ) dη dζ
# Polar: η = ρ cos φ, ζ = ρ sin φ, dη dζ = ρ dρ dφ
# Domain split at φ = atan(c/b):
#   Region 1 (0 <= phi <= atan(c/b)):  rho_max = b/cos(phi)  (limited by η = b)
#   Region 2 (atan(c/b) <= phi <= π/2): rho_max = c/sin(phi)  (limited by ζ = c)

phi_split = mp.atan(c / b)

def integrand_polar_log(phi, rho):
    eta = rho * mp.cos(phi)
    zeta = rho * mp.sin(phi)
    if rho == 0: return 0
    return mp.log(rho) * g_y(eta) * g_z(zeta) * rho


print("Computing INT_log via polar (region 1: 0 < phi < atan(c/b))...")
t0 = time.time()
INT_log_R1 = mp.quad(
    lambda phi: mp.quad(lambda rho: integrand_polar_log(phi, rho),
                        [0, b / mp.cos(phi)]),
    [0, phi_split]
)
print(f"  R1 = {mp.nstr(INT_log_R1, 30)}, time {time.time()-t0:.2f} s")

print("Computing INT_log via polar (region 2: atan(c/b) < phi < π/2)...")
t0 = time.time()
INT_log_R2 = mp.quad(
    lambda phi: mp.quad(lambda rho: integrand_polar_log(phi, rho),
                        [0, c / mp.sin(phi)]),
    [phi_split, mp.pi / 2]
)
print(f"  R2 = {mp.nstr(INT_log_R2, 30)}, time {time.time()-t0:.2f} s")

INT_log = INT_log_R1 + INT_log_R2
print(f"  INT_log (polar) = {mp.nstr(INT_log, 30)}")
print(f"  MATLAB ref      = -1.52143398669749040857e-19")
print()

# === Compute INT_smooth (Cartesian, mp.quad) ===
print("Computing INT_smooth (Cartesian mp.quad)...")
t0 = time.time()
def integrand_smooth_cart(eta, zeta):
    R = mp.sqrt(eta*eta + zeta*zeta)
    return ixi_smooth(R) * g_y(eta) * g_z(zeta)

INT_smooth = mp.quad(integrand_smooth_cart, [0, b], [0, c])
print(f"  INT_smooth = {mp.nstr(INT_smooth, 30)}, time {time.time()-t0:.2f} s")
print()

# === Combine ===
K_log_part = -8 * a * INT_log
K_smooth_part = 8 * INT_smooth
D_yy = K_log_part + K_smooth_part

print("=" * 64)
print("Result")
print("=" * 64)
print(f"  D_yy (Python polar log removal) = {mp.nstr(D_yy, 35)}")

ref = mp.mpf("5.24918270035279297821091928151251136119220629893627828114e-21")
print(f"  Reference (Math WP=50)          = {mp.nstr(ref, 35)}")
err = abs((D_yy - ref) / ref)
print(f"  Relative error = {mp.nstr(err, 6)}")

if err < mp.mpf("1e-30"):
    n_digits = -mp.log10(err)
    print(f"  SUCCESS: {mp.nstr(n_digits, 4)} digit precision (>= 30)")
elif err < mp.mpf("1e-20"):
    n_digits = -mp.log10(err)
    print(f"  PARTIAL: {mp.nstr(n_digits, 4)} digit precision (better than 20)")
else:
    n_digits = -mp.log10(err)
    print(f"  FAIL: only {mp.nstr(n_digits, 4)} digit precision")
