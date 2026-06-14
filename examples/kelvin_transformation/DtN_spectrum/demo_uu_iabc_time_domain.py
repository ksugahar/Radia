# -*- coding: utf-8 -*-
r"""
demo_uu_iabc_time_domain.py  (Track A -- kelvin branch)
=======================================================
TIME-DOMAIN representation of the exact spherical exterior DtN that IABC
(nested isotropic absorbing shells) approximates.

WHY A CAS (the user's instinct that "Mathematica's power" is needed):
the per-mode exterior Dirichlet-to-Neumann symbol

      Lambda_n(z) = z * h_n^{(1)'}(z) / h_n^{(1)}(z),     z = k R = omega R / c

is a RATIONAL function of z, because the spherical Hankel function is a degree-n
polynomial times e^{i z} / z^{n+1}:

      h_n^{(1)}(z) = (-i)^{n+1} * e^{i z} * z^{-(n+1)} * theta_n(-i z),

theta_n = reverse Bessel polynomial.  Symbolic partial fraction + inverse Fourier
transform turns the *frequency-domain* boundary operator into a *local-in-time*
auxiliary-ODE system.  That symbolic step is what a CAS is for; Mathematica is not
installed on this host, so we use sympy and VERIFY every result numerically.

PROVEN HERE (all asserted):
  (1) Closed form  h_n^{(1)}(z) = (-i)^{n+1} e^{iz} z^{-(n+1)} theta_n(-iz)
      matches scipy's spherical Hankel  (rel.err ~ machine).
  (2) Pole-residue DtN:  Lambda_n(z) = i z - 1 + sum_j z_j/(z - z_j),
      z_j = i*(roots of theta_n) = zeros of h_n^{(1)}.  Matches scipy DtN.
  (3) THREE INDEPENDENT pole sets agree:  sympy roots(theta_n),  the zeros of the
      scipy spherical Hankel, and scipy.signal.besselap (Bessel/Thomson filter).
      => the time-domain spherical open boundary is a BESSEL-FILTER NETWORK per mode.
  (4) Time-domain realization (Grote-Keller form, here R=c=1):
          g(t) = -du/dt - u + sum_j psi_j           (g = R d_r u, the DtN output)
          dpsi_j/dt = -i z_j (psi_j + u)             (one auxiliary ODE per pole)
      reproduces Lambda_n(z) for a steady tone -- checked BOTH algebraically and by
      transient time integration from rest.
  (5) The auxiliary relaxation rates lambda_j = -i z_j EQUAL the Bessel-filter poles
      and all have Re < 0  => the realization is STABLE / causal.

HONEST PRIOR ART: the rational/pole-residue exact sphere DtN and its local-in-time
auxiliary realization is Grote & Keller (1995) / Hagstrom; the reverse-Bessel-
polynomial == Bessel/Thomson-filter-pole identity is classical network synthesis.
The new angle (recorded in the knowledge module, NOT asserted as a theorem here) is
reading IABC's N shells AS this N-pole set (datasheet: N shells -> multipoles to
order N) and the Kelvin/DtN-spectral linkage.  This file supplies only the VERIFIED
time-domain core; it does NOT claim IABC's exact shell constants (those live in
Sugahara ICEAA-2015, not in hand).

No overclaim: every printed 'ok' is gated on an executed numerical assertion.
"""
import os
os.environ.setdefault("PYTHONIOENCODING", "utf-8")
import numpy as np
import sympy as sp
from scipy.special import hankel1
from scipy.signal import besselap
from scipy.integrate import solve_ivp

np.set_printoptions(precision=6, suppress=False)

# ===========================================================================
# scipy reference: exact spherical Hankel of the first kind and the exact DtN
#   h_n^{(1)}(z) = sqrt(pi/(2z)) * H^{(1)}_{n+1/2}(z)   (hankel1 accepts complex z)
# ===========================================================================
def sph_h1(n, z):
    z = np.asarray(z, dtype=complex)
    return np.sqrt(np.pi / (2.0 * z)) * hankel1(n + 0.5, z)

def sph_h1_prime(n, z):
    # spherical-Bessel derivative recurrence: f_n' = f_{n-1} - (n+1)/z f_n
    z = np.asarray(z, dtype=complex)
    return sph_h1(n - 1, z) - (n + 1) / z * sph_h1(n, z)

def dtn_exact(n, z):
    z = np.asarray(z, dtype=complex)
    return (z * sph_h1_prime(n, z) / sph_h1(n, z))

# ===========================================================================
# (1)+(2) SYMBOLIC: reverse Bessel polynomial, closed-form Hankel, pole-residue DtN
# ===========================================================================
xsym = sp.symbols('x')

def reverse_bessel_poly_expr(n):
    # theta_n(x) = sum_{k=0}^n (n+k)! / ((n-k)! k! 2^k) * x^{n-k}
    return sp.expand(sum(
        sp.factorial(n + k) / (sp.factorial(n - k) * sp.factorial(k) * 2**k) * xsym**(n - k)
        for k in range(n + 1)))

def theta_roots(n):
    """Roots of theta_n(x) as complex numpy array (the 'x' variable)."""
    p = sp.Poly(reverse_bessel_poly_expr(n), xsym)
    return np.array([complex(r) for r in p.nroots(n=30)], dtype=complex)

def hankel_via_theta(n, zval):
    """Closed-form spherical Hankel built from the reverse Bessel polynomial.
    h_n^{(1)}(z) = -i * e^{iz} * z^{-(n+1)} * theta_n(-iz).
    The leading constant is -i for EVERY n: the naive (-i)^{n+1} prefactor of the
    descending-series form is cancelled by theta_n(-iz) = (-i)^n * Q_n(z), where
    Q_n is the actual polynomial factor of h_n; (-i)^{n+1} * i^n = -i."""
    th = reverse_bessel_poly_expr(n)
    val = complex(th.subs(xsym, -sp.I * zval))
    return (-1j) * np.exp(1j * zval) * zval**(-(n + 1)) * val

def dtn_pole_residue(n, z, zpoles):
    """Lambda_n(z) = i z - 1 + sum_j z_j/(z - z_j)."""
    z = complex(z)
    return 1j * z - 1.0 + np.sum(zpoles / (z - zpoles))

NMAX = 6
# test points (complex), kept away from any pole
ZTEST = [0.5, 1.0, 2.0, 4.0, 0.3 + 0.0j, 1.0 + 0.4j, 2.5 - 0.6j, 3.0 + 1.0j]

print("=" * 78)
print(" demo_uu : TIME-DOMAIN representation of the exact spherical (IABC-target) DtN")
print("=" * 78)

print("\n[1] closed form  h_n = -i * e^{iz} z^{-(n+1)} theta_n(-iz)  vs scipy:")
max_h = 0.0
for n in range(1, NMAX + 1):
    for zt in ZTEST:
        a = hankel_via_theta(n, zt)
        b = complex(sph_h1(n, zt))
        rel = abs(a - b) / abs(b)
        max_h = max(max_h, rel)
print(f"    n=1..{NMAX}, {len(ZTEST)} pts each : max rel.err = {max_h:.2e}")
assert max_h < 1e-10, "closed-form Hankel mismatch"
print("    ok  (theta_n closed form reproduces the spherical Hankel)")

print("\n[2] poles z_j = i*roots(theta_n) ARE zeros of h_n^{(1)}  (so DtN poles):")
poles_z = {}
max_zero = 0.0
for n in range(1, NMAX + 1):
    zj = 1j * theta_roots(n)
    poles_z[n] = zj
    # h_n(z_j) should vanish; normalize by a nearby non-zero scale to test 'zero-ness'
    hv = np.array([complex(sph_h1(n, z)) for z in zj])
    scale = abs(complex(sph_h1(n, 1.0)))
    err = np.max(np.abs(hv)) / scale
    max_zero = max(max_zero, err)
print(f"    n=1..{NMAX} : max |h_n(z_j)|/|h_n(1)| = {max_zero:.2e}")
assert max_zero < 1e-8, "theta-roots are not Hankel zeros"
print("    ok  (the n reverse-Bessel roots are exactly the n Hankel zeros)")

print("\n[3] pole-residue DtN   Lambda_n = iz - 1 + sum z_j/(z-z_j)   vs scipy DtN:")
max_dtn = 0.0
for n in range(1, NMAX + 1):
    for zt in ZTEST:
        a = dtn_pole_residue(n, zt, poles_z[n])
        b = complex(dtn_exact(n, zt))
        max_dtn = max(max_dtn, abs(a - b) / abs(b))
print(f"    n=1..{NMAX}, {len(ZTEST)} pts each : max rel.err = {max_dtn:.2e}")
assert max_dtn < 1e-9, "pole-residue DtN mismatch"
print("    ok  (exact per-mode DtN is the degree-n rational  iz-1+sum z_j/(z-z_j))")

print("\n[3b] CROSS-CHECK: roots(theta_n)  vs  scipy.signal.besselap (Bessel filter):")
max_bes = 0.0
for n in range(1, NMAX + 1):
    rx = np.sort_complex(theta_roots(n))
    _, pb, _ = besselap(n, norm='delay')   # poles of the analog Bessel/Thomson filter
    pb = np.sort_complex(np.asarray(pb, dtype=complex))
    err = np.max(np.abs(rx - pb))
    max_bes = max(max_bes, err)
print(f"    n=1..{NMAX} : max |roots(theta_n) - besselap poles| = {max_bes:.2e}")
assert max_bes < 1e-6, "reverse-Bessel roots != Bessel filter poles"
print("    ok  (=> the spherical open boundary's poles ARE Bessel/Thomson filter poles)")

# ===========================================================================
# (4) TIME-DOMAIN realization (R=c=1):  g = -u' - u + sum psi_j ;
#                                       psi_j' = -i z_j (psi_j + u)
# (a) algebraic steady-state for a tone u=e^{-i w t}  =>  psi_j = z_j/(w - z_j) u
#     and  g/u = i w - 1 + sum z_j/(w - z_j) = Lambda_n(w).
# (b) transient integration from rest must converge to the same DtN.
# ===========================================================================
print("\n[4a] time-domain auxiliary-ODE realization, ALGEBRAIC steady state (tone):")
max_alg = 0.0
for n in range(1, NMAX + 1):
    zj = poles_z[n]
    for w in [0.5, 1.0, 2.0, 3.5]:
        psi = zj / (w - zj)                 # steady-state amplitude (per pole)
        g_over_u = 1j * w - 1.0 + np.sum(psi)
        max_alg = max(max_alg, abs(g_over_u - complex(dtn_exact(n, w))) / abs(complex(dtn_exact(n, w))))
print(f"    n=1..{NMAX} : max rel.err vs exact DtN = {max_alg:.2e}")
assert max_alg < 1e-9, "algebraic steady-state mismatch"
print("    ok  (the auxiliary ODEs reproduce the DtN for every tone)")

print("\n[4b] TRANSIENT: integrate psi_j' = -i z_j (psi_j+u) from rest, read g/u at steady state:")
def transient_dtn(n, w, t_end=400.0):
    zj = poles_z[n]
    m = zj.size
    def rhs(t, y):
        psi = y[:m] + 1j * y[m:]
        u = np.exp(-1j * w * t)
        dpsi = -1j * zj * (psi + u)
        return np.concatenate([dpsi.real, dpsi.imag])
    y0 = np.zeros(2 * m)
    sol = solve_ivp(rhs, [0.0, t_end], y0, rtol=1e-9, atol=1e-11,
                    t_eval=[t_end], method='RK45')
    psi = sol.y[:m, -1] + 1j * sol.y[m:, -1]
    u = np.exp(-1j * w * t_end)
    up = -1j * w * u                          # du/dt
    g = -up - u + np.sum(psi)
    return g / u
max_tr = 0.0
for n in [1, 2, 3]:
    for w in [1.0, 2.0]:
        g_over_u = transient_dtn(n, w)
        max_tr = max(max_tr, abs(g_over_u - complex(dtn_exact(n, w))) / abs(complex(dtn_exact(n, w))))
print(f"    n=1..3, w in {{1,2}} : max rel.err vs exact DtN = {max_tr:.2e}")
assert max_tr < 1e-3, "transient realization did not converge to the DtN"
print("    ok  (time integration of the auxiliary network converges to the exact DtN)")

# ===========================================================================
# (5) STABILITY / causality:  relaxation rate lambda_j = -i z_j = root(theta_n)
#     == Bessel filter pole, Re < 0.
# ===========================================================================
print("\n[5] stability: auxiliary relaxation rates lambda_j = -i z_j  (== theta roots):")
worst_re = -np.inf
for n in range(1, NMAX + 1):
    lam = -1j * poles_z[n]                    # = roots(theta_n)
    worst_re = max(worst_re, np.max(lam.real))
    rmatch = np.max(np.abs(np.sort_complex(lam) - np.sort_complex(theta_roots(n))))
    assert rmatch < 1e-9
print(f"    n=1..{NMAX} : worst Re(lambda_j) = {worst_re:.3e}  (must be < 0)")
assert worst_re < 0.0, "an auxiliary mode is unstable"
print("    ok  (all auxiliary ODEs decay => stable, causal time-domain boundary)")

# ---- explicit small-n formulas (for the manuscript) -----------------------
print("\n[info] explicit pole/rate tables (z_j = pole in z=wR/c ; lambda_j = -i z_j = ODE rate):")
for n in [1, 2, 3]:
    zj = poles_z[n]
    lam = -1j * zj
    print(f"    n={n}:  z_j = {np.round(zj, 4)}")
    print(f"           lambda_j (=Bessel poles) = {np.round(lam, 4)}")

print("\n[interpretation]")
print("  * Exact spherical open boundary, mode n: DtN is the degree-n rational")
print("    Lambda_n(z) = iz - 1 + sum_{j=1..n} z_j/(z - z_j).")
print("  * Time domain = n LOCAL auxiliary ODEs (Grote-Keller form); their rates")
print("    are the Bessel/Thomson filter poles => a Bessel-filter NETWORK per mode.")
print("  * n poles  <=>  multipole order n  (the DtN datasheet, in the TIME domain).")
print("  * IABC's N nested shells approximate this N-pole operator => 'time-domain")
print("    IABC' is an N-pole Bessel-filter network (analysis: knowledge module).")
print("\nALL CHECKS PASSED.")
