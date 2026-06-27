"""
Dowell F_R(s) and F_L(s) Pade (Rational Function) Expansion

F_R and F_L are transcendental functions of xi = sqrt(-j*tau*s).
We approximate them as rational functions of s for use in PRIMA.

Pade approximation: P(s)/Q(s) where P, Q are polynomials
"""

import numpy as np
import matplotlib.pyplot as plt

MU_0 = 4 * np.pi * 1e-7

print("="*70)
print("Dowell F_R(s) and F_L(s) Pade (Rational Function) Expansion")
print("="*70)

# ============================================================
# Exact Dowell formulas (coth form)
# ============================================================

def dowell_F_R_exact(Delta):
    """Exact Dowell F_R: Re[(1+j)*Delta * coth((1+j)*Delta)]"""
    if Delta < 1e-10:
        return 1.0
    z = (1 + 1j) * Delta
    coth_z = 1 / np.tanh(z)
    return np.real(z * coth_z)

def dowell_F_L_exact(Delta):
    """Exact Dowell F_L: (3/2) * Im[(1+j)*Delta * coth((1+j)*Delta)] / Delta^2"""
    if Delta < 1e-10:
        return 1.0
    z = (1 + 1j) * Delta
    coth_z = 1 / np.tanh(z)
    return 1.5 * np.imag(z * coth_z) / Delta**2

# ============================================================
# Pade approximation derivation
# ============================================================
print("\n" + "="*70)
print("Pade Approximation Derivation")
print("="*70)

print("""
The key function is: z * coth(z) where z = (1+j)*Delta

coth(z) has a well-known continued fraction expansion:
  coth(z) = 1/z + z/3 - z^3/45 + 2*z^5/945 - ...

Therefore:
  z * coth(z) = 1 + z^2/3 - z^4/45 + 2*z^6/945 - ...

With z = (1+j)*Delta and z^2 = 2j*Delta^2:
  z * coth(z) = 1 + (2j*Delta^2)/3 - (2j)^2*Delta^4/45 + ...
              = 1 + (2j/3)*Delta^2 + (4/45)*Delta^4 + ...

F_R = Re[z*coth(z)] = 1 + (4/45)*Delta^4 + O(Delta^8)
F_L = (3/2)*Im[z*coth(z)]/Delta^2 = (3/2)*(2/3) + O(Delta^4) = 1 + O(Delta^4)

For rational approximation, we use the continued fraction form of coth(z):

coth(z) = 1/z + 1/(3z + 1/(5z + 1/(7z + ...)))

This is a Stieltjes continued fraction!
""")

# ============================================================
# Continued fraction of coth(z)
# ============================================================
print("\n" + "="*70)
print("Continued Fraction of coth(z)")
print("="*70)

print("""
coth(z) has the Stieltjes continued fraction:

coth(z) = 1/z + 1/(3/z + 1/(5/z + 1/(7/z + ...)))

Or equivalently:
  z * coth(z) = 1 + z^2/(3 + z^2/(5 + z^2/(7 + ...)))

This is a K-type continued fraction.

Truncating at n terms gives a rational [n,n] Pade approximant.
""")

def coth_cf(z, n_terms=5):
    """Continued fraction approximation of coth(z)"""
    # coth(z) = 1/z + 1/(3/z + 1/(5/z + ...))
    # We compute z*coth(z) = 1 + z^2/(3 + z^2/(5 + z^2/(7 + ...)))
    z2 = z * z
    result = 2 * n_terms + 1  # Starting from innermost: (2n+1)
    for k in range(n_terms - 1, 0, -1):
        result = (2 * k + 1) + z2 / result
    # z*coth(z) = 1 + z^2/result
    return 1 + z2 / result

print("\nVerification of continued fraction for z*coth(z):")
print(f"{'z':>12} {'Exact':>14} {'CF(5 terms)':>14} {'CF(10 terms)':>14}")
print("-"*60)
for z_test in [0.1, 0.5, 1.0, 2.0, 3.0]:
    z = complex(z_test, z_test)  # z = (1+j)*Delta/sqrt(2)
    exact = z * (1 / np.tanh(z))
    cf5 = coth_cf(z, 5)
    cf10 = coth_cf(z, 10)
    print(f"  {z_test:>10.2f} {np.abs(exact):>14.8f} {np.abs(cf5):>14.8f} {np.abs(cf10):>14.8f}")

# ============================================================
# F_R and F_L as rational functions
# ============================================================
print("\n" + "="*70)
print("F_R and F_L as Rational Functions of Delta^2")
print("="*70)

print("""
Let u = Delta^2.  Then z^2 = (1+j)^2 * Delta^2 = 2j*u

z*coth(z) = 1 + 2j*u/(3 + 2j*u/(5 + 2j*u/(7 + ...)))

This is a continued fraction in the complex variable 2j*u.

Truncating at N terms:
  z*coth(z) ~= P_N(2j*u) / Q_N(2j*u)

where P_N, Q_N are polynomials of degree N in u.

F_R = Re[P_N/Q_N], F_L = (3/2)*Im[P_N/Q_N]/u

Since we want rational functions of s (not sqrt(s)):
  Delta^2 = xi^2/2 = tau*omega/2 = -j*tau*s/2

  u = Delta^2 = -j*tau*s/2
  2j*u = 2j*(-j*tau*s/2) = tau*s

So: z*coth(z) = f(tau*s) where f is a rational function!
""")

# ============================================================
# Explicit Pade coefficients
# ============================================================
print("\n" + "="*70)
print("Explicit Pade Coefficients")
print("="*70)

print("""
From the continued fraction:
  z*coth(z) = 1 + w/(3 + w/(5 + w/(7 + ...)))

where w = 2j*Delta^2 = tau*s

N=1 (Pade [1,1]):
  z*coth(z) ~= 1 + w/3 = (3 + w)/3 = (3 + tau*s)/3

N=2 (Pade [2,2]):
  z*coth(z) ~= 1 + w/(3 + w/5) = 1 + 5w/(15 + w)
             = (15 + w + 5w)/(15 + w) = (15 + 6w)/(15 + w)
             = (15 + 6*tau*s)/(15 + tau*s)

N=3 (Pade [3,3]):
  Inner: 5 + w/7 = (35 + w)/7
  Middle: 3 + w/((35+w)/7) = 3 + 7w/(35+w) = (105 + 3w + 7w)/(35+w) = (105 + 10w)/(35+w)
  Outer: 1 + w/((105+10w)/(35+w)) = 1 + w*(35+w)/(105+10w)
       = (105 + 10w + 35w + w^2)/(105 + 10w)
       = (105 + 45w + w^2)/(105 + 10w)
""")

def zcoth_pade(w, order=3):
    """Pade approximation of z*coth(z) where w = 2j*Delta^2 = tau*s"""
    if order == 1:
        return (3 + w) / 3
    elif order == 2:
        return (15 + 6*w) / (15 + w)
    elif order == 3:
        return (105 + 45*w + w**2) / (105 + 10*w)
    elif order == 4:
        # Continue the pattern...
        # Inner: 7 + w/9 = (63+w)/9
        # ...
        # This gets complex, use numerical continued fraction instead
        return coth_cf(np.sqrt(w/2) * (1+1j), order)
    else:
        return coth_cf(np.sqrt(w/2) * (1+1j), order)

print("\nVerification of Pade approximants:")
print(f"{'Delta':>8} {'Exact':>12} {'Pade[1,1]':>12} {'Pade[2,2]':>12} {'Pade[3,3]':>12}")
print("-"*60)
for Delta in [0.1, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0]:
    z = (1 + 1j) * Delta
    w = 2j * Delta**2
    exact = z * (1 / np.tanh(z))
    p1 = zcoth_pade(w, 1)
    p2 = zcoth_pade(w, 2)
    p3 = zcoth_pade(w, 3)
    print(f"  {Delta:>6.2f} {np.abs(exact):>12.6f} {np.abs(p1):>12.6f} {np.abs(p2):>12.6f} {np.abs(p3):>12.6f}")

# ============================================================
# F_R and F_L from Pade
# ============================================================
print("\n" + "="*70)
print("F_R and F_L from Pade Approximation")
print("="*70)

def F_R_pade(Delta, order=3):
    """F_R from Pade approximation"""
    if Delta < 1e-10:
        return 1.0
    w = 2j * Delta**2
    zcoth = zcoth_pade(w, order)
    return np.real(zcoth)

def F_L_pade(Delta, order=3):
    """F_L from Pade approximation"""
    if Delta < 1e-10:
        return 1.0
    w = 2j * Delta**2
    zcoth = zcoth_pade(w, order)
    return 1.5 * np.imag(zcoth) / Delta**2

print("\nF_R comparison:")
print(f"{'Delta':>8} {'Exact':>12} {'Pade[1,1]':>12} {'Pade[2,2]':>12} {'Pade[3,3]':>12} {'Err[3,3]':>10}")
print("-"*70)
for Delta in [0.1, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0]:
    exact = dowell_F_R_exact(Delta)
    p1 = F_R_pade(Delta, 1)
    p2 = F_R_pade(Delta, 2)
    p3 = F_R_pade(Delta, 3)
    err = abs(exact - p3) / exact * 100
    print(f"  {Delta:>6.2f} {exact:>12.6f} {p1:>12.6f} {p2:>12.6f} {p3:>12.6f} {err:>9.4f}%")

print("\nF_L comparison:")
print(f"{'Delta':>8} {'Exact':>12} {'Pade[1,1]':>12} {'Pade[2,2]':>12} {'Pade[3,3]':>12} {'Err[3,3]':>10}")
print("-"*70)
for Delta in [0.1, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0]:
    exact = dowell_F_L_exact(Delta)
    p1 = F_L_pade(Delta, 1)
    p2 = F_L_pade(Delta, 2)
    p3 = F_L_pade(Delta, 3)
    err = abs(exact - p3) / max(exact, 1e-10) * 100
    print(f"  {Delta:>6.2f} {exact:>12.6f} {p1:>12.6f} {p2:>12.6f} {p3:>12.6f} {err:>9.4f}%")

# ============================================================
# Conversion to s-domain
# ============================================================
print("\n" + "="*70)
print("Conversion to s-domain (Laplace Variable)")
print("="*70)

print("""
Key relationship:
  w = 2j*Delta^2 = 2j*(xi/sqrt(2))^2 = j*xi^2 = j*tau*omega = tau*s

where tau = d^2*mu*sigma/2

So w = tau*s, and:

z*coth(z) = P(tau*s) / Q(tau*s)  (rational function of s!)

For Pade[3,3]:
  z*coth(z) = (105 + 45*tau*s + (tau*s)^2) / (105 + 10*tau*s)

F_R = Re[z*coth(z)]
F_L = (3/2) * Im[z*coth(z)] / Delta^2
    = (3/2) * Im[z*coth(z)] / (tau*omega/2)
    = 3 * Im[z*coth(z)] / (tau*omega)

At s = j*omega:
  tau*s = j*tau*omega

For real omega > 0:
  w = tau*s = j*tau*omega (purely imaginary)

Let's denote x = tau*omega (real):
  w = j*x

Pade[3,3]:
  z*coth(z) = (105 + 45*j*x - x^2) / (105 + 10*j*x)
            = (105 - x^2 + j*45*x) / (105 + j*10*x)

Multiply by conjugate:
  = [(105-x^2 + j*45*x)(105 - j*10*x)] / [(105)^2 + (10*x)^2]
  = [(105-x^2)*105 + 45*10*x^2 + j*(45*105*x - (105-x^2)*10*x)] / [105^2 + 100*x^2]
  = [105*(105-x^2) + 450*x^2 + j*x*(45*105 - 10*(105-x^2))] / [11025 + 100*x^2]
  = [11025 - 105*x^2 + 450*x^2 + j*x*(4725 - 1050 + 10*x^2)] / [11025 + 100*x^2]
  = [11025 + 345*x^2 + j*x*(3675 + 10*x^2)] / [11025 + 100*x^2]

F_R = [11025 + 345*x^2] / [11025 + 100*x^2]  (x = tau*omega)
F_L = 3 * x*(3675 + 10*x^2) / [(11025 + 100*x^2) * tau*omega]
    = 3 * (3675 + 10*x^2) / (11025 + 100*x^2)
    = (11025 + 30*x^2) / (11025 + 100*x^2)
""")

def F_R_rational(x):
    """F_R as rational function of x = tau*omega (Pade[3,3])"""
    x2 = x**2
    return (11025 + 345*x2) / (11025 + 100*x2)

def F_L_rational(x):
    """F_L as rational function of x = tau*omega (Pade[3,3])"""
    x2 = x**2
    return (11025 + 30*x2) / (11025 + 100*x2)

print("\nVerification of rational form (Pade[3,3]):")
print(f"{'x=tau*w':>10} {'F_R(exact)':>12} {'F_R(rat)':>12} {'F_L(exact)':>12} {'F_L(rat)':>12}")
print("-"*65)
for x in [0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0]:
    # Delta^2 = tau*omega/2 = x/2
    Delta = np.sqrt(x/2)
    F_R_ex = dowell_F_R_exact(Delta)
    F_L_ex = dowell_F_L_exact(Delta)
    F_R_rat = F_R_rational(x)
    F_L_rat = F_L_rational(x)
    print(f"  {x:>8.2f} {F_R_ex:>12.6f} {F_R_rat:>12.6f} {F_L_ex:>12.6f} {F_L_rat:>12.6f}")

# ============================================================
# Continued Fraction Form (for PRIMA)
# ============================================================
print("\n" + "="*70)
print("Continued Fraction Form for PRIMA")
print("="*70)

print("""
The rational functions can be converted to continued fractions:

F_R(x) = (11025 + 345*x^2) / (11025 + 100*x^2)
       = (a0 + a1*x^2) / (b0 + b1*x^2)

This is a [2,2] rational function in x^2, which corresponds to a
2-stage continued fraction.

For a general [n,n] Pade from z*coth(z):
  z*coth(z) = 1 + w/(3 + w/(5 + w/(7 + ...)))

This IS already in continued fraction form in w = tau*s!

For PRIMA implementation with N stages:
  Z(s) = R_dc * F_R(tau*s) + s * L_dc * F_L(tau*s)

Using Pade[N,N] approximations gives a (2N)-order rational function
that can be realized as a ladder network.
""")

# ============================================================
# Plot comparison
# ============================================================
print("\n" + "="*70)
print("Generating comparison plot...")
print("="*70)

x_range = np.logspace(-2, 2, 200)

F_R_exact_arr = np.array([dowell_F_R_exact(np.sqrt(x/2)) for x in x_range])
F_L_exact_arr = np.array([dowell_F_L_exact(np.sqrt(x/2)) for x in x_range])
F_R_pade1_arr = np.array([F_R_pade(np.sqrt(x/2), 1) for x in x_range])
F_R_pade2_arr = np.array([F_R_pade(np.sqrt(x/2), 2) for x in x_range])
F_R_pade3_arr = np.array([F_R_pade(np.sqrt(x/2), 3) for x in x_range])
F_R_rat_arr = F_R_rational(x_range)
F_L_rat_arr = F_L_rational(x_range)

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

ax1 = axes[0]
ax1.semilogx(x_range, F_R_exact_arr, 'k-', linewidth=2, label='F_R exact')
ax1.semilogx(x_range, F_R_pade1_arr, 'b--', linewidth=1.5, label='Pade[1,1]')
ax1.semilogx(x_range, F_R_pade2_arr, 'g--', linewidth=1.5, label='Pade[2,2]')
ax1.semilogx(x_range, F_R_pade3_arr, 'r--', linewidth=1.5, label='Pade[3,3]')
ax1.set_xlabel('x = tau * omega')
ax1.set_ylabel('F_R')
ax1.set_title('F_R: Exact vs Pade Approximations')
ax1.legend()
ax1.grid(True, which='both', alpha=0.3)
ax1.set_ylim([0.9, 3])

ax2 = axes[1]
ax2.semilogx(x_range, F_L_exact_arr, 'k-', linewidth=2, label='F_L exact')
ax2.semilogx(x_range, F_L_rat_arr, 'r--', linewidth=1.5, label='Pade[3,3] rational')
ax2.set_xlabel('x = tau * omega')
ax2.set_ylabel('F_L')
ax2.set_title('F_L: Exact vs Pade Approximation')
ax2.legend()
ax2.grid(True, which='both', alpha=0.3)
ax2.set_ylim([0, 1.2])

plt.tight_layout()
plt.savefig('derive_dowell_pade.png', dpi=150)
print("Saved: derive_dowell_pade.png")
plt.close()

# ============================================================
# Conclusion
# ============================================================
print("\n" + "="*70)
print("Conclusion")
print("="*70)
print("""
1. z*coth(z) has a natural continued fraction expansion:
   z*coth(z) = 1 + w/(3 + w/(5 + w/(7 + ...)))
   where w = tau*s, tau = d^2*mu*sigma/2

2. Truncating at N terms gives Pade[N,N] approximants:
   - Pade[1,1]: (3 + tau*s) / 3
   - Pade[2,2]: (15 + 6*tau*s) / (15 + tau*s)
   - Pade[3,3]: (105 + 45*tau*s + (tau*s)^2) / (105 + 10*tau*s)

3. F_R and F_L as rational functions of x = tau*omega (Pade[3,3]):
   F_R(x) = (11025 + 345*x^2) / (11025 + 100*x^2)
   F_L(x) = (11025 + 30*x^2) / (11025 + 100*x^2)

4. These rational functions CAN be converted to continued fractions
   and implemented as PRIMA ladder networks.

5. Higher-order Pade gives better accuracy at high frequencies,
   but Pade[3,3] is already very accurate for most applications.

6. The continued fraction form is natural because coth(z) itself
   has a Stieltjes continued fraction expansion.
""")
