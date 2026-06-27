"""
Dowell F_R(s) and F_L(s) - Continued Fraction via Viskovatov Algorithm

Viskovatov (or quotient-difference) algorithm converts a power series
to a continued fraction.

Given: f(x) = c0 + c1*x + c2*x^2 + c3*x^3 + ...
Find:  f(x) = b0 + a1*x / (1 + a2*x / (1 + a3*x / (1 + ...)))
"""

import numpy as np
import matplotlib.pyplot as plt

MU_0 = 4 * np.pi * 1e-7

print("="*70)
print("Viskovatov Algorithm for Continued Fraction Expansion")
print("="*70)

# ============================================================
# Exact Dowell formulas
# ============================================================

def dowell_F_R_exact(Delta):
    """Exact Dowell F_R"""
    if Delta < 1e-10:
        return 1.0
    z = (1 + 1j) * Delta
    coth_z = 1 / np.tanh(z)
    return np.real(z * coth_z)

def dowell_F_L_exact(Delta):
    """Exact Dowell F_L"""
    if Delta < 1e-10:
        return 1.0
    z = (1 + 1j) * Delta
    coth_z = 1 / np.tanh(z)
    return 1.5 * np.imag(z * coth_z) / Delta**2

# ============================================================
# Taylor coefficients of z*coth(z)
# ============================================================
print("\n" + "="*70)
print("Taylor Series of z*coth(z)")
print("="*70)

print("""
z*coth(z) = 1 + z^2/3 - z^4/45 + 2*z^6/945 - z^8/4725 + 2*z^10/93555 - ...

The coefficients are related to Bernoulli numbers:
  z*coth(z) = sum_{n=0}^inf (2^(2n) * B_{2n} / (2n)!) * z^(2n)

where B_0 = 1, B_2 = 1/6, B_4 = -1/30, B_6 = 1/42, ...
""")

# Taylor coefficients of z*coth(z) for z^0, z^2, z^4, z^6, ...
# These are: 1, 1/3, -1/45, 2/945, -1/4725, 2/93555, ...
taylor_zcoth = [
    1.0,           # z^0
    1/3,           # z^2
    -1/45,         # z^4
    2/945,         # z^6
    -1/4725,       # z^8
    2/93555,       # z^10
    -1382/638512875,  # z^12
]

print("\nTaylor coefficients of z*coth(z):")
for i, c in enumerate(taylor_zcoth):
    print(f"  z^{2*i}: {c:.10f}")

# ============================================================
# Viskovatov Algorithm
# ============================================================
print("\n" + "="*70)
print("Viskovatov Algorithm")
print("="*70)

print("""
For a power series: f(x) = c0 + c1*x + c2*x^2 + ...

The algorithm computes continued fraction:
  f(x) = b0 + a1*x / (1 + a2*x / (1 + a3*x / ...))

Step 1: b0 = c0, then work with (f(x) - c0) / x
Step 2: Invert: x / (f(x) - c0) = 1 / (c1 + c2*x + ...)
Step 3: Repeat...

For z*coth(z), we have a series in z^2:
  z*coth(z) = 1 + (1/3)*w - (1/45)*w^2 + (2/945)*w^3 - ...
  where w = z^2

So we expand in w and get a continued fraction in w.
""")

def viskovatov_cf(coeffs, n_terms):
    """
    Viskovatov algorithm to convert power series to continued fraction.

    Input: coeffs = [c0, c1, c2, ...] for f(x) = c0 + c1*x + c2*x^2 + ...
    Output: (b0, [a1, a2, a3, ...]) for f(x) = b0 + a1*x/(1 + a2*x/(1 + ...))

    Note: This computes the S-fraction form.
    """
    # Work with coefficient arrays
    c = np.array(coeffs, dtype=float)
    n = len(c)

    # Result arrays
    b0 = c[0]
    a_list = []

    # Remove constant term and divide by x
    # f(x) - c0 = c1*x + c2*x^2 + ...
    # (f(x) - c0)/x = c1 + c2*x + c3*x^2 + ...
    c = c[1:]  # Now c = [c1, c2, c3, ...]

    for _ in range(n_terms):
        if len(c) == 0 or abs(c[0]) < 1e-15:
            break

        # a_k = 1/c[0]
        a_k = 1.0 / c[0]
        a_list.append(a_k)

        # Invert the series: 1/(c0 + c1*x + c2*x^2 + ...)
        # Result: d0 + d1*x + d2*x^2 + ... where d0 = 1/c0
        # This is series inversion
        c_inv = invert_series(c)

        # Subtract 1/c[0] and divide by x
        # (1/(c0 + c1*x + ...) - 1/c0) / x = ...
        c_inv[0] = 0  # Remove constant term
        c = c_inv[1:]  # Divide by x

    return b0, a_list

def invert_series(c, max_terms=20):
    """
    Invert power series: given f(x) = c0 + c1*x + c2*x^2 + ...
    Find g(x) = d0 + d1*x + d2*x^2 + ... such that f(x)*g(x) = 1
    """
    n = min(len(c), max_terms)
    d = np.zeros(n)

    if abs(c[0]) < 1e-15:
        return d  # Cannot invert if c0 = 0

    d[0] = 1.0 / c[0]

    for k in range(1, n):
        s = 0.0
        for j in range(1, min(k+1, len(c))):
            if j < len(c):
                s += c[j] * d[k-j]
        d[k] = -s / c[0]

    return d

# Apply to z*coth(z)
print("\nApplying Viskovatov to z*coth(z) (in w = z^2):")
b0, a_list = viskovatov_cf(taylor_zcoth, 10)
print(f"  b0 = {b0}")
for i, a in enumerate(a_list):
    print(f"  a{i+1} = {a:.10f}")

# ============================================================
# Known analytical result
# ============================================================
print("\n" + "="*70)
print("Analytical Continued Fraction of z*coth(z)")
print("="*70)

print("""
The well-known continued fraction for z*coth(z) is:

z*coth(z) = 1 + z^2/(3 + z^2/(5 + z^2/(7 + z^2/(9 + ...))))

This is a J-fraction (Jacobi continued fraction).

In the form: b0 + w/(b1 + w/(b2 + w/(b3 + ...)))
where w = z^2:
  b0 = 1, b1 = 3, b2 = 5, b3 = 7, b4 = 9, ...
  b_n = 2n + 1

This is NOT the same form as Viskovatov's S-fraction!
""")

def zcoth_jfrac(w, n_terms=10):
    """J-fraction: z*coth(z) = 1 + w/(3 + w/(5 + w/(7 + ...)))"""
    result = 2*n_terms + 1  # Start from innermost
    for k in range(n_terms-1, 0, -1):
        result = (2*k + 1) + w / result
    return 1 + w / result

print("\nVerification of J-fraction:")
print(f"{'w':>10} {'Exact':>14} {'J-frac(5)':>14} {'J-frac(10)':>14}")
print("-"*55)
for w_val in [0.01, 0.1, 1.0, 4.0, 10.0]:
    z = np.sqrt(w_val)
    z_complex = z * (1 + 1j) / np.sqrt(2)  # z = (1+j)*Delta, w = z^2 = 2j*Delta^2
    exact = z_complex * (1 / np.tanh(z_complex))
    jf5 = zcoth_jfrac(w_val, 5)
    jf10 = zcoth_jfrac(w_val, 10)
    print(f"  {w_val:>8.2f} {np.abs(exact):>14.8f} {np.abs(jf5):>14.8f} {np.abs(jf10):>14.8f}")

# ============================================================
# Convert J-fraction to PRIMA form
# ============================================================
print("\n" + "="*70)
print("J-fraction to PRIMA Ladder Network")
print("="*70)

print("""
J-fraction: f(w) = 1 + w/(3 + w/(5 + w/(7 + ...)))

where w = tau*s = d^2*mu*sigma*s/2

For impedance Z(s) = R_dc * F_R + s*L_dc * F_L:

At s = j*omega:
  w = j*tau*omega

The impedance can be written as a ladder network:
  Z(s) = R_dc + s*L_dc * f(tau*s)

The J-fraction in w naturally corresponds to:
  Stage 1: Series R, Shunt L
  Stage 2: Series R, Shunt L
  ...

But the exact correspondence depends on the impedance formula structure.
""")

# ============================================================
# PRIMA parameters from continued fraction
# ============================================================
print("\n" + "="*70)
print("PRIMA I Parameters (from 1D diffusion eigenvalue expansion)")
print("="*70)

print("""
The PRIMA I form comes from eigenvalue expansion of 1D diffusion:

R[n] = (4n-5) * 4 / (sigma*d)  for n >= 2, R[1] ~ 0
L[n] = d*mu / (4n-3)

The continued fraction for PRIMA I is different from z*coth(z).

PRIMA I represents:
  Z_prima(s) = sL[1] + R[1] || (sL[2] + R[2] || (sL[3] + R[3] || ...))

which matches the 1D diffusion equation solution:
  Z_diff(s) = s*mu*(2/(k*d))*tan(k*d/2)*d

NOT Dowell's formula: Z_dowell = R_dc*F_R + s*L_dc*F_L

So:
- PRIMA I matches 1D Diffusion (exact)
- PRIMA I does NOT match Dowell (different model)
- Dowell has its own continued fraction (J-fraction)
""")

def skin_effect_prima_params(d, sigma, mu=MU_0, n_stages=7):
    """PRIMA I parameters for 1D diffusion skin effect"""
    R_prima = np.zeros(n_stages)
    L_prima = np.zeros(n_stages)

    for n in range(1, n_stages + 1):
        if n == 1:
            R_prima[n-1] = 1e-16  # Nearly zero
        else:
            R_prima[n-1] = (4*n - 5) * 4.0 / (sigma * d)
        L_prima[n-1] = d * mu / (4*n - 3)

    return R_prima, L_prima

# Example parameters
d = 0.1e-3  # 0.1 mm
sigma = 5.8e7  # Copper

print(f"\nPRIMA I parameters for d = {d*1e3} mm, sigma = {sigma:.1e} S/m:")
R_prima, L_prima = skin_effect_prima_params(d, sigma, MU_0, 7)
for i in range(7):
    print(f"  Stage {i+1}: R = {R_prima[i]:.4e} Ohm*m^2, L = {L_prima[i]*1e9:.4f} nH*m^2")

# ============================================================
# Dowell PRIMA form
# ============================================================
print("\n" + "="*70)
print("Dowell-based PRIMA (Different from 1D Diffusion PRIMA)")
print("="*70)

print("""
Dowell's formula:
  Z(s) = R_dc * F_R(w) + s * L_dc * F_L(w)
  where w = tau*s, tau = d^2*mu*sigma/2

From J-fraction:
  z*coth(z) = 1 + w/(3 + w/(5 + w/(7 + ...)))
  where z^2 = w*2j

We can write F_R and F_L in terms of continued fractions.

However, the structure is:
  Z = R_dc * (rational function of tau*s)

This is different from the PRIMA I structure:
  Z = s*L + parallel RC network

The "DC PRIMA + Dowell correction" approach avoids this mismatch:
1. Build PRIMA with DC parameters
2. Apply F_R(omega), F_L(omega) at each frequency
3. No need for continued fraction conversion
""")

# ============================================================
# Conclusion
# ============================================================
print("\n" + "="*70)
print("Conclusion")
print("="*70)
print("""
1. z*coth(z) has a well-known J-fraction (Jacobi continued fraction):
   z*coth(z) = 1 + z^2/(3 + z^2/(5 + z^2/(7 + ...)))

2. This can also be written in w = z^2 = tau*s:
   1 + w/(3 + w/(5 + w/(7 + ...)))

   Coefficients: b_n = 2n + 1 (i.e., 3, 5, 7, 9, ...)

3. The Viskovatov algorithm can convert Taylor series to S-fraction,
   but the J-fraction form is more natural for z*coth(z).

4. The PRIMA I form (from 1D diffusion) has DIFFERENT structure:
   - PRIMA I matches 1D diffusion equation
   - PRIMA I does NOT match Dowell formula

5. For Dowell-based PEEC, use "DC PRIMA + Dowell correction":
   - No continued fraction conversion needed
   - Apply F_R(omega), F_L(omega) directly per frequency
   - This is EXACT (no approximation error)
""")

# ============================================================
# Plot
# ============================================================
print("\n" + "="*70)
print("Generating comparison plot...")
print("="*70)

Delta_range = np.linspace(0.01, 3, 200)
w_range = 2 * Delta_range**2  # w = 2*Delta^2 for z = (1+j)*Delta

F_R_exact = np.array([dowell_F_R_exact(D) for D in Delta_range])
F_L_exact = np.array([dowell_F_L_exact(D) for D in Delta_range])

# J-fraction approximations
def F_R_jfrac(Delta, n):
    z = (1 + 1j) * Delta
    w = z**2  # = 2j*Delta^2
    return np.real(zcoth_jfrac(w, n))

def F_L_jfrac(Delta, n):
    if Delta < 1e-10:
        return 1.0
    z = (1 + 1j) * Delta
    w = z**2
    return 1.5 * np.imag(zcoth_jfrac(w, n)) / Delta**2

F_R_jf3 = np.array([F_R_jfrac(D, 3) for D in Delta_range])
F_R_jf5 = np.array([F_R_jfrac(D, 5) for D in Delta_range])
F_R_jf10 = np.array([F_R_jfrac(D, 10) for D in Delta_range])
F_L_jf3 = np.array([F_L_jfrac(D, 3) for D in Delta_range])
F_L_jf5 = np.array([F_L_jfrac(D, 5) for D in Delta_range])
F_L_jf10 = np.array([F_L_jfrac(D, 10) for D in Delta_range])

fig, axes = plt.subplots(2, 2, figsize=(12, 10))

ax1 = axes[0, 0]
ax1.plot(Delta_range, F_R_exact, 'k-', linewidth=2, label='F_R exact')
ax1.plot(Delta_range, F_R_jf3, 'b--', linewidth=1.5, label='J-frac (3 terms)')
ax1.plot(Delta_range, F_R_jf5, 'g--', linewidth=1.5, label='J-frac (5 terms)')
ax1.plot(Delta_range, F_R_jf10, 'r--', linewidth=1.5, label='J-frac (10 terms)')
ax1.set_xlabel('Delta = d/(sqrt(2)*delta)')
ax1.set_ylabel('F_R')
ax1.set_title('F_R: Exact vs J-fraction')
ax1.legend()
ax1.grid(True, alpha=0.3)

ax2 = axes[0, 1]
ax2.plot(Delta_range, F_L_exact, 'k-', linewidth=2, label='F_L exact')
ax2.plot(Delta_range, F_L_jf3, 'b--', linewidth=1.5, label='J-frac (3 terms)')
ax2.plot(Delta_range, F_L_jf5, 'g--', linewidth=1.5, label='J-frac (5 terms)')
ax2.plot(Delta_range, F_L_jf10, 'r--', linewidth=1.5, label='J-frac (10 terms)')
ax2.set_xlabel('Delta = d/(sqrt(2)*delta)')
ax2.set_ylabel('F_L')
ax2.set_title('F_L: Exact vs J-fraction')
ax2.legend()
ax2.grid(True, alpha=0.3)

# Error plots
ax3 = axes[1, 0]
err_jf3 = np.abs(F_R_jf3 - F_R_exact) / F_R_exact * 100
err_jf5 = np.abs(F_R_jf5 - F_R_exact) / F_R_exact * 100
err_jf10 = np.abs(F_R_jf10 - F_R_exact) / F_R_exact * 100
ax3.semilogy(Delta_range, err_jf3, 'b-', linewidth=1.5, label='J-frac (3 terms)')
ax3.semilogy(Delta_range, err_jf5, 'g-', linewidth=1.5, label='J-frac (5 terms)')
ax3.semilogy(Delta_range, err_jf10, 'r-', linewidth=1.5, label='J-frac (10 terms)')
ax3.set_xlabel('Delta')
ax3.set_ylabel('F_R Error [%]')
ax3.set_title('F_R Approximation Error')
ax3.legend()
ax3.grid(True, alpha=0.3)

ax4 = axes[1, 1]
err_jf3_L = np.abs(F_L_jf3 - F_L_exact) / np.maximum(F_L_exact, 1e-10) * 100
err_jf5_L = np.abs(F_L_jf5 - F_L_exact) / np.maximum(F_L_exact, 1e-10) * 100
err_jf10_L = np.abs(F_L_jf10 - F_L_exact) / np.maximum(F_L_exact, 1e-10) * 100
ax4.semilogy(Delta_range, err_jf3_L, 'b-', linewidth=1.5, label='J-frac (3 terms)')
ax4.semilogy(Delta_range, err_jf5_L, 'g-', linewidth=1.5, label='J-frac (5 terms)')
ax4.semilogy(Delta_range, err_jf10_L, 'r-', linewidth=1.5, label='J-frac (10 terms)')
ax4.set_xlabel('Delta')
ax4.set_ylabel('F_L Error [%]')
ax4.set_title('F_L Approximation Error')
ax4.legend()
ax4.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('derive_dowell_cf_algorithm.png', dpi=150)
print("Saved: derive_dowell_cf_algorithm.png")
plt.close()
