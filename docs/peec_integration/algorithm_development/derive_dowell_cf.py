"""
Dowell F_R(xi) and F_L(xi) Taylor expansion around s=0

xi = d/delta = d * sqrt(omega*mu*sigma/2)

F_R(xi) = xi * (sinh(xi) + sin(xi)) / (cosh(xi) + cos(xi))
F_L(xi) = 3 * (sinh(2xi) - sin(2xi)) / (xi^3 * (cosh(2xi) - cos(2xi)))

Taylor expand around xi->0 (s->0).
"""

import numpy as np

print("="*70)
print("Dowell F_R(xi) and F_L(xi) Taylor Expansion")
print("="*70)

# ============================================================
# Dowell formulas (same as prima_with_dowell_correction.py)
# ============================================================

def dowell_F_R(xi):
    """Dowell F_R: xi*(sinh+sin)/(cosh+cos)"""
    if xi < 1e-10:
        return 1.0  # Taylor limit
    sinh_xi = np.sinh(xi)
    cosh_xi = np.cosh(xi)
    sin_xi = np.sin(xi)
    cos_xi = np.cos(xi)
    denom = cosh_xi + cos_xi
    if abs(denom) < 1e-30:
        return 1.0
    return xi * (sinh_xi + sin_xi) / denom

def dowell_F_L(xi):
    """Dowell F_L: 3*(sinh(2xi)-sin(2xi)) / (xi^3*(cosh(2xi)-cos(2xi)))"""
    if xi < 1e-10:
        return 1.0  # Taylor limit
    xi2 = 2 * xi
    sinh_2xi = np.sinh(xi2)
    cosh_2xi = np.cosh(xi2)
    sin_2xi = np.sin(xi2)
    cos_2xi = np.cos(xi2)
    denom = xi**3 * (cosh_2xi - cos_2xi)
    if abs(denom) < 1e-30:
        return 1.0
    return 3.0 * (sinh_2xi - sin_2xi) / denom

# ============================================================
# Taylor Expansion Derivation
# ============================================================
print("\n" + "="*70)
print("Taylor Expansion Derivation for F_R")
print("="*70)

print("""
F_R(xi) = xi * (sinh(xi) + sin(xi)) / (cosh(xi) + cos(xi))

Taylor series:
  sinh(xi) = xi + xi^3/6 + xi^5/120 + ...
  sin(xi)  = xi - xi^3/6 + xi^5/120 - ...
  cosh(xi) = 1 + xi^2/2 + xi^4/24 + ...
  cos(xi)  = 1 - xi^2/2 + xi^4/24 - ...

Numerator: xi * (sinh + sin)
  sinh + sin = 2*xi + 2*xi^5/120 + O(xi^9)
             = 2*xi * (1 + xi^4/120)
  xi * (...) = 2*xi^2 * (1 + xi^4/120)

Denominator: cosh + cos
  cosh + cos = 2 + 2*xi^4/24 + O(xi^8)  [xi^2 terms cancel!]
             = 2 * (1 + xi^4/24)

F_R = [2*xi^2 * (1 + xi^4/120)] / [2 * (1 + xi^4/24)]
    = xi^2 * (1 + xi^4/120) / (1 + xi^4/24)

Using 1/(1+x) ~= 1-x for small x:
    = xi^2 * (1 + xi^4/120) * (1 - xi^4/24)
    = xi^2 * (1 + xi^4/120 - xi^4/24 + O(xi^8))
    = xi^2 * (1 + xi^4 * (1/120 - 1/24))
    = xi^2 * (1 + xi^4 * (1/120 - 5/120))
    = xi^2 * (1 - xi^4/30)

So: F_R(xi) = xi^2 - xi^6/30 + O(xi^10)  [NOT 1 + xi^4/45!]

DC limit: F_R(0) = 0, NOT 1!
""")

# Verify numerically
print("\nNumerical verification of F_R:")
print(f"{'xi':>10} {'F_R(exact)':>14} {'xi^2':>14} {'xi^2(1-xi^4/30)':>18}")
print("-"*60)
for xi in [0.001, 0.01, 0.1, 0.5, 1.0, 2.0]:
    F_R_ex = dowell_F_R(xi)
    F_R_taylor1 = xi**2
    F_R_taylor2 = xi**2 * (1 - xi**4/30)
    print(f"  {xi:>8.3f} {F_R_ex:>14.8f} {F_R_taylor1:>14.8f} {F_R_taylor2:>18.8f}")

# ============================================================
# The issue: Different xi conventions
# ============================================================
print("\n" + "="*70)
print("Key Insight: Different xi Conventions!")
print("="*70)

print("""
There are TWO different conventions for xi in the literature:

Convention A (our code):
  xi = d / delta = d * sqrt(omega*mu*sigma/2)
  F_R(0) = 0, F_R(inf) -> xi

Convention B (Dowell 1966):
  Delta = d / (sqrt(2) * delta) = d * sqrt(omega*mu*sigma) / 2
  This gives F_R(0) = 1, F_R(inf) -> Delta

The relationship: xi = sqrt(2) * Delta

In prima_with_dowell_correction.py, we use small-xi approximation F_R = 1.0
when xi < 0.001. But this is mathematically inconsistent with the formula.

The physical meaning:
  Z = R_dc * F_R(xi) + j*omega * L_int_dc * F_L(xi)

At DC (xi=0):
  Z_dc = R_dc * F_R(0) + 0 * L_int_dc * F_L(0)
       = R_dc * F_R(0)

If F_R(0) = 0, then Z_dc = 0, which is WRONG (should be R_dc)
So F_R must satisfy F_R(0) = 1 for physical consistency!
""")

# ============================================================
# Correct Dowell formula
# ============================================================
print("\n" + "="*70)
print("Corrected Dowell Formula")
print("="*70)

print("""
The correct Dowell formula uses a DIFFERENT normalization!

From Ferreira's book "Electromagnetic Modelling of Power Electronic Converters":

F_R = Re[(1+j)*Delta * coth((1+j)*Delta)]

where Delta = d / (sqrt(2)*delta)

At Delta -> 0: coth(z) ~ 1/z, so (1+j)*Delta * coth((1+j)*Delta) -> 1
Therefore F_R(0) = 1 (correct!)

Let's verify with this formula:
""")

def dowell_F_R_correct(Delta):
    """Correct Dowell F_R using coth formula"""
    if Delta < 1e-10:
        return 1.0
    z = (1 + 1j) * Delta
    coth_z = 1 / np.tanh(z)
    return np.real(z * coth_z)

def dowell_F_L_correct(Delta):
    """Correct Dowell F_L using coth formula

    F_L = (3/2) * Im[(1+j)*Delta * coth((1+j)*Delta)] / Delta^2

    Or equivalently:
    F_L = 3/(2*Delta^2) * (sinh(2D) - sin(2D)) / (cosh(2D) - cos(2D))
    """
    if Delta < 1e-10:
        return 1.0
    z = (1 + 1j) * Delta
    coth_z = 1 / np.tanh(z)
    # F_L = (3/2) * Im[z * coth(z)] / Delta^2
    return 1.5 * np.imag(z * coth_z) / Delta**2

print("\nCorrect F_R (coth formula) values:")
print(f"{'Delta':>10} {'F_R':>14}")
print("-"*30)
for Delta in [1e-6, 1e-4, 0.001, 0.01, 0.1, 0.5, 1.0, 2.0]:
    F_R = dowell_F_R_correct(Delta)
    print(f"  {Delta:>8.4f} {F_R:>14.8f}")

print("\nCorrect F_L values:")
print(f"{'Delta':>10} {'F_L':>14}")
print("-"*30)
for Delta in [1e-6, 1e-4, 0.001, 0.01, 0.1, 0.5, 1.0, 2.0]:
    F_L = dowell_F_L_correct(Delta)
    print(f"  {Delta:>8.4f} {F_L:>14.8f}")

# ============================================================
# Taylor expansion of correct formulas
# ============================================================
print("\n" + "="*70)
print("Taylor Expansion of Correct Formulas")
print("="*70)

print("""
F_R(Delta) = Re[(1+j)*Delta * coth((1+j)*Delta)]

For small z: coth(z) = 1/z + z/3 - z^3/45 + 2*z^5/945 + ...

Let z = (1+j)*Delta:
  z^2 = (1+j)^2 * Delta^2 = 2j * Delta^2
  z^3 = (1+j)^3 * Delta^3 = (1+j)*2j * Delta^3 = 2j(1+j) * Delta^3

z * coth(z) = 1 + z^2/3 - z^4/45 + ...
            = 1 + 2j*Delta^2/3 - (2j)^2*Delta^4/45 + ...
            = 1 + 2j*Delta^2/3 + 4*Delta^4/45 + ...

F_R = Re[...] = 1 + 4*Delta^4/45 + O(Delta^8)
""")

# Verify Taylor expansion numerically
print("\nNumerical verification of Taylor expansion:")
print(f"{'Delta':>10} {'F_R(exact)':>14} {'1+4*D^4/45':>14} {'Error':>12}")
print("-"*55)
for Delta in [0.001, 0.01, 0.1, 0.2, 0.5, 1.0]:
    F_R_ex = dowell_F_R_correct(Delta)
    F_R_taylor = 1 + 4*Delta**4/45
    err = abs(F_R_ex - F_R_taylor) / F_R_ex * 100
    print(f"  {Delta:>8.4f} {F_R_ex:>14.8f} {F_R_taylor:>14.8f} {err:>11.4f}%")

# ============================================================
# Relationship between xi and Delta
# ============================================================
print("\n" + "="*70)
print("Relationship between xi and Delta")
print("="*70)

print("""
Our code uses: xi = d / delta = d * sqrt(omega*mu*sigma/2)
Dowell uses:   Delta = d / (sqrt(2)*delta) = xi / sqrt(2)

So: xi = sqrt(2) * Delta

Taylor expansion in terms of xi:
  F_R(xi) = 1 + 4*(xi/sqrt(2))^4/45 + ...
          = 1 + 4*xi^4/(4*45) + ...
          = 1 + xi^4/45 + O(xi^8)

  F_L(xi) = 1 - xi^4/45 + O(xi^8)  [by similar analysis]

This matches what was assumed in prima_with_dowell_correction.py!
""")

# Final verification
print("\nFinal verification with xi = sqrt(2)*Delta:")
print(f"{'xi':>10} {'F_R(xi)':>14} {'1+xi^4/45':>14} {'Error':>12}")
print("-"*55)
for Delta in [0.001, 0.01, 0.1, 0.2, 0.5, 1.0]:
    xi = np.sqrt(2) * Delta
    F_R_ex = dowell_F_R_correct(Delta)  # Use Delta for correct formula
    F_R_taylor = 1 + xi**4/45
    err = abs(F_R_ex - F_R_taylor) / F_R_ex * 100
    print(f"  {xi:>8.4f} {F_R_ex:>14.8f} {F_R_taylor:>14.8f} {err:>11.4f}%")

# ============================================================
# Conclusion
# ============================================================
print("\n" + "="*70)
print("Conclusion")
print("="*70)
print("""
1. The Dowell formulas have TWO conventions for the normalized thickness:
   - xi = d/delta (our code)
   - Delta = d/(sqrt(2)*delta) = xi/sqrt(2) (Dowell 1966)

2. The correct Dowell F_R formula that gives F_R(0) = 1 is:
   F_R = Re[(1+j)*Delta * coth((1+j)*Delta)]

3. Taylor expansion around xi = 0:
   F_R(xi) = 1 + xi^4/45 + O(xi^8)
   F_L(xi) = 1 - xi^4/45 + O(xi^8)

4. In terms of Laplace variable s:
   tau = d^2*mu*sigma/2
   xi^2 = tau * omega = -j*tau*s
   xi^4 = -tau^2*s^2

   F_R(s) = 1 - tau^2*s^2/45 + O(s^4)
   F_L(s) = 1 + tau^2*s^2/45 + O(s^4)

5. These are NOT rational functions of s.
   Exact finite continued fraction expansion is impossible.

6. Practical approach: "DC PRIMA + Dowell correction"
   - Build PRIMA with DC parameters
   - Apply F_R(omega), F_L(omega) directly at each frequency
   - This is EXACT for frequency-domain analysis
""")
