# Radia PEEC Solver: Surface Impedance and SPICE Export

## Overview

This document describes Radia's integrated electromagnetic solver framework combining:

- **PEEC (Loop-Star)**: Conductor current distribution, inductive/capacitive coupling
- **HDiv-VIM / reduced FEM**: Magnetic-material response
- **Surface Impedance Zs**: Skin effect via Dowell, ESIM, or PyKAN learning

```
+-------------------------------------------------------------+
|                    Radia Integrated Solver                   |
+-------------------------------------------------------------+
|  +-----------+   +-----------+   +-----------+              |
|  |   PEEC    |   | HDiv-VIM |   | Surface   |              |
|  |(Loop-Star)| + | (Magnetic)| + |    Zs     |              |
|  +-----------+   +-----------+   +-----------+              |
|       |               |               |                     |
|   L, R, C         M(r), B(r)      R(f), L(f)               |
|   matrices        distribution    freq-dependent            |
+-------------------------------------------------------------+
|                    SPICE Output (Verilog-A)                 |
+-------------------------------------------------------------+
```

**Note**: Surface impedance is also **field-dependent**: Zs(H) = Re{Zs(H)} + j·Im{Zs(H)}

**Target Applications**: Power electronics (DC - 10 MHz)
- Wireless Power Transfer (WPT)
- Induction Heating
- High-Frequency Transformers
- EMI/EMC Filters

### Framework Strengths

| Component | Role | Key Feature |
|-----------|------|-------------|
| **PEEC Loop** | Inductive coupling | Handles μ(ω) via L matrix |
| **PEEC Star** | Capacitive coupling | Handles ε(ω) via C matrix |
| **HDiv-VIM / reduced FEM** | Magnetic materials | Nonlinear B-H, saturation |
| **Surface Zs** | Skin effect | Arbitrary shapes via PyKAN |

### Frequency Range and Physics

| Frequency | Dominant Effect | Solver Component |
|-----------|-----------------|------------------|
| DC - 1 kHz | DC resistance | R_dc only |
| 1 kHz - 100 kHz | Skin effect begins | Dowell/ESIM for Zs |
| 100 kHz - 1 MHz | Strong skin effect | Zs dominates R |
| 1 MHz - 10 MHz | Capacitive effects | PEEC Star (ε coupling) |

---

## PEEC Basic Formulation

### Standard PEEC Impedance Matrix

For a conductor segment, the PEEC method constructs an impedance matrix:

```
V = Z · I

Z = R + jωL
```

where:
- **R**: Resistance matrix (diagonal for self, off-diagonal for mutual resistance)
- **L**: Inductance matrix (self and mutual inductance from Neumann formula)

### DC Resistance

For a conductor segment with length l, cross-section area A, and conductivity σ:

```
R_dc = ρ · l / A = l / (σ · A)
```

### External Inductance (Neumann Formula)

The mutual inductance between two conductor segments:

```
L_ij = (μ_0 / 4π) ∫∫ (dl_i · dl_j) / |r_i - r_j|
```

For self-inductance (i = j), the GMD (Geometric Mean Distance) is used to avoid singularity.

### Internal Inductance

For a round wire at DC:
```
L_int,dc = μ_0 · l / (8π)
```

For a rectangular conductor at DC:
```
L_int,dc = μ · d / 3   [H·m] (per unit area)
```

where d is the thickness in the direction of current skin effect.

### Full PEEC Matrix Assembly

```
Z_PEEC = R_dc + Z_s(ω) + jω·L_ext + jω·L_int(ω)
       = R_dc·F(ξ) + jω·(L_ext + L_int,dc·G(ξ))
```

where F(ξ) and G(ξ) are the Dowell factors for skin effect.

---

## Incorporating Surface Impedance

### Physical Model

At high frequencies, current concentrates near the conductor surface (skin effect). The surface impedance Z_s relates the tangential electric field E_t to the tangential magnetic field H_t at the surface:

```
E_t = Z_s · H_t
```

For a semi-infinite conductor:
```
Z_s = (1 + j) / (σ · δ) = √(jωμ/σ)
```

where δ = √(2ρ/(ωμ)) is the skin depth.

### Skin Depth Formula

The skin depth δ is the characteristic penetration distance of electromagnetic fields into a conductor:

```
δ = √(2ρ / (ωμ)) = √(2 / (ωμσ)) = √(ρ / (πfμ))
```

| Material | σ [S/m] | μ_r | δ at 50 Hz | δ at 1 kHz | δ at 100 kHz | δ at 1 MHz |
|----------|---------|-----|------------|------------|--------------|------------|
| Copper | 5.8×10^7 | 1 | 9.4 mm | 2.1 mm | 0.21 mm | 66 μm |
| Aluminum | 3.5×10^7 | 1 | 12 mm | 2.7 mm | 0.27 mm | 85 μm |
| Steel (μ_r=200) | 2×10^6 | 200 | 0.5 mm | 0.11 mm | 11 μm | 3.5 μm |
| Ferrite | 0.01 | 2000 | 503 mm | 113 mm | 11.3 mm | 3.6 mm |

### From Surface Impedance to Conductor Resistance

For a conductor with perimeter P and length l:

```
Z_conductor = V / I = Z_s · l / P
R_ac = Re(Z_conductor) = Re(Z_s) · l / P
```

### Impedance Matrix with Skin Effect

The PEEC impedance matrix becomes:

```
Z_ii = R_ac,i + jωL_ii + jωL_int,i
```

where:
- **R_ac,i**: AC resistance of segment i
- **L_ii**: External self-inductance (from Neumann formula)
- **L_int,i**: Internal inductance (from skin effect)

### SIBC (Surface Impedance Boundary Condition)

In 3D PEEC, SIBC replaces volumetric current distribution with surface current:

```
Traditional PEEC:
  J(r) = 3D current density distribution (expensive mesh inside conductor)

PEEC + SIBC:
  J_s(r) = Surface current distribution
  Z_s = √(jωμ/σ)  (replaces conductor interior)
```

**Advantages**:
- No mesh inside conductor (only surface mesh)
- Frequency-dependent skin effect captured analytically
- Orders of magnitude faster for high-conductivity materials

---

## Dowell's Formula for Rectangular Conductors

### Physical Basis

Dowell's formula derives from solving the 1D diffusion equation for a rectangular conductor with half-thickness $a$:

$$
\frac{\partial^2 H}{\partial z^2} = j\omega\mu\sigma H
$$

with boundary conditions $H(z=\pm a) = H_0$.

The solution gives the complex current density distribution, from which resistance and inductance ratios are computed.

### Resistance Ratio $F(\xi)$

$$
F(\xi) = \xi \cdot \frac{\sinh(2\xi) + \sin(2\xi)}{\cosh(2\xi) - \cos(2\xi)}
$$

Alternative form using complex propagation constant $\gamma = (1+j)/\delta$:

$$
F(\xi) = \text{Re}[\gamma a \cdot \coth(\gamma a)]
$$

where $\xi = a/\delta$ is the normalized thickness (half-thickness divided by skin depth).

### Internal Inductance Ratio $G(\xi)$

$$
G(\xi) = \frac{3}{2\xi} \cdot \frac{\sinh(2\xi) - \sin(2\xi)}{\cosh(2\xi) - \cos(2\xi)}
$$

The internal inductance represents energy stored in the magnetic field **inside** the conductor.

### Taylor Series Expansions (Low Frequency)

For $\xi < 0.5$:

$$
F(\xi) \approx 1 + \frac{\xi^4}{45} - \frac{2\xi^8}{4725} + O(\xi^{12})
$$

$$
G(\xi) \approx 1 - \frac{\xi^4}{15} + \frac{2\xi^8}{945} + O(\xi^{12})
$$

### High Frequency Asymptotic (Large $\xi$)

For $\xi > 5$:

$$
F(\xi) \approx \xi
$$

$$
G(\xi) \approx \frac{3}{2\xi}
$$

### Numerical Values

| ξ | F(ξ) R_ac/R_dc | G(ξ) L_int/L_int,dc |
|---|----------------|---------------------|
| 0.1 | 1.000 | 1.000 |
| 0.5 | 1.003 | 0.975 |
| 1.0 | 1.086 | 0.806 |
| 2.0 | 1.932 | 0.476 |
| 3.0 | 3.010 | 0.329 |
| 5.0 | 5.004 | 0.200 |
| 10.0 | 10.001 | 0.100 |

**Asymptotic Behavior**:
- DC (ξ → 0): F → 1, G → 1
- High frequency (ξ → ∞): F → ξ, G → 3/(2ξ)

### Python Implementation

```python
import numpy as np

MU_0 = 4 * np.pi * 1e-7  # Permeability of free space [H/m]

def calc_skin_depth(freq, sigma, mu=MU_0):
    """Calculate skin depth delta = sqrt(2 / (omega * mu * sigma))."""
    omega = 2 * np.pi * freq
    omega = np.maximum(omega, 1e-10)  # Avoid division by zero
    return np.sqrt(2.0 / (omega * mu * sigma))

def dowell_F(xi):
    """Dowell's formula for R_ac/R_dc.

    F(xi) = xi * (sinh(2*xi) + sin(2*xi)) / (cosh(2*xi) - cos(2*xi))

    Args:
        xi: Normalized thickness = half_thickness / skin_depth
    Returns:
        Ratio R_ac/R_dc
    """
    if np.isscalar(xi):
        if xi < 0.01:
            return 1.0 + xi**4 / 45
        sh2, sn2 = np.sinh(2*xi), np.sin(2*xi)
        ch2, cs2 = np.cosh(2*xi), np.cos(2*xi)
        return xi * (sh2 + sn2) / (ch2 - cs2)
    else:
        result = np.ones_like(xi, dtype=float)
        mask = xi >= 0.01
        xi_m = xi[mask]
        sh2, sn2 = np.sinh(2*xi_m), np.sin(2*xi_m)
        ch2, cs2 = np.cosh(2*xi_m), np.cos(2*xi_m)
        result[mask] = xi_m * (sh2 + sn2) / (ch2 - cs2)
        return result

def dowell_G(xi):
    """Dowell's formula for L_int/L_int,dc.

    G(xi) = (3/2xi) * (sinh(2*xi) - sin(2*xi)) / (cosh(2*xi) - cos(2*xi))

    Args:
        xi: Normalized thickness = half_thickness / skin_depth
    Returns:
        Ratio L_int/L_int,dc
    """
    if np.isscalar(xi):
        if xi < 0.01:
            return 1.0 - xi**4 / 15
        sh2, sn2 = np.sinh(2*xi), np.sin(2*xi)
        ch2, cs2 = np.cosh(2*xi), np.cos(2*xi)
        return (3 / (2*xi)) * (sh2 - sn2) / (ch2 - cs2)
    else:
        result = np.ones_like(xi, dtype=float)
        mask = xi >= 0.01
        xi_m = xi[mask]
        sh2, sn2 = np.sinh(2*xi_m), np.sin(2*xi_m)
        ch2, cs2 = np.cosh(2*xi_m), np.cos(2*xi_m)
        result[mask] = (3 / (2*xi_m)) * (sh2 - sn2) / (ch2 - cs2)
        return result

def dowell_impedance(R_dc, L_int_dc, thickness, sigma, mu, freq):
    """Calculate impedance using Dowell's formula.

    Z(s) = R_dc * F(xi) + j*omega * L_int_dc * G(xi)

    Args:
        R_dc: DC resistance [Ohm]
        L_int_dc: DC internal inductance [H]
        thickness: Conductor thickness [m]
        sigma: Conductivity [S/m]
        mu: Permeability [H/m]
        freq: Frequency [Hz]
    Returns:
        Complex impedance [Ohm]
    """
    omega = 2 * np.pi * freq
    delta = calc_skin_depth(freq, sigma, mu)
    xi = (thickness / 2) / delta  # Half-thickness / skin_depth
    F = dowell_F(xi)
    G = dowell_G(xi)
    return R_dc * F + 1j * omega * L_int_dc * G
```

### Usage Example

```python
# Copper foil: 0.1 mm thick, 10 mm wide, 100 mm long
thickness = 0.1e-3  # m
width = 10e-3       # m
length = 100e-3     # m
sigma = 5.8e7       # S/m (copper)
mu = MU_0

# DC parameters
R_dc = length / (sigma * width * thickness)  # Ohm
L_int_dc = mu * length * thickness / (3 * width)  # H (approximate)

# Impedance at 100 kHz
Z = dowell_impedance(R_dc, L_int_dc, thickness, sigma, mu, 100e3)
print(f"Z at 100 kHz: {Z.real:.6f} + j{Z.imag:.6f} Ohm")
print(f"R_ac/R_dc = {Z.real/R_dc:.3f}")
```

---

## Dowell vs ESIM: Method Comparison

### Boundary Condition Difference

**ESIM (dH/dz = 0 BC) and Dowell's formula (H = 0 BC) solve DIFFERENT problems!**

| Formula | Function | BC at Center | DC Limit | Use Case |
|---------|----------|--------------|----------|----------|
| Dowell | **coth** | H(a) = 0 | F = 1 | R_ac/R_dc ratio |
| ESIM | **tanh** | dH/dz(a) = 0 | Z_s = 0 | Surface impedance |

Note: coth(x) = 1/tanh(x), reflecting different boundary conditions.

### Comparison Table

| Aspect | Dowell | ESIM |
|--------|--------|------|
| **Output** | F, G (real ratios) | Z_s (complex) |
| **Cross-section** | Rectangular only | Any (with 2D FEM) |
| **Nonlinear μ(H)** | Not supported | Supported |
| **DC limit** | F=1, G=1 (correct) | Z_s→0 (need homogenization) |
| **Computation** | Analytical (fast) | Numerical (slower) |

### When to Use Which

**Use Dowell when**:
- Linear material (constant μ)
- Rectangular cross-section
- Speed is important

**Use ESIM when**:
- Nonlinear material μ(H) (saturable iron)
- Need Z_s for SIBC boundary conditions
- Arbitrary cross-section (future 2D FEM)

---

## Cross-Section Specific Formulas

### Rectangular Cross-Section (Dowell)

For width w >> height h (thin strip), use Dowell's formula with:
- Half-thickness: a = h / 2
- ξ = a / δ

### Round Wire (Kelvin Functions)

For circular cross-section, use Kelvin functions (ber, bei):

```python
from scipy.special import kelvin

def kelvin_rac_ratio(xi_sqrt2):
    """R_ac/R_dc for round wire. xi_sqrt2 = sqrt(2) * radius / delta"""
    ber, bei, _, _ = kelvin(xi_sqrt2)
    h = 1e-6
    ber_h, bei_h, _, _ = kelvin(xi_sqrt2 + h)
    berp, beip = (ber_h - ber) / h, (bei_h - bei) / h
    return (xi_sqrt2 / 2) * (complex(berp, beip) / complex(ber, bei)).real
```

### Geometry Summary

| Cross-Section | Formula | Parameter |
|--------------|---------|-----------|
| Rectangular | Dowell (planar 1D) | ξ = half-thickness / δ |
| Round | Kelvin functions | ξ = √2 × radius / δ |
| Square | ~1.05 × Dowell | Approximate correction |
| Litz wire | Complex formula | Depends on strand arrangement |
| Arbitrary | 2D FEM required | No analytical formula |

---

## ESIM Homogenization for Nonlinear Materials

### Why Homogenization?

For linear materials (constant μ), ESIM and Dowell give equivalent results.

**ESIM value is for NONLINEAR materials** where μ(H) varies with field strength.

### Homogenization Algorithm

1. **Solve 1D cell problem** with BC: H(0)=H0, dH/dz(a)=0
   - Get H(z) and mu(z) distributions

2. **Compute effective permeability** (|H|²-weighted average):
   ```
   mu_eff = integral{|H|^2 * mu(z)} dz / integral{|H|^2} dz
   ```

3. **Compute effective xi**:
   ```
   delta_eff = sqrt(2*rho / (omega * mu_eff))
   xi_eff = a / delta_eff
   ```

4. **Apply Dowell's formula**:
   ```
   R_ac/R_dc = F(xi_eff)
   ```

### Validation Results (Linear Material)

| ξ | Dowell | Direct Z_s | Homogenization |
|----|--------|------------|----------------|
| 0.3 | 1.0007 | 1.3335 | **1.0007** |
| 1.0 | 1.0856 | 1.3513 | **1.0856** |
| 3.0 | 3.0101 | 2.4779 | **3.0101** |

The homogenization approach matches Dowell exactly for linear materials.

### Nonlinear Material Example (Steel BH curve)

| H0 [A/m] | μ_r_eff | ξ_eff | R_ac/R_dc |
|----------|---------|-------|-----------|
| 10 | 1588 | 25.0 | 25.0 |
| 1000 | 1276 | 22.4 | 22.4 |
| 10000 | 217 | 9.3 | 9.3 |

As material saturates (high H0), μ_r decreases, ξ_eff decreases, and skin effect reduces.

---

## Loop-Star Decomposition and Z_s

In PEEC with Loop-Star decomposition:

```
Z_PEEC = [Z_LL  Z_LS]
         [Z_SL  Z_SS]
```

**Z_s affects only Loop components**:

| Component | Z_s Effect | Physical Reason |
|-----------|------------|-----------------|
| Z_LL (Loop-Loop) | **Modified** | Current flows through conductor |
| Z_LS, Z_SL (Coupling) | Minimal | Geometric coupling |
| Z_SS (Star-Star) | **None** | Capacitive (charge accumulation) |

Star components represent **charge accumulation** and **displacement current**, which are not affected by conductor skin effect.

---

## Frequency Dependence Analysis

### Fundamental Frequency-Dependent Quantities

All frequency dependence originates from ω = 2πf:

```
ω = 2πf                                    [rad/s]
δ(ω) = √(2ρ / (ωμ)) = √(ρ / (πfμ))        ∝ ω^(-1/2)  [m]
ξ(ω) = a / δ(ω) = a√(πfμ/ρ)               ∝ ω^(+1/2)  [-]
γ(ω) = (1+j) / δ(ω) = (1+j)√(ωμσ/2)       ∝ ω^(+1/2)  [1/m]
```

### PEEC Matrix Element Frequency Dependence

```
Z_PEEC(ω) = R_dc · F(ξ(ω)) + jω(L_ext + L_int(ω))
```

| Component | Symbol | Low Frequency | High Frequency |
|-----------|--------|---------------|----------------|
| DC Resistance | R_dc | constant | constant |
| AC/DC Ratio | F(ξ) | ≈ 1 | ∝ ω^(1/2) |
| **AC Resistance** | R_ac = R_dc·F | **≈ R_dc** | **∝ ω^(1/2)** |
| External Inductance | jωL_ext | ∝ ω | ∝ ω |
| Internal Inductance | jωL_int | ≈ 0 | ∝ ω^(1/2) |

### Practical Implications

1. **Low frequency (f < f_transition)**:
   - Skin effect negligible
   - R_ac ≈ R_dc (use DC resistance)

2. **Transition frequency**:
   ```
   f_transition ≈ ρ / (πμa²)

   Example (copper, a = 1mm):
   f_transition ≈ 4.3 kHz
   ```

3. **High frequency (f > f_transition)**:
   - Strong skin effect
   - R_ac increases as √f

---

## Laplace Domain Representation

### Fundamental Relationship

The Dowell functions F and G can be expressed in the Laplace domain using the magnetic diffusion time constant:

```
τ = a²μσ    [s]  (magnetic diffusion time constant)
s = jω      (Laplace variable)
```

### F(s) and G(s) in Laplace Domain

```
Z(s)/R_dc = √(τs) · coth(√(τs))
```

### Continued Fraction Expansion

The function √(τs)·coth(√(τs)) has the elegant continued fraction:

```
√(τs)·coth(√(τs)) = 1 + τs/(3 + τs/(5 + τs/(7 + τs/(9 + ...))))
```

**Truncated approximations**:

| Order | Approximation | Accuracy |
|-------|---------------|----------|
| 0th | F₀ = 1 | DC exact |
| 1st | F₁ = 1 + τs/3 | Good for τω < 1 |
| 2nd | F₂ = (15 + 8τs)/(15 + 3τs) | Good for τω < 4.5 (ξ < 1.5) |

### PRIMA Ladder Network

The continued fraction maps to a ladder circuit:

```
     R_dc     L₁      L₂      L₃
  o──/\/\/──○──⊃⊃⊃──○──⊃⊃⊃──○──⊃⊃⊃──○
             │       │       │
            ═╧═     ═╧═     ═╧═
            R₁      R₂      R₃

where:
  L₁ = τ·R_dc/3,  R₁ = 3R_dc
  L₂ = τ·R_dc/5,  R₂ = 5R_dc/3
  L₃ = τ·R_dc/7,  R₃ = 7R_dc/5
```

### Time Constant Examples

| Material | a | σ [S/m] | μ_r | τ = a²μσ |
|----------|---|---------|-----|----------|
| Copper, 1mm | 1mm | 5.8×10⁷ | 1 | **73 ns** |
| Copper, 5mm | 5mm | 5.8×10⁷ | 1 | **1.8 μs** |
| Steel, 1mm | 1mm | 2×10⁶ | 1000 | **2.5 ms** |

### Python Implementation

```python
def dowell_F_continued_fraction(tau_s, order=3):
    """Continued fraction approximation of F(s)."""
    denominators = [2*k + 1 for k in range(order, 0, -1)]
    result = denominators[0]
    for d in denominators[1:]:
        result = d + tau_s / result
    return 1 + tau_s / result

def dowell_F_2nd_order(tau_s):
    """2nd order rational approximation of F(s)."""
    return (15 + 8*tau_s) / (15 + 3*tau_s)
```

---

## Proximity Effect in PEEC + Dowell

### Overview

When combining PEEC with Dowell's formula, proximity effects are **partially included**:

| Effect | Dowell Alone | PEEC + Dowell |
|--------|:------------:|:-------------:|
| Self skin effect | ✓ | ✓ |
| External proximity effect | ✗ | **✓** |
| Internal proximity effect | ✓ (m-layer formula) | △ (approximate) |

### External Proximity Effect

Current redistribution due to magnetic fields from **other conductors**.

**PEEC + Dowell**: Automatically included through PEEC mutual inductance terms.

### Practical Guidelines

| Situation | Recommended Approach |
|-----------|---------------------|
| Conductor spacing > 3δ | PEEC + Dowell is sufficient |
| Conductor spacing ≈ δ | Consider filament subdivision |
| Closely wound coils | Use Dowell's m-layer formula |
| High accuracy required | Filament subdivision or FEM |

### Improvement Methods

#### Method 1: Filament Subdivision

Subdivide one conductor into multiple filaments. Each filament is an independent PEEC element, and proximity effect is automatically included.

#### Method 2: ESIM with Boundary Magnetic Field

For nonlinear materials or complex geometries:
- Include external field in ESIM boundary conditions
- Solve with asymmetric boundary: H(+a) ≠ H(-a)

---

## PRIMA Model Order Reduction for PEEC

### Overview

PRIMA (Passive Reduced-order Interconnect Macromodeling Algorithm) reduces large PEEC systems while preserving terminal impedance. This is essential for SPICE simulation where full PEEC matrices are too large.

**Key Assumption**: **Uniform conductor** (same material and cross-section for all segments)

This assumption ensures:
- tau = a^2*mu*sigma is identical for all segments
- F(s) and G(s) become scalar multipliers
- Diagonal structure is **exactly preserved** after Lanczos transformation

### Why Model Order Reduction?

| Scenario | Full PEEC DOF | SPICE Simulation |
|----------|---------------|------------------|
| Single coil (20 turns) | 20 | Manageable |
| Multi-layer transformer | 200-500 | Slow |
| Complex PCB traces | 1000+ | Impractical |
| WPT system (TX+RX) | 100+ | Needs reduction |

PRIMA reduces DOF by 10-100x while preserving terminal behavior.

### Problem Statement

For N conductor segments:
```
V = Z(s) * I
Z(s) = R(s) + sL
```

The goal is to find a reduced-order model:
```
V_r = Z_r(s) * I_r
```
where Z_r has dimension k << N but matches terminal impedance.

### PRIMA Reduction Strategy: DC Basis with Diagonal s-Correction

1. **Generate basis at DC** (s=0)
2. **Apply Dowell correction as diagonal s-dependent terms**

This preserves the simple structure where **only diagonal elements have s-dependence**.

### Step 1: DC Basis Generation

```
K = L_ext (external inductance matrix, dense with mutual inductances)
N = R_dc  (DC resistance matrix, diagonal)
v0 = terminal excitation vector

Lanczos transformation at DC:
  U, V = lanczos_tridiagonalize(K, N, v0)
```

### Step 2: Transform to Reduced Coordinates

```
L' = U^H * L_ext * V  -> Tridiagonal (constant)
R' = U^H * R_dc * V   -> Diagonal (constant)
```

### Step 3: Final Reduced System Structure

```
Z'(s) = s * L'_tridiag + R'_diag + Delta_diag(s)

where:
  L'_tridiag: Constant tridiagonal (from DC Lanczos)
  R'_diag:    Constant diagonal (from DC Lanczos)
  Delta_diag(s): s-dependent diagonal correction from Dowell
```

### Accuracy Analysis

**Lanczos Truncation Error**:

The Lanczos reduction preserves terminal impedance moments:
```
Z(s) = m0 + m1*s + m2*s^2 + ...

k_reduced stages preserve first 2*k_reduced moments
```

| k_reduced | Bandwidth | Notes |
|-----------|-----------|-------|
| 5 | Narrowband | Single resonance |
| 10 | Wideband | Multiple resonances |
| 20 | Very wide | Complex structures |

### Applications

**WPT Coil Design**:
```
Multi-turn WPT coil:
- N = 20 turns
- k_reduced = 10

Full PEEC: 20 DOF, slow frequency sweep
PRIMA reduced: 10 DOF, fast frequency sweep

Speedup: ~8x for impedance calculation
```

---

## ACA Low-Rank Approximation for Large-Scale PEEC

### Overview

For large PEEC systems (N > 1000 segments), the full inductance matrix $L_{ext}$ becomes memory-prohibitive ($O(N^2)$ storage). **ACA (Adaptive Cross Approximation)** provides efficient low-rank approximation for off-diagonal blocks.

### H-Matrix Structure

The PEEC inductance matrix can be hierarchically partitioned:

$$
L = \begin{pmatrix}
L_{11} & L_{12} & \cdots \\
L_{21} & L_{22} & \cdots \\
\vdots & \vdots & \ddots
\end{pmatrix}
$$

For well-separated clusters, off-diagonal blocks $L_{ij}$ are numerically low-rank:

$$
L_{ij} \approx U_k V_k^T \quad \text{where } k \ll \min(n_i, n_j)
$$

### ACA Algorithm

**Adaptive Cross Approximation (ACA+)** builds low-rank approximation without computing full matrix:

```
Input: Matrix block A (m x n), tolerance eps
Output: Low-rank factors U (m x k), V (n x k)

1. Initialize: R = A, k = 0
2. While ||R||_F > eps * ||A||_F:
   a. Find pivot (i*, j*) = argmax |R_ij|
   b. u_k = R[:, j*]          # Column
   c. v_k = R[i*, :] / R[i*, j*]  # Row (normalized)
   d. R = R - u_k * v_k^T     # Update residual
   e. k = k + 1
3. Return U = [u_1, ..., u_k], V = [v_1, ..., v_k]
```

### Admissibility Condition

Blocks are approximated if clusters are **well-separated**:

$$
\text{dist}(C_i, C_j) \geq \eta \cdot \max(\text{diam}(C_i), \text{diam}(C_j))
$$

where $\eta \approx 2.0$ is the admissibility parameter.

### Complexity Comparison

| Operation | Dense | H-Matrix (ACA) |
|-----------|-------|----------------|
| Storage | $O(N^2)$ | $O(N \log N)$ |
| Matrix-Vector | $O(N^2)$ | $O(N \log N)$ |
| LU Factorization | $O(N^3)$ | $O(N \log^2 N)$ |

### Integration with PRIMA

ACA acceleration can be combined with PRIMA reduction:

1. **Build H-matrix** for $L_{ext}$ using ACA
2. **Apply PRIMA** using H-matrix-vector products
3. **Result**: Reduced model from large-scale PEEC

```python
# Pseudo-code for H-matrix accelerated PRIMA
from hacapk import HMatrix

# Build H-matrix for L_ext
L_hmat = HMatrix(segments, eps=1e-4, eta=2.0)

# Lanczos with H-matrix MVP
def matvec(v):
    return L_hmat @ v

U, T = lanczos_tridiagonalize(matvec, R_dc, v0, k_reduced)
```

### Radia Implementation

Radia uses **HACApK** library for H-matrix operations:

```python
import radia as rad

# Enable H-matrix solver
rad.SolverConfig(hacapk_eps=1e-4, hacapk_leaf=10, hacapk_eta=2.0)

# Solve with H-matrix acceleration
rad.Solve(container, precision=0.0001, max_iter=1000, method=2)  # method=2 = HACApK
```

### When to Use ACA

| Problem Size | Recommended Solver |
|--------------|-------------------|
| N < 500 | Dense LU (fast for small N) |
| 500 < N < 2000 | BiCGSTAB (iterative) |
| N > 2000 | **HACApK (ACA)** |

---

## Magnetic-Core Coupling Policy

PEEC in Radia is now the conductor / circuit-reduction path. Magnetic material
coupling is not implemented through a PEEC-owned magnetic moment subsystem.
Use HDiv-VIM or reduced FEM for ferrite, soft iron, saturation, and open-boundary
magnetic-material response, then exchange only the required terminal quantities
or source fields with the PEEC circuit model.

The old monolithic PEEC-plus-magnetic-core derivation has been removed from the
public docs to avoid reviving the retired moment-path API. Keep PEEC examples
conductor-only unless a separate HDiv-VIM / reduced-FEM artifact supplies the
magnetic response.
## Application Examples

### 1. Wireless Power Transfer (WPT)

```
Configuration: TX coil + RX coil + Ferrite shield

PEEC: Coil mutual inductance
HDiv-VIM: Flux concentration by ferrite
Zs: Litz wire skin/proximity effect
```

### 2. Induction Heating

```
Configuration: Heating coil + Workpiece (conductive)

PEEC: Coil self-inductance
HDiv-VIM: Workpiece nonlinear mu (temperature-dependent)
Zs: Workpiece surface frequency-dependent impedance (ESIM)
```

### 3. High-Frequency Transformer

```
Configuration: Primary winding + Secondary winding + Magnetic core

PEEC: Winding coupling, leakage inductance
HDiv-VIM: Core magnetization, saturation
Zs: Winding skin/proximity effect
```

### 4. EMI/EMC Filter

```
Configuration: Common-mode choke + Wiring

PEEC: Wiring impedance, capacitive coupling
HDiv-VIM: Ferrite core frequency characteristics mu(omega)
Zs: Wiring cross-section high-frequency loss
```

---

## Future Roadmap

### Current Status (2026-01)

| Feature | Status | Notes |
|---------|--------|-------|
| 2D cross-section Zs (Dowell) | Implemented | Rectangular conductors |
| 2D cross-section Zs (ESIM) | Implemented | Nonlinear mu(H) |
| 2D cross-section Zs (PyKAN) | Design complete | Arbitrary shapes |
| PEEC Loop-Star | Implemented | Low-frequency stable |
| PRIMA MOR | Design complete | SPICE export |
| PEEC magnetic-core coupling | Retired | Use HDiv-VIM / reduced FEM |

### Planned: 3D Volumetric Effects

Current surface impedance assumes **1D or 2D** current distribution (depth direction only). For complex 3D geometries:

**3D Effects Not Yet Covered**:
- Corner effects in PCB vias
- 3D current crowding at bends
- Non-uniform cross-section along conductor length
- Proximity effect from multiple directions

**Planned Approach**:

$$
Z_{3D}(s) = Z_{2D}(s) \cdot K_{3D}(\text{geometry})
$$

where $K_{3D}$ is a geometry-dependent correction factor.

**Implementation Options**:

1. **3D FEM + PyKAN**: Extend NGSolve training to 3D eddy current problems
   - Input: [omega, sigma, geometry_params_3D...]
   - Output: [Re(Z), Im(Z)]
   - Computational cost: Higher training time, same inference time

2. **Correction Factors**: Analytical/semi-analytical corrections for common 3D effects
   - Corner correction: $K_{corner} = f(r/\delta)$ for bend radius r
   - Via correction: $K_{via} = f(d_{via}/\delta, AR)$ for via diameter and aspect ratio

3. **Sub-segmentation**: Divide conductor into smaller 2D segments
   - Each segment uses 2D Zs
   - 3D effects captured through PEEC mutual coupling

**Priority Applications for 3D**:

| Application | 3D Effect | Priority |
|-------------|-----------|----------|
| PCB vias | Current crowding at via barrel | High |
| Connector pins | Corner effects at bends | Medium |
| Busbars | Current redistribution at corners | Medium |
| Spiral inductors | 3D proximity in multi-layer | High |

### Planned: Temperature-Dependent Properties

Power electronics applications require temperature effects:

$$
\sigma(T) = \sigma_0 / (1 + \alpha (T - T_0))
$$

$$
\mu(T, H) = \mu(H) \cdot f(T)  \text{(Curie temperature effects)}
$$

**Integration with thermal solvers** (future):
- Coupled electro-thermal simulation
- Loss distribution -> Temperature rise -> Updated sigma, mu
- Iterative convergence

---

## Comparison with Other Tools

| Feature | Radia Integrated | ANSYS Maxwell | FastHenry | FEKO |
|---------|------------------|---------------|-----------|------|
| Magnetic materials | HDiv-VIM / reduced FEM | FEM | None | MoM |
| Conductors | PEEC+Zs | FEM | PEEC | MoM |
| Skin effect | PyKAN (arbitrary) | FEM | Analytical only | SIBC |
| SPICE output | Verilog-A | Limited | RLCK | None |
| Low-freq stability | Loop-Star | OK | Unstable | Unstable |
| Computation speed | Fast | Slow | Fast | Medium |
| License | Open source | Commercial | Academic | Commercial |

---

## Implementation Files

| File | Purpose |
|------|---------|
| `esim_cell_problem.py` | Base ESIM solvers |
| `esim_coupled_solver.py` | PEEC+ESIM coupled solver |
| `lanczos_reduction.py` | PRIMA/Lanczos MOR |

---

## References

### PEEC and Surface Impedance

1. **A. Ruehli**, "Equivalent Circuit Models for Three-Dimensional Multiconductor Systems," *IEEE Trans. Microwave Theory and Techniques*, vol. 22, no. 3, pp. 216-221, 1974.
   - Original PEEC formulation

2. **P.L. Dowell**, "Effects of eddy currents in transformer windings," *Proc. IEE*, vol. 113, no. 8, pp. 1387-1394, 1966.
   - Analytical skin effect formulas F(xi) and G(xi) for rectangular conductors

3. **J.A. Ferreira**, "Improved analytical modeling of conductive losses in magnetic components," *IEEE Trans. Power Electronics*, vol. 9, no. 1, pp. 127-131, 1994.
   - Extended Dowell analysis for transformer windings

### ESIM and Nonlinear Materials

4. **K. Hollaus, V. Hanser, and M. Schobinger**, "A Nonlinear Effective Surface Impedance in a Magnetic Scalar Potential Formulation," *IEEE Trans. Magnetics*, 2025.
   - ESIM homogenization for nonlinear mu(H) materials
   - Key reference for ESIM implementation in Radia

### Loop-Star Decomposition

5. **G. Vecchi**, "Loop-Star decomposition of basis functions in the discretization of the EFIE," *IEEE Trans. Antennas and Propagation*, vol. 47, no. 2, pp. 339-346, 1999.
   - Low-frequency stabilization for EFIE

6. **J.S. Zhao, W.C. Chew**, "Integral equation solution of Maxwell's equations from zero frequency to microwave frequencies," *IEEE Trans. Antennas and Propagation*, vol. 48, no. 10, pp. 1635-1645, 2000.
   - Augmented EFIE formulation

### Model Order Reduction

7. **A. Odabasioglu, M. Celik, L.T. Pileggi**, "PRIMA: Passive Reduced-order Interconnect Macromodeling Algorithm," *IEEE Trans. Computer-Aided Design*, vol. 17, no. 8, pp. 645-654, 1998.
   - Passivity-preserving model order reduction

8. **B. Gustavsen, A. Semlyen**, "Rational approximation of frequency domain responses by vector fitting," *IEEE Trans. Power Delivery*, vol. 14, no. 3, pp. 1052-1061, 1999.
   - Vector Fitting for rational approximation

### PyKAN (Kolmogorov-Arnold Networks)

9. **Z. Liu, Y. Wang, S. Vaidya, F. Ruehle, J. Halverson, M. Soljacic, T.Y. Hou, M. Tegmark**, "KAN: Kolmogorov-Arnold Networks," *arXiv:2404.19756*, 2024.
   - Original KAN paper with symbolic formula extraction
   - Key reference for PyKAN-based surface impedance learning

10. **PyKAN GitHub Repository**: https://github.com/KindXiaoming/pykan
    - Official implementation of Kolmogorov-Arnold Networks
    - Used for learning arbitrary cross-section surface impedance

### Complex Permeability and Ferrite Materials

11. **E.C. Snelling**, "Soft Ferrites: Properties and Applications," *Butterworth-Heinemann*, 1988.
    - Complex permeability mu' - j*mu" for ferrite materials

12. **C.R. Sullivan**, "Optimal choice for number of strands in a Litz-wire transformer winding," *IEEE Trans. Power Electronics*, vol. 14, no. 2, pp. 283-291, 1999.
    - Litz wire skin and proximity effect analysis

---

**Date**: 2026-01-16
**Author**: Claude Code
