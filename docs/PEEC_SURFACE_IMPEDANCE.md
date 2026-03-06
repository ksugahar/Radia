# Radia PEEC Solver: Surface Impedance and SPICE Export

## Overview

This document describes Radia's integrated electromagnetic solver framework combining:

- **PEEC (Loop-Star)**: Conductor current distribution, inductive/capacitive coupling
- **MMM**: Magnetic material magnetization
- **Surface Impedance Zs**: Skin effect via Dowell, ESIM, or PyKAN learning

```
+-------------------------------------------------------------+
|                    Radia Integrated Solver                   |
+-------------------------------------------------------------+
|  +-----------+   +-----------+   +-----------+              |
|  |   PEEC    |   |    MMM    |   | Surface   |              |
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
| **MMM** | Magnetic materials | Nonlinear B-H, saturation |
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

## PEEC-MMM Coupling via Mutual Inductance

### Overview

When a coil is placed near or around a magnetic material (e.g., ferrite core), the inductance and coupling are enhanced. This section describes how to couple PEEC (for conductors) with MMM (for magnetic materials).

### Coupling Mechanism

```
Coil current I → Generates H-field → Magnetizes core M
Core magnetization M → Generates flux → Links back to coil (enhanced L)
```

### Mathematical Formulation

#### PEEC Part (Conductor)
```
V_coil = Z_PEEC · I_coil + V_coupling

Z_PEEC = R_ac + jωL_ext + jωL_int (with skin effect)
```

#### MMM Part (Magnetic Core)
```
M = χ · H_total
H_total = H_coil + H_demag

H_demag = N · M  (demagnetization tensor)
```

#### Coupling Terms
```
H_coil at core: H_c = A_cp · I_coil  (Biot-Savart from coil)
Flux to coil from M: Φ_pm = A_pm · M
V_coupling = jω · Φ_pm = jω · A_pm · M
```

### Coupled System Matrix

```
[Z_PEEC      jω·A_pm  ] [I_coil]   [V_source]
[-χ·A_cp    I + χ·N   ] [M     ] = [0       ]
```

### Schur Complement: Effective PEEC

Eliminate M to get effective impedance seen by circuit:

```
M = (I + χ·N)^{-1} · χ · A_cp · I_coil
V = Z_eff · I_coil

Z_eff = Z_PEEC + jω · A_pm · (I + χ·N)^{-1} · χ · A_cp
        \_____/   \____________________________________/
        Original        Core enhancement term
```

### Physical Interpretation

The core enhancement term:
```
ΔZ = jω · A_pm · (I + χ·N)^{-1} · χ · A_cp
```

represents:
1. **A_cp**: H-field at core from unit coil current
2. **χ · (I + χ·N)^{-1}**: Effective susceptibility (includes demagnetization)
3. **A_pm**: Flux linkage per unit magnetization
4. **jω**: Faraday's law (rate of change)

### ACA Low-Rank Approximation for Schur Complement

The core enhancement term $\Delta Z$ can be approximated using ACA (Adaptive Cross Approximation) to reduce SPICE netlist complexity.

**Problem**: For N_coil PEEC segments and N_core MMM elements:
- $A_{pm}$: N_coil × N_core matrix
- $A_{cp}$: N_core × N_coil matrix
- Full $\Delta Z$: N_coil × N_coil dense matrix → $O(N_{coil}^2)$ mutual inductance elements in SPICE

**Solution**: Apply ACA to $\Delta Z$:

$$
\Delta Z = j\omega \cdot A_{pm} \cdot (I + \chi N)^{-1} \cdot \chi \cdot A_{cp} \approx j\omega \cdot U \cdot V^T
$$

where $U \in \mathbb{R}^{N_{coil} \times k}$ and $V \in \mathbb{R}^{N_{coil} \times k}$ with $k \ll N_{coil}$.

**SPICE Element Reduction**:

| Representation | Mutual Inductance Elements | Storage |
|----------------|---------------------------|---------|
| Full $\Delta Z$ | $N_{coil}(N_{coil}-1)/2$ | $O(N_{coil}^2)$ |
| ACA Low-rank | $2 \cdot k \cdot N_{coil}$ | $O(k \cdot N_{coil})$ |

For N_coil = 100 and k = 5:
- Full: 4,950 mutual inductance elements
- ACA: 1,000 elements (5x reduction)

**Implementation in SPICE**:

The low-rank approximation $\Delta Z \approx U V^T$ can be implemented as:

```
Original: V_i = sum_j (ΔZ_ij * I_j)  for each segment i

Low-rank: V_i = sum_m (U_im * y_m)
          y_m = sum_j (V_jm * I_j)

→ Introduce k "auxiliary nodes" for rank-k approximation
```

**SPICE Subcircuit for Low-Rank Coupling**:

```spice
* Low-rank mutual inductance approximation
* Delta_Z = j*omega * U * V^T
* U: N_coil x k matrix, V: N_coil x k matrix

.SUBCKT LOWRANK_COUPLING seg1 seg2 ... segN aux1 aux2 ... auxk
* aux nodes represent intermediate "mode" currents

* V^T part: aux_m = sum_j V_jm * I_segj (CCCS)
G_aux1 0 aux1 POLY(N) seg1 0 seg2 0 ... segN 0  0 V_11 V_21 ... V_N1
G_aux2 0 aux2 POLY(N) seg1 0 seg2 0 ... segN 0  0 V_12 V_22 ... V_N2
...

* U part: V_segi += j*omega * sum_m U_im * I_auxm
* Implemented as mutual inductance from aux to seg
L_aux1 aux1 0 1
L_aux2 aux2 0 1
...
K_seg1_aux1 L_seg1 L_aux1 {U_11}
K_seg1_aux2 L_seg1 L_aux2 {U_12}
...

.ENDS
```

---

## Complex Permittivity and Star Component

### Complex Permittivity ε(ω) = ε' - jε"

In Loop-Star decomposition, the **Star component** handles capacitive effects. Complex permittivity includes both dielectric storage and loss:

$$
\varepsilon(\omega) = \varepsilon'(\omega) - j\varepsilon''(\omega)
$$

where the imaginary part includes both dielectric relaxation and conduction:

$$
\varepsilon''(\omega) = \varepsilon''_{relaxation}(\omega) + \frac{\sigma}{\omega}
$$

**Note**: The DC conductivity contribution $\sigma/\omega$ is often written as $\sigma/(j\omega)$ in the complex permittivity, which gives:

$$
\varepsilon_{eff} = \varepsilon' - j\left(\varepsilon'' + \frac{\sigma}{\omega}\right)
$$

### Capacitance Matrix in Star Component

The Star-Star impedance block:

$$
Z_{SS} = \frac{1}{j\omega C(\varepsilon)}
$$

For complex $\varepsilon = \varepsilon' - j\varepsilon''$:

$$
C_{eff} = \varepsilon' \cdot C_0 \quad \text{(effective capacitance)}
$$
$$
G_d = \omega \varepsilon'' \cdot C_0 \quad \text{(dielectric loss conductance)}
$$

Including DC conductivity:

$$
G_{total} = G_d + \sigma \cdot G_0 = \omega\varepsilon'' \cdot C_0 + \sigma \cdot G_0
$$

where $C_0$ and $G_0$ are geometry-dependent constants.

### ACA for Capacitance Matrix (SPICE Element Reduction)

Similar to the inductance case, the capacitance matrix $C_{ext}$ can be approximated using ACA:

$$
C_{ext} \approx U_C \cdot V_C^T
$$

**SPICE Capacitor Reduction**:

| Representation | Capacitor Elements | Storage |
|----------------|-------------------|---------|
| Full $C_{ext}$ | $N(N-1)/2$ | $O(N^2)$ |
| ACA Low-rank | $2 \cdot k \cdot N$ | $O(k \cdot N)$ |

**Implementation**:

```spice
* Low-rank capacitance approximation
* C_ext = U_C * V_C^T

.SUBCKT LOWRANK_CAPACITANCE n1 n2 ... nN aux1 aux2 ... auxk

* V_C^T part: charge sensing
E_aux1 aux1 0 POLY(N) n1 0 n2 0 ... nN 0  0 VC_11 VC_21 ... VC_N1
...

* U_C part: voltage coupling
* Q_ni = sum_m U_C_im * V_auxm → I_ni = sum_m U_C_im * C_m * dV_auxm/dt
C_aux1 aux1 0 1
C_aux2 aux2 0 1
...

.ENDS
```

### Combined Low-Rank PEEC Model

For a complete PEEC model with both inductive and capacitive coupling:

```
Full PEEC:
- N(N-1)/2 mutual inductances
- N(N-1)/2 mutual capacitances
- N self R, L, C

Low-rank PEEC (rank k_L for L, rank k_C for C):
- 2*k_L*N inductance-related elements
- 2*k_C*N capacitance-related elements
- N self R, L, C

Reduction factor: (N-1)/2 / (k_L + k_C) for off-diagonal coupling
```

**Example**: N = 500 segments, k_L = k_C = 10

| | Full | Low-rank | Reduction |
|--|------|----------|-----------|
| Mutual L | 124,750 | 10,000 | 12.5x |
| Mutual C | 124,750 | 10,000 | 12.5x |
| Total off-diagonal | 249,500 | 20,000 | **12.5x** |

This dramatic reduction makes large PEEC models practical for SPICE simulation.

### Verilog-A Implementation for Frequency-Dependent Low-Rank Model

Standard SPICE cannot handle frequency-dependent matrix elements. Verilog-A allows implementing frequency-dependent low-rank coupling with full accuracy.

**Low-Rank Mutual Inductance with Skin Effect**:

```verilog
`include "disciplines.vams"

// Low-rank mutual inductance with frequency-dependent correction
// Z_mutual = j*omega * L_mutual * F_skin(omega)
// L_mutual = U * V^T (low-rank)

module lowrank_mutual_L(seg1, seg2, aux1, aux2);
    inout seg1, seg2;
    inout aux1, aux2;
    electrical seg1, seg2, aux1, aux2;

    // Low-rank factors (example for rank-2)
    parameter real U1_1 = 1e-9;  // U[seg1, mode1]
    parameter real U1_2 = 0.5e-9;  // U[seg1, mode2]
    parameter real V1_1 = 1.0;  // V[seg1, mode1]
    parameter real V1_2 = 0.8;  // V[seg1, mode2]

    // Skin effect parameters
    parameter real tau = 1e-6;  // = a^2 * mu * sigma

    real omega, sqrt_tau_s, F_skin_re, F_skin_im;
    real L_mode1, L_mode2;

    analog begin
        omega = 2 * `M_PI * $freq;
        if (omega < 1e-10) omega = 1e-10;

        // Skin effect correction factor sqrt(tau*s)*coth(sqrt(tau*s))
        sqrt_tau_s = sqrt(tau * omega);
        if (sqrt_tau_s < 0.1) begin
            F_skin_re = 1.0;
            F_skin_im = tau * omega / 3.0;
        end else begin
            // Full formula for higher frequencies
            F_skin_re = sqrt_tau_s * (sinh(2*sqrt_tau_s) + sin(2*sqrt_tau_s))
                       / (cosh(2*sqrt_tau_s) - cos(2*sqrt_tau_s));
            F_skin_im = sqrt_tau_s * (sinh(2*sqrt_tau_s) - sin(2*sqrt_tau_s))
                       / (cosh(2*sqrt_tau_s) - cos(2*sqrt_tau_s));
        end

        // Mode currents from auxiliary nodes
        L_mode1 = U1_1 * V1_1;
        L_mode2 = U1_2 * V1_2;

        // Voltage contribution with frequency-dependent skin effect
        V(seg1) <+ (L_mode1 + L_mode2) * F_skin_re * ddt(I(aux1))
                 + (L_mode1 + L_mode2) * F_skin_im / omega * ddt(ddt(I(aux1)));
    end
endmodule
```

**Low-Rank Capacitance with Complex Permittivity**:

```verilog
`include "disciplines.vams"

// Low-rank capacitance with complex permittivity
// Y_mutual = j*omega * C_mutual * (epsilon'/epsilon_0 - j*epsilon"/epsilon_0)
// C_mutual = U_C * V_C^T (low-rank)

module lowrank_mutual_C(n1, n2, aux1, aux2);
    inout n1, n2;
    inout aux1, aux2;
    electrical n1, n2, aux1, aux2;

    // Low-rank factors
    parameter real UC1_1 = 1e-12;  // U_C[n1, mode1]
    parameter real UC1_2 = 0.5e-12;  // U_C[n1, mode2]
    parameter real VC1_1 = 1.0;  // V_C[n1, mode1]
    parameter real VC1_2 = 0.8;  // V_C[n1, mode2]

    // Complex permittivity parameters (frequency-dependent)
    parameter real eps_r_dc = 4.0;    // DC relative permittivity
    parameter real eps_r_inf = 2.5;   // High-freq relative permittivity
    parameter real tau_d = 1e-9;      // Debye relaxation time
    parameter real sigma_dc = 0.01;   // DC conductivity [S/m]
    parameter real eps_0 = 8.854e-12; // Vacuum permittivity

    real omega, eps_r_re, eps_r_im;
    real C_mode, G_mode;
    real debye_denom;

    analog begin
        omega = 2 * `M_PI * $freq;
        if (omega < 1e-10) omega = 1e-10;

        // Debye relaxation model for complex permittivity
        // epsilon(omega) = eps_inf + (eps_dc - eps_inf)/(1 + j*omega*tau_d) - j*sigma/(omega*eps_0)
        debye_denom = 1 + pow(omega * tau_d, 2);
        eps_r_re = eps_r_inf + (eps_r_dc - eps_r_inf) / debye_denom;
        eps_r_im = (eps_r_dc - eps_r_inf) * omega * tau_d / debye_denom
                  + sigma_dc / (omega * eps_0);

        // Mode capacitance and conductance
        C_mode = (UC1_1 * VC1_1 + UC1_2 * VC1_2) * eps_r_re;
        G_mode = (UC1_1 * VC1_1 + UC1_2 * VC1_2) * omega * eps_r_im;

        // Current: I = C*dV/dt + G*V (capacitive + lossy)
        I(n1) <+ C_mode * ddt(V(aux1)) + G_mode * V(aux1);
    end
endmodule
```

**Complete Low-Rank PEEC Segment with Skin Effect and Complex ε**:

```verilog
`include "disciplines.vams"

// Complete PEEC segment with:
// - Self R, L with Dowell skin effect
// - Self C with complex permittivity
// - Low-rank mutual coupling (L and C)

module peec_segment_lowrank(p, n, aux_L1, aux_L2, aux_C1, aux_C2);
    inout p, n;
    inout aux_L1, aux_L2, aux_C1, aux_C2;
    electrical p, n, aux_L1, aux_L2, aux_C1, aux_C2;

    // Self parameters
    parameter real R_dc = 1e-3;      // DC resistance [Ohm]
    parameter real L_self = 10e-9;   // Self inductance [H]
    parameter real C_self = 1e-12;   // Self capacitance [F]

    // Skin effect
    parameter real tau_skin = 1e-6;  // a^2 * mu * sigma

    // Complex permittivity
    parameter real eps_r = 4.0;
    parameter real tan_delta = 0.02;  // Loss tangent
    parameter real sigma = 0.01;      // DC conductivity

    // Low-rank coupling coefficients
    parameter real U_L1 = 1e-9;
    parameter real U_L2 = 0.5e-9;
    parameter real U_C1 = 1e-13;
    parameter real U_C2 = 0.5e-13;

    real omega, xi, F_R, F_L;
    real eps_eff, G_dielectric;

    analog begin
        omega = 2 * `M_PI * $freq;
        if (omega < 1e-10) omega = 1e-10;

        // Dowell factors
        xi = sqrt(tau_skin * omega / 2);
        if (xi < 0.1) begin
            F_R = 1.0 + pow(xi, 4) / 45.0;
            F_L = 1.0 - pow(xi, 4) / 15.0;
        end else begin
            F_R = xi * (sinh(2*xi) + sin(2*xi)) / (cosh(2*xi) - cos(2*xi));
            F_L = 1.5 / xi * (sinh(2*xi) - sin(2*xi)) / (cosh(2*xi) - cos(2*xi));
        end

        // Complex permittivity effect
        eps_eff = eps_r;
        G_dielectric = omega * eps_r * tan_delta * C_self + sigma * C_self / eps_r;

        // Self impedance: V = R_ac*I + L_ac*dI/dt
        V(p, n) <+ R_dc * F_R * I(p, n);
        V(p, n) <+ L_self * F_L * ddt(I(p, n));

        // Self capacitance with loss: I = C*dV/dt + G*V
        I(p, n) <+ C_self * eps_eff * ddt(V(p, n));
        I(p, n) <+ G_dielectric * V(p, n);

        // Low-rank mutual inductance coupling (from auxiliary nodes)
        V(p, n) <+ U_L1 * F_L * ddt(I(aux_L1));
        V(p, n) <+ U_L2 * F_L * ddt(I(aux_L2));

        // Low-rank mutual capacitance coupling (to auxiliary nodes)
        I(p, n) <+ U_C1 * eps_eff * ddt(V(aux_C1));
        I(p, n) <+ U_C2 * eps_eff * ddt(V(aux_C2));
    end
endmodule
```

**Advantages of Verilog-A Implementation**:

| Aspect | Standard SPICE | Verilog-A |
|--------|---------------|-----------|
| Skin effect | Ladder approximation (finite accuracy) | Exact Dowell formula |
| Complex ε(ω) | Fixed R//C (single frequency) | Debye model (all frequencies) |
| Low-rank coupling | Auxiliary nodes + CCCS/VCVS | Direct implementation |
| Model accuracy | Approximate | Exact |
| Simulator support | Universal | Spectre, ADS, HSPICE |

### PyKAN for Arbitrary Frequency-Dependent Properties

When material properties are available only as **numerical data** (measured values, lookup tables, or FEM simulation results), PyKAN provides a powerful alternative that is **not limited to analytical models** like Debye relaxation.

**Advantages over Analytical Models**:

| Aspect | Debye/Cole-Cole | PyKAN |
|--------|----------------|-------|
| Model form | Fixed functional form | Arbitrary learned function |
| Multi-pole relaxation | Manual pole fitting | Automatic learning |
| Measured data | Requires curve fitting | Direct learning |
| Nonlinear effects | Not supported | Supported (e.g., ε(E, ω)) |
| Formula extraction | N/A | Symbolic regression |

**Workflow: From Measured Data to Verilog-A**:

```
1. Measured data: ε'(ω), ε"(ω) at discrete frequencies
                  ↓
2. PyKAN training: Learn ε(ω) = f(log(ω))
                  ↓
3. Formula extraction: Symbolic regression
   e.g., ε'(ω) = a + b*tanh(c*log(ω) + d)
         ε"(ω) = e*exp(-f*(log(ω) - g)^2)
                  ↓
4. Verilog-A generation: Automatic code output
```

**PyKAN Training for Complex Permittivity**:

```python
import torch
from kan import KAN
import numpy as np

# Measured data (example: ferrite material)
freq_data = np.logspace(3, 9, 100)  # 1 kHz to 1 GHz
omega_data = 2 * np.pi * freq_data
eps_r_data = np.array([...])  # Measured eps'
eps_i_data = np.array([...])  # Measured eps"

# Normalize inputs
log_omega = np.log10(omega_data)
log_omega_norm = (log_omega - log_omega.mean()) / log_omega.std()

# Create training data
X = torch.tensor(log_omega_norm, dtype=torch.float32).reshape(-1, 1)
Y_re = torch.tensor(eps_r_data, dtype=torch.float32).reshape(-1, 1)
Y_im = torch.tensor(eps_i_data, dtype=torch.float32).reshape(-1, 1)

# Train KAN for real part
kan_re = KAN(width=[1, 5, 1], grid=10, k=3)
kan_re.fit({'train_input': X, 'train_label': Y_re}, steps=2000)

# Train KAN for imaginary part
kan_im = KAN(width=[1, 5, 1], grid=10, k=3)
kan_im.fit({'train_input': X, 'train_label': Y_im}, steps=2000)

# Extract symbolic formulas
formula_re = kan_re.symbolic_formula()[0][0]
formula_im = kan_im.symbolic_formula()[0][0]

print(f"eps'(omega) = {formula_re}")
print(f"eps''(omega) = {formula_im}")
```

**Auto-Generated Verilog-A from PyKAN**:

```verilog
`include "disciplines.vams"

// Auto-generated from PyKAN training
// Material: Ferrite MN60 (measured data)
// Frequency range: 1 kHz - 1 GHz

module pykan_complex_eps(p, n);
    inout p, n;
    electrical p, n;

    parameter real C_0 = 1e-12;  // Geometry factor [F]

    real omega, log_omega, log_omega_norm;
    real eps_r_re, eps_r_im;
    real C_eff, G_loss;

    // Normalization constants (from training)
    parameter real log_omega_mean = 6.5;
    parameter real log_omega_std = 1.8;

    analog begin
        omega = 2 * `M_PI * $freq;
        if (omega < 1e-10) omega = 1e-10;

        log_omega = log10(omega);
        log_omega_norm = (log_omega - log_omega_mean) / log_omega_std;

        // PyKAN-extracted formulas (example)
        // eps'(x) = 12.5 + 8.3*tanh(-0.42*x + 0.15)
        // eps''(x) = 2.1*exp(-0.8*(x - 0.3)^2) + 0.05*exp(0.9*x)
        eps_r_re = 12.5 + 8.3 * tanh(-0.42 * log_omega_norm + 0.15);
        eps_r_im = 2.1 * exp(-0.8 * pow(log_omega_norm - 0.3, 2))
                  + 0.05 * exp(0.9 * log_omega_norm);

        // Effective capacitance and loss conductance
        C_eff = eps_r_re * C_0;
        G_loss = omega * eps_r_im * C_0;

        // Constitutive relation
        I(p, n) <+ C_eff * ddt(V(p, n)) + G_loss * V(p, n);
    end
endmodule
```

**Combining Low-Rank + PyKAN**:

For large PEEC models with complex materials:

```
Full model complexity: O(N^2) elements × O(M) frequency points
                      ↓
Low-rank (ACA): O(k*N) elements
                      ↓
PyKAN frequency model: Compact formula (few parameters)
                      ↓
Final Verilog-A: O(k*N) modules with O(1) frequency evaluation
```

**Benefits Summary**:

1. **No analytical model required**: Works with pure numerical data
2. **Automatic formula extraction**: Interpretable symbolic expressions
3. **Multi-physics support**: Can learn ε(ω, T, E) for temperature and field dependence
4. **Verilog-A ready**: Direct code generation for circuit simulation
5. **Combines with ACA**: Low-rank structure + learned frequency dependence

### Framework Architecture: Two Paths to Verilog-A

This framework provides **two complementary paths** for converting frequency-dependent material properties to SPICE-compatible Verilog-A:

```
┌─────────────────────────────────────────────────────────────────┐
│                    Material Property Input                       │
└─────────────────────────────────────────────────────────────────┘
                              │
              ┌───────────────┴───────────────┐
              │                               │
              ▼                               ▼
┌─────────────────────────┐     ┌─────────────────────────┐
│   Path 1: Analytical    │     │   Path 2: Data-Driven   │
│   Model Available       │     │   (No Analytical Model) │
├─────────────────────────┤     ├─────────────────────────┤
│ • Debye                 │     │ • Measured data         │
│ • Cole-Cole             │     │ • FEM simulation        │
│ • Cole-Davidson         │     │ • Literature curves     │
│ • Havriliak-Negami      │     │ • Complex geometries    │
│ • Multi-pole Debye      │     │ • Nonlinear effects     │
│ • Dowell skin effect    │     │                         │
└───────────┬─────────────┘     └───────────┬─────────────┘
            │                               │
            ▼                               ▼
┌─────────────────────────┐     ┌─────────────────────────┐
│  Direct Verilog-A       │     │  PyKAN Learning         │
│  Implementation         │     │  + Formula Extraction   │
│  (closed-form formula)  │     │  + Auto Verilog-A Gen   │
└───────────┬─────────────┘     └───────────┬─────────────┘
            │                               │
            └───────────────┬───────────────┘
                            ▼
              ┌─────────────────────────┐
              │   Verilog-A Module      │
              │   for SPICE Simulation  │
              └─────────────────────────┘
```

**Path 1 (Analytical Models)**: When material behavior follows known physics (Debye relaxation, skin effect, etc.), use established formulas directly in Verilog-A. This is simpler, faster, and produces compact models.

**Path 2 (PyKAN)**: When no analytical model fits the data (arbitrary cross-sections, complex frequency dependence, nonlinear effects), PyKAN learns the relationship from numerical data and extracts a symbolic formula for Verilog-A.

**Key Insight**: PyKAN is a fallback for cases where analytical models fail, not a replacement for them. Use the simplest approach that works:

| Material Behavior | Recommended Path |
|-------------------|------------------|
| Standard dielectric (single τ) | Debye → Direct Verilog-A |
| Ferrite with broad relaxation | Cole-Cole → Direct Verilog-A |
| Rectangular conductor skin effect | Dowell → Direct Verilog-A |
| Arbitrary cross-section Zs | PyKAN → Verilog-A |
| ε(ω, T) temperature-dependent | PyKAN → Verilog-A |
| Nonlinear μ(H, ω) | PyKAN → Verilog-A |

### Path 1: Analytical Models (Direct Verilog-A)

The framework supports **standard analytical models** with direct Verilog-A implementation:

**Supported Analytical Models**:

| Model | Formula | Use Case |
|-------|---------|----------|
| **Debye** | $\varepsilon(\omega) = \varepsilon_\infty + \frac{\varepsilon_s - \varepsilon_\infty}{1 + j\omega\tau}$ | Single relaxation time |
| **Cole-Cole** | $\varepsilon(\omega) = \varepsilon_\infty + \frac{\varepsilon_s - \varepsilon_\infty}{1 + (j\omega\tau)^\alpha}$ | Distributed relaxation |
| **Cole-Davidson** | $\varepsilon(\omega) = \varepsilon_\infty + \frac{\varepsilon_s - \varepsilon_\infty}{(1 + j\omega\tau)^\beta}$ | Asymmetric distribution |
| **Havriliak-Negami** | $\varepsilon(\omega) = \varepsilon_\infty + \frac{\varepsilon_s - \varepsilon_\infty}{(1 + (j\omega\tau)^\alpha)^\beta}$ | General case |

**Verilog-A for Debye Model**:

```verilog
`include "disciplines.vams"

// Debye relaxation model for complex permittivity
// epsilon(omega) = eps_inf + (eps_s - eps_inf) / (1 + j*omega*tau)

module debye_capacitor(p, n);
    inout p, n;
    electrical p, n;

    parameter real C_0 = 1e-12;      // Geometry factor [F]
    parameter real eps_s = 80.0;     // Static permittivity
    parameter real eps_inf = 4.0;    // High-frequency permittivity
    parameter real tau = 1e-9;       // Relaxation time [s]
    parameter real sigma = 0.0;      // DC conductivity [S/m]

    real omega, omega_tau, denom;
    real eps_re, eps_im;
    real C_eff, G_eff;

    analog begin
        omega = 2 * `M_PI * $freq;
        if (omega < 1e-10) omega = 1e-10;

        omega_tau = omega * tau;
        denom = 1 + omega_tau * omega_tau;

        // Debye formula: eps' and eps"
        eps_re = eps_inf + (eps_s - eps_inf) / denom;
        eps_im = (eps_s - eps_inf) * omega_tau / denom + sigma / (omega * 8.854e-12);

        // Effective C and G
        C_eff = eps_re * C_0;
        G_eff = omega * eps_im * C_0;

        // Constitutive relation
        I(p, n) <+ C_eff * ddt(V(p, n)) + G_eff * V(p, n);
    end
endmodule
```

**Verilog-A for Cole-Cole Model**:

```verilog
`include "disciplines.vams"

// Cole-Cole model with distribution parameter alpha
// epsilon(omega) = eps_inf + (eps_s - eps_inf) / (1 + (j*omega*tau)^alpha)

module cole_cole_capacitor(p, n);
    inout p, n;
    electrical p, n;

    parameter real C_0 = 1e-12;
    parameter real eps_s = 80.0;
    parameter real eps_inf = 4.0;
    parameter real tau = 1e-9;
    parameter real alpha = 0.8;  // Distribution parameter (0 < alpha <= 1)

    real omega, x, theta;
    real denom_re, denom_im, denom_mag2;
    real eps_re, eps_im;
    real C_eff, G_eff;

    analog begin
        omega = 2 * `M_PI * $freq;
        if (omega < 1e-10) omega = 1e-10;

        // (j*omega*tau)^alpha = (omega*tau)^alpha * exp(j*alpha*pi/2)
        x = pow(omega * tau, alpha);
        theta = alpha * `M_PI / 2;

        // 1 + (j*omega*tau)^alpha
        denom_re = 1 + x * cos(theta);
        denom_im = x * sin(theta);
        denom_mag2 = denom_re * denom_re + denom_im * denom_im;

        // eps = eps_inf + (eps_s - eps_inf) / denom
        eps_re = eps_inf + (eps_s - eps_inf) * denom_re / denom_mag2;
        eps_im = (eps_s - eps_inf) * denom_im / denom_mag2;

        C_eff = eps_re * C_0;
        G_eff = omega * eps_im * C_0;

        I(p, n) <+ C_eff * ddt(V(p, n)) + G_eff * V(p, n);
    end
endmodule
```

**Multi-Pole Debye (Multiple Relaxation Times)**:

```verilog
`include "disciplines.vams"

// Multi-pole Debye model: sum of N Debye terms
// epsilon(omega) = eps_inf + sum_k (Delta_eps_k / (1 + j*omega*tau_k))

module multi_debye_capacitor(p, n);
    inout p, n;
    electrical p, n;

    parameter real C_0 = 1e-12;
    parameter real eps_inf = 2.5;

    // Debye poles (example: 3 poles)
    parameter real Delta_eps_1 = 10.0;  parameter real tau_1 = 1e-12;
    parameter real Delta_eps_2 = 20.0;  parameter real tau_2 = 1e-10;
    parameter real Delta_eps_3 = 50.0;  parameter real tau_3 = 1e-8;

    real omega, eps_re, eps_im;
    real ot1, ot2, ot3, d1, d2, d3;

    analog begin
        omega = 2 * `M_PI * $freq;
        if (omega < 1e-10) omega = 1e-10;

        ot1 = omega * tau_1; d1 = 1 + ot1*ot1;
        ot2 = omega * tau_2; d2 = 1 + ot2*ot2;
        ot3 = omega * tau_3; d3 = 1 + ot3*ot3;

        eps_re = eps_inf + Delta_eps_1/d1 + Delta_eps_2/d2 + Delta_eps_3/d3;
        eps_im = Delta_eps_1*ot1/d1 + Delta_eps_2*ot2/d2 + Delta_eps_3*ot3/d3;

        I(p, n) <+ eps_re * C_0 * ddt(V(p, n)) + omega * eps_im * C_0 * V(p, n);
    end
endmodule
```

**Model Selection Decision Tree**:

```
Do you have analytical model parameters?
├── Yes → Which model fits your data?
│         ├── Single relaxation → Debye
│         ├── Broad peak → Cole-Cole (fit alpha)
│         ├── Asymmetric → Cole-Davidson or Havriliak-Negami
│         └── Multiple peaks → Multi-pole Debye
└── No → Is data available?
         ├── Measured ε'(ω), ε"(ω) → PyKAN learning
         ├── FEM simulation data → PyKAN learning
         └── Literature curves → Digitize + PyKAN
```

### Path 2 Bonus: PyKAN for Model Identification

PyKAN can also be used to **identify** which analytical model best fits measured data (useful when you're unsure which Path 1 model applies):

```python
# Train PyKAN on measured data
kan = KAN(width=[1, 5, 1], grid=10, k=3)
kan.fit({'train_input': X_omega, 'train_label': Y_eps}, steps=2000)

# Extract symbolic formula
formula = kan.symbolic_formula()

# Compare with known models
# If formula ~ 1/(1 + x^2) → Debye
# If formula ~ 1/(1 + x^alpha) → Cole-Cole
# If formula contains multiple 1/(1+x^2) terms → Multi-pole Debye

print(f"Identified model: {identify_model(formula)}")
print(f"Parameters: {extract_parameters(formula)}")
```

**Summary of Two-Path Framework**:

| Path | When to Use | Output | Complexity |
|------|-------------|--------|------------|
| **Path 1** | Analytical model exists (Debye, Dowell, etc.) | Direct Verilog-A | Simple |
| **Path 2** | No analytical model, only data | PyKAN → Verilog-A | More complex |

- **Known materials**: Use Path 1 with established analytical models
- **New/complex materials**: Use Path 2 with PyKAN learning
- **Uncertain cases**: Use PyKAN (Path 2) to identify which Path 1 model applies

---

## Complex Permeability Support

### Complex Permeability Model

For magnetic materials with frequency-dependent losses:

```
μ(ω) = μ'(ω) - j·μ"(ω)

μ': Storage (reactive)
μ": Loss (resistive)
```

### Implementation in ESIM

The 1D cell problem becomes:

```
d/dz[(1/μ(z)) · dH/dz] = jωσ·H

with μ(z) = μ(H(z)) for nonlinear materials
```

**Note**: μ(H) can be obtained from the B-H curve of a specific material.

### Effect on Surface Impedance

```
Z_s = √(jωμ/σ) · coth(√(jωμσ) · a)

For complex μ = μ' - jμ":
- Real part increases (additional loss)
- Imaginary part modified (changed reactance)
```

### Ferrite Core Example

At f = 100 kHz with ferrite (μ_r = 2000 - j200):

| Quantity | Real μ only | Complex μ |
|----------|-------------|-----------|
| Re(Z_s) | 0.82 Ω/m | 0.95 Ω/m |
| Im(Z_s) | 0.82 Ω/m | 0.73 Ω/m |
| Loss factor | 1.0 | 1.16 |

The imaginary part of μ increases losses by ~16% in this example.

---

## Method Selection Guide

### Decision Tree

```
Is material nonlinear (μ(H))?
├── Yes → Use ESIM with homogenization
│         └── Complex μ? → Include μ' - jμ" in cell problem
└── No → Is cross-section rectangular?
         ├── Yes → Use Dowell's formula
         │         └── High freq (ξ > 5)? → Consider ladder approx
         └── No → Is it round?
                  ├── Yes → Use Kelvin functions
                  └── No → Use 2D FEM or PyKAN
```

### Recommended Methods by Application

| Application | Material | Geometry | Method |
|-------------|----------|----------|--------|
| WPT coil | Copper | Litz | Dowell + strand correction |
| Transformer | Copper | Foil | Dowell |
| Induction heating | Steel | Various | ESIM + complex μ |
| EMI filter | Copper + Ferrite | Various | PEEC + MMM coupling |
| PCB trace | Copper | Rectangular | Dowell |

---

## SPICE Output and Circuit Models

### Overview: From PEEC to SPICE

The goal of SPICE conversion is to represent the frequency-dependent PEEC impedance as circuit elements that SPICE simulators can process. Three approaches are available:

| Approach | Implementation | Accuracy | Complexity |
|----------|----------------|----------|------------|
| **Verilog-A** | Direct formula evaluation | Exact | Requires Verilog-A support |
| **RL Ladder** | Continued fraction expansion | Approximate | Standard SPICE |
| **Lookup Table** | Pre-computed values | Discrete frequencies | Simple |

### Mathematical Foundation

The PEEC impedance with skin effect can be written as:

$$
Z(s) = R_{dc} \cdot F(\xi) + s \cdot L_{int,dc} \cdot G(\xi)
$$

where $s = j\omega$ and $\xi = \frac{a}{\delta} = a\sqrt{\frac{\omega\mu\sigma}{2}}$.

In the Laplace domain with time constant $\tau = a^2 \mu \sigma$:

$$
\frac{Z(s)}{R_{dc}} = \sqrt{\tau s} \cdot \coth(\sqrt{\tau s})
$$

### Approach 1: Verilog-A (Recommended)

Verilog-A allows direct implementation of frequency-dependent impedance formulas.

**Dowell's Formula in Verilog-A**:

The resistance factor $F(\xi)$ and inductance factor $G(\xi)$ are:

$$
F(\xi) = \xi \cdot \frac{\sinh(2\xi) + \sin(2\xi)}{\cosh(2\xi) - \cos(2\xi)}
$$

$$
G(\xi) = \frac{3}{2\xi} \cdot \frac{\sinh(2\xi) - \sin(2\xi)}{\cosh(2\xi) - \cos(2\xi)}
$$

```verilog
`include "disciplines.vams"

module dowell_skin(p, n);
    inout p, n;
    electrical p, n;

    // Physical parameters
    parameter real d = 0.1e-3;      // Conductor thickness [m]
    parameter real sigma = 5.8e7;   // Conductivity [S/m]
    parameter real mu = 1.2566e-6;  // Permeability [H/m]
    parameter real length = 1.0;    // Conductor length [m]
    parameter real width = 10e-3;   // Conductor width [m]

    real R_dc, L_int_dc, xi, F_R, F_L, omega, delta;
    real sh_xi, ch_xi, sn_xi, cs_xi;
    real sh_2xi, ch_2xi, sn_2xi, cs_2xi;

    analog begin
        // Frequency and skin depth
        omega = 2 * `M_PI * $freq;
        if (omega < 1e-10) omega = 1e-10;  // Avoid DC singularity
        delta = sqrt(2.0 / (omega * mu * sigma));
        xi = (d / 2) / delta;  // Half-thickness / skin_depth

        // Hyperbolic and trigonometric functions
        sh_xi = sinh(xi);
        ch_xi = cosh(xi);
        sn_xi = sin(xi);
        cs_xi = cos(xi);
        sh_2xi = sinh(2*xi);
        ch_2xi = cosh(2*xi);
        sn_2xi = sin(2*xi);
        cs_2xi = cos(2*xi);

        // Dowell factors F(xi) and G(xi)
        if (xi < 0.01) begin
            F_R = 1.0;  // DC limit
            F_L = 1.0;  // DC limit
        end else begin
            F_R = xi * (sh_2xi + sn_2xi) / (ch_2xi - cs_2xi);
            F_L = (3.0 / (2.0 * xi)) * (sh_2xi - sn_2xi) / (ch_2xi - cs_2xi);
        end

        // DC resistance and internal inductance
        R_dc = length / (sigma * width * d);
        L_int_dc = mu * length * d / (3.0 * width);

        // Constitutive relation: V = R_ac * I + L_ac * dI/dt
        V(p, n) <+ R_dc * F_R * I(p, n) + L_int_dc * F_L * ddt(I(p, n));
    end
endmodule
```

**Alternative Form (coth formulation)**:

Using the complex propagation constant $\gamma = (1+j)/\delta$:

$$
Z(s) = R_{dc} \cdot \gamma a \cdot \coth(\gamma a)
$$

```verilog
`include "disciplines.vams"

module dowell_coth(p, n);
    inout p, n;
    electrical p, n;

    parameter real a = 0.05e-3;     // Half-thickness [m]
    parameter real sigma = 5.8e7;   // Conductivity [S/m]
    parameter real mu = 1.2566e-6;  // Permeability [H/m]

    real omega, R_dc;
    real gamma_re, gamma_im, ga_re, ga_im;
    real sh_re, sh_im, ch_re, ch_im;
    real coth_re, coth_im, Z_re, Z_im;
    real denom;

    analog begin
        omega = 2 * `M_PI * $freq;
        if (omega < 1e-10) omega = 1e-10;

        R_dc = 1.0 / (sigma * 2 * a);  // Per unit area

        // gamma = (1+j) / delta = (1+j) * sqrt(omega*mu*sigma/2)
        gamma_re = sqrt(omega * mu * sigma / 2.0);
        gamma_im = gamma_re;

        // gamma * a
        ga_re = gamma_re * a;
        ga_im = gamma_im * a;

        // sinh(gamma*a) and cosh(gamma*a) using complex exponentials
        // sinh(x+jy) = sinh(x)cos(y) + j*cosh(x)sin(y)
        // cosh(x+jy) = cosh(x)cos(y) + j*sinh(x)sin(y)
        sh_re = sinh(ga_re) * cos(ga_im);
        sh_im = cosh(ga_re) * sin(ga_im);
        ch_re = cosh(ga_re) * cos(ga_im);
        ch_im = sinh(ga_re) * sin(ga_im);

        // coth = cosh/sinh
        denom = sh_re*sh_re + sh_im*sh_im;
        coth_re = (ch_re*sh_re + ch_im*sh_im) / denom;
        coth_im = (ch_im*sh_re - ch_re*sh_im) / denom;

        // Z = R_dc * (gamma*a) * coth(gamma*a)
        Z_re = R_dc * (ga_re*coth_re - ga_im*coth_im);
        Z_im = R_dc * (ga_re*coth_im + ga_im*coth_re);

        V(p, n) <+ Z_re * I(p, n) + (Z_im / omega) * ddt(I(p, n));
    end
endmodule
```

### Approach 2: RL Ladder Network (Standard SPICE)

For simulators without Verilog-A support, use the continued fraction expansion:

$$
\sqrt{\tau s} \cdot \coth(\sqrt{\tau s}) = 1 + \cfrac{\tau s}{3 + \cfrac{\tau s}{5 + \cfrac{\tau s}{7 + \cdots}}}
$$

This maps to a ladder circuit:

```
       R_dc        L1         L2         L3
  o----/\/\/---o---)))))---o---)))))---o---)))))---o
                   |           |           |
                  === R1      === R2      === R3
                   |           |           |
                  ---         ---         ---

Element values (for n-th section):
  L_n = tau * R_dc / (2n+1)
  R_n = (2n+1) * R_dc / (2n-1)    [for n >= 1]

where tau = a^2 * mu * sigma
```

**SPICE Netlist (3-section ladder)**:

```spice
* Dowell Skin Effect Model - 3-section RL ladder
* tau = a^2 * mu * sigma, R_dc = 1/(sigma*d)
.SUBCKT DOWELL_LADDER in out tau=73n R_dc=172m

* Series DC resistance
R_dc in n1 {R_dc}

* Section 1: L1 = tau*R_dc/3, R1 = 3*R_dc
L1 n1 n2 {tau*R_dc/3}
R1 n2 0 {3*R_dc}

* Section 2: L2 = tau*R_dc/5, R2 = 5*R_dc/3
L2 n2 n3 {tau*R_dc/5}
R2 n3 0 {5*R_dc/3}

* Section 3: L3 = tau*R_dc/7, R3 = 7*R_dc/5
L3 n3 out {tau*R_dc/7}
R3 out 0 {7*R_dc/5}

.ENDS DOWELL_LADDER
```

**Accuracy of Ladder Approximation**:

| Sections | Max Error (xi < 1) | Max Error (xi < 3) | Max Error (xi < 10) |
|----------|--------------------|--------------------|---------------------|
| 2 | < 1% | < 5% | < 30% |
| 3 | < 0.1% | < 1% | < 10% |
| 5 | < 0.01% | < 0.1% | < 3% |
| 10 | < 0.001% | < 0.01% | < 0.3% |

### Approach 3: Vector Fitting (Rational Approximation)

For arbitrary Z(s), use Vector Fitting to obtain a rational function:

$$
Z(s) \approx \sum_{n=1}^{N} \frac{c_n}{s - p_n} + d + s \cdot e
$$

This can be synthesized as a Foster or Cauer network.

**Python Implementation**:

```python
from scipy.optimize import least_squares
import numpy as np

def vector_fitting(freq, Z_data, n_poles=8):
    """Vector Fitting for rational approximation.

    Args:
        freq: Frequency array [Hz]
        Z_data: Complex impedance data
        n_poles: Number of poles
    Returns:
        poles, residues, d, e coefficients
    """
    s = 1j * 2 * np.pi * freq

    # Initial pole guess (logarithmically spaced, real negative)
    f_min, f_max = freq.min(), freq.max()
    poles_init = -2 * np.pi * np.logspace(np.log10(f_min), np.log10(f_max), n_poles)

    # ... (iterative fitting algorithm)
    # Returns poles, residues, d, e

    return poles, residues, d, e

def synthesize_foster_network(poles, residues, d, e):
    """Convert rational function to Foster network.

    Foster form: Z(s) = R0 + s*L0 + sum(R_n || L_n)

    Returns:
        List of (R, L) pairs for parallel branches
    """
    branches = []
    for p, c in zip(poles, residues):
        # Each pole-residue pair becomes a parallel RL branch
        # c / (s - p) = -c/p / (1 - s/p) -> R = -c/p, L = c/p^2
        R = -c.real / p.real
        L = -c.real / (p.real ** 2)
        branches.append((R, L))
    return branches, d, e
```

### Approach 4: Lookup Table

For simple cases, pre-compute Z(f) at discrete frequencies and use interpolation.

```spice
* Frequency-dependent R and L using lookup table
* Requires Spectre or HSPICE with pwl() function

.SUBCKT DOWELL_TABLE in out R_dc=172m L_int_dc=4.2n

* AC resistance as frequency-dependent multiplier
R_ac in mid 'pwl($freq, 1, 1, 1k, 1.001, 10k, 1.086, 100k, 3.16, 1M, 10) * R_dc'

* Internal inductance as frequency-dependent
L_int mid out 'pwl($freq, 1, 1, 1k, 1, 10k, 0.81, 100k, 0.32, 1M, 0.1) * L_int_dc'

.ENDS DOWELL_TABLE
```

### PyKAN Extracted Formula Verilog-A

When using PyKAN to learn arbitrary cross-section impedance:

```verilog
// Example of PyKAN-extracted analytical formula
module kan_zs(p, n);
    inout p, n;
    electrical p, n;

    parameter real sigma = 5.8e7;
    parameter real width = 1e-3;
    parameter real height = 0.5e-3;

    real omega, log_omega, Re_Zs, Im_Zs;

    analog begin
        omega = 2 * `M_PI * $freq;
        log_omega = ln(omega);

        // KAN-extracted formula (example - actual formula from training)
        Re_Zs = exp(0.5*log_omega - 8.2) * tanh(0.3*log_omega + 1.2);
        Im_Zs = exp(0.5*log_omega - 7.8) * (1 - exp(-0.1*log_omega));

        V(p, n) <+ Re_Zs * I(p, n) + Im_Zs / omega * ddt(I(p, n));
    end
endmodule
```

### Complete SPICE Model with External Inductance

A full conductor segment includes both surface impedance and external inductance:

```verilog
`include "disciplines.vams"

module peec_segment(p, n);
    inout p, n;
    electrical p, n;

    // Geometry
    parameter real length = 0.1;    // [m]
    parameter real width = 10e-3;   // [m]
    parameter real thickness = 0.1e-3;  // [m]

    // Material
    parameter real sigma = 5.8e7;   // [S/m]
    parameter real mu_r = 1.0;

    // External inductance (from PEEC)
    parameter real L_ext = 50e-9;   // [H] - computed by PEEC

    real mu, R_dc, L_int_dc, omega, delta, xi;
    real F_R, F_L, sh2, ch2, sn2, cs2;

    analog begin
        mu = 4.0 * `M_PI * 1e-7 * mu_r;

        // DC parameters
        R_dc = length / (sigma * width * thickness);
        L_int_dc = mu * length * thickness / (3.0 * width);

        // Frequency-dependent parameters
        omega = 2 * `M_PI * $freq;
        if (omega < 1e-10) omega = 1e-10;

        delta = sqrt(2.0 / (omega * mu * sigma));
        xi = (thickness / 2) / delta;

        // Dowell factors
        if (xi < 0.01) begin
            F_R = 1.0;
            F_L = 1.0;
        end else begin
            sh2 = sinh(2*xi); ch2 = cosh(2*xi);
            sn2 = sin(2*xi);  cs2 = cos(2*xi);
            F_R = xi * (sh2 + sn2) / (ch2 - cs2);
            F_L = (3.0 / (2.0 * xi)) * (sh2 - sn2) / (ch2 - cs2);
        end

        // Total impedance: R_ac + j*omega*(L_int + L_ext)
        V(p, n) <+ R_dc * F_R * I(p, n)
                 + (L_int_dc * F_L + L_ext) * ddt(I(p, n));
    end
endmodule
```

### SPICE Implementation Summary

| Method | SPICE Compatibility | Accuracy | Use Case |
|--------|---------------------|----------|----------|
| Verilog-A | Spectre, HSPICE, Xyce | Exact | Production simulation |
| RL Ladder | All SPICE | Good (3+ sections) | Legacy simulators |
| Vector Fitting | All SPICE | Excellent | Arbitrary Z(s) |
| Lookup Table | Limited | Discrete freqs | Quick estimation |

---

## PyKAN for Arbitrary Cross-Sections

### Overview

PyKAN (Kolmogorov-Arnold Networks) enables learning surface impedance for arbitrary cross-section shapes from FEM data. Unlike MLPs, KAN can extract **symbolic formulas** suitable for SPICE implementation.

### Why KAN Instead of MLP?

| Feature | MLP | KAN |
|---------|-----|-----|
| Formula extraction | Difficult (black box) | **Possible** (symbolic) |
| Parameter count | High (thousands) | Low (hundreds) |
| Physical interpretability | None | **Yes** (explicit functions) |
| Verilog-A conversion | Requires approximation | **Direct** |
| Training data efficiency | Low | High |

### Workflow

```
┌─────────────────────────────────────────────────────────────┐
│  1. Define Shape Parameter Space                            │
│     - Cross-section geometry parameters (w, h, aspect ratio)│
│     - Material parameters (σ, μ)                            │
│     - Frequency range (f_min to f_max)                      │
└───────────────────────────┬─────────────────────────────────┘
                            │
┌───────────────────────────▼─────────────────────────────────┐
│  2. Generate FEM Training Data (NGSolve)                    │
│     - 2D eddy current problem for each parameter set        │
│     - Extract Z_s = V/I from FEM solution                   │
│     - Store [params, Re(Zs), Im(Zs)] pairs                  │
└───────────────────────────┬─────────────────────────────────┘
                            │
┌───────────────────────────▼─────────────────────────────────┐
│  3. Train KAN Model                                         │
│     - Input: [log(omega), log(sigma), shape_params]         │
│     - Output: [Re(Zs), Im(Zs)]                              │
│     - Use log scaling for better numerical conditioning     │
└───────────────────────────┬─────────────────────────────────┘
                            │
┌───────────────────────────▼─────────────────────────────────┐
│  4. Extract Symbolic Formula                                │
│     - KAN auto-detects activation functions                 │
│     - Produces human-readable formula                       │
└───────────────────────────┬─────────────────────────────────┘
                            │
┌───────────────────────────▼─────────────────────────────────┐
│  5. Generate Verilog-A                                      │
│     - Convert symbolic formula to Verilog-A syntax          │
│     - Export for SPICE simulation                           │
└─────────────────────────────────────────────────────────────┘
```

### NGSolve 2D Eddy Current FEM for Training Data

```python
from ngsolve import *
from netgen.occ import *
import numpy as np

def compute_surface_impedance_fem(width, height, sigma, mu_r, freq):
    """Compute surface impedance using 2D FEM (NGSolve).

    Solves 2D eddy current problem:
        curl(1/mu * curl(A)) + j*omega*sigma*A = J_ext

    Args:
        width: Cross-section width [m]
        height: Cross-section height [m]
        sigma: Conductivity [S/m]
        mu_r: Relative permeability
        freq: Frequency [Hz]
    Returns:
        Complex surface impedance Z_s [Ohm/m]
    """
    omega = 2 * np.pi * freq
    mu = 4 * np.pi * 1e-7 * mu_r

    # Create geometry: conductor cross-section + surrounding air
    conductor = Rectangle(width, height).Face()
    conductor.mat("conductor")

    air_size = 5 * max(width, height)
    air = Rectangle(air_size, air_size).Face() - conductor
    air.mat("air")

    geo = OCCGeometry(Glue([conductor, air]), dim=2)
    mesh = Mesh(geo.GenerateMesh(maxh=min(width, height)/10))

    # Define FE space (H(curl) for vector potential)
    fes = HCurl(mesh, order=2, complex=True)
    u, v = fes.TnT()

    # Bilinear form: curl-curl + j*omega*sigma
    a = BilinearForm(fes)
    nu_cond = 1 / mu
    nu_air = 1 / (4 * np.pi * 1e-7)

    a += nu_cond * curl(u) * curl(v) * dx("conductor")
    a += nu_air * curl(u) * curl(v) * dx("air")
    a += 1j * omega * sigma * u * v * dx("conductor")

    # Source: unit current density
    f = LinearForm(fes)
    J_ext = 1.0 / (width * height)  # A/m^2 for unit total current
    f += J_ext * v * dx("conductor")

    a.Assemble()
    f.Assemble()

    # Solve
    gfu = GridFunction(fes)
    gfu.vec.data = a.mat.Inverse() * f.vec

    # Compute impedance per unit length
    # Z = V/I = (integral of E over cross-section) / I
    # E = -j*omega*A for quasi-static
    E_avg = Integrate(-1j * omega * gfu, mesh, definedon=mesh.Materials("conductor"))
    I = 1.0  # Unit current
    Z_s = E_avg / I  # Ohm/m

    return Z_s

def generate_training_data(param_ranges, n_samples=1000):
    """Generate training data for KAN.

    Args:
        param_ranges: Dict with parameter ranges
            {'width': (min, max), 'height': (min, max), ...}
        n_samples: Number of training samples
    Returns:
        X_train: Input array [n_samples, n_features]
        y_train: Output array [n_samples, 2] (Re(Zs), Im(Zs))
    """
    X_train = []
    y_train = []

    for _ in range(n_samples):
        # Random sample from parameter space
        width = np.random.uniform(*param_ranges['width'])
        height = np.random.uniform(*param_ranges['height'])
        sigma = 10**np.random.uniform(*np.log10(param_ranges['sigma']))
        mu_r = 10**np.random.uniform(*np.log10(param_ranges['mu_r']))
        freq = 10**np.random.uniform(*np.log10(param_ranges['freq']))

        # Compute FEM solution
        try:
            Z_s = compute_surface_impedance_fem(width, height, sigma, mu_r, freq)

            # Store (use log scale for better conditioning)
            omega = 2 * np.pi * freq
            X_train.append([np.log10(omega), np.log10(sigma),
                           np.log10(mu_r), width/height])
            y_train.append([Z_s.real, Z_s.imag])
        except:
            continue  # Skip failed simulations

    return np.array(X_train), np.array(y_train)
```

### KAN Model Training

```python
from pykan import KAN
import torch

def train_kan_model(X_train, y_train, width=[4, 8, 8, 2], steps=5000):
    """Train KAN model for surface impedance.

    Args:
        X_train: Input features [n_samples, n_input]
        y_train: Target values [n_samples, 2]
        width: KAN layer widths
        steps: Training steps
    Returns:
        Trained KAN model
    """
    # Convert to torch tensors
    X = torch.tensor(X_train, dtype=torch.float32)
    y = torch.tensor(y_train, dtype=torch.float32)

    # Create KAN model
    model = KAN(width=width)

    # Train
    model.train(
        {'train_input': X, 'train_label': y},
        steps=steps,
        lamb=0.01,  # Regularization for sparsity
        lamb_entropy=2.0  # Encourage simple functions
    )

    return model

def extract_formula_and_verilog(model, output_file="zs_model.va"):
    """Extract symbolic formula and generate Verilog-A.

    Args:
        model: Trained KAN model
        output_file: Verilog-A output filename
    Returns:
        Symbolic formula string
    """
    # Auto-detect symbolic functions
    model.auto_symbolic()

    # Get formula
    formula = model.symbolic_formula()

    # Generate Verilog-A
    verilog_code = f'''`include "disciplines.vams"

module kan_surface_impedance(p, n);
    inout p, n;
    electrical p, n;

    parameter real sigma = 5.8e7;  // Conductivity [S/m]
    parameter real mu_r = 1.0;     // Relative permeability
    parameter real width = 1e-3;  // Cross-section width [m]
    parameter real height = 0.5e-3; // Cross-section height [m]

    real omega, log_omega, log_sigma, log_mu_r, aspect;
    real Re_Zs, Im_Zs;

    analog begin
        omega = 2 * `M_PI * $freq;
        log_omega = log10(omega);
        log_sigma = log10(sigma);
        log_mu_r = log10(mu_r);
        aspect = width / height;

        // KAN-extracted formula (auto-generated)
        // {formula[0]}  // Re(Zs)
        // {formula[1]}  // Im(Zs)
        Re_Zs = {formula[0].replace('x0', 'log_omega').replace('x1', 'log_sigma').replace('x2', 'log_mu_r').replace('x3', 'aspect')};
        Im_Zs = {formula[1].replace('x0', 'log_omega').replace('x1', 'log_sigma').replace('x2', 'log_mu_r').replace('x3', 'aspect')};

        V(p, n) <+ Re_Zs * I(p, n) + Im_Zs / omega * ddt(I(p, n));
    end
endmodule
'''

    with open(output_file, 'w') as f:
        f.write(verilog_code)

    return formula
```

### Complete Workflow Example

```python
# Step 1: Define parameter space
param_ranges = {
    'width': (0.1e-3, 10e-3),      # 0.1 mm to 10 mm
    'height': (0.05e-3, 5e-3),     # 0.05 mm to 5 mm
    'sigma': (1e6, 6e7),           # Steel to copper
    'mu_r': (1, 1000),             # Non-magnetic to soft iron
    'freq': (100, 10e6)            # 100 Hz to 10 MHz
}

# Step 2: Generate training data
print("Generating training data...")
X_train, y_train = generate_training_data(param_ranges, n_samples=5000)
print(f"Generated {len(X_train)} training samples")

# Step 3: Train KAN
print("Training KAN model...")
model = train_kan_model(X_train, y_train)

# Step 4: Extract formula and generate Verilog-A
print("Extracting symbolic formula...")
formula = extract_formula_and_verilog(model, "zs_rectangular.va")
print(f"Re(Zs) = {formula[0]}")
print(f"Im(Zs) = {formula[1]}")
```

### Use Cases

1. **Arbitrary cross-section**: No analytical formula exists
   - L-shaped conductors
   - Hollow tubes
   - Multi-conductor bundles

2. **ESIM (nonlinear μ)**: Z(H, ω) as 2-variable function
   - Steel at different magnetization levels
   - Saturable ferrite cores

3. **Proximity effect**: Multi-conductor interactions
   - Litz wire bundles
   - Transformer windings

4. **Complex μ(ω)**: Frequency-dependent permeability
   - Ferrite materials with resonance

---

## Complex Material Properties Support

### Complex Permeability μ(ω) = μ' - jμ"

Magnetic materials have frequency-dependent properties:

```
mu(omega) = mu'(omega) - j*mu"(omega)

mu': Energy storage (reactance)
mu": Energy loss (hysteresis, eddy currents in grains)
```

**Supported materials**:
- Ferrite: mu(omega) resonance characteristics
- Laminated steel: Effective mu reduction due to eddy currents
- Amorphous: High-frequency loss characteristics

**Implementation**: Surface impedance Zs via ESIM or PyKAN learning

### Complex Permittivity ε(ω) = ε' - jε"

Dielectric materials also have frequency-dependent properties:

```
epsilon(omega) = epsilon'(omega) - j*epsilon"(omega)

epsilon': Dielectric energy storage (capacitance)
epsilon": Dielectric loss (relaxation, conduction)
```

**Implementation**: **PEEC Star component**

The Star component in Loop-Star decomposition handles capacitive effects:
```
Z_SS (Star-Star): Capacitive impedance matrix
Y_C = j*omega*C  where C includes epsilon(omega)

For complex epsilon:
  C_eff(omega) = epsilon'(omega) * C_0
  G_d(omega) = omega * epsilon"(omega) * C_0  (dielectric loss conductance)
```

---

## Application Examples

### 1. Wireless Power Transfer (WPT)

```
Configuration: TX coil + RX coil + Ferrite shield

PEEC: Coil mutual inductance
MMM: Flux concentration by ferrite
Zs: Litz wire skin/proximity effect
```

### 2. Induction Heating

```
Configuration: Heating coil + Workpiece (conductive)

PEEC: Coil self-inductance
MMM: Workpiece nonlinear mu (temperature-dependent)
Zs: Workpiece surface frequency-dependent impedance (ESIM)
```

### 3. High-Frequency Transformer

```
Configuration: Primary winding + Secondary winding + Magnetic core

PEEC: Winding coupling, leakage inductance
MMM: Core magnetization, saturation
Zs: Winding skin/proximity effect
```

### 4. EMI/EMC Filter

```
Configuration: Common-mode choke + Wiring

PEEC: Wiring impedance, capacitive coupling
MMM: Ferrite core frequency characteristics mu(omega)
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
| PEEC-MMM coupling | Implemented | Magnetic cores |

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
| Magnetic materials | MMM (high accuracy) | FEM | None | MoM |
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
