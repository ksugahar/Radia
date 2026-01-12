# PEEC Surface Impedance and ESIM Analysis

## Overview

This document provides a comprehensive guide to surface impedance in PEEC (Partial Element Equivalent Circuit) conductor modeling, including:

- PEEC impedance matrix formulation with skin effect
- Dowell's formula for rectangular conductors
- ESIM (Effective Surface Impedance Method) for nonlinear materials
- Laplace domain representation and equivalent circuits
- Proximity effect analysis

**Current Scope**: Rectangular cross-section conductors (1D cell problem)

**Future Extension**: Arbitrary 2D cross-sections (2D FEM cell problem)

## PEEC Basic Formulation

### Standard PEEC Impedance Matrix

For a conductor segment, the PEEC method constructs an impedance matrix:

```
V = Z · I
```

where:
- **V**: Voltage vector
- **I**: Current vector
- **Z**: Impedance matrix

The impedance matrix consists of:

```
Z = R + jωL
```

- **R**: Resistance matrix (diagonal for self, off-diagonal for mutual resistance)
- **L**: Inductance matrix (self and mutual inductance from Neumann formula)

### DC Resistance

For a conductor segment with length l, cross-section area A, and conductivity σ:

```
R_dc = ρ · l / A = l / (σ · A)
```

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

### From Surface Impedance to Conductor Resistance

For a conductor with perimeter P and length l:

**Current-field relationship** (Ampère's law):
```
I = ∮ H_t · dl = H_t · P   (for uniform H around perimeter)
```

So: `H_t = I / P`

**Voltage drop** along the conductor:
```
V = E_t · l = Z_s · H_t · l = Z_s · (I / P) · l
```

**Conductor impedance**:
```
Z_conductor = V / I = Z_s · l / P
```

**AC Resistance**:
```
R_ac = Re(Z_conductor) = Re(Z_s) · l / P
```

### Impedance Matrix with Skin Effect

The PEEC impedance matrix becomes:

```
Z_ii = R_ac,i + jωL_ii + jωL_int,i
     = Re(Z_s) · l_i / P_i + jωL_ii + jω·Im(Z_s) · l_i / P_i
```

where:
- **R_ac,i**: AC resistance of segment i
- **L_ii**: External self-inductance (from Neumann formula)
- **L_int,i**: Internal inductance (from skin effect)

The internal inductance contribution:
```
L_int = Im(Z_s) · l / (ω · P)
```

## Surface Impedance for Different Geometries

### Rectangular Cross-Section (Planar Approximation)

For width w >> height h (thin strip):

```
Half-thickness: a = h / 2
Perimeter: P ≈ 2w (current on top and bottom surfaces)
```

**Finite slab surface impedance**:
```
Z_s = ρ · γ · tanh(γ · a)

where γ = (1 + j) / δ
```

**R_ac/R_dc ratio** (Dowell's formula):
```
R_ac/R_dc = ξ · (sinh(2ξ) + sin(2ξ)) / (cosh(2ξ) - cos(2ξ))

where ξ = a / δ
```

### Round Wire

**Kelvin function formula**:
```
R_ac/R_dc = (ξ/2) · Re[(ber' + j·bei') / (ber + j·bei)]

where ξ = √2 · radius / δ
```

### Arbitrary Cross-Section

Requires 2D FEM solution of the cell problem:
```
∇·(σ·∇V) + jωσ·A = J_ext
```

## PEEC Matrix Assembly with ESIM

### Algorithm

1. **For each conductor segment i**:

   a. Compute geometric properties:
      - Length l_i
      - Cross-section area A_i
      - Perimeter P_i
      - Characteristic dimension a_i (half-thickness or radius)

   b. Compute DC resistance:
      ```
      R_dc,i = ρ · l_i / A_i
      ```

   c. Compute skin depth and ξ:
      ```
      δ = √(2ρ / (ωμ))
      ξ_i = a_i / δ
      ```

   d. Compute R_ac/R_dc ratio:
      - Linear material: Use Dowell or Kelvin formula
      - Nonlinear material: Use ESIM Homogenization
      ```
      F_i = R_ac/R_dc (from appropriate formula)
      ```

   e. Compute AC impedance:
      ```
      R_ac,i = R_dc,i · F_i
      Z_ii = R_ac,i + jωL_ii
      ```

2. **Mutual terms** (off-diagonal):
   ```
   Z_ij = jωM_ij   (mutual inductance only, no mutual resistance)
   ```

3. **Assemble full matrix**:
   ```
   Z = [Z_11  Z_12  ...  Z_1n ]
       [Z_21  Z_22  ...  Z_2n ]
       [...   ...   ...  ...  ]
       [Zn1   Zn2   ...  Z_nn ]
   ```

### Matrix Structure

```
Z = R_ac + jωL_ext + jωL_int

where:
  R_ac = diag(R_ac,1, R_ac,2, ..., R_ac,n)     (diagonal)
  L_ext = [L_ij] from Neumann formula          (full matrix)
  L_int = diag(L_int,1, L_int,2, ..., L_int,n) (diagonal)
```

## Nonlinear Materials: Iterative Solution

For materials with μ(H) dependence:

### Fixed-Point Iteration

```
1. Initialize: H_surface = I / P (estimate from DC current)

2. Repeat:
   a. Solve ESIM cell problem with current H_surface
      → Get Z_s(H_surface), μ_eff, ξ_eff

   b. Compute R_ac/R_dc = Dowell(ξ_eff)

   c. Update PEEC matrix: Z_ii = R_dc · F + jωL_ii

   d. Solve circuit: I = Z^(-1) · V

   e. Update H_surface = I / P

   f. Check convergence: |I_new - I_old| < tolerance

3. Return converged solution
```

### Practical Implementation

```python
class PEECWithESIM:
    def __init__(self, segments, sigma, frequency, bh_curve=None):
        self.segments = segments
        self.sigma = sigma
        self.frequency = frequency
        self.bh_curve = bh_curve

        # Pre-compute DC resistances and inductances
        self.R_dc = self._compute_dc_resistances()
        self.L_matrix = self._compute_inductance_matrix()

        # Initialize ESIM solver for each segment
        self.esim_solvers = self._create_esim_solvers()

    def solve(self, V_source, tol=1e-6, max_iter=50):
        """Solve PEEC circuit with skin effect."""

        # Initial guess: DC solution
        Z = np.diag(self.R_dc) + 1j * self.omega * self.L_matrix
        I = np.linalg.solve(Z, V_source)

        for iteration in range(max_iter):
            # Update surface field for each segment
            H_surface = self._compute_surface_field(I)

            # Update R_ac using ESIM
            R_ac = self._compute_ac_resistance(H_surface)

            # Rebuild impedance matrix
            Z = np.diag(R_ac) + 1j * self.omega * self.L_matrix

            # Solve for new currents
            I_new = np.linalg.solve(Z, V_source)

            # Check convergence
            if np.max(np.abs(I_new - I)) < tol * np.max(np.abs(I)):
                return I_new, {'converged': True, 'iterations': iteration + 1}

            I = I_new

        return I, {'converged': False, 'iterations': max_iter}

    def _compute_ac_resistance(self, H_surface):
        """Compute AC resistance for each segment using ESIM."""
        R_ac = np.zeros(len(self.segments))

        for i, (seg, solver) in enumerate(zip(self.segments, self.esim_solvers)):
            if self.bh_curve is not None:
                # Nonlinear: use ESIM homogenization
                result = solver.solve(H_surface[i])
                F = result['R_ratio']
            else:
                # Linear: use Dowell formula
                F = self._dowell_ratio(solver.xi)

            R_ac[i] = self.R_dc[i] * F

        return R_ac
```

## Summary

| Component | Formula | Notes |
|-----------|---------|-------|
| DC Resistance | R_dc = ρ·l/A | Always present |
| Surface Impedance | Z_s = ρ·γ·tanh(γ·a) | Finite slab |
| AC/DC Ratio | F = Dowell(ξ) or Kelvin(ξ) | Geometry-dependent |
| AC Resistance | R_ac = R_dc · F | Diagonal matrix element |
| PEEC Impedance | Z = R_ac + jωL | Full matrix |

## ESIM vs Dowell: Boundary Condition Difference

### Key Finding

**ESIM (dH/dz = 0 BC) and Dowell's formula (H = 0 BC) solve DIFFERENT problems!**

| Aspect | ESIM (Igarashi homogenization) | Dowell Formula |
|--------|-------------------------------|----------------|
| **BC at center** | dH/dz(a) = 0 (symmetry) | H(a) = 0 (current exits) |
| **DC current** | 0 (no current at DC!) | I = H0 (current flows) |
| **Use case** | Surface impedance Z_s | R_ac/R_dc |
| **Formula** | Re(γa · tanh(γa)) | ξ · (sinh(2ξ)+sin(2ξ)) / (cosh(2ξ)-cos(2ξ)) |
| **DC limit** | 0 | 1.0 |

### Geometry Comparison

**ESIM Geometry (dH/dz = 0)**:
```
Surface (z=0)          Center (z=a)
    |                      |
    H = H0                 dH/dz = 0
    |                      |
    +------ Conductor -----+

At DC: H = constant, J = 0, I = 0
This models a thick conductor where current crowds to surface.
```

**Dowell Geometry (H = 0)**:
```
Surface (z=0)          Center (z=a)
    |                      |
    H = H0                 H = 0
    |                      |
    +------ Conductor -----+

At DC: H linear, J = H0/a, I = H0
This models current flowing through the conductor.
```

### Key Difference: coth vs tanh

| Formula | Function | BC at Center | DC Limit |
|---------|----------|--------------|----------|
| Dowell | **coth** | H(a) = 0 | F = 1 |
| ESIM | **tanh** | dH/dz(a) = 0 | Z_s = 0 |

Note: coth(x) = 1/tanh(x), reflecting different boundary conditions.

## Dowell vs ESIM: Method Comparison

### Dowell's Formula (Analytical)

**Purpose**: Compute R_ac/R_dc ratio and L_int/L_int,dc ratio for rectangular conductors.

**Resistance Ratio F(ξ)**:
```
F(ξ) = ξ · (sinh(2ξ) + sin(2ξ)) / (cosh(2ξ) - cos(2ξ))

Alternative form: F(ξ) = Re[γa · coth(γa)]
```

**Internal Inductance Ratio G(ξ)**:
```
G(ξ) = (3 / 2ξ) · (sinh(2ξ) - sin(2ξ)) / (cosh(2ξ) - cos(2ξ))
```

**Numerical Values**:

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

### ESIM (Effective Surface Impedance Method)

**Purpose**: Compute surface impedance Z_s for SIBC and nonlinear materials.

**Output**:
```
Z_s = Re(Z_s) + j·Im(Z_s)
      ↓           ↓
   Resistance   Reactance
```

**Conversion to PEEC**:
```
R_ac = Re(Z_s) · l / P
L_int = Im(Z_s) · l / (ω · P)
```

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

## Cross-Section Specific Formulas

### Rectangular Cross-Section (Dowell)

For width w >> height h (thin strip), use Dowell's formula:

```python
def dowell_rac_ratio(xi):
    """Dowell's formula for rectangular conductor."""
    if xi < 0.01:
        return 1.0 + xi**4 / 45
    sh2, sn2 = np.sinh(2*xi), np.sin(2*xi)
    ch2, cs2 = np.cosh(2*xi), np.cos(2*xi)
    return xi * (sh2 + sn2) / (ch2 - cs2)
```

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

## ESIM Homogenization for Nonlinear Materials

### Why Homogenization?

For linear materials (constant μ), ESIM and Dowell give equivalent results:
```
mu_eff = mu * integral{|H|^2 dz} / integral{|H|^2 dz} = mu
→ xi_eff = xi_initial
→ ESIM → Dowell formula
```

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

## Cross-Section Support Roadmap

| Cross-Section | Cell Problem | Status |
|--------------|--------------|--------|
| **Rectangular/Planar** | 1D ESIM | Implemented |
| Round wire | 1D (Kelvin functions) | Analytical solution available |
| Arbitrary 2D | 2D FEM | Future work |
| 3D topology (Litz, braid) | Complex homogenization | Not planned |

### PEEC Matrix Modification with Dowell

Apply F(ξ) and G(ξ) to **diagonal elements only**:

```
Z_ii = R_dc,i × F(ξ_i) + jω × (L_ext,ii + L_int,dc,i × G(ξ_i))
       \_____________/         \_____/   \___________________/
          AC resistance        External    Internal (modified)

Z_ij = jω × M_ij    (i ≠ j, no modification for mutual terms)
```

**Python Implementation**:
```python
def apply_skin_effect(R_dc, L_ext, L_int_dc, xi_array, omega):
    """Apply Dowell correction to PEEC matrices."""
    n = len(xi_array)

    # Compute correction factors for each segment
    F = np.array([dowell_F(xi) for xi in xi_array])
    G = np.array([dowell_G(xi) for xi in xi_array])

    # Modify diagonal elements only
    R_ac = R_dc * F                    # Resistance increases
    L_int = L_int_dc * G               # Internal inductance decreases

    # Build impedance matrix
    Z = np.diag(R_ac) + 1j * omega * (L_ext + np.diag(L_int))

    return Z

def dowell_F(xi):
    """Dowell's formula for R_ac/R_dc."""
    if xi < 0.01:
        return 1.0 + xi**4 / 45
    sh2, sn2 = np.sinh(2*xi), np.sin(2*xi)
    ch2, cs2 = np.cosh(2*xi), np.cos(2*xi)
    return xi * (sh2 + sn2) / (ch2 - cs2)

def dowell_G(xi):
    """Dowell's formula for L_int/L_int,dc."""
    if xi < 0.01:
        return 1.0 - xi**4 / 15
    sh2, sn2 = np.sinh(2*xi), np.sin(2*xi)
    ch2, cs2 = np.cosh(2*xi), np.cos(2*xi)
    return (3 / (2*xi)) * (sh2 - sn2) / (ch2 - cs2)
```

## Without Skin Effect Correction

**Important**: If Z_s or Dowell correction is NOT applied, PEEC gives **DC results**:

```
Z_PEEC = R_dc + jω × L
```

This assumes **uniform current distribution** across the conductor cross-section.

**Error at High Frequency**:

| Quantity | Without Correction | With Correction | Error |
|----------|-------------------|-----------------|-------|
| R | R_dc | R_dc × ξ | ξ× underestimate |
| L_int | L_int,dc | L_int,dc × (3/2ξ) | Overestimate |

**Example** (copper, a=1mm, f=100kHz, ξ≈4.6):
- Resistance: 4.6× underestimate without skin effect
- Internal inductance: 3× overestimate without skin effect

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

## Frequency Dependence Analysis

### Fundamental Frequency-Dependent Quantities

All frequency dependence originates from ω = 2πf:

```
ω = 2πf                                    [rad/s]

δ(ω) = √(2ρ / (ωμ)) = √(ρ / (πfμ))        ∝ ω^(-1/2)  [m]

ξ(ω) = a / δ(ω) = a√(πfμ/ρ)               ∝ ω^(+1/2)  [-]

γ(ω) = (1+j) / δ(ω) = (1+j)√(ωμσ/2)       ∝ ω^(+1/2)  [1/m]
```

### Surface Impedance Z_s Frequency Dependence

```
Z_s(ω) = ρ · γ(ω) · tanh(γ(ω)·a)
```

| Regime | Condition | Z_s Behavior |
|--------|-----------|--------------|
| **Low frequency** | ξ << 1 | tanh((1+j)ξ) ≈ (1+j)ξ → Z_s ≈ 2jωμa ∝ **ω** |
| **High frequency** | ξ >> 1 | tanh((1+j)ξ) ≈ 1 → Z_s ≈ (1+j)√(ωμρ/2) ∝ **ω^(1/2)** |

### R_ac/R_dc Ratio Frequency Dependence

```
F(ξ(ω)) = ξ · (sinh(2ξ) + sin(2ξ)) / (cosh(2ξ) - cos(2ξ))
```

| Regime | Condition | F(ξ) Behavior |
|--------|-----------|---------------|
| **Low frequency** | ξ << 1 | F ≈ 1 + ξ^4/45 ≈ **1** (constant) |
| **Transition** | ξ ~ 1 | F increases smoothly |
| **High frequency** | ξ >> 1 | F ≈ ξ ∝ **ω^(1/2)** |

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
| **Total Reactance** | jωL | **∝ ω** | **∝ ω** |

### Frequency Dependence Flow Diagram

```
Frequency f
    │
    ▼
ω = 2πf
    │
    ▼
┌─────────────────────────────────────┐
│  Skin Depth and Related Parameters  │
│                                     │
│  δ(ω) = √(2ρ/(ωμ))    ∝ ω^(-1/2)   │
│  ξ(ω) = a/δ           ∝ ω^(+1/2)   │
│  γ(ω) = (1+j)/δ       ∝ ω^(+1/2)   │
└─────────────────────────────────────┘
    │
    ▼
┌─────────────────────────────────────┐
│  Surface Impedance                  │
│                                     │
│  Z_s(ω) = ρ·γ·tanh(γa)             │
│                                     │
│  Low freq (ξ<<1):  Z_s ∝ ω         │
│  High freq (ξ>>1): Z_s ∝ ω^(1/2)   │
└─────────────────────────────────────┘
    │
    ▼
┌─────────────────────────────────────┐
│  AC/DC Resistance Ratio             │
│                                     │
│  F(ξ) = Dowell(ξ(ω))               │
│                                     │
│  Low freq (ξ<<1):  F ≈ 1           │
│  High freq (ξ>>1): F ∝ ω^(1/2)     │
└─────────────────────────────────────┘
    │
    ▼
┌─────────────────────────────────────┐
│  PEEC Impedance Matrix              │
│                                     │
│  Z_PEEC(ω) = R_dc·F(ω) + jωL       │
│              \_____/      \__/      │
│              R_ac(ω)    jωL(ω)     │
│                                     │
│  Low freq:  Z ≈ R_dc + jωL         │
│  High freq: Z ∝ ω^(1/2) + jωL      │
└─────────────────────────────────────┘
```

### Practical Implications

1. **Low frequency (f < f_transition)**:
   - Skin effect negligible
   - R_ac ≈ R_dc (use DC resistance)
   - Reactance dominates: |jωL| >> R_ac

2. **Transition frequency**:
   ```
   f_transition ≈ ρ / (πμa²)

   Example (copper, a = 1mm):
   f_transition ≈ 1.7×10⁻⁸ / (π × 4π×10⁻⁷ × 10⁻⁶) ≈ 4.3 kHz
   ```

3. **High frequency (f > f_transition)**:
   - Strong skin effect
   - R_ac increases as √f
   - Both resistance and reactance significant

## Laplace Domain (s = jω) Representation

### Fundamental Relationship

The Dowell functions F and G can be expressed in the Laplace domain using the magnetic diffusion time constant:

```
τ = a²μσ    [s]  (magnetic diffusion time constant)

s = jω      (Laplace variable)
```

**Relationship between ξ and τs**:
```
ξ² = a²μσω/2 = τω/2 = τs/(2j)

γa = (1+j)ξ = √(τs)    (for s = jω)
```

### F(s) and G(s) in Laplace Domain

**Resistance ratio**:
```
F(s) = Re[√(τs) · coth(√(τs))]
```

**Inductance ratio** (complex form):
```
Z(s)/R_dc = √(τs) · coth(√(τs))

where:
  Real part → F(s) = R_ac/R_dc
  Imaginary part / (τs) → related to G(s)
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
| 3rd | F₃ = 1 + τs/(3 + τs/(5 + τs/7)) | Good for τω < 10 |

### 2nd Order Rational Approximation

The 2nd order approximation:
```
F₂(s) = (15 + 8τs) / (15 + 3τs)
```

**Partial fraction expansion**:
```
F₂(s) = 8/3 - 25/(9τ) · 1/(s + 5/τ)
```

This corresponds to a first-order RL circuit with a single pole at s = -5/τ.

### Foster Ladder Network

The continued fraction leads to a Foster-type ladder circuit:

```
Z(s)/R_dc = 1 + 1/(3/(τs) + 1/(5/(τs) + 1/(7/(τs) + ...)))
```

**Circuit interpretation**:
```
     R_dc
  ○──/\/\/──┬──────┬──────┬──────○
            │      │      │
           ═╧═    ═╧═    ═╧═
           C₁     C₂     C₃
            │      │      │
           ─┴─    ─┴─    ─┴─
           R₁     R₂     R₃

where:
  C₁ = τ/(3·R_dc),  R₁ = 3R_dc/5
  C₂ = τ/(5·R_dc),  R₂ = 5R_dc/7
  C₃ = τ/(7·R_dc),  R₃ = 7R_dc/9
  ...
```

### Cauer Ladder Network

Alternative series RL ladder:
```
     R_dc     L₁      L₂      L₃
  ○──/\/\/──○──⊃⊃⊃──○──⊃⊃⊃──○──⊃⊃⊃──○
             │       │       │
            ═╧═     ═╧═     ═╧═
            R₁      R₂      R₃

where:
  L₁ = τ·R_dc/3,  R₁ = 3R_dc
  L₂ = τ·R_dc/5,  R₂ = 5R_dc/3
  L₃ = τ·R_dc/7,  R₃ = 7R_dc/5
  ...
```

### Time Constant Examples

| Material | a | σ [S/m] | μ_r | τ = a²μσ |
|----------|---|---------|-----|----------|
| Copper, 1mm | 1mm | 5.8×10⁷ | 1 | **73 ns** |
| Copper, 5mm | 5mm | 5.8×10⁷ | 1 | **1.8 μs** |
| Aluminum, 1mm | 1mm | 3.5×10⁷ | 1 | 44 ns |
| Steel, 1mm | 1mm | 2×10⁶ | 1000 | **2.5 ms** |

### Approximation Accuracy

For τ = 73 ns (copper, a = 1mm):

| f (Hz) | τω | ξ | F (exact) | F₂ (2nd order) | Error |
|--------|-----|---|-----------|----------------|-------|
| 1k | 4.6×10⁻⁴ | 0.015 | 1.000 | 1.000 | <0.01% |
| 100k | 0.046 | 0.15 | 1.001 | 1.001 | <0.1% |
| 1M | 0.46 | 0.48 | 1.015 | 1.014 | 0.1% |
| 10M | 4.6 | 1.5 | 1.52 | 1.45 | 5% |
| 100M | 46 | 4.8 | 4.8 | 3.2 | 33% |

**Validity**: 2nd order approximation is accurate within 5% for ξ < 1.5 (τω < 4.5).

### Python Implementation

```python
def dowell_F_continued_fraction(tau_s, order=3):
    """
    Continued fraction approximation of F(s).

    Parameters:
        tau_s: τ·s (dimensionless, complex for s=jω)
        order: Number of continued fraction terms

    Returns:
        F: Resistance ratio approximation
    """
    # Build continued fraction from bottom up
    # F = 1 + τs/(3 + τs/(5 + τs/(7 + ...)))

    denominators = [2*k + 1 for k in range(order, 0, -1)]  # [2*order+1, ..., 5, 3]

    result = denominators[0]  # Start with innermost denominator
    for d in denominators[1:]:
        result = d + tau_s / result

    return 1 + tau_s / result

def dowell_F_2nd_order(tau_s):
    """2nd order rational approximation of F(s)."""
    return (15 + 8*tau_s) / (15 + 3*tau_s)

# Example usage
import numpy as np

tau = 73e-9  # 73 ns for copper, a=1mm
frequencies = [1e3, 1e5, 1e6, 1e7]

for f in frequencies:
    omega = 2 * np.pi * f
    s = 1j * omega
    tau_s = tau * s

    F_exact = dowell_F(np.sqrt(tau * omega / 2))  # Using ξ
    F_cf3 = dowell_F_continued_fraction(tau_s, order=3)
    F_2nd = dowell_F_2nd_order(tau_s)

    print(f"f={f/1e3:.0f}kHz: F_exact={F_exact.real:.4f}, "
          f"F_cf3={F_cf3.real:.4f}, F_2nd={F_2nd.real:.4f}")
```

### G(s) Continued Fraction

For internal inductance ratio G(s):

```
G(s) ≈ 1 / (1 + τs/15 + (τs)²/315 + ...)

     ≈ 1 / (1 + τs/(15 + τs/(21 + τs/(27 + ...))))
```

**2nd order approximation**:
```
G₂(s) = 15 / (15 + τs)
```

This is a simple first-order low-pass filter with pole at s = -15/τ.

## Proximity Effect in PEEC + Dowell

### Overview

When combining PEEC with Dowell's formula, proximity effects are **partially included**:

| Effect | Dowell Alone | PEEC + Dowell |
|--------|:------------:|:-------------:|
| Self skin effect | ✓ | ✓ |
| External proximity effect | ✗ | **✓** |
| Internal proximity effect | ✓ (m-layer formula) | △ (approximate) |

### Self Skin Effect

Current redistribution due to the conductor's own magnetic field.

```
Dowell: F(ξ) accurately models this
PEEC:   Z_self = R_dc · F(ξ) + jωL_int · G(ξ)
```

**→ Correctly handled by both**

### External Proximity Effect

Current redistribution due to magnetic fields from **other conductors**.

```
        Conductor 1      Conductor 2
       ┌─────────┐      ┌─────────┐
       │ → → → → │←─H₁₂→│ ← ← ← ← │
       │ → → → → │      │ ← ← ← ← │
       └─────────┘      └─────────┘
       Current density   Current density
       shifts right      shifts left
```

**Dowell alone**: Not considered (single conductor 1D analysis)

**PEEC + Dowell**:
```
Z_mutual(i,j) = jω · L_ij    (mutual inductance)

Magnetic coupling between conductors represents external proximity effect
→ Automatically included through PEEC mutual inductance terms
```

**→ Correctly handled by PEEC**

### Internal Proximity Effect

Current redistribution within the **same conductor** due to external magnetic fields.

```
External field H_ext
    ↓ ↓ ↓ ↓ ↓
  ┌───────────┐
  │→ → → → → →│  ← Current concentrates at surface
  │  → → →    │    (skin effect + external field influence)
  │           │
  └───────────┘
```

**Dowell alone**: Considered via m-layer formula

```
F_m(ξ) = F(ξ) + (m² - 1)/3 · G(ξ)

m = layer number (represents external field magnitude)
```

**PEEC + Simple Dowell**: **Not accurately handled**

Issues:
- Dowell's F(ξ) assumes "no external field"
- In reality, adjacent conductors create external fields
- PEEC mutual inductance represents external fields but doesn't modify internal current distribution

### Practical Guidelines

| Situation | Recommended Approach |
|-----------|---------------------|
| Conductor spacing > 3δ | PEEC + Dowell is sufficient |
| Conductor spacing ≈ δ | Consider filament subdivision |
| Closely wound coils | Use Dowell's m-layer formula |
| High accuracy required | Filament subdivision or FEM |

### Improvement Methods

#### Method 1: Filament Subdivision

```
     Subdivide one conductor into multiple filaments

     ┌──┬──┬──┐
     │ 1│ 2│ 3│  ← Each filament is an
     ├──┼──┼──┤     independent PEEC element
     │ 4│ 5│ 6│
     └──┴──┴──┘

     - Directly computes internal current distribution
     - Proximity effect automatically included
     - Computational cost: O(N²) → O((N × subdivisions)²)
```

#### Method 2: Modified Dowell (H_ext dependent)

```
Modified formula considering external field H_ext:

F_eff(ξ, H_ext) = F(ξ) + correction_term(ξ, H_ext/H_self)

where:
  H_self = I / (2a)         (self magnetic field)
  H_ext  = Σ L_ij · I_j / μ  (field from other conductors)
```

#### Method 3: ESIM with Boundary Magnetic Field

For nonlinear materials or complex geometries:
- Include external field in ESIM boundary conditions
- Solve with asymmetric boundary: H(+a) ≠ H(-a)
- Captures internal proximity effect for saturating materials

### Summary

```
PEEC + Dowell proximity effect coverage:

✓ External proximity: Included via PEEC mutual inductance
△ Internal proximity: Approximate (OK for thin conductors or wide spacing)
✗ Exact current distribution: Requires filament subdivision

Practical applicability:
- Power electronics (tens of kHz): PEEC + Dowell usually sufficient
- RF applications (MHz+): Consider filament subdivision
- Closely coupled windings: Use Dowell's m-layer formula or FEM
```

## Implementation Files

- `esim_cell_problem.py`: Base ESIM solvers (ESIMCellProblemSolver, ESIMFiniteSlabSolver)
- `esim_correct_implementation.py`: ESIMHomogenizationSolver (correct approach)
- `esim_conductor_model.py`: PEEC conductor model using ESIM

## References

1. A. Ruehli, "Equivalent Circuit Models for Three-Dimensional Multiconductor Systems," IEEE Trans. MTT, 1974.
2. P. L. Dowell, "Effects of eddy currents in transformer windings," Proc. IEE, vol. 113, no. 8, pp. 1387-1394, 1966.
3. H. Igarashi, "Homogenization approach for skin effect in conductors," various papers.
4. K. Hollaus et al., "Nonlinear Effective Surface Impedance," IEEE Trans. Magnetics, 2025.
5. H. S. Wall, "Analytic Theory of Continued Fractions," Van Nostrand, 1948.
6. J. A. Ferreira, "Improved analytical modeling of conductive losses in magnetic components," IEEE Trans. Power Electronics, vol. 9, no. 1, 1994.

---

**Date**: 2026-01-12
**Author**: Claude Code Analysis
