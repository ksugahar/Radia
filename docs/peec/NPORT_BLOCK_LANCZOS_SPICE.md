# N-Port Block Lanczos SPICE Generation

## Overview

This document describes the N-port Block Lanczos algorithm for SPICE-compatible circuit extraction from coupled electromagnetic systems (conductors + magnetic materials + dielectrics).

## Scalar vs Block Lanczos: When to Use Which

| Ports | Lanczos Type | Rationale |
|-------|--------------|-----------|
| **1-port** | Scalar Lanczos OK | Single starting vector captures port behavior |
| **N-port (N>=2)** | **Block Lanczos required** | Need p starting vectors for p ports |

**WARNING**: Using scalar Lanczos for N-port systems is **meaningless** - it cannot capture multi-port interactions correctly.

## Key Insight: Block Lanczos with DC Starting Vectors

For N-port systems, **Block Lanczos** (not scalar Lanczos) must be used for magnetic and dielectric subsystems. The critical insight is:

**Starting vectors = DC solution from each port**

```
At DC: Z_L = R_L, so I_L = R_L^{-1} @ V_port

For magnetic: B_M = K_LM.T @ (R_L^{-1} @ B_L)   [n_M x p]
For dielectric: B_D = K_LD.T @ (R_L^{-1} @ B_L) [n_D x p]
```

This captures how port excitation propagates through coupling to each subsystem.

## Convergence Comparison

### Scalar Lanczos (old approach)
```
k_D=2:  49.20% error
k_D=4:  49.20% error
k_D=6:  15.22% error
k_D=8:   0.00% error (no compression)
```

### Block Lanczos with DC starting vectors (new approach)
```
k_D=1:   0.01% error
k_D=2:   0.00% error  <- Full accuracy with only 2 stages!
```

Block Lanczos converges **dramatically faster** because it captures the port-relevant modes directly.

## System Structure

### LoopStar + Magnetic (L + M)
```
[sL_L + R_L    K_LM  ] [I_L]   [V_port]
[K_ML      sL_M + R_M] [M  ] = [0     ]
```

### LoopStar + Dielectric (L + D)
```
[sL_L + R_L    K_LD  ] [I_L]   [V_port]
[K_DL      sC_D + G_D] [V_D] = [0     ]
```

### Full Coupled (L + M + D)
```
[sL_L + R_L    K_LM       K_LD     ] [I_L]   [V_port]
[K_ML      sL_M + R_M     0        ] [M  ] = [0     ]
[K_DL         0       sC_D + G_D   ] [V_D]   [0     ]
```

## Port Impedance via Schur Complement

After Block Lanczos reduction, port DOFs are at the first p indices.

```python
# Schur complement eliminates non-port DOFs
Z_eff = Z_L - K_LM @ inv(Y_M) @ K_LM.T - K_LD @ inv(Y_D) @ K_LD.T

# Port impedance: extract first p x p submatrix
Z_port = Z_eff[:p, :p]
```

## Complex Material Support

### Complex Permittivity
```
eps(f) = eps_0 * (eps'_r - j*eps''_r)

Circuit equivalent:
  C_eff = C * eps'_r
  G_loss = omega * C * eps''_r  (parallel conductance)
```

### Complex Permeability
```
mu(f) = mu_0 * (mu'_r - j*mu''_r)

Circuit equivalent:
  L_eff = L * mu'_r
  R_loss = omega * L * mu''_r  (series resistance)
```

### Loss Tangent Conversion
```
tan(delta_e) = eps''/eps'  -> eps'' = eps' * tan(delta_e)
tan(delta_m) = mu''/mu'    -> mu'' = mu' * tan(delta_m)
```

## API Classes

### NPortBlockLanczosSPICE
Basic conductor-only N-port reduction.

```python
solver = NPortBlockLanczosSPICE(n_ports=2, n_stages=5)
result = solver.reduce(L, R, port_indices)
netlist = solver.to_spice(result, "COIL_2PORT")
```

### NPortCoupledMagneticSPICE
LoopStar + magnetic-material reduced model.

```python
solver = NPortCoupledMagneticSPICE(n_ports=2, n_stages_L=5, n_stages_M=3)
result = solver.reduce(L_L, R_L, L_M, R_M, K_LM, port_indices)
netlist = solver.to_spice(result, "COIL_WITH_CORE")
```

### NPortCoupledDielectricSPICE
LoopStar + Dielectric coupled system.

```python
solver = NPortCoupledDielectricSPICE(n_ports=2, n_stages_L=4, n_stages_D=2)
result = solver.reduce(L_L, R_L, C_D, G_D, K_LD, port_indices)
netlist = solver.to_spice(result, "LC_FILTER")
```

### NPortFullCoupledSPICE
Full L + M + D coupled system with complex materials.

```python
ferrite = ComplexMaterialParams(mu_r_real=2000, tan_delta_m=0.01)
pcb = ComplexMaterialParams(eps_r_real=4.5, tan_delta_e=0.02)

solver = NPortFullCoupledSPICE(n_ports=2, n_stages_L=4, n_stages_M=2, n_stages_D=2)
result = solver.reduce(L_L, R_L, L_M, R_M, C_D, G_D, K_LM, K_LD, port_indices,
                       material_M=ferrite, material_D=pcb)
netlist = solver.to_spice(result, "TRANSFORMER_ON_PCB", f_ref=100e3)
```

## SPICE Output Format

Generated netlists use standard SPICE subcircuit format:

```spice
* N-Port Full Coupled (L + M + D) SPICE Netlist
* Conductor stages: 4, Magnetic stages: 2, Dielectric stages: 2
* Reference frequency for loss: 100.0 kHz
* Magnetic: mu'_r = 2000, tan(delta_m) = 0.01
* Dielectric: eps'_r = 4.5, tan(delta_e) = 0.02

.SUBCKT TRANSFORMER_ON_PCB p1_in p1_out p2_in p2_out

* === Conductor Inductance Ladder ===
L1 p1_in cond_s1_p1 1.000000e-05
...

* === Magnetic Core Ladder ===
Lmag1 mag_in mag_n1 1.000000e-01
Rmag_loss1 mag_in mag_loss1 1.256637e+01
...

* === Dielectric Capacitor Ladder ===
Cdiel1 cap_in cap_n1 2.250000e-10
Gdiel_loss1 cap_in 0 2.827433e-05
...

.ENDS TRANSFORMER_ON_PCB
```

## Frequency-Dependent Materials

For materials with frequency-dependent properties (Debye relaxation, Cole-Cole, etc.), use:

1. **Verilog-A**: Direct frequency-dependent impedance modeling
2. **Continued Fraction Expansion (CFE)**: Rational approximation -> RL/RC ladder

The current implementation uses a reference frequency `f_ref` for loss elements. For wideband accuracy, CFE or Verilog-A is recommended.

## Test Results

### LoopStar + Dielectric (demo_coupled_dielectric)
```
System: 15 conductor DOF, 8 dielectric DOF, 2 ports
Reduction: 15 -> 8 conductor, 8 -> 4 dielectric (2 Block Lanczos stages x 2 ports)
Compression: 56.5%
Z_11 error: 0.00%
Z_12 error: 0.03%
```

### Full L + M + D (demo_full_coupled_lmd)
```
System: 12 conductor, 6 magnetic, 4 dielectric DOF, 2 ports
Reduction: 22 -> 12 total DOF (45.5% compression)
Complex materials: ferrite (mu_r=2000, tan_delta_m=1%), PCB (eps_r=4.5, tan_delta_e=2%)
```

## Design Decision: PRIMA Only, No ACA

**Policy (2026-01-17)**: For PEEC interaction matrices, use **PRIMA Lanczos only**. ACA (Adaptive Cross Approximation) is **not used**.

**Rationale**:
1. **PRIMA directly reduces element count**: Lanczos produces a tridiagonal matrix that maps directly to an RL/RC ladder in SPICE
2. **ACA only compresses storage**: Low-rank approximation doesn't reduce circuit complexity
3. **Port impedance preservation**: PRIMA with DC starting vectors preserves port behavior exactly
4. **Simplicity**: Single reduction method, no need for ACA tolerance tuning

**Comparison**:
| Method | Storage | SPICE Elements | Port Accuracy |
|--------|---------|----------------|---------------|
| Full | O(N^2) | N elements | Exact |
| ACA | O(N log N) | Still N elements | Approximate |
| **PRIMA** | O(k) | **k elements** | Exact at ports |

## PyKAN for Continued Fraction: Evaluation (2026-01-17)

**Question**: Can PyKAN learn continued fraction expansions and extract RL ladder parameters?

**Test Results**:
| Test | Result | Notes |
|------|--------|-------|
| z*coth(z) function learning | PASS | 1.34% max error |
| Symbolic regression to CF | FAIL | Cannot extract rational form |
| Ladder coefficient extraction | FAIL | Derivative-based extraction unreliable |

**Conclusion**: KAN can **approximate** skin effect functions accurately, but **cannot extract** the continued fraction structure [3, 5, 7, ...].

**Recommended Approach**:
1. **Standard skin effect**: Use analytical Dowell formula (exact CF coefficients)
2. **Nonlinear Z(H,f) steady-state**: Use KAN to learn from ESIM data, export as **PWL lookup table**
3. **Nonlinear Z(H,f) transient**: **DEFERRED** - requires Verilog-A with dynamic H tracking
4. **DO NOT** attempt to extract ladder structure from KAN

```python
# RECOMMENDED: Analytical Dowell for standard skin effect
dowell_coeffs = [3, 5, 7, 9, 11]  # Exact continued fraction

# FOR NONLINEAR: KAN with PWL export (not CF extraction)
kan_model.train_from_esim_data(H_values, f_values, Zs_data)
spice_pwl = kan_model.to_spice_pwl(H_export, f_ref)
```

## CLN Type-I Coordinate Transform Verification

This section documents the verification of the CLN I (Cauer Ladder Network Type I) coordinate transform used in PEEC model order reduction. The CLN I transform is the s=0 (DC) expansion variant of the Lanczos algorithm that converts dense PEEC matrices into sparse tridiagonal/diagonal form.

### CLN I Transform Definition

Given two Hermitian matrices R (resistance, diagonal) and L (inductance, dense), the CLN I transform finds a transformation matrix Q such that:

```
R_diag    = Q^T * R * Q    (diagonal matrix)
L_tridiag = Q^T * L * Q    (tridiagonal matrix)
```

```python
from cln import lanczos

result = lanczos(K=R, N=L)
# result.R_diag:    diagonal matrix (resistance)
# result.L_tridiag: tridiagonal matrix (inductance)
# result.U, result.V: transformation matrices (U = V = Q for symmetric case)
```

### Key Property: U = V = Q

For symmetric PEEC matrices, the Lanczos algorithm produces identical left and right transformation matrices:

```
||U - V||_F ~ 1e-16  (machine precision)
```

This means a single transformation matrix Q fully characterizes the coordinate change.

### K-Inner Product Orthogonality

Q is orthogonal under the K-inner product but **not** under the standard inner product:

```
Q^T * K * Q = R_diag     (diagonal = K-inner product orthogonal)
Q^T * Q != I              (NOT standard-orthogonal)
||Q^T * Q - I||_F = 1.105e+09  (very large)
```

This distinction is critical when adding terms like surface impedance Z_s to the transformed system.

### Verification Test 1: Without Surface Impedance

**Original PEEC circuit:**
```
Z(s) = R + s*L  (dense matrices)
```

**After CLN I transform:**
```
Z'(s) = R_diag + s*L_tridiag  (sparse matrices)
```

**Result:**
- Max relative error: 8.23e-16
- Mean relative error: 2.09e-16
- **PASS** (machine precision)

### Verification Test 2: With Surface Impedance Z_s

When surface impedance (skin effect) is included:

**Original:**
```
Z(s) = R + s*L + Z_s(s)*I  (dense matrices)
```

**After CLN I transform:**
```
Z'(s) = R_diag + s*L_tridiag + Z_s(s)*(Q^T * Q)  (sparse + transformed identity)
```

Since Q^T * Q != I, the transformed identity matrix `Q^T * Q` must be used (not I).

**Surface impedance Z_s:**
```
Z_s = (1+j) * sqrt(pi * f * mu / sigma)
```

Copper conductor (sigma = 5.8e7 S/m) examples:

| Frequency | Z_s | Skin depth delta |
|-----------|-----|------------------|
| 10 Hz | 8.25e-7 + j8.25e-7 | 20.9 mm |
| 10 kHz | 2.70e-5 + j2.70e-5 | 0.64 mm |
| 10 MHz | 8.25e-4 + j8.25e-4 | 0.021 mm |

**Result:**
- Max relative error: 9.53e-16
- Mean relative error: 2.10e-16
- **PASS** (machine precision)

### Test Conditions

```python
n = 10  # number of loops

# Dense inductance matrix (with mutual inductance)
L0 = 1e-6  # 1 uH
decay = 3.0
L_dense[i,j] = L0 * exp(-|i-j| / decay)

# Diagonal resistance matrix
R0 = 0.01  # 10 mOhm
R_diag = diag([R0, R0, ..., R0])

# Frequency range
frequencies = logspace(1, 7, 100)  # 10 Hz to 10 MHz
```

### Coupling Matrix Transformation Rules for Retired PEEC-Magnetic/STAR Systems

When coupling the CLN I-transformed LoopStar subsystem with other physics
(magnetic-material or STAR subsystems), the coupling matrices require one-sided
transformation:

| Coupling Term | Original | CLN I Coordinates | Note |
|---------------|----------|-------------------|------|
| Z_s * I | Z_s * I | Z_s * (Q^T * Q) | Surface impedance |
| Z_LM | Z_LM | Q^T * Z_LM | Loop to magnetic subsystem (left multiply) |
| Z_ML | Z_ML | Z_ML * Q | Magnetic subsystem to Loop (right multiply) |
| Z_LS | Z_LS | Q^T * Z_LS | Loop to STAR |
| Z_SL | Z_SL | Z_SL * Q | STAR to Loop |

**Important rules:**
1. **Magnetic/STAR subsystem matrices are NOT transformed** - Z_MAG, Z_STAR remain unchanged
2. **Coupling matrices are one-sided** - Z_LM gets Q^T on the left; Z_ML gets Q on the right
3. **Excitation vectors are also transformed** - Loop excitation V_L becomes Q^T * V_L

### Coupled System in CLN I Coordinates

Original coupled system:
```
[R + sL    Z_LM  ] [I_L]   [V_L]
[Z_ML     Z_MAG ] [I_M] = [V_M]
```

After CLN I transform:
```
[R_diag + s*L_tridiag    Q^T*Z_LM ] [I_L']   [Q^T*V_L]
[Z_ML*Q                  Z_MAG    ] [I_M ] = [V_M    ]
```

where `I_L' = Q^{-1} * I_L` is the transformed Loop current.

### Verification Artifacts

- `cln_frequency_response_verification.png` - Test 1 results
- `cln_frequency_response_with_zs.png` - Test 2 results (with Z_s)
- `verify_cln_frequency_response.py` - Frequency response verification script
- `verify_lanczos_uv_equality.py` - U, V matrix and CLN I structure verification script

## References

1. A. Odabasioglu, M. Celik, L.T. Pileggi, "PRIMA: Passive Reduced-order Interconnect Macromodeling Algorithm," IEEE TCAD, 1998.

2. P.J. Dowell, "Effects of eddy currents in transformer windings," Proc. IEE, 1966.

3. van Oosterom & Strackee, "The Solid Angle of a Plane Triangle," IEEE Trans. BME, 1983.

4. Lanczos algorithm for generalized eigenvalue problems and Cauer circuit model order reduction of RLC equivalent circuits.

5. PEEC (Partial Element Equivalent Circuit) method.

---
Last Updated: 2026-02-22
