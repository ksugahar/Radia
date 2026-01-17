# N-Port Block Lanczos SPICE Generation

## Overview

This document describes the N-port Block Lanczos algorithm for SPICE-compatible circuit extraction from coupled electromagnetic systems (conductors + magnetic materials + dielectrics).

## Scalar vs Block Lanczos: When to Use Which

| Ports | Lanczos Type | Rationale |
|-------|--------------|-----------|
| **1-port** | Scalar Lanczos OK | Single starting vector captures port behavior |
| **N-port (N≥2)** | **Block Lanczos required** | Need p starting vectors for p ports |

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
LoopStar + MMM (magnetic) coupled system.

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

## References

1. A. Odabasioglu, M. Celik, L.T. Pileggi, "PRIMA: Passive Reduced-order Interconnect Macromodeling Algorithm," IEEE TCAD, 1998.

2. P.J. Dowell, "Effects of eddy currents in transformer windings," Proc. IEE, 1966.

3. van Oosterom & Strackee, "The Solid Angle of a Plane Triangle," IEEE Trans. BME, 1983.

---
Last Updated: 2026-01-17
