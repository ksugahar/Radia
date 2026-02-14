# PEEC-MSC Validation Plan

**Status**: Testing Strategy (2026-02-12)
**Purpose**: Systematic validation of PEEC-MSC coupled solver

---

## Table of Contents

1. [Validation Philosophy](#validation-philosophy)
2. [Phase 1: PEEC Single Conductor](#phase-1-peec-single-conductor)
3. [Phase 2: MSC with External Field](#phase-2-msc-with-external-field)
4. [Phase 3: PEEC-MSC Linear Coupling](#phase-3-peec-msc-linear-coupling)
5. [Phase 4: PEEC-MSC Nonlinear Coupling](#phase-4-peec-msc-nonlinear-coupling)
6. [Phase 5: Complex Geometries](#phase-5-complex-geometries)
7. [Acceptance Criteria](#acceptance-criteria)

---

## Validation Philosophy

### Pyramid of Validation

```
         ┌─────────────────────┐
         │  Complex Geometries │  Phase 5: Real applications
         │  (Transformers, WPT)│
         └─────────────────────┘
               ▲
               │
         ┌─────────────────────┐
         │ Nonlinear Coupling  │  Phase 4: μ_r(H), Newton-Raphson
         │  PEEC + MSC(B-H)    │
         └─────────────────────┘
               ▲
               │
         ┌─────────────────────┐
         │  Linear Coupling    │  Phase 3: PEEC + MSC(μ_r=const)
         │  PEEC + MSC         │
         └─────────────────────┘
               ▲
               │
      ┌────────┴────────┐
      │                 │
┌─────────────┐   ┌─────────────┐
│ PEEC Alone  │   │  MSC Alone  │  Phase 1-2: Individual solvers
│ (No MSC)    │   │  (No PEEC)  │
└─────────────┘   └─────────────┘
```

### Validation Principles

1. **Bottom-up**: Validate simple cases before complex
2. **Analytical first**: Compare to closed-form solutions when available
3. **Cross-verification**: Compare to existing tools (CplMag, ELF/MAGIC, FEM)
4. **Convergence study**: Mesh refinement, solver tolerance
5. **Energy conservation**: Verify power balance

---

## Phase 1: PEEC Single Conductor

**Status**: 🔄 Next Priority
**Goal**: Verify PEEC solver without magnetic materials

### Test 1.1: Circular Loop Self-Inductance

**Geometry**:
- Circular wire loop, radius **a = 0.05 m**, wire radius **r_wire = 1 mm**
- Material: Copper (σ = 5.8e7 S/m)
- Frequency: DC, 1 kHz, 10 kHz, 100 kHz

**Analytical solution** (Neumann formula):
```
L_self ≈ μ_0 · a · [ln(8a/r_wire) - 2]  (for r_wire << a)
```

For **a = 0.05 m**, **r_wire = 0.001 m**:
```
L_self ≈ 4π×10^-7 × 0.05 × [ln(8×0.05/0.001) - 2]
       ≈ 4π×10^-7 × 0.05 × [ln(400) - 2]
       ≈ 4π×10^-7 × 0.05 × [5.99 - 2]
       ≈ 2.51e-7 H = 251 nH
```

**PEEC mesh**:
- Discretize loop into **N_segments** (try N = 8, 16, 32, 64)
- Use filamentary approximation or rectangular cross-section

**Validation**:
```python
import radia as rad
import numpy as np

rad.FldUnits('m')

# Create PEEC loop
a = 0.05  # radius
r_wire = 0.001  # wire radius
sigma = 5.8e7  # copper
N_seg = 32

loop = rad.PEECLoop(center=[0,0,0], radius=a, axis=[0,0,1],
                    wire_radius=r_wire, sigma=sigma, n_segments=N_seg)

# Compute self-inductance at DC
freq = 0  # DC
Z = rad.PEECImpedance(loop, freq)
L_peec = np.imag(Z) / (2 * np.pi * freq) if freq > 0 else Z.imag

# Analytical
L_analytical = 4e-7 * np.pi * a * (np.log(8*a/r_wire) - 2)

# Compare
error = abs(L_peec - L_analytical) / L_analytical * 100
print(f"L_peec = {L_peec*1e9:.2f} nH")
print(f"L_analytical = {L_analytical*1e9:.2f} nH")
print(f"Error = {error:.2f}%")
```

**Acceptance**: Error < 5% for N_seg ≥ 32

**Convergence study**:
| N_seg | L (nH) | Error (%) |
|-------|--------|-----------|
| 8     | ?      | ?         |
| 16    | ?      | ?         |
| 32    | ?      | < 5%      |
| 64    | ?      | < 2%      |

---

### Test 1.2: Rectangular Conductor DC Resistance

**Geometry**:
- Rectangular conductor: length **L = 0.1 m**, width **w = 10 mm**, thickness **t = 1 mm**
- Material: Copper (σ = 5.8e7 S/m)
- Frequency: DC

**Analytical solution**:
```
R_dc = ρ · L / A = (1/σ) · L / (w · t)
```

For **L = 0.1 m**, **w = 10 mm**, **t = 1 mm**, **σ = 5.8e7 S/m**:
```
R_dc = (1/5.8e7) · 0.1 / (0.01 × 0.001)
     = 1.724e-8 · 0.1 / 1e-5
     = 1.724e-4 Ω = 0.1724 mΩ
```

**PEEC mesh**:
- Discretize conductor surface into rectangular panels
- Use filamentary model or volume discretization

**Validation**:
```python
# Create rectangular conductor
conductor = rad.PEECRectangle(
    corner=[0, 0, 0],
    length=0.1, width=0.01, thickness=0.001,
    sigma=5.8e7, n_length=10, n_width=5
)

# Compute DC resistance
R_peec = np.real(rad.PEECImpedance(conductor, freq=0))

# Analytical
R_analytical = (1/5.8e7) * 0.1 / (0.01 * 0.001)

# Compare
error = abs(R_peec - R_analytical) / R_analytical * 100
print(f"R_peec = {R_peec*1e3:.4f} mΩ")
print(f"R_analytical = {R_analytical*1e3:.4f} mΩ")
print(f"Error = {error:.2f}%")
```

**Acceptance**: Error < 1%

---

### Test 1.3: Skin Effect (Dowell's Formula)

**Geometry**:
- Rectangular conductor: length **L = 0.1 m**, width **w = 10 mm**, thickness **t = 0.1 mm** (thin foil)
- Material: Copper (σ = 5.8e7 S/m, μ_r = 1)
- Frequency: 1 kHz, 10 kHz, 100 kHz, 1 MHz

**Analytical solution** (Dowell's formula):
```
ξ = t / (2δ)  where δ = sqrt(2 / (ω μ_0 σ))
F_R(ξ) = ξ · [sinh(2ξ) + sin(2ξ)] / [cosh(2ξ) - cos(2ξ)]
R_ac = R_dc · F_R(ξ)
```

**Validation**:
```python
# Compute impedance vs frequency
freqs = [1e3, 1e4, 1e5, 1e6]
R_peec = []
R_dowell = []

for f in freqs:
    # PEEC with SIBC
    Z = rad.PEECImpedance(conductor, freq=f, sibc=True)
    R_peec.append(np.real(Z))

    # Dowell's formula
    delta = np.sqrt(2 / (2*np.pi*f * 4e-7*np.pi * 5.8e7))
    xi = 0.0001 / (2*delta)
    F_R = xi * (np.sinh(2*xi) + np.sin(2*xi)) / (np.cosh(2*xi) - np.cos(2*xi))
    R_dowell.append(R_dc * F_R)

# Plot comparison
import matplotlib.pyplot as plt
plt.loglog(freqs, R_peec, 'o-', label='PEEC')
plt.loglog(freqs, R_dowell, 's--', label='Dowell')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Resistance (Ω)')
plt.legend()
plt.grid(True)
plt.savefig('validation_skin_effect.pdf')
```

**Acceptance**: Error < 5% at all frequencies

---

## Phase 2: MSC with External Field

**Status**: ✅ Complete (MSC solver validated)
**Goal**: Verify MSC solver with external background field

### Test 2.1: Single Hex in Uniform Field (Linear)

**Geometry**:
- Single hexahedron: 1 cm cube
- Material: Linear (μ_r = 1000)
- External field: **H_ext = [0, 0, 10000]** A/m (uniform)

**Analytical solution**:
```
M = χ · H_ext = (μ_r - 1) · H_ext = 999 · [0, 0, 10000] = [0, 0, 9.99e6] A/m
```

**Validation**:
```python
import radia as rad
rad.FldUnits('m')

# Create hex
vertices = [[-0.005,-0.005,-0.005], [0.005,-0.005,-0.005],
            [0.005,0.005,-0.005], [-0.005,0.005,-0.005],
            [-0.005,-0.005,0.005], [0.005,-0.005,0.005],
            [0.005,0.005,0.005], [-0.005,0.005,0.005]]
hex = rad.ObjHexahedron(vertices, [0, 0, 0])

mat = rad.MatLin(1000)
rad.MatApl(hex, mat)

# Background field
bkg = rad.ObjBckg(lambda p: [0, 0, 10000/796000])  # T (mu_0*H)

container = rad.ObjCnt([hex, bkg])
rad.Solve(container, 0.0001, 100, 0)

# Get magnetization
M = rad.ObjM(hex)
M_z = M[2]

# Analytical
M_analytical = 999 * 10000

# Compare
error = abs(M_z - M_analytical) / M_analytical * 100
print(f"M_z = {M_z:.2e} A/m")
print(f"M_analytical = {M_analytical:.2e} A/m")
print(f"Error = {error:.2f}%")
```

**Acceptance**: Error < 1%

**Status**: ✅ Already validated (MSC linear solver works)

---

### Test 2.2: Two Hex Mutual Interaction

**Goal**: Verify MSC interaction matrix (hex-hex coupling)

**Geometry**:
- Two hexahedra: 1 cm cubes, separated by 2 cm (center-to-center)
- Material: Linear (μ_r = 500)
- External field: **H_ext = [0, 0, 10000]** A/m

**Expected**:
- Cube 1 magnetization → field at Cube 2
- Cube 2 magnetization → field at Cube 1
- Mutual coupling reduces total M compared to isolated case

**Validation**: Compare to analytical dipole-dipole interaction or FEM

**Status**: ✅ Covered by existing benchmarks

---

### Test 2.3: Nonlinear Material (B-H Curve)

**Goal**: Verify Newton-Raphson solver for nonlinear μ_r(H)

**Geometry**:
- Single hexahedron: 1 cm cube
- Material: Steel (B-H curve)
- External field: **H_ext = [0, 0, H_range]** where H_range = 100 to 50000 A/m

**B-H curve** (example):
```python
BH_DATA = [
    [0, 0.0],
    [100, 0.1],
    [500, 0.8],
    [1000, 1.2],
    [5000, 1.7],
    [50000, 2.0],
]
```

**Validation**:
```python
mat = rad.MatSatIsoTab(BH_DATA)
rad.MatApl(hex, mat)

H_range = [100, 500, 1000, 5000, 10000, 50000]
M_computed = []

for H_ext in H_range:
    bkg = rad.ObjBckg(lambda p: [0, 0, H_ext/796000])
    container = rad.ObjCnt([hex, bkg])
    rad.Solve(container, 0.0001, 100, 0)
    M = rad.ObjM(hex)[2]
    M_computed.append(M)

# Compare to B-H curve
# M = B/μ_0 - H
```

**Acceptance**: M vs H follows B-H curve within 5%

**Status**: ✅ Already validated (nonlinear benchmarks complete)

---

## Phase 3: PEEC-MSC Linear Coupling

**Status**: 🔄 After PEEC validation
**Goal**: Couple PEEC conductor + MSC core (linear materials)

### Test 3.1: Loop + Single Hex Core

**Geometry**:
- Circular loop: radius **a = 5 cm**, wire radius **r_wire = 1 mm**
- Hex core: 3 cm cube at loop center
- Material: Linear (μ_r = 500)
- Frequency: 1 kHz
- Excitation: 1 A current

**Expected**:
1. Loop current → field at core
2. Core magnetization → increased loop inductance
3. L_with_core > L_air

**Analytical estimate** (for small core):
```
L_with_core ≈ L_air · (1 + k · μ_r)
```

where **k** is geometric filling factor (≈ V_core / V_loop)

**Validation**:
```python
# Without core
loop_only = rad.PEECLoop([0,0,0], 0.05, [0,0,1], 0.001, 5.8e7, 32)
Z_air = rad.PEECImpedance(loop_only, 1000)
L_air = np.imag(Z_air) / (2*np.pi*1000)

# With core
core = rad.ObjHexahedron(vertices_3cm_cube, [0, 0, 0])
mat = rad.MatLin(500)
rad.MatApl(core, mat)

coupled = rad.CoupledPEECMSC(loop_only, core)
Z_coupled = rad.SolveCoupled(coupled, freq=1000, I_source=1.0)
L_core = np.imag(Z_coupled) / (2*np.pi*1000)

# Compare
ratio = L_core / L_air
print(f"L_air = {L_air*1e9:.2f} nH")
print(f"L_core = {L_core*1e9:.2f} nH")
print(f"Ratio = {ratio:.2f}")
```

**Acceptance**:
- L_core / L_air > 1 (core increases inductance)
- Compare to CplMag

---

### Test 3.2: Loop + Multi-Element Core

**Geometry**:
- Circular loop: radius **a = 5 cm**
- Core: 3×3×3 hex mesh (27 elements), total 3 cm cube
- Material: Linear (μ_r = 500)

**Goal**: Verify multi-element MSC works with PEEC

**Validation**: Compare to Test 3.1 (single element) → should give similar L_core

**Acceptance**: |L_multi - L_single| / L_single < 5%

---

### Test 3.3: Frequency Sweep (Impedance vs Frequency)

**Geometry**: Same as Test 3.1

**Frequency range**: 100 Hz to 100 kHz (logarithmic sweep)

**Expected**:
```
Z(ω) = R(ω) + jω L(ω)
```

- **Low frequency**: R ≈ R_dc, L ≈ L_dc
- **High frequency**: R increases (skin effect), L may decrease (core losses)

**Validation**:
```python
freqs = np.logspace(2, 5, 20)  # 100 Hz to 100 kHz
Z = [rad.SolveCoupled(coupled, freq=f, I_source=1.0) for f in freqs]

R = [np.real(z) for z in Z]
L = [np.imag(z) / (2*np.pi*f) for z, f in zip(Z, freqs)]

# Plot
plt.subplot(2,1,1)
plt.semilogx(freqs, R, 'o-')
plt.ylabel('R (Ω)')
plt.grid(True)

plt.subplot(2,1,2)
plt.semilogx(freqs, np.array(L)*1e9, 'o-')
plt.xlabel('Frequency (Hz)')
plt.ylabel('L (nH)')
plt.grid(True)
plt.savefig('validation_freq_sweep.pdf')
```

**Acceptance**: Smooth curves, consistent with physics

---

### Test 3.4: Comparison with CplMag

**Goal**: Verify PEEC-MSC gives same results as existing CplMag

**Geometry**: Same as CplMag test cases

**Validation**:
```python
# CplMag (existing)
coil = rad.CndLoop([0,0,0], 0.05, [0,0,1], 'r', 0.002, 0.002, 5.8e7, 8, 36)
core = rad.ObjCnt([hex1, hex2, ...])

solver_cplmag = rad.CplMagCreate(coil, core)
rad.CplMagSetFrequency(solver_cplmag, 1000)
rad.CplMagSetMu(solver_cplmag, 500, 0)
result_cplmag = rad.CplMagSolve(solver_cplmag)

# PEEC-MSC (new)
peec_mesh = rad.PEECLoop([0,0,0], 0.05, [0,0,1], 0.001, 5.8e7, 32)
msc_mesh = core  # same core

coupled = rad.CoupledPEECMSC(peec_mesh, msc_mesh)
result_peecmsc = rad.SolveCoupled(coupled, freq=1000, I_source=1.0)

# Compare impedance
error = abs(result_peecmsc - result_cplmag) / abs(result_cplmag) * 100
print(f"Error vs CplMag = {error:.2f}%")
```

**Acceptance**: Error < 5%

---

## Phase 4: PEEC-MSC Nonlinear Coupling

**Status**: ⏳ After Phase 3
**Goal**: Nonlinear core materials (μ_r(H))

### Test 4.1: Loop + Nonlinear Core (B-H Curve)

**Geometry**: Same as Test 3.1, but with nonlinear material

**Material**: Steel (B-H curve)

**Excitation**: Current sweep **I = 0.1 to 10 A**

**Expected**:
- Low current: High μ_r → high inductance
- High current: Core saturates → μ_r drops → inductance drops

**Validation**:
```python
I_range = [0.1, 0.5, 1.0, 2.0, 5.0, 10.0]
L_computed = []

mat = rad.MatSatIsoTab(BH_DATA)
rad.MatApl(core, mat)

for I in I_range:
    Z = rad.SolveCoupled(coupled, freq=1000, I_source=I)
    L = np.imag(Z) / (2*np.pi*1000)
    L_computed.append(L)

# Plot L vs I
plt.plot(I_range, np.array(L_computed)*1e9, 'o-')
plt.xlabel('Current (A)')
plt.ylabel('Inductance (nH)')
plt.grid(True)
plt.savefig('validation_nonlinear_L_vs_I.pdf')
```

**Acceptance**: L decreases with increasing I (saturation effect)

---

### Test 4.2: Convergence with Newton Damping

**Goal**: Verify Newton-Raphson converges for difficult cases

**Geometry**: Loop + highly saturated core

**Validation**:
- Monitor convergence: `|ΔM| / |M|` vs iteration
- Try different damping factors: α = 1.0, 0.5, 0.3

**Acceptance**: Converges within 20 iterations for α = 0.3

---

## Phase 5: Complex Geometries

**Status**: ⏳ Future
**Goal**: Real-world applications

### Test 5.1: Planar Transformer

**Geometry**:
- Primary winding: PCB spiral (PEEC mesh)
- Secondary winding: PCB spiral (PEEC mesh)
- Core: Ferrite E-core (MSC hex mesh from Cubit)

**Validation**: Compare to FEM or measurement

---

### Test 5.2: Wireless Power Transfer (WPT)

**Geometry**:
- Transmitter coil: Litz wire (PEEC mesh)
- Receiver coil: Litz wire (PEEC mesh)
- Ferrite shields: MSC hex mesh

**Frequency**: 6.78 MHz or 13.56 MHz

**Validation**: Compare to measurement or ngbem

---

## Acceptance Criteria

### Error Tolerances

| Test Category | Acceptance Criterion |
|---------------|---------------------|
| **Analytical comparison** | Error < 5% |
| **CplMag comparison** | Error < 5% |
| **FEM comparison** | Error < 10% |
| **Measurement** | Error < 15% |

### Convergence Requirements

| Solver | Max Iterations | Residual |
|--------|---------------|----------|
| **Linear** | < 100 | < 1e-4 |
| **Nonlinear (easy)** | < 20 | < 1e-4 |
| **Nonlinear (hard)** | < 50 | < 1e-3 |

### Energy Conservation

**For lossless case**:
```
|P_source - jω(W_peec + W_msc)| / |P_source| < 5%
```

### Performance Benchmarks

| Problem Size | Target Time |
|--------------|-------------|
| 100 PEEC elements + 100 MSC elements | < 1 s |
| 1000 + 1000 | < 10 s |
| 10000 + 10000 | < 100 s (with HACApK) |

---

## Summary

**Validation strategy**: Bottom-up, analytical first, cross-verification

**Current status**:
- ✅ Phase 2: MSC validated
- 🔄 Phase 1: PEEC validation (next)
- ⏳ Phase 3: Linear coupling
- ⏳ Phase 4: Nonlinear coupling
- ⏳ Phase 5: Complex geometries

**Next action**: Implement and validate **Test 1.1** (circular loop self-inductance)

---

**References**:
1. Grover, F.W. (2004). *Inductance Calculations*. Dover Publications.
2. Dowell, P.L. (1966). "Effects of eddy currents in transformer windings." Proc. IEE.
3. Radia Documentation: `PEEC_MSC_COUPLING.md`
