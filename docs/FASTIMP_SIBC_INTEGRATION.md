# FastImp + SIBC Integration Design Document

**Date**: 2026-01-08
**Status**: Planning

## Overview

This document describes the design for integrating FastImp-based conductor modeling and Nonlocal SIBC (Surface Impedance Boundary Condition) into Radia for wide-band electromagnetic analysis.

## Goals

1. **Coil/Conductor modeling**: Import FastImp formulation for surface current analysis
2. **Magnetic material**: Use existing ELF/Radia MSC formulation
3. **Conductive magnetic material**: Implement Nonlocal SIBC for materials with both conductivity and permeability

## Target Applications

- **Accelerator magnets**: Insertion devices, undulators, wigglers
- **Kicker magnets**: Fast pulsed magnets for beam injection/extraction (eddy currents in vacuum chamber)
- **WPT (Wireless Power Transfer)**: Coil impedance and resonance analysis
- **On-chip inductors**: High-frequency parasitic extraction
- **Transformers**: Core losses and winding impedance
- **Electromagnetic shielding**: Eddy current shielding effectiveness

## Architecture

```
+-------------------------------------------------------------+
|                    Radia Unified Solver                      |
+-------------------------------------------------------------+
|                                                             |
|  +------------------+  +------------------+  +--------------+|
|  |  Coil/Conductor  |  | Magnetic (s=0)   |  | Magnetic     ||
|  |                  |  |                  |  | Conductor    ||
|  |  FastImp         |  |  ELF/Radia MSC   |  | (s!=0,ur!=1) ||
|  |  Surface K, s    |  |  Surface sm, M   |  | Nonlocal SIBC||
|  +--------+---------+  +--------+---------+  +------+-------+|
|           |                     |                   |        |
|           +----------+----------+-------------------+        |
|                      v                                       |
|              +-------------------+                            |
|              |  Coupled Solver   |                            |
|              +-------------------+                            |
|                      |                                       |
|           +----------+----------+                            |
|           v                     v                            |
|     +-----------+         +-----------+                      |
|     |  HACApK   |         | pFFT(MKL) |                      |
|     | Low freq  |         | High freq |                      |
|     +-----------+         +-----------+                      |
+-------------------------------------------------------------+
```

## Three Modules

### 1. Coil/Conductor (FastImp Formulation)

**Target**: Copper coils, aluminum conductors, wiring
**Properties**: s ~ 10^7 S/m, ur = 1

**Unknowns**:
- K: Surface current density [A/m]
- s: Surface charge density [C/m^2]

**Formulation**: FastImp Full-wave IE
```
A = u * integral{ g(r,r') * K dF' }
Phi = (1/e) * integral{ g(r,r') * s dF' }
g = exp(-jkr) / (4*pi*r)
```

**Acceleration**: pFFT with MKL FFT backend (GPL-free)

### 2. Magnetic Material (ELF/Radia MSC)

**Target**: Permanent magnets, ferrite (high resistivity), soft iron (low frequency)
**Properties**: s ~ 0, ur >> 1

**Unknowns**:
- sm: Magnetic surface charge density [Wb/m^2]
- M: Magnetization vector [A/m]

**Formulation**: MSC (existing Radia)
```
H = -(1/4pi) * integral{ sm * (r-r')/|r-r'|^3 dF' }
sm = M . n_hat
```

**Acceleration**: HACApK (ACA+)

### 3. Conductive Magnetic Material (Nonlocal SIBC)

**Target**: Electrical steel sheets, iron yoke (high frequency)
**Properties**: s ~ 10^6 S/m, ur ~ 1000-10000

**Unknowns**:
- K: Surface current density [A/m]
- Internal: Solved by 2D cross-section FEM

**Formulation**: Nonlocal SIBC (Bilicz et al., 2023)
```
Et = Z{n_hat x K}
Z: Nonlocal impedance operator from 2D FEM

2D cross-section problem:
nabla^2 Eu - jw*u(H)*k*Eu = 0  in Omega (cross-section)
dEu/dn = -jw*u*Hv              on Gamma (boundary)
```

**Reference**:
- Bilicz S, Badics Z, Pavo J. "Wide-band nonlocal impedance boundary condition model for high-conductivity regions in integral equation framework", December 2023.

**Acceleration**: Surface integrals with pFFT/HACApK, 2D FEM solved once per cross-section

## Interaction Matrix

```
[Z_cc  Z_cm  Z_cs] [K_coil ]   [V_ext]
[Z_mc  Z_mm  Z_ms] [sm     ] = [H_ext]
[Z_sc  Z_sm  Z_ss] [K_sibc ]   [E_ext]

Z_cc: Coil-Coil (FastImp)
Z_mm: Magnetic-Magnetic (MSC)
Z_ss: SIBC-SIBC (Nonlocal SIBC)
Z_cm, Z_mc: Coil-Magnetic (cross terms)
Z_cs, Z_sc: Coil-SIBC (cross terms)
Z_ms, Z_sm: Magnetic-SIBC (cross terms)
```

## Implementation Phases

### Phase 1a: FastImp Port (FFTW -> MKL)

**Tasks**:
1. Analyze FastImp source code (https://github.com/ediloren/FastImp)
2. Identify FFTW dependencies in pfft++
3. Replace FFTW calls with Intel MKL FFT (DftiCreateDescriptor, etc.)
4. Build with MSVC + MKL (same as Radia)
5. Validate against original FastImp results

**Key Files**:
- `pfft++/src/` - pFFT implementation
- `fastImp/src/surf/formulation.cc` - Core formulation

**Dependencies**: Intel MKL (already used by Radia)

### Phase 1b: MSC Verification

**Tasks**:
1. Verify existing Radia MSC implementation
2. Ensure compatibility with FastImp surface element format
3. Document interface requirements

### Phase 2: Coil + Magnetic Material Coupling

**Tasks**:
1. Define cross-term computation (Z_cm, Z_mc)
2. Implement coupled solver (iterative or direct)
3. Test with simple coil + iron core problem
4. Benchmark against reference solutions

**Cross-term Physics**:
```
Coil -> Magnetic: B field from coil currents induces magnetization
Magnetic -> Coil: H field from magnetization affects coil impedance
```

### Phase 3: Nonlocal SIBC Implementation

**Tasks**:
1. Implement 2D cross-section FEM solver
   - Option A: Use NGSolve (already integrated)
   - Option B: Simple FEM for Helmholtz equation
2. Implement nonlocal SIBC operator Z{.}
3. Couple with FastImp surface formulation
4. Test with conductive magnetic material
5. Validate against 3D FEM (Ansys HFSS or similar)

**2D FEM Choice**:
- Recommend NGSolve (same build configuration as radia_ngsolve)
- Cross-section problem is simple: Helmholtz equation with Robin BC

### Phase 4: Full Integration

**Tasks**:
1. Unified API design
2. Automatic material type detection
3. Frequency sweep support
4. Resonance finding
5. Documentation and examples

### Phase 5: Transient Analysis (Future Work)

**Goal**: Convert frequency-domain Z(s) to time-domain response

**Candidate Methods**:

1. **Cauer Ladder Network (CLN) Method**
   - Z(s) -> continued fraction expansion -> L-C ladder
   - Physically meaningful equivalent circuit
   - Passive and stable by construction
   - **Requirement: Symmetric impedance matrix**
   - Challenge: May fail for coupled systems with non-symmetric cross terms

2. **Arnoldi-based Model Order Reduction** (Likely choice)
   - Krylov subspace projection
   - PRIMA (Passive Reduced-order Interconnect Macromodeling Algorithm)
   - **Works with non-symmetric matrices**
   - Mathematically robust for coupled multi-physics systems
   - Challenge: Less physical interpretation

3. **Vector Fitting**
   - Rational function approximation of Z(s)
   - Widely used in signal integrity
   - Can be converted to state-space or SPICE model

**Open Questions**:
- Which method works best for electromagnetic systems with eddy currents?
- Where does each method break down?
- How to handle nonlinear magnetic materials in time domain?

**Matrix Symmetry Consideration**:
- Single-physics blocks (Z_cc, Z_mm, Z_ss): Symmetric (reciprocity)
- Cross-coupling blocks (Z_cm, Z_mc, etc.): May be non-symmetric
- For coupled multi-physics systems, **Arnoldi-based methods are likely required**

**Note**: This phase requires further research to determine the optimal approach.
The choice between CLN, Arnoldi, or other methods depends on the specific
characteristics of the impedance function Z(s) obtained from the IE solver.

## API Design (Preliminary)

```python
import radia as rad

# Set units
rad.FldUnits('m')

# Create coil (FastImp)
coil = rad.ObjCoil(vertices, conductivity=5.8e7)

# Create magnetic core (MSC)
core_vertices = [...]
core = rad.ObjHexahedron(core_vertices, [0, 0, 0])
mat_core = rad.MatLin(1000)  # ur = 1000
rad.MatApl(core, mat_core)

# Create conductive magnetic shield (SIBC)
shield = rad.ObjConductiveMagnetic(
    vertices,
    conductivity=1e6,
    permeability=1000,
    cross_section='rectangular'
)

# Assemble and solve
assembly = rad.ObjCnt([coil, core, shield])

# Frequency domain analysis
freq = 1e6  # 1 MHz
Z = rad.SolveImpedance(assembly, freq)

# Frequency sweep
freqs = np.logspace(3, 9, 100)  # 1 kHz to 1 GHz
Z_sweep = rad.ImpedanceSweep(assembly, freqs)

# Find resonances
resonances = rad.FindResonances(assembly, freq_range=[1e6, 1e9])
```

## Build Configuration

All components use MSVC + Intel MKL (same as current Radia):

```
Compiler: MSVC (Visual Studio 2022)
Libraries:
  - Intel MKL (BLAS/LAPACK + FFT)
  - NGSolve (for 2D FEM in SIBC, optional)

No new dependencies:
  - FFTW (GPL) is NOT used - replaced by MKL FFT
  - Intel Compiler is NOT used - only MKL library
```

## Performance Considerations

### pFFT vs HACApK Selection

| Frequency Range | Kernel | Recommended |
|-----------------|--------|-------------|
| DC / Low freq | 1/r | HACApK |
| MQS | 1/r, jw/r | HACApK or pFFT |
| Full-wave | exp(-jkr)/r | pFFT |

### Memory Estimates

For N surface elements:
- Dense matrix: O(N^2)
- HACApK: O(N log N)
- pFFT: O(N)

### Expected Performance (from Bilicz paper)

| Problem | IE + SIBC | 3D FEM |
|---------|-----------|--------|
| Loop (720 elem) | 1 s/freq | 46 s/freq |
| Spiral (1800 elem) | 11 s/freq | 65 s/freq |

## References

### Nonlocal SIBC (Primary Reference)

[1] S. Bilicz, Z. Badics, and J. Pávó, "Wide-band nonlocal impedance boundary condition model for high-conductivity regions in integral equation framework," presented at ISEM 2023 (International Symposium on Electromagnetic Fields in Mechatronics, Electrical and Electronic Engineering), December 2023.
- Affiliation: Budapest University of Technology and Economics, Hungary; Tensor Research, LLC, USA
- ORCID: S. Bilicz (0000-0003-4995-6698), Z. Badics (0000-0001-6176-3675), J. Pávó (0000-0002-9501-7176)
- Funding: Hungarian Scientific Research Fund, Grant K-135307

[2] S. Bilicz, Z. Badics, S. Gyimóthy, and J. Pávó, "A Full-Wave Integral Equation Method Including Accurate Wide-Frequency-Band Wire Models for WPT Coils," IEEE Transactions on Magnetics, vol. 54, no. 3, pp. 1-4, March 2018.
- DOI: 10.1109/TMAG.2017.2771366

### FastImp

[3] Z. Zhu, B. Song, and J. K. White, "Algorithms in FastImp: a fast and wide-band impedance extraction program for complicated 3-D geometries," IEEE Transactions on Computer-Aided Design of Integrated Circuits and Systems, vol. 24, no. 7, pp. 981-998, July 2005.
- DOI: 10.1109/TCAD.2005.850814
- Affiliation: Massachusetts Institute of Technology

[4] FastImp source code: https://github.com/ediloren/FastImp (MIT License)

### Related Work

[5] M. Al-Qedra, J. Aronsson, and V. Okhmatovski, "A Novel Skin-Effect Based Surface Impedance Formulation for Broadband Modeling of 3-D Interconnects With Electric Field Integral Equation," IEEE Transactions on Microwave Theory and Techniques, vol. 58, no. 12, pp. 3872-3881, December 2010.

[6] W. C. Gibson, "The Method of Moments in Electromagnetics," Boca Raton: Chapman & Hall/CRC, 2008.

## License Considerations

- FastImp: MIT license (source code)
- FFTW: GPL (NOT used - replaced by MKL)
- Intel MKL: Intel EULA (redistributable, already used by Radia)
- NGSolve: LGPL (dynamic linking OK)

**Result**: No GPL contamination, compatible with Radia's license.
