# Magnetic Aftereffect Research Plan

Research plan for quantitative evaluation of magnetic aftereffect,
conducted at Sugahara Lab (Kindai University) in collaboration with Prof. Hane's group.

This document describes the experimental and theoretical framework for separating
magnetic aftereffect from eddy current effects in soft magnetic materials, and its
connection to Radia's energy-based hysteresis model (see `B_INPUT_PLAY_MODEL.md`).

---

## 1. Starting Point: Fourier Separation of B-H Loops

### Taka's Formula (SA-26-005, IEEJ Static Apparatus / Rotating Machinery Joint Workshop @ Kawagoe, 2026/03/05)

**Paper**: "Proposal of Iron Loss Calculation Formula for PWM Excitation of Electrical Steel Sheets"
**Authors**: Taka Enki (Oita Univ.), Watanabe Masashi et al. (Toyota Industries)

By applying DFT to the magnetic field h under sinusoidal B excitation and using
the substitution `sin(wt) -> b/Bmax`, `cos(wt) -> db/(w*Bmax*dt)`:

```
h = (a1 + a3*b^2 + a5*b^4 + ...) * (db/dt)
  + (c1 + c3*b^2 + c5*b^4 + ...) * b              ... Eq.(4)
```

| Term | Coefficient | Physical meaning |
|------|-------------|-----------------|
| `a_i * (db/dt)` | a1, a3, a5, ... | **Dynamic component** (eddy current + aftereffect) |
| `c_i * b` | c1, c3, c5, ... | **Static component** (hysteresis loss) |

- Coefficients a_i, c_i are computed from harmonic amplitudes H_i_max and phase angles phi_i
- Sinusoidal excitation: odd harmonics only (symmetric loop). Paper considers up to **7th harmonic**
- PWM extension achieves < 1% total iron loss error

### Unverified Points in Taka's Work

- Only **one condition** tested: 50 Hz sinusoidal, Bmax = 1 T
- **No frequency variation** -> frequency invariance of c_i unconfirmed
- **No Bmax variation** -> smoothness of a_i(Bmax), c_i(Bmax) unconfirmed

---

## 2. Core Research Thesis

### Theme

**"Quantitative evaluation of magnetic aftereffect using Fourier-based magnetic property
separation and reinforcement-learning-based adaptive low-frequency measurement"**

### Central Claim: "DC Hysteresis Does Not Exist"

The conventional "static hysteresis" term (c_i * b) may itself contain magnetic aftereffect.

**If c_i(f) does not converge to a constant as f -> 0**, then DC hysteresis is merely
a convenient concept inseparable from magnetic aftereffect.

### Logical Chain

```
Physical claim    -> DC hysteresis doesn't exist (inseparable from aftereffect)
Evidence needed   -> Frequency dependence of c_i(f) at ultra-low frequency
Enabling tech     -> dSPACE + reinforcement learning for adaptive measurement
```

Both physics and engineering are required; neither alone is sufficient.

---

## 3. Separation of Eddy Current and Magnetic Aftereffect

The a_i * (db/dt) term contains both eddy current and aftereffect contributions:

| | Eddy current | Magnetic aftereffect |
|---|---|---|
| Frequency dependence | Proportional to f (classical) / f^1.5 (anomalous) | Depends on domain wall relaxation time (proportional to ln(t)) |
| Low-frequency behavior | Vanishes as f -> 0 | **Persists** |

**Separation method**:
- Measure a_i(f) at multiple frequencies
- Decompose into f-proportional component (eddy current) and residual (aftereffect)
- **Intercept of a_i(f) extrapolated to f -> 0 = aftereffect contribution**
- Lower frequency -> smaller eddy current -> aftereffect more visible, but S/N degrades
  -> reinforcement learning needed

---

## 4. Comparison with Prior Methods

### Morishita's Thesis (2025, Sugahara Lab)

Gaussian basis function H-axis warping for B-H loop smoothing, with warping
parameters w interpolated vs Bmax using PCHIP. 5 Hz measurement, < 1% relative error
(35A360, grain-oriented Si steel, powder core).

| | Morishita's method | Taka's separation + interpolation |
|---|---|---|
| Interpolation target | Entire B-H loop | Fourier coefficients a_i, c_i |
| Interpolation difficulty | Steep changes -> Gaussian warping needed | 1D functions -> PCHIP/spline sufficient |
| Aftereffect separation | Not possible (mixed) | **Possible** |
| Frequency extension | Independent re-measurement needed | a_i frequency interpolation |

### Doshisha 2-Frequency Method

- Coarse separation based on f^2 dependence vs independence (iron loss "values")
- Taka's method: per-harmonic separation (B-H loop "structure") -> far more information

---

## 5. Open Technical Issues

### Smoothness of Fourier Coefficients vs Bmax

- Taka tested Bmax = 1 T only -> completely unexplored territory
- B-H loops vary continuously with Bmax, so a_i(Bmax), c_i(Bmax) should be continuous
- **Near saturation knee**: rate of change may be large
- Higher harmonics are more sensitive to fine B-H shape changes -> noise-like behavior possible
- **First task**: Plot a_i(Bmax), c_i(Bmax) experimentally (Muramatsu's initial work)
- If not smooth, Morishita's Gaussian warping may be needed as fallback

### Harmonic Truncation Order

- Taka uses 7th harmonic. Sinusoidal excitation -> odd harmonics only
- Strongly nonlinear materials may need higher orders
- Start with 7th, evaluate truncation error, adjust

### Ring Core Radial Non-Uniformity

- Ring specimen: H = N1*I/(2*pi*r), non-uniform in radial direction
- Mean path length approximation fails for wide rings
- Errors prominent at low Bmax (Rayleigh regime) and near saturation
- **Must correct this to detect the subtle aftereffect signal**
- Without correction, geometric B-H distortion is indistinguishable from aftereffect

---

## 6. Matsuo's Ring Core Correction Method

### Reference

- **Author**: Matsuo Tetsuji (Kyoto University)
- **Paper**: "Accurate Identification of BH-Loops From Measurement Using Wide Ring Specimen"
- **Journal**: IET Science, Measurement & Technology, 2025

### Method Overview

Analytically corrects measured i-phi loops for radial field non-uniformity to
reconstruct the true B-H loop.

Core formula (initial magnetization curve):
```
f(Ha)/Ha = mu_1 + sum_{n=0}^inf F(c^n * i)
F(i) = dphi/di - phi/i,  c = a/b (inner/outer radius ratio)
```

Hysteresis loop extension (alpha method):
```
f(alpha*Ha, Ha)/Ha = mu_alpha_1 + sum_{n=0}^inf F_alpha(c^n * im)
```
Alpha varies from -1 to 1 to reconstruct the full B-H loop.

### Characteristics

- Analytical formulation (simpler than prior numerical methods)
- Accurate even for wide rings (c ~ 0.5)
- Compatible with Play model for hysteresis loop correction

---

## 7. Measurement System Design

### Hardware

| Component | Model | Role |
|-----------|-------|------|
| Real-time controller | **dSPACE** | High-precision A/D, real-time control |
| Bipolar power supply | **NF 4520** | Low-frequency, high-current excitation |
| Multi-function supply | **NF 4505** | CC mode (relatively clean excitation) |

### NF 4520 Control Strategy

CC mode (current control) produces PWM-like noise at low currents.
-> Use **CV mode (voltage control) + dSPACE feedback**:

1. dSPACE generates target b(t) waveform
2. Secondary coil dB/dt signal integrated to obtain B
3. Voltage command based on B_target - B_measured error

### Two-Stage Reinforcement Learning

**Stage 1: B waveform tracking control (excitation control)**

- Nonlinear material -> fixed PI gains insufficient across conditions
- SAC learns adaptive voltage commands for nonlinear materials
- Reward: minimize tracking error between B_target(t) and B_measured(t)
- Clean sinusoidal B -> better Fourier separation -> better aftereffect detection
- This control work is publishable independently

**Stage 2: Adaptive measurement strategy (S/N improvement)**

- Address S/N degradation at low frequency
- Reward: estimation accuracy of a_i, c_i coefficients

---

## 8. Research Steps

### Step 0: Theory Understanding (Miura & Muramatsu, shared)

Reproduce Taka's Eq.(4) derivation using MATLAB Symbolic.

### Step 1: Measurement System (Miura)

dSPACE + NF 4520 CV mode + RL-based B waveform tracking.
Material: 35A360 (same as Morishita's thesis -> direct comparison possible).

### Step 2: Multi-Frequency / Multi-Bmax Measurement (Miura)

- Frequencies: 5 Hz, 50 Hz, 200 Hz, 1 kHz
- Multiple Bmax values

### Step 2.5: Ring Core Correction (Muramatsu)

Apply Matsuo's correction to measured i-phi loops.
Remove radial non-uniformity to reconstruct true B-H loops.

### Step 3: Fourier Separation and Verification (Muramatsu)

- Compute a_i, c_i from corrected B-H loops at each condition
- Verify frequency invariance of c_i:
  - If invariant -> separation validated
  - If **not invariant** -> evidence for magnetic aftereffect in "static" term
- Check smoothness of a_i(Bmax), c_i(Bmax)
- Separate eddy current and aftereffect from a_i(f)

### Step 4: Ultra-Low Frequency Challenge (Miura & Muramatsu)

Below 1 Hz -> S/N problems -> Stage 2 reinforcement learning needed.

---

## 9. Data Processing Pipeline

```
Ring core measurement (Miura)
  |
  v
Raw i-phi loops at multiple (f, Bmax)
  |
  v
Matsuo's radial correction (Muramatsu)        <- Remove geometric artifact
  |
  v
True B-H loops
  |
  v
Taka's Fourier separation                     <- Extract a_i(f, Bmax), c_i(f, Bmax)
  |
  v
Frequency extrapolation                       <- a_i(f->0) = aftereffect
  |                                              c_i(f) dependence = "DC hyst" test
  v
Aftereffect parameters
  |
  v
Radia energy-based model (Section 12 of B_INPUT_PLAY_MODEL.md)
  +-- Arrhenius rate: dJ_k/dt ~ exp(-Delta_U_k / kT)
  +-- Pinning landscape from measured shape functions
```

---

## 10. Connection to Radia's Energy-Based Hysteresis

### Current State

Radia implements the Francois-Lavet / Egger energy-based vector hysteresis model
with O(K) forward and inverse operators (see `B_INPUT_PLAY_MODEL.md`). This model
is **rate-independent** -- it captures static hysteresis but not time-dependent effects.

### Future Extension: Aftereffect in Energy Framework

The energy-based model provides a natural framework for magnetic aftereffect:

| Component | Rate-independent (current) | With aftereffect (planned) |
|-----------|--------------------------|---------------------------|
| State evolution | J_k jumps to energy minimum | J_k evolves via thermal activation |
| Pinning | chi_k gives sharp threshold | chi_k defines barrier height |
| Time scale | Instantaneous | Arrhenius: tau_k = f_0^{-1} exp(Delta_U_k / kT) |
| Congruency | H-axis congruent (exact) | Modified by relaxation (non-congruent) |

### Parameter Identification from Measurement

The Fourier separation provides a direct path to aftereffect parameters:

1. **a_i(f -> 0) intercept** -> rate-dependent component magnitude
   -> constrains attempt frequency f_0 and barrier distribution
2. **c_i(f) dependence** -> how "static" hysteresis changes with rate
   -> constrains barrier height distribution (related to chi_k)
3. **Temperature dependence** (if measured) -> directly gives barrier heights Delta_U_k

---

## 11. Team and Timeline

### Students (current 3rd year -> 2026 4th year)

| Student | Role | Focus |
|---------|------|-------|
| **Miura** | Measurement | dSPACE + CV mode control, RL excitation, data acquisition |
| **Muramatsu** | Theory | MATLAB Symbolic derivation, Fourier coefficient analysis |

### Research Component Dependencies

```
Taka's theory (SA-26-005)        -> Fourier separation foundation
Matsuo's correction (IET 2025)   -> Measurement accuracy foundation
Ito's SAC research (SA-26-019)   -> Reinforcement learning expertise
Hane Lab collaboration           -> Hysteresis measurement data
Morishita's thesis (2025)        -> Bmax interpolation prior work
dSPACE + NF 4520/4505           -> Real-time measurement platform
```

### Publications

- Cite Taka (SA-25-001, 2025; SA-26-005, 2026)
- Cite Matsuo (IET Sci. Meas. Technol., 2025)
- Can develop independently; joint work with Taka/Matsuo is an option

---

## 12. Paper Concept: Unified Aftereffect Identification and Simulation

### Three-Layer Novelty

The research combines three independently novel contributions into a closed loop:

```
Layer 1: Measurement           Layer 2: Model            Layer 3: Identification
─────────────────              ────────────               ──────────────────────
Taka's Fourier separation      Energy-based hysteresis    Fourier coefficients
+ Matsuo's ring core           U_k, chi_k from B-H       -> energy parameters
  correction                   + Arrhenius aftereffect    No assumed parameters
      |                              |                          |
      v                              v                          v
  a_i(f,Bmax), c_i(f,Bmax)    dJ_k/dt = f(U_k, chi_k, T)    f_0 from a_i(f) slope
                                                               Delta_U from chi_k
      |                              |                          |
      +──────────────────────────────+──────────────────────────+
                                     |
                                     v
                              Closed validation loop:
                              Model predicts B-H at new (f, Bmax, T)
                              -> Compare with measurement
```

### Key Claim

**"DC hysteresis does not exist as an independent phenomenon"**
-- it is inseparable from magnetic aftereffect.

Evidence: c_i(f) (the "static" Fourier component) shows frequency dependence
even at ultra-low frequencies, meaning the static B-H loop itself changes
with excitation rate.

### Why This Is a Strong Paper

1. **Measurement + Model + Identification form a closed loop** -- each validates the others
2. **No assumed parameters**: f_0 is measured (not assumed from literature)
3. **Physical claim with experimental evidence**: c_i(f) non-convergence
4. **Tool exists**: Radia's energy-based model already implements the rate-independent
   part; aftereffect extension is a natural addition (rate equation on existing U_k, chi_k)
5. **Reviewer-proof**: hard to attack because the loop is self-consistent

---

## 13. Accelerator Application

### Why Aftereffect Matters for Accelerators

Accelerator magnets require ppm-level field stability. Magnetic aftereffect
causes time-dependent field drift that is currently not predictable by
simulation tools.

```
Accelerator operation:
  Beam energy stabilization -> field stability at ppm level
  Post-magnetization yoke  -> M(t) relaxes as ln(t) -> beam orbit drift
  Temperature fluctuation  -> Delta_U/(k_B*T) changes -> field reproducibility
```

### Radia + Aftereffect for Accelerator Design

| Design Task | Current Radia | With aftereffect extension |
|-------------|--------------|---------------------------|
| PM undulator field quality | Static M distribution | **M(t) drift after installation** |
| Iron yoke magnetization | Reversible B-H curve | **Irreversible M after magnet cycling** |
| Beamline field stability | Snapshot only | **Time evolution over hours/days** |
| Thermal compensation | Not possible | **T-dependent M(t) at operating conditions** |
| Accelerated aging test design | Empirical | **Physics-based: barrier height -> lifetime** |

### Why Radia (Not FEM)

Radia's BEM is uniquely suited for this combination:

1. **No air mesh**: field evaluation anywhere along beamline without meshing air
2. **Open boundary**: natural for accelerator geometries (long beamlines)
3. **Far-field accuracy**: exact analytical integrals, no truncation
4. **Hysteresis + BEM**: FEM hysteresis exists, but FEM + open boundary + far-field
   is extremely difficult (requires PML/infinite elements)

The combination of energy-based aftereffect + BEM open boundary is not
available in any existing simulation tool.

### Target Venues

- **IPAC** (International Particle Accelerator Conference)
- **IEEE Trans. Magn.** (modeling + identification)
- **Compumag / CEFC** (computational electromagnetics)
- **J. Synchrotron Rad.** (if undulator-specific application)

---

**Last Updated:** 2026-03-06
