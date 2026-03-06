# Energy-Based Vector Hysteresis Model

Radia's hysteresis module implements the Francois-Lavet / Egger energy-based formulation,
providing both forward (H -> B) and inverse (B -> H) operators with O(K) complexity.

**References**:
1. Egger, Engertsberger, Schafelner: "Efficient evaluation of forward and inverse energy-based magnetic hysteresis operators", MAGCON-25-07-0171 (Compumag 2025)
2. Francois-Lavet et al.: J. Comp. Appl. Math., 246:243-250, 2013

---

## 1. Why Energy-Based (Not Raw Play)

| Criterion | B-input Play (Matsuo) | Energy-Based (Egger) |
|-----------|----------------------|----------------------|
| Thermodynamic consistency | Sacrificed for efficiency | **Strictly satisfied** |
| Iron loss definition | Approximate (loop area) | **Direct from local energy dissipation** |
| Inverse operator | Requires fminsearch | **O(K) via Schur complement** |
| Parameters | Shape functions + eta | U_k (internal energy) + chi_k (pinning) |
| Congruency | **H-axis congruent** (see Section 2-5) | **Structurally guaranteed** |
| Vector extension | 2D -> 3D natural | **Arbitrary dimension (d <= 3)** |

The B-input Play shape functions f_k are reused as U_k' (derivative of internal energy),
so measured data is fully preserved.

---

## 2. Mathematical Formulation

### 2-1. Gibbs Free Energy

```
G*(B, {J_k}, {J_{k,p}}) = W_field + sum_k W_int,k + sum_k D_k
```

| Term | Formula | Physical meaning |
|------|---------|-----------------|
| W_field | (nu_0/2) \|B - sum J_k\|^2 = (mu_0/2)\|H\|^2 | Vacuum field energy |
| W_int,k | U_k(\|J_k\|) | Anisotropy + exchange (reversible) |
| D_k | chi_k \|J_k - J_{k,p}\|_eps | Pinning dissipation (irreversible) |

### 2-2. Internal Energy

Analytical form (per partial polarization k):

```
U_k(r) = -(2 A_{s,k} J_{s,k}) / pi * log(cos(pi/2 * r / J_{s,k}))
U_k'(r) = A_{s,k} * tan(pi/2 * r / J_{s,k})
```

Parameters:
- `A_{s,k}`: saturation slope (controls initial permeability)
- `J_{s,k}`: saturation polarization per operator (total J_sat = sum J_{s,k})
- `chi_k`: pinning strength [A/m] (chi_k = 0 = reversible, chi_k large = wide hysteresis)

### 2-3. Regularization

Non-smooth |x| replaced by smooth approximation:

```
|x|_eps = sqrt(|x|^2 + eps)
```

- eps > 0 makes everything C-infinity -> standard Newton
- Error: |J_k - J_k^eps| = O(sqrt(eps)) -> eps=1e-8 is sufficient
- Forward and inverse are **exactly invertible** under regularization (Assertion 6)

### 2-4. Relationship to Play Model

| B-input Play | Energy-based | H-input Play |
|-------------|--------------|-------------|
| eta_k (B-space threshold) | **chi_k = eta_k** (pinning) | zeta_k (H-space threshold) |
| f_k(\|p\|) (B -> H shape function) | **U_k' = f_k** | g_k = (U_k')^{-1} |

Key identity: **f_k (B-input shape function) = U_k' (derivative of internal energy)**.

Conversion from B-input Play data:
```
U_k(r) = integral_0^r f_k(s) ds
chi_k = eta_k
```

### 2-5. Congruency Property (合同性)

**B-input Play model satisfies H-axis congruency** (H軸方向の合同性).

The Play operators work in B-space with thresholds eta_k defined in B. For minor loops
with the same B-amplitude at different DC bias points:

```
Same delta-B  →  Same delta-p_k  →  Same H output shape
```

This is because the play operator state update depends only on the B excursion:
```
p_k = B - eta_k * (B - p_k_prev) / max(eta_k, |B - p_k_prev|)
```
When |B - B_dc| = const (same amplitude), the delta-p_k trajectory is identical
regardless of B_dc, producing congruent H output loops.

| Model | Input | Output | Congruency Direction |
|-------|-------|--------|---------------------|
| H-input Play | H | B (M) | Not congruent (non-physical) |
| **B-input Play** | **B** | **H** | **H-axis congruent** |
| Energy-based (Egger) | H <-> B | both | Structurally guaranteed |

**Verification procedure**: Drive the material at different DC bias points with the
same AC amplitude, compare the resulting H minor loops. They should be congruent
(same shape, only translated in B).

```python
import radia as rad
import numpy as np

MU_0 = 4e-7 * np.pi

rad.UtiDelAll()
mat = rad.MatEnergyHysteresis(K, As, Js, chi, eps)

H_dc_list = [0, 200, 500]  # Different DC bias [A/m]
H_ac = 100.0                # Same AC amplitude

minor_loops = {}
for H_dc in H_dc_list:
    # Drive to bias point
    for H_val in np.linspace(0, H_dc, 50):
        rad.MatMvsH(mat, 'm', [H_val, 0, 0])

    # Trace minor loop
    t = np.linspace(0, 2 * np.pi, 200)
    H_loop = H_dc + H_ac * np.sin(t)
    B_loop = []
    for H_val in H_loop:
        M = rad.MatMvsH(mat, 'm', [H_val, 0, 0])
        B_loop.append(MU_0 * (H_val + M[0]))
    minor_loops[H_dc] = (H_loop, np.array(B_loop))

# Compare: H_loop shapes at different bias should be congruent
# (same delta-B produces same delta-H)
```

**Collaboration note**: Experimental congruency measurements are ongoing at
Sugahara Lab (Kindai University) in collaboration with Prof. Hane's group,
including amorphous materials where non-congruency has been observed.

---

## 3. Forward Operator: H -> B, O(K)

Each J_k is solved **independently** (H given, no coupling):

```
J_k = argmin_J  U_k(J) - <H, J> + chi_k |J - J_{k,p}|_eps
B = mu_0 H + sum_k J_k
```

Newton's method with Armijo backtracking, 3-5 iterations to tol=1e-12.

**This is the operator used in Radia's MMM/MSC solver** (H-field formulation).

### Hessian Structure

```
H_k = hess(U_k) + chi_k * hess(|.|_eps)

hess(U_k)(J) = (U_k'' - U_k'/|J|) * e e^T + (U_k'/|J|) * I    where e = J/|J|
hess(|x|_eps) = (I - x x^T / |x|_eps^2) / |x|_eps
```

At J=0: `H_k = U_k''(0) * I + chi_k / sqrt(eps) * I` (well-conditioned).

### Jacobian dB/dH (for Newton on Field Problem)

```
dB/dH = mu_0 I + sum_k H_k^{-1}
chi_diff = trace(dB/dH) / (d * mu_0) - 1
```

Available via `jacobian_dBdH()`. No separate derivative computation needed.

---

## 4. Inverse Operator: B -> H, O(K) via Schur Complement

All J_k are **coupled** through sum J_k:

```
min_{J_k} (nu_0/2)|B - sum J_k|^2 + sum [U_k(J_k) + chi_k |J_k - J_{k,p}|_eps]
H = nu_0 (B - sum J_k)
```

### Schur Complement Trick (Egger et al.)

Full Newton system has K blocks of size d. Introduce auxiliary variable delta = sum dJ_k:

```
H_k_priv = hess_Uk + chi_k * hess_norm     (d x d, private to each k)
Full block: (H_k_priv + nu_0 I) dJ_k + nu_0 sum_{l!=k} dJ_l = -g_k

Rearrange: H_k_priv * dJ_k + nu_0 * delta = -g_k
  => dJ_k = H_k_priv^{-1} (-g_k - nu_0 delta)

Sum over k: [I + nu_0 sum H_k_priv^{-1}] delta = -sum H_k_priv^{-1} g_k
```

The Schur matrix `S = I + nu_0 sum H_k_priv^{-1}` is d x d (d=2 or 3).
Total cost: K inversions of d x d + one d x d solve = **O(K)**.

### Performance (Egger Table I, K=20, 500 time steps)

| Method | Time | Iterations/step |
|--------|------|-----------------|
| Forward (Newton) | 10 ms | 6.3 |
| **Inverse (Schur)** | **10 ms** | **6.4** |
| Inverse (Standard) | 40 ms | 6.4 |

K=100: Standard=1071 ms vs **Efficient=59 ms** (18x speedup).

---

## 5. Implementation (C++ Material Class)

### Source Files

| File | Description |
|------|-------------|
| `src/core/rad_material_def.h` | `radTEnergyHysteresisMaterial` class definition (Type 5) |
| `src/core/rad_material_impl.cpp` | Forward/inverse operators, Newton+Armijo, Jacobian |
| `src/core/rad_relaxation_methods.cpp` | Solver integration (`UpdateChiAndCheckConvergence`) |
| `src/radia/hysteresis_io.py` | .hys/.mat readers, Play-to-Energy conversion |
| `examples/hysteresis/verify_cpp_hysteresis.py` | C++ verification suite (4 tests) |

### Python API

```python
import radia as rad

# Create material from analytical parameters
mat = rad.MatEnergyHysteresis(K=10, As=As, Js=Js, chi=chi, eps=1e-8)

# Or from .hys file via hysteresis_io
from hysteresis_io import hys_to_radia
params = hys_to_radia('material.hys', K=20)
mat = rad.MatEnergyHysteresis(**params)

# Apply to element and solve
iron = rad.ObjRecMag([0,0,0], [0.01,0.01,0.01], [0,0,0])
rad.MatApl(iron, mat)
rad.Solve(container, 0.001, 100, 1)

# Direct M(H) evaluation
M = rad.MatMvsH(mat, 'm', [500, 0, 0])  # returns M in A/m
```

### I/O Utilities

```python
from hysteresis_io import load_hys, build_shape_functions, hys_to_radia

# One-step: .hys file -> ready-to-use model
model = hys_to_radia('material.hys', K=20, eps=1e-8, d=3)

# Manual pipeline
loops = load_hys('material.hys')
eta, f_k_tables, Bplay = build_shape_functions(loops)
model = EnergyBasedHysteresis.from_shape_functions(eta, f_k_tables)

# MATLAB .mat format (Potter-Schmulian)
from hysteresis_io import mat_to_radia
model = mat_to_radia('B_input.mat')
```

### .hys File Format (JMAG Compatible)

```
Line 1: Model type ("Jiles Atherton")
Line 2: 0  nx  ny          (nx=points/loop, ny=number of loops)
Line 3: 1  nx  0           (flags)
Line 4: Units ("tesla;A/m" or "A/m;tesla")
Line 5+: col1  col2        (nx * ny data rows, descending branch)
```

### Verification Results (Phase 1)

All 5 tests PASS:

| Test | Description | Result |
|------|-------------|--------|
| Forward-Inverse roundtrip | forward(H) -> B -> inverse(B) -> H' | \|H-H'\|/\|H\| < 3e-10 |
| 1D B-H loop | Sinusoidal H drive, hysteresis detection | Width = 0.44 T |
| O(K) scaling | K = 5, 10, 20, 50 timing | 8.5x for 10x K ratio |
| Regularization | eps = 1e-4 to 1e-12 convergence | B converges to limit |
| 2D rotating field | Circular H trajectory | Smooth \|B\| trajectory |

### Key Bug Fix: Armijo Condition in Inverse

The inverse Newton with Schur complement requires the **directional derivative**
for the Armijo line search condition:

```python
# WRONG (stalls after 2 iterations):
if G_trial <= G_curr - 1e-4 * tau * max_grad:   # max_grad is not the directional derivative

# CORRECT:
dir_deriv = sum(dot(g_k, dJ_k) for k in range(K))
if G_trial <= G_curr + 0.1 * tau * dir_deriv:    # proper Armijo with directional derivative
```

---

## 6. Data Flow: Measured Data -> Radia Solver

```
[Measured data]                       [Model]                    [Radia]

B-input Play fit               EnergyBasedHysteresis
(.hys, .mat, CSV)                   (As, Js, chi)
       |                              |
       v                              |
  load_hys() / load_mat()             |
       |                              |
  build_shape_functions()             |
       |                              |
  convert_play_to_energy()            |
       |                              |
       +---> As[], Js[], chi[] -----> |
                                      |
                                      +-- forward(H): O(K)  <-- Radia H
                                      |   Each J_k independent    (recommended)
                                      |
                                      +-- inverse(B): O(K)  <-- (for A-formulation)
                                      |   Schur complement
                                      |
                                      +-- jacobian_dBdH()   <-- (for Newton on field)
```

---

## 7. MMM/MSC C++ Integration Plan (Phase 2-4)

### Integration Architecture

```
Current flow (MatSatIsoTab):
  AutoRelax_Unified()
    +-- SolveLinearStep()              ... LU/BiCGSTAB/HACApK
    +-- ComputeActualHFieldFromSigma()
    +-- UpdateChiAndCheckConvergence()
          +-- NonlinMater->ComputeChiDualMethod(H_mag, mu_old)
                +-- B-H curve interpolation -> chi = B/(mu_0*H) - 1

New flow (EnergyHysteresis):
  AutoRelax_Unified()                  ... loop structure unchanged
    +-- SolveLinearStep()              ... unchanged
    +-- ComputeActualHFieldFromSigma() ... unchanged
    +-- UpdateChiAndCheckConvergence()
          +-- HystMater->ComputeChiFromHysteresis(H_vec, elem_idx)
                +-- Forward(H_vec) -> B_vec -> chi = |B|/(mu_0*|H|) - 1
                    + J_k state update
```

### Integration Points (4 files)

#### (A) Material Class: `radTEnergyHysteresisMaterial` in `rad_material_def.h`

```cpp
class radTEnergyHysteresisMaterial : public radTMaterial {
    int m_K;                              // Number of play operators
    std::vector<double> m_As, m_Js, m_chi;
    double m_eps;
    std::vector<TVector3d> m_Jk_prev;     // Current J_{k,p}
    std::vector<TVector3d> m_Jk_pinning;  // Pinning reference
    std::vector<TVector3d> m_Jk_current;  // Last computed J_k

    int Type_Material() { return 5; }     // New type ID

    // Core operators
    TVector3d Forward(const TVector3d& H);
    TVector3d Inverse(const TVector3d& B);

    // Solver interface
    double ComputeChiFromH(const TVector3d& H);
    double ComputeChiDualMethod(double H_mag, double mu_old, double relax);
    double ComputeDifferentialChi(double H_mag);
};
```

#### (B) `UpdateChiAndCheckConvergence()` Change

Add `dynamic_cast<radTEnergyHysteresisMaterial*>` branch in
`rad_relaxation_methods.cpp:706`. Existing `ComputeChiDualMethod` path unchanged.

#### (C) `NonlinearContext` Extension

Add `radTEnergyHysteresisMaterial* hystMaterial` and `hystStatesSaved` flag
to `rad_relaxation_methods.h`.

#### (D) Python Binding

```python
mat_hyst = rad.MatEnergyHyst(K=20, As=[...], Js=[...], chi=[...], eps=1e-8)
rad.MatApl(element, mat_hyst)
rad.Solve(assembly, 0.0001, 100, 0)
rad.HystCommit(mat_hyst)
```

### Quasi-Static State Management

```
Step n -> n+1:
  SaveState()                    // Freeze J_{k,p}
  for iter = 0, 1, 2, ...:      // Nonlinear iteration
    RestoreState()               // Rewind J_{k,p} at each retry
    SolveLinearStep()
    ComputeH()
    Forward(H) -> B, {J_k}      // H -> B, update J_k
    chi = |B|/(mu_0*|H|) - 1
    if converged: break
  CommitState()                  // Confirm converged J_k as new J_{k,p}
```

### MMM vs MSC: Same Constitutive Law

| | MMM (tet, 3 DOF) | MSC (hex, 5/6 DOF) |
|--|---|---|
| DOF | Mx, My, Mz | sigma_1, ..., sigma_6 |
| H source | FlatField array | NewFieldArray |
| Forward(H) | **Same**: B = mu_0 H + J -> chi | **Same** |
| J_k state | Per-element 3D vector x K | **Same** |

### Memory Estimate

```
Per element: K x 3 x 8 bytes x 2 (prev + saved) = 48K bytes
K=20, N=10,000:  48 x 20 x 10,000 = 9.6 MB
K=50, N=100,000: 48 x 50 x 100,000 = 240 MB
```

---

## 8. Magnetization (Pulse Magnetization) Application

With hysteresis in the MMM/MSC solver, the following become possible:

| Capability | Current (MatSatIsoTab) | With Hysteresis |
|------------|----------------------|-----------------|
| Magnetization process | Not possible (reversible B-H) | **Possible** (irreversible Play) |
| Remanence prediction | Manual `ObjHexahedron(v, [0,0,Mr])` | **Automatic** (M at H=0) |
| Partial magnetization | Not possible | **Possible** (pulse-dependent) |
| Demagnetization risk | Not possible | **Possible** (irreversible demagnetization) |

### API Sketch

```python
import radia as rad
from hysteresis_io import hys_to_radia

# 1. Create hysteresis material from .hys data
model = hys_to_radia('NdFeB.hys', K=20)
mat = rad.MatEnergyHyst(model.K, model.As, model.Js, model.chi)

# 2. Unmagnetized blank
magnet = rad.ObjHexahedron(vertices, [0, 0, 0])  # M=0
rad.MatApl(magnet, mat)

# 3. Magnetization pulse
coil = rad.ObjBckg(lambda p: magnetizing_field(p, I_peak))
assembly = rad.ObjCnt([magnet, coil])
rad.Solve(assembly, 0.0001, 100, 0)
rad.HystCommit(mat)

# 4. Remove external field -> remanence
coil_off = rad.ObjBckg(lambda p: [0, 0, 0])
assembly_off = rad.ObjCnt([magnet, coil_off])
rad.Solve(assembly_off, 0.0001, 100, 0)
rad.HystCommit(mat)

# 5. Remanent magnetization pattern
B_ext = rad.Fld(magnet, 'b', observation_point)
```

---

## 9. Anisotropy Extension Roadmap

| Step | Internal energy | Pinning | Use case |
|------|----------------|---------|----------|
| **Current** | U_k(\|J\|) isotropic | scalar chi_k | Isotropic soft iron |
| **Step 1** | U_k(\|J\|) isotropic | tensor X_k | Directional pinning |
| **Step 2** | U_k(\|J\|) + W_a(J, e) | tensor X_k | Grain-oriented Si steel |
| **Step 3** | Polycrystalline averaging | ODF-based | Textured materials |

Step 1 (tensor pinning) is simplest and most impactful:
- Shape function data unchanged
- chi_k becomes a positive-definite matrix X_k
- Newton Hessian changes from scalar to 3x3 tensor

Egger's paper (Assumption 1 note):
> "anisotropic variants |chi_k(J_k - J_{k,p})| with matrix-valued parameters chi_k could be considered"

---

## 10. Unified Picture

```
                    Gibbs Free Energy G*(B, {J_k})
                    +-------------------------------+
                    | (mu_0/2)|H|^2 + sum U_k + sum D_k |
                    +---------------+---------------+
                                    |
                    +---------------+---------------+
                    |               |               |
                    v               v               v
            [B-input Play]  [Egger inverse]  [Egger forward]
            Each p_k indep  Schur O(K)       Each J_k indep O(K)
            B -> H (approx) B -> H (exact)   H -> B (exact)
            Measured fit    Energy exact      Radia integration
                    |               |               |
                    +-------+-------+-------+-------+
                            |               |
                            v               v
                    [Isotropic]     [Anisotropic]
                    U_k(|J|)       U_k(|J|) + W_a
                    scalar chi_k   tensor X_k
```

**Unifying key**: f_k (B-input shape function) = U_k' (internal energy derivative).
Through this identity, measured data -> energy parameters -> forward/inverse/anisotropy
all derive from **one energy function G***.

---

## 11. Implementation Phases

| Phase | Content | Files | Status |
|-------|---------|-------|--------|
| **Phase 1** | Python prototype (reference, now removed) | ~~`energy_hysteresis.py`~~ | **DONE** |
| **Phase 2** | C++ material class (Type 5) | `rad_material_def.h`, `rad_material_impl.cpp` | **DONE** |
| **Phase 3** | Solver integration | `rad_relaxation_methods.cpp` | **DONE** |
| **Phase 4** | Python binding | `radia_pybind.cpp`, `radentry.cpp/h` | **DONE** |
| **Phase 5** | Play data -> parameter conversion | `hysteresis_io.py` | **DONE** |
| **Phase 6** | Validation (C-type electromagnet + magnetization) | `examples/hysteresis/` | Planned |

C++ is **1000x+ faster** than the Python prototype (0.07ms vs 78ms for K=50).

---

## 12. Future: Magnetic Aftereffect (磁気余効)

The Sugahara laboratory (Kindai University) plans to extend the energy-based
hysteresis model to include **magnetic aftereffect** (magnetic viscosity).

### Physical Background

Magnetic aftereffect is the time-dependent relaxation of magnetization after
a change in applied field. It is an intrinsic material property originating
from thermal activation over pinning energy barriers -- distinct from eddy
current effects (which are extrinsic/geometric).

After a field change, M(t) relaxes logarithmically:
```
M(t) = M_0 + S * ln(t / t_0)
```
where S is the magnetic viscosity coefficient and t_0 is a microscopic time constant.

### Energy-Based Formulation

The pinning potential U_k naturally provides the energy barrier landscape for
thermal activation. Each partial polarization J_k can evolve via Arrhenius-type
rate equations:

```
dJ_k/dt = f_0 * exp(-Delta_U_k(J_k) / (k_B * T)) * direction
```

where:
- `f_0`: attempt frequency (~10^9 Hz for domain wall processes)
- `Delta_U_k(J_k)`: energy barrier from U_k + chi_k pinning
- `k_B * T`: thermal energy (~4.1e-21 J at 300K)

The energy-based framework is well-suited for this extension because:
1. Energy barriers are **already computed** (U_k and chi_k define the landscape)
2. Thermodynamic consistency is **preserved** (rates satisfy detailed balance)
3. The O(K) Schur complement structure is **unchanged** (each J_k still evolves independently in the aftereffect regime)

### Target Applications

| Application | Description |
|-------------|-------------|
| Post-magnetization relaxation (着磁後の緩和) | M decay after pulse magnetization |
| Long-term stability of PM (経年減磁予測) | Irreversible flux loss over years |
| Temperature-dependent demagnetization (温度減磁) | Accelerated aging at elevated T |
| Accommodation / reptation | Minor loop drift under repeated cycling |

### Accelerator Magnet Applications

Radia originates from accelerator magnet design (ESRF). Magnetic aftereffect
directly impacts accelerator operation where ppm-level field stability is required:

| Capability | Current (rate-independent) | With aftereffect |
|-----------|--------------------------|-----------------|
| PM undulator design | Static M only | **Long-term flux loss prediction** |
| Iron yoke design | Reversible B-H only | **Post-magnetization drift M(t)** |
| Field stability | Cannot evaluate | **ppm-level temporal variation** |
| Temperature compensation | Manual estimation | **T-dependent M(t) simulation** |

Radia's BEM (no air mesh) computes fields along long beamlines with natural
open boundary -- this combination (hysteresis + open boundary + far-field
accuracy) is extremely difficult with FEM.

### Parameter Identification via Fourier Separation

The energy-based model has a key advantage: **all aftereffect parameters can
be identified experimentally** using Taka's Fourier separation method
(see `MAGNETIC_AFTEREFFECT_RESEARCH.md`):

```
Measured B-H loops at multiple frequencies
  |
  v
Taka's Fourier separation: h = (a_i * db/dt) + (c_i * b)
  |
  v
Static parameters:                    Dynamic parameters:
  c_i(Bmax) -> shape functions f_k      a_i(f) decomposition:
            -> U_k' = f_k                 f-proportional = eddy current
            -> chi_k = eta_k              f->0 intercept = aftereffect
                                          |
                                          v
                                        f_0 from a_i(f) slope
                                        Delta_U_k from chi_k (already known)
```

| Parameter | Conventional | With Fourier separation |
|-----------|-------------|------------------------|
| U_k | B-H loop fitting | Same (from c_i static component) |
| chi_k | B-H loop fitting | Same + independent verification via a_i intercept |
| f_0 | Literature value (~10^9, assumed) | **Experimentally determined from a_i(f) slope** |

No assumed parameters -- the complete model is identified from measurement.

### Relationship to Experimental Work

Congruency measurements with Prof. Hane's group (amorphous materials) have
revealed non-congruency behavior that may be explained by thermal aftereffect:
the time-dependent barrier crossing modifies the effective shape functions,
breaking ideal congruency.

For detailed research plan, see `MAGNETIC_AFTEREFFECT_RESEARCH.md`.

**Status**: Planned. Not yet implemented.

---

## References

1. Egger, Engertsberger, Schafelner: "Efficient evaluation of forward and inverse energy-based magnetic hysteresis operators", MAGCON-25-07-0171 (Compumag 2025)
2. Francois-Lavet et al.: J. Comp. Appl. Math., 246:243-250, 2013
3. Prigozhin et al.: IEEE Trans. Magn., 52:1-11, 2016
4. Egger, Engertsberger et al.: IEEE Trans. Magn., 61:7300207, 2025
5. Matsuo, T.: Play model for vector hysteresis modeling
6. MATLAB reference: `lib/PlayHysteron.m`, `lib/ShapeFunction.m`

---

**Last Updated:** 2026-03-06
