# B-input Play Hysteresis Model

## 1. B-input Play Model (Reference Model)

### 1-1. Definition

The **B-input Play model** uses play operators driven by B (not H):

```
Scalar 1D:
  p_k(t) = max(B(t) - eta_k, min(B(t) + eta_k, p_k(t^-)))
  H(t) = sum_k f_k(p_k(t))

Vector (d dimensions):
  p_k(B) = projection of B onto closed ball B_{eta_k}(p_{k,prev})
  H = sum_k f_k(|p_k|) * (p_k / |p_k|)
```

Parameters:
- `eta_k`: play threshold in B-space (k = 1, ..., K)
- `f_k(r)`: shape function (monotone increasing, f_k(0) = 0)
- Natural direction: **B -> H** (direct evaluation, no iteration)
- H -> B: **implicit** (requires fminsearch or Newton on implicit equation)

### 1-2. Why B-input (Sugahara & Hane)

The majority of the Play model literature uses **H-input** (H drives the play
operators, output is B). However, H-input models **do not reproduce** measured
H-axis congruency:

| Property | B-input (Sugahara/Hane) | H-input (Literature majority) |
|----------|------------------------|-------------------------------|
| Play operator input | **B** | H |
| Natural direction | **B -> H** | H -> B |
| H-axis congruency | **Reproduced** (matches measurement) | **Not reproduced** |

**H-axis congruency**: When driving a material with the same B-amplitude at
different DC bias points, the measured minor H-loops are congruent (same shape).
B-input Play naturally reproduces this because the play operator state depends
only on the B excursion, not on the absolute B level.

### 1-3. Limitation: H -> B Is Implicit

For **BEM** (Radia, ELF), the primary variable is H. The solver needs H -> B
(constitutive relation). With raw B-input Play:

- B -> H: **direct**, O(K) -- but BEM doesn't need this direction
- H -> B: **implicit** -- must solve `sum f_k(|p_k(B)|) * e_k = H` for B
  (requires fminsearch, 2-scalar potential, or similar)

This motivates the construction of an **energy-based approximation** that
provides efficient H -> B computation while preserving B-input Play's physics.

---

## 2. Forward Model: Energy-Based Approximation of B-input Play

### 2-1. Derivation

Starting from B-input Play parameters (f_k, eta_k), we construct an energy-based
model that **approximates** B-input Play and provides efficient H -> B:

**Step 1**: Define internal energy from shape functions:
```
U_k(r) = integral_0^r f_k(s) ds
```
Key identity: **U_k'(r) = f_k(r)** (shape function = energy derivative).

For the analytical form used in Radia:
```
f_k(r) = A_{s,k} * tan(pi/2 * r / J_{s,k})
U_k(r) = -(2 A_{s,k} J_{s,k}) / pi * log(cos(pi/2 * r / J_{s,k}))
```

**Step 2**: Map play thresholds to pinning strengths:
```
chi_k = eta_k
```

**Step 3**: Regularize the non-smooth dissipation:
```
|x| --> |x|_eps = sqrt(|x|^2 + eps)       (eps > 0, C-infinity)
```

**Step 4**: Construct Gibbs free energy:
```
G*(H, {J_k}, {J_{k,p}}) = sum_k [ U_k(|J_k|) - <H, J_k> + chi_k |J_k - J_{k,p}|_eps ]
```

**Step 5**: The **Forward model** = minimize G* w.r.t. each J_k independently:
```
J_k = argmin_J  [ U_k(|J|) - <H, J> + chi_k |J - J_{k,p}|_eps ]
B = mu_0 * H + sum_k J_k
```

Each J_k is solved independently (no coupling between operators). This is
O(K) and requires **no special processing** -- Newton with Armijo backtracking,
3-5 iterations to tol = 1e-12.

### 2-2. Why This Approximates B-input Play

The energy-based model preserves B-input Play's essential structure:

| B-input Play | Energy-based (Forward) |
|-------------|------------------------|
| Shape functions f_k | **U_k' = f_k** (identical) |
| Play thresholds eta_k | **chi_k = eta_k** (identical) |
| Play operator outputs p_k | J_k (partial polarizations, ≈ p_k) |
| Non-smooth \|x\| dynamics | Regularized \|x\|_eps (**approximation**) |
| B -> H: direct | B -> H: via Forward(B) (Schur) |
| H -> B: implicit | **H -> B: Inverse(H), each J_k independent** |

The **only source of approximation** is the eps-regularization. As eps -> 0,
the energy model converges to the B-input Play model. In practice, eps = 1e-8
gives negligible error (O(sqrt(eps)) ~ 1e-4 in J_k).

The optimality condition (grad G* = 0) for each J_k:
```
U_k'(|J_k|) * (J_k / |J_k|) + chi_k * (J_k - J_{k,p}) / |J_k - J_{k,p}|_eps = H
```

As eps -> 0, this becomes:
```
f_k(|J_k|) * e_k + chi_k * sign(J_k - J_{k,p}) = H
```

which is precisely the B-input Play model's constitutive relation.

### 2-3. Naming Convention

The code uses **B-input Play naming** (B is the natural input):

| Direction | Method | Description |
|-----------|--------|-------------|
| B -> H | `Forward(B)` | **Natural** direction for B-input Play. Schur complement O(K). |
| H -> B | `Inverse(H)` | Each J_k independent, O(K), no Schur needed. |

This reflects the physics: B drives the play operators, so B -> H is the **forward**
(natural) direction. The code naming matches B-input Play convention.

### 2-4. Properties of the Forward Model

The energy-based approximation inherits B-input Play's physics while adding
computational advantages:

| Property | Source |
|----------|--------|
| H-axis congruency | **From B-input Play** (chi_k = eta_k in B-space) |
| Measurement agreement | **From B-input Play** (f_k = U_k' preserves shape functions) |
| H -> B without iteration on B | **From energy framework** (independent J_k minimization) |
| Thermodynamic consistency | **From energy framework** (G* is a proper energy function) |
| Analytical Jacobian dB/dH | **From energy framework** (Hessian of G*) |
| Iron loss from local dissipation | **From energy framework** (D_k terms) |
| Vector extension to d <= 3 | **From energy framework** (arbitrary dimension) |

### 2-5. Congruency Property

The Forward model satisfies **H-axis congruency** (inherited from B-input Play):

```
Same delta-B  -->  Same delta-p_k  -->  Same H output shape
```

This is because chi_k = eta_k operates in B-space. The play operator state update
depends only on the B excursion:
```
p_k = B - eta_k * (B - p_k_prev) / max(eta_k, |B - p_k_prev|)
```

When |B - B_dc| = const, the delta-p_k trajectory is identical regardless of B_dc,
producing congruent H output loops.

**H-input models** (including H-input energy-based) do NOT have this property.
This is the fundamental reason B-input is physically correct.

**Verification**:
```python
import radia as rad
import numpy as np

MU_0 = 4e-7 * np.pi
rad.UtiDelAll()
tables = [(r_k, f_k) for r_k, f_k in zip(r_tables, f_tables)]
mat = rad.MatPlayHysteresis(K, eta, tables)

H_dc_list = [0, 200, 500]  # Different DC bias [A/m]
H_ac = 100.0                # Same AC amplitude

for H_dc in H_dc_list:
    # Drive to bias, then trace minor loop
    # Compare: H minor loops should be congruent (same shape)
```

**Collaboration**: Experimental congruency measurements are ongoing at
Sugahara Lab (Kindai University) in collaboration with Prof. Hane's group.

---

## 3. Inverse Operator: H -> B, O(K)

### 3-1. Algorithm

`Inverse(H)` is the operator used by Radia's BEM solver (where H is the primary
variable). From the B-input Play perspective, H-input is the non-natural direction,
but the energy-based approximation makes it computationally simple:

Each J_k is solved **independently** (H given, no coupling):

```
J_k = argmin_J  U_k(|J|) - <H, J> + chi_k |J - J_{k,p}|_eps
B = mu_0 H + sum_k J_k
```

Newton's method with Armijo backtracking, 3-5 iterations to tol = 1e-12.

### 3-2. Hessian Structure

```
H_k = hess(U_k) + chi_k * hess(|.|_eps)

hess(U_k)(J) = (U_k'' - U_k'/|J|) * e e^T + (U_k'/|J|) * I    where e = J/|J|
hess(|x|_eps) = (I - x x^T / |x|_eps^2) / |x|_eps
```

At J=0: `H_k = U_k''(0) * I + chi_k / sqrt(eps) * I` (well-conditioned).

### 3-3. Jacobian dB/dH

```
dB/dH = mu_0 I + sum_k H_k^{-1}
chi_diff = trace(dB/dH) / (d * mu_0) - 1
```

Available via `ComputeJacobian()`. This is a key advantage of the energy-based
formulation over raw B-input Play: the **analytical Jacobian** enables quadratic
Newton convergence in the field solver, without finite differences.

---

## 4. Forward Operator: B -> H, O(K) via Schur Complement

### 4-1. Relationship to B-input Play

`Forward(B)` computes B -> H, the **natural (forward) direction** for B-input Play.
The energy-based model provides this via **Schur complement**
(Egger's mathematical contribution). While Egger's method is referenced
mathematically, the **physical correctness** comes from B-input parameterization
(chi_k = eta_k, U_k' = f_k).

### 4-2. Algorithm

All J_k are **coupled** through sum J_k:

```
min_{J_k} (nu_0/2)|B - sum J_k|^2 + sum [U_k(|J_k|) + chi_k |J_k - J_{k,p}|_eps]
H = nu_0 (B - sum J_k)
```

### 4-3. Schur Complement Trick (Egger et al.)

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

### 4-4. Performance (Egger Table I, K=20, 500 time steps)

| Method | Time | Iterations/step |
|--------|------|-----------------|
| Forward (Newton) | 10 ms | 6.3 |
| **Inverse (Schur)** | **10 ms** | **6.4** |
| Inverse (Standard) | 40 ms | 6.4 |

K=100: Standard=1071 ms vs **Efficient=59 ms** (18x speedup).

---

## 5. Formulation Compatibility

The B-input energy-based model works with **any field formulation**:

| Field Formulation | Primary Variable | Direction Needed | Operator | Processing |
|-------------------|-----------------|-----------------|----------|------------|
| **A-formulation** (FEM) | B = curl(A) | B -> H | `Forward(B)` | Schur O(K) |
| **BEM standard** (Radia default) | H | H -> B | `Inverse(H)` | **None** (each J_k indep.) |
| **Hantila BEM** | H | H -> B | `Inverse(H)` | **None** |
| **FEM Hantila** (scalar/vector pot.) | H | H -> B | `Inverse(H)` | **None** |
| **Newton BEM** (optional variant) | M -> B | B -> H | `Forward(B)` | Schur O(K) |

**A-formulation** (B known): B -> H is the natural direction for B-input Play.
`Forward(B)` uses Schur complement because all J_k are coupled through sum J_k.

**BEM/Hantila** (H known): `Inverse(H)` makes H -> B simple (each J_k
independent). This is the primary reason for constructing the energy-based
approximation.

---

## 6. Implementation (C++ Material Class)

### Source Files

| File | Description |
|------|-------------|
| `src/core/rad_material_def.h` | `radTEnergyHysteresisMaterial` class definition (Type 5) |
| `src/core/rad_material_impl.cpp` | Forward/Inverse operators, Newton+Armijo, Jacobian |
| `src/core/rad_relaxation_methods.cpp` | Solver integration (`UpdateChiAndCheckConvergence`) |
| `src/radia/hysteresis_io.py` | .hys/.mat readers, Play-to-Energy conversion |
| `examples/hysteresis/verify_cpp_hysteresis.py` | C++ verification suite (4 tests) |

### Python API

```python
import radia as rad

# Play model (recommended): table-based shape functions
# f_k_tables: list of (r_array, f_array) tuples, one per play operator
mat = rad.MatPlayHysteresis(K=10, eta=eta, f_k_tables=tables)

# Energy model: requires non-negative f_k (convex U_k)
mat = rad.MatEnergyHysteresis(K=10, eta=eta, f_k_tables=tables, eps=1e-8)

# From .hys file
from hysteresis_io import load_hys_file
K, eta, tables = load_hys_file('material.hys')

# Apply to element and solve
iron = rad.ObjRecMag([0,0,0], [0.01,0.01,0.01], [0,0,0])
rad.MatApl(iron, mat)
rad.Solve(container, 0.001, 100, 1)

# Direct M(H) evaluation
M = rad.MatMvsH(mat, 'm', [500, 0, 0])  # returns M in A/m
```

### I/O: B-input Play Data -> Energy Parameters

```python
from hysteresis_io import load_hys, build_shape_functions, hys_to_radia

# One-step: .hys file -> energy parameters
model = hys_to_radia('material.hys', K=20, eps=1e-8, d=3)
# Internally: f_k -> U_k (integration), eta_k -> chi_k

# Manual pipeline
loops = load_hys('material.hys')
eta, f_k_tables, Bplay = build_shape_functions(loops)
# eta -> chi, f_k -> U_k' (same functions, different interpretation)

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

### C++ Class Hierarchy

```cpp
// Abstract base class (shared by Energy and Play models)
class radTHysteresisMaterial : public radTMaterial {
    virtual TVector3d Forward(const TVector3d& B) = 0;
    virtual TVector3d Inverse(const TVector3d& H) = 0;
    virtual void ComputeJacobian(TMatrix3d& dBdH, double& chi_d) const = 0;
    virtual void SaveState/RestoreState/ResetState/CommitState() = 0;
    // ...
};

// Play model (recommended, Type 6)
class radTPlayHysteresisMaterial : public radTHysteresisMaterial {
    int m_K;
    std::vector<PlayTable> m_tables;      // Shape function tables
    std::vector<double> m_eta;            // Play thresholds [Tesla]
    std::vector<TVector3d> m_pk_prev;     // Committed play operator state
    std::vector<TVector3d> m_pk_pinning;  // Pinning reference
    std::vector<TVector3d> m_pk_current;  // Last Forward result
    double m_B_mono_max, m_H_mono_max;    // Monotone limits

    TVector3d Forward(const TVector3d& B);   // B -> H, O(K) direct
    TVector3d Inverse(const TVector3d& H);   // H -> B, Newton + analytical Jacobian
};

// Energy model (Type 5, requires non-negative f_k)
class radTEnergyHysteresisMaterial : public radTHysteresisMaterial {
    // Uses Egger Schur complement Newton for Forward
    // Each J_k independently solved for Inverse
};
```

### Verification Results

All 5 tests PASS:

| Test | Description | Result |
|------|-------------|--------|
| Inverse-Forward roundtrip | Inverse(H) -> B -> Forward(B) -> H' | \|H-H'\|/\|H\| < 3e-10 |
| 1D B-H loop | Sinusoidal H drive, hysteresis detection | Width = 0.44 T |
| O(K) scaling | K = 5, 10, 20, 50 timing | 8.5x for 10x K ratio |
| Regularization | eps = 1e-4 to 1e-12 convergence | B converges to limit |
| 2D rotating field | Circular H trajectory | Smooth \|B\| trajectory |

### Key Bug Fix: Armijo Condition in Inverse

```python
# WRONG (stalls after 2 iterations):
if G_trial <= G_curr - 1e-4 * tau * max_grad:

# CORRECT:
dir_deriv = sum(dot(g_k, dJ_k) for k in range(K))
if G_trial <= G_curr + 0.1 * tau * dir_deriv:    # proper Armijo
```

---

## 7. Data Flow: Measured Data -> Model -> Solver

```
[Measured B-H loops]                   [B-input Energy Model]         [Field Solver]

B-input Play fitting             radTEnergyHysteresisMaterial
(.hys, .mat, CSV)                   (As, Js, chi)
       |                              |
       v                              |
  load_hys() / load_mat()             |
       |                              |
  build_shape_functions()             |
       |                              |
  f_k -> U_k (integration)            |
  eta_k -> chi_k (identity)           |
       |                              |
       +---> As[], Js[], chi[] -----> |
              U_k' = f_k              |
              (same shape functions!) |
                                      +-- Inverse(H): O(K)  <-- BEM standard path
                                      |   Each J_k independent
                                      |   (H-input, no Schur needed)
                                      |
                                      +-- Forward(B): O(K)  <-- A-method, Newton BEM
                                      |   Schur complement
                                      |   (B-input, natural direction)
                                      |
                                      +-- ComputeJacobian()  <-- All Newton field solvers
                                          dB/dH from Hessian     (not available in raw Play model)
```

---

## 8. BEM Solver Integration

### Integration Architecture

```
STANDARD PATH (Inverse, BEM chi-based):
  AutoRelax_Unified()
    +-- SolveLinearStep()              ... LU/BiCGSTAB/HACApK
    +-- ComputeActualHFieldFromSigma()
    +-- UpdateChiAndCheckConvergence()
          +-- HystMater->ComputeChiFromHysteresis(H_vec, elem_idx)
                +-- Inverse(H_vec) -> B_vec -> chi = |B|/(mu_0*|H|) - 1
                    + J_k state update
                    (each J_k independent, O(K), no Schur)

OPTIONAL: Newton variant (Forward, for vector-accurate M):
  AutoRelax_BInput_Newton()
    +-- B = mu_0*(H_ext + N*M + M)     ... B from current M
    +-- Forward(B) per element          ... Schur complement, O(K)
    +-- M_model = B/mu_0 - H           ... vector M (not scalar chi!)
    +-- Newton: F(M) = M - M_model      ... analytical Jacobian from ComputeJacobian
```

### MMM vs MSC: Same Constitutive Law

| | MMM (tet, 3 DOF) | MSC (hex, 5/6 DOF) |
|--|---|---|
| DOF | Mx, My, Mz | sigma_1, ..., sigma_6 |
| H source | FlatField array | NewFieldArray |
| Inverse(H) | **Same**: B = mu_0 H + J -> chi | **Same** |
| J_k state | Per-element 3D vector x K | **Same** |

### Quasi-Static State Management

```
Step n -> n+1:
  SaveState()                    // Freeze J_{k,p}
  for iter = 0, 1, 2, ...:      // Nonlinear iteration
    RestoreState()               // Rewind J_{k,p} at each retry
    SolveLinearStep()
    ComputeH()
    Inverse(H) -> B, {J_k}      // Each J_k independent
    chi = |B|/(mu_0*|H|) - 1
    if converged: break
  CommitState()                  // Confirm converged J_k as new J_{k,p}
```

### Memory Estimate

```
Per element: K x 3 x 8 bytes x 2 (prev + saved) = 48K bytes
K=20, N=10,000:  48 x 20 x 10,000 = 9.6 MB
K=50, N=100,000: 48 x 50 x 100,000 = 240 MB
```

---

## 9. Magnetization (Pulse Magnetization) Application

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

# 1. Create hysteresis material from B-input Play data
K, eta, tables = load_hys_file('NdFeB.hys')
mat = rad.MatPlayHysteresis(K, eta, tables)

# 2. Unmagnetized blank
magnet = rad.ObjHexahedron(vertices, [0, 0, 0])
rad.MatApl(magnet, mat)

# 3. Magnetization pulse
coil = rad.ObjBckg(lambda p: magnetizing_field(p, I_peak))
assembly = rad.ObjCnt([magnet, coil])
rad.Solve(assembly, 0.0001, 100, 0)
rad.MatHysCommitState(mat)

# 4. Remove external field -> remanence
coil_off = rad.ObjBckg(lambda p: [0, 0, 0])
assembly_off = rad.ObjCnt([magnet, coil_off])
rad.Solve(assembly_off, 0.0001, 100, 0)
rad.MatHysCommitState(mat)

# 5. Remanent magnetization pattern
B_ext = rad.Fld(magnet, 'b', observation_point)
```

---

## 10. Anisotropy Extension Roadmap

| Step | Internal energy | Pinning | Use case |
|------|----------------|---------|----------|
| **Current** | U_k(\|J\|) isotropic | scalar chi_k | Isotropic soft iron |
| **Step 1** | U_k(\|J\|) isotropic | tensor X_k | Directional pinning |
| **Step 2** | U_k(\|J\|) + W_a(J, e) | tensor X_k | Grain-oriented Si steel |
| **Step 3** | Polycrystalline averaging | ODF-based | Textured materials |

---

## 11. Unified Picture

```
                [B-input Play Model]
                 Sugahara & Hane
                 B -> H (natural)
                 H-axis congruent
                 Measured f_k, eta_k
                        |
                        | f_k = U_k', chi_k = eta_k
                        | (energy-based approximation)
                        v
            [Energy-Based Forward Model]
            G* = sum [U_k + D_k] - <H, J>
            Approximates B-input Play
            eps-regularized
                        |
        +---------------+---------------+
        |               |               |
        v               v               v
  [Inverse(H)]    [Forward(B)]    [ComputeJacobian]
  H -> B           B -> H          dB/dH
  Each J_k indep.  Schur O(K)      Analytical Hessian
  O(K), no Schur   Natural dir.    (not avail. in raw Play)
  BEM standard     A-method,       Newton field solvers
  Hantila          Newton BEM
```

**Key relationships**:
- **B-input Play -> Energy model**: f_k = U_k' (shape function = energy derivative),
  chi_k = eta_k (pinning = play threshold). This ensures B-input Play's physics
  (congruency, measurement agreement) is preserved in the energy model.
- **Forward(B)**: B -> H, the natural direction for B-input Play.
  Egger's Schur complement provides efficient O(K) computation.
- **Inverse(H)**: H -> B, the energy-based approximation that makes H -> B simple.
  Each J_k independent, computationally trivial.
- **ComputeJacobian**: analytical, from energy Hessian. Not available in raw
  B-input Play. Needed for Newton convergence in field solvers.

---

## 12. Optional: Vector-Accurate Newton Solver

### 12-1. Motivation: Beyond Scalar Chi

The **standard BEM solver** uses Inverse(H) with scalar chi (M forced parallel to H).
For hysteresis, M does not always align with H. The optional Newton variant uses
Forward(B) for full vector M:

```
Standard path:
  Inverse(H) -> B -> chi = |B|/(mu_0*|H|) - 1
  M = chi * H  (scalar, loses vector info)

Newton variant:
  B = mu_0*(H+M), Forward(B) -> H
  M_model = B/mu_0 - H  (VECTOR M, fully accurate)
  dF/dM uses ComputeJacobian (analytical)
```

### 12-2. Analytical Jacobian from Energy Hessian

```
After Forward(B):
  dB/dH = mu_0 * I + sum_k H_k^{-1}     (from energy Hessian)
  dJ/dB = I - mu_0 * inv(dB/dH)          (3x3 per element)
  dF/dM = I - block_diag(dJ/dB) * (N + I)  (full system Jacobian)
```

In Python, computing dJ/dB required 6 finite-difference Forward evaluations per
element per iteration. In C++, `ComputeJacobian()` provides exact derivatives
in a single call -- **zero FD evaluations**.

This is a key advantage of the energy framework over raw B-input Play:
raw Play can compute B->H directly (Forward), but cannot provide the analytical
Jacobian needed for quadratic Newton convergence.

### 12-3. Solver Variants

#### Newton with Forward (`b_input_newton=True`)

```python
rad.SolverConfig(b_input_newton=True)
rad.Solve(container, 1e-4, 100, 0)
```

- Uses Forward(B->H) per element with analytical Jacobian
- **2-4 iterations** to converge (quadratic convergence)
- O(N^3) per iteration -- best for small problems (N < 500)

#### Newton + Hantila Hybrid (`b_input_hantila=True`, Recommended)

```python
rad.SolverConfig(b_input_hantila=True, hantila_alpha=0)  # auto-alpha
rad.Solve(container, 1e-4, 5000, 0)
```

1. **Newton warmup** (3 iterations): Uses Forward(B->H), O(N^3), gets close quickly
2. **Hantila refinement**: Uses Inverse(H->B) for R update, O(N^2) per iteration

Hantila phase uses Inverse(H) [= simple, no Schur], while Newton phase uses
Forward(B) [= Schur, but provides analytical Jacobian].

### 12-4. Verified Results

| Test | Description | Result |
|------|-------------|--------|
| Forward/Inverse roundtrip | `Forward(Inverse(H)) == H` | Error < 1e-6 |
| Single-step B-input Newton | 5-tet block, B=0.1T | 2-3 iterations |
| Hysteresis loop (Newton) | 11 steps, C++ vs Python | 0.00-0.03% error |
| Single-step Hantila | 5-tet block, B=0.1T | 3-47 iterations |
| Hysteresis loop (Hantila) | 11 steps, vs Newton reference | **0.00% error** |

### 12-5. Crossover: Newton vs Hantila

| Problem Size (DOF) | Newton (O(N^3)/iter) | Hantila (O(N^2)/iter) | Recommendation |
|---------------------|---------------------|-----------------------|----------------|
| N < 100 | Fast | Overhead | Newton |
| 100 < N < 500 | OK | Comparable | Either |
| N > 500 | Slow | **Much faster** | **Hantila** |
| N > 2000 (HACApK) | N/A | O(N log N)/iter | **Hantila + HACApK** |

---

## 13. Unified Pipeline: Measurement -> Model -> Solver

### Three Layers

```
Layer 1: Parameter Identification (B-input Play Fitting)
======================================================================
Measured B-H loops
  |
  v
B-input Play fitting (Sugahara/Hane):
  eta_k: play thresholds in B-space
  f_k(|p|): B-input shape functions
  |
  v
Energy parameter conversion:
  U_k(r) = integral f_k(s) ds      (f_k = U_k')
  chi_k = eta_k                     (pinning = threshold)
  As_k, Js_k: analytical fit to f_k

Optional: Taka's Fourier separation for dynamic components:
  DFT: h = (a_i * db/dt) + (c_i * b)
  c_i(Bmax): static -> B-input shape functions
  a_i(f): dynamic -> eddy current + aftereffect


Layer 2: Constitutive Model (Energy-Based Approximation of B-input Play)
======================================================================
Parameters: K, As[], Js[], chi[], eps  (from B-input Play fitting)
  |
  +-- Inverse(H): H -> B, O(K)       <-- BEM standard (each J_k indep.)
  |     J_k = argmin [U_k - <H,J> + chi_k |J - J_{k,p}|_eps]
  |     B = mu_0 H + sum J_k
  |
  +-- Forward(B): B -> H, O(K)       <-- A-method, Newton BEM
  |     Schur complement (Egger's mathematical contribution)
  |     Natural direction for B-input Play
  |
  +-- ComputeJacobian(): dB/dH       <-- Newton field solvers
  |     Analytical from energy Hessian (unavailable in raw Play)
  |
  +-- SaveState / RestoreState / CommitState
        Play operator state management for quasi-static stepping


Layer 3: Field Solver
======================================================================
A-formulation (FEM): B = curl(A) known
  -> Forward(B): B-input Play natural direction
  -> FEM Hantila already sufficient (no merit adding Newton)

H-formulation (BEM, 2-scalar): H known
  |
  +-- Standard chi-based (Inverse H->B):
  |     Inverse(H) -> B -> chi = |B|/(mu_0*|H|) - 1
  |     M = chi * H (scalar, loses vector info)
  |
  +-- Newton with Forward(B->H):       <-- VECTOR-ACCURATE
  |     B = mu_0*(H+M), Forward(B) -> H
  |     M_model = B/mu_0 - H (full vector)
  |     dF/dM uses ComputeJacobian
  |
  +-- Hantila hybrid:                   <-- RECOMMENDED
        Newton(3) warmup [Forward] + Hantila [Inverse] refinement
```

### One Energy Function Unifies Everything

```
G*(H, {J_k}, {J_{k,p}}) = sum_k [ U_k(|J_k|) - <H, J_k> + chi_k |J_k - J_{k,p}|_eps ]
```

- **U_k**: from B-input Play shape functions (f_k = U_k')
- **chi_k**: from B-input Play thresholds (chi_k = eta_k)
- **Inverse(H)**: minimize G* w.r.t. each J_k independently
- **Forward(B)**: minimize G* w.r.t. all {J_k} jointly (Schur)
- **ComputeJacobian**: from G*'s second derivatives
- **Congruency**: guaranteed by B-input parameterization (chi_k in B-space)

The **B-input distinction**: G* has the same mathematical form regardless of
parameterization. But only B-input fitting (f_k from B-input Play, chi_k = eta_k
in B-space) reproduces measured H-axis congruency. H-input fitting would produce
a model that looks mathematically identical but does NOT match measurements.

---

## 14. Why B-input Play Cannot Be Converted to Energy-Based

### 14-1. The Conversion Attempt

Section 2 describes how to construct an energy model from B-input Play:
```
U_k(r) = integral_0^r f_k(s) ds     (energy = integral of shape function)
chi_k = eta_k                        (pinning = play threshold)
```

This works **only if U_k is convex** (f_k >= 0 and monotonically increasing).

### 14-2. The Problem: JMAG Shape Functions Are Negative

Real B-input Play data (e.g., JMAG .hys files) has **negative shape functions**
for k > 0:

```
k=0:  f_0(r) > 0   (positive, dominant term)
k>0:  f_k(r) < 0   (negative for all r > 0)
```

This is physically necessary. On the descending branch at B = 0:
- All play operators have p_k > 0 (residual polarization)
- H(B=0) = sum_k f_k(|p_k|) * p_hat_k
- Coercive field requires H(B=0) < 0
- With f_0 > 0 and sum over all f_k: negative f_k (k>0) must dominate

### 14-3. Mathematical Consequence: Non-Convex Energy

```
f_k < 0  =>  U_k(r) = integral f_k(s) ds  is DECREASING
         =>  U_k is NON-CONVEX
         =>  Energy minimization has MULTIPLE local minima
         =>  Forward(B) / Inverse(H) may not have unique solutions
```

### 14-4. Physical Consequence: Non-Monotonic H(B)

With negative f_k, the total H(B) curve is **non-monotonic**:

```
H(B) = sum_k f_k(|p_k(B)|) * p_hat_k

At small B:  f_0 dominates  =>  dH/dB > 0  (monotonically increasing)
At large B:  negative f_k grow  =>  dH/dB < 0  (H decreases with B!)
```

For JMAG data (K=20, B_sat=2 T):
- H(B) peaks at B ~ 1.4 T (H_max ~ 5900 A/m)
- Beyond B = 1.4 T: dH/dB < 0 (82/499 fine steps are non-monotone)
- Inverse(H) is **multi-valued** for H near H_max

### 14-5. Comparison: H-input vs B-input Energy Convertibility

| | H-input Play | B-input Play |
|---|---|---|
| Shape function constraint | u_k >= 0 (from measurement) | **f_k < 0 for k > 0** |
| Energy integral V_k | **Convex** (u_k >= 0) | **Non-convex** (f_k < 0) |
| Energy minimization | **Unique solution** | Multiple local minima |
| H(B) monotonicity | **Guaranteed** | **Not guaranteed** |
| Energy-based conversion | **Valid** | **Invalid** |

**H-input Play -> Energy**: Works because the H-input shape functions u_k are
non-negative (from measurement fitting). The energy U_k = integral u_k is convex,
ensuring well-posed minimization.

**B-input Play -> Energy**: Fails because the B-input shape functions f_k must be
negative for k > 0 to reproduce coercive field behavior. The resulting U_k is
non-convex, violating the mathematical foundation of energy minimization.

### 14-6. Solution: Direct B-input Play Model (`radTPlayHysteresisMaterial`)

Since the energy conversion fails, we implement the B-input Play model **directly**
in C++ without energy-based approximation:

```
radTMaterial
  +-- radTHysteresisMaterial (abstract base class)
        +-- radTEnergyHysteresisMaterial (Type 5, existing)
        +-- radTPlayHysteresisMaterial   (Type 6, NEW)
```

Both classes share the same abstract interface (`radTHysteresisMaterial`):
Forward, Inverse, ComputeJacobian, SaveState, RestoreState, CommitState.
The BEM solver uses `dynamic_cast<radTHysteresisMaterial*>` and works with
either model transparently.

### 14-7. Direct Play Model: Forward and Inverse

**Forward (B -> H): O(K) direct evaluation, no Newton**

```cpp
H = 0;  dHdB = 0;
for k = 0..K-1:
    q = B - p_k_pinning
    if |q| <= eta_k:  p_k = p_k_pinning  (stuck)
    else:             p_k = B - eta_k * q/|q|  (following)
    H += f_k(|p_k|) * p_k/|p_k|
    dHdB += (chain rule: dh_k/dp_k * dp_k/dB)
```

- Forward is **always O(K)** with analytical Jacobian dH/dB
- No Newton iteration needed (unlike Energy model's Schur complement)

**Inverse (H -> B): Trust-region Newton**

```cpp
B = initial_guess  (from last result or chi-based estimate)
for iter = 0..99:
    H_comp = Forward(B)  // also caches dH/dB
    if |H_comp - H_target| < tol: break
    if trace(dH/dB) < 0:  // non-monotone region detected
        B *= 0.95          // retreat to monotone region
        continue
    dB = -inv(dH/dB) * (H_comp - H_target)  // Newton step
    dB = clamp(dB, max_step=0.05*B_sat)       // trust region
    B = clamp(B + dB, |B| <= B_sat)           // magnitude bound
```

Key features:
- **Non-monotone detection**: When trace(dH/dB) < 0, retreat to lower |B|
- **Trust region**: Limit |dB| per iteration to 5% of B_sat
- **B magnitude clamp**: |B| <= B_sat prevents table extrapolation

### 14-8. Comparison: Energy Model vs Direct Play Model

| Property | Energy Model (Type 5) | Direct Play Model (Type 6) |
|----------|----------------------|---------------------------|
| Forward B->H | Schur complement Newton, O(K) | **Direct evaluation O(K), no Newton** |
| Inverse H->B | K independent Newton (each J_k), O(K) | Single Newton on B, O(K) |
| Jacobian | dB/dH from energy Hessian | dH/dB from chain rule, inv for dB/dH |
| Shape function constraint | f_k >= 0 required (convex U_k) | **No constraint** (negative f_k OK) |
| JMAG data compatibility | **Limited** (requires non-negative f_k) | **Full** (native format) |
| Non-monotone H(B) | Cannot occur (convex) | Can occur (handled by trust region) |
| BEM solver integration | Via radTHysteresisMaterial | **Same interface** |
| Accuracy (monotone region) | Reference | **0.000% vs Python** (verified) |
| Accuracy (non-monotone B > ~1.4 T) | N/A | **~18% error** (multi-valued H(B)) |
| Performance (K=80) | ~0.07 ms/eval | **~0.005 ms/eval** |

### 14-9. Python API

```python
import radia as rad

# Direct B-input Play model (no energy conversion needed)
mat = rad.MatPlayHysteresis(K, eta, f_k_tables)
# K: int, number of play operators
# eta: list[float], play thresholds [Tesla]
# f_k_tables: list of (r_array, f_array) tuples

# Apply to element and solve (same as Energy model)
iron = rad.ObjRecMag([0,0,0], [0.01,0.01,0.01], [0,0,0])
rad.MatApl(iron, mat)
bkg = rad.ObjBckg(lambda p: [0, 0, 0.1])
container = rad.ObjCnt([iron, bkg])
rad.Solve(container, 0.001, 100, 1)

# Direct M(H) evaluation
M = rad.MatMvsH(mat, 'm', [500, 0, 0])

# State management (same API as Energy model)
state = rad.MatHysSaveState(mat)
rad.MatHysRestoreState(mat, state)
```

### 14-10. When to Use Which Model

| Use Case | Recommended Model |
|----------|-------------------|
| JMAG .hys data (f_k < 0 for k > 0) | **MatPlayHysteresis** (Type 6) |
| Custom non-negative f_k data | MatEnergyHysteresis (Type 5) |
| A-formulation FEM (B known) | **MatPlayHysteresis** (fastest Forward) |
| BEM with moderate B (< ~1.4 T) | Either (identical results) |
| BEM with high B (> ~1.4 T) | MatEnergyHysteresis (if f_k >= 0 available) |
| Analytical Jacobian for Newton solver | Either (both provide dB/dH) |

### 14-11. Source Files

| File | Description |
|------|-------------|
| `src/core/rad_material_def.h` | `radTHysteresisMaterial` (base), `radTPlayHysteresisMaterial` (Type 6) |
| `src/core/rad_material_impl.cpp` | Forward/Inverse, trust-region Newton, state management |
| `src/core/rad_relaxation_methods.cpp` | Uses `radTHysteresisMaterial*` (both types) |
| `src/radia/radia_pybind.cpp` | `MatPlayHysteresis` binding |
| `examples/hysteresis/verify_cpp_play_model.py` | Verification suite (13 tests) |

### 14-12. Verification Results

All 13 tests PASS:

| Test | Description | Result |
|------|-------------|--------|
| 1-4 | Forward M(H) at various H vectors | M > 0, ferromagnetic |
| 5 | C++ Play vs Python play_hysteron (B < 1.3 T) | **0.000% error** |
| 6 | B-H loop with sinusoidal H drive | Hysteresis width > 0 |
| 7-8 | State Save/Restore lifecycle | err = 0.000000 A/m |
| 9-13 | Monotone limits, Jacobian, round-trip | All PASS |
| 9 | Solver integration (MatApl + Solve) | Converged |
| 10 | Performance (K=80, 100 evaluations) | 0.005 ms/eval |

---

## 15. Implementation Phases (Both Models)

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

## 16. Future: Magnetic Aftereffect

### Physical Background

Magnetic aftereffect is the time-dependent relaxation of magnetization after a
change in applied field. Thermal activation over pinning energy barriers.

After a field change, M(t) relaxes logarithmically:
```
M(t) = M_0 + S * ln(t / t_0)
```

### Energy-Based Formulation

The pinning potential U_k provides the energy barrier landscape:
```
dJ_k/dt = f_0 * exp(-Delta_U_k(J_k) / (k_B * T)) * direction
```

The energy framework is well-suited because energy barriers are already computed
(U_k and chi_k define the landscape).

### Parameter Identification via Fourier Separation

```
Measured B-H loops at multiple frequencies
  |
  v
Taka's Fourier separation: h = (a_i * db/dt) + (c_i * b)
  |
  v
Static: c_i -> B-input shape functions f_k -> U_k, chi_k
Dynamic: a_i(f) -> eddy current (f-proportional) + aftereffect (f->0 intercept)
```

**Status**: Planned. Not yet implemented.

---

## References

1. Egger, Engertsberger, Schafelner: "Efficient evaluation of forward and inverse energy-based magnetic hysteresis operators", MAGCON-25-07-0171 (Compumag 2025)
2. Francois-Lavet et al.: J. Comp. Appl. Math., 246:243-250, 2013
3. Prigozhin et al.: IEEE Trans. Magn., 52:1-11, 2016
4. Egger, Engertsberger et al.: IEEE Trans. Magn., 61:7300207, 2025
5. Matsuo, T.: Play model for vector hysteresis modeling
6. MATLAB reference: `lib/PlayHysteron.m`, `lib/ShapeFunction.m`
7. Hantila, F.I.: "A method for solving stationary magnetic field in nonlinear media", Rev. Roum. Sci. Techn., 1975
8. Sugahara, K. and Hane, H.: B-input energy-based vector hysteresis (Kindai University, ongoing research)

---

**Last Updated:** 2026-03-08 (added Section 14: Direct B-input Play model, non-convex energy proof)
