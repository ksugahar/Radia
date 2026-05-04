# Cauer Ladder Network (CLN) — Tanimoto-Kameari Method

Comprehensive theory + implementation guide for Cauer ladder network
construction from FEM eddy current analysis. Source material: Tanimoto
master's thesis (W:/00_CAE/NGSolve/谷本/) + Kameari iterative
orthogonalization. Connection to Kelvin transformation (recent
2026-05-04 findings) is documented in §6.

---

## 1. Overview

Cauer Ladder Network (CLN) is a model-order reduction technique that
expresses the eddy-current admittance of a conductor as a series of
{R_n, L_n} circuit elements, providing a compact equivalent circuit
suitable for time-domain or frequency-domain coupling with external
circuits.

### 1.1 Method Genealogy

- **Foster decomposition**: any rational admittance can be written as
  Σ a_n / (1 + s τ_n) (partial fractions on the imaginary axis).
- **Cauer-II ladder**: continued-fraction realization at s = 0; gives
  ladder of (R, L) pairs forming a physical circuit topology.
- **Kameari iteration** (2008): obtain {R_n, L_n} by FEM-based
  iterative orthogonalization in σ-weighted L² inner product on the
  conductor. Converges to the Foster modes order-by-order.
- **Tanimoto formulation** (修論, 2024): three FE-space variants
  (A-T, T-Ω, A-Φ) for 3D HCurl, all yielding the same ladder.

### 1.2 Use Cases

- **Eddy current Cauer model** for transformer / inductor / motor
  windings → port impedance vs frequency.
- **Time-domain SPICE / circuit simulation** of induction heating,
  switching transients.
- **Sensitivity analysis**: ladder parameters depend smoothly on
  geometry / material, suitable for optimization.

---

## 2. Mathematical Foundation

### 2.1 Iterative Orthogonalization (Kameari)

Given a conductor Ω_cond with conductivity σ and an impressed
electric source E_imp on the boundary (or impressed J inside), the
Kameari iteration generates the Cauer ladder via:

```
J_0 = σ · A_imp                       (initial impressed current)
A_pot_0 = 0
for n in 0, 1, 2, ...:
    Solve curl-curl for A_n with source J_n in Ω_cond:
        a(A_n, v) = (J_n, v)_cond            # bilinear form
    R_n = 1 / ⟨J_n, J_n / σ⟩_cond            # σ-weighted norm
    A_pot_{n+1} = A_pot_n + R_n · A_n
    L_n = R_n · ⟨J_n, A_pot_{n+1}⟩_cond     (Tanimoto accumulated form)
    J_{n+1} = J_n - σ · A_pot_{n+1} / L_n   # Schmidt orthogonalization
```

The σ-weighted inner product `⟨·,·⟩_cond` on the conductor is the
natural Hilbert structure for eddy current dissipation (resistive
power).

### 2.2 Cauer-II Ladder Synthesis

After N stages, the admittance is reconstructed via Cauer-II
continued fraction at s = 0:

```
Y(s) = 1 / (s·L_0 + 1 / (R_0 + 1 / (s·L_1 + 1 / (R_1 + ...))))
```

Equivalently, in Foster partial-fraction form:

```
Y_Foster(s) = Σ_n  a_n / (1 + s τ_n)
```

where τ_n = L_n · R_n (NOT L_n / R_n — the Cauer ladder gives time
constants from the products) and a_n are residues at the poles.

### 2.3 Foster Pole Identification

The dominant Foster pole τ_lead = MAX_n τ_n corresponds to the
slowest decay mode. For typical geometries this is the lowest
TE/TM eigenmode of the conductor.

**Important**: τ_lead is NOT necessarily τ_0 (the first Cauer stage).
Cauer stages probe different "directions" in the eigenspace; the
dominant Foster mode emerges by superposition.

### 2.4 Validation: Cylinder Analytical

For a long cylindrical conductor of radius r, height h:

```
R_theory[n] = (2n + 1) / (π · r² · σ · h)
L_theory[n] = μ_0 · (2n² + 2n + 1) / (8 · (n + 1) · π · h)
```

(Bessel-mode decomposition; Tanimoto 修論 §3.) These provide
machine-precision benchmarks for code verification.

For a closed PEC cuboid a × b × c (uniform B_z applied):

```
τ_lead = μ_0 · σ · a² · b² / (π² · (a² + b²))      (TE_z(1,1,0) mode)
```

For 5 × 2 × 1 mm Cu: τ_lead = 25.46 μs (analytical) vs Kameari max
τ_n over 12 stages = 25.33 μs (0.5% match).

---

## 3. Three Formulations (Equal Physics, Different FE Spaces)

### 3.1 A-T Formulation (primary 3D)

Two HCurl spaces, decoupled:

```python
fesA = HCurl(mesh, order, nograds=True, dirichlet="conductorBND|in|out")
fesT = HCurl(mesh, order, nograds=True)
```

**T solve** (curl-curl in conductor with boundary E source):
```python
a_T(T, W) = ∫_cond (1/σ) curl T · curl W dx
f_T(W)   = -∫_∂cond (E_s × n) · W ds       # E source on boundary
J = curl T                                    # current density
```

**A solve** (curl-curl with R·J in conductor):
```python
a_A(A, N) = ∫_full (1/μ) curl A · curl N dx
f_A(N)   = ∫_cond R_n · J · N dx
```

**Distinctive**: T resides in the conductor only (since 1/σ = ∞ in
air would make the bilinear form singular elsewhere); A is global.

### 3.2 T-Ω Formulation (hybrid)

```python
fesT = HCurl(mesh, order, nograds=True)
fesΩ = H1(mesh, order, definedon=mesh.Materials("conductor"))
```

Scalar Ω confined to conductor for current flux constraint:
∫_∂cond H_t = 0 across non-source boundaries (curl-free in air, T = -∇Ω).

Coupled solve via Lagrange multiplier or augmented bilinear form.

**Advantage**: smaller global DOF count (Ω is conductor-only).
**Disadvantage**: implementation complexity (mixed space coupling).

### 3.3 A-Φ Formulation (volume current)

```python
fesA = HCurl(mesh, order, dirichlet="...")
fesΦ = H1(mesh, order, definedon="conductor", dirichlet="in|out")
```

Body current J = σ ∇Φ (implicitly via weak form). Lukas variant
explores order = 3 with mixed-space `V = HCurl × H1`.

**Use case**: when scalar Φ is more natural (e.g., end-effect studies
in machines).

### 3.4 Equivalence

All three formulations produce the same ladder {R_n, L_n} (within
discretization error). Choice depends on:
- DOF count (T-Ω smallest, A-T largest)
- Coupling complexity (A-T simplest, T-Ω most complex)
- Boundary condition convenience (A-Φ for end-effect, A-T for general)

---

## 4. 2D CLN (Reference / Validation)

For 2D problems (plane geometry, B_z applied perpendicular to plane),
the formulation reduces to a single H1 unknown.

### 4.1 Strong / Weak Form

Magnetic vector potential `A = A_z ẑ` only (out of plane). The 2D
eddy current equation:

```
σ ∂A_z/∂t - ∇·(ν ∇A_z) = -σ E_z       (E_z impressed)
```

Weak form (steady-state Kameari):
```
∫_Ω ν ∇A · ∇v dx = ∫_cond J · v dx
```

### 4.2 Implementation (Tanimoto 2D notebook)

```python
fes = H1(mesh, order=order, dirichlet="outer")
J = sigma * E         # E impressed (constant)

for nStage in range(N):
    R = 1 / Integrate(J*J/sigma * dx, mesh)
    a = BilinearForm(fes)
    a += 1/mu * grad(u) * grad(v) * dx
    f = LinearForm(fes)
    f += v * J * dx
    # solve, update Apot, L, J
    L = R * Integrate(J * Apot * dx, mesh)
    J = J - sigma * Apot / L
```

### 4.3 Cylinder Analytical (2D)

For a cylindrical conductor radius r:
```
R_theory[n] = (2n + 1) / (π · r² · σ)              (per unit length)
L_theory[n] = μ_0 · (2n² + 2n + 1) / (8 · (n + 1) · π)
```

The 2D form is the **per-unit-length** version of the 3D cylinder
analytics (cf. §2.4). Used as a sanity check before scaling to 3D.

### 4.4 When to Use 2D vs 3D

| Problem | 2D CLN | 3D CLN |
|---|---|---|
| Long uniform conductor (translation symmetry) | ✅ ideal | overkill |
| End-effect study | ❌ misses | ✅ required |
| Rotational machine cross-section | ✅ for one slice | ✅ for full machine |
| Small benchmark / validation | ✅ fast | slow |

---

## 5. Constraint and Gauge Variants

In HCurl, the curl-curl operator has a non-trivial nullspace
(gradients), requiring gauge fixing to ensure solvability.

### 5.1 nograds=True (Tanimoto canonical)

NGSolve `HCurl(mesh, order, nograds=True)` removes basis functions
that are pure gradients — eliminates the nullspace by construction.
The Kameari iteration is well-posed with this choice.

### 5.2 A-Penalty

Add a small mass-matrix term:
```
a_penalty(A, N) = ∫ (1/μ) curl A · curl N dx + ε ∫ (1/μ) A · N dx
```
with ε ~ 1e-6. Stabilizes the curl-curl operator without explicit
gauge correction. Loses ∇·A = 0 exactly but gives equivalent ladder
within tolerance.

**Pros**: simpler than explicit gauge.
**Cons**: ε tuning; loses divergence-free property.

### 5.3 A-Gauge (Coulomb explicit)

Two-step:
1. Solve A via HCurl (nograds=True) — gives A modulo gradient
2. Project to Coulomb gauge: solve H1 problem ∇²φ = ∇·A,
   then update A := A - ∇φ

Maintains ∇·A = 0 to machine precision but adds an H1 solve per
stage (overhead).

### 5.4 Comparison

| Variant | Notebook | Gauge precision | Overhead |
|---|---|---|---|
| nograds=True | CLN_AT (修論) | exact | none |
| A-Penalty | 20231211_A_(Penalty) | ε-tolerance | mass term |
| A-Gauge | 20231221_A_gauge_CLN | machine | H1 solve / stage |

Production recommends **nograds=True** (cleanest, no tuning).

---

## 6. Connection to Kelvin Transformation

CLN treats a CONDUCTOR in finite domain (typically truncated air box).
Kelvin transformation extends this to UNBOUNDED exterior (open
boundary). Combining the two enables CLN for problems where the
eddy current decays into an unbounded region (e.g., isolated coil
in vacuum).

### 6.1 Status (2026-05-04)

**Recent finding**: CLN + Kelvin had been an open research problem.
Combining standard reduced-A formulation with Kelvin pullback gave
+43% inductance error on a benchmark torus problem until the root
cause was identified.

**Root cause**: the popular `(ν - ν₀) curl(A_s) · curl(v)` reduced-A
weak form is **invalid when A_s in the Kelvin region is the proper
1-form pullback**. The simplification requires A_s to satisfy
ν₀ Maxwell globally; the pullback satisfies ν' Maxwell instead
(metric-dependent).

**Correct form**: drop the (ν - ν₀) middle step and use directly:
```
a(A_r, v) = -∫_kext ν' · (∇×A_s_pullback) · (∇×v) dV
```
which gives +6% (matching the J-source baseline).

See [`docs/kelvin/KELVIN_TRANSFORMATION.md`](../kelvin/KELVIN_TRANSFORMATION.md)
§7.5 for the full derivation.

### 6.2 Implications for CLN + Kelvin

The Kameari iteration for CLN involves solving curl-curl problems
each stage with a specific impressed source. If Kelvin is used:

1. **A-T formulation + Kelvin**: A solve has the same (ν - ν₀)
   pitfall — must use direct form on kext.
2. **T-Ω + Kelvin**: T resides in conductor only (inner) — no Kelvin
   pullback needed for T. Ω in air + kext — uses scalar pullback (0-form,
   trivial). **T-Ω is structurally cleaner for CLN + Kelvin** than A-T.
3. **A-Φ + Kelvin**: similar to A-T, has (ν - ν₀) pitfall.

### 6.3 Specific Benchmark: Cuboid 5×2×1 in Vacuum

A 5 × 2 × 1 mm Cu cuboid with uniform B_z applied at infinity is a
canonical test for CLN + Kelvin. Without Kelvin (with absorbing BC at
finite distance), τ_lead converges to 25.46 μs (analytical TE_z(1,1,0)).

Recent attempt (2026-05-04, NGSolve A-formulation reduced-A + Kelvin):
- Initial naive (ν-ν₀) form: τ_0 = 1.95 × 10⁹ μs (broken)
- Various A_s pullback variants: stage-0 sign flip, +43% errors
- Diagnosis: same (ν - ν₀) pitfall as above + uniform B_z producing
  unbounded A_s at infinity (incompatible with Kelvin pullback)

**Status**: H-formulation or T-Ω with reduced-Ω = -H_0 z + Ω_r (Ω_r
decays at infinity, no unbounded A) is the recommended path for
applied uniform field. PEEC source (decaying field) + reduced A-form
+ Kelvin direct-form works (proven on torus benchmark).

### 6.4 Recommended Path Forward

For CLN of an isolated conductor in **applied uniform field**:
1. Use H-formulation (or T-Ω) with reduced potential approach.
2. Apply Convention B background field transformation:
   `H_s' = -(ρ'/R)² H_s` (vanishes at offset, no singularity).
3. Use direct form (avoid (ν - ν₀) middle step).
4. Iterate Kameari for ladder construction.

For CLN of conductor with **localized PEEC coil source** (filament
bundle with finite extent):
1. Use A-formulation reduced (Convention A pullback for A_s in kext).
2. Direct form: `-ν' (∇×A_s) (∇×v) dx("kelvin")` (NOT (ν-ν₀)).
3. Iterate Kameari.

---

## 7. Solver Variants and Production Code

### 7.1 Solver Choices

| Solver | Notebook | Distinctive |
|---|---|---|
| Direct sparse (Pardiso/MKL) | A_direct, CASE_*_direct | Reference baseline, robust |
| NGSolve CG + local pre | A_CG, CASE_*_CG | Standard Krylov |
| SparseSolvPy ICCG | CLN_AT (修論) | JP-MARs research backend |
| accICCG | CASE_*_accICCG, A_ICCG_最新版 | Acceleration param tuning |

**Production**: A + ICCG (`20240917_A_ICCG_最新版.ipynb`), includes
inline gauge correction.

### 7.2 bonus_intorder for Higher-Order Elements

For order ≥ 3 HCurl with quadratic A_imp source (e.g., (B₀/2)(-y, x, 0)),
default integration order is too low and Schmidt orthogonality drift
appears prematurely. Use `bonus_intorder=8` in all `dx()` calls to
keep drift at machine precision through 11+ stages.

```python
a += (1/mu) * curl(A) * curl(N) * dx(bonus_intorder=8)
```

### 7.3 Schmidt Orthogonality Drift Diagnostic

```
drift_n = max_{m<n} |⟨J_n, J_m / σ⟩_cond| / ⟨J_n, J_n / σ⟩_cond
```

For 3D HCurl order=3 with bonus_intorder=8 (closed PEC cuboid):
| Stage range | Drift | Status |
|---|---|---|
| N ≤ 11 | ≤ 1e-12 | machine precision |
| N = 12 | 3e-12 | exponential growth onset (~×6/stage) |
| N = 25 | 5.5% | 1% breakdown threshold |
| N ≥ 26 | corrupted | regenerate basis |

---

## 8. Benchmark Cases

| Case | Geometry | Applied | Reference | Validation |
|---|---|---|---|---|
| 2D cylinder | r = 0.01 m | uniform E_z | Bessel analytics | machine precision |
| 3D cylinder | r × h = 0.01 m × 0.01 m | uniform E_z | Bessel × geometric | <1% match |
| 3D cuboid 5×2×1 (closed PEC) | Cu cuboid in PEC box | uniform B_z | TE_z(1,1,0) τ = μσa²b²/(π²(a²+b²)) | 0.5% match |
| 3D cuboid 5×2×1 (vacuum) | Cu cuboid in air | uniform B_z at infinity | (Kelvin or absorbing) | open problem, see §6.3 |
| TEAM-28 | scaled benchmark | (per spec) | TEAM database | validated |

---

## 9. Implementation File Index

### Knowledge Base (mcp-server)

```
packages/radia-mcp/src/radia_mcp/radia_ngsolve/knowledge/
├── cln_3d.py                     # high-level overview
└── cln_notebooks/
    ├── __init__.py               # registry
    ├── CLN_AT.py                 # Tanimoto 修論 A-T (primary 3D)
    ├── CLN_T_Omega.py            # T-Ω formulation
    ├── CLN_APhi.py               # A-Φ formulation
    ├── CLN_2D.py                 # 2D scalar reference
    └── A_ICCG_production.py      # latest 2024-09-17 production
```

### Source Notebooks

```
W:/00_CAE/NGSolve/谷本/
├── 修論/                         # Master's thesis (canonical)
│   ├── 2次元CLN.ipynb
│   ├── CLN_AT.ipynb              ← primary 3D
│   ├── CLN_T-Omega.ipynb
│   ├── CLN_APhi.ipynb
│   └── メッシュ.ipynb
├── 定式_誤差検証/                  # Formulation / error verification
│   ├── 20231211_A_(Penalty)_CLN.ipynb
│   ├── 20231221_A_gauge_CLN.ipynb
│   ├── 20240108_gauge_test.ipynb
│   └── CASE_*.ipynb              # ablations
└── 20240910_静止器回転機用/        # Production for machines
    ├── 20240917_A_ICCG_最新版.ipynb  ⭐ Latest stable
    ├── A-T_formulation.ipynb
    ├── A_CG_TEAM28size.ipynb
    └── Lukas_A-Φ_test.ipynb
```

### Research Project

`W:/30_CauerLadderNetwork/` — CLN + Kelvin research (2026):
- 2026_04_01_長方形CLN/ — current rectangular CLN + Kelvin work
- prior phases for various conductor shapes (sphere, cuboid, etc.)

---

## 10. References

1. **A. Kameari**, "Calculation of transient 3D eddy current using
   edge-elements," IEEE Trans. Magn. (foundational paper for FEM eddy
   current).
2. **谷本** (Tanimoto), Master's Thesis, Kindai University, 2024.
   3D HCurl CLN formulations (A-T, T-Ω, A-Φ) and verification.
3. **K. Sugahara**, "Electromagnetic analysis of eddy current testing
   with Kelvin transformation," IEEE Trans. Magn. 58(9), 2022.
4. **H. Nagamine, T. Yamaguchi, K. Sugahara**, "A Pullback-Based
   Formulation of Kelvin Transformation in Electromagnetic Field
   Analysis," CEFC 2026 (Thessaloniki), id 350.
5. *Foster network synthesis* — classical circuit theory references
   (Brune, Foster, Cauer).
6. **TEAM Workshop Problem 28** — induction levitation benchmark
   used by `A_CG_TEAM28size.ipynb`.

---

## 11. Open Questions / Future Work

1. **CLN + Kelvin for applied uniform field** (§6.3): A-formulation
   reduced + Kelvin breaks down due to A_s unboundedness; H-form or
   T-Ω with Convention B background field is the recommended path.
   Implementation TODO.
2. **Convergence acceleration**: Schmidt drift limits N ≤ 25; could
   re-orthogonalize periodically (Gram-Schmidt-Modified) to extend.
3. **Adaptive mesh + CLN**: standard adaptive refinement via ZZ
   estimator could reduce ladder construction cost; not yet
   integrated.
4. **CLN + Kelvin + adaptive**: combining all three for production
   open-boundary CLN. Major future work.

---

## 12. Related Documentation

- [`docs/kelvin/KELVIN_TRANSFORMATION.md`](../kelvin/KELVIN_TRANSFORMATION.md) — Kelvin transformation theory + (ν-ν₀) form pitfall
- [`docs/peec/`](../peec/) — PEEC integration with FEM (filament source for A-formulation)
- [mcp-server `cln_3d.py`](../../packages/radia-mcp/src/radia_mcp/radia_ngsolve/knowledge/cln_3d.py) — programmatic knowledge access
