# Cauer Ladder Network (CLN) — Tanimoto-Kameari Method

Comprehensive theory + implementation guide for Cauer ladder network
construction from FEM eddy current analysis. Source material: Tanimoto's
2024 master's thesis and Kameari's iterative
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
- **High-order response compression**: use NGSolve `HCurl(p)` as a parent
  space, keep only the Eddy-Visible Response Space (EVRS), eliminate
  Eddy-Invisible DoFs (EIDs / eddy bubbles), and connect the resulting
  `curl(T)` current basis to HDiv-MMM/VIM/BEM/SIBC branches.  See
  [`HCURL_RESPONSE_COMPRESSION.md`](HCURL_RESPONSE_COMPRESSION.md).

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

### 6.1 Status (2026-05-05)

**Important context**: Kelvin transformation in NGSolve is a **proven,
high-accuracy** technique when properly formulated. Kameari's 2025/10/14
presentation demonstrates:
- Magnetic sphere in uniform B with COARSE mesh (Order 3, μ_r=1000):
  **A-Ω_r gives 0.001% error**, **Ω-Ω_r Order 4 gives 0.029% error**
- **Independent of Kelvin radius rk** — even rk=100 (huge) gives same
  accuracy with adaptive refinement (slide 18)
- Three canonical reductions: Ω-Ω_r, A-Ω_r, A-A_r (slide 4)
- Eddy current uses A-φ-A_r (slide 25, TEAM Workshop Problem 7)

**Canonical Kameari weak-form pattern** (slide 7-9):
- Ω-Ω_r: `∫_Ω μ ∇ω·∇Ω = -∫_∂Ω ω B_s·n ds  ∀ω ∈ H¹⁰(order N)`
- A-Ω_r: `∫_Ω μ⁻¹ ∇×N·∇×A = ∫_∂Ω N×H_s·n ds  ∀N ∈ H_curl⁰(order N)`
- A-φ-A_r (eddy): adds `jωσ(N+∇ψ)·(A+∇φ)` and same boundary integral

The **source field is injected via a boundary integral at the inner-Kelvin
interface ∂Ω**, NOT via a `(ν-ν₀)` bulk term. This is structurally
different from the volume-source reduced-A patterns we tried in v11-v14
on the cuboid problem.

**v11-v14 issue revisited (2026-05-05)**: the historical `(ν-ν₀)` bulk
form was NOT Kameari's canonical method. The correct pattern is the
boundary-integral injection. Re-examining v14 against Kameari reference
is on the open work-list.

See [`docs/kelvin/KELVIN_TRANSFORMATION.md`](../kelvin/KELVIN_TRANSFORMATION.md)
§7.5 for the (ν-ν₀) pitfall analysis.

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

### 6.4.1 Kameari's Three Reduction Methods (2025/10/14 reference)

From Kameari, "Electromagnetic Analyses Using Higher Order Hierarchic
Finite Elements" (presentation, 2025/10/14):

| Method | Inner region | Outer (Kelvin) region | Use case |
|---|---|---|---|
| Ω-Ω_r | H = -∇Ω_t | H = -∇Ω_r + H_s | Linear magnetic, H source |
| A-Ω_r | B = ∇×A_t | H = -∇Ω_r + H_s | Mixed, recommended for sphere/cuboid |
| A-A_r | B = ∇×A_t | B = ∇×(A_r + A_s) | Vector source / coil |
| A-φ-A_r | A in air, A+φ in cond | A_r in Kelvin | Eddy current (TEAM7) |

**Convergence (magnetic sphere a=1m, μ_r=1000, B_0=1T, theory=2.9940 T)**:
| Method | Order | Coarse-mesh Bz0 | Error |
|---|---|---|---|
| Ω-Ω_r | 2 | 3.3813 | 12.9% |
| Ω-Ω_r | 3 | 2.9995 | 0.184% |
| Ω-Ω_r | 4 | 2.9985 | 0.029% |
| A-A_r | 3 | 2.9928 | 0.040% |
| **A-Ω_r** | **3** | **2.9940** | **0.001%** ⭐ |

A-Ω_r with Order 3 is the recommended winner.

Source injection is **boundary integral at inner-Kelvin interface**, not
bulk volume term:
```python
# A-Ω_r weak form (slide 8)
a += 1/mu * curl(N) * curl(A) * dx
f += N.Trace() * Cross(H_s, n) * ds(kelvin_interface)  # H_s on boundary
```
where `H_s` is the applied uniform H on the inner-Kelvin boundary
(NOT inside the bulk).

### 6.5 Open-Boundary CLN — When Truncation Works vs When It Doesn't

**Numerically the Kameari iteration is healthy** with the canonical recipe
(nograds=True, NO penalty, tree-cotree, box-shaped outer PEC). Verified
2026-05-03 (`cuboid_521_vacuum_kameari_breakdown.py`, AIR_SCALE=5 box,
ORDER=2, 12 stages):

- All R_n > 0, L_n > 0 (no sign flip)
- Schmidt drift `|⟨J_n, J_m/σ⟩|/⟨J_n, J_n/σ⟩` ranges 3.5e-17 (n=1) to
  9.8e-10 (n=11) — machine precision throughout
- Gram matrix conditioning bounded

The classical "Kameari breaks down at high N" claim is **incorrect for the
iteration itself**; that breakdown is inner-solver origin.

**Physical accuracy depends on geometry / current pattern**:

| Setup | Field decay | Dirichlet truncation? |
|---|---|---|
| 2 parallel cylinders, opposite I (Sugahara 2017 Compumag) | 1/r² (dipole) | ✅ Works at moderate AIR_SCALE |
| Closed magnetic circuit (DC-DC core, 2017 §III.B) | localized | ✅ Works |
| **Isolated conductor + uniform B_z (no return)** | **1/r³ (induced dipole)** | ❌ τ_0 wrong by 9× at AIR_SCALE=5 |
| Coil + open-air return | depends | depends |

For "isolated conductor in applied uniform B_z" specifically (our cuboid
5×2×1 case), AIR_SCALE-sweep showed:
- AIR_SCALE = 3 / 5 / 8 / 12 → τ_0 = 42.5 / 104.4 / 148.3 / 204.6 μs
  (**diverges** as box grows — PEC reflection generates progressively
  longer-wavelength spurious cavity modes)
- vs ELF τ_lead = 11.51 μs → never converges

This is the use case that **Kelvin transformation specifically solves**
(unbounded exterior with proper radiation-like BC). The recent v14
work fixed the historical `(ν-ν₀)` pitfall in reduced-A + Kelvin (see
§6.1); applying the corrected formulation to the cuboid problem is the
ongoing research thread.

**Anti-patterns that look like structural failure but are bugs**:
- `GAUGE_EPS / mu * u * v * dx` penalty term → numerical sign flip
- Sphere-shaped outer air with tetrahedral mesh near cuboid corners →
  irregular elements near sharp interfaces
- Use BOX outer (matches conductor symmetry) and avoid all penalty terms.

**Conclusion**: open-boundary 3D CLN is **not unfixable** — the choice
of method must match the geometry:
- Closed-loop / fast-decay fields → Dirichlet truncation works
  (Sugahara 2017 Compumag is the canonical demonstration)
- Isolated conductor + slow-decay field → Kelvin (or BEM) is required.
  **Kelvin works extremely well** when properly implemented (Kameari
  2025 demos: 0.001% error on magnetic sphere with coarse mesh, no
  dependence on Kelvin radius rk). The key is using the canonical
  weak-form pattern (boundary integral at ∂Ω, NOT (ν-ν₀) bulk term).
  See §6.4.1.
- For canonical analytical benchmarks, Nagamine-style infinite-series
  CLN (CEFC 2026) is a complementary validation tool.

### 6.6 The "AJ vs B²" Inductance Formula Choice

Tanimoto canonical implementations use BOTH formulas:
- **AJ form**: `L_n = R_n × ∫_cond J_n · A_pot dV`
- **B² form**: `L_n = ∫_full curl(A_pot)·curl(A_pot)/μ dV` = magnetic energy

For closed PEC (Dirichlet on conductor surface), they are EQUAL via
integration by parts (no boundary term). For air-box / Kelvin, they
DIFFER by the conductor-surface boundary term.

| Property | AJ form | B² form |
|---|---|---|
| Compute domain | conductor only | full + Kelvin region |
| Sign | can go negative (air-box) | always positive |
| Use in Schmidt update | gives meaningful J update | denominator too large → J_n+1 ≈ J_n (drift→1) |
| With Kelvin region | trivial (no Kelvin in cond) | must integrate B² in pulled-back metric |

**Tradeoff**: B² is physically correct (magnetic energy) but expensive
to compute over Kelvin region. AJ is cheap but can fail (sign or
unphysical Schmidt) when conductor is surrounded by σ=0 region.

The **canonical 20240917 production** (closed PEC) uses AJ. For air-box
+ Kelvin, neither formula has been demonstrated to give clean Cauer
ladder for "uniform B_z applied" — open research question.

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

### 7.1.1 Canonical A-form Recipe (Tanimoto, do not deviate)

Audited from Tanimoto's 2024 A-ICCG and CLN A-T thesis notebooks:

```python
fes = HCurl(mesh, order=order, dirichlet="...", nograds=True)   # OR type1=True
gauge = H1(mesh, order=order, dirichlet="...")

# Bilinear forms — NO penalty term, NO gauge_eps
a_HC += (1/mu) * curl(u) * curl(v) * dx
a_HH += grad(uu) * grad(vv) * dx

for nStage in range(N):
    R = 1/Integrate(J*J/sigma*dx, mesh)

    # Solve A from J
    f += v * J * dx
    solve gfA = inv_HC * f

    # Helmholtz-Hodge gauge correction (REQUIRED, even with nograds=True)
    ff += grad(vv) * gfA * dx
    solve gfu = inv_HH * (-ff)        # rhs sign flip

    # Accumulate (CF expression)
    if nStage == 0:
        Apot = R * (gfA + grad(gfu))   # A is corrected, R-weighted
        B = R * curl(gfA)              # B unaffected by gauge
    else:
        Apot = Apot + R * (gfA + grad(gfu))
        B = B + R * curl(gfA)

    # L (Tanimoto canonical 20240917): AJ form
    L = Integrate(R * J * Apot * dx, mesh)
    # OR Tanimoto修論 (CLN_AT): B^2 form
    # L = Integrate(B*B/mu * dx, mesh)

    # Schmidt orthogonalization
    J = J - sigma * Apot / L
```

**Anti-patterns (avoid)**:
- `a += GAUGE_EPS / mu * u * v * dx` — penalty perturbs ladder values
- Skipping Helmholtz-Hodge correction — divergence accumulates over stages
- Tree-cotree mask + nograds=True simultaneously — over-constrains
- Mixing seed types between stages

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

### Tracked Validation and Provenance

- [`validation_test/cln/`](../../validation_test/cln/) contains the maintained
  CLN smoke cases and committed JSON evidence.
- [`validation_test/maglev/research_cln/`](../../validation_test/maglev/research_cln/)
  preserves the imported research notebooks, Mathematica derivations, and
  historical cross-validation programs used to reconstruct the method.
- [`validation_test/axifem/`](../../validation_test/axifem/) owns the maintained
  axisymmetric element evidence and generated-matrix checks.

The thesis and presentation listed in the references below remain the source
provenance. Local laboratory archive paths are intentionally not part of the
public execution contract.

---

## 10. References

### CLN method

1. **A. Kameari**, "Calculation of transient 3D eddy current using edge
   elements," *IEEE Trans. Magn.*, vol. 26, no. 2, pp. 466-469, 1990.
   Foundational FEM eddy-current paper introducing edge-element
   formulations on which the CLN iterative orthogonalisation builds.
2. **A. Kameari, J. Ebrahimi, K. Fujiwara, Y. Takahashi, N. Takahashi**,
   "Cauer Ladder Network Representation of Eddy-Current Fields for
   Model Order Reduction Using Finite-Element Method," *IEEE Trans.
   Magn.*, vol. 54, no. 3, 7201804, 2018. CLN method introduction
   (the basis of Tanimoto's thesis below).
3. **谷本** (Tanimoto), Master's Thesis, Kindai University, 2024.
   3D HCurl CLN formulations (A-T, T-Ω, A-Φ) and verification.
4. **A. Kameari**, "Electromagnetic Analyses Using Higher Order
   Hierarchic Finite Elements", presentation 2025/10/14. Canonical demonstration of
   Ω-Ω_r / A-Ω_r / A-A_r reductions with Kelvin transformation in
   NGSolve. Magnetic sphere reference: A-Ω_r Order 3 with coarse mesh
   gives 0.001% error, independent of Kelvin radius `r_k`.

### Classical circuit-theory background

5. **W. Cauer**, *Synthesis of Linear Communication Networks*,
   McGraw-Hill, 1958.
6. **R.M. Foster**, "A reactance theorem," *Bell Syst. Tech. J.*,
   vol. 3, pp. 259-267, 1924.
7. **O. Brune**, "Synthesis of a finite two-terminal network whose
   driving point impedance is a prescribed function of frequency,"
   *J. Math. Phys.*, vol. 10, pp. 191-236, 1931.

### Open-boundary side (Kelvin)

8. **K. Sugahara**, "Electromagnetic analysis of eddy current testing
   with Kelvin transformation," *IEEE Trans. Magn.*, vol. 58, no. 9,
   1-6, 2022. A-formulation Kelvin truncation used in §6.
9. **H. Nagamine, T. Yamaguchi, K. Sugahara**, "A Pullback-Based
   Formulation of Kelvin Transformation in Electromagnetic Field
   Analysis," CEFC 2026 (Thessaloniki), id 350. Pullback derivation of
   the (ν − ν₀) reduced-A form used jointly with CLN.
10. **Q. Chen**, "A Review of Finite Element Open Boundary Techniques
    for Static and Quasi-Static Electromagnetic Field Problems,"
    *IEEE Trans. Magn.*, vol. 58, no. 9, 2022. 3D Kelvin theory cited
    by Kameari (Ref [4]) for the open-boundary side.

### Validation

11. **TEAM Workshop Problem 28** — induction-levitation benchmark used
    by `A_CG_TEAM28size.ipynb` for 3D CLN validation.

---

## 11. Open Questions / Future Work

1. **CLN + Kelvin for applied uniform field** (§6.3, §6.5): A-formulation
   reduced + Kelvin breaks down due to A_s unboundedness; H-form or
   T-Ω with Convention B background field is the recommended path.
   Step 0a (no Kelvin, just air-box, 2026-05-05) showed that even WITHOUT
   Kelvin, Kameari struggles when conductor is surrounded by σ=0 region:
   sign flip (AJ form) or Schmidt collapse (B² form) at stage 1.
   Implementation TODO.
2. **AJ vs B² inductance choice** (§6.6): which form converges to clean
   Cauer ladder when conductor is in air/Kelvin? Closed PEC: equivalent.
   Air-box: AJ goes negative (boundary term), B² stays positive but
   produces drift→1 in Schmidt. Hybrid (B² for L_n display, AJ for
   Schmidt update with proper sign control) untested.
3. **Convergence acceleration**: Schmidt drift limits N ≤ 25; could
   re-orthogonalize periodically (Gram-Schmidt-Modified) to extend.
4. **Adaptive mesh + CLN**: standard adaptive refinement via ZZ
   estimator could reduce ladder construction cost; not yet
   integrated.
5. **CLN + Kelvin + adaptive**: combining all three for production
   open-boundary CLN. Major future work.
6. **Nagamine-style analytical infinite-domain CLN as alternative**:
   For canonical geometries (cuboid in vacuum, sphere, etc.), an
   analytical infinite-series CLN (Nagamine CEFC 2026 derivation) avoids
   the FEM+Kelvin+Kameari structural difficulties documented in §6.5.
   When 3D FEM hits open-boundary Kameari barriers, falling back to the
   analytical Nagamine approach for benchmark validation is recommended.

---

## 12. Related Documentation

- [`docs/kelvin/KELVIN_TRANSFORMATION.md`](../kelvin/KELVIN_TRANSFORMATION.md) — Kelvin transformation theory + (ν-ν₀) form pitfall
- [`docs/peec/`](../peec/) — PEEC integration with FEM (filament source for A-formulation)
- [mcp-server `cln_3d.py`](../../packages/radia-mcp/src/radia_mcp/radia_ngsolve/knowledge/cln_3d.py) — programmatic knowledge access
