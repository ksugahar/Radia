# HACApK-accelerated PEEC with PRIMA — paper plan

**Status**: Paper idea (2026-04-16). Not started.

## Co-authorship

- Lead: Sugahara Lab (Radia, PEEC application side)
- **Co-author: Akihiro Ida (伊田明弘)** — HACApK library author, H-matrix
  / BEM-QR / lattice H-matrix expert. Personal collaborator of lead
  (always available for co-authorship). Collaboration unlocks:
  - Direct access to HACApK internals for extensions
  - Credible path to follow-up work on H-QR factorization (below)
  - Top-tier reviewer fit in the H-matrix / BEM community

## Factorization choice: matrix-free PRIMA (corrected)

**HACApK does not provide full-matrix H-LU or H-QR factorization.**
Its public API (`cHACApK_base.h`) offers only:
- `cHACApK_acaplus` — low-rank block compression via ACA+
- `cHACApK_RRQR`, `cHACApK_SVD` — *within-block* rank-revealing helpers
  (NOT hierarchical QR of the full matrix in the Bebendorf-Hackbusch sense)
- `cHACApK_calc_vec` — matvec

Therefore linear solves must be iterative (BiCGSTAB, GMRES) with HACApK
matvec as the operator. There is no factor-once / back-subst fast path.

Our approach:

- **PEEC structure saves us**: `Z(s) = R_diag + sL_H` with R diagonal.
  Lanczos / Arnoldi needs only `L_H @ v` (HACApK matvec) and
  `R^(-1) @ v` (diagonal inverse, O(N) scalar divide).
- **Matrix-free block-Lanczos PRIMA** runs purely on matvec, preserving
  passivity (by Arnoldi / Gram-Schmidt structure) and causality.
- **Multi-frequency rational Krylov**: each expansion point s₀ needs an
  iterative solve of `(R + s₀ L_H) x = b`, accelerated by HACApK matvec
  inside BiCGSTAB. Per-s₀ cost ~ iterations × O(N log N).
- **Passivity-margin bound**: provable from ACA tolerance ε via the
  spectral perturbation bound on the reduced Z_k(s). Q-orthogonality is
  NOT available here, so the bound is slightly weaker than what a
  hypothetical H-QR factorization would give.

## Strategic recommendation (2026-04-16)

Start with **matrix-free PRIMA** (Phase 1). H-matrix QR / LU (Phase 2)
is a contingent extension, not a committed next step.

Rationale, derived from IH usage patterns:

| Usage | Frequency points | Winner |
|-------|------------------|--------|
| Single-frequency design (10 kHz typical IH) | 1 | **matrix-free clearly** (factor-once irrelevant) |
| Skin-effect verification sweep | 10–50 | matrix-free acceptable (~10× slower than H-QR at most) |
| Broadband WPT / resonance search | 100–1000 | H-QR meaningfully faster |
| SPICE netlist extraction | — (reduced model evaluation is cheap) | matrix-free sufficient |

Most IH work is the first row. Matrix-free is the right starting
point; it carries ~1/4 the implementation cost (~1000 LOC vs ~4000
for full H-matrix arithmetic) and already admits a strong paper angle.

### Angles for Phase 1 (matrix-free)

The Phase 1 paper can be framed along any of these, or combined:

- **Angle A (system integration)**: First partitioned Radia-NGSolve
  framework for IH coil-workpiece coupling — HACApK-accelerated PEEC
  plus BEM-SIBC, no air mesh, with guaranteed-passive PRIMA output.
  **Weight**: highest novelty for TMAG / COMPEL venue.
- **Angle B (algorithmic)**: Matrix-free block-Lanczos PRIMA for
  H-matrix operators — structural passivity preservation without
  requiring H-factorization, exploiting the diagonal-R structure of
  PEEC. **Weight**: standalone algorithmic contribution.
- **Angle C (application case study)**: HACApK-PEEC applied to an
  industrial IH coil design — broadband inductance extraction for a
  realistic workpiece. **Weight**: demonstration of practical impact,
  stronger if joined to a collaborator's measurement data.

Recommended composition: **A + B** (methodology + algorithm) as a
single full paper, with C sketched as a case study section.

### When to trigger Phase 2

Only revisit the H-matrix QR / LU investment if Phase 1 benchmarks
show:
- IH users routinely run frequency sweeps with ≥ 100 points, OR
- iterative solve per-s₀ convergence is poor (≥ 200 BiCGSTAB iters
  consistently), OR
- a theory-venue opportunity (SIAM SISC / Numer. Math.) specifically
  values the H-QR analysis angle.

Otherwise, matrix-free remains the deployed path.

## Phased plan

- **Phase 1 (this paper)**: HACApK compression + matrix-free block
  Lanczos PRIMA + BEM-SIBC partitioned coupling + SPICE output.
  Passivity demonstrated numerically; perturbation bound from ACA ε
  stated. Target: CEFC 2026 / IEEE TMAG.
- **Phase 2 (follow-up, implementation-level contribution):
  full H-matrix QR from scratch**

  Context: Ida-san's public "H-QR" work is actually BLR-QR (flat
  2-level, O(N^1.5) memory), not full recursive H-matrix QR. He has
  stated that BLR is a **practical compromise chosen for
  implementation ease**, not the theoretical ideal. Full H-matrix QR
  (Bebendorf-Hackbusch recursive form) retains the O(N log N) scaling
  of H-matrix compression throughout the factorization.

  The ambition of Phase 2 is to implement **full H-matrix QR directly
  on top of HACApK 1.0.0's primitives**, without falling back to BLR.

  ### Algorithmic skeleton (recursive)

  ```
  function H_QR(A: H-matrix):    # block-column QR
    [A1 | A2] ← split columns
    (Q1, R11) ← H_QR(A1)                   # recurse
    R12 ← Q1^T H@ A2                        # H-matmul, re-compress
    A2⊥ ← A2 H- Q1 H@ R12                   # H-add, re-compress
    (Q2, R22) ← H_QR(A2⊥)                   # recurse
    return [Q1 | Q2], [[R11, R12]; [0, R22]]
  ```

  ### Primitives already in HACApK 1.0.0 (we reuse, do not reinvent)

  - `cHACApK_acaplus` — ACA+ compression (for re-compressing products
    and sums of H-blocks)
  - `cHACApK_RRQR` — rank-revealing QR inside blocks
    (orthogonalization + truncation steps)
  - `cHACApK_SVD` — alternate rank truncation
  - `cHACApK_calc_vec` — matvec for residual checks
  - `st_cHACApK_leafmtxp`, cluster tree infrastructure

  ### What must be built on top

  - H × H matrix multiplication with re-compression after each product
  - H + H, H − H with re-compression
  - Modified Block Gram-Schmidt at each recursion level for numerical
    stability
  - Triangular solve (H-triangular × vector) for back-substitution
  - Error accumulation control via tolerance tuning across recursion
    depth

  ### Scope estimate

  | Item | C++ LOC (rough) |
  |------|-----------------|
  | H × H multiplication + re-compression | 600 |
  | H ± H with re-compression | 300 |
  | Block Gram-Schmidt driver | 400 |
  | Recursive H_QR core | 500 |
  | Triangular H-solve (back-sub) | 400 |
  | OpenMP parallelization | 500 |
  | Tests (analytical problems, convergence) | 1000 |
  | **Total** | **~3500-4000 LOC** |

  Feasible as a focused 2–3 month effort. AI-assisted implementation
  (Claude) is realistic given the algorithmic clarity and primitive
  reuse; the mathematical recipe is public (Bebendorf 2008 textbook
  Ch. 4, Hackbusch 2015 Ch. 7 both give pseudocode).

  ### Add H-LU in parallel (~15% extra effort)

  Once the H-matrix arithmetic backbone exists (H×H, H±H with
  re-compression, H-triangular solve), adding H-LU costs only the
  recursive block LU driver with Schur complement updates and
  block-level partial pivoting — ~500 LOC on top of ~3500 for H-QR.

  **Implementing both enables a head-to-head comparison** for the
  PEEC shifted system `(R_diag + s₀ L_H)`:

  | | H-LU (predicted) | H-QR (predicted) |
  |--|------------------|------------------|
  | Speed | Faster (simpler arithmetic) | 2-3× slower |
  | Stability | Pivot may fragment block hierarchy on near-singular cases | Unconditionally stable via orthogonality |
  | Passivity-margin bound | Weaker (product bound with pivoting) | Tighter (Q-orthogonality) |
  | PEEC-specific behavior | May be fine since R is diagonally dominant | Fine by construction |

  **Which actually wins on real IH coil geometries is unknown.**
  The experimental comparison is itself a paper-worthy contribution:
  no published study compares H-LU vs H-QR on PEEC shifted systems.

  ### Papers unlocked

  - **Algorithmic paper A**: *H-matrix QR on HACApK primitives:
    passivity-margin bound and PEEC application.* SIAM SISC / Numer.
    Math.
  - **Algorithmic paper B** (or extended A): *H-LU vs H-QR for
    complex-symmetric shifted systems: a PEEC benchmark study.*
    IEEE TMAG / COMPEL algorithmic slot.
  - **Application paper** (Phase 1, extended): broad-band IH coil
    sweep with H-matrix factorization — O(factor + f_points × N log² N),
    strictly better than BLR's O(factor + f_points × N^1.5).

  ### Open questions for Ida-san

  - Would he collaborate / co-author on the algorithmic paper(s)?
  - Is contributing the new H-QR / H-LU back to upstream HACApK
    welcome, or does his lab have strategic reasons to keep advanced
    factorizations out of the public library?
  - Are there existing unpublished attempts in his lab we should
    build on rather than starting from scratch?

## Note on upstream HACApK 1.0.0

Verified 2026-04-16 on `S:\<external>\ppOpenHPC-MATH-HACApK\src\HACApK_1.0.0\`
(external upstream) and `src/ext/HACApK/` (Radia bundled): the public
upstream release provides:
- `HACApK_adot_*` (matvec, serial + hybrid MPI+OpenMP)
- `HACApK_bicgstab_lfmtx` (BiCGSTAB with H-matvec)
- `HACApK_gcrm_lfmtx` (GCR(m) with H-matvec)
- `HACApK_RRQR`, `HACApK_SVD` (within-block low-rank helpers — NOT
  full-matrix H-QR in the Bebendorf-Hackbusch sense)

Neither H-LU nor full-matrix H-QR is in 1.0.0. BLR-QR (above) is the
path through Ida-san's research codebase / HACApK-MAGMA. The public
1.0.0 is iterative-only as shipped.

## Novelty claim

We present the first integrated PEEC inductance extractor combining:
1. **HACApK H-matrix compression** of the partial-inductance matrix L
   (O(N log N) assembly and storage),
2. **Matrix-free rational Krylov PRIMA** exploiting the fact that the
   PEEC resistance matrix R is diagonal, so the reduction runs on
   HACApK matvec + diagonal inverse without requiring H-LU,
3. A SPICE-ready reduced rational Z(s) = R_k + sL_k at the port level,
4. **Coupling to NGSolve-BEM-SIBC** for the workpiece, so that the IH
   coil + workpiece system is solved without an air mesh.

Although H-LU is known in the literature (Bebendorf 2005, Grasedyck
2003), it is **not implemented in the HACApK library** we rely on.
Sidestepping the H-LU dependency via the diagonal-R structure of PEEC
is itself a methodological contribution.

## Gap analysis

| Tool / paper | Compression | PRIMA | SPICE | Workpiece coupling |
|--------------|-------------|-------|-------|--------------------|
| **FastHenry** | FMM | no | no | no (air mesh via FastModel) |
| **FastCap / pFFT++** | FMM (precorrected) | no | no | — (capacitance only) |
| **Ansys Q3D** | FMM + SDM | likely (undocumented) | yes | lumped |
| **ngsolve.bem** | **dense** (no H-matrix) | no | no | yes (BEM only) |
| **Radia PEEC (this)** | **H-matrix (HACApK)** | **yes (matrix-free block-Lanczos)** | **yes** | **BEM-SIBC** |

Key novelty: **H-matrix is preferable to FMM in combination with PRIMA.**
FMM supplies only matvec via tree walk and cannot persist state between
calls; every PRIMA Lanczos iteration repeats the full O(N log N) walk.
HACApK's low-rank block storage stays in memory across the k Krylov
iterations, so the compressed matrix is paid for once and reused. For
IH coil design with broadband (DC–MHz) sweeps of 10²–10³ frequency
points, this matters.

For the BEM-SIBC side, ngsolve.bem is currently dense. Because the
workpiece surface mesh is small (O(10³)–O(10⁴) DOFs after SIBC reduces
the volume to a surface), dense BEM is acceptable. The H-matrix
acceleration matters only on the coil side, where nwinc × nhinc
subdivision can push filament count to 10⁴–10⁵.

## Passivity and causality: structural guarantees

PRIMA is not just "a reduction method" — it is the **Passive Reduced-order
Interconnect Macromodeling Algorithm** (Odabasioglu et al. 1998). Built-in
guarantees distinguish it from ad-hoc alternatives:

| Property | PRIMA | Vector fitting | SVD truncation | Balanced truncation |
|----------|-------|----------------|----------------|---------------------|
| **Passivity preserved** | yes (by Arnoldi structure) | no (needs post-hoc fix) | no | yes (if dissipative) |
| **Causality preserved** | yes (analytic in RHP) | conditional | no | yes |
| **Stability preserved** | yes (no RHP poles) | conditional | no | yes |
| **Cost per reduction** | O(N log N) × k (matvec) | Fit per s-point | O(N³) SVD | O(N³) Lyapunov |
| **SPICE-ready** | directly | after passivity fix | no | directly (complex) |

Why this matters for IH / WPT papers:

1. **Lossless systems, lossless output**: a passive original (RLC with
   symmetric SPD matrices) produces a passive reduced model. The reduced
   impedance Z_k(s) satisfies Re(Z_k(jω)) ≥ 0 for all ω, by construction.
2. **Direct time-domain co-simulation**: the reduced SPICE netlist can
   drive LTspice / ngspice transient simulation without oscillations or
   numerical blow-up — a failure mode common with ad-hoc rational fits.
3. **Extrapolation safety**: PRIMA Krylov subspace matches moments of the
   original transfer function at the expansion point s₀. Away from s₀,
   deviation grows gracefully (power-law), not catastrophically (poles
   wandering into RHP).
4. **Certifiable designs**: for safety-critical applications
   (medical IH tempering, automotive WPT charging), a passivity proof
   on the reduced model is a regulatory asset.

We claim in the paper:
- The full pipeline — HACApK compression → matrix-free Arnoldi (PEEC
  structure) → PRIMA projection → SPICE — preserves passivity and
  causality at every stage.
- HACApK compression error (ACA tolerance ε) bounds the reduced-model
  passivity margin; an explicit bound `Re(Z_k) ≥ -C·ε` is derived
  (Lemma target).
- Empirical validation: on test cases 1–6 below, the reduced Z_k(s) is
  strictly passive (positive-real) within numerical tolerance.

This "素性の良さ" (structural soundness) is the headline methodological
contribution, beyond raw speedup numbers.

## Implementation architecture (code reuse)

The current Radia `RadHACApKManager` (`src/core/rad_hacapk.h`) is bound
to `radTInteraction*`, i.e., tightly coupled to the MSC/MMM kernel.
~70% of the workflow is nevertheless kernel-agnostic. The cost-effective
refactor is a base class that exposes a single virtual `ComputeEntry`:

| Component | Kernel-agnostic? |
|-----------|------------------|
| Cluster tree (geometric) | **yes** |
| Admissibility check | **yes** |
| ACA+ block compression | **yes** |
| H-matrix matvec | **yes** |
| BiCGSTAB + Block-Jacobi preconditioning | **yes** |
| `ComputeEntry(i,j)` | **no (kernel-specific)** |
| Element metadata (face/filament geometry) | **no** |

Refactor target:

```cpp
class RadHACApKBase {                 // shared ~700 LOC
    virtual double ComputeEntry(int i, int j) = 0;
    // everything else shared
};

class RadHACApKMSCManager : RadHACApKBase { ... }   // existing, lift from radTInteraction
class RadHACApKPEECManager : RadHACApKBase { ... }  // new, ~400 LOC
```

**Revised implementation cost for Phase 1**: ~400 LOC of new PEEC code,
plus a one-time ~300 LOC refactor of the existing MSC manager. The
refactor is valuable independently: it also enables future applications
(BEM-SIBC acceleration, domain-decomposition coupling) to plug into
the same H-matrix backbone.

## PRIMA is a feature, not the main solver

For Radia PEEC, the primary solver path is **direct iterative solve
(BiCGSTAB + HACApK matvec) per frequency**. PRIMA layers on top for
specific use cases:

| Use case | Primary solver |
|----------|----------------|
| Single-frequency design (10–500 kHz IH, 6.78 MHz WPT) | **Direct iterative** |
| Frequency-sweep verification (10–100 points) | Direct iterative × points |
| SPICE netlist extraction | **PRIMA** (main value) |
| Transient / time-domain co-simulation | **PRIMA** (low-order required) |
| Broadband WPT resonance (100–1000 points) | PRIMA (or H-factor in Phase 2) |
| BEM-SIBC coupled solve | Direct iterative inside Picard outer loop |

Implementation order:

1. **HACApK-PEEC adapter + BiCGSTAB direct solve** (~400 LOC new,
   ~300 LOC refactor). Covers 80% of day-to-day use cases. Gives
   10× speedup at N = 18k filaments.
2. **PRIMA with HACApK matvec** (~200 LOC, wires `lanczos_reduction.py`
   into the H-matrix operator). Adds SPICE / transient support.
3. Phase 2 H-factorization (contingent on Phase 1 measurements).

## Method outline

### 1. PEEC filament discretization
- Existing in Radia: `PEECBuilder`, `fasthenry_parser`, nwinc/nhinc
  cross-section subdivision. No change.

### 2. H-matrix compression of L
- Kernel: Biot-Savart `A = (μ₀/4π) ∫ J/r dV`, so `L_ij ∝ ∫∫ 1/r` —
  smooth kernel, ACA-friendly.
- Reuse `src/core/rad_hacapk.cpp` HACApK wrapper (already integrated for
  MMM/MSC). Add a PEEC adapter: cluster tree on filament centroids,
  admissibility by filament-pair distance vs. cluster diameter.
- Target compression: < 5% memory vs. dense at ACA ε = 10⁻⁴.

### 3. H-matrix LU and back-substitution
- HACApK supports H-LU (confirm via library audit).
- At Ref frequency f₀: factor `Z(s₀) = R + j·2πf₀·L_H` into H-LU once.
- For all other frequencies: back-substitute (perturbed-LU or
  direct-LU on each f if H-LU memory is acceptable).

### 4. Block-Lanczos / PRIMA
- Reuse `src/radia/lanczos_reduction.py` (8937 lines, already has
  `LCResonantPRIMA`, `HierarchicalReducer`, `LanczosReducer`).
- Input: HACApK matvec `L_H @ v` and `R @ v`.
- Output: `Z_reduced(s)` of order k ≪ N, with passivity preserved.
- k depends on port count p and target bandwidth (k ~ 10p typically).

### 5. BEM-SIBC workpiece coupling (partitioned, NOT monolithic)

The coil side (HACApK H-matrix) and the workpiece side (dense
ngsolve.bem) are kept **block-decoupled**. No monolithic mixed
H-matrix + dense system is assembled.

```
Outer Picard loop (k = 0, 1, 2, ...):
  1. Coil PEEC solve:   (R + jωL_H) I_coil^{k+1} = V_drive − V_back^{k}
                        [HACApK H-QR, fast]
  2. Biot–Savart:       φ_inc = BiotSavart(I_coil^{k+1}) on workpiece surface
                        [HACApK matvec on filament kernel]
  3. BEM-SIBC solve:    (DL − SL·M^{-1}·K) σ = φ_inc
                        [ngsolve.bem dense, PARDISO, small N_surf]
  4. Back-reaction:     V_back^{k+1} = induced voltage on coil from σ
                        [BiotSavart from surface current]
  5. Converged?         ‖I_coil^{k+1} − I_coil^{k}‖ < tol → exit
```

Convergence:
- Non-magnetic workpiece (copper, aluminum, brass): **3–5 Picard
  iterations** typical (weak feedback).
- Magnetic workpiece (steel, stainless): may need Anderson
  acceleration or relaxation (up to ~10 iterations).
- Per memory `bem_coupled_solver_existing.md`, the existing
  `CoupledBEMSolver` validated against `calc_fem_kelvin` within +0.3%
  (copper, 50 kHz) already uses an equivalent partitioned structure.

Advantages of partitioned vs. monolithic:
- Each block stays in its natural solver (HACApK vs. dense PARDISO).
- No mixed-format linear algebra to invent.
- Reuses existing `CoupledBEMSolver` plumbing — the paper's new
  contribution is on the coil side (HACApK + PRIMA), not on the
  coupling itself.
- Passivity preservation inherits cleanly from each block.

PRIMA is applied on the coil block (coil PEEC reduced order k ≪ N_coil).
The workpiece is not reduced — it is already small (N_surf) and its
transfer function feeds into the coil block as a modified source term
V_back(jω), evaluated per frequency point.

## Validation benchmarks

| # | Configuration | Reference | Target |
|---|---------------|-----------|--------|
| 1 | Single round loop, analytical | Biot-Savart + Rosa formulas | < 0.1% on L at DC |
| 2 | Two concentric coils | Mutual inductance series | < 0.5% on M |
| 3 | Air-core multi-turn solenoid | FastHenry (small N) | < 0.5%, benchmark speedup 10×+ at N > 10⁴ |
| 4 | Pancake coil over billet (IH) | `calc_fem_kelvin` (FEM-SIBC) | < 1% on L, R, P at 10 kHz–1 MHz |
| 5 | Full frequency sweep (DC–10 MHz, 200 points) | HACApK single-f repeated | > 20× speedup with PRIMA |
| 6 | SPICE time-domain transient | LTspice ngspice | Waveform ~ 1% envelope |

## Venues

- **CEFC 2026** (June, Thessaloniki) — already on project roadmap per
  memory. MMM/MSC H-matrix paper was the original slot; this work is
  larger scope, may need a separate follow-up paper at IEEE COMPEL or
  TMAG.
- **NGSolve User Meeting** — fitting venue for the BEM-SIBC coupling
  component; live demo possible.
- **IEEE TMAG** — inductance extraction methods paper.
- **ICS Newsletter / NGSolve forum** — short announcement once open-source.

## Implementation sub-phases

1. **HACApK-PEEC adapter** (~1-2 weeks): C++ wrapper for PEEC L kernel.
   Unit tests vs. dense dense L for small N.
2. **H-LU audit** (~1 day): confirm Radia's HACApK wrapper exposes LU;
   if not, add it (or use iterative + H-matrix precond as fallback).
3. **PRIMA integration** (~1 week): wire `lanczos_reduction.py` to
   consume HACApK matvec.
4. **Benchmark suite** (~2 weeks): 6 tests above, generate paper figures.
5. **Paper draft** (~3 weeks): 8-10 pages, double column.

Total: ~2-3 months effort for one competent implementer + reviewer.

## Risks and open questions

- **H-LU accuracy vs. rank**: H-LU can lose accuracy with deep
  factorization. Fallback: use HACApK matvec with iterative solver
  (BiCGSTAB + sparse preconditioner) per frequency; PRIMA still provides
  the factor-once benefit via Krylov.
- **nwinc / nhinc scaling**: for aggressive subdivision (5×5), N grows
  25× per turn — does HACApK compression still help? Needs measurement.
- **Proximity effect accuracy**: filament count must be high enough;
  cross-check against Dowell's analytical formula on simple bars.
- **Comparison with Q3D**: accessing a Q3D license for apples-to-apples
  benchmark is hard. Plan B: compare against published Q3D numbers
  (e.g., the IEEE 2020 power-module paper) on the same geometry.

## Relationship to other Radia work

- Builds on existing: `rad_hacapk.cpp`, `peec_matrices.py`,
  `peec_coupled.py`, `lanczos_reduction.py`, `bem_coupled_solver.py`.
- Independent of the `FORCE_COMPUTATION_DESIGN.md` line of work — they
  share the "Radia analytical kernel + NGSolve high-order" philosophy
  but address different quantities (force vs. impedance).
- **Not started.** This document is a research plan, not an
  implementation checklist.
