# Paper outline — SF inverse design with monotone Path-A compensation

The session 2026-05-30 produced enough material for a journal-grade
paper.  This page is the skeleton for the manuscript draft.

## Working title

> *Open-source kernel-agnostic stream-function coil design with monotone
> Path-A compensation, regularisation choices, and surface-deformation
> bilevel optimisation*

Authors: K. Sugahara et al. (Sugahara Lab, Kindai University).

Target venues (in order of fit):

1. **IEEE Transactions on Magnetics** — methodology paper, broad
   community fit
2. **IEEE Transactions on Medical Imaging** — if we lead with MRI
   gradient/shim coil application
3. **Compumag 2027** — workshop version of the methodology
4. **IEEJ Joint Tech Meeting Static Apparatus / Rotating Machinery** —
   SA-25-020 follow-on (already established lineage)

## Candidate contributions — calibrated by what is genuinely novel

The list below is the set of points worth investigating for the
manuscript.  Each needs a literature search to confirm before final
wording.  ACA+TSVD itself is NOT novel — Bebendorf 2000 (original ACA),
Bebendorf & Rjasanow 2003 (ACA+), Hackbusch 2008 *Hierarchical
Matrices* (TSVD recompression / "round" operation) are all standard
H-matrix algorithm primitives.  The SA-25-020 paper (2025) put ACA+TSVD
together with the SF coil-design pipeline and published it; even that
combination probably has international precedent (the MRI gradient
coil community has used H-matrix acceleration; the specific
combination needs to be checked against e.g. Kurz-Rain-Rjasanow IEEE
TMag and the MRI literature before claiming originality).

What the current session adds on TOP of SA-25-020 is what may have
real novelty:

1. **Path-A compensated iteration on FE-direct ψ** with empirical
   monotone convergence on the planar uniform Bz benchmark (iter 40–47,
   residual 0.62 % → 0.49 %).  Extends the Kuijpers 2023 "deviation =
   field error" observation from a SELECTION criterion to a SOLVE
   update — Kuijpers et al. observe the coupling and pick the
   least-bad chain; we fold the chain's parasitic field back into the
   inverse-problem RHS and iterate.  Worth checking literature for
   prior compensated-iteration work in coil design.

2. **Representation-dependence of Path-A** — the same fixed-point map
   only converges when ψ is a continuous FE function; on the grid-
   sampled / matplotlib-contoured basis-loop representation it
   oscillates.  This is a methodological observation that has not (to
   our knowledge) been pointed out in print.

3. **Empirical complexity-tier framework** (Easy / Medium / Hard /
   Harder) for SF coil design with measured RMS bounds per tier.
   Probably an articulation of community-implicit knowledge rather than
   a new discovery; still useful as a practical design guide.

4. **Material-kernel extension via the kernel-agnostic callback
   contract**, demonstrated on a shielded coil with iron back plate
   using Radia MMM through the same callback interface.  *(Pending the
   shielded-coil benchmark.)*  The callback abstraction itself is a
   standard H-matrix library pattern (HTool, HLib, H2Lib all support
   it), so the contribution here is the SF-coil-design APPLICATION on
   magnetic materials, not the contract itself.  Originality against
   commercial FEM equivalents to be verified.

5. **Regularisation choice catalogue + bilevel deformation** — L²
   (Euclidean), H¹ (current density), 1/σ H¹ (true ohmic dissipation),
   inductance-diagonal proxy, L∞ peak (experimental); composed with
   Optuna CMA-ES outer loop on surface geometry.  Each individual
   piece has precedent (Liu-Hennig-Korvink H² regularisation, Forbes-
   Crozier inductance min, SA-25-020 CMA-ES); the integrated pipeline
   that lets the user pick + compose is a software contribution, not a
   numerical-methods novelty.

6. **Architectural alignment with NGSolve 6.2.2604+ `ngsolve.bem`
   HTool bridge** — same block-entry callback contract, so the matrix-
   assembly layer can be swapped for the H-matrix-accelerated operator
   when applicable (surface-to-surface only; off-surface targets stay
   on our path).  Compatibility observation, not a contribution.

### Honest priority of which to lead with

The strongest candidate for a methods-paper headline is **(1) Path-A
compensated iteration with monotone convergence demonstration on
FE-direct ψ**, because the Kuijpers 2023 paper at Compumag is recent
and our extension is a specific algorithmic update to it.  The
representation-dependence observation (2) is the supporting result.

The OTHER items (4 / 5 / 6) are best framed as **implementation
contributions**, not methods contributions — see next section.

## Implementation contribution (= separate from algorithmic novelty)

Even though the underlying numerical methods (ACA+TSVD, FE Galerkin
BEM, CMA-ES) are standard, the integrated software stack we ship is
a contribution in its own right.  This is the framing for a software-
focused paper (JOSS, SoftwareX) or for the Implementation section of
a methods paper.

**The integration combines three libraries that are valuable
separately but have not been combined for SF coil design (pending
literature check)**:

  - **NGSolve `ngsolve.bem`** (Schöberl et al. 2024+; HTool bridge by
    Pierre Marchand 2026-04-20): surface BEM operators (Laplace SL/DL,
    Maxwell SL/DL, Biot-Savart with FMM), arbitrary-order H1 / HCurl /
    HDiv / HDivSurface FE basis on surfaces, Sauter-Schwab quadrature,
    block-entry callback API (`CalcSubMatrix`, `CalcSubMatrixCapsule`).
  - **HACApK** (Ida et al., ppOpen-HPC, MIT 2015+): battle-tested
    ACA+ + H-matrix solver, kernel-agnostic via `HACApK_set_entry_func`
    callback.  Already integrated into Radia for the MMM/MSC system
    matrix.
  - **Radia** (ESRF + Sugahara Lab): MMM/MSC material kernels
    (permanent magnet, soft iron, BH curves, SIBC), PEEC inductance
    extraction, single-stroke chain construction, Path-A iteration,
    Biot-Savart segment evaluator, kernel-agnostic
    `radia.stream_function.aca_tsvd` entry callback.

**The integration value**:

  1. **Callback contract alignment** — `ngsolve.bem.CalcSubMatrix`
     (Marchand 2026), `HACApK_set_entry_func` (Ida 2015), and
     `radia.stream_function.aca_tsvd(entry)` (Sugahara 2025) are
     all `entry(i, j) -> float` (or block-form) compatible.  This
     alignment was lucky-or-foresighted upstream; we provide the
     glue that uses it for SF coil design.
  2. **Kernel-swap transparency** — a user can run the same
     pipeline on free-space Biot-Savart (default), Radia MMM with
     iron yoke, SIBC workpiece (Karl iteration), or future
     application-specific kernels by replacing the `entry` function.
     The single-stroke chain, Path-A iteration, regularisation
     selection, and deformation outer loop are unchanged.
  3. **End-to-end OSS pipeline** — to our knowledge no other OSS
     SF coil designer combines kernel-agnostic SF inverse design +
     single-stroke chain construction + material-kernel support
     in one package (CoilGen is MRI-only and free-space-only;
     commercial Comsol AC/DC + Opt Module is comparable in scope
     but closed source and licence-gated).
  4. **Active upstream alignment** — NGSolve 6.2.2604+ shipped the
     `ngsolve.bem` HTool bridge that our `aca_tsvd` callback
     contract is structurally compatible with.  Replacement of the
     interim per-target `LinearForm` matrix assembly with a true
     H-matrix-backed operator is a contained change in our
     pipeline (the rest stays).

**Where this gets credited**:

  - Methods paper: short Implementation section near the end,
    crediting Schöberl, Marchand, Ida et al., Sugahara Lab; framing
    as "we integrate" not "we invent".
  - Software paper (JOSS / SoftwareX): the integrated pipeline
    is the contribution, with the demos and validation as the
    evidence of usefulness.  Methods novelty (Path-A) is
    secondary; software completeness + reproducibility +
    documentation are primary.

The two framings are not mutually exclusive — a methods paper +
companion JOSS paper is a common pattern.

## Section structure

### 1. Introduction (1.5 pages)

  - SF method history (Turner 1986 → Forbes-Crozier → Peeren → Liu-
    Hennig-Korvink → Kuijpers 2023).
  - The single-stroke manufacturability problem.
  - The deviation → field-error coupling (Kuijpers 2023).
  - Our contribution: FOLD the deviation INTO the solve, validated as
    monotone convergence with FE-direct ψ.

### 2. SF method recap (1 page)

Standard material — surface current ⊥ ∇ψ, Bz target = least-norm
inverse problem, ACA+ + TSVD pseudo-inverse.

Reference our existing [theory.md](theory.md) for full detail.

### 3. FE-direct ψ formulation (1 page)

  - H¹(Γ) FES on Netgen triangulation.
  - Per-target `LinearForm` for matrix entry.
  - Closed-form Lagrangian for min-seminorm regularisation.
  - Why this is essential: continuous ψ → smooth contour family →
    Path-A converges.

### 4. Path-A compensated iteration (1.5 pages)

  - Fixed-point ψ ← ψ + α · pseudo_inverse(B_target − I_w · Bz_chain_unit)
  - Convergence theorem (informal): if `B_c(ψ)` is Lipschitz in ψ,
    Picard contracts for `α` small enough.
  - Empirical: basis-loop ψ → topology jumps → not Lipschitz → no
    convergence (best-effort).  FE-direct ψ → smooth → monotone
    convergence iter 40-47.
  - Figure: residual vs iter for both representations.

### 5. Regularisation choices (1 page)

Table from [regularization.md](regularization.md).  Note the
non-monotone behaviour in FE polynomial order p (sweet spot at p=3
for our mesh density) — interpret as discretisation-matched
resolution.

### 6. Surface deformation outer loop (1 page)

  - Bilevel: inner SF + Path-A (cached factorisation, fast); outer
    CMA-ES via Optuna over geometry params.
  - Measured: bump-only -63 % RMS in 20 trials, 23 s wall.
  - Design rule: turn deform OFF when baseline already sub-0.5 %.

### 7. Complexity tier framework (0.5 page)

Tier table from [single_stroke.md](single_stroke.md).  Recommended
escalation paths.

### 8. Validation (2–3 pages)

**This is what needs the 2-week campaign before submission**.

Benchmarks (from [benchmarks.md](benchmarks.md)):

  - Bilac et al. planar shim (TODO)
  - Turner cylindrical Gz analytical (TODO)
  - Lemdiasov-Ludwig 2005 (TODO)
  - CoilGen head-to-head (TODO)
  - Shielded coil with iron back plate via Radia MMM (TODO,
    material-kernel demo)

Each benchmark → row in a comparison table → discussion.

### 9. Open-source release (0.5 page)

  - PyPI: `pip install radia[cubit] radia-mcp`
  - GitHub: ksugahar/Radia
  - Documentation: this `docs/stream_function/` folder
  - Reproducibility: every benchmark + figure has a script + JSON
    result.

### 10. Discussion / future work (0.5 page)

  - True inductance min via ngsolve.bem MaxwellSL on a thin shell.
  - Bilevel σ as design variable (let CMA-ES choose σ(x, y) too).
  - Multi-objective Pareto (RMS / inductance / peak current / wire
    length) via Optuna NSGA-II.
  - Multivalued-potential reformulation (D-path in
    `streamfunction(topic=single_stroke)`).
  - GPU acceleration of the per-entry kernel (CUDA Biot-Savart for
    massive-M targets).

### 11. Conclusions (0.25 page)

## Figures (must-have)

  1. Pipeline schematic: target → SF solve → ψ → contours → single-
     stroke spiral → CAD → PEEC → field eval.
  2. Path-A monotone convergence residual plot (FE-direct vs basis-
     loop on the same problem).
  3. Complexity tier table (graphic).
  4. Single-stroke method comparison: greedy vs lobe vs kuijpers
     chain visualisation on Gx fingerprint.
  5. Regularisation choice comparison (ψ contours for each mode).
  6. Deformation outer-loop CMA-ES progression (cost vs trial).
  7. Benchmark table (Bilac, Turner, Lemdiasov-Ludwig, CoilGen,
     shielded coil) — the validation table.

## References (preliminary)

  - Turner, R. *J. Phys. D: Appl. Phys.* 19, L147 (1986).
  - Forbes, L.K., Crozier, S., various 2002-2010 inverse design.
  - Peeren, G., *Stream Function Approach for Determining Optimal
    Surface Currents* (2003 PhD thesis, TU/e).
  - Lemdiasov, R.A., Ludwig, R. *Concepts Magn Reson Part B* 26B(1),
    67-80 (2005).
  - Liu, F., Hennig, J., Korvink, J.G. *IEEE Trans Magn* 48(4), 1179
    (2012) — H² smoothness regularisation.
  - Kuijpers, B.J.A., Jansen, J.W., Lomonova, E.A. *Compumag 2023* [525].
  - Schwartz, A., et al. *CoilGen*
    (github.com/Philipp-MR/CoilGen).
  - HACApK (ppOpen-HPC, MIT).
  - NGSolve (Schöberl et al.).
  - Sugahara Lab, *SA-25-020* — (ACA+)+TSVD + CMA-ES lineage.

## Submission gate

Before submission:

  - [ ] 5 validation benchmarks promoted from the target list in
        `docs/stream_function/benchmarks.md`, each with JSON results
  - [ ] Comparison table assembled
  - [ ] All figures generated by scripts + checked-in PNGs
  - [ ] CHANGES.md or release notes for the relevant Radia version
        bump
  - [ ] PyPI release tagged
  - [ ] arXiv preprint posted (proof-of-priority — concurrent with
        journal submission)

## Timing note

NGSolve 6.2.2604 (April 2026) shipped the H-matrix bridge that other
groups working in SF + BEM will likely build on over 6-12 months.  If
proof-of-priority on any specific contribution is intended, the
manuscript draft + arXiv preprint are worth doing sooner rather than
later (target: ~1 month for arXiv, 2-3 months for journal submission).

## Cross-reference

  - All sub-docs in this folder
  - MCP topic: `streamfunction(topic=session_2026_05_30)` section 11 +
    `streamfunction(topic=single_stroke)` for full background
