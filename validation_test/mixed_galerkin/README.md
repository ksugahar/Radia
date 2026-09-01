# Mixed Galerkin formulation (CLN bulk + HOIBC surface)

Research scripts for the **mixed Galerkin** formulation of eddy-current
admittance Y(s): bulk Cauer Ladder Network (CLN) at low frequency +
Higher-Order Impedance Boundary Condition (HOIBC) at high frequency,
coupled via the Schur complement (= discrete Steklov-Poincaré operator =
discrete DtN map at the algebraic level).

> **History note (2026-06-12)**: This directory was previously
> `examples/hierarchical_cauer_sibc/` and centred on the **Warburg-Schur
> termination** (`Y_R = Y_CLN + K_SIBC √s / (s + d)`, with `d` tuned).
> That approach was superseded by the Mixed Galerkin framework below,
> which removes the `d` parameter entirely and improves wall-band
> accuracy by 1–4 orders of magnitude. The Warburg-Schur code was
> hard-deleted; see
> `memory/project_warburg_schur_deprecated_2026_06_12.md` for the
> history and lessons.

## Framework summary

The variational problem on a conductor of arbitrary shape is

  (-∇² + s μ σ) v = -s μ σ    on Ω,  v = 0 on ∂Ω

with the volumetric admittance

  Y(s) = V σ (1 + ⟨v⟩)

decaying from Y(0) = Vσ at DC to Y(s) ~ K_SIBC/√s at deep skin.

The two endpoints have different natural bases:

| Endpoint | Basis | Captures |
|---|---|---|
| s → 0 (DC, low-freq) | bulk CLN Foster modes (sine eigenfunctions) | volume diffusion |
| s → ∞ (deep skin)    | HOIBC Senior tower of fractional-power envelopes | surface skin effect |

The **mixed Galerkin** combines both:

  v(r, s) ≈ Σ ξ_k^{bulk} φ_k(r) + Σ ξ_k^{surf} ψ_k(r, s)

with bulk φ_k from CLN Krylov-at-s=0 (frequency-independent) and surface
ψ_k(r, s) the planar SIBC envelope plus Senior tower curvature
corrections. The Galerkin system has size (N_bulk + N_surf) and is
solved per-frequency. Intermediate-frequency accuracy is bounded by the
**two-point Padé** approximation in t = √s.

## Directory layout

| Subdir | Contents |
|---|---|
| `_references/` | Remaining square/cube analytic references pending src/API promotion |
| `cylinder/`    | Mixed Galerkin on infinite-z cylinder cross-section |
| `sphere/`      | Mixed Galerkin on solid sphere |
| `square2d/`    | Mixed Galerkin on infinite-z square cross-section |
| `cube3d/`      | Mixed Galerkin on solid cube (rank-N + closed K_ss + NGSolve FEM ground truth verified) |
| `time_domain/` | Time-domain realization: AAA rational fit → 21 stable poles (cube). No `d` tuning. |
| `ngsolve_validation/` | Framework-agnostic NGSolve FEM cross-validation (cube / cuboid Kelvin) |
| `lshape3d_ngsolve_mellin.py` | Non-tensor 3D shape (L-shape) probe — Mellin universality check, mixed Galerkin extension is future work |

Each geometry directory contains numbered scripts (`01_*.py`, `02_*.py`,
...) implementing successive refinements: baseline 1-DOF surface,
Senior tower corrections, rank-N bulk sweep, etc.

## Headline results (post Phase 8b correction)

After cross-checking against properly-implemented analytic references
(`radia.maglev.mixed_galerkin.references` for cylinder/sphere; `_references/`
for remaining square/cube candidates):

| Geometry | 1-DOF (planar SIBC) | + γ_1 | + γ_2 | + γ_3 |
|---|---|---|---|---|
| Cylinder    | 0.04% wall band  | 2.4e-4% | 2e-5% | 1e-5% |
| Sphere      | 0.11%            | 0.001% | (Senior tower terminates at γ_1) | |
| 2D square   | 0.03 – 0.26%     | — (flat face: Senior trivial) | corner Mellin needed | |
| 3D cube     | 0.33% vs NGSolve FEM (rank 20, closed K_ss) | — | — | |

All results obtained with **zero free parameters**: the bulk-surface
crossover frequency is determined by the Galerkin system, not a
user-supplied `d` (cf. the Warburg-terminated CLN of the IGTE digest
where `d` is fit). For the smooth-boundary bodies (sphere, cylinder),
the Senior tower coefficients γ_k = -H_mean, (K_gauss - H_mean²)/2,
... are the canonical curvature corrections (Senior 1962 / Mitzner
1967 / Yuferev-Ida 2010).

## Phase history

The development went through ~11 numbered "Phases" within the
2026-05-28 → 2026-06-12 research sprint. The condensed lineage:

| Phase | Discovery |
|---|---|
| 1   | Reproduce digest's CLN3 + Warburg with `d` tuned (17% wall band). |
| 2   | Replace Warburg block with s-dependent SIBC envelope → 0.04% wall band, no `d` needed. |
| 3   | Bulk rank-N sweep — does NOT help (envelope-limited). |
| 4   | Extend to sphere — 1-DOF gives 0.11%. |
| 5   | Extend to 2D square — naive v1 fails (89%), corner-aware v2 envelope ψ = f(x)f(y) gives apparent 0.99%. |
| 6   | Extend to 3D cube — corner envelope ψ = f(x)f(y)f(z), apparent 20%. |
| 7   | Sphere HOIBC γ_1 = -1/a curvature correction → 100× improvement (0.11% → 0.001%). |
| 8   | Cylinder HOIBC γ_1 = -1/(2a) — appeared NOT to help at wall band. |
| **8b** | **Discovered Y_exact bug**: the cylinder reference switched from full Bessel to K_SIBC/√s asymptote at \|γa\| = 50. The "0.93% saturation" and "γ_1 ineffective" findings of Phases 2/3/8 were ARTIFACTS. True cylinder accuracy: 0.04% baseline, with Senior tower truncation giving 100× per added DOF. |
| **8c** | **Foster sum reference audit for 2D square**: N=1999 had 5% bias at wall band; the "0.99% peak" was 97% Foster truncation. True 2D square accuracy: 0.03–0.26% baseline. |

The Phase 8b/8c corrections are the reason this directory exists as a
clean restart — the scripts archived in
`../../digest/supplement/2026_05_28_*.py` reflect the original (now
known-buggy) versions and are kept frozen for publication-history
purposes. See `digest/supplement/PHASE_8B_ARTIFACT_NOTE.md`.

## Open questions

- **Task #183** — 3D cube Foster reference audit. CLOSED by replacing the reference rather than converging it: `cube3d/06_ngsolve_ground_truth.py` measures against an NGSolve FEM solution, giving 0.33% at rank 20 with the closed `K_ss`. The Phase 6 "20% wall band" was indeed a Foster N=99 truncation artifact, the same failure the 2D square showed in Phase 8c. `cube3d/01_corner_envelope_uncertain.py` still carries the provisional Foster-referenced number in its docstring and is superseded by 06.
- **Two-point Padé Theorem** — formalize the bound err_intermediate ~ (δ/L)^{2N} for rank-N bulk + N-DOF Senior tower.
- **Time-domain Cauer realization** — each Senior tower correction has a fractional-power impedance signature; the diffusive Foster quantization technique (digest §IV) can realize each as a finite RC ladder. This is the Paper 2 direction.

## Related references in this repo

- `../README.md` — top-level Hierarchical Cauer + Warburg-Schur framework
- `../frequency_domain/circle_warburg_plot.py` — IGTE digest cylinder demo
- `../frequency_domain/3d_sphere.py` — sphere Foster sum
- `../../digest/supplement/` — Phase 2/3 archive (2026-05-28 frozen versions, with Y_exact artifact)
- `../docs/hierarchical_cauer_sibc.md` — overall documentation

## Bibliography (selected)

- Schur, I. (1917). On power series bounded in the unit disc. *J. Reine Angew. Math.* — **Schur complement**.
- Guyan, R. J. (1965). Reduction of stiffness and mass matrices. *AIAA J.* 3(2):380 — **static condensation**.
- Kron, G. (1939). *Tensor Analysis of Networks*. Wiley — **Kron reduction**.
- Senior, T. B. A. (1962). Impedance boundary conditions for imperfectly conducting surfaces. *Appl. Sci. Res. B* 8:418 — **γ_k tower**.
- Mitzner, K. M. (1967). An integral equation approach to scattering from a body of finite conductivity. *Radio Sci.* 2:1459 — **curvature correction**.
- Quarteroni, A. & Valli, A. (1999). *Domain Decomposition Methods for Partial Differential Equations*. Oxford — **Steklov-Poincaré operator**.
- Yuferev, S. V. & Ida, N. (2010). *Surface Impedance Boundary Conditions: A Comprehensive Approach*. CRC Press — **Rytov-form HOIBC**.
- Kameari, A. *et al.* (2018). Cauer ladder network representation of eddy-current fields for MOR using FEM. *IEEE TMag* 54(11):7202804 — **CLN bulk basis**.
