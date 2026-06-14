# Track A: performance of the Kelvin transformation, via the DtN spectrum (the SA / Hachinohe paper)

> **SCOPE (2026-06-15): TRACK A — the core paper of this session.** Evaluate/quantify the *performance*
> of the Kelvin-transformation open boundary using the exterior Dirichlet-to-Neumann (DtN) spectrum as
> the metric. Track B (stream-function coil design with iron) is a SEPARATE paper, handed off on the
> `streamfunction` branch (`HANDOFF_sommerfeld_dtn_kelvin_streamfunction.md`). This document is Track A.

## Thesis (one line)
The exterior **DtN (Steklov-Poincaré) spectrum** is the natural performance metric for the Kelvin open
boundary: the closed-form ladder `λ_n = −(n+1)/R` (3D; `−n/R` in 2D) is the exact target, and **how high
in degree `n` a discretization climbs it before "peel-off" quantifies the accuracy**. The order `p` sets
the multipole reach (`n ≈ p`), so the spectrum is a **predictive datasheet** that fixes the required
element order / surface resolution BEFORE solving — and puts Kelvin, BEM, PML and infinite elements on a
single comparison axis.

## What "performance" means here (the axes the paper measures)
1. **Accuracy** — per-degree DtN error vs `n`; the peel-off knee moves right with order `p` (controllability).
2. **Cost** — sparse SPD volume (~tens nnz/row, local grad-grad, **no Green function, no singular
   quadrature**) vs dense BEM (`nnz = DoF²`, `G` + singular quadrature). The "**DoF ↑ but cost ↓**" result.
3. **Coarseness** — how coarse the truncation mesh may be: the spectrum says **order, not h** (a single
   coarse curved shell at `p ≥ n` beats graded multi-layer/h-refinement and PML for a decaying field).
4. **Scope/robustness** — it is the *exact* exterior+open DtN (∞ baked in), and generalises to arbitrary
   exterior **material** (demo_t) and arbitrary **body geometry** (demo_w/bb) — the Sommerfeld surrogate.

## Evidence chain (Track A, all committed + verified)
| theme | demos |
|---|---|
| the ladder / datasheet | `demo_d`, `demo_m` (multipole ceiling), `floor_vs_curve`, `p_vs_h_study`, `sufficient_mesh` |
| coarse-mesh / minimal ball | `demo_l`, `kelvin_exterior_mesh`(3), `demo_g` (closure hierarchy: Dirichlet/Robin/Kelvin) |
| Kelvin vs BEM (cost & equivalence) | `demo_k`, `demo_n`, `demo_o` (vs H-matrix), `demo_r` (Schur = dense kernel), `demo_s` (refined-limit equivalence) |
| the matrix & what it IS | `demo_v` (assemble material-aware DtN matrix), `demo_w`/`demo_bb` (arbitrary body / non-layered), `demo_cc` (it is condensed FEM, not BEM — the Green-function line), `demo_dd` (WHEN to form it) |
| Sommerfeld performance/scope | `demo_x` (static isomorphism), `demo_y`/`demo_z`/`demo_aa` (multilayer / all-frequency / low-freq eddy) |
| certified quantities | `inductance_dtn`, complementary bracket `demo2`, A-formulation `demo3` |

## The performance story (comparisons to make explicit)
- **vs BEM / H-matrix**: same exterior DtN operator (method-agnostic), but sparse-SPD-no-G vs dense-G;
  measured fill/time gap (`demo_k`/`demo_n`/`demo_r`); H-matrix compresses the dense block, Kelvin has
  structural zeros (`demo_o`). The arbitrary-Γ parity is met by condensing on the body surface (`demo_w/bb`).
- **vs PML / infinite elements**: PML is for radiating (wave) problems; a magnetostatic field decays
  algebraically and the Kelvin inversion is the *exact* infinite element — one coarse `p ≥ n` shell beats
  graded layers (knowledge `dtn_coarse_mesh`, section 3).
- **vs analytic**: the ladder `(n+1)/R` and its material/multilayer shifts are the closed-form references.

## Honest novelty (state precisely; the discipline that held all along)
"Lighten open-boundary BEM with transformed/Kelvin FE, DoF↑/cost↓" is **Remacle ~1995 / Lowther 1989**
(free-space) — CITE, do NOT claim. The genuinely-new contributions of THIS paper:
1. the **DtN-spectral performance reinterpretation** — recasting the empirical "it is faster" as a
   closed-form `−(n+1)/R` datasheet (`p` = multipole reach) that *predicts* accuracy/cost before solving;
2. the **material / Sommerfeld generalisation** of the open-boundary operator (exterior material, arbitrary
   body, the Sommerfeld isomorphism/surrogate) on top of the authors' own Sugahara 2013/2017/2022/2025.

## Paper plan (`W:\02_学会資料\...\Kelvin_DtN@...\原稿\kelvin_dtn.tex`)
1. Open-boundary problem + the DtN/Steklov-Poincaré operator as the unifying object.
2. Kelvin-FEM formulation; the closed-form ladder; `p` = multipole reach (the datasheet).
3. Performance: accuracy (peel-off), cost (sparse vs dense, DoF↑/cost↓), coarse-mesh-suffices.
4. Comparisons: BEM/H-matrix, PML/infinite elements (one-axis via the DtN spectrum).
5. Scope: exterior material + arbitrary body + Sommerfeld isomorphism (the differentiator).
6. Honest related work (Remacle/Lowther free-space; SBFEM; Demarcke/Knockaert; Sugahara own chain).
7. むすび. (Figures via `radia_mcp.figure`; the O_h split (demo_w) and the Sommerfeld isomorphism
   (demo_x) are strong new-figure candidates alongside the existing `fig_dtn_overlay`/`fig_exterior_material`.)
