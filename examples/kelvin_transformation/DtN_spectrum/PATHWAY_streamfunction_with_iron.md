# New pathway: stream-function coil design with magnetic material via a Kelvin-FEM material-aware DtN matrix

> **SCOPE (2026-06-15): this is TRACK B — a SEPARATE paper from the DtN+Kelvin core.**
> Decision: keep two distinct tracks so the core is not diluted by the application.
> - **Track A — DtN + Kelvin (core; the SA / Hachinohe paper):** the DtN-spectrum datasheet, the sparse
>   Kelvin open-boundary, the Sommerfeld isomorphism/surrogate, the directly-assembled material-aware DtN
>   matrix and what it IS (FEM-condensed, not BEM). Demos d…dd + x/y/z/aa/cc.
> - **Track B — stream-function coil design with iron (this document; a separate paper):** uses the
>   Track-A operator as the *material-aware design kernel*. Demos ee/ff + a future general-iron design.
> The two share machinery but are written up independently.

*(Consolidation of the stream-function track, 2026-06-15. Honest novelty status at the bottom — a targeted check is running; phrase any claim as "to our knowledge".)*

## The idea (one line)
The **transfer / DtN matrix** that a stream-function (surface-current / current-potential) coil design
inverts — the linear map `psi -> field` — is, **with magnetic material present**, the system's
*material Green operator*. Generate it **sparsely and Green-function-free** as the Schur complement of a
**Kelvin-transformed FEM** (which carries arbitrary `mu(x)` in the inverted exterior). Then coil design
is the same clean linear inverse as in free space, but with the correct **material-aware** kernel.

## Why it was hard (the user's observation: "流れ関数法は磁性体があると楽じゃない")
- Free space: `psi -> field` kernel = **Biot-Savart** (analytic, easy). Design = a clean linear solve.
- With iron (yoke / shield / core): total field = coil field **+ iron reaction**; the kernel becomes the
  **material Green operator**. For planar/cylindrical iron that is the (hard) layered/**Sommerfeld** Green
  function; for **arbitrary** iron there is **no closed-form Green function at all** (a volume integral
  equation revives the dense volume unknown). So the clean "psi x kernel" structure is lost.

## The mechanism (what is new)
1. Mesh **iron + coil surface + open exterior** once; the Kelvin inversion makes the unbounded exterior a
   bounded sparse SPD volume (no Green function, infinity exact).
2. **Condense (Schur)** onto {coil surface, target} -> a small **dense, material-aware transfer/DtN matrix
   M** (the operator BEM would need a Green function to build).
3. **Design = invert M**: `psi = M^+ B_target`. The matrix is the deliverable (this is precisely the
   "operator is the deliverable" case — you do NOT just solve one field; the inverse design consumes M).

## Evidence chain (all committed, verified)
| demo | establishes |
|---|---|
| `demo_t` | FEM-Kelvin carries arbitrary exterior material (layered shell matches analytic ~1e-4) |
| `demo_v` | the material-loaded exterior DtN/Green **matrix** is assembled directly (Schur), spectrum = analytic |
| `demo_w`,`demo_bb` | the matrix for **arbitrary geometry** (cube O_h split; non-layered on-axis inclusion C∞v split) — where **no Sommerfeld Green function exists** |
| `demo_x`,`demo_y`,`demo_z`,`demo_aa` | Kelvin-FEM is the Sommerfeld operator (static isomorphism; multilayer kernel; works DC->wave; low-freq eddy-current = Bannister complex image) |
| `demo_cc` | it is still **FEM** (condensed substructure / SBFEM), not BEM — the Green-function criterion |
| `demo_dd` | **when** to form the matrix: only when the operator is the deliverable (not to solve one field) |
| `demo_ee` | coil + iron shield: free-space kernel off by up to **~16x**; material-aware matches ~1e-4 |
| `demo_ff` | **design inverts M** (concentric, MODAL): with the material-aware M the target is hit (2e-16); the free-space-designed coil misses by **77%** in the iron system |
| `demo_hh` | **the GENERAL next step** (this track's deliverable): a REAL winding-surface stream function `psi` (the order-p H1 trace, 378 nodal DoFs -- not 3 modes) + ARBITRARY (NON-concentric) iron blob, where no Sommerfeld Green function exists. M[target, psi-dof] assembled directly (one factorisation + a back-sub per coil DoF). Concentric sub-case ANCHORS M to the analytic layered transfer (rel 3e-3..2e-2; M@psi == a fresh FEM solve to 1e-15). Non-concentric: the iron-aware design HITS the target in a fresh full Kelvin-FEM solve (1e-14), the free-space design MISSES by ~43% (stable 31-39% across mesh refinement -- a real shield effect, not a coarse artifact) |
| `demo_jj` | **can ACA-TSVD be applied to M?** YES, and the two ingredients answer two questions. **TSVD** -> the GLOBAL inverse design (`psi=M^+ B`), unchanged + necessary (M = compact forward map, SVs decay, measured cond(1e-6)~9.5e3). **ACA** -> BLOCK-WISE as an H-matrix, NOT a global low-rank factor: global M is near-full-rank (189/195, near field) but the material Green kernel `G_mu(x_t,y_c)` is asymptotically smooth for well-separated clusters even with iron between them, so an ADMISSIBLE block (compact source cap + far targets) is low rank (verified 6/20, sigma 1->4e-3; rSVD 10 oracle-calls -> 5e-6) while a NEAR block stays dense (18/20) -- the **same near/far split HACApK already runs** on the MMM/MSC matrix. Only the entry oracle changes: 1 column = 1 Kelvin-FEM back-substitution (a fast matvec of M), 1 row = 1 adjoint solve -> a **randomized SVD/Lanczos is even more natural than entrywise ACA**. Answers the "wire through `radia.streamfunction` ACA-TSVD" open item: yes, with the Kelvin-FEM-solve oracle |
| `demo_kk` | **the wiring made real: DtN M -> the EXISTING `radia.stream_function` (ACA+)+TSVD + ridge solver.** The cheap DtN material-aware M is handed to the kernel-agnostic `aca_tsvd(M_rows, N_cols, entry)` (entry = array lookup -- the cheap-build path the benchmark justifies; ACA+ reproduces the formed M to **1e-14**), then designed by `RegularizedTSVD.from_stiffness(res, S)` with `S` = the coil **surface current-density seminorm** (tangential trace gradient; `psi^T S psi = ||K||^2`). Verified (non-concentric blob, mu_r=50, order 2, 378 coil DoFs, 24 targets): TSVD (min-L2) and the **ridge (min current density)** both HIT the iron-system target in a fresh Kelvin-FEM solve, the ridge with **23% less current density** (||K|| 3.37 vs 4.39); an alpha L-curve trades field accuracy (9.7e-10->4.8e-1) for ||K|| (3.37->1.69) re-solving only the k x k core; the **free-space-designed coil MISSES by 22.5%** through the same solver. **No DtN-specific design code** -- only the `entry` callback changes |

## To turn into a contribution (next steps)
1. **General (arbitrary iron) coil-design demo** — **DONE: `demo_hh_general_iron_design.py`** (2026-06,
   on the `streamfunction` branch). Non-concentric iron blob -> a genuinely coupled M; a REAL surface-
   current `psi` (the order-p H1 nodal trace on the coil, 378 DoFs, not modal amplitudes); the designed
   coil forward-verified in a fresh full Kelvin-FEM solve (hits 1e-14) while the free-space design MISSES
   by ~43% (stable across mesh). Concentric sub-case anchors M to the analytic transfer (rel 3e-3..2e-2).
   *Open refinements:* a non-spherical coil former (cylinder), a m≠0 tesseral target, and wiring the
   inverse through `radia.streamfunction`'s ACA-TSVD rather than a dense numpy TSVD. **The
   ACA-TSVD applicability is now settled by `demo_jj_aca_tsvd_on_dtn_matrix.py`** (TSVD applies to the
   global inverse unchanged; ACA applies block-wise as an H-matrix — admissible block rank 6/20, near
   block 18/20, global 189/195 — the same near/far machinery as HACApK, oracle = one Kelvin-FEM solve
   per column). **The wiring itself is now DONE: `demo_kk_streamfunction_ridge_with_dtn.py`** hands the
   formed DtN M to the existing `radia.stream_function.aca_tsvd` (entry = array lookup; ACA+ reproduces M
   to 1e-14) and designs with `RegularizedTSVD` (the ridge / current-density-min solver) -- both TSVD and
   ridge hit the iron target, the ridge with 23% less current density; the free-space kernel through the
   same solver misses by 22.5%. So the cylinder former and m≠0 tesseral target are the only remaining
   refinements of step 1.
2. **Benchmark** M-build (sparse Kelvin-FEM Schur) vs the dense layered-Green / FE-BEM baseline —
   **DONE: `bench_dtn_mbuild.py`** (JSON + figure committed). Three measured contrasts: **(C1
   sparsity/scaling)** the Kelvin-FEM volume matrix has CONSTANT nnz/row (~15, fill→7e-4) so its
   storage is linear in ndof, while the dense operator a Green/BEM route forms is O(ndof²) — **941×
   larger at 25k DoF** (3388 MB vs 3.6 MB), sparse factor 0.6 s; **(C2 accuracy)** at order-2 the FEM
   transfer R_n reproduces the analytic layered-sphere Green transfer to **rel_max 2.5e-2** (n=1
   2.9e-3, n=2 1.6e-2, n=3 2.5e-2) — Kelvin-FEM builds the SAME operator the Green route gives;
   **(C3 generality)** a non-concentric iron blob builds at the identical sparse cost where **no
   closed-form Green function exists at all**. Selling point confirmed: sparse, material-aware, no
   Green function.
3. **Manuscript**: position as a *formulation* contribution; cite the free-space transformed-FE prior art
   (Brunotte 1992, Meeker 2013), the stream-function/target-field design lineage, and the author's own
   Sugahara 2022 (uniform specimen) as the foundation extended.

## Honest novelty status — VERDICT (targeted check, 2026-06-15)
**NOVEL, confidence 0.83** — 5 search agents + 4 adversarial refuters, ALL "not preempted". The specific
combination (Kelvin/shell-transformed FE open boundary → Schur-condensed DtN/transfer matrix carrying
finite-permeability iron in the transformed exterior → used AS the stream-function coil-design
sensitivity kernel) was not found anywhere.
- **It FUSES Sugahara's OWN two threads** (which exist only separately): his Kelvin open-boundary FEM
  (Extended Kelvin, IEICE E108-C 2024/25; ECT-with-Kelvin, IEEE TMAG 58(9) 2022 — forward solvers) vs
  his free-space stream-function coil design (Koiso/Sugahara/Ida, TSVD+ACA CEFC 2024; ACA+CMA-ES 3D coil
  IEEJ 2025). Frame the contribution as the FUSION, not as new over either half (defuses self-citation).
- Where coil design DOES carry iron elsewhere, the kernel is categorically DIFFERENT — BEM/μ→∞
  equipotential (bfieldtools, Mäkinen–Zetter 2020), image/modified-Green, saturated dipole (Landreman
  2025), linearized magnetization-response (passive shimming) — never a transformed-FE DtN material
  operator. Closest (non-preempting): Wang et al., *Measurement* 2024 (gradient coil + ferromagnetic
  shield via Green/image); bfieldtools.
- **Recommended phrasing (NOT a bare "first"):** "To the best of our knowledge, this is the first method
  to use a Kelvin/shell-transformed FE open boundary, condensed by a Schur complement into a DtN
  (transfer) matrix that carries the finite-permeability iron in the transformed exterior, directly as
  the sensitivity kernel of a stream-function coil-design inverse problem. We are not aware of any prior
  work fusing these two ingredients ..."
- **Residual checks before a 'first' claim:** (1) full texts of paywalled shield/coil papers (esp. Wang
  2024 *Measurement* S0263224124008339); (2) Japanese grey-lit (IEEJ 静止器/マグネティックス研究会,
  J-STAGE, CEFC/COMPUMAG 2022-26, in Japanese: Sugahara/Koiso/Sato/Ida + ケルビン変換 + 電流ポテンシャル);
  (3) 2025-26 preprints + the stellarator REGCOIL/current-potential line; (4) confirm no in-press Sugahara
  paper already fuses them; (5) match wording to the exact construction (total vs reduced scalar potential;
  transfer vs DtN; finite vs infinite μ); (6) accelerator-magnet (ROXIE/CERN) field-quality design re-scan.
