# HDiv-type VIM — productionization roadmap (the path to retire yano-type)

The HDiv-type VIM method is validated (see [README.md](README.md), feec suite 85/85). This document is
the roadmap to make it a **shipped production solver that retires the yano-type distortion elements** —
i.e. a drop-in MMM/MSC replacement that is at-least-as-good as the shipped solver on **every** case
yano-type handles, behind a clean Radia API. Honest scope, milestone-based, with a hard
**definition-of-done** (the parity gate).

## Current state (inventory, 2026-06-08)

| Layer | C++ (compiled, `_radia_pybind.pyd`) | Python-only (prototype) |
|---|---|---|
| Charge map B + HDiv mass | structured **hex** (`rad_hdiv_vim.cpp`) | unstructured **tet** (NGSolve HDiv extraction, `examples/feec_vim/`) |
| Coulomb Gram G | monopole + sub-point (`CoulombGramEntry`); **analytic Wilton `TriPotential` + `PhiTet` WIRED into the `_ChargeGramHMatrix` analytic entry (M2a+M2b, golden-locked == dense `analytic_gram` ~1e-9)** | Wilton surface + `phi_tet` volume (analytic), dense `build_demag` |
| Scalable Gram H-matrix | **`_ChargeGramHMatrix`** monopole **+ analytic mode (M2b)**, `_HDivVimHMatrix` (hex) | (assembly driven from Python) |
| H-LDLᵀ factorization | cHACApK (`_hldlt_self_test*`, `factor_solve_hldlt`) — **TO DELETE (消去, see M3)** | (not wired into the material-system apply) |
| Linear solve | **`SolveLinearMaterial`: Jacobi-PCG for ((1/chi)M_mass + B^T G B) in C++ (M3, golden-locked vs scipy MINRES + dense)** | scipy MINRES / CG |
| Nonlinear damped Newton | — (the per-element constitutive + nonsymmetric tangent GMRES + Newton loop still Python) | `solve_nonlinear_newton` (dense + scalable; scalable uses the analytic C++ Gram, M2b) |
| Curved + high-order Gram | — (uses NGSolve) | `ngsolve.bem` single-layer (sphere/spheroid/ellipsoid validation) |
| Symmetry image method | — | `hdiv_demag_symmetry_image.py` (sphere 1/2,1/4,1/8 validation, crude Gram) |
| **Public Radia API** | analytic `_ChargeGramHMatrix` + `SolveLinearMaterial` (scalable demag + linear solve) | `radia.hdiv_vim` package (build_demag / solve_nonlinear_newton[_scalable]) |

Progress (2026-06-08): **M2 DONE** (the accurate analytic charge Gram is in the C++ scalable path,
golden-locked; the scalable nonlinear Newton rewired onto it). **M3 partial** (the iterative linear
solve is in C++; the full nonlinear Newton-in-C++ remains). The remaining gap to seal yano-type: the
full C++ Newton (M3), the M4 production pieces, and the **speed** parity number (M0, never measured).

## Definition of done — the parity gate (M0)

yano-type is sealed **only when** HDiv-VIM is at-least-as-good on the full case matrix the shipped MSC
handles, measured head-to-head:

| Case | Accuracy criterion | Speed criterion |
|---|---|---|
| Linear soft iron (mu_r 10–1e5), convex + non-convex | match shipped MSC to its mesh-converged value | wall-clock within ~2x of MSC (target: faster) |
| Nonlinear BH iron (cube, C-yoke, real table) | match MSC volume-avg `<1%` (done at prototype level) | within ~2x |
| Permanent magnet + soft iron (mixed) | match MSC | within ~2x |
| IMA symmetry 1/4, 1/8 | match the full model | fewer DOF → faster |
| Distorted meshes (high distortion) | mu_r-independent (done) + correct values | bounded iters (done) |
| Standard validation set (the lab's MSC golden problems) | parity | parity |

**M0 deliverable**: (1) enumerate the exact golden problems above into a `tests/feec/parity_vs_msc/`
suite; (2) **measure the speed gap honestly** — isolate the *algorithm* from the Python orchestration:
time the C++ `_ChargeGramHMatrix` demag MatVec + the iteration count vs the shipped MSC interaction
MatVec at matched N (the orchestration becomes C++ in production, so the per-MatVec + iters are the
production-representative numbers). This sizes the C++ lift and fixes the "done" bar.

## Milestones

- **M0 — parity gate + speed-gap measurement** *(START HERE; mostly measurement, low risk).* The
  definition-of-done above + the honest speed number. Until M0, "retire yano-type" is conference-ready
  (validated method) but the production gap is unquantified.
- **M1 — production module + public Radia API.** Move the validated solve out of `examples/feec_vim/`
  into `src/radia/` with a clean entry (e.g. `rad.hdiv_demag_solve(mesh, materials, source)`), driving
  the existing C++ `_ChargeGramHMatrix` + the Newton. Golden-test against the examples' validated
  numbers. Makes HDiv-VIM a usable Radia feature (first shippable step).
- **M2 — DONE (2026-06-08): accurate Gram in the C++ scalable path.** `RadHACApKChargeGram` gained an
  ANALYTIC mode so its entry is the exact Wilton/`phi_tet` charge Gram (was pure centroid-monopole),
  and the scalable nonlinear Newton was rewired onto it. Both sub-steps landed + golden-locked:
  - **M2a — DONE (commit 8f25788a).** Ported `tri_potential` (Wilton `log`/`atan2` exact triangle 1/r
    surface potential) + `phi_tet` (divergence-theorem tet Newtonian volume potential) to C++ as
    `rad_hdiv::TriPotential` / `PhiTet`, exposed as `_hdiv_tri_potential` / `_hdiv_phi_tet` and
    **golden-locked vs Python to ~machine precision** (test_hdiv_vim_cpp_potentials.py: 7.6e-16 / 2.3e-15).
  - **M2b — DONE (commits c31e1501 core + 6e8aa494 nonlinear).** Threaded per-charge vertices + type into
    the `RadHACApKChargeGram` analytic ctor + `_ChargeGramHMatrix(cell_verts,face_verts,n_el)` pybind; the
    analytic entry `G[a][b] = 0.5(QuadDot(a,b)+QuadDot(b,a))` (PhiTet/TriPotential inner × tet-subpoint /
    Dunavant-5 outer) matches dense `build_demag(analytic_gram=True)` entry-by-entry (all-dense <1e-9;
    test_hdiv_vim_chargegram_analytic.py). The scalable nonlinear Newton now applies `N v = B^T(H.matvec(B v))`
    with the analytic H-matrix (no Python near-correction) → reproduces the dense ANALYTIC Newton on the
    NON-uniform cube (test_hdiv_vim_newton_scalable.py).
  GOTCHA (learned M2a/M2b, keep): on the NAS-mounted `S:` source, `Build.ps1`'s freshness check
  (`src.LastWriteTime -le dst.LastWriteTime`) can skip the `.pyd` copy AND a fresh `import radia` can
  briefly read the pre-copy `.pyd` (NAS read-cache lag). Workaround: `Copy-Item -Force
  build-msvc\_radia_pybind.cp312-win_amd64.pyd src\radia\_radia_pybind.pyd` then verify in a NEW process.
- **M3 — nonlinear Newton in C++. PARTIAL (2026-06-08).**
  - **DONE — the iterative-solve hot kernel in C++ (commit 029236d8).** `RadHACApKChargeGram::SolveLinearMaterial`
    solves the SPD material system `((1/chi)M_mass + B^T G B) m = rhs` by Jacobi-PCG with G applied as the
    analytic H-matvec — the linear demag solve + the symmetric Picard warmstart, no Python per-iteration
    glue (pybind `solve_linear_material`; golden test_hdiv_vim_cpp_linear_solve.py: vs scipy MINRES on the
    identical operator <1e-6, vs dense analytic <5e-4).
  - **REMAINING — the full nonlinear Newton-in-C++.** The per-element tensor-tangent constitutive
    (M(H) + dM/dH from the B-H curve) + the RT0 field reconstruction + the NONSYMMETRIC tangent GMRES +
    the damped Newton loop, all in C++ (reusing the `SolveLinearMaterial` Krylov + the analytic H-matvec).
    This is the main speed lever; it is a substantial port and its value is gated by the M0 speed number,
    so it is the right next deliberate (clean-context) C++ task — not a session-tail rush.
  > **H-LDLᵀ TO BE DELETED — 消去 confirmed 2026-06-08.** The direct symmetric H-LDLᵀ factorization of
  > the H-matrix is **not** on the production path: the HDiv-VIM is μr-INDEPENDENT by construction, so
  > the iterative solve (MINRES linear / GMRES nonlinear) is cheap and well-conditioned — a direct
  > factorization buys little and adds a C++ component to mature. Decision: **remove it from the
  > codebase entirely** (not just off-path). This is its own focused cleanup commit with a rebuild +
  > suite (feec drops ~13 tests, 88 → 75), best done with clean context (it touches an external-lib
  > header). Scope (8 code files + 3 doc files):
  > - C++: the `Symmetric H-LDL^T factorization` section of `src/ext/HACApK/cHACApK_harith.h`
  >   (`cHACApK_hldlt_decomp` / `_solve_vec` / `_factor_leafmtxp` / `_apply` / `_free_factors` /
  >   `_get_storage` + self-test) **and its `.cpp` implementations**; `factor_solve_hldlt` in
  >   `src/core/rad_hacapk_hdiv.{h,cpp}`; the `HLDLTSelfTest*` helpers + `_hldlt_self_test` / `_rk`
  >   registrations in `src/lib/radia_pybind.cpp`.
  > - tests: `tests/feec/test_hldlt_factorization.py`, `test_hldlt_real_gram.py`.
  > - examples: `examples/feec_vim/hldlt_block_flops.py`; strip the H-LDLᵀ paths from
  >   `hdiv_demag_speedup.py` + `hdiv_demag_twolevel.py`.
  > - docs/knowledge: this note → "deleted"; `radia_mcp …/knowledge/hdiv_vim.py` + `server.py`.
- **M4 — curved/high-order + symmetry production + the curved-nonlinear-volume gap.** Wire the
  `ngsolve.bem` single-layer (curved Gram) + the symmetry image method + a true reduced-DOF symmetry-BC
  solve into the API. Build the genuine method gap: the curved nonlinear **volume** charge (`phi_tet` on
  curved cells; `ngsolve.bem` is boundary-only).
- **M5 — the seal.** HDiv-VIM passes the M0 parity gate on the full matrix → flip the production
  default, deprecate + seal yano-type (preserve in the private ELF repo, per the original plan).

## Honest unknowns / risks

- **Speed** is the one unmeasured quantity and the gating risk. At lowest order HDiv-VIM ≈ shipped MSC
  algorithmically (same demag spectrum); the raw-speed *win* needs high-order. The *loss* risk is the
  H-matrix constant factors + the Newton orchestration. M0 measures it before committing M2/M3 effort.
- **Curved nonlinear volume charge** has no clean analytic reference (Radia facets curved geometry), so
  M4's curved-nonlinear validation is mesh-convergence / self-consistency, not vs-truth.
- The Wilton/`phi_tet` Gram may be adequate to keep **NGSolve-side / Python** (the loops are not the hot
  path); M2's C++ port is justified only if the speed-gap (M0) demands it.

## Sequencing note

M0 → M1 are the immediate, low-risk, high-value steps (measure the gap; ship a usable API). M2/M3 (the
heavy C++) are justified *by M0's speed number*. Do M0 first — it converts "retire yano-type" from a
goal into a quantified, executable plan.
