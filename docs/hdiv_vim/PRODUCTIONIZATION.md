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
| Linear solve | **`SolveLinearMaterial`: Jacobi-PCG for ((1/chi)M_mass + B^T G B) in C++ (M3, golden-locked vs scipy MINRES + dense)** | scipy MINRES / CG |
| Nonlinear demag | **scalar-chi Picard `SolveNonlinearPicard` in C++ (M3): isotropic nonlinear demag M=Mof(H0−Dscal·M), golden vs Python Picard <1e-5 + analytic fixed point <1%**; per-element tensor-tangent Newton (non-uniform M) still NGSolve | `solve_nonlinear_newton` (dense + scalable; scalable uses the analytic C++ Gram, M2b) |
| Curved + high-order Gram | — (uses NGSolve) | `ngsolve.bem` single-layer (sphere/spheroid/ellipsoid validation) |
| Symmetry image method | — | `hdiv_demag_symmetry_image.py` (sphere 1/2,1/4,1/8 validation, crude Gram) |
| **Public Radia API** | analytic `_ChargeGramHMatrix` + `SolveLinearMaterial` (scalable demag + linear solve) | `radia.hdiv_vim` package (build_demag / solve_nonlinear_newton[_scalable]) |

Progress (2026-06-08): **M2 DONE** (the accurate analytic charge Gram is in the C++ scalable path,
golden-locked; the scalable nonlinear Newton rewired onto it). **M3 DONE** (the warranted C++ work:
SolveLinearMaterial Jacobi-PCG + SolveNonlinearPicard scalar-χ nonlinear; the per-element tensor-tangent
Newton is **measured-not-warranted** — the production Newton converges in 5–6 iters so a C++ port buys
little). The remaining gap to seal yano-type: the M4 production pieces + the broader **speed** parity
matrix (M0; the nonlinear-Newton case is now measured — 5–6 iters, orchestration negligible).

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

**The comparison headline — LOOPS (not A_SS, not H-LDLᵀ).** Frame the HDiv-VIM advantage at the LOOP
level: HDiv-VIM's loops are ker(B), FIELD-NULL BY CONSTRUCTION (de Rham / Piola) — no loop-star, no
cohomology, μr-INDEPENDENT iterations. In the conventional MMM/MSC path the loops are NOT field-null,
so they are the OBSTRUCTION: they need explicit loop-star handling, which on the solver side becomes
the A_SS "star block" that must be H-LU/H-ILU-preconditioned. (Keep the layers distinct: yano-type is
the distortion-ELEMENT formulation; A_SS is the loop-star SOLVER block — do not conflate them.) The
payoff to MEASURE (M0): because the HDiv-VIM tangent stays well-conditioned, Newton-Raphson should
converge in far fewer iterations than the conventional path whose loops can ill-condition / stall the
tangent — the head-to-head Newton-iteration count is the headline parity number. (Precise: the
HDiv-VIM *Newton* is the fast part, <30 iters at saturation; the 0.5-damped scalar Picard warmstart is
the slow-but-robust part, ~230 iters at deep saturation — do not confuse the two.)

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
- **M3 — nonlinear solve in C++. DONE (2026-06-08): the WARRANTED C++ work is complete.** (linear-solve
  kernel + scalar-χ Picard in C++; the per-element tensor-tangent Newton is measured-not-warranted —
  the Python Newton is the production path at 5–6 iters, see below.)
  - **DONE — the iterative-solve hot kernel in C++ (commit 029236d8).** `RadHACApKChargeGram::SolveLinearMaterial`
    solves the SPD material system `((1/chi)M_mass + B^T G B) m = rhs` by Jacobi-PCG with G applied as the
    analytic H-matvec — the linear demag solve + the symmetric Picard warmstart, no Python per-iteration
    glue (pybind `solve_linear_material`; golden test_hdiv_vim_cpp_linear_solve.py: vs scipy MINRES on the
    identical operator <1e-6, vs dense analytic <5e-4).
  - **DONE — the scalar-chi nonlinear solve in C++ (commit ae81ec23).** `RadHACApKChargeGram::SolveNonlinearPicard`
    solves the isotropic nonlinear demag `M = Mof(H0 − Dscal·M)` entirely in C++: each Picard step is a
    `SolveLinearMaterial` solve + the closed-form `chi_sec` update (no NGSolve per iteration). Golden
    test_hdiv_vim_cpp_nonlinear.py: == Python scalar-chi Picard on the identical operator <1e-5, == the
    analytic scalar fixed point <1% (near-saturation sphere). The nonlinear PHYSICS for an isotropic body
    is now in C++.
  - **NOT WARRANTED (M0-measured 2026-06-08) — the per-element tensor-tangent Newton in C++.** The
    tensor-tangent Newton-Raphson is the PRODUCTION nonlinear solver (`solve_nonlinear_newton`, Python +
    NGSolve assembly) and it converges in **5–6 iterations** on sphere AND non-uniform cube at genuine
    saturation (H0 = 3e5–8e5), matching the analytic fixed point to 0–0.8% and Radia MMM/MSC to <0.05%
    (test_hdiv_vim_newton_vs_radia.py). Because the Newton does only 5–6 iterations, the Python
    orchestration overhead is NEGLIGIBLE (each iteration's cost is the C++ H-matvec, already C++) — so
    porting the Newton loop to C++ buys almost nothing. (Contrast the scalar-χ Picard ~230 iters, where
    C++ DID help — which is why the Picard, not the Newton, was C++-ported in M3.) So the full C++ Newton
    is NOT built: the algorithm is the win and it is already fast. The per-element constitutive assembly
    (b_M, T) stays in NGSolve — the standard FE stack, not a limitation. (This re-scoping is the M0-type
    measurement for the nonlinear case: 5–6 iters ⇒ C++ port not justified.)
  > **H-LDLᵀ DELETED — 消去 done 2026-06-08.** The direct symmetric H-LDLᵀ factorization was removed from
  > the codebase entirely: the HDiv-VIM is μr-INDEPENDENT by construction, so the iterative solve
  > (MINRES linear / GMRES/Picard nonlinear) is cheap and well-conditioned — a direct factorization
  > bought little and was unused (only self-tests + HDiv-VIM research called it). Removed: the
  > `Symmetric H-LDL^T` section of `cHACApK_harith.{c,h}` (1188 + 119 lines), the `_hldlt_self_test*` +
  > `factor_solve_hldlt` pybind, `tests/feec/test_hldlt_*.py` (feec 94 → 81), the hldlt examples.
  > **H-LU / H-ILU were NOT touched** — they are the SAME `cHACApK_hlu_*` subsystem (the rk-truncation
  > tol switches accurate-H-LU ↔ incomplete-H-ILU), load-bearing as the MMM/MSC solver's A_SS
  > preconditioner (a separate solver layer from the yano-type element formulation).
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
