# HDiv-type VIM — productionization roadmap (the path to retire yano-type)

The HDiv-type VIM method is validated (see [README.md](README.md), feec suite 85/85). This document is
the roadmap to make it a **shipped production solver that retires the yano-type distortion elements** —
i.e. a drop-in MMM/MSC replacement that is at-least-as-good as the shipped solver on **every** case
yano-type handles, behind a clean Radia API. Honest scope, milestone-based, with a hard
**definition-of-done** (the parity gate).

## Current state (inventory, 2026-06-08)

| Layer | C++ (compiled, `_radia_pybind.pyd`) | Python-only (prototype) |
|---|---|---|
| Charge map B + HDiv mass | structured **hex** (`rad_hdiv_vim.cpp`) | unstructured **tet** (NGSolve HDiv extraction, `examples/vim/`) |
| Coulomb Gram G | monopole + sub-point (`CoulombGramEntry`); **analytic Wilton `TriPotential` + `PhiTet` WIRED into the `_ChargeGramHMatrix` analytic entry (M2a+M2b, golden-locked == dense `analytic_gram` ~1e-9)** | Wilton surface + `phi_tet` volume (analytic), dense `build_demag` |
| Scalable Gram H-matrix | **`_ChargeGramHMatrix`** monopole **+ analytic mode (M2b)**, `_HDivVimHMatrix` (hex) | (assembly driven from Python) |
| Linear solve | **`SolveLinearMaterial`: Jacobi-PCG for ((1/chi)M_mass + B^T G B) in C++ (M3, golden-locked vs scipy MINRES + dense)** | scipy MINRES / CG |
| Nonlinear demag | **scalar-chi Picard `SolveNonlinearPicard` in C++ (M3): isotropic nonlinear demag M=Mof(H0−Dscal·M), golden vs Python Picard <1e-5 + analytic fixed point <1%**; per-element tensor-tangent Newton (non-uniform M) still NGSolve | `solve_nonlinear_newton` (dense + scalable; scalable uses the analytic C++ Gram, M2b) |
| Curved + high-order Gram | — (uses NGSolve) | `ngsolve.bem` single-layer (sphere/spheroid/ellipsoid validation) |
| Symmetry image method | — | `hdiv_demag_symmetry_image.py` (sphere 1/2,1/4,1/8 validation, crude Gram) |
| **Public Radia API** | analytic `_ChargeGramHMatrix` + `SolveLinearMaterial` (scalable demag + linear solve) | `radia.vim` package (build_demag / solve_nonlinear_newton[_scalable]) |

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

**HACApK scalability — DEMONSTRATED for the charge Gram (M0 partial, 2026-06-08).** The charge-Gram
1/r H-matrix (`_ChargeGramHMatrix`) genuinely compresses: on spheres n_charge 322 → 3560, the ACA
low-rank blocks grow 0 → 1780 and the H-matrix/dense memory ratio falls 1.00 → 0.37, with H-matrix
memory growing ~N^1.6 (vs dense N²) — **sub-quadratic, trending O(N log N)**. This is BETTER than the
compact MMM/MSC "materialize-fallback" caveat: the charge Gram is a cleaner far-field 1/r kernel, so
ACA works on it. Benchmark: `examples/vim/hdiv_demag_hacapk_scaling.py` (+ `.json`). Honest scope:
shown to n~3560 (build_demag's dense-G reference is O(N²) and caps N); the trend is clear + favorable;
larger-N (10k+) confirmation needs a dense-G-free charge extraction — the remaining M0 scalability item.

**Head-to-head target — the saved yano-type C-type benchmarks (M0).** `examples/c_type_electromagnet/
nonlinear/quarter/hacapk/*.json` record the SHIPPED Radia MMM/MSC (HACApK, nonlinear Newton) on the
C-type electromagnet: at **165600 DOF, 214 nonlinear iterations, t_solve = 2607 s** (= 582 s H-matrix
build [22%, once] + 1953 s linear solve over 2686 linear iters); at 18900 DOF, 174 iters, 99 s. The
HDiv-VIM Newton converges in **5–6 iters** (sphere/cube) → ~35–40× fewer nonlinear iterations, so the
solve phase should win big; the H-matrix BUILD (582 s) is the cost to beat ("h-matrix化が遅い" is the
real risk, not the matvec). The μr-independent C++ material solver for the head-to-head is
`SolveMaterialMINRES` (commit 87c6591f; iters FLAT ~120 for μr 100–1e6). REMAINING for a real wall-clock
head-to-head: set up the C-type geometry for HDiv-VIM (mesh + charges) and time build + solve against
these JSONs — the core M0 deliverable.

**Build-time measured (2026-06-08, `examples/vim/hdiv_demag_buildtime_scaling.py` + .json).** The
analytic charge-Gram H-matrix build (charges straight from a tet mesh, KELVIN-LESS — iron only, the 1/r
Gram is the open boundary, no air/Kelvin): n_charge 281→7278 → t_build 0.19→24 s, compr 1.0→0.21,
matvec 8.5 ms @ 7278 (O(N log N)). Build scales ~N^1.1–1.3 at large N → extrapolated to C-type scale
(~30k–150k charges) ~150–1200 s, the **same order as yano's 582 s** (comparable, NOT a clear win). So
**the BUILD is the bottleneck, not the solve**: the SOLVE wins big (5–6 Newton iters vs 214 → ~tens of s
vs 1953 s). Total estimate ~**2–6× faster, build-limited**. The build cost is the ALWAYS-analytic entry
(every pair pays PhiTet/TriPotential); the **lever is a near/far split** (cheap monopole for ACA-far
pairs, analytic only for near). (KELVIN-LESS throughout: the
HDiv-VIM is a volume integral method like MMM/MSC; only the iron is meshed — Kelvin is a FEM-only need.)

**Near/far split DONE + measured (2026-06-09, `_ChargeGramHMatrix(..., near_factor)` + golden
`tests/feec/test_hdiv_vim_chargegram_nearfar.py` + demo `C:/temp/demo_nearfar_speedup.py`).** The
analytic entry `0.5(QuadDot(a,b)+QuadDot(b,a))` (PhiTet/TriPotential) is now used ONLY for NEAR pairs
(`|c_a−c_b| ≤ near_factor·(size_a+size_b)`, size = `cbrt(vol)`/`sqrt(area)`); FAR pairs use the cheap
centroid monopole `meas_a·meas_b/(4π r)`. `near_factor` default **1e30 = all-analytic** (preserves the
M2b golden); `near_factor=2.0` = the split. Measured build speedup on the non-uniform cube:
**n_charge 901→3.6×, 1541→5.2×, 4225→5.0×, 7278→6.7×** (20.95 s → 3.13 s @ 7278), with **identical
compression** (0.206 vs 0.206) — the H-matrix STRUCTURE is unchanged, the speedup is purely the cheaper
far-pair entry. Accuracy preserved: the split demag (near_factor=2.0) matches dense FULL-analytic within
**3 % on BOTH the uniform sphere AND the non-uniform cube** (the cube being where the old surface-only
near-correction failed — so near_factor=2.0 does NOT revert to that). Extrapolated: at C-type scale the
~150–1200 s all-analytic build drops to ~**25–250 s**, now **below yano's 582 s** → HDiv-VIM becomes a
clear win on BOTH build AND solve. feec 84/84.

**Scalable path is now dense-N²-FREE (2026-06-09, `build_demag(..., skip_dense_gram=True)` returns sparse
M_mass/B; golden `tests/feec/test_hdiv_vim_build_sparse.py`).** A hidden non-scalability remained: the
"scalable" Newton (`solve_nonlinear_newton_scalable`) called `build_demag(mesh)` with DEFAULT flags,
which built the dense O(N²) Gram G, the dense demag N, AND the O(N²) loop SVD, then ignored almost all
of it (the demag apply is driven by the analytic C++ charge-Gram H-matvec). At C-type scale a 150k-charge
dense `M_mass`/`N` is ~180 GB → OOM. Fix: (1) `build_demag`'s `skip_dense_gram` branch now returns
`M_mass`, `B`, `B_csr` as scipy SPARSE (the RT0 HDiv mass + the L2/SurfaceL2 charge map are LOCAL → sparse)
and `N`/`G`/`loops`/`n_loop` as `None` — no dense N² object is ever formed; (2) `solve_nonlinear_newton_scalable`
now requests `skip_dense_gram=True` and its two M_mass uses are sparse-safe (`mu @ (M_mass @ mu)`,
`M_mass.diagonal()`). The C++ solvers already took B as CSR + M_mass as COO, so no C++ change was needed.
Bit-for-bit identical to the dense reference (sparse M_mass == dense, Rayleigh denom to machine precision);
the dense REFERENCE path (small-N `demag_factor` + dense Newton) is unchanged. feec 88/88. With this, the
scalable Newton's BUILD is genuinely O(N log N) analytic-Gram + sparse FE assembly — no dense N² anywhere.

**HONEST CORRECTION (2026-06-09): the scalable nonlinear Newton is NOT yet mesh-robust — the "5–6 iters →
clear win on SOLVE" above held only at COARSE mesh.** The first real C-yoke wall-clock head-to-head
(`examples/vim/hdiv_cyoke_headtohead.py`) measured the SOLVE degrading sharply with refinement:
iters 6 (h=0.008) → 27 (h=0.006) → 37 (h=0.005), with Mz appearing to "drift" 589k → 509k. A full
instrumented diagnosis (per-iter ‖F‖/λ/Mavg trajectory) found:
  1. **The method + tangent are CORRECT** — once in the basin the Newton converges QUADRATICALLY
     (relF 0.035 → 6e-8, λ=1.0) to the right Mz=587566 (matches the mesh trend); verified at h=0.006.
  2. **The "drift" was a FALSE-CONVERGENCE BUG** — the old break test `|M_now-M_prev|<1e-8` (Mavg
     stagnation) fires during the slow globalization phase while relF is still ~0.1, silently returning
     an under-converged M. **Fix A (committed): break on relF=‖F‖/‖M_mass m‖ < newton_tol, and RAISE
     if maxit is hit (fail loud, CLAUDE.md "No Fallbacks") — never silently return.** golden
     `tests/feec/test_hdiv_vim_newton_scalable.py::test_scalable_newton_fails_loud_*`.
  3. **The slowness (27 vs 5–6 iters) is the INEFFECTIVE WARMSTART — and the root cause is a 2.5–5%
     ASYMMETRIC charge-Gram H-matrix, NOT the loop conditioning.** The warmstart's MINRES on the +N
     system hit its 2000-iter cap (`info=2000`). The diagnosis went deeper than "loop near-null modes":
     - **B1 (−N material formulation) is RULED OUT.** The −N MINRES (`(1/chi)M_mass − N`, the golden
       `test_hdiv_vim_solve.py` mu_r-independent solve) is **non-physical** — measured: for the sphere it
       gives `m_avg` of the WRONG sign/value (−4.503 at mu_r=10 vs the physical +2.250). The **physical**
       magnetization system is **+N** (`((1/chi)M_mass + N)m = M_mass h_ext`, matches the analytic
       demag-limited M exactly); only its CG count grows with mu_r (36→133). You cannot substitute the
       −N solver — different operator. (`C:/temp/check_pmN.py`.)
     - **With the FULL `M_mass⁻¹` preconditioner, +N CG is already mu_r-bounded at coarse mesh (37 iters
       flat to mu_r=1e6)** — so the loop conditioning is NOT the blocker (`C:/temp/proto_B2.py`).
     - **The real blocker: the ACA charge-Gram `Hg` is 2.5–5% ASYMMETRIC at scale** (n_charge=6812:
       `vᵀNw` vs `wᵀNv` rel 5.2e-2). The dense analytic entry `0.5(QuadDot(a,b)+QuadDot(b,a))` is
       symmetric, and the M2b golden checks Hg==dense `<1e-9` — but only at COARSE mesh; the ACA
       symmetric-part stays accurate (demag factor right) while a spurious ANTISYMMETRIC part **grows
       with N** (independent ACA pivots on block (I,J) vs (J,I)). That asymmetry makes **CG/MINRES
       DIVERGE** (residual 8e1–1e3 at 3000 iters) while **GMRES converges** (92–151 iters) — which is
       exactly why the Newton step (GMRES) always worked but the warmstart (MINRES) failed.
       (`C:/temp/char_plusN.py`, `char_sym.py`.)
     - **Fix B (partial, committed): warmstart MINRES → GMRES.** The warmstart now converges
       (`info=0`), Newton starts at relF≈0.5 (not 1.4), iters drop **27 → 19**, solve **113 → 82 s** on
       cyoke h=0.006. NOT a full mesh-independence fix — the globalization still costs ~14 iters.
     - **Fix B RESOLVED (the asymmetry was just a too-LOOSE ACA tolerance).** Measured: the Hg asymmetry
       scales directly with `eps` (`1e-6`→15%, `1e-8`→1.6%, `1e-10`→3e-5, `1e-12`→1.5e-6). The default
       `gram_eps=1e-7` gave the ~5% asymmetric/inaccurate operator → the Newton tangent was inconsistent
       → slow LINEAR convergence (19–37 iters) + the 0.6% wrong Mz. **Fix: default `gram_eps` 1e-7 →
       1e-10** in `solve_nonlinear_newton_scalable`. End-to-end on the C-yoke this gives a **mesh-INDEPENDENT
       6-iter quadratic convergence** and the drift vanishes:

       | maxh | ndof | Mz | iters (1e-7→1e-10) | solve (1e-7→1e-10) |
       |---|---|---|---|---|
       | 0.008 | 5949 | 584635 | 6 → 6 | — → 4.1 s |
       | 0.006 | 10759 | 585538 | 27 → **6** | 87 → **10.4 s** (8×) |
       | 0.005 | 15726 | 586032 | 37 → **6** | — → 23.9 s |

       Mz is now monotonic + converging (584.6k→585.5k→586.0k, the dense trend), no false-convergence
       drift. (The GMRES warmstart from the prior commit stays as robustness; the near/far `near_factor`
       offsets the modestly-higher tight-eps build cost.) The ACA tolerance was the cause of the
       OUTER-Newton (tangent-consistency) degradation — NOT B1 (−N material, non-physical).

     - **"SCALE WALL at 44k" — RESOLVED (2026-06-09): it was a GMRES `restart=50` artifact, NOT a wall,
       NOT the loops.** A cautious intermediate point (cyoke h=0.004, ndof 38383) first showed the
       warmstart not converging (`gmres_info=20`). Investigation (`C:/temp/validate_calderon.py`,
       `resolve_44k.py`, clean per-inner confirm): (i) at `gram_eps=1e-10` the inner +N GMRES with plain
       `M_mass⁻¹` is mesh/mu_r-INDEPENDENT — the loops are deflated by `M_mass⁻¹` (a loop-aware/Calderón
       preconditioner gives **ZERO** benefit, measured: identical iters → **B2 ruled out**); (ii) the
       inner-iter count grows only MILDLY (star-space demag spectrum): 31 @ ndof 8573 → **115 @ 38383**;
       (iii) the stall was purely that `restart=50` < 115 needed → restarted GMRES STAGNATES (at 38383:
       restart=50 → 1000 iters res 1.4e-3 fail; restart=200 → 115 iters res 7e-9 converged). **Fix:
       `gmres_restart` default 50 → 400** (warmstart + Newton step). With it, ndof 38383 completes in
       **7 Newton iters** (Mz 586732, on the trend), build 19.4 s + solve 119.6 s. So Fix A + B + the
       restart fix scale to **ndof ≈ 38383 confirmed**; for very large N the inner iters keep growing
       (~star-space), so the lever there is a **STAR-space preconditioner** (NOT the loops) or a larger
       restart — but no auxiliary-space machinery is needed at ≤ ~100k.

     - **CORRECTED COMPARISON (2026-06-09): the yano-type reference is NO-loop-star + BLOCK JACOBI, NOT
       H-ILU.** H-ILU is the SEPARATE loop-star `A_SS` "star block" solver (`RadHACApKMSCManager` mode 2)
       — conflating it with yano-type is the exact mistake this doc earlier warned against. yano-type
       scaled to 165600 at 2686 linear iters with **Block Jacobi** + the robust **修正反復法 / Picard**
       outer loop, which TOLERATES a loose inner solve (it does not need the inner solve to converge
       tightly). So my earlier "HDiv Newton + M_mass⁻¹ vs yano + H-ILU" comparison was apples-to-oranges
       on BOTH axes (outer method Newton-vs-Picard AND preconditioner M⁻¹-vs-Block-Jacobi — and H-ILU was
       the wrong solver anyway). A growing inner-iter count is EXPECTED for a simple preconditioner —
       yano-type has it too (2686 @ 165600); it just absorbs it via Picard. HDiv-VIM's Newton is LESS
       forgiving (the warmstart needs a converged inner solve), which together with my too-low maxiter
       cap (1000) is why it stalled at ndof~44k — NOT a fundamental wall. **The FAIR head-to-head is the
       SAME solver — 修正反復法/Picard + Block Jacobi for both — where HDiv-VIM's de-Rham loops (EXACTLY
       field-null on ANY mesh, ker B) should beat yano-type's distortion-element loops (field-null only
       on affine hexes; they carry field on distorted ones).** That comparison has NOT been run yet; it
       is the right next experiment.

So the corrected status: **BUILD is scalable + a clear win; SOLVE is mesh-INDEPENDENT (6-iter) + accurate
up to ndof ≈ 15726 (Fix A + Fix B), but the INNER +N solve needs a real preconditioner to go past ~20–40k
ndof toward the 165600 scale.** The head-to-head JSON is honest at the measured scale; do NOT read it as a
165600-DOF result.
## Milestones

- **M0 — parity gate + speed-gap measurement** *(START HERE; mostly measurement, low risk).* The
  definition-of-done above + the honest speed number. Until M0, "retire yano-type" is conference-ready
  (validated method) but the production gap is unquantified.
- **M1 — production module + public Radia API.** Move the validated solve out of `examples/vim/`
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
