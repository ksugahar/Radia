# HDiv-MMM Topology Optimization for Isochronous Magnet Design

Status: **design document** (2026-07-28). The enabling structure is verified
numerically (Sec. 3, dated measurements); the objective/adjoint and the design
loop are planned (Sec. 7). Per the documentation policy this stays Markdown
until the method runs end-to-end, at which point a result-bearing companion
notebook is added next to the showcase notebooks in this directory.

---

## 1. Mission: why the objective is `gL`, not `BL`

The design driver is an **isochronous** machine. Isochronism fixes the average
bending field along radius to `<B>(r) = B0 gamma(r)`, which written as a field
index is

    n = (r/B) dB/dr = gamma^2 beta^2 = gamma^2 - 1  >  0,

so **isochronism itself forces vertical defocusing**. The vertical tune obeys

    nu_z^2 = -n + F (1 + 2 tan^2 zeta),

where `F` is the flutter and `zeta` the spiral angle of the hill/valley
boundary. Every unit of vertical focusing must be manufactured from `F` and
`zeta`; the **integrated gradient `gL` (focusing) is the scarce resource** the
optimizer buys, while `BL` (bending) is a *constraint* handed down by
isochronism, not a free objective. Because `tan(zeta)` enters `nu_z^2`
squared, a few percent of effective-edge error is 10-15 % of the focusing
budget -- the objective evaluation chain must resolve real fringe physics
(Sec. 2), not hard-edge surrogates.

## 2. Why HDiv-MMM is the right forward engine

HDiv-MMM (the RT1/BDM1 charge-Gram volume integral method, `radia.vim`) meshes
**only the iron**; the open boundary is exact through the Laplace kernel.

* The focusing objective is a **fringe integral extending into the field
  tail** -- exactly the quantity an air-box FEM truncates. HDiv-MMM has no air
  mesh at all.
* Measured head-to-head (2026-07-13, `docs/clebsch_hodograph/
  edge_focusing_fem_results.json`, key `hdiv_vim_cross_check`): reproduces the
  reduced-Omega FEM edge-focusing `dK_in` to **0.8 %** at **~10x lower cost
  per case**, from an iron-only tet mesh and one batch `rad.Fld` map.
* **Amortization fits topology optimization perfectly**: the geometry operator
  `N = B^T G B` is built once (1.35 s on the 270-tet check problem; minutes at
  design scale), and each design iterate costs only a weighted-mass
  reassembly (measured 1-10 ms) plus one SPD solve.
* The field anywhere is an exact analytic functional of the solved
  magnetization -- the tracking-chain objective consumes it directly.

## 3. Verified foundations (dated; do not re-derive)

System form (from `src/radia/vim/_vim.py` and `src/core/rad_hacapk_hdiv.cpp`):

    A(s) x = M_s x + N x,     N = B^T G B  (geometry only),
    M_s    = HDiv mass weighted by the per-element  s_e = 1/chi_e.

| # | Property | Measurement (2026-07-26/28) |
|---|----------|------------------------------|
| 1 | Separation / reuse | `N` built once, reused for every design; uniform-chi matches the exact sphere demag to < 0.05 % |
| 2 | Per-element material hook | one `L2(order=0)` weight on the mass form; `N` untouched |
| 3 | Affine in the design variable | `\|\|M_(a+b)x - (M_a+M_b)x\|\| = 1.45e-16` -> `dA/ds_e` = the element mass matrix (exact local sensitivity) |
| 4 | Self-adjoint | `<r2, A r1>` vs `<r1, A r2>` bit-identical -> the adjoint reuses the same operator and preconditioner |
| 5 | Void is clean | per-element void response is EXACTLY the physical `chi_v * H` over 5 decades (zero anomalous leak); the historical "1.54 % residual" was a test artifact of a spatial-IfPos split |
| 6 | `chi_min` is free | with the `s`-weighted mass-Riesz preconditioner CG iterations DECREASE as `chi_v -> 0` (74 at 1e-1 -> 46 at 1e-6); use `chi_min = 1e-6` |
| 7 | Interface bias | embedded designs carry an `O(h)` iron under-magnetization at iron/void boundaries (charge smearing into the last iron layer: -7.6/-5.5/-3.6 % at maxh .35/.22/.14 on the smooth-aligned sphere) -> Stage-3 verification rule |
| 8 | Multi-region meshes | legal since the charge-layer fix (internal interface faces carry no single-sided charge; commit 327ce4aa8, locked by `tests/test_hdiv_vim_multiregion_interface.py`) |

Solver behavior: SPD CG + mass-Riesz, kappa ~ 5, iteration counts flat in `N`;
HACApK build ~ `N^1.23` and build-dominated -- which is exactly why the
build-once/iterate-many topology loop is the method's best-case workload.

## 4. Formulation

**Design variable.** `s_e = 1/chi_e` per element (an `L2(order=0)` weight),
bounded by `s_iron = 1/chi_iron` and `s_void = 1/chi_min`, `chi_min = 1e-6`.
Design boundaries are element-aligned by construction (no straddling-element
artifacts). Standard density filtering (Helmholtz/convolution on `s_e`) plus
projection supplies minimum-feature control; the filter radius is a
manufacturability input, not a numerical crutch.

**State.**  `(M_s + N) m = P H_ext`, SPD; `H_ext` from the coil model
(Biot-Savart `CoilBuilder` source, no coil mesh).

**Objective.** The field is LINEAR in the magnetization: `B(x) = F(x) m +
B_coil(x)` with `F` the analytic magnetization-to-field map. The optimization
objective is the linearized focusing functional along the reference orbit,

    J(m) = c^T F m + const,      c = quadrature weights of  dB_z/dx  on orbit
                                     points of the verified tracking chain,

so `dJ/dm = F^T c` is a constant vector. The **verification-level** objective
remains the full RK4 linearized Hill integral of
`docs/clebsch_hodograph/edge_focusing_tracking.{ipynb,py}` (golden:
`tests/feec/test_edge_focusing_tracking.py`); the linear functional is its
optimization-grade surrogate evaluated on the same orbit discretization.

**Isochronism constraint.** `<B>(r_i) = B0 gamma(r_i)` on a family of orbit
radii -- also linear functionals of `m`. Enter as equality constraints (or a
penalized quadratic) alongside an iron volume/mass budget.

**Adjoint and sensitivity.** Self-adjointness makes the adjoint mechanical:

    A lambda = F^T c        (same operator, same preconditioner as the state),
    dJ/ds_e  = - lambda^T M^(e) m      (exact, local; M^(e) = element mass).

Forward and adjoint share `A` -> solve as one block of right-hand sides.

**Update.** Candidate drivers already in the repository:
`radia.topology_optimization.optimize_vim_lp` / `solve_lp_update` (LP with
move limits) or a standard OC/MMA step on the filtered variable. The choice is
a Stage-0 reconnaissance decision, not new research.

## 5. Existing assets — Stage-0 reconnaissance findings (2026-07-28, DONE)

Full read of `radia.topology_optimization` (651 lines, single module) plus its
consumers (`tests/test_topology_optimization.py`,
`radia.magnetic_shield_optimization`).

**EXISTS — reuse as-is:**

| Asset | What it is | Role here |
|---|---|---|
| Shape route, complete | GetTrafo mode sampling -> analytic `dG/dM/dB` (C++ kernels, tet/hex/wedge, dense or streaming directional H-matrices) -> `dA` -> `production_vim_rms_adjoint_gradient_streaming` (state CG + adjoint CG on the same operator + per-mode contraction) | Sec.-6 SHAPE polish stage; also the PROOF PATTERN for adjoint solves in this codebase |
| `solve_lp_update` | bounded LP (HiGHS), move limits, volume budget, **extra `A_ub`/`b_ub` rows** | update step; the isochronism constraints go into the existing `A_ub` slot (equalities as +/- row pairs) |
| `optimize_vim_lp` | sequential linearize -> LP loop with history/convergence | reusable AS-IS by an adjoint-backed callback: with a scalar objective, `weights=[1.0]` and `response_jacobian` = the 1 x n_el gradient row satisfy its interface |
| `write_cubit_density_journal` | thresholded density -> Cubit block journal | converged-design export toward the Stage-3 iron-only remesh |
| `linearize_vim_system` | dense FORWARD sensitivities (one solve per design cell) | small/sheet problems only; NOT the m-scale route (kept, not extended) |
| `radia.magnetic_shield_optimization` | an application module composing the above | the COMPOSITION PRECEDENT this project mirrors |
| `radia.vim.build_charge_gram` / `DemagOperator` | build-once `N`, weighted mass, mass-Riesz solves | forward engine |
| `docs/clebsch_hodograph/edge_focusing_tracking.*` + golden | verified orbit/tracking chain | orbit points/weights; Stage-3 verification objective |
| `radia.coil_builder` | Biot-Savart coil source | `H_ext` without a coil mesh |

**MISSING — the Stage-1 build list (verified non-duplicating):**

1. **Per-element DENSITY adjoint at scale.**  The existing density path is
   forward-mode (O(n_el) solves, dense `A`); the existing adjoint is
   shape-mode only.  To build: one state + one adjoint solve of
   `A(s) = M_s + N` (weighted-mass-Riesz preconditioner, warm starts,
   block rhs), then ALL element sensitivities in one call:
   `dJ/ds_e = -lambda^T M^(e) m = -Integrate(gf_lambda * gf_m,
   element_wise=True)` -- the element mass contraction IS the element-wise
   integral of the product of the two solution fields.  No per-mode
   machinery, no P-matrix.
2. **`gL` objective vector `F^T c` by kernel reciprocity.**  The dipole
   kernel is symmetric, so the adjoint load is the mass projection of the
   field of a weighted DIPOLE-PAIR array on the orbit points
   (`dB_z/dx` realized as +/- pairs): an analytic CoefficientFunction /
   `rad.Fld` evaluation.  The dense magnetization-to-field matrix `F` is
   never formed.
3. **Density -> `s` mapping, filter, projection.**  `s(rho)` interpolation
   with the validated floor `chi_min = 1e-6`; no density filter exists
   anywhere yet (Helmholtz/convolution + projection to build).
4. Stage-3 iron-only remesh/verification driver.

**Placement decision (mutex):** `topology_optimization.py` is under active
co-agent development (recent commits: streaming adjoint contractions, Cubit
remeshing).  Stage 1 therefore lands as a separate APPLICATION module
`radia/isochronous_topopt.py` -- mirroring the `magnetic_shield_optimization`
precedent, importing the LP driver, editing nothing in the shared module.

**Stage-1 status (2026-07-28): build items 1 and 2 SHIPPED** in
`radia.isochronous_topopt` (plus the `density_to_s` floor of item 3 and the
warm-start lever of Sec. 8).

**Stage-2 status (2026-07-28): item 3 SHIPPED** (`HelmholtzFilter` with the
self-adjoint chain and the clipped-with-exact-chain-rule composition inside
`optimize_density`; smooth tanh PROJECTION remains future polish -- the
converged surrogate designs carry ~45 % intermediate-density elements at
verification-mesh coarseness).

**Stage-3 status (2026-07-28): item 4 SHIPPED** (`iron_only_mesh` +
`verify_design_iron_only` + SIMP `penalty`).  The build list of Sec. 5 is
complete; what remains is the STUDY itself (Sec. 6.2: real isochronous
profile, penalty/projection continuation, study-scale mdx runs, RK4 Hill
re-measure on the verified iron-only design).

Knowledge homes unchanged: `radia_mcp.topology_optimization`,
`radia_mcp.accelerator`; validation lanes land under `validation_test/`.

## 6. Two design-variable routes (both exist; ordered)

1. **Density route (this document's primary path)**: per-element `s_e` on a
   fixed mesh -- the verified foundations of Sec. 3 apply verbatim; topology
   is free (holes, hills/valleys, spiral edges emerge).
2. **Shape route (final polish)**: the existing GetTrafo deformation-mode
   linearizations (gram/operator/rhs Jacobians per tet/hex/wedge family) give
   exact shape derivatives on a body-fitted mesh. Use AFTER the density
   topology freezes, to sharpen boundaries beyond the `O(h)` embedded-design
   accuracy -- ending exactly where Sec. 7 Stage 3 verifies.

## 7. Execution plan

| Stage | Content | Gate |
|---|---|---|
| 0 | **DONE 2026-07-28.** Reconnaissance of `radia.topology_optimization` (LP semantics, reuse points, filter status) | Sec. 5 findings: adjoint exists only shape-side; density adjoint + `F^T c` + filter are the build list; LP driver and `A_ub` constraint slot reused as-is; new work lands in `radia/isochronous_topopt.py` |
| 1 | **DONE 2026-07-28.** `radia.isochronous_topopt` shipped: `DensityAdjointVIM` (one weighted-mass assembly + factorization, state+adjoint CG, ALL element sensitivities from one element-wise `Integrate`), `field_functional_load` (`F^T c` by dipole-pair reciprocity; dense `F` never formed), `gradient_pair_points`, `density_to_s` (validated `chi_min = 1e-6` floor) + chain rule, warm-started CG with the rhs-anchored ABSOLUTE tolerance (krylovspace's relative tol is anchored to the run's first residual, which makes warm restarts harder, not cheaper -- anchor to `\|\|rhs\|\|_pre`; upstream `abstol=` kwarg is broken, use `atol=`) | gate PASSED (unit ball, log-uniform `s` in [1e-2, 1]): adjoint == central FD to 8.1e-10 directional / 3.9e-7 worst per-element (the smallest-gradient element, FD noise floor); reciprocity load == the independent C++ charge evaluator to 1.1e-10 at `bonus_intorder=10`; locked by `tests/test_isochronous_topopt.py` (9 tests) |
| 2 | **DONE 2026-07-28 (verification scale).** `optimize_density` shipped in `radia.isochronous_topopt`: trust-region SLP (LP via `solve_lp_update` with normalized `A_ub` constraint rows -- Tesla-scale rows sit below HiGHS's ABSOLUTE feasibility tolerance and must be scaled to O(1)), absolute engineering bands (default 0.5 % of target) with restore/hold modes, acceptance = monotone J AND (violation <= 1.25 band OR strict geometric decrease) -- no ratchet path; `HelmholtzFilter` (self-adjoint chain, FD-locked; filtered density clipped to [0,1] with the piecewise-exact chain rule -- the P1 realization undershoots ~1e-2 at bang-bang transitions), `DensityAdjointVIM.linearize` (state + K adjoints on ONE factorization, warm-started), per-point-direction `gradient_pair_points` + `orbit_arc_points` | gate MET at verification scale: sector-pole surrogate (194 tets; arc-orbit radial `dB_z/dr` objective, two mean-`B_z` arc constraints, volume budget) J strictly MONOTONE **+16.1 %** riding the active band at 1.06 x band (ball: +0.7 %, peak 1.24 x band <= the 1.25 cap); median 115-120 ms/iterate incl. trial evaluations (informal LAB; ~50 warm CG iters/solve, LP 3-5 ms); locked by `tests/test_isochronous_topopt.py` (15 tests). STUDY-SCALE timing = the planned mdx job (Sec. 8); rejected schemes recorded in the docstring (fixed-move SLP limit-cycles; violation-relative shrink bands lock into a ratcheting hold mode) |
| 3 | **DONE 2026-07-28 (protocol shipped; full RK4 Hill re-measure joins the Sec.-6.2 design study).** `iron_only_mesh` (kept-element extraction into a NEW straight-tet netgen mesh: exact void removal at the design's own discretization; volume identity 4e-16; NGSolve-vs-netgen element-order verified per call, fail-loud; SURFACE-ORIENTATION TRAP documented -- netgen stores boundary triangles right-hand-outward for `domin=1/domout=0`, the Radia `TETRA_FACES` handedness is opposite and flips every charge into runaway magnetization, `<Mz>` 18-34 instead of 2.2) + `verify_design_iron_only` (threshold -> matched-0/1 embedded vs exact-void bands per functional, loads rebuilt from BUILDERS on both spaces) + SIMP `penalty` in `density_to_s`/`optimize_density` | gate MET: the protocol MEASURES and STATES the band per design. Hemisphere machinery validation: matched-shape ersatz `<Mz>` band -8.3/-8.8 % (maxh .35/.22, pre-asymptotic jagged interface; the smooth-interface case is the established O(h)); jagged-vs-smooth +4.6/+7.7 % (STAIRCASE boundaries do not converge to the smooth body -- manufacturing-shape numbers need the Sec.-6 shape route). Sector surrogate: matched-0/1 FUNCTIONAL band **+0.69 %** (filtered design) / +10.2 % (bang-bang, more interface); the protocol EXPOSED the gray-design gap (continuous-vs-threshold -96 %, constraints destroyed) -> `penalty=3` shipped (intermediate 0.73 -> 0.28, gap -91 -> -33 %, iron-only constraint shift -81 -> -23 % at 30 move-limited iterations); threshold-READY designs need penalty/projection continuation at study scale = the Sec.-6.2 design study. Locked by `tests/test_isochronous_topopt.py` (19 tests) |
| 4 | **DONE 2026-07-29.** Promotion complete per the ladder: `validation_test/isochronous_topopt/` (adjoint-gate lane + design-loop/verification lane with golden bands from the measured values, committed record JSON `results_design_loop_lane.json`; 6 lane tests, ~1 min); result-bearing companion notebook `isochronous_topopt.ipynb` in this directory (executed: FD 1.4e-8 on its direction / reciprocity 1.1e-10 / +16.10 % monotone / ersatz +0.69 %, with history+density figures; sidecar `isochronous_topopt_result.json`, `radia.notebook_result.v1`); knowledge synced to `radia_mcp.topology_optimization` (`topology_opt_applications(topic="isochronous")`: API map + the six measured traps) | ladder rules of CLAUDE.md; the remaining WORK is the STUDY itself (real isochronism profile, penalty/projection continuation, study-scale mdx runs, RK4 Hill re-measure) -- tracked as the Sec.-6.2 design study, not a stage |

In parallel: the design-scale timing of the Python assembled route runs on
idle mdx (the same job doubles as the CUDA-lane baseline), and the CUDA lane
proceeds through its Phase 0-1 (Sec. 8).

## 8. Performance plan

Per-iterate cost = mass reweight (ms) + one SPD solve with two right-hand
sides (state + adjoint). Levers, in order of cost:

1. **Warm starts across design iterates** (small design steps -> CG from the
   previous solution; expected 3-10x fewer iterations).
2. **Block rhs** (state + adjoint together; ~2x).
3. **C++ vector `inv_chi`** so the fast native path takes per-element
   material (today scalar, `rad_hacapk_hdiv.cpp`; codex-owned -- coordinate).
4. **Krylov recycling / deflation** across iterates (optional).
5. **CUDA lane** (active, staged): the hot kernel is the HACApK H-matrix
   matvec (batched low-rank GEMM -- higher arithmetic intensity than FEM
   SpMV, i.e. GPU-friendlier than the FEM alternative). Measured GPU facts
   (2026-07-28): LAB Quadro RTX 5000 FP64 0.15 / FP32 2.49 TFLOPS (dev box
   only); mdx/hibino have no GPU; an A100 exists dormant (FP64-native ~65x
   the Quadro). Phase 0 = mdx baseline + `RadHACApKBase::MatVec` backend
   interface design; Phase 1 = CuPy block-matvec prototype on real ACA block
   distributions -> the "CG seconds on A100" prediction table that decides
   waking the A100; Phase 2 = opt-in C++ CUDA build (PyPI wheel stays
   CPU-only); Phase 3 = production on the woken A100. CPU levers 1-4 compose
   multiplicatively with the GPU and stay in the plan.

## 9. Known limits and recorded negatives (do not re-walk)

* **Field-EFB slope** for edge focusing: characterized negative (wrong sign,
  many times `tan(beta)`); the RK4 Hill integral chain is the measurement.
* **Partial `(x, A_z)` hodograph as a solving formulation**: evaluated and
  declined (no Green's function/BIE; compactification geometrically
  restricted); retained only as a design-reading chart.
* **Pyramid elements**: NGSolve HDiv NOT_IMPLEMENTED (checked through the
  6.2.2606 nightlies) -> tet/hex/wedge meshes only; no mixed hex/tet
  transition layers.
* **Curved hexes**: the accepted 0.78 % demag-spectrum leak on strongly
  curved hexes; tet `Curve(2)` is the supported curved path.
* **Embedded-design accuracy** is `O(h)` at design boundaries (Sec. 3.7);
  never report final numbers from the embedded model -- Stage 3 exists for
  that.
* **Staircase boundaries are their own limit**: the iron-only extraction
  keeps the design's jagged element boundary, and jagged does NOT converge
  to the smooth body (hemisphere `<Mz>`: +4.6/+7.7 % vs the smooth-meshed
  hemisphere at maxh .35/.22 -- the staircase surface has its own
  homogenized boundary response, Schwarz-lantern style).  The exact-void
  verification is the gold standard AT THE DESIGN'S OWN DISCRETIZATION;
  manufacturing-shape numbers come from the Sec.-6 shape route on a
  smoothed body-fitted remesh.
* **Netgen surface-element handedness**: boundary triangles are stored
  right-hand-OUTWARD for `FaceDescriptor(domin=1, domout=0)`; the Radia
  `TETRA_FACES` ordering is the opposite handedness in this context and
  produces runaway magnetization if used for surface reconstruction
  (`iron_only_mesh` carries the correct `_TET_BOUNDARY_FACES`).
* **Ownership**: `src/core/rad_hacapk_hdiv.*`, `src/radia/vim/**` are the
  co-agent's active area -- coordinate before editing; heavy timings run on
  idle mdx/hibino only.

## 10. References

* Le-Van T. et al., IEEE Trans. Magn. 51(7) (2015) -- facet-element magnetic
  moment VIM (formulation positioning); IEEE Trans. Magn. 50(2) (2014) --
  ACA-compressed VIM.
* Bendsoe M. P., Sigmund O., "Topology Optimization: Theory, Methods and
  Applications", Springer (2003) -- density method, filtering, OC/MMA.
* The workshop paper of this method family (HDiv-MMM construction,
  amortization, loop-free hysteresis) and its committed goldens
  (`validation_test/feec/`, `validation_test/hysteresis/`).
