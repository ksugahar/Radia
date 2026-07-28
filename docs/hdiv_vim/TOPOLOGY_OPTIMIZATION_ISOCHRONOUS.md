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

## 5. Existing assets (reuse, do not rebuild)

| Asset | Location | Role here |
|---|---|---|
| VIM linearization + LP + shape tangents | `radia.topology_optimization` (`linearize_vim_system`, `optimize_vim_lp`, `solve_lp_update`, per-family `production_*_derivatives`, `write_cubit_density_journal`; exercised by `tests/test_topology_optimization.py`) | update driver; the SHAPE route of Sec. 6; Cubit export of converged densities |
| Operator API | `radia.vim.DemagOperator` / `HDivSolver` / `build_charge_gram` | build-once `N`, weighted mass, solves |
| `gL` measurement chain | `docs/clebsch_hodograph/edge_focusing_tracking.*` + its golden | orbit points, weights, verification objective |
| Coil source | `radia.coil_builder` | `H_ext` without a coil mesh |
| Knowledge | `radia_mcp.topology_optimization`, `radia_mcp.accelerator` (`edge_focusing_tracking`) | formulation and lane conventions |
| Validation homes | `validation_test/feec/`, `validation_test/vim_coupled/` | where the Stage gates land |

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
| 0 | Reconnaissance of `radia.topology_optimization` (LP semantics, reuse points, filter status) | design memo; no duplicated machinery |
| 1 | `J`, `F^T c`, adjoint on a sanity geometry (sphere/bar in a uniform + gradient drive) | adjoint gradient == finite differences to 1e-6 class on every tested `s_e` |
| 2 | Isochronous sector case: orbit, `BL` constraints, volume budget; design loop with warm-started CG and block forward+adjoint rhs | monotone constrained descent; per-iterate wall time ms-to-seconds class at study scale |
| 3 | **Final verification protocol**: remesh the converged design iron-only (void REMOVED -- exact-void gold standard, discharging the `O(h)` embedded bias), re-measure with the full RK4 Hill chain | removed-void `gL` inside a stated band of the embedded prediction; band quantifies the ersatz error honestly |
| 4 | Promotion: `validation_test/` lane with golden bands; companion result-bearing notebook in this directory; knowledge sync | ladder rules of CLAUDE.md |

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
