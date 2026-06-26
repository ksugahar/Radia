# MSC Operator Null Space, Conditioning, and Loop-Mode Deflation

This note documents the eigenvalue structure of the MMM/MSC interaction
operator, why the discrete solution is contaminated by spurious "loop"
modes at high permeability, and how those modes can be removed. It is the
solver-side foundation of the CEFC 2026 eigenvalue study.

Reproducible scripts: REMOVED 2026-06-27. The `examples/mmm_eigenvalue_study/`
corpus (and its docs notebook) was deleted: the loop-deflation runtime API is
gone and **MMMM is the official MMM H-matrix route**, so the loop-deflation
exploration has no forward value to ship. This note is retained as the theory
record only; the inline `Script: *.py` references below are historical.

> **Status (updated 2026-06-23): the RUNTIME loop-mode handling once described here
> -- matrix-free deflation (`SetHACApKDeflation` / `SetDeflateNullspace`), the
> alpha-free loop-star gauge (`SolveLoopStar`, the `A_SS = S^T A S` star block),
> and the Helmholtz-Hodge loop projection (`SetLoopProjection`) -- was REMOVED
> from the solver.** The live surface-charge MSC backend has since consolidated on
> **multipole-moment MMM** (hex/wedge/pyramid), while **HDiv-VIM** remains the FEEC
> complementary backend: its charge map `B` makes the loop space exactly `ker(B)`,
> field-null BY CONSTRUCTION via the de Rham complex, `mu_r`-independent, with no
> runtime deflation / gauge / projection. The null-space THEORY and conditioning
> analysis below remain valid, but the API-level "how to remove the modes at
> runtime" material is historical; class names such as `RadHACApKMSCManager` are
> pre-2026-06-23 names (`RadHACApKMMMManager` is the live MMM manager).

---

## 1. The MSC operator and its system matrix

`rad.BuildMatrix(obj)` + `rad.GetInteractMatrix(handle)` returns the
**geometric** interaction matrix `N` (stored `+N`, the demagnetization
tensor). It is **independent of permeability** (building it at `mu_r=10`
vs `mu_r=5000` gives the identical matrix) and **non-symmetric** (MSC
collocates the normal field at face centers).

The system matrix actually solved is

```
A = diag(1/(mu_r - 1)) - N        (= -N + diag(1/chi), chi = mu_r - 1)
```

In the digest/paper notation `A_ij = delta_ij/(mu_r-1) + G_ij`, the "G"
there equals `-N`; the `+N/-N` flip is a storage convention.
`dof = sum of per-element DOF` (tet 3, wedge 5, hex 6).

---

## 2. The null space = cycle space of the element-adjacency graph

`N` has an **exact null space**: surface-charge distributions that produce
zero normal field at every collocation point. Its dimension equals the
**first Betti number (cycle rank) of the element-adjacency graph** (nodes
= elements, edges = shared internal faces):

```
null_dim = F_internal - N_elem + 1
```

verified exactly (machine precision) for a cube-element block:

| block   | 1x1x1 | 2x1x1 | 2x2x1 | 2x2x2 | 3x3x3 | 4x4x2 | 4x4x4 | 5x5x5 |
|---------|------:|------:|------:|------:|------:|------:|------:|------:|
| N_elem  | 1     | 2     | 4     | 8     | 27    | 32    | 64    | 125   |
| null    | 0     | 0     | 1     | 5     | 28    | 33    | 81    | 176   |

A single element and a 1-D chain have NO null space; it first appears at a
`2x2` arrangement (one mode = the minimal circulating "loop" of four
elements).

**Physical derivation.** A null mode puts opposite charges on the two
sides of each shared internal face (the field of the coincident faces
cancels) AND zero net charge per element (so the compensating centroid
point charge `m_e = -sum_f sigma_{e,f}` vanishes). Defining one "flux"
variable per shared face, "zero net charge per element" is a discrete
divergence-free condition, whose solution space is exactly the cycle
(circulation) space of the graph. Minimal cycle = `2x2` plaquette.

Scripts: `block_spectrum.py`, `nullspace_scan.py`, `nullmode_cycle_basis.py`.

---

## 3. Conditioning: high permeability is the worst case

Because `A V = (1/chi) V` for any null vector `V` of `N`, the null modes
become eigenvalues of `A` located **exactly at `1/(mu_r - 1)`**. They are
the smallest eigenvalues of `A`, so

```
cond(A) ~ mu_r          (largest eigenvalue is O(1))
```

As `mu_r -> infinity` the cluster collapses to zero and `A` becomes
singular. The corresponding eigenvectors are circulating ("loop")
magnetization patterns that radiate almost no external field.

---

## 4. "Beautiful -> ugly" iterative behavior

For a high-`mu_r` soft-iron block in a uniform field, BiCGSTAB from a zero
initial guess passes a physically aligned solution at few iterations
("beautiful") and then drifts to a loop-dominated converged / LU solution
("ugly"), because Krylov resolves the large-eigenvalue (physical)
components first and the near-null (loop) components last. The effect grows
with `mu_r`: mean alignment of `M` with the applied field falls from

- `0.94 (few iter) -> 0.92 (converged)` at `mu_r = 5000`
- `0.94 (few iter) -> 0.14 (converged)` at `mu_r = 1e5`

Reproducible with Radia APIs only (`Solve` + `ObjM` + `SolverConfig`,
using `bicgstab_tol` as the iteration knob). Script: `beautiful_ugly_viz.py`.

This is a **formulation / conditioning property, not a solver bug** — LU
shows the same loop-dominated converged solution.

---

## 5. Does the H-matrix remove the bad modes? (No.)

A controlled low-rank (SVD) truncation of the far-field blocks does **not**
regularize the near-null spectrum. At the HACApK working tolerance
`eps = 1e-4` the smallest eigenvalues / condition number are essentially
unchanged (`3.9e4 -> 4.0e4` on a 6x6x6 block); more aggressive truncation
moves them **closer** to zero. So the data-sparse representation buys
memory and time, not spectral regularization. Script:
`hmatrix_truncation_spectrum.py`.

To validate this against the **actual** HACApK ACA+ operator (not the SVD
proxy), the C++ binding `rad.HMatrixDensify(handle)` densifies the real
HACApK system matrix `A = -N + diag(1/chi)` by applying its `MatVec` to
unit vectors (original DOF ordering), so its full spectrum can be compared
with the exact dense `A`. Script: `hmatrix_validate.py`.

**Result** (8x8x8, `mu_r=1e5`): the densified real ACA+ operator matches the
exact dense `A` to the ACA tolerance (relative error `5e-7` at `eps=1e-4`,
diagonals match exactly -- which also confirms the `A = diag(1/chi) - N`
sign convention). Its spectrum, however, perturbs the exact-zero cluster
into a *spread* of near-zero eigenvalues lying BELOW `1/(mu_r-1)`, worsening
the condition number: `7.7e5 -> 1.08e6` at `eps=1e-4`, and `-> 8.9e6` at
`eps=1e-2`. So the real ACA+ confirms the SVD proxy -- it does NOT regularize
the near-null subspace; it mildly worsens it. (The `rad.HMatrixDensify`
binding: `radTApplication::HMatrixDensify` -> `RadHACApKMSCManager::MatVec`.)

---

## 6. Removing the loop modes (deflation)

Two textbook techniques, validated on the dense system
(`nullmode_removal.py`):

1. **Projection deflation** of the solution: `sigma_clean = sigma - V(V^T sigma)`
   removes the loop content. Because `N V = 0`, the collocation field
   `N*sigma` is preserved exactly.
2. **Eigenvalue shift** (the recommended solver fix):
   `A_s = A + alpha * V (W^T V)^{-1} W^T` moves the near-zero cluster from
   `1/(mu_r-1)` up to `1/(mu_r-1) + alpha`. With `alpha ~ O(1)` this drops
   `cond(A)` from `7.7e5` to `46` on a 4x4x2 / `mu_r=1e5` block and makes a
   fully converged BiCGSTAB solve loop-free.

Here `V` / `W` are the right / left null spaces of `N` (SVD).

---

## 7. Scalable deflation: a LOCAL loop basis

The full null space is ~17-21% of the DOF, so a dense null basis `V` is
`O(N^2)` and would defeat HACApK. However the null projector `P = V V^T`
is **local** (decays ~2x per cell of element distance; ~50% of its mass is
within face-neighbors), so a **local basis exists**:

- **SCDM** (selected columns of `P` by QR-with-pivoting) gives a local
  basis that spans the null space exactly and is 97% localized within 1.5
  cells (`nullmode_local_deflation.py`).
- **Explicit 2x2 plaquette cycles** give the local basis **without any
  SVD**, directly from mesh topology: each plaquette circulation is an
  exact null mode (`||N*loop||/||loop|| ~ 1e-15`), and the plaquette set
  spans the null space (`rank == null_dim`, span error `~1e-15`). Each
  basis vector touches only 4 elements (8 shared-face DOF) -> `O(1)`
  support (`nullmode_cycle_basis.py`).

This makes a **matrix-free deflated / shifted BiCGSTAB** feasible inside
the HACApK solver: build the local cycle basis from mesh topology and
deflate without ever forming a dense null basis.

**Shipped (C++).** `rad.SetHACApKDeflation(offsets, dofs, signs, alpha)`
passes the sparse plaquette basis `L` (CSR) to the solver;
`RadHACApKMSCManager::MatVec` then applies the matrix-free symmetric shift
`A x + alpha * L (L^T x)` -- RAW `L`, NO orthonormalization, NO inverse, NO
SVD (`L L^T` acts only on `span(L)` = the null subspace, `O(nnz)=O(N)`).
`alpha~1` reconditions; a final projection `sigma - L(L^T L)^{-1}(L^T sigma)`
clears any residual. Validated end-to-end: the converged HACApK solution's
mean magnetization alignment rises from `0.22` (undeflated, loop-dominated)
to `0.94` (deflated) on a 4x4x2 / `mu_r=1e5` block, and from `-0.02` to
`0.97` on 4x4x4. Default `alpha=0` leaves the solver unchanged (no
regression). Tests: `test_deflated_hacapk.py` (manual basis),
`test_deflate_auto.py` (auto).

**General mesh (C++).** `rad.SetDeflateNullspace(True, alpha)` makes the
HACApK solve auto-build the cycle basis IN C++ from the mesh topology:
`RadHACApKMSCManager::BuildLoopBasis` finds adjacency + shared-face DOF from
coincident face centroids (`poly->FaceCenter[f]`, DOF `= dof_offset(e)+f`),
enumerates the element-graph's short cycles (3- and 4-cycles), and (for
multiply-connected bodies) completes them with the global belt loops -- see
the **belted-tree** paragraph in Section 8 -- then installs the loop CSR. This
is geometry- and topology-robust (no grid assumption, any genus): validated to
clean the solution on a SHEARED (non-axis-aligned) block (alignment
`0.07 -> 0.88`) where the Python column-cosine shared-DOF heuristic fails, and
on closed RINGS (see below). No Python pre-processing -- self-contained for any
conforming mesh. The number of installed cycles is reported in
`rad.GetSolveStats()['deflation_cycles']` (diagnostic).

---

## 8. Theoretical backbone: tree-cotree gauge (div side) / loop-star

The null space is not an accident -- it is the discrete **tree-cotree gauge**
structure of the magnetic-charge (div-side) formulation, dual to the gauge of
the FEM vector-potential (A) formulation.

| | FEM A-formulation (curl side) | MSC charge (div side, here) |
|---|---|---|
| variable | A on edges (1-form) | sigma on faces (~2-form) |
| null space | ker(rot)=im(grad) = **gradients** | **circulating (solenoidal) magnetization** |
| graph | primal (vertex-edge) | dual (element-face) |
| dim | N_node - 1 (+ b1) | F_int - N_elem + 1 |
| gauge fix | A=0 on a spanning tree of edges | a spanning tree of FACES; cotree faces = loops |
| reference | Albanese-Rubinacci 1988 | the div-side dual |

By the Helmholtz/Hodge decomposition `M = grad(phi) + curl(T) + harmonic`,
the charge-producing part is `grad(phi)` while the **solenoidal part
`curl(T)` produces no charge and no field** -- exactly the null space. The
loop (cotree) modes ARE the discrete `curl(T)` subspace; our plaquette /
short-cycle basis is a local basis for it. Deflation = removing this div-side
gauge freedom.

**Difference from a pure gauge.** In the A-formulation `ker(rot-rot)` is
exactly singular (a true gauge). Here the *field* operator `N` has the exact
null space (loops are invisible to the field), but the *system* matrix
`A = diag(1/(mu_r-1)) - N` lifts those modes to eigenvalue `1/(mu_r-1)`: the
constitutive term weakly breaks the gauge. So it is a NEAR-singularity whose
conditioning degrades as `mu_r -> infinity` -- the magnetostatic analog of the
**low-frequency breakdown** of MoM/EFIE, whose standard cure is the
**loop-star (loop-tree) decomposition**. Our deflation is the loop-star /
tree-cotree gauge removal applied matrix-free inside HACApK.

**Topology completion (belted tree) -- IMPLEMENTED.** The short-cycle
(plaquette / edge-ring) basis spans only the *locally generated* (contractible)
loops. For a **multiply-connected** iron body (first Betti number `b1 > 0`,
e.g. a closed O-core / ring) there are `b1` additional GLOBAL loops (homology
generators) that local cycles miss. `BuildLoopBasis` now completes the basis
with exactly those `b1` belt loops (a **belted tree**), so it spans the FULL
null space for any topology:

1. Index the element-graph edges and build a BFS spanning forest
   (parent/depth, component count `C`); the graph cycle rank is
   `b1_graph = E - V_active + C` (`= dim ker N`).
2. Seed a GF(2) basis (sparse rows keyed by pivot edge) with the short cycles;
   `rankShort` is the contractible dimension. The reduction stops early once
   `rankShort == b1_graph` (simply-connected fast exit: no cotree pass, and the
   installed DOF vectors are identical to the short-cycle-only path -- zero
   regression).
3. `beltNeeded = b1_graph - rankShort`. If `> 0`, walk the cotree edges; each
   one's fundamental cycle (tree path + the cotree edge) that is GF(2)-
   independent of the basis is a genuine belt loop -- install its alternating
   +/-1 DOF vector (a divergence-free circulation = exact null mode) and add it
   to the GF(2) basis. Exactly `beltNeeded` are added.

Verified (`examples/mmm_eigenvalue_study/belt_loop_validation.py`): a thin
1-element-thick ring has a single `C_N` adjacency cycle (ZERO short cycles), so
the belt completion is the ONLY source of its one null mode -- it installs
`deflation_cycles == 1 == null_dim` and a short-cycle-only basis would install
0. A tube (short cycles + 1 belt) installs `9 == null_dim` and restores the
converged alignment `0.19 -> 0.89`. Simply-connected controls (solid block,
cube `4^3 .. 8^3`) are unchanged (belt count 0; iteration counts and alignment
match the short-cycle-only results).

References: Albanese-Rubinacci (tree-cotree gauging, 1988); Bossavit (Whitney
forms, primal/dual de Rham); loop-star / low-frequency EFIE (Wilton, Vecchi,
Andriulli). See radia-mcp `differential_forms_homology` (tree_cotree),
`fem_gauge_open_boundary`, `bem_low_freq`.

## 8.5. Loop removal as a post-solve Hodge projection (`SetLoopProjection`, DEFAULT ON)

Sections 6-8 deflate the loop modes *during* the iteration (to recondition the
solve). But the loops are also **non-physical in the answer itself**: they are the
solenoidal `curl(T)` part of the magnetization (Section 8) -- circulating surface
charges that produce **no field**. If you want a physically clean `sigma` / `M`
with zero loop content, the cleanest route is a **Helmholtz-Hodge projection of the
converged solution** off the cycle space, independent of how the system was solved:

```
c solves (L^T L) c = L^T sigma        (CG; L = topological cycle basis = ker N)
sigma  <-  sigma - L c                 (remove the span(L) = loop component)
```

This is **default ON** (`rad.SetLoopProjection(True)` is the default; pass `False`
to keep the raw loop-containing `sigma`). Key properties:

- **Field-exact / transparent.** `N L = 0`, so `N*sigma` -- the on-element field --
  is unchanged: removing the loop does NOT change `rad.Fld()` (`dB/B ~ 1e-15`).
  Only the non-physical circulating part of the magnetization distribution changes.
  Verified: **102/102** core field + golden tests are unchanged with the default on.
- **Cheap, no slowdown.** `L^T L` is the loop **Gram** matrix -- geometric, sparse,
  `mu_r`-INDEPENDENT, and WELL-CONDITIONED (`cond ~ 1`, the cycle basis is
  near-orthonormal), the OPPOSITE of the system matrix `A` (`cond ~ mu_r`). The CG
  converges in a handful of iters (cube ~6-8, C-type ~54, `mu_r`-independent). Cost
  is negligible vs the solve (`wall_proj ~ wall_plain`).
- **Scalable.** Uses ONLY the `O(N)` topological cycle CSR (`BuildLoopBasis`) plus
  matrix-free `L`/`L^T` -- no H-matrix, no `MatVec`, no dense null basis. NOTE this
  uses the **topological** cycles (`E-V+C`), which is what "loop" means here; it is
  NOT the full numerical `ker(N)` (on some geometries `ker N` exceeds the cycle
  space by a few accidental near-null modes that GROW with `N` and are a separate,
  non-scalable issue -- those are NOT loops and are deliberately left alone).
- **All solver paths, consistently.** Method 2 (HACApK) projects in-manager after
  convergence; methods 0 (LU) and 1 (dense BiCGSTAB) call the static
  `RadHACApKMSCManager::ProjectOutLoopsStandalone` (a temporary manager bound to the
  interaction runs the same pure-sparse projection on the element `sigma`). So LU,
  dense BiCGSTAB, and HACApK all return the SAME loop-free `sigma`.
- **Nonlinear-safe.** `L = ker(N)` is `chi`-independent (pure geometry), so the same
  `L` projects at any `chi`. The projection is applied **ONCE after the nonlinear
  iteration converges**, NOT per Picard/Newton step -- so the `chi(H)` iteration is
  driven by the true (loop-included) `sigma` and only the final answer is made
  loop-free. (Projecting every step instead fed a loop-free `sigma` into the `chi`
  update and perturbed the field ~1%, and crashed on the multiply-connected 1/4
  C-type -- both fixed by the post-convergence placement.) Verified on a C-type
  nonlinear `MatSatIsoTab` solve (57 Picard iters): no crash, loop removed 8 orders,
  field unchanged `dB/B = 8.7e-15`; `tests/test_hysteresis.py` (25 cases) green.
- **Interaction with the gauges.** Auto-skipped when `SetLoopStarGauge` (keeps loops
  by design) or the loop-deflated block-Jacobi gauge is active -- those already
  handle the loop content, so the post-projection would double-process.

Diagnostics: `rad.GetLoopProjStats()` -> `{n_loop, cg_iters, loop_before,
loop_after, loop_frac}`. `loop_frac` (`||L c|| / ||sigma||`) reports how much of the
solved `sigma` was circulating (non-physical): typically `~0.26` at `mu_r=2` rising
to `~0.99` at `mu_r >= 1e4` -- i.e. at high permeability almost ALL of the raw
`sigma` is loop content, which is exactly why the field is loop-robust but the raw
magnetization distribution is not. Test: `tests/test_loop_projection.py`.

## 9. Practical guidance

- The converged `sigma` is loop-free by default (`SetLoopProjection`, Section 8.5);
  `rad.Fld()` is unchanged either way (`N L = 0`). Pass `SetLoopProjection(False)`
  only if you specifically need the raw loop-containing magnetization.
- At very high `mu_r` the *iteration* can still be slow / ill-conditioned (the
  near-singular `1/(mu_r-1)` loop eigenvalues); the post-solve projection cleans the
  ANSWER but does not speed the SOLVE. For solve speed at high `mu_r` use a deflated
  / eigenvalue-shifted iteration (Sections 6-8) or keep `mu_r` physical.
- HACApK's ACA does not regularize the spectrum — in-solve deflation must be added
  explicitly (a local cycle basis keeps it `O(N)`).

---

## References

- Hackbusch, "A Sparse Matrix Arithmetic Based on H-Matrices", Computing 62 (1999).
- Ida-Iwashita-Mifune-Takahashi, J. Inf. Process. 22(4) (2014) — HACApK.
- Yano, J. Magn. Soc. Jpn. 47 (2023) — MSC formulation.
- Loop-star / cycle-space view: standard tree-cotree decomposition of
  computational electromagnetics.
