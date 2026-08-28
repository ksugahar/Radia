# HDiv-MMM implementation review

Review of the HDiv-MMM charge-Gram implementation, and of the ESRF example-5
C-type magnet cross-route discrepancy that it turned into.

- Reviewed: `src/core/rad_hdiv_vim.cpp` (2,495 lines, analytic kernels) and
  `src/core/rad_hacapk_hdiv.cpp` (11,834 lines, H-matrix and nonlinear solvers).
- Lenses: correctness, policy conformance, design/maintainability, performance.
- Branches: `claude/hdiv-mmm-review-fixes` (`ef986f0a6`),
  `claude/cubit-vol-main-path` (`2c939a8ed`, `4b30479b6`).

Every number below is measured. Wall-clock figures are deliberately absent:
lab policy puts timing on an idle mdx/hibino, never on LAB.

---

## 1. C++ implementation review

### 1.1 Fixed

| # | Finding | Site |
|---|---|---|
| B1 | Fill-time state (`cHACApK_set_sym_fill`, `m_fillNormalized`, `m_sigmaActive`) had no RAII guard | `rad_hacapk_hdiv.cpp` |
| B2 | CG had no breakdown detection on the `p^T A p` denominator | `SolveLinearMaterial` |
| B3 | `ComputeChargeSigma` silently kept sigma = 1 on a pathological diagonal | `ComputeChargeSigma` |
| B4 | `TriMoment1`, `TriMoment2`, `LinTriField`, `QuadTriFieldBasis` divided by an unguarded normal length | `rad_hdiv_vim.cpp` |
| B6 | The fixed 84-entry affine scratch had no bound check | `HexPolyMulLinear*` |
| Perf1 | Hot kernels zeroed worst-case (degree-18) stack blocks on every call | `rad_hdiv_vim.cpp` |
| Perf2 | `getenv` per matrix entry and per block lookup inside the parallel fill | `rad_hacapk_hdiv.cpp` |

**B1 — the entry oracle really does throw.** `BuildHMatrix` asserted that "no
exceptions cross the C fill", but `GetInteractionMatrixElementRaw` raises
`std::out_of_range` on an out-of-range index, `GetHexBlock` raises
`std::invalid_argument` past 63 images, and the affine map helpers raise
`std::logic_error` on singular geometry. Those unwind the pure-C frames of
`cHACApK_base.c`, whose ACA fill loops call `cHACApK_entry_ij`, so the
post-call reset was skipped. A stranded sym-fill flag makes the *next*
H-matrix build in the process — PEEC or BEM, which are not symmetric fills —
silently drop its lower leaves. This fix has since landed on main together
with a regression test that checks `_TestPEECHACApKSanity` after a failed
fill.

**B2 — breakdown was indistinguishable from slow convergence.** A
non-positive or non-finite `p^T A p` means the operator lost positive
definiteness or produced NaN. Previously the iterate silently became NaN and
the run surfaced as "did not converge in maxit", sending the reader after a
tuning problem that did not exist. The guard fires in practice; see 3.3.

**Perf1 — scratch sizing.** `POLY_MAX_MOMENTS` is 1,330. At the degree-2
charges that dominate production only a small prefix is used:

| Scratch | Was | Now (degree 2) |
|---|---|---|
| `face_moments[4][1330]` | 42,560 B | 320 B |
| per-moment `poly[19][19]` | 2,888 B | 48 B |
| `poly2_mul_linear` `tmp` | 2,888 B | 80 B |
| `TetMomentMemo` | 61,731 B | 27 B |

These sit inside the H-matrix fill.

**B6 — zero headroom, guarded only from Python.** `HEX_AFFINE_POLY_N` is 84,
exactly the total-degree-6 moment count, which is exactly what order-2
charges reach: `m_hexAffinePolyCount = (3*order+1)(3*order+2)(3*order+3)/6`
is 84 at order 2. `HexPolyIdx(7,0,0)` is already 119. The only bound was the
Python wrapper's order-above-2 guard, which the MATLAB/MEX entry points do
not go through.

### 1.2 Reported, not fixed

- **Determinism is inconsistent inside one function.** `SolveLinearMaterial`
  has a bit-deterministic dot product (fixed 4,096-element blocks,
  compensated summation, ordered reduction), but its transpose scatter falls
  back to `ngcore::AtomicAdd` when no charge map is configured, and the whole
  nonlinear path uses atomics. A deterministic transposed CSR already exists
  in the class and is used when a charge map *is* configured. The
  same-process determinism test only exercises the linear path.
- **God class.** `RadHACApKChargeGram`: 1,295-line header, about 202 method
  declarations, 12 constructors, 12 mode flags, and an eight-way flag cascade
  in `GetInteractionMatrixElementRaw`.
- **Exception-message string matching for control flow** in the block-PCG
  breakdown handler.
- **Thirteen environment variables** silently change numerical behaviour, and
  not all of them reach the run artifacts.

### 1.3 Worth preserving

The deterministic compensated dot product; the true-residual re-check at the
convergence exit; the sigma normalization of the fill (both born of real
incidents); the reference-coordinate moment expansion that avoids
global-monomial cancellation; and the documented "read the value before the
next fetch" discipline that prevents dangling block references.

---

## 2. The C-type magnet discrepancy

The campaign reported BDM2 versus Kelvin FEM at **1.409849%** full-vector B
RMS, with the pole edge holding **95.0905%** of the vector-error energy while
the central flat agreed to 0.110757%. Mesh refinement, quadrature order,
`leaf_size`, the Kelvin gap and the IMA count had all been excluded.

### 2.1 Pole geometry, the dominant defect

`_example5_iron()` built the pole-face chamfer with a single three-section
`netgen.occ.ThruSections`. netgen.occ does not expose OCC's `ruled` flag, so
that loft is a smooth spline surface, not the planar-faced solid
`ObjMltExtRtg` builds. The section list 34x24 to 50x40 to 50x40 has a kink at
z = 13 mm that a smooth loft cannot reproduce, so it overshoots:

| z (mm) | Actual area | Ruled reference | Ratio |
|---|---|---|---|
| 5.0 (gap-facing pole face) | 886.6 mm2 | 816.0 | **+8.7%** |
| 8.0 (chamfer mid) | 1471.6 mm2 | 1200.0 | **+22.6%** |
| 13.0 | 1999.4 mm2 | 2000.0 | 1.000 |
| 18.0 (should be constant) | 2146.2 mm2 | 2000.0 | **+7.3%** |

Pole volume +7.65%; bounding box 56.93 x 46.93 mm instead of 50 x 40. All of
it silently.

The chamfer is genuinely in the original: `RADIA_Example05.py` sets
`chamfer = 8` and builds `k1`/`k2`/`k3` for `ObjMltExtRtg`, and the reference
harness sets `ex.chamfer = 8`.

**Field impact, geometry isolated.** Both geometries meshed with the same
mesher at the same 25 mm size, same coil, material and solver, order 1:

- full-vector B RMS difference **7.2429%**
- maximum local difference **22.1482%** at (30, 0, 0) mm, which is 5 mm
  outside the pole edge, in the fringe
- magnitude of B at the origin, -0.5224%

The wrong CAD also loses iron when hex-meshed: CAD 1,451,442 against mesh
1,430,863 mm3, -1.42%, because hexes chord the convex spline bulge. The
corrected CAD is preserved exactly at every refinement.

**Why a shared wrong geometry still separates the two routes.** Both routes
consume the same OCC compound, so the defect alone cannot produce a
cross-route difference. But the defective surface is *curved*, and the two
routes discretize it differently: the FEM calls `mesh.Curve(order)` and
follows the spline closely, while BDM2 uses affine facets. The disagreement
therefore localizes exactly where the error was observed, and the measured
BDM2 refinement correction was anti-correlated (cosine -0.920220) with the
FEM-required correction, which is the signature of two representations
converging to different effective geometries rather than a resolution
deficit. A correctly planar chamfer removes this mechanism entirely.

### 2.2 B-H interpolation, a direct cross-route difference

Both routes are handed the same `get_esrf_bh_table(5)` and interpolate it
differently:

| Route | Representation |
|---|---|
| FEM | `ng.BSpline(2, ...)`, piecewise linear |
| HDiv-MMM | `scipy PchipInterpolator`, monotone cubic |

They agree at the table nodes to 1.7e-16 and diverge between them. Against
the analytic `MatSatIsoFrm` truth over the operating knee, B = 1.7 to 2.1 T:

| Route | Max relative error | RMS |
|---|---|---|
| FEM (piecewise linear) | **1.460e-3** | 4.285e-4 |
| HDiv-MMM (PCHIP) | 5.778e-5 | 1.603e-5 |

The FEM material law is 25 times less accurate, and only **20 of 221** table
points cover the entire B = 1.7 to 2.0 T knee, because the sampling is
geometric in mu0*H across eleven decades. Differential permeability differs
by up to 70% near H = 339 A/m.

Unlike the geometry, this is a *direct* cross-route difference: the two
routes literally solve with different material laws.

### 2.3 Coil and CoilBuilder, clean

| Item | Value | Original |
|---|---|---|
| Current | -2000.0 A | `ex.current = -2000` |
| Closure | `is_closed=True`, gap 0.0 mm | - |
| radius 22.5 / width 35 | r1 = 5, r2 = 40 mm | consistent with `ObjRaceTrk` |

Both routes evaluate the same `CoilBuilder.to_radia` solid-current field
exactly inside the iron, so the source cannot cause a cross-route difference.

---

## 3. Cubit-first CAD and the HEX/TET cross-check

### 3.1 Policy

Adopted in CLAUDE.md: **Cubit (ACIS) over build123d over netgen.occ** for
authoring; netgen.occ is I/O and meshing only; **Cubit into `export netgen`
into `.vol` is the lab's main mesh path**. Python never calls Cubit; a
journal emits the artifacts and Python reads them, which preserves the
Layer-4 separation policy.

Three measurements drove this past a preference:

1. `ThruSections` has no `ruled` option and silently produced the +7.65%
   solid above.
2. Cubit/ACIS reproduces the intended geometry exactly: total volume
   **1,446,095.333333 mm3** against the analytic 1,446,095.333333, relative
   error 0.0, bounding box exact, all ten faces planar. Through
   `export netgen` into NGSolve it is +1.5e-12% with `check-vol` reporting
   zero errors.
3. `netgen.occ.OCCGeometry(<step>)` returns **1 solid from a 14-solid STEP**,
   dropping thirteen without a word. The file carries 14
   `MANIFOLD_SOLID_BREP` entities and OCP's reader returns all 14. A
   Cubit-STEP-netgen route is therefore unsafe; Cubit into `.vol` is not.

Also recorded: Cubit's APREPRO evaluates brace expressions inside `#` comment
lines.

### 3.2 The HEX/TET tension

A C-shaped yoke needs decomposition to hex-mesh; a united solid fails `auto`,
`sweep`, `webcut_cyl_auto` and `polyhedron` alike. But decomposition leaves
same-material internal interfaces, and the exporter writes those with
DomainIn equal to DomainOut equal to 1. It has no omission logic, contrary to
`build_esrf_cubit_hdiv_iron`'s docstring. TET needs no decomposition and
passes `check-vol` cleanly.

A chamfered pole cannot be meshed into affine hexes. The mapping residual
does not converge under refinement: 0.745 at 25 mm, 0.666 at 12.5 mm, 0.405
at 6.25 mm, against a 1e-10 gate, because cells straddling the 45-degree face
are trapezoidal prisms by construction. Removing the chamfer gives zero
non-affine cells at 2.2e-14, which is the comparison floor rather than a
production option while the benchmark keeps its chamfer.

### 3.3 Cross-check result

Corrected Cubit/ACIS quarter yoke, `image='+x-z'`, order 1:

| Case | ne | ndof | B at origin | B RMS vs reference |
|---|---|---|---|---|
| HEX 8.75 mm (reference) | 372 | 9,232 | 0.249965 T | - |
| HEX 12.5 mm | 176 | 4,432 | 0.249908 T | 0.1898% |
| TET 12.5 mm | 1,314 | 8,208 | 0.249708 T | **0.4788%** |

At comparable DoF the two families agree to **0.48% RMS and 0.10% at the
centre**. The HEX internal-interface warning therefore does not corrupt the
result, which lowers the priority of the exporter change.

Finer TET meshes break down: `p^T A p` is -1.90e5 at 8.75 mm and -7.51e8 at
7.0 mm, four orders worse under refinement. This is an operator defect, not
slow convergence, and it was the B2 guard that separated the two. The TET
lane currently cannot go below 12.5 mm.

---

## 4. Regression found and fixed

The `image='+x-z'` quarter route was **completely stopped on main**. The
image-folded self-energy gate accepted only an exact zero diagonal, which no
finite-precision cancellation reaches: a charge on an antisymmetric mirror
plane has `G_self + sign*G_refl(a,a)` equal to zero in exact arithmetic, but
in double precision the residue lands within roundoff of zero with a
geometry-dependent sign. All five example-5 quarter meshes, TET at 12.5, 8.75
and 7.0 mm and HEX at 12.5 and 8.75 mm, failed at about -5e-7 against an
order-one diagonal scale. The existing roundoff test covered only the
positive side of the same cancellation, so it stayed green.

Fixed in `4b30479b6`: the accept band is now 1e-12 of the largest diagonal
actually present, the strict test is kept when no image folding can cause a
cancellation, and a new test locks the negative side with an annihilated
diagonal of -2.8e-17 against a 0.19 companion, verified negative rather than
merely small.

---

## 5. Corrections to earlier statements in this review

Recorded because the record should show what was wrong, not only what was
found.

1. **"The non-convergence contract is missing" was wrong.** `_solve.py`
   enforces it and says so. The real gap was only the missing `p^T A p`
   breakdown detection.
2. **"`ComputeChargeSigma`'s permissive guard is a silent fallback" was
   wrong.** It is deliberate and load-bearing for antisymmetric mirror-plane
   charges. The first fix broke `test_hdiv_vim_2d_ima.py`; the corrected
   version distinguishes image-folded from non-folded operators.
3. **"All-affine hex meshing of the chamfer is impossible" was asserted
   before it was measured.** Parallelepipeds need not be axis-aligned, so the
   assertion was unsound as stated. The refinement measurement in 3.2 is the
   evidence that actually supports the conclusion.
4. **"The spline surface causes the non-affine cells" was wrong.** The
   chamfer taper causes them in both geometries, 4 of 88 and 4 of 92.
5. **A comparison script written for this review used a silent fallback** and
   reported a 0.000000% difference by measuring the coil field alone, which
   is exactly the failure mode the No-Fallbacks policy exists to prevent.
6. **The MCP disconnection was misdiagnosed twice**, first as a possible
   `export netgen` crash and then as a healthy server behind a dead
   transport. The actual cause was an expired Cubit license, which the
   daemon-free checks could not see.

---

## 6. Open items

| # | Item | State |
|---|---|---|
| a | TET operator indefiniteness under refinement | Open; worsens four orders from 8.75 to 7.0 mm |
| b | Unify the B-H interpolant, then re-measure BDM2 against FEM on corrected geometry | Blocked on the untracked corpus |
| c | Exporter omission of same-material internal interfaces | Deprioritized; no measured harm |
| d | Determinism of the nonlinear path | Reported, not fixed |
| e | `RadHACApKChargeGram` decomposition | Reported, not fixed |

The example-5 corpus, `src/radia/esrf_examples.py` and
`validation_test/esrf_radia_examples_cubit/`, is untracked and exists only in
the shared checkout, which blocks item b. It is being committed separately.
