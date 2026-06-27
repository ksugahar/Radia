# Why scalar BIE-SIBC + ESIM is the right combination

**Audience.** Reviewers / readers who ask "why scalar BEM, why not vector
EFIE / MFIE, and why not just stay in FEM?".

**TL;DR.** For an MQS induction-heating workpiece with a Leontovich
surface impedance, the *minimal* compatible variational triple is

```
   scalar magnetic potential phi      <-- only 1 DOF / surface node
+  Lagrange P2 basis on Tri6 surface  <-- O(h^3) approximation
+  isoparametric curved Tri6 geometry <-- O(h^3) geometric error
```

and the **three orders match**.  No Piola transformation, no edge-element
ABI translation, no flat-element area loss, no air mesh.  This document
explains why every adjacent choice (vector BEM-A, flat Tri3, FEM-Kelvin)
is either redundantly more expensive or geometrically lossy at the same
target accuracy.

---

## 1. The four combinations on the table

| # | Workpiece formulation | Coil source path | DOFs | Geometric truncation |
|---|---|---|---|---|
| **A** | **scalar BIE-SIBC on phi, Tri6 P2** | PEEC filament + Biot-Savart | ndof_surface | none (BEM) |
| B | Vector BEM-A on J_s, RWG | PEEC filament + Biot-Savart | 2 x nedge | none (BEM) |
| C | HCurl A on volume + Robin SIBC | PEEC filament line-integral RHS | ndof_volume | **Kelvin transformation** |
| D | HCurl A-V compound (coil volume + workpiece) | volumetric coil mesh | ndof_volume + n_coil | **Kelvin transformation** |

All four are implemented in Radia ([`calc_inductance.py`](../../src/radia/panels/calc_inductance.py),
[`calc_fem_kelvin.py`](../../src/radia/panels/calc_fem_kelvin.py),
[`calc_fem_coilmesh.py`](../../src/radia/panels/calc_fem_coilmesh.py)) — A and B
both via `calc_inductance.py --coil-solver {peec,bem-a}`.

The internal-consistency [§3 of CROSS_VALIDATION.md](CROSS_VALIDATION.md)
shows A ≡ C ≡ D within the SIBC's own validity window.  This document
explains why **path A is the right one to publish**.

---

## 2. The physics permits a scalar — and only a scalar

Inside the MQS air domain, all current is on the coil (Biot-Savart
known) and the workpiece (Leontovich SIBC encapsulates it).  The
air-side magnetic field decomposes as

$$
\mathbf{H}(\mathbf{r}) = -\nabla \varphi(\mathbf{r})
+ \mathbf{H}_{\mathrm{inc}}(\mathbf{r}),
\qquad \mathbf{r} \in \Omega_{\mathrm{air}},
$$

where `H_inc` is the coil's Biot-Savart contribution (curl-free
in air, exact) and `phi` is the workpiece-scattered scalar potential.
Single-valuedness of `phi` requires only that `Omega_air` is
simply connected; this is satisfied for a closed workpiece + closed
coil topology (no flux linkage between air-side circuits).

**The vector A field is redundant here.**  Once one accepts the
SIBC model (which already presumes the workpiece is "thin-skin-thin"
relative to its tangential extent), the air-side problem has
*scalar* number of DOFs per spatial point, not vector.  Solving
for a 3-component A field on Tri6 surface DOFs is **6x as many DOFs
as the underlying physics needs**, and those extra DOFs are
constrained by `n·curl A = B_n` and gauge conditions to give
exactly the same `phi` after reconstruction.

This is the classical "EFIE vs MFIE vs single-layer/double-layer"
distinction: for closed surfaces in source-free regions, the
scalar Laplace BIE (Green's function of `nabla^2`) is the minimal
formulation.  Calderón calculus references: Sauter-Schwab 2011
*Boundary Element Methods*, §3.4.

---

## 3. Why Tri6 P2 pairs naturally with scalar phi (and not with RWG)

### 3.1 The scalar P2 path (path A)

Tri6 = 6 isoparametric Lagrange nodes (3 vertices + 3 edge midpoints).
The same 6 nodes carry:

- **Geometry** : `x(u, v) = sum_i x_i * N_i^{P2}(u, v)` — quadratic
  approximation of the workpiece surface.  Curved.  Sourced from the
  Cubit `export netgen "f.vol" order 2` companion JSON via
  [`bem_sibc_solver.py`](../../src/radia/bem_sibc_solver.py).
- **Basis** : `phi(u, v) = sum_i phi_i * N_i^{P2}(u, v)` — scalar
  Lagrange P2 finite element.
- **DOFs** : `phi_i` at each of the 6 nodes per element.

Result: **the basis function and the geometry use the same 6 nodes
with the same shape functions** — perfect isoparametric setup.

Quadrature:

- Singular pairs (`i = j` or shared vertex/edge): Sauter-Schwab Duffy
  transformation at degree 6.
- Regular pairs: Gauss-Legendre on the parametric triangle at degree 7.

Error rates (for a smooth surface and smooth solution):

| Quantity | Rate |
|---|---|
| `\|\|\varphi - \varphi_h\|\|_{L^2}` | O(h^3) (P2 basis order) |
| `\|\|x - x_h\|\|_\infty` on Gamma | O(h^3) (Tri6 geometric order) |
| End-to-end `P_wp` error | O(h^3) (matching rates) |

This **rate match is the unique feature of the scalar P2 path**.
The next three subsections show why no other reasonable
discretisation gets all three rates aligned.

### 3.2 The vector BEM-A path (path B): Piola transformation required

Vector unknown `J_s = n × H_t` lives in `H(div_Γ; Γ)` — the
surface div-conforming space.  The canonical basis is RWG (= RT_0
on triangles).  Under the parametric mapping `x: T_hat → T`,
contravariant Piola transformation is needed:

$$
\mathbf{J}_s(x) = \frac{1}{|\det J|}\,J\,\hat{\mathbf{J}}_s(\hat{x})
$$

The Piola transformation:

1. **Introduces `det J`** — for a curved Tri6 mapping `det J` varies
   across the element.  Quadrature must integrate
   `1/det J × kernel × Piola J × Piola J` exactly, demanding a
   higher quadrature degree.
2. **Couples geometry to basis non-trivially** — the basis "twists"
   as the element curves, and the discrete divergence `div_Γ J_h` on
   the curved geometry is no longer piecewise-constant per element.
   Loss of the conservation property that makes RWG attractive on
   flat triangles.
3. **For RT_0** (the only widely-implemented RWG basis): degree-1
   vector approximation that gives `O(h)` in `||J - J_h||` on a flat
   geometry, `O(h^2)` only with quadratic curving plus a
   higher-order RT_k basis (uncommon in standard BEM packages).

So the Tri6 + RWG path either:

- (i) gives `O(h)` accuracy in `J_s` regardless of geometric order
  — flat-class accuracy paid for at curved-class cost, OR
- (ii) requires raising to RT_1 / RT_2 surface elements (rare
  off-the-shelf), at which point the DOF count balloons further.

The Radia BEM-A path implements (i) via [NGSolve `ngsbem`](https://docu.ngsolve.org/latest/i-tutorials/unit-12.html)
HDivSurface with Weggler low-frequency stabilisation; it works for
linear-SIBC coil source but uses **2× the surface DOFs of the
scalar BIE** without delivering matched `O(h^3)`.

### 3.3 Flat Tri3 + P1 path: geometry caps accuracy

The cheapest variant — Tri3 flat geometry + Lagrange P1 phi — has

- basis-order error : `O(h^2)`
- geometric error   : `O(h^2)` (chordal flat-triangle on a curved surface)

Both rates are quadratic, *matched*, so the path is *internally
consistent*.  But: on a curved workpiece (cylinder, gear root,
fillet) the flat-triangle area error converges as `O(h^2)` only in
the asymptotic regime.  For engineering-mesh sizes (h/R ≈ 0.05) the
constant pre-factor often dominates and `P_wp` can be off by 1-2 %
just from surface area.  Section [§6 of CROSS_VALIDATION.md](CROSS_VALIDATION.md)
shows that on a smooth Cu cylinder this effect is small (+0.3 % from
flat→curved), but on a high-curvature workpiece (gear teeth) the
flat-Tri3 path is provably worse.

The scalar P2 path **captures both `O(h^3)` rates with the same DOF
count** as flat P1 on a refined mesh (h_P2 ≈ 2 × h_P1) — i.e. P2 is
typically cheaper than P1 refinement at a fixed accuracy target.

### 3.4 FEM-Kelvin path (C): same DOF order on volume, plus truncation

Going volumetric with HCurl A on a Kelvin-transformed exterior
domain:

- **DOFs** = O(N_volume) on the workpiece + Kelvin sphere interior.
  Typically `10^4`-`10^5` DOFs vs `10^2`-`10^3` for the BEM surface
  path.
- **Kelvin transformation** is exact for the homogeneous Laplace
  exterior, but the discretised Kelvin mesh introduces its own
  truncation error (geometric and basis).
- **Curvature handling** for HCurl: NGSolve `Curve(p)` works on
  tetrahedra, so the workpiece SURFACE gets curved correctly via
  the boundary projection — same `O(h^3)` rate as the scalar BIE on
  the workpiece surface.  But the workpiece *interior* mesh and the
  Kelvin exterior consume the bulk of the DOFs.

So path C is mathematically clean and, like Hollaus et al., it is
a **volume FEM** that needs an (unbounded) air mesh.  NB the
potential differs: our path C discretises the vector potential
**HCurl A** on the Kelvin-transformed exterior, whereas Hollaus
et al. use a **magnetic scalar potential** ("A Nonlinear Effective
Surface Impedance in a Magnetic Scalar Potential Formulation",
IEEE TMag 2025) — closer in spirit to our scalar BIE, but still
volumetric.  Either way the volume FEM pays a **100x DOF count**
for the same surface-side accuracy.  Worth it only when the air
mesh genuinely needs to resolve something the BEM cannot — e.g.,
a magnetic core or a thin-gap structure inside the air.  For pure
IH (coil + workpiece + vacuum), it is wasteful.

### 3.5 Volumetric coil A-V (D): a separate problem

Path D is "no SIBC on the coil side either" — useful as a reference
solver (workpiece SIBC, coil volumetric) and as a way to measure
the Joule loss on the coil itself (P_coil).  It is **not** the
right comparison for "BEM vs FEM ESIM" because it changes the coil
representation too.  It serves as a high-fidelity cross-check
([§3 of CROSS_VALIDATION.md](CROSS_VALIDATION.md): at 10 kHz the
agreement is <1 %; at higher frequency the coil mesh resolution
becomes the limiting factor and not the workpiece formulation).

---

## 4. ESIM is invariant under the workpiece-formulation choice

The cell problem itself

$$
(\rho / r)\,\partial_r[r \partial_r H] + j \omega \mu(|H|)\,H = 0
$$

is a 1-D PDE solved per surface DOF (or per scalar `H_t_rms`).
**It does not care whether the outer solve is BEM or FEM.**  All
four paths in §1 plug the same cell solver into the same Karl
fixed-point.  The cell problem returns `Z_s(|H_t|)`, which is then
inserted into:

- BIE: row scaling of `(1/2 M - DL + jω/Z_s · SL M⁻¹ K)`
  ([`bem_sibc_solver.py:401`](../../src/radia/bem_sibc_solver.py#L401))
- HCurl: Robin coefficient `jω/Z_s` on the workpiece BND
  ([`calc_fem_kelvin.py:565`](../../src/radia/panels/calc_fem_kelvin.py#L565))

The **scalar BIE is the cheapest carrier**; the cell-problem
calculation is the dominant cost only for per-element Karl on very
large surface DOFs (>5000), where it can match the BIE solve time.

## 5. Headline numerical impact (the marketing line)

The IGTE 2026 digest reports:

> For a steel cylinder driven at 50 kHz, `I_port = 100 A` (through
> the BH knee), the per-element scalar BIE-SIBC reports
> `P_wp = 18.75 W`, **38.5 % below** the scalar-Z_s formulation's
> `30.51 W`, because local saturation reduces `Z_s` at hot-spot DOFs
> that the uniform scalar model averages away.

This gap is the dense-sweep source-of-truth value for the digest
representative cell.  Similar per-element vs scalar-Z_s comparisons on
path C or D are useful cross-checks, but the digest headline should cite
the committed 108-case scalar-BIE sweep.  The formulation argument remains:

- Path A (scalar BIE) makes per-element Z_s the **natural row-wise
  scaling** of an existing BIE matrix — 3-line code change in the
  assembler.
- Path C (HCurl Kelvin) requires the Robin coefficient `jω/Z_s` to
  become a `CoefficientFunction` of the per-BND-DOF Z_s array — a
  full re-assembly of the Robin term per Karl iteration (vs
  diagonal-only re-scaling for BIE).
- Path D (HCurl A-V) inherits path C's costs *and* re-assembles the
  coil-side conductivity term per iter.

**Per-iteration Karl cost ratio (production gapped-torus benchmark)**:

| Path | Per-iter cost (re-assembly + solve) |
|---|---|
| A (scalar BIE, dense LU) | 0.2 s (166 DOFs) — diagonal re-scaling only |
| C (FEM-Kelvin, pardiso)  | 12 s (12k DOFs) — full re-assembly |
| D (FEM-coilmesh, pardiso) | 25 s (38k DOFs) — full re-assembly |

The scalar BIE-SIBC path is **50-100x cheaper per Karl iteration**
on top of being cheaper to *assemble* in the first place.

---

## 6. Summary — what to claim in publication

The IGTE digest abstract states three contributions; this document
unpacks WHY they are not orthogonal but a single tightly-coupled
choice:

1. **Scalar BIE on the workpiece surface only** — the minimal MQS
   formulation when air is source-free.  Eliminates the air mesh
   (FEM Kelvin) and the redundant DOFs of vector BEM (EFIE/MFIE).
2. **Curved isoparametric Tri6 P2** — the geometric and basis
   orders **match** at `O(h^3)`.  No Piola, no edge-element
   ABI translation, no flat-element area loss.
3. **Element-by-element nonlinear Z_s** — natural diagonal scaling
   of the BIE assembly; 3-line code change vs full re-assembly in
   FEM paths.  Captures saturation patterns that scalar mesh-RMS
   cannot resolve.

Each of (1)-(3) is implementable independently in other frameworks,
but only their combination on **path A** is *uniquely cheap and
geometrically lossless at O(h^3)*.  Path A is the SUBJECT of the
paper; the BEM-A / FEM-Kelvin / FEM-coilmesh paths are
*consistency-check references*, **not competitors to be
benchmarked against**.

---

## 7. Suggested reviewer Q&A

| Q | A |
|---|---|
| "Why not vector EFIE/MFIE for completeness?" | MQS air domain is curl-free; vector BEM solves 6× the DOFs of the underlying scalar physics.  We use vector BEM-A only as a coil-source variant when the coil topology requires it (§§ R_MISMATCH_PEEC_VS_BEMA.md). |
| "Why scalar BIE in 2026 — isn't this old?" | The contribution is not the BIE.  It is the demonstration that scalar BIE + curved Tri6 + per-element ESIM is the **uniquely well-matched** discretisation for nonlinear SIBC IH problems, with measurable downstream cost (this document §5). |
| "Don't you lose generality by assuming closed surfaces?" | Yes; the scalar BIE requires `Omega_air` simply connected.  This is satisfied for closed-conductor IH workpieces but not for, e.g., open-strip transformer windings.  For those cases path B (BEM-A) or path C (FEM-Kelvin) remains appropriate.  Documented as a limitation in §6 of [`MATHEMATICAL_ANALYSIS.md`](MATHEMATICAL_ANALYSIS.md). |
| "Why is the FEM-coilmesh path in your figures if you say BEM is enough?" | As a high-fidelity reference, not a competitor.  FEM-coilmesh resolves coil-side ohmic loss explicitly; the BEM path uses Dowell SIBC on the coil.  Both agree on `P_wp` to <2 % at 10 kHz; the BEM path is 50× cheaper per Karl iteration. |
| "Can per-element Z_s be done in vector BEM-A too?" | Yes, in principle, but the row-scaling is per-edge (RWG) rather than per-node (P2) and the row pattern interacts with the Piola transformation; the analytical simplicity of "one Z_s per surface DOF" is specific to scalar P2.  This is exactly the matching-orders argument of §3. |

---

## 8. Cross-references

- [`MATHEMATICAL_ANALYSIS.md`](MATHEMATICAL_ANALYSIS.md) §§ 1-4 — scalar BIE
  derivation, weak form, curvature handling.
- [`IMPLEMENTATION.md`](IMPLEMENTATION.md) §§ 3-4 — Karl loop, Robin BC
  implementation per path.
- [`CROSS_VALIDATION.md`](CROSS_VALIDATION.md) § 3 — three-path numerical
  agreement.
- [`R_MISMATCH_PEEC_VS_BEMA.md`](R_MISMATCH_PEEC_VS_BEMA.md) — when path B
  (vector BEM-A) is preferred over path A's PEEC filament coil
  (n_peri-free, but R disagreement).

---

**Document version**: 2026-05-30 (radia v4.67.0+ dense-sweep baseline).
