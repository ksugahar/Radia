# Coarse-Mesh Accuracy of the Kelvin Transformation as a DtN-Spectrum Property

Kameari's classic demonstration shows, by mesh **refinement**, that the Kelvin
transformation already gives good accuracy on a relatively coarse mesh. This
document reframes that *empirical* observation as a *property of the
Dirichlet-to-Neumann (DtN) matrix* measured across mesh sizes — the accuracy
becomes a **spectral** fact, readable off the boundary operator **before any
field is solved**. It is the academic/reference companion to
[KELVIN_TRANSFORMATION.md](KELVIN_TRANSFORMATION.md); the runnable layer is the
`radia-mcp` tool `dtn_coarse_mesh(topic=...)` and the code/tests/example listed
in §10.

> **Venue.** Prepared as the supporting analysis for the IEEJ Static
> Apparatus / Rotating Machinery joint technical meeting (静止器・回転機合同研究会),
> 2026-06. The 2-D circle-inversion result (§7) is the static-apparatus /
> rotating-machine cross-section case.

---

## 1. The question (Kameari, reframed)

Kameari demonstrated the Kelvin transformation's accuracy the empirical way:
take a problem with a known answer, solve it on a sequence of meshes, and show
the error is already small on a relatively **coarse** mesh (and merely polishes
as you refine). The question answered here:

> Can that coarse-mesh accuracy be stated directly from the **properties of the
> DtN matrix at various mesh sizes**, instead of from a refinement experiment?

Yes. The accuracy is a **spectral** fact about the discrete DtN operator,
visible without ever solving the field problem, and it **isolates the
open-boundary error from the interior FEM error** — which a field-only
refinement study conflates.

---

## 2. One operator behind every open-boundary method

Truncate the unbounded exterior at a surface `Γ` and impose the exact
transparent condition there:

```
    ∂u/∂n |_Γ  =  Λ_ext u |_Γ ,
```

where `Λ_ext` is the exterior Steklov–Poincaré (Dirichlet-to-Neumann) operator:
it returns, for any boundary trace `u|_Γ`, the outward normal derivative of the
unique exterior harmonic function that decays at infinity. Every practical
open-boundary technique realises some discrete operator `S_h` on `Γ` targeting
this single `Λ_ext`, in two distinct fidelity classes:

| Method                 | How it realises an `S_h` on `Γ`        | Fidelity            |
|------------------------|----------------------------------------|---------------------|
| Kelvin transformation  | FEM on the inverted ball; trace on `Γ` | → `Λ_ext` as h→0    |
| BEM coupling           | `S_h = Λ_h = V⁻¹(−½M+K)` (boundary)    | → `Λ_ext` as h→0    |
| PML                    | absorbing layer; trace on `Γ`          | fixed-error model   |
| Infinite elements      | decaying shape functions; trace on `Γ` | → `Λ_ext` (basis)   |
| Asymptotic Robin BC    | `S_h = −(1/R) I`                       | exact `n=0` only    |

Kelvin / BEM / (rich-enough) infinite elements are **discretizations** that
converge to the exact `Λ_ext` for every degree as `h→0`; PML and the asymptotic
Robin BC are **model surrogates** with a fixed, mesh-independent modal error.
This document concerns the first class: there the open-BC error is governed by
how well the discrete `S_h` reproduces `Λ_ext` **on the modes the solution
contains**.

**The Kelvin transformation IS a spherical inversion** (`r ↦ R²/r` about a
centre), so its `Γ` is necessarily the inversion sphere of radius `R` — there is
no non-spherical Kelvin. Hence the spherical eigenvalue ladder below is not a
special case for Kelvin; it is the **complete, general** spectral structure of
any Kelvin coupling. (The interior physical domain may be any shape; it sits
inside the Kelvin sphere and the coupling is on that sphere.) BEM / PML /
infinite elements admit any closed `Γ` — on a non-sphere the relevant spectrum
is that surface's own Steklov spectrum.

---

## 3. The continuous DtN spectrum is known in closed form

On a sphere of radius `R` the spherical harmonics `Y_n^m` diagonalise `Λ_ext`:

```
    Λ_ext Y_n^m = −(n+1)/R · Y_n^m ,     m = −n … n   (multiplicity 2n+1)

    n = 0 monopole   : λ₀ = −1/R
    n = 1 dipole     : λ₁ = −2/R
    n = 2 quadrupole : λ₂ = −3/R
    n = 3 octupole   : λ₃ = −4/R
```

These eigenvalues are properties of the **continuous** operator — no mesh, no
`h`. They are the yardstick: a discrete DtN matrix `Λ_h` is "good for mode `n`"
exactly when its corresponding eigenvalue matches `−(n+1)/R`. (In 2-D, the
circle inversion has no conformal prefactor and the ladder is `−n/R`; see §7.)

---

## 4. The matrix property that explains coarse-mesh accuracy

Assemble the dense DtN matrix `Λ_h = V⁻¹(−½M+K)` on a sphere at several mesh
sizes and read off its eigenvalues. Three facts fall out.

**Measured — BEM `Λ_h` spectrum, sphere R=1, SurfaceL2 order=1, intorder=10.**
The sphere mesh floors at `ndof=336` (every `maxh ≥ 0.5` gives that coarsest
mesh), so the three columns are the three genuinely distinct refinement levels:

```
   per-degree relative eigenvalue error  |λ_h,n − λ_n| / |λ_n|

   degree n      | ndof=336    ndof=564    ndof=924      refines →
                 | (maxh 0.5)  (maxh 0.4)  (maxh 0.3)
   --------------|--------------------------------------------------
    0  monopole  | 2.6e-04     9.9e-05     3.6e-05
    1  dipole    | 7.4e-04     2.9e-04     1.0e-04
    2  quadrupole| 2.1e-03     8.0e-04     3.0e-04
    3  octupole  | 5.2e-03     2.1e-03     7.8e-04
    4            | 1.1e-02     4.6e-03     1.8e-03
   --------------|--------------------------------------------------
   accurate band | n ≤ 2       n ≤ 4       n ≤ 4
   (rel_err<0.5%)|
```

1. **The low modes are already accurate on the coarsest mesh.** On the coarsest
   admissible sphere mesh the `n = 0,1,2` eigenvalues match `−(n+1)/R` to
   `< ~0.25 %` (the quadrupole is the worst at 0.21 %), because the
   eigenfunctions (low spherical harmonics) are smooth: a coarse surface mesh
   interpolates them with tiny error, and the Galerkin eigenvalue (which for an
   analytic eigenfunction superconverges as the **squared** L2 trace error)
   inherits an even smaller error. **This is the coarse-mesh accuracy itself —
   visible in the matrix, before any solve.**

2. **At every mesh size the error is ordered by degree `n`** (the spectral
   signature). `rel_err` increases smoothly and monotonically with `n` — not a
   cliff. A fixed mesh is "accurate up to some low degree, progressively worse
   above."

3. **Refinement lowers every mode steeply and widens the accurate band.** Each
   refinement cuts the per-degree error `≈ ×2.6` (measured rate `p ≈ 3.9`, an
   order-1 Galerkin eigenvalue superconvergence `~O(h⁴)`), pushing the
   resolution limit up so more degrees fall below any fixed tolerance. But the
   low modes were already below engineering tolerance on the coarse mesh —
   refining buys **higher modes**, not materially better low ones.

The honest distinction (not "the low modes are mesh-independent"): the low-mode
error *does* fall under refinement, but it starts from such a small floor on the
coarse mesh that there is nothing of engineering value to gain there. What
refinement adds is **bandwidth** — accuracy for higher multipoles.

---

## 5. Why this matches Kameari's observation

For a **compact** source inside `Γ` the exterior field is a multipole series
whose `n`-th term decays like `r^{−(n+1)}`. At a truncation radius a few times
the source size the boundary data is dominated by the lowest harmonics.
Therefore the open-BC part of the error depends on `(S_h − Λ_ext)` only through
its **low-degree block** — the very block that is already accurate on a coarse
mesh (fact 1) and is the most accurate part of the spectrum (fact 2). Kameari's
field-refinement curve is the per-degree convergences of fact 3 superposed and
weighted by the source's multipole content; because that content is
low-degree-dominated and those modes are already accurate coarse, the open-BC
contribution is already small on the coarse mesh and merely polishes.

This is an **explanation** (implication), not a strict equivalence; it rests on
three premises, all checkable from the spectrum:

- **(i)** the boundary data is low-degree dominated (truncation a few
  source-radii out; no sharply-structured source sitting near `Γ`);
- **(ii)** the open-BC error is not the bottleneck — the total field error also
  carries the **interior FEM** discretisation error, typically the larger term
  on a coarse mesh (≈ 5 % L2 from the interior solve while the dipole DtN
  eigenvalue is already 0.07 %). The spectrum isolates the boundary
  contribution from this interior error (§6);
- **(iii)** for the Kelvin transform specifically, that its `S_h` shares the
  low-mode fidelity — now **measured** directly (§6 bridge).

With those premises the refinement experiment and the spectral statement are two
views of the same continuous operator. The spectral one is stronger: it says
the answer is good **before** you refine, says exactly which degrees a given
mesh can be trusted for, and isolates the open-BC error from the interior FEM
error — without solving the field problem at all.

---

## 6. The Kelvin bridge: measured, not just argued

The Kelvin inversion `r ↦ R²/r'` maps the slowly-decaying low-`n` exterior modes
to **bounded solid harmonics** on the Kelvin ball. With the 3-D conformal
weight `(R/r')^{d−2} = R/r'`:

```
    u*_n(r') = (R/r')·u_n(R²/r') = (R/r')·(R²/r')^{−(n+1)} Y_n  ∝  r'^{n} Y_n
```

so the exterior `r^{−(n+1)} Y_n` becomes the **solid harmonic** `r'^{n} Y_n` on
the ball — a **degree-`n` polynomial**. Order-`p` Lagrange FEM represents a
degree-`n` polynomial **exactly iff `p ≥ n`** (in the reference space; on the
curved 3-D sphere the realized accuracy then floors at the geometry + conformal-
weight error, ~5–6 digits — exactly Kameari's observation, not machine-exact).
That is the mechanism, measured
directly by `kelvin_dtn_eigenvalue` (a volume-FEM solve on the inverted ball,
reading `λ_eff = −1/R − ∫_Ω|∇u*|² / ∮_Γ u*²` and comparing to `−(n+1)/R`):

```
   Kelvin closure effective DtN, rel_err vs −(n+1)/R  (volume FEM, R=1):

   mode          inverts to | order 1            order 2     order 3
   --------------------------|---------------------------------------
   dipole    n=1   linear    | 0.5% -> 0.03%      --          --
                             | (maxh 0.6->0.25, converges ~O(h³))
   quadrupole n=2  quadratic | 18%  (maxh 0.4)    0.03%       --
   octupole  n=3   cubic     | 40%  (maxh 0.4)    --          4e-5
```

The **dominant dipole inverts to a linear field**, so even order-1 FEM nails its
DtN eigenvalue on the coarsest mesh — the residual is only **geometry**
(curved-sphere) error, converging fast; nothing about the dipole needs refining.
The quadrupole needs order ≥ 2, the octupole order ≥ 3. So the Kelvin closure is
coarse-mesh accurate for exactly the low modes a compact source radiates — a
**measured** fact, matching the BEM `Λ_h` table via the shared continuous
`Λ_ext`.

**Two complementary mechanisms, one conclusion:**

- **BEM `Λ_h` (boundary):** all low surface harmonics are smooth, so a coarse
  surface mesh resolves them; error grows **smoothly** with degree `n`.
- **Kelvin (volume):** mode `n` inverts to a degree-`n` polynomial; the closure
  is exact up to FEM order — a **sharp threshold** at `n = p`.

### The floor IS geometry — proof by raising only the Curve order

The 5–6 digit floor that remains at `order ≥ n` is the **curved-sphere geometry**,
not the multipole or the method. Proof: hold the FE order `p ≥ n` **and** the
coarse mesh fixed, and raise ONLY the isoparametric geometry order `k`
(`mesh.Curve(k)`):

```
   Curve (geometry) order k :   1 (flat)    2          3
   n=2 (p=3) rel_err        :   1.33e-2     3.8e-4     1.30e-5
   n=3 (p=4) rel_err        :   1.48e-2     4.4e-4     3.42e-5
```

Raising only the geometry order from `k=1` (flat polyhedron, ~1 % faceting) to
`k=3` drops the error **~1000×** to the 5–6 digit floor — the polynomial image and
the FE order untouched. A flat truncation is ~1 % off; an isoparametric (curved)
one reaches 5–6 digits. (Past `k ≥ 3` it plateaus ~1e-5: the residual
conformal-weight quadrature / energy-quotient limit.) Script:
`examples/kelvin_transformation/DtN_spectrum/floor_vs_curve.py`.

### Isolating the Kelvin open-BC error from the interior FEM error

On **one shared shell mesh**, swap only the `Γ` operator — exact-DtN Robin
(`λ = −(n+1)/R`, zero open-BC error), Kelvin-DtN Robin (`λ` from
`kelvin_dtn_eigenvalue`), and the BEM-DtN Schur — so the interior FEM error
**cancels** and `‖u_method − u_exactDtN‖` is the method's **pure** open-boundary
error:

```
   interior_fem_error   ≈ 5.3e-2   (shared; open BC is EXACT here)
   kelvin_openbc_error  ≈ 1.2e-3   (ISOLATED Kelvin open-BC error)
   bem_openbc_error     ≈ 1.7e-3   (ISOLATED BEM open-BC error)
```

The Kelvin closure's own error is **≈ 0.1 %, ≈ 45× below** the ≈ 5 % interior
FEM error, even on a coarse Kelvin mesh — Kameari's coarse-mesh accuracy as a
**separated error budget** rather than a refinement curve.

### Accuracy vs the exterior (Kelvin-ball) mesh size

In the Kelvin transform the unbounded exterior **is** the Kelvin ball, so the
"exterior-region mesh size" is the Kelvin-ball `maxh`. Sweeping it with the
interior mesh held fixed, the (fixed) interior FEM error cancels and the
reported open-BC error is purely the exterior-discretisation contribution:

```
   open-BC error: 1.2e-3 → 7.5e-5  (≈ ×4/level, kelvin_ndof 58 → 301)
   interior FEM error: fixed 5.3e-2   (dominates 45× → 709×)
```

A **coarse exterior mesh already suffices** — Kameari's exterior-refinement
check, isolated: refine the exterior, watch the open-BC error converge below the
fixed interior floor. (The Kelvin-ball mesh floors at `ndof = 58` for
`maxh ≳ 0.5`, the volume analogue of the sphere-surface floor.)

---

## 7. 2-D cross-section: static apparatus / rotating machinery

The Kelvin transform is a spherical inversion, so in 2-D the truncation surface
is a **circle** and there is no conformal prefactor (offset `0`): the eigenvalue
ladder is **`−n/R`** (not `−(n+1)/R`), while the order-threshold mechanism (order
`≥ n` removes the polynomial error; the dipole inverts to a **linear** field →
order-1 coarse-accurate) is the **same** as 3-D — but with no conformal weight the
realized geometry floor is **deeper** than 3-D (~`1e-7`…`1e-9` vs the 3-D sphere's
~5–6 digits).

```
   2-D circle inversion, effective DtN rel_err vs −n/R  (R=1):
   dipole     n=1  → linear     : 0.06 %  at order 1 / coarse mesh
   quadrupole n=2  → quadratic  : needs order ≥ 2
```

This is the directly relevant case for the planar cross-section of a static
apparatus or rotating machine: the same coarse-mesh accuracy, with the
2-D ladder.

### The lab's real two-sphere periodic Kelvin (end-to-end validation)

`kelvin_twosphere_shell_dipole` solves the genuine two-offset-sphere periodic-BC
Kelvin (Nagamine convention: material `μ' = (R/r')²` in the inverted ball, GND
at the Kelvin centre, periodic identification by translation) — **not** the
single-ball effective-DtN equivalent:

```
   shell L2 vs exact dipole u = R_in² z / r³ :  ≈ 5.2e-2  (maxh 0.25, order 2)
```

This matches the analytic-DtN shell solve, so **the real Kelvin coupling = the
analytic exterior = the single-ball effective DtN** — the whole picture
confirmed on the actual implementation.

---

## 8. Using the DtN-spectrum view

1. **Sizing the air box / truncation radius without a convergence sweep.**
   Decide which multipoles the source carries (a compact magnet/coil is
   dipole-dominated; a symmetric one starts at quadrupole). The DtN need only be
   accurate up to that degree; `dtn_spectrum_vs_mesh` reports the coarsest mesh
   whose `accurate_band` already covers it. This replaces "refine until the
   field stops moving" with "resolve the DtN modes the source excites."

2. **Choosing an open-BC method.** The §2 table ranks methods by which DtN modes
   they capture: a plain fixed-potential / asymptotic-Robin truncation keeps
   only `n=0` (or `n=0,1`) and needs a large air box; Kelvin and BEM DtN capture
   the whole low-mode ladder and stay accurate on a compact, coarse mesh.

3. **Trusting a coarse Kelvin model.** As long as the low-degree DtN eigenvalues
   sit on the `−(n+1)/R` ladder (they do, even on coarse meshes), the open
   boundary is faithful for the multipoles your source radiates. Do **not**
   refine reflexively — refine only if the source carries high-degree content,
   which the spectrum flags as an under-resolved band.

4. **Diagnosing an open-BC bug.** Low-mode eigenvalues far off `−(n+1)/R`
   indicate a real operator/sign/scaling bug (e.g. the `PᵀΛP` sign, the `−½M`
   exterior jump, or a wrong `R`), **not** a mesh-resolution problem — because
   mesh resolution only polishes the low-mode eigenvalues from an already-tiny
   coarse-mesh floor. This separates "the operator is wrong" from "the mesh is
   too coarse," which a pure field-refinement study cannot.

### Relation to the radia-ngsolve DtN hierarchy

| DtN type          | Approximates          | Coarse-mesh behaviour (measured)   |
|-------------------|-----------------------|------------------------------------|
| Asymptotic Robin  | exterior, `n=0` only  | exact `n=0`, wrong `n≥1`           |
| Kelvin transform  | exterior, all `n`     | order `≥ n` kills poly error, then geom floor |
| BEM DtN (`Λ_h`)   | exterior, all `n`     | low-`n` accurate, error smooth in `n` |
| SIBC / GIBC       | conductor interior    | curvature-corrected (separate)     |
| AGE               | annular air gap       | analytic per harmonic              |

The Kelvin transform and BEM DtN share the **same** continuous operator `Λ_ext`,
and both are measured here to be coarse-mesh accurate for the low modes — by
complementary mechanisms: BEM by surface-harmonic smoothness, Kelvin by
polynomial representability (order `≥ n` removes the polynomial error and drops
onto the curved-geometry floor — ~5–6 digits in 3-D = Kameari's result; the
dominant dipole inverts to a linear field, accurate at order 1).

---

## 9. The headline result, in one line

> Kameari's coarse-mesh accuracy of the Kelvin transformation is a **spectral
> property of the discrete DtN matrix**: its low-degree eigenvalues already lie
> on the `−(n+1)/R` ladder on the coarsest mesh (the dipole to 0.07 %), and the
> isolated Kelvin open-boundary error (≈ 0.1 %) sits ≈ 45× below the interior
> FEM error — readable off the operator, before any field is solved, and
> separated from the interior discretisation that a field-refinement study
> conflates.

### Computational cost — measured by the DoF increment (not solve time)

The cost of an open-boundary closure should be read off the **DoF increment** `ΔDoF`
it adds to the FE system — machine-independent and reproducible — **not** wall-clock
solve time (solver/hardware dependent). Measured that way the Kelvin closure is
**cheap**:

- The closure adds the inverted Kelvin ball = `ΔDoF` unknowns. Meshed as fine as the
  interior it is `ΔDoF ≈ N_interior` (total ≈ 2×; measured 1.72–2.03×). But the
  exterior is **volume-irrelevant**, so the ball can be a **coarse, Γ-scale** sphere.
- **Measured** (`kelvin_openbc_error_vs_exterior_mesh`, shell dipole): the coarsest
  ball — `ΔDoF = 58`, about the size of the truncation surface `Γ` — already leaves
  the closure error at `1.2e-3`, **≈ 1/45 of the interior FE error** (`5.3e-2`).
  Refining it (`ΔDoF 58 → 301`) lowers the closure error (`1.2e-3 → 7.5e-5`) but it is
  already a non-bottleneck. So the open boundary costs **~a Γ-scale coarse ball of
  DoF and never limits accuracy** — the precise sense in which Kelvin is "cheap".
- The **DtN spectrum is the measure**: the order threshold `p ≥ n_src` + the defect
  law `defect_n ~ n²(h/R)⁴` give the **minimum `ΔDoF`** for a target accuracy *a
  priori* — size the cheapest admissible exterior without a convergence sweep.

**Keep Kelvin sparse.** Condensing the exterior into the explicit DtN
(Steklov–Poincaré) operator removes the ball's DoF but turns `Γ` into a **dense
`N_Γ²` clique**: measured **nnz 10–20× the sparse extension, growing as `N^{4/3}`**
(vs the sparse extension's `N`). Condensation pays off only when `N_Γ` is small, the
exterior is reused across many solves (factor once), or the dense block is
H-matrix/FMM-compressed. *(Measured & settled 2026-06-14, hex vs tet: **a wash — neither wins for the Kelvin
sphere.** The full sphere hexes easily via `volume scheme sphere` (an earlier "impractical"
note was an error; only 1/4 & 1/8 sectors aren't covered, where tet wins by default). The
geometry floor is **two-regime**: at low curving order it tracks the curved-mesh **volume
(geometry) accuracy**; at high curving order (export order 4+) the volume error keeps
falling (hex `1.4e-5`, tet `8e-7`) while the floor **plateaus** at the FE/mesh-discretization
level (tet `~1e-5`, hex `~9e-4`) — so `floor = max(geometry-curving error, FE-discretization
error)`, the latter lowered only by mesh refinement. At **matched DoF** the two are
**comparable** (order 2, `N≈2.3k`: hex `1.6e-5` vs tet `1.9e-5`); hex's higher floor was
mostly its coarser element count, not worse geometry. So "hex lowers the floor / cuts ΔDoF"
is **not** supported. **Tet stays the practical default** (simpler, handles symmetry sectors);
high-order hex's strength stays sweepable bodies.)*

### Connection to the Cauer Ladder Network (CLN)

The `−(n+1)/R` DtN eigenvalue ladder is a **spectral** object of the same kind as
the Cauer-Ladder-Network (CLN) eigenmode decomposition: each characterises a
closure by how it transmits the modes a source excites. CLN folds the **interior**
response by circuit order `{R_n, L_n}`; the Kelvin closure folds the **exterior**
response by element order `p`. Same idea — decompose the physics into eigenmodes
and resolve only the modes that matter — applied to the interior network vs the
open boundary.

### The two scalar readouts: capacitance (n=0) and external inductance (n=1)

Capacitance `C` and external inductance `L_ext` are each **one Steklov mode of the
same exterior scalar Laplace DtN**, on different rungs of the `−(n+1)/R` ladder — so
each inherits the datasheet directly:

| quantity | exterior multipole | datasheet rung | measured open-BC defect |
|---|---|---|---|
| `C` | **monopole `n=0`** (isolated charged conductor) | `defect_0` | **0** (constant image, exact at every order; sphere `C=4π` machine-precision) |
| `L_ext` | **dipole `n=1`** (a current loop has no magnetic monopole) | `defect_1` | `1.4e-3 (p=1) → 2.4e-5 (p=2) → 7.6e-6 (p=3)` |

The bridge to inductance is the identity (hand-checked, reconfirmed numerically to
machine zero) that the DtN eigenvalue **is** the exterior-energy coefficient:

```
W_ext = ½ μ₀ · (n+1)/R · ∮_Γ φ² dS        (decaying mode r^−(n+1))
      ⇒  ½·(2/R)∮φ² = ∫_{r>R}|H|² = m²/(6πR³)   (dipole, n=1)
```

so `L_ext = 2 W_ext / I²` inherits the **dipole** defect exactly: captured at order
`p ≥ 1`, **floored at curved-Γ geometry** (Curve `k=1→3`: `4.7e-3 → 2.4e-5`),
mesh-independent, and **exterior-volume-irrelevant** (open-BC error stays ≈ 45× below
the interior FEM error at every exterior mesh). This is the precise **dual of
capacitance** — `C ↔ n=0`, `L_ext ↔ n=1`.

**Which operator is certified.** This is the open-boundary *truncation* accuracy of a
**field-energy** inductance (Kelvin / air-box), via the magnetic *potential* exterior
— scalar `Ω` (single-valued for a magnetisation source; a free-current loop needs a
cohomology cut) or the vector potential `A` (no cut; same `−(n+1)/R` gradient block).
It is **not** the `ngsolve.bem` vector single-layer energy `L = μ₀ Jᵀ(LaplaceSL)J`,
which is a *different* integral operator (the single-layer potential, no `−(n+1)/R`
ladder, order-0 current basis). So "C and L are both DtN-certified" is only **half**
true — keep the two operators distinct. (Scope: `L_ext` is the *external* energy share;
a thin loop's full self-inductance is near-field/log-dominated = interior FEM accuracy,
not a DtN question. Verified by
[`inductance_dtn.py`](../../examples/kelvin_transformation/DtN_spectrum/inductance_dtn.py).)

---

## 10. Runnable layer (source of truth)

The numbers and mechanisms above are produced by the `radia-mcp` code; query the
MCP knowledge tool for the live recipe.

- **MCP tool:** `dtn_coarse_mesh(topic="overview" | "numerics" | "api" | "applications" | "p_method" | "formulation" | "datasheet")`
  (server `mcp-server-radia-ngsolve`). Companion: `fem_bem_schur(...)`,
  `kelvin_transformation(...)` — including `topic="mesh_control"`, the
  consolidated "where to spend elements" guide (the Γ-conforming mesh
  constraint + the six measured pillars above: exterior-volume-free,
  floor = Curve order, `p ≥ n` & p-vs-h, optimal `R/a ≈ 3`, corner `hp`,
  DoF-cost `1/45`; confirmed in 2D on the `−n/R` ladder).
- **BEM side** — `radia_mcp.radia_ngsolve.bem_integral`:
  `laplace_exterior_dtn()` (assemble `Λ_h`), `exterior_dtn_spectrum()`
  (eigenvalues matched to `−(n+1)/R`), `dtn_spectrum_vs_mesh()` (per-degree
  error vs mesh).
- **Kelvin side** — `radia_mcp.radia_ngsolve.fem_bem_coupling`:
  `kelvin_dtn_eigenvalue()` (effective DtN, `dim=3` sphere `−(n+1)/R` /
  `dim=2` circle `−n/R`), `kelvin_vs_exact_open_bc_error()` (isolated open-BC
  error), `kelvin_openbc_error_vs_exterior_mesh()` (vs exterior mesh),
  `kelvin_twosphere_shell_dipole()` (real two-sphere periodic Kelvin),
  `laplace_fem_bem_schur()` (BEM DtN as an exact open BC).
- **Example:**
  [`examples/dtn_spectrum_coarse_mesh_demo.py`](../../packages/radia-mcp/examples/dtn_spectrum_coarse_mesh_demo.py)
  — Part A (BEM spectrum), Part B (Kelvin effective DtN), Part C (exterior-mesh
  sweep). Runs end-to-end.
- **Tests:**
  [`tests/test_dtn_spectrum_coarse.py`](../../packages/radia-mcp/tests/test_dtn_spectrum_coarse.py)
  (BEM spectrum + Kelvin + 2-D Kelvin),
  [`tests/test_fem_bem_coupling.py`](../../packages/radia-mcp/tests/test_fem_bem_coupling.py)
  (analytic shell, full BEM Schur, error-isolation, exterior-mesh sweep,
  two-sphere periodic Kelvin).

---

## References

1. A. Kameari, open-boundary / Kelvin-transformation magnetostatics (the
   coarse-mesh accuracy demonstration reframed here).
2. See [KELVIN_TRANSFORMATION.md](KELVIN_TRANSFORMATION.md) for the Kelvin map,
   Jacobian, material modulation, and full FEM workflow, and
   [../cln/CAUER_LADDER_NETWORK.md](../cln/CAUER_LADDER_NETWORK.md) for Kelvin
   transformation coupling within the Cauer Ladder Network.
