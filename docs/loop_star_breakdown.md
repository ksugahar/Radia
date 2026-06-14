# The Loop-Star Breakdown: high-μ magnetostatics ↔ low-frequency MoM

**One sentence:** the high-`μ_r` conditioning breakdown of the magnetization
integral equation (Radia MMM/MSC) is the *exact* magnetostatic analogue of the
classic **low-frequency breakdown** of the Electric Field Integral Equation
(EFIE / Method of Moments with RWG basis) — same cause (a "loop" subspace
pinned only by a term that vanishes in the limit), same first remedy
(**Loop-Star** decomposition).

---

## 1. Low-frequency EFIE breakdown (the classic MoM story)

The EFIE for a surface current `J` on a conductor, tested with RWG functions,
assembles to

```
Z J = E_inc ,     Z = jωμ₀ L  +  (1/(jωε₀)) S
```

- `L` — the magnetic **vector-potential** operator; acts on `J` directly
  ("inductive").
- `S` — the **scalar-potential** operator; acts on the charge `ρ = (1/jω)∇·J`
  ("capacitive").

Split the current into **loop** (solenoidal, `∇·J = 0`) and **star/tree**
(`∇·J ≠ 0`, charge-carrying):

- `S` annihilates the loops: `S·J_loop = 0` (no charge → no scalar potential).
- A loop current is therefore seen **only** by the inductive term:
  `Z·J_loop = jωμ₀ L·J_loop  ~ ω`.
- A star current is dominated by the charge term: `Z·J_star ~ 1/ω`.

As `ω → 0` the loop block scales like `ω` and the star block like `1/ω`; their
ratio is `~ω²` and `cond(Z) ~ 1/ω²` → the matrix is numerically singular. The
loop currents are constrained **only** by the vanishing `jωL` term, so they are
nearly free and the solve (iterative or direct) breaks down.

**Fix — Loop-Star / Loop-Tree:** change basis to `{loops, stars}`, which
block-separates the `ω`-scaling, and rescale each block to an `O(1)` dynamic
range. (Wilton–Glisson RWG 1982; Vecchi, *IEEE TAP* 1999; Zhao–Chew 2000;
Andriulli et al. Calderón, *IEEE TAP* 2008.)

---

## 2. High-μ magnetostatic breakdown (Radia MMM/MSC)

The magnetization integral equation for a linear soft magnet
(susceptibility `χ = μ_r − 1`) is

```
( (1/χ) I  −  N ) M = H_ext ,     N = the demagnetising operator (M ↦ H_demag)
```

where `N` acts **through the magnetic charges** of `M`:
`ρ = −div M` (volume), `σ = M·n` (surface). Split `M` into **loop**
(charge-free: `div M = 0` **and** `M·n = 0`) and **star** (charge-carrying):

- `N` annihilates the loops: `N·M_loop = 0` (no charge → no field). The loops
  are exactly `ker N`.
- A loop magnetization is therefore seen **only** by the regularising term:
  `A·M_loop = (1/χ) M_loop`.
- A star magnetization is dominated by the demag term:
  `A·M_star ≈ −N·M_star ~ O(1)`.

As `μ_r → ∞`, `1/χ → 0`: the loop block has eigenvalue `1/χ → 0` while the star
block is `O(1)`; `cond(A) ~ χ ~ μ_r` → the matrix becomes singular. The loop
magnetizations are constrained **only** by the vanishing `(1/χ)I` term — exactly
as the EFIE loops are constrained only by the vanishing `jωL`.

---

## 3. The correspondence (the dictionary)

| | EFIE, `ω → 0` | MMM/MSC, `μ_r → ∞` |
|---|---|---|
| unknown | surface current `J` | magnetization `M` |
| **loop** modes | solenoidal `J`, `∇·J = 0` | charge-free `M`, `div M = 0` & `M·n = 0` |
| loops produce | no charge → no scalar potential | no charge → no field (`M_loop ∈ ker N`) |
| operator | `Z = jωμ L + (1/jωε) S` | `A = (1/χ) I − N` |
| term that **pins the loops** | `jωμ L`  (`~ω`) | `(1/χ) I`  (`~1/χ`) |
| dominant term | `(1/jωε) S`  (`~1/ω`) | `−N`  (`~O(1)`) |
| conditioning | `cond ~ 1/ω²` | `cond ~ μ_r` |
| null operator | charge op `S` (`S J_loop = 0`) | demag op `N` (`N M_loop = 0`) |
| first remedy | **Loop-Star / Loop-Tree** | **Loop-Star** (de-Rham loops) |

The dictionary is **(solenoidal current, frequency ω) ↔ (charge-free
magnetization, 1/χ)**. The EFIE charge operator `S` plays the role of the
magnetostatic demag operator `N`; both have the loops as their exact null space,
and in both the loop-pinning term vanishes in the limit.

---

## 4. Why the loops must be EXACTLY field-null — yano-type vs HDiv-type element

Loop-Star only works if the **discrete** loops are exactly in the null space
(`S J_loop = 0`; `N M_loop = 0`). For RWG on a flat triangulation this is
automatic — the RWG loop functions are exactly solenoidal. For the magnetostatic
**volume** problem on **distorted hexes** it is *not* automatic:

> the naive constant-magnetization / `±1`-cycle "loops" are field-null only on
> **affine** elements; on a **non-affine (distorted)** hex they carry residual
> charge → `N M_loop ≠ 0` → the Loop-Star separation leaks and the field is
> wrong (Radia's "C-type" magnet: ~28 % `B_z` error at `μ_r = 1e5`).

Two ways to restore exactly field-null loops on distorted elements:

- **yano-type** — engineer the element shape functions so the loop component
  produces no field (Yano's distortion-handling MMM/MSC elements).
- **HDiv-type** — use the NGSolve **H(div) (FEEC)** basis, where the loops
  `= curl(interior H(curl))` are charge-free (`div = 0` **and** `M·n = 0`) **by
  construction**: the contravariant **Piola** map preserves *both* the
  divergence and the normal trace, so the loops stay field-null under *any*
  distortion. This is the discrete **de Rham** analogue of RWG's exactly
  solenoidal loops. (Verified to machine precision on distorted hexes —
  `examples/vim/`.)

Either way, the H-matrix acceleration is unchanged: the demag operator is the
same dense `1/r` (Laplace) integral operator, so it is compressed by **HACApK**
(ACA⁺) exactly as the present MSC is.

---

## 5. The honest limit — Loop-Star fixes the loops, not the whole conditioning

Loop-Star removes the loop-induced singularity (the `1/χ → 0` mode), which fixes
field **correctness** — the loops no longer pollute the solution (the
magnetostatic *Problem A*). It does **not**, by itself, bound the conditioning of
the remaining **star** block: the demag operator `N` carries a near-continuum of
small singular values (weak-demag, charge-carrying modes) *above* the exact loop
kernel, so the star-block conditioning can still grow with `μ_r` (the
magnetostatic *Problem B*; empirically deflating only `ker N` barely moves
`cond`). This is the magnetostatic counterpart of the EFIE **dense-mesh
(h-refinement) breakdown** that survives Loop-Star and needs a second-kind /
**Calderón** preconditioner. In Radia the residual high-`μ_r` conditioning is
carried by the **H-ILU** preconditioner (the HACApK H-matrix `A_SS` factor),
validated for `μ_r ≤ 1e4`; a scalable spectral coarse space (GenEO-type) is open.

### Does it have to be H-ILU? — the loop-quality hypothesis

How strong a preconditioner the **star** block needs depends on **how exactly
the loops were removed**:

- With **inexact** loops (the yano `±1` cycles on distorted hexes), the
  un-removed *near-loop* modes — charge-carrying magnetizations with a tiny,
  nonzero demag eigenvalue — stay **inside** the star block as a weak-demag
  near-null cluster. That cluster is what forces a strong preconditioner
  (**H-ILU**). Empirically, deflating only the exact `ker N` of that matrix
  barely moved `cond` — consistent with the near-null mass being *near*-loops,
  not the exact kernel.
- With **exactly** field-null loops (de-Rham / HDiv, `N M_loop = 0` to machine
  precision), the Loop-Star separation removes that cluster **as part of the
  loops**, so the star block may be well-conditioned enough for a **plain
  diagonal (Jacobi)** preconditioner — the same way EFIE Loop-Star + diagonal
  rescaling tames the low-frequency breakdown without a Calderón operator at
  moderate mesh density.

**This is a concrete, testable hypothesis**, not yet established: measure the
HDiv-type star-block iteration count under **Jacobi** vs `μ_r` and mesh size. If
Jacobi keeps it bounded for `μ_r ≤ 1e4` at moderate meshes, the HDiv-type needs
**no H-matrix factor** — which would make it not only correctness-equivalent but
*cheaper* than the yano-type + H-ILU, i.e. a genuine **compute-time superset**
(the switch criterion). If a weak-demag continuum survives even with exact loops
(thin/elongated geometry, very fine meshes), a stronger rung is needed.

### The preconditioner ladder — where Jacobi, AMS, and Calderón sit

**Loop-Star is the lowest-order auxiliary-space decomposition.** Its systematic
generalization is **AMS / ADS** (Hiptmair–Xu 2007 for H(curl); Kolev–Vassilevski
for H(div)) = "loop-star + algebraic multigrid on the auxiliary gradient/curl
spaces". The FEEC/HDiv basis here *is* that decomposition — the loops are
`curl(H(curl))`, the curl-auxiliary space.

Two features make the VIM **easier** than the sparse-FEM curl-curl problem AMS
was built for:

1. The loops are **field-null** (`N M_loop = 0`). In sparse FEM the gradient
   kernel is *coupled* and AMS must actively correct it; here the loop /
   curl-auxiliary block is **decoupled** — the demag operator ignores it,
   `(1/χ)I` handles it trivially. **Only the star (charge) block needs a real
   preconditioner.**
2. That star block is a **dense** `1/r` integral operator stored as an
   **H-matrix** (HACApK).

AMS/ADS are formulated for **sparse** FEM operators. An *"H-matrix AMS/ADS"* —
the auxiliary-space construction applied to a dense H-matrix integral operator —
**does not appear to be in the literature**: the dense-operator world
preconditions with **Calderón / operator preconditioning** instead (Hiptmair
2006; Andriulli et al. 2008), and even those are for *boundary* integral
equations — *volume*-IE preconditioning is less settled. So an H-matrix-AMS for
this volume VIM would be genuinely new.

Ladder for the HDiv-type star block, cheapest first — **repo-first: use the
simplest rung that works; the novelty is not the goal, a conditioned solver is:**

1. **Jacobi** (diagonal) — try first; may suffice once the exact de-Rham loops
   are removed (the weak-demag cluster leaves with them).
2. **Auxiliary-space / operator preconditioning** — the principled next step; the
   lab's HYPRE-free **Compact AMS** (`radia.sparsesolv_ngsolve`) already provides
   the sparse-auxiliary machinery to build on, and a dense (H-matrix) variant
   would be the novel "H-matrix AMS".
3. **H-ILU** (HACApK H-matrix factor) — the current fallback, validated
   `μ_r ≤ 1e4`.

### H-matrix preconditioner roadmap

The HACApK H-matrix (`O(N log N)` operator apply) is the shared substrate; the
preconditioners on it form a suite with graceful escalation (repo-first: use the
lightest that works — Jacobi where the evidence below allows):

- **H-LU / H-ILU** — *keep* (current factor preconditioner + scalable fallback;
  validated `μ_r ≤ 1e4`). Not deleted even if Jacobi handles the main path.
- **H-QR** (A. Ida's proposal) — *future*. An H-matrix QR factorization: more
  numerically stable than H-LU (no pivot growth on ill-conditioned high-`μ`
  systems), enables least-squares / overdetermined solves, and — being
  **rank-revealing** — can expose the operator's null space (the **loops**)
  directly from the factorization, i.e. a factorization-level loop/star split.
- **H-AMS** — *future*. The auxiliary-space (Hiptmair–Xu) preconditioner applied
  to the dense H-matrix operator, with the de-Rham loop/star as its lowest-order
  auxiliary decomposition. Appears unpublished for dense integral operators (see
  "Does it have to be H-ILU?" above) — a genuine novelty.

### Evidence so far (favours Jacobi for `μ_r ≤ 1e4`)

`examples/vim/demag_spectrum_jacobi.py` measures the demag operator's spectrum
with Radia's **exact** `ObjRecMag` field on a compact body (3×3×3 cube grid,
constant-M = an *all-star*, loop-free operator):

| grid | smallest `\|μ\|` | `cond(A)` @ `μ_r=1e4` | Jacobi-GMRES iters @ `μ_r=1e4` |
|---|---|---|---|
| regular | 0.051 | 17 | 10 |
| distorted (30 % jitter) | 0.0076 | 1.1e2 | 84 |

In both, `cond(A)` is **`μ_r`-independent** (no `~μ_r` blowup) and Jacobi-GMRES is
**bounded**, because the smallest charge `|μ| (≈0.01–0.05) ≫ 1/χ (=1e-4` at
`μ_r=1e4)` — a large spectral gap. So **the charge-carrying (star) demag operator
of a moderate compact body is Jacobi-solvable at `μ_r ≤ 1e4`**; the `~μ_r`
breakdown needs `μ_min < 1/χ`, which on the C-type came from the **yano `±1` loops
not being exactly field-null** (spurious `μ≈1/χ` modes) rather than a physical
weak-demag continuum. This **supports the hypothesis**: with the HDiv-type's
*exactly* field-null loops removed, the star block should be Jacobi-solvable for
`μ_r ≤ 1e4` (no H-matrix factor → also the compute-time superset).

A stress sweep (`demag_jacobi_stress.py`: cube / rods up to 1×1×16 / 4×4×1 slab /
C-yoke, regular **and** 30 %-distorted) confirms the regime: in **every** case
`|μ|_min ≈ 0.007–0.12 ≫ 1/χ`, `cond(A)` stays **bounded (≤ ~270)** with **no `~μ_r`
blowup**, and Jacobi-GMRES converges in **≤ ~200 iters** at `μ_r = 1e4`. A subtle
point: the conditioning is *worst* where `1/χ ≈ μ_min` (gap ≈ 1) — at *moderate*
`μ_r` (~1e2–1e3) — and *improves* as `μ_r → ∞` (the gap opens). So for the
charge-carrying star block of a moderate body, Jacobi is fine across `μ_r`.

Caveats: this is the constant-M proxy on moderate bodies; the **stress test that
remains** is the HDiv-type star block on the genuinely distorted C-type with a
*thin* magnet gap, and an ultra-weak-demag mode (extreme thin features, very fine
meshes — where `μ_min` drops below `1/χ`) would still need the H-ILU fallback.

---

## Takeaways

1. The **high-μ breakdown IS the low-frequency breakdown**, with
   `(charge-free M, 1/χ)` in the role of `(solenoidal J, ω)`.
2. **Loop-Star** restores field correctness — *provided the discrete loops are
   exactly field-null* (de-Rham/HDiv basis, or Yano-engineered elements). On
   distorted hexes the naive `±1` loops are **not** exact, which is the whole
   reason the distortion-robust element exists.
3. The **residual conditioning** is a separate, preconditioner-level problem
   (H-ILU / HACApK), analogous to EFIE's post-Loop-Star dense-mesh breakdown —
   not something the loop basis alone can fix.

---

*See also:* `examples/vim/` (the HDiv-type loop/star split, loops field-null
on distorted hexes), `src/core/rad_hacapk.cpp` (`BuildLoopBasis` /
`SolveLoopStar`, the de-Rham-exact loop basis + H-ILU `A_SS` solve),
`docs/HMATRIX_EVALUATION.md` (HACApK ACA⁺).
