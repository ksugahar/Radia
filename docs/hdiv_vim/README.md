# HDiv-type VIM — the FEEC Volume Integral Method for MMM/MSC demagnetization

**One sentence:** the magnetization "loop" modes that break the constant-M MMM/MSC integral
equation on distorted meshes (see [`../loop_star_breakdown.md`](../loop_star_breakdown.md)) are
**field-null by construction** when the magnetization lives in NGSolve's H(div) (RT) finite-element
space — so the **HDiv-type VIM** complements the canonical **multipole-moment MMM** MSC backend with a
de-Rham-exact operator for **tetrahedral RT1 soft-iron demag**, including curved P2 tetrahedral
geometry.  Non-tet bodies, symmetry-reduced image models, and permanent-magnet mixing are intentionally
routed to the collocation MMMM backend.

The tradeoff is deliberate: HDiv-VIM gives a symmetric Galerkin matrix and high-order extensibility,
but the charge-Coulomb Gram integrals dominate matrix construction.  Multipole-moment MMM gives up
that full Galerkin symmetry and instead uses Mathematica-derived moment rows, so the 3-DOF and
surface-charge MMM systems keep cheap local element functionals and scalable H-matrix matvecs.

This is the canonical technical reference. The runnable examples + their numbers live in
[`examples/vim/README.md`](../../examples/vim/README.md); the decision/narrative record is
the radia-mcp `hdiv_vim` knowledge (MCP tool `mcp__radia-ngsolve__hdiv_vim`).

---

## 1. The operator

The demagnetizing operator is assembled as a symmetric Galerkin **charge–Coulomb** form:

```
N = Bᵀ G B
```

- **B** — the charge map of the magnetization `M ∈ HDiv(order p)`:
  `B M = ( ρ = −div M  on L2(p) ,  σ = M·n  on SurfaceL2(p) )`  (volume charge ρ, surface charge σ).
- **G** — the Coulomb Gram between charges, `G_ij = ∫∫ q_i q_j / (4π r)`.
- **M_mass** — the HDiv mass; the physical demag factors are the eigenvalues of `M_mass⁻¹ N`
  (basis-invariant), and `D_axis = (mᵀ N m)/(mᵀ M_mass m)` for the uniform mode m.

**Why the loops are field-null by construction.** A magnetization is field-null (a "loop", in
`ker N`) **iff it is charge-free** — `ρ = 0` *and* `σ = 0`. The loops are exactly `ker B`, and since
`N = Bᵀ G B`, `N·loop = Bᵀ G (B·loop) = 0` for **any** Gram G, on **any** mesh (affine, distorted, or
cut for symmetry). The contravariant Piola map preserves both `div` and the normal trace (the de Rham
commuting diagram), so distortion cannot break it. This is the strong (everywhere) field-null property,
versus the constant-M basis's fragile (collocation-only) field-null that breaks under distortion and
forces the cohomology-aware `installCycle` loop-star construction in multipole-moment MMM MSC.

## 2. The Gram — operator choices, pick by problem

| Gram | Used for | Accuracy |
|---|---|---|
| centroid-monopole + sub-point self | quick probes, near-field correction | crude (~2–3% on the demag factor) |
| **Wilton analytic surface** (`wilton_surface=True`) | **uniform-M linear demag** (div M = 0 → surface charge only) | sphere/cube → 1/3 to `<0.15%` |
| **full analytic volume** (`analytic_gram=True`, `phi_tet`) | **NON-uniform / nonlinear** (div M ≠ 0 → volume charge) | required — see §4 |
| **Retired Gauss-point H-matrix** | not a live public backend | it was an RT0 build-speed experiment and is now fail-loud; RT1 uses the analytic charge Gram |
| **Curved P2 tetrahedral charge Gram** | curved RT1 tet demag | matched curved cell/face geometry; no straight-Gram drift on curved tetrahedra |

The surface single-layer is the key architectural unlock: the uniform-M surface demag Gram **is** the
Laplace single-layer of σ = M·n.  The current shipped route is the analytic charge Gram used by the
RT1 tetrahedral solver; earlier `ngsolve.bem` / Gauss-point experiments remain research history rather
than public backends.

## 3. Validation matrix (golden-locked) — with reference honesty

What each number is measured against, precisely:

| Result | Validated against | Status |
|---|---|---|
| Loops field-null on distorted hex | exact | `4e-16` (machine zero) |
| Linear demag, sphere/cube → 1/3 | **ANALYTIC** | `<0.15%` (Wilton surface Gram) |
| Spheroid tensor (prolate+oblate, polar+transverse) | **ANALYTIC** Osborn | curved+o2 exact; sum rule `N_x+N_y+N_z=1` to `~1e-6` |
| General triaxial ellipsoid (3 distinct factors) | **ANALYTIC** Osborn integral | all exact (`hdiv_demag_bem_singlelayer.py`) |
| Nonlinear, sphere deep-saturation | **ANALYTIC** fixed point | `<0.05%` (damped Newton) |
| Nonlinear cross-check + real BH table | Radia / ANALYTIC | sphere `<0.05%`; BH table `<0.2%` |
| Nonlinear cube / C-yoke (non-convex) | **Radia** (both flat → valid) | volume-avg M_z `<1%` at every mesh, 6 iters |
| Curved + high-order demag | **ANALYTIC** | curved+o2 exact; flat floored |
| Curved field vs shipped Radia (flat) | **ANALYTIC** dipole | accuracy-per-DOF **~10–30×** |
| Curved × nonlinear field | **ANALYTIC** dipole | flat `~9%` → curved `<0.4%` (~23×) |
| Symmetry models 1/2, 1/4, 1/8 | retired from HDiv-VIM | use collocation MMMM for reduced image/symmetry models |

**Reference rules:** sphere / spheroid / ellipsoid / dipole are validated vs **analytic truth** (real
errors). Cube / C-yoke have **no analytic solution** — validated vs **shipped Radia** (a trusted
solver, *valid here because both are flat*); stated as "agreement", `<1%` volume-avg. **Radia cannot
referee curved geometry** (its `ObjHex/Tet` facet the body), so curved results are validated vs
analytic only.

## 4. Nonlinear — the Gram requirement (a trap, now fail-loud)

A NON-uniform-M nonlinear solve (`div M ≠ 0`: cube, C-yoke, any non-ellipsoid) **requires
`analytic_gram=True`** (the full volume Gram). The surface-only `wilton_surface` Gram leaves the volume
(cell) blocks crude → wrong per-element fields → damped Newton does **not** converge (stalls at maxit,
M_avg drifts). Only uniform-M nonlinear (sphere/spheroid, div M = 0) converges with `wilton_surface`.

`solve_nonlinear_newton` now **raises** on non-convergence (`require_convergence=True` default), with a
message pointing at `analytic_gram` — no silent wrong result. With the right Gram the C-yoke converges
in **6 Newton iters** and matches Radia to `<1%`.

## 5. Symmetry models (1/2, 1/4, 1/8)

Image / IMA symmetry is retired from HDiv-VIM.  HDiv-VIM is now the tet/RT1 soft-iron path; reduced
symmetry models, hex/wedge bodies, and permanent-magnet mixing are the responsibility of collocation
MMMM.  This keeps the public HDiv contract small and fail-loud instead of carrying experimental image
rules in the Galerkin path.

## 6. HDiv-type vs multipole-moment MMM

| Capability | multipole-moment MMM MSC | HDiv-type VIM |
|---|---|---|
| Linear demag (sphere/spheroid/triaxial) | ✓ | ✓ exact vs analytic |
| Nonlinear (cube / C-yoke) | ✓ | ✓ `<1%` vs Radia, 6 iters |
| Distorted-mesh μr-independence | ✓ (hand-crafted) | ✓ **by construction** (`4e-16`) |
| Symmetry 1/4, 1/8 | ✓ (loop-star / image handling) | retired; use collocation MMMM |
| **Curved tetrahedral geometry** | ✗ (flat elements) | ✓ RT1 + curved P2 tet geometry |
| Hand-crafted elements | required | **not needed for tet RT1** (de-Rham-exact) |

## 7. Honest status & open work

> **The concrete milestone roadmap for production HDiv-VIM is [PRODUCTIONIZATION.md](PRODUCTIONIZATION.md)**
> — current C++/Python inventory, the parity-gate definition-of-done (incl. the unmeasured speed gap),
> and milestones M0 (parity + speed measurement) -> M5 (production seal).

The HDiv-type VIM is a **validated research prototype** (Python + NGSolve) with a quantified
accuracy-per-DOF win over the shipped flat solver on curved problems, and parity on the flat cases. The
remaining lift to make it a **production backend alongside multipole-moment MMM MSC**:

1. **C++ productionization** — keep the live analytic charge Gram + RT1 Newton path behind a Radia API
   and measure it honestly against collocation MMMM.  The retired Gauss-point and H-LU paths are not
   part of the public HDiv-VIM contract.
2. **Curved nonlinear volume charge** — `ngsolve.bem` is boundary-only, so non-uniform nonlinear on
   curved cells still needs the Newtonian volume potential (`phi_tet`) on curved geometry.

## 8. Code map

| Concern | Example | Golden test |
|---|---|---|
| Loops field-null / Hodge split | `ngsolve_loopfree_verify.py`, `hdiv_loop_star_split.py` | `test_hdiv_vim_solve.py` |
| Operator + Wilton / volume Gram | `hdiv_demag_tet.py` | `test_hdiv_vim_{wilton,volume}_gram.py` |
| Nonlinear Newton (+ BH table, scalable) | `hdiv_demag_tet_nonlinear.py` | `test_hdiv_vim_{tet_newton,newton_vs_radia,newton_table,newton_scalable}.py` |
| Curved + high-order demag (single-layer) | `hdiv_demag_bem_singlelayer.py` | `test_hdiv_vim_bem_demag.py` |
| Curved geometry / field win | `hdiv_demag_curved.py`, `hdiv_curved_nonlinear_field.py` | `test_hdiv_vim_curved{,_nonlinear,_nonlinear_field}.py` |
| Head-to-head vs shipped Radia | `compare_curved_vs_radia_field.py` | `test_curved_vs_radia_field.py` |
| C-yoke nonlinear (non-convex) | `hdiv_cyoke_nonlinear.py` | `test_hdiv_vim_cyoke_nonlinear.py` |
| RT1 / tet-only public contract | `hdiv_demag_solve(..., order=1)` | `test_hdiv_vim_rt1_contract.py` |

All under `examples/vim/` and `validation_test/feec/` (the FEEC suite is split out of lightweight CI).

## 9. Research plan — the eddy-current VIM (future directions)

HDiv-VIM above solves the **magnetostatic** demag operator: magnetization in H(div)
splits into charge-carrying modes (which drive demag) and field-null **loop** modes
(`ker B` — charge-free, divergence-free, zero normal trace). The **eddy-current**
problem is the natural next VIM, and its unknown is exactly the *solenoidal* part:
eddy current is divergence-free (`∇·J = 0`), i.e. it lives in the loop space.

**Motivating negative result.** The eddy-current VIM route taken so far — a
Newton-kernel volume-Galerkin **Nagamine–Foster–Born series** (the `radia_vim`
prototype, deleted 2026-06-14) — is impractical. The obstacle is the **Foster-series
summation itself**: slow convergence at the wall band and the high-N Hankel/QD
breakdown in float64 (see memory `foster-convergence-central-obstacle`,
`cln-high-stage-degrades-below-foster`). An efficient summation would help, but none
is in hand, so the series route is set aside. (An extended-precision "DD" port of the
hex VIM lives under `examples/maglev/research_cln/ngsolve_validation/dd_*` — a
separate line attacking the float64 breakdown directly; it is NOT part of the deleted
engine.)

Two **matrix-free, non-series** routes to the same solenoidal eddy-current VIM, both
built the way the production HDiv-VIM is built (FEEC, de-Rham-exact, analytic field
operator — the `N = Bᵀ G B` machinery), are **unverified research directions**:

- **(A) HCurl-VIM** — eddy current as `J = curl T` with the current vector potential
  `T` in NGSolve's H(curl) (Nédélec) space (curl-conforming, natural for curl-curl /
  the A-formulation). *May be revived* as a sibling to HDiv-VIM.
- **(B) loop-basis-only VIM** — expand the unknown **directly in the loop subspace
  `ker B`** (the divergence-free, field-null modes that HDiv-VIM already constructs
  automatically on any mesh). Since the demag VIM already builds and validates `ker B`
  (loops field-null to `4e-16`, §1/§5), restricting a VIM to that basis is a small
  step from the accumulated HDiv-VIM work — it just *uses* what HDiv-VIM discards.

(A) and (B) target the **same** space by the de Rham complex: `ker(div)` in H(div)
(the loops) **=** `range(curl)` (curls of H(curl) potentials). They are two
representations of the solenoidal eddy-current VIM; (B) reuses the existing `ker B`
construction, (A) uses the curl-conforming potential.

**Prerequisite / sequencing.** Do HDiv-VIM productionization first
([PRODUCTIONIZATION.md](PRODUCTIONIZATION.md) M0–M5); both routes reuse its operator /
Gram / Newton machinery. All **unverified** — directions, not results.
