# HDiv-type VIM — the FEEC Volume Integral Method for MMM/MSC demagnetization

**One sentence:** the magnetization "loop" modes that break the constant-M MMM/MSC integral
equation on distorted meshes (see [`../loop_star_breakdown.md`](../loop_star_breakdown.md)) are
**field-null by construction** when the magnetization lives in NGSolve's H(div) (RT) finite-element
space — so the **HDiv-type VIM** is the **primary accurate soft-iron demag route** (decision
2026-06-30; updated 2026-07-05 toward Radia HDiv-only after the MMMM migration to
`ELF_MAGIC@研究室版`), with a de-Rham-exact operator for **RT1 demag
on tetrahedra (flat + curved P2), pure-hex meshes (flat + curved Q2, Piola-exact charges,
H-matrix build), and flat pure-wedge meshes**, plus a **2D planar tri/quad layer** (log
kernel, motor cross-sections). Permanent-magnet regions mix directly into `Solve` (`pm_M=`);
flat pure-TET / pure-HEX / pure-WEDGE symmetry-reduced image models are supported by the
charge-Gram IMA path. Mixed/pyramid bodies are temporary migration gaps rather than a reason to keep
MMMM as the Radia production backend.

The tradeoff is deliberate, and the main product motive is **Reduced FEM coupling**: HDiv-VIM lives in
NGSolve's finite-element world (`Mesh`, `GridFunction`, `CoefficientFunction`, `BilinearForm`,
`TaskManager`), so VIM iron can be coupled weakly to conductor / motor FEM without translating through a
Radia-only collocation object layer.  It also gives a symmetric Galerkin matrix, high-order / curved
geometry, and 2D planar support; the charge-Coulomb Gram integrals dominate matrix construction, but the
HACApK charge-Gram makes the build scalable at engineering size.  Multipole-moment MMMM remains useful in
`ELF_MAGIC@研究室版` and as a transitional Radia cross-check, but Radia's production direction is
HDiv-only.  Very small DoF timing is treated as "both fast", not as a backend selection criterion.

This is the canonical technical reference. The runnable legacy corpus that the
FEEC goldens still import now lives under
[`validation_test/feec/vim_legacy`](../../validation_test/feec/vim_legacy);
the retired source-tree prototype inventory is recorded in
[`vim_examples_retirement.ipynb`](vim_examples_retirement.ipynb) with
[`vim_examples_retirement_results.json`](vim_examples_retirement_results.json).
The decision/narrative record is the radia-mcp `hdiv_vim` knowledge (MCP tool
`mcp__radia-ngsolve__hdiv_vim`).

**Executed showcase notebooks** (code + embedded results + `_result.json` sidecars, this directory):
[`hdiv_curved_showcase.ipynb`](hdiv_curved_showcase.ipynb) (curved P2 tet win),
[`polynomial_charge_field.ipynb`](polynomial_charge_field.ipynb) (RT1 polynomial charges),
[`hex_rt1_and_2d_showcase.ipynb`](hex_rt1_and_2d_showcase.ipynb) (2026-07-03: hex RT1 + H-matrix
build timing + 2D planar closed-form gates + the production `Solve` one-call), and
[`hex_vs_mmmm_crossvalidation.ipynb`](hex_vs_mmmm_crossvalidation.ipynb) (2026-07-04: HDiv-VIM hex
RT1 vs collocation-MMMM hex on the same cube and C-yoke meshes).
The HDiv-only transition benchmark plan is
[`HDiv_vs_MMMM_benchmark_plan.md`](HDiv_vs_MMMM_benchmark_plan.md).

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
| Pure-hex RT1 (affine / distorted / real-warp hex) | physical bound + cube 1/3 | eig ∈ [0,1] incl. the worst real Cubit cylinder hex (pre-Piola leaked to 1.0105); `hex_rt1_and_2d_showcase.ipynb` |
| Curved hex cylinder (Q2 geometry) | **ANALYTIC** moment (volume capture) | moment err 0.13%; max eig 1.0078 — the one documented bound violation, open item |
| 2D disk / ellipse / Clausius–Mossotti | **ANALYTIC** closed forms | disk demag 1/2 exact (0.50000); ellipse 2:1 → 0.33438/0.66562; CM `2–3e-4` |
| Symmetry models 1/2, 1/4, 1/8 | explicit full model + parity | flat pure-TET/HEX/WEDGE IMA; `M_avg` is full-domain, `M_avg_reduced` is diagnostic |

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

Image / IMA symmetry is live for **flat pure-TET, pure-HEX, and pure-WEDGE RT1** HDiv-VIM.
The reduced mesh solve folds mirror charges into the charge-Gram H-matrix; when called through
`rad.Solve(..., image=...)`, the Radia bridge also materializes explicit mirror polyhedra after
write-back.  Therefore `rad.Fld(iron, ...)` on the solved `MeshSoftIron` container evaluates the
**materialized full field of that reduced solution**, suitable for probes and visualization.  This is a
roundoff-level contract against `res["field_object"]`, not a license to accept percent-level differences
against an unconstrained explicit full solve.  The latter still exposes a hex RT1 write-back gap: small
transverse components in the full solve follow the hex ChargeGram quadrature details, so the 10-eps
explicit-full `rad.Fld` target is tracked as an xfail validation item until the reflection symmetry of the
hex RT1 operator is fixed.  The result dict reports:

- `M_avg`: physical full-domain average magnetization; odd components across an image plane cancel.
- `M_avg_reduced`: raw reduced-domain average, kept for diagnostics and debugging only.
- `image_field_handles`: explicit mirror element handles added for `rad.Fld` when `image=` is used.

Curved reduced models, mixed-element bodies, and pyramid-containing meshes fail loud or route to the
collocation MMMM bridge; do not silently drop image symmetry in those cases.

## 6. HDiv-type vs multipole-moment MMM

| Capability | multipole-moment MMM MSC | HDiv-type VIM |
|---|---|---|
| Linear demag (sphere/spheroid/triaxial) | ✓ | ✓ exact vs analytic |
| Nonlinear (cube / C-yoke) | ✓ | ✓ `<1%` vs Radia, 6 iters |
| Distorted-mesh μr-independence | ✓ (hand-crafted) | ✓ **by construction** (`4e-16`) |
| Symmetry 1/4, 1/8 | ✓ (image handling) | ✓ flat pure-TET/HEX/WEDGE; `rad.Fld` materializes images |
| **Curved tetrahedral geometry** | ✗ (flat elements) | ✓ RT1 + curved P2 tet geometry |
| **Pure-hex meshes** | ✓ (legacy / ELF / transition cross-check) | ✓ RT1, Piola-exact charges, flat + curved Q2, H-matrix build |
| **2D planar tri/quad (motor cross-sections)** | ✗ | ✓ log-kernel Gram, closed-form gated |
| Wedge / mixed-element bodies | ✓ | pure-wedge RT1 live; mixed/pyramid still bridge/open |
| Hand-crafted elements | required | **not needed** (de-Rham-exact) |

## 7. Honest status & open work

> **The milestone record for production HDiv-VIM is [PRODUCTIONIZATION.md](PRODUCTIONIZATION.md)**
> — C++/Python inventory, the parity-gate definition-of-done, milestones M0–M5, and the dated
> progress log (latest: 2026-07-03 — hex RT1 + 2D layers shipped; build speed measured).

**Status 2026-07-03: the HDiv-type VIM is a production C++ backend**, not a Python prototype — the
analytic charge Gram is an all-C++ HACApK H-matrix (symmetric leaf fill + static-site radial inner
quadrature, ~10× vs the bring-up build on the cylinder benchmarks), the linear solve is the
symmetric mass-Riesz CG, the nonlinear tet path is the all-C++ energy-Newton (deep-saturation
robust), and the 2026-06-30 / 2026-07-05 role decision makes it the **primary accurate Radia
soft-iron route**.  Collocation MMMM is now a transition cross-check in Radia and the active
continuation belongs in `ELF_MAGIC@研究室版`.  Shipped and golden-locked: tet RT1 flat + curved P2
(linear/nonlinear/PM-mixed/per-region/holed bodies), pure-hex RT1 flat + curved Q2 (Piola-exact
charges; the ~20k-charge use-after-free is fixed, commit `20e6e9e2`), and the 2D planar tri/quad
layer (commit `a9999dd7`).  Evidence: `validation_test/feec/` (40+ HDiv goldens) + the executed
showcase notebooks above.

**Open work** (honest list, tracked in memory `hdiv-tet-hex-coupling-pyramid-gated` /
`hdiv-vim-tri-quad-motor`):

1. **Curved-hex eigenvalue 1.0078** — the curved cylinder's demag spectrum still slightly exceeds
   the [0,1] bound (halved from 1.0166 by the site anchors); remaining self/touching curved
   quadrature refinement.
2. **Hex public API — DONE (2026-07-04).**  `Solve(hexmesh, mu_r=/bh_table=)` and
   `rad.Solve(demag_backend='hdiv')` now solve a pure-hex mesh **LINEAR + NONLINEAR** (the C++
   energy-Newton was already Gram-agnostic — it takes the hex `(H, B, M_mass)` unchanged), matching
   collocation MMMM to ~1% (golden `test_hdiv_vim_hex_public_solve.py`).  The `auto` default routes
   mesh-backed pure TET/HEX/WEDGE iron to HDiv-VIM; explicit
   `demag_backend='collocation_mmmm'` is a transitional cross-check while the MMMM migration/deletion is
   audited.  The legacy `solve_nonlinear_newton_scalable`
   (`tet.build_demag` head-to-head path) stays tet-only. **2D nonlinear** iron is the 2D planar layer's
   own track.
3. **`rad.Fld` / application contract hardening (2026-07-05):** `rad.Solve` now has validation gates
   for HDiv write-back into Radia field objects, IMA full-domain `rad.Fld`, HACApK charge-Gram build
   stats, and a minimal planar motor saliency angle sweep.  Important correction: IMA `rad.Fld` is exact
   for the materialized reduced solution, while unconstrained explicit-full hex RT1 `rad.Fld` is not yet
   a machine-precision oracle because the hex ChargeGram still has a small reflection-symmetry defect.
   The remaining application work is to fix that explicit-full field parity, then promote the validation
   gates into the radia-motor operating API and larger mdx benchmarks.
4. **VIM ↔ reduced-FEM weak coupling — verified end-to-end + promoted to docs (2026-07-04).**  The 2D
   maglev-pattern weak coupling (open-boundary VIM iron + reduced complex A_z FEM on the conductor) is
   shipped as the executed showcase [`../electric_machine/planar_vim_motor.ipynb`](../electric_machine/planar_vim_motor.ipynb)
   (+ `em_reference_audit.ipynb`): the nonlinear-iron ↔ eddy stagger converges in ~4 iters through deep
   saturation (plate loss within ~2–3% of an all-in-one FEM), the salient-bar motor torque-angle matches
   an exact-Newton FEM (mean 0.58%), and two induction machines pass (rotating cylinder vs Bessel 0.19%;
   mini cage 0.57%).  Remaining: a PRODUCTION `radia.` coupling API (the docs layer is a research helper),
   and the multi-harmonic AGE comparison bench (the linear-cylinder dual-lane bench already passes: VIM
   0.01% @ 502 dof vs AGE 7e-4% @ 31.6k dof).
5. The Sauter–Schwab 6D inner quadrature is a **negative result** so far (prototype plateaus at
   2e-3; needs a rigorous per-shuffle CPS second transform).  The retired Gauss-point and H-LU
   paths remain outside the public contract.

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
| RT1 pure-TET / pure-HEX public contract | `Solve(..., order=1)` | `test_hdiv_vim_rt1_contract.py`, `test_hdiv_vim_hex_public_solve.py` |
| Radia field write-back / IMA field contract | `rad.Solve(..., image=...)`, `rad.Fld(...)` | `test_hdiv_radfld_contract.py` |
| HACApK Gram build stats / speed guard | `_ChargeGramHMatrix`, `res["hmat_stats"]` | `test_hdiv_hacapk_gram_performance.py` |
| Planar motor minimal contract | `PlanarDemagBody`, angle sweep | `test_hdiv_motor_minimal_contract.py` |

The live executable checks are under `validation_test/feec/`; the legacy helper
scripts imported by those checks are under `validation_test/feec/vim_legacy/`
after the prototype retirement.

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
hex VIM lives under `validation_test/maglev/research_cln/ngsolve_validation/dd_*` — a
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
