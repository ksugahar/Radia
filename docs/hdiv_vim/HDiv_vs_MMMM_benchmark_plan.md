# HDiv-VIM vs MMMM Benchmark Plan -- Radia HDiv-Only Transition

Decision date: 2026-07-05.

Radia's production soft-iron demagnetization direction is **HDiv-VIM first**.
The MMMM implementation is being moved to `ELF_MAGIC@研究室版`, where it remains
useful as a laboratory / legacy / cross-check implementation. The Radia-side
benchmark should therefore answer a narrower question:

> Is HDiv-VIM the better Radia production backend when accuracy, Reduced-FEM
> coupling, and engineering-scale runtime are considered together?

It does **not** need to prove that HDiv is faster on tiny cases. At very low
DoF both methods are fast enough, so small-N timing is recorded only as a
latency sanity check, not as a method-selection criterion.

## Primary Motivation: Reduced FEM Coupling

The strongest reason to adopt HDiv-VIM in Radia is not that it wins every tiny
timing comparison.  It is that HDiv-VIM lives in the same NGSolve finite-element
world as the Reduced FEM workflows Radia needs for motors, conductors, and
open-boundary multiphysics.

| Coupling concern | Why HDiv / NGSolve is better than MMMM |
|---|---|
| Shared data model | HDiv magnetization is a NGSolve `GridFunction` / `CoefficientFunction` on a NGSolve `Mesh`; Reduced FEM uses the same mesh, material labels, integration rules, and `BilinearForm` machinery. MMMM lives in Radia object / collocation elements and must be sampled or projected before FEM can consume it. |
| Weak-form handoff | VIM iron fields can enter conductor / motor FEM as source `CoefficientFunction`s or projected fields, so coupling can be written as integrals. MMMM naturally gives element moments / `rad.Fld` samples, which are good for postprocessing but less natural as weak-form inputs. |
| FEEC compatibility | H(div), H(curl), H1, and L2 are part of one de Rham sequence. Divergence, normal trace, curl, and gradient constraints line up with Reduced FEM spaces; loop modes are `ker(B)` by construction. MMMM needs hand-crafted loop/co-loop logic outside the FE complex. |
| Curved / high-order geometry | NGSolve `mesh.Curve`, Piola maps, and `GetTrafo` are shared by HDiv and FEM. MMMM's flat-element / custom-element representation creates a geometry translation layer and loses the clean curved high-order path. |
| Solver/runtime integration | `TaskManager`, sparse FE assembly, preconditioners, material `GridFunction`s, and notebook/webgui inspection are already NGSolve-native. MMMM can be fast as a standalone VIM, but coupling it to Reduced FEM adds glue code and another ownership boundary. |
| Existing evidence | `docs/electric_machine/planar_vim_motor.ipynb` already demonstrates the intended split: laminated nonlinear iron via VIM, conducting bars / rotor via reduced complex FEM, and staggered coupling against all-in-one FEM references. |

Therefore the benchmark should treat Reduced FEM coupling as a first-class
result: report field handoff quality, torque/loss agreement, stagger
convergence, and rebuild/remesh counts, not only matrix build time.

## Current Evidence

| Axis | Existing artifact | What it says |
|---|---|---|
| Curved accuracy | `hdiv_curved_showcase.ipynb`, `compare_curved_vs_radia_field.json` | Curved HDiv beats flat Radia/MMMM-style faceting by roughly 10-30x in field accuracy per resolution. |
| Curved nonlinear field | `hdiv_curved_nonlinear_field.json` | Flat geometry leaves about 9% field error; curved HDiv drops this below 0.4% in the showcased case. |
| Hex head-to-head | `hex_vs_mmmm_crossvalidation.ipynb` | On the same hex cube / C-yoke meshes, HDiv and collocation MMMM converge to the same continuum field; HDiv pins the cube demag factor at 1/3 and the gaps shrink under refinement. |
| 2D and Reduced-FEM link | `../electric_machine/planar_vim_motor.ipynb` | Planar HDiv-VIM supports motor cross-sections, nonlinear saturation, torque extraction, and reduced conductor coupling; representative errors are 0.19-1.28% on the documented checks. |
| Build scalability | `build_scaling_hdiv_vs_mmmm.ipynb`, `build_scaling_mdx_data.json` | HACApK charge-Gram build scales as about `n_charge^1.23`; HDiv build reaches MMMM build parity near 34k HDiv DoF and keeps compressing at larger N. |
| Application contract | `validation_test/feec/test_hdiv_radfld_contract.py` | `rad.Solve -> ObjSetM -> rad.Fld` is locked for HDiv write-back; materialized IMA fields are roundoff-consistent with the reduced solution. |

## Benchmark Matrix

| Lane | Purpose | Required metrics | Current / next artifact |
|---|---|---|---|
| Accuracy, analytic | Prove HDiv accuracy against cases with true solutions. | demag tensor error, field max/RMS error, p- or mesh-convergence. | Existing curved sphere / spheroid / ellipsoid notebooks and `validation_test/feec/test_hdiv_vim_*`. |
| Accuracy, same mesh | Show MMMM and HDiv approach the same continuum solution where both are valid. | volume-average M gap, external B gap, monotonic refinement trend. | Existing `hex_vs_mmmm_crossvalidation.ipynb`; next mdx refresh should extend hex sizes. |
| Reduced FEM coupling | Show HDiv is the natural Radia backend for motor / conductor workflows. | torque/loss error, stagger iterations, reusable field handoff, no remesh/rebuild per angle when possible. | Existing `docs/electric_machine/planar_vim_motor.ipynb`; next larger mdx nonlinear cube / motor run. |
| Engineering-scale speed | Check total runtime, not just build. | charge count, HDiv DoF, build time, solve time, iterations, memory, compression. | Existing tet mdx build scaling; next required run is pure-hex HDiv vs historical MMMM on mdx. |
| Small-N latency | Make sure the method is not awkward for interactive use. | wall time only, broad pass/fail. | Treat as "both fast"; do not use this lane to block HDiv-only Radia. |

## Interpretation Rules

1. **Accuracy wins before small-N speed.** If curved geometry, high-order RT
   charges, or Reduced-FEM coupling matter, HDiv is the Radia production route.
2. **Small-N speed is not decisive.** If both runs finish interactively, the
   result is "both fast"; the benchmark report should not promote MMMM because
   it saves milliseconds or sub-seconds on tiny meshes.
3. **Engineering scale is mdx-only.** Timing claims for build / solve / memory
   are measured on mdx, with host, Radia version, thread count, and JSON sidecar.
4. **MMMM is a transition reference in Radia.** Radia tests may keep MMMM only as
   a temporary cross-check while the deletion / ELF migration is audited. New
   user-facing Radia docs should lead with HDiv and analytic / Reduced-FEM
   evidence.
5. **Do not overclaim IMA.** HDiv materialized reduced-solution fields are
   roundoff-consistent. Explicit-full hex parity remains a separate reflection
   symmetry fix and must stay labelled until the xfail is closed.

## Exit Criteria For Radia HDiv-Only

- `rad.Solve(auto)` routes mesh-backed pure TET / HEX / WEDGE and 2D planar
  soft iron through HDiv-VIM.
- User-facing Radia docs present MMMM only as migrated / legacy / cross-check
  material, not as the recommended Radia backend.
- MMMM-specific tests, validation, and notebooks are either moved to
  `ELF_MAGIC@研究室版`, converted to HDiv/analytic references, or explicitly
  labelled as temporary migration gates.
- The mdx benchmark report contains at least one engineering-scale nonlinear
  HDiv run with build time, solve time, iterations, memory, and result JSON.
- The pure-hex build-scaling refresh is run on mdx after the current HDiv hex
  reflection-symmetry work is closed.
