# Mesh-quality study: netgen vs cubit

Promoted into the validation lane on 2026-08-06. One STEP per geometry
(lab `build123d` helpers), meshed by netgen (tet) and the Cubit
batch daemon (tet + hex where sweepable), all judged by ONE referee:
gmsh **minSICN** via `radia_mcp.gmsh.msh_inspect.mesh_quality` on
`.msh v4.1`. First-order elements throughout (shape quality is a
corner-geometry property — 2026-08-06 decision).

Every script writes its committed results JSON beside itself (Data
Persistence Policy) and its scratch meshes into `artifacts/`
(gitignored). All five require Cubit + netgen + gmsh + build123d; all are
quality-class runs (no timing), so LAB execution is allowed.

| Script | Question | Result JSON | Runtime |
|---|---|---|---|
| `run_study.py` | which mesher wins, equal-h and equal-budget? | `results_mesh_quality_study.json` | ~80 s |
| `run_nonmonotone_sweep.py` | does refinement always improve the worst element? | `results_nonmonotone_sweep.json` | ~45 s |
| `run_size_target_sensitivity.py` | is `min` a reproducible mesher property at all? | `results_size_target_sensitivity.json` | ~60 s |
| `run_solver_impact.py` | does any of it change the solver? | `results_solver_impact.json` | ~10 s |
| `run_accuracy_per_dof.py` | which mesh is more accurate for the same **dof**? tet vs tet vs hex | `results_accuracy_per_dof.json` | ~60 s |
| `run_vector_elements.py` | does any of it hold for **HCurl / HDiv**? | `results_vector_elements.json` | ~110 s |
| `run_rotation_control.py` | is the hex win just **axis alignment**? | `results_rotation_control.json` | ~20 s |

`run_study.py` compares **equal-h** (same target size; netgen produces
~half the elements) and **equal-budget** (netgen `maxh` calibrated by
n ∝ 1/h³ until its element count matches cubit_tet within ±10 %).
`run_solver_impact.py` poses one manufactured-solution Poisson problem —
u = sin(kx)sin(ky)exp(kz), H1 order 1, Dirichlet on every boundary — on
each mesh and reports CG iteration counts (unpreconditioned and Jacobi)
plus relative L2/H1 error against the exact solution.
`run_accuracy_per_dof.py` turns that single point into a **convergence
curve**: the same manufactured problem over ~1 decade of ndof for every
route, compared at matched ndof. Element count is not the cost that
matters — the linear system is sized by **dof**, so that is the axis.

`results_mesh_quality_study.json` is the committed result (Data
Persistence Policy); `run_study.py` regenerates it on an idle compute
host with Cubit + netgen + gmsh + build123d. Scratch meshes land in
`artifacts/`.

## Measured findings (LAB, Cubit 2025.12 / netgen 6.2.2604)

1. **Cubit tet holds the higher worst-element BAND across the whole size
   range**, not just at one point. Over 24-point maxh sweeps:
   c_core netgen 0.272–0.546 vs cubit 0.538–0.593; sphere netgen
   0.238–0.612 vs cubit 0.571–0.616. `run_study.py` confirms it case by
   case: cubit_tet's `min` beats netgen's in **28/28** comparisons
   (14 geometry/size cases × equal-h and equal-budget). Equalizing the
   element budget narrows the MEAN gap (both ~0.82–0.86) but not the
   tail.

2. **`min` quality is CHAOTIC in the size target — for BOTH meshers.**
   At a *fixed* target both are bit-exactly reproducible (5/5 identical
   meshes, spread 0.0000). But moving the target by **0.125 %** — far
   below any meaningful refinement step — moves `min` by **0.118**
   (netgen, c_core 0.2817→0.3997) and **0.086** (cubit_tet, c_core
   0.5155→0.6014). Two consequences:
   * refinement does **not** monotonically improve the worst element —
     8 of 23 c_core steps and 7 of 23 sphere steps make it worse, with
     the largest single drop 0.132. This is a symptom of the chaos, not
     a separate netgen pathology;
   * **a single `min` number is one sample, not a mesher property.**
     Compare bands over several sizes; never rank two meshers on one
     call each. (This caveat is now in the
     `cubit_netgen_quality_compare` docstring.)

3. **At these levels `min` does not predict solver behaviour.** Order-1
   Poisson, identical problem on each mesh:

   | pair | mesher | min | ndof | CG | CG+Jac | L2 rel | H1 rel |
   |---|---|---|---|---|---|---|---|
   | sphere equal-budget | cubit_tet | 0.567 | 329 | 34 | 31 | 3.75 % | 23.73 % |
   | | netgen | 0.470 | 502 | 34 | 32 | 3.66 % | 23.38 % |
   | c_core equal-budget | cubit_tet | 0.593 | 2236 | 42 | 40 | 0.20 % | 5.73 % |
   | | netgen | 0.391 | 2571 | 37 | 35 | 0.28 % | 6.82 % |
   | c_core, netgen refining | netgen | 0.512 | 313 | 14 | 14 | 1.15 % | 14.01 % |
   | | netgen | 0.374 | 449 | 16 | 15 | 0.98 % | 13.18 % |
   | | netgen | 0.272 | 573 | 17 | 17 | 0.88 % | 12.68 % |

   A 0.20 deficit in `min` left CG *unchanged* on the sphere and made it
   *faster* on c_core; and watching one mesher's worst element collapse
   0.512→0.272 under refinement, the error improved monotonically
   throughout. The worst element left no trace in either observable.
   **The metric that actually matters is `negative > 0`** (inverted
   elements); `min` is a floor/safety indicator, not a solver predictor.

4. **Cubit's real advantage is accuracy per dof — same convergence
   RATE, better constant.** Over ~1 decade of ndof the two meshers
   converge at the same order (L2: sphere −0.735 vs −0.733, c_core
   −0.730 vs −0.724; the P1 ideal is −2/3), so neither is
   asymptotically superior. What differs is the constant: at matched
   ndof netgen's L2 error is **1.29×** (c_core, ndof 1481) to **1.39×**
   (sphere, ndof 700) cubit's, and point-wise interpolation over the
   whole range gives 1.14–1.89×. Netgen's economy is per *element*
   (~2× fewer at equal h), never per *dof* — and dof is what sizes the
   linear system.

   **Mechanism: netgen spends a much larger share of its nodes on the
   boundary.** Interior-node fraction over the sweeps — cubit_tet
   19–77 %, netgen 2–49 %; on the thin plate netgen bottoms out at
   13 %. Fewer interior unknowns per dof is fewer unknowns doing work.
   (Beware: an *equal-element-count* comparison hides this and makes
   netgen look better — that is exactly how the earlier equal-budget
   sphere pair appeared to favour netgen, while it was quietly using
   1.53× the dofs.)

5. **Hex: wins the energy norm and the conditioning, not L2.** On the
   thin plate (8×8×1, sizes chosen so the thickness is resolved),
   cubit_hex vs cubit_tet interpolated to matched ndof:

   | ndof | CG hex/tet | L2 hex/tet | H1 hex/tet |
   |---|---|---|---|
   | 1764 | 0.54× | 0.75× | 0.75× |
   | 2704 | 0.54× | 1.10× | 0.90× |
   | 5445 | 0.55× | 1.05× | 0.88× |
   | 10086 | 0.56× | 1.03× | 0.87× |

   So at equal dof hex needs **~1.8× fewer CG iterations**, uniformly
   across the range, and is **10–25 % better in H1** (the gradient/energy
   norm — the one that matters when the answer is a field derived from
   gradients). **In L2 it is a wash** (±10 %). Its 5–8× element-count
   advantage does not shrink the linear system, whose size is set by dof.
   Its minSICN dominance (plate 0.93, gap pair 1.000 vs tet ~0.6) is
   therefore NOT an accuracy claim — it buys conditioning and energy-norm
   accuracy, not L2.

   netgen is clearly last on the thin plate: at ndof ≈ 2011 its L2 is
   4.36 % against cubit_tet 2.17 % and cubit_hex 1.89 % at the same dof
   (log-log interpolated) — roughly **2× worse than either**. Its
   convergence rate there could not even be fitted — only 2 of 5 points
   had a ≥ 20 % interior.

6. **No inverted elements anywhere** — all routes, all cases,
   `negative = 0`.

7. **Methodological trap, now fail-loud.** The first thin-plate attempt
   (50×50×1, sizes 4.0→1.0) produced meshes whose every node lay on a
   Dirichlet face: **zero free dofs**, no solve at all, and a flat ~28 %
   "error" curve that was pure interpolation error — it looked like three
   meshers tying. `run_accuracy_per_dof.py` now raises on `free_dofs == 0`
   and excludes points with < 20 % interior nodes from rate fits (two
   coarse netgen sphere points had silently corrupted an earlier fitted
   slope). Always check the interior fraction before believing a
   convergence rate.

8. **Vector spaces (HCurl / HDiv) — the verdict mostly carries, but hex
   gains and the honest metric changes.** Manufactured curl-curl+mass
   and div-div+mass problems at lowest order, on the same three
   geometries, compared at matched dof:

   * **cubit_tet still leads netgen, but by less and no longer
     everywhere**: ahead in **29 of 32** matched-dof comparisons,
     0.91–1.27× (vs 1.14–1.89× and 28/28 on scalar H1). All three
     exceptions are the same one — the HCurl **curl** seminorm on the
     thin plate, where netgen edges ahead by 2–9 %. The boundary-node
     penalty is diluted because HCurl/HDiv dofs live on edges and faces
     rather than nodes, though netgen's free-dof fraction is still the
     lower one (HCurl 30–59 % vs cubit 51–83 %).
   * **hex helps more here than in H1**, and differently per space —
     matched-dof, alignment-free: HCurl L2 0.28×, **curl seminorm
     0.66×**, CG 0.45×; HDiv L2 0.56×, **div seminorm 0.97× (a dead
     heat)**, CG 0.47×. So for eddy-current (HCurl) work hex buys a
     real ~1.5× in the curl seminorm and ~2.2× in iterations; for
     flux-space (HDiv) work it buys iterations and L2, but nothing in
     the div seminorm.

9. **The metric that was lying: L2 on an axis-aligned benchmark.**
   Before the control, hex looked 30× better than tet in HCurl L2 and
   10–18× in HDiv L2 — while the curl/div seminorms on the *same runs*
   showed near parity. That asymmetry was the tell. Rotating the FIELD
   30° about z and 20° about y (rotation commutes with curl and div, so
   `curl curl u = k²u` still holds exactly — same geometry, same meshes,
   same exact-solution machinery) collapsed hex/tet L2 from 0.017–0.031×
   to 0.39–0.44× (HCurl) and from 0.055–0.100× to 0.93–1.01× (HDiv, a
   dead heat), and inflated hex's HDiv CG advantage by ~7×. The
   derivative seminorms did not move at all across the rotation
   (0.658→0.657, 0.973→0.967).

   **Rule: never verify or benchmark a structured mesh with an
   axis-aligned separable exact solution, and treat an order-of-magnitude
   mesh-comparison gap as an artifact until a control says otherwise.**
   Catalogued as bug pattern
   `axis-aligned-manufactured-solution-fake-superconvergence`.

## Practical rule for the lab pipeline

* **Cost is dof, not elements.** Rank meshes on error-vs-ndof, never on
  element count — element count is what makes netgen look economical and
  hex look revolutionary, and neither survives the dof axis.
* **Cubit tetmesh over netgen for accuracy**: same convergence rate,
  1.14–1.89× lower L2 error at matched dof, because more of its nodes are
  interior. netgen's strengths here are speed and element economy.
* **Cubit hex where the geometry sweeps cheaply**: at equal dof, ~1.8×
  fewer CG iterations and 10–25 % better H1 on scalar problems, and
  more on vector ones — ~2.2× fewer iterations plus ~1.5× in the curl
  seminorm for HCurl. The exception is the HDiv div seminorm, where hex
  and tet tie. Worth the meshing effort when the solve is
  iteration-bound or the answer is a gradient/curl quantity.
* **Gate on inverted elements** (`negative > 0`), not on a single `min`
  value — and when comparing `min`, compare a **band** of sizes.
* **Check the interior-node fraction** before trusting any convergence
  number (finding 7), and **check a rotated control** before trusting an
  order-of-magnitude win (finding 9).
* `mesh_quality` now reports these axes directly — `mesh_stats.
  dof_estimate`, `mesh_stats.interior_node_fraction`, and per-type
  `aspect_ratio` — so an agent can apply this rule without re-deriving
  it, and `cubit_netgen_quality_compare` carries them per route.

## Scope of these conclusions

Findings 1–7 are order-1 elements on an isotropic scalar (Poisson)
problem; findings 8–9 extend the comparison to lowest-order HCurl and
HDiv, which is what closes the scope hole that findings 1–7 originally
carried. Still untested: anisotropic materials, high-order curved
elements, and nonlinear problems, where element shape enters the
inf-sup / interpolation constants differently. The c_core vector cases
are mass-dominated (k·L ≈ 1.6 gives k² ≪ 1), so their CG counts say
little about curl-curl conditioning — the sphere is the balanced case.
Timing was deliberately not measured anywhere (LAB is a contended
host); every observable here is deterministic.
