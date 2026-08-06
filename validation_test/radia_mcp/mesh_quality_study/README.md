# Mesh-quality study: netgen vs cubit

Promoted into the validation lane on 2026-08-06. One STEP per geometry
(lab `build123d` helpers), meshed by netgen (tet) and the Cubit
batch daemon (tet + hex where sweepable), all judged by ONE referee:
gmsh **minSICN** via `radia_mcp.gmsh.msh_inspect.mesh_quality` on
`.msh v4.1`. First-order elements throughout (shape quality is a
corner-geometry property — 2026-08-06 decision).

Every script writes its committed results JSON beside itself (Data
Persistence Policy) and its scratch meshes into `artifacts/`
(gitignored). All four require Cubit + netgen + gmsh + build123d; all are
quality-class runs (no timing), so LAB execution is allowed.

| Script | Question | Result JSON | Runtime |
|---|---|---|---|
| `run_study.py` | which mesher wins, equal-h and equal-budget? | `results_mesh_quality_study.json` | ~80 s |
| `run_nonmonotone_sweep.py` | does refinement always improve the worst element? | `results_nonmonotone_sweep.json` | ~45 s |
| `run_size_target_sensitivity.py` | is `min` a reproducible mesher property at all? | `results_size_target_sensitivity.json` | ~60 s |
| `run_solver_impact.py` | does any of it change the solver? | `results_solver_impact.json` | ~10 s |

`run_study.py` compares **equal-h** (same target size; netgen produces
~half the elements) and **equal-budget** (netgen `maxh` calibrated by
n ∝ 1/h³ until its element count matches cubit_tet within ±10 %).
`run_solver_impact.py` poses one manufactured-solution Poisson problem —
u = sin(kx)sin(ky)exp(kz), H1 order 1, Dirichlet on every boundary — on
each mesh and reports CG iteration counts (unpreconditioned and Jacobi)
plus relative L2/H1 error against the exact solution.

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

4. **Cubit's real advantage is accuracy per dof, and it comes from the
   whole distribution rather than the tail.** On c_core cubit reached
   0.20 % L2 with 2236 dofs where netgen needed 2571 dofs for 0.28 %
   (1.15× the dofs, 1.43× the error). On the sphere the errors are equal
   but netgen spends 1.53× the dofs. Netgen's economy shows up per
   *element* (~2× fewer at equal h), not per *dof*.

5. **Hex dominates wherever the geometry sweeps**: thin plate min 0.93,
   stacked gap pair min = mean = **1.000** (perfect lattice), with 5–8×
   fewer elements than any tet route.

6. **No inverted elements anywhere** — all routes, all cases,
   `negative = 0`.

## Practical rule for the lab pipeline

Prefer Cubit hex where the geometry sweeps. For tet, Cubit tetmesh gives
the better accuracy-per-dof and a higher worst-element floor; netgen is
the cheaper mesh per element at equal h. Judge with the shared minSICN
referee, over a **band** of sizes — and gate on inverted elements, not on
a single `min` value.

## Scope of these conclusions

Everything above is order-1 elements on an isotropic scalar (Poisson)
problem. Shape quality is expected to matter more for HCurl/HDiv vector
problems, anisotropic materials, and high-order curved elements, where it
enters the inf-sup / interpolation constants differently — none of that
is tested here. Findings 2 and 3 should not be extrapolated past that
boundary without measuring.
