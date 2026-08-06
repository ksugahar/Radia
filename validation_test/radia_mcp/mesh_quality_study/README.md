# Mesh-quality study: netgen vs cubit (equal-h and equal-budget)

Promoted 2026-08-06 from `C:/temp/mesh_quality_study`. One STEP per
geometry (lab `build123d` helpers), meshed by netgen (tet) and the Cubit
batch daemon (tet + hex where sweepable), all judged by ONE referee:
gmsh **minSICN** via `radia_mcp.gmsh.msh_inspect.mesh_quality` on
`.msh v4.1`. First-order elements throughout (shape quality is a
corner-geometry property — 2026-08-06 decision).

Two comparisons per case:

* **equal-h** — both meshers get the same target size (netgen produces
  ~half the elements at equal h);
* **equal-budget** — netgen `maxh` calibrated (n ∝ 1/h³ iteration) until
  its element count matches cubit_tet within ±10 %, answering *"same
  element budget, which mesher is better?"*

`results_mesh_quality_study.json` is the committed result (Data
Persistence Policy); `run_study.py` regenerates it (~80 s on LAB;
requires Cubit + netgen + gmsh + build123d; scratch meshes land in
`artifacts/`).

## Measured findings (LAB, Cubit 2025.12 / netgen 6.2.2604, 14 cases)

1. **Cubit tet wins the worst-element metric in every case, equal-h AND
   equal-budget** (e.g. sphere 0.59 vs 0.41; c_core 0.54 vs 0.12;
   holed cylinder 0.56 vs 0.41). Equalizing the element budget narrows
   the MEAN quality gap (both ~0.82–0.86) but not the tail: the
   difference lives in the worst elements, and cubit_tet's minimum
   stays in a tight 0.50–0.65 band across all geometries and sizes.
2. **Netgen's minimum quality is NOT monotone under refinement**:
   calibrating c_core from maxh 8.0 to 6.97 *dropped* min from 0.47 to
   0.12 (new slivers appear as h shrinks). Do not assume "finer netgen
   mesh ⇒ better worst element".
3. **Netgen is ~2× more economical at equal h** (half the elements for
   the same size target) — its strength is economy and speed, not the
   worst-element tail.
4. **Hex dominates wherever the geometry sweeps**: thin plate min 0.93,
   stacked gap pair min = mean = **1.000** (perfect lattice), with
   5–8× fewer elements than any tet route.
5. **No inverted elements anywhere** (all routes, all 14 cases,
   `negative = 0`).

Practical rule distilled for the lab pipeline: prefer Cubit hex where
sweepable; for tet, Cubit tetmesh gives the safer worst-element floor,
netgen the cheaper mesh at equal h — and always judge by the shared
minSICN referee, not by element count.
