# Mesh Fusion in NGSolve -- proof of concept

Demonstrates the academic foundations of Ansys HFSS Mesh Fusion using
NGSolve: independently mesh two regions at different densities, then
couple them weakly across the shared interface.

## Background

See `radia_mcp.fem.fem_nonconforming_mesh_coupling(topic=...)` for
the full theoretical context (mortar / Nitsche / FETI-DP / DG / etc.).

## Examples (build incrementally)

### `phase1_nitsche_h1_poisson.py` -- Nitsche coupling, 2D Poisson

Two side-by-side sub-rectangles meshed at DIFFERENT `maxh`.  The
interface between them is non-conformal in the SENSE that the two
sides have different element sizes (one side has more interface
vertices than the other).  This is the simplest tractable example
that demonstrates the core challenge.

Solves: -nabla.(mu nabla u) = f  on a unit square domain split at
x = 0.5 into left/right sub-rectangles.  Each sub-rectangle has a
different mu (heterogeneous material).

Output:
- Convergence rates (h vs L2 error) per side
- Comparison against a uniform fine reference mesh
- JSON results -> `results_phase1.json`

### `phase2_mortar_h1_poisson.py` (TODO)

Same problem as Phase 1, but with proper MORTAR Lagrange multiplier
formulation.  Demonstrates LBB-stable saddle-point assembly + the
mortar projection between non-matching interface meshes.

### `phase3_hcurl_maxwell_mortar.py` (TODO)

The actually-on-point variant: 2D TM eddy current with HCurl
trial space.  Mortar tangential trace on the interface.  This is
the analog of what HFSS Mesh Fusion does for high-frequency EM.

## Folder convention

Mirror of `examples/hiruma_xfem_comparison/` -- phaseN_<topic>.py
+ results JSON, README for orientation.

## References (most relevant)

- Buffa-Maday-Rapetti 2001 (M2AN) -- 2D EM mortar standard
- Becker-Hansbo-Stenberg 2003 (M2AN) -- Nitsche-style coupling
- Egger et al. 2020 (arXiv 2005.12020) -- harmonic mortar (modern)
- Zhang-Liang 2020 (arXiv 2009.04400) -- transfinite mortar (HFSS-style)

All PDFs at W:/03_文献・論文/00_電磁界解析/15_非接合/.
