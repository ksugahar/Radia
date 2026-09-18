# Fin-capable PEEC: first conservative topology slice (2026-09-19)

This is an internal implementation/validation ledger, not a user workflow.
The production IH PEEC path remains the existing series-filament bundle.
No CLN code or CLN-derived snapshots are used by this work.

## Implemented

- `radia.peec_fin_topology` accepts CAD-derived perimeter station rings.
- Longitudinal branches exist between every pair of stations. Transverse
  branches exist only at selected interior stations. Every lane retains its
  own transition node; no artificial equipotential junction is introduced.
- The STEP adapter accepts only direct perimeter extraction
  (`step_uv`, `step_per_station`, `step_section_planes`). Equivalent-circle
  fallbacks are rejected. For an unknown section, the first-pass auto rule
  meshes every interior station; selective mesh removal is not yet claimed.
- A small dense MNA reference solve and a C++ PEEC assembly experiment both
  preserve node KCL on the hybrid graph. The terminal is a logical common
  node while each branch retains its physical CAD endpoints for L assembly.

## Not yet accepted as physical IH fin solver

The C++ experiment approximates a surface branch as a thin rectangular
segment. Its partial self/mutual inductances, internal-reactance split, and
surface resistance at a sharp fin have **not** been validated. There is no
SIBC coupling, no frequency/temperature-dependent skin-depth selection, no
thin-fin two-face coupling, and no branch-current back-reaction to the
workpiece BEM. It must not be enabled in `calc_inductance.py` or the Simulink
IH operator assembler on the strength of topology/KCL tests alone.

## Required next gates

1. Obtain a representative beak-fin STEP and preserve its tip radius. Check
   direct CAD extraction, station alignment, branch geometry, and convergence
   under perimeter/axial refinement. Reject missing tips and zero-area cells.
2. Implement a consistent surface-current discretization with passive SIBC
   Gram matrix and validated nonorthogonal partial elements. Avoid adding the
   existing isolated-wire Dowell/Bessel or proximity correction on top.
3. Compare local current, tip/root loss, terminal impedance, and field at the
   workpiece against an independent 3-D A-phi/HCurl reference over frequency,
   conductivity, geometry, and mesh sweeps. KCL alone is not an accuracy gate.
4. Extend the workpiece weak and strong coupling APIs to accept per-branch
   currents and per-branch induced EMFs. The present K-by-K bundle reduction
   and one-current-per-filament field projection cannot represent transverse
   branches or longitudinal redistribution.
5. Only after those gates, add an explicit `fin-surface` IH backend and
   document its contract through the owning Radia MCP manual. Existing
   rectangle/circle production behavior must remain unchanged.
