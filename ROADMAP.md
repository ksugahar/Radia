# Radia Roadmap

This file is tracked. Keep it aligned with the public repository policy and
with the agent-facing guidance in `AGENTS.md` / `CLAUDE.md`. Every item states
its status: **Achieved** (merged and backed by tests or validation records),
**In progress** (work under way, not yet accepted) or **Planned** (target
fiscal year). Do not mark an item Achieved without that evidence.

Updated 2026-10-24.

---

## Toolchain Position

### Post-processing: Gmsh 5 + MSH 4.2

- Gmsh is the post-processing tool. It renders fields well, including curved
  high-order elements, and every view is reproducible from options and scripts,
  so it can be driven from CLI and MCP prompts.
- Radia writes MSH 4.2 only, and only for post-processing. A `.msh` file is
  never a solver input.
- Status: **In progress.** The radia-mcp Gmsh tools, cubit-mesh-export and the
  Simulink spatial artifacts move to MSH 4.2 and Gmsh 5. Removal of the
  remaining Gmsh-to-solver routes and the radia-mcp knowledge update are under
  review.
- Next: write Gmsh `InterpolationScheme` data so high-order (p >= 2) fields are
  drawn adaptively inside each element, instead of element-node values only.
  Status: **Planned.**
- Viewing: `.vol` opens in the Netgen GUI and `.msh` in the Gmsh GUI through one
  radia-mcp launcher, usable from MCP prompts and from Windows file
  associations; `.vol` can also be converted to MSH 4.2 for viewing in Gmsh.
  Status: **In progress.**

### Pre-processing: Cubit, or build123d + Netgen

- Coreform Cubit (commercial, so some restrictions apply) is the strongest
  pre-processor. Its ACIS geometry lets cubit-mesh-export build high-quality
  curved high-order elements. Status: **Achieved.**
- On the current benchmark, cubit-mesh-export's geometric refit gives more
  accurate high-order geometry than Netgen's own high-order coefficients
  (`validation_test/cubit_mesh_export/geometric_refit_benchmark.py`).
- build123d + Netgen is the open-source route. It produces the checked `.vol`
  solver mesh without Cubit. Status: **Achieved.**
- The solver boundary is always a checked `.vol` with a versioned label
  contract, whichever route created it.

---

## Solver Capabilities

### HDiv-MMM

- HDiv-MMM uses NGSolve HDiv spaces, so high-order elements, curved geometry
  and Piola maps come from NGSolve. Status: **Achieved.**
- Topology optimization with Clebsch-hodograph design representation on top of
  HDiv-MMM (`docs/clebsch_hodograph/`). Completing and validating it is the
  first repository priority. Status: **In progress.**

### Fast Interactions

- HACApK compresses the core compact interactions; NGSolve BEM and its
  multipole facilities serve smooth free-space BEM and Biot-Savart work, so
  `GridFunction` sources are evaluated quickly. Status: **Achieved.**

### GPU And Supercomputer Execution

- GPU kernels built with MATLAB mexcuda, validated against the Python CPU
  reference. Status: **In progress** (first radia-simulink release pending).
- Supercomputer execution, CPU and GPU, on the Flow2 system. Target: October
  2026. Status: **In progress.**

### Control, Optimization And Simulink

- radia-simulink: masked blocks in one Simulink library for control and
  optimization studies. Deployment target includes the mdx platform, working in
  October 2026. Status: **In progress.**
- radia-optuna: a MATLAB implementation of Optuna, matched against
  upstream-generated fixtures. Status: **In progress.**

---

## Application Domains

| Domain | Target | Status |
|---|---|---|
| Induction heating (radia-ih) | FY2026 | In progress |
| Electromagnets (radia-electromagnet) | FY2026 | In progress |
| PCB (radia-pcb) | FY2027 | Planned |
| Motors (radia-motor) | FY2027 | Planned |
