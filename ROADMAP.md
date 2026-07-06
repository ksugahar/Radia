# Radia Research Roadmap

This file is tracked. Keep it aligned with the public repository policy and
with the agent-facing guidance in `AGENTS.md` / `CLAUDE.md`.

---

## Current Direction (2026-07)

### Magnetic Materials: HDiv-VIM First

- The production soft-iron and nonlinear magnetic-material path is HDiv-VIM.
- The main engineering reason is reduced-FEM and NGSolve interoperability:
  mesh labels, FE spaces, TaskManager parallelism, high-order/curved meshes,
  and open-boundary coupling all stay in one mathematical vocabulary.
- Future acceleration work should target HDiv-VIM matrix assembly,
  ChargeGram/HACApK, preconditioning, and reduced-FEM coupling.
- Retired moment-demagnetization prototypes are not roadmap items. Keep the
  source tree focused on the HDiv-VIM path and its validations.

### PEEC / BEM / ESIM

- PEEC remains the circuit-extraction path for conductor geometry, skin and
  proximity effects, and SPICE-facing reduced models.
- BEM and ESIM remain surface-integral / surface-impedance methods used by the
  PCB and induction-heating application domains.
- Conductor-core coupling should connect PEEC/BEM source fields to HDiv-VIM or
  reduced-FEM magnetic-material solves.

### Hysteresis And Nonlinear Materials

- Extend the current material model around the HDiv-VIM solve path.
- Keep the Play-model and rate-dependent loss terms compatible with energy-based
  formulations where symmetry and solver robustness matter.
- Record nonlinear benchmark claims under `validation_test/` with mdx timing
  provenance when the run is heavy.

### Field Reconstruction For Beam And Motion Workflows

- Use NGSolve HDiv / HCurl spaces to reconstruct fields that satisfy Maxwell
  constraints, especially for beam trajectory, accelerator magnet, motor, and
  MagLev workflows.
- Preferred route:

  ```text
  Radia source or HDiv-VIM solve
      -> NGSolve CoefficientFunction / GridFunction
      -> divergence-compatible interpolation / trajectory computation
  ```

- Keep particle and motion workflows batch-friendly with NGSolve TaskManager and
  durable JSON result records.

### Cubit / Mesh Export / Visualization

- Cubit export focuses on Coreform Cubit 2025.12 and high-order curved mesh
  quality.
- `.geo` is the standard human-facing Gmsh post route; `.msh v4.1` remains the
  durable headless validation artifact.
- Mesh evaluation demonstrations belong in result-bearing docs notebooks, not
  in transient GUI menu commands.

---

## Near-Term Backlog

- HDiv-VIM API naming cleanup to match NGSolve style.
- HDiv-VIM 2D support and focused tests for dispatch, `rad.solve`, and `rad.Fld`.
- ChargeGram/HACApK review: symmetry, image-method consistency, and assembly
  timing across small, medium, and high-order/curved meshes.
- PEEC plus HDiv-VIM conductor-core coupling examples with reproducible JSON.
- Panel migration cleanup: keep notebook workbenches as the current operating
  surface, with Cubit-only menu/toolbars handled by `cubit-mesh-export`.
- `examples/` retirement sweep: promote useful material to `src/`, `tests/`,
  `validation_test/`, or result-bearing `docs/*.ipynb`.
