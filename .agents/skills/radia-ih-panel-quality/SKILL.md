---
name: radia-ih-panel-quality
description: Compatibility entry point for Radia induction-heating interface QA. The retired PySide panel must not be revived; verify the production IH Simulink block with simulink-app-health, its DesignSpec/headless contract, native MEX handles, artifacts, and deployment.
---

# Radia IH Simulink Quality

This historical skill name is retained only so old prompts route to the current
production interface. Radia IH is a masked block in the single Radia Simulink
library. There is no standalone PySide application or notebook workbench.

Run `simulink-app-health` with IH in scope. Also use:

- `panel-cli-diff` for `IHDesignSpec` and headless calculation parity.
- `verify-deploy` to prove the loaded MATLAB, MEX, Python, and Cubit assets come
  from the intended checkout or release.
- `gmsh-verify` for the required spatial post-processing artifact.
- `radia-plugin-check` only when the issue concerns Cubit's independent toolbar
  or mesh export path.

The acceptance surface is the tracked IH block reopened through the official
MathWorks model lane. Verify mask parameters, typed ports, sample-time and
lifecycle behavior, independent Eddy/Thermal handles, `.vol` label checks,
`result.json`, `run.log`, and GMSH `.msh v4.1` output. Run focused unit tests in
`tests/`; keep solver-heavy, long-run, MATLAB/Simulink, and multi-machine
evidence in `validation_test/` with result JSON.

Do not recreate `src/radia/radia_ih.py`, `radia_gui_base.py`, Qt tests, a PySide
dependency, or a notebook launcher. Coreform Cubit's bundled PySide runtime is a
separate exception owned by the Cubit toolbar.
