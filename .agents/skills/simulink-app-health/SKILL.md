---
name: simulink-app-health
description: Verify Radia's production Simulink application blocks for EM, PCB, Motor, Stream Function, and IH. Use after changing buildLibrary, applicationSFunction, runApplication, application.py, a DesignSpec/calc script, the application interface manifest, or the tracked Simulink library; also use when a block does not run or leaves incomplete artifacts.
---

# Simulink Application Health

Radia's human production interface is one Simulink library. EM, PCB, Motor,
Stream Function, and IH are Simulink-only. Notebook workbenches are retired.

## Backend Boundary

- A block launches one batch solve only on an explicit rising-edge trigger.
- The initial production backend is the validated Python `DesignSpec` plus
  `calc_*.py` path through `radia.simulink.application`.
- Never call Python at every simulation step.
- MEX/ROM is optional. Do not make it the default until numerical parity,
  errors, state lifecycle, repeated runs, and long-run stability are tested.
- IH is the native-only exception. Its first standalone release is a preview
  runtime for checked preassembled operators. Do not call it production-complete
  until Cubit `.vol` to PEEC/BEM-A/BIM/FEM operator assembly and strict label
  contracts pass on the supported machines.
- Solver logic stays in tested Radia APIs; the block is orchestration only.

## Required Gates

Run from the repository root:

```powershell
python -m pytest tests/test_simulink_application.py -q
python -m pytest tests/test_application_interface_manifest.py -q
python tools/audit_new_panel_contract.py
```

Use the official MATLAB MCP Server for MATLAB checks:

1. Run Code Analyzer on `buildLibrary.m`, `applicationSFunction.m`,
   `runApplication.m`, and `writeApplicationConfig.m`.
2. Run the focused application-block tests in
   `tests/matlab/test_simulink_workflow.m`.
3. Build the library in `C:\temp`, load it, and update a model containing one
   copied application block connected to a Boolean trigger source.
4. Verify Electromagnet, PCB, Motor, and Stream Function have masks and use
   `radia_application_sfun`. Verify IH is a masked subsystem with separate
   `radia_ih_eddy_sfun` and `radia_ih_thermal_sfun` native blocks, closed
   temperature feedback, and no Python fallback.
5. For an IH preview release, run the extracted-package verifier through
   `release_qud simulink-candidate` on LAB, 100号機, mdx, and hibino.

## Artifact Contract

For success and failure, require a run directory containing:

```text
launcher_command.txt
command.txt        # when DesignSpec command construction succeeds
run.log
result.json
solver_result.json  # when the solver reached output generation
<application>_fields.msh  # required for a spatial-field mode
```

`result.json` must identify the application, backend, status, timing,
versions/platform, configuration hash, primary scalar or error, and GMSH
artifact list. A spatial-field mode must expose `--msh-output`; the runner
redirects it into the run directory and rejects anything except `.msh v4.1`.
Scalar/circuit-only modes record GMSH as not applicable.

Verify every generated GMSH file before opening it:

```powershell
python tests/panels/verify_gmsh_output.py <run-dir>\<application>_fields.msh
```

After the gates pass, regenerate the tracked library:

```matlab
radia.simulink.buildLibrary(OutputDirectory=fullfile(pwd, "matlab"));
```

Check Cubit's export toolbar separately with its deploy/smoke tools.
