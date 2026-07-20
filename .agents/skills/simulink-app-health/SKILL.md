---
name: simulink-app-health
description: Verify Radia's production Simulink application blocks for EM, PCB, Motor, Stream Function, and IH. Use after changing buildLibrary, applicationSFunction, runApplication, application.py, a DesignSpec/calc script, the application interface manifest, or the tracked Simulink library; also use when a block does not run or leaves incomplete artifacts.
---

# Simulink Application Health

Radia's human production interface is one Simulink library. EM, PCB, Motor,
and Stream Function are Simulink-only. IH temporarily keeps its notebook
workbench as a comparison surface.

## Backend Boundary

- A block launches one batch solve only on an explicit rising-edge trigger.
- The initial production backend is the validated Python `DesignSpec` plus
  `calc_*.py` path through `radia.simulink.application`.
- Never call Python at every simulation step.
- MEX/ROM is optional. Do not make it the default until numerical parity,
  errors, state lifecycle, repeated runs, and long-run stability are tested.
- Solver logic stays in tested Radia APIs; the block is orchestration only.

## Required Gates

Run from the repository root:

```powershell
python -m pytest tests/test_simulink_application.py -q
python -m pytest validation_test/panels/test_notebook_workbench.py -q
python tools/audit_new_panel_contract.py
```

The notebook test is an IH-only compatibility gate. It must also prove that
the retired EM, PCB, Motor, and Stream Function workbenches are absent.

Use the official MATLAB MCP Server for MATLAB checks:

1. Run Code Analyzer on `buildLibrary.m`, `applicationSFunction.m`,
   `runApplication.m`, and `writeApplicationConfig.m`.
2. Run the focused application-block tests in
   `tests/matlab/test_simulink_workflow.m`.
3. Build the library in `C:\temp`, load it, and update a model containing one
   copied application block connected to a Boolean trigger source.
4. Verify all five blocks have masks and use `radia_application_sfun`.

## Artifact Contract

For success and failure, require a run directory containing:

```text
launcher_command.txt
command.txt        # when DesignSpec command construction succeeds
run.log
result.json
solver_result.json  # when the solver reached output generation
```

`result.json` must identify the application, backend, status, timing,
versions/platform, configuration hash, and primary scalar or error.

After the gates pass, regenerate the tracked library:

```matlab
radia.simulink.buildLibrary(OutputDirectory=fullfile(pwd, "matlab"));
```

Check Cubit's export toolbar separately with its deploy/smoke tools.
