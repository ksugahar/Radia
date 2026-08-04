# Radia IH Simulink Preview

`radia_ih.slx` is the first native induction-heating runtime preview. It contains
separate readable Level-2 MATLAB Eddy and Thermal S-Functions backed by checked
`radia_mex` object handles, explicit current/angle/ambient source blocks, a
visible temperature-feedback loop, and an `IH Parameters` block. Initialization
uses the configured Python installation only to locate NGSolve's shared runtime
DLLs. Eddy and Thermal steps execute in MEX and never call Python per step.
This package targets 64-bit Windows and is generated and gated with MATLAB and
Simulink R2026a.

After extracting the release, start MATLAB in the extracted directory and run:

```matlab
model = install_radia_ih();
```

The saved model contains a one-DOF diagnostic configuration so installation
can be checked immediately. It is not a physical design result. For an actual
case, set `IH Parameters > IH configuration MAT/JSON` to a file containing
`config` or `radia_ih_config`. Use `radia.simulink.makeIHNativeConfig` when
constructing the checked native configuration. That function requires
preassembled Eddy/Thermal operators and successful strict-label `.vol` reports.

The `Geometry Update` block provides `Browse...` controls for the workpiece
Netgen `.vol`/`.vol.gz` mesh and the coil `.step`/`.stp` (PEEC) or
`.vol`/`.vol.gz` (BEM-A) input. After
the assemble command and its output MAT/JSON configuration are set once,
selecting a replacement file is enough: the next diagram update or simulation
start detects path and content changes, rebuilds the operators, and reloads the
configuration. `Rebuild now` forces the same operation explicitly.

This preview does not yet construct PEEC, BEM-A, BIM, or FEM operators from a
Cubit `.vol` file. That native assembly boundary is required before the model
can be called production-complete. LUT and lumped-state-space IH helpers are
not included in this package. Preparing a physical configuration currently also
requires the separately installed `check-vol` command for strict mesh-label
reports; the diagnostic model does not require it.

The discrete update order is fixed:

```text
Eddy at T(t)
  -> conservative workpiece transport from theta_prev to theta_now
  -> Thermal update to T(t + dt)
```

For a linear magnetic law, a current-only change rescales heat density without
refactorizing the Eddy operator. A temperature-dependent operator causes the
Eddy solve to update. Nonlinear BH iteration is not present in this preview and
`bh_mode="nonlinear"` fails during configuration instead of silently using a
linear solve.
