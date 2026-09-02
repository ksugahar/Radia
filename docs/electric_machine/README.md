# Electric-Machine Docs

This docs topic now owns the human-facing electric-machine demo material that
used to live under `examples/electric_machine`.

| Artifact | Purpose |
|---|---|
| `cogging_skew_demo.ipynb` | Result-saved cogging / skew notebook. |
| `cogging_skew_demo.py` | Notebook-coupled helper for the finite-element torque sweep. |
| `cogging_skew_demo_results.json` | Computed torque, skew, checks, and runtime versions. |
| `planar_vim_motor.ipynb` | Executed planar HDiv-VIM machine showcase (2026-07-03): nonlinear iron deep-saturation vs the analytic fixed point; ellipse reluctance torque three ways with the demag operator built once per rotation sweep; salient-bar motor torque-angle vs an exact-Newton FEM (mean 0.58%); rotating-field conducting cylinder vs the Bessel closed form (0.19%); mini cage induction machine (VIM core + reduced-FEM bars) torque-slip vs an all-in-one FEM (0.57%). |
| `planar_vim_motor_helpers.py` | Notebook-coupled reference/coupling layer (exact-Newton nonlinear FEM reference, Bessel closed form, single-valued polar conjugate potential, cage stagger, frozen-secant all-in-one). |
| `planar_vim_motor_result.json` | Sidecar with the executed outputs, versions, and notebook sha. |
| `em_reference_audit.ipynb` | Executed reference-audit methodology (2026-07-04): the diagnostics that exonerate or convict a FEM cross-validation reference — coil-disk polygon current deficit (−5.4%) + the drive-equivalence probe (uniform 4.9%); the conjugate-potential gradient gate (1e-10) and atan2 branch cuts (72% jumps) vs the single-valued polar construction (closure 4e-15); the finite-Dirichlet dipole image matched by its closed form (0.00%) and the exact n=1 open Robin; the secant-Picard plateau vs the exact Newton (9 iters). Cross-linked to `bug_patterns_lookup(topic="validation")` and MCP `hdiv_vim(topic="reference_audit")`. |
| `em_reference_audit_result.json` | Sidecar for the audit notebook. |
| `angle_periodic_motor_rom.ipynb` | Executed angle-periodic motor-ROM study: 8-pole/24-slot curved PMSM cross-section, interlaced holdout angles, Maxwell-stress/virtual-work/ROM torque audit, and the production coupling contract. |
| `angle_periodic_motor_rom_result.json` | Synchronized benchmark metrics, runtime provenance, and notebook SHA-256. |
| `angle_periodic_motor_rom_torque.pdf` / `.png` | Publication and preview forms of the torque comparison. |

The solver behind `planar_vim_motor.ipynb` is the promoted 2D layer in `radia.vim`
(`PlanarDemagBody` / `Solve` on a 2D mesh / `maxwell_torque_circle`), golden-locked in
`validation_test/feec/test_hdiv_vim_2d_solve.py`.

The production reduced reluctance-motor API is `radia.motor_hdiv.HDivReducedMotor`.
It keeps the rotor mesh in a local frame, builds the symmetric BDM1 charge Gram
once, and reuses it over a mechanical-angle sweep.  The public torque contract
compares the air-gap Maxwell stress, the magnetization-volume coupling, and the
fixed-current coenergy derivative. The Motor Simulink block exposes the same
path as the **HDiv Reduced** study through
`src/radia/panels/calc_motor_hdiv_reduced.py`; its input mesh is a rotor-only 2D
`.vol`, not the full-motor mesh used by the transient A-formulation.

The production dynamic-machine API is `radia.motor_rom.AnglePeriodicMotorROM`.
It combines angle-periodic phase and reduced eddy-current ports with the
HDiv-MMM hysteresis restart contract, rigid-rotation `v x B` flux, cogging
coenergy, skew, end-winding corrections, temperature-dependent resistance,
thermal state, and an implicit electromechanical step.  The model bundle writer
`radia.motor_rom_export.SaveMotorROMBundle` emits NPZ, MAT, JSON, a MATLAB
loader, and an FMI model-variable fragment; `src/core/rad_motor_rom_c.h`
provides the deterministic C ABI used by external motion solvers.  The
Simulink boundary in `src/radia/simulink/` builds a fixed-step C-MEX
`radia_motor_rom_sfun` block from that same C ABI.  Its input is the phase
voltage/load/ambient-temperature vector and its output is the phase current,
eddy-current, flux, electromechanical torque, loss, temperature, and energy
diagnostic vector described in the adapter README.

The production block policy is deliberately layered: MagLev reduced models
export MATLAB `ss` state-space objects for the native Simulink State-Space
block, while dynamic machine ROMs use the C-MEX adapter or FMI Co-Simulation.
The HCurl eddy-bubble and HDiv-MMM transient adapters will use the same
physical-port contract after their time-domain state interfaces are finalized.

The curved 2D qualification uses 33 training and 33 interlaced holdout angles.
Its saved run reports a maximum holdout flux error of `1.67e-15`, a
Maxwell-stress versus ROM torque relative RMSE of `1.90e-3`, and a minimum
inductance eigenvalue of `0.137 H`.  These results establish the periodic-ROM
contract for that machine-equivalent benchmark; machine-specific 3D end-region
qualification remains a separate production acceptance task.

`PlanarDemagBody.field_cf(...)` is the native C++/NGSolve source-field bridge.
It rotates a solved source from its local frame into a target body's local
frame without a Python point loop.  This is the interface primitive for the
next fixed-stator reduced-FEM/AGE coupling increment; that full coupled machine
path is not claimed by the current single-rotor reduced solver.

The executable validation corpus is `validation_test/electric_machine/`; this
directory is the rendered, result-bearing docs layer.
