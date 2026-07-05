# C-Type Dipole Panel Fixtures

Small panel-owned fixtures rescued from the former
`examples/cubit_panels/accel_magnet` staging area.

- `coil_wire.step` is a STEP fixture used by the EM panel rejection test:
  normal EM coil input must be a Python `build_coil()` script; STEP coil input
  is reserved for PEEC-style workflows and should fail with a clear message.
- `yoke.step` is the historical C-type dipole yoke STEP used by the
  `step_mesh_builder.py` source-check path.

The old exploratory Python scripts were removed from `examples/`; their lesson
is recorded in `memory/cubit_panels_accel_magnet_prune.md`.  Reusable APIs
belong in `src/`, while panel-only fixtures belong under this `panels/` tree.
