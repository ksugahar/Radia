# Electric-Machine Validation Suite

Validation-class rotating-machine and electromagnetic force/torque checks.
These scripts were promoted from the retired electric-machine examples topic
after the docs archive existed; each script writes a matching `*_summary.json`
beside itself.

Run individual gates from the repository root, for example:

```powershell
$env:PYTHONPATH = (Get-Location).Path
python validation_test/electric_machine/validation_machine_torque_scaling.py
python validation_test/electric_machine/validation_planar_lorentz_block_force.py
python validation_test/electric_machine/validation_pmsm_angle_periodic_rom.py
python validation_test/electric_machine/validate_motor_rom_c_abi.py --library build-msvc/radia_motor_rom.dll
```

`validation_pmsm_angle_periodic_rom.py` builds an 8-pole/24-slot curved PMSM
cross-section and locks angle interpolation, positive inductance, and three
torque routes on interlaced holdout angles.  The solver-heavy reference run is
executed on an idle compute host; the adjacent summary records the actual host
and elapsed time.

`validate_motor_rom_c_abi.py` advances the same passive ROM through the Python
and C implementations for 1000 implicit steps.  It locks generalized currents,
torque including cogging and `v x B`, and the discrete power-balance residual.

The cogging/skew readable demo is now the result-saved
[`docs/electric_machine/cogging_skew_demo.ipynb`](../../docs/electric_machine/cogging_skew_demo.ipynb)
with its notebook-coupled helper beside it. The helper's standalone entry point
writes `cogging_skew_demo_results.json` here; the notebook embeds the public
presentation without creating a docs-local result sidecar.
