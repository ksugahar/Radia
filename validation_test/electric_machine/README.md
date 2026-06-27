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
```

The docs-facing archive is
[`docs/electric_machine/electric_machine_validation_archive.ipynb`](../../docs/electric_machine/electric_machine_validation_archive.ipynb).

The cogging/skew readable demo is now the result-saved
[`docs/electric_machine/cogging_skew_demo.ipynb`](../../docs/electric_machine/cogging_skew_demo.ipynb)
with its notebook-coupled helper beside it.
