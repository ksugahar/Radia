# Force validation scripts

This directory contains the heavy electromagnetic force/energy
cross-validation suite for `radia_mcp.radia_ngsolve.force`.

Run the heavy suite explicitly from `packages/radia-mcp`:

```powershell
python validation/force/validate_force_xval.py
```

The full force suite may take minutes because it runs real 3D FEM solves.

For loop-learning target extraction, first run the lightweight public-safe
electromagnetic-force pass:

```powershell
python validation/force/electromagnetic_force_target.py --source-json <autonomous-basic-learning.json> --out-dir <out-dir>
```

This selects `force_torque_motor` slots, attaches analytic rows for Lorentz
force, air-gap pressure, PM force-gap scaling, and dq torque decomposition, and
leaves commercial/source-tool solver runs as follow-up candidates.
