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

For axisymmetric-to-3D force checks, use the public-safe gate in
`radia_mcp.radia_ngsolve.axisymmetric_3d_validation`.  An axisymmetric result
from `eggshell_force_axi` is already a full `2*pi*r` revolution quantity, so a
3D force artifact should record its vector basis (`full_revolution` or
`symmetry_sector`), sector angle, axis convention, solver versions, and timing
before comparing the axial force and transverse cancellation.  The executable
validation example is:

```powershell
python validation_test/force_validation/validation_axisymmetric_to_3d_force_gate.py
python validation_test/force_validation/validation_axisymmetric_to_3d_vol_force.py
```

The second script writes a real solver-run artifact and a solver-ready
artifact: `validation_axisymmetric_to_3d_vol_force_summary.json` and
`validation_axisymmetric_to_3d_vol_force_ready.json`.  It generates or reuses a
Netgen `.vol`, reloads it with NGSolve, integrates `J x B` over the 3D target
torus, and then runs the same axisymmetric-to-3D gate.  It also writes the
equivalent Cubit journal so a Coreform Cubit `export netgen` mesh can replace
the default Netgen/OCC mesh when a Cubit license is available.
