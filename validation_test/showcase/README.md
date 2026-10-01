# Selected application reruns

The six application notebooks read records in `results/`. Running those
notebooks needs Python, NumPy, Matplotlib and Jupyter; it does not start Radia,
NGSolve, CAD or an MCP server. Drivers here own numerical execution. Each
record names its runtime, execution date, checks and physical scope.

## Recompute on an available compute host

From the repository root, in an environment with Radia and NGSolve:

```powershell
python validation_test/showcase/run_selected.py complex_coil --output validation_test/showcase/results/complex_coil.json
python validation_test/showcase/run_selected.py motor --output validation_test/showcase/results/motor.json
python validation_test/showcase/run_physics.py winding --output validation_test/showcase/results/winding.json
python validation_test/showcase/run_physics.py lift --output validation_test/showcase/results/lift.json
python validation_test/showcase/run_physics.py particles --output validation_test/showcase/results/particles.json
python validation_test/showcase/run_heating.py
```

The motor driver imports the numerical helpers from `radia-mcp`; this is an
isolated Python calculation, not a deployed MCP service. The heating driver
also requires CAD STEP support and the declared ESIM input files. The recorded
rerun used an isolated cadquery-ocp 7.8.1.1.post1 / VTK 9.3.1 environment.
Follow [compute routing](../../AGENTS.md#compute-host-routing); run one heavy
job at a time. Drivers can also need SciPy and Matplotlib through their helpers.

| Case | Check | Limits |
| --- | --- | --- |
| Complex coil | Finite, nonzero field and reversal under reversed current | Consistency check, not independent accuracy certification |
| Motor ripple | Dominant order, near-zero mean and analytical skew factor | Two-dimensional example with declared stack scaling; no drive controller |
| Spherical winding | Continuous fit and connected-wire RMS | Fixed Z2 target and geometry; no universal manufacturing claim |
| Sphere lift | Frequency limits and modal reduction | Prescribed-gradient induced-dipole approximation; no stable levitation claim |
| Particle fan | Kinetic-energy conservation in a magnetic field | Saddle-coil dispersion; no quadrupole or fringe qualification |
| Heating | Positive integrated power and nonlinear convergence | Weak coupling; local loss density and transient temperature need separate checks |

The result readers intentionally show these focused cases. Wider historical
research studies remain accessible through the validation lanes and repository
history. Re-executing a reader redraws the recorded result; it does not renew
the numerical acceptance date. Failed numerical runs must not replace a passed
record. Publishing a new runtime version requires its own verification.
