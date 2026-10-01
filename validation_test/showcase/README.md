# Selected application reruns

The six application notebooks read records in `results/`. Running those
notebooks needs Python, NumPy, Matplotlib and Jupyter; it does not start Radia,
NGSolve, CAD or an MCP server. Drivers here own numerical execution. Each
record names its runtime, execution date, checks and physical scope.

## Recompute on an available compute host

From the repository root, in an environment with Radia and NGSolve:

```powershell
python validation_test/showcase/run_selected.py complex_coil --output C:/temp/showcase-candidate/complex_coil.json
python validation_test/showcase/run_selected.py motor --output C:/temp/showcase-candidate/motor.json
python validation_test/showcase/run_physics.py winding --output C:/temp/showcase-candidate/winding.json
python validation_test/showcase/run_physics.py lift --output C:/temp/showcase-candidate/lift.json
python validation_test/showcase/run_physics.py particles --output C:/temp/showcase-candidate/particles.json
python validation_test/showcase/run_heating.py --output C:/temp/showcase-candidate/heating.json
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

## Recorded source and runtime

The six saved results are historical execution records from Radia 4.95.90 /
NGSolve 6.2.2606. Their timestamps and hashes have not been updated merely
because documentation or a driver changed. They do not validate the current
release. `campaign.json` records hashes of the sources actually executed,
not a promise that today's checkout has identical bytes.

An audit against the recovered execution archive found a substantive difference
between the executed `radia_mcp/radia_ngsolve/solve.py` and the PR #302 source:
the executed helper included residual checks and nonconvergence failures absent
from that published revision. The sphere helper and `machine_scaling.py` had
newline-only differences. Therefore exact reproduction of the motor result from
that public revision is **not established**. The recorded checks remain evidence
for the recorded run only; rerun and validate on the selected public revision
before using it as release evidence. No execution-source hashes were rewritten
to hide this difference.

Some older case records do not include Python, helper hashes or thread counts.
The readers label missing metadata instead of inferring it. Heating's additional
recorded dependencies are in `campaign.json`. Updated drivers record a common
runtime header for future runs; unknown thread counts remain null. The accepted
JSON files are unchanged by these documentation fixes.

## Redraw the published figures

Install NumPy, Matplotlib, nbformat, nbclient and ipykernel in the viewer
environment, then run from the repository root:

```powershell
python docs/application_cases/render_previews.py
```

This executes the six JSON-only readers, updates their saved notebook outputs,
and extracts the gallery PNGs from those same outputs. It also extracts the
historical correction plot from the retained study and redraws the explicitly
labelled analytical thermal reference. It does not run a field solver.
`--output-dir` changes the PNG destination; notebook outputs are still updated.

## Isolate future computation outputs

The physics and heating drivers accept `--workdir` for intermediate files.
It must be outside the repository; the Windows default is a unique directory
under `C:/temp`. Only the explicitly selected final `--output` record is written
into the repository. The sphere driver passes its scratch directory to the
helper, so it does not overwrite another validation lane's accepted record or
the original docs plot. Review scratch logs and results before removing them.
These driver revisions have contract checks; the historical solver runs have
not been repeated with the revised drivers.

Execution hashes describe the original byte streams, including their line endings.
Source files may be normalized by Git; compare with the recorded execution bytes,
not a differently normalized checkout. Accepted `results/*.json` are marked
`-text` in `.gitattributes`, preserving those bytes on every checkout. New runs
should go to a candidate directory; adoption into `results/` is a separate review
that updates the accepted campaign intentionally.
