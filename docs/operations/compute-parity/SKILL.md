---
name: compute-parity
description: Deploy and verify identical Radia wheel environments on mdx1, mdx2 and hibino, then execute committed electromagnetic analysis jobs with hash-verified recovery and cleanup. Use before cross-host numerical comparisons or after solver workflow improvements.
---

# Identical compute environments

Use `tools/compute_runtime.py` and `tools/cleanup_compute_job.ps1` from the
same reviewed local Radia commit. Customer CAD, material tables and results
stay in the customer's private case repository; the shared repository owns
reusable solver implementations, operation helpers and independent validation.
Do not publish customer data or push merely to make a compute host current.

## Establish one numerical implementation

Inspect local WIP and preserve other sessions' edits. Commit the owned workflow
changes in an isolated worktree. Fix solver defects in their owning repository,
with an appropriate regression, before deployment. Do not inject alternate
functions into `radia.vim._solve`, execute AST snapshots as solver replacements,
or patch installed Python/native files during a customer comparison. An explicitly
declared experimental formulation is separate evidence, not a production result.

If solver code changes, read the `build` skill and build one candidate wheel
from the reviewed clean commit with current `Build.ps1`. Preserve its source
and native build provenance. Deploy those same wheel bytes everywhere; a
version string or editable checkout alone is insufficient. If only the operation
workflow changes, an already verified released wheel can remain the numerical
implementation. Do not claim an unknown wheel source commit from its version.
This procedure does not tag, publish or perform a release-quad release.

Resolve a complete Windows/Python-compatible wheelhouse using package metadata.
Select the intended dependency versions once, then retain all resolved wheel
bytes. Include pytest only when the declared acceptance job needs it. Do not
install Radia MCP, Cubit, Git or developer tools on compute hosts for this flow.

Create the lock after committing the execution helper:

```powershell
python tools/compute_runtime.py lock --wheelhouse C:/temp/compute-stage/wheelhouse --output C:/temp/compute-stage/lock.json --procedure-commit <reviewed-SHA>
```

The lock records the Python patch version, every wheel hash, dependency version,
critical source/native hashes, thread settings and operation commit. Transfer
the same helper, lock and wheelhouse to each host's owned `C:\temp` staging
directory. Verify the staging archive hash on arrival. Use foreground SSH as in
the `mdx-compute` skill. Never start a detached computation.

Deploy on each host with its bootstrap Python of the required patch version:

```powershell
python C:/temp/compute-stage/compute_runtime.py deploy --lock C:/temp/compute-stage/lock.json --wheelhouse C:/temp/compute-stage/wheelhouse --runtime C:/ProgramData/Radia/compute-runtimes/<first-16-lock-hash>
```

These versioned virtual environments are maintained runtime deployments,
separate from temporary jobs and existing CI/default Python. The command uses
offline, hash-checked wheels, checks dependency closure, verifies installed
RECORD hashes and imports, and writes `acceptance.json`. Compare all hosts'
lock hashes, package/module hashes, shared-file winners and loaded numerical
DLL hashes. NGSolve and Netgen share a delvewheel initialization file; its
installed bytes must be one of the hashes in the locked wheel records and
must match across hosts. The launcher explicitly sets MKLROOT to the selected
runtime's absolute Library directory and rejects numerical DLLs loaded outside
that runtime. Never substitute default `python` for
the recorded runtime's `Scripts/python.exe`.

## Admit and execute a job

Probe hibino via SSH, not ping; use it first only when already started and idle.
Otherwise select an idle mdx host. Inspect CI runner/queue state from the
controller and Python/MATLAB jobs and free memory on the host immediately before
execution. CI/preflight has priority. The launcher also rejects active CI workers
and, for analysis jobs, other Python/MATLAB processes. A small runtime acceptance
job may declare `workload: smoke`; it is not a heavy performance measurement.

Read [job-spec.md](references/job-spec.md) for bundle fields. `pack` exports
program files from their declared Git commits and hashes every input. Use a
fresh owned case directory for each run. Never stage uncommitted numerical
patches or stale model inputs. For Netgen meshes, record the `.vol` generation
source and versioned label contract; the case entrypoint must pass check-vol
before solver initialization. Mesh checking is Cubit-independent and does not
authorize installing or executing licensed CAD software.

Declare SI units, current orientation/mirror signs, material interpolation,
outer boundary, mesh/order/quadrature and linear/nonlinear tolerances in the
physical contract and in the actual entrypoint. Reject unsupported choices;
do not smooth a B-H table or change the constitutive route silently. For HDiv-MMM,
use public `vim.HDivSolver`/`Solve` with caller-owned `ngsolve.TaskManager` and
SparseCholesky. Confirm true linear residual <=1e-6 and the declared nonlinear
residual gate. A process exit code or a settled Newton step is not convergence.

```powershell
python tools/compute_runtime.py pack --lock C:/temp/compute-stage/lock.json --spec C:/temp/compute-stage/job-spec.json --output C:/temp/compute-stage/job.zip
```

For runtime acceptance, bundle `tools/compute_smoke.py` and the two existing
test modules named there. Run exactly this bundle on all three hosts. Preserve
JSON and logs and require both tests to pass. This proves installed-runtime
acceptance, not customer-model accuracy or a formulation's superiority.

After hash verification and unpacking, run in foreground:

```powershell
C:/ProgramData/Radia/compute-runtimes/<lock-id>/Scripts/python.exe -I C:/ProgramData/Radia/compute-runtimes/<lock-id>/compute_runtime.py run --lock C:/ProgramData/Radia/compute-runtimes/<lock-id>/lock.json --job-root C:/temp/<owned-job>
```

The launcher rechecks the environment and input/code hashes, clears Python
source overrides, fixes thread settings and records runtime, command, exit,
time and memory. The case must save its own numerical result JSON including
convergence and comparison criteria. Keep host hardware/context distinct from
software identity. Do not compare accuracy until discretization, source,
material law and boundary differences are explicitly accounted for.

## Recover before deleting

After the foreground process has ended, archive the owned job with `recover`.
Copy the archive to durable private storage, then use `verify-recovery` to check
the archive SHA and every recovered file. Retain failed-run diagnostics too.
Commit intended numerical JSON/code/provenance in the private case repository;
raw logs, meshes and binary staging remain durable outside tracked source.

Run `cleanup_compute_job.ps1` only after verifying that evidence commit. Supply
the literal job root, recovery archive, confirmed archive hash, full evidence
commit SHA and any owned staging files. The helper rejects escaping targets,
links/junctions and active owned Python jobs, deletes only those targets, and
verifies absence. Recover the cleanup certificate to private storage. Remove
deployment staging after the matching wheel/lock hashes and acceptance evidence
are durably preserved. Keep only the intentionally deployed runtime on hosts,
never case data, disposable job environments or completed job archives.

Report the reviewed commit, runtime lock, three-host acceptance, evidence
locations and cleanup state. If the environments differ, stop admission and
repair/deploy the common candidate before any comparative computation.
