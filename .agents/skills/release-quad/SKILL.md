---
name: release-quad
description: Prepare, verify, and hand off CPU, Python GPU, and MATLAB release artifacts from one fixed source.
---

# Radia release

Start with CPU development. Advance the fixed source through the GPU and MATLAB
builds, verify the exact artifacts on the assigned machines, and hand them to
the release owner for publication.

| Order | Role | Work and evidence |
| --- | --- | --- |
| 1 | CPU build host | Develop and debug the CPU implementation. Run the shared precheck and the affected CPU checks. Freeze the source and build the CPU release candidate. Record its artifact hashes. |
| 2 | Linux GPU host | Build Linux Python GPU artifacts from the frozen source. Verify the artifacts on the intended device architecture. |
| 2 | Windows GPU host | Build Windows Python GPU artifacts and the MATLAB/Simulink package's CUDA MEX artifacts from the frozen source. Use the reviewed Windows CUDA, compiler, and MATLAB ABI. |
| 3 | Designated validation roles | Verify the exact hashes before execution. Run the applicable native, numerical, and package acceptance checks in isolated runtimes. |
| 4 | Assigned release owner | Verify the evidence against the publication scope. Publish and distribute the accepted artifacts under each package's release workflow. |

Resolve the machine mapping privately. The CPU build host stops at CPU
candidate preparation. Schedule CUDA compilation, GPU CI, and GPU execution
on the assigned GPU hosts. Run isolated exact-wheel CPU acceptance on the
validation host. Keep editable runtime maintenance separate. After freezing the CPU candidate source and dependency identities, the Linux and
Windows GPU builds may proceed in parallel. The coordinated handoff requires all
four lane receipts for the same source SHA.

## Acceptance and handoff

- Record the exact source and dependency identities, platform, compiler,
  runtime, architecture, artifact hashes, and the checks run on each artifact.
  Verify package contents and wheel RECORD where applicable.
- Retain the applicable exact-source CI results and adopted numerical checks.
  Require explicit backend selection, the existing CPU/GPU parity checks and
  the existing transfer policy. Record performance as diagnostic evidence.
- Keep MATLAB acceptance in the owning package. Check the MATLAB ABI, native
  MEX execution, shared-kernel parity and applicable Simulink model checks.
  Keep MATLAB and Simulink artifacts outside the Radia Python release gates.
- Give the manager a release manifest with source, artifact identities,
  evidence links, acceptance status and pending items. Public evidence uses
  neutral roles. Keep machine identities and private paths outside it.
- MCP, Cubit, Simulink and Optuna retain independent publication ownership and
  package gates. Use live MCP catalog discovery; do not restore a tracked generated `docs/TOOLS.md` gate.

## Controller coverage

Run `python tools/release_quad.py done --source-sha SHA --manifest FILE
--release-source CHECKOUT` on the reviewed candidate. The controller verifies
the supplied receipts and artifacts without deploying any runtime. Build and
validate each lane using its owning package workflow.

The manifest uses `schema: 1`, the full `source_sha`, an exact-source
`native_ci` record with `passed: true` and `evidence`, and a `lanes` object:

| Lane | Build role | Validation role |
| --- | --- | --- |
| cpu | cpu-build | windows-validation |
| linux-python-gpu | linux-gpu | linux-gpu |
| windows-python-gpu | windows-gpu | windows-gpu |
| matlab | windows-gpu | windows-gpu |

Each lane has a `build` and `validation` record with `role`, `source_sha`,
and `evidence`. The build lists `artifacts`. Validation declares `passed: true`
and an `artifact_sha256` mapping from artifact paths to hashes for exactly
the files it accepted. GPU Python lanes require a `checks.cpu_gpu_parity`
entry. The MATLAB lane requires `checks.simulink_acceptance` and
`checks.ltspice_acceptance`. Every evidence or artifact entry has a manifest-relative
`path` and lowercase SHA-256 `sha256`. Keep the manifest and evidence private.
Use explicit deployment tools only as part of a separately authorized operation.
