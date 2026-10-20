---
name: release-cubit-mesh-export
description: Release and verify the CAD exporter on its sole licensed build target.
---

# CAD exporter release

`tools/release_cubit_dual.py` deploys cubit-mesh-export only on the build target.
Solver releases must not install, repoint or version-gate it. The compute
validation target has no Cubit and is forbidden. Resolve machines and paths privately.

| Role | CLI target key | Deployment tier |
|---|---|---|
| Build | `100` | editable |

Build reviewed source with the current Cubit build script. Verify native
provenance and VC runtime imports against the supported Cubit runtime, curver
import and curving, all six export formats, check-vol and unchanged-case parity.
Run the owning package tests and `tests/test_release_cubit_dual.py`. Missing
required native payloads fail acceptance.

After verified publication, download and hash-check the exact wheel. Create a
clean exact-tag checkout whose packaged bytes equal that wheel. Run:

```powershell
python tools/release_cubit_dual.py --action deploy --run-from 100 `
  --wheel $Wheel --source-sha $SourceSha --source-root-100 $ReleaseRoot `
  --evidence-100 $EvidenceRoot
```

Repeat with `--action done`. Preflight precedes installation. Require fresh imports,
registration, export handlers and MCP checks. Use the shared
[cubit-toolbar-acceptance](../cubit-toolbar-acceptance/SKILL.md) contract on the
sole target; mouse dispatch is optional. MCP stays batch/nographics.

Preserve active work. Force-closing Cubit requires explicit authorization for that
operation. Recover and hash-check evidence before scoped scratch cleanup; retained
unfinished work needs an owner, reason and cleanup trigger. Completion requires
release-owner review, exact artifact evidence and cleanup accounting.
