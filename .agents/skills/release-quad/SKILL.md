---
name: release-quad
description: Verify exact CPU solver artifacts on the build and validation roles.
---

# Solver release

The build target (`100`) owns preflight, CI, release gates and editable development.
The NGSolve compute validation target (`validation`) consumes verified radia and
radia-mcp wheels. Cubit is forbidden there. Evidence storage is independent of
deployment routing. Resolve SSH aliases, machine identities and paths privately.

Use `tools/release_quad.py preflight`, `phase8 --target validation,100`, `phase9`
and `done`. `all` runs deployment and consistency checks. Set
`RADIA_RELEASE_EDITABLE_REPO_100` to the clean exact-release source. A repaired
controller may use `done --release-source <exact-tag-checkout>`.

Before publication require exact-SHA native CI, matching wheel hashes, unchanged
numerical gates and machine-readable evidence. Configure private
`RADIA_ACCEPTANCE_HOSTNAME_VALIDATION` and `RADIA_ACCEPTANCE_HOSTNAME_100` for
raw acceptance export. Public evidence uses neutral roles and opaque machine
tokens; raw metadata and batch keys stay private.

MCP releases independently through its package checks and runtime policy.
Use live MCP catalog discovery; do not restore a tracked generated `docs/TOOLS.md` gate.
Cubit releases independently on the build target only. Simulink and Optuna checks
run on the build target. GPU artifacts have independent release gates.
Never stash, reset or clean shared work, interrupt active jobs or silently
substitute a backend. Completion requires `done` to exit zero.
