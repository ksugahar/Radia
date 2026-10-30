---
name: release-quad
description: Build, accept, publish and distribute exact Radia solver wheels to the build and production consumer roles.
---

# Solver release

The build target (`100`) owns preflight, CI, release gates and editable development.
The `validation` role is also a permanent production consumer of the complete
Radia distribution, not a temporary or HDiv-only acceptance target. It consumes
verified wheels built on the build target. Cubit is forbidden there. Resolve
SSH aliases, machine identities and runtime paths privately; evidence storage
is independent of deployment routing.

Use `tools/release_quad.py preflight`, `phase8 --target validation,100`, `phase9`
and `done`. `all` runs deployment and consistency checks. Set
`RADIA_RELEASE_EDITABLE_REPO_100` to the clean exact-release source. A repaired
controller may use `done --release-source <exact-tag-checkout>`.

Before publication require exact-SHA native CI, matching wheel hashes, unchanged
numerical gates and machine-readable evidence. Configure private
`RADIA_ACCEPTANCE_HOSTNAME_VALIDATION` and `RADIA_ACCEPTANCE_HOSTNAME_100` for
raw acceptance export. Public evidence uses neutral roles and opaque machine
tokens; raw metadata and batch keys stay private.

## Production wheel distribution

Keep development and normal distribution runtimes distinct. Development uses an
owned venv; normal distribution may use the explicitly selected standard Python
or a dedicated runtime venv. A venv is not a prerequisite for normal pip use.
Preserve existing environments and editable sources; the assigned executor owns
installation and deployment.

Build the formal artifact on the build target, record its source SHA and wheel
SHA256, and transfer the same hash-verified wheel to the production consumer.
Before publication, install the accepted candidate using that runtime's absolute
Python executable and `-m pip install <verified-wheel-path>`. After publication,
normal distribution uses `-m pip install "radia==<release-version>"`, or
`"radia[gpu]==<release-version>"` when the released GPU feature/runtime is needed.
Resolve/download the published wheel and verify it against the accepted wheel
hash; never rebuild or substitute another artifact during promotion. Unsupported
platform/ABI or missing native payload fails rather than falling back.

Verify the actual Python executable, installed/imported Radia version and source
identity, wheel RECORD/native hashes, dependency ABI and `pip check` in each
selected runtime. Run the required CPU acceptance and the relevant GPU acceptance
on the consumer for each supported released capability. Record hardware, backend,
artifact hashes and machine-readable results. Import success alone does not
certify numerical accuracy, GPU execution, architecture compatibility or speed.
Retain experimental/unqualified labels and all user-approved acceptance criteria.

Use one explicit private SSH/runtime configuration for deploy, wheel verification,
phase9 and done; a PATH alias, WindowsApps launcher or legacy endpoint must not
silently select a different Python. If the current controller cannot target the
selected runtime, fix and review the controller before deployment; do not bypass
it with wrappers or SSH shims. No reinstall is required merely because this skill
was read. Completion includes verified production distribution and `done` exit0.

MCP releases independently through its package checks and runtime policy.
Use live MCP catalog discovery; do not restore a tracked generated `docs/TOOLS.md` gate.
Cubit releases independently on the build target only. Simulink and Optuna checks
run on the build target. GPU artifacts have independent release gates.
Never stash, reset or clean shared work, interrupt active jobs or silently
substitute a backend. Completion requires `done` to exit zero.
