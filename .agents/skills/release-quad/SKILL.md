---
name: release-quad
description: Coordinate the four Radia release cases, preserving exact-artifact and numerical acceptance gates.
---

# Radia release wrapper

Coordinate all four release cases through the dedicated, privately configured skills:

- Build-role editable installation.
- Validation-role editable installation.
- CPU compute build.
- GPU compute build for the qualified device architecture.

`release-dual` and the interim `release-flow` are retired; do not create aliases
or route new work through them. Host-specific skills and host/runtime mappings
are private configuration. Public evidence uses neutral roles.

Normal Radia operation uses editable sources on both Windows roles. Keep
isolated exact-wheel acceptance separate on both roles. Editable installation
does not waive wheel RECORD, native artifact provenance, dependency ABI,
CPU/GPU numerical acceptance or exact-source CI gates. Build cases alone do
not certify GPU runtime, numerical accuracy or performance.

Use the reviewed `tools/release_quad.py` controller. Resolve explicit private
host, interpreter, source and commit bindings before mutation. Preserve other
packages' release ownership and existing environments. If the controller does
not implement all four cases, report the missing gates and keep wrapper
completion pending. Never equate a two-role deployment result with four-case
completion.

Publication requires exact-release native CI, accepted artifact hashes and all
mandatory adopted numerical cases. Preserve user-owned acceptance criteria and
experimental/unqualified labels. Require evidence for all four cases before
claiming wrapper completion. MCP, Cubit and Simulink release independently;
their artifacts and unrelated gates are not folded into the Radia wrapper.

MCP discovery uses the live catalog; do not restore a tracked generated `docs/TOOLS.md` gate.
