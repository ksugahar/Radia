---
name: deploy
description: Current Radia deployment entry point. Use release-quad for releases, verify-deploy for editable source identity, and cubit-plugin-install for the Cubit plugin.
---

# Deploy

Use this skill for deployment triage that does not yet require a new release.
For a release candidate, load `.agents/skills/release-quad/SKILL.md` and follow
that gate instead.

## Machine roles

- LAB and 100号機 are development hosts with editable `radia`,
  `cubit-mesh-export`, and `radia-mcp` installs.
- mdx is the priority CI runner and a MATLAB compute host. CI builds in an
  isolated per-run environment.
- hibino runs long optimization and validation jobs when available.

## Supported operations

1. Confirm editable imports with the `verify-deploy` skill. A running MCP
   server may need code reload or one client reconnect after source changes.
2. Install and verify the independently packaged Cubit backend and embedded
   toolbar with `cubit-plugin-install` and `cubit-smoke-test`.
3. Build native extensions on the machine that will validate them, using the
   selected Python environment and pip `mkl-devel`.
4. Publish only through `release-quad`; its `done` gate is the release decision.

## Forbidden legacy routes

- Do not copy `.pyd`, `.dll`, or `.mex*` files directly between machines.
- Do not restore binary-download scripts or pre-push artifact uploads.
- Do not install PySide or Qt into normal Radia Python. Cubit's private PySide
  runtime remains valid for the Cubit toolbar.
- Do not use retired `.ccl`, Cubit 2025.3, notebook workbench, or standalone
  desktop-panel recipes.
