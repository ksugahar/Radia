---
name: deploy
description: Current Radia deployment entry point. Use release-quad for releases, verify-deploy for editable source identity, and cubit-plugin-install for the Cubit plugin.
---

# Deploy

Use this skill for deployment triage that does not yet require a new release.
For a release candidate, use the owning release skill: `release-quad` for the
numerical Radia solver/Simulink package, `release-radia-mcp` for radia-mcp, and
`release-cubit-mesh-export` for the exporter and Cubit MCP.

## Machine roles

- LAB and 100号機 are development hosts with editable `radia`,
  `cubit-mesh-export`, and `radia-mcp` installs.
- mdx1 and mdx2 are the priority CI runners and MATLAB compute hosts. CI builds
  in isolated per-run environments.
- hibino runs long optimization and validation jobs when available.
- No workflow selects a LAB runner. The EqnEdit64 release lane is GitHub-hosted
  end to end: the signed executable reaches it as an asset of the
  eqnedit64-staging release, staged by sync_to_o.ps1 alongside the O: update. A
  runner service runs as NETWORK SERVICE and can read neither O: nor the
  workgroup share, and a LAB runner would also inherit Radia's heavy builds.

## Supported operations

1. Confirm editable imports with the `verify-deploy` skill. A running MCP
   server may need code reload or one client reconnect after source changes.
2. Install and verify the independently packaged Cubit backend and embedded
   toolbar with `cubit-plugin-install` and `cubit-smoke-test`.
3. Build native extensions on the machine that will validate them, using the
   selected Python environment and pip `mkl-devel`.
4. Publish only through the owning release lane. `release-quad done` decides
   the numerical solver/Simulink release; it does not publish or deploy
   radia-mcp or cubit-mesh-export.

## Forbidden legacy routes

- Do not copy `.pyd`, `.dll`, or `.mex*` files directly between machines.
- Do not restore binary-download scripts or pre-push artifact uploads.
- Do not install PySide or Qt into normal Radia Python. Cubit's private PySide
  runtime remains valid for the Cubit toolbar.
- Do not use retired `.ccl`, Cubit 2025.3, notebook workbench, or standalone
  desktop-panel recipes.
