---
name: cubit-toolbar-acceptance
description: Accept the radia-mcp Cubit Vol Export toolbar/menu using APREPRO/Python execution, the same GUI-handler dispatch path and validated .vol artifacts on the sole build target. Mouse-driven end-to-end GUI tests are optional diagnostics.
---

# Cubit toolbar acceptance

This repository file is the canonical contract for Codex and Claude. Local skill entrypoints reference it rather than maintain separate acceptance rules.

Scope: the **Vol Export** toolbar and menu of `radia_mcp.cubit.gui` (one action, `Export .vol...`), installed by `python -m radia_mcp.cubit.install` into Coreform Cubit 2026.8 on the sole build target. LLM/MCP Cubit execution stays batch/nographics and never runs this acceptance.

## Acceptance policy (human decision, 2026-10-16)

On **the sole build target**, successful APREPRO/Python backend execution plus verification that the GUI action invokes that same tested path is sufficient. Mouse clicks add no mandatory acceptance gate when they invoke the same handler with equivalent options. Do not require mouse-driven exports, Cancel testing, repeated dialog opening, screenshots or a dedicated GUI account merely to close an acceptance.

Historical blocked GUI attempts remain historical evidence; do not rewrite them as successful clicks. A capture/RDP/mouse-automation limitation alone does not block acceptance. Actual registration, dispatch, native identity or export failures still block their affected checks. Do not call backend-only success GUI dispatch evidence.

Human scope update (2026-10-30): the build target is the only Cubit qualification target. Do not require a retired host, another machine/account, or mandatory mouse automation.

## Required evidence

1. Runtime identity: radia-mcp source SHA (or exact wheel SHA-256), interpreter and import paths, Cubit version, and the loaded native plugin path/hash. Disk-file equality alone does not prove the loaded path.
2. Registration: `python -m radia_mcp.cubit.install --verify-only` checks that the deployed plugin matches the built one (size and SHA-256) and that the toolbar is registered. `python -m radia_mcp.cubit.toolbar_smoke` cold-starts the real Cubit GUI and records what Cubit's embedded PySide6 actually shows; run it when registration changed or its evidence is missing.
3. Handler dispatch: `packages/radia-mcp/tests/cubit_mcp/test_cubit_export_menu.py` exercises the production menu helpers (`export_command`, `kelvin_options`, `verification_command`). A headless test may replace dialog inputs through a test seam, but it must execute the production handler and command builder; replacing the handler with a fake or inspecting its source alone is insufficient. Record whether this is programmatic dispatch or an observed GUI invocation; never claim mouse operation that did not occur.
4. Backend: execute the resulting `export vol "<path>" order <n> overwrite` command (with the Kelvin step when selected) through licensed Cubit APREPRO/Python on the sole build target, using independently authored small geometry and explicit settings. Verify the fresh `.vol` with `python -I -B -m radia.vol_check` and its versioned label contract. Reuse already passing evidence for the same bytes and contracts rather than repeat work.

## Allowed test seams and explicit limits

A harness may replace dialog input/exec results, the parent-window lookup, settings persistence, OS-associated-file opening and result presentation. It must retain the installed toolbar payload, the production menu handler and command builder, licensed `cubit.cmd` execution, the native exporter and the real `radia.vol_check` subprocess. Suppressing a warning or backend failure to obtain PASS is prohibited. Record every substituted seam.

Handler equivalence does not certify mouse dispatch, visible dialog rendering, Cancel or dispatch after restart. Such optional observations stay null when untested.

## Verdict and durable record

For the build target record runtime identity, registration result, input source/contract hashes, units/options and the action record:

- action label and registered script/handler;
- `dispatch_pass`, `backend_pass`, `artifact_validation_pass` (true/false/null), dispatch route and evidence paths;
- generated and executed command/options, output paths/sizes/SHA-256 and diagnostics;
- optional `mouse_dispatch_pass`, `cancel_pass`, `restart_dispatch_pass` (true/false/null; untested values must be null).

PASS requires runtime identity, correct registration and the dispatch/backend/artifact checks. Never infer unmeasured checks from green CI or a menu screenshot. Numerical validation gates remain unchanged.

## Evidence and cleanup

Recover result JSON, artifacts, commands and hashes to private durable evidence, and preserve failed diagnostics. Delete only this run's enumerated task-owned scratch after hash-verified recovery; preserve deployed runtimes, evidence, other sessions' work and active jobs. Close only owned sessions normally; force-closing Cubit requires explicit authorization for that operation. Report retained paths with owner, reason and cleanup trigger.

## Optional real-desktop diagnostics

Use mouse/UI testing only for a reported GUI-specific problem or explicitly requested additional coverage. Read the computer-use skill, preserve active human sessions, use an authorized task-owned session and confirm a nonblack capture and detectable dialog before further input. If these fail, report the diagnostic limitation without substituting it for product failure. Never turn such diagnostics back into a universal acceptance gate.
