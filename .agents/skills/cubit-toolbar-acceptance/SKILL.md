---
name: cubit-toolbar-acceptance
description: Accept cubit-mesh-export toolbar releases using APREPRO/Python execution, the same GUI-handler dispatch path, native identity and validated artifacts on the sole build target, with distinct exact-wheel and editable evidence. Mouse-driven end-to-end GUI tests are optional diagnostics.
---

# Cubit toolbar acceptance

This repository file is the canonical contract for Codex and Claude. Local skill entrypoints reference it rather than maintain separate acceptance rules.

## Acceptance policy (human decision, 2026-10-16)

On **the sole build target**, keeping exact-wheel qualification and normal editable identity/registration as explicit records, successful APREPRO/Python backend execution plus verification that the GUI actions invoke that same tested path is sufficient. Mouse clicks add no mandatory acceptance gate when they invoke the same handler with equivalent options. Do not require six mouse-driven exports, Cancel testing, repeated dialog opening, screenshots, restart dispatch of every button, or a dedicated GUI account merely to close a release.

This supersedes this skill's previous mandatory real-desktop six-button policy. Historical blocked GUI attempts remain historical evidence; do not rewrite them as successful clicks. A capture/RDP/mouse-automation limitation alone does not block acceptance under this policy. Actual registration, dispatch, native identity or export failures still block their affected checks. Do not call backend-only success GUI dispatch evidence.

Codex task owners implement and self-validate, and submit exact-source evidence to Codex repository management. Codex repository management owns GitHub execution, release and deployment; Claude provides reciprocal independent review. Creating or changing this skill does not itself complete the tests.

Human scope update (2026-10-30): the build target is the only Cubit deployment/qualification target. Do not require a retired host, another machine/account, or mandatory mouse automation. Prior two-host attempts remain historical evidence. The versioned v2 test seams and six production handler/backend/artifact gates below remain required.

## Required evidence

1. Establish the selected release version, source/tag SHA, published wheel SHA-256, deployment receipt, interpreter/import paths, Cubit version and selected runtime identity. Verify the actual loaded plugin path/hash against native_payloads.json (.ccm may be deployed as plugins/cubit_mesh_export.dll), and the matching pyd hash. Disk-file equality alone does not prove the loaded path. Preserve the declared editable / wheel identity for each runtime record.
2. Verify all six installed toolbar script registrations resolve to the intended package's GUI handlers. Exercise the same handler/command-building path programmatically with explicit options and record the resulting APREPRO/Python invocation. A headless test may control dialog inputs through a test seam, but it must execute the production handler and command builder; replacing the handler with a fake or inspecting its source alone is insufficient. Record whether this is programmatic GUI-handler dispatch or an observed GUI invocation; never claim mouse operation that did not occur. Existing automated dispatch tests may be reused when they match the release and runtime contract.
3. Execute the resulting backend commands through licensed Cubit APREPRO/Python on the sole build target using independently authored small geometry and explicit settings. Match dispatch command/options to the real execution and validate fresh output artifacts. Reuse already passing release evidence for the same bytes and contracts rather than repeat work. Missing format coverage stays unverified.
4. Confirm installed startup/menu/toolbar registration using deployment verification and existing startup/display evidence. A separate restart/visibility probe is appropriate when registration changed or evidence is missing; it is not an automatic demand for six mouse clicks. Preserve supported batch/nographics-only behavior for MCP execution.

| Action | Backend/output contract |
|---|---|
| Netgen Vol (.vol) | `.vol`, checked with the versioned label contract and check-vol |
| GMSH (.msh) | `.msh`, supported format-reader/header validation |
| Nastran (.bdf) | `.bdf`, supported reader/header validation |
| VTK (.vtk) | `.vtk`, supported reader/header validation |
| FEMEEM | directory with in.dat, sin.dat.B, sina.dat, d3; sizes/hashes and contents checked |
| MEG (ELF/MAGIC) | `.meg`, required block labels/options explicitly recorded and validated |

Read the selected release's package-owned cubit_gui/cubit_export_menu.py, toolbar scripts/template and C++ writers when establishing the dispatch/output contract. FEMEEM's quoted argument denotes a directory (ExportFemeemCommand.cpp::write_femeem), despite its filename parser label. Use first-order tetrahedra for FEMEEM/MEG and record scale/block mapping. No external commercial application is required.

The bbnd_brick.py fixture defaults to order 2 and includes the curve label test_edge. Select BBND_ORDER=1 for the baseline or use a recorded elementary first-order brick with its own material/surface contract. Netgen order-2/BBND validation is separate affected-feature coverage; other formats need not reproduce codimension-2 labels. Never reset another session's model or use private commercial reference data.

## Allowed test seams and explicit limits

The v2 harness may replace dialog input/exec results, the parent-window lookup, settings persistence, OS-associated-file opening and result presentation. It must retain the installed toolbar scripts, production launch_export/format handler, production ExportDialog constructor and cubit_command builder, licensed cubit.cmd execution, native exporter and real Netgen verification subprocess. Suppressing a warning or backend failure to obtain PASS is prohibited. Record every substituted seam.

This v2 acceptance does not certify mouse dispatch, visible dialog rendering, Cancel, dispatch after restart, automatic discovery during normal GUI startup, or the BBND GUI route. Such optional observations stay null when untested, even if not required by policy. Existing startup/display or BBND backend evidence may be cited separately for its own scope; handler equivalence must not be described as measuring those GUI behaviors.

## Verdict and durable record

Use schema `cubit-mesh-export.toolbar-acceptance.v2` and `acceptance_policy: aprepro-python-gui-handler-equivalence-2026-10-16`. Keep v1 records unchanged. For each host record release/runtime identity, registration result, input source/contract hashes, units/options and six action records containing:

- action label and registered script/handler;
- `dispatch_pass`, `backend_pass`, `artifact_validation_pass` (true/false/null), dispatch route and evidence paths;
- generated and executed command/options, output paths/sizes/SHA-256 and diagnostics;
- optional `mouse_dispatch_pass`, `cancel_pass`, `restart_dispatch_pass` (true/false/null; untested values must be null, with mouse_scope or outside_acceptance_scope recording the policy scope separately).

Host PASS requires release/native identity, correct registration and all six required dispatch/backend/artifact checks. Overall PASS requires exact-wheel qualification and normal editable identity/registration on the sole build target. Same-byte/options/native execution evidence may be reused explicitly; it does not replace the distinct import/source and registration records. Example: `editable identity PASS; wheel qualification null; overall incomplete` reflects missing required execution, while `mouse_capture blocked; required checks PASS on the sole build target; overall PASS under v2` is valid. Never infer unmeasured checks from green CI or a menu screenshot. Numerical validation gates remain unchanged.

Recover result JSON, artifacts, commands and hashes to private durable release evidence. Preserve failed diagnostics. Remove only verified task-owned scratch within authorized scope, preserving deployed releases, evidence and other sessions' WIP. Report remaining session/process cleanup separately; acceptance-policy changes do not authorize force-killing a process.

## Standard post-release cleanup

Acceptance and cleanup are separate recorded gates for every release. Before starting validation, assign one cleanup executor through management and enumerate only this run's absolute scratch paths and owned processes. The selected human authorization applies to that scope; creating this skill does not authorize unrelated deletion or forced shutdown.

After tests, recover inputs, scripts, commands, artifacts, successful and failed logs, runtime/source identity and result JSON to durable private release storage. Build a SHA-256 inventory and verify every destination against the original before removing the original. Preserve the selected deployed checkout/wheel and canonical evidence.

The assigned executor checks exact resolved paths, active/queued consumers, unique untracked/ignored work and symlinks/junctions/reparse points. Close only owned sessions normally. Delete only approved task-owned paths using literal path operations, then verify their absence. Preserve another session's work and active jobs. Retire obsolete Git branches/worktrees only after main inclusion (including equivalent rewritten commits), consumer and unique-work checks, through the appropriate Git/app-managed workflow.

Record a cleanup receipt with executor, timestamp/timezone, exact targets, recovery manifest hash, prechecks, deleted/absent outcomes and retained owner/reason/trigger. If approval review rejects deletion, record the rejection and zero/partial execution accurately; do not retry through another shell, API or agent to bypass it. Management resolves the applicable restriction before further action. Changing executor is not itself resolution.

The management closure record references the separate wheel/editable acceptance records, Claude's actual independent review and the cleanup receipt, and verifies the final inventory after any review edits. A queued packet is not a confirmation. A cleanup blocker remains a separate open item even if numerical/handler acceptance passes; name retained paths and their owner instead of claiming full completion.

## Optional real-desktop diagnostics

Use mouse/UI testing only for a reported GUI-specific problem or explicitly requested additional coverage. Read the computer-use skill, preserve active human sessions, use an authorized task-owned session and confirm a nonblack capture and detectable dialog before further input. If these fail, report the diagnostic limitation without substituting it for product failure. Account creation/licensing/session switching require their own authorization. Never turn such diagnostics back into a universal release gate.

## Handoff

Return the sole-target wheel/editable verdicts, required six-action coverage, exact evidence/hashes, unresolved findings and cleanup status to Codex management. Ask management to reconcile the selected release evidence against v2 and close only genuinely satisfied remaining gates; do not automatically label an old blocked test PASS.
