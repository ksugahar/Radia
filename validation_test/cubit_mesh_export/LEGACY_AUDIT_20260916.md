# Cubit independent-runtime legacy audit

Scope: independent distribution ownership, executable MCP operating guidance,
recipe artifacts, and the old GUI transport. This is a source audit, not a new
release acceptance or proof that every numerical/native route is defect-free.

## Corrected in the first audit change

- The operating manual required Radia installation, attributed the toolbar to
  Radia, and prescribed Cubit deployment to compute hosts. It now names the
  independent package and LAB/100 release boundary; combined Radia acceptance
  remains explicit and optional.
- Recipe curation emitted JSON directly as Python source. Boolean/null values
  made the generated module fail at import. It now loads an escaped JSON literal;
  an executable round-trip regression covers true/false/null values.
- Generated recipe text asserted contributor consent merely from using Radia
  MCP and named the wrong distribution for publication. It now requires explicit
  permission/confidentiality review and identifies cubit-mesh-export.
- Validation: 759 MCP tests passed, including two new audit regressions.

## Second pass: retired implementations and contradictory contracts

- Deleted the GUI file-drop transport and Qt bootstrap runner, shared PID
  attachment and foreign-PID termination. Both the client constructor and
  daemon reject non-batch execution. Human toolbar/rendering release probes
  remain separate and unchanged.
- Removed the always-failing `cubit_snapshot` tool, underlying hardcopy RPC,
  and ignored `include_failed` journal argument. Native recording and APREPRO
  preservation remain the provenance contract.
- Replaced old GUI logs/session-mode status with current owned-process stderr.
  Removed the lockstep 4.7.0 release synopsis and shared-session topic aliases;
  corrected warmup and independently owned cross-server support guidance.
- Removed three no-op lint compatibility functions and their stub-only tests;
  the real deleted-API detection and lint golden still run.
- Fixed latent lifecycle faults uncovered during retirement: unbounded startup,
  competing stderr reads, automatic replay after an unknown command outcome,
  silent continuation after geometry loss, and concurrent duplicate startup.
  Startup is bounded, calls serialized, failed children reaped, and Windows
  children tied to their owner by a kill-on-close Job Object.
- Fixed the constructor ignoring `set_cubit_bin_dir` and discovery ignoring the
  documented `CUBIT_INSTALL_DIR`; moved native journals to Cubit-owned scratch.

### Validation

- MCP tests: **763 passed** (759 baseline, 15 retired checks, 19 new lifecycle
  regressions). No skips in this run.
- Installer/menu/standalone-GUI-assets/toolbar-smoke/release-dual tests:
  **78 passed**. These are regression tests, not a fresh two-host GUI acceptance.
- Real LAB Cubit: headless ready/ping/brick command/native journal/shutdown pass;
  only the newly owned child was stopped. No interactive GUI was started.
- Fresh-build local wheel: native provenance gate passed; installed into an
  isolated test environment without Radia/radia-mcp/cae-mcp-core. MCP stdio
  discovery/status and selftest passed, **54 tools**, no bootstrap module or
  snapshot tool. This local test wheel retains version 1.0.4 and is NOT a new
  published artifact; do not distribute it as the released 1.0.4 wheel.

## Third pass: legacy removal

- Removed the old native command alias and its command class. The rebuilt plugin
  exposes only `export nastran_bdf`; binary inspection confirms the old class,
  command, plugin identity, and product header are absent.
- Renamed the native plugin, message filter, C++ namespace, GUI module, toolbar,
  icon, logs, eager-start environment variable, and exporter-owned settings root
  to Cubit Mesh Export names. Deleted the obsolete hard-coded batch builder.
- The installer removes exact exporter-owned legacy startup, toolbar, archive,
  settings, and native-plugin paths. These names are one-way deletion targets,
  not compatibility relays; unrelated files and parent directories are retained.
- Removed the obsolete MCP Netgen-export compatibility gate. Current guidance
  and active validation journals use only supported native commands.
- Removed the Radia-dependent native coil command, obsolete GUI-transport E2E,
  legacy test journals, and exporter-owned `RADIA_*` environment variables.
  Radia's coil-generation API remains on the Radia side of the dependency.
- Added a clean-package build command and wheel-content CI gate. This prevents
  deleted Python, toolbar, icon, or bootstrap files from leaking out of a stale
  incremental `build/lib` directory.
- Solver-consumed Kelvin boundary labels and optional Radia integration checks
  remain because they are current interoperability contracts, not legacy paths.
- Removed programmatic license-cache warmup. `--setup` is now a non-destructive
  doctor; license activation stays in Coreform's official UI.
- Rejected unsafe race-history identifiers before file lookup and moved every
  synchronous MCP tool invocation to a worker thread so a long Sculpt or Cubit
  call does not freeze the MCP event loop.
- Made missing material domains fail mesh-quality acceptance, made malformed
  managed startup markers fatal, and removed transient staging paths from the
  official toolbar archive.
- Deleted the synthetic v29-v56 mixed-transition identity ladder: 28 generated
  test generations, 12 recursive identity modules, and 9,400+ lines of repeated
  predicates had no package-side evidence producer. The production gate retains
  shared-face ownership, two-sided manifold, family/quality inventory, Gmsh
  connectivity, independent volume closure, and headless-process classification.

### Third-pass validation

- Rebuilt the `.ccm` and native curver; provenance manifest updated.
- Current package source suite: **488 passed** after deleting the synthetic
  identity generations. Security/runtime, mesh-quality and installer review
  focus: **59 passed**; GUI independence, native registry, Kelvin and
  release-dual focus: **84 passed**. No skips in these runs.
- Built-wheel standalone checks passed without Radia, radia-mcp, or
  cae-mcp-core: **54 MCP tools**, **469 installed-wheel MCP tests**, current
  toolbar paths, and no retired bootstrap, snapshot tool, GUI module, icon,
  toolbar template, license warmup, or versioned identity modules. The wheel
  contains **83 members**; Sculpt/vfrac contracts run in this clean environment
  with the declared extras and no skips.
- A real headless Cubit run loaded the rebuilt plugin and exported a canonical
  BDF whose header identifies cubit-mesh-export. The installed live plugin was
  not changed during this source audit.

Published 1.0.4 and both live editable installations are unchanged by this audit.
The corrected source must be reviewed/integrated and released before describing
these changes as deployed. No compatibility bridge or separate support package
is added by this change.
