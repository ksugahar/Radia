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

## Boundaries still requiring a separate migration/validation decision

- Native compatibility command aliases, including `export jmag_nastran`, were
  inspected but not changed. Removing them requires rebuilding and validating
  the matching `.ccm` payload, not just editing Python guidance.
- Existing human GUI settings and Learn credential environment names still
  contain historical Radia naming. They are not Radia runtime dependencies;
  changing them needs explicit settings/launcher migration, not silent deletion
  of user preferences or credentials. No legacy user data was deleted.
- Solver-consumed Kelvin boundary labels and optional Radia integration checks
  remain actual interoperability contracts, not unused compatibility bridges.

Published 1.0.4 and both live editable installations are unchanged by this audit.
The corrected source must be reviewed/integrated and released before describing
these changes as deployed. No compatibility bridge or separate support package
is added by this change.
