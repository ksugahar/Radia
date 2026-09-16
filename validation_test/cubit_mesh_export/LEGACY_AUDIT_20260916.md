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

## Remaining retirement work (not complete)

1. `mcp/session.py` still implements a GUI file-drop transport and
   `mcp/bootstrap.py` supplies its Qt polling runner. Public MCP helpers request
   batch explicitly, but direct `CubitSession(mode="gui")` remains possible.
   Retire the transport, not the separate human-facing GUI release probes.
2. `cubit_doctor`, error-log pointers and `cubit_session_status` still describe
   GUI drop files, attached PIDs and `RADIA_CUBIT_SESSION_MODE`. Replace these
   diagnostics with actual process-owned stdio state. Preserve native journal
   provenance and owned-child cleanup; do not use historical PID files to kill
   an unrelated process.
3. Remove file-drop-only tests and replace relevant ownership/lifecycle checks
   with batch regressions. Update scripting knowledge about persistent GUI
   attachment, warmup and old editor restart instructions at the same time.
4. Audit backend compatibility command aliases and retired manual parameters
   separately against current callers. Do not change a native command contract
   without rebuilding and validating the matching native wheel payload.

Published 1.0.4 and both live editable installations are unchanged by this audit.
The corrected source must be reviewed/integrated and released before describing
these changes as deployed. No compatibility bridge or separate support package
is added by this change.
