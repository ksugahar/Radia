# radia-mcp runtime review reconciliation

Scope: review supplied on 2026-09-13, targeting main at PR #224
(`73a56b49c`). This record separates reproduced defects from proposals that
still need verification. It does not certify all package maintenance complete.

## Reproduced and repaired in this change

- Presentation bedrock lint called an undefined `_scan_hedges`. Reuse the
  existing `_shared.hedges.scan_hedges`; do not merge grant, paper and slide
  policy into one undifferentiated checker.
- Three public slot gates lacked `re` for comma/semicolon/space-separated
  coordinates. Add the import and compare their string and numeric input paths.
  Incomplete artifact packets in these regression tests are not solver evidence.
- Two legacy gate tests incorrectly required the shared dispatcher to advertise
  read-only and idempotent behavior. Keep conservative dispatcher annotations;
  its other operations can write artifacts. Existing stdio tests remain enabled.
- Worker `SystemExit` or `KeyboardInterrupt` could leave AsyncRunner RUNNING.
  Catch these at the worker boundary, record FAILED/end_time, and test restart.
  This does not intercept interrupts in the caller's main thread.

## Remaining review lanes (not closed by this change)

1. SDK dependency floor: reproduce registration, metadata and tool-call behavior
   in isolated environments, then synchronize package metadata and a minimum-SDK
   CI lane. Importing FastMCP alone is insufficient compatibility evidence.
2. Hot reload: assess explicit disable/registration controls and newly discovered
   callable exposure. Preserve experimental editable development permissions;
   do not introduce an implicit source freeze or terminate unrelated clients.
3. Coarse dispatcher: preserve public signatures while reusing SDK argument
   validation; test defaults, coercion, invalid inputs and async behavior.
4. Test distribution: distinguish wheel-only contracts from monorepo integration
   tests. Missing required evidence must not become a blanket successful skip.
5. Static findings: verify duplicate dictionary keys, duplicate definitions and
   remaining undefined names individually before deleting or consolidating code.
6. Documentation counts, stdout encoding and large-module refactoring: derive
   inventory counts from the catalog, check actual CLI behavior, and preserve
   domain-specific policy and public contracts during extraction.
7. Full-audit coverage: decide and test an appropriate CI trigger separately;
   do not claim targeted CI is an all-tests or installed-wheel audit.

No shared editable installation, live MCP process, solver binary, release tag,
or other contributor's uncommitted work is changed by this maintenance patch.
