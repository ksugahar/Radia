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

## Follow-up audit at main b3bfeed5a

This is a focused source and contract audit, not a complete-package or ver5
numerical acceptance result.

- SDK floor: replace the unsupported `mcp>=1.0,<2` declaration with the tested
  `mcp>=1.20.0,<2` range. The minimum-SDK CI job derives its exact pin from
  package metadata and tests Python 3.10 separately from the normal resolver.
  Registration with annotations and metadata, schema listing, default/coerced
  arguments, coarse dispatch, and real stdio reload notifications are checked.
  Local Python 3.12 runs pass all 20 selected tests with SDK 1.20.0 and 1.27.0.
  This establishes a supported floor, not a claim to have tested every SDK
  release or every optional application dependency. It is not a wheel-only test.
- Reload controls remain open: `register_reload_tool` is unconditional and
  new callable discovery still uses a shared name prefix plus module ownership.
  Rollback and conservative annotations do not close exposure-policy questions.
- Grouped argument validation remains open: `CoarseToolRegistry.run` calls
  `entry.function(**kwargs)` directly, outside the individual SDK tool schema.
- Full-audit infrastructure already exists through `workflow_dispatch`; the
  missing item is current acceptance evidence, not another duplicate workflow.
- The ordinary tests intentionally add this checkout's source to `sys.path`.
  Passing them is not evidence of installed-wheel resource completeness or
  correct live editable ownership.

### MATLAB / NGSolve: knowledge versus acceptance

`matlab_radia_mex_contract("ngsolve")` and the MATLAB skill expose the native
`radia.ngsolve` mesh, space, form, matrix, vector and field-handle boundary.
The skill distinguishes native local HCurl/CLN reduction from the full
HCurl-VIM/BEM route and from exported reduced-state families. This is substantial
executable guidance, not a full clone of upstream NGSolve or model retraining.

Durable evidence is in `validation_test/ngsolve_matlab_parity/`:

- Core 100 linear cases, with sparse entries, matvec and native solves.
- HIBINO extended run dated 2026-09-01: breadth 500, scale 20 (21,515 to
  134,130 DoFs), and manufactured-solution 15 cases passed, with zero remaining
  native handles. The result identifies MATLAB R2026a, NGSolve 6.2.2606 and
  MEX SHA-256 `9149f254621d30128a661e1d4159e5939459eaa49f82427204fd1ebb78e07beb`.
- The eight fast core/extended artifact contracts pass in this audit. They
  validate the saved evidence; they do not rerun MATLAB or certify a new binary.
- The MATLAB MCP contract suite has 14 passes and one existing failure:
  `test_root_readme_publishes_native_topology_mex_parity` finds a 365-command
  current inventory but `matlab/README.md` still says 364. Do not suppress it
  or equate the stale documentation number with a numerical solver failure.

Before ver5 acceptance, select the affected numerical parity lane against the
exact candidate MEX/source provenance. Production geometry, nonlinear/moving
applications and the complete upstream NGSolve API are not proven by these
bounded linear parity results. No new heavy numerical run was launched here.
