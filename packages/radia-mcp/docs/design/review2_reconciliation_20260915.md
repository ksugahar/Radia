# Second MCP review reconciliation

Review baseline: main `afb8c621c`, radia-mcp 1.4.54. This records the supplied
review summary, not an independent reproduction of every finding in the full
external report. No shared checkout, editable installation or live client is
changed by this maintenance patch.

| Finding | Disposition and evidence |
| --- | --- |
| Non-Windows `RADIA_MCP_TEMP` ignored | Fixed both native journal and headless batch paths through one helper. Mocked Windows/Linux/macOS cases check configured roots and platform defaults. The native journal suite passed twice consecutively (20 cases each); production Cubit was not started. Existing journal overwrite refusal remains in place. |
| README `49 servers` | Replaced current fixed counts with catalog discovery; historical changelog counts are not rewritten. |
| Missing 1.4.54 changes | Added SDK floor, argument binding/validation and restricted reload entries. New scratch-path changes are under Unreleased, not falsely described as published. |
| Blanket aliases and CBCR | Naming guidance explicitly excludes retired GUI tools; CBCR carries a historical-design notice and no longer claims current production status. Supported headless checkpoint operations remain available. |
| NumPy in minimal tests | Clarified runtime versus test dependencies. Mandatory runtime remains MCP; scoped CSV tests install NumPy, while minimum-SDK registration/dispatch does not add it. This does not establish all tools work without their optional dependencies. |
| 29 package tests with monorepo dependencies (review count) | OPEN. MATLAB/policy/validation artifact checks need classification as repository integration tests or conversion to synthetic local fixtures. Static path inspection confirms examples, but no fresh exhaustive count or package-only full-suite acceptance is claimed. Missing coverage must not be closed by adding skips. |
| Old S: source | An old checkout alone does not identify any live client. Use the original client's provenance; do not pull/repoint the shared tree based on its name. Integration and live acceptance remain separate. |

Remaining work is the package-test/integration boundary, including the reported
29 cases. This patch does not claim package-wide completion or a release.
