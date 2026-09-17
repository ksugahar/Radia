# `radia_ngsolve` ownership and growth boundaries

Status: adopted for incremental migration; this is a source-layout design, not a
new MCP server or a change to numerical acceptance criteria.

## Problem

`radia_ngsolve` combines numerical workflows, MCP tool registration, long-form
guidance, and historical validation cases. Splitting files by size alone would
leave that coupling intact. In particular, `ngsolve_usage("all")` concatenated
every alias value, repeating shared documents. The default response was about
1.9 million characters even though the distinct document bodies total roughly
0.3 million. Token cost and ambiguous ownership matter more than source lines.

## Ownership model

| Plane | Owns | Must not own |
| --- | --- | --- |
| Domain | NGSolve/Radia method, typed inputs, numerical result and artifact contract | MCP registration or model-facing prose |
| MCP | Public tool names, argument schemas, routing, bounded discovery and guidance | A second implementation of a numerical method |
| Validation | Independent reference routes, scientific acceptance tests, case evidence | Runtime behavior merely because a historical test needs it |

Keep the `radia-ngsolve` entry point. A source package is not a server: moving a
method into a cohesive module does not require another process or another public
tool name. Existing `force` and `radia-acoustic` servers own their respective
public discovery surfaces; shared domain implementations should have one owner
and be called rather than copied. Cross-server relocation requires a separate
public-contract and dependency review.

### Knowledge contract

- Default discovery returns a compact index of unique documents. Aliases point
  to one canonical document; they are not duplicate documents.
- A named topic returns its existing body. The explicit `all` route remains for
  compatibility during migration but is too large for ordinary MCP context.
- Future migration should give each document a stable topic ID, move long text
  to packaged resources, and verify those resources from an installed wheel.
  Do not add `.md` resources without package-data and wheel tests.
- Before changing or removing `all`, inventory consumers, publish a migration
  path, and test the versioned public contract. Avoid silently truncating text.

### Gate contract

`slot_gates.py` is legacy aggregation, not the destination for new work. Its
runtime imports are currently concentrated in `em_force_target.py`,
`loop_autolearn.py`, and `regularized_trace_inverse_gate.py`. This observation
does **not** prove that every other function is unused: tests and external Python
users may import functions directly.

For each gate family, record one owner, input/output contract, evidence role
(`runtime acceptance`, `validation-only`, or `historical`), public exposure,
and lifecycle (`active`, `deprecated`, or `retired`). Move a cohesive family
only after direct imports, tool registration, docs, and installed-wheel behavior
are checked. Preserve shims for supported direct imports until a declared
deprecation boundary. Keep scientific thresholds and independent comparison
routes unchanged during structural moves.

## Migration order and acceptance

1. Inventory public tool names and schemas, runtime imports, response sizes,
   package contents, and the tests that exercise each gate family.
2. Make default documentation discovery bounded while preserving named-topic
   lookups and explicit legacy `all`. Test alias resolution and payload size.
3. Assign ownership and evidence roles to gate families; move one connected
   family at a time. Require numerical parity and installed-wheel import tests.
4. Separate server registration into capability-owned modules only after the
   public inventory is stable. Keep the CLI entry point and tool names stable.
5. Replace temporary no-new-definitions guards after `slot_gates.py` has shrunk
   and an owner/lifecycle check covers all exposed gates.

Track default MCP response bytes, alias duplication, public-tool ownership,
cross-plane imports, wheel contents, and numerical parity. Report source lines
as a diagnostic, not a pass/fail quota: a smaller file is not proof of a better
boundary.
