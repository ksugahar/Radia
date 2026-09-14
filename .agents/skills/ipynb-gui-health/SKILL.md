---
name: ipynb-gui-health
description: Compatibility note for auditing result-bearing Radia docs notebooks and their saved parameterized WebGUI scenes. Notebook workbenches, including IH, are retired; use simulink-app-health for production applications.
---

# Docs Notebook Capability and WebGUI Health

No Radia application has a notebook production workbench. `docs/**/*.ipynb`
is the public capability-showcase, reproduction, and field-inspection layer
only. The `radia-mcp` MCP servers are the primary AI-facing entrypoint and
canonical manual for current operating instructions. Masked Simulink blocks
are the formal human UI but are not documented through competing notebook
workbenches. A separate standalone MATLAB edition remains undecided; `.mlx`
and `.m` files remain implementation assets under `matlab/`, not `docs/`
content. Public notebooks in `docs/` are executed Python capability showcases.
The Python shown there is an LLM-driven implementation behind Radia MCP, not a
separate direct-user product. Simulink operation requires MathWorks' official
MATLAB MCP Server and fails fast when that connection is unavailable.

## Gates

```powershell
python -m pytest packages/radia-mcp/tests/test_document_meta_notebook_audit.py -q
```

Verify that:

- the opening states the engineering problem, the Radia capability/route, the
  expected result or artifact, and the owning MCP capability pack;
- method/application notebooks include governing equations, sufficient
  derivation, and relevant literature citations declared in
  `metadata.radia.citation_keys` and resolved against canonical
  `references.bib`; a citation-free operational notebook records a specific
  `metadata.radia.citation_audit_exempt_reason`;
- the notebook is executed and has no saved error output;
- public examples save WebGUI rich output;
- a primary field uses `Draw(field, mesh, name=..., ...)` with explicit display
  arguments and sets `metadata.radia.webgui_field_required=true`;
- no notebook contains application control widgets or a workbench adapter;
- migration ledgers, cleanup routing, release evidence, and implementation
  handoffs are moved outside `docs/` or rewritten as genuine capability
  showcases;
- procedural instructions are linked to the owning MCP status/usage/recipe
  tool instead of being maintained as a competing notebook manual.

Saved notebook outputs demonstrate the capability claim. Machine-readable
validation thresholds and benchmark evidence belong under `validation_test/`;
an adjacent JSON sidecar is not required merely because a notebook is public.

Use `simulink-app-health` for all five production blocks.
