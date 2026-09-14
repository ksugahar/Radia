---
name: ipynb-gui-health
description: Compatibility note for auditing result-bearing Radia docs notebooks and their saved parameterized WebGUI scenes. Notebook workbenches, including IH, are retired; use simulink-app-health for production applications.
---

# Docs Notebook Capability and WebGUI Health

No Radia application has a notebook production workbench. `docs/**/*.ipynb`
is the public capability-showcase, reproduction, and field-inspection layer
only. The `radia-mcp` MCP servers are the primary manual for current operating
instructions and the only currently supported user entrypoint. MATLAB/MEX and
Simulink remain integration/parity surfaces with an undecided product role;
`.mlx` files are never the sole public showcase because GitHub does not render
them as readable pages.

## Gates

```powershell
python -m pytest packages/radia-mcp/tests/test_document_meta_notebook_audit.py -q
```

Verify that:

- the opening states the engineering problem, the Radia capability/route, the
  expected result or artifact, and the owning MCP capability pack;
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
