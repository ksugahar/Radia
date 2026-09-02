---
name: ipynb-gui-health
description: Compatibility note for auditing result-bearing Radia docs notebooks and their saved parameterized WebGUI scenes. Notebook workbenches, including IH, are retired; use simulink-app-health for production applications.
---

# Docs Notebook WebGUI Health

No Radia application has a notebook production workbench. `docs/**/*.ipynb`
is the public example, reproduction, and field-inspection layer only.

## Gates

```powershell
python -m pytest packages/radia-mcp/tests/test_document_meta_notebook_audit.py -q
```

Verify that:

- the notebook is executed and has no saved error output;
- no adjacent JSON is required; benchmark/evidence JSON belongs in
  `validation_test/`;
- public examples save WebGUI rich output;
- a primary field uses `Draw(field, mesh, name=..., ...)` with explicit display
  arguments and sets `metadata.radia.webgui_field_required=true`;
- no notebook contains application control widgets or a workbench adapter.

Use `simulink-app-health` for all five production blocks.
