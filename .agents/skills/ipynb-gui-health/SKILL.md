---
name: ipynb-gui-health
description: Verify the temporary Radia IH notebook workbench used alongside the IH Simulink block. Use after changes to ih_notebook.py, ih_design.py, radia_ih.ipynb, notebook_workbench.py, or the IH dual-interface manifest. This is not the production gate for EM, PCB, Motor, or Stream Function; use simulink-app-health for those applications.
---

# IH Notebook Comparison Health

Only IH retains a notebook workbench, temporarily, so its operating quality
can be compared with the IH Simulink block. EM, PCB, Motor, and Stream Function
must not regain notebook production surfaces.

## Gates

```powershell
python -m pytest validation_test/panels/test_notebook_workbench.py -q
python -m pytest tests/test_simulink_application.py -q
python tools/audit_new_panel_contract.py
```

Verify that:

- `radia_ih.ipynb` delegates to `IHDesignSpec` and `IHWorkbench`.
- The notebook and Simulink block reach the same validated headless CLI.
- Runs retain timeout, cancellation, `run.log`, and structured result artifacts.
- The interface manifest marks IH `active-dual-comparison`.
- No EM, PCB, Motor, or Stream Function notebook/adapter exists.

Use `simulink-app-health` for the five production blocks. Retire this skill
when IH has met its documented Simulink migration gates and its notebook is
removed.
