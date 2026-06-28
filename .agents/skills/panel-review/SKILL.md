---
name: panel-review
description: Retired PySide6 panel review checklist. Current panel review should target Jupyter notebook workbenches, DesignSpec-to-calc wiring, result artifacts, and no-PySide regression via ipynb-gui-health.
---

# panel-review

This skill is retired for the old PySide6 Layer-3 desktop panels.

For current panel review, inspect:

- `src/radia/*_design.py` for the canonical `DesignSpec`
- `src/radia/*_notebook.py` for `CommandWorkbench` wiring
- `src/radia/panels/notebooks/radia_*.ipynb` for result-bearing notebook cells
- `src/radia/panels/calc_*.py` for the validated headless solver surface
- adjacent `result.json` / `radia_result.v2` artifacts for durable run records

Run:

```powershell
python -m pytest validation_test/panels/test_notebook_workbench.py -q
```

Normal Radia Python should not depend on PySide6.  Do not revive the old GUI
extra as a production fix; use `ipynb-gui-health` and the notebook workbench
contract instead.
