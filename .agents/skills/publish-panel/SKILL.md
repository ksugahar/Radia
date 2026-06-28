---
name: publish-panel
description: Retired PySide panel release checklist. The canonical Radia panel release gate is now the Jupyter notebook workbench path; use ipynb-gui-health for notebook panels and deploy/release-qud for Cubit plugin deployment. Do not install PySide6 into normal Radia Python for production.
---

# publish-panel

This skill is retired for the old Layer-3 PySide6 desktop panel era.

The current panel surface is:

```text
src/radia/panels/notebooks/radia_<app>.ipynb
radia.<app>_design.<App>DesignSpec
radia.<app>_notebook.<App>Workbench
src/radia/panels/calc_<app>.py
```

## Current Release Gate

Use `ipynb-gui-health` after editing notebook panel code:

```powershell
python -m pytest validation_test/panels/test_notebook_workbench.py -q
```

Then verify the Cubit export/plugin surface separately:

```powershell
cubit-plugin-install --verify-only
cubit-smoke-test
```

Production installs use `radia[cubit]`, not the old GUI extra.  Normal Radia
Python should not install PySide6.  Coreform Cubit's bundled PySide6 is
protected and must not be uninstalled or deleted.

Use git history for the old PySide screenshot checklist if you must maintain a
legacy adapter intentionally.
