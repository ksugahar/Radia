---
name: panel-qt-test
description: Retired PySide6 panel test skill. For current Radia panels use ipynb-gui-health and validation_test/panels/test_notebook_workbench.py; keep PySide6 out of normal Radia Python. Use this only as a pointer for intentional legacy adapter archaeology.
---

# panel-qt-test

This skill is retired.  The old headless PySide6 widget tests were for
`src/radia/radia_*.py` desktop panels.  The canonical Radia panel surface is now
the Jupyter notebook workbench.

Run the current panel gate instead:

```powershell
python -m pytest validation_test/panels/test_notebook_workbench.py -q
```

If the Cubit export toolbar is involved, verify the Cubit plugin separately:

```powershell
cubit-plugin-install --verify-only
cubit-smoke-test
```

Do not install PySide6 into normal Radia Python for production.  Cubit's
embedded PySide6 remains protected because it belongs to Coreform Cubit.
