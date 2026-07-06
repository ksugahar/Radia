# Panel Validation Lane

Current Radia analysis panels are Jupyter notebook workbenches:

- `src/radia/panels/notebooks/radia_*.ipynb`
- `src/radia/*_design.py`
- `src/radia/*_notebook.py`
- `src/radia/notebook_workbench.py`
- headless `src/radia/panels/calc_*.py`

The primary health gate is:

```powershell
python -m pytest validation_test/panels/test_notebook_workbench.py -q
```

This locks the `DesignSpec` initial-value contract, the
`CommandWorkbench` run artifact schema, background/cancel/timeout behavior,
manifest states, and the no-PySide rule for notebook-facing panels.

The old desktop `radia_*.py` PySide6 analysis panels were removed. Tests that
instantiate those windows are retired and ignored from collection in
`conftest.py`, even on machines where PySide6 is installed. Keep them only as
historical notes until their lessons are distilled into notebook-workbench
tests or deleted.

Cubit is the exception: `src/radia/panels/radia_export_menu.py` and
`register_toolbar.py` run inside Coreform Cubit's embedded PySide6 runtime.
Their local test (`test_radia_export_menu.py`) is skipped when normal Radia
Python has no PySide6.
