# Application Interface Validation Lane

Current Radia human analysis interfaces are Simulink blocks backed by
`DesignSpec` and headless `calc_*.py` contracts. IH temporarily also keeps
`radia_ih.ipynb` and `ih_notebook.py` for the comparison trial.

The primary health gate is:

```powershell
python -m pytest validation_test/panels/test_notebook_workbench.py -q
```

This locks the IH `CommandWorkbench` artifact/cancel/timeout behavior and the
application manifest's Simulink-first states. The common block runner is locked
by `tests/test_simulink_application.py`; the library is locked by
`tests/matlab/test_simulink_workflow.m`.

The old desktop `radia_*.py` PySide6 analysis panels were removed. Tests that
instantiate those windows are retired and ignored from collection in
`conftest.py`, even on machines where PySide6 is installed. Keep them only as
historical notes until their lessons are distilled into current tests or deleted.

Cubit is the exception: `src/radia/panels/radia_export_menu.py` and
`register_toolbar.py` run inside Coreform Cubit's embedded PySide6 runtime.
Their local test (`test_radia_export_menu.py`) is skipped when normal Radia
Python has no PySide6.
