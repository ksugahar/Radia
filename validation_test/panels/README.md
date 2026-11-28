# Application Interface Validation Lane

Current Radia human analysis interfaces are Simulink blocks backed by
`DesignSpec` and headless `calc_*.py` contracts. Notebook workbenches are
retired for every application, including IH.

The primary health gate is:

```powershell
python -m pytest tests/test_application_interface_manifest.py tests/test_simulink_application.py -q
```

This locks the application manifest's Simulink-only states and the common
block runner. The library is locked by
`tests/matlab/test_simulink_workflow.m`.

The old desktop `radia_*.py` PySide6 analysis panels and their ignored GUI
tests have been removed. Interface behavior belongs to Simulink tests; solver
evidence remains here as result-bearing validation.

Cubit has no Radia menu or toolbar; its export runs headless in the external
worker.
