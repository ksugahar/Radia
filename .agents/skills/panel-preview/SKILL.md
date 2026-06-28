---
name: panel-preview
description: Retired PySide6 desktop panel screenshot skill. Current Radia panels are Jupyter notebook workbenches; verify notebooks with ipynb-gui-health and saved result-bearing notebooks, not desktop Qt screenshots.
---

# panel-preview

This skill is retired for the old PySide6 desktop panel screenshots.

The current previewable artifact is the executed notebook itself:

```text
src/radia/panels/notebooks/radia_<app>.ipynb
```

For health checks, run:

```powershell
python -m pytest validation_test/panels/test_notebook_workbench.py -q
```

For human-facing visualization inside notebooks, use `netgen.webgui`.  For
headless/LLM validation, use durable `.msh v4.1`, JSON, and saved notebook
outputs.  Do not install PySide6 into normal Radia Python for production.
