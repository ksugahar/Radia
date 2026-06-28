"""Retired PySide panel review knowledge, redirected to notebook panels."""

PANEL_REVIEW_NOTEBOOK = """
# Radia panel review after notebook migration

The old Layer-3 PySide6 desktop panel review chain is retired.  Current Radia
panel review targets the Jupyter notebook workbench:

```
src/radia/panels/notebooks/radia_<app>.ipynb
src/radia/*_design.py
src/radia/*_notebook.py
src/radia/notebook_workbench.py
src/radia/panels/calc_<app>.py
validation_test/panels/test_notebook_workbench.py
```

Review checklist:

1. `DesignSpec(...)` is the canonical initial-value store.
2. JSON files are run artifacts, not preset storage.
3. `Workbench.build_command()` matches the target `calc_*.py` argparse surface.
4. `run_local()` writes `radia_result.v2` with timing/version/runtime context.
5. Active notebooks do not import PySide6 / PyQt.
6. Human visualization uses `netgen.webgui`; headless validation uses durable
   `.msh v4.1`, JSON, and saved notebook outputs.

Run:

```powershell
python -m pytest validation_test/panels/test_notebook_workbench.py -q
```

Production deploy uses `radia[cubit]` plus `cubit-plugin-install`.  Do not
install PySide6 into normal Radia Python.  Coreform Cubit's bundled PySide6 is
protected and must not be removed.
"""


TOPICS = {
    "overview": PANEL_REVIEW_NOTEBOOK,
    "5_skills_chain": PANEL_REVIEW_NOTEBOOK,
    "13_checks": PANEL_REVIEW_NOTEBOOK,
    "bug_catalogue": PANEL_REVIEW_NOTEBOOK,
    "val_checkbox_trap": PANEL_REVIEW_NOTEBOOK,
    "map_value_reject": PANEL_REVIEW_NOTEBOOK,
    "widget_calc_gap": PANEL_REVIEW_NOTEBOOK,
    "smoke_scenarios": PANEL_REVIEW_NOTEBOOK,
    "red_flags": PANEL_REVIEW_NOTEBOOK,
    "workflow": PANEL_REVIEW_NOTEBOOK,
}


def get_panel_review_documentation(topic: str = "overview") -> str:
    """Return notebook-panel review guidance; keep old topic names stable."""
    if topic == "all":
        return "\n\n".join(f"# Topic: {k}\n{v}" for k, v in TOPICS.items())
    if topic in TOPICS:
        return TOPICS[topic]
    return (
        f"Unknown topic: {topic!r}\n\n"
        f"Available topics:\n"
        + "\n".join(f"  - {k}" for k in TOPICS)
    )
