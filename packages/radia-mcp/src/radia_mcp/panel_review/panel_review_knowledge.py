"""Radia notebook panel review and construction knowledge."""

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


NOTEBOOK_GUI_BUILD = """
# Build a Radia Jupyter notebook GUI

Use this when promoting an `examples/<topic>/*.py` script, a former PySide
panel, or a Cubit-panel prototype into the current browser-native Radia panel
surface.  The construction target is always:

```
panels/calc_<app>.py                  # target layout for panel CLI scripts
panels/notebooks/radia_<app>.ipynb    # target layout for notebook panels
panels/samples/<app>/...              # panel-owned samples and fixtures
src/radia/<app>_design.py
src/radia/<app>_notebook.py
src/radia/notebook_workbench.py
validation_test/panels/test_notebook_workbench.py
```

The notebook GUI is thin: it maps `DesignSpec` fields into ipywidgets and
launches the validated `calc_*.py` command.  It does not re-implement solver
logic.

Construction checklist:

1. Move reusable solver, geometry, parser, or mesh-helper code to `src/`.
2. Keep heavy numerical checks in `validation_test/<topic>/`; add a docs
   notebook only as the human-facing showcase layer.
3. Define `<App>DesignSpec` as the canonical initial-value store.  JSON is a
   run artifact, not a preset store.
4. Implement `DesignSpec.build_command()`, `missing_required_inputs()`, and
   `visible_fields()` where needed.  The command must target a headless
   `panels/calc_*.py` script with argparse flags.  During the staged migration,
   existing `src/radia/panels/calc_*.py` scripts remain legacy-compatible.
5. In `<app>_notebook.py`, declare `NotebookFieldSpec(...)` rows, subclass
   `CommandWorkbench`, set `title`, `field_specs`, `section_order`, and an
   app-specific `run_root` such as `runs/radia_ih`.
6. In `radia_<app>.ipynb`, include concise Markdown notes plus a code cell:

```python
from radia.<app>_design import <App>DesignSpec
from radia.<app>_notebook import <App>Workbench

spec = <App>DesignSpec()
workbench = <App>Workbench(spec)
workbench.display()
```

7. Register the notebook in `panels/notebooks/panel_notebook_manifest.json`
   once the root panel tree owns it.  During the staged migration, keep the
   legacy manifest in `src/radia/panels/notebooks/` synchronized.
8. Extend `validation_test/panels/test_notebook_workbench.py` so the workbench
   builds the expected headless command and the notebook has no PySide/PyQt
   imports.
9. Run:

```powershell
python -m pytest validation_test/panels/test_notebook_workbench.py -q
```

Run artifacts:

- `CommandWorkbench.run_local()` writes a timestamped run directory with
  `command.txt`, `run.log`, and `result.json`.
- `result.json` starts with `radia_result`, schema
  `radia.notebook_panel_run.v2`, Radia/Python/platform versions, status,
  command, wall time, and the four heaviest timing stages discoverable from
  solver JSON outputs.
- Docs/showcase notebooks must be committed with outputs and synchronized
  adjacent JSON containing `generated_at_utc`, version/runtime data, and a
  matching `notebook_sha256`.

Presentation template:

The NGSolve User Meeting draft
`public-safe curated corpus`
is a presentation shell for the IH workbench.  Its pattern is:

1. Markdown title explaining that the notebook is the panel.
2. Optional dark/presentation CSS.
3. `IHDesignSpec` + `IHWorkbench` display cell.
4. Short tips about required inputs and saved run artifacts.

Use that shell for talks and operator-facing demos; keep the repository
notebook contract in `panels/notebooks/radia_ih.ipynb` once migrated.  The
current legacy source-of-truth path is
`src/radia/panels/notebooks/radia_ih.ipynb`.

CSS / JupyterLab guard:

Presentation CSS may restyle colors, typography, and widget backgrounds, but
must not change pointer events, z-index, positioning, or layout ownership of
Jupyter cells and input areas.  If hovering the notebook run/play button makes
it translucent or unclickable, a cell-selection layer is covering the toolbar;
scope the CSS to visual properties and verify that both the cell run button and
the workbench `Run` button remain clickable.

Cubit boundary:

The Cubit Export Mesh toolbar is still a Cubit-embedded plugin surface and may
use Coreform Cubit's bundled PySide6.  Normal Radia Python panels are notebook
workbenches and must not depend on PySide6/PyQt.
"""


CUBIT_PANELS_MIGRATION = """
# `examples/cubit_panels` migration plan

As of 2026-06-29, the induction-heating side has been drained from
`examples/cubit_panels/inductance` into
`validation_test/induction_heating/cubit_panels_legacy`.  The accel-magnet
staging scripts were later pruned from `examples`; only rescued panel fixtures
belong under `src/radia/panels/samples/em/c_type_dipole/`.  Treat any remaining
`examples/cubit_panels` path as reference debt to eliminate, not as a
long-lived examples tier.

## Destination rules

- Reusable geometry/model builders -> `src/radia/...` API plus focused tests.
- Numerical checks, golden locks, and solver comparisons -> `validation_test`.
- User-facing demonstrations -> result-saved docs notebooks plus synchronized
  JSON.
- Presentation/operator GUI and panel-owned samples -> repo-root `panels/`
  using `DesignSpec` + `CommandWorkbench`.
- Cubit journals, `.geo`, BH tables, and mesh/CAD assets remain protected only
  while an owning API/notebook/test still needs them.  Do not keep stale
  development logs as examples.

## Accel magnet side

The former `examples/cubit_panels/accel_magnet` tree was source material for
the EM notebook/panel track:

- `coil_dipole.py` duplicated the panel coil sample and was pruned.
- `experiment_mmm_ima.py`, `experiment_occ_dipole.py`, and STEP/FEM helper
  probes were development diagnostics; their lesson lives in memory, not in
  the source tree.
- `coil_wire.step` and `yoke.step` were rescued to
  `src/radia/panels/samples/em/c_type_dipole/`.
- `BH.txt` is embedded in `src/radia/em_material.py` and shipped as
  `src/radia/panels/samples/em_sample_bh.txt`; the examples copy was removed.

## Induction heating side

The original `examples/cubit_panels/inductance` Python and display-asset corpus
now lives under `validation_test/induction_heating/cubit_panels_legacy`.  This
is a protected legacy validation corpus, not a final public docs surface:

- Validation-first: `compare_bem_coupled_vs_fem_kelvin.py`, `verify_*.py`, and
  `test_*.py` belong in `validation_test/induction_heating` or a more specific
  validation subtree.
- Docs/showcase: `scalar_bie_sibc.py`, `bem_sibc_workpiece.py`,
  `efie_sibc.py`, `pmchwt_sibc.py`, `mfie_sphere_demo.py`,
  `fem_esim_3d.py`, `fem_esim_kelvin.py`, `fem_total_field.py`, and
  `impedance_esim.py` should become result-saved notebooks only after reusable
  kernels are in `src`.
- API candidates: `create_induction_model.py`, `fem_esim_3d_cubit.py`,
  `inductance_hodge.py`, `inductance_source_sink.py`, and shared torus/coil
  builders should move from the legacy corpus to `src/radia` only when they
  become reusable APIs.  Panel-only wiring belongs under `panels/`.
- Display `.geo` files are visualization assets; keep or regenerate them next
  to the notebook/test that owns them.

## Original 35-script inventory

Route each Python script before deleting the corresponding legacy copy:

| Script | Target after unblock |
|--------|----------------------|
| `accel_magnet/coil_dipole.py` | pruned; duplicate of panel coil sample |
| `accel_magnet/experiment_mmm_ima.py` | pruned after memory distillation |
| `accel_magnet/experiment_occ_dipole.py` | pruned; superseded by panel calc/sample |
| `accel_magnet/experiment_step_fem.py` | pruned; superseded by `step_mesh_builder.py` and panel calc |
| `accel_magnet/experiment_step_fem_full.py` | pruned after memory distillation |
| `accel_magnet/experiment_step_fem_nokelvin.py` | pruned after memory distillation |
| `inductance/bem_sibc_workpiece.py` | IH docs notebook plus src kernel |
| `inductance/compare_bem_coupled_vs_fem_kelvin.py` | `validation_test` |
| `inductance/create_induction_model.py` | `src` Cubit/IH model API |
| `inductance/efie_sibc.py` | IH docs notebook plus validation |
| `inductance/experiment_bem_sibc_solver.py` | `validation_test` |
| `inductance/experiment_coupled_bem.py` | IH docs notebook plus validation |
| `inductance/experiment_coupled_bem_steel.py` | IH docs notebook plus validation |
| `inductance/experiment_surface_J_accuracy.py` | `validation_test` |
| `inductance/fem_esim_3d.py` | IH docs notebook plus src kernel |
| `inductance/fem_esim_3d_cubit.py` | `src` Cubit/IH API plus validation |
| `inductance/fem_esim_kelvin.py` | IH docs notebook plus validation |
| `inductance/fem_scattered_coil.py` | `src` incident-field/coil helper |
| `inductance/fem_total_field.py` | IH docs notebook plus src helper |
| `inductance/impedance_esim.py` | IH docs notebook plus src kernel |
| `inductance/inductance_hodge.py` | `src` BEM inductance API |
| `inductance/inductance_source_sink.py` | `src` BEM inductance API |
| `inductance/inductance_torus.py` | `src`/panel calc or Cubit fixture |
| `inductance/mfie_sphere_demo.py` | IH/BEM docs notebook plus validation |
| `inductance/pmchwt_sibc.py` | IH/BEM docs notebook plus validation |
| `inductance/pmchwt_sibc_test.py` | `validation_test` (renamed to `validation_pmchwt_sibc.py` to avoid pytest auto-collection) |
| `inductance/scalar_bie_sibc.py` | IH docs notebook plus validation |
| `inductance/test_interp_quality.py` | `validation_test` (renamed to `validation_interp_quality.py`) |
| `inductance/test_nxH_rhs.py` | `validation_test` (renamed to `validation_nxH_rhs.py`) |
| `inductance/verify_esim.py` | `validation_test` |
| `inductance/verify_laplace_bem.py` | `validation_test` |
| `inductance/verify_per_node_sibc_sphere.py` | `validation_test` |
| `inductance/verify_per_node_sibc_spheroid.py` | `validation_test` |
| `inductance/verify_per_node_sibc_torus.py` | `validation_test` |
| `inductance/verify_sphere_sibc.py` | `validation_test` |

## Deletion gate

After each promotion batch, run a repository reference search for
`examples/cubit_panels`.  Remaining references must name a concrete
`target_after_unblock` (`src`, `validation_test`, `docs`, `panels`, or delete
after distillation).  New long-lived references to examples are not allowed.
"""


TOPICS = {
    "overview": PANEL_REVIEW_NOTEBOOK,
    "build_notebook_gui": NOTEBOOK_GUI_BUILD,
    "presentation_template": NOTEBOOK_GUI_BUILD,
    "cubit_panels_migration": CUBIT_PANELS_MIGRATION,
    "5_skills_chain": PANEL_REVIEW_NOTEBOOK,
    "13_checks": PANEL_REVIEW_NOTEBOOK,
    "bug_catalogue": PANEL_REVIEW_NOTEBOOK,
    "val_checkbox_trap": PANEL_REVIEW_NOTEBOOK,
    "map_value_reject": PANEL_REVIEW_NOTEBOOK,
    "widget_calc_gap": PANEL_REVIEW_NOTEBOOK,
    "smoke_scenarios": PANEL_REVIEW_NOTEBOOK,
    "red_flags": PANEL_REVIEW_NOTEBOOK,
    "workflow": NOTEBOOK_GUI_BUILD,
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
