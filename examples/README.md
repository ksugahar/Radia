# Radia Examples Migration Area

`examples/` is no longer the long-term home for public demos, reusable APIs, or
validation corpora. It is a temporary research/teaching tier being drained into
the canonical destinations:

| Final destination | Use |
|---|---|
| `docs/<topic>/*.ipynb` | Public method showcase with Markdown, saved outputs, and synchronized JSON. |
| `src/` | Reusable computation, parsers, mesh readers, formulas, and solver helpers. |
| `validation_test/` | Golden locks, convergence sweeps, benchmarks, and heavier validation surfaces. |
| `panels/` or `src/radia/panels/` | Panel-owned operating assets and samples during the staged panel migration. |
| `memory/` or docs note, then delete | Superseded experiments and failed development iterations. |

Do not add new long-lived references to `examples/`. If a test, notebook, MCP
knowledge file, or panel still points here, treat that as migration debt and
record the intended `target_after_unblock`.

## Current Large Topics

The current worktree snapshot is tracked in
[`docs/examples_classification/examples_classification.ipynb`](../docs/examples_classification/examples_classification.ipynb).
As of the latest migration batch, the largest remaining topics are:

| Topic | Current lane |
|---|---|
| `maglev` | Split into `validation_test`, `docs`, `src`, and distill-delete. |
| `peec_integration` | Move validation corpus first, then promote reusable readers/generators to `src`. |
| `vim` | Move `validation_test/feec` direct imports first; promote reusable HDiv/VIM helpers to `src`. |
| `clebsch_hodograph` | Use the validation-test research harness as migration driver; keep docs notebooks as showcase. |
| `cube_uniform_field` | Move benchmark drivers and result corpus under validation/docs ownership. |
| `cubit_panels` | IH inductance scripts moved to `validation_test/induction_heating/cubit_panels_legacy`; remaining accel-magnet sources still need `src`/docs/panels ownership. |
| `stream_function` | Extract reusable code to `src`, keep result-saved docs notebooks as the public layer. |
| `induction_heating` | Promote BEM helpers to `src`; move checks to `validation_test`; keep ESIM demos in docs. |

## Policy

Run the `promote-examples-to-docs` workflow before deleting anything:

1. Inventory with `rg --files`.
2. Search references in `docs`, `tests`, `validation_test`, `src`, and
   `packages`.
3. Create or refresh result-saved docs notebooks plus JSON sidecars.
4. Move reusable or validation code to its canonical home.
5. Update docs/MCP/panel references away from `examples/`.
6. Delete from `examples/` only after the owning artifact exists elsewhere.
