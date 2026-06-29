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

## Current State

The current worktree snapshot is tracked in
[`docs/examples_classification/examples_classification.ipynb`](../docs/examples_classification/examples_classification.ipynb).
As of the 2026-06-29 migration check, the tracked non-`vim` contents of
`examples/` have been drained. The only remaining tracked files are this README
and the active `examples/vim/` research corpus.

| Topic | Current lane |
|---|---|
| `vim` | Active HDiv/VIM and FEEC research corpus. Move reusable helpers to `src`, validation locks to `validation_test/feec`, and public demonstrations to result-saved docs notebooks before pruning. |

Former non-`vim` topics have been promoted to their owning `docs/`, `src/`,
`validation_test/`, `panels/`, or `memory/` locations. If a deleted topic must
be recovered, use git history as the archive and recreate only the canonical
artifact in its final destination.

## Policy

Run the `promote-examples-to-docs` workflow before deleting anything:

1. Inventory with `rg --files`.
2. Search references in `docs`, `tests`, `validation_test`, `src`, and
   `packages`.
3. Create or refresh result-saved docs notebooks plus JSON sidecars.
4. Move reusable or validation code to its canonical home.
5. Update docs/MCP/panel references away from `examples/`.
6. Delete from `examples/` only after the owning artifact exists elsewhere.
