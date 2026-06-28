# Demoted IH samples

These files used to live in `src/radia/panels/samples/` and shipped in the
`radia` PyPI wheel as panel Browse-dialog choices. On **2026-04-23** they
were demoted out of the shipped panel samples because they violate the Panel
Samples Quality Policy (AGENTS.md § "Sample Promotion Ladder"). On
**2026-06-29** they were moved here from the loose examples tier so validation
history no longer depends on examples paths:

- **Non-canonical**: no golden test under `tests/panels/test_*_golden.py`
  references them, so their numerical correctness is not automatically
  verified on every release.
- **Misleading "auto-Kelvin" comments**: several of these `.jou` files
  claimed Kelvin open-boundary was added automatically by the panel.
  That feature was never implemented; running the `.jou` alone produced
  a `.vol` *without* a `kelvin` material, and `calc_fem_coilmesh.py` then
  silently fell back to a Dirichlet A=0 truncation (~3 % L bias,
  ~20 % P_wp bias on the gapped torus geometry).

## Contents

| File | Was used for | Canonical replacement |
|---|---|---|
| `ih_fem_kelvin_sample.jou` / `.py` | Small-mesh FEM+Kelvin demo | Use `ih_fem_kelvin_skin_fine.jou` (golden-tested) |
| `ih_fem_sample.jou`                | No-Kelvin FEM baseline (truncation demo) | — (educational only) |
| `ih_fem_kelvin_skin.jou` / `.py`   | Coarse `_skin` variant | Use `ih_fem_kelvin_skin_fine.jou` |
| `ih_closed_torus.jou` / `.py`      | Closed (no-gap) torus research variant | — (research-stage) |

## How to use these now

They still work for archaeology and validation experiments — they're just no
longer automatically installed into the Cubit panel's Browse dialog. To run one:

1. From the Radia source tree, copy the `.jou` / `.py` to a working
   directory of your choice (e.g. `C:\temp\`).
2. If you want Kelvin open-boundary, **always pick the `.py`** — it
   explicitly calls `add_kelvin_cubit()`.  The `.jou` files in this
   directory do *not* add Kelvin.
3. In Cubit: `Play → <file.py>` (or `<file.jou>`).
4. Export: `export netgen "model.vol" order 1 overwrite`.

The `.py` loader finds `add_kelvin.py` via two paths: repo-relative
(`../../../src/radia/panels/add_kelvin.py`) or the pip-installed radia
package (`radia.panels.add_kelvin`).  Either location works.

## Why the `.jou` vs `.py` split exists

Cubit `.jou` is a journal (text command stream).  It cannot import
Python modules directly — only via the `python "..."` escape.  The
companion `.py` files are Cubit-side Python scripts that call
`cubit.cmd(...)` plus `add_kelvin_cubit()`, giving a one-shot
"geometry + Kelvin + export" workflow.

For panel use, Cubit users prefer `.jou` (familiar, replayable), so any
sample that needs Kelvin must either (a) inline the Kelvin Cubit Python
or (b) ship as a `.py`.  These demoted samples took the second route
and then ended up confusing users who double-clicked the `.jou` by
mistake — hence the demotion.
