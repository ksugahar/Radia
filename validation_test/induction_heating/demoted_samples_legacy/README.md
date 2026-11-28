# Demoted IH samples

These files used to live in `src/radia/panels/samples/` and shipped in the
`radia` PyPI wheel as panel Browse-dialog choices. On **2026-04-23** they
were demoted out of the shipped panel samples because they violate the Panel
Samples Quality Policy (AGENTS.md § "Sample Promotion Ladder"). On
**2026-06-29** they were moved here from the loose examples tier so validation
history no longer depends on examples paths:

- **Non-canonical**: no evidence test under `validation_test/panels/`
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
| `ih_fem_kelvin_sample.jou`        | Small-mesh FEM+Kelvin demo | Use `ih_fem_kelvin_skin_fine.jou` (golden-tested) |
| `ih_fem_sample.jou`                | No-Kelvin FEM baseline (truncation demo) | — (educational only) |
| `ih_fem_kelvin_skin.jou`          | Coarse `_skin` variant | Use `ih_fem_kelvin_skin_fine.jou` |
| `ih_closed_torus.jou`             | Closed (no-gap) torus research variant | — (research-stage) |

## How to use these now

They still work for archaeology and validation experiments; they are just no
longer automatically installed into the Cubit panel's Browse dialog. The
journals build the physical domain only and do *not* add Kelvin. To export
one with the Kelvin open boundary, pass it to the MCP tool
`cubit_export_vol(<file.jou>, "model.vol", order=1, add_kelvin=True)` or use
the Kelvin option of the Vol Export menu; both run
`radia_mcp.cubit.kelvin` before `export vol`.
