# Kelvin Archive Retirement Map

This file is the maintained routing table for the remaining Kelvin source
archives.  The goal is to keep MCP knowledge, docs, source APIs, and validation
tests pointed at canonical maintained artifacts instead of depending on full
historical source ledgers.

## Maintained Surfaces

| Former archive | Maintained destination | Status |
|---|---|---|
| `kelvin_dtn_spectrum_archive.*` | `src/radia/open_boundary/{dtn_cln.py,kelvin_dtn.py}`, `src/radia/infinite_element.py`, `validation_test/open_boundary/`, `validation_test/kelvin_dtn_spectrum/` | Deleted after maintained API/validation routing was verified; use these destinations instead. |
| `kelvin_adaptive_mesh_archive.*` | `docs/kelvin/Supplement/{CG-smoother.md,ErrorEstimator.md}`, `docs/kelvin/Supplement/cg_smoother_demo.ipynb`, memory notes | Deleted after repetitive adaptive-mesh runners were collapsed to distilled method notes. |
| `kelvin_remaining_examples_archive.*` | `src/radia/kelvin_source.py`, `validation_test/kelvin_source/`, `validation_test/cubit/kelvin_1_4_p_convergence/`, `docs/kelvin/TEAM7_ADAPTIVE_RETIREMENT.md`, `memory/kelvin_remaining_examples_retired_2026_06_28.md` | Deleted after A-formulation pullbacks, Omega-Reduced Omega p-convergence, adaptive-method notes, and TEAM7 retirement rationale were routed to maintained destinations. |
| `kelvin_examples_migration.*`, `kelvin_classic_demos.*` | `docs/kelvin/{CONVENTION.md,KELVIN_TRANSFORMATION.md,kelvin_exterior_source_and_aphi.ipynb}`, `src/radia/kelvin_source.py`, `validation_test/kelvin_source/`, `validation_test/cubit/kelvin_1_4_p_convergence/` | Deleted after the completed migration ledger and full-source archive were replaced by maintained theory, API, notebook, and executable validation routes. |
| `packages/radia-mcp/examples/dtn_spectrum_coarse_mesh_demo.py` | `docs/kelvin/dtn_spectrum_coarse_mesh_demo.py`, `docs/kelvin/DTN_SPECTRUM_COARSE_MESH.md`, `radia_mcp.radia_ngsolve.knowledge.dtn_coarse_mesh` | Promoted out of the retired `examples/` tree during the radia-mcp diet; keep future runnable docs helpers next to their maintained docs. |

## MCP Rule

MCP knowledge should cite the maintained destination above.  It may mention that
the full historical source was once archived, but MCP tools should not require
those archive JSON files to answer ordinary Kelvin/Open-boundary questions.

## Deletion Gate

Delete a Kelvin archive only after all of these are true:

1. Public docs and MCP no longer point to the archive as the primary reference.
2. Any runnable claim has a maintained `validation_test/` or `src/` API path.
3. Debug-only lessons are distilled into `memory/` or concise docs.
4. The final reference scan confirms no non-archive file needs the deleted path.
