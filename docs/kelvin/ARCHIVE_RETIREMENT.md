# Kelvin Archive Retirement Map

This file is the maintained routing table for the remaining Kelvin source
archives.  The goal is to keep MCP knowledge, docs, source APIs, and validation
tests pointed at canonical maintained artifacts instead of depending on full
historical source ledgers.

## Maintained Surfaces

| Former archive | Maintained destination | Status |
|---|---|---|
| `kelvin_dtn_spectrum_archive.*` | `src/radia/open_boundary/{dtn_cln.py,kelvin_dtn.py}`, `src/radia/infinite_element.py`, `validation_test/open_boundary/`, `validation_test/kelvin_dtn_spectrum/` | Production API and validation already exist; archive is only historical source. |
| `kelvin_adaptive_mesh_archive.*` | `docs/kelvin/Supplement/{CG-smoother.md,ErrorEstimator.md}`, `docs/kelvin/Supplement/cg_smoother_demo.ipynb`, memory notes | Repetitive adaptive-mesh runners collapse to distilled method notes; promote only compact validation if a maintained regression is needed. |
| `kelvin_remaining_examples_archive.*` | `validation_test/cubit/kelvin_1_4_p_convergence/`, `src/radia/kelvin_*`, topic-specific future validation tests | Cubit p-convergence has a validation lane; the remaining A/Omega/TEAM7 scripts need one-by-one triage before deleting this archive. |

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
