# axifem validation corpus

This directory holds the runnable validation-class axifem checks promoted out
of the former axifem example tree.

- `research/verification/` keeps the Hiruma/Cauer and per-element validation
  scripts plus their committed JSON result records.
- `research/validate_q2_codegen.py` checks the closed-form Q2 generated matrix
  values against `research/q2_henrotte_test_values.json`.
- `axifem_element_evidence.json` is the checked aggregate consumed by the
  result-bearing public AXIFEM notebook.

The human-facing theory and executed notebooks live under `docs/axifem/`.
Development-history snapshots are intentionally not retained in docs.
