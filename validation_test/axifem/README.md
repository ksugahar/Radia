# axifem validation corpus

This directory holds the runnable validation-class axifem checks promoted out
of the former axifem example tree.

- `research/verification/` keeps the Hiruma/Cauer and per-element validation
  scripts plus their committed JSON result records.
- `research/validate_q2_codegen.py` checks the closed-form Q2 generated matrix
  values against `research/q2_henrotte_test_values.json`.
- `axifem_element_evidence.json` is the checked aggregate consumed by the
  element evidence notebook in this directory.
- `test_heat_axisym_e2e.py` exercises the production thermal CLI at Q1 and Q2.
  `validate_heat_axisym_q2.py` owns the corresponding Bessel-series comparison
  and regenerates `results_heat_axisym_q2_near_axis_20260903.json`.

The [element evidence notebook](AXIFEM_ELEMENT_EVIDENCE.ipynb) lives beside its
checked JSON here. It displays stored results; opening it does not rerun the
numerical acceptance suite. Human-facing theory remains under `docs/axifem/`.
Development-history snapshots are intentionally not retained in docs.
