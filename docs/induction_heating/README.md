# Induction Heating Examples

This directory is the docs promotion layer for the remaining
`examples/induction_heating` scripts.

Start with the result-saved catalog notebook:

- `induction_heating_examples_catalog.ipynb`
- `induction_heating_examples_catalog_results.json`
- `induction_heating_examples_catalog_result.json`

The catalog stores the full source text, SHA-256 hashes, existing result/media
artifact hashes, protected references, and a migration lane for every current
example script.

## Current Routing

- ESIM tutorials (`esim_demo.py`, `esim_induction_heating_demo.py`,
  `demo_esim_impedance.py`) should become human-facing docs notebooks.
- `bem_reference/` is protected by panel tests and should be split into
  reusable `src` API plus runnable `validation_test` checks before deletion.
- `scattered_rhs_clean_test/` is a validation-test candidate because it already
  carries a `.vol` fixture and `results.json`.
- `demoted_samples/` is panel/sample history. Do not delete it until the
  corresponding notebook-panel or validation fixture replaces each sample.
