# Induction Heating Examples

This directory is the docs promotion layer for the induction-heating example
cleanup. The old loose examples topic is closed; public demonstrations now live
as executed, result-bearing notebooks, while numerical evidence and its JSON
live under `validation_test/induction_heating/`.

Start with the public result-saved notebook. The adjacent JSON files are
historical migration records and are not required by the current docs policy:

- `public_demo.ipynb`
- `public_demo_results.json`
- `public_demo_result.json`

The full source/hash catalog is:

- `induction_heating_examples_catalog.ipynb`
- `induction_heating_examples_catalog_results.json`
- `induction_heating_examples_catalog_result.json`

The closed public demo showcase is:

- `induction_heating_demo_showcase.ipynb`
- `induction_heating_demo_showcase_results.json`
- `induction_heating_demo_showcase_result.json`

The catalogs store source hashes, route decisions, existing result/media
artifact hashes, protected references, and the final migration lane.

## Current Routing

- ESIM/WPT/RWG tutorials are represented by
  `induction_heating_demo_showcase.ipynb` and its JSON source/hash archive.
- Legacy `bem_reference/` has been split: reusable solver modules now live as
  `radia.bem_inductance`, `radia.bem_coupled_solver`, and `radia.ngsbem_*`;
  runnable reference scripts and sweep data live under
  `validation_test/induction_heating/bem_reference/`.
- `scattered_rhs_clean_test/` has been promoted to
  `validation_test/induction_heating/scattered_rhs_clean_test/` because it
  already carries a `.vol` fixture and `results.json`.
- Former demoted Cubit samples live under
  `validation_test/induction_heating/demoted_samples_legacy/`.
