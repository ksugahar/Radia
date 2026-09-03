# Induction Heating Examples

This directory is the public notebook layer for induction heating. The old
loose examples topic is closed. Numerical evidence and its checked JSON remain
under `validation_test/`.

Start with the executed ESIM showcase:

- [`induction_heating_demo_showcase.ipynb`](induction_heating_demo_showcase.ipynb)

## Current Routing

- The notebook reproduces the analytical ESIM/Bessel cross-check and embeds
  selected publication figures.
- `validation_test/ih_esim_benchmark/` owns the benchmark scripts, checked JSON,
  sweep corpus, and source figures.
- Legacy `bem_reference/` has been split: reusable solver modules now live as
  `radia.bem_inductance`, `radia.bem_coupled_solver`, and `radia.ngsbem_*`;
  runnable reference scripts and sweep data live under
  `validation_test/induction_heating/bem_reference/`.
- `scattered_rhs_clean_test/` has been promoted to
  `validation_test/induction_heating/scattered_rhs_clean_test/` because it
  already carries a `.vol` fixture and `results.json`.
- Former demoted Cubit samples live under
  `validation_test/induction_heating/demoted_samples_legacy/`.
