# Induction Heating Validation

This directory holds induction-heating checks promoted out of `examples/`.

## Scattered RHS Clean Test

`scattered_rhs_clean_test/` is a historical validation fixture for the
scattered-field Robin RHS investigation. It keeps the script, the generated
Netgen `.vol` mesh fixture, and the measured `results.json` together.

## BEM Reference

`bem_reference/` contains the executable reference scripts and sweep result
that used to live under the induction-heating examples tree.  The
reusable solver code was promoted to `src/radia` as `radia.bem_inductance`,
`radia.bem_coupled_solver`, `radia.cubit_bem_extractor`, and
`radia.ngsbem_*`; validation scripts should import those APIs instead of
reaching back into `examples/`.

The public docs catalog is `docs/induction_heating/`.

## Cubit Panels Legacy

`cubit_panels_legacy/` contains the former IH Cubit-panel Python and
Cubit/display assets.  Keep it as a protected validation corpus while reusable
kernels are promoted to `src/` and public explanations move into result-saved
docs notebooks.  Do not add new examples-path references back to the old IH
Cubit-panel examples location.

## Demoted Sample Legacy Fixtures

`demoted_samples_legacy/` contains the old IH Cubit journals and Cubit-side
Python generators that are no longer shipped as panel samples and no longer
belong in `examples/`. They are kept here for reproducible archaeology of
Kelvin open-boundary setup, SIBC hole tagging, and closed-torus edge cases.
