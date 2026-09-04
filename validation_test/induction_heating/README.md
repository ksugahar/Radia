# Induction Heating Validation

This directory holds induction-heating checks promoted out of `examples/`.

## ESIM Cross-Formulation Consistency

`esim_cross_formulation/` checks that the nonlinear per-element ESIM
`P_wp` agrees between the scalar-potential BIE
(`calc_inductance.py --impedance-model esim`) and the HCurl A + Kelvin
FEM (`calc_fem_kelvin.py --impedance esim`) implementations -- both
driving the same 1-D cell solver and Karl loop -- on two cylinders
(BH-knee dia-20 and deep-saturation dia-50). `make_meshes.py`
regenerates the meshes (Cubit); `run_cross_formulation.py` runs the
comparison matrix and asserts golden bands (no Cubit import). Backs the
validation tier (vi) of the SA-26-070 (八戸) and IGTE 2026 papers.

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

## Simulink Operator Physics Golden

`test_ih_operator_physics_golden.py` runs the real PEEC/BEM-SIBC unit-current
solve and verifies thermal power closure, assembled-operator identities, and
the exact single-current response contract used by the IH Simulink block.

```powershell
python -m pytest validation_test/induction_heating/test_ih_operator_physics_golden.py -q
python validation_test/induction_heating/generate_ih_operator_physics_results.py
```

The generator owns `ih_operator_physics_results.json`; do not edit its
numerical values by hand.

## Legacy BEM/ESIM Validation

`cubit_panels_legacy/` is a migration corpus from the former IH Cubit panel.
It retains research scripts and indispensable Cubit input fixtures only while
their numerical claims are promoted to checked JSON evidence. Reusable kernels
and application behavior belong to `src/radia`; generated GMSH companion files
are run artifacts and are not tracked here. Do not add new production imports
or examples-path references to this directory.

## Demoted Sample Legacy Fixtures

`demoted_samples_legacy/` contains the old IH Cubit journals and Cubit-side
Python generators that are no longer shipped as panel samples and no longer
belong in `examples/`. They are kept here for reproducible archaeology of
Kelvin open-boundary setup, SIBC hole tagging, and closed-torus edge cases.
