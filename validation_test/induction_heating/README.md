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

## BEM Validation

[`test_coil_bem_a_volume_vol.py`](../bem/test_coil_bem_a_volume_vol.py)
contains executable BEM validation. Reusable solver code is provided by
`radia.bem_inductance`, `radia.bem_coupled_solver`, and `radia.ngsbem_*`;
validation scripts should import those APIs directly. Retired investigation
scripts and unreproducible sweep records remain available in Git history.

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

Rerun record (2026-11-27, regenerated on NGSolve 6.2.2608 and radia 5.4.6;
previous record NGSolve 6.2.2607, radia native 5.1.0). Thresholds unchanged,
PASS. Every difference is below its threshold by orders of magnitude:

- Thermal entries (row-sum error, stiffness maximum and nullspace, mass):
  cause identified. Netgen 6.2.2608 places one interior box node 2.1e-11 m
  away from its 6.2.2607 position, and the summation order changes with it.
  An A/B swap of element order and node coordinate selects the old bit
  pattern by the coordinate alone; the remaining spread is run-to-run.
- Electromagnetic power (and the thermal source power and
  `relative_power_error` derived from it): 5.771847527664258e-06 W before,
  5.771847527664283e-06 W now (25 ulp, 4.3e-15 relative). Cause not
  identified. Twelve runs on both meshes and both swapped meshes stay within
  5.771847527664280e-06 to 5.771847527664285e-06 W, so the mesh does not
  explain the old value. The remaining candidates (radia native 5.1.0 vs
  5.4.6, the old uncommitted source, the NGSolve runtime under the EM route)
  cannot be separated, because radia native 5.1.0 is no longer available.
  The new value is recorded as measured; no threshold was changed.

## Demoted Sample Legacy Fixtures

`demoted_samples_legacy/` contains the old IH Cubit journals and Cubit-side
Python generators that are no longer shipped as panel samples and no longer
belong in `examples/`. They are kept here for reproducible archaeology of
Kelvin open-boundary setup, SIBC hole tagging, and closed-torus edge cases.
