# Kelvin Source Validation

This directory is the maintained validation lane for Kelvin source and
material helper APIs in `radia.kelvin_source`.

The old `examples/kelvin_transformation/A-formulation` scripts were not kept
as runnable examples. Their durable checks are split here:

- `test_kelvin_material_factors.py` pins the 3D, axisymmetric, and 2D Kelvin
  material factors.
- `test_kelvin_pullback_forms.py` pins the 1-form A pullback, 2-form B
  pullback, local energy-density identity, and Kelvin-aware Biot-Savart
  source evaluation.

Use this directory, not the retired docs archive, when changing Kelvin
source-factor code.
