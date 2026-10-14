# Cubit Kelvin 1/4 P-Convergence

This validation fixture is the maintained home for the former
`examples/kelvin_transformation/Cubit_1_4_p_convergence` scripts.

It exercises the Omega-Reduced Omega Kelvin workflow on a Cubit-generated
quarter-sector magnetic sphere:

- `mesh_and_export.py` builds the Cubit mesh and exports `.vol` files at
  p-orders 1, 2, and 3.
- `solve_p_convergence.py` solves the Omega-Reduced Omega Kelvin problem and
  checks the field at the sphere center against the analytical sphere result.
- `validation_test/cubit/test_kelvin_1_4_p_convergence.py` is the runnable
  regression wrapper.

Use this fixture for high-order Kelvin periodic coupling and Cubit export
regression work. Do not restore the deleted example copy.
