# Cubit Mesh Export Validation Corpus

Runnable validation scripts for the Cubit `.vol` export live here.

Presentation/demo material is kept under `docs/cubit_mesh_export`; these files
are the regression-oriented checks for Netgen `.vol` boundary inventory,
surface quality, pressure/traction resultants, and FEM/BEM topology handling.
`cubit_mesh_export_showcase_results.json` aggregates that checked evidence for
the executed public showcase notebook without placing a duplicate result JSON
under `docs/`.
`p_convergence_demo_results.json` records the Cubit-to-`.vol` p-convergence
run consumed by the executed docs notebook of the same name.

## Geometric curving refit

`refit_map_derivatives.py` measures the higher derivatives of the curved
element map on `.vol` files exported from the same Cubit meshes with Netgen's
coefficients and with the geometric refit, and writes
`refit_map_derivatives_results.json` next to the script.
