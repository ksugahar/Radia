# PEEC Integration Docs

The public entry point is the executed, result-bearing
[`public_demo.ipynb`](public_demo.ipynb). Adjacent JSON files are retained as
optional historical records, not as the current docs contract.

The numerical showcase is [`peec_showcase.ipynb`](peec_showcase.ipynb).
Checked numerical evidence belongs under `validation_test/peec_integration/`.

The full examples source/hash catalog is the historical
[`examples_catalog.ipynb`](examples_catalog.ipynb).

The cleanup route for the 63 non-public-demo scripts is recorded in the
historical [`cleanup_routing.ipynb`](cleanup_routing.ipynb).

The completed cleanup of the remaining examples tree is recorded in the
historical [`post_examples_migration.ipynb`](post_examples_migration.ipynb).
It records the 97 Python scripts and result assets that were moved from
`docs/peec_integration/demos` into `docs/peec_integration/demos` or
`validation_test/peec_integration`, plus the one-off scripts that were
distilled and removed.

The completed verification migration is recorded in the historical
[`verification_migration.ipynb`](verification_migration.ipynb). It
moved the PEEC verification corpus to
`validation_test/peec_integration/verification`, promoted the reusable GMSH
centerline reader to `radia.peec_mesh_import`, and removed the distilled
`check_funcs.py` scratch file from source.

The `docs/peec_integration/demos/` tree is now the docs-coupled public demo
corpus. Benchmark, comparison, smoke, and solver-regression scripts live under
`validation_test/peec_integration`; reusable behavior belongs in `src/radia`.
