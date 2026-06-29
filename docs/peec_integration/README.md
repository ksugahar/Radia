# PEEC Integration Docs

The public entry point is
[`public_demo.ipynb`](public_demo.ipynb), synchronized with
[`public_demo_results.json`](public_demo_results.json) and
[`public_demo_result.json`](public_demo_result.json).

The verified numerical showcase is
[`peec_showcase.ipynb`](peec_showcase.ipynb), synchronized with
[`peec_showcase_result.json`](peec_showcase_result.json).

The full examples source/hash catalog is
[`examples_catalog.ipynb`](examples_catalog.ipynb), synchronized with
[`examples_catalog_results.json`](examples_catalog_results.json) and
[`examples_catalog_result.json`](examples_catalog_result.json).

The cleanup route for the 63 non-public-demo scripts is
[`cleanup_routing.ipynb`](cleanup_routing.ipynb), synchronized with
[`cleanup_routing_results.json`](cleanup_routing_results.json) and
[`cleanup_routing_result.json`](cleanup_routing_result.json).

The completed cleanup of the remaining examples tree is
[`post_examples_migration.ipynb`](post_examples_migration.ipynb), synchronized
with [`post_examples_migration_results.json`](post_examples_migration_results.json)
and [`post_examples_migration_result.json`](post_examples_migration_result.json).
It records the 97 Python scripts and result assets that were moved from
`docs/peec_integration/demos` into `docs/peec_integration/demos` or
`validation_test/peec_integration`, plus the one-off scripts that were
distilled and removed.

The completed verification migration is
[`verification_migration.ipynb`](verification_migration.ipynb), synchronized with
[`verification_migration_results.json`](verification_migration_results.json) and
[`verification_migration_result.json`](verification_migration_result.json). It
moved the PEEC verification corpus to
`validation_test/peec_integration/verification`, promoted the reusable GMSH
centerline reader to `radia.peec_mesh_import`, and removed the distilled
`check_funcs.py` scratch file from source.

The `docs/peec_integration/demos/` tree is now the docs-coupled public demo
corpus. Benchmark, comparison, smoke, and solver-regression scripts live under
`validation_test/peec_integration`; reusable behavior belongs in `src/radia`.
