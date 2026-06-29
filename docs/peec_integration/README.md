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

The completed verification migration is
[`verification_migration.ipynb`](verification_migration.ipynb), synchronized with
[`verification_migration_results.json`](verification_migration_results.json) and
[`verification_migration_result.json`](verification_migration_result.json). It
moved the PEEC verification corpus to
`validation_test/peec_integration/verification`, promoted the reusable GMSH
centerline reader to `radia.peec_mesh_import`, and removed the distilled
`check_funcs.py` scratch file from source.

The remaining `examples/peec_integration/` tree is still a mixed staging area:
maintained public demos, source-API candidates, benchmark scripts, experiments,
smoke checks, and analysis helpers. Do not treat every script as a public demo.
Use the cleanup-review table in `public_demo.ipynb` and the migration notebook
above before rerouting any remaining examples.
