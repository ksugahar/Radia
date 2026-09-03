# PEEC Integration Docs

The public entry point is the executed numerical showcase
[`peec_showcase.ipynb`](peec_showcase.ipynb). It contains the maintained Dowell,
ngsolve.bem, and publication-figure demonstrations without a docs-side result
JSON. Checked numerical evidence belongs under
`validation_test/peec_integration/`.

The `docs/peec_integration/demos/` tree is now the docs-coupled public demo
corpus. Benchmark, comparison, smoke, and solver-regression scripts live under
`validation_test/peec_integration`; reusable behavior belongs in `src/radia`.
Completed source-routing ledgers are retained by Git history, not as public
notebooks.
