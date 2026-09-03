# PEEC Solver Benchmark Validation

This directory owns long-running PEEC solver benchmarks that are too expensive
for the default `tests/` loop. The executable drivers write JSON records here;
those records, rather than a validation notebook, are the durable numerical and
performance evidence.

## Dense vs HACApK

- `bench_peec_dense.py`
- `bench_peec_hacapk.py`
- `bench_peec_mna_crossover.py`
- `validate_peec_circuit_hacapk.py`

The committed records are `results_bench_peec_dense.json`,
`results_bench_peec_hacapk.json`, `results_bench_peec_mna_crossover.json`, and
the aggregate `peec_solver_benchmarks_results.json`. These benchmarks are
manual validation work and are not part of routine pull-request CI.

The retired mesh-less magnetostatic `rad.Solve` method comparisons do not
belong to this lane. Legacy methods 1 and 2 are no longer public solver routes.
Current magnetostatic HDiv-MMM and HACApK evidence lives under
`validation_test/feec/`, with coupled HDiv-MMM/HCurl evidence under
`validation_test/vim_coupled/`.
