# Solver Benchmark Validation Scripts

This directory contains long-running solver benchmark and validation drivers
that are too expensive for the default `tests/` loop.

## PEEC Dense vs HACApK

- `bench_peec_dense.py`
- `bench_peec_hacapk.py`
- `bench_peec_mna_crossover.py`
- `validate_peec_circuit_hacapk.py`

These scripts write their committed benchmark JSON to
`docs/solver_benchmarks/`, where the result-bearing notebook renders the
human-facing tables and plots.

## Legacy H-Matrix Solver Benchmarks

The older magnetostatic H-matrix benchmark scripts that used to live in
`examples/solver_benchmarks/` are kept here as validation drivers, not as
public examples:

- `benchmark_solver*.py`
- `benchmark_*construction.py`
- `benchmark_field_evaluation.py`
- `benchmark_lu_vs_hmatrix.py`
- `verify_field_accuracy.py`
- `run_all_benchmarks.py`
- `plot_benchmark_results.py`

Their historical notes are `legacy_hmatrix_README.md` and
`legacy_hmatrix_BENCHMARK_RESULTS.md`. Prefer the PEEC notebook in
`docs/solver_benchmarks/` for current public solver benchmark presentation.
