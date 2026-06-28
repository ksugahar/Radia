# PEEC Solver Benchmarks

This directory is the public, result-bearing documentation layer for the
PEEC dense-Ruehli vs HACApK solver benchmarks.

## Canonical Artifacts

- `peec_solver_benchmarks.ipynb`: rendered notebook with saved outputs.
- `peec_solver_benchmarks_results.json`: synchronized debug record generated
  from the notebook and benchmark JSON.
- `results_bench_peec_*.json`: committed benchmark measurements used by the
  notebook.
- `findings_peec_mna_crossover.md`: written interpretation of the MNA
  crossover result.
- `comparison_peec_dense_vs_hacapk.{md,tex}`: table exports for papers and
  slide material.

## Runnable Layer

The executable benchmark drivers live in
`validation_test/solver_benchmarks/`. Running them refreshes the JSON in this
docs directory; re-run the notebook afterward so the saved outputs and
sidecar JSON stay synchronized.
