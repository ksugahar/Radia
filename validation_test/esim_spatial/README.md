# Retained ESIM spatial evidence

The [documentation example](../../docs/esim_spatial/esim_spatial_demo.ipynb)
uses these three weak-coupled, 100 A results, computed on mdx2 on 2026-09-15.
`compute_record.json` retains runtime versions, original source hashes and
input hashes. Recovery replaced machine-local paths with repository-relative
paths and removed embedded source-code duplicates. Original result hashes are
retained separately from the hashes of the portable copies.

The whole-file mesh check used 10% tolerance because the unused air and coil
volume discrepancies exceed 1%. The workpiece volume and SIBC area each satisfy
1%. This boundary-only demonstration does not relax volume-FEM acceptance.

`spatial_metrics.json` distinguishes the local surface diagnostic from the
power-normalized transfer artifact. Saved outputs are historical evidence;
they are not a new benchmark or proof of spatial-loss accuracy. Rendered visual
QA remains unverified.

`test_saved_evidence.py` checks input/result integrity, convergence records and
the declared mesh/observable scope without rerunning a solver.
