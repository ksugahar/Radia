# HDiv-MMM Cube Uniform-Field Validation

This corpus validates the current mesh-backed HDiv-VIM formulation on a
soft-iron cube in a uniform applied field. The maintained drivers and evidence
live in [`hdiv/`](hdiv/).

The former six-face Radia object benchmarks compared legacy `rad.Solve`
methods 0, 1, and 2. Those solver routes cannot represent the current
mesh-backed HDiv-MMM implementation and have been retired from this corpus.
Their source and measured outputs remain available in Git history.

Useful problem sizes are validation workloads, not pull-request tests. Run them
on hibino when available, or on mdx while its CI queue is idle, and preserve the
result as JSON under `validation_test/`.
