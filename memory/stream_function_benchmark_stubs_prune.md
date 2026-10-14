# Stream-Function Benchmark Stubs Prune

## Lesson

`examples/stream_function/benchmarks` should contain runnable benchmark code,
not TODO-only stubs that intentionally fail. Planned benchmark targets belong
in docs or memory until they are implemented and emit JSON results.

## Removed TODO Stubs

- `bench_bilac_planar_shim.py`
- `bench_turner_cylinder_gz.py`
- `bench_lemdiasov_ludwig_2005.py`
- `bench_coilgen_headtohead.py`
- `bench_shielded_iron_yoke.py`

## Kept

- `bench_helmholtz_pair.py`
- `benchmark_framework.py`
- `result_helmholtz.json`

When a target is implemented, add it back as a runnable benchmark with a
sidecar JSON and update `docs/stream_function/benchmarks.md`.
