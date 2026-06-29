# Stream-function benchmark suite

Reproducible benchmark cases for validating the Radia SF inverse-design
pipeline against published references and other open-source coil design
tools. Each benchmark loads a literature-defined target spec, runs the
full pipeline (SF solve → contour → single-stroke → field), and reports
the key metrics in a single JSON file (so the table can be regenerated
unattended).

## What's here

| Status | File | Reference |
|--------|------|-----------|
| shipped | `bench_helmholtz_pair.py` | Analytical Helmholtz pair (uniform Bz, Maxwell) |

This directory was promoted out of the examples tier.
Previously that tier also held TODO-only benchmark stubs. Those
were removed because they were not runnable examples. The target list is kept
in `docs/stream_function/benchmarks.md` and
`memory/stream_function_benchmark_stubs_prune.md` until a target is implemented
with JSON output.

## How to run

```bash
# Single benchmark
python bench_helmholtz_pair.py --json out_helmholtz.json

# All benchmarks
python -m pytest validation_test/stream_function/benchmarks -v
```

## How to add a new benchmark

1. Subclass `Benchmark` from `benchmark_framework.py`.
2. Define `target_spec`, `solve()`, `metrics()`.
3. Add to the README table above.
4. Run, eyeball the numbers, commit `expected_*.json` for regression.

## Reporting format

Each benchmark writes a JSON file with at minimum:

```json
{
  "name": "...",
  "reference": "Author et al., Journal Vol (Year), DOI",
  "target_spec": {...},
  "our_result": {
    "rms": ..., "p2p_mean": ..., "wire_length_m": ...,
    "n_contours": ..., "inductance_H": ...,
    "elapsed_s": ...
  },
  "published_baseline": {...},
  "ratio_better_than": "field error ratio (us / them)"
}
```
