# Stream-function benchmark suite

Reproducible benchmark cases for validating the Radia SF inverse-design
pipeline against published references and other open-source coil design
tools. Each benchmark loads a literature-defined target spec, runs the
full pipeline (SF solve → contour → single-stroke → field), and reports
the key metrics in a single JSON file (so the table can be regenerated
unattended).

## What's here vs what's TODO

| Status | File | Reference |
|--------|------|-----------|
| ✅ shipped | `bench_helmholtz_pair.py` | Analytical Helmholtz pair (uniform Bz, Maxwell) |
| 🚧 stub | `bench_bilac_planar_shim.py` | Bilac et al., MRI shim coil design |
| 🚧 stub | `bench_turner_cylinder_gz.py` | Turner, IEEE TMI 5 (1986), cylindrical Gz |
| 🚧 stub | `bench_lemdiasov_ludwig_2005.py` | Lemdiasov & Ludwig 2005, target-field method |
| 🚧 stub | `bench_coilgen_headtohead.py` | Schwartz et al. CoilGen (OSS comparison) |
| 🚧 stub | `bench_shielded_iron_yoke.py` | Material kernel via Radia MMM (material-kernel demo) |

The TODO stubs document the literature target spec + expected output and
fail explicitly with `NotImplementedError` so they're visible in the test
suite. Each is ~1-2 days of careful work to bring up; see the validation
strategy in MCP topic `aca_tsvd(session_2026_05_30)` section 11.

## How to run

```bash
# Single benchmark
python bench_helmholtz_pair.py --json out_helmholtz.json

# All benchmarks (skips stubs)
python -m pytest benchmarks/ -v
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
