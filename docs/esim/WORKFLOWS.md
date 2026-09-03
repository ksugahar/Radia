# ESIM Workflows

Catalog of the public notebook and checked validation workflows for ESIM.
The notebook is the human-facing demonstration; scripts and JSON evidence live
under `validation_test/` and are reproducible from a clean
`pip install radia[cubit]` checkout.

## Cell-problem level

| Script | What it does | Validation tier |
|---|---|---|
| [`validation_test/ih_esim_benchmark/analytical_bessel_baseline.py`](../../validation_test/ih_esim_benchmark/analytical_bessel_baseline.py) | Compares the 1-D cell solver against the closed-form Bessel $I_0/I_1$ reference for a Cu cylinder with imposed $\mu_r{=}100$, $\xi\!\in\![4,140]$.  Publication-ready Table I. | **(i)** Bessel linear-$\mu$ |
| [`validation_test/ih_esim_benchmark/plot_cell_envelope.py`](../../validation_test/ih_esim_benchmark/plot_cell_envelope.py) | Sweeps the cell solver across $|H_t|\!\in\![1,10^5]$ A/m at 50 kHz on the IH benchmark steel; overlays low-H Bessel and high-H thin-skin asymptotes; highlights the IH headline $|H_t|$ band.  Publication-ready Fig. 1. | **(ii)** Nonlinear envelope |

## End-to-end IH workpiece

| Script | What it does | Notes |
|---|---|---|
| [`docs/induction_heating/induction_heating_demo_showcase.ipynb`](../induction_heating/induction_heating_demo_showcase.ipynb) | Executed ESIM/Bessel public showcase with embedded figures; checked numerical JSON remains in `validation_test/ih_esim_benchmark/`. | Tutorial entry point |
| [`validation_test/ih_esim_benchmark/benchmark.py`](../../validation_test/ih_esim_benchmark/benchmark.py) | Drives the 3-path (PEEC-BEM / FEM-Kelvin / FEM-coilmesh) Karl benchmark at 10/50/100/500 kHz.  Emits `results.json` and `benchmark_plot.png`. | Validation tier **(iv)** consistency table |
| [`validation_test/ih_esim_benchmark/plot_zs_per_dof_map.py`](../../validation_test/ih_esim_benchmark/plot_zs_per_dof_map.py) | 3-panel side-wall map of per-DOF $\mathrm{Re}\,Z_s$, $\mathrm{Im}\,Z_s$, $|H_t|$ from a `--esim-per-panel` JSON.  Publication-ready Fig. 2 (hotspot pattern). | Requires `esim_per_panel_H_t` array (radia $\geq$ 4.55.x) |
| [`validation_test/ih_esim_benchmark/plot_karl_history.py`](../../validation_test/ih_esim_benchmark/plot_karl_history.py) | Karl-iteration convergence diagnostic.  3-panel plot of $dZ$ (log), $\|Z_s\|$ with per-DOF min/max band, and $\|H_t\|$ per iteration.  Accepts both scalar- and per-panel-Karl JSONs.  Use to distinguish *convergence*, *plateau-at-max-iter* and *divergence* (see [`IMPLEMENTATION.md`](IMPLEMENTATION.md) § 3.4 decision table). | Diagnostic / triage |
| [`validation_test/induction_heating/bem_reference/`](../../validation_test/induction_heating/bem_reference/) + `radia.bem_inductance` / `radia.ngsbem_*` | Reference BEM scripts and reusable BEM APIs retained for research cross-check. | Research / validation |

## Tests

| Script | What it does |
|---|---|
| [`tests/test_esim_integration.py`](../../tests/test_esim_integration.py) | Pytest suite covering cell-problem solver, ESI table interface, and coupled solver invocation. |
| [`validation_test/panels/golden/`](../../validation_test/panels/golden/) | Checked JSON evidence for `calc_inductance.py` ESIM mode (Cu and steel benchmarks at fixed frequency / current). |

## Reproducer for the IGTE 2026 paper headline

The 48% per-element-vs-scalar gap reported in the IGTE 2026 paper
is reproduced by these three commands in sequence:

```bash
# 1) Cell-solver linear-mu Bessel validation (Table I)
python validation_test/ih_esim_benchmark/analytical_bessel_baseline.py

# 2) Per-element + scalar Karl runs at 50 kHz / 100 A
python src/radia/panels/calc_inductance.py \
    --coil-step src/radia/panels/samples/ih_fem_kelvin_demo_coil.step \
    --coil-solver peec \
    --vol src/radia/panels/samples/ih_bem_sample_p1.vol --wp-label sibc \
    --sigma 2e6 --mu-r 100 --half-thickness 0.005 \
    --frequency 50000 --current 100.0 --coil-sigma 5.8e7 \
    --impedance-model esim --bh-file src/radia/panels/samples/em_sample_bh.txt \
    --esim-max-iter 15 --esim-tol 1e-3 --esim-relax 0.5 \
    --esim-per-panel \
    --h1-order 1 --wp-bem-backend intree-dense \
    --output C:/temp/I100_per_panel.json
# (drop --esim-per-panel for the scalar comparator)

# 3) Visualisation (Fig. 2 of the paper)
python validation_test/ih_esim_benchmark/plot_zs_per_dof_map.py \
    C:/temp/I100_per_panel.json
```

For the full numerical results JSON used in the paper see
[`validation_test/ih_esim_benchmark/results.json`](../../validation_test/ih_esim_benchmark/results.json).

## Cross-references

- [`README.md`](README.md) — top-level index of `docs/esim/`.
- [`USAGE.md`](USAGE.md) — production CLI flag reference.
- [`SCALAR_BIE_VS_VECTOR_BEM.md`](SCALAR_BIE_VS_VECTOR_BEM.md) — the methodological argument tying the scripts above into a coherent paper-grade workflow.

---

**Document version**: 2026-09-03 (docs/validation ownership aligned).
