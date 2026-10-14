# Universal Relaxation Network (URN)

This directory contains public documentation, result-bearing notebooks, and
reproduction drivers for the Universal Relaxation Network (URN), a KAN-inspired
approach for automatic discovery of physical relaxation mechanisms from
impedance data. The maintained implementation is the `radia.urn` package under
`src/radia/urn/`.

The canonical operating manual is the radia-ngsolve MCP `urn(topic="overview")`
tool; use `urn(topic="method")` for the model families. This page provides
discovery and historical evidence, not an independent workflow policy.

## Validation Results (2026-01-19)

### Comprehensive Real-World Data Performance

Historical results below are not renewed experimental acceptance. In particular,
the NASA frequency axis was model-assigned: its row and aggregate improvement
must not be cited as measured frequency-domain accuracy.

| Dataset | VF NRMSE | URN NRMSE | Improvement | URN Time |
|---------|----------|-----------|-------------|----------|
| NASA 18650 Battery | 0.2700 | **0.2454** | 9.1% | 178s |
| TDK PC47 Ferrite | 0.0146 | **0.0088** | 39.4% | 175s |
| TDK PC50 Ferrite | 0.0288 | **0.0098** | 65.9% | 176s |
| TDK PC95 Ferrite | **0.0080** | 0.0120 | -48.9% | 176s |
| TDK PC200 Ferrite | 0.0108 | **0.0056** | 48.4% | 177s |
| **Average** | --- | --- | **22.8%** | 176s |

**Key Findings**:
- **Overall**: URN outperforms Vector Fitting on 4/5 datasets (average 22.8% improvement)
- **Ferrite (PC47, PC50, PC200)**: URN achieves 39-66% lower error on materials with Cole-Cole relaxation
- **Ferrite (PC95)**: VF outperforms URN (-49%) on near-ideal Debye behavior
- **Legacy Attention Study**: older validation notebooks include an attention
  ablation; SA/RM work now uses the attention-free Y-domain dictionary.
- **Honest Assessment**: URN's advantage emerges when fractional-order dynamics dominate

## Directory Structure

```
universal_relaxation_network/
  data/
    real_world/                   # Real measurement datasets
      nasa_battery/               # NASA Li-ion Battery Aging Dataset
        README.md                # Acquisition/provenance; measurements kept private
      tdk_ferrite/                # TDK MnZn ferrite datasheet data
        tdk_pc50_impedance.csv    # PC50 impedance (included)
  demo_spice_timedomain.py        # Time-domain SPICE simulation demo
  cq_urn_bridge.py                # Passive URN H(s) -> BDF2 CQ teaching artifact
  cq_urn_bridge.ipynb             # Result-bearing CQ bridge notebook
```

Benchmark and validation JSON records are stored under
`validation_test/universal_relaxation_network/` together with their runnable
drivers; `docs/` retains public notebooks, narrative, immutable input data, and
embedded demonstration results.

## Quick Start

```python
from radia.urn import (
    UniversalRelaxationNetwork, URNConfig, train_urn, generate_spice_netlist
)
import numpy as np
import torch

# Load caller-owned data only after establishing the sample-aligned frequency axis.
# Set RADIA_NASA_EIS_CSV to its path outside the repository.
from generate_paper_figures import load_nasa_battery_data
freq, Z = load_nasa_battery_data()

# Configure and train the legacy URN path
config = URNConfig(n_debye=3, n_cole_cole=2, n_warburg=2, sparsity_weight=0.01)
model = train_urn(freq, Z, config)  # Returns trained model

# View discovered mechanisms
mechanisms = model.get_active_components()
for name, components in mechanisms.items():
    for comp in components:
        print(f"{name}: {comp}")

# Generate SPICE netlist (uses learned parameters)
netlist = generate_spice_netlist(model, "BATTERY")
with open("battery_model.sp", "w") as f:
    f.write(netlist)
print("SPICE netlist saved to battery_model.sp")
```

## SA/RM-2026 Y-Domain Variant

### Noise-aware basis selection

Include the richer basis candidate in comparisons, but do not increase basis
count merely to force fitting error below an arbitrary value such as 1e-3.
For noisy measurements, prefer a smaller active set when it explains the data
within independently estimated measurement uncertainty. Report the starting
basis count and the retained count separately, along with residuals, held-out
performance when available, and stability of relaxation times/weights under
repeated fits or noise perturbations. If noise is unknown, report sensitivity
rather than inventing a noise floor. Parsimony alone does not establish physical
truth; check passivity and consistency with the measurement conditions too.
The undocumented NASA frequency axis remains a separate blocker to interpreting
relaxation times, irrespective of the number of bases.

The current SA/RM-2026 research route starts from a single-layer, attention-free
**34-basis Y-domain dictionary**, not the earlier 22-basis route.
Its engineering objective is a compact passive response within a caller-owned
measurement error budget, rather than the smallest possible training residual.
Ten or twelve retained bases are comparison points, not lower bounds: even ten
bases can contain repeated response shapes. Basis count is not mechanism count.
Parameter fitting and basis-count selection are separate decisions. A better
fit alone does not justify extra bases that explain measurement error rather
than a reproducible response. A likelihood-based fit requires a declared noise
model; the current Huber fit/error-budget selection is not calibrated
maximum-likelihood inference. Parsimony remains an explicit selection objective.

The owning MCP `urn(topic="overview")` and `urn(topic="method")` describe the
workflow and acceptance contract. The implementation composition is:

```python
from radia.urn import (
    YAdmittanceURNConfig,
    reduce_y_admittance_urn,
    train_y_admittance_urn,
)

cfg = YAdmittanceURNConfig.research_34_basis()
model = train_y_admittance_urn(freq_hz, Z_measured, cfg)
compact, trace = reduce_y_admittance_urn(
    model, freq_hz, Z_measured, uncertainty_ohm=measurement_error_budget_ohm,
)
# compact is None if no tested candidate meets the supplied measurement budget.
```

The error budget must come from the caller's measurement assumptions; it is
not inferred as `1e-3`. The trace reports S-domain and log-component errors,
duplicate-response pairs and all candidate decisions, without raw data rows.
This greedy path selects the smallest tested acceptable model, not a globally
minimal circuit or an independently validated physical interpretation.

When checking a conference figure rather than the original measurement CSV,
treat digitized points as approximate. Use original measurement data for
quantitative claims; figure-extracted points are suitable only for qualitative
workflow checks.


See [`model_inventory.md`](model_inventory.md) for candidate models beyond the
original 22-basis dictionary, including parallel-RLC anti-resonance branches,
skin/proximity ladders, Havriliak-Negami relaxation, DRT diagnostics, and
passive rational macromodeling.

## Convolution Quadrature Bridge

`cq_urn_bridge.ipynb` records the compact path from an identified passive URN
relaxation model to a causal time-domain operator:

1. Fit a non-negative Debye ladder on a candidate relaxation-time grid.
2. Expose the fit as a Laplace-domain evaluator `H(s)`.
3. Generate BDF2 convolution-quadrature weights from `H(delta(zeta)/dt)`.
4. Compare the causal CQ response with a deliberately naive periodic IFFT
   contrast.

The builder writes `validation_test/universal_relaxation_network/cq_urn_bridge_results.json`,
while the notebook embeds the executed figure and summary.
This is the educational contract for later acoustic FEM/BEM and time-domain
Maxwell examples: replace only the scalar teaching `H(s)` with the solver's
passive boundary/material/operator response.

## Data Sources

### Bundled measured data

The repository includes extracted NASA 18650 EIS measurements and TDK PC47,
PC50, PC95, and PC200 impedance tables under `data/real_world/`. Their source
notes and conversion provenance are stored beside the CSV files. Synthetic
benchmark cases are generated by the validation drivers rather than maintained
as undocumented CSV fixtures.

### Additional real-world datasets

For broader validation, the following public datasets can be supplied to the
drivers under `validation_test/universal_relaxation_network/`:

#### NASA Li-ion Battery Aging Dataset
- **URL**: https://c3.ndc.nasa.gov/dashlink/resources/133/
- **Description**: Li-ion batteries run through charge, discharge, and EIS at different temperatures
- **Download**: http://ti.arc.nasa.gov/c/5/ (Dataset 1), http://ti.arc.nasa.gov/c/9/ (Dataset 2)
- **Format**: MATLAB .mat files
- **License**: Public Domain (NASA)

#### Mendeley SoC EIS Dataset (2024)
- **URL**: https://data.mendeley.com/datasets/cb887gkmxw/2
- **DOI**: 10.17632/cb887gkmxw.2
- **Description**: EIS measurements on 11 LiFePO4 batteries at 19 SoC levels
- **Format**: CSV files
- **License**: CC BY 4.0
- **Citation**: Mingant, R., Petit, M. (2024). SoC estimation on Li-ion batteries: A new EIS-based dataset for data-driven applications. Data in Brief, 56, 110807.

Pass downloaded NASA MAT or Mendeley CSV files explicitly with `--nasa-path`
or `--mendeley-path`. External downloads are not silently mixed with the
bundled validation inputs.

## Validation Drivers

The runnable validation and benchmark drivers are maintained under
`validation_test/universal_relaxation_network/`. Run them from the repository
root so their generated evidence remains in that validation directory.

### Real-World Data Validation
```bash
# With NASA dataset
python validation_test/universal_relaxation_network/validate_real_data.py --nasa-path data/real_world/nasa_battery/B0005.mat

# With Mendeley dataset
python validation_test/universal_relaxation_network/validate_real_data.py --mendeley-path data/real_world/mendeley_eis/cell1_soc50.csv

# Run all theoretical analyses (convergence, sensitivity, noise)
python validation_test/universal_relaxation_network/validate_real_data.py --bundled-nasa --all-tests
```

### Time-Domain Stability Verification
```bash
# Demonstrates URN vs Vector Fitting stability comparison
python validation_test/universal_relaxation_network/verify_timedomain_stability.py
```

### Ablation Study
```bash
# Feature contribution analysis
python validation_test/universal_relaxation_network/ablation_study.py --dataset battery --n-trials 3
```

### scikit-rf Vector Fitting Comparison
```bash
# Benchmark against industry-standard VF (requires scikit-rf)
pip install scikit-rf
python validation_test/universal_relaxation_network/benchmark_urn_vs_skrf_vf.py --dataset battery
```

### Actual SPICE Verification (LTspice)
```bash
# Run actual LTspice simulation (requires LTspice + PyLTSpice)
pip install PyLTSpice
python validation_test/universal_relaxation_network/run_ltspice_verification.py --dataset battery

# Specify custom LTspice path if needed
python validation_test/universal_relaxation_network/run_ltspice_verification.py --ltspice-path "C:/Program Files/ADI/LTspice/LTspice.exe"
```

## Related Paper

This implementation accompanied the following **draft**, not a verified
published or accepted IEEE Access article:

> K. Sugahara and Y. Sato, "KAN-inspired Universal Relaxation Network for Automatic Discovery of Physical Relaxation Mechanisms with Direct Circuit Synthesis," IEEE Access-format draft, 2026.

The located `paper/urn_paper.tex` has placeholder publication dates and DOI
`10.1109/ACCESS.2026.XXXXXXX`. Its adjacent `REVIEW.md` records an AI-agent
review, not publisher acceptance. Do not cite this as an accepted article.

The manuscript source is maintained **outside** the code repository in the
lab's conference-materials archive. The reproducibility scripts, result JSON,
and showcase notebooks in this directory remain in the repo.

## Requirements

```
numpy
scipy
torch
matplotlib
```

Optional for comparison and verification:
```
scikit-rf  # Industry-standard Vector Fitting (recommended)
PyLTSpice  # LTspice automation via Python
```

Optional system tools:
```
LTspice   # Analog Devices SPICE simulator (free, Windows/macOS)
```

## License

MIT License (same as main Radia project)
