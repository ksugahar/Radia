# Universal Relaxation Network (URN) Examples

This directory contains the implementation and examples for the Universal Relaxation Network (URN), a KAN-inspired approach for automatic discovery of physical relaxation mechanisms from impedance data.

## Validation Results (2026-01-19)

### Comprehensive Real-World Data Performance

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
- **Attention Mechanism**: 79-83% accuracy improvement on real data (ablation study)
- **Honest Assessment**: URN's advantage emerges when fractional-order dynamics dominate

## Directory Structure

```
universal_relaxation_network/
  data/
    synthetic/                    # Physics-based synthetic benchmark data
      liion_battery_eis.csv       # Synthetic Li-ion battery EIS
      mnzn_ferrite_impedance.csv  # Synthetic MnZn ferrite impedance
    real_world/                   # Real measurement datasets
      nasa_battery/               # NASA Li-ion Battery Aging Dataset
        nasa_18650_eis.csv        # Extracted EIS data (included)
      tdk_ferrite/                # TDK MnZn ferrite datasheet data
        tdk_pc50_impedance.csv    # PC50 impedance (included)
  universal_relaxation_network.py # Main URN implementation (3900+ lines)
  relaxation_basis_library.py     # Basis function library
  validate_urn_vs_vf.py           # URN vs Vector Fitting comparison
  demo_spice_timedomain.py        # Time-domain SPICE simulation demo
  urn_benchmark_suite.py          # Full benchmark suite
  validate_real_data.py           # Real-world data validation (NASA/Mendeley)
  verify_timedomain_stability.py  # URN vs VF time-domain stability test
  ablation_study.py               # Feature contribution analysis
  benchmark_urn_vs_skrf_vf.py     # scikit-rf Vector Fitting comparison
  run_ltspice_verification.py     # Actual LTspice simulation (PyLTSpice)
  cq_urn_bridge.py                # Passive URN H(s) -> BDF2 CQ teaching artifact
  cq_urn_bridge.ipynb             # Result-bearing CQ bridge notebook
  cq_urn_bridge_results.json      # Machine-readable CQ bridge checks
```

## Quick Start

```python
from universal_relaxation_network import (
    UniversalRelaxationNetwork, URNConfig, train_urn, generate_spice_netlist
)
import numpy as np
import torch

# Load impedance data (use real NASA battery data if available)
data = np.loadtxt('data/real_world/nasa_battery/nasa_18650_eis.csv',
                  delimiter=',', skiprows=24)
freq = data[:, 0]
Z = data[:, 1] + 1j * data[:, 2]

# Configure and train URN (with attention enabled by default)
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

## Convolution Quadrature Bridge

`cq_urn_bridge.ipynb` records the compact path from an identified passive URN
relaxation model to a causal time-domain operator:

1. Fit a non-negative Debye ladder on a candidate relaxation-time grid.
2. Expose the fit as a Laplace-domain evaluator `H(s)`.
3. Generate BDF2 convolution-quadrature weights from `H(delta(zeta)/dt)`.
4. Compare the causal CQ response with a deliberately naive periodic IFFT
   contrast.

The notebook writes `cq_urn_bridge_results.json` and embeds the executed figure.
This is the educational contract for later acoustic FEM/BEM and time-domain
Maxwell examples: replace only the scalar teaching `H(s)` with the solver's
passive boundary/material/operator response.

## Data Sources

### Synthetic Benchmark Data (included)

The `data/synthetic/` directory contains physics-based synthetic impedance data:

- **liion_battery_eis.csv**: Synthetic Li-ion battery EIS modeled using standard Randles circuit parameters from literature (Barsoukov & Macdonald, "Impedance Spectroscopy", 2018). Represents typical 18650 cell at 50% SOC.

- **mnzn_ferrite_impedance.csv**: Synthetic MnZn ferrite impedance modeled using Cole-Cole relaxation with parameters typical of power ferrites (Snelling, "Soft Ferrites", 1988).

**IMPORTANT**: These are synthetic data for algorithm validation, not real measurements.

### Real-World Datasets (download required)

For validation on real measured data, we recommend the following publicly available datasets:

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

To use real-world data:
1. Download the dataset from the links above
2. Place CSV/MAT files in `data/real_world/nasa_battery/` or `data/real_world/mendeley_eis/`
3. Run the validation scripts with `--real-data` flag

## Validation Scripts

### Real-World Data Validation
```bash
# With NASA dataset
python validate_real_data.py --nasa-path data/real_world/nasa_battery/B0005.mat

# With Mendeley dataset
python validate_real_data.py --mendeley-path data/real_world/mendeley_eis/cell1_soc50.csv

# Run all theoretical analyses (convergence, sensitivity, noise)
python validate_real_data.py --synthetic --all-tests
```

### Time-Domain Stability Verification
```bash
# Demonstrates URN vs Vector Fitting stability comparison
python verify_timedomain_stability.py
```

### Ablation Study
```bash
# Feature contribution analysis
python ablation_study.py --dataset battery --n-trials 3
```

### scikit-rf Vector Fitting Comparison
```bash
# Benchmark against industry-standard VF (requires scikit-rf)
pip install scikit-rf
python benchmark_urn_vs_skrf_vf.py --dataset battery
```

### Actual SPICE Verification (LTspice)
```bash
# Run actual LTspice simulation (requires LTspice + PyLTSpice)
pip install PyLTSpice
python run_ltspice_verification.py --dataset battery

# Specify custom LTspice path if needed
python run_ltspice_verification.py --ltspice-path "C:/Program Files/ADI/LTspice/LTspice.exe"
```

## Related Paper

This implementation accompanies the paper:

> K. Sugahara and Y. Sato, "KAN-inspired Universal Relaxation Network for Automatic Discovery of Physical Relaxation Mechanisms with Direct Circuit Synthesis," IEEE Access, 2026.

The paper is located at `docs/paper/urn_paper.tex`.

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
