# Universal Relaxation Network (URN) Examples

This directory contains the implementation and examples for the Universal Relaxation Network (URN), a KAN-inspired approach for automatic discovery of physical relaxation mechanisms from impedance data.

## Directory Structure

```
Universal_Relaxation_Network/
  data/
    synthetic/                    # Physics-based synthetic benchmark data
      liion_battery_eis.csv       # Synthetic Li-ion battery EIS
      mnzn_ferrite_impedance.csv  # Synthetic MnZn ferrite impedance
    real_world/                   # Publicly available real datasets
      nasa_battery/               # NASA Li-ion Battery Aging Dataset
      mendeley_eis/               # Mendeley SoC EIS Dataset
  universal_relaxation_network.py # Main URN implementation
  relaxation_basis_library.py     # Basis function library
  validate_urn_vs_vf.py           # URN vs Vector Fitting comparison
  demo_spice_timedomain.py        # Time-domain SPICE simulation demo
  urn_benchmark_suite.py          # Full benchmark suite
  validate_real_data.py           # Real-world data validation (NASA/Mendeley)
  verify_timedomain_stability.py  # URN vs VF time-domain stability test
  ablation_study.py               # Feature contribution analysis
  benchmark_urn_vs_skrf_vf.py     # scikit-rf Vector Fitting comparison
  run_ltspice_verification.py     # Actual LTspice simulation (PyLTSpice)
```

## Quick Start

```python
from universal_relaxation_network import (
    UniversalRelaxationNetwork, URNConfig, train_urn, generate_spice_netlist
)
import numpy as np

# Load impedance data
data = np.loadtxt('data/synthetic/liion_battery_eis.csv', delimiter=',', skiprows=7)
freq = data[:, 0]
Z = data[:, 1] + 1j * data[:, 2]

# Configure and train URN
config = URNConfig(n_debye=3, n_cole_cole=2, n_warburg=2, sparsity_weight=0.01)
model, history = train_urn(freq, Z, config)

# View discovered mechanisms
mechanisms = model.get_discovered_mechanisms()
for name, params in mechanisms.items():
    print(f"{name}: {params}")

# Generate SPICE netlist
netlist = generate_spice_netlist(model, "BATTERY")
print(netlist)
```

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
