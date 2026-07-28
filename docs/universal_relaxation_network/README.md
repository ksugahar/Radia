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
- **Legacy Attention Study**: older validation notebooks include an attention
  ablation; SA/RM work now focuses on attention-free CLN peeling.
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

The SA/RM-2026 research-meeting manuscript uses a newer attention-free Y-base
formulation: 22 physical basis functions are
summed as a parallel admittance network, fitted through an S-domain Huber loss,
and selected by output ablation.  The Radia package exposes this variant
without replacing the older Z-domain implementation:

```python
from radia.urn import (
    YAdmittanceURNConfig,
    refit_y_admittance_active_bases,
    s_domain_rmse,
    train_y_admittance_urn,
)

cfg = YAdmittanceURNConfig.paper_22_basis()
model = train_y_admittance_urn(freq_hz, Z_measured, cfg)
active = model.active_bases()  # output-ablation ranking
Z_fit = model.predict(freq_hz)
rmse_s = s_domain_rmse(Z_fit, Z_measured)

realizable = refit_y_admittance_active_bases(freq_hz, Z_measured, active, cfg)
Z_realizable = realizable.predict(freq_hz)
```

Use this API when reproducing the SA-26/RM-26 draft model
(Debye/magnetic-Debye/Cole-Cole/magnetic-Cole-Cole/inductive-CPE/
capacitive-CPE/series-RLC in Y space).  Use `train_urn` for the original
Radia URN model used by the NASA/TDK validation notebooks.

When checking a conference figure rather than the original measurement CSV,
treat digitized points as approximate. Use original measurement data for
quantitative claims; figure-extracted points are suitable only for qualitative
workflow checks.

For time-domain review, this attention-free Y-branch is now treated mainly as a
branch model for CLN peeling.  A single parallel 22-basis sum is not rich enough
to claim VF-level accuracy on the extracted SA/RM curves; the preferred path is
to use it inside a continued-fraction residual peeling topology.

See [`model_inventory.md`](model_inventory.md) for candidate models beyond the
current 22-basis dictionary, including parallel-RLC anti-resonance branches,
skin/proximity ladders, Havriliak-Negami relaxation, DRT diagnostics, and
passive rational macromodeling.

## Cauer-Ladder Direction

For SA/RM review work, Radia also includes an experimental differentiable
Cauer-ladder candidate:

```python
from radia.urn import (
    CauerLadderURNConfig,
    train_cauer_ladder_alternating,
    train_cauer_ladder_progressive,
)

cfg = CauerLadderURNConfig.twenty_two_parameter_candidate()
model = train_cauer_ladder_alternating(freq_hz, Z_measured, cfg)
Z_fit = model.predict(freq_hz)
sections = model.parameter_summary()
```

The ladder evaluates a positive continued fraction,
`Z_k = R_k + s L_k + 1/(G_k + s C_k + 1/Z_{k+1})`, using PyTorch autograd.
A six-section ladder has 24 positive parameters, close to the 22-basis
Y-domain dictionary, but it can represent pole-zero/anti-resonance behavior
through series/parallel nesting rather than by adding many parallel basis
functions.

`train_cauer_ladder_alternating` alternates direct impedance-domain updates for
the series elements (`R,L`) with direct admittance-domain updates for the shunt elements
(`G,C`).  Blocks that make the combined Z/Y objective unstable are rolled back
and retried with a lower learning rate.  This follows the review idea that some
parameters are better identified in impedance form while others are better
identified in admittance form.

There is also an experimental `use_peeling_initialization=True` mode that tries
to initialize the ladder by alternately peeling series impedance and shunt
admittance terms.  Early checks on the PDF-extracted SA/RM curves show that
naive peeling is not yet reliable; direct alternating optimization is the safer
baseline for now.

For harder resonance/anti-resonance data, use the pole-zero assisted path:

```python
from radia.urn import (
    CauerLadderURNConfig,
    fit_rational_pole_zero,
    train_cauer_ladder_tail_then_polish,
)

cfg = CauerLadderURNConfig.twenty_two_parameter_candidate(
    use_rational_initialization=True,
    use_least_squares_polish=True,
    frozen_outer_sections=2,
)
teacher = fit_rational_pole_zero(freq_hz, Z_measured, order=6)
model = train_cauer_ladder_tail_then_polish(freq_hz, Z_measured, cfg)
```

This first builds a small pole-zero rational teacher, distills that response
into positive Cauer parameters by nonlinear least squares, trains the inner
ladder with the outer sections frozen, and finally polishes all sections.

## CLN Peeling Direction

`train_cln_peeling_urn` implements paired-basis CLN peeling.  Stage ``n``
represents the current driving-point impedance as an even series branch plus an
odd shunt branch loaded by the next tail:

```text
R_n = Z_2n + (Z_2n+1 || R_n+1)
```

One 22-basis composite fit supplies the physical basis shapes of the pair; the
fitted coefficients are split continuously between the even and odd branches
(soft split ``a_2n,k = a_k p_k``, ``a_2n+1,k = a_k (1 - p_k)``) while a fresh
22-basis lookahead model represents ``R_n+1``.  Accepted pairs are frozen and
the measured tail is peeled by the exact inverse map
``R_n+1 = [1/(R_n - Z_2n) - 1/Z_2n+1]^-1``; earlier stages are never
re-trained (no global polish).

```python
from radia.urn import CLNPeelingConfig, train_cln_peeling_urn

cfg = CLNPeelingConfig(n_stages=2)
model = train_cln_peeling_urn(freq_hz, Z_measured, cfg)
Z_lookahead = model.predict_terminated(freq_hz, termination="lookahead")
audit = model.audit_passivity()  # dense-grid positive-real audit
```

Evaluation policy: report the ``termination="lookahead"`` S-domain RMSE (the
learned physical tail) together with the dense-grid ``audit_passivity()``
report.  ``termination="stored"`` re-inserts the exactly peeled measurement
residue, reconstructs the training grid to machine precision by construction,
and is only defined on the training grid -- it is an identity check and must
never be quoted as fit accuracy.  All terminations except ``stored`` accept an
arbitrary frequency grid, which is what the audit uses for interpolation and
extrapolation checks.

### Stage-wise trust region (2026-07-27)

The exact peel amplifies measurement error wherever ``R_n - Z_2n`` cancels,
and it produces sign-unstable spikes wherever the peeled tail admittance
``1/(R_n - Z_2n) - 1/Z_2n+1`` nearly vanishes: in such parallel-resonance
bands the tail barely loads the ladder, so its peeled value is amplified noise.
On the SA/RM PCB coil the 1.36--1.49 MHz self-resonance band produced a
-59.7 kOhm negative-real spike in ``R_1`` (median ``|R_1|`` is 56 Ohm), which
no passive dictionary can represent; the second stage was rejected with
``min_parallel_real_normalized = -1063`` no matter what the series branch did.

Each stage therefore assigns a per-frequency trust weight to its peeled tail,
built from the relative series cancellation and the peeled-tail admittance
magnitude (config fields ``denominator_margin_relative`` and
``tail_admittance_margin_relative``).  The next stage fits its composite seed
and pair split with those weights, and its acceptance gates use trusted points
only (``min_parallel_real_trusted``, ``min_tail_admittance_real_trusted``,
``seed/parent_s_rmse_trusted``) plus a ``trusted_fraction >=
min_trusted_fraction`` gate; inherited weights multiply stage by stage.  A
``denominator_margin_relative`` hinge in the pair loss additionally discourages
the current stage from manufacturing new cancellation bands.  The stored exact
tail is never modified, so identity reconstruction is untouched.

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
