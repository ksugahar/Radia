# PEEC Integration Examples

PEEC (Partial Element Equivalent Circuit) method implementation and verification,
combined with PRIMA model order reduction and Dowell continued-fraction expansion.

> **Note**: We use "PRIMA" (Passive Reduced-order Interconnect Macromodeling Algorithm)
> instead of "CLN" (Cauer Ladder Network). Both are mathematically equivalent
> (Lanczos tridiagonalization for ladder network representation), but PRIMA is
> the widely recognized term from the 1998 IEEE paper and avoids patent issues.

**Theory**: See [docs/PEEC_SURFACE_IMPEDANCE.md](../../docs/peec/PEEC_SURFACE_IMPEDANCE.md).

## Directory Structure

```
peec_integration/
├── ngbem/                   # ngsbem (NGSolve BEM) based PEEC demos
├── ngsbem_peec_demo/        # ngsbem PEEC Jupyter notebook
├── basic_peec/              # Basic PEEC demos (loop, DC, mesh input)
├── coil/                    # Coil impedance analysis
├── coupled/                 # Coupled analysis (mutual inductance, eddy currents)
├── wpt/                     # Wireless power transfer (WPT)
├── applications/            # Applications (dielectric, magnetic core, induction heating, NMR)
├── spice/                   # SPICE / Verilog-A export
├── benchmarks/              # Benchmarks and comparisons
├── verification/            # Verification, validation, and debug scripts
├── analysis/                # Analysis and paper figure generation
├── gmsh_models/             # Historical GMSH fixtures and .vol migration demos
├── data/                    # Measurement data (CSV)
└── model/                   # Model cache
```

Each folder contains its own `README.md`.
Dowell continued-fraction derivation scripts have been promoted to
`docs/peec_integration/algorithm_development/` with result-bearing notebook
archives preserving the original examples source hashes.

### Demo Folders

| Folder | Description |
|--------|-------------|
| `ngbem/` | ngsbem PEEC Loop-Star / FEM-BEM eddy current solver |
| `ngsbem_peec_demo/` | **Jupyter notebook**: ngsbem PEEC impedance extraction |
| `basic_peec/` | Core PEEC: loop analysis, DC, mesh input, Loop-Star |
| `coil/` | Coil impedance, coil on magnetic core |
| `coupled/` | Coupled inductance, eddy current inductance reduction |
| `wpt/` | WPT 85 kHz coils, shielded models |
| `applications/` | Dielectric, FastHenry comparison, induction heating, NMR magnet |
| `spice/` | PRIMA reduced-order model SPICE/Verilog-A export |

## Quick Start

```bash
cd examples/peec_integration

# PRIMA + Dowell correction verification
python spice/prima_with_dowell_correction.py

# PEEC + PRIMA integration demo
python spice/demo_peec_prima_reduction.py

# ngsbem PEEC notebook
jupyter notebook ngsbem_peec_demo/peec_impedance.ipynb
```

## Key Verification Results

### 1. PRIMA(DC) + Dowell correction = Dowell formula (exact match)

Verified by `spice/prima_with_dowell_correction.py`:

| Frequency | |Z_Dowell| | |Z_PRIMA+Dowell| | Error |
|-----------|-----------|----------------|-------|
| 1 Hz | 1.724e-04 | 1.724e-04 | 0% |
| 1 kHz | 1.724e-04 | 1.724e-04 | 0% |
| 100 kHz | 2.044e-04 | 2.044e-04 | 0% |
| 10 MHz | 1.229e-03 | 1.229e-03 | 0% |
| 100 MHz | 3.867e-03 | 3.867e-03 | 0% |

Maximum error < 1e-10% (numerical precision limit).

### 2. Continued-fraction expansion of z*coth(z)

```
z*coth(z) = 1 + w/(3 + w/(5 + w/(7 + w/(9 + ...))))
w = z^2 = tau * s,   tau = d^2 * mu * sigma / 2
```

This expansion is **exact** (no approximation error).

### 3. Ladder network conversion

The continued fraction maps to an RL ladder circuit:

```
        R1        R2        R3
  o----[  ]------[  ]------[  ]------...
         |         |         |
        L1        L2        L3
         |         |         |
        ===       ===       ===
```

## Theory

### Dowell formula

Impedance with skin-effect correction:

```
Z(s) = R_dc * F_R(xi) + s * L_int_dc * F_L(xi)
```

where `R_dc = 1/(sigma*d)`, `L_int_dc = mu*d/3`, `xi = d/delta`,
`delta = sqrt(2/(omega*mu*sigma))`.

### Extension with magnetic core

```
Z_total(s) = Z_cond(s) + Z_mag(s)

Z_cond(s) = R_dc * F_R + s * L_int_dc * F_L + s * L_ext_air   [conductor: Dowell applies]
Z_mag(s)  = R_mag(omega) + s * L_mag                           [magnetic: Dowell does not apply]
```

See `coil/coil_on_magnetic_core_peec.py` for implementation.

## SPICE / Verilog-A Export

```bash
python spice/demo_prima_spice_export.py              # SPICE netlist
python spice/demo_prima_spice_export.py --verilog-a   # Verilog-A
python spice/demo_prima_spice_export.py --lanczos 10  # Set Lanczos order
```

| Format | Extension | Compatibility | Use case |
|--------|-----------|---------------|----------|
| SPICE netlist | `.sp` | High (LTspice, PSpice, ngspice, etc.) | Circuit simulation |
| Verilog-A | `.va` | Limited (Spectre, ADS, HSPICE) | Advanced modeling |

## Paper Figure Generation

```bash
python analysis/generate_paper_figures.py
```

| Figure | Content | Output |
|--------|---------|--------|
| fig1 | PEEC matrix verification | `analysis/figures/fig1_peec_matrix_verification.pdf` |
| fig2 | PRIMA accuracy comparison | `analysis/figures/fig2_prima_accuracy.pdf` |
| fig3 | Adaptive MSC error | `analysis/figures/fig3_adaptive_msc_error.pdf` |
| fig4 | WPT coil (85 kHz) | `analysis/figures/fig4_wpt_coil_85khz.pdf` |

## References

1. P.L. Dowell, "Effects of eddy currents in transformer windings," Proc. IEE, 1966.
2. J.A. Ferreira, "Electromagnetic Modelling of Power Electronic Converters," 1989.
3. A. Ruehli, "Equivalent Circuit Models for Three-Dimensional Multiconductor Systems," IEEE Trans. MTT, 1974.
4. A. Odabasioglu et al., "PRIMA: Passive Reduced-Order Interconnect Macromodeling Algorithm," IEEE Trans. CAD, 1998.
