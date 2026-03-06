# TDK MnZn Ferrite Real Measurement Data

## Data Source

**Document**: TDK Corporation "Mn-Zn Ferrite Material characteristics" (May 2022)
**Filename**: `ferrite_mn-zn_material_characteristics_en.pdf`
**Publisher**: TDK Corporation
**Date**: May 2022 (Document revision: 20220510)

## Data Files

### Permeability Data (Extracted from Datasheet Graphs)

| File | Material | Source Page | Description |
|------|----------|-------------|-------------|
| `tdk_pc50_permeability.csv` | PC50 | Page 11 | Complex permeability vs frequency |
| `tdk_pc200_permeability.csv` | PC200 | Page 12 | Complex permeability vs frequency |

### Impedance Data (Converted)

| File | Material | Description |
|------|----------|-------------|
| `tdk_pc50_impedance.csv` | PC50 | Impedance for 10-turn toroidal inductor |
| `tdk_pc200_impedance.csv` | PC200 | Impedance for 10-turn toroidal inductor |

## Material Specifications

### PC50 (High-Frequency Power Supply Ferrite)

| Property | Value |
|----------|-------|
| Initial permeability (μi) | 1400 ± 25% |
| Saturation flux density (Bs) | 470 mT @ 25°C |
| Curie temperature (Tc) | > 240°C |
| Density | 4.8 × 10³ kg/m³ |
| Electrical resistivity | 30 Ω·m |

### PC200 (High-Frequency Power Supply Ferrite)

| Property | Value |
|----------|-------|
| Initial permeability (μi) | 800 ± 25% |
| Saturation flux density (Bs) | 485 mT @ 25°C |
| Curie temperature (Tc) | > 280°C |
| Density | 4.9 × 10³ kg/m³ |
| Electrical resistivity | 22 Ω·m |

## Test Conditions

Data extracted from "Magnetic permeability vs. frequency characteristics (Typ.)" graphs:

- **Temperature**: 23°C
- **Test signal**: Hm = 0.4 A/m
- **Frequency range**: 10 kHz - 10 MHz

## Conversion to Impedance

The permeability data was converted to impedance using the following parameters:

- **Core geometry**: Toroidal (OD=31mm, ID=19mm, TH=8mm) - standard TDK test core
- **Effective cross-section (Ae)**: 48 mm²
- **Effective path length (le)**: 78.54 mm
- **Number of turns (N)**: 10 (assumed for demonstration)
- **Air-core inductance (L0)**: 76.8 nH

Conversion formula:
```
Z = jωL = jωL₀(μ' - jμ")
  = ωL₀μ" + jωL₀μ'
  = R(f) + jX(f)
```

## Data Extraction Method

The complex permeability values (μ' and μ") were manually extracted from the
graphical data in the TDK datasheet PDF. The graphs show typical values for
each material grade.

**Note**: These are "Typ." (typical) values, not specification limits.
Individual cores may vary within the stated tolerance (±25% for μi).

## Usage for URN Validation

This data is used to validate the Universal Relaxation Network (URN) algorithm
with REAL manufacturer measurement data, addressing the "Lack of real-world
validation" critique in the peer review.

Run validation:
```bash
python validate_tdk_ferrite.py --material PC50
python validate_tdk_ferrite.py --material PC200
python validate_tdk_ferrite.py --material both
```

## License and Attribution

The original data is from TDK Corporation's publicly available datasheet.
This extraction is for academic research purposes.

**Citation**:
> TDK Corporation, "Mn-Zn Ferrite Material characteristics,"
> Document No. ferrite_mn-zn_material_characteristics_en, May 2022.
> Available: https://www.tdk-electronics.tdk.com/

## Files in This Directory

```
tdk_ferrite/
├── README.md                          # This file
├── tdk_pc50_permeability.csv          # PC50 μ'(f), μ"(f) data
├── tdk_pc50_impedance.csv             # PC50 Z(f) data (converted)
├── tdk_pc200_permeability.csv         # PC200 μ'(f), μ"(f) data
├── tdk_pc200_impedance.csv            # PC200 Z(f) data (converted)
├── convert_permeability_to_impedance.py  # Conversion script
└── ferrite_mn-zn_material_characteristics_en.pdf  # Original datasheet
```
