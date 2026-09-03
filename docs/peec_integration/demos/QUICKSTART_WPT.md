# PEEC for Wireless Power Transfer (WPT) - Quick Start

## Overview

This guide demonstrates using Radia PEEC for WPT coil design at 85 kHz (SAE J2954 standard).

## Prerequisites

```bash
pip install radia numpy matplotlib
```

## 1. Basic Coil Impedance Analysis

```python
import numpy as np

# Import PEEC builder
from radia.peec_matrices import PEECBuilder

# Physical parameters
LENGTH = 0.1       # 100mm coil segment
WIDTH = 2e-3       # 2mm wire width
HEIGHT = 0.5e-3    # 0.5mm wire height
SIGMA = 5.8e7      # Copper conductivity

# Create PEEC model
builder = PEECBuilder()
builder.create_wire([0, 0, 0], [LENGTH, 0, 0], WIDTH, HEIGHT, 20, SIGMA)
L, R, P, M_LS = builder.build()

# Analyze at 85 kHz
f = 85000
omega = 2 * np.pi * f
Z_matrix = np.diag(R) + 1j * omega * L
Z_total = np.sum(Z_matrix)

print(f"Frequency: {f/1000} kHz")
print(f"Impedance: {abs(Z_total)*1000:.3f} mOhm")
print(f"Phase: {np.angle(Z_total, deg=True):.1f} deg")
```

## 2. Generate SPICE Model

```bash
# Generate SPICE netlist
python spice/demo_prima_spice_export.py

# Output files:
#   wire_full.sp  - Full PEEC model
#   wire_prima.sp - PRIMA-reduced model (recommended)
```

## 3. Simulate in LTspice/ngspice

```spice
* Include generated subcircuit
.include wire_prima.sp

* Test circuit
Xcoil port_in port_out WIRE_PRIMA
Vin port_in 0 AC 1

* AC analysis
.AC DEC 100 1k 10MEG
.END
```

## Key Parameters for WPT Design

| Parameter | Typical Value | Notes |
|-----------|---------------|-------|
| Frequency | 85 kHz | SAE J2954 standard |
| Coil Q | 100-300 | Higher is better |
| Coupling k | 0.1-0.3 | Depends on gap |
| Efficiency target | >90% | Requires k*Q > 10 |

## Validation Results

Tested against analytical formulas:

| Metric | PEEC | Analytical | Error |
|--------|------|------------|-------|
| L | 89.65 nH | 97.41 nH | 8% |
| R_DC | 1.724 mOhm | 1.724 mOhm | 0% |
| Phase@85kHz | 84.4° | - | - |

## Next Steps

1. **Multi-coil systems**: Use `CplMag` for TX/RX coil coupling
2. **Magnetic core**: Add ferrite using `MatLin(mu_r)`
3. **Shielding**: Add aluminum shield with PEEC conductor model

## References

- `verification/verify_frequency_response_85khz.py` - Full validation script
- `wpt/demo_wpt.py` - WPT system demo with multiple coils
- [PEEC_SURFACE_IMPEDANCE.md](../../docs/peec/PEEC_SURFACE_IMPEDANCE.md) - Theory documentation
