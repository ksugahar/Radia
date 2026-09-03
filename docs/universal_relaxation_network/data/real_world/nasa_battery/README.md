# NASA 18650 Li-ion Battery EIS Data

## Data Source

**Dataset**: Li-ion Battery Aging Datasets
**Provider**: NASA Ames Prognostics Center of Excellence (PCoE)
**URL**: https://data.nasa.gov/dataset/li-ion-battery-aging-datasets

## Reference

```bibtex
@misc{nasa_battery,
    author = {Saha, B. and Goebel, K.},
    title = {{Battery Data Set}},
    year = {2007},
    publisher = {NASA Ames Prognostics Data Repository},
    institution = {NASA Ames Research Center, Moffett Field, CA},
    url = {https://ti.arc.nasa.gov/tech/dash/groups/pcoe/prognostic-data-repository/}
}
```

## Battery Specifications

| Property | Value |
|----------|-------|
| Form factor | 18650 cylindrical |
| Chemistry | Li-ion (LiCoO2 cathode) |
| Rated capacity | 2 Ah |
| Rated voltage | 3.7 V |
| Test temperature | 24 C |

## EIS Measurement Conditions

From NASA dataset documentation:
- **Frequency range**: 0.1 Hz - 5 kHz
- **Batteries tested**: B0005, B0006, B0007, B0018
- **Charge protocol**: CC-CV at 1.5 A to 4.2 V
- **Discharge**: Constant current at various levels

## Equivalent Circuit Model

The 18650 Li-ion battery impedance follows a modified Randles circuit:

```
Z(s) = R_s + R_ct / (1 + (j*omega*tau)^alpha) + sigma_w / sqrt(j*omega)
```

Parameters:
- R_s = 35 mOhm (series/ohmic resistance)
- R_ct = 25 mOhm (charge transfer resistance)
- tau = R_ct * C_dl = 50 us (RC time constant)
- alpha = 0.85 (Cole-Cole/CPE exponent)
- sigma_w = 15 Ohm*s^-0.5 (Warburg coefficient)

## Usage for URN Validation

This electrochemical impedance data demonstrates URN's ability to fit:
1. Cole-Cole relaxation (charge transfer)
2. Warburg diffusion (45-degree line)
3. Multi-element RC networks

```python
from radia.urn import UniversalRelaxationNetwork

# Load data
import pandas as pd
df = pd.read_csv('nasa_18650_eis.csv', comment='#')

# Create URN model
urn = UniversalRelaxationNetwork(n_rc_elements=8)
urn.fit(df['frequency_Hz'].values,
        df['Z_real_Ohm'].values + 1j * df['Z_imag_Ohm'].values)
```

## License

The NASA battery dataset is publicly available for research purposes.
See NASA's data policy at https://data.nasa.gov/
