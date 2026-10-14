# NASA 18650 Li-ion Battery EIS Data

## Acceptance status: historical hybrid, not frequency-domain validation

The 47 complex impedance values in the historical `nasa_18650_eis.csv` match NASA B0005's
first impedance operation (zero-based cycle 40), after dropping its first
sample. The frequency column does **not** come from that MAT file: it is a
47-point model-assigned logarithmic axis. The original README states sweep
endpoints only; it does not establish sample spacing or ordering.

Do not use that historical CSV to claim measured Bode accuracy, identified
relaxation times, or experimental time-domain/SPICE validation. The fitting
and SPICE loaders reject its undocumented axis. The CSVs are now retained
privately, not distributed in the current public tree; do not replace the axis
with another guessed sequence.

To enable measured fitting, obtain an independent, sample-aligned instrument
frequency record and regenerate using `extract_real_eis.py` with explicit
`frequency_hz` and `frequency_source`. Its provenance header is required by
both consumers. A synthetic example belongs in the separate `synthetic_*`
files and must remain labelled synthetic.

Archive identities and the exact complex-value comparison are recorded in
`docs/design/notebook_evidence_closeout_20260914.md` at the repository root.

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

## Historical fitting example (disabled until the axis is verified)

The model families of interest are:
1. Cole-Cole relaxation (charge transfer)
2. Warburg diffusion (45-degree line)
3. Multi-element RC networks

```python
from radia.urn import UniversalRelaxationNetwork

# Do not run this against the historical CSV. First regenerate it
# from an independently documented sample-frequency record as described above.
import pandas as pd
df = pd.read_csv('nasa_18650_eis.csv', comment='#')

# Create URN model
urn = UniversalRelaxationNetwork(n_rc_elements=8)
urn.fit(df['frequency_Hz'].values,
        df['Z_real_Ohm'].values + 1j * df['Z_imag_Ohm'].values)
```

## License

The official [NASA dataset listing](https://data.nasa.gov/dataset/groups/li-ion-battery-aging-datasets)
reports **License not specified** (checked 2026-09-14). Public download access
does not by itself establish redistribution permission. Do not add downloaded
MAT/ZIP files or newly extracted measurement CSVs to this public repository
without verifying the applicable permission. The owner chose private retention
on 2026-09-14; the three pre-existing NASA-named CSVs were backed up and removed
from the current tree. Git history was not rewritten. This README does not
grant a license to NASA data.

## Private input workflow

Obtain the data from the official provider under its applicable terms and keep
downloads and extracted CSVs outside the repository. Both consumers accept an
explicit `data_path` argument or the `RADIA_NASA_EIS_CSV` environment variable.
The latter also works when running either complete script. Local CSV/MAT/ZIP
files in this directory are ignored to prevent accidental re-addition.
Frequency provenance is still required; private storage does not make a guessed
axis valid. Automated loader tests use small self-generated values only.
