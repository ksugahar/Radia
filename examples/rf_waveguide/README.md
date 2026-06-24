# RF waveguide examples

Analytic-gated RF and microwave examples for rectangular waveguide,
S-parameters, and cavity-style post-processing.

| Example | Shows | Capabilities used |
|---|---|---|
| [`validation_waveguide_vna_calibration.py`](validation_waveguide_vna_calibration.py) | Validation-class WR-90 offset-short calibration sweep: S11 phase -> group delay -> physical short offset, plus TE/TM impedance duality and dielectric-slab reflection metrics | `waveguide_offset_short_s11`, `sparameter_group_delay`, `waveguide_offset_short_length_from_group_delay`, `waveguide_wave_impedance`, `waveguide_dielectric_slab_sparams`, `reflection_metrics` |

```powershell
python validation_waveguide_vna_calibration.py
```

The examples are self-contained and use closed-form transmission-line /
waveguide formulas. Commercial RF solvers can be used internally as independent
references, but no solver-derived benchmark values are required to run or
validate these public examples.
