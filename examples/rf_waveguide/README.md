# RF waveguide examples

Analytic-gated RF and microwave examples for rectangular waveguide,
S-parameters, and cavity-style post-processing.

| Example | Shows | Capabilities used |
|---|---|---|
| [`validation_waveguide_vna_calibration.py`](validation_waveguide_vna_calibration.py) | Validation-class WR-90 offset-short calibration sweep: S11 phase -> group delay -> physical short offset, plus TE/TM impedance duality and dielectric-slab reflection metrics | `waveguide_offset_short_s11`, `sparameter_group_delay`, `waveguide_offset_short_length_from_group_delay`, `waveguide_wave_impedance`, `waveguide_dielectric_slab_sparams`, `reflection_metrics` |
| [`validation_waveguide_multisection_filter.py`](validation_waveguide_multisection_filter.py) | Validation-class TE10 multi-section dielectric filter: quarter-wave slab consistency, Bragg reflection, unitarity, and ABCD determinant | `waveguide_cascade_sparams`, `waveguide_dielectric_slab_sparams`, `reflection_metrics` |
| [`validation_waveguide_mode_table_sweep.py`](validation_waveguide_mode_table_sweep.py) | Validation-class rectangular-guide mode table: TE/TM cutoffs, single-mode band, multimode onset, guide wavelength, and below-cutoff attenuation | `rectangular_waveguide_mode_table`, `rectangular_waveguide_band_summary`, `waveguide_dispersion`, `waveguide_evanescent_attenuation` |
| [`validation_circular_waveguide_mode_sweep.py`](validation_circular_waveguide_mode_sweep.py) | Validation-class circular-guide mode table: Bessel-zero TE/TM cutoffs, degeneracy, single-mode band, and cutoff dispersion | `circular_waveguide_mode_table`, `circular_waveguide_band_summary`, `circular_waveguide_cutoff`, `waveguide_dispersion` |
| [`validation_tem_line_geometry_sweep.py`](validation_tem_line_geometry_sweep.py) | Validation-class TEM line geometry sweep for coax, two-wire, wire-plane, and microstrip quasi-static RF checks | `coaxial_line_parameters`, `two_wire_line_parameters`, `microstrip_line_parameters`, `tem_lc_identity_summary` |
| [`validation_waveguide_conductor_loss.py`](validation_waveguide_conductor_loss.py) | Validation-class WR-90 TE10 conductor-loss sweep: near-cutoff attenuation, length-linear insertion loss, conductivity scaling, and matched-line power balance | `rectangular_waveguide_te10_conductor_loss`, `rectangular_waveguide_cutoff` |
| [`validation_waveguide_te10_port_normalization.py`](validation_waveguide_te10_port_normalization.py) | Validation-class WR-90 TE10 port normalization: 1 W field amplitudes, Poynting power integral, longitudinal/transverse H ratio, and power scaling | `rectangular_waveguide_te10_port_normalization`, `waveguide_wave_impedance`, `waveguide_dispersion` |

```powershell
python validation_waveguide_vna_calibration.py
python validation_waveguide_multisection_filter.py
python validation_waveguide_mode_table_sweep.py
python validation_circular_waveguide_mode_sweep.py
python validation_tem_line_geometry_sweep.py
python validation_waveguide_conductor_loss.py
python validation_waveguide_te10_port_normalization.py
```

The examples are self-contained and use closed-form transmission-line /
waveguide formulas. Commercial RF solvers can be used internally as independent
references, but no solver-derived benchmark values are required to run or
validate these public examples.
