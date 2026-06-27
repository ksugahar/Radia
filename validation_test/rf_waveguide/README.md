# RF waveguide validation corpus

Analytic-gated RF and microwave validation scripts for rectangular waveguide,
S-parameters, and cavity-style post-processing.

This directory is the executable validation surface. Each script refreshes its
adjacent summary JSON with `schema`, `generated_at_utc`, and runtime version
metadata, so no separate docs-layer source archive is required.

| Example | Shows | Capabilities used |
|---|---|---|
| [`validation_waveguide_vna_calibration.py`](validation_waveguide_vna_calibration.py) | Validation-class WR-90 offset-short calibration sweep: S11 phase -> group delay -> physical short offset, plus TE/TM impedance duality and dielectric-slab reflection metrics | `waveguide_offset_short_s11`, `sparameter_group_delay`, `waveguide_offset_short_length_from_group_delay`, `waveguide_wave_impedance`, `waveguide_dielectric_slab_sparams`, `reflection_metrics` |
| [`validation_waveguide_multisection_filter.py`](validation_waveguide_multisection_filter.py) | Validation-class TE10 multi-section dielectric filter: quarter-wave slab consistency, Bragg reflection, unitarity, and ABCD determinant | `waveguide_cascade_sparams`, `waveguide_dielectric_slab_sparams`, `reflection_metrics` |
| [`validation_waveguide_mode_table_sweep.py`](validation_waveguide_mode_table_sweep.py) | Validation-class rectangular-guide mode table: TE/TM cutoffs, single-mode band, multimode onset, guide wavelength, and below-cutoff attenuation | `rectangular_waveguide_mode_table`, `rectangular_waveguide_band_summary`, `waveguide_dispersion`, `waveguide_evanescent_attenuation` |
| [`validation_circular_waveguide_mode_sweep.py`](validation_circular_waveguide_mode_sweep.py) | Validation-class circular-guide mode table: Bessel-zero TE/TM cutoffs, degeneracy, single-mode band, and cutoff dispersion | `circular_waveguide_mode_table`, `circular_waveguide_band_summary`, `circular_waveguide_cutoff`, `waveguide_dispersion` |
| [`validation_tem_line_geometry_sweep.py`](validation_tem_line_geometry_sweep.py) | Validation-class TEM line geometry sweep for coax, two-wire, wire-plane, and microstrip quasi-static RF checks | `coaxial_line_parameters`, `two_wire_line_parameters`, `microstrip_line_parameters`, `tem_lc_identity_summary` |
| [`validation_waveguide_conductor_loss.py`](validation_waveguide_conductor_loss.py) | Validation-class WR-90 TE10 conductor-loss sweep: near-cutoff attenuation, length-linear insertion loss, conductivity scaling, and matched-line power balance | `rectangular_waveguide_te10_conductor_loss`, `rectangular_waveguide_cutoff` |
| [`validation_waveguide_te10_port_normalization.py`](validation_waveguide_te10_port_normalization.py) | Validation-class WR-90 TE10 port normalization: 1 W field amplitudes, Poynting power integral, longitudinal/transverse H ratio, and power scaling | `rectangular_waveguide_te10_port_normalization`, `waveguide_wave_impedance`, `waveguide_dispersion` |
| [`validation_radiation_pressure_poynting.py`](validation_radiation_pressure_poynting.py) | Validation-class RF radiation pressure: TE10 Poynting power to absorber/reflector force via momentum flux | `radiation_pressure_summary`, `radiation_force_from_power`, `plane_wave_intensity_from_electric_field`, `rectangular_waveguide_te10_port_normalization` |
| [`validation_oblique_radiation_pressure.py`](validation_oblique_radiation_pressure.py) | Validation-class oblique-incidence radiation pressure: normal `cos^2(theta)` force plus absorbed tangential momentum | `oblique_radiation_pressure_summary`, `radiation_pressure_summary` |
| [`validation_poynting_patch_force_vector.py`](validation_poynting_patch_force_vector.py) | Validation-class vector radiation force from a 3D Poynting vector, checked against oblique-incidence pressure components | `poynting_patch_force_summary`, `oblique_radiation_pressure_summary` |
| [`validation_time_harmonic_maxwell_stress.py`](validation_time_harmonic_maxwell_stress.py) | Validation-class complex phasor Maxwell stress: peak/RMS conventions, plane-wave momentum flux, and local traction sign | `time_average_maxwell_stress_tensor`, `time_average_maxwell_traction_summary`, `radiation_pressure_from_intensity` |
| [`validation_scattering_radiation_force.py`](validation_scattering_radiation_force.py) | Validation-class normal-incidence scattering force from reflectance/transmittance, matching `F=(1+R-T)P/c=(A+2R)P/c` | `radiation_scattering_force_summary`, `radiation_force_from_normal_scattering` |
| [`validation_two_port_scattering_momentum_force.py`](validation_two_port_scattering_momentum_force.py) | Validation-class vector momentum balance for a two-port scatterer: straight through-line, short, absorber, lossy line, and 90-degree bend | `two_port_scattering_momentum_force_summary` |
| [`validation_one_port_reflection_momentum_force.py`](validation_one_port_reflection_momentum_force.py) | Validation-class one-port S11 momentum force: matched load, perfect short, and phase-independent partial reflection | `one_port_reflection_momentum_force_summary` |
| [`validation_one_port_reflection_sweep_force.py`](validation_one_port_reflection_sweep_force.py) | Validation-class one-port S11 sweep force audit: max/min force frequencies, mean force, and passivity flags | `one_port_reflection_sweep_momentum_force_summary` |
| [`validation_two_port_scattering_sweep_force.py`](validation_two_port_scattering_sweep_force.py) | Validation-class two-port S11/S21 sweep force audit: extrema, mean force, absorbed fraction, and passivity flags | `two_port_scattering_sweep_momentum_force_summary` |
| [`validation_two_port_sparameter_health.py`](validation_two_port_sparameter_health.py) | Validation-class two-port S-parameter health: passivity, reciprocity, return symmetry, and momentum-force extrema | `two_port_sparameter_sweep_health_summary` |

```powershell
python validation_test/rf_waveguide/validation_waveguide_vna_calibration.py
python validation_test/rf_waveguide/validation_waveguide_multisection_filter.py
python validation_test/rf_waveguide/validation_waveguide_mode_table_sweep.py
python validation_test/rf_waveguide/validation_circular_waveguide_mode_sweep.py
python validation_test/rf_waveguide/validation_tem_line_geometry_sweep.py
python validation_test/rf_waveguide/validation_waveguide_conductor_loss.py
python validation_test/rf_waveguide/validation_waveguide_te10_port_normalization.py
python validation_test/rf_waveguide/validation_radiation_pressure_poynting.py
python validation_test/rf_waveguide/validation_oblique_radiation_pressure.py
python validation_test/rf_waveguide/validation_poynting_patch_force_vector.py
python validation_test/rf_waveguide/validation_time_harmonic_maxwell_stress.py
python validation_test/rf_waveguide/validation_scattering_radiation_force.py
python validation_test/rf_waveguide/validation_two_port_scattering_momentum_force.py
python validation_test/rf_waveguide/validation_one_port_reflection_momentum_force.py
python validation_test/rf_waveguide/validation_one_port_reflection_sweep_force.py
python validation_test/rf_waveguide/validation_two_port_scattering_sweep_force.py
python validation_test/rf_waveguide/validation_two_port_sparameter_health.py
```

The validations are self-contained and use closed-form transmission-line /
waveguide formulas. Commercial RF solvers can be used internally as independent
references, but no solver-derived benchmark values are required to run or
validate these public artifacts.
