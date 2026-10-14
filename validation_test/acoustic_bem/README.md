# Acoustic FEM/BEM validation corpus

These validation scripts keep acoustic radiation checks outside the fast unit-test
suite. They are intended as readable gates for scalar Helmholtz FEM/BEM teaching
code: closed-form first, then mesh or boundary-integral implementations can be
checked against the same numbers.

The cross-validation registry advertises
`validation_test/acoustic_bem/validation_*_summary.json` as the reusable
machine-readable artifact family. Keep runnable validation here; promote only a
polished user-facing explanation to `docs/` when there is one.

## Pulsating sphere

`validation_pulsating_sphere_radiation.py` validates the outgoing spherical-wave
solution for a uniformly pulsating sphere. It records low-frequency radiation
resistance scaling, high-frequency efficiency, exact spherical power
conservation, and near-to-far pressure decay.

## Low-frequency Helmholtz kernel

`validation_low_frequency_helmholtz_kernel.py` validates the outgoing 3D
Helmholtz Green function split into the static Laplace singular term and smooth
low-frequency corrections. It is a kernel-level gate for stable acoustic BEM
assembly before moving to full FEM/BEM coupling examples.

## Spherical DtN modes

`validation_spherical_dtn_modes.py` validates exact outgoing spherical
Helmholtz Dirichlet-to-Neumann eigenvalues and radiation impedances for
low-order modes. It is the smallest readable open-boundary FEM/BEM coupling
gate: FEM pressure trace in, exterior normal derivative out.

## Planar DtN symbol

`validation_planar_dtn_symbol.py` validates the exact outgoing half-space
Helmholtz DtN symbol for planar trace modes. It records normal/oblique
propagating modes and evanescent modes, making the FEM/BEM sign convention and
low-frequency radiation/near-field split visible without a mesh.

## Baffled piston radiation

`validation_baffled_piston_radiation.py` validates the closed-form radiation
impedance of a uniformly vibrating circular piston in an infinite baffle. It
records low-frequency resistance/reactance asymptotes, high-frequency
plane-wave limit, and radiated active power.

## Impedance to DtN bridge

`validation_acoustic_impedance_dtn_bridge.py` validates the sign convention
between acoustic specific impedance/admittance and Helmholtz DtN/Robin
coefficients. It checks planar and spherical radiation modes and a baffled
piston average impedance round trip.

## Boundary power

`validation_acoustic_boundary_power.py` validates active/reactive acoustic
boundary power from complex pressure and outward normal velocity phasors. It
checks peak/RMS conventions, a plane-wave resistive load, a purely reactive
near-field load, and a mixed impedance boundary.

## Impedance reflection and absorption

`validation_acoustic_impedance_reflection.py` validates the plane-wave
reflection coefficient, absorption coefficient, and power balance of a local
acoustic impedance boundary. It checks matched, mismatched, purely reactive,
pressure-release, and oblique-incidence limits.

## Impedance radiation pressure

`validation_acoustic_impedance_radiation_pressure.py` validates the normal
momentum pressure from an acoustic impedance reflection summary. It checks the
matched absorber, partial reflector, and lossless reactive reflector limits.

## Impedance sweep

`validation_acoustic_impedance_sweep.py` validates a frequency-indexed acoustic
impedance table. It records reflection, absorption, normal momentum pressure,
and passive-load diagnostics in the same sweep summary.

```powershell
python validation_pulsating_sphere_radiation.py
python validation_low_frequency_helmholtz_kernel.py
python validation_spherical_dtn_modes.py
python validation_planar_dtn_symbol.py
python validation_baffled_piston_radiation.py
python validation_acoustic_impedance_dtn_bridge.py
python validation_acoustic_boundary_power.py
python validation_acoustic_impedance_reflection.py
python validation_acoustic_impedance_radiation_pressure.py
python validation_acoustic_impedance_sweep.py
```
