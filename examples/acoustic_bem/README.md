# Acoustic FEM/BEM validation examples

These examples keep acoustic radiation checks outside the fast unit-test suite.
They are intended as readable gates for scalar Helmholtz FEM/BEM teaching code:
closed-form first, then mesh or boundary-integral implementations can be checked
against the same numbers.

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

```powershell
python validation_pulsating_sphere_radiation.py
python validation_low_frequency_helmholtz_kernel.py
python validation_spherical_dtn_modes.py
```
