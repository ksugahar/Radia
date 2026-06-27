# TEAM Problem 7: Asymmetrical Conductor with a Hole

## Problem Description

TEAM (Testing Electromagnetic Analysis Methods) Problem 7 is a benchmark problem for 3D eddy current analysis.

- **Geometry**: Aluminum plate (294 x 294 x 19 mm) with an eccentric square hole (108 x 108 mm)
- **Excitation**: Racetrack coil with 2742 A-turns
- **Frequencies**: 50 Hz and 200 Hz
- **Material**: Aluminum (σ = 3.526e7 S/m)

## References

- [COMSOL Blog: Solving TEAM Problem 7](https://www.comsol.com/blogs/solving-team-problem-7-acdc-module/)
- [FreeFEM TEAM7 Module](https://modules.freefem.org/modules/team7/)
- [Compumag TEAM Problem 7 PDF](https://www.compumag.org/jsite/images/stories/TEAM/problem7.pdf)

## Files

The historical TEAM7 Python scripts (`team7_geometry.py`, `team7_solver.py`,
`team7_coil_current.py`, `team7_A_method.py`, and
`experiment_weighted_average.py`) were promoted out of `examples` and archived
with full source and SHA-256 in
`docs/kelvin/kelvin_remaining_examples_archive_results.json`.

Future maintained TEAM7 behavior should be reintroduced either as a
`validation_test/` regression or a result-bearing notebook, not as standalone
example scripts.

## Physical Parameters

| Parameter | Value |
|-----------|-------|
| Plate dimensions | 294 x 294 x 19 mm |
| Hole dimensions | 108 x 108 mm at (18, 18) mm |
| Coil turns | 2742 |
| Coil current | 1 A/turn |
| Aluminum conductivity | 3.526e7 S/m |
| Skin depth (50 Hz) | ~12 mm |
| Skin depth (200 Hz) | ~6 mm |

## Formulation

Time-harmonic A-method (vector potential):

```
curl(1/μ · curl(A)) + jωσA = J₀
```

where:
- A: magnetic vector potential
- μ: permeability (μ₀ for all regions)
- σ: conductivity (non-zero only in aluminum plate)
- ω: angular frequency
- J₀: source current density in coil

## Measurement Points

Validation is performed along line A1-B1:
- y = 72 mm
- z = 34 mm (above the plate)
- x varies from 0 to 288 mm

## Notes

- The 50 Hz solution gives reasonable results with coarse mesh
- The 200 Hz solution requires finer mesh near the plate surface (skin depth ~6 mm)
- Direct solver (UMFPACK) is used for problems with < 150k DOFs
- GMRES with Jacobi preconditioner for larger problems
