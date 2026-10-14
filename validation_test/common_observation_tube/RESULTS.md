# Retained result — 2026-08-21

The retained four-formulation comparison passed all seven gates in
`results/comparison_report.json`.

## Discretization

- common tube: 41 longitudinal stations from -55 mm to +55 mm;
- transverse sampling: 24 points on a 6 mm-radius circle;
- legacy Radia reference: rectangular segmentation 5 and 7, with level 7
  selected (12,348 interaction degrees of freedom);
- HDiv-MMM: 1,954 tetrahedra, 13,332 BDM1 degrees of freedom, 31 iterations;
- reduced-A: order 2, 13,800 tetrahedra, 73,564 HCurl degrees of freedom;
- Omega-reduced-Omega: order 2, 17,758 tetrahedra, 26,437 H1 degrees of
  freedom, periodic Kelvin exterior.

The legacy-Radia field changed by 2.617% RMS from segmentation 5 to 7, below
the 3% reference-convergence gate.

## Error relative to the legacy Radia reference

| Formulation | pointwise B RMS | integrated-field error | maximum vector error |
|---|---:|---:|---:|
| HDiv-MMM | 2.99% | 4.04% | 2.21 mT |
| reduced-A | 4.94% | 6.61% | 9.81 mT |
| Omega-reduced-Omega | 4.00% | 4.42% | 7.63 mT |

The integrated main-field coefficients were 0.00343309 T m (Radia),
0.00330948 T m (HDiv-MMM), 0.00322988 T m (reduced-A), and 0.00333116 T m
(Omega-reduced-Omega).

## Integrated multipoles at the 6 mm reference radius

The square yoke and vertical excitation permit odd normal components.  The
Radia reference gave normal sextupole and decapole components of -211.33 and
-35.81 units.  The corresponding values were:

| Formulation | normal sextupole | normal decapole |
|---|---:|---:|
| HDiv-MMM | -157.32 units | -24.25 units |
| reduced-A | -182.88 units | -23.68 units |
| Omega-reduced-Omega | -122.22 units | -55.18 units |

All complex coefficients through order six and their differences from Radia
are retained in the JSON report. `results/integrated_multipoles.csv` is
regenerated as a plotting convenience and is not separate validation truth.
Higher multipoles are more sensitive than the integrated dipole to mesh and
open-boundary discretization, so these values are comparison baselines rather
than final high-accuracy multipole certificates.

## Scope

This is the common observation-tube infrastructure proof on a finite linear
square yoke with a beam aperture.  The next retained magnet models are the
C-type dipole and quadrupoles.  The hybrid undulator is explicitly outside the
requested validation scope.
