# Surface modes under saturation (1-D, 2026-09-27)

Does deliberately adding surface modes to a bulk basis still pay off when the
material saturates?  `nonlinear_surface_modes_1d.py` answers it on the
simplest case and writes `nonlinear_surface_modes_1d.json` next to itself.

Model: `sigma dB(H)/dt = d2H/dx2` on a steel half-slab of 10 mm,
`H(0,t) = H0 sin(wt)` at 50 Hz, symmetry at the centre plane,
`B = mu0 H + Js tanh(mu0 mur H / Js)` with `mur = 1000`, `Js = 1.6 T`,
`sigma = 5 MS/m`; linear skin depth 1.0 mm.  Backward Euler and Newton on
2,000 cells, 800 steps per period, four periods.

Metric: the relative space-time L2 error of projecting the last period of
`H(x,t)` onto each basis.  This is the best any Galerkin model on that basis
could reach, not the error of a solved reduced model (a nonlinear ROM also
needs hyper-reduction).

| H0 [A/m] | surface B | bulk 7 | bulk 5 + linear surface 2 | bulk 5 + surface POD 4 | bulk 5 + linear 2 + POD 2 |
|---|---|---|---|---|---|
| 10 | 0.01 T | 6.1e-2 | 2.1e-4 | 1.7e-3 | 8.3e-5 |
| 1,000 | 1.05 T | 5.9e-2 | 3.9e-3 | 2.0e-3 | 1.9e-3 |
| 3,000 (train) | 1.58 T | 5.5e-2 | 2.1e-2 | 8.6e-3 | 8.9e-3 |
| 10,000 | 1.61 T | 5.8e-2 | 4.2e-2 | 1.5e-2 | 1.6e-2 |
| 100,000 | 1.73 T | 2.2e-2 | 2.3e-2 | 1.5e-2 | 1.5e-2 |

Bulk: Legendre polynomials.  Linear surface modes: real and imaginary parts
of `exp(-(1+i) x/delta)` with the linear depth.  Surface POD: SVD of training
snapshots at 30, 300, 3,000 and 30,000 A/m, each normalized to unit norm,
after removing their bulk (and, for the last column, linear-surface) content.

Findings:

- The linear surface modes beat the same number of bulk modes by 300x in the
  linear range, 2.6x at 1.58 T and not at all at 1.73 T: under saturation the
  field penetrates as a front whose depth and shape the linear exponential no
  longer match.
- At equal dimension nine, bulk 5 + POD 4 improves on bulk 9 by about
  2.3x at 3,000 A/m and 2.1x at 10,000 A/m. At 100,000 A/m, bulk 9
  is slightly better (0.01447 versus 0.01478); there is no uniform gain.
- Linear surface 2 + POD 2 is best in the low-amplitude rows. POD 4 is
  better from 3,000 to 30,000 A/m, and bulk 9 is best at 100,000 A/m.
- Unnormalized POD training lets the large amplitudes swamp the linear shapes
  (8.2e-2 at 10 A/m in an earlier run); normalize per amplitude.

The four-period reference has no recorded cycle-to-cycle or mesh/time-step
convergence study; “periodic-steady” in the original script is an assumption.
1-D only; corners and edges of 3-D conductors are not represented.
