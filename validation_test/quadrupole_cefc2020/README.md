# CEFC 2020 quadrupole ("Q-mag"): HDiv-MMM on an iron-dominated quadrupole

The laboratory's CEFC 2020 test quadrupole is the second iron-dominated quadrupole
validation of HDiv-MMM after ESRF example 6.  It exercises the same conforming
all-HEX charge-Gram route on a compact yoke with two racetrack coils per pole, over
a linear permeability sweep and a nonlinear current sweep.

## Geometry and sources

* Iron: `qmag_iron.step`, the laboratory's ACIS CAD of the yoke exported by Coreform
  Cubit 2025.12 after scaling to metres (24 solids, all extrusions along the beam
  axis z, length 60 mm, yoke 176 mm square, pole face at 20 mm on the +-x / +-y
  axes, total volume 1.0247373e-3 m^3; the STEP re-imports with the same volume).
* Coils: two Biot-Savart racetracks per pole (`qmag_case.build_qmag_coils`), the
  `radObjRaceTrk` definition of the original design notebook -- inner radii 2..15 mm,
  straights 60 x 26 mm, height 40 mm on the pole coordinate 30..70 mm; outer radii
  15..31 mm, height 20 mm on 50..70 mm.  The +-x coils carry the opposite azimuthal
  sign of the +-y coils, giving `B_x = G x`, `B_y = -G y` with `G > 0`
  (`tests/test_qmag_cefc2020_case.py` locks the pattern).  Design current density
  3 A/mm^2; coil-only gradient 1.03 T/m.
* Iron law: `iron_bh_table.json` (49 monotone `[H, B]` samples) for the nonlinear
  cases; linear `mu_r` (100 / 1000 / 10000) for the linear cases.

## Observable

`B_perp = -(B_x - B_y)/sqrt(2)` on the diagonal `x = y`, `z = 0`, 31 points at 1 mm
spacing (`r = -15..15 mm`); the summary value is `B_perp(15 mm)` and the gradient
`B_perp / r`.  Internal checks: finite field, odd symmetry of `B_perp` along the
diagonal (quadrupole), converged nonlinear loop; mesh convergence is judged across
the `qmag_h*.vol` series (15 / 10 / 6 / 4 mm: 980 / 2352 / 6040 / 16408 HEX).

## Acceptance: agreement with the mixed total/reduced Omega FEM

HDiv-MMM is accepted on this magnet when it agrees with the repository's own
mixed total/reduced Omega formulation (and optionally HCurl reduced-A) on the
same coils and the same iron law -- the same three-engine contract as the
C-type and ESRF lanes.  `build_qmag_fem_mesh.py` builds the full-domain FEM
mesh (iron + physical air sphere of 0.16 m + periodic Kelvin exterior, tets,
order 2, bore refinement; ~183k elements) with the ESRF coil-yoke builder, and
`run_qmag_three_engine.py` runs both formulations on the 31 diagonal points
(the FEM field cube-averaged around each point because the points lie on the
z = 0 mesh seam) and reports the relative RMS of the vector field difference
and the `B_perp(15 mm)` gap against a 3 % band.

## Running

```powershell
python build_qmag_cubit_mesh.py --output-dir C:\temp\radia-qmag\meshes --sizes 0.010 0.006 0.004
python run_qmag_hdiv.py --mesh C:\temp\radia-qmag\meshes\qmag_h10.vol --case mu1000 --output results\qmag_h10_mu1000_hibino.json
python run_qmag_hdiv.py --mesh C:\temp\radia-qmag\meshes\qmag_h10.vol --case J3.0 --output results\qmag_h10_J3.0_hibino.json
python build_qmag_fem_mesh.py --output-dir C:\temp\radia-qmag\fem
python run_qmag_three_engine.py --fem-mesh C:\temp\radia-qmag\fem\coil_yoke_kelvin.vol --fem-mesh-report C:\temp\radia-qmag\fem\qmag_fem_kelvin.mesh.json --case mu1000 --hdiv-result results\qmag_h10_mu1000_hibino.json --output results\three_engine_mu1000_hibino.json
```

The HEX mesh builder writes conforming z-swept HEX (imprint/merge, per-volume
explicit sweep, `check-vol` and `radia.vim.mesh_conformity_report` gates).  The
solves are heavy (one charge Gram per HDiv case, a 183k-element Kelvin FEM per
formulation) and run on hibino, one job at a time; results are committed under
`results/` with the host name in the file name.

## Results (hibino, 2026-09-07)

Linear cases, HDiv-MMM on the `h = 10 mm` conforming HEX mesh (2352 BDM1
elements, 60,816 unknowns, mass-Riesz CG to `1e-8`) against the mixed
total/reduced Omega FEM on the 183k-element order-2 Kelvin mesh.  The RMS is
the relative RMS of the vector field difference over the 31 diagonal points;
the last column is the relative gap of the diagonal field at 15 mm.  All
three cases pass the 3 % band by a wide margin; the JSONs under `results/`
carry the full point-wise fields and the mesh hashes.

| case | `mu_r` | CG iterations | field RMS | `B_perp(15 mm)` HDiv-MMM [T] | mixed Omega FEM [T] | gap |
|---|---|---|---|---|---|---|
| `mu100` | 100 | 30 | 0.08 % | -0.17875 | -0.17858 | 0.09 % |
| `mu1000` | 1000 | 33 | 0.06 % | -0.22722 | -0.22725 | 0.01 % |
| `mu10000` | 10000 | 36 | 0.07 % | -0.23394 | -0.23401 | 0.03 % |
| `J3.0` (B(H) table, 3 A/mm^2) | nonlinear | 2 Newton / 866 CG | 0.05 % | -0.23401 | -0.23393 | 0.03 % |

The nonlinear row is the tabulated B(H) law at the design current: HDiv-MMM
converged its energy Newton loop (mass-Riesz CG inner solves); the mixed Omega
Picard loop needed the constrained Anderson mixing (`--mixed-anderson-depth 2`,
52 iterations, 126 min on hibino) after a first attempt with
plain damped Picard stalled at the 80-iteration cap with a relative B change of
7.8e-4 -- the runner now saves that partial state and resumes from it.

## Multipole convergence: curved against straight pole faces, BDM1 against BDM2

`B_perp(15 mm)` is dominated by the gradient and cannot tell a curved pole face
from a faceted one (the curved and straight TET BDM2 solves agree to five
digits on it).  The observable that can is the allowed-harmonic content of the
field: `run_qmag_multipoles.py` expands `B_y + i B_x` on a circle of radius
15 mm (75 % of the pole radius) by FFT and reports `b_6, b_10, b_14` in units
of `1e-4 B_2`, plus the quadrupole-forbidden `b_3, b_4, b_5` as a mesh-symmetry
check.  Linear `mu_r = 1000`, design current; every mesh from the same Cubit
import (`build_qmag_cubit_mesh.py`, HEX swept or `--scheme tet`, exported at
curve order 2 or 1).  Correctness study on LAB, no timings:

| route | h [mm] | unknowns | `b_6` | `b_10` | `b_14` | forbidden `b_3` |
|---|---|---|---|---|---|---|
| HEX curved Q2, BDM1 | 15 | 25,984 | +21.32 | -4.01 | +0.37 | 0.000 |
| HEX curved Q2, BDM1 | 10 | 60,816 | +26.04 | -2.08 | +0.76 | 0.000 |
| HEX curved Q2, BDM1 | 6 | 153,296 | +26.17 | -1.97 | +0.80 | 0.000 |
| HEX curved Q2, BDM1 | 4 | 410,128 | +26.14 | -1.97 | +0.80 | 0.000 |
| TET straight, BDM1 | 10 | 80,202 | +26.20 | -1.84 | +0.85 | 0.024 |
| TET straight, BDM1 | 7 | 172,644 | +26.25 | -1.82 | +0.89 | 0.004 |
| TET straight, BDM1 | 5 | 397,047 | +26.10 | -1.85 | +0.83 | 0.003 |
| TET straight, BDM2 | 10 | 233,892 | +25.69 | -1.86 | +0.80 | 0.003 |
| TET curved Q2, BDM2 | 10 | 233,892 | +24.88 | -2.31 | +0.61 | 0.015 |
| TET curved Q2, BDM2 | 7 | 507,210 | +26.05 | -1.73 | +0.82 | 0.022 |

Both the curved HEX route and the straight TET BDM1 route converge to
`b_6 = 26.1 +- 0.1`; the faceted pole face at 10 mm and below costs less than
0.1 unit in `b_6` and about 0.12 unit in `b_10` (HEX `-1.97` against straight
TET `-1.85`, the curvature showing first in the higher order).  The 15 mm HEX
mesh is off by 5 units through field resolution, not geometry: at 10 mm the
Q2 pole face is already exact to the pole sag of order `h^3 / R^2`, so a
higher curve order buys nothing here, and the curved touching-block
quadrature would only get dearer.  The swept HEX meshes keep the quadrupole
symmetry exactly (forbidden harmonics at round-off); the TET meshes leak up
to 0.04 unit.  The curved TET BDM2 solve at 10 mm sits 1.2 units below the
converged `b_6` -- farther than the straight BDM2 on the same tets -- and
reaches 26.05 at 7 mm while its `b_10` overshoots to `-1.73`: the curved TET
route converges to the same limit but from further away, so at a given mesh
its curved touching-block rule, not the geometry, sets the error.  That rule
is also 87 % of the curved TET build time (section 8.14 of the review): the
next TET-side lever is a cheaper and more accurate curved touching family.
