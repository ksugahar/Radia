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
