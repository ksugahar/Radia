# Four-formulation common-observation-tube comparison

This validation compares the same physical magnetic flux density on one beam
observation tube.  It deliberately keeps the four field engines independent:

1. the unmodified non-MPI ESRF Radia 2023 Python 3.8 BEM extension, with
   two-level reference-discretization convergence;
2. BDM1 HDiv-MMM on an iron-only mesh;
3. reduced-A on an HCurl finite-element space;
4. Omega-reduced-Omega on an H1 space with a periodic Kelvin exterior.

The model is a finite 40 mm square, 60 mm long, linear iron yoke with a 20 mm
square beam aperture.  The observation tube follows the aperture axis from
-55 mm to +55 mm and samples a 6 mm-radius circle.  The source is a uniform
vertical 0.1 T field and the yoke has relative permeability 100.

Run a quick local lane:

```powershell
python validation_test/common_observation_tube/compare_four_formulations.py `
  --output-dir C:/temp/radia_common_tube_smoke `
  --stations 17 --circle-points 16 --radia-segmentation 2 --fe-order 1 `
  --gate-profile smoke
```

The smoke profile checks execution and finite fields only. It still records
the retained accuracy diagnostics, but coarse meshes do not fail the run for
missing the publication comparison bands.

Run the retained comparison:

```powershell
python validation_test/common_observation_tube/compare_four_formulations.py `
  --output-dir validation_test/common_observation_tube/results `
  --stations 41 --circle-points 24 --radia-segmentation 7 --fe-order 2
```

The legacy oracle needs Python 3.8 and the unmodified ESRF Radia extension.
The runner finds Python through `uv python find 3.8`; set
`RADIA_LEGACY_PYTHON38` when that interpreter is elsewhere. Set
`RADIA_LEGACY_EXTENSION` when the extension is not at the LAB default path.

`comparison_report.json` is the machine-readable verdict.
It records the Radia, NGSolve, and Python versions plus the SHA-256 identity of
the independent legacy-Radia binary; workstation-specific binary paths are not
part of the retained artifact.
`longitudinal_main_field.csv` stores the complex dipole coefficient along the
tube.  `integrated_multipoles.csv` stores the integrated normal and skew
multipoles through order six, including accelerator units at the 6 mm
reference radius. The CSV files are derived plotting outputs; the JSON files
are the retained validation evidence.
