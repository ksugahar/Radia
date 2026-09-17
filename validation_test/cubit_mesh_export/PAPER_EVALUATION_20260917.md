# Cubit high-order export: comparative screening result (2026-09-17)

This is a **screening table, not yet a publication-ready superiority claim**.
It compares `cubit-mesh-export`'s Cubit-to-Netgen `.vol` route (TET and HEX)
with an independent, established Netgen/OCC native TET route. Both describe
the unit sphere, but their meshes and degrees of freedom differ. The native
route is not an alternate converter for the *same Cubit mesh*. Do not use
ratios between routes as a speed or accuracy advantage without a matched-DOF
study. Netgen/OCC has no corresponding HEX row in this comparison.

The electrostatic manufactured problem is
`-Delta(phi) = 3*pi^2*phi`, with
`phi = sin(pi*x)cos(pi*y)cos(pi*z)` prescribed on the entire surface.
The table reports the NGSolve H1 potential L2 error and the electric-field
`E=-grad(phi)` L2 error. Geometry error is integrated volume against `4*pi/3`.
The polynomial field order and Cubit curved-export order are the same. Cubit
meshing is batch/nographics; no GUI is launched.

| Route | Element | Order | Elements | DOF | Volume error | Potential L2 error | Electric-field L2 error | Labels | `check-vol` |
|---|---|---:|---:|---:|---:|---:|---:|---:|---|
| Cubit → exporter → `.vol` | TET | 1 | 231 | 75 | −10.008% | 2.539e−1 | 2.385 | 1/1 | CAD tolerance fails, as expected |
| Cubit → exporter → `.vol` | TET | 2 | 231 | 434 | −0.212% | 5.116e−2 | 8.962e−1 | 1/1 | pass |
| Cubit → exporter → `.vol` | TET | 3 | 231 | 1309 | +0.019% | 1.040e−2 | 1.753e−1 | 1/1 | pass |
| Cubit → exporter → `.vol` | HEX | 1 | 32 | 53 | −23.359% | 3.009e−1 | 2.146 | 1/1 | CAD tolerance fails, as expected |
| Cubit → exporter → `.vol` | HEX | 2 | 32 | 321 | −0.211% | 5.303e−2 | 1.011 | 1/1 | pass |
| Cubit → exporter → `.vol` | HEX | 3 | 32 | 997 | +0.131% | 2.877e−2 | 2.312e−1 | 1/1 | pass |
| Netgen/OCC native | TET | 1 | 107 | 58 | −9.376% | 6.847e−1 | 4.409 | 1/1 | not applicable¹ |
| Netgen/OCC native | TET | 2 | 107 | 278 | −0.187% | 3.427e−1 | 2.705 | 1/1 | not applicable¹ |
| Netgen/OCC native | TET | 3 | 107 | 768 | +0.016% | 2.430e−1 | 1.553 | 1/1 | not applicable¹ |

¹ `check-vol` checks an exported `.vol` and its Cubit CAD-reference sidecar;
the native OCC meshes were not exported through Cubit. All six Cubit sphere
rows had **zero invalid curved-Jacobian samples**, including the two linear
baselines whose geometry failed the 1% CAD tolerance.

The one-label sphere is too weak to assess boundary semantics. Therefore the
bundled gapped-coil/workpiece/air sample was separately exported at order 2:
**5/5 distinct boundary labels** (`source`, `sink`, `sibc`, `coil_surface`,
`outer`) had positive measured area, **3/3 materials** (`coil`, `workpiece`,
`air`) were retained, and strict-label `check-vol` passed with zero invalid
Jacobian samples. This is label preservation, not an electromagnetic solve of
the coil device.

Two independent full Cubit batch generations and two independent OCC
generations gave the same element counts, labels and field errors. The maximum
absolute repeat difference across volume, potential error and field error was
`4.44e-16`. This checks numerical repeatability on this host, not portability
across machines or Cubit releases.

## Reproduce

Run from the repository root on a licensed Cubit host with the released
`cubit-mesh-export==2.0.0`, NGSolve/Netgen `6.2.2606`, Python `3.12.10`,
Coreform Cubit `2025.12`, and the bundled exporter installed:

```powershell
python validation_test/cubit_mesh_export/paper_sphere_benchmark.py --repeats 2
```

The machine-readable [result snapshot](paper_sphere_benchmark_results.json)
is tracked beside this report. The raw per-run `.vol`, sidecar, journals,
Cubit logs and `results.json` are retained at
`C:\temp\cubit-paper-sphere-benchmark\` on LAB. Cubit returned
exit code 2 after producing each mesh; its log reported plugin/startup file
errors and `Errors found during session`. We did **not** reinterpret that exit
code as success. Independently reloaded `.vol` files, verified expected labels,
solved the FE problem, and applied `check-vol`; process-exit cleanliness remains
an operational issue to resolve before manuscript submission.

## Remaining work before a paper claim

- Compare matched-DOF or matched-element-size sequences, not just one fixed
  size per route; separate mesh-density effects from curving effects.
- Run a physically configured electromagnetic device case (coil/current and
  material interfaces), not only a manufactured electrostatic solution, and
  check field observables against an independent reference.
- Repeat on another licensed machine and retain machine/version provenance.
- Establish a clean Cubit batch exit, or document and isolate its startup
  error with a defensible artifact-acceptance rule.
