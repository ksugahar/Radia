# Cubit mesh export: final comparison table (2026-09-18)

This is the final table for the declared two-case protocol, not a universal
claim that one mesher is more accurate. Cubit-to-Netgen .vol is compared
with independently generated Netgen/OCC meshes. The paths start from the same
analytic CAD shapes but not the same mesh connectivity. TET degrees of
freedom were brought into the same range (sphere: 75/434/1309 versus
81/395/1099; capacitor: 977/7007/22783 versus 1030/6823/21536).
Each route keeps its own topology fixed across polynomial orders.
Netgen/OCC has no corresponding curved-HEX route in this protocol.

Case S is the unit sphere. Geometry error is integrated volume against 4π/3.
The electrostatic manufactured problem is -Δφ=3π²φ, with
φ=sin(πx)cos(πy)cos(πz) as exact full-boundary Dirichlet data. Potential
and electric-field E=-grad(φ) errors are L2 norms on the represented mesh.

Case C is a spherical capacitor: inner/outer radii a=0.3, b=1, unit
permittivity, electrode voltages 1/0. Exact capacitance is
C=4πab/(b-a)=5.385587406 in normalized units. Numerical capacitance
is field energy for the 1-V difference. This is a physical Laplace boundary
value problem with two distinct electrode labels, not a manufactured source.

| Case | Route | Element | p | Elements | DOF | Geometry / capacitance error¹ | Potential L2 error | E-field L2 error | Boundary labels | check-vol² |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---|
| S | Cubit → exporter | TET | 1 | 231 | 75 | −10.008% | 2.539e−1 | 2.385 | 1/1 | CAD tolerance fail |
| S | Cubit → exporter | TET | 2 | 231 | 434 | −0.212% | 5.116e−2 | 8.962e−1 | 1/1 | pass |
| S | Cubit → exporter | TET | 3 | 231 | 1309 | +0.019% | 1.040e−2 | 1.753e−1 | 1/1 | pass |
| S | Cubit → exporter | HEX | 1 | 32 | 53 | −23.359% | 3.009e−1 | 2.146 | 1/1 | CAD tolerance fail |
| S | Cubit → exporter | HEX | 2 | 32 | 321 | −0.211% | 5.303e−2 | 1.011 | 1/1 | pass |
| S | Cubit → exporter | HEX | 3 | 32 | 997 | +0.131% | 2.877e−2 | 2.312e−1 | 1/1 | pass |
| S | Netgen/OCC native | TET | 1 | 156 | 81 | −6.876% | 5.412e−1 | 3.275 | 1/1 | n/a |
| S | Netgen/OCC native | TET | 2 | 156 | 395 | −0.099% | 1.629e−1 | 1.671 | 1/1 | n/a |
| S | Netgen/OCC native | TET | 3 | 156 | 1099 | +0.009% | 3.605e−2 | 3.859e−1 | 1/1 | n/a |
| C | Cubit → exporter | TET | 1 | 4691 | 977 | +2.510% | 1.782e−2 | 6.308e−1 | 2/2 | CAD tolerance fail |
| C | Cubit → exporter | TET | 2 | 4691 | 7007 | +0.018% | 1.461e−3 | 7.973e−2 | 2/2 | pass |
| C | Cubit → exporter | TET | 3 | 4691 | 22783 | +0.011% | 1.456e−4 | 1.027e−2 | 2/2 | pass |
| C | Netgen/OCC native | TET | 1 | 4155 | 1030 | +4.973% | 3.327e−2 | 7.198e−1 | 2/2 | n/a |
| C | Netgen/OCC native | TET | 2 | 4155 | 6823 | +0.135% | 2.590e−3 | 1.094e−1 | 2/2 | n/a |
| C | Netgen/OCC native | TET | 3 | 4155 | 21536 | +0.014% | 2.906e−4 | 1.695e−2 | 2/2 | n/a |

¹ Case S reports signed sphere-volume error; Case C reports signed
capacitance error. These quantities must not be compared across cases.
² The exporter gate consumes the Cubit CAD-reference sidecar and applies its
1% geometry tolerance. Linear p=1 meshes fail that geometry threshold as
expected; all nine Cubit rows have zero invalid Jacobian samples. The OCC route
does not produce a Cubit sidecar, so this gate does not apply.

A separate gapped-coil/workpiece/air export tests realistic naming: all 5/5
boundary labels (source, sink, sibc, coil_surface, outer) have positive
integrated area, all 3/3 materials (coil, workpiece, air) are retained, and
strict-label check-vol passes with zero invalid Jacobian samples. This label
case is not itself a device-field solve.

## Reproducibility and provenance

- LAB, Windows Server 2022; Cubit 2025.12 student license; Python 3.12.10;
  cubit-mesh-export 2.0.0; Netgen/NGSolve 6.2.2606. All Cubit calls used
  -batch -nographics -nojournal, with its command-plugin directory explicitly
  supplied. No Cubit GUI was launched.
- Two independent full Cubit generations and two independent OCC generations
  for each case produced identical topology, labels, and reported values to
  roundoff. The maximum absolute repeat difference across monitored volume,
  capacitance, potential and field errors was 2.66e−15.
- Cubit ended with code 0 on every benchmark generation. A missing
  -commandplugindir caused the previous code-2 startup diagnostic;
  explicit configuration eliminated it. The packaged cubit-smoke-test
  command now uses the same argument and passes its full solver-ready gate.
- Raw journals, logs, sidecars, and meshes are retained on LAB under
  C:\temp\cubit-paper-sphere-benchmark\ and
  C:\temp\cubit-paper-capacitor-benchmark\. The committed numerical
  snapshots are [sphere results](paper_sphere_benchmark_results.json) and
  [capacitor results](paper_capacitor_benchmark_results.json).

Run from the repository root on a licensed Cubit host:

    python validation_test/cubit_mesh_export/paper_sphere_benchmark.py --repeats 2
    python validation_test/cubit_mesh_export/paper_capacitor_benchmark.py --repeats 2

This table establishes high-order geometry transfer, label preservation,
solver convergence, and local rerun reproducibility for these cases. It does
not establish statistical performance superiority, cross-machine agreement,
or applicability to all element families and industrial geometries.
