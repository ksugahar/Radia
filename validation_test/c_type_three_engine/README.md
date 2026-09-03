# C-type formulation validation

The primary comparison is the production HDiv-MMM magnetostatic route against
the NGSolve TOSCA-style mixed total/reduced Omega route. HCurl reduced-A is
retained as an independent third-formulation audit:

1. HDiv-MMM (BDM1 or BDM2), with the Coulomb charge Gram as the exact open
   boundary operator;
2. HCurl reduced-A;
3. H1 TOSCA mixed total/reduced Omega.

## Canonical geometry route

`cad/c_type_iron.jou` is the C-yoke CAD authority. It creates the ruled pole
chamfer in Cubit/ACIS and preserves the analytic Example-5 iron volume. The
validation must not reconstruct the pole with `netgen.occ`.

`build_cubit_meshes.py` writes two checked Netgen `.vol` artifacts:

- `iron.vol`: exact iron only, for HDiv-MMM. Adding finite air would replace
  the method's Coulomb open boundary with an unrelated truncation.
- `kelvin_domain.vol`: the same exact iron, a locally refined physical-air
  sphere, and a translated Kelvin sphere. Reduced-A and
  mixed total/reduced Omega share its one-to-one periodic identification. A finite
  outer air box is forbidden.

The engines share the same solid `CoilBuilder`, B-H table, and physical
observation points. The acceptance quantity is B, not a gauge-dependent
potential. Fixed-mesh machine equality is not claimed across different FE
spaces and different open-boundary treatments; mesh and outer-domain
convergence must tighten the pairwise B discrepancy.

For nonlinear comparisons, HDiv-MMM and mixed total/reduced Omega use the same
monotone PCHIP B(H) interpolation and continue beyond the table with vacuum
slope. Sharing only the table samples is not considered a shared material law.

The acceptance gate uses the median-plane-projected B field in the useful gap
core (`|x| <= 10 mm` by default). The artifact also stores the raw full-tube
comparison, projected full-fringe comparison, and each engine's symmetry
defect. A passing core result therefore cannot hide an unconverged fringe or
an asymmetric mesh; those remain explicit convergence diagnostics.

## Run

Build meshes on a Cubit 2025.12 host:

```powershell
python validation_test/c_type_three_engine/build_cubit_meshes.py `
  --output-dir C:/temp/radia_ctype_three_engine/meshes
```

Run a fast linear preflight, then the nonlinear production comparison:

```powershell
python validation_test/c_type_three_engine/run_three_engine.py `
  --mesh-dir C:/temp/radia_ctype_three_engine/meshes `
  --mode linear `
  --output C:/temp/radia_ctype_three_engine/linear.json

python validation_test/c_type_three_engine/run_three_engine.py `
  --mesh-dir C:/temp/radia_ctype_three_engine/meshes `
  --mode nonlinear `
  --reduced-a-solver direct `
  --output C:/temp/radia_ctype_three_engine/nonlinear.json
```

Use `--primary-only` for the faster direct HDiv-MMM versus mixed-Omega
development comparison. Omitting it also runs reduced-A as the independent
third route. A full three-engine run is accepted only if every selected pair
passes the requested gap-core `B` tolerance and, for nonlinear runs, every
engine converges. The default `--source-trace-tolerance 0.05` is a separate
gate for both the physical iron/air and physical-air/Kelvin source-potential
traces. A failed trace projection requires an explicit cut/cohomology
representation rather than a relaxed numerical tolerance.

Add `--resume` for remote production runs. The runner writes a hash-checked
checkpoint after each of the HDiv, reduced-A, and Omega engines and emits one
JSON progress event at each engine boundary. A checkpoint with different
mesh, material, order, mode, ChargeGram tolerance, or observation points is
rejected rather than silently reused. The comparison defaults to
`--hdiv-gram-eps 1e-14`; this keeps H-matrix truncation below the symmetry and
cross-formulation accuracy being measured.

Nonlinear runs also record convergence, iteration count, final relative
change, tolerance, and iteration limit for every engine. Pairwise agreement is
not a pass when any engine is unconverged. Use `--nonlinear-verbose` for remote
progress logs.

Heavy nonlinear runs belong on hibino first, with mdx allowed only when hibino
is unavailable and its CI queue is idle. The mesh build remains a
foreground Cubit job, and the solver run remains a foreground Python job so a
result cannot be mistaken for a completed validation while a detached process
is still running.

For an accuracy certificate, build the default four-level geometric family and
run the convergence driver. The levels use scales `1.25, 1.0, 0.8, 0.64`; the
last three define the observed order and Richardson estimate. The fourth level
is intentional: one unusually small increment between two independent
unstructured Cubit meshes must not be mistaken for asymptotic convergence.

```powershell
python validation_test/c_type_three_engine/build_mesh_family.py `
  --output-dir C:/temp/radia_ctype_accuracy/meshes

python validation_test/c_type_three_engine/run_mesh_convergence.py `
  --mesh-family C:/temp/radia_ctype_accuracy/meshes/mesh_family.json `
  --output C:/temp/radia_ctype_accuracy/mdx_certificate.json `
  --replicate-final-result C:/temp/radia_ctype_accuracy/hibino_finer.json `
  --resume
```

The independent-host result is mandatory. Its mesh, software version,
implementation hashes, comparison contract, and observation points must match
the mdx finest result. The certificate bounds discretization and
cross-formulation spread; it does not rename their agreement as an unavailable
analytic solution.

## Tracked evidence

`results/lab_20260829_mesh.json` records the passing Cubit/Kelvin topology
contract. The `20260829` and `20260830` field artifacts remain tracked as
historical evidence for the former global reduced-Omega route only. They do
not certify the current TOSCA mixed formulation, because they omit its required
physical-air/Kelvin source-potential jump.

The current pre-release evidence is explicit about native binary provenance:
the current Python mixed-formulation source was overlaid on the installed
`radia 4.95.77` wheel on Hibino; no `.pyd` was copied. The native HDiv kernel,
CoilBuilder, Radia source evaluation, and NGSolve assembly were therefore the
installed wheel's components.

- `results/hibino_20260902_linear_order3_tosca_mixed_v3.json` is the linear
  order-3 full three-engine run. At a 1% all-pair gate its HDiv/mixed,
  HDiv/reduced-A, and mixed/reduced-A gap-core RMS differences are 0.41955%,
  0.38661%, and 0.45977%. The physical source-trace residuals are 1.26933%
  on iron/air and 1.75381% on `kelvin_int`, both below the 5% cut gate.
- `results/hibino_20260902_nonlinear_order2_tosca_mixed_v3.json` is the
  nonlinear full three-engine run with the shared monotone PCHIP B(H) table.
  All engines converge; the three respective gap-core RMS differences are
  0.12324%, 0.10674%, and 0.16023%. Its trace residuals are 2.98511% and
  4.53077%, both below the same 5% gate. The HDiv, reduced-A, and mixed-Omega
  runtimes are 11.15 s, 221.56 s, and 135.76 s.

These two artifacts complete the fixed-mesh three-formulation acceptance for
the current TOSCA split. A new multi-level, independent-host convergence
certificate must be generated from the released binary before stating an
absolute numerical envelope for this formulation; the historic global-Omega
certificate must not be relabelled as that result.
