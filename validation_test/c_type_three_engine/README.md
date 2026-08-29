# C-type three-formulation validation

This suite compares the production HDiv-MMM magnetostatic route with two
independent NGSolve finite-element formulations:

1. HDiv-MMM (BDM1 or BDM2), with the Coulomb charge Gram as the exact open
   boundary operator;
2. HCurl reduced-A;
3. H1 Omega-reduced-Omega.

## Canonical geometry route

`cad/c_type_iron.jou` is the C-yoke CAD authority. It creates the ruled pole
chamfer in Cubit/ACIS and preserves the analytic Example-5 iron volume. The
validation must not reconstruct the pole with `netgen.occ`.

`build_cubit_meshes.py` writes two checked Netgen `.vol` artifacts:

- `iron.vol`: exact iron only, for HDiv-MMM. Adding finite air would replace
  the method's Coulomb open boundary with an unrelated truncation.
- `kelvin_domain.vol`: the same exact iron, a locally refined physical-air
  sphere, and a translated Kelvin sphere. Reduced-A and
  Omega-reduced-Omega share its one-to-one periodic identification. A finite
  outer air box is forbidden.

The three engines share the same solid `CoilBuilder`, B-H table, and physical
observation points. The acceptance quantity is B, not a gauge-dependent
potential. Fixed-mesh machine equality is not claimed across different FE
spaces and different open-boundary treatments; mesh and outer-domain
convergence must tighten the pairwise B discrepancy.

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

Heavy nonlinear runs belong on mdx or hibino. The mesh build remains a
foreground Cubit job, and the solver run remains a foreground Python job so a
result cannot be mistaken for a completed validation while a detached process
is still running.

## Tracked evidence

`results/lab_20260829_mesh.json` records the passing Cubit/Kelvin topology
contract. `results/lab_20260829_linear_order2.json` is the passing order-2
linear comparison. `results/lab_20260829_nonlinear_order1.json` deliberately
records a failed accuracy gate: all three engines converged, but the Omega
result remains 5.85% from HDiv in the gap core. A failed artifact is retained
so the mdx/hibino order-2 follow-up cannot silently replace or conceal it.
