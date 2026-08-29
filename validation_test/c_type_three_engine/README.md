# C-type formulation validation

The primary comparison is the production HDiv-MMM magnetostatic route against
NGSolve Omega-reduced-Omega. HCurl reduced-A is retained as an independent
third-formulation audit:

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

The engines share the same solid `CoilBuilder`, B-H table, and physical
observation points. The acceptance quantity is B, not a gauge-dependent
potential. Fixed-mesh machine equality is not claimed across different FE
spaces and different open-boundary treatments; mesh and outer-domain
convergence must tighten the pairwise B discrepancy.

For nonlinear comparisons, HDiv-MMM and Omega-reduced-Omega use the same
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

Use `--primary-only` for the faster direct HDiv-MMM versus
Omega-reduced-Omega production comparison. Omitting it also runs reduced-A as
the independent third route. The pass/fail accuracy metric is always the
primary pair; every selected nonlinear engine must also converge.

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
result remains 5.85% from HDiv in the gap core.
`results/lab_20260829_nonlinear_order2_primary.json` records the passing direct
comparison after unifying the B-H interpolation: HDiv-MMM and
Omega-reduced-Omega differ by 0.18032%. The failed order-1 artifact is retained
so the order-convergence result cannot silently replace or conceal it. LAB
timings are correctness evidence only; publish performance after repeating the
same release on idle mdx and hibino.
