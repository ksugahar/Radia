# C-type formulation validation

This comparison lane is three-way: all three engines are mandatory to claim
that this lane passed. It is not an unconditional HDiv-MMM release prerequisite;
see [production acceptance](../../docs/hdiv_vim/PRODUCTION_ACCEPTANCE.md) for
scope-specific acceptance when a comparator is unvalidated or resource-limited.
Such runs retain their failed or incomplete comparison status.

The physical C-type model
is evaluated by HDiv-MMM, HCurl reduced-A, and the NGSolve H1 mixed
total/reduced Omega route:

1. HDiv-MMM (BDM1 or BDM2), with the Coulomb charge Gram as the exact open
   boundary operator;
2. HCurl reduced-A;
3. H1 mixed total/reduced Omega.

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

Every run executes all three formulations. It is accepted only if every
gap-core `B` pair passes the requested tolerance and, for nonlinear runs,
every engine converges. The default `--source-trace-tolerance 0.05` is a separate
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
not certify the current mixed total/reduced Omega formulation, because they omit its required
physical-air/Kelvin source-potential jump.

The current pre-release evidence is explicit about native binary provenance:
the current Python mixed-formulation source was overlaid on the installed
`radia 4.95.77` wheel on Hibino; no `.pyd` was copied. The native HDiv kernel,
CoilBuilder, Radia source evaluation, and NGSolve assembly were therefore the
installed wheel's components.

- `results/hibino_20260903_linear_order3_mixed_omega_v4.json` is the linear
  order-3 full three-engine run. At a 1% all-pair gate its HDiv/mixed,
  HDiv/reduced-A, and mixed/reduced-A gap-core RMS differences are 0.41955%,
  0.38661%, and 0.45977%. The physical source-trace residuals are 1.26933%
  on iron/air and 1.75381% on `kelvin_int`, both below the 5% cut gate.
- `results/hibino_20260903_nonlinear_order2_mixed_omega_v4.json` is the
  nonlinear full three-engine run with the shared monotone PCHIP B(H) table.
  All engines converge; the three respective gap-core RMS differences are
  0.12324%, 0.10674%, and 0.16023%. Its trace residuals are 2.98511% and
  4.53077%, both below the same 5% gate. The HDiv, reduced-A, and H1
  mixed total/reduced Omega
  runtimes are 11.27 s, 216.93 s, and 134.71 s.

These two artifacts complete the fixed-mesh three-formulation acceptance for
the current mixed total/reduced Omega split. The released-binary four-level campaign is recorded
in `results/c_type_20260903_nonlinear_bdm2_mesh_convergence_certificate.json`
with its portable level artifacts:

- finest all-pair gap-core RMS: `0.27714%`;
- maximum discretisation uncertainty: `0.17601%`;
- conservative combined numerical envelope: `0.35399%`;
- mdx/Hibino independent-host replay: `5.25e-14` relative RMS.

All certificate checks pass for `radia 4.95.77`, including nonlinear
convergence, mesh contraction, cross-formulation agreement, and independent
host reproducibility. It certifies numerical agreement, not analytic absolute
truth. A later implementation hash must rerun the campaign; the historical
global-Omega certificate must not be relabelled as this result.

## Kelvin exterior of the mixed total/reduced Omega engine (2026-09-08)

Until `ad5d9d553` the mixed total/reduced Omega engine's Kelvin ball carried no
unknown at all: its total-potential space was restricted to the total
materials, but the material touching the identified physical sphere is the
reduced source enclosure, so that space owned no degree of freedom there and
the periodic identification paired nothing. The recorded open boundary was in
fact a Dirichlet truncation at the sphere. The four `results/
c_type_20260908_linear_kelvin_*_mdx1.json` artifacts are the same linear
order-2 run before and after the fix on the same two meshes, measured against
HDiv-MMM whose open boundary is exact.

| gap-core relative RMS | coarse before | coarse after | medium before | medium after |
|---|---|---|---|---|
| HDiv-MMM vs mixed Omega | 0.53990% | 0.15236% | 0.45521% | 0.11700% |
| HDiv-MMM vs reduced-A | 0.09819% | 0.09819% | 0.10879% | 0.10879% |
| reduced-A vs mixed Omega | 0.57760% | 0.16134% | 0.50904% | 0.11569% |

The HDiv/reduced-A pair is the control and is bit-identical across the pair,
as it must be. Before the fix the Omega error barely moved under refinement,
which is the signature of a truncation floor; after it converges like the other
two and reaches the reduced-A level. The `after` artifacts carry a different
`radia.kelvin_solver` implementation hash, which is how they are told apart.

## Gap-thickness mesh family (2026-09-11)

`--gap-size` is a pole-face SURFACE size: it refines the gap air in plane
while a tetrahedron stays free to span the whole half-gap in z.  On the scale
family above, lines through the gap crossed 6 to 18 elements but the longest
stay inside one element was 2.7 to 2.9 mm of the 10 mm gap at every level, so
the family did not control the field's resolution through the gap thickness.

`build_cubit_meshes.py` now measures what a line meets on its way across the
gap (`gap_inventory.line_profile`): the exact straight-tetrahedron clipping
on five representative lines, the same measurement on the curved cells by
bisection, two sampled counts, and coverage.  The acceptance requires a
minimum crossing count `--gap-elements-across N` on both geometries and a
longest stay of at most `--gap-segment-factor` x gap_height/N.  Three probe
traps are handled rather than trusted: NGSolve's point locator accepts a point
up to about 1e-4 reference barycentric OUTSIDE the returned element (a
tolerance-band answer is never read as a membership); pieces shorter than
1 nm are grazing contacts, not crossings; and both the curved walk and the
sampled counts are seeded with the straight-segment midpoints so thin elements
are not stepped over.  `results/lab_20260911_gap_probe_locator_regression.json`
records the probe on the four scale-family meshes: all twenty lines agree
between the four counts and the coverage deviation is zero.

`--gap-layers N` (even) cuts the positive-z gap air under a 60 x 50 mm column
into N/2 parallel slabs, meshes them at gap_height/N, and reflects the meshed
slabs with the rest of the half.  Two Cubit facts shape the construction:
the default merge tolerance is 5e-4 MODEL units (0.5 mm in this metre model),
and the imprint has an absolute tolerance of about the same size that no
setting changes, so slabs thinner than 0.5 mm are collapsed.  A slabbed half is
therefore built, imprinted, merged, meshed and reflected in a working frame
scaled by 1000 and scaled back to metres; the unlayered path keeps the
historical commands and rebuilds byte-identical.  Every slab face except the
symmetry plane must be merged, and every volume meshed, before export.
`add_kelvin_cubit` (cubit-mesh-export) was made to select only faces on the
air sphere: it used to take the largest face of every air volume, which with
slabs put interior faces into `kelvin_int`.

`build_mesh_family.py --gap-layers 6 12 24` writes the separately named
manifest `gap_family.json` (schema `c-type-cubit-gap-family.v1`): the gap is
refined alone, the iron, air and Kelvin sizes are held at the medium base
level, and the manifest must not be read by `run_mesh_convergence.py`.
`analyze_gap_family.py` reads one three-engine result per level, requires
HDiv-MMM (iron mesh only) to be identical across levels, and reports the two
FEM routes' gap-core increments, observed order in N and Richardson estimate,
which estimate gap-refinement sensitivity under an asymptotic-convergence
assumption, not a rigorous error bound for the gap or whole model.

`results/lab_20260912_gap_family.json` and the three
`lab_20260912_gap_family_n{06,12,24}_mesh.json` contracts (LAB, curve order
2, Cubit 2025.12):

| N | elements | gap elements | minimum crossings | longest stay / limit |
|---|---|---|---|---|
| 6 | 104,794 | 11,874 | 16 | 1.427 / 2.000 mm |
| 12 | 528,094 | 96,842 | 34 | 0.633 / 1.000 mm |
| 24 | 3,501,632 | 801,988 | 72 | 0.211 / 0.500 mm |

All three pass the exact reflection, periodic Kelvin (462 pairs, trace ratio
1 to 1e-15) and on-sphere `kelvin_int` gates.

### The direct reduced-A solve, and why N=24 first failed

The first attempts at the family stopped in PARDISO: N=6 order 3 and N=12
order 2 exhausted the 57 GB of mdx at about 1.5 M HCurl unknowns, and N=24
order 1 (4.1 M) exhausted the 220 GB of hibino.  The iterative routes do not
apply (AMS is refused on periodic Kelvin HCurl; BDDC at order 1 keeps every
edge in the coarse space).  The cause was the pip NGSolve PARDISO wrapper:
it hard-codes the minimum-degree ordering and ignores the ordering it is
handed, and that ordering fills far more than METIS on three-dimensional
HCurl systems.  Peak process memory of the order-1 reduced-A solve, same
field to 1e-11 (`probe_n24_direct.py`, a research probe kept in `C:/temp`):

| unknowns | shipped default | METIS (+SPD) | sparsecholesky | umfpack |
|---|---|---|---|---|
| 124 k (N=6) | 1.40 GB | 0.81 (0.88) GB | 1.56 GB | 6.7 GB, 47 s |
| 618 k (N=12) | 12.7 GB, 62 s | 5.0 (5.1) GB, 20 s | 16.5 GB, 79 s | fails at 29 GB |
| 4.09 M (N=24) | out of memory, 220 GB | 54 GB, 200 s | MemoryError | -- |

`radia.vector_potential_solver` now registers its own SPD PARDISO with
METIS ordering (`DIRECT_INVERSE_TYPE`) through NGSolve's own hooks and uses
it on its direct paths; the solutions equal the shipped wrapper's to 2e-14
on every level below.  N=24 at order 2 (about 21 M unknowns) is still out
of reach, so the order-2 family stops at N=12.

### Three-engine results on the gap family

hibino, one contract per family: the branch's Python source over the
installed `radia 4.95.90` binaries (`results/hibino_20260912_gap_family_
native_overlay_manifest.json` records their hashes; the driver's
`radia_version` field is the installed distribution's `4.95.90` while the
source is `4.95.91`), linear, HDiv BDM2, direct reduced-A with METIS,
72 threads.  Gap-core relative RMS of the median-plane-projected B:

| level | FEM order | HDiv vs reduced-A | HDiv vs mixed Omega | reduced-A vs mixed | reduced-A unknowns |
|---|---|---|---|---|---|
| N=6 | 1 | 0.200% | 0.192% | 0.226% | 123,777 |
| N=12 | 1 | 0.101% | 0.155% | 0.189% | 617,859 |
| N=24 | 1 | 0.091% | 0.163% | 0.189% | 4,087,820 |
| N=6 | 2 | 0.092% | 0.117% | 0.051% | 544,793 |
| N=12 | 2 | 0.088% | 0.107% | 0.035% | 2,732,075 |

(`results/hibino_20260912_gap_family_*_metis.json`; the `mdx*_20260911_*`
and `hibino_20260912_*_order{1,2}.json` files without `_metis` are the
same levels on the shipped solver, identical fields to 4e-14.)

`results/hibino_20260912_gap_family_report_linear_order1.json` is the
`analyze_gap_family.py` report on the order-1 family.  It passes: HDiv-MMM
is identical across the three levels (it never sees the gap air), and both
FEM routes' gap-core increments contract in N -- reduced-A with ratio 0.31
(observed order 1.7, Richardson estimate 1.9e-4), mixed Omega with ratio
0.47 (order 1.1, Richardson 2.9e-4).  What the two routes converge to
differs.  Reduced-A approaches HDiv-MMM at every level (0.200, 0.101,
0.091%); the mixed Omega route does not (0.192, 0.155, 0.163%) and stays
0.189% from reduced-A at N=12 and N=24.  Its Kelvin source-trace tangential
residual is 4.53% at every level (gate 5%), also unchanged by the gap.
Gap resolution is therefore not what limits the mixed Omega route on this
base mesh; that is consistent with the source-load quadrature localised on
ESRF #6, and is not a proof of it.  At order 2 both routes sit within
0.12% of HDiv-MMM on two levels; three levels are not available.

The report bounds the gap-resolution error of this base mesh only.  It says
nothing about the iron, outer-air or Kelvin discretisation, and nothing
about a family built with other base sizes.

