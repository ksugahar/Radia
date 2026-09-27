# C-type source-load quadrature and response order (2026-09-27)

Diagnostic evidence, not three-engine acceptance. Linear C-type model,
`mu_r = 1000`, coarse Cubit Kelvin mesh (`mesh_result_sha256` in every JSON;
24,134 elements), Radia 5.0.1 release tree, NGSolve 6.2.2606, INTEL11,
8 threads, PARDISO. Each JSON records the host, versions and the hash of the
script that produced it; the scripts are committed unmodified so the hashes
still match.

Question: does quadrature of the coil source load near the coil (the
reciprocal of near-field potential evaluation) explain the ~0.15% gap-core
difference between mixed total/reduced Omega and HDiv-MMM?

## Load-only interventions (`ctype_rhs_quadrature.py`, `result.json`)

The stiffness, interface and Kelvin terms stay at production
`bonus_intorder = 4`; only one source load changes. The harness reproduces the
production route to 1.1e-14, and the externally assembled air load (absolute
tetrahedral order 8, the production rule) to 9.8e-15.

| variant | gap-core B change vs production |
|---|---|
| air volume load, order 12 / 16 / 20 | 1.0e-6 / 2.1e-7 / 5.2e-7 |
| air load as enclosing-face flux (`surface_flux`) | 4.3e-7 |
| iron total-Hodge projection bonus 8 / 12 / 16 | 1.2e-14 / 2.0e-15 / 1.8e-14 |

The Hodge projection is B-invariant by construction: a change of `Phi_s` in
the order-p space is absorbed by `phi_total` through the interface jump, so its
quadrature cannot reach B even though its harmonic norm moves
(0.0292141863 at bonus 4, 0.0292143428 at bonus 8).

Element-wise source-integral quadrature error against order 20, by signed
distance `d` from the element to the coil conductor over the longest edge `h`:

| air, production order 8 | elements | median rel. error | share of error^2 |
|---|---|---|---|
| straddles conductor surface | 534 | 6.5e-4 | 52% |
| inside conductor | 58 | 1.8e-4 | 45% |
| d/h < 0.5 | 172 | 4.6e-6 | 4% |
| d/h >= 1 | 13,529 | <= 1e-11 | 0% |

The air error converges only algebraically (global 5.3e-3, 2.9e-3, 1.7e-3,
1.3e-3 at orders 6, 8, 12, 16) because `grad H_s` jumps on the conductor
surface. Iron has no straddling element; its global error is 1.9e-6 at order
8, concentrated in the 104 elements with d/h < 0.5, and converges quickly
(4.4e-8 at order 16). A near-element subdivision would target the ~600 air
elements at the coil, but they do not affect the gap-core field.

## Whole mixed-solver bonus (`ctype_mixed_bonus.py`, `mixed_bonus.json`)

Air load frozen at the production vector; bonus 8 / 12 / 16 for stiffness,
interface, Kelvin and the iron source remainder together: 3.6e-7 / 4.4e-7 /
4.6e-7.

## Exact Kelvin exterior source (`ctype_exact_exterior.py`, `exact_exterior.json`)

`exact_exterior_source=True` (pulled-back `KelvinRadiaFieldStrength`) instead
of the projected `kelvin_int` trace: 2.8e-7.

## Response order (`ctype_order3_pair.py`, `order3_pair.json`)

All engines recomputed on 5.0.1 (the stored `20260908` fields came from an
older release; drift 1.7e-4 HDiv, 2.4e-4 reduced-A). `vim.Solve` accepts only
BDM1/BDM2, so HDiv-MMM stays at BDM2.

| gap-core relative RMS | p = 2 | p = 3 |
|---|---|---|
| reduced-A vs mixed Omega | 0.162% | 0.052% |
| HDiv BDM2 vs reduced-A | 0.102% | 0.091% |
| HDiv BDM2 vs mixed Omega | 0.152% | 0.122% |

Mixed Omega moves 0.205% from p = 2 to 3 and reduced-A 0.049%. The two FEM
routes contract towards each other by a factor of about three, and at p = 3
HDiv BDM2 is the outlier by about 0.1%.

## Near-coil field with near-element subdivision (`ctype_near_coil_subdivision.py`, `near_coil_subdivision.json`)

The air load keeps the production order-8 rule on far elements and, on the
1,325 near elements (534 straddling the conductor, 58 inside it, the rest
within one element size), splits the reference tetrahedron into 8^L Bey
children with the order-8 rule on each.  Near + far with the unchanged rule
reproduces the whole load to 1.3e-16.  B is sampled on three lines through
the +x straight leg (z = 0 and z = 0.06 m across x, and along z at the
conductor centre) and at the gap core; the reference is L = 2.

| air load | near-coil lines, relative RMS | gap core | load assembly |
|---|---|---|---|
| production, order 8 everywhere | 3.2e-4 to 4.7e-4 (max 9.6e-6 T) | 4.3e-7 | 114 s |
| near elements split L = 1 | 3.7e-5 to 7.0e-5 | 9.0e-8 | + 114 s |
| order 20 everywhere | 2.9e-5 to 4.9e-5 | 8.5e-8 | 1,486 s |

Splitting only the near elements once gives comparable order-20 accuracy
on these near-coil lines. The reported 114 s is the near L1 assembly only;
it excludes the far-element assembly and is not the complete load or solve
cost. The 1/13 comparison therefore applies only to that partial stage.
Assembly times were recorded in the run log, not the committed JSON.
The gap-core field changes by less than 1e-6 relative.  The largest pointwise relative values (up to 2.2e-2) sit where
|B| is nearly zero and are 9e-6 T in absolute terms.

## Conclusion

On this model the ~0.15% mixed-Omega discrepancy is response discretization
(dominantly mixed Omega at p = 2), not source-load quadrature, source
projection or the Kelvin trace, all of which move gap-core B by less than
1e-6. The remaining ~0.1% at p = 3 cannot be assigned uniquely to HDiv BDM2
without further independent mesh/order convergence evidence. One coarse mesh only; this is not a mesh-convergence
certificate.

## Reproduction

The scripts expect one staging directory holding `mesh/` (`kelvin_domain.vol`,
`iron.vol`, `mesh_result.json` of the coarse family),
`validation_test/c_type_three_engine/run_three_engine.py`, and
`results/c_type_20260908_linear_kelvin_after_coarse_mdx1.json`.
`ctype_rhs_quadrature.py` and `ctype_mixed_bonus.py` take `--mesh-dir`,
`--three-engine`, `--reference`, `--output`; run them with
`--inverse pardiso` as recorded (their `sparsecholesky` default is not usable
for this saddle-point system: it returned an all-NaN field on mdx2).
`ctype_exact_exterior.py` and `ctype_order3_pair.py` read the staging
directory next to the script and take the thread count as the only argument.
`ctype_near_coil_subdivision.py` takes `--staging` (the same directory, which
must also hold `ctype_rhs_quadrature.py`) and `--output`.
