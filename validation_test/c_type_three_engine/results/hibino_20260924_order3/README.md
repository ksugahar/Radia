# Coarse-mesh, third-order FEM three-engine diagnostic

2026-09-24, hibino, eight threads, retained Radia 5.0.0 wheel and NGSolve
6.2.2606. Same checked coarse mesh, coil, linear relative permeability 1000
and observation stencil as earlier runs. HDiv BDM2; both FEM routes order
three, mixed source projection order three, bonus quadrature four.

All three engines completed, the fixed-mesh 1% projected gap-core pairwise
gate passed. This is NOT an absolute accuracy or matched-error speed certificate.

| engine | main DOF | build/solve/field seconds | sampled process RSS maximum bytes |
|---|---:|---:|---:|
| HDiv | 8784 | 14.0756531 | 334524416 |
| reduced-A | 345810 | 79.1419698 | 3162120192 |
| mixed Omega | 128071 | 1259.9262756 | 4653805568 |

HDiv: 41 iterations, native convergence true, residual 9.654505507e-9;
Gram 10.8830043 s, H-matrix max rank 294, compression 0.6282090,
H-matrix storage 64.9068 MiB. RSS is whole process sampled every 50 ms,
can retain earlier allocations and miss peaks; not isolated engine memory.
Mixed auxiliary Hodge and interface DOFs are 34919 and 2657, not simultaneous
main-system unknowns.

Projected core differences: HDiv/A 0.09135337%, HDiv/Omega 0.12160546%,
A/Omega 0.05177438%. Raw and full-fringe metrics are retained in the JSON.
Kelvin trace residual 1.62912846%, below the existing 5% trace gate;
iron harmonic norm 0.61028360% is not a pure projection error.

Source cache: 1478304 points, preparation 635.2897611 s (included in total),
Hodge inclusive 728.1971906 s. Hodge reported 4434912 hits, 13688 misses,
hit rate 99.6923077%. Unlike order two this is not a zero-miss cache.
Actual re-evaluation remains in the measured time; cache correctness is not
inferred from an assumed 100% hit rate.

Same-coarse-mesh p2-to-p3 field changes: A 0.04917762%, Omega 0.20487956%.
Coarse p3 versus finer p2: A 0.04085300%, Omega 0.09354510%.
These use each comparator's projected-core vector norm as denominator.
The postprocessor verifies common geometry authority, coil and observation
points. It does not certify absolute error. p-refinement also changes
quadrature; this does not isolate one error source. The prior coarse
full-three HDiv was BDM1, so do not call the pairwise change a controlled
HDiv comparison. The finer result also has a different HDiv mesh.

Recovery on LAB:
S:/Radia/validation_artifacts/hdiv_order3_20260924/hdiv-order3-20260924-recovery.zip
SHA256 f885e8ff245b576127776b6ce900171be03fc38dbe0a54470cce7f2f0425b1a2
matched remote and LAB. Archive includes logs, scripts, wheel and checked inputs.
Wheel SHA256 41b5f5ec02c5c26ad47147da7e4496525e3e2066216d85d787487ecc253023e6.
No solver source was changed, no detached launch, no other job stopped.

Next: retain speed-claim caveats; select a controlled refinement family with
fixed HDiv order and source/response compatibility before equal-error timing.
