# Cached source three-engine linear smoke

All three fixed-mesh field agreement gates passed. This is not matched-error
performance certification or evidence about commercial software performance.
Linear mu_r=1000; 8 configured threads; same coil and observation points;
HDiv BDM1, FEM order 2; exact point-value caching inside mixed Omega timing.

| Engine | DOF | Build + solve + field evaluation (s) |
|---|---:|---:|
| HDiv-MMM | 3090 | 3.1265337 |
| reduced-A | 126834 | 19.8621060 |
| mixed total/reduced Omega | 41305 | 671.0836538 |

Maximum pairwise median-projected gap-core relative RMS difference:
0.004535457592290168 (0.453546%). This is inter-method agreement, not an
absolute field error certificate. Global process peak: 939.875 MiB; not
per-engine memory. HDiv reports 40 linear iterations; no final algebraic
residual is emitted by this driver, so the passing field gate is not a
separate residual audit. No nonlinear convergence claim applies.

## Accuracy scope audit

The declared pass is ONLY the symmetry-projected useful gap core. Raw full
81-point tube discrepancies are 1.29777% (HDiv/reduced-A), 1.85978%
(HDiv/mixed Omega), and 2.03945% (reduced-A/mixed Omega). Do not describe
this run as below 1% everywhere. Source-trace diagnostics (2.92142% iron
harmonic ratio; 3.96789% Kelvin tangential residual) are distinct from B
accuracy and use their own 5% contract. Off-plane reflection defects are
about 2e-10; the reduced-A one-sided median-plane tangential trace is a
different diagnostic, not a reflection-symmetry failure.

Reported DOF are the primary spaces returned by each engine, not a sum of
every auxiliary projection space. A performance table calling these total
workflow DOF would need additional source-projection space accounting.

Before a matched-error claim, refine both discretizations against a shared
reference/declared observable, report raw and projected core/fringe errors,
explicit algebraic residuals, primary and auxiliary DOF, per-engine peak
memory, and repeated timings including source preparation. Do not substitute
pairwise agreement for mesh convergence or use a cross-host old result as
the performance baseline.

855500 exact source points prepared in 364.3916166 s. Hodge projection
including preparation: 365.2732088 s; 4277500 cache hits, zero misses during
that phase. The remaining mixed Omega phases remain in the total time.
The source wrapper changes no quadrature, material, solver or geometry.
It is an experimental implementation, not a production release change.
The wrapper patches kelvin_solver.project_source_total_hodge, imported
locally by static_electromagnet; the initial incorrect patch target failed
before numerical work and its log is in the recovery archive.

Wheel source: 1b9f6769bbd402a9859e0a892d2a729e13294d6c, native CI 35863320903.
Wheel SHA256: 41b5f5ec02c5c26ad47147da7e4496525e3e2066216d85d787487ecc253023e6.
Driver SHA256: 7baa0feaee8bf3c9b921944d4b5c4102dec0ac072748f0c29e287f5fb1f0c334.
Python 3.12.10, NGSolve 6.2.2606, Radia 5.0.0 in a task-owned venv;
system Radia 4.95.90 was preserved. HIBINO: 38 logical CPUs, 61680418816 B RAM.
OMP_NUM_THREADS, MKL_NUM_THREADS, OPENBLAS_NUM_THREADS were all 8.
Foreground invocation: venv/Scripts/python.exe -u run_cached_three.py.

Full recovery (inputs, logs, scripts, wheel, checkpoints):
S:/Radia/validation_artifacts/hdiv_mixed_20260924/recovery.zip
SHA256: 44ea1fc3d3ff002680f2b629cac790563d5075b0ee5caa28aa50560ae1d1c636.
Remote/source and LAB/archive hashes match. Remote task scratch cleanup
follows this evidence commit; shared C:/temp/radia-ctype-family is preserved.

Next: field-error convergence and repeated common-scope timing, plus explicit
linear residual and per-engine memory records, before manuscript speed claims.
