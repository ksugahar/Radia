# Fine linear three-engine run

All three engines completed on hibino, Python 3.12.10, Radia 5.0.0 wheel, NGSolve 6.2.2606, eight threads, same settings as coarse/medium. Exact-point source caching remains inside mixed Omega total runtime. The added auxiliary counters and sampled phase RSS do not change equations; measurement overhead is included.

| Engine | Primary DOFs | Build/solve/observation seconds |
|---|---:|---:|
| HDiv BDM1 | 5766 | 4.9254894 |
| reduced-A | 288902 | 40.8551605 |
| mixed total/reduced Omega | 90863 | 1421.5401162 |

All-pair projected core agreement passes the 1% gate; the maximum is HDiv/mixed Omega, 0.8998604033%. This is fixed-mesh agreement, not an accuracy or algebraic convergence certificate.

The three-level increments expose a warning: HDiv core medium-to-fine change is 1.00149%, versus coarse-to-medium 0.016278%, both normalized by fine. The contraction ratio is 61.52, not convergence. Reduced-A and Omega ratios are 0.6224 and 0.8116. Investigate the HDiv fine-mesh change before further scaling or asserting matched-error speed superiority.

Mixed Omega actual auxiliary counts: iron Hodge 23033, Kelvin trace 2966. These are separate spaces, not a simultaneous linear system dimension. Source cache preparation 798.0158 s, inclusive Hodge 799.8329 s; 9233750 hits and zero misses during Hodge. This does not imply every later source evaluation is cached.

Phase sampled whole-process RSS maxima: HDiv 246497280 bytes, reduced-A 2233356288 bytes, Omega 2140250112 bytes. These may include retained earlier-engine allocations and miss transient peaks. They are NOT isolated per-engine peak memory. Sampling errors were empty. Algebraic residual is not provided by the driver.

Recovery: S:/Radia/validation_artifacts/hdiv_fine_20260924/recovery.zip; remote and LAB SHA256 verified as b3c449dd0f3a999a578f0f7d838bf6318848e1564cef5a3a4fa2162d8f911ce3. Archive includes inputs, wheel, scripts, checkpoints, events and run.log. Command/runtime identity is in run_manifest.json. Remote task scratch cleanup follows evidence commit; shared mesh family is preserved.

Commercial software was not run. No manuscript speed-superiority claim is justified by this result.
