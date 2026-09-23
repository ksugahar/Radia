# Finer three-engine comparison, hibino 2026-09-24

Completed foreground SSH run, exit 0; all three mandatory engines pass the
declared 1% parity-projected useful-gap-core comparison and trace contracts.
This is fixed-case cross-formulation agreement, not a continuum-error
certificate or matched-error performance claim. No commercial solver ran.

Same checked physical geometry, coil, linear relative permeability 1000,
observation points, eight threads, and construction/solve/field timing scope.
HDiv uses BDM2 and Gram tolerance 1e-14; both FEM routes use order 2.
Native HDiv convergence is true with residual 9.530603207840468e-9.

| Engine | Total seconds | Sampled whole-process RSS maximum, bytes |
|---|---:|---:|
| HDiv | 59.5483153 | 1014005760 |
| reduced-A | 62.6141462 | 3728904192 |
| mixed total/reduced Omega | 2246.6865511 | 3959480320 |

HDiv/Omega projected core relative RMS difference is 0.002477312180679713;
maximum three-way difference is 0.0026417559039416225. Raw/full-fringe and
reflection diagnostics remain in the complete result. Memory sampling is
whole-process, may retain prior allocations, and can miss peaks.

Omega exact source-cache preparation takes 1268.6466925 s and is INCLUDED in
the total; Hodge projection inclusive time is 1271.6526311 s, with 14677500
hits and zero misses. Auxiliary Hodge and interface spaces have 35999 and
4546 DOFs respectively; these are not simultaneous-system DOFs.

Runtime/import provenance, complete field arrays, mesh hashes, checkpoint
contracts, execution command, driver and instrumentation are retained here.
The wheel identifies native source 1b9f6769bbd402a9859e0a892d2a729e13294d6c;
its SHA256 is 41b5f5ec02c5c26ad47147da7e4496525e3e2066216d85d787487ecc253023e6.

Full inputs, wheel and raw log are recovered outside tracked source at
S:/Radia/validation_artifacts/hdiv_fullfiner_20260924/. Recovery ZIP SHA256
5c7c569fcb91df97c3aaf96f98ffbe8d0cb6ad44cbd3b0bf852e591d35cab12b
matches the remote archive. Cleanup follows this evidence commit; the
shared C:/temp/radia-ctype-family remains untouched.

The prior BDM2 refinement diagnostics show contraction at the last step,
but independent-host replication and comparable error control remain open.
Do not advertise the observed time ratio as a matched-accuracy speedup.
