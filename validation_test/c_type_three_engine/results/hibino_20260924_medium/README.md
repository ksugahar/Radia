# Medium linear fixed-mesh validation

Same wheel, driver, wrapper and 8-thread configuration as the adjacent
hibino_20260924_cached coarse evidence; only checked mesh level changes.
All three engines ran from scratch, foreground, with no other compute jobs.
Linear mu_r=1000; HDiv BDM1 and FEM order 2. All declared core pairwise
field gates passed. This is NOT matched-error performance certification.

HDiv: 3.492863 s, 3804 primary DOF. Reduced-A: 27.9125853 s, 190396
primary DOF. Mixed Omega: 975.0445808 s, 60848 primary DOF. Timings include
build, solve and observation evaluation; mixed Omega includes source caching.
Auxiliary projection DOF and separate per-engine peaks remain unrecorded.
Source cache preparation: 529.2549417 s, 1218000 points; Hodge phase total
530.4572528 s, 6090000 hits, zero misses.

Maximum projected gap-core pairwise relative RMS is 0.412377922% (coarse:
0.453545759%). Coarse-to-medium projected-core increments: HDiv 0.01615%,
reduced-A 0.07186%, mixed Omega 0.13263%. Raw full-tube increments: HDiv
0.01656%, reduced-A 0.97807%, mixed Omega 1.71586%. Medium is not an exact
reference; these increments are not absolute errors or observed orders.
Two levels are insufficient for a convergence certificate. Retain raw
fringe data, and do not generalize the useful-core pass to all points.

Full recovery: S:/Radia/validation_artifacts/hdiv_medium_20260924/recovery.zip
SHA256 81d01706714c9bb14a1c12b7d32f739598a5f16afb553e8a3e6241beb5e35dd8.
LAB and remote archive hashes match. Archive contains scripts, input meshes,
wheel, command provenance, logs and checkpoints. Following this commit,
delete only C:/temp/hdiv-medium-20260924 and its task recovery ZIP; preserve
the shared C:/temp/radia-ctype-family. No manuscript superiority claim yet.
