# Persisted native TPE history — mdx2, 2026-09-10

## Change and scope

Table/MAT is still authoritative. Eligible persisted scalar TPE studies now
keep the same computation-only native history used by in-memory studies.
Attaching to a restored/populated Study rebuilds completed observations and
distribution/intersection/group metadata once, preserving parameter order.
Rebuild does not consume random numbers. The existing RNG snapshot is restored
separately. No Rust runtime, new MEX ABI, storage schema change, or reduced
AutoSave frequency was introduced.

General-history eligibility guards remain. A public FrozenTrial after_trial
hook invalidates the computation cache rather than duplicating an existing
completed observation. Arbitrary mid-session imports/edits are not claimed to
receive incremental cache repair. Large-history resume latency is not measured.

## Measurements

Same benchmark and mdx2 host as schema-6 baseline commit `27a0e1c82`:
MATLAB R2026a Update 5, seed 37, one float, scalar TPE, StoragePath enabled,
AutoSave=false. Seven timed repetitions, first two discarded. Each trial probe
is followed by table materialization and a full save; profiling is separate.
Trial time below is the **sum of component medians**, not median total latency.

| History | Baseline trial | Candidate trial | Reduction | Baseline dirty tables | Candidate dirty tables | Baseline save | Candidate save |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 100 | 3.568 ms | 3.190 ms | 10.6% | 3.016 ms | 2.964 ms | 0.084 s | 0.080 s |
| 1,000 | 4.236 ms | 3.639 ms | 14.1% | 3.034 ms | 3.140 ms | 0.103 s | 0.102 s |
| 10,000 | 15.520 ms | 12.925 ms | 16.7% | 10.750 ms | 21.247 ms | 0.264 s | 0.324 s |

The 1,000-to-10,000 fill segment (8,986 trials, no per-trial materialize/save)
fell from 61.575 to 42.445 s. This is not the full benchmark elapsed time.
At 10,000 rows, profiling recorded seven sampleNativeGroupedHistory calls and
21 optuna_mex calls over seven trials, confirming the intended cached path.
MEX time includes native work plus its gateway; pure transfer cost is not isolated.

The dirty-table and save measurements worsened at 10,000 rows. Therefore this
is evidence of reduced sampler-loop cost, **not** of faster per-trial table
monitoring, AutoSave throughput, or all-workflow performance. Measurements are
non-interleaved sessions and include runtime variation. Investigate repeated
table materialization and RNG-state table delete/append costs next, and specify
crash recovery before changing full-checkpoint durability. Rustuna was not run.

## Evidence and reproduction

- Candidate CI: https://github.com/ksugahar/Radia/actions/runs/34453133245
  (success; focused normalized-snapshot roundtrip test passed before timing).
- Baseline: `results_optuna_normalized_storage_mdx2_20260910.json`, run
  `34451777406`.
- Candidate raw artifact: `results_optuna_persisted_native_mdx2_20260910.json`.
  SHA256: `ab94b75697209f5ada7ee2374d2f39274a27931bd72344698418dce4ec0907ed`.
- Executed TPESampler.m SHA256 (checked against the service run copy):
  `d9d5e5a63538fe95d9124f5051725b6c01e9a83869d4deb0eaf01369c55c6857`.
- Unchanged Study.m SHA256:
  `1f72b4c826926ec2cf6407606e90793724158370b334ccd6ac5397c7ad5ed20a`.
- Unchanged MEX SHA256:
  `f2bbd68bb3e427ddf75fdc0880d790cd24f44336b237812aba812927358f7c3e`.
- Run `benchmark_optuna_storage_scaling(outputPath,true)` in a dedicated
  Engine on an idle compute host. Use the authenticated mdx service context.
  The experimental service workflow is not a production workflow change.
- The owned MATLAB Engine shut down; no python/MATLAB process remained on mdx2.
- 76 pinned upstream-oracle tests passed locally, including repeated persisted
  resume for intersection transitions and grouped categorical branches.
  Expected proposals remain upstream-generated, not handwritten values.
- 26 table/reliability/core MATLAB tests and 11 standalone package tests passed.

Long performance evidence stays under validation_test; short seeded oracle
checks stay under tests. This change is not a release/merge acceptance claim.
