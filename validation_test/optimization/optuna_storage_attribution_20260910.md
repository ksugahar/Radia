# Optuna table/storage attribution — mdx2, 2026-09-10

## Finding

MATLAB table is not an inherent ceiling on sampler performance. The current
implementation has two separable costs: storage-path eligibility disables the
incremental native history path even with AutoSave=false, and full persistence
serializes a large number of nested table objects and validates two copies.
No production behavior was changed by this investigation.

## Measurements

MATLAB R2026a Update 5, seed 37, scalar TPE with one float parameter,
100/1,000/10,000 completed trials. Both modes use AutoSave=false. Each level
has seven timed probes; the first two are discarded. The trial-time columns
below sum the separate ask/suggest/tell medians (not the median of their sum).
Table materialization follows every timed trial; these are diagnostic probes,
not the earlier uninterrupted 100-trial throughput benchmark.

| Initial history | Trial, no storage path | Trial, storage path | Dirty tables, storage path | Full save |
|---:|---:|---:|---:|---:|
| 100 | 1.84 ms | 3.88 ms | 4.91 ms | 0.208 s |
| 1,000 | 2.19 ms | 5.26 ms | 3.51 ms | 1.328 s |
| 10,000 | 10.11 ms | 18.73 ms | 11.01 ms | 9.902 s |

At 10,000 trials, reading already-cached tables takes approximately 7 us in
the storage-path mode. The MAT file is only 2,508,755 bytes. File size alone
does not explain the save latency. No OS-level disk attribution was performed.

## Source and profiler evidence

1. `TPESampler.canUseNativeGroupedHistory` requires an empty StoragePath;
   `attach` sets NativeHistoryValid=false for a nonempty path. The persisted
   route also records RNG state in SamplerStateTable. Therefore the path-only
   comparison includes both cache eligibility and state bookkeeping, not
   cache eligibility alone. It does not measure disk writes in ask/tell.
2. `Study.save` serializes StudyData, validates the temporary MAT, atomically
   replaces the primary file, then copies and validates the backup. The
   10,000-trial save profile records 40,096 table.saveobj calls and 20,048
   table.loadobj calls. TrialTable contains an IntermediateValues table in
   every trial cell, including empty tables, while the normalized
   IntermediateTable is also stored. This creates per-trial object work.
3. `Study.ask`, `recordParameter`, and `finishTrial` each invoke save when
   AutoSave and StoragePath are enabled. Thus an ordinary one-parameter
   ask/suggest/tell loop requests three saves per trial. This is source
   evidence, not an end-to-end AutoSave timing; do not present 3x9.9 s as
   a directly measured trial latency. Multiple parameters/report calls add
   more potential persistence operations.
4. Without a storage path, the 10,000-trial profile confirms the native-history
   route: optuna_mex has 14 calls across seven trials (proposal and history
   update). Its 39 ms profiled total includes gateway and native work. Pure
   transfer time is not isolated, and profiler timings are not wall timings.

## Recommended implementation order (not implemented here)

1. Keep the public table API and authoritative table/MAT history, but normalize
   the durable representation so empty/nested per-trial tables are not
   redundantly serialized. Preserve reconstruction, schema migration,
   attributes, intermediate results, and corruption detection in tests.
2. Separate mutation durability from full checkpoints. Do not simply disable
   saves or reduce frequency silently: specify the crash-recovery contract
   first. An incremental durable record plus checked table/MAT checkpoints
   is a candidate, not a selected implementation.
3. Permit a rebuildable native computation cache for eligible persisted
   studies, with explicit invalidation/rebuild on restore/import/state edits.
   Preserve RNG state and proposal sequence; removing the path guard alone
   is not a sufficient fix.
4. Then investigate long-history native computation and array growth/copying.
   Keep pure transfer attribution separate until measured with a suitable
   narrow native diagnostic. Rustuna was not run in this investigation.

## Reproduction and evidence

Run `benchmark_optuna_storage_scaling(outputPath,true)` and
`benchmark_optuna_storage_scaling(outputPath,false)` in separate dedicated
MATLAB Engine sessions on an otherwise idle host. Use the bundled Engine in
the authenticated mdx service context; the SSH startup failure is separate.
All raw timings and top-40 inclusive/self-time profiles are in
`results_optuna_storage_attribution_mdx2_20260910.json`.
Both CI runs succeeded: `34445398310` and `34445704876`. Owned MATLAB sessions
were closed and no Python/MATLAB process remained after validation.
