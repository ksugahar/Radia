# Sampler-state row update — mdx2, 2026-09-10

## Implementation

`Study.recordSamplerState` previously removed all rows matching sampler/schema
and appended a replacement on every update. The usual case has one matching
row at the end of the table. That case now replaces the row in place, avoiding
the deletion and regrowth of seven table columns. Non-tail matches and duplicate
matches retain remove-and-append behavior. Revision, generation, trial number,
timestamp, payload, schema and public row ordering are unchanged.

There is no new state cache, native ABI, storage format or dependency. AutoSave
frequency, exact RNG snapshots, read-back validation, atomic replacement and
verified backup are unchanged. Table remains the authoritative public storage.

The new short MATLAB-only integration test checks tail replacement, previously
returned snapshot independence, multiple sampler/schema keys, public ordering,
revision/generation updates, and exact table equality after save/reload. Its
classification is recorded in the checked test manifest; it is not an upstream
behavioral oracle. Existing seeded upstream comparisons remain the algorithm gate.

## Measurement conditions

`benchmark_optuna_storage_scaling(outputPath,true)`, MATLAB R2026a Update 5,
mdx2 authenticated service Engine, seed 37, scalar TPE with one float,
StoragePath enabled, AutoSave=false. Each history level has seven probes with
the first two discarded; profiles are separate. Trial figures sum the
ask/suggest/tell component medians rather than measuring a median total.

Candidate run: https://github.com/ksugahar/Radia/actions/runs/34453844259
(success, including the focused normalized-snapshot storage test).
Raw result: `results_optuna_state_row_mdx2_20260910.json`.
Its SHA256 is `88756763d6f143eb534e7365e097e82042da4e05812b6af0a0196cb82501de28`.
Candidate Study.m SHA256:
`eb503047c2a8fd51a6c453d312591259a686d7e3e8c741aa2d3df0f3c918cd8f`.
Baseline Study.m SHA256:
`1f72b4c826926ec2cf6407606e90793724158370b334ccd6ac5397c7ad5ed20a`.
TPESampler.m and MEX are unchanged from commit `3c6ea278c`.

The original baseline is run `34453133245`, retained in
`results_optuna_persisted_native_mdx2_20260910.json`.

| History | Original baseline trial | Candidate trial | Candidate dirty tables | Candidate save |
|---:|---:|---:|---:|---:|
| 100 | 3.190 ms | 2.828 ms | 3.013 ms | 0.083 s |
| 1,000 | 3.639 ms | 3.310 ms | 3.087 ms | 0.106 s |
| 10,000 | 12.925 ms | 12.721 ms | 11.060 ms | 0.299 s |

The 8,986-trial fill segment from the 1,000-row probes to 10,000 rows took
39.355 s versus the original baseline's 42.445 s. No table materialization or
disk save occurs in this fill loop. This does not measure AutoSave throughput.
Profiled recordSamplerState inclusive time over seven trials decreased from
0.1008 to 0.0759 s; these instrumented values are not wall-clock trial timings.

The 10,000-row trial-probe difference is small. Runtime variation and repeated
baseline evidence must be considered; do not extrapolate a universal speedup.
Dirty-table timing returned near the earlier schema-6 baseline without any
table-construction change, so the previous 21 ms observation alone does not
justify another cache layer. Pure transfer cost, peak memory, large-history
resume, and end-to-end AutoSave performance are not measured here.

## Reverse-order baseline check and acceptance

After the candidate, the unchanged baseline Study.m was staged again and
verified by its SHA256 in the actual service run. Run
https://github.com/ksugahar/Radia/actions/runs/34454038390 passed and produced
`results_optuna_state_row_baseline_mdx2_20260910.json`.
Its trial component sums were 3.047/3.645/12.901 ms at
100/1,000/10,000 rows. The long fill took 42.209 s, compared with the
candidate's 39.355 s: a 6.8% observed reduction. The earlier baseline took
42.445 s. Thus the baseline/candidate/baseline sessions support a modest
loop-cost benefit, not a large or universal speedup. At 10,000 rows the
individual trial probe only improves about 1.4% over the repeated baseline.

The repeated baseline's dirty-table and save medians were 10.900 ms and
0.300 s, essentially the candidate's 11.060 ms and 0.299 s. Those changes
are not credited to this optimization. Its recordSamplerState profile was
0.0993 s including seven table deletions (0.0248 s); the candidate removes
that deletion work. No algorithm or save-frequency change is involved.

Acceptance: 76 upstream-oracle MATLAB tests, 27 table/reliability/core MATLAB
tests, and 11 package Python tests passed. Both mdx2 runs passed their focused
storage test and shut down their owned Engine; no Python/MATLAB process
remained. The staging source was restored to the validated candidate after
the reverse check. No production release or merge was performed.
