# Single-layer Y-URN reduction sensitivity

## Scope

The current MCP operating manual is `urn(topic="overview")` / `urn(topic="method")`.
This record is numerical evidence, not a second operating manual.

The private SA/RM PCB and NL87 34-basis and compact 10-basis log-component
checkpoints were replayed on mdx2, CPU, one thread. Saved predictions reproduced
to maximum relative differences below 4.1e-15. Each curve has 2001 samples.
No raw curves, checkpoints or fitted parameter arrays are published here.

Reduction starts from each saved 10-basis checkpoint, not a new 34-basis fit.
Each candidate receives 6000 Adam epochs and one restart. The path is greedy,
down through one basis; it does not enumerate all subsets. All samples are
training samples. No held-out, repeated-start or physical-uniqueness claim is made.

## What changed after the first numerical run

The initial reducer ranked removals using raw impedance output-ablation
importance. That can hide component errors. The corrected reducer minimizes
unrefitted error against measurements in the configured S/log-component metric,
within the same-family duplicate pool when one exists. It then refits the
selected reduced structure. A regression test makes the raw-magnitude ordering
deliberately wrong and verifies the configured fitting metric is used.

The initial full run is retained in `yurn_reduction_initial_20260914.json`;
the corrected run is `yurn_reduction_sensitivity_20260914.json`. Both carry
numerical source hashes, private input hashes, runtime identity and full traces.
The corrected run took 304.6 seconds (PCB) and 274.5 seconds (NL87).

## Results: component tradeoffs remain visible

| Dataset | Bases | S-domain RMSE | Log-component RMSE | Maximum relative complex error |
|---|---:|---:|---:|---:|
| PCB | 10 | 0.00404969 | 0.0435159 | 18.7285% |
| PCB | 9 | 0.00407201 | 0.0439035 | 18.7416% |
| PCB | 8 | 0.00405167 | 0.0452624 | 18.8208% |
| NL87 | 10 | 0.00155136 | 0.0916337 | 3.91294% |
| NL87 | 9 | 0.00174472 | 0.0814743 | 4.00802% |
| NL87 | 8 | 0.000949545 | 0.207068 | 4.08309% |

PCB's saved 10-basis model contains a parallel-RLC pair with normalized complex
response distance 1.27e-7. Removing one and refitting largely preserves both
metrics. Eight bases remain close in S-domain error, with about 4% larger
log-component RMSE. This supports investigating reduction below ten, not claiming
that eight is the unique physical circuit.

NL87's corrected nine-basis path improves log-component RMSE while slightly
worsening maximum complex error. Eight bases improve S-domain RMSE but worsen
log-component RMSE. Neither metric should be hidden to promote a smaller model.

## Assumed error-budget sensitivity, not instrument calibration

The scenarios use `uncertainty_ohm = relative_budget * abs(Z_measured)` at every
sample, with acceptance requiring every complex error to be within that budget.
These are explicit assumptions, not measured instrument precision or confidence
intervals. Their purpose is to quantify what simplification costs.

| Assumed relative budget | Smallest tested PCB candidate | Smallest tested NL87 candidate |
|---|---:|---:|
| 0.5%, 1%, 2% | none | none |
| 5% | none | 8 |
| 10% | none | 5 |
| 20% | 4 | 5 |

Under 5% and 10%, the initial NL87 path retained ten; the corrected path finds
eight and five respectively. This is evidence that deletion order matters,
not proof of global optimality. PCB's four-basis candidate at the loose 20%
complex budget has much worse log-component error (0.4576); it must not be
presented as preserving loss-component accuracy. The roughly 19% PCB maximum
residual is not attributed to measurement noise without separate evidence.

The trace's `accepted` field refers to the reference **1%** budget; the
`scenarios` field reevaluates that same candidate path at the other budgets.
Parameter fitting and model selection are separate: better fitted likelihood
alone would not justify extra bases. A calibrated maximum-likelihood claim
requires a declared measurement-error model, which this experiment does not
provide. Parsimony is nevertheless an explicit, useful objective.

## Verification

28 focused URN/MCP tests passed locally with the pure-Python URN package loaded
without the unrelated native Radia root. This is not a native solver validation.
The actual numerical computation uses the checked source modules and an isolated
mdx2 environment (versions in the JSON); no shared editable install is changed.
Private inputs and operational evidence are retained on LAB, with SHA-256
verification before removal of the task-owned remote scratch directory.
