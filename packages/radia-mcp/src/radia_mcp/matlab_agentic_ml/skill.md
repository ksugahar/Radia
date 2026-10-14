---
name: matlab-agentic-ml
description: Use the official MATLAB MCP execution layer for machine-learning and reinforcement-learning workflows, with explicit split, evaluation, release, seed, environment, and result-artifact gates.
---

# MATLAB ML/RL Workflow

Use the official MATLAB MCP server as the execution substrate. This package
adds small, solver-neutral artifact gates; it is not a second MATLAB runtime.

## Public video learning catalogs

For AIcia Solid Project, inventory all public channel tabs rather than only the
Videos tab. The 2026-07-20 snapshot contains 397 distinct IDs: 245 videos, two
Shorts, and 150 streams. Public titles and metadata may nominate a CAE lesson,
but transcripts and media are not copied. A candidate is not numerical ground
truth and cannot be promoted without deterministic seeds, units/schema,
provenance, and an independent forward-solver verification.

## Workflow

1. Pin MATLAB release, toolbox availability, session owner, and random seed.
2. Define the data split or RL environment before training.
3. Train with recorded options and elapsed-time stages.
4. Evaluate on data/episodes disjoint from training and record the evaluation
   metric separately from the training metric.
5. Store a result artifact with schema, units, release, seed, split/environment,
   digest, and one to four dominant timing stages.
6. Replay the artifact through `validate_matlab_ml_rl_artifact` and, for the
   v2 split/evaluation contract, `validate_matlab_ml_rl_v44_identity` before
   treating it as a learning result.

For the v46 generalization contract, also bind nonfinite filtering, worker
seed, restart state, RL episode timeout, environment-defined termination,
training-only exploration, checkpoint lineage, structured tool-argument shape,
existing-session attachment, timeout/error class, and monotone checkpoint
order. Treat any mismatch as `needs_attention` even when the numeric output
looks plausible.

For supervised learning, cross-validation or a holdout is required; a
resubstitution score is not a validation result. For reinforcement learning,
the environment, training episode count, evaluation episode count, and mean
evaluation return are required. Evaluation must use the pinned environment, a
fresh evaluation seed, and no exploration. A training curve alone is
insufficient.

Use MathWorks' official MATLAB MCP/Agentic Toolkit for tool discovery,
execution, testing, and code analysis. Keep public knowledge generic: do not
embed lab paths, licensed solver results, or private MATLAB MCP implementation.
The legacy Gypsilab acoustic lane is replay-only; new MATLAB learning belongs
to this official ML/RL lane.
