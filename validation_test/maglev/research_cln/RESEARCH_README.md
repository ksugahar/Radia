# radia-maglev CLN validation corpus

This directory preserves numerical verification sweeps, Mathematica
derivations, iteration histories, and reference probes that support the CLN
theory. These files belong to `validation_test/` because they are durable
engineering evidence; this directory is not a scratch area or a package-demo
tier.

Fast implementation regressions belong in `tests/`. Human-facing explanations
belong in executed notebooks under `docs/`. New exploratory work stays in
`C:\temp` until it is distilled into one of those maintained destinations.

## Contents

| Subdir | What |
|--------|------|
| `cln/` | The absorbed radia-cln validation corpus (CLN = Cauer Ladder Network for eddy-current transients). See `cln/README.md`. |

`cln/` was migrated here from `examples/cln/` (2026-06-14) when the
radia-cln scope was absorbed into radia-maglev: it is eddy-current
CLN research, the natural domain of CLN, but the bulk of it (NGSolve
validation sweeps, Schur-F derivations, BEM-CLN iteration history) is
durable validation evidence rather than a public demo.
