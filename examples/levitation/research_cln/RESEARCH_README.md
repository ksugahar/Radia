# research/ — radia-levitation research corpus

This directory holds **research scratch**: raw verification sweeps,
Mathematica derivations, iteration history, and one-off probe scripts
that back the theory but are **not** README-backed package examples.

Per the lab Sample Promotion Ladder (tests -> examples -> panels), a
sample lives in exactly one tier:

- `tests/` — golden fixtures (machine-readable, minimal)
- `examples/` — README-backed, standalone-runnable demonstrations
- `research/` — **this dir**: exploratory corpus, no standalone-run or
  README guarantee

A script earns promotion `research/ -> examples/` only after it gains a
README, runs standalone (`python <script>.py` with no hidden state), and
demonstrates one clear concept.

## Contents

| Subdir | What |
|--------|------|
| `cln/` | The absorbed radia-cln research corpus (CLN = Cauer Ladder Network for eddy-current transients). See `cln/README.md`. |

`cln/` was migrated here from `examples/cln/` (2026-06-14) when the
radia-cln scope was absorbed into radia-levitation: it is eddy-current
CLN research, the natural domain of CLN, but the bulk of it (NGSolve
validation sweeps, Schur-F derivations, BEM-CLN iteration history) is
verification scratch rather than curated examples.
