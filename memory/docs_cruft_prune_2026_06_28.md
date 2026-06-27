# Docs Cruft Prune 2026-06-28

## Lesson

Bulk promotion ledgers and copied legacy source trees do not belong in `docs/`
once their lessons are captured elsewhere. `docs/` should show maintained
theory, result-bearing notebooks, synchronized JSON, and canonical APIs; it
should not become a warehouse for development chronology.

## Pruned

- `docs/examples_consolidation/`: batch-by-batch migration scratch space for
  "next 100" example sweeps. Keep those decisions in memory notes and topic
  notebooks, not as public docs.
- `docs/stream_function/stream_function_examples_archive.*`: source-only topic
  archive that preserved deleted TODO benchmark stubs. The maintained docs now
  live in `docs/stream_function/{theory,regularization,deformation,benchmarks}`.
- `docs/hdiv_vim/vim_examples_archive.*`: full-source inventory of the
  `examples/vim` corpus. The public docs keep `README.md`, productionization
  notes, and result-bearing showcase notebooks; the validation corpus remains
  in `examples/vim` + `validation_test/feec`.
- Topic `*_examples_archive.*` triples for topics that already have maintained
  result-bearing notebooks or validation surfaces. Public docs should point at
  those maintained artifacts, not source-only archive ledgers.
- `docs/kelvin/legacy_assets/kelvin_transformation/`: old examples mirror and
  debug notes. The one live document, the Kelvin convention, was promoted to
  `docs/kelvin/CONVENTION.md`; old source-level history remains recoverable
  from git if needed.

## Rule

When a docs artifact is only a source archive, migration ledger, failed
attempt log, or old path mirror, prune it before promoting more examples.
Promote the distilled rule/API/result instead.
