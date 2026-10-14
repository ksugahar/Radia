# MCP quality validation

This directory stores public-safe validation summaries for radia-mcp server
fleet quality gates.

These records are intentionally lighter than physics validation scripts: they
capture catalog health, generated documentation drift checks, publish-boundary
lint, per-server selftests, and representative loop-lane regression tests.

Use them as release evidence and as a checklist for future self-driving loop
slots. They must not contain private solver paths, internal campaign paths, or
commercial-tool provenance.

Records:

- `golden_gate_2026-06-25.json` — catalog health, generated TOOLS.md drift,
  public-boundary lint, and representative selftest evidence.
- `release_candidate_review_2026-06-26.json` — review result after the
  radia-mcp matrix reached 229 passed: strong public-facing health, with
  operational quality still gated on PyPI entry-point smoke and release-quad
  machine checks.
