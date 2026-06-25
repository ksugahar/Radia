# MCP quality validation

This directory stores public-safe validation summaries for radia-mcp server
fleet quality gates.

These records are intentionally lighter than physics validation scripts: they
capture catalog health, generated documentation drift checks, publish-boundary
lint, per-server selftests, and representative loop-lane regression tests.

Use them as release evidence and as a checklist for future self-driving loop
slots. They must not contain private solver paths, internal campaign paths, or
commercial-tool provenance.
