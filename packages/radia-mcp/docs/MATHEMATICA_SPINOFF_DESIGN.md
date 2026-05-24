# mathematica subpackage spinoff — design discussion

Status: **design discussion only**, no code change yet (2026-05-25).
Owner: Sugahara lab.

## Question

Should `radia_mcp.mathematica` move out of the radia-mcp wheel into
its own PyPI package (`mcp-server-mathematica`)?

Raised in the 2026-05-24 thorough review: "`mathematica` is
essentially Radia-unrelated. Same justification as why
`mcp-server-elf` is a separate repo — wolframscript bridge has
nothing to do with Radia's C++ core or NGSolve integration."

## Current state

- Located at `packages/radia-mcp/src/radia_mcp/mathematica/`
- Promoted 2026-05-20 from `s:/mcp-server/src/mcp_server_document/mathematica/`
- 11 MCP tools: `mathematica_evaluate`, `mathematica_status`,
  `mathematica_simplify`, `mathematica_to_tex`, `mathematica_check_identity`,
  `mathematica_vector_calc`, `mathematica_unit_convert`,
  `mathematica_solve`, `mathematica_integrate`,
  `mathematica_differentiate`, `mathematica_server_status`
- Subprocess wrapper around `wolframscript`
- Cross-references: `differential-forms` (symbolic verification),
  `radia-ngsolve` (formula derivation)
- Standalone — does NOT depend on radia, NGSolve, or any other lab
  code. Pure wolframscript bridge.

## Pros of spinoff (`mcp-server-mathematica`)

1. **Useful to non-Radia users**. Anyone with Mathematica + Claude
   Code would benefit. Independent PyPI = much broader audience than
   the radia-mcp[radia] users.
2. **Lighter dep for radia-mcp**. Currently the mathematica
   subpackage is in the radia-mcp wheel even for users who never use
   Mathematica.
3. **Versioning independence**. Mathematica-specific bug fixes can
   ship without bumping radia-mcp version.
4. **Matches the `mcp-server-elf` precedent** (similar argument:
   ELF600 has nothing to do with Radia internals, so it's a
   separate PyPI dist).

## Cons of spinoff

1. **Cross-references break**. The Kelvin-transform demo (radia-mcp
   README's killer-demo) chains `differential-forms` + `mathematica`
   in one session. Both being in radia-mcp means a single
   `pip install radia-mcp` lights up the demo. If mathematica spins
   off, the demo doc must mention an extra install.
2. **Repository overhead**. Another GitHub repo to maintain CI,
   releases, version drift checks.
3. **Smaller community**. Lab maintains 4 packages already
   (radia-mcp, mcp-server-elf, mcp-server-document, COMSOL fork).
   Adding a 5th = more cognitive load.

## Comparison with mcp-server-elf split

The elf split happened because:
- elf was 1980 lines of dense ELF600-specific knowledge
- elf had its own user base (ELF customers)
- elf was 100% self-contained (no Radia / NGSolve coupling)
- elf grew first, then split (organic separation)

For mathematica:
- 11 tools, ~700 lines of wolframscript wrapping
- No clear "Mathematica user base" outside of Radia for the lab
- Self-contained (yes)
- Promoted inward to radia-mcp 2026-05-20 (only 5 days ago) — has
  not yet had time to grow an identity

## Recommendation

**Stay in radia-mcp for now**. Reasons:

- The recent promotion means the structure is still settling. A
  spinoff this soon is premature.
- The killer-demo chain (`differential-forms` + `mathematica`)
  benefits from co-location in one wheel.
- 11 tools is a small budget; the maintenance overhead of a new
  repo isn't justified yet.
- The `EXTERNAL_PACKAGES` mechanism (added 2026-05-24) handles the
  discovery side — when a user calls `radia_mcp_overview()` they
  see all 37 catalog servers; mathematica is discoverable.

**Revisit if**:
- Mathematica grows past 30 tools / 2000 LOC
- An identifiable non-Radia user base appears (e.g. someone asks
  "can I use just the mathematica tools without Radia?")
- The killer-demo chain is replaced or supplemented (so co-location
  matters less)

## Action

No code change. This design doc is the lock — if a future session
proposes spinoff, point here.
