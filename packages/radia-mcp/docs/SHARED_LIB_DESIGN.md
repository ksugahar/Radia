# Shared lib design — radia-mcp / mcp-server-elf / COMSOL fork

Status: **design discussion**, no code change yet (2026-05-24).
Owner: Sugahara lab.

## Problem statement

3 lab MCP packages are now active:
- `radia-mcp` (37 servers, includes the meta + status/topics factories)
- `mcp-server-elf` (16 tools, standalone)
- `COMSOL_Multiphysics_MCP` fork (~80 tools, fork of upstream)

The 2026-05-24 thorough review identified that the **same MCP
patterns** are being reinvented in each: status tool, topics tool,
markdown prompts loader, meta-overview. radia-mcp `common/` has them
(`status.py`, `topics.py`, `prompts_loader.py`, `async_runner.py`,
`chroma_retriever.py`). elf and comsol-fork don't, and the review
recommends adopting them.

**Question**: do we extract `radia_mcp.common.*` into a shared lib
all 3 packages depend on, or copy-paste / vendor / hand-port?

## Options

### A. Status quo (copy/port manually)

elf + comsol-fork copy code from `radia_mcp.common.status` (~80 LOC)
into their own tree. Each maintains independently.

**Pros**: zero coordination. No new package. No dep version dance.
**Cons**: drift inevitable. Bug fixes / improvements in one don't
propagate. Three implementations of the same thing.

### B. New PyPI package `sugahara-mcp-common`

Extract `status.py + topics.py + prompts_loader.py` to a new wheel.
All 3 repos `pip install sugahara-mcp-common`.

**Pros**: single source of truth. Standard Python dependency model.
**Cons**:
- 4-repo coordination (need to release common before consumers).
- elf + comsol users pay a transitive-dep install cost.
- Bumping the common lib for one consumer forces all to upgrade.
- Versioning policy needed (semver discipline).
- 3 packages' CI / release-quad now has a 4-machine gate.

### C. Vendor radia-mcp.common into elf + comsol-fork via sync script

Each consumer ships a `vendored/sugahara_mcp_common/` directory
that's an exact copy of `radia_mcp/common/{status,topics,prompts_
loader}.py`. A `tools/sync_vendored.py` script downloads the latest
from `radia-mcp` (via GitHub raw URL or local NAS path) on demand.

```python
# tools/sync_vendored.py (pseudo)
SOURCES = [
    ("https://raw.githubusercontent.com/ksugahar/Radia/main/"
     "packages/radia-mcp/src/radia_mcp/common/status.py",
     "src/<pkg>/vendored/status.py"),
    ...
]
for url, dest in SOURCES: download(url) -> dest
```

**Pros**: no extra wheel. No install dep. Each consumer can pin a
specific snapshot. Sync is opt-in (run the script when you want
fresh).
**Cons**:
- Vendoring is unusual in Python (more idiomatic in Go).
- Each consumer has duplicate code in its tree (showing on GitHub
  as ~3x the LOC).
- Vendor copies can still drift if sync is forgotten.

### D. radia-mcp as the lib (consumers `pip install radia-mcp`)

elf + comsol-fork add `radia-mcp>=0.70.0` as a dep and import
`from radia_mcp.common import register_status_tool`.

**Pros**: zero new package. No copy-paste. Single source of truth.
**Cons**:
- radia-mcp wheel pulls in 37 server subpackages + their knowledge
  payloads (~7-8 MB knowledge text) just to use 80 lines of helper.
- The 3 repos become tightly coupled — bumping radia-mcp forces a
  re-test cycle on all 3.
- Conceptually backwards: radia-mcp is the BIG package; comsol fork
  is smaller. Depending on the big from the small is heavy.

## Recommendation

**Option C (vendor + sync script)** with a small tweak:

- Vendor `status.py + topics.py + prompts_loader.py` (~150 LOC
  combined) into `<consumer>/<pkg>/vendored/` in elf and comsol-fork.
- Header comment in each vendored file points at the upstream:
  ```
  # VENDORED from radia-mcp v0.70.0 (commit b36ea560).
  # Source: packages/radia-mcp/src/radia_mcp/common/<file>.py
  # Sync: python tools/sync_vendored.py
  # Do not edit by hand — edits will be overwritten by the next sync.
  ```
- Provide a one-shot `tools/sync_vendored.py` script (~30 LOC, stdlib
  only, uses urllib to fetch from GitHub raw).
- Each repo's CI runs the sync script (or its --check mode) as a
  drift gate, parallel to radia-mcp's TOOLS.md drift gate.

**Why not B (PyPI lib)**: at this scale (3 consumers, 150 LOC), the
overhead of running a 4th release lifecycle is higher than the
sometimes-drift cost of vendoring. If we ever grow to 5+ consumers
or the common lib grows past 500 LOC, reconsider B.

**Why not D (radia-mcp as lib)**: the size asymmetry is wrong. elf
is 16 tools, comsol-fork is ~80, radia-mcp is 281. Depending the
small-on-big locks the small to release-coordinate with the big.

## Action plan (if recommendation accepted)

1. Create `tools/sync_vendored.py` template in `radia-mcp/scripts/`
   (so it's the canonical reference implementation).
2. Apply to elf (lower priority — currently 0 of these patterns in
   use).
3. Apply to comsol-fork (after upstream PR for multilingual lands).
4. Each consumer's CI: `python tools/sync_vendored.py --check`
   (exit 1 if vendored files differ from upstream).

## Alternative: do nothing, accept drift

Given that elf and comsol-fork have a much smaller surface and
their MCP usage patterns differ from radia-mcp's broad-spectrum
discovery problem (elf is essentially one big knowledge dump; comsol
fork is upstream-derived with its own conventions), the case for
heavy unification is weaker than it first looks.

**Pragmatic stance**: copy `status.py` once into each consumer when
they actually need it, accept that they'll drift slightly, and
revisit if maintenance pain shows up.

## Decision

**TBD** — recommend filing as a GitHub Discussion on `Radia` repo
once at least one of elf / comsol actually adopts a status tool,
to gauge real demand vs theoretical.

For now: leave radia-mcp.common as the working reference
implementation. Patterns documented in this file are the contract;
any future consumer of these patterns should ping back here.
