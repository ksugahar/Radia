# Contributing to radia-mcp

`radia-mcp` is the Sugawara Lab (菅原研究室) MCP suite for the Radia
CAE ecosystem. Contributions are welcome — particularly in the
directions outlined below — but the lab maintains the architectural
conventions described in [README.md](README.md), so please skim that
first.

## Quick orientation

- **Python ≥ 3.10**, pure-Python except for OCCT (via build123d),
  gmsh, and the optional Radia core wheel.
- Source layout:
  ```
  packages/radia-mcp/
    src/radia_mcp/
      cubit/          — mcp-server-cubit
      build123d/      — mcp-server-build123d
      gmsh/           — mcp-server-gmsh
      common/         — failure_log, web_docs, examples (scrape index)
    pyproject.toml    — entry points + extras
    CHANGELOG.md      — append a release entry
  ```
- Primary workflows are `@mcp.tool()` functions with docstrings. Fine-grained
  validation/identity helpers use a server-local `CoarseToolRegistry` and
  `@_validation.tool()` so production clients see one searchable catalog and
  one runner instead of dozens of near-identical schemas. Keep the underlying
  function directly importable for focused tests.
- Do not import optional numerical, RAG/Chroma, or validation-only modules at
  server import time. Use `common.lazy_call.lazy_callable` for call-only helper
  dependencies and add a cold-import validation result when changing a large
  server's dependency graph.
- The `cubit_session` Plan A architecture (PySide6 QTimer + file-drop
  IPC, no sockets) is non-negotiable for the live-Cubit path; if you
  need a different transport, justify it in the PR description.

## High-value contribution areas

### 1. New scrape sub-sources in `radia_mcp.common.examples`

Pattern: write `refresh_<source>(...) -> dict` that fetches, indexes
to `_examples_dir(<source>) / *.md`, and writes the index.json.
Then add to `REFRESH_FUNCS` and the appropriate `FAMILIES` list.
Examples already shipped: GitHub Issues / Discussions (REST + GraphQL),
Discourse forum search, GitLab issues, StackExchange API, YouTube
transcripts via `youtube-transcript-api`. Targets: mailing list
archives (mailman / Pipermail), more YouTube channels, blog posts via
RSS.

### 2. Cookbook KB topics

Both `radia_mcp.build123d.build123d_knowledge._TOPICS` and
`radia_mcp.gmsh.gmsh_knowledge.TOPICS` are dictionaries of markdown strings
keyed by topic name. Add an entry, register it in the dict, and the
unified retrieval picks it up automatically (`build123d_lookup` /
`gmsh_ask`). Topics that would help: more Plane/Axis recipes,
boundary-layer mesh recipes, sideset/nodeset conventions for specific
solvers (Abaqus / Sierra / OpenFOAM).

### 3. Lint rules

`radia_mcp.cubit.rules` (cubit) and `radia_mcp.build123d.rules`
(build123d) both follow the same shape: a function taking
`(filepath, lines: list[str]) -> list[dict]` that returns
`{line, severity, rule, message}` records. Add the function and
append it to `ALL_RULES`. Keep `severity` to one of
`CRITICAL | HIGH | MODERATE | LOW`.

### 4. New starter templates in `generate_build123d_script`

Add a key to `_B3D_TEMPLATES` (raw Python source string) and update
the docstring on `generate_build123d_script`. Lab-relevant patterns
particularly welcome: more Radia magnet topologies (Halbach 3D, ring
quadrupole, …), motor / generator components, accelerator devices.

## Common pitfalls

### `subprocess.Popen` MUST set `stdin=subprocess.DEVNULL` for child GUI apps

When an MCP tool spawns a long-running external GUI process
(Cubit, Gmsh GUI, FreeCAD, etc.) with `subprocess.Popen`, set
**all three** stdio streams explicitly. Specifically `stdin` —
omitting it inherits the MCP server's stdin, which is the
JSON-RPC pipe to Claude Code / the MCP host.

Symptom of the bug: the spawned GUI starts but is `Responding=False`
from the very first poll, threads decline (e.g. 23 → 15 → 12), and
the process never reaches an interactive state. The same binary
launched from PowerShell `Start-Process` works fine, masking the
issue during local manual testing.

```python
# BAD — Cubit inherits MCP stdin (JSON-RPC pipe) and hangs.
proc = subprocess.Popen(
    [exe, "-nojournal", "-input", str(wrapper)],
    stdout=subprocess.DEVNULL,
    stderr=subprocess.DEVNULL,
    creationflags=DETACHED_PROCESS | CREATE_NEW_PROCESS_GROUP,
)

# GOOD — stdin explicitly redirected.
proc = subprocess.Popen(
    [exe, "-nojournal", "-input", str(wrapper)],
    stdin=subprocess.DEVNULL,    # <-- crucial
    stdout=subprocess.DEVNULL,
    stderr=subprocess.DEVNULL,
    creationflags=DETACHED_PROCESS | CREATE_NEW_PROCESS_GROUP,
)
```

This is **not** fixed by `mcp` 1.27.0 (current as of 2026-04-25);
the framework cannot guess which child processes a tool will spawn.
Confirmed root cause of the long-standing "MCP-launched Cubit hangs
intermittently" issue (fix landed in `cubit/server.py` `open_in_cubit`
and `cubit/cubit_session.py` daemon launcher).

The same caution applies whenever a tool runs a long-lived child
that does **not** itself read from MCP stdin. For short-lived
`subprocess.run` of CLI tools that finish in milliseconds, stdin
inheritance usually does not matter, but explicit `stdin=DEVNULL`
is still recommended for clarity.

## Workflow

1. Open an issue first if the change is more than ~50 lines or
   touches the public tool surface — let's align on tool name /
   shape before you write the code.
2. Fork → branch → PR. Normal CI runs compact contracts, tests related to the
   changed package family, and affected server selftests. Use the explicit
   `radia-mcp` full-audit workflow before a release candidate; run only the
   relevant `validation_test/` lane when numerical or runtime evidence is
   needed, and retain its result JSON.
3. Bump the version in `pyproject.toml` and `radia_mcp/__init__.py`,
   add a CHANGELOG entry, and the maintainer will tag + publish.

## Releasing (maintainer notes)

Before publishing, treat a green matrix as a release candidate, not as
operational completion.  The minimum pre-push health evidence is:

- radia-mcp pytest matrix green (latest review evidence: 229 passed)
- package policy lint green
- package version consistency green
- generated `docs/TOOLS.md` drift check green
- top-level pytest collection green

After PyPI upload, the release is not "done" until wheel-installed MCP
entry points smoke and the release-quad deploy checks pass on their
intended machines.  See
`release_workflow(topic="mcp_quality_review")` and
`validation/mcp_quality/release_candidate_review_2026-06-26.json`.

```bash
cd packages/radia-mcp
rm -rf build/ dist/radia_mcp-*
python -m build && python -m build --wheel
PYTHONIOENCODING=utf-8 \
  TWINE_USERNAME=__token__ TWINE_PASSWORD=$PYPI_TOKEN \
  twine upload --disable-progress-bar dist/radia_mcp-X.Y.Z*
```

After upload, follow the monorepo release-quad deployment policy:
LAB and 100号機 remain editable installs, mdx is a PyPI consumer for
`radia` / `cubit-mesh-export` without `radia-mcp`, and hibino is the
PyPI MCP consumer.  Do not claim operational release quality until the
PyPI entry-point smoke and Phase 8/9 machine checks are green.

## License

By contributing you agree your contribution is BSD-3-Clause licensed
(matching the project license).

## Stance reminder

The lab actively pushes **build123d + Cubit** as the primary CAD /
mesh pair. CadQuery, OpenSCAD, FreeCAD adapters are welcome and
maintained as **interop / compat** layers — they are not first-class
authoring targets. FreeCAD specifically is a `friendly` source: the
lab respects the FreeCAD community and maintains first-class interop
through `freecad_to_cubit_hex`, but new lab work is authored in
build123d.
