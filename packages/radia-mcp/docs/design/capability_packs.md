# MCP capability packs

Capability packs combine existing domain tools in one stdio process. The domain
modules remain the implementation owners. Tool names, schemas, annotations,
structured results, images, prompts and static resources retain their original
contracts. Calls use the SDK's public APIs and retain the domain call-log and
reload wrappers. There are no child MCP transports or generic tool-execution
proxies in the client-visible tool list.

## Available entry points

| Command | Profiles (in addition to `all`) |
| --- | --- |
| `mcp-server-radia-design` | `optimization`, `learning` |
| `mcp-server-radia-electrical` | `all` (PEEC, PCB, litz/transmission) |
| `mcp-server-radia-motion` | `motor`, `maglev`, `ih`, `accelerator` |
| `mcp-server-radia-analysis` | `field`, `verification`, `theory` |
| `mcp-server-radia-acoustic-workflows` | `production`, `education` |
| `mcp-server-document-ops` | `inspect`, `convert`, `project` |
| `mcp-server-radia-publication` | `paper`, `poster`, `figures` |

The default is `all` for compatibility with the complete member surface. Use a
narrow profile in everyday sessions, e.g. `mcp-server-radia-motion --profile
maglev`. Selection happens before domain imports. Discovery is fixed at startup;
restart after adding/removing tools or changing a profile. Existing body-only
reload behavior belongs to each domain. A duplicate tool, prompt or resource
identifier aborts startup rather than choosing an implementation silently.

`capability_pack_status` reports the active profile, source, version and exact
tool owners. Source-domain status tools still describe the domain, not the whole
pack. The acoustic production and educational routes remain explicitly named;
neither substitutes for the other. Motion integration groups existing workflows
and verification tools; it does not claim a newly implemented autonomous solver.

## Writing consolidation

`mcp-server-paper-writing` now also registers `grant_writing_*` and `poster_*`,
alongside the already integrated presentation and figure tools. Grant diagnostics
retain read-only annotations, document-type requirements and separate scoring.
Manuscript and funding-review scores must not be averaged. The publication pack
adds bibliography and, in the poster/all profiles, chart tools to that surface.
`references.bib` management and `.bbl` production remain bibliography-owned.

Start a review with `paper_writing_review_route(document_type="paper" |
"grant_proposal", language="ja" | "en", venue_or_program=...)`. Language and
genre are independent choices. Japanese and English papers share scientific
standards but have different prose diagnostics. Grants use funding-review
criteria and proposed-work feasibility. No validated English grant readability
score is claimed. The route never supplies an aggregate score or marks venue
requirements verified; the agent must consult the actual submission rules.

The old `mcp-server-grant-writing` and `mcp-server-poster` commands remain available
as narrow compatibility endpoints. The catalog marks them `integrated_into:
paper-writing`. Other old commands likewise remain valid, with `capability_packs`
listing their replacements. Keep compatibility for at least two published
releases and verify client migration before removing commands.

## Client migration

After installing the wheel, replace selected old client entries with the pack
command and optional `--profile` arguments. For example, replace the six design
server entries with one `mcp-server-radia-design` entry. Do not configure both an
old member and its replacement pack in the same session. Reconnect the client
after editing its configuration. Client configuration changes are separate from
package installation and are not performed by this implementation.

Cubit, build123d, Gmsh, official MATLAB/Optuna, Radia MATLAB, Mathematica and
literature/RAG remain independent process/installation boundaries. Choose only
the capability packs needed for a task. `all` intentionally preserves the large
publication surface; combining servers alone does not reduce schema size.

## Verification

Run `pytest tests/test_capability_packs.py` from the package, plus the affected
writing/catalog tests. Tests compare complete live schemas and cover real stdio
initialization, discovery, domain calls, prompts, resource reads and errors.
Every new command supports `--selftest` (profile import and discovery only).
Wheel checks must run outside the checkout without `PYTHONPATH`; use a separate
environment so another session's editable installation is not redirected.

New pack members must have no server-owned startup/shutdown lifecycle. Resource
templates currently fail explicitly until a routing adapter is supplied. Do not
add session-owning external servers to the manifest. Keep numerical kernels and
independent validation methods in their owning modules.
