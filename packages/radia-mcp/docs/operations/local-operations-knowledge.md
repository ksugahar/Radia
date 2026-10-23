# Local operations knowledge

`radia_mcp_operations()` lists locally registered operations topics.
`radia_mcp_operations(topic="<id>")` reads an exact topic with its source
provenance. The meta server catalog advertises this read-only tool. It does not
connect to hosts, execute scripts, submit/cancel jobs or grant resource admission.

Site manuals, connection details and derived site instructions stay outside
public source and wheels. Register UTF-8 JSON at the per-user
`radia_mcp/operations/knowledge.json` below `LOCALAPPDATA` (or `~/.cache`), or
set `RADIA_OPERATIONS_KNOWLEDGE` to a private file. The tool rereads the file on
each call so data updates need no server restart. New tool registration requires
normal MCP startup/reconnection; a package source edit alone is not live proof.

Keep credentials out of this private JSON too: never store private keys,
passwords, tokens or recovery codes. Public-key fingerprints may be recorded
when needed; key material is not needed for operations knowledge.

The document schema is `radia.operations-knowledge.v1` with a `topics` list.
Each topic has nonempty `id`, `title`, `summary`, `content`, and a nonempty
`sources` list. Sources should identify the local manual path, file SHA-256,
PDF page numbers and printed page/section where applicable. Live observations
should retain their command, capture time and hashed output separately.
Missing registration returns `not_configured`; malformed or duplicate topic
data fails loudly. Unknown ids list the available topics and never become
filesystem paths.

Distinguish vendor PBS syntax, site-specific wrapper/resource restrictions,
historical trial notes and current read-only scheduler observations. Queue
visibility is not entitlement or permission to submit. Queue limits can change:
refresh them before an admitted job. Preserve Windows-to-Linux shell quoting
and LF job-script line endings; do not mistake a local Windows PBS example for
a Linux cluster submission recipe.
