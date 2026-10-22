# Historical validation records

This directory holds one-time evidence, not a numerical release lane. Package-owned records belong in `packages/<package>/validation_records/`.

Each record has `record.json`, a hash-checked result JSON and a local entry script with SHA256. Use `tools/validation_record.schema.json`. Preserve original source, neutral runtime role, units, limits and timings. Keep failures and never rewrite raw results or claim stale evidence as current.

The stdlib checker supports only the fixed keywords used by this repository-owned schema and rejects unsupported keywords; it is not a general JSON Schema validator.

Run `python tools/check_validation_records.py --report <owned-output>/records.json` for read-only schema, path, hash, covered-content staleness checks. `repo-contracts` runs the same gate without executing recorded physics. Privacy patterns are a guard, not publication clearance; source/data still require confidentiality review. Creating this directory authorizes no record migration.

Each record is a development milestone. Required metadata includes `milestone`, a one-line `claim`, calendar `date` (YYYY-MM-DD), `radia_version`, the exact `source.repo_sha`, and hash-checked `result` and entry artifacts. `conditions` declares the problem, mesh and neutral hardware description; `environment.threads` records the thread count. `comparison` names the reference and its independence (including any limitation, such as a same-solver regression lock). This metadata does not establish oracle independence by itself.

`status` is `current` or `stale`; `covered_paths` declares the source paths whose changes invalidate current status. Preserve the stale reason and original evidence. `timing_s` values are descriptive records only and must never determine acceptance; the checker does not compare recorded times with limits.

Staleness uses `covered_sha256`, a nonempty map from covered repository paths to source-time SHA256. UTF-8 text normalizes CRLF and CR to LF; other bytes are hashed unchanged. Store a hash per matched file, with dependency paths selected by the record author when relevant to the claim; no root or package pyproject is forced. New, removed or changed covered files invalidate `current`; history ancestry is not required. Source SHA is provenance and may be unavailable after squash/rewrite. Optional `environment.ngsolve` records the actual experiment version, including pre-pin trials; it is not forced on unrelated packages. Only explicitly covered dependency metadata changes invalidate current status; the author is responsible for selecting a sufficient covered scope. Bundle the local script with the record so checking it needs no historical Git objects.

`source.source_only=true` explicitly identifies a source-only milestone. Otherwise at least one native or wheel SHA256 is required, in addition to the always-required result and entry artifact hashes. `superseded_by`, when set, must resolve to an existing different record in the same distribution. The checker never executes entry scripts or uses recorded timings to accept/reject a numerical result.
