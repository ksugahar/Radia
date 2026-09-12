# Bibliography finalization contract

## Citation decisions

Use `bibliography_make_bbl(tex_path, style="", out_path=None)` for ordinary
static citation/input syntax. For macros, conditional text or multicite commands,
compile the manuscript using its actual TeX engine and pass that run's main aux:

```python
bibliography_make_bbl(tex_path="paper.tex", aux_path="paper.aux")
```

The aux reader consumes only citation, bibstyle and nested @input records.
It never executes aux commands or launches a TeX engine. BibTeX is still run in
an isolated temporary directory against a snapshot of the single canonical
references.bib. Explicit `style` overrides the compiled style; conflicting styles
inside the aux chain, missing inputs, cycles and malformed records fail closed.
This is the BibTeX aux contract, not a biber/bcf implementation.

The caller must regenerate aux after changing the source or compile options.
Hash checks detect edits during generation, not a stale aux supplied before the
operation. The result explicitly labels compiled-aux freshness as caller-owned.
This is an artifact handoff to the actual TeX interpreter, not a second partial
macro interpreter with a claim of universal TeX support.

## Writer exclusion

Source edits and bbl generation acquire a fail-fast operating-system advisory
lock for their destination. Competing cooperating processes cannot publish to
that target simultaneously. After acquiring the lock, snapshot comparisons still
prevent overwriting a source changed since the operation read it.

The stable sibling `.<filename>.radia-write-lock` is deliberately retained.
The OS releases its lock when the process exits, including abnormal termination;
an existing sidecar does not indicate a busy writer. Do not delete the sidecar
while writers may be active: replacing its inode breaks mutual exclusion.
Unrelated editors that ignore this protocol are not completely excluded, and
the final compare/replace window remains optimistic for them. Hard-link aliases
and filesystem-specific network locking semantics need separate acceptance.

## Live deployment

Source integration, interpreter registration and live-client verification are
different acceptance steps. A fresh Python selftest is not live acceptance.
An added optional tool parameter changes the advertised schema and requires a
client-owned reconnect, not only code reload.

At the finalization check on 2026-09-12, this task's live radia-publication
bibliography status reported version 1.4.53 from the shared
`release-quad/radia-4.95.91-runtime` tree, not the maintenance worktree.
No idle-state evidence or callable client reconnect control was available.
The runtime was not repointed or reloaded. Integrate reviewed changes into the
chosen shared source, establish a quiet boundary, restart only radia-publication
through the client, then verify its schema/provenance and a harmless affected
tool through that same connection. This remains an operational handoff, not a
code defect to mask with process termination.
