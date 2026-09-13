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

## Notebook references

For `.ipynb`, the same `bibliography_make_bbl` entry point reads only
`metadata.radia.bibliography.keys` (a nonempty ordered list of unique canonical
keys) and optional `style`. It never executes cells or guesses identities from
author/year prose. No keys means no bibliography request: do not add dummy entries
to citation-free notebooks. A notebook cannot use `aux_path`.

```python
bibliography_make_bbl(tex_path="docs/open_boundary/open_boundary_demo.ipynb")
```

This writes the sibling `open_boundary_demo.bbl`, not a local `.bib`. Keep that
generated file with the notebook. For readable display, use the installed TeX4ht
renderer, not a hand-maintained second bibliography. In a scratch directory,
copy the generated bbl as `references.bbl` and create `references.tex`:

```latex
\documentclass{article}
\usepackage[T1]{fontenc}
\usepackage[utf8]{inputenc}
\usepackage{url}
\begin{document}
\input{references.bbl}
\end{document}
```

Run `make4ht -u references.tex`. Preserve the generated HTML body in the notebook's
single reference Markdown cell, tagged with `metadata.radia_bibliography` holding
the sibling `bbl` filename, its `bbl_sha256_lf` (SHA-256 after CRLF-to-LF
normalization for portable Git checkouts), and the `renderer` version. Retain
the bbl link and all generated entry anchors. Keep scratch TeX/HTML/CSS out of
docs; existing calculation cells and outputs are not rerun by this operation.
If BibTeX or TeX4ht is unavailable or rendering fails, report the missing stage;
do not hand-author replacement reference text. This adds no Python dependency.

After a parent entry, citation key or style changes, regenerate both bbl and
display. The fast docs contract checks declared keys, generated item membership and
display/bbl hash association; it does not prove scholarly citation completeness
or independently re-run TeX4ht. Preserve an explicit unresolved list for legacy
notebooks rather than calling absent declarations compliant.

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
