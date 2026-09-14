# Existing citation workflow integration

The user approved using this worktree as the eventual parent bibliography on
2026-09-14, while preserving the established citation workflow.

## Reused implementation

The active source is the sparse `mcp-development-main` worktree, at
`df452eaa4` when inspected. Its bibliography work is substantially newer than
the corresponding files in this checkout. The following implementation was
integrated, rather than replaced with another citation parser or renderer:

- source-span-preserving BibTeX parsing;
- cooperative writer exclusion and source snapshot checks;
- static/compiled-aux citation handling and TeX lexical boundaries;
- canonical entry lookup and notebook-aware `bibliography_make_bbl`;
- canonical/compiled-aux regression tests.

This checkout's additional DOI verification, unpublished-reference refresh,
stray-bibliography discovery, key checks, search and landmark tools are retained.
The imported regression expectations were migrated from old citation keys to
the DOI/title-verified identities already used by this parent. The selected
regression run passed **102 tests** before the whole-corpus display check.

The generated-display contract was taken from the existing
`tests/test_docs_notebook_contract.py` on that branch. Only its citation
contract is integrated here; unrelated document removals and numerical changes
are not silently imported. The local renderer is an orchestration helper for
BibTeX and TeX4ht, not a new bibliography formatter. Math packages are loaded
for mathematical article titles such as Nedelec's R-cubed title.
Japanese bibliography notes use TeX4ht's CJK support without changing their
meaning or hand-editing the generated `.bbl`.

All 53 tracked public notebooks now declare `metadata.radia.bibliography`,
retain their existing `citation_keys`, and have a generated sibling `.bbl` plus
one generated reference-display cell. This replaces the earlier hand-formatted
reference lists. The schema checks compare key membership, selected parent
entry fingerprints, style, `.bbl` hash and display hash. Original numerical
code/output cells remain unchanged.

Final focused verification: **105 tests passed**, including the whole tracked
corpus's generated-reference contract. All 53 generated displays contain no
Unicode replacement glyph. This is structural/text integrity, not a claim of
complete visual or scientific review.
The citation audit now recognizes DOI links in generated HTML without treating
quote delimiters or closing tags as part of the DOI; a focused regression
protects both single- and double-quoted links.

## Identity reconciliation

The initial [identity comparison](bibliography_reconciliation_20260914.json)
records both source hashes and candidate mappings. It is a review artifact,
not an approved alias table: generated-key or conflicting-metadata matches
must not be silently substituted in manuscripts.

Twelve missing DOI records were fetched through the **existing live**
`bibliography_doi_to_bibtex` tool; the
[lookup receipt](bibliography_crossref_20260914.json) preserves its responses.
Eleven verified records were added to this parent using their existing MCP
keys and protected surnames. Ortner's online-2022 versus issue-2023 date remains
for explicit reconciliation rather than an automatic year rewrite.

`egger2025tmag` is again an explicitly unverified legacy identity. The verified
2025 five-author article is now the distinct `egger2025forward`, and the
hysteresis notebook cites that key. A successful DOI lookup for a similar title
is not permission to reassign an existing citation key to a different paper.

The live parent also contains material not yet incorporated here, including
the nine IEEJ formula-report records and generic `bib1` through `bib13` entries.
Some other unmatched records are spelling/alias candidates. Do not append all
of them merely to make the entry counts agree; reconcile provenance and actual
source identities. The live file remains intact, so none of these records are
lost while review continues.

## Deployment boundary

No editable installation, client configuration or running MCP process was
changed in this step. The original connection still reads its own 501-entry
parent. Source integration and offline tests are not live deployment.
Complete data/key reconciliation before switching the runtime parent, then
use the scoped reconnect and live-provenance verification procedure. Do not
restart Windows, terminate Codex, or interrupt other MATLAB/CAD jobs for this.

Notebook code and saved numerical outputs are preserved. Bibliography
typesetting on LAB is not a numerical recomputation; solver reruns still belong
on mdx1/mdx2 and remain a separate acceptance task.
