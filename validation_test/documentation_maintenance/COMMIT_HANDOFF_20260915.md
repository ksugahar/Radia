# Documentation maintenance commit handoff

This commit records the documentation work as a coherent checkpoint, not a
claim of complete rendered/publication acceptance or a numerical solver release.

Included:

- Method/equation/implementation/limitation narratives for 53 public notebooks.
- Canonical citation declarations, generated BibTeX reference displays and BBLs.
- Verified bibliography corrections/additions and reused bibliography tooling.
- Eight internal notebook ledgers relocated with their companion records.
- Twelve added saved WebGUI scenes across four notebooks, with distinct
  historical-campaign and new mdx2 execution provenance.
- Corresponding MCP discovery/caveat updates and focused regression tests.

Pre-commit verification: 111 tests passed across canonical/finalized
bibliography, document-meta auditing, generated references, and both saved
WebGUI suites. `git diff --cached --check` passed. The latest citation audit
resolved all 53 notebooks against the 520-entry worktree parent without gaps.
No solver kernels, MATLAB implementation changes, release versions, or unrelated
README/package changes are included in this docs checkpoint.

Remaining acceptance work:

- Rendered visual QA in a permitted WebGUI-capable viewer; saved scene/state
  checks alone do not prove useful camera framing, vector density or clipping.
- Broader field/mesh-view applicability review and additional examples.
- Remaining bibliography identity reconciliation and live MCP source alignment;
  committing this parent does not repoint any running MCP process.
- A genuine AXIFEM public example and other previously recorded evidence gaps.

The dated WEBGUI and CITATION_INTEGRATION notes describe intermediate states;
their statements that docs were uncommitted are historical, superseded by this
checkpoint. Numerical evidence was already committed separately as `68f685fe3`
and `c39d56490`. Both mdx2 task directories were removed after hash-verified
recovery and evidence commits. LAB scene recovery is retained for pending
visual QA under the ownership/cleanup triggers in those notes.
