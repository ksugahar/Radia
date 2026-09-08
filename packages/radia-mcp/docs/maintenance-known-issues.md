# Maintenance issue ledger

Owner is the named subsystem, not an assumed student or external maintainer.
This is a curated backlog, not a generated tool inventory. Update evidence and
status together; a test file merely existing is not a successful run.

| ID / priority | Owner / concern | Evidence or required check | Status |
| --- | --- | --- | --- |
| M01 / high | deployment: editable and loaded source drift | `maintenance doctor` with approved root/version/full SHA, then real transport status for each client/user | Automated checks implemented; user-specific acceptance required at each update |
| M02 / high | deployment: Windows locked entry points | pip must exit successfully after controlled client shutdown; reconnect and probe | Operational limitation; no automatic process killing |
| M03 / high | paper: absent inputs and detector failures incorrectly passed | `test_paper_writing_review_regressions.py`: supplied missing files and detector-status mapping | Regression available; limited to asserted cases |
| M04 / high | bibliography: nested braces, false title matches, transient DOI errors | Same suite: balanced parser, exact titles, HTTP-status classification | Regression available; not certification of every bibliography parser |
| M05 / high | paper: destructive encoding/newline normalization | Same suite: line-ending preservation and no replacement-character writes | Regression available; no blanket approval of other writing tools |
| M06 / medium | paper: PDF cropbox, float proximity, whitespace false reports | Same suite: cropbox, aux-page references, missing-aux behavior, whitespace threshold | Regression available; representative PDFs still required |
| M07 / medium | paper: remaining approximately 40-item external review | Reconcile the complete 2026-09-03 review against current code, one finding per test/closure | Open: pasted summary is not a complete closure checklist |
| M08 / high | packaging: wheel missing skills/canonical bibliography | Existing wheel-content verifier and installed-wheel lane; editable-only pass is insufficient | Required release gate; not replaced by doctor |
| M09 / medium | publication: JA/EN/grant have distinct review objectives | `test_capability_packs.py` checks distinct routes and unsupported English grant scoring | Covered route contract; no unified quality score claimed |
| M10 / medium | CAD: actual meshing vs available MCP tools | Cubit conformal gate regression, plus licensed live CAD/mesh lane on affected changes | Contract regression available; connectivity is not numerical acceptance |

Completion requires matching the approved revision, zero unexpected failures
or skips in the relevant lane, and retained local machine-readable evidence.
Do not commit transient test logs or mark the package globally "complete" from
one maintenance pass.
