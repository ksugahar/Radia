# Maintenance issue ledger

Owner is the named subsystem, not an assumed student or external maintainer.
This is a curated backlog, not a generated tool inventory. Update evidence and
status together; a test file merely existing is not a successful run.

The [2026-09-11 paper review reconciliation](design/paper_review_reconciliation_20260911.md)
maps the complete original report to numbered finding groups and current evidence.
Use its statuses: 修正済み (named reproduction verified), 未修正 (reproduced or
static-confirmed), 検証待ち (not established either way), 継続検証 (operational
acceptance), and 方針変更 (old recommendation superseded). Release gates and
client verification are not automatically open code defects.

| ID / priority | Owner / concern | Evidence or required check | Status |
| --- | --- | --- | --- |
| M01 / high | development: editable and loaded source mismatch | `maintenance doctor` with selected root/version/full SHA; separately inspect the affected original clients | Intentional repointing is allowed; fresh-process evidence is not live-client verification |
| M02 / high | deployment: Windows locked entry points | pip must exit successfully after controlled client shutdown; reconnect and probe | Operational limitation; no automatic process killing |
| M03 / high | paper: absent inputs and detector failures incorrectly passed | Reconciliation P09–P14 | Input/status cases and all 24 gate detector adapters have focused fault-injection coverage (P12); lower-detector accuracy and nested consumers remain separate 検証待ち work |
| M04 / high | bibliography: nested braces, false title matches, transient DOI errors | Reconciliation P28, P60–P68 | Named parser/title/HTTP cases 修正済み; metadata/author/escaping, failure paths and HTTP resource closure have scoped evidence. All three S2 routes now share identifier normalization/path encoding with offline regressions (905-test lane). Other producers and rendered bibliography acceptance remain 検証待ち; unsupported math/markup requires manual verification |
| M05 / high | paper: destructive encoding/newline normalization | Reconciliation P75; both regressions passed in the 170-test focused lane | 修正済み for those two behaviors; no blanket approval of other writing tools |
| M06 / medium | paper: PDF diagnostics | Reconciliation P39–P50 | Named crop/float/threshold cases, four detector failure paths and PNG/thumbnail resource closure have focused evidence; visual accuracy, vector figures and other PDF consumers remain 検証待ち |
| M07 / high | paper: external review follow-up | Full report mapped in reconciliation P01–P76; groups may contain related findings | Inventory complete, not fixes complete. Gate/composite execution-status fixes have scoped evidence. Publisher landing/PDF failure stages and resource closure now have 27 offline cases in the 795-test lane. Missing pattern producers, semantic rules, other external-input/PDF consumers and live publisher acceptance remain open |
| M08 / high | packaging: wheel missing skills/canonical bibliography | Existing wheel-content verifier and installed-wheel lane; editable-only pass is insufficient | Required release gate; not replaced by doctor |
| M09 / medium | publication: JA/EN/grant have distinct review objectives | `test_capability_packs.py` checks distinct routes and unsupported English grant scoring | Covered route contract; no unified quality score claimed |
| M10 / medium | CAD: actual meshing vs available MCP tools | Cubit conformal gate regression, plus licensed live CAD/mesh lane on affected changes | Contract regression available; connectivity is not numerical acceptance |

M04 follow-up (2026-09-12): bibliography T2 arXiv identity, explicit version,
required metadata and primary-category preservation have 23 offline regression
cases. A further 12 cases cover T2 shared text escaping and unstructured author
order/separator preservation; the six-file lane passes 940 tests. Unsupported
TeX/math requires manual verification. Citation-key surname heuristics and
rendered bibliography acceptance remain open (P63/P67).

M04 T1 follow-up (2026-09-12): standalone DOI-to-BibTeX now has DOI identity,
required metadata, corporate-author/text escaping and page-range regressions.
The seven-file lane passes 960 tests. T3 Crossref search response handling and
rendered bibliography acceptance remain open; this is not package-wide closure.

Completion of a scoped validation requires matching the selected source, zero unexpected failures
or skips in the relevant lane, and retained local machine-readable evidence.
Do not commit transient test logs or mark the package globally "complete" from
one maintenance pass.
