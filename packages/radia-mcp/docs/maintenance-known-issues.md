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
| M03 / high | paper: absent inputs and detector failures incorrectly passed | Reconciliation P09–P14 | Named input/status cases 修正済み; full detector fault injection and nested gate consumers 検証待ち |
| M04 / high | bibliography: nested braces, false title matches, transient DOI errors | Reconciliation P28, P60–P68 | Named parser/title/HTTP cases 修正済み; other consumers and external-input handling 検証待ち |
| M05 / high | paper: destructive encoding/newline normalization | Reconciliation P75; both regressions passed in the 170-test focused lane | 修正済み for those two behaviors; no blanket approval of other writing tools |
| M06 / medium | paper: PDF diagnostics | Reconciliation P39–P50 | Named crop/float/threshold cases 修正済み; positive clipping, vector figures and extraction failures 検証待ち |
| M07 / high | paper: external review follow-up | Full report mapped in reconciliation P01–P76; groups may contain related findings | Inventory complete, not fixes complete. P31/P32/P37/P38 未修正; prioritize gate/composite failure evidence P12/P51/P52 |
| M08 / high | packaging: wheel missing skills/canonical bibliography | Existing wheel-content verifier and installed-wheel lane; editable-only pass is insufficient | Required release gate; not replaced by doctor |
| M09 / medium | publication: JA/EN/grant have distinct review objectives | `test_capability_packs.py` checks distinct routes and unsupported English grant scoring | Covered route contract; no unified quality score claimed |
| M10 / medium | CAD: actual meshing vs available MCP tools | Cubit conformal gate regression, plus licensed live CAD/mesh lane on affected changes | Contract regression available; connectivity is not numerical acceptance |

Completion of a scoped validation requires matching the selected source, zero unexpected failures
or skips in the relevant lane, and retained local machine-readable evidence.
Do not commit transient test logs or mark the package globally "complete" from
one maintenance pass.
