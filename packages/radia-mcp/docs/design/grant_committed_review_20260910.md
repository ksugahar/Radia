# Committed grant-writing integration review — 2026-09-10

## Scope and ownership

This is a focused integration of committed grant-writing corrections, based on
main `d6db6ad780868b4d3949dd89707e2b83a269ef4c`. It is not a merge of the
divergent `backup/main-pre-release-20260821` branch or a review of all its history.
The shared working tree, its mixed-owner uncommitted changes, and the separate
A/B/C patches remain untouched. Use a fresh branch from the accepted main to
review and apply those patches; do not assume they apply without conflict.

## Disposition

- Integrated `6c89576a5`: ignore researcher role titles as integration-claim
  triggers; support grouped budget headings; correct stale figure-server and
  corpus guidance; remove the obsolete software name from test fixtures.
- Preserved main's wrapped-paragraph regression when resolving the overlapping
  test insertion. No main implementation was replaced with an old file snapshot.
- The stale-feasibility and trailing-budget-row behaviors from `ecf5a9b2e`
  already exist on main; no duplicate implementation was added.
- Poster tools already compose into paper-writing on main with compatibility
  routes. The older `6440c2900` standalone-removal patch was not blindly applied
  over that newer composition contract.
- Presentation commit `48c9a37f6` is handled by the separate presentation task;
  this integration deliberately does not touch its files.
- Reader-reconstruction rules (`9102324ed`) and co-listing rules (`7b135627b`)
  are not included in this focused merge. The old branch includes prerequisite
  readability/paired-object behavior absent from main; automatic three-way
  application pulled that ancestry into conflicts. Review that feature set as
  a coherent follow-up, not as part of the mixed-owner working-tree patches.

## Finding corrected before integration

The grouped-budget helper skipped member categories whenever the ledger already
contained a group heading. Distinct expense rows labelled `物品費=200`, `A=300`,
and `B=50` therefore produced a grouped total of 200 instead of 550. An expected
total of 200 could pass while excluding 350 of accepted expenditures.

The corrected helper sums the group-labelled expenditures and all member-labelled
expenditures exactly once. It still expands only actual ledger totals, never
declared totals; raw grand totals are unchanged. Summary rows are the ledger
reader's responsibility, not something inferred from a category name. A focused
regression verifies both accepting 550 and rejecting 200. Existing code/Japanese
heading, grouped-heading, zero-category, role-title and paragraph tests remain.

## Verification

The original port passed 239 grant-writing tests. The added mixed-category
regression then reproduced the incorrect total and passed after the fix.
Combined grant-writing and paper-writing tests passed: 380 tests. Run the
paper-writing real-stdio probe and scoped CI before merging. No editable
installation, running MCP connection, package version or release is changed.
