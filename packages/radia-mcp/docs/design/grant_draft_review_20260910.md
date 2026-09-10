# Draft-check review follow-up — 2026-09-10

Based on main `498409b5a`, including presentation PR #134. This change ports
the separately supplied A/B/C draft work, not the shared mixed-owner working
tree. Original patches and the release editable installation remain untouched.

## Disposition

- Added the three draft tools in a private module and retained public
  `grant_writing_*` registration through `tools.py` and paper-writing composition.
- Question punctuation and heading boundaries have regression coverage.
  Separate announced questions cannot borrow each other's following questions.
  Repeated enumerations remain observations, not a prohibition on multiple aims.
- No draft tool assigns scientific-quality points. Form vocabulary is only a
  candidate index; unmatched vocabulary does not prove missing prose. Length
  capacity is an estimate with explicit status and validated numeric inputs.
- Health-report question candidates remain visible without contributing scores.
- Ported the kanji-cause observations, but avoided claiming the cause is proven
  by a regex and avoided recommending ambiguous replacements of defined terms.
- Latest-main Optuna quality tests pass. The old 47/45/95 expectation failure
  does not reproduce here; no Optuna source or oracle policy was changed.
- Latest-main MEX source contains 364 commands, not the old review's 367.
  Corrected the two stale 361-command documentation statements and replaced
  the fixed test count with comparison against the general C++ command table,
  including uniqueness. Existing named-capability checks remain in place.
- Generated TOOLS.md is not an inventory gate on current main; no snapshot added.
- Shared junction/lock cleanup and release deployment are outside this isolated
  integration. No active runtime was restarted or uninstalled.

The follow-up's 34 regression cases cover punctuation, short and wrapped questions,
applicability, question preservation in health reports, limits and invalid
inputs, lexical false certainty, file input, and actual MCP registration.
Impacted grant/paper/MATLAB tests pass: 433 tests. Run real paper-writing stdio
and scoped CI before merging. This is not a full-package audit or a live-client refresh.
