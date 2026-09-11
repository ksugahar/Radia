# Paper-writing review reconciliation — 2026-09-11

## Baseline and interpretation

Source inspected: `f5e083375f82233b69fe64097ca570bc671db45e`.
Original: `C:/temp/pw_review/REVIEW_paper_writing_2026-09-03.md`, read in full.
Its baseline was `ddd50dfe4`, compared with `d88fbf9d2`; its line numbers,
installed versions and tool counts are historical, not current assertions.
This document maps that report rather than copying its old recommendations.

Statuses:

- **修正済み / verified fix**: the named reproduction is covered by a test that
  passed in this audit. This is not certification of the entire tool.
- **未修正 / reproduced open**: current execution reproduces the finding, or
  the precise defective operation is still present (explicitly marked static).
- **検証待ち / pending**: no current reproduction/closure proof was completed.
  This does not assert that the old defect still exists.
- **継続検証 / ongoing acceptance**: environment, artifact or release checks,
  not an open source-code bug.
- **方針変更 / superseded**: the old proposed action conflicts with current policy.

Related occurrences are grouped below with original section references. Row
counts must not be advertised as a count of independent defects.

## Executed evidence

From the full repository checkout, process-local `PYTHONPATH` selected
`packages/radia-mcp/src`; no editable installation or client was changed.

```powershell
python -X utf8 -m pytest -q -c packages/radia-mcp/pyproject.toml packages/radia-mcp/tests/test_paper_writing_review_regressions.py packages/radia-mcp/tests/test_paper_writing.py packages/radia-mcp/tests/test_paper_writing_conclusion_first_use.py --basetemp C:/temp/paper-ledger-a
```

Result: **170 passed, no failures or skips**, 6.98 seconds. This is a focused
three-file lane, NOT the complete package suite or live-client acceptance.
`R::name` below abbreviates `tests/test_paper_writing_review_regressions.py::name`.
Paths in tables are relative to `src/radia_mcp/paper_writing` unless specified.

## Gate, source, bibliography and registration

| ID | Original | Status | Evidence / next closure check |
| --- | --- | --- | --- |
| P01 | 1.1; 1.2 | 継続検証 | Old checkout drift/divergence does not define today's intended source. Check interpreter and original clients when updating. |
| P02 | 1.1 uninstall-first; 9.1 | 方針変更 | No uninstall-first or automatic canonical-tree restoration. MCP development permits coordinated editable repointing. |
| P03 | 1.3; 3.7; 6.1; 6.12; 7.3 filename failures | 検証待ち | Audit all default entry points and instructions against the single parent `references.bib` and generated `.bbl` policy; do not restore local bibliography copies. |
| P04 | 1.4; 6 canonical-path paragraph | 継続検証 | Wheel-content and installed-wheel checks are release gates. No wheel was built/installed in this audit. |
| P05 | 1.5; 3.8 | 修正済み | R::test_computational_electromagnetics_venues_are_supported (IGTE/COMPUMAG/CEFC). IEEE PAC alias precedence remains pending separately. |
| P06 | 2 missing five tools; 7.1 counts | 検証待ち | Check current `tools/list` against the five named functions and current advice; old fixed counts are not an oracle. |
| P07 | 2 duplicate acronym and shared packages | 検証待ち | Trace current consumers before consolidating; prove both public routes retain behavior. |
| P08 | 2 stale docstrings/selftest thresholds | 検証待ち | Reconcile descriptions and registration assertions without regenerating a tool inventory. |
| P09 | 3 missing inputs | 修正済み | R::test_submission_gate_rejects_supplied_missing_files. |
| P10 | 3.1; 3.2 | 修正済み | R::test_submission_gate_maps_detector_failures_to_status checks figure references and underlines. |
| P11 | 3.4 equation/reference/citation/abstract length | 修正済み | Same test checks each result status, not just check-name presence. |
| P12 | 3.3; remaining 3.4 self-citation/background | 修正済み | test_paper_writing_fault_injection.py covers all 24 detector adapters with exception/error/status=error/ok=false/skip/not-applicable/empty/None/malformed payloads. Includes self-citation OK/WARN/HIGH and background 0.25/0.40 boundaries, nonfinite/invalid measurements, preprocessing failures and an all-clean control. Missing evidence fails; skipped checks cannot claim submission-ready. This certifies adapter failure propagation, not each detector's linguistic or PDF accuracy. |
| P13 | 3.5 | 検証待ち | Gate consumers must see nested citations consistently; resolver-only success does not close this. |
| P14 | 3.6 | 検証待ち | Observe merged-source temporary-file cleanup on success and exception paths. |

## Text and citation diagnostics

| ID | Original | Status | Evidence / next closure check |
| --- | --- | --- | --- |
| P15 | 4.1 | 修正済み | R::test_display_math_scanner_does_not_treat_line_spacing_as_math protects plain-text preservation. Other three consumers still need parameterized coverage. |
| P16 | 4.2 IMRaD/conclusion aliases | 検証待ち | Static inspection: plural Conclusions and Results and Discussion patterns now exist. Exercise section boundaries, combined headings and bibliography exclusion before closure. |
| P17 | 4.2 References page boundary | 検証待ち | Reproduce a body page mentioning “Cross References” before the bibliography. |
| P18 | 4.3 Japanese abstract splitting | 修正済み | R::test_japanese_abstract_sentences_split_without_spaces. |
| P19 | 4.3 ending variety/decimal/shadowed endings | 検証待ち | Decimal-preserving sentence split and precedence tests required. |
| P20 | 4.3 misleading ratio/Fig. 3 | 検証待ち | Keep both source values with the ratio across abbreviations. |
| P21 | 4.3 English stop words/passive denominator | 検証待ち | Separate English token repetition and passive-voice ratio fixtures. |
| P22 | 4.3 percent in conclusion | 修正済み | R::test_conclusion_plain_text_preserves_percent_suffix. |
| P23 | 4.4 citation usage @string/@comment/natbib | 検証待ち | Balanced parser test alone does not cover all citation-usage consumers or two optional arguments. |
| P24 | 4.4 reference format @/bare year | 検証待ち | Exercise email/URL values and unbraced numeric year in the public lint. |
| P25 | 4.4 self-citation Han/Hane | 検証待ち | Exact author-token behavior needs reproduction. |
| P26 | 4.4 local bib search | 検証待ち | Verify parent bibliography defaults and explicit bibliography resolution without local copies. |
| P27 | 4.4 cref lists | 修正済み | R::test_cref_comma_list_is_split_into_individual_keys covers figure references; equation cref/vref remains pending. |
| P28 | 4.4; 6.3; 6.4 | 修正済み | R::test_balanced_bib_parser_and_title_match_are_exact covers nested title, @string exclusion and exact title matching. Other parser implementations/arXiv reuse remain pending. |
| P29 | 4.5 three-letter acronym exemption | 検証待ち | Reproduce POD without definition; do not infer closure from a different acronym route. |
| P30 | 4.5 lowercase expansion/CJK boundary/whitelist | 検証待ち | Test LPOD expansion, FEMを, IEEE/SI/Section II in both routes. |
| P31 | 4.6 reviewer regex | 修正済み | R::test_reviewer_figure_caption_pattern_and_priority covers Fig./Figure, case and spacing, negative word boundaries and D/C/B precedence over minor-caption detection. |
| P32 | 4.6 D response schema | 修正済み | R::test_reviewer_d_returns_matched_triggers verifies D results include the matched triggers. This does not certify semantic classification accuracy. |
| P33 | 4.6 tense section argument | 検証待ち | Compare identical prose under Methods/Results/Introduction and documented contract. |
| P34 | 4.6 missing-conclusion response schema | 検証待ち | Assert all public keys on not-found and update consumers consistently. |
| P35 | 4.6 concept-drop exact text | 検証待ち | Test emphasized DtN text against original-source spans. |
| P36 | 4.6 typography scope/dead warning branch | 検証待ち | Exercise thanks/TikZ size commands and warning counts, not just ordinary body sizes. |
| P37 | 4.6; 8.6 personal author default | 修正済み | R::test_cover_letter_does_not_invent_author covers empty, whitespace, None and explicit names; omitted author becomes a neutral placeholder. |
| P38 | 4.6 overfull log forms | 修正済み | R::test_overfull_hbox_context_variants covers paragraph/alignment/detected/output-active, wrapped locations and absent context. R::test_overfull_hbox_count_not_limited_to_detail_cap preserves total count and excludes underfull/vbox. Unknown source locations remain empty, not invented. |
| P39 | 4.6 PDF edge overflow/docstring/headers | 検証待ち | Test drawings/images and header/footer cases separately from page-boundary overflow. |

## PDF and composite plans

| ID | Original | Status | Evidence / next closure check |
| --- | --- | --- | --- |
| P40 | 5.1 floats always clean | 修正済み | R::test_float_distance_without_aux_is_not_reported_as_pass and test_float_distance_uses_aux_page_and_pdf_prose_reference. |
| P41 | 5.1 whitespace threshold | 修正済み | R::test_whitespace_detector_default_does_not_use_old_75_percent_cutoff. Real-document sensitivity remains 継続検証, not certified by a signature check. |
| P42 | 5.1 invalid range/dependency errors | 検証待ち | Check malformed page_range and missing PyMuPDF contracts. |
| P43 | 5.2 crop origin | 修正済み | R::test_cropbox_with_nonzero_pdf_origin_does_not_create_false_overflow. Actual right-edge clipping needs a separate positive fixture. |
| P44 | 5.2 math text-block overlap | 検証待ち | Representative display-math PDFs needed. |
| P45 | 5.2; 5.3 vector-only figures | 検証待ち | Check vector drawings are not blank pages and participate in relevant geometry checks. |
| P46 | 5.3 equation number/body reference | 検証待ち | Distinguish standalone number from prose ending “as in (3)”. |
| P47 | 5.3 caption geometry | 検証待ち | Check image top/bottom, horizontal overlap and prose “Figure 1 shows”. |
| P48 | 5.3 duplicate image xrefs | 検証待ち | Repeated placements must not self-overlap. |
| P49 | 5.3 INFO score/docstrings | 検証待ち | Base-14 font information alone must not imply a defect. |
| P50 | 5.3 extraction errors/resource closure; 8.2 | 検証待ち | Four PDF detectors have failure/closure evidence. PNG renderer and thumbnail helpers now close PDF, tiles and canvas on parsing/render/decode/composition/save/stat failures; test_page_png_renderer_closes_on_every_exit and test_thumbnail_strip_closes_pdf_tiles_and_canvas include successful and empty-selection controls. Other PDF consumers and visual accuracy remain unverified. |
| P51 | 5.4 T18 phase4 | 修正済み | test_paper_writing_fault_injection.py::test_workflow_keeps_active_triggers_and_forwards_bib verifies a nonempty phase4, bibliography and author forwarding. Incomplete/invalid reviewer results no longer produce a completion summary. This does not certify other phase payloads. |
| P52 | 5.4 T9 detector exception | 修正済み | test_paper_writing_fault_injection.py injects exceptions, error/ok=false/skip, empty and non-dict results into each of seven detectors, plus complete failure and invalid counts. Unknown checks suppress score/risk (None) while retaining confirmed findings; a complete clean fixture retains score 10. |
| P53 | 5.4 T10 English triangle | 修正済み | R::test_triangle_uses_lowercase_english_content_words. |
| P54 | 5.4 T13 bib entry | 修正済み | R::test_citation_health_parses_single_line_nested_entry. |
| P55 | 5.4 T11 coil/Table metadata | 修正済み | R::test_reproducibility_does_not_read_coil_or_table_as_metadata. |
| P56 | 5.4 T12 power/p-value | 修正済み | R::test_physical_unit_is_not_misread_as_p_value. |
| P57 | 5.4 T3 unit boundaries/T4 absent discussion | 検証待ち | Reproduce “2 settings”, “Sec. 3.” and missing-section behavior. |
| P58 | 5.4 T12/T13/T15 unreachable rules | 検証待ち | T16 now lists unknown_patterns and missing_signals instead of inferring a full pattern from partial evidence. test_root_cause_requires_all_pattern_signals covers this. Missing T8 producers (including T9/T12/T13/T15) are not added by this fix and remain open. |
| P59 | 5.4 sort/zero/substring/year/schema nits | 検証待ち | T8/T20 propagate partial/unavailable into T18 phase2. T16/T17 now propagate invalid/incomplete health and T18 phase3/5/6 reject invalid results. test_actions_separate_unknown_from_real_zero closes the T17 zero-versus-None case. T8/T20 ordering, T20 substring, T5 year and other schema cases remain open. |

## External data, documentation and cross-cutting work

| ID | Original | Status | Evidence / next closure check |
| --- | --- | --- | --- |
| P60 | 6.2 transient DOI failure | 修正済み | R::test_temporary_doi_failure_is_not_called_fabrication and test_crossref_http_status_classification. Retry-After/UA policy remains pending. |
| P61 | 6.5 nested input resolution | 修正済み | R::test_input_chain_resolves_from_main_compile_directory. Explicit .pgf suffix and unreadable includes still need checks. |
| P62 | 6.6 arXiv escaping/error feed | 修正済み | test_paper_writing_external_inputs.py verifies Requests-prepared query roundtrips for &, #, C++, percent and Japanese text. Atom error entries, unexpected XML roots, missing title/invalid identifier, HTTP errors and genuine empty feeds are distinguished. Citation verification preserves search failure as verdict=error rather than no_candidate_found. All network calls are mocked. |
| P63 | 6.7 metadata/BibTeX fidelity | 検証待ち | Missing published/issued, organization authors, escaped entities and complete arXiv authors. |
| P64 | 6.8 damaged PDF persistence | 修正済み | test_pdf_magic_alone_does_not_validate_corruption and test_pdf_verification_requires_parser reject magic-only validation. _save_verified_pdf streams into a same-directory temporary file, validates, then replaces the destination; failure removes only its temporary file. Tests cover interruption/size limits and all three publisher adapters preserving an existing file before a valid retry. No live publisher access or visual-fidelity certification. |
| P65 | 6.9 DOI normalization | 検証待ち | Test #/?/%, doi: and dx.doi.org consistently across consumers. |
| P66 | 6.10 archive resource bounds | 修正済み | test_arxiv_source_resource_limits covers streamed download, gzip expansion, member count, tar-member/single-source size limits, interruption and successful tar/gzip controls. Limits apply before unbounded download/decompression; tar members are iterated with a cap rather than getmembers(). No archive files are extracted to disk. |
| P67 | 6.11 fallback paths | 検証待ち | arXiv search-to-citation fallback now rejects error/malformed/inconsistent payloads and preserves verdict=error (external-input tests). Unreadable bib, resolver exceptions and publisher landing errors still need their own closure evidence. |
| P68 | 6 duplicate implementations/docstring promise | 検証待ち | Trace Crossref/arXiv/BibTeX/DOI consumers before deduplication; verify promised search routes. |
| P69 | 7.2 missing/stale tool guidance | 検証待ち | Correct names, routes, supported callable advice and guidance organization; old counts are historical. |
| P70 | 7.1; 7.2 generated inventory proposal | 方針変更 | TOOLS.md inventory gate is retired. Live discovery/catalog is the source of truth, not a new committed generated index. |
| P71 | 7.3 collection and coverage counts | 継続検証 | This checkout's focused lane collected and passed 170 tests. Shared checkout/root collection, full suite and live registration were not tested. |
| P72 | 7.3 related_servers/catalog | 検証待ち | Reconcile current aliases and routes against catalog. |
| P73 | 8.1 language policy | 検証待ち | Decide current writing-tool language policy before bulk translation; old English-only claim is not authorization. |
| P74 | 8.3; 8.4 duplication/source ownership | 検証待ち | Map current shared lint consumers and presentation copies; do not mass-delete apparent duplicates. |
| P75 | 8.5 destructive normalization | 修正済み | R::test_terminology_normalizer_preserves_line_endings and test_terminology_normalizer_never_writes_replacement_char. Backup policy is a separate decision. |
| P76 | 8.6 other personal-name heuristics | 検証待ち | P37 proves only the cover-letter default, not every public heuristic. |

## Next bounded work

Follow-up fault-injection verification (2026-09-11): 481 passed across
`test_paper_writing_fault_injection.py`, `test_paper_writing_review_regressions.py`,
`test_paper_writing.py`, and `test_paper_writing_conclusion_first_use.py`.
This is a focused source-isolated test run, not a full-package or live-MCP check.
No editable installation or live client was changed.

The subsequent P31/P32/P37/P38 fix adds 21 regression cases in the same four-file
lane. These are source-level checks, not a live-client or full-package acceptance.

The T8 -> T20 -> T18 phase2 follow-up preserves unknown/skipped checks, suppresses
incomplete overall scores, retains confirmed findings, and distinguishes detector
failure from a confirmed critical manuscript defect. Context filtering does not
erase UNKNOWN issues. The four-file focused lane now passes 593 tests; no live
client or editable source was changed.

The subsequent T16/T17 and T18 phase3/5/6 follow-up adds structured unavailable
results, preserves unresolved checks separately from rewrite actions, requires all
root-pattern signals, and rejects failed rewrite/response generation. Focused
fault-injection, real-health propagation and successful-generation controls pass
in the same four-file lane (626 tests). This is not a semantic accuracy audit;
missing pattern producers remain tracked in P58. No live client was changed.

The PDF-analysis follow-up verifies four detectors with mocked page extraction
and raster failures, successful positive controls, document closure, and a real
detector-to-gate error path. The same four-file lane passes 642 tests. No PDF
artifact was authored or visually reviewed, and no live MCP client was changed.

The PNG/thumbnail resource-lifecycle follow-up adds 14 focused cases (656 passed
in the same four-file lane). Mocked rendering verifies exception propagation,
resource closure and unchanged success metadata, not actual image fidelity.
Partial output files may remain after failure; cleanup does not delete outputs.
No live MCP client or editable installation was changed.

The arXiv search/citation follow-up adds offline external-response tests. Search
parameter encoding was already correct; error-feed/root validation and fallback
error propagation were strengthened. The five-file focused lane (the prior four
files plus test_paper_writing_external_inputs.py) passes 680 tests. No live API,
download, bibliography modification or MCP reconnection occurred.

The download/archive follow-up passes 697 tests in the same five-file lane.
PDF publication is staged and verified (PyMuPDF required); the source archive
lane caps streamed bytes, expanded bytes, member count and per-file reads. Tests
use synthetic in-memory archives and mocked HTTP, not live publisher downloads.
No live MCP client, editable installation or parent bibliography was changed.

1. Continue metadata/DOI fidelity and remaining fallback paths (P63/P65/P67), then language heuristics and
   deduplication; other PDF consumers remain separate P50 follow-up work.
2. Review missing pattern producers (P58) and remaining ordering/context rules
   (P59); do not confuse execution-status coverage with semantic correctness.

Update a row's evidence and status in the same change. Do not use passing tests
for a neighboring consumer, a changed regex, or a previous release's live status
as closure evidence. This audit intentionally makes no package-wide readiness claim.
