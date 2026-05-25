# Changelog

All notable changes to `radia-mcp`. Format: each release lists **what
shipped** + **why** in compact form. Older releases (≤ 0.4) are
omitted; the 0.5 → 0.6 jump is when the standalone `radia-mcp` wheel
crystallized as its own package.

## 0.80.0 — graph: GitHub-MCP/tueplots/Wong-2011 absorption batch

Released 2026-05-26.

After v0.79.0 user prompted to survey existing GitHub MCP plotting
servers + scientific-figure libraries.  Two parallel general-purpose
agents scanned:

  PYTHON MCP:
    StacklokLabs/plotting-mcp, newsbubbles/matplotlib_mcp,
    xlisp/visualization-mcp-server, LindseyyyLi/MCP-Server,
    arshlibruh/plotly-mcp-cursor, antvis/mcp-server-chart (4.1k★),
    isaacwasserman/mcp-vegalite-server (97★),
    stephaneberle9/mcp-server-vegalite-viewer,
    jjsantos01/jupyter-notebook-mcp (130★), Vizro MCP

  MATLAB / scientific-figure libraries:
    matlab/matlab-mcp-core-server (official, 771★),
    Tsuchijo/matlab-mcp, garrettj403/SciencePlots,
    pnkraemer/tueplots, proplot/Ultraplot, BayesWatch/mpl_sizes,
    masumhabib/PlotPub, altmany/export_fig

Distilled 5 high-value patterns this lab DIDN'T already have and
absorbed them this release.  Bigger-picture finding: no existing MCP
server treats figures as first-class publication artifacts; our
emit_paper_figure(...) gate-stack remains genuinely novel.

What shipped:

1. **NEW: Okabe-Ito CVD-safe color palette as default for every
   profile** (Wong 2011 Nature Methods 8:441 = SciencePlots `bright`).
     OKABE_ITO = ['#000000', '#E69F00', '#56B4E9', '#009E73',
                  '#F0E442', '#0072B2', '#D55E00', '#CC79A7']
   - PaperProfile.color_cycle field (default = OKABE_ITO)
   - paper_figure() sets rcParams['axes.prop_cycle'] from the cycle
   - emit_paper_figure() lints every Line2D color via
     _check_colors_are_cvd_safe (Okabe-Ito + greyscale exception),
     raises on violation.  Matplotlib's default tab10 (red+orange
     confusable in deuteranopia) no longer leaks into lab figures.

2. **NEW: PaperProfile.from_base() tueplots-style derivation** —
   declare ONE `base_pt` font; legend_pt = base_pt, tick_pt =
   base_pt - small_offset auto-computed.  Pins the IEEE/IEEJ
   "tick 1 pt below body" convention so future profile edits
   cannot desync.

3. **NEW: paper_figure(rel_width=1.0)** — fraction of profile width
   (tueplots `rel_width` pattern).  rel_width=0.5 of double-column
   gives a half-column inset; rel_width=1.5 of single-column gives
   a 1.5-column figure.  Custom column widths without forking a
   profile.

4. **NEW: paper_figure(panel_labels='auto')** — default is now
   "auto", which applies (a)(b)(c) iff nrows*ncols > 1.  Mirrors
   ultraplot `abc=True` ergonomics.  panel_labels=True forces even
   for 1x1; False suppresses.

5. **NEW: Post-save PDF font-embedding verifier**
   (_check_pdf_fonts_embedded).  Scans the output PDF binary for
   /Subtype /Type3 (raster glyph, IEEE/Elsevier reject) and reports
   Type-1-only embeds.  Catches the case where rcParams['pdf.fonttype']
   got reset to 3 between paper_figure() and savefig().  Returned
   in emit_paper_figure result as `font_violations: list`.

emit_paper_figure() gate stack now (in order):
  pre-flight 1:  no in-figure titles
  pre-flight 1b: NEW colorblind-safe palette (raise on violation)
  pre-flight 2:  no legend overlapping data lines
  pre-flight 3:  resize to profile width
  measure:       axes_area_fraction >= min_axes_fraction (auto_tighten optional)
  post-save:     NEW PDF Type-42 font embedding (raise on Type-3 found)

`paper_figure_quality_rules` (MCP tool) gains 2 new topics:
  `colorblind_safe` -- Okabe-Ito palette + why tab10 is harmful
  `font_embedding`  -- Type-42 requirement + pdffonts verification

tests/test_graph_paper_figure.py: 37 -> 53 tests
  + test_default_color_cycle_is_okabe_ito
  + test_paper_figure_sets_okabe_ito_rcparams
  + test_emit_raises_on_non_cvd_safe_color
  + test_emit_greyscale_passes
  + test_emit_cvd_check_can_be_disabled
  + test_from_base_derives_font_sizes_correctly
  + test_from_base_with_different_base
  + test_rel_width_scales_figure
  + test_rel_width_1_5_for_1_5_column_figure
  + test_panel_labels_auto_applied_on_multipanel
  + test_panel_labels_not_applied_on_single_panel
  + test_panel_labels_can_be_forced_on_single_panel
  + test_panel_labels_false_suppresses_on_multipanel
  + test_emit_pdf_has_no_type3_fonts
  + test_emit_pdf_check_can_be_disabled
  + test_emit_raises_on_type3_fonts

Suite: **124/124 pytest pass** (was 108 after 0.79.0; +16 new).

Patterns surveyed but DEFERRED (not absorbed, with reasons):
  - MCP `Image` content type return (LindseyyyLi, isaacwasserman):
    radia_mcp.graph returns RECIPES (text), not rendered images --
    keeps the server purely informational and lets the user execute
    locally.  Image-return is wrong for our model.
  - `usetex=True` switch (tueplots): adds LaTeX dependency to CI,
    cm/stix mathtext is good enough for IEEE compliance already.
  - Pydantic-typed request models: overkill for our flat kwargs.
  - Out-of-process rendering (antvis): added complexity not justified
    for a local-dev MCP server.
  - generate_sample_data tool (arshlibruh): nice-to-have but
    overlaps with mcp-server-mathematica + Python REPL.
  - DISABLED_TOOLS filter (antvis): catalog is small (5 tools), no
    pruning needed.
  - export_fig CMYK conversion (MATLAB): RGB-output is the dominant
    publisher requirement now (IEEE since 2020); CMYK is print-only.

These can be revisited in v0.81+ if user demand surfaces.

## 0.79.0 — graph: 10pt-at-8cm absolute rule + title/legend gates

Released 2026-05-26.

User corrections to v0.78.0 design:

  1. **Font is ABSOLUTE 10 pt @ 8 cm**, not 8 pt and not relative to
     column width.  Wider columns keep the same 10 pt -- the axes box
     grows, the text doesn't.
  2. **"余白" (waste) = white space between AXES OUTER EDGE and
     FIGURE BOUNDING BOX**, not inside-axes whitespace.  Principle:
     "情報がなく無駄はやめる" -- every mm of figure bbox should be
     axes interior, axis label, tick label, tick mark, or legend.
  3. **Titles go in the LaTeX `\caption{}`, NEVER in the figure.**
     gate raises ValueError on `ax.set_title()` or `fig.suptitle()`.
  4. **Legends MUST NOT overlap data lines.**  gate raises on detected
     overlap.

What shipped:

- **6 PaperProfiles updated**: font_pt 8 → 10, legend_pt 7 → 10,
  tick_pt 7 → 9 (1pt below body per IEEE/IEEJ convention) for every
  profile.  Margins recomputed for the wider labels:
    IEEE_SINGLE_COLUMN: margin_left 0.155 → 0.165, margin_bottom 0.180 → 0.200
    IGTE_DIGEST_SINGLE: margin_left 0.165 → 0.175, margin_bottom 0.190 → 0.210
    (double-column profiles tightened slightly: 8.5% left → 8% etc.,
     since 10 pt absolute on 18 cm = 0.019 of width vs 0.038 on 8 cm,
     so labels eat proportionally less of the figure)

- **NEW `_check_no_in_figure_title(fig)`** -- walks every `ax.title`
  and the figure-level `_suptitle`, returns the list of non-empty
  titles found.

- **NEW `_check_legend_no_overlap(fig)`** -- for every axis-legend
  pair, samples 200 points along every Line2D and reports any line
  with >= 1 sample inside the legend bbox.

- **`emit_paper_figure(...)` extended**:
    new arg `check_title_in_figure=True` (raises if titles present)
    new arg `check_legend_overlap=True` (raises if legend overlaps)
    Both checks run BEFORE the efficiency gate so the user gets the
    actionable single-line fix first (delete `set_title`, move the
    legend) before being told to also tighten margins.

- **`paper_figure_quality_rules` extended**: 3 new topics
    `font_rule`         -- the 10pt-at-8cm absolute rule + why wider
                            columns don't scale the font
    `no_title_in_figure` -- titles → LaTeX caption, why, override
    `no_legend_overlap`  -- detection + 3 fix recipes ranked by lab
                            preference (direct labels > best_loc >
                            outside-axes)
  `efficiency` topic now defines 余白 = white between axes outer edge
  and figure bbox (was previously fuzzy).

- **`tests/test_graph_paper_figure.py` extended**: 23 → 37 tests
    + test_profile_uses_10pt_body_font (parameterized over 6 profiles)
    + test_paper_figure_rcparams_have_10pt_body_font
    + test_emit_raises_on_ax_set_title
    + test_emit_raises_on_fig_suptitle
    + test_emit_allows_empty_title
    + test_emit_title_check_can_be_disabled
    + test_emit_raises_on_legend_overlap
    + test_emit_passes_when_legend_in_safe_corner
    + test_emit_legend_check_can_be_disabled

Suite: **108/108 pytest pass** (was 94 after 0.78.0 -- +14 new tests).

Why this matters: a default-matplotlib figure exported by an
inexperienced author typically has (a) 8-pt-or-smaller text that the
reviewer can't read at print scale, (b) a title duplicating the
caption, and (c) a legend covering the most interesting curve.  The
three gates now refuse-to-ship each of these.

## 0.78.0 — graph: paper-grade figure scaffolds + efficiency gate

Released 2026-05-26.

The graph subpackage gains a serious paper-quality figure pipeline:
profile-based scaffolds at the journal's EXACT column width, a
measurement gate that refuses to ship wasteful figures, and an
iterative auto-tighten loop that shrinks margins until labels would
clip.

What shipped:

- **NEW: `radia_mcp.graph.paper_figure(profile, nrows, ncols, ...)`**
  -- one-shot scaffold that returns `(fig, axes_2d)` at the journal's
  exact width in mm with pre-tuned subplots_adjust per (R, C) layout
  delta.  Always returns axes as a 2D ndarray so the same loop
  works for any layout.

- **NEW: 6 `PaperProfile`s** (`dataclass(frozen=True)`):
  | Profile | Width | Note |
  |---|---|---|
  | `ieee_single_column` | 88.9 mm | 3.5 in IEEE Transactions single |
  | `ieee_double_column` | 181 mm | 7.16 in IEEE Transactions \figure* |
  | `ieej_single_column` | 88 mm  | IEEJ-D / IEEJ-B 単欄 |
  | `ieej_double_column` | 180 mm | IEEJ-D / IEEJ-B 両欄 |
  | `igte_digest_double` | 170 mm | IGTE / Compumag digest A4 2-col |
  | `igte_digest_single` | 82 mm  | IGTE / Compumag digest single |

  All use 8 pt body font (IEEE-recommended figure-text minimum) +
  Times New Roman serif + Type-42 (TrueType) PDF embedding +
  `xtick.direction='in'` + `units (in parentheses)` lab convention.

- **NEW: `measure_figure_efficiency(fig)`** -- returns the
  axes_area/total_area fraction + per-margin (L/R/T/B) breakdown +
  estimated wspace/hspace.  The metric for the gate.

- **NEW: `auto_tighten(fig, target_axes_fraction=0.80)`** -- iterative
  per-side subplots_adjust shrinker.  Snapshots baseline per-side
  overhang of text artists past `fig.bbox`, then shrinks each side
  by 0.005-step increments and rejects only when overhang grows past
  baseline + (2% width / 3% height) tolerance.  Multi-side per
  iteration (not first-success-only) so a single iter can tighten
  L + B + wspace + hspace together.  Empirically: IEEE 2-col 1x2
  baseline 0.687 → 0.776 (+8.9 pts) without label clipping.

- **NEW: `add_panel_labels(axes, ...)`** -- places (a), (b), (c)...
  at consistent in-axes positions for multi-panel figures.  Bold,
  IEEE convention, with optional bbox.

- **NEW: `emit_paper_figure(fig, path, profile, ...)`** -- the
  validation gate.  Resizes to profile width if needed, measures
  efficiency, then per `on_fail`:
    `'raise'` (default): ValueError + per-margin suggestion of which
                        margin is the biggest waste
    `'warn'`: warnings.warn() and save anyway
    `'auto_tighten'`: run auto_tighten once, re-measure, save
  Saves PDF + PNG at 600 DPI at the profile's exact width.

- **NEW: 3 MCP tools** on `mcp-server-graph`:
    `paper_figure_profiles(query)` -- list profiles with exact mm specs
    `paper_figure_recipe(profile, nrows, ncols, panel_labels)` --
      returns a ready-to-paste Python recipe ending in the
      `emit_paper_figure(..., on_fail='raise')` gate
    `paper_figure_quality_rules(query)` -- the WHY: efficiency,
      margins, units, font_embedding, multipanel
  Total `mcp-server-graph` tool count: 2 -> 5.

- **NEW: `tests/test_graph_paper_figure.py`** (23 tests, ~24s):
  per-profile baseline-in-band locks, auto_tighten gains >= 5 pt on
  IEEE 1x2, no-new-clipping invariant, gate raise/warn/auto-tighten
  behaviour, panel-label placement, profile-width enforcement,
  Type-42 font verification.

- Catalog: `graph` primary_tools extended to 5; meta_health remains
  9/9 PASS, per-server selftest now hits the extended graph
  `--selftest` (paper-quality profiles + recipes + quality_rules +
  runtime smoke of paper_figure / measure_figure_efficiency).

Suite: **94/94 pytest pass** (was 71 after 0.77.0 — +23 new
paper-figure tests).

Why this matters: a typical default-matplotlib figure exported
straight to PDF wastes 30-40% of its area on default outer margins.
For an 8 cm IEEE single-column figure that's ~12 mm of lost axes
width — visible to reviewers as "the curves are tiny".  The gate
makes the wastage refuse-to-ship.

## 0.77.0 — graph subpackage + 4 housekeeping items from review

Released 2026-05-26.

Outcome of an in-conversation `radia-mcp` review (2026-05-26).  The
review found the package healthy (37/37 servers import OK, meta_health
9/9 PASS) and surfaced 4 small housekeeping gaps; all four are
addressed in this release plus the `graph` subpackage migration that
follows.

What shipped:

1. **NEW: `radia_mcp.graph` subpackage** — promoted from
   `s:/mcp-server/src/mcp_server_document/graph/`.  Sugahara Lab
   publication-figure style guide: IEEE / IEEJ font/size profiles,
   MATLAB + Matplotlib snippets, lab style rules (units in parentheses,
   no in-figure title, Times New Roman serif).  Two MCP tools
   (`graph_style_guide`, `graph_size_for_target`) + 10 Python helpers
   (apply_lab_style / lab_figsize / lab_savefig / tighten_margins /
   label_curve_endpoints / add_slope_guide / check_legend_overlap /
   find_best_legend_loc / plot_asymptote_ratio_sweep /
   plot_basis_size_convergence) for direct import.  938-line tools.py
   carried over verbatim with `mcp_server_document` → `radia_mcp` path
   updates in 2 docstring examples.  Catalog count: **37 → 38**.

2. **NEW: LICENSE file (BSD-3-Clause)**.  Closes the SPDX-compliance
   gap where pyproject.toml declared the license but no LICENSE file
   was on disk.  PyPI Warehouse / pip-licenses now see the full text.

3. **NEW: `tests/test_each_server_selftest.py`** (40 tests, ~95s).
   Subprocess-launches every `mcp-server-<x> --selftest` script,
   complementing `test_meta_health.py` which only import-tests.
   Catches: broken `pyproject.toml [project.scripts]`, stale editable
   install after rename, `if __name__ == "__main__"` bugs, cp932
   decode failures in selftest output (decode as UTF-8 with
   errors='replace' on the harness side).  Auto-parameterized from
   the meta catalog so new servers are tested automatically.
   Includes a floor invariant (`test_at_least_30_servers_runnable`)
   that detects the LAB editable-install drift incident pattern from
   CLAUDE.md 2026-05-19.

4. **Catalog alias resolution** (`radia_mcp.meta.catalog`).
   `catalog.get('radia-meta')`, `('radia_meta')`, `('mcp-server-radia-meta')`,
   `('magnetic_materials')`, etc. now all resolve.  Eliminates the trap
   where a user typing the CLI script name into `radia_mcp_get(...)`
   would see "Unknown server".  Adds `_ALIASES` map (auto-populated
   with underscore variants of every hyphenated key) + `_resolve()`
   helper.  `find_related()` is also alias-aware.

5. **Review finding: `panel_review.review_a_panel` is NOT a missing
   `@mcp.tool()`** — it is correctly declared as `@mcp.prompt()`.
   MCP Prompts are a distinct protocol surface from MCP Tools and
   don't show up in a `@mcp.tool` grep.  No change needed; documented
   for future audits.

Suite: 71/71 pytest pass (9 meta_health + 22 chroma_multilingual +
38 server selftests + 2 selftest-harness invariants).
mcp-server-graph --selftest verifies all 7 figure profiles (digest /
paper / presentation / matlab-oversized).

## 0.76.0 — optuna: 5 advanced lab BBO recipes (no Gurobi spinoff)

Released 2026-05-25.

User decision: Gurobi (white-box LP/MIP/QP) is structurally
inappropriate for the lab's black-box FEM-as-objective EM design
problems. `radia_mcp.optuna` stays as the canonical optimization
MCP. Reinforces the locked decision recorded at
`memory/decision_gurobi_dropped_optuna_only.md`.

What shipped:

- NEW `radia_mcp.optuna.recipes_advanced_knowledge` (~625 lines /
  5 topics / ~25k chars) -- complements the existing 5 pattern-level
  recipes in `lab_applications_knowledge` with production-grade
  deep dives that wire Optuna onto an existing Stage-2 calc_*.py
  script.

  | Recipe | Drives | Headline |
  |---|---|---|
  | `pmsm_cogging` | calc_motor_transient.py | NSGA-II multi-obj cogging T_pp + T_avg over magnet alpha_p + slot b_s + skew |
  | `wpt_misalignment` | calc_inductance.py --coil-solver peec | Worst-case eta across 5x3 lateral/vertical offset grid; MedianPruner intermediate reporting |
  | `shielding_layout` | calc_shielding.py | mu-metal / Cu sheet 1-4 placement; Pareto |B| at sensor vs shield mass |
  | `litz_strand_design` | calc_inductance.py --coil-solver peec | n_strands x strand_d x twist_pitch with cost+DC_R pre-filter |
  | `karl_multifidelity` | calc_fem_kelvin.py | Karl iter intermediate_value reporting kills bad geometry in seconds |

- NEW `@mcp.tool() optuna_recipes_advanced(topic)` in
  `optuna/server.py`. Wired into --selftest (6 explicit topics
  + "all" verified > 500 chars each).
- TOOLS.md regenerated: mcp-server-optuna now lists 5 tools
  (was 4: usage / algorithm / lab_applications / status).
- Aliases supported: pmsm, cogging, wpt, misalignment, robustness,
  shielding, shield, litz, strand, karl, pruning_recipe,
  multifidelity, multi_fidelity.

Suite: 9/9 meta_health pytest pass, optuna --selftest PASSED
(27 topics x 4 tools, 25k+ chars on the new tool).

## 0.74.0 — Full CLN corpus absorption (W:\30_CauerLadderNetwork)

Released 2026-05-25.

Complete absorption of the Sugahara lab's **Cauer Ladder Network**
practice corpus at `W:\30_CauerLadderNetwork\` -- ~500 .m / .mph /
.docx / .pdf files across 16 topic folders + 6 root references -- into
the `radia_mcp.mor` subpackage. CLN is the lab's signature MOR
method; Sugahara is co-author on the canonical Kameari-Ebrahimi-
Sugahara-Shindo-Matsuo 2018 IEEE TMAG paper. The corpus was
previously only accessible via direct filesystem inspection; this
release makes it queryable via 5 grouped MCP tools.

**Total new content**: **5238 lines / 57 topics across 5 modules /
213,065 chars of CLN-specific knowledge** -- 5 parallel agents, one
per theme group.

**New MCP tools** (all on `mcp-server-mor`):

| Tool | Source folders | Lines | Topics | Headline content |
|---|---|---:|---:|---|
| `mor_cln_practice` | 01, 02, 09, 2020_11_04 + A-phi.pdf + 2D-rethink + Bessel | 1263 | 12 | Full 71-line `CLN.m` MATLAB class verbatim; COMSOL `HelmholtzEquation(c=0)+withsol('sol2',...)` recursion idiom; Legendre analytical formulas to n=9; Robin/Infinite/Kelvin BC comparison incl. the `Kelvin_NG.m` documented failure mode |
| `mor_cln_multiport` | 03, 04, 10, 11 | 801 | 10 | Kuriyama 2019 multi-expansion `K = C^T nu C + s_0 sigma` with 4 variants (A/T/3D/AK); FreeFEM++ `Multi-turnLadderSeries.edp` quoted; 3D HCurl/H1 saddle-point via `A_phi_Gridap.jl` |
| `mor_cln_advanced` | 05, 06, 07, 14, 16, 2020_12_07 | 1172 | 12 | **FP-CLN** (Fixed-Point CLN; CEFC 2024 Sugahara-Tobita-Matsuo-Takahashi); 4-generation nonlinear lineage 2017-2023 culminating in Tobita's jw method; CLN-as-SPICE-block via Shindo electromagnet 437-line FreeFEM++ driver |
| `mor_cln_specialty` | 08, 12, 13, 15 | 1148 | 11 | **Hiruma method** (Shingo Hiruma, Hokkaido Igarashi -> Kyoto Matsuo): non-symmetric Lanczos producing Cauer ladder from algebraic `(G+sC)x=b`, unifying CLN with PVL/SyPVL/PRIMA; **Nagamine error theory** (Hideaki Nagamine, Kyoto Matsuo): mesh-adequacy rule `delta_n >= 10*Delta_x` from Foster cut-off; BEM+FEM TSVD coupling reducing `O(M*N_m)` to `K=5-15` ports |
| `mor_cln_collab` | 2021_CauerI_to_II, 2022_遠藤, 2023_松本, 2026_長方形, 2017_inverter | 854 | 10 | CauerI vs CauerII (continued-fraction expansion of Z(s) around s=0 vs s=infty); two-matrix Lanczos in K-inner-product (N<=7 stability); Endo @ Hosei 4-square+1-cylinder COMSOL LiveLink sweep; CLN-as-inverter-subcircuit 2017 design memo |

**Wired through**:
- `mor/server.py` -- 5 new `@mcp.tool()` entries + `--selftest`
  exercises each (all 213k chars produced + each `overview` > 200 chars).
- `docs/TOOLS.md` regenerated -- mcp-server-mor now lists 9 tools
  (3 original + 5 CLN deep-dive + 1 status).

**Suite**: 31/31 pytest pass + `mor --selftest` PASSED.

## 0.72.0 — COMSOL fork multilingual RAG absorption

Released 2026-05-25.

Cross-pollination from the upstream wjc9011/COMSOL_Multiphysics_MCP fork
this lab maintains (`ksugahar/COMSOL_Multiphysics_MCP`). The fork added
Japanese / Chinese support to its ChromaDB RAG layer for the COMSOL PDF
manual corpus; radia-mcp's lab corpus at W:/03_文献・論文/00_電磁界解析
is **more multilingual** (roughly 50/50 Japanese textbooks + English IEEE
papers), so the same infrastructure pays off bigger here.

**What shipped**:

- `radia_mcp.common.chroma_retriever.detect_filename_language()` —
  CJK Unicode-range heuristic (cheap, no langdetect / fasttext
  dependency). Returns "ja" / "zh" / "en" / None.
- `radia_mcp.common.chroma_retriever.find_chapters()` +
  `CHAPTER_PATTERNS` constant — multilingual chapter detection:
  English (`Chapter N`, `N.M`), Japanese (`第N章`, `N章`, `第N節`),
  Chinese kanji-numeral (`第一章`, `第十二章`).
- `ChromaRetriever.search(..., language_filter="ja")` — restrict
  semantic hits to chunks tagged with the given language.
- `extract_pdf_chunks(..., default_language=, auto_detect_language=)`
  — tag every chunk's metadata with a language code at index time.
- `literature_index.literature_semantic_search(..., language_filter=)`
  — exposes the filter to LLM clients.
- `literature_index.literature_build_vector_index(...,
  default_language=, auto_detect_language=True)` — defaults to
  filename-based auto-detect for the bilingual lab corpus.
- 22 new tests in `tests/test_chroma_multilingual.py` (filename
  heuristic edge cases + JA/ZH chapter regex + re-export sanity).

**Why this matters**: previously a query like "ヒステリシス測定"
against the full index returned mostly English IEEE papers (more
numerous so they dominate the top-K hits). Adding
`language_filter="ja"` lets a Japanese-language search hit the lab's
Japanese textbook content directly. No re-indexing required for
existing chunks; new builds with `auto_detect_language=True`
(default) populate the metadata tag.

**Total tests**: 31/31 PASS (22 new + 9 existing meta_health).

## 0.69.0 — meta server + uniform tooling + 5 thin-server PDF enrichments

Released 2026-05-24.

**Discovery infrastructure** (the headline change):
- NEW `radia_mcp.meta` subpackage — 36-server cross-server catalog
  (★ recommended first call). Tools: `radia_mcp_overview`,
  `radia_mcp_get(name)`, `radia_mcp_by_tag(tag)`,
  `radia_mcp_related(name)`, `radia_mcp_health`. Entry point
  `mcp-server-radia-meta`. Solves the "which server has knowledge X"
  discovery problem with 3-call lookup instead of guess-and-error.
- NEW `radia_mcp.common.register_status_tool` factory — uniform
  `<server>_status()` introspection (tool list + dep probe + related
  servers). Wired into all 36 servers.
- NEW `radia_mcp.common.register_topics_tool` factory — uniform
  `<short>_topics()` enum for dispatcher-style servers. Wired into 11
  dispatchers (accelerator/bayesian-opt/data-assimilation/electromagnet/
  evolutionary/fusion/gnn/litz-transmission/maglev-linear/pinn/rna-mec).
- NEW `radia_mcp.common` modules: `prompts_loader` (.md knowledge
  loading), `async_runner` (long-running command wrapper),
  `chroma_retriever` (optional ChromaDB+sentence-transformers RAG).

**New subpackages from W:/04_機械学習と最適化 + 99_アプリケーション**:
- 8 ML/optimization: `bayesian_opt`, `evolutionary`, `gnn`,
  `data_assimilation`, `mcmc` (Hokkaido Sato/Yin MCTS lineage +
  Saotome SPM), `optuna` (Sano-Akiba-Imamura textbook), `pinn`,
  `topology_optimization`.
- 19 application + theory: `motor` (ONELAB + Liu Xinyao + Hollaus +
  Wakao + Hane Cauer), `accelerator`, `fusion`, `maglev_linear`,
  `nmr_mri`, `ndt`, `wpt`, `metamaterial`, `magnetic_materials`,
  `litz_transmission`, `rna_mec`, `team_benchmark`, `mor`,
  `matrix_solvers`, `fem`, `bem`, `differential_forms`,
  `mathematica`, `literature_index`.

**Thin-server enrichments** (5 of 6 batch-promoted servers got
substantive PDF-sourced content; +7400 lines total):
- `fusion`: 142→1380 lines, 12 topics (ITER coil system, W7-X
  modular, LHD helical, NbTi/Nb3Sn/HTS CICC, NESCOIL→FOCUS coil
  design, error field, transient eddy, RMP for ELM control).
- `ndt`: 207→1642 lines, 14 topics + 36 aliases (probe types,
  defect models, FEM A-V/T-Omega, MFL pipeline PIG, JSAEM
  benchmarks, ML for NDT).
- `litz_transmission`: 165→1551 lines, 14 topics (Dowell/Wojda/
  Ferreira/Bartoli/Tourkhani M1-M4 taxonomy, Umetani multi-level
  twisting, Igarashi homogenization, Rosskopf FEM+PEEC coupling,
  multiconductor TL).
- `rna_mec`: 177→1591 lines, 12 topics (Derbas 2009 nodal-vs-mesh,
  Lee 2005 TEAM-28 reduced model, Kameari-Ebrahimi-Sugahara-Shindo-
  Matsuo 2018 canonical 3D-FEM CLN, Hane 2020 dynamic hysteresis +
  Cauer MEC, Janet 2004-2005 RNA-MMM mixed method).
- `metamaterial`: 90→1244 lines, 12 topics (Veselago/Pendry/Smith
  LH materials, SRR Pendry LC model, transformation optics with
  explicit Kelvin-inversion cross-link to electromagnet subpackage,
  Sadatgol Bi:YIG+Au 9x Faraday enhancement, Toyota CRLH).
- `maglev_linear`: 173→937 lines, 10 topics, 32k chars
  sources (4 small PDFs only, under 25MB budget per the agent's
  per-request limit): Murata eddy-current demo, Sumitomo Heavy
  patents JP 7-327337 + JP 2007-215264 (PM bearing + planar mover),
  Saiki 2021 Kansai Univ PM maglev thesis. Topics:
  pm_maglev_zero_power (Earnshaw workaround), eddy_current_maglev
  (Hsu-Hill), sumitomo_heavy_industrial, kansai_research,
  lim_lsm_propulsion (Yamamura 1972), scmaglev_eds (Post-Ryutov
  2000), halbach_arrays, end_effects. The 4 open-literature topics
  flag "(open literature, not lab PDF)" — cross-check available
  once the 10-44MB 09_リニアドライブ year-PDFs become accessible.

**Now 6/6 of the originally thin servers got real enrichment**:
total +8300 lines across fusion / ndt / litz_transmission / rna_mec
/ metamaterial / maglev_linear.

**README** (`packages/radia-mcp/README.md`):
- New "## ★ Discovery — start here" section
- Added meta + literature-index rows to the Standalone server table
- Updated JSON config example

**CI**:
- NEW `.github/workflows/radia-mcp-matrix.yml` — Python 3.10/3.11/3.12
  matrix on ubuntu-latest, runs in minutes. Complements existing
  self-hosted Windows `build-test.yml` (45-min full integration).
  Steps: compileall + meta_health + pytest + 36-server --selftest.

**Tests**:
- NEW `packages/radia-mcp/tests/test_meta_health.py` (6 cases):
  importability of all 36 subpackages, catalog floor (≥30),
  status-tool-policy gate (every server must wire
  `register_status_tool`), overview shape, by_tag('optimization')
  finds ≥4, related('mcmc') includes 'optuna'.
- NEW `tests/conftest.py` — resolves `radia_mcp` from this checkout's
  src/ regardless of editable install state.

## 0.55.0 — coordinated bump for radia 4.55.0 (cap-centroid endpoint anchoring)

Released 2026-05-16.  No knowledge changes; coordinated version
release with radia 4.55.0 (rim-end kink fix at lead caps).
See radia CHANGELOG 4.55.0.

## 0.54.0 — coordinated bump for radia 4.54.0 (RMF + corner densification)

Released 2026-05-16.  No knowledge changes; coordinated version
release with radia 4.54.0 (Wang-Joe RMF + corner-densification for
filament viz smoothing).  See radia CHANGELOG 4.54.0.

## 0.53.0 — coordinated bump for radia 4.53.0 (keiko's "arc + leads" coil now PEEC-solvable)

Released 2026-05-16.  No knowledge changes; coordinated version
release with radia 4.53.0 which integrates keiko's CCW winding
fix verbatim and replaces her spine-thinning workaround with
adaptive resampling (policy-compliant).  See radia CHANGELOG 4.53.0.

## 0.52.0 — coordinated bump for radia 4.52.0 (magic-number audit complete)

Released 2026-05-16.  No knowledge changes; coordinated version
release with radia 4.52.0 which adds 4 magic-number pin tests + 2
negative-confidence tests, closing the PEEC STEP-loading audit
started in v4.48.2.  See radia CHANGELOG 4.52.0.

## 0.51.0 — coordinated bump for radia 4.51.0 (Strong Tier C: per-point distance check)

Released 2026-05-16.  No knowledge changes; coordinated version
release with radia 4.51.0 which adds `_check_centerline_near_solid_surface`
(BRepExtrema_DistShapeShape sub-sampled per-point distance check)
as the third orthogonal positive proof in the centerline-verification
chain.  See radia CHANGELOG 4.51.0.

## 0.50.1 — coordinated bump for radia 4.50.1 (PEEC pipeline polish)

Released 2026-05-16.  No knowledge changes; coordinated version
release with radia 4.50.1 (doc lies cleanup, peec_bundle.py
readability fix, 4 magic-number pin tests).  See radia CHANGELOG
4.50.1.

## 0.50.0 — coordinated bump for radia 4.50.0 (Tier C: PEEC STEP-loading sweep complete)

Released 2026-05-16.  No knowledge changes; coordinated version
release with radia 4.50.0 which adds `_check_centerline_inside_solid`
(bbox-containment positive proof) wired into all 5 predicates and
4 filament-construction paths.  Completes the PEEC STEP-loading
weakness sweep started in v4.48.2.  See radia CHANGELOG 4.50.0.

## 0.49.0 — coordinated bump for radia 4.49.0 (Tier A+B+D+E weakness sweep)

Released 2026-05-16.  No knowledge changes; coordinated version
release with radia 4.49.0 which removes silent fallbacks in the
open-spine extractor, adds spine-corner detection at the filament
construction layer (covers HACApK), adds entry guards
(multi-solid STEP raise, cad_to_m silent-1.0/0 bug fix), and adds
8 boundary tests pinning the magic numbers.  See radia CHANGELOG
4.49.0 for the full breakdown.

## 0.48.5 — coordinated bump for radia 4.48.2 (PEEC L fail-fast on NaN)

Released 2026-05-16.  No knowledge changes; coordinated version
release with radia 4.48.2 which adds `_assert_solver_L_finite` to
`peec_bundle.build_bundle_solver` so that silent NaN in the PEEC
mutual-inductance matrix is converted to a hard `ValueError` with
a HINT pointing at vertex-aligned-loft fix path.  See radia
CHANGELOG 4.48.2.

## 0.48.4 — peec_inductance knowledge updated for v4.48.1 STEP-only centerline

Released 2026-05-16.  Pairs with radia 4.48.1 which replaced the
spine-extractor try/except cascades in `coil_from_cad.py` with
classification-based single dispatch + removed the `path_points_m`
parameter ("STEP-Only Centerline: Auto-Detect or Fail" policy in
CLAUDE.md).  The 0.48.3 bump was a version-coordination only; this
0.48.4 ships the actual knowledge updates:

### What shipped

- **`PEEC_IND_FILAMENT_DISPATCH`** (topic `filament_dispatch`): rewritten
  from "3-tier fallback chain" language to **classification-based
  single dispatch**.  Documents Path 1 (UV-map; predicate now
  includes the UV-closure check so downstream sampling MUST succeed
  -- no try/except in Path 1), Path 2 (per-station faces), Path 2b
  (CIRCLE-edge stations), Path 2c (section-planes), Path 3
  (equivalent-circle catch-all with the new fail-fast sanity check).
- **`PEEC_IND_CENTERLINE`** (topic `centerline`): expanded from 3
  paths to **5 classification predicates** (Loft / Circle-edge /
  Revolution-sweep / OPEN longest-edge / CLOSED full-revolution).
  Documents the CLOSED-only guard in `_centerline_from_topology_spine`
  and the keiko `1turn_coil_loft_outsideline.step` lesson (OPEN
  geometries with leads must route to Predicate 4, not 5).
- **New `PEEC_IND_STEP_AUTHORING`** (topic `step_authoring` +
  aliases `cubit_recipe`, `build123d_recipe`, `anti_patterns`):
  concrete recipes for authoring auto-detect-friendly STEPs.
  Quick-decision table mapping Cubit/build123d operations to
  predicate hits, full Cubit `.jou` recipes for gapped torus and
  multi-turn pancake, build123d `sweep()` recipe for curved
  spine + circular profile, anti-patterns (lateral split into 2
  halves, pairwise loft chain, hardcoded IDs, non-manifold,
  self-intersecting), and a 10-line build123d probe script for
  verifying a STEP is auto-detect-friendly BEFORE running the panel.

### Why

radia-mcp 0.48.3 (released 2026-05-15 alongside radia 4.48.1) only
bumped versions for release-triple coordination -- the knowledge
documents still described the obsolete try/except cascade.  Users
asking the `peec_inductance(topic=...)` MCP tool got stale guidance.
0.48.4 reconciles the knowledge layer with the v4.48.1 dispatcher.

## 0.40.0 — 3D CLN (Tanimoto-Kameari) knowledge module

New `radia_ngsolve.knowledge.cln_3d` module captures Tanimoto's 3D
Cauer Ladder Network (CLN) methods from W:/00_CAE/NGSolve/谷本/
master's thesis + production code (~25 notebooks). Covers:

  - **A-T**, **T-Ω**, **A-Φ** formulations (mathematical foundation,
    iteration pseudocode, common boilerplate)
  - **Constraint variants**: penalty stabilization, explicit Coulomb gauge
  - **Solver variants**: SparseSolvPy ICCG, accICCG, NGSolve CG, direct
  - **Validation**: cylindrical TM-mode analytical R/L, Schmidt drift
    diagnostic, bonus_intorder=8 critical setting
  - **Open research note**: Kameari + Kelvin combination remains
    unsolved (3D HCurl A-formulation gives ~25× discrepancy with
    mpmath BEM Foster target due to A_ext gauge unboundedness)

Five canonical notebooks embedded as `cln_notebooks/*.py` resources:
  - `CLN_AT.py` (primary 修論 reference, 7.4 KB)
  - `CLN_T_Omega.py` (T-Ω formulation, 7.6 KB)
  - `CLN_APhi.py` (A-Φ formulation, 8.6 KB)
  - `CLN_2D.py` (2D scalar reference, 2.7 KB)
  - `A_ICCG_production.py` (latest 2024-09-17 production, 6.9 KB)

New MCP tools:
  - `cln_3d(topic="all"|"overview"|"notebooks"|"formulas")`:
    structured documentation
  - `cln_3d_notebook(name="list"|"AT"|"T_Omega"|"APhi"|"2D"|"production")`:
    raw Python code retrieval


## 0.33.5 — Sync with radia 4.10.0 (PEEC-inductance Window merged into IH)

`radia_ngsolve.peec_inductance_knowledge` Source list updated: the
standalone `radia_peec_inductance.py` wrapper was merged into IHWindow
in radia 4.10.0; the analysis is now reached via Method dropdown.
Knowledge text re-points new users at the IHWindow path so MCP
suggestions stay accurate.

No behavioural changes to any MCP tool.

## 0.33.4 — Kelvin knowledge maturity pass (republished)

Same content as 0.33.3 but with a shortened pyproject `description`
field (PyPI's 512-char `summary` limit rejected 0.33.3's metadata
upload at 596 chars, so the wheel never made it to PyPI).  No
behavioural / knowledge changes vs the unreleased 0.33.3; see below
for the actual changes.

## 0.33.3 — Kelvin knowledge maturity pass

Knowledge-only release across 3 subpackages, capturing the
2026-04-26 1/2 + 1/4 Kelvin Benchmark debug session and clarifying
why the 1/8 case has two completely different answers depending on
which panel mode is asking.

- **`radia_ngsolve.kelvin_transformation` (`benchmark_panel` topic)**:
  - Why 1/8 is unsupported for the magnetic-sphere-in-uniform-Hz BVP
    (the source `H0 z_hat` reverses sign under z=0 mirror -- a
    physical limitation, not a Cubit/NGSolve bug).
  - **rho_min sweep diagnostic**: setting rho_min = R collapses
    Mu = mu_0 *(R/rho')^2 to uniform mu_0; if the answer becomes
    correct, the bug is in the Mu coefficient; if still wrong,
    the bug is in BCs / Periodic / mesh.  One solve isolates the
    layer.
  - Surprise: for compact geometry (Kelvin offset = 3*R), even
    Mu = mu_0 in the Kelvin region gives 1/2 +0.34% / 1/4 -0.02% --
    Periodic + sym BCs do most of the open-boundary work.
  - **Cubit-meshed Kelvin needs `-specialcf.normal`** in the
    reduced-Omega Neumann correction term (Cubit assigns surface
    normals with opposite sign to NGSolve's WorkPlane OCC; sign-
    flip A/B test takes 30 seconds and catches it).

- **`cubit` (new `kelvin_reduction_traps` topic)**:
  - Trap 1: `subtract A from B keep` is a silent no-op in Cubit
    2025.3 -- workaround is to drop `keep` and re-create A as a
    fresh primitive.
  - Trap 2: 1/8 octant copy-mesh anchor curve picking is non-
    deterministic (3 equal-length quarter-arcs); fix is
    `min(curves, key=(centroid_z, y, x))` -- 143/143 pairs at
    machine precision.
  - Trap 3: surface normal sign convention differs between Cubit
    and OCC (cross-ref to `radia_ngsolve.kelvin_transformation`).

- **`electromagnet` (new `symmetry_reductions` topic)**:
  - Two distinct Kelvin panel paths -- "Kelvin Benchmark" sphere
    (1/2, 1/4 only) vs "EM panel FEM/MSC" C-yoke (1/1, 1/2, 1/4,
    1/8).  Don't conflate.
  - C-yoke 1/8 sample paths and ELF CEFC 2020 convention
    `ht=0_x, ht=0_y, bn=0_z`.
  - "Don't add a 1/8 sphere benchmark" -- multi-hour debug trail
    capture so the next session doesn't re-investigate.

## 0.32.0 — PEEC-inductance public topic + Cubit daemon speedup

- **`peec_inductance` tool** in `mcp-server-radia-ngsolve`: 5 sub-topics
  (overview / centerline / jou / sibling_jou / japanese_path) promoted
  from LAB-private `mcp-server-ih` after the feature stabilised.
- **Cubit daemon license warmup**: `cubit_license_warmup.py` mirrors
  `coreform_cubit.ps1` renewals cache logic (3-day cache + 7-day
  expiry).  Cold daemon start 30 – 60 s → 3 s.
- **Cubit daemon Phase 1 attach**: per-user stable drop-dir
  (`%LOCALAPPDATA%\radia-mcp\cubit-session\`) + `pid.lock` discovery.
  VSCode restart → new MCP server attaches to living Cubit in
  **0.01 s** instead of re-spawning (6 s cold).
- `open_in_cubit`: same license warmup applied so one-shot GUI
  launches from VSCode also get the speedup.
- `cubit_session_status` reports `mode = owned | attached`.
- New MCP knowledge placement policy in `CLAUDE.md`: stable /
  general → public `radia-mcp` (PyPI), research-stage / lab-only →
  `S:\mcp-server\mcp-server-ih`.

## 0.23.x — YouTube + training pack + GitHub `.jou` search

- **0.23.1** (planned, docs-only): full README rewrite with badges /
  multi-server table / quickstart / lab stance / acknowledgments;
  CHANGELOG.md + CONTRIBUTING.md added. (You're reading it.)
- **0.23.0**: YouTube tutorial transcript scraping for
  `cubit_youtube` / `build123d_youtube` / `gmsh_youtube`
  sub-sources (`youtube-transcript-api` extra). Coreform training
  `examples_only.zip` (24 MB / 30 .jou) auto-folded into
  `cubit_local`. PAT-gated `gmsh_post_jou_github` GitHub-wide `.jou`
  code search. New optional extra `radia-mcp[youtube]`.

## 0.22.x — Universal CAD-MCP mesh backend

- **0.22.4**: lab stance refinement — FreeCAD marked `friendly` /
  `compat — Sugawara Lab respects the FreeCAD community`; build123d
  + Cubit explicitly tagged `主力 (push)` in `lab_policy` topic.
- **0.22.3**: Sugawara Lab primary-pair stance reflected in
  `lab_policy` KB topic + `list_cad_mcp_interop` payload (`lab`,
  `primary_pair` fields) + memory.
- **0.22.2**: build123d marked `PREFERRED` in adapter list, others
  flagged `compat`; `note` clarifies "new lab work should be
  authored in build123d".
- **0.22.1**: expanded CAD detection — `_find_openscad` /
  `_find_freecad` walk Windows `Program Files\FreeCAD*\bin\` and
  macOS `/Applications` so installed-but-not-on-PATH FreeCAD is
  auto-discovered.
- **0.22.0**: new server `mcp-server-radia-interop` —
  `any_step_to_cubit_hex` (universal STEP receiver) +
  `openscad_to_cubit_hex` (CLI) + `freecad_to_cubit_hex`
  (FreeCADCmd subprocess) + `list_cad_mcp_interop`. Position:
  "the mesh backend any CAD MCP can dispatch to."

## 0.21.0 — gmsh community scrape

- New `gmsh_examples(query)` + `gmsh_examples_refresh` MCP tools.
- Sub-sources `gmsh_issues` (gitlab.onelab.info, 3000+ tickets)
  and `gmsh_stackoverflow` (StackOverflow + SciComp.SE `[gmsh]`).
- FAMILIES["gmsh"] union for ranked retrieval.

## 0.20.0 — gmsh post-processing forged

- mcp-server-gmsh-post: bundled auto-generated **gmsh API
  reference** (651 entries across `model` / `view` / `option` /
  `fltk` / …, 2 008 lines, via `_gen_api_reference.py`).
- New cookbooks: `view_data_cookbook`
  (`$NodeData`/`$ElementData`/`$ElementNodeData` decision tree)
  and `physical_groups_cookbook` (dim/tag, downstream solver
  conventions).
- New tools: `gmsh_post_api` (focused tf-idf), `gmsh_post_quality`
  (min Jacobian / skew histogram), `gmsh_post_extract_physical`,
  `gmsh_post_boundary`, `gmsh_post_add_view_from_csv` (most-frequent
  post workflow).

## 0.19.0 — build123d depth gaps closed

- Bundled auto-generated **build123d API reference** (142 classes /
  65 functions / 1 673 lines, via `_gen_api_reference.py`).
- New cookbooks: `plane_axis_location_cookbook` (the 3 most-
  confused classes, 20+ worked recipes) and
  `builder_vs_algebra_rosetta` (side-by-side conversion table).
- New tool `build123d_api(query)` for API-focused tf-idf.

## 0.18.0 — Radia-specific build123d templates + STEP gating

- 7 new templates in `generate_build123d_script`: `magnet_ring`,
  `halbach_array`, `c_core`, `e_core`, `pole_piece`,
  `stator_lamination`, `racetrack_coil`.
- `build123d_inspect_step(path)` — OCCT validity / bbox /
  micro-edge ratio / labels report; gates external STEPs before
  Cubit.
- `build123d_heal(step_in, step_out)` — `OCP.ShapeFix_Shape`
  auto-repair (small edges / face orientation / degenerate fixes).

## 0.17.0 — build123d parity with Cubit

- `lint_build123d_script` + `lint_build123d_directory` (7 rules:
  `missing-buildpart-context`, `sweep-no-path`,
  `polyline-not-closed`, `buildsketch-ambiguous-arg`,
  `missing-export`, `cadquery-in-build123d`,
  `micro-fillet-radius`).
- `build123d_suggest_next(goal, script)` — state-aware (5 goals).
- `generate_build123d_script(pattern)` — 6 starter templates
  (helix_coil, l_bracket, cae_block, gear_bd_warehouse,
  fastener_assembly, sweep_square_path).
- `build123d_try(script)` — fresh subprocess; OCCT segfault
  containment + clean namespace.
- `build123d_to_cubit_hex(script, target_size)` — one-call
  pipeline (build123d → STEP → cubit_mesh_auto → live GUI replay).
- 3 new KB topics: `joints_and_mates`, `assemblies_and_compounds`,
  `cae_workflow_tips`.

## 0.16.0 — Unified search + safety gate

- GitHub PAT auto-discovery (`GITHUB_TOKEN` / `GH_TOKEN` /
  `gh auth token`); 60 → 5000 req/h on GitHub API + GraphQL access.
- Threaded Coreform forum walk (300 topics, ~30 s on 8 threads).
- `build123d_github_discussions` via GraphQL (PAT-gated, 50
  discussions).
- `cubit_ask` / `build123d_ask` unified retrieval across
  bundled KB + scraped examples + optional live web (`include_web`).
- Pre-flight check: `cubit_exec` / `execute_build123d` scan
  failure log for similar inputs (token Jaccard ≥ 0.6) and
  surface the past hint non-blockingly.
- `cubit_mesh_auto` geometry-split rung — auto-detects compound
  bodies (`vol ≤ 3 ∧ surf/vol ≥ 7`) and `webcut volume all with
  cylinder axis z` before retrying scheme auto.

## 0.15.0 — build123d community scrape

- `build123d_discussions(query)` — `gumyr/build123d` GitHub Issues
  + comments (anonymous REST, 60 issues default).

## 0.14.x — gmsh-post lab v4.1 standardization

- **0.14.1**: lint rule `gmsh-v22-deprecated` (HIGH) — flags
  `export mesh "...msh"` without `mesh_version 4.1`. Lab policy
  is v4.1 only; `.vol` (NETGEN native) is the sole exception for
  HO curved meshes.
- **0.14.0**: new server `mcp-server-gmsh-post` —
  `gmsh_post_inspect`, `gmsh_post_validate`, `gmsh_post_convert`
  (lifts any older .msh to v4.1), `gmsh_post_write_node_data` /
  `_element_data` (append `$NodeData` / `$ElementData` blocks
  while keeping the file v4.1-compliant), `gmsh_post_spec`.
  `cubit_exec_safely` — auto-checkpoint to `.cub5`, batch dry-run
  on the snapshot, replay on live GUI on success; silent-error
  detection via `cubit.get_error_count()` delta.

## 0.13.0 — CadQuery interop

- `execute_cadquery(script)` (sibling OCCT lib) +
  `cadquery_to_cubit_hex(script)` one-call pipeline.
- `radia-mcp[cadquery]` extra; integration with cadquery-mcp
  community.

## 0.12.0 — Multi-source example unions

- FAMILY mapping: `cubit` = `[cubit, cubit_local]`; `build123d`
  = `[build123d, bd_warehouse]`.
- `cubit_local` indexer walks `S:\CoreformCubit` (lab archive of
  ~145 .jou) + `S:\Radia\01_GitHub\examples` (~400 files); 753
  files indexed.
- `bd_warehouse` (15 modules: gear, bearing, fastener, flange,
  pipe, …).
- Forum seed queries 5 → 15.

## 0.11.0 — Scraped example libraries

- `build123d_examples(query)` — `gumyr/build123d/examples` (65
  curated scripts).
- `cubit_examples(query)` — Coreform forum (Discourse search.json,
  triple-backtick code-fence extract).

## 0.10.0 — Batch ladder safety pattern

- `cubit_batch_try(commands)` — disposable headless Cubit.
- `cubit_mesh_auto(step_path)` — scheme ladder
  (auto → sweep → polyhedron → tetmesh) batch-validated, winning
  recipe replayed in live GUI. 4-turn spiral coil yielded 1668
  hex on first run.

## 0.9.0 — Failure log + tf-idf retrieval + live web docs

- Persistent jsonl failure log per kind (`cubit` / `build123d`),
  fed into every `*_lookup`.
- tf-idf retrieval with heading boost replaces substring counter.
- `cubit_web_docs` (Discourse JSON for forum.coreform.com) +
  `build123d_web_docs` (readthedocs).

## 0.8.0 — Standalone wheel crystallized

- Plan A established (Cubit GUI + PyQt5 QTimer + file-drop IPC).
- `cubit_session.py` dual-mode (gui / batch) + auto-restart on
  RPC failure.
- `cubit_checkpoint(label)` / `cubit_restore(label)` — `.cub5`
  snapshot undo.
- `cubit_mesh_diagnose` (per-volume scheme alternatives),
  `cubit_suggest_next(goal)` (state-aware), `cubit_lookup(query)`
  (heading-chunk retrieval over 8000-line knowledge).
- 4-turn coil + KEIKO 6-letter text both produced pure hex
  meshes via the build123d → Cubit pipeline.

## 0.5 / 0.6 / 0.7 — Initial wheel

- Standalone `radia-mcp` package extracted from the `radia` core
  repo (Option Y restructure).
- `mcp-server-cubit` and `mcp-server-build123d` as the first two
  entry points.
- OCP CAD Viewer retired in favor of the persistent Cubit GUI.
