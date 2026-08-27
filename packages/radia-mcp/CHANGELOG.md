# Changelog

All notable changes to `radia-mcp`. Format: each release lists **what
shipped** + **why** in compact form. Older releases (≤ 0.4) are
omitted; the 0.5 → 0.6 jump is when the standalone `radia-mcp` wheel
crystallized as its own package.

## [Unreleased]

## [1.4.46] - 2026-08-27

- cubit: added a reusable Nastran consumer gate that verifies the exact BDF
  digest, dimension, element order, counts, bounding box, retained imported
  mesh, material assignment, and optional set semantics without confusing a
  downstream importer limitation with an exporter defect.
- cubit: synchronized executable knowledge with the solver-neutral
  `export nastran_bdf` command and retained `export jmag_nastran` only as an
  explicit legacy-journal compatibility path.

## [1.4.45] - 2026-08-27

- streamfunction: documented and exposed the checked ACA QR-TSVD and Optuna
  outer-loop workflow used by the production Stream Function Simulink block.
- release: synchronized the executable Stream Function knowledge and tests with
  the Radia 4.95.64 H-matrix accuracy and tuning contract.

## [1.4.44] - 2026-08-26

- matlab: kept standalone `optuna_mex` commands outside the regular Radia MEX
  inventory while preserving complete checked coverage of the Radia gateway.
- release: documented the narrow `radia-optuna` recovery lane that reuses only
  an immutable successful push-CI artifact for the exact tagged SHA and
  rechecks both required CI jobs before publication.

## [1.4.43] - 2026-08-25

- radia-ngsolve: migrated deterministic two-dimensional Netgen mesh writers
  to the 6.2.2606 `EdgeDescriptor` contract and taught the strict `.vol`
  inventory and execution gates to accept `edgesegmentsgi3` plus
  `edgedescriptors` while retaining legacy `edgesegmentsgi2` input.
- matlab: classified the private directional-Schur proposal helper explicitly
  so the checked pybind/MEX parity inventory remains complete without
  presenting a Python optimizer implementation detail as a public MATLAB API.

## [1.4.42] - 2026-08-22

- matlab: narrowed all Optuna claims to the checked 4.9.0 subset and exposed
  that boundary through the compatibility contract. The direct oracle now
  covers `TrialPruned`, trial-number `tell`, FrozenTrial-returning collections,
  FrozenTrial public behavior, and the snake-case `create_trial` factory.

- grant-writing: added the non-scoring
  `grant_writing_adjacent_reviewer_readability_check` for prose that remains
  hard for an adjacent-domain reviewer despite short sentences. It reports
  compressed concept density, method-name piles, paragraphs that mix science,
  decision, and infrastructure layers, and assurance language distributed
  across the draft. It also catches a research answer collapsed into an
  implementation representation, vague relation or decision objects, and a
  required scope that names applications but no deliverable. Reviewer-
  vocabulary checks now reject descriptions of MCP as a storage location, and
  the argument-evidence map locates whether preparation evidence is connected
  to a research item the team can start. The integrated report exposes the
  readability and traceability diagnostics without adding them to the defect
  score or corpus finding baseline. Subject-predicate distance checks now also
  recognize the full-width Japanese comma `，` used by Word documents.

- grant-writing: generalized a project-specific wording prohibition into a
  central-question semantic contract. Reviews now distinguish the scientific
  answer (condition), its operational decision criterion, and its calibrated
  application scope or limit; `boundary` is not banned when a variable, two
  states, and an estimation method are defined. Exact wording and one- versus
  two-sentence locks remain local to each proposal instead of becoming
  universal lint rules. The central-claim checker and guidance now state this
  scope explicitly.

- grant-writing: added the source-grounded `grant_writing_kaken_review_axes`
  reference for current KAKENHI Scientific Research B/C (General). It keeps
  the three research-plan elements, the separately rated internationality
  element, and budget validity distinct; records the two-stage review and
  three-reviewer structure for Scientific Research C; and carries official
  JSPS URLs, revision/access dates, budget itemization requirements, and the
  consequence of multiple budget-problem ratings. `kaken_generic` now checks
  internationality as well as the three research-plan axes.

- grant-writing: corrected two corpus-exposed budget false positives against
  the official FY2027 entry guide. A future conference whose venue is not yet
  announced is no longer treated as an unfilled form field. Budget necessity
  prose may state amounts when it gives a recomputable basis such as unit
  price times quantity/duration, a quotation, or an official tariff; only a
  bare amount is now reported. This replaces the overly broad 2019 rule that
  all amounts must live in the table.

- grant-writing: review-axis presence now runs on applicant prose after form
  instructions are removed. The funder's own request to discuss academic
  importance, method, capability, or internationality can no longer satisfy
  the corresponding applicant-content check.

- grant-writing: expanded the private outcome corpus from 8 to 19 documents
  (7 adopted, 11 rejected, 1 unsubmitted) across ordinary KAKENHI, Go-Tech,
  JSPS invitation, Power Academy, foundation, and internal-grant genres.
  Outcome-labelled entries can now retain immutable submitted sources, a
  written classification basis, and award/review evidence. Added
  `sweep.py --compare-outcomes`, which normalizes findings by prose length and
  counts document-level prevalence while explicitly refusing to treat outcome
  as a score. In the 10 ordinary KAKENHI documents, detector density was 5.55
  per 10k characters for adopted drafts and 4.37 for rejected drafts: another
  direct negative result for adoption prediction.

- grant-writing: separated ordinary KAKENHI review (`kaken_generic`) from the
  current OSS-platform theme (`kaken_oss`). Historical KAKENHI applications
  now receive the general review axes without being judged against
  project-specific GitHub/MCP vocabulary. Wrapped ethics examples in older PDF
  forms are also removed before the human-subjects check, preventing the
  funder's own アンケート/インタビュー instructions from becoming an
  applicant finding.

- grant-writing: added a non-scoring `argument_evidence_map` for locating the
  central question, prior gap, method operation, decision rule, knowledge
  output, preliminary evidence, responsibilities, and negative-result plan.
  It provides excerpts and manual review prompts rather than a keyword score.
  Public Japanese-lint wrappers now all apply the same form/LaTeX prose filter;
  acronym parsing accepts English-gloss definitions, Japanese particles after
  acronyms, and hyphenated project names such as JP-MARs.

- grant-writing: the private corpus lane can now read an ordered `paths` list
  for a live multi-file proposal, and records a SHA-256 source fingerprint in
  its adjudicated baseline. The current 科研費 draft therefore follows its
  four TeX sources directly instead of a stale extracted-text snapshot, and
  any later prose edit forces a new review even when finding counts happen to
  stay unchanged. While connecting the live draft, the review-format check
  was also fixed to ignore `\newcommand` metadata in the final-year field;
  its repeated 「該当なし」 values are not bare human-rights rationales.

- grant-writing: audited which checks ever say anything about a real
  proposal (`sweep.py --audit`). Eight were silent on all eight documents;
  adjudicating them found one genuinely broken.
  `literature_gap_evidence_check` was applicable to **nothing**: it was built
  for a search report (確認できなかった / 見当たらなかった) while real
  proposals assert absence outright -- 「統合的なマルチスケールモデル縮約法が
  存在しない」, 「直接的な競合製品は存在しない」 -- with no account of how
  they looked. It now reports `absence_claimed_without_search` for that, the
  more common and more punctured form: a reviewer needs one counterexample.
  This is the second time a vocabulary here was written from assumption and
  matched zero real proposals, after `_GAP_MARKERS`; derive the words from
  documents. Three of the four instances are in adopted proposals, so it
  predicts nothing about adoption either.
- grant-writing: `check_misuse_japanese` had no test at all. It inherits the
  shared Japanese table aimed at speech and email (よろしかったでしょうか,
  のほう, こんにちわ), which no research proposal trips, so its silence is
  the genre and not a fault -- now recorded, and covered by a test that
  proves it still fires on the text it was built for.

## [1.4.41] - 2026-08-20

- presentation: `presentation_script_vs_slide_coverage` could not read
  Japanese prose, so every acronym in a deck came back as never spoken --
  a report that names specific words, looks precise, and is confidently
  false. `\b` does not work next to a CJK character (Python's `\w` matches
  Japanese, so "のHACApK" has no word boundary before the H): a slide puts
  the acronym in a box of its own where a boundary exists, a script says it
  mid-sentence where one does not. Also fixed: deck furniture such as the
  affiliation in the corner of every slide was scored as unspoken content
  (tokens on 80%+ of slides are now set aside, and listed in the report so
  the exclusion is auditable); exact segmentation was required, so saying
  "Gram二重積分" where the slide says "二重積分" counted as saying less; and
  the cover slide, whose formal title need not be read aloud, is no longer
  scored. Measured on a real deck: average coverage 46.9% -> 58.9%, with the
  false "HACApK / Jacobi / Richardson are never spoken" verdict gone and the
  remaining flags genuine.

## [1.4.40] - 2026-08-19

- presentation: recorded the script-first authoring format. The skill has
  said since Phase 2 that the script comes first and the slide is built to
  support it; there was nowhere to put the script, so the talk lived in one
  file and the deck in another and the checks that compare them
  (`presentation_script_vs_slide_coverage`,
  `presentation_estimate_per_slide_time`, `presentation_speaker_note_ratio`)
  had nothing reliable to read. One Markdown file now holds both halves: a
  heading starts a slide, a `>` blockquote is what you SAY and becomes the
  speaker notes, everything else is what the audience sees, and
  `radia.equation.markdown_to_pptx` builds the deck with the equations native
  rather than pictures. The header's "AI cannot make the slide body" caveat is
  retired with it -- that was true when a deck could only be built by hand.

- meta: catalogued `new-src-module-unclassified-in-parity-manifest`. A module
  added under `src/radia` without a rule in the MATLAB parity manifest leaves
  `ci_preflight` green and CI red, because the preflight's top-level gate is
  collect-only: the new file imports fine, it is simply unaccounted for.
  Adding a module is exactly the change that gate cannot see.

## [1.4.39] - 2026-08-19

- matlab: mapped the four beam Lie-map / orbit pybind names onto their
  `beam.*` MEX commands so the parity audit sees the coverage that the
  shared C++ kernels already provide (98/98 public names, 362 gateway
  commands), and recorded the one real gap: the MATLAB orbit tracker
  drives Radia-object sources only while the HDiv iron evaluator remains
  a pybind-owned handle.

- grant-writing: encoded the R9 (FY2027) in-house KAKENHI call briefing.
  New `grant_writing_kaken_review_format_check` catches color-only figure
  discrimination (some categories are reviewed as monochrome prints),
  missing safeguards when surveys / animal experiments / personal data
  appear, a bare 「該当なし」 without a rationale in the human-rights/legal
  box (the box reviewers flag most often), unidentifiable publication
  mentions in the researchmap era, and an incomplete funding-overlap box
  (相違点・応募理由・所属組織役職); full drafts are additionally checked
  for three-review-criteria coverage and emphasis/figure use (reviewers
  read up to ~100 proposals in about a month). The integrated health
  report now runs the check for every program, `skill.md` gained
  review-reality / program-strategy (充足率, 重複応募, 開拓の2段階書面
  審査化, 補助金 vs 基金) / compliance-layer (research integrity via
  e-Rad, funding-overlap disclosure incl. foreign and private funds,
  effort definition, DMP, immediate open access, equipment sharing,
  cost-item rules) sections, and the budget policy now states the ~70%
  award rate of 基盤 categories versus full funding of 挑戦的研究.

- md2html: replaced the hand-rolled math-protection regexes with
  `pymdownx.arithmatex` (new `pymdown-extensions` dependency in the
  `[md2html]` extra) and repaired seven defect classes found by a full
  audit: bare `\begin{align}` blocks lost their `\\` row separators;
  `$` inside inline code was mistaken for a math opener; `[N]` inside
  code/math became citation links; `||A|B||` stayed unconverted and
  split table rows; the References-section detector anchored on the
  first prose occurrence of the word (numbering procedure lists instead
  of the bibliography); and `<` inside math was injected raw, so the
  browser swallowed the following text as a tag. Pipes inside math
  spans are now neutralized as `&#124;` before the table parser runs,
  and the citation pass skips code, math, and existing anchors.
  `tests/test_md2html.py` grew from 12 to 21 goldens locking each
  repaired behaviour.

## [1.4.38] - 2026-08-18

- Added PPTX figure auditing for authored width, paste scale, final visible
  text size, clipping, and duplicate slide images, plus deterministic PDF-page
  rasterization for presentation workflows.
- Added presentation integrity checks for raw math markup and embedded figure
  text, together with a repair path that preserves placement and reports
  unsupported crops or rotations explicitly.
- Synchronized accelerator and MATLAB parity knowledge with the exact
  equivalent-current vector-potential source, reflection-conditioned field
  diagnostics, and the CanonicalHCurl initialization boundary.
- Added the headless OpenCV dependency to the document extra used by figure and
  presentation inspection without introducing a desktop GUI dependency.

## [1.4.37] - 2026-08-16

- Corrected the presentation message hierarchy: slide titles use the shortest
  practical concrete "target + viewpoint" phrase (for example, "Model 1
  accuracy evaluation"), while the bottom line states what was learned from
  the slide evidence. Added `presentation_check_slide_title_specificity` and
  retained the former title-verb tool name as a compatibility alias. The title
  check now uses a five-axis, 10-point acceptance rubric: concrete target,
  explicit viewpoint, concise one-line form, no result claim, and deck-level
  uniqueness.
- Re-audited the 12-book writing corpus under the lab's `09_作文技術`
  shelf and documented routing boundaries between paper, figure,
  presentation, and grant guidance.  Added the Kinoshita goal-statement
  workflow and `paper_writing_check_conclusion_first_use`, which distinguishes
  legitimate new synthesis in a conclusion from technical terms, symbols,
  numbers, and citations introduced there without body support.
- Operationalized Miyano's Japanese slide-copy guidance with
  `presentation_check_japanese_copy_style`: the PPTX audit now identifies
  hard line breaks that split dependent phrases, estimates the body-text
  nominal-ending ratio, and flags screen-only lead-ins such as
  `そこで本研究では`.
- Strengthened presentation message hierarchy: each content slide should state
  its most important claim in the title and place a one-sentence interpretation
  of the evidence above the footer.  Added
  `presentation_check_slide_message_hierarchy` to audit both requirements.
- Raised the presentation audience-content floor to 24 pt for body text,
  captions, annotations, tables, and chart labels.  Titles remain 32 pt or
  larger; footer, page-number, date, and source chrome are explicit exceptions.
- Strengthened the 16:9 figure profiles: `beamer_169_full`,
  `beamer_169_half`, and `presentation_slide` now use a 24 pt minimum for
  labels, ticks, legends, panel labels, and annotations.  The export gate
  rejects smaller per-artist overrides before saving.
- Added the electromagnetic force/torque MCP family and synchronized its
  Lorentz, Maxwell-stress, air-gap, coenergy, and virtual-work contracts with
  the Radia Python and MATLAB APIs.
- Extended accelerator knowledge for the HCurl EarlyTimes Bishop/RMF loft
  chain, fifth-degree canonical Hamiltonian, independent reference-curvature
  control, and separately reported Lie, A-map, and direct B-map errors.

## [1.4.36] - 2026-08-14

- Added accelerator-magnet knowledge for distributed transfer-matrix and
  nonlinear-map attribution directly over solved electromagnetic fields,
  including Enge fringe integrals, edge focusing, and multi-momentum FFAG
  diagnostics.
- Added Gmsh electromagnetic post-processing for charged-particle trajectories
  and flying-beam views, with bounded inputs and public-safe artifact metadata.
  Field comparison, flux and line integrals, Maxwell force, gap harmonics, and
  Poincare reductions are now registered first-class tools.
- Extended Cubit/Gmsh contracts for topology-optimization material fractions
  and the integrated Radia LTspice operating surface. A closure miss now names
  both the Sculpt lattice and the source design mesh without prescribing a
  misleading automatic refinement.

- Added an electromagnet knowledge contract for Enge fringe-field form factors:
  `I1` reference-orbit displacement, `I2` effective edge focusing, pole-face
  geometry, and momentum-indexed soft-edge diagnostics for FFAG magnets.  The
  guidance distinguishes these integrals from Enge profile coefficients and
  scopes HDiv-MMM as air-volume-mesh-free rather than vacuum-computation-free;
  it also names the executable multi-momentum target, native-row, and
  whole-element optimization APIs.

## [1.4.35] - 2026-08-08

- Added six Gmsh workflow tools for cross-file field ranges, comparable panel
  rendering, compound element selection, slice-stack volume views,
  dense-streamline flow textures, and file-series temporal statistics.
- Added named camera, colour, glyph, clipping, axis, and annotation controls to
  PNG and animation rendering, with shared transfer functions propagated to
  generated slice views.
- Replaced compound-selection `eval` with a bounded AST parser/evaluator and
  protected coordinate/function names from field-name collisions.
- Made file-series statistics verify node coordinates, connectivity, view
  section/component contracts, finite values, and valid scalar output; unsafe
  ElementNodeData reduction now fails with explicit conversion guidance.
- Made multi-panel comparisons require one physical quantity and fail loudly
  when an automatic shared range or shared zoom cannot be honored.
- Documented the isosurface diagnostic as an outer-boundary contact check
  rather than a general topology proof.

## [1.4.34] - 2026-08-07

- Strengthened grant-writing budget review so external costs must be traceable
  to dated official prices or quotations, including tax, purchase-unit,
  expiry, currency, rounding, direct-cost ceiling, and future-travel status.
- Added grant-writing checks that require platform proposals to end in a
  field-specific, falsifiable knowledge product; require proposal-specific
  metrics to separate calibration from held-out validation; and require
  cross-organization pilots to identify the transferred artifact, independent
  action, observed result, and remaining gap.
- Added a grant-writing reviewer-vocabulary check that expands OSS/AI terms by
  their Japanese role, prefers readable institution and domain terminology,
  and prevents named benchmarks from substituting for engineering significance.
- Added a grant-writing persuasion check for self-negating evidence, equation
  on-ramps and symbol definitions, post-equation interpretation, defensive
  prose, optional side branches, acronym density, and internal memo shorthand.
- Extended the persuasion check to flag short inline conditions such as a bare
  symbol-equals-zero claim before the proposal explains its physical meaning
  and engineering consequence.
- Added a grant-writing abstraction check that keeps named software out of
  titles, questions, aims, novelty, and impact while preserving concrete names
  for reproducible methods, preliminary evidence, collaboration, rights, and
  costs.
- Added a grant-writing evidence-scope check that flags field-wide adoption,
  causal-barrier, and academic-gap claims inferred from non-detection in a
  bounded literature corpus, with health-report integration and rewrite
  guidance.

## [1.4.33] - 2026-08-06

- Added a headless Netgen/Cubit mesh-quality comparison using one Gmsh minSICN
  referee, with bounded inputs and per-route failure isolation.
- Added STEP/BREP and post-data overlays to Gmsh PNG/GIF tools, including
  shaded CAD controls and robust line-view rendering.
- Added executable mesh-quality studies and contracts for output orientation,
  high-order connectivity, PEEC current views, and Cubit session cleanup.
- Regenerated the tool inventory from the integrated MCP source so unrelated
  editable-install work does not leak into the published catalog.

## [1.4.32] - 2026-08-06

- Added the executed electromagnetic field-line and adaptive-isosurface
  notebook as the public Gmsh post-processing showcase.
- Linked Gmsh knowledge to saved WebGUI scenes and synchronized quantitative
  evidence without publishing generated heavyweight render artifacts.

## [1.4.31] - 2026-08-06

- Fixed Cubit GUI bootstrap imports under Cubit's string-execution path and
  consolidated probes and snapshots across GUI and batch transports.
- Made private sessions fail loudly unless their process is protected by a
  checked Windows kill-on-close Job Object, with deterministic handle and
  private-directory cleanup.
- Made snapshots write through a unique temporary PNG, verify non-empty output,
  and atomically replace the destination so stale captures cannot pass.
- Added 15 ParaView-parity Gmsh tools covering derived fields, thresholds,
  symmetry expansion, transforms, skins, regular-grid and curve resampling,
  CSV/histogram/history export, montages, and camera orbits.
- Added adaptive RK4 field tracing with closed-loop detection, exact equal-flux
  contours from planar or axisymmetric potentials, and Jobard-Lefer evenly
  spaced streamlines on arbitrary plane patches.
- Added bounded adaptive isosurface extraction for high-order views and smooth
  surface normals by default; PNG/GIF MCP tools now expose numeric/string view
  options, clipping, adaptation, smoothing, and linked-view controls.
- Made field-line controls reject non-finite or non-physical values before
  launching Gmsh, enforce the 2D integration-step and spacing contract, and
  stop exactly at the configured total-step budget.
- Established physical 1:1:1 axis scaling for spatial Gmsh figures, with an
  explicit warning whenever a render requests an exaggerated axis.
- Hardened the pure-Python MSH data reader against malformed section counts,
  duplicate or missing tags, invalid row widths, tag ranges, and connectivity;
  post plots also follow the repository's no-in-figure-title rule.
- Completed build123d/Cubit structured error-kind coverage and made Cubit block
  element counts deduplicate direct and geometry-based membership.
- Added tool-classification, call-log rotation, Job Object, malformed-MSH,
  stale-snapshot, concurrency, and live GUI plus headless Cubit coverage,
  including explicit Popen stream and private-directory cleanup.

## [1.4.30] - 2026-08-06

- Added a persistent headless Gmsh session with MATLAB MCP-style execution
  verbs, option control, post-processing, probing, profiles, integration,
  derived views, cuts, isosurfaces, harmonics, and streamlines.
- Added structural and field-aware MSH audits, strict NaN/Inf and parser
  gates, Gmsh `minSICN` shape-quality checks, and overall plus per-step field
  diffs that detect changes hidden by unchanged extrema.
- Expanded Cubit session setup, doctor, journal, call-log, and lifecycle
  diagnostics, including collision-free file-drop requests for multiple MCP
  clients attached to one daemon.
- Hardened PEEC MSH v4.1 handling and synchronized the generated public tool
  inventory with the current Gmsh and HDiv-VIM surfaces.

## [1.4.29] - 2026-08-05

- Added headless GMSH v4.1 inspection, validation, PNG/GIF rendering, and
  invalid-option lint tools, plus Cubit `check-vol` and entity/label probes.
- Ported the official MathWorks MATLAB MCP execution patterns and documented
  the MagLev and electromagnet-topology Simulink production interfaces.
- Kept solver-backed `check-vol` numerical checks in the validation lane and
  refreshed the generated public tool inventory.

## [1.4.28] - 2026-08-04

- Added a public-safe visual differential-geometry knowledge and validation
  layer connecting intrinsic/extrinsic geometry, Gauss-Bonnet, holonomy,
  differential forms, and the H1/HCurl/HDiv/L2 de Rham sequence to Radia
  electromagnetic workflows.
- Expanded the generated MATLAB/MEX contract to 334 commands with complete
  121/121 numerical-class coverage, including configured ChargeGram material
  element and candidate-cluster inspection.
- Enforced the package's provenance and internal-path boundary in the local
  pre-push matrix, matching the GitHub policy-lint workflow before publication.

## [1.4.27] - 2026-08-02

- Added fail-closed Motor evidence tools for BDM1 HEX torque, armature
  reaction, and absolute demagnetization attribution, with strict mesh,
  angle-grid, metric-range, and cross-artifact identity checks.
- Updated MATLAB/MEX release knowledge for the Level-2 induction-heating
  object-handle lifecycle and expanded the generated command inventory to the
  complete current native surface.

## [1.4.26] - 2026-08-01

- Added a public-safe accelerator-fundamentals layer to
  `radia_mcp.electromagnet`, covering magnetic rigidity, Twiss/dispersion and
  integrated-field handoff, magnet-family selection, normal-conducting and
  rapid-cycling design, superconducting conductor/quench constraints, field
  measurement, commissioning, and explicit beam-physics scope boundaries.
- Added a searchable bibliographic guide for the 12-source, 1,844-page
  accelerator textbook corpus used by that layer, with topical/page locators
  and no machine-local paths or redistributed PDF text.

## [1.4.25] - 2026-07-31

- Aligned the release with Radia 4.95.31's native ChargeGram streaming
  lifecycle fix, preventing concurrent deformation modes from reentering one
  native parent after successive NGSolve TaskManager regions.

## [1.4.24] - 2026-07-31

- Published the native Motor angle-family MEX/Simulink lifecycle contract,
  including split output/update semantics, Custom SimState, source hashes, and
  live MATLAB evidence.
- Added content-addressed Motor dual-lane artifacts for an honestly scoped
  smooth annular first-harmonic fixture, with aligned transverse excitation and
  fail-loud shared-identity checks.
- Extended release knowledge and uninstall-safety evidence to bind the Motor
  solver artifacts, native Simulink tests, and fresh dependency scans before a
  legacy solver can be removed.
- Raised the MATLAB MEX command inventory to 316 and documented the standalone
  MEX debugging surface alongside the production S-Function path.
- Recorded compute-host provenance for field and native Motor evidence while
  keeping optional Optimization Toolbox diagnostics explicit on lean hosts.

## [1.4.23] - 2026-07-30

- Added executable MCP contracts for scalar, thermal, harmonic magnetic,
  circuit, force, transient, post-processing, and periodic air-gap field
  studies, including owned workers and checked production artifacts.
- Added reversible FEMM/AGE and axisymmetric-signature migration evidence,
  uninstall-safety checks, and protocol-level validation without making a
  legacy solver installation part of the production contract.
- Published the typed MATLAB/Simulink material, winding, and Field Study
  interfaces, with strict `.vol` boundary ranges, portable region labels, and
  a collision-free Simulink library knowledge index.
- Corrected harmonic-study observables to the documented RMS-phasor convention
  and expanded numerical power-closure coverage.
- Made production-evidence hashes independent of platform line endings and
  restored dependency-light FEM/Motor server imports through lazy solver loads.

## [1.4.22] - 2026-07-29

- Added the verified hodograph free-boundary design guidance for saturable
  pole faces and flux-concentrator horns, including the constructive cap
  guarantee, practical applicability conditions, and explicit limits against
  treating the construction as a general optimizer.
- Added exact axisymmetric ring-current and point-potential source contracts,
  with fail-loud vertex and Dirichlet checks, NumPy/generator normalization,
  preserved nonlinear positional compatibility, and linear/nonlinear
  constant-reluctivity validation.
- Added signed-DOF reduction, mixed-boundary, residual, and dual-boundary
  helpers for axisymmetric models, with persistent nonlinear constraints and
  strict finite/integer input validation.

## [1.4.21] - 2026-07-29

- Constrained the MCP Python SDK to the supported 1.x FastMCP contract after
  the incompatible 2.0 package removed `mcp.server.fastmcp`; the minimal CI
  matrix now installs the same declared dependency range.
- Published the isochronous HDiv-MMM topology-optimization contract, including
  staged feasibility restoration, exact-void/SIMP promotion, and measured
  machine-level interpretation without machine-local provenance.
- Added the conducting-sphere analytical reference and the named MATLAB
  fallback ownership required by the checked Python-to-MATLAB parity surface.

## [1.4.20] - 2026-07-28

- Added executable URN guidance for Y-admittance, Cauer-ladder, and CLN-peeling
  workflows, including passivity checks and the distinction between stored-grid
  reconstruction and look-ahead termination.

## [1.4.19] - 2026-07-24

- Added JSON-only TEAM Workshop Problem 36 contract and validation tools. The
  cross-validation gate now resolves named 250 s observables from the supplied
  Radia artifact and requires an identity-matched independent reference.
- Published the analytic-adjoint MMA/SQP and native HCurl topology contracts,
  including the rule that directional finite differences are QA only.
- Updated the executable MATLAB guidance for automatic multivariate TPE,
  table-persistent full-covariance CMA-ES, and the Stream Function Simulink
  optimization interface.

## [1.4.18] - 2026-07-22

- Published the Simulink-only application policy, MATLAB-only IH sample-object
  contract, checked GMSH result-artifact rules, and strict `.vol` preflight in
  the executable MATLAB, Cubit, and panel-review guidance.
- Extended document-meta auditing for public CAE examples to require executed
  WebGUI scenes and explicit `Draw(field, mesh, name=..., ...)` field views.
- Retired the remaining IH notebook-workbench guidance while preserving
  result-bearing docs notebooks as the public reproduction and visualization
  layer.

## [1.4.17] - 2026-07-21

- Published the acoustic, axifem, and HCurl-topology Python/MEX equivalence
  gates in the MATLAB executable contract and locked the documented MEX
  inventory to 311 commands and 232 covered pybind entries.
- Strengthened the minimal-dependency acoustic server test across BDF1 and
  BDF2 while explicitly proving that neither NumPy nor Radia is imported.
- Added the checked Python-to-MATLAB capability manifest and native-promotion
  backlog to the MCP-visible MATLAB workflow guidance.

## [1.4.16] - 2026-07-21

- Added the Radia Simulink application-library contract and synchronized MCP
  guidance with the Simulink-first production policy and temporary IH dual
  operation.
- Documented the native NGSolve MATLAB MEX surface, HCurl topology workflow,
  Python/MEX parity gates, and official MATLAB MCP Server execution boundary.
- Added Radia acoustic and acoustic FEM-BEM server contracts, strengthened
  grant-writing and catalog validation, and removed machine-local test-data
  dependencies from the public package checks.
- Moved 68 solver-backed FEM/BEM, geometry, convergence, and application
  checks to `validation_test/radia_mcp/`; the package test lane now enforces
  fast API/MCP contracts without direct Netgen/NGSolve imports.
- Kept the Radia acoustic MCP import and selftest usable in the minimal MCP
  environment by evaluating its BDF1/BDF2 CQ grid gate with the Python
  standard library instead of importing NumPy or Radia at startup.

- Strengthened `mcp-server-mathematica` for reusable validation workflows:
  tracked `.wls` / `.wl` / `.m` files can now run as one-kernel batches with
  elapsed time and parsed JSON failure reports, and named identities can be
  checked together without paying one Wolfram kernel cold start per formula.
  The server now also exposes an electromagnetics, differential-forms, and
  paper verification guide, and the meta catalog points to the real tools.

## [1.4.15] - 2026-07-17

- Added motor guidance for the angle-periodic ROM, reduced HCurl/HDiv basis
  construction, native C ABI export, energy accounting, and stateful
  hysteresis integration.
- **Merged the figure server into paper-writing** (Sugahara: figure was a
  shared middle-layer server for paper-writing + presentation; since
  presentation merged into paper-writing on 2026-07-17, figure no longer
  needs to be a shared middle layer). `mcp-server-paper-writing` now serves
  all `figure_*` / `paper_figure_*` tools; the `mcp-server-figure` entry
  point and the standalone-server scaffolding in `radia_mcp/figure/server.py`
  (module-level FastMCP, auto-register loop, `register_status_tool`, `main`)
  were retired -- `radia_mcp.figure` remains the implementation home and now
  exposes a `register(mcp)` used by the paper-writing server. The `figure`
  catalog entry was folded into `paper-writing` with `figure` /
  `mcp-server-figure` -> `paper-writing` discovery aliases; `docs/TOOLS.md`
  regenerated. The standalone `figure_status` tool is dropped (paper-writing's
  status tool covers the merged tools).

- **Merged the presentation server into paper-writing** (Sugahara: slide
  decks cannot yet be authored end-to-end by AI, so the slide lint / PPTX
  toolset does not warrant a standalone server). `mcp-server-paper-writing`
  now serves all `presentation_*` tools; the `mcp-server-presentation`
  entry point and `radia_mcp/presentation/server.py` were retired (the
  `radia_mcp.presentation` module remains the implementation home, still
  consumed by `document_meta` cross-lint). Catalog entry folded into
  `paper-writing` with a `presentation -> paper-writing` discovery alias;
  skills cross-linked (`paper_writing/skill.md` gained 発表スライド +
  `radia_mcp.figure` operational sections; `presentation/skill.md` header
  notes the merge).

## [1.4.14] - 2026-07-17

- Documented genus-1 loop-current closure across weak and strong induction-
  heating coupling, including screening diagnostics and the canonical current-
  phase convention.
 - Added method-selection guidance based on the weak-route coil back-reaction
   ratio and synchronized the topology references with the production API.

## [1.4.13] - 2026-07-17

- Added executable Mathematica studies for RT/BDM simplex families,
  parent-order admissibility, conductor-topology reduction, and the discrete
  de Rham bridge from HCurl T modes to solenoidal HDiv currents.
- Documented the topology-aware eddy-bubble HCurl response basis, SIBC and
  conductor-cycle retention, and the shared-mesh BDM-MMM coupling API.
- Corrected HDiv family terminology to BDM1/BDM2 throughout the live API
  knowledge and added focused family-dimension and production-contract tests.
- Strengthened CQ contour, ACA pivot, LTspice event-window, and group-delay
  evidence gates, and synchronized IH guidance with the removed invalid
  incident-potential and loop-current routes.

## [1.4.12] - 2026-07-17

- Added genus-aware induction-heating guidance for surface orientation,
  cohomology loop-current closure, surface-Poisson incident potentials, and
  strong-coupling current normalization, grounded in analytic sphere and ring
  contracts.
- Strengthened evidence-lineage gates across build123d, Cubit, COMSOL, FEMM,
  MATLAB, NGSolve, motor, RF, inverse-problem, and SPICE workflows so reported
  quantities remain bound to their active mesh, basis, frame, units, source,
  operating point, and artifact generation.
- Extended counterfactual coverage for mixed-mesh exports, nonlinear
  inductance, force/coenergy, motor loss, regularized traces, rotational eddy
  braking, PWM loss, and CAD assembly gates.
- Documented the NGSolve BDM/RT family distinction and the topology-aware
  HCurl/BDM hybrid VIM, EVRS/T-method, and reduced-response workflow.

## [1.4.11] - 2026-07-16

- **HDiv-VIM production completeness**: documented the persistent solver,
  BDM2 nonlinear and quadrature-state hysteresis paths, independent-body
  coupling, multi-body Radia dispatch, and focused validation contracts.
- Updated HDiv capability guidance and public API knowledge for BDM1/BDM2 field
  reconstruction, restart state, coupled permanent-magnet/iron workflows, and
  accelerator workbench order selection.

## [1.4.10] - 2026-07-16

- **Cubit mixed-mesh evidence gates**: require two-sided manifold transition
  interfaces when recorded, verify every required export artifact is fresh,
  nonempty, digest-bearing, and uniquely named, and lock replayed journal/model
  identities to the pinned source evidence.
- Preserve older summaries as readable evidence with explicit warnings when
  per-artifact freshness or replay identity metadata was not recorded.

## [1.4.9] - 2026-07-16

- **Mapped BDM2 HDiv-VIM**: documented the supported 2D, tetrahedral,
  hexahedral, and wedge order/geometry combinations and the curved-mesh
  compute-host validation route.
- **Validation-learning gates**: retained counterfactual evidence and added
  focused curriculum gates for geometry handoff, nonlinear energy balance,
  force-method selection, inverse regularization, and minor-loop behavior.

## [1.4.8] - 2026-07-15

- **Motor HDiv-VIM lane**: documented the planar reduced-motor analysis,
  notebook workbench route, torque cross-checks, and focused validation targets.

## [1.4.7] - 2026-07-15

- **HDiv-VIM production stack**: documented the fixed/given magnetization,
  linear-recoil permanent-magnet, native energy-stop hysteresis, persistent
  matrix, field-evaluator, BDM2, and IMA interfaces with focused knowledge tests.

## [1.4.5] - 2026-07-13

- **radia-ngsolve gates**: added heterogeneous P1 current-flow,
  transformer inductance identity, and sphere `.vol` convergence gates
  with focused pytest coverage.
- **motor gates**: added electrothermal result-handoff and replayed
  demagnetization-history validation gates.
- **build123d gate**: strengthened external solid validation by rejecting
  zero-volume curved shell STEP inputs before Cubit handoff.
- **FEEC / hysteresis validation**: locked HDiv charge-Gram
  loop-eigenvalue protection and B-input loop-pollution regression
  records.

- **motor / AGE validation routing**: added public-safe `motor_age_quality`
  and `motor_age_validation_plan` tools, plus the ELF/MAGIC motor bridge
  and field quick-check router.  Motor prompts now route from SPM/IPM/IM/SRM/
  SynRM/hysteresis intent to explicit NGSolve AGE gates, physical quantities,
  pytest targets, and publication labels without exposing private solver
  provenance.
- **CI / optional dependencies**: the baffled-piston acoustics regression now
  skips when scipy is unavailable, matching the package's lightweight
  selftest policy.  This patch supersedes the failed 1.4.2 tag run.
- **BREAKING**: removed the in-package `radia_mcp.optuna` server and
  `mcp-server-optuna` entry point.  Optuna operation now belongs to the
  official public `optuna/optuna-mcp`; radia-mcp keeps the CAE objective,
  Bayesian/surrogate, evolutionary, and topology-optimization knowledge
  layers that pair with that external server.  Install the official
  server separately with `pip install --upgrade optuna optuna-mcp` and
  run it as `optuna-mcp` (optionally with
  `--storage sqlite:///C:/temp/optuna_mcp.db` for persistence).
- **radia-ngsolve (dtn_coarse_mesh, accuracy)**: tightened the `p`-method
  claim after an adversarial methodology review. "mode `n` exact iff order
  `p ≥ n`" is precise only in the **reference space** / 2-D; on the curved
  3-D sphere `p ≥ n` removes the *polynomial* error but the realized
  accuracy then **floors at the curved-geometry + conformal-weight error**
  (`n=2`→`2.5e-3`, `n=3`→`8.4e-4`, both →`~1e-5` as order rises = ~5–6
  digits — *exactly Kameari's observation*, not machine-exact; 2-D floors
  deeper, `~1e-7…1e-9`). New "REALIZED floor" section + softened headline
  wording across `p_method` / `api` / overview, `fem_bem_coupling` and
  `server` docstrings, the docs page, and the manuscript. Also added the
  missing source-factor numbers (`c_5/c_1 = (4/15)(a/R)^4`, `n≡1 mod 4`;
  optimal `R/a ~ 3`) and explicit **domain-of-validity caveats**: the
  `Σ c_n·defect_n` factorisation needs spherical truncation + linear
  materials + interior source + interior/closure error separation; Kelvin's
  `defect_n` is a *discretisation* error (→0) vs PML/Robin's *model* floor;
  "formulation-independent" qualifies the operator eigenvalue, not the
  discrete defect. The 8 verified experiments are now routed through
  `docs/kelvin/ARCHIVE_RETIREMENT.md` to the maintained docs/API/validation
  surfaces.
- **radia-ngsolve**: new `dtn_coarse_mesh` tool + `bem_integral` /
  `fem_bem_coupling` measurements that reframe Kameari's coarse-mesh
  accuracy of the Kelvin transformation as a **DtN-matrix spectral
  property**. The exterior Dirichlet-to-Neumann operator `Λ_ext` has the
  closed-form ladder `−(n+1)/R` (3D) / `−n/R` (2D); the discrete `Λ_h`
  lands the low multipoles on that ladder on the coarsest mesh (dipole
  0.07 %), and the **isolated** Kelvin open-BC error (~0.1 %) sits ~45×
  below the interior FEM error. New code: `exterior_dtn_spectrum`,
  `dtn_spectrum_vs_mesh` (`bem_integral`); `kelvin_dtn_eigenvalue`
  (`dim=3`/`dim=2`), `kelvin_vs_exact_open_bc_error`,
  `kelvin_openbc_error_vs_exterior_mesh`, `kelvin_twosphere_shell_dipole`
  (`fem_bem_coupling`). Knowledge module `dtn_coarse_mesh`, example
  `dtn_spectrum_coarse_mesh_demo.py`, tests `test_dtn_spectrum_coarse.py`
  + extended `test_fem_bem_coupling.py`, and the academic reference
  `docs/kelvin/DTN_SPECTRUM_COARSE_MESH.md`. Supports the IEEJ
  static-apparatus / rotating-machinery joint technical meeting.
- **radia-ngsolve (fix)**: `laplace_fem_bem_schur` corrected — the
  boundary Schur term is `K_FEM − Pᵀ M_bnd Λ_ext P` (three bugs fixed:
  the `K +` → `K −` sign, the missing `SurfaceL2` mass `M_bnd` that turns
  the coefficient→coefficient `Λ` into a weak-form contribution, and the
  free-DOF restriction of the densified `V` solve). Also fixed a
  nondeterministic NGSolve teardown crash: a `definedon`-restricted
  `SurfaceL2` used as the **test** space corrupts the heap on a volume
  mesh — the coupling matrix now uses it as the **trial** space.
  Validated `< 6 %` L2 on the spherical shell at order 2.
- **BREAKING**: renamed the `radia_mcp.graph` subpackage to
  `radia_mcp.figure` (server `mcp-server-graph` -> `mcp-server-figure`;
  tools `graph_style_guide` / `graph_size_for_target` /
  `graph_matlab2tikz_recipe` -> `figure_style_guide` /
  `figure_size_for_target` / `figure_matlab2tikz_recipe`).  Clean
  rename, no backward-compat shim.
- **figure**: added two beamer 16:9 talk profiles
  (`beamer_169_full` 150 mm, `beamer_169_half` 72 mm) — 10 pt body
  (on-page 10 pt @ 8 cm, same as paper), heavier strokes for projector
  legibility, NO in-figure title (the beamer frametitle carries it).
  Times New Roman + Okabe-Ito remain the defaults.
- **presentation**: added reference-citation tools for talks —
  `presentation_cite_format` (talk / bracket / numeric / full-IEEE from
  a BibTeX entry or explicit fields), `presentation_references_slide`
  (beamer `thebibliography` frame + plain numbered list for a PPTX text
  box), `presentation_add_citation_footer` (insert a Times-New-Roman
  footnote textbox on a .pptx slide), and `presentation_citation_audit`
  (flag dangling `[N]` and never-cited references).
- **figure (policy)**: codified two figure conventions in the quality
  rules. (1) The no-in-figure-title rule applies to ALL figures —
  including the `beamer_169_*` slide profiles (`emit_paper_figure`'s
  no-title gate fires for them too; the slide's title is the beamer
  frametitle / caption). The beamer slide profiles use **10 pt** (not
  11) to honour the on-page rule. (2) Made the font rule explicit:
  on-page **10 pt at 8 cm = 1.25 pt/cm**; matplotlib authors at the
  embed size (→ 10 pt @ 8 cm), MATLAB authors oversized (→ **20 pt @
  16 cm**, a 2× downscale to the 8 cm column).
- **figure (API)**: misuse-proof slide/paper figure API, forged from the
  CEFC-2026 incident (figures authored at 150 mm then
  `\linewidth`-downscaled to ~6 pt, carrying in-figure titles, at risk of
  a silent DejaVu fallback). `lab_figure(embed_width_cm)` authors AT the
  on-page width with verified Times New Roman; `save_lab_figure` runs
  FAIL-LOUD gates (no in-figure title, TNR actually used pre+post, font
  embedding, no CJK) and returns the exact `\includegraphics[width=<cm>]`
  snippet (embed at 100% → on-page 10 pt). Adds one-call builders
  (`scaling_loglog` / `grouped_bars` / `convergence` / `quiver_pair` /
  `bh_curve`), `_assert_times_new_roman` (catches the silent DejaVu
  fallback), and a `figure_audit_embeds(tex)` MCP tool that lints a .tex
  for height-constrained / `\linewidth` / DejaVu-embedding figures.
- **figure (layout)**: `legend_no_overlap(ax)` places the legend by the
  standard escalation and never on top of the data — `loc='best'`, then the
  exhaustive best-of-six in-axes search (`find_best_legend_loc`), then (if a
  curve still passes under it) OUTSIDE-right with the axes shrunk to fit the
  SAME figure width.  `save_lab_figure(tighten=...)` now pushes the axes box
  to the limit within the FIXED figure size (`auto_tighten`, overhang-safe so
  labels never clip) and reports `axes_fraction` — the supported
  "bbox ぎりぎり" path, NOT `bbox_inches='tight'` (which would change the
  width and break 10 pt @ width).  Builders use `legend_no_overlap` by default.

## 0.99.1 — radia-ih knowledge sync + 3 new bug-patterns

Released 2026-06-02.

- `ih` / `radia_ngsolve` knowledge: removed the dangling `calc_heating.py`
  references (the orphan was deleted in radia 4.89.1); reframed to the
  EM -> q_surf -> Thermal flow.
- `meta.bug_patterns`: +3 entries
  (pardiso-mkl-thread-dll-fails-in-pytest-subprocess,
  stale-index-lock-in-shared-clone,
  policy-lint-helmholtz-hodge-false-positive); 20 total.

## 0.99.0 — version bump alongside radia 4.89.0 (loop-free-by-default)

Released 2026-06-02.  Bumped in lockstep with the radia 4.89.0 triple
(Helmholtz-Hodge loop removal); see the radia CHANGELOG.

## 0.98.0 — Cubit knowledge updated for the `export` command verb

Released 2026-06-01.

The cubit / radia_ngsolve knowledge now documents the renamed Cubit
plugin commands (`export netgen / gmsh / vtk / femeem / meg` + `export
jmag_nastran`) instead of the old `radia_export` verb, matching
cubit-mesh-export 0.11.0 / radia 4.88.0.

## 0.97.0 — SF-coil RegularizedTSVD + single-stroke field_aware + sheet-metal distort + bug-pattern catalog

Released 2026-05-31.

Accumulated `radia_ngsolve` knowledge + meta tooling since 0.96.1, shipped
alongside `radia` v4.86.0:

- **`aca_tsvd`**: RegularizedTSVD closed form (`ψ = S⁻¹V·W⁻¹·Σ⁻¹·UᵀB`) + Path-A
  cache + Optuna 3-mode (RMS / constrained / Pareto); single-stroke chain
  `field_aware` (beats kuijpers: 9.3% vs 16.2%); and the new single-current
  **sheet-metal coil distortion (bankin-ho)** technique — `--distort` bends one
  series wire in 3D (control-grid VectorH1-style deformation) to drive a planar
  uniform-Bz coil from ~12000 ppm to **340–2015 ppm** with ONE current.
- **`bug_patterns`** (radia-meta): learned bug-pattern catalog + `bug_patterns_lookup`
  / `bug_patterns_stats` MCP tools (learn-once-from-incidents policy).
- knowledge cleanup: scrub stale hardcoded `2025.3` / `6.2.2603` paths.

## 0.93.0 — multi-file .tex resolver + abstract auto-extract + E2E real-paper validation

Released 2026-05-26.

User directive: "残課題をクリアしてからpypi公開＋100号機とmdxにデプロイ
だね。" (clear the 3 residual paper-writing gaps then PyPI publish +
deploy to 100号機/mdx).

Closes 3 known gaps that surfaced during v0.92.0 selftest review:

  1. `em_submission_gate` could only see a SINGLE .tex file -- any
     paper that used `\input{chapter1}` style modular structure had
     half its checks running on an essentially-empty main file.
  2. The gate forced the caller to pass `abstract_text=...` manually
     because there was no way to extract the abstract body from a
     .tex source.  Real submissions have the abstract in the source
     already; manual re-typing was pure friction.
  3. The strengthening passes (v0.88-0.92) never ran end-to-end on
     an accepted real lab paper -- only on synthetic test fixtures.

### NEW module `_tex_resolver.py` (~340 LOC) -- 2 tools, 2 helpers

  1. `resolve_input_chain(main_tex_path, max_depth=8, encoding="utf-8")`
        - Recursively inlines `\input{X}` and `\include{X}` into a
          single merged string.
        - Walks the LaTeX-canonical search semantics (relative paths
          resolved against INCLUDING file's directory, not main file's).
        - max_depth=8 cap prevents pathological cyclic includes from
          blowing the stack.
        - Detects + skips duplicates (LaTeX would re-typeset; for
          static analysis a single pass is correct).
        - Detects + skips `% \input{X}` commented lines.
        - Returns metadata: files_resolved, files_missing,
          files_duplicate, total_chars.

  2. `extract_abstract_from_tex(tex_source)` -- 4-pattern extraction:
        - `\begin{abstract}...\end{abstract}` (most common)
        - `\begin{IEEEabstract}...\end{IEEEabstract}` (IEEE conf)
        - `\abstract{...}` (IEEEtran legacy)
        - `\textbf{Abstract}\\\\ ... \end{quote}` (AAAI-style)

  3. `paper_writing_resolve_input_chain(main_tex_path, max_depth=8,
      return_merged_text=False, encoding="utf-8")` -- MCP tool.

  4. `paper_writing_extract_abstract(tex_path, encoding="utf-8",
      auto_resolve_inputs=True)` -- MCP tool.  Falls back to
      `resolve_input_chain` if abstract not found in main file --
      handles the case where the abstract is in a separate
      `\input{abstract}` subfile.

### `em_submission_gate` extended (transparent automation)

  Two new parameters (both default True, fully backward-compatible):
  - `auto_extract_abstract=True`: if `abstract_text` is empty AND
    `tex_path` is supplied, automatically extract via the 4-pattern
    matcher and run the abstract-only checks on it.
  - `auto_resolve_inputs=True`: if `tex_path` is supplied, run
    `resolve_input_chain` first.  When `\input` subfiles are found,
    write a merged temp file and route the 4 tex-based checks
    (forward_reference / equation_numbering / count_underlines /
    undefined_variables) through the MERGED text.

  Two new check categories appear in the gate report:
  - `multifile_resolved`: pass + "resolved N \input subfiles" OR
    skip + reason.
  - `abstract_extracted`: pass + character count + source file OR
    skip + reason.

### Real-paper E2E validation

  Ran the full gate on a real 18-year-old Japanese lab thesis
  (cp932 encoding, 7 chapter subfiles, 71KB merged):
  - multifile_resolved: pass (7 \input chains inlined, 0 missing)
  - figure_forward_reference: pass
  - equation_numbering: pass
  - count_underlines: pass
  - undefined_variables: caught 40 real undefined math symbols
  - text_image_overlap: caught 4 real text-on-figure overlaps
  - page_whitespace_anomalies: warned about 80 mostly-blank pages
  All checks ran cross-file, the resolver+extract architecture works
  on a real (not toy) document.

### Tool count

  v0.92.0: 89 paper_writing_* tools + 1 prompt
  v0.93.0: 89 + 2 = **91 paper_writing_* tools** + 1 prompt

### Tests

  +10 tests covering:
  - resolve_input_chain: inlining, missing-subfile graceful,
    commented \input skipped
  - extract_abstract: standard env, IEEEabstract env, not-found
  - paper_writing_extract_abstract: direct path, via input-resolution
  - em_submission_gate: auto-extracts abstract, resolves multifile tex

### Verification

  - Full in-process pytest: 258 + 10 = **268 tests pass**
  - mcp-server-paper-writing --selftest: all sections OK
  - Real-paper E2E: gate ran 12 checks on cp932 7-chapter thesis

## 0.92.0 — pixel-accurate overlap/overflow detection + undefined-variable check

Released 2026-05-26.

User directives:
  1. "web 検索で、agentic paper writing の場合に、画像出力 or PDF 出力
     の品質チェック（重なりやはみだし）をどうやって検出しているかを
     調べて。"
  2. "未定義の変数をチェックすることも追加してください。"

Web-research basis (2026-05-26):
  * arXiv:2106.00676 -- VILA (Shen-Lo-Wang) visual layout groups
  * arXiv:2604.05018 -- PaperOrchestra multi-agent paper writing
  * arXiv:2512.02589 -- PaperDebugger plugin in-editor checker
  * arXiv:2604.01128 -- Paper Reconstruction Evaluation (Presentation
    + Hallucination axes)
  * pymupdf Rect.intersects() + Rect.intersect() bbox primitives
  * pdfminer.six char_margin / line_overlap grouping algorithm
  * misc0110/paper-linter -- generic LaTeX lint (for symbol-check
    comparison)

### NEW module `_pdf_overlap_detection.py` (~450 LOC) -- 4 tools

  1. `paper_writing_detect_text_image_overlap(pdf, iou_threshold=0.02,
      min_intersection_area_pt2=50.0, ignore_full_page_images=True)`
        - Walks every text block bbox vs every image bbox via
          pymupdf.get_image_info + get_text("dict").
        - Flags IoU + intersection-area dual threshold.
        - Ignores full-page journal-template images by default.
        - Catches: caption-on-figure, full-image-reappearing-under-text
          (clip+viewport trap), page-number-on-header-logo.

  2. `paper_writing_detect_text_overflow_page(pdf, use_cropbox=True,
      overflow_tolerance_pt=1.0)`
        - Compares every text bbox to CropBox (or MediaBox).
        - Flags > 1pt past any side.
        - Catches: wide \\verbatim, wide table, full-width
          \\includegraphics in single-column figure, long inline eq.
        - Complements paper_writing_check_overfull_hbox (reads .log) --
          this catches OVERFLOW that LaTeX itself didn't warn about.

  3. `paper_writing_detect_overlapping_text_blocks(pdf, iou_threshold=0.05)`
        - Pairwise text-text IoU on every page.
        - Catches: caption straddling page break, two-column
          collision, negative \\vspace pushing paragraphs together.

  4. `paper_writing_pdf_overlap_recipe()` -- 4-step decision tree +
     web-research citations.

### NEW module `_undefined_variables.py` (~350 LOC) -- 1 tool

  `paper_writing_check_undefined_variables(tex_path, extra_whitelist,
   report_first_use=True, max_reported=50)`

  Detection algorithm (3-stage regex):
    1. EXTRACT symbols from 8 math environments (equation / align /
       eqnarray / gather / $...$ / \\(...\\) / $$...$$ / \\[...\\])
    2. BUILD defined-set from:
       * Nomenclature blocks (\\begin{nomenclature} /
         \\begin{IEEEdescription} / \\section{Notation|Symbols})
       * "where ..." clauses after equations
       * Inline definitions ("Let $X$ denote ..." / "$X$ is the ...")
    3. SUBTRACT a UNIVERSAL_WHITELIST of 80+ universally-recognized
       symbols (\\pi, \\omega, \\mu_0, t, f, x, y, z, E, B, H, ...)

  Returns ranked list of symbols never defined, with first-occurrence
  source context (~100 chars).

  Lab POLICY (matches em_paper_style 'reviewer_patterns' #11):
  every math symbol used MUST be defined in Nomenclature OR "where"
  OR inline.  Universal whitelist exempts conventional symbols.
  extra_whitelist parameter for project-specific generic indices.

### em_submission_gate extended with 3 new checks

  `paper_writing_em_submission_gate()` now runs (when respective
  inputs supplied):
  - `undefined_variables` (tex_path): FAIL if any math symbol undefined.
  - `text_image_overlap` (pdf_path): FAIL if any overlap.
  - `text_overflow_page` (pdf_path): FAIL if any text past page edge.

  All three default to status="fail" since these are reviewer-visible
  defects.

### Tool count

  v0.91.0: 84 paper_writing_* tools + 1 prompt
  v0.92.0: 84 + 5 = **89 paper_writing_* tools** + 1 prompt

### Tests

  +17 tests covering:
  - PDF overlap: text-on-image flagged (synthetic PDF), clean-PDF
    case, text-overflow flagged (text past right edge), clean case,
    text-text overlap clean, recipe content, IoU math (4 corner cases)
  - Undefined variables: clean paper with where-clauses (0 undefined),
    \\eta_s undefined flagged, IEEEdescription nomenclature exempts
    defined symbols, universal whitelist exempts \\pi/\\omega/\\sigma/f,
    extra_whitelist runtime arg, missing tex returns error,
    first-occurrence context retrieved

### Verification

  - mcp-server-paper-writing --selftest: 14 sections OK
  - Full in-process pytest: 244 + 14 = **258 tests pass**

## 0.91.0 — citation-verification policy: reference.bib + search-and-verify enforcement

Released 2026-05-26.

User directive: "paper-writing では、reference.bib を使うことと検索
して裏を取ることを忘れずに。"

The #1 AI-assisted-paper failure mode is **citation hallucination** --
a plausible-looking DOI / author / year that does NOT correspond to
any real paper.  Reviewers catch every single one; one hallucinated
cite taints the whole submission.  This release ENFORCES the lab
policy at three levels:

### Level 1: NEW behavioral recipe

  `paper_writing_citation_workflow_recipe()` -- 6-step mandatory
  workflow with strong language ("NEVER invent", "ALWAYS verify").
  Read this BEFORE generating ANY \\cite{} or BibTeX entry.

  The 6 steps:
  1. READ the user's reference.bib first (single source of truth).
  2. SEARCH for grounding via Crossref / Semantic Scholar / arXiv.
  3. MATCH candidates against reference.bib (already-cited? skip).
  4. VERIFY the DOI resolves (defensive double-check).
  5. INSERT into reference.bib only if "ready_to_insert" verdict.
  6. LINT after insertion (paper_writing_lint_reference_format).

  Plus explicit "when verification cannot complete" guidance:
  emit \\cite{TODO: verify -- ...} marker, ask the user; do NOT
  fabricate a placeholder DOI / author / year.

### Level 2: NEW composite verification tool

  `paper_writing_verify_citation(claim, bib_path, candidate_doi,
  candidate_arxiv_id, candidate_title, search_arxiv_if_no_doi=True)`

  Returns dict with verdict one of:
  - **`"found_in_bib"`** -- already cited; reuse `matching_key`.
    DOI matching is case-insensitive + handles `https://doi.org/`
    prefix.  Title matching is whitespace-normalised + first 80
    chars substring (~90% precision in practice).
  - **`"ready_to_insert"`** -- DOI resolved via Crossref;
    `suggested_bibtex` is the ready-to-paste entry.
  - **`"needs_disambiguation"`** -- arXiv search returned multiple
    candidates; user picks.
  - **`"no_candidate_found"`** -- DOI did not resolve OR arXiv
    search empty.  Emits TODO advice; does NOT fabricate.
  - **`"error"`** -- bib read failed / bib_path missing.

  Composes the existing tools (resolve_doi + doi_to_bibtex +
  semantic_scholar_lookup + arxiv_search) into a single check.
  Lightweight bib parser (regex; no pybtex dep) reads existing
  citations.

### Level 3: submission_gate enforcement

  `paper_writing_em_submission_gate()` now FAILS when bib_path is
  not supplied -- previously this was a soft skip.  The new
  `bib_policy` check status="fail" message:

    "reference.bib was not supplied.  Lab POLICY: every citation
     must come from the user's actual .bib and be verified via
     paper_writing_verify_citation BEFORE insertion.  Re-call
     with bib_path=/path/to/reference.bib OR explicitly justify
     why no .bib check applies."

### NEW MCP prompt

  `cite_a_claim(claim, bib_path="reference.bib")` -- surfaces the
  policy as a prompt the AI sees when asked to insert a citation.

### Tests

  +8 tests covering:
  - recipe contains the 8 policy keywords (NEVER invent, ALWAYS,
    reference.bib, Crossref, Semantic Scholar, arXiv, verify, TODO)
  - verify_citation no bib_path → error verdict
  - verify_citation missing bib file → error verdict
  - DOI already in bib → found_in_bib verdict (case + URL handling)
  - Title already in bib → found_in_bib verdict
  - Insufficient info → no_candidate_found verdict
  - DOI resolve failure → no_candidate_found (does NOT fabricate)
  - arXiv search fallback (mocked) → needs_disambiguation verdict

### Tool count

  v0.90.0: 82 paper_writing_* tools + 0 prompts
  v0.91.0: 82 + 2 = **84 paper_writing_* tools + 1 prompt**

### Selftest update

  mcp-server-paper-writing --selftest now exercises:
  - citation_workflow_recipe size + key-strings
  - verify_citation no-bib refusal (offline)
  - em_submission_gate FAILS on missing bib_path (enforces policy)

## 0.90.0 — paper_writing production hardening (test suite + EM style + submission gate)

Released 2026-05-26.

User directive: "paper-writing の mcp-server は徹底的に鍛えないと。"
(thoroughly strengthen paper-writing to production grade).

**3-pass hardening**:

### Pass 1: Comprehensive test suite

  Until v0.89.0 there was ZERO dedicated test coverage of the 80
  paper_writing_* tools (only the server's --selftest smoke check).
  **NEW `tests/test_paper_writing.py`** locks 65+ test cases across:

  * Text-style (count_weak_expressions, analyze_sentences,
    paragraph_opener / paragraph_length / sentence_ending_variety,
    word_repetition, tense_consistency, prose_density, passive_voice)
  * Abstract / IMRaD (validate_abstract_length, background_ratio,
    imrad_balance, abstract_strength)
  * Figure / equation (forward_reference, equation_numbering,
    figure_caption_showing on minimal IEEEtran .tex)
  * Citation / bibliography (lint_reference_format,
    check_citation_usage, self_citation_ratio with mini .bib)
  * JA-lint (kanji_ratio, notation_variants, find_undefined_acronyms,
    acronym_usage_audit, lint_bedrock, misuse_japanese,
    suggest_redundancy_fixes, subject_predicate_distance)
  * PDF (validate_pdf_pages + check_pdf_edge_overflow with a
    pymupdf-generated minimal valid PDF)
  * v0.88.0 layout (tex_figure_placement aliases, layout_visual_recipe)
  * v0.89.0 external sources (normalize_arxiv_id round-trip,
    extract_equations on synthetic 6-env input, mocked
    arxiv_fetch_latex_source + arxiv_search + semantic_scholar_lookup
    + references + citations via monkeypatched requests.get)
  * Health / orchestrator (health_report, adaptive_health_report,
    run_full_workflow, root_cause_diagnosis, next_5_actions,
    rewrite_suggest)
  * Reviewer (classify_reviewer_comment, generate_response_letter,
    generate_cover_letter, reviewer_2_trigger_summary)
  * Plan-B tier scores (contribution_clarity, claim_quantification,
    limitation_statement, related_work_density,
    figure_referencing_coverage, journal_fit_assessment)

### Pass 1 bonus: real bug uncovered + fixed

  `_ja_lint.py::grant_writing_lint_bedrock` referenced `_scan_hedges`
  but the import was stripped during the v0.88.0 AST-extraction
  inline-copy from grant_writing.tools.  Restored
  `from ._shared.hedges import scan_hedges as _scan_hedges` --
  lint_bedrock now works.

### Pass 2: NEW `_em_paper_style.py` (~700 LOC) -- EM-domain knowledge

  Generic text/lint is journal-agnostic; this module captures the
  conventions EM reviewers (IEEE TAP/TMag/TMTT, IEEJ Trans D/B, IGTE)
  actually enforce.  6 topic sections + index:

  * `sign_conventions` -- engineering exp(+jwt) vs physics exp(-iwt)
    + conjugation rule for porting between conventions.
  * `vector_tensor_notation` -- bold-italic E vs arrow vs underline;
    tensor styles; operators upright vs italic; differential-forms
    cross-link.
  * `b_vs_h` -- "magnetic flux density B" vs "magnetic field
    intensity H" terminology; M (A/m) vs J (T).
  * `si_units` -- 4 ironclad rules; EM-specific unit checklist;
    siunitx package recommendation.
  * `equation_typesetting` -- equation vs equation*, align (NOT
    eqnarray), `\\ref{eq:foo}~(5)` with tilde, where-clause rule.
  * `reviewer_patterns` -- **12 common EM-paper reviewer comments**
    from lab 2019-2025 paper-review history + the pre-check tool
    for each (sign inconsistency, B-vs-H, mesh independence,
    CPU/memory cost, contribution clarity, missing recent
    citations, error bars, eq-symbol collision, ...).

  **NEW tool** `paper_writing_em_paper_style(topic)` with 17 topic
  keys + aliases.

### Pass 3: NEW pre-submission gate orchestrator

  **NEW tool** `paper_writing_em_submission_gate(tex, pdf, bib,
  abstract, author_last_names, page_limit, ...)` runs ALL relevant
  checks in one pass and returns a single verdict:

  * `"verdict"`: "pass" | "warn" | "fail"
  * `"checks"`: per-check status (pass / warn / fail / skip)
  * `"advice"`: actionable next step

  Chains: figure_forward_reference + equation_numbering +
  count_underlines (tex) + lint_reference_format +
  check_citation_usage + self_citation_ratio (bib) +
  validate_abstract_length + abstract_background_ratio +
  abstract_weak_expressions (abstract) + validate_pdf_pages +
  detect_page_whitespace_anomalies + check_floats_far_from_reference
  (pdf).  Skips gracefully on missing inputs.

### Tool count

  v0.89.0: 80 tools
  v0.90.0: 80 + 2 = **82 paper_writing_* tools total**

### Verification

  - mcp-server-paper-writing --selftest -- 10 sections OK
  - tests/test_paper_writing.py -- 75 tests pass (65 from Pass 1 +
    10 from Pass 2/3)
  - Full in-process pytest suite -- 162 + 75 = **237 tests pass**

## 0.89.0 — GitHub survey absorption: arXiv LaTeX source + Semantic Scholar graph

Released 2026-05-26.

User directives (2 of them):

  1. "mcp-server-document には、paper-writing は残さない。世界に公開
     する以上。"  (the LAB-private copy at
     mcp-server-document.paper_writing was deleted in the
     accompanying mcp-server-document v3.3.0 release.)
  2. "世の中の github を探し電磁場にとってよいものがあれば吸収して
     進化させる。"

Surveyed 2026-05-26 (3 web searches, 20+ projects):

  * takashiishida/arxiv-latex-mcp  -- arXiv .tex source fetch
  * blazickjp/arxiv-mcp-server     -- Semantic Scholar ref/citation graph
  * openags/paper-search-mcp        -- multi-source arXiv search
  * Tejas242/arxiv-mcp             -- arXiv ID-based search
  * ScienceAIHub/PaperMCP           -- 32 unified tools (over-scope)
  * PaperDebugger/PaperDebugger    -- Research->Critique->Revision (already
                                       covered by run_full_workflow)
  * citecheck (arXiv:2603.17339)    -- workspace-aware bib repair (already
                                       partially in lint_reference_format)

Absorbed the 2 most useful patterns for EM papers (math-heavy,
single-author-LaTeX-source matters):

### NEW module `_arxiv_source.py` (~600 LOC) + 7 NEW MCP tools

  * `paper_writing_arxiv_fetch_latex_source(arxiv_id)` -- download the
    e-print tarball + extract all .tex files.  Identifies main .tex
    by \\documentclass.  Returns up to N files capped to M chars
    each (defaults 5 / 200k).  Accepts any common arXiv ID form
    (DOI-style, arXiv:, URL, old physics/0501123, with or without
    vN suffix).  Inspired by takashiishida/arxiv-latex-mcp.

  * `paper_writing_arxiv_extract_equations(latex_source)` -- regex-
    scan 6 displayed-math environments (equation / align / eqnarray
    / gather / $$..$$ / \\[...\\]) and return all eqs with body +
    char offset.  Truncates per-eq at 600 chars by default.

  * `paper_writing_arxiv_search(query, categories="", sort=)` -- arXiv
    Atom XML API search.  EM-tuned default `categories` filter
    (eess.SP, physics.class-ph, physics.app-ph, physics.comp-ph,
    physics.plasm-ph, physics.space-ph, cs.CE, math.NA).  Pass
    `categories="all"` to disable the filter.  Inspired by
    openags/paper-search-mcp.

  * `paper_writing_semantic_scholar_lookup(paper_id, fields=)` --
    Semantic Scholar API metadata lookup.  Accepts DOI, arXiv ID,
    Corpus ID, PubMed ID, ACL ID, or S2 URL.

  * `paper_writing_semantic_scholar_references(paper_id, limit=50)`
    -- list of papers CITED BY the given paper.

  * `paper_writing_semantic_scholar_citations(paper_id, limit=50)`
    -- list of papers CITING the given paper.  Useful for finding
    follow-up work to a foundational paper.

  * `paper_writing_external_sources_recipe()` -- knowledge module
    with the 7-project survey writeup + decision tree (which tool
    when) + EM-default category list + credits.

### Lazy import

`requests` is imported only when an arXiv/S2 tool is called.  Users
without `requests` see a clear ImportError pointing at pip install.

### Selftest

mcp-server-paper-writing --selftest now exercises:
  - external_sources_recipe size + key strings
  - _normalize_arxiv_id() on 3 input forms (offline)
  - arxiv_extract_equations() on synthetic 2-equation LaTeX (offline)

### Tool count

  67 (v0.88.0) + 5 (layout-visual) + 1 (tex-figure-placement) +
  7 (arxiv-source) = **80 paper_writing_* tools total**.

### Companion release (mcp-server-document)

mcp-server-document v3.3.0 (LAB-private) deleted the now-redundant
`mcp_server_document/paper_writing/` directory.  The public PyPI
`radia-mcp.paper_writing` is the canonical implementation going
forward.

## 0.88.0 — paper_writing migration + image-based layout + LaTeX figure-placement

Released 2026-05-26.

User directive: "S:\\mcp-server の paper writing は、radia に完全移植
しよう。paper-writing では、ちゃんと画像でレイアウトを確認するスキル
もつけよう。 どうやって確認するかは、web 検索でノウハウを吸収しよう。
TeX の figure の配置スキルを paper writing に反映してほしい"

3-part release:

### Part 1: complete migration

- **NEW subpackage** `radia_mcp.paper_writing` (~5000 LOC, 67 tools)
  promoted from `s:/mcp-server/src/mcp_server_document/paper_writing/`:
    * tools.py (33 paper_writing_* tools + 20 Plan-B helpers from
      plans/T1-T20)  -- byte-identical migration
    * cross_lint.py + _ja_lint.py (8 JA-lint wrappers; original
      cross-package import to ..grant_writing.tools replaced by
      inline-copy in _ja_lint.py so radia-mcp is self-contained;
      30 KB extracted via AST from grant_writing/tools.py)
    * _shared/{hedges,language,sentence}.py (3 small shared
      modules inline-copied so the cross-package _shared dep is
      severed)
    * paper_download.py (IEEE Xplore / ScienceDirect / Emerald
      cookie-seeded PDF download + paper_writing_fetch_and_cite
      orchestration -- byte-identical)
    * skill.md (712-line writing-skill knowledge)
    * plans/T1-T20 (20 Plan-B composite-score helpers --
      byte-identical)

- **NEW entry point** `mcp-server-paper-writing`.
- **NEW catalog entry** `paper-writing` (tag=`meta`), cross-linked
  bidirectionally with literature-index, graph, chart2d, md2html.

### Part 2: image-based PDF layout verification (NEW skill)

User directive: "画像でレイアウトを確認するスキル ... web 検索で
ノウハウを吸収しよう".  Researched 2026-05-26:
  * pymupdf get_pixmap(dpi=N) for per-page PNG render
  * hgustafsson/skillnad (SyncTeX visual-diff inspiration)
  * Overleaf "Understanding underfull and overfull box warnings"
    -- the text-only checks are insufficient

- **NEW module** `_pdf_layout_visual.py` (~400 LOC) +
  **5 new MCP tools**:
    * `paper_writing_render_pages_to_png(pdf, out_dir, dpi,
      page_range)` -- per-page PNG dump for vision-agent
      inspection
    * `paper_writing_detect_page_whitespace_anomalies(pdf,
      whitespace_threshold=0.75, dpi=100)` -- algorithmically
      flag mostly-white pages (sign of bad [!h] / [H] float
      placement forcing figure to own page)
    * `paper_writing_layout_thumbnail_strip(pdf, out_path,
      dpi=80, cols=4)` -- composite PNG of ALL pages as tiles
      for one-glance visual scan
    * `paper_writing_check_floats_far_from_reference(tex, pdf,
      max_pages_apart=1)` -- heuristic detection of figures
      whose \\ref{} appears > 1 page from the float
    * `paper_writing_layout_visual_recipe()` -- 4-step lab
      workflow combining algorithmic + vision-agent inspection

- Lazy imports of pymupdf + PIL with clear ImportError hints
  (radia-mcp ships without these by default; users install via
  pip install pymupdf Pillow).

### Part 3: LaTeX figure-placement knowledge (NEW skill)

User directive: "TeX の figure の配置スキルを paper writing に
反映してほしい".  Web-researched 2026-05-26:
  * LaTeX/Floats, Figures and Captions (Wikibooks canonical ref)
  * Hyndman blog "Controlling figure and table placement"
  * Overleaf "!h float specifier changed to !ht" (silent rewrite)
  * IEEE / IEEJ / IGTE journal LaTeX templates

- **NEW module** `_tex_figure_placement.py` (~600 LOC, 20 topic
  keys + aliases) +
  **1 new MCP tool** `paper_writing_tex_figure_placement(topic)`:
    * OVERVIEW: TL;DR + lab reviewer-pattern
    * FLOAT_SPECIFIERS: h/t/b/p/H + ! semantics; why [H] is
      usually wrong; why [h!] silently becomes [!ht]
    * PLACEINS_FLOATBARRIER: placeins.sty + [section] option as
      lab default
    * WIDTH_FRACTIONS: \\columnwidth vs \\textwidth, figure vs
      figure*, subcaption/subfigure widths
    * ANTIPATTERNS: 8 common mistakes (h! trap, center vs
      centering, blank line in figure body, missing \\label
      after \\caption, mixing absolute and relative widths, ...)
    * JOURNAL_PROFILES: IEEEtran / jiee / igte_digest / scrbook
      defaults
    * CROSS_REFERENCES: which paper_writing_check_* tool maps
      to which placement defect + write-compile-check-fix loop

### Tests + release verification

- mcp-server-paper-writing --selftest: clean (67 tools + 6 new)
- meta-health: 9/9 (bidirectional related-edge invariant
  maintained: 4 reverse edges added to literature-index, graph,
  chart2d, md2html, all pointing back to paper-writing)
- LAB editable install at `S:\\Radia\\01_GitHub\\packages\\radia-mcp`
  refreshed to 0.88.0
- pyproject.toml + __init__.py bumped to 0.88.0

Catalog: **39 → 40 servers**.

NOTE: mcp-server-document.paper_writing source remains for now
(unlike md2html / graph / mathematica which were deleted during
their migrations).  Per migration-policy convention, the LAB-private
copy can be retired in a follow-up commit once the user confirms the
public version is the canonical one.

## 0.87.0 — Sommerfeld layered-medium Green's function knowledge

Released 2026-05-26.

User directive: "W:\\03_文献・論文\\00_電磁界解析\\11_BEM_モーメント法\\
10_sommerfeld_layered ここを学ばせてください" — absorb 3 Sommerfeld
integral PDFs into the BEM subpackage.

Sources absorbed:

  1. **Chew Lecture 35** ("Sommerfeld Integral, Weyl Identity", from
     "Lectures on Electromagnetic Field Theory", Purdue ECE 604/618
     post-2018 revisions, 22 pages).  Pedagogical derivation: scalar
     Helmholtz -> 3D Fourier -> contour deformation -> Weyl identity
     (35.1.11) -> polar reduction -> Sommerfeld identity (35.1.14) ->
     layered-medium dipoles (VED 35.2.6, HED 35.2.12-13) + branch
     cuts / Riemann sheets / Sommerfeld Integration Path.

  2. **Koh & Yook 2006**, IEEE Trans. Antennas Propagat. 54(9),
     2568-2576, DOI 10.1109/TAP.2006.880747.  EXACT closed-form
     Sommerfeld integral for a dipole over an impedance half-plane
     via incomplete Weber integrals / incomplete Lipschitz-Hankel
     integrals / Lommel functions.  Three-layer construction:
     exponential-integral series (eq. 8), Lommel series (eq. 10),
     exact closed form (eq. 17).  Valid for ANY complex eta_s and
     ANY (rho, z+z'), no patchwork like Banos / Wait / Norton.

  3. **Spectral Expansions of Source Fields** (same Chew lecture 35
     from a different revision, 16 pages).  Same derivation +
     stationary-phase asymptotic.

What shipped:

- **NEW knowledge module**
  `radia_mcp/bem/sommerfeld_layered_knowledge.py` (~600 LOC, 22
  topic keys + aliases).  9 main sections:
    * OVERVIEW (when to use, identity summary, exp(-i omega t) vs
      exp(+j omega t) conversion warning)
    * WEYL_DERIVATION (3D Fourier + contour-deformation construction)
    * SOMMERFELD_DERIVATION (polar -> Bessel J_0 -> Hankel form)
    * LAYERED_MEDIUM_DIPOLES (VED TM-only, HED TM+TE,
      tilde R^{TM/TE} recursion as black box, sanity-check recipe)
    * BRANCH_CUTS_AND_SIP (log branch at origin, algebraic branches
      at +/- k_0, surface-wave poles, lossless regularisation)
    * KOH_YOOK_2006_CLOSED_FORM (the headline result + all 3 layers
      + special-function inventory + validation Figs. 2-8)
    * NUMERICAL_EVALUATION (gap-filler since the absorbed PDFs do
      NOT cover this: DCIM Chow 1991 / Aksun two-level 1996 / MPIE
      Michalski-Mosig 1997 / GPOF Hua-Sarkar 1989 / VECTFIT
      Gustavsen-Semlyen 1999 + open-source impl pointers)
    * LAB_USAGE_NOTES (where lab DOES and DOES NOT need Sommerfeld;
      cross-references to bem_low_freq, bem_mmm_msc, bem_h_matrix,
      analytical_formulas, matrix_solvers, ndt, litz_transmission)
    * REFERENCES (the absorbed bundle + Chew 1995 + Kong 1972 +
      Brekhovskikh + Felsen-Marcuvitz + Banos + acceleration refs)

- **NEW `@mcp.tool() bem_sommerfeld_layered(topic)`** in
  `mcp-server-bem`.  22 topic keys with aliases:
    overview / weyl / weyl_derivation / sommerfeld / layered /
    ved_hed / branch_cuts / sip / riemann / koh_yook /
    impedance_plane / closed_form / numerical / dcim / gpof /
    acceleration / lab_usage / lab_notes / references / bibliography.

- **NEW prompt entry** `pick_a_bem_method("dipole_over_layered_medium")`
  with 4-step guidance (overview -> Koh-Yook closed form for single
  interface -> DCIM for multi-layer -> SIP contour theory).

- **bem `--selftest` extended** to exercise the new tool with
  size + unknown-topic assertions.

POLICY: Radia core (HDiv-VIM) is Laplace-kernel only.  This module
is the THEORY POINTER for users who need a full layered-medium
Green's function — production answer is usually either Koh-Yook
2006 closed form (single interface) or ngsolve.bem + DCIM
(multi-layer).  The knowledge module documents the path; it does
not ship runnable code, because lab production rarely needs it.

Catalog: **39 servers** unchanged (new tool inside existing
mcp-server-bem).

## 0.86.0 — FEM flux-line tracing + 3 lab-canonical post-processing plots

Released 2026-05-26.

User directive: "S:\\FEMM\\2020_01_06_磁束線 ここもradia-mcpのgraphに
反映" — reflect the 2020-01-06 lab MATLAB/FEMM workflow for tracing
magnetic flux lines and visualising the field along them into the
graph subpackage.

Source absorbed:

- `S:\\FEMM\\2020_01_06_磁束線\\main.m`         (MATLAB driver)
- `S:\\FEMM\\2020_01_06_磁束線\\plot_trajectory.m`
- `S:\\FEMM\\2020_01_06_磁束線\\velocity.m`     (rhs for ode45)
- `S:\\FEMM\\2020_01_06_磁束線\\磁束線の方程式.docx` (dx/ds = B equations)
- `S:\\FEMM\\2020_01_06_磁束線\\陽的_陰的シンプレクティック.jpg`
                                                 (symplectic Euler note)
- 4 reference PNGs: flux_line, B_in_elements, s-B, s-Az

What shipped:

- **NEW module `radia_mcp/graph/_fem_postprocess.py`** (~470 LOC, 4
  public functions + knowledge text):

    * `trace_flux_line(B_func, x0, s_span, method='RK45', max_step,
      n_eval, ...)` — scipy.integrate.solve_ivp wrapper that
      integrates `dx/ds = B_x(x, y)`, `dy/ds = B_y(x, y)` from a seed
      point.  Replaces MATLAB `ode45(@velocity, [0,10000], [5,0])`.
      Returns `(s, xy)` where `xy.shape == (N, 2)`.  Lazy scipy
      import so module loads without scipy.

    * `plot_flux_line_trajectory(ax, traj_x, traj_y, mesh_nodes,
      mesh_elements, source_xy, ...)` — reproduces `flux_line.png`:
      thin blue FEM-mesh wireframe + black trajectory + red filled
      source marker, equal-aspect axes.

    * `plot_field_probe_line(ax, y, Bx, By, element_change_indices,
      ...)` — reproduces `B_in_elements.png`: Bx (vermillion) +
      By (blue) along a probe line with **vertical dotted lines
      at element-boundary crossings**.  Makes the piecewise-linear
      FEM honest in the figure.

    * `plot_field_vs_arclength(ax, s, Bx, By, Az, ...)` — reproduces
      `s-B.png` (Bx, By vs s in m/T) and `s-Az.png` (A_z drift
      diagnostic).  Automatic twin-y when Az is provided alongside
      Bx/By.

    * `get_flux_line_knowledge()` — full recipe text (equations,
      symplectic-integrator note, units gotchas, complete example
      pipeline using `paper_figure_8cm(panels=2)`).

- **NEW `@mcp.tool() flux_line_recipe(query='knowledge'|'recipe'|'api')`**
  in `mcp-server-graph`:
    * `'knowledge'` — full recipe + equations + symplectic note
    * `'recipe'` — ready-to-paste Python (uses paper_figure_8cm panels=2)
    * `'api'` — function signatures

- Lab colour convention upgraded from MATLAB pure red/blue to the
  CVD-safe Okabe-Ito pair:
    * Bx = `#D55E00` vermillion
    * By = `#0072B2` blue
  Distinguishable in greyscale + passes the lab colorblind-safe gate.

- Knowledge captures the SYMPLECTIC INTEGRATOR note (陽的/陰的) for
  Hamiltonian charged-particle tracking, distinguishing the flux-line
  ODE (NOT Hamiltonian, RK45 is fine) from particle tracking
  (Hamiltonian, needs symplectic).

- **NEW tests** (16 in `tests/test_graph_paper_figure.py`):
    * uniform B → straight trajectory; circular B → orbit closes to <1e-4
    * scipy absent → clean ImportError with pip-install hint
    * trajectory draws mesh + traj + source artists; skips mesh when None
    * sets equal aspect (circular orbit looks circular)
    * field-probe draws bx_line + by_line + N boundary verticals
    * field-probe shape validation
    * arclength: B-only / Az-only / B+Az creates twin y-axis
    * arclength shape validation
    * recipe tool: knowledge / recipe / api / unknown query

**Suite total: 187 → 203 pytest pass.**

Catalog: **39 servers** unchanged (new tool lives inside the existing
mcp-server-graph).

## 0.85.0 — paper_figure_8cm: lab 8-cm-column anchor (1 vs 2 panels)

Released 2026-05-26.

User directive: "8cmに2横に並べて2枚置く場合には極力axesオブジェクトを
大きくして、情報量を増やす。8cmで1枚の場合にはaxesは5cm程度にして
場所を取りすぎないようにする" — the 8 cm column has TWO opposite
sizing regimes depending on whether one or two panels share it.

What shipped:

- **NEW `paper_figure_8cm(panels=1|2, profile, target_axes_cm=5.0)`**
  in `radia_mcp.graph._paper_figure`, exported from
  `radia_mcp.graph`.  Single helper covers both regimes:

  | panels | Figure width | Axes box (each) | Strategy |
  |--------|------------|-----------------|----------|
  | 1 | ~6.5 cm (rel_width derived from target) | **~5.0 × 3.3 cm** | Leave whitespace on the column; figure does NOT dominate |
  | 2 | full 8.89 cm | **~3.6 × 5.3 cm** | `sharey=True` drops inner ylabel + ticks; `wspace=0.10` (~3 mm inner gap); each panel claims every mm |

  Both keep the LAB FONT RULE (10 pt absolute @ 8 cm; text size
  identical between regimes — only axes box scales).

- **NEW `paper_figure_8cm_recipe(panels, profile, target_axes_cm,
  panel_labels)`** MCP tool in `radia_mcp.graph.server`.  Emits
  self-contained Python with the derived geometry comment block
  (figure width, margins, per-panel axes dimensions).

- Default profile is `ieee_single_column` (88.9 mm); also accepts
  `ieej_single_column` (88 mm) and `igte_digest_single` (82 mm).
  Other profiles raise with the allowed list.

- `panels` must be 1 or 2 — for 1x3 / 2x2 / etc., users keep using
  `paper_figure()` directly with its per-(R,C) margin deltas.

- `target_axes_cm` for panels=1 derives `rel_width` from the user's
  target axes width; validates the derived `rel_width` is in
  [0.50, 1.00] and raises with an actionable range hint otherwise.

- **NEW pytest coverage** (13 tests in `tests/test_graph_paper_figure.py`):
    * panels=1 axes width ~5 cm (±0.3); figure < 7.5 cm (doesn't dominate)
    * panels=2 figure = full 8.89 cm; each panel axes 3.3-3.9 cm
    * panels=2 default `sharey=True`; explicit `sharey=False` works
    * `target_axes_cm` scales axes width linearly (±0.4 cm at 4/5/6 cm)
    * font is 10 pt in BOTH regimes (lab font rule unbroken)
    * panels=1 figure < panels=2 × 0.85 (regression guard against
      rel_width drift to 1.0)
    * panels=2 total axes area > panels=1 × 1.5 (max info density check)
    * `panels` invalid (0, 3, 4, -1) → raises with policy hint
    * profile not 8 cm column → raises with allowed list
    * `target_axes_cm` outside derivable range → raises with hint
    * recipe tool returns runnable Python with correct
      `min_axes_fraction` (0.55 for panels=1, 0.72 for panels=2)
    * recipe tool rejects bad panels / profile

**Suite total: 174 → 187 pytest pass.**

Catalog: **39 servers** unchanged (new tool lives inside the existing
mcp-server-graph).

Quick-start:
```python
from radia_mcp.graph import paper_figure_8cm, emit_paper_figure

# Single figure on a column — don't dominate
fig, axes = paper_figure_8cm(panels=1)            # axes ~ 5 × 3.3 cm
axes[0, 0].plot(x, y)
emit_paper_figure(fig, 'single', 'ieee_single_column',
                  min_axes_fraction=0.55)

# Two figures sharing a column — max info density
fig, axes = paper_figure_8cm(panels=2)            # each panel ~ 3.6 × 5.3 cm
for ax in axes.flat:
    ax.plot(...)
axes[0, 0].set_ylabel('|Z| (Ω)')                   # only LEFT panel
emit_paper_figure(fig, 'pair', 'ieee_single_column')
```

## 0.84.0 — NGSolve mirror of COMSOL topology-optimization knowledge

Released 2026-05-26.

User directive: "COMSOLのmcp-serverの知見は、radia-mcpにはNGSolveにも
反映させた上で継続学習。S:\\NGSolve\\03_TolologyOptimization にあるPDF
も追加学習" — i.e. the topology-optimization knowledge that lives in
COMSOL MCP (`docs/TOPOLOGY_OPTIMIZATION.md` + RAG prompt) must also
be reflected in the radia-mcp NGSolve subpackage as the
implementation-flavoured mirror, plus the 3 Gangl PDFs from
`S:\NGSolve\03_TolologyOptimization` should be absorbed.

What shipped:

- **NEW knowledge module**
  `radia_mcp/radia_ngsolve/knowledge/topology_optimization.py`
  (~1900 LOC, 13 super-topics with primary keys + aliases). Mirrors
  the COMSOL MCP long-form tutorial in an NGSolve-implementation
  flavour: weak forms, HCurl + eps*mass regulariser pattern,
  Newton-Raphson tangent reluctivity assembly, periodic / Kelvin BC
  interplay, Amstutz-Andra fixed-point recipe.

- **NEW `@mcp.tool() topology_optimization(topic)`** wired into
  `radia_mcp.radia_ngsolve.server`.  Canonical topics + aliases:
  `overview/framework`, `ch1_2/foundations`, `ch3/evolution_equation`,
  `ch4_5/heat_elasticity`, `ch6_8/fluid`, `ch9/lbm`,
  `appendix_a/objectives`, `appendix_c/helmholtz_filter`,
  `appendix_d/kkt`, `gangl_part1/gangl_sensitivity/averaged_adjoint`,
  `gangl_part2/ngsolve_implementation`, `gangl_motor/ipm_motor`,
  `ngsolve_recipes/recipes`.

- Source bundles absorbed (continuous-learning batch, 2026-05):
    * **Nishiwaki-Kondoh-Yachi "Foundations of Topology Optimization"**
      (Corona Publishing 2024, 267 pages, lab's canonical
      Japanese-language reference): ch.1 history; ch.2 4-method-family
      + OC/SCP; ch.3 unified reaction-diffusion engine;
      ch.4 heat conduction; ch.5 elasticity; ch.6 Stokes; ch.7 laminar
      NS; ch.8 conjugate HT; ch.9 LBM; App A objective cookbook;
      App C Helmholtz filter; App D KKT dual evolution.
    * **Gangl-Sturm 3-paper bundle** (TU Graz + TU Wien, 2015-2019,
      `S:\NGSolve\03_TolologyOptimization` 3 PDFs):
        - Part I = sensitivity analysis for nonlinear curl-curl
        - Part II = NGSolve implementation pattern
        - IPM motor case study (27% cogging-torque-surrogate
          reduction; **corrected** from the previously-cited 87% --
          captured in `GANGL_IPM_MOTOR_CASE` knowledge).

- Pairs with the existing standalone `radia_mcp.topology_optimization`
  server (Gangl theory: shape_optimization / topology_derivative /
  applications) -- the new NGSolve tool is the implementation-side
  reflection, not a replacement.

- `radia_ngsolve._selftest()` extended to exercise the new tool:
  asserts `topology_optimization("all")` > 10k chars,
  `topology_optimization("overview")` > 500 chars, unknown-topic
  raises with `"Unknown topic"` substring.

Cross-reference:
  - COMSOL MCP: `docs/TOPOLOGY_OPTIMIZATION.md` (long-form tutorial)
  - COMSOL MCP: `src/knowledge/prompts/topology_optimization.md`
    (RAG-indexed prompt for retrieval)

Catalog: **39 servers** unchanged (no new entry point -- the new
tool lives inside an existing server).

## 0.83.0 — catalog cleanup: graph shim + mcmc removal

Released 2026-05-26.

Bookkeeping release — no new MCP servers, just hygiene around the
catalog after the v0.82.0 mass-migration batch.

What shipped:

- `mcmc` removed from catalog + entry-points + pyproject scripts
  (subpackage source dir was deleted upstream by the user).
- `tests/test_meta_health.py::test_meta_related_to_mcmc_includes_optuna`
  renamed to `test_meta_related_to_chart2d_includes_graph` to use a
  stable cross-link pair (both chart2d and graph shipped 2026-05 and
  are not at risk of removal).
- Companion cleanup at `s:/mcp-server/`:
    - `mcp_server_document/graph/` directory **deleted** (959 LOC
      backward-compat shim).  `mcp-server-graph` is canonical at
      `radia_mcp.graph` since v0.77.0; the shim was no longer needed.
    - `mcp_server_document/server.py::_build_mcp()` no longer
      imports/registers `graph` (was 12 subpackages → 11; 326 → 324
      tools).
    - `mcp-server-document` bumped to v3.2.0.

Catalog: **40 → 39 servers** (mcmc removal).

Suite: **174/174 pytest pass**.

## 0.82.0 — md2html migration + mcp-server-document cleanup

Released 2026-05-26.

User decision (after s:\mcp-server enumeration):
  - DISCARD: `mcp_server_document/diagram/` (raster → Excalidraw
             pipeline) + sibling `s:/mcp-server/excalidraw/` Node.js
             tree.  Not LAB-aligned; bulky deps; better third-party
             options exist.
  - COMPLETE MIGRATION: `mathematica/` → `radia_mcp.mathematica`
             (already done in v0.69; this release removes the stale
             source copy from `mcp-server-document`).
  - COMPLETE MIGRATION: `md2html/` → `radia_mcp.md2html` (NEW this
             release).

What shipped:

- **NEW subpackage `radia_mcp.md2html`** (393 LOC, verbatim from
  `s:/mcp-server/src/mcp_server_document/md2html/`):
    converter.py (302 LOC) + tools.py (69 LOC) carried over
    byte-identical; only the registration layer (server.py, 95 LOC) is
    new and follows the standard radia_mcp.<topic>.server pattern
    (--selftest, status_tool wired, FastMCP).
  - Tool: `md2html_convert(md_file, output_file=None, title=None)`
  - Features (preserved): MathJax v3 CDN, $..$/$$..$$/```math```
    blocks, ||x|| → \\Vert x \\Vert, [N] reference linking +
    auto <li id="refN"> under References / 参考文献, base64 image embed,
    UTF-8 + cp932 legacy fallback.
  - Optional dep: `pip install radia-mcp[md2html]` brings `markdown>=3.0`.

- **NEW catalog entry** `md2html` (tags=['meta']).  Catalog grows
  **39 → 40 servers**.  Related: mathematica, graph, literature-index
  (all bidirectionally cross-linked).

- **NEW `mcp-server-md2html`** entry point in pyproject.toml.

- **NEW tests/test_md2html.py** -- 12 golden tests covering:
    basic .md → .html round-trip; MathJax script emission; $$..$$
    block preservation; ||v|| normalisation; [N] linking; Japanese
    参考文献 section; base64 image embed; cp932 fallback; missing
    file rejection; non-.md rejection; lower-level md_to_html dict
    shape; explicit title override.

- **mcp-server-document v3.0.0 → v3.1.0** (LAB-private):
    Removed `diagram`, `mathematica`, `md2html` from
    `server.py::_build_mcp()` (12 subpackages, 326 tools remain).
    Removed `diagram` keyword + extras + `diagram/skill.md`
    package-data + `md2html` extras from pyproject.toml.
    Description string updated to drop "Excalidraw diagram pipeline".

Suite: **175/175 pytest pass** (was 162 after 0.81.0; +12 md2html
+1 bidirectional-link fix for md2html cross-links).

Catalog evolution:
  v0.76.0   37 servers
  v0.77.0   38 (+ graph)
  v0.80.0   38 (graph upgrades)
  v0.81.0   39 (+ chart2d)
  v0.82.0  *40* (+ md2html)  ← this release

Survey verdict (4 candidates from s:\mcp-server, 2 chosen):
  ✅ MOVED: md2html (393 LOC, generic utility, low dep)
  ✅ STALE COPY DELETED: mathematica (already in radia-mcp)
  ⏸  DEFERRED: pdf (1055 LOC PDF toolkit) -- still recommended for
                v0.83+, would absorb the Type-42 verifier we built
                in v0.80.0 into a richer toolkit
  ❌ KEPT IN mcp-server-document (LAB-private, scope mismatch):
     grant_writing, paper_writing, presentation, poster, circuit,
     bibliography, doc_convert, ocr, document_meta, research_project,
     graph (already migrated, OLD source intentionally kept as
     backward-compat shim for now), eqnedt32, eqnedt32_native
  ❌ DELETED outright (per user "捨てよう"):
     mcp_server_document/diagram/ + sibling excalidraw/ Node.js tree

## 0.81.0 — chart2d: 22 paper-quality 2D chart MCP tools

Released 2026-05-26.

After v0.80.0 user requested a dedicated MCP server for 2D chart
rendering, citing the GitHub survey results (StacklokLabs/plotting-mcp,
xlisp/visualization-mcp-server, LindseyyyLi/MCP-Server,
antvis/mcp-server-chart 4.1k★, isaacwasserman/mcp-vegalite-server,
arshlibruh/plotly-mcp-cursor) as prior art.  The architectural
finding: no existing MCP server combines (a) paper-quality styling
with (b) recipe + image dual return, so chart2d fills both gaps.

What shipped:

- **NEW subpackage `radia_mcp.chart2d`** -- companion to
  radia_mcp.graph.  graph = styling / profiles / gates.
  chart2d = data → rendered 2D chart.  Catalog: **38 → 39 servers**.

- **NEW `mcp-server-chart2d`** with **22 chart-type tools**:

    LINE / TIME-SERIES (8):
      chart2d_line, chart2d_loglog, chart2d_semilogx, chart2d_semilogy,
      chart2d_step, chart2d_errorbar, chart2d_fill_between, chart2d_bode

    DISTRIBUTION / STATS (5):
      chart2d_histogram, chart2d_bar, chart2d_box, chart2d_violin,
      chart2d_ecdf

    SCIENTIFIC 2D FIELDS (7):
      chart2d_contour, chart2d_contourf, chart2d_pcolormesh,
      chart2d_quiver, chart2d_streamplot, chart2d_imshow, chart2d_polar

    POINT (2):
      chart2d_scatter, chart2d_phase  (Nyquist / impedance locus,
        aliases nyquist / complex_plane / argand / smith_like)

  + `chart2d_catalog(group)` introspection tool (23 total).

- **DUAL RETURN MODE** -- every chart tool accepts
  `return_mode='recipe'|'image'|'both'`:

    'recipe' (default):  returns Python code string ending in
                          emit_paper_figure().  User pastes into a
                          script and runs locally.  No matplotlib
                          dependency on the server.
    'image':             server-side render via matplotlib + paper_figure
                          profile, returns PNG bytes as MCP Image
                          content type (LindseyyyLi/MCP-Server +
                          isaacwasserman/mcp-vegalite-server pattern)
                          -- inline preview in Claude Desktop /
                          claude.ai.
    'both':              returns {'recipe': str, 'image': Image}.

- **PROFILE INHERITANCE** -- every render path calls
  radia_mcp.graph.paper_figure() so the chart inherits 10 pt @ 8 cm
  font, Okabe-Ito CVD-safe palette, Type-42 PDF embed, and IEEE/IEEJ
  margins automatically.

- **SPECIAL-CASE RENDERERS**:
    `_render_bode`  -- 2x1 panel pair (magnitude/phase) with sharex,
                       top panel's x-tick labels hidden, lab-style
                       labels (|H| (dB), arg(H) (deg)) by default.
    `_render_polar` -- handles projection='polar' axes (which
                       paper_figure() does not directly expose).

- **`tests/test_chart2d.py`** -- 37 tests, ~3 s:
    pin: catalog has exactly 22 types in 4 groups; aliases resolve;
         every type has a non-empty paper-figure-gated recipe (5×22
         parameterised tests); 6 representative types render to
         valid PNG; Bode produces (2,1) layout; phase plot has
         equal aspect; polar has projection='polar'; first-line color
         is OKABE_ITO[0] (palette inheritance lock); unknown type
         raises with hint; chart2d_catalog tool shape.

Suite: **162/162 pytest pass** (was 124 after 0.80.0; +38 new).

Survey patterns DEFERRED (with reasons):
  - Vega-Lite spec return (isaacwasserman): consider for chart2d v2
    if a user requests inline-renderable specs that don't require
    matplotlib server-side; current MCP Image path is good enough.
  - 26-chart antvis catalog: their counts include flow / sankey /
    treemap which are out-of-scope for 2D scientific charting; our
    22 covers every chart we've used in lab papers since 2019.
  - generate_sample_data tool (arshlibruh): the recipe-mode default
    already gives the AI a working sample script with synthetic data
    -- the inline `data = np.random.normal(...)` lines serve the
    same purpose.

Catalog evolution:
  v0.76.0   37 servers
  v0.77.0   38 (+ graph)
  v0.80.0   38 (graph upgrades)
  v0.81.0  *39* (+ chart2d)  ← this release

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

## 0.76.0 — optuna: 5 advanced lab BBO recipes

Released 2026-05-25.

`radia_mcp.optuna` stays as the canonical optimization MCP for the
lab's black-box FEM-as-objective EM design problems.

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
- 7 ML/optimization: `bayesian_opt`, `evolutionary`, `gnn`,
  `data_assimilation`, `optuna` (Sano-Akiba-Imamura textbook), `pinn`,
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
  Cauer MEC, Janet 2004-2005 RNA calibrated leakage mixed method).
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
  finds ≥4, related('optuna') includes 'bayesian-opt'.
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
    (1/2, 1/4 only) vs "EM panel FEM/HDiv-VIM" C-yoke (1/1, 1/2, 1/4,
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
