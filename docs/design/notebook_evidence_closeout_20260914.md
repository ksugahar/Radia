# Notebook evidence closeout

Base: `b79c281cd`, 2026-09-14. This investigation separates reproducible
mathematics, repository-generated figures, unsupported historical comparisons
and live-client deployment. It does not rename uncertainty as acceptance.

## Cylinder coefficients: constructive derivation

Let `R(z) = I1(z)/I0(z)`. The modified-Bessel equation
([NIST DLMF 10.25.1](https://dlmf.nist.gov/10.25.E1)) and derivative identity
([10.29.3](https://dlmf.nist.gov/10.29.E3)) give

```math
R'(z) = 1 - R(z)/z - R(z)^2.
```

Substitute `R(z) ~ sum(n>=0) b_n z^(-n)`, `b_0=1`. Coefficient matching gives

```math
2b_n=(n-2)b_{n-1}-\sum_{i=1}^{n-1}b_i b_{n-i}.
```

Thus `b_1..b_6 = (-1/2, -1/8, -1/8, -25/128, -13/32, -1073/1024)`, exactly the
tuple in `src/radia/maglev/mixed_galerkin/references.py`. The existing fast test
module now derives it with rational arithmetic and separately compares the
production continuation against scaled Bessel functions across the crossover
and deep-skin points on the positive-real and positive/negative frequency rays.
This closes the coefficient origin without claiming an unlocated Senior passage.
It does not prove that a finite enriched Galerkin basis automatically recovers
all these coefficients, or revalidate the saved wall-band accuracy table.

For the implemented spherical reference, the normalized surface factor is
`coth(z)-1/z`. In the right half-plane it is `1-1/z` plus exponentially small
terms. The algebraic inverse-power expansion terminates; this does **not** mean
that every further finite-frequency basis function adds no value. That stronger
claim in the historical script has no such proof.

## PEEC figure lineage

The executable source is the `fig1_peec_matrix_verification` through
`fig4_wpt_coil_characteristics` definitions in the saved code cell of
`docs/peec_integration/peec_showcase.ipynb`. No external measurement file is
read by these four functions.

| Figure | Actual source | Acceptance limit |
| --- | --- | --- |
| 1, matrix properties | `PEECBuilder.create_loop`, matrix eigenvalues, symmetry norms | Internal numerical diagnostic, not external measurement. |
| 2, reduction | `LanczosReducer.lanczos_symmetric` and projected R/L matrices | Not demonstrated to be the complete PRIMA algorithm; catches singular solves by substituting a diagonal entry. Plot labels do not establish algorithm equivalence. |
| 3, adaptive error/speedup | Assumed `(a/r)^3` curve, assumed near-field fraction and cost ratio 10 | Conceptual model, **not measured error or benchmark**. |
| 4, WPT | PEEC loop matrices, `sum(R + j*omega*L)` | Simulated model only; the marked 85 kHz sweep point uses the nearest grid sample, not necessarily exactly 85 kHz. |

The original code and outputs remain intact; prominent figure-section text
qualifies their use. No source publication or measurement record is inferred
from a filename. Replacing these figures with benchmark/measurement evidence
requires actual run/data artifacts, not a new bibliography citation.

The `analytical_square_loop` function in
`validation_test/peec_integration/ngsbem_peec_demo/verify_loop_peec_vs_ngsbem.py`
was evaluated in isolation (only that AST function, no solver imports): its
10 mm side / 1 mm width inputs give **29.760609107860034 nH**, not the 24 nH
printed reference. This proves a mismatch between the historical constant and
the displayed formula, not the physical validity of either. Original FastHenry
run/model evidence and a Grover equation locator remain required.

## Remaining-43 screening: source questions are now explicit

The Markdown source of all 43 undeclared notebooks was rescreened for named
methods, references, papers, data providers and benchmark attributions. This
is source-oriented screening, not a claim of full scientific peer review.
Do not add dummy bibliographies to API examples or internal regression galleries.
The following substantive source families need reconciliation before the set
can be called citation-complete:

- Clebsch/hodograph and particle-orbit notebooks: SCOFF/Enge and classical
  fringe/edge-focusing claims; existing internal derivations are not publication
  identities. The cut-selection notebook also names a Takahashi workpiece.
- Hysteresis and NGSolve integration: Egger, measured B-H and Simkin benchmark
  source identities and original measurement lineage.
- Maglev: TEAM 28 benchmark specification. Induction heating: the named IGTE
  digest versus internal results. Neither a benchmark name nor a conference year
  identifies a source by itself.
- Dowell and stream-function companions: the linked theory Markdown owns
  literature today; reconcile the actual Dowell/Kuijpers claims with those
  sources rather than inventing notebook-local bibliographies.
- URN: the additional concrete data-provenance findings below take priority
  over adding generic NASA/TDK citations.

### URN data source correction

`extract_nasa_eis.py:create_representative_eis_csv` generates model data but
formerly labelled it `REAL MEASUREMENT DATA` and overwrote `nasa_18650_eis.csv`.
It now writes `synthetic_18650_eis.csv`, explicitly labelled synthetic; its
documentation writes `SYNTHETIC_DATA.md`, never the measured-data README.
The fresh/aged synthetic generators in `download_nasa_data.py` likewise have
synthetic names and labels. Fast isolated function tests protect measured-file
sentinels and verify the actual generated headers and row counts.

The existing tracked `nasa_18650_eis.csv` claims B0005 cycle 40 extraction and
contains 47 valid points, so the mere existence of a synthetic generator is
**not proof that this saved curve is synthetic**. However,
`extract_real_eis.py` assumes a log-spaced 0.1 Hz--5 kHz frequency vector instead
of reading instrument frequency metadata. Original MAT identity, exact sample
ordering, frequency mapping and CSV hashes must be reconciled before claiming
independent measured-spectrum validation. The TDK generator converts permeability
to impedance with an explicitly assumed ten-turn geometry: the resulting Z is
derived, not directly measured impedance. The notebook now states these limits.
No existing CSV numeric value, fit output or figure was rewritten or retrained.

## Live publication server

The current task's actual `radia-publication` connection was queried through
`capability_pack_status`, `bibliography_status`, `bibliography_canonical_path`
and `bibliography_get_entries`. It reports editable radia-mcp 1.4.53 from
`release-quad/radia-4.95.91-runtime`, Python 3.12.10 and MCP SDK 1.27.0.
Its bibliography has 466 entries; `ieej2007integrals9`, `senior1962note` and
`grover2004inductance` are absent. Therefore main integration is **not live
deployment**. A server-module unchanged hash does not attest to bibliography
currency. No install, source repoint, reload or process termination was performed
by this initial read-only check.

## Primary NASA archive verification (2026-09-14 follow-up)

The [NASA PCoE repository](https://www.nasa.gov/intelligent-systems-division/discovery-and-systems-health/pcoe/pcoe-data-set-repository/)
links the [battery archive](https://phm-datasets.s3.amazonaws.com/NASA/5.+Battery+Data+Set.zip).
The original-source availability question is now closed; a user copy is not
needed. The downloaded archive was inspected in memory without extracting into
the repository. Reproducible identities (SHA-256):

- Outer archive: `82302a7db4fc1b34e0b6676326610438d43b816bdf11a69d1d012a464ef2f92e`.
- Nested `1. BatteryAgingARC-FY08Q4.zip`: `f3a4fcd7c3434e06906f4222778445ce6ad42f49c5c78a785bdd676b7c223a2e`.
- `B0005.mat`: `0eae4585baf3f200c09fe24c5ab884f1889679fc75206ca1aa19da704104f0b0`.

The first impedance operation is zero-based cycle 40. Its 48 complex
`Battery_impedance` samples are finite. The tracked CSV's 47 complex samples
match `Battery_impedance[1:]` **exactly** (maximum absolute difference zero).
The excluded first sample has negative real impedance. This establishes the
saved impedance lineage, not the validity of that filtering choice.

The MAT fields are `Sense_current`, `Battery_current`, `Current_ratio`,
`Battery_impedance`, `Rectified_Impedance`, `Re`, and `Rct`; no frequency field
is present. The archive README provides sweep endpoints but not sample ordering
or spacing. The saved frequency column matches the generated 47-point
`logspace(-1, log10(5000), 47)` exactly. It is therefore a model-assigned axis,
not a verified instrument axis. Frequency-dependent fitting accuracy remains
unqualified until instrument records identify that axis; no retraining can
repair missing independent frequency evidence.

`extract_real_eis.py` now refuses measured export without an explicit,
sample-aligned positive finite frequency vector and source locator. It preserves
the caller's sample order and does not fabricate or silently sort frequencies.
The ZIP member is read in memory, removing a separate overwrite/delete hazard
for pre-existing MAT files. Fast tests cover missing/invalid axes, sample order,
real MAT parsing, unknown cycle selectors, and preservation of caller files.
Historical CSV values and notebook outputs remain unchanged.

## Focused companion citation reconciliation

The stream-function and Dowell executable companions now declare selected
canonical keys and carry generated BibTeX/TeX4ht reference displays. This is two
additional notebook migrations, not a claim that all 43 screened notebooks or
all methods in these two companions are citation-complete.

- The [official COMPUMAG program](https://compumag2023.com/program/www.conftool.pro/compumag2023/index43b4.html?form_session=29&ismobile=false&page=browseSessions&print=export)
  identifies contribution 525 (PC-A1:6, May 25, 2023) as *Comparison of
  Discretization Methods for Continuous Stream-Function Distributions* by
  Kuijpers, Jansen and Lomonova. The parent gains the distinct key
  `kuijpers2023discretization`; `kuijpers2023` continues to identify their
  different journal article. The university repository confirms the conference
  title but lists only Kuijpers, so the author list follows the official program.
  Program identity does not verify the former spatial-collocation attribution.
  That attribution is withdrawn; Path-A is labelled Radia's own residual
  correction. Continuous FE representation alone is no contraction proof, and
  the companion now distinguishes saved empirical convergence from a theorem.
- `dowell1966eddy` is reused for the winding-loss background. Publisher-deposited
  [Crossref metadata](https://api.crossref.org/works/10.1049/piee.1966.0236)
  confirms author, title, volume 113, issue 8, year 1966 and first page 1387.
  This check does not certify the companion's boundary interpretation, scaling,
  continued-fraction/PRIMA composition or saved numerical agreement. Kernel
  agreement is explicitly an internal check, not experimental validation.

No calculation cells, saved numerical outputs or figures were rerun or changed.
The generated-reference contract checks both new notebooks and the earlier
migrations against the parent. ACA/TSVD/DUCAS source reconciliation, the other
named source families above, and missing instrument/external-run evidence remain
separate follow-up work; no dummy citations were added to API-only notebooks.

### Live status correction

The earlier live-server paragraph is a historical pre-update observation.
Following the user's Restart, the original task directly verified the new
`mcp-development-main` source, 349 publication tools and all three previously
absent bibliography keys. That task's publication reconnect is closed. It does
not establish the state of other clients, users or the 100 host.

## Remaining 41 notebooks: full editorial triage

All 41 notebooks without a bibliography declaration at the start of this pass
received an individually written source/evidence scope cell. This is an editorial
assessment, not 41 completed literature reviews or numerical validations.
All original code cells, saved outputs, execution counts and existing notebook
metadata were preserved; only cut selection gains bibliography metadata.

The decisions are **5 API examples, 1 navigation page, 10 internal evidence
records, and 25 with open source or scientific-acceptance work**. An `internal`
classification supports the particular recorded comparison, not a blanket
citation exemption for a future methods paper. No dummy references were added.

Cut selection now cites existing canonical `kotiuga1987cuts` and
`pellikka2013homology` as general background through generated BibTeX/TeX4ht
artifacts. Kotiuga's [author publication list](https://people.bu.edu/prk/Publications.htm)
confirms the 1987 title, journal, volume and pages. These references do not
identify the unnamed Takahashi benchmark or prove the one-handle algorithm.
There are now **11 generated bibliographies among 51 notebooks**. The other
40 absent declarations neither prove citation completeness nor automatically
require a bibliography: follow each notebook's explicit review scope.

### Work retained explicitly

- **CG and Kelvin:** the unconditional Prager-Synge paragraph was rewritten.
  Boundary/source admissibility, constitutive energy norm and theorem
  correspondence remain open. Kelvin retained/lost fractions and growth
  conventions need reconciliation rather than a guessed correction.
- **Fringe/hodograph theory:** resolve Enge/SCOFF and Chaplygin/von Mises
  sources, focusing units and tanh-profile assumptions. Angular-drift reduction
  is not a measured saturation-range extension. Historical output remains saved
  but its overstatements are explicitly superseded by the scope note.
- **Named algorithms and benchmarks:** Love/Stratton-Chu, streamline/LIC,
  geometry processing, hysteresis, TEAM 28, Simkin, optimization and CQ require
  source-to-implementation correspondence. Finding a paper title alone does
  not establish a numerical acceptance certificate.
- **Model/data provenance:** closed-torus BEM behavior is case-specific;
  analytical hysteresis fixtures are not measured samples; slab/cylinder
  comparisons require a geometry argument. URN TDK captions now say derived
  impedance. NASA frequency and external-run gaps above remain open.
- **Optimization:** sampled feasible/nondominated designs are not guaranteed
  optima or complete Pareto fronts. Current norm, magnetic energy and dissipation
  are distinct objectives. Low-rank Tikhonov comparisons must state the retained
  operator and metric/nullspace assumptions.

The next pass should close each precise source/assumption gap and regenerate
references where applicable. It should not remove an open qualification merely
because a bibliography entry exists, or rerun expensive notebooks without a
specific validation question.

Verification: all 41 code/output/metadata preservation checks passed, and the
existing notebook plus canonical-bibliography regression selection passed
**103 tests**. No solver code, parent bibliography, editable installation or
live-client connection was changed in this pass.

## Source-resolution follow-up for the 25 open notebooks

Every one of the 25 received a claim-to-source or claim-to-implementation
resolution, replacing the initial screening note. **19 are now qualified as
bounded demonstrations; 6 retain specific evidence gaps.** Qualified means the
claims were narrowed to the supported demonstration, not that a broad theorem
was proved or a historical benchmark rerun. In particular, unlocated classical
attributions in the hodograph/edge examples are not silently promoted to verified
citations; the displayed local derivation or explicit proxy is the asserted model.

Eighteen additional notebooks now carry canonical BibTeX/TeX4ht artifacts;
cut selection reuses its existing generated references. There are 29 declared
bibliographies in the 51-notebook corpus. Six of these 25 use repository-local
derivations/roadmap claims without invented external references. Code cells,
execution counts, saved outputs and non-bibliography metadata were preserved
in all 25, checked against the start commit.

### Corrections supported by implementation and algebra

- The Kelvin tests check field ratios Z/B1 = 2/3 and B0/B1 = 4/3, hence a
  one-third deficit/excess, not a loss of two thirds. The uniform potential
  grows as inverse radius after inversion; its gradient grows as inverse
  radius squared. Historical diagnostic strings inside saved code/output are
  retained but explicitly superseded in the prose.
- The edge proxy at zero tilt is `-tan(w/(2*rho))/rho`; `-w/(2*rho^2)` is its
  first-order approximation. The finite-z polynomial fringe is a local
  expansion with an order-z-squared curl residual, not an exact vacuum field
  over an arbitrary aperture.
- The IH Bessel helper explicitly passes `geometry='cylinder'`. The class name
  `ESIMFiniteSlabSolver` did not imply a slab/cylinder comparison.
- The CG notebook computes an unweighted H reconstruction difference (with
  radial integration weight), not the constitutive-energy norm or error against
  an exact solution. [Bertrand and Boffi, section 2](https://arxiv.org/html/1907.00440v1)
  supplies the relevant admissibility/equilibrium context, not a blanket bound
  for this saved indicator.
- The regularized formula applies to the retained operator A_k and an SPD free-DOF
  metric. An independent small NumPy comparison with the dense A_k Tikhonov
  equations agreed to 2.33e-16. This is an algebra check, not solver validation.
- A current seminorm is not magnetic stored energy; sampled nondominated trials
  are not a complete Pareto front. Those interpretations were corrected.

### Parent bibliography reconciliation

New entries identify Lubich CQ, Bertrand/Boffi hypercircle reconstruction,
Jobard/Lefer streamline placement, Cabral/Leedom LIC, Lewiner marching cubes,
Trimesh smoothing API and the Optuna framework. Verification used the
[Springer CQ record](https://doi.org/10.1007/BF01398686),
[Springer streamline chapter](https://link.springer.com/chapter/10.1007/978-3-7091-6876-9_5),
[Lewiner's author page](https://thomas.lewiner.org/pub/marching_cubes_jgt.html),
[Trimesh API](https://trimesh.org/trimesh.smoothing.html), author arXiv records,
and publisher-deposited Crossref metadata. References are background only where
the implementation or experimental correspondence has not been demonstrated.

Publisher metadata for [10.1109/TMAG.2025.3621134](https://api.crossref.org/works/10.1109/TMAG.2025.3621134)
identifies Herbert Egger, Felix Engertsberger and Andreas Schafelner, the
forward/inverse title, volume 62(7), article 7300706 (2026).
The misleading conference/journal hybrid record is corrected under its retained
key; the separate [arXiv version](https://arxiv.org/abs/2507.15289) retains a
separate key and correct title. The old `egger2025tmag` title/DOI cannot be
confidently identified: a 404 is not proof of fabrication, so it remains flagged
unverified rather than being silently reassigned to another paper. None of the
new references cite it. Potter's [publisher metadata](https://api.crossref.org/works/10.1109/TMAG.1971.1067251)
restores coauthor R. Schmulian and pages 873-880. Existing generated references
were checked for affected selected-source hashes.

### Six evidence gaps that references cannot close

1. **BEM extractor:** independently specify the current/internal-inductance
   convention and investigate the actual matrix rank/spectrum before assigning
   a general cause to the closed-torus observation.
2. **Shape regeneration:** `topopt_cad.py` passes `nu=-0.53` to Trimesh's
   subtractive dilation pass. Trimesh 4.12.2 independently reproduces this sign
   problem on a one-subdivision icosphere: two-pass volume ratios are 0.57641
   for negative nu versus 0.99187 for positive nu. The code and saved geometry
   are not modified in this documentation pass. A focused implementation fix
   plus regenerated shape evidence is required; merely citing Taubin is not a fix.
3. **Hysteresis gallery:** recover dataset/identification/preprocessing lineage
   for measured-B-H claims. The analytical Potter fixture is not a measured specimen.
4. **Maglev:** match inputs and observations to the current official TEAM 28
   specification, not just its preliminary-description citation.
5. **NGSolve integration gallery:** identify the Simkin-labelled hysteresis image
   and underlying data. An unrelated scalar-potential paper is not its provenance.
6. **URN:** obtain the original sample-to-frequency record and accepted manuscript
   identity; NASA samples and derived TDK impedance remain as qualified earlier.

The 103 notebook/canonical-bibliography tests pass. This pass changes no solver,
installation, running client, or saved numerical result. Parent edits take effect
for a manuscript only after its generated bibliography is refreshed; this is not
a claim of live-client synchronization.
