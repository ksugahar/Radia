"""EM-paper-specific style knowledge for radia_mcp.paper_writing.

Added 2026-05-26 (radia-mcp v0.90.0).  The generic text/lint tools in
``tools.py`` are journal-agnostic; this module captures the
EM-DOMAIN-specific conventions that the lab's reviewers care about
(sign-convention, vector/tensor notation, B vs H usage, SI units,
EM-paper-specific reviewer-comment patterns).

Sources:
  * IEEE Editorial Style Manual (10th edition, IEEE Publication
    Services and Products Board)
  * IEEE Author Center: Manuscript Templates and Tools
  * IEEJ 論文投稿の手引き (D-section / B-section)
  * IGTE Symposium digest template
  * Sugahara lab CLN / IH / motor paper review history
    (reviewer comments collected from accepted IEEE TAP, IEEE TMag,
    IEEE TMTT, IEEJ Trans. D, IEEJ Trans. B papers, 2019-2025)
"""

from __future__ import annotations

import os
import re


SIGN_CONVENTIONS = r"""
# Time / phase sign conventions for EM papers

## The two conventions

EM papers use ONE of two time-harmonic conventions:

  * **Engineering convention**: `exp(+j omega t)`  (IEEE / IEEJ default)
  * **Physics convention**:     `exp(-i omega t)`  (Chew, Brekhovskikh,
                                                    Felsen-Marcuvitz)

These are mathematically equivalent (complex conjugate of each other),
but reviewers HATE mixing them within one paper.

## Translation table (engineering <-> physics)

  imaginary unit       j   <->   -i
  outgoing sph wave    exp(-j k r) / r   <->   exp(+i k r) / r
  Hankel outgoing      H_0^(2)  <->   H_0^(1)
  Fourier kernel       F(omega) = int f(t) exp(-j omega t) dt
                       <->   F(omega) = int f(t) exp(+i omega t) dt
  Drude conductivity   sigma = sigma_0 / (1 + j omega tau)
                       <->   sigma = sigma_0 / (1 - i omega tau)
  Skin depth           delta = sqrt(2 / (mu sigma omega))
                       same (real-valued)
  Complex permittivity epsilon = epsilon' - j epsilon''  (lossy)
                       <->   epsilon = epsilon' + i epsilon''  (lossy)

## Lab default: engineering convention

The lab (Sugahara) uses **`exp(+j omega t)` engineering convention**
throughout Radia C++, NGSolve panels, and all lab papers since 2018.

When porting equations from Chew, Brekhovskikh, Felsen-Marcuvitz, or
other physics-side sources: CONJUGATE EVERY COMPLEX QUANTITY.  When
porting equations from any modern IEEE / IEEJ paper: use as-is.

## Why this matters (lab reviewer pattern)

  Reviewer (IEEE TMag 2023):
    "In eq. (5) the author writes exp(-j omega t) on the left side and
     uses H_0^(1) on the right side as the outgoing wave.  These are
     INCONSISTENT.  Please choose either engineering (exp(+j omega t) +
     H_0^(2)) or physics (exp(-i omega t) + H_0^(1)) and apply it
     consistently throughout."

This is the SINGLE MOST COMMON theory-paper rejection cause.  When you
cite a paper from the "wrong" convention, you must EXPLICITLY state the
conversion in a footnote.

## How to verify in your paper

  1. Search the .tex for ``exp(j`` and ``exp(-j`` (engineering) AND
     ``exp(i`` and ``exp(-i`` (physics).  Mixing the two ⇒ likely bug.
  2. Search for ``H_0^{(1)}`` and ``H_0^{(2)}``.  If you use both,
     verify they label DIFFERENT physical objects (e.g. outgoing vs
     incoming), not the same one in two equations.
  3. Search for ``omega tau`` followed by `+` or `-` `j` / `i`.  All
     occurrences must share the same sign on the imaginary unit.
"""


VECTOR_TENSOR_NOTATION = r"""
# Vector and tensor notation in EM papers

## Vector notation -- 3 common styles

| Style | Example | Used by |
|-------|---------|---------|
| **Bold italic** | $\\boldsymbol{E}$, $\\boldsymbol{B}$ | IEEE (modern), IEEJ |
| Bold upright | $\\mathbf{E}$, $\\mathbf{B}$ | Physics (Griffiths) |
| Arrow over | $\\vec{E}$, $\\vec{B}$ | Older texts, blackboard |
| Underline | $\\underline{E}$, $\\underline{B}$ | German engineering |

Lab default: **bold italic** ($\\boldsymbol{E}$) per IEEE Author Center.
Stickier rule: WHATEVER you pick, use it CONSISTENTLY in the same paper.

## Tensor notation

| Style | Example | Used by |
|-------|---------|---------|
| Bold italic capital | $\\boldsymbol{T}$, $\\boldsymbol{\\sigma}$ | IEEE / IEEJ |
| Double-underline | $\\underline{\\underline{T}}$ | German engineering |
| Sans-serif bold | $\\mathsf{T}$ | Differential geometry |
| Component form | $T_{ij}$, $T^{ij}$ | All (preferred for index notation) |

Lab default: **bold italic capital** for symbolic tensor; component
form ($T_{ij}$, $\\nu_{ij}$ for the reluctivity tensor) for indexed
operations.  The CLN module uses ${{R}}_n$, ${{L}}_n$ (italic
non-bold) for SCALAR rungs and ${{\\boldsymbol{R}}}_n$,
${{\\boldsymbol{L}}}_n$ for MATRIX rungs.

## Math operators: italic VS upright

  upright: differentials (d/dx, $\\mathrm{d}x$), constants ($\\pi$,
           $\\mathrm{e}$, $\\mathrm{j}$ -- though IEEE allows italic j),
           operator names (det, tr, rank, sin, cos, exp, log)
  italic:  variables ($x, y, t, \\omega, \\sigma$ as variable),
           function names ($f(x)$, $g(\\omega)$)

LaTeX shortcuts:
  $\\sin$, $\\cos$, $\\exp$, $\\log$, $\\det$, $\\tr$, $\\rank$ -- all
  auto-upright via LaTeX kernel macros.  Do NOT write $sin$ (no
  backslash) -- LaTeX renders that as 3 italic letters $s i n$, which
  reviewers will flag.

## Differential forms / coordinate-independent notation

For the lab's differential-forms-flavoured work (Bossavit, FEEC,
``differential_forms`` MCP), use:

  * 1-forms (E, A):     $\\boldsymbol{E}$ (just the bold italic; the
                         form-degree is implicit from context)
  * 2-forms (B, J, D):   same notation; equations like
                         $\\mathrm{d}\\boldsymbol{B} = 0$ make the
                         degree explicit.
  * Hodge star:          $\\star$ (NOT $*$)
  * Wedge product:       $\\wedge$
  * Exterior derivative: $\\mathrm{d}$ (UPRIGHT to distinguish from
                         "delta" or variable d)

See `differential_forms_basics` for fuller convention.
"""


B_VS_H_USAGE = r"""
# B vs H: when to use which in EM papers

## Strict IEEE / IEEJ usage

| Quantity | SI unit | Strict term | Common (slightly imprecise) term |
|----------|---------|-------------|----------------------------------|
| $\\boldsymbol{B}$ | T (tesla) | magnetic flux density | "magnetic field" |
| $\\boldsymbol{H}$ | A/m | magnetic field intensity | "magnetic field" |
| $\\boldsymbol{M}$ | A/m | magnetization | "magnetization" |
| $\\boldsymbol{J}$ | T (tesla) | magnetic polarization | "magnetization" (in some communities) |

Lab convention:
  * **Always write B explicitly as "magnetic flux density"** in
    headings, abstract, captions.  Body text may use "magnetic field
    B" once context is established.
  * **Always write H explicitly as "magnetic field intensity"** OR
    "magnetic field strength" in headings.  Just "H" is OK in
    equations.
  * **DO NOT use "magnetization" for both M and J**.  Lab uses M
    (A/m) per the Radia C++ POLICY (see CLAUDE.md "Magnetization
    Units: A/m").

## Reviewer pattern

  Reviewer (IEEE TMag 2024, rejecting a lab draft):
    "The author calls B 'magnetic field' throughout the paper.  This
     is informal usage; the IEEE Editorial Style Manual requires
     'magnetic flux density' for B.  Please update headings, captions,
     and the abstract."

  Reviewer 2 (same paper):
    "Eq. (7) uses M with units of tesla, but eq. (12) treats M as A/m.
     The author appears to confuse magnetization M (A/m) with magnetic
     polarization J (T).  Please clarify."

The 2nd reviewer is right: M is A/m by SI convention.  Lab Radia uses
M=A/m consistently (Br=1.2 T → M=954930 A/m via M = Br / mu_0).

## How to lint

  Grep the .tex for:
    "magnetic field B"      -> flag (use "magnetic flux density")
    "magnetic field H"      -> flag (use "magnetic field intensity")
    "magnetization J"       -> flag (use "magnetic polarization")
    "polarization M"        -> flag (use "magnetization")

  Grep equations for:
    "M = .* [Tt]\\b"        -> M in tesla? almost certainly a J/M mixup
    "B = .* A/m"            -> B in A/m? definitely an H/B mixup
"""


SI_UNITS_NOTATION = r"""
# SI units notation in EM papers

## The 4 ironclad SI rules

  1. **Space between number and unit** (except %, °):
       "5 mm" -- correct
       "5mm" -- WRONG (catches "100kHz" too)
       "20 percent" or "20\\,\\%" -- correct
       "10\\,$^\\circ$" (angle) or "10$^\\circ$C" (temperature)

  2. **Unit symbols are UPRIGHT**, variables are italic:
       "$f = 50$ Hz"          -- correct (f italic, Hz upright)
       "$f = 50$ \\textit{Hz}" -- WRONG
       "$f = 50 Hz$"           -- WRONG (renders Hz in italic math)

  3. **Compound units use space, NOT center dot in body text**:
       "kg m/s$^2$" or "kg\\,m/s$^2$"  -- correct
       "kg$\\cdot$m/s$^2$"             -- formal (allowed)
       "kgm/s$^2$" without space        -- WRONG

  4. **Unit prefixes** are part of the symbol:
       "5 mm" -- the "m" is the milli prefix; "mm" is the symbol
       "5 nA/m" -- "n" is nano prefix, "A" is ampere, "m" is meter
       Never write "5 m m" (space inside the unit).

## EM-specific units checklist

  Quantity                Symbol      Unit
  ---------------------------------------------
  Frequency               f, omega    Hz, rad/s
  Magnetic flux density   B           T  (tesla)
  Magnetic field intens.  H           A/m  (NOT "Oersted" outside CGS)
  Magnetization           M           A/m
  Magnetic polarization   J           T
  Permeability            mu          H/m  (NOT "G/Oe")
  Permittivity            epsilon     F/m
  Conductivity            sigma       S/m
  Resistivity             rho         Ohm m  (NOT $\\Omega$.m)
  Inductance              L           H  (henry)
  Capacitance             C           F  (farad)
  Impedance               Z           Ohm
  Vector potential        A           Wb/m (or T*m)
  Current density         J           A/m^2

## Common errors caught by reviewers

  * Mixing CGS (Gauss, Oersted) with SI -- only do this if the
    paper is explicitly a comparison with old CGS literature.
  * Writing "Oersted" or "Oe" anywhere -- forbidden in modern
    IEEE papers; use A/m.
  * Lab convention: "tesla", "henry", "weber" lowercase when written
    out (units are not capitalised); symbols T, H, Wb uppercase.

## LaTeX siunitx macros (lab-recommended)

```latex
\\usepackage{siunitx}
\\sisetup{detect-all}        % use surrounding font style
% Then:
\\SI{50}{\\hertz}            % renders "50 Hz"
\\SI{5.8e7}{\\siemens\\per\\meter}   % "5.8 x 10^7 S/m"
\\SI{1.2}{\\tesla}           % "1.2 T"
\\num{5.8e7}                % "5.8 x 10^7"
```

`siunitx` enforces the space-between-number-and-unit rule automatically
and produces consistent output (e.g. \\num{5.8e7} renders as
$5.8 \\times 10^{7}$, not "5.8e+07").  Highly recommended -- a
single missing siunitx import is a reviewer-visible defect.
"""


EQUATION_TYPESETTING = r"""
# Equation typesetting for EM papers (IEEE / IEEJ rules)

## Number every displayed equation that is referenced

```latex
\\begin{equation}
\\nabla \\times \\boldsymbol{H} = \\boldsymbol{J} + \\frac{\\partial \\boldsymbol{D}}{\\partial t}
\\label{eq:ampere}
\\end{equation}
```

  * Use `equation` (auto-numbered) when you intend to \\ref it.
  * Use `equation*` (unnumbered) when the equation is purely
    illustrative and never referenced.
  * NEVER write a long equation inline as `$ ... $` if you intend to
    \\ref it -- displayed form only.

## Multi-line: align, NOT eqnarray

  ```latex
  \\begin{align}
  E_x &= \\sin(\\omega t - kz)  \\\\
  E_y &= 0  \\\\
  E_z &= 0
  \\end{align}
  ```

  `eqnarray` is DEPRECATED (poor spacing around `=`).  Always use
  `align`.

## Equation references

  IEEE convention: "(5)" with parens, NOT "Eq. (5)" or "equation 5".

  ```latex
  As shown in~(\\ref{eq:ampere}), the curl of H equals...
  ```

  Add tilde `~` between the word and the cite to prevent line-break
  separation between "in" and "(5)".

## Variable definitions BELOW the equation

  IEEE convention:
  ```latex
  \\begin{equation}
  P_{\\rm loss} = R I^2
  \\end{equation}
  where $P_{\\rm loss}$ is the dissipated power (W), $R$ is the
  resistance ($\\Omega$), and $I$ is the rms current (A).
  ```

  All variables that appear in an equation should be defined either
  (a) immediately before it ("Let $X$ denote ...") OR
  (b) immediately after it ("where $X$ is ...").
  Reviewers WILL flag undefined variables.

  For reduced-order / Schur-complement papers, this also applies to
  block matrices and model-order symbols.  Write frequency-domain blocks
  as `K_{bb}(s)`, `K_{bs}(s)`, `K_{sb}(s)`, `K_{ss}(s)` when they depend
  on `s`, and define them in words as bulk, surface, and coupling
  Galerkin blocks.  If the result is really "two-rung CLN + SIBC",
  say that directly; do not make the reader parse whether a uniform dc
  term is included in `N_b`.  Define `N_b` explicitly only when the
  symbol is needed: what kind of rung/mode/basis function it counts.
  If `N_b` and `N` both appear, state the distinction instead of relying
  on convention; e.g. `N_b` counts mixed-model bulk correction functions,
  whereas `N` is the order of the bulk-only CLN ladder.
  For surface-impedance expansions, define the surface order explicitly.
  In particular, `p_H=0` is the leading SIBC term, not a high-order
  impedance-boundary result.  If the manuscript also mentions HOIBC, make
  the separation clear: the present benchmark uses the leading SIBC term,
  while higher-order impedance terms are the extension needed in some
  curved or three-dimensional applications.
  If L/R CLN terminations appear in a figure or caption, define them:
  L-terminated means the last rung is closed inductively, and
  R-terminated means it is closed resistively.  State whether `N=10`
  is a rung, mode, or basis count.

  Do not equate a Schur complement with a Dirichlet-to-Neumann /
  Steklov-Poincare map in a single unexplained sentence.  The Schur
  complement is the algebraic result of eliminating unknowns; the DtN
  map is the interpretation that the resulting block maps surface
  Dirichlet data to Neumann flux.  Write those two steps explicitly.

  The core method equation should be numbered, labeled, and cited in
  prose.  Do not leave the main block system in `\\[ ... \\]` with no
  equation number, and do not number an equation that is never mentioned
  by `\\eqref{...}` / `\\ref{...}`.

## Units inside equations: AVOID

  Wrong:  $f = 50~\\mathrm{Hz}$
  Right:  $f = 50$ Hz   (move the unit OUTSIDE the math)
  Right:  $f = 50~\\mathrm{Hz}$ ONLY when the equation is showing
                                 the result of a unit-aware calculation

  The reason: variables in equations are pure NUMBERS in implicit SI
  unless the equation is symbolic.  Mixing numerical values with
  unit symbols inside math mode obscures dimensional analysis.

## Lab-specific: complex impedance / admittance

  Lab default:
  ```latex
  Z = R + j\\omega L   \\quad \\text{(engineering convention)}
  Y = G + j\\omega C
  ```

  $j$ NOT $i$.  $j$ italic is accepted by IEEE Author Center.
  $\\mathrm{j}$ (upright) is the formal alternative but visually
  busy; lab uses italic $j$.
"""


FIGURE_FORMAT = r"""
# Figure file format: import via PDF (vector), NOT PNG (raster)

## The rule

When including a generated graph / plot in the LaTeX manuscript, use
the **PDF (vector)** output, NOT the PNG (raster) output:

```latex
% CORRECT — vector PDF, scales losslessly, sharp at any zoom
\\includegraphics[width=\\columnwidth]{fig/result.pdf}

% WRONG — raster PNG, blurs when the typesetter scales it, fails
%         publisher 600/1200-dpi pre-flight on line art
\\includegraphics[width=\\columnwidth]{fig/result.png}
```

The lab figure pipeline (`mcp-server-figure` `emit_paper_figure`) emits
BOTH `<name>.pdf` and `<name>.png` at the journal's exact width.  The
`.png` is for quick preview / slides / GitHub README ONLY.  The `.pdf`
is the one that goes into the paper.

## Why PDF, not PNG, for line-art figures

| Aspect | PDF (vector) | PNG (raster) |
|--------|--------------|--------------|
| Scaling | lossless at any size | blurs / pixelates when scaled |
| Text in figure | stays selectable + sharp | fixed pixel grid, fuzzy |
| Thin strokes (0.7 pt axes) | crisp at print | alias / disappear |
| Publisher pre-flight | passes (Type-42 fonts embed) | flagged "low-res raster" on line art |
| File size (line art) | small | large at 600+ dpi |

## The ONE exception: inherently raster content

Use PNG / TIFF only for content that IS a pixel image:
- field maps / contour heatmaps rendered as a bitmap
- photographs of the experimental rig
- microscopy / camera images

For those, embed at >= 300 dpi (>= 600 dpi for mixed line+image) and
still prefer wrapping them so the AXES / LABELS are vector (matplotlib:
keep the colormesh raster but let the axis frame + ticks + labels be
vector by saving the whole figure as PDF — `rasterized=True` on just
the `pcolormesh`/`imshow` artist gives a small vector-frame + raster-
data hybrid PDF, the best of both).

## How to verify in your manuscript

  1. Grep the .tex for ``\\includegraphics`` and check every path ends
     in ``.pdf`` (or ``.eps``), not ``.png`` / ``.jpg`` — except the
     deliberate raster exceptions above.
  2. If a figure looks fuzzy in the compiled PDF, it is almost always
     a PNG that should have been the PDF sibling.
  3. The lab `emit_paper_figure(...)` already produced the `.pdf`; just
     point `\\includegraphics` at it.

## Cross-reference

- `mcp-server-figure` `emit_paper_figure` — emits .pdf + .png at exact
  journal width; the .pdf is the manuscript figure.
- The companion `lab-graph-style` skill documents the full figure
  pipeline.
"""


FIGURE_CAPTION_BODY_POLICY = r"""
# Figure caption and body-text consistency

## The rule

Everything important enough to write in a figure or table caption must
also be stated in the main text.  Captions should help the reader identify
what is plotted, the essential conditions, and the meaning of symbols or
series, but they are not the primary place for claims, interpretation, or
conclusions.

If a caption becomes long, prioritize the body text:

  * keep the caption short enough to identify the plotted quantity,
    geometry, parameter sweep, and series;
  * move comparison logic, physical interpretation, numerical scaling
    claims, and conclusions into the paragraph that cites the figure;
  * do not leave a number, trend, or "therefore" statement only in the
    caption;
  * for dense digests, one body sentence after the figure reference is
    usually better than a caption that tries to be a mini Results section.

## Good manuscript pattern

Body:
  "Fig. 5 summarizes the MMPM-HACApK scaling.  The large-DOF-side fit
  gives memory scaling of about O(n^1.2) and solve-time scaling of about
  O(n^1.6)."

Caption:
  "Scaling of MMPM solvers for the fixed-boundary cube benchmark.  Cases
  above 100,000 DOF were evaluated only with HACApK."

## Bad manuscript pattern

Body:
  "Fig. 5 shows the scaling."

Caption:
  "Scaling of MMPM solvers... The large-DOF-side fit gives memory scaling
  of O(n^1.2) and solve-time scaling of O(n^1.6), demonstrating scalable
  performance for practical models."

This fails because the main claim appears only in the caption.  Reviewers
who read the body first miss the conclusion; readers who skim figures get
a result that is not integrated into the argument.

## Cross-reference

- `paper_writing_check_figure_caption_showing` reports the caption-body
  policy in its recommendation.
- `paper_writing_check_digest_human_review_triggers` flags overloaded
  captions and points interpretation back to the body.
"""


ITERATION_COMPARISON_POLICY = r"""
# Fixed-iteration reporting and fair solver comparison

## Disclose a fixed iteration budget

A fixed linear or nonlinear iteration count is a legitimate experimental
protocol, but it must not be hidden.  State in the Methods section:

  * the update sequence and count (for example, quasi-Newton followed by
    two Picard blocks);
  * the number of inner linear iterations, damping/relaxation, and
    preconditioner;
  * the residual or physical-observable acceptance threshold; and
  * how the reported state is selected when the best state seen so far is
    retained.

Describe the count as a **fixed computational budget**, not as a
convergence count.  Report an acceptance check as well: a fixed budget by
itself does not show that the physical solution is sufficiently accurate.

Observable-based selection can support a claim of iteration stability or
mesh-to-mesh consistency when the residual gate and observable change are
reported.  It does not establish absolute accuracy without an analytical
solution, independent reference solution, or measurement.  Do not turn
"the observable is less iteration-sensitive than the internal state" into
"the observable is exact."

## Compare iteration counts only under identical semantics

Iteration counts are directly comparable only when the methods use the
same update definition, initial state, residual norm, tolerance, stopping
rule, and maximum-iteration policy.  If one method stops on a residual and
another uses a fixed budget or an observable-based best-state rule:

  * do not place the counts in one performance-comparison column;
  * do not claim that the smaller count is faster or more robust;
  * report the protocols separately in Methods or supplementary material;
  * compare physical quantities, accuracy, measured time, and measured
    memory under the declared protocols.

## Interpret scaling with a fixed outer count carefully

When every problem size uses the same outer iteration budget, say why:
the experiment removes problem-size-dependent nonlinear iteration growth
to isolate matrix assembly, factorization, or matrix-vector-product
scaling.  A fitted exponent from that experiment is the scaling of the
measured fixed-budget workload.  It is not evidence that the complete
nonlinear convergence cost has the same asymptotic exponent.

## Recommended manuscript pattern

Methods:
  "To isolate matrix-operation scaling, the same outer-update budget was
  used for all problem sizes.  Among states satisfying the nonlinear
  residual threshold, the state with the smallest physical-observable
  change was retained."

Results table:
  Show the physical quantities and accuracy.  Omit an iteration-count
  column when the compared formulations use different stopping semantics.
"""


BILINGUAL_WORKFLOW = r"""
# Bilingual workflow: English paper + Japanese translation

## The rule (Sugahara Lab)

When writing an **English** paper, ALSO produce a **Japanese version**.
The Japanese version is **generated AS the translation of the English
manuscript** (English is the master; Japanese is derived from it), and
serves as a lab deliverable (internal review, co-author readability,
archival).

**Key point — the Japanese version is INDIFFERENT to the page limit.**
The page limit is a property of the TARGET VENUE (IEEE TMag 8 pp, IGTE
digest 2 pp, etc.) and applies ONLY to the English submission.  For the
Japanese translation you simply **do not care about page count at all**
— do not measure it, do not enforce it, do not adjust the Japanese to
hit any page number.  It runs to whatever length the translation
naturally takes.

This is NOT a rule "make the Japanese fit / not fit one page".  Page
count is just not a concern for the Japanese version, in either
direction.

## What this means in practice

| Artifact | Page limit? | Purpose |
|----------|-------------|---------|
| English manuscript (submitted) | **YES** — venue page limit enforced | the actual submission |
| Japanese translation (derived) | **N/A — indifferent** — any length | lab review, co-author readability, archive |

- Run the page-limit check (`paper_writing_validate_pdf_pages`,
  `paper_writing_em_submission_gate(page_limit=...)`) **only on the
  English PDF**.
- For the Japanese translation, simply **do not pass `page_limit`** —
  page count is not evaluated for it at all (its length is whatever the
  translation produces; that is neither good nor bad).
- Figures, equations, and tables are SHARED between the two versions
  (same .pdf figure files); only the prose is translated.  This keeps
  the two in sync and means the figure-format rule
  (`figure_format` topic) applies identically to both.

## Why keep a Japanese translation at all

- Co-authors / students who read Japanese faster catch logic errors
  the English draft hides.
- The Japanese version is reusable for a domestic-conference
  (IEEJ 全国大会 / 研究会) submission later.
- Archival: the lab keeps both language versions of every paper.

## How to apply

  1. Draft + finalize the English manuscript to the venue page limit
     (English is the master document).
  2. Generate the Japanese version AS a translation of the FINAL
     English text (translate the prose; reuse the identical figure
     .pdf files and equations verbatim).  Let it run to whatever
     length results — do not trim or pad it toward any page count.
  3. Gate the English PDF WITH `page_limit`; gate the Japanese PDF
     WITHOUT `page_limit` (page count is simply not checked for the
     Japanese version).  All OTHER checks — undefined variables,
     ref/label consistency, figure format — still apply to both.

## Cross-reference

- `paper_writing_em_submission_gate(pdf_path=..., page_limit=N)` —
  pass page_limit ONLY for the English PDF.
- `figure_format` topic — figures are shared; both versions import
  the same vector .pdf figures.
"""


ABSTRACT_AND_CONCISENESS = r"""
# Abstract content rules + page-fitting conciseness philosophy

## Abstract rules (Sugahara Lab / IEEE / IEEJ)

1. **NO equations / math in the abstract.**
   The abstract is plain prose.  Do NOT put displayed equations,
   inline math (`$...$`), or symbol-heavy expressions in it.  State
   results in words ("the eddy-current loss is reduced by 18 %"), not
   formulas.  Reason: the abstract is indexed / displayed by databases
   (IEEE Xplore, Scopus, Google Scholar) as plain text; math renders as
   garbage there, and reviewers expect a math-free abstract.

2. **NO citations in the abstract.**
   Do NOT put `\cite{}` / `[1]` reference markers in the abstract.
   The abstract must stand alone; a reference number is meaningless to
   a reader who hasn't reached the bibliography.  If prior work must be
   acknowledged, name it in words ("unlike conventional Preisach
   models", "the Hollaus et al. method") without a citation marker.

   **No exception — not even for the foundational paper the work
   re-casts, nor for a co-author's paper.**  The tempting loophole is
   "keep `... of X et al.~\cite{X}` so the foundational reference is
   numbered [1]".  Resist it: name X et al. in words in the abstract
   and move the `\cite{}` to the first body mention (Introduction /
   Method).  The reference still renders — it is simply numbered by its
   body order of appearance (so it may become [4] rather than [1], which
   is correct and fine).  Do NOT carve out a "keep [1] in the abstract"
   exception.  (Real 2026 IGTE-digest slip: a co-author's foundational
   ref was left as `\cite{}`/[1] in the abstract "as a deliberate
   exception" and had to be stripped on review — the renumbering to a
   body-order number is the intended outcome, not a problem to avoid.)

3. **NO domain-specific acronyms in the abstract.**
   Avoid abbreviations such as FEM, BEM, MCP, or LLM in the
   abstract, even if they are familiar inside the lab.  Write the
   abstract in plain words, then introduce acronyms at their first body
   use: `Finite Element Method (FEM)`,
   etc.  Universal acronyms such as PDF, USB, or CPU are ignored by the
   checker.

These rules are enforced by
`paper_writing_check_abstract_no_math_no_citation` (and the
`em_submission_gate` runs it automatically).  Violations are a
reviewer-visible, database-visible defect.

## Page-fitting philosophy: don't OVER-compress

When fitting a hard page limit (IGTE / Compumag 2-page digest, a
1-page extended abstract, IEEE TMag 8-page full paper), there are two
ways to make content fit:

| Approach | Verdict |
|----------|---------|
| **Cram**: shrink fonts below the 10 pt rule, squeeze margins, kill whitespace, pack every sentence | ✗ AVOID — hurts readability, breaks the 10 pt-at-8 cm figure font rule, reviewers complain |
| **Select**: reduce SCOPE — cut secondary results, move detail to a future full paper, keep fewer points but keep them readable | ✓ PREFERRED |

**Rule**: to fit one page (or any page limit), you do NOT have to
over-compress the information.  It is equally valid — and usually
better — to **reduce the amount of content** (be selective) rather than
cram everything in at the cost of readability.  Keep the surviving
content at full readability (10 pt fonts, adequate whitespace, figures
at proper size); drop or defer the rest.

Why this matters:
- A digest is a TEASER, not a compressed full paper.  It should leave
  the reader wanting the full paper, not exhaust every detail.
- Over-compression triggers the reviewer comment "the figures/text are
  too small to read at print scale" (see `reviewer_patterns` #7).
- The 10 pt-at-8 cm figure font rule (graph pipeline) is a HARD floor;
  if content won't fit at 10 pt, cut content, do not shrink the font.

## How to apply

  1. Draft the abstract in plain prose; remove every `$...$`,
     `\cite{}` / `[1]`, and domain-specific acronym — verify with
     `paper_writing_check_abstract_no_math_no_citation`.
  2. If over the page limit, FIRST cut scope (secondary results,
     redundant figures), THEN check fit — do not reach for smaller
     fonts / tighter margins as the primary tool.
  3. Re-gate; the figure font rule and whitespace checks confirm you
     did not over-compress.

## Cross-reference

- `paper_writing_check_abstract_no_math_no_citation` — enforces rules 1-3.
- `paper_writing_validate_abstract_length` — abstract length.
- `bilingual` topic — page limit applies to the EN submission only.
- `reviewer_patterns` #7 — "figure font too small" over-compression flag.
"""


HUMAN_REVIEW_LEARNING_POLICY = r"""
# Human-review learning policy

When a coauthor or human review finds a problem that is generic rather
than manuscript-specific, do not treat it as a one-off edit.  Convert it
into the paper-writing knowledge layer.

## What counts as generic

- Abstract conventions that will recur across venues: no math, no
  citation markers, no local acronyms, no lists of exact percent errors.
- Digest structure problems: main result figures floating before the
  section that explains them, a numerical example hidden inside a bold
  paragraph instead of a numbered section, or a one-page digest promising
  too many future geometries.
- Figure placement policy: use the standard LaTeX `figure` environment
  whenever possible and place it in the source immediately after the
  relevant section heading / first reference.  The `float` package with
  `[H]`, manual `\refstepcounter{figure}` captions, non-float minipage
  figures, and forced `\clearpage` / `\newpage` placement are last
  resorts.  First reduce text, figure width, caption length, or move the
  `figure` environment to the correct source position.
- Figure readability problems: over-compressed aspect ratio, result
  labels such as `f_N` that do not state what they mark, or captions that
  carry interpretation that belongs in the text.  If a caption contains
  a factual claim, numerical result, trend, or conclusion, write the same
  point in the main text too.  If the caption becomes long, prioritize the
  body and shorten the caption.
- Reproducibility omissions in standard benchmarks: missing material
  properties, reference solutions, model order, or boundary-condition
  order.
- Proposal / research-plan framing problems: section titles written as
  the author's drafting process ("notes", "plain explanation") instead
  of the reader-facing theoretical or design role; immature future work
  shown as a result-like integrated graph instead of a workflow or
  hypothesis; internal tool or MCP implementation notes left in the
  public manuscript.
- Proposal / research-plan figure selection: at the idea stage, figures
  should either ground the proposal in a real calculation / measurement
  or explain the iteration that will test an unverified hypothesis.  Do
  not write defensively that a combined figure is "not ready"; state what
  each shown figure contributes.  Keep unverified implementation ideas in
  an explicitly marked inventory table or hypothesis list, not in a
  result-like conclusion figure.
- Coordinate-transform explanations that overclaim "the nonlinearity is
  gone" without saying what is invariant and what changed.  State the
  unchanged operator / conservation law and the changed metric, Hodge
  operator, coefficient, or unknown immediately after the equation.
- Claim strength and antecedent problems: a finite numerical sweep should
  normally "support", "suggest", or "provide a finding", not serve as
  "proof/証拠" of a broader mechanism.  At paragraph boundaries, repeat the
  named method in phrases such as "the derived MMPM moment conditions";
  do not make the reader infer what "the derived conditions" belong to.
- Iteration-reporting fairness: disclose fixed linear/nonlinear budgets,
  relaxation and accepted-state criteria.  Do not compare a fixed budget
  with a residual-based convergence count in one iteration column.  For a
  fixed-budget scaling sweep, state that the design isolates matrix-work
  scaling and does not establish end-to-end nonlinear convergence
  complexity.

## Required action

For each generic lesson, update at least one durable artifact:

1. `paper_writing/skill.md` for prose policy.
2. `_em_paper_style.py` when the rule should be returned by
   `paper_writing_em_paper_style(...)`.
3. `paper_writing_check_*` plus bad/clean tests when the pattern is
   mechanically detectable.

Do not generalize manuscript-specific numbers, author preferences, or
unconfirmed future-work promises.  A full-paper outlook such as
"extension to three-dimensional conductors will be examined" is allowed
only when it is a real plan and is clearly framed as future scope, not as
a current digest result.
"""


REFERENCE_BIB_POLICY = r"""
# references.bib policy: shared across digest + full paper, web-verified

## The three rules (Sugahara Lab)

1. **ONE shared `references.bib`** for BOTH the digest (IGTE / Compumag
   2-page) AND the full paper (IEEE TMag / IEEJ Trans).  Do NOT keep a
   separate `.bib` per document.  The digest and the full paper of the
   same research cite the same literature; a single `references.bib` is
   the single source of truth, kept in the shared project folder and
   `\bibliography{reference}`-d from both `.tex` files.

2. **EVERY entry MUST be web-verified** before it goes into the paper.
   No citation is inserted from memory / a hand-typed guess.  Each
   entry is confirmed against an authoritative web source (Crossref
   DOI, IEEE Xplore, Semantic Scholar, arXiv) so the authors, title,
   venue, year, volume, pages, and DOI are EXACT.

3. **Numbered reference lists follow first citation order.**  For IEEE
   styles and manual `thebibliography`, order entries by the first
   appearance of their `\cite{}` key in the text.  After moving related
   work sections, tables, or citations, re-check the bibliography order;
   do not leave a mixture of author order, year order, and drafting
   order.

## Why one shared .bib

| Benefit | Detail |
|---------|--------|
| No drift | digest + full paper can't disagree on a reference's year/pages |
| Less work | verify once, cite in both documents |
| Promotion path | a digest grows into a full paper; the .bib carries over verbatim |
| Audit | one file to lint (`paper_writing_lint_reference_format`) |

Layout:
```
project/
  references.bib          <- THE shared bibliography (single source)
  digest/igte_digest.tex     \bibliography{../reference}
  full/ieee_tmag.tex         \bibliography{../reference}
```

## Why web-verification is mandatory

Fabricated or misremembered citations are a top reviewer-trust killer
(wrong year, wrong volume, author initials swapped, a DOI that resolves
to a different paper, or an entirely hallucinated entry).  The lab rule
is **strict-in, no fabrication**: if a reference cannot be confirmed on
the web, it does NOT go in the .bib.

Verification workflow (per entry):
1. Search the web for the exact title.
2. Resolve the DOI via Crossref (`paper_writing_verify_citation` /
   `paper_writing_resolve_doi` / `paper_writing_doi_to_bibtex`).
3. Confirm authors / venue / year / volume / pages match the
   authoritative record (IEEE Xplore, publisher page, Semantic
   Scholar, arXiv).
4. Only THEN add the verified BibTeX entry to the shared
   `references.bib`.

The `paper_writing_verify_citation` tool **refuses to fabricate** — if
given no bib / no resolvable source it returns `verdict=error` rather
than inventing an entry.  This is by design.

## Enforcement in the submission gate

`paper_writing_em_submission_gate` already has a **bib policy gate**:
calling it WITHOUT `bib_path` returns `verdict=fail` with the message
that every citation must trace back to the real `references.bib` and be
verified.  Pass the SAME shared `references.bib` when gating BOTH the
digest and the full paper.

## How to apply

  1. Keep one `references.bib` at the project root; both `.tex` point to it.
  2. For every new reference: web-search -> verify DOI -> add entry.
     Never paste an unverified entry.
  3. For numbered citation styles, confirm the bibliography follows the
     first citation order in the compiled text.
  4. Gate digest AND full paper with the SAME `bib_path=references.bib`.
  5. Lint with `paper_writing_lint_reference_format(references.bib)` and
     `paper_writing_check_citation_keys_exist(tex, references.bib)` for
     each document — every `\cite{}` key must exist in the shared bib.

## Cross-reference

- `paper_writing_verify_citation` — web/DOI verification (refuses to
  fabricate).
- `paper_writing_citation_workflow_recipe` — full per-entry recipe.
- `paper_writing_check_citation_keys_exist` — every \cite resolves in
  the shared bib.
- `bilingual` topic — the Japanese translation reuses the SAME shared
  `references.bib` too (references are language-independent).
"""


EM_REVIEWER_PATTERNS = r"""
# Common EM-reviewer comments (and how to preempt them)

Collected from lab paper-review history 2019-2025 across IEEE TAP,
IEEE TMag, IEEE TMTT, IEEJ Trans. D / B.  For each pattern, the
recommended PRE-SUBMISSION check.

## 1. "Sign convention is inconsistent"

  See `paper_writing_em_paper_style('sign_conventions')`.
  PRE-CHECK: grep .tex for ``exp(-?j`` and ``H_0\\^\\{(1|2)\\}``.

## 2. "Magnetic field B vs magnetic flux density B"

  See `paper_writing_em_paper_style('b_vs_h')`.
  PRE-CHECK: grep ``magnetic field B`` (should be "magnetic flux
  density"); grep ``magnetic field H`` (should be "magnetic field
  intensity" or "magnetic field strength").

## 3. "Mesh independence not demonstrated"

  Lab pattern: include a convergence table (DOF vs result) in the
  Results section.  Three rows is the minimum: coarse / medium /
  fine.
  PRE-CHECK: search for "convergence", "DOF", "mesh refinement";
  if all three are absent, add the convergence table.

## 4. "What is the computational cost (CPU time / memory) of the
       proposed method vs the baseline?"

  Lab pattern: include a single table with N_DOF, T_solve, T_assemble,
  peak memory.  Reviewers especially hate FEM/BEM papers that don't
  report cost.
  PRE-CHECK: search for ``CPU`` and ``memory`` in Results section;
  if absent, add the cost table.

## 5. "How does this differ from prior work [X]?"

  Lab pattern: include an explicit comparison table in Introduction
  (rows: method axes; columns: lab vs prior X vs prior Y).
  PRE-CHECK: ensure Introduction has a sentence "The proposed method
  differs from [X] in that...".

## 6. "Why is the proposed method preferred over [closed-form
       analytical solution]?"

  EM papers often have an analytical closed-form competitor (image
  theory, infinite-plate analytical, Kelvin-image method, etc.).
  Reviewers WILL ask "why didn't you just use that?".
  Answer in the paper: cite the analytical solution + state its
  range-of-validity limitation that the proposed method overcomes.
  PRE-CHECK: use `radia_mcp.radia_ngsolve.analytical_formulas
  ('validation_use_cases')` to find which closed forms apply, and
  cite them explicitly.

## 7. "The figure font is too small to read at print scale"

  Lab pattern: enforce the 10pt-at-8cm rule (see
  `mcp-server-figure.paper_figure_quality_rules('font_rule')`).
  PRE-CHECK: use `paper_writing_layout_thumbnail_strip` and visually
  scan; or compute via `paper_writing_detect_page_whitespace_anomalies`.

## 8. "Citations to recent (last 5 years) literature are missing"

  Lab pattern: many lab papers cite the 2000s-era classical
  references (Bossavit, Chew, Wakao) without recent extensions.
  PRE-CHECK: use `paper_writing_semantic_scholar_citations` on the
  most-cited prior-art paper to find recent extensions you may have
  missed.

## 9. "Uncertainty / error bars are not shown"

  Lab pattern: experimental validation figures should include error
  bars (or shaded uncertainty regions).
  PRE-CHECK: search figure captions for "uncertainty", "error",
  "standard deviation".  If experimental data is plotted without
  error bars, add them or explain why.

## 10. "Reviewer 2: The proposed method is incremental, not novel"

  Lab pattern: this is the #1 rejection cause for IEEE TMag.
  PRE-CHECK: `paper_writing_contribution_clarity_score` and
  `paper_writing_reviewer_2_trigger_summary`.  Strengthen the
  Introduction's "Contribution" paragraph (3-5 bullets, each
  starting "We propose..." or "We derive..." or "We demonstrate...").

## 11. "Equations (5) and (12) appear contradictory"

  When equation (5) and (12) use the same SYMBOL for different
  PHYSICAL quantities (e.g. F = force in (5) but F = matrix in (12)).
  PRE-CHECK: build a symbol table (a "Nomenclature" section
  preceding Section I, with EVERY symbol used in the paper).  IEEE
  templates support this via `\\begin{IEEEdescription}`.

## 12. "Reference [X] is not the correct citation for this claim"

  Reviewer expertise often catches when a lab paper cites a vaguely
  related work instead of the specific paper that originated the
  claim.
  PRE-CHECK: use `paper_writing_semantic_scholar_references` on
  the cited paper to check it actually contains the cited result.
"""


TOPICS = {
    "overview":                  None,  # placeholder; resolved in dispatcher
    "sign":                      SIGN_CONVENTIONS,
    "sign_conventions":          SIGN_CONVENTIONS,
    "phase_convention":          SIGN_CONVENTIONS,
    "vector":                    VECTOR_TENSOR_NOTATION,
    "tensor":                    VECTOR_TENSOR_NOTATION,
    "notation":                  VECTOR_TENSOR_NOTATION,
    "vector_tensor_notation":    VECTOR_TENSOR_NOTATION,
    "b_vs_h":                    B_VS_H_USAGE,
    "magnetic_field_terminology": B_VS_H_USAGE,
    "units":                     SI_UNITS_NOTATION,
    "si_units":                  SI_UNITS_NOTATION,
    "siunitx":                   SI_UNITS_NOTATION,
    "equations":                 EQUATION_TYPESETTING,
    "equation_typesetting":      EQUATION_TYPESETTING,
    "displayed_math":            EQUATION_TYPESETTING,
    "figure_format":             FIGURE_FORMAT,
    "figures":                   FIGURE_FORMAT,
    "pdf_not_png":               FIGURE_FORMAT,
    "includegraphics":           FIGURE_FORMAT,
    "caption_body":              FIGURE_CAPTION_BODY_POLICY,
    "caption_body_policy":       FIGURE_CAPTION_BODY_POLICY,
    "caption_body_consistency":  FIGURE_CAPTION_BODY_POLICY,
    "figure_caption":            FIGURE_CAPTION_BODY_POLICY,
    "captions":                  FIGURE_CAPTION_BODY_POLICY,
    "iteration":                 ITERATION_COMPARISON_POLICY,
    "iterations":                ITERATION_COMPARISON_POLICY,
    "iteration_fairness":        ITERATION_COMPARISON_POLICY,
    "fixed_iteration":           ITERATION_COMPARISON_POLICY,
    "fixed_iterations":          ITERATION_COMPARISON_POLICY,
    "nonlinear_solver_reporting": ITERATION_COMPARISON_POLICY,
    "stopping_criteria":         ITERATION_COMPARISON_POLICY,
    "abstract":                  ABSTRACT_AND_CONCISENESS,
    "abstract_rules":            ABSTRACT_AND_CONCISENESS,
    "conciseness":               ABSTRACT_AND_CONCISENESS,
    "page_fitting":              ABSTRACT_AND_CONCISENESS,
    "no_math_abstract":          ABSTRACT_AND_CONCISENESS,
    "no_citation":               ABSTRACT_AND_CONCISENESS,
    "no_citation_abstract":      ABSTRACT_AND_CONCISENESS,
    "no_cite":                   ABSTRACT_AND_CONCISENESS,
    "citation_in_abstract":      ABSTRACT_AND_CONCISENESS,
    "human_review_learning":     HUMAN_REVIEW_LEARNING_POLICY,
    "review_learning":           HUMAN_REVIEW_LEARNING_POLICY,
    "policy_learning":           HUMAN_REVIEW_LEARNING_POLICY,
    "bilingual":                 BILINGUAL_WORKFLOW,
    "bilingual_workflow":        BILINGUAL_WORKFLOW,
    "japanese_translation":      BILINGUAL_WORKFLOW,
    "translation":               BILINGUAL_WORKFLOW,
    "page_limit":                BILINGUAL_WORKFLOW,
    "reference_bib":             REFERENCE_BIB_POLICY,
    "bib":                       REFERENCE_BIB_POLICY,
    "bib_policy":                REFERENCE_BIB_POLICY,
    "shared_bib":                REFERENCE_BIB_POLICY,
    "citation_verification":     REFERENCE_BIB_POLICY,
    "reviewer_patterns":         EM_REVIEWER_PATTERNS,
    "common_reviewer_comments":  EM_REVIEWER_PATTERNS,
    "reviewer_comments":         EM_REVIEWER_PATTERNS,
}


_OVERVIEW = r"""
# EM-paper style knowledge layer for radia_mcp.paper_writing

The generic text/lint tools (kanji ratio, hedge count, IMRaD balance,
etc.) are journal-agnostic.  THIS module is EM-paper-specific:

  sign / sign_conventions / phase_convention
        Engineering exp(+j omega t) vs physics exp(-i omega t).
        Conjugation rule for porting between conventions.  Lab default
        is engineering.

  notation / vector / tensor
        Bold-italic E vs arrow-over E vs underline E.  Tensor
        conventions.  Operators upright vs italic.  Differential-forms
        notation cross-link.

  b_vs_h
        "magnetic flux density B" vs "magnetic field intensity H"
        terminology.  Magnetization M (A/m) vs polarization J (T).

  units / si_units / siunitx
        Space between number and unit, upright units vs italic
        variables, compound-unit dot, EM-specific unit checklist
        (T, A/m, H/m, S/m, etc.), siunitx package recommendation.

  equations / equation_typesetting / displayed_math
        equation vs equation*, align (NOT eqnarray), \\ref{eq:foo}
        with tilde, variable-definition-after-equation rule.

  figure_format / figures / pdf_not_png / includegraphics
        Import figures via PDF (vector), NOT PNG (raster).  The lab
        graph pipeline emits both; the .pdf is the manuscript figure.
        PNG only for inherently-raster content (field maps, photos).

  caption_body / figure_caption / captions
        Anything important enough to write in a caption must also be in
        the main text.  If the caption grows long, keep the caption short
        and move interpretation, numerical claims, and conclusions to the
        paragraph that cites the figure.

  iteration_fairness / fixed_iteration / stopping_criteria
        Disclose fixed linear/nonlinear budgets and accepted-state rules.
        Compare iteration counts only under identical stopping semantics;
        interpret fixed-budget scaling as matrix-work scaling, not complete
        nonlinear convergence complexity.

  abstract / abstract_rules / conciseness / page_fitting
        Abstract carries NO math and NO citations.  To fit a page
        limit, reduce SCOPE rather than over-compressing (don't shrink
        fonts / cram); keep surviving content fully readable.

  bilingual / japanese_translation / translation / page_limit
        English paper + Japanese translation lab workflow.  The page
        limit applies ONLY to the English submission, NOT the Japanese
        translation (figures/equations shared between the two).

  reference_bib / bib / bib_policy / shared_bib / citation_verification
        ONE shared references.bib for digest + full paper; EVERY entry
        web-verified (Crossref / IEEE Xplore / S2 / arXiv) before
        insertion -- no fabricated citations.

  reviewer_patterns / reviewer_comments
        12 common EM-paper reviewer-comment patterns (sign
        inconsistency, B-vs-H, mesh independence, CPU/memory,
        contribution clarity, missing recent citations, etc.) and
        the pre-check tool that catches each.

  all
        Concatenate all topic sections.

Use as a CHECKLIST before submission, alongside the algorithmic checks
in tools.py and the new v0.88-0.89 layout / arxiv tools.
"""


def paper_writing_em_paper_style(topic: str = "overview") -> str:
    """EM-domain-specific style/notation/convention knowledge.

    Companion to the generic checks in paper_writing.tools.  Captures
    the rules that EM reviewers (IEEE TAP/TMag/TMTT, IEEJ Trans D/B,
    and accelerator-community venues) actually enforce, plus the
    12-pattern catalogue of common
    reviewer comments collected from lab paper-review history.

    Args:
        topic: one of TOPICS keys.  Aliases accepted.

    Available topics:
        overview                              -- this index
        sign / sign_conventions               -- exp(+jwt) vs exp(-iwt)
        vector / tensor / notation            -- bold italic vs arrow
        b_vs_h / magnetic_field_terminology   -- B "flux density" vs H
        units / si_units / siunitx            -- SI unit notation
        equations / equation_typesetting      -- align, \\ref, units
        figure_format / pdf_not_png           -- import figures as PDF
        caption_body / figure_caption         -- claims in captions also
                                                  belong in the body
        iteration_fairness / fixed_iteration  -- fixed-budget disclosure
                                                  and fair count comparison
        abstract / conciseness / page_fitting -- no math/cite in abstract;
                                                 cut scope, don't over-compress
        bilingual / japanese_translation      -- EN paper + JA translation,
                                                 page limit EN-only
        reference_bib / bib_policy            -- shared references.bib for
                                                 digest+full, web-verified
        reviewer_patterns / reviewer_comments -- 12-pattern catalogue
        all                                    -- concatenate everything
    """
    key = topic.strip().lower()
    if key == "overview":
        return _OVERVIEW
    if key == "all":
        return "\n\n---\n\n".join([
            _OVERVIEW,
            SIGN_CONVENTIONS,
            VECTOR_TENSOR_NOTATION,
            B_VS_H_USAGE,
            SI_UNITS_NOTATION,
            EQUATION_TYPESETTING,
            FIGURE_FORMAT,
            FIGURE_CAPTION_BODY_POLICY,
            ITERATION_COMPARISON_POLICY,
            ABSTRACT_AND_CONCISENESS,
            BILINGUAL_WORKFLOW,
            REFERENCE_BIB_POLICY,
            EM_REVIEWER_PATTERNS,
        ])
    if key not in TOPICS or TOPICS[key] is None:
        return (f"Unknown topic '{topic}'. Available: "
                f"{', '.join(sorted(set(k for k, v in TOPICS.items() if v is not None) | {'overview', 'all'}))}.")
    return TOPICS[key]


# ============================================================
# Pass 3: pre-submission gate orchestrator
# ============================================================


def paper_writing_em_submission_gate(
    tex_path: str = "",
    pdf_path: str = "",
    bib_path: str = "",
    abstract_text: str = "",
    author_last_names: str = "",
    target_venue: str = "",
    page_limit: int = 0,
    whitespace_threshold: float = 0.75,
    layout_max_pages_apart: int = 1,
    auto_extract_abstract: bool = True,
    auto_resolve_inputs: bool = True,
) -> dict:
    """One-shot EM-paper pre-submission gate.

    Runs ALL relevant checks in sequence on the user's manuscript +
    typeset PDF + bibliography, then returns a single pass/fail
    verdict + a ranked-by-severity issue list.

    Args:
        tex_path: path to the main .tex file (required for most checks).
        pdf_path: path to the typeset PDF (required for layout checks).
        bib_path: path to the .bib file (required for citation checks).
        abstract_text: the abstract verbatim (required for abstract
            checks; can be auto-extracted from tex_path in a future
            extension).
        author_last_names: comma-separated lastnames for self-citation
            check (e.g. "Sugahara,Nagamine").
        target_venue: concrete target venue. CAE-AI Lab policy accepts
            only IEEJ, IEEE, and accelerator-community venues.
        page_limit: page limit for the target venue (0 = skip the check).
        whitespace_threshold: threshold for whitespace flag (default 0.75).
        layout_max_pages_apart: max pages a float can drift from its
            first reference (default 1).

    Returns:
        dict with:
          "verdict": "pass" | "warn" | "fail"
          "n_checks_run": int
          "n_critical": int
          "n_warning": int
          "checks": list of {"name": str, "status": "pass"|"warn"|"fail"|"skip",
                              "summary": str, "detail": dict}
          "advice": text recommendation
    """
    from . import tools as _t
    from ._pdf_layout_visual import (
        paper_writing_detect_page_whitespace_anomalies,
        paper_writing_check_floats_far_from_reference,
    )
    # v0.92.0: PDF overlap + undefined-variable checks
    from ._pdf_overlap_detection import (
        paper_writing_detect_text_image_overlap,
        paper_writing_detect_text_overflow_page,
    )
    from ._undefined_variables import (
        paper_writing_check_undefined_variables,
    )
    # v0.93.0: multi-file .tex resolver + abstract auto-extract
    from ._tex_resolver import (
        resolve_input_chain,
        extract_abstract_from_tex,
    )
    # 2026-05-27: digest-review additions (acronym + citation key existence
    # + ref/label + IEEE keywords + PDF unresolved markers + abstract math/cite)
    from ._undefined_acronyms import (
        paper_writing_check_undefined_acronyms,
    )
    from ._citation_verify import (
        paper_writing_check_citation_keys_exist,
    )
    from ._digest_lints import (
        paper_writing_check_ref_label_consistency,
        paper_writing_check_ieee_keywords,
        paper_writing_check_pdf_unresolved_markers,
    )

    checks: list[dict] = []

    def _add(name, status, summary, detail=None):
        checks.append({
            "name": name,
            "status": status,
            "summary": summary,
            "detail": detail or {},
        })

    target_profile = _t.paper_writing_target_venue_policy(target_venue)
    target_category = target_profile.get("target_category")
    target_status = target_profile.get("status")
    _add(
        "target_venue_policy",
        (
            "pass" if target_status == "pass"
            else "fail" if target_status == "reject"
            else "warn"
        ),
        (
            f"target={target_venue or '(not selected)'}; "
            f"category={target_category or 'none'}; policy={target_status}"
        ),
        target_profile,
    )

    # v0.93.0: multi-file .tex resolution -- if tex_path uses
    # \input{...}, merge into a single file for downstream tools
    # that operate on a single .tex string.
    merged_tex_path = tex_path   # default = passed-in tex_path
    if tex_path and auto_resolve_inputs:
        try:
            r = resolve_input_chain(tex_path)
            if r.get("ok") and r.get("files_resolved"):
                n_inputs = len(r["files_resolved"]) - 1   # exclude main
                if n_inputs > 0:
                    # Write merged tex to a sibling temp file so
                    # downstream file-path-based tools can read it.
                    import tempfile
                    tf = tempfile.NamedTemporaryFile(
                        mode="w", suffix="_merged.tex", delete=False,
                        encoding="utf-8")
                    tf.write(r["merged_tex"])
                    tf.close()
                    merged_tex_path = tf.name
                    _add(
                        "multifile_resolved",
                        "pass",
                        (f"resolved {n_inputs} \\input subfiles -> "
                         f"merged to {tf.name}"),
                        {"files_resolved": r["files_resolved"],
                         "files_missing": r["files_missing"]},
                    )
        except Exception as e:  # noqa: BLE001
            _add("multifile_resolved", "skip",
                 f"resolver error: {e}")

    # v0.93.0: abstract auto-extraction
    if tex_path and not abstract_text and auto_extract_abstract:
        try:
            with open(merged_tex_path, encoding="utf-8",
                      errors="replace") as fh:
                src = fh.read()
            extracted = extract_abstract_from_tex(src)
            if extracted:
                abstract_text = extracted
                _add(
                    "abstract_extracted",
                    "pass",
                    (f"abstract auto-extracted ({len(extracted)} chars) "
                     f"from {os.path.basename(merged_tex_path)}"),
                )
        except Exception as e:  # noqa: BLE001
            _add("abstract_extracted", "skip",
                 f"extraction error: {e}")

    # ------ TEX-based checks (use merged tex if multi-file) ------
    if tex_path:
        _scan_tex = merged_tex_path
        try:
            r = _t.paper_writing_check_figure_forward_reference(_scan_tex)
            n_fwd = len(r.get("forward_refs", r.get("violations", [])))
            status = "fail" if n_fwd > 0 else "pass"
            _add("figure_forward_reference",
                 status,
                 f"{n_fwd} forward-referenced figures" if n_fwd
                 else "no forward references",
                 r)
        except Exception as e:  # noqa: BLE001
            _add("figure_forward_reference", "skip", f"tool error: {e}")

        try:
            r = _t.paper_writing_check_equation_numbering(_scan_tex)
            _add("equation_numbering", "pass", "equation numbering ok", r)
        except Exception as e:  # noqa: BLE001
            _add("equation_numbering", "skip", f"tool error: {e}")

        try:
            r = _t.paper_writing_count_underlines(_scan_tex)
            n_under = r.get("count", r.get("n_underlines", 0))
            status = "warn" if n_under > 0 else "pass"
            _add("count_underlines",
                 status,
                 f"{n_under} \\underline detected (consider \\emph)"
                 if n_under else "no underlines",
                 r)
        except Exception as e:  # noqa: BLE001
            _add("count_underlines", "skip", f"tool error: {e}")

        # v0.92.0: undefined-variable check
        try:
            r = paper_writing_check_undefined_variables(_scan_tex)
            n_undef = r.get("n_undefined", 0)
            status = "fail" if n_undef > 0 else "pass"
            _add("undefined_variables",
                 status,
                 f"{n_undef} undefined math symbols found" if n_undef
                 else "every math symbol is defined",
                 r)
        except Exception as e:  # noqa: BLE001
            _add("undefined_variables", "skip", f"tool error: {e}")

        # 2026-05-27: undefined acronyms (IH/MQS/FEM not spelled out)
        try:
            r = paper_writing_check_undefined_acronyms(_scan_tex)
            n_und = r.get("n_undefined", 0)
            status = "fail" if n_und > 0 else "pass"
            _add("undefined_acronyms",
                 status,
                 f"{n_und} acronyms not spelled out on first use" if n_und
                 else "all acronyms spelled out / in nomenclature",
                 r)
        except Exception as e:  # noqa: BLE001
            _add("undefined_acronyms", "skip", f"tool error: {e}")

        # 2026-05-27: \ref{} ↔ \label{} consistency (catches "[??]")
        try:
            r = paper_writing_check_ref_label_consistency(_scan_tex)
            st = r.get("status", "skip")
            status = "fail" if st == "fail" else ("warn" if st == "warning" else "pass")
            _add("ref_label_consistency",
                 status,
                 (f"{r.get('n_dangling_refs', 0)} dangling refs, "
                  f"{r.get('n_orphan_labels', 0)} orphan labels"),
                 r)
        except Exception as e:  # noqa: BLE001
            _add("ref_label_consistency", "skip", f"tool error: {e}")

        # Keyword blocks are venue-specific. In particular, IEEJ's
        # jkeyword/ekeyword pair must not be rejected for lacking
        # IEEEkeywords. Accelerator templates vary by venue and year.
        try:
            if target_category == "ieej":
                with open(_scan_tex, encoding="utf-8", errors="replace") as fh:
                    keyword_src = fh.read()
                ja_match = re.search(
                    r"\\begin\{jkeyword\}(.+?)\\end\{jkeyword\}",
                    keyword_src,
                    re.DOTALL,
                )
                en_match = re.search(
                    r"\\begin\{ekeyword\}(.+?)\\end\{ekeyword\}",
                    keyword_src,
                    re.DOTALL,
                )

                def _keyword_count(match, separators):
                    if not match:
                        return 0
                    items = re.split(separators, match.group(1))
                    return sum(1 for item in items if item.strip())

                ja_count = _keyword_count(ja_match, r"[，,]")
                en_count = _keyword_count(en_match, r"[,;]")
                both_present = bool(ja_match and en_match)
                enough = ja_count >= 3 and en_count >= 3
                status = "pass" if both_present and enough else "fail"
                r = {
                    "status": status,
                    "found_blocks": [
                        name for name, present in (
                            ("jkeyword", ja_match), ("ekeyword", en_match)
                        ) if present
                    ],
                    "japanese_keyword_count": ja_count,
                    "english_keyword_count": en_count,
                    "policy": "IEEJ bilingual keyword blocks",
                }
                _add(
                    "target_keywords",
                    status,
                    f"IEEJ keywords: JA {ja_count}, EN {en_count}",
                    r,
                )
            elif target_category == "accelerator":
                _add(
                    "target_keywords",
                    "skip",
                    "accelerator keyword rules are template-specific; verify the official venue template",
                    target_profile,
                )
            else:
                r = paper_writing_check_ieee_keywords(_scan_tex)
                st = r.get("status", "skip")
                status = ("fail" if st == "missing"
                          else "warn" if st == "warning"
                          else "pass")
                _add("ieee_keywords",
                      status,
                      (f"{r.get('n_keywords', 0)} keywords found "
                       f"(block: {r.get('found_block')})"),
                      r)
        except Exception as e:  # noqa: BLE001
            _add("target_keywords", "skip", f"tool error: {e}")

        # 2026-06-24: one-page digest human-review triggers
        try:
            r = _t.paper_writing_check_digest_human_review_triggers(_scan_tex)
            st = r.get("status", "skip")
            status = "warn" if st == "warning" else "pass"
            _add("digest_human_review_triggers",
                 status,
                 (f"{r.get('n_findings', 0)} digest-specific "
                  "human-review triggers"),
                 r)
        except Exception as e:  # noqa: BLE001
            _add("digest_human_review_triggers", "skip",
                 f"tool error: {e}")

        # Sentence length alone misses structurally valid but cognitively
        # dense manuscripts. Japanese and English need different units and
        # thresholds, so this check keeps both results separate.
        try:
            r = _t.paper_writing_bilingual_readability_check(_scan_tex)
            st = r.get("status", "not_applicable")
            status = (
                "fail" if st == "fail"
                else "warn" if st == "warning"
                else "pass" if st == "pass"
                else "skip"
            )
            ja = r.get("japanese", {})
            en = r.get("english", {})
            _add(
                "bilingual_adjacent_reviewer_readability",
                status,
                (
                    f"JA {ja.get('heuristic_score', 0)}/100 "
                    f"({ja.get('risk_count', 0)} risks), "
                    f"EN {en.get('heuristic_score', 0)}/100 "
                    f"({en.get('risk_count', 0)} risks)"
                ),
                r,
            )
        except Exception as e:  # noqa: BLE001
            _add(
                "bilingual_adjacent_reviewer_readability",
                "skip",
                f"tool error: {e}",
            )
    else:
        _add("tex_checks", "skip", "no tex_path supplied")

    # ------ BIB POLICY GATE (v0.91.0) ------
    # The lab policy is that EVERY citation must trace back to
    # references.bib + be verified via Crossref / S2 / arXiv.  Running
    # the gate without bib_path is a soft "fail" -- the gate emits a
    # critical warning and skips the bib-dependent checks.
    if not bib_path:
        _add(
            "bib_policy",
            "fail",
            ("references.bib was not supplied.  Lab POLICY: every "
             "citation must come from the user's actual .bib and be "
             "verified via paper_writing_verify_citation BEFORE "
             "insertion.  Re-call with bib_path=/path/to/references.bib "
             "OR explicitly justify why no .bib check applies."),
        )

    # ------ BIB-based checks ------
    if bib_path:
        try:
            r = _t.paper_writing_lint_reference_format(bib_path)
            _add("lint_reference_format", "pass",
                 "bib reference format OK", r)
        except Exception as e:  # noqa: BLE001
            _add("lint_reference_format", "skip", f"tool error: {e}")

        if tex_path:
            try:
                r = _t.paper_writing_check_citation_usage(tex_path, bib_path)
                _add("check_citation_usage", "pass",
                     "citations resolve", r)
            except Exception as e:  # noqa: BLE001
                _add("check_citation_usage", "skip", f"tool error: {e}")

            if author_last_names:
                try:
                    r = _t.paper_writing_check_self_citation_ratio(
                        tex_path, bib_path,
                        author_last_names=author_last_names,
                    )
                    _add("self_citation_ratio", "pass",
                         "self-citation ratio computed", r)
                except Exception as e:  # noqa: BLE001
                    _add("self_citation_ratio", "skip",
                         f"tool error: {e}")

            # 2026-05-27: \cite{} key ↔ .bib entry existence
            try:
                r = paper_writing_check_citation_keys_exist(
                    tex_path, bib_path)
                st = r.get("status", "skip")
                status = ("fail" if st == "fail"
                          else "warn" if st == "warning"
                          else "pass")
                _add("citation_keys_exist",
                     status,
                     (f"{r.get('n_missing_in_bib', 0)} cite keys missing in bib, "
                      f"{r.get('n_unused_in_tex', 0)} bib entries unused"),
                     r)
            except Exception as e:  # noqa: BLE001
                _add("citation_keys_exist", "skip", f"tool error: {e}")
    else:
        _add("bib_checks", "skip", "no bib_path supplied")

    # ------ Abstract-based checks ------
    if abstract_text:
        try:
            r = _t.paper_writing_validate_abstract_length(abstract_text)
            _add("validate_abstract_length", "pass",
                 "abstract length checked", r)
        except Exception as e:  # noqa: BLE001
            _add("validate_abstract_length", "skip", f"tool error: {e}")

        try:
            r = _t.paper_writing_check_abstract_background_ratio(
                abstract_text)
            _add("abstract_background_ratio", "pass",
                 "abstract background ratio checked", r)
        except Exception as e:  # noqa: BLE001
            _add("abstract_background_ratio", "skip",
                 f"tool error: {e}")

        try:
            r = _t.paper_writing_count_weak_expressions(abstract_text)
            n = r.get("total_weak_expressions", 0)
            status = "warn" if n > 0 else "pass"
            _add("abstract_weak_expressions", status,
                 f"{n} hedge expressions in abstract" if n
                 else "abstract has no weak expressions",
                 r)
        except Exception as e:  # noqa: BLE001
            _add("abstract_weak_expressions", "skip",
                 f"tool error: {e}")

        # 2026-05-27: abstract must not contain math or \cite{}
        try:
            r = _t.paper_writing_check_abstract_no_math_no_citation(
                abstract_text)
            st = r.get("status", "skip")
            status = ("fail" if st == "fail"
                      else "warn" if st == "warning"
                      else "pass")
            _add("abstract_no_math_no_citation",
                 status,
                 (f"display_math={r.get('n_math_display', 0)}, "
                  f"inline_math={r.get('n_math_inline', 0)}, "
                  f"citations={r.get('n_citations', 0)}, "
                  f"acronyms={r.get('n_acronyms', 0)}"),
                 r)
        except Exception as e:  # noqa: BLE001
            _add("abstract_no_math_no_citation", "skip",
                 f"tool error: {e}")
    else:
        _add("abstract_checks", "skip", "no abstract_text supplied")

    # ------ PDF-based checks (need both tex AND pdf for floats check) ------
    if pdf_path:
        if page_limit > 0:
            try:
                r = _t.paper_writing_validate_pdf_pages(
                    pdf_path, page_limit=page_limit)
                # Find the page-count key (varies by tool version)
                n_pages = (r.get("page_count")
                            or r.get("n_pages")
                            or r.get("count")
                            or 0)
                status = "fail" if n_pages > page_limit else "pass"
                _add("validate_pdf_pages", status,
                     f"{n_pages} pages (limit {page_limit})", r)
            except Exception as e:  # noqa: BLE001
                _add("validate_pdf_pages", "skip", f"tool error: {e}")

        try:
            r = paper_writing_detect_page_whitespace_anomalies(
                pdf_path, whitespace_threshold=whitespace_threshold)
            n_flag = r.get("flagged_count", 0)
            status = "warn" if n_flag > 0 else "pass"
            _add("page_whitespace_anomalies", status,
                 f"{n_flag} mostly-white pages flagged" if n_flag
                 else "no mostly-white pages",
                 r)
        except Exception as e:  # noqa: BLE001
            _add("page_whitespace_anomalies", "skip",
                 f"tool error: {e}")

        if tex_path:
            try:
                r = paper_writing_check_floats_far_from_reference(
                    tex_path, pdf_path,
                    max_pages_apart=layout_max_pages_apart)
                n_flag = r.get("flagged_count", 0)
                status = "warn" if n_flag > 0 else "pass"
                _add("floats_far_from_reference", status,
                     f"{n_flag} floats drift > {layout_max_pages_apart}"
                     f" pages from their refs" if n_flag
                     else "floats stay near their refs",
                     r)
            except Exception as e:  # noqa: BLE001
                _add("floats_far_from_reference", "skip",
                     f"tool error: {e}")

        # v0.92.0: pixel-accurate overlap checks
        try:
            r = paper_writing_detect_text_image_overlap(pdf_path)
            n_ovl = r.get("n_overlaps", 0)
            status = "fail" if n_ovl > 0 else "pass"
            _add("text_image_overlap", status,
                 f"{n_ovl} text-on-image overlaps detected" if n_ovl
                 else "no text-image overlap",
                 r)
        except Exception as e:  # noqa: BLE001
            _add("text_image_overlap", "skip", f"tool error: {e}")

        try:
            r = paper_writing_detect_text_overflow_page(pdf_path)
            n_ofl = r.get("n_overflows", 0)
            status = "fail" if n_ofl > 0 else "pass"
            _add("text_overflow_page", status,
                 f"{n_ofl} text blocks overflow past page edge" if n_ofl
                 else "no text overflow",
                 r)
        except Exception as e:  # noqa: BLE001
            _add("text_overflow_page", "skip", f"tool error: {e}")

        # 2026-05-27: post-compile [?] / [??] marker scan
        try:
            r = paper_writing_check_pdf_unresolved_markers(pdf_path)
            st = r.get("status", "skip")
            status = "fail" if st == "fail" else "pass"
            _add("pdf_unresolved_markers",
                 status,
                 (f"{r.get('n_markers_total', 0)} '[?]'/'[??]' markers "
                  f"in rendered PDF"),
                 r)
        except Exception as e:  # noqa: BLE001
            _add("pdf_unresolved_markers", "skip", f"tool error: {e}")
    else:
        _add("pdf_checks", "skip", "no pdf_path supplied")

    # ------ Verdict aggregation ------
    n_critical = sum(1 for c in checks if c["status"] == "fail")
    n_warning = sum(1 for c in checks if c["status"] == "warn")
    n_passed = sum(1 for c in checks if c["status"] == "pass")
    n_skipped = sum(1 for c in checks if c["status"] == "skip")

    if n_critical > 0:
        verdict = "fail"
        advice = (f"BLOCKING: {n_critical} critical issue(s) must be "
                  f"fixed before submission.  Fix the 'fail' rows then "
                  f"re-run paper_writing_em_submission_gate.")
    elif n_warning > 0:
        verdict = "warn"
        advice = (f"REVIEWABLE: {n_warning} warning(s).  Each may be "
                  f"acceptable depending on context; review the 'warn' "
                  f"rows and silence by editing OR by passing tighter "
                  f"thresholds (whitespace_threshold, "
                  f"layout_max_pages_apart) if the defaults are wrong "
                  f"for your venue.")
    else:
        verdict = "pass"
        advice = (f"SUBMISSION-READY: all {n_passed} checks passed "
                  f"({n_skipped} skipped due to missing inputs).  "
                  f"Also run the EM-style knowledge checklist via "
                  f"paper_writing_em_paper_style('reviewer_patterns').")

    return {
        "verdict": verdict,
        "target_venue": target_venue,
        "target_category": target_category,
        "n_checks_run": len(checks),
        "n_critical": n_critical,
        "n_warning": n_warning,
        "n_passed": n_passed,
        "n_skipped": n_skipped,
        "checks": checks,
        "advice": advice,
    }
