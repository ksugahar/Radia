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
  `mcp-server-graph.paper_figure_quality_rules('font_rule')`).
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

  reviewer_patterns / reviewer_comments
        12 common EM-paper reviewer-comment patterns (sign
        inconsistency, B-vs-H, mesh independence, CPU/memory,
        contribution clarity, missing recent citations, etc.) and
        the pre-check tool that catches each.

  all
        Concatenate all 6 sections (~25 KB total).

Use as a CHECKLIST before submission, alongside the algorithmic checks
in tools.py and the new v0.88-0.89 layout / arxiv tools.
"""


def paper_writing_em_paper_style(topic: str = "overview") -> str:
    """EM-domain-specific style/notation/convention knowledge.

    Companion to the generic checks in paper_writing.tools.  Captures
    the rules that EM reviewers (IEEE TAP/TMag/TMTT, IEEJ Trans D/B,
    IGTE) actually enforce, plus the 12-pattern catalogue of common
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

    checks: list[dict] = []

    def _add(name, status, summary, detail=None):
        checks.append({
            "name": name,
            "status": status,
            "summary": summary,
            "detail": detail or {},
        })

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
    else:
        _add("tex_checks", "skip", "no tex_path supplied")

    # ------ BIB POLICY GATE (v0.91.0) ------
    # The lab policy is that EVERY citation must trace back to
    # reference.bib + be verified via Crossref / S2 / arXiv.  Running
    # the gate without bib_path is a soft "fail" -- the gate emits a
    # critical warning and skips the bib-dependent checks.
    if not bib_path:
        _add(
            "bib_policy",
            "fail",
            ("reference.bib was not supplied.  Lab POLICY: every "
             "citation must come from the user's actual .bib and be "
             "verified via paper_writing_verify_citation BEFORE "
             "insertion.  Re-call with bib_path=/path/to/reference.bib "
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
        "n_checks_run": len(checks),
        "n_critical": n_critical,
        "n_warning": n_warning,
        "n_passed": n_passed,
        "n_skipped": n_skipped,
        "checks": checks,
        "advice": advice,
    }
