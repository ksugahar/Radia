"""LaTeX figure placement knowledge for radia_mcp.paper_writing.

Added 2026-05-26 (radia-mcp v0.88.0).  Companion to the existing
text-style tools (figure_forward_reference, citation_usage); this
module captures the LAYOUT side of figure handling that the lab's
papers historically get wrong.

Web-research basis (2026-05-26):
  * LaTeX/Floats, Figures and Captions (Wikibooks): the canonical
    open reference.
  * Hyndman R.J., "Controlling figure and table placement in LaTeX",
    blog 2014.
  * Overleaf docs: "!h float specifier changed to !ht" -- the
    famous LaTeX silent rewrite that catches many authors.
  * LaTeX Cloud Studio docs: positioning [h], [t], [b], [p], [H]
    semantics + placeins/FloatBarrier strategy.
  * IEEE journal LaTeX templates: \\figure*{} vs \\figure{} +
    \\columnwidth vs \\textwidth conventions for two-column papers.

What this module is NOT: it does NOT include the lab's existing
text-based figure checks (forward references, caption density)
which already live in paper_writing.tools as
``paper_writing_check_figure_forward_reference`` and
``paper_writing_check_figure_caption_showing``.  This module is the
PLACEMENT-and-WIDTH knowledge layer those tools couldn't carry.
"""

from __future__ import annotations


OVERVIEW = r"""
# LaTeX figure placement -- lab quick reference

The single most common reviewer-visible defect in lab papers is
figures placed on the wrong page (or, worse, on a wasteful
mostly-blank page).  The fix is almost always to change the float
specifier ``[...]`` argument.  This file is the cookbook.

## TL;DR

  * USE ``[htbp]`` (or ``[!htbp]`` for last resort) as the DEFAULT.
  * AVOID ``[h]`` alone (LaTeX silently rewrites it to ``[ht]``).
  * AVOID ``[H]`` (from the ``float`` package) unless you genuinely
    need NON-FLOATING -- in 90% of cases what you actually want is
    ``[!ht]`` plus a ``\\FloatBarrier``.
  * NEVER place a figure FORWARD of its first ``\\ref{fig:X}``;
    LaTeX cannot pull a float backward.
  * WIDTH: use ``\\columnwidth`` inside ``figure``, use
    ``\\textwidth`` inside ``figure*`` (two-column spanning).
  * CHECK the rendered PDF with the image-based recipe (see
    ``paper_writing_layout_visual_recipe``).

## Why this matters (lab pattern)

Reviewer comment, anonymous IEEE TAP 2025 (lab paper):
  "Figure 4 is referenced on page 3 but appears on page 6, after
   the conclusions.  Please fix."

Root cause: the author wrote ``[h!]`` thinking it would force
LaTeX to put the figure "here".  LaTeX's actual behaviour: ``h``
means "try here", ``!`` waives the float-fraction constraint, but
neither overrides the running-text balance algorithm.  Result:
LaTeX bumped the float to the next valid spot it could find, which
happened to be after the bibliography.

The correct construct was ``[!t]`` + earlier source placement OR
``[H]`` (from ``\\usepackage{float}``).
"""


FLOAT_SPECIFIERS = r"""
# Float specifiers: what each letter means

LaTeX's ``[htbp]`` letters in any order tell LaTeX WHERE it's
ALLOWED to put the float (not commanded -- ALLOWED).  The order in
the brackets does NOT matter; the order of attempts is hardwired
in the LaTeX kernel: first ``h`` (only if requested), then ``t``,
then ``b``, then ``p``.

| Letter | Meaning | When LaTeX uses it |
|--------|---------|--------------------|
| ``h`` | here, in place | only IF the float fits in remaining space on current page |
| ``t`` | top of next page | rule-of-thumb default for IEEE/IEEJ |
| ``b`` | bottom of current/next page | works well for small figures |
| ``p`` | float-only page | LAST RESORT; collects "pX" deferred floats |
| ``H`` | here, FORCED, non-float | requires ``\\usepackage{float}``; treats fig as inline image |

## The ``!`` modifier

``!`` waives the float-fraction constraints (\\floatpagefraction,
\\textfraction, \\topfraction, \\bottomfraction).  These are
internal LaTeX limits like "no more than 70% of a text page can be
floats".  ``!`` says "ignore those limits and place the float
anyway".

  * ``[!h]`` -- LaTeX SILENTLY rewrites this to ``[!ht]`` since
    LaTeX2e (Overleaf-documented behaviour).  If you really mean
    "here, no flexibility, no float", use ``[H]`` from the float
    package.
  * ``[!htbp]`` -- the maximum-flexibility default.

## Lab default recommendation

For every figure UNLESS you have a specific reason otherwise:

```latex
\\begin{figure}[!htbp]
    \\centering
    \\includegraphics[width=\\columnwidth]{fig.pdf}
    \\caption{...}
    \\label{fig:foo}
\\end{figure}
```

For two-column papers (IEEE Transactions, IEEJ Trans):
  * Single-column figure  -> ``\\begin{figure}[!t]`` + ``width=\\columnwidth``
  * Two-column-span fig   -> ``\\begin{figure*}[!t]`` + ``width=\\textwidth``
  * Top of NEXT page is the lab default because two-column papers
    really only fit floats at top of page.

## Why ``[H]`` is usually wrong

``[H]`` (from ``\\usepackage{float}``) treats the figure as an
inline element with NO floating.  This guarantees position but at a
huge cost:

  * Can leave gigantic blank space below the figure if it doesn't
    fit on the current page.
  * Forces LaTeX into awkward page breaks.
  * Defeats LaTeX's whole text-flow optimisation.

Use ``[H]`` ONLY when the figure MUST be at the exact source
position -- e.g. inside a step-by-step proof, or as a side-by-side
comparison that must stay aligned with specific paragraphs.  For
normal papers, ``[!ht]`` + ``\\FloatBarrier`` does the same job
without the blank space.
"""


PLACEINS_FLOATBARRIER = r"""
# placeins.sty + \\FloatBarrier -- the lab tactic for keeping floats
# in their section

The ``placeins`` package gives you ``\\FloatBarrier`` -- a barrier
that LaTeX floats CANNOT cross.

```latex
\\usepackage{placeins}
% ... later ...
\\section{Results}
... text and figures ...
\\FloatBarrier
\\section{Discussion}
```

The barrier forces any unplaced floats to appear BEFORE the next
section starts.  Without it, a figure from the Results section can
drift into the Discussion (or further), separating it from the
text that introduced it.

## ``[section]`` option -- the one-line fix

```latex
\\usepackage[section]{placeins}
```

This option redefines ``\\section`` to AUTOMATICALLY include a
``\\FloatBarrier``.  Add this single line to your preamble and you
never have to manually scatter ``\\FloatBarrier`` commands.

**Lab default**: always use ``\\usepackage[section]{placeins}``
unless the journal template forbids it.  No downside; large benefit.

## Sub-section barrier

```latex
\\usepackage[subsection]{placeins}
```

Barrier at every ``\\subsection`` too.  Use this for review
articles or thesis chapters where each subsection has its own
figures.  For standard 6-page conference papers, ``[section]`` is
enough.

## Manual barriers

For finer control (e.g. inside a long Results section):

```latex
\\subsection{Validation against analytical solution}
... text ...
\\begin{figure}[!ht]
    \\includegraphics{validation.pdf}
    \\caption{...}
\\end{figure}
... text ...
\\FloatBarrier      % force the figure above to settle here
\\subsection{Parametric study}
```

## Common pitfalls

  * ``\\FloatBarrier`` inside a two-column ``figure*`` environment
    is a no-op (the spanning float plays by different rules).
  * If a figure is too tall to fit before the barrier, LaTeX will
    complain ``! LaTeX Error: Float(s) lost.`` -- the figure won't
    appear.  Fix by making the figure shorter (or by removing the
    barrier and letting the float drift).
"""


WIDTH_FRACTIONS = r"""
# Figure widths: \\columnwidth, \\textwidth, \\linewidth

In a two-column journal layout (IEEE Trans, IEEJ Trans, IGTE
digest), the document has THREE conceptual widths:

  * ``\\columnwidth`` -- width of ONE column (~88.9 mm = 3.5 in for IEEE).
  * ``\\textwidth``   -- width of BOTH columns together (~181 mm).
  * ``\\linewidth``   -- equal to ``\\columnwidth`` in 2-col body,
                          equal to ``\\textwidth`` in single-col body.

## The figure / figure* distinction

```latex
\\begin{figure}[!t]
    \\centering
    \\includegraphics[width=\\columnwidth]{onecol.pdf}
    \\caption{Single-column figure.}
\\end{figure}

\\begin{figure*}[!t]
    \\centering
    \\includegraphics[width=\\textwidth]{twocol.pdf}
    \\caption{Two-column-spanning figure.}
\\end{figure*}
```

The ``*`` form spans both columns.  In two-column journals, ``*``
floats can ONLY appear at the TOP of a page (``[!t]`` is the only
sensible specifier).

## When to span both columns

  * Wide plots that lose too much information at column width.
  * Multi-panel figures where each panel is ~column-wide.
  * Schematic diagrams of complex apparatus.

For most line graphs and convergence plots: stay single-column,
use ``\\columnwidth``.  The lab default is single-column because
two-column-span eats vertical space disproportionately.

## What NOT to do

  * ``width=8cm`` -- ABSOLUTE width breaks when the journal
    template changes.  Always use relative width
    (``\\columnwidth``, ``0.5\\textwidth``, etc.).
  * ``width=1.0\\columnwidth`` and ``\\linewidth`` MIXED -- pick
    one and use it everywhere.  ``\\linewidth`` is safer because it
    auto-adapts inside subfigures and minipages.
  * Forgetting ``\\centering`` -- figure gets left-aligned;
    looks wrong in tight columns.
  * Forgetting the trailing units in ``width=`` -- ``width=8``
    silently becomes ``width=8pt`` (4 px wide image, often
    invisible in the PDF).

## Subfigure widths

Use the ``subcaption`` package (NOT the older ``subfigure``):

```latex
\\usepackage{subcaption}
% ...
\\begin{figure}[!t]
    \\centering
    \\begin{subfigure}{0.48\\columnwidth}
        \\includegraphics[width=\\linewidth]{left.pdf}
        \\caption{Left panel.}
        \\label{fig:foo-left}
    \\end{subfigure}
    \\hfill
    \\begin{subfigure}{0.48\\columnwidth}
        \\includegraphics[width=\\linewidth]{right.pdf}
        \\caption{Right panel.}
        \\label{fig:foo-right}
    \\end{subfigure}
    \\caption{Two-panel comparison.}
    \\label{fig:foo}
\\end{figure}
```

Note:
  * Sum of subfigure widths + ``\\hfill`` should NOT exceed
    ``\\columnwidth`` (or ``\\textwidth`` for ``figure*``).  Use
    ``0.48`` each + ``\\hfill`` to leave 4% gap between panels.
  * ``\\includegraphics[width=\\linewidth]`` inside a subfigure
    uses the SUBFIGURE width, not the figure width.  That's why
    ``\\linewidth`` is safer than ``\\columnwidth`` here.
  * Each subfigure should have its OWN ``\\caption`` AND
    ``\\label``; the outer ``\\caption`` is the figure-level title.
  * Refer to panels as ``Fig.~\\ref{fig:foo}(a)`` or
    ``Fig.~\\subref{fig:foo-left}`` (with the right
    ``subcaption`` options).

## Cross-reference with radia_mcp.figure

The lab's matplotlib ``paper_figure_8cm(panels=1)`` already
produces a single-column-friendly figure at ~5 cm axes.  When you
``\\includegraphics{out.pdf}`` the result, use ``width=\\columnwidth``
in the figure environment so the absolute-size font scaling stays
correct.  See ``mcp-server-figure`` ``paper_figure_8cm_recipe``.
"""


ANTIPATTERNS = r"""
# Anti-patterns -- what NOT to do

## 1. The "h!" trap

```latex
\\begin{figure}[h!]       % WRONG: LaTeX silently rewrites to [!ht]
    ...
\\end{figure}
```

The ``!`` belongs BEFORE the letters, not after:
  * ``[!ht]`` -- correct (force flex, allow here-or-top)
  * ``[h!]``  -- silently becomes ``[!ht]`` per LaTeX2e kernel

If you really wanted "h forced", you wanted ``[H]`` from the
``float`` package.

## 2. ``\\centering`` vs ``\\begin{center}``

```latex
\\begin{figure}
    \\begin{center}        % WRONG: adds extra vertical space
        \\includegraphics{...}
    \\end{center}
    \\caption{...}
\\end{figure}
```

The ``center`` environment adds vertical spacing above and below
that breaks the figure's vertical box.  Always use ``\\centering``:

```latex
\\begin{figure}
    \\centering            % CORRECT
    \\includegraphics{...}
    \\caption{...}
\\end{figure}
```

## 3. Caption before \\includegraphics for figures

```latex
\\begin{figure}
    \\caption{...}         % WRONG for figure (right for table)
    \\includegraphics{...}
\\end{figure}
```

Convention:
  * ``figure``: caption AFTER ``\\includegraphics``.
  * ``table``:  caption BEFORE ``\\tabular``.

Putting the caption in the wrong place doesn't break compilation
but breaks lab style and looks unfamiliar to reviewers.

## 4. Forgetting \\label after \\caption

```latex
\\begin{figure}
    \\centering
    \\includegraphics{...}
    \\caption{...}
    \\label{fig:foo}       % must come AFTER \\caption
\\end{figure}
```

Order matters: ``\\label`` reads the LAST ``\\refstepcounter`` to
know which counter to label.  Before ``\\caption``, the figure
counter hasn't been stepped yet, so ``\\ref{fig:foo}`` resolves to
the previous figure's number (off-by-one bug).

## 5. Blank line inside the float environment

```latex
\\begin{figure}[!ht]

    \\centering            % WRONG: blank line before this
    \\includegraphics{...}
    \\caption{...}
\\end{figure}
```

A blank line inside ``figure`` is treated as a paragraph break,
which can confuse LaTeX's vertical alignment.  Symptom: caption
shifts downward by one ``\\baselineskip``.  Always keep the
``figure`` body contiguous.

## 6. Mixing absolute and relative widths in one paper

```latex
% Fig. 3
\\includegraphics[width=\\columnwidth]{fig3.pdf}
% Fig. 4
\\includegraphics[width=8cm]{fig4.pdf}        % WRONG: 8cm
```

Use relative width EVERYWHERE in the same paper.  Otherwise,
switching journals (one-column to two-column, or letter to A4)
makes some figures the right size and others jarringly wrong.

## 7. \\FloatBarrier inside figure*

```latex
\\begin{figure*}[!t]
    ...
    \\FloatBarrier         % WRONG: no-op for two-column-spanning floats
\\end{figure*}
```

Two-column floats use the ``\\dbltop`` queue, which placeins does
NOT manage.  To barrier a ``figure*``, use ``\\afterpage{\\clearpage}``
from the ``afterpage`` package or just put the float earlier in
source.

## 8. Forward-referencing a figure

```latex
... see Fig.~\\ref{fig:results} for details ...   % page 3

\\begin{figure}[!t]
    \\caption{Results.}
    \\label{fig:results}                          % page 7 (after [t] dump)
\\end{figure}
```

LaTeX CANNOT pull a float backward.  If ``\\begin{figure}`` is on
source line 1000 and ``\\ref{fig:results}`` is on line 500, the
figure will appear AFTER line 500 in the output -- never before.

Fix: move ``\\begin{figure}`` to BEFORE the first reference.  The
mechanical rule: ``\\begin{figure}`` source order should match
intended figure order.

This is also why ``paper_writing_check_figure_forward_reference``
exists -- it catches this pattern automatically.
"""


JOURNAL_PROFILES = r"""
# Journal-specific placement profiles

## IEEE Transactions (IEEEtran.cls, two-column)

```latex
\\documentclass[journal,onecolumn]{IEEEtran}
% OR (default):
\\documentclass[journal]{IEEEtran}  % two-column

\\usepackage[section]{placeins}    % LAB DEFAULT
\\usepackage{subcaption}
\\usepackage{graphicx}
```

  * Single-column figure: ``\\begin{figure}[!t]`` + ``width=\\columnwidth``
  * Two-column-spanning:  ``\\begin{figure*}[!t]`` + ``width=\\textwidth``
  * Top placement ``[!t]`` is strongly preferred; IEEE typesetting
    rejects mid-column figures in production.

## IEEJ Transactions (jiee.cls, two-column)

Same conventions as IEEE.  IEEJ explicitly uses ``\\columnwidth``
= 88 mm in their template; absolute-sized figures will break the
production pipeline.

## IGTE / Compumag digest (two-column A4 wide)

```latex
\\documentclass{igte_digest}     % template

\\begin{figure}[!htbp]            % digest tolerates inline floats
    \\centering
    \\includegraphics[width=\\columnwidth]{...}
    ...
\\end{figure}
```

Digest papers are 4-6 pages; ``[!htbp]`` is preferred over ``[!t]``
because the brevity makes per-page balance more flexible.

## Thesis / book class (single column, scrbook/memoir)

```latex
\\usepackage[section]{placeins}    % barrier at every chapter too
\\begin{figure}[!htbp]
    \\centering
    \\includegraphics[width=0.7\\textwidth]{...}
    ...
\\end{figure}
```

Single-column thesis: use ``width=0.7\\textwidth`` (or ``0.85``)
to leave margin space; full ``\\textwidth`` overwhelms a single
column.

## Submission / review version with endfloat

```latex
\\usepackage{endfloat}
% body: figures and tables stay where you wrote them in source
% PDF output: all floats are collected at the END of the document
```

Some journals (e.g. older biology / physics) REQUIRE figures
collected at the end for review.  ``\\usepackage{endfloat}`` does
this transparently.  Remove for camera-ready.
"""


CROSS_REFERENCES = r"""
# Cross-references to other paper_writing tools

This module is the KNOWLEDGE side of figure placement.  The
ENFORCEMENT side already lives in paper_writing.tools.  Both work
together:

  * ``paper_writing_check_figure_forward_reference(tex_path)``
        -> textual check: \\ref{fig:X} BEFORE \\begin{figure}.
    Symptom of bad placement source order.

  * ``paper_writing_check_pdf_edge_overflow(pdf_path)``
        -> PDF binary check: text or image overflows printable area.
    Often caused by ``width=\\textwidth`` inside ``figure`` (should
    be ``\\columnwidth``).

  * ``paper_writing_check_overfull_hbox(log_path)``
        -> parse .log for "Overfull \\hbox" warnings.  These almost
    always come from too-wide \\includegraphics or too-long
    \\caption.

  * ``paper_writing_layout_thumbnail_strip(pdf_path)``
        -> visual: render every page as a thumbnail tile.  Bad
    figure placement = obvious at a glance (wasteful white pages,
    figure-on-wrong-page, etc.).

  * ``paper_writing_detect_page_whitespace_anomalies(pdf_path)``
        -> algorithmic: flag pages > 75% white.  These are usually
    [!h] or [H] forcing a figure to its own page.

  * ``paper_writing_check_floats_far_from_reference(tex, pdf)``
        -> heuristic: figure label refs spanning > 1 page indicates
    the float drifted from its first reference.

## Recommended workflow (writing -> compile -> verify)

```
1. WRITE   -> follow this module's defaults: [!htbp] + \\columnwidth +
              \\centering + \\label after \\caption
2. COMPILE -> pdflatex paper.tex
3. CHECK   -> paper_writing_check_overfull_hbox(paper.log)
           -> paper_writing_check_figure_forward_reference(paper.tex)
           -> paper_writing_layout_thumbnail_strip(paper.pdf)
           -> paper_writing_detect_page_whitespace_anomalies(paper.pdf)
4. FIX     -> based on flagged issues, edit float specifier OR move
              \\begin{figure} earlier in source
5. RE-COMPILE + RE-CHECK until clean.
```
"""


TOPICS = {
    "overview":               OVERVIEW,
    "specifiers":             FLOAT_SPECIFIERS,
    "float_specifiers":       FLOAT_SPECIFIERS,
    "htbp":                   FLOAT_SPECIFIERS,
    "placeins":               PLACEINS_FLOATBARRIER,
    "floatbarrier":           PLACEINS_FLOATBARRIER,
    "widths":                 WIDTH_FRACTIONS,
    "columnwidth":            WIDTH_FRACTIONS,
    "subcaption":             WIDTH_FRACTIONS,
    "subfigure":              WIDTH_FRACTIONS,
    "antipatterns":           ANTIPATTERNS,
    "bad_patterns":           ANTIPATTERNS,
    "h_exclamation_trap":     ANTIPATTERNS,
    "journals":               JOURNAL_PROFILES,
    "journal_profiles":       JOURNAL_PROFILES,
    "ieee":                   JOURNAL_PROFILES,
    "ieej":                   JOURNAL_PROFILES,
    "igte":                   JOURNAL_PROFILES,
    "cross_refs":             CROSS_REFERENCES,
    "workflow":               CROSS_REFERENCES,
    "related_tools":          CROSS_REFERENCES,
}


def paper_writing_tex_figure_placement(topic: str = "overview") -> str:
    """LaTeX figure placement knowledge: float specifiers, placeins,
    widths, anti-patterns, per-journal profiles, related tools.

    The KNOWLEDGE-side companion to paper_writing.tools' existing
    text-based figure checks (forward_reference, edge_overflow,
    overfull_hbox).  Captures the placement-rule expertise that
    those check tools enforce.

    Args:
        topic: one of TOPICS keys.  Aliases accepted.

    Available topics:
      overview                 -- TL;DR + lab pattern + why this matters
      specifiers / htbp        -- h/t/b/p/H + ! modifier semantics
      placeins / floatbarrier  -- placeins package + [section] option
      widths / columnwidth / subcaption
                                -- columnwidth/textwidth, figure vs figure*,
                                   subfigure widths
      antipatterns             -- 8 common mistakes incl. [h!] trap,
                                   center vs centering, blank line in figure,
                                   forward reference, etc.
      journals / ieee / ieej / igte
                                -- per-journal defaults (IEEEtran, jiee,
                                   igte_digest, scrbook/memoir, endfloat)
      cross_refs / workflow    -- which paper_writing tool to call when,
                                   and the write-compile-check-fix loop
    """
    key = topic.strip().lower()
    if key == "all":
        return "\n\n---\n\n".join([
            OVERVIEW,
            FLOAT_SPECIFIERS,
            PLACEINS_FLOATBARRIER,
            WIDTH_FRACTIONS,
            ANTIPATTERNS,
            JOURNAL_PROFILES,
            CROSS_REFERENCES,
        ])
    if key not in TOPICS:
        return (f"Unknown topic '{topic}'. Available: "
                f"{', '.join(sorted(set(TOPICS)))}, or 'all'.")
    return TOPICS[key]
