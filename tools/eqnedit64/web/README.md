# Eqnedit Web

This directory is the canonical source of the browser/JavaScript edition of
the equation editor. The native application and Web edition are maintained in
the same Radia change and follow [`../docs/PRODUCT_PARITY.md`](../docs/PRODUCT_PARITY.md).

The Sugahara Laboratory homepage is the current publication surface, not the
source repository. Its site builder reads this directory from the Radia
checkout (normally `S:\Radia\01_GitHub`, or the checkout named by the
`RADIA_REPOSITORY` environment variable), expands
`equation-editor.fragment.html` into the teaching page, copies
`equation-editor.js` into the generated site, and verifies the copied
JavaScript with SHA-256. Do not retain or edit an independent homepage source
copy.

The homepage release gate is intentionally isolated from the other 3D teaching
material: run `site_builder/tools/run_eqnedit64_release_qa.ps1` in the homepage
workspace. It builds only the editor page and asset, then runs the two-viewport
browser contract and hidden PowerPoint native-equation render test. An
Eqnedit64-only change does not require the Mathematica, NGSolve, or all-page 3D
curriculum suite.

The host page must load MathJax 3 and include the markup in
`equation-editor.fragment.html`. The script is deliberately dependency-free
apart from that host-provided MathJax runtime.

The TeX source pane and the most-recent insertion display share the
`--eqed-source-font` CSS variable. Its browser fallback stack includes a
monospace system face and Japanese-capable UI faces; do not hard-code Consolas
or assume that one locally installed font resolves correctly. This is the Web
counterpart of the native source-legibility contract, not a JavaScript port of
the Win32 physical-face, cmap, and raster-ink probes.

`Tab` and `Shift+Tab` move to the next and previous empty `{}` slot, and
`Enter` is the same structural row break the native editor performs in its TeX
source pane: it writes the row separator `\\`, keeps the alignment column by
starting the new row with `&`, and wraps a bare expression in `aligned` when no
row environment is open yet. It never splits a group such as `\frac{}{}` or
`\text{}`; inside one it opens an empty row after the current one instead.
`Shift+Enter` keeps the plain source newline, and an IME confirmation `Enter`
(`isComposing`, key code 229) is left to the input method. Inserted environment
templates are laid out with the native source rules — break after `\begin`,
before `\end`, and after each `\\`, indented by environment depth — and their
cells are `{}` slots so `Tab` reaches every one of them.

Every insertion and row break is applied through the browser's own editing
command, so `Ctrl+Z` undoes a palette click exactly like typed text. Assigning
`textarea.value` would discard the undo history; do not reintroduce it.

The hint line under the source pane and the recent-insertion display show bare
TeX on purpose and are marked `data-tex-literal-ok` so the homepage rendering QA
does not read them as unrendered math.

The `R x` / `I x` / `B x` math-alphabet group is always visible beside the
category tabs. It inserts `\mathrm{}`, `\mathit{}`, or `\mathbf{}`; a source
selection is wrapped in the chosen command and an empty selection leaves the
caret inside the new braces.
Less common alphabets stay in Decoration: `\mathsf`, `\mathtt`, `\mathcal`,
`\mathbb`, `\mathfrak`, `\bm`, and the `\mathnormal` reset. The Web source
keeps `\bm`; only the MathJax boundary expands it to `\boldsymbol`.

Inputs and saved/copied source are TeX. MTEF and `.eqn` are not supported
formats. Office copy emits inline MathML in `text/html` plus plain TeX; image
copy is a separate action so PowerPoint cannot accidentally prefer PNG over an
editable equation.
