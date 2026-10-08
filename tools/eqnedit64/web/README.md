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

A formal Eqnedit64 release includes both the native distribution and this Web
publication. After GitHub/PyPI publication, build the homepage from the same
release tag, update its versionless Windows download link, publish the editor
page and changed assets, and verify public JS SHA-256 and browser behavior.
The release is incomplete until this homepage phase passes. The complete
procedure is in `.agents/skills/release-eqnedit64/SKILL.md` at the repository root.
The editor lives at `elemag/equation-editor.php`; the learning portal links to it.

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
formats. Office copy candidate emits conditional RichEdit HTML OMML in `text/html`
plus plain TeX and a non-Office MathML fallback. PNG copy remains a separate
action. Native and Web need not use identical transports: UXP-0029 prioritizes
editable limits and row structure. Current production Web MathML was observed
to flatten bounds and rows; the OMML candidate requires the user's PowerPoint
H5/H6 re-test before 3.1.2 is released. Browser payload tests are not Office
acceptance. Candidate BUILD is `3.1.2 (2026-10-08 OMML candidate)`.

The candidate exporter handles token runs, normal/italic/bold/bold-italic,
colour inheritance, Unicode mathematical alphanumerics (including script,
fraktur, double-struck, sans-serif and monospace), fractions and binomials, roots, sub/sup scripts and
large-operator limits, accents/bars, arrays/matrices, phantom and box/strike
enclosures. Unsupported MathML elements, dimensional padding (e.g. smash),
other mathvariant values and other enclosures stop copy with a diagnostic;
they are not flattened silently. Preview and saved TeX remain available.
This support boundary must be reviewed before final acceptance.

The OMML candidate approximates each nonzero mspace with one NBSP; quad,
thin and negative spacing therefore do not preserve exact widths. The user's
H5/H6 re-test must also reject duplicated fallback content or visible
conditional markup. Unicode alphabet payloads are browser-tested, not a
claim of independently measured PowerPoint font rendering.


## Shared Office clipboard contract (UXP-0031, 2026-10-08)

This section supersedes earlier candidate route descriptions above. The user
reported H5/H6 Web OK on 8361ec916. Both editions now publish conditional
RichEdit-HTML OMML through CF_HTML as the primary Office format. Native copy
publishes neither registered `MathML` nor `MathML Presentation`, so PowerPoint
cannot prefer an untested competing registered equation format. Raw TeX, Office
TeX, EMF and DIBV5 remain available in the native clipboard.

`eqnedit64.office-omml-bits.v1` compares the literal UTF-8 bytes inside the
`<!--[if gte msEquation 12]>...<![endif]-->` branch. It excludes CF_HTML byte
offsets, HTML document envelopes, SourceURL and the non-Office MathML fallback.
It does not normalize the compared bytes. Producers canonicalize named colours
to uppercase hex and merge adjacent equal-style runs before publication; text,
namespace, scripts, colour, styles and row structure remain in the comparison.
The fixed corpus covers H1/H5/H6, cases, pmatrix, colour and math alphabets.

The mandatory Windows CI gate reads the actual native and browser clipboard,
rejects registered MathML formats, and checks both wire payloads against their
producers. Missing conditional OMML or any byte difference fails the gate.
`test_office_bit_parity.cjs` producer-only mode is useful locally but cannot
qualify the shared route. Evidence records source SHA, Web/corpus hashes, mode,
per-item byte lengths and SHA-256 values. This gate proves corpus payload identity,
not Office rendering or equivalence of all accepted TeX.

The changed EXE route was tested on signed ea787664d (H1/H5/H6 OK; H7
colour loss accepted for 3.1.2).
After this qualification, a passing actual-wire gate on the release candidate
plus an EXE hand test permits omitting the separate Web hand test for covered
fixtures. A missing/failing gate, a changed shared transport contract or different
clipboard/browser/Office environment requires renewed Web qualification.
Historical run 37692902458 on ea787664d used headless Chromium: the Web side provided in-memory browser-clipboard (producer-equivalent) evidence, not Web OS-clipboard evidence. The first genuine Web OS-clipboard evidence is 3.1.3 run 37718032042 on b06b8d89a (headed full Chromium).

The wire gate runs headed full Chromium. On the hosted Windows runner, both old-headless and new-headless modes do not reach the OS clipboard. The harness independently proves the Web write reached Windows before comparing captured payloads. Hosted run 37718032042 on b06b8d89a established real Windows clipboard identity for both editions: 36/36 covered fixtures; boxed/cancel/binom remain non-covered.
User H1/H5/H6 EXE passed; H7 colour loss is accepted for 3.1.2 as documented below.


## 3.1.2 PowerPoint colour limitation (UXP-0030/0031)

PowerPoint import drops colour on the shared conditional OMML CF_HTML route in
both editions. The user measured all r10 colour variants A–G as not red on
2026-10-08 and approved releasing 3.1.2 with this documented limitation. Colour
support remains enabled in the native/Web canvas, SVG, saved TeX, and native
EMF/DIB output. This decision does not strip colour from the source or renderer.
PowerPoint colour retention is not an acceptance requirement for this release.

User hand tests: H1/H5/H6 EXE OK on ea787664d; H7 EXE canvas red, PowerPoint
not red. Web H5/H6 OK on 8361ec916. Historical run 37692902458 on ea787664d used headless Chromium: the Web side provided in-memory browser-clipboard (producer-equivalent) evidence, not Web OS-clipboard evidence. The first genuine Web OS-clipboard evidence is 3.1.3 run 37718032042 on b06b8d89a (headed full Chromium).
The user hand tests remain valid; the headed gate establishes the covered corpus
wire contract, not every TeX input or every Office/browser environment.
Changed payload behaviour requires renewed qualification. The release metadata
commit still requires its own CI. Formal publication follows the release gates.


## 3.1.3 candidate coverage limits

The expanded Office gate distinguishes 36 matched fixtures from three non-covered
constructs: native `\boxed` and `\binom` omit their structure, and native
`\cancel` emits a bottom bar instead of Web's diagonal strike. These cases do
not qualify for Web hand-test exemption. Hosted wire/review and changed-boundary
qualification remain pending; 3.1.3 is not released. See PRODUCT_PARITY's fixture
matrix. PowerPoint colour loss on the shared OMML route remains a known limitation.
