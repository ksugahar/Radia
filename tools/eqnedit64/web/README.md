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

Inputs and saved/copied source are TeX. MTEF and `.eqn` are not supported
formats. Office copy emits inline MathML in `text/html` plus plain TeX; image
copy is a separate action so PowerPoint cannot accidentally prefer PNG over an
editable equation.
