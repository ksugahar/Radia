# TeX equation services

Equation editing and conversion belong to the independent
[`eqnedit64` package](../../packages/eqnedit64/README.md) and
[`native/Web product`](../../tools/eqnedit64/README.md).
The Radia solver no longer builds or distributes an equation module.

TeX is the source format. For an existing Office session, use
`radia_mcp.presentation.presentation_copy_equation` with the Eqnedit64 backend
to copy editable Office math. For explicitly selected image output, use
`presentation_render_equation`; images are not editable Office math.

For Markdown-to-PPTX generation with editable equations, use
`presentation_markdown_to_pptx` with `eqnedit64>=3.1.4`. The independent
Eqnedit64 backend exports checked OOXML math; the presentation adapter packages
it into the deck. Unsupported input raises without an image or legacy backend
fallback. `presentation_equation_policy()` reports the implementation and
required version. The `presentation_pptx_*` tools also edit presentation documents.

Standalone product contracts and tests remain under `tools/eqnedit64`;
the MCP adapter contracts are in
`packages/radia-mcp/tests/test_presentation_equation_cli.py` and
`packages/radia-mcp/tests/test_presentation_markdown_math.py`.
