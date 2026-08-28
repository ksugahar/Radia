# TeX equation services

`radia.equation` is the headless TeX-to-document layer used by Radia and
`radia_mcp.presentation`. It parses TeX directly and provides:

- editable Office Math (OMML) for Word and PowerPoint;
- MathML, RTF, SVG, EMF, PNG, and DIB output;
- native-equation `.docx` and `.pptx` generation from Markdown;
- a structural editing model and Markdown math scanner for automated workflows.

The human-facing native and browser editors are maintained together in
[`tools/eqnedit64`](../../tools/eqnedit64). The standalone signed executable
is released separately and is not bundled into the Radia Python wheel.

## Supported source contract

TeX is the only equation source format. MTEF and `.eqn` readers, writers,
converters, fixtures, and the old duplicate GUI were retired when Eqnedit64
became the canonical editor. Git history preserves the former migration code
without keeping it in the supported product.

## Presentation use

For a complete deck, prefer one-process generation:

```python
from radia.equation import markdown_to_pptx

markdown_to_pptx(markdown, "deck.pptx")
```

This inserts OMML rather than equation screenshots. For an existing interactive
PowerPoint/Word session, Google Slides, or a standalone PNG/EMF, use the
`radia_mcp.presentation` Eqnedit64 tools documented by
`presentation_equation_policy()`.

## Verification

Focused tests live in `tests/equation`. The product-level native/Web parity
contract is
[`tools/eqnedit64/docs/PRODUCT_PARITY.md`](../../tools/eqnedit64/docs/PRODUCT_PARITY.md).
