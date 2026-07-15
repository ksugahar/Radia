# HDiv-VIM axis-aligned closed-form charge Gram rejected

2026-07-05.  The semi-closed charge-Gram entry idea for HDiv-VIM hex/wedge
blocks was retired from public docs after scope review.

What was tested:

- Inner closed Newtonian box potential plus symmetric tensor-Gauss outer
  quadrature for Q1 polynomial charges.
- The formula works as an analytic oracle for axis-aligned rectangular boxes.
- The arithmetic is real.  Removable singularities and branch choices need real
  limits / real quadrant handling; promoting the expression to complex
  arithmetic is not the fix.

Why it was rejected as a production route:

- The derivation assumes Cartesian-product, axis-aligned boxes with monomial
  charges in box-local axis coordinates.
- It does not cover warped, curved, swept, or general Cubit hex/wedge cells
  without a separate transformation theory.
- Publishing it as "the" HDiv-VIM hex/wedge near-block replacement would invite
  an invalid C++ port for the actual production mesh path.

Current rule:

- Keep `QuadBlockHex` / `QuadBlockWedge` as the production charge-Gram path for
  general HDiv-VIM hex/wedge cells.
- Use the box closed forms only as a small local oracle for axis-aligned
  fixtures if needed.
- Do not revive the public implementation plan unless the missing geometric
  transformation theory and validation for non-axis-aligned cells are added.

## Reflection-invariance lesson (2026-07-15)

- Mirroring only the stored upper triangle guarantees an algebraically
  symmetric matrix, but it does not make a one-sided finite quadrature rule
  invariant under replacement by an explicitly reflected mesh.
- The former FAR one-sided HEX rule left a roughly `1e-5` reduced/full field
  defect; increasing quadrature only hid it.  Explicitly averaging directed
  `AB` and `BA` rules restored the multicell full-vs-image field comparison to
  roundoff at the normal quadrature order.
- Production HEX, WEDGE, and high-order TET FAR blocks therefore average both
  directions.  The one-sided environment switches are diagnostic-only and
  must not be used for release results.
