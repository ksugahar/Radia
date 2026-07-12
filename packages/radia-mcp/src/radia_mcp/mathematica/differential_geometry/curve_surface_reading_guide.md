# Curve / Surface Differential Geometry Reading Guide

This guide connects the local Mathematica textbook
`Mathematica 曲線と曲面の微分幾何.pdf` to the Radia MCP symbolic
differential-geometry suite.  The PDF itself is a local literature asset and
is not copied into this package; this file records the study map, search terms,
and executable Radia anchors.

Local use:

```powershell
$env:LAB_MATHEMATICA_CURVE_SURFACE_DG_PDF = "<local path to the PDF>"
```

Do not commit extracted book text or machine-local absolute paths.  Use the
book as background reading, then encode Radia-facing claims as original notes,
small symbolic Mathematica checks, or self-tested `.wls` assertions.

## Why this book belongs here

Hodograph and stream-function workflows both depend on the same geometric
language:

| book theme | Radia meaning | executable anchor |
|---|---|---|
| Parametric curves, tangent, normal, curvature | stream-function contours are current paths; single-stroke wiring needs curve geometry | `hodograph.wls`, stream-function contour validation |
| Surface parametrization and surface normal | winding surfaces, current sheets, and curved former geometry | `surface_derham.wls`, `dtn_geometry.wls` |
| First fundamental form / metric tensor | surface gradients, area weights, and Hodge-star material weights | `weakform_hodge.wls` |
| Second fundamental form / curvature | curved elements, geometric error, and whether a surface model is faithful | curved-mesh validation and HOIBC geometry notes |
| Coordinate transforms and pullbacks | Kelvin maps, hodograph charts, and transformation-optics style material modulation | `weakform_hodge.wls`, `hodograph.wls` |
| Differential operators on surfaces | `grad_Gamma`, `div_Gamma`, Laplace-Beltrami, and surface DtN spectra | `surface_derham.wls`, `dtn_geometry.wls` |

## Hodograph reading route

Read the book's curve/surface chapters with these questions in mind:

1. What is the parameter space, and what is the physical space?
2. Which quantities are topological or differential-form data (`d`, pullback,
   closedness), and which quantities are metric data (inner product, area,
   Hodge star)?
3. When a coordinate chart changes, which part of the weak form changes?

Then run:

```powershell
wolframscript -file hodograph.wls
wolframscript -file canonical.wls
wolframscript -file weakform_hodge.wls
```

Interpretation:

- `hodograph.wls` is the coordinate-transform backbone: Kelvin, Clebsch,
  Chaplygin, and the `A_z` potential.
- `canonical.wls` is the Hamiltonian / Legendre reading of flux lines.
- `weakform_hodge.wls` explains why metric changes enter as Hodge/material
  weights rather than as changes to `dB = 0`.

## Stream-function reading route

For stream-function coil design, the key surface identity is

```text
K = n x grad_Gamma psi
```

where `psi` is a scalar on the winding surface, `grad_Gamma` is the surface
gradient, and `n` is the oriented surface normal.  Equal increments of `psi`
define current contours; connecting those contours into one manufacturable wire
is a curve-geometry problem on the surface.

Use the book to keep three implementation contracts straight:

- The surface metric controls `grad_Gamma psi` and `|K|`; it is not optional on
  curved formers.
- The surface normal orientation controls the sign of `K`.
- A contour is a level set on the surface, not a planar polyline unless the
  surface chart is actually planar.

Relevant Radia anchors:

- `radia_mcp.streamfunction.streamfunction("single_stroke")`
- `radia_mcp.streamfunction.streamfunction("regularized")`
- `radia_mcp.streamfunction.streamfunction("fusion")`
- `surface_derham.wls` for the surface de Rham split
- `dtn_geometry.wls` for the surface Laplace-Beltrami / DtN viewpoint

## What to promote from reading into MCP knowledge

Promote only distilled, Radia-specific claims:

- a Mathematica snippet that verifies an identity used in Radia;
- a `.wls` self-test when the identity becomes policy;
- a short note that maps a geometry concept to a Radia API, validation test,
  or notebook;
- a search term list for local PDF lookup.

Do not promote:

- scanned page images;
- copied book prose;
- long formula lists that are not tied to Radia behavior;
- local absolute paths.

## Search terms for local study

Use these terms when searching the OCR text locally:

```text
曲線, 曲率, 接ベクトル, 法線, 曲面, 接平面, 第一基本形式,
第二基本形式, 測地線, ガウス曲率, 平均曲率, パラメータ表示,
面積要素, 座標変換, Mathematica
```

Map any useful result back to one of the executable anchors above before
turning it into MCP guidance.
