# Hodograph Differential-Geometry `.wls` Suite

This note is the human-readable index for the 2026-06-22 hodograph
differential-geometry commits.  The executable sources live in
`packages/radia-mcp/src/radia_mcp/mathematica/differential_geometry/`;
the MCP mirror is
`differential_forms_mathematica_recipes(topic=...)` in
`mcp-server-differential-forms`.

The purpose is narrow: before a FEM implementation changes a Kelvin weight,
hodograph chart, saturation linearization, or HOIBC boundary model, run the
symbolic file that fixes the corresponding geometric identity.

## Coverage

| topic | script | settled statement | MCP topic |
|---|---|---|---|
| Weak form / Hodge | `weakform_hodge.wls` | The weak form is a Hodge pairing.  A coordinate map moves the material/Hodge weight, not the exterior derivative.  Pullback Kelvin weights and transformation-optics material tensors are the same cofactor rule. | `weakform_hodge` |
| Hodograph backbone | `hodograph.wls` | Kelvin, Clebsch, potential-plane, field-plane, and Chaplygin transforms are cells of one diagram.  Chaplygin linearizes the 2-D saturable field-plane problem; ordinary magnetic saturation has no limiting line. | `hodograph` |
| Canonical structure | `canonical.wls` | Flux lines are Hamiltonian trajectories with `A_z` as Hamiltonian.  `H = dw/dB` is a Legendre transform, and the hodograph is a canonical transform of the potential. | `canonical` |
| Surface de Rham / HOIBC | `surface_derham.wls` | HOIBC has a topological part (surface harmonic 1-forms counted by `b1`) and a geometric part (local Pade/Laplace-Beltrami approximation of the analytic exterior DtN); vector HOIBC splits into two de Rham Steklov ladders. | `surface_derham` |
| DtN / Steklov geometry | `dtn_geometry.wls` | The DtN operator is the exterior Hodge star condensed to the boundary: Schur complement, radial Riccati fixed point, shifted square-root of the surface Laplacian, self-adjoint positive boundary metric. | `dtn_geometry` |

The suite index is exposed as MCP topic `differential_geometry`.  The topology
companion is `../basis_functions/cohomology.wls`, which locks the fact that
the harmonic-loop count is Betti data and cannot be changed by a Hodge/material
metric.

## Run

```powershell
cd packages/radia-mcp/src/radia_mcp/mathematica/differential_geometry
wolframscript -file weakform_hodge.wls
wolframscript -file hodograph.wls
wolframscript -file canonical.wls
wolframscript -file surface_derham.wls
wolframscript -file dtn_geometry.wls
wolframscript -file ../basis_functions/cohomology.wls
```

Each file is self-contained and prints `ALL PASS` when its symbolic assertions
hold.  If Mathematica is being driven through MCP, send the file content to
`mathematica_evaluate(code=<file text>, timeout=120)`.

## Where It Fits

`HODOGRAPH_BACKBONE.md` is the conceptual map.  This suite is the executable
symbolic guardrail for that map:

- `dB = 0` and the de Rham/cochain structure are metric-free.
- The Hodge star carries metric, material, and coordinate-map weights.
- Convex magnetic energy makes the nonlinear material tangent SPD, so
  saturation deforms the metric but does not fold it.
- 3-D global Clebsch variables are obstructed by helicity.
- HOIBC global loop counts are surface topology; impedance values are the
  surface Hodge/Steklov geometry.
- DtN is nonlocal exterior geometry condensed to the boundary; local HOIBC can
  only approximate its square-root spectrum.

That split is the rule of thumb: if a proposed implementation changes topology,
look at de Rham / cohomology first; if it changes weights, look at the Hodge
script first.
