# Rectangular Loft: Prescribed Stream Tubes

`CoilBuilder.to_radia_loft_filaments(nw, nh)` exports straight rectangular
segments and straight linear rectangular lofts into native `ObjFlmCur` objects.
The existing solid `to_radia()` API does not silently select this approximation.

Positive-angle rectangular arc lofts are also supported, with an explicit
`n_arc` chord count (default 64). Their centerline and interpolated section
are sampled on the arc; each chord uses the closed-form straight-filament
kernel. This is not an analytic curved-filament integral. Refine `n_arc`
independently of `nw` and `nh`. Radii that let the section reach the arc axis
and angles outside (0, 360] are rejected before allocation.

For a straight loft and canonical section coordinates a,b in [-1/2,1/2], use

    X(a,b,s) = origin + a*w(s)*ex + s*L*ey + b*h(s)*ez

where w and h interpolate linearly between positive endpoint dimensions.
Each path is a straight line. Its Biot-Savart line integral is evaluated by
the native closed-form filament kernel, without longitudinal Gaussian
quadrature. The section integral uses a uniform midpoint grid, assigning
current I/(nw*nh) to every path. Doubling nw and nh is still required for
section convergence; this is not an exact finite-volume field formula.

The continuous current model is the flux-preserving push-forward

    J(X) = I * (dX/ds) / det(dX/d(a,b,s))

using a consistently oriented local frame. The longitudinal component is
I/(w*h); transverse components follow expanding/contracting stream tubes.
This satisfies internal current continuity but does not solve Ohm's law
with equipotential terminals. A complete magnetostatic circuit still needs
a return path; an open segment is only a source-field contribution.

## Limits

- External-field approximation only: individual filaments are singular.
- Not suitable for internal fields, self-energy, self-force, or Joule loss.
- No skin effect, proximity effect, or resistance-weighted current model.
- Negative-angle arc lofts, arbitrary profiles, ordinary ArcSegment objects,
  and discontinuous joins are rejected by this explicit approximation API.
- Cross-section corners are checked even with a 1x1 sampling grid.

`tests/test_coil_builder_loft_native.py` checks an independent volume integral,
section convergence at three exterior points, the constant-section solid
limit, total current, reversal, and pre-allocation rejection. The independent
test oracle uses Gaussian quadrature; the exported model does not.
`tests/test_coil_builder_curved_loft.py` additionally checks a curved-volume
oracle, chord convergence, rigid transforms, and a bend-to-straight join.

## CAD Export

Matching rectangular or circular straight lofts use centered endpoint wires. Arc
lofts use interpolated section wires on the same circular centerline and
linearly varying dimensions as the current model. The curved CAD sides
are approximations: refine the segment's `n_sub` independently of field
sampling and check geometry convergence. CAD accepts `0 < angle < 360`,
positive dimensions, a clear inner radius, and `n_sub >= 4`; negative and
closed bends and cross-type profiles remain unsupported. Circular lofts
interpolate radius linearly. OCC section interpolation can have nonmonotone
volume errors; check the analytic volume and successive refined geometries,
not just a presumed convergence rate. Circular CAD support does not imply
circular support in `to_radia_loft_filaments`. This does not change the native-solid
restrictions of `to_radia()`.

`tests/test_coil_builder_loft_cad.py` checks analytic volume, refinement,
rigid poses, STEP round trips and unsupported-geometry rejection.
