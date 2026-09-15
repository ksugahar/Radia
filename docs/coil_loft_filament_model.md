# Rectangular Straight Loft: Prescribed Stream Tubes

`CoilBuilder.to_radia_loft_filaments(nw, nh)` exports straight rectangular
segments and straight linear rectangular lofts into native `ObjFlmCur` objects.
The existing solid `to_radia()` API does not silently select this approximation.

For canonical section coordinates a,b in [-1/2,1/2], use

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
- Arc lofts, arbitrary profiles, and discontinuous joins are rejected.
- Cross-section corners are checked even with a 1x1 sampling grid.

`tests/test_coil_builder_loft_native.py` checks an independent volume integral,
section convergence at three exterior points, the constant-section solid
limit, total current, reversal, and pre-allocation rejection. The independent
test oracle uses Gaussian quadrature; the exported model does not.
