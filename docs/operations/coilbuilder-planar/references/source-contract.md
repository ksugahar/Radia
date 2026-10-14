# Excitation identity for planar coil conversion

A planar CAD ribbon and its extrusion determine a winding-pack volume. They do not determine turns, electrical connection, signed current or the homogenized current distribution. These remain declared input data.

## Supported conversion

One planar ring face (two boundary wires), constant width, perpendicular extrusion, matched inner/outer line and concentric-circle segments. Full circles are split into four arcs. Adjacent sections must be tangent. The input may lie in any plane, use explicit CAD units, or be a known half-loop completed by a declared symmetry. Unsupported splitting, splines, nonconcentric bends or extra holes are rejected rather than silently approximated. More general profiles need a separately validated loft/filament route.

The helper preserves the specified width and measures it independently from paired boundaries. Centerline endpoints within the declared geometry tolerance may be snapped to make a closed chain; the certificate reports the change. Arc fitting and junction tangent checks expose rounded input coordinates.

## Current models are different

`CoilBuilder.to_radia()` constructs constant rectangular straight and arc sources. Arc primitives use uniform azimuthal current density. A resistive scalar-potential current in a curved solid generally has a radial distribution (approximately inverse radius in a concentric bend). Equal geometry and ampere-turns do not prove equal source fields. Actual insulated turns can follow another declared distribution.

Compare coil-only vector fields before attributing any full-model discrepancy to reduced-A, Omega or HDiv-MMM. Record core and fringe differences, geometric/source refinement, and current closure. If the source difference exceeds the engineering limit, keep the original conservative volume-current source for solver parity, or adopt a documented filament/current-distribution representation and validate it. Do not switch representations silently.

## Formulation identity

The shared inputs are SI geometry, material B(H)/H(B) law, declared source, boundary/symmetry convention and physical observation points. An original total-A calculation may use a different constitutive interpolant or finite exterior approximation: certify their differences explicitly.

An FE-vacuum split of total A is a useful algebraic check on a common mesh but is not evidence of a CoilBuilder-based reduced-A calculation. State which route each result used. Report true linear and nonlinear residuals separately from magnetic-field agreement.

## CAD Boolean diagnostics

For tangent arc/straight compounds, OCC Boolean masses can be inconsistent. A zero difference is not a certificate when union volume or bounds are invalid. Check primitive extents, reconstructed boundary samples and analytic swept volume independently; retain failed Boolean diagnostics as unverified. Do not infer a valid watertight solver solid from those diagnostics.

## CoilBuilder current and identity manifests

The converter sets `builder.current = current_A * turns`: native straight/arc sources use the total ampere-turns of the homogenized winding. The generated module override is therefore named `ampere_turns_A`. When using `to_identity_manifest`, pass `turns=1` for this already aggregated builder and attach the separately declared physical current and turn count to the conversion certificate; passing the physical turns again would double-count the excitation.
