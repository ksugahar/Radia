# Three-dimensional CoilBuilder workflow

## Current implementation boundary

`src/radia/coil_builder.py` provides `CoilBuilder.set_start(position, orientation)`, straight and circular-arc segments with `tilt`, loft segments, rigid transforms and a `saddle_coil` constructor. These can describe spatial chains when their geometry is explicitly supplied. Review the actual implementation in the selected wheel/source before choosing a route. A tilted planar loop remains planar; it is not evidence of arbitrary nonplanar winding support.

The bundled `face_to_coilbuilder` converter supports planar ribbons only. There is no certified automatic extractor of a general spatial centerline and section frame from an arbitrary Netgen/OCC solid in this skill. Spline, helical and general nonplanar CAD conversion require an implementation and independent validation before being advertised as supported. Never silently replace them with a planar loop or an unverified polyline.

## Required physical and geometric inputs

Record the centerline with its parameterization, cross-section shape and dimensions, section orientation along the path, CAD-to-metre scale, current, turns and circuit connectivity. Declare whether the source is a homogenized winding pack or explicitly resolved conductors. Specify joints, leads and return paths; an open curve is not a closed magnetostatic current circuit. For a fusion assembly, distinguish each circuit and its excitation; rotational copies do not imply identical polarity or electrical connection.

For an explicitly supplied segment chain, use the existing public APIs. Check right-handed orthonormal frames, position continuity, tangent continuity and section alignment at every junction. A closed centerline alone does not certify closure of a twisted rectangular section. Do not call `close()` to alter the intended geometry without recording and accepting the measured changes. Report bend inner radius, section distortion, self-intersection, segment overlaps and clearance to other coils and structures. A loft or filament approximation must preserve its declared current distribution and convergence contract.

For a future free-form converter, define how the section frame is transported, how deliberate twist is represented, and how inflections or near-zero curvature are handled. Do not assume a Frenet frame is defined everywhere. Bound deviations of both centerline and swept boundary; report closure of position, tangent and frame separately. Unsupported topology or a failed reconstruction must raise, rather than returning a visually plausible source.

## Acceptance sequence

1. Validate synthetic spatial coils before using customer geometry: a closed saddle, a chain with different bend planes, and a nonplanar closed curve with a declared frame. The last case is a development target until its representation and reconstruction are implemented and tested.
2. Check exported/reloaded excitation identity and CAD geometry, including sections and joints. Retain numerical tolerances, source hashes and failures. Do not treat a successful CAD Boolean alone as acceptance.
3. Compare vector fields against an independent Biot–Savart or resolved-conductor reference at declared core and fringe points. Refine path, section and quadrature independently; check current conservation and polarity. A filament reference does not automatically validate a finite winding-pack distribution.
4. Use the accepted identical source in the selected formulations. Separate source representation, material law, exterior boundary and solver discretization errors. A fusion-device accuracy claim needs device-specific acceptance evidence; the existing planar validation does not provide it.

Run numerical acceptance in the selected verified compute runtime; do not compare hosts using different editable or installed implementations. Store project CAD and results privately. Repository examples must use public synthetic geometry and must not require Cubit or a commercial license.
