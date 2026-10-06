---
name: coilbuilder-planar
description: Build checked CoilBuilder excitation models for planar and three-dimensional coils, including fusion-magnet workflows. Convert supported planar Netgen/OCC faces with the bundled converter; use declared spatial segment geometry for 3D coils and identify unsupported free-form conversion explicitly. Preserve current, turns, units and source physics across reduced-A, mixed Omega and HDiv-MMM.
---

# CoilBuilder: planar and spatial coils

Three-dimensional coils are part of the intended scope, including fusion toroidal-field, poloidal-field and nonplanar windings. Read [the spatial-coil contract](references/spatial-coils.md) before constructing a 3D model. Distinguish an existing spatial segment API from automatic conversion of an arbitrary 3D CAD winding: the latter is not implemented by the planar helper. Do not flatten a spatial winding into a plane to make the converter accept it.

Use the executable [planar converter](scripts/planar_face_to_coilbuilder.py), not a fresh hand-written centerline approximation, for a constant-width planar ribbon with one hole and a perpendicular rectangular extrusion. It recognizes matched line/circle inner and outer boundaries and writes a portable `build_coil(ampere_turns_A=None)` Python module. Unsupported curves, ambiguous pairs, variable widths and non-tangent corners fail explicitly. Read [the source contract](references/source-contract.md) before using a converted source in a solver.

## Netgen Python inputs

Inspect the original Python file before executing it. Extract only its geometry constructor into a job-owned module; avoid its solver loop, GUI calls, filesystem writes and meshing side effects. Retain its path/hash and parameter provenance. The helper takes a Netgen OCC `Face`; it does not parse arbitrary Python into a winding or infer electrical data.

Collect the declared CAD-to-metre scale, extrusion vector, current, turns, circuit connectivity and winding orientation. A known symmetry half may be joined to its declared mirror to create a closed planar face. Do not interpret its open half-loop as a closed coil or invent a symmetry. An existing extruded winding solid can supply its planar base face when the extrusion direction and thickness are independently known.

Call `face_to_coilbuilder(face, current_A=..., turns=..., extrusion_vector_m=..., section_width_m=..., path_normal=..., geometry_tolerance_m=...)`. The signed current multiplied by turns becomes CoilBuilder.current. `path_normal` specifies the positive loop orientation. Lengths in a CAD face use `length_scale_to_m`; the extrusion vector is always in metres. Preserve an explicitly specified section width; report the independently measured width.

## Gates and artifacts

- Retain the conversion certificate: planarity, paired boundary error, coordinate snapping, reconstructed centerline error, closure, tangent continuity and section dimensions. Select a tolerance from the source's coordinate precision and the engineering acceptance; never increase it simply to hide a mismatch.
- Export `build_coil()` using `write_builder_module`, reload it, and verify the current and closed path. Where available, use CoilBuilder's identity manifest and segment-overlap check. Confirm reconstructed boundary samples, analytic swept volume and primitive extents against the original; a matching centerline alone is insufficient. Treat inconsistent OCC Boolean masses as unverified, as described in the source contract.
- Compare the coil-only **vector** field at the solver's declared core and fringe points. Compare source refinement independently. Record the current-density convention and report every observation region; do not pass a core-only gate while claiming the entire field matches.
- Only then use the same accepted CoilBuilder excitation in reduced-A, mixed total/reduced Omega and HDiv-MMM. Keep an FE-vacuum A split explicitly labeled as discrete equivalence; it is not a CoilBuilder source validation.

Use Netgen directly when requested or when commercial CAD tools are disallowed. No Cubit dependency is needed. Keep project-specific geometry and results in the user's private project; skill examples and validation contain synthetic geometry only.

Validate helper changes with `python scripts/validate_planar_converter.py --output C:/temp/coilbuilder_planar_validation`. The script checks reconstructed OCC volume, arbitrary-plane handling, current/turn semantics, module reload and rejection of unsupported contracts.
