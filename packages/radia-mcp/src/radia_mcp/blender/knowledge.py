"""User-facing guidance for the Blender animation domain."""

GUIDE = """\
# Blender 3-D animation in radia-mcp

## Roles
- Gmsh: mesh and field inspection, cut planes and sections of Radia results.
- Blender: 3-D model visualisation and animation (blender_field_animation).
- The solver (Radia main) is the source of truth; Blender only maps its data.

## Input contract: radia.field-animation.v1 (HDF5, SI, index base 0)
/time_s (n)                      physical time of each solver step
/solver_step (n)                 solver step index (optional, provenance)
/magnet/center_m (n,3)           magnet centre of that step
/magnet/rotation (n,3,3)         columns = magnet body u, v, w axes (orthonormal)
/magnet/size_m (3), /magnet/magnetization_A_per_m (3, body frame)
/conductor/box_min_m, box_max_m (3)
/samples/id (m) unique stable ids, /samples/point_m (m,3)
/field/J_A_per_m2 (n,m,3)        signed solved current density
Root attributes: schema, units, source script and its SHA-256, solver settings.
Example producer: docs/maglev/demos/rotating_magnet_eddy_animation_export.py
(rotating 1 mm NdFeB-like cube over a copper plate, full-FEM A-phi reference).

## Fidelity rules
- Magnet pose and J change only at solver steps: CONSTANT keyframes, no
  interpolation between samples.  One frame = one solver step.
- Fixed |J| colour scale for the whole clip (default: global max), shown in
  an on-screen legend with SI units.
- Arrow direction = in-plane direction of the computed J; length and colour
  are linear in |J| on the fixed scale.  Glyphs are drawn a small display
  offset above the plate top; the stored attribute keeps the solved point.
- Camera is static.  The playback rate (physical seconds per video second)
  is written on screen and in the scene properties.
- Nothing uncomputed is added: no field lines, sparks, levitation or extra motion.

## Worker setup (optional, outside the wheel)
- Python 3.13 venv with bpy: pip install -r packages/radia-mcp/tools/blender/requirement.txt
- RADIA_BLENDER_PYTHON = that interpreter.  It is started with -I -B.
- RADIA_BLENDER_WORKER = tools/blender/field_animation_worker.py; this may be
  omitted in an editable source checkout.  Missing settings give a clear error.
- Rendering uses Cycles on the CPU (no GPU); set "threads" in spec_json.

## Verification performed by the tool
A second process reopens the .blend without a UI and reports frame range,
fps, physical time per frame, |J| scale, camera animation, the magnet's world
pose and the stored (sample_id, J) of selected steps.  The server compares
these with the HDF5 dataset: J relative error <= 1e-6 (float32 storage),
centre <= 1e-8 m, rotation <= 1e-6, exactly one glyph object visible per
frame and the time label of that step.  Review the preview PNGs and MP4 for
readability; that judgement is separate from numerical validation.

## Design references (public Blender MCP servers, MIT; ideas reimplemented, no code copied)
- ahujasid/mcp-for-blender: scene summary fields (frame range, fps, engine,
  resolution, units, world bounding boxes), per-key interpolation choice
  (CONSTANT for stepped data), Blender 4.4+/5 slotted actions (channelbags),
  multi-frame camera-view checks.
- HoldMyBeer-gg/blenderwright, seehiong/blender-mcp-bridge, djeada/blender-mcp-server:
  typed animation/render operations, one-call multi-view/frame previews.
Deliberately NOT adopted: live-session sockets, execute-code tools, asset
download, generative content.  Here the reopened file is checked by
numbers (interpolation must be CONSTANT only; conductor bbox equals the
dataset), not by a vision judgement.

## Analytic mapping check
For a field-mapping test independent of the FEM, export a dataset whose J is
known in closed form (for example a rigidly rotating uniform field) and
confirm the readback values equal the closed form at the stored points.
"""
