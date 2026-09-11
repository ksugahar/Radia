# Reflection audit follow-up

The dictionary optimization in `4aa3e3d28` retained the rounded-coordinate
membership gate but changed nearest-distance semantics: `setdefault` selected
the first point in a rounding cell. Coincident rounded keys and neighboring
rounding cells can contain a closer point. A three-point counterexample gives
zero in the original audit and 4e-15 m in the first-key implementation.

The follow-up retains the original rounded-key gate. Exact reflected points
use set membership and need no spatial tree. For non-exact matches, a cKDTree
returns the nearest distance without an all-pairs scan. The exact-mirror path
has expected linear lookup cost; the general tree path is not advertised as
strictly O(N). Vertex-order permutations do not change the result.

## Verification

- 33 focused tests passed: the new reflection tests and existing gap acceptance
  tests. Cases include rounding collisions, a nearest point in an adjacent
  rounding cell, duplicates, missing keys, perturbed symmetric clouds, and
  invalid coordinates. The reference is the original all-pairs algorithm.
- The N=6 `kelvin_domain.vol` from `C:/temp/ctype_gapfamily2/n06` was checked
  again. All five output fields match the saved `mesh_result.json`:
  18,065 vertices, 104,794 volume elements, zero missing reflected vertices,
  zero maximum reflected vertex error, and zero missing reflected elements.
  The full audit took 9.94 seconds on LAB (one observation, not a benchmark).

## Acceptance scope

This audit checks vertices and material-labelled element vertex signatures.
It does not certify the complete high-order geometry map or magnetic-field
reflection parity. Gap line measurements certify only the named probe lines.

At this review, N=24 had a passing independent early check but no official
`n24/mesh_result.json`. Do not synthesize that result or silently substitute the
early check. Before completing the family, compare the official checks with
the early check and record the hashes of each mesh and result in the manifest.

The N=6/12/24 family holds the other requested mesh sizes fixed. A subsequent
field study measures sensitivity to this gap-refinement sequence. Neither mesh
acceptance nor observed convergence alone is a rigorous discretization-error
bound or a substitute for whole-model convergence and three-engine validation.
