# Kelvin consumer API audit — 2026-10-28

This targeted audit covers producer-persisted periodic contracts. It does not
certify an electromagnetic formulation, curved-map accuracy or a field result.

| Family | Classification | Current action |
| --- | --- | --- |
| Post-load point matching and identification insertion | plumbing-delete | Remove the consumer mutator. CAD/mesh producers own correspondence and write persistent IDs. |
| Panel repair wrappers and aliases | deprecated-drop | Remove; migrate callers to explicit read-only checks. |
| `verify_kelvin_identification` | user-intent | Verify an explicit directed native relation, type and complete labelled coverage. It cannot prove a persisted name through the pinned Python mesh API. |
| `load_verified_kelvin_mesh` | user-intent | Parse named persisted metadata, load that exact immutable file with NGSolve, and cross-check actual native rows and boundary coverage. Return mesh and evidence together. |
| `detect_kelvin_offset` | method-demote | Retain as an observation-only coordinate summary used by existing physical-configuration callers. It never matches, repairs or accepts periodic topology. |

The public top-level family replaces three former entries with three current
entries; there is no compatibility alias. An internal structural verifier serves
solver constructors with in-memory/OCC meshes, reports native direction and
explicitly makes no persisted-name or producer-orientation claim. Solver
constructor signatures remain unchanged.

All consumers use the backend's own point IDs, face descriptors, BC properties
and periodic type. No replacement mesher, coordinate matcher or solver is added.
Names and directed persistence are checked at the file-load boundary; detached
metadata is never proof for an unrelated in-memory mesh. Case-specific pair
oracles, geometry fingerprints and space/field acceptance stay separate.
