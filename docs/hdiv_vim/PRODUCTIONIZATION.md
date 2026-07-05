# HDiv-VIM Production Note

Radia soft iron is standardized on HDiv-VIM.  This note is the current
production checklist, not a migration archive.

## Done

- Mesh-backed `radia.vim.MeshSoftIron` integrates with `rad.Solve`.
- The charge-Gram backend is C++/HACApK based and TaskManager-aware.
- TET/HEX/WEDGE and 2D planar paths have validation coverage under
  `validation_test/feec/`.
- `rad.Fld` is part of the public contract after HDiv write-back.
- MCP `hdiv_vim` documents the live API and reduced-FEM coupling policy.

## Release Gate

Before release or mdx deployment:

- run focused HDiv smoke tests on LAB;
- run heavy validation/benchmark sweeps on mdx when idle;
- record charge count, HDiv DoF, H-matrix compression, build time, solve time,
  iteration count, and machine label;
- verify image symmetry with an explicit full model when the mesh is truly
  symmetric;
- keep public docs free of obsolete backend names and local validation
  provenance.

## Open Work

- tighten image-symmetry roundoff contracts, especially for full-model hex
  comparisons;
- extend `rad.Fld` and force/energy tests around motor workflows;
- keep Cubit/GMSH mesh export aligned with the HDiv API;
- continue mdx scaling measurements for charge-Gram build and solve time.
