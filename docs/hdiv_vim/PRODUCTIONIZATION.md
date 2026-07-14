# HDiv-VIM Production Note

Radia soft iron is standardized on HDiv-VIM.  This note is the current
production checklist, not a migration archive.

## Done

- Mesh-backed `radia.vim.MeshSoftIron` integrates with `rad.Solve`.
- The charge-Gram backend is C++/HACApK based and TaskManager-aware.
- TET/HEX/WEDGE and 2D planar paths have validation coverage under
  `validation_test/feec/`.
- `rad.Fld` is part of the public contract after HDiv write-back.
- RT1 solve results own a persistent C++ field evaluator: NumPy target buffers,
  one-pass IMA, TaskManager observation parallelism, analytic tet near kernels,
  and a guarded large-map quadrupole tree.  IMA auto evaluation remains direct
  to preserve the reduced/full roundoff contract.
- On geometrically and topologically symmetric full/reduced hex meshes,
  `rad.Fld` image parity is locked to the explicit full solve at the
  roundoff-level contract (`< 10 eps` relative error).
- RT1 is public for pure TET/HEX/WEDGE, planar 2D, IMA, and field evaluation.
- RT2 is public for flat pure-TET linear/nonlinear material solves and the NGSolve
  `ChargeGram`/`DemagOperator` surface.  RT2 topology/image/field extensions
  and curved RT2 remain fail-loud until separately implemented and validated.
- MCP `hdiv_vim` documents the live API and reduced-FEM coupling policy.

## Release Gate

Before release or `mdx`/`hibino` deployment:

- run focused HDiv smoke tests on LAB/100号機;
- run heavy validation/benchmark sweeps on an idle `mdx` or `hibino` host
  (`mdx` by default, `hibino` for MATLAB, large-memory, long-running, or
  mdx-occupied jobs);
- record the actual validation host in the result JSON/log;
- record charge count, HDiv DoF, H-matrix compression, build time, solve time,
  iteration count, and machine label;
- record field source count, evaluator build time, selected direct/tree route,
  observation count, direct-reference error, and public `rad.Fld` wall time;
- verify image symmetry with an explicit full model when the mesh is truly
  symmetric and enforce the roundoff contract;
- keep public docs free of obsolete backend names and local validation
  provenance.

## Open Work

- extend force/energy tests around motor workflows;
- keep Cubit/GMSH mesh export aligned with the HDiv API;
- continue mdx scaling measurements for charge-Gram build and solve time.
