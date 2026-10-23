# Fixed post-processing meshes

These ASCII MSH 4.2 fixtures have explicit, hand-authored connectivity; no Gmsh mesher is invoked.
The unit hex has volume 1. Tetrahedron boundary counts are four faces (4 or 10 nodes).
The adjacent hex/prism share one quadrilateral: 6 + 5 - 2 = 9 boundary faces, 10 nodes.
They exercise the production Gmsh post-processing reader only; they do not create solver-bound .vol files.

| File | SHA256 |
|---|---|
| hex_prism.msh | c5524ccf527f8b04da11f4e5acfb41b0c05941bc8d70c2b127fc54526057a46f |
| tet_p1.msh | d846c6d540ac891afdba5a4c25662636a69e9118638175b9a3bc06c9adce6175 |
| tet_p2.msh | 919148fac5df882d17f73eace6f3432ca761861a8591e046dd06c2634aa3b9e6 |
| unit_cube.msh | 1867f228287a73cd742c57d34f54bfdd088609857cb16a51c446d9f54e0166e9 |
