# Frequency-domain verification

Hierarchical Cauer Y(s) verification across 2D and 3D geometries.

| Script | Geometry | Reference | Result |
|---|---|---|---|
| `circle_cpe_schur_plot.py` | 2D Cu cylinder, $a\!=\!5$ mm | exact Bessel | basis 4 → 17% wall-band, sub-0.01% asymptote |
| `circle_multik_schur_plot.py` | 2D cylinder, Multi-K vs canonical | exact Bessel | M0=3+4pt+Schur basis 8 floor still 17.5% (no improvement over basis 4) |
| `3d_sphere.py` | 3D sphere | Bessel-like closed form | Verifies $K_{\rm SIBC} = 4\pi R^2 \sqrt{\sigma/\mu}$ |
| `3d_cuboid.py` | 3D cuboid | Mellin asymptote + 64-corner sum | Verifies $K_{\rm SIBC} = (2(ab+bc+ca))\sqrt{\sigma/\mu}$ + edge corrections |
| `ngsolve_cuboid_Y.py` | 3D cuboid via NGSolve FEM | itself (production) | Full FEM Y(s) sweep for ROM extraction |
| `ngsolve_cube_minimal.py` | Minimal NGSolve cube | sanity | Small-scale FEM test |
