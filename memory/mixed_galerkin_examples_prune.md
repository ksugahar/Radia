# mixed_galerkin examples prune

2026-06-28: The cylinder/sphere analytic references were promoted from the
example `_references` tree into `radia.maglev.mixed_galerkin.references` so
docs notebooks and lightweight examples no longer import example-local helpers.

The deliberately broken 2D square v1 envelope was deleted from source. Lesson:
`psi = f(x) sin(pi y/L) + sin(pi x/L) f(y)` vanishes on the boundary but misses
the orthogonal boundary-layer structure. The canonical square example uses
`psi = f(x) f(y)`.

Do not recreate source-code archive notebooks for mixed_galerkin. Keep result
showcases in docs, reusable analytic formulas in src, runnable heavy checks in
validation_test, and lessons in memory.
