# High-order curved HEX mesh in NGSolve

Demonstrates the `cubit-mesh-export` capability that Netgen alone does not provide:
**curved, high-order hexahedral meshes consumed by NGSolve.**

A sphere meshed with 56 hexes (`volume scheme sphere` in Coreform Cubit) is exported to a
high-order Netgen `.vol` (`block 1 add hex all; export netgen "hexsph_oN.vol" order N`).
NGSolve loads each `.vol` and integrates the volume; the curved hex boundary makes it
converge to the analytic `4/3 pi r^3`:

| order | element | volume error |
|:-----:|:-------:|:------------:|
| 1 (straight hex) | HEX × 56 | **-23.4 %** |
| 2 (curved hex)   | HEX × 56 | **-0.2 %** |
| 3 (curved hex)   | HEX × 56 | **+0.1 %** |

## Run

The `.vol` meshes are generated artifacts and are not committed. Regenerate
`hexsph_o1.vol`, `hexsph_o2.vol`, and `hexsph_o3.vol` with Cubit first, then
run the NGSolve checks from this directory:

```bash
python hex_sphere_curved_ngsolve.py
```

## Gotcha

A high-order `.vol` already carries its curved mid-side nodes. Load it with plain
`Mesh(path)` and **do not** call `mesh.Curve()` — `mesh.Curve(p)` re-curves from the
underlying CAD geometry (which a loaded `.vol` does not have) and resets every element
to straight-sided (the volume jumps back to the −23 % linear value). `mesh.Curve()` is
the right call only for meshes built from an in-memory CAD object (e.g. `netgen.occ` /
`SplineGeometry`).

## Regenerating the meshes

```
volume 1 size 0.012 ; volume 1 scheme sphere ; mesh volume 1
block 1 add hex all
export netgen "hexsph_oN.vol" order N overwrite
```
