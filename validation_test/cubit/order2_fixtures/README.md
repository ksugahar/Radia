# Order-2 Netgen fixtures

`mixed_cylinder_o2.jou` produces `tests/fixtures/cubit/mixed_cylinder_o2.vol.gz`; `tet_sphere_o2.jou` is a
second order-2 case with no committed fixture.
Run each journal in a separate empty directory using Coreform Cubit 2026.8 with the radia-mcp
`export netgen` command plugin installed. Use Cubit's batch console:

```powershell
& "$cubit\bin\cubitx.exe" -batch -nographics -nojournal -noinitfile -workingdir "$PWD" -input mixed_cylinder_o2.jou
```

For the sphere, substitute `tet_sphere_o2.jou`. Each journal saves the input as `model.cub5` and writes
`order2.vol` with geometry order 2. Keep the generated `order2.vol.json` with the validation record.
Compress the .vol bytes directly for consumers that expect .vol.gz.

| Fixture | Geometry | Volume elements recorded when the fixture was made |
| --- | --- | --- |
| Mixed cylinder | Radius 1 m, height 2 m | 111 hex, 729 tet, 37 pyramid |
| Tet sphere | Radius 1 m | 231 tet |

Read the exported geometry at its stored order. Compare the cylinder with volume 2π and area 6π,
and the sphere with volume 4π/3 and area 4π.
