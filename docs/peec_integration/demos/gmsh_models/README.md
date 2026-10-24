# Solver meshes and display artifacts

PEEC solver inputs are checked Cubit/Netgen `.vol` files with versioned label
contracts and SI coordinates. Gmsh 5.0.0 consumes derived MSH 4.2 solely for
post-processing. This directory contains checked `.vol` demo helpers; no
`.msh` mesh assets remain here. Derived display output is never imported
into PEEC or converted into `.vol`.

- `demo_gmsh_surface_mesh.py` constructs the existing Netgen torus and can save
  its `.vol` using `--write-vol`. Its existing area criterion is unchanged.
- `demo_vol_to_peec.py` reads an explicitly supplied checked `.vol` and builds
  the existing PEEC model.
- `demo_vol_ports_peec.py` uses the same solver-input boundary and the demo's
  declared geometric port selection.
- `solver_vol.py` uses Netgen's own mesh data and the canonical label checker;
  it implements no Gmsh reader or compatibility format.

```shell
python demo_vol_to_peec.py --mesh conductor.vol --contract conductor.labels.json
```

The legacy MSH solver reader and MSH solver-fixture generators are deleted.
Supply a mesh and contract from the responsible mesher; no mesh generation or
numerical substitution is performed while loading a demo.

The canonical label-contract check requires `cubit-mesh-export` (available via
`radia[cubit]`); it does not invoke licensed Cubit. Mesh access uses Netgen/NGSolve.

The versioned label contract must declare SI units `m` (top-level `units`
for a surface mesh, or `centerline.units` for a centerline). Other units
are rejected before native mesh loading; these demos do not convert units.
