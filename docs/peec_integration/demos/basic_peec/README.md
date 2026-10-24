# basic_peec — Basic PEEC Demos

Demonstrations of core PEEC functionality.

## Files

| File | Description |
|------|-------------|
| `demo_peec_simple_loop.py` | Simple loop PEEC analysis |
| `demo_peec_dc.py` | DC analysis (resistance computation) |
| `demo_peec_from_vol.py` | Build PEEC model from a checked Netgen .vol |
| `demo_peec_from_dual_mesh.py` | Build from dual mesh |
| `demo_peec_from_1d_mesh.py` | Build from 1D centerline mesh |
| `demo_frequency_response.py` | Frequency response computation |
| `demo_loop_star.py` | Loop-Star decomposition demo |
| `demo_port_extraction.py` | Port impedance extraction |

## Outputs

- `demo_peec_frequency_response.png` — Frequency response
- `demo_loop_star_coil.png`, `demo_loop_star_simple.png` — Loop-Star results

Mesh-backed demos require explicit SI solver inputs and their versioned label
contract. For example:

```shell
python demo_peec_from_vol.py --mesh conductor.vol --contract conductor.labels.json
```

The mesh must originate in Cubit or Netgen. These demos never read or convert
Gmsh output. Material properties, cross-sections, excitation and solve settings
remain defined by each demo. No bundled display mesh is used as solver input.

The canonical label-contract check requires `cubit-mesh-export` (available via
`radia[cubit]`); it does not invoke licensed Cubit. Mesh access uses Netgen/NGSolve.

The versioned label contract must declare SI units `m` (top-level `units`
for a surface mesh, or `centerline.units` for a centerline). Other units
are rejected before native mesh loading; these demos do not convert units.
