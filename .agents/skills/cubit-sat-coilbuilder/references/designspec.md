# SAT CoilBuilder DesignSpec

Use this contract for every SAT-to-CoilBuilder migration.  Values marked
`REQUIRED` must be decided from design data or confirmed by the model owner;
they must never be inferred from a solid-coil CAD body.

```yaml
schema: radia.sat-coilbuilder-design/v1
name: ffag_main_bending_magnet

cad:
  source_sat: W:/path/to/model.SAT                 # REQUIRED
  sha256: <audit SHA-256>                           # REQUIRED
  source_length_unit: mm                            # REQUIRED
  metres_per_source_unit: 0.001                     # REQUIRED
  cubit_version: 2025.12                            # REQUIRED
  cubit_model: W:/path/to/imported_model.cub5       # REQUIRED

volumes:
  iron:                                             # Cubit volume IDs
    - 1
  coils:
    f_main:
      cubit_volume_ids: [2]                         # Geometry evidence only
      circuit: f_series
    d_main:
      cubit_volume_ids: [3]
      circuit: d_series
  discard: []                                       # Legacy air / construction bodies

circuits:
  f_series:
    current_a: 0.0                                  # REQUIRED, signed
    turns: 1                                        # REQUIRED, positive integer
    path:
      representation: explicit_polyline_m            # REQUIRED for first delivery
      points_m:                                      # Closed: first equals last
        - [0.0, 0.0, 0.0]
        - [0.0, 0.0, 0.0]
    conductor:
      profile: rectangle
      width_m: 0.0                                  # REQUIRED
      height_m: 0.0                                 # REQUIRED
    direction_reference: "positive current follows the listed point order"

open_region:
  exterior: kelvin
  physical_air_gaps:
    - name: pole_gap
      required_mesh_resolution: "declare before meshing"

acceptance:
  closure_tolerance_m: 1.0e-9
  source_probe_file: W:/path/to/jmag_probes.json    # REQUIRED when available
  jmag_reference_file: W:/path/to/jmag_reference.json
  required_components: [Bx, By, Bz]
  periodic_sector:
    enabled: true                                  # REQUIRED decision
    fold: 12                                       # REQUIRED if enabled
    sector_angle_deg: 30.0                         # 360 / fold
    rotation_axis: [0.0, 0.0, 1.0]                 # REQUIRED
    field_phase: periodic                          # periodic | antiperiodic
    body_crosses_periodic_planes: false            # REQUIRED
    periodic_boundaries: []                        # [periodic_min, periodic_max]
    hdiv_trace_identified: false                   # required if body crosses
    q2_geometry_relative_tolerance: 1.0e-9
    coil_source: full_ring_from_one_sector         # REQUIRED if enabled
```

## Rules

- `points_m` must describe the electrical centerline after the stated unit
  conversion.  Use a higher-level `CoilBuilder` recipe only when it regenerates
  those same points deterministically.
- Multiply either the current or the number of identical paths by `turns`, but
  document which representation is used.  Do not multiply both.
- A `cubit_volume_id` proves only the enclosing CAD region.  It does not define
  the electrical path.
- Keep imported external-air and construction volumes under `discard` for
  Kelvin HDiv solves unless a physical local air gap explicitly requires them.
- Save the SAT audit, Cubit volume inventory, DesignSpec, probe set, and every
  numerical comparison next to the validation result.
- A periodic source has two distinct contracts.  Rotate one declared coil
  source to the full ring for both the HDiv RHS and native tracking.  Use
  `image_cyclic=fold` only for the repeated solved iron charges.
- When `body_crosses_periodic_planes` is false, do not manufacture periodic
  CAD faces; a disjoint sector cell uses cyclic images without an HDiv trace
  identification.  When it is true, Cubit must create exact rotationally
  matched cut faces, Netgen must identify every cut vertex, and curved-Q2
  geometry must meet the stated relative tolerance before solving.
