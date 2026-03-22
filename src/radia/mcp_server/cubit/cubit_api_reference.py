"""
Cubit Python API Reference (CubitInterface namespace).

Organized reference of the cubit.* Python API functions and entity
class methods, sourced from the official Coreform Cubit 2025.8 docs.

Source: https://coreform.com/cubit_help/cubithelp.htm
  appendix/python/namespace_cubit_interface.htm
"""

API_CORE = """
# Core Functions (cubit.*)

## Command Execution

```python
cubit.cmd(command_string) -> bool
    # Execute a Cubit command. Returns True on success.
    cubit.cmd("brick x 10")

cubit.silent_cmd(command_string) -> bool
    # Execute without echoing to console.
    cubit.silent_cmd("volume all size 1")
```

## Initialization

```python
cubit.init(args_list)
    # Initialize Cubit engine.
    cubit.init(["cubit", "-nojournal", "-batch"])

cubit.destroy()
    # Shutdown Cubit engine. Only one init/destroy cycle per process.
```

## Entity Parsing

```python
cubit.parse_cubit_list(entity_type, list_string) -> tuple[int]
    # Parse a Cubit entity selection string into IDs.
    # entity_type: "volume", "surface", "curve", "vertex",
    #              "body", "group", "block", "nodeset", "sideset",
    #              "hex", "tet", "tri", "quad", "face", "edge", "node",
    #              "wedge", "pyramid", "sphere"

    # Examples:
    cubit.parse_cubit_list("volume", "all")
    cubit.parse_cubit_list("surface", "all in volume 1")
    cubit.parse_cubit_list("hex", "all in block 1")
    cubit.parse_cubit_list("node", "all in surface 3")
    cubit.parse_cubit_list("surface", "all with is_merged=false")
    cubit.parse_cubit_list("tet", "all with x_coord > 0")
    cubit.parse_cubit_list("volume", "all except 1 3 5")

cubit.get_id_string(entity_ids) -> str
    # Convert list of IDs to space-separated string.
    ids = cubit.parse_cubit_list("volume", "all")
    cubit.cmd(f"mesh volume {cubit.get_id_string(ids)}")

cubit.string_from_id_list(ids) -> str
    # Convert list of IDs to compact range string (e.g., "1 to 5 7 9").
```

## Entity ID Management

```python
cubit.get_last_id(entity_type) -> int
    # Get the highest ID of the given entity type.
    cubit.cmd("brick x 10")
    vol_id = cubit.get_last_id("volume")  # e.g., 1

cubit.entity_exists(entity_type, id) -> bool
    # Check if entity exists.
    if cubit.entity_exists("volume", 5):
        print("Volume 5 exists")

cubit.get_entities(entity_type) -> tuple[int]
    # Get all entity IDs of a type.
    all_vols = cubit.get_entities("volume")

cubit.get_id_from_name(name) -> int
    # Get entity ID from its name.
    vid = cubit.get_id_from_name("my_block")

cubit.get_all_ids_from_name(geo_type, name) -> tuple[int]
    # Get all IDs matching a name.
    ids = cubit.get_all_ids_from_name("volume", "part_*")
```

## Entity Naming

```python
cubit.get_entity_name(entity_type, entity_id) -> str
cubit.set_entity_name(entity_type, entity_id, new_name) -> bool
    vol_name = cubit.get_entity_name("volume", 1)
    cubit.set_entity_name("volume", 1, "matrix")
```

## APREPRO Variables

```python
cubit.get_aprepro_vars() -> list[str]
    # List all APREPRO variable names.

cubit.get_aprepro_numeric_value(variable_name) -> float
    # Get numeric APREPRO variable value.
    cubit.cmd("{radius = 5.0}")
    r = cubit.get_aprepro_numeric_value("radius")  # 5.0

cubit.get_aprepro_value_as_string(variable_name) -> str
```

## Version Info

```python
cubit.get_version() -> str        # e.g., "2025.8"
cubit.get_python_version() -> str # e.g., "3.10.12"
cubit.get_exodus_version() -> str
cubit.get_acis_version() -> str
```

Source: coreform.com/cubit_help - CubitInterface namespace
"""

API_GEOMETRY_QUERIES = """
# Geometry Query Functions

## Coordinates and Positions

```python
cubit.get_center_point(entity_type, entity_id) -> tuple[float, float, float]
    # Works for any entity: volume, surface, curve, vertex
    center = cubit.get_center_point("volume", 1)
    # Returns (x, y, z)

cubit.is_point_contained(geometry_type, entity_id, xyz_point) -> int
    # Check if a point is inside a geometry entity.
    # Returns: 0=outside, 1=inside, 2=on boundary
    result = cubit.is_point_contained("volume", 1, (0.5, 0.5, 0.5))
```

## Measurements

```python
cubit.get_curve_length(curve_id) -> float
cubit.get_surface_area(surface_id) -> float
cubit.get_volume_volume(vol_id) -> float
cubit.get_total_volume(volume_list) -> float

cubit.get_distance_between(vertex_id_1, vertex_id_2) -> float
cubit.get_distance_between_entities(type1, id1, type2, id2) -> float
```

## Bounding Box

```python
cubit.get_bounding_box(geometry_type, entity_id) -> tuple[10 floats]
    # Returns: (min_x, max_x, range_x, min_y, max_y, range_y,
    #           min_z, max_z, range_z, diagonal)

cubit.get_total_bounding_box(geometry_type, entity_list) -> tuple[10 floats]
cubit.get_tight_bounding_box(geometry_type, entity_list) -> tuple[15 floats]
```

## Surface Properties

```python
cubit.get_surface_normal(surface_id) -> tuple[float, float, float]
cubit.get_surface_normal_at_coord(surface_id, xyz) -> tuple[float, float, float]
cubit.get_surface_type(surface_id) -> str  # "plane", "cylinder", etc.
cubit.is_surface_planar(surface_id) -> bool
cubit.is_cylinder_surface(surface_id) -> bool
cubit.is_cone_surface(surface_id) -> bool

cubit.get_surface_principal_curvatures(surface_id) -> tuple[float, float]
```

## Curve Properties

```python
cubit.get_curve_type(curve_id) -> str        # "line", "arc", "spline", etc.
cubit.get_curve_radius(curve_id) -> float    # For arcs/circles
cubit.get_curve_center(curve_id) -> tuple[float, float, float]
```

## Topology Navigation

```python
cubit.get_relatives(source_type, source_id, target_type) -> tuple[int]
    # Get related entities (parent/child in topology).
    surfaces = cubit.get_relatives("volume", 1, "surface")
    curves = cubit.get_relatives("surface", 3, "curve")
    vertices = cubit.get_relatives("curve", 5, "vertex")

cubit.get_adjacent_surfaces(geometry_type, entity_id) -> tuple[int]
cubit.get_adjacent_volumes(geometry_type, entity_id) -> tuple[int]
cubit.get_owning_volume(geometry_type, entity_id) -> int
cubit.get_owning_body(geometry_type, entity_id) -> int

cubit.get_common_curve_id(surface_1_id, surface_2_id) -> int
cubit.get_common_vertex_id(curve_1_id, curve_2_id) -> int
```

## Entity State

```python
cubit.is_meshed(geometry_type, entity_id) -> bool
cubit.is_merged(geometry_type, entity_id) -> bool
cubit.is_visible(geometry_type, entity_id) -> bool
cubit.is_virtual(geometry_type, entity_id) -> bool
cubit.is_sheet_body(volume_id) -> bool
cubit.is_periodic(geometry_type, entity_id) -> bool
cubit.is_volume_meshable(volume_id) -> bool
cubit.is_surface_meshable(surface_id) -> bool
```

## Entity Counts

```python
cubit.get_volume_count() -> int
cubit.get_surface_count() -> int
cubit.get_curve_count() -> int
cubit.get_vertex_count() -> int
cubit.get_body_count() -> int
cubit.get_block_count() -> int
cubit.get_nodeset_count() -> int
cubit.get_sideset_count() -> int
```

## Geometry Cleanup Detection

```python
cubit.get_small_curves(volume_ids, mesh_size) -> tuple[int]
cubit.get_small_surfaces(volume_ids, mesh_size) -> tuple[int]
cubit.get_narrow_surfaces(volume_ids, mesh_size) -> tuple[int]
cubit.get_small_volumes(volume_ids, mesh_size) -> tuple[int]
cubit.get_blend_surfaces(volume_ids) -> tuple[int]
cubit.get_cone_surfaces(volume_ids) -> tuple[int]
cubit.get_chamfer_surfaces(volume_ids, thickness) -> list
cubit.get_overlapping_volumes(volume_ids) -> tuple[int]
cubit.get_bad_geometry(volume_ids, ...) -> (bodies, volumes, surfaces, curves)
cubit.get_similar_volumes(volume_ids, tol) -> tuple[int]
cubit.get_similar_surfaces(surface_ids, tol) -> tuple[int]
```

Source: coreform.com/cubit_help - CubitInterface namespace
"""

API_MESH_ACCESS = """
# Mesh Element Access Functions

## Node Coordinates

```python
cubit.get_nodal_coordinates(node_id) -> tuple[float, float, float]
    coord = cubit.get_nodal_coordinates(42)
    x, y, z = coord

cubit.get_closest_node(x, y, z) -> int
    # Find node closest to given coordinates.
    node_id = cubit.get_closest_node(0.0, 0.0, 0.0)
```

## Element Connectivity

```python
cubit.get_connectivity(entity_type, entity_id) -> tuple[int]
    # Get 1st-order (corner) node IDs for an element.
    # TET4: returns 4 nodes, HEX8: returns 8 nodes
    nodes = cubit.get_connectivity("tet", tet_id)

cubit.get_expanded_connectivity(entity_type, entity_id) -> tuple[int]
    # Get ALL node IDs including mid-edge/mid-face nodes.
    # TET10: returns 10 nodes, HEX20: returns 20 nodes, HEX27: returns 27
    # Node ordering follows Exodus II convention.
    nodes = cubit.get_expanded_connectivity("hex", hex_id)
```

## Elements by Region

```python
# Volume elements
cubit.get_volume_hexes(volume_id) -> tuple[int]
cubit.get_volume_tets(volume_id) -> tuple[int]
cubit.get_volume_wedges(volume_id) -> tuple[int]
cubit.get_volume_pyramids(volume_id) -> tuple[int]

# Surface elements
cubit.get_surface_quads(surface_id) -> tuple[int]
cubit.get_surface_tris(surface_id) -> tuple[int]

# Curve elements
cubit.get_curve_edges(curve_id) -> tuple[int]

# Nodes by geometry
cubit.get_volume_nodes(volume_id) -> tuple[int]
cubit.get_surface_nodes(surface_id) -> tuple[int]
cubit.get_curve_nodes(curve_id) -> tuple[int]
cubit.get_vertex_node(vertex_id) -> int
```

## Element Properties

```python
cubit.get_element_type(element_id) -> str
    # Returns type by global element ID: "HEX8", "TET10", "TRI3", etc.

cubit.get_global_element_id(element_type, local_id) -> int
    # Convert from Cubit-type-specific ID to global Exodus element ID.
    global_id = cubit.get_global_element_id("hex", hex_id)

cubit.get_hex_global_element_id(hex_id) -> int
cubit.get_tet_global_element_id(tet_id) -> int

cubit.get_node_global_id(node_id) -> int
cubit.get_element_block(element_id) -> int
    # Get block ID containing this element.
```

## Element Counts

```python
cubit.get_hex_count() -> int
cubit.get_tet_count() -> int
cubit.get_tri_count() -> int
cubit.get_quad_count() -> int
cubit.get_wedge_count() -> int
cubit.get_pyramid_count() -> int
cubit.get_node_count() -> int
cubit.get_element_count() -> int
cubit.get_volume_element_count(volume_id) -> int
cubit.get_surface_element_count(surface_id) -> int
```

## Mesh Edge

```python
cubit.get_mesh_edge_length(edge_id) -> float
```

## Node Relationships

```python
cubit.get_node_faces(node_id) -> tuple[int]
cubit.get_node_tris(node_id) -> tuple[int]
cubit.get_node_edges(node_id) -> tuple[int]
cubit.get_node_exists(node_id) -> bool
cubit.get_element_exists(element_id) -> bool
```

Source: coreform.com/cubit_help - CubitInterface namespace
"""

API_BLOCKS_SETS = """
# Block, Nodeset, Sideset, and Group Functions

## Block Access

```python
cubit.get_block_count() -> int
cubit.get_block_id_list() -> tuple[int]
cubit.get_next_block_id() -> int

# Elements in blocks
cubit.get_block_hexes(block_id) -> tuple[int]
cubit.get_block_tets(block_id) -> tuple[int]
cubit.get_block_tris(block_id) -> tuple[int]
cubit.get_block_faces(block_id) -> tuple[int]
cubit.get_block_wedges(block_id) -> tuple[int]
cubit.get_block_pyramids(block_id) -> tuple[int]
cubit.get_block_edges(block_id) -> tuple[int]
cubit.get_block_nodes(block_id) -> tuple[int]

# Geometry in blocks
cubit.get_block_volumes(block_id) -> tuple[int]
cubit.get_block_surfaces(block_id) -> tuple[int]
cubit.get_block_curves(block_id) -> tuple[int]
cubit.get_block_vertices(block_id) -> tuple[int]

# Block properties
cubit.get_block_element_type(block_id) -> str
    # e.g., "HEX8", "TETRA10", "TRI6", "QUAD9"
cubit.get_valid_block_element_types(block_id) -> list[str]
cubit.get_block_material(block_id) -> int
cubit.get_block_attribute_count(block_id) -> int
cubit.get_block_attribute_value(block_id, index) -> float
cubit.get_block_attribute_name(block_id, index) -> str

# Exodus naming
cubit.get_exodus_entity_name(entity_type, entity_id) -> str
cubit.get_exodus_entity_type(entity_type, entity_id) -> str
cubit.get_exodus_element_count(entity_id, entity_type) -> int
cubit.get_exodus_id(entity_type, entity_id) -> int
```

## Nodeset Access

```python
cubit.get_nodeset_count() -> int
cubit.get_nodeset_id_list() -> tuple[int]
cubit.get_next_nodeset_id() -> int

cubit.get_nodeset_nodes(nodeset_id) -> tuple[int]
cubit.get_nodeset_nodes_inclusive(nodeset_id) -> tuple[int]
cubit.get_nodeset_node_count(nodeset_id) -> int
cubit.get_nodeset_volumes(nodeset_id) -> tuple[int]
cubit.get_nodeset_surfaces(nodeset_id) -> tuple[int]
cubit.get_nodeset_curves(nodeset_id) -> tuple[int]
cubit.get_nodeset_vertices(nodeset_id) -> tuple[int]
```

## Sideset Access

```python
cubit.get_sideset_count() -> int
cubit.get_sideset_id_list() -> tuple[int]
cubit.get_next_sideset_id() -> int

cubit.get_sideset_quads(sideset_id) -> tuple[int]
cubit.get_sideset_surfaces(sideset_id) -> tuple[int]
cubit.get_sideset_curves(sideset_id) -> tuple[int]
cubit.get_sideset_edges(sideset_id) -> tuple[int]
cubit.get_sideset_element_type(sideset_id) -> str
```

## Group Access

```python
cubit.group_names_ids() -> list[tuple[str, int]]
    # Returns list of (name, id) pairs.
    for name, gid in cubit.group_names_ids():
        print(f"Group '{name}' = ID {gid}")

cubit.create_new_group() -> int
cubit.get_next_group_id() -> int
cubit.delete_group(group_id)
cubit.delete_all_groups()

cubit.add_entity_to_group(group_id, entity_id, entity_type)
cubit.add_entities_to_group(group_id, entity_ids, entity_type)
cubit.remove_entity_from_group(group_id, entity_id, entity_type)

cubit.get_group_volumes(group_id) -> tuple[int]
cubit.get_group_surfaces(group_id) -> tuple[int]
cubit.get_group_curves(group_id) -> tuple[int]
cubit.get_group_vertices(group_id) -> tuple[int]
cubit.get_group_hexes(group_id) -> tuple[int]
cubit.get_group_tets(group_id) -> tuple[int]
cubit.get_group_tris(group_id) -> tuple[int]
cubit.get_group_nodes(group_id) -> tuple[int]
cubit.get_group_bodies(group_id) -> tuple[int]
cubit.get_group_groups(group_id) -> tuple[int]
```

## Material Access

```python
cubit.get_material_name(material_id) -> str
cubit.get_material_name_list() -> list[str]
cubit.get_material_property(property_enum, entity_id) -> float
cubit.get_blocks_with_materials() -> list[list[int]]
```

Source: coreform.com/cubit_help - CubitInterface namespace
"""

API_MESH_SETTINGS = """
# Mesh Settings and Quality Functions

## Mesh Scheme and Size

```python
cubit.get_mesh_scheme(geometry_type, entity_id) -> str
    # Returns: "tetmesh", "map", "submap", "pave", "sweep", etc.

cubit.get_mesh_scheme_firmness(geometry_type, entity_id) -> str
    # "soft", "hard"

cubit.get_mesh_size(geometry_type, entity_id) -> float
cubit.get_requested_mesh_size(geometry_type, entity_id) -> float
cubit.get_mesh_intervals(geometry_type, entity_id) -> int
cubit.get_requested_mesh_intervals(geometry_type, entity_id) -> int
cubit.get_auto_size(geometry_type, entity_ids, size) -> float

cubit.get_smooth_scheme(geometry_type, entity_id) -> str
```

## Sizing Estimation

```python
cubit.estimate_curve_mesh_size(curve_id, percent_capture) -> float
cubit.get_element_budget(element_type, entity_ids, auto_factor) -> int
    # Estimate element count for given auto factor.
    count = cubit.get_element_budget("tet", [1], 5)
```

## Quality Assessment

```python
cubit.get_quality_value(mesh_type, mesh_id, metric_name) -> float
    # Metrics: "scaled jacobian", "shape", "aspect ratio",
    #          "jacobian", "condition", "distortion", "stretch"
    q = cubit.get_quality_value("hex", 42, "scaled jacobian")

cubit.get_quality_values(mesh_type, mesh_ids, metric_name) -> tuple[float]
    # Batch query for multiple elements.
    ids = cubit.get_volume_hexes(1)
    qualities = cubit.get_quality_values("hex", list(ids), "scaled jacobian")

cubit.get_quality_stats(entity_type, id_list, metric_name,
    single_threshold, use_low_threshold, low_threshold, high_threshold,
    ...) -> (min, max, mean, std, min_id, max_id, mesh_list, elem_type, group_id)
    # Comprehensive quality statistics.
```

## Tetmesh Settings

```python
cubit.get_tetmesh_parallel() -> bool
cubit.get_tetmesh_growth_factor(volume_id) -> float
cubit.get_tetmesh_proximity_flag(volume_id) -> bool
cubit.get_tetmesh_proximity_layers(volume_id) -> int
cubit.get_tetmesh_optimization_level() -> int
cubit.get_tetmesh_insert_mid_nodes() -> bool
```

## Trimesh Settings

```python
cubit.get_trimesh_surface_gradation() -> float
cubit.get_trimesh_volume_gradation() -> float
cubit.get_trimesh_tiny_edge_length() -> float
```

## Curve Bias Settings

```python
cubit.get_curve_bias_type(curve_id) -> str
cubit.get_curve_bias_fine_size(curve_id) -> float
cubit.get_curve_bias_coarse_size(curve_id) -> float
cubit.get_curve_bias_geometric_factor(curve_id) -> float
cubit.get_curve_bias_start_vertex_id(curve_id) -> int
```

## Sweep Source/Target

```python
cubit.get_source_surfaces(volume_id) -> tuple[int]
cubit.get_target_surfaces(volume_id) -> tuple[int]
```

## Boundary Layer

```python
cubit.get_next_boundary_layer_id() -> int
cubit.is_boundary_layer_id_available(id) -> bool
cubit.get_boundary_layer_algorithm(id) -> str
```

## Node Constraint

```python
cubit.get_node_constraint() -> bool
    # Whether mid-edge nodes are constrained to geometry.
cubit.get_node_constraint_value() -> int
```

Source: coreform.com/cubit_help - CubitInterface namespace
"""

API_ENTITY_CLASSES = """
# Entity Class Methods

Access entity objects via: `cubit.volume(id)`, `cubit.surface(id)`,
`cubit.curve(id)`, `cubit.vertex(id)`, `cubit.body(id)`.

## Common Methods (all GeomEntity classes)

```python
entity = cubit.volume(1)   # or cubit.surface(3), cubit.curve(5), etc.

entity.id() -> int                     # Entity ID
entity.entity_name() -> str            # Get name
entity.entity_name(name)               # Set name (setter overload)
entity.set_entity_name(name)           # Set name
entity.entity_names() -> list[str]     # All names
entity.center_point() -> (x, y, z)     # Center coordinates
entity.bounding_box() -> (x_min, x_max, y_min, y_max, z_min, z_max)
entity.is_meshed() -> bool
entity.mesh()                          # Mesh this entity
entity.remove_mesh()                   # Delete mesh
entity.smooth()                        # Smooth mesh
entity.is_visible() -> bool
entity.set_visible(flag)
entity.is_transparent() -> int
entity.set_transparent(flag)

# Navigate topology
entity.bodies() -> list[Body]
entity.volumes() -> list[Volume]
entity.surfaces() -> list[Surface]
entity.curves() -> list[Curve]
entity.vertices() -> list[Vertex]
```

## Volume-Specific Methods

```python
vol = cubit.volume(1)
vol.volume() -> float              # Volume measure
vol.centroid() -> (x, y, z)        # Centroid coordinates
vol.principal_axes() -> tuple[9]   # 3x3 principal axes
vol.principal_moments() -> (I1, I2, I3)
vol.color() -> (r, g, b, a)
vol.set_color((r, g, b, a))
```

## Surface-Specific Methods

```python
surf = cubit.surface(1)
surf.area() -> float
surf.normal_at(location) -> (nx, ny, nz)
    # Normal vector at a point on the surface.
    n = surf.normal_at((1.0, 2.0, 3.0))

surf.closest_point_trimmed(location) -> (x, y, z)
    # Project point onto surface.

surf.closest_point_along_vector(location, direction) -> (x, y, z)
    # Ray-surface intersection.

surf.is_planar() -> bool
surf.is_cylindrical() -> bool

# Parametric (u,v) operations
surf.get_param_range_U() -> (u_min, u_max)
surf.get_param_range_V() -> (v_min, v_max)
surf.position_from_u_v(u, v) -> (x, y, z)
    # Map (u,v) parameter to 3D position.
surf.u_v_from_position(location) -> (u, v)
    # Map 3D position to (u,v) parameter.

surf.principal_curvatures(point) -> (k1, k2)
    # Principal curvature values at a point.

surf.point_containment(point) -> int
    # 0=outside, 1=inside, 2=on boundary

surf.ordered_loops() -> list[list[Curve]]
    # Get boundary loops as ordered curve lists.

surf.color() -> (r, g, b, a)
surf.set_color((r, g, b, a))
```

## Curve-Specific Methods

```python
crv = cubit.curve(1)
crv.length() -> float
crv.curve_center() -> (x, y, z)  # For arcs: center of circle

# Parameterization
crv.start_param() -> float
crv.end_param() -> float
crv.position_from_u(u) -> (x, y, z)
    # Map parameter u to 3D position.
crv.u_from_position(position) -> float
    # Map 3D position to parameter u.
crv.position_from_fraction(fraction) -> (x, y, z)
    # Position at fraction (0-1) along curve.

# Arc length
crv.length_from_u(u1, u2) -> float
crv.u_from_arc_length(root_param, arc_length) -> float
crv.point_from_arc_length(root_param, arc_length) -> (x, y, z)
crv.fraction_from_arc_length(root_vertex, length) -> float

# Differential geometry
crv.closest_point(point) -> (x, y, z)
crv.closest_point_trimmed(point) -> (x, y, z)
crv.tangent(point) -> (tx, ty, tz)
crv.curvature(point) -> (kx, ky, kz)

crv.is_periodic() -> bool
crv.color() -> (r, g, b, a)
crv.set_color((r, g, b, a))
```

## Vertex-Specific Methods

```python
vtx = cubit.vertex(1)
vtx.coordinates() -> (x, y, z)
    # Get vertex position.
    x, y, z = vtx.coordinates()

vtx.color() -> (r, g, b, a)
vtx.set_color((r, g, b, a))
```

Source: coreform.com/cubit_help - Entity class pages
"""

API_GRAPHICS_SELECTION = """
# Graphics and Selection Functions

## Graphics Control

```python
cubit.flush_graphics()
cubit.reset_camera()

cubit.get_view_distance() -> float
cubit.get_view_at() -> (x, y, z)
cubit.get_view_from() -> (x, y, z)
cubit.get_view_up() -> (x, y, z)

cubit.is_perspective_on() -> bool
cubit.is_mesh_visibility_on() -> bool
cubit.is_geometry_visibility_on() -> bool
cubit.get_rendering_mode() -> int
cubit.set_rendering_mode(mode)
```

## Selection and Highlighting

```python
cubit.highlight(entity_type, entity_id)
cubit.clear_highlight()
cubit.clear_preview()

cubit.get_selected_ids() -> tuple[int]
cubit.get_selected_id(index) -> int
cubit.get_selected_type(index) -> str
cubit.current_selection_count() -> int
cubit.clear_picked_list()

cubit.get_pick_type() -> str
cubit.set_pick_type(pick_type, silent=False)

# Entity label settings
cubit.set_label_type(entity_type, label_flag)
cubit.get_label_type(entity_type) -> int
```

## Journaling

```python
cubit.journal_commands(state)          # Enable/disable journaling
cubit.is_command_journaled() -> bool
cubit.write_to_journal(words)
cubit.get_current_journal_file() -> str
```

## Undo

```python
cubit.get_undo_enabled() -> bool
cubit.number_undo_commands() -> int
cubit.was_last_cmd_undoable() -> bool
```

Source: coreform.com/cubit_help - CubitInterface namespace
"""

API_ADVANCED = """
# Advanced / Specialized Functions

## Merge and Overlap Detection

```python
cubit.get_merge_tolerance() -> float
cubit.get_mergeable_surfaces(volume_ids) -> list[list[int]]
cubit.get_mergeable_curves(volume_ids) -> list[list[int]]
cubit.get_mergeable_vertices(volume_ids) -> list[list[int]]
cubit.estimate_merge_tolerance(volume_ids, accurate=False) -> float

cubit.get_overlapping_volumes(volume_ids) -> tuple[int]
cubit.get_overlapping_surfaces_in_bodies(body_ids) -> list[list[int]]
```

## Geometry Analysis

```python
cubit.evaluate_exterior_angle_at_curve(curve_id, volume_id) -> float
cubit.evaluate_surface_angle_at_vertex(surf_id, vert_id) -> float

cubit.find_floating_volumes(volume_ids) -> list[int]
cubit.find_nonmanifold_curves(volume_ids) -> list[int]
cubit.find_nonmanifold_vertices(volume_ids) -> list[int]

cubit.is_blend_surface(surface_id) -> bool
cubit.is_chamfer_surface(surface_id, thickness) -> bool
cubit.is_cavity_surface(surface_id) -> bool
cubit.is_hole_surface(surface_id, radius) -> bool
cubit.is_narrow_surface(surface_id, mesh_size) -> bool

cubit.get_continuous_surfaces(surface_id, angle_tol) -> tuple[int]
cubit.get_continuous_curves(curve_id, angle_tol) -> tuple[int]
cubit.are_adjacent_surfaces(surface_ids) -> bool
cubit.are_adjacent_curves(curve_ids) -> bool
```

## Solutions / Auto-Repair Suggestions

```python
cubit.get_solutions_for_volumes(vol_id, small_curve, mesh_size) -> list
cubit.get_solutions_for_small_surfaces(surface_id, ...) -> list
cubit.get_solutions_for_small_curves(curve_id, ...) -> list
cubit.get_solutions_for_decomposition(volume_list, angle, ...) -> list
cubit.get_solutions_for_blends(surface_id) -> list
cubit.get_solutions_for_forced_sweepability(volume_id, ...) -> list
cubit.get_solutions_for_imprint_merge(surface_id1, surface_id2) -> list
cubit.get_solutions_for_overlapping_volumes(vol1, vol2, ...) -> list
cubit.get_mesh_error_solutions(error_code) -> list[str]
```

## Sculpt Defaults

```python
cubit.get_dbl_sculpt_default(variable) -> float
cubit.get_int_sculpt_default(variable) -> int
cubit.get_bool_sculpt_default(variable) -> bool
cubit.get_string_sculpt_default(variable) -> str
```

## Coordinate Systems

```python
cubit.get_coordinate_systems_id_list() -> tuple[int]
```

## Measurement

```python
cubit.measure_between_entities(type1, id1, type2, id2) -> list[float]
    # Returns distance and closest points.
cubit.snap_locations_to_geometry(locations, entity_type, entity_id, tol) -> list
cubit.calculate_timestep_estimate(entity_type, entity_ids) -> float
```

## Exodus File Reading

```python
cubit.get_all_exodus_times(filename) -> list[float]
cubit.get_all_exodus_variable_names(filename, variable_type) -> list[str]
cubit.get_block_ids(mesh_file_name) -> tuple[int]
```

Source: coreform.com/cubit_help - CubitInterface namespace
"""


def get_api_reference(topic: str = "all") -> str:
	"""Return Cubit Python API reference by topic."""
	topics = {
		"core": API_CORE,
		"geometry_queries": API_GEOMETRY_QUERIES,
		"mesh_access": API_MESH_ACCESS,
		"blocks_sets": API_BLOCKS_SETS,
		"mesh_settings": API_MESH_SETTINGS,
		"entity_classes": API_ENTITY_CLASSES,
		"graphics_selection": API_GRAPHICS_SELECTION,
		"advanced": API_ADVANCED,
	}

	topic = topic.lower().strip()
	if topic == "all":
		return "\n\n".join(topics.values())
	elif topic in topics:
		return topics[topic]
	else:
		return (
			f"Unknown topic: '{topic}'. "
			f"Available: all, {', '.join(topics.keys())}"
		)
