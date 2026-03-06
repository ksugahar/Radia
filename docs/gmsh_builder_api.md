# GmshBuilder API Reference

GmshBuilder is a Gmsh-based mesh generation library built on the GMSH OCC (OpenCASCADE) kernel. It provides approximately 719 public methods across 16 modules, covering geometry creation, boolean operations, transforms, mesh control, mesh generation, mesh quality evaluation, mesh data access, physical group management, geometry queries, export, entity wrappers, geometry analysis, selection groups, post-processing views, and GMSH options.

All coordinates are in meters. GmshBuilder must be used as a context manager to handle GMSH initialization and finalization.

## Quick Start

```python
from radia.gmsh_builder import GmshBuilder

with GmshBuilder(verbose=False) as cm:
    # Create geometry
    box = cm.add_box([0, 0, 0], [0.1, 0.02, 0.03])

    # Mesh control
    cm.set_mesh_size(0.005)

    # Generate mesh
    cm.generate()

    # Export
    cm.export('output.msh')

    # Or convert directly to Radia objects
    radia_obj = cm.to_radia(mu_r=1000)
```

**Constructor**: `GmshBuilder(model_name='gmsh_model', verbose=True)`

---

## Table of Contents

1. [CoreMixin (core.py)](#1-coremixin-corepy)
2. [GeometryMixin (geometry.py)](#2-geometrymixin-geometrypy)
3. [BooleanMixin (boolean_ops.py)](#3-booleanmixin-boolean_opspy)
4. [TransformMixin (transforms.py)](#4-transformmixin-transformspy)
5. [MeshControlMixin (mesh_control.py)](#5-meshcontrolmixin-mesh_controlpy)
6. [MeshGenerateMixin (mesh_generate.py)](#6-meshgeneratemixin-mesh_generatepy)
7. [MeshQualityMixin (mesh_quality.py)](#7-meshqualitymixin-mesh_qualitypy)
8. [MeshAccessMixin (mesh_access.py)](#8-meshaccessmixin-mesh_accesspy)
9. [PhysicalGroupMixin (physical_groups.py)](#9-physicalgroupmixin-physical_groupspy)
10. [QueryMixin (query.py)](#10-querymixin-querypy)
11. [ExportMixin (export.py)](#11-exportmixin-exportpy)
12. [EntityWrapperMixin (entity_wrappers.py)](#12-entitywrappermixin-entity_wrapperspy)
13. [GeometryAnalysisMixin (geometry_analysis.py)](#13-geometryanalysismixin-geometry_analysispy)
14. [GroupMixin (groups.py)](#14-groupmixin-groupspy)
15. [PostProcessingMixin (post_processing.py)](#15-postprocessingmixin-post_processingpy)
16. [OptionsMixin (options.py)](#16-optionsmixin-optionspy)

---

## 1. CoreMixin (core.py)

Context manager, state management, and volume ID allocation.

| Method | Parameters | Returns | Description |
|--------|-----------|---------|-------------|
| `__enter__` | -- | `self` | Initialize GMSH and enter context manager |
| `__exit__` | exc_type, exc_val, exc_tb | `bool` | Finalize GMSH and exit context manager |
| `reset` | -- | `None` | Reset all geometry and mesh state |
| `compress_ids` | -- | `dict` | Re-number volume IDs to be contiguous from 1 |
| `set_verbose` | verbose: bool | `None` | Enable or disable verbose output |
| `get_verbose` | -- | `bool` | Return current verbose setting |
| `get_volume_count` | -- | `int` | Return number of registered volumes |
| `get_all_volume_ids` | -- | `list[int]` | Return sorted list of all volume IDs |
| `has_volume` | vol_id: int | `bool` | Check if volume ID exists |
| `get_gmsh_volume_tags` | vol_id: int | `list[int]` | Return raw GMSH tags for a volume ID |
| `set_model_name` | name: str | `None` | Set/rename GMSH model name |
| `get_model_name` | -- | `str` | Return current model name |
| `is_initialized` | -- | `bool` | Return True if GMSH is initialized |
| `is_meshed` | -- | `bool` | Return True if mesh has been generated |
| `snapshot` | -- | `dict` | Return dict capturing current state |
| `log` | message: str | `None` | Print message if verbose mode is on |

---

## 2. GeometryMixin (geometry.py)

Geometry creation: primitives, sweeps, imports, and advanced OCC shapes. Cubit equivalents: create brick, cylinder, sphere, cone, torus, wedge, sweep, revolve, loft, fillet, chamfer, import step/iges/brep/stl.

### Solid Primitives

| Method | Parameters | Returns | Description |
|--------|-----------|---------|-------------|
| `add_box` | center: list, size: list | `int` | Create a box from center and dimensions |
| `add_box_at` | x0, y0, z0, dx, dy, dz: float | `int` | Create a box from corner point and dimensions |
| `add_cylinder` | center, axis: str, radius, height | `int` | Create a cylinder along x/y/z axis |
| `add_cylinder_vector` | center, direction: list, radius, height | `int` | Create a cylinder with arbitrary direction vector |
| `add_cylinder_between` | p1: list, p2: list, radius | `int` | Create a cylinder between two arbitrary points |
| `add_sphere` | center: list, radius: float | `int` | Create a sphere |
| `add_cone` | center, axis: str, r_base, r_top, height | `int` | Create a cone or truncated cone |
| `add_cone_vector` | center, direction: list, r_base, r_top, height | `int` | Create a cone with arbitrary direction vector |
| `add_torus` | center, axis: str, R, r | `int` | Create a torus (major and minor radii) |
| `add_wedge` | center: list, size: list, ltx: float | `int` | Create a wedge (0 ltx = triangular) |
| `add_ellipsoid` | center, rx, ry, rz | `int` | Create an ellipsoid by scaling a sphere |
| `add_half_space` | axis: str, position: float | `int` | Create a half-space box for boolean ops |
| `add_hollow_cylinder` | center, axis, r_inner, r_outer, height | `int` | Create a hollow cylinder (tube) |
| `add_hollow_sphere` | center, r_inner, r_outer | `int` | Create a hollow sphere (shell) |
| `add_c_shape` | center, axis, outer_r, inner_r, height, gap_angle | `int` | Create C-shaped magnet with angular gap |
| `add_prism` | base_points, height, axis | `int` | Create prism by extruding a polygon |
| `add_volume` | shell_tags: list | `int` | Create volume from shell(s) |
| `add_pipe` | wire_tag, radius | `int` | Create pipe (tube) along a wire path |
| `add_thick_solid` | vol_id, exclude_face_tags, offset | `int` | Create thick solid (shell) by offsetting faces |

### Sweep / Loft / Revolve

| Method | Parameters | Returns | Description |
|--------|-----------|---------|-------------|
| `extrude` | vol_id, dx, dy, dz, num_elements, heights, recombine | `list[int]` | Extrude surfaces of a volume (structured with num_elements) |
| `extrude_surface` | surf_tag, dx, dy, dz, num_elements, heights, recombine | `list[int]` | Extrude a surface to create volume(s) (structured with num_elements) |
| `revolve` | vol_id, center, axis, angle, num_elements, heights, recombine | `list[int]` | Revolve surfaces of a volume (structured with num_elements) |
| `revolve_surface` | surf_tag, center, axis, angle, num_elements, heights, recombine | `list[int]` | Revolve a surface around an axis (structured with num_elements) |
| `loft` | wire_tags: list | `int` | Create volume by lofting through sections |
| `add_thru_sections` | wire_tags, make_solid, make_ruled | `int` | Create volume by lofting wire cross-sections |
| `add_general_pipe` | section_dimtags, wire_tag | `int` | Create pipe by sweeping section along wire |
| `extrude_rotate` | surf_tags, center, axis, angle, n_layers | `list[int]` | Extrude surfaces with rotation |
| `extrude_along_wire` | surf_tag, wire_tag | `list[int]` | Extrude surface along a wire path |
| `extrude_with_twist` | surf_tag, dx, dy, dz, angle, n_layers | `list[int]` | Extrude surface with twist |
| `revolve_partial` | surf_tag, center, axis, angle | `list[int]` | Revolve surface partially |
| `sweep_surface` | surf_tag, wire_tag, make_solid | `list[int]` | Sweep surface along wire |
| `extrude_to_point` | surf_tag, point | `list[int]` | Extrude surface to a point |

### Low-Level Geometry Primitives

| Method | Parameters | Returns | Description |
|--------|-----------|---------|-------------|
| `add_point` | x, y, z, mesh_size | `int` (GMSH tag) | Create a point (vertex) |
| `add_line` | p1, p2 | `int` (GMSH tag) | Create a line between two points |
| `add_circle_arc` | start, center, end | `int` (GMSH tag) | Create a circle arc |
| `add_circular_arc` | center, axis, radius, start_angle, end_angle | `int` (GMSH tag) | Create a circular arc with angles |
| `add_circle` | center, axis, radius | `int` (GMSH tag) | Create a full circle curve |
| `add_full_ellipse` | center, rx, ry, axis | `int` (GMSH tag) | Create a full ellipse curve |
| `add_ellipse_arc` | start, center, major_point, end | `int` (GMSH tag) | Create an ellipse arc |
| `add_spline` | point_tags: list | `int` (GMSH tag) | Create a spline through points |
| `add_bspline` | point_tags: list | `int` (GMSH tag) | Create a B-spline through points |
| `add_bezier` | point_tags: list | `int` (GMSH tag) | Create a Bezier curve |
| `add_wire` | curve_tags: list | `int` (GMSH tag) | Create a wire (curve loop) |
| `add_curve_loop` | curve_tags: list | `int` (GMSH tag) | Create a curve loop (alias for add_wire) |
| `add_rectangle` | center, dx, dy | `int` (GMSH tag) | Create a rectangular surface |
| `add_disk` | center, rx, ry | `int` (GMSH tag) | Create a disk surface |
| `add_polygon` | points: list, z: float | `int` (GMSH tag) | Create a polygon surface from 2D points |
| `add_plane_surface` | wire_tags: list | `int` (GMSH tag) | Create plane surface from wire(s) |
| `add_surface_filling` | wire_tag | `int` (GMSH tag) | Create surface filling a wire boundary |
| `add_surface_loop` | surface_tags: list | `int` (GMSH tag) | Create a surface loop (shell) |
| `add_offset_curve` | curve_tag, offset, dx, dy, dz | `int` (GMSH tag) | Create curve offset from existing curve |

### Advanced Surfaces

| Method | Parameters | Returns | Description |
|--------|-----------|---------|-------------|
| `add_bspline_surface` | point_tags, n_u, n_v, degree_u, degree_v | `int` (GMSH tag) | B-spline surface from control grid |
| `add_bspline_filling` | wire_tag, type_str | `int` (GMSH tag) | B-spline surface filling a wire boundary |
| `add_bezier_surface` | point_tags, n_u, n_v | `int` (GMSH tag) | Bezier surface from control grid |
| `add_bezier_filling` | wire_tag | `int` (GMSH tag) | Bezier surface filling a wire boundary |
| `add_trimmed_surface` | surf_tag, wire_tags | `int` (GMSH tag) | Create trimmed surface |

### Geometry Modification

| Method | Parameters | Returns | Description |
|--------|-----------|---------|-------------|
| `fillet` | vol_id, curve_ids, radius | `int` | Apply fillet to edges |
| `chamfer` | vol_id, curve_ids, distances | `int` | Apply chamfer to edges |
| `heal_geometry` | vol_id, tolerance | `None` | Heal geometry issues |
| `remove_small_edges` | vol_id, tolerance | `None` | Remove small edges |
| `convert_to_nurbs` | vol_id | `None` | Convert geometry to NURBS |

### Import

| Method | Parameters | Returns | Description |
|--------|-----------|---------|-------------|
| `import_step` | filename: str | `list[int]` | Import STEP file geometry |
| `import_iges` | filename: str | `list[int]` | Import IGES file geometry |
| `import_brep` | filename: str | `list[int]` | Import BREP file geometry |
| `import_stl` | filename: str | `list[int]` | Import STL file |
| `import_mesh` | filename: str | `list[int]` | Import .msh mesh file |
| `import_sat` | filename: str | `list[int]` | Import ACIS SAT file |
| `import_x_t` | filename: str | `list[int]` | Import Parasolid X_T file |
| `import_from_string` | format_str, data, ... | `list[int]` | Import geometry from string |
| `export_iges` | filename: str | `None` | Export geometry to IGES |

### Domain-Specific Shapes

| Method | Parameters | Returns | Description |
|--------|-----------|---------|-------------|
| `add_multi_box` | centers: list, sizes: list | `list[int]` | Create multiple boxes at once |
| `add_coil_helix` | center, axis, radius, pitch, n_turns, wire_radius | `int` | Create helical coil geometry |
| `add_sector` | center, axis, r_inner, r_outer, height, angle | `int` | Create an angular sector volume |
| `add_racetrack` | center, straight_length, radius, height | `int` | Create a racetrack coil shape |
| `add_e_shape` | center, size, slot_depths | `int` | Create E-shaped core geometry |
| `add_u_shape` | center, size, slot_depth, slot_width | `int` | Create U-shaped core geometry |

### Miscellaneous

| Method | Parameters | Returns | Description |
|--------|-----------|---------|-------------|
| `copy_geometry` | vol_id | `int` | Copy a volume (convenience alias) |

---

## 3. BooleanMixin (boolean_ops.py)

Boolean operations and webcut. Cubit equivalents: unite, subtract, intersect, webcut, imprint, merge.

### Basic Boolean Operations

| Method | Parameters | Returns | Description |
|--------|-----------|---------|-------------|
| `fuse` | vol_ids: list | `int` | Boolean union of volumes |
| `cut` | target_id, tool_id | `int` | Boolean subtraction (removes tool) |
| `intersect` | vol_id1, vol_id2 | `int` | Boolean intersection |
| `fragment` | vol_ids: list | `list[int]` | Boolean fragment (split at interfaces) |
| `fragment_tracked` | vol_ids: list | `list[int]` | Fragment with tracking of parent-child |

### Extended Boolean Operations

| Method | Parameters | Returns | Description |
|--------|-----------|---------|-------------|
| `cut_keep_tool` | target_id, tool_id | `int` | Boolean subtraction keeping tool volume |
| `fuse_keep` | vol_ids: list | `int` | Boolean union keeping originals |
| `subtract_list` | target_id, tool_ids: list | `int` | Subtract multiple tools from target |
| `fuse_all` | -- | `int` | Fuse all volumes into one |
| `intersect_list` | vol_ids: list | `int` | Intersect multiple volumes |
| `split_by_volume` | vol_id, tool_id | `list[int]` | Split volume by another volume |
| `chop` | vol_id, tool_ids: list | `list[int]` | Chop volume with multiple tools |
| `separate_volumes` | vol_id | `list[int]` | Separate disconnected volumes |
| `cut_with_plane` | vol_id, point, normal | `list[int]` | Cut volume with an arbitrary plane |

### Webcut Operations

| Method | Parameters | Returns | Description |
|--------|-----------|---------|-------------|
| `webcut_plane` | vol_id, axis, position | `list[int]` | Cut volume with coordinate plane |
| `webcut_cylinder` | vol_id, center, axis, radius | `list[int]` | Cut volume with cylinder |
| `webcut_sphere` | vol_id, center, radius | `list[int]` | Cut volume with sphere |
| `webcut_box` | vol_id, center, size | `list[int]` | Cut volume with box |
| `webcut_cone` | vol_id, center, axis, r_base, r_top, height | `list[int]` | Cut volume with cone |
| `webcut_torus` | vol_id, center, axis, R, r | `list[int]` | Cut volume with torus |
| `webcut_general` | vol_id, tool_dimtags | `list[int]` | Cut volume with arbitrary tool shape |
| `section` | vol_id, axis, position | `list[int]` | Section volume (alias for webcut_plane) |

### Imprint, Merge, and Removal

| Method | Parameters | Returns | Description |
|--------|-----------|---------|-------------|
| `imprint_and_merge` | -- | `None` | Imprint and merge all volumes |
| `imprint_volumes` | vol_ids: list | `None` | Imprint specific volumes |
| `merge_volumes` | vol_ids, tolerance | `None` | Merge specific volumes |
| `remove_volume` | vol_id | `None` | Remove a volume |
| `remove_volumes` | vol_ids: list | `None` | Remove multiple volumes |
| `remove_all_duplicates` | -- | `None` | Remove duplicate entities |
| `keep_only` | vol_ids: list | `None` | Keep only specified volumes, remove others |

---

## 4. TransformMixin (transforms.py)

Transform operations. Cubit equivalents: move, rotate, copy, reflect, scale, pattern.

### Basic Transforms

| Method | Parameters | Returns | Description |
|--------|-----------|---------|-------------|
| `translate` | vol_id, dx, dy, dz | `None` | Translate a volume |
| `rotate` | vol_id, center, axis, angle | `None` | Rotate a volume (radians) |
| `rotate_degrees` | vol_id, center, axis, angle_deg | `None` | Rotate a volume (degrees) |
| `copy` | vol_id | `int` | Copy a volume |
| `mirror` | vol_id, plane_axis: str | `None` | Mirror/reflect about x/y/z plane |
| `mirror_plane` | vol_id, point, normal | `None` | Mirror about an arbitrary plane |
| `scale` | vol_id, sx, sy, sz, center | `None` | Scale a volume |
| `affine_transform` | vol_id, matrix: list | `None` | Apply 4x4 affine transformation |
| `translate_to` | vol_id, target_point | `None` | Move volume center to target point |
| `center_at_origin` | vol_id | `None` | Move volume center to origin |
| `align_to_axis` | vol_id, from_axis, to_axis | `None` | Rotate to align direction vectors |

### Copy-Transform Combos

| Method | Parameters | Returns | Description |
|--------|-----------|---------|-------------|
| `copy_translate` | vol_id, dx, dy, dz | `int` | Copy and translate in one step |
| `copy_rotate` | vol_id, center, axis, angle | `int` | Copy and rotate in one step |
| `copy_mirror` | vol_id, plane_axis | `int` | Copy and mirror in one step |
| `copy_scale` | vol_id, factor, center | `int` | Copy and scale in one step |
| `symmetrize` | vol_id, plane_axis | `int` | Mirror copy and fuse to make symmetric |

### Array / Pattern

| Method | Parameters | Returns | Description |
|--------|-----------|---------|-------------|
| `array_linear` | vol_id, direction, n, spacing | `list[int]` | Create linear array of copies |
| `array_circular` | vol_id, center, axis, n, angle | `list[int]` | Create circular array of copies |
| `array_grid` | vol_id, nx, ny, nz, sx, sy, sz | `list[int]` | Create 3D grid array of copies |
| `array_along_curve` | vol_id, curve_tag, n | `list[int]` | Distribute copies along a curve |

### Bulk Transforms

| Method | Parameters | Returns | Description |
|--------|-----------|---------|-------------|
| `translate_all` | dx, dy, dz | `None` | Translate all volumes |
| `rotate_all` | center, axis, angle | `None` | Rotate all volumes |
| `mirror_all` | plane_axis | `None` | Mirror all volumes |
| `scale_all` | factor, center | `None` | Scale all volumes |
| `copy_all` | -- | `dict` | Copy all volumes, return {old: new} map |
| `swap_volumes` | vol_id1, vol_id2 | `None` | Swap GMSH tags between two volume IDs |

---

## 5. MeshControlMixin (mesh_control.py)

Mesh size, scheme, bias, boundary layer, and sizing fields. Cubit equivalents: size, interval, scheme, bias, boundary layer.

### Size and Division Control

| Method | Parameters | Returns | Description |
|--------|-----------|---------|-------------|
| `set_mesh_size` | size: float | `None` | Set global mesh element size |
| `set_mesh_size_at` | vol_id, size | `None` | Set mesh size at a specific volume |
| `set_mesh_size_on_surface` | surf_tag, size | `None` | Set mesh size on a surface |
| `set_mesh_size_on_curve` | curve_tag, size | `None` | Set mesh size on a curve |
| `set_mesh_size_on_point` | point_tag, size | `None` | Set mesh size at a point |
| `set_mesh_size_factor` | factor: float | `None` | Set global mesh size factor |
| `set_mesh_size_from_boundary` | extend: int | `None` | Control size extension from boundary |
| `set_min_max_size` | min_size, max_size | `None` | Set global min/max element size |
| `set_divisions` | vol_id, nx, ny, nz | `None` | Set transfinite divisions for a volume |
| `set_divisions_all` | nx, ny, nz | `None` | Set transfinite divisions for all volumes |

### Element Order and Scheme

| Method | Parameters | Returns | Description |
|--------|-----------|---------|-------------|
| `set_order` | order: int | `None` | Set mesh element order (1 or 2) |
| `set_scheme` | vol_id, scheme: str | `None` | Set meshing scheme for a volume |
| `set_algorithm` | algo_2d, algo_3d | `None` | Set meshing algorithm |
| `set_algorithm_on_surface` | surf_tag, algo | `None` | Set algorithm for a specific surface |
| `set_algorithm_on_volume` | vol_tag, algo | `None` | Set algorithm for a specific volume |

### Curve Control

| Method | Parameters | Returns | Description |
|--------|-----------|---------|-------------|
| `set_curve_divisions` | curve_tag, n, mesh_type, coef | `None` | Set curve divisions with distribution |
| `set_bias` | curve_tag, n, ratio | `None` | Set biased mesh grading on a curve |
| `set_double_bias` | curve_tag, n, ratio | `None` | Set double-biased grading on a curve |

### Transfinite Control

| Method | Parameters | Returns | Description |
|--------|-----------|---------|-------------|
| `set_transfinite_surface` | surf_tag, arrangement, corners | `None` | Set transfinite surface meshing |
| `set_transfinite_volume` | vol_tag, corners | `None` | Set transfinite volume meshing |
| `set_transfinite_automatic` | vol_tag, corner_angle, recombine | `None` | Automatic transfinite setup |
| `set_transfinite_automatic_all` | corner_angle, recombine | `None` | Automatic transfinite for all volumes |

### Boundary Layer

| Method | Parameters | Returns | Description |
|--------|-----------|---------|-------------|
| `set_boundary_layer` | surf_tags, thickness, n_layers, ratio | `None` | Set boundary layer mesh |

### Sizing Fields

| Method | Parameters | Returns | Description |
|--------|-----------|---------|-------------|
| `add_field_distance` | entity_tags, dim | `int` | Add distance field |
| `add_field_threshold` | distance_field, dist_min, dist_max, size_min, size_max | `int` | Add threshold field |
| `add_field_box` | xmin, xmax, ymin, ymax, zmin, zmax, size_in, size_out | `int` | Add box sizing field |
| `add_field_ball` | xc, yc, zc, radius, size_in, size_out | `int` | Add ball sizing field |
| `add_field_cylinder` | xc, yc, zc, xa, ya, za, radius, size_in, size_out | `int` | Add cylinder sizing field |
| `add_field_frustum` | x1, y1, z1, x2, y2, z2, r1_in, r1_out, ... | `int` | Add frustum sizing field |
| `add_field_math` | expression: str | `int` | Add math expression field |
| `add_field_min` | field_tags: list | `int` | Add min of fields |
| `add_field_max` | field_tags: list | `int` | Add max of fields |
| `add_field_mean` | field_tags: list | `int` | Add mean of fields |
| `add_field_restrict` | field_tag, vol_ids | `int` | Restrict field to volumes |
| `add_field_attractor` | entity_tags, dim, n_per_radius | `int` | Add attractor field |
| `add_field_laplacian` | field_tag | `int` | Add Laplacian smoothed field |
| `add_field_gradient` | field_tag | `int` | Add gradient field |
| `add_field_curvature_field` | -- | `int` | Add curvature-based field |
| `add_field_octree` | size_in, size_out | `int` | Add octree sizing field |
| `add_field_external_process` | command: str | `int` | Add external process field |
| `add_field_structured` | filename: str | `int` | Add structured field from file |
| `add_field_post_view` | view_tag | `int` | Add field from post-processing view |
| `add_field_param` | field_tag, key, value | `None` | Set field parameter |
| `set_background_field` | field_tag | `None` | Set active background field |
| `remove_field` | field_tag | `None` | Remove a sizing field |
| `remove_all_fields` | -- | `None` | Remove all sizing fields |
| `get_field_count` | -- | `int` | Return number of fields |
| `get_field_type` | field_tag | `str` | Return field type string |
| `get_field_list` | -- | `list` | Return all fields info |

### CFD Convenience Methods

| Method | Parameters | Returns | Description |
|--------|-----------|---------|-------------|
| `add_size_gradient` | entity_tags, dim, lc_near, lc_far, dist_near, dist_far, set_as_background | `int` | Distance+Threshold in one step (most common CFD sizing pattern) |
| `set_background_mesh_from_fields` | field_tags, operator ('min'/'max'/'mean') | `int` | Combine fields and set as background mesh |
| `set_field_numbers` | field_tag, key, values (list) | `None` | Set list-valued field parameter (setNumbers) |
| `add_boundary_layer_field` | entity_tags, dim, lc_min, lc_max, dist_min, dist_max, n_layers, ratio, fan_points, intersect_metrics | `int` | BoundaryLayer field for CFD wall resolution |
| `add_field_auto_mesh_size` | -- | `int` | Automatic geometry-based sizing (GMSH>=4.11) |
| `set_mesh_size_at_points_direct` | point_tags (list), size | `None` | Batch mesh size at multiple points |
| `add_wake_refinement` | curve_tags, direction, length, lc_min, lc_max, width | `int` | Anisotropic wake region refinement (airfoil/bluff body) |
| `add_terrain_field` | filename, format ('structured'/'post_view') | `int` | Terrain/elevation data as sizing field |

### Mesh Optimization and Recombination

| Method | Parameters | Returns | Description |
|--------|-----------|---------|-------------|
| `set_curvature_mesh` | n_per_2pi: int | `None` | Set curvature-based mesh sizing |
| `set_recombine` | surf_tag, angle (default 45.0) | `None` | Set recombine on a surface (tri -> quad) |
| `set_recombine_all` | -- | `None` | Set recombine on all surfaces |
| `set_recombination_algorithm` | algo: int or str | `None` | Set recombination algorithm (0=simple, 1=blossom, 2=full_quad, 3=blossom_full_quad) |
| `set_smoothing` | n_steps: int | `None` | Set global smoothing steps |
| `set_smoothing_on_surface` | surf_tag, n_steps | `None` | Set smoothing on a surface |
| `set_smoothing_on_volume` | vol_tag, n_steps | `None` | Set smoothing on a volume |
| `set_optimization` | method: str | `None` | Set optimization method |
| `set_growth_rate` | rate: float | `None` | Set mesh growth rate |
| `set_angle_threshold` | angle_deg: float | `None` | Set angle threshold for mesh |
| `set_subdivision_algorithm` | algo: int | `None` | Set subdivision algorithm |
| `get_mesh_options` | -- | `dict` | Return current mesh options |

### Compound and Misc

| Method | Parameters | Returns | Description |
|--------|-----------|---------|-------------|
| `set_compound` | dim, tags | `None` | Set compound meshing |
| `set_compound_surfaces` | surf_tags: list | `None` | Set compound surface meshing |
| `set_compound_curves` | curve_tags: list | `None` | Set compound curve meshing |
| `set_reverse_mesh` | dim, tag | `None` | Reverse mesh element orientation |
| `set_size_callback` | callback | `None` | Set size callback function |
| `remove_size_callback` | -- | `None` | Remove size callback |
| `set_size_at_parametric_points` | dim, tag, parametric_coords, sizes | `None` | Set sizes at parametric points |
| `set_mesh_visibility` | dim, tag, visible | `None` | Set mesh visibility |
| `remove_mesh_constraints` | dim, tag | `None` | Remove mesh constraints |
| `set_outward_orientation` | tag | `None` | Set outward surface orientation |
| `remove_embedded` | dim, tag | `None` | Remove embedded entities |
| `get_mesh_sizes` | -- | `dict` | Return current mesh size info |
| `get_embedded_entities` | dim, tag | `list` | Get embedded entities |

---

## 6. MeshGenerateMixin (mesh_generate.py)

Mesh generation and optimization. Cubit equivalents: mesh vol, mesh surface, delete mesh, refine, smooth.

### Core Generation

| Method | Parameters | Returns | Description |
|--------|-----------|---------|-------------|
| `generate` | element_type: str ('hex'/'tet') | `None` | Generate 3D mesh for all volumes |
| `generate_surface` | mesh_size | `None` | Generate 2D surface mesh |
| `generate_curve` | curve_tags | `None` | Generate 1D curve mesh |
| `generate_volume` | vol_ids, element_type | `None` | Generate mesh for specific volumes |
| `generate_tet` | mesh_size | `None` | Generate tetrahedral mesh |
| `generate_hex` | mesh_size | `None` | Generate hexahedral mesh |
| `generate_mixed` | -- | `None` | Generate mixed hex/tet mesh |
| `generate_surface_quad` | mesh_size | `None` | Generate quad surface mesh |
| `generate_and_optimize` | element_type, method | `None` | Generate and immediately optimize |
| `generate_adaptive` | field_tag, n_iter, target_quality | `None` | Generate with adaptive refinement |
| `generate_boundary_layer_mesh` | surf_tags, thickness, n_layers, ratio | `None` | Generate boundary layer mesh |
| `generate_prism_layers` | surf_tags, n_layers, total_thickness, ratio | `None` | Generate prism layer mesh |
| `generate_surface_on_entity` | surf_tag | `None` | Generate surface mesh on single entity |
| `generate_curve_on_entity` | curve_tag | `None` | Generate curve mesh on single entity |
| `import_stl_mesh` | filename, angle | `None` | Import and mesh an STL file |

### Mesh Modification

| Method | Parameters | Returns | Description |
|--------|-----------|---------|-------------|
| `delete_mesh` | vol_id | `None` | Delete mesh for a volume |
| `delete_all_mesh` | -- | `None` | Delete all meshes |
| `refine` | n_refine: int | `None` | Refine mesh n times |
| `optimize_mesh` | method: str | `None` | Optimize mesh quality |
| `smooth_laplacian` | n_iterations: int | `None` | Laplacian smoothing |
| `convert_to_second_order` | -- | `None` | Convert mesh to 2nd order |
| `convert_to_first_order` | -- | `None` | Convert mesh to 1st order |
| `recombine_mesh` | -- | `None` | Recombine triangles to quads |
| `split_quadrangles` | -- | `None` | Split quads into triangles |
| `reverse_elements` | dim, tag | `None` | Reverse element orientations |
| `partition_mesh` | n_parts: int | `None` | Partition mesh into n parts |
| `unpartition_mesh` | -- | `None` | Remove mesh partitioning |

### Node and Element Management

| Method | Parameters | Returns | Description |
|--------|-----------|---------|-------------|
| `renumber_nodes` | -- | `None` | Renumber nodes contiguously |
| `renumber_elements` | -- | `None` | Renumber elements contiguously |
| `reorder_elements` | method, dim, tag | `None` | Reorder elements (e.g. RCMK) |
| `compute_renumbering` | method, dim, tag | `list` | Compute renumbering without applying |
| `reclassify_nodes` | -- | `None` | Reclassify nodes on boundaries |
| `relocate_nodes` | dim, tag | `None` | Relocate nodes to geometry |
| `set_mesh_coherence` | value: bool | `None` | Set mesh coherence mode |
| `rebuild_node_cache` | -- | `None` | Rebuild internal node cache |
| `rebuild_element_cache` | -- | `None` | Rebuild internal element cache |

### Embedding

| Method | Parameters | Returns | Description |
|--------|-----------|---------|-------------|
| `embed_point_in_volume` | point_tag, vol_tag | `None` | Embed point in volume mesh |
| `embed_curve_in_surface` | curve_tag, surf_tag | `None` | Embed curve in surface mesh |
| `embed_surface_in_volume` | surf_tag, vol_tag | `None` | Embed surface in volume mesh |
| `embed_entities` | dim, tags, in_dim, in_tag | `None` | General-purpose entity embedding |

### Topology and Homology

| Method | Parameters | Returns | Description |
|--------|-----------|---------|-------------|
| `classify_surfaces` | -- | `None` | Classify surface elements |
| `classify_surfaces_parametric` | angle, boundary, for_reparametrize, curve_angle | `None` | Classify surfaces with configurable thresholds |
| `adapt_mesh_iterative` | field_tag, n_iterations, size_factor | `None` | Field-driven iterative mesh adaptation |
| `create_geometry_from_mesh` | -- | `None` | Create geometry from mesh |
| `create_topology` | make_simply_connected | `None` | Create topology from mesh |
| `add_homology_request` | type_str, domain_tags, subdomain_tags, ... | `None` | Add homology computation request |
| `compute_homology` | -- | `None` | Compute homology |
| `compute_cross_field` | -- | `None` | Compute cross field |
| `clear_homology_requests` | -- | `None` | Clear homology requests |
| `get_mesh_statistics` | -- | `dict` | Return detailed mesh statistics |
| `get_ghost_elements` | dim, tag | `list` | Get ghost elements for partition |
| `get_partition_count` | -- | `int` | Get number of mesh partitions |
| `create_edges` | -- | `None` | Create edge data structures |
| `create_faces` | order | `None` | Create face data structures |

---

## 7. MeshQualityMixin (mesh_quality.py)

Mesh quality evaluation. Cubit equivalents: quality vol, quality histogram, scaled jacobian.

| Method | Parameters | Returns | Description |
|--------|-----------|---------|-------------|
| `get_quality` | metric: str | `np.array` | Get quality array for all elements |
| `get_quality_stats` | metric: str | `dict` | Get quality statistics (min/max/mean/std) |
| `get_quality_histogram` | metric, n_bins | `dict` | Get quality histogram |
| `get_scaled_jacobian` | -- | `np.array` | Get scaled Jacobian for all elements |
| `get_aspect_ratio` | -- | `np.array` | Get aspect ratio for all elements |
| `get_min_quality` | metric: str | `float` | Get minimum quality value |
| `get_max_quality` | metric: str | `float` | Get maximum quality value |
| `get_quality_range` | metric: str | `tuple` | Return (min, max) quality range |
| `get_quality_percentile` | metric, percentile | `float` | Quality value at given percentile |
| `get_worst_elements` | metric, n | `list[dict]` | Get n worst elements by quality |
| `get_element_quality` | elem_tag, metric | `float` | Get quality for a single element |
| `get_quality_by_type` | metric, elem_type | `np.array` | Get quality for specific element type |
| `get_element_quality_map` | metric: str | `dict` | Map of element tag to quality value |
| `get_volume_quality` | vol_id, metric | `dict` | Quality stats for a specific volume |
| `get_quality_summary` | -- | `dict` | Multiple quality metrics summary |
| `get_distortion` | -- | `dict` | Get mesh distortion metrics |
| `get_edge_length_range` | -- | `tuple` | Return (min, max) edge lengths |
| `get_volume_ratio_range` | -- | `tuple` | Return (min, max, ratio) element volumes |
| `get_face_quality` | metric: str | `np.array` | Get quality of 2D (surface) elements |
| `count_elements_below_threshold` | metric, threshold | `int` | Count poor quality elements |
| `check_inverted_elements` | -- | `dict` | Check for negative Jacobian elements |
| `improve_quality` | method: str | `None` | Improve mesh quality |
| `fix_negative_jacobians` | -- | `int` | Fix inverted elements |
| `smooth_mesh` | n_iterations, method | `None` | Smooth mesh with repeated optimization |
| `untangle_mesh` | -- | `int` | Fix tangled (inverted) elements |
| `optimize_high_order` | -- | `None` | Optimize high-order element quality |
| `print_quality_report` | metric: str | `None` | Print formatted quality report |

---

## 8. MeshAccessMixin (mesh_access.py)

Mesh data access: nodes, elements, connectivity. Cubit equivalents: get_node_count, get_connectivity, get_nodal_coordinates.

### Node Access

| Method | Parameters | Returns | Description |
|--------|-----------|---------|-------------|
| `get_nodes` | -- | `np.array` | Get all node tags |
| `get_node_count` | -- | `int` | Get total node count |
| `get_node_coordinates` | node_tag: int | `list` | Get [x,y,z] of a node |
| `get_all_node_coordinates` | -- | `np.array` | Get all node coords as (N,3) array |
| `get_node_coordinates_batch` | node_tags: list | `np.array` | Get coords for multiple nodes |
| `get_closest_node` | x, y, z | `dict` | Find node closest to a point |
| `get_nodes_in_box` | bbox: dict | `list` | Get node tags within bounding box |
| `get_boundary_nodes` | -- | `list` | Get nodes on model boundary |
| `get_nodes_on_surface` | surf_tag | `list` | Get nodes on a surface |
| `get_nodes_on_curve` | curve_tag | `list` | Get nodes on a curve |
| `get_internal_nodes` | -- | `list` | Get internal (non-boundary) nodes |
| `get_volume_node_count` | vol_id | `int` | Get node count for a volume |
| `get_surface_node_count` | surf_tag | `int` | Get node count for a surface |
| `get_duplicate_nodes` | -- | `list` | Find duplicate nodes |
| `set_node_coordinates` | node_tag, x, y, z | `None` | Set coordinates of a node |
| `add_nodes_manual` | dim, tag, coords, tags | `None` | Manually add nodes |
| `remove_duplicate_nodes` | tolerance | `None` | Remove duplicate nodes |
| `relocate_nodes` | dim, tag | `None` | Relocate nodes to geometry |
| `get_nodes_for_physical_group` | dim, tag | `list` | Get nodes for a physical group |
| `get_nodes_by_element_type` | elem_type, tag | `list` | Get nodes for element type |

### Element Access

| Method | Parameters | Returns | Description |
|--------|-----------|---------|-------------|
| `get_elements` | dim: int | `dict` | Get all elements for a dimension |
| `get_element_nodes` | elem_tag | `list` | Get node tags of an element |
| `get_element_type` | elem_tag | `str` | Get element type name |
| `get_element_centroid` | elem_tag | `list` | Get centroid [x,y,z] of element |
| `get_element_volume` | elem_tag | `float` | Get volume of a 3D element |
| `get_element_face_count` | elem_tag | `int` | Get face count for an element |
| `get_hex_count` | -- | `int` | Count hexahedral elements |
| `get_tet_count` | -- | `int` | Count tetrahedral elements |
| `get_tri_count` | -- | `int` | Count triangular elements |
| `get_quad_count` | -- | `int` | Count quadrilateral elements |
| `get_wedge_count` | -- | `int` | Count wedge elements |
| `get_pyramid_count` | -- | `int` | Count pyramid elements |
| `get_prism_count` | -- | `int` | Count prism elements |
| `get_element_count` | dim: int | `int` | Get total element count for dimension |
| `get_volume_elements` | vol_id | `dict` | Get elements for a specific volume |
| `get_surface_elements` | surf_tag | `dict` | Get elements for a surface |
| `get_volume_element_count` | vol_id | `int` | Element count for a volume |
| `get_volume_hex_count` | vol_id | `int` | Hex count for a volume |
| `get_volume_tet_count` | vol_id | `int` | Tet count for a volume |
| `get_elements_by_type` | elem_type_name, dim | `list` | Get elements by type name |
| `get_elements_by_type_raw` | elem_type, tag | `tuple` | Get raw element data by type |
| `get_elements_in_sphere` | center, radius, dim | `list` | Get elements within a sphere |
| `get_element_by_coordinates` | x, y, z, dim | `int` | Find element at a point |
| `get_elements_by_coordinates` | x, y, z, dim | `list` | Find elements near a point |
| `get_element_properties` | elem_type | `dict` | Get element type properties |
| `get_max_element_tag` | -- | `int` | Get maximum element tag |
| `get_max_node_tag` | -- | `int` | Get maximum node tag |
| `get_element_types` | dim, tag | `list` | Get element types present |
| `add_elements_manual` | dim, tag, elem_type, node_tags, tags | `None` | Manually add elements |
| `remove_duplicate_elements` | -- | `None` | Remove duplicate elements |

### Connectivity and Topology

| Method | Parameters | Returns | Description |
|--------|-----------|---------|-------------|
| `get_elements_as_numpy` | dim: int | `dict` | Get elements as numpy arrays |
| `get_node_to_element_map` | -- | `dict` | Map node tag to element tags |
| `get_element_neighbors` | elem_tag | `list` | Get neighboring elements |
| `get_edge_nodes` | -- | `list` | Get edge node pairs |
| `get_face_nodes` | dim: int | `list` | Get face node lists |
| `get_connectivity_matrix` | dim: int | `np.array` | Get connectivity matrix |
| `get_element_edge_nodes` | elem_type, tag, primary | `list` | Get edge nodes for element type |
| `get_element_face_nodes` | elem_type, tag, primary | `list` | Get face nodes for element type |
| `get_all_edges` | -- | `list` | Get all edges |
| `get_all_faces` | order | `list` | Get all faces |
| `get_edges_for_element` | elem_tag | `list` | Get edges for a specific element |

### FEM-Related

| Method | Parameters | Returns | Description |
|--------|-----------|---------|-------------|
| `get_local_coords_in_element` | elem_tag, x, y, z | `list` | Get local coords in element |
| `get_jacobian` | elem_type, local_coords, tag | `np.array` | Get Jacobian at a point |
| `get_jacobians` | elem_type, local_coords, tag | `tuple` | Get Jacobians for multiple elements |
| `get_basis_functions` | elem_type, local_coords, function_type | `tuple` | Get basis function values |
| `get_integration_points` | elem_type, order | `tuple` | Get Gauss integration points |
| `get_element_barycenters` | elem_type, tag, fast, primary | `np.array` | Get element barycenters |
| `get_element_sizes` | dim, tag | `np.array` | Get element sizes |
| `get_element_type_for_name` | name: str | `int` | Get GMSH element type for name |
| `get_nodes_per_element` | elem_type: int | `int` | Get nodes per element for type |

### Periodic

| Method | Parameters | Returns | Description |
|--------|-----------|---------|-------------|
| `set_periodic_mesh` | dim, tags, source_tags, affine_transform | `None` | Set periodic mesh constraint |
| `get_periodic_nodes` | dim, tag | `dict` | Get periodic node mapping |
| `get_periodic_keys` | elem_type, function_type | `list` | Get periodic keys |

---

## 9. PhysicalGroupMixin (physical_groups.py)

Physical group, block, sideset, nodeset, and material management. Cubit equivalents: block, sideset, nodeset, material assignment.

### Physical Groups

| Method | Parameters | Returns | Description |
|--------|-----------|---------|-------------|
| `set_physical_group` | vol_id, name | `None` | Assign physical group to volume |
| `get_physical_groups` | dim | `list` | Get all physical groups |
| `get_physical_group_name` | tag, dim | `str` | Get name of physical group |
| `get_physical_group_entities` | tag, dim | `list` | Get entities in physical group |
| `remove_physical_group` | tag, dim | `None` | Remove a physical group |
| `remove_all_physical_groups` | -- | `None` | Remove all physical groups |
| `rename_physical_group` | tag, dim, new_name | `None` | Rename a physical group |

### Blocks (3D)

| Method | Parameters | Returns | Description |
|--------|-----------|---------|-------------|
| `add_block` | vol_ids, block_id, name | `int` | Add a block |
| `add_block_by_name` | vol_ids, name | `int` | Add a block with required name |
| `auto_assign_blocks` | -- | `None` | Auto-assign blocks for all volumes |
| `get_block_count` | -- | `int` | Count blocks |
| `get_all_blocks` | -- | `list[dict]` | Get info for all blocks |
| `get_block_name` | block_tag | `str` | Get name of a block |
| `rename_block` | block_tag, new_name | `None` | Rename a block |
| `get_block_elements` | block_tag | `dict` | Get elements in block |
| `get_block_hex_count` | block_tag | `int` | Count hexes in block |
| `get_block_tet_count` | block_tag | `int` | Count tets in block |
| `get_block_node_count` | block_tag | `int` | Count nodes in block |
| `get_block_nodes` | block_tag | `list` | Get node tags in block |
| `get_block_element_count` | block_tag | `int` | Count elements in block |
| `get_block_volumes` | block_tag | `list` | Get volume tags in block |
| `get_block_surfaces` | block_tag | `list` | Get surface tags in block |
| `get_block_bounding_box` | block_tag | `dict` | Get block bounding box |
| `get_block_center` | block_tag | `list` | Get block center |
| `set_block_attribute` | block_tag, key, value | `None` | Set block attribute |
| `get_block_attribute` | block_tag, key | value | Get block attribute |
| `get_block_attribute_names` | block_tag | `list` | Get block attribute names |
| `remove_block_attribute` | block_tag, key | `None` | Remove block attribute |

### Sidesets (2D)

| Method | Parameters | Returns | Description |
|--------|-----------|---------|-------------|
| `add_sideset` | surf_tags, sideset_id, name | `int` | Add a sideset |
| `add_sideset_by_name` | surf_tags, name | `int` | Add sideset with name |
| `auto_assign_sidesets` | -- | `None` | Auto-assign sidesets |
| `get_sideset_count` | -- | `int` | Count sidesets |
| `get_all_sidesets` | -- | `list[dict]` | Get all sidesets info |
| `get_sideset_name` | ss_tag | `str` | Get name of a sideset |
| `get_sideset_elements` | ss_tag | `dict` | Get elements in sideset |
| `get_sideset_tri_count` | ss_tag | `int` | Count triangles in sideset |
| `get_sideset_quad_count` | ss_tag | `int` | Count quads in sideset |
| `get_sideset_node_count` | ss_tag | `int` | Count nodes in sideset |
| `get_sideset_nodes` | ss_tag | `list` | Get node tags in sideset |
| `get_sideset_surfaces` | ss_tag | `list` | Get surface tags in sideset |
| `get_sideset_area` | ss_tag | `float` | Get total sideset area |
| `get_sideset_bounding_box` | ss_tag | `dict` | Get sideset bounding box |

### Nodesets (0D)

| Method | Parameters | Returns | Description |
|--------|-----------|---------|-------------|
| `add_nodeset` | point_tags, nodeset_id, name | `int` | Add a nodeset |
| `add_nodeset_by_name` | point_tags, name | `int` | Add nodeset with name |
| `get_nodeset_count` | -- | `int` | Count nodesets |
| `get_all_nodesets` | -- | `list[dict]` | Get info for all nodesets |
| `get_nodeset_name` | ns_tag | `str` | Get name of a nodeset |
| `get_nodeset_nodes` | ns_tag | `list` | Get node tags in nodeset |
| `get_nodeset_node_count` | ns_tag | `int` | Count nodes in nodeset |
| `get_nodeset_coordinates` | ns_tag | `list` | Get coords of nodeset nodes |
| `get_nodeset_bounding_box` | ns_tag | `dict` | Get nodeset bounding box |

### Materials

| Method | Parameters | Returns | Description |
|--------|-----------|---------|-------------|
| `set_material` | vol_id, name, mu_r | `None` | Set material for a volume |
| `set_material_by_name` | vol_ids, name, properties | `None` | Set material by name with properties |
| `set_material_conductivity` | vol_id, sigma | `None` | Set material conductivity |
| `set_material_density` | vol_id, density | `None` | Set material density |
| `get_material` | vol_id | `dict` | Get material info for volume |
| `get_material_property` | vol_id, key | value | Get specific material property |
| `get_materials` | -- | `dict` | Get all materials |
| `has_material` | vol_id | `bool` | Check if volume has material |
| `clear_materials` | -- | `None` | Clear all materials |
| `copy_material` | from_vol, to_vol | `None` | Copy material from one volume to another |
| `export_material_table` | -- | `dict` | Export material table |

### Misc

| Method | Parameters | Returns | Description |
|--------|-----------|---------|-------------|
| `assign_boundary_names` | mapping: dict | `None` | Assign boundary names from mapping |
| `auto_assign_all` | -- | `None` | Auto-assign blocks, sidesets, nodesets |

---

## 10. QueryMixin (query.py)

Geometry and topology queries. Cubit equivalents: get_bounding_box, get_center, get_volume, list entities.

### Volume Info

| Method | Parameters | Returns | Description |
|--------|-----------|---------|-------------|
| `get_volumes` | -- | `list[int]` | Get all managed volume IDs |
| `get_gmsh_tags` | vol_id | `list[int]` | Get underlying GMSH volume tags |
| `get_mesh_info` | -- | `dict` | Get mesh statistics |
| `get_bounding_box` | vol_id | `dict` | Get axis-aligned bounding box |
| `get_center` | vol_id | `list` | Get center [x,y,z] of volume |
| `get_center_of_mass` | vol_id | `list` | Get center of mass |
| `get_volume_measure` | vol_id | `float` | Get volume measure (m^3) |
| `get_mass` | vol_id, density | `float` | Get mass (volume * density) |
| `get_inertia` | vol_id | `list` | Get inertia matrix (9 elements) |
| `get_volume_topology` | vol_id | `dict` | Get topology summary |
| `is_meshed` | vol_id | `bool` | Check if volume/model is meshed |
| `get_adjacent_volumes` | vol_id | `list` | Get volumes sharing a surface |
| `get_shared_surfaces` | vol_id1, vol_id2 | `list` | Get surfaces shared between volumes |
| `is_point_inside_volume` | vol_tag, point | `bool` | Check if point is inside volume |

### Surface/Curve/Point Info

| Method | Parameters | Returns | Description |
|--------|-----------|---------|-------------|
| `get_surface_area` | surf_tag | `float` | Get surface area |
| `get_curve_length` | curve_tag | `float` | Get curve length |
| `get_surfaces` | vol_id | `list` | Get surfaces bounding a volume |
| `get_curves` | surf_tag | `list` | Get curves bounding a surface |
| `get_vertices_of` | curve_tag | `list` | Get vertices of a curve |
| `get_surface_normal` | surf_tag, uv | `list` | Get surface normal at point |
| `get_surface_type` | surf_tag | `str` | Get surface type name |
| `is_planar` | surf_tag | `bool` | Check if surface is planar |
| `get_point_on_curve` | curve_tag, t | `list` | Get point on curve at parameter t |
| `get_curve_tangent` | curve_tag, t | `list` | Get tangent on curve at t |
| `get_surface_curvature` | surf_tag, u, v | `float` | Get curvature at surface point |
| `get_surface_derivative` | surf_tag, u, v | `list` | Get surface derivatives |
| `get_second_derivative` | dim, tag, params | `list` | Get second derivatives |
| `get_point_coordinates` | point_tag | `list` | Get point [x,y,z] |
| `is_curve_closed` | curve_tag | `bool` | Check if curve is closed |

### Entity Listing

| Method | Parameters | Returns | Description |
|--------|-----------|---------|-------------|
| `get_all_surfaces` | -- | `list` | Get all surface tags |
| `get_all_curves` | -- | `list` | Get all curve tags |
| `get_all_points` | -- | `list` | Get all point tags |
| `get_entities` | dim | `list` | Get all entities of a dimension |
| `list_entities` | dim | `None` | Print entity listing |
| `get_entity_count` | dim | `int` | Count entities of a dimension |
| `get_model_summary` | -- | `dict` | Get model summary |
| `get_named_entities` | -- | `dict` | Get entities with names |
| `find_volumes_by_name` | pattern | `list` | Find volumes matching name pattern |

### Distance and Search

| Method | Parameters | Returns | Description |
|--------|-----------|---------|-------------|
| `get_distance` | vol_id1, vol_id2 | `float` | Distance between volume centers |
| `find_closest_entity` | point, dim | `dict` | Find closest entity to a point |
| `measure_distance_point_to_surface` | point, surf_tag | `dict` | Distance from point to surface |
| `distance_point_to_curve` | point, curve_tag | `dict` | Distance from point to curve |
| `distance_point_to_entity` | point, dim, tag | `dict` | Distance from point to entity |
| `distance_between_volumes` | vol_id1, vol_id2 | `dict` | Distance between two volumes |
| `closest_point_on_curve` | curve_tag, point | `dict` | Closest point on curve |
| `closest_point_on_entity` | dim, tag, point | `dict` | Closest point on entity |

### Naming and Attributes

| Method | Parameters | Returns | Description |
|--------|-----------|---------|-------------|
| `set_name` | vol_id, name | `None` | Set volume name |
| `get_name` | vol_id | `str` | Get volume name |
| `set_attribute` | dim, tag, name, values | `None` | Set entity attribute |
| `get_attribute` | dim, tag, name | `list` | Get entity attribute |
| `get_attribute_names` | dim, tag | `list` | Get entity attribute names |
| `remove_attribute` | dim, tag, name | `None` | Remove entity attribute |
| `set_model_filename` | filename | `None` | Set model filename |
| `get_model_filename` | -- | `str` | Get model filename |

### Aggregates

| Method | Parameters | Returns | Description |
|--------|-----------|---------|-------------|
| `get_bounding_box_all` | -- | `dict` | Get bounding box of all volumes |
| `get_total_volume` | -- | `float` | Get total volume of all volumes |
| `get_total_surface_area` | -- | `float` | Get total surface area |
| `get_total_curve_length` | -- | `float` | Get total curve length |
| `get_volumes_in_bounding_box` | xmin, ymin, zmin, xmax, ymax, zmax | `list` | Get volumes in bounding box |
| `get_surfaces_in_bounding_box` | xmin, ymin, zmin, xmax, ymax, zmax | `list` | Get surfaces in bounding box |
| `get_curves_in_bounding_box` | xmin, ymin, zmin, xmax, ymax, zmax | `list` | Get curves in bounding box |

### Adjacency and Topology

| Method | Parameters | Returns | Description |
|--------|-----------|---------|-------------|
| `get_upward_adjacencies` | dim, tag | `list` | Get higher-dim adjacent entities |
| `get_downward_adjacencies` | dim, tag | `list` | Get lower-dim adjacent entities |
| `get_parent` | dim, tag | `list` | Get parent entities |
| `get_parametric_bounds` | dim, tag | `dict` | Get parametric bounds |
| `point_on_surface` | surf_tag, u, v | `list` | Get point at parametric coords |
| `get_parametrization` | dim, tag, point | `list` | Get parametric coords for point |
| `reparametrize_on_surface` | dim, tag, surf_tag, params | `list` | Reparametrize on surface |
| `get_surface_loops_query` | vol_tag | `list` | Get surface loops of volume |
| `get_curve_loops_query` | surf_tag | `list` | Get curve loops of surface |
| `get_surface_count` | vol_id | `int` | Count surfaces of volume |
| `get_curve_count` | surf_tag | `int` | Count curves of surface |
| `get_entity_partitions` | dim, tag | `list` | Get entity partition info |
| `get_number_of_partitions` | -- | `int` | Get total partition count |
| `get_changed_entities` | dim | `list` | Get entities changed by last geometry operation |

---

## 11. ExportMixin (export.py)

Mesh and geometry export. Cubit equivalents: export mesh, export step, export stl.

### Mesh Export

| Method | Parameters | Returns | Description |
|--------|-----------|---------|-------------|
| `export` | filename, version (default 2.2) | `None` | Export to GMSH .msh (Radia/NGSolve compatible) |
| `export_msh` | filename, version (default 4.1) | `None` | Export to MSH with specified version |
| `export_mesh2d` | filename, version | `None` | Export 2D surface mesh only |
| `export_stl` | filename | `None` | Export to STL |
| `export_vtk` | filename | `None` | Export to VTK |
| `export_vtu` | filename | `None` | Export to VTU |
| `export_ply` | filename | `None` | Export to PLY |
| `export_pos` | filename | `None` | Export to Gmsh POS |
| `export_unv` | filename | `None` | Export to IDEAS Universal (.unv) |
| `export_abaqus` | filename | `None` | Export to Abaqus INP |
| `export_cgns` | filename | `None` | Export to CGNS |
| `export_med` | filename | `None` | Export to MED |
| `export_auto` | filename | `None` | Auto-detect format from extension |

### Geometry Export

| Method | Parameters | Returns | Description |
|--------|-----------|---------|-------------|
| `export_step` | filename | `None` | Export to STEP |
| `export_brep` | filename | `None` | Export to BREP |
| `export_geo` | filename | `None` | Export to GEO (Gmsh script) |
| `save_session` | filename | `None` | Save as BREP |
| `show` | dim | `None` | Show model in GMSH GUI |

### Import / Merge

| Method | Parameters | Returns | Description |
|--------|-----------|---------|-------------|
| `merge_file` | filename: str | `None` | Merge a file (mesh or geometry) into the current model |

### Radia Integration

| Method | Parameters | Returns | Description |
|--------|-----------|---------|-------------|
| `to_radia` | mu_r, magnetization | `int` | Convert to Radia container object |
| `to_radia_per_volume` | materials: dict | `int` | Convert with per-volume materials |
| `to_radia_with_material` | vol_id, mu_r, magnetization | `int` | Export single volume to Radia |

### NGSolve Integration

| Method | Parameters | Returns | Description |
|--------|-----------|---------|-------------|
| `to_ngsolve_surface` | label, mesh_size | `ngsolve.Mesh` | Convert to NGSolve surface mesh |
| `to_ngsolve_volume` | -- | `ngsolve.Mesh` | Convert to NGSolve volume mesh |

### Data Conversion

| Method | Parameters | Returns | Description |
|--------|-----------|---------|-------------|
| `to_numpy` | -- | `dict` | Export mesh as numpy arrays |
| `to_pyvista` | -- | `pyvista.UnstructuredGrid` | Export as PyVista grid |
| `to_dict` | -- | `dict` | Export mesh as plain Python dict |
| `to_meshio` | -- | `meshio.Mesh` | Export via meshio |
| `to_scipy_sparse` | -- | `csr_matrix` | Export connectivity as sparse matrix |
| `get_export_formats` | -- | `list[str]` | List supported export formats |

---

## 12. EntityWrapperMixin (entity_wrappers.py)

Gmsh-based entity wrapper queries operating on raw GMSH entity tags (not GmshBuilder volume IDs).

### Volume Queries

| Method | Parameters | Returns | Description |
|--------|-----------|---------|-------------|
| `volume_surfaces` | vol_tag | `list[int]` | Get surfaces bounding a volume |
| `volume_curves` | vol_tag | `list[int]` | Get curves bounding a volume |
| `volume_vertices` | vol_tag | `list[int]` | Get vertices of a volume |
| `volume_bounding_box` | vol_tag | `dict` | Get bounding box of a volume |
| `volume_center` | vol_tag | `list` | Get center of mass of a volume |
| `volume_measure` | vol_tag | `float` | Get volume (3D measure) |
| `volume_inertia` | vol_tag | `list` | Get inertia matrix (9 elements) |
| `volume_type` | vol_tag | `str` | Get volume type name |
| `volume_adjacent_volumes` | vol_tag | `list[int]` | Get adjacent volumes |
| `volume_physical_groups` | vol_tag | `list[int]` | Get physical groups for volume |
| `volume_is_inside` | vol_tag, point | `bool` | Check if point is inside volume |
| `volume_surface_area` | vol_tag | `float` | Get total surface area of volume |
| `volume_curve_count` | vol_tag | `int` | Count curves bounding volume |
| `volume_vertex_count` | vol_tag | `int` | Count vertices of volume |

### Surface Queries

| Method | Parameters | Returns | Description |
|--------|-----------|---------|-------------|
| `surface_curves` | surf_tag | `list[int]` | Get curves bounding a surface |
| `surface_vertices` | surf_tag | `list[int]` | Get vertices of a surface |
| `surface_volumes` | surf_tag | `list[int]` | Get volumes adjacent to surface |
| `surface_bounding_box` | surf_tag | `dict` | Get bounding box of a surface |
| `surface_center` | surf_tag | `list` | Get center of mass of a surface |
| `surface_area` | surf_tag | `float` | Get area of a surface |
| `surface_type` | surf_tag | `str` | Get surface type name |
| `surface_normal_at` | surf_tag, u, v | `list` | Get outward normal at parametric point |
| `surface_closest_point` | surf_tag, point | `dict` | Find closest point on surface |
| `surface_is_planar` | surf_tag | `bool` | Check if surface is planar |
| `surface_is_cylindrical` | surf_tag | `bool` | Check if surface is cylindrical |
| `surface_is_conical` | surf_tag | `bool` | Check if surface is conical |
| `surface_is_spherical` | surf_tag | `bool` | Check if surface is spherical |
| `surface_physical_groups` | surf_tag | `list[int]` | Get physical groups for surface |

### Curve Queries

| Method | Parameters | Returns | Description |
|--------|-----------|---------|-------------|
| `curve_vertices` | curve_tag | `list[int]` | Get endpoint vertices of a curve |
| `curve_surfaces` | curve_tag | `list[int]` | Get surfaces adjacent to curve |
| `curve_bounding_box` | curve_tag | `dict` | Get bounding box of a curve |
| `curve_center` | curve_tag | `list` | Get center of mass of a curve |
| `curve_length` | curve_tag | `float` | Get length of a curve |
| `curve_type` | curve_tag | `str` | Get curve type name |
| `curve_point_at` | curve_tag, t | `list` | Get point on curve at parameter |
| `curve_tangent_at` | curve_tag, t | `list` | Get tangent at parameter |
| `curve_curvature_at` | curve_tag, t | `float` | Get curvature at parameter |
| `curve_closest_point` | curve_tag, point | `dict` | Find closest point on curve |
| `curve_is_line` | curve_tag | `bool` | Check if curve is a line |
| `curve_is_circle` | curve_tag | `bool` | Check if curve is circular |
| `curve_is_bspline` | curve_tag | `bool` | Check if curve is B-spline |
| `curve_physical_groups` | curve_tag | `list[int]` | Get physical groups for curve |

### Vertex Queries

| Method | Parameters | Returns | Description |
|--------|-----------|---------|-------------|
| `vertex_coordinates` | point_tag | `list` | Get [x,y,z] coordinates |
| `vertex_curves` | point_tag | `list[int]` | Get curves adjacent to vertex |
| `vertex_surfaces` | point_tag | `list[int]` | Get surfaces reachable from vertex |
| `vertex_volumes` | point_tag | `list[int]` | Get volumes reachable from vertex |
| `vertex_physical_groups` | point_tag | `list[int]` | Get physical groups for vertex |
| `vertex_distance_to` | point_tag, point | `float` | Distance from vertex to point |
| `vertex_on_curve` | point_tag, curve_tag | `bool` | Check if vertex is on curve |
| `vertex_on_surface` | point_tag, surf_tag | `bool` | Check if vertex is on surface |
| `vertex_on_volume` | point_tag, vol_tag | `bool` | Check if vertex is on volume boundary |
| `vertex_bounding_box` | point_tag | `dict` | Get bounding box (degenerate) |
| `get_vertex_count` | -- | `int` | Count all point entities |
| `get_all_vertices` | -- | `list[int]` | Get all point tags |
| `vertices_in_bounding_box` | xmin, ymin, zmin, xmax, ymax, zmax | `list[int]` | Get vertices in bounding box |

---

## 13. GeometryAnalysisMixin (geometry_analysis.py)

Geometry cleanup, distance/overlap detection, topology analysis, and differential geometry. Cubit equivalents: heal, defeature, validate, analysis.

### Cleanup

| Method | Parameters | Returns | Description |
|--------|-----------|---------|-------------|
| `heal_all` | tolerance | `dict` | Heal all shapes in model |
| `defeature_volume` | vol_tag, tolerance | `list[int]` | Remove small features from volume |
| `remove_small_faces` | vol_tag, tolerance | `int` | Remove surfaces below area threshold |
| `find_small_curves` | min_length | `list[tuple]` | Find curves shorter than minimum |
| `find_small_surfaces` | min_area | `list[tuple]` | Find surfaces below area minimum |
| `find_small_volumes` | min_volume | `list[tuple]` | Find volumes below volume minimum |
| `find_narrow_surfaces` | aspect_ratio_max | `list[tuple]` | Find sliver surfaces |
| `remove_duplicates` | -- | `int` | Remove duplicate entities |
| `find_orphan_entities` | dim | `list[int]` | Find unreferenced entities |
| `clean_geometry` | vol_ids, tolerance | `dict` | Combined geometry cleanup |

### Distance and Overlap Detection

| Method | Parameters | Returns | Description |
|--------|-----------|---------|-------------|
| `distance_between_entities` | dim1, tag1, dim2, tag2 | `float` | Min distance between entities |
| `closest_entities` | dim1, tag1, dim2, tag2 | `dict` | Closest points between entities |
| `find_overlapping_volumes` | tolerance | `list[tuple]` | Find volumes with overlapping bboxes |
| `find_touching_surfaces` | vol_tag1, vol_tag2, tolerance | `list[tuple]` | Find close surfaces between volumes |
| `entities_in_bounding_box` | xmin, ymin, zmin, xmax, ymax, zmax, dim | `list[int]` | Find entities in bounding box |
| `find_coincident_vertices` | tolerance | `list[list]` | Find overlapping vertices |
| `find_coincident_curves` | tolerance | `list[tuple]` | Find coincident curves |
| `check_watertight` | vol_tag | `bool` | Check if volume is watertight |

### Topology Analysis

| Method | Parameters | Returns | Description |
|--------|-----------|---------|-------------|
| `get_curve_loops` | surf_tag | `list[list]` | Get curve loops of a surface |
| `get_surface_loops` | vol_tag | `list[list]` | Get surface loops of a volume |
| `get_adjacencies` | dim, tag | `dict` | Get upward and downward adjacencies |
| `get_entity_type` | dim, tag | `str` | Get entity type name |
| `get_entity_name_raw` | dim, tag | `str` | Get raw entity name |
| `set_entity_name_raw` | dim, tag, name | `None` | Set raw entity name |
| `is_entity_orphan` | dim, tag | `bool` | Check if entity is orphaned |
| `count_topology` | vol_tag | `dict` | Count surfaces/curves/vertices |
| `find_nonmanifold_edges` | -- | `list[int]` | Find non-manifold edges |
| `find_free_surfaces` | -- | `list[int]` | Find free (one-sided) surfaces |

### Differential Geometry

| Method | Parameters | Returns | Description |
|--------|-----------|---------|-------------|
| `curvature_at_curve` | curve_tag, t | `float` | Get curvature at curve point |
| `curvature_at_surface` | surf_tag, u, v | `float` | Get mean curvature at surface point |
| `principal_curvatures` | surf_tag, u, v | `dict` | Get principal curvatures and directions |
| `derivative_curve` | curve_tag, t | `list` | Get first derivative of curve |
| `derivative_surface` | surf_tag, u, v | `list` | Get first derivatives of surface |
| `second_derivative_curve` | curve_tag, t | `list` | Get second derivative of curve |
| `second_derivative_surface` | surf_tag, u, v | `list` | Get second derivatives of surface |

---

## 14. GroupMixin (groups.py)

Selection group management and set operations. Cubit equivalents: group create/delete/add/remove, set operations.

### Group Lifecycle

| Method | Parameters | Returns | Description |
|--------|-----------|---------|-------------|
| `create_group` | name: str | `int` | Create a new selection group |
| `delete_group` | group_id | `None` | Delete a group |
| `rename_group` | group_id, new_name | `None` | Rename a group |
| `get_group_name` | group_id | `str` | Get group name |
| `get_group_count` | -- | `int` | Count groups |
| `get_all_groups` | -- | `list[dict]` | Get summary for all groups |
| `find_group_by_name` | name | `int or None` | Find group by name |
| `clear_all_groups` | -- | `None` | Remove all groups |
| `group_exists` | group_id | `bool` | Check if group exists |
| `get_group_ids` | -- | `list[int]` | Get sorted group IDs |

### Entity Management

| Method | Parameters | Returns | Description |
|--------|-----------|---------|-------------|
| `group_add_volumes` | group_id, vol_ids | `None` | Add volumes to group |
| `group_add_surfaces` | group_id, surf_tags | `None` | Add surfaces to group |
| `group_add_curves` | group_id, curve_tags | `None` | Add curves to group |
| `group_add_vertices` | group_id, point_tags | `None` | Add vertices to group |
| `group_remove_volumes` | group_id, vol_ids | `None` | Remove volumes from group |
| `group_remove_surfaces` | group_id, surf_tags | `None` | Remove surfaces from group |
| `group_remove_curves` | group_id, curve_tags | `None` | Remove curves from group |
| `group_remove_vertices` | group_id, point_tags | `None` | Remove vertices from group |
| `group_get_entities` | group_id, dim | `list` | Get entities in group |
| `group_entity_count` | group_id, dim | `int` | Count entities in group |

### Curvesets and Set Operations

| Method | Parameters | Returns | Description |
|--------|-----------|---------|-------------|
| `add_curveset` | curve_tags, name | `int` | Create 1D physical group |
| `add_curveset_by_name` | curve_tags, name | `int` | Create named curveset |
| `get_all_curvesets` | -- | `list[dict]` | Get all curveset info |
| `get_curveset_count` | -- | `int` | Count curvesets |
| `get_curveset_curves` | pg_tag | `list[int]` | Get curves in curveset |
| `get_pg_for_entity` | dim, tag | `list[int]` | Get physical groups for entity |
| `get_entities_by_name` | name | `list[tuple]` | Find entities by physical group name |
| `merge_groups` | group_id1, group_id2 | `int` | Union of two groups |
| `subtract_groups` | group_id1, group_id2 | `int` | Difference of two groups |
| `intersect_groups` | group_id1, group_id2 | `int` | Intersection of two groups |

---

## 15. PostProcessingMixin (post_processing.py)

Post-processing view management. Wraps GMSH's gmsh.view API for creating, populating, querying, and configuring views.

### View Lifecycle

| Method | Parameters | Returns | Description |
|--------|-----------|---------|-------------|
| `add_view` | name: str | `int` | Create a new view |
| `remove_view` | view_tag | `None` | Remove a view |
| `get_all_views` | -- | `list[int]` | Get all view tags |
| `get_view_count` | -- | `int` | Count views |
| `clear_all_views` | -- | `None` | Remove all views |
| `write_view` | view_tag, filename | `None` | Write view to file |
| `combine_views` | what, how, remove_original | `None` | Combine existing views |
| `add_alias_view` | view_tag, copy_options | `int` | Create alias of view |

### Add Data

| Method | Parameters | Returns | Description |
|--------|-----------|---------|-------------|
| `add_scalar_node_data` | view_tag, step, time, node_tags, values | `None` | Add scalar data at nodes |
| `add_vector_node_data` | view_tag, step, time, node_tags, values | `None` | Add vector data at nodes |
| `add_tensor_node_data` | view_tag, step, time, node_tags, values | `None` | Add tensor data at nodes |
| `add_scalar_element_data` | view_tag, step, time, elem_tags, values | `None` | Add scalar data on elements |
| `add_vector_element_data` | view_tag, step, time, elem_tags, values | `None` | Add vector data on elements |
| `add_list_data` | view_tag, data_type, n_elem, data | `None` | Add list-based (non-model) data |
| `add_list_data_string` | view_tag, coord, data, style | `None` | Add string annotation |
| `add_homogeneous_data` | view_tag, step, time, data_type, n_components, data | `None` | Add homogeneous model data |
| `add_field_to_view` | view_tag, field_callable, step, time | `None` | Evaluate callable at all nodes |

### Get / Query Data

| Method | Parameters | Returns | Description |
|--------|-----------|---------|-------------|
| `get_view_model_data` | view_tag, step | `dict` | Retrieve model data from view |
| `get_view_list_data` | view_tag | `dict` | Retrieve list data from view |
| `probe_view` | view_tag, x, y, z, step | `list` | Probe view at spatial location |
| `get_view_data_range` | view_tag | `tuple` | Get (min, max) value range |

### View Options

| Method | Parameters | Returns | Description |
|--------|-----------|---------|-------------|
| `set_view_option_number` | view_tag, name, value | `None` | Set numeric view option |
| `set_view_option_string` | view_tag, name, value | `None` | Set string view option |
| `get_view_option_number` | view_tag, name | `float` | Get numeric view option |
| `get_view_option_string` | view_tag, name | `str` | Get string view option |

---

## 16. OptionsMixin (options.py)

GMSH global options. Cubit equivalents: set/get options for geometry, mesh, colors, visibility.

### Generic Options

| Method | Parameters | Returns | Description |
|--------|-----------|---------|-------------|
| `set_option_number` | name, value | `None` | Set a GMSH numeric option |
| `get_option_number` | name | `float` | Get a GMSH numeric option |
| `set_option_string` | name, value | `None` | Set a GMSH string option |
| `get_option_string` | name | `str` | Get a GMSH string option |
| `restore_default_options` | -- | `None` | Restore all defaults |
| `get_all_mesh_options` | -- | `dict` | Read common Mesh.* options |

### Geometry Options

| Method | Parameters | Returns | Description |
|--------|-----------|---------|-------------|
| `set_geometry_tolerance` | tolerance: float | `None` | Set OCC geometry tolerance |
| `get_geometry_tolerance` | -- | `float` | Get geometry tolerance |
| `set_occ_scaling` | factor: float | `None` | Set OCC scaling factor |
| `set_occ_parallel` | enabled: bool | `None` | Enable/disable parallel OCC |

### Mesh Behavior Options

| Method | Parameters | Returns | Description |
|--------|-----------|---------|-------------|
| `set_mesh_only_visible` | enabled: bool | `None` | Mesh only visible entities |
| `set_mesh_save_all` | enabled: bool | `None` | Save all elements on export |
| `set_mesh_save_groups_of_elements` | value: int | `None` | Control element group saving |
| `set_mesh_save_groups_of_nodes` | value: int | `None` | Control node group saving |
| `set_mesh_binary` | enabled: bool | `None` | Enable binary mesh output |
| `set_mesh_random_factor` | factor: float | `None` | Set Delaunay random factor |

### Color and Visibility

| Method | Parameters | Returns | Description |
|--------|-----------|---------|-------------|
| `set_entity_color` | dim, tag, r, g, b, a | `None` | Set entity color (RGBA 0-255) |
| `get_entity_color` | dim, tag | `tuple` | Get entity color (r,g,b,a) |
| `set_entity_visibility` | dim, tag, visible | `None` | Set entity visibility |
| `get_entity_visibility` | dim, tag | `bool` | Get entity visibility |

---

## Appendix: CFD Meshing Techniques

GmshBuilder uses OCC kernel only. `extrudeBoundaryLayer` (geo kernel) is not available;
use field-based boundary layers and structured extrusion instead.

### Technique 1: Wall-Normal Grading

```python
# Single-call Distance+Threshold for boundary layer resolution
f = gb.add_size_gradient(wall_surfs, dim=2,
                         lc_near=0.0002,  # first cell (y+ ~ 1)
                         lc_far=0.005,    # far field
                         dist_near=0.001, dist_far=0.02)
```

Y+ estimation: `y_wall(y+=1) ~ 6 * nu / u_tau`, where `u_tau = sqrt(tau_w/rho)`.

### Technique 2: Structured BL via Extrusion

```python
# Graded layers: 5 thin (20% height) + 5 thick (80% height)
gb.extrude_surface(rect, 0, 0, 0.05,
                   num_elements=[5, 5], heights=[0.2, 1.0], recombine=True)
```

### Technique 3: Multi-Zone Refinement (Field Composition)

```python
f1 = gb.add_size_gradient(body_surfs, 2, 0.002, 0.05, 0.005, 0.1, set_as_background=False)
f2 = gb.add_wake_refinement(trailing_curves, [1,0,0], 1.0, 0.005, 0.05)
f3 = gb.add_field_box(xmin, xmax, ymin, ymax, zmin, zmax, 0.01, 0.05)
gb.set_background_mesh_from_fields([f1, f2, f3], operator='min')
```

### Technique 4: Airfoil 2D Mesh

Create NACA profile from spline points, apply wall grading + wake refinement.

### Technique 5: Pipe Flow (Structured Hex)

Revolve 2D rectangle with `num_elements` for structured axisymmetric mesh.

### Technique 6: STL Terrain

```python
gb.import_stl_mesh('terrain.stl', angle=30)
gb.classify_surfaces_parametric(angle=25, curve_angle=120)
gb.create_geometry_from_mesh()
f = gb.add_terrain_field('terrain_size.pos', format='structured')
```

### Technique 7: Adaptive Refinement

```python
f = gb.add_field_math("0.005 + 0.015 * sqrt(x*x + y*y)")
gb.adapt_mesh_iterative(f, n_iterations=3, size_factor=0.7)
```

### Technique 8: BoundaryLayer Field Composition

```python
f_bl = gb.add_boundary_layer_field(obs_surfs, 2, 0.0005, 0.01, 0.001, 0.02,
                                    ratio=1.3, n_layers=8)
f_far = gb.add_field_box(...)
gb.set_background_mesh_from_fields([f_bl, f_far], operator='min')
```

### OCC Kernel Limitations

| Feature | Workaround |
|---------|------------|
| `extrudeBoundaryLayer` (geo only) | Field-based BL + structured extrusion |
| Prism layers | `extrude_surface` + `num_elements` + graded `heights` |
| Anisotropic stretching | `add_wake_refinement` Box field |
