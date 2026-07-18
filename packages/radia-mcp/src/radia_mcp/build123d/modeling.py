# -*- coding: utf-8 -*-
r"""Parametric solid-modeling operations for CAE/EM geometry in build123d.

A small, **tested** library of the parametric modeling operations a solid-modeller user reaches for
when authoring electromagnetic structures -- the building blocks that build123d has as raw kernel
calls but that are fiddly to get clean and meshable by hand: a robust **annular wedge/sector**, a
hollow **tube**, a **racetrack** (rounded-rectangle) coil, and the **array-with-copies** transforms
(polar / linear / mirror) plus a labelled multi-region **assembly**.  Every helper returns a clean,
``is_valid`` build123d solid (or a labelled :class:`~build123d.Compound`), CAE-safe (primitives +
boolean, no micro-edges) and ready for the Netgen -> Radia / NGSolve tet pipeline.

Design rules (see build123d_usage "cae_guidelines"):
  * primitives + boolean only -- no tiny fillets, no over-constrained sketches;
  * every solid carries a ``.label`` (maps to a material region / Gmsh physical group);
  * z-centred solids (span -h/2..h/2) so stacking and symmetry planes are simple;
  * pure build123d -- no radia_mcp dependency, so the helpers are equally usable inside an
    ``execute_build123d`` subprocess or copy-pasted into a user script.

These are intentionally generic CAD operations (sector, array, tube, racetrack); the EM *device*
archetypes that compose them (Halbach ring, dipole/quadrupole yokes, C-core, solenoid ...) live in
:mod:`radia_mcp.build123d.archetypes`.
"""
from __future__ import annotations

from collections import Counter
import math

try:
    from build123d import (Axis, Box, BuildLine, BuildSketch, CenterArc, Circle, Compound, Cylinder,
                           Helix, Keep, Line, Mode, Plane, Pos, Rectangle, RectangleRounded, Spline,
                           chamfer, extrude, fillet, loft, make_face, offset, revolve, split, sweep)
except ImportError as exc:
    _BUILD123D_IMPORT_ERROR = exc

    def _missing_build123d(*_args, **_kwargs):
        raise ImportError(
            "build123d geometry helpers require the optional 'build123d' dependency"
        ) from _BUILD123D_IMPORT_ERROR

    class _MissingBuild123dType:
        def __new__(cls, *_args, **_kwargs):
            _missing_build123d()

    class _MissingBuild123dNamespace:
        def __getattr__(self, _name):
            _missing_build123d()

    Axis = Keep = Mode = Plane = _MissingBuild123dNamespace()
    Compound = _MissingBuild123dType
    Box = BuildLine = BuildSketch = CenterArc = Circle = Cylinder = Helix = _missing_build123d
    Line = Pos = Rectangle = RectangleRounded = Spline = _missing_build123d
    chamfer = extrude = fillet = loft = make_face = offset = revolve = split = sweep = _missing_build123d
else:
    _BUILD123D_IMPORT_ERROR = None


def _require_build123d():
    if _BUILD123D_IMPORT_ERROR is not None:
        _missing_build123d()

__all__ = ["annular_segment", "tube", "racetrack_coil", "polar_array", "linear_array",
           "mirrored", "assembly", "shape_envelope_row", "enclosing_box",
           "enclosure_clearance_row", "enclosure_difference_region",
           "shape_measurement_row", "shape_measurement_rows",
           "box_through_cylinder_reference_row", "mounting_plate_boss_reference_row",
           "keyed_terminal_plate_reference_row", "flanged_sleeve_reference_row",
           "coax_annular_sleeve_reference_row",
           "ribbed_busbar_heat_sink_reference_row",
           "three_phase_busbar_snubber_plate_reference_row",
           "v_type_ipm_rotor_coupon_reference_row",
           "rcd_snubber_heat_spreader_reference_row",
           "rcd_snubber_capacitance_sweep_rows",
           "thermal_robin_cooling_plate_reference_row",
           "motor_housing_radial_fin_reference_row",
           "box_face_vector_area_rows", "box_face_pressure_force_rows",
           "box_face_pressure_moment_rows", "box_face_pressure_resultant_summary",
           "box_face_traction_moment_rows",
           "compare_boundary_vector_area_rows",
           "compare_shape_measurement_rows", "shape_measurement_comparison_summary",
           "compare_shape_volume_rows", "shape_volume_crosscheck_summary",
           "shape_perforated_prism_roundtrip_gate",
           "shape_volume_crosscheck_source_coverage_gate",
           "shape_volume_crosscheck_source_identity_gate",
           "shape_external_cad_volume_evidence_package_gate",
           "shape_cad_route_source_contract_gate",
           "cst_cad_volume_export_manifest_gate",
           "shape_name_identity_gate",
           "shape_role_metadata_gate",
           "shape_transition_role_metadata_gate",
           "shape_cubit_meshing_scheme_intent_gate",
           "shape_mass_property_crosscheck_summary",
           "shape_cubit_export_package_handoff_gate",
           "shape_cubit_quality_package_handoff_gate",
           "shape_cubit_quality_ledger_handoff_gate",
           "shape_cubit_solver_route_handoff_gate",
           "shape_cad_handoff_manifest_gate",
           "shape_submodel_cad_handoff_gate",
           "shape_curvilinear_mesh_intent_gate",
           "shape_mesh_environment_handoff_gate",
           "shape_measurement_inventory_summary", "worst_shape_measurement_comparison_rows",
           "shape_measurement_health_summary", "shape_bbox_pair_clearance_summary",
           "shape_parameter_sweep_summary",
           # generic solid-modelling operations (constructors / local mods / arrays)
           "swept", "revolved", "lofted", "coil", "helix_centerline_length",
           "round_wire_helix_metrics", "strut", "thicken", "draft_extrude",
           "shell", "fillet_edges", "chamfer_edges", "grid_array", "path_array",
           # boolean / slice / sheet-metal
           "fuse", "cut", "common", "slice_solid", "bend_sheet",
           "fillet_varied", "chamfer_varied", "slice_array"]


def annular_segment(r_in, r_out, h, start_angle=0.0, end_angle=90.0, label="segment"):
    r"""A radial **wedge of an annulus**: ``r_in < r < r_out``, angular span
    ``[start_angle, end_angle]`` degrees, height ``h`` (z-centred).

    Robust for any span ``< 360`` deg: the full annulus is INTERSECTED with a pie wedge swept from the
    centre, so there are no self-intersecting sketch boundaries (the failure mode of a hand-built
    sector polyline).  The canonical building block for segmented PM arrays (Halbach), pole sectors,
    and angular slices.  Volume ``= (span/360) * pi (r_out^2 - r_in^2) * h``.
    """
    if not (0.0 < r_in < r_out):
        raise ValueError("require 0 < r_in < r_out")
    span = end_angle - start_angle
    if not (0.0 < span < 360.0):
        raise ValueError("require 0 < (end_angle - start_angle) < 360")
    annulus = Cylinder(radius=r_out, height=h) - Cylinder(radius=r_in, height=h)
    big = r_out * 1.5
    a0, a1 = math.radians(start_angle), math.radians(end_angle)
    with BuildSketch(Plane.XY) as sk:
        with BuildLine():
            Line((0, 0), (big * math.cos(a0), big * math.sin(a0)))
            CenterArc((0, 0), big, start_angle, span)
            Line((big * math.cos(a1), big * math.sin(a1)), (0, 0))
        make_face()
    wedge = Pos(0, 0, -h / 2) * extrude(sk.sketch, amount=h)
    seg = (annulus & wedge).solid()        # single labelled Solid (label survives .solids() downstream)
    seg.label = label
    return seg


def tube(r_in, r_out, h, label="tube"):
    r"""Hollow cylinder (annulus extruded), z-centred.  Volume ``= pi (r_out^2 - r_in^2) h``."""
    if not (0.0 <= r_in < r_out):
        raise ValueError("require 0 <= r_in < r_out")
    t = (Cylinder(radius=r_out, height=h) - Cylinder(radius=r_in, height=h)).solid()
    t.label = label
    return t


def racetrack_coil(length, width, band, h, corner_radius, label="coil"):
    r"""A **racetrack** (rounded-rectangle) coil: an outer rounded rectangle minus an inner one,
    extruded by ``h`` (z-centred).  ``length`` x ``width`` is the outer envelope, ``band`` the
    conductor wall thickness, ``corner_radius`` the outer corner radius; the conductor cross-section is
    ``band`` (in-plane) x ``h`` (axial).  The standard flat racetrack / pancake coil outline used for
    accelerator and MRI windings.
    """
    if band <= 0 or band * 2 >= min(length, width):
        raise ValueError("require 0 < 2*band < min(length, width)")
    with BuildSketch(Plane.XY) as sk:
        RectangleRounded(length, width, corner_radius)
        RectangleRounded(length - 2 * band, width - 2 * band,
                         max(corner_radius - band, 0.1), mode=Mode.SUBTRACT)
    coil = (Pos(0, 0, -h / 2) * extrude(sk.sketch, amount=h)).solid()
    coil.label = label
    return coil


def polar_array(part, count, total_angle=360.0, axis=None, label=None, label_fmt="{base}_{k:02d}"):
    r"""``count`` rotated copies of ``part`` about ``axis``, returned as a labelled
    :class:`~build123d.Compound`.  A full ``360`` deg ring spaces the copies by ``360/count``; a
    partial fan (``total_angle < 360``) spaces them by ``total_angle/(count-1)`` so the first and last
    copies sit at the fan ends.  The "rotate with copies" / segmented-array verb.
    """
    if axis is None:
        _require_build123d()
        axis = Axis.Z
    if count < 1:
        raise ValueError("count must be >= 1")
    base = label or (part.label or "part")
    full = abs(total_angle) >= 360.0 - 1e-9
    step = total_angle / count if full else (total_angle / (count - 1) if count > 1 else 0.0)
    children = []
    for k in range(count):
        c = part.rotate(axis, k * step)
        c.label = label_fmt.format(base=base, k=k)
        children.append(c)
    return Compound(children=children, label=base + "_array")


def linear_array(part, count, spacing, direction=(1, 0, 0), label=None, label_fmt="{base}_{k:02d}"):
    r"""``count`` translated copies of ``part`` along (normalised) ``direction`` at ``spacing``,
    returned as a labelled :class:`~build123d.Compound`.  The "translate with copies" verb.
    """
    if count < 1:
        raise ValueError("count must be >= 1")
    base = label or (part.label or "part")
    d = direction
    norm = math.sqrt(d[0] ** 2 + d[1] ** 2 + d[2] ** 2) or 1.0
    ux, uy, uz = d[0] / norm, d[1] / norm, d[2] / norm
    children = []
    for k in range(count):
        c = Pos(k * spacing * ux, k * spacing * uy, k * spacing * uz) * part
        c.label = label_fmt.format(base=base, k=k)
        children.append(c)
    return Compound(children=children, label=base + "_array")


def mirrored(part, about=None, keep_original=True, label=None):
    r"""The mirror image of ``part`` across plane ``about`` (default ``Plane.XZ``); with
    ``keep_original`` (default) returns a :class:`~build123d.Compound` of the original + its mirror --
    the symmetry-completion verb (build a quarter/half model, then mirror to whole).
    """
    _require_build123d()
    if about is None:
        about = Plane.XZ
    from build123d import mirror as _mirror
    base = label or (part.label or "part")
    mir = _mirror(part, about=about)
    mir.label = base + "_mir"
    if not keep_original:
        return mir
    orig = part
    try:
        orig.label = base + "_orig"
    except Exception:
        pass
    return Compound(children=[orig, mir], label=base + "_mirrored")


def assembly(*parts, label="assembly"):
    r"""Group ``parts`` (labelled solids and/or Compounds) into a single multi-region
    :class:`~build123d.Compound`, KEEPING each part as a separate **labelled child** (do NOT fuse) so
    it maps to its own material region / Gmsh physical group downstream.  Flattens one level of nested
    Compounds via ``.children`` (NOT ``.solids()`` -- the latter drops the region labels).
    """
    children = []
    for p in parts:
        if isinstance(p, Compound) and p.children:
            children.extend(p.children)          # labelled sub-parts: keep them (solids() strips labels)
        elif isinstance(p, Compound):
            sols = p.solids()                    # raw Part (Box / extrude / ... has no children): its solids
            if p.label and len(sols) == 1:
                sols[0].label = p.label          # carry the Part's own label onto its (stripped) solid
            children.extend(sols)
        else:
            children.append(p)                   # a bare Solid
    return Compound(children=children, label=label)


def _shape_list(shapes):
    if isinstance(shapes, Compound):
        return list(shapes.children) if shapes.children else list(shapes.solids())
    try:
        return list(shapes)
    except TypeError:
        return [shapes]


def _margin_xyz(margin):
    try:
        values = [float(v) for v in margin]
    except TypeError:
        values = [float(margin)] * 3
    if len(values) != 3:
        raise ValueError("margin must be a scalar or a 3-value iterable")
    if any(v < 0.0 for v in values):
        raise ValueError("margin values must be >= 0")
    return values


def shape_envelope_row(shapes, margin=0.0, name="envelope"):
    """Return a JSON-friendly bounding-box envelope for one or more shapes.

    This is a pre-mesh guard for multi-region CAE models: compute the union
    bounding box of labelled inner regions, expand it by a scalar or xyz
    ``margin``, and expose centre/size/volume in the same simple vocabulary as
    :func:`shape_measurement_row`.
    """

    items = _shape_list(shapes)
    if not items:
        raise ValueError("shape_envelope_row needs at least one shape")
    boxes = [shape.bounding_box() for shape in items]
    margin_xyz = _margin_xyz(margin)
    mins = [
        min(bb.min.X for bb in boxes) - margin_xyz[0],
        min(bb.min.Y for bb in boxes) - margin_xyz[1],
        min(bb.min.Z for bb in boxes) - margin_xyz[2],
    ]
    maxs = [
        max(bb.max.X for bb in boxes) + margin_xyz[0],
        max(bb.max.Y for bb in boxes) + margin_xyz[1],
        max(bb.max.Z for bb in boxes) + margin_xyz[2],
    ]
    size = [maxs[i] - mins[i] for i in range(3)]
    center = [(mins[i] + maxs[i]) / 2.0 for i in range(3)]
    return {
        "name": str(name),
        "n_shapes": len(items),
        "margin": margin_xyz,
        "min": [float(v) for v in mins],
        "max": [float(v) for v in maxs],
        "size": [float(v) for v in size],
        "center": [float(v) for v in center],
        "volume": float(size[0] * size[1] * size[2]),
    }


def enclosing_box(shapes, margin=0.0, label="enclosure"):
    """Build a labelled box enclosing one or more shapes with the given margin."""

    row = shape_envelope_row(shapes, margin=margin, name=label)
    sx, sy, sz = row["size"]
    cx, cy, cz = row["center"]
    box = (Pos(cx, cy, cz) * Box(sx, sy, sz)).solid()
    box.label = label
    return box


def enclosure_clearance_row(enclosure, inner_shapes, name=None):
    """Summarise bbox clearances and nominal void volume for an enclosure.

    ``nominal_void_volume`` assumes the inner shapes are disjoint.  It is a
    readable first check before a STEP/Cubit/Netgen round trip, not a
    replacement for boolean validation.
    """

    inner = _shape_list(inner_shapes)
    if not inner:
        raise ValueError("enclosure_clearance_row needs at least one inner shape")
    enc_bb = enclosure.bounding_box()
    inner_row = shape_envelope_row(inner, margin=0.0, name="inner")
    clearances = {
        "xmin": inner_row["min"][0] - float(enc_bb.min.X),
        "xmax": float(enc_bb.max.X) - inner_row["max"][0],
        "ymin": inner_row["min"][1] - float(enc_bb.min.Y),
        "ymax": float(enc_bb.max.Y) - inner_row["max"][1],
        "zmin": inner_row["min"][2] - float(enc_bb.min.Z),
        "zmax": float(enc_bb.max.Z) - inner_row["max"][2],
    }
    inner_volume = sum(float(shape.volume) for shape in inner)
    enclosure_volume = float(enclosure.volume)
    min_clearance = min(clearances.values())
    return {
        "name": str(name or getattr(enclosure, "label", "") or "enclosure"),
        "n_inner_shapes": len(inner),
        "inner_envelope": inner_row,
        "clearances": {key: float(value) for key, value in clearances.items()},
        "min_clearance": float(min_clearance),
        "contained_by_bbox": bool(min_clearance >= -1.0e-12),
        "enclosure_volume": enclosure_volume,
        "inner_volume_sum": inner_volume,
        "nominal_void_volume": enclosure_volume - inner_volume,
        "inner_volume_fraction": inner_volume / enclosure_volume if enclosure_volume else math.inf,
    }


def enclosure_difference_region(enclosure, inner_shapes, label="air"):
    """Return ``enclosure - inner_shapes`` as a labelled region.

    This is the build123d-side companion to the Netgen/Gmsh multi-region
    contract: keep inner material solids as labelled children, and make the
    surrounding region explicitly disjoint by subtracting them from the outer
    enclosure before meshing.
    """

    region = enclosure
    for shape in _shape_list(inner_shapes):
        region = region - shape
    region = region.solid()
    region.label = label
    return region


def shape_measurement_row(shape, name=None, index=1):
    """Return a JSON-friendly build123d measurement row for one shape.

    The row is intentionally close to Cubit's geometry API vocabulary:
    volume, surface area, topology counts, and bounding-box coordinates.  It is
    useful both as a quick sanity report for generated CAD and as the build123d
    side of a STEP round-trip cross validation against an external geometry
    kernel.
    """

    label = name or getattr(shape, "label", "") or f"shape_{index}"
    bb = shape.bounding_box()
    faces = shape.faces()
    edges = shape.edges()
    vertices = shape.vertices()
    solids = shape.solids()
    is_valid = bool(all(s.is_valid for s in solids)) if isinstance(shape, Compound) else bool(shape.is_valid)
    bbox_min = [float(bb.min.X), float(bb.min.Y), float(bb.min.Z)]
    bbox_max = [float(bb.max.X), float(bb.max.Y), float(bb.max.Z)]
    bbox_size = [float(bb.size.X), float(bb.size.Y), float(bb.size.Z)]
    bbox_center = [(lo + hi) / 2.0 for lo, hi in zip(bbox_min, bbox_max)]
    bbox_diagonal = math.sqrt(sum(value * value for value in bbox_size))
    return {
        "index": int(index),
        "name": str(label),
        "type": type(shape).__name__,
        "is_valid": is_valid,
        "volume": float(shape.volume),
        "area": float(shape.area),
        "faces": len(faces),
        "edges": len(edges),
        "vertices": len(vertices),
        "solids": len(solids),
        "bounding_box": {
            "min": bbox_min,
            "max": bbox_max,
            "center": bbox_center,
            "size": bbox_size,
            "diagonal": bbox_diagonal,
        },
        "characteristic_length": max(bbox_size),
    }


def shape_measurement_rows(shapes):
    """Return measurement rows for labelled shapes or ``(shape, name)`` tuples.

    If a single labelled :class:`~build123d.Compound` with children is passed,
    the children are measured as separate rows.  That mirrors multi-region CAE
    assemblies where each child maps to a material/block.
    """

    if isinstance(shapes, Compound):
        items = list(shapes.children) if shapes.children else list(shapes.solids())
    else:
        items = list(shapes)

    rows = []
    for index, item in enumerate(items, start=1):
        if isinstance(item, tuple) and len(item) == 2:
            shape, name = item
        else:
            shape, name = item, None
        rows.append(shape_measurement_row(shape, name=name, index=index))
    return rows


def box_face_vector_area_rows(size, center=(0.0, 0.0, 0.0), names=None):
    """Return analytic oriented area vectors for the six faces of an axis-aligned box.

    The row vocabulary mirrors
    ``radia_mcp.radia_ngsolve.netgen_vol.NetgenTriTetVolMesh.boundary_normal_summary_rows``:
    each face has a scalar ``surface_area``, an oriented ``vector_area`` equal
    to ``normal * area``, a unit normal, and the face center.  This gives CAD
    scripts a compact reference for checking boundary orientation before
    Maxwell-stress or pressure loads are integrated on a surface mesh.
    """

    sx, sy, sz = [float(value) for value in size]
    cx, cy, cz = [float(value) for value in center]
    if sx <= 0.0 or sy <= 0.0 or sz <= 0.0:
        raise ValueError("box size values must be positive")

    default_names = ["xmin", "xmax", "ymin", "ymax", "zmin", "zmax"]
    if names is None:
        labels = default_names
    elif isinstance(names, dict):
        labels = [str(names.get(name, name)) for name in default_names]
    else:
        labels = [str(name) for name in names]
        if len(labels) != 6:
            raise ValueError("names must have exactly six entries")

    specs = [
        ((-1.0, 0.0, 0.0), sy * sz, (cx - 0.5 * sx, cy, cz)),
        ((1.0, 0.0, 0.0), sy * sz, (cx + 0.5 * sx, cy, cz)),
        ((0.0, -1.0, 0.0), sx * sz, (cx, cy - 0.5 * sy, cz)),
        ((0.0, 1.0, 0.0), sx * sz, (cx, cy + 0.5 * sy, cz)),
        ((0.0, 0.0, -1.0), sx * sy, (cx, cy, cz - 0.5 * sz)),
        ((0.0, 0.0, 1.0), sx * sy, (cx, cy, cz + 0.5 * sz)),
    ]

    rows = []
    for index, (label, (unit_normal, area, face_center)) in enumerate(zip(labels, specs), start=1):
        vector_area = tuple(component * area for component in unit_normal)
        vector_norm = math.sqrt(sum(component * component for component in vector_area))
        rows.append({
            "index": index,
            "name": label,
            "surface_area": float(area),
            "vector_area": vector_area,
            "vector_area_norm": float(vector_norm),
            "vector_area_norm_over_area": vector_norm / area,
            "unit_normal": tuple(float(component) for component in unit_normal),
            "face_center": tuple(float(component) for component in face_center),
            "box_center": (cx, cy, cz),
            "box_size": (sx, sy, sz),
        })
    return rows


def box_face_pressure_force_rows(
    size,
    pressure_by_face,
    center=(0.0, 0.0, 0.0),
    names=None,
    default_pressure=0.0,
):
    """Return analytic box face force rows from scalar pressure values.

    Pressure is positive along each face's outward normal:
    ``force = pressure * vector_area``.  This is the build123d-side analytic
    companion to ``NetgenTriTetVolMesh.boundary_pressure_force_rows``.
    """

    rows = []
    for area_row in box_face_vector_area_rows(size, center=center, names=names):
        name = area_row["name"]
        if name in pressure_by_face:
            pressure = float(pressure_by_face[name])
            source = "name"
        elif area_row["index"] in pressure_by_face:
            pressure = float(pressure_by_face[area_row["index"]])
            source = "index"
        elif default_pressure is not None:
            pressure = float(default_pressure)
            source = "default"
        else:
            raise KeyError(f"missing pressure for face {name}")
        force = tuple(pressure * component for component in area_row["vector_area"])
        rows.append({
            **area_row,
            "pressure_Pa": pressure,
            "pressure_source": source,
            "force_N": force,
            "force_magnitude_N": math.sqrt(sum(component * component for component in force)),
        })
    return rows


def _cross3(a, b):
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def box_face_pressure_moment_rows(
    size,
    pressure_by_face,
    center=(0.0, 0.0, 0.0),
    names=None,
    default_pressure=0.0,
    pivot_m=(0.0, 0.0, 0.0),
):
    """Return analytic box face pressure force and pivot-moment rows.

    This is the build123d-side analytic companion to
    ``NetgenTriTetVolMesh.boundary_pressure_force_moment_rows``.  Each face is
    planar with a constant pressure, so its resultant acts at ``face_center``:

    ``F = pressure * vector_area`` and ``M = (face_center - pivot) x F``.
    """

    pivot = tuple(float(value) for value in pivot_m)
    if len(pivot) != 3:
        raise ValueError("pivot_m must have three components")
    rows = []
    for row in box_face_pressure_force_rows(
        size,
        pressure_by_face,
        center=center,
        names=names,
        default_pressure=default_pressure,
    ):
        lever = tuple(row["face_center"][axis] - pivot[axis] for axis in range(3))
        moment = _cross3(lever, row["force_N"])
        rows.append({
            **row,
            "pivot_m": pivot,
            "lever_arm_m": lever,
            "moment_about_pivot_Nm": moment,
            "moment_magnitude_Nm": math.sqrt(sum(component * component for component in moment)),
        })
    return rows


def box_face_pressure_resultant_summary(
    size,
    pressure_by_face,
    center=(0.0, 0.0, 0.0),
    names=None,
    default_pressure=0.0,
    pivot_m=(0.0, 0.0, 0.0),
):
    """Return analytic box-face pressure rows plus net force/moment metrics.

    The return vocabulary mirrors
    ``NetgenTriTetVolMesh.boundary_pressure_resultant_summary`` so a CAD-side
    build123d box reference and a Cubit/Coreform-exported `.vol` boundary mesh
    can be compared without adapter code.
    """

    rows = box_face_pressure_moment_rows(
        size,
        pressure_by_face,
        center=center,
        names=names,
        default_pressure=default_pressure,
        pivot_m=pivot_m,
    )
    total_force = tuple(sum(row["force_N"][axis] for row in rows) for axis in range(3))
    total_moment = tuple(
        sum(row["moment_about_pivot_Nm"][axis] for row in rows)
        for axis in range(3)
    )
    total_force_norm = math.sqrt(sum(component * component for component in total_force))
    total_moment_norm = math.sqrt(sum(component * component for component in total_moment))
    absolute_force_sum = sum(row["force_magnitude_N"] for row in rows)
    absolute_moment_sum = sum(row["moment_magnitude_Nm"] for row in rows)
    surface_vector_area = tuple(
        sum(row["vector_area"][axis] for row in rows)
        for axis in range(3)
    )
    surface_vector_area_norm = math.sqrt(
        sum(component * component for component in surface_vector_area)
    )
    total_area = sum(row["surface_area"] for row in rows)

    return {
        "boundary_count": len(rows),
        "rows": rows,
        "pivot_m": tuple(float(value) for value in pivot_m),
        "total_force_N": total_force,
        "total_force_magnitude_N": total_force_norm,
        "total_moment_about_pivot_Nm": total_moment,
        "total_moment_magnitude_Nm": total_moment_norm,
        "absolute_force_sum_N": absolute_force_sum,
        "absolute_moment_sum_Nm": absolute_moment_sum,
        "force_balance_ratio": (
            total_force_norm / absolute_force_sum
            if absolute_force_sum > 0.0
            else 0.0
        ),
        "moment_balance_ratio": (
            total_moment_norm / absolute_moment_sum
            if absolute_moment_sum > 0.0
            else 0.0
        ),
        "surface_vector_area": surface_vector_area,
        "surface_vector_area_norm": surface_vector_area_norm,
        "surface_vector_area_norm_over_area": (
            surface_vector_area_norm / total_area if total_area > 0.0 else None
        ),
        "box_center": tuple(float(value) for value in center),
        "box_size": tuple(float(value) for value in size),
    }


def _face_vector_value(vector_by_face, area_row, default_vector, value_name):
    name = area_row["name"]
    if name in vector_by_face:
        value = vector_by_face[name]
        source = "name"
    elif area_row["index"] in vector_by_face:
        value = vector_by_face[area_row["index"]]
        source = "index"
    elif default_vector is not None:
        value = default_vector
        source = "default"
    else:
        raise KeyError(f"missing {value_name} for face {name}")
    vector = tuple(float(component) for component in value)
    if len(vector) != 3:
        raise ValueError(f"{value_name} values must have three components")
    return vector, source


def box_face_traction_moment_rows(
    size,
    traction_by_face,
    center=(0.0, 0.0, 0.0),
    names=None,
    default_traction=(0.0, 0.0, 0.0),
    pivot_m=(0.0, 0.0, 0.0),
):
    """Return analytic box face vector-traction force and pivot-moment rows.

    The traction vector is a constant global vector [N/m2] on each planar box
    face, so ``F = traction * area`` and ``M = (face_center - pivot) x F``.
    Unlike scalar pressure, the force direction is not tied to the face normal.
    This is the build123d-side analytic companion to
    ``NetgenTriTetVolMesh.boundary_traction_force_moment_rows``.
    """

    pivot = tuple(float(value) for value in pivot_m)
    if len(pivot) != 3:
        raise ValueError("pivot_m must have three components")
    rows = []
    for area_row in box_face_vector_area_rows(size, center=center, names=names):
        traction, source = _face_vector_value(
            traction_by_face,
            area_row,
            default_traction,
            "traction",
        )
        force = tuple(float(area_row["surface_area"]) * component for component in traction)
        lever = tuple(area_row["face_center"][axis] - pivot[axis] for axis in range(3))
        moment = _cross3(lever, force)
        rows.append({
            **area_row,
            "traction_N_per_m2": traction,
            "traction_source": source,
            "force_N": force,
            "force_magnitude_N": math.sqrt(sum(component * component for component in force)),
            "pivot_m": pivot,
            "lever_arm_m": lever,
            "moment_about_pivot_Nm": moment,
            "moment_magnitude_Nm": math.sqrt(sum(component * component for component in moment)),
        })
    return rows


def _surface_area_from_row(row):
    if "surface_area" in row:
        return float(row["surface_area"])
    if "area" in row:
        return float(row["area"])
    return None


def _row_vector(row, key):
    values = row.get(key)
    if values is None:
        return None
    vector = list(values)
    if len(vector) != 3:
        return None
    return [float(value) for value in vector]


def _vector_norm(values):
    return math.sqrt(sum(float(value) * float(value) for value in values))


def compare_boundary_vector_area_rows(
    reference_rows,
    measured_rows,
    vector_atol=1.0e-9,
    area_rtol=1.0e-9,
    measured_label="measured",
):
    """Compare oriented boundary-area rows by name.

    ``reference_rows`` can come from :func:`box_face_vector_area_rows`, while
    ``measured_rows`` can come from a triangle/tetrahedral ``.vol`` boundary
    summary.  The comparison checks scalar area, oriented area vector, and unit
    normal when both sides provide one.
    """

    measured_by_name = {row["name"]: row for row in measured_rows}
    rows = []
    for ref in reference_rows:
        name = ref["name"]
        measured = measured_by_name.get(name)
        ref_area = _surface_area_from_row(ref)
        ref_vector = _row_vector(ref, "vector_area")
        ref_normal = _row_vector(ref, "unit_normal")
        if measured is None:
            rows.append({
                "name": name,
                "measured_label": measured_label,
                "reference_surface_area": ref_area,
                "measured_surface_area": None,
                "area_rel_error": None,
                "reference_vector_area": ref_vector,
                "measured_vector_area": None,
                "vector_abs_error": None,
                "reference_unit_normal": ref_normal,
                "measured_unit_normal": None,
                "unit_normal_abs_error": None,
                "vector_atol": float(vector_atol),
                "area_rtol": float(area_rtol),
                "passed": False,
                "reason": "missing measured row",
            })
            continue

        measured_area = _surface_area_from_row(measured)
        measured_vector = _row_vector(measured, "vector_area")
        measured_normal = _row_vector(measured, "unit_normal")
        area_rel_error = (
            _relative_measurement_error(ref_area, measured_area)
            if ref_area is not None and measured_area is not None
            else None
        )
        vector_abs_error = (
            _vector_norm(float(a) - float(b) for a, b in zip(ref_vector, measured_vector))
            if ref_vector is not None and measured_vector is not None
            else None
        )
        unit_normal_abs_error = (
            _vector_norm(float(a) - float(b) for a, b in zip(ref_normal, measured_normal))
            if ref_normal is not None and measured_normal is not None
            else None
        )
        area_ok = area_rel_error is not None and area_rel_error <= area_rtol
        vector_ok = vector_abs_error is not None and vector_abs_error <= vector_atol
        normal_ok = unit_normal_abs_error is None or unit_normal_abs_error <= vector_atol
        passed = bool(area_ok and vector_ok and normal_ok)
        rows.append({
            "name": name,
            "measured_label": measured_label,
            "reference_surface_area": ref_area,
            "measured_surface_area": measured_area,
            "area_rel_error": area_rel_error,
            "reference_vector_area": ref_vector,
            "measured_vector_area": measured_vector,
            "vector_abs_error": vector_abs_error,
            "reference_unit_normal": ref_normal,
            "measured_unit_normal": measured_normal,
            "unit_normal_abs_error": unit_normal_abs_error,
            "vector_atol": float(vector_atol),
            "area_rtol": float(area_rtol),
            "passed": passed,
            "reason": "ok" if passed else "outside tolerance",
        })
    return rows


def _relative_measurement_error(reference, measured):
    return abs(float(reference) - float(measured)) / max(
        abs(float(reference)),
        abs(float(measured)),
        1.0e-300,
    )


def _xyz(values):
    if values is None:
        return None
    xyz = list(values)
    if len(xyz) != 3:
        return None
    return [float(value) for value in xyz]


def _row_bounding_box(row):
    if not isinstance(row, dict):
        return None
    box = row.get("bounding_box")
    if isinstance(box, dict):
        bbox_min = _xyz(box.get("min"))
        bbox_max = _xyz(box.get("max"))
        bbox_center = _xyz(box.get("center"))
        bbox_size = _xyz(box.get("size"))
        bbox_diagonal = box.get("diagonal")
    else:
        bbox_min = _xyz(row.get("bbox_min"))
        bbox_max = _xyz(row.get("bbox_max"))
        bbox_center = _xyz(row.get("bbox_center"))
        bbox_size = _xyz(row.get("bbox_size"))
        bbox_diagonal = row.get("bbox_diagonal")

    if bbox_min is None and bbox_center is not None and bbox_size is not None:
        bbox_min = [center - 0.5 * size for center, size in zip(bbox_center, bbox_size)]
    if bbox_max is None and bbox_center is not None and bbox_size is not None:
        bbox_max = [center + 0.5 * size for center, size in zip(bbox_center, bbox_size)]
    if bbox_min is None or bbox_max is None:
        return None
    if bbox_size is None:
        bbox_size = [hi - lo for lo, hi in zip(bbox_min, bbox_max)]
    if bbox_center is None:
        bbox_center = [(lo + hi) / 2.0 for lo, hi in zip(bbox_min, bbox_max)]
    if bbox_diagonal is None:
        bbox_diagonal = math.sqrt(sum(value * value for value in bbox_size))

    return {
        "min": bbox_min,
        "max": bbox_max,
        "center": bbox_center,
        "size": bbox_size,
        "diagonal": float(bbox_diagonal),
    }


def _max_abs_delta(reference, measured):
    return max(abs(float(a) - float(b)) for a, b in zip(reference, measured))


def shape_measurement_inventory_summary(rows):
    """Return assembly-level inventory statistics from measurement rows.

    ``shape_measurement_rows`` is intentionally row-oriented.  This companion
    answers the next question a CAE user asks: how much volume each labelled
    region contributes, how much area is present, and how tightly the labelled
    regions fill their union bounding box.
    """

    rows = list(rows)
    if not rows:
        return {
            "n_shapes": 0,
            "n_valid": 0,
            "total_volume": 0.0,
            "total_area": 0.0,
            "total_solids": 0,
            "total_faces": 0,
            "total_edges": 0,
            "total_vertices": 0,
            "bounding_box": None,
            "bbox_volume": None,
            "bbox_fill_fraction": None,
            "volume_fraction_rows": [],
            "largest_volume_name": None,
            "smallest_volume_name": None,
        }

    total_volume = sum(float(row["volume"]) for row in rows)
    total_area = sum(float(row["area"]) for row in rows)
    boxes = [_row_bounding_box(row) for row in rows]
    complete_bbox = all(box is not None for box in boxes)
    bounding_box = None
    bbox_volume = None
    fill_fraction = None
    if complete_bbox:
        mins = [min(box["min"][axis] for box in boxes) for axis in range(3)]
        maxs = [max(box["max"][axis] for box in boxes) for axis in range(3)]
        size = [maxs[axis] - mins[axis] for axis in range(3)]
        center = [(mins[axis] + maxs[axis]) / 2.0 for axis in range(3)]
        bbox_volume = size[0] * size[1] * size[2]
        fill_fraction = total_volume / bbox_volume if bbox_volume > 0.0 else None
        bounding_box = {
            "min": [float(value) for value in mins],
            "max": [float(value) for value in maxs],
            "center": [float(value) for value in center],
            "size": [float(value) for value in size],
            "diagonal": math.sqrt(sum(value * value for value in size)),
        }

    volume_fraction_rows = []
    for row in rows:
        volume = float(row["volume"])
        area = float(row["area"])
        volume_fraction_rows.append({
            "name": str(row["name"]),
            "volume": volume,
            "volume_fraction": volume / total_volume if total_volume > 0.0 else None,
            "area": area,
            "area_fraction": area / total_area if total_area > 0.0 else None,
        })

    by_volume = sorted(volume_fraction_rows, key=lambda row: row["volume"])
    return {
        "n_shapes": len(rows),
        "n_valid": sum(1 for row in rows if bool(row.get("is_valid", True))),
        "total_volume": total_volume,
        "total_area": total_area,
        "total_solids": sum(int(row.get("solids", 0)) for row in rows),
        "total_faces": sum(int(row.get("faces", 0)) for row in rows),
        "total_edges": sum(int(row.get("edges", 0)) for row in rows),
        "total_vertices": sum(int(row.get("vertices", 0)) for row in rows),
        "bounding_box": bounding_box,
        "bbox_volume": bbox_volume,
        "bbox_fill_fraction": fill_fraction,
        "volume_fraction_rows": volume_fraction_rows,
        "largest_volume_name": by_volume[-1]["name"],
        "smallest_volume_name": by_volume[0]["name"],
    }


def _bbox_pair_clearance_row(row_a, box_a, row_b, box_b, tolerance):
    axes = ("x", "y", "z")
    axis_gaps = {}
    axis_overlaps = {}
    separated_axes = []
    touching_axes = []
    for axis, label in enumerate(axes):
        lo_a = float(box_a["min"][axis])
        hi_a = float(box_a["max"][axis])
        lo_b = float(box_b["min"][axis])
        hi_b = float(box_b["max"][axis])
        if hi_a < lo_b:
            gap = lo_b - hi_a
        elif hi_b < lo_a:
            gap = lo_a - hi_b
        else:
            gap = 0.0
        overlap = min(hi_a, hi_b) - max(lo_a, lo_b)
        if gap > tolerance:
            separated_axes.append(label)
        elif abs(overlap) <= tolerance:
            touching_axes.append(label)
        axis_gaps[label] = float(gap if gap > 0.0 else 0.0)
        axis_overlaps[label] = float(overlap if overlap > 0.0 else 0.0)

    bbox_separated = bool(separated_axes)
    if bbox_separated:
        status = "separated"
    elif touching_axes:
        status = "touching_bbox"
    else:
        status = "bbox_overlap_needs_precise_check"
    intersection_size = [axis_overlaps[label] for label in axes]
    intersection_volume = (
        intersection_size[0] * intersection_size[1] * intersection_size[2]
        if not bbox_separated else 0.0
    )
    positive_gaps = [value for value in axis_gaps.values() if value > tolerance]
    center_distance = math.sqrt(sum(
        (float(box_a["center"][axis]) - float(box_b["center"][axis])) ** 2
        for axis in range(3)
    ))
    name_a = str(row_a.get("name", "shape_a"))
    name_b = str(row_b.get("name", "shape_b"))
    return {
        "pair": f"{name_a}::{name_b}",
        "name_a": name_a,
        "name_b": name_b,
        "status": status,
        "bbox_separated": bbox_separated,
        "bbox_intersects_or_touches": not bbox_separated,
        "separated_axes": separated_axes,
        "touching_axes": touching_axes,
        "axis_gaps": axis_gaps,
        "axis_overlaps": axis_overlaps,
        "separation_distance": min(positive_gaps) if positive_gaps else 0.0,
        "center_distance": float(center_distance),
        "bbox_intersection_size": intersection_size,
        "bbox_intersection_volume": float(intersection_volume),
    }


def shape_bbox_pair_clearance_summary(rows, clearance_tolerance=1.0e-12):
    """Audit pairwise bounding-box clearances from measurement rows.

    A positive gap on any one axis is a cheap proof that two shapes are
    separated before STEP/mesh export.  If all three bbox intervals overlap the
    pair is not declared intersecting; it is marked for precise geometry or
    boolean checking.
    """

    rows = list(rows)
    tolerance = float(clearance_tolerance)
    if tolerance < 0.0:
        raise ValueError("clearance_tolerance must be >= 0")
    complete = []
    missing = []
    for index, row in enumerate(rows, start=1):
        box = _row_bounding_box(row)
        if box is None:
            missing.append(str(row.get("name", f"shape_{index}") if isinstance(row, dict) else f"shape_{index}"))
        else:
            complete.append((row, box))

    pair_rows = []
    for i, (row_a, box_a) in enumerate(complete):
        for row_b, box_b in complete[i + 1:]:
            pair_rows.append(_bbox_pair_clearance_row(row_a, box_a, row_b, box_b, tolerance))

    separated_count = sum(1 for row in pair_rows if row["status"] == "separated")
    touching_count = sum(1 for row in pair_rows if row["status"] == "touching_bbox")
    overlap_count = sum(1 for row in pair_rows if row["status"] == "bbox_overlap_needs_precise_check")
    gaps = [row["separation_distance"] for row in pair_rows if row["separation_distance"] > tolerance]
    overlap_volumes = [
        row["bbox_intersection_volume"]
        for row in pair_rows
        if row["status"] == "bbox_overlap_needs_precise_check"
    ]
    ok = not missing and touching_count == 0 and overlap_count == 0
    return {
        "policy": "build123d_bbox_pair_clearance_pre_mesh_audit",
        "status": "ok" if ok else "needs_attention",
        "ok_for_bbox_clearance": ok,
        "clearance_tolerance": tolerance,
        "n_shapes": len(rows),
        "n_complete_bbox_rows": len(complete),
        "missing_bbox_names": missing,
        "n_pairs": len(pair_rows),
        "separated_pair_count": separated_count,
        "touching_pair_count": touching_count,
        "bbox_overlap_pair_count": overlap_count,
        "min_positive_gap": min(gaps) if gaps else None,
        "max_bbox_intersection_volume": max(overlap_volumes) if overlap_volumes else 0.0,
        "pair_rows": pair_rows,
    }


def compare_shape_measurement_rows(
    reference_rows,
    measured_rows,
    rtol=1.0e-5,
    measured_label="measured",
    bbox_atol=1.0e-6,
):
    """Compare build123d measurement rows with external measurement rows.

    ``reference_rows`` are normally from :func:`shape_measurement_rows`.
    ``measured_rows`` only need ``name``, ``volume`` and ``area`` keys.  If both
    sides also provide ``bounding_box`` data, the bounding box is compared with
    ``bbox_atol``.  Rows can come from Cubit, another CAD kernel, a mesher, or
    an analytic table.  The return value is a list of JSON-friendly pass/fail
    rows.
    """

    measured_by_name = {row["name"]: row for row in measured_rows}
    rows = []
    for ref in reference_rows:
        name = ref["name"]
        measured = measured_by_name.get(name)
        ref_bbox = _row_bounding_box(ref)
        if measured is None:
            rows.append({
                "name": name,
                "measured_label": measured_label,
                "reference_volume": ref["volume"],
                "reference_area": ref["area"],
                "reference_bounding_box": ref_bbox,
                "measured_volume": None,
                "measured_area": None,
                "measured_bounding_box": None,
                "volume_rel_error": None,
                "area_rel_error": None,
                "bbox_compared": False,
                "bbox_min_abs_error": None,
                "bbox_max_abs_error": None,
                "bbox_center_abs_error": None,
                "bbox_size_abs_error": None,
                "bbox_abs_error": None,
                "rtol": float(rtol),
                "bbox_atol": float(bbox_atol),
                "passed": False,
                "reason": "missing measured row",
            })
            continue

        volume_rel_error = _relative_measurement_error(ref["volume"], measured["volume"])
        area_rel_error = _relative_measurement_error(ref["area"], measured["area"])
        measured_bbox = _row_bounding_box(measured)
        bbox_compared = ref_bbox is not None and measured_bbox is not None
        bbox_min_abs_error = None
        bbox_max_abs_error = None
        bbox_center_abs_error = None
        bbox_size_abs_error = None
        bbox_abs_error = None
        bbox_passed = True
        if bbox_compared:
            bbox_min_abs_error = _max_abs_delta(ref_bbox["min"], measured_bbox["min"])
            bbox_max_abs_error = _max_abs_delta(ref_bbox["max"], measured_bbox["max"])
            bbox_center_abs_error = _max_abs_delta(ref_bbox["center"], measured_bbox["center"])
            bbox_size_abs_error = _max_abs_delta(ref_bbox["size"], measured_bbox["size"])
            bbox_abs_error = max(
                bbox_min_abs_error,
                bbox_max_abs_error,
                bbox_center_abs_error,
                bbox_size_abs_error,
            )
            bbox_passed = bbox_abs_error <= bbox_atol

        volume_area_passed = volume_rel_error <= rtol and area_rel_error <= rtol
        passed = volume_area_passed and bbox_passed
        if passed:
            reason = "ok"
        elif not volume_area_passed:
            reason = "outside tolerance"
        else:
            reason = "bbox outside tolerance"
        rows.append({
            "name": name,
            "measured_label": measured_label,
            "reference_volume": ref["volume"],
            "reference_area": ref["area"],
            "reference_bounding_box": ref_bbox,
            "measured_volume": measured["volume"],
            "measured_area": measured["area"],
            "measured_bounding_box": measured_bbox,
            "volume_rel_error": volume_rel_error,
            "area_rel_error": area_rel_error,
            "bbox_compared": bbox_compared,
            "bbox_min_abs_error": bbox_min_abs_error,
            "bbox_max_abs_error": bbox_max_abs_error,
            "bbox_center_abs_error": bbox_center_abs_error,
            "bbox_size_abs_error": bbox_size_abs_error,
            "bbox_abs_error": bbox_abs_error,
            "rtol": float(rtol),
            "bbox_atol": float(bbox_atol),
            "passed": passed,
            "reason": reason,
        })
    return rows


def _volume_row_name(row, index):
    if isinstance(row, dict):
        return str(row.get("name") or row.get("label") or f"shape_{index}")
    return f"shape_{index}"


def _normalize_volume_rows(rows):
    if isinstance(rows, dict):
        if "rows" in rows:
            rows = rows["rows"]
        elif "measurements" in rows:
            rows = rows["measurements"]
        elif "volume" in rows:
            rows = [rows]
        else:
            rows = [{"name": key, "volume": value} for key, value in rows.items()]
    normalized = []
    for index, row in enumerate(list(rows), start=1):
        if isinstance(row, dict):
            if "volume" not in row:
                raise KeyError(f"missing volume for row {index}")
            volume = float(row["volume"])
        else:
            volume = float(row)
        if not math.isfinite(volume):
            raise ValueError(f"volume for row {index} must be finite")
        normalized.append({
            "index": index,
            "name": _volume_row_name(row, index),
            "volume": volume,
        })
    return normalized


def compare_shape_volume_rows(
    reference_rows,
    measured_rows,
    rtol=1.0e-5,
    measured_label="external_cad",
):
    """Compare shape volumes by name against an external CAD/kernel table.

    This is intentionally weaker than :func:`compare_shape_measurement_rows`:
    it only requires ``name`` and ``volume``.  That makes it the common-denominator
    check for build123d -> STEP round trips through Cubit or another CAD system
    whose scripting API can report body volumes but not identical topology/area
    rows.  Use the full measurement comparison when area and bounding boxes are
    also available.
    """

    reference = _normalize_volume_rows(reference_rows)
    measured = _normalize_volume_rows(measured_rows)
    measured_by_name = {row["name"]: row for row in measured}
    rows = []
    for ref in reference:
        measured_row = measured_by_name.get(ref["name"])
        if measured_row is None:
            rows.append({
                "name": ref["name"],
                "measured_label": measured_label,
                "reference_volume": ref["volume"],
                "measured_volume": None,
                "volume_rel_error": None,
                "rtol": float(rtol),
                "passed": False,
                "reason": "missing measured row",
            })
            continue
        rel_error = _relative_measurement_error(ref["volume"], measured_row["volume"])
        passed = rel_error <= float(rtol)
        rows.append({
            "name": ref["name"],
            "measured_label": measured_label,
            "reference_volume": ref["volume"],
            "measured_volume": measured_row["volume"],
            "volume_rel_error": rel_error,
            "rtol": float(rtol),
            "passed": passed,
            "reason": "ok" if passed else "outside tolerance",
        })
    return rows


def _normalize_volume_measurement_sets(measured_sets):
    if isinstance(measured_sets, dict):
        if "rows" in measured_sets or "volume" in measured_sets:
            meta = {key: value for key, value in measured_sets.items() if key != "rows"}
            return [("external_cad", measured_sets, meta)]
        return [
            (
                str(label),
                rows,
                (
                    {key: value for key, value in rows.items() if key != "rows"}
                    if isinstance(rows, dict) and "rows" in rows
                    else {}
                ),
            )
            for label, rows in measured_sets.items()
        ]
    normalized = []
    for index, item in enumerate(list(measured_sets), start=1):
        if isinstance(item, dict) and "rows" in item:
            label = str(item.get("source") or item.get("label") or f"external_cad_{index}")
            rows = item["rows"]
            meta = {key: value for key, value in item.items() if key not in {"rows", "source", "label"}}
        else:
            label = f"external_cad_{index}"
            rows = item
            meta = {}
        normalized.append((label, rows, meta))
    return normalized


def shape_volume_crosscheck_summary(reference_rows, measured_sets, rtol=1.0e-5):
    """Return a volume-only crosscheck summary for one or more CAD sources.

    ``reference_rows`` are usually build123d measurement rows. ``measured_sets``
    can be ``{"cubit": rows, "external_cad": rows}`` or a list of
    ``{"source": label, "rows": rows}`` dictionaries.  The summary is designed
    for cross-validation artifacts where volume is the stable common contract
    across CAD kernels; it does not publish or assume any private-tool
    provenance.
    """

    reference = _normalize_volume_rows(reference_rows)
    sets = []
    all_rows = []
    for label, rows, metadata in _normalize_volume_measurement_sets(measured_sets):
        compared = compare_shape_volume_rows(
            reference,
            rows,
            rtol=rtol,
            measured_label=label,
        )
        passed = sum(1 for row in compared if row["passed"])
        errors = [row["volume_rel_error"] or 0.0 for row in compared]
        source_summary = {
            "source": label,
            "n_cases": len(compared),
            "n_passed": passed,
            "max_volume_rel_error": max(errors) if errors else 0.0,
            "status": "ok" if passed == len(compared) else "needs_attention",
            "rows": compared,
        }
        for key in (
            "source_artifact_id",
            "source_artifact_digest",
            "measurement_method",
            "body_identity_key",
            "volume_unit",
            "cad_kernel",
            "parameter_set_artifact_id",
            "parameter_set_digest",
            "parameter_set_path",
            "objective_observable_id",
            "objective_observable_family",
        ):
            if key in metadata and metadata[key] is not None:
                source_summary[key] = str(metadata[key]).strip()
        sets.append(source_summary)
        all_rows.extend(compared)

    failed = [row for row in all_rows if not row["passed"]]
    all_errors = [row["volume_rel_error"] or 0.0 for row in all_rows]
    return {
        "policy": "build123d_external_cad_volume_crosscheck",
        "status": "ok" if not failed else "needs_attention",
        "ok_for_cad_roundtrip_volume": not failed,
        "rtol": float(rtol),
        "n_reference_rows": len(reference),
        "n_sources": len(sets),
        "n_failed_rows": len(failed),
        "max_volume_rel_error": max(all_errors) if all_errors else 0.0,
        "sources": [item["source"] for item in sets],
        "comparison_sets": sets,
    }


def shape_perforated_prism_roundtrip_gate(
    *,
    reference_volume,
    imported_volume,
    hole_count,
    hole_side_count,
    imported_surface_count,
    imported_body_count=1,
    outer_side_count=4,
    volume_rtol=1.0e-9,
):
    """Gate a perforated prism by volume and expected boundary topology."""

    reference = float(reference_volume)
    imported = float(imported_volume)
    holes = int(hole_count)
    hole_sides = int(hole_side_count)
    imported_surfaces = int(imported_surface_count)
    imported_bodies = int(imported_body_count)
    outer_sides = int(outer_side_count)
    tolerance = float(volume_rtol)
    if reference <= 0.0 or imported < 0.0:
        raise ValueError("reference_volume must be > 0 and imported_volume must be >= 0")
    if holes < 0 or hole_sides < 3 or outer_sides < 3:
        raise ValueError("hole_count must be >= 0 and polygon side counts must be >= 3")
    if imported_surfaces < 0 or imported_bodies < 0:
        raise ValueError("imported counts must be non-negative")
    if tolerance < 0.0:
        raise ValueError("volume_rtol must be >= 0")

    expected_surfaces = 2 + outer_sides + holes * hole_sides
    volume_relative_error = abs(imported - reference) / reference
    checks = {
        "volume_agrees": volume_relative_error <= tolerance,
        "single_body_preserved": imported_bodies == 1,
        "surface_topology_preserved": imported_surfaces == expected_surfaces,
    }
    return {
        "policy": "build123d_perforated_prism_roundtrip_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "reference_volume": reference,
        "imported_volume": imported,
        "volume_relative_error": volume_relative_error,
        "volume_rtol": tolerance,
        "hole_count": holes,
        "hole_side_count": hole_sides,
        "outer_side_count": outer_sides,
        "expected_surface_count": expected_surfaces,
        "imported_surface_count": imported_surfaces,
        "imported_body_count": imported_bodies,
        "checks": checks,
        "lesson": (
            "Volume agreement is necessary but not sufficient for a perforated CAD "
            "roundtrip; require the expected cap, outer-wall, and hole-wall faces."
        ),
    }


def shape_volume_crosscheck_source_coverage_gate(
    volume_summary,
    required_sources=("cubit", "cst_import"),
    max_allowed_volume_rel_error=None,
):
    """Check that external CAD volume crosscheck used the required sources.

    build123d can validate analytic CAD volume internally, but CAE handoff is
    stronger when at least Cubit and external CAD rows agree on the
    same named bodies.  This gate sits on top of
    ``shape_volume_crosscheck_summary`` and prevents a build123d-only or
    single-kernel replay from masquerading as a multi-source CAD crosscheck.
    """

    summary = dict(volume_summary)
    required = [str(source).strip() for source in required_sources if str(source).strip()]
    sources = [str(source).strip() for source in summary.get("sources", [])]
    source_set = set(sources)
    missing_sources = [source for source in required if source not in source_set]
    max_error = float(summary.get("max_volume_rel_error", math.inf))
    allowed_error = (
        None if max_allowed_volume_rel_error is None else float(max_allowed_volume_rel_error)
    )
    checks = {
        "summary_policy_known": summary.get("policy") == "build123d_external_cad_volume_crosscheck",
        "summary_status_ok": summary.get("status") == "ok"
        and summary.get("ok_for_cad_roundtrip_volume") is True,
        "required_sources_present": not missing_sources,
        "source_count_sufficient": len(source_set) >= len(set(required)),
        "volume_error_within_limit": allowed_error is None or max_error <= allowed_error,
    }
    issues = [name for name, ok in checks.items() if not ok]
    return {
        "policy": "build123d_volume_crosscheck_source_coverage_gate",
        "status": "ok" if not issues else "needs_attention",
        "required_sources": required,
        "sources": sources,
        "missing_sources": missing_sources,
        "max_volume_rel_error": max_error,
        "max_allowed_volume_rel_error": allowed_error,
        "checks": checks,
        "issues": issues,
        "notes": [
            "Use this after shape_volume_crosscheck_summary when Cubit and external CAD rows are both required.",
            "Volume is the common CAD currency; solver-ready promotion still needs labels, files, and mesh/package gates.",
        ],
    }


def shape_volume_crosscheck_source_identity_gate(
    volume_summary,
    *,
    expected_measurement_methods=None,
    expected_body_identity_keys=None,
    expected_source_artifact_ids=None,
    expected_parameter_set_artifact_ids=None,
    expected_parameter_set_digests=None,
    expected_parameter_set_paths=None,
    expected_objective_observable_ids=None,
    expected_objective_observable_families=None,
):
    """Check source identity metadata for external CAD volume crosschecks.

    Volume values alone cannot prove which CAD body list, export artifact, or
    measurement API produced each source.  This gate sits after
    :func:`shape_volume_crosscheck_summary` and requires per-source method,
    body identity key, and artifact identity when the caller knows them.
    """

    summary = dict(volume_summary)
    sets = list(summary.get("comparison_sets", []) or [])

    def norm_map(value):
        if value is None:
            return {}
        return {str(key).strip(): str(val).strip() for key, val in dict(value).items()}

    expected_methods = norm_map(expected_measurement_methods)
    expected_body_keys = norm_map(expected_body_identity_keys)
    expected_artifacts = norm_map(expected_source_artifact_ids)
    expected_parameter_ids = norm_map(expected_parameter_set_artifact_ids)
    expected_parameter_digests = norm_map(expected_parameter_set_digests)
    expected_parameter_paths = norm_map(expected_parameter_set_paths)
    expected_objective_ids = norm_map(expected_objective_observable_ids)
    expected_objective_families = norm_map(expected_objective_observable_families)
    by_source = {str(item.get("source", "")).strip(): dict(item) for item in sets}

    def field_matches(expected, field):
        missing = []
        mismatched = []
        for source, expected_value in expected.items():
            row = by_source.get(source, {})
            actual = str(row.get(field, "") or "").strip()
            if not actual:
                missing.append(source)
            elif actual != expected_value:
                mismatched.append(source)
        return missing, mismatched

    missing_methods, mismatched_methods = field_matches(expected_methods, "measurement_method")
    missing_body_keys, mismatched_body_keys = field_matches(expected_body_keys, "body_identity_key")
    missing_artifacts, mismatched_artifacts = field_matches(expected_artifacts, "source_artifact_id")
    missing_parameter_ids, mismatched_parameter_ids = field_matches(
        expected_parameter_ids, "parameter_set_artifact_id"
    )
    missing_parameter_digests, mismatched_parameter_digests = field_matches(
        expected_parameter_digests, "parameter_set_digest"
    )
    missing_parameter_paths, mismatched_parameter_paths = field_matches(
        expected_parameter_paths, "parameter_set_path"
    )
    missing_objective_ids, mismatched_objective_ids = field_matches(
        expected_objective_ids, "objective_observable_id"
    )
    missing_objective_families, mismatched_objective_families = field_matches(
        expected_objective_families, "objective_observable_family"
    )
    expected_sources = (
        set(expected_methods)
        | set(expected_body_keys)
        | set(expected_artifacts)
        | set(expected_parameter_ids)
        | set(expected_parameter_digests)
        | set(expected_parameter_paths)
        | set(expected_objective_ids)
        | set(expected_objective_families)
    )
    missing_expected_sources = sorted(expected_sources - set(by_source))
    checks = {
        "summary_policy_known": summary.get("policy") == "build123d_external_cad_volume_crosscheck",
        "summary_status_ok": summary.get("status") == "ok",
        "expected_sources_present": not missing_expected_sources,
        "measurement_methods_recorded_when_expected": not missing_methods,
        "expected_measurement_methods_match": not mismatched_methods,
        "body_identity_keys_recorded_when_expected": not missing_body_keys,
        "expected_body_identity_keys_match": not mismatched_body_keys,
        "source_artifact_ids_recorded_when_expected": not missing_artifacts,
        "expected_source_artifact_ids_match": not mismatched_artifacts,
        "parameter_set_artifact_ids_recorded_when_expected": not missing_parameter_ids,
        "expected_parameter_set_artifact_ids_match": not mismatched_parameter_ids,
        "parameter_set_digests_recorded_when_expected": not missing_parameter_digests,
        "expected_parameter_set_digests_match": not mismatched_parameter_digests,
        "parameter_set_paths_recorded_when_expected": not missing_parameter_paths,
        "expected_parameter_set_paths_match": not mismatched_parameter_paths,
        "objective_observable_ids_recorded_when_expected": not missing_objective_ids,
        "expected_objective_observable_ids_match": not mismatched_objective_ids,
        "objective_observable_families_recorded_when_expected": not missing_objective_families,
        "expected_objective_observable_families_match": not mismatched_objective_families,
    }
    issues = [name for name, ok in checks.items() if not ok]
    return {
        "policy": "build123d_volume_crosscheck_source_identity_gate",
        "status": "ok" if not issues else "needs_attention",
        "sources": sorted(by_source),
        "expected_measurement_methods": expected_methods,
        "expected_body_identity_keys": expected_body_keys,
        "expected_source_artifact_ids": expected_artifacts,
        "expected_parameter_set_artifact_ids": expected_parameter_ids,
        "expected_parameter_set_digests": expected_parameter_digests,
        "expected_parameter_set_paths": expected_parameter_paths,
        "expected_objective_observable_ids": expected_objective_ids,
        "expected_objective_observable_families": expected_objective_families,
        "missing_expected_sources": missing_expected_sources,
        "missing_measurement_methods": missing_methods,
        "mismatched_measurement_methods": mismatched_methods,
        "missing_body_identity_keys": missing_body_keys,
        "mismatched_body_identity_keys": mismatched_body_keys,
        "missing_source_artifact_ids": missing_artifacts,
        "mismatched_source_artifact_ids": mismatched_artifacts,
        "missing_parameter_set_artifact_ids": missing_parameter_ids,
        "mismatched_parameter_set_artifact_ids": mismatched_parameter_ids,
        "missing_parameter_set_digests": missing_parameter_digests,
        "mismatched_parameter_set_digests": mismatched_parameter_digests,
        "missing_parameter_set_paths": missing_parameter_paths,
        "mismatched_parameter_set_paths": mismatched_parameter_paths,
        "missing_objective_observable_ids": missing_objective_ids,
        "mismatched_objective_observable_ids": mismatched_objective_ids,
        "missing_objective_observable_families": missing_objective_families,
        "mismatched_objective_observable_families": mismatched_objective_families,
        "checks": checks,
        "issues": issues,
        "notes": [
            "Run this after volume crosscheck when Cubit, CST, or another CAD row set supplies source metadata.",
            "A volume match is not reusable if the source artifact id, measurement method, or body identity key drifted.",
            "For parametric CAD or optimization loops, bind the parameter-set artifact and objective observable identity before reusing volume rows.",
        ],
    }


def shape_cad_route_source_contract_gate(
    shape_rows,
    external_crosscheck_summary,
    *,
    required_source_groups=(("cubit", "coreform_cubit", "coreform"), ("external_cad", "cst_import")),
    expected_route="cubit_hex_or_mixed_path",
    expected_authoring_sources=("build123d", "build123d_occt", "occt"),
    disallowed_routes=("netgen_tri_tet_path", "tet_only"),
    expected_length_unit=None,
    expected_area_unit=None,
    expected_volume_unit=None,
    required_metadata_fields=(),
):
    """Check build123d CAD rows before Cubit/CAD cross-validation handoff.

    This gate ties together three pieces of evidence that otherwise drift apart
    in loop artifacts: build123d/OCCT authored shape rows, an external CAD
    volume or mass-property crosscheck, and the downstream Cubit hex/mixed mesh
    route.  It is deliberately a policy gate; the numerical tolerance still
    lives in ``shape_volume_crosscheck_summary`` or
    ``shape_mass_property_crosscheck_summary``.
    """

    rows = [dict(row) for row in shape_rows]
    if not rows:
        raise ValueError("shape_rows must not be empty")
    summary = dict(external_crosscheck_summary or {})
    route = str(expected_route or "").strip()
    accepted_authors = {
        str(source).strip().lower()
        for source in expected_authoring_sources
        if str(source).strip()
    }
    blocked_routes = {
        str(item).strip().lower()
        for item in disallowed_routes
        if str(item).strip()
    }

    def has_bbox(row):
        bbox = row.get("bounding_box") or row.get("bbox")
        return isinstance(bbox, dict) and bbox.get("size") is not None

    def _norm(value):
        return str(value or "").strip()

    def _norm_lower(value):
        return _norm(value).lower()

    def _norm_unit(value):
        return _norm_lower(value).replace(" ", "")

    def _row_unit(row, kind):
        key_groups = {
            "length": ("length_unit", "unit_length", "lengthUnit"),
            "area": ("area_unit", "unit_area", "areaUnit"),
            "volume": ("volume_unit", "unit_volume", "volumeUnit"),
        }
        units = row.get("units")
        for key in key_groups[kind]:
            value = row.get(key)
            if value:
                return _norm_unit(value)
        if isinstance(units, dict):
            for key in key_groups[kind]:
                value = units.get(key) or units.get(kind)
                if value:
                    return _norm_unit(value)
        return ""

    names = [_norm(row.get("name")) for row in rows]
    geometry_ids = [_norm(row.get("geometry_id")) for row in rows]
    authoring_sources = [
        _norm(row.get("authoring_source") or row.get("source_kind") or row.get("source"))
        for row in rows
    ]
    routes = [
        _norm(row.get("mesh_route") or row.get("routing_hint"))
        for row in rows
    ]
    length_units = sorted({unit for unit in (_row_unit(row, "length") for row in rows) if unit})
    area_units = sorted({unit for unit in (_row_unit(row, "area") for row in rows) if unit})
    volume_units = sorted({unit for unit in (_row_unit(row, "volume") for row in rows) if unit})
    expected_length_unit_norm = _norm_unit(expected_length_unit)
    expected_area_unit_norm = _norm_unit(expected_area_unit)
    expected_volume_unit_norm = _norm_unit(expected_volume_unit)
    route_set = {item for item in routes if item}
    summary_policy = summary.get("policy")
    summary_sources = [_norm(source) for source in summary.get("sources", [])]
    summary_source_lut = {_norm_lower(source): source for source in summary_sources}
    required_groups = [
        tuple(_norm(source) for source in group if _norm(source))
        for group in required_source_groups
    ]
    required_groups = [group for group in required_groups if group]
    required_metadata = [
        _norm(field)
        for field in required_metadata_fields
        if _norm(field)
    ]
    missing_metadata_by_shape = {
        name or f"shape_{index}": [
            field for field in required_metadata
            if not _norm(row.get(field))
        ]
        for index, (name, row) in enumerate(zip(names, rows), start=1)
    }
    missing_metadata_by_shape = {
        name: missing
        for name, missing in missing_metadata_by_shape.items()
        if missing
    }
    matched_groups = []
    missing_groups = []
    for group in required_groups:
        match = next(
            (summary_source_lut[_norm_lower(source)] for source in group if _norm_lower(source) in summary_source_lut),
            None,
        )
        if match is None:
            missing_groups.append(list(group))
        else:
            matched_groups.append({"accepted": list(group), "matched": match})

    summary_known = summary_policy in {
        "build123d_external_cad_volume_crosscheck",
        "build123d_external_cad_volume_area_bbox_crosscheck",
    }
    summary_ok = summary.get("status") == "ok" and (
        summary.get("ok_for_cad_roundtrip_volume") is True
        or summary.get("ok_for_cad_roundtrip_mass_properties") is True
    )
    route_lower = {_norm_lower(item) for item in route_set}
    checks = {
        "shape_rows_present": bool(rows),
        "shape_names_recorded": all(bool(name) for name in names),
        "shape_geometry_ids_recorded": all(bool(value) for value in geometry_ids),
        "shape_geometry_ids_unique": len(set(geometry_ids)) == len(geometry_ids),
        "shape_rows_have_volume_area_bbox": all(
            row.get("volume") is not None and row.get("area") is not None and has_bbox(row)
            for row in rows
        ),
        "authoring_source_recorded": all(bool(source) for source in authoring_sources),
        "authoring_source_is_build123d": all(
            _norm_lower(source) in accepted_authors for source in authoring_sources
        ),
        "mesh_route_recorded": all(bool(item) for item in routes),
        "mesh_route_matches_expected": route_set == {route},
        "disallowed_routes_absent": not (route_lower & blocked_routes),
        "external_crosscheck_policy_known": summary_known,
        "external_crosscheck_ok": summary_ok,
        "external_crosscheck_sources_recorded": bool(summary_sources),
        "required_source_groups_present": not missing_groups,
    }
    if required_metadata:
        checks["required_shape_metadata_present"] = not missing_metadata_by_shape
    if length_units or area_units or volume_units:
        checks["shape_unit_metadata_unique_when_present"] = (
            len(length_units) <= 1 and len(area_units) <= 1 and len(volume_units) <= 1
        )
    if expected_length_unit_norm:
        checks["shape_length_unit_expected_ok"] = length_units == [expected_length_unit_norm]
    if expected_area_unit_norm:
        checks["shape_area_unit_expected_ok"] = area_units == [expected_area_unit_norm]
    if expected_volume_unit_norm:
        checks["shape_volume_unit_expected_ok"] = volume_units == [expected_volume_unit_norm]
    issues = []
    if not checks["shape_rows_have_volume_area_bbox"]:
        issues.append("shape rows must carry volume, area, and bounding_box")
    if not checks["authoring_source_is_build123d"]:
        issues.append("shape rows must record a build123d/OCCT authoring source")
    if not checks["mesh_route_matches_expected"]:
        issues.append("shape rows must route to the expected Cubit hex/mixed mesh lane")
    if not checks["external_crosscheck_ok"]:
        issues.append("external CAD volume or mass-property crosscheck is missing or not ok")
    if missing_groups:
        issues.append("external crosscheck is missing one or more required source groups")
    if missing_metadata_by_shape:
        issues.append("shape rows are missing required CAD provenance metadata")
    if (
        expected_length_unit_norm
        and checks.get("shape_length_unit_expected_ok") is False
        or expected_area_unit_norm
        and checks.get("shape_area_unit_expected_ok") is False
        or expected_volume_unit_norm
        and checks.get("shape_volume_unit_expected_ok") is False
    ):
        issues.append("shape rows must carry the expected CAD units before cross-source volume reuse")

    return {
        "policy": "build123d_cad_route_source_contract_gate",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "shape_names": names,
        "geometry_ids": sorted(set(geometry_ids)),
        "authoring_sources": sorted(set(source for source in authoring_sources if source)),
        "routes": sorted(route_set),
        "units": {"length": length_units, "area": area_units, "volume": volume_units},
        "expected_units": {
            "length": expected_length_unit_norm or None,
            "area": expected_area_unit_norm or None,
            "volume": expected_volume_unit_norm or None,
        },
        "expected_route": route,
        "external_crosscheck_policy": summary_policy,
        "external_crosscheck_sources": summary_sources,
        "required_source_groups": [list(group) for group in required_groups],
        "matched_source_groups": matched_groups,
        "missing_source_groups": missing_groups,
        "required_metadata_fields": required_metadata,
        "missing_required_metadata_by_shape": missing_metadata_by_shape,
        "checks": checks,
        "issues": issues,
        "notes": [
            "Use after volume or mass-property crosscheck, before Cubit/CST/solver-ready promotion.",
            "build123d owns CAD intent and authoring identity; Cubit owns hex/mixed mesh evidence.",
            "Tet-only Netgen routes remain useful, but they should not be mixed into this Cubit route gate.",
            "Record length/area/volume units when this row will be compared with Cubit, CST, or external CAD volume evidence.",
            "When solver-ready promotion needs provenance, require CAD row fields such as recipe_id, cad_kernel, cad_kernel_version, script_path, and export_id.",
        ],
    }


def shape_external_cad_volume_evidence_package_gate(
    shape_rows,
    volume_summary,
    *,
    required_sources=("cubit", "cst_import"),
    max_allowed_volume_rel_error=None,
    expected_measurement_methods=None,
    expected_body_identity_keys=None,
    expected_source_artifact_ids=None,
    expected_parameter_set_artifact_ids=None,
    expected_parameter_set_digests=None,
    expected_parameter_set_paths=None,
    expected_objective_observable_ids=None,
    expected_objective_observable_families=None,
    expected_route="cubit_hex_or_mixed_path",
    expected_length_unit=None,
    expected_area_unit=None,
    expected_volume_unit=None,
    required_metadata_fields=(),
):
    """Bundle external CAD volume evidence before Cubit/CST reuse.

    This is the build123d-side package gate for slots where a geometry claim is
    stronger than "OCCT volume matched a number": the same build123d shape rows
    must have a valid multi-source volume summary, explicit source identity for
    each CAD exporter/importer, and a route contract for the Cubit hex/mixed
    lane.  It composes the narrower gates instead of replacing them.
    """

    rows = [dict(row) for row in shape_rows]
    summary = dict(volume_summary or {})
    required = tuple(str(source).strip() for source in required_sources if str(source).strip())
    coverage_gate = shape_volume_crosscheck_source_coverage_gate(
        summary,
        required_sources=required,
        max_allowed_volume_rel_error=max_allowed_volume_rel_error,
    )
    identity_gate = shape_volume_crosscheck_source_identity_gate(
        summary,
        expected_measurement_methods=expected_measurement_methods,
        expected_body_identity_keys=expected_body_identity_keys,
        expected_source_artifact_ids=expected_source_artifact_ids,
        expected_parameter_set_artifact_ids=expected_parameter_set_artifact_ids,
        expected_parameter_set_digests=expected_parameter_set_digests,
        expected_parameter_set_paths=expected_parameter_set_paths,
        expected_objective_observable_ids=expected_objective_observable_ids,
        expected_objective_observable_families=expected_objective_observable_families,
    )
    route_gate = shape_cad_route_source_contract_gate(
        rows,
        summary,
        required_source_groups=tuple((source,) for source in required),
        expected_route=expected_route,
        expected_length_unit=expected_length_unit,
        expected_area_unit=expected_area_unit,
        expected_volume_unit=expected_volume_unit,
        required_metadata_fields=required_metadata_fields,
    )

    comparison_sets = [dict(item) for item in summary.get("comparison_sets", []) or []]
    summary_sources = sorted(str(source).strip() for source in summary.get("sources", []) if str(source).strip())
    coverage_sources = sorted(str(source).strip() for source in coverage_gate.get("sources", []) if str(source).strip())
    identity_sources = sorted(str(source).strip() for source in identity_gate.get("sources", []) if str(source).strip())
    route_sources = sorted(
        str(source).strip()
        for source in route_gate.get("external_crosscheck_sources", [])
        if str(source).strip()
    )

    def _has(item, key):
        return bool(str(item.get(key, "") or "").strip())

    artifact_ids = [
        str(item.get("source_artifact_id", "") or "").strip()
        for item in comparison_sets
        if str(item.get("source_artifact_id", "") or "").strip()
    ]
    checks = {
        "shape_rows_present": bool(rows),
        "volume_summary_policy_known": summary.get("policy") == "build123d_external_cad_volume_crosscheck",
        "volume_summary_ok": summary.get("status") == "ok"
        and summary.get("ok_for_cad_roundtrip_volume") is True,
        "coverage_gate_ok": coverage_gate.get("status") == "ok",
        "source_identity_gate_ok": identity_gate.get("status") == "ok",
        "route_contract_gate_ok": route_gate.get("status") == "ok",
        "reference_row_count_matches_shape_rows": int(summary.get("n_reference_rows", -1)) == len(rows),
        "source_identity_metadata_complete": bool(comparison_sets)
        and all(
            _has(item, "measurement_method")
            and _has(item, "body_identity_key")
            and _has(item, "source_artifact_id")
            for item in comparison_sets
        ),
        "source_artifact_ids_unique": len(set(artifact_ids)) == len(artifact_ids),
        "gate_source_sets_consistent": summary_sources == coverage_sources == identity_sources
        and set(route_sources) == set(summary_sources),
    }
    issues = [name for name, ok in checks.items() if not ok]
    return {
        "policy": "build123d_external_cad_volume_evidence_package_gate",
        "status": "ok" if not issues else "needs_attention",
        "required_sources": list(required),
        "sources": summary_sources,
        "coverage_gate_status": coverage_gate.get("status"),
        "source_identity_gate_status": identity_gate.get("status"),
        "route_contract_gate_status": route_gate.get("status"),
        "max_volume_rel_error": summary.get("max_volume_rel_error"),
        "source_artifact_ids": artifact_ids,
        "coverage_gate": coverage_gate,
        "source_identity_gate": identity_gate,
        "route_contract_gate": route_gate,
        "checks": checks,
        "issues": issues,
        "notes": [
            "Use after shape_volume_crosscheck_summary when Cubit and external CAD rows are both part of the claim.",
            "This package gate keeps volume values, source identity, body identity, and Cubit hex/mixed route intent together.",
            "It is CAD evidence only; Cubit still owns mesh quality, order, and .vol inventory evidence.",
        ],
    }


def cst_cad_volume_export_manifest_gate(
    manifest,
    *,
    expected_geometry_id=None,
    expected_export_id=None,
    required_shape_names=(),
):
    """Check CST CAD volume rows before using them in build123d/Cubit crosschecks.

    CST CAD/solid volume exports are useful as an external-CAD source for
    ``shape_volume_crosscheck_summary``, but only after source, identity, unit,
    and body names are explicit.  This gate returns normalized metre-cubed rows
    that can be passed on as a ``cst_import`` measured set.
    """

    data = dict(manifest)
    rows = list(data.get("volume_rows") or data.get("rows") or data.get("cad_volume_rows") or [])
    if not rows:
        raise ValueError("manifest must include volume_rows/rows")

    def _norm(value):
        return str(value or "").strip()

    def _status_ok(value):
        return str(value or "ok").strip().lower().replace("-", "_") in {
            "ok",
            "pass",
            "passed",
            "verified",
            "exported",
        }

    unit_scales = {
        "m3": 1.0,
        "m^3": 1.0,
        "meter^3": 1.0,
        "metre^3": 1.0,
        "mm3": 1.0e-9,
        "mm^3": 1.0e-9,
        "cm3": 1.0e-6,
        "cm^3": 1.0e-6,
    }

    manifest_unit = _norm(data.get("volume_unit") or data.get("unit") or "m^3").lower().replace(" ", "")
    normalized_rows = []
    unknown_unit_rows = []
    missing_name_rows = []
    nonpositive_rows = []
    bad_status_rows = []
    for index, row in enumerate(rows, start=1):
        name = _norm(row.get("name") or row.get("solid_name") or row.get("body_name"))
        if not name:
            missing_name_rows.append(index)
        row_unit = _norm(row.get("volume_unit") or row.get("unit") or manifest_unit).lower().replace(" ", "")
        if row_unit not in unit_scales:
            unknown_unit_rows.append({"row": index, "unit": row_unit})
            scale = math.nan
        else:
            scale = unit_scales[row_unit]
        if "volume_m3" in row:
            volume_m3 = float(row["volume_m3"])
        else:
            volume = row.get("volume", row.get("volume_value"))
            if volume is None:
                volume = row.get("cad_volume")
            volume_m3 = float(volume) * scale if volume is not None and math.isfinite(scale) else math.nan
        if not math.isfinite(volume_m3) or volume_m3 <= 0.0:
            nonpositive_rows.append(index)
        if not _status_ok(row.get("status", data.get("status", "ok"))):
            bad_status_rows.append(index)
        normalized_rows.append(
            {
                "name": name,
                "volume": volume_m3,
                "volume_m3": volume_m3,
                "source_unit": row_unit,
                "source_tool": "CST",
                "row": index,
            }
        )

    names = [row["name"] for row in normalized_rows if row["name"]]
    duplicate_names = sorted(name for name, count in Counter(names).items() if count > 1)
    required = [_norm(name) for name in required_shape_names if _norm(name)]
    missing_required = [name for name in required if name not in set(names)]
    source_key = str(data.get("source_tool", "")).strip().lower()
    geometry_id = _norm(data.get("geometry_id"))
    export_id = _norm(data.get("export_id"))
    project_id = _norm(data.get("project_id"))
    run_id = _norm(data.get("run_id"))
    checks = {
        "source_tool_is_cst": "cst" in source_key,
        "manifest_identity_present": bool(project_id and run_id and export_id and geometry_id),
        "expected_geometry_id_matches": expected_geometry_id is None or geometry_id == str(expected_geometry_id),
        "expected_export_id_matches": expected_export_id is None or export_id == str(expected_export_id),
        "volume_rows_present": bool(rows),
        "volume_units_known": not unknown_unit_rows,
        "shape_names_present": not missing_name_rows,
        "shape_names_unique": not duplicate_names,
        "required_shape_names_present": not missing_required,
        "positive_volumes": not nonpositive_rows,
        "row_status_ok": not bad_status_rows,
    }
    issues = [name for name, ok in checks.items() if not ok]
    return {
        "policy": "cst_cad_volume_export_manifest_gate",
        "status": "ok" if not issues else "needs_attention",
        "project_id": project_id,
        "run_id": run_id,
        "export_id": export_id,
        "geometry_id": geometry_id,
        "source_tool": data.get("source_tool"),
        "volume_unit": manifest_unit,
        "required_shape_names": required,
        "shape_names": names,
        "duplicate_shape_names": duplicate_names,
        "missing_required_shape_names": missing_required,
        "unknown_unit_rows": unknown_unit_rows,
        "missing_name_rows": missing_name_rows,
        "nonpositive_rows": nonpositive_rows,
        "bad_status_rows": bad_status_rows,
        "normalized_rows": normalized_rows,
        "checks": checks,
        "issues": issues,
        "notes": [
            "Use normalized_rows as the CST measured set for shape_volume_crosscheck_summary.",
            "Volume is a CAD handoff gate; solver-ready status still needs mesh, ports, or physics-specific manifests.",
        ],
    }


def _normalize_shape_measurement_sets(measured_sets):
    def _rows_payload(value):
        if isinstance(value, dict):
            rows = value.get("rows", value.get("measurements", value))
            if isinstance(rows, dict) and ("volume" in rows or "area" in rows):
                return [rows]
            return rows
        return value

    if isinstance(measured_sets, dict):
        if "rows" in measured_sets or "measurements" in measured_sets or "volume" in measured_sets:
            return [("external_cad", _rows_payload(measured_sets))]
        return [
            (str(label), _rows_payload(rows))
            for label, rows in measured_sets.items()
        ]
    normalized = []
    for index, item in enumerate(list(measured_sets), start=1):
        if isinstance(item, dict) and ("rows" in item or "measurements" in item):
            label = str(item.get("source") or item.get("label") or f"external_cad_{index}")
            rows = _rows_payload(item)
        else:
            label = f"external_cad_{index}"
            rows = _rows_payload(item)
        normalized.append((label, rows))
    return normalized


def shape_name_identity_gate(reference_rows, measured_rows, measured_label="measured"):
    """Check that a CAD round trip preserved the same named shape multiset.

    Volume and area comparisons iterate over reference rows, so an imported
    CAD source can otherwise carry an extra solid without being noticed.  Run
    this gate before trusting assembly-level mass-property comparisons.
    """

    reference = list(reference_rows)
    measured = list(measured_rows)
    reference_names = [str(row.get("name", "")).strip() for row in reference]
    measured_names = [str(row.get("name", "")).strip() for row in measured]
    reference_counts = Counter(reference_names)
    measured_counts = Counter(measured_names)
    reference_missing_name_count = reference_counts.pop("", 0)
    measured_missing_name_count = measured_counts.pop("", 0)
    duplicate_reference_names = sorted(
        name for name, count in reference_counts.items() if count > 1
    )
    duplicate_measured_names = sorted(
        name for name, count in measured_counts.items() if count > 1
    )
    missing_names = sorted(
        name for name, count in reference_counts.items()
        if measured_counts.get(name, 0) < count
    )
    extra_names = sorted(
        name for name, count in measured_counts.items()
        if reference_counts.get(name, 0) < count
    )
    count_mismatches = [
        {
            "name": name,
            "reference_count": reference_counts.get(name, 0),
            "measured_count": measured_counts.get(name, 0),
        }
        for name in sorted(set(reference_counts) | set(measured_counts))
        if reference_counts.get(name, 0) != measured_counts.get(name, 0)
    ]
    checks = {
        "reference_names_present": reference_missing_name_count == 0,
        "measured_names_present": measured_missing_name_count == 0,
        "reference_names_unique": not duplicate_reference_names,
        "measured_names_unique": not duplicate_measured_names,
        "same_name_multiset": not count_mismatches,
    }
    return {
        "policy": "build123d_cad_roundtrip_named_shape_identity_gate",
        "measured_label": str(measured_label),
        "n_reference_rows": len(reference),
        "n_measured_rows": len(measured),
        "reference_names": sorted(reference_counts),
        "measured_names": sorted(measured_counts),
        "reference_missing_name_count": reference_missing_name_count,
        "measured_missing_name_count": measured_missing_name_count,
        "duplicate_reference_names": duplicate_reference_names,
        "duplicate_measured_names": duplicate_measured_names,
        "missing_names": missing_names,
        "extra_names": extra_names,
        "count_mismatches": count_mismatches,
        "checks": checks,
        "status": "ok" if all(checks.values()) else "needs_attention",
        "version_note": (
            "Use this before volume/area/bbox gates for multi-body build123d "
            "STEP round trips through Cubit, CST, or another CAD kernel."
        ),
    }


def shape_role_metadata_gate(
    rows,
    required_names=(),
    required_roles=(),
    required_materials=(),
    source_label="build123d",
):
    """Check that named CAD bodies carry solver-handoff role/material metadata.

    Geometry gates can prove that a STEP round trip preserved bodies, volumes
    and areas, but a solver also needs to know what each body *means*.  Keep a
    tiny row contract next to build123d assemblies before exporting them:
    ``{"name": "core", "role": "magnetic_core", "material": "steel"}``.
    """

    normalized = []
    for row in list(rows):
        item = dict(row)
        name = str(item.get("name", "")).strip()
        role = str(
            item.get("role")
            or item.get("solver_role")
            or item.get("region_role")
            or item.get("boundary_role")
            or ""
        ).strip()
        material = str(
            item.get("material")
            or item.get("material_name")
            or item.get("mat")
            or ""
        ).strip()
        normalized.append(
            {
                "name": name,
                "role": role,
                "material": material,
                "source_row": item,
            }
        )

    name_counts = Counter(row["name"] for row in normalized)
    missing_name_count = name_counts.pop("", 0)
    duplicate_names = sorted(name for name, count in name_counts.items() if count > 1)
    names = set(name_counts)
    roles = {row["role"] for row in normalized if row["role"]}
    materials = {row["material"] for row in normalized if row["material"]}
    rows_missing_role = [row["name"] for row in normalized if row["name"] and not row["role"]]
    rows_missing_material = [
        row["name"] for row in normalized if row["name"] and not row["material"]
    ]
    required_name_set = {str(name).strip() for name in required_names if str(name).strip()}
    required_role_set = {str(role).strip() for role in required_roles if str(role).strip()}
    required_material_set = {
        str(material).strip()
        for material in required_materials
        if str(material).strip()
    }
    missing_required_names = sorted(required_name_set - names)
    missing_required_roles = sorted(required_role_set - roles)
    missing_required_materials = sorted(required_material_set - materials)
    checks = {
        "names_present": missing_name_count == 0,
        "names_unique": not duplicate_names,
        "all_rows_have_role": not rows_missing_role,
        "all_rows_have_material": not rows_missing_material,
        "required_names_present": not missing_required_names,
        "required_roles_present": not missing_required_roles,
        "required_materials_present": not missing_required_materials,
    }
    return {
        "policy": "build123d_solver_handoff_role_material_metadata_gate",
        "source_label": str(source_label),
        "n_rows": len(normalized),
        "names": sorted(names),
        "roles": sorted(roles),
        "materials": sorted(materials),
        "missing_name_count": missing_name_count,
        "duplicate_names": duplicate_names,
        "rows_missing_role": sorted(rows_missing_role),
        "rows_missing_material": sorted(rows_missing_material),
        "missing_required_names": missing_required_names,
        "missing_required_roles": missing_required_roles,
        "missing_required_materials": missing_required_materials,
        "checks": checks,
        "status": "ok" if all(checks.values()) else "needs_attention",
        "version_note": (
            "Run this after shape_name_identity_gate and before meshing so "
            "build123d bodies retain solver role/material intent across STEP, "
            "Cubit, CST, or Netgen handoff."
        ),
    }


def shape_transition_role_metadata_gate(
    rows,
    *,
    required_roles=("hex_region", "mesh_transition", "tet_region"),
    transition_role="mesh_transition",
    required_transition_kind="pyramid",
    required_connected_roles=("hex_region", "tet_region"),
    required_surface_kinds=(),
    required_interface_roles=(),
    expected_downstream_material_names=(),
    allowed_zero_downstream_material_names=(),
    require_positive_volume=True,
    source_label="build123d",
):
    """Check CAD-side metadata for a future hex-to-tet transition handoff.

    build123d does not create Cubit pyramid elements, but it can preserve the
    solver intent before STEP/Cubit handoff: which body is the hex-led region,
    which body is the tet region, which body is the transition envelope, and
    which downstream Cubit material/block labels should appear in sidecars.
    """

    def as_list(value):
        if value is None:
            return []
        if isinstance(value, str):
            return [part.strip() for part in value.replace(";", ",").split(",") if part.strip()]
        return [str(item).strip() for item in value if str(item).strip()]

    normalized = []
    for row in list(rows):
        item = dict(row)
        name = str(item.get("name", "")).strip()
        role = str(
            item.get("role")
            or item.get("solver_role")
            or item.get("region_role")
            or ""
        ).strip()
        material = str(
            item.get("material")
            or item.get("material_name")
            or item.get("mat")
            or ""
        ).strip()
        transition_kind = str(
            item.get("transition_kind")
            or item.get("mesh_transition_kind")
            or ""
        ).strip()
        connected_roles = as_list(
            item.get("connects_roles")
            or item.get("connected_roles")
            or item.get("transition_between_roles")
            or item.get("transition_between")
        )
        surface_kinds = as_list(
            item.get("surface_kinds")
            or item.get("expected_surface_kinds")
            or item.get("mesh_surface_kinds")
            or item.get("surface_family_intent")
        )
        interface_roles = as_list(
            item.get("interface_roles")
            or item.get("expected_interface_roles")
            or item.get("transition_interface_roles")
            or item.get("mesh_interface_roles")
            or item.get("interface_role_intent")
        )
        downstream_material_name = str(
            item.get("downstream_material_name")
            or item.get("cubit_material_name")
            or item.get("material_block_name")
            or item.get("sidecar_material_name")
            or item.get("expected_material_name")
            or ""
        ).strip()
        volume = item.get("volume")
        normalized.append({
            "name": name,
            "role": role,
            "material": material,
            "transition_kind": transition_kind,
            "connected_roles": connected_roles,
            "surface_kinds": surface_kinds,
            "interface_roles": interface_roles,
            "downstream_material_name": downstream_material_name,
            "volume": None if volume is None else float(volume),
            "source_row": item,
        })

    names = [row["name"] for row in normalized if row["name"]]
    roles = {row["role"] for row in normalized if row["role"]}
    required_role_set = {str(role).strip() for role in required_roles if str(role).strip()}
    connected_role_set = {str(role).strip() for role in required_connected_roles if str(role).strip()}
    required_surface_kind_set = {
        str(kind).strip() for kind in required_surface_kinds if str(kind).strip()
    }
    required_interface_role_set = {
        str(role).strip()
        for role in required_interface_roles
        if str(role).strip()
    }
    expected_downstream_material_set = {
        str(name).strip()
        for name in expected_downstream_material_names
        if str(name).strip()
    }
    allowed_zero_downstream_material_set = {
        str(name).strip()
        for name in allowed_zero_downstream_material_names
        if str(name).strip()
    }
    transition_rows = [row for row in normalized if row["role"] == transition_role]
    transition_kinds = {row["transition_kind"] for row in transition_rows if row["transition_kind"]}
    connected_roles_union = {
        role for row in transition_rows for role in row["connected_roles"] if role
    }
    transition_surface_kinds = {
        kind for row in transition_rows for kind in row["surface_kinds"] if kind
    }
    transition_interface_roles = {
        role for row in transition_rows for role in row["interface_roles"] if role
    }
    downstream_material_names = {
        row["downstream_material_name"] for row in normalized if row["downstream_material_name"]
    }
    rows_missing_material = [
        row["name"] for row in normalized if row["name"] and not row["material"]
    ]
    rows_missing_downstream_material_name = [
        row["name"] for row in normalized if row["name"] and not row["downstream_material_name"]
    ]
    rows_missing_interface_roles = [
        row["name"]
        for row in transition_rows
        if row["name"] and required_interface_role_set and not row["interface_roles"]
    ]
    rows_missing_volume = [
        row["name"] for row in normalized if row["name"] and row["volume"] is None
    ]
    rows_nonpositive_volume = [
        row["name"]
        for row in normalized
        if row["name"] and row["volume"] is not None and row["volume"] <= 0.0
    ]
    missing_required_roles = sorted(required_role_set - roles)
    missing_expected_downstream_material_names = sorted(
        expected_downstream_material_set - downstream_material_names
    )
    missing_required_interface_roles = sorted(
        required_interface_role_set - transition_interface_roles
    )
    missing_allowed_zero_downstream_material_names = sorted(
        allowed_zero_downstream_material_set - downstream_material_names
    )
    duplicate_names = sorted(
        name for name, count in Counter(names).items() if count > 1
    )
    checks = {
        "names_present": len(names) == len(normalized),
        "names_unique": not duplicate_names,
        "required_roles_present": not missing_required_roles,
        "transition_row_present": bool(transition_rows),
        "transition_kind_recorded": bool(transition_kinds),
        "transition_kind_matches": transition_kinds == {str(required_transition_kind)},
        "transition_connects_required_roles": connected_role_set.issubset(connected_roles_union),
        "required_surface_kinds_present": (
            not required_surface_kind_set
            or required_surface_kind_set.issubset(transition_surface_kinds)
        ),
        "interface_roles_recorded": (
            not required_interface_role_set
            or bool(transition_interface_roles)
        ),
        "required_interface_roles_present": not missing_required_interface_roles,
        "downstream_material_names_recorded": (
            not expected_downstream_material_set
            or not rows_missing_downstream_material_name
        ),
        "expected_downstream_material_names_present": not missing_expected_downstream_material_names,
        "allowed_zero_downstream_material_names_declared": (
            not allowed_zero_downstream_material_set
            or not missing_allowed_zero_downstream_material_names
        ),
        "all_rows_have_material": not rows_missing_material,
        "all_rows_have_volume": not rows_missing_volume,
        "all_rows_have_positive_volume": not rows_nonpositive_volume if require_positive_volume else True,
    }
    return {
        "policy": "build123d_hex_tet_transition_role_metadata_gate",
        "source_label": str(source_label),
        "n_rows": len(normalized),
        "names": sorted(names),
        "duplicate_names": duplicate_names,
        "roles": sorted(roles),
        "required_roles": sorted(required_role_set),
        "missing_required_roles": missing_required_roles,
        "transition_role": str(transition_role),
        "transition_kinds": sorted(transition_kinds),
        "required_transition_kind": str(required_transition_kind),
        "connected_roles": sorted(connected_roles_union),
        "required_connected_roles": sorted(connected_role_set),
        "surface_kinds": sorted(transition_surface_kinds),
        "required_surface_kinds": sorted(required_surface_kind_set),
        "interface_roles": sorted(transition_interface_roles),
        "required_interface_roles": sorted(required_interface_role_set),
        "missing_required_interface_roles": missing_required_interface_roles,
        "downstream_material_names": sorted(downstream_material_names),
        "expected_downstream_material_names": sorted(expected_downstream_material_set),
        "allowed_zero_downstream_material_names": sorted(allowed_zero_downstream_material_set),
        "missing_expected_downstream_material_names": missing_expected_downstream_material_names,
        "missing_allowed_zero_downstream_material_names": missing_allowed_zero_downstream_material_names,
        "rows_missing_material": sorted(rows_missing_material),
        "rows_missing_downstream_material_name": sorted(rows_missing_downstream_material_name),
        "rows_missing_interface_roles": sorted(rows_missing_interface_roles),
        "rows_missing_volume": sorted(rows_missing_volume),
        "rows_nonpositive_volume": sorted(rows_nonpositive_volume),
        "checks": checks,
        "status": "ok" if all(checks.values()) else "needs_attention",
        "version_note": (
            "Run this on build123d assembly metadata before handing a future "
            "hex+tet model to Cubit; the pyramid is a mesh transition contract, "
            "not a build123d primitive requirement.  When the downstream Cubit "
            "path will verify mixed surface families, preserve that surface "
            "family intent here before STEP export.  If the downstream Cubit "
            "sidecar will verify material/block labels, record those expected "
            "downstream material names here rather than relying on CAD body "
            "names to match by accident.  If the downstream Cubit package will "
            "verify hex-to-transition and transition-to-tet interface roles, "
            "record those role names in the CAD-side transition manifest too."
        ),
    }


def shape_cubit_meshing_scheme_intent_gate(
    shape_rows,
    *,
    scheme_trace_gate=None,
    required_roles=("hex_region", "mesh_transition", "tet_region"),
    expected_scheme_by_role=None,
    required_command_fragments=("imprint all", "merge all", "export netgen"),
    expected_export_order=None,
    expected_trace_id=None,
    expected_export_output_artifact_id=None,
    expected_export_output_digest=None,
    expected_export_output_path=None,
    require_downstream_export_output_artifact=False,
    expected_route="cubit_hex_or_mixed_path",
    source_label="build123d",
):
    """Check CAD-side Cubit meshing-scheme intent before Cubit owns the mesh.

    build123d can say which CAD roles are intended for mapped hex, tetmesh, or
    transition meshing, but Cubit still owns the actual volume ids and exported
    mesh.  This gate keeps the CAD-side role-to-scheme intent tied to the
    downstream Cubit scheme trace id and, when available, to
    ``cubit_meshing_scheme_trace_gate`` evidence.
    """

    def text(value):
        return str(value or "").strip()

    def as_list(value):
        if value is None:
            return []
        if isinstance(value, str):
            return [part.strip() for part in value.replace(";", ",").split(",") if part.strip()]
        return [str(item).strip() for item in value if str(item).strip()]

    def norm_scheme(value):
        value = text(value).lower()
        aliases = {
            "mapped": "map",
            "mapped_hex": "map",
            "hex_map": "map",
            "tet": "tetmesh",
            "tetrahedral": "tetmesh",
            "tet_mesh": "tetmesh",
        }
        return aliases.get(value, value)

    rows = [dict(row) for row in shape_rows]
    if not rows:
        raise ValueError("shape_rows must not be empty")

    expected_role_scheme = {
        text(role): norm_scheme(scheme)
        for role, scheme in (expected_scheme_by_role or {}).items()
        if text(role) and text(scheme)
    }
    if not expected_role_scheme:
        expected_role_scheme = {
            "hex_region": "map",
            "mesh_transition": "tetmesh",
            "tet_region": "tetmesh",
        }
    required_role_set = {text(role) for role in required_roles if text(role)}
    required_fragments = [fragment.lower() for fragment in as_list(required_command_fragments)]
    expected_order = None if expected_export_order is None else int(expected_export_order)
    expected_trace = None if expected_trace_id is None else text(expected_trace_id)
    route = text(expected_route)

    normalized = []
    command_fragments = set()
    trace_ids = set()
    export_orders = set()
    route_values = set()
    for row in rows:
        name = text(row.get("name"))
        role = text(row.get("role") or row.get("solver_role") or row.get("region_role"))
        scheme = norm_scheme(
            row.get("expected_cubit_scheme")
            or row.get("cubit_scheme")
            or row.get("mesh_scheme")
            or row.get("expected_mesh_scheme")
            or row.get("downstream_mesh_scheme")
        )
        route_value = text(row.get("mesh_route") or row.get("routing_hint"))
        trace_id = text(
            row.get("downstream_meshing_trace_id")
            or row.get("cubit_meshing_trace_id")
            or row.get("scheme_trace_id")
            or row.get("meshing_trace_id")
        )
        export_order = row.get("expected_cubit_export_order", row.get("export_order"))
        if export_order is not None and text(export_order):
            export_orders.add(int(export_order))
        if route_value:
            route_values.add(route_value)
        if trace_id:
            trace_ids.add(trace_id)
        command_fragments.update(fragment.lower() for fragment in as_list(
            row.get("expected_cubit_command_fragments")
            or row.get("required_command_fragments")
            or row.get("command_fragments")
        ))
        normalized.append(
            {
                "name": name,
                "role": role,
                "expected_cubit_scheme": scheme,
                "mesh_route": route_value,
                "downstream_meshing_trace_id": trace_id,
                "expected_cubit_export_order": export_order,
            }
        )

    names = [row["name"] for row in normalized if row["name"]]
    roles = {row["role"] for row in normalized if row["role"]}
    role_scheme = {
        row["role"]: row["expected_cubit_scheme"]
        for row in normalized
        if row["role"] and row["expected_cubit_scheme"]
    }
    missing_required_roles = sorted(required_role_set - roles)
    missing_expected_scheme_roles = sorted(
        role for role in expected_role_scheme if role_scheme.get(role) != expected_role_scheme[role]
    )
    missing_fragments = sorted(set(required_fragments) - command_fragments)
    duplicate_names = sorted(name for name, count in Counter(names).items() if count > 1)

    downstream = dict(scheme_trace_gate or {})
    downstream_checks = downstream.get("checks", {})
    if not isinstance(downstream_checks, dict):
        downstream_checks = {}
    downstream_commands = "\n".join(
        str(command).lower()
        for command in downstream.get("commands", [])
    )
    downstream_trace_id = text(downstream.get("trace_id"))
    downstream_export_order = downstream.get("export_order")
    downstream_export_order_int = (
        None if downstream_export_order is None or text(downstream_export_order) == "" else int(downstream_export_order)
    )
    downstream_output_artifact_id = text(downstream.get("export_output_artifact_id"))
    downstream_output_digest = text(downstream.get("export_output_digest"))
    downstream_output_path = text(downstream.get("export_output_path"))
    expected_output_artifact = (
        None if expected_export_output_artifact_id is None else text(expected_export_output_artifact_id)
    )
    expected_output_digest = (
        None if expected_export_output_digest is None else text(expected_export_output_digest)
    )
    expected_output_path = (
        None if expected_export_output_path is None else text(expected_export_output_path)
    )
    output_required = bool(
        require_downstream_export_output_artifact
        or expected_output_artifact is not None
        or expected_output_digest is not None
        or expected_output_path is not None
    )

    checks = {
        "shape_rows_present": bool(normalized),
        "shape_names_recorded": len(names) == len(normalized),
        "shape_names_unique": not duplicate_names,
        "required_roles_present": not missing_required_roles,
        "mesh_route_recorded": all(bool(row["mesh_route"]) for row in normalized),
        "mesh_route_matches_expected": route_values == {route},
        "scheme_intent_recorded": all(bool(row["expected_cubit_scheme"]) for row in normalized),
        "expected_scheme_by_role_matches": not missing_expected_scheme_roles,
        "downstream_trace_id_recorded": bool(trace_ids),
        "expected_trace_id_matches": expected_trace is None or trace_ids == {expected_trace},
        "command_fragments_recorded": bool(command_fragments),
        "required_command_fragments_present": not missing_fragments,
        "export_order_recorded": expected_order is None or export_orders == {expected_order},
        "expected_export_order_matches": expected_order is None or export_orders == {expected_order},
    }
    if downstream:
        checks["downstream_scheme_trace_gate_ok"] = (
            downstream.get("policy") == "cubit_meshing_scheme_trace_gate"
            and downstream.get("status") == "ok"
        )
        checks["downstream_trace_id_matches"] = (
            not trace_ids or downstream_trace_id in trace_ids
        )
        checks["downstream_export_order_matches"] = (
            expected_order is None or downstream_export_order_int == expected_order
        )
        checks["downstream_command_fragments_present"] = all(
            fragment in downstream_commands for fragment in required_fragments
        )
        checks["downstream_export_output_artifact_id_recorded_when_required"] = (
            not output_required or bool(downstream_output_artifact_id)
        )
        checks["downstream_export_output_digest_recorded_when_required"] = (
            not output_required or bool(downstream_output_digest)
        )
        checks["downstream_export_output_path_recorded_when_required"] = (
            not output_required or bool(downstream_output_path)
        )
        checks["downstream_export_output_artifact_id_matches"] = (
            expected_output_artifact is None or downstream_output_artifact_id == expected_output_artifact
        )
        checks["downstream_export_output_digest_matches"] = (
            expected_output_digest is None or downstream_output_digest == expected_output_digest
        )
        checks["downstream_export_output_path_matches"] = (
            expected_output_path is None or downstream_output_path == expected_output_path
        )

    issues = [name for name, ok in checks.items() if not ok]
    return {
        "policy": "build123d_cubit_meshing_scheme_intent_gate",
        "source_label": str(source_label),
        "status": "ok" if not issues else "needs_attention",
        "shape_names": sorted(names),
        "duplicate_names": duplicate_names,
        "roles": sorted(roles),
        "required_roles": sorted(required_role_set),
        "missing_required_roles": missing_required_roles,
        "role_scheme_intent": role_scheme,
        "expected_scheme_by_role": expected_role_scheme,
        "missing_expected_scheme_roles": missing_expected_scheme_roles,
        "mesh_routes": sorted(route_values),
        "expected_route": route,
        "downstream_meshing_trace_ids": sorted(trace_ids),
        "expected_trace_id": expected_trace,
        "command_fragments": sorted(command_fragments),
        "required_command_fragments": sorted(required_fragments),
        "missing_required_command_fragments": missing_fragments,
        "export_orders": sorted(export_orders),
        "expected_export_order": expected_order,
        "downstream_scheme_trace_policy": downstream.get("policy"),
        "downstream_scheme_trace_status": downstream.get("status"),
        "downstream_export_output_artifact_id": downstream_output_artifact_id or None,
        "downstream_export_output_digest": downstream_output_digest or None,
        "downstream_export_output_path": downstream_output_path or None,
        "expected_export_output_artifact_id": expected_output_artifact,
        "expected_export_output_digest": expected_output_digest,
        "expected_export_output_path": expected_output_path,
        "require_downstream_export_output_artifact": output_required,
        "checks": checks,
        "issues": issues,
        "notes": [
            "Use this after build123d transition/material/interface intent and before Cubit meshing.",
            "build123d records role-to-scheme intent; Cubit records actual volume ids and export commands.",
            "When a downstream .vol is already exported, bind its artifact id, digest, and path before promoting the CAD handoff.",
            "Do not treat this as proof that the .vol inventory is solver-ready; run Cubit inventory gates separately.",
        ],
    }


def _assembly_location_boolean_mass_identity(row):
    value = row.get(
        "assembly_location_boolean_operand_revision_mass_property_generation_identity"
    )
    if not isinstance(value, dict):
        return None
    generation = str(value.get("assembly_generation", "")).strip()
    part_ids = [str(item).strip() for item in value.get("part_ids", [])]
    evaluated_part_ids = [
        str(item).strip() for item in value.get("evaluated_part_ids", [])
    ]
    locations = [
        str(item).lower() for item in value.get("part_location_sha256", [])
    ]
    evaluated_locations = [
        str(item).lower()
        for item in value.get("evaluated_part_location_sha256", [])
    ]
    operands = [
        str(item).strip() for item in value.get("boolean_operand_revisions", [])
    ]
    evaluated_operands = [
        str(item).strip()
        for item in value.get("evaluated_boolean_operand_revisions", [])
    ]
    try:
        members = [int(item) for item in value.get("compound_member_ids", [])]
        evaluated_members = [
            int(item) for item in value.get("evaluated_compound_member_ids", [])
        ]
    except (TypeError, ValueError):
        return None

    def valid_digest(digest):
        return len(digest) == 64 and all(
            character in "0123456789abcdef" for character in digest
        )

    density_digest = str(value.get("density_map_sha256", "")).lower()
    mass_digest = str(value.get("mass_property_sha256", "")).lower()
    if (
        not generation
        or any(
            value.get(key) != generation
            for key in (
                "location_assembly_generation",
                "boolean_operand_assembly_generation",
                "compound_membership_assembly_generation",
                "density_map_assembly_generation",
                "mass_property_assembly_generation",
            )
        )
        or not part_ids
        or any(not item for item in part_ids)
        or len(set(part_ids)) != len(part_ids)
        or evaluated_part_ids != part_ids
        or len(locations) != len(part_ids)
        or not all(valid_digest(digest) for digest in locations)
        or evaluated_locations != locations
        or len(operands) < 2
        or any(not item for item in operands)
        or len(set(operands)) != len(operands)
        or evaluated_operands != operands
        or not members
        or any(item <= 0 for item in members)
        or len(set(members)) != len(members)
        or evaluated_members != members
        or not valid_digest(density_digest)
        or value.get("evaluated_density_map_sha256") != density_digest
        or not valid_digest(mass_digest)
        or value.get("evaluated_mass_property_sha256") != mass_digest
    ):
        return None
    return (
        generation,
        tuple(part_ids),
        tuple(locations),
        tuple(operands),
        tuple(members),
        density_digest,
        mass_digest,
    )


def _loft_spline_watertight_volume_identity(row):
    value = row.get(
        "loft_spline_tessellation_watertight_volume_generation_identity"
    )
    if not isinstance(value, dict):
        return None
    generation = str(value.get("shape_generation", "")).strip()
    spline_digest = str(value.get("spline_sha256", "")).lower()
    shell_digest = str(value.get("tessellated_shell_sha256", "")).lower()
    volume_digest = str(value.get("volume_result_sha256", "")).lower()
    try:
        chord = float(value.get("chord_tolerance"))
        evaluated_chord = float(value.get("evaluated_chord_tolerance"))
        angle = float(value.get("angular_tolerance_deg"))
        evaluated_angle = float(value.get("evaluated_angular_tolerance_deg"))
        volume = float(value.get("volume"))
        evaluated_volume = float(value.get("evaluated_volume"))
    except (TypeError, ValueError):
        return None

    def valid_digest(digest):
        return len(digest) == 64 and all(
            character in "0123456789abcdef" for character in digest
        )

    unit = str(value.get("length_unit", "")).strip()
    if (
        not generation
        or any(
            value.get(key) != generation
            for key in (
                "spline_shape_generation",
                "tessellation_shape_generation",
                "watertight_shape_generation",
                "volume_shape_generation",
            )
        )
        or not valid_digest(spline_digest)
        or value.get("evaluated_spline_sha256") != spline_digest
        or not math.isfinite(chord)
        or chord <= 0.0
        or not math.isclose(evaluated_chord, chord, rel_tol=0.0, abs_tol=1.0e-18)
        or not math.isfinite(angle)
        or not 0.0 < angle <= 180.0
        or not math.isclose(evaluated_angle, angle, rel_tol=0.0, abs_tol=1.0e-12)
        or unit not in {"m", "cm", "mm"}
        or value.get("evaluated_length_unit") != unit
        or value.get("watertight") is not True
        or value.get("evaluated_watertight") is not True
        or not valid_digest(shell_digest)
        or value.get("evaluated_tessellated_shell_sha256") != shell_digest
        or not math.isfinite(volume)
        or volume <= 0.0
        or not math.isclose(evaluated_volume, volume, rel_tol=1.0e-12, abs_tol=1.0e-18)
        or not valid_digest(volume_digest)
        or value.get("evaluated_volume_result_sha256") != volume_digest
    ):
        return None
    return generation, spline_digest, chord, angle, unit, shell_digest, volume, volume_digest


def _transformed_assembly_mass_identity(row):
    value = row.get(
        "transformed_assembly_com_inertia_axis_density_unit_generation_identity"
    )
    if not isinstance(value, dict):
        return None
    generation = str(value.get("assembly_generation", "")).strip()
    part_ids = [str(item).strip() for item in value.get("part_ids", [])]
    result_part_ids = [
        str(item).strip() for item in value.get("result_part_ids", [])
    ]
    transforms = [
        str(item).lower()
        for item in value.get("local_to_global_transform_sha256", [])
    ]
    result_transforms = [
        str(item).lower()
        for item in value.get("result_local_to_global_transform_sha256", [])
    ]
    try:
        densities = [float(item) for item in value.get("density_kg_m3", [])]
        result_densities = [
            float(item) for item in value.get("result_density_kg_m3", [])
        ]
        center = [float(item) for item in value.get("center_of_mass_m", [])]
        result_center = [
            float(item) for item in value.get("result_center_of_mass_m", [])
        ]
        inertia = [
            [float(component) for component in matrix_row]
            for matrix_row in value.get("inertia_tensor_kg_m2", [])
        ]
        result_inertia = [
            [float(component) for component in matrix_row]
            for matrix_row in value.get("result_inertia_tensor_kg_m2", [])
        ]
    except (TypeError, ValueError):
        return None

    def valid_digest(digest):
        return len(digest) == 64 and all(
            character in "0123456789abcdef" for character in digest
        )

    axes_digest = str(value.get("principal_axes_sha256", "")).lower()
    mass_digest = str(value.get("mass_property_sha256", "")).lower()
    unit = str(value.get("length_unit", "")).strip()
    if (
        not generation
        or any(
            value.get(key) != generation
            for key in (
                "transform_assembly_generation",
                "density_assembly_generation",
                "unit_assembly_generation",
                "mass_property_assembly_generation",
                "result_assembly_generation",
            )
        )
        or not part_ids
        or any(not item for item in part_ids)
        or len(set(part_ids)) != len(part_ids)
        or result_part_ids != part_ids
        or len(transforms) != len(part_ids)
        or not all(valid_digest(digest) for digest in transforms)
        or result_transforms != transforms
        or len(densities) != len(part_ids)
        or not all(math.isfinite(item) and item > 0.0 for item in densities)
        or result_densities != densities
        or unit not in {"m", "cm", "mm"}
        or value.get("result_length_unit") != unit
        or len(center) != 3
        or not all(math.isfinite(item) for item in center)
        or result_center != center
        or len(inertia) != 3
        or any(len(matrix_row) != 3 for matrix_row in inertia)
        or not all(
            math.isfinite(component)
            for matrix_row in inertia
            for component in matrix_row
        )
        or any(inertia[index][index] <= 0.0 for index in range(3))
        or any(
            not math.isclose(
                inertia[row_index][column_index],
                inertia[column_index][row_index],
                rel_tol=1.0e-12,
                abs_tol=1.0e-18,
            )
            for row_index in range(3)
            for column_index in range(3)
        )
        or result_inertia != inertia
        or not valid_digest(axes_digest)
        or value.get("result_principal_axes_sha256") != axes_digest
        or not valid_digest(mass_digest)
        or value.get("result_mass_property_sha256") != mass_digest
    ):
        return None
    return (
        generation,
        tuple(part_ids),
        tuple(transforms),
        tuple(densities),
        unit,
        tuple(center),
        tuple(tuple(matrix_row) for matrix_row in inertia),
        axes_digest,
        mass_digest,
    )


def _fillet_chamfer_topology_identity(row):
    value = row.get(
        "fillet_chamfer_topology_naming_edge_selection_fingerprint_identity"
    )
    if not isinstance(value, dict):
        return None
    generation = str(value.get("build_generation", "")).strip()
    operation_order = [
        str(item).strip() for item in value.get("operation_order", [])
    ]
    result_operation_order = [
        str(item).strip() for item in value.get("result_operation_order", [])
    ]
    try:
        edge_ids = [int(item) for item in value.get("selected_edge_ids", [])]
        result_edge_ids = [
            int(item) for item in value.get("result_selected_edge_ids", [])
        ]
    except (TypeError, ValueError):
        return None
    names = [
        str(item).strip() for item in value.get("persistent_edge_names", [])
    ]
    result_names = [
        str(item).strip()
        for item in value.get("result_persistent_edge_names", [])
    ]

    def valid_digest(digest):
        return len(digest) == 64 and all(
            character in "0123456789abcdef" for character in digest
        )

    pre_digest = str(value.get("pre_operation_topology_sha256", "")).lower()
    final_digest = str(value.get("final_topology_sha256", "")).lower()
    fingerprint = str(value.get("build_fingerprint_sha256", "")).lower()
    if (
        not generation
        or any(
            value.get(key) != generation
            for key in (
                "selection_build_generation",
                "fillet_build_generation",
                "chamfer_build_generation",
                "naming_build_generation",
                "result_build_generation",
            )
        )
        or operation_order != ["fillet", "chamfer"]
        or result_operation_order != operation_order
        or not edge_ids
        or any(item <= 0 for item in edge_ids)
        or len(set(edge_ids)) != len(edge_ids)
        or result_edge_ids != edge_ids
        or len(names) != len(edge_ids)
        or any(not item for item in names)
        or len(set(names)) != len(names)
        or result_names != names
        or not valid_digest(pre_digest)
        or value.get("result_pre_operation_topology_sha256") != pre_digest
        or not valid_digest(final_digest)
        or value.get("result_final_topology_sha256") != final_digest
        or not valid_digest(fingerprint)
        or value.get("result_build_fingerprint_sha256") != fingerprint
    ):
        return None
    return (
        generation,
        tuple(operation_order),
        tuple(edge_ids),
        tuple(names),
        pre_digest,
        final_digest,
        fingerprint,
    )


def _boolean_tolerance_healing_topology_volume_identity(row):
    value = row.get(
        "boolean_tolerance_healing_topology_volume_generation_identity"
    )
    if not isinstance(value, dict):
        return None
    generation = str(value.get("boolean_generation", "")).strip()
    operand_ids = [str(item).strip() for item in value.get("operand_ids", [])]
    result_operand_ids = [
        str(item).strip() for item in value.get("result_operand_ids", [])
    ]
    healing_policy = str(value.get("healing_policy", "")).strip()
    topology_keys = ("solids", "shells", "faces", "edges")
    topology = value.get("topology_signature")
    result_topology = value.get("result_topology_signature")
    try:
        linear_tolerance = float(value.get("linear_tolerance"))
        result_linear_tolerance = float(value.get("result_linear_tolerance"))
        volume = float(value.get("volume_m3"))
        result_volume = float(value.get("result_volume_m3"))
        normalized_topology = tuple(int(topology[key]) for key in topology_keys)
        normalized_result_topology = tuple(
            int(result_topology[key]) for key in topology_keys
        )
    except (KeyError, TypeError, ValueError):
        return None

    def valid_digest(digest):
        return len(digest) == 64 and all(
            character in "0123456789abcdef" for character in digest
        )

    operand_digest = str(value.get("operand_shape_sha256", "")).lower()
    boolean_digest = str(value.get("boolean_shape_sha256", "")).lower()
    if (
        not generation
        or any(
            value.get(key) != generation
            for key in (
                "operand_boolean_generation",
                "tolerance_boolean_generation",
                "healing_boolean_generation",
                "topology_boolean_generation",
                "volume_boolean_generation",
                "result_boolean_generation",
            )
        )
        or len(operand_ids) < 2
        or any(not item for item in operand_ids)
        or len(set(operand_ids)) != len(operand_ids)
        or result_operand_ids != operand_ids
        or not math.isfinite(linear_tolerance)
        or linear_tolerance <= 0.0
        or not math.isclose(
            result_linear_tolerance,
            linear_tolerance,
            rel_tol=0.0,
            abs_tol=1.0e-18,
        )
        or not healing_policy
        or value.get("result_healing_policy") != healing_policy
        or set(topology) != set(topology_keys)
        or set(result_topology) != set(topology_keys)
        or normalized_topology[0] < 1
        or any(item < 0 for item in normalized_topology)
        or normalized_result_topology != normalized_topology
        or not math.isfinite(volume)
        or volume <= 0.0
        or not math.isclose(result_volume, volume, rel_tol=1.0e-12, abs_tol=1.0e-18)
        or not valid_digest(operand_digest)
        or value.get("result_operand_shape_sha256") != operand_digest
        or not valid_digest(boolean_digest)
        or value.get("result_boolean_shape_sha256") != boolean_digest
    ):
        return None
    return (
        generation,
        tuple(operand_ids),
        linear_tolerance,
        healing_policy,
        normalized_topology,
        volume,
        operand_digest,
        boolean_digest,
    )


def _assembly_mate_transform_dof_loop_closure_identity(row):
    value = row.get("assembly_mate_transform_dof_loop_closure_generation_identity")
    if not isinstance(value, dict):
        return None
    generation = str(value.get("assembly_generation", "")).strip()
    mate_ids = [str(item).strip() for item in value.get("mate_ids", [])]
    result_mate_ids = [
        str(item).strip() for item in value.get("result_mate_ids", [])
    ]
    transforms = [
        str(item).lower() for item in value.get("part_transform_sha256", [])
    ]
    result_transforms = [
        str(item).lower()
        for item in value.get("result_part_transform_sha256", [])
    ]
    try:
        remaining_dof = int(value.get("remaining_dof"))
        result_remaining_dof = int(value.get("result_remaining_dof"))
        closure_residual = float(value.get("loop_closure_residual_m"))
        result_closure_residual = float(
            value.get("result_loop_closure_residual_m")
        )
    except (TypeError, ValueError):
        return None

    def valid_digest(digest):
        return len(digest) == 64 and all(
            character in "0123456789abcdef" for character in digest
        )

    solver_digest = str(value.get("kinematic_solver_sha256", "")).lower()
    pose_digest = str(value.get("assembly_pose_sha256", "")).lower()
    if (
        not generation
        or any(
            value.get(key) != generation
            for key in (
                "mate_assembly_generation",
                "transform_assembly_generation",
                "dof_assembly_generation",
                "closure_assembly_generation",
                "solver_assembly_generation",
                "result_assembly_generation",
            )
        )
        or not mate_ids
        or any(not item for item in mate_ids)
        or len(set(mate_ids)) != len(mate_ids)
        or result_mate_ids != mate_ids
        or len(transforms) != len(mate_ids)
        or not all(valid_digest(digest) for digest in transforms)
        or result_transforms != transforms
        or remaining_dof < 0
        or result_remaining_dof != remaining_dof
        or not math.isfinite(closure_residual)
        or closure_residual < 0.0
        or closure_residual > 1.0e-8
        or not math.isclose(
            result_closure_residual,
            closure_residual,
            rel_tol=0.0,
            abs_tol=1.0e-18,
        )
        or not valid_digest(solver_digest)
        or value.get("result_kinematic_solver_sha256") != solver_digest
        or not valid_digest(pose_digest)
        or value.get("result_assembly_pose_sha256") != pose_digest
    ):
        return None
    return (
        generation,
        tuple(mate_ids),
        tuple(transforms),
        remaining_dof,
        closure_residual,
        solver_digest,
        pose_digest,
    )


def _valid_identity_digest(value):
    value = str(value).lower()
    return len(value) == 64 and all(character in "0123456789abcdef" for character in value)


def _fillet_chamfer_selector_generation_identity(row):
    value = row.get(
        "fillet_chamfer_edge_selector_topology_naming_tolerance_shape_generation_identity"
    )
    if not isinstance(value, dict):
        return None
    generation = str(value.get("feature_generation", "")).strip()
    names = [str(item).strip() for item in value.get("edge_selector_names", [])]
    result_names = [str(item).strip() for item in value.get("result_edge_selector_names", [])]
    try:
        edge_ids = [int(item) for item in value.get("persistent_edge_ids", [])]
        result_edge_ids = [int(item) for item in value.get("result_persistent_edge_ids", [])]
        radius = float(value.get("feature_radius_m"))
        result_radius = float(value.get("result_feature_radius_m"))
        tolerance = float(value.get("linear_tolerance_m"))
        result_tolerance = float(value.get("result_linear_tolerance_m"))
        topology = tuple(int(value["topology_signature"][key]) for key in ("solids", "faces", "edges"))
        result_topology = tuple(int(value["result_topology_signature"][key]) for key in ("solids", "faces", "edges"))
    except (KeyError, TypeError, ValueError):
        return None
    input_digest = str(value.get("input_shape_sha256", "")).lower()
    feature_digest = str(value.get("feature_shape_sha256", "")).lower()
    if (
        not generation
        or any(value.get(key) != generation for key in (
            "selector_feature_generation", "topology_feature_generation",
            "tolerance_feature_generation", "shape_feature_generation",
            "result_feature_generation"))
        or value.get("feature_type") not in {"fillet", "chamfer"}
        or value.get("result_feature_type") != value.get("feature_type")
        or not names or any(not item for item in names) or len(set(names)) != len(names)
        or result_names != names
        or len(edge_ids) != len(names) or any(item <= 0 for item in edge_ids)
        or len(set(edge_ids)) != len(edge_ids) or result_edge_ids != edge_ids
        or not math.isfinite(radius) or radius <= 0.0 or result_radius != radius
        or not math.isfinite(tolerance) or tolerance <= 0.0 or result_tolerance != tolerance
        or topology[0] < 1 or any(item < 0 for item in topology) or result_topology != topology
        or not _valid_identity_digest(input_digest) or value.get("result_input_shape_sha256") != input_digest
        or not _valid_identity_digest(feature_digest) or value.get("result_feature_shape_sha256") != feature_digest
    ):
        return None
    return generation, tuple(names), tuple(edge_ids), radius, tolerance, topology, input_digest, feature_digest


def _mass_density_frame_generation_identity(row):
    value = row.get(
        "mass_density_center_inertia_reference_frame_assembly_generation_identity"
    )
    if not isinstance(value, dict):
        return None
    generation = str(value.get("mass_generation", "")).strip()
    try:
        density = float(value.get("density_kg_m3"))
        result_density = float(value.get("result_density_kg_m3"))
        volume = float(value.get("volume_m3"))
        result_volume = float(value.get("result_volume_m3"))
        mass = float(value.get("mass_kg"))
        result_mass = float(value.get("result_mass_kg"))
        center = tuple(float(item) for item in value.get("center_of_mass_m", []))
        result_center = tuple(float(item) for item in value.get("result_center_of_mass_m", []))
        inertia = tuple(tuple(float(item) for item in row) for row in value.get("inertia_tensor_kg_m2", []))
        result_inertia = tuple(tuple(float(item) for item in row) for row in value.get("result_inertia_tensor_kg_m2", []))
    except (TypeError, ValueError):
        return None
    frame = str(value.get("reference_frame", "")).strip()
    transform_digest = str(value.get("assembly_transform_sha256", "")).lower()
    result_digest = str(value.get("mass_property_sha256", "")).lower()
    if (
        not generation
        or any(value.get(key) != generation for key in (
            "density_mass_generation", "center_mass_generation", "inertia_mass_generation",
            "frame_mass_generation", "assembly_mass_generation", "result_mass_generation"))
        or not all(math.isfinite(item) and item > 0.0 for item in (density, volume, mass))
        or result_density != density or result_volume != volume or result_mass != mass
        or not math.isclose(mass, density * volume, rel_tol=1.0e-12, abs_tol=1.0e-15)
        or len(center) != 3 or any(not math.isfinite(item) for item in center) or result_center != center
        or len(inertia) != 3 or any(len(row) != 3 for row in inertia)
        or any(not math.isfinite(item) for row in inertia for item in row) or result_inertia != inertia
        or any(not math.isclose(inertia[i][j], inertia[j][i], rel_tol=0.0, abs_tol=1.0e-15) for i in range(3) for j in range(3))
        or not frame or value.get("result_reference_frame") != frame
        or not _valid_identity_digest(transform_digest) or value.get("result_assembly_transform_sha256") != transform_digest
        or not _valid_identity_digest(result_digest) or value.get("result_mass_property_sha256") != result_digest
    ):
        return None
    return generation, density, volume, mass, center, inertia, frame, transform_digest, result_digest


def _boolean_imprint_generation_identity(row):
    value = row.get(
        "boolean_imprint_interface_owner_topology_name_tolerance_mass_generation_identity"
    )
    if not isinstance(value, dict):
        return None
    generation = str(value.get("boolean_generation", "")).strip()
    faces = [str(item).strip() for item in value.get("interface_face_names", [])]
    result_faces = [str(item).strip() for item in value.get("result_interface_face_names", [])]
    owners = [tuple(str(item).strip() for item in pair) for pair in value.get("interface_owner_pairs", [])]
    result_owners = [tuple(str(item).strip() for item in pair) for pair in value.get("result_interface_owner_pairs", [])]
    names = [str(item).strip() for item in value.get("persistent_topology_names", [])]
    result_names = [str(item).strip() for item in value.get("result_persistent_topology_names", [])]
    try:
        tolerance = float(value.get("linear_tolerance_m"))
        result_tolerance = float(value.get("result_linear_tolerance_m"))
        solids = int(value.get("solid_count"))
        result_solids = int(value.get("result_solid_count"))
        volume = float(value.get("total_volume_m3"))
        result_volume = float(value.get("result_total_volume_m3"))
    except (TypeError, ValueError):
        return None
    interface_digest = str(value.get("interface_map_sha256", "")).lower()
    mass_digest = str(value.get("mass_property_sha256", "")).lower()
    if (
        not generation
        or any(value.get(key) != generation for key in (
            "interface_boolean_generation", "owner_boolean_generation",
            "topology_boolean_generation", "tolerance_boolean_generation",
            "mass_boolean_generation", "result_boolean_generation"))
        or value.get("operation") != "imprint"
        or value.get("result_operation") != "imprint"
        or not faces or any(not item for item in faces) or len(set(faces)) != len(faces)
        or result_faces != faces
        or len(owners) != len(faces)
        or any(len(pair) != 2 or not all(pair) or pair[0] == pair[1] for pair in owners)
        or result_owners != owners
        or not names or any(not item for item in names) or len(set(names)) != len(names)
        or result_names != names
        or not math.isfinite(tolerance) or tolerance <= 0.0 or result_tolerance != tolerance
        or solids < 1 or result_solids != solids
        or not math.isfinite(volume) or volume <= 0.0 or result_volume != volume
        or not _valid_identity_digest(interface_digest)
        or value.get("result_interface_map_sha256") != interface_digest
        or not _valid_identity_digest(mass_digest)
        or value.get("result_mass_property_sha256") != mass_digest
    ):
        return None
    return generation, tuple(faces), tuple(owners), tuple(names), tolerance, solids, volume, interface_digest, mass_digest


def _loft_section_generation_identity(row):
    value = row.get("loft_section_wire_seam_continuity_solid_volume_generation_identity")
    if not isinstance(value, dict):
        return None
    generation = str(value.get("loft_generation", "")).strip()
    sections = [str(item).strip() for item in value.get("section_names", [])]
    result_sections = [str(item).strip() for item in value.get("result_section_names", [])]
    seams = [str(item).strip() for item in value.get("seam_vertex_names", [])]
    result_seams = [str(item).strip() for item in value.get("result_seam_vertex_names", [])]
    try:
        orientations = [int(item) for item in value.get("wire_orientation_signs", [])]
        result_orientations = [int(item) for item in value.get("result_wire_orientation_signs", [])]
        volume = float(value.get("volume_m3"))
        result_volume = float(value.get("result_volume_m3"))
    except (TypeError, ValueError):
        return None
    digest = str(value.get("loft_shape_sha256", "")).lower()
    continuity = str(value.get("continuity", ""))
    if (
        not generation
        or any(value.get(key) != generation for key in (
            "section_loft_generation", "wire_loft_generation", "seam_loft_generation",
            "continuity_loft_generation", "solid_loft_generation",
            "volume_loft_generation", "result_loft_generation"))
        or len(sections) < 2 or any(not item for item in sections)
        or len(set(sections)) != len(sections) or result_sections != sections
        or len(orientations) != len(sections) or any(item != 1 for item in orientations)
        or result_orientations != orientations
        or len(seams) != len(sections) or any(not item for item in seams)
        or len(set(seams)) != len(seams) or result_seams != seams
        or continuity not in {"C1", "C2"} or value.get("result_continuity") != continuity
        or value.get("is_valid_solid") is not True
        or value.get("result_is_valid_solid") is not True
        or not math.isfinite(volume) or volume <= 0.0 or result_volume != volume
        or not _valid_identity_digest(digest)
        or value.get("result_loft_shape_sha256") != digest
    ):
        return None
    return generation, tuple(sections), tuple(orientations), tuple(seams), continuity, volume, digest


def _shell_offset_generation_identity(row):
    value = row.get(
        "shell_offset_face_normal_thickness_join_self_intersection_mass_generation_identity"
    )
    if not isinstance(value, dict):
        return None
    generation = str(value.get("shell_generation", "")).strip()
    faces = [str(item).strip() for item in value.get("selected_face_names", [])]
    result_faces = [
        str(item).strip() for item in value.get("result_selected_face_names", [])
    ]
    try:
        thickness = float(value.get("thickness_m"))
        result_thickness = float(value.get("result_thickness_m"))
        volume = float(value.get("volume_m3"))
        result_volume = float(value.get("result_volume_m3"))
    except (TypeError, ValueError):
        return None
    direction = str(value.get("normal_direction", "")).strip()
    join_mode = str(value.get("join_mode", "")).strip()
    shape_generation = str(value.get("shape_generation", "")).strip()
    digest = str(value.get("mass_property_sha256", "")).lower()
    if (
        not generation
        or any(value.get(key) != generation for key in (
            "face_shell_generation", "normal_shell_generation",
            "thickness_shell_generation", "join_shell_generation",
            "intersection_shell_generation", "shape_shell_generation",
            "mass_shell_generation", "result_shell_generation"))
        or not faces or any(not item for item in faces) or len(set(faces)) != len(faces)
        or result_faces != faces
        or direction not in {"outward", "inward"}
        or value.get("result_normal_direction") != direction
        or not math.isfinite(thickness) or thickness <= 0.0 or result_thickness != thickness
        or join_mode not in {"arc", "intersection", "tangent"}
        or value.get("result_join_mode") != join_mode
        or value.get("self_intersection") is not False
        or value.get("result_self_intersection") is not False
        or not shape_generation
        or value.get("result_shape_generation") != shape_generation
        or value.get("is_valid_solid") is not True
        or value.get("result_is_valid_solid") is not True
        or not math.isfinite(volume) or volume <= 0.0 or result_volume != volume
        or not _valid_identity_digest(digest)
        or value.get("result_mass_property_sha256") != digest
    ):
        return None
    return (
        generation, tuple(faces), direction, thickness, join_mode,
        shape_generation, volume, digest,
    )


def _path_sweep_generation_identity(row):
    value = row.get(
        "path_sweep_frame_transition_profile_orientation_solid_volume_generation_identity"
    )
    if not isinstance(value, dict):
        return None
    generation = str(value.get("sweep_generation", "")).strip()
    paths = [str(item).strip() for item in value.get("path_edge_names", [])]
    result_paths = [str(item).strip() for item in value.get("result_path_edge_names", [])]
    profiles = [str(item).strip() for item in value.get("profile_wire_names", [])]
    result_profiles = [str(item).strip() for item in value.get("result_profile_wire_names", [])]
    try:
        orientations = [float(item) for item in value.get("profile_orientation_deg", [])]
        result_orientations = [float(item) for item in value.get("result_profile_orientation_deg", [])]
        volume = float(value.get("volume_m3"))
        result_volume = float(value.get("result_volume_m3"))
    except (TypeError, ValueError):
        return None
    frame = str(value.get("moving_frame", "")).strip()
    transition = str(value.get("transition_mode", "")).strip()
    digest = str(value.get("sweep_shape_sha256", "")).lower()
    if (
        not generation
        or any(value.get(key) != generation for key in (
            "path_sweep_generation", "frame_sweep_generation",
            "transition_sweep_generation", "profile_sweep_generation",
            "orientation_sweep_generation", "solid_sweep_generation",
            "volume_sweep_generation", "result_sweep_generation"))
        or not paths or any(not item for item in paths) or len(set(paths)) != len(paths)
        or result_paths != paths
        or frame not in {"parallel_transport", "frenet", "fixed"}
        or value.get("result_moving_frame") != frame
        or transition not in {"round", "right", "transformed"}
        or value.get("result_transition_mode") != transition
        or not profiles or any(not item for item in profiles) or len(set(profiles)) != len(profiles)
        or result_profiles != profiles
        or len(orientations) != len(paths)
        or any(not math.isfinite(item) for item in orientations)
        or result_orientations != orientations
        or value.get("is_valid_solid") is not True
        or value.get("result_is_valid_solid") is not True
        or not math.isfinite(volume) or volume <= 0.0 or result_volume != volume
        or not _valid_identity_digest(digest)
        or value.get("result_sweep_shape_sha256") != digest
    ):
        return None
    return generation, tuple(paths), frame, transition, tuple(profiles), tuple(orientations), volume, digest


def _sheet_metal_flat_pattern_generation_identity(row):
    value = row.get(
        "sheet_metal_bend_allowance_kfactor_neutral_axis_relief_thickness_flat_pattern_area_generation_identity"
    )
    if not isinstance(value, dict):
        return None
    generation = str(value.get("sheet_generation", "")).strip()
    try:
        radius = float(value.get("bend_radius_m"))
        result_radius = float(value.get("result_bend_radius_m"))
        angle = float(value.get("bend_angle_deg"))
        result_angle = float(value.get("result_bend_angle_deg"))
        k_factor = float(value.get("k_factor"))
        result_k_factor = float(value.get("result_k_factor"))
        neutral_radius = float(value.get("neutral_axis_radius_m"))
        result_neutral_radius = float(value.get("result_neutral_axis_radius_m"))
        allowance = float(value.get("bend_allowance_m"))
        result_allowance = float(value.get("result_bend_allowance_m"))
        relief_width = float(value.get("relief_width_m"))
        result_relief_width = float(value.get("result_relief_width_m"))
        thickness = float(value.get("thickness_m"))
        result_thickness = float(value.get("result_thickness_m"))
        width = float(value.get("strip_width_m"))
        result_width = float(value.get("result_strip_width_m"))
        straight = [float(item) for item in value.get("straight_lengths_m", [])]
        result_straight = [
            float(item) for item in value.get("result_straight_lengths_m", [])
        ]
        area = float(value.get("flat_pattern_area_m2"))
        result_area = float(value.get("result_flat_pattern_area_m2"))
    except (TypeError, ValueError):
        return None
    relief = str(value.get("relief_type", "")).strip()
    digest = str(value.get("flat_pattern_sha256", "")).lower()
    expected_neutral_radius = radius + k_factor * thickness
    expected_allowance = math.radians(abs(angle)) * expected_neutral_radius
    expected_area = width * (sum(straight) + expected_allowance)
    if (
        not generation
        or any(value.get(key) != generation for key in (
            "bend_sheet_generation", "neutral_axis_sheet_generation",
            "relief_sheet_generation", "thickness_sheet_generation",
            "pattern_sheet_generation", "area_sheet_generation",
            "result_sheet_generation"))
        or not all(math.isfinite(item) for item in (
            radius, angle, k_factor, neutral_radius, allowance, relief_width,
            thickness, width, area))
        or radius <= 0.0 or thickness <= 0.0 or width <= 0.0
        or not 0.0 < abs(angle) <= 180.0 or not 0.0 <= k_factor <= 1.0
        or relief not in {"rectangular", "round", "tear"}
        or relief_width <= 0.0
        or not straight or any(item <= 0.0 or not math.isfinite(item) for item in straight)
        or not math.isclose(neutral_radius, expected_neutral_radius, rel_tol=1.0e-12, abs_tol=1.0e-15)
        or not math.isclose(allowance, expected_allowance, rel_tol=1.0e-12, abs_tol=1.0e-15)
        or not math.isclose(area, expected_area, rel_tol=1.0e-12, abs_tol=1.0e-15)
        or result_radius != radius or result_angle != angle
        or result_k_factor != k_factor or result_neutral_radius != neutral_radius
        or result_allowance != allowance or value.get("result_relief_type") != relief
        or result_relief_width != relief_width or result_thickness != thickness
        or result_width != width or result_straight != straight or result_area != area
        or value.get("flat_pattern_wire_closed") is not True
        or value.get("result_flat_pattern_wire_closed") is not True
        or not _valid_identity_digest(digest)
        or value.get("result_flat_pattern_sha256") != digest
    ):
        return None
    return (
        generation, radius, angle, k_factor, neutral_radius, allowance, relief,
        relief_width, thickness, width, tuple(straight), area, digest,
    )


def _joint_kinematic_loop_generation_identity(row):
    value = row.get(
        "joint_kinematic_loop_graph_dof_limit_connector_frame_closure_configuration_swept_volume_generation_identity"
    )
    if not isinstance(value, dict):
        return None
    generation = str(value.get("joint_generation", "")).strip()
    edges = [tuple(str(item).strip() for item in edge) for edge in value.get("joint_graph_edges", [])]
    result_edges = [
        tuple(str(item).strip() for item in edge)
        for edge in value.get("result_joint_graph_edges", [])
    ]
    names = [str(item).strip() for item in value.get("dof_names", [])]
    result_names = [str(item).strip() for item in value.get("result_dof_names", [])]
    types = [str(item).strip() for item in value.get("dof_types", [])]
    result_types = [str(item).strip() for item in value.get("result_dof_types", [])]
    try:
        lower = [float(item) for item in value.get("lower_limits", [])]
        result_lower = [float(item) for item in value.get("result_lower_limits", [])]
        upper = [float(item) for item in value.get("upper_limits", [])]
        result_upper = [float(item) for item in value.get("result_upper_limits", [])]
        configuration = [float(item) for item in value.get("configuration_values", [])]
        result_configuration = [
            float(item) for item in value.get("result_configuration_values", [])
        ]
        closure = float(value.get("loop_closure_error_m"))
        result_closure = float(value.get("result_loop_closure_error_m"))
        tolerance = float(value.get("loop_closure_tolerance_m"))
        result_tolerance = float(value.get("result_loop_closure_tolerance_m"))
        swept_volume = float(value.get("swept_volume_m3"))
        result_swept_volume = float(value.get("result_swept_volume_m3"))
    except (TypeError, ValueError):
        return None
    nodes = {item for edge in edges if len(edge) == 3 for item in (edge[0], edge[2])}
    joint_names = [edge[1] for edge in edges if len(edge) == 3]
    graph = {node: set() for node in nodes}
    for edge in edges:
        if len(edge) == 3 and edge[0] in graph and edge[2] in graph:
            graph[edge[0]].add(edge[2])
            graph[edge[2]].add(edge[0])
    visited = set()
    pending = [next(iter(nodes))] if nodes else []
    while pending:
        node = pending.pop()
        if node in visited:
            continue
        visited.add(node)
        pending.extend(graph[node] - visited)
    frame_digest = str(value.get("connector_frame_sha256", "")).lower()
    swept_digest = str(value.get("swept_shape_sha256", "")).lower()
    configuration_id = str(value.get("configuration_id", "")).strip()
    if (
        not generation
        or any(value.get(key) != generation for key in (
            "graph_joint_generation", "dof_joint_generation", "limit_joint_generation",
            "frame_joint_generation", "closure_joint_generation",
            "configuration_joint_generation", "swept_joint_generation",
            "result_joint_generation"))
        or len(edges) < 3 or any(len(edge) != 3 or not all(edge) for edge in edges)
        or len(set(edges)) != len(edges) or len(set(joint_names)) != len(joint_names)
        or visited != nodes or len(edges) < len(nodes)
        or result_edges != edges
        or not names or len(set(names)) != len(names) or result_names != names
        or len(types) != len(names) or any(item not in {"revolute", "prismatic"} for item in types)
        or result_types != types
        or not (len(lower) == len(upper) == len(configuration) == len(names))
        or any(not all(math.isfinite(item) for item in values) for values in (lower, upper, configuration))
        or any(lo >= hi or not lo <= current <= hi for lo, current, hi in zip(lower, configuration, upper))
        or result_lower != lower or result_upper != upper or result_configuration != configuration
        or not _valid_identity_digest(frame_digest)
        or value.get("result_connector_frame_sha256") != frame_digest
        or not all(math.isfinite(item) for item in (closure, tolerance, swept_volume))
        or tolerance <= 0.0 or closure < 0.0 or closure > tolerance
        or result_closure != closure or result_tolerance != tolerance
        or not configuration_id or value.get("result_configuration_id") != configuration_id
        or swept_volume <= 0.0 or result_swept_volume != swept_volume
        or not _valid_identity_digest(swept_digest)
        or value.get("result_swept_shape_sha256") != swept_digest
    ):
        return None
    return (
        generation, tuple(edges), tuple(names), tuple(types), tuple(lower),
        tuple(upper), tuple(configuration), frame_digest, closure, tolerance,
        configuration_id, swept_volume, swept_digest,
    )


def _helical_sweep_generation_identity(row):
    value = row.get(
        "helical_sweep_pitch_handedness_profile_frame_turn_self_intersection_volume_centroid_shape_generation_identity"
    )
    if not isinstance(value, dict):
        return None
    generation = str(value.get("helical_generation", "")).strip()
    try:
        pitch = float(value.get("pitch_m"))
        result_pitch = float(value.get("result_pitch_m"))
        frame = [
            [float(item) for item in row]
            for row in value.get("profile_frame_matrix", [])
        ]
        result_frame = [
            [float(item) for item in row]
            for row in value.get("result_profile_frame_matrix", [])
        ]
        turns = float(value.get("turn_count"))
        result_turns = float(value.get("result_turn_count"))
        volume = float(value.get("volume_m3"))
        result_volume = float(value.get("result_volume_m3"))
        centroid = [float(item) for item in value.get("centroid_m", [])]
        result_centroid = [
            float(item) for item in value.get("result_centroid_m", [])
        ]
    except (TypeError, ValueError):
        return None
    handedness = str(value.get("handedness", "")).strip()
    profile_digest = str(value.get("profile_frame_sha256", "")).lower()
    shape_digest = str(value.get("helical_shape_sha256", "")).lower()
    intersection = value.get("self_intersection")
    result_intersection = value.get("result_self_intersection")
    if (
        not generation
        or any(value.get(key) != generation for key in (
            "pitch_helical_generation", "handedness_helical_generation",
            "profile_helical_generation", "turn_helical_generation",
            "intersection_helical_generation", "mass_helical_generation",
            "shape_helical_generation", "result_helical_generation"))
        or not math.isfinite(pitch) or pitch <= 0.0 or result_pitch != pitch
        or handedness not in {"right", "left"}
        or value.get("result_handedness") != handedness
        or len(frame) != 4
        or any(len(items) != 4 or any(not math.isfinite(item) for item in items) for items in frame)
        or frame[3] != [0.0, 0.0, 0.0, 1.0]
        or result_frame != frame
        or not _valid_identity_digest(profile_digest)
        or value.get("result_profile_frame_sha256") != profile_digest
        or not math.isfinite(turns) or turns <= 0.0 or result_turns != turns
        or not isinstance(intersection, bool)
        or not isinstance(result_intersection, bool)
        or intersection is not False or result_intersection != intersection
        or not math.isfinite(volume) or volume <= 0.0 or result_volume != volume
        or len(centroid) != 3 or any(not math.isfinite(item) for item in centroid)
        or result_centroid != centroid
        or not _valid_identity_digest(shape_digest)
        or value.get("result_helical_shape_sha256") != shape_digest
    ):
        return None
    return (
        generation, pitch, handedness, tuple(tuple(items) for items in frame),
        profile_digest, turns, intersection, volume, tuple(centroid), shape_digest,
    )


def _boolean_history_generation_identity(row):
    value = row.get(
        "boolean_tolerance_operand_order_history_volume_centroid_inertia_shape_generation_identity"
    )
    if not isinstance(value, dict):
        return None
    generation = str(value.get("boolean_generation", "")).strip()
    operands = [str(item).strip() for item in value.get("operand_order", [])]
    result_operands = [
        str(item).strip() for item in value.get("result_operand_order", [])
    ]
    try:
        tolerance = float(value.get("model_tolerance_m"))
        result_tolerance = float(value.get("result_model_tolerance_m"))
        volume = float(value.get("volume_m3"))
        result_volume = float(value.get("result_volume_m3"))
        centroid = [float(item) for item in value.get("centroid_m", [])]
        result_centroid = [
            float(item) for item in value.get("result_centroid_m", [])
        ]
        inertia = [
            [float(item) for item in row]
            for row in value.get("inertia_tensor_kg_m2", [])
        ]
        result_inertia = [
            [float(item) for item in row]
            for row in value.get("result_inertia_tensor_kg_m2", [])
        ]
    except (TypeError, ValueError):
        return None
    operation = str(value.get("operation", "")).strip()
    history_digest = str(value.get("subshape_history_sha256", "")).lower()
    shape_digest = str(value.get("boolean_shape_sha256", "")).lower()
    if (
        not generation
        or any(value.get(key) != generation for key in (
            "tolerance_boolean_generation", "operand_boolean_generation",
            "history_boolean_generation", "mass_boolean_generation",
            "inertia_boolean_generation", "shape_boolean_generation",
            "result_boolean_generation"))
        or operation not in {"cut", "fuse", "intersect"}
        or value.get("result_operation") != operation
        or not math.isfinite(tolerance) or tolerance <= 0.0
        or result_tolerance != tolerance
        or len(operands) != 2 or not all(operands) or operands[0] == operands[1]
        or result_operands != operands
        or not _valid_identity_digest(history_digest)
        or value.get("result_subshape_history_sha256") != history_digest
        or not math.isfinite(volume) or volume <= 0.0 or result_volume != volume
        or len(centroid) != 3 or any(not math.isfinite(item) for item in centroid)
        or result_centroid != centroid
        or len(inertia) != 3
        or any(len(items) != 3 or any(not math.isfinite(item) for item in items) for items in inertia)
        or any(inertia[index][index] <= 0.0 for index in range(3))
        or any(not math.isclose(inertia[i][j], inertia[j][i], rel_tol=0.0, abs_tol=1.0e-18) for i in range(3) for j in range(3))
        or result_inertia != inertia
        or not _valid_identity_digest(shape_digest)
        or value.get("result_boolean_shape_sha256") != shape_digest
    ):
        return None
    return (
        generation, operation, tolerance, tuple(operands), history_digest,
        volume, tuple(centroid), tuple(tuple(items) for items in inertia), shape_digest,
    )


def _assembly_occurrence_mass_generation_identity(row):
    value = row.get(
        "assembly_occurrence_location_density_unit_suppression_mass_center_inertia_parallel_axis_shape_generation_identity"
    )
    if not isinstance(value, dict):
        return None
    generation = str(value.get("assembly_generation", "")).strip()
    occurrence_ids = [str(item).strip() for item in value.get("occurrence_ids", [])]
    result_occurrence_ids = [
        str(item).strip() for item in value.get("result_occurrence_ids", [])
    ]
    suppressed = [
        str(item).strip() for item in value.get("suppressed_occurrence_ids", [])
    ]
    result_suppressed = [
        str(item).strip()
        for item in value.get("result_suppressed_occurrence_ids", [])
    ]
    try:
        locations = [
            [[float(item) for item in matrix_row] for matrix_row in matrix]
            for matrix in value.get("location_matrices", [])
        ]
        result_locations = [
            [[float(item) for item in matrix_row] for matrix_row in matrix]
            for matrix in value.get("result_location_matrices", [])
        ]
        densities = [float(item) for item in value.get("densities_kg_m3", [])]
        result_densities = [
            float(item) for item in value.get("result_densities_kg_m3", [])
        ]
        volumes = [float(item) for item in value.get("part_volumes_m3", [])]
        result_volumes = [
            float(item) for item in value.get("result_part_volumes_m3", [])
        ]
        mass = float(value.get("assembly_mass_kg"))
        result_mass = float(value.get("result_assembly_mass_kg"))
        center = [float(item) for item in value.get("center_of_mass_m", [])]
        result_center = [
            float(item) for item in value.get("result_center_of_mass_m", [])
        ]
        inertia = [
            [float(item) for item in matrix_row]
            for matrix_row in value.get("assembly_inertia_kg_m2", [])
        ]
        result_inertia = [
            [float(item) for item in matrix_row]
            for matrix_row in value.get("result_assembly_inertia_kg_m2", [])
        ]
    except (TypeError, ValueError):
        return None
    parallel_axis = value.get("parallel_axis_applied")
    result_parallel_axis = value.get("result_parallel_axis_applied")
    shape_digest = str(value.get("assembly_shape_sha256", "")).lower()
    part_masses = [density * volume for density, volume in zip(densities, volumes)]
    calculated_mass = sum(part_masses)
    calculated_center = [
        sum(
            part_mass * locations[index][axis][3]
            for index, part_mass in enumerate(part_masses)
        )
        / calculated_mass
        for axis in range(3)
    ] if calculated_mass > 0.0 else []
    if (
        not generation
        or any(value.get(key) != generation for key in (
            "occurrence_assembly_generation", "location_assembly_generation",
            "density_assembly_generation", "suppression_assembly_generation",
            "mass_assembly_generation", "inertia_assembly_generation",
            "shape_assembly_generation", "result_assembly_generation"))
        or not occurrence_ids or len(set(occurrence_ids)) != len(occurrence_ids)
        or result_occurrence_ids != occurrence_ids
        or not (len(locations) == len(densities) == len(volumes) == len(occurrence_ids))
        or any(len(matrix) != 4 or any(len(matrix_row) != 4 for matrix_row in matrix) for matrix in locations)
        or any(not math.isfinite(item) for matrix in locations for matrix_row in matrix for item in matrix_row)
        or any(matrix[3] != [0.0, 0.0, 0.0, 1.0] for matrix in locations)
        or result_locations != locations
        or any(not math.isfinite(item) or item <= 0.0 for item in densities + volumes)
        or result_densities != densities or result_volumes != volumes
        or value.get("density_unit") != "kg/m^3"
        or value.get("result_density_unit") != value.get("density_unit")
        or any(item not in occurrence_ids for item in suppressed)
        or result_suppressed != suppressed
        or suppressed
        or not math.isfinite(mass) or mass <= 0.0 or result_mass != mass
        or not math.isclose(mass, calculated_mass, rel_tol=1.0e-12, abs_tol=1.0e-15)
        or len(center) != 3 or any(not math.isfinite(item) for item in center)
        or result_center != center
        or any(not math.isclose(a, b, rel_tol=1.0e-12, abs_tol=1.0e-15) for a, b in zip(center, calculated_center))
        or value.get("inertia_reference_frame") != "assembly_global_center_of_mass"
        or value.get("result_inertia_reference_frame") != value.get("inertia_reference_frame")
        or parallel_axis is not True or result_parallel_axis != parallel_axis
        or len(inertia) != 3
        or any(len(matrix_row) != 3 or any(not math.isfinite(item) for item in matrix_row) for matrix_row in inertia)
        or any(inertia[index][index] <= 0.0 for index in range(3))
        or any(not math.isclose(inertia[i][j], inertia[j][i], rel_tol=0.0, abs_tol=1.0e-18) for i in range(3) for j in range(3))
        or result_inertia != inertia
        or not _valid_identity_digest(shape_digest)
        or value.get("result_assembly_shape_sha256") != shape_digest
    ):
        return None
    return (
        generation, tuple(occurrence_ids), tuple(tuple(tuple(items) for items in matrix) for matrix in locations),
        tuple(densities), tuple(volumes), mass, tuple(center),
        tuple(tuple(items) for items in inertia), shape_digest,
    )


def _loft_face_lineage_generation_identity(row):
    value = row.get(
        "loft_sweep_profile_order_seam_guide_orientation_face_lineage_shell_volume_shape_generation_identity"
    )
    if not isinstance(value, dict):
        return None
    generation = str(value.get("loft_generation", "")).strip()
    profiles = [str(item).strip() for item in value.get("profile_ids", [])]
    result_profiles = [str(item).strip() for item in value.get("result_profile_ids", [])]
    lineage = [tuple(str(item).strip() for item in pair) for pair in value.get("face_lineage_pairs", [])]
    result_lineage = [tuple(str(item).strip() for item in pair) for pair in value.get("result_face_lineage_pairs", [])]
    try:
        parameters = [float(item) for item in value.get("profile_parameters", [])]
        result_parameters = [float(item) for item in value.get("result_profile_parameters", [])]
        seams = [float(item) for item in value.get("seam_parameters", [])]
        result_seams = [float(item) for item in value.get("result_seam_parameters", [])]
        volume = float(value.get("volume_m3"))
        result_volume = float(value.get("result_volume_m3"))
    except (TypeError, ValueError):
        return None
    shell_closed = value.get("shell_closed")
    result_shell_closed = value.get("result_shell_closed")
    lineage_digest = str(value.get("face_lineage_sha256", "")).lower()
    shape_digest = str(value.get("loft_shape_sha256", "")).lower()
    if (
        not generation
        or any(value.get(key) != generation for key in (
            "profile_loft_generation", "seam_loft_generation",
            "guide_loft_generation", "lineage_loft_generation",
            "shell_loft_generation", "mass_loft_generation",
            "shape_loft_generation", "result_loft_generation"))
        or len(profiles) < 2 or len(set(profiles)) != len(profiles)
        or result_profiles != profiles
        or len(parameters) != len(profiles)
        or any(not math.isfinite(item) for item in parameters)
        or any(right <= left for left, right in zip(parameters, parameters[1:]))
        or result_parameters != parameters
        or len(seams) != len(profiles)
        or any(not math.isfinite(item) or not 0.0 <= item < 1.0 for item in seams)
        or result_seams != seams
        or value.get("guide_orientation") != "start_to_end_right_handed"
        or value.get("result_guide_orientation") != value.get("guide_orientation")
        or len(lineage) < len(profiles)
        or any(len(pair) != 2 or not all(pair) for pair in lineage)
        or result_lineage != lineage
        or not _valid_identity_digest(lineage_digest)
        or value.get("result_face_lineage_sha256") != lineage_digest
        or shell_closed is not True or result_shell_closed != shell_closed
        or not math.isfinite(volume) or volume <= 0.0 or result_volume != volume
        or not _valid_identity_digest(shape_digest)
        or value.get("result_loft_shape_sha256") != shape_digest
    ):
        return None
    return (
        generation, tuple(profiles), tuple(parameters), tuple(seams),
        tuple(lineage), lineage_digest, volume, shape_digest,
    )


def _boolean_topology_ancestry_generation_identity(row):
    value = row.get(
        "boolean_fuzzy_tolerance_topology_name_face_ancestry_count_volume_centroid_shape_generation_identity"
    )
    if not isinstance(value, dict):
        return None
    generation = str(value.get("boolean_generation", "")).strip()
    names = [str(item).strip() for item in value.get("surviving_topology_names", [])]
    result_names = [
        str(item).strip() for item in value.get("result_surviving_topology_names", [])
    ]
    ancestry = [
        tuple(str(item).strip() for item in pair)
        for pair in value.get("face_ancestry", [])
    ]
    result_ancestry = [
        tuple(str(item).strip() for item in pair)
        for pair in value.get("result_face_ancestry", [])
    ]
    try:
        tolerance = float(value.get("fuzzy_tolerance_m"))
        result_tolerance = float(value.get("result_fuzzy_tolerance_m"))
        solid_count = int(value.get("solid_count"))
        result_solid_count = int(value.get("result_solid_count"))
        volume = float(value.get("volume_m3"))
        result_volume = float(value.get("result_volume_m3"))
        centroid = tuple(float(item) for item in value.get("centroid_m", []))
        result_centroid = tuple(
            float(item) for item in value.get("result_centroid_m", [])
        )
    except (TypeError, ValueError):
        return None
    shape_digest = str(value.get("boolean_shape_sha256", "")).lower()
    if (
        not generation
        or any(
            value.get(key) != generation
            for key in (
                "tolerance_generation",
                "topology_generation",
                "ancestry_generation",
                "mass_generation",
                "shape_generation",
                "result_generation",
            )
        )
        or not math.isfinite(tolerance)
        or tolerance <= 0.0
        or result_tolerance != tolerance
        or not names
        or len(set(names)) != len(names)
        or result_names != names
        or not ancestry
        or any(len(pair) != 2 or not all(pair) for pair in ancestry)
        or result_ancestry != ancestry
        or solid_count != 1
        or result_solid_count != solid_count
        or not math.isfinite(volume)
        or volume <= 0.0
        or result_volume != volume
        or len(centroid) != 3
        or any(not math.isfinite(item) for item in centroid)
        or result_centroid != centroid
        or not _valid_identity_digest(shape_digest)
        or value.get("result_boolean_shape_sha256") != shape_digest
    ):
        return None
    return (
        generation,
        tolerance,
        tuple(names),
        tuple(ancestry),
        solid_count,
        volume,
        centroid,
        shape_digest,
    )


def _sweep_frame_transition_generation_identity(row):
    value = row.get(
        "sweep_frenet_frame_twist_transition_profile_orientation_self_intersection_volume_owner_shape_generation_identity"
    )
    if not isinstance(value, dict):
        return None
    generation = str(value.get("sweep_generation", "")).strip()
    try:
        twist = tuple(float(item) for item in value.get("twist_parameters_rad", []))
        result_twist = tuple(
            float(item) for item in value.get("result_twist_parameters_rad", [])
        )
        orientation = tuple(
            int(item) for item in value.get("profile_orientation_signs", [])
        )
        result_orientation = tuple(
            int(item) for item in value.get("result_profile_orientation_signs", [])
        )
        volume = float(value.get("sweep_volume_m3"))
        result_volume = float(value.get("result_sweep_volume_m3"))
    except (TypeError, ValueError):
        return None
    shape_digest = str(value.get("sweep_shape_sha256", "")).lower()
    if (
        not generation
        or any(
            value.get(key) != generation
            for key in (
                "frame_generation",
                "twist_generation",
                "transition_generation",
                "orientation_generation",
                "intersection_generation",
                "mass_generation",
                "owner_generation",
                "result_generation",
            )
        )
        or value.get("frame_convention") != "corrected_frenet"
        or value.get("result_frame_convention") != value.get("frame_convention")
        or len(twist) < 2
        or any(not math.isfinite(item) for item in twist)
        or any(right < left for left, right in zip(twist, twist[1:]))
        or result_twist != twist
        or value.get("transition_mode") not in {"right_corner", "round_corner"}
        or value.get("result_transition_mode") != value.get("transition_mode")
        or len(orientation) != len(twist)
        or any(item != 1 for item in orientation)
        or result_orientation != orientation
        or value.get("self_intersection") is not False
        or value.get("result_self_intersection") is not False
        or not math.isfinite(volume)
        or volume <= 0.0
        or result_volume != volume
        or not str(value.get("shape_owner", "")).strip()
        or value.get("result_shape_owner") != value.get("shape_owner")
        or not _valid_identity_digest(shape_digest)
        or value.get("result_sweep_shape_sha256") != shape_digest
    ):
        return None
    return (
        generation,
        value.get("frame_convention"),
        twist,
        value.get("transition_mode"),
        orientation,
        volume,
        value.get("shape_owner"),
        shape_digest,
    )


def _guided_loft_generation_identity(row):
    value = row.get(
        "loft_section_guide_parameterization_seam_orientation_mode_intersection_volume_shape_generation_identity"
    )
    if not isinstance(value, dict):
        return None
    generation = str(value.get("loft_generation", "")).strip()
    sections = tuple(str(item).strip() for item in value.get("section_order", []))
    result_sections = tuple(
        str(item).strip() for item in value.get("result_section_order", [])
    )
    intersections = tuple(
        tuple(str(item).strip() for item in pair)
        for pair in value.get("guide_intersections", [])
    )
    result_intersections = tuple(
        tuple(str(item).strip() for item in pair)
        for pair in value.get("result_guide_intersections", [])
    )
    try:
        seams = tuple(int(item) for item in value.get("seam_orientation_signs", []))
        result_seams = tuple(
            int(item) for item in value.get("result_seam_orientation_signs", [])
        )
        volume = float(value.get("loft_volume_m3"))
        result_volume = float(value.get("result_loft_volume_m3"))
    except (TypeError, ValueError):
        return None
    parameter_digest = str(value.get("wire_parameterization_sha256", "")).lower()
    shape_digest = str(value.get("loft_shape_sha256", "")).lower()
    if (
        not generation
        or any(
            value.get(key) != generation
            for key in (
                "section_generation",
                "parameter_generation",
                "guide_generation",
                "seam_generation",
                "mode_generation",
                "intersection_generation",
                "mass_generation",
                "shape_generation",
                "result_generation",
            )
        )
        or len(sections) < 2
        or len(set(sections)) != len(sections)
        or result_sections != sections
        or not _valid_identity_digest(parameter_digest)
        or value.get("result_wire_parameterization_sha256") != parameter_digest
        or not intersections
        or len(set(intersections)) != len(intersections)
        or any(len(pair) != 2 or not all(pair) or pair[1] not in sections for pair in intersections)
        or {pair[1] for pair in intersections} != set(sections)
        or result_intersections != intersections
        or len(seams) != len(sections)
        or any(item != 1 for item in seams)
        or result_seams != seams
        or value.get("loft_mode") != "smooth"
        or value.get("result_loft_mode") != value.get("loft_mode")
        or value.get("self_intersection") is not False
        or value.get("result_self_intersection") is not False
        or not math.isfinite(volume)
        or volume <= 0.0
        or result_volume != volume
        or not _valid_identity_digest(shape_digest)
        or value.get("result_loft_shape_sha256") != shape_digest
    ):
        return None
    return (
        generation,
        sections,
        parameter_digest,
        intersections,
        seams,
        value.get("loft_mode"),
        volume,
        shape_digest,
    )


def _mass_inertia_parallel_axis_generation_identity(row):
    value = row.get(
        "mass_property_density_unit_origin_center_principal_axis_degeneracy_parallel_axis_owner_shape_generation_identity"
    )
    if not isinstance(value, dict):
        return None
    generation = str(value.get("mass_property_generation", "")).strip()
    try:
        density = float(value.get("density_kg_m3"))
        result_density = float(value.get("result_density_kg_m3"))
        mass = float(value.get("mass_kg"))
        result_mass = float(value.get("result_mass_kg"))
        origin = tuple(float(item) for item in value.get("inertia_origin_m", []))
        result_origin = tuple(
            float(item) for item in value.get("result_inertia_origin_m", [])
        )
        center = tuple(float(item) for item in value.get("center_of_mass_m", []))
        result_center = tuple(
            float(item) for item in value.get("result_center_of_mass_m", [])
        )
        inertia_com = tuple(
            tuple(float(item) for item in matrix_row)
            for matrix_row in value.get("inertia_at_com_kg_m2", [])
        )
        result_inertia_com = tuple(
            tuple(float(item) for item in matrix_row)
            for matrix_row in value.get("result_inertia_at_com_kg_m2", [])
        )
        inertia_origin = tuple(
            tuple(float(item) for item in matrix_row)
            for matrix_row in value.get("inertia_at_origin_kg_m2", [])
        )
        result_inertia_origin = tuple(
            tuple(float(item) for item in matrix_row)
            for matrix_row in value.get("result_inertia_at_origin_kg_m2", [])
        )
        moments = tuple(
            float(item) for item in value.get("principal_moments_kg_m2", [])
        )
        result_moments = tuple(
            float(item) for item in value.get("result_principal_moments_kg_m2", [])
        )
        axes = tuple(
            tuple(float(item) for item in matrix_row)
            for matrix_row in value.get("principal_axes", [])
        )
        result_axes = tuple(
            tuple(float(item) for item in matrix_row)
            for matrix_row in value.get("result_principal_axes", [])
        )
    except (TypeError, ValueError):
        return None

    def matrix3(matrix):
        return len(matrix) == 3 and all(
            len(matrix_row) == 3
            and all(math.isfinite(item) for item in matrix_row)
            for matrix_row in matrix
        )

    def close(left, right):
        return math.isclose(left, right, rel_tol=1.0e-10, abs_tol=1.0e-12)

    axes_orthonormal = matrix3(axes) and all(
        close(
            sum(axes[row][index] * axes[column][index] for index in range(3)),
            1.0 if row == column else 0.0,
        )
        for row in range(3)
        for column in range(3)
    )
    axes_determinant = (
        axes[0][0] * (axes[1][1] * axes[2][2] - axes[1][2] * axes[2][1])
        - axes[0][1] * (axes[1][0] * axes[2][2] - axes[1][2] * axes[2][0])
        + axes[0][2] * (axes[1][0] * axes[2][1] - axes[1][1] * axes[2][0])
        if matrix3(axes)
        else 0.0
    )
    reconstructed_com = (
        tuple(
            tuple(
                sum(axes[index][row] * moments[index] * axes[index][column] for index in range(3))
                for column in range(3)
            )
            for row in range(3)
        )
        if matrix3(axes) and len(moments) == 3
        else ()
    )
    displacement = (
        tuple(center[index] - origin[index] for index in range(3))
        if len(center) == len(origin) == 3
        else ()
    )
    displacement_sq = sum(item * item for item in displacement)
    parallel_axis = (
        tuple(
            tuple(
                inertia_com[row_index][column_index]
                + mass
                * (
                    displacement_sq * (1.0 if row_index == column_index else 0.0)
                    - displacement[row_index] * displacement[column_index]
                )
                for column_index in range(3)
            )
            for row_index in range(3)
        )
        if matrix3(inertia_com) and len(displacement) == 3
        else ()
    )
    shape_digest = str(value.get("mass_shape_sha256", "")).lower()
    if (
        not generation
        or any(
            value.get(key) != generation
            for key in (
                "density_generation",
                "origin_generation",
                "center_generation",
                "principal_generation",
                "axis_generation",
                "parallel_axis_generation",
                "owner_generation",
                "shape_generation",
                "result_generation",
            )
        )
        or not math.isfinite(density)
        or density <= 0.0
        or result_density != density
        or value.get("density_unit") != "kg/m^3"
        or value.get("result_density_unit") != value.get("density_unit")
        or not math.isfinite(mass)
        or mass <= 0.0
        or result_mass != mass
        or len(origin) != 3
        or any(not math.isfinite(item) for item in origin)
        or result_origin != origin
        or len(center) != 3
        or any(not math.isfinite(item) for item in center)
        or result_center != center
        or not matrix3(inertia_com)
        or not matrix3(inertia_origin)
        or result_inertia_com != inertia_com
        or result_inertia_origin != inertia_origin
        or len(moments) != 3
        or any(not math.isfinite(item) or item <= 0.0 for item in moments)
        or tuple(sorted(moments)) != moments
        or result_moments != moments
        or not axes_orthonormal
        or not close(axes_determinant, 1.0)
        or result_axes != axes
        or not matrix3(reconstructed_com)
        or any(
            not close(reconstructed_com[row_index][column_index], inertia_com[row_index][column_index])
            for row_index in range(3)
            for column_index in range(3)
        )
        or not matrix3(parallel_axis)
        or any(
            not close(parallel_axis[row_index][column_index], inertia_origin[row_index][column_index])
            for row_index in range(3)
            for column_index in range(3)
        )
        or value.get("degeneracy_convention") != "right_handed_sorted_moments"
        or value.get("result_degeneracy_convention") != value.get("degeneracy_convention")
        or not str(value.get("shape_owner", "")).strip()
        or value.get("result_shape_owner") != value.get("shape_owner")
        or not _valid_identity_digest(shape_digest)
        or value.get("result_mass_shape_sha256") != shape_digest
    ):
        return None
    return (
        generation,
        density,
        mass,
        origin,
        center,
        inertia_com,
        inertia_origin,
        moments,
        axes,
        value.get("shape_owner"),
        shape_digest,
    )


def _assembly_mate_mass_inertia_generation_identity(row):
    value = row.get(
        "assembly_mate_transform_cycle_frame_handedness_mass_center_inertia_owner_shape_result_generation_identity"
    )
    if not isinstance(value, dict):
        return None

    def matrix(raw, size):
        parsed = tuple(tuple(float(item) for item in matrix_row) for matrix_row in raw)
        return parsed if len(parsed) == size and all(len(matrix_row) == size for matrix_row in parsed) else ()

    def multiply(left, right):
        return tuple(
            tuple(
                sum(left[row][inner] * right[inner][column] for inner in range(3))
                for column in range(3)
            )
            for row in range(3)
        )

    def transpose(raw):
        return tuple(tuple(raw[column][row] for column in range(3)) for row in range(3))

    def determinant(raw):
        return (
            raw[0][0] * (raw[1][1] * raw[2][2] - raw[1][2] * raw[2][1])
            - raw[0][1] * (raw[1][0] * raw[2][2] - raw[1][2] * raw[2][0])
            + raw[0][2] * (raw[1][0] * raw[2][1] - raw[1][1] * raw[2][0])
        )

    def close(left, right):
        return math.isclose(left, right, rel_tol=1.0e-9, abs_tol=1.0e-12)

    generation = str(value.get("assembly_generation", "")).strip()
    try:
        part_ids = tuple(str(item).strip() for item in value.get("part_ids", []))
        result_part_ids = tuple(str(item).strip() for item in value.get("result_part_ids", []))
        mate_edges = tuple(tuple(str(item).strip() for item in pair) for pair in value.get("mate_edges", []))
        result_mate_edges = tuple(
            tuple(str(item).strip() for item in pair)
            for pair in value.get("result_mate_edges", [])
        )
        cycle = matrix(value.get("mate_cycle_transform", []), 4)
        result_cycle = matrix(value.get("result_mate_cycle_transform", []), 4)
        frame_determinants = tuple(float(item) for item in value.get("frame_determinants", []))
        result_frame_determinants = tuple(
            float(item) for item in value.get("result_frame_determinants", [])
        )
        masses = tuple(float(item) for item in value.get("part_masses_kg", []))
        result_masses = tuple(float(item) for item in value.get("result_part_masses_kg", []))
        centers = tuple(tuple(float(item) for item in row) for row in value.get("part_centers_m", []))
        result_centers = tuple(
            tuple(float(item) for item in row)
            for row in value.get("result_part_centers_m", [])
        )
        assembly_mass = float(value.get("assembly_mass_kg"))
        result_assembly_mass = float(value.get("result_assembly_mass_kg"))
        assembly_center = tuple(float(item) for item in value.get("assembly_center_of_mass_m", []))
        result_assembly_center = tuple(
            float(item) for item in value.get("result_assembly_center_of_mass_m", [])
        )
        rotations = tuple(matrix(raw, 3) for raw in value.get("part_rotation_matrices", []))
        result_rotations = tuple(
            matrix(raw, 3) for raw in value.get("result_part_rotation_matrices", [])
        )
        inertia_local = tuple(matrix(raw, 3) for raw in value.get("part_inertia_local_kg_m2", []))
        result_inertia_local = tuple(
            matrix(raw, 3) for raw in value.get("result_part_inertia_local_kg_m2", [])
        )
        inertia_global = tuple(matrix(raw, 3) for raw in value.get("part_inertia_global_kg_m2", []))
        result_inertia_global = tuple(
            matrix(raw, 3) for raw in value.get("result_part_inertia_global_kg_m2", [])
        )
    except (TypeError, ValueError):
        return None
    part_count = len(part_ids)
    expected_identity = tuple(
        tuple(1.0 if row == column else 0.0 for column in range(4))
        for row in range(4)
    )
    edge_sources = [pair[0] for pair in mate_edges if len(pair) == 2]
    edge_targets = [pair[1] for pair in mate_edges if len(pair) == 2]
    weighted_center = (
        tuple(
            sum(mass * center[axis] for mass, center in zip(masses, centers))
            / assembly_mass
            for axis in range(3)
        )
        if assembly_mass > 0.0 and len(centers) == part_count and all(len(center) == 3 for center in centers)
        else ()
    )
    inertia_rotations_ok = len(rotations) == len(inertia_local) == len(inertia_global) == part_count
    if inertia_rotations_ok:
        for rotation, local, global_tensor in zip(rotations, inertia_local, inertia_global):
            if not rotation or not local or not global_tensor:
                inertia_rotations_ok = False
                break
            orthogonality = multiply(rotation, transpose(rotation))
            rotated = multiply(multiply(rotation, local), transpose(rotation))
            if (
                not close(determinant(rotation), 1.0)
                or any(
                    not close(orthogonality[row][column], 1.0 if row == column else 0.0)
                    for row in range(3)
                    for column in range(3)
                )
                or any(
                    not close(rotated[row][column], global_tensor[row][column])
                    for row in range(3)
                    for column in range(3)
                )
            ):
                inertia_rotations_ok = False
                break
    shape_digest = str(value.get("assembly_shape_sha256", "")).lower()
    if (
        not generation
        or any(
            value.get(key) != generation
            for key in (
                "mate_generation",
                "cycle_generation",
                "frame_generation",
                "mass_generation",
                "center_generation",
                "inertia_generation",
                "owner_generation",
                "shape_generation",
                "result_generation",
            )
        )
        or part_count < 2
        or len(set(part_ids)) != part_count
        or not all(part_ids)
        or result_part_ids != part_ids
        or len(mate_edges) != part_count
        or any(len(pair) != 2 or pair[0] == pair[1] for pair in mate_edges)
        or set(edge_sources) != set(part_ids)
        or set(edge_targets) != set(part_ids)
        or result_mate_edges != mate_edges
        or cycle != expected_identity
        or result_cycle != cycle
        or len(frame_determinants) != part_count
        or any(not close(item, 1.0) for item in frame_determinants)
        or result_frame_determinants != frame_determinants
        or len(masses) != part_count
        or any(not math.isfinite(item) or item <= 0.0 for item in masses)
        or result_masses != masses
        or len(centers) != part_count
        or any(len(center) != 3 or any(not math.isfinite(item) for item in center) for center in centers)
        or result_centers != centers
        or not close(assembly_mass, sum(masses))
        or result_assembly_mass != assembly_mass
        or len(assembly_center) != 3
        or len(weighted_center) != 3
        or any(not close(item, expected) for item, expected in zip(assembly_center, weighted_center))
        or result_assembly_center != assembly_center
        or not inertia_rotations_ok
        or result_rotations != rotations
        or result_inertia_local != inertia_local
        or result_inertia_global != inertia_global
        or not str(value.get("assembly_owner", "")).strip()
        or value.get("result_assembly_owner") != value.get("assembly_owner")
        or not _valid_identity_digest(shape_digest)
        or value.get("accepted_assembly_shape_sha256") != shape_digest
    ):
        return None
    return (
        generation,
        part_ids,
        mate_edges,
        cycle,
        masses,
        centers,
        assembly_mass,
        assembly_center,
        rotations,
        inertia_local,
        inertia_global,
        value.get("assembly_owner"),
        shape_digest,
    )


def _shell_fillet_topology_generation_identity(row):
    value = row.get(
        "shell_fillet_topology_euler_manifold_thickness_volume_area_inertia_convergence_brep_result_generation_identity"
    )
    if not isinstance(value, dict):
        return None
    generation = str(value.get("shell_fillet_generation", "")).strip()
    try:
        vertices = int(value.get("vertex_count"))
        edges = int(value.get("edge_count"))
        faces = int(value.get("face_count"))
        euler = int(value.get("euler_characteristic"))
        incidence = tuple(int(item) for item in value.get("edge_face_incidence_counts", []))
        thickness = float(value.get("nominal_wall_thickness_m"))
        thickness_samples = tuple(float(item) for item in value.get("wall_thickness_samples_m", []))
        original_volume = float(value.get("original_volume_m3"))
        removed_volume = float(value.get("removed_volume_m3"))
        shell_volume = float(value.get("shell_volume_m3"))
        area = float(value.get("surface_area_m2"))
        inertia = tuple(
            tuple(float(item) for item in matrix_row)
            for matrix_row in value.get("inertia_tensor_kg_m2", [])
        )
        tolerances = tuple(float(item) for item in value.get("convergence_tolerances_m", []))
        convergence_volumes = tuple(float(item) for item in value.get("convergence_volumes_m3", []))
    except (TypeError, ValueError):
        return None
    mirrored_fields = (
        "vertex_count",
        "edge_count",
        "face_count",
        "euler_characteristic",
        "edge_face_incidence_counts",
        "nominal_wall_thickness_m",
        "wall_thickness_samples_m",
        "original_volume_m3",
        "removed_volume_m3",
        "shell_volume_m3",
        "surface_area_m2",
        "inertia_tensor_kg_m2",
        "convergence_tolerances_m",
        "convergence_volumes_m3",
    )
    matrix3 = len(inertia) == 3 and all(len(matrix_row) == 3 for matrix_row in inertia)
    symmetric = matrix3 and all(
        math.isclose(inertia[row][column], inertia[column][row], rel_tol=1.0e-9, abs_tol=1.0e-12)
        for row in range(3)
        for column in range(3)
    )
    positive_definite = (
        symmetric
        and inertia[0][0] > 0.0
        and inertia[0][0] * inertia[1][1] - inertia[0][1] ** 2 > 0.0
        and (
            inertia[0][0] * (inertia[1][1] * inertia[2][2] - inertia[1][2] * inertia[2][1])
            - inertia[0][1] * (inertia[1][0] * inertia[2][2] - inertia[1][2] * inertia[2][0])
            + inertia[0][2] * (inertia[1][0] * inertia[2][1] - inertia[1][1] * inertia[2][0])
        )
        > 0.0
    )
    shape_digest = str(value.get("shell_brep_sha256", "")).lower()
    if (
        not generation
        or any(
            value.get(key) != generation
            for key in (
                "topology_generation",
                "thickness_generation",
                "volume_generation",
                "area_generation",
                "inertia_generation",
                "convergence_generation",
                "brep_generation",
                "result_generation",
            )
        )
        or min(vertices, edges, faces) <= 0
        or euler != 2
        or vertices - edges + faces != euler
        or len(incidence) != edges
        or any(item != 2 for item in incidence)
        or not math.isfinite(thickness)
        or thickness <= 0.0
        or not thickness_samples
        or any(
            not math.isfinite(item)
            or not math.isclose(item, thickness, rel_tol=1.0e-3, abs_tol=1.0e-12)
            for item in thickness_samples
        )
        or any(not math.isfinite(item) or item <= 0.0 for item in (original_volume, removed_volume, shell_volume, area))
        or removed_volume >= original_volume
        or not math.isclose(shell_volume, original_volume - removed_volume, rel_tol=1.0e-9, abs_tol=1.0e-12)
        or not positive_definite
        or len(tolerances) < 3
        or len(convergence_volumes) != len(tolerances)
        or any(item <= 0.0 or not math.isfinite(item) for item in tolerances + convergence_volumes)
        or any(right >= left for left, right in zip(tolerances, tolerances[1:]))
        or not math.isclose(convergence_volumes[-1], shell_volume, rel_tol=1.0e-9, abs_tol=1.0e-12)
        or not math.isclose(convergence_volumes[-2], convergence_volumes[-1], rel_tol=1.0e-5, abs_tol=1.0e-12)
        or any(value.get(f"result_{field}") != value.get(field) for field in mirrored_fields)
        or not _valid_identity_digest(shape_digest)
        or value.get("accepted_shell_brep_sha256") != shape_digest
    ):
        return None
    return (
        generation,
        vertices,
        edges,
        faces,
        euler,
        incidence,
        thickness,
        thickness_samples,
        original_volume,
        removed_volume,
        shell_volume,
        area,
        inertia,
        tolerances,
        convergence_volumes,
        shape_digest,
    )


def shape_mass_property_crosscheck_summary(
    reference_rows,
    measured_sets,
    rtol=1.0e-5,
    bbox_atol=1.0e-6,
    worst_limit=5,
):
    """Return volume/area/bbox crosscheck summary for one or more CAD sources.

    This is the stronger companion to :func:`shape_volume_crosscheck_summary`.
    Use it when the external CAD kernel can report at least ``name``,
    ``volume`` and ``area`` rows, and optionally a ``bounding_box`` compatible
    with :func:`shape_measurement_row`.  It is intended for build123d -> STEP
    round trips through Cubit, CST, or another CAD kernel before meshing or
    solver setup.
    """

    if worst_limit < 0:
        raise ValueError("worst_limit must be non-negative")
    reference = list(reference_rows)
    normalized_sets = [
        (label, list(rows))
        for label, rows in _normalize_shape_measurement_sets(measured_sets)
    ]
    identity_rows = [*reference]
    for _, rows in normalized_sets:
        identity_rows.extend(rows)

    assembly_location_evidence_present = any(
        row.get(
            "assembly_location_boolean_operand_revision_mass_property_generation_identity"
        )
        is not None
        for row in identity_rows
    )
    reference_assembly_locations = {
        str(row.get("name", "")): _assembly_location_boolean_mass_identity(row)
        for row in reference
    }
    assembly_location_identity_ok = not assembly_location_evidence_present
    if assembly_location_evidence_present:
        assembly_location_identity_ok = bool(reference_assembly_locations) and all(
            value is not None for value in reference_assembly_locations.values()
        )
        for _, rows in normalized_sets:
            assembly_location_identity_ok = assembly_location_identity_ok and all(
                _assembly_location_boolean_mass_identity(row)
                == reference_assembly_locations.get(str(row.get("name", "")))
                for row in rows
            )

    loft_volume_evidence_present = any(
        row.get("loft_spline_tessellation_watertight_volume_generation_identity")
        is not None
        for row in identity_rows
    )
    reference_loft_volumes = {
        str(row.get("name", "")): _loft_spline_watertight_volume_identity(row)
        for row in reference
    }
    loft_volume_identity_ok = not loft_volume_evidence_present
    if loft_volume_evidence_present:
        loft_volume_identity_ok = bool(reference_loft_volumes) and all(
            value is not None for value in reference_loft_volumes.values()
        )
        for _, rows in normalized_sets:
            loft_volume_identity_ok = loft_volume_identity_ok and all(
                _loft_spline_watertight_volume_identity(row)
                == reference_loft_volumes.get(str(row.get("name", "")))
                for row in rows
            )

    transformed_assembly_evidence_present = any(
        row.get(
            "transformed_assembly_com_inertia_axis_density_unit_generation_identity"
        )
        is not None
        for row in identity_rows
    )
    reference_transformed_assemblies = {
        str(row.get("name", "")): _transformed_assembly_mass_identity(row)
        for row in reference
    }
    transformed_assembly_identity_ok = not transformed_assembly_evidence_present
    if transformed_assembly_evidence_present:
        transformed_assembly_identity_ok = bool(
            reference_transformed_assemblies
        ) and all(
            value is not None for value in reference_transformed_assemblies.values()
        )
        for _, rows in normalized_sets:
            transformed_assembly_identity_ok = (
                transformed_assembly_identity_ok
                and all(
                    _transformed_assembly_mass_identity(row)
                    == reference_transformed_assemblies.get(
                        str(row.get("name", ""))
                    )
                    for row in rows
                )
            )

    fillet_chamfer_evidence_present = any(
        row.get(
            "fillet_chamfer_topology_naming_edge_selection_fingerprint_identity"
        )
        is not None
        for row in identity_rows
    )
    reference_fillet_chamfer = {
        str(row.get("name", "")): _fillet_chamfer_topology_identity(row)
        for row in reference
    }
    fillet_chamfer_identity_ok = not fillet_chamfer_evidence_present
    if fillet_chamfer_evidence_present:
        fillet_chamfer_identity_ok = bool(reference_fillet_chamfer) and all(
            value is not None for value in reference_fillet_chamfer.values()
        )
        for _, rows in normalized_sets:
            fillet_chamfer_identity_ok = fillet_chamfer_identity_ok and all(
                _fillet_chamfer_topology_identity(row)
                == reference_fillet_chamfer.get(str(row.get("name", "")))
                for row in rows
            )

    boolean_operation_evidence_present = any(
        row.get("boolean_tolerance_healing_topology_volume_generation_identity")
        is not None
        for row in identity_rows
    )
    reference_boolean_operations = {
        str(row.get("name", "")): _boolean_tolerance_healing_topology_volume_identity(
            row
        )
        for row in reference
    }
    boolean_operation_identity_ok = not boolean_operation_evidence_present
    if boolean_operation_evidence_present:
        boolean_operation_identity_ok = bool(reference_boolean_operations) and all(
            value is not None for value in reference_boolean_operations.values()
        )
        for _, rows in normalized_sets:
            boolean_operation_identity_ok = boolean_operation_identity_ok and all(
                _boolean_tolerance_healing_topology_volume_identity(row)
                == reference_boolean_operations.get(str(row.get("name", "")))
                for row in rows
            )

    assembly_kinematics_evidence_present = any(
        row.get("assembly_mate_transform_dof_loop_closure_generation_identity")
        is not None
        for row in identity_rows
    )
    reference_assembly_kinematics = {
        str(row.get("name", "")): _assembly_mate_transform_dof_loop_closure_identity(
            row
        )
        for row in reference
    }
    assembly_kinematics_identity_ok = not assembly_kinematics_evidence_present
    if assembly_kinematics_evidence_present:
        assembly_kinematics_identity_ok = bool(reference_assembly_kinematics) and all(
            value is not None for value in reference_assembly_kinematics.values()
        )
        for _, rows in normalized_sets:
            assembly_kinematics_identity_ok = assembly_kinematics_identity_ok and all(
                _assembly_mate_transform_dof_loop_closure_identity(row)
                == reference_assembly_kinematics.get(str(row.get("name", "")))
                for row in rows
            )

    feature_selection_evidence_present = any(
        row.get("fillet_chamfer_edge_selector_topology_naming_tolerance_shape_generation_identity")
        is not None for row in identity_rows
    )
    reference_feature_selection = {
        str(row.get("name", "")): _fillet_chamfer_selector_generation_identity(row)
        for row in reference
    }
    feature_selection_identity_ok = not feature_selection_evidence_present
    if feature_selection_evidence_present:
        feature_selection_identity_ok = bool(reference_feature_selection) and all(
            value is not None for value in reference_feature_selection.values()
        )
        for _, rows in normalized_sets:
            feature_selection_identity_ok = feature_selection_identity_ok and all(
                _fillet_chamfer_selector_generation_identity(row)
                == reference_feature_selection.get(str(row.get("name", "")))
                for row in rows
            )

    mass_frame_evidence_present = any(
        row.get("mass_density_center_inertia_reference_frame_assembly_generation_identity")
        is not None for row in identity_rows
    )
    reference_mass_frame = {
        str(row.get("name", "")): _mass_density_frame_generation_identity(row)
        for row in reference
    }
    mass_frame_identity_ok = not mass_frame_evidence_present
    if mass_frame_evidence_present:
        mass_frame_identity_ok = bool(reference_mass_frame) and all(
            value is not None for value in reference_mass_frame.values()
        )
        for _, rows in normalized_sets:
            mass_frame_identity_ok = mass_frame_identity_ok and all(
                _mass_density_frame_generation_identity(row)
                == reference_mass_frame.get(str(row.get("name", "")))
                for row in rows
            )

    boolean_imprint_evidence_present = any(
        row.get("boolean_imprint_interface_owner_topology_name_tolerance_mass_generation_identity")
        is not None for row in identity_rows
    )
    reference_boolean_imprints = {
        str(row.get("name", "")): _boolean_imprint_generation_identity(row)
        for row in reference
    }
    boolean_imprint_identity_ok = not boolean_imprint_evidence_present
    if boolean_imprint_evidence_present:
        boolean_imprint_identity_ok = bool(reference_boolean_imprints) and all(
            value is not None for value in reference_boolean_imprints.values()
        )
        for _, rows in normalized_sets:
            boolean_imprint_identity_ok = boolean_imprint_identity_ok and all(
                _boolean_imprint_generation_identity(row)
                == reference_boolean_imprints.get(str(row.get("name", "")))
                for row in rows
            )

    loft_section_evidence_present = any(
        row.get("loft_section_wire_seam_continuity_solid_volume_generation_identity")
        is not None for row in identity_rows
    )
    reference_loft_sections = {
        str(row.get("name", "")): _loft_section_generation_identity(row)
        for row in reference
    }
    loft_section_identity_ok = not loft_section_evidence_present
    if loft_section_evidence_present:
        loft_section_identity_ok = bool(reference_loft_sections) and all(
            value is not None for value in reference_loft_sections.values()
        )
        for _, rows in normalized_sets:
            loft_section_identity_ok = loft_section_identity_ok and all(
                _loft_section_generation_identity(row)
                == reference_loft_sections.get(str(row.get("name", "")))
                for row in rows
            )

    shell_offset_evidence_present = any(
        row.get("shell_offset_face_normal_thickness_join_self_intersection_mass_generation_identity")
        is not None for row in identity_rows
    )
    reference_shell_offsets = {
        str(row.get("name", "")): _shell_offset_generation_identity(row)
        for row in reference
    }
    shell_offset_identity_ok = not shell_offset_evidence_present
    if shell_offset_evidence_present:
        shell_offset_identity_ok = bool(reference_shell_offsets) and all(
            value is not None for value in reference_shell_offsets.values()
        )
        for _, rows in normalized_sets:
            shell_offset_identity_ok = shell_offset_identity_ok and all(
                _shell_offset_generation_identity(row)
                == reference_shell_offsets.get(str(row.get("name", "")))
                for row in rows
            )

    path_sweep_evidence_present = any(
        row.get("path_sweep_frame_transition_profile_orientation_solid_volume_generation_identity")
        is not None for row in identity_rows
    )
    reference_path_sweeps = {
        str(row.get("name", "")): _path_sweep_generation_identity(row)
        for row in reference
    }
    path_sweep_identity_ok = not path_sweep_evidence_present
    if path_sweep_evidence_present:
        path_sweep_identity_ok = bool(reference_path_sweeps) and all(
            value is not None for value in reference_path_sweeps.values()
        )
        for _, rows in normalized_sets:
            path_sweep_identity_ok = path_sweep_identity_ok and all(
                _path_sweep_generation_identity(row)
                == reference_path_sweeps.get(str(row.get("name", "")))
                for row in rows
            )

    sheet_metal_evidence_present = any(
        row.get(
            "sheet_metal_bend_allowance_kfactor_neutral_axis_relief_thickness_flat_pattern_area_generation_identity"
        ) is not None
        for row in identity_rows
    )
    reference_sheet_metal = {
        str(row.get("name", "")): _sheet_metal_flat_pattern_generation_identity(row)
        for row in reference
    }
    sheet_metal_identity_ok = not sheet_metal_evidence_present
    if sheet_metal_evidence_present:
        sheet_metal_identity_ok = bool(reference_sheet_metal) and all(
            value is not None for value in reference_sheet_metal.values()
        )
        for _, rows in normalized_sets:
            sheet_metal_identity_ok = sheet_metal_identity_ok and all(
                _sheet_metal_flat_pattern_generation_identity(row)
                == reference_sheet_metal.get(str(row.get("name", "")))
                for row in rows
            )

    joint_loop_evidence_present = any(
        row.get(
            "joint_kinematic_loop_graph_dof_limit_connector_frame_closure_configuration_swept_volume_generation_identity"
        ) is not None
        for row in identity_rows
    )
    reference_joint_loops = {
        str(row.get("name", "")): _joint_kinematic_loop_generation_identity(row)
        for row in reference
    }
    joint_loop_identity_ok = not joint_loop_evidence_present
    if joint_loop_evidence_present:
        joint_loop_identity_ok = bool(reference_joint_loops) and all(
            value is not None for value in reference_joint_loops.values()
        )
        for _, rows in normalized_sets:
            joint_loop_identity_ok = joint_loop_identity_ok and all(
                _joint_kinematic_loop_generation_identity(row)
                == reference_joint_loops.get(str(row.get("name", "")))
                for row in rows
            )

    helical_sweep_evidence_present = any(
        row.get(
            "helical_sweep_pitch_handedness_profile_frame_turn_self_intersection_volume_centroid_shape_generation_identity"
        ) is not None
        for row in identity_rows
    )
    reference_helical_sweeps = {
        str(row.get("name", "")): _helical_sweep_generation_identity(row)
        for row in reference
    }
    helical_sweep_identity_ok = not helical_sweep_evidence_present
    if helical_sweep_evidence_present:
        helical_sweep_identity_ok = bool(reference_helical_sweeps) and all(
            value is not None for value in reference_helical_sweeps.values()
        )
        for _, rows in normalized_sets:
            helical_sweep_identity_ok = helical_sweep_identity_ok and all(
                _helical_sweep_generation_identity(row)
                == reference_helical_sweeps.get(str(row.get("name", "")))
                for row in rows
            )

    boolean_history_evidence_present = any(
        row.get(
            "boolean_tolerance_operand_order_history_volume_centroid_inertia_shape_generation_identity"
        ) is not None
        for row in identity_rows
    )
    reference_boolean_histories = {
        str(row.get("name", "")): _boolean_history_generation_identity(row)
        for row in reference
    }
    boolean_history_identity_ok = not boolean_history_evidence_present
    if boolean_history_evidence_present:
        boolean_history_identity_ok = bool(reference_boolean_histories) and all(
            value is not None for value in reference_boolean_histories.values()
        )
        for _, rows in normalized_sets:
            boolean_history_identity_ok = boolean_history_identity_ok and all(
                _boolean_history_generation_identity(row)
                == reference_boolean_histories.get(str(row.get("name", "")))
                for row in rows
            )

    assembly_occurrence_evidence_present = any(
        row.get(
            "assembly_occurrence_location_density_unit_suppression_mass_center_inertia_parallel_axis_shape_generation_identity"
        ) is not None
        for row in identity_rows
    )
    reference_assembly_occurrences = {
        str(row.get("name", "")): _assembly_occurrence_mass_generation_identity(row)
        for row in reference
    }
    assembly_occurrence_identity_ok = not assembly_occurrence_evidence_present
    if assembly_occurrence_evidence_present:
        assembly_occurrence_identity_ok = bool(reference_assembly_occurrences) and all(
            value is not None for value in reference_assembly_occurrences.values()
        )
        for _, rows in normalized_sets:
            assembly_occurrence_identity_ok = assembly_occurrence_identity_ok and all(
                _assembly_occurrence_mass_generation_identity(row)
                == reference_assembly_occurrences.get(str(row.get("name", "")))
                for row in rows
            )

    loft_lineage_evidence_present = any(
        row.get(
            "loft_sweep_profile_order_seam_guide_orientation_face_lineage_shell_volume_shape_generation_identity"
        ) is not None
        for row in identity_rows
    )
    reference_loft_lineages = {
        str(row.get("name", "")): _loft_face_lineage_generation_identity(row)
        for row in reference
    }
    loft_lineage_identity_ok = not loft_lineage_evidence_present
    if loft_lineage_evidence_present:
        loft_lineage_identity_ok = bool(reference_loft_lineages) and all(
            value is not None for value in reference_loft_lineages.values()
        )
        for _, rows in normalized_sets:
            loft_lineage_identity_ok = loft_lineage_identity_ok and all(
                _loft_face_lineage_generation_identity(row)
                == reference_loft_lineages.get(str(row.get("name", "")))
                for row in rows
            )

    boolean_topology_evidence_present = any(
        row.get(
            "boolean_fuzzy_tolerance_topology_name_face_ancestry_count_volume_centroid_shape_generation_identity"
        ) is not None
        for row in identity_rows
    )
    reference_boolean_topologies = {
        str(row.get("name", "")): _boolean_topology_ancestry_generation_identity(row)
        for row in reference
    }
    boolean_topology_identity_ok = not boolean_topology_evidence_present
    if boolean_topology_evidence_present:
        boolean_topology_identity_ok = bool(reference_boolean_topologies) and all(
            value is not None for value in reference_boolean_topologies.values()
        )
        for _, rows in normalized_sets:
            boolean_topology_identity_ok = boolean_topology_identity_ok and all(
                _boolean_topology_ancestry_generation_identity(row)
                == reference_boolean_topologies.get(str(row.get("name", "")))
                for row in rows
            )

    sweep_frame_transition_evidence_present = any(
        row.get(
            "sweep_frenet_frame_twist_transition_profile_orientation_self_intersection_volume_owner_shape_generation_identity"
        ) is not None
        for row in identity_rows
    )
    reference_sweep_frame_transitions = {
        str(row.get("name", "")): _sweep_frame_transition_generation_identity(row)
        for row in reference
    }
    sweep_frame_transition_identity_ok = not sweep_frame_transition_evidence_present
    if sweep_frame_transition_evidence_present:
        sweep_frame_transition_identity_ok = bool(reference_sweep_frame_transitions) and all(
            value is not None for value in reference_sweep_frame_transitions.values()
        )
        for _, rows in normalized_sets:
            sweep_frame_transition_identity_ok = (
                sweep_frame_transition_identity_ok
                and all(
                    _sweep_frame_transition_generation_identity(row)
                    == reference_sweep_frame_transitions.get(str(row.get("name", "")))
                    for row in rows
                )
            )

    guided_loft_evidence_present = any(
        row.get(
            "loft_section_guide_parameterization_seam_orientation_mode_intersection_volume_shape_generation_identity"
        )
        is not None
        for row in identity_rows
    )
    reference_guided_lofts = {
        str(row.get("name", "")): _guided_loft_generation_identity(row)
        for row in reference
    }
    guided_loft_identity_ok = not guided_loft_evidence_present
    if guided_loft_evidence_present:
        guided_loft_identity_ok = bool(reference_guided_lofts) and all(
            value is not None for value in reference_guided_lofts.values()
        )
        for _, rows in normalized_sets:
            guided_loft_identity_ok = guided_loft_identity_ok and all(
                _guided_loft_generation_identity(row)
                == reference_guided_lofts.get(str(row.get("name", "")))
                for row in rows
            )

    mass_inertia_evidence_present = any(
        row.get(
            "mass_property_density_unit_origin_center_principal_axis_degeneracy_parallel_axis_owner_shape_generation_identity"
        )
        is not None
        for row in identity_rows
    )
    reference_mass_inertias = {
        str(row.get("name", "")): _mass_inertia_parallel_axis_generation_identity(row)
        for row in reference
    }
    mass_inertia_identity_ok = not mass_inertia_evidence_present
    if mass_inertia_evidence_present:
        mass_inertia_identity_ok = bool(reference_mass_inertias) and all(
            value is not None for value in reference_mass_inertias.values()
        )
        for _, rows in normalized_sets:
            mass_inertia_identity_ok = mass_inertia_identity_ok and all(
                _mass_inertia_parallel_axis_generation_identity(row)
                == reference_mass_inertias.get(str(row.get("name", "")))
                for row in rows
            )

    assembly_mate_evidence_present = any(
        row.get(
            "assembly_mate_transform_cycle_frame_handedness_mass_center_inertia_owner_shape_result_generation_identity"
        )
        is not None
        for row in identity_rows
    )
    reference_assembly_mates = {
        str(row.get("name", "")): _assembly_mate_mass_inertia_generation_identity(row)
        for row in reference
    }
    assembly_mate_mass_inertia_identity_ok = not assembly_mate_evidence_present
    if assembly_mate_evidence_present:
        assembly_mate_mass_inertia_identity_ok = bool(reference_assembly_mates) and all(
            value is not None for value in reference_assembly_mates.values()
        )
        for _, rows in normalized_sets:
            assembly_mate_mass_inertia_identity_ok = (
                assembly_mate_mass_inertia_identity_ok
                and all(
                    _assembly_mate_mass_inertia_generation_identity(row)
                    == reference_assembly_mates.get(str(row.get("name", "")))
                    for row in rows
                )
            )

    shell_fillet_evidence_present = any(
        row.get(
            "shell_fillet_topology_euler_manifold_thickness_volume_area_inertia_convergence_brep_result_generation_identity"
        )
        is not None
        for row in identity_rows
    )
    reference_shell_fillets = {
        str(row.get("name", "")): _shell_fillet_topology_generation_identity(row)
        for row in reference
    }
    shell_fillet_identity_ok = not shell_fillet_evidence_present
    if shell_fillet_evidence_present:
        shell_fillet_identity_ok = bool(reference_shell_fillets) and all(
            value is not None for value in reference_shell_fillets.values()
        )
        for _, rows in normalized_sets:
            shell_fillet_identity_ok = shell_fillet_identity_ok and all(
                _shell_fillet_topology_generation_identity(row)
                == reference_shell_fillets.get(str(row.get("name", "")))
                for row in rows
            )

    revision_evidence_present = any(
        row.get("brep_identity") is not None
        or row.get("mass_property_identity") is not None
        for row in identity_rows
    )

    def mass_property_revision_ok(row):
        brep = row.get("brep_identity")
        measured = row.get("mass_property_identity")
        if not isinstance(brep, dict) or not isinstance(measured, dict):
            return False
        revision = str(brep.get("revision", "")).strip()
        digest = str(brep.get("sha256", "")).strip().lower()
        return (
            bool(revision)
            and len(digest) == 64
            and all(character in "0123456789abcdef" for character in digest)
            and measured.get("brep_revision") == revision
            and str(measured.get("brep_sha256", "")).lower() == digest
        )

    revision_identity_ok = not revision_evidence_present or all(
        mass_property_revision_ok(row) for row in identity_rows
    )

    assembly_evidence_present = any(
        row.get("assembly_identity") is not None for row in identity_rows
    )

    def assembly_identity(row):
        value = row.get("assembly_identity")
        if not isinstance(value, dict):
            return None
        generation = str(value.get("generation", "")).strip()
        children = value.get("child_revisions")
        if not generation or not isinstance(children, dict) or not children:
            return None
        normalized = {str(name): str(revision) for name, revision in children.items()}
        if not all(normalized.values()):
            return None
        return generation, normalized

    reference_assembly = {
        str(row.get("name", "")): assembly_identity(row) for row in reference
    }
    assembly_identity_ok = not assembly_evidence_present
    if assembly_evidence_present:
        reference_values = list(reference_assembly.values())
        assembly_identity_ok = (
            bool(reference_values)
            and all(value is not None for value in reference_values)
            and len({repr(value) for value in reference_values}) == 1
        )
        for _, rows in normalized_sets:
            measured_values = [assembly_identity(row) for row in rows]
            assembly_identity_ok = assembly_identity_ok and (
                bool(measured_values)
                and all(value is not None for value in measured_values)
                and len({repr(value) for value in measured_values}) == 1
                and all(
                    assembly_identity(row)
                    == reference_assembly.get(str(row.get("name", "")))
                    for row in rows
                )
            )

    frame_evidence_present = any(
        row.get("mass_property_frame_identity") is not None for row in identity_rows
    )

    def frame_identity(row):
        value = row.get("mass_property_frame_identity")
        if not isinstance(value, dict):
            return None
        frame_id = str(value.get("frame_id", "")).strip()
        transform_generation = str(value.get("transform_generation", "")).strip()
        return (frame_id, transform_generation) if frame_id and transform_generation else None

    reference_frames = {
        str(row.get("name", "")): frame_identity(row) for row in reference
    }
    frame_identity_ok = not frame_evidence_present
    if frame_evidence_present:
        frame_identity_ok = bool(reference_frames) and all(
            value is not None for value in reference_frames.values()
        )
        for _, rows in normalized_sets:
            frame_identity_ok = frame_identity_ok and all(
                frame_identity(row)
                == reference_frames.get(str(row.get("name", "")))
                for row in rows
            )

    topology_evidence_present = any(
        row.get("topology_identity") is not None for row in identity_rows
    )

    def topology_identity(row):
        value = row.get("topology_identity")
        brep = row.get("brep_identity")
        if not isinstance(value, dict) or not isinstance(brep, dict):
            return None
        adjacency_digest = str(value.get("face_adjacency_sha256", "")).lower()
        identity = (
            str(value.get("brep_revision", "")),
            str(value.get("brep_sha256", "")).lower(),
            adjacency_digest,
        )
        if (
            identity[0] != str(brep.get("revision", ""))
            or identity[1] != str(brep.get("sha256", "")).lower()
            or len(adjacency_digest) != 64
        ):
            return None
        return identity

    reference_topology = {
        str(row.get("name", "")): topology_identity(row) for row in reference
    }
    topology_identity_ok = not topology_evidence_present
    if topology_evidence_present:
        topology_identity_ok = bool(reference_topology) and all(
            value is not None for value in reference_topology.values()
        )
        for _, rows in normalized_sets:
            topology_identity_ok = topology_identity_ok and all(
                topology_identity(row)
                == reference_topology.get(str(row.get("name", "")))
                for row in rows
            )

    compound_evidence_present = any(
        row.get("compound_volume_identity") is not None for row in identity_rows
    )

    def compound_volume_identity(row):
        value = row.get("compound_volume_identity")
        if not isinstance(value, dict):
            return None
        try:
            overlap_volume = float(value.get("overlap_volume"))
        except (TypeError, ValueError):
            return None
        generation = str(value.get("topology_generation", "")).strip()
        volume_generation = str(value.get("volume_generation", "")).strip()
        if (
            value.get("topology_kind") != "physical_union_solid"
            or value.get("reported_volume_basis") != "physical_union"
            or overlap_volume < 0.0
            or not generation
            or volume_generation != generation
        ):
            return None
        return generation, volume_generation

    reference_compounds = {
        str(row.get("name", "")): compound_volume_identity(row) for row in reference
    }
    compound_volume_identity_ok = not compound_evidence_present
    if compound_evidence_present:
        compound_volume_identity_ok = bool(reference_compounds) and all(
            value is not None for value in reference_compounds.values()
        )
        for _, rows in normalized_sets:
            compound_volume_identity_ok = compound_volume_identity_ok and all(
                compound_volume_identity(row)
                == reference_compounds.get(str(row.get("name", "")))
                for row in rows
            )

    placement_evidence_present = any(
        row.get("placement_transform_identity") is not None for row in identity_rows
    )

    def placement_transform_identity(row):
        value = row.get("placement_transform_identity")
        if not isinstance(value, dict):
            return None
        frame = str(value.get("center_of_mass_frame", "")).strip()
        center_generation = str(
            value.get("center_of_mass_transform_generation", "")
        ).strip()
        final_generation = str(
            value.get("final_placement_transform_generation", "")
        ).strip()
        if not frame or not center_generation or center_generation != final_generation:
            return None
        return frame, final_generation

    reference_placements = {
        str(row.get("name", "")): placement_transform_identity(row)
        for row in reference
    }
    placement_transform_identity_ok = not placement_evidence_present
    if placement_evidence_present:
        placement_transform_identity_ok = bool(reference_placements) and all(
            value is not None for value in reference_placements.values()
        )
        for _, rows in normalized_sets:
            placement_transform_identity_ok = placement_transform_identity_ok and all(
                placement_transform_identity(row)
                == reference_placements.get(str(row.get("name", "")))
                for row in rows
            )

    healing_evidence_present = any(
        row.get("shape_healing_identity") is not None for row in identity_rows
    )

    def valid_sha256(value):
        digest = str(value or "").lower()
        return len(digest) == 64 and all(
            character in "0123456789abcdef" for character in digest
        )

    def shape_healing_identity(row):
        value = row.get("shape_healing_identity")
        if not isinstance(value, dict):
            return None
        pre_heal = str(value.get("pre_heal_brep_sha256", "")).lower()
        healed = str(value.get("healed_brep_sha256", "")).lower()
        final = str(value.get("final_brep_sha256", "")).lower()
        mass_property = str(value.get("mass_property_brep_sha256", "")).lower()
        generation = str(value.get("final_shape_generation", "")).strip()
        mass_generation = str(
            value.get("mass_property_shape_generation", "")
        ).strip()
        if (
            not all(valid_sha256(digest) for digest in (pre_heal, healed, final))
            or pre_heal == healed
            or healed != final
            or mass_property != final
            or not generation
            or mass_generation != generation
        ):
            return None
        return final, generation

    reference_healing = {
        str(row.get("name", "")): shape_healing_identity(row) for row in reference
    }
    shape_healing_identity_ok = not healing_evidence_present
    if healing_evidence_present:
        shape_healing_identity_ok = bool(reference_healing) and all(
            value is not None for value in reference_healing.values()
        )
        for _, rows in normalized_sets:
            shape_healing_identity_ok = shape_healing_identity_ok and all(
                shape_healing_identity(row)
                == reference_healing.get(str(row.get("name", "")))
                for row in rows
            )

    inertia_evidence_present = any(
        row.get("inertia_tensor_identity") is not None for row in identity_rows
    )

    def inertia_tensor_identity(row):
        value = row.get("inertia_tensor_identity")
        frame = frame_identity(row)
        placement = placement_transform_identity(row)
        if not isinstance(value, dict) or frame is None or placement is None:
            return None
        tensor_frame = str(value.get("tensor_frame_id", "")).strip()
        center_frame = str(value.get("center_of_mass_frame_id", "")).strip()
        tensor_generation = str(
            value.get("tensor_transform_generation", "")
        ).strip()
        final_generation = str(
            value.get("final_placement_transform_generation", "")
        ).strip()
        try:
            determinant = float(value.get("mirror_transform_determinant"))
        except (TypeError, ValueError):
            return None
        if (
            tensor_frame != center_frame
            or tensor_frame != frame[0]
            or tensor_generation != final_generation
            or tensor_generation != placement[1]
            or value.get("mirror_transform_applied") is not True
            or not math.isclose(determinant, -1.0, rel_tol=0.0, abs_tol=1.0e-12)
            or value.get("tensor_basis_handedness") != "right_handed"
        ):
            return None
        return tensor_frame, tensor_generation, determinant

    reference_inertia = {
        str(row.get("name", "")): inertia_tensor_identity(row) for row in reference
    }
    inertia_tensor_identity_ok = not inertia_evidence_present
    if inertia_evidence_present:
        inertia_tensor_identity_ok = bool(reference_inertia) and all(
            value is not None for value in reference_inertia.values()
        )
        for _, rows in normalized_sets:
            inertia_tensor_identity_ok = inertia_tensor_identity_ok and all(
                inertia_tensor_identity(row)
                == reference_inertia.get(str(row.get("name", "")))
                for row in rows
            )

    assembly_coordinate_evidence_present = any(
        row.get("assembly_mass_property_coordinate_identity") is not None
        for row in identity_rows
    )

    def assembly_mass_property_coordinate_identity(row):
        value = row.get("assembly_mass_property_coordinate_identity")
        if not isinstance(value, dict):
            return None
        assembly_generation = str(value.get("assembly_generation", "")).strip()
        placement_generation = str(
            value.get("placement_matrix_generation", "")
        ).strip()
        coordinate_frame = str(value.get("coordinate_frame_id", "")).strip()
        placement_digest = str(value.get("placement_matrix_sha256", "")).lower()
        if (
            not assembly_generation
            or not placement_generation
            or not coordinate_frame
            or value.get("centroid_transform_generation") != placement_generation
            or value.get("inertia_transform_generation") != placement_generation
            or value.get("centroid_coordinate_frame_id") != coordinate_frame
            or value.get("inertia_coordinate_frame_id") != coordinate_frame
            or not valid_sha256(placement_digest)
            or value.get("centroid_placement_matrix_sha256") != placement_digest
            or value.get("inertia_placement_matrix_sha256") != placement_digest
        ):
            return None
        return assembly_generation, placement_generation, coordinate_frame, placement_digest

    reference_assembly_coordinates = {
        str(row.get("name", "")): assembly_mass_property_coordinate_identity(row)
        for row in reference
    }
    assembly_coordinate_identity_ok = not assembly_coordinate_evidence_present
    if assembly_coordinate_evidence_present:
        assembly_coordinate_identity_ok = bool(reference_assembly_coordinates) and all(
            value is not None for value in reference_assembly_coordinates.values()
        )
        for _, rows in normalized_sets:
            assembly_coordinate_identity_ok = assembly_coordinate_identity_ok and all(
                assembly_mass_property_coordinate_identity(row)
                == reference_assembly_coordinates.get(str(row.get("name", "")))
                for row in rows
            )

    final_shape_evidence_present = any(
        row.get("boolean_final_shape_identity") is not None for row in identity_rows
    )

    def boolean_final_shape_identity(row):
        value = row.get("boolean_final_shape_identity")
        if not isinstance(value, dict):
            return None
        pre_heal = str(value.get("pre_heal_brep_sha256", "")).lower()
        final = str(value.get("final_brep_sha256", "")).lower()
        final_generation = str(value.get("final_shape_generation", "")).strip()
        if (
            not value.get("boolean_result_generation")
            or not value.get("healing_generation")
            or not final_generation
            or not valid_sha256(pre_heal)
            or not valid_sha256(final)
            or pre_heal == final
            or value.get("mass_property_shape_generation") != final_generation
            or value.get("validity_shape_generation") != final_generation
            or value.get("topology_shape_generation") != final_generation
            or value.get("mass_property_brep_sha256") != final
            or value.get("validity_brep_sha256") != final
            or value.get("topology_brep_sha256") != final
        ):
            return None
        return final_generation, final

    reference_final_shapes = {
        str(row.get("name", "")): boolean_final_shape_identity(row)
        for row in reference
    }
    final_shape_identity_ok = not final_shape_evidence_present
    if final_shape_evidence_present:
        final_shape_identity_ok = bool(reference_final_shapes) and all(
            value is not None for value in reference_final_shapes.values()
        )
        for _, rows in normalized_sets:
            final_shape_identity_ok = final_shape_identity_ok and all(
                boolean_final_shape_identity(row)
                == reference_final_shapes.get(str(row.get("name", "")))
                for row in rows
            )

    tessellation_unit_evidence_present = any(
        row.get("tessellation_tolerance_unit_identity") is not None for row in identity_rows
    )

    def tessellation_tolerance_unit_identity(row):
        value = row.get("tessellation_tolerance_unit_identity")
        if not isinstance(value, dict):
            return None
        units = {"m": 1.0, "mm": 1.0e-3, "um": 1.0e-6}
        unit = str(value.get("linear_deflection_unit", "")).strip()
        evaluation_unit = str(value.get("area_evaluation_deflection_unit", "")).strip()
        try:
            scale = float(value.get("linear_deflection_scale_to_m"))
            evaluation_scale = float(value.get("area_evaluation_deflection_scale_to_m"))
            deflection = float(value.get("linear_deflection_value"))
        except (TypeError, ValueError):
            return None
        generation = str(value.get("tessellation_generation", "")).strip()
        if (
            unit not in units or evaluation_unit != unit or deflection <= 0.0
            or not math.isclose(scale, units[unit], rel_tol=0.0, abs_tol=0.0)
            or not math.isclose(evaluation_scale, scale, rel_tol=0.0, abs_tol=0.0)
            or not generation or value.get("surface_area_generation") != generation
        ):
            return None
        return generation, deflection * scale

    reference_tessellation_units = {str(row.get("name", "")): tessellation_tolerance_unit_identity(row) for row in reference}
    tessellation_unit_identity_ok = not tessellation_unit_evidence_present
    if tessellation_unit_evidence_present:
        tessellation_unit_identity_ok = bool(reference_tessellation_units) and all(value is not None for value in reference_tessellation_units.values())
        for _, rows in normalized_sets:
            tessellation_unit_identity_ok = tessellation_unit_identity_ok and all(
                tessellation_tolerance_unit_identity(row) == reference_tessellation_units.get(str(row.get("name", ""))) for row in rows
            )

    label_topology_evidence_present = any(row.get("compound_label_topology_identity") is not None for row in identity_rows)

    def compound_label_topology_identity(row):
        value = row.get("compound_label_topology_identity")
        if not isinstance(value, dict):
            return None
        generation = str(value.get("boolean_generation", "")).strip()
        digest = str(value.get("final_shape_sha256", "")).lower()
        try:
            topology_index = int(value.get("topology_index"))
            label_index = int(value.get("label_topology_index"))
        except (TypeError, ValueError):
            return None
        if (
            not generation or value.get("label_table_boolean_generation") != generation
            or value.get("selector_boolean_generation") != generation or not value.get("label")
            or topology_index < 0 or label_index != topology_index or not valid_sha256(digest)
            or value.get("selected_subshape_parent_sha256") != digest
        ):
            return None
        return generation, str(value.get("label")), topology_index, digest

    reference_label_topologies = {str(row.get("name", "")): compound_label_topology_identity(row) for row in reference}
    label_topology_identity_ok = not label_topology_evidence_present
    if label_topology_evidence_present:
        label_topology_identity_ok = bool(reference_label_topologies) and all(value is not None for value in reference_label_topologies.values())
        for _, rows in normalized_sets:
            label_topology_identity_ok = label_topology_identity_ok and all(
                compound_label_topology_identity(row) == reference_label_topologies.get(str(row.get("name", ""))) for row in rows
            )

    mass_unit_evidence_present = any(
        row.get("center_of_mass_density_length_unit_identity") is not None
        for row in identity_rows
    )

    def center_of_mass_density_length_unit_identity(row):
        value = row.get("center_of_mass_density_length_unit_identity")
        if not isinstance(value, dict):
            return None
        length_units = {"m": 1.0, "cm": 1.0e-2, "mm": 1.0e-3}
        geometry_unit = str(value.get("geometry_length_unit", "")).strip()
        center_unit = str(value.get("center_of_mass_length_unit", "")).strip()
        volume_unit = str(value.get("volume_length_unit", "")).strip()
        generation = str(value.get("mass_property_generation", "")).strip()
        try:
            geometry_scale = float(value.get("geometry_length_scale_to_m"))
            center_scale = float(value.get("center_of_mass_length_scale_to_m"))
            volume_scale = float(value.get("volume_length_scale_to_m"))
            density_scale = float(value.get("density_scale_to_kg_per_m3"))
        except (TypeError, ValueError):
            return None
        expected_scale = length_units.get(geometry_unit)
        if (
            expected_scale is None
            or center_unit != geometry_unit
            or volume_unit != geometry_unit
            or not math.isclose(geometry_scale, expected_scale, rel_tol=0.0, abs_tol=0.0)
            or not math.isclose(center_scale, geometry_scale, rel_tol=0.0, abs_tol=0.0)
            or not math.isclose(volume_scale, geometry_scale, rel_tol=0.0, abs_tol=0.0)
            or value.get("density_unit") != "kg/m^3"
            or not math.isclose(density_scale, 1.0, rel_tol=0.0, abs_tol=0.0)
            or value.get("reported_mass_unit") != "kg"
            or not generation
            or value.get("geometry_generation") != generation
            or value.get("density_generation") != generation
        ):
            return None
        return generation, geometry_unit, geometry_scale

    reference_mass_units = {
        str(row.get("name", "")): center_of_mass_density_length_unit_identity(row)
        for row in reference
    }
    mass_unit_identity_ok = not mass_unit_evidence_present
    if mass_unit_evidence_present:
        mass_unit_identity_ok = bool(reference_mass_units) and all(
            value is not None for value in reference_mass_units.values()
        )
        for _, rows in normalized_sets:
            mass_unit_identity_ok = mass_unit_identity_ok and all(
                center_of_mass_density_length_unit_identity(row)
                == reference_mass_units.get(str(row.get("name", "")))
                for row in rows
            )

    periodic_selector_evidence_present = any(
        row.get("periodic_face_selector_fillet_topology_identity") is not None
        for row in identity_rows
    )

    def periodic_face_selector_fillet_topology_identity(row):
        value = row.get("periodic_face_selector_fillet_topology_identity")
        if not isinstance(value, dict):
            return None
        fillet_generation = str(value.get("final_fillet_generation", "")).strip()
        topology_generation = str(value.get("final_topology_generation", "")).strip()
        source_ids = list(value.get("source_face_ids") or [])
        selected_source_ids = list(value.get("selected_source_face_ids") or [])
        target_ids = list(value.get("target_face_ids") or [])
        selected_target_ids = list(value.get("selected_target_face_ids") or [])
        digest = str(value.get("final_shape_sha256", "")).lower()
        if (
            not fillet_generation
            or not topology_generation
            or value.get("selector_topology_generation") != topology_generation
            or value.get("periodic_pair_topology_generation") != topology_generation
            or not source_ids
            or not target_ids
            or len(set(source_ids)) != len(source_ids)
            or len(set(target_ids)) != len(target_ids)
            or selected_source_ids != source_ids
            or selected_target_ids != target_ids
            or not valid_sha256(digest)
            or value.get("selector_parent_shape_sha256") != digest
        ):
            return None
        return fillet_generation, topology_generation, tuple(source_ids), tuple(target_ids), digest

    reference_periodic_selectors = {
        str(row.get("name", "")): periodic_face_selector_fillet_topology_identity(row)
        for row in reference
    }
    periodic_selector_identity_ok = not periodic_selector_evidence_present
    if periodic_selector_evidence_present:
        periodic_selector_identity_ok = bool(reference_periodic_selectors) and all(
            value is not None for value in reference_periodic_selectors.values()
        )
        for _, rows in normalized_sets:
            periodic_selector_identity_ok = periodic_selector_identity_ok and all(
                periodic_face_selector_fillet_topology_identity(row)
                == reference_periodic_selectors.get(str(row.get("name", "")))
                for row in rows
            )

    boolean_tolerance_evidence_present = any(
        row.get("boolean_tolerance_length_unit_identity") is not None
        for row in identity_rows
    )

    def boolean_tolerance_length_unit_identity(row):
        value = row.get("boolean_tolerance_length_unit_identity")
        if not isinstance(value, dict):
            return None
        units = {"m": 1.0, "cm": 1.0e-2, "mm": 1.0e-3, "um": 1.0e-6}
        model_unit = str(value.get("model_length_unit", "")).strip()
        tolerance_unit = str(value.get("tolerance_unit", "")).strip()
        generation = str(value.get("boolean_generation", "")).strip()
        input_digest = str(value.get("input_shape_sha256", "")).lower()
        try:
            tolerance_value = float(value.get("tolerance_value"))
            tolerance_scale = float(value.get("tolerance_scale_to_m"))
            kernel_tolerance = float(value.get("kernel_tolerance_m"))
        except (TypeError, ValueError):
            return None
        expected_scale = units.get(tolerance_unit)
        tolerance_m = tolerance_value * tolerance_scale
        if (
            model_unit not in units
            or expected_scale is None
            or not math.isclose(
                tolerance_scale, expected_scale, rel_tol=0.0, abs_tol=0.0
            )
            or not math.isfinite(tolerance_m)
            or tolerance_m <= 0.0
            or not math.isclose(
                kernel_tolerance,
                tolerance_m,
                rel_tol=1.0e-12,
                abs_tol=1.0e-30,
            )
            or not generation
            or value.get("result_geometry_generation") != generation
            or not valid_sha256(input_digest)
            or value.get("boolean_input_shape_sha256") != input_digest
        ):
            return None
        return generation, model_unit, tolerance_m, input_digest

    reference_boolean_tolerances = {
        str(row.get("name", "")): boolean_tolerance_length_unit_identity(row)
        for row in reference
    }
    boolean_tolerance_identity_ok = not boolean_tolerance_evidence_present
    if boolean_tolerance_evidence_present:
        boolean_tolerance_identity_ok = bool(reference_boolean_tolerances) and all(
            value is not None for value in reference_boolean_tolerances.values()
        )
        for _, rows in normalized_sets:
            boolean_tolerance_identity_ok = boolean_tolerance_identity_ok and all(
                boolean_tolerance_length_unit_identity(row)
                == reference_boolean_tolerances.get(str(row.get("name", "")))
                for row in rows
            )

    nested_placement_evidence_present = any(
        row.get("nested_assembly_placement_order_identity") is not None
        for row in identity_rows
    )

    def nested_assembly_placement_order_identity(row):
        value = row.get("nested_assembly_placement_order_identity")
        if not isinstance(value, dict):
            return None
        generation = str(value.get("assembly_generation", "")).strip()
        placement_digest = str(value.get("placement_chain_sha256", "")).lower()
        order = str(value.get("multiplication_order", "")).strip()
        if (
            not generation
            or value.get("parent_placement_generation") != generation
            or value.get("child_placement_generation") != generation
            or value.get("world_placement_generation") != generation
            or order != "parent_then_child"
            or value.get("applied_multiplication_order") != order
            or not valid_sha256(placement_digest)
            or value.get("world_transform_sha256") != placement_digest
        ):
            return None
        return generation, order, placement_digest

    reference_nested_placements = {
        str(row.get("name", "")): nested_assembly_placement_order_identity(row)
        for row in reference
    }
    nested_placement_identity_ok = not nested_placement_evidence_present
    if nested_placement_evidence_present:
        nested_placement_identity_ok = bool(reference_nested_placements) and all(
            value is not None for value in reference_nested_placements.values()
        )
        for _, rows in normalized_sets:
            nested_placement_identity_ok = nested_placement_identity_ok and all(
                nested_assembly_placement_order_identity(row)
                == reference_nested_placements.get(str(row.get("name", "")))
                for row in rows
            )

    mass_inertia_frame_evidence_present = any(
        row.get("mass_inertia_reference_frame_placement_identity") is not None
        for row in identity_rows
    )

    def mass_inertia_reference_frame_placement_identity(row):
        value = row.get("mass_inertia_reference_frame_placement_identity")
        if not isinstance(value, dict):
            return None
        shape_generation = str(value.get("shape_generation", "")).strip()
        placement_generation = str(value.get("placement_generation", "")).strip()
        frame = str(value.get("mass_reference_frame", "")).strip()
        shape_digest = str(value.get("placed_shape_sha256", "")).lower()
        if (
            not shape_generation
            or value.get("mass_property_shape_generation") != shape_generation
            or not placement_generation
            or value.get("mass_property_placement_generation")
            != placement_generation
            or frame != "world"
            or value.get("inertia_reference_frame") != frame
            or value.get("center_of_mass_reference_frame") != frame
            or not valid_sha256(shape_digest)
            or value.get("mass_property_shape_sha256") != shape_digest
        ):
            return None
        return shape_generation, placement_generation, frame, shape_digest

    reference_mass_inertia_frames = {
        str(row.get("name", "")): mass_inertia_reference_frame_placement_identity(row)
        for row in reference
    }
    mass_inertia_frame_identity_ok = not mass_inertia_frame_evidence_present
    if mass_inertia_frame_evidence_present:
        mass_inertia_frame_identity_ok = bool(reference_mass_inertia_frames) and all(
            value is not None for value in reference_mass_inertia_frames.values()
        )
        for _, rows in normalized_sets:
            mass_inertia_frame_identity_ok = mass_inertia_frame_identity_ok and all(
                mass_inertia_reference_frame_placement_identity(row)
                == reference_mass_inertia_frames.get(str(row.get("name", "")))
                for row in rows
            )

    loft_seam_evidence_present = any(
        row.get("loft_wire_correspondence_seam_identity") is not None
        for row in identity_rows
    )

    def loft_wire_correspondence_seam_identity(row):
        value = row.get("loft_wire_correspondence_seam_identity")
        if not isinstance(value, dict):
            return None
        loft_generation = str(value.get("loft_generation", "")).strip()
        seam_generation = str(
            value.get("seam_normalization_generation", "")
        ).strip()
        wire_ids = list(value.get("section_wire_ids") or [])
        digest = str(value.get("wire_correspondence_sha256", "")).lower()
        if (
            not loft_generation
            or value.get("section_wire_loft_generation") != loft_generation
            or not seam_generation
            or value.get("wire_correspondence_seam_generation") != seam_generation
            or not wire_ids
            or len(set(wire_ids)) != len(wire_ids)
            or list(value.get("loft_section_wire_ids") or []) != wire_ids
            or not valid_sha256(digest)
            or value.get("loft_wire_correspondence_sha256") != digest
        ):
            return None
        return loft_generation, seam_generation, tuple(wire_ids), digest

    reference_loft_seams = {
        str(row.get("name", "")): loft_wire_correspondence_seam_identity(row)
        for row in reference
    }
    loft_seam_identity_ok = not loft_seam_evidence_present
    if loft_seam_evidence_present:
        loft_seam_identity_ok = bool(reference_loft_seams) and all(
            value is not None for value in reference_loft_seams.values()
        )
        for _, rows in normalized_sets:
            loft_seam_identity_ok = loft_seam_identity_ok and all(
                loft_wire_correspondence_seam_identity(row)
                == reference_loft_seams.get(str(row.get("name", "")))
                for row in rows
            )

    boolean_unit_generation_evidence_present = any(
        row.get("boolean_tolerance_model_length_unit_generation_identity")
        is not None
        for row in identity_rows
    )

    def boolean_tolerance_model_length_unit_generation_identity(row):
        value = row.get(
            "boolean_tolerance_model_length_unit_generation_identity"
        )
        if not isinstance(value, dict):
            return None
        units = {"m": 1.0, "cm": 1.0e-2, "mm": 1.0e-3, "um": 1.0e-6}
        generation = str(value.get("model_length_unit_generation", "")).strip()
        unit = str(value.get("model_length_unit", "")).strip()
        digest = str(value.get("boolean_tolerance_sha256", "")).lower()
        try:
            tolerance_value = float(value.get("boolean_tolerance_value"))
            tolerance_si = float(value.get("boolean_tolerance_si_m"))
            result_tolerance_si = float(
                value.get("boolean_result_tolerance_si_m")
            )
        except (TypeError, ValueError):
            return None
        expected_scale = units.get(unit)
        if (
            not generation
            or value.get("boolean_tolerance_unit_generation") != generation
            or value.get("boolean_result_unit_generation") != generation
            or expected_scale is None
            or value.get("tolerance_length_unit") != unit
            or value.get("boolean_result_length_unit") != unit
            or not math.isfinite(tolerance_value)
            or tolerance_value <= 0.0
            or not math.isclose(
                tolerance_value * expected_scale,
                tolerance_si,
                rel_tol=1.0e-12,
                abs_tol=1.0e-30,
            )
            or not math.isclose(
                result_tolerance_si,
                tolerance_si,
                rel_tol=1.0e-12,
                abs_tol=1.0e-30,
            )
            or not valid_sha256(digest)
            or value.get("boolean_result_tolerance_sha256") != digest
        ):
            return None
        return generation, unit, tolerance_si, digest

    reference_boolean_unit_generations = {
        str(row.get("name", "")): (
            boolean_tolerance_model_length_unit_generation_identity(row)
        )
        for row in reference
    }
    boolean_unit_generation_identity_ok = not (
        boolean_unit_generation_evidence_present
    )
    if boolean_unit_generation_evidence_present:
        boolean_unit_generation_identity_ok = bool(
            reference_boolean_unit_generations
        ) and all(
            value is not None
            for value in reference_boolean_unit_generations.values()
        )
        for _, rows in normalized_sets:
            boolean_unit_generation_identity_ok = (
                boolean_unit_generation_identity_ok
                and all(
                    boolean_tolerance_model_length_unit_generation_identity(row)
                    == reference_boolean_unit_generations.get(
                        str(row.get("name", ""))
                    )
                    for row in rows
                )
            )

    assembly_density_evidence_present = any(
        row.get("assembly_center_of_mass_part_density_mapping_identity")
        is not None
        for row in identity_rows
    )

    def assembly_center_of_mass_part_density_mapping_identity(row):
        value = row.get("assembly_center_of_mass_part_density_mapping_identity")
        if not isinstance(value, dict):
            return None
        generation = str(
            value.get("assembly_configuration_generation", "")
        ).strip()
        names = [str(name) for name in value.get("part_names", [])]
        density_names = [
            str(name) for name in value.get("density_part_names", [])
        ]
        digest = str(value.get("density_mapping_sha256", "")).lower()
        try:
            densities = [
                float(number) for number in value.get("part_densities_kg_m3", [])
            ]
            used_densities = [
                float(number)
                for number in value.get(
                    "center_of_mass_density_values_kg_m3", []
                )
            ]
        except (TypeError, ValueError):
            return None
        if (
            not generation
            or value.get("part_density_mapping_generation") != generation
            or value.get("center_of_mass_configuration_generation") != generation
            or not names
            or len(names) != len(densities)
            or len(set(names)) != len(names)
            or density_names != names
            or not all(math.isfinite(number) and number > 0.0 for number in densities)
            or used_densities != densities
            or not valid_sha256(digest)
            or value.get("center_of_mass_density_mapping_sha256") != digest
        ):
            return None
        return generation, tuple(names), tuple(densities), digest

    reference_assembly_density_maps = {
        str(row.get("name", "")): (
            assembly_center_of_mass_part_density_mapping_identity(row)
        )
        for row in reference
    }
    assembly_density_identity_ok = not assembly_density_evidence_present
    if assembly_density_evidence_present:
        assembly_density_identity_ok = bool(reference_assembly_density_maps) and all(
            value is not None for value in reference_assembly_density_maps.values()
        )
        for _, rows in normalized_sets:
            assembly_density_identity_ok = assembly_density_identity_ok and all(
                assembly_center_of_mass_part_density_mapping_identity(row)
                == reference_assembly_density_maps.get(str(row.get("name", "")))
                for row in rows
            )

    nested_transform_evidence_present = any(
        row.get("nested_assembly_location_transform_composition_identity") is not None
        for row in identity_rows
    )

    def nested_assembly_location_transform_composition_identity(row):
        value = row.get("nested_assembly_location_transform_composition_identity")
        if not isinstance(value, dict):
            return None
        assembly_generation = str(value.get("assembly_generation", "")).strip()
        transform_order_generation = str(
            value.get("transform_order_generation", "")
        ).strip()
        child_paths = [str(path) for path in value.get("child_paths", [])]
        location_paths = [
            str(path) for path in value.get("location_child_paths", [])
        ]
        local_digests = [
            str(digest).lower() for digest in value.get("local_transform_sha256", [])
        ]
        composed_digests = [
            str(digest).lower()
            for digest in value.get("composed_transform_sha256", [])
        ]
        resolved_digests = [
            str(digest).lower()
            for digest in value.get("resolved_composed_transform_sha256", [])
        ]
        tree_digest = str(value.get("location_tree_sha256", "")).lower()
        order = str(value.get("composition_order", "")).strip()
        if (
            not assembly_generation
            or value.get("location_tree_assembly_generation")
            != assembly_generation
            or value.get("composition_assembly_generation") != assembly_generation
            or not transform_order_generation
            or value.get("location_transform_order_generation")
            != transform_order_generation
            or not child_paths
            or len(set(child_paths)) != len(child_paths)
            or location_paths != child_paths
            or len(local_digests) != len(child_paths)
            or len(composed_digests) != len(child_paths)
            or resolved_digests != composed_digests
            or not all(valid_sha256(digest) for digest in local_digests)
            or not all(valid_sha256(digest) for digest in composed_digests)
            or order != "parent_then_child"
            or value.get("resolved_composition_order") != order
            or not valid_sha256(tree_digest)
            or value.get("resolved_location_tree_sha256") != tree_digest
        ):
            return None
        return (
            assembly_generation,
            transform_order_generation,
            tuple(child_paths),
            tuple(local_digests),
            tuple(composed_digests),
            order,
            tree_digest,
        )

    reference_nested_transform_maps = {
        str(row.get("name", "")): (
            nested_assembly_location_transform_composition_identity(row)
        )
        for row in reference
    }
    nested_transform_identity_ok = not nested_transform_evidence_present
    if nested_transform_evidence_present:
        nested_transform_identity_ok = bool(reference_nested_transform_maps) and all(
            value is not None for value in reference_nested_transform_maps.values()
        )
        for _, rows in normalized_sets:
            nested_transform_identity_ok = nested_transform_identity_ok and all(
                nested_assembly_location_transform_composition_identity(row)
                == reference_nested_transform_maps.get(str(row.get("name", "")))
                for row in rows
            )

    retained_face_history_evidence_present = any(
        row.get("boolean_retained_face_name_history_refine_identity") is not None
        for row in identity_rows
    )

    def boolean_retained_face_name_history_refine_identity(row):
        value = row.get("boolean_retained_face_name_history_refine_identity")
        if not isinstance(value, dict):
            return None
        boolean_generation = str(value.get("boolean_generation", "")).strip()
        refine_generation = str(value.get("refine_generation", "")).strip()
        names = [str(name) for name in value.get("retained_face_names", [])]
        resolved_names = [str(name) for name in value.get("resolved_face_names", [])]
        try:
            face_ids = [int(face_id) for face_id in value.get("retained_face_ids", [])]
            resolved_face_ids = [
                int(face_id) for face_id in value.get("resolved_face_ids", [])
            ]
        except (TypeError, ValueError):
            return None
        history_digest = str(value.get("topology_history_sha256", "")).lower()
        if (
            not boolean_generation
            or value.get("retained_name_boolean_generation") != boolean_generation
            or not refine_generation
            or value.get("retained_name_refine_generation") != refine_generation
            or value.get("topology_history_refine_generation") != refine_generation
            or not names
            or len(set(names)) != len(names)
            or any(not name.strip() for name in names)
            or resolved_names != names
            or len(face_ids) != len(names)
            or len(set(face_ids)) != len(face_ids)
            or any(face_id <= 0 for face_id in face_ids)
            or resolved_face_ids != face_ids
            or not valid_sha256(history_digest)
            or value.get("resolved_topology_history_sha256") != history_digest
        ):
            return None
        return (
            boolean_generation,
            refine_generation,
            tuple(names),
            tuple(face_ids),
            history_digest,
        )

    reference_retained_face_histories = {
        str(row.get("name", "")): (
            boolean_retained_face_name_history_refine_identity(row)
        )
        for row in reference
    }
    retained_face_history_identity_ok = not retained_face_history_evidence_present
    if retained_face_history_evidence_present:
        retained_face_history_identity_ok = bool(
            reference_retained_face_histories
        ) and all(
            value is not None for value in reference_retained_face_histories.values()
        )
        for _, rows in normalized_sets:
            retained_face_history_identity_ok = (
                retained_face_history_identity_ok
                and all(
                    boolean_retained_face_name_history_refine_identity(row)
                    == reference_retained_face_histories.get(str(row.get("name", "")))
                    for row in rows
                )
            )

    boolean_subshape_evidence_present = any(
        row.get("boolean_history_subshape_label_fillet_order_identity") is not None
        for row in identity_rows
    )

    def boolean_history_subshape_label_fillet_order_identity(row):
        value = row.get("boolean_history_subshape_label_fillet_order_identity")
        if not isinstance(value, dict):
            return None
        boolean_generation = str(value.get("boolean_generation", "")).strip()
        fillet_generation = str(value.get("fillet_generation", "")).strip()
        labels = [str(label) for label in value.get("subshape_labels", [])]
        resolved_labels = [
            str(label) for label in value.get("resolved_subshape_labels", [])
        ]
        try:
            face_ids = [int(item) for item in value.get("history_face_ids", [])]
            resolved_face_ids = [
                int(item) for item in value.get("resolved_face_ids", [])
            ]
            edge_ids = [int(item) for item in value.get("fillet_edge_ids", [])]
            resolved_edge_ids = [
                int(item) for item in value.get("resolved_fillet_edge_ids", [])
            ]
        except (TypeError, ValueError):
            return None
        digest = str(value.get("subshape_history_sha256", "")).lower()
        if (
            not boolean_generation
            or value.get("history_boolean_generation") != boolean_generation
            or not fillet_generation
            or value.get("label_fillet_generation") != fillet_generation
            or value.get("edge_order_fillet_generation") != fillet_generation
            or not labels
            or len(set(labels)) != len(labels)
            or any(not label.strip() for label in labels)
            or resolved_labels != labels
            or len(face_ids) != len(labels)
            or len(set(face_ids)) != len(face_ids)
            or resolved_face_ids != face_ids
            or not edge_ids
            or len(set(edge_ids)) != len(edge_ids)
            or resolved_edge_ids != edge_ids
            or not valid_sha256(digest)
            or value.get("resolved_subshape_history_sha256") != digest
        ):
            return None
        return (
            boolean_generation,
            fillet_generation,
            tuple(labels),
            tuple(face_ids),
            tuple(edge_ids),
            digest,
        )

    reference_boolean_subshape_histories = {
        str(row.get("name", "")): boolean_history_subshape_label_fillet_order_identity(
            row
        )
        for row in reference
    }
    boolean_subshape_identity_ok = not boolean_subshape_evidence_present
    if boolean_subshape_evidence_present:
        boolean_subshape_identity_ok = bool(
            reference_boolean_subshape_histories
        ) and all(
            value is not None
            for value in reference_boolean_subshape_histories.values()
        )
        for _, rows in normalized_sets:
            boolean_subshape_identity_ok = boolean_subshape_identity_ok and all(
                boolean_history_subshape_label_fillet_order_identity(row)
                == reference_boolean_subshape_histories.get(str(row.get("name", "")))
                for row in rows
            )

    assembly_mate_evidence_present = any(
        row.get("assembly_mate_frame_unit_location_generation_identity") is not None
        for row in identity_rows
    )

    def assembly_mate_frame_unit_location_generation_identity(row):
        value = row.get("assembly_mate_frame_unit_location_generation_identity")
        if not isinstance(value, dict):
            return None
        generation = str(value.get("assembly_generation", "")).strip()
        length_unit = str(value.get("length_unit", "")).strip()
        names = [str(name) for name in value.get("mate_names", [])]
        resolved_names = [str(name) for name in value.get("resolved_mate_names", [])]
        local_digests = [
            str(item).lower() for item in value.get("local_frame_sha256", [])
        ]
        resolved_local_digests = [
            str(item).lower()
            for item in value.get("resolved_local_frame_sha256", [])
        ]
        parent_digests = [
            str(item).lower() for item in value.get("parent_location_sha256", [])
        ]
        resolved_parent_digests = [
            str(item).lower()
            for item in value.get("resolved_parent_location_sha256", [])
        ]
        digest = str(value.get("mate_resolution_sha256", "")).lower()
        if (
            not generation
            or value.get("mate_frame_assembly_generation") != generation
            or value.get("parent_location_assembly_generation") != generation
            or value.get("unit_assembly_generation") != generation
            or length_unit not in {"m", "cm", "mm"}
            or value.get("mate_frame_length_unit") != length_unit
            or value.get("parent_location_length_unit") != length_unit
            or not names
            or len(set(names)) != len(names)
            or any(not name.strip() for name in names)
            or resolved_names != names
            or len(local_digests) != len(names)
            or resolved_local_digests != local_digests
            or len(parent_digests) != len(names)
            or resolved_parent_digests != parent_digests
            or not all(valid_sha256(item) for item in local_digests)
            or not all(valid_sha256(item) for item in parent_digests)
            or not valid_sha256(digest)
            or value.get("resolved_mate_resolution_sha256") != digest
        ):
            return None
        return (
            generation,
            length_unit,
            tuple(names),
            tuple(local_digests),
            tuple(parent_digests),
            digest,
        )

    reference_assembly_mate_maps = {
        str(row.get("name", "")): (
            assembly_mate_frame_unit_location_generation_identity(row)
        )
        for row in reference
    }
    assembly_mate_identity_ok = not assembly_mate_evidence_present
    if assembly_mate_evidence_present:
        assembly_mate_identity_ok = bool(reference_assembly_mate_maps) and all(
            value is not None for value in reference_assembly_mate_maps.values()
        )
        for _, rows in normalized_sets:
            assembly_mate_identity_ok = assembly_mate_identity_ok and all(
                assembly_mate_frame_unit_location_generation_identity(row)
                == reference_assembly_mate_maps.get(str(row.get("name", "")))
                for row in rows
            )

    mass_density_location_evidence_present = any(
        row.get("mass_properties_density_unit_location_generation_identity")
        is not None
        for row in identity_rows
    )

    def mass_properties_density_unit_location_generation_identity(row):
        value = row.get(
            "mass_properties_density_unit_location_generation_identity"
        )
        if not isinstance(value, dict):
            return None
        assembly_generation = str(
            value.get("assembly_generation", "")
        ).strip()
        density_generation = str(
            value.get("density_mapping_generation", "")
        ).strip()
        location_generation = str(
            value.get("part_location_generation", "")
        ).strip()
        density_unit = str(value.get("density_unit", "")).strip()
        part_names = [str(name) for name in value.get("part_names", [])]
        resolved_names = [
            str(name) for name in value.get("resolved_part_names", [])
        ]
        location_digests = [
            str(item).lower()
            for item in value.get("part_location_sha256", [])
        ]
        resolved_location_digests = [
            str(item).lower()
            for item in value.get("resolved_part_location_sha256", [])
        ]
        table_digest = str(
            value.get("mass_property_table_sha256", "")
        ).lower()
        if (
            not assembly_generation
            or any(
                value.get(key) != assembly_generation
                for key in (
                    "mass_assembly_generation",
                    "center_of_mass_assembly_generation",
                    "inertia_assembly_generation",
                )
            )
            or not density_generation
            or any(
                value.get(key) != density_generation
                for key in (
                    "mass_density_mapping_generation",
                    "center_of_mass_density_mapping_generation",
                    "inertia_density_mapping_generation",
                )
            )
            or not location_generation
            or any(
                value.get(key) != location_generation
                for key in (
                    "mass_part_location_generation",
                    "center_of_mass_part_location_generation",
                    "inertia_part_location_generation",
                )
            )
            or density_unit not in {"kg/m^3", "g/cm^3"}
            or any(
                value.get(key) != density_unit
                for key in (
                    "mass_density_unit",
                    "center_of_mass_density_unit",
                    "inertia_density_unit",
                )
            )
            or not part_names
            or len(set(part_names)) != len(part_names)
            or any(not name.strip() for name in part_names)
            or resolved_names != part_names
            or len(location_digests) != len(part_names)
            or not all(valid_sha256(item) for item in location_digests)
            or resolved_location_digests != location_digests
            or not valid_sha256(table_digest)
            or value.get("resolved_mass_property_table_sha256") != table_digest
        ):
            return None
        return (
            assembly_generation,
            density_generation,
            location_generation,
            density_unit,
            tuple(part_names),
            tuple(location_digests),
            table_digest,
        )

    reference_mass_density_location_maps = {
        str(row.get("name", "")): (
            mass_properties_density_unit_location_generation_identity(row)
        )
        for row in reference
    }
    mass_density_location_identity_ok = not mass_density_location_evidence_present
    if mass_density_location_evidence_present:
        mass_density_location_identity_ok = bool(
            reference_mass_density_location_maps
        ) and all(
            value is not None
            for value in reference_mass_density_location_maps.values()
        )
        for _, rows in normalized_sets:
            mass_density_location_identity_ok = (
                mass_density_location_identity_ok
                and all(
                    mass_properties_density_unit_location_generation_identity(row)
                    == reference_mass_density_location_maps.get(
                        str(row.get("name", ""))
                    )
                    for row in rows
                )
            )

    sweep_frame_evidence_present = any(
        row.get(
            "sweep_path_frame_twist_profile_orientation_generation_identity"
        )
        is not None
        for row in identity_rows
    )

    def sweep_path_frame_twist_profile_orientation_generation_identity(row):
        value = row.get(
            "sweep_path_frame_twist_profile_orientation_generation_identity"
        )
        if not isinstance(value, dict):
            return None
        sweep_generation = str(value.get("sweep_generation", "")).strip()
        path_generation = str(value.get("path_generation", "")).strip()
        profile_generation = str(value.get("profile_generation", "")).strip()
        try:
            parameters = [float(item) for item in value.get("path_parameters", [])]
            frame_parameters = [
                float(item) for item in value.get("frame_path_parameters", [])
            ]
            twists = [float(item) for item in value.get("twist_degrees", [])]
            solid_twists = [
                float(item) for item in value.get("solid_twist_degrees", [])
            ]
        except (TypeError, ValueError):
            return None
        frame_digests = [
            str(item).lower() for item in value.get("path_frame_sha256", [])
        ]
        solid_frame_digests = [
            str(item).lower()
            for item in value.get("solid_path_frame_sha256", [])
        ]
        orientation_digest = str(
            value.get("profile_orientation_sha256", "")
        ).lower()
        solid_digest = str(value.get("swept_solid_sha256", "")).lower()
        if (
            not sweep_generation
            or value.get("solid_sweep_generation") != sweep_generation
            or not path_generation
            or value.get("frame_path_generation") != path_generation
            or value.get("twist_path_generation") != path_generation
            or not profile_generation
            or value.get("orientation_profile_generation") != profile_generation
            or value.get("solid_profile_generation") != profile_generation
            or len(parameters) < 2
            or not all(math.isfinite(item) for item in parameters)
            or not all(left < right for left, right in zip(parameters, parameters[1:]))
            or frame_parameters != parameters
            or len(twists) != len(parameters)
            or not all(math.isfinite(item) for item in twists)
            or solid_twists != twists
            or len(frame_digests) != len(parameters)
            or not all(valid_sha256(item) for item in frame_digests)
            or solid_frame_digests != frame_digests
            or not valid_sha256(orientation_digest)
            or value.get("solid_profile_orientation_sha256")
            != orientation_digest
            or not valid_sha256(solid_digest)
            or value.get("resolved_swept_solid_sha256") != solid_digest
        ):
            return None
        return (
            sweep_generation,
            path_generation,
            profile_generation,
            tuple(parameters),
            tuple(twists),
            tuple(frame_digests),
            orientation_digest,
            solid_digest,
        )

    reference_sweep_frame_maps = {
        str(row.get("name", "")): (
            sweep_path_frame_twist_profile_orientation_generation_identity(row)
        )
        for row in reference
    }
    sweep_frame_identity_ok = not sweep_frame_evidence_present
    if sweep_frame_evidence_present:
        sweep_frame_identity_ok = bool(reference_sweep_frame_maps) and all(
            value is not None for value in reference_sweep_frame_maps.values()
        )
        for _, rows in normalized_sets:
            sweep_frame_identity_ok = sweep_frame_identity_ok and all(
                sweep_path_frame_twist_profile_orientation_generation_identity(row)
                == reference_sweep_frame_maps.get(str(row.get("name", "")))
                for row in rows
            )

    boolean_result_evidence_present = any(
        row.get(
            "boolean_result_solid_orientation_location_label_generation_identity"
        )
        is not None
        for row in identity_rows
    )

    def boolean_result_orientation_location_label_identity(row):
        value = row.get(
            "boolean_result_solid_orientation_location_label_generation_identity"
        )
        if not isinstance(value, dict):
            return None
        generation = str(value.get("boolean_generation", "")).strip()
        operand_generations = [
            str(item).strip() for item in value.get("operand_generations", [])
        ]
        result_operands = [
            str(item).strip()
            for item in value.get("result_operand_generations", [])
        ]
        orientation = str(value.get("solid_orientation", "")).strip()
        labels = [str(item).strip() for item in value.get("semantic_labels", [])]
        resolved_labels = [
            str(item).strip()
            for item in value.get("resolved_semantic_labels", [])
        ]
        location_digest = str(
            value.get("result_location_sha256", "")
        ).lower()
        result_digest = str(value.get("boolean_result_sha256", "")).lower()
        if (
            not generation
            or any(
                value.get(key) != generation
                for key in (
                    "result_boolean_generation",
                    "orientation_boolean_generation",
                    "location_boolean_generation",
                    "label_boolean_generation",
                )
            )
            or len(operand_generations) < 2
            or any(not item for item in operand_generations)
            or result_operands != operand_generations
            or orientation not in {"forward", "reversed"}
            or value.get("resolved_solid_orientation") != orientation
            or not labels
            or len(set(labels)) != len(labels)
            or any(not item for item in labels)
            or resolved_labels != labels
            or not valid_sha256(location_digest)
            or value.get("resolved_result_location_sha256") != location_digest
            or not valid_sha256(result_digest)
            or value.get("resolved_boolean_result_sha256") != result_digest
        ):
            return None
        return (
            generation,
            tuple(operand_generations),
            orientation,
            tuple(labels),
            location_digest,
            result_digest,
        )

    reference_boolean_result_maps = {
        str(row.get("name", "")): boolean_result_orientation_location_label_identity(
            row
        )
        for row in reference
    }
    boolean_result_identity_ok = not boolean_result_evidence_present
    if boolean_result_evidence_present:
        boolean_result_identity_ok = bool(reference_boolean_result_maps) and all(
            value is not None for value in reference_boolean_result_maps.values()
        )
        for _, rows in normalized_sets:
            boolean_result_identity_ok = boolean_result_identity_ok and all(
                boolean_result_orientation_location_label_identity(row)
                == reference_boolean_result_maps.get(str(row.get("name", "")))
                for row in rows
            )

    tessellation_generation_evidence_present = any(
        row.get("tessellation_chord_angle_unit_location_generation_identity")
        is not None
        for row in identity_rows
    )

    def tessellation_chord_angle_unit_location_identity(row):
        value = row.get(
            "tessellation_chord_angle_unit_location_generation_identity"
        )
        if not isinstance(value, dict):
            return None
        shape_generation = str(value.get("shape_generation", "")).strip()
        tessellation_generation = str(
            value.get("tessellation_generation", "")
        ).strip()
        try:
            chord = float(value.get("chord_tolerance"))
            evaluated_chord = float(value.get("evaluated_chord_tolerance"))
            angle = float(value.get("angular_tolerance_deg"))
            evaluated_angle = float(value.get("evaluated_angular_tolerance_deg"))
        except (TypeError, ValueError):
            return None
        length_unit = str(value.get("length_unit", "")).strip()
        location_digest = str(
            value.get("object_location_sha256", "")
        ).lower()
        tessellation_digest = str(
            value.get("tessellation_sha256", "")
        ).lower()
        if (
            not shape_generation
            or value.get("tessellation_shape_generation") != shape_generation
            or not tessellation_generation
            or value.get("metric_tessellation_generation")
            != tessellation_generation
            or value.get("location_tessellation_generation")
            != tessellation_generation
            or not math.isfinite(chord)
            or chord <= 0.0
            or not math.isclose(
                evaluated_chord, chord, rel_tol=0.0, abs_tol=1.0e-18
            )
            or not math.isfinite(angle)
            or not 0.0 < angle <= 180.0
            or not math.isclose(
                evaluated_angle, angle, rel_tol=0.0, abs_tol=1.0e-12
            )
            or length_unit not in {"m", "cm", "mm"}
            or value.get("evaluated_length_unit") != length_unit
            or not valid_sha256(location_digest)
            or value.get("evaluated_object_location_sha256") != location_digest
            or not valid_sha256(tessellation_digest)
            or value.get("evaluated_tessellation_sha256")
            != tessellation_digest
        ):
            return None
        return (
            shape_generation,
            tessellation_generation,
            chord,
            angle,
            length_unit,
            location_digest,
            tessellation_digest,
        )

    reference_tessellation_generation_maps = {
        str(row.get("name", "")): tessellation_chord_angle_unit_location_identity(
            row
        )
        for row in reference
    }
    tessellation_generation_identity_ok = not tessellation_generation_evidence_present
    if tessellation_generation_evidence_present:
        tessellation_generation_identity_ok = bool(
            reference_tessellation_generation_maps
        ) and all(
            value is not None
            for value in reference_tessellation_generation_maps.values()
        )
        for _, rows in normalized_sets:
            tessellation_generation_identity_ok = (
                tessellation_generation_identity_ok
                and all(
                    tessellation_chord_angle_unit_location_identity(row)
                    == reference_tessellation_generation_maps.get(
                        str(row.get("name", ""))
                    )
                    for row in rows
                )
            )

    connector_evidence_present = any(
        row.get("joint_connector_frame_labeled_face_subshape_generation_identity")
        is not None
        for row in identity_rows
    )

    def connector_frame_subshape_identity(row):
        value = row.get(
            "joint_connector_frame_labeled_face_subshape_generation_identity"
        )
        if not isinstance(value, dict):
            return None
        generation = str(value.get("shape_generation", "")).strip()
        face_id = str(value.get("labeled_face_subshape_id", "")).strip()
        face_digest = str(value.get("labeled_face_geometry_sha256", "")).lower()
        location_digest = str(value.get("parent_location_sha256", "")).lower()
        frame_digest = str(value.get("connector_frame_sha256", "")).lower()
        try:
            origin = tuple(float(item) for item in value.get("connector_origin", []))
            evaluated_origin = tuple(
                float(item) for item in value.get("evaluated_connector_origin", [])
            )
            axis = tuple(float(item) for item in value.get("connector_axis", []))
            evaluated_axis = tuple(
                float(item) for item in value.get("evaluated_connector_axis", [])
            )
        except (TypeError, ValueError):
            return None
        if (
            not generation
            or any(
                value.get(key) != generation
                for key in (
                    "label_table_shape_generation",
                    "connector_shape_generation",
                    "location_shape_generation",
                )
            )
            or not face_id
            or value.get("resolved_labeled_face_subshape_id") != face_id
            or not valid_sha256(face_digest)
            or value.get("resolved_labeled_face_geometry_sha256") != face_digest
            or len(origin) != 3
            or not all(math.isfinite(item) for item in origin)
            or evaluated_origin != origin
            or len(axis) != 3
            or not all(math.isfinite(item) for item in axis)
            or math.isclose(sum(item * item for item in axis), 0.0)
            or evaluated_axis != axis
            or not valid_sha256(location_digest)
            or value.get("evaluated_parent_location_sha256") != location_digest
            or not valid_sha256(frame_digest)
            or value.get("evaluated_connector_frame_sha256") != frame_digest
        ):
            return None
        return generation, face_id, face_digest, origin, axis, location_digest, frame_digest

    reference_connector_maps = {
        str(row.get("name", "")): connector_frame_subshape_identity(row)
        for row in reference
    }
    connector_identity_ok = not connector_evidence_present
    if connector_evidence_present:
        connector_identity_ok = bool(reference_connector_maps) and all(
            value is not None for value in reference_connector_maps.values()
        )
        for _, rows in normalized_sets:
            connector_identity_ok = connector_identity_ok and all(
                connector_frame_subshape_identity(row)
                == reference_connector_maps.get(str(row.get("name", "")))
                for row in rows
            )

    inertia_evidence_present = any(
        row.get(
            "inertia_tensor_principal_axes_density_unit_location_generation_identity"
        )
        is not None
        for row in identity_rows
    )

    def inertia_density_location_identity(row):
        value = row.get(
            "inertia_tensor_principal_axes_density_unit_location_generation_identity"
        )
        if not isinstance(value, dict):
            return None
        generation = str(value.get("shape_generation", "")).strip()
        unit = str(value.get("density_unit", "")).strip()
        location_digest = str(value.get("shape_location_sha256", "")).lower()
        tensor_digest = str(value.get("inertia_tensor_sha256", "")).lower()
        axes_digest = str(value.get("principal_axes_sha256", "")).lower()
        mass_digest = str(value.get("mass_property_sha256", "")).lower()
        try:
            density = float(value.get("density_value"))
            evaluated_density = float(value.get("evaluated_density_value"))
            center = tuple(float(item) for item in value.get("center_of_mass", []))
            evaluated_center = tuple(
                float(item) for item in value.get("evaluated_center_of_mass", [])
            )
        except (TypeError, ValueError):
            return None
        if (
            not generation
            or any(
                value.get(key) != generation
                for key in (
                    "density_shape_generation",
                    "mass_property_shape_generation",
                    "location_shape_generation",
                    "principal_axis_shape_generation",
                )
            )
            or not math.isfinite(density)
            or density <= 0.0
            or not math.isclose(evaluated_density, density, rel_tol=0.0, abs_tol=1.0e-12)
            or unit not in {"kg/m^3", "g/cm^3"}
            or value.get("evaluated_density_unit") != unit
            or not valid_sha256(location_digest)
            or value.get("evaluated_shape_location_sha256") != location_digest
            or len(center) != 3
            or not all(math.isfinite(item) for item in center)
            or evaluated_center != center
            or not valid_sha256(tensor_digest)
            or value.get("evaluated_inertia_tensor_sha256") != tensor_digest
            or not valid_sha256(axes_digest)
            or value.get("evaluated_principal_axes_sha256") != axes_digest
            or not valid_sha256(mass_digest)
            or value.get("evaluated_mass_property_sha256") != mass_digest
        ):
            return None
        return (
            generation,
            density,
            unit,
            location_digest,
            center,
            tensor_digest,
            axes_digest,
            mass_digest,
        )

    reference_inertia_maps = {
        str(row.get("name", "")): inertia_density_location_identity(row)
        for row in reference
    }
    inertia_density_location_identity_ok = not inertia_evidence_present
    if inertia_evidence_present:
        inertia_density_location_identity_ok = bool(reference_inertia_maps) and all(
            value is not None for value in reference_inertia_maps.values()
        )
        for _, rows in normalized_sets:
            inertia_density_location_identity_ok = (
                inertia_density_location_identity_ok
                and all(
                    inertia_density_location_identity(row)
                    == reference_inertia_maps.get(str(row.get("name", "")))
                    for row in rows
                )
            )

    inventory = shape_measurement_inventory_summary(reference)
    sets = []
    all_rows = []
    identity_gates = []
    for label, rows_list in normalized_sets:
        identity_gate = shape_name_identity_gate(
            reference,
            rows_list,
            measured_label=label,
        )
        identity_gates.append(identity_gate)
        compared = shape_measurement_comparison_summary(
            reference,
            rows_list,
            rtol=rtol,
            measured_label=label,
            bbox_atol=bbox_atol,
        )
        comparison_rows = compared["rows"]
        failed = [row for row in comparison_rows if not row["passed"]]
        source_summary = {
            "source": label,
            "status": "ok"
            if not failed and identity_gate["status"] == "ok"
            else "needs_attention",
            "name_identity_gate": identity_gate,
            "n_cases": compared["n_cases"],
            "n_passed": compared["n_passed"],
            "n_bbox_compared": compared["n_bbox_compared"],
            "max_volume_rel_error": compared["max_volume_rel_error"],
            "max_area_rel_error": compared["max_area_rel_error"],
            "max_bbox_abs_error": compared["max_bbox_abs_error"],
            "worst_comparisons": worst_shape_measurement_comparison_rows(
                comparison_rows,
                limit=worst_limit,
            ),
            "rows": comparison_rows,
        }
        sets.append(source_summary)
        all_rows.extend(comparison_rows)

    failed_rows = [row for row in all_rows if not row["passed"]]
    failed_identity_gates = [gate for gate in identity_gates if gate["status"] != "ok"]
    checks = {
        "all_reference_shapes_valid": inventory["n_valid"] == inventory["n_shapes"],
        "all_sources_present_and_within_tolerance": not failed_rows,
        "all_sources_preserve_named_shape_identity": not failed_identity_gates,
        "mass_properties_bind_current_brep_revision": revision_identity_ok,
        "assembly_children_match_reference_revision_map": assembly_identity_ok,
        "mass_property_centers_share_reference_frames": frame_identity_ok,
        "face_adjacency_matches_current_brep_revision": topology_identity_ok,
        "compound_volume_uses_physical_union_not_child_sum": (
            compound_volume_identity_ok
        ),
        "center_of_mass_uses_final_placement_transform": (
            placement_transform_identity_ok
        ),
        "mass_properties_follow_final_healed_brep": shape_healing_identity_ok,
        "mirrored_inertia_tensor_uses_final_global_frame": (
            inertia_tensor_identity_ok
        ),
        "assembly_mass_properties_use_final_coordinate_frame": (
            assembly_coordinate_identity_ok
        ),
        "boolean_validity_topology_and_mass_share_final_healed_shape": (
            final_shape_identity_ok
        ),
        "tessellated_area_uses_one_length_unit_tolerance": tessellation_unit_identity_ok,
        "compound_labels_resolve_on_final_boolean_topology": label_topology_identity_ok,
        "center_of_mass_density_and_volume_share_length_unit_covariance": mass_unit_identity_ok,
        "periodic_face_selectors_follow_final_fillet_topology": periodic_selector_identity_ok,
        "boolean_tolerance_uses_one_physical_model_length_basis": (
            boolean_tolerance_identity_ok
        ),
        "nested_assembly_placements_use_parent_then_child_order": (
            nested_placement_identity_ok
        ),
        "mass_inertia_uses_final_world_placement_frame": (
            mass_inertia_frame_identity_ok
        ),
        "loft_sections_use_current_seam_normalized_correspondence": (
            loft_seam_identity_ok
        ),
        "boolean_tolerance_uses_current_model_length_unit_generation": (
            boolean_unit_generation_identity_ok
        ),
        "assembly_center_of_mass_uses_current_part_density_mapping": (
            assembly_density_identity_ok
        ),
        "nested_assembly_locations_use_current_transform_composition_order": (
            nested_transform_identity_ok
        ),
        "boolean_retained_face_names_follow_post_refine_topology_history": (
            retained_face_history_identity_ok
        ),
        "boolean_subshape_labels_follow_current_fillet_edge_order": (
            boolean_subshape_identity_ok
        ),
        "assembly_mates_share_current_frame_unit_and_parent_location": (
            assembly_mate_identity_ok
        ),
        "assembly_mass_properties_share_density_unit_and_part_locations": (
            mass_density_location_identity_ok
        ),
        "swept_solid_uses_current_path_frames_twist_and_profile_orientation": (
            sweep_frame_identity_ok
        ),
        "boolean_result_uses_current_orientation_location_and_labels": (
            boolean_result_identity_ok
        ),
        "tessellation_uses_current_tolerances_units_and_object_location": (
            tessellation_generation_identity_ok
        ),
        "joint_connectors_use_current_labeled_face_frame_and_parent_location": (
            connector_identity_ok
        ),
        "inertia_uses_current_density_unit_location_and_principal_axes": (
            inertia_density_location_identity_ok
        ),
        "assembly_mass_properties_use_current_locations_operands_members_and_density": (
            assembly_location_identity_ok
        ),
        "loft_volume_uses_current_spline_tolerances_and_watertight_shell": (
            loft_volume_identity_ok
        ),
        "transformed_assembly_mass_properties_use_current_transforms_density_units_and_axes": (
            transformed_assembly_identity_ok
        ),
        "fillet_chamfer_topology_uses_current_selection_names_order_and_fingerprint": (
            fillet_chamfer_identity_ok
        ),
        "boolean_result_uses_current_operands_tolerance_healing_topology_and_volume": (
            boolean_operation_identity_ok
        ),
        "assembly_mates_use_current_transforms_dof_solver_and_loop_closure": (
            assembly_kinematics_identity_ok
        ),
        "fillet_chamfer_features_use_current_selectors_topology_names_tolerance_and_shape": (
            feature_selection_identity_ok
        ),
        "mass_properties_use_current_density_center_inertia_frame_and_assembly": (
            mass_frame_identity_ok
        ),
        "boolean_imprints_use_current_interfaces_owners_topology_names_tolerance_and_mass": (
            boolean_imprint_identity_ok
        ),
        "lofts_use_current_sections_wire_orientation_seams_continuity_solid_and_volume": (
            loft_section_identity_ok
        ),
        "shell_offsets_use_current_faces_normals_thickness_join_intersection_shape_and_mass": (
            shell_offset_identity_ok
        ),
        "path_sweeps_use_current_path_frame_transition_profile_orientation_solid_and_volume": (
            path_sweep_identity_ok
        ),
        "sheet_metal_flat_patterns_use_current_bends_neutral_axis_relief_thickness_and_area": (
            sheet_metal_identity_ok
        ),
        "joint_loops_use_current_graph_dofs_limits_frames_closure_configuration_and_swept_volume": (
            joint_loop_identity_ok
        ),
        "helical_sweeps_use_current_pitch_handedness_profile_frame_turns_intersection_mass_and_shape": (
            helical_sweep_identity_ok
        ),
        "boolean_results_use_current_operation_tolerance_operands_history_mass_inertia_and_shape": (
            boolean_history_identity_ok
        ),
        "assembly_occurrences_use_current_locations_densities_units_suppression_mass_center_parallel_axis_inertia_and_shape": (
            assembly_occurrence_identity_ok
        ),
        "lofts_use_current_profile_order_seams_guide_face_lineage_shell_volume_and_shape": (
            loft_lineage_identity_ok
        ),
        "boolean_results_use_current_fuzzy_tolerance_topology_names_face_ancestry_count_volume_centroid_and_shape": (
            boolean_topology_identity_ok
        ),
        "swept_solids_use_current_frenet_frame_twist_transition_orientation_intersection_volume_owner_and_shape": (
            sweep_frame_transition_identity_ok
        ),
        "guided_lofts_use_current_sections_parameterization_guides_seams_mode_intersection_volume_and_shape": (
            guided_loft_identity_ok
        ),
        "mass_properties_use_current_density_units_origin_center_principal_axes_parallel_axis_owner_and_shape": (
            mass_inertia_identity_ok
        ),
        "assembly_mates_close_current_cycles_frames_mass_centers_rotated_inertia_owner_and_shape": (
            assembly_mate_mass_inertia_identity_ok
        ),
        "shell_fillets_use_current_euler_manifold_thickness_volume_area_inertia_convergence_and_brep": (
            shell_fillet_identity_ok
        ),
    }
    issues = []
    if not checks["all_reference_shapes_valid"]:
        issues.append("at least one reference build123d shape is invalid")
    if not checks["all_sources_present_and_within_tolerance"]:
        issues.append("at least one CAD source row is missing or outside tolerance")
    if not checks["all_sources_preserve_named_shape_identity"]:
        issues.append("at least one CAD source has missing, extra, duplicate, or unnamed shapes")
    if not checks["mass_properties_bind_current_brep_revision"]:
        issues.append("at least one mass-property row belongs to a different BREP revision")
    if not checks["assembly_children_match_reference_revision_map"]:
        issues.append("assembly child revisions differ across rows or from the reference")
    if not checks["mass_property_centers_share_reference_frames"]:
        issues.append("mass-property centers use a different reference frame")
    if not checks["face_adjacency_matches_current_brep_revision"]:
        issues.append("face adjacency belongs to another BREP revision")
    if not checks["compound_volume_uses_physical_union_not_child_sum"]:
        issues.append("compound volume is not bound to the physical union topology")
    if not checks["center_of_mass_uses_final_placement_transform"]:
        issues.append("center of mass belongs to a previous placement transform")
    if not checks["mass_properties_follow_final_healed_brep"]:
        issues.append("mass properties belong to the pre-heal or another BREP generation")
    if not checks["mirrored_inertia_tensor_uses_final_global_frame"]:
        issues.append("mirrored inertia tensor is not expressed in the final global frame")
    if not checks["assembly_mass_properties_use_final_coordinate_frame"]:
        issues.append("assembly mass properties are not expressed in the final placement frame")
    if not checks["boolean_validity_topology_and_mass_share_final_healed_shape"]:
        issues.append("boolean validity, topology, and mass belong to different heal generations")
    if not checks["center_of_mass_density_and_volume_share_length_unit_covariance"]:
        issues.append("center of mass, density, volume, and mass use inconsistent unit covariance")
    if not checks["periodic_face_selectors_follow_final_fillet_topology"]:
        issues.append("periodic face selectors belong to a pre-fillet topology")
    if not checks["boolean_tolerance_uses_one_physical_model_length_basis"]:
        issues.append("boolean tolerance is not bound to the model length basis")
    if not checks["nested_assembly_placements_use_parent_then_child_order"]:
        issues.append("nested assembly placements use a stale or reversed transform order")
    if not checks["mass_inertia_uses_final_world_placement_frame"]:
        issues.append("mass and inertia properties do not use the final world placement frame")
    if not checks["loft_sections_use_current_seam_normalized_correspondence"]:
        issues.append("loft wire correspondence predates seam normalization")
    if not checks["boolean_tolerance_uses_current_model_length_unit_generation"]:
        issues.append("boolean tolerance predates the active model length-unit generation")
    if not checks["assembly_center_of_mass_uses_current_part_density_mapping"]:
        issues.append("assembly center of mass uses a stale part-density mapping")
    if not checks[
        "nested_assembly_locations_use_current_transform_composition_order"
    ]:
        issues.append("nested assembly locations use a stale transform composition order")
    if not checks[
        "boolean_retained_face_names_follow_post_refine_topology_history"
    ]:
        issues.append("Boolean retained face names use topology history from before refine")
    if not checks["boolean_subshape_labels_follow_current_fillet_edge_order"]:
        issues.append("Boolean subshape labels use history from before fillet reordering")
    if not checks[
        "assembly_mates_share_current_frame_unit_and_parent_location"
    ]:
        issues.append("assembly mates mix stale frames, units, or parent locations")
    if not checks[
        "assembly_mass_properties_share_density_unit_and_part_locations"
    ]:
        issues.append(
            "assembly mass properties mix density units or part locations from stale generations"
        )
    if not checks[
        "swept_solid_uses_current_path_frames_twist_and_profile_orientation"
    ]:
        issues.append(
            "swept solid uses stale path frames, twist samples, or profile orientation"
        )
    if not checks[
        "joint_connectors_use_current_labeled_face_frame_and_parent_location"
    ]:
        issues.append("joint connector uses a stale labeled face, frame, or parent location")
    if not checks[
        "inertia_uses_current_density_unit_location_and_principal_axes"
    ]:
        issues.append("inertia properties mix stale density, unit, location, or principal axes")
    if not checks[
        "assembly_mass_properties_use_current_locations_operands_members_and_density"
    ]:
        issues.append(
            "assembly mass properties mix stale locations, Boolean operands, compound members, or density maps"
        )
    if not checks[
        "loft_volume_uses_current_spline_tolerances_and_watertight_shell"
    ]:
        issues.append(
            "loft volume belongs to a stale spline, tessellation tolerance, or non-watertight shell"
        )
    if not checks[
        "transformed_assembly_mass_properties_use_current_transforms_density_units_and_axes"
    ]:
        issues.append(
            "transformed assembly mass properties mix stale transforms, densities, units, or principal axes"
        )
    if not checks[
        "fillet_chamfer_topology_uses_current_selection_names_order_and_fingerprint"
    ]:
        issues.append(
            "fillet/chamfer topology uses stale edge selections, names, order, or build fingerprints"
        )
    if not checks[
        "boolean_result_uses_current_operands_tolerance_healing_topology_and_volume"
    ]:
        issues.append(
            "Boolean result mixes stale operands, tolerances, healing, topology, or volume"
        )
    if not checks[
        "assembly_mates_use_current_transforms_dof_solver_and_loop_closure"
    ]:
        issues.append(
            "assembly mates mix stale transforms, DOF state, solver, or loop closure"
        )
    volume_errors = [row["volume_rel_error"] or 0.0 for row in all_rows]
    area_errors = [row["area_rel_error"] or 0.0 for row in all_rows]
    bbox_errors = [row["bbox_abs_error"] or 0.0 for row in all_rows]
    ok = all(checks.values())
    return {
        "policy": "build123d_external_cad_volume_area_bbox_crosscheck",
        "status": "ok" if ok else "needs_attention",
        "ok_for_cad_roundtrip_mass_properties": ok,
        "rtol": float(rtol),
        "bbox_atol": float(bbox_atol),
        "n_reference_rows": len(reference),
        "n_sources": len(sets),
        "n_failed_rows": len(failed_rows),
        "max_volume_rel_error": max(volume_errors) if volume_errors else 0.0,
        "max_area_rel_error": max(area_errors) if area_errors else 0.0,
        "max_bbox_abs_error": max(bbox_errors) if bbox_errors else 0.0,
        "sources": [item["source"] for item in sets],
        "checks": checks,
        "issues": issues,
        "inventory": inventory,
        "comparison_sets": sets,
    }


def shape_cubit_export_package_handoff_gate(
    shape_rows,
    package_gate,
    *,
    geometry_id_key="geometry_id",
    expected_export_id=None,
    require_sidecar_inventory_counts=False,
    require_sidecar_schema=False,
) -> dict:
    """Check that build123d CAD rows hand off to the intended Cubit package.

    This is the CAD-side companion to ``cubit_export_package_identity_gate``.
    build123d keeps the geometry intent and mass properties; Cubit later owns
    the hex-led ``.vol`` package.  The two should meet on an explicit
    ``geometry_id`` instead of only a human-readable filename.
    """

    rows = [dict(row) for row in shape_rows]
    if not rows:
        raise ValueError("shape_rows must not be empty")
    if not isinstance(package_gate, dict):
        raise ValueError("package_gate must be a mapping")

    geometry_ids = [str(row.get(geometry_id_key, "")).strip() for row in rows]
    names = [str(row.get("name", "")).strip() for row in rows]
    package_geometry_ids = [str(value).strip() for value in package_gate.get("geometry_ids", [])]
    package_export_ids = [str(value).strip() for value in package_gate.get("export_ids", [])]
    package_checks = package_gate.get("checks", {})
    if not isinstance(package_checks, dict):
        package_checks = {}
    expected_export = None if expected_export_id is None else str(expected_export_id).strip()
    sidecar_counts_required = bool(require_sidecar_inventory_counts)
    sidecar_schema_required = bool(require_sidecar_schema)

    def has_bbox(row):
        bbox = row.get("bounding_box") or row.get("bbox")
        return isinstance(bbox, dict) and bbox.get("size") is not None

    checks = {
        "shape_rows_present": bool(rows),
        "shape_names_recorded": all(bool(name) for name in names),
        "shape_geometry_ids_recorded": all(bool(value) for value in geometry_ids),
        "shape_geometry_ids_unique": len(set(geometry_ids)) == len(geometry_ids),
        "shape_rows_have_volume_area_bbox": all(
            row.get("volume") is not None and row.get("area") is not None and has_bbox(row)
            for row in rows
        ),
        "package_gate_ok": package_gate.get("status") == "ok",
        "package_has_export_id": bool(package_export_ids),
        "package_has_geometry_id": bool(package_geometry_ids),
        "package_vol_sidecar_pairs_vol": package_checks.get("vol_sidecar_pairs_vol") is True,
        "package_sidecar_schema_recorded": (
            not sidecar_schema_required
            or package_checks.get("vol_sidecar_schema_id_recorded_when_required") is True
        ),
        "package_sidecar_schema_matches_expected": (
            not sidecar_schema_required
            or package_checks.get("expected_vol_sidecar_schema_id_matches") is True
        ),
        "package_sidecar_inventory_counts_recorded": (
            not sidecar_counts_required
            or package_checks.get("vol_sidecar_inventory_counts_recorded_when_required") is True
        ),
        "package_sidecar_element_count_matches_inventory": (
            not sidecar_counts_required
            or package_checks.get("vol_sidecar_element_count_matches_inventory") is True
        ),
        "package_sidecar_point_count_matches_inventory": (
            not sidecar_counts_required
            or package_checks.get("vol_sidecar_point_count_matches_inventory") is True
        ),
        "package_sidecar_order_matches_expected": (
            not sidecar_counts_required
            or package_checks.get("vol_sidecar_order_matches_expected") is True
        ),
        "package_raw_result_present": "raw_result" in package_gate.get("kinds", []),
        "geometry_id_matches_package": bool(geometry_ids)
        and bool(package_geometry_ids)
        and set(geometry_ids) == set(package_geometry_ids),
        "export_id_matches_expected": expected_export is None or set(package_export_ids) == {expected_export},
    }
    issues = []
    if not checks["shape_geometry_ids_recorded"]:
        issues.append("build123d rows need an explicit geometry_id before Cubit handoff")
    if not checks["geometry_id_matches_package"]:
        issues.append("build123d geometry_id does not match the Cubit export package")
    if not checks["package_vol_sidecar_pairs_vol"]:
        issues.append("Cubit package did not prove the .vol/.vol.json pairing")
    if sidecar_counts_required and not (
        checks["package_sidecar_element_count_matches_inventory"]
        and checks["package_sidecar_point_count_matches_inventory"]
        and checks["package_sidecar_order_matches_expected"]
    ):
        issues.append("Cubit package sidecar counts/order do not match the parsed .vol inventory")
    if sidecar_schema_required and not (
        checks["package_sidecar_schema_recorded"]
        and checks["package_sidecar_schema_matches_expected"]
    ):
        issues.append("Cubit package sidecar schema id is missing or stale")
    return {
        "policy": "build123d_cubit_export_package_handoff_gate",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "geometry_id_key": geometry_id_key,
        "shape_names": names,
        "shape_geometry_ids": sorted(set(geometry_ids)),
        "package_geometry_ids": sorted(set(package_geometry_ids)),
        "package_export_ids": sorted(set(package_export_ids)),
        "expected_export_id": expected_export,
        "require_sidecar_inventory_counts": sidecar_counts_required,
        "require_sidecar_schema": sidecar_schema_required,
        "package_policy": package_gate.get("policy"),
        "checks": checks,
        "issues": issues,
        "notes": [
            "Use after build123d mass-property rows pass and before Cubit .vol packages enter notebooks or solver-ready runs.",
            "The CAD intent row and the mesh export package should share geometry_id; filenames alone are not enough.",
            "When Cubit .vol.json schema identity is available, require the package schema id before accepting the CAD-to-mesh handoff.",
            "When Cubit .vol.json metadata is available, require sidecar element/point counts and order to match the parsed .vol inventory before CAD-to-mesh handoff.",
        ],
    }


def shape_cubit_quality_package_handoff_gate(
    shape_rows,
    quality_package_gate,
    *,
    geometry_id_key="geometry_id",
    expected_export_id=None,
    require_export_inventory=False,
) -> dict:
    """Check that build123d CAD rows match a Cubit headless quality package.

    This is the CAD-side companion to
    ``cubit_headless_batch_quality_package_gate``.  It is used when Cubit has
    already proven a headless mesh-quality run.  When a ``.vol`` package is
    being consumed, set ``require_export_inventory=True`` so the CAD row also
    demands the parsed export inventory checks from the Cubit package gate.
    The CAD intent and mesh-quality evidence still need to meet on explicit
    ``geometry_id`` and ``export_id`` values.
    """

    rows = [dict(row) for row in shape_rows]
    if not rows:
        raise ValueError("shape_rows must not be empty")
    if not isinstance(quality_package_gate, dict):
        raise ValueError("quality_package_gate must be a mapping")

    geometry_ids = [str(row.get(geometry_id_key, "")).strip() for row in rows]
    names = [str(row.get("name", "")).strip() for row in rows]
    package_geometry_id = str(quality_package_gate.get("geometry_id", "")).strip()
    package_export_id = str(quality_package_gate.get("export_id", "")).strip()
    package_checks = quality_package_gate.get("checks", {})
    if not isinstance(package_checks, dict):
        package_checks = {}
    expected_export = None if expected_export_id is None else str(expected_export_id).strip()
    inventory_required = bool(require_export_inventory)

    def has_bbox(row):
        bbox = row.get("bounding_box") or row.get("bbox")
        return isinstance(bbox, dict) and bbox.get("size") is not None

    checks = {
        "shape_rows_present": bool(rows),
        "shape_names_recorded": all(bool(name) for name in names),
        "shape_geometry_ids_recorded": all(bool(value) for value in geometry_ids),
        "shape_geometry_ids_unique": len(set(geometry_ids)) == len(geometry_ids),
        "shape_rows_have_volume_area_bbox": all(
            row.get("volume") is not None and row.get("area") is not None and has_bbox(row)
            for row in rows
        ),
        "quality_package_gate_ok": quality_package_gate.get("status") == "ok",
        "quality_package_policy_known": (
            quality_package_gate.get("policy") == "cubit_headless_batch_quality_package_gate"
        ),
        "quality_package_has_export_id": bool(package_export_id),
        "quality_package_has_geometry_id": bool(package_geometry_id),
        "quality_package_headless": package_checks.get("headless_command_recorded") is True,
        "quality_package_count_positive": package_checks.get("quality_count_positive") is True
        or int(quality_package_gate.get("quality_count", 0) or 0) > 0,
        "quality_package_export_inventory_present": (
            not inventory_required
            or package_checks.get("export_inventory_recorded") is True
        ),
        "quality_package_export_inventory_volume_positive": (
            not inventory_required
            or package_checks.get("export_inventory_volume_elements_positive") is True
        ),
        "quality_package_export_inventory_routing_ok": (
            not inventory_required
            or package_checks.get("export_inventory_routing_hint_matches_expected") is True
        ),
        "quality_package_export_inventory_count_matches": (
            not inventory_required
            or package_checks.get("export_inventory_count_matches_quality") is True
        ),
        "quality_package_export_inventory_contains_quality_element": (
            not inventory_required
            or package_checks.get("export_inventory_contains_quality_element") is True
        ),
        "quality_package_export_inventory_not_tri_tet_only_for_cubit_route": (
            not inventory_required
            or package_checks.get("export_inventory_not_tri_tet_only_for_cubit_hex_route") is True
        ),
        "geometry_id_matches_quality_package": bool(geometry_ids)
        and bool(package_geometry_id)
        and set(geometry_ids) == {package_geometry_id},
        "export_id_matches_expected": expected_export is None or package_export_id == expected_export,
    }
    issues = []
    if not checks["shape_geometry_ids_recorded"]:
        issues.append("build123d rows need an explicit geometry_id before Cubit quality handoff")
    if not checks["geometry_id_matches_quality_package"]:
        issues.append("build123d geometry_id does not match the Cubit quality package")
    if not checks["quality_package_headless"]:
        issues.append("Cubit quality package did not prove headless execution")
    if not checks["quality_package_count_positive"]:
        issues.append("Cubit quality package has no positive element-quality count")
    if inventory_required and not checks["quality_package_export_inventory_present"]:
        issues.append("Cubit quality package did not include parsed .vol inventory evidence")
    if inventory_required and not checks["quality_package_export_inventory_count_matches"]:
        issues.append("Cubit .vol inventory count does not match the quality package count")
    if inventory_required and not checks["quality_package_export_inventory_contains_quality_element"]:
        issues.append("Cubit .vol inventory does not contain the replayed quality element kind")
    if inventory_required and not checks["quality_package_export_inventory_not_tri_tet_only_for_cubit_route"]:
        issues.append("Cubit hex-led quality package is paired with tri/tet-only inventory")
    return {
        "policy": "build123d_cubit_quality_package_handoff_gate",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "geometry_id_key": geometry_id_key,
        "shape_names": names,
        "shape_geometry_ids": sorted(set(geometry_ids)),
        "quality_package_geometry_id": package_geometry_id,
        "quality_package_export_id": package_export_id,
        "expected_export_id": expected_export,
        "quality_package_policy": quality_package_gate.get("policy"),
        "quality_count": int(quality_package_gate.get("quality_count", 0) or 0),
        "require_export_inventory": inventory_required,
        "quality_package_export_inventory_source": quality_package_gate.get("export_inventory_source"),
        "quality_package_export_inventory_volume_kind_counts": quality_package_gate.get("export_inventory_volume_kind_counts"),
        "quality_package_export_inventory_routing_hint": quality_package_gate.get("export_inventory_routing_hint"),
        "quality_package_export_inventory_is_tri_tet_only": quality_package_gate.get("export_inventory_is_tri_tet_only"),
        "checks": checks,
        "issues": issues,
        "notes": [
            "Use after build123d mass-property rows pass and Cubit has a headless quality package.",
            "Mesh-quality evidence is only reusable when it matches the CAD row geometry_id and export_id.",
            "When a .vol is being consumed, require the Cubit package gate's export_inventory checks too.",
            "Cubit hex-led quality evidence must stay on the Cubit hex/mixed route and must not be paired with a tri/tet-only inventory.",
        ],
    }


def shape_cubit_quality_ledger_handoff_gate(
    shape_rows,
    quality_ledger_gate,
    *,
    geometry_id_key="geometry_id",
    expected_quality_artifact_id=None,
    expected_quality_digest=None,
    expected_metric_set_id=None,
    expected_export_id=None,
    expected_mesh_artifact_id=None,
    expected_mesh_digest=None,
    expected_route="cubit_hex_or_mixed_path",
    min_scaled_jacobian_threshold=0.2,
    require_hex_or_mixed_route=True,
    require_quality_execution_metadata=False,
) -> dict:
    """Check that build123d CAD rows match a Cubit mesh-quality ledger.

    This is the build123d-side companion to
    ``cubit_mesh_quality_ledger_identity_gate``.  It is stricter than the
    older quality-package handoff: a reusable quality row must travel with its
    own artifact id, digest, metric-set id, exported mesh id/digest, geometry
    id, and Cubit hex/mixed routing hint.  That keeps a fresh CAD package from
    accidentally reusing a stale quality JSON or a quality row from another
    mesh.
    """

    rows = [dict(row) for row in shape_rows]
    if not rows:
        raise ValueError("shape_rows must not be empty")
    if not isinstance(quality_ledger_gate, dict):
        raise ValueError("quality_ledger_gate must be a mapping")

    def as_text(value):
        return str(value or "").strip()

    def as_float(value):
        try:
            return float(value)
        except (TypeError, ValueError):
            return None

    def as_int(value):
        try:
            return int(value)
        except (TypeError, ValueError):
            return None

    def has_bbox(row):
        bbox = row.get("bounding_box") or row.get("bbox")
        return isinstance(bbox, dict) and bbox.get("size") is not None

    def counts_from(row, *names):
        for name in names:
            value = row.get(name)
            if isinstance(value, dict):
                counts = {}
                for key, count in value.items():
                    try:
                        counts[str(key).strip().lower().replace("-", "_").replace(" ", "_")] = int(count)
                    except (TypeError, ValueError):
                        counts[str(key).strip().lower().replace("-", "_").replace(" ", "_")] = -1
                return counts
        return {}

    geometry_ids = [as_text(row.get(geometry_id_key)) for row in rows]
    names = [as_text(row.get("name")) for row in rows]
    routes = {
        as_text(row.get("mesh_route", row.get("routing_hint")))
        for row in rows
        if as_text(row.get("mesh_route", row.get("routing_hint")))
    }
    ledger_checks = quality_ledger_gate.get("checks", {})
    if not isinstance(ledger_checks, dict):
        ledger_checks = {}
    quality_artifact_id = as_text(quality_ledger_gate.get("quality_artifact_id"))
    quality_digest = as_text(quality_ledger_gate.get("quality_digest"))
    metric_set_id = as_text(quality_ledger_gate.get("metric_set_id"))
    export_id = as_text(quality_ledger_gate.get("export_id"))
    ledger_geometry_id = as_text(quality_ledger_gate.get("geometry_id"))
    mesh_artifact_id = as_text(quality_ledger_gate.get("mesh_artifact_id"))
    mesh_digest = as_text(quality_ledger_gate.get("mesh_digest"))
    routing_hint = as_text(quality_ledger_gate.get("routing_hint"))
    element_counts = counts_from(quality_ledger_gate, "element_type_counts", "volume_kind_counts")
    min_j = as_float(quality_ledger_gate.get("min_scaled_jacobian"))
    negative_j = as_int(quality_ledger_gate.get("negative_jacobian_count"))
    expected_quality_artifact = (
        None if expected_quality_artifact_id is None else as_text(expected_quality_artifact_id)
    )
    expected_quality_hash = None if expected_quality_digest is None else as_text(expected_quality_digest)
    expected_metric_set = None if expected_metric_set_id is None else as_text(expected_metric_set_id)
    expected_export = None if expected_export_id is None else as_text(expected_export_id)
    expected_mesh_artifact = None if expected_mesh_artifact_id is None else as_text(expected_mesh_artifact_id)
    expected_mesh_hash = None if expected_mesh_digest is None else as_text(expected_mesh_digest)
    expected_route_text = as_text(expected_route)
    threshold = float(min_scaled_jacobian_threshold)
    if threshold <= 0.0:
        raise ValueError("min_scaled_jacobian_threshold must be > 0")
    route_required = bool(require_hex_or_mixed_route)
    tri_tet_only = bool(quality_ledger_gate.get("inventory_is_tri_tet_only")) or (
        bool(element_counts) and set(element_counts).issubset({"tet", "tri", "triangle"})
    )
    hex_or_mixed_present = any(element_counts.get(kind, 0) > 0 for kind in ("hex", "pyramid", "wedge"))
    quality_execution_checks = (
        "created_at_utc_recorded_when_required",
        "created_at_utc_parseable_when_present",
        "version_recorded_when_required",
        "elapsed_s_recorded_when_required",
        "elapsed_s_finite_nonnegative_when_present",
        "timing_breakdown_recorded_when_required",
        "timing_breakdown_has_required_stage_count",
        "timing_breakdown_values_finite_nonnegative",
        "timing_breakdown_total_within_elapsed_when_present",
    )
    quality_execution_metadata_ok = all(
        ledger_checks.get(name) is True for name in quality_execution_checks
    )

    checks = {
        "shape_rows_present": bool(rows),
        "shape_names_recorded": all(bool(name) for name in names),
        "shape_geometry_ids_recorded": all(bool(value) for value in geometry_ids),
        "shape_geometry_ids_unique": len(set(geometry_ids)) == len(geometry_ids),
        "shape_rows_have_volume_area_bbox": all(
            row.get("volume") is not None and row.get("area") is not None and has_bbox(row)
            for row in rows
        ),
        "mesh_route_recorded": (not route_required) or bool(routes),
        "mesh_route_matches_expected": (not route_required) or routes == {expected_route_text},
        "quality_ledger_gate_ok": quality_ledger_gate.get("status") == "ok",
        "quality_ledger_policy_known": (
            quality_ledger_gate.get("policy") == "cubit_mesh_quality_ledger_identity_gate"
        ),
        "quality_ledger_execution_metadata_ok": (
            (not bool(require_quality_execution_metadata)) or quality_execution_metadata_ok
        ),
        "quality_artifact_id_recorded": bool(quality_artifact_id),
        "quality_digest_recorded": bool(quality_digest),
        "metric_set_id_recorded": bool(metric_set_id),
        "export_id_recorded": bool(export_id),
        "quality_ledger_geometry_id_recorded": bool(ledger_geometry_id),
        "mesh_artifact_id_recorded": bool(mesh_artifact_id),
        "mesh_digest_recorded": bool(mesh_digest),
        "routing_hint_recorded": bool(routing_hint),
        "element_type_counts_recorded": bool(element_counts),
        "element_type_counts_positive": bool(element_counts) and all(count > 0 for count in element_counts.values()),
        "min_scaled_jacobian_recorded": min_j is not None,
        "min_scaled_jacobian_above_threshold": min_j is not None and min_j >= threshold,
        "negative_jacobian_count_recorded": negative_j is not None,
        "negative_jacobian_count_zero": negative_j == 0,
        "geometry_id_matches_quality_ledger": bool(geometry_ids)
        and bool(ledger_geometry_id)
        and set(geometry_ids) == {ledger_geometry_id},
        "routing_hint_matches_expected": (not route_required) or routing_hint == expected_route_text,
        "hex_or_mixed_volume_family_present": (not route_required) or hex_or_mixed_present,
        "not_tri_tet_only_for_cubit_quality_ledger": (not route_required) or not tri_tet_only,
        "ledger_not_tri_tet_only_check_not_failed": (
            ledger_checks.get("not_tri_tet_only_for_cubit_quality_ledger") is not False
        ),
        "expected_quality_artifact_id_matches": (
            expected_quality_artifact is None or quality_artifact_id == expected_quality_artifact
        ),
        "expected_quality_digest_matches": expected_quality_hash is None or quality_digest == expected_quality_hash,
        "expected_metric_set_id_matches": expected_metric_set is None or metric_set_id == expected_metric_set,
        "expected_export_id_matches": expected_export is None or export_id == expected_export,
        "expected_mesh_artifact_id_matches": (
            expected_mesh_artifact is None or mesh_artifact_id == expected_mesh_artifact
        ),
        "expected_mesh_digest_matches": expected_mesh_hash is None or mesh_digest == expected_mesh_hash,
    }
    issues = []
    if not checks["geometry_id_matches_quality_ledger"]:
        issues.append("build123d geometry_id does not match the Cubit quality ledger")
    if not checks["quality_ledger_gate_ok"]:
        issues.append("Cubit quality ledger identity gate is not ok")
    if not checks["quality_ledger_execution_metadata_ok"]:
        issues.append("Cubit quality ledger is missing required version/date/elapsed/timing execution metadata")
    if not checks["quality_digest_recorded"] or not checks["mesh_digest_recorded"]:
        issues.append("Cubit quality ledger must record both quality digest and mesh digest")
    if not checks["min_scaled_jacobian_above_threshold"]:
        issues.append("Cubit quality ledger minimum scaled Jacobian is below threshold")
    if not checks["negative_jacobian_count_zero"]:
        issues.append("Cubit quality ledger reports negative Jacobians")
    if not checks["mesh_route_matches_expected"] or not checks["routing_hint_matches_expected"]:
        issues.append("build123d CAD route and Cubit quality ledger route do not match")
    if route_required and not checks["not_tri_tet_only_for_cubit_quality_ledger"]:
        issues.append("Cubit quality ledger is tri/tet-only but this CAD handoff requested the Cubit hex/mixed route")
    if not checks["expected_quality_digest_matches"]:
        issues.append("Cubit quality ledger digest does not match the expected quality artifact digest")
    if not checks["expected_mesh_digest_matches"]:
        issues.append("Cubit quality ledger mesh digest does not match the expected mesh artifact digest")

    return {
        "policy": "build123d_cubit_quality_ledger_handoff_gate",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "geometry_id_key": geometry_id_key,
        "shape_names": names,
        "shape_geometry_ids": sorted(set(geometry_ids)),
        "mesh_routes": sorted(routes),
        "expected_route": expected_route_text,
        "quality_ledger_policy": quality_ledger_gate.get("policy"),
        "quality_artifact_id": quality_artifact_id,
        "quality_digest": quality_digest,
        "metric_set_id": metric_set_id,
        "export_id": export_id,
        "quality_ledger_geometry_id": ledger_geometry_id,
        "mesh_artifact_id": mesh_artifact_id,
        "mesh_digest": mesh_digest,
        "routing_hint": routing_hint,
        "min_scaled_jacobian": min_j,
        "min_scaled_jacobian_threshold": threshold,
        "negative_jacobian_count": negative_j,
        "element_type_counts": element_counts,
        "inventory_is_tri_tet_only": tri_tet_only,
        "require_hex_or_mixed_route": route_required,
        "require_quality_execution_metadata": bool(require_quality_execution_metadata),
        "expected_quality_artifact_id": expected_quality_artifact,
        "expected_quality_digest": expected_quality_hash,
        "expected_metric_set_id": expected_metric_set,
        "expected_export_id": expected_export,
        "expected_mesh_artifact_id": expected_mesh_artifact,
        "expected_mesh_digest": expected_mesh_hash,
        "checks": checks,
        "issues": issues,
        "notes": [
            "Use after build123d mass-property rows pass and Cubit has produced a mesh-quality ledger identity gate.",
            "This gate binds CAD geometry_id to quality artifact id/digest, metric-set id, exported mesh id/digest, route, and element counts.",
            "For Cubit hex-led lanes, reject tri/tet-only ledgers even when the stale quality metric itself looks healthy.",
            "When required, propagate the Cubit quality-ledger execution metadata checks before CAD-to-mesh handoff reuse.",
        ],
    }


def shape_cubit_solver_route_handoff_gate(
    shape_rows,
    solver_route_gate,
    *,
    geometry_id_key="geometry_id",
    expected_route="cubit_hex_or_mixed_path",
    expected_solver_route_package_id=None,
    expected_solver_contract_artifact_id=None,
    expected_solver_contract_digest=None,
    expected_solver_contract_path=None,
    expected_solver_route_convention_schema_id=None,
    require_solver_contract_artifact=False,
    require_solver_route_convention_schema=False,
    require_no_implicit_tetization=True,
) -> dict:
    """Check that build123d CAD rows match a Cubit solver-route manifest.

    build123d owns CAD intent, while Cubit owns the mixed hex+pyramid+tet mesh
    route.  This gate binds the CAD rows to the downstream solver-route gate so
    a STEP/measurement package cannot claim a Cubit mixed route while the
    actual solver route silently tetizes pyramids or omits transition roles.
    """

    rows = [dict(row) for row in shape_rows]
    if not rows:
        raise ValueError("shape_rows must not be empty")
    if not isinstance(solver_route_gate, dict):
        raise ValueError("solver_route_gate must be a mapping")

    def norm(value):
        return str(value or "").strip().lower().replace("-", "_").replace(" ", "_")

    geometry_ids = [str(row.get(geometry_id_key, "")).strip() for row in rows]
    names = [str(row.get("name", "")).strip() for row in rows]
    routes = {
        str(row.get("mesh_route", row.get("routing_hint", ""))).strip()
        for row in rows
        if str(row.get("mesh_route", row.get("routing_hint", ""))).strip()
    }
    roles = {
        norm(row.get("role", row.get("mesh_role", row.get("region_role", ""))))
        for row in rows
        if norm(row.get("role", row.get("mesh_role", row.get("region_role", ""))))
    }
    route_checks = solver_route_gate.get("checks", {})
    if not isinstance(route_checks, dict):
        route_checks = {}
    expected_package = (
        None
        if expected_solver_route_package_id is None
        else str(expected_solver_route_package_id).strip()
    )
    solver_route_package = str(
        solver_route_gate.get("solver_route_package_id", solver_route_gate.get("package_id", ""))
    ).strip()
    solver_contract_artifact_id = str(solver_route_gate.get("solver_contract_artifact_id", "")).strip()
    solver_contract_digest = str(solver_route_gate.get("solver_contract_digest", "")).strip()
    solver_contract_path = str(solver_route_gate.get("solver_contract_path", "")).strip()
    solver_route_convention_schema_id = str(
        solver_route_gate.get("solver_route_convention_schema_id", "")
    ).strip()
    expected_contract_id = (
        None
        if expected_solver_contract_artifact_id is None
        else str(expected_solver_contract_artifact_id).strip()
    )
    expected_contract_digest = (
        None if expected_solver_contract_digest is None else str(expected_solver_contract_digest).strip()
    )
    expected_contract_path = (
        None if expected_solver_contract_path is None else str(expected_solver_contract_path).strip()
    )
    expected_route_convention_schema = (
        None
        if expected_solver_route_convention_schema_id is None
        else str(expected_solver_route_convention_schema_id).strip()
    )
    expected_route_text = str(expected_route or "").strip()
    no_implicit_required = bool(require_no_implicit_tetization)
    solver_contract_required = bool(require_solver_contract_artifact)
    route_convention_required = bool(
        require_solver_route_convention_schema or expected_route_convention_schema is not None
    )

    checks = {
        "shape_rows_present": bool(rows),
        "shape_names_recorded": all(bool(name) for name in names),
        "shape_geometry_ids_recorded": all(bool(value) for value in geometry_ids),
        "shape_geometry_ids_unique": len(set(geometry_ids)) == len(geometry_ids),
        "mesh_route_recorded": bool(routes),
        "mesh_route_matches_expected": routes == {expected_route_text},
        "solver_route_gate_ok": solver_route_gate.get("status") == "ok",
        "solver_route_policy_known": (
            solver_route_gate.get("policy") == "cubit_mixed_solver_route_manifest_gate"
        ),
        "solver_route_package_id_recorded": bool(solver_route_package),
        "expected_solver_route_package_id_matches": (
            expected_package is None or solver_route_package == expected_package
        ),
        "solver_contract_artifact_id_recorded_when_required": (
            not solver_contract_required
            or (
                bool(solver_contract_artifact_id)
                and route_checks.get("solver_contract_artifact_id_recorded_when_required") is not False
            )
        ),
        "solver_contract_digest_recorded_when_required": (
            not solver_contract_required
            or (
                bool(solver_contract_digest)
                and route_checks.get("solver_contract_digest_recorded_when_required") is not False
            )
        ),
        "solver_contract_path_recorded_when_required": (
            not solver_contract_required
            or (
                bool(solver_contract_path)
                and route_checks.get("solver_contract_path_recorded_when_required") is not False
            )
        ),
        "expected_solver_contract_artifact_id_matches": (
            expected_contract_id is None
            or (
                solver_contract_artifact_id == expected_contract_id
                and route_checks.get("expected_solver_contract_artifact_id_matches") is not False
            )
        ),
        "expected_solver_contract_digest_matches": (
            expected_contract_digest is None
            or (
                solver_contract_digest == expected_contract_digest
                and route_checks.get("expected_solver_contract_digest_matches") is not False
            )
        ),
        "expected_solver_contract_path_matches": (
            expected_contract_path is None
            or (
                solver_contract_path == expected_contract_path
                and route_checks.get("expected_solver_contract_path_matches") is not False
            )
        ),
        "solver_route_convention_schema_id_recorded_when_required": (
            not route_convention_required
            or (
                bool(solver_route_convention_schema_id)
                and route_checks.get("solver_route_convention_schema_id_recorded_when_required") is not False
            )
        ),
        "expected_solver_route_convention_schema_id_matches": (
            expected_route_convention_schema is None
            or (
                solver_route_convention_schema_id == expected_route_convention_schema
                and route_checks.get("expected_solver_route_convention_schema_id_matches") is not False
            )
        ),
        "solver_route_routing_hint_matches_expected": (
            solver_route_gate.get("routing_hint") == expected_route_text
        ),
        "solver_route_volume_routes_cover_inventory": (
            route_checks.get("volume_route_kinds_cover_inventory") is True
        ),
        "solver_route_surface_routes_cover_inventory": (
            route_checks.get("surface_route_kinds_cover_inventory") is True
        ),
        "solver_route_pyramid_transition_role_recorded": (
            route_checks.get("pyramid_transition_role_recorded") is True
        ),
        "solver_route_no_implicit_tetization": (
            not no_implicit_required
            or route_checks.get("no_implicit_tetization_recorded") is True
        ),
        "solver_route_tet_only_owner_is_netgen": (
            route_checks.get("tet_only_owner_is_netgen") is True
        ),
        "cad_role_not_tet_only": "tet_only" not in roles,
    }
    issues = []
    if not checks["mesh_route_matches_expected"]:
        issues.append("build123d mesh_route does not match the expected Cubit mixed route")
    if not checks["solver_route_gate_ok"]:
        issues.append("Cubit solver-route manifest gate is not ok")
    if not checks["solver_route_pyramid_transition_role_recorded"]:
        issues.append("Cubit solver route did not keep pyramid cells as transition bridge cells")
    if not checks["solver_route_no_implicit_tetization"]:
        issues.append("Cubit solver route allows implicit tetization")
    if not checks["solver_contract_path_recorded_when_required"]:
        issues.append("Cubit solver route did not record the downstream solver contract path")
    if not checks["expected_solver_contract_digest_matches"]:
        issues.append("Cubit solver route contract digest does not match the expected downstream reader contract")
    if not checks["solver_route_convention_schema_id_recorded_when_required"]:
        issues.append("Cubit solver route did not record the solver-route convention schema id")
    if not checks["expected_solver_route_convention_schema_id_matches"]:
        issues.append("Cubit solver route convention schema id does not match the expected route convention")
    return {
        "policy": "build123d_cubit_solver_route_handoff_gate",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "geometry_id_key": geometry_id_key,
        "shape_names": names,
        "shape_geometry_ids": sorted(set(geometry_ids)),
        "mesh_routes": sorted(routes),
        "expected_route": expected_route_text,
        "cad_roles": sorted(roles),
        "solver_route_policy": solver_route_gate.get("policy"),
        "solver_route_package_id": solver_route_package,
        "expected_solver_route_package_id": expected_package,
        "solver_contract_artifact_id": solver_contract_artifact_id,
        "solver_contract_digest": solver_contract_digest,
        "solver_contract_path": solver_contract_path,
        "solver_route_convention_schema_id": solver_route_convention_schema_id,
        "expected_solver_contract_artifact_id": expected_contract_id,
        "expected_solver_contract_digest": expected_contract_digest,
        "expected_solver_contract_path": expected_contract_path,
        "expected_solver_route_convention_schema_id": expected_route_convention_schema,
        "require_solver_contract_artifact": solver_contract_required,
        "require_solver_route_convention_schema": route_convention_required,
        "require_no_implicit_tetization": no_implicit_required,
        "checks": checks,
        "issues": issues,
        "notes": [
            "Use after build123d CAD rows request a Cubit mixed route and before Cubit/solver notebooks consume the package.",
            "The CAD row owns mesh intent; the Cubit route gate owns hex primary, pyramid transition, tet compatibility, and no implicit tetization.",
            "When the Cubit route is promoted as solver-ready, propagate the downstream solver/reader contract artifact id, digest, and path from the Cubit route gate.",
            "When the route convention is part of the claim, propagate solver_route_convention_schema_id so CAD handoffs do not reuse value-only mixed-route manifests.",
        ],
    }


def shape_cad_handoff_manifest_gate(
    shape_rows,
    *,
    file_manifest=(),
    external_volume_summary=None,
    cubit_export_handoff=None,
    cubit_quality_handoff=None,
    cubit_quality_ledger_handoff=None,
    cubit_solver_route_handoff=None,
    cubit_meshing_scheme_handoff=None,
    required_file_kinds=("step", "build123d_measurement_json"),
    expected_geometry_ids=None,
    expected_cad_output_artifact_id=None,
    expected_cad_output_digest=None,
    require_cad_output_artifact=False,
    expected_cad_observable_id=None,
    expected_cad_observable_family=None,
    require_cad_observable=False,
    expected_length_unit=None,
    expected_area_unit=None,
    expected_volume_unit=None,
    expected_measurement_convention=None,
    expected_measurement_postprocess_row_convention_schema_id=None,
    expected_measurement_component_basis_schema_id=None,
    require_measurement_postprocess_row_convention_schema=False,
    require_measurement_component_basis_schema=False,
) -> dict:
    """Check the final build123d CAD handoff manifest before meshing/solver work.

    build123d is the geometry authoring lane.  Before a shape package is handed
    to Cubit, CST, or a solver notebook, keep the build123d mass-property rows,
    exported file list, external CAD volume check, and optional Cubit package
    handoff gates together.  This prevents a good volume row from travelling
    with the wrong STEP file, stale mesh-quality package, or stale mixed solver
    route manifest.
    """

    rows = [dict(row) for row in shape_rows]
    files = [dict(item) for item in (file_manifest or ())]
    if not rows:
        raise ValueError("shape_rows must not be empty")

    def has_bbox(row):
        bbox = row.get("bounding_box") or row.get("bbox")
        return isinstance(bbox, dict) and bbox.get("size") is not None

    def collect_file_values(*names):
        values = []
        for item in files:
            for name in names:
                if name in item and item[name] is not None:
                    text = str(item[name]).strip()
                    if text:
                        values.append(text)
        return values

    def unique_strings(values):
        return sorted(set(values))

    def norm_unit(value):
        return str(value or "").strip().lower().replace(" ", "")

    def norm_label(value):
        return str(value or "").strip().lower().replace("-", "_").replace(" ", "_")

    def collect_row_values(*names):
        values = []
        for row in rows:
            units = row.get("units")
            for name in names:
                value = row.get(name)
                if value is not None:
                    text = str(value).strip()
                    if text:
                        values.append(text)
            if isinstance(units, dict):
                for name in names:
                    unit_key = {
                        "length_unit": "length",
                        "unit_length": "length",
                        "area_unit": "area",
                        "unit_area": "area",
                        "volume_unit": "volume",
                        "unit_volume": "volume",
                    }.get(name)
                    for key in (name, unit_key):
                        if not key:
                            continue
                        value = units.get(key)
                        if value is not None:
                            text = str(value).strip()
                            if text:
                                values.append(text)
        return values

    geometry_ids = [str(row.get("geometry_id", "")).strip() for row in rows]
    names = [str(row.get("name", "")).strip() for row in rows]
    expected_ids = None if expected_geometry_ids is None else {
        str(value).strip() for value in expected_geometry_ids if str(value).strip()
    }
    file_kinds = [str(item.get("kind", "")).strip() for item in files]
    file_paths = [str(item.get("path", "")).strip() for item in files]
    required_kinds = [str(kind).strip() for kind in required_file_kinds if str(kind).strip()]
    cad_output_artifact_ids = unique_strings(collect_file_values(
        "cad_output_artifact_id",
        "step_output_artifact_id",
        "export_output_artifact_id",
        "output_artifact_id",
    ))
    cad_output_digests = unique_strings(collect_file_values(
        "cad_output_digest",
        "step_output_digest",
        "export_output_digest",
        "output_digest",
        "cad_output_sha256",
        "step_sha256",
    ))
    cad_output_paths = unique_strings(collect_file_values(
        "cad_output_path",
        "step_output_path",
        "export_output_path",
        "output_path",
    ))
    cad_observable_ids = unique_strings(collect_file_values(
        "cad_observable_id",
        "handoff_observable_id",
        "output_observable_id",
        "observable_id",
    ))
    cad_observable_families = unique_strings(collect_file_values(
        "cad_observable_family",
        "handoff_observable_family",
        "output_observable_family",
        "observable_family",
    ))
    cad_output_artifact_id = cad_output_artifact_ids[0] if cad_output_artifact_ids else ""
    cad_output_digest = cad_output_digests[0] if cad_output_digests else ""
    cad_output_path = cad_output_paths[0] if cad_output_paths else ""
    cad_observable_id = cad_observable_ids[0] if cad_observable_ids else ""
    cad_observable_family = cad_observable_families[0] if cad_observable_families else ""
    length_units = unique_strings(
        norm_unit(value)
        for value in (
            collect_row_values("length_unit", "unit_length")
            + collect_file_values("length_unit", "unit_length")
        )
        if str(value).strip()
    )
    area_units = unique_strings(
        norm_unit(value)
        for value in (
            collect_row_values("area_unit", "unit_area")
            + collect_file_values("area_unit", "unit_area")
        )
        if str(value).strip()
    )
    volume_units = unique_strings(
        norm_unit(value)
        for value in (
            collect_row_values("volume_unit", "unit_volume")
            + collect_file_values("volume_unit", "unit_volume")
        )
        if str(value).strip()
    )
    measurement_conventions = unique_strings(
        norm_label(value)
        for value in (
            collect_row_values(
                "cad_measurement_convention",
                "measurement_convention",
                "mass_property_convention",
            )
            + collect_file_values(
                "cad_measurement_convention",
                "measurement_convention",
                "mass_property_convention",
            )
        )
        if str(value).strip()
    )
    measurement_postprocess_row_convention_schema_ids = unique_strings(
        value
        for value in (
            collect_row_values(
                "cad_measurement_postprocess_row_convention_schema_id",
                "cadMeasurementPostprocessRowConventionSchemaId",
                "measurement_postprocess_row_convention_schema_id",
                "measurementPostprocessRowConventionSchemaId",
                "postprocess_row_convention_schema_id",
                "postprocessRowConventionSchemaId",
            )
            + collect_file_values(
                "cad_measurement_postprocess_row_convention_schema_id",
                "cadMeasurementPostprocessRowConventionSchemaId",
                "measurement_postprocess_row_convention_schema_id",
                "measurementPostprocessRowConventionSchemaId",
                "postprocess_row_convention_schema_id",
                "postprocessRowConventionSchemaId",
            )
        )
        if str(value).strip()
    )
    measurement_component_basis_schema_ids = unique_strings(
        value
        for value in (
            collect_row_values(
                "cad_measurement_component_basis_schema_id",
                "cadMeasurementComponentBasisSchemaId",
                "measurement_component_basis_schema_id",
                "measurementComponentBasisSchemaId",
                "mass_property_component_basis_schema_id",
                "massPropertyComponentBasisSchemaId",
                "component_basis_schema_id",
                "componentBasisSchemaId",
            )
            + collect_file_values(
                "cad_measurement_component_basis_schema_id",
                "cadMeasurementComponentBasisSchemaId",
                "measurement_component_basis_schema_id",
                "measurementComponentBasisSchemaId",
                "mass_property_component_basis_schema_id",
                "massPropertyComponentBasisSchemaId",
                "component_basis_schema_id",
                "componentBasisSchemaId",
            )
        )
        if str(value).strip()
    )
    length_unit = length_units[0] if length_units else ""
    area_unit = area_units[0] if area_units else ""
    volume_unit = volume_units[0] if volume_units else ""
    measurement_convention = measurement_conventions[0] if measurement_conventions else ""
    measurement_postprocess_row_convention_schema_id = (
        measurement_postprocess_row_convention_schema_ids[0]
        if measurement_postprocess_row_convention_schema_ids
        else ""
    )
    measurement_component_basis_schema_id = (
        measurement_component_basis_schema_ids[0]
        if measurement_component_basis_schema_ids
        else ""
    )
    expected_cad_output_artifact = (
        None if expected_cad_output_artifact_id is None else str(expected_cad_output_artifact_id).strip()
    )
    expected_cad_output_hash = (
        None if expected_cad_output_digest is None else str(expected_cad_output_digest).strip()
    )
    expected_cad_observable = (
        None if expected_cad_observable_id is None else str(expected_cad_observable_id).strip()
    )
    expected_cad_observable_kind = (
        None if expected_cad_observable_family is None else str(expected_cad_observable_family).strip()
    )
    expected_length_unit_norm = norm_unit(expected_length_unit)
    expected_area_unit_norm = norm_unit(expected_area_unit)
    expected_volume_unit_norm = norm_unit(expected_volume_unit)
    expected_measurement_convention_norm = norm_label(expected_measurement_convention)
    expected_measurement_postprocess_row_convention_schema = (
        ""
        if expected_measurement_postprocess_row_convention_schema_id is None
        else str(expected_measurement_postprocess_row_convention_schema_id).strip()
    )
    expected_measurement_component_basis_schema = (
        ""
        if expected_measurement_component_basis_schema_id is None
        else str(expected_measurement_component_basis_schema_id).strip()
    )
    measurement_postprocess_row_convention_schema_required = bool(
        require_measurement_postprocess_row_convention_schema
        or expected_measurement_postprocess_row_convention_schema
    )
    measurement_component_basis_schema_required = bool(
        require_measurement_component_basis_schema
        or expected_measurement_component_basis_schema
    )

    volume_summary = dict(external_volume_summary or {})
    export_handoff = dict(cubit_export_handoff or {})
    quality_handoff = dict(cubit_quality_handoff or {})
    quality_ledger_handoff = dict(cubit_quality_ledger_handoff or {})
    solver_route_handoff = dict(cubit_solver_route_handoff or {})
    meshing_scheme_handoff = dict(cubit_meshing_scheme_handoff or {})
    checks = {
        "shape_rows_present": bool(rows),
        "shape_names_recorded": all(bool(name) for name in names),
        "shape_geometry_ids_recorded": all(bool(value) for value in geometry_ids),
        "shape_geometry_ids_unique": len(set(geometry_ids)) == len(geometry_ids),
        "shape_rows_have_volume_area_bbox": all(
            row.get("volume") is not None and row.get("area") is not None and has_bbox(row)
            for row in rows
        ),
        "geometry_ids_match_expected": expected_ids is None or set(geometry_ids) == expected_ids,
        "file_manifest_present": bool(files),
        "file_kinds_recorded": all(bool(kind) for kind in file_kinds),
        "file_paths_recorded": all(bool(path) for path in file_paths),
        "required_file_kinds_present": all(kind in file_kinds for kind in required_kinds),
        "external_volume_summary_ok": volume_summary.get("status") == "ok",
        "external_volume_sources_recorded": bool(volume_summary.get("sources")),
        "cubit_export_handoff_ok": not export_handoff or export_handoff.get("status") == "ok",
        "cubit_quality_handoff_ok": not quality_handoff or quality_handoff.get("status") == "ok",
        "cubit_quality_ledger_handoff_ok": (
            not quality_ledger_handoff or quality_ledger_handoff.get("status") == "ok"
        ),
        "cubit_solver_route_handoff_ok": (
            not solver_route_handoff or solver_route_handoff.get("status") == "ok"
        ),
        "cubit_meshing_scheme_handoff_ok": (
            not meshing_scheme_handoff or meshing_scheme_handoff.get("status") == "ok"
        ),
        "cad_output_artifact_id_consistent_when_present": len(cad_output_artifact_ids) <= 1,
        "cad_output_digest_consistent_when_present": len(cad_output_digests) <= 1,
        "cad_output_path_consistent_when_present": len(cad_output_paths) <= 1,
        "cad_output_artifact_id_recorded_when_required": (
            not require_cad_output_artifact or bool(cad_output_artifact_id)
        ),
        "cad_output_digest_recorded_when_required": (
            not require_cad_output_artifact or bool(cad_output_digest)
        ),
        "cad_output_path_recorded_when_required": (
            not require_cad_output_artifact or bool(cad_output_path)
        ),
        "cad_output_artifact_id_recorded_when_expected": (
            expected_cad_output_artifact is None or bool(cad_output_artifact_id)
        ),
        "expected_cad_output_artifact_id_matches": (
            expected_cad_output_artifact is None or cad_output_artifact_id == expected_cad_output_artifact
        ),
        "cad_output_digest_recorded_when_expected": (
            expected_cad_output_hash is None or bool(cad_output_digest)
        ),
        "expected_cad_output_digest_matches": (
            expected_cad_output_hash is None or cad_output_digest == expected_cad_output_hash
        ),
        "cad_observable_id_consistent_when_present": len(cad_observable_ids) <= 1,
        "cad_observable_family_consistent_when_present": len(cad_observable_families) <= 1,
        "cad_unit_metadata_consistent_when_present": (
            len(length_units) <= 1 and len(area_units) <= 1 and len(volume_units) <= 1
        ),
        "cad_measurement_convention_consistent_when_present": len(measurement_conventions) <= 1,
        "cad_measurement_postprocess_row_convention_schema_id_consistent_when_present": (
            len(measurement_postprocess_row_convention_schema_ids) <= 1
        ),
        "cad_measurement_component_basis_schema_id_consistent_when_present": (
            len(measurement_component_basis_schema_ids) <= 1
        ),
        "cad_observable_id_recorded_when_required": (
            not require_cad_observable or bool(cad_observable_id)
        ),
        "cad_observable_family_recorded_when_required": (
            not require_cad_observable or bool(cad_observable_family)
        ),
        "cad_observable_id_recorded_when_expected": (
            expected_cad_observable is None or bool(cad_observable_id)
        ),
        "expected_cad_observable_id_matches": (
            expected_cad_observable is None or cad_observable_id == expected_cad_observable
        ),
        "cad_observable_family_recorded_when_expected": (
            expected_cad_observable_kind is None or bool(cad_observable_family)
        ),
        "expected_cad_observable_family_matches": (
            expected_cad_observable_kind is None or cad_observable_family == expected_cad_observable_kind
        ),
        "cad_length_unit_recorded_when_expected": (
            not expected_length_unit_norm or bool(length_unit)
        ),
        "expected_cad_length_unit_matches": (
            not expected_length_unit_norm or length_units == [expected_length_unit_norm]
        ),
        "cad_area_unit_recorded_when_expected": (
            not expected_area_unit_norm or bool(area_unit)
        ),
        "expected_cad_area_unit_matches": (
            not expected_area_unit_norm or area_units == [expected_area_unit_norm]
        ),
        "cad_volume_unit_recorded_when_expected": (
            not expected_volume_unit_norm or bool(volume_unit)
        ),
        "expected_cad_volume_unit_matches": (
            not expected_volume_unit_norm or volume_units == [expected_volume_unit_norm]
        ),
        "cad_measurement_convention_recorded_when_expected": (
            not expected_measurement_convention_norm or bool(measurement_convention)
        ),
        "expected_cad_measurement_convention_matches": (
            not expected_measurement_convention_norm
            or measurement_conventions == [expected_measurement_convention_norm]
        ),
        "cad_measurement_postprocess_row_convention_schema_id_recorded_when_required": (
            not measurement_postprocess_row_convention_schema_required
            or bool(measurement_postprocess_row_convention_schema_id)
        ),
        "cad_measurement_postprocess_row_convention_schema_id_recorded_when_expected": (
            not expected_measurement_postprocess_row_convention_schema
            or bool(measurement_postprocess_row_convention_schema_id)
        ),
        "expected_cad_measurement_postprocess_row_convention_schema_id_matches": (
            not expected_measurement_postprocess_row_convention_schema
            or measurement_postprocess_row_convention_schema_ids
            == [expected_measurement_postprocess_row_convention_schema]
        ),
        "cad_measurement_component_basis_schema_id_recorded_when_required": (
            not measurement_component_basis_schema_required
            or bool(measurement_component_basis_schema_id)
        ),
        "cad_measurement_component_basis_schema_id_recorded_when_expected": (
            not expected_measurement_component_basis_schema
            or bool(measurement_component_basis_schema_id)
        ),
        "expected_cad_measurement_component_basis_schema_id_matches": (
            not expected_measurement_component_basis_schema
            or measurement_component_basis_schema_ids
            == [expected_measurement_component_basis_schema]
        ),
    }
    issues = []
    if not checks["shape_rows_have_volume_area_bbox"]:
        issues.append("shape rows must carry volume, area, and bounding_box before CAD handoff")
    if not checks["required_file_kinds_present"]:
        issues.append("file manifest is missing one or more required export/result kinds")
    if not checks["external_volume_summary_ok"]:
        issues.append("external CAD volume crosscheck is missing or not ok")
    if not checks["cubit_export_handoff_ok"]:
        issues.append("Cubit export package handoff is not ok")
    if not checks["cubit_quality_handoff_ok"]:
        issues.append("Cubit quality package handoff is not ok")
    if not checks["cubit_quality_ledger_handoff_ok"]:
        issues.append("Cubit quality ledger handoff is not ok")
    if not checks["cubit_solver_route_handoff_ok"]:
        issues.append("Cubit solver-route handoff is not ok")
    if not checks["cubit_meshing_scheme_handoff_ok"]:
        issues.append("Cubit meshing-scheme handoff is not ok")
    if (
        expected_length_unit_norm
        and checks["expected_cad_length_unit_matches"] is False
        or expected_area_unit_norm
        and checks["expected_cad_area_unit_matches"] is False
        or expected_volume_unit_norm
        and checks["expected_cad_volume_unit_matches"] is False
    ):
        issues.append("CAD handoff rows/files must carry the expected length/area/volume units")
    if (
        expected_measurement_convention_norm
        and checks["expected_cad_measurement_convention_matches"] is False
    ):
        issues.append("CAD handoff rows/files must carry the expected mass-property measurement convention")
    if (
        expected_measurement_postprocess_row_convention_schema
        and checks["expected_cad_measurement_postprocess_row_convention_schema_id_matches"] is False
    ):
        issues.append("CAD handoff rows/files must carry the expected mass-property postprocess-row convention schema")
    if (
        expected_measurement_component_basis_schema
        and checks["expected_cad_measurement_component_basis_schema_id_matches"] is False
    ):
        issues.append("CAD handoff rows/files must carry the expected mass-property component-basis schema")
    return {
        "policy": "build123d_cad_handoff_manifest_gate",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "shape_names": names,
        "geometry_ids": sorted(set(geometry_ids)),
        "expected_geometry_ids": None if expected_ids is None else sorted(expected_ids),
        "file_kinds": file_kinds,
        "file_paths": file_paths,
        "required_file_kinds": required_kinds,
        "external_volume_policy": volume_summary.get("policy"),
        "external_volume_sources": volume_summary.get("sources", []),
        "cubit_export_handoff_policy": export_handoff.get("policy"),
        "cubit_quality_handoff_policy": quality_handoff.get("policy"),
        "cubit_quality_ledger_handoff_policy": quality_ledger_handoff.get("policy"),
        "cubit_solver_route_handoff_policy": solver_route_handoff.get("policy"),
        "cubit_meshing_scheme_handoff_policy": meshing_scheme_handoff.get("policy"),
        "cad_output_artifact_id": cad_output_artifact_id or None,
        "cad_output_digest": cad_output_digest or None,
        "cad_output_path": cad_output_path or None,
        "cad_output_artifact_ids": cad_output_artifact_ids,
        "cad_output_digests": cad_output_digests,
        "cad_output_paths": cad_output_paths,
        "expected_cad_output_artifact_id": expected_cad_output_artifact,
        "expected_cad_output_digest": expected_cad_output_hash,
        "require_cad_output_artifact": bool(require_cad_output_artifact),
        "cad_observable_id": cad_observable_id or None,
        "cad_observable_family": cad_observable_family or None,
        "cad_observable_ids": cad_observable_ids,
        "cad_observable_families": cad_observable_families,
        "expected_cad_observable_id": expected_cad_observable,
        "expected_cad_observable_family": expected_cad_observable_kind,
        "require_cad_observable": bool(require_cad_observable),
        "units": {"length": length_units, "area": area_units, "volume": volume_units},
        "expected_units": {
            "length": expected_length_unit_norm or None,
            "area": expected_area_unit_norm or None,
            "volume": expected_volume_unit_norm or None,
        },
        "cad_measurement_convention": measurement_convention or None,
        "cad_measurement_conventions": measurement_conventions,
        "expected_cad_measurement_convention": expected_measurement_convention_norm or None,
        "cad_measurement_postprocess_row_convention_schema_id": (
            measurement_postprocess_row_convention_schema_id or None
        ),
        "cad_measurement_postprocess_row_convention_schema_ids": (
            measurement_postprocess_row_convention_schema_ids
        ),
        "expected_cad_measurement_postprocess_row_convention_schema_id": (
            expected_measurement_postprocess_row_convention_schema or None
        ),
        "require_measurement_postprocess_row_convention_schema": (
            measurement_postprocess_row_convention_schema_required
        ),
        "cad_measurement_component_basis_schema_id": (
            measurement_component_basis_schema_id or None
        ),
        "cad_measurement_component_basis_schema_ids": (
            measurement_component_basis_schema_ids
        ),
        "expected_cad_measurement_component_basis_schema_id": (
            expected_measurement_component_basis_schema or None
        ),
        "require_measurement_component_basis_schema": (
            measurement_component_basis_schema_required
        ),
        "checks": checks,
        "issues": issues,
        "notes": [
            "Use this as the last build123d-side preflight before Cubit/CST/solver notebooks consume CAD exports.",
            "Volume is the common CAD currency, but the manifest must also keep area, bbox, file identity, and downstream package gates together.",
            "When a STEP or manifest package is consumed downstream, bind the CAD output artifact id, digest, and path explicitly.",
            "Bind the CAD observable id/family too, so volume/area/bbox manifests are not confused with mesh, field, or solver-ready observables.",
            "Bind length/area/volume units and mass-property measurement convention so Cubit/CST/build123d volume cross-checks cannot mix model-unit, mm, and SI interpretations.",
            "Bind the mass-property postprocess-row convention schema so volume/area/bbox rows cannot silently switch selected solids, compound aggregation, or objective reduction semantics.",
            "Bind the mass-property component-basis schema so volume, area, bbox, center, and derived components cannot be reinterpreted as one scalar value.",
            "When a Cubit mixed route is requested, bind the Cubit solver-route handoff too so pyramid transition cells cannot be silently tetized downstream.",
            "When Cubit scheme intent is available, bind the meshing-scheme handoff too so a fresh CAD package cannot reuse a stale exported .vol artifact.",
            "When Cubit quality-ledger identity is available, bind it too so a CAD package cannot reuse stale mesh-quality digests.",
        ],
    }


def shape_submodel_cad_handoff_gate(
    shape_rows,
    *,
    recipe_id,
    parent_model_id,
    submodel_region_id,
    crop_box,
    export_id="",
    unit="",
    file_manifest=(),
    boundary_handoff=None,
    transition_handoff=None,
    expected_geometry_ids=None,
    bbox_atol=1.0e-12,
) -> dict:
    """Check build123d local CAD identity before submodel mesh handoff.

    This is the CAD-side companion to solver/Cubit submodel boundary-handoff
    gates.  It keeps the local recipe, crop box, shape mass properties, export
    files, and optional downstream boundary-handoff gate together.
    """

    rows = [dict(row) for row in shape_rows]
    if not rows:
        raise ValueError("shape_rows must not be empty")
    recipe_text = str(recipe_id or "").strip()
    parent_text = str(parent_model_id or "").strip()
    submodel_text = str(submodel_region_id or "").strip()
    export_text = str(export_id or "").strip()
    unit_text = str(unit or "").strip()
    files = [dict(item) for item in (file_manifest or ())]
    boundary = dict(boundary_handoff or {})
    transition = dict(transition_handoff or {})
    expected_ids = None if expected_geometry_ids is None else {
        str(value).strip() for value in expected_geometry_ids if str(value).strip()
    }

    def has_bbox(row):
        bbox = row.get("bounding_box") or row.get("bbox")
        return isinstance(bbox, dict) and bbox.get("min") is not None and bbox.get("max") is not None

    def normalize_box(box):
        data = dict(box or {})
        if data.get("min") is not None and data.get("max") is not None:
            mins = [float(v) for v in data["min"]]
            maxs = [float(v) for v in data["max"]]
        elif data.get("center") is not None and data.get("size") is not None:
            center = [float(v) for v in data["center"]]
            size = [float(v) for v in data["size"]]
            mins = [center[i] - 0.5 * size[i] for i in range(3)]
            maxs = [center[i] + 0.5 * size[i] for i in range(3)]
        else:
            return None
        if len(mins) != 3 or len(maxs) != 3:
            return None
        size = [maxs[i] - mins[i] for i in range(3)]
        if any(not math.isfinite(v) for v in mins + maxs + size):
            return None
        if any(v < 0.0 for v in size):
            return None
        center = [(mins[i] + maxs[i]) * 0.5 for i in range(3)]
        return {"min": mins, "max": maxs, "size": size, "center": center}

    crop = normalize_box(crop_box)
    names = [str(row.get("name", "")).strip() for row in rows]
    geometry_ids = [str(row.get("geometry_id", "")).strip() for row in rows]
    row_bboxes = [_row_bounding_box(row) if has_bbox(row) else None for row in rows]
    tol = float(bbox_atol)

    def row_inside_crop(bbox):
        if crop is None or bbox is None:
            return False
        return all(
            bbox["min"][i] >= crop["min"][i] - tol
            and bbox["max"][i] <= crop["max"][i] + tol
            for i in range(3)
        )

    file_kinds = [str(item.get("kind", "")).strip() for item in files]
    file_paths = [str(item.get("path", "")).strip() for item in files]
    boundary_checks = boundary.get("checks", {})
    if not isinstance(boundary_checks, dict):
        boundary_checks = {}
    boundary_submodel = str(boundary.get("submodel_region_id", "") or "").strip()
    transition_checks = transition.get("checks", {})
    if not isinstance(transition_checks, dict):
        transition_checks = {}
    boundary_volume_counts = boundary.get("volume_kind_counts", {})
    if not isinstance(boundary_volume_counts, dict):
        boundary_volume_counts = {}
    boundary_surface_counts = boundary.get("surface_kind_counts", {})
    if not isinstance(boundary_surface_counts, dict):
        boundary_surface_counts = {}
    boundary_has_pyramid = int(boundary_volume_counts.get("pyramid", 0) or 0) > 0
    boundary_surface_kinds = {
        str(kind).strip()
        for kind, count in boundary_surface_counts.items()
        if str(kind).strip() and int(count or 0) > 0
    }

    def as_name_set(value):
        if value is None:
            return set()
        if isinstance(value, str):
            return {
                part.strip()
                for part in value.replace(";", ",").split(",")
                if part.strip()
            }
        if isinstance(value, dict):
            return {str(item).strip() for item in value.values() if str(item).strip()}
        return {str(item).strip() for item in value if str(item).strip()}

    boundary_material_names = as_name_set(
        boundary.get("material_names")
        or boundary.get("materials")
        or boundary.get("expected_material_names")
        or boundary.get("sidecar_material_names")
        or boundary.get("row_names")
    )
    boundary_allowed_zero_material_names = as_name_set(
        boundary.get("allowed_zero_measurement_names")
        or boundary.get("allowed_zero_material_names")
        or boundary.get("zero_volume_material_names")
    )
    transition_kinds_raw = transition.get("transition_kinds", ())
    if isinstance(transition_kinds_raw, str):
        transition_kinds = {transition_kinds_raw.strip()}
    else:
        transition_kinds = {
            str(value).strip() for value in transition_kinds_raw if str(value).strip()
        }
    transition_surface_kinds_raw = (
        transition.get("surface_kinds")
        or transition.get("transition_surface_kinds")
        or transition.get("required_surface_kinds")
        or ()
    )
    if isinstance(transition_surface_kinds_raw, str):
        transition_surface_kinds = {
            part.strip()
            for part in transition_surface_kinds_raw.replace(";", ",").split(",")
            if part.strip()
        }
    else:
        transition_surface_kinds = {
            str(value).strip()
            for value in transition_surface_kinds_raw
            if str(value).strip()
        }
    transition_material_names = as_name_set(
        transition.get("downstream_material_names")
        or transition.get("material_names")
    )
    transition_allowed_zero_material_names = as_name_set(
        transition.get("allowed_zero_downstream_material_names")
        or transition.get("allowed_zero_material_names")
    )
    boundary_interface_roles = as_name_set(
        boundary.get("interface_roles")
        or boundary.get("roles_present")
        or boundary.get("required_roles")
        or boundary.get("interface_adjacency_roles")
    )
    transition_interface_roles = as_name_set(
        transition.get("interface_roles")
        or transition.get("required_interface_roles")
        or transition.get("transition_interface_roles")
    )

    checks = {
        "shape_rows_present": bool(rows),
        "shape_names_recorded": all(bool(name) for name in names),
        "shape_geometry_ids_recorded": all(bool(value) for value in geometry_ids),
        "shape_geometry_ids_unique": len(set(geometry_ids)) == len(geometry_ids),
        "geometry_ids_match_expected": expected_ids is None or set(geometry_ids) == expected_ids,
        "shape_rows_have_volume_area_bbox": all(
            row.get("volume") is not None and row.get("area") is not None and has_bbox(row)
            for row in rows
        ),
        "recipe_id_recorded": bool(recipe_text),
        "parent_model_id_recorded": bool(parent_text),
        "submodel_region_id_recorded": bool(submodel_text),
        "crop_box_recorded": crop is not None,
        "shape_bboxes_inside_crop": bool(row_bboxes) and all(row_inside_crop(bbox) for bbox in row_bboxes),
        "export_id_recorded": bool(export_text),
        "unit_recorded": bool(unit_text),
        "file_manifest_present": bool(files),
        "file_kinds_recorded": all(bool(kind) for kind in file_kinds),
        "file_paths_recorded": all(bool(path) for path in file_paths),
        "step_file_present": "step" in file_kinds,
        "measurement_json_present": "build123d_measurement_json" in file_kinds,
        "boundary_handoff_ok": (
            not boundary
            or (
                boundary.get("status") == "ok"
                and "boundary_handoff" in str(boundary.get("policy", ""))
            )
        ),
        "boundary_handoff_submodel_matches": (
            not boundary
            or (bool(boundary_submodel) and boundary_submodel == submodel_text)
        ),
        "transition_handoff_ok": (
            not transition
            or (
                transition.get("status") == "ok"
                and transition.get("policy") == "build123d_hex_tet_transition_role_metadata_gate"
            )
        ),
    }
    if boundary:
        checks["boundary_handoff_has_zoom_boundary"] = bool(boundary.get("zoom_boundary_id"))
        checks["boundary_handoff_error_recorded"] = (
            boundary.get("boundary_transfer_error_estimate") is not None
            or boundary_checks.get("boundary_transfer_error_estimate_recorded") is True
        )
    if boundary_has_pyramid or transition:
        checks["transition_handoff_present_for_pyramid_boundary"] = bool(transition)
        checks["transition_handoff_kind_matches_boundary"] = (
            not boundary_has_pyramid or "pyramid" in transition_kinds
        )
        checks["transition_handoff_connects_hex_tet"] = (
            not boundary_has_pyramid
            or transition_checks.get("transition_connects_required_roles") is True
        )
    if boundary_surface_kinds or transition_surface_kinds:
        checks["transition_handoff_surface_kinds_recorded"] = bool(transition_surface_kinds)
        checks["transition_handoff_surface_kinds_match_boundary"] = (
            not boundary_surface_kinds
            or boundary_surface_kinds.issubset(transition_surface_kinds)
        )
    if boundary_material_names or transition_material_names:
        checks["boundary_material_names_recorded"] = bool(boundary_material_names)
        checks["transition_handoff_material_names_recorded"] = bool(transition_material_names)
        checks["transition_handoff_material_names_match_boundary"] = (
            not boundary_material_names
            or boundary_material_names.issubset(transition_material_names)
        )
    if boundary_allowed_zero_material_names or transition_allowed_zero_material_names:
        checks["transition_handoff_zero_material_names_match_boundary"] = (
            not boundary_allowed_zero_material_names
            or boundary_allowed_zero_material_names.issubset(transition_allowed_zero_material_names)
        )
    if boundary_interface_roles or transition_interface_roles:
        checks["boundary_interface_roles_recorded"] = bool(boundary_interface_roles)
        checks["transition_handoff_interface_roles_recorded"] = bool(transition_interface_roles)
        checks["transition_handoff_interface_roles_match_boundary"] = (
            not boundary_interface_roles
            or boundary_interface_roles.issubset(transition_interface_roles)
        )

    issues = [name for name, ok in checks.items() if not ok]
    return {
        "policy": "build123d_submodel_cad_handoff_gate",
        "status": "ok" if not issues else "needs_attention",
        "shape_names": names,
        "geometry_ids": sorted(set(geometry_ids)),
        "expected_geometry_ids": None if expected_ids is None else sorted(expected_ids),
        "recipe_id": recipe_text,
        "parent_model_id": parent_text,
        "submodel_region_id": submodel_text,
        "crop_box": crop,
        "export_id": export_text,
        "unit": unit_text,
        "file_kinds": file_kinds,
        "file_paths": file_paths,
        "boundary_handoff_policy": boundary.get("policy"),
        "boundary_handoff_status": boundary.get("status"),
        "boundary_volume_kind_counts": boundary_volume_counts,
        "boundary_surface_kind_counts": boundary_surface_counts,
        "boundary_material_names": sorted(boundary_material_names),
        "boundary_allowed_zero_material_names": sorted(boundary_allowed_zero_material_names),
        "boundary_interface_roles": sorted(boundary_interface_roles),
        "transition_handoff_policy": transition.get("policy"),
        "transition_handoff_status": transition.get("status"),
        "transition_handoff_kinds": sorted(transition_kinds),
        "transition_handoff_surface_kinds": sorted(transition_surface_kinds),
        "transition_handoff_material_names": sorted(transition_material_names),
        "transition_handoff_allowed_zero_material_names": sorted(transition_allowed_zero_material_names),
        "transition_handoff_interface_roles": sorted(transition_interface_roles),
        "checks": checks,
        "issues": issues,
        "notes": [
            "Use this before a local build123d CAD crop is promoted to Cubit/Netgen meshing.",
            "The CAD recipe, crop box, shape measurements, files, and boundary-handoff gate must travel together.",
            "A good local volume is not enough if the parent model or inherited boundary contract is missing.",
            "When the downstream Cubit boundary handoff contains a pyramid bridge, also attach the build123d hex-to-tet transition intent gate.",
            "If the downstream Cubit gate has surface-family evidence, keep the expected quad/triangle surface-family intent in the build123d transition handoff too.",
            "If the downstream Cubit gate has material/block labels or sidecar row names, keep matching downstream material names in the build123d transition handoff.",
            "If the downstream Cubit gate has interface-adjacency roles, keep matching hex-to-transition and transition-to-tet intent in the build123d transition handoff.",
        ],
    }


def shape_curvilinear_mesh_intent_gate(
    shape_rows,
    *,
    downstream_manifest_gate=None,
    downstream_order_series_gate=None,
    required_roles=("hex_region",),
    expected_route="cubit_hex_or_mixed_path",
    expected_handoff="cubit_curvilinear_handoff",
) -> dict:
    """Check build123d CAD rows before Cubit curvilinear mesh handoff.

    build123d owns CAD shape intent and mass properties.  Cubit owns high-order
    hex/mixed meshing.  This gate keeps that split explicit by requiring CAD
    rows to carry volume/area/bbox plus a mesh role, route, and downstream
    handoff label before the Cubit curvilinear manifest is allowed to consume
    the STEP package.
    """

    rows = [dict(row) for row in shape_rows]
    if not rows:
        raise ValueError("shape_rows must not be empty")
    manifest_gate = dict(downstream_manifest_gate or {})
    order_series_gate = dict(downstream_order_series_gate or {})
    required = {str(role).strip() for role in required_roles if str(role).strip()}
    route = str(expected_route or "").strip()
    handoff = str(expected_handoff or "").strip()

    def has_bbox(row):
        bbox = row.get("bounding_box") or row.get("bbox")
        return isinstance(bbox, dict) and bbox.get("size") is not None

    names = [str(row.get("name", "")).strip() for row in rows]
    geometry_ids = [str(row.get("geometry_id", "")).strip() for row in rows]
    roles = {str(row.get("role", row.get("mesh_role", ""))).strip() for row in rows if str(row.get("role", row.get("mesh_role", ""))).strip()}
    routes = {str(row.get("mesh_route", row.get("routing_hint", ""))).strip() for row in rows if str(row.get("mesh_route", row.get("routing_hint", ""))).strip()}
    handoffs = {str(row.get("downstream_handoff", row.get("handoff", ""))).strip() for row in rows if str(row.get("downstream_handoff", row.get("handoff", ""))).strip()}
    manifest_checks = manifest_gate.get("checks", {})
    if not isinstance(manifest_checks, dict):
        manifest_checks = {}
    order_series_checks = order_series_gate.get("checks", {})
    if not isinstance(order_series_checks, dict):
        order_series_checks = {}
    checks = {
        "shape_rows_present": bool(rows),
        "shape_names_recorded": all(bool(name) for name in names),
        "shape_geometry_ids_recorded": all(bool(value) for value in geometry_ids),
        "shape_geometry_ids_unique": len(set(geometry_ids)) == len(geometry_ids),
        "shape_rows_have_volume_area_bbox": all(
            row.get("volume") is not None and row.get("area") is not None and has_bbox(row)
            for row in rows
        ),
        "required_roles_present": required.issubset(roles),
        "route_is_cubit_hex_or_mixed": routes == {route},
        "handoff_label_recorded": handoffs == {handoff},
        "not_tet_only_route": "netgen_tri_tet_path" not in routes and "tet_only" not in roles,
        "downstream_manifest_ok": not manifest_gate or manifest_gate.get("status") == "ok",
        "downstream_projection_error_ok": not manifest_gate or (
            manifest_checks.get("projection_error_recorded") is True
            and manifest_checks.get("projection_error_within_tolerance") is True
        ),
        "downstream_negative_jacobian_zero": not manifest_gate or (
            manifest_checks.get("negative_jacobian_count_recorded") is True
            and manifest_checks.get("negative_jacobian_count_zero") is True
        ),
        "downstream_order_series_ok": not order_series_gate or order_series_gate.get("status") == "ok",
        "downstream_order_series_policy_known": not order_series_gate or (
            order_series_gate.get("policy") == "cubit_mixed_order_series_inventory_gate"
        ),
        "downstream_order_series_topology_invariant": not order_series_gate or (
            order_series_checks.get("volume_kind_counts_invariant") is True
            and order_series_checks.get("surface_kind_counts_invariant") is True
        ),
        "downstream_order_series_route_matches": not order_series_gate or (
            order_series_checks.get("routing_hint_is_cubit_mixed") is True
        ),
        "downstream_order_series_first_order_inventory": not order_series_gate or (
            order_series_checks.get("first_order_inventory_present") is True
            and order_series_checks.get("first_order_inventory_not_curved") is True
        ),
    }
    issues = []
    if not checks["shape_rows_have_volume_area_bbox"]:
        issues.append("CAD rows must carry volume, area, and bounding_box before curvilinear handoff")
    if not checks["required_roles_present"]:
        issues.append("CAD rows do not include all required mesh roles")
    if not checks["route_is_cubit_hex_or_mixed"]:
        issues.append("curvilinear mesh intent must route to Cubit hex/mixed, not tet-only Netgen")
    if not checks["downstream_manifest_ok"]:
        issues.append("downstream Cubit curvilinear manifest gate is not ok")
    if not checks["downstream_projection_error_ok"]:
        issues.append("downstream Cubit curvilinear manifest must record projection error within tolerance")
    if not checks["downstream_negative_jacobian_zero"]:
        issues.append("downstream Cubit curvilinear manifest must record zero negative-Jacobian elements")
    if not checks["downstream_order_series_ok"]:
        issues.append("downstream Cubit order-series inventory gate is not ok")
    if not checks["downstream_order_series_policy_known"]:
        issues.append("downstream Cubit order-series gate policy is not cubit_mixed_order_series_inventory_gate")
    if not checks["downstream_order_series_topology_invariant"]:
        issues.append("downstream Cubit order series must keep volume and surface topology invariant")
    if not checks["downstream_order_series_route_matches"]:
        issues.append("downstream Cubit order series must keep cubit_hex_or_mixed_path routing")
    if not checks["downstream_order_series_first_order_inventory"]:
        issues.append("downstream Cubit order series must include non-curved first-order inventory")

    return {
        "policy": "build123d_curvilinear_mesh_intent_gate",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "shape_names": names,
        "geometry_ids": sorted(set(geometry_ids)),
        "roles": sorted(roles),
        "routes": sorted(routes),
        "handoffs": sorted(handoffs),
        "required_roles": sorted(required),
        "expected_route": route,
        "expected_handoff": handoff,
        "downstream_manifest_policy": manifest_gate.get("policy"),
        "downstream_manifest_checks": manifest_checks,
        "downstream_order_series_policy": order_series_gate.get("policy"),
        "downstream_order_series_checks": order_series_checks,
        "checks": checks,
        "issues": issues,
        "notes": [
            "Use this on build123d CAD rows before Cubit owns high-order hex or mixed meshing.",
            "The first-order tri/tet education path remains separate from Cubit curvilinear handoff.",
            "Volume/area/bbox are CAD evidence; order, projection error, and Jacobian quality remain downstream Cubit evidence.",
            "If a Cubit order series is attached, first-order topology/routing must remain invariant when curvedelements grow with order.",
        ],
    }


def shape_mesh_environment_handoff_gate(
    shape_rows,
    mesh_environment_gate,
    *,
    expected_route="cubit_hex_or_mixed_path",
    expected_environment_policy="cubit_headless_installation_route_gate",
) -> dict:
    """Check that CAD intent is handed to a verified mesh environment.

    build123d owns geometry and mass properties, but it should not imply that a
    downstream mesher version was actually available.  This gate binds CAD rows
    to a mesh-environment replay gate so release-note watchlists cannot be
    mistaken for installed headless execution evidence.
    """

    rows = [dict(row) for row in shape_rows]
    if not rows:
        raise ValueError("shape_rows must not be empty")
    environment = dict(mesh_environment_gate or {})
    route = str(expected_route or "").strip()
    expected_policy = str(expected_environment_policy or "").strip()

    def has_bbox(row):
        bbox = row.get("bounding_box") or row.get("bbox")
        return isinstance(bbox, dict) and bbox.get("size") is not None

    names = [str(row.get("name", "")).strip() for row in rows]
    geometry_ids = [str(row.get("geometry_id", "")).strip() for row in rows]
    routes = {
        str(row.get("mesh_route", row.get("routing_hint", ""))).strip()
        for row in rows
        if str(row.get("mesh_route", row.get("routing_hint", ""))).strip()
    }
    env_checks = environment.get("checks", {})
    if not isinstance(env_checks, dict):
        env_checks = {}
    license_status = str(environment.get("license_status", "")).strip()
    version_probe_command = str(environment.get("version_probe_command", "")).strip()
    version_probe_summary = environment.get("version_probe_summary", {})
    if version_probe_summary is None:
        version_probe_summary = {}
    if not isinstance(version_probe_summary, dict):
        version_probe_summary = {}
    checks = {
        "shape_rows_present": bool(rows),
        "shape_names_recorded": all(bool(name) for name in names),
        "shape_geometry_ids_recorded": all(bool(value) for value in geometry_ids),
        "shape_geometry_ids_unique": len(set(geometry_ids)) == len(geometry_ids),
        "shape_rows_have_volume_area_bbox": all(
            row.get("volume") is not None and row.get("area") is not None and has_bbox(row)
            for row in rows
        ),
        "mesh_route_matches_expected": routes == {route},
        "mesh_environment_gate_ok": environment.get("status") == "ok",
        "mesh_environment_policy_known": environment.get("policy") == expected_policy,
        "installed_version_recorded": bool(environment.get("installed_version")),
        "binary_path_recorded": bool(environment.get("binary_path")),
        "binary_exists": environment.get("binary_exists") is True
        or env_checks.get("binary_exists") is True,
        "binary_path_is_console_com": env_checks.get("binary_path_is_console_com") is True,
        "headless_flags_present": env_checks.get("required_headless_flags_present") is True,
        "live_claim_matches_installed": env_checks.get("live_claim_matches_installed_version") is True,
        "release_note_watchlist_not_live_claim": env_checks.get("release_note_watchlist_not_live_claim") is True,
        "license_status_allows_headless_probe": (
            not license_status or env_checks.get("license_status_allows_headless_probe") is True
        ),
        "version_probe_is_synchronous_console": (
            not version_probe_command or env_checks.get("version_probe_is_synchronous_console") is True
        ),
        "version_probe_uses_recorded_binary": (
            not version_probe_command or env_checks.get("version_probe_uses_recorded_binary") is True
        ),
        "version_probe_summary_records_installed_version": (
            not version_probe_summary
            or env_checks.get("version_probe_summary_records_installed_version") is True
        ),
        "version_probe_summary_records_license_status": (
            not version_probe_summary
            or env_checks.get("version_probe_summary_records_license_status") is True
        ),
    }
    issues = []
    if not checks["shape_rows_have_volume_area_bbox"]:
        issues.append("build123d CAD rows must carry volume, area, and bbox before mesh-environment handoff")
    if not checks["mesh_route_matches_expected"]:
        issues.append("CAD mesh_route does not match the expected downstream mesh lane")
    if not checks["mesh_environment_gate_ok"]:
        issues.append("mesh environment gate is not ok")
    if not checks["binary_path_is_console_com"]:
        issues.append("mesh environment binary must be the console coreform_cubit.com executable")
    if not checks["version_probe_uses_recorded_binary"]:
        issues.append("mesh environment version probe must use the recorded binary path")
    if not checks["live_claim_matches_installed"]:
        issues.append("mesh environment live claim does not match installed version")
    if not checks["version_probe_summary_records_installed_version"]:
        issues.append("mesh environment version-probe summary does not record the installed version")
    return {
        "policy": "build123d_mesh_environment_handoff_gate",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "shape_names": names,
        "geometry_ids": sorted(set(geometry_ids)),
        "routes": sorted(routes),
        "expected_route": route,
        "mesh_environment_policy": environment.get("policy"),
        "installed_version": environment.get("installed_version"),
        "binary_path": environment.get("binary_path"),
        "license_status": license_status,
        "version_probe_command": version_probe_command,
        "version_probe_summary": version_probe_summary,
        "release_note_version": environment.get("release_note_version"),
        "live_claimed_release_version": environment.get("live_claimed_release_version"),
        "checks": checks,
        "issues": issues,
        "notes": [
            "Use this before sending build123d CAD rows into Cubit/Coreform headless meshing.",
            "A release-note watchlist is useful learning, but installed-version evidence must support live mesh claims.",
            "Slot354-style Cubit handoff requires coreform_cubit.com, not the GUI launcher executable, and the version probe must use that same recorded binary.",
            "Volume/area/bbox remain build123d evidence; headless flags and mesher version remain downstream environment evidence.",
            "When provided, sanitized license status and synchronous version-probe summary remain downstream environment evidence, not CAD evidence.",
        ],
    }


def shape_measurement_comparison_summary(
    reference_rows,
    measured_rows,
    rtol=1.0e-5,
    measured_label="measured",
    bbox_atol=1.0e-6,
):
    """Return compact summary statistics for a measurement comparison."""

    rows = compare_shape_measurement_rows(
        reference_rows,
        measured_rows,
        rtol=rtol,
        measured_label=measured_label,
        bbox_atol=bbox_atol,
    )
    volume_errors = [row["volume_rel_error"] or 0.0 for row in rows]
    area_errors = [row["area_rel_error"] or 0.0 for row in rows]
    bbox_errors = [row["bbox_abs_error"] or 0.0 for row in rows]
    return {
        "measured_label": measured_label,
        "rtol": float(rtol),
        "bbox_atol": float(bbox_atol),
        "n_cases": len(rows),
        "n_passed": sum(1 for row in rows if row["passed"]),
        "n_bbox_compared": sum(1 for row in rows if row["bbox_compared"]),
        "max_volume_rel_error": max(volume_errors) if volume_errors else 0.0,
        "max_area_rel_error": max(area_errors) if area_errors else 0.0,
        "max_bbox_abs_error": max(bbox_errors) if bbox_errors else 0.0,
        "rows": rows,
    }


def _shape_measurement_comparison_score(row):
    if row.get("volume_rel_error") is None or row.get("area_rel_error") is None:
        return math.inf
    score = max(float(row["volume_rel_error"]), float(row["area_rel_error"]))
    bbox_error = row.get("bbox_abs_error")
    if bbox_error is not None:
        bbox_atol = float(row.get("bbox_atol") or 0.0)
        score = max(score, float(bbox_error) / bbox_atol if bbox_atol > 0.0 else float(bbox_error))
    return score


def worst_shape_measurement_comparison_rows(comparison_rows, limit=5):
    """Return the largest measurement mismatches first."""

    if limit < 0:
        raise ValueError("limit must be non-negative")
    rows = sorted(
        list(comparison_rows),
        key=lambda row: (
            0 if bool(row.get("passed")) else 1,
            _shape_measurement_comparison_score(row),
            str(row.get("name", "")),
        ),
        reverse=True,
    )
    return rows[:limit]


def shape_measurement_health_summary(
    reference_rows,
    measured_rows,
    rtol=1.0e-5,
    measured_label="measured",
    bbox_atol=1.0e-6,
    worst_limit=5,
):
    """Return a readable health report for CAD measurement cross validation."""

    if worst_limit < 0:
        raise ValueError("worst_limit must be non-negative")
    reference_rows = list(reference_rows)
    inventory = shape_measurement_inventory_summary(reference_rows)
    comparison = shape_measurement_comparison_summary(
        reference_rows,
        measured_rows,
        rtol=rtol,
        measured_label=measured_label,
        bbox_atol=bbox_atol,
    )
    checks = {
        "all_reference_shapes_valid": inventory["n_valid"] == inventory["n_shapes"],
        "all_measurements_present_and_within_tolerance": comparison["n_passed"] == comparison["n_cases"],
    }
    issues = []
    if not checks["all_reference_shapes_valid"]:
        issues.append("at least one reference build123d shape is invalid")
    if not checks["all_measurements_present_and_within_tolerance"]:
        issues.append("at least one external measurement is missing or outside tolerance")
    ok = all(checks.values())
    return {
        "policy": "build123d_shape_measurement_volume_area_bbox_health",
        "status": "ok" if ok else "needs_attention",
        "ok_for_geometry_roundtrip": ok,
        "measured_label": measured_label,
        "checks": checks,
        "issues": issues,
        "inventory": inventory,
        "comparison_summary": {
            key: value
            for key, value in comparison.items()
            if key != "rows"
        },
        "worst_comparisons": worst_shape_measurement_comparison_rows(
            comparison["rows"],
            limit=worst_limit,
        ),
    }


def box_through_cylinder_reference_row(x, y, z, radius, axis="z", label="box_hole"):
    """Analytic mass-property row for a centered box with a through cylindrical hole.

    This is a compact build123d slot gate: create the CAD with
    ``Box(x, y, z) - Cylinder(radius, height > hole_length)`` and compare the
    resulting OCCT mass properties with this closed form before meshing.  The
    helper returns the same row shape expected by
    :func:`compare_shape_measurement_rows`.
    """

    x = float(x)
    y = float(y)
    z = float(z)
    radius = float(radius)
    axis_key = str(axis).lower().strip()
    dims = {"x": x, "y": y, "z": z}
    if axis_key not in dims:
        raise ValueError("axis must be 'x', 'y', or 'z'")
    if x <= 0.0 or y <= 0.0 or z <= 0.0:
        raise ValueError("box dimensions must be positive")
    if radius <= 0.0:
        raise ValueError("hole radius must be positive")
    perpendicular = [value for key, value in dims.items() if key != axis_key]
    if 2.0 * radius >= min(perpendicular):
        raise ValueError("hole diameter must be smaller than both perpendicular box dimensions")

    hole_length = dims[axis_key]
    box_volume = x * y * z
    box_area = 2.0 * (x * y + y * z + x * z)
    circle_area = math.pi * radius * radius
    lateral_area = 2.0 * math.pi * radius * hole_length
    volume = box_volume - circle_area * hole_length
    area = box_area - 2.0 * circle_area + lateral_area
    bbox = {
        "min": [-x / 2.0, -y / 2.0, -z / 2.0],
        "max": [x / 2.0, y / 2.0, z / 2.0],
        "center": [0.0, 0.0, 0.0],
        "size": [x, y, z],
        "diagonal": math.sqrt(x * x + y * y + z * z),
    }
    return {
        "name": str(label),
        "volume": volume,
        "area": area,
        "bounding_box": bbox,
        "axis": axis_key,
        "hole_radius": radius,
        "hole_length": hole_length,
        "policy": "analytic_box_through_cylinder_mass_property_reference",
    }


def mounting_plate_boss_reference_row(
    base_x,
    base_y,
    base_h,
    boss_r,
    boss_h,
    central_hole_r,
    corner_hole_r,
    corner_hole_x,
    corner_hole_y,
    label="mounting_plate_boss_five_holes",
):
    """Analytic volume row for a plate + cylindrical boss + five vertical holes.

    The geometry is a centered rectangular plate with a cylindrical boss on
    its top face, one central through-hole through plate and boss, and four
    corner holes through the plate only.  It is a compact CAD-kernel
    round-trip gate for build123d -> STEP -> Cubit/CST volume checks.
    """

    base_x = float(base_x)
    base_y = float(base_y)
    base_h = float(base_h)
    boss_r = float(boss_r)
    boss_h = float(boss_h)
    central_hole_r = float(central_hole_r)
    corner_hole_r = float(corner_hole_r)
    corner_hole_x = abs(float(corner_hole_x))
    corner_hole_y = abs(float(corner_hole_y))
    dims = (base_x, base_y, base_h, boss_r, boss_h, central_hole_r, corner_hole_r)
    if any(value <= 0.0 for value in dims):
        raise ValueError("all dimensions and radii must be positive")
    if 2.0 * boss_r >= min(base_x, base_y):
        raise ValueError("boss diameter must fit on the plate")
    if central_hole_r >= boss_r:
        raise ValueError("central hole radius must be smaller than boss radius")
    if corner_hole_x + corner_hole_r >= base_x / 2.0:
        raise ValueError("corner hole x location must fit inside the plate")
    if corner_hole_y + corner_hole_r >= base_y / 2.0:
        raise ValueError("corner hole y location must fit inside the plate")
    if math.hypot(corner_hole_x, corner_hole_y) - corner_hole_r <= boss_r:
        raise ValueError("corner holes must not intersect the boss footprint")

    base_volume = base_x * base_y * base_h
    boss_volume = math.pi * boss_r * boss_r * boss_h
    central_hole_volume = math.pi * central_hole_r * central_hole_r * (base_h + boss_h)
    corner_hole_volume = 4.0 * math.pi * corner_hole_r * corner_hole_r * base_h
    volume = base_volume + boss_volume - central_hole_volume - corner_hole_volume
    bbox = {
        "min": [-base_x / 2.0, -base_y / 2.0, -base_h / 2.0],
        "max": [base_x / 2.0, base_y / 2.0, base_h / 2.0 + boss_h],
        "center": [0.0, 0.0, boss_h / 2.0],
        "size": [base_x, base_y, base_h + boss_h],
        "diagonal": math.sqrt(base_x * base_x + base_y * base_y + (base_h + boss_h) ** 2),
    }
    return {
        "name": str(label),
        "volume": volume,
        "bounding_box": bbox,
        "terms": {
            "base": base_volume,
            "boss": boss_volume,
            "central_hole": -central_hole_volume,
            "four_corner_holes": -corner_hole_volume,
        },
        "policy": "analytic_mounting_plate_boss_volume_reference",
    }


def keyed_terminal_plate_reference_row(
    base_x,
    base_y,
    base_h,
    boss_r,
    boss_h,
    boss_x,
    boss_hole_r,
    window_x,
    window_y,
    key_slot_x,
    key_slot_y,
    mount_hole_r,
    mount_hole_x,
    mount_hole_y,
    label="keyed_terminal_plate_two_bosses",
):
    """Analytic volume row for a keyed terminal plate with two bosses.

    The body is a centered rectangular base plate, two cylindrical bosses on
    the top face, a centered rectangular window through the base, a centered
    edge key-slot cut through the base at ``+Y``, two boss holes through
    plate+boss, and two mirrored mounting holes through the base only.  This
    is a compact motor-fixture/terminal-plate CAD gate for build123d -> STEP
    -> external CAD-kernel volume checks.
    """

    base_x = float(base_x)
    base_y = float(base_y)
    base_h = float(base_h)
    boss_r = float(boss_r)
    boss_h = float(boss_h)
    boss_x = abs(float(boss_x))
    boss_hole_r = float(boss_hole_r)
    window_x = float(window_x)
    window_y = float(window_y)
    key_slot_x = float(key_slot_x)
    key_slot_y = float(key_slot_y)
    mount_hole_r = float(mount_hole_r)
    mount_hole_x = abs(float(mount_hole_x))
    mount_hole_y = float(mount_hole_y)

    positive = (
        base_x, base_y, base_h, boss_r, boss_h, boss_x, boss_hole_r,
        window_x, window_y, key_slot_x, key_slot_y, mount_hole_r, mount_hole_x,
    )
    if any(value <= 0.0 for value in positive):
        raise ValueError("all dimensions and radii except mount_hole_y must be positive")
    if boss_hole_r >= boss_r:
        raise ValueError("boss hole radius must be smaller than boss radius")
    if boss_x + boss_r >= base_x / 2.0:
        raise ValueError("boss footprint must fit on the plate")
    if 2.0 * boss_r >= base_y:
        raise ValueError("boss diameter must fit inside the plate width")
    if window_x >= base_x or window_y >= base_y:
        raise ValueError("window must fit inside the plate")
    if key_slot_x >= base_x or key_slot_y >= base_y:
        raise ValueError("key slot must fit inside the plate")
    if base_y / 2.0 - key_slot_y - window_y / 2.0 <= 0.0:
        raise ValueError("center window and edge key slot must not overlap")
    if boss_x - boss_r <= max(window_x, key_slot_x) / 2.0:
        raise ValueError("boss footprint must not overlap the centered window or key slot")
    if mount_hole_x + mount_hole_r >= base_x / 2.0:
        raise ValueError("mount hole x location must fit inside the plate")
    if abs(mount_hole_y) + mount_hole_r >= base_y / 2.0:
        raise ValueError("mount hole y location must fit inside the plate")
    if math.hypot(mount_hole_x - boss_x, mount_hole_y) <= boss_r + mount_hole_r:
        raise ValueError("mount holes must not intersect the boss footprints")

    base_volume = base_x * base_y * base_h
    boss_volume = 2.0 * math.pi * boss_r * boss_r * boss_h
    window_volume = window_x * window_y * base_h
    key_slot_volume = key_slot_x * key_slot_y * base_h
    boss_hole_volume = 2.0 * math.pi * boss_hole_r * boss_hole_r * (base_h + boss_h)
    mount_hole_volume = 2.0 * math.pi * mount_hole_r * mount_hole_r * base_h
    volume = base_volume + boss_volume - window_volume - key_slot_volume - boss_hole_volume - mount_hole_volume
    bbox = {
        "min": [-base_x / 2.0, -base_y / 2.0, -base_h / 2.0],
        "max": [base_x / 2.0, base_y / 2.0, base_h / 2.0 + boss_h],
        "center": [0.0, 0.0, boss_h / 2.0],
        "size": [base_x, base_y, base_h + boss_h],
        "diagonal": math.sqrt(base_x * base_x + base_y * base_y + (base_h + boss_h) ** 2),
    }
    return {
        "name": str(label),
        "volume": volume,
        "bounding_box": bbox,
        "terms": {
            "base": base_volume,
            "two_bosses": boss_volume,
            "rectangular_window": -window_volume,
            "edge_key_slot": -key_slot_volume,
            "two_boss_holes": -boss_hole_volume,
            "two_mount_holes": -mount_hole_volume,
        },
        "policy": "analytic_keyed_terminal_plate_volume_reference",
    }


def flanged_sleeve_reference_row(
    flange_r,
    flange_h,
    hub_r,
    hub_h,
    bore_r,
    bolt_circle_r,
    bolt_r,
    bolt_count=4,
    label="flanged_sleeve_four_bolt_holes",
):
    """Analytic volume row for a coaxial flanged sleeve with bolt holes.

    The body is a lower flange annulus plus a raised hub annulus, with a
    central bore through both regions and equally spaced vertical bolt holes
    through the flange only.  It is a compact motor-fixture / bearing-seat CAD
    gate for build123d -> STEP -> external CAD-kernel volume checks.
    """

    flange_r = float(flange_r)
    flange_h = float(flange_h)
    hub_r = float(hub_r)
    hub_h = float(hub_h)
    bore_r = float(bore_r)
    bolt_circle_r = float(bolt_circle_r)
    bolt_r = float(bolt_r)
    bolt_count = int(bolt_count)

    if any(value <= 0.0 for value in (flange_r, flange_h, hub_r, hub_h, bore_r, bolt_circle_r, bolt_r)):
        raise ValueError("all dimensions and radii must be positive")
    if bolt_count < 1:
        raise ValueError("bolt_count must be positive")
    if not (bore_r < hub_r < flange_r):
        raise ValueError("require bore_r < hub_r < flange_r")
    if bolt_circle_r + bolt_r >= flange_r:
        raise ValueError("bolt holes must fit inside the flange")
    if bolt_circle_r - bolt_r <= hub_r:
        raise ValueError("bolt holes must not intersect the raised hub footprint")
    if bolt_circle_r - bolt_r <= bore_r:
        raise ValueError("bolt holes must not intersect the central bore")
    if bolt_count > 1:
        chord = 2.0 * bolt_circle_r * math.sin(math.pi / bolt_count)
        if chord <= 2.0 * bolt_r:
            raise ValueError("adjacent bolt holes must not overlap")

    flange_volume = math.pi * (flange_r * flange_r - bore_r * bore_r) * flange_h
    hub_volume = math.pi * (hub_r * hub_r - bore_r * bore_r) * hub_h
    bolt_hole_volume = bolt_count * math.pi * bolt_r * bolt_r * flange_h
    volume = flange_volume + hub_volume - bolt_hole_volume
    total_h = flange_h + hub_h
    bbox = {
        "min": [-flange_r, -flange_r, -flange_h / 2.0],
        "max": [flange_r, flange_r, flange_h / 2.0 + hub_h],
        "center": [0.0, 0.0, hub_h / 2.0],
        "size": [2.0 * flange_r, 2.0 * flange_r, total_h],
        "diagonal": math.sqrt((2.0 * flange_r) ** 2 + (2.0 * flange_r) ** 2 + total_h ** 2),
    }
    return {
        "name": str(label),
        "volume": volume,
        "bounding_box": bbox,
        "terms": {
            "flange_annulus": flange_volume,
            "hub_annulus": hub_volume,
            "bolt_holes": -bolt_hole_volume,
        },
        "bolt_count": bolt_count,
        "policy": "analytic_flanged_sleeve_volume_reference",
    }


def coax_annular_sleeve_reference_row(
    inner_radius,
    outer_radius,
    height,
    label="coax_annular_sleeve",
):
    """Analytic row for a hollow coaxial sleeve.

    This is the CAD companion of the coaxial C/R field gate: before a motor
    drive, cable shield, or winding-insulation model promotes a capacitance or
    resistance into an equivalent circuit, verify that the annular solid volume
    survived the build123d -> STEP -> external CAD-kernel round trip.
    """

    r_in = float(inner_radius)
    r_out = float(outer_radius)
    h = float(height)
    if not (0.0 < r_in < r_out):
        raise ValueError("require 0 < inner_radius < outer_radius")
    if h <= 0.0:
        raise ValueError("height must be positive")

    volume = math.pi * (r_out * r_out - r_in * r_in) * h
    area = (
        2.0 * math.pi * r_out * h
        + 2.0 * math.pi * r_in * h
        + 2.0 * math.pi * (r_out * r_out - r_in * r_in)
    )
    bbox = {
        "min": [-r_out, -r_out, -h / 2.0],
        "max": [r_out, r_out, h / 2.0],
        "center": [0.0, 0.0, 0.0],
        "size": [2.0 * r_out, 2.0 * r_out, h],
        "diagonal": math.sqrt((2.0 * r_out) ** 2 + (2.0 * r_out) ** 2 + h * h),
    }
    return {
        "name": str(label),
        "volume": volume,
        "area": area,
        "bounding_box": bbox,
        "terms": {
            "outer_cylinder": math.pi * r_out * r_out * h,
            "inner_void": -math.pi * r_in * r_in * h,
        },
        "parameters": {
            "inner_radius": r_in,
            "outer_radius": r_out,
            "height": h,
        },
        "policy": "analytic_coax_annular_sleeve_volume_reference",
    }


def ribbed_busbar_heat_sink_reference_row(
    base_x,
    base_y,
    base_h,
    fin_count,
    fin_w,
    fin_h,
    fin_pitch,
    hole_r,
    hole_x,
    hole_y,
    label="ribbed_busbar_heat_sink_four_holes",
):
    """Analytic volume row for a finned busbar / heat-sink plate.

    The body is a centered rectangular base plate with evenly spaced straight
    ribs on the top face and four vertical bolt holes through the base only.
    Holes are kept outside the fin band on purpose, so the volume decomposes
    into three readable terms: base, ribs, and bolt-hole removal.  This gives a
    compact CAD-kernel round-trip gate for motor terminals, heat-spreader
    plates, and conduction-cooled busbar fixtures before meshing.
    """

    base_x = float(base_x)
    base_y = float(base_y)
    base_h = float(base_h)
    fin_count = int(fin_count)
    fin_w = float(fin_w)
    fin_h = float(fin_h)
    fin_pitch = float(fin_pitch)
    hole_r = float(hole_r)
    hole_x = abs(float(hole_x))
    hole_y = abs(float(hole_y))

    if any(value <= 0.0 for value in (base_x, base_y, base_h, fin_w, fin_h, fin_pitch, hole_r)):
        raise ValueError("all dimensions, pitches, and radii must be positive")
    if fin_count < 1:
        raise ValueError("fin_count must be positive")
    if fin_w >= fin_pitch:
        raise ValueError("fin_w must be smaller than fin_pitch")

    fin_band_half_width = 0.5 * ((fin_count - 1) * fin_pitch + fin_w)
    if fin_band_half_width >= base_y / 2.0:
        raise ValueError("fin band must fit on the base")
    if hole_x + hole_r >= base_x / 2.0:
        raise ValueError("bolt hole x location must fit inside the base")
    if hole_y + hole_r >= base_y / 2.0:
        raise ValueError("bolt hole y location must fit inside the base")
    if hole_y - hole_r <= fin_band_half_width:
        raise ValueError("bolt holes must stay outside the fin footprint")

    base_volume = base_x * base_y * base_h
    fin_volume = fin_count * base_x * fin_w * fin_h
    bolt_hole_volume = 4.0 * math.pi * hole_r * hole_r * base_h
    volume = base_volume + fin_volume - bolt_hole_volume
    total_h = base_h + fin_h
    bbox = {
        "min": [-base_x / 2.0, -base_y / 2.0, -base_h / 2.0],
        "max": [base_x / 2.0, base_y / 2.0, base_h / 2.0 + fin_h],
        "center": [0.0, 0.0, fin_h / 2.0],
        "size": [base_x, base_y, total_h],
        "diagonal": math.sqrt(base_x * base_x + base_y * base_y + total_h * total_h),
    }
    return {
        "name": str(label),
        "volume": volume,
        "bounding_box": bbox,
        "terms": {
            "base": base_volume,
            "straight_ribs": fin_volume,
            "four_base_holes": -bolt_hole_volume,
        },
        "fin_count": fin_count,
        "clearances": {
            "hole_to_x_edge": base_x / 2.0 - hole_x - hole_r,
            "hole_to_y_edge": base_y / 2.0 - hole_y - hole_r,
            "hole_to_fin_band": hole_y - hole_r - fin_band_half_width,
            "fin_band_to_y_edge": base_y / 2.0 - fin_band_half_width,
        },
        "policy": "analytic_ribbed_busbar_heat_sink_volume_reference",
    }


def three_phase_busbar_snubber_plate_reference_row(
    base_x,
    base_y,
    base_h,
    phase_count,
    phase_tab_x,
    phase_tab_y,
    phase_tab_h,
    phase_pitch,
    snubber_count,
    snubber_pad_x,
    snubber_pad_y,
    snubber_pad_h,
    snubber_pitch,
    mount_hole_r,
    mount_hole_x,
    mount_hole_y,
    phase_tab_y0=0.0,
    snubber_pad_y0=-1.0,
    label="three_phase_busbar_snubber_plate",
):
    """Analytic volume row for a motor-drive busbar with snubber pads.

    The body is a rectangular busbar plate, three raised phase terminal tabs,
    two raised snubber/component pads, and four mounting holes through the base
    only.  It is intentionally decomposed into readable terms for public
    motor-drive CAD validation before STEP handoff to Cubit/CST or downstream
    meshing.
    """

    base_x = float(base_x)
    base_y = float(base_y)
    base_h = float(base_h)
    phase_count = int(phase_count)
    phase_tab_x = float(phase_tab_x)
    phase_tab_y = float(phase_tab_y)
    phase_tab_h = float(phase_tab_h)
    phase_pitch = float(phase_pitch)
    snubber_count = int(snubber_count)
    snubber_pad_x = float(snubber_pad_x)
    snubber_pad_y = float(snubber_pad_y)
    snubber_pad_h = float(snubber_pad_h)
    snubber_pitch = float(snubber_pitch)
    mount_hole_r = float(mount_hole_r)
    mount_hole_x = abs(float(mount_hole_x))
    mount_hole_y = abs(float(mount_hole_y))
    phase_tab_y0 = float(phase_tab_y0)
    snubber_pad_y0 = float(snubber_pad_y0)

    positives = (
        base_x, base_y, base_h, phase_tab_x, phase_tab_y, phase_tab_h,
        phase_pitch, snubber_pad_x, snubber_pad_y, snubber_pad_h,
        snubber_pitch, mount_hole_r, mount_hole_x, mount_hole_y,
    )
    if any(value <= 0.0 for value in positives):
        raise ValueError("all dimensions, pitches, and radii must be positive")
    if phase_count < 1 or snubber_count < 1:
        raise ValueError("phase_count and snubber_count must be positive")

    phase_span = (phase_count - 1) * phase_pitch + phase_tab_x
    snubber_span = (snubber_count - 1) * snubber_pitch + snubber_pad_x
    if phase_span >= base_x:
        raise ValueError("phase terminal tabs must fit across the base")
    if snubber_span >= base_x:
        raise ValueError("snubber pads must fit across the base")
    if abs(phase_tab_y0) + phase_tab_y / 2.0 >= base_y / 2.0:
        raise ValueError("phase terminal tabs must fit inside the base width")
    if abs(snubber_pad_y0) + snubber_pad_y / 2.0 >= base_y / 2.0:
        raise ValueError("snubber pads must fit inside the base width")
    if mount_hole_x + mount_hole_r >= base_x / 2.0:
        raise ValueError("mount hole x location must fit inside the base")
    if mount_hole_y + mount_hole_r >= base_y / 2.0:
        raise ValueError("mount hole y location must fit inside the base")

    phase_centers = [
        ((index - (phase_count - 1) / 2.0) * phase_pitch, phase_tab_y0)
        for index in range(phase_count)
    ]
    snubber_centers = [
        ((index - (snubber_count - 1) / 2.0) * snubber_pitch, snubber_pad_y0)
        for index in range(snubber_count)
    ]
    hole_centers = [
        (sx * mount_hole_x, sy * mount_hole_y)
        for sx in (-1.0, 1.0)
        for sy in (-1.0, 1.0)
    ]

    def _hole_clears_rectangle(hole_x, hole_y, rect_x, rect_y, rect_w, rect_h):
        return (
            abs(hole_x - rect_x) > rect_w / 2.0 + mount_hole_r
            or abs(hole_y - rect_y) > rect_h / 2.0 + mount_hole_r
        )

    for hole_x, hole_y in hole_centers:
        for rect_x, rect_y in phase_centers:
            if not _hole_clears_rectangle(hole_x, hole_y, rect_x, rect_y, phase_tab_x, phase_tab_y):
                raise ValueError("mount holes must not intersect phase terminal tabs")
        for rect_x, rect_y in snubber_centers:
            if not _hole_clears_rectangle(hole_x, hole_y, rect_x, rect_y, snubber_pad_x, snubber_pad_y):
                raise ValueError("mount holes must not intersect snubber pads")

    base_volume = base_x * base_y * base_h
    phase_tab_volume = phase_count * phase_tab_x * phase_tab_y * phase_tab_h
    snubber_pad_volume = snubber_count * snubber_pad_x * snubber_pad_y * snubber_pad_h
    mount_hole_volume = 4.0 * math.pi * mount_hole_r * mount_hole_r * base_h
    top_h = max(phase_tab_h, snubber_pad_h)
    volume = base_volume + phase_tab_volume + snubber_pad_volume - mount_hole_volume
    bbox = {
        "min": [-base_x / 2.0, -base_y / 2.0, -base_h / 2.0],
        "max": [base_x / 2.0, base_y / 2.0, base_h / 2.0 + top_h],
        "center": [0.0, 0.0, top_h / 2.0],
        "size": [base_x, base_y, base_h + top_h],
        "diagonal": math.sqrt(base_x * base_x + base_y * base_y + (base_h + top_h) ** 2),
    }
    return {
        "name": str(label),
        "volume": volume,
        "bounding_box": bbox,
        "terms": {
            "base": base_volume,
            "three_phase_tabs": phase_tab_volume,
            "two_snubber_pads": snubber_pad_volume,
            "four_mount_holes": -mount_hole_volume,
        },
        "counts": {
            "phase_tabs": phase_count,
            "snubber_pads": snubber_count,
            "mount_holes": 4,
        },
        "clearances": {
            "mount_hole_to_x_edge": base_x / 2.0 - mount_hole_x - mount_hole_r,
            "mount_hole_to_y_edge": base_y / 2.0 - mount_hole_y - mount_hole_r,
            "phase_tab_gap": phase_pitch - phase_tab_x,
            "snubber_pad_gap": snubber_pitch - snubber_pad_x,
            "phase_span_to_x_edge": (base_x - phase_span) / 2.0,
            "snubber_span_to_x_edge": (base_x - snubber_span) / 2.0,
        },
        "policy": "analytic_three_phase_busbar_snubber_plate_volume_reference",
    }


def rcd_snubber_heat_spreader_reference_row(
    base_x,
    base_y,
    base_h,
    rib_count,
    rib_x,
    rib_y,
    rib_h,
    rib_pitch,
    snubber_count,
    snubber_pad_x,
    snubber_pad_y,
    snubber_pad_h,
    snubber_pitch,
    mount_hole_r,
    mount_hole_x,
    mount_hole_y,
    snubber_pad_y0=-1.55,
    label="rcd_snubber_heat_spreader",
):
    """Analytic volume row for a readable RCD-snubber heat-spreader plate.

    This is a public-safe motor-drive hardware gate: a base plate, straight
    cooling ribs, two snubber/component pads, and four mounting holes through
    the base only.  Holes are required to stay outside ribs and pads so the
    volume decomposition remains a transparent pre-mesh CAD check.
    """

    base_x = float(base_x)
    base_y = float(base_y)
    base_h = float(base_h)
    rib_count = int(rib_count)
    rib_x = float(rib_x)
    rib_y = float(rib_y)
    rib_h = float(rib_h)
    rib_pitch = float(rib_pitch)
    snubber_count = int(snubber_count)
    snubber_pad_x = float(snubber_pad_x)
    snubber_pad_y = float(snubber_pad_y)
    snubber_pad_h = float(snubber_pad_h)
    snubber_pitch = float(snubber_pitch)
    mount_hole_r = float(mount_hole_r)
    mount_hole_x = abs(float(mount_hole_x))
    mount_hole_y = abs(float(mount_hole_y))
    snubber_pad_y0 = float(snubber_pad_y0)

    positives = (
        base_x, base_y, base_h, rib_x, rib_y, rib_h, rib_pitch,
        snubber_pad_x, snubber_pad_y, snubber_pad_h,
        snubber_pitch, mount_hole_r, mount_hole_x, mount_hole_y,
    )
    if any(value <= 0.0 for value in positives):
        raise ValueError("all dimensions, pitches, and radii must be positive")
    if rib_count < 1 or snubber_count < 1:
        raise ValueError("rib_count and snubber_count must be positive")
    if rib_x >= base_x:
        raise ValueError("ribs must fit inside the base length")
    rib_span = (rib_count - 1) * rib_pitch + rib_y
    if rib_span >= base_y:
        raise ValueError("rib band must fit inside the base width")
    snubber_span = (snubber_count - 1) * snubber_pitch + snubber_pad_x
    if snubber_span >= base_x:
        raise ValueError("snubber pads must fit across the base")
    if snubber_pad_y >= base_y:
        raise ValueError("snubber pads must fit inside the base width")
    if abs(snubber_pad_y0) + snubber_pad_y / 2.0 >= base_y / 2.0:
        raise ValueError("snubber pads must fit inside the base width at snubber_pad_y0")
    if mount_hole_x + mount_hole_r >= base_x / 2.0:
        raise ValueError("mount hole x location must fit inside the base")
    if mount_hole_y + mount_hole_r >= base_y / 2.0:
        raise ValueError("mount hole y location must fit inside the base")

    rib_band_half_width = rib_span / 2.0
    if mount_hole_y - mount_hole_r <= rib_band_half_width:
        raise ValueError("mount holes must stay outside the rib band")
    if mount_hole_x - mount_hole_r <= snubber_span / 2.0:
        raise ValueError("mount holes must stay outside the snubber pad span")
    if abs(snubber_pad_y0) - snubber_pad_y / 2.0 <= rib_band_half_width:
        raise ValueError("snubber pads must stay outside the rib band for additive volume")

    base_volume = base_x * base_y * base_h
    rib_volume = rib_count * rib_x * rib_y * rib_h
    snubber_pad_volume = snubber_count * snubber_pad_x * snubber_pad_y * snubber_pad_h
    mount_hole_volume = 4.0 * math.pi * mount_hole_r * mount_hole_r * base_h
    top_h = max(rib_h, snubber_pad_h)
    volume = base_volume + rib_volume + snubber_pad_volume - mount_hole_volume
    bbox = {
        "min": [-base_x / 2.0, -base_y / 2.0, -base_h / 2.0],
        "max": [base_x / 2.0, base_y / 2.0, base_h / 2.0 + top_h],
        "center": [0.0, 0.0, top_h / 2.0],
        "size": [base_x, base_y, base_h + top_h],
        "diagonal": math.sqrt(base_x * base_x + base_y * base_y + (base_h + top_h) ** 2),
    }
    return {
        "name": str(label),
        "volume": volume,
        "bounding_box": bbox,
        "terms": {
            "base": base_volume,
            "straight_ribs": rib_volume,
            "snubber_pads": snubber_pad_volume,
            "four_mount_holes": -mount_hole_volume,
        },
        "counts": {
            "ribs": rib_count,
            "snubber_pads": snubber_count,
            "mount_holes": 4,
        },
        "clearances": {
            "mount_hole_to_x_edge": base_x / 2.0 - mount_hole_x - mount_hole_r,
            "mount_hole_to_y_edge": base_y / 2.0 - mount_hole_y - mount_hole_r,
            "mount_hole_to_rib_band": mount_hole_y - mount_hole_r - rib_band_half_width,
            "mount_hole_to_snubber_span": mount_hole_x - mount_hole_r - snubber_span / 2.0,
            "snubber_pad_to_rib_band": abs(snubber_pad_y0) - snubber_pad_y / 2.0 - rib_band_half_width,
            "rib_gap": rib_pitch - rib_y,
            "snubber_pad_gap": snubber_pitch - snubber_pad_x,
            "rib_span_to_y_edge": (base_y - rib_span) / 2.0,
            "snubber_span_to_x_edge": (base_x - snubber_span) / 2.0,
            "snubber_pad_to_y_edge": base_y / 2.0 - abs(snubber_pad_y0) - snubber_pad_y / 2.0,
        },
        "parameters": {
            "rib_pitch": rib_pitch,
            "snubber_pad_y0": snubber_pad_y0,
        },
        "policy": "analytic_rcd_snubber_heat_spreader_volume_reference",
    }


def rcd_snubber_capacitance_sweep_rows(
    capacitance_uF_values,
    snubber_pad_x_values,
    *,
    base_x=10.0,
    base_y=4.0,
    base_h=0.32,
    rib_count=5,
    rib_x=8.0,
    rib_y=0.12,
    rib_h=0.45,
    rib_pitch=0.45,
    snubber_count=2,
    snubber_pad_y=0.70,
    snubber_pad_h=0.38,
    snubber_pitch=2.4,
    mount_hole_r=0.16,
    mount_hole_x=4.2,
    mount_hole_y=1.55,
    snubber_pad_y0=-1.55,
    label_prefix="rcd_snubber_heat_spreader",
):
    """Return CAD design-table rows for RCD snubber capacitance variants.

    The electrical capacitance value is not inferred from geometry here.  It is
    carried as explicit provenance next to the pad dimensions and volume terms
    so a motor-drive notebook can compare overshoot rows without losing the CAD
    variant that produced each component footprint.
    """

    capacitances = [float(value) for value in capacitance_uF_values]
    pad_lengths = [float(value) for value in snubber_pad_x_values]
    if not capacitances:
        raise ValueError("capacitance_uF_values must not be empty")
    if len(capacitances) != len(pad_lengths):
        raise ValueError("capacitance_uF_values and snubber_pad_x_values must have the same length")
    if any(value <= 0.0 for value in capacitances):
        raise ValueError("capacitance_uF_values must be positive")
    if any(value <= 0.0 for value in pad_lengths):
        raise ValueError("snubber_pad_x_values must be positive")

    rows = []
    for capacitance_uF, pad_x in zip(capacitances, pad_lengths):
        label = f"{label_prefix}_{capacitance_uF:g}uF"
        row = rcd_snubber_heat_spreader_reference_row(
            base_x, base_y, base_h,
            rib_count, rib_x, rib_y, rib_h, rib_pitch,
            snubber_count, pad_x, snubber_pad_y, snubber_pad_h, snubber_pitch,
            mount_hole_r, mount_hole_x, mount_hole_y,
            snubber_pad_y0=snubber_pad_y0,
            label=label,
        )
        row["capacitance_uF"] = capacitance_uF
        row["snubber_pad_x"] = pad_x
        row["snubber_pad_volume"] = row["terms"]["snubber_pads"]
        row["design_table_role"] = "RCD snubber capacitance-to-CAD-footprint handoff"
        rows.append(row)
    return rows


def thermal_robin_cooling_plate_reference_row(
    base_x,
    base_y,
    base_h,
    fin_count,
    fin_x,
    fin_y,
    fin_h,
    fin_pitch,
    device_pad_count,
    device_pad_x,
    device_pad_y,
    device_pad_h,
    device_pad_pitch,
    mount_hole_r,
    mount_hole_x,
    mount_hole_y,
    fin_y0=0.45,
    device_pad_y0=-1.45,
    label="thermal_robin_cooling_plate",
):
    """Analytic volume row for a readable convection-cooled drive plate.

    The geometry is a base plate, straight cooling fins, raised device pads,
    and four base-only mounting holes.  It is meant as a public-safe bridge
    between CAD volume checks and thermal Robin-boundary examples: the cooling
    surface can change later, but the solid volume should already be measurable
    before meshing or applying convection data.
    """

    base_x = float(base_x)
    base_y = float(base_y)
    base_h = float(base_h)
    fin_count = int(fin_count)
    fin_x = float(fin_x)
    fin_y = float(fin_y)
    fin_h = float(fin_h)
    fin_pitch = float(fin_pitch)
    device_pad_count = int(device_pad_count)
    device_pad_x = float(device_pad_x)
    device_pad_y = float(device_pad_y)
    device_pad_h = float(device_pad_h)
    device_pad_pitch = float(device_pad_pitch)
    mount_hole_r = float(mount_hole_r)
    mount_hole_x = abs(float(mount_hole_x))
    mount_hole_y = abs(float(mount_hole_y))
    fin_y0 = float(fin_y0)
    device_pad_y0 = float(device_pad_y0)

    positives = (
        base_x, base_y, base_h, fin_x, fin_y, fin_h, fin_pitch,
        device_pad_x, device_pad_y, device_pad_h, device_pad_pitch,
        mount_hole_r, mount_hole_x, mount_hole_y,
    )
    if any(value <= 0.0 for value in positives):
        raise ValueError("all dimensions, pitches, and radii must be positive")
    if fin_count < 1 or device_pad_count < 1:
        raise ValueError("fin_count and device_pad_count must be positive")
    if fin_x >= base_x:
        raise ValueError("fins must fit inside the base length")
    if device_pad_y >= base_y:
        raise ValueError("device pads must fit inside the base width")
    if mount_hole_x + mount_hole_r >= base_x / 2.0:
        raise ValueError("mount hole x location must fit inside the base")
    if mount_hole_y + mount_hole_r >= base_y / 2.0:
        raise ValueError("mount hole y location must fit inside the base")

    fin_span_y = (fin_count - 1) * fin_pitch + fin_y
    fin_min_y = fin_y0 - fin_span_y / 2.0
    fin_max_y = fin_y0 + fin_span_y / 2.0
    if fin_min_y <= -base_y / 2.0 or fin_max_y >= base_y / 2.0:
        raise ValueError("fin band must fit inside the base width")
    device_span_x = (device_pad_count - 1) * device_pad_pitch + device_pad_x
    if device_span_x >= base_x:
        raise ValueError("device pads must fit across the base")
    device_min_y = device_pad_y0 - device_pad_y / 2.0
    device_max_y = device_pad_y0 + device_pad_y / 2.0
    if device_min_y <= -base_y / 2.0 or device_max_y >= base_y / 2.0:
        raise ValueError("device pads must fit inside the base width at device_pad_y0")

    if not (device_max_y < fin_min_y or fin_max_y < device_min_y):
        raise ValueError("device pads must stay outside the fin band for additive volume")

    hole_centers = [
        (sx * mount_hole_x, sy * mount_hole_y)
        for sx in (-1.0, 1.0)
        for sy in (-1.0, 1.0)
    ]

    def _hole_clears_rectangle(hole_x, hole_y, rect_x, rect_y, rect_w, rect_h):
        return (
            abs(hole_x - rect_x) > rect_w / 2.0 + mount_hole_r
            or abs(hole_y - rect_y) > rect_h / 2.0 + mount_hole_r
        )

    fin_centers = [
        (0.0, fin_y0 + (index - (fin_count - 1) / 2.0) * fin_pitch)
        for index in range(fin_count)
    ]
    device_centers = [
        ((index - (device_pad_count - 1) / 2.0) * device_pad_pitch, device_pad_y0)
        for index in range(device_pad_count)
    ]
    for hole_x, hole_y in hole_centers:
        for rect_x, rect_y in fin_centers:
            if not _hole_clears_rectangle(hole_x, hole_y, rect_x, rect_y, fin_x, fin_y):
                raise ValueError("mount holes must not intersect fins")
        for rect_x, rect_y in device_centers:
            if not _hole_clears_rectangle(
                hole_x, hole_y, rect_x, rect_y, device_pad_x, device_pad_y
            ):
                raise ValueError("mount holes must not intersect device pads")

    base_volume = base_x * base_y * base_h
    fin_volume = fin_count * fin_x * fin_y * fin_h
    device_pad_volume = device_pad_count * device_pad_x * device_pad_y * device_pad_h
    mount_hole_volume = 4.0 * math.pi * mount_hole_r * mount_hole_r * base_h
    top_h = max(fin_h, device_pad_h)
    volume = base_volume + fin_volume + device_pad_volume - mount_hole_volume
    bbox = {
        "min": [-base_x / 2.0, -base_y / 2.0, -base_h / 2.0],
        "max": [base_x / 2.0, base_y / 2.0, base_h / 2.0 + top_h],
        "center": [0.0, 0.0, top_h / 2.0],
        "size": [base_x, base_y, base_h + top_h],
        "diagonal": math.sqrt(base_x * base_x + base_y * base_y + (base_h + top_h) ** 2),
    }
    fin_device_gap = fin_min_y - device_max_y if device_max_y < fin_min_y else device_min_y - fin_max_y
    return {
        "name": str(label),
        "volume": volume,
        "bounding_box": bbox,
        "terms": {
            "base": base_volume,
            "straight_cooling_fins": fin_volume,
            "device_pads": device_pad_volume,
            "four_mount_holes": -mount_hole_volume,
        },
        "counts": {
            "fins": fin_count,
            "device_pads": device_pad_count,
            "mount_holes": 4,
        },
        "clearances": {
            "mount_hole_to_x_edge": base_x / 2.0 - mount_hole_x - mount_hole_r,
            "mount_hole_to_y_edge": base_y / 2.0 - mount_hole_y - mount_hole_r,
            "fin_gap": fin_pitch - fin_y,
            "device_pad_gap": device_pad_pitch - device_pad_x,
            "fin_band_to_negative_y_edge": fin_min_y + base_y / 2.0,
            "fin_band_to_positive_y_edge": base_y / 2.0 - fin_max_y,
            "device_span_to_x_edge": (base_x - device_span_x) / 2.0,
            "device_pad_to_y_edge": base_y / 2.0 - max(abs(device_min_y), abs(device_max_y)),
            "device_pad_to_fin_band": fin_device_gap,
        },
        "parameters": {
            "fin_y0": fin_y0,
            "device_pad_y0": device_pad_y0,
            "fin_pitch": fin_pitch,
            "device_pad_pitch": device_pad_pitch,
        },
        "policy": "analytic_thermal_robin_cooling_plate_volume_reference",
    }


def motor_housing_radial_fin_reference_row(
    inner_radius,
    outer_radius,
    length,
    fin_count,
    fin_height,
    fin_thickness,
    label="motor_housing_radial_fins",
):
    """Analytic multi-body row for a cylindrical motor housing and radial fins.

    The sleeve and fins are intentionally separate tangent bodies.  This keeps
    the component basis readable for thermal contacts and makes the STEP
    round-trip check independent of boolean healing.  Total area includes all
    surfaces of the separate bodies; a later tied/contact thermal model decides
    how shared interfaces are treated.
    """
    inner_radius = float(inner_radius)
    outer_radius = float(outer_radius)
    length = float(length)
    fin_count = int(fin_count)
    fin_height = float(fin_height)
    fin_thickness = float(fin_thickness)
    if any(value <= 0.0 for value in (
        inner_radius, outer_radius, length, fin_height, fin_thickness
    )):
        raise ValueError("radii, length, fin height, and fin thickness must be positive")
    if outer_radius <= inner_radius:
        raise ValueError("outer_radius must be larger than inner_radius")
    if fin_count < 3:
        raise ValueError("fin_count must be at least 3")
    chord = 2.0 * outer_radius * math.sin(math.pi / fin_count)
    if fin_thickness >= chord:
        raise ValueError("fin_thickness must leave separated radial fins")

    sleeve_volume = math.pi * (outer_radius**2 - inner_radius**2) * length
    one_fin_volume = fin_height * fin_thickness * length
    volume = sleeve_volume + fin_count * one_fin_volume
    sleeve_area = (
        2.0 * math.pi * (outer_radius + inner_radius) * length
        + 2.0 * math.pi * (outer_radius**2 - inner_radius**2)
    )
    one_fin_area = 2.0 * (
        fin_height * fin_thickness
        + fin_height * length
        + fin_thickness * length
    )
    area = sleeve_area + fin_count * one_fin_area
    envelope_radius = outer_radius + fin_height
    return {
        "name": str(label),
        "volume": volume,
        "area": area,
        "body_count": fin_count + 1,
        "units": {"length": "mm", "area": "mm^2", "volume": "mm^3"},
        "terms": {
            "sleeve_volume": sleeve_volume,
            "one_fin_volume": one_fin_volume,
            "all_fins_volume": fin_count * one_fin_volume,
            "sleeve_surface_area": sleeve_area,
            "one_fin_surface_area": one_fin_area,
            "all_fins_surface_area": fin_count * one_fin_area,
        },
        "parameters": {
            "inner_radius": inner_radius,
            "outer_radius": outer_radius,
            "length": length,
            "fin_count": fin_count,
            "fin_height": fin_height,
            "fin_thickness": fin_thickness,
        },
        "bounding_box": {
            "min": [-envelope_radius, -envelope_radius, 0.0],
            "max": [envelope_radius, envelope_radius, length],
            "size": [2.0 * envelope_radius, 2.0 * envelope_radius, length],
        },
        "roundtrip_tolerances": {
            "curved_step_volume_rtol": 1.0e-4,
            "surface_area_rtol": 1.0e-10,
            "body_count_exact": True,
        },
        "policy": "analytic_motor_housing_radial_fin_multibody_reference",
    }


def v_type_ipm_rotor_coupon_reference_row(
    coupon_x,
    coupon_y,
    thickness,
    magnet_slot_length,
    magnet_slot_width,
    magnet_slot_angle_deg,
    magnet_slot_center_x,
    magnet_slot_center_y,
    bore_radius,
    label="v_type_ipm_rotor_coupon",
):
    """Analytic volume row for a readable V-type IPM rotor coupon.

    The public CAD gate is a rectangular rotor-lamination coupon with two
    through rectangular magnet pockets mirrored about the y-axis and rotated
    into a V, plus one central through bore.  The pockets and bore are required
    to stay inside the coupon and not overlap, so the volume is exactly
    ``coupon - 2*pocket - bore``.  This is a compact pre-FEM geometry contract
    before a full motor rotor sector is meshed.
    """

    coupon_x = float(coupon_x)
    coupon_y = float(coupon_y)
    thickness = float(thickness)
    magnet_slot_length = float(magnet_slot_length)
    magnet_slot_width = float(magnet_slot_width)
    angle = float(magnet_slot_angle_deg)
    magnet_slot_center_x = abs(float(magnet_slot_center_x))
    magnet_slot_center_y = float(magnet_slot_center_y)
    bore_radius = float(bore_radius)
    positives = (
        coupon_x, coupon_y, thickness, magnet_slot_length, magnet_slot_width,
        magnet_slot_center_x, bore_radius,
    )
    if any(value <= 0.0 for value in positives):
        raise ValueError("all dimensions except magnet_slot_center_y must be positive")
    if not (0.0 < abs(angle) < 90.0):
        raise ValueError("magnet_slot_angle_deg must be between 0 and 90 degrees")
    half_x = coupon_x / 2.0
    half_y = coupon_y / 2.0
    theta = math.radians(abs(angle))
    pocket_half_x = 0.5 * (
        magnet_slot_length * math.cos(theta)
        + magnet_slot_width * math.sin(theta)
    )
    pocket_half_y = 0.5 * (
        magnet_slot_length * math.sin(theta)
        + magnet_slot_width * math.cos(theta)
    )
    if magnet_slot_center_x + pocket_half_x >= half_x:
        raise ValueError("magnet pockets must fit inside coupon_x")
    if abs(magnet_slot_center_y) + pocket_half_y >= half_y:
        raise ValueError("magnet pockets must fit inside coupon_y")
    if bore_radius >= min(half_x, half_y):
        raise ValueError("bore must fit inside the coupon")
    pocket_inner_x = magnet_slot_center_x - pocket_half_x
    pocket_inner_y_clearance = abs(magnet_slot_center_y) - pocket_half_y
    if pocket_inner_x <= bore_radius:
        raise ValueError("magnet pockets must not overlap the bore")
    if pocket_inner_y_clearance < -bore_radius:
        raise ValueError("magnet pockets must clear the bore in y")
    if 2.0 * magnet_slot_center_x <= 2.0 * pocket_half_x:
        raise ValueError("mirrored magnet pockets must not overlap each other")

    base_volume = coupon_x * coupon_y * thickness
    pocket_volume = 2.0 * magnet_slot_length * magnet_slot_width * thickness
    bore_volume = math.pi * bore_radius * bore_radius * thickness
    volume = base_volume - pocket_volume - bore_volume
    bbox = {
        "min": [-half_x, -half_y, -thickness / 2.0],
        "max": [half_x, half_y, thickness / 2.0],
        "center": [0.0, 0.0, 0.0],
        "size": [coupon_x, coupon_y, thickness],
        "diagonal": math.sqrt(coupon_x * coupon_x + coupon_y * coupon_y + thickness * thickness),
    }
    return {
        "name": str(label),
        "volume": volume,
        "bounding_box": bbox,
        "terms": {
            "coupon": base_volume,
            "two_v_magnet_pockets": -pocket_volume,
            "central_bore": -bore_volume,
        },
        "counts": {
            "magnet_pockets": 2,
            "bore": 1,
        },
        "clearances": {
            "pocket_to_x_edge": half_x - magnet_slot_center_x - pocket_half_x,
            "pocket_to_y_edge": half_y - abs(magnet_slot_center_y) - pocket_half_y,
            "pocket_to_bore_x": pocket_inner_x - bore_radius,
            "mirrored_pocket_gap": 2.0 * (magnet_slot_center_x - pocket_half_x),
        },
        "parameters": {
            "magnet_slot_angle_deg": angle,
            "magnet_slot_center_x": magnet_slot_center_x,
            "magnet_slot_center_y": magnet_slot_center_y,
            "magnet_slot_length": magnet_slot_length,
            "magnet_slot_width": magnet_slot_width,
            "bore_radius": bore_radius,
        },
        "policy": "analytic_v_type_ipm_rotor_coupon_volume_reference",
    }


def shape_parameter_sweep_summary(
    rows,
    parameter_key,
    metric_keys=("volume", "area"),
    limits_by_metric=None,
    monotonic_tolerance=1.0e-12,
):
    """Summarize a CAD parameter sweep from measurement rows.

    Each row should contain ``parameter_key`` and the requested metric keys
    (for example ``volume`` and ``area`` from :func:`shape_measurement_row`).
    The summary sorts by the parameter value and reports monotonicity, extrema,
    spans, and optional min/max constraint violations for each metric.  It is a
    small pre-mesh design table before geometry rows are sent to meshing,
    validation, or optimization.
    """

    rows = [dict(row) for row in rows]
    if not rows:
        raise ValueError("rows must not be empty")
    tolerance = float(monotonic_tolerance)
    if tolerance < 0.0:
        raise ValueError("monotonic_tolerance must be >= 0")
    metric_keys = tuple(str(key) for key in metric_keys)
    if not metric_keys:
        raise ValueError("metric_keys must not be empty")
    limits_by_metric = limits_by_metric or {}

    for row in rows:
        if parameter_key not in row:
            raise KeyError(f"missing parameter {parameter_key!r}")
        row[parameter_key] = float(row[parameter_key])
        if not math.isfinite(row[parameter_key]):
            raise ValueError("parameter values must be finite")
        for key in metric_keys:
            if key not in row:
                raise KeyError(f"missing metric {key!r}")
            row[key] = float(row[key])
            if not math.isfinite(row[key]):
                raise ValueError("metric values must be finite")

    rows.sort(key=lambda row: row[parameter_key])
    parameter_values = [row[parameter_key] for row in rows]
    duplicate_parameters = len(set(parameter_values)) != len(parameter_values)
    parameter_strict = all(parameter_values[i] < parameter_values[i + 1] for i in range(len(rows) - 1))
    metric_rows = []
    violations = []

    for key in metric_keys:
        values = [row[key] for row in rows]
        deltas = [values[i + 1] - values[i] for i in range(len(values) - 1)]
        min_index = min(range(len(values)), key=lambda i: values[i])
        max_index = max(range(len(values)), key=lambda i: values[i])
        limits = limits_by_metric.get(key, {})
        lower = limits.get("min") if isinstance(limits, dict) else None
        upper = limits.get("max") if isinstance(limits, dict) else None
        lower = None if lower is None else float(lower)
        upper = None if upper is None else float(upper)
        if lower is not None and not math.isfinite(lower):
            raise ValueError(f"lower limit for metric {key!r} must be finite")
        if upper is not None and not math.isfinite(upper):
            raise ValueError(f"upper limit for metric {key!r} must be finite")
        metric_violations = []
        for row, value in zip(rows, values):
            if lower is not None and value < lower - tolerance:
                metric_violations.append({
                    "parameter": row[parameter_key],
                    "metric": key,
                    "value": value,
                    "limit": lower,
                    "kind": "below_min",
                })
            if upper is not None and value > upper + tolerance:
                metric_violations.append({
                    "parameter": row[parameter_key],
                    "metric": key,
                    "value": value,
                    "limit": upper,
                    "kind": "above_max",
                })
        violations.extend(metric_violations)
        min_value = values[min_index]
        max_value = values[max_index]
        metric_rows.append({
            "metric": key,
            "min": min_value,
            "min_parameter": parameter_values[min_index],
            "max": max_value,
            "max_parameter": parameter_values[max_index],
            "span": max_value - min_value,
            "relative_span": (
                (max_value - min_value) / abs(min_value)
                if min_value != 0.0
                else (0.0 if max_value == 0.0 else math.inf)
            ),
            "first": values[0],
            "last": values[-1],
            "delta_first_to_last": values[-1] - values[0],
            "monotonic_non_decreasing": all(delta >= -tolerance for delta in deltas),
            "monotonic_non_increasing": all(delta <= tolerance for delta in deltas),
            "min_step_delta": min(deltas) if deltas else 0.0,
            "max_step_delta": max(deltas) if deltas else 0.0,
            "limits": {"min": lower, "max": upper},
            "constraint_violation_count": len(metric_violations),
        })

    ok = not duplicate_parameters and not violations
    issues = []
    if duplicate_parameters:
        issues.append("duplicate parameter values")
    if violations:
        issues.append("at least one metric is outside requested limits")
    return {
        "policy": "build123d_parameter_sweep_measurements_are_sorted_and_audited",
        "parameter_key": str(parameter_key),
        "n_cases": len(rows),
        "parameter_values": parameter_values,
        "parameter_min": parameter_values[0],
        "parameter_max": parameter_values[-1],
        "parameter_strictly_increasing": parameter_strict,
        "duplicate_parameter_values": duplicate_parameters,
        "metric_rows": metric_rows,
        "constraint_violations": violations,
        "constraint_violation_count": len(violations),
        "status": "ok" if ok else "needs_attention",
        "ok_for_design_table": ok,
        "issues": issues,
        "rows": rows,
    }


# =====================================================================================================
# Generic solid-modelling operations -- the constructor / local-modification / array verbs a full 3D
# modeller exposes, wrapped as clean, labelled, CAE-safe build123d helpers.
# =====================================================================================================
def swept(profile, path, label="swept"):
    r"""**Sweep** a 2D ``profile`` (a face / sketch such as ``Circle(r)`` or ``Rectangle(w, h)``) along a
    ``path`` (a wire / edge -- a line, arc, spline or :class:`~build123d.Helix`), keeping the profile
    perpendicular to the path tangent at the start.  The "sweep curve" verb -- pipes, coil conductors,
    swept beams.  Returns a single labelled :class:`~build123d.Solid`."""
    sec = Plane(origin=path @ 0.0, z_dir=path % 0.0) * profile
    s = sweep(sec, path=path).solid()
    s.label = label
    return s


def revolved(profile, axis=None, angle=360.0, label="revolved"):
    r"""**Revolve** (spin) a 2D ``profile`` about ``axis`` by ``angle`` degrees -- bodies of revolution
    (shafts, vases, toroidal cores, pulleys).  Returns a labelled :class:`~build123d.Solid`."""
    if axis is None:
        _require_build123d()
        axis = Axis.Y
    s = revolve(profile, axis, angle).solid()
    s.label = label
    return s


def lofted(sections, label="lofted"):
    r"""**Loft** (blend) a list of 2D ``sections`` (faces, positioned at their stations) into a solid that
    interpolates between them -- transitions, blades, ducts, frusta.  Returns a labelled
    :class:`~build123d.Solid`."""
    s = loft(sections).solid()
    s.label = label
    return s


def coil(profile, pitch, height, radius, label="coil"):
    r"""A **helical coil**: sweep a 2D ``profile`` (e.g. ``Circle(wire_radius)``) along a helix of the
    given ``pitch`` (axial advance per turn), ``height`` and ``radius``.  Solenoids, springs, helical
    conductors.  Returns a labelled :class:`~build123d.Solid` (volume ``= area(profile) * helix length``)."""
    return swept(profile, Helix(pitch=pitch, height=height, radius=radius), label=label)


def helix_centerline_length(radius, pitch, height):
    r"""Centreline length of a constant-radius helix.

    ``pitch`` is the axial advance per turn and ``height`` is the axial span, so
    ``turns = height / pitch`` and

        ``length = sqrt(height^2 + (2*pi*radius*turns)^2)``.

    This is the analytic length behind :func:`coil`, useful before CAD generation for resistance,
    copper volume and mesh-size estimates.
    """
    radius = float(radius)
    pitch = float(pitch)
    height = float(height)
    if radius < 0.0:
        raise ValueError("radius must be >= 0")
    if pitch <= 0.0:
        raise ValueError("pitch must be > 0")
    if height < 0.0:
        raise ValueError("height must be >= 0")
    turns = height / pitch
    return math.hypot(height, 2.0 * math.pi * radius * turns)


def round_wire_helix_metrics(radius, wire_radius, pitch, height):
    r"""Pre-CAD metrics for a round wire swept on a constant-radius helix.

    Returns turns, centreline length, circular cross-section area, conductor volume, and
    ``resistance_per_resistivity = length / area``.  Multiply the last value by material resistivity to
    get DC resistance before skin/proximity corrections.
    """
    wire_radius = float(wire_radius)
    if wire_radius <= 0.0:
        raise ValueError("wire_radius must be > 0")
    length = helix_centerline_length(radius, pitch, height)
    area = math.pi * wire_radius * wire_radius
    return {
        "radius": float(radius),
        "wire_radius": wire_radius,
        "pitch": float(pitch),
        "height": float(height),
        "turns": float(height) / float(pitch),
        "centerline_length": length,
        "cross_section_area": area,
        "conductor_volume": area * length,
        "resistance_per_resistivity": length / area,
    }


def strut(p0, p1, radius, label="strut"):
    r"""A round **strut** (cylinder) of given ``radius`` between two arbitrary 3D points ``p0`` and
    ``p1`` -- the member primitive for trusses, space frames, spokes and lattices."""
    return swept(Circle(radius), Line(p0, p1), label=label)


def thicken(profile, thickness, label="sheet"):
    r"""**Thicken** a planar ``profile`` (face / sketch) into a solid of the given ``thickness``,
    symmetric about the profile plane -- the "thicken sheet" verb (plates, laminations, membranes)."""
    s = extrude(profile, amount=thickness / 2.0, both=True).solid()
    s.label = label
    return s


def draft_extrude(profile, height, taper_deg, label="draft"):
    r"""**Extrude with draft**: extrude a 2D ``profile`` by ``height`` while tapering its walls by
    ``taper_deg`` (a positive angle shrinks the top) -- moulded / cast features, frusta."""
    s = extrude(profile, amount=height, taper=taper_deg).solid()
    s.label = label
    return s


def shell(solid, thickness, openings=None, label="shell"):
    r"""**Shell** (hollow out) a ``solid`` to a wall of ``thickness``, removing material inward.  With
    ``openings`` (a face or list of faces of ``solid``) the shell is OPEN on those faces (a cup / box);
    without, it is a CLOSED hollow shell.  The "shell solid" verb."""
    if openings is None:
        sh = (solid - offset(solid, amount=-abs(thickness)))
    else:
        sh = offset(solid, amount=-abs(thickness), openings=openings)
    sh = sh.solid()
    sh.label = label
    return sh


def fillet_edges(solid, radius, edge_filter=None, label="filleted"):
    r"""**Blend (fillet)** the edges of ``solid`` with the given ``radius``.  ``edge_filter`` is an
    optional callable ``edges -> edges`` to pick a subset (e.g. ``lambda e: e.filter_by(Axis.Z)`` for the
    vertical edges only); without it, every edge is rounded."""
    edges = solid.edges() if edge_filter is None else edge_filter(solid.edges())
    f = fillet(edges, radius)
    try:
        f.label = label
    except Exception:
        pass
    return f


def chamfer_edges(solid, length, edge_filter=None, label="chamfered"):
    r"""**Chamfer** the edges of ``solid`` by the given ``length``.  ``edge_filter`` is an optional
    callable ``edges -> edges`` to pick a subset; without it, every edge is chamfered."""
    edges = solid.edges() if edge_filter is None else edge_filter(solid.edges())
    c = chamfer(edges, length)
    try:
        c.label = label
    except Exception:
        pass
    return c


def grid_array(part, nx, ny, dx, dy, label=None, label_fmt="{base}_{i:02d}_{j:02d}"):
    r"""An ``nx`` x ``ny`` **rectangular grid** of copies of ``part`` at pitches ``dx`` (x) and ``dy``
    (y), returned as a labelled :class:`~build123d.Compound` -- the 2D "translate with copies" verb
    (pin grids, hole patterns, fin / standoff arrays)."""
    if nx < 1 or ny < 1:
        raise ValueError("nx, ny must be >= 1")
    base = label or (part.label or "part")
    children = []
    for i in range(nx):
        for j in range(ny):
            c = Pos(i * dx, j * dy, 0) * part
            c.label = label_fmt.format(base=base, i=i, j=j)
            children.append(c)
    return Compound(children=children, label=base + "_grid")


def path_array(part, path, count, label=None, label_fmt="{base}_{k:02d}"):
    r"""``count`` copies of ``part`` placed at equal parameter steps along ``path`` (a wire / edge) --
    the "array along a curve" verb (beads on a wire, bolts along a slot, stations along a spline)."""
    if count < 1:
        raise ValueError("count must be >= 1")
    base = label or (part.label or "part")
    children = []
    for k in range(count):
        t = k / (count - 1) if count > 1 else 0.0
        c = Pos(*tuple(path @ t)) * part
        c.label = label_fmt.format(base=base, k=k)
        children.append(c)
    return Compound(children=children, label=base + "_patharray")


# ---- boolean / slice / sheet-metal verbs -----------------------------------------------------------
def fuse(*parts, label="fused"):
    r"""**Boolean union** (add): fuse ``parts`` into one body, relabelled.  Overlapping inputs merge into
    a single solid; disjoint inputs remain separate solids in the returned body."""
    if not parts:
        raise ValueError("fuse needs at least one part")
    r = parts[0]
    for p in parts[1:]:
        r = r + p
    try:
        r.label = label
    except Exception:
        pass
    return r


def cut(base, *tools, label="cut"):
    r"""**Boolean subtract**: ``base`` minus each of ``tools`` (drill holes / remove material),
    relabelled."""
    r = base
    for t in tools:
        r = r - t
    try:
        r.label = label
    except Exception:
        pass
    return r


def common(*parts, label="common"):
    r"""**Boolean intersect**: the common (overlapping) volume of ``parts``, relabelled."""
    if not parts:
        raise ValueError("common needs at least one part")
    r = parts[0]
    for p in parts[1:]:
        r = r & p
    try:
        r.label = label
    except Exception:
        pass
    return r


def slice_solid(solid, plane=None, keep="top", label="slice"):
    r"""**Slice** a solid by a ``plane`` and keep one side -- ``keep`` is ``"top"`` (the +normal side),
    ``"bottom"`` (the -normal side) or ``"both"`` (returns both halves).  The plane-cut / split verb
    (half models on a symmetry plane, sectioning)."""
    if plane is None:
        _require_build123d()
        plane = Plane.XY
    keepmap = {"top": Keep.TOP, "bottom": Keep.BOTTOM, "both": Keep.BOTH}
    if keep not in keepmap:
        raise ValueError("keep must be 'top', 'bottom' or 'both'")
    r = split(solid, bisect_by=plane, keep=keepmap[keep])
    try:
        r.label = label
    except Exception:
        pass
    return r


def bend_sheet(leg1, leg2, width, thickness, angle_deg, radius, n_arc=24, label="bent"):
    r"""A **sheet-metal bend**: a strip of cross-section ``thickness`` (in the bend plane) x ``width``
    (out of plane) that runs a flat leg of length ``leg1``, a cylindrical bend of inner ``radius``
    through ``angle_deg``, then a flat leg of length ``leg2``.  Built by sweeping the cross-section along
    the neutral-axis path (straight -> arc -> straight), so it is one clean labelled solid -- the
    sheet-metal "bend" verb (brackets, flanges, folded chassis)."""
    a = math.radians(angle_deg)
    rn = radius + thickness / 2.0
    pts = [(-leg1, 0.0, 0.0), (-leg1 * 0.5, 0.0, 0.0), (0.0, 0.0, 0.0)]
    for k in range(1, n_arc + 1):
        t = a * k / n_arc
        pts.append((rn * math.sin(t), rn * (1.0 - math.cos(t)), 0.0))
    pa = pts[-1]
    pts.append((pa[0] + leg2 * math.cos(a), pa[1] + leg2 * math.sin(a), 0.0))
    return swept(Rectangle(thickness, width), Spline(*pts), label=label)


def fillet_varied(solid, specs, label="filleted"):
    r"""**Variable-radius fillet**: ``specs`` is a list of ``(edge_filter, radius)`` applied in turn, so
    different edge groups take different blend radii (e.g. a large radius on the load-bearing edges, a
    small one elsewhere).  Each ``edge_filter`` is a callable ``edges -> edges`` (e.g.
    ``lambda e: e.filter_by(Axis.Z)``); the fillets are applied sequentially (each operates on the
    already-filleted solid, so order edge-groups outer-to-inner)."""
    g = solid
    for edge_filter, radius in specs:
        g = fillet(edge_filter(g.edges()), radius)
    try:
        g.label = label
    except Exception:
        pass
    return g


def chamfer_varied(solid, specs, label="chamfered"):
    r"""**Variable chamfer**: ``specs`` is a list of ``(edge_filter, length)`` applied in turn -- the
    chamfer analogue of :func:`fillet_varied`."""
    g = solid
    for edge_filter, length in specs:
        g = chamfer(edge_filter(g.edges()), length)
    try:
        g.label = label
    except Exception:
        pass
    return g


def slice_array(solid, n, thickness, axis="z", label="slice"):
    r"""**Slice** a solid into ``n`` slabs of ``thickness`` along ``axis`` (``'x'`` / ``'y'`` / ``'z'``),
    centred on the solid, returned as a labelled :class:`~build123d.Compound` of the ``n`` pieces -- the
    lamination-stack / parting verb (split an iron core into laminations, a bar into segments).  Each
    piece is the solid intersected with one slab; choose ``n * thickness`` to span the solid's extent."""
    if n < 1:
        raise ValueError("n must be >= 1")
    if axis not in ("x", "y", "z"):
        raise ValueError("axis must be 'x', 'y' or 'z'")
    bb = solid.bounding_box()
    big = max(bb.size.X, bb.size.Y, bb.size.Z) * 3.0
    c = solid.center()
    total = n * thickness
    pieces = []
    for k in range(n):
        off = -total / 2.0 + (k + 0.5) * thickness
        if axis == "z":
            slab = Pos(c.X, c.Y, c.Z + off) * Box(big, big, thickness)
        elif axis == "x":
            slab = Pos(c.X + off, c.Y, c.Z) * Box(thickness, big, big)
        else:
            slab = Pos(c.X, c.Y + off, c.Z) * Box(big, thickness, big)
        piece = solid & slab
        try:
            piece = piece.solid()
        except Exception:
            pass
        piece.label = f"{label}_{k:02d}"
        pieces.append(piece)
    return Compound(children=pieces, label=label + "_stack")
