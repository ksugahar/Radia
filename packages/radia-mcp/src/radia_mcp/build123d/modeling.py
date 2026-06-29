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

from build123d import (Axis, Box, BuildLine, BuildSketch, CenterArc, Circle, Compound, Cylinder, Helix,
                       Keep, Line, Mode, Plane, Pos, Rectangle, RectangleRounded, Spline, chamfer,
                       extrude, fillet, loft, make_face, offset, revolve, split, sweep)

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
           "box_face_vector_area_rows", "box_face_pressure_force_rows",
           "box_face_pressure_moment_rows", "box_face_pressure_resultant_summary",
           "box_face_traction_moment_rows",
           "compare_boundary_vector_area_rows",
           "compare_shape_measurement_rows", "shape_measurement_comparison_summary",
           "compare_shape_volume_rows", "shape_volume_crosscheck_summary",
           "shape_name_identity_gate",
           "shape_role_metadata_gate",
           "shape_transition_role_metadata_gate",
           "shape_mass_property_crosscheck_summary",
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


def polar_array(part, count, total_angle=360.0, axis=Axis.Z, label=None, label_fmt="{base}_{k:02d}"):
    r"""``count`` rotated copies of ``part`` about ``axis``, returned as a labelled
    :class:`~build123d.Compound`.  A full ``360`` deg ring spaces the copies by ``360/count``; a
    partial fan (``total_angle < 360``) spaces them by ``total_angle/(count-1)`` so the first and last
    copies sit at the fan ends.  The "rotate with copies" / segmented-array verb.
    """
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


def mirrored(part, about=Plane.XZ, keep_original=True, label=None):
    r"""The mirror image of ``part`` across plane ``about`` (default ``Plane.XZ``); with
    ``keep_original`` (default) returns a :class:`~build123d.Compound` of the original + its mirror --
    the symmetry-completion verb (build a quarter/half model, then mirror to whole).
    """
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
            return [("external_cad", measured_sets)]
        return [(str(label), rows) for label, rows in measured_sets.items()]
    normalized = []
    for index, item in enumerate(list(measured_sets), start=1):
        if isinstance(item, dict) and "rows" in item:
            label = str(item.get("source") or item.get("label") or f"external_cad_{index}")
            rows = item["rows"]
        else:
            label = f"external_cad_{index}"
            rows = item
        normalized.append((label, rows))
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
    for label, rows in _normalize_volume_measurement_sets(measured_sets):
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
    require_positive_volume=True,
    source_label="build123d",
):
    """Check CAD-side metadata for a future hex-to-tet transition handoff.

    build123d does not create Cubit pyramid elements, but it can preserve the
    solver intent before STEP/Cubit handoff: which body is the hex-led region,
    which body is the tet region, and which body is the transition envelope.
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
        volume = item.get("volume")
        normalized.append({
            "name": name,
            "role": role,
            "material": material,
            "transition_kind": transition_kind,
            "connected_roles": connected_roles,
            "volume": None if volume is None else float(volume),
            "source_row": item,
        })

    names = [row["name"] for row in normalized if row["name"]]
    roles = {row["role"] for row in normalized if row["role"]}
    required_role_set = {str(role).strip() for role in required_roles if str(role).strip()}
    connected_role_set = {str(role).strip() for role in required_connected_roles if str(role).strip()}
    transition_rows = [row for row in normalized if row["role"] == transition_role]
    transition_kinds = {row["transition_kind"] for row in transition_rows if row["transition_kind"]}
    connected_roles_union = {
        role for row in transition_rows for role in row["connected_roles"] if role
    }
    rows_missing_material = [
        row["name"] for row in normalized if row["name"] and not row["material"]
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
        "rows_missing_material": sorted(rows_missing_material),
        "rows_missing_volume": sorted(rows_missing_volume),
        "rows_nonpositive_volume": sorted(rows_nonpositive_volume),
        "checks": checks,
        "status": "ok" if all(checks.values()) else "needs_attention",
        "version_note": (
            "Run this on build123d assembly metadata before handing a future "
            "hex+tet model to Cubit; the pyramid is a mesh transition contract, "
            "not a build123d primitive requirement."
        ),
    }


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
    inventory = shape_measurement_inventory_summary(reference)
    sets = []
    all_rows = []
    identity_gates = []
    for label, rows in _normalize_shape_measurement_sets(measured_sets):
        rows_list = list(rows)
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
    }
    issues = []
    if not checks["all_reference_shapes_valid"]:
        issues.append("at least one reference build123d shape is invalid")
    if not checks["all_sources_present_and_within_tolerance"]:
        issues.append("at least one CAD source row is missing or outside tolerance")
    if not checks["all_sources_preserve_named_shape_identity"]:
        issues.append("at least one CAD source has missing, extra, duplicate, or unnamed shapes")
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


def revolved(profile, axis=Axis.Y, angle=360.0, label="revolved"):
    r"""**Revolve** (spin) a 2D ``profile`` about ``axis`` by ``angle`` degrees -- bodies of revolution
    (shafts, vases, toroidal cores, pulleys).  Returns a labelled :class:`~build123d.Solid`."""
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


def slice_solid(solid, plane=Plane.XY, keep="top", label="slice"):
    r"""**Slice** a solid by a ``plane`` and keep one side -- ``keep`` is ``"top"`` (the +normal side),
    ``"bottom"`` (the -normal side) or ``"both"`` (returns both halves).  The plane-cut / split verb
    (half models on a symmetry plane, sectioning)."""
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
