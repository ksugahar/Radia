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

import math

from build123d import (Axis, Box, BuildLine, BuildSketch, CenterArc, Circle, Compound, Cylinder, Helix,
                       Keep, Line, Mode, Plane, Pos, Rectangle, RectangleRounded, Spline, chamfer,
                       extrude, fillet, loft, make_face, offset, revolve, split, sweep)

__all__ = ["annular_segment", "tube", "racetrack_coil", "polar_array", "linear_array",
           "mirrored", "assembly", "shape_envelope_row", "enclosing_box",
           "enclosure_clearance_row", "enclosure_difference_region",
           "shape_measurement_row", "shape_measurement_rows",
           "box_face_vector_area_rows", "box_face_pressure_force_rows",
           "box_face_pressure_moment_rows", "box_face_pressure_resultant_summary",
           "box_face_traction_moment_rows",
           "compare_boundary_vector_area_rows",
           "compare_shape_measurement_rows", "shape_measurement_comparison_summary",
           "shape_measurement_inventory_summary", "worst_shape_measurement_comparison_rows",
           "shape_measurement_health_summary",
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
