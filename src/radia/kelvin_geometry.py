"""kelvin_geometry.py

Layer 1 of the Kelvin helper API (api_plan.md). Encapsulates the
boilerplate to assemble a Sugahara two-sphere Kelvin geometry: a set of
inner physical sub-shapes glued together + an offset Kelvin exterior
sphere, with their boundaries identified via a periodic BC and a GND
vertex at the Kelvin exterior sphere center (image of physical
infinity).
"""

from __future__ import annotations

from netgen.occ import (Sphere, Vertex, Glue, OCCGeometry, Pnt,
                         IdentificationType)


def _iter_inner_shapes(inner):
    """Normalize ``inner`` to a non-empty list of OCC shapes."""
    if isinstance(inner, (list, tuple)):
        shapes = list(inner)
    else:
        shapes = [inner]
    if not shapes:
        raise ValueError("add_kelvin_exterior_domain: inner shape list is empty")
    return shapes


def _find_named_face(shapes, face_name):
    """Return (shape_index, face) for the first face matching ``face_name``."""
    for i, s in enumerate(shapes):
        for f in s.faces:
            if f.name == face_name:
                return i, f
    return None, None


def add_kelvin_exterior_domain(inner_shape, offset, R_K,
                                 inner_kelvin_face_name="kelvin_int",
                                 outer_kelvin_face_name="kelvin_ext",
                                 kelvin_mat="kelvin",
                                 gnd_vertex_name="GND",
                                 outer_maxh_factor=2.0,
                                 inner_maxh=None):
    """Glue ``inner_shape`` (or a list of inner sub-shapes) with an
    offset Kelvin exterior sphere.

    The Sugahara two-sphere convention:
      - inner physical sphere at ORIGIN, radius R_K (carrier of the
        physical model, e.g. coil + air). The caller supplies this as
        ``inner_shape``: either a single OCC shape (already confined
        to the sphere |r| <= R_K) OR a list/tuple of sub-shapes (for
        example ``[inner_air, coil_torus]``) which will be glued
        together here. Exactly one of those sub-shapes must carry an
        outward face named ``inner_kelvin_face_name`` (typically the
        air sub-shape). The helper locates that face and uses it for
        the periodic identification.
      - outer Kelvin sphere at ``offset``, radius R_K, material
        ``kelvin_mat``.
      - Periodic BC identifies the two sphere surfaces (Sugahara 2022).
      - GND vertex at the Kelvin sphere center maps to physical infinity.

    Args:
        inner_shape: OCC shape or list/tuple of OCC shapes.
        offset: 3-tuple, center of the outer Kelvin sphere.
        R_K: Kelvin sphere radius (same for both inner and outer).
        inner_kelvin_face_name, outer_kelvin_face_name: names used for
            the periodic identification and bcname tagging.
        kelvin_mat: material name on the outer sphere (default "kelvin").
        gnd_vertex_name: BBND name for the Dirichlet GND vertex.
        outer_maxh_factor: multiplier for the Kelvin sphere maxh
            relative to ``inner_maxh``. Default 2.0 (coarser).
        inner_maxh: maxh hint for deriving the outer sphere maxh. If
            None the outer sphere keeps whatever maxh the caller set.

    Returns:
        ``(geometry, info)`` where ``geometry`` is the glued OCC
        compound ready to feed to ``OCCGeometry(...)`` and ``info`` is
        a dict with keys ``inner_face``, ``outer_face``, ``outer_shape``,
        ``gnd_vertex``, ``offset``, ``R_K`` and ``inner_shapes``
        (the normalized list).
    """
    inner_shapes = _iter_inner_shapes(inner_shape)

    # --- Locate the inner Kelvin-interface face (must already be named).
    _, inner_kelvin_face = _find_named_face(inner_shapes,
                                              inner_kelvin_face_name)
    if inner_kelvin_face is None:
        raise ValueError(
            f"add_kelvin_exterior_domain: no face named "
            f"{inner_kelvin_face_name!r} found on the provided inner "
            f"shape(s); tag the outward sphere face before calling "
            f"(e.g. ``for f in inner_sphere.faces: f.name = "
            f"{inner_kelvin_face_name!r}``)")

    # --- Build the outer Kelvin sphere at offset. --------------------
    outer = Sphere(Pnt(*offset), R_K)
    if inner_maxh is not None:
        outer.maxh = inner_maxh * outer_maxh_factor
    for f in outer.faces:
        f.name = outer_kelvin_face_name
    outer.name = kelvin_mat

    # First (and only) face of a raw Sphere is the sphere surface.
    outer_kelvin_face = None
    for f in outer.faces:
        if f.name == outer_kelvin_face_name:
            outer_kelvin_face = f
            break
    if outer_kelvin_face is None:  # pragma: no cover - defensive
        raise RuntimeError("outer Kelvin sphere has no tagged face")

    # --- Periodic identification: inner outer face <-> outer sphere face.
    # Done BEFORE Glue so the identification is baked into the OCC
    # topology that Netgen sees.
    inner_kelvin_face.Identify(outer_kelvin_face, "kelvin_periodic",
                                IdentificationType.PERIODIC)

    # --- GND vertex at the Kelvin center (image of physical infinity).
    gnd = Vertex(Pnt(*offset))
    gnd.name = gnd_vertex_name

    geometry = Glue(list(inner_shapes) + [outer, gnd])
    info = {
        "inner_shapes": inner_shapes,
        "inner_face": inner_kelvin_face,
        "outer_face": outer_kelvin_face,
        "outer_shape": outer,
        "gnd_vertex": gnd,
        "offset": offset,
        "R_K": R_K,
    }
    return geometry, info
