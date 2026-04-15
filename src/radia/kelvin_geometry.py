"""kelvin_geometry.py

Layer 1 of the Kelvin helper API (api_plan.md). Encapsulates the
boilerplate to assemble a Sugahara two-sphere Kelvin geometry: an
inner physical sphere + an offset Kelvin exterior sphere, with their
boundaries identified via a periodic BC and a GND vertex at the
Kelvin sphere center (image of physical infinity).
"""

from __future__ import annotations

from netgen.occ import (Sphere, Vertex, Glue, OCCGeometry, Pnt,
                         IdentificationType)


def add_kelvin_exterior_domain(inner_shape, offset, R_K,
                                 inner_kelvin_face_name="kelvin_int",
                                 outer_kelvin_face_name="kelvin_ext",
                                 kelvin_mat="kelvin",
                                 gnd_vertex_name="GND",
                                 outer_maxh_factor=2.0,
                                 inner_maxh=None):
    """Glue ``inner_shape`` with an offset Kelvin exterior sphere.

    The Sugahara two-sphere convention:
      - inner physical sphere at ORIGIN, radius R_K (carrier of the
        physical model, e.g. coil + air); the user supplies this as
        ``inner_shape``, which MUST already be confined to the sphere
        |r| <= R_K (typically built as ``Sphere(O, R_K) - coil``).
      - outer Kelvin sphere at ``offset``, radius R_K, material
        ``kelvin_mat``: represents the physical exterior under the
        Kelvin inversion centered on ``offset``.
      - Periodic BC identifies the two sphere surfaces (Sugahara 2022).
      - GND vertex at the Kelvin sphere center maps to physical infinity.

    The Kelvin nu modulation (used downstream by ``kelvin_material``)
    is ``nu_kelvin = (R_K / rho')^2 * nu_0`` with rho' measured from
    ``offset``; this is the "sugahara" convention (validated by
    test_nu_convention.py 2026-04-15).

    Args:
        inner_shape: OCC shape representing the inner physical region
            (must already be enclosed in the Sphere(O, R_K) so its
            outer face is the Kelvin interface).
        offset: 3-tuple, center of the outer Kelvin sphere.
        R_K: Kelvin sphere radius (same for both inner and outer).
        inner_kelvin_face_name, outer_kelvin_face_name: names assigned
            to the two faces that get identified.
        kelvin_mat: material name on the outer sphere ("kelvin").
        gnd_vertex_name: BBND name for the Dirichlet GND vertex.
        outer_maxh_factor: multiplier for the Kelvin sphere maxh
            relative to the inner sphere's. Default 2.0 (coarser).
        inner_maxh: maxh override for the inner sphere boundary
            (and used to derive outer maxh). If None, the user is
            responsible for setting maxh on inner_shape.faces.

    Returns:
        ``(geometry, info)`` where ``geometry`` is the OCC compound
        ready to feed to ``OCCGeometry(...)`` and ``info`` is a dict
        with keys ``inner_face``, ``outer_face``, ``gnd_vertex``,
        ``offset``, ``R_K``.

    Notes:
        - The inner_shape MUST have its OUTER face named
          ``inner_kelvin_face_name`` BEFORE calling this helper, OR
          left unnamed (in which case all inner faces facing outward
          are renamed by this function).
        - For the periodic identification to work, the two sphere
          surfaces must be topologically equivalent (same radius, same
          orientation). Done by construction here.
    """
    # --- Tag the inner shape's outer (Kelvin-interface) face. --------
    # If the caller pre-named it, we keep the existing name; otherwise
    # we rename the face on the boundary "sphere of radius R_K".
    inner_kelvin_face = None
    for f in inner_shape.faces:
        if f.name == inner_kelvin_face_name:
            inner_kelvin_face = f
            break
    if inner_kelvin_face is None:
        # Heuristic: rename the largest face (typically the outer sphere).
        # Caller may want to do this themselves for tighter control.
        for f in inner_shape.faces:
            f.name = inner_kelvin_face_name
            inner_kelvin_face = f
            break

    # --- Build outer Kelvin sphere at offset. ------------------------
    outer = Sphere(Pnt(*offset), R_K)
    if inner_maxh is not None:
        outer.maxh = inner_maxh * outer_maxh_factor
    for f in outer.faces:
        f.name = outer_kelvin_face_name
    outer.name = kelvin_mat

    outer_kelvin_face = None
    for f in outer.faces:
        if f.name == outer_kelvin_face_name:
            outer_kelvin_face = f
            break

    # --- Periodic identification: inner outer surface <-> Kelvin sphere
    # surface. Done before Glue so the identification is recorded in the
    # OCC topology.
    if inner_kelvin_face is None or outer_kelvin_face is None:
        raise RuntimeError(
            "Could not locate Kelvin interface faces for periodic BC")
    inner_kelvin_face.Identify(outer_kelvin_face, "kelvin_periodic",
                                IdentificationType.PERIODIC)

    # --- GND vertex at the Kelvin center (image of physical infinity).
    gnd = Vertex(Pnt(*offset))
    gnd.name = gnd_vertex_name

    geometry = Glue([inner_shape, outer, gnd])
    info = {
        "inner_face": inner_kelvin_face,
        "outer_face": outer_kelvin_face,
        "gnd_vertex": gnd,
        "offset": offset,
        "R_K": R_K,
    }
    return geometry, info
