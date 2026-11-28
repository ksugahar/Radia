"""Kelvin open-boundary helpers: Cubit 3D, OCC 3D and OCC 2D axisymmetric.

The module is standard-library only at import time; each entry
imports ``cubit`` or ``netgen.occ`` inside its body, so the Cubit menu loads
this file by path inside Cubit's embedded Python.

* ``add_kelvin_cubit`` / ``auto_add_kelvin_from_current_model`` /
  ``apply_options``: Cubit Python, after the journal is played and before
  ``export vol``; called by the external curved-export worker and by the
  Vol Export menu (:mod:`radia_mcp.cubit.gui.cubit_export_menu`).
* ``add_kelvin_occ`` / ``add_kelvin_2d_axisym``: NGSolve/OCC geometry, before
  meshing, in a normal Python with ``netgen``.

Cubit 3D semantics: create an exterior sphere of the air radius at an
offset, handle the symmetry cuts (mesh-seam or domain-reduction mode), copy
the air outer-surface mesh 1:1, tet-mesh the Kelvin volume, assign the
``kelvin`` block, sidesets and the GND nodeset.  The plugin labels the shared
air|kelvin faces ``kelvin_int`` and the outer Kelvin faces ``kelvin_ext``, and
writes their ``kelvin`` periodic group when the export command carries
``export_arguments(info)`` (the plugin refuses a Kelvin mesh without them).

Change from the retired plugin: a missing or empty air block, an invalid per-axis
symmetry value or an unusable air surface now raises instead of warning and
exporting without the requested open boundary.
"""

from __future__ import annotations

import math


# ====================================================================
# Shared logic
# ====================================================================

# Map symmetry plane name -> normal vector
_SYM_NORMALS = {"x": (1, 0, 0), "y": (0, 1, 0), "z": (0, 0, 1)}

# Valid symmetry-plane BC labels in reduction mode.
#   "bn=0"  -- B . n = 0 on the plane (flux parallel to plane).
#              Omega: natural (do nothing).  A: Dirichlet (A x n = 0).
#              Magnetic image sign: '+' (mirror-symmetric).
#   "ht=0"  -- H x n = 0 on the plane (flux perpendicular to plane).
#              Omega: Dirichlet (Omega=const).  A: natural.
#              Magnetic image sign: '-' (mirror-antisymmetric).
_VALID_BC = ("bn=0", "ht=0")


def _cubit_module():
    """Cubit's Python API; only the Cubit-route helpers need it.

    The OCC and axisymmetric helpers in this module work without Cubit, so
    the import stays inside the Cubit functions and fails with a clear message.
    """
    try:
        import cubit
    except ImportError as exc:
        raise RuntimeError(
            "this Kelvin helper drives Cubit and needs its Python API: run it "
            "inside Cubit or after cubit.init()") from exc
    return cubit


def sym_sideset_name(axis, bc):
    """Return the canonical sideset label for a symmetry plane.

    Args:
        axis: one of "x", "y", "z" (lowercase).
        bc: one of "bn=0", "ht=0".

    Returns:
        e.g. "sym_bn=0_x", "sym_ht=0_z".  This string is the single
        source of truth -- the Radia EM headless backends
        parse it with `parse_sym_label`.
    """
    axis = axis.lower()
    if axis not in ("x", "y", "z"):
        raise ValueError(f"axis must be x|y|z, got {axis!r}")
    if bc not in _VALID_BC:
        raise ValueError(
            f"bc must be one of {_VALID_BC}, got {bc!r}")
    return f"sym_{bc}_{axis}"


def parse_sym_label(label):
    """Inverse of `sym_sideset_name`.

    Returns:
        (axis, bc) tuple if label is a sym_*_* label, else None.

    Recognises both the canonical "sym_bn=0_x" / "sym_ht=0_z" names and
    the LEGACY "sym_tangential" / "sym_normal" labels (no axis).  For
    the legacy names we return ("", "bn=0") / ("", "ht=0") so callers
    can still apply the physically-correct BC without knowing which
    axis it lives on.  Legacy = pre-2026-04-25; new code should emit
    the canonical names.
    """
    if not isinstance(label, str):
        return None
    # Legacy back-compat.
    if label == "sym_tangential":
        return ("", "bn=0")   # tangential B = B parallel to plane = Bn=0
    if label == "sym_normal":
        return ("", "ht=0")   # normal B = B perpendicular = Ht=0 on plane
    # Canonical: sym_<bc>_<axis>
    for bc in _VALID_BC:
        prefix = f"sym_{bc}_"
        if label.startswith(prefix):
            axis = label[len(prefix):].lower()
            if axis in ("x", "y", "z"):
                return (axis, bc)
    return None


def auto_offset_direction(symmetry):
    """Return the offset axis that is compatible with all symmetry planes.

    Rule: the offset direction must lie in the **intersection** of all
    symmetry planes.  That means it must be perpendicular to every
    symmetry-plane normal -- i.e. along an axis whose plane is NOT a
    symmetry plane.

    Examples:
        symmetry=[]        -> "x"   (arbitrary default)
        symmetry=["z"]     -> "x"   (offset in x-y plane)
        symmetry=["x","z"] -> "y"   (only free axis)
        symmetry=["x","y","z"] -> ValueError (1/8 impossible)

    Returns:
        One of "x", "y", "z".
    """
    blocked = set(s.lower() for s in symmetry)
    free = [d for d in ("x", "y", "z") if d not in blocked]
    if not free:
        # 1/8 symmetry: all axes are symmetry planes.
        # Offset along diagonal (1,1,1).  The Kelvin sphere center
        # is at (d/sqrt3, d/sqrt3, d/sqrt3); all three symmetry
        # planes cut through it asymmetrically.
        return "diag"
    return free[0]


def compute_offset_vector(R, symmetry, offset_dir=None, offset_dist=None):
    """Compute the (ox, oy, oz) offset vector.

    Args:
        R: Kelvin sphere radius [m].
        symmetry: list of symmetry plane names, e.g. ["z"] or ["x","z"].
        offset_dir: explicit axis ("x"/"y"/"z"), or None for auto.
        offset_dist: explicit distance [m], or None for 3*R default.

    Returns:
        (ox, oy, oz) tuple [m].
    """
    if offset_dir is None:
        offset_dir = auto_offset_direction(symmetry)
    offset_dir = offset_dir.lower()

    if offset_dir != "diag" and offset_dir in set(s.lower() for s in symmetry):
        raise ValueError(
            "Offset direction '%s' conflicts with symmetry plane %s=0.  "
            "The offset must be perpendicular to all symmetry normals."
            % (offset_dir, offset_dir))

    if offset_dist is None:
        offset_dist = 3.0 * R  # >= 2R required, 3R is comfortable
    if offset_dist < 2.0 * R:
        raise ValueError(
            "Offset distance %.4f must be >= 2*R = %.4f to avoid overlap."
            % (offset_dist, 2.0 * R))

    if offset_dir == "diag":
        # 1/8 symmetry: diagonal offset (1,1,1)/sqrt(3) * dist
        s3 = 3.0 ** 0.5
        d = offset_dist / s3
        return (d, d, d)

    idx = {"x": 0, "y": 1, "z": 2}[offset_dir]
    vec = [0.0, 0.0, 0.0]
    vec[idx] = offset_dist
    return tuple(vec)


def _sweep_params(symmetry):
    """Determine arc sweep angle and seed direction from symmetry planes.

    The Kelvin sphere is created by sweeping an arc around the z-axis.
    Symmetry planes reduce the sweep angle.

    Returns:
        (sweep_angle_deg, seed_direction, use_quarter_arc)
        seed_direction is (sx, sy, sz) on the equator at distance R=1.
    """
    sym = set(s.lower() for s in symmetry)

    if "x" in sym and "y" in sym:
        # 1/4 about x=0 and y=0 -> 90 deg sweep
        return 90, (1, 0, 0), "z" in sym
    elif "y" in sym:
        # 1/2 about y=0 -> 180 deg, seed on +x
        return 180, (1, 0, 0), "z" in sym
    elif "x" in sym:
        # 1/2 about x=0 -> 180 deg, seed on +y
        return 180, (0, 1, 0), "z" in sym
    else:
        # full sphere -> 360 deg
        return 360, (1, 0, 0), "z" in sym


# ====================================================================
# Cubit path (Layer 2 -- Cubit Python 3.10)
# ====================================================================

# ====================================================================
# Which faces are the Kelvin interface (non-reduction path)
# ====================================================================

_SPHERE_TOLERANCE = 1.0e-6      # relative to R
_AXIS_PROBES = ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1))


def _faces_on_sphere(cubit, volume, centre, R):
    """Faces of ``volume`` lying on the sphere of radius ``R`` about ``centre``."""
    faces = []
    for surface in cubit.get_relatives("volume", volume, "surface"):
        if cubit.get_surface_type(surface) != "sphere surface":
            continue
        near = cubit.surface(surface).closest_point_trimmed(list(centre))
        if abs(math.dist(near, centre) - R) <= _SPHERE_TOLERANCE * R:
            faces.append(surface)
    if len(faces) > 1:
        raise RuntimeError(
            "Volume %d has %d faces on the sphere of radius %g about %s; "
            "expected at most one" % (volume, len(faces), R, tuple(centre)))
    return faces


def _air_sphere_faces(cubit, air_vols, R, air_block):
    """Faces of the air domain on the physical sphere, which must sit at the origin.

    Only these belong to kelvin_int.  The largest face of every air volume
    used to be taken instead, so the faces of gap slabs inside the air went
    in too.  The rest of add_kelvin_cubit -- the exterior sphere placed at an
    offset from the origin, the zplane webcut, the periodic identification --
    assumes an air sphere centred at the origin, so a shifted air domain is
    rejected here rather than half supported.
    """
    boxes = [cubit.volume(v).bounding_box() for v in air_vols]
    centre = [0.5 * (min(b[i] for b in boxes) + max(b[i + 3] for b in boxes))
              for i in range(3)]
    if math.dist(centre, (0.0, 0.0, 0.0)) > _SPHERE_TOLERANCE * R:
        raise NotImplementedError(
            "add_kelvin_cubit supports an air sphere centred at the origin "
            "only; block '%s' is centred at %s" % (air_block, tuple(centre)))
    faces = [face for volume in air_vols
             for face in _faces_on_sphere(cubit, volume, (0.0, 0.0, 0.0), R)]
    if not faces:
        raise RuntimeError(
            "No face of block '%s' lies on the sphere of radius %g about the "
            "origin; R must match the air sphere" % (air_block, R))
    return faces


def _sphere_signature(cubit, face, centre, R):
    """Which of the six points ``centre +/- R e_i`` lie on ``face``."""
    signature = []
    for axis in _AXIS_PROBES:
        point = [centre[i] + R * axis[i] for i in range(3)]
        near = cubit.surface(face).closest_point_trimmed(point)
        signature.append(math.dist(near, point) <= _SPHERE_TOLERANCE * R)
    return tuple(signature)


def _pair_sphere_faces(cubit, air_faces, kelvin_faces, offset, R):
    """Match every air sphere face to the Kelvin face that will carry its mesh.

    Returns the copy-mesh pairs ``(air, kelvin)``.
    Faces are matched by which axis points of their own sphere they contain,
    so a Kelvin face cut on the wrong plane does not match.  A face without a
    partner raises; the former ``min(len(air), len(kelvin))`` dropped it.
    """
    origin = (0.0, 0.0, 0.0)
    air = {face: _sphere_signature(cubit, face, origin, R) for face in air_faces}
    kelvin = {face: _sphere_signature(cubit, face, offset, R) for face in kelvin_faces}
    if len(air) != len(kelvin):
        raise RuntimeError("%d air sphere faces but %d Kelvin faces"
                           % (len(air), len(kelvin)))
    pairs, remaining = [], dict(kelvin)
    for face, signature in air.items():
        match = [k for k, s in remaining.items() if s == signature]
        if len(match) != 1:
            raise RuntimeError("Air sphere face %d has %d matching Kelvin faces"
                               % (face, len(match)))
        pairs.append((face, match[0]))
        del remaining[match[0]]
    return pairs


def _surface_tri_nodes(cubit, surface):
    return [tuple(cubit.get_connectivity("tri", t))
            for t in cubit.get_surface_tris(surface)]


def _check_translated_copy(cubit, src, dst, offset, R):
    """Raise unless surface ``dst`` carries ``src``'s mesh moved by ``offset``.

    Every node and triangle of the air face must map one to one onto the
    Kelvin face, each node within 1e-12 R of its translated position.
    Returns the largest node deviation.
    """
    tol = 1e-12 * R
    cell = 1e-9 * R
    src_tris = _surface_tri_nodes(cubit, src)
    dst_tris = _surface_tri_nodes(cubit, dst)
    src_nodes = sorted({n for t in src_tris for n in t})
    dst_nodes = sorted({n for t in dst_tris for n in t})
    if len(src_tris) != len(dst_tris) or len(src_nodes) != len(dst_nodes):
        raise RuntimeError(
            "Kelvin surface %d has %d nodes / %d triangles; air surface %d has "
            "%d / %d" % (dst, len(dst_nodes), len(dst_tris), src,
                         len(src_nodes), len(src_tris)))
    grid = {}
    for n in dst_nodes:
        p = cubit.get_nodal_coordinates(n)
        grid.setdefault(tuple(int(round(c / cell)) for c in p), []).append((n, p))
    mapping, worst = {}, 0.0
    for n in src_nodes:
        q = [c + o for c, o in zip(cubit.get_nodal_coordinates(n), offset)]
        key = [int(round(c / cell)) for c in q]
        best = None
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    cand = grid.get((key[0] + dx, key[1] + dy, key[2] + dz), ())
                    for m, p in cand:
                        d = max(abs(a - b) for a, b in zip(p, q))
                        if best is None or d < best[0]:
                            best = (d, m)
        if best is None or best[0] > tol:
            raise RuntimeError(
                "Air node %d on surface %d has no Kelvin node within %.3g m "
                "of its translated position on surface %d" % (n, src, tol, dst))
        worst = max(worst, best[0])
        mapping[n] = best[1]
    if len(set(mapping.values())) != len(dst_nodes):
        raise RuntimeError("Kelvin surface %d node map is not one to one" % dst)
    dst_set = {frozenset(t) for t in dst_tris}
    for t in src_tris:
        if frozenset(mapping[n] for n in t) not in dst_set:
            raise RuntimeError("Air triangle %s on surface %d has no translated "
                               "Kelvin triangle on surface %d" % (t, src, dst))
    return worst


def add_kelvin_cubit(R, air_block="air", symmetry=None, reduction=None,
                     offset_dir=None, offset_dist=None, mesh_size=None,
                     kelvin_block="kelvin", gnd_nodeset=100):
    """Add Kelvin open-boundary sphere to the current Cubit model.

    Call **after** meshing the physical domain (coil, workpiece, air).
    Creates a fresh exterior sphere at the offset position, webcutted
    for symmetry, with 1:1 mesh copy from the air sphere outer surface.

    Strategy:
      1. Find the existing air volumes and their outer spherical surfaces
      2. Create a fresh ACIS sphere at the offset position
      3. Webcut the exterior sphere for z-symmetry (equator curves)
      4. Imprint + merge the exterior sphere halves (NOT with air!)
      5. Copy mesh from air outer surfaces to exterior sphere surfaces
      6. Tet-mesh the exterior sphere volumes
      7. Assign blocks, sidesets, GND nodeset

    Args:
        R: Kelvin sphere radius [m] (must match air sphere radius).
        air_block: Name of the existing air block (to find outer surface).
        symmetry: List of webcut-seam axes (mesh-seam mode), e.g.
            ["z"] = full sphere with z=0 mesh seam kept on both sides.
            Mutually exclusive with `reduction`.  Default None.
        reduction: Dict {axis: bc} for domain-reduction mode, e.g.
            {"x": "ht=0", "z": "bn=0"}.  Air block is expected to
            ALREADY be reduced; the Kelvin sphere is cut on the same
            planes and each cut face gets a `sym_<bc>_<axis>` sideset.
            Supports 1/2 (one axis) and 1/4 (two axes).  1/8 is a
            geometric blocker -- raises NotImplementedError.
            Mutually exclusive with `symmetry`.
        offset_dir: Explicit offset axis ("x"/"y"/"z") or None (auto).
        offset_dist: Explicit offset distance [m] or None (3*R).
        mesh_size: Kelvin tet size [m] or None (auto from copy mesh).
        kelvin_block: Name for the Kelvin material block.
        gnd_nodeset: Nodeset ID for GND vertex.

    Returns:
        dict with keys:
            'R', 'center', 'symmetry', 'reduction', 'offset_dir',
            'air_vols', 'outer_vols', 'sym_sidesets'
    """
    cubit = _cubit_module()

    if symmetry is not None and reduction is not None:
        raise ValueError(
            "symmetry= (webcut-seam mode) and reduction= (domain-"
            "reduction mode) are mutually exclusive.  Pick one.")

    # Domain-reduction path delegates to the dedicated implementation.
    if reduction is not None:
        return _add_kelvin_cubit_reduction(
            R=R, air_block=air_block, reduction=reduction,
            offset_dir=offset_dir, offset_dist=offset_dist,
            mesh_size=mesh_size, kelvin_block=kelvin_block,
            gnd_nodeset=gnd_nodeset)

    if symmetry is None:
        symmetry = []

    ox, oy, oz = compute_offset_vector(R, symmetry, offset_dir, offset_dist)
    actual_dir = auto_offset_direction(symmetry) if offset_dir is None else offset_dir

    # ---- 1. Find existing air volumes ----
    air_bid = None
    for bid in cubit.parse_cubit_list("block", "all"):
        try:
            n = cubit.get_block_name(bid) or ""
        except Exception:
            n = ""
        if n.lower() == air_block.lower():
            air_bid = bid
            break

    if air_bid is None:
        raise RuntimeError(
            "Block '%s' not found.  Create and mesh the air sphere "
            "before calling add_kelvin_cubit()." % air_block)

    air_vols = list(cubit.parse_cubit_list("volume", "in block %d" % air_bid))
    if not air_vols:
        raise RuntimeError("Block '%s' has no volumes." % air_block)

    # ---- 1b. Faces of the air domain on the physical sphere ----
    # Selected before any geometry is created, so a rejected model is left
    # untouched.
    air_outer_surfs = _air_sphere_faces(cubit, air_vols, R, air_block)

    # ---- 2. Create exterior sphere (clean ACIS sphere) ----
    cubit.cmd("create sphere radius %g" % R)
    kelvin_sphere = cubit.get_last_id("volume")
    cubit.cmd("move volume %d x %g y %g z %g" % (kelvin_sphere, ox, oy, oz))

    # ---- 3. Webcut for symmetry ----
    sym = set(s.lower() for s in symmetry)
    kelvin_vols = [kelvin_sphere]

    if "z" in sym:
        # Both hemispheres are kept and each receives a copy of the matching
        # air hemisphere mesh.  A meshed-and-reflected lower hemisphere is
        # not used: it matches kelvin_int only when the air mesh itself is
        # mirror-symmetric, which an independently meshed air domain is not.
        before_cut = set(cubit.parse_cubit_list("volume", "all"))
        cubit.cmd("webcut volume %d with plane zplane" % kelvin_sphere)
        halves = sorted({kelvin_sphere} | (
            set(cubit.parse_cubit_list("volume", "all")) - before_cut))
        if len(halves) != 2:
            raise RuntimeError(
                "Kelvin z-webcut produced %d volumes instead of two: %s"
                % (len(halves), halves))
        kelvin_vols = halves

    # Imprint + merge kelvin halves (NOT with air volumes!)
    if len(kelvin_vols) > 1:
        k_str = " ".join(str(v) for v in kelvin_vols)
        cubit.cmd("imprint volume %s" % k_str)
        cubit.cmd("merge volume %s" % k_str)

    for vid in kelvin_vols:
        cubit.cmd('volume %d rename "kelvin_ext"' % vid)

    # ---- 4/5. Sphere faces of the Kelvin exterior ----
    # (The air side was selected in step 1b.)
    kelvin_outer_surfs = [face for vid in kelvin_vols
                          for face in _faces_on_sphere(cubit, vid, (ox, oy, oz), R)]

    # ---- 6. Copy mesh: air outer surface -> kelvin surface ----
    # Every air sphere face is matched to the Kelvin face that will carry its
    # mesh; after the volume mesh each pair is checked to be an exact
    # translated copy.
    pairs = _pair_sphere_faces(
        cubit, air_outer_surfs, kelvin_outer_surfs, (ox, oy, oz), R)
    for src, dst in pairs:

        src_curves = list(cubit.get_relatives("surface", src, "curve"))
        dst_curves = list(cubit.get_relatives("surface", dst, "curve"))
        if not src_curves or not dst_curves:
            raise RuntimeError(
                "Surface %d or %d has no curves -- cannot copy mesh"
                % (src, dst))

        src_c = max(src_curves, key=lambda c: cubit.curve(c).length())
        dst_c = max(dst_curves, key=lambda c: cubit.curve(c).length())
        src_v = cubit.get_relatives("curve", src_c, "vertex")[0]
        dst_v = cubit.get_relatives("curve", dst_c, "vertex")[0]

        cubit.cmd(
            "copy mesh surface %d onto surface %d "
            "source curve %d source vertex %d "
            "target curve %d target vertex %d"
            % (src, dst, src_c, src_v, dst_c, dst_v))
        if not cubit.is_meshed("surface", dst):
            raise RuntimeError("Cubit could not copy the mesh of air surface "
                               "%d onto Kelvin surface %d" % (src, dst))

    # ---- 7. Mesh exterior kelvin volumes ----
    # When mesh_size is None (default), let the copied surface mesh
    # constrain the interior — this matches adjacent air-sphere density.
    # When mesh_size is given, impose it so the Kelvin exterior can be
    # coarser than the physical domain (typical case, since Kelvin is
    # just an open-boundary trick and doesn't need fine physics).
    k_str = " ".join(str(v) for v in kelvin_vols)
    for vid in kelvin_vols:
        cubit.cmd("volume %d scheme tetmesh" % vid)
        if mesh_size is not None:
            cubit.cmd("volume %d size %g" % (vid, float(mesh_size)))
    cubit.cmd("mesh volume %s" % k_str)

    for src, dst in pairs:
        _check_translated_copy(cubit, src, dst, (ox, oy, oz), R)

    # ---- 8. Block assignment ----
    existing_blocks = set(cubit.parse_cubit_list("block", "all"))
    nb = (max(existing_blocks) + 1) if existing_blocks else 1

    existing_block_names = {}
    for bid in cubit.parse_cubit_list("block", "all"):
        try:
            n = cubit.get_block_name(bid) or ""
        except Exception:
            n = ""
        existing_block_names[n.lower()] = bid

    if kelvin_block.lower() not in existing_block_names:
        cubit.cmd("block %d add volume %s" % (nb, k_str))
        cubit.cmd('block %d name "%s"' % (nb, kelvin_block))
        nb += 1

    # ---- 9. Sideset assignment ----
    existing_ss = set(cubit.parse_cubit_list("sideset", "all"))
    ns = (max(existing_ss) + 1) if existing_ss else 1

    existing_ss_names = set()
    for sid in cubit.parse_cubit_list("sideset", "all"):
        try:
            n = cubit.get_sideset_name(sid) or ""
        except Exception:
            n = ""
        existing_ss_names.add(n.lower())

    if "kelvin_int" not in existing_ss_names:
        cubit.cmd("sideset %d add surface %s"
                  % (ns, " ".join(str(s) for s in air_outer_surfs)))
        cubit.cmd('sideset %d name "kelvin_int"' % ns)
        ns += 1

    if "kelvin_ext" not in existing_ss_names:
        # Add surfaces of kelvin volumes; skip merged-away surfaces.
        # Use parse_cubit_list which only returns live (non-merged) IDs.
        all_kelvin_surfs = set()
        live_surfs = set(cubit.parse_cubit_list("surface", "all"))
        for vid in kelvin_vols:
            for s in cubit.get_relatives("volume", vid, "surface"):
                # Port change: skip the merged seam face shared by two
                # Kelvin volumes; the export plugin rejects an explicit
                # label on a same-material internal surface.
                if s in live_surfs and len(cubit.get_relatives("surface", s, "volume")) == 1:
                    all_kelvin_surfs.add(s)
        if all_kelvin_surfs:
            cubit.cmd("sideset %d add surface %s"
                      % (ns, " ".join(str(s)
                                      for s in sorted(all_kelvin_surfs))))
            cubit.cmd('sideset %d name "kelvin_ext"' % ns)

    # ---- 10. GND nodeset at Kelvin center ----
    cubit.cmd("create vertex %g %g %g" % (ox, oy, oz))
    gnd_vid = cubit.get_last_id("vertex")
    cubit.cmd("nodeset %d add vertex %d" % (gnd_nodeset, gnd_vid))
    cubit.cmd('nodeset %d name "GND"' % gnd_nodeset)

    # ---- 11. Visual feedback ----
    for vid in kelvin_vols:
        cubit.cmd("volume %d visibility off" % vid)

    n_sym = len(symmetry)
    frac = "1/%d" % (2 ** n_sym) if n_sym > 0 else "full"

    print("")
    print("=== add_kelvin_cubit ===")
    print("  R = %g m" % R)
    print("  offset = (%g, %g, %g) m" % (ox, oy, oz))
    print("  symmetry = %s (%s)" % (symmetry, frac))
    print("  air volumes = %s (kelvin_int surfaces = %s)"
          % (air_vols, air_outer_surfs))
    print("  kelvin volumes = %s" % kelvin_vols)
    print("  GND nodeset %d at (%g, %g, %g)" % (gnd_nodeset, ox, oy, oz))
    print("  block '%s', sidesets 'kelvin_int'/'kelvin_ext'" % kelvin_block)

    return {
        "R": R,
        "center": (ox, oy, oz),
        "symmetry": list(symmetry),
        "reduction": None,
        "offset_dir": actual_dir,
        "air_vols": air_vols,
        "outer_vols": kelvin_vols,
        "sym_sidesets": {},
    }


# ====================================================================
# Cubit path -- reduction mode (2026-04-25)
# ====================================================================

def _add_kelvin_cubit_reduction(R, air_block, reduction,
                                 offset_dir, offset_dist, mesh_size,
                                 kelvin_block, gnd_nodeset):
    """Reduction-mode implementation (1/2 or 1/4 domain).

    Contract:
      - The air block already contains only the reduced domain
        (vertices all >= 0 on each reduction axis).
      - The Kelvin sphere is created at the offset location, then
        webcut on the same symmetry planes so it also occupies the
        reduced region.
      - Each symmetry plane gets a `sym_<bc>_<axis>` sideset that
        includes the air's cut face AND the Kelvin's cut face (both
        sides of the plane - solver sees a single boundary).
      - 1/8 (all three planes) raises NotImplementedError: the Kelvin
        offset lands at (d/sqrt3, d/sqrt3, d/sqrt3) which does NOT
        lie on the x=0/y=0/z=0 planes, so the Kelvin cut faces
        cannot coincide with the air cut faces.  Requires a different
        Kelvin mapping (e.g. 7 reflected copies) -- deferred.
    """
    cubit = _cubit_module()

    # Validate inputs.
    if not isinstance(reduction, dict):
        raise TypeError(
            f"reduction must be a dict {{axis: bc}}, got {type(reduction).__name__}")
    axes = []
    for axis, bc in reduction.items():
        axis = axis.lower()
        if axis not in ("x", "y", "z"):
            raise ValueError(f"reduction axis must be x|y|z, got {axis!r}")
        if bc not in _VALID_BC:
            raise ValueError(
                f"reduction[{axis!r}] must be one of {_VALID_BC}, got {bc!r}")
        axes.append(axis)
    if len(axes) == 0:
        raise ValueError(
            "reduction dict is empty.  For full domain, pass "
            "symmetry=None (no reduction).")
    if len(axes) > 3:
        raise ValueError(
            f"reduction has {len(axes)} axes; only x, y, z are valid.")
    is_eighth = (len(axes) == 3)

    # Physically impossible: 3 mutually perpendicular planes ALL set to
    # bn=0 means B parallel to all 3 planes, which forces B = 0.  Reject
    # with a clear message so the user knows to pick at least one ht=0.
    if is_eighth and all(reduction[a] == "bn=0" for a in axes):
        raise ValueError(
            "1/8 reduction with all 3 planes set to 'bn=0' is "
            "physically impossible: B parallel to three mutually "
            "perpendicular planes forces B = 0 everywhere.  "
            "At least one axis must be 'ht=0' (antimirror).  "
            "Common ELF-style choice: {'x': 'ht=0', 'y': 'ht=0', "
            "'z': 'bn=0'} (Radia '-x-y+z').")

    # Offset direction.
    #   1/2 / 1/4: must lie along a FREE axis (perpendicular to all
    #              reduction-axis normals); the reduction sym planes
    #              are sym BCs and the Kelvin sits "outside" them.
    #   1/8     : no free axis exists, so the offset MUST be along one
    #              of the reduction axes.  Two of the three Kelvin cut
    #              planes (those of the non-offset reduction axes) end
    #              up coinciding with air sym planes; the third (the
    #              offset-axis cut, through the Kelvin centre) is the
    #              "infinity plane" and gets a fresh `kelvin_far`
    #              sideset for a Dirichlet anchor.
    if offset_dir is None:
        free = [d for d in ("x", "y", "z") if d not in axes]
        if free:
            offset_dir = free[0]
        else:
            # 1/8 default: offset along the first reduction axis.
            offset_dir = axes[0]
    offset_dir = offset_dir.lower()
    if offset_dir not in ("x", "y", "z"):
        raise ValueError(
            f"offset_dir must be x|y|z, got {offset_dir!r}")
    if (not is_eighth) and offset_dir in axes:
        raise ValueError(
            f"offset_dir={offset_dir!r} conflicts with reduction on "
            f"the same axis (Kelvin sphere would poke through the "
            f"symmetry plane).  Pick a free axis from "
            f"{[d for d in ('x','y','z') if d not in axes]}.")

    if offset_dist is None:
        offset_dist = 3.0 * R
    if offset_dist < 2.0 * R:
        raise ValueError(
            f"offset_dist {offset_dist:.4f} must be >= 2*R = {2*R:.4f}.")

    idx = {"x": 0, "y": 1, "z": 2}[offset_dir]
    center = [0.0, 0.0, 0.0]
    center[idx] = offset_dist
    ox, oy, oz = center

    # ---- 1. Locate the air block ----
    air_bid = None
    for bid in cubit.parse_cubit_list("block", "all"):
        n = cubit.get_exodus_entity_name("block", bid) or ""
        if n.lower() == air_block.lower():
            air_bid = bid
            break
    if air_bid is None:
        raise RuntimeError(
            f"Block {air_block!r} not found.  Create and mesh the "
            f"reduced air domain before calling add_kelvin_cubit().")
    air_vols = list(cubit.parse_cubit_list("volume", f"in block {air_bid}"))
    if not air_vols:
        raise RuntimeError(f"Block {air_block!r} has no volumes.")

    # ---- 2. Sanity-check reduction: all air vertices on the reduced side ----
    tol = 1e-6
    for axis in axes:
        ai = {"x": 0, "y": 1, "z": 2}[axis]
        min_coord = float("inf")
        for vid in air_vols:
            for v in cubit.get_relatives("volume", vid, "vertex"):
                c = cubit.vertex(v).coordinates()[ai]
                if c < min_coord:
                    min_coord = c
        if min_coord < -tol * max(R, 1.0):
            raise RuntimeError(
                f"reduction[{axis!r}] requested, but air block has "
                f"vertices with {axis}={min_coord:.4e} (expected "
                f">= 0).  The air domain must be reduced to the "
                f"positive side of each reduction plane before "
                f"calling add_kelvin_cubit(reduction=...).")

    # ---- 3. Create exterior sphere, then cut on each reduction plane ----
    cubit.cmd(f"create sphere radius {R:g}")
    kelvin_sphere = cubit.get_last_id("volume")
    cubit.cmd(f"move volume {kelvin_sphere} x {ox:g} y {oy:g} z {oz:g}")

    # Webcut Kelvin at three (or fewer) planes.  Each cut plane passes
    # through the corresponding Kelvin-centre coordinate:
    #   - reduction axis NOT == offset_dir  -> coord = 0  (sym plane,
    #     coincides with the air's sym plane, sym_<bc>_<axis> sideset)
    #   - reduction axis == offset_dir       -> coord = offset_dist
    #     (1/8 only: "infinity plane" through Kelvin centre, kelvin_far
    #     sideset)
    # Always keep the half with centroid > cut_at; air points (in the
    # positive octant by convention) translate to Kelvin points which
    # are also > cut_at on every cut axis (sym cuts at 0 keep > 0;
    # offset cut at offset_dist keeps > offset_dist, since translated
    # air y goes from [0, R] to [offset_dist, offset_dist + R]).
    kelvin_vol = kelvin_sphere
    far_axis = offset_dir if is_eighth else None
    for axis in axes:
        ai = {"x": 0, "y": 1, "z": 2}[axis]
        cut_at = offset_dist if axis == far_axis else 0.0
        # Cubit's `xplane offset <X>` takes a signed offset along the axis.
        plane_cmd = f"{axis}plane offset {cut_at:g}" if cut_at != 0.0 \
                    else f"{axis}plane"
        before = set(cubit.parse_cubit_list("volume", "all"))
        cubit.cmd(f"webcut volume {kelvin_vol} with plane {plane_cmd}")
        after = set(cubit.parse_cubit_list("volume", "all"))
        new_vols = list(after - before)
        keep, drop = None, None
        candidates = [kelvin_vol] + new_vols
        for vid in candidates:
            try:
                c = cubit.volume(vid).centroid()[ai]
            except Exception:
                continue
            if c > cut_at + tol * max(R, 1.0):
                keep = vid
            elif c < cut_at - tol * max(R, 1.0):
                drop = vid
        if keep is None or drop is None:
            raise RuntimeError(
                f"webcut on {axis}={cut_at:g} did not produce a clear "
                f"keep/drop pair (candidates={candidates}, "
                f"centroids={[cubit.volume(v).centroid()[ai] for v in candidates]}).")
        cubit.cmd(f"delete volume {drop}")
        kelvin_vol = keep

    cubit.cmd(f'volume {kelvin_vol} rename "kelvin_ext"')

    # ---- 4. Find curved outer face of the reduced Kelvin volume ----
    # After reduction, the Kelvin volume has 1+len(axes) outer faces:
    # the curved sphere cap (we want this one) and the flat cut faces
    # (which will become sym_* sidesets).  The curved one is the one
    # whose bounding box spans all three coordinates when measured
    # relative to the Kelvin center; the flat cut faces lie in one of
    # the x=0 / y=0 / z=0 planes.
    # Each cut produces a flat face on a known plane:
    #   sym cuts (axis != far_axis): plane at coord = 0
    #   far cut (axis == far_axis):  plane at coord = offset_dist
    # Detect by bounding-box: zero span on axis A AND coord matches the
    # expected plane for A.
    k_surfs = list(cubit.get_relatives("volume", kelvin_vol, "surface"))
    k_curved = None
    k_flat_by_axis = {}
    k_far = None             # 1/8 only: face on the offset axis at offset_dist
    for sid in k_surfs:
        bb = cubit.surface(sid).bounding_box()
        # bb returns (xmin, ymin, zmin, xmax, ymax, zmax, ...)
        flat_axis = None
        flat_role = None     # 'sym' or 'far'
        for ai, ax_name in enumerate(("x", "y", "z")):
            lo, hi = bb[ai], bb[3 + ai]
            if abs(hi - lo) > 1e-5 * R:
                continue   # not flat in this axis
            plane_coord = 0.5 * (lo + hi)
            if ax_name in axes and ax_name != far_axis \
                    and abs(plane_coord) < 1e-5 * max(R, 1.0):
                flat_axis = ax_name
                flat_role = 'sym'
                break
            if ax_name == far_axis \
                    and abs(plane_coord - offset_dist) < 1e-5 * max(R, 1.0):
                flat_axis = ax_name
                flat_role = 'far'
                break
        if flat_role == 'sym':
            k_flat_by_axis[flat_axis] = sid
        elif flat_role == 'far':
            k_far = sid
        else:
            # Curved surface is the largest-area non-flat face.
            if (k_curved is None
                or cubit.surface(sid).area()
                   > cubit.surface(k_curved).area()):
                k_curved = sid
    if k_curved is None:
        raise RuntimeError(
            f"Cannot find curved face on Kelvin volume {kelvin_vol}.  "
            f"Surfaces = {k_surfs}.")
    for axis in axes:
        if axis == far_axis:
            if k_far is None:
                raise RuntimeError(
                    f"1/8 reduction: Kelvin offset-axis cut face "
                    f"(at {axis}={offset_dist:g}) not found among "
                    f"surfaces {k_surfs}.")
        else:
            if axis not in k_flat_by_axis:
                raise RuntimeError(
                    f"Cannot find flat sym face on axis {axis!r} for "
                    f"the Kelvin volume (expected face with "
                    f"{axis}=0 plane).  Surfaces = {k_surfs}.")

    # ---- 5. Find air's curved outer face + flat cut faces ----
    a_curved = None
    a_flat_by_axis = {axis: [] for axis in axes}
    for vid in air_vols:
        for sid in cubit.get_relatives("volume", vid, "surface"):
            bb = cubit.surface(sid).bounding_box()
            flat_axis = None
            for ai, ax_name in enumerate(("x", "y", "z")):
                lo, hi = bb[ai], bb[3 + ai]
                if ax_name in axes and abs(hi - lo) < 1e-5 * R and abs(lo) < 1e-5 * R:
                    flat_axis = ax_name
                    break
            if flat_axis is not None:
                a_flat_by_axis[flat_axis].append(sid)
            else:
                # Largest-area non-flat surface is the curved outer.
                if (a_curved is None
                    or cubit.surface(sid).area()
                       > cubit.surface(a_curved).area()):
                    a_curved = sid
    if a_curved is None:
        raise RuntimeError(
            "Cannot find curved outer face on air block.")

    # ---- 6. Copy-mesh the Kelvin curved face from the air curved face ----
    a_curves = list(cubit.get_relatives("surface", a_curved, "curve"))
    k_curves = list(cubit.get_relatives("surface", k_curved, "curve"))
    if not a_curves or not k_curves:
        raise RuntimeError(
            f"Air curved surface {a_curved} has {len(a_curves)} curves, "
            f"Kelvin curved surface {k_curved} has {len(k_curves)} "
            f"curves -- cannot copy mesh.")
    # Anchor curve selection: pick the arc with the SMALLEST centroid
    # z-coordinate.  For 1/4 / 1/8 reductions in the +octant the
    # boundary arcs of a sphere cap are equal-length quarter-arcs;
    # `max(..., key=length)` would tie and Cubit's internal listing
    # order can pick DIFFERENT arcs on source vs target, producing a
    # ~120-deg rotational mis-projection that lands ~5% of vertex
    # pairs ~2 cm off (verified 2026-04-25).  The arc on the z=0 plane
    # has centroid z=0; arcs on y=0 / x=0(or kelvin_far) planes have
    # centroid z>0.  Picking min-z is therefore deterministic AND
    # translation-equivalent on source and target.
    def _pick_anchor(curves):
        # Tuple-key for stable sorting: (centroid z, centroid y, centroid x)
        def _key(c):
            cp = cubit.curve(c).center_point()
            return (cp[2], cp[1], cp[0])
        return min(curves, key=_key)
    a_c = _pick_anchor(a_curves)
    k_c = _pick_anchor(k_curves)

    # Vertex anchor: pick the endpoint with the smallest (z, y, x) tuple.
    def _pick_anchor_vertex(curve_id):
        verts = list(cubit.get_relatives("curve", curve_id, "vertex"))
        def _key(v):
            cp = cubit.vertex(v).coordinates()
            return (cp[2], cp[1], cp[0])
        return min(verts, key=_key)
    a_v = _pick_anchor_vertex(a_c)
    k_v = _pick_anchor_vertex(k_c)
    cubit.cmd(
        f"copy mesh surface {a_curved} onto surface {k_curved} "
        f"source curve {a_c} source vertex {a_v} "
        f"target curve {k_c} target vertex {k_v}")

    if not cubit.is_meshed("surface", k_curved):
        raise RuntimeError("Cubit could not copy the mesh of air surface "
                           "%d onto Kelvin surface %d" % (a_curved, k_curved))

    # ---- 7. Mesh the Kelvin volume ----
    cubit.cmd(f"volume {kelvin_vol} scheme tetmesh")
    if mesh_size is not None:
        cubit.cmd(f"volume {kelvin_vol} size {float(mesh_size):g}")
    cubit.cmd(f"mesh volume {kelvin_vol}")
    _check_translated_copy(cubit, a_curved, k_curved, (ox, oy, oz), R)

    # ---- 8. Block assignment (kelvin) ----
    existing_blocks = set(cubit.parse_cubit_list("block", "all"))
    nb = (max(existing_blocks) + 1) if existing_blocks else 1
    existing_block_names = {
        (cubit.get_exodus_entity_name("block", bid) or "").lower()
        for bid in cubit.parse_cubit_list("block", "all")}
    if kelvin_block.lower() not in existing_block_names:
        cubit.cmd(f"block {nb} add volume {kelvin_vol}")
        cubit.cmd(f'block {nb} name "{kelvin_block}"')

    # ---- 9. Sidesets: kelvin_int / kelvin_ext + sym_<bc>_<axis> ----
    existing_ss = set(cubit.parse_cubit_list("sideset", "all"))
    ns = (max(existing_ss) + 1) if existing_ss else 1
    existing_ss_names = {
        (cubit.get_exodus_entity_name("sideset", sid) or "").lower()
        for sid in cubit.parse_cubit_list("sideset", "all")}

    def _next_ss():
        nonlocal ns
        cur = ns
        ns += 1
        return cur

    if "kelvin_int" not in existing_ss_names:
        sid = _next_ss()
        cubit.cmd(f"sideset {sid} add surface {a_curved}")
        cubit.cmd(f'sideset {sid} name "kelvin_int"')

    if "kelvin_ext" not in existing_ss_names:
        sid = _next_ss()
        cubit.cmd(f"sideset {sid} add surface {k_curved}")
        cubit.cmd(f'sideset {sid} name "kelvin_ext"')

    # Per-axis sym sidesets.  Each includes the AIR's flat cut face on
    # the same plane.  For axes != far_axis, the Kelvin's cut face on
    # that plane (at coord = 0) is also added -- both sides share the
    # sideset so the solver sees one boundary.  For the far_axis (1/8
    # only), the Kelvin's cut on that axis is at offset_dist (not 0)
    # and gets its own `kelvin_far` sideset; the air's cut on far_axis
    # at coord = 0 still becomes a sym sideset (single-side, air only).
    sym_sidesets = {}
    for axis in axes:
        bc = reduction[axis]
        label = sym_sideset_name(axis, bc)
        if label.lower() in existing_ss_names:
            sym_sidesets[axis] = None
            continue
        faces = list(a_flat_by_axis.get(axis, []))
        if axis != far_axis:
            kf = k_flat_by_axis.get(axis)
            if kf is not None:
                faces.append(kf)
        if not faces:
            raise RuntimeError(
                f"No faces found for sym axis={axis!r}.")
        sid = _next_ss()
        cubit.cmd(f"sideset {sid} add surface {' '.join(str(f) for f in faces)}")
        cubit.cmd(f'sideset {sid} name "{label}"')
        sym_sidesets[axis] = sid

    # 1/8 only: the Kelvin offset-axis cut face becomes the "kelvin_far"
    # sideset.  This is the "infinity plane" that passes through the
    # Kelvin centre, perpendicular to the offset axis -- not a sym
    # plane of the air.  The electromagnet calc paths treat it as a
    # Dirichlet boundary regardless of formulation (the same way GND works).
    if k_far is not None and "kelvin_far" not in existing_ss_names:
        sid = _next_ss()
        cubit.cmd(f"sideset {sid} add surface {k_far}")
        cubit.cmd(f'sideset {sid} name "kelvin_far"')

    # ---- 10. GND nodeset at Kelvin center ----
    cubit.cmd(f"create vertex {ox:g} {oy:g} {oz:g}")
    gnd_vid = cubit.get_last_id("vertex")
    cubit.cmd(f"nodeset {gnd_nodeset} add vertex {gnd_vid}")
    cubit.cmd(f'nodeset {gnd_nodeset} name "GND"')

    cubit.cmd(f"volume {kelvin_vol} visibility off")

    frac = "1/%d" % (2 ** len(axes))
    print("")
    print("=== add_kelvin_cubit (reduction mode) ===")
    print(f"  R = {R:g} m,  offset = ({ox:g}, {oy:g}, {oz:g}) m  "
          f"(offset_dir = {offset_dir})")
    print(f"  reduction = {reduction} ({frac} model)")
    print(f"  air_vols = {air_vols}, kelvin_vol = {kelvin_vol}")
    print(f"  sideset 'kelvin_int' (air curved) / 'kelvin_ext' (kelvin curved)")
    for axis, sid in sym_sidesets.items():
        label = sym_sideset_name(axis, reduction[axis])
        if axis == far_axis:
            note = f"(air cut at {axis}=0; kelvin face is kelvin_far)"
        else:
            note = f"(air cut + kelvin cut at {axis}=0)"
        print(f"  sideset {sid} '{label}'  {note}")
    if k_far is not None:
        print(f"  sideset 'kelvin_far'  (kelvin cut at {far_axis}={offset_dist:g})")
    print(f"  GND nodeset {gnd_nodeset} at Kelvin center")

    return {
        "R": R,
        "center": (ox, oy, oz),
        "symmetry": [],
        "reduction": dict(reduction),
        "offset_dir": offset_dir,
        "air_vols": air_vols,
        "outer_vols": [kelvin_vol],
        "sym_sidesets": sym_sidesets,
        "has_kelvin_far": k_far is not None,
    }


def auto_add_kelvin_from_current_model(air_block="air",
                                        kelvin_block="kelvin",
                                        mesh_size=None,
                                        reduction=None):
    """Detect air sphere + symmetry, then call add_kelvin_cubit().

    Called by the export worker after the journal and before
    `export vol`, in the worker's Cubit Python session.

    Steps (matching the 2026-04-14 c60a6007 implementation):
      1. If a `<kelvin_block>` block already exists, skip (idempotent).
      2. Find the `<air_block>` block; raise when it is missing or empty.
      3. Across all volumes in that block, pick the surface with the
         largest area -> the outer sphere boundary.
      4. R = max vertex distance from the origin on that surface.
      5. Detect symmetry axes: for each of x/y/z, if the air block's
         vertices all have coord >= 0 AND at least one vertex sits on
         the axis plane (coord ~ 0), that axis is a symmetry plane.
      6. Call ``add_kelvin_cubit(R, symmetry=[...], mesh_size=...)``.

    ``mesh_size`` (float in meters, or None): tet edge length for the
    Kelvin exterior.  None lets add_kelvin_cubit inherit from the air
    outer surface via copy-mesh (the usual sensible default).  The
    Kelvin region can usually be coarser than the physical domain —
    pass an explicit value (e.g. 2-3x air surface size) to override.

    Returns the info dict from ``add_kelvin_cubit``, or ``None`` when the
    Kelvin block already exists.
    Once an air sphere is detected, construction failures raise and abort the
    export; silently replacing the requested open boundary with Dirichlet
    truncation is forbidden.

    Args:
        air_block: block name holding the air volumes.
        kelvin_block: block name to assign to the Kelvin volume.
        mesh_size: override Kelvin tet size (m), or None to inherit.
        reduction: dict {axis: bc} forwarded verbatim to
            add_kelvin_cubit(reduction=...).  When provided, the
            function skips the auto-detection of mesh-seam symmetry
            and takes the reduction path directly.  E.g.
            reduction={"x": "ht=0", "z": "bn=0"} for a 1/4 xz model.
    """
    import math as _m
    cubit = _cubit_module()

    # --- Step 1: idempotent skip ---
    for bid in cubit.get_block_id_list():
        bn = cubit.get_exodus_entity_name("block", bid)
        if bn and bn.lower() == kelvin_block.lower():
            print("Auto-Kelvin: '%s' block already present — skipping."
                  % kelvin_block)
            return None

    # --- Step 2: locate the air block ---
    air_bid = None
    for bid in cubit.get_block_id_list():
        bn = cubit.get_exodus_entity_name("block", bid)
        if bn and bn.lower() == air_block.lower():
            air_bid = bid
            break
    if air_bid is None:
        raise RuntimeError("Auto-Kelvin needs an '%s' block; none found." % air_block)

    try:
        # --- Step 3: largest surface area among air volumes ---
        air_vols = list(cubit.parse_cubit_list(
            "volume", "in block %d" % air_bid))
        if not air_vols:
            raise RuntimeError("Auto-Kelvin: '%s' block is empty." % air_block)

        best_sid, best_area = 0, 0.0
        for vid in air_vols:
            for sid in cubit.get_relatives("volume", vid, "surface"):
                a = cubit.surface(sid).area()
                if a > best_area:
                    best_area = a
                    best_sid = sid
        if best_sid == 0:
            raise RuntimeError("Auto-Kelvin: no surfaces found on air volumes.")

        # --- Step 4: R from outer-surface vertices, with bbox fallback ---
        # A FULL closed sphere (1/1 case, no symmetry webcut) has zero
        # explicit vertices on its outer surface in Cubit's CAD topology
        # -- the surface wraps around without seams.  Fall back to the
        # surface bounding box in that case.
        vids_outer = list(cubit.get_relatives("surface", best_sid, "vertex"))
        if vids_outer:
            R = max(
                _m.sqrt(sum(c * c for c in cubit.vertex(v).coordinates()))
                for v in vids_outer)
        else:
            bb = cubit.surface(best_sid).bounding_box()
            # bounding_box returns (xmin, ymin, zmin, xmax, ymax, zmax).
            # A sphere's box corner is sqrt(3)*R away and must not be used.
            R = 0.5 * max(bb[3] - bb[0], bb[4] - bb[1], bb[5] - bb[2])
            print("Auto-Kelvin: closed outer surface (no vertices); "
                  "R=%.4f from bounding box." % R)

        # --- Step 5: short-circuit for explicit reduction= ---
        if reduction is not None:
            print("Auto-Kelvin: air R=%.4f m, air_vols=%d, "
                  "reduction=%s, mesh_size=%s"
                  % (R, len(air_vols), reduction, mesh_size))
            info = add_kelvin_cubit(R=R, air_block=air_block,
                                    reduction=reduction,
                                    kelvin_block=kelvin_block,
                                    mesh_size=mesh_size)
            ox, oy, oz = info["center"]
            print("Auto-Kelvin: added at offset=(%.3f, %.3f, %.3f), "
                  "reduction=%s" % (ox, oy, oz, reduction))
            return info

        # --- Step 6: symmetry detection (mesh-seam mode) ---
        # Two cases end up as the same axis-in-symmetry:
        #   (a) half-domain: all vertices on one side of the plane, AND
        #       at least one vertex sits on the plane (min == 0).
        #   (b) webcut-kept-both: the air is split into >=2 volumes along
        #       this specific axis plane (centroids straddle 0).
        # Case (b) used to be detected via "multi-vol AND vertex range
        # straddles 0", but that triggered SPURIOUSLY for any axis on a
        # multi-vol sphere -- e.g. a 1/1 model with only a z-webcut would
        # see x detected too, because the air's vertex range spans x in
        # [-R, +R].  Fix (2026-04-25): require centroids of distinct
        # volumes to straddle 0 on this axis -- this is true ONLY for
        # the actually-webcut axis.
        all_verts = set()
        for vid in air_vols:
            for v in cubit.get_relatives("volume", vid, "vertex"):
                all_verts.add(v)
        # Volume centroids: used to discriminate which axis was webcut.
        try:
            centroids = [cubit.volume(vid).centroid() for vid in air_vols]
        except Exception:
            centroids = []
        # Thresholds scale with R.  Vertex coordinates are exact geometry,
        # so the on-plane test uses a roundoff-level 1e-9 R.  Volume
        # centroids come from Cubit's numerical integration and carry
        # errors near 1e-5 R; an air piece on one side of a seam plane has
        # its centroid a sizeable fraction of R off the plane, so 1e-3 R
        # separates the two without reacting to integration noise.
        plane_tol = 1e-9 * R
        centroid_tol = 1e-3 * R
        symmetry = []
        for axis, name in enumerate(("x", "y", "z")):
            coords = [cubit.vertex(v).coordinates()[axis]
                      for v in all_verts]
            cmin, cmax = min(coords), max(coords)
            on_plane = any(abs(c) < plane_tol for c in coords)
            # Case (a): half-domain
            half_domain = cmin >= -plane_tol and on_plane
            # Case (b): centroids straddle this axis plane.
            cent_coords = [c[axis] for c in centroids]
            has_pos = any(cc > centroid_tol for cc in cent_coords)
            has_neg = any(cc < -centroid_tol for cc in cent_coords)
            webcut_full = on_plane and has_pos and has_neg
            if half_domain or webcut_full:
                symmetry.append(name)

        print("Auto-Kelvin: air R=%.4f m, air_vols=%d, symmetry=%s, "
              "mesh_size=%s" % (R, len(air_vols), symmetry, mesh_size))

        # --- Step 7: delegate to add_kelvin_cubit ---
        info = add_kelvin_cubit(R=R, air_block=air_block,
                                symmetry=symmetry,
                                kelvin_block=kelvin_block,
                                mesh_size=mesh_size)
        ox, oy, oz = info["center"]
        print("Auto-Kelvin: added at offset=(%.3f, %.3f, %.3f), "
              "symmetry=%s" % (ox, oy, oz, symmetry))
        return info
    except Exception as e:
        raise RuntimeError("Auto-Kelvin failed; open-boundary export aborted") from e


SYM_VALUES = {"off": None, "bn": "bn=0", "ht": "ht=0"}


def reduction_from_sym(sym_x="off", sym_y="off", sym_z="off"):
    """Map the per-axis ``off|bn|ht`` options to a ``reduction`` dict or None.

    All ``off`` keeps the automatic mesh-seam symmetry detection (None).
    Any other value raises.
    """
    reduction = {}
    for axis, value in (("x", sym_x), ("y", sym_y), ("z", sym_z)):
        key = str(value).strip().lower()
        if key not in SYM_VALUES:
            raise ValueError(f"kelvin_sym_{axis} must be off|bn|ht, got {value!r}")
        if SYM_VALUES[key] is not None:
            reduction[axis] = SYM_VALUES[key]
    return reduction or None


def validate_options(options):
    """Normalize and check the Kelvin export options; raise on any invalid value."""
    allowed = {"add_kelvin", "kelvin_air", "kelvin_block", "kelvin_mesh",
               "kelvin_sym_x", "kelvin_sym_y", "kelvin_sym_z"}
    unknown = set(options) - allowed
    if unknown:
        raise ValueError(f"unknown Kelvin option(s): {sorted(unknown)}")
    out = {"add_kelvin": bool(options.get("add_kelvin", False)),
           "kelvin_air": str(options.get("kelvin_air") or "air"),
           "kelvin_block": str(options.get("kelvin_block") or "kelvin"),
           "kelvin_mesh": options.get("kelvin_mesh"),
           "kelvin_sym_x": str(options.get("kelvin_sym_x") or "off").lower(),
           "kelvin_sym_y": str(options.get("kelvin_sym_y") or "off").lower(),
           "kelvin_sym_z": str(options.get("kelvin_sym_z") or "off").lower()}
    if out["kelvin_mesh"] is not None:
        mesh = float(out["kelvin_mesh"])
        if not (mesh > 0 and math.isfinite(mesh)):
            raise ValueError(f"kelvin_mesh must be a positive length in m, got {out['kelvin_mesh']!r}")
        out["kelvin_mesh"] = mesh
    out["reduction"] = reduction_from_sym(out["kelvin_sym_x"], out["kelvin_sym_y"], out["kelvin_sym_z"])
    return out


def export_arguments(info, match=None):
    """``export vol`` arguments that make the plugin write the Kelvin pairs.

    ``info`` is the dict ``add_kelvin_cubit`` returns; its ``center`` is the
    Kelvin sphere offset, the translation from each ``kelvin_int`` point to its
    ``kelvin_ext`` copy.  ``match`` is the pairing distance in m (default
    ``1e-6 * R``).  The plugin rejects a Kelvin mesh exported without these
    fields, so the ``.vol`` always carries its ``kelvin`` periodic group.
    """
    offset, tol = pair_fields(info, match)
    return ("kelvin_tx %.17g kelvin_ty %.17g kelvin_tz %.17g kelvin_match %.17g"
            % (offset[0], offset[1], offset[2], tol))


def pair_fields(info, match=None):
    """The Kelvin pairing of :func:`export_arguments` as ``(translation, match)``."""
    if not info or "center" not in info or "R" not in info:
        raise ValueError("Kelvin export arguments need the add_kelvin_cubit info (R, center)")
    R = float(info["R"])
    offset = [float(c) for c in info["center"]]
    tol = 1e-6 * R if match is None else float(match)
    if not (R > 0 and math.isfinite(R) and tol > 0 and math.isfinite(tol)
            and len(offset) == 3 and all(math.isfinite(c) for c in offset)):
        raise ValueError(f"invalid Kelvin export geometry: R={info['R']!r}, "
                         f"center={info['center']!r}, match={match!r}")
    return offset, tol


def apply_options(options):
    """Run Auto-Kelvin on the current Cubit model when ``add_kelvin`` is set."""
    opts = validate_options(options)
    if not opts["add_kelvin"]:
        return None
    return auto_add_kelvin_from_current_model(air_block=opts["kelvin_air"],
                                              kelvin_block=opts["kelvin_block"],
                                              mesh_size=opts["kelvin_mesh"],
                                              reduction=opts["reduction"])


# ====================================================================
# OCC path (Layer 4 -- system Python 3.12 + NGSolve)
# ====================================================================

def add_kelvin_occ(air_shape, R, symmetry=None,
                   offset_dir=None, offset_dist=None, maxh_kelvin=None):
    """Add Kelvin open-boundary sphere to an OCC geometry.

    Call **before** meshing.  Returns the modified compound shape with
    the Kelvin sphere, periodic identification, and GND vertex.

    The caller is responsible for creating the interior air sphere
    (radius R, centered at origin) and subtracting the physical objects
    from it.  This function adds the exterior Kelvin sphere.

    Args:
        air_shape: OCC shape for the air domain (sphere with holes).
            Its outer spherical face (area ~ 4*pi*R^2 or fraction)
            will be identified with the Kelvin sphere surface.
        R: Kelvin sphere radius [m].
        symmetry: List of symmetry planes, e.g. ["z"] or ["x","z"].
        offset_dir: Explicit offset axis or None (auto).
        offset_dist: Explicit offset distance [m] or None (3*R).
        maxh_kelvin: Max mesh size for Kelvin volume [m] or None.

    Returns:
        (compound_shape, info_dict)
        The compound_shape is ready for OCCGeometry().GenerateMesh().
        No scalar point constraint is fabricated from an isolated OCC vertex.
        For Radia scalar H1 use, pass the generated mesh through
        radia.kelvin_geometry.ground_kelvin_scalar_mesh(mesh, info["center"])
        before curving or creating spaces. This fixes an additive gauge, not
        a prescribed physical value at infinity; HCurl needs its own gauge.
        info_dict has keys: 'R', 'center', 'symmetry', 'offset_dir'.
    """
    from netgen.occ import (Sphere, Pnt, Dir, Glue, HalfSpace,
                             IdentificationType)

    if symmetry is None:
        symmetry = []

    ox, oy, oz = compute_offset_vector(R, symmetry, offset_dir, offset_dist)
    actual_dir = auto_offset_direction(symmetry) if offset_dir is None else offset_dir

    sym = sorted(set(s.lower() for s in symmetry))
    if actual_dir == "diag" or offset_dir == "diag":
        raise ValueError(
            "add_kelvin_occ: no offset lies in all three symmetry planes, so "
            "the 1/8 Kelvin sphere cannot be a translated copy of the air "
            "domain; use at most two symmetry planes on the OCC route.")

    # ---- 1. Exterior Kelvin sphere = translated copy of the cut air sphere.
    # The offset lies in every symmetry plane, so cutting at the origin and
    # then moving equals cutting at the offset; building it this way keeps
    # the OCC face split (seams) identical on both sides.
    # netgen.occ.HalfSpace(p, n) keeps the side n points away from, so the
    # inward normal keeps the positive side of each plane.
    ext_sphere = Sphere(Pnt(0, 0, 0), R)
    for plane_name in sym:
        normal = _SYM_NORMALS[plane_name]
        ext_sphere = ext_sphere * HalfSpace(Pnt(0, 0, 0), Dir(*(-c for c in normal)))
    ext_sphere = ext_sphere.Move((ox, oy, oz))
    ext_sphere.name = "kelvin"
    if maxh_kelvin is not None:
        ext_sphere.maxh = maxh_kelvin

    # ---- 2. Pair the sphere faces one by one -------------------------
    normals = [_SYM_NORMALS[n] for n in sym]
    tol = 1e-6 * R

    def sphere_faces(shape, centre):
        found = []
        for f in shape.faces:
            if not all(abs(math.dist(tuple(v.p), centre) - R) < tol for v in f.vertices):
                continue
            c = f.center
            rel = (c.x - centre[0], c.y - centre[1], c.z - centre[2])
            if any(abs(sum(r * n for r, n in zip(rel, nn))) < tol for nn in normals):
                continue  # planar cut face whose rim lies on the sphere
            found.append(f)
        return found

    int_faces = sphere_faces(air_shape, (0.0, 0.0, 0.0))
    ext_faces = sphere_faces(ext_sphere, (ox, oy, oz))
    if not int_faces:
        raise RuntimeError(
            "Cannot find the interior sphere face(s) of radius R = %g at the "
            "origin; air_shape must have a spherical outer boundary on the "
            "positive side of each symmetry plane." % R)
    pairs = []
    for fi in int_faces:
        ci = fi.center
        match = [fe for fe in ext_faces
                 if math.dist((fe.center.x, fe.center.y, fe.center.z),
                              (ci.x + ox, ci.y + oy, ci.z + oz)) < tol]
        if len(match) != 1:
            raise RuntimeError(
                "kelvin_int face at %r has %d translated copies by the offset "
                "on the Kelvin sphere (need exactly 1); the air domain must "
                "lie on the positive side of each symmetry plane."
                % ((ci.x, ci.y, ci.z), len(match)))
        pairs.append((fi, match[0]))
    if len(pairs) != len(ext_faces):
        raise RuntimeError(
            "Kelvin sphere has %d faces but only %d were paired with kelvin_int"
            % (len(ext_faces), len(pairs)))
    for fi, fe in pairs:
        fi.name = "kelvin_int"
        fe.name = "kelvin_ext"
        fi.Identify(fe, "kelvin", IdentificationType.PERIODIC)

    # A detached OCC Vertex does not constrain a volume-mesh unknown.
    # Scalar users must explicitly attach a connected point gauge after meshing.
    shape = Glue([air_shape, ext_sphere])

    info = {
        "R": R,
        "center": (ox, oy, oz),
        "symmetry": list(symmetry),
        "offset_dir": actual_dir,
    }

    print("")
    print("=== add_kelvin_occ ===")
    print("  R = %g m" % R)
    print("  offset = (%g, %g, %g) m" % (ox, oy, oz))
    n_sym = len(sym)
    print("  symmetry = %s (1/%d)" % (symmetry, 2 ** n_sym) if n_sym else
          "  symmetry = [] (full)")
    print("  Identify: kelvin_int <-> kelvin_ext")
    print("  Scalar point gauge: attach to a connected vertex after meshing")

    return shape, info


# ====================================================================
# OCC 2D axisymmetric path (Layer 4 -- system Python 3.12 + NGSolve)
# ====================================================================

def add_kelvin_2d_axisym(interior_shape, R, z_offset=None, maxh_kelvin=None):
    """Add z-offset Kelvin half-circle for 2D axisymmetric FEM.

    Call **before** meshing.  The caller must provide `interior_shape`:
    a 2D OCC face representing the physical domain on x >= 0 (the (r, z)
    projection of the axisymmetric 3D geometry), with its outer curved
    edge (radius R arc centered at origin) ready to be identified.

    Edges of `interior_shape` should already be named:
      - "axis"       on the x = 0 boundary (r = 0)
      - "kelvin_int" on the outer arc (to be identified with exterior)
      - anything else for physical boundaries (coil, wp, etc.)

    Args:
        interior_shape: OCC 2D face for the interior half-circle with
            any physical inclusions subtracted. Must live on x >= 0.
        R: Kelvin arc radius [m]. Must enclose all physical objects.
        z_offset: Vertical offset of exterior half-circle center [m].
            Default: 2.5*R (ensures interior and exterior are disjoint).
        maxh_kelvin: Max mesh size for Kelvin face [m] or None.

    Returns:
        (compound_shape, info_dict)
        The compound_shape is ready for OCCGeometry(shape, dim=2).GenerateMesh().
        info_dict has keys: 'R', 'z_offset', 'kelvin_factor' (string expr),
        'axis_labels' (Dirichlet edges to pass to FES).

    Weak form hint (A-formulation, u = r*A_phi, mesh coordinates x=r, y=z):
        fes = Periodic(H1(mesh, order=p, complex=True,
                          dirichlet="axis|axis_ext",
                          dirichlet_bbnd="GND"))
        from ngsolve import x, y, sqrt, IfPos
        r_safe = IfPos(x - 1e-10, x, 1e-10)
        rho_prime = sqrt(x**2 + (y - z_offset)**2)
        rho_safe = IfPos(rho_prime - 1e-10, rho_prime, 1e-10)
        kelvin_fac = (rho_safe / R)**2
        nu_dict = {m: nu_0 * kelvin_fac if "kelvin" in m else nu_0
                   for m in mesh.GetMaterials()}
        nu_cf = mesh.MaterialCF(nu_dict, default=nu_0)
        a_bf += nu_cf / r_safe * grad(u) * grad(v) * dx
    """
    from netgen.occ import (WorkPlane, MoveTo, Axes, Pnt, X, Z,
                             Glue, Vertex, IdentificationType)

    if z_offset is None:
        z_offset = 2.5 * R
    if z_offset < 2.0 * R:
        raise ValueError(
            "z_offset %.4f must be >= 2*R = %.4f to avoid overlap."
            % (z_offset, 2.0 * R))

    # ---- 1. Build exterior half-circle centered at (0, z_offset) ----
    wp_ext = WorkPlane(Axes((0, z_offset, 0), n=Z, h=X))
    outer_full = wp_ext.Circle(R).Face()
    cutter_ext = MoveTo(-R - 0.1, z_offset - R - 0.1).Rectangle(
        R + 0.1, 2 * R + 0.2).Face()
    exterior = outer_full - cutter_ext
    exterior.name = "kelvin"
    if maxh_kelvin is not None:
        exterior.maxh = maxh_kelvin

    # ---- 2. Name exterior edges ----
    kelvin_ext_edges = []
    for edge in exterior.edges:
        cx = edge.center.x
        try:
            v0, v1 = edge.vertices
            d0 = math.sqrt(v0.p.x ** 2 + (v0.p.y - z_offset) ** 2)
            d1 = math.sqrt(v1.p.x ** 2 + (v1.p.y - z_offset) ** 2)
            is_arc = (abs(d0 - R) < 0.01 * R and abs(d1 - R) < 0.01 * R
                      and cx > 0.01 * R)
        except Exception:
            is_arc = False
        if cx < 1e-4:
            edge.name = "axis_ext"
        elif is_arc:
            edge.name = "kelvin_ext"
            kelvin_ext_edges.append(edge)
        else:
            edge.name = "default"

    # ---- 3. Find kelvin_int edges in interior (expected named already) ----
    kelvin_int_edges = [e for e in interior_shape.edges
                        if getattr(e, 'name', '') == "kelvin_int"]
    if not kelvin_int_edges:
        raise RuntimeError(
            "No edges named 'kelvin_int' found in interior_shape. "
            "The caller must name the outer arc(s) of the interior "
            "half-circle 'kelvin_int' before calling add_kelvin_2d_axisym.")
    if not kelvin_ext_edges:
        raise RuntimeError(
            "Failed to identify exterior Kelvin arc edges. "
            "Check z_offset and R values.")

    # ---- 4. Periodic identification (match by y-sign) ----
    matched = 0
    for int_e in kelvin_int_edges:
        iy = int_e.center.y
        for ext_e in kelvin_ext_edges:
            ey = ext_e.center.y - z_offset
            # Match same-sign halves (Netgen splits full arc into top/bottom)
            if (iy > 0 and ey > 0) or (iy < 0 and ey < 0) or \
               (abs(iy) < 1e-8 and abs(ey) < 1e-8):
                int_e.Identify(ext_e, "kelvin", IdentificationType.PERIODIC)
                matched += 1
                break
    if matched == 0:
        raise RuntimeError(
            "Failed to match any kelvin_int edge with a kelvin_ext edge. "
            "Verify that the interior arc radius matches R = %g." % R)

    # ---- 5. GND vertex at exterior center (= image of infinity) ----
    gnd = Vertex(Pnt(0, z_offset, 0))
    gnd.name = "GND"

    # ---- 6. Glue ----
    shape = Glue([interior_shape, exterior, gnd])

    info = {
        "R": R,
        "z_offset": z_offset,
        "n_periodic_pairs": matched,
        "axis_labels": "axis|axis_ext",
        "kelvin_factor": "(rho'/R)**2 with rho' = sqrt(x**2 + (y-z_offset)**2)",
    }

    print("")
    print("=== add_kelvin_2d_axisym ===")
    print("  R        = %g m" % R)
    print("  z_offset = %g m" % z_offset)
    print("  Periodic pairs matched: %d (kelvin_int <-> kelvin_ext)" % matched)
    print("  GND vertex at (0, %g)" % z_offset)
    print("  Reluctivity factor (exterior): nu_0 * (rho'/R)^2")

    return shape, info
