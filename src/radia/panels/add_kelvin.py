"""Add Kelvin open-boundary sphere to an existing model.

Provides three entry points -- Cubit 3D, OCC 3D, OCC 2D axisym -- with
consistent semantics:

    1.  Create an exterior sphere/half-circle (same radius R, offset in space)
    2.  Handle symmetry cuts (1/2, 1/4)  [3D only]
    3.  Establish 1:1 surface/edge mesh correspondence (Periodic BC)
    4.  Mesh the Kelvin volume/face
    5.  Assign labels (block/sideset or OCC face/edge names)
    6.  Create GND vertex at Kelvin center (= physical infinity)

The caller supplies only the physical-domain geometry and the Kelvin
radius.  Everything else is automatic.

## Two symmetry modes (3D Cubit path)

**webcut (existing, default)** -- the air sphere has a z-plane webcut
but BOTH halves are retained in the mesh.  The webcut is a mesh seam
for Periodic copy-mesh only; the geometry is still a full sphere.
Symmetry of the physical problem (e.g. mirror coils about z=0) is
handled by the source definition (coil + coil-mirror), not by the
mesh.  Use this for field-line-continuation validation.  API::

    add_kelvin_cubit(R=0.06, symmetry=["z"])

**reduction (new, 2026-04-25)** -- the air block already contains only
the reduced domain (x>=0 and/or y>=0 and/or z>=0).  The Kelvin sphere
is cut on the same planes.  Each symmetry plane is given a sideset
label that encodes the boundary-condition *physics* (not the math -
A-formulation Dirichlet and Omega-formulation Dirichlet mean opposite
things for the same symmetry).  The solver then chooses Dirichlet or
Natural per formulation::

    add_kelvin_cubit(R=0.06,
                     reduction={"x": "ht=0", "z": "bn=0"})

    # sideset labels produced:
    #   sym_ht=0_x    -- H x n = 0 on x=0 plane (field perpendicular)
    #                    -> Omega: Dirichlet (Omega=const)
    #                    -> A:     Natural
    #   sym_bn=0_z    -- B . n  = 0 on z=0 plane (field parallel)
    #                    -> Omega: Natural
    #                    -> A:     Dirichlet (A x n = 0)

The two modes are mutually exclusive.  Passing both `symmetry=...` and
`reduction=...` raises.

Usage (Cubit Python -- Layer 2, inside Cubit):
    >>> from add_kelvin import add_kelvin_cubit
    >>> info = add_kelvin_cubit(R=0.06, symmetry=["z"])
    >>> info = add_kelvin_cubit(R=0.06, reduction={"x": "ht=0", "z": "bn=0"})

Usage (OCC 3D -- Layer 4, system Python 3.12):
    >>> from add_kelvin import add_kelvin_occ
    >>> shape, info = add_kelvin_occ(air_shape, R=0.06, symmetry=["z"])

Usage (OCC 2D axisymmetric -- Layer 4, system Python 3.12):
    >>> from add_kelvin import add_kelvin_2d_axisym
    >>> shape, info = add_kelvin_2d_axisym(interior_face, R=0.10, z_offset=0.25)

## 2D Axisymmetric Kelvin (z-offset strategy)

In 2D axisymmetric FEM (r, z plane representation of an axisymmetric
3D problem), the physical domain lives on x >= 0 with r = x, z = y.
The Kelvin transformation maps the infinite exterior (r^2 + z^2 > R^2)
into a bounded fictitious half-circle centered at (0, z_offset):

    interior:  half-circle centered (0, 0),        radius R, x >= 0
    exterior:  half-circle centered (0, z_offset), radius R, x >= 0
    link:      Periodic BC along the curved arcs (kelvin_int, kelvin_ext)
    GND:       vertex at (0, z_offset) (= image of infinity, Dirichlet)

Reluctivity modulation (A-formulation, u = r*A_phi):
    nu_exterior(x, y) = nu_0 * (rho' / R)^2
    where rho' = sqrt(x^2 + (y - z_offset)^2)

This factor is 0 at rho'=0 (infinity) and continuous nu_0 at rho'=R
(Kelvin boundary), providing the "image of infinity" degeneracy that
requires GND Dirichlet for uniqueness.

IMPORTANT: The 2D axisym Kelvin radius must enclose ALL physical
objects (coil + workpiece + sources). In practice R >= 3*R_coil is
recommended. A recent convergence test (2026-04-14) verified that L
is stable within 0.07% across R in [0.08, 0.25]m for a R_coil=30mm
coil + R_wp=25mm workpiece.
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
#              Radia image sign: '+' (mirror-symmetric).
#   "ht=0"  -- H x n = 0 on the plane (flux perpendicular to plane).
#              Omega: Dirichlet (Omega=const).  A: natural.
#              Radia image sign: '-' (mirror-antisymmetric).
_VALID_BC = ("bn=0", "ht=0")


def sym_sideset_name(axis, bc):
    """Return the canonical sideset label for a symmetry plane.

    Args:
        axis: one of "x", "y", "z" (lowercase).
        bc: one of "bn=0", "ht=0".

    Returns:
        e.g. "sym_bn=0_x", "sym_ht=0_z".  This string is the single
        source of truth -- `calc_accel_magnet`/`calc_accel_msc` both
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

def add_kelvin_cubit(R, air_block="air", symmetry=None, reduction=None,
                     offset_dir=None, offset_dist=None, mesh_size=None,
                     kelvin_block="kelvin", gnd_nodeset=100):
    """Add Kelvin open-boundary sphere to the current Cubit model.

    Call **after** meshing the physical domain (coil, workpiece, air).
    Creates a fresh exterior sphere at the offset position, webcutted
    for symmetry, with 1:1 mesh copy from the air sphere outer surface.

    Strategy (matches the verified examples/induction_heating/
    demoted_samples/ih_fem_kelvin_sample.py):
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
    import cubit  # Layer 2 only

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

    # ---- 2. Create exterior sphere (clean ACIS sphere) ----
    cubit.cmd("create sphere radius %g" % R)
    kelvin_sphere = cubit.get_last_id("volume")
    cubit.cmd("move volume %d x %g y %g z %g" % (kelvin_sphere, ox, oy, oz))

    # ---- 3. Webcut for symmetry ----
    sym = set(s.lower() for s in symmetry)
    kelvin_vols = [kelvin_sphere]

    if "z" in sym:
        # Webcut with zplane to match the air sphere's equator webcut
        id_before = cubit.get_last_id("volume")
        cubit.cmd("webcut volume %d with plane zplane" % kelvin_sphere)
        kelvin_top = kelvin_sphere
        kelvin_bot = cubit.get_last_id("volume")
        kelvin_vols = [kelvin_top, kelvin_bot]

    # Imprint + merge kelvin halves (NOT with air volumes!)
    if len(kelvin_vols) > 1:
        k_str = " ".join(str(v) for v in kelvin_vols)
        cubit.cmd("imprint volume %s" % k_str)
        cubit.cmd("merge volume %s" % k_str)

    for vid in kelvin_vols:
        cubit.cmd('volume %d rename "kelvin_ext"' % vid)

    # ---- 4. Find outer spherical surface on each air volume ----
    # After subtract of coil/workpiece, the air volumes still have the
    # original sphere outer surface as their LARGEST surface.
    air_outer_surfs = []
    for vid in air_vols:
        surfs = list(cubit.get_relatives("volume", vid, "surface"))
        if surfs:
            outer = max(surfs, key=lambda s: cubit.surface(s).area())
            air_outer_surfs.append(outer)

    # ---- 5. Find hemisphere surface on each kelvin volume ----
    kelvin_outer_surfs = []
    for vid in kelvin_vols:
        surfs = list(cubit.get_relatives("volume", vid, "surface"))
        if surfs:
            outer = max(surfs, key=lambda s: cubit.surface(s).area())
            kelvin_outer_surfs.append(outer)

    # ---- 6. Copy mesh: air outer surface -> kelvin surface ----
    # Pair by matching air/kelvin volume order (both webcutted at z=0).
    n_pairs = min(len(air_outer_surfs), len(kelvin_outer_surfs))
    for i in range(n_pairs):
        src = air_outer_surfs[i]
        dst = kelvin_outer_surfs[i]

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
                if s in live_surfs:
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
    import cubit

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
    if len(axes) == 3:
        raise NotImplementedError(
            "1/8 reduction (x+y+z) is not yet supported.  The Kelvin "
            "offset along the (1,1,1) diagonal does not intersect the "
            "x=0/y=0/z=0 planes, so the Kelvin cut faces cannot be "
            "coincident with the air cut faces.  Needs a different "
            "Kelvin mapping (e.g. 7 reflected copies).  Track at "
            "src/radia/panels/samples/em/README.md.")
    if len(axes) not in (1, 2):
        raise ValueError(
            f"reduction supports 1 or 2 axes (1/2 or 1/4), got {len(axes)}")

    # Offset direction: must be perpendicular to all reduction-axis
    # normals, i.e. along a FREE axis not listed in reduction.
    if offset_dir is None:
        free = [d for d in ("x", "y", "z") if d not in axes]
        if not free:
            raise NotImplementedError("1/8 case -- see above")
        offset_dir = free[0]
    elif offset_dir.lower() in axes:
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

    idx = {"x": 0, "y": 1, "z": 2}[offset_dir.lower()]
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

    # Webcut keeps the "positive" side of each plane by discarding the
    # negative-side volume produced by `webcut ... with plane <axis>plane`.
    kelvin_vol = kelvin_sphere
    for axis in axes:
        plane = f"{axis}plane"
        before = set(cubit.parse_cubit_list("volume", "all"))
        cubit.cmd(f"webcut volume {kelvin_vol} with plane {plane}")
        after = set(cubit.parse_cubit_list("volume", "all"))
        new_vols = list(after - before)
        # Identify the two halves: one has centroid coord > 0 (keep),
        # the other < 0 (delete).  After webcut, kelvin_vol is one of
        # {kelvin_vol (unchanged id), the new id}.
        keep, drop = None, None
        candidates = [kelvin_vol] + new_vols
        ai = {"x": 0, "y": 1, "z": 2}[axis]
        for vid in candidates:
            try:
                c = cubit.volume(vid).centroid()[ai]
            except Exception:
                continue
            if c > 0:
                keep = vid
            elif c < 0:
                drop = vid
        if keep is None or drop is None:
            raise RuntimeError(
                f"webcut on {plane} did not produce a positive/negative "
                f"pair (candidates={candidates}).")
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
    k_surfs = list(cubit.get_relatives("volume", kelvin_vol, "surface"))
    k_curved = None
    k_flat_by_axis = {}
    for sid in k_surfs:
        bb = cubit.surface(sid).bounding_box()
        # bounding_box returns (xmin, ymin, zmin, xmax, ymax, zmax, ...)
        xmin, ymin, zmin = bb[0], bb[1], bb[2]
        xmax, ymax, zmax = bb[3], bb[4], bb[5]
        spans = ((xmax - xmin), (ymax - ymin), (zmax - zmin))
        # Flat face on axis=0 plane: the span along that axis is ~0 AND
        # the coordinate is ~0 on that plane.
        flat_axis = None
        for ai, ax_name in enumerate(("x", "y", "z")):
            lo, hi = bb[ai], bb[3 + ai]
            if ax_name in axes and abs(hi - lo) < 1e-5 * R and abs(lo) < 1e-5 * R:
                flat_axis = ax_name
                break
        if flat_axis is not None:
            k_flat_by_axis[flat_axis] = sid
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
        if axis not in k_flat_by_axis:
            raise RuntimeError(
                f"Cannot find flat sym face on axis {axis!r} for the "
                f"Kelvin volume (expected face with {axis}=0 plane).  "
                f"Surfaces = {k_surfs}.")

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
    a_c = max(a_curves, key=lambda c: cubit.curve(c).length())
    k_c = max(k_curves, key=lambda c: cubit.curve(c).length())
    a_v = cubit.get_relatives("curve", a_c, "vertex")[0]
    k_v = cubit.get_relatives("curve", k_c, "vertex")[0]
    cubit.cmd(
        f"copy mesh surface {a_curved} onto surface {k_curved} "
        f"source curve {a_c} source vertex {a_v} "
        f"target curve {k_c} target vertex {k_v}")

    # ---- 7. Mesh the Kelvin volume ----
    cubit.cmd(f"volume {kelvin_vol} scheme tetmesh")
    if mesh_size is not None:
        cubit.cmd(f"volume {kelvin_vol} size {float(mesh_size):g}")
    cubit.cmd(f"mesh volume {kelvin_vol}")

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

    # Per-axis sym sidesets -- each includes BOTH the air cut face(s)
    # and the Kelvin cut face.  The solver sees a single boundary
    # labelled `sym_<bc>_<axis>` and applies the appropriate BC.
    sym_sidesets = {}
    for axis in axes:
        bc = reduction[axis] if axis in reduction else reduction[axis.lower()]
        label = sym_sideset_name(axis, bc)
        if label.lower() in existing_ss_names:
            sym_sidesets[axis] = None
            continue
        faces = list(a_flat_by_axis.get(axis, []))
        faces.append(k_flat_by_axis[axis])
        if not faces:
            raise RuntimeError(
                f"No faces found for sym axis={axis!r}.")
        sid = _next_ss()
        cubit.cmd(f"sideset {sid} add surface {' '.join(str(f) for f in faces)}")
        cubit.cmd(f'sideset {sid} name "{label}"')
        sym_sidesets[axis] = sid

    # ---- 10. GND nodeset at Kelvin center ----
    cubit.cmd(f"create vertex {ox:g} {oy:g} {oz:g}")
    gnd_vid = cubit.get_last_id("vertex")
    cubit.cmd(f"nodeset {gnd_nodeset} add vertex {gnd_vid}")
    cubit.cmd(f'nodeset {gnd_nodeset} name "GND"')

    cubit.cmd(f"volume {kelvin_vol} visibility off")

    frac = "1/%d" % (2 ** len(axes))
    print("")
    print("=== add_kelvin_cubit (reduction mode) ===")
    print(f"  R = {R:g} m,  offset = ({ox:g}, {oy:g}, {oz:g}) m")
    print(f"  reduction = {reduction} ({frac} model)")
    print(f"  air_vols = {air_vols}, kelvin_vol = {kelvin_vol}")
    print(f"  sideset 'kelvin_int' (air curved) / 'kelvin_ext' (kelvin curved)")
    for axis, sid in sym_sidesets.items():
        label = sym_sideset_name(axis, reduction[axis])
        print(f"  sideset {sid} '{label}'  (air cut + kelvin cut on {axis}=0)")
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
    }


def auto_add_kelvin_from_current_model(air_block="air",
                                        kelvin_block="kelvin",
                                        mesh_size=None,
                                        reduction=None):
    """Detect air sphere + symmetry, then call add_kelvin_cubit().

    Meant to be invoked by the Radia-NGSolve launcher just before
    `radia_export netgen`.  Runs inside Cubit's embedded Python.

    Steps (matching the 2026-04-14 c60a6007 implementation):
      1. If a `<kelvin_block>` block already exists, skip (idempotent).
      2. Find the `<air_block>` block.  Abort with a warning if missing.
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

    Returns the info dict from ``add_kelvin_cubit``, or None on skip /
    failure.  Never raises — the launcher should continue (and fall
    back to Dirichlet truncation) if auto-detection fails.

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
    try:
        import cubit
    except ImportError:
        print("WARNING: auto_add_kelvin_from_current_model requires the "
              "cubit Python module (must run inside Cubit).")
        return None

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
        print("WARNING: Auto-Kelvin needs an '%s' block.  "
              "None found — skipping." % air_block)
        return None

    try:
        # --- Step 3: largest surface area among air volumes ---
        air_vols = list(cubit.parse_cubit_list(
            "volume", "in block %d" % air_bid))
        if not air_vols:
            print("WARNING: Auto-Kelvin: '%s' block is empty." % air_block)
            return None

        best_sid, best_area = 0, 0.0
        for vid in air_vols:
            for sid in cubit.get_relatives("volume", vid, "surface"):
                a = cubit.surface(sid).area()
                if a > best_area:
                    best_area = a
                    best_sid = sid
        if best_sid == 0:
            print("WARNING: Auto-Kelvin: no surfaces found on air volumes.")
            return None

        # --- Step 4: R from outer-surface vertices ---
        vids_outer = cubit.get_relatives("surface", best_sid, "vertex")
        R = max(
            _m.sqrt(sum(c * c for c in cubit.vertex(v).coordinates()))
            for v in vids_outer)

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
        #       the plane but both halves are retained (full sphere).
        #       Detected by: multi-volume air AND vertex range straddles
        #       0 symmetrically AND a vertex sits on the plane.
        # In case (b), passing the axis as "symmetry" to add_kelvin_cubit
        # is what triggers the matching webcut on the Kelvin side, so
        # the copy-mesh step can pair up the equator curves.
        all_verts = set()
        for vid in air_vols:
            for v in cubit.get_relatives("volume", vid, "vertex"):
                all_verts.add(v)
        symmetry = []
        multi_vol = len(air_vols) > 1
        for axis, name in enumerate(("x", "y", "z")):
            coords = [cubit.vertex(v).coordinates()[axis]
                      for v in all_verts]
            cmin, cmax = min(coords), max(coords)
            on_plane = any(abs(c) < 1e-6 for c in coords)
            # Case (a): half-domain
            half_domain = cmin >= -1e-6 and on_plane
            # Case (b): full-domain with webcut equator kept
            scale = max(abs(cmin), abs(cmax), 1e-12)
            webcut_full = (multi_vol and on_plane
                           and abs(cmin + cmax) < 1e-6 * scale
                           and cmin < -1e-6)
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
        print("WARNING: Auto-Kelvin failed: %s" % e)
        print("Proceeding without Kelvin (Dirichlet truncation on outer "
              "boundary).")
        return None


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
        info_dict has keys: 'R', 'center', 'symmetry', 'offset_dir'.
    """
    from netgen.occ import (Sphere, Pnt, Dir, Glue, Vertex, HalfSpace,
                             IdentificationType)

    if symmetry is None:
        symmetry = []

    ox, oy, oz = compute_offset_vector(R, symmetry, offset_dir, offset_dist)
    actual_dir = auto_offset_direction(symmetry) if offset_dir is None else offset_dir

    # ---- 1. Create exterior Kelvin sphere ----
    ext_sphere = Sphere(Pnt(ox, oy, oz), R)
    ext_sphere.name = "kelvin"
    if maxh_kelvin is not None:
        ext_sphere.maxh = maxh_kelvin

    # ---- 2. Cut with symmetry planes ----
    sym = set(s.lower() for s in symmetry)
    for plane_name in sorted(sym):
        normal = _SYM_NORMALS[plane_name]
        # Keep the positive side (normal points outward from cut)
        hs = HalfSpace(Pnt(0, 0, 0), Dir(*normal))
        ext_sphere = ext_sphere * hs
        ext_sphere.name = "kelvin"

    # ---- 3. Identify periodic faces ----
    # Interior sphere outer face: area ~ expected fraction of 4*pi*R^2
    n_sym = len(symmetry)
    expected_area = 4.0 * math.pi * R ** 2 / (2 ** n_sym)

    int_face = None
    for f in air_shape.faces:
        if abs(f.mass - expected_area) / expected_area < 0.15:
            f.name = "kelvin_int"
            int_face = f
            break

    ext_face = None
    for f in ext_sphere.faces:
        if abs(f.mass - expected_area) / expected_area < 0.15:
            f.name = "kelvin_ext"
            ext_face = f
            break

    if int_face is None:
        raise RuntimeError(
            "Cannot find interior sphere face (expected area ~ %.4f m^2).  "
            "Make sure air_shape has a spherical outer boundary with "
            "radius R = %g." % (expected_area, R))
    if ext_face is None:
        raise RuntimeError(
            "Cannot find exterior sphere face (expected area ~ %.4f m^2)."
            % expected_area)

    int_face.Identify(ext_face, "kelvin", IdentificationType.PERIODIC)

    # ---- 4. GND vertex at Kelvin center ----
    gnd = Vertex(Pnt(ox, oy, oz))
    gnd.name = "GND"

    # ---- 5. Glue everything ----
    shape = Glue([air_shape, ext_sphere, gnd])

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
    print("  symmetry = %s (1/%d)" % (symmetry, 2 ** n_sym) if n_sym else
          "  symmetry = [] (full)")
    print("  Identify: kelvin_int <-> kelvin_ext")
    print("  GND vertex at (%g, %g, %g)" % (ox, oy, oz))

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
