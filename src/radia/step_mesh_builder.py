"""Build FEM mesh from STEP file with air domain and Kelvin transformation.

Workflow:
    STEP (yoke) -> OCC import -> air sphere + Kelvin exterior domain
    -> auto-label faces -> Netgen all-tet mesh -> mesh.Curve(order)

Designed for accelerator electromagnet analysis where:
    - Yoke geometry comes from any CAD tool as STEP
    - Air domain and Kelvin exterior domain are added automatically
    - Symmetry planes are auto-detected or specified
    - All-tet mesh avoids pyramid element issues with high-order

Usage:
    from step_mesh_builder import build_mesh_from_step
    mesh = build_mesh_from_step("yoke.step", symmetry="quarter_xz")
"""

import numpy as np
from netgen.occ import OCCGeometry, Glue, Pnt, Vec, Sphere, Box, HalfSpace


def build_mesh_from_step(step_file, symmetry="quarter_xz",
                         kelvin_radius=None, kelvin_factor=2.0,
                         mesh_size=None, mesh_size_yoke=None,
                         mesh_size_air=None, mesh_size_kelvin=None,
                         curve_order=2, fes_order=None,
                         add_kelvin=True):
    """Build NGSolve mesh from yoke STEP file.

    Parameters
    ----------
    step_file : str
        Path to STEP file containing yoke geometry.
    symmetry : str
        Symmetry mode:
        - "full": no symmetry
        - "half_x": symmetry at x=0 (yoke in x<0)
        - "half_z": symmetry at z=0 (yoke in z<0)
        - "quarter_xz": symmetry at x=0 and z=0 (yoke in x<0, z<0)
    kelvin_radius : float or None
        Inner Kelvin sphere radius in meters. If None, auto from bounding box.
    kelvin_factor : float
        Kelvin outer/inner radius ratio (default 2.0).
    mesh_size : float or None
        Global mesh size. If None, auto from bounding box.
    mesh_size_yoke : float or None
        Mesh size for yoke. Overrides global.
    mesh_size_air : float or None
        Mesh size for air. Overrides global.
    mesh_size_kelvin : float or None
        Mesh size for Kelvin exterior domain. Overrides global.
    curve_order : int
        Geometry curving order (default 2).
    fes_order : int or None
        Not used here, but stored in mesh metadata.
    add_kelvin : bool
        If True, add Kelvin transformation shell.

    Returns
    -------
    mesh : ngsolve.Mesh
        NGSolve mesh with labeled boundaries and materials.
    info : dict
        Mesh info: kelvin_radius, kelvin_center, symmetry, etc.
    """
    from ngsolve import Mesh

    # Load STEP
    occ_geo = OCCGeometry(step_file)
    yoke_shape = occ_geo.shape

    # Get yoke bounding box
    bb = yoke_shape.bounding_box
    bb_min = np.array([bb[0].x, bb[0].y, bb[0].z])
    bb_max = np.array([bb[1].x, bb[1].y, bb[1].z])
    bb_center = 0.5 * (bb_min + bb_max)
    bb_size = bb_max - bb_min

    print(f"Yoke BBox: min={bb_min}, max={bb_max}")
    print(f"Yoke center: {bb_center}, size: {bb_size}")

    # Determine symmetry planes
    sym_planes = _parse_symmetry(symmetry, bb_min, bb_max)

    # Determine Kelvin sphere parameters
    if kelvin_radius is None:
        # Auto: sphere must enclose yoke with margin
        max_extent = np.max(np.abs(np.concatenate([bb_min, bb_max])))
        kelvin_radius = max_extent * 1.5
    kelvin_center = _kelvin_center(sym_planes)
    r_inner = kelvin_radius
    r_outer = kelvin_radius * kelvin_factor

    print(f"Kelvin: center={kelvin_center}, r_inner={r_inner:.4f}, "
          f"r_outer={r_outer:.4f}")

    # Auto mesh sizes
    if mesh_size is None:
        mesh_size = np.min(bb_size) / 3
    if mesh_size_yoke is None:
        mesh_size_yoke = mesh_size
    if mesh_size_air is None:
        mesh_size_air = mesh_size * 2.0
    if mesh_size_kelvin is None:
        mesh_size_kelvin = r_inner * 0.5

    # Build OCC geometry with domains
    geo = _build_occ_geometry(
        yoke_shape, sym_planes, kelvin_center,
        r_inner, r_outer, add_kelvin,
        mesh_size_yoke, mesh_size_air, mesh_size_kelvin,
        bb_min, bb_max)

    # Generate mesh
    print(f"Meshing (all tet)...")
    ngmesh = geo.GenerateMesh(maxh=mesh_size_kelvin if add_kelvin else mesh_size_air)
    mesh = Mesh(ngmesh)

    # Curve
    if curve_order > 1:
        print(f"Curving order {curve_order}...")
        mesh.Curve(curve_order)

    # Print mesh info
    mats = set(mesh.GetMaterials())
    bnds = set(mesh.GetBoundaries())
    ne = mesh.ne
    print(f"Mesh: {ne} elements, materials={mats}, boundaries={bnds}")

    info = {
        'kelvin_radius': r_inner,
        'kelvin_center': kelvin_center.tolist(),
        'kelvin_factor': kelvin_factor,
        'symmetry': symmetry,
        'sym_planes': {k: v for k, v in sym_planes.items()},
        'ne': ne,
        'curve_order': curve_order,
        'bb_min': bb_min.tolist(),
        'bb_max': bb_max.tolist(),
    }

    return mesh, info


def _parse_symmetry(symmetry, bb_min, bb_max):
    """Parse symmetry string into plane definitions.

    Returns dict of {axis: coordinate} for symmetry planes.
    E.g., {"x": 0.0, "z": 0.0} for quarter_xz.
    """
    sym = {}
    if symmetry == "full":
        pass
    elif symmetry == "half_x":
        sym["x"] = _detect_sym_plane("x", bb_min, bb_max)
    elif symmetry == "half_z":
        sym["z"] = _detect_sym_plane("z", bb_min, bb_max)
    elif symmetry == "quarter_xz":
        sym["x"] = _detect_sym_plane("x", bb_min, bb_max)
        sym["z"] = _detect_sym_plane("z", bb_min, bb_max)
    else:
        raise ValueError(f"Unknown symmetry: {symmetry}")
    return sym


def _detect_sym_plane(axis, bb_min, bb_max):
    """Detect symmetry plane coordinate.

    For a full model symmetric about 0: bb_min ≈ -bb_max → plane at 0.
    For a quarter model: one end is near 0 → plane there.
    """
    idx = {"x": 0, "y": 1, "z": 2}[axis]
    lo, hi = bb_min[idx], bb_max[idx]

    # Check if model is symmetric about 0 (full model)
    if abs(lo + hi) < 0.1 * abs(hi - lo):
        return 0.0

    # Quarter model: plane at the end closer to zero
    if abs(lo) < abs(hi):
        val = lo
    else:
        val = hi
    if abs(val) < 1e-6:
        return 0.0
    return float(val)


def _kelvin_center(sym_planes):
    """Kelvin sphere center at intersection of symmetry planes."""
    center = np.array([0.0, 0.0, 0.0])
    for axis, val in sym_planes.items():
        idx = {"x": 0, "y": 1, "z": 2}[axis]
        center[idx] = val
    return center


def _build_occ_geometry(yoke, sym_planes, kelvin_center,
                        r_inner, r_outer, add_kelvin,
                        maxh_yoke, maxh_air, maxh_kelvin,
                        yoke_bb_min, yoke_bb_max):
    """Build OCC compound geometry with labeled domains.

    Creates: yoke + air + kelvin (optional), with face labels.

    Strategy: build full compound first, then cut to quarter.
    This avoids Boolean failures when cutting yoke and sphere separately.
    """
    kc = Pnt(*kelvin_center)

    # Create spheres (full, not yet cut)
    sphere_inner = Sphere(kc, r_inner)

    # Air = full inner sphere - full yoke
    air = sphere_inner - yoke

    # Name domains
    yoke.name = "yoke"
    for s in yoke.solids:
        s.name = "yoke"
    yoke.maxh = maxh_yoke

    air.name = "air"
    for s in air.solids:
        s.name = "air"
    air.maxh = maxh_air

    parts = [yoke, air]

    if add_kelvin:
        sphere_outer = Sphere(kc, r_outer)

        # Kelvin = outer - inner
        kelvin = sphere_outer - sphere_inner
        kelvin.name = "kelvin"
        for s in kelvin.solids:
            s.name = "kelvin"
        kelvin.maxh = maxh_kelvin
        parts.append(kelvin)

    # Glue parts together (shared faces) on full model
    compound = Glue(parts)

    # Cut entire compound to quarter model
    if sym_planes:
        compound = _cut_symmetry(compound, sym_planes,
                                  yoke_bb_min, yoke_bb_max)

    # Label faces
    _label_faces(compound, sym_planes, kelvin_center, r_inner, r_outer,
                 add_kelvin)

    geo = OCCGeometry(compound)
    return geo


def _cut_symmetry(shape, sym_planes, yoke_bb_min, yoke_bb_max):
    """Cut shape at symmetry planes to keep the side where the yoke is.

    Detects which side of each symmetry plane the yoke is on,
    and removes the opposite half.
    """
    for axis, val in sym_planes.items():
        idx = {"x": 0, "y": 1, "z": 2}[axis]

        # Determine which side of the plane the yoke is on
        yoke_center = 0.5 * (yoke_bb_min[idx] + yoke_bb_max[idx])
        yoke_on_positive = yoke_center > val

        # OCC HalfSpace(origin, normal) represents the region OPPOSITE to normal.
        # shape - HalfSpace(+x) removes x<0, keeps x>0.
        # So normal should point TOWARD the side we want to keep.
        if axis == "x":
            normal = Vec(1, 0, 0) if yoke_on_positive else Vec(-1, 0, 0)
            origin = Pnt(val, 0, 0)
        elif axis == "y":
            normal = Vec(0, 1, 0) if yoke_on_positive else Vec(0, -1, 0)
            origin = Pnt(0, val, 0)
        elif axis == "z":
            normal = Vec(0, 0, 1) if yoke_on_positive else Vec(0, 0, -1)
            origin = Pnt(0, 0, val)

        hp = HalfSpace(origin, normal)
        shape = shape - hp
    return shape


def _label_faces(compound, sym_planes, kelvin_center, r_inner, r_outer,
                 add_kelvin):
    """Label faces based on geometric properties.

    Labels:
        sym_normal: B perpendicular to face (Omega Dirichlet)
        sym_tangential: B parallel to face (A-form Dirichlet)
        kelvin_int: inner Kelvin sphere surface
        kelvin_ext: outer Kelvin sphere surface (= outer boundary)
        outer: outer boundary (if no Kelvin)
    """
    kc = np.array(kelvin_center)
    tol_plane = max(r_inner * 0.01, 1e-4)  # tolerance for plane detection
    tol_r = 0.05  # relative tolerance for sphere radius

    for face in compound.faces:
        c = np.array([face.center.x, face.center.y, face.center.z])
        bb = face.bounding_box
        fb_min = np.array([bb[0].x, bb[0].y, bb[0].z])
        fb_max = np.array([bb[1].x, bb[1].y, bb[1].z])
        fb_extent = fb_max - fb_min

        # Check symmetry planes: face center near plane coordinate
        labeled = False
        for axis, val in sym_planes.items():
            idx = {"x": 0, "y": 1, "z": 2}[axis]
            face_on_plane = abs(c[idx] - val) < tol_plane
            if face_on_plane:
                # Default labeling by axis:
                #   x=0: sym_tangential (Omega: natural, dOmega/dx=0)
                #   z=0: sym_normal (Omega: Dirichlet, Omega=0)
                # This matches IMA '+x-z' convention for dipole magnets
                if axis == "z":
                    face.name = "sym_normal"
                elif axis == "x":
                    face.name = "sym_tangential"
                else:
                    face.name = f"sym_{axis}"
                labeled = True
                break

        if labeled:
            continue

        # Check sphere surfaces using vertex distances
        verts = face.vertices
        if len(verts) > 0:
            dists = []
            for v in verts:
                vp = np.array([v.p.x, v.p.y, v.p.z])
                dists.append(np.linalg.norm(vp - kc))
            mean_dist = np.mean(dists)
            dist_var = np.std(dists) / mean_dist if mean_dist > 0 else 1.0
        else:
            mean_dist = np.linalg.norm(c - kc)
            dist_var = 0.0

        # Spherical face: all vertices at similar distance from center
        is_spherical = dist_var < 0.1

        if add_kelvin:
            if is_spherical and abs(mean_dist - r_inner) / r_inner < tol_r:
                face.name = "kelvin_int"
                continue
            if is_spherical and abs(mean_dist - r_outer) / r_outer < tol_r:
                face.name = "kelvin_ext"
                continue
        else:
            if is_spherical and abs(mean_dist - r_inner) / r_inner < tol_r:
                face.name = "outer"
                continue

        # Remaining faces: yoke-air interface (no special label needed)


if __name__ == "__main__":
    import sys
    import os

    if len(sys.argv) < 2:
        # Default test with yoke.step
        step_file = os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            '..', '..', 'panels', 'samples', 'em', 'c_type_dipole',
            'yoke.step')
    else:
        step_file = sys.argv[1]

    print(f"Building mesh from: {step_file}")
    mesh, info = build_mesh_from_step(step_file, symmetry="quarter_xz",
                                      curve_order=2)

    print(f"\nResult:")
    print(f"  Elements: {info['ne']}")
    print(f"  Kelvin radius: {info['kelvin_radius']:.4f} m")
    print(f"  Kelvin center: {info['kelvin_center']}")
    print(f"  Materials: {set(mesh.GetMaterials())}")
    print(f"  Boundaries: {set(mesh.GetBoundaries())}")
