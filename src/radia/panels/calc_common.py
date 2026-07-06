"""
Shared utilities for panel calc scripts (BEM + FEM).

Provides:
  - Cubit initialization with path cleanup
  - OCC fallback mesh generation (no Cubit required)
  - ESIM/SIBC solver creation
  - BH curve constants and loading
  - Subprocess output protocol (progress / calc_main)
"""

import io
import json
import math
import os
import sys
import traceback

import numpy as np

MU_0 = 4e-7 * np.pi
NU_0 = 1.0 / MU_0

# BH curve data and material presets are canonical in em_material.py.
# Re-export here for backward compatibility.
_radia_src = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
if os.path.abspath(_radia_src) not in sys.path:
    sys.path.insert(0, os.path.abspath(_radia_src))
from em_material import STEEL_BH, EMMaterial, add_material_args  # noqa: E402


def setup_paths():
    """Add radia src to sys.path (for esim_cell_problem, etc.)."""
    radia_src = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
    abs_src = os.path.abspath(radia_src)
    if abs_src not in sys.path:
        sys.path.insert(0, abs_src)


def setup_cubit(cub5_path=None, load_plugins=False):
    """Initialize Cubit and optionally open a .cub5 file.

    Handles path cleanup to avoid scipy/numpy DLL conflicts.

    Args:
        cub5_path: Optional .cub5 file to open after init.
        load_plugins: If True, set CUBIT_PLUGIN_DIR so radia .ccm
            plugin commands (export gmsh/nastran/vtk) are available.

    Returns:
        cubit module (or None if Cubit unavailable)
    """
    setup_paths()
    try:
        from install_panels import find_cubit_bin
    except ImportError:
        return None

    cubit_path = find_cubit_bin()
    if not cubit_path:
        return None

    if cubit_path not in sys.path:
        sys.path.append(cubit_path)

    # Set CUBIT_PLUGIN_DIR before cubit.init() so .ccm plugins are loaded
    if load_plugins:
        plugin_dir = os.path.join(cubit_path, "plugins")
        if os.path.isdir(plugin_dir):
            os.environ["CUBIT_PLUGIN_DIR"] = plugin_dir

    # Remove Cubit's bundled site-packages (scipy/numpy conflicts)
    _clean_cubit_paths()

    import contextlib
    import cubit
    with contextlib.redirect_stderr(io.StringIO()):
        try:
            cubit.init(["cubit", "-nojournal", "-batch", "-noinit"])
        except Exception:
            pass  # Already initialized

    _clean_cubit_paths()

    if cub5_path:
        cubit.cmd(f'open "{cub5_path}"')

    return cubit


def _clean_cubit_paths():
    """Remove Cubit's bundled site-packages from sys.path."""
    for p in list(sys.path):
        if "site-packages" in p and ("cubit" in p.lower() or "Cubit" in p):
            sys.path.remove(p)


def export_mesh(cubit_mod, order=2, surface_only=False, split_quads=None):
    """Export Cubit mesh to NGSolve via C++ export netgen command.

    Args:
        cubit_mod: Cubit module (already initialized with model open)
        order: Curve order (1-5)
        surface_only: True for BEM surface mesh (auto-extracts from volume)
        split_quads: Unused (kept for API compatibility)

    Returns:
        ngsolve.Mesh
    """
    import tempfile
    setup_paths()
    from ngsolve import Mesh as NGMesh
    vol_path = tempfile.mktemp(suffix='.vol')
    cubit_mod.cmd(f'export netgen "{vol_path}" order {order} overwrite')
    return NGMesh(vol_path)


# ============================================================
# OCC Fallback Mesh Generation (no Cubit required)
# ============================================================

def build_occ_ih_mesh_3d(R_coil=0.030, a_coil=0.003, gap_deg=5,
                         R_wp=0.010, H_wp=0.020,
                         a_kelvin=0.060, kelvin_offset=0.300,
                         maxh_coil=0.005, maxh_air=0.015,
                         order=2):
    """Build 3D IH mesh via OCC with Periodic Kelvin (no Cubit required).

    Creates gapped torus coil + workpiece cylinder + air sphere
    + exterior Kelvin sphere (periodic identification).
    Source/sink faces on coil gap.

    The 2-sphere Periodic Kelvin approach:
      - Interior sphere at origin (radius a_kelvin): physical domain
      - Exterior sphere at (kelvin_offset, 0, 0) (same radius): mapped domain
      - Periodic BC on sphere surfaces (Identify)
      - GND vertex at exterior sphere center (Dirichlet = 0 = infinity)
      - nu_kelvin = nu0 * (a/r')^2, r' from exterior center

    Args:
        R_coil, a_coil: Coil major/minor radius [m]
        gap_deg: Gap angle [degrees]
        R_wp, H_wp: Workpiece radius/height [m]
        a_kelvin: Kelvin boundary sphere radius [m]
        kelvin_offset: Exterior sphere center offset [m]
        maxh_coil, maxh_air: Mesh sizes [m]
        order: Mesh curve order

    Returns:
        (ngsolve.Mesh, dict) where dict contains:
            'a_kelvin': Kelvin radius
            'kelvin_center': (offset, 0, 0) exterior sphere center
            'R_coil', 'a_coil', 'R_wp', 'H_wp': geometry params
    """
    from ngsolve import Mesh, TaskManager
    from netgen.occ import (WorkPlane, Axes, Axis, Pnt, Dir, Sphere,
                             Cylinder, OCCGeometry, Glue, Vertex,
                             IdentificationType)

    # Gapped torus coil
    wp_plane = WorkPlane(Axes(p=Pnt(R_coil, 0, 0),
                              n=Dir(0, 1, 0), h=Dir(0, 0, 1)))
    circle = wp_plane.Circle(a_coil).Face()
    torus = circle.Revolve(Axis(Pnt(0, 0, 0), Dir(0, 0, 1)), 360 - gap_deg)
    torus.name = "coil"
    torus.maxh = maxh_coil

    # Label gap faces as source/sink
    expected_gap_area = math.pi * a_coil**2
    for f in torus.faces:
        area = f.mass
        if abs(area - expected_gap_area) / expected_gap_area < 0.3:
            angle = math.atan2(f.center.y, f.center.x)
            if abs(angle) < math.radians(gap_deg):
                f.name = "source"
            else:
                f.name = "sink"

    # Workpiece cylinder (hole: subtracted from air, NOT meshed)
    wp_cyl = Cylinder(Pnt(0, 0, -H_wp / 2), Dir(0, 0, 1), R_wp, H_wp)

    # Interior sphere (physical domain)
    inner_sphere = Sphere(Pnt(0, 0, 0), a_kelvin)

    # Boolean: air = inner_sphere - coil - wp (wp is a hole)
    air = inner_sphere - torus - wp_cyl
    air.name = "air"

    # Label cavity faces (from wp subtraction) as "sibc"
    # These are the faces whose area matches the workpiece surface
    wp_lateral_area = 2 * math.pi * R_wp * H_wp
    wp_cap_area = math.pi * R_wp**2
    for f in air.faces:
        area = f.mass
        if (abs(area - wp_lateral_area) / wp_lateral_area < 0.3 or
                abs(area - wp_cap_area) / wp_cap_area < 0.3):
            f.name = "sibc"
            f.maxh = min(R_wp / 3, maxh_coil)

    # Exterior sphere (Kelvin-mapped domain, same radius, offset center)
    ext_sphere = Sphere(Pnt(kelvin_offset, 0, 0), a_kelvin)
    ext_sphere.name = "kelvin"
    ext_sphere.maxh = maxh_air * 2

    # GND vertex at exterior sphere center (= physical infinity)
    gnd = Vertex(Pnt(kelvin_offset, 0, 0))
    gnd.name = "GND"

    # Periodic identification: match sphere surfaces
    # After boolean, air has multiple faces; the sphere surface is the one
    # with area close to 4*pi*a^2
    expected_sphere_area = 4 * math.pi * a_kelvin**2
    int_face = None
    for f in air.faces:
        if abs(f.mass - expected_sphere_area) / expected_sphere_area < 0.1:
            f.name = "kelvin_int"
            int_face = f
            break

    ext_face = None
    for f in ext_sphere.faces:
        if abs(f.mass - expected_sphere_area) / expected_sphere_area < 0.1:
            f.name = "kelvin_ext"
            ext_face = f
            break

    if int_face is not None and ext_face is not None:
        int_face.Identify(ext_face, "kelvin", IdentificationType.PERIODIC)

    # Hole approach: wp_cyl is NOT included (it's just a hole in air)
    shape = Glue([air, torus, ext_sphere, gnd])
    geo = OCCGeometry(shape)
    with TaskManager():
        ngmesh = geo.GenerateMesh(maxh=maxh_air, grading=0.3)
    mesh = Mesh(ngmesh)
    mesh.Curve(order)

    info = {
        'a_kelvin': a_kelvin,
        'kelvin_center': (kelvin_offset, 0, 0),
        'R_coil': R_coil, 'a_coil': a_coil,
        'R_wp': R_wp, 'H_wp': H_wp,
    }
    return mesh, info


def build_occ_ih_mesh_simple(R_coil=0.030, a_coil=0.003, gap_deg=5,
                             R_wp=0.010, H_wp=0.020,
                             R_air=0.080,
                             maxh_coil=0.005, maxh_air=0.015,
                             order=2):
    """Build 3D IH mesh via OCC WITHOUT Kelvin (truncated air sphere).

    For quick FEM testing / BEM-vs-FEM comparison. Dirichlet on outer sphere.
    Same coil+workpiece geometry as build_occ_ih_mesh_3d.

    Args:
        R_coil, a_coil: Coil major/minor radius [m]
        gap_deg: Gap angle [degrees]
        R_wp, H_wp: Workpiece radius/height [m]
        R_air: Outer air sphere radius [m] (Dirichlet truncation)
        maxh_coil, maxh_air: Mesh sizes [m]
        order: Mesh curve order

    Returns:
        (ngsolve.Mesh, dict)
    """
    from ngsolve import Mesh, TaskManager
    from netgen.occ import (WorkPlane, Axes, Axis, Pnt, Dir, Sphere,
                             Cylinder, OCCGeometry, Glue)

    # Gapped torus coil
    wp_plane = WorkPlane(Axes(p=Pnt(R_coil, 0, 0),
                              n=Dir(0, 1, 0), h=Dir(0, 0, 1)))
    circle = wp_plane.Circle(a_coil).Face()
    torus = circle.Revolve(Axis(Pnt(0, 0, 0), Dir(0, 0, 1)), 360 - gap_deg)
    torus.name = "coil"
    torus.maxh = maxh_coil

    # Label gap faces as source/sink
    expected_gap_area = math.pi * a_coil**2
    for f in torus.faces:
        area = f.mass
        if abs(area - expected_gap_area) / expected_gap_area < 0.3:
            angle = math.atan2(f.center.y, f.center.x)
            if abs(angle) < math.radians(gap_deg):
                f.name = "source"
            else:
                f.name = "sink"

    # Workpiece cylinder (hole: subtracted from air, NOT meshed)
    wp_cyl = Cylinder(Pnt(0, 0, -H_wp / 2), Dir(0, 0, 1), R_wp, H_wp)

    # Air sphere (truncated boundary)
    air_sphere = Sphere(Pnt(0, 0, 0), R_air)

    # Boolean: air = sphere - coil - workpiece (wp is a hole)
    air = air_sphere - torus - wp_cyl
    air.name = "air"
    # Label outer sphere as Dirichlet
    expected_sphere_area = 4 * math.pi * R_air**2
    wp_lateral_area = 2 * math.pi * R_wp * H_wp
    wp_cap_area = math.pi * R_wp**2
    for f in air.faces:
        area = f.mass
        if abs(area - expected_sphere_area) / expected_sphere_area < 0.1:
            f.name = "outer"
        elif (abs(area - wp_lateral_area) / wp_lateral_area < 0.3 or
              abs(area - wp_cap_area) / wp_cap_area < 0.3):
            f.name = "sibc"
            f.maxh = min(R_wp / 3, maxh_coil)

    # Hole approach: wp_cyl is NOT included
    shape = Glue([air, torus])
    geo = OCCGeometry(shape)
    with TaskManager():
        ngmesh = geo.GenerateMesh(maxh=maxh_air, grading=0.3)
    mesh = Mesh(ngmesh)
    mesh.Curve(order)

    info = {
        'R_air': R_air,
        'R_coil': R_coil, 'a_coil': a_coil,
        'R_wp': R_wp, 'H_wp': H_wp,
    }
    return mesh, info


def build_occ_ih_mesh_hole(R_coil=0.030, a_coil=0.003, gap_deg=5,
                           R_wp=0.010, H_wp=0.020,
                           R_air=0.080,
                           maxh_coil=0.005, maxh_air=0.015,
                           order=2):
    """Build 3D IH mesh with workpiece as HOLE (scattered-field FEM-SIBC).

    Workpiece volume is NOT meshed. wp_surface is an external boundary
    of the air domain. This is required for scattered-field formulation
    where H_t is extracted from A_total = A_inc + A_scat on BND.

    Faces are named BEFORE boolean to ensure OCC propagates labels.

    Returns:
        (ngsolve.Mesh, dict)
    """
    from ngsolve import Mesh, TaskManager
    from netgen.occ import (WorkPlane, Axes, Axis, Pnt, Dir, Sphere,
                             Cylinder, OCCGeometry, Glue)

    # Gapped torus coil
    wp_plane = WorkPlane(Axes(p=Pnt(R_coil, 0, 0),
                              n=Dir(0, 1, 0), h=Dir(0, 0, 1)))
    circle = wp_plane.Circle(a_coil).Face()
    torus = circle.Revolve(Axis(Pnt(0, 0, 0), Dir(0, 0, 1)), 360 - gap_deg)
    torus.name = "coil"
    torus.maxh = maxh_coil

    expected_gap_area = math.pi * a_coil**2
    for f in torus.faces:
        if abs(f.mass - expected_gap_area) / expected_gap_area < 0.3:
            angle = math.atan2(f.center.y, f.center.x)
            if abs(angle) < math.radians(gap_deg):
                f.name = "source"
            else:
                f.name = "sink"

    # Workpiece cylinder: name faces BEFORE boolean (OCC propagates names)
    wp_cyl = Cylinder(Pnt(0, 0, -H_wp / 2), Dir(0, 0, 1), R_wp, H_wp)
    for f in wp_cyl.faces:
        f.name = "wp_surface"
        f.maxh = min(R_wp / 3, maxh_coil)

    # Air sphere with hole (no workpiece volume)
    air = Sphere(Pnt(0, 0, 0), R_air) - torus - wp_cyl
    air.name = "air"

    expected_sphere_area = 4 * math.pi * R_air**2
    for f in air.faces:
        if f.name is None and abs(f.mass - expected_sphere_area) / expected_sphere_area < 0.1:
            f.name = "outer"

    shape = Glue([air, torus])  # NO wp_cyl (hole)
    geo = OCCGeometry(shape)
    with TaskManager():
        ngmesh = geo.GenerateMesh(maxh=maxh_air, grading=0.3)
    mesh = Mesh(ngmesh)
    mesh.Curve(order)

    info = {
        'R_air': R_air,
        'R_coil': R_coil, 'a_coil': a_coil,
        'R_wp': R_wp, 'H_wp': H_wp,
        'approach': 'hole',
    }
    return mesh, info


def _add_kelvin_point_idents_manual(mesh, kelvin_offset, idnr=1,
                                     point_tolerance=None):
    """Compatibility shim -- delegates to `kelvin_identify_ngsolve`.

    The full implementation lives in
    `radia.kelvin_identify_ngsolve.add_kelvin_identification`.  This
    function is preserved only because it appears as a private
    helper name in older docstrings / code review threads.  New code
    should use the public API.
    """
    from radia.kelvin_identify_ngsolve import add_kelvin_identification
    info = add_kelvin_identification(
        mesh, kelvin_offset=kelvin_offset, idnr=idnr,
        point_tolerance=point_tolerance, skip_if_existing=False)
    return info["n_pairs"]


def add_periodic_kelvin(mesh, kelvin_offset):
    """Ensure Kelvin Periodic identification exists on the mesh.

    Thin compatibility shim around the public
    `radia.kelvin_identify_ngsolve.add_kelvin_identification`.  Kept here
    so existing panel code (`calc_fem_kelvin.py`, `calc_accel_magnet.py`,
    ...) continues to work without import changes.  New callers should
    prefer the public API:

        from radia import add_kelvin_identification
        info = add_kelvin_identification(mesh)   # auto-detects offset

    Behaviour matches the pre-2026-05-13 implementation:
      - returns True iff some identification ended up on the mesh
      - skips silently when kelvin_int / kelvin_ext labels are missing
      - skips silently when ngmesh already has identifications
      - delegates the actual point-pair matching + AddPointIdentification
        calls to the public module
    """
    from radia.kelvin_identify_ngsolve import (
        add_kelvin_identification, has_kelvin_identification,
    )

    boundaries = mesh.GetBoundaries()
    if "kelvin_int" not in boundaries or "kelvin_ext" not in boundaries:
        return False

    if has_kelvin_identification(mesh):
        return True

    info = add_kelvin_identification(
        mesh, kelvin_offset=kelvin_offset, skip_if_existing=False)
    if info["n_pairs"] > 0:
        print(f"[add_periodic_kelvin/manual] {info['n_pairs']} point pairs "
              f"(max_dist={info['max_dist']:.2e}, "
              f"unmatched={info['n_unmatched']}, "
              f"tol={info['point_tolerance']:.2e})")
    return info["n_pairs"] > 0


def detect_outer_boundary(mesh):
    """Detect the outer boundary of the air domain (largest BND area).

    Works for full, 1/2, 1/4, 1/8 symmetry models: the outer spherical
    surface is ALWAYS the largest BND face because symmetry planes have
    holes (coil, workpiece subtracted).

    Returns:
        (name, area, R) where:
            name: boundary label string
            area: total area of the boundary [m^2]
            R: estimated sphere radius from max vertex distance [m]
        or (None, 0, 0) if no boundary found.
    """
    from ngsolve import Integrate, CF, BND
    bnd_names = sorted(set(mesh.GetBoundaries()))

    best_name, best_area = None, 0.0
    for name in bnd_names:
        if not name or name == "default":
            continue
        try:
            area = float(Integrate(CF(1), mesh, BND,
                                   definedon=mesh.Boundaries(name)).real)
        except Exception:
            continue
        if area > best_area:
            best_area = area
            best_name = name

    if best_name is None:
        return (None, 0.0, 0.0)

    # Estimate R from vertex positions on the outer boundary
    outer_verts = set()
    boundaries = mesh.GetBoundaries()
    for el in mesh.Elements(BND):
        if boundaries[el.index] == best_name:
            for v in el.vertices:
                outer_verts.add(v.nr)

    if outer_verts:
        coords = np.array([mesh.vertices[v].point for v in outer_verts])
        R = float(np.max(np.linalg.norm(coords, axis=1)))
    else:
        R = 0.0

    return (best_name, best_area, R)


def detect_symmetry(mesh, tol=1e-6):
    """Detect symmetry planes from mesh vertex positions.

    If ALL vertices have x >= -tol AND some vertices have x < tol
    (i.e., vertices lie ON the x=0 plane), then x=0 is a symmetry plane.

    Works for 1/2 (1 plane), 1/4 (2 planes), 1/8 (3 planes) models.

    Returns:
        list of plane names, e.g. ["z"], ["x", "z"], ["x", "y", "z"]
    """
    pts = np.array([v.point for v in mesh.vertices])
    sym = []
    for axis, name in enumerate(["x", "y", "z"]):
        all_positive = pts[:, axis].min() >= -tol
        has_on_plane = (pts[:, axis] < tol).any()
        if all_positive and has_on_plane:
            sym.append(name)
    return sym


def estimate_kelvin_maxh(mesh, outer_bnd_name):
    """Estimate mesh size for Kelvin domain from outer boundary edges.

    Returns the median edge length on the outer boundary, multiplied
    by 1.5 (Kelvin domain can be coarser than the physical domain).
    """
    from ngsolve import BND
    boundaries = mesh.GetBoundaries()
    outer_edges = set()
    for el in mesh.Elements(BND):
        if boundaries[el.index] == outer_bnd_name:
            for e in el.edges:
                outer_edges.add(e.nr)

    if not outer_edges:
        return 0.02  # fallback

    lengths = []
    for e_nr in outer_edges:
        edge = mesh.edges[e_nr]
        v0, v1 = edge.vertices
        p0 = np.array(mesh.vertices[v0.nr].point)
        p1 = np.array(mesh.vertices[v1.nr].point)
        lengths.append(float(np.linalg.norm(p1 - p0)))

    return float(np.median(lengths)) * 1.5


def detect_kelvin_offset(mesh):
    """Compatibility shim -- delegates to `kelvin_identify_ngsolve`.

    The full implementation lives in
    `radia.kelvin_identify_ngsolve.detect_kelvin_offset`.  Kept here for
    backward compatibility with existing panel callers
    (`calc_fem_kelvin.py`, `calc_fem_coilmesh.py`, `calc_kelvin_benchmark.py`,
    `calc_accel_magnet.py`).  New code should `from radia import
    detect_kelvin_offset` directly.
    """
    from radia.kelvin_identify_ngsolve import detect_kelvin_offset as _impl
    return _impl(mesh)


# ============================================================
# ESIM / SIBC Solver
# ============================================================

def get_bh_curve(material="steel", bh_file=None):
    """Get BH curve data for material.

    Delegates to :class:`EMMaterial`.  Kept for backward compatibility.
    """
    mat = EMMaterial.from_name(material, bh_file=bh_file or "")
    return mat.get_bh_curve()


def create_esim_solver(material="steel", frequency=7000, sigma=2e6,
                       half_thickness=0.010, geometry='cylinder',
                       bh_file=None):
    """Create ESIMFiniteSlabSolver with appropriate material settings.

    Delegates to :class:`EMMaterial`.  Kept for backward compatibility.
    """
    mat = EMMaterial.from_name(material, bh_file=bh_file or "")
    # Override sigma if the caller passed an explicit value (legacy API).
    if sigma != 2e6 or material not in ("steel", "elf_steel"):
        mat = EMMaterial(name=mat.name, sigma=sigma, mu_r=mat.mu_r,
                         bh_curve=mat.bh_curve, hys_file=mat.hys_file)
    return mat.create_esim_solver(frequency, half_thickness, geometry)


def sigma_for_material(material):
    """Default conductivity [S/m] for common materials.

    Delegates to :class:`EMMaterial`.  Kept for backward compatibility.
    """
    try:
        return EMMaterial.from_name(material).sigma
    except ValueError:
        return 2e6


def write_surface_only_vol(src_ngmesh, output_path, keep_material=""):
    """Extract surface elements from a volume mesh and save as surface-only .vol.

    BEM (HDivSurface) requires surface-only mesh. On a volume mesh,
    HDivSurface includes all interior edges as DOFs, making the SL
    matrix singular.

    Args:
        src_ngmesh: netgen.meshing.Mesh with volume + surface elements
        output_path: path to write surface-only .vol
        keep_material: if set, only keep surface elements adjacent to this
            material (by domin/domout in facedescriptors).  Empty string
            keeps all surfaces (original behaviour).
    """
    import tempfile

    tmp = os.path.join(tempfile.gettempdir(), "_vol_tmp.vol")
    src_ngmesh.Save(tmp)

    with open(tmp, "r") as f:
        lines = f.readlines()

    # --- When keep_material is specified, find the material index ---
    mat_idx = 0  # 0 = no filter
    if keep_material:
        mat_idx = _find_material_index(lines, keep_material)

    # Build output lines in memory (avoids Windows text-mode seek issues)
    out = []
    section = ""
    se_insert_idx = -1  # where to insert filtered count
    kept_elems = []

    for line in lines:
        s = line.strip()

        # --- skip volumeelements entirely ---
        if s == "volumeelements":
            section = "volumeelements"
            out.append("volumeelements\n")
            out.append("0\n")
            continue
        if section == "volumeelements":
            if _is_section_header(s):
                section = ""
            else:
                continue

        # --- filter surfaceelements / surfaceelementsuv by domin/domout ---
        if s in ("surfaceelements", "surfaceelementsuv"):
            section = "surfaceelements"
            out.append(line)
            if mat_idx > 0:
                se_insert_idx = len(out)  # count placeholder
                out.append("")  # will be replaced
            continue

        if section == "surfaceelements":
            if _is_section_header(s):
                section = ""
            elif s.startswith("#") or not s:
                continue
            elif mat_idx > 0:
                parts = s.split()
                # Format: surfid bcprop domin domout np p1 p2 p3 ...
                if len(parts) >= 5:
                    try:
                        domin = int(parts[2])
                        domout = int(parts[3])
                        if domin == mat_idx or domout == mat_idx:
                            kept_elems.append(line)
                    except (ValueError, IndexError):
                        pass
                continue
            # mat_idx == 0: pass through (count + elements)

        out.append(line)

    # Insert filtered surface elements
    if mat_idx > 0 and se_insert_idx >= 0:
        out[se_insert_idx] = f"{len(kept_elems)}\n"
        for i, el in enumerate(kept_elems):
            out.insert(se_insert_idx + 1 + i, el)

    with open(output_path, "w") as f:
        f.writelines(out)

    # When filtering by material, also drop FaceDescriptors and bcnames that
    # have no surface elements left. Otherwise NGSolve registers spurious
    # boundary regions and the BEM saddle-point matrix becomes singular.
    if mat_idx > 0:
        _strip_unused_face_descriptors(output_path)

    os.remove(tmp)


def _strip_unused_face_descriptors(vol_path):
    """Remove FaceDescriptors and bcnames that have zero surface elements.

    Required after `write_surface_only_vol` filters by material: the
    surface element list is reduced but the FD/bcname tables are not.
    NGSolve treats orphan FDs as separate boundary regions and the BEM
    saddle-point matrix becomes rank-deficient.

    Strategy: parse the .vol into sections, find which 1-based FD indices
    still appear in surface elements, build an old->new index map, then
    rewrite facedescriptors / surfaceelements / bcnames consistently.
    """
    with open(vol_path, "r") as f:
        raw = f.read()

    # Split into "section blocks": each block starts at a section header
    # line (a single keyword) and runs until the next section header.
    section_keywords = {
        "mesh3d", "dimension", "geomtype", "facedescriptors",
        "surfaceelements", "surfaceelementsuv", "surfaceelementsgi",
        "volumeelements", "edgesegmentsgi2", "points", "pointelements",
        "identifications", "identificationtypes", "identificationnames",
        "materials", "bcnames", "cd2names", "cd3names",
        "singular_points", "singular_lines", "singular_faces",
        "curvedelements", "endmesh",
    }

    lines = raw.splitlines(keepends=True)
    # Identify section boundaries.
    section_starts = []  # list of (line_index, keyword)
    for i, line in enumerate(lines):
        s = line.strip()
        if s in section_keywords:
            section_starts.append((i, s))
    section_starts.append((len(lines), "EOF"))

    # Build a dict: keyword -> list of lines AFTER the keyword and count
    # line, until the next section.
    blocks = {}
    for j in range(len(section_starts) - 1):
        i0, kw = section_starts[j]
        i1, _ = section_starts[j + 1]
        # Skip the keyword line itself
        body = lines[i0 + 1 : i1]
        blocks[kw] = (i0, body)

    # --- Pass 1: find used FD indices in surfaceelements ---
    # Surface element format: surfid bcprop domin domout np p1 p2 p3 ...
    # surfid is the original Cubit surface ID; bcprop is the contiguous
    # 1-based FD index that NGSolve uses to look up the FaceDescriptor.
    # We must collect bcprop values (parts[1]), NOT surfid (parts[0]).
    used_fds = set()
    for kw in ("surfaceelements", "surfaceelementsuv", "surfaceelementsgi"):
        if kw not in blocks:
            continue
        body = blocks[kw][1]
        count_skipped = False
        for ln in body:
            s = ln.strip()
            if not s or s.startswith("#"):
                continue
            if not count_skipped:
                count_skipped = True
                continue
            parts = s.split()
            if len(parts) >= 5:
                try:
                    used_fds.add(int(parts[1]))  # bcprop = 1-based FD index
                except ValueError:
                    pass

    if not used_fds:
        return

    # --- Pass 2: parse facedescriptors to get total nfd ---
    if "facedescriptors" not in blocks:
        return
    fd_body = blocks["facedescriptors"][1]
    fd_data = []
    count_line = None
    for ln in fd_body:
        s = ln.strip()
        if not s or s.startswith("#"):
            continue
        if count_line is None:
            count_line = ln
            continue
        fd_data.append(ln)

    nfd = int(count_line.strip())
    if len(fd_data) != nfd:
        # Defensive: if parsing failed, leave file untouched.
        return

    # Build old (1-based) -> new (1-based) FD map.
    old_to_new = {}
    new_idx = 1
    for k in range(1, nfd + 1):
        if k in used_fds:
            old_to_new[k] = new_idx
            new_idx += 1
    if len(old_to_new) == nfd:
        return  # nothing to strip

    # --- Pass 3: rewrite the file ---
    out = []
    skip_until_next_section = False
    j = 0
    while j < len(lines):
        line = lines[j]
        s = line.strip()
        if s == "facedescriptors":
            out.append(line)
            new_fds = []
            for k, fl in enumerate(fd_data, 1):
                if k not in old_to_new:
                    continue
                parts = fl.split()
                if len(parts) >= 5:
                    # facedescriptors: surfnr domin domout tlosurf bcprop
                    #   surfnr / tlosurf = original geometric surface id
                    #     (Cubit surface ID) — KEEP as-is
                    #   bcprop = the 1-based FD index used by surface elements
                    #     to look up this FD; renumber consecutively.
                    parts[4] = str(old_to_new[k])
                new_fds.append(" " + " ".join(parts) + "\n")
            out.append(f"{len(new_fds)}\n")
            out.extend(new_fds)
            # Skip the original FD body in the input lines
            j += 1
            # Skip count line + nfd data lines (preserving comments)
            consumed = 0
            target = 1 + nfd  # count + entries
            while j < len(lines) and consumed < target:
                t = lines[j].strip()
                if t and not t.startswith("#"):
                    consumed += 1
                j += 1
            continue

        if s in ("surfaceelements", "surfaceelementsuv", "surfaceelementsgi"):
            out.append(line)
            j += 1
            # Read count + nse element lines from input
            se_body_lines = []
            count_consumed = False
            nse = 0
            while j < len(lines):
                t = lines[j].strip()
                if t in section_keywords:
                    break
                if not t or t.startswith("#"):
                    se_body_lines.append(lines[j])
                    j += 1
                    continue
                if not count_consumed:
                    nse = int(t)
                    se_body_lines.append(lines[j])
                    count_consumed = True
                    j += 1
                    continue
                parts = lines[j].split()
                if len(parts) >= 5:
                    try:
                        # Update bcprop (parts[1]) only; surfid (parts[0])
                        # is the original Cubit surface ID and should be
                        # preserved.
                        old_bc = int(parts[1])
                        if old_bc in old_to_new:
                            parts[1] = str(old_to_new[old_bc])
                    except ValueError:
                        pass
                se_body_lines.append(" " + " ".join(parts) + "\n")
                j += 1
            out.extend(se_body_lines)
            continue

        if s == "bcnames":
            out.append(line)
            j += 1
            # Read count + nbc lines
            count_consumed = False
            nbc = 0
            bc_data = []
            while j < len(lines):
                t = lines[j].strip()
                if t in section_keywords:
                    break
                if not t:
                    j += 1
                    continue
                if not count_consumed:
                    nbc = int(t)
                    count_consumed = True
                    j += 1
                    continue
                bc_data.append(lines[j])
                j += 1
            new_bcs = []
            for k, bl in enumerate(bc_data, 1):
                if k not in old_to_new:
                    continue
                parts = bl.split(None, 1)
                rest = parts[1].rstrip("\r\n") if len(parts) > 1 else "default"
                new_bcs.append(f"{old_to_new[k]}\t{rest}\n")
            out.append(f"{len(new_bcs)}\n")
            out.extend(new_bcs)
            continue

        out.append(line)
        j += 1

    with open(vol_path, "w") as f:
        f.writelines(out)


def _find_material_index(vol_lines, name):
    """Return 1-based material index for *name* in a .vol file, or 0."""
    in_mat = False
    for line in vol_lines:
        s = line.strip()
        if s == "materials":
            in_mat = True
            continue
        if in_mat:
            if _is_section_header(s):
                break
            parts = s.split(None, 1)
            if len(parts) == 2 and parts[1] == name:
                return int(parts[0])
    return 0


def _find_adjacent_fds(vol_lines, mat_idx):
    """Return set of 1-based face descriptor indices adjacent to *mat_idx*."""
    in_fd = False
    fd_set = set()
    idx = 0
    for line in vol_lines:
        s = line.strip()
        if s == "facedescriptors":
            in_fd = True
            continue
        if in_fd:
            if _is_section_header(s):
                break
            parts = s.split()
            if len(parts) >= 4 and parts[0][0].isdigit():
                idx += 1
                # facedescriptors: surfnr domin domout tlosurf bcprop
                domin = int(parts[1])
                domout = int(parts[2])
                if domin == mat_idx or domout == mat_idx:
                    fd_set.add(idx)
    return fd_set


def _is_section_header(s):
    """True if *s* looks like a .vol section keyword (not a data line)."""
    if not s or s.startswith("#"):
        return False
    return (s.isalpha() or s in (
        "edgesegmentsgi2", "facedescriptors", "face_transparencies",
        "surfaceelements", "surfaceelementsuv",
        "volumeelements", "curvedelements",
        "points", "materials", "bcnames", "endmesh",
    ))


# ============================================================
# Subprocess Output Protocol
# ============================================================

def progress(tag, msg):
    """Write progress message to stderr AND the panel debug log file.

    Protocol: stderr lines matching ``TAG:detail`` are parsed by
    register_toolbar.py ``_on_stderr`` handlers.  Tags should be
    short uppercase identifiers.

    History (2026-05-23): historically this only wrote to stderr.
    Keiko's 27-minute BEM-A + BEM run on mdx finished cleanly but the
    panel log file had ZERO intermediate progress lines between start
    and end.  For long runs the user had no post-hoc record of WHERE
    the time was spent.  Mirror every ``progress()`` call to
    ``panel_log`` so successful long runs leave the same timeline in
    the file that the notebook/workbench shows live.  panel_log gates
    the file write through its own
    handle (only opens once per process), so the cost is one
    formatted write per progress line -- negligible vs the work
    those lines are summarising.

    Examples::

        progress("MESH_READY", "42 elements exported")
        progress("BEM-SIBC", "wp R=10.0mm, sigma=2e6")
        progress("FEM", "354 DOFs, solving...")
    """
    sys.stderr.write(f"{tag}:{msg}\n")
    sys.stderr.flush()
    _progress_to_panel_log(f"{tag}:{msg}")


# Cached panel_log function -- resolved lazily on first progress() call
# so we don't pay the import cost on every line, and so calc_*.py
# modules that import this file BEFORE calc_main() runs still work
# (panel_log's init_panel_log is called from calc_main).  When
# panel_log is unavailable (e.g. pytest direct import without going
# through calc_main), the no-op shim keeps progress() side-effect-free.
_progress_panel_log_fn = None


def _progress_to_panel_log(line):
    """Mirror a single progress line to the shared panel debug log."""
    global _progress_panel_log_fn
    if _progress_panel_log_fn is None:
        try:
            from panel_log import panel_log as _pl
            _progress_panel_log_fn = _pl
        except Exception:
            # No-op shim: panel_log unreachable from this subprocess.
            _progress_panel_log_fn = lambda _msg: None  # noqa: E731
    try:
        _progress_panel_log_fn(line)
    except Exception:
        # Log writes must never break the actual computation.
        pass


def calc_main(solve_func, parser):
    """Common main() wrapper for all calc_*.py scripts.

    Handles:
      - stdout suppression (NGSolve/Cubit banners)
      - Exception -> ``{"error": ...}`` JSON
      - ``--output FILE`` JSON dump + stdout JSON print
      - Traceback on stderr for debugging

    Usage in calc_*.py::

        def main():
            parser = argparse.ArgumentParser(...)
            parser.add_argument(...)

            def run(args):
                return solve_something(args.cub5, args.order, ...)

            calc_main(run, parser)

    Args:
        solve_func: callable(args) -> dict.  Receives parsed argparse
            Namespace; must return a JSON-serializable dict.
        parser: argparse.ArgumentParser (must include ``--output``
            argument if file output is desired; added automatically
            if missing).
    """
    # Ensure --output is available
    known_dests = {a.dest for a in parser._actions}
    if "output" not in known_dests:
        parser.add_argument("--output", default="",
                            help="JSON output file (optional)")

    args = parser.parse_args()

    # Tag this subprocess in the shared panel debug log so its lines
    # are interleaved with Cubit / IH window lines under one Cubit
    # session.
    try:
        from panel_log import init_panel_log, panel_log
        # Use the entry-point script name (e.g. calc_inductance.py)
        # as the component tag, truncated to 12 chars.
        _entry = os.path.basename(sys.argv[0])
        _tag = os.path.splitext(_entry)[0].replace("calc_", "")
        init_panel_log(_tag, truncate=False, banner=False)
        panel_log(f"calc_main: {_entry} args={vars(args)}")
    except Exception:
        def panel_log(_msg):
            pass

    import time as _time
    _t0 = _time.perf_counter()

    real_stdout = sys.stdout
    sys.stdout = io.StringIO()
    exception_caught = False
    try:
        result = solve_func(args)
    except Exception as e:
        result = {"error": str(e)}
        exception_caught = True
        sys.stderr.write(traceback.format_exc())
        sys.stderr.flush()
        try:
            panel_log(f"calc_main: EXCEPTION {e}")
            for line in traceback.format_exc().rstrip().splitlines():
                panel_log(f"  {line}")
        except Exception:
            pass
    sys.stdout = real_stdout

    # Stamp the radia version that actually ran this calc into the
    # result dict (auditability: the persisted <base>_<suffix>.json
    # next to the geometry records which radia produced the numbers).
    # Prefer the already-imported module's __version__ (authoritative,
    # free if solve_func imported radia); fall back to dist metadata.
    if isinstance(result, dict) and "radia_version" not in result:
        _rv = None
        try:
            _m = sys.modules.get("radia")
            _rv = getattr(_m, "__version__", None) if _m is not None else None
            if not _rv:
                from importlib.metadata import version as _pkgver
                _rv = _pkgver("radia")
        except Exception:
            _rv = None
        if _rv:
            result["radia_version"] = _rv

    _elapsed = _time.perf_counter() - _t0
    try:
        if isinstance(result, dict):
            keys = sorted(result.keys())[:15]
            panel_log(f"calc_main: done {_elapsed:.2f}s keys={keys}")
            if "error" in result:
                panel_log(f"  ERROR: {result['error']}")
    except Exception:
        pass

    output_path = getattr(args, "output", "")
    if output_path:
        with open(output_path, "w", encoding="utf-8") as f:
            json.dump(result, f)
    print(json.dumps(result))

    # CLAUDE.md "No Fallbacks -- Fail Fast, Fail Loud": a calc_*.py run
    # that caught an exception OR returned {"error": ...} / {"status":
    # "error"} MUST exit non-zero so subprocess callers (sweep scripts,
    # panel UI, CI) see the failure unambiguously.  The error JSON is
    # still written for audit (panel Output window + sweep_results.json
    # error entry) but exit code reflects success / failure.
    if exception_caught or (
        isinstance(result, dict)
        and ("error" in result or result.get("status") == "error")
    ):
        sys.exit(1)
