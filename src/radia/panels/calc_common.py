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

# BH curve (100 points, from ELF_MAGIC CEFC 2020 nonlinear model)
# Source: examples/c_type_electromagnet/nonlinear/BH.txt
# Units: H [A/m], B [T]
STEEL_BH = [
    [0.0, 0.0],
    [13.898, 0.22296], [15.397, 0.25304], [17.058, 0.28380],
    [18.898, 0.31552], [20.936, 0.34852], [23.194, 0.38323],
    [25.696, 0.42011], [28.467, 0.45974], [31.538, 0.50272],
    [34.939, 0.54965], [38.708, 0.60110], [42.883, 0.65744],
    [47.508, 0.71868], [52.632, 0.78437], [58.309, 0.85340],
    [64.598, 0.92403], [71.565, 0.99400], [79.284, 1.06090],
    [87.836, 1.12254], [97.310, 1.17738], [107.806, 1.22465],
    [119.433, 1.26440], [132.315, 1.29727], [146.587, 1.32427],
    [162.397, 1.34654], [179.913, 1.36518], [199.319, 1.38116],
    [220.817, 1.39530], [244.634, 1.40821], [271.020, 1.42039],
    [300.252, 1.43217], [332.636, 1.44381], [368.514, 1.45547],
    [408.262, 1.46728], [452.296, 1.47930], [501.081, 1.49157],
    [555.127, 1.50410], [615.002, 1.51691], [681.335, 1.52999],
    [754.823, 1.54332], [836.238, 1.55689], [926.433, 1.57068],
    [1026.357, 1.58467], [1137.059, 1.59883], [1259.701, 1.61315],
    [1395.571, 1.62761], [1546.096, 1.64220], [1712.856, 1.65688],
    [1897.603, 1.67166], [2102.276, 1.68651], [2329.025, 1.70142],
    [2580.231, 1.71638], [2858.532, 1.73137], [3166.850, 1.74640],
    [3508.423, 1.76144], [3886.837, 1.77649], [4306.067, 1.79154],
    [4770.514, 1.80658], [5285.057, 1.82162], [5855.097, 1.83664],
    [6486.621, 1.85165], [7186.261, 1.86663], [7961.362, 1.88158],
    [8820.066, 1.89651], [9771.388, 1.91142], [10825.319, 1.92629],
    [11992.926, 1.94114], [13286.469, 1.95596], [14719.532, 1.97077],
    [16307.164, 1.98555], [18066.037, 2.00033], [20014.619, 2.01510],
    [22173.373, 2.02987], [24564.968, 2.04466], [27214.517, 2.05948],
    [30149.844, 2.07434], [33401.772, 2.08927], [37004.449, 2.10428],
    [40995.707, 2.11941], [45417.458, 2.13468], [50316.133, 2.15014],
    [55743.174, 2.16583], [61755.569, 2.18181], [68416.455, 2.19815],
    [75795.776, 2.21493], [83971.022, 2.23225], [93028.041, 2.25023],
    [103061.940, 2.26901], [114178.084, 2.28877], [126493.203, 2.30972],
    [140136.616, 2.33211], [155251.592, 2.35623], [171996.852, 2.38247],
    [190548.237, 2.41122], [211100.554, 2.44300], [233869.620, 2.47834],
    [259094.531, 2.51786], [287040.173, 2.56214], [318000.0, 2.61173],
]


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
            plugin commands (radia_export gmsh/nastran/vtk) are available.

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
    """Export Cubit mesh to NGSolve via C++ radia_export netgen command.

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
    cubit_mod.cmd(f'radia_export netgen "{vol_path}" order {order} overwrite')
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

    # Workpiece cylinder (meshed as air, SIBC on interface)
    wp_cyl = Cylinder(Pnt(0, 0, -H_wp / 2), Dir(0, 0, 1), R_wp, H_wp)
    wp_cyl.name = "workpiece"
    for f in wp_cyl.faces:
        f.name = "wp_surface"
        f.maxh = min(R_wp / 3, maxh_coil)

    # Interior sphere (physical domain)
    inner_sphere = Sphere(Pnt(0, 0, 0), a_kelvin)

    # Boolean: air = inner_sphere - coil - wp
    air = inner_sphere - torus - wp_cyl
    air.name = "air"

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

    shape = Glue([air, wp_cyl, torus, ext_sphere, gnd])
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

    # Workpiece cylinder
    wp_cyl = Cylinder(Pnt(0, 0, -H_wp / 2), Dir(0, 0, 1), R_wp, H_wp)
    wp_cyl.name = "workpiece"
    for f in wp_cyl.faces:
        f.name = "wp_surface"
        f.maxh = min(R_wp / 3, maxh_coil)

    # Air sphere (truncated boundary)
    air_sphere = Sphere(Pnt(0, 0, 0), R_air)

    # Boolean: air = sphere - coil - workpiece
    air = air_sphere - torus - wp_cyl
    air.name = "air"
    # Label outer sphere as Dirichlet
    expected_sphere_area = 4 * math.pi * R_air**2
    for f in air.faces:
        if abs(f.mass - expected_sphere_area) / expected_sphere_area < 0.1:
            f.name = "outer"
            break

    shape = Glue([air, wp_cyl, torus])
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


def add_periodic_kelvin(mesh, kelvin_offset):
    """Add periodic identification to Cubit-exported mesh for Kelvin transform.

    Matches vertices on 'kelvin_int' boundary to corresponding points
    on 'kelvin_ext' (at position + kelvin_offset) using translation mapping.
    Modifies the ngmesh in place.

    For Cubit meshes, the sphere must be webcut first to create curves,
    then `copy mesh surface` ensures matching node topology.

    Args:
        mesh: ngsolve.Mesh (from export_NGSolveCurvedMesh)
        kelvin_offset: (x, y, z) offset of exterior sphere center

    Returns:
        True if identification was added, False otherwise
    """
    boundaries = mesh.GetBoundaries()
    if "kelvin_int" not in boundaries and "kelvin_ext" not in boundaries:
        return False

    from netgen.meshing import Trafo, Vec3d
    trafo = Trafo(Vec3d(*kelvin_offset))

    ngmesh = mesh.ngmesh
    identnr = ngmesh.IdentifyPeriodicBoundaries(
        "kelvin", "kelvin_int", trafo, point_tolerance=1e-3)

    return identnr > 0


def detect_kelvin_offset(mesh):
    """Detect Kelvin sphere center offset from mesh vertex positions.

    Finds the centroid of 'kelvin' material region as the offset.

    Returns:
        (x, y, z) offset or (0.3, 0, 0) default
    """
    from ngsolve import VOL
    try:
        materials = mesh.GetMaterials()
        if "kelvin" not in materials:
            return (0.3, 0, 0)

        # Average vertex positions in kelvin region
        kelvin_verts = set()
        for el in mesh.Elements(VOL):
            if materials[el.mat] == "kelvin" if hasattr(el, 'mat') else False:
                for v in el.vertices:
                    kelvin_verts.add(v.nr)

        if not kelvin_verts:
            # Fallback: scan all elements
            for i, el in enumerate(mesh.Elements(VOL)):
                mat_idx = el.index if hasattr(el, 'index') else i
                if mat_idx < len(materials) and "kelvin" in materials[mat_idx].lower():
                    for v in el.vertices:
                        kelvin_verts.add(v.nr)

        if kelvin_verts:
            import numpy as np
            coords = np.array([mesh.vertices[v].point for v in kelvin_verts])
            center = coords.mean(axis=0)
            return tuple(center)
    except Exception:
        pass

    return (0.3, 0, 0)


# ============================================================
# ESIM / SIBC Solver
# ============================================================

def get_bh_curve(material="steel", bh_file=None):
    """Get BH curve data for material.

    Args:
        material: "steel", "copper", or "aluminum"
        bh_file: Optional path to custom BH file (2 columns: H B)

    Returns:
        (bh_curve, mu_r): bh_curve is list or None, mu_r is float or None
    """
    if bh_file and bh_file != "(built-in Steel)" and os.path.exists(bh_file):
        data = np.loadtxt(bh_file)
        if data.ndim == 2 and data.shape[1] >= 2:
            return data[:, :2].tolist(), None
        raise ValueError(f"Invalid BH file format: {bh_file}")

    if material in ("steel", "elf_steel"):
        return STEEL_BH, None
    else:
        return None, 1.0  # non-magnetic


def create_esim_solver(material="steel", frequency=7000, sigma=2e6,
                       half_thickness=0.010, geometry='cylinder',
                       bh_file=None):
    """Create ESIMFiniteSlabSolver with appropriate material settings.

    Args:
        material: "steel", "copper", or "aluminum"
        frequency: Operating frequency [Hz]
        sigma: Conductivity [S/m]
        half_thickness: Workpiece half-thickness or radius [m]
        geometry: 'cylinder' (Bessel) or 'slab' (cosh/sinh)
        bh_file: Optional path to custom BH file

    Returns:
        ESIMFiniteSlabSolver instance
    """
    setup_paths()
    from esim_cell_problem import ESIMFiniteSlabSolver

    bh_curve, mu_r = get_bh_curve(material, bh_file)

    return ESIMFiniteSlabSolver(
        half_thickness=half_thickness,
        bh_curve=bh_curve,
        sigma=sigma,
        frequency=frequency,
        mu_r=mu_r if bh_curve is None else None,
        n_nodes=200,
        geometry=geometry)


def sigma_for_material(material):
    """Default conductivity [S/m] for common materials."""
    return {'steel': 2e6, 'copper': 5.8e7, 'aluminum': 3.5e7}.get(material, 2e6)


def write_surface_only_vol(src_ngmesh, output_path):
    """Extract surface elements from a volume mesh and save as surface-only .vol.

    BEM (HDivSurface) requires surface-only mesh. On a volume mesh,
    HDivSurface includes all interior edges as DOFs, making the SL
    matrix singular.

    Args:
        src_ngmesh: netgen.meshing.Mesh with volume + surface elements
        output_path: path to write surface-only .vol
    """
    import tempfile

    tmp = os.path.join(tempfile.gettempdir(), "_vol_tmp.vol")
    src_ngmesh.Save(tmp)

    with open(tmp, "r") as f:
        lines = f.readlines()

    with open(output_path, "w") as f:
        skip = False
        for line in lines:
            s = line.strip()
            if s == "volumeelements":
                skip = True
                f.write("volumeelements\n0\n")
                continue
            if skip:
                if (s and not s[0].isdigit() and s != "edgesegmentsgi2"
                        and not s.startswith("#")):
                    skip = False
                elif s in ("edgesegmentsgi2", "points", "surfaceelements",
                           "materials", "bcnames", "curvedelements",
                           "endmesh", "face_transparencies",
                           "facedescriptors"):
                    skip = False
                else:
                    continue
            f.write(line)

    os.remove(tmp)


# ============================================================
# T0 Source/Sink Current Injection
# ============================================================

def compute_T0_source(mesh, order, I_total=1.0):
    """Compute current density J via T0 technique (source/sink faces).

    Solves Laplace equation in coil volume:
        -div(grad(T0)) = 0 in "coil"
        T0 = 1 on "source", T0 = 0 on "sink"
        dT0/dn = 0 on lateral surface (natural BC)

    Returns:
        CoefficientFunction J = I_total / Phi * grad(T0)
        where Phi = total flux through source face
    """
    from ngsolve import (H1, BilinearForm, LinearForm, GridFunction,
                         Integrate, grad, dx, ds, CF, BND)

    fes_T0 = H1(mesh, order=order, dirichlet="source|sink")
    u, v = fes_T0.TnT()

    a = BilinearForm(fes_T0, check_unused=False)
    a += grad(u) * grad(v) * dx("coil")
    a.Assemble()

    gf_T0 = GridFunction(fes_T0)
    gf_T0.Set(CF(1), definedon=mesh.Boundaries("source"))

    # Solve with Dirichlet lift
    r = gf_T0.vec.CreateVector()
    r.data = -a.mat * gf_T0.vec
    gf_T0.vec.data += a.mat.Inverse(fes_T0.FreeDofs(), inverse="pardiso") * r

    # Normalize: Phi = integral of grad(T0) . n over source face
    # By Gauss's theorem, Phi = integral of |grad(T0)|^2 in coil
    # (since div(grad T0)=0 and T0 drops from 1 to 0)
    Phi = Integrate(grad(gf_T0) * grad(gf_T0), mesh,
                    definedon=mesh.Materials("coil"))
    if abs(Phi) < 1e-30:
        raise RuntimeError("T0 flux is zero - check source/sink blocks")

    J_source = (I_total / Phi) * grad(gf_T0)
    return J_source, gf_T0


def compute_J_theta(I_total, a_coil):
    """Compute J = J0 * e_theta for axisymmetric torus coil (fallback).

    Returns:
        CoefficientFunction J = J0 * (-y/r, x/r, 0)
    """
    from ngsolve import CF, x, y, sqrt, IfPos

    J0 = I_total / (math.pi * a_coil**2)
    r_xy = sqrt(x * x + y * y)
    r_safe = IfPos(r_xy - 1e-10, r_xy, 1e-10)
    return J0 * CF((-y / r_safe, x / r_safe, 0))


# ============================================================
# Subprocess Output Protocol
# ============================================================

def progress(tag, msg):
    """Write progress message to stderr for panel GUI consumption.

    Protocol: stderr lines matching ``TAG:detail`` are parsed by
    register_toolbar.py ``_on_stderr`` handlers.  Tags should be
    short uppercase identifiers.

    Examples::

        progress("MESH_READY", "42 elements exported")
        progress("BEM-SIBC", "wp R=10.0mm, sigma=2e6")
        progress("FEM", "354 DOFs, solving...")
    """
    sys.stderr.write(f"{tag}:{msg}\n")
    sys.stderr.flush()


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

    real_stdout = sys.stdout
    sys.stdout = io.StringIO()
    try:
        result = solve_func(args)
    except Exception as e:
        result = {"error": str(e)}
        sys.stderr.write(traceback.format_exc())
        sys.stderr.flush()
    sys.stdout = real_stdout

    output_path = getattr(args, "output", "")
    if output_path:
        with open(output_path, "w", encoding="utf-8") as f:
            json.dump(result, f)
    print(json.dumps(result))
