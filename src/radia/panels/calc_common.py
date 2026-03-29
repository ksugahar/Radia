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

STEEL_BH = [
    [0, 0], [50, 0.1], [100, 0.25], [200, 0.55],
    [500, 0.95], [1000, 1.2], [2000, 1.4], [5000, 1.55],
    [10000, 1.65], [20000, 1.75], [50000, 1.9], [100000, 2.0],
]


def setup_paths():
    """Add radia src to sys.path (for esim_cell_problem, etc.)."""
    radia_src = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
    abs_src = os.path.abspath(radia_src)
    if abs_src not in sys.path:
        sys.path.insert(0, abs_src)


def setup_cubit(cub5_path=None):
    """Initialize Cubit and optionally open a .cub5 file.

    Handles path cleanup to avoid scipy/numpy DLL conflicts.

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
    """Export Cubit mesh to NGSolve via export_NGSolveCurvedMesh.

    Args:
        cubit_mod: Cubit module (already initialized with model open)
        order: Curve order (1-3)
        surface_only: True for BEM surface mesh
        split_quads: Split quads into tris (default: True for surface_only)

    Returns:
        ngsolve.Mesh
    """
    setup_paths()
    import cubit_mesh_export
    if split_quads is None:
        split_quads = surface_only
    return cubit_mesh_export.export_NGSolveCurvedMesh(
        cubit_mod, order=order, surface_only=surface_only,
        split_quads=split_quads)


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

    if material == "steel":
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
