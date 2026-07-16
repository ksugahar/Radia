"""
Accelerator magnet HDiv-VIM solver for the electromagnet panel.

Pipeline:
  Netgen .vol mesh -> NGSolve read -> radia.vim.MeshSoftIron -> rad.Solve

User workflow:
  1. Write coil Python script (defines build_coil() -> CoilBuilder)
  2. Write .jou to create iron yoke hex mesh + named blocks, export .vol
     via `export netgen`
  3. Run this script with --vol <file>

Mesh materials (user sets in .jou before export):
  yoke  - volume material, hex/tet/wedge elements (iron)

Symmetry (auto-detected from `sym_<bc>_<axis>` sidesets in the .vol):

  sym_bn=0_<axis> -- B . n = 0 on axis=0 plane (field parallel)
                     -> Radia image sign '+' (mirror-symmetric)
  sym_ht=0_<axis> -- H x n = 0 on axis=0 plane (field perpendicular)
                     -> Radia image sign '-' (mirror-antisymmetric)

E.g. a 1/4 xz model with sym_ht=0_x + sym_bn=0_z auto-produces the
IMA string '-x+z'.  Pass --ima <string> to override.  Pre-2026-04-25
samples that omit the sym_*_* sidesets work unchanged -- auto-detect
yields '' which is "no IMA".

Coil is NOT meshed -- Biot-Savart analytical source via CoilBuilder.to_radia().

Mesh interface:
  This script takes a Netgen .vol file via --vol (not .cub5).  The .vol
  is the sole interface between Cubit and NGSolve per the Cubit/NGSolve
  Complete Separation Policy in CLAUDE.md -- this script does NOT
  import cubit.  Element types (hex/tet/wedge) are extracted from the
  NGSolve mesh directly.
"""

import argparse
import importlib.util
import math
import os
import sys
import time

import numpy as np

# Shared utilities
_this_dir = os.path.dirname(os.path.abspath(__file__))
if _this_dir not in sys.path:
    sys.path.insert(0, _this_dir)

from calc_common import (MU_0, setup_paths,
                          get_bh_curve, progress, calc_main,
                          EMMaterial, add_material_args)


def _log(msg):
    """Write progress to stderr (panel reads these)."""
    progress("HDiv", msg)


def ima_from_mesh_labels(boundaries):
    """Auto-assemble a Radia IMA string from sym_<bc>_<axis> labels.

    Scans the list of mesh boundary labels for the canonical
    `sym_bn=0_<axis>` / `sym_ht=0_<axis>` patterns created by
    `add_kelvin.add_kelvin_cubit(reduction=...)` and maps them:

        sym_bn=0_<axis> -> '+<axis>'   (mirror-symmetric, Radia '+')
        sym_ht=0_<axis> -> '-<axis>'   (antisymmetric,    Radia '-')

    Legacy labels "sym_tangential" / "sym_normal" have no axis, so they
    are ignored here (the user should pass --ima explicitly for
    pre-2026-04-25 samples).

    The output axes come in canonical x/y/z order, e.g.
    {'+x', '-z'} -> '+x-z'.  Empty string if no sym_*_* labels found.
    """
    # Accept either iterable of labels or a single string.
    if isinstance(boundaries, str):
        boundaries = [boundaries]
    sign_by_axis = {}
    for label in boundaries:
        if not isinstance(label, str):
            continue
        for bc, sign in (("bn=0", "+"), ("ht=0", "-")):
            prefix = f"sym_{bc}_"
            if label.startswith(prefix):
                axis = label[len(prefix):].lower()
                if axis in ("x", "y", "z"):
                    sign_by_axis[axis] = sign
                break
    return "".join(f"{sign_by_axis[a]}{a}"
                   for a in ("x", "y", "z") if a in sign_by_axis)


def _load_coil_script(script_path):
    """Load coil input -> CoilBuilder.

    Accepts either a CAD STEP file (`.step` / `.stp`) or a Python module
    that exposes ``build_coil() -> CoilBuilder``.  The canonical source
    of truth is the STEP; the Python path is a legacy escape hatch.
    """
    if not script_path or not os.path.exists(script_path):
        raise FileNotFoundError(f"Coil script not found: {script_path}")

    setup_paths()

    # EM panel contract (2026-04-25): .py only.  See
    # calc_accel_magnet._load_coil_script for the rationale.
    ext = os.path.splitext(script_path)[1].lower()
    if ext in ('.step', '.stp'):
        raise ValueError(
            f"EM panel does not accept STEP coil inputs.  Use a Python "
            f"module that defines `build_coil() -> CoilBuilder`.  "
            f"PEEC workflows (IH panel) accept STEP via peec_step.  "
            f"Got: {script_path}")

    script_dir = os.path.dirname(os.path.abspath(script_path))
    if script_dir not in sys.path:
        sys.path.insert(0, script_dir)

    spec = importlib.util.spec_from_file_location("coil_module", script_path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)

    if not hasattr(mod, 'build_coil'):
        raise AttributeError(
            f"Coil script must define build_coil() function: {script_path}")

    return mod.build_coil()


def _extract_elements_from_mesh(mesh, material_name="yoke"):
    """Extract hex/tet/wedge elements for a given material from an NGSolve mesh.

    Replaces the pre-policy `_extract_elements_from_cubit` helper.  The
    `export netgen` C++ command emits coordinates in meters (the
    lab-wide unit) so no coordinate scaling is needed.

    Args:
        mesh: ngsolve.Mesh
        material_name: Material name whose elements should be extracted
                       (matches `el.mat` from the .vol).

    Returns:
        list of dicts: [{'type': 'hex'|'tet'|'wedge', 'vertices': [[x,y,z],...]}]
    """
    from ngsolve import VOL
    from ngsolve.fem import ET

    elements = []
    target = material_name.lower()
    for el in mesh.Elements(VOL):
        if el.mat.lower() != target:
            continue
        # NGSolve element type -> Radia polyhedron type
        if el.type == ET.TET:
            etype = 'tet'
        elif el.type == ET.HEX:
            etype = 'hex'
        elif el.type == ET.PRISM:
            etype = 'wedge'  # triangular prism
        else:
            # PYRAMID and other exotic types are not supported here yet.
            continue
        verts = [list(mesh.vertices[v.nr].point) for v in el.vertices]
        elements.append({'type': etype, 'vertices': verts})

    return elements


def solve_hdiv(coil_script="", vol_file="",
              mat=None, ima="", solver=0,
              max_iter=100, tol=1e-3, relax=0.0,
              msh_output="", demag_backend="hdiv", hdiv_order=1):
    """HDiv-VIM solver: Netgen .vol mesh -> MeshSoftIron -> rad.Solve.

    Args:
        coil_script: Path to Python script with build_coil()
        vol_file: Netgen .vol file with hex/tet/wedge elements for a
            material named 'yoke' (sole interface between Cubit and
            NGSolve; this script does NOT import cubit).  Coordinates
            are in meters (export netgen policy).
        mat: EMMaterial instance (iron yoke properties)
        ima: IMA string e.g. '+x-z' for quarter model ('' = no IMA)
        solver: 0=LU, 1=BiCGSTAB, 2=HACApK
        max_iter: Max nonlinear iterations
        tol: Convergence tolerance
        relax: Under-relaxation (0=full step)
        msh_output: Optional GMSH .msh output path
        demag_backend: only 'hdiv' is exposed by this panel
            (the FEEC HDiv-VIM via radia.vim.MeshSoftIron + rad.Solve).
            The 'hdiv' backend is KELVIN-LESS and iron-only: the .vol must
            contain ONLY the 'yoke' volume material (no air / kelvin elements),
            and IMA symmetry is not supported there yet (both -> fail-loud).
        hdiv_order: HDiv finite-element order, 1 (RT1) or 2 (RT2).

    Returns:
        dict with B_center, n_elem, ndof, iterations, converged, etc.
    """
    if mat is None:
        mat = EMMaterial.from_name("steel")
    material = mat.name
    mu_r = mat.mu_r
    bh_file = ""  # BH data now lives in mat.bh_curve

    setup_paths()
    import radia as rad

    t_total_start = time.perf_counter()

    # Clear ALL prior Radia state ONCE, up front -- before building the coil.
    # (Bug fix 2026-06-17: this UtiDelAll used to live at the start of "Step 3",
    # AFTER the coil was built, so it DESTROYED the coil_container -- rad.UtiDelAll
    # invalidates every existing handle.  model = ObjCnt([yoke, coil_container]) then
    # referenced a stale coil index, silently dropping the coil's field.  This solver had
    # no end-to-end golden (only a static --help smoke), so it went unnoticed.)
    rad.UtiDelAll()

    # ============================================================
    # Step 1: Load coil -> Radia objects (exact analytical)
    # ============================================================
    _log("COIL:loading script")
    coil = _load_coil_script(coil_script)
    current = coil.current
    radia_coil_objs = coil.to_radia()
    coil_container = rad.ObjCnt(radia_coil_objs)
    _log(f"COIL:{len(radia_coil_objs)} Radia objects, I={current}A")

    # ============================================================
    # Step 2: Extract 'yoke' elements from the .vol via NGSolve
    # ============================================================
    _log("MESH:loading .vol")

    if not vol_file or not os.path.exists(vol_file):
        return {"error": f".vol file required: {vol_file}"}

    from ngsolve import Mesh
    mesh = Mesh(vol_file)

    if "yoke" not in set(mesh.GetMaterials()):
        return {"error": "Material 'yoke' not found in .vol "
                         f"(materials: {sorted(set(mesh.GetMaterials()))})"}

    # Auto-assemble IMA string from sym_*_* sidesets when the caller
    # left --ima empty.  Explicit --ima overrides (precedence: CLI >
    # auto-detected labels).
    if not ima:
        mesh_bnds = list(mesh.GetBoundaries())
        ima_auto = ima_from_mesh_labels(mesh_bnds)
        if ima_auto:
            _log(f"IMA:auto-detected {ima_auto!r} from sym_*_* labels "
                 f"(pass --ima to override)")
            ima = ima_auto

    elements = _extract_elements_from_mesh(mesh, material_name="yoke")

    n_hex = sum(1 for e in elements if e['type'] == 'hex')
    n_tet = sum(1 for e in elements if e['type'] == 'tet')
    n_wedge = sum(1 for e in elements if e['type'] == 'wedge')
    n_total = len(elements)
    if n_total == 0:
        return {"error": "No volume elements found with material 'yoke'"}

    _log(f"MESH:{n_hex} hex + {n_tet} tet + {n_wedge} wedge = {n_total} elem")

    # ============================================================
    # Step 3: Create Radia objects with material  (+ Step 4: solve)
    # ============================================================
    _log("RADIA:creating objects")

    is_linear = (material == "linear")
    bh_data = None
    if not is_linear:
        bh_data, _ = mat.get_bh_curve()
        if bh_data is None:
            return {"error": f"No BH curve for material: {material}"}

    solver_names = {0: "LU", 1: "BiCGSTAB", 2: "HACApK"}
    n_dof = n_hex * 6 + n_tet * 3 + n_wedge * 5

    # --- FEEC HDiv-VIM backend (radia.vim) -- the only backend exposed by this panel.  KELVIN-less, IRON-ONLY:
    # the VIM solves the WHOLE registered mesh as iron, so the .vol must contain ONLY 'yoke' volume
    # elements (no air / kelvin).  Production supports one pure TET / HEX / WEDGE topology.
    if demag_backend != "hdiv":
        return {"error": "demag_backend=%r is not available in this panel. Use demag_backend='hdiv' "
                         "(the FEEC HDiv-VIM)."
                         % (demag_backend,)}
    from ngsolve import VOL as _VOL, TaskManager as _TM
    n_nonyoke = sum(1 for el in mesh.Elements(_VOL) if el.mat != "yoke")
    if n_nonyoke > 0:
        return {"error": "the HDiv-VIM needs an IRON-ONLY .vol (it is KELVIN-less; air/kelvin regions "
                         "are not used).  This .vol has %d non-'yoke' volume elements (materials %s).  "
                         "Export the yoke block alone." % (n_nonyoke, sorted(set(mesh.GetMaterials())))}
    if sum(count > 0 for count in (n_tet, n_hex, n_wedge)) != 1:
        return {"error": "the HDiv-VIM panel requires one pure TET, HEX, or WEDGE topology; "
                         "mixed 3D meshes require HDiv pyramid transition elements.  "
                         "This mesh has %d tet, %d hex, and %d wedge elements."
                         % (n_tet, n_hex, n_wedge)}
    hdiv_order = int(hdiv_order)
    if hdiv_order not in (1, 2):
        return {"error": "hdiv_order must be 1 (RT1) or 2 (RT2)"}
    if ima:
        return {"error": "IMA/image symmetry is not exposed by this panel yet. Run the full "
                         "HDiv-VIM model with --ima empty."}
    import radia.vim as _vim
    if is_linear:
        iron = _vim.MeshSoftIron(mesh, mu_r=float(mu_r), order=hdiv_order)
        _log(f"MAT:linear mu_r={mu_r} (HDiv-VIM)")
    else:
        iron = _vim.MeshSoftIron(mesh, bh_table=bh_data, order=hdiv_order)
        _log(f"MAT:nonlinear BH ({len(bh_data)} pts) (HDiv-VIM)")
    model = rad.ObjCnt([iron, coil_container])
    topology = "tet" if n_tet else ("hex" if n_hex else "wedge")
    _log(f"SOLVE:backend=hdiv (FEEC HDiv-VIM RT{hdiv_order}/{topology}), "
         f"tol={tol}, maxiter={max_iter}")
    t_solve_start = time.perf_counter()
    with _TM():
        res = rad.Solve(model, tol, max_iter, solver, demag_backend="hdiv")
    t_solve = time.perf_counter() - t_solve_start
    # the HDiv dispatch returns the radia.vim.Solve dict (raises on non-convergence)
    n_iter = int(res.get("iters", 0)) if isinstance(res, dict) else 0
    n_dof = int(res.get("ndof", n_dof)) if isinstance(res, dict) else n_dof
    M_avg = [float(c) for c in res["M_avg"]] if isinstance(res, dict) else [0.0, 0.0, 0.0]
    residual = 0.0
    converged = True
    solver_label = "HDiv-VIM"

    _log(f"SOLVE:done in {t_solve:.2f}s, {n_iter} iter, "
         f"residual={residual:.2e}")

    # ============================================================
    # Step 5: Evaluate field at origin
    # ============================================================
    B = rad.Fld(model, 'b', [0, 0, 0])
    B_mag = math.sqrt(B[0]**2 + B[1]**2 + B[2]**2)
    _log(f"FIELD:B(0)=({B[0]*1e3:.2f}, {B[1]*1e3:.2f}, {B[2]*1e3:.2f}) mT, "
         f"|B|={B_mag*1e3:.2f} mT")

    # ============================================================
    # Step 6: Optional GMSH output (element-center B field)
    # Write .vol + .sol first, then vol2msh for GMSH visualization
    # ============================================================
    gmsh_file = ""
    if msh_output:
        try:
            gmsh_file = _write_hdiv_vol_sol_msh(
                msh_output, elements, model, rad)
            _log(f"GMSH:{gmsh_file}")
        except Exception as e:
            _log(f"GMSH:export failed: {e}")

    t_total = time.perf_counter() - t_total_start

    result_dict = {
        "B_center": [float(B[0]), float(B[1]), float(B[2])],
        "B_center_mag": float(B_mag),
        "M_avg": [float(M_avg[0]), float(M_avg[1]), float(M_avg[2])],
        "n_elem": n_total,
        "n_hex": n_hex,
        "n_tet": n_tet,
        "n_wedge": n_wedge,
        "ndof": n_dof,
        "solver": solver_label,
        "demag_backend": demag_backend,
        "hdiv_order": hdiv_order,
        "ima": ima,
        "iterations": n_iter,
        "residual": residual,
        "converged": converged,
        "material": material,
        "mu_r": mu_r if is_linear else None,
        "current": current,
        "t_solve": t_solve,
        "t_total": t_total,
        "msh_file": gmsh_file,
    }

    rad.UtiDelAll()
    return result_dict


def _build_ngmesh_from_elements(elements):
    """Build a netgen Mesh from the Radia hex/tet/wedge element list.

    Shares vertices by coordinate rounding so adjacent elements
    reference the same PointIndex.
    """
    from netgen.meshing import (Mesh as NGMesh, MeshPoint, Element3D,
                                 Pnt, FaceDescriptor)
    m = NGMesh(dim=3)
    m.Add(FaceDescriptor(surfnr=1, domin=1, bc=1))
    m.SetMaterial(1, 'yoke')

    vmap = {}  # (rx,ry,rz) -> PointIndex
    for elem in elements:
        pids = []
        for v in elem['vertices']:
            key = (round(v[0], 10), round(v[1], 10), round(v[2], 10))
            pid = vmap.get(key)
            if pid is None:
                pid = m.Add(MeshPoint(Pnt(v[0], v[1], v[2])))
                vmap[key] = pid
            pids.append(pid)
        m.Add(Element3D(1, pids))
    return m


def _write_hdiv_vol_sol_msh(msh_path, elements, model, rad):
    """.vol + .sol first, then vol2msh to GMSH v4.1 with ElementData.

    Outputs:
      <base>.vol      netgen mesh (hex/tet/wedge)
      <base>_B.sol    L2 order=0 dim=3 element-centered B [T]
      <base>_Bmag.sol L2 order=0 scalar element-centered |B| [T]
      <base>.msh      GMSH v4.1 with $ElementData views (via vol2msh)
    """
    from ngsolve import Mesh as NMesh, L2, GridFunction
    from ngsolve import TaskManager

    base, _ = os.path.splitext(msh_path)
    vol_path = base + ".vol"
    b_sol = base + "_B.sol"
    bmag_sol = base + "_Bmag.sol"

    ngm = _build_ngmesh_from_elements(elements)
    ngm.Save(vol_path)

    # Reload so FES/GridFunction build on the saved-and-loaded mesh
    # (keeps element ordering deterministic between save and later vol2msh).
    from netgen.meshing import Mesh as _NGMesh
    ng2 = _NGMesh()
    ng2.Load(vol_path)
    nsm = NMesh(ng2)

    with TaskManager():
        fes_v = L2(nsm, order=0, dim=3)
        gf_b = GridFunction(fes_v)
        fes_s = L2(nsm, order=0)
        gf_bm = GridFunction(fes_s)

        # Compute B at NGSolve element centroids directly — same centroids
        # the source element list used, since we built the mesh from them.
        for ngs_idx, el in enumerate(nsm.Elements()):
            pts = [nsm[v].point for v in el.vertices]
            n = len(pts)
            cx = sum(p[0] for p in pts) / n
            cy = sum(p[1] for p in pts) / n
            cz = sum(p[2] for p in pts) / n
            B = rad.Fld(model, 'b', [cx, cy, cz])
            gf_b.vec.FV()[ngs_idx * 3 + 0] = float(B[0])
            gf_b.vec.FV()[ngs_idx * 3 + 1] = float(B[1])
            gf_b.vec.FV()[ngs_idx * 3 + 2] = float(B[2])
            gf_bm.vec.FV()[ngs_idx] = math.sqrt(
                float(B[0]) ** 2 + float(B[1]) ** 2 + float(B[2]) ** 2)

        gf_b.Save(b_sol)
        gf_bm.Save(bmag_sol)

        from gmsh_post_export import vol2msh
        vol2msh(
            msh_path, vol_path,
            [
                {"name": "B [T]", "sol": b_sol, "fes": "L2",
                 "fes_order": 0, "fes_dim": 3, "ncomp": 3, "cell_data": True},
                {"name": "|B| [T]", "sol": bmag_sol, "fes": "L2",
                 "fes_order": 0, "ncomp": 1, "cell_data": True},
            ],
        )
        return msh_path


def build_argparser():
    """argparse factory shared by main() and notebook DesignSpec callers."""
    parser = argparse.ArgumentParser(
        description="Accelerator magnet HDiv-VIM solver")
    parser.add_argument("--coil-script", required=True,
                        help="Python script with build_coil() -> CoilBuilder")
    parser.add_argument("--vol", default="",
                        help="Netgen .vol file with a 'yoke' material "
                             "(sole interface between Cubit and NGSolve; "
                             "this script no longer accepts .cub5)")
    # include_custom=True so the EM panel's HDiv-VIM build_command
    # (which sends `--material custom --mu-r <user_value>` for the
    # "mu_r (Linear)" radio) is accepted at argparse.  Without this
    # the panel command fails silently with
    # `error: argument --material: invalid choice: 'custom'` --
    # caught 2026-04-26.  See validation_test/panels/test_em_hdiv_smoke.py.
    add_material_args(parser, default_material="steel",
                      include_custom=True, include_hys=True)
    parser.add_argument("--ima", default="",
                        help="IMA string: '+x-z' (quarter), '+x' (half-x), "
                             "'-z' (half-z), '' (full)")
    parser.add_argument("--solver", type=int, default=0,
                        choices=[0, 1, 2],
                        help="0=LU, 1=BiCGSTAB, 2=HACApK")
    parser.add_argument("--demag-backend", default="hdiv",
                        choices=["hdiv"],
                        help="Demag backend: 'hdiv' (the FEEC HDiv-VIM; KELVIN-less, requires an "
                             "iron-only .vol, does not support IMA symmetry yet).")
    parser.add_argument("--hdiv-order", type=int, choices=[1, 2], default=1,
                        help="HDiv finite-element order: 1=RT1, 2=RT2")
    parser.add_argument("--max-iter", type=int, default=100,
                        help="Max nonlinear iterations")
    parser.add_argument("--tol", type=float, default=1e-3,
                        help="Convergence tolerance")
    parser.add_argument("--relax", type=float, default=0.0,
                        help="Under-relaxation (0=full step)")
    parser.add_argument("--msh-output", default="",
                        help="GMSH .msh output path")
    parser.add_argument("--output", default="",
                        help="JSON output file")
    return parser


def main():
    parser = build_argparser()

    def run(args):
        return solve_hdiv(
            coil_script=args.coil_script,
            vol_file=args.vol,
            mat=EMMaterial.from_args(args),
            ima=args.ima,
            solver=args.solver,
            max_iter=args.max_iter,
            tol=args.tol,
            relax=args.relax,
            msh_output=args.msh_output,
            demag_backend=args.demag_backend,
            hdiv_order=args.hdiv_order,
        )

    calc_main(run, parser)


if __name__ == "__main__":
    main()
