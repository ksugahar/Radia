"""
Inductance extractor using ngsolve.bem LaplaceSL BEM.

Called as subprocess from Cubit panel:
    python calc_inductance.py --cub5 model.cub5 --order 1

Workflow (STEP reimport + CalcSurfacesOfNode):
  1. Import NGSolve FIRST (before cubit)
  2. Open cub5, get CAD info, export STEP
  3. Reset + reimport STEP heal (ACIS -> OCC seam fix)
  4. Remesh + export_NetgenMesh(geometry=OCC)
  5. CalcSurfacesOfNode() (CRITICAL for BEM edge orientation)
  6. SetGeomInfo + Curve(order)
  7. HDivSurface + LaplaceSL -> L extraction

IMPORTANT: NGSolve must be imported BEFORE cubit.
Outputs JSON to stdout (suppresses all other print output).
"""

import argparse
import json
import os
import sys
import tempfile


def _setup_cubit():
    """Import cubit with path cleanup (NGSolve already imported)."""
    radia_src = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
    if os.path.abspath(radia_src) not in sys.path:
        sys.path.insert(0, os.path.abspath(radia_src))
    from install_panels import find_cubit_bin

    cubit_path = find_cubit_bin()
    if cubit_path and cubit_path not in sys.path:
        sys.path.append(cubit_path)

    for p in list(sys.path):
        if "site-packages" in p and ("cubit" in p.lower() or "Cubit" in p):
            sys.path.remove(p)

    import cubit
    import io, contextlib
    with contextlib.redirect_stderr(io.StringIO()):
        cubit.init(["cubit", "-nojournal", "-batch", "-noinit"])

    for p in list(sys.path):
        if "site-packages" in p and ("cubit" in p.lower() or "Cubit" in p):
            sys.path.remove(p)

    return cubit


def _get_mesh_size(cubit, vol_ids):
    """Estimate mesh size from existing mesh."""
    for vid in vol_ids:
        ne = len(cubit.get_volume_tets(vid)) + len(cubit.get_volume_hexes(vid))
        if ne > 0:
            bb = cubit.get_bounding_box("volume", vid)
            diag = bb[9] if len(bb) > 9 else 1.0
            return diag / max(ne ** (1.0 / 3.0), 1.0)
    return 0.1


def extract_inductance(cub5_file, order):
    """Extract self-inductance using BEM LaplaceSL.

    Uses the full STEP reimport workflow:
      STEP export -> reimport heal -> remesh -> export_NetgenMesh(geometry=OCC)
      -> CalcSurfacesOfNode -> SetGeomInfo -> Curve -> HDivSurface + LaplaceSL
    """
    import numpy as np
    import math
    from ngsolve import Mesh, HDivSurface, TaskManager, ds, Integrate, CF, BND, GridFunction, Norm
    from ngsolve.bem import LaplaceSL
    from netgen.occ import OCCGeometry

    MU_0 = 4.0 * math.pi * 1e-7

    cubit = _setup_cubit()
    cubit.cmd(f'open "{cub5_file}"')

    vol_ids = list(cubit.get_entities("volume"))
    if not vol_ids:
        return {"error": "No volumes found in model"}

    # Check if meshed
    total_elems = sum(
        len(cubit.get_volume_tets(vid)) + len(cubit.get_volume_hexes(vid))
        for vid in vol_ids
    )
    if total_elems == 0:
        return {"error": "Volumes are not meshed."}

    # Export STEP (for OCC geometry reference only, no reimport)
    tmpdir = tempfile.mkdtemp(prefix="radia_ind_")
    step_file = os.path.join(tmpdir, "geometry.step").replace("\\", "/")
    vol_list = " ".join(str(v) for v in vol_ids)
    cubit.cmd(f'export step "{step_file}" volume {vol_list} overwrite')

    # Register blocks (use existing mesh, no remesh)
    cubit.cmd("delete block all")
    for i, vid in enumerate(vol_ids):
        cubit.cmd(f"block {i + 1} add volume {vid}")
    cubit.cmd(f"block {len(vol_ids) + 1} add tri all")
    cubit.cmd(f'block {len(vol_ids) + 1} name "conductor"')

    # Export to Netgen with OCC geometry
    radia_src = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
    if os.path.abspath(radia_src) not in sys.path:
        sys.path.insert(0, os.path.abspath(radia_src))
    import cubit_mesh_export

    geo = OCCGeometry(step_file)
    ngmesh = cubit_mesh_export.export_NetgenMesh(cubit, geometry=geo)
    # CalcSurfacesOfNode is now called inside export_NetgenMesh

    mesh = Mesh(ngmesh)

    # Curve (with SetGeomInfo already applied if available)
    if order > 1:
        try:
            mesh.Curve(order)
        except Exception:
            pass  # Curve may fail without SetGeomInfo for specific geometry

    # HDivSurface + LaplaceSL
    fes = HDivSurface(mesh, order=0)
    u, v = fes.TnT()
    ndof = fes.ndof
    fd = fes.FreeDofs()
    free_idx = np.array([i for i in range(ndof) if fd[i]])
    n_free = len(free_idx)

    bnd_label = list(set(mesh.GetBoundaries()))
    label = bnd_label[0] if bnd_label else None

    # NOTE: Do NOT use TaskManager() with LaplaceSL - it is not thread-safe
    # and produces non-deterministic results (sign flips in L matrix).
    if label:
        L_op = LaplaceSL(u.Trace() * ds(label)) * v.Trace() * ds(label)
    else:
        L_op = LaplaceSL(u.Trace() * ds) * v.Trace() * ds

    # Extract free-DOF L matrix (O(N^2), may take minutes for large meshes)
    L_free = np.zeros((n_free, n_free))
    ei = L_op.mat.CreateColVector()
    col = L_op.mat.CreateColVector()
    for jl, jg in enumerate(free_idx):
        ei[:] = 0; ei[jg] = 1.0; L_op.mat.Mult(ei, col)
        for il, ig in enumerate(free_idx):
            L_free[il, jl] = col[ig]
    L_matrix = MU_0 * L_free

    # Check for negative diagonals
    neg_count = int(np.sum(np.diag(L_matrix) < 0))

    # Uniform current excitation -> self-inductance
    e = np.ones(n_free) / n_free
    try:
        x = np.linalg.solve(L_matrix, e)
        L_total = 1.0 / (e @ x)
    except np.linalg.LinAlgError:
        x = np.zeros(n_free)
        L_total = 0.0

    # Surface area for verification
    area = float(Integrate(CF(1), mesh, VOL_or_BND=BND))

    # Per-node |J| for visualization
    # Solve L*x = e -> x is the current coefficient vector
    gf = GridFunction(fes)
    for il, ig in enumerate(free_idx):
        gf.vec[ig] = x[il]

    # Element-wise |J| integral / area -> average |J| per element
    elem_J = Integrate(Norm(gf), mesh, VOL_or_BND=BND, element_wise=True)
    elem_A = Integrate(CF(1), mesh, VOL_or_BND=BND, element_wise=True)

    # Distribute to nodes (average over adjacent elements)
    num_nodes = mesh.nv
    node_sum = np.zeros(num_nodes)
    node_cnt = np.zeros(num_nodes)
    for el in mesh.Elements(BND):
        avg_J = abs(elem_J[el.nr]) / max(abs(elem_A[el.nr]), 1e-30)
        for v in el.vertices:
            node_sum[v.nr] += avg_J
            node_cnt[v.nr] += 1
    node_J = np.where(node_cnt > 0, node_sum / node_cnt, 0.0)

    # Export GMSH .msh with per-node |J| as NodeData
    # Place .msh next to input cub5 file
    cub5_dir = os.path.dirname(os.path.abspath(cub5_file))
    gmsh_file = os.path.join(cub5_dir, "inductance_J.msh").replace("\\", "/")
    try:
        from gmsh_post_export import GmshPostExport
        post = GmshPostExport(mesh)
        # Use pre-computed node_J array (per-node scalar from BND integration)
        post.add_field("|J|", node_J, ncomp=1)
        post.write(gmsh_file)
    except Exception:
        gmsh_file = ""

    return {
        "inductance_H": float(L_total),
        "n_free_dofs": n_free,
        "neg_diag": neg_count,
        "surface_area": area,
        "order": order,
        "node_J": node_J.tolist(),
        "gmsh_file": gmsh_file,
    }


def main():
    parser = argparse.ArgumentParser(description="ngsolve.bem inductance extractor")
    parser.add_argument("--cub5", required=True, help="Cubit .cub5 model file")
    parser.add_argument("--order", type=int, default=1, help="Curve order (1-5)")
    parser.add_argument("--source", default="", help="Source block name (future)")
    parser.add_argument("--sink", default="", help="Sink block name (future)")
    parser.add_argument("--sigma", type=float, default=5.8e7, help="Conductivity (future)")
    parser.add_argument("--freq", type=float, default=0.0, help="Frequency (future)")
    parser.add_argument("--output", default="", help="Output JSON file (optional)")
    args = parser.parse_args()

    # Redirect stdout during computation (suppress Cubit/BEM print output)
    import io as _io
    real_stdout = sys.stdout
    sys.stdout = _io.StringIO()
    try:
        result = extract_inductance(args.cub5, args.order)
    except Exception as e:
        result = {"error": str(e)}
    sys.stdout = real_stdout

    # Output: file if specified, otherwise stdout
    if args.output:
        with open(args.output, "w", encoding="utf-8") as f:
            json.dump(result, f)
    print(json.dumps(result))


if __name__ == "__main__":
    main()
