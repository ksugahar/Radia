"""
Inductance extractor using ngsolve.bem LaplaceSL BEM.

Called as subprocess from Cubit panel:
    python calc_inductance.py --cub5 model.cub5 --order 1

Workflow:
  1. Import NGSolve FIRST (before cubit)
  2. Open cub5, register blocks
  3. export_NGSolveCurvedMesh(order=N, surface_only=True) -> NGSolve Mesh
  4. HDivSurface + LaplaceSL -> L extraction

IMPORTANT: NGSolve must be imported BEFORE cubit.
Outputs JSON to stdout (suppresses all other print output).
"""

import argparse
import json
import os
import sys


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


def extract_inductance(cub5_file, order, msh_output=""):
    """Extract self-inductance using BEM LaplaceSL.

    Uses export_NGSolveCurvedMesh() for mesh export with automatic curving:
      export_NGSolveCurvedMesh(order=N, surface_only=True) -> HDivSurface + LaplaceSL
    """
    import numpy as np
    import math
    from ngsolve import HDivSurface, TaskManager, ds, Integrate, CF, BND, GridFunction, Norm
    from ngsolve.bem import LaplaceSL

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

    # Register blocks (use existing mesh, no remesh)
    # Support both tri (tet mesh) and quad (hex mesh) surface elements
    cubit.cmd("set duplicate block elements on")
    for bid in list(cubit.get_block_id_list()):
        cubit.cmd(f"delete block {bid}")
    for i, vid in enumerate(vol_ids):
        cubit.cmd(f"block {i + 1} add volume {vid}")
    bnd_id = len(vol_ids) + 1
    # Always try both tri and face(quad) - Cubit silently skips if none exist
    cubit.cmd(f"block {bnd_id} add tri all")
    cubit.cmd(f"block {bnd_id} add face all")
    cubit.cmd(f'block {bnd_id} name "conductor"')

    # Export to NGSolve mesh with curving
    radia_src = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
    if os.path.abspath(radia_src) not in sys.path:
        sys.path.insert(0, os.path.abspath(radia_src))
    import cubit_mesh_export

    mesh = cubit_mesh_export.export_NGSolveCurvedMesh(
        cubit, order=order, surface_only=True,
        split_quads=True)

    # HDivSurface + LaplaceSL
    fes = HDivSurface(mesh, order=0)
    u, v = fes.TnT()
    ndof = fes.ndof
    fd = fes.FreeDofs()
    free_idx = np.array([i for i in range(ndof) if fd[i]])
    n_free = len(free_idx)

    bnd_label = list(set(mesh.GetBoundaries()))
    label = bnd_label[0] if bnd_label else None

    # --- Phase 1: Export mesh to GMSH immediately (before solve) ---
    gmsh_file = msh_output if msh_output else os.path.join(
        os.path.dirname(os.path.abspath(cub5_file)), "inductance_J.msh"
    ).replace("\\", "/")
    try:
        from gmsh_post_export import GmshPostExport
        post = GmshPostExport(mesh, boundary=True)
        post.write_mesh(gmsh_file)
        sys.stderr.write(f"MESH_READY:{gmsh_file}\n")
        sys.stderr.flush()
    except Exception:
        pass

    # --- Phase 2: BEM solve ---
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

    # --- Phase 3: Rewrite .msh with field data ---
    # Per-node |J| for visualization
    gf = GridFunction(fes)
    for il, ig in enumerate(free_idx):
        gf.vec[ig] = x[il]

    elem_J = Integrate(Norm(gf), mesh, VOL_or_BND=BND, element_wise=True)
    elem_A = Integrate(CF(1), mesh, VOL_or_BND=BND, element_wise=True)

    num_nodes = mesh.nv
    node_sum = np.zeros(num_nodes)
    node_cnt = np.zeros(num_nodes)
    for el in mesh.Elements(BND):
        avg_J = abs(elem_J[el.nr]) / max(abs(elem_A[el.nr]), 1e-30)
        for v in el.vertices:
            node_sum[v.nr] += avg_J
            node_cnt[v.nr] += 1
    node_J = np.where(node_cnt > 0, node_sum / node_cnt, 0.0)

    # Overwrite mesh-only .msh with mesh + field data
    try:
        post = GmshPostExport(mesh, boundary=True)
        post.add_field("|J|", node_J, ncomp=1)
        post.write(gmsh_file)
        sys.stderr.write(f"FIELD_READY:{gmsh_file}\n")
        sys.stderr.flush()
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
    parser.add_argument("--msh-output", default="", help="GMSH .msh output path")
    parser.add_argument("--output", default="", help="Output JSON file (optional)")
    args = parser.parse_args()

    # Redirect stdout during computation (suppress Cubit/BEM print output)
    import io as _io
    real_stdout = sys.stdout
    sys.stdout = _io.StringIO()
    try:
        result = extract_inductance(args.cub5, args.order, args.msh_output)
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
