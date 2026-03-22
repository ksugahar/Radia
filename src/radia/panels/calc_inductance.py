"""
Inductance extractor using ngsolve.bem LaplaceSL BEM (energy method).

Called as subprocess from Cubit panel:
    python calc_inductance.py --cub5 model.cub5 --order 1

Energy method: L = mu_0 * J^T * SL * J  (for I=1 toroidal current)
where SL is the LaplaceSL single-layer operator on the conductor surface.

Workflow:
  1. Import NGSolve FIRST (before cubit)
  2. Open cub5, register blocks
  3. export_NGSolveCurvedMesh(order=N, surface_only=True) -> NGSolve Mesh
  4. HDivSurface + LaplaceSL (use_fmm=False, ToDense()) -> SL matrix
  5. Project toroidal current onto HDivSurface
  6. L = mu_0 * J^T * SL * J

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


def extract_inductance(cub5_file, order, msh_output=""):
    """Extract self-inductance using BEM LaplaceSL energy method.

    Energy method: L = mu_0 * J^T * SL * J
    where J is the toroidal surface current for I=1, projected onto HDivSurface.

    Uses export_NGSolveCurvedMesh() for mesh export with automatic curving:
      export_NGSolveCurvedMesh(order=N, surface_only=True) -> HDivSurface + LaplaceSL
    """
    import numpy as np
    import math
    from ngsolve import (HDivSurface, TaskManager, ds, Integrate, CF, BND,
                         GridFunction, Norm, sqrt, x, y, z)
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
    cubit.cmd("set duplicate block elements on")
    for bid in list(cubit.get_block_id_list()):
        cubit.cmd(f"delete block {bid}")
    for i, vid in enumerate(vol_ids):
        cubit.cmd(f"block {i + 1} add volume {vid}")
    bnd_id = len(vol_ids) + 1
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

    # HDivSurface
    fes = HDivSurface(mesh, order=0)
    u, v = fes.TnT()
    ndof = fes.ndof

    # Toroidal current for I=1: J = e_phi / (2*pi*a)
    # Estimate minor radius from bounding box
    bb = cubit.get_bounding_box("volume", vol_ids[0])
    z_span = abs(bb[7] - bb[6])  # Z extent = 2*a for a torus in XY plane
    a_est = z_span / 2.0

    r_cf = sqrt(x*x + y*y)
    J_toroidal = CF((-y/r_cf, x/r_cf, 0)) / (2 * math.pi * a_est)

    gf_J = GridFunction(fes)
    gf_J.Set(J_toroidal, definedon=mesh.Boundaries(".*"), dual=True)
    J_vec = gf_J.vec.FV().NumPy().copy()

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

    # --- Phase 2: BEM (energy method) ---
    # use_fmm=False for reproducibility and faster dense extraction.
    # TaskManager must be used for BOTH setup and extraction, or neither.
    # See: Joachim Schoeberl feedback (2026-03-22)
    with TaskManager():
        L_op = LaplaceSL(u.Trace() * ds, use_fmm=False) * v.Trace() * ds
        SL = L_op.mat.ToDense().NumPy()

    # Energy method: L = mu_0 * J^T * SL * J  (for I=1)
    L_total = MU_0 * float(J_vec @ SL @ J_vec)

    # Surface area for verification
    area = float(Integrate(CF(1), mesh, VOL_or_BND=BND))

    # --- Phase 3: Rewrite .msh with field data ---
    elem_J = Integrate(Norm(gf_J), mesh, VOL_or_BND=BND, element_wise=True)
    elem_A = Integrate(CF(1), mesh, VOL_or_BND=BND, element_wise=True)

    num_nodes = mesh.nv
    node_sum = np.zeros(num_nodes)
    node_cnt = np.zeros(num_nodes)
    for el in mesh.Elements(BND):
        avg_J = abs(elem_J[el.nr]) / max(abs(elem_A[el.nr]), 1e-30)
        for vtx in el.vertices:
            node_sum[vtx.nr] += avg_J
            node_cnt[vtx.nr] += 1
    node_J = np.where(node_cnt > 0, node_sum / node_cnt, 0.0)

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
        "n_dofs": ndof,
        "surface_area": area,
        "order": order,
        "gmsh_file": gmsh_file,
    }


def main():
    parser = argparse.ArgumentParser(description="ngsolve.bem inductance extractor (energy method)")
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
