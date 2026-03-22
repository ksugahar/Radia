"""
Inductance extractor using Hodge decomposition + ngsolve.bem LaplaceSL.

Called as subprocess from Cubit panel:
    python calc_inductance.py --cub5 model.cub5 --order 2

Method (EFIE-based, NOT energy method):
  1. export_NGSolveCurvedMesh(surface_only=True) -> surface mesh
  2. Hodge decomposition: harmonic = ker([D; C^T M_J])
  3. LaplaceSL projected onto harmonic subspace
  4. Generalized eigenvalue -> toroidal mode -> L
  5. GMSH v2.2 export with |J| field

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
    """Extract self-inductance via Hodge decomposition + LaplaceSL."""
    import numpy as np
    import math
    from scipy.linalg import null_space, eigh
    from ngsolve import (HDivSurface, SurfaceL2, TaskManager, ds,
                         Integrate, CF, BND, GridFunction, Norm,
                         BilinearForm, InnerProduct, div)
    from ngsolve.bem import LaplaceSL

    MU_0 = 4.0 * math.pi * 1e-7

    cubit = _setup_cubit()
    cubit.cmd(f'open "{cub5_file}"')

    vol_ids = list(cubit.get_entities("volume"))
    if not vol_ids:
        return {"error": "No volumes found in model"}

    total_elems = sum(
        len(cubit.get_volume_tets(vid)) + len(cubit.get_volume_hexes(vid))
        for vid in vol_ids
    )
    if total_elems == 0:
        return {"error": "Volumes are not meshed."}

    # Register blocks
    cubit.cmd("set duplicate block elements on")
    for bid in list(cubit.get_block_id_list()):
        cubit.cmd(f"delete block {bid}")
    for i, vid in enumerate(vol_ids):
        cubit.cmd(f"block {i + 1} add volume {vid}")
    bnd_id = len(vol_ids) + 1
    cubit.cmd(f"block {bnd_id} add tri all")
    cubit.cmd(f"block {bnd_id} add face all")
    cubit.cmd(f'block {bnd_id} name "conductor"')

    # Estimate R and a from bounding box
    bb = cubit.get_bounding_box("volume", vol_ids[0])
    x_span = abs(bb[1] - bb[0])
    z_span = abs(bb[7] - bb[6])
    R_est = x_span / 2.0
    a_est = z_span / 2.0

    # Export surface mesh
    radia_src = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
    if os.path.abspath(radia_src) not in sys.path:
        sys.path.insert(0, os.path.abspath(radia_src))
    import cubit_mesh_export

    mesh = cubit_mesh_export.export_NGSolveCurvedMesh(
        cubit, order=order, surface_only=True, split_quads=True)

    nse = mesh.GetNE(BND)
    nv = mesh.nv
    ne = sum(1 for e in mesh.edges)
    euler = nv - ne + nse
    area = float(Integrate(CF(1), mesh, VOL_or_BND=BND))

    # --- Phase 1: Export mesh to GMSH (before solve) ---
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

    # --- Phase 2: Hodge decomposition ---
    def _extract_rect(mat, nrow, ncol):
        M = np.zeros((nrow, ncol))
        ei = mat.CreateRowVector(); col = mat.CreateColVector()
        for j in range(ncol):
            ei[:] = 0; ei[j] = 1.0; mat.Mult(ei, col)
            for i in range(nrow): M[i, j] = col[i]
        return M

    def _extract_dense(mat, n):
        M = np.zeros((n, n))
        ei = mat.CreateColVector(); col = mat.CreateColVector()
        for j in range(n):
            ei[:] = 0; ei[j] = 1.0; mat.Mult(ei, col)
            for i in range(n): M[i, j] = col[i]
        return M

    fes_J = HDivSurface(mesh, order=0)
    fes_L2 = SurfaceL2(mesh, order=0)
    n_J, n_f = fes_J.ndof, fes_L2.ndof

    u_J, q = fes_J.TrialFunction(), fes_L2.TestFunction()
    bf_D = BilinearForm(trialspace=fes_J, testspace=fes_L2)
    bf_D += div(u_J.Trace()) * q * ds
    bf_D.Assemble()
    D = _extract_rect(bf_D.mat, n_f, n_J)

    C = np.zeros((n_J, nv))
    for i, e in enumerate(mesh.edges):
        vv = list(e.vertices)
        if i < n_J:
            C[i, vv[0].nr] = -1
            C[i, vv[1].nr] = +1

    u, v = fes_J.TnT()
    bf_M = BilinearForm(fes_J)
    bf_M += InnerProduct(u.Trace(), v.Trace()) * ds
    bf_M.Assemble()
    M_J = _extract_dense(bf_M.mat, n_J)

    c_h_all = null_space(np.vstack([D, C.T @ M_J]), rcond=1e-10)
    n_harm = c_h_all.shape[1]

    if n_harm == 0:
        return {"error": "No harmonic modes found (check Euler characteristic)"}

    # --- Phase 3: LaplaceSL + eigenvalue ---
    jt, jv = fes_J.TnT()
    with TaskManager():
        V_op = LaplaceSL(jt.Trace() * ds, use_fmm=False) * jv.Trace() * ds
        SL = V_op.mat.ToDense().NumPy()

    eigvals, eigvecs = eigh(c_h_all.T @ SL @ c_h_all,
                            c_h_all.T @ M_J @ c_h_all)

    idx = 1 if n_harm >= 2 else 0
    L_total = MU_0 * eigvals[idx] * R_est / a_est

    # --- Phase 4: GMSH with |J| field ---
    c_tor = c_h_all @ eigvecs[:, idx]
    gf_J = GridFunction(fes_J)
    gf_J.vec.FV().NumPy()[:] = c_tor

    elem_J = Integrate(Norm(gf_J), mesh, VOL_or_BND=BND, element_wise=True)
    elem_A = Integrate(CF(1), mesh, VOL_or_BND=BND, element_wise=True)
    ns, nc = np.zeros(nv), np.zeros(nv)
    for el in mesh.Elements(BND):
        aj = abs(elem_J[el.nr]) / max(abs(elem_A[el.nr]), 1e-30)
        for vtx in el.vertices:
            ns[vtx.nr] += aj
            nc[vtx.nr] += 1
    node_J = np.where(nc > 0, ns / nc, 0.0)

    try:
        post = GmshPostExport(mesh, boundary=True)
        post.add_field("|J_toroidal|", node_J, ncomp=1)
        post.write(gmsh_file)
        sys.stderr.write(f"FIELD_READY:{gmsh_file}\n")
        sys.stderr.flush()
    except Exception:
        gmsh_file = ""

    return {
        "inductance_H": float(L_total),
        "n_dofs": n_J,
        "n_harmonic": n_harm,
        "euler": euler,
        "surface_area": area,
        "R_est": R_est,
        "a_est": a_est,
        "order": order,
        "gmsh_file": gmsh_file,
    }


def main():
    parser = argparse.ArgumentParser(description="BEM inductance (Hodge decomposition)")
    parser.add_argument("--cub5", required=True, help="Cubit .cub5 model file")
    parser.add_argument("--order", type=int, default=2, help="Curve order (1-2)")
    parser.add_argument("--msh-output", default="", help="GMSH .msh output path")
    parser.add_argument("--output", default="", help="Output JSON file (optional)")
    args = parser.parse_args()

    import io as _io
    real_stdout = sys.stdout
    sys.stdout = _io.StringIO()
    try:
        result = extract_inductance(args.cub5, args.order, args.msh_output)
    except Exception as e:
        result = {"error": str(e)}
    sys.stdout = real_stdout

    if args.output:
        with open(args.output, "w", encoding="utf-8") as f:
            json.dump(result, f)
    print(json.dumps(result))


if __name__ == "__main__":
    main()
