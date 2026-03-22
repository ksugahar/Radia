"""
Inductance extractor using source/sink saddle point EFIE + ngsolve.bem.

Called as subprocess from Cubit panel:
    python calc_inductance.py --cub5 model.cub5 --order 2 --source source --sink sink

Method (saddle point EFIE):
  1. export_NGSolveCurvedMesh(surface_only=True) -> surface mesh
     Block names become NGSolve boundary labels ("source", "sink", etc.)
  2. Build SL (LaplaceSL) and D (divergence) matrices
  3. Solve saddle point: [SL D^T; D 0] [J; p] = [0; g]
     where g encodes unit current injection at source/sink
  4. L = mu_0 * J^T @ SL @ J
  5. GMSH v2.2 export with |J| field

Prerequisites (Cubit journal):
    block N name "source"    # tri/quad on source face
    block M name "sink"      # tri/quad on sink face
    block K name "boundary"  # all surface tri/quad (optional)

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


def _to_dense(mat):
    """Extract dense NumPy array from NGSolve BaseMatrix via ToDense()."""
    return mat.ToDense().NumPy()


def extract_inductance(cub5_file, order, source_label="source",
                       sink_label="sink", fes_order=0, msh_output=""):
    """Extract self-inductance via source/sink saddle point EFIE.

    Args:
        cub5_file: Path to Cubit .cub5 model
        order: Curve order (1-2)
        source_label: Block name for source face
        sink_label: Block name for sink face
        fes_order: HDivSurface basis order (0=RWG, 1+=higher-order)
        msh_output: Optional path for GMSH .msh output

    Returns:
        dict with inductance_H, diagnostics
    """
    import numpy as np
    import math
    from scipy.linalg import solve as scipy_solve
    from ngsolve import (HDivSurface, SurfaceL2, TaskManager, ds,
                         Integrate, CF, BND, GridFunction, Norm,
                         BilinearForm, LinearForm, div)
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

    # Check that source/sink blocks exist
    block_names = []
    for bid in cubit.get_block_id_list():
        block_names.append(cubit.get_exodus_entity_name("block", bid))

    has_source = source_label in block_names
    has_sink = sink_label in block_names

    if not has_source or not has_sink:
        return {"error": f"Missing blocks: "
                f"{'source' if not has_source else ''} "
                f"{'sink' if not has_sink else ''}. "
                f"Available blocks: {block_names}. "
                f"Define source/sink blocks in Cubit journal."}

    # Export surface mesh (block names -> boundary labels)
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

    # Verify boundary labels
    bnd_labels = list(set(mesh.GetBoundaries()))
    if source_label not in bnd_labels or sink_label not in bnd_labels:
        return {"error": f"Boundary labels {bnd_labels} do not contain "
                f"'{source_label}' and/or '{sink_label}'."}

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

    # --- Phase 2: Saddle point EFIE ---
    # HDivSurface order can be 0 (RWG) or higher; constraint always order=0
    fes_J = HDivSurface(mesh, order=fes_order)
    fes_L2 = SurfaceL2(mesh, order=0)  # element-wise current conservation
    n_J, n_f = fes_J.ndof, fes_L2.ndof

    # Divergence matrix D: n_f x n_J
    u_J = fes_J.TrialFunction()
    q = fes_L2.TestFunction()
    bf_D = BilinearForm(trialspace=fes_J, testspace=fes_L2)
    bf_D += div(u_J.Trace()) * q * ds
    bf_D.Assemble()
    D = _to_dense(bf_D.mat)

    # LaplaceSL: n_J x n_J
    jt, jv = fes_J.TnT()
    with TaskManager():
        V_op = LaplaceSL(jt.Trace() * ds, use_fmm=False) * jv.Trace() * ds
        SL = V_op.mat.ToDense().NumPy()

    # Source/sink RHS
    f_src = LinearForm(fes_L2)
    f_src += q * ds(source_label)
    f_src.Assemble()
    g_src = f_src.vec.FV().NumPy().copy()
    A_src = np.sum(g_src)

    f_snk = LinearForm(fes_L2)
    f_snk += q * ds(sink_label)
    f_snk.Assemble()
    g_snk = f_snk.vec.FV().NumPy().copy()
    A_snk = np.sum(g_snk)

    if A_src < 1e-30 or A_snk < 1e-30:
        return {"error": f"Source/sink faces empty (A_src={A_src}, A_snk={A_snk})."}

    g = g_src / A_src - g_snk / A_snk

    # Saddle point solve (remove last constraint for regularity)
    D_red = D[:-1, :]
    g_red = g[:-1]
    n_c = n_f - 1

    K = np.block([
        [SL,              D_red.T],
        [D_red, np.zeros((n_c, n_c))]
    ])
    rhs = np.zeros(n_J + n_c)
    rhs[n_J:] = g_red

    x = scipy_solve(K, rhs)
    J = x[:n_J]

    # Inductance: L = mu_0 * J^T @ SL @ J
    L_total = MU_0 * J @ SL @ J
    residual = float(np.max(np.abs(D @ J - g)))

    # --- Phase 3: GMSH with |J| field ---
    gf_J = GridFunction(fes_J)
    gf_J.vec.FV().NumPy()[:] = J

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
        post.add_field("|J|", node_J, ncomp=1)
        post.write(gmsh_file)
        sys.stderr.write(f"FIELD_READY:{gmsh_file}\n")
        sys.stderr.flush()
    except Exception:
        gmsh_file = ""

    return {
        "inductance_H": float(L_total),
        "n_dofs": n_J,
        "n_faces": n_f,
        "euler": euler,
        "surface_area": area,
        "source_area": float(A_src),
        "sink_area": float(A_snk),
        "constraint_residual": residual,
        "curve_order": order,
        "fes_order": fes_order,
        "gmsh_file": gmsh_file,
    }


def main():
    parser = argparse.ArgumentParser(
        description="BEM inductance (source/sink saddle point EFIE)")
    parser.add_argument("--cub5", required=True, help="Cubit .cub5 model file")
    parser.add_argument("--order", type=int, default=2, help="Curve order (1-2)")
    parser.add_argument("--fes-order", type=int, default=0, help="HDivSurface order (0=RWG)")
    parser.add_argument("--source", default="source", help="Source block name")
    parser.add_argument("--sink", default="sink", help="Sink block name")
    parser.add_argument("--msh-output", default="", help="GMSH .msh output path")
    parser.add_argument("--output", default="", help="Output JSON file (optional)")
    args = parser.parse_args()

    import io as _io
    real_stdout = sys.stdout
    sys.stdout = _io.StringIO()
    try:
        result = extract_inductance(args.cub5, args.order,
                                    args.source, args.sink,
                                    args.fes_order, args.msh_output)
    except Exception as e:
        result = {"error": str(e)}
    sys.stdout = real_stdout

    if args.output:
        with open(args.output, "w", encoding="utf-8") as f:
            json.dump(result, f)
    print(json.dumps(result))


if __name__ == "__main__":
    main()
