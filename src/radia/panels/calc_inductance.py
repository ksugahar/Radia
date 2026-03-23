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
  5. GMSH v2.2 export: J-distribution (surface) and B-distribution (volume)

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
        msh_output: Optional path for GMSH .msh output (J-distribution)

    Returns:
        dict with inductance_H, gmsh_file_J, gmsh_file_B, diagnostics
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
    # Output dir: same as msh_output (all files must be co-located for .geo Merge)
    if msh_output:
        base_dir = os.path.dirname(os.path.abspath(msh_output))
        gmsh_file_J = msh_output
    else:
        base_dir = os.path.dirname(os.path.abspath(cub5_file))
        gmsh_file_J = os.path.join(base_dir, "inductance_J.msh").replace("\\", "/")
    try:
        from gmsh_post_export import GmshPostExport
        post = GmshPostExport(mesh, boundary=True)
        post.write_mesh(gmsh_file_J)
        sys.stderr.write(f"MESH_READY:{gmsh_file_J}\n")
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

    # --- Phase 3: J-distribution GMSH (surface |J| + J vector) ---
    gf_J = GridFunction(fes_J)
    gf_J.vec.FV().NumPy()[:] = J

    # Per-node |J| (scalar) and J (vector) via vertex averaging
    elem_J_mag = Integrate(Norm(gf_J), mesh, VOL_or_BND=BND, element_wise=True)
    elem_A = Integrate(CF(1), mesh, VOL_or_BND=BND, element_wise=True)

    # J vector components at element centroids
    elem_Jx = Integrate(gf_J[0], mesh, VOL_or_BND=BND, element_wise=True)
    elem_Jy = Integrate(gf_J[1], mesh, VOL_or_BND=BND, element_wise=True)
    elem_Jz = Integrate(gf_J[2], mesh, VOL_or_BND=BND, element_wise=True)

    ns_mag = np.zeros(nv)
    ns_vec = np.zeros((nv, 3))
    nc_count = np.zeros(nv)
    for el in mesh.Elements(BND):
        a = max(abs(elem_A[el.nr]), 1e-30)
        mag = abs(elem_J_mag[el.nr]) / a
        jvec = [elem_Jx[el.nr] / a, elem_Jy[el.nr] / a, elem_Jz[el.nr] / a]
        for vtx in el.vertices:
            ns_mag[vtx.nr] += mag
            ns_vec[vtx.nr] += jvec
            nc_count[vtx.nr] += 1

    node_J_mag = np.where(nc_count > 0, ns_mag / nc_count, 0.0)
    for k in range(3):
        ns_vec[:, k] = np.where(nc_count > 0, ns_vec[:, k] / nc_count, 0.0)

    try:
        post = GmshPostExport(mesh, boundary=True)
        post.add_field("|J|", node_J_mag, ncomp=1)
        post.add_field("J", ns_vec, ncomp=3)
        post.write(gmsh_file_J)
        sys.stderr.write(f"FIELD_READY:{gmsh_file_J}\n")
        sys.stderr.flush()
    except Exception:
        gmsh_file_J = ""

    # --- Phase 4: B-distribution (volume B via Biot-Savart) ---
    gmsh_file_B = ""
    try:
        gmsh_file_B = _compute_B_distribution(
            mesh, fes_J, gf_J, jt, base_dir, MU_0)
    except Exception as e:
        sys.stderr.write(f"B_FIELD_ERROR:{e}\n")
        sys.stderr.flush()

    # --- Phase 5: Combined .geo (all views in one GMSH window) ---
    gmsh_geo = os.path.join(base_dir, "inductance.geo").replace("\\", "/")
    _write_combined_geo(gmsh_geo, gmsh_file_J, gmsh_file_B, base_dir)

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
        "gmsh_file": gmsh_geo,
    }


def _compute_B_distribution(mesh_surf, fes_J, gf_J, jt, base_dir, MU_0):
    """Compute B field in air volume via Biot-Savart and export to GMSH.

    Creates an air volume mesh (box around conductor), evaluates the
    vector potential A via Biot-Savart quadrature over surface elements,
    then projects into HCurl and computes B = curl(A).

    Args:
        mesh_surf: Surface-only BEM mesh
        fes_J: HDivSurface space on mesh_surf
        gf_J: Solved surface current GridFunction
        jt: Trial function of fes_J
        base_dir: Output directory for .msh file
        MU_0: Vacuum permeability

    Returns:
        Path to B-distribution .msh file
    """
    import numpy as np
    from ngsolve import Mesh, BND, CF, Integrate
    from netgen.occ import Box, Pnt, OCCGeometry
    from netgen.meshing import MeshingParameters

    INV_4PI = 1.0 / (4.0 * np.pi)

    # --- Step 1: Extract per-element J and geometry from surface mesh ---
    # For each BND element: centroid, area, average J vector
    elem_data = []
    elem_A = Integrate(CF(1), mesh_surf, VOL_or_BND=BND, element_wise=True)
    elem_Jx = Integrate(gf_J[0], mesh_surf, VOL_or_BND=BND, element_wise=True)
    elem_Jy = Integrate(gf_J[1], mesh_surf, VOL_or_BND=BND, element_wise=True)
    elem_Jz = Integrate(gf_J[2], mesh_surf, VOL_or_BND=BND, element_wise=True)

    for el in mesh_surf.Elements(BND):
        area = abs(elem_A[el.nr])
        if area < 1e-30:
            continue
        jvec = np.array([elem_Jx[el.nr], elem_Jy[el.nr], elem_Jz[el.nr]]) / area
        verts = [mesh_surf.vertices[v.nr].point for v in el.vertices]
        centroid = np.mean([(v[0], v[1], v[2]) for v in verts], axis=0)
        elem_data.append((centroid, area, jvec))

    n_elem = len(elem_data)
    sys.stderr.write(f"B_PROGRESS:Biot-Savart over {n_elem} elements\n")
    sys.stderr.flush()

    # Pack into arrays for vectorized computation
    centroids = np.array([d[0] for d in elem_data])  # (n_elem, 3)
    areas = np.array([d[1] for d in elem_data])       # (n_elem,)
    J_vecs = np.array([d[2] for d in elem_data])       # (n_elem, 3)

    # --- Step 2: Create air volume mesh ---
    coords = np.array([(v.point[0], v.point[1], v.point[2])
                       for v in mesh_surf.vertices])
    bbox_min = coords.min(axis=0)
    bbox_max = coords.max(axis=0)
    extent = bbox_max - bbox_min
    margin = max(extent) * 0.3

    box = Box(Pnt(*(bbox_min - margin)), Pnt(*(bbox_max + margin)))
    for f in box.faces:
        f.name = "outer"

    maxh_vol = max(extent) * 0.2  # coarse volume mesh (visualization only)
    geo = OCCGeometry(box)
    mesh_vol = Mesh(geo.GenerateMesh(
        mp=MeshingParameters(maxh=maxh_vol)))

    nv_vol = mesh_vol.nv
    obs_pts = np.array([(v.point[0], v.point[1], v.point[2])
                        for v in mesh_vol.vertices])  # (nv_vol, 3)

    sys.stderr.write(f"B_PROGRESS:Biot-Savart {n_elem} elems -> {nv_vol} pts\n")
    sys.stderr.flush()

    # --- Step 3: Biot-Savart B directly (fully vectorized) ---
    # B(x) = mu_0/(4*pi) * sum_e (J_e x (x - c_e)) * A_e / |x - c_e|^3
    # Vectorized over all observation points at once.
    # obs_pts: (nv, 3), centroids: (ne, 3) -> dx: (nv, ne, 3)
    dx = obs_pts[:, None, :] - centroids[None, :, :]  # (nv, ne, 3)
    r2 = np.sum(dx**2, axis=2)  # (nv, ne)
    r = np.sqrt(np.maximum(r2, 1e-30))  # (nv, ne)
    r3_inv = areas[None, :] / (r * r * r)  # (nv, ne)

    # Cross product J x dx for all (obs, elem) pairs
    # J_vecs: (ne, 3) broadcast to (nv, ne, 3)
    Jx = J_vecs[None, :, :]  # (1, ne, 3)
    cross = np.cross(Jx, dx)  # (nv, ne, 3)

    B_nodes = MU_0 * INV_4PI * np.sum(cross * r3_inv[:, :, None], axis=1)  # (nv, 3)

    # --- Step 4: Compute |B| per node ---
    node_B_mag = np.sqrt(np.sum(B_nodes**2, axis=1))

    from gmsh_post_export import GmshPostExport
    gmsh_file_B = os.path.join(base_dir, "inductance_B.msh").replace("\\", "/")
    post = GmshPostExport(mesh_vol, boundary=False)
    post.add_field("|B|", node_B_mag, ncomp=1)
    post.add_field("B", B_nodes, ncomp=3)
    post.write(gmsh_file_B)

    # Write coil wireframe as 1D line elements (always visible on top of volume)
    coil_file = os.path.join(base_dir, "inductance_coil.msh").replace("\\", "/")
    _write_coil_wireframe(mesh_surf, coil_file)

    # Write .geo that merges B field + coil wireframe
    geo_file = os.path.splitext(gmsh_file_B)[0] + '.geo'
    with open(geo_file, 'w', encoding='utf-8') as f:
        f.write('// B-distribution with coil wireframe overlay\n')
        f.write(f'Merge "{os.path.basename(gmsh_file_B)}";\n')
        if os.path.exists(coil_file):
            f.write(f'Merge "{os.path.basename(coil_file)}";\n')
        f.write('Mesh.NumSubEdges = 4;\n')

    sys.stderr.write(f"B_FIELD_READY:{gmsh_file_B}\n")
    sys.stderr.flush()

    return gmsh_file_B


def _write_coil_wireframe(mesh_surf, filename):
    """Write coil surface edges as 1D line elements in GMSH v2.2 format.

    1D elements are always visible on top of 3D volume elements in GMSH,
    providing a clear coil outline overlay on the B-distribution.
    """
    import numpy as np
    from ngsolve import BND

    # Collect unique edges from boundary elements
    edges = set()
    for el in mesh_surf.Elements(BND):
        verts = [v.nr for v in el.vertices]
        nv = len(verts)
        for i in range(nv):
            a, b = verts[i], verts[(i + 1) % nv]
            edges.add((min(a, b), max(a, b)))

    # Get vertex coordinates
    nodes = [(v.point[0], v.point[1], v.point[2]) for v in mesh_surf.vertices]
    n_nodes = len(nodes)
    n_edges = len(edges)

    with open(filename, 'w', encoding='utf-8') as f:
        f.write('$MeshFormat\n2.2 0 8\n$EndMeshFormat\n')
        f.write('$PhysicalNames\n1\n1 1 "coil"\n$EndPhysicalNames\n')
        f.write(f'$Nodes\n{n_nodes}\n')
        for i, (x, y, z) in enumerate(nodes):
            f.write(f'{i + 1} {x:.15e} {y:.15e} {z:.15e}\n')
        f.write('$EndNodes\n')
        f.write(f'$Elements\n{n_edges}\n')
        for idx, (a, b) in enumerate(sorted(edges)):
            f.write(f'{idx + 1} 1 2 1 1 {a + 1} {b + 1}\n')  # type 1 = 2-node line
        f.write('$EndElements\n')


def _write_combined_geo(geo_file, gmsh_file_J, gmsh_file_B, base_dir):
    """Write a .geo that merges all result files into one GMSH window.

    GMSH tree will show:
      Post-processing
        [0] |B|   (volume)
        [1] B     (volume vectors)
        [2] |J|   (surface, visible on coil)
        [3] J     (surface vectors)
    Plus coil wireframe as 1D mesh overlay.
    """
    coil_file = os.path.join(base_dir, "inductance_coil.msh").replace("\\", "/")
    with open(geo_file, 'w', encoding='utf-8') as f:
        f.write('// Combined inductance visualization\n')
        # B-distribution first (volume mesh establishes 3D context)
        if gmsh_file_B and os.path.exists(gmsh_file_B):
            f.write(f'Merge "{os.path.basename(gmsh_file_B)}";\n')
        # J-distribution (surface mesh + field views overlay)
        if gmsh_file_J and os.path.exists(gmsh_file_J):
            f.write(f'Merge "{os.path.basename(gmsh_file_J)}";\n')
        # Coil wireframe (1D lines, always visible on top)
        if os.path.exists(coil_file):
            f.write(f'Merge "{os.path.basename(coil_file)}";\n')
        f.write('Mesh.NumSubEdges = 4;\n')
        # Hide volume mesh edges for cleaner view
        f.write('Mesh.VolumeEdges = 0;\n')


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
