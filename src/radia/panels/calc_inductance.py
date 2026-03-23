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

    # --- Phase 3: Per-node J vector via vertex averaging ---
    gf_J = GridFunction(fes_J)
    gf_J.vec.FV().NumPy()[:] = J

    elem_A = Integrate(CF(1), mesh, VOL_or_BND=BND, element_wise=True)
    elem_Jx = Integrate(gf_J[0], mesh, VOL_or_BND=BND, element_wise=True)
    elem_Jy = Integrate(gf_J[1], mesh, VOL_or_BND=BND, element_wise=True)
    elem_Jz = Integrate(gf_J[2], mesh, VOL_or_BND=BND, element_wise=True)

    J_nodes = np.zeros((nv, 3))
    nc_count = np.zeros(nv)
    for el in mesh.Elements(BND):
        a = max(abs(elem_A[el.nr]), 1e-30)
        jvec = [elem_Jx[el.nr] / a, elem_Jy[el.nr] / a, elem_Jz[el.nr] / a]
        for vtx in el.vertices:
            J_nodes[vtx.nr] += jvec
            nc_count[vtx.nr] += 1
    for k in range(3):
        J_nodes[:, k] = np.where(nc_count > 0, J_nodes[:, k] / nc_count, 0.0)

    sys.stderr.write(f"FIELD_READY:J computed\n")
    sys.stderr.flush()

    # --- Phase 4: B-distribution (volume B via Biot-Savart) ---
    B_nodes, vol_nodes, vol_elems = _compute_B_field(
        mesh, gf_J, elem_A, MU_0)

    sys.stderr.write(f"B_FIELD_READY:B computed\n")
    sys.stderr.flush()

    # --- Phase 5: Write single combined .msh (no node ID collision) ---
    gmsh_file = os.path.join(base_dir, "inductance.msh").replace("\\", "/")
    _write_combined_msh(gmsh_file, mesh, J_nodes, vol_nodes, vol_elems, B_nodes)

    # Companion .geo
    geo_file = os.path.join(base_dir, "inductance.geo").replace("\\", "/")
    with open(geo_file, 'w', encoding='utf-8') as f:
        f.write(f'Merge "{os.path.basename(gmsh_file)}";\n')
        f.write('Mesh.NumSubEdges = 4;\n')
        f.write('Mesh.VolumeEdges = 0;\n')
        f.write('View[0].ArrowSizeMin = 20;\n')
        f.write('View[0].ArrowSizeMax = 20;\n')
        f.write('View[1].ArrowSizeMin = 20;\n')
        f.write('View[1].ArrowSizeMax = 20;\n')

    # COMSOL-compatible text interpolation files
    _write_comsol_txt(os.path.join(base_dir, "B_field.txt"), vol_nodes, B_nodes,
                      "B", ["Bx", "By", "Bz"])
    _write_comsol_txt(os.path.join(base_dir, "J_surface.txt"),
                      [(mesh.vertices[i].point[0], mesh.vertices[i].point[1],
                        mesh.vertices[i].point[2]) for i in range(nv)],
                      J_nodes, "J", ["Jx", "Jy", "Jz"])

    # Nastran BDF with CTRIA6 (2nd order surface mesh for COMSOL geometry import)
    _write_nastran_ctria6(os.path.join(base_dir, "coil_surface.bdf"), mesh)

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
        "gmsh_file": geo_file,
    }


def _write_nastran_ctria6(filename, mesh_surf):
    """Write coil surface as Nastran BDF with CTRIA6 (2nd order triangles).

    Uses NGSolve GetTrafo to compute mid-edge node positions on the
    curved surface (from mesh.Curve(2)). Edge nodes are cached to avoid
    duplicates at shared edges.

    Nastran CTRIA6 node ordering: v0 v1 v2 m01 m12 m20
    (corners then mid-edge, counterclockwise)
    """
    import numpy as np
    from ngsolve import BND, IntegrationRule

    # Vertex nodes
    nodes = [(v.point[0], v.point[1], v.point[2]) for v in mesh_surf.vertices]

    # Compute mid-edge nodes via GetTrafo
    edge_cache = {}  # (min_v, max_v) -> (start_vertex, mid_node_index)

    elements = []
    for el in mesh_surf.Elements(BND):
        verts = [v.nr for v in el.vertices]
        if len(verts) != 3:
            continue

        trafo = mesh_surf.GetTrafo(el)

        # Reference triangle mid-edge points
        # Edge 0->1: (0.5, 0.0), Edge 1->2: (0.5, 0.5), Edge 2->0: (0.0, 0.5)
        ref_mids = [(0.5, 0.0), (0.5, 0.5), (0.0, 0.5)]
        edge_pairs = [(0, 1), (1, 2), (2, 0)]

        # Match ref corners to physical vertices
        ref_corners = [(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)]
        corner_mapped = []
        for u, v in ref_corners:
            ir = IntegrationRule([(u, v)], [1.0])
            for ip in ir:
                mip = trafo(ip)
                corner_mapped.append(
                    np.array([mip.point[0], mip.point[1], mip.point[2]]))

        ref_to_vert = []
        for ci in range(3):
            dists = [np.linalg.norm(corner_mapped[ci] -
                     np.array(mesh_surf.vertices[verts[vi]].point))
                     for vi in range(3)]
            ref_to_vert.append(verts[np.argmin(dists)])

        mid_indices = []
        for ei, (rc0, rc1) in enumerate(edge_pairs):
            va = ref_to_vert[rc0]
            vb = ref_to_vert[rc1]
            edge_key = (min(va, vb), max(va, vb))

            if edge_key in edge_cache:
                mid_indices.append(edge_cache[edge_key])
            else:
                u, v = ref_mids[ei]
                ir = IntegrationRule([(u, v)], [1.0])
                for ip in ir:
                    mip = trafo(ip)
                    pt = (mip.point[0], mip.point[1], mip.point[2])
                    nidx = len(nodes)
                    nodes.append(pt)
                    edge_cache[edge_key] = nidx
                    mid_indices.append(nidx)

        # CTRIA6: v0 v1 v2 m01 m12 m20 (Nastran convention)
        elements.append([ref_to_vert[0], ref_to_vert[1], ref_to_vert[2],
                         mid_indices[0], mid_indices[1], mid_indices[2]])

    # Write Nastran BDF
    with open(filename, 'w', encoding='utf-8') as f:
        f.write('$ Nastran BDF - Coil surface mesh (CTRIA6)\n')
        f.write('$ Generated by Radia BEM inductance panel\n')
        f.write('BEGIN BULK\n')

        # GRID cards
        for i, (x, y, z) in enumerate(nodes):
            nid = i + 1
            f.write(f'GRID    {nid:8d}        {x:8.5e}{y:8.5e}{z:8.5e}\n')

        # CTRIA6 cards
        for i, conn in enumerate(elements):
            eid = i + 1
            n = [c + 1 for c in conn]  # 1-indexed
            f.write(f'CTRIA6  {eid:8d}{1:8d}'
                    f'{n[0]:8d}{n[1]:8d}{n[2]:8d}{n[3]:8d}{n[4]:8d}{n[5]:8d}\n')

        f.write('ENDDATA\n')


def _write_comsol_txt(filename, nodes, field_data, field_name, comp_names):
    """Write COMSOL-compatible text interpolation file.

    COMSOL import: Global Definitions > Interpolation > File
    Format: % x y z Fx Fy Fz (space-separated, one header line)

    Args:
        filename: Output file path
        nodes: list of (x, y, z) or array (n, 3)
        field_data: array (n, 3) field values
        field_name: e.g. "B"
        comp_names: e.g. ["Bx", "By", "Bz"]
    """
    import numpy as np
    nodes = np.asarray(nodes)
    field_data = np.asarray(field_data)
    with open(filename, 'w', encoding='utf-8') as f:
        f.write(f'% x y z {" ".join(comp_names)}\n')
        for i in range(len(nodes)):
            x, y, z = nodes[i]
            vals = ' '.join(f'{field_data[i, c]:.15e}' for c in range(field_data.shape[1]))
            f.write(f'{x:.15e} {y:.15e} {z:.15e} {vals}\n')


def _compute_B_field(mesh_surf, gf_J, elem_A, MU_0):
    """Compute B field in air volume via direct Biot-Savart.

    Returns:
        B_nodes: (nv_vol, 3) B field at volume vertices
        vol_nodes: list of (x, y, z) volume vertex coordinates
        vol_elems: list of (v0, v1, v2, v3) tet connectivity (0-indexed)
    """
    import numpy as np
    from ngsolve import Mesh, BND, CF, Integrate, VOL
    from netgen.occ import Box, Pnt, OCCGeometry
    from netgen.meshing import MeshingParameters

    INV_4PI = 1.0 / (4.0 * np.pi)

    # Extract per-element J and geometry
    elem_Jx = Integrate(gf_J[0], mesh_surf, VOL_or_BND=BND, element_wise=True)
    elem_Jy = Integrate(gf_J[1], mesh_surf, VOL_or_BND=BND, element_wise=True)
    elem_Jz = Integrate(gf_J[2], mesh_surf, VOL_or_BND=BND, element_wise=True)

    centroids, areas, J_vecs = [], [], []
    for el in mesh_surf.Elements(BND):
        area = abs(elem_A[el.nr])
        if area < 1e-30:
            continue
        jvec = np.array([elem_Jx[el.nr], elem_Jy[el.nr], elem_Jz[el.nr]]) / area
        verts = [mesh_surf.vertices[v.nr].point for v in el.vertices]
        c = np.mean([(v[0], v[1], v[2]) for v in verts], axis=0)
        centroids.append(c); areas.append(area); J_vecs.append(jvec)

    centroids = np.array(centroids)
    areas = np.array(areas)
    J_vecs = np.array(J_vecs)

    # Create air volume mesh
    coords = np.array([(v.point[0], v.point[1], v.point[2])
                       for v in mesh_surf.vertices])
    bbox_min, bbox_max = coords.min(axis=0), coords.max(axis=0)
    extent = bbox_max - bbox_min
    margin = max(extent) * 0.3

    box = Box(Pnt(*(bbox_min - margin)), Pnt(*(bbox_max + margin)))
    maxh_vol = max(extent) * 0.1
    mesh_vol = Mesh(OCCGeometry(box).GenerateMesh(
        mp=MeshingParameters(maxh=maxh_vol)))

    nv_vol = mesh_vol.nv
    obs_pts = np.array([(v.point[0], v.point[1], v.point[2])
                        for v in mesh_vol.vertices])

    sys.stderr.write(f"B_PROGRESS:{len(centroids)} elems -> {nv_vol} pts\n")
    sys.stderr.flush()

    # Biot-Savart B (vectorized)
    dx = obs_pts[:, None, :] - centroids[None, :, :]
    r = np.sqrt(np.maximum(np.sum(dx**2, axis=2), 1e-30))
    r3_inv = areas[None, :] / (r ** 3)
    cross = np.cross(J_vecs[None, :, :], dx)
    B_nodes = MU_0 * INV_4PI * np.sum(cross * r3_inv[:, :, None], axis=1)

    # Extract volume mesh topology
    vol_nodes = [(v.point[0], v.point[1], v.point[2]) for v in mesh_vol.vertices]
    vol_elems = []
    for el in mesh_vol.Elements(VOL):
        vol_elems.append([v.nr for v in el.vertices])

    return B_nodes, vol_nodes, vol_elems


def _write_combined_msh(filename, mesh_surf, J_nodes, vol_nodes, vol_elems, B_nodes):
    """Write single .msh v2.2 with volume B + surface J + coil wireframe.

    Node numbering:
      1..nv_vol                  : volume mesh vertices (B field)
      nv_vol+1..nv_vol+nv_surf   : surface mesh vertices (J field)
    Elements:
      1..ne_vol                  : volume tets (physical group "air")
      ne_vol+1..ne_vol+ne_surf   : surface tris (physical group "coil_surface")
      after that                 : coil wireframe lines (physical group "coil_wire")
    """
    import numpy as np
    from ngsolve import BND

    nv_vol = len(vol_nodes)
    nv_surf = mesh_surf.nv
    surf_offset = nv_vol  # surface node IDs start after volume nodes

    # Surface nodes
    surf_nodes = [(v.point[0], v.point[1], v.point[2]) for v in mesh_surf.vertices]

    # Surface triangles
    surf_elems = []
    for el in mesh_surf.Elements(BND):
        surf_elems.append([v.nr for v in el.vertices])

    # Coil wireframe edges
    edges = set()
    for el in mesh_surf.Elements(BND):
        verts = [v.nr for v in el.vertices]
        nv = len(verts)
        for i in range(nv):
            a, b = verts[i], verts[(i + 1) % nv]
            edges.add((min(a, b), max(a, b)))
    edges = sorted(edges)

    n_nodes = nv_vol + nv_surf
    ne_vol = len(vol_elems)
    ne_surf = len(surf_elems)
    ne_wire = len(edges)
    n_elems = ne_vol + ne_surf + ne_wire

    with open(filename, 'w', encoding='utf-8') as f:
        f.write('$MeshFormat\n2.2 0 8\n$EndMeshFormat\n')

        # Physical groups
        f.write('$PhysicalNames\n3\n')
        f.write('3 1 "air"\n')
        f.write('2 2 "coil_surface"\n')
        f.write('1 3 "coil_wire"\n')
        f.write('$EndPhysicalNames\n')

        # Nodes: volume then surface
        f.write(f'$Nodes\n{n_nodes}\n')
        for i, (x, y, z) in enumerate(vol_nodes):
            f.write(f'{i + 1} {x:.15e} {y:.15e} {z:.15e}\n')
        for i, (x, y, z) in enumerate(surf_nodes):
            f.write(f'{surf_offset + i + 1} {x:.15e} {y:.15e} {z:.15e}\n')
        f.write('$EndNodes\n')

        # Elements: volume tets + surface tris + wireframe lines
        f.write(f'$Elements\n{n_elems}\n')
        eid = 1
        for verts in vol_elems:
            ns = ' '.join(str(v + 1) for v in verts)
            f.write(f'{eid} 4 2 1 1 {ns}\n')  # type 4 = tet
            eid += 1
        for verts in surf_elems:
            ns = ' '.join(str(v + surf_offset + 1) for v in verts)
            f.write(f'{eid} 2 2 2 2 {ns}\n')  # type 2 = tri
            eid += 1
        for a, b in edges:
            f.write(f'{eid} 1 2 3 3 {a + surf_offset + 1} {b + surf_offset + 1}\n')
            eid += 1
        f.write('$EndElements\n')

        # NodeData: B on volume nodes
        f.write('$NodeData\n1\n"B"\n1\n0.0\n3\n0\n3\n')
        f.write(f'{nv_vol}\n')
        for i in range(nv_vol):
            bx, by, bz = B_nodes[i]
            f.write(f'{i + 1} {bx:.15e} {by:.15e} {bz:.15e}\n')
        f.write('$EndNodeData\n')

        # NodeData: J on surface nodes (offset IDs)
        f.write('$NodeData\n1\n"J"\n1\n0.0\n3\n0\n3\n')
        f.write(f'{nv_surf}\n')
        for i in range(nv_surf):
            jx, jy, jz = J_nodes[i]
            f.write(f'{surf_offset + i + 1} {jx:.15e} {jy:.15e} {jz:.15e}\n')
        f.write('$EndNodeData\n')


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
