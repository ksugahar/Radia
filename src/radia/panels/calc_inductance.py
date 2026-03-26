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


def extract_inductance(cub5_file, order, source_label="source",
                       sink_label="sink", fes_order=0, msh_output="",
                       workpiece="", impedance_model="esim",
                       frequency=50000, sigma=2e6, half_thickness=0.005,
                       material="steel"):
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
    import time as _time
    from ngsolve import Integrate, CF, BND, GridFunction

    radia_src = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
    if os.path.abspath(radia_src) not in sys.path:
        sys.path.insert(0, os.path.abspath(radia_src))
    from bem_inductance import compute_inductance_source_sink, MU_0

    t_total_start = _time.perf_counter()

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

    t0_export = _time.perf_counter()
    mesh = cubit_mesh_export.export_NGSolveCurvedMesh(
        cubit, order=order, surface_only=True, split_quads=True)
    t_export = _time.perf_counter() - t0_export

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

    # --- Phase 2: Saddle point EFIE (shared solver) ---
    sol = compute_inductance_source_sink(mesh, source_label, sink_label, fes_order)
    if 'error' in sol:
        return sol

    L_total = sol['L']
    n_J = sol['n_J']
    n_f = sol['n_f']
    A_src = sol['A_source']
    A_snk = sol['A_sink']
    residual = sol['residual']
    t_solve = sol['t_total']

    sys.stderr.write(f"SOLVE_DONE:{t_solve:.1f}s\n")
    sys.stderr.flush()

    # Save solve results for post-processing (separate step)
    j_npy = os.path.join(base_dir, "J_coeffs.npy").replace("\\", "/")
    np.save(j_npy, sol['J'])
    mesh_vol = os.path.join(base_dir, "surface_mesh.vol").replace("\\", "/")
    mesh.ngmesh.Save(mesh_vol)

    # Conductor bounding box for default post volume
    coords = np.array([(v.point[0], v.point[1], v.point[2])
                       for v in mesh.vertices])
    bbox_min, bbox_max = coords.min(axis=0), coords.max(axis=0)
    extent = bbox_max - bbox_min
    default_half = float(max(extent) * 0.65)
    default_maxh = float(max(extent) * 0.1)

    result = {
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
        "t_solve": round(t_solve, 2),
        "t_export": round(t_export, 2),
        "t_assembly": sol['t_assembly'],
        "t_lu": sol['t_solve'],
        "j_npy": j_npy,
        "mesh_vol": mesh_vol,
        "default_lxyz": round(default_half, 4),
        "default_maxh": round(default_maxh, 4),
    }

    # --- Workpiece surface impedance (ESIM or Dowell) ---
    if workpiece:
        sys.stderr.write("ESIM_START:computing workpiece impedance\n")
        sys.stderr.flush()
        wp_result = _compute_workpiece_impedance(
            mesh, sol['gf_J'], workpiece, cubit,
            impedance_model=impedance_model,
            frequency=frequency, sigma=sigma,
            half_thickness=half_thickness, material=material)
        result.update(wp_result)
        sys.stderr.write(f"ESIM_DONE:R={wp_result.get('wp_R_effective', 0):.4e}\n")
        sys.stderr.flush()

    return result


def _compute_workpiece_impedance(mesh_coil, gf_J, workpiece_label, cubit_mod,
                                  impedance_model="esim", frequency=50000,
                                  sigma=2e6, half_thickness=0.005,
                                  material="steel"):
    """Compute workpiece surface impedance via ESIM or Dowell.

    Pipeline:
      1. Extract workpiece surface panels from Cubit block geometry
      2. Biot-Savart from coil J -> H at workpiece panels
      3. ESIM cell problem or Dowell formula -> Z_s, P, Q per panel
      4. Integrate over surface

    Returns dict with wp_* keys to merge into main result.
    """
    import math
    from ngsolve import Integrate, CF, BND

    rho = 1.0 / sigma
    omega = 2 * math.pi * frequency

    # --- Step 1: Workpiece surface panels from Cubit block ---
    # Get workpiece volume IDs from block
    wp_vol_ids = []
    for bid in cubit_mod.get_block_id_list():
        try:
            name = cubit_mod.get_exodus_entity_name("block", bid)
            if name == workpiece_label:
                # Block may contain volumes
                for vid in cubit_mod.get_block_volumes(bid):
                    wp_vol_ids.append(vid)
                break
        except Exception:
            pass

    if not wp_vol_ids:
        return {"wp_error": f"Workpiece block '{workpiece_label}' has no volumes"}

    # Get workpiece surface panels: centroid + normal + area from Cubit surfaces
    panels = []
    for vid in wp_vol_ids:
        surfaces = cubit_mod.get_relatives("volume", vid, "surface")
        for sid in surfaces:
            s = cubit_mod.surface(sid)
            cx, cy, cz = s.center_point()
            # Normal from Cubit (outward)
            nx, ny, nz = s.normal_at(s.center_point())
            area = s.area()
            panels.append({
                'center': np.array([cx, cy, cz]),
                'normal': np.array([nx, ny, nz]),
                'area': area,
            })

    if not panels:
        return {"wp_error": "No workpiece surfaces found"}

    n_panels = len(panels)

    # --- Step 2: Biot-Savart H at workpiece panels ---
    INV_4PI = 1.0 / (4.0 * np.pi)

    elem_A = Integrate(CF(1), mesh_coil, VOL_or_BND=BND, element_wise=True)
    elem_Jx = Integrate(gf_J[0], mesh_coil, VOL_or_BND=BND, element_wise=True)
    elem_Jy = Integrate(gf_J[1], mesh_coil, VOL_or_BND=BND, element_wise=True)
    elem_Jz = Integrate(gf_J[2], mesh_coil, VOL_or_BND=BND, element_wise=True)

    centroids, areas_c, J_vecs = [], [], []
    for el in mesh_coil.Elements(BND):
        area = abs(elem_A[el.nr])
        if area < 1e-30:
            continue
        jvec = np.array([elem_Jx[el.nr], elem_Jy[el.nr], elem_Jz[el.nr]]) / area
        verts = [mesh_coil.vertices[v.nr].point for v in el.vertices]
        c = np.mean([(v[0], v[1], v[2]) for v in verts], axis=0)
        centroids.append(c)
        areas_c.append(area)
        J_vecs.append(jvec)

    centroids = np.array(centroids)
    areas_c = np.array(areas_c)
    J_vecs = np.array(J_vecs)

    obs_pts = np.array([p['center'] for p in panels])
    dx = obs_pts[:, None, :] - centroids[None, :, :]
    r = np.sqrt(np.maximum(np.sum(dx**2, axis=2), 1e-30))
    r3_inv = areas_c[None, :] / (r ** 3)
    cross = np.cross(J_vecs[None, :, :], dx)
    H_panels = INV_4PI * np.sum(cross * r3_inv[:, :, None], axis=1)

    # Tangential H
    normals = np.array([p['normal'] for p in panels])
    H_n = np.sum(H_panels * normals, axis=1, keepdims=True)
    H_t_vec = H_panels - H_n * normals
    H_t_mag = np.linalg.norm(H_t_vec, axis=1)

    # --- Step 3: Surface impedance per panel ---
    # BH curve for steel
    STEEL_BH = [
        [0.0, 0.0], [50.0, 0.10], [100.0, 0.25], [200.0, 0.55],
        [500.0, 0.95], [1000.0, 1.20], [2000.0, 1.40], [5000.0, 1.55],
        [10000.0, 1.65], [20000.0, 1.75], [50000.0, 1.90], [100000.0, 2.00],
    ]

    bh_curve = STEEL_BH if material == "steel" else None
    mu_r = 1.0 if material in ("copper", "aluminum") else None

    from esim_cell_problem import ESIMFiniteSlabSolver

    P_total = 0.0
    Q_total = 0.0
    Z_sum = 0.0 + 0.0j
    delta_min = float('inf')
    delta_max = 0.0

    if impedance_model == "esim":
        solver = ESIMFiniteSlabSolver(
            half_thickness=half_thickness, bh_curve=bh_curve,
            sigma=sigma, frequency=frequency,
            mu_r=mu_r if mu_r else 1.0, n_nodes=200)

        for i in range(n_panels):
            H0 = max(float(H_t_mag[i]), 1e-3)
            sol = solver.solve(H0)
            P_total += sol['P_prime'] * panels[i]['area']
            Q_total += sol['Q_prime'] * panels[i]['area']
            Z_sum += sol['Z'] * H_t_mag[i]**2 * panels[i]['area']
            delta_min = min(delta_min, sol['delta'])
            delta_max = max(delta_max, sol['delta'])

    elif impedance_model == "dowell":
        # Dowell: Z_s = (rho/a) * gamma*a * tanh(gamma*a)
        mu_eff = MU_0 * (mu_r if mu_r else 1.0)
        delta = math.sqrt(2 * rho / (omega * mu_eff)) if omega > 0 else 1e10
        xi = half_thickness / delta
        gamma_a = complex(1, 1) * xi
        try:
            Z_s = (rho / half_thickness) * gamma_a * np.tanh(gamma_a)
        except OverflowError:
            Z_s = (rho / half_thickness) * gamma_a

        delta_min = delta_max = delta
        for i in range(n_panels):
            H0 = float(H_t_mag[i])
            P_panel = Z_s.real * H0**2 / 2 * panels[i]['area']
            Q_panel = Z_s.imag * H0**2 / 2 * panels[i]['area']
            P_total += P_panel
            Q_total += Q_panel
            Z_sum += Z_s * H0**2 * panels[i]['area']

    R_eff = float(Z_sum.real)
    X_eff = float(Z_sum.imag)

    return {
        "wp_model": impedance_model,
        "wp_material": material,
        "wp_frequency": frequency,
        "wp_sigma": sigma,
        "wp_half_thickness": half_thickness,
        "wp_n_panels": n_panels,
        "wp_H_t_max": float(H_t_mag.max()),
        "wp_H_t_min": float(H_t_mag.min()),
        "wp_P_total": float(P_total),
        "wp_Q_total": float(Q_total),
        "wp_R_effective": R_eff,
        "wp_X_effective": X_eff,
        "wp_delta_min": float(delta_min),
        "wp_delta_max": float(delta_max),
    }


def post_process(mesh_vol_path, fes_order=0, msh_output="", j_npy="",
                 lx=0, ly=0, lz=0, maxh_vol=0):
    """Post-process: B-field Biot-Savart, GMSH/Nastran/COMSOL export.

    Loads Solve results (mesh .vol + J_coeffs.npy). No Cubit needed.

    Args:
        mesh_vol_path: Path to surface_mesh.vol from solve step
        fes_order: HDivSurface basis order
        msh_output: Optional .msh output path
        j_npy: Path to J_coeffs.npy from solve step
        lx, ly, lz: Box half-sizes [m].
        maxh_vol: Volume mesh element size [m].

    Returns:
        dict with gmsh_file and diagnostics
    """
    import numpy as np
    import time as _time
    from ngsolve import Mesh, Integrate, CF, BND, HDivSurface, GridFunction
    from netgen.meshing import Mesh as NetgenMesh

    radia_src = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
    if os.path.abspath(radia_src) not in sys.path:
        sys.path.insert(0, os.path.abspath(radia_src))
    from bem_inductance import MU_0

    t0 = _time.perf_counter()

    # Load mesh from Solve step (no Cubit needed)
    ngmesh = NetgenMesh()
    ngmesh.Load(mesh_vol_path)
    mesh = Mesh(ngmesh)
    nv = mesh.nv

    # Output directory
    base_dir = os.path.dirname(os.path.abspath(j_npy))
    if msh_output:
        base_dir = os.path.dirname(os.path.abspath(msh_output))

    # Load J coefficients
    J = np.load(j_npy)
    fes_J = HDivSurface(mesh, order=fes_order)
    gf_J = GridFunction(fes_J)
    gf_J.vec.FV().NumPy()[:] = J

    # Per-node J vector via vertex averaging
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

    # B-distribution (volume B via Biot-Savart)
    B_nodes, vol_nodes, vol_elems = _compute_B_field(
        mesh, gf_J, elem_A, MU_0,
        lx=lx, ly=ly, lz=lz, maxh_vol=maxh_vol)

    sys.stderr.write(f"B_FIELD_READY:B computed\n")
    sys.stderr.flush()

    # Write single combined .msh
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

    # Nastran BDF with CTRIA6
    _write_nastran_ctria6(os.path.join(base_dir, "coil_surface.bdf"), mesh)

    t_post = _time.perf_counter() - t0

    return {
        "gmsh_file": geo_file,
        "t_post": round(t_post, 2),
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

    # Write Nastran BDF with strict fixed-format fields
    #
    # GRID* (long format): 2 lines
    #   Line 1: columns 1-8="GRID*   ", 9-24=ID, 25-40=CP, 41-56=X1, 57-72=X2
    #   Line 2: columns 1-8="*       ", 9-24=X3, 25-40=CD
    #   Each field is exactly 16 characters (right-justified)
    #
    # CTRIA6 (short format): 1 line, 10 fields of 8 characters each
    #   columns 1-8="CTRIA6  ", 9-16=EID, 17-24=PID, 25-32=G1, ...
    #   Each field is exactly 8 characters (right-justified)

    def _grid_star(nid, x, y, z):
        """Format GRID* card (long format, 16-char fields)."""
        # Coordinate formatting: fit within 16 chars
        def _fmt16(val):
            s = f'{val:.8e}'
            if len(s) > 16:
                s = f'{val:.6e}'
            if len(s) > 16:
                s = f'{val:.5e}'
            return f'{s:>16}'
        line1 = f'{"GRID*":8s}{nid:>16d}{0:>16d}{_fmt16(x)}{_fmt16(y)}\n'
        line2 = f'{"*":8s}{_fmt16(z)}\n'
        return line1 + line2

    def _ctria6(eid, pid, n):
        """Format CTRIA6* card (long format, 16-char fields)."""
        line1 = (f'{"CTRIA6*":8s}{eid:>16d}{pid:>16d}'
                 f'{n[0]:>16d}{n[1]:>16d}\n')
        line2 = (f'{"*":8s}{n[2]:>16d}{n[3]:>16d}'
                 f'{n[4]:>16d}{n[5]:>16d}\n')
        return line1 + line2

    with open(filename, 'w', encoding='utf-8') as f:
        f.write('$ Nastran BDF - Coil surface mesh (CTRIA6)\n')
        f.write('$ Generated by Radia BEM inductance panel\n')
        f.write('BEGIN BULK\n')
        f.write('$\n')
        f.write('$ Grid cards (long format)\n')
        f.write('$\n')

        for i, (x, y, z) in enumerate(nodes):
            f.write(_grid_star(i + 1, x, y, z))

        f.write('$\n')
        f.write('$ Element cards\n')
        f.write('$\n')

        for i, conn in enumerate(elements):
            n = [c + 1 for c in conn]  # 1-indexed
            f.write(_ctria6(i + 1, 1, n))

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


def _compute_B_field(mesh_surf, gf_J, elem_A, MU_0,
                     lx=0, ly=0, lz=0, maxh_vol=0):
    """Compute B field in air volume via direct Biot-Savart.

    Args:
        lx, ly, lz: Box half-sizes [m]. 0 = auto (conductor bbox + 30%).
        maxh_vol: Volume mesh element size [m]. 0 = auto (10% of extent).

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
    center = (bbox_min + bbox_max) / 2

    half = np.array([lx, ly, lz])
    box = Box(Pnt(*(center - half)), Pnt(*(center + half)))
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
    parser.add_argument("--mode", default="solve", choices=["solve", "post"],
                        help="solve: BEM solve only; post: B-field + GMSH export")
    # Solve mode args
    parser.add_argument("--cub5", default="", help="Cubit .cub5 model file (solve mode)")
    parser.add_argument("--order", type=int, default=2, help="Curve order (1-2)")
    parser.add_argument("--fes-order", type=int, default=0, help="HDivSurface order (0=RWG)")
    parser.add_argument("--source", default="source", help="Source block name")
    parser.add_argument("--sink", default="sink", help="Sink block name")
    parser.add_argument("--msh-output", default="", help="GMSH .msh output path")
    parser.add_argument("--output", default="", help="Output JSON file (optional)")
    # Workpiece ESIM/Dowell args
    parser.add_argument("--workpiece", default="", help="Workpiece block name")
    parser.add_argument("--impedance-model", default="esim",
                        choices=["esim", "dowell"], help="Surface impedance model")
    parser.add_argument("--frequency", type=float, default=50000, help="Frequency [Hz]")
    parser.add_argument("--sigma", type=float, default=2e6, help="Conductivity [S/m]")
    parser.add_argument("--half-thickness", type=float, default=0.005,
                        help="Slab half-thickness [m]")
    parser.add_argument("--material", default="steel",
                        choices=["steel", "copper", "aluminum"])
    # Post mode args
    parser.add_argument("--j-npy", default="", help="J_coeffs.npy path (post mode)")
    parser.add_argument("--mesh-vol", default="", help="surface_mesh.vol path (post mode)")
    parser.add_argument("--lx", type=float, default=0.07, help="Box half-size X [m]")
    parser.add_argument("--ly", type=float, default=0.07, help="Box half-size Y [m]")
    parser.add_argument("--lz", type=float, default=0.07, help="Box half-size Z [m]")
    parser.add_argument("--maxh-vol", type=float, default=0.01, help="Volume mesh size [m]")
    args = parser.parse_args()

    import io as _io
    real_stdout = sys.stdout
    sys.stdout = _io.StringIO()
    try:
        if args.mode == "solve":
            result = extract_inductance(args.cub5, args.order,
                                        args.source, args.sink,
                                        args.fes_order, args.msh_output,
                                        workpiece=args.workpiece,
                                        impedance_model=args.impedance_model,
                                        frequency=args.frequency,
                                        sigma=args.sigma,
                                        half_thickness=args.half_thickness,
                                        material=args.material)
        else:
            result = post_process(args.mesh_vol, args.fes_order,
                                  args.msh_output, args.j_npy,
                                  args.lx, args.ly, args.lz,
                                  args.maxh_vol)
    except Exception as e:
        result = {"error": str(e)}
    sys.stdout = real_stdout

    if args.output:
        with open(args.output, "w", encoding="utf-8") as f:
            json.dump(result, f)
    print(json.dumps(result))


if __name__ == "__main__":
    main()
