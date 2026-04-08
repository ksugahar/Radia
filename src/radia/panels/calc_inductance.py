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


_this_dir = os.path.dirname(os.path.abspath(__file__))
if _this_dir not in sys.path:
    sys.path.insert(0, _this_dir)

from calc_common import (setup_cubit as _setup_cubit_common, setup_paths, MU_0,
                          progress, calc_main, write_surface_only_vol)


def _setup_cubit():
    """Import cubit with path cleanup (NGSolve already imported)."""
    cubit = _setup_cubit_common()
    if cubit is None:
        raise RuntimeError("Cubit not available")
    return cubit


def _write_surface_only_vol(src_ngmesh, output_path):
    """Legacy wrapper -- delegates to calc_common.write_surface_only_vol."""
    write_surface_only_vol(src_ngmesh, output_path)


def extract_inductance_vol(vol_file, source_label="source",
                           sink_label="sink", fes_order=0, msh_output=""):
    """Extract self-inductance from .vol file (no Cubit needed).

    Args:
        vol_file: Path to Netgen .vol file (surface mesh with boundary labels)
        source_label: Boundary label for source face (sideset name)
        sink_label: Boundary label for sink face (sideset name)
        fes_order: HDivSurface basis order (0=RWG, 1+=higher-order)
        msh_output: Optional path for GMSH .msh output

    Returns:
        dict with inductance_H, diagnostics
    """
    import numpy as np
    import time as _time
    from ngsolve import Mesh, Integrate, CF, BND

    setup_paths()
    from bem_inductance import compute_inductance_source_sink

    t_total_start = _time.perf_counter()

    # Load mesh from .vol
    mesh = Mesh(vol_file)

    # BEM requires surface-only mesh (HDivSurface on volume mesh
    # includes interior edges as DOFs -> singular SL matrix).
    if mesh.ne > 0:
        import tempfile
        surf_vol = os.path.join(tempfile.gettempdir(), "bem_surface.vol")
        # Re-export as surface-only: read .vol, remove volume elements, save
        ng = mesh.ngmesh
        # Delete all volume elements by creating a fresh mesh
        # that only has surface elements, points, and face descriptors
        _write_surface_only_vol(ng, surf_vol)
        mesh = Mesh(surf_vol)

    t_export = _time.perf_counter() - t_total_start

    nse = mesh.GetNE(BND)
    nv = mesh.nv
    ne = sum(1 for e in mesh.edges)
    euler = nv - ne + nse
    area = float(Integrate(CF(1), mesh, VOL_or_BND=BND))

    # Verify boundary labels
    bnd_labels = list(set(mesh.GetBoundaries()))
    if source_label not in bnd_labels or sink_label not in bnd_labels:
        return {"error": f"Boundary labels {bnd_labels} do not contain "
                f"'{source_label}' and/or '{sink_label}'. "
                f"Define source/sink as sidesets in Cubit before export."}

    # Export mesh to GMSH (before solve)
    if msh_output:
        base_dir = os.path.dirname(os.path.abspath(msh_output))
        gmsh_file_J = msh_output
    else:
        base_dir = os.path.dirname(os.path.abspath(vol_file))
        gmsh_file_J = os.path.join(base_dir, "inductance_J.msh").replace("\\", "/")
    try:
        from gmsh_post_export import GmshPostExport
        post = GmshPostExport(mesh, boundary=True)
        post.write_mesh(gmsh_file_J)
        sys.stderr.write(f"MESH_READY:{gmsh_file_J}\n")
        sys.stderr.flush()
    except Exception:
        pass

    # Saddle point EFIE
    sol = compute_inductance_source_sink(mesh, source_label, sink_label, fes_order)
    if 'error' in sol:
        return sol

    L_total = sol['L']
    t_solve = sol['t_total']

    sys.stderr.write(f"SOLVE_DONE:{t_solve:.1f}s\n")
    sys.stderr.flush()

    # Save results
    j_npy = os.path.join(base_dir, "J_coeffs.npy").replace("\\", "/")
    np.save(j_npy, sol['J'])
    mesh_vol = os.path.join(base_dir, "surface_mesh.vol").replace("\\", "/")
    mesh.ngmesh.Save(mesh_vol)

    coords = np.array([(v.point[0], v.point[1], v.point[2])
                       for v in mesh.vertices])
    bbox_min, bbox_max = coords.min(axis=0), coords.max(axis=0)
    extent = bbox_max - bbox_min
    default_half = float(max(extent) * 0.65)
    default_maxh = float(max(extent) * 0.1)

    return {
        "inductance_H": float(L_total),
        "n_dofs": sol['n_J'],
        "n_faces": sol['n_f'],
        "euler": euler,
        "surface_area": area,
        "source_area": float(sol['A_source']),
        "sink_area": float(sol['A_sink']),
        "constraint_residual": sol['residual'],
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


def extract_inductance(cub5_file, order, source_label="source",
                       sink_label="sink", fes_order=0, msh_output="",
                       workpiece="", impedance_model="esim",
                       frequency=50000, sigma=2e6, half_thickness=0.005,
                       material="steel", esim_geometry="cylinder"):
    """Extract self-inductance via source/sink saddle point EFIE (legacy cub5 path).

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

    setup_paths()
    from bem_inductance import compute_inductance_source_sink

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
    from calc_common import export_mesh

    t0_export = _time.perf_counter()
    mesh = export_mesh(cubit, order=order, surface_only=True)
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
            half_thickness=half_thickness, material=material,
            esim_geometry=esim_geometry)
        result.update(wp_result)
        sys.stderr.write(f"ESIM_DONE:R={wp_result.get('wp_R_effective', 0):.4e}\n")
        sys.stderr.flush()

    return result


def _compute_workpiece_impedance(mesh_coil, gf_J, workpiece_label, cubit_mod,
                                  impedance_model="esim", frequency=50000,
                                  sigma=2e6, half_thickness=0.005,
                                  material="steel", esim_geometry="cylinder"):
    """Compute workpiece surface impedance via ESIM or Dowell.

    Pipeline:
      1. Extract workpiece surface panels from Cubit block geometry
      2. Biot-Savart from coil J -> H at workpiece panels
      3. ESIM cell problem or Dowell formula -> Z_s, P, Q per panel
      4. Integrate over surface

    Returns dict with wp_* keys to merge into main result.
    """
    import math
    import numpy as np
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
    from calc_common import create_esim_solver, get_bh_curve

    bh_curve, mu_r = get_bh_curve(material)

    P_total = 0.0
    Q_total = 0.0
    Z_sum = 0.0 + 0.0j
    delta_min = float('inf')
    delta_max = 0.0

    if impedance_model == "esim":
        solver = create_esim_solver(
            material=material, frequency=frequency, sigma=sigma,
            half_thickness=half_thickness, geometry=esim_geometry)

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

    elif impedance_model == "bem-sibc":
        return _compute_workpiece_bem_sibc(
            mesh_coil, gf_J, panels, cubit_mod, wp_vol_ids,
            frequency=frequency, sigma=sigma, half_thickness=half_thickness,
            material=material, esim_geometry=esim_geometry)

    elif impedance_model in ("fem-esim", "fem-dowell", "fem-sibc"):
        _zs_mode = {"fem-esim": "esim", "fem-dowell": "dowell",
                     "fem-sibc": "sibc"}[impedance_model]
        return _compute_workpiece_fem_esim(
            mesh_coil, gf_J, panels, cubit_mod, wp_vol_ids,
            frequency=frequency, sigma=sigma, half_thickness=half_thickness,
            material=material, esim_geometry=esim_geometry,
            zs_mode=_zs_mode)

    R_eff = float(Z_sum.real)
    X_eff = float(Z_sum.imag)
    wp_area = sum(p['area'] for p in panels)

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
        "wp_P_density": float(P_total / wp_area) if wp_area > 0 else 0,
        "wp_R_effective": R_eff,
        "wp_X_effective": X_eff,
        "wp_delta_min": float(delta_min),
        "wp_delta_max": float(delta_max),
    }


def _compute_workpiece_bem_sibc(mesh_coil, gf_J, panels, cubit_mod, wp_vol_ids,
                                 frequency=50000, sigma=2e6, half_thickness=0.005,
                                 material="steel", esim_geometry="cylinder"):
    """Compute workpiece impedance via Scalar BIE + SIBC (BEM).

    Uses ScalarBIESIBCSolver from bem_sibc_solver.py.
    Workpiece surface mesh auto-generated from Cubit geometry.
    phi_inc from coil computed via Biot-Savart (coil center line approximation).

    Returns dict with wp_* keys.
    """
    import math
    import time as _time
    from ngsolve import (Mesh, Integrate, CF, BND, TaskManager)
    from netgen.occ import Cylinder, Pnt, Dir, OCCGeometry, Glue
    from bem_sibc_solver import (ScalarBIESIBCSolver,
                                 compute_phi_inc_from_surface_J)
    from esim_cell_problem import ESIMFiniteSlabSolver

    omega = 2 * math.pi * frequency

    # --- Auto-detect workpiece geometry from Cubit surfaces ---
    centers = np.array([p['center'] for p in panels])
    wp_center = np.mean(centers, axis=0)
    z_min, z_max = centers[:, 2].min(), centers[:, 2].max()
    wp_height = z_max - z_min
    rho = np.sqrt((centers[:, 0] - wp_center[0])**2 +
                  (centers[:, 1] - wp_center[1])**2)
    wp_radius = float(np.max(rho))
    if wp_radius < 1e-6:
        wp_radius = half_thickness

    sys.stderr.write(f"BEM-SIBC: wp R={wp_radius*1e3:.1f}mm, "
                     f"H={wp_height*1e3:.1f}mm\n")
    sys.stderr.flush()

    # --- Extract per-element J from EFIE solution (actual surface current) ---
    from ngsolve import Integrate, CF, BND
    elem_A = Integrate(CF(1), mesh_coil, VOL_or_BND=BND, element_wise=True)
    elem_Jx = Integrate(gf_J[0], mesh_coil, VOL_or_BND=BND, element_wise=True)
    elem_Jy = Integrate(gf_J[1], mesh_coil, VOL_or_BND=BND, element_wise=True)
    elem_Jz = Integrate(gf_J[2], mesh_coil, VOL_or_BND=BND, element_wise=True)

    src_centroids, src_areas, src_J_vecs = [], [], []
    for el in mesh_coil.Elements(BND):
        area = abs(elem_A[el.nr])
        if area < 1e-30:
            continue
        jvec = np.array([elem_Jx[el.nr], elem_Jy[el.nr], elem_Jz[el.nr]]) / area
        verts = [mesh_coil.vertices[v.nr].point for v in el.vertices]
        c = np.mean([(v[0], v[1], v[2]) for v in verts], axis=0)
        src_centroids.append(c)
        src_areas.append(area)
        src_J_vecs.append(jvec)

    src_centroids = np.array(src_centroids)
    src_areas = np.array(src_areas)
    src_J_vecs = np.array(src_J_vecs)
    sys.stderr.write(f"BEM-SIBC: {len(src_areas)} coil elements for Biot-Savart\n")
    sys.stderr.flush()

    # --- Workpiece surface mesh ---
    # Try Cubit export first (ACIS curved, order=2), fallback to OCC
    t0 = _time.perf_counter()
    import tempfile as _tmpfile
    try:
        # Mesh workpiece volumes in Cubit if not already meshed
        for vid in wp_vol_ids:
            n_tets = len(cubit_mod.get_volume_tets(vid))
            n_hexes = len(cubit_mod.get_volume_hexes(vid))
            if n_tets + n_hexes == 0:
                cubit_mod.cmd(f'volume {vid} scheme tetmesh')
                cubit_mod.cmd(f'volume {vid} size {wp_radius / 4}')
                cubit_mod.cmd(f'mesh volume {vid}')
                sys.stderr.write(f"BEM-SIBC: auto-meshed wp volume {vid}\n")
                sys.stderr.flush()

        # Export via C++ plugin
        vol_path = _tmpfile.mktemp(suffix='.vol')
        cubit_mod.cmd(f'radia_export netgen "{vol_path}" order 2 overwrite')
        wp_mesh = Mesh(vol_path)
        sys.stderr.write(f"BEM-SIBC: Cubit export OK (ACIS order=2)\n")
        sys.stderr.flush()
    except Exception as e:
        sys.stderr.write(f"BEM-SIBC: Cubit export failed ({e}), OCC fallback\n")
        sys.stderr.flush()
        from netgen.occ import (Cylinder as _Cyl, Pnt as _P, Dir as _D,
                                 OCCGeometry as _G, Glue as _Gl)
        wp_cyl = _Cyl(_P(wp_center[0], wp_center[1], z_min),
                       _D(0, 0, 1), wp_radius, wp_height)
        for f in wp_cyl.faces:
            f.name = "wp_surface"
        wp_mesh = Mesh(_G(_Gl(wp_cyl.faces)).GenerateMesh(maxh=wp_radius / 4))
        wp_mesh.Curve(3)
    t_mesh = _time.perf_counter() - t0

    ne_wp = wp_mesh.GetNE(BND)
    sys.stderr.write(f"BEM-SIBC: wp mesh {ne_wp} elems ({t_mesh:.1f}s)\n")
    sys.stderr.flush()

    # --- Assemble BEM operators ---
    solver = ScalarBIESIBCSolver(wp_mesh, order=1)
    sys.stderr.write(f"BEM-SIBC: assembled {solver.ndof} DOFs ({solver.t_assembly:.1f}s)\n")
    sys.stderr.flush()

    # --- Compute phi_inc ---
    t0 = _time.perf_counter()
    node_coords = np.array([[wp_mesh.vertices[i].point[j] for j in range(3)]
                            for i in range(wp_mesh.nv)])

    # Use EFIE-solved surface J for phi_inc (accounts for proximity effects)
    phi_inc = compute_phi_inc_from_surface_J(
        node_coords, src_centroids, src_areas, src_J_vecs, n_quad=20)
    t_phi = _time.perf_counter() - t0

    sys.stderr.write(f"BEM-SIBC: phi_inc [{phi_inc.min():.4f}, {phi_inc.max():.4f}] "
                     f"({t_phi:.1f}s)\n")
    sys.stderr.flush()

    # --- ESIM + Karl iteration ---
    STEEL_BH = [
        [0.0, 0.0], [50.0, 0.10], [100.0, 0.25], [200.0, 0.55],
        [500.0, 0.95], [1000.0, 1.20], [2000.0, 1.40], [5000.0, 1.55],
        [10000.0, 1.65], [20000.0, 1.75], [50000.0, 1.90], [100000.0, 2.00],
    ]
    bh_curve = STEEL_BH if material == "steel" else None
    mu_r = 1.0 if material in ("copper", "aluminum") else None

    esim_solver = ESIMFiniteSlabSolver(
        half_thickness=half_thickness, bh_curve=bh_curve,
        sigma=sigma, frequency=frequency,
        mu_r=mu_r if bh_curve is None else None,
        n_nodes=200, geometry=esim_geometry)

    Z_s = esim_solver.solve(5.0)['Z']
    max_iter = 15
    tol = 1e-3
    relax = 0.5

    for iteration in range(max_iter):
        result = solver.solve(phi_inc, Z_s=Z_s, omega=omega)
        H_t_rms = result['H_t_rms']

        Z_s_old = Z_s
        sol = esim_solver.solve(max(H_t_rms, 1e-3))
        Z_s_new = sol['Z']
        Z_s = relax * Z_s_new + (1 - relax) * Z_s_old

        dZ = abs(Z_s - Z_s_old) / max(abs(Z_s_old), 1e-30)
        if dZ < tol and iteration > 0:
            break

    # Final solve with converged Z_s
    result = solver.solve(phi_inc, Z_s=Z_s, omega=omega)
    H_t_rms = result['H_t_rms']
    P_density = result['P_density']
    P_total = P_density * result['area']
    Q_density = 0.5 * Z_s.imag * H_t_rms**2 if Z_s != 0 else 0
    Q_total = Q_density * result['area']

    mu_eff = 4e-7 * math.pi * (mu_r if mu_r else 100.0)
    delta = math.sqrt(2.0 / (omega * mu_eff * sigma)) if omega > 0 else 1e10

    return {
        "wp_model": "bem-sibc",
        "wp_material": material,
        "wp_frequency": frequency,
        "wp_sigma": sigma,
        "wp_half_thickness": half_thickness,
        "wp_n_panels": ne_wp,
        "wp_H_t_max": float(H_t_rms),
        "wp_H_t_min": float(H_t_rms),
        "wp_P_total": float(P_total),
        "wp_Q_total": float(Q_total),
        "wp_P_density": float(P_density),
        "wp_R_effective": float(0.5 * Z_s.real * H_t_rms**2 * result['area']),
        "wp_X_effective": float(0.5 * Z_s.imag * H_t_rms**2 * result['area']),
        "wp_delta_min": float(delta),
        "wp_delta_max": float(delta),
        "wp_bem_ndof": solver.ndof,
        "wp_bem_t_assembly": solver.t_assembly,
        "wp_bem_n_coil_elems": len(src_areas),
        "wp_bem_wp_radius": wp_radius,
    }


def _compute_workpiece_fem_esim(mesh_coil, gf_J, panels, cubit_mod, wp_vol_ids,
                                 frequency=50000, sigma=2e6, half_thickness=0.005,
                                 material="steel", esim_geometry="cylinder",
                                 zs_mode="esim"):
    """Compute workpiece impedance via 3D FEM with SIBC Robin condition.

    Creates 3D volume mesh (air + coil + workpiece cavity),
    solves HCurl curl-curl + SIBC penalty, Karl iteration for Z_s.

    Returns dict with wp_* keys (same format as BEM-SIBC).
    """
    import math
    import time as _time
    from ngsolve import (Mesh, HCurl, BilinearForm, LinearForm, GridFunction,
                         Integrate, curl, dx, ds, CF, BND, VOL, TaskManager,
                         InnerProduct, x, y, z, sqrt, IfPos)
    from netgen.occ import (WorkPlane, Axes, Axis, Pnt, Dir, Sphere,
                             Cylinder, OCCGeometry, Glue)
    from esim_cell_problem import ESIMFiniteSlabSolver

    MU_0 = 4e-7 * np.pi
    omega = 2 * math.pi * frequency
    nu0 = 1.0 / MU_0

    # --- Auto-detect geometry from panels + coil mesh ---
    centers = np.array([p['center'] for p in panels])
    wp_center = np.mean(centers, axis=0)
    z_min, z_max = centers[:, 2].min(), centers[:, 2].max()
    wp_height = z_max - z_min
    rho = np.sqrt((centers[:, 0] - wp_center[0])**2 +
                  (centers[:, 1] - wp_center[1])**2)
    R_wp = float(np.max(rho))
    if R_wp < 1e-6:
        R_wp = half_thickness
    H_wp = wp_height if wp_height > 1e-6 else 2 * R_wp

    coil_coords = np.array([(v.point[0], v.point[1], v.point[2])
                            for v in mesh_coil.vertices])
    coil_rho = np.sqrt(coil_coords[:, 0]**2 + coil_coords[:, 1]**2)
    R_coil = float(np.mean(coil_rho))
    coil_extent = np.max(np.abs(coil_coords))
    a_coil = float(np.std(coil_rho)) * 2  # rough cross-section estimate

    sys.stderr.write(f"FEM-ESIM: coil R={R_coil*1e3:.1f}mm, "
                     f"wp R={R_wp*1e3:.1f}mm H={H_wp*1e3:.1f}mm\n")
    sys.stderr.flush()

    # --- 3D Geometry ---
    t0 = _time.perf_counter()
    gap_deg = 5
    wp_plane = WorkPlane(Axes(p=Pnt(R_coil, 0, 0),
                              n=Dir(0, 1, 0), h=Dir(0, 0, 1)))
    circle = wp_plane.Circle(a_coil).Face()
    torus = circle.Revolve(Axis(Pnt(0, 0, 0), Dir(0, 0, 1)), 360 - gap_deg)
    torus.name = "coil"
    torus.maxh = a_coil

    R_air = max(R_coil * 4, coil_extent * 3)
    air_sphere = Sphere(Pnt(0, 0, 0), R_air)
    air_sphere.name = "air"
    for f in air_sphere.faces:
        f.name = "outer"

    wp_cyl = Cylinder(Pnt(wp_center[0], wp_center[1], z_min),
                      Dir(0, 0, 1), R_wp, H_wp)
    for f in wp_cyl.faces:
        f.name = "wp_surface"
        f.maxh = min(R_wp / 4, a_coil)

    air = air_sphere - torus - wp_cyl
    air.name = "air"
    shape = Glue([air, torus])

    geo = OCCGeometry(shape)
    maxh_air = R_air / 8
    with TaskManager():
        ngmesh = geo.GenerateMesh(maxh=maxh_air, grading=0.3)
    mesh = Mesh(ngmesh)
    mesh.Curve(2)
    t_mesh = _time.perf_counter() - t0

    ne = mesh.ne
    bnds = mesh.GetBoundaries()
    has_wp = "wp_surface" in bnds

    sys.stderr.write(f"FEM-ESIM: {ne} tets, {mesh.nv} verts ({t_mesh:.1f}s)\n")
    sys.stderr.flush()

    if not has_wp:
        return {"wp_error": "wp_surface boundary not found in FEM mesh"}

    # --- FEM-SIBC solve with Karl iteration ---
    t0 = _time.perf_counter()
    fes = HCurl(mesh, order=1, dirichlet="outer", complex=True)
    u, v = fes.TnT()

    I_total = 1.0
    J0 = I_total / (math.pi * a_coil**2)
    r_xy = sqrt(x * x + y * y)
    r_safe = IfPos(r_xy - 1e-10, r_xy, 1e-10)
    J_source = J0 * CF((-y / r_safe, x / r_safe, 0))

    f_lf = LinearForm(fes)
    f_lf += J_source * v * dx("coil")
    f_lf.Assemble()

    wp_region = mesh.Boundaries("wp_surface")
    A_wp = Integrate(CF(1), mesh, BND, definedon=wp_region).real

    STEEL_BH = [
        [0, 0], [50, 0.1], [100, 0.25], [200, 0.55],
        [500, 0.95], [1000, 1.2], [2000, 1.4], [5000, 1.55],
        [10000, 1.65], [20000, 1.75], [50000, 1.9], [100000, 2.0],
    ]
    bh_curve = STEEL_BH if material == "steel" else None
    mu_r = 1.0 if material in ("copper", "aluminum") else None

    # Z_s computation: ESIM (nonlinear, Karl iteration), Dowell or SIBC (linear, no iteration)
    rho = 1.0 / sigma
    mu_eff = MU_0 * (mu_r if mu_r else 100.0)
    delta = math.sqrt(2.0 / (omega * mu_eff * sigma)) if omega > 0 else 1e10
    _label = f"FEM-{zs_mode.upper()}"
    esim_solver = None

    if zs_mode == "sibc":
        # Classical SIBC: Z_s = (1+j)/(sigma*delta), semi-infinite half-space
        Z_s = complex(1, 1) / (sigma * delta)
        max_iter = 1  # no iteration (linear, no thickness correction)
    elif zs_mode == "dowell":
        # Dowell: Z_s = (rho/a)*gamma_a*tanh(gamma_a), finite slab
        xi = half_thickness / delta
        gamma_a = complex(1, 1) * xi
        Z_s = (rho / half_thickness) * gamma_a * np.tanh(gamma_a)
        max_iter = 1  # no iteration (linear)
    else:
        # ESIM: 1D cell problem, supports nonlinear BH + curvature
        esim_solver = ESIMFiniteSlabSolver(
            half_thickness=R_wp, bh_curve=bh_curve, sigma=sigma,
            frequency=frequency,
            mu_r=mu_r if bh_curve is None else None, n_nodes=200,
            geometry=esim_geometry)
        Z_s = esim_solver.solve(5.0)['Z']
        max_iter = 15  # Karl iteration for nonlinear

    gfu = GridFunction(fes)
    tol = 1e-3
    relax = 0.5
    H_t_rms = 5.0
    P_total = 0.0
    Q_total = 0.0

    # Per-element Z_s on workpiece surface (for ESIM spatially varying Z_s)
    from ngsolve import SurfaceL2, NumberSpace
    fes_zs = SurfaceL2(mesh, order=0, complex=True)
    gf_zs_inv = GridFunction(fes_zs)  # stores 1/Z_s per element (for Robin coeff)
    # Initialize with uniform Z_s
    gf_zs_inv.Set(1.0 / Z_s, definedon=wp_region)

    sys.stderr.write(f"{_label}: {fes.ndof} DOFs, solving...\n")
    sys.stderr.flush()

    for iteration in range(max_iter):
        a_bf = BilinearForm(fes)
        a_bf += nu0 * curl(u) * curl(v) * dx
        # SIBC Robin: (jw/Z_s)*A_t*v_t, with Z_s varying per element
        if esim_solver:
            sibc_cf = 1j * omega * gf_zs_inv  # per-element
        else:
            sibc_cf = 1j * omega / Z_s  # uniform scalar
        a_bf += sibc_cf * u.Trace() * v.Trace() * ds("wp_surface")
        a_bf.Assemble()

        with TaskManager():
            gfu.vec.data = a_bf.mat.Inverse(
                fes.FreeDofs(), inverse="pardiso") * f_lf.vec

        if esim_solver:
            # Per-element H_t and Z_s update
            At_sq = sum(gfu[i].real * gfu[i].real +
                        gfu[i].imag * gfu[i].imag for i in range(3))
            elem_At2 = Integrate(At_sq, mesh, BND,
                                 definedon=wp_region, element_wise=True)
            elem_area = Integrate(CF(1), mesh, BND,
                                  definedon=wp_region, element_wise=True)

            Z_s_old_arr = gf_zs_inv.vec.FV().NumPy().copy()
            P_total = 0.0
            Q_total = 0.0
            max_dZ = 0.0
            n_wp_elems = 0

            for el in mesh.Elements(BND):
                if mesh.GetBCName(el.index) != "wp_surface":
                    continue
                area_e = abs(elem_area[el.nr])
                if area_e < 1e-30:
                    continue
                At2_e = abs(elem_At2[el.nr])
                At_rms_e = math.sqrt(At2_e / area_e)
                # H_t from SIBC: H_t = |jw/Z_s|*|A_t|
                Zs_old_e = 1.0 / Z_s_old_arr[el.nr] if abs(Z_s_old_arr[el.nr]) > 1e-30 else Z_s
                H_t_e = abs(1j * omega / Zs_old_e) * At_rms_e

                sol = esim_solver.solve(max(float(H_t_e), 1e-3))
                Zs_new = sol['Z']
                Zs_e = relax * Zs_new + (1 - relax) * Zs_old_e
                gf_zs_inv.vec[el.nr] = 1.0 / Zs_e

                P_total += sol['P_prime'] * area_e
                Q_total += sol['Q_prime'] * area_e
                dZ_e = abs(Zs_e - Zs_old_e) / max(abs(Zs_old_e), 1e-30)
                max_dZ = max(max_dZ, dZ_e)
                n_wp_elems += 1

            dZ = max_dZ
            # Compute overall H_t_rms
            At_int = Integrate(At_sq, mesh, BND, definedon=wp_region).real
            At_rms = math.sqrt(max(At_int, 0) / max(A_wp, 1e-30))
            # Mean |Z_s| for reporting
            Zs_mean = A_wp / max(abs(np.sum(gf_zs_inv.vec.FV().NumPy()[:n_wp_elems])), 1e-30)
            H_t_rms = abs(1j * omega * np.mean(np.abs(
                gf_zs_inv.vec.FV().NumPy()[:fes_zs.ndof]))) * At_rms
        else:
            # Uniform Z_s (Dowell or SIBC): single scalar
            At_sq = sum(gfu[i].real * gfu[i].real +
                        gfu[i].imag * gfu[i].imag for i in range(3))
            At_int = Integrate(At_sq, mesh, BND, definedon=wp_region).real
            At_rms = math.sqrt(max(At_int, 0) / max(A_wp, 1e-30))
            H_t_rms = abs(1j * omega / Z_s) * At_rms
            P_total = 0.5 * Z_s.real * H_t_rms**2 * A_wp
            Q_total = 0.5 * Z_s.imag * H_t_rms**2 * A_wp
            dZ = 0.0

        sys.stderr.write(f"{_label}: iter {iteration} "
                         f"H_t={H_t_rms:.2f} P={P_total:.4e} dZ={dZ:.4e}\n")
        sys.stderr.flush()
        if dZ < tol and iteration > 0:
            break

    t_solve = _time.perf_counter() - t0

    return {
        "wp_model": _label.lower(),
        "wp_material": material,
        "wp_frequency": frequency,
        "wp_sigma": sigma,
        "wp_half_thickness": half_thickness,
        "wp_n_panels": ne,
        "wp_H_t_max": float(H_t_rms),
        "wp_H_t_min": float(H_t_rms),
        "wp_P_total": float(P_total),
        "wp_Q_total": float(Q_total),
        "wp_P_density": float(P_total / A_wp) if A_wp > 0 else 0,
        "wp_R_effective": float(P_total),
        "wp_X_effective": float(Q_total),
        "wp_delta_min": float(delta),
        "wp_delta_max": float(delta),
        "wp_fem_ndof": fes.ndof,
        "wp_fem_ne": ne,
        "wp_fem_t_mesh": round(t_mesh, 1),
        "wp_fem_t_solve": round(t_solve, 1),
    }


def _compute_B_field_cf(mesh_surf, gf_J, elem_A, lx=0, ly=0, lz=0, maxh_vol=0):
    """Compute B field CoefficientFunction via BiotSavartCF.

    Converts EFIE-solved surface current J into wire segments,
    then builds a BiotSavartCF (multipole-accelerated) for fast evaluation
    on any volume mesh.

    Args:
        mesh_surf: Surface mesh (from BEM solve)
        gf_J: HDivSurface GridFunction (solved J)
        elem_A: Per-element areas
        lx, ly, lz: Box half-sizes [m] for volume mesh
        maxh_vol: Volume mesh element size [m]

    Returns:
        (B_cf, mesh_vol): CoefficientFunction for B [T] and volume mesh
    """
    import numpy as np
    import math
    from ngsolve import Mesh, BND, Integrate
    from netgen.occ import Box, Pnt, OCCGeometry
    from netgen.meshing import MeshingParameters

    setup_paths()
    from biot_savart import biot_savart_wire

    # Extract per-element J and geometry
    elem_Jx = Integrate(gf_J[0], mesh_surf, VOL_or_BND=BND, element_wise=True)
    elem_Jy = Integrate(gf_J[1], mesh_surf, VOL_or_BND=BND, element_wise=True)
    elem_Jz = Integrate(gf_J[2], mesh_surf, VOL_or_BND=BND, element_wise=True)

    # Convert surface current elements to wire segments:
    # Each triangle with J -> short filament at centroid with current I = |J|*A/(2*eps)
    # Filament endpoints: centroid +/- eps * J_hat, length = 2*eps
    # This preserves the magnetic dipole moment: m = I * area = J * A * eps
    segments = []
    currents = []
    for el in mesh_surf.Elements(BND):
        area = abs(elem_A[el.nr])
        if area < 1e-30:
            continue
        jvec = np.array([elem_Jx[el.nr], elem_Jy[el.nr], elem_Jz[el.nr]]) / area
        j_mag = np.linalg.norm(jvec)
        if j_mag < 1e-30:
            continue

        verts = [mesh_surf.vertices[v.nr].point for v in el.vertices]
        centroid = np.mean([(v[0], v[1], v[2]) for v in verts], axis=0)

        # Filament length ~ element size
        eps = math.sqrt(area) * 0.5
        j_hat = jvec / j_mag
        p1 = centroid - eps * j_hat
        p2 = centroid + eps * j_hat
        I_fil = j_mag * area / (2 * eps)

        segments.append((tuple(p1), tuple(p2)))
        currents.append(I_fil)

    progress("BIOT_SAVART", f"{len(segments)} filaments from surface J")

    # Build BiotSavartCF with per-segment current
    from ngsolve.bem import BiotSavartCF
    from ngsolve.bla import Vec3D

    # Compute center and radius from all segment endpoints
    all_pts = np.array([p for seg in segments for p in seg])
    center = np.mean(all_pts, axis=0)
    dists = np.linalg.norm(all_pts - center[np.newaxis, :], axis=1)
    rad = float(max(np.max(dists) * 1.5, 1e-10))

    bs = BiotSavartCF(order=15, kappa=1e-10,
                      center=Vec3D(*center), rad=rad)
    for (p1, p2), I in zip(segments, currents):
        bs.AddCurrent(Vec3D(*p1), Vec3D(*p2), complex(I), 4)

    # BiotSavartCF returns complex CF; take real part for DC problems
    B_cf = MU_0 * bs.real  # H [A/m] -> B [T], real part only

    # Create volume mesh for visualization
    coords = np.array([(v.point[0], v.point[1], v.point[2])
                       for v in mesh_surf.vertices])
    bbox_min, bbox_max = coords.min(axis=0), coords.max(axis=0)
    vol_center = (bbox_min + bbox_max) / 2
    half = np.array([lx, ly, lz])

    box = Box(Pnt(*(vol_center - half)), Pnt(*(vol_center + half)))
    mesh_vol = Mesh(OCCGeometry(box).GenerateMesh(
        mp=MeshingParameters(maxh=maxh_vol)))

    progress("BIOT_SAVART", f"vol mesh {mesh_vol.nv} verts")
    return B_cf, mesh_vol


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
    from ngsolve import (Mesh, Integrate, CF, BND, HDivSurface, GridFunction,
                         Norm, HDiv)
    from netgen.meshing import Mesh as NetgenMesh

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

    progress("FIELD_READY", "J computed")

    # B-distribution via BiotSavartCF + GridFunction.Set()
    B_cf, mesh_vol = _compute_B_field_cf(
        mesh, gf_J, elem_A, lx=lx, ly=ly, lz=lz, maxh_vol=maxh_vol)

    progress("B_FIELD_READY", "B computed")

    # Write output via GmshPostExport
    gmsh_file_B = os.path.join(base_dir, "inductance_B.msh").replace("\\", "/")
    gmsh_file_J = os.path.join(base_dir, "inductance_J.msh").replace("\\", "/")

    from gmsh_post_export import GmshPostExport

    # Volume B-field via BiotSavartCF -> GridFunction.Set()
    fes_B = HDiv(mesh_vol, order=1)
    gf_B = GridFunction(fes_B)
    gf_B.Set(B_cf)
    post_B = GmshPostExport(mesh_vol)
    post_B.add_vector_field("B", gf_B)
    post_B.add_scalar_field("|B|", Norm(gf_B))
    post_B.write(gmsh_file_B)
    progress("GMSH", f"B written: {gmsh_file_B}")

    # Surface J-field (high-order curved elements)
    post_J = GmshPostExport(mesh, boundary=True)
    post_J.add_field("|J|", np.linalg.norm(J_nodes, axis=1), ncomp=1)
    post_J.add_field("J", J_nodes, ncomp=3)
    post_J.write(gmsh_file_J)
    progress("GMSH", f"J written: {gmsh_file_J}")

    # Companion .geo (Merge both files)
    geo_file = os.path.join(base_dir, "inductance.geo").replace("\\", "/")
    with open(geo_file, 'w', encoding='utf-8') as f:
        f.write(f'Merge "{os.path.basename(gmsh_file_B)}";\n')
        f.write(f'Merge "{os.path.basename(gmsh_file_J)}";\n')
        f.write('Mesh.NumSubEdges = 4;\n')
        f.write('Mesh.VolumeEdges = 0;\n')

    # COMSOL-compatible text interpolation files
    vol_nodes_list = [(v.point[0], v.point[1], v.point[2])
                      for v in mesh_vol.vertices]
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




def main():
    parser = argparse.ArgumentParser(
        description="BEM inductance (source/sink saddle point EFIE)")
    parser.add_argument("--mode", default="solve", choices=["solve", "post"],
                        help="solve: BEM solve only; post: B-field + GMSH export")
    # Solve mode args
    parser.add_argument("--vol", default="", help="Netgen .vol file (preferred, no Cubit needed)")
    parser.add_argument("--cub5", default="", help="Cubit .cub5 model file (legacy)")
    parser.add_argument("--order", type=int, default=2, help="Curve order (1-2)")
    parser.add_argument("--fes-order", type=int, default=0, help="HDivSurface order (0=RWG)")
    parser.add_argument("--source", default="source", help="Source block name")
    parser.add_argument("--sink", default="sink", help="Sink block name")
    parser.add_argument("--msh-output", default="", help="GMSH .msh output path")
    parser.add_argument("--output", default="", help="Output JSON file (optional)")
    # Workpiece ESIM/Dowell args
    parser.add_argument("--workpiece", default="", help="Workpiece block name")
    parser.add_argument("--impedance-model", default="esim",
                        choices=["esim", "dowell", "bem-sibc",
                                 "fem-esim", "fem-dowell", "fem-sibc"],
                        help="esim/dowell=BEM per-panel, bem-sibc=BEM scalar BIE, "
                             "fem-esim/fem-dowell/fem-sibc=3D FEM with Z_s model")
    parser.add_argument("--frequency", type=float, default=50000, help="Frequency [Hz]")
    parser.add_argument("--sigma", type=float, default=2e6, help="Conductivity [S/m]")
    parser.add_argument("--half-thickness", type=float, default=0.005,
                        help="Slab half-thickness [m]")
    parser.add_argument("--material", default="steel",
                        choices=["steel", "copper", "aluminum"])
    parser.add_argument("--mu-r", type=float, default=0,
                        help="Relative permeability for SIBC/Dowell (0=auto from material)")
    parser.add_argument("--esim-geometry", default="local_curvature",
                        choices=["local_curvature", "none"],
                        help="ESIM curvature: local_curvature=Bessel (R=half_thickness), "
                             "none=flat slab (cosh/sinh)")
    parser.add_argument("--bh-file", default="",
                        help="BH curve file (2-column: H[A/m] B[T], overrides --material)")
    # Post mode args
    parser.add_argument("--j-npy", default="", help="J_coeffs.npy path (post mode)")
    parser.add_argument("--mesh-vol", default="", help="surface_mesh.vol path (post mode)")
    parser.add_argument("--lx", type=float, default=0.07, help="Box half-size X [m]")
    parser.add_argument("--ly", type=float, default=0.07, help="Box half-size Y [m]")
    parser.add_argument("--lz", type=float, default=0.07, help="Box half-size Z [m]")
    parser.add_argument("--maxh-vol", type=float, default=0.01, help="Volume mesh size [m]")

    def run(args):
        if args.mode == "solve":
            # Prefer --vol (no Cubit needed)
            if args.vol:
                return extract_inductance_vol(args.vol,
                                              args.source, args.sink,
                                              args.fes_order, args.msh_output)
            elif args.cub5:
                _geom_map = {"local_curvature": "cylinder", "none": "slab"}
                esim_geom = _geom_map.get(args.esim_geometry, "cylinder")
                bh_file = getattr(args, 'bh_file', '')
                if bh_file and os.path.isfile(bh_file):
                    import numpy as _np
                    bh_data = _np.loadtxt(bh_file).tolist()
                    material_override = "steel"
                else:
                    bh_data = None
                    material_override = args.material

                return extract_inductance(args.cub5, args.order,
                                          args.source, args.sink,
                                          args.fes_order, args.msh_output,
                                          workpiece=args.workpiece,
                                          impedance_model=args.impedance_model,
                                          frequency=args.frequency,
                                          sigma=args.sigma,
                                          half_thickness=args.half_thickness,
                                          material=material_override,
                                          esim_geometry=esim_geom)
            else:
                return {"error": "Either --vol or --cub5 is required"}
        else:
            return post_process(args.mesh_vol, args.fes_order,
                                args.msh_output, args.j_npy,
                                args.lx, args.ly, args.lz,
                                args.maxh_vol)

    calc_main(run, parser)


if __name__ == "__main__":
    main()
