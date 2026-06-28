"""
BEM-SIBC workpiece eddy current calculator.

Called as subprocess from Cubit panel:
    python calc_heating_bem.py --vol workpiece.vol
        --coil-radius 0.030 --frequency 7000 --sigma 2e6

Pipeline:
  1. Load surface mesh from .vol (auto-extract BND from volume mesh)
  2. Assemble BEM operators (ScalarBIESIBCSolver)
  3. Compute phi_inc from Biot-Savart coil (filamentary loop)
  4. Karl iteration (nonlinear Z_s via ESIM)
  5. Return JSON with P_total, H_t_rms, Z_s

No coil mesh needed -- coil is analytical (Biot-Savart).

IMPORTANT: NGSolve must be imported BEFORE cubit.
Outputs JSON to stdout (suppresses all other print output).
"""

import argparse
import math
import os
import sys

import numpy as np

_this_dir = os.path.dirname(os.path.abspath(__file__))
_repo = os.path.abspath(os.path.join(_this_dir, "..", "..", ".."))
_src = os.path.join(_repo, "src")
_panels = os.path.join(_src, "radia", "panels")
for _p in (_this_dir, _src, _panels):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from calc_common import (setup_paths, MU_0, progress, calc_main,
                          create_esim_solver,
                          EMMaterial, add_material_args)


def _extract_surface_mesh_filtered(vol_mesh, keep_label="",
                                    return_vertex_map=False):
    """Extract a clean 2D surface mesh containing only boundary elements
    that touch the specified material.

    Args:
        vol_mesh: NGSolve Mesh (3D volume mesh)
        keep_label: material name to keep (e.g. "coil"). If empty, keeps all.
        return_vertex_map: if True, also return the
            ``new_to_old`` mapping (extracted-mesh vertex nr -> volume-mesh
            vertex nr) so the caller can re-evaluate per-vertex fields on
            the parent (curved) volume mesh.

    Returns:
        NGSolve Mesh (surface only, no orphan vertices, all boundary labels
        renumbered consecutively). If ``return_vertex_map`` is True a
        ``(mesh, new_to_old)`` tuple is returned instead.
    """
    from ngsolve import Mesh, BND, VOL
    import netgen.meshing as ngm

    # Determine which volume index corresponds to keep_label.
    # An empty keep_label means "no filter" — extract all boundaries.
    # A non-empty keep_label that does not match any material is an
    # error: the caller asked us to filter by a label that does not
    # exist, and silently switching to all-boundaries would produce a
    # mesh containing the air sphere outer surface and other unrelated
    # boundaries that ruin the BEM solve.
    keep_dom = 0  # 0 = no filter
    materials = vol_mesh.GetMaterials()
    if keep_label:
        for i, m in enumerate(materials, 1):
            if m == keep_label:
                keep_dom = i
                break
        if keep_dom == 0:
            raise ValueError(
                f"Surface mesh extractor: requested keep_label "
                f"{keep_label!r} is not in the .vol's materials list "
                f"({sorted(set(materials))}). Fix the .jou (block name) "
                f"or pass an existing material name.")

    ngmesh_new = ngm.Mesh(dim=3)
    bnd_labels = list(vol_mesh.GetBoundaries())

    # Pre-scan: which boundary labels actually have elements adjacent to keep_dom?
    used_labels = []
    label_to_fd = {}
    for el in vol_mesh.Elements(BND):
        if keep_dom > 0:
            # Check FaceDescriptor's domin/domout for keep_dom adjacency
            fd = vol_mesh.ngmesh.FaceDescriptor(el.index + 1)
            if fd.domin != keep_dom and fd.domout != keep_dom:
                continue
        lbl = bnd_labels[el.index]
        if lbl not in label_to_fd:
            new_idx = len(used_labels) + 1
            fd_new = ngm.FaceDescriptor(bc=new_idx)
            fd_idx = ngmesh_new.Add(fd_new)
            ngmesh_new.SetBCName(new_idx - 1, lbl)
            label_to_fd[lbl] = fd_idx
            used_labels.append(lbl)

    if not used_labels:
        # Material exists but has no adjacent boundary elements. The
        # .vol is malformed (volume tagged but no surface elements).
        raise ValueError(
            f"Surface mesh extractor: material {keep_label!r} has no "
            f"adjacent boundary elements. The .vol export is "
            f"inconsistent — re-run the Cubit .jou and re-export.")

    # Add boundary elements (and their vertices) for the filtered set.
    old_to_new = {}
    for el in vol_mesh.Elements(BND):
        if keep_dom > 0:
            fd = vol_mesh.ngmesh.FaceDescriptor(el.index + 1)
            if fd.domin != keep_dom and fd.domout != keep_dom:
                continue
        lbl = bnd_labels[el.index]
        fd_idx = label_to_fd.get(lbl)
        if fd_idx is None:
            continue
        new_verts = []
        for v in el.vertices:
            if v.nr not in old_to_new:
                pt = vol_mesh.vertices[v.nr].point
                old_to_new[v.nr] = ngmesh_new.Add(
                    ngm.MeshPoint(ngm.Pnt(pt[0], pt[1], pt[2])))
            new_verts.append(old_to_new[v.nr])
        se = ngm.Element2D(fd_idx, new_verts)
        ngmesh_new.Add(se)

    surf_mesh = Mesh(ngmesh_new)
    if return_vertex_map:
        # ngmesh.Add returns a netgen ``PointId`` (1-indexed). The
        # NGSolve Mesh exposes vertices via 0-indexed ``v.nr``.
        # PointId is NOT directly castable to int — use its ``.nr``
        # attribute. (Discovered the hard way 2026-04-12 when the
        # original ``int(new_id)`` cast crashed BEM with TypeError.)
        new_to_old = {new_id.nr - 1: int(old_nr)
                       for old_nr, new_id in old_to_new.items()}
        return surf_mesh, new_to_old
    return surf_mesh


def _extract_surface_mesh(vol_mesh, order=2):
    """Extract a clean 2D surface mesh from a 3D volume mesh.

    Creates a new mesh containing only boundary surface elements
    and their vertices (no interior points). This avoids singular
    mass matrices when using H1 on the surface.

    Args:
        vol_mesh: NGSolve Mesh (3D volume mesh)
        order: Curve order for the surface mesh

    Returns:
        NGSolve Mesh (2D surface-only)
    """
    from ngsolve import Mesh, BND
    import netgen.meshing as ngm

    ngmesh_new = ngm.Mesh(dim=3)

    # Collect boundary vertices and remap indices
    old_to_new = {}

    bnd_labels = list(vol_mesh.GetBoundaries())

    # Add face descriptors for each boundary label
    label_to_fd = {}
    seen_labels = []
    for lbl in bnd_labels:
        if lbl not in label_to_fd:
            fd = ngm.FaceDescriptor(bc=len(seen_labels) + 1)
            fd_idx = ngmesh_new.Add(fd)
            ngmesh_new.SetBCName(len(seen_labels), lbl)
            label_to_fd[lbl] = fd_idx
            seen_labels.append(lbl)

    for el in vol_mesh.Elements(BND):
        verts = el.vertices
        new_verts = []
        for v in verts:
            if v.nr not in old_to_new:
                pt = vol_mesh.vertices[v.nr].point
                old_to_new[v.nr] = ngmesh_new.Add(
                    ngm.MeshPoint(ngm.Pnt(pt[0], pt[1], pt[2])))
            new_verts.append(old_to_new[v.nr])

        # Get face descriptor index for this element's boundary label
        el_label = bnd_labels[el.index]
        fd_idx = label_to_fd.get(el_label, 1)

        se = ngm.Element2D(fd_idx, new_verts)
        ngmesh_new.Add(se)

    mesh = Mesh(ngmesh_new)
    with TaskManager():
        mesh.Curve(order)
        return mesh


def compute_heating_bem(vol_file, coil_radius=0.030, coil_current=1.0,
                        gap_deg=5, frequency=7000, mat=None,
                        h1_order=1, wp_label="wp_surface",
                        half_thickness=0.010, esim_geometry="cylinder",
                        impedance_model="esim",
                        max_iter=15, tol=1e-3,
                        msh_output="",
                        coil_vol="", coil_source="source", coil_sink="sink",
                        coil_label="coil"):
    """BEM-SIBC workpiece eddy current from .vol file.

    Args:
        vol_file: Path to Netgen .vol (volume or surface mesh)
        coil_radius: Coil loop radius [m]
        coil_current: Total coil current [A]
        gap_deg: Coil gap angle [degrees] (0 = closed loop)
        frequency: Operating frequency [Hz]
        mat: EMMaterial instance (workpiece properties)
        h1_order: H1 polynomial order on surface
        wp_label: Boundary label for workpiece surface
        half_thickness: Workpiece radius for ESIM [m]
        esim_geometry: "cylinder" or "planar"
        max_iter: Karl iteration maximum
        tol: Karl convergence tolerance
        msh_output: Optional GMSH output path

    Returns:
        dict with P_total, H_t_rms, Z_s, etc.
    """
    if mat is None:
        mat = EMMaterial.from_name("steel")
    sigma = mat.sigma
    mu_r = mat.mu_r
    material = mat.name

    import time as _time
    from ngsolve import Mesh, Integrate, CF, BND

    setup_paths()
    from radia.bem_sibc_solver import ScalarBIESIBCSolver, compute_phi_inc_from_loop

    omega = 2 * np.pi * frequency
    t_total_start = _time.perf_counter()

    # === 1. Load mesh and extract workpiece surface ===
    progress("MESH", "Loading .vol file...")
    mesh_full = Mesh(vol_file)

    if mesh_full.ne > 0:
        # Try to extract the workpiece volume only — keeps the BEM problem
        # to the conducting body, not the entire mesh boundary.
        if wp_label and wp_label in mesh_full.GetMaterials():
            mesh = _extract_surface_mesh_filtered(mesh_full,
                                                  keep_label=wp_label)
        else:
            mesh = _extract_surface_mesh(mesh_full, order=2)
    else:
        mesh = mesh_full

    nse = mesh.GetNE(BND)
    nv = mesh.nv
    area = float(Integrate(CF(1), mesh, BND))
    t_mesh = _time.perf_counter() - t_total_start

    progress("MESH", f"{nse} surface elements, {nv} vertices, "
             f"area={area:.4e} m^2 ({t_mesh:.1f}s)")

    # Verify boundary label exists
    bnd_labels = list(set(mesh.GetBoundaries()))
    if wp_label and wp_label not in bnd_labels:
        sys.stderr.write(f"WARNING: boundary '{wp_label}' not found in "
                         f"{bnd_labels}. Using all surfaces.\n")
        sys.stderr.flush()

    # === 2. Assemble BEM operators ===
    progress("BEM_ASSEMBLY", "Assembling SL, DL, M, K...")
    solver = ScalarBIESIBCSolver(mesh, order=h1_order)
    progress("BEM_ASSEMBLY", f"{solver.ndof} DOFs, {solver.t_assembly:.1f}s")

    # === 3. Compute phi_inc from Biot-Savart coil ===
    progress("PHI_INC", "Computing phi_inc (Biot-Savart)...")
    t0 = _time.perf_counter()

    node_coords = np.array([[mesh.vertices[i].point[j] for j in range(3)]
                            for i in range(mesh.nv)])

    if coil_vol:
        # --- Mesh coil path: BEM EFIE for J, then phi_inc from surface J ---
        from radia.bem_sibc_solver import compute_phi_inc_from_surface_J
        from radia.bem_inductance import compute_inductance_source_sink
        from ngsolve import Mesh as NGSolveMesh
        from ngsolve import Integrate, CF, BND

        progress("COIL_BEM", f"Loading coil mesh: {coil_vol}")
        coil_mesh_full = NGSolveMesh(coil_vol)
        # Extract clean surface mesh (no orphan vertices) — coil only.
        if coil_mesh_full.ne > 0:
            coil_mesh = _extract_surface_mesh_filtered(
                coil_mesh_full, keep_label=coil_label)
        else:
            coil_mesh = coil_mesh_full

        progress("COIL_BEM",
                 f"coil surface: nse={coil_mesh.GetNE(BND)}, nv={coil_mesh.nv}")

        # Solve EFIE saddle point system for J on coil surface
        sol_coil = compute_inductance_source_sink(
            coil_mesh, coil_source, coil_sink, fes_order=0)
        if 'error' in sol_coil:
            return {"error": f"Coil BEM failed: {sol_coil['error']}"}

        gf_J = sol_coil['gf_J']
        L_coil = sol_coil['L']
        progress("COIL_BEM",
                 f"L_coil={L_coil*1e9:.2f} nH, n_J={sol_coil['n_J']}")

        # Extract per-element J vector and area for compute_phi_inc_from_surface_J
        elem_A = Integrate(CF(1), coil_mesh, VOL_or_BND=BND, element_wise=True)
        elem_Jx = Integrate(gf_J[0], coil_mesh, VOL_or_BND=BND, element_wise=True)
        elem_Jy = Integrate(gf_J[1], coil_mesh, VOL_or_BND=BND, element_wise=True)
        elem_Jz = Integrate(gf_J[2], coil_mesh, VOL_or_BND=BND, element_wise=True)

        centroids, areas, Jvecs = [], [], []
        for el in coil_mesh.Elements(BND):
            a = abs(elem_A[el.nr])
            if a < 1e-30:
                continue
            jvec = np.array([elem_Jx[el.nr], elem_Jy[el.nr],
                             elem_Jz[el.nr]]) / a
            verts = [coil_mesh.vertices[v.nr].point for v in el.vertices]
            c = np.mean([(v[0], v[1], v[2]) for v in verts], axis=0)
            centroids.append(c)
            areas.append(a)
            Jvecs.append(jvec)

        centroids = np.array(centroids)
        areas = np.array(areas)
        Jvecs = np.array(Jvecs)
        # Scale to physical current
        Jvecs *= coil_current

        progress("PHI_INC", f"computing from {len(centroids)} surface panels")
        phi_inc_nodes = compute_phi_inc_from_surface_J(
            node_coords, centroids, areas, Jvecs, n_quad=20)
    else:
        # --- Analytical loop path (original) ---
        phi_inc_nodes = compute_phi_inc_from_loop(
            node_coords, loop_center=[0, 0, 0], loop_radius=coil_radius,
            current=coil_current, n_quad=30, gap_deg=gap_deg)

    t_phi = _time.perf_counter() - t0
    progress("PHI_INC", f"range=[{phi_inc_nodes.min():.4f}, "
             f"{phi_inc_nodes.max():.4f}] ({t_phi:.1f}s)")

    # === 4. Surface impedance + Karl iteration ===
    # Determine mu_r for skin depth reporting
    if mu_r is None:
        mu_r_report = 100.0 if material == "steel" else 1.0
    else:
        mu_r_report = mu_r
    mu_eff = MU_0 * mu_r_report
    delta = math.sqrt(2.0 / (omega * mu_eff * sigma)) if omega * sigma > 0 else float('inf')

    if impedance_model == "linear":
        # Linear SIBC: Z_s = (1+j) * rho / delta, fixed (no iteration)
        rho = 1.0 / sigma
        mu_lin = MU_0 * (mu_r if mu_r is not None else mu_r_report)
        delta_lin = math.sqrt(2 * rho / (omega * mu_lin)) if omega > 0 else float('inf')
        Z_s = complex(1, 1) * rho / delta_lin
        delta = delta_lin

        progress("SIBC", f"Linear Z_s={Z_s:.4e}, mu_r={mu_r_report}, "
                 f"delta={delta*1e3:.4f}mm")

        result = solver.solve(phi_inc_nodes, Z_s=Z_s, omega=omega)
        H_t_rms = result['H_t_rms']
        P_density = result['P_density']
        P_total = P_density * result['area']
        n_converged = 1
        mu_r_eff = mu_r_report

        progress("SOLVE_DONE", f"P={P_total:.4e}W, H_t_rms={H_t_rms:.2f}A/m")

    else:
        # ESIM: nonlinear Z_s(H_t), Karl iteration
        progress("KARL", "Starting Karl iteration...")

        esim_solver = mat.create_esim_solver(
            frequency, half_thickness, geometry=esim_geometry)

        H_t_init = 5.0
        Z_s = esim_solver.solve(H_t_init)['Z']
        relax = 0.5

        n_converged = 0
        for iteration in range(max_iter):
            result = solver.solve(phi_inc_nodes, Z_s=Z_s, omega=omega)
            H_t_rms = result['H_t_rms']
            P_density = result['P_density']
            P_total = P_density * result['area']

            Z_s_old = Z_s
            sol = esim_solver.solve(max(H_t_rms, 1e-3))
            Z_s_new = sol['Z']
            Z_s = relax * Z_s_new + (1 - relax) * Z_s_old

            dZ = abs(Z_s - Z_s_old) / max(abs(Z_s_old), 1e-30)
            progress("KARL", f"iter {iteration}: |Z_s|={abs(Z_s):.4e}, "
                     f"H_t_rms={H_t_rms:.2f}, P={P_total:.4e}W, dZ/Z={dZ:.4e}")

            if dZ < tol and iteration > 0:
                n_converged = iteration
                break
        else:
            n_converged = max_iter

        # Final result with converged Z_s
        result = solver.solve(phi_inc_nodes, Z_s=Z_s, omega=omega)
        H_t_rms = result['H_t_rms']
        P_density = result['P_density']
        P_total = P_density * result['area']

        # Effective mu_r from ESIM (for reporting)
        esim_info = esim_solver.solve(max(H_t_rms, 1e-3))
        mu_final = esim_info.get('mu_final', MU_0 * mu_r_report)
        mu_r_eff = float(abs(mu_final) / MU_0) if mu_final else mu_r_report

        t_total_check = _time.perf_counter() - t_total_start
        progress("SOLVE_DONE", f"P={P_total:.4e}W, H_t_rms={H_t_rms:.2f}A/m "
                 f"({t_total_check:.1f}s)")

    Q_density = 0.5 * Z_s.imag * H_t_rms**2 if Z_s != 0 else 0
    Q_total = Q_density * result['area']

    t_total = _time.perf_counter() - t_total_start

    progress("SOLVE_DONE", f"P={P_total:.4e}W, H_t_rms={H_t_rms:.2f}A/m "
             f"({t_total:.1f}s)")

    # === 4b. Save .vol + .sol (workpiece mesh + heating density) ===
    base_dir = os.path.dirname(os.path.abspath(vol_file))
    if msh_output:
        base_dir = os.path.dirname(os.path.abspath(msh_output))
    try:
        from ngsolve import GridFunction, grad, Norm, H1

        # q = sigma * |H_t|^2 / 2 per vertex.
        # Compute per-element |grad phi|^2 via Integrate, then vertex-average.
        phi_vec = result['phi_vec']
        gf_phi_re = GridFunction(solver.fes)
        gf_phi_im = GridFunction(solver.fes)
        gf_phi_re.vec.FV().NumPy()[:] = phi_vec.real
        gf_phi_im.vec.FV().NumPy()[:] = phi_vec.imag

        H_t_sq_cf = Norm(grad(gf_phi_re))**2 + Norm(grad(gf_phi_im))**2
        elem_area = Integrate(CF(1), solver.mesh, BND, element_wise=True)
        elem_Ht2 = Integrate(H_t_sq_cf, solver.mesh, BND, element_wise=True)

        # Per-vertex q [W/m^2] via element averaging.
        #
        # SIBC heating density: q'' = 0.5 * Re(Z_s) * |H_t|^2
        # (= 0.5 * |J_s|^2 / (sigma * delta) since J_s = n × H_t).
        # Re(Z_s) comes directly from the solver's nonlinear ESIM
        # (varies with |H_t|) or from the linear sibc fallback, matching
        # the GMSH P_loss viz formula below at line ~505. The previous
        # `0.5 * sigma * |H_t|^2` was dimensionally W/m^3 (volumetric)
        # not W/m^2 — saved values were ~10^8 too large (sugahara 2026-
        # 04-15 GUI feedback, same root cause as calc_inductance v4.5.9).
        Re_Zs = float(Z_s.real)
        q_arr = np.zeros(solver.mesh.nv)
        nc = np.zeros(solver.mesh.nv)
        for el in solver.mesh.Elements(BND):
            a = max(abs(float(elem_area[el.nr])), 1e-30)
            ht2 = float(elem_Ht2[el.nr]) / a
            qe = 0.5 * Re_Zs * ht2
            for v in el.vertices:
                q_arr[v.nr] += qe
                nc[v.nr] += 1
        q_arr = np.where(nc > 0, q_arr / nc, 0.0)

        # Phase B: compute + save the four scalar fields displayed in GMSH
        # as individual .vol + .sol pairs, then vol2msh() to combine.
        # Values all use the same formula and Re_Zs so the saved file
        # exactly matches what GMSH renders.
        Im_Zs = float(Z_s.imag)

        def _elemwise_cf_to_vertex_array(cf):
            """Integrate `cf` per element, divide by area, vertex-average."""
            arr = np.zeros(solver.mesh.nv)
            nc_local = np.zeros(solver.mesh.nv)
            elem_val = Integrate(cf, solver.mesh, BND, element_wise=True)
            for el in solver.mesh.Elements(BND):
                a = max(abs(float(elem_area[el.nr])), 1e-30)
                v_elem = float(elem_val[el.nr]) / a
                for vtx in el.vertices:
                    arr[vtx.nr] += v_elem
                    nc_local[vtx.nr] += 1
            return np.where(nc_local > 0, arr / nc_local, 0.0)

        from ngsolve import sqrt
        from ngsolve import TaskManager
        H_t_cf = sqrt(H_t_sq_cf)
        P_cf = 0.5 * Re_Zs * H_t_sq_cf
        Q_cf = 0.5 * Im_Zs * H_t_sq_cf

        Ht_arr = _elemwise_cf_to_vertex_array(H_t_cf)
        P_arr = _elemwise_cf_to_vertex_array(P_cf)
        Q_arr = _elemwise_cf_to_vertex_array(Q_cf)

        # All four fields live on solver.mesh (the BEM surface). Save as
        # H1 order=1 scalar GridFunctions so vol2msh can reload them.
        from gmsh_post_export import save_vol_sol_pair, vol2msh

        def _save_scalar(arr, sol_path):
            gf = GridFunction(H1(solver.mesh, order=1))
            gf.vec.FV().NumPy()[:len(arr)] = arr
            gf.Save(sol_path)

        wp_vol_path = os.path.join(base_dir, "workpiece.vol").replace("\\", "/")
        q_sol_path = os.path.join(base_dir, "q_heating.sol").replace("\\", "/")
        ht_sol_path = os.path.join(base_dir, "H_t.sol").replace("\\", "/")
        p_sol_path = os.path.join(base_dir, "P_loss.sol").replace("\\", "/")
        q_react_sol_path = os.path.join(base_dir, "Q_reactive.sol").replace("\\", "/")

        gf_q = GridFunction(H1(solver.mesh, order=1))
        gf_q.vec.FV().NumPy()[:len(q_arr)] = q_arr
        save_vol_sol_pair(wp_vol_path, q_sol_path,
                           solver.mesh.ngmesh, gf_q)
        _save_scalar(Ht_arr, ht_sol_path)
        _save_scalar(P_arr, p_sol_path)
        _save_scalar(Q_arr, q_react_sol_path)
        progress("SAVE", "workpiece.vol + q_heating.sol + H_t.sol + "
                          "P_loss.sol + Q_reactive.sol")
    except Exception as e:
        wp_vol_path = ""
        q_sol_path = ""
        sys.stderr.write(f"SOL save failed: {e}\n")
        sys.stderr.flush()

    # === 5. GMSH output via the shared vol2msh converter ===
    gmsh_file = ""
    if msh_output and wp_vol_path:
        try:
            vol2msh(msh_output, wp_vol_path, [
                {"sol": ht_sol_path, "fes": "H1", "fes_order": 1,
                 "name": "|H_t| [A/m]", "ncomp": 1},
                {"sol": p_sol_path, "fes": "H1", "fes_order": 1,
                 "name": "P_loss [W/m^2]", "ncomp": 1},
                {"sol": q_react_sol_path, "fes": "H1", "fes_order": 1,
                 "name": "Q_reactive [var/m^2]", "ncomp": 1},
            ])
            gmsh_file = msh_output
            progress("MSH", gmsh_file)
        except Exception as e:
            sys.stderr.write(f"GMSH export failed: {e}\n")
            sys.stderr.flush()

    return {
        "P_total_W": float(P_total),
        "Q_total_var": float(Q_total),
        "H_t_rms_Am": float(H_t_rms),
        "Z_s_abs_Ohm": float(abs(Z_s)),
        "Z_s_phase_deg": float(np.degrees(np.angle(Z_s))),
        "Z_s_real_Ohm": float(Z_s.real),
        "Z_s_imag_Ohm": float(Z_s.imag),
        "skin_depth_mm": float(delta * 1e3),
        "ndof": solver.ndof,
        "n_elements": nse,
        "area_m2": float(area),
        "n_iter": n_converged,
        "frequency_Hz": frequency,
        "wp_sigma_Sm": sigma,
        "wp_mu_r": float(mu_r_eff),
        "material": material,
        "coil_radius_m": coil_radius,
        "coil_current_A": coil_current,
        "coil_sigma_Sm": 5.8e7,  # copper (filamentary, not used in BEM)
        "t_mesh_s": round(t_mesh, 2),
        "t_assembly_s": round(solver.t_assembly, 2),
        "t_phi_inc_s": round(t_phi, 2),
        "t_total_s": round(t_total, 2),
        "gmsh_file": gmsh_file,
        "msh_output": gmsh_file,
    }


def main():
    parser = argparse.ArgumentParser(
        description="BEM-SIBC workpiece eddy current calculator")
    parser.add_argument("--vol", required=True, help="Path to .vol file")
    parser.add_argument("--coil-radius", type=float, default=0.030,
                        help="Coil loop radius [m]")
    parser.add_argument("--coil-current", type=float, default=1.0,
                        help="Coil current [A]")
    parser.add_argument("--gap-deg", type=float, default=5,
                        help="Coil gap angle [degrees]")
    parser.add_argument("--frequency", type=float, default=7000,
                        help="Frequency [Hz]")
    add_material_args(parser, include_custom=False)
    parser.add_argument("--h1-order", type=int, default=1,
                        help="H1 polynomial order")
    parser.add_argument("--wp-label", default="wp_surface",
                        help="Boundary label for workpiece surface")
    parser.add_argument("--half-thickness", type=float, default=0.010,
                        help="Workpiece radius for ESIM [m]")
    parser.add_argument("--esim-geometry", default="cylinder",
                        choices=["cylinder", "planar"])
    parser.add_argument("--impedance-model", default="esim",
                        choices=["esim", "linear"],
                        help="esim (nonlinear BH) or linear (fixed mu_r)")
    parser.add_argument("--max-iter", type=int, default=15,
                        help="Karl iteration max")
    parser.add_argument("--tol", type=float, default=1e-3,
                        help="Karl convergence tolerance")
    parser.add_argument("--msh-output", default="",
                        help="GMSH output path")
    parser.add_argument("--output", default="",
                        help="JSON output file")
    parser.add_argument("--coil-vol", default="",
                        help="Optional coil mesh .vol file. When set, the "
                             "coil current is solved from a BEM EFIE "
                             "(source/sink saddle point) on the coil mesh "
                             "and used as the source for phi_inc instead "
                             "of the analytical loop. Coil shape can be "
                             "arbitrary.")
    parser.add_argument("--coil-source", default="source",
                        help="Source boundary label on coil mesh "
                             "(only used with --coil-vol)")
    parser.add_argument("--coil-sink", default="sink",
                        help="Sink boundary label on coil mesh")
    parser.add_argument("--coil-label", default="coil",
                        help="Material name of the coil volume "
                             "(used to extract coil-only surface mesh)")

    def run(args):
        return compute_heating_bem(
            vol_file=args.vol,
            coil_radius=args.coil_radius,
            coil_current=args.coil_current,
            gap_deg=args.gap_deg,
            frequency=args.frequency,
            mat=EMMaterial.from_args(args),
            impedance_model=args.impedance_model,
            h1_order=args.h1_order,
            wp_label=args.wp_label,
            half_thickness=args.half_thickness,
            esim_geometry=args.esim_geometry,
            max_iter=args.max_iter,
            tol=args.tol,
            msh_output=args.msh_output,
            coil_vol=args.coil_vol,
            coil_source=args.coil_source,
            coil_sink=args.coil_sink,
            coil_label=args.coil_label,
        )

    calc_main(run, parser)


if __name__ == "__main__":
    main()
