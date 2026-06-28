"""
Inductance extractor using source/sink saddle point EFIE + ngsolve.bem.

Usage:
    python calc_inductance.py --vol model.vol --source source --sink sink

Input: Netgen .vol file with boundary labels (prefix-matched):
  - source*: current injection face(s)
  - sink*:   current extraction face(s)
  - Optional: coil material label for surface filtering
  - Optional: workpiece boundary labels for impedance calculation

Method (saddle point EFIE):
  1. Load .vol, extract surface mesh (filter by coil material if set)
  2. Build SL (LaplaceSL) and D (divergence) matrices
  3. Solve saddle point: [SL D^T; D 0] [J; p] = [0; g]
     where g encodes unit current injection at source/sink
  4. L = mu_0 * J^T @ SL @ J
  5. Optional: workpiece impedance via Biot-Savart + ESIM/Dowell

Outputs JSON to stdout.
"""

import argparse
import json
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
                          write_surface_only_vol,
                          EMMaterial, add_material_args)
from gmsh_post_export import write_combined_msh_v41 as _write_combined_msh_v41


def _resolve_bnd_labels(query, bnd_labels):
    """Look up *query* in *bnd_labels*. Exact match only.

    No case-insensitive matching, no prefix matching, no fuzzy fallback.
    The .jou contractually defines the labels; if the exact name is not
    present, error out so the user fixes the .jou rather than getting a
    silently-wrong result.
    """
    return [query] if query in bnd_labels else []


def extract_inductance_vol(vol_file, source_label="source",
                           sink_label="sink", fes_order=0, msh_output="",
                           coil_label="",
                           workpiece="", impedance_model="esim",
                           frequency=50000, mat=None,
                           half_thickness=0.005,
                           esim_geometry="cylinder",
                           coil_sigma=0, solver="lu",
                           field_air=False, use_local_curvature=False):
    """Extract self-inductance from .vol file (no Cubit needed).

    Args:
        vol_file: Path to Netgen .vol file (surface mesh with boundary labels)
        source_label: Boundary label for source face (prefix match)
        sink_label: Boundary label for sink face (prefix match)
        fes_order: HDivSurface basis order (0=RWG, 1+=higher-order)
        msh_output: Optional path for GMSH .msh output
        coil_label: Material name of the coil block.  When set, only surface
            elements adjacent to this material are kept (filters out
            workpiece/air surfaces so BEM solves on coil only).
        workpiece: Workpiece boundary label keyword (prefix match).
            Empty string disables workpiece computation.
        impedance_model: 'esim' or 'dowell'
        frequency: Operating frequency [Hz]
        mat: EMMaterial instance (workpiece properties: sigma, mu_r, BH)
        half_thickness: Slab half-thickness [m]
        esim_geometry: 'cylinder' (local curvature) or 'slab' (flat)
        coil_sigma: Coil conductivity [S/m] (stored in results, not used in BEM)

    Returns:
        dict with inductance_H, diagnostics, and wp_* keys if workpiece enabled
    """
    if mat is None:
        mat = EMMaterial.from_name("steel")
    sigma = mat.sigma
    mu_r = mat.mu_r

    import numpy as np
    import time as _time
    from ngsolve import Mesh, Integrate, CF, BND

    setup_paths()
    from radia.bem_inductance import compute_inductance_source_sink

    t_total_start = _time.perf_counter()

    # Load mesh from .vol (keep full mesh for workpiece panel extraction)
    mesh_full = Mesh(vol_file)
    # Activate high-order curving from the .vol's curvedelements
    # section. The Cubit exporter writes order=1..5; calling Curve(2)
    # is enough to make GmshPostExport detect the mesh as curved and
    # emit Tri6 elements with correct mid-edge nodes for the coil
    # visualization. NGSolve clips to the actual stored order, so
    # this is safe even when the .vol is order=1.
    try:
        mesh_full.Curve(2)
    except Exception:
        pass

    # BEM requires surface-only mesh (HDivSurface on volume mesh
    # includes interior edges as DOFs -> singular SL matrix).
    # Use the in-memory extractor that prunes orphan vertices and
    # renumbers FaceDescriptors. The .vol-text path
    # (write_surface_only_vol) preserves all 7000+ vertices from the
    # parent volume mesh, leaving thousands of unused vertices and
    # making HDivSurface DOFs / SL singular.
    if mesh_full.ne > 0:
        from calc_heating_bem import _extract_surface_mesh_filtered
        mesh, surf_new_to_old = _extract_surface_mesh_filtered(
            mesh_full, keep_label=coil_label, return_vertex_map=True)
    else:
        mesh = mesh_full
        surf_new_to_old = {i: i for i in range(mesh.nv)}

    t_export = _time.perf_counter() - t_total_start

    nse = mesh.GetNE(BND)
    nv = mesh.nv
    ne = sum(1 for e in mesh.edges)
    euler = nv - ne + nse
    area = float(Integrate(CF(1), mesh, VOL_or_BND=BND))

    # Resolve source/sink labels by keyword matching
    bnd_labels = list(set(mesh.GetBoundaries()))
    src_matches = _resolve_bnd_labels(source_label, bnd_labels)
    sink_matches = _resolve_bnd_labels(sink_label, bnd_labels)
    if not src_matches or not sink_matches:
        missing = []
        if not src_matches:
            missing.append(f"source='{source_label}'")
        if not sink_matches:
            missing.append(f"sink='{sink_label}'")
        return {"error":
                f"No boundary labels match {', '.join(missing)}. "
                f"Available boundary labels: {sorted(bnd_labels)[:20]}. "
                f"Fix: in your Cubit .jou, define sidesets named 'source' "
                f"and 'sink' on the coil terminal faces, e.g.:\n"
                f"  sideset 1 add surface <source_face_id>\n"
                f"  sideset 1 name \"source\"\n"
                f"  sideset 2 add surface <sink_face_id>\n"
                f"  sideset 2 name \"sink\"\n"
                f"Then re-export with: export netgen \"model.vol\" order N"}
    source_label = "|".join(src_matches)
    sink_label = "|".join(sink_matches)

    # Export coil surface mesh (no field yet) to GMSH so the user can
    # at least see the mesh while the BEM solve runs. The J / B fields
    # are added later (or in a separate post pass) — see field_air.
    # We use ``mesh_full`` (with Curve(2) activated above) so the
    # exporter can emit Tri6 elements with the correct mid-edge nodes
    # — the linear extracted mesh ``mesh`` would only render flat
    # triangles. The exporter writes BND elements grouped by
    # material, so the user can hide non-coil entities in GMSH.
    if msh_output:
        base_dir = os.path.dirname(os.path.abspath(msh_output))
        gmsh_file_J = msh_output
    else:
        base_dir = os.path.dirname(os.path.abspath(vol_file))
        gmsh_file_J = os.path.join(base_dir, "inductance_J.msh").replace("\\", "/")
    gmsh_file_J = gmsh_file_J.replace("\\", "/")
    try:
        from gmsh_post_export import GmshPostExport
        post = GmshPostExport(mesh_full, boundary=True)
        post.write_mesh(gmsh_file_J)
        sys.stderr.write(f"MESH_READY:{gmsh_file_J}\n")
        sys.stderr.flush()
    except Exception:
        gmsh_file_J = ""

    # Saddle point EFIE
    sol = compute_inductance_source_sink(
        mesh, source_label, sink_label, fes_order, solver=solver)
    if 'error' in sol:
        return sol

    L_total = sol['L']
    t_solve = sol['t_total']

    sys.stderr.write(f"SOLVE_DONE:{t_solve:.1f}s\n")
    sys.stderr.flush()

    # Save results (.vol + .sol — NGSolve standard pair)
    progress("POST_SOLVE", "saving surface_mesh.vol + J_coeffs.sol")
    mesh_vol = os.path.join(base_dir, "surface_mesh.vol").replace("\\", "/")
    mesh.ngmesh.Save(mesh_vol)

    from ngsolve import HDivSurface, GridFunction
    fes_save = HDivSurface(mesh, order=fes_order)
    gf_save = GridFunction(fes_save)
    gf_save.vec.FV().NumPy()[:] = sol['J']
    j_sol = os.path.join(base_dir, "J_coeffs.sol").replace("\\", "/")
    gf_save.Save(j_sol)

    # Legacy .npy (backward compat with existing post_process callers)
    j_npy = os.path.join(base_dir, "J_coeffs.npy").replace("\\", "/")
    np.save(j_npy, sol['J'])

    # Add the actual J vector field to the surface .msh so the GMSH
    # button has something meaningful to show even when --field-air
    # is off. We re-export the same file (overwrites the mesh-only
    # version written before the solve) with the J VECTOR view —
    # |J| is intentionally not exported because GMSH already shows
    # the magnitude when the user picks the vector view.
    #
    # The J coefficients live on the linear extracted mesh, so we
    # vertex-average the per-element J there and then map the result
    # back onto mesh_full's parent vertices via surf_new_to_old. The
    # write step uses mesh_full so the curved Tri6 elements survive.
    # J_nodes_full is also reused later by the workpiece coupled
    # path to re-write the same .msh with material-restricted views.
    J_nodes_full = None
    if gmsh_file_J:
        progress("POST_SOLVE", "writing coil J vector field to .msh")
        try:
            from gmsh_post_export import GmshPostExport
            from ngsolve import (Integrate, CF, BND, HDivSurface,
                                 GridFunction)
            fes_J_post = HDivSurface(mesh, order=fes_order)
            gf_J_post = GridFunction(fes_J_post)
            gf_J_post.vec.FV().NumPy()[:] = sol['J']

            elem_A = Integrate(CF(1), mesh, VOL_or_BND=BND,
                               element_wise=True)
            elem_Jx = Integrate(gf_J_post[0], mesh, VOL_or_BND=BND,
                                element_wise=True)
            elem_Jy = Integrate(gf_J_post[1], mesh, VOL_or_BND=BND,
                                element_wise=True)
            elem_Jz = Integrate(gf_J_post[2], mesh, VOL_or_BND=BND,
                                element_wise=True)

            J_nodes_extracted = np.zeros((mesh.nv, 3))
            nc_count = np.zeros(mesh.nv)
            for el in mesh.Elements(BND):
                a = max(abs(elem_A[el.nr]), 1e-30)
                jvec = (elem_Jx[el.nr] / a,
                        elem_Jy[el.nr] / a,
                        elem_Jz[el.nr] / a)
                for vtx in el.vertices:
                    J_nodes_extracted[vtx.nr, 0] += jvec[0]
                    J_nodes_extracted[vtx.nr, 1] += jvec[1]
                    J_nodes_extracted[vtx.nr, 2] += jvec[2]
                    nc_count[vtx.nr] += 1
            for k in range(3):
                J_nodes_extracted[:, k] = np.where(
                    nc_count > 0,
                    J_nodes_extracted[:, k] / nc_count, 0.0)

            # Lift J onto mesh_full's vertices via the new->old map.
            # Vertices not adjacent to a coil boundary stay zero.
            J_nodes_full = np.zeros((mesh_full.nv, 3))
            for new_nr, old_nr in surf_new_to_old.items():
                if 0 <= new_nr < mesh.nv and 0 <= old_nr < mesh_full.nv:
                    J_nodes_full[old_nr] = J_nodes_extracted[new_nr]

            post_J = GmshPostExport(mesh_full, boundary=True)
            post_J.add_field("J", J_nodes_full, ncomp=3)
            post_J.write(gmsh_file_J)
        except Exception as _e:
            sys.stderr.write(f"WARN: J post export failed: {_e}\n")
            sys.stderr.flush()

    coords = np.array([(v.point[0], v.point[1], v.point[2])
                       for v in mesh.vertices])
    bbox_min, bbox_max = coords.min(axis=0), coords.max(axis=0)
    extent = bbox_max - bbox_min
    default_half = float(max(extent) * 0.65)
    default_maxh = float(max(extent) * 0.1)

    result = {
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
        # GMSH button payload: the coil surface mesh file written
        # before the solve. The field_air path overwrites this with
        # the merged .geo (B + J + air-field) when --field-air is set.
        "msh_output": gmsh_file_J,
        "default_lxyz": round(default_half, 4),
        "default_maxh": round(default_maxh, 4),
    }
    if coil_sigma > 0:
        result["coil_sigma_Sm"] = coil_sigma

    # Workpiece (optional). Two paths:
    #   - dowell -> CoupledBEMSolver (iterative back-reaction with
    #     per-DOF f_back[i] = int v_i . A_wp dS_coil; gives the
    #     physically correct sign of Delta_L: <0 for non-magnetic
    #     conductors with strong screening, >0 for ferromagnetic
    #     workpieces). See bem_coupled_solver.py.
    #   - esim   -> one-way SIBC post-processing (nonlinear cell
    #     problem; the coupled solver does not yet support nonlinear
    #     Z_s, so this path keeps the legacy uncoupled estimator with
    #     P/Q reported but no Delta L).
    if workpiece:
        # "sibc" is the GUI name for analytical surface-impedance with
        # back-reaction. Internally it routes through _run_coupled_bem
        # (which uses the analytical Dowell Z_s and writes wp_J / wp_q
        # to the .msh for visualization). "dowell" is the legacy alias.
        if impedance_model == "sibc":
            impedance_model = "dowell"
        progress("WORKPIECE", f"impedance_model={impedance_model}, "
                              f"workpiece={workpiece}")
        if impedance_model == "dowell":
            try:
                coupled_result = _run_coupled_bem(
                    mesh_full=mesh_full,
                    workpiece_label=workpiece,
                    source_label=source_label,
                    sink_label=sink_label,
                    fes_order=fes_order,
                    frequency=frequency,
                    mat=mat,
                    half_thickness=half_thickness,
                    use_local_curvature=use_local_curvature,
                    output_dir=base_dir)
                # Pull out private wp_J/wp_q arrays before merging the
                # rest of the dict (they are NumPy arrays and would
                # break json.dumps later).
                wp_J_full_arr = coupled_result.pop("_wp_J_full", None)
                wp_q_full_arr = coupled_result.pop("_wp_q_full", None)
                result.update(coupled_result)
                # Headline inductance becomes the coupled L_total.
                # Keep the air-only value as L_air_H so the panel can
                # display both.
                result["L_air_H"] = result["inductance_H"]
                result["inductance_H"] = coupled_result[
                    "coupled_L_total_H"]
                # Re-write the J .msh with workpiece views appended
                # (J vector + heating density), restricted to the
                # "workpiece" material so they only render on the
                # workpiece surface.
                if (gmsh_file_J and wp_J_full_arr is not None
                        and J_nodes_full is not None):
                    try:
                        from gmsh_post_export import GmshPostExport
                        post_J2 = GmshPostExport(mesh_full, boundary=True)
                        post_J2.add_field(
                            "J_coil", J_nodes_full, ncomp=3,
                            material="coil")
                        post_J2.add_field(
                            "J_workpiece", wp_J_full_arr, ncomp=3,
                            material="workpiece")
                        post_J2.add_field(
                            "q_workpiece [W/m^2]", wp_q_full_arr, ncomp=1,
                            material="workpiece")
                        post_J2.write(gmsh_file_J)
                    except Exception as _e:
                        sys.stderr.write(
                            f"WARN: wp J/q export failed: {_e}\n")
                        sys.stderr.flush()
            except Exception as _exc:
                # Coupled solve failed — surface the error rather than
                # silently falling back to a one-way (wrong) result.
                import traceback as _tb
                result["coupled_error"] = (
                    f"CoupledBEMSolver failed: {_exc}\n{_tb.format_exc()}")
        else:  # esim (nonlinear cell problem)
            wp_panels = _extract_panels_from_mesh(mesh_full, workpiece)
            if not wp_panels:
                result["wp_error"] = (
                    f"No boundary labels starting with '{workpiece}' "
                    f"found. Available: "
                    f"{sorted(set(mesh_full.GetBoundaries()))[:20]}...")
            else:
                wp_result = _compute_wp_impedance_from_panels(
                    mesh, sol['gf_J'], wp_panels,
                    impedance_model=impedance_model, frequency=frequency,
                    mat=mat, half_thickness=half_thickness,
                    esim_geometry=esim_geometry,
                    mesh_wp=mesh_full,
                    use_local_curvature=use_local_curvature)
                result.update(wp_result)
                wp_P = wp_result.get("wp_P_total", 0.0)
                omega = 2.0 * math.pi * frequency
                if omega > 0:
                    result["coupled_R_effective_Ohm"] = float(2.0 * wp_P)

    # Optional: B-field post-processing in the air domain.
    # Reuses the post_process() pipeline (BiotSavartCF -> volume mesh
    # via Set() -> GMSH .msh + .geo). Runs inline after solve so the
    # user gets one-shot inductance + field in a single GUI Run.
    if field_air:
        progress("FIELD", "Computing B-field in air domain...")
        # default_lxyz/default_maxh come from the solve step's bbox
        lxyz = result.get("default_lxyz", 0.05)
        maxh = result.get("default_maxh", 0.01)
        post_result = post_process(
            mesh_vol_path=result["mesh_vol"],
            fes_order=fes_order,
            msh_output=result.get("msh_output", msh_output),
            j_npy=result["j_npy"],
            lx=lxyz, ly=lxyz, lz=lxyz, maxh_vol=maxh)
        # Merge field-eval keys with `field_` prefix to keep them
        # separate from inductance results.
        for k, v in post_result.items():
            result[f"field_{k}"] = v

    return result


def _dowell_Zs(R, rho, omega, mu_eff):
    """Linear SIBC Z_s from the Dowell tanh formula.

    Z_s = (rho / R) * gamma * tanh(gamma)
    gamma = (1 + j) * R / delta
    delta = sqrt(2 * rho / (omega * mu_eff))
    """
    if omega <= 0:
        return complex(0.0)
    delta = math.sqrt(2.0 * rho / (omega * mu_eff))
    xi = R / delta
    gamma_a = complex(1, 1) * xi
    # tanh saturates to 1 at |arg| > ~20 (double precision); avoid the
    # NumPy RuntimeWarning "overflow encountered in tanh" that otherwise
    # floods the panel Output window. PEC-limit Z_s = (rho/R) * gamma_a.
    if xi > 20.0:
        return (rho / R) * gamma_a
    return (rho / R) * gamma_a * np.tanh(gamma_a)


def _build_per_node_Zs(mesh_wp, panels, R_panel, rho, omega, mu_eff,
                        ndof_phi, fes=None):
    """Build a per-node Z_s array on the workpiece H1 mesh.

    For each workpiece panel, compute Z_s_panel from the Dowell tanh
    formula at the panel's local radius R_panel[i]. Project the
    per-panel values onto the H1 DOFs.

    Two projection paths:

      1. **Order = 1** (the production case): vertex averaging is
         exact for piecewise-linear data. The H1 ordering is
         vertex-first so node index == vertex index. This is the
         fast path.

      2. **Order > 1** (or fes supplied): assemble a piecewise-constant
         L2 source from the per-panel Z_s values, then **L2-project**
         onto the H1 space via ``GridFunction.Set`` of a
         CoefficientFunction built per material region. This properly
         distributes the per-panel data to all H1 DOFs (vertex + edge
         + face) instead of leaving the higher-order DOFs at the
         panel-mean fallback. Slower but correct for any order.

    The default ``fes=None`` keeps the legacy P1 vertex averaging
    behaviour for backward compatibility with the production
    coupled BEM solver where order=1 is the only supported case.

    Returns ``(Zs_node, Zs_panel)`` where ``Zs_node`` has shape
    ``(ndof_phi,)`` and ``Zs_panel`` has shape ``(n_panels,)``.
    """
    Zs_panel = np.array(
        [_dowell_Zs(float(R_panel[i]), rho, omega, mu_eff)
         for i in range(len(panels))], dtype=complex)

    # Path 2: order > 1 — L2-project a piecewise-constant CF onto fes
    # via SurfaceL2 -> H1 mass-matrix solve. This handles edge / face
    # DOFs correctly for any polynomial order.
    if fes is not None and getattr(fes, 'globalorder', 1) > 1:
        try:
            from ngsolve import (SurfaceL2, GridFunction, BilinearForm,
                                 LinearForm, ds, BND, TaskManager)

            # Build SurfaceL2 of the same dimensionality, attach the
            # per-panel Z_s as element-wise data via a GridFunction.
            fes_l2 = SurfaceL2(mesh_wp, order=0)
            gf_l2_re = GridFunction(fes_l2)
            gf_l2_im = GridFunction(fes_l2)
            re_vec = gf_l2_re.vec.FV().NumPy()
            im_vec = gf_l2_im.vec.FV().NumPy()
            re_vec[:] = 0
            im_vec[:] = 0
            # SurfaceL2 ordering matches mesh.Elements(BND) iteration
            # which is the same iteration order used by
            # _extract_panels_from_mesh; we have to map our panels list
            # back to the SurfaceL2 element index. Use el.nr.
            from ngsolve import BND as _BND
            el_nr_to_panel = {}
            j = 0
            for el in mesh_wp.Elements(_BND):
                el_nr_to_panel[el.nr] = j
                j += 1
            for el in mesh_wp.Elements(_BND):
                # Find the panel index that corresponds to this element
                pi = el_nr_to_panel.get(el.nr)
                if pi is None or pi >= len(Zs_panel):
                    continue
                z = Zs_panel[pi]
                # SurfaceL2 order 0: one DOF per element. el.nr is
                # the element index but the DOF mapping uses the
                # local SurfaceL2 dof for this element.
                dofs = fes_l2.GetDofNrs(el)
                for d in dofs:
                    if 0 <= d < len(re_vec):
                        re_vec[d] = z.real
                        im_vec[d] = z.imag

            # L2-project onto fes (the H1 space) by solving the
            # H1 mass equation: (u, v)_H1 = (gf_l2, v)_H1 for u in fes.
            # GridFunction.Set with VOL_or_BND=BND does this.
            gf_re = GridFunction(fes)
            gf_im = GridFunction(fes)
            with TaskManager():
                gf_re.Set(gf_l2_re, definedon=mesh_wp.Boundaries(".*"))
                gf_im.Set(gf_l2_im, definedon=mesh_wp.Boundaries(".*"))
            re = gf_re.vec.FV().NumPy().copy()
            im = gf_im.vec.FV().NumPy().copy()
            return (re + 1j * im).astype(complex), Zs_panel
        except Exception:
            # Fall through to vertex averaging if the L2 projection
            # fails for any reason (NGSolve API surface drift,
            # missing SurfaceL2 etc.). The user gets the order-1
            # answer instead of an error.
            pass

    # Path 1: order = 1 vertex averaging (default, fast)
    Zs_node = np.zeros(ndof_phi, dtype=complex)
    count = np.zeros(ndof_phi, dtype=int)
    for i, p in enumerate(panels):
        for vid in p['vert_ids']:
            if vid < ndof_phi:
                Zs_node[vid] += Zs_panel[i]
                count[vid] += 1

    # Nodes that no panel touched (e.g. orphan H1 dof) get the
    # mean Z_s as a safe default; in practice this should be empty
    # because every BND vertex belongs to >=1 panel.
    Zs_mean = complex(np.mean(Zs_panel)) if len(Zs_panel) > 0 else 0+0j
    for j in range(ndof_phi):
        if count[j] > 0:
            Zs_node[j] /= count[j]
        else:
            Zs_node[j] = Zs_mean
    return Zs_node, Zs_panel


def _run_coupled_bem(mesh_full, workpiece_label, source_label, sink_label,
                     fes_order, frequency, mat, half_thickness,
                     use_local_curvature=False, output_dir=None):
    """Iterative coil-workpiece BEM via bem_coupled_solver.CoupledBEMSolver.

    Builds clean coil + workpiece surface meshes from the same .vol
    (filtered by material name with orphan vertices pruned), then runs
    the iterative coupled solve with per-DOF back-reaction RHS.

    When ``use_local_curvature`` is True (Phase 5, 2026-04-12), the
    workpiece SIBC uses **per-node Z_s** built from per-panel local
    curvature instead of a single global Z_s value. This is the most
    accurate variant of the coupled BEM solver: each workpiece node
    gets a Z_s computed from its local sphere/cylinder radius
    extracted from the mesh.

    Sign of Delta_L (verified 2026-04-12):
      - mu_r = 1, copper, 1 kHz: -0.77 nH (Lenz screening)
      - mu_r = 1, copper, 1 MHz: -1.02 nH (PEC limit)
      - mu_r = 100, steel, 50 kHz: +0.34 nH (flux concentration starts)
      - mu_r = 1000, steel, 50 kHz: +1.30 nH (flux concentration dominates)

    Returns dict of coupled_* keys ready to merge into the main result.
    """
    setup_paths()
    from radia.bem_coupled_solver import CoupledBEMSolver
    from calc_heating_bem import _extract_surface_mesh_filtered
    import numpy as np

    sigma = mat.sigma
    mu_r = mat.mu_r
    material = mat.name

    progress("COUPLED", "extracting coil + workpiece surfaces")
    mesh_coil = _extract_surface_mesh_filtered(mesh_full, keep_label="coil")
    # The user-facing --workpiece argument is actually a BOUNDARY name
    # (the sideset, e.g. "sibc"), but _extract_surface_mesh_filtered
    # filters by VOLUME material. The .jou convention is to name the
    # workpiece block "workpiece" and the surface sideset "sibc", so
    # always use "workpiece" here. The new_to_old vertex map lets us
    # lift per-element/per-vertex workpiece quantities back onto
    # mesh_full so the GMSH visualization keeps the curved Tri6
    # rendering of the parent mesh.
    mesh_wp, wp_new_to_old = _extract_surface_mesh_filtered(
        mesh_full, keep_label="workpiece", return_vertex_map=True)

    rho = mat.rho
    omega = 2.0 * math.pi * frequency
    mu_eff = MU_0 * mu_r
    delta = mat.skin_depth(frequency)

    # Global Z_s (always computed for the legacy reporting fields)
    Zs_global = _dowell_Zs(half_thickness, rho, omega, mu_eff)

    progress("COUPLED", f"building coupled solver: coil {mesh_coil.nv}v "
                        f"/ wp {mesh_wp.nv}v, f={frequency} Hz")
    solver = CoupledBEMSolver(
        mesh_coil, mesh_wp,
        source_label=source_label, sink_label=sink_label,
        fes_order=fes_order)

    extra_keys = {}
    if use_local_curvature:
        # Per-panel R from the workpiece mesh, then per-node Z_s by
        # vertex averaging onto the SIBC solver's H1 nodes.
        wp_panels = _extract_panels_from_mesh(mesh_wp, workpiece_label)
        if not wp_panels:
            # Fall through to global Z_s; the user will see the same
            # answer as the non-local-curvature path. Errors out only
            # if the workpiece label is wrong.
            Z_s_arg = Zs_global
        else:
            # The extractor's adaptive panel-diameter floor handles
            # sliver triangles internally; the only global cap we
            # still apply is the upper bound = 1 m (= "essentially
            # flat"), so an unbounded default_R cannot leak out.
            R_panel = _compute_panel_local_radii(
                mesh_wp, wp_panels, default_R=1.0)
            Zs_node, Zs_panel = _build_per_node_Zs(
                mesh_wp, wp_panels, R_panel, rho, omega, mu_eff,
                solver.wp_solver.ndof,
                fes=solver.wp_solver.fes)
            Z_s_arg = Zs_node
            extra_keys.update({
                "coupled_use_local_curvature": True,
                "coupled_R_local_min_m": float(R_panel.min()),
                "coupled_R_local_max_m": float(R_panel.max()),
                "coupled_Zs_local_re_min": float(Zs_panel.real.min()),
                "coupled_Zs_local_re_max": float(Zs_panel.real.max()),
                "coupled_Zs_local_im_min": float(Zs_panel.imag.min()),
                "coupled_Zs_local_im_max": float(Zs_panel.imag.max()),
                "coupled_n_wp_panels": int(len(wp_panels)),
            })
    else:
        Z_s_arg = Zs_global
        extra_keys["coupled_use_local_curvature"] = False

    progress("COUPLED", "iterating coil-workpiece back-reaction...")
    sol = solver.solve(Z_s=Z_s_arg, omega=omega,
                       max_iter=10, tol=1e-3, relax=0.5)
    progress("COUPLED", f"converged in {int(sol['iterations'])} iter, "
                        f"L_total={sol['L_total']*1e9:.2f} nH")

    # ---- Lift wp_J onto mesh_full vertices for GMSH visualization ----
    # The CoupledBEMSolver returns wp_J as per-element complex vectors
    # on the *extracted* mesh_wp. To get a per-vertex |J|^2 on the
    # parent mesh_full we (1) average per-element |J| onto mesh_wp
    # vertices, (2) map mesh_wp vertices back to mesh_full vertices,
    # and (3) compute the heating density sigma * |J_complex|^2 / 2.
    # The result arrays span all of mesh_full.nv with zeros outside
    # the workpiece surface; the GMSH writer restricts the view to
    # the "workpiece" material so the zeros never show up.
    wp_c = sol.get('wp_c')
    wp_a = sol.get('wp_a')
    wp_J_re_arr = sol.get('wp_J_re')
    wp_J_im_arr = sol.get('wp_J_im')
    wp_J_full = None
    wp_q_full = None
    if (wp_c is not None and wp_a is not None
            and wp_J_re_arr is not None and len(wp_c) > 0):
        from ngsolve import BND
        wp_J_re_arr = np.asarray(wp_J_re_arr)
        wp_J_im_arr = np.asarray(wp_J_im_arr)
        # Per-vertex accumulators on the extracted workpiece mesh.
        Jre_v = np.zeros((mesh_wp.nv, 3))
        Jim_v = np.zeros((mesh_wp.nv, 3))
        cnt_v = np.zeros(mesh_wp.nv)
        for el_idx, el in enumerate(mesh_wp.Elements(BND)):
            if el_idx >= len(wp_J_re_arr):
                break
            for v in el.vertices:
                Jre_v[v.nr] += wp_J_re_arr[el_idx]
                Jim_v[v.nr] += wp_J_im_arr[el_idx]
                cnt_v[v.nr] += 1
        for k in range(3):
            Jre_v[:, k] = np.where(cnt_v > 0, Jre_v[:, k] / cnt_v, 0.0)
            Jim_v[:, k] = np.where(cnt_v > 0, Jim_v[:, k] / cnt_v, 0.0)
        # Lift to mesh_full vertices.
        wp_J_full = np.zeros((mesh_full.nv, 3))
        wp_q_full = np.zeros(mesh_full.nv)
        Jmag2_v = np.sum(Jre_v ** 2 + Jim_v ** 2, axis=1)
        # Heating density [W/m^2] for a SIBC conductor carrying surface
        # current J_s [A/m]:
        #
        #   Z_s  = (1 + j) / (sigma * delta)   [Ω/sq]  (SIBC)
        #   q''  = 0.5 * Re(Z_s) * |J_s|^2
        #        = 0.5 * |J_s|^2 / (sigma * delta)     [W/m^2]
        #
        # Previous formula `sigma * |J_s|^2 / 2` was dimensionally wrong
        # (A^2 * S / m^3, not W/m^2) and off by a factor of ~delta*sigma
        # ~ 1e-1, giving q ~10^8 W/m^2 where the physics expects ~1 W/m^2
        # for the IH BEM sample. Fixed 2026-04-15 after sugahara GUI
        # feedback ("発熱密度が一様で気持ち悪い" - values saturated the
        # log-scale legend because they were 8 orders of magnitude too
        # large). Jre/Jim here are the real / imaginary parts of the
        # surface current phasor (peak amplitudes), so the 0.5 factor
        # converts peak^2 to time-averaged power.
        q_v = 0.5 * Jmag2_v / (sigma * delta)
        for new_nr, old_nr in wp_new_to_old.items():
            if 0 <= new_nr < mesh_wp.nv and 0 <= old_nr < mesh_full.nv:
                wp_J_full[old_nr] = Jre_v[new_nr]
                wp_q_full[old_nr] = q_v[new_nr]

    result = {
        "coupled_L_air_H": float(sol['L_air']),
        "coupled_L_total_H": float(sol['L_total']),
        "coupled_dL_H": float(sol['Delta_L']),
        "coupled_P_total_W": float(sol['P_total']),
        "coupled_R_effective_Ohm": float(2.0 * sol['P_total']),
        "coupled_H_t_rms_Am": float(sol['H_t_rms']),
        "coupled_iterations": int(sol['iterations']),
        "coupled_n_J_coil": int(sol['n_J_coil']),
        "coupled_n_phi_wp": int(sol['n_phi_wp']),
        "coupled_Z_s_real_Ohm": float(Zs_global.real),
        "coupled_Z_s_imag_Ohm": float(Zs_global.imag),
        "coupled_delta_skin_m": float(delta),
        "coupled_mu_r": float(mu_r),
        "coupled_sigma_Sm": float(sigma),
        "coupled_half_thickness_m": float(half_thickness),
        "coupled_material": material,
    }
    result.update(extra_keys)
    # Pass per-vertex workpiece quantities back to the caller as
    # private keys (underscore prefix); the caller pops them before
    # returning the JSON-serializable result dict to the panel.
    if wp_J_full is not None:
        result["_wp_J_full"] = wp_J_full
        result["_wp_q_full"] = wp_q_full
        # Save workpiece mesh + per-vertex J / q as .vol + .sol pair
        # (NGSolve native format). The unified vol2msh converter at the
        # end of extract_inductance_vol reads these and produces the
        # combined GMSH .msh — no inline GMSH writes here.
        try:
            from ngsolve import H1, GridFunction
            wp_dir = output_dir if output_dir else "."
            wp_vol_path = os.path.join(wp_dir, "workpiece.vol").replace("\\", "/")
            wp_J_sol = os.path.join(wp_dir, "wp_J.sol").replace("\\", "/")
            wp_q_sol = os.path.join(wp_dir, "wp_q.sol").replace("\\", "/")

            # Save extracted workpiece mesh.
            mesh_wp.ngmesh.Save(wp_vol_path)

            # wp_J as H1 vector-3 GridFunction on extracted wp mesh.
            fes_J = H1(mesh_wp, order=1, dim=3)
            gf_J_wp = GridFunction(fes_J)
            # H1 dim=3 layout: interleaved (x0,y0,z0, x1,y1,z1, ...) with
            # vertex-first DOF ordering for order=1.
            for v_idx in range(mesh_wp.nv):
                for k in range(3):
                    gf_J_wp.vec.FV()[v_idx * 3 + k] = float(Jre_v[v_idx, k])
            gf_J_wp.Save(wp_J_sol)

            # wp_q as H1 scalar GridFunction on extracted wp mesh.
            fes_q = H1(mesh_wp, order=1)
            gf_q_wp = GridFunction(fes_q)
            for v_idx in range(mesh_wp.nv):
                gf_q_wp.vec.FV()[v_idx] = float(q_v[v_idx])
            gf_q_wp.Save(wp_q_sol)

            result["_wp_vol_path"] = wp_vol_path
            result["_wp_J_sol_path"] = wp_J_sol
            result["_wp_q_sol_path"] = wp_q_sol
            from calc_common import progress as _progress
            _progress("SAVE",
                      f"workpiece.vol + wp_J.sol + wp_q.sol -> {wp_dir}")
        except Exception as _e:
            sys.stderr.write(f"WARN: workpiece .vol/.sol save failed: {_e}\n")
            sys.stderr.flush()
    return result


def _extract_panels_from_mesh(mesh, workpiece_label):
    """Extract workpiece surface panels from mesh boundary elements.

    Finds boundary elements whose label starts with *workpiece_label*
    (prefix match, case-insensitive) and computes per-element centroid,
    outward normal, area, and the **vertex indices** (used downstream
    by ``_compute_panel_local_radii`` for per-panel local sphere fit).

    Returns list of dicts with keys ``center``, ``normal``, ``area``,
    ``vert_ids``.
    """
    import numpy as np
    from ngsolve import BND

    bnd_labels = list(set(mesh.GetBoundaries()))
    wp_labels = _resolve_bnd_labels(workpiece_label, bnd_labels)
    if not wp_labels:
        return []

    wp_label_set = set(wp_labels)
    panels = []
    for el in mesh.Elements(BND):
        if el.mat not in wp_label_set:
            continue
        vert_ids = [v.nr for v in el.vertices]
        verts = [mesh.vertices[v.nr].point for v in el.vertices]
        pts = np.array([(v[0], v[1], v[2]) for v in verts])
        center = pts.mean(axis=0)
        # Triangle normal via cross product
        if len(pts) >= 3:
            e1 = pts[1] - pts[0]
            e2 = pts[2] - pts[0]
            cross = np.cross(e1, e2)
            area = 0.5 * np.linalg.norm(cross)
            normal = cross / (2.0 * area) if area > 1e-30 else np.array([0, 0, 1.0])
        else:
            area = 0.0
            normal = np.array([0, 0, 1.0])
        if area > 1e-30:
            panels.append({
                'center': center,
                'normal': normal,
                'area': area,
                'vert_ids': vert_ids,
            })
    return panels


def _compute_panel_local_radii(mesh, panels, default_R=1e30):
    """Compute a per-panel local radius via discrete mean curvature.

    Uses the standard discrete-differential-geometry formula based on
    the **angle between adjacent panel normals**. For each panel ``i``
    we look at its **edge-adjacent** panels ``j`` (sharing 2 vertices,
    i.e. an edge) and compute::

        angle_ij  = arccos(n_i . n_j)
        dist_ij   = |c_j - c_i|
        kappa_ij  ~ angle_ij / dist_ij           (local curvature)
        R_ij      = 1 / kappa_ij = dist_ij / angle_ij

    The panel's local radius is the **median** of the per-edge
    estimates (median is robust to outliers from sliver triangles).

    This formula is the limiting case of the integrated geodesic
    curvature in 1D and gives the right answer for the three
    canonical local shapes:

      - **Sphere** R: all neighbors are at the same arc-length and
        the normal rotates by the same angle -> R_ij = R exactly.
      - **Cylinder** R: edge-adjacent panels around the
        circumferential direction give R; edge-adjacent panels in
        the axial direction give parallel normals (angle ~ 0,
        radius -> infinity). The median picks up the smaller R
        from circumferential neighbors -> R_local = R.
      - **Flat plate**: parallel normals everywhere, angle ~ 0,
        radius -> infinity -> clipped to ``default_R``.

    Returns ``np.ndarray[len(panels)]`` of local radii [m]. Panels
    with no edge neighbors fall back to ``default_R``.
    """
    import numpy as np

    if not panels:
        return np.array([])

    # vertex -> list of panel indices that touch it
    vert_to_panels = {}
    for pi, p in enumerate(panels):
        for vid in p['vert_ids']:
            vert_to_panels.setdefault(vid, []).append(pi)

    # Pre-compute per-panel diameter from the panel's vertex coordinates.
    # The local radius cannot meaningfully be smaller than the panel
    # diameter (the chord-vs-arc approximation breaks down) — this is
    # the **adaptive sliver-suppression bound**, more meaningful than
    # the previous global ``0.5 * half_thickness`` clamp.
    panel_diam = np.zeros(len(panels))
    for pi, p in enumerate(panels):
        coords = np.array(
            [mesh.vertices[v].point for v in p['vert_ids']], dtype=float)
        # Max pairwise distance among the (3) triangle vertices
        d = 0.0
        for i in range(len(coords)):
            for j in range(i + 1, len(coords)):
                dij = float(np.linalg.norm(coords[i] - coords[j]))
                if dij > d:
                    d = dij
        panel_diam[pi] = d

    R_local = np.full(len(panels), float(default_R))
    MIN_ANGLE = 1e-4  # below this, treat as flat (numerical noise)

    for pi, p in enumerate(panels):
        n_i = np.asarray(p['normal'], dtype=float)
        n_i /= max(float(np.linalg.norm(n_i)), 1e-30)
        c_i = np.asarray(p['center'], dtype=float)

        # Find directly edge-adjacent panels (share exactly 2 vertices).
        own_vids = set(p['vert_ids'])
        candidate_neighbors = set()
        for vid in own_vids:
            for qi in vert_to_panels.get(vid, []):
                if qi != pi:
                    candidate_neighbors.add(qi)
        edge_neighbors = [
            qi for qi in candidate_neighbors
            if len(own_vids & set(panels[qi]['vert_ids'])) == 2]

        if not edge_neighbors:
            continue

        # Adaptive sliver lower bound: a panel of diameter h cannot
        # resolve a curvature radius below h/2 (the chord-vs-arc
        # approximation breaks). Drop any R_est below this floor.
        # Above the floor, the discrete formula is still meaningful.
        R_floor = 0.5 * float(panel_diam[pi])

        radii = []
        for qi in edge_neighbors:
            n_q = np.asarray(panels[qi]['normal'], dtype=float)
            n_q /= max(float(np.linalg.norm(n_q)), 1e-30)
            c_q = np.asarray(panels[qi]['center'], dtype=float)

            cos_angle = float(np.clip(np.dot(n_i, n_q), -1.0, 1.0))
            angle = float(np.arccos(cos_angle))
            if angle < MIN_ANGLE:
                continue  # parallel normals -> flat -> infinite R

            dist = float(np.linalg.norm(c_q - c_i))
            if dist < 1e-30:
                continue

            R_est = dist / angle
            # Filter out spurious tiny R from sliver triangles, fold-
            # over panels (e.g. gap end-caps), or cross-mesh edge
            # neighbors. Anything below the panel-diameter floor is
            # geometrically suspicious.
            if R_est < R_floor:
                continue
            if R_est > default_R or R_est <= 0:
                continue
            radii.append(R_est)

        if not radii:
            continue
        # We want the **maximum local curvature direction** = the
        # smallest R, because the SIBC cylinder cell problem is solved
        # in the most-curved 1D radial coordinate. For:
        #   - sphere R: every neighbor gives R, percentile-10 = R
        #   - cylinder R: circumferential neighbors give R,
        #     axial neighbors give large R (or are filtered out by
        #     MIN_ANGLE), percentile-10 picks up the small R = R
        #   - flat: all R large, percentile-10 large -> default_R clip
        # We take percentile 10 (not the strict min) for robustness
        # against sliver triangles that produce spurious tiny R.
        radii_arr = np.asarray(radii)
        R_local[pi] = float(np.percentile(radii_arr, 10))

    return R_local


def _compute_wp_impedance_from_panels(mesh_coil, gf_J, panels,
                                       impedance_model="esim", frequency=50000,
                                       mat=None, half_thickness=0.005,
                                       esim_geometry="cylinder",
                                       mesh_wp=None,
                                       use_local_curvature=False):
    """Compute workpiece surface impedance from mesh panels.

    Pipeline:
      1. Biot-Savart from coil J -> H at workpiece panels
      2. ESIM cell problem or Dowell formula -> Z_s, P, Q per panel
         If ``use_local_curvature`` is True, each panel's R is taken
         from the mesh-driven local curvature instead of the global
         ``half_thickness``.
      3. Integrate over surface

    Args:
        mesh_wp: workpiece NGSolve mesh (required when
            use_local_curvature=True; used to compute per-panel R)
        use_local_curvature: when True, the SIBC R is per-panel,
            extracted from the mesh by ``_compute_panel_local_radii``.
            Sphere / ellipsoid / curved-shell workpieces benefit
            because the cell problem now matches the local geometry
            instead of using a single global radius.

    Returns dict with wp_* keys to merge into main result.
    """
    import math
    import numpy as np
    from ngsolve import Integrate, CF, BND

    if mat is None:
        mat = EMMaterial.from_name("steel")
    sigma = mat.sigma
    mu_r = mat.mu_r
    material = mat.name

    if not panels:
        return {"wp_error": "No workpiece surfaces found"}

    n_panels = len(panels)
    rho = mat.rho
    omega = 2 * math.pi * frequency

    # --- Biot-Savart H at workpiece panels ---
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

    # --- Surface impedance per panel ---
    P_total = 0.0
    Q_total = 0.0
    Z_sum = 0.0 + 0.0j
    delta_min = float('inf')
    delta_max = 0.0

    # Per-panel local radii from the mesh (for the curvature-aware
    # SIBC variant). All radii default to ``half_thickness`` so when
    # use_local_curvature=False the legacy global-R behaviour is
    # exactly preserved.
    R_per_panel = np.full(n_panels, float(half_thickness))
    R_local_min = float(half_thickness)
    R_local_max = float(half_thickness)
    if use_local_curvature and mesh_wp is not None:
        # The extractor's adaptive panel-diameter floor handles
        # sliver suppression internally. Only the upper bound (=1 m,
        # essentially flat) needs to be enforced here so unbounded
        # default_R values cannot leak through.
        R_per_panel = _compute_panel_local_radii(
            mesh_wp, panels, default_R=1.0)
        R_local_min = float(R_per_panel.min())
        R_local_max = float(R_per_panel.max())

    # "sibc" is the user-facing name for the linear analytical Dowell SIBC.
    if impedance_model == "sibc":
        impedance_model = "dowell"
    if impedance_model == "esim":
        solver = mat.create_esim_solver(frequency, half_thickness,
                                        geometry=esim_geometry)

        for i in range(n_panels):
            H0 = max(float(H_t_mag[i]), 1e-3)
            R_i = float(R_per_panel[i])
            if use_local_curvature:
                sol = solver.solve(H0, R_local=R_i)
            else:
                sol = solver.solve(H0)
            P_total += sol['P_prime'] * panels[i]['area']
            Q_total += sol['Q_prime'] * panels[i]['area']
            Z_sum += sol['Z'] * H_t_mag[i]**2 * panels[i]['area']
            delta_min = min(delta_min, sol['delta'])
            delta_max = max(delta_max, sol['delta'])

    elif impedance_model == "dowell":
        mu_eff = MU_0 * mu_r
        delta = math.sqrt(2 * rho / (omega * mu_eff)) if omega > 0 else 1e10
        delta_min = delta_max = delta

        # Per-panel Z_s using per-panel R (Dowell tanh)
        for i in range(n_panels):
            R_i = float(R_per_panel[i])
            xi = R_i / delta
            gamma_a = complex(1, 1) * xi
            try:
                Z_s_i = (rho / R_i) * gamma_a * np.tanh(gamma_a)
            except OverflowError:
                Z_s_i = (rho / R_i) * gamma_a
            H0 = float(H_t_mag[i])
            area_i = panels[i]['area']
            P_total += 0.5 * Z_s_i.real * H0 ** 2 * area_i
            Q_total += 0.5 * Z_s_i.imag * H0 ** 2 * area_i
            Z_sum += Z_s_i * H0 ** 2 * area_i

    else:
        return {"wp_error": f"Unknown impedance model: {impedance_model}"}

    R_eff = float(Z_sum.real)
    X_eff = float(Z_sum.imag)
    wp_area = sum(p['area'] for p in panels)

    return {
        "wp_model": impedance_model,
        "wp_material": material,
        "wp_frequency": frequency,
        "wp_sigma": sigma,
        "wp_half_thickness": half_thickness,
        "wp_use_local_curvature": bool(use_local_curvature),
        "wp_R_local_min": R_local_min,
        "wp_R_local_max": R_local_max,
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

    Loads Solve results (.vol + .sol or .npy). No Cubit needed.

    Args:
        mesh_vol_path: Path to surface_mesh.vol from solve step
        fes_order: HDivSurface basis order
        msh_output: Optional .msh output path
        j_npy: Path to J_coeffs.sol or J_coeffs.npy from solve step
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

    # Load J coefficients (.sol or .npy)
    fes_J = HDivSurface(mesh, order=fes_order)
    gf_J = GridFunction(fes_J)
    if j_npy.endswith('.sol'):
        gf_J.Load(j_npy)
    else:
        J = np.load(j_npy)
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

    # Write a SINGLE v4.1 .msh with:
    #   - volume tets (air box) carrying B NodeData
    #   - Tri6 curved surface (coil) carrying J NodeData
    #   - coil wireframe (1D edges) for visual reference
    #
    # Coil surface uses Tri6 (curved) elements from mesh.Curve(2)
    # for the smooth "kirei" rendering. .msh.opt auto-loaded by GMSH
    # sets arrow sizes and hides air mesh edges.
    # 2026-04-14 Phase A.3: write to user's msh_output (single file with
    # B + coil J + workpiece overlay) instead of a hardcoded
    # "inductance.msh" that left the panel-output and post-output files
    # diverged.
    if msh_output:
        gmsh_file = msh_output.replace("\\", "/")
    else:
        gmsh_file = os.path.join(
            base_dir, "inductance.msh").replace("\\", "/")

    # B per vertex on the volume mesh
    fes_B = HDiv(mesh_vol, order=1)
    gf_B = GridFunction(fes_B)
    gf_B.Set(B_cf)

    # Save B (.vol + .sol pair for the air volume mesh)
    air_vol = os.path.join(base_dir, "air_mesh.vol").replace("\\", "/")
    b_sol = os.path.join(base_dir, "B_field.sol").replace("\\", "/")
    mesh_vol.ngmesh.Save(air_vol)
    gf_B.Save(b_sol)
    progress("SAVE", f"air_mesh.vol + B_field.sol")

    vol_nodes_xyz = [(v.point[0], v.point[1], v.point[2])
                     for v in mesh_vol.vertices]
    B_nodes = np.zeros((len(vol_nodes_xyz), 3))
    for vi, v in enumerate(mesh_vol.vertices):
        bval = gf_B(mesh_vol(*v.point))
        B_nodes[vi] = [float(bval[k]) for k in range(3)]

    # Volume elements (tets)
    vol_elems = [[v.nr for v in el.vertices]
                 for el in mesh_vol.Elements()]

    # Phase A.2: load workpiece .vol + wp_q.sol if present and pass to
    # the combined .msh writer. This is how the GMSH viewer shows
    # workpiece + heating density alongside coil J + air B in ONE file
    # (single source of truth = the .vol+.sol pairs the solver wrote).
    wp_nodes = wp_elems = q_nodes = wp_J_nodes = None
    wp_vol_path = os.path.join(base_dir, "workpiece.vol").replace("\\", "/")
    wp_q_sol_path = os.path.join(base_dir, "wp_q.sol").replace("\\", "/")
    wp_J_sol_path = os.path.join(base_dir, "wp_J.sol").replace("\\", "/")
    if os.path.isfile(wp_vol_path) and os.path.isfile(wp_q_sol_path):
        try:
            from netgen.meshing import Mesh as NgMesh
            from ngsolve import H1, GridFunction, Mesh as NgsMesh
            wp_ngmesh = NgMesh()
            wp_ngmesh.Load(wp_vol_path)
            wp_mesh_obj = NgsMesh(wp_ngmesh)
            wp_nodes = [(v.point[0], v.point[1], v.point[2])
                         for v in wp_mesh_obj.vertices]
            wp_elems = [[v.nr for v in el.vertices]
                         for el in wp_mesh_obj.Elements(BND)]
            fes_q = H1(wp_mesh_obj, order=1)
            gf_q_load = GridFunction(fes_q)
            gf_q_load.Load(wp_q_sol_path)
            q_nodes = np.array(
                [float(gf_q_load.vec.FV()[i]) for i in range(wp_mesh_obj.nv)])

            # Eddy-current vector on wp surface (optional companion)
            if os.path.isfile(wp_J_sol_path):
                fes_J_wp = H1(wp_mesh_obj, order=1, dim=3)
                gf_J_wp = GridFunction(fes_J_wp)
                gf_J_wp.Load(wp_J_sol_path)
                wp_J_nodes = np.zeros((wp_mesh_obj.nv, 3))
                for v_idx in range(wp_mesh_obj.nv):
                    for k in range(3):
                        wp_J_nodes[v_idx, k] = float(
                            gf_J_wp.vec.FV()[v_idx * 3 + k])

            progress("GMSH", f"loaded workpiece overlay: "
                              f"{len(wp_nodes)} nodes, {len(wp_elems)} tris"
                              + (", J_workpiece" if wp_J_nodes is not None else ""))
        except Exception as _e:
            sys.stderr.write(f"WARN: could not load workpiece overlay: {_e}\n")
            sys.stderr.flush()
            wp_nodes = wp_elems = q_nodes = wp_J_nodes = None

    _write_combined_msh_v41(
        gmsh_file, mesh, J_nodes, vol_nodes_xyz, vol_elems, B_nodes,
        wp_nodes=wp_nodes, wp_elems=wp_elems, q_nodes=q_nodes,
        wp_J_nodes=wp_J_nodes)
    progress("GMSH", f"B + J + wireframe + workpiece (v4.1): {gmsh_file}")

    # .msh.opt is written by _write_combined_msh_v41 (replaces .geo)

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
        "gmsh_file": gmsh_file,
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




def run_bem_pipeline(vol_file, source_label="source", sink_label="sink",
                     coil_label="coil", wp_label="workpiece",
                     fes_order=0, frequency=10000, mat=None,
                     impedance_model="linear",
                     half_thickness=0.010, coil_radius=0.030,
                     coil_current=100.0, field_air=True,
                     lx=0.06, ly=0.06, lz=0.03, maxh_vol=0.015):
    """Full BEM pipeline: .vol -> .sol files + .msh visualization.

    Input:
        {name}.vol       Cubit-exported mesh (export netgen order 2)

    Output (all in same directory as .vol):
        {name}_coil.vol  Coil surface mesh
        {name}_J.sol     Coil current density J (HDivSurface)
        {name}_air.vol   Air volume mesh (Biot-Savart box)
        {name}_B.sol     Magnetic flux density B (HDiv)
        {name}_wp.vol    Workpiece surface mesh
        {name}_q.sol     Heating density q = sigma*|H_t|^2/2
        {name}.msh       GMSH visualization (B + J + q + wireframe)
        {name}.msh.opt   GMSH display options

    Args:
        vol_file: Path to .vol file from Cubit export.
        source_label: Boundary label for current source face.
        sink_label: Boundary label for current sink face.
        coil_label: Material name of the coil block.
        wp_label: Material name of the workpiece block.
        fes_order: HDivSurface basis order (0=RWG).
        frequency: Operating frequency [Hz].
        mat: EMMaterial instance (workpiece properties).
        impedance_model: "linear" or "esim".
        half_thickness: Workpiece thickness for ESIM [m].
        coil_radius: Coil loop radius [m] (for workpiece phi_inc).
        coil_current: Coil current [A].
        field_air: Compute B field in air box.
        lx, ly, lz: Air box half-sizes [m].
        maxh_vol: Air box element size [m].

    Returns:
        dict with file paths and solver results.
    """
    base = os.path.splitext(vol_file)[0]
    base_dir = os.path.dirname(os.path.abspath(vol_file))

    # === 1. Coil BEM: J ===
    progress("PIPELINE", "Step 1: Coil BEM solve")
    result = extract_inductance_vol(
        vol_file=vol_file,
        source_label=source_label, sink_label=sink_label,
        fes_order=fes_order, coil_label=coil_label,
        field_air=field_air,
        msh_output=base + '.msh',
    )
    if 'error' in result:
        return result

    # Rename to {name}_xxx convention
    _rename_if_exists(base_dir, 'surface_mesh.vol', base + '_coil.vol')
    _rename_if_exists(base_dir, 'J_coeffs.sol', base + '_J.sol')
    _rename_if_exists(base_dir, 'J_coeffs.npy', base + '_J.npy')
    _rename_if_exists(base_dir, 'air_mesh.vol', base + '_air.vol')
    _rename_if_exists(base_dir, 'B_field.sol', base + '_B.sol')
    progress("PIPELINE", f"Coil: {base}_coil.vol + _J.sol")

    # === 2. Workpiece BEM: q ===
    wp_result = {}
    wp_vol = base + '_wp.vol'
    q_sol = base + '_q.sol'
    try:
        sys.path.insert(0, os.path.dirname(__file__))
        from calc_heating_bem import compute_heating_bem
        progress("PIPELINE", "Step 2: Workpiece BEM heating")
        wp_result = compute_heating_bem(
            vol_file=vol_file,
            coil_radius=coil_radius, coil_current=coil_current,
            frequency=frequency, mat=mat,
            impedance_model=impedance_model,
            half_thickness=half_thickness,
            wp_label=wp_label,
        )
        _rename_if_exists(base_dir, 'workpiece.vol', wp_vol)
        _rename_if_exists(base_dir, 'q_heating.sol', q_sol)
        progress("PIPELINE", f"WP: P={wp_result.get('P_total_W', 0):.4e} W")
    except Exception as e:
        progress("PIPELINE", f"Workpiece skipped: {e}")
        wp_vol = None
        q_sol = None

    # === 3. Convert .sol -> .msh ===
    progress("PIPELINE", "Step 3: convert_sol_to_msh")
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
    from gmsh_post_export import convert_sol_to_msh
    coil_vol = base + '_coil.vol'
    j_sol = base + '_J.sol'
    air_vol = base + '_air.vol' if field_air else None
    b_sol = base + '_B.sol' if field_air else None

    convert_sol_to_msh(
        output_msh=base + '.msh',
        surface_vol=coil_vol, j_sol=j_sol,
        air_vol=air_vol, b_sol=b_sol,
        wp_vol=wp_vol if wp_vol and os.path.exists(wp_vol) else None,
        q_sol=q_sol if q_sol and os.path.exists(q_sol) else None,
        fes_order=fes_order,
    )
    progress("PIPELINE", f"Done: {base}.msh")

    return {
        "vol": vol_file,
        "coil_vol": coil_vol, "j_sol": j_sol,
        "air_vol": air_vol, "b_sol": b_sol,
        "wp_vol": wp_vol, "q_sol": q_sol,
        "msh": base + '.msh',
        "L_nH": result.get('L_nH'),
        "P_W": wp_result.get('P_total_W'),
        "H_t_rms": wp_result.get('H_t_rms_Am'),
    }


def _rename_if_exists(base_dir, old_name, new_path):
    """Rename a file from base_dir/old_name to new_path."""
    old = os.path.join(base_dir, old_name)
    if os.path.exists(old):
        if os.path.exists(new_path):
            os.remove(new_path)
        os.rename(old, new_path)


def main():
    parser = argparse.ArgumentParser(
        description="BEM inductance (source/sink saddle point EFIE)")
    parser.add_argument("--mode", default="solve", choices=["solve", "post"],
                        help="solve: BEM solve only; post: B-field + GMSH export")
    # Solve mode args
    parser.add_argument("--vol", required=True, help="Netgen .vol file")
    parser.add_argument("--fes-order", type=int, default=0, help="HDivSurface order (0=RWG)")
    parser.add_argument("--source", default="source",
                        help="Source boundary label (default 'source' by .jou convention)")
    parser.add_argument("--sink", default="sink",
                        help="Sink boundary label (default 'sink' by .jou convention)")
    parser.add_argument("--msh-output", default="", help="GMSH .msh output path")
    parser.add_argument("--output", default="", help="Output JSON file (optional)")
    # Coil / Workpiece args
    parser.add_argument("--coil", default="coil",
                        help="Coil material name (default 'coil' by .jou convention; filters BEM to coil surfaces)")
    parser.add_argument("--coil-sigma", type=float, default=0,
                        help="Coil conductivity [S/m] (stored in results)")
    parser.add_argument("--workpiece", default="", help="Workpiece boundary label (prefix match)")
    parser.add_argument("--impedance-model", default="esim",
                        choices=["esim", "sibc", "dowell"],
                        help="esim=ESIM cell problem, sibc=linear surface "
                             "impedance (per-panel curvature). 'dowell' is "
                             "kept as a deprecated alias for sibc.")
    parser.add_argument("--frequency", type=float, default=50000, help="Frequency [Hz]")
    parser.add_argument("--current", type=float, default=1.0,
                        help="Coil terminal current [A] (1 A reference; "
                             "BEM solves a unit-current problem and the "
                             "user can rescale post hoc - accepted here "
                             "for symmetry with the FEM driver)")
    parser.add_argument("--half-thickness", type=float, default=0.005,
                        help="Slab half-thickness [m]")
    add_material_args(parser, include_custom=False)
    parser.add_argument("--esim-geometry", default="local_curvature",
                        choices=["local_curvature", "none"],
                        help="ESIM curvature: local_curvature=Bessel (R=half_thickness), "
                             "none=flat slab (cosh/sinh)")
    parser.add_argument("--solver", default="lu",
                        choices=["lu", "minres", "gmres"],
                        help="BEM saddle-point solver: lu=dense direct (default, "
                             "best for ndof<5000), minres=Krylov symmetric "
                             "indefinite (large problems), gmres=general Krylov")
    parser.add_argument("--field-air", action="store_true",
                        help="After BEM solve, also compute B-field in the "
                             "air domain via Biot-Savart and write a "
                             "GMSH .msh + .geo for visualization.")
    parser.add_argument("--use-local-curvature", action="store_true",
                        help="Use mesh-driven per-panel local curvature "
                             "for the SIBC. Each workpiece panel gets its "
                             "own local radius from a discrete normal-angle "
                             "fit. Replaces the global half_thickness "
                             "assumption - required for correct sphere or "
                             "doubly-curved workpiece SIBC.")
    # Post mode args
    parser.add_argument("--j-npy", default="", help="J_coeffs.npy path (post mode)")
    parser.add_argument("--mesh-vol", default="", help="surface_mesh.vol path (post mode)")
    parser.add_argument("--lx", type=float, default=0.07, help="Box half-size X [m]")
    parser.add_argument("--ly", type=float, default=0.07, help="Box half-size Y [m]")
    parser.add_argument("--lz", type=float, default=0.07, help="Box half-size Z [m]")
    parser.add_argument("--maxh-vol", type=float, default=0.01, help="Volume mesh size [m]")

    def run(args):
        _geom_map = {"local_curvature": "cylinder", "none": "slab"}
        esim_geom = _geom_map.get(args.esim_geometry, "cylinder")

        if args.mode == "solve":
            # SIBC and the legacy 'dowell' label both map to the
            # linear surface impedance path; the rest of this module
            # still uses "dowell" as the internal key.
            imp_model = args.impedance_model
            if imp_model == "sibc":
                imp_model = "dowell"
            return extract_inductance_vol(
                args.vol, args.source, args.sink,
                args.fes_order, args.msh_output,
                coil_label=args.coil,
                workpiece=args.workpiece,
                impedance_model=imp_model,
                frequency=args.frequency,
                mat=EMMaterial.from_args(args),
                half_thickness=args.half_thickness,
                esim_geometry=esim_geom,
                coil_sigma=args.coil_sigma,
                solver=args.solver,
                field_air=args.field_air,
                use_local_curvature=args.use_local_curvature)
        else:
            return post_process(args.mesh_vol, args.fes_order,
                                args.msh_output, args.j_npy,
                                args.lx, args.ly, args.lz,
                                args.maxh_vol)

    calc_main(run, parser)


if __name__ == "__main__":
    main()
