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
if _this_dir not in sys.path:
    sys.path.insert(0, _this_dir)

from calc_common import (setup_paths, MU_0, progress, calc_main,
                          create_esim_solver)


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
    mesh.Curve(order)
    return mesh


def compute_heating_bem(vol_file, coil_radius=0.030, coil_current=1.0,
                        gap_deg=5, frequency=7000, sigma=2e6,
                        material="steel", h1_order=1, wp_label="wp_surface",
                        half_thickness=0.010, esim_geometry="cylinder",
                        bh_file=None, impedance_model="esim", mu_r=None,
                        max_iter=15, tol=1e-3,
                        msh_output=""):
    """BEM-SIBC workpiece eddy current from .vol file.

    Args:
        vol_file: Path to Netgen .vol (volume or surface mesh)
        coil_radius: Coil loop radius [m]
        coil_current: Total coil current [A]
        gap_deg: Coil gap angle [degrees] (0 = closed loop)
        frequency: Operating frequency [Hz]
        sigma: Workpiece conductivity [S/m]
        material: "steel", "copper", or "aluminum"
        h1_order: H1 polynomial order on surface
        wp_label: Boundary label for workpiece surface
        half_thickness: Workpiece radius for ESIM [m]
        esim_geometry: "cylinder" or "planar"
        bh_file: Optional custom BH curve file
        max_iter: Karl iteration maximum
        tol: Karl convergence tolerance
        msh_output: Optional GMSH output path

    Returns:
        dict with P_total, H_t_rms, Z_s, etc.
    """
    import time as _time
    from ngsolve import Mesh, Integrate, CF, BND

    setup_paths()
    from bem_sibc_solver import ScalarBIESIBCSolver, compute_phi_inc_from_loop

    omega = 2 * np.pi * frequency
    t_total_start = _time.perf_counter()

    # === 1. Load mesh and extract surface ===
    progress("MESH", "Loading .vol file...")
    mesh = Mesh(vol_file)

    if mesh.ne > 0:
        mesh = _extract_surface_mesh(mesh, order=2)

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

        esim_solver = create_esim_solver(
            material=material, frequency=frequency, sigma=sigma,
            half_thickness=half_thickness, geometry=esim_geometry,
            bh_file=bh_file)

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

    # === 5. GMSH output (optional) ===
    gmsh_file = ""
    if msh_output:
        try:
            from gmsh_post_export import GmshPostExport
            post = GmshPostExport(mesh, boundary=(mesh.ne > 0))

            phi_vec = result['phi_vec']
            from ngsolve import GridFunction, grad, Norm, sqrt

            gf_phi_re = GridFunction(solver.fes)
            gf_phi_re.vec.FV().NumPy()[:] = phi_vec.real
            gf_phi_im = GridFunction(solver.fes)
            gf_phi_im.vec.FV().NumPy()[:] = phi_vec.imag

            # |H_t|^2 = |grad phi_re|^2 + |grad phi_im|^2
            H_t_sq = Norm(grad(gf_phi_re))**2 + Norm(grad(gf_phi_im))**2
            H_t_cf = sqrt(H_t_sq)

            # P_loss [W/m^2] = 0.5 * Re(Z_s) * |H_t|^2
            Re_Zs = Z_s.real
            P_cf = 0.5 * Re_Zs * H_t_sq

            # Q_reactive [var/m^2] = 0.5 * Im(Z_s) * |H_t|^2
            Im_Zs = Z_s.imag
            Q_cf = 0.5 * Im_Zs * H_t_sq

            post.add_field("|H_t| [A/m]", H_t_cf, ncomp=1)
            post.add_field("P_loss [W/m^2]", P_cf, ncomp=1)
            post.add_field("Q_reactive [var/m^2]", Q_cf, ncomp=1)

            post.write(msh_output)

            # Overwrite .geo companion with SIBC parameters
            geo_path = os.path.splitext(msh_output)[0] + ".geo"
            with open(geo_path, "w") as gf:
                gf.write("// BEM-SIBC workpiece result\n")
                gf.write(f"// Material: {material}\n")
                gf.write(f"// WP sigma = {sigma:.4e} S/m\n")
                gf.write(f"// WP mu_r = {mu_r_eff:.1f}\n")
                gf.write(f"// Coil sigma = 5.8e+07 S/m (copper)\n")
                gf.write(f"// frequency = {frequency:.0f} Hz\n")
                gf.write(f"// |Z_s| = {abs(Z_s):.4e} Ohm "
                         f"(phase = {np.degrees(np.angle(Z_s)):.1f} deg)\n")
                gf.write(f"// Re(Z_s) = {Re_Zs:.4e}, "
                         f"Im(Z_s) = {Im_Zs:.4e}\n")
                gf.write(f"// skin_depth = {delta*1e3:.3f} mm\n")
                gf.write(f"// half_thickness = {half_thickness*1e3:.1f} mm\n")
                gf.write(f"// coil_radius = {coil_radius*1e3:.1f} mm\n")
                gf.write(f"// coil_current = {coil_current:.2f} A\n")
                gf.write(f"// P_total = {P_total:.4e} W\n")
                gf.write(f"// Q_total = {Q_total:.4e} var\n")
                gf.write(f"// H_t_rms = {H_t_rms:.4f} A/m\n")
                gf.write(f"// Karl iterations = {n_converged}\n")
                gf.write(f'Merge "{os.path.basename(msh_output)}";\n')
                gf.write("Mesh.NumSubEdges = 4;\n")
                gf.write("Mesh.VolumeEdges = 0;\n")
                gf.write("Mesh.SurfaceEdges = 0;\n")
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
    parser.add_argument("--sigma", type=float, default=2e6,
                        help="Workpiece conductivity [S/m]")
    parser.add_argument("--material", default="steel",
                        choices=["steel", "copper", "aluminum"])
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
    parser.add_argument("--mu-r", type=float, default=None,
                        help="Relative permeability (linear mode)")
    parser.add_argument("--bh-file", default="",
                        help="Custom BH curve file")
    parser.add_argument("--max-iter", type=int, default=15,
                        help="Karl iteration max")
    parser.add_argument("--tol", type=float, default=1e-3,
                        help="Karl convergence tolerance")
    parser.add_argument("--msh-output", default="",
                        help="GMSH output path")
    parser.add_argument("--output", default="",
                        help="JSON output file")

    def run(args):
        return compute_heating_bem(
            vol_file=args.vol,
            coil_radius=args.coil_radius,
            coil_current=args.coil_current,
            gap_deg=args.gap_deg,
            frequency=args.frequency,
            sigma=args.sigma,
            impedance_model=args.impedance_model,
            mu_r=args.mu_r,
            material=args.material,
            h1_order=args.h1_order,
            wp_label=args.wp_label,
            half_thickness=args.half_thickness,
            esim_geometry=args.esim_geometry,
            bh_file=args.bh_file or None,
            max_iter=args.max_iter,
            tol=args.tol,
            msh_output=args.msh_output,
        )

    calc_main(run, parser)


if __name__ == "__main__":
    main()
