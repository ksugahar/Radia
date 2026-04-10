"""
3D FEM-SIBC solver for Cubit panel (IHFEMDialog).

Two approaches for workpiece:
  Hole (preferred): workpiece subtracted from mesh, SIBC on hole boundary "sibc".
  Interface (legacy): workpiece meshed as volume, SIBC on internal "wp_surface".

User workflow:
  1. Write .jou to create geometry + mesh + named blocks
  2. Run journal in Cubit panel
  3. Click Solve -> panel saves temp .cub5 -> this script runs as subprocess

Cubit blocks (user sets in .jou):
  coil        - volume block, J source
  air         - volume block (includes where workpiece was, for hole approach)
  sibc        - sideset on hole boundary (hole approach, preferred)
  wp_surface  - surface on air/workpiece interface (legacy interface approach)
  kelvin      - volume block (optional, Periodic Kelvin)
  source      - surface block, T0=1 (gap face, optional)
  sink        - surface block, T0=0 (gap face, optional)

IMPORTANT: NGSolve must be imported BEFORE cubit.
"""

import argparse
import json
import math
import os
import sys
import time

# Pre-import scipy BEFORE Cubit (Cubit bundles broken scipy)
try:
    import scipy  # noqa: F401
    import scipy.interpolate  # noqa: F401
except ImportError:
    pass

import numpy as np

# Shared utilities
_this_dir = os.path.dirname(os.path.abspath(__file__))
if _this_dir not in sys.path:
    sys.path.insert(0, _this_dir)

from calc_common import (MU_0, NU_0, setup_paths, setup_cubit, export_mesh,
                          add_periodic_kelvin, detect_kelvin_offset,
                          create_esim_solver, get_bh_curve,
                          compute_T0_source, compute_J_theta,
                          progress, calc_main)


def _log(msg):
    """Write progress to stderr (panel reads these)."""
    progress("FEM", msg)


def solve_fem(cub5_file="", vol_file="", order=2, fes_order=1,
              frequency=7000, sigma=2e6,
              impedance_model="sibc", mu_r=100.0,
              bh_file=None, material="steel",
              I_total=1.0, a_coil=0.003, R_wp=0.010,
              max_iter=15, tol=1e-3, relax=0.5,
              solver="pardiso", reg=1e-6, shift_eps=1e-6,
              nthreads=0,
              msh_output=""):
    """3D FEM-SIBC solver for Cubit mesh.

    Args:
        cub5_file: Cubit .cub5 file (temp file from panel dialog)
        order: Mesh curve order (1-3)
        fes_order: HCurl polynomial order (1-3)
        frequency: Operating frequency [Hz]
        sigma: Workpiece conductivity [S/m]
        impedance_model: "sibc" (linear) or "esim" (nonlinear)
        mu_r: Relative permeability (SIBC mode)
        bh_file: BH curve file path (ESIM mode)
        material: "steel", "copper", "aluminum"
        I_total: Coil current [A]
        a_coil: Coil wire radius [m] (for J0*e_theta fallback)
        R_wp: Workpiece radius [m] (for ESIM half-thickness)
        max_iter: Max Karl iterations
        tol: Z_s convergence tolerance
        relax: Under-relaxation (0-1)
        solver: "pardiso" (direct), "bddc" (iterative), "iccg" (shifted IC+CG),
                "ams" (Compact AMS+COCR)
        reg: Gauge regularization parameter (add reg*nu0*u*v*dx to system)
        shift_eps: Shifted preconditioner eps (AMS/ICCG, add eps*nu*u*v*dx to prec)
        nthreads: TaskManager thread count (0=auto, 1=single thread)
        msh_output: Optional GMSH .msh output path

    Returns:
        dict with P_total, L, Z_s, H_t_rms, etc.
    """
    # NGSolve must be imported BEFORE Cubit
    import ngsolve  # noqa: F401
    from ngsolve import (H1, HCurl, Periodic, BilinearForm, LinearForm,
                         GridFunction, Integrate, Conj, curl, dx, ds, CF,
                         BND, VOL, TaskManager, Preconditioner)
    from ngsolve import x, y, z, sqrt, IfPos

    setup_paths()

    t_total_start = time.perf_counter()

    # ============================================================
    # Step 1: Load mesh (.vol direct or Cubit .cub5)
    # ============================================================
    _log("MESH:loading")
    t0 = time.perf_counter()

    if vol_file and os.path.exists(vol_file):
        # Direct .vol path: skip Cubit entirely
        from ngsolve import Mesh as NGMesh
        mesh = NGMesh(vol_file)
        _log(f"MESH:loaded .vol directly ({vol_file})")
    elif cub5_file and os.path.exists(cub5_file):
        cubit = setup_cubit(cub5_file)
        if cubit is None:
            return {"error": "Cubit not available"}
        mesh = export_mesh(cubit, order=order)
    else:
        return {"error": f"Either --vol or --cub5 required"}

    # Kelvin: detect offset, estimate radius, add periodic identification
    # Only attempt if "kelvin" material exists in mesh
    a_kelvin = 0.060
    kelvin_center = (0, 0, 0)
    mats_pre = mesh.GetMaterials()
    if "kelvin" in mats_pre:
        kelvin_center = detect_kelvin_offset(mesh)
        try:
            # Estimate a_kelvin from max radius of kelvin_int boundary vertices
            kelvin_verts = set()
            for el in mesh.Elements(VOL):
                if "kelvin" in str(mats_pre[el.nr]).lower():
                    for v in el.vertices:
                        kelvin_verts.add(v.nr)
            if kelvin_verts:
                kc = np.array(kelvin_center)
                coords = np.array([mesh.vertices[v].point for v in kelvin_verts])
                dists = np.linalg.norm(coords - kc[np.newaxis, :], axis=1)
                a_kelvin = float(np.min(dists))  # inner radius of kelvin shell
        except Exception:
            pass
        if add_periodic_kelvin(mesh, kelvin_center):
            mesh = ngsolve.Mesh(mesh.ngmesh)
            _log(f"PERIODIC:added (a={a_kelvin:.4f}, offset={kelvin_center})")

    t_mesh = time.perf_counter() - t0

    materials = mesh.GetMaterials()
    boundaries = mesh.GetBoundaries()
    ne = mesh.GetNE(VOL)
    _log(f"MESH:done ({ne} elems, {t_mesh:.1f}s)")

    # ============================================================
    # Step 2: Source current (T0 or J0*e_theta)
    # ============================================================
    has_source = "source" in boundaries
    has_sink = "sink" in boundaries
    has_kelvin = "kelvin" in materials

    # SIBC boundary: "sibc" (hole approach, preferred) or "wp_surface" (legacy)
    if "sibc" in boundaries:
        sibc_bnd = "sibc"
        has_wp = True
    elif "wp_surface" in boundaries:
        sibc_bnd = "wp_surface"
        has_wp = True
    else:
        sibc_bnd = ""
        has_wp = False

    # Hole approach: workpiece NOT in materials (subtracted from mesh)
    is_hole = has_wp and "workpiece" not in materials
    if has_wp:
        _log(f"SIBC:boundary={sibc_bnd}, approach={'hole' if is_hole else 'interface'}")

    if has_source and has_sink:
        _log("SOURCE:T0 source/sink technique")
        J_source, gf_T0 = compute_T0_source(mesh, fes_order, I_total)
    else:
        _log("SOURCE:J0*e_theta (torus fallback)")
        J_source = compute_J_theta(I_total, a_coil)
        gf_T0 = None

    # ============================================================
    # Step 3: Material properties
    # ============================================================
    # Kelvin weight: nu = nu0 * (a/r')^2 in kelvin domain
    # r' = distance from exterior sphere center (kelvin_center)
    # All other domains: nu = nu0 (workpiece treated as air for SIBC)
    kx, ky, kz = kelvin_center
    nu_dict = {}
    for m in materials:
        if "kelvin" in m.lower():
            dx_k = x - kx
            dy_k = y - ky
            dz_k = z - kz
            rp_sq = dx_k * dx_k + dy_k * dy_k + dz_k * dz_k + 1e-20
            kelvin_fac = a_kelvin**2 / rp_sq  # (a/r')^2
            nu_dict[m] = NU_0 * kelvin_fac
        else:
            nu_dict[m] = NU_0
    nu_cf = mesh.MaterialCF(nu_dict, default=NU_0)
    if has_kelvin:
        _log(f"KELVIN:a={a_kelvin}, center=({kx},{ky},{kz})")

    # ============================================================
    # Step 4: ESIM / SIBC solver
    # ============================================================
    omega = 2 * math.pi * frequency

    if impedance_model == "esim":
        esim = create_esim_solver(
            material=material, frequency=frequency, sigma=sigma,
            half_thickness=R_wp, geometry='cylinder', bh_file=bh_file)
        Z_s = esim.solve(5.0)['Z']
    else:
        # Linear SIBC: Z_s = (1+j) * rho / delta
        rho = 1.0 / sigma
        mu_eff = MU_0 * mu_r
        delta = math.sqrt(2 * rho / (omega * mu_eff)) if omega > 0 else 1e10
        Z_s = complex(1, 1) * rho / delta
        esim = None

    _log(f"SIBC:Z_s={Z_s:.4e}, model={impedance_model}")

    # ============================================================
    # Step 5: FE space + Karl iteration
    # ============================================================
    # Periodic Kelvin: GND vertex at exterior sphere center = Dirichlet
    # No "outer" boundary needed (periodic BC handles far field)
    dirichlet_bnd = ""
    if "GND" in boundaries:
        dirichlet_bnd = "GND"
    elif "outer" in boundaries:
        dirichlet_bnd = "outer"

    base_fes = HCurl(mesh, order=fes_order, complex=True,
                     dirichlet=dirichlet_bnd)
    # Wrap with Periodic if Kelvin identification exists
    if has_kelvin:
        fes = Periodic(base_fes)
        _log("FES:Periodic HCurl (Kelvin)")
    else:
        fes = base_fes
    u, v = fes.TnT()
    ndof = fes.ndof
    _log(f"FES:ndof={ndof}")

    _log(f"SOLVER:{solver}, reg={reg}, shift_eps={shift_eps}, threads={nthreads}")

    # TaskManager thread count
    if nthreads > 0:
        import ngsolve
        ngsolve.SetNumThreads(nthreads)

    # RHS (constant across Karl iterations)
    f_lf = LinearForm(fes)
    f_lf += J_source * v * dx("coil")
    f_lf.Assemble()

    # SIBC surface area for H_t estimation
    if has_wp:
        wp_region = mesh.Boundaries(sibc_bnd)
        A_wp = Integrate(CF(1), mesh, BND, definedon=wp_region).real
    else:
        A_wp = 2 * math.pi * R_wp * 0.020  # rough estimate

    # Karl iteration
    gfu = GridFunction(fes)
    history = []

    for iteration in range(max_iter):
        t0_iter = time.perf_counter()

        # Robin coefficient:
        #   Hole approach: +jw/Z_s (SIBC on external boundary of air domain)
        #   Interface: sign depends on which side NGSolve evaluates from
        #              +jw/Z_s works for both (tested with tangential projection)
        robin = 1j * omega / Z_s if has_wp else 0

        # System bilinear form
        a_bf = BilinearForm(fes)
        a_bf += nu_cf * curl(u) * curl(v) * dx(bonus_intorder=4)
        if reg > 0:
            a_bf += reg * NU_0 * u * v * dx  # gauge regularization

        if has_wp and abs(robin) > 0:
            a_bf += robin * u.Trace() * v.Trace() * ds(sibc_bnd)

        # Solver-specific setup
        if solver == "ams":
            # Shifted preconditioner (AMS+COCR)
            a_shifted = BilinearForm(fes)
            a_shifted += nu_cf * curl(u) * curl(v) * dx(bonus_intorder=4)
            a_shifted += shift_eps * NU_0 * u * v * dx  # shift for prec only
            if reg > 0:
                a_shifted += reg * NU_0 * u * v * dx
            if has_wp and abs(robin) > 0:
                a_shifted += robin * u.Trace() * v.Trace() * ds(sibc_bnd)

            from ngsolve.la import CompactAMSPreconditioner, COCRSolver
            pre_ams = CompactAMSPreconditioner(a_shifted, fes)
            a_shifted.Assemble()
            a_bf.Assemble()

            with TaskManager():
                cocr = COCRSolver(a_bf.mat, pre_ams, maxiter=500, tol=1e-8,
                                  freedofs=fes.FreeDofs())
                gfu.vec.data = cocr * f_lf.vec
        elif solver == "iccg":
            # Shifted ICCG (IC preconditioned CG, single-thread efficient)
            a_bf.Assemble()

            import sparsesolv_ngsolve as ssn
            iccg = ssn.SparseSolvSolver(
                a_bf.mat, method="ICCG",
                freedofs=fes.FreeDofs(),
                tol=1e-8, maxiter=500, shift=shift_eps,
                save_best_result=True, printrates=False,
                use_abmc=True, abmc_block_size=4, abmc_num_colors=4)
            iccg.auto_shift = True
            gfu.vec.data = iccg * f_lf.vec
        elif solver == "bddc":
            pre = Preconditioner(a_bf, "bddc")
            a_bf.Assemble()
            with TaskManager():
                from ngsolve import solvers
                solvers.BVP(bf=a_bf, lf=f_lf, gf=gfu,
                            pre=pre, maxsteps=500, tol=1e-8)
        else:
            # pardiso (direct)
            a_bf.Assemble()
            with TaskManager():
                gfu.vec.data = a_bf.mat.Inverse(
                    fes.FreeDofs(), inverse="pardiso") * f_lf.vec

        t_solve_iter = time.perf_counter() - t0_iter

        # H_t from SIBC relation: H_t = |jw/Z_s| * |A_t|
        #
        # Tangential projection: |A_t|^2 = |A|^2 - (A . n)^2
        # n = specialcf.normal(3) = outward from mesh domain (air)
        #
        # Energy-balance P: P_sibc = P_input - P_curlcurl
        #   P_input = Re(int J . A* dV)  (volume integral, no BND issues)
        #   P_curlcurl = int nu |curl(A)|^2 dV + reg*nu0*|A|^2 dV
        if has_wp and abs(Z_s) > 0:
            from ngsolve import specialcf
            n_bnd = specialcf.normal(3)

            # |A|^2
            A_sq = sum(gfu[i].real * gfu[i].real +
                       gfu[i].imag * gfu[i].imag for i in range(3))
            # (A . n)^2 = (Re(A).n)^2 + (Im(A).n)^2
            Adn_re = sum(gfu[i].real * n_bnd[i] for i in range(3))
            Adn_im = sum(gfu[i].imag * n_bnd[i] for i in range(3))
            An_sq = Adn_re * Adn_re + Adn_im * Adn_im
            At_sq = A_sq - An_sq  # |A_t|^2

            int_At2 = Integrate(At_sq, mesh, BND,
                                definedon=wp_region).real
            At_rms = math.sqrt(max(int_At2, 0) / max(A_wp, 1e-30))
            H_t_bnd = abs(1j * omega / Z_s) * At_rms

            # BND integral with tangential projection works for both approaches:
            #   Hole: single-side trace (air external boundary), fully correct
            #   Interface: tangential projection removes normal artifacts
            # Energy-balance is unreliable for interface (Robin leaks into workpiece)
            H_t_rms = H_t_bnd
        else:
            H_t_rms = 5.0  # fallback

        # Update Z_s
        Z_s_old = Z_s
        if esim is not None:
            sol_new = esim.solve(max(float(H_t_rms), 1e-3))
            Z_s = relax * sol_new['Z'] + (1 - relax) * Z_s_old
        # else: linear SIBC, Z_s constant (no iteration needed)

        dZ = abs(Z_s - Z_s_old) / max(abs(Z_s_old), 1e-30)
        history.append({
            'iteration': iteration,
            'Z_s_abs': abs(Z_s),
            'H_t_rms': float(H_t_rms),
            'dZ': float(dZ),
            't_solve': t_solve_iter,
        })
        _log(f"ITER:{iteration} |Z_s|={abs(Z_s):.4e} H_t={H_t_rms:.2f} "
             f"dZ={dZ:.4e} t={t_solve_iter:.1f}s")

        if dZ < tol and iteration > 0:
            _log(f"CONVERGED:iter={iteration}")
            break

        if esim is None:
            break  # linear SIBC: single iteration

    # ============================================================
    # Step 6: Post-process
    # ============================================================
    # Time-averaged power: P = 0.5 * Re(Z_s) * H_t_rms^2 * A_wp
    if has_wp:
        P_total = 0.5 * Z_s.real * H_t_rms**2 * A_wp
        Q_total = 0.5 * Z_s.imag * H_t_rms**2 * A_wp
    else:
        P_total = 0.0
        Q_total = 0.0

    # Inductance from magnetic energy
    W_mag = Integrate(
        0.5 * nu_cf * curl(gfu) * Conj(curl(gfu)), mesh, order=10).real
    L = 2 * W_mag / I_total**2

    t_total = time.perf_counter() - t_total_start
    _log(f"DONE:P={P_total:.4e} L={L*1e9:.2f}nH t={t_total:.1f}s")

    # ============================================================
    # Step 7: GMSH export
    # ============================================================
    gmsh_file = ""
    if msh_output:
        try:
            from gmsh_post_export import GmshPostExport

            # Volume field: |B|
            curl_A = curl(gfu)
            B_mag = sqrt(sum(curl_A[i].real * curl_A[i].real +
                             curl_A[i].imag * curl_A[i].imag
                             for i in range(3)))

            post = GmshPostExport(mesh, boundary=False)
            # Compute per-vertex |B| values
            fes_h1 = ngsolve.H1(mesh, order=1)
            gf_B = GridFunction(fes_h1)
            gf_B.Set(B_mag)
            node_B = np.array([gf_B(mesh(*mesh.vertices[v.nr].point))
                               for v in mesh.vertices])
            post.add_field("|B|", node_B, ncomp=1)
            post.write(msh_output)
            gmsh_file = msh_output
            _log(f"GMSH:{msh_output}")
        except Exception as e:
            _log(f"GMSH_ERROR:{e}")

    # ============================================================
    # Step 8: Result JSON
    # ============================================================
    delta_skin = math.sqrt(2.0 / (omega * MU_0 * (mu_r if esim is None else 100) * sigma)) if omega > 0 else 0

    result = {
        "P_total": float(P_total),
        "Q_total": float(Q_total),
        "L": float(L),
        "Z_s": str(Z_s),
        "H_t_rms": float(H_t_rms),
        "delta": float(delta_skin),
        "ndof": ndof,
        "ne": ne,
        "iterations": len(history),
        "t_mesh": round(t_mesh, 2),
        "t_total": round(t_total, 2),
        "frequency": frequency,
        "sigma": sigma,
        "impedance_model": impedance_model,
        "source_type": "T0" if gf_T0 is not None else "J_theta",
        "has_kelvin": has_kelvin,
        "msh_file": gmsh_file,
    }
    return result


def main():
    parser = argparse.ArgumentParser(
        description="3D FEM-SIBC with source/sink + optional Kelvin")
    parser.add_argument("--cub5", default="", help="Cubit .cub5 file")
    parser.add_argument("--vol", default="", help="Netgen .vol file (skip Cubit)")
    parser.add_argument("--order", type=int, default=2,
                        help="Curve order (1-3)")
    parser.add_argument("--fes-order", type=int, default=1,
                        help="HCurl polynomial order (1-3)")
    parser.add_argument("--frequency", type=float, default=7000,
                        help="Frequency [Hz]")
    parser.add_argument("--sigma", type=float, default=2e6,
                        help="Conductivity [S/m]")
    parser.add_argument("--impedance", default="sibc",
                        choices=["sibc", "esim"],
                        help="Impedance model")
    parser.add_argument("--mu-r", type=float, default=100,
                        help="mu_r (SIBC mode)")
    parser.add_argument("--bh-file", default="",
                        help="BH curve file (ESIM mode)")
    parser.add_argument("--material", default="steel",
                        choices=["steel", "copper", "aluminum"])
    parser.add_argument("--current", type=float, default=1.0,
                        help="Coil current [A]")
    parser.add_argument("--a-coil", type=float, default=0.003,
                        help="Coil wire radius [m] (for J_theta fallback)")
    parser.add_argument("--r-wp", type=float, default=0.010,
                        help="Workpiece radius [m]")
    parser.add_argument("--max-iter", type=int, default=15,
                        help="Max Karl iterations")
    parser.add_argument("--solver", default="pardiso",
                        choices=["pardiso", "bddc", "iccg", "ams"],
                        help="pardiso (direct), bddc (iterative), "
                             "iccg (shifted IC+CG), ams (Compact AMS+COCR)")
    parser.add_argument("--reg", type=float, default=1e-6,
                        help="Gauge regularization (add reg*nu0*u*v*dx)")
    parser.add_argument("--shift-eps", type=float, default=1e-6,
                        help="Shifted preconditioner eps (AMS/ICCG)")
    parser.add_argument("--nthreads", type=int, default=0,
                        help="TaskManager threads (0=auto, 1=single)")
    parser.add_argument("--msh-output", default="",
                        help="GMSH .msh output path")
    parser.add_argument("--output", default="",
                        help="JSON output file")

    def run(args):
        return solve_fem(
            cub5_file=args.cub5,
            vol_file=args.vol,
            order=args.order,
            fes_order=args.fes_order,
            frequency=args.frequency,
            sigma=args.sigma,
            impedance_model=args.impedance,
            mu_r=args.mu_r,
            bh_file=args.bh_file if args.bh_file else None,
            material=args.material,
            I_total=args.current,
            a_coil=args.a_coil,
            R_wp=args.r_wp,
            max_iter=args.max_iter,
            solver=args.solver,
            reg=args.reg,
            shift_eps=args.shift_eps,
            nthreads=args.nthreads,
            msh_output=args.msh_output,
        )

    calc_main(run, parser)


if __name__ == "__main__":
    main()
