"""
3D FEM-SIBC solver for Radia-NGSolve panel.

Usage:
    python calc_fem_kelvin.py --vol model.vol --frequency 7000 --sigma 2e6

Input: Netgen .vol file with material/boundary labels:
  coil        - volume material, J source
  air         - volume material (includes where workpiece was, for hole approach)
  sibc        - boundary on hole surface (hole approach, preferred)
  wp_surface  - boundary on air/workpiece interface (legacy interface approach)
  kelvin      - volume material (optional, Periodic Kelvin)
  source      - boundary, T0=1 (gap face, optional)
  sink        - boundary, T0=0 (gap face, optional)

Two approaches for workpiece:
  Hole (preferred): workpiece subtracted from mesh, SIBC on hole boundary "sibc".
  Interface (legacy): workpiece meshed as volume, SIBC on internal "wp_surface".

Outputs JSON to stdout.
"""

import argparse
import json
import math
import os
import sys
import time

import numpy as np

# Shared utilities
_this_dir = os.path.dirname(os.path.abspath(__file__))
if _this_dir not in sys.path:
    sys.path.insert(0, _this_dir)

from calc_common import (MU_0, NU_0, setup_paths,
                          add_periodic_kelvin, detect_kelvin_offset,
                          create_esim_solver, get_bh_curve,
                          compute_T0_source, compute_J_theta,
                          progress, calc_main)


def _log(msg):
    """Write progress to stderr (panel reads these)."""
    progress("FEM", msg)


def solve_fem(vol_file="", fes_order=1,
              frequency=7000, sigma=2e6,
              impedance_model="sibc", mu_r=100.0,
              bh_file=None, material="steel",
              I_total=1.0, half_thickness=0.005,
              max_iter=15, tol=1e-3, relax=0.5,
              solver="pardiso", reg=1e-6, shift_eps=1e-6,
              nthreads=0,
              msh_output=""):
    """3D FEM-SIBC solver for .vol mesh.

    Args:
        vol_file: Netgen .vol file (with material/boundary labels)
        fes_order: HCurl polynomial order (1-3)
        frequency: Operating frequency [Hz]
        sigma: Workpiece conductivity [S/m]
        impedance_model: "sibc" (linear) or "esim" (nonlinear)
        mu_r: Relative permeability (SIBC mode)
        bh_file: BH curve file path (ESIM mode)
        material: "steel", "copper", "aluminum"
        I_total: Coil current [A]
        half_thickness: Workpiece characteristic radius / half
            thickness [m]. R in the Dowell tanh formula for the
            linear SIBC path AND cell-problem domain length for
            the ESIM nonlinear path.
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
    # Step 1: Load mesh (.vol)
    # ============================================================
    _log("MESH:loading")
    t0 = time.perf_counter()

    if not vol_file or not os.path.exists(vol_file):
        return {"error": f"--vol is required (got: {vol_file!r})"}

    from ngsolve import Mesh as NGMesh
    mesh = NGMesh(vol_file)
    _log(f"MESH:loaded {vol_file}")

    # POLICY: FES order must not exceed mesh curve order.
    # curve_order=1 + fes_order>=2 produces large geometry/basis mismatch error
    # because high-order HCurl basis functions assume curved edges that
    # don't exist on a piecewise-linear mesh. Inductance diverges with order.
    try:
        curve_order = mesh.ngmesh.GetCurveOrder()
    except Exception:
        curve_order = 1
    if fes_order > curve_order:
        return {
            "error": f"fes_order ({fes_order}) > mesh curve_order ({curve_order}). "
                     f"Re-export .vol with 'radia_export netgen ... order {fes_order}' "
                     f"or use --fes-order {curve_order}."
        }
    _log(f"MESH:curve_order={curve_order}, fes_order={fes_order}")

    # Kelvin: detect offset, measure radius, require periodic identification
    # to be already embedded in the .vol by the C++ export. No fallback —
    # if Kelvin is present in the mesh but periodic ID is missing, abort.
    # (A Kelvin-transformed mesh without periodic identification gives a
    # silently wrong answer that looks like a Dirichlet truncation result.)
    kelvin_center = (0, 0, 0)
    a_kelvin = 0.0
    mats_pre = mesh.GetMaterials()
    has_kelvin_periodic = False
    if "kelvin" in mats_pre:
        kelvin_center = detect_kelvin_offset(mesh)

        # Estimate a_kelvin = sphere radius from MAX vertex distance from
        # the kelvin centroid. The 2-sphere model uses a solid sphere
        # (not a shell), so np.min would return ~0 (interior vertex), which
        # is wrong. The correct value is np.max = sphere radius R.
        kelvin_verts = set()
        for el in mesh.Elements(VOL):
            if el.mat == "kelvin":
                for v in el.vertices:
                    kelvin_verts.add(v.nr)
        if not kelvin_verts:
            return {"error":
                    "Mesh has 'kelvin' material in the materials list but "
                    "no volume elements are tagged 'kelvin'. The .vol export "
                    "is inconsistent — re-run the Cubit .jou and re-export."}
        kc = np.array(kelvin_center)
        coords = np.array([mesh.vertices[v].point for v in kelvin_verts])
        dists = np.linalg.norm(coords - kc[np.newaxis, :], axis=1)
        a_kelvin = float(np.max(dists))  # sphere radius R

        # Periodic identification MUST come from the C++ export. The
        # Python add_periodic_kelvin path was a debug crutch that produced
        # subtly different identification tolerances and is now removed.
        n_ident = mesh.ngmesh.GetNrIdentifications()
        if n_ident == 0:
            return {"error":
                    "Mesh has 'kelvin' material but no periodic "
                    "identification in the .vol file. The C++ export "
                    "(radia_export netgen ... order N) is responsible for "
                    "writing the inner-sphere/outer-sphere triangle pairs "
                    "into the .vol. Without it the FEM solve falls back to "
                    "Dirichlet truncation and silently produces a wrong "
                    "answer for the open-boundary problem. "
                    "Re-export the .vol with the latest plugin and ensure "
                    "the .jou uses the 'subtract -> imprint -> merge' "
                    "pattern for the kelvin sphere."}
        has_kelvin_periodic = True
        _log(f"PERIODIC:from_vol ({n_ident} ident(s), a={a_kelvin:.4f})")

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

    # SIBC boundary: must be named "sibc" by .jou convention.
    # The legacy "wp_surface" name is no longer accepted (single supported
    # convention, no fallback).
    if "sibc" in boundaries:
        sibc_bnd = "sibc"
        has_wp = True
    else:
        sibc_bnd = ""
        has_wp = False

    # Hole approach: workpiece NOT in materials (subtracted from mesh)
    is_hole = has_wp and "workpiece" not in materials
    if has_wp:
        _log(f"SIBC:boundary={sibc_bnd}, approach={'hole' if is_hole else 'interface'}")

    if not (has_source and has_sink):
        return {"error":
                "Mesh has no 'source' and/or 'sink' boundary labels. "
                "These are required for the T0 source/sink technique. "
                f"Available boundaries: {sorted(set(boundaries))}. "
                "Fix the .jou: define sidesets named 'source' and 'sink' "
                "on the coil terminal faces and re-export."}
    _log("SOURCE:T0 source/sink technique")
    J_source, gf_T0 = compute_T0_source(mesh, fes_order, I_total)

    # ============================================================
    # Step 3: Material properties
    # ============================================================
    # Kelvin weight for HCurl A-formulation:  nu_eff = (r'/R)^2 * nu0
    #
    # Physical rule (Kelvin transform under inversion r -> R^2/r'):
    #   - mu, eps, sigma -> infinity at r'=0 (image of physical infinity)
    #   - nu = 1/mu, rho_e = 1/sigma -> 0 at r'=0
    # The PDE is unchanged; only the constitutive coefficient is distorted.
    #
    # r' = distance from exterior sphere center (kelvin_center)
    # R = a_kelvin (sphere radius, same for both spheres)
    # All other domains: nu = nu0 (workpiece treated as air for SIBC)
    kx, ky, kz = kelvin_center
    nu_dict = {}
    for m in materials:
        if "kelvin" in m.lower():
            dx_k = x - kx
            dy_k = y - ky
            dz_k = z - kz
            rp_sq = dx_k * dx_k + dy_k * dy_k + dz_k * dz_k + 1e-20
            kelvin_fac = rp_sq / a_kelvin**2  # (r'/R)^2
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

    if sigma <= 0 or omega <= 0:
        # DC or no conductor: pure inductance, no SIBC
        Z_s = complex(0, 0)
        esim = None
        has_wp = False
    elif impedance_model == "esim":
        esim = create_esim_solver(
            material=material, frequency=frequency, sigma=sigma,
            half_thickness=half_thickness, geometry='cylinder',
            bh_file=bh_file)
        Z_s = esim.solve(5.0)['Z']
    else:
        # Linear SIBC: Dowell tanh formula with R = half_thickness.
        #   Z_s = (rho / R) * (1+j)*xi * tanh((1+j)*xi)
        #   xi  = R / delta,  delta = sqrt(2 rho / (omega mu_eff))
        # Reduces to the thick-limit Z_s = (1+j) rho / delta when
        # xi >> 1 (typical Cu workpiece at 50 kHz: xi ~ 30) but gives
        # the correct lower magnitude when the workpiece is on the
        # order of one skin depth (e.g. magnetic materials at high
        # frequency, where mu_r raises delta).
        rho = 1.0 / sigma
        mu_eff = MU_0 * mu_r
        if omega <= 0 or half_thickness <= 0:
            delta = 1e10
            Z_s = complex(0, 0)
        else:
            delta = math.sqrt(2 * rho / (omega * mu_eff))
            xi = half_thickness / delta
            gamma_a = complex(1, 1) * xi
            try:
                Z_s = ((rho / half_thickness)
                       * gamma_a * np.tanh(gamma_a))
            except OverflowError:
                # tanh saturates to 1 for large xi -> thick-limit
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
    elif not has_kelvin and "outer" in boundaries:
        dirichlet_bnd = "outer"
    # 2-sphere Kelvin: Periodic BC handles far field, no Dirichlet on kelvin_ext.
    # GND (vertex at kelvin sphere center) provides uniqueness if present.

    is_dc = (omega <= 0 or sigma <= 0)
    # POLICY: nograds=True is REQUIRED for HCurl curl-curl at any order.
    # Reference: NGSolve Maxwell tutorial unit-2.4
    #   HCurl(mesh, order=p, dirichlet="outer", nograds=True)
    # The HCurl high-order basis includes gradients of H1 basis functions
    # which form the curl-curl gauge kernel. nograds=True removes them
    # so the system is well-conditioned for any p (including p>=2).
    # complex flag follows AC vs DC (eddy term needs complex arithmetic).
    use_complex = not is_dc
    base_fes = HCurl(mesh, order=fes_order,
                     complex=use_complex,
                     nograds=True,
                     dirichlet=dirichlet_bnd)
    # 2-sphere Kelvin: Periodic BC couples interior and exterior sphere DOFs
    if has_kelvin and has_kelvin_periodic:
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
        # Workpiece is not in the .vol — use a rough placeholder
        # (the H_t estimate is only used for the iteration printout
        # in the no-workpiece case, not for the inductance result).
        A_wp = 4 * math.pi * half_thickness ** 2

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
        # Gauge regularization. nograds=True removes the H1 gradient
        # subspace but residual harmonic forms (cohomology of nontrivial
        # topology) still leave the curl-curl operator singular. Add a
        # tiny mass term to fix the kernel.  This is needed in BOTH DC
        # and AC mode unless an SIBC Robin term covers the air domain.
        #
        # IMPORTANT: in the Kelvin domain nu_cf = (r'/R)^2 * nu_0 vanishes
        # near r'=0, so a constant NU_0 mass term would dominate the curl
        # term there and corrupt the magnetic energy integral. Restrict the
        # regularization to non-kelvin materials when Periodic Kelvin is on;
        # otherwise apply it everywhere.
        if reg > 0:
            if has_kelvin and has_kelvin_periodic:
                non_kelvin_mats = [m for m in materials if "kelvin" not in m.lower()]
                if non_kelvin_mats:
                    a_bf += reg * NU_0 * u * v * dx("|".join(non_kelvin_mats))
            else:
                a_bf += reg * NU_0 * u * v * dx

        if has_wp and abs(robin) > 0:
            a_bf += robin * u.Trace() * v.Trace() * ds(sibc_bnd)

        # Solver-specific setup
        if solver == "ams":
            # Compact AMS preconditioner + COCR solver (AC eddy-current only).
            # Reference: src/ext/sparsesolv/examples/hiruma/bench_compact_ams.py
            #
            # Constraints:
            #   - fes_order = 1 (CompactAMS auxiliary space is order-1 HCurl)
            #   - nograds = True (H1 auxiliary space matches vertex count)
            #   - AC mode (sigma > 0): AMS is designed for K + jw*sigma*M
            #   - No Periodic BC (DOF dimensions don't match)
            if fes_order > 1:
                return {"error": f"solver=ams requires fes_order=1 (got {fes_order})"}
            if is_dc:
                return {"error": "solver=ams requires AC (frequency>0, sigma>0). "
                                 "For DC use solver=bddc or pardiso."}
            if has_kelvin and has_kelvin_periodic:
                return {"error": "solver=ams is not supported with Periodic Kelvin BC."}

            # Rebuild the system fes/bilinear form with nograds=True
            # (overrides the use_complex code path which uses nograds=False)
            import sparsesolv_ngsolve as ssn

            fes_ams = HCurl(mesh, order=1, complex=True, nograds=True,
                            dirichlet=dirichlet_bnd)
            u_a, v_a = fes_ams.TnT()
            a_bf = BilinearForm(fes_ams)
            a_bf += nu_cf * curl(u_a) * curl(v_a) * dx(bonus_intorder=4)
            a_bf += reg * NU_0 * u_a * v_a * dx
            a_bf += 1j * omega * sigma * u_a * v_a * dx  # AC eddy term
            if has_wp and abs(robin) > 0:
                a_bf += robin * u_a.Trace() * v_a.Trace() * ds(sibc_bnd)
            a_bf.Assemble()

            f_ams = LinearForm(fes_ams)
            f_ams += J_source * v_a * dx("coil")
            f_ams.Assemble()
            f_lf = f_ams  # used by H_t evaluation later

            # Real auxiliary fes for the AMS preconditioner
            fes_real = HCurl(mesh, order=1, complex=False, nograds=True,
                             dirichlet=dirichlet_bnd)
            ur, vr = fes_real.TnT()
            a_real = BilinearForm(fes_real)
            a_real += nu_cf * curl(ur) * curl(vr) * dx(bonus_intorder=4)
            a_real += shift_eps * NU_0 * ur * vr * dx
            a_real += abs(omega) * sigma * ur * vr * dx
            if has_wp and abs(robin) > 0:
                a_real += abs(robin) * ur.Trace() * vr.Trace() * ds(sibc_bnd)
            a_real.Assemble()

            grad_mat, fes_h1 = fes_real.CreateGradient()
            nv = mesh.nv
            coord_x = [0.0] * nv; coord_y = [0.0] * nv; coord_z = [0.0] * nv
            for i in range(nv):
                pt = mesh.vertices[i].point
                coord_x[i] = pt[0]; coord_y[i] = pt[1]; coord_z[i] = pt[2]

            pre_ams = ssn.ComplexCompactAMSPreconditioner(
                a_real_mat=a_real.mat, grad_mat=grad_mat,
                freedofs=fes_real.FreeDofs(),
                coord_x=coord_x, coord_y=coord_y, coord_z=coord_z,
                ndof_complex=fes_ams.ndof, cycle_type=1, print_level=0)

            gfu = GridFunction(fes_ams)
            with TaskManager():
                cocr = ssn.COCRSolver(a_bf.mat, pre_ams,
                                      freedofs=fes_ams.FreeDofs(),
                                      maxiter=500, tol=1e-8, printrates=False)
                gfu.vec.data = cocr * f_ams.vec
            fes = fes_ams  # downstream code uses `fes`
            _log(f"AMS:iters={cocr.iterations}")
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
        # Only meaningful when there is a workpiece. For the inductance-
        # only path (no workpiece) the Karl loop already breaks at the
        # esim-is-None check below; H_t_rms is not used downstream.
        H_t_rms = None
        if has_wp:
            if not (abs(Z_s) > 0):
                raise RuntimeError(
                    "Z_s == 0 inside Karl iteration with workpiece present. "
                    "This is an uninitialized state — check the SIBC setup "
                    "(material, sigma, frequency, half-thickness, esim) "
                    "before the iteration starts.")
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
            # BND integral with tangential projection works for both approaches:
            #   Hole: single-side trace (air external boundary), fully correct
            #   Interface: tangential projection removes normal artifacts
            # Energy-balance is unreliable for interface (Robin leaks into workpiece)
            H_t_rms = abs(1j * omega / Z_s) * At_rms

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
            'H_t_rms': float(H_t_rms) if H_t_rms is not None else None,
            'dZ': float(dZ),
            't_solve': t_solve_iter,
        })
        h_t_str = f"{H_t_rms:.2f}" if H_t_rms is not None else "n/a"
        _log(f"ITER:{iteration} |Z_s|={abs(Z_s):.4e} H_t={h_t_str} "
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
    # Step 7: GMSH export — B vector + J vector + companion .msh.opt
    # ============================================================
    # Match the BEM viz convention: VECTORS only (GMSH shows the
    # magnitude under the vector view automatically). The companion
    # .msh.opt hides volume mesh elements so arrows are visible.
    gmsh_file = ""
    if msh_output:
        _log("GMSH:start")
        try:
            from gmsh_post_export import GmshPostExport

            curl_A = curl(gfu)

            def _eval_vec_at_vertices(cf):
                """Evaluate a vector CF at every mesh vertex.

                Returns an (nv, 3) numpy array. Uses NGSolve's
                MeshPoint evaluation directly because building a
                vector H1 GridFunction + Set() can fail on the FEM
                volume mesh when the source is a curl of a complex
                AC GridFunction (NGSolve internals try to allocate a
                complex vector space). Per-vertex point evaluation
                of the CF avoids that path entirely.
                """
                arr = np.zeros((mesh.nv, 3))
                for vi, v in enumerate(mesh.vertices):
                    pt = mesh(*v.point)
                    val = cf(pt)
                    # cf(pt) returns a tuple/list/scalar; coerce
                    # to a length-3 vector.
                    if hasattr(val, "__len__"):
                        arr[vi, :len(val)] = [
                            float(getattr(x, "real", x))
                            for x in val[:3]]
                    else:
                        arr[vi, 0] = float(getattr(val, "real", val))
                return arr

            # Real-valued vector for the GMSH view: take the real
            # part of curl_A in the AC case (GMSH view files do not
            # carry phase — write a second .msh later for the
            # imaginary part if needed).
            if is_dc:
                B_cf = curl_A
            else:
                B_cf = ngsolve.CoefficientFunction(
                    tuple(curl_A[i].real for i in range(3)))

            node_B = _eval_vec_at_vertices(B_cf)
            _log(f"GMSH:B node values done ({mesh.nv} verts)")

            # J on the coil — restrict via IfPos on the material
            # indicator function (J_source is already zero outside
            # the coil region, so simple per-vertex evaluation works).
            node_J = _eval_vec_at_vertices(J_source)
            _log(f"GMSH:J node values done")

            post = GmshPostExport(mesh, boundary=False)
            post.add_field("B", node_B, ncomp=3)
            post.add_field("J", node_J, ncomp=3)
            post.write(msh_output)
            _log(f"GMSH:wrote {os.path.basename(msh_output)}")

            # Companion .msh.opt: GMSH auto-loads this when opening
            # the .msh. Replaces the old .geo companion approach.
            from gmsh_post_export import write_companion_opt
            opt_file = write_companion_opt(msh_output, n_vector_views=2)
            _log(f"GMSH:wrote {os.path.basename(opt_file)}")
            gmsh_file = msh_output
        except Exception as e:
            import traceback
            tb_text = traceback.format_exc()
            _log(f"GMSH_ERROR:{type(e).__name__}: {e}")
            # Surface the last 3 lines of the traceback so the panel
            # output window shows where the failure happened (the
            # one-line catch we had before was useless).
            for line in tb_text.splitlines()[-4:]:
                _log(f"GMSH_ERROR:  {line}")
            # Even on failure, try to write a minimal mesh-only .msh
            # so the user can at least open the geometry in GMSH.
            try:
                from gmsh_post_export import GmshPostExport
                post_min = GmshPostExport(mesh, boundary=False)
                post_min.write_mesh(msh_output)
                gmsh_file = msh_output
                _log(f"GMSH:fallback mesh-only -> {msh_output}")
            except Exception as e2:
                _log(f"GMSH_ERROR:fallback also failed: {e2}")

    # ============================================================
    # Step 8: Result JSON
    # ============================================================
    delta_skin = math.sqrt(2.0 / (omega * MU_0 * (mu_r if esim is None else 100) * sigma)) if (omega > 0 and sigma > 0) else 0

    result = {
        "P_total": float(P_total),
        "Q_total": float(Q_total),
        "L": float(L),
        "Z_s": str(Z_s),
        # H_t_rms is None when there is no workpiece (coil-only run);
        # serialize as null in JSON instead of crashing.
        "H_t_rms": float(H_t_rms) if H_t_rms is not None else None,
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
    parser.add_argument("--vol", required=True, help="Netgen .vol file")
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
                        choices=["steel", "copper", "aluminum", "custom"],
                        help="Workpiece material; 'custom' takes mu_r "
                             "and sigma from the explicit --mu-r / --sigma "
                             "arguments instead of the built-in BH table.")
    parser.add_argument("--current", type=float, default=1.0,
                        help="Coil current [A]")
    parser.add_argument("--half-thickness", type=float, default=0.005,
                        help="Workpiece half-thickness / characteristic "
                             "radius [m]. Used as R in the Dowell tanh "
                             "formula for the linear SIBC path AND as "
                             "the cell-problem domain length for the "
                             "ESIM nonlinear path.")
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
            vol_file=args.vol,
            fes_order=args.fes_order,
            frequency=args.frequency,
            sigma=args.sigma,
            impedance_model=args.impedance,
            mu_r=args.mu_r,
            bh_file=args.bh_file if args.bh_file else None,
            material=args.material,
            I_total=args.current,
            half_thickness=args.half_thickness,
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
