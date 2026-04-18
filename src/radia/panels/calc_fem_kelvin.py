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
  (T0 source/sink labels retired 2026-04-18; PEEC filament source only)

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
                          detect_kelvin_offset,
                          progress, calc_main,
                          EMMaterial, add_material_args)


def _log(msg):
    """Write progress to stderr (panel reads these)."""
    progress("FEM", msg)


def solve_fem(vol_file="", fes_order=1,
              frequency=7000, mat=None,
              impedance_model="sibc",
              formulation="total",
              I_total=1.0, half_thickness=0.005,
              max_iter=15, tol=1e-3, relax=0.5,
              solver="pardiso", reg=1e-6, shift_eps=1e-6,
              nthreads=0,
              msh_output="",
              peec_step="",
              peec_sigma=5.8e7,
              peec_nwinc=1,
              peec_nhinc=1):
    """3D FEM-SIBC solver for .vol mesh.

    Args:
        vol_file: Netgen .vol file (with material/boundary labels)
        fes_order: HCurl polynomial order (1-3)
        frequency: Operating frequency [Hz]
        mat: EMMaterial instance (workpiece properties: sigma, mu_r, BH)
        impedance_model: "sibc" (linear) or "esim" (nonlinear)
        formulation: "total" (standard) or "scattered" (two-step:
            solve free-space A_inc first, then A_scat with Robin RHS)
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
    if mat is None:
        mat = EMMaterial.from_name("steel")
    sigma = mat.sigma
    mu_r = mat.mu_r
    # NGSolve must be imported BEFORE Cubit
    import ngsolve  # noqa: F401
    from ngsolve import (H1, HCurl, HDiv, Periodic, BilinearForm, LinearForm,
                         GridFunction, Integrate, Conj, curl, grad, dx, ds, CF,
                         BND, VOL, TaskManager, Preconditioner, InnerProduct)
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

        # Periodic identification comes from C++ export (translation-based).
        # The .vol must contain identification pairs from ident.Add().
        n_ident = mesh.ngmesh.GetNrIdentifications()
        if n_ident == 0:
            return {"error":
                    "Mesh has 'kelvin' material but no periodic "
                    "identification in the .vol file. Ensure the .jou/.py "
                    "sets sideset 'kelvin_int' (interior sphere outer "
                    "boundary) and 'kelvin_ext' (exterior sphere boundary), "
                    "and uses 'copy mesh surface' for 1:1 correspondence. "
                    "The C++ exporter writes identification via translation "
                    "offset matching."}
        has_kelvin_periodic = True
        _log(f"PERIODIC:from_vol ({n_ident} ident(s), a={a_kelvin:.4f})")

    t_mesh = time.perf_counter() - t0

    materials = mesh.GetMaterials()
    boundaries = mesh.GetBoundaries()
    ne = mesh.GetNE(VOL)
    _log(f"MESH:done ({ne} elems, {t_mesh:.1f}s)")

    # ============================================================
    # Step 2: Source current -- PEEC filament only
    # ============================================================
    # T0 (Dirichlet-at-source/sink Laplace) was retired 2026-04-18 after
    # being shown unreliable on gapped geometries (1/r cusps at gap
    # corners). The remaining FEM source is the PEEC line-integral
    # total-field excitation; for a scattered A_r variant use
    # solve_fem_biot_savart().
    has_kelvin = "kelvin" in materials

    if "sibc" in boundaries:
        sibc_bnd = "sibc"
        has_wp = True
    else:
        sibc_bnd = ""
        has_wp = False

    is_hole = has_wp and "workpiece" not in materials
    if has_wp:
        _log(f"SIBC:boundary={sibc_bnd}, approach={'hole' if is_hole else 'interface'}")

    use_peec = bool(peec_step and os.path.isfile(peec_step))

    if use_peec:
        _log(f"SOURCE:PEEC filament from {peec_step}")
    else:
        return {"error":
                "--peec-step is required. The T0 source/sink technique "
                "was retired 2026-04-18 (unreliable on gapped geometries). "
                "Use --peec-step to pass a coil STEP file for the PEEC "
                "filament source."}

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
        esim = mat.create_esim_solver(frequency, half_thickness,
                                      geometry='cylinder')
        Z_s = esim.solve(5.0)['Z']
    else:
        # Linear SIBC via Dowell tanh formula
        Z_s = mat.dowell_Zs(frequency, half_thickness)
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
    use_complex = not is_dc

    # HCurl space for A (nograds=True always). T0-era A-V compound was
    # retired with T0; the PEEC filament source does not need an in-coil
    # scalar potential gauge.
    base_fes = HCurl(mesh, order=fes_order,
                     complex=use_complex,
                     nograds=True,
                     dirichlet=dirichlet_bnd)
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
    peec_A_s_cf = None
    peec_nu_cf_for_L = None
    f_lf = LinearForm(fes)
    if use_peec:
        # ---- PEEC filament source: total-field via line integrals ----
        # Strategy: assemble RHS f[j] = sum_k I_k * integral_{path_k} phi_j . dl
        # directly from filament paths and HCurl basis functions.
        # This is the total-field formulation: FEM solves for A_total,
        # not the scattered field A_r. Robin BC is on the LHS (A_total).
        #
        # No Biot-Savart vertex sampling, no HCurl projection, no
        # reduced-A cancellation error.
        from coil_from_cad import filaments_from_step
        from ngsolve import ElementId

        _log("PEEC:building topology")
        topo = filaments_from_step(peec_step, sigma=peec_sigma,
                                    nwinc=peec_nwinc, nhinc=peec_nhinc)

        if "filament_paths" in topo:
            fil_paths = topo["filament_paths"]
            n_fil = topo.get("n_loop", len(fil_paths))
            if n_fil == 1:
                fil_currents = [I_total]
            else:
                peec_solver = topo["solver"]
                I_branch = peec_solver.compute_branch_currents(
                    frequency, [I_total])
                seg_of_fil = topo["seg_of_filament"]
                fil_currents = []
                for segs in seg_of_fil:
                    avg_I = np.mean([I_branch[s] for s in segs])
                    fil_currents.append(complex(avg_I))
            _log(f"PEEC:{n_fil} filaments, "
                 f"{sum(len(p) for p in fil_paths)} total segments")
        elif "filament_path" in topo:
            fp = topo["filament_path"]
            fil_paths = [[(fp[i], fp[i + 1])
                          for i in range(len(fp) - 1)]]
            fil_currents = [I_total]
            _log(f"PEEC:legacy centerline {len(fp)-1} segments")
        else:
            return {"error": "filaments_from_step returned unknown format: "
                             + str(list(topo.keys()))}

        # Assemble line-integral RHS:
        #   f[j] = sum_k I_k * sum_seg integral phi_j(x) . dl dx
        # using Gauss-Legendre quadrature per segment.
        n_quad_line = 4
        t_gl, w_gl = np.polynomial.legendre.leggauss(n_quad_line)
        t01 = 0.5 * (t_gl + 1.0)
        w01 = 0.5 * w_gl

        _log("PEEC:assembling line-integral RHS")
        t0_li = time.perf_counter()
        f_vec = f_lf.vec.CreateVector()
        f_vec[:] = 0

        # Scratch GF for per-DOF basis function evaluation
        gf_scratch = GridFunction(fes)
        n_skip = 0

        for k, path in enumerate(fil_paths):
            Ik = complex(fil_currents[k])
            for (p1, p2) in path:
                p1a = np.asarray(p1, dtype=float)
                p2a = np.asarray(p2, dtype=float)
                dl = p2a - p1a
                seg_len = np.linalg.norm(dl)
                if seg_len < 1e-15:
                    continue
                for q in range(n_quad_line):
                    pt = p1a + t01[q] * dl
                    try:
                        mp = mesh(*pt)
                        elnr = mp.nr
                        if elnr < 0:
                            raise ValueError("outside mesh")
                        el = mesh[ElementId(VOL, elnr)]
                    except Exception:
                        n_skip += 1
                        continue
                    dofs = fes.GetDofNrs(el)
                    for d in dofs:
                        if d < 0:
                            continue  # constrained DOF
                        gf_scratch.vec[:] = 0
                        gf_scratch.vec[d] = 1.0
                        val = gf_scratch(mp)
                        dot = sum(complex(val[i]) * dl[i]
                                  for i in range(3))
                        f_vec[d] += Ik * w01[q] * dot

        t_li = time.perf_counter() - t0_li
        _log(f"PEEC:line-integral RHS assembled ({t_li:.1f}s, "
             f"{n_skip} skipped)")
        # Mark as total-field (no A_s to add later)
        peec_A_s_cf = None
        # Assemble empty LinearForm first, then overwrite with f_vec.
        f_lf.Assemble()
        f_lf.vec.data = f_vec

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

        # PEEC total-field: Robin BC is on the LHS (a_bf already has it).
        # No surface RHS from A_s needed.
        rhs_vec = f_lf.vec

        # Solver-specific setup
        if solver == "ams":
            # AMS was built around the T0 / J_source volume source; with
            # PEEC line-integral RHS (the only remaining path) it has
            # never been wired up. Fail loud rather than run the wrong
            # configuration silently.
            return {"error": "solver=ams is not supported with the PEEC "
                             "filament source. Use pardiso or bddc."}
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
            gfu.vec.data = iccg * rhs_vec
        elif solver == "bddc":
            pre = Preconditioner(a_bf, "bddc")
            a_bf.Assemble()
            # BVP reads f_lf.vec; for PEEC, rhs_vec includes the Robin
            # surface term from A_s. Copy into f_lf.vec if they differ.
            if rhs_vec is not f_lf.vec:
                f_lf.vec.data = rhs_vec
            with TaskManager():
                from ngsolve import solvers
                solvers.BVP(bf=a_bf, lf=f_lf, gf=gfu,
                            pre=pre, maxsteps=500, tol=1e-8)
        else:
            # pardiso (direct)
            if formulation == "scattered" and has_wp and abs(robin) > 0:
                # Scattered-field two-step solve:
                #   Step 1: A_inc = free-space solution (no Robin)
                #   Step 2: A_scat with Robin, RHS = -robin * A_inc on SIBC
                #   Total: A = A_inc + A_scat
                #
                # Volume RHS terms cancel exactly: A_inc satisfies
                #   int nu curl(A_inc) . curl(v) dx = int J . v dx
                # so the scattered-field RHS is surface-only.

                # Step 1: free-space (no Robin)
                a_free = BilinearForm(fes)
                a_free += nu_cf * curl(u) * curl(v) * dx(bonus_intorder=4)
                if reg > 0:
                    if has_kelvin and has_kelvin_periodic:
                        if non_kelvin_mats:
                            a_free += reg * NU_0 * u * v * dx(
                                "|".join(non_kelvin_mats))
                    else:
                        a_free += reg * NU_0 * u * v * dx
                a_free.Assemble()
                gfu_inc = GridFunction(fes)
                with TaskManager():
                    gfu_inc.vec.data = a_free.mat.Inverse(
                        fes.FreeDofs(), inverse="pardiso") * rhs_vec
                _log(f"SCATTERED:A_inc solved")

                # Step 2: A_scat with Robin
                # a_bf already has curl-curl + Robin + reg
                a_bf.Assemble()
                f_scat = LinearForm(fes)
                f_scat += -robin * gfu_inc * v.Trace() * ds(sibc_bnd)
                f_scat.Assemble()
                gfu_scat = GridFunction(fes)
                with TaskManager():
                    gfu_scat.vec.data = a_bf.mat.Inverse(
                        fes.FreeDofs(), inverse="pardiso") * f_scat.vec
                _log(f"SCATTERED:A_scat solved")

                # Total field
                gfu.vec.data = gfu_inc.vec + gfu_scat.vec
            else:
                # Standard total-field solve
                a_bf.Assemble()
                with TaskManager():
                    gfu.vec.data = a_bf.mat.Inverse(
                        fes.FreeDofs(), inverse="pardiso") * rhs_vec

        t_solve_iter = time.perf_counter() - t0_iter

        # H_t from SIBC relation: H_t = |jw/Z_s| * |A_t|
        #
        # Tangential projection: |A_t|^2 = |A|^2 - (A . n)^2
        # n = specialcf.normal(3) = outward from mesh domain (air)
        #
        # For PEEC (reduced-A): gfu is A_r, the total field is A_r + A_s.
        # H_t must use A_total, not just A_r.
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

            # A_eval: the field whose tangential trace drives H_t.
            if peec_A_s_cf is not None:
                A_eval = CF(tuple(gfu[i] + peec_A_s_cf[i]
                                  for i in range(3)))
            else:
                A_eval = gfu

            # |A|^2
            A_sq = sum(A_eval[i].real * A_eval[i].real +
                       A_eval[i].imag * A_eval[i].imag for i in range(3))
            # (A . n)^2 = (Re(A).n)^2 + (Im(A).n)^2
            Adn_re = sum(A_eval[i].real * n_bnd[i] for i in range(3))
            Adn_im = sum(A_eval[i].imag * n_bnd[i] for i in range(3))
            An_sq = Adn_re * Adn_re + Adn_im * Adn_im
            At_sq = A_sq - An_sq  # |A_t|^2

            int_At2 = Integrate(At_sq, mesh, BND,
                                definedon=wp_region).real
            At_rms = math.sqrt(max(int_At2, 0) / max(A_wp, 1e-30))
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
    if peec_A_s_cf is not None:
        B_total = curl(gfu) + curl(peec_A_s_cf)
        nu_for_L = peec_nu_cf_for_L if peec_nu_cf_for_L is not None else nu_cf
        W_mag = Integrate(
            0.5 * nu_for_L * B_total * Conj(B_total),
            mesh, order=10).real
    else:
        W_mag = Integrate(
            0.5 * nu_cf * curl(gfu) * Conj(curl(gfu)),
            mesh, order=10).real
    L = 2 * W_mag / I_total**2

    t_total = time.perf_counter() - t_total_start
    _log(f"DONE:P={P_total:.4e} L={L*1e9:.2f}nH t={t_total:.1f}s")

    # ============================================================
    # Step 7: GMSH export — B vector + J vector + companion .msh.opt
    # ============================================================
    # Match the BEM viz convention: VECTORS only (GMSH shows the
    # magnitude under the vector view automatically). The companion
    # .msh.opt hides volume mesh elements so arrows are visible.
    # Phase B: save .vol + .sol (NGSolve native) first, then convert to
    # GMSH via the shared vol2msh helper. Gives a single source of truth
    # (the .vol+.sol pair) so the .msh can be regenerated / swapped
    # later without re-solving.
    gmsh_file = ""
    vol_path = ""
    sol_paths = {}
    if msh_output:
        _log("SAVE:start .vol + .sol")
        try:
            from gmsh_post_export import save_vol_sol_pair, vol2msh
            base_dir = os.path.dirname(os.path.abspath(msh_output))
            name_stem = os.path.splitext(os.path.basename(msh_output))[0]

            # 1. Mesh + B field as a real-valued HDiv GridFunction.
            #    For the AC case we save the real component (phase is
            #    not carried in .sol; callers needing complex B should
            #    reconstruct from gfu.Load + curl on the loaded mesh).
            curl_A = curl(gfu)
            if is_dc:
                B_cf = curl_A
            else:
                B_cf = ngsolve.CoefficientFunction(
                    tuple(curl_A[i].real for i in range(3)))
            fes_B = HDiv(mesh, order=1)
            gf_B = GridFunction(fes_B)
            try:
                gf_B.Set(B_cf)
            except Exception:
                # Some AC problems resist Set on HDiv — fall back to
                # vertex-eval + H1 dim=3 GridFunction.
                fes_B = H1(mesh, order=1, dim=3)
                gf_B = GridFunction(fes_B)
                for vi, v in enumerate(mesh.vertices):
                    pt = mesh(*v.point)
                    val = B_cf(pt)
                    if hasattr(val, "__len__"):
                        for k in range(3):
                            gf_B.vec.FV()[vi * 3 + k] = float(
                                getattr(val[k], "real", val[k]))
            vol_B = os.path.join(base_dir, f"{name_stem}_fem.vol").replace("\\", "/")
            sol_B = os.path.join(base_dir, f"{name_stem}_B.sol").replace("\\", "/")
            save_vol_sol_pair(vol_B, sol_B, mesh.ngmesh, gf_B)
            vol_path = vol_B
            sol_paths["B"] = sol_B

            # 2. Source field as a separate GridFunction.
            sol_entries = [
                {"sol": sol_B, "fes": type(fes_B).__name__,
                 "fes_order": 1,
                 "fes_dim": getattr(fes_B, "dim", 1),
                 "name": "B", "ncomp": 3},
            ]
            if peec_A_s_cf is not None:
                # PEEC source: export A_s as H1(dim=3) for visualization
                sol_As = os.path.join(base_dir,
                                      f"{name_stem}_As.sol").replace("\\", "/")
                peec_A_s_cf.Save(sol_As)
                sol_paths["A_s"] = sol_As
                sol_entries.append(
                    {"sol": sol_As, "fes": "H1", "fes_order": 1,
                     "fes_dim": 3, "name": "A_s", "ncomp": 3})
            _log(f"SAVE:{os.path.basename(vol_B)} + {len(sol_paths)} .sol")

            # 3. Convert .vol + .sol to .msh via the shared helper.
            vol2msh(msh_output, vol_B, sol_entries)
            gmsh_file = msh_output
            _log(f"GMSH:wrote {os.path.basename(msh_output)}")
        except Exception as e:
            import traceback
            tb_text = traceback.format_exc()
            _log(f"GMSH_ERROR:{type(e).__name__}: {e}")
            for line in tb_text.splitlines()[-4:]:
                _log(f"GMSH_ERROR:  {line}")

    # ============================================================
    # Step 8: Result JSON
    # ============================================================
    delta_skin = mat.skin_depth(frequency)

    # PEEC: L from port impedance (mesh-independent, validated)
    L_peec = None
    if use_peec and 'topo' in dir():
        peec_solver_obj = topo.get("solver")
        if peec_solver_obj is not None:
            Z_port = peec_solver_obj.compute_port_impedance(frequency)
            L_peec = float(Z_port.imag) / (2 * math.pi * frequency)
            _log(f"PEEC:L_peec={L_peec*1e9:.3f}nH (Z_port)")

    result = {
        "P_total": float(P_total),
        "Q_total": float(Q_total),
        "L": float(L),
        "L_peec": L_peec,
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
        "formulation": formulation,
        "source_type": "PEEC",
        "has_kelvin": has_kelvin,
        "msh_file": gmsh_file,
    }
    return result


def solve_fem_biot_savart(vol_file="", fes_order=1,
                           frequency=7000, mat=None,
                           I_total=1.0, half_thickness=0.005,
                           peec_step="",
                           peec_sigma=5.8e7,
                           peec_nwinc=3, peec_nhinc=3,
                           solver="pardiso", nthreads=0):
    """Scattered A_r formulation with Biot-Savart A_s excitation (no T0).

    Solves for the scattered vector potential A_r on the given .vol mesh,
    using the PEEC filament Biot-Savart A_s as the excitation (entering
    through the SIBC Robin boundary term). The coil volume is treated
    as air from the FEM side; PEEC captures the coil skin / proximity
    effect through the nwinc/nhinc filament subdivision.

    L is returned as ``L_peec + Delta_L`` where L_peec comes from the
    PEEC circuit port impedance and Delta_L is the line integral of
    A_r along the filament paths. This decomposition is numerically
    well-behaved because neither term relies on a source/sink T0
    Dirichlet lift, so there are no gap-corner 1/r cusps.

    Returns a dict with keys matching solve_fem: P_total, L, Z_s, ...
    """
    if mat is None:
        mat = EMMaterial.from_name("steel")
    if not peec_step or not os.path.exists(peec_step):
        return {"error": f"--peec-step is required for Biot-Savart mode "
                         f"(got: {peec_step!r})"}

    import ngsolve  # noqa: F401
    from ngsolve import (Mesh, HCurl, Periodic, BilinearForm, LinearForm,
                         GridFunction, Integrate, curl, dx, ds, CF,
                         BND, VOL, TaskManager)
    from ngsolve import x, y, z
    from coil_from_cad import filaments_from_step
    from kelvin_source import (project_A_s_to_hcurl,
                                line_integral_A_filaments,
                                compute_back_reaction)

    setup_paths()
    t_total_start = time.perf_counter()

    # --- 1. PEEC topology + port impedance + per-filament currents ---
    _log("MESH:loading")
    t0 = time.perf_counter()
    mesh = Mesh(vol_file)
    t_mesh = time.perf_counter() - t0

    materials = mesh.GetMaterials()
    boundaries = set(mesh.GetBoundaries())

    has_kelvin = 'kelvin' in materials
    if 'sibc' not in boundaries:
        return {"error": "Mesh has no 'sibc' boundary — SIBC workpiece required"}

    if has_kelvin:
        kelvin_center = np.array(detect_kelvin_offset(mesh))
        kelvin_verts = set()
        for el in mesh.Elements(VOL):
            if el.mat == 'kelvin':
                for v in el.vertices:
                    kelvin_verts.add(v.nr)
        coords = np.array([mesh.vertices[v].point for v in kelvin_verts])
        dists = np.linalg.norm(coords - kelvin_center[None, :], axis=1)
        a_kelvin = float(np.max(dists))
        n_ident = mesh.ngmesh.GetNrIdentifications()
        has_kelvin_periodic = n_ident > 0
    else:
        kelvin_center = np.array([0.0, 0.0, 0.0])
        a_kelvin = 0.0
        has_kelvin_periodic = False
    _log(f"MESH:ne={mesh.GetNE(VOL)} kelvin={has_kelvin} "
         f"periodic={has_kelvin_periodic} t={t_mesh:.1f}s")

    # --- 2. PEEC filaments + currents ---
    t0 = time.perf_counter()
    topo = filaments_from_step(peec_step, sigma=peec_sigma,
                                nwinc=peec_nwinc, nhinc=peec_nhinc)
    peec_solver = topo['solver']
    Z_port = peec_solver.compute_port_impedance(frequency)
    L_peec = float(Z_port.imag) / (2 * math.pi * frequency)
    R_peec = float(Z_port.real)

    fil_paths = topo['filament_paths']
    seg_of_fil = topo['seg_of_filament']
    I_branch = peec_solver.compute_branch_currents(frequency, [I_total])
    fil_currents = np.array([
        complex(np.mean([I_branch[s] for s in segs]))
        for segs in seg_of_fil
    ])
    _log(f"PEEC:n_fil={len(fil_paths)} L_peec={L_peec*1e9:.3f}nH "
         f"t={time.perf_counter()-t0:.1f}s")

    # --- 3. nu_cf: Kelvin weight, coil treated as air ---
    kx, ky, kz = kelvin_center
    nu_dict = {}
    for m in materials:
        if 'kelvin' in m.lower():
            dx_k = x - kx
            dy_k = y - ky
            dz_k = z - kz
            rp_sq = dx_k * dx_k + dy_k * dy_k + dz_k * dz_k + 1e-20
            nu_dict[m] = NU_0 * rp_sq / a_kelvin ** 2
        else:
            nu_dict[m] = NU_0
    nu_cf = mesh.MaterialCF(nu_dict, default=NU_0)

    # --- 4. Biot-Savart A_s (Kelvin pullback in exterior) ---
    t0 = time.perf_counter()
    gf_As = project_A_s_to_hcurl(
        mesh, fil_paths, fil_currents,
        kelvin_center=kelvin_center if has_kelvin else None,
        R_kelvin=a_kelvin if has_kelvin else None,
        factor_mode='pullback', is_complex=True, order=fes_order)
    _log(f"A_s:proj t={time.perf_counter()-t0:.1f}s")

    # --- 5. FES + SIBC Robin system for scattered A_r ---
    Z_s = mat.dowell_Zs(frequency, half_thickness)
    omega = 2 * math.pi * frequency
    robin = 1j * omega / Z_s
    _log(f"SIBC:Z_s={Z_s:.4e}, Robin={robin:.4e}")

    dirichlet_bnd = 'GND' if 'GND' in boundaries else ''
    base_fes = HCurl(mesh, order=fes_order, nograds=True, complex=True,
                     dirichlet=dirichlet_bnd)
    fes = Periodic(base_fes) if has_kelvin_periodic else base_fes

    u, v_ = fes.TnT()
    a = BilinearForm(fes)
    a += nu_cf * curl(u) * curl(v_) * dx(bonus_intorder=4)

    non_kelvin_mats = [m for m in materials if 'kelvin' not in m.lower()]
    if non_kelvin_mats:
        a += 1e-6 * NU_0 * u * v_ * dx('|'.join(non_kelvin_mats))

    a += robin * u.Trace() * v_.Trace() * ds('sibc')

    f_lf = LinearForm(fes)
    f_lf += -robin * gf_As * v_.Trace() * ds('sibc')

    t0 = time.perf_counter()
    a.Assemble()
    f_lf.Assemble()
    t_asm = time.perf_counter() - t0

    gfu = GridFunction(fes)
    t0 = time.perf_counter()
    if nthreads > 0:
        ngsolve.SetNumThreads(nthreads)
    with TaskManager():
        gfu.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse=solver) * f_lf.vec
    t_solve = time.perf_counter() - t0
    _log(f"FEM:ndof={fes.ndof} asm={t_asm:.1f}s solve={t_solve:.1f}s")

    # --- 6. Back-reaction Delta_L, Delta_R from line integral ---
    delta_phi = line_integral_A_filaments(gfu, mesh, fil_paths, n_quad=4)
    back = compute_back_reaction(fil_currents, delta_phi, I_total)
    Delta_L = back['Delta_L']
    Delta_R = omega * back['Delta_R_over_omega']
    P_wp = 0.5 * Delta_R * abs(I_total) ** 2

    # --- 6b. Power + H_t_rms from SIBC surface integral ---
    # (fem_esim_3d convention -- validated against 2D axisym to <1%)
    # A_total on the wp surface = gf_As + gfu (scattered).
    # Tangential A: A_t = A - (A.n) n. H_t_rms uses |A_t|^2 integrated
    # over the sibc surface divided by wp area.
    from ngsolve import specialcf
    n_bnd = specialcf.normal(3)
    A_total = CF(tuple(gf_As[i] + gfu[i] for i in range(3)))
    A_sq = sum(A_total[i].real * A_total[i].real
               + A_total[i].imag * A_total[i].imag for i in range(3))
    Adn_re = sum(A_total[i].real * n_bnd[i] for i in range(3))
    Adn_im = sum(A_total[i].imag * n_bnd[i] for i in range(3))
    An_sq = Adn_re * Adn_re + Adn_im * Adn_im
    At_sq = A_sq - An_sq  # |A_t|^2

    wp_region = mesh.Boundaries('sibc')
    A_wp = Integrate(CF(1), mesh, BND, definedon=wp_region).real
    At_int = Integrate(At_sq, mesh, BND, definedon=wp_region).real
    At_rms = math.sqrt(max(At_int, 0) / max(A_wp, 1e-30))
    H_t_rms = abs(1j * omega / Z_s) * At_rms
    P_surf = 0.5 * Z_s.real * H_t_rms**2 * A_wp
    Q_surf = 0.5 * Z_s.imag * H_t_rms**2 * A_wp

    # --- 6c. Skin-layer stored energy term (per sibc_skin_energy_fix) ---
    # L_skin = omega * Im(Z_s)/|Z_s|^2 * int |A_t|^2 dS / |I|^2
    # Captures the magnetic energy inside the workpiece skin layer that
    # SIBC hides from the FEM mesh. For Cu (mu_r=1) this is small; for
    # steel (mu_r=100) it dominates.
    L_skin = omega * Z_s.imag / (abs(Z_s) ** 2) * At_int / abs(I_total) ** 2

    L_total = L_peec + Delta_L
    # Total L = PEEC circuit L + flux-linkage back-reaction + skin-layer stored energy.
    L_total_with_skin = L_total + L_skin

    _log(f"DONE:L={L_total_with_skin*1e9:.2f}nH "
         f"(L_peec={L_peec*1e9:.2f} + DL={Delta_L*1e9:+.2f} "
         f"+ L_skin={L_skin*1e9:+.2f}) "
         f"P_line={P_wp:.4e} P_surf={P_surf:.4e} H_t={H_t_rms:.2f} "
         f"t={time.perf_counter()-t_total_start:.1f}s")

    return {
        "P_total": P_surf,
        "P_line": P_wp,
        "Q_total": Q_surf,
        "L": L_total_with_skin,
        "L_no_skin": L_total,
        "L_peec": L_peec,
        "R_peec": R_peec,
        "Delta_L": Delta_L,
        "Delta_R": Delta_R,
        "L_skin": L_skin,
        "Z_s_real": Z_s.real,
        "Z_s_imag": Z_s.imag,
        "H_t_rms": H_t_rms,
        "A_wp": A_wp,
        "iterations": 0,
        "converged": True,
        "t_mesh": round(t_mesh, 2),
        "t_total": round(time.perf_counter() - t_total_start, 2),
        "frequency": frequency,
        "sigma": mat.sigma,
        "impedance_model": "sibc",
        "formulation": "scattered-biot-savart",
        "source_type": "biot-savart",
        "has_kelvin": has_kelvin,
        "msh_file": "",
    }


def main():
    parser = argparse.ArgumentParser(
        description="3D FEM-SIBC with PEEC filament source + optional Kelvin")
    parser.add_argument("--vol", required=True, help="Netgen .vol file")
    parser.add_argument("--fes-order", type=int, default=1,
                        help="HCurl polynomial order (1-3)")
    parser.add_argument("--frequency", type=float, default=7000,
                        help="Frequency [Hz]")
    add_material_args(parser, include_custom=True)
    parser.add_argument("--impedance", default="sibc",
                        choices=["sibc", "esim"],
                        help="Impedance model")
    parser.add_argument("--formulation", default="total",
                        choices=["total", "scattered"],
                        help="total: standard single solve; "
                             "scattered: two-step (A_inc free-space "
                             "+ A_scat with Robin RHS)")
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
    parser.add_argument("--peec-step", default="",
                        help="STEP file for PEEC coil (filament source)")
    parser.add_argument("--peec-sigma", type=float, default=5.8e7,
                        help="PEEC coil conductivity [S/m]")
    parser.add_argument("--peec-nwinc", type=int, default=1,
                        help="PEEC width subdivision (default 1)")
    parser.add_argument("--peec-nhinc", type=int, default=1,
                        help="PEEC height subdivision (default 1)")
    parser.add_argument("--source-mode", default="scattered",
                        choices=["scattered", "total"],
                        help="scattered (default): Biot-Savart A_s excitation "
                             "with FEM scattered A_r (fast, PEEC captures "
                             "coil skin via nwinc). "
                             "total: FEM line-integral PEEC RHS (total-field, "
                             "supports ESIM Karl iteration and GMSH viz).")

    def run(args):
        if not args.peec_step:
            return {"error": "--peec-step is required. T0/AV paths were "
                             "retired 2026-04-18; all FEM source injection "
                             "goes through PEEC filaments."}
        if args.source_mode == "scattered":
            return solve_fem_biot_savart(
                vol_file=args.vol,
                fes_order=args.fes_order,
                frequency=args.frequency,
                mat=EMMaterial.from_args(args),
                I_total=args.current,
                half_thickness=args.half_thickness,
                peec_step=args.peec_step,
                peec_sigma=args.peec_sigma,
                peec_nwinc=args.peec_nwinc,
                peec_nhinc=args.peec_nhinc,
                solver=args.solver,
                nthreads=args.nthreads,
            )
        return solve_fem(
            vol_file=args.vol,
            fes_order=args.fes_order,
            frequency=args.frequency,
            mat=EMMaterial.from_args(args),
            impedance_model=args.impedance,
            formulation=args.formulation,
            I_total=args.current,
            half_thickness=args.half_thickness,
            max_iter=args.max_iter,
            solver=args.solver,
            reg=args.reg,
            shift_eps=args.shift_eps,
            nthreads=args.nthreads,
            msh_output=args.msh_output,
            peec_step=args.peec_step,
            peec_sigma=args.peec_sigma,
            peec_nwinc=args.peec_nwinc,
            peec_nhinc=args.peec_nhinc,
        )

    calc_main(run, parser)


if __name__ == "__main__":
    main()
