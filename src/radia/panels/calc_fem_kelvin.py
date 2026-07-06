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
              peec_nhinc=1,
              peec_n_peri=0,
              esim_per_panel=False):
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
    # Resolve solver=auto: ams for fes_order=1 (Compact AMS+COCR is the
    # cheapest preconditioner for HCurl order-1 AC), bddc for higher
    # orders (CompactAMS auxiliary space requires order-1).
    if solver == "auto":
        solver = "ams" if fes_order == 1 else "bddc"
        _log(f"SOLVER:auto -> {solver} (fes_order={fes_order})")

    if fes_order > curve_order:
        return {
            "error": f"fes_order ({fes_order}) > mesh curve_order ({curve_order}). "
                     f"Re-export .vol with 'export netgen ... order {fes_order}' "
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
    # total-field excitation. The alternative Biot-Savart scattered
    # formulation (solve_fem_biot_savart) was retired 2026-04-24 because
    # of an unresolved ~3.4x P_wp under-prediction.
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
        # Seed Z_s with a small-H Picard solve; cap inner iter at 5
        # because the outer Karl loop will refresh Z_s immediately.
        esim = mat.create_esim_solver(frequency, half_thickness,
                                      geometry='cylinder')
        Z_s = esim.solve(5.0, max_iter=5)['Z']
    else:
        # Linear SIBC via Dowell tanh formula
        Z_s = mat.dowell_Zs(frequency, half_thickness)
        esim = None

    _log(f"SIBC:Z_s={Z_s:.4e}, model={impedance_model}")

    # ============================================================
    # Step 5: FE space + Karl iteration
    # ============================================================
    # Far-field truncation for HCurl A.
    #
    # NB (2026-04-28): an earlier version also accepted "GND" here as
    # a Dirichlet target, on the assumption that pinning A at a single
    # vertex would lock the gauge.  Empirically this is a no-op:
    # HCurl DOFs live on EDGES, not vertices, so HCurl(dirichlet="GND")
    # leaves every DOF free (verified: free DOFs unchanged before /
    # after on the IH coarse sample).  Vertex Dirichlet works in the
    # A-phi compound only on the H1 phi side.  Removed here to avoid
    # a misleading "Dirichlet on GND vertex" log line and a panel
    # warning that promised non-existent gauge fixing.
    #
    # Three real cases for the HCurl A solve:
    #   (A) kelvin material + Periodic identifications
    #          -> Periodic Kelvin closes the far field.  Gauge null-
    #             space is locked by nograds=True + reg*nu0*u*v
    #             (+ j*omega*sigma*A in AC), NOT by Dirichlet.
    #   (B) "outer" boundary face, no kelvin
    #          -> Dirichlet A=0 on outer face (the only Dirichlet
    #             that actually constrains HCurl edge DOFs).
    #   (C) neither
    #          -> reg-only gauge fix.  Result depends on reg unless
    #             an SIBC face absorbs the radiating field.
    dirichlet_bnd = ""
    if not has_kelvin and "outer" in boundaries:
        dirichlet_bnd = "outer"

    if has_kelvin and has_kelvin_periodic:
        _log("FARFIELD:Periodic Kelvin (open boundary; gauge by "
             "nograds + reg + j*omega*sigma)")
    elif dirichlet_bnd == "outer":
        _log("FARFIELD:Dirichlet A=0 on 'outer' boundary "
             "(truncation error if outer is too close to coil)")
    else:
        _log(f"FARFIELD:gauge regularisation only (reg={reg:.1e}). "
             "No 'kelvin' material AND no 'outer' face -- result is "
             "gauge-dependent unless an SIBC face absorbs the "
             "radiating field.")

    is_dc = (omega <= 0 or sigma <= 0)
    use_complex = not is_dc

    # HCurl space for A (nograds=True always). T0-era A-V compound was
    # retired with T0; the PEEC filament source does not need an in-coil
    # scalar potential gauge.
    #
    # Design split (2026-04-20, confirmed by user):
    #   - calc_fem_kelvin (this file):    HCurl A-only + PEEC filament
    #     source (total-field RHS via line integrals). Production-fast
    #     path, the coil is NOT meshed as material.
    #   - calc_fem_coilmesh (sibling):    compound HCurl(A) x H1(phi)
    #     with voltage-driven conductor source (--coil-sigma). "FEM truth"
    #     for PEEC validation. Requires a `coil` material + source/sink
    #     boundaries in the .vol.
    #
    # FEM reduced-A path (Biot-Savart B_s on integration points, HCurl
    # projection) is intentionally NOT provided: it cancellation-errors
    # against the coil-excluded A_r solution, and the PEEC line-integral
    # RHS covers the same use case without the projection error. See
    # project_fem_vs_peec_skin_2026_04_18.md for the validation
    # comparison (A-V vs PEEC+BEM-SIBC: H_t +3.0%, P +4.4%).
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
        from radia.coil_from_cad import filaments_from_step
        from ngsolve import ElementId

        _log("PEEC:building topology")
        if peec_n_peri > 0:
            # Perimeter-only placement (thin-skin, dispatch path 2b
            # with circle-edge centers).  Use this for 3turn helix
            # / racetrack STEPs where the walker fails on seed plane.
            topo = filaments_from_step(peec_step, sigma=peec_sigma,
                                        n_peri=peec_n_peri,
                                        use_coil_builder=True)
        else:
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
        # Total quadrature points attempted: sum over filaments and
        # segments of n_quad_line.  If 100% of them were skipped
        # (every filament point lies outside the volume mesh), the
        # RHS will be exactly zero and the downstream solve will
        # silently return gfu=0 → all P_total, L, q_surf evaluate to
        # zero.  This is a frequent setup error: the workpiece .vol
        # mesh doesn't include the coil region.
        # Detect and abort BEFORE solving so the user sees a clear,
        # actionable error instead of "AMS broken".
        n_seg_total = sum(len(p) for p in fil_paths)
        n_quad_total = n_seg_total * n_quad_line
        if n_quad_total > 0 and n_skip == n_quad_total:
            return {
                "error": (
                    f"PEEC line-integral RHS is zero: every one of "
                    f"the {n_skip}/{n_quad_total} filament "
                    f"quadrature points falls OUTSIDE the workpiece "
                    f".vol mesh.  The coil STEP "
                    f"({os.path.basename(peec_step)}) must lie "
                    f"inside the air region of the workpiece .vol "
                    f"({os.path.basename(vol_file)}).  Either: (a) extend "
                    f"the .vol's air region (Cubit .jou) to enclose "
                    f"the coil, OR (b) use a different method "
                    f"(PEEC inductance / weak-coupled BEM) that "
                    f"does NOT require the coil to be inside the "
                    f"workpiece mesh.")
            }
        if n_quad_total > 0 and n_skip > 0:
            pct = 100.0 * n_skip / n_quad_total
            _log(f"PEEC:WARN {pct:.1f}% of quadrature points are "
                 f"outside the workpiece mesh ({n_skip}/"
                 f"{n_quad_total}).  The coil may extend beyond the "
                 f"air region; check .vol geometry.")
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

    # Per-panel ESIM (v4.48.0+): build a surface H1 GridFunction on the
    # workpiece BND so the Robin coefficient becomes a per-DOF CF instead
    # of a scalar.  Iteration math is per-element variational projection;
    # see docs/esim/MATHEMATICAL_ANALYSIS.md § 3.3-3.4 for the BEM
    # equivalent and the rationale.
    if esim_per_panel and esim is not None and has_wp:
        # Boundary H1 space at the SAME order as the FES so the trace
        # integral has compatible dofs.  L2 boundary projection of Z_s
        # is sound for piecewise smoothly-varying surface impedances
        # (the typical case for IH workpieces).
        fes_Zs = H1(mesh, order=fes_order, complex=True)
        gf_Zs = GridFunction(fes_Zs)
        gf_Zs.vec[:] = 0
        gf_Zs.Set(CF(Z_s), definedon=mesh.Boundaries(sibc_bnd))
        # Mask: which DOFs got set (i.e. live on the wp BND surface).
        # Mirror the q_surf precedent for picking out BND dofs.
        mask_gf = GridFunction(H1(mesh, order=fes_order))
        mask_gf.vec[:] = 0
        mask_gf.Set(CF(1.0), definedon=mesh.Boundaries(sibc_bnd))
        bnd_mask = np.asarray(mask_gf.vec.FV().NumPy()) > 0.5
        n_bnd_dofs = int(bnd_mask.sum())
        _log(f"PER_PANEL:wp surface H1 dofs={n_bnd_dofs}")
    else:
        fes_Zs = None
        gf_Zs = None
        bnd_mask = None
        n_bnd_dofs = 0

    for iteration in range(max_iter):
        t0_iter = time.perf_counter()

        # Robin coefficient:
        #   Hole approach: +jw/Z_s (SIBC on external boundary of air domain)
        #   Interface: sign depends on which side NGSolve evaluates from
        #              +jw/Z_s works for both (tested with tangential projection)
        if gf_Zs is not None:
            # Per-DOF Z_s; CF arithmetic broadcasts naturally.
            robin = 1j * omega / gf_Zs
        else:
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

        # robin is either a complex scalar (legacy) or a CF (per-DOF).
        # Skip the surface term only when has_wp is False or scalar robin
        # is identically zero; CF case always has the term.
        if has_wp:
            if gf_Zs is not None:
                a_bf += robin * u.Trace() * v.Trace() * ds(sibc_bnd)
            elif abs(robin) > 0:
                a_bf += robin * u.Trace() * v.Trace() * ds(sibc_bnd)

        # PEEC total-field: Robin BC is on the LHS (a_bf already has it).
        # No surface RHS from A_s needed.
        rhs_vec = f_lf.vec

        # Solver-specific setup
        if solver == "ams":
            # Compact AMS preconditioner + COCR solver (AC eddy-current).
            # Reference: src/ext/sparsesolv/examples/hiruma/bench_compact_ams.py
            #
            # Re-wired 2026-05-08 for the PEEC line-integral RHS path:
            # the bilinear form a_bf above already has curl-curl + reg +
            # SIBC Robin; the PEEC RHS is in f_lf.vec.  We only need to
            # build the REAL auxiliary bilinear form a_real (CompactAMS
            # internal preconditioner space) + grab vertex coords.
            #
            # Constraints carried over from the pre-fec3465f impl:
            #   - fes_order = 1 (CompactAMS auxiliary space is order-1 HCurl)
            #   - nograds = True (already set on base_fes; H1 aux space
            #                     matches vertex count)
            #   - AC mode (omega > 0)
            #   - No Periodic Kelvin (DOF dimensions don't match)
            if fes_order > 1:
                return {"error": f"solver=ams requires fes_order=1 "
                                 f"(got {fes_order}); use solver=bddc for p>=2"}
            if not (omega > 0):
                return {"error": "solver=ams requires AC (frequency>0). "
                                 "For DC use solver=bddc or pardiso."}
            if has_kelvin and has_kelvin_periodic:
                return {"error": "solver=ams is not supported with "
                                 "Periodic Kelvin BC. Use bddc."}

            import radia.sparsesolv_ngsolve as ssn

            # Real auxiliary fes for the AMS preconditioner.  Same edge
            # topology as the complex `fes` (HCurl order=1 nograds=True
            # dirichlet=dirichlet_bnd).
            fes_real = HCurl(mesh, order=1, complex=False, nograds=True,
                              dirichlet=dirichlet_bnd)
            ur, vr = fes_real.TnT()
            a_real = BilinearForm(fes_real)
            a_real += nu_cf * curl(ur) * curl(vr) * dx(bonus_intorder=4)
            if reg > 0:
                if has_kelvin and has_kelvin_periodic:
                    if non_kelvin_mats:
                        a_real += shift_eps * NU_0 * ur * vr * dx(
                            "|".join(non_kelvin_mats))
                else:
                    a_real += shift_eps * NU_0 * ur * vr * dx
            if has_wp and abs(robin) > 0:
                a_real += abs(robin) * ur.Trace() * vr.Trace() * ds(sibc_bnd)
            # AIR-REGION REGULARIZATION (preconditioner only — NOT a_bf):
            # CompactAMS auxiliary form expects a substantial volumetric
            # mass `|omega| * sigma * u*v*dx("cond")` that dominates
            # curl-curl (cf. bench_compact_ams.py:80-83).  With SIBC the
            # workpiece is a HOLE — no volumetric sigma anywhere — so
            # the auxiliary form a_real has only curl-curl + tiny
            # shift_eps mass + boundary Robin trace.  The CompactAMS
            # internal factorization then segfaults at COCR matvec.
            # Fix: add an air-region "synthetic AC mass" |omega| * NU_0
            # to a_real ONLY (a_bf is unchanged → physics preserved).
            # Scale matches the bench's conductor mass within ~2 orders
            # of magnitude and gives the AMS algebra something to grip.
            a_real += abs(omega) * NU_0 * ur * vr * dx
            a_real.Assemble()

            grad_mat, fes_h1 = fes_real.CreateGradient()
            nv = mesh.nv
            coord_x = [0.0] * nv
            coord_y = [0.0] * nv
            coord_z = [0.0] * nv
            for i in range(nv):
                pt = mesh.vertices[i].point
                coord_x[i] = float(pt[0])
                coord_y[i] = float(pt[1])
                coord_z[i] = float(pt[2])

            t0_pre = time.perf_counter()
            pre_ams = ssn.ComplexCompactAMSPreconditioner(
                a_real_mat=a_real.mat, grad_mat=grad_mat,
                freedofs=fes_real.FreeDofs(),
                coord_x=coord_x, coord_y=coord_y, coord_z=coord_z,
                ndof_complex=fes.ndof, cycle_type=1, print_level=0)
            _log(f"AMS:preconditioner setup {time.perf_counter()-t0_pre:.1f}s")

            a_bf.Assemble()
            with TaskManager():
                cocr = ssn.COCRSolver(a_bf.mat, pre_ams,
                                       freedofs=fes.FreeDofs(),
                                       maxiter=500, tol=1e-8,
                                       printrates=False)
                gfu.vec.data = cocr * rhs_vec
            _log(f"AMS:COCR iters={cocr.iterations}")

            # Sanity check: COCR can spuriously "converge" in 0
            # iterations when the synthetic-air-mass preconditioner
            # (line 582 above) is so dominant that pre_ams * b is
            # below the tolerance threshold.  The result is gfu == 0
            # everywhere — physics broken, downstream P_total / L /
            # q_surf all evaluate to 0.  Detect this case and refuse
            # rather than emit silently-wrong numbers.
            # (kubota mdx report 2026-05-12: 3turnCoil_work, no
            # Kelvin, SIBC workpiece; iters=0, all stats zero.)
            import numpy as _np
            sol_norm = float(_np.linalg.norm(
                _np.asarray(gfu.vec.FV().NumPy())))
            rhs_norm = float(_np.linalg.norm(
                _np.asarray(rhs_vec.FV().NumPy())))
            if cocr.iterations == 0 or (
                    rhs_norm > 1e-30 and sol_norm / rhs_norm < 1e-30):
                return {
                    "error": (
                        f"solver=ams produced a null solution "
                        f"(iters={cocr.iterations}, "
                        f"|gfu|={sol_norm:.3e}, |rhs|={rhs_norm:.3e}).  "
                        f"The Compact AMS preconditioner is not "
                        f"well-conditioned for SIBC-only problems "
                        f"(workpiece is a hole, no volumetric "
                        f"conductor mass).  Re-run with "
                        f"--solver bddc (iterative, supports any "
                        f"order) or --solver pardiso (direct, "
                        f"memory-heavy).  Tracked: AMS works for "
                        f"problems with a volumetric conductor "
                        f"(non-SIBC).")
                }
        elif solver == "iccg":
            # Shifted ICCG (IC preconditioned CG, single-thread efficient)
            a_bf.Assemble()

            import radia.sparsesolv_ngsolve as ssn
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
            if gf_Zs is None and not (abs(Z_s) > 0):
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
            if gf_Zs is not None:
                # Per-DOF mode: report area-averaged |Z_s| for the
                # scalar H_t_rms diagnostic.
                Z_s_arr = np.asarray(gf_Zs.vec.FV().NumPy())
                Z_s_avg_abs = float(np.mean(np.abs(Z_s_arr[bnd_mask])))
                H_t_rms = abs(1j * omega / max(Z_s_avg_abs, 1e-30)) * At_rms
            else:
                H_t_rms = abs(1j * omega / Z_s) * At_rms

        # Update Z_s
        if gf_Zs is not None and esim is not None and has_wp:
            # Per-DOF Karl: project |A_t|^2 onto a real surface-H1 field,
            # extract per-DOF values, run ESIM per BND DOF, blend with
            # under-relaxation, write back into gf_Zs.
            fes_real = H1(mesh, order=fes_order)
            gf_At_re = GridFunction(fes_real)
            gf_At_im = GridFunction(fes_real)
            gf_At_re.vec[:] = 0
            gf_At_im.vec[:] = 0
            A_sq_re = sum(A_eval[i].real * A_eval[i].real
                          for i in range(3))
            A_sq_im = sum(A_eval[i].imag * A_eval[i].imag
                          for i in range(3))
            Adn_re_sq = Adn_re * Adn_re
            Adn_im_sq = Adn_im * Adn_im
            At_sq_re = A_sq_re - Adn_re_sq  # |Re(A_t)|^2
            At_sq_im = A_sq_im - Adn_im_sq  # |Im(A_t)|^2
            gf_At_re.Set(At_sq_re, definedon=wp_region)
            gf_At_im.Set(At_sq_im, definedon=wp_region)
            At_re_arr = np.asarray(gf_At_re.vec.FV().NumPy())
            At_im_arr = np.asarray(gf_At_im.vec.FV().NumPy())
            # Per-DOF |A_t| amplitude (matches the scalar mesh-RMS
            # convention: amplitude of phasor, not time-RMS).
            At_amp_per_dof = np.sqrt(
                np.maximum(At_re_arr, 0.0) + np.maximum(At_im_arr, 0.0))
            # Per-DOF H_t = |jw / Z_s| * |A_t|
            Z_s_old_arr = np.asarray(gf_Zs.vec.FV().NumPy()).copy()
            Z_s_new_arr = Z_s_old_arr.copy()
            n_called = 0
            for i in np.flatnonzero(bnd_mask):
                Zsi_old = Z_s_old_arr[i]
                if abs(Zsi_old) <= 1e-30:
                    continue
                Ht_i = abs(1j * omega / Zsi_old) * float(At_amp_per_dof[i])
                sol_new = esim.solve(max(Ht_i, 1e-3))
                Z_s_new_arr[i] = complex(sol_new['Z'])
                n_called += 1
            # Under-relax per-DOF.
            Z_s_blend = (relax * Z_s_new_arr
                         + (1 - relax) * Z_s_old_arr)
            gf_Zs.vec.FV().NumPy()[:] = Z_s_blend
            dZ_per = (np.abs(Z_s_blend[bnd_mask] - Z_s_old_arr[bnd_mask])
                      / np.maximum(np.abs(Z_s_old_arr[bnd_mask]),
                                   1e-30))
            dZ = float(np.max(dZ_per)) if dZ_per.size else 0.0
            Zabs = np.abs(Z_s_blend[bnd_mask])
            Z_s_old = Z_s  # legacy: keep scalar var, used only for logging
            Z_s = complex(np.mean(Z_s_blend[bnd_mask]))  # mean for log/scalar fallback
            Ht_amp = abs(1j * omega) * At_amp_per_dof[bnd_mask] / np.maximum(
                np.abs(Z_s_old_arr[bnd_mask]), 1e-30)
            history.append({
                'iteration': iteration,
                'Z_s_abs_mean': float(np.mean(Zabs)),
                'Z_s_abs_min': float(np.min(Zabs)),
                'Z_s_abs_max': float(np.max(Zabs)),
                'H_t_per_dof_mean': float(np.mean(Ht_amp)),
                'H_t_per_dof_max': float(np.max(Ht_amp)),
                'dZ_max': dZ,
                't_solve': t_solve_iter,
            })
            _log(f"ITER:{iteration} per-panel "
                 f"<|Z_s|>={np.mean(Zabs):.4e} "
                 f"<H_t>={np.mean(Ht_amp):.2f} "
                 f"max(H_t)={np.max(Ht_amp):.2f} "
                 f"max(dZ)={dZ:.4e} t={t_solve_iter:.1f}s")
        else:
            # Scalar Karl (legacy path)
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

        # Convergence: dZ below tol.  Require iteration > 0 to avoid
        # spurious convergence on the seed Z_s, EXCEPT when max_iter
        # <= 1: the user explicitly asked for one iteration, so the
        # dZ we just computed is what they get -- report it as
        # converged if it passes tol, NOT-CONVERGED otherwise.
        if dZ < tol and (iteration > 0 or max_iter <= 1):
            _log(f"CONVERGED:iter={iteration}")
            break

        if esim is None:
            break  # linear SIBC: single iteration

    # ============================================================
    # Step 6: Post-process
    # ============================================================
    # Time-averaged power: P = 0.5 * Re(Z_s) * |H_t|_amp^2 * A_wp.
    #
    # Naming note: `H_t_rms` here is **spatial RMS of the amplitude**
    # of H_t over the workpiece surface (= sqrt(2) * time-RMS for a
    # single-frequency phasor).  The factor 0.5 in the power formula
    # corresponds to <cos^2(omega t)> = 1/2 and is correct for this
    # amplitude convention.  Do NOT replace H_t_rms with a true
    # time-RMS without dropping the 0.5 -- the same misnomer is used
    # consistently in calc_inductance, calc_fem_coilmesh, and the
    # underlying bem_sibc_solver, ESIM cell solver.
    if has_wp:
        if gf_Zs is not None:
            # Per-DOF mode: integrate the actual local Re/Im(Z_s) over
            # the workpiece surface via a CF: P = 0.5 * int Re(Z_s) *
            # (omega/|Z_s|)^2 * |A_t|^2 dS = 0.5 * int Re(Z_s) * |H_t|^2
            # dS.  We approximate by area-averaging Re(Z_s) at the BND
            # H1 DOFs (variational mean), then multiply by the global
            # H_t_rms^2 * A_wp the scalar diagnostic returns.  For Karl
            # solutions where Z_s varies <2x across the surface this
            # is within 1 % of the true integral; for stiffer variations
            # the per-DOF cumulative integral is the rigorous form (TODO).
            Z_s_arr = np.asarray(gf_Zs.vec.FV().NumPy())
            Z_s_avg_re = float(np.mean(Z_s_arr[bnd_mask].real))
            Z_s_avg_im = float(np.mean(Z_s_arr[bnd_mask].imag))
            P_total = 0.5 * Z_s_avg_re * H_t_rms**2 * A_wp
            Q_total = 0.5 * Z_s_avg_im * H_t_rms**2 * A_wp
        else:
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

    # ----------------------------------------------------------------
    # Time-averaged EM force on the SIBC workpiece [N]
    # ----------------------------------------------------------------
    # The workpiece is a HOLE (no volume mesh), so the only force handle is
    # the time-averaged Maxwell stress over its "sibc" surface.  With
    # B = curl(A_total) (PEEC: A_r + A_s) and the time-average of the
    # quadratic stress (n real, B complex phasor):
    #   <F_k> = oint_S [ (1/(2 mu0)) Re(B_k conj(B.n))
    #                     - (1/(4 mu0)) |B|^2 n_k ] dS .
    # specialcf.normal here points OUT of the air mesh = INTO the hole, the
    # OPPOSITE of the workpiece outward normal, so the force ON the workpiece
    # is the NEGATED surface integral.  The 0.5/0.25 factors (vs the static
    # 1/0.5) are the cos^2 time-average; the identity is validated by
    # reduction to the static Maxwell force in
    # packages/radia-mcp/tests/test_maxwell_surface_harmonic.py.
    F_sibc = None
    if has_wp and use_complex:
        try:
            from ngsolve import specialcf
            n_f = specialcf.normal(3)
            B_f = (curl(gfu) + curl(peec_A_s_cf)) if peec_A_s_cf is not None \
                else curl(gfu)
            Bn_f = sum(B_f[k] * n_f[k] for k in range(3))
            B2_f = sum((B_f[k] * Conj(B_f[k])).real for k in range(3))
            F_sibc = [
                -float(Integrate(
                    ((0.5 / MU_0) * (B_f[k] * Conj(Bn_f)).real
                     - (0.25 / MU_0) * B2_f * n_f[k]),
                    mesh, BND, definedon=wp_region).real)
                for k in range(3)]
            _log(f"FORCE:F_sibc=({F_sibc[0]:+.3e},{F_sibc[1]:+.3e},"
                 f"{F_sibc[2]:+.3e}) N (time-avg Maxwell stress on workpiece)")
        except Exception as _fe:  # noqa: BLE001 -- force is a diagnostic add-on
            _log(f"FORCE:skipped ({_fe})")
            F_sibc = None

    t_total = time.perf_counter() - t_total_start
    _log(f"DONE:P={P_total:.4e} L={L*1e9:.2f}nH t={t_total:.1f}s")

    # ============================================================
    # Step 6b: Surface heat flux q_surf(x,y,z) [W/m^2]
    # ============================================================
    # The SIBC boundary radiates time-averaged Poynting power per area
    #   q_surf = 0.5 * Re(Z_s) * |H_t|^2
    # where |H_t|^2 is the converged tangential field on the workpiece
    # surface.  Building it as a boundary CoefficientFunction lets us:
    #   - persist the spatial distribution as a .sol companion (so the
    #     downstream thermal solver can apply it as a Neumann BC),
    #   - emit hotspot statistics (max / mean / p95) for optimization
    #     filtering, and
    #   - cross-check that the surface integral matches P_total.
    # No-op when there is no workpiece (At_sq does not exist in scope).
    q_surf_max = None
    q_surf_mean = None
    q_surf_p95 = None
    P_total_check = None
    qsurf_sol_path = ""
    if has_wp:
        H_t_sq_cf = (omega / abs(Z_s)) ** 2 * At_sq
        q_surf_cf = 0.5 * Z_s.real * H_t_sq_cf

        # Surface integral as the cross-check for P_total above.
        P_total_check = float(
            Integrate(q_surf_cf, mesh, BND, definedon=wp_region).real)
        q_surf_mean = P_total_check / max(A_wp, 1e-30)

        # Project to a global H1 GridFunction for stats + .sol storage.
        # Set on the workpiece boundary only; the rest of the field is
        # zero and tells the loader (Phase B thermal) which DOFs are
        # active.  H1 order matches fes_order so the projection has the
        # same accuracy as the EM solve itself.
        try:
            fes_q = H1(mesh, order=fes_order)
            gf_q = GridFunction(fes_q)
            gf_q.vec[:] = 0
            gf_q.Set(q_surf_cf, definedon=wp_region)

            # Stats on the boundary DOFs.  Workpiece-surface vertices
            # are exactly the H1 DOFs whose value Set populated; we
            # tag them via the scalar 1 mask trick so the stats do not
            # depend on a numerical zero threshold.
            mask_gf = GridFunction(fes_q)
            mask_gf.vec[:] = 0
            mask_gf.Set(CF(1.0), definedon=wp_region)
            mask_arr = np.asarray(mask_gf.vec.FV().NumPy())
            q_arr = np.asarray(gf_q.vec.FV().NumPy())
            on_wp = mask_arr > 0.5
            if np.any(on_wp):
                vals = q_arr[on_wp]
                q_surf_max = float(np.max(vals))
                q_surf_p95 = float(np.percentile(vals, 95))
            else:
                q_surf_max = 0.0
                q_surf_p95 = 0.0
            _log(f"Q_SURF:max={q_surf_max:.3e} mean={q_surf_mean:.3e} "
                 f"p95={q_surf_p95:.3e} W/m^2 "
                 f"(P_check={P_total_check:.4e} vs P_total={P_total:.4e})")
        except Exception as e:
            # Stats are best-effort.  A failure here must not break the
            # main solve / JSON output.
            _log(f"Q_SURF:stats failed: {type(e).__name__}: {e}")
            gf_q = None

        # Surface current density on the SIBC face.
        #
        # |J_s| [A/m] = |H_t| is saved as a scalar for the colormap.
        #
        # Re(J_s) is saved as a 3D vector for the arrow plot.  The
        # SIBC relation E_t = Z_s * J_s combined with E_t = -j*omega*A_t
        # gives J_s = -j*omega/Z_s * A_t.  We save the real part as
        # an instantaneous-phase snapshot (the "電流分布" arrows
        # Kubota wants); GMSH renders complex AC fields by phase
        # animation only when both Re and Im are saved as separate
        # views, so for v1 the time-snapshot is enough.
        gf_J = None
        gf_J_vec = None
        try:
            J_mag_cf = sqrt(At_sq) * (omega / abs(Z_s))
            fes_J = H1(mesh, order=fes_order)
            gf_J = GridFunction(fes_J)
            gf_J.vec[:] = 0
            gf_J.Set(J_mag_cf, definedon=wp_region)
        except Exception as e:
            _log(f"J_SURF:scalar gen failed: {type(e).__name__}: {e}")
        try:
            from ngsolve import specialcf as _scf
            n_bnd_v = _scf.normal(3)
            # Tangential A as complex 3D CF.
            A_dot_n = sum(A_eval[i] * n_bnd_v[i] for i in range(3))
            A_t_vec = CF(tuple(
                A_eval[i] - A_dot_n * n_bnd_v[i] for i in range(3)))
            # -j*omega/Z_s = -j*omega*conj(Z_s)/|Z_s|^2
            # Re(-j*omega/Z_s) = -omega*Im(Z_s)/|Z_s|^2
            # Im(-j*omega/Z_s) = -omega*Re(Z_s)/|Z_s|^2
            zs_sq = Z_s.real ** 2 + Z_s.imag ** 2
            re_coef = -omega * Z_s.imag / zs_sq
            im_coef = -omega * Z_s.real / zs_sq
            # Re(J_s) = re_coef*Re(A_t) - im_coef*Im(A_t)
            J_s_re_cf = CF(tuple(
                re_coef * A_t_vec[i].real - im_coef * A_t_vec[i].imag
                for i in range(3)))
            fes_J_vec = H1(mesh, order=fes_order, dim=3)
            gf_J_vec = GridFunction(fes_J_vec)
            gf_J_vec.vec[:] = 0
            gf_J_vec.Set(J_s_re_cf, definedon=wp_region)
        except Exception as e:
            _log(f"J_SURF:vector gen failed: {type(e).__name__}: {e}")
            gf_J_vec = None
    else:
        gf_q = None
        gf_J = None
        gf_J_vec = None

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

            # Surface heat-flux distribution q_surf [W/m^2] on the
            # workpiece SIBC face.  Phase B (calc_heat.py) loads this
            # .sol and applies it as the Neumann BC for the heat solve.
            # Saved at the same H1 order as the EM solve so the
            # spatial detail is preserved in the round trip.
            if gf_q is not None:
                sol_Q = os.path.join(base_dir,
                                     f"{name_stem}_qsurf.sol").replace("\\", "/")
                gf_q.Save(sol_Q)
                qsurf_sol_path = sol_Q
                sol_paths["q_surf"] = sol_Q
                sol_entries.append(
                    {"sol": sol_Q, "fes": "H1", "fes_order": fes_order,
                     "fes_dim": 1, "name": "q_surf", "ncomp": 1})
            if gf_J is not None:
                sol_J = os.path.join(base_dir,
                                     f"{name_stem}_Jsurf.sol").replace("\\", "/")
                gf_J.Save(sol_J)
                sol_paths["J_surf"] = sol_J
                sol_entries.append(
                    {"sol": sol_J, "fes": "H1", "fes_order": fes_order,
                     "fes_dim": 1, "name": "J_surf_Am", "ncomp": 1})
            if gf_J_vec is not None:
                sol_Jv = os.path.join(base_dir,
                                      f"{name_stem}_Jvec.sol").replace("\\", "/")
                gf_J_vec.Save(sol_Jv)
                sol_paths["J_surf_vec"] = sol_Jv
                sol_entries.append(
                    {"sol": sol_Jv, "fes": "H1", "fes_order": fes_order,
                     "fes_dim": 3, "name": "J_surf_vec",
                     "ncomp": 3})
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
        # Time-averaged EM force on the SIBC workpiece [N] (Maxwell stress
        # over the "sibc" surface).  None for DC / no-workpiece runs.
        "F_sibc": F_sibc,
        # Surface heat-flux statistics for optimization filtering.
        # All keys present (None when no workpiece) so downstream
        # consumers can rely on the schema.
        "q_surf_max": q_surf_max,
        "q_surf_mean": q_surf_mean,
        "q_surf_p95": q_surf_p95,
        "P_total_check": P_total_check,
        "qsurf_sol": qsurf_sol_path,
        "delta": float(delta_skin),
        "ndof": ndof,
        "ne": ne,
        "iterations": len(history),
        "esim_history": history,
        "esim_per_panel": bool(gf_Zs is not None),
        "t_mesh_s": round(t_mesh, 2),
        # Aggregated per-iteration solve time.  Each iteration's t_solve is
        # already captured in history[i]['t_solve']; sum into a single
        # top-level t_solve_s so the notebook/workbench summary shows the
        # solve breakdown alongside t_mesh_s + t_total_s.  Per-iter detail
        # is still available in esim_history for users who want it.
        "t_solve_s": round(sum(h.get("t_solve", 0.0) for h in history), 2),
        "t_total_s": round(t_total, 2),
        "frequency": frequency,
        "sigma": sigma,
        "impedance_model": impedance_model,
        "formulation": formulation,
        "source_type": "PEEC",
        "has_kelvin": has_kelvin,
        "msh_file": gmsh_file,
    }
    if gf_Zs is not None:
        Z_s_arr = np.asarray(gf_Zs.vec.FV().NumPy())
        Z_s_bnd = Z_s_arr[bnd_mask]
        result["esim_per_panel_Z_s_real"] = Z_s_bnd.real.tolist()
        result["esim_per_panel_Z_s_imag"] = Z_s_bnd.imag.tolist()
        result["esim_per_panel_n_dof"] = int(n_bnd_dofs)
    return result


def build_argparser():
    """argparse factory shared by main() and notebook DesignSpec callers."""
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
    parser.add_argument("--esim-per-panel", action="store_true",
                        help="Per-BND-DOF ESIM Karl (v4.48.0+): one ESIM "
                             "cell solve per H1 workpiece-boundary DOF "
                             "using a per-DOF |H_t| extracted from the "
                             "Karl-iter solution.  Resolves spatial "
                             "saturation pattern; ~N_bnd_dofs × cell-solve "
                             "cost per Karl iter.")
    parser.add_argument("--solver", default="auto",
                        choices=["auto", "pardiso", "bddc", "iccg", "ams"],
                        help="auto: ams for fes-order=1, bddc otherwise. "
                             "pardiso (direct), bddc (iterative), "
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
    parser.add_argument("--peec-n-peri", type=int, default=0,
                        help="If >0, use perimeter-only filament "
                             "placement (path 2b: circle-edge centres) "
                             "instead of the volume-grid walker.  "
                             "Required for 3-turn helix / radia "
                             "STEPs whose seed plane is hard to find.")
    parser.add_argument("--require-kelvin", action="store_true",
                        help="Fail-fast if the .vol does not have a "
                             "'kelvin' material with periodic "
                             "identification.  The IH panel's "
                             "PEEC+FEM+Kelvin method sets this to "
                             "prevent silent fallback to "
                             "reg-only gauge truncation (which "
                             "produces a physically wrong open-BC "
                             "result that looks like a successful "
                             "solve).  Per CLAUDE.md 'No Fallbacks'.")
    return parser


def main():
    parser = build_argparser()

    def run(args):
        if not args.peec_step:
            return {"error": "--peec-step is required. T0/AV paths were "
                             "retired 2026-04-18; all FEM source injection "
                             "goes through PEEC filaments."}
        if args.require_kelvin:
            # Sniff the .vol's materials list before the solver gets a
            # chance to silently downgrade to reg-only gauge truncation.
            # Cheap O(1) check using Netgen's mesh loader; full mesh
            # construction happens inside solve_fem anyway.
            from ngsolve import Mesh
            try:
                probe_mesh = Mesh(args.vol)
            except Exception as e:
                return {"error":
                        f"--require-kelvin: cannot open {args.vol!r}: {e}"}
            mats = list(probe_mesh.GetMaterials())
            if "kelvin" not in mats:
                return {"error":
                        f"--require-kelvin: .vol has no 'kelvin' material "
                        f"(materials = {mats}). The IH panel's "
                        f"PEEC+FEM+Kelvin method requires a Kelvin-extended "
                        f"mesh: re-run the Cubit .jou with "
                        f"add_kelvin enabled, OR switch to a method "
                        f"that doesn't need the open boundary."}
            if probe_mesh.ngmesh.GetNrIdentifications() == 0:
                return {"error":
                        f"--require-kelvin: .vol has 'kelvin' material but "
                        f"no periodic identifications. The Kelvin pair "
                        f"sidesets (kelvin_int, kelvin_ext) need 'copy "
                        f"mesh surface' in the Cubit .jou so the C++ "
                        f"exporter can write Netgen identification "
                        f"records."}
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
            peec_n_peri=args.peec_n_peri,
            esim_per_panel=args.esim_per_panel,
        )

    calc_main(run, parser)


if __name__ == "__main__":
    main()
