"""calc_fem_coilmesh.py -- A-V formulation: coil meshed + wp SIBC + Kelvin.

Layer 4 subprocess calc for the IH panel "FEM (coil meshed + SIBC +
Kelvin)" method.  Gapped torus only (real IH coils have physical
port terminations; closed-torus is a topological abstraction).

Formulation (following MCP `INDUCTION_HEATING_AV_COIL_SIGMA`):

  FES: HCurl(A) x H1(phi)   phi defined only on coil material.
  Dirichlet: phi=1 on 'source', phi=0 on 'sink'.  Solve, then scale.
  Current extraction: VOLUME integral I_out = int J . grad(psi_n) dV
     where psi_n is a scalar H1 test function with psi_n=1 on source,
     psi_n=0 on sink.  By Gauss' theorem this equals the surface flux
     int_source J . n dS and is FEM-consistent (preferred over direct
     cut-plane surface integral).

Required .vol:
  materials: 'coil' + air (+ 'kelvin' for open boundary)
  boundaries: 'source' (one gap face), 'sink' (other gap face),
              'sibc' (workpiece hole).
  coil mesh size <= delta_coil / 3 (~0.26 mm for Cu at 7 kHz).

Output: JSON to stdout.  L, P_total, P_coil, P_wp, I_out (pre-scale).
"""
from __future__ import annotations

import argparse
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SRC_RADIA = os.path.abspath(os.path.join(HERE, ".."))
if SRC_RADIA not in sys.path:
    sys.path.insert(0, SRC_RADIA)
if HERE not in sys.path:
    sys.path.insert(0, HERE)

from calc_common import calc_main, progress
from em_material import EMMaterial, MU_0

NU_0 = 1.0 / MU_0


def solve_fem_coilmesh(vol, frequency, I_target,
                       coil_sigma,
                       wp_sigma, wp_mu_r, half_thickness,
                       fes_order=1, solver="pardiso",
                       sibc_bnd="sibc",
                       source_bnd="source", sink_bnd="sink",
                       coil_mat="coil",
                       msh_output="",
                       impedance_model="sibc",
                       bh_file="",
                       esim_max_iter=15,
                       esim_tol=1e-3,
                       esim_relax=0.5,
                       esim_per_panel=False):
    """A-V formulation for volumetric coil + SIBC workpiece + Kelvin.

    ``impedance_model="esim"`` runs a Karl iteration: per outer iteration
    the BilinearForm is re-assembled with the current ``Z_s_wp`` (the
    Robin SIBC term ``robin_wp = s / Z_s_wp`` depends on it), the A-V
    system is solved, ``H_t_rms_wp`` is read off the workpiece surface,
    and ``Z_s_wp`` is refreshed by ``esim.solve(H_t_rms)``.  ``bh_file``
    is required for ESIM (a 2-column ``H[A/m] B[T]`` table).
    """
    import radia  # noqa: F401  DLL path setup

    from ngsolve import (Mesh, HCurl, H1, Periodic, BilinearForm, LinearForm,
                         GridFunction, Integrate, Conj, curl, grad, dx, ds,
                         CF, BND, VOL, TaskManager, sqrt as ngsqrt, IfPos,
                         x, y, z, specialcf)
    from calc_fem_kelvin import detect_kelvin_offset

    omega = 2 * math.pi * frequency
    s = 1j * omega

    bh_curve = None
    if bh_file:
        from em_material import load_bh_file
        bh_curve = load_bh_file(bh_file)
    if impedance_model == "esim" and bh_curve is None:
        raise ValueError(
            "impedance_model='esim' requires bh_file with a 2-column "
            "BH curve (H[A/m], B[T]).")

    wp_mat = EMMaterial(name="wp", sigma=wp_sigma, mu_r=wp_mu_r,
                         bh_curve=bh_curve)
    delta_wp = wp_mat.skin_depth(frequency)
    if impedance_model == "esim":
        # Seed Z_s with a small-H Picard solve; cap inner iter at 5
        # because the outer Karl loop will refresh Z_s immediately.
        esim_solver = wp_mat.create_esim_solver(
            frequency, half_thickness, geometry='cylinder')
        Z_s_wp = complex(esim_solver.solve(5.0, max_iter=5)['Z'])
    else:
        esim_solver = None
        Z_s_wp = wp_mat.dowell_Zs(frequency, half_thickness)
    robin_wp = s / Z_s_wp

    # Non-magnetic coil convention (Cu/Al only -- coil_mu_r=1).
    coil_delta = math.sqrt(2.0 / (omega * MU_0 * coil_sigma))

    progress("FEM", f"load {os.path.basename(vol)}")
    t0 = time.perf_counter()
    mesh = Mesh(vol)
    materials = mesh.GetMaterials()
    boundaries = set(mesh.GetBoundaries())
    t_load = time.perf_counter() - t0
    progress("FEM", f"ne={mesh.ne} mats={set(materials)} ({t_load:.1f}s)")

    # Validate required labels
    if coil_mat not in materials:
        raise ValueError(
            f"A-V FEM requires material {coil_mat!r} (meshed coil).  "
            f"Available: {sorted(set(materials))}.")
    for req, kind in ((source_bnd, "source port"),
                       (sink_bnd, "sink port"),
                       (sibc_bnd, "workpiece SIBC")):
        if req not in boundaries:
            raise ValueError(
                f"A-V FEM requires boundary {req!r} ({kind}).  "
                f"Real IH coils have physical port terminations; "
                f"the .jou must tag the gap faces as '{source_bnd}' / "
                f"'{sink_bnd}'.  Available: {sorted(boundaries)}.")

    # Coil mesh resolution check (advisory warning)
    # Sample a few coil tets to get a ballpark edge length.
    n_coil = 0
    coil_h_max = 0.0
    for el in mesh.Elements(VOL):
        if el.mat == coil_mat:
            pts = np.array([mesh.vertices[v.nr].point for v in el.vertices])
            # max pairwise distance
            d = float(np.linalg.norm(pts[None] - pts[:, None], axis=-1).max())
            coil_h_max = max(coil_h_max, d)
            n_coil += 1
            if n_coil >= 200:
                break
    if coil_h_max > coil_delta:
        progress("FEM",
                 f"WARN coil h_max={coil_h_max*1e3:.2f} mm > "
                 f"delta={coil_delta*1e3:.2f} mm; skin effect under-resolved")

    # Kelvin detection + weight
    has_kelvin = "kelvin" in materials
    has_kelvin_periodic = False
    a_kelvin = 0.0
    kelvin_center = np.array([0.0, 0.0, 0.0])
    if has_kelvin:
        kelvin_center = np.array(detect_kelvin_offset(mesh))
        kelvin_verts = set()
        for el in mesh.Elements(VOL):
            if el.mat == "kelvin":
                for v in el.vertices:
                    kelvin_verts.add(v.nr)
        coords = np.array([mesh.vertices[v].point for v in kelvin_verts])
        dists = np.linalg.norm(coords - kelvin_center[None, :], axis=1)
        a_kelvin = float(np.max(dists))
        has_kelvin_periodic = mesh.ngmesh.GetNrIdentifications() > 0
        progress("FEM",
                 f"Kelvin R={a_kelvin*1e3:.1f}mm periodic={has_kelvin_periodic}")

    # nu per material: NU_0 everywhere; Kelvin exterior gets (rho'/a_k)^2.
    kx, ky, kz = kelvin_center
    nu_dict = {}
    for m in materials:
        if "kelvin" in m.lower():
            dxk, dyk, dzk = x - kx, y - ky, z - kz
            rp_sq = dxk * dxk + dyk * dyk + dzk * dzk + 1e-20
            nu_dict[m] = NU_0 * rp_sq / a_kelvin ** 2
        else:
            nu_dict[m] = NU_0  # non-magnetic coil + air share NU_0
    nu_cf = mesh.MaterialCF(nu_dict, default=NU_0)

    # Compound FES: HCurl(A) x H1(phi on coil).
    #
    # NB: HCurl A has no Dirichlet here.  An earlier version passed
    # dirichlet="GND" when a "GND" boundary was present; empirically
    # (2026-04-28) HCurl(dirichlet=<vertex_tag>) is a no-op because
    # HCurl DOFs live on edges, so the path was silently dead.  The
    # gauge null-space of HCurl A in the A-phi system is locked by:
    #   - phi's source/sink Dirichlet drives the coil current;
    #   - nograds=True removes pure-gradient modes;
    #   - the coil's j*omega*sigma*A term provides positive-definite
    #     gauge fixing in the conducting region (AC);
    #   - reg=1e-6 on the air mass term covers DC and harmonic-form
    #     residuals.
    fesA_base = HCurl(mesh, order=fes_order, nograds=True, complex=True)
    fesA = Periodic(fesA_base) if has_kelvin_periodic else fesA_base
    fesPhi = H1(mesh, order=fes_order, complex=True,
                definedon=mesh.Materials(coil_mat),
                dirichlet=f"{source_bnd}|{sink_bnd}")
    fes = fesA * fesPhi
    (A, phi), (N, psi) = fes.TnT()
    progress("FEM", f"ndof={fes.ndof} (A:{fesA.ndof} + phi:{fesPhi.ndof})")

    # Bilinear form pieces (Robin term depends on Z_s_wp and is rebuilt
    # per Karl iteration; everything else is geometry/material-only and
    # is reassembled trivially because BilinearForm composition is cheap).
    non_kelvin = [m for m in materials if "kelvin" not in m.lower()]

    def _assemble_a_bf(robin_coeff):
        """(Re)build + assemble the A-V BilinearForm at given Robin coef.

        Re-built per Karl iteration because robin_wp = s / Z_s_wp varies
        with the updated SIBC impedance.
        """
        a_bf = BilinearForm(fes, symmetric=True)
        a_bf += nu_cf * curl(A) * curl(N) * dx(bonus_intorder=4)
        if non_kelvin:
            a_bf += 1e-6 * NU_0 * A * N * dx("|".join(non_kelvin))
        a_bf += robin_coeff * A.Trace() * N.Trace() * ds(sibc_bnd)
        a_bf += s * coil_sigma * (A + grad(phi)) * (N + grad(psi)) * \
                dx(coil_mat)
        with TaskManager():
            a_bf.Assemble()
        return a_bf

    # Test function psi_n for I_out volume-integral extraction
    # (geometry-only, independent of Z_s_wp).
    fes_psi_n = H1(mesh, order=1, complex=True,
                    definedon=mesh.Materials(coil_mat),
                    dirichlet=sink_bnd)
    gf_psi_n = GridFunction(fes_psi_n)
    gf_psi_n.Set(CF(1), definedon=mesh.Boundaries(source_bnd))

    # Pre-compute wp surface (geometry-only); used inside the loop for
    # H_t_rms_wp and in the post-step.
    wp_region = mesh.Boundaries(sibc_bnd)
    A_wp = float(Integrate(CF(1), mesh, BND, definedon=wp_region).real)

    gfu = GridFunction(fes)
    gf_A, gf_phi = gfu.components
    n_bnd = specialcf.normal(3)

    # Per-panel ESIM (v4.48.0+): build a surface H1 GridFunction on the
    # workpiece BND so the Robin coefficient becomes a per-DOF CF.  Same
    # pattern as calc_fem_kelvin.py.
    if esim_per_panel and esim_solver is not None:
        fes_Zs = H1(mesh, order=fes_order, complex=True)
        gf_Zs = GridFunction(fes_Zs)
        gf_Zs.vec[:] = 0
        gf_Zs.Set(CF(Z_s_wp), definedon=wp_region)
        mask_gf = GridFunction(H1(mesh, order=fes_order))
        mask_gf.vec[:] = 0
        mask_gf.Set(CF(1.0), definedon=wp_region)
        bnd_mask = np.asarray(mask_gf.vec.FV().NumPy()) > 0.5
        n_bnd_dofs = int(bnd_mask.sum())
        progress("FEM",
            f"PER_PANEL: wp surface H1 dofs={n_bnd_dofs}")
    else:
        fes_Zs = None
        gf_Zs = None
        bnd_mask = None
        n_bnd_dofs = 0

    esim_history = []
    esim_converged = (esim_solver is None)
    max_iter = max(int(esim_max_iter), 1) if esim_solver is not None else 1
    t_asm = 0.0
    t_solve = 0.0
    n_iter_done = 0
    dZ = float("inf")
    H_t_rms_iter = 0.0
    a_bf = None

    for iteration in range(max_iter):
        n_iter_done = iteration + 1
        if gf_Zs is not None:
            robin_wp = s / gf_Zs  # CF
            progress("FEM",
                f"assemble A-V (iter {iteration}, per-panel, "
                f"<|Z_s|>={float(np.mean(np.abs(np.asarray(gf_Zs.vec.FV().NumPy())[bnd_mask]))):.3e})")
        else:
            robin_wp = s / Z_s_wp
            progress("FEM",
                f"assemble A-V (iter {iteration}, Z_s={Z_s_wp:.3e})")
        t0 = time.perf_counter()
        a_bf = _assemble_a_bf(robin_wp)
        t_asm_iter = time.perf_counter() - t0
        t_asm += t_asm_iter

        # Reset Dirichlet lift each iteration
        gfu.vec[:] = 0
        gf_phi.Set(CF(1), definedon=mesh.Boundaries(source_bnd))

        progress("FEM", f"solve ({solver}, iter {iteration})")
        t0 = time.perf_counter()
        with TaskManager():
            r = gfu.vec.CreateVector()
            r.data = -a_bf.mat * gfu.vec
            gfu.vec.data += a_bf.mat.Inverse(fes.FreeDofs(),
                                               inverse=solver) * r
        t_solve_iter = time.perf_counter() - t0
        t_solve += t_solve_iter

        # Scale gfu to physical I_target so H_t_rms reflects the actual
        # operating point ESIM is supposed to see.  Without this Karl
        # iterates around the Dirichlet-phi=1 system, whose H_t depends
        # on workpiece back-reaction (I_out) and Z_s_wp simultaneously
        # — the Karl fixed point would no longer be physical.
        J_coil = -s * coil_sigma * (gf_A + grad(gf_phi))
        I_out_pre = complex(Integrate(
            J_coil * grad(gf_psi_n), mesh,
            definedon=mesh.Materials(coil_mat)))
        if abs(I_out_pre) < 1e-20:
            raise RuntimeError(
                f"FEM I_out is zero at Karl iter {iteration}; A-V setup "
                f"failed.  Check source/sink labels and coil mesh.")
        scale = complex(I_target) / I_out_pre
        gfu.vec.data = complex(scale) * gfu.vec

        if esim_solver is None:
            progress("FEM",
                f"solved ({t_solve_iter:.1f}s, I_out={abs(I_out_pre):.4e})")
            break

        # H_t_rms_wp from the just-scaled gfu (physical magnitude).
        A_sq = sum(gf_A[i].real ** 2 + gf_A[i].imag ** 2 for i in range(3))
        Adn_re = sum(gf_A[i].real * n_bnd[i] for i in range(3))
        Adn_im = sum(gf_A[i].imag * n_bnd[i] for i in range(3))
        An_sq = Adn_re ** 2 + Adn_im ** 2
        At_sq_loop = A_sq - An_sq
        At_int_wp_loop = float(Integrate(
            At_sq_loop, mesh, BND, definedon=wp_region).real)
        At_rms_loop = math.sqrt(max(At_int_wp_loop, 0.0) / max(A_wp, 1e-30))

        # Karl update with under-relaxation (matches calc_fem_kelvin
        # default 0.5; damps oscillation near saturation).
        if gf_Zs is not None:
            # Per-DOF Karl: project |A_t|^2 per-DOF, ESIM per-DOF.
            fes_real = H1(mesh, order=fes_order)
            gf_At_re = GridFunction(fes_real)
            gf_At_im = GridFunction(fes_real)
            gf_At_re.vec[:] = 0
            gf_At_im.vec[:] = 0
            A_sq_re = sum(gf_A[i].real ** 2 for i in range(3))
            A_sq_im = sum(gf_A[i].imag ** 2 for i in range(3))
            Adn_re_sq = Adn_re ** 2
            Adn_im_sq = Adn_im ** 2
            At_sq_re = A_sq_re - Adn_re_sq
            At_sq_im = A_sq_im - Adn_im_sq
            gf_At_re.Set(At_sq_re, definedon=wp_region)
            gf_At_im.Set(At_sq_im, definedon=wp_region)
            At_re_arr = np.asarray(gf_At_re.vec.FV().NumPy())
            At_im_arr = np.asarray(gf_At_im.vec.FV().NumPy())
            At_amp_per_dof = np.sqrt(
                np.maximum(At_re_arr, 0.0) + np.maximum(At_im_arr, 0.0))
            Z_s_old_arr = np.asarray(gf_Zs.vec.FV().NumPy()).copy()
            Z_s_new_arr = Z_s_old_arr.copy()
            n_called = 0
            for i in np.flatnonzero(bnd_mask):
                Zsi_old = Z_s_old_arr[i]
                if abs(Zsi_old) <= 1e-30:
                    continue
                Ht_i = abs(s / Zsi_old) * float(At_amp_per_dof[i])
                sol_new = esim_solver.solve(max(Ht_i, 1e-3))
                Z_s_new_arr[i] = complex(sol_new['Z'])
                n_called += 1
            Z_s_blend = (esim_relax * Z_s_new_arr
                         + (1 - esim_relax) * Z_s_old_arr)
            gf_Zs.vec.FV().NumPy()[:] = Z_s_blend
            dZ_per = (np.abs(Z_s_blend[bnd_mask] - Z_s_old_arr[bnd_mask])
                      / np.maximum(np.abs(Z_s_old_arr[bnd_mask]), 1e-30))
            dZ = float(np.max(dZ_per)) if dZ_per.size else 0.0
            Zabs = np.abs(Z_s_blend[bnd_mask])
            Z_s_avg = complex(np.mean(Z_s_blend[bnd_mask]))
            Z_s_wp = Z_s_avg  # for downstream legacy uses + log
            Ht_amp = abs(s) * At_amp_per_dof[bnd_mask] / np.maximum(
                np.abs(Z_s_old_arr[bnd_mask]), 1e-30)
            H_t_rms_iter = float(np.mean(Ht_amp))
            esim_history.append({
                "iteration": iteration,
                "Z_s_abs_mean": float(np.mean(Zabs)),
                "Z_s_abs_min": float(np.min(Zabs)),
                "Z_s_abs_max": float(np.max(Zabs)),
                "H_t_per_dof_mean": float(np.mean(Ht_amp)),
                "H_t_per_dof_max": float(np.max(Ht_amp)),
                "dZ_max": dZ,
                "t_solve": float(t_solve_iter),
            })
            progress("FEM",
                f"ESIM:ITER {iteration} per-panel "
                f"<|Z_s|>={np.mean(Zabs):.4e} "
                f"<H_t>={np.mean(Ht_amp):.2f} "
                f"max(H_t)={np.max(Ht_amp):.2f} "
                f"max(dZ)={dZ:.4e} "
                f"t={t_asm_iter+t_solve_iter:.1f}s")
        else:
            H_t_rms_iter = abs(s / Z_s_wp) * At_rms_loop
            Z_s_old = Z_s_wp
            sol_new = esim_solver.solve(max(float(H_t_rms_iter), 1e-3))
            Z_s_new = complex(sol_new['Z'])
            Z_s_wp = esim_relax * Z_s_new + (1 - esim_relax) * Z_s_old
            dZ = abs(Z_s_wp - Z_s_old) / max(abs(Z_s_old), 1e-30)
            esim_history.append({
                "iteration": iteration,
                "Z_s_abs": float(abs(Z_s_wp)),
                "H_t_rms": float(H_t_rms_iter),
                "dZ": float(dZ),
                "t_solve": float(t_solve_iter),
            })
            progress("FEM",
                f"ESIM:ITER {iteration} |Z_s|={abs(Z_s_wp):.4e} "
                f"H_t={H_t_rms_iter:.2f} dZ={dZ:.4e} "
                f"t={t_asm_iter+t_solve_iter:.1f}s")
        # Require iteration > 0 to avoid spurious convergence on the
        # seed Z_s (= esim.solve(5.0)), EXCEPT when esim_max_iter
        # <= 1: the user explicitly asked for one iteration so the
        # dZ we just computed is the only result available.
        if dZ < esim_tol and (iteration > 0 or esim_max_iter <= 1):
            esim_converged = True
            progress("FEM", f"ESIM:CONVERGED iter={iteration}")
            break
    else:
        if esim_solver is not None:
            esim_converged = False
            progress("FEM",
                f"ESIM:NOT-CONVERGED after {max_iter} iter "
                f"(dZ={dZ:.4e} > tol={esim_tol:.1e})")
    # After ESIM convergence, re-solve once at the converged Z_s so
    # downstream post-processing (energy, P_wp, q_surf) sees the same
    # gfu that produced the final Karl Z_s.
    if esim_solver is not None:
        if gf_Zs is not None:
            progress("FEM", f"final re-solve at Z_s (per-panel)")
            robin_final = s / gf_Zs
        else:
            progress("FEM", f"final re-solve at Z_s={Z_s_wp:.3e}")
            robin_final = s / Z_s_wp
        t0 = time.perf_counter()
        a_bf = _assemble_a_bf(robin_final)
        t_asm += time.perf_counter() - t0
        gfu.vec[:] = 0
        gf_phi.Set(CF(1), definedon=mesh.Boundaries(source_bnd))
        t0 = time.perf_counter()
        with TaskManager():
            r = gfu.vec.CreateVector()
            r.data = -a_bf.mat * gfu.vec
            gfu.vec.data += a_bf.mat.Inverse(fes.FreeDofs(),
                                               inverse=solver) * r
        t_solve += time.perf_counter() - t0
        J_coil = -s * coil_sigma * (gf_A + grad(gf_phi))
        I_out_pre = complex(Integrate(
            J_coil * grad(gf_psi_n), mesh,
            definedon=mesh.Materials(coil_mat)))
        if abs(I_out_pre) < 1e-20:
            raise RuntimeError("FEM I_out is zero on final re-solve.")
        scale = complex(I_target) / I_out_pre
        gfu.vec.data = complex(scale) * gfu.vec

    progress("FEM",
        f"FEM total (asm {t_asm:.1f}s + solve {t_solve:.1f}s "
        f"over {n_iter_done} iter)")

    # Recompute J_coil CF at final physical-scale gfu for downstream
    # post-processing (P_coil, J_surf export, etc.).
    J_coil = -s * coil_sigma * (gf_A + grad(gf_phi))

    # Post-solve quantities
    # L from volumetric curl energy
    W_vol = float(Integrate(0.5 * nu_cf * curl(gf_A) * Conj(curl(gf_A)),
                              mesh, order=10).real)
    L_vol = 2 * W_vol / I_target ** 2

    # WP dissipation via SIBC surface integral
    n_bnd = specialcf.normal(3)
    A_sq = sum(gf_A[i].real ** 2 + gf_A[i].imag ** 2 for i in range(3))
    Adn_re = sum(gf_A[i].real * n_bnd[i] for i in range(3))
    Adn_im = sum(gf_A[i].imag * n_bnd[i] for i in range(3))
    An_sq = Adn_re ** 2 + Adn_im ** 2
    At_sq = A_sq - An_sq
    wp_region = mesh.Boundaries(sibc_bnd)
    A_wp = float(Integrate(CF(1), mesh, BND, definedon=wp_region).real)
    At_int_wp = float(Integrate(At_sq, mesh, BND, definedon=wp_region).real)
    H_t_rms_wp = abs(s / Z_s_wp) * math.sqrt(max(At_int_wp, 0.0) / A_wp)
    P_wp = 0.5 * Z_s_wp.real * H_t_rms_wp ** 2 * A_wp
    L_skin_wp = omega * Z_s_wp.imag / (abs(Z_s_wp) ** 2) * At_int_wp \
                 / I_target ** 2

    # Coil dissipation via volumetric |J|^2/sigma.  Requires the coil
    # mesh to resolve the skin depth (h_coil <= delta_coil / 3).  Coarser
    # meshes over-estimate P_coil (seen ~1.8x on h=1mm / delta=0.79mm).
    # Trust A-V + volumetric — mesh resolution is the engineering
    # responsibility of the .jou author, not a solver workaround.
    J_sq = (J_coil[0].real ** 2 + J_coil[0].imag ** 2
            + J_coil[1].real ** 2 + J_coil[1].imag ** 2
            + J_coil[2].real ** 2 + J_coil[2].imag ** 2)
    P_coil = 0.5 / coil_sigma * float(Integrate(
        J_sq, mesh, definedon=mesh.Materials(coil_mat),
        order=10).real)
    progress("FEM",
             f"P_coil (volumetric |J|^2/sigma) = {P_coil:.3e} W")

    P_total = P_coil + P_wp
    L_total = L_vol + L_skin_wp

    # ============================================================
    # Field export to GMSH .msh: B + q_surf + J_surf (mirror of
    # calc_fem_kelvin.py save block).  Without this, the panel's
    # Open GMSH button shows a mesh with NO views attached -- kubota
    # reported "qsurf not in view list" for FEM-full because that
    # mode had never produced any field views (calc_fem_coilmesh
    # used to call vol2msh with an empty fields list).
    # ============================================================
    q_surf_max = None
    q_surf_p95 = None
    q_surf_mean = None
    P_total_check = None
    qsurf_sol_path = ""
    gmsh_file = ""
    msh_export_error = None
    if msh_output:
        from ngsolve import HDiv
        try:
            H_t_sq_cf = (omega / abs(Z_s_wp)) ** 2 * At_sq
            q_surf_cf = 0.5 * Z_s_wp.real * H_t_sq_cf
            P_total_check = float(
                Integrate(q_surf_cf, mesh, BND, definedon=wp_region).real)
            q_surf_mean = P_total_check / max(A_wp, 1e-30)
            # H1 scalar field; non-wp DOFs remain at 0 by .Set's
            # definedon=wp_region restriction.  Thermal Phase B loads
            # this same .sol as the Neumann BC source.
            fes_q = H1(mesh, order=fes_order)
            gf_q = GridFunction(fes_q)
            gf_q.vec[:] = 0
            gf_q.Set(q_surf_cf, definedon=wp_region)
            # Stats via mask trick (matches calc_fem_kelvin).
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
            progress("FEM",
                     f"Q_SURF max={q_surf_max:.3e} mean={q_surf_mean:.3e} "
                     f"p95={q_surf_p95:.3e} W/m^2 "
                     f"(P_check={P_total_check:.4e} vs P_total={P_total:.4e})")
        except Exception as e:
            progress("FEM", f"Q_SURF stats failed: {type(e).__name__}: {e}")
            gf_q = None
        # |J_s| = |H_t| on wp surface (scalar) and Re(J_s) (vector).
        gf_J = None
        gf_J_vec = None
        try:
            J_mag_cf = ngsqrt(At_sq) * (omega / abs(Z_s_wp))
            fes_J = H1(mesh, order=fes_order)
            gf_J = GridFunction(fes_J)
            gf_J.vec[:] = 0
            gf_J.Set(J_mag_cf, definedon=wp_region)
        except Exception as e:
            progress("FEM",
                     f"J_SURF scalar gen failed: {type(e).__name__}: {e}")
        try:
            n_bnd_v = specialcf.normal(3)
            A_dot_n = sum(gf_A[i] * n_bnd_v[i] for i in range(3))
            A_t_vec = CF(tuple(
                gf_A[i] - A_dot_n * n_bnd_v[i] for i in range(3)))
            zs_sq = Z_s_wp.real ** 2 + Z_s_wp.imag ** 2
            re_coef = -omega * Z_s_wp.imag / zs_sq
            im_coef = -omega * Z_s_wp.real / zs_sq
            J_s_re_cf = CF(tuple(
                re_coef * A_t_vec[i].real - im_coef * A_t_vec[i].imag
                for i in range(3)))
            fes_J_vec = H1(mesh, order=fes_order, dim=3)
            gf_J_vec = GridFunction(fes_J_vec)
            gf_J_vec.vec[:] = 0
            gf_J_vec.Set(J_s_re_cf, definedon=wp_region)
        except Exception as e:
            progress("FEM",
                     f"J_SURF vector gen failed: {type(e).__name__}: {e}")
            gf_J_vec = None
        # Save .vol + .sol pairs and build the GMSH .msh via vol2msh.
        try:
            import sys as _sys
            import os as _os
            radia_src = _os.path.dirname(_os.path.abspath(__file__)) + "/.."
            if _os.path.abspath(radia_src) not in _sys.path:
                _sys.path.insert(0, _os.path.abspath(radia_src))
            from gmsh_post_export import save_vol_sol_pair, vol2msh
            base_dir = _os.path.dirname(_os.path.abspath(msh_output))
            name_stem = _os.path.splitext(_os.path.basename(msh_output))[0]
            # Real-component B field (HDiv order=1 with H1(dim=3) fallback).
            B_cf = CF(tuple(curl(gf_A)[i].real for i in range(3)))
            fes_B = HDiv(mesh, order=1)
            gf_B = GridFunction(fes_B)
            try:
                gf_B.Set(B_cf)
            except Exception:
                fes_B = H1(mesh, order=1, dim=3)
                gf_B = GridFunction(fes_B)
                gf_B.Set(B_cf)
            vol_B = _os.path.join(base_dir, f"{name_stem}_fem.vol").replace("\\", "/")
            sol_B = _os.path.join(base_dir, f"{name_stem}_B.sol").replace("\\", "/")
            save_vol_sol_pair(vol_B, sol_B, mesh.ngmesh, gf_B)
            sol_entries = [
                {"sol": sol_B, "fes": type(fes_B).__name__,
                 "fes_order": 1,
                 "fes_dim": getattr(fes_B, "dim", 1),
                 "name": "B", "ncomp": 3},
            ]
            if gf_q is not None:
                sol_Q = _os.path.join(base_dir,
                                      f"{name_stem}_qsurf.sol").replace("\\", "/")
                gf_q.Save(sol_Q)
                qsurf_sol_path = sol_Q
                sol_entries.append(
                    {"sol": sol_Q, "fes": "H1", "fes_order": fes_order,
                     "fes_dim": 1, "name": "q_surf", "ncomp": 1})
            if gf_J is not None:
                sol_J = _os.path.join(base_dir,
                                      f"{name_stem}_Jsurf.sol").replace("\\", "/")
                gf_J.Save(sol_J)
                sol_entries.append(
                    {"sol": sol_J, "fes": "H1", "fes_order": fes_order,
                     "fes_dim": 1, "name": "J_surf_Am", "ncomp": 1})
            if gf_J_vec is not None:
                sol_Jv = _os.path.join(base_dir,
                                       f"{name_stem}_Jvec.sol").replace("\\", "/")
                gf_J_vec.Save(sol_Jv)
                sol_entries.append(
                    {"sol": sol_Jv, "fes": "H1", "fes_order": fes_order,
                     "fes_dim": 3, "name": "J_surf_vec", "ncomp": 3})
            vol2msh(msh_output, vol_B, sol_entries)
            gmsh_file = msh_output
            progress("FEM",
                     f"GMSH wrote {_os.path.basename(msh_output)} "
                     f"({len(sol_entries)} views: "
                     f"{', '.join(e['name'] for e in sol_entries)})")
        except Exception as e:
            msh_export_error = f"{type(e).__name__}: {e}"
            progress("FEM", f"GMSH export failed: {msh_export_error}")

    result = {
        "status": "ok",
        "method": "FEM A-V (coil meshed + wp SIBC + Kelvin)",
        "frequency_hz": float(frequency),
        "current_A": float(I_target),
        "ndof": int(fes.ndof),
        "ne": int(mesh.ne),
        "fes_order": int(fes_order),
        "has_kelvin": bool(has_kelvin),
        "kelvin_periodic": bool(has_kelvin_periodic),
        # Current extraction diagnostic
        "I_out_pre_scale_abs": float(abs(I_out_pre)),
        "scale_factor_abs": float(abs(scale)),
        # Inductance
        "L_total_nH": L_total * 1e9,
        "L_vol_nH": L_vol * 1e9,
        "L_skin_wp_nH": L_skin_wp * 1e9,
        # Dissipation
        "P_total_W": P_total,
        "P_coil_W": P_coil,
        "P_wp_W": P_wp,
        "R_total_ohm": 2 * P_total / I_target ** 2,
        # Diagnostics
        "H_t_rms_wp_Am": H_t_rms_wp,
        "wp_area_m2": A_wp,
        "coil_delta_mm": coil_delta * 1e3,
        "coil_h_max_mm": coil_h_max * 1e3,
        "Z_s_wp_real": float(Z_s_wp.real),
        "Z_s_wp_imag": float(Z_s_wp.imag),
        "impedance_model": impedance_model,
        "esim_iterations": int(n_iter_done),
        "esim_converged": bool(esim_converged),
        "esim_history": esim_history,
        "esim_per_panel": bool(gf_Zs is not None),
        # Surface heat-flux stats (Phase B thermal input)
        "q_surf_max": q_surf_max,
        "q_surf_mean": q_surf_mean,
        "q_surf_p95": q_surf_p95,
        "P_total_check": P_total_check,
        # File outputs
        "msh_file": gmsh_file,
        "qsurf_sol": qsurf_sol_path,
        "msh_export_error": msh_export_error,
        # Timings
        "t_load_s": float(t_load),
        "t_assembly_s": float(t_asm),
        "t_solve_s": float(t_solve),
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
        description="FEM A-V coil + wp SIBC + Kelvin (gapped torus, "
                    "source/sink ports)")
    parser.add_argument("--vol", required=True, help=".vol file")
    parser.add_argument("--frequency", type=float, required=True)
    parser.add_argument("--current", type=float, default=1.0,
                        help="Port current I_target [A]")
    parser.add_argument("--coil-sigma", type=float, default=5.8e7)
    parser.add_argument("--sigma", type=float, required=True,
                        help="Workpiece conductivity [S/m]")
    parser.add_argument("--mu-r", type=float, default=1.0,
                        help="Workpiece relative permeability")
    parser.add_argument("--half-thickness", type=float, default=0.01)
    parser.add_argument("--fes-order", type=int, default=1)
    parser.add_argument("--solver", default="pardiso",
                        choices=["pardiso", "bddc", "iccg", "ams",
                                  "shifted_ams"])
    parser.add_argument("--impedance-model", default="sibc",
                        choices=["sibc", "esim"],
                        help="sibc: linear Dowell (production). "
                             "esim: nonlinear Karl iteration (requires "
                             "--bh-file; re-assembles A-V per outer iter; "
                             "converges in a few iters, verified 2026-06-02).")
    parser.add_argument("--bh-file", default="",
                        help="BH table (2-col H[A/m] B[T]); required "
                             "for --impedance-model esim.")
    parser.add_argument("--esim-max-iter", type=int, default=15,
                        help="ESIM Karl iteration max.")
    parser.add_argument("--esim-tol", type=float, default=1e-3,
                        help="ESIM Karl iteration relative tolerance "
                             "(|dZ_s| / |Z_s|).")
    parser.add_argument("--esim-relax", type=float, default=0.5,
                        help="Karl under-relaxation (0,1]; 1.0 = full "
                             "step, 0.5 = half-step (default, matches "
                             "calc_fem_kelvin).  Lower if Karl "
                             "oscillates near saturation.")
    parser.add_argument("--esim-per-panel", action="store_true",
                        help="Per-BND-DOF ESIM Karl (v4.48.0+): one ESIM "
                             "cell solve per H1 workpiece-boundary DOF "
                             "using a per-DOF |H_t| from the FEM "
                             "solution.  Resolves spatial saturation "
                             "patterns; ~N_bnd_dofs cell-solve cost "
                             "per Karl iter.")
    parser.add_argument("--sibc-bnd", default="sibc")
    parser.add_argument("--source-bnd", default="source")
    parser.add_argument("--sink-bnd", default="sink")
    parser.add_argument("--coil-mat", default="coil")
    parser.add_argument("--msh-output", default="",
                        help="Optional GMSH .msh output path. When set, "
                             "the .vol mesh is converted to .msh after "
                             "solve so the application result action can "
                             "view the geometry.")
    parser.add_argument("--require-kelvin", action="store_true",
                        help="Fail-fast if the .vol does not have a "
                             "'kelvin' material with periodic "
                             "identification. The IH panel's FEM-full "
                             "method sets this; without it a typo'd "
                             "kelvin block name would silently demote "
                             "the open BC to reg-only gauge truncation. "
                             "Per CLAUDE.md 'No Fallbacks'.")
    return parser


def main():
    parser = build_argparser()

    def run(args):
        if args.require_kelvin:
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
                        f"(materials = {mats}). The IH panel's FEM-full "
                        f"method requires a Kelvin-extended mesh: re-run "
                        f"the Cubit .jou with add_kelvin enabled."}
            if probe_mesh.ngmesh.GetNrIdentifications() == 0:
                return {"error":
                        f"--require-kelvin: .vol has 'kelvin' material but "
                        f"no periodic identifications. Ensure the .jou "
                        f"sidesets kelvin_int / kelvin_ext are paired via "
                        f"'copy mesh surface' in the Cubit script."}
        if args.impedance_model == "esim" and not args.bh_file:
            return {
                "error": "--impedance-model esim requires --bh-file "
                         "with a 2-column BH curve (H[A/m], B[T])."
            }
        if args.solver == "shifted_ams":
            return {
                "error": "shifted_ams solver is not yet wired into "
                         "calc_fem_coilmesh. Use pardiso (direct) or "
                         "bddc (iterative p>=2) for now."
            }
        result = solve_fem_coilmesh(
            vol=args.vol,
            frequency=args.frequency,
            I_target=args.current,
            coil_sigma=args.coil_sigma,
            wp_sigma=args.sigma,
            wp_mu_r=args.mu_r,
            half_thickness=args.half_thickness,
            fes_order=args.fes_order,
            solver=args.solver,
            sibc_bnd=args.sibc_bnd,
            source_bnd=args.source_bnd,
            sink_bnd=args.sink_bnd,
            coil_mat=args.coil_mat,
            msh_output=args.msh_output,
            impedance_model=args.impedance_model,
            bh_file=args.bh_file,
            esim_max_iter=args.esim_max_iter,
            esim_tol=args.esim_tol,
            esim_relax=args.esim_relax,
            esim_per_panel=args.esim_per_panel,
        )
        return result

    calc_main(run, parser)


if __name__ == "__main__":
    main()
