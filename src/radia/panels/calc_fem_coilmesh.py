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
                       coil_sigma, coil_mu_r,
                       wp_sigma, wp_mu_r, half_thickness,
                       fes_order=1, solver="pardiso",
                       sibc_bnd="sibc",
                       source_bnd="source", sink_bnd="sink",
                       coil_mat="coil"):
    """A-V formulation for volumetric coil + SIBC workpiece + Kelvin."""
    import radia  # noqa: F401  DLL path setup

    from ngsolve import (Mesh, HCurl, H1, Periodic, BilinearForm, LinearForm,
                         GridFunction, Integrate, Conj, curl, grad, dx, ds,
                         CF, BND, VOL, TaskManager, sqrt as ngsqrt, IfPos,
                         x, y, z, specialcf)
    from calc_fem_kelvin import detect_kelvin_offset

    omega = 2 * math.pi * frequency
    s = 1j * omega

    wp_mat = EMMaterial(name="wp", sigma=wp_sigma, mu_r=wp_mu_r)
    delta_wp = wp_mat.skin_depth(frequency)
    Z_s_wp = wp_mat.dowell_Zs(frequency, half_thickness)
    robin_wp = s / Z_s_wp

    coil_delta = math.sqrt(2.0 / (omega * coil_mu_r * MU_0 * coil_sigma))

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
            d = 0.0
            for i in range(len(pts)):
                for j in range(i + 1, len(pts)):
                    d = max(d, float(np.linalg.norm(pts[i] - pts[j])))
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
            nu_dict[m] = 1.0 / (coil_mu_r * MU_0) if m == coil_mat else NU_0
    nu_cf = mesh.MaterialCF(nu_dict, default=NU_0)

    # Compound FES: HCurl(A) x H1(phi on coil).
    dirichlet_A = "GND" if "GND" in boundaries else ""
    fesA_base = HCurl(mesh, order=fes_order, nograds=True, complex=True,
                      dirichlet=dirichlet_A)
    fesA = Periodic(fesA_base) if has_kelvin_periodic else fesA_base
    fesPhi = H1(mesh, order=fes_order, complex=True,
                definedon=mesh.Materials(coil_mat),
                dirichlet=f"{source_bnd}|{sink_bnd}")
    fes = fesA * fesPhi
    (A, phi), (N, psi) = fes.TnT()
    progress("FEM", f"ndof={fes.ndof} (A:{fesA.ndof} + phi:{fesPhi.ndof})")

    # Bilinear form (MCP pattern)
    a_bf = BilinearForm(fes, symmetric=True)
    a_bf += nu_cf * curl(A) * curl(N) * dx(bonus_intorder=4)
    # Tiny mass regularisation on non-Kelvin materials (HCurl gauge)
    non_kelvin = [m for m in materials if "kelvin" not in m.lower()]
    if non_kelvin:
        a_bf += 1e-6 * NU_0 * A * N * dx("|".join(non_kelvin))
    # Workpiece SIBC Robin
    a_bf += robin_wp * A.Trace() * N.Trace() * ds(sibc_bnd)
    # A-V compound eddy term on coil
    a_bf += s * coil_sigma * (A + grad(phi)) * (N + grad(psi)) * dx(coil_mat)

    progress("FEM", "assemble A-V")
    t0 = time.perf_counter()
    with TaskManager():
        a_bf.Assemble()
    t_asm = time.perf_counter() - t0
    progress("FEM", f"assembled ({t_asm:.1f}s)")

    # Dirichlet lift: phi=1 on source, phi=0 on sink
    gfu = GridFunction(fes)
    gf_A, gf_phi = gfu.components
    gf_phi.Set(CF(1), definedon=mesh.Boundaries(source_bnd))

    progress("FEM", f"solve ({solver})")
    t0 = time.perf_counter()
    with TaskManager():
        r = gfu.vec.CreateVector()
        r.data = -a_bf.mat * gfu.vec
        gfu.vec.data += a_bf.mat.Inverse(fes.FreeDofs(),
                                           inverse=solver) * r
    t_solve = time.perf_counter() - t0
    progress("FEM", f"solved ({t_solve:.1f}s)")

    # Volume-integral current extraction (Gauss-consistent).
    # I_out = int_coil J . grad(psi_n) dV, with psi_n scalar H1
    # Dirichlet 1 on source, 0 on sink.  Equals port flux int J . n dS.
    J_coil = -s * coil_sigma * (gf_A + grad(gf_phi))
    fes_psi_n = H1(mesh, order=1, complex=True,
                    definedon=mesh.Materials(coil_mat),
                    dirichlet=sink_bnd)
    gf_psi_n = GridFunction(fes_psi_n)
    gf_psi_n.Set(CF(1), definedon=mesh.Boundaries(source_bnd))
    I_out_pre = complex(Integrate(
        J_coil * grad(gf_psi_n), mesh,
        definedon=mesh.Materials(coil_mat)))
    progress("FEM", f"I_out (pre-scale) = {abs(I_out_pre):.4e}")

    if abs(I_out_pre) < 1e-20:
        raise RuntimeError("FEM I_out is zero; A-V setup failed.")

    scale = complex(I_target) / I_out_pre
    gfu.vec.data = complex(scale) * gfu.vec
    # Recompute J_coil expression after scaling (CF re-evaluates).
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

    return {
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
        # Timings
        "t_load_s": float(t_load),
        "t_assembly_s": float(t_asm),
        "t_solve_s": float(t_solve),
    }


def main():
    parser = argparse.ArgumentParser(
        description="FEM A-V coil + wp SIBC + Kelvin (gapped torus, "
                    "source/sink ports)")
    parser.add_argument("--vol", required=True, help=".vol file")
    parser.add_argument("--frequency", type=float, required=True)
    parser.add_argument("--current", type=float, default=1.0,
                        help="Port current I_target [A]")
    parser.add_argument("--coil-sigma", type=float, default=5.8e7)
    parser.add_argument("--coil-mu-r", type=float, default=1.0)
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
                             "esim: Karl iteration (WIP, raises).")
    parser.add_argument("--bh-file", default="",
                        help="BH table for ESIM (WIP).")
    parser.add_argument("--esim-max-iter", type=int, default=15,
                        help="ESIM Karl iteration max (WIP).")
    parser.add_argument("--esim-tol", type=float, default=1e-3,
                        help="ESIM Karl iteration tol (WIP).")
    parser.add_argument("--sibc-bnd", default="sibc")
    parser.add_argument("--source-bnd", default="source")
    parser.add_argument("--sink-bnd", default="sink")
    parser.add_argument("--coil-mat", default="coil")

    def run(args):
        if args.impedance_model == "esim":
            return {
                "error": "ESIM is not implemented in calc_fem_coilmesh "
                         "yet. Use --impedance-model sibc (linear SIBC "
                         "Dowell) for now. Karl iteration + BH-curve "
                         "ESIM coupling is WIP."
            }
        if args.solver == "shifted_ams":
            return {
                "error": "shifted_ams solver is not yet wired into "
                         "calc_fem_coilmesh. Use pardiso (direct) or "
                         "bddc (iterative p>=2) for now."
            }
        return solve_fem_coilmesh(
            vol=args.vol,
            frequency=args.frequency,
            I_target=args.current,
            coil_sigma=args.coil_sigma,
            coil_mu_r=args.coil_mu_r,
            wp_sigma=args.sigma,
            wp_mu_r=args.mu_r,
            half_thickness=args.half_thickness,
            fes_order=args.fes_order,
            solver=args.solver,
            sibc_bnd=args.sibc_bnd,
            source_bnd=args.source_bnd,
            sink_bnd=args.sink_bnd,
            coil_mat=args.coil_mat,
        )

    calc_main(run, parser)


if __name__ == "__main__":
    main()
