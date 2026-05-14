"""A1_AC_sweep.py — Direct NGSolve eddy-current AC sweep on A1 cuboid.

Solves the eddy-current frequency-domain problem at multiple ω, extracts the
induced magnetic moment m_z(jω) per unit applied B0_z, and exposes that as
the "ground truth" Y(jω) for cross-checking the Cauer + SIBC unified
prediction emitted by A1_sibc_unified_compare.py.

Formulation
-----------
For uniform applied B0_z and applied vector potential A_ext = ½B0(-y, x, 0),
we solve in conductor-only HCurl complex space:

    ∫ (1/μ) curl(A) · curl(W) + jω σ A · W  dV  =  -jω σ ∫ A_ext · W dV

(i.e., A is the **reduced** complex potential such that the total potential
inside the conductor is A_ext + A.)  The eddy current density at AC steady
state is J = -jωσ(A_ext + A).  The induced magnetic dipole moment is

    m_z(jω) = ½ ∫_V (x J_y - y J_x) dV.

For uniform B excitation, m_z(jω)/B0 is the {\em complex magnetic
polarizability} of the conductor — at DC it equals zero (no current for
static B) but at finite ω it tracks the same Foster spectrum used in the
CLN admittance Y(s) up to an overall geometric prefactor.  Normalizing both
the AC ground-truth and the CLN prediction by their maximum |m_z| / |Y|
gives shape-independent comparison.

Usage
-----
    python A1_AC_sweep.py --order 3 --h-cond-mm 0.5 \
        --f-min 10 --f-max 1e7 --n-freq 15 --out-tag A1_AC_sweep
"""
from __future__ import annotations

import argparse
import json
import math
import sys
import time
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8", errors="replace", line_buffering=True)

from netgen.occ import Box, Pnt, OCCGeometry
from ngsolve import (
    Mesh, HCurl, BilinearForm, LinearForm, GridFunction,
    CoefficientFunction, curl, dx, x, y, z,
    Integrate, TaskManager, ngsglobals,
)

mu0 = 4 * math.pi * 1e-7
sigma_Cu = 5.8e7
B0 = 1.0   # Tesla
Lx, Ly, Lz = 17.72e-3, 17.72e-3, 2.0e-3   # A1 geometry


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--f-min", type=float, default=10.0)
    ap.add_argument("--f-max", type=float, default=1e7)
    ap.add_argument("--n-freq", type=int, default=15)
    ap.add_argument("--order", type=int, default=3)
    ap.add_argument("--h-cond-mm", type=float, default=0.5)
    ap.add_argument("--bonus-int", type=int, default=8)
    ap.add_argument("--out-tag", type=str, default="A1_AC_sweep")
    args = ap.parse_args()

    ngsglobals.msg_level = 0
    H_COND = args.h_cond_mm * 1e-3

    print("=" * 76)
    print(" A1 cuboid eddy-current AC sweep (NGSolve, complex HCurl)")
    print("=" * 76)
    print(f"  geometry = {Lx*1000}×{Ly*1000}×{Lz*1000} mm Cu")
    print(f"  σ = {sigma_Cu:.3g} S/m, B0 = {B0} T")
    print(f"  HCurl order = {args.order}, h_cond = {args.h_cond_mm} mm")
    print(f"  f range = [{args.f_min:g}, {args.f_max:g}] Hz, "
          f"n_freq = {args.n_freq}")
    print()

    # Build conductor-only geometry / mesh
    cond = Box(Pnt(-Lx/2, -Ly/2, -Lz/2), Pnt(Lx/2, Ly/2, Lz/2))
    cond.faces.name = "conductorBND"
    cond.mat("conductor")
    cond.maxh = H_COND
    geo = OCCGeometry(cond)
    t0 = time.time()
    with TaskManager():
        ngmesh = geo.GenerateMesh(maxh=H_COND, grading=0.5)
    mesh = Mesh(ngmesh)
    mesh.Curve(args.order)
    print(f"  Mesh: ne={mesh.ne}, generation+curve {time.time()-t0:.1f} s")

    # Complex HCurl with A_red·t = 0 on conductor boundary, equivalent to
    # A_total = A_ext outside the conductor (eddy current confined inside).
    fes = HCurl(mesh, order=args.order, complex=True, nograds=True,
                dirichlet="conductorBND")
    print(f"  fes.ndof = {fes.ndof} (free {sum(fes.FreeDofs())})")
    A_tr, A_te = fes.TnT()

    # External source potential: A_ext = ½ B0_z × r = ½B0 (-y, x, 0)
    A_ext = CoefficientFunction((-y, x, 0)) * (B0 / 2.0)

    # Sweep
    log_f_min = math.log10(args.f_min)
    log_f_max = math.log10(args.f_max)
    freqs = [10 ** (log_f_min + (log_f_max - log_f_min) * k / (args.n_freq - 1))
             for k in range(args.n_freq)]

    # also expose analytical |A_ext|² and m_z_naive for downstream extraction
    A_ext_norm2 = Integrate(A_ext * A_ext * dx(
        bonus_intorder=args.bonus_int), mesh)
    print(f"  |A_ext|^2_L2 (volume integral) = {float(A_ext_norm2):.6e}")
    print()
    print(f" {'f (Hz)':>11}  {'Re(m_z)':>14}  {'Im(m_z)':>14}  "
          f"{'Re<Ar,Ae>':>14}  {'Im<Ar,Ae>':>14}  {'t (s)':>9}")
    results = []
    for f in freqs:
        omega = 2 * math.pi * f
        # ∫ (1/μ)curl(A)·curl(W) + jωσ A·W = -jωσ ∫ A_ext · W
        a = BilinearForm(fes)
        a += (1.0/mu0) * curl(A_tr) * curl(A_te) * dx(
            bonus_intorder=args.bonus_int)
        a += 1j * omega * sigma_Cu * A_tr * A_te * dx(
            bonus_intorder=args.bonus_int)

        rhs = LinearForm(fes)
        rhs += (-1j) * omega * sigma_Cu * A_ext * A_te * dx(
            bonus_intorder=args.bonus_int)

        t_s = time.time()
        with TaskManager():
            a.Assemble()
            rhs.Assemble()
            try:
                inv = a.mat.Inverse(fes.FreeDofs(), inverse="pardiso")
            except Exception:
                inv = a.mat.Inverse(fes.FreeDofs(),
                                     inverse="sparsecholesky")
            gfA = GridFunction(fes, name=f"A_red_{f:.0e}")
            gfA.vec.data = inv * rhs.vec
        t_solve = time.time() - t_s

        # Eddy current J = -jωσ (A_ext + A_red), where A_red := gfA
        A_total = A_ext + gfA
        J = (-1j * omega * sigma_Cu) * A_total

        # Induced magnetic dipole m_z = ½ ∫(x J_y - y J_x) dV  (complex)
        m_z_integrand = (x * J[1] - y * J[0]) * 0.5
        m_z_val = Integrate(m_z_integrand * dx(
            bonus_intorder=args.bonus_int), mesh)

        m_re = float(m_z_val.real)
        m_im = float(m_z_val.imag)
        m_abs = math.hypot(m_re, m_im)
        m_phase_deg = math.degrees(math.atan2(m_im, m_re))

        # Cauer-Y-equivalent observable: ⟨A_red, A_ext⟩ (complex, L²) where
        # A_red = gfA is the induced (reduced) potential.
        # In a Foster decomposition with M-orthonormal basis e_n and
        # g_n = ⟨A_ext, e_n⟩_M, the AC solution satisfies
        #   c_n(jω) = -jωτ_n g_n / (1 + jω τ_n)
        # so  ⟨A_red, A_ext⟩ = Σ c_n g_n = -Σ g_n²·jωτ_n/(1+jωτ_n)
        #                    = -[Y(0) - Y(jω)]  =  Y(jω) - Y(0)
        # The complex shape of ⟨A_red, A_ext⟩(ω)+|A_ext|² therefore tracks
        # the CLN Y(jω) directly.
        ared_dot_aext = Integrate(gfA * A_ext * dx(
            bonus_intorder=args.bonus_int), mesh)
        ar_re = float(ared_dot_aext.real)
        ar_im = float(ared_dot_aext.imag)
        ar_abs = math.hypot(ar_re, ar_im)
        ar_phase_deg = math.degrees(math.atan2(ar_im, ar_re))

        results.append({
            "f_Hz": float(f),
            "omega_rad_per_s": float(omega),
            "m_z_real": m_re,
            "m_z_imag": m_im,
            "m_z_abs": m_abs,
            "m_z_phase_deg": m_phase_deg,
            "ared_dot_aext_real": ar_re,
            "ared_dot_aext_imag": ar_im,
            "ared_dot_aext_abs": ar_abs,
            "ared_dot_aext_phase_deg": ar_phase_deg,
            "t_solve_s": float(t_solve),
        })
        print(f"  {f:>9.3e}  {m_re:>14.6e}  {m_im:>14.6e}  "
              f"{ar_re:>14.6e}  {ar_im:>14.6e}  {t_solve:>9.2f}")

    # Save
    out = {
        "method": "Direct NGSolve eddy-current AC sweep (complex HCurl)",
        "geometry_mm": [Lx * 1000, Ly * 1000, Lz * 1000],
        "params": {
            "sigma_S_per_m": sigma_Cu, "B0_T": B0,
            "order": args.order, "h_cond_mm": args.h_cond_mm,
            "bonus_int": args.bonus_int,
            "ne": mesh.ne, "ndof": fes.ndof,
        },
        "A_ext_norm2_L2": float(A_ext_norm2),
        "results": results,
    }
    out_path = Path(__file__).parent / f"{args.out_tag}_results.json"
    out_path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print()
    print(f"Saved: {out_path}")


if __name__ == "__main__":
    main()
