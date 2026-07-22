"""Hollaus Effective Material homogenization for laminated iron cores.

Stage-2 CLI script (per Panel Design Workflow Policy) implementing
the Effective Material (EM) approach from Hanser-Schoebinger-Hollaus
2024 (IEEE TMAG 60(10):7402307).

Knowledge prerequisite: `radia_mcp.motor.hollaus_eddy_knowledge`.

### Algorithm summary

Two-stage workflow:

**Offline (run once per iron grade + frequency range):**
1. Build a 1D cell problem covering ONE lamination period:
     iron strip (thickness d_iron) + insulation (thickness d_ins).
2. Apply uniform tangential H_t at the top boundary.
3. Solve 1D nonlinear time-harmonic eddy current PDE:
     d/dz [(1/(jw mu(B^2))) dH_t/dz] = sigma jw H_t (in iron)
     dH_t/dz = 0                                    (in insulation)
4. Integrate eddy-current losses (P_ECL) and reactive power (P_RP)
   over the cell volume.
5. Define effective material:
     mu_eff(B_avg, f) = J0_avg / B_avg          (real part)
     sigma_eff(B_avg, f) = P_ECL / (V_cell * |H_t|^2)
6. Store (mu_eff, sigma_eff) as a 2D lookup table over (|B|, f).

**Online (global FE with homogenized material):**
1. Replace the laminated stator iron region with a SINGLE
   homogeneous volume.
2. In the global 2D/3D FE, use `nu_eff(|B|, f) = 1/(mu_0 mu_eff)`.
3. Post-process: integrate `sigma_eff |E|^2` over the homogenized
   region → total ECL.
4. JSON stdout: total ECL, RP, B distribution.

### Usage

Offline cell problem only:
    python calc_motor_lamination.py --mode cell \\
        --d-iron 0.35e-3 --d-ins 0.05e-3 \\
        --sigma 4e6 --mu-r-iron 5000 \\
        --B-list 0.5,1.0,1.5 --freq-list 50,500,5000 \\
        --output cell_em_table.json

Online global FE:
    python calc_motor_lamination.py --mode global \\
        --vol stator_yoke.vol --em-table cell_em_table.json \\
        --H-amplitude 1000 --freq 1000 \\
        --output global_results.json

Full workflow (cell + global in one command):
    python calc_motor_lamination.py --mode full \\
        --vol stator_yoke.vol \\
        --d-iron 0.35e-3 --d-ins 0.05e-3 \\
        --sigma 4e6 --mu-r-iron 5000 \\
        --H-amplitude 1000 --freq 1000 \\
        --output full_results.json
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys

import numpy as np

# Direct solver for the FE Inverse(); default "pardiso" (MKL, fastest) is
# overridden by --linear-solver in main().  Golden tests pass "sparsecholesky"
# (ngsolve built-in, no MKL) -- same solution, but immune to the MKL PARDISO
# threading-layer load that is flaky under the pytest subprocess.
_LINEAR_SOLVER = "pardiso"

# Panel-common boilerplate
_this_dir = os.path.dirname(os.path.abspath(__file__))
if _this_dir not in sys.path:
    sys.path.insert(0, _this_dir)

from calc_common import (MU_0, NU_0, setup_paths, progress, calc_main)


def _log(msg):
    progress("LAMINATION", msg)


def solve_cell_problem(d_iron=0.35e-3, d_ins=0.05e-3,
                       sigma=4e6, mu_r_iron=5000.0,
                       B_avg=1.0, freq=1000.0,
                       n_elements=100,
                       drive="meanB", J_s_iron=0.0):
    """1D cell problem for one lamination period (iron + insulation).

    Solves the time-harmonic eddy-current equation in the iron layer
    (linear here for the first cut; nonlinear extension straightforward
    by Picard iteration on mu(|B|)).

    Returns dict with:
      mu_eff_real, mu_eff_imag : complex effective permeability
      sigma_eff                 : effective conductivity
      ECL_density               : eddy current loss per volume (W/m^3)
      RP_density                : reactive power per volume (VAR/m^3)

    ### Formulation (v2)

    We use the A-formulation for the in-plane vector potential A_y(z)
    (the 2D scalar equivalent of the magnetic vector potential).  The
    cell extends from z = 0 (iron bottom) to z = d_total (insulation top).
    Boundary conditions:

      A(z=0)        = 0
      A(z=d_total)  = B_avg . d_total          (linear-in-z applied A)

    so that the **cell-averaged flux density** equals B_avg.

    The weak form in iron:
      int [ nu (dA/dz)^2 + jw sigma A v ] dx = 0     (no source)
    in insulation:
      int [ nu (dA/dz)^2 ] dx = 0

    With uniform B_avg imposed at top/bottom, the field penetrates the
    lamination and the eddy currents close in loops, reproducing the
    classical Bertotti scaling P_ecl ~ d^2 * sigma * B_avg^2 * omega^2
    in the thin-sheet regime (d << skin depth).

    ECL density (per cell volume) = (sigma * omega^2 / V_total) *
                                     int_{iron} |A|^2 dz

    derived from J = -jw sigma A.
    """
    setup_paths()
    from ngsolve import (Mesh, H1, BilinearForm, LinearForm, GridFunction,
                         CoefficientFunction, grad, x, y, dx, ds, Integrate, Conj, TaskManager)

    # Build a thin 2D strip mesh as a 1D-in-z proxy.
    from netgen.geom2d import SplineGeometry
    geo = SplineGeometry()
    d_total = d_iron + d_ins
    width = max(d_total, 1e-4)
    pts = [
        geo.AppendPoint(0, 0),
        geo.AppendPoint(width, 0),
        geo.AppendPoint(width, d_iron),
        geo.AppendPoint(width, d_total),
        geo.AppendPoint(0, d_total),
        geo.AppendPoint(0, d_iron),
    ]
    geo.Append(["line", pts[0], pts[1]], leftdomain=1, rightdomain=0, bc="bottom")
    geo.Append(["line", pts[1], pts[2]], leftdomain=1, rightdomain=0, bc="ironside")
    geo.Append(["line", pts[2], pts[5]], leftdomain=1, rightdomain=2, bc="interface")
    geo.Append(["line", pts[5], pts[0]], leftdomain=1, rightdomain=0, bc="ironside")
    geo.Append(["line", pts[2], pts[3]], leftdomain=2, rightdomain=0, bc="insside")
    geo.Append(["line", pts[3], pts[4]], leftdomain=2, rightdomain=0, bc="top")
    geo.Append(["line", pts[4], pts[5]], leftdomain=2, rightdomain=0, bc="insside")
    geo.SetMaterial(1, "iron")
    geo.SetMaterial(2, "insulation")
    maxh = d_iron / max(n_elements, 5)
    ngmesh = geo.GenerateMesh(maxh=maxh)
    mesh = Mesh(ngmesh)

    omega = 2.0 * math.pi * freq
    nu_iron_lin = 1.0 / (mu_r_iron * MU_0)
    nu_ins = 1.0 / MU_0

    # Dispatch on drive mode (Hanser-Schobinger-Hollaus 2025 extension):
    #   meanB    -- impose uniform mean B via Dirichlet on A (v2.1, default)
    #   current  -- impose uniform J_s in iron, A free at boundaries
    #               (Neumann, symmetric strip → measures iron's intrinsic
    #               response to a current-source drive)
    #   voltage  -- impose terminal voltage (= jω flux), Lagrange multiplier
    #               binds the integrated A difference (not yet implemented
    #               -- raise)
    if drive not in ("meanB", "current", "voltage"):
        raise ValueError(f"unknown drive mode: {drive!r}")
    if drive == "voltage":
        raise NotImplementedError(
            "drive='voltage' requires Lagrange multiplier coupling -- "
            "planned for Hanser 2025 follow-up.  Use drive='current' "
            "or drive='meanB' for now."
        )

    if drive == "meanB":
        dirichlet_bnd = "top|bottom"
    else:  # current
        dirichlet_bnd = ""  # no Dirichlet -- free boundary
    fes = H1(mesh, order=2, complex=True, dirichlet=dirichlet_bnd)
    A_sol = GridFunction(fes)

    if drive == "meanB":
        # Dirichlet: A(y=0)=0, A(y=d_total)=B_avg*d_total
        A_top_val = complex(B_avg * d_total, 0.0)
        A_sol.Set(CoefficientFunction(0.0),
                  definedon=mesh.Boundaries("bottom"))
        A_sol.Set(CoefficientFunction(A_top_val),
                  definedon=mesh.Boundaries("top"))

    u, v = fes.TrialFunction(), fes.TestFunction()
    # Iron: nu (dA/dy)^2 + jw sigma A v  (sigma > 0)
    # Insulation: nu (dA/dy)^2 only
    a = BilinearForm(fes, symmetric=False)
    a += (nu_iron_lin * grad(u)[1] * grad(v)[1]
          + 1j * omega * sigma * u * v) * dx(definedon=mesh.Materials("iron"))
    a += (nu_ins * grad(u)[1] * grad(v)[1]) * dx(definedon=mesh.Materials("insulation"))
    # Penalize the constant mode in current-drive mode (no Dirichlet → 1
    # null vector).  A tiny mass term anchors the gauge.
    if drive == "current":
        a += 1e-12 * u * v * dx
    with TaskManager():
        a.Assemble()

        f = LinearForm(fes)
        if drive == "current":
            # Apply uniform J_s in iron: integrand = J_s * v
            f += complex(J_s_iron, 0.0) * v * dx(definedon=mesh.Materials("iron"))
        f.Assemble()

        inv = a.mat.Inverse(fes.FreeDofs(), inverse=_LINEAR_SOLVER)
        res = f.vec - a.mat * A_sol.vec
        A_sol.vec.data = A_sol.vec + inv * res

        V_iron = Integrate(CoefficientFunction(1.0), mesh,
                            definedon=mesh.Materials("iron"))
        V_total = Integrate(CoefficientFunction(1.0), mesh)

        # GAUGE FIX (v2.1): subtract the iron-region mean of A so that
        # <A>_iron = 0.  This corresponds to placing the zero-EMF reference
        # at the sheet center, matching the classical Bertotti derivation
        # (where the induced E field circulates around the sheet midplane).
        # Without this gauge fix, ECL gets a constant offset that makes the
        # absolute magnitudes wrong even though the d-dependence is correct.
        A_mean_iron_complex = (
            Integrate(A_sol, mesh, definedon=mesh.Materials("iron"))
            / float(V_iron)
        )

        # Build a gauge-shifted CF: A_shifted = A_sol - <A>_iron (complex constant)
        A_shifted = A_sol - CoefficientFunction(A_mean_iron_complex)
        A_shifted_conj = Conj(A_shifted)
        A_sq = (A_shifted * A_shifted_conj).real

        grad_A = grad(A_sol)
        Bx_real = -grad_A[1].real
        Bx_imag = -grad_A[1].imag
        B_sq = Bx_real * Bx_real + Bx_imag * Bx_imag

        # ECL = (1/2) * sigma * omega^2 * integral |A_shifted|^2 in iron
        int_A_sq_iron = Integrate(A_sq, mesh,
                                   definedon=mesh.Materials("iron")).real
        ECL_iron = 0.5 * sigma * omega * omega * int_A_sq_iron
        ECL_density = float(ECL_iron / V_total)

        # Reactive power: omega * (1/(2mu)) * integral |B|^2 dV averaged
        int_B_sq_iron = Integrate(B_sq, mesh,
                                   definedon=mesh.Materials("iron")).real
        mu_iron = mu_r_iron * MU_0
        RP_iron = 0.5 * omega * (1.0 / mu_iron) * int_B_sq_iron
        int_B_sq_ins = Integrate(B_sq, mesh,
                                  definedon=mesh.Materials("insulation")).real
        RP_ins = 0.5 * omega * NU_0 * int_B_sq_ins
        RP_density = float((RP_iron + RP_ins) / V_total)

        # Effective (complex) permeability of the LINEAR cell -- the in-plane
        # lamination homogenization
        #   mu_eff = mu0 [ fill*mu_r*tanh(b)/b + (1-fill) ],
        #   b = (d_iron/2) sqrt(j w mu0 mu_r sigma).
        # This is the canonical radia_mcp.radia_ngsolve.solve.laminated_mu_eff
        # (validated < 2.2 % vs 1D FE from |b|=0.5 to ~10 --
        # validation_test/radia_mcp/test_laminated_mu_eff_freqsweep.py).  It is INLINED here rather
        # than imported, because the GUI panel runs against the deployed radia_mcp
        # which may predate that helper.  It replaces the old volume-average, which
        # ignored the high-frequency skin reduction of the REAL part (mu' drops as
        # the field is expelled from the sheet).
        #
        # The eddy LOSS stays the separate, FE-computed ECL above (resistive,
        # 0.5 sigma |E|^2).  Do NOT derive the loss from Im(mu_eff): mu_eff is the
        # flux-averaging <B>/H_s permeability, and 0.5 w B^2 Im(1/mu_eff)
        # under-predicts the true eddy loss by ~70 % at |b|~5 (keep RP/mu_eff and
        # ECL as independent outputs).  For a future NONLINEAR cell, replace this
        # closed form with mu_eff = <B>/H_s read from the FE on a symmetric cell.
        import cmath
        fill_eff = float(V_iron / V_total)
        if omega > 0.0 and sigma > 0.0:
            _b = (d_iron / 2.0) * cmath.sqrt(1j * omega * mu_iron * sigma)
            _mu_eff_c = MU_0 * (fill_eff * mu_r_iron * (cmath.tanh(_b) / _b)
                                + (1.0 - fill_eff))
        else:
            _mu_eff_c = MU_0 * (fill_eff * mu_r_iron + (1.0 - fill_eff))
        mu_eff_real = float(_mu_eff_c.real)
        mu_eff_imag = float(_mu_eff_c.imag)
        sigma_eff = float((V_iron / V_total) * sigma)

        return {
            "B_avg": B_avg,
            "freq": freq,
            "mu_eff_real": mu_eff_real,
            "mu_eff_imag": mu_eff_imag,   # complex effective-mu (skin + loss-equiv)
            "sigma_eff": sigma_eff,
            "ECL_density_W_per_m3": ECL_density,
            "RP_density_VAR_per_m3": RP_density,
            "V_iron_fraction": float(V_iron / V_total),
            "n_dof": fes.ndof,
        }


def build_em_table(d_iron, d_ins, sigma, mu_r_iron,
                   B_list, freq_list, n_elements=100,
                   drive="meanB", J_s_iron=0.0):
    """Build a (|B|, freq) -> (mu_eff, sigma_eff, ECL, RP) lookup table.

    Args:
        drive: "meanB" (canonical Hollaus, default) or "current"
               (Hanser 2025 J_s-driven extension).  "voltage" is
               planned but not yet implemented.
        J_s_iron: A/m^2, uniform current density in iron (drive='current').

    Returns list of dicts (one per (B, freq) pair) suitable for JSON.
    """
    table = []
    for B in B_list:
        for f in freq_list:
            _log(f"cell problem: B={B}T, f={f}Hz, drive={drive}")
            res = solve_cell_problem(
                d_iron=d_iron, d_ins=d_ins,
                sigma=sigma, mu_r_iron=mu_r_iron,
                B_avg=B, freq=f, n_elements=n_elements,
                drive=drive, J_s_iron=J_s_iron,
            )
            res["drive"] = drive
            table.append(res)
    return table


def solve_global_fe_with_em(vol_file, em_table, H_amplitude, freq,
                            fes_order=2, msh_output=""):
    """Online global FE using effective-material properties on
    'stator_iron' (or similar laminated) region.

    For this first-cut implementation: linear magnetostatic 2D FE with
    nu = 1/(mu_eff * mu_0).  Eddy current losses are post-processed
    using sigma_eff * |B|^2 / mu_eff.
    """
    setup_paths()
    from ngsolve import (Mesh, H1, BilinearForm, LinearForm, GridFunction,
                         CoefficientFunction, grad, x, y, dx, ds, Integrate,
                         InnerProduct)
    from ngsolve import TaskManager
    if not os.path.exists(vol_file):
        return {"error": f"vol_file not found: {vol_file}"}

    mesh = Mesh(vol_file)
    _log(f"Loaded {vol_file}: ne={mesh.ne}")
    _log(f"  materials: {sorted(set(mesh.GetMaterials()))}")

    # Interpolate em_table to (current H_amplitude, freq)
    # For first cut: just average the table entries (no interpolation)
    if not em_table:
        return {"error": "em_table is empty"}
    mu_eff_avg = sum(r["mu_eff_real"] for r in em_table) / len(em_table)
    sigma_eff_avg = sum(r["sigma_eff"] for r in em_table) / len(em_table)
    ECL_density_lookup = sum(r["ECL_density_W_per_m3"] for r in em_table) / len(em_table)

    # Global 2D linear magnetostatic
    fes = H1(mesh, order=fes_order, dirichlet="outer")
    A = GridFunction(fes)
    u, v = fes.TrialFunction(), fes.TestFunction()

    nu_dict = {}
    sigma_dict = {}
    for m in set(mesh.GetMaterials()):
        if m.startswith("stator_iron") or m.startswith("rotor_iron"):
            nu_dict[m] = 1.0 / (mu_eff_avg * MU_0)  # effective
            sigma_dict[m] = sigma_eff_avg
        else:
            nu_dict[m] = NU_0
            sigma_dict[m] = 0.0
    nu_cf = mesh.MaterialCF(nu_dict, default=NU_0)
    sigma_cf = mesh.MaterialCF(sigma_dict, default=0.0)

    a = BilinearForm(fes, symmetric=True)
    a += nu_cf * grad(u) * grad(v) * dx
    a.Assemble()

    # No real source for a passive test — apply Neumann H_t on outer or
    # uniform source J_s as a placeholder
    f = LinearForm(fes)
    f.Assemble()

    if A.vec.Norm() == 0:
        # Apply test source: uniform current density 1 A/m^2 in
        # any "coil_*" material if present, otherwise return
        # the zero-source consistency check
        coil_mats = [m for m in set(mesh.GetMaterials())
                     if m.startswith("coil") or m.startswith("stator_ind")]
        if coil_mats:
            for cm in coil_mats:
                f += 1.0 * v * dx(definedon=mesh.Materials(cm))
            f.Assemble()

    inv = a.mat.Inverse(fes.FreeDofs(), inverse=_LINEAR_SOLVER)
    A.vec.data = inv * f.vec

    # Post-process: average |B|^2 in stator_iron, multiply by ECL density
    gA = grad(A)
    B2 = gA[0]*gA[0] + gA[1]*gA[1]
    iron_volume = 0.0
    B2_avg = 0.0
    if "stator_iron" in set(mesh.GetMaterials()):
        iron_volume = Integrate(CoefficientFunction(1.0), mesh,
                                 definedon=mesh.Materials("stator_iron"))
        if iron_volume > 0:
            B2_avg = Integrate(B2, mesh,
                                definedon=mesh.Materials("stator_iron")) / iron_volume
    total_ECL = ECL_density_lookup * iron_volume

    gmsh_file = ""
    if msh_output:
        from radia.gmsh_post_export import GmshPostExport

        gmsh_file = os.path.abspath(msh_output)
        os.makedirs(os.path.dirname(gmsh_file), exist_ok=True)
        magnetic_flux_density = CoefficientFunction(
            (gA[1], -gA[0], 0.0)
        )
        with TaskManager():
            post = GmshPostExport(mesh)
            post.add_scalar_field("A_z_Wb_per_m", A)
            post.add_vector_field(
                "B_T", magnetic_flux_density, cell_data=True
            )
            post.write(gmsh_file)

    return {
        "vol_file": vol_file,
        "fes_order": fes_order,
        "ndof": fes.ndof,
        "n_em_table_entries": len(em_table),
        "mu_eff_used": mu_eff_avg,
        "sigma_eff_used": sigma_eff_avg,
        "ECL_density_lookup_W_per_m3": ECL_density_lookup,
        "iron_volume_m3": iron_volume,
        "B2_avg_T2": B2_avg,
        "total_ECL_W": total_ECL,
        "gmsh_file": gmsh_file,
        "method_paper": "Hanser-Schöbinger-Hollaus, IEEE TMAG 60(10) 2024",
        "knowledge": "radia_mcp.motor.hollaus_eddy_knowledge.effective_material",
    }


def build_argparser():
    parser = argparse.ArgumentParser(
        description="Hollaus Effective Material homogenization for laminated iron")
    parser.add_argument("--mode", choices=["cell", "global", "full"],
                        default="cell",
                        help="cell=offline cell problem only, "
                             "global=use existing em-table on .vol, "
                             "full=cell+global in one run")
    # Cell problem args
    parser.add_argument("--d-iron", type=float, default=0.35e-3,
                        help="m, lamination iron thickness")
    parser.add_argument("--d-ins", type=float, default=0.05e-3,
                        help="m, insulation layer thickness")
    parser.add_argument("--sigma", type=float, default=4e6,
                        help="S/m, iron conductivity")
    parser.add_argument("--mu-r-iron", type=float, default=5000.0,
                        help="relative permeability of iron (linear approx)")
    parser.add_argument("--B-list", default="0.5,1.0,1.5",
                        help="comma-separated B amplitudes (T)")
    parser.add_argument("--freq-list", default="50,500,5000",
                        help="comma-separated frequencies (Hz)")
    parser.add_argument("--cell-n-elements", type=int, default=100,
                        help="resolution of 1D cell problem")
    parser.add_argument("--drive", default="meanB",
                        choices=["meanB", "current", "voltage"],
                        help="cell driving: meanB (canonical Hollaus 2024), "
                             "current (Hanser 2025 J_s), voltage (planned)")
    parser.add_argument("--J-s-iron", type=float, default=0.0,
                        help="A/m^2, uniform current density in iron "
                             "(drive='current' only)")
    # Global FE args
    parser.add_argument("--vol", default="",
                        help="Netgen .vol mesh for global FE")
    parser.add_argument("--em-table", default="",
                        help="JSON file with em_table from prior cell run")
    parser.add_argument("--H-amplitude", type=float, default=1000.0,
                        help="A/m, peak applied H field")
    parser.add_argument("--freq", type=float, default=1000.0,
                        help="Hz, operating frequency for global FE")
    parser.add_argument("--fes-order", type=int, default=2,
                        help="HCurl order")
    parser.add_argument("--linear-solver", default="pardiso",
                        choices=["pardiso", "sparsecholesky", "umfpack"],
                        help="FE direct solver. pardiso=MKL (fastest, default); "
                             "sparsecholesky=ngsolve built-in (no MKL) -- use on "
                             "machines with a conflicting MKL on PATH (e.g. CST).")
    parser.add_argument(
        "--msh-output", default="",
        help="GMSH .msh v4.1 output for global/full FE modes",
    )
    return parser


def main():
    parser = build_argparser()

    def run(args):
        global _LINEAR_SOLVER
        _LINEAR_SOLVER = args.linear_solver
        B_list = [float(x) for x in args.B_list.split(",")]
        freq_list = [float(x) for x in args.freq_list.split(",")]

        if args.mode == "cell":
            if args.msh_output:
                raise ValueError(
                    "cell mode produces an effective-material table, not a spatial field"
                )
            table = build_em_table(
                args.d_iron, args.d_ins, args.sigma, args.mu_r_iron,
                B_list, freq_list, args.cell_n_elements,
                drive=args.drive, J_s_iron=args.J_s_iron,
            )
            return {"mode": "cell", "em_table": table,
                    "params": {"d_iron": args.d_iron, "d_ins": args.d_ins,
                                "sigma": args.sigma,
                                "mu_r_iron": args.mu_r_iron}}

        elif args.mode == "global":
            if not args.em_table or not os.path.exists(args.em_table):
                return {"error": "global mode requires --em-table FILE"}
            with open(args.em_table, "r", encoding="utf-8") as f:
                table_data = json.load(f)
            table = table_data.get("em_table", [])
            res = solve_global_fe_with_em(
                args.vol, table, args.H_amplitude, args.freq,
                args.fes_order, args.msh_output,
            )
            res["mode"] = "global"
            return res

        else:   # full
            _log("Step 1/2: building effective-material table from cell problem")
            table = build_em_table(
                args.d_iron, args.d_ins, args.sigma, args.mu_r_iron,
                B_list, freq_list, args.cell_n_elements,
                drive=args.drive, J_s_iron=args.J_s_iron,
            )
            _log("Step 2/2: solving global FE with effective material")
            res = solve_global_fe_with_em(
                args.vol, table, args.H_amplitude, args.freq,
                args.fes_order, args.msh_output,
            )
            res["mode"] = "full"
            res["em_table"] = table
            return res

    calc_main(run, parser)


if __name__ == "__main__":
    main()
