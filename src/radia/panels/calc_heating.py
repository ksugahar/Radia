"""
FEM-ESIM workpiece heating calculator (2D axisymmetric).

Called as subprocess from Cubit panel:
    python calc_heating.py --frequency 50000 --sigma 2e6 --material steel
        --r-coil 0.03 --a-coil 0.003 --r-wp 0.01 --h-wp 0.02

Stage 2 of induction heating design:
  Stage 1: BEM -> coil inductance L (quick, surface mesh only)
  Stage 2: FEM-ESIM -> workpiece heating P (this script)

Pipeline (independent from BEM):
  1. Auto-generate 2D axi mesh: air + coil + workpiece hole + Kelvin
  2. Solve static FEM (A-formulation, coil = J0 source)
  3. Sample H_t on workpiece surface from grad(u)
  4. ESIM cell problem per surface segment -> P', Q'
  5. Return P distribution as JSON

IMPORTANT: NGSolve must be imported BEFORE cubit.
Outputs JSON to stdout.
"""

import argparse
import json
import math
import os
import sys
import time
import numpy as np


MU_0 = 4e-7 * np.pi
NU_0 = 1.0 / MU_0

# Import canonical material definitions
_radia_src = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
if os.path.abspath(_radia_src) not in sys.path:
    sys.path.insert(0, os.path.abspath(_radia_src))
from em_material import STEEL_BH, EMMaterial, add_material_args


def build_axi_mesh(R_coil, a_coil, R_wp, H_wp, a_kelvin, z_offset, maxh):
    """Auto-generate 2D axisymmetric mesh: air + coil + workpiece hole + Kelvin.

    No Cubit needed. Pure OCC geometry.
    """
    from ngsolve import Mesh
    from netgen.occ import (WorkPlane, Axes, MoveTo, Pnt, Dir, Z, X,
                             OCCGeometry, Glue, Vertex, IdentificationType)

    # Interior half-circle
    inner = WorkPlane().Circle(a_kelvin).Face()
    cutter = MoveTo(-a_kelvin - 0.1, -a_kelvin - 0.1).Rectangle(
        a_kelvin + 0.1, 2 * a_kelvin + 0.2).Face()
    inner = inner - cutter

    # Coil cross-section
    w = 2 * a_coil
    coil = MoveTo(R_coil - a_coil, -a_coil).Rectangle(w, w).Face()
    coil.name = "coil"
    coil.maxh = maxh / 2

    # Workpiece hole (name edges before subtraction)
    wp = MoveTo(0, -H_wp / 2).Rectangle(R_wp, H_wp).Face()
    for e in wp.edges:
        e.name = "wp_sibc" if e.center.x > 1e-4 else "axis"

    air = inner - coil - wp
    air.name = "air_inner"

    # Kelvin exterior
    outer = WorkPlane(Axes((0, z_offset, 0), n=Z, h=X)).Circle(a_kelvin).Face()
    cutter_ext = MoveTo(-a_kelvin - 0.1, z_offset - a_kelvin - 0.1).Rectangle(
        a_kelvin + 0.1, 2 * a_kelvin + 0.2).Face()
    outer = outer - cutter_ext
    outer.name = "air_outer"
    gnd = Vertex(Pnt(0, z_offset, 0))
    gnd.name = "GND"

    # Edge naming
    ki, ke = [], []
    for e in air.edges:
        if hasattr(e, 'name') and e.name == 'wp_sibc':
            continue
        cx = e.center.x
        try:
            v0, v1 = e.vertices
            d0 = math.sqrt(v0.p.x**2 + v0.p.y**2)
            d1 = math.sqrt(v1.p.x**2 + v1.p.y**2)
            is_arc = abs(d0 - a_kelvin) < 0.01 and abs(d1 - a_kelvin) < 0.01 and cx > 0.01
        except:
            is_arc = False
        if is_arc:
            e.name = "kelvin_int"; ki.append(e)
        elif cx < 1e-4:
            e.name = "axis"
        else:
            e.name = "default"

    for e in outer.edges:
        cx = e.center.x
        try:
            v0, v1 = e.vertices
            d0 = math.sqrt(v0.p.x**2 + (v0.p.y - z_offset)**2)
            d1 = math.sqrt(v1.p.x**2 + (v1.p.y - z_offset)**2)
            is_arc = abs(d0 - a_kelvin) < 0.01 and abs(d1 - a_kelvin) < 0.01 and cx > 0.01
        except:
            is_arc = False
        if cx < 1e-4:
            e.name = "axis_ext"
        elif is_arc:
            e.name = "kelvin_ext"; ke.append(e)
        else:
            e.name = "default"

    for ie in ki:
        for ee in ke:
            if (ie.center.y > 0 and ee.center.y - z_offset > 0) or \
               (ie.center.y < 0 and ee.center.y - z_offset < 0):
                ie.Identify(ee, "kelvin", IdentificationType.PERIODIC)
                break

    shape = Glue([air, coil, outer, gnd])
    mesh = Mesh(OCCGeometry(shape, dim=2).GenerateMesh(maxh=maxh, grading=0.4))
    return mesh


def compute_heating(R_coil, a_coil, R_wp, H_wp,
                    mat=None, frequency=50000,
                    I_total=1.0, maxh=0.003, order=3,
                    esim_geometry='cylinder'):
    """Compute workpiece heating via FEM-ESIM.

    Returns dict with P_total, P_distribution, L_coil, diagnostics.
    """
    if mat is None:
        mat = EMMaterial.from_name("steel")
    sigma = mat.sigma

    from ngsolve import (H1, Periodic, BilinearForm, LinearForm, GridFunction,
                         Integrate, Conj, IfPos, sqrt, grad, dx, ds, x, y)

    radia_src = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
    if os.path.abspath(radia_src) not in sys.path:
        sys.path.insert(0, os.path.abspath(radia_src))
    from esim_cell_problem import ESIMFiniteSlabSolver

    omega = 2 * np.pi * frequency
    J0 = I_total / (2 * a_coil)**2
    rho = mat.rho
    delta = mat.skin_depth(frequency)

    bh_curve = mat.bh_curve
    mu_r = mat.mu_r if bh_curve is None else None

    # Kelvin parameters (auto-sized)
    a_kelvin = max(R_coil * 3, R_wp * 5, 0.15)
    z_offset = a_kelvin * 2.5

    sys.stderr.write("HEATING_MESH:building auto mesh\n")
    sys.stderr.flush()
    t0 = time.perf_counter()
    mesh = build_axi_mesh(R_coil, a_coil, R_wp, H_wp, a_kelvin, z_offset, maxh)
    t_mesh = time.perf_counter() - t0

    r_cf = IfPos(x - 1e-10, x, 1e-10)
    yl = y - z_offset
    rp = sqrt(x**2 + yl**2)
    rs = IfPos(rp - 1e-10, rp, 1e-10)
    kf = (rs / a_kelvin)**2

    mats = mesh.GetMaterials()
    nu_d = {m: NU_0 * kf if "outer" in m.lower() else NU_0 for m in mats}
    from ngsolve import CoefficientFunction
    nu_cf = mesh.MaterialCF(nu_d, default=NU_0)

    # Static solve.
    #
    # POLICY VIOLATION FLAG (2026-05-10): this script uses standard
    # `H1(mesh, order=p)` which CLAUDE.md "Axisymmetric FE: Henrotte
    # Basis Only" forbids for axisymmetric solves.  Migration to
    # `radia.radia_axifemm.H1Henrotte + AxiHenrotteStiffnessBFI` was
    # attempted in radia 4.33.0 development and **deferred** because
    # the Kelvin transformation here uses a spatially-varying nu_cf
    # (kf = (rs/a_kelvin)^2 in the air_outer region) and the existing
    # AxiHenrotteStiffnessBFI samples mu_cf only at the element
    # centroid.  Centroid sampling cannot represent the kf variation
    # within an element, and the migration produced numerically wrong
    # P_total (~3580x off) and L (~70x off) on the steel cylinder
    # benchmark.  Proper fix needs per-quadrature-point mu evaluation
    # in the BFI (loses closed-form benefit) OR rederiving the Kelvin
    # transform on the Henrotte (s = r^2, z) basis; both are larger
    # tasks than a calc-script migration.
    #
    # This script remains on standard H1 + 2 pi r weighting until
    # the Henrotte BFI gains spatially-varying mu support.
    sys.stderr.write("HEATING_SOLVE:FEM static solve (legacy H1 -- "
                      "Henrotte migration blocked on Kelvin nu_cf)\n")
    sys.stderr.flush()
    t0 = time.perf_counter()
    fes = Periodic(H1(mesh, order=order, complex=False,
                      dirichlet="axis|axis_ext", dirichlet_bbnd="GND"))
    u, v = fes.TnT()
    bf = BilinearForm(fes)
    bf += nu_cf / r_cf * grad(u) * grad(v) * dx
    lf = LinearForm(fes)
    lf += J0 * v * dx("coil")
    bf.Assemble()
    lf.Assemble()
    gfu = GridFunction(fes)
    gfu.vec.data = bf.mat.Inverse(fes.FreeDofs(), inverse="pardiso") * lf.vec
    t_solve = time.perf_counter() - t0

    # Inductance
    W = 2 * np.pi * float(Integrate(
        0.5 * nu_cf / r_cf * grad(gfu) * grad(gfu), mesh))
    L = 2 * W / I_total**2

    # --- ESIM on workpiece surface ---
    sys.stderr.write("HEATING_ESIM:computing surface impedance\n")
    sys.stderr.flush()
    t0 = time.perf_counter()

    esim_solver = ESIMFiniteSlabSolver(
        half_thickness=R_wp, bh_curve=bh_curve, sigma=sigma,
        frequency=frequency,
        mu_r=mu_r if bh_curve is None else None, n_nodes=200,
        geometry=esim_geometry)

    grad_u = grad(gfu)
    eps = maxh * 0.2
    n_lat, n_cap = 30, 12

    P_total = 0.0
    Q_total = 0.0
    panels = []  # per-panel results for distribution

    # Lateral surface
    for j in range(n_lat):
        z_j = -H_wp / 2 + (j + 0.5) * H_wp / n_lat
        dA = 2 * np.pi * R_wp * (H_wp / n_lat)
        try:
            mip = mesh(R_wp + eps, z_j)
            H_t = abs(float(grad_u[0](mip)) / (MU_0 * (R_wp + eps)))
        except:
            H_t = 0.0
        sol = esim_solver.solve(max(H_t, 1e-3))
        P_panel = sol['P_prime'] * dA
        P_total += P_panel
        Q_total += sol['Q_prime'] * dA
        panels.append({
            "type": "lateral", "z": round(z_j, 6),
            "H_t": round(H_t, 4),
            "P_density": round(sol['P_prime'], 6),
            "P": round(float(P_panel), 10),
        })

    # Top and bottom caps
    for sign, label in [(+1, "top"), (-1, "bottom")]:
        for k in range(n_cap):
            r_k = (k + 0.5) * R_wp / n_cap
            if r_k < 1e-6:
                continue
            dA = 2 * np.pi * r_k * (R_wp / n_cap)
            try:
                mip = mesh(r_k, sign * (H_wp / 2 + eps))
                H_t = abs(float(grad_u[1](mip)) / (MU_0 * r_k))
            except:
                H_t = 0.0
            sol = esim_solver.solve(max(H_t, 1e-3))
            P_panel = sol['P_prime'] * dA
            P_total += P_panel
            Q_total += sol['Q_prime'] * dA
            panels.append({
                "type": label, "r": round(r_k, 6),
                "H_t": round(H_t, 4),
                "P_density": round(sol['P_prime'], 6),
                "P": round(float(P_panel), 10),
            })

    t_esim = time.perf_counter() - t0

    return {
        "P_total": float(P_total),
        "Q_total": float(Q_total),
        "L_coil": float(L),
        "delta": float(delta),
        "n_panels": len(panels),
        "ndof": fes.ndof,
        "ne": mesh.ne,
        "t_mesh": round(t_mesh, 2),
        "t_solve": round(t_solve, 2),
        "t_esim": round(t_esim, 2),
        "frequency": frequency,
        "sigma": sigma,
        "material": mat.name,
        "panels": panels,
    }


def main():
    _panels_dir = os.path.dirname(os.path.abspath(__file__))
    if _panels_dir not in sys.path:
        sys.path.insert(0, _panels_dir)
    from calc_common import calc_main

    parser = argparse.ArgumentParser(
        description="FEM-ESIM workpiece heating (2D axisymmetric)")
    parser.add_argument("--r-coil", type=float, required=True, help="Coil radius [m]")
    parser.add_argument("--a-coil", type=float, required=True, help="Coil wire radius [m]")
    parser.add_argument("--r-wp", type=float, required=True, help="Workpiece radius [m]")
    parser.add_argument("--h-wp", type=float, required=True, help="Workpiece height [m]")
    parser.add_argument("--frequency", type=float, default=50000, help="Frequency [Hz]")
    add_material_args(parser, include_custom=False)
    parser.add_argument("--maxh", type=float, default=0.003, help="Mesh size [m]")
    parser.add_argument("--esim-geometry", default="local_curvature",
                        choices=["local_curvature", "none"],
                        help="ESIM curvature: local_curvature or none (flat)")
    parser.add_argument("--output", default="", help="Output JSON file")

    def run(args):
        _geom_map = {"local_curvature": "cylinder", "none": "slab"}
        esim_geom = _geom_map.get(args.esim_geometry, "cylinder")
        return compute_heating(
            args.r_coil, args.a_coil, args.r_wp, args.h_wp,
            mat=EMMaterial.from_args(args), frequency=args.frequency,
            maxh=args.maxh, esim_geometry=esim_geom)

    calc_main(run, parser)


if __name__ == "__main__":
    main()
