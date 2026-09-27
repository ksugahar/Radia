"""Temperature-dependent induction heating of a cylinder through the Curie band.

Route A (production, radia.ih_axisym_coupled): axisymmetric volumetric
eddy-current FE for A_phi on an (r, z) mesh of workpiece + coil + air,
re-solved from the element temperatures with sigma(T) and mu_r(T); the
Joule density heats the workpiece volume (enthalpy Newton integrator).
It is run on a ladder of time steps and EM-update intervals to show
convergence of the staggered coupling.

Route B (surface model): the wall source of route A at the initial
temperature -- q(z) and |H_t|(z) -- written on a 3D surface mesh and fed to
calc_heat_axisym, (i) as a fixed source and (ii) through the frozen-|H_t|
surface-impedance table (--em-table --allow-frozen-ht) built from the same
sigma(T), mu_r(T) (planar SIBC).  The comparison documents why the frozen
|H_t| model is refused by default: for a current-driven coil around a
ferromagnetic part |H_t| is not frozen.

Material curves are illustrative (a carbon-steel-like shape), not alloy
data.  Output: results/coupled_curie_cylinder_frozen_ht.json.

    python validation_test/induction_heating/coupled_curie_cylinder.py
"""
from __future__ import annotations

import argparse
import json
import math
import os
import platform
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, os.path.join(ROOT, "src", "radia"))
sys.path.insert(0, os.path.join(ROOT, "src", "radia", "panels"))

import ih_axisym_coupled as C  # noqa: E402  (no native Radia needed)

MU0 = C.MU0
R_WP, H_WP = 0.020, 0.040
COIL_R0, COIL_R1, COIL_H = 0.026, 0.030, 0.012
R_AIR = 0.15
RHO, CP, K = 7800.0, 600.0, 35.0
T0 = 20.0
T_TAB = [0.0, 400.0, 700.0, 770.0, 800.0, 1200.0, 2000.0]
SIGMA_TAB = [6.0e6, 2.4e6, 1.3e6, 1.15e6, 1.0e6, 0.85e6, 0.8e6]
MU_TAB = [100.0, 90.0, 80.0, 1.0, 1.0, 1.0, 1.0]


def em_material():
    return C.EMMaterialTable(T=T_TAB, sigma=SIGMA_TAB, mu_r=MU_TAB,
                             source="illustrative carbon-steel-like")


def build_mesh(maxh_skin=6e-5, maxh_wp=1.5e-3, maxh_air=0.01, skin=1.5e-3):
    from netgen.occ import Glue, MoveTo, OCCGeometry, X, Y
    from ngsolve import Mesh

    core = MoveTo(0, -H_WP / 2).Rectangle(R_WP - skin, H_WP).Face()
    core.faces.name = "wp"
    core.edges.Min(X).name = "axis"
    core.maxh = maxh_wp
    shell = MoveTo(R_WP - skin, -H_WP / 2).Rectangle(skin, H_WP).Face()
    shell.faces.name = "wp"
    shell.maxh = maxh_skin
    coil = MoveTo(COIL_R0, -COIL_H / 2).Rectangle(COIL_R1 - COIL_R0,
                                                   COIL_H).Face()
    coil.faces.name = "coil"
    coil.maxh = 1.0e-3
    air = MoveTo(0, -R_AIR).Rectangle(R_AIR, 2 * R_AIR).Face()
    air.faces.name = "air"
    air.edges.Max(X).name = "outer"
    air.edges.Min(X).name = "axis"
    air.edges.Max(Y).name = "top"
    air.edges.Min(Y).name = "bot"
    air = air - core - shell - coil
    return Mesh(OCCGeometry(Glue([air, core, shell, coil]), dim=2)
                .GenerateMesh(maxh=maxh_air))


def wall_source(em, zs):
    """q(z) [W/m^2] into the side wall and |H_t|(z) [A/m] just outside."""
    from ngsolve import IfPos, Integrate, x as r, y
    Q = em.heat_density()
    out_q, out_h = [], []
    dz = 0.5e-3
    for zi in zs:
        band = IfPos(dz / 2 - (y - zi), 1.0, 0.0) * \
            IfPos(dz / 2 + (y - zi), 1.0, 0.0)
        P = float(Integrate(Q * band * 2 * math.pi * r, em.mesh,
                            definedon=em.mesh.Materials("wp"), order=8).real)
        out_q.append(P / (2 * math.pi * R_WP * dz))
        e = 2e-5
        p1, p2 = R_WP + e, R_WP + 2 * e
        A1 = complex(em.gfA(em.mesh(p1, zi)))
        A2 = complex(em.gfA(em.mesh(p2, zi)))
        Bz = ((p2 * A2 - p1 * A1) / (p2 - p1)) / p1
        out_h.append(abs(Bz) / MU0)
    return np.asarray(out_q), np.asarray(out_h)


def route_a(mesh, frequency, current, dt, t_end, em_every):
    import ih_heat_transient as iht
    import ih_thermal_material as itm
    th = itm.ThermalMaterial.constant(RHO, CP, K)
    t0 = time.time()
    result, gfT, em = C.run_coupled(
        mesh, frequency=frequency, workpiece="wp", coils={"coil": current},
        dirichlet="axis|outer|top|bot", em_material=em_material(),
        thermal_material=th, boundaries=iht.HeatBoundaryTerms(), dt=dt,
        t_end=t_end, t_initial=T0, em_every=em_every)
    return {"dt_s": dt, "em_every": em_every,
            "runtime_s": time.time() - t0,
            "t_s": [h["t_s"] for h in result["history"]],
            "P_W": [h["P_W"] for h in result["history"]],
            "T_max_C": [h["T_max_nodal_C"] for h in result["history"]],
            "T_surface_mid_final_C": float(gfT(mesh(R_WP - 1e-6, 0.0))),
            "energy_balance_relative_error":
                result["nonlinear_transient"]["energy_balance_relative_error"],
            "em_power_balance_max":
                max(abs(h["power_balance_relative_error"])
                    for h in result["em_history"]),
            "skin_depth_final_max_m":
                result["em_history"][-1]["skin_depth_max_m"]}


def impedance_table(path, frequency):
    omega = 2 * math.pi * frequency
    H = np.logspace(2, 6.5, 40)
    T = np.linspace(0.0, 1300.0, 261)
    sig = np.interp(T, T_TAB, SIGMA_TAB)
    mu = np.interp(T, T_TAB, MU_TAB)
    reZ = np.sqrt(omega * MU0 * mu / (2 * sig))
    q = 0.5 * reZ[None, :] * H[:, None] ** 2
    np.savez(path, H_grid=H, T_grid=T, Zs_re=np.repeat(reZ[None, :], 40, 0),
             Zs_im=np.repeat(reZ[None, :], 40, 0), q_surf=q,
             meta=np.asarray(json.dumps({"frequency": frequency,
                                         "material": "illustrative-steel",
                                         "model": "linear planar SIBC"})))


def route_b(workdir, zs, q0, h0, frequency, dt, t_end):
    import calc_heat_axisym
    import ih_thermal
    from netgen.geom2d import SplineGeometry
    from netgen.occ import Axes, Cylinder, OCCGeometry, Pnt, Z
    from ngsolve import GridFunction, H1, Mesh

    solid = Cylinder(Axes(Pnt(0, 0, -H_WP / 2), Z), r=R_WP, h=H_WP)
    solid.faces.name = "side"
    solid.faces.Max(Z).name = "top"
    solid.faces.Min(Z).name = "bottom"
    vol = os.path.join(workdir, "em_cylinder.vol")
    OCCGeometry(solid).GenerateMesh(maxh=0.001).Save(vol)
    em = Mesh(vol)
    pts = ih_thermal.mesh_vertices(em)
    side = ih_thermal.boundary_vertex_numbers(em, ["side"])
    q = GridFunction(H1(em, order=1))
    h = GridFunction(H1(em, order=1))
    q.vec[:] = 0.0
    h.vec[:] = 0.0
    q.vec.FV().NumPy()[side] = np.interp(pts[side, 2], zs, q0)
    h.vec.FV().NumPy()[side] = np.interp(pts[side, 2], zs, h0)
    qs = os.path.join(workdir, "q.sol")
    hs = os.path.join(workdir, "q_Ht.sol")
    q.Save(qs)
    h.Save(hs)
    for sol, quantity, unit in ((qs, ih_thermal.QSURF_QUANTITY,
                                 ih_thermal.QSURF_UNIT),
                                (hs, ih_thermal.HT_QUANTITY,
                                 ih_thermal.HT_UNIT)):
        ih_thermal.write_field_sidecar(
            sol, mesh_path=vol, mesh=em, fes_order=1, quantity=quantity,
            unit=unit, boundaries=["side"],
            extra={"frequency_Hz": frequency, "producer": "route A wall"})
    table = os.path.join(workdir, "sibc_table.npz")
    impedance_table(table, frequency)
    geo = SplineGeometry()
    geo.AddRectangle((0, -H_WP / 2), (R_WP, H_WP / 2),
                     bcs=("bottom", "side", "top", "axis"))
    m2 = Mesh(geo.GenerateMesh(maxh=0.0005))
    common = dict(material="custom", rho=RHO, cp=CP, k=K, h_conv=0.0,
                  heat_flux_boundaries="side", qsurf_sol=qs, em_vol=vol,
                  dt=dt, t_end=t_end, t_initial=T0, fes_order=2,
                  _wp_mesh=m2, _write_solution=False, n_phi_samples=16,
                  power_tolerance=0.03)
    out = {}
    for label, extra in (("fixed_source", {}),
                         ("frozen_ht_table", {"em_table": table, "ht_sol": hs,
                                              "em_table_azimuths": 8,
                                              "allow_frozen_ht": True})):
        res = calc_heat_axisym.solve_heat_axisym("<meridian>", **common,
                                                 **extra)
        if "error" in res:
            raise RuntimeError(f"route B {label}: {res['error']}")
        out[label] = {"t_s": res["t_history_s"][1:],
                      "T_max_C": res["T_max_history_C"][1:],
                      "Q_input_J": res["Q_input_J"]}
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--frequency", type=float, default=7000.0)
    ap.add_argument("--current", type=float, default=6000.0)
    ap.add_argument("--t-end", type=float, default=4.0)
    ap.add_argument("--ladder", default="0.1:1,0.05:1,0.025:1,0.1:2,0.1:5",
                    help="dt:em_every pairs for route A")
    ap.add_argument("--out", default=os.path.join(
        HERE, "results", "coupled_curie_cylinder_frozen_ht.json"))
    ap.add_argument("--work", default=os.path.join(
        os.environ.get("TEMP", "C:\\temp"), "coupled_curie_cylinder"))
    a = ap.parse_args()
    os.makedirs(a.work, exist_ok=True)
    mesh = build_mesh()
    ladder = [(float(p.split(":")[0]), int(p.split(":")[1]))
              for p in a.ladder.split(",")]
    runs = []
    for dt, every in ladder:
        run = route_a(mesh, a.frequency, a.current, dt, a.t_end, every)
        runs.append(run)
        print(f"route A dt={dt} em_every={every}: T_max(end)="
              f"{run['T_max_C'][-1]:.1f} C  P(end)={run['P_W'][-1]:.0f} W  "
              f"({run['runtime_s']:.0f}s)", flush=True)
    em = C.AxisymEddyCurrent(mesh, frequency=a.frequency, workpiece="wp",
                             coils={"coil": a.current},
                             dirichlet="axis|outer|top|bot")
    s0, m0 = em_material().evaluate(np.full(len(em.workpiece_elements), T0))
    rec0 = em.solve(s0, m0)
    zs = np.linspace(-H_WP / 2 * 0.95, H_WP / 2 * 0.95, 41)
    q0, h0 = wall_source(em, zs)
    rb = route_b(a.work, zs, q0, h0, a.frequency, ladder[0][0], a.t_end)
    finest = min(runs, key=lambda r: r["dt_s"] * r["em_every"])
    import ngsolve
    result = {
        "case": "coaxial coil around a steel-like cylinder through the "
                "Curie band",
        "geometry_m": {"R_wp": R_WP, "H_wp": H_WP,
                       "coil_r": [COIL_R0, COIL_R1], "coil_h": COIL_H,
                       "air_R": R_AIR},
        "frequency_Hz": a.frequency, "current_A": a.current,
        "t_end_s": a.t_end,
        "thermal": {"rho": RHO, "cp": CP, "k": K, "T0_C": T0},
        "em_material": {"T_C": T_TAB, "sigma_S_m": SIGMA_TAB, "mu_r": MU_TAB},
        "P0_W": rec0["P_joule_W"],
        "route_A_ladder": runs,
        "route_A_reference": {"dt_s": finest["dt_s"],
                              "em_every": finest["em_every"]},
        "route_B_surface": rb,
        "wall_source_T0": {"z_m": zs.tolist(), "q_W_m2": q0.tolist(),
                           "Ht_A_m": h0.tolist()},
        "host": platform.node(), "python": platform.python_version(),
        "ngsolve": ngsolve.__version__,
    }
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    with open(a.out, "w", encoding="utf-8") as fh:
        json.dump(result, fh, indent=1)
    ref = finest
    tB = np.asarray(rb["fixed_source"]["t_s"])
    for t in (1.0, 2.0, 3.0, a.t_end):
        iA = int(np.argmin(np.abs(np.asarray(ref["t_s"]) - t)))
        iB = int(np.argmin(np.abs(tB - t)))
        print(f"t={t:.1f}s  route A T_max={ref['T_max_C'][iA]:.0f} C "
              f"(P={ref['P_W'][iA]:.0f} W)  fixed q: "
              f"{rb['fixed_source']['T_max_C'][iB]:.0f} C  frozen |H_t|: "
              f"{rb['frozen_ht_table']['T_max_C'][iB]:.0f} C")
    print("wrote", a.out)


if __name__ == "__main__":
    main()
