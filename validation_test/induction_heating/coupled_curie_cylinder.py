"""Temperature-dependent induction heating: two independent routes.

Route A (reference): axisymmetric volumetric eddy-current FE for A_phi on an
(r, z) mesh of workpiece + coil + air, re-solved every time step with the
element temperature setting sigma(T) and mu_r(T); the time-averaged Joule
density Q = sigma |omega A|^2 / 2 heats the workpiece volume.

Route B (production): the surface source of route A at the initial
temperature -- q(z) and |H_t|(z) on the heated cylinder wall -- is written
on a 3D surface mesh as a q_surf / _Ht.sol pair and consumed by
calc_heat_axisym --em-table with a surface-impedance table built from the
same sigma(T), mu_r(T) (planar SIBC, q = |H_t|^2 sqrt(omega mu / (2 sigma)) / 2).

Both routes share mesh-independent thermal properties and boundary
conditions, so their difference is the error of route B's model: |H_t|
frozen at the initial state, the heat deposited at the surface instead of
over the skin depth, and the planar impedance.  Below the Curie band the
skin depth is a few tenths of a millimetre and the routes must agree;
above it the skin depth reaches millimetres and the comparison measures
how far the surface model can be trusted.

The material curves are illustrative (a carbon-steel-like shape), not
alloy data.  Results go to ``results/coupled_curie_cylinder.json``.

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

MU0 = 4e-7 * math.pi

# Geometry [m]
R_WP, H_WP = 0.020, 0.040
COIL_R0, COIL_R1, COIL_H = 0.026, 0.030, 0.012
R_AIR = 0.15
# Thermal (constant, shared by both routes)
RHO, CP, K = 7800.0, 600.0, 35.0
T0 = 20.0


def sigma_of_T(T):
    """Illustrative electrical conductivity [S/m]."""
    return np.interp(T, [0.0, 400.0, 800.0, 1200.0],
                     [6.0e6, 2.4e6, 1.0e6, 0.85e6])


def mu_r_of_T(T):
    """Illustrative linear permeability collapsing over 700..770 C."""
    return np.interp(T, [0.0, 700.0, 770.0, 2000.0], [100.0, 80.0, 1.0, 1.0])


def build_mesh(maxh_skin, maxh_wp, maxh_air, skin):
    from netgen.occ import Glue, MoveTo, OCCGeometry, X, Y
    from ngsolve import Mesh

    core = MoveTo(0, -H_WP / 2).Rectangle(R_WP - skin, H_WP).Face()
    core.faces.name = "wp"
    core.edges.Min(X).name = "axis"     # A_phi = 0 on the whole axis
    core.maxh = maxh_wp
    shell = MoveTo(R_WP - skin, -H_WP / 2).Rectangle(skin, H_WP).Face()
    shell.faces.name = "wp"
    shell.maxh = maxh_skin
    coil = MoveTo(COIL_R0, -COIL_H / 2).Rectangle(COIL_R1 - COIL_R0,
                                                   COIL_H).Face()
    coil.faces.name = "coil"
    coil.maxh = 1.0e-3
    air = MoveTo(0, -R_AIR).Rectangle(R_AIR, 2 * R_AIR).Face()
    air = air - core - shell - coil
    air.faces.name = "air"
    air.edges.Max(X).name = "outer"
    air.edges.Min(X).name = "axis"
    air.edges.Max(Y).name = "top"
    air.edges.Min(Y).name = "bot"
    shape = Glue([air, core, shell, coil])
    mesh = Mesh(OCCGeometry(shape, dim=2).GenerateMesh(maxh=maxh_air))
    return mesh


class EddyCurrent:
    """Time-harmonic A_phi with element-wise sigma, mu_r (P0 fields)."""

    def __init__(self, mesh, frequency, current, order=2):
        from ngsolve import H1, L2, GridFunction
        self.mesh = mesh
        self.omega = 2 * math.pi * frequency
        self.fes = H1(mesh, order=order, complex=True,
                      dirichlet="axis|outer|top|bot")
        self.p0 = L2(mesh, order=0)
        self.gf_sigma = GridFunction(self.p0)
        self.gf_mu = GridFunction(self.p0)
        self.gfA = GridFunction(self.fes)
        self.J = current / ((COIL_R1 - COIL_R0) * COIL_H)
        self.wp_elems = np.asarray([el.nr for el in mesh.Elements()
                                    if el.mat == "wp"])

    def solve(self, T_elem):
        """T_elem: temperature per workpiece element (order of wp_elems)."""
        from ngsolve import (BilinearForm, LinearForm, TaskManager, dx, grad,
                             x as r)
        s = np.zeros(self.mesh.ne)
        m = np.ones(self.mesh.ne)
        s[self.wp_elems] = sigma_of_T(T_elem)
        m[self.wp_elems] = mu_r_of_T(T_elem)
        self.gf_sigma.vec.FV().NumPy()[:] = s
        self.gf_mu.vec.FV().NumPy()[:] = m
        nu = 1.0 / (MU0 * self.gf_mu)
        u, v = self.fes.TnT()
        a = BilinearForm(self.fes)
        a += nu / r * (r * grad(u)[0] + u) * (r * grad(v)[0] + v) * dx
        a += nu * r * grad(u)[1] * grad(v)[1] * dx
        a += 1j * self.omega * self.gf_sigma * r * u * v * dx("wp")
        f = LinearForm(self.fes)
        f += self.J * r * v * dx("coil")
        with TaskManager():
            a.Assemble()
            f.Assemble()
            self.gfA.vec.data = a.mat.Inverse(self.fes.FreeDofs(),
                                              inverse="pardiso") * f.vec
        from ngsolve import InnerProduct
        return 0.5 * self.gf_sigma * self.omega ** 2 * InnerProduct(
            self.gfA, self.gfA)

    def wall_source(self, z):
        """q(z) [W/m^2] into the wall and |H_t|(z) [A/m] just outside it."""
        from ngsolve import Integrate, x as r
        # q(z): revolved Joule power of the thin horizontal slab around z,
        # per unit wall area, from the element-wise density.
        Q = self.last_Q
        out_q, out_h = [], []
        dz = 0.5e-3
        for zi in z:
            from ngsolve import IfPos, y
            band = IfPos(dz / 2 - (y - zi), 1.0, 0.0) * IfPos(dz / 2 + (y - zi), 1.0, 0.0)
            P = float(Integrate(Q.real * band * 2 * math.pi * r, self.mesh,
                                definedon=self.mesh.Materials("wp"),
                                order=8).real)
            out_q.append(P / (2 * math.pi * R_WP * dz))
            e = 2e-5
            p1, p2 = R_WP + e, R_WP + 2 * e
            A1 = complex(self.gfA(self.mesh(p1, zi)))
            A2 = complex(self.gfA(self.mesh(p2, zi)))
            Bz = ((p2 * A2 - p1 * A1) / (p2 - p1)) / p1
            out_h.append(abs(Bz) / MU0)
        return np.asarray(out_q), np.asarray(out_h)


def route_a(mesh, frequency, current, dt, t_end, depth_T):
    """Staggered EM + volumetric heat on the workpiece elements."""
    from ngsolve import (BilinearForm, CF, GridFunction, H1, Integrate,
                         LinearForm, TaskManager, dx, grad, InnerProduct,
                         x as r)
    em = EddyCurrent(mesh, frequency, current)
    fesT = H1(mesh, order=2, definedon="wp")
    u, v = fesT.TnT()
    gfT = GridFunction(fesT)
    gfT.Set(CF(T0), definedon=mesh.Materials("wp"))
    w = 2 * math.pi * r
    m = BilinearForm(fesT, symmetric=True)
    m += RHO * CP * u * v * w * dx("wp")
    k = BilinearForm(fesT, symmetric=True)
    k += K * InnerProduct(grad(u), grad(v)) * w * dx("wp")
    with TaskManager():
        m.Assemble()
        k.Assemble()
    mstar = m.mat.CreateMatrix()
    mstar.AsVector().data = m.mat.AsVector() + dt * k.mat.AsVector()
    inv = mstar.Inverse(fesT.FreeDofs(), inverse="sparsecholesky")
    # element-centroid temperatures for the material update
    from ngsolve import VOL
    cent = [(np.mean([mesh.vertices[vv.nr].point[0] for vv in el.vertices]),
             np.mean([mesh.vertices[vv.nr].point[1] for vv in el.vertices]))
            for el in mesh.Elements(VOL) if el.mat == "wp"]
    cent = np.asarray(cent)
    mips = mesh(np.ascontiguousarray(cent[:, 0]),
                np.ascontiguousarray(cent[:, 1]))
    hist = []
    t = 0.0
    em.last_Q = em.solve(np.full(len(cent), T0))
    q0, h0 = None, None
    zs = np.linspace(-H_WP / 2 * 0.95, H_WP / 2 * 0.95, 41)
    q0, h0 = em.wall_source(zs)
    P0 = float(Integrate(em.last_Q.real * w, mesh,
                         definedon=mesh.Materials("wp"), order=8).real)
    n = int(round(t_end / dt))
    for step in range(1, n + 1):
        T_el = np.asarray(gfT(mips)).reshape(-1)
        em.last_Q = em.solve(T_el)
        f = LinearForm(fesT)
        f += em.last_Q.real * v * w * dx("wp", bonus_intorder=4)
        with TaskManager():
            f.Assemble()
            rhs = gfT.vec.CreateVector()
            rhs.data = m.mat * gfT.vec + dt * f.vec
            gfT.vec.data = inv * rhs
        t = step * dt
        P = float(Integrate(em.last_Q.real * w, mesh,
                            definedon=mesh.Materials("wp"), order=8).real)
        hist.append({"t_s": t, "P_W": P,
                     "T_surface_mid_C": float(gfT(mesh(R_WP - 1e-6, 0.0))),
                     "T_axis_mid_C": float(gfT(mesh(1e-6, 0.0)))})
    depth = _radial_depth(lambda rr: float(gfT(mesh(rr, 0.0))), depth_T)
    return {"P0_W": P0, "history": hist, "depth_mid_m": depth,
            "wall_z_m": zs.tolist(), "wall_q0_W_m2": q0.tolist(),
            "wall_Ht0_A_m": h0.tolist()}


def _radial_depth(T_at, thr):
    rr = np.linspace(R_WP - 1e-6, 1e-6, 4000)
    T = np.asarray([T_at(x) for x in rr])
    if T[0] < thr:
        return 0.0
    below = np.flatnonzero(T < thr)
    if below.size == 0:
        return float(R_WP)
    j = below[0]
    d0, d1 = R_WP - rr[j - 1], R_WP - rr[j]
    return float(d0 + (d1 - d0) * (T[j - 1] - thr) / (T[j - 1] - T[j]))


def impedance_table(path, frequency):
    omega = 2 * math.pi * frequency
    H = np.logspace(2, 6.5, 40)
    T = np.linspace(0.0, 1300.0, 261)
    reZ = np.sqrt(omega * MU0 * mu_r_of_T(T) / (2 * sigma_of_T(T)))
    q = 0.5 * reZ[None, :] * H[:, None] ** 2
    np.savez(path, H_grid=H, T_grid=T, Zs_re=np.repeat(reZ[None, :], 40, 0),
             Zs_im=np.repeat(reZ[None, :], 40, 0), q_surf=q,
             meta=np.asarray(json.dumps({"frequency": frequency,
                                         "material": "illustrative-steel",
                                         "model": "linear planar SIBC"})))


def route_b(workdir, a, frequency, dt, t_end, depth_T):
    import calc_heat_axisym
    import ih_thermal
    import ih_thermal_post
    from netgen.geom2d import SplineGeometry
    from netgen.occ import Axes, Cylinder, OCCGeometry, Pnt, Z
    from ngsolve import GridFunction, H1, Mesh, y as _y, z as zc

    solid = Cylinder(Axes(Pnt(0, 0, -H_WP / 2), Z), r=R_WP, h=H_WP)
    solid.faces.name = "side"
    solid.faces.Max(Z).name = "top"
    solid.faces.Min(Z).name = "bottom"
    vol = os.path.join(workdir, "em_cylinder.vol")
    OCCGeometry(solid).GenerateMesh(maxh=0.001).Save(vol)
    em = Mesh(vol)
    zs = np.asarray(a["wall_z_m"])
    pts = ih_thermal.mesh_vertices(em)
    side = ih_thermal.boundary_vertex_numbers(em, ["side"])
    q = GridFunction(H1(em, order=1))
    h = GridFunction(H1(em, order=1))
    q.vec[:] = 0.0
    h.vec[:] = 0.0
    q.vec.FV().NumPy()[side] = np.interp(pts[side, 2], zs, a["wall_q0_W_m2"])
    h.vec.FV().NumPy()[side] = np.interp(pts[side, 2], zs, a["wall_Ht0_A_m"])
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
    for label, extra in (("frozen", {}),
                         ("coupled", {"em_table": table, "ht_sol": hs,
                                      "em_table_azimuths": 8,
                                      "allow_frozen_ht": True})):
        res = calc_heat_axisym.solve_heat_axisym("<meridian>", **common,
                                                 **extra)
        if "error" in res:
            raise RuntimeError(f"route B {label}: {res['error']}")
        out[label] = {k: res[k] for k in ("T_max_C", "Q_input_J",
                                          "T_max_history_C", "t_history_s")}
        out[label]["qsurf_projection"] = {
            k: v for k, v in res["qsurf_projection"].items()
            if k in ("power_balance", "temperature_dependent_source")}
        out[label]["nonlinear_transient"] = res.get("nonlinear_transient")
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--frequency", type=float, default=7000.0)
    ap.add_argument("--current", type=float, default=6000.0,
                    help="coil ampere-turns (peak)")
    ap.add_argument("--dt", type=float, default=0.1)
    ap.add_argument("--t-end", type=float, default=4.0)
    ap.add_argument("--depth-T", type=float, default=850.0)
    ap.add_argument("--out", default=os.path.join(
        HERE, "results", "coupled_curie_cylinder_frozen_ht.json"))
    ap.add_argument("--work", default=os.path.join(
        os.environ.get("TEMP", "C:\\temp"), "coupled_curie_cylinder"))
    a = ap.parse_args()
    os.makedirs(a.work, exist_ok=True)
    t0 = time.time()
    mesh = build_mesh(maxh_skin=6e-5, maxh_wp=1.5e-3, maxh_air=0.01,
                      skin=1.5e-3)
    ra = route_a(mesh, a.frequency, a.current, a.dt, a.t_end, a.depth_T)
    t_a = time.time() - t0
    rb = route_b(a.work, ra, a.frequency, a.dt, a.t_end, a.depth_T)
    t_b = time.time() - t0 - t_a
    import ngsolve
    result = {
        "case": "coaxial coil around a steel-like cylinder through the "
                "Curie band",
        "geometry_m": {"R_wp": R_WP, "H_wp": H_WP, "coil_r": [COIL_R0,
                                                             COIL_R1],
                       "coil_h": COIL_H, "air_R": R_AIR},
        "frequency_Hz": a.frequency, "current_A": a.current,
        "dt_s": a.dt, "t_end_s": a.t_end,
        "thermal": {"rho": RHO, "cp": CP, "k": K, "T0_C": T0},
        "skin_depth_mm": {
            "T0": 1e3 * math.sqrt(2 / (2 * math.pi * a.frequency * MU0 *
                                       mu_r_of_T(T0) * sigma_of_T(T0))),
            "800C": 1e3 * math.sqrt(2 / (2 * math.pi * a.frequency * MU0 *
                                         mu_r_of_T(800.0) *
                                         sigma_of_T(800.0)))},
        "route_A_volumetric": ra,
        "route_B_surface": rb,
        "runtime_s": {"route_A": t_a, "route_B": t_b},
        "host": platform.node(), "python": platform.python_version(),
        "ngsolve": ngsolve.__version__,
    }
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    with open(a.out, "w", encoding="utf-8") as fh:
        json.dump(result, fh, indent=1)
    hA = ra["history"]
    print(f"P0: route A {ra['P0_W']:.1f} W")
    for i in range(0, len(hA), max(1, len(hA) // 8)):
        tA = hA[i]
        tb_f = rb["frozen"]["T_max_history_C"][i + 1]
        tb_c = rb["coupled"]["T_max_history_C"][i + 1]
        print(f"t={tA['t_s']:.2f}s  A: P={tA['P_W']:.0f} W "
              f"T_surf={tA['T_surface_mid_C']:.0f} C   "
              f"B frozen T_max={tb_f:.0f}  B coupled T_max={tb_c:.0f}")
    print("wrote", a.out)


if __name__ == "__main__":
    main()
