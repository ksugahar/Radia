"""3D FEM-SIBC for induction heating: a brute-force axisymmetric check, then
rotor states of a part with a through cross hole.

Part 1 (``axisym``).  A copper tube (axial bore) in a coaxial loop coil is
solved by the 3D FEM-SIBC production solver (air meshed, the tube a hole
with a surface impedance) and, independently, by brute force: the 2D
axisymmetric volumetric eddy-current solver, which resolves the skin in the
metal and assumes nothing about it.  Both truncate at the same sphere with
A = 0.  Compared: the total loss and its split between the outer wall, the
bore and the end faces, at each frequency and on two meshes of each route.

Part 2 (``rotation``).  A cylinder with a radial through hole (genus 1 about
x, which the BEM-SIBC loop extension refuses) turns in front of a side coil.
The FEM-SIBC solver gives one EM solution per rotor angle by two routes:
the part re-meshed turned (world frame) and the part fixed with the coil
turned back (body frame).  Recorded: route agreement of the body-frame
source, convergence in the number of states, the angle-resolution gate, the
power identity, and an end-to-end rotating heat solve for both routes.

Output: results/fem_sibc_through_hole.json

    python validation_test/induction_heating/fem_sibc_through_hole.py
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, os.path.join(ROOT, "src", "radia", "panels"))
sys.path.insert(0, HERE)

import fem_sibc_geometry as G  # noqa: E402
from radia import ih_thermal  # noqa: E402

SIGMA = 5.8e7
R, H = 0.025, 0.025
BORE = 0.008                     # part 1: axial bore of the tube
LOOP_R = 0.035                   # part 1: coaxial loop radius
CROSS = 0.005                    # part 2: through cross hole along x
SIDE_R, SIDE_X = 0.012, 0.036    # part 2: side loop in the plane x = SIDE_X
R_AIR, RING_R, RING_H = 0.15, 0.052, 0.06
PERIOD = math.pi


def skin_depth(f):
    return math.sqrt(2.0 / (2 * math.pi * f * 4e-7 * math.pi * SIGMA))


def fem_solve(vol, stem, f, paths, current, work):
    import calc_fem_kelvin as K
    from radia.em_material import EMMaterial
    t0 = time.perf_counter()
    res = K.solve_fem(vol_file=vol, fes_order=2, frequency=f,
                      mat=EMMaterial(name="copper", sigma=SIGMA, mu_r=1.0),
                      half_thickness=R, solver="sparsecholesky",
                      msh_output=os.path.join(work, stem + ".msh"),
                      filaments=([paths], [current]))
    if "error" in res:
        raise RuntimeError(f"{stem}: {res['error']}")
    res["_runtime_s"] = time.perf_counter() - t0
    return res


def wall_split_3d(res):
    """Loss on the outer wall, the bore and the end faces from the saved
    P1 q_surf (triangles classified by centroid)."""
    from ngsolve import GridFunction, H1, Mesh
    mesh = Mesh(res["qsurf_em_vol"])
    gf = GridFunction(H1(mesh, order=1))
    gf.Load(res["qsurf_sol"])
    field = ih_thermal.SurfaceP1Field.from_gridfunction(mesh, gf, ["sibc"])
    P, T, V = field.points, field.tris, field.values
    c = P[T].mean(axis=1)
    r, z = np.hypot(c[:, 0], c[:, 1]), c[:, 2]
    area = ih_thermal.triangle_areas(P, T)
    pw = area * V[T].mean(axis=1)
    tol = 0.25 * field.h_median
    ends = np.abs(np.abs(z) - H / 2) < tol
    outer = ~ends & (np.abs(r - R) < tol)
    bore = ~ends & (np.abs(r - BORE) < tol)
    if np.any(~(ends | outer | bore)):
        raise RuntimeError("unclassified sibc triangles")
    return {"outer_W": float(pw[outer].sum()), "bore_W": float(pw[bore].sum()),
            "ends_W": float(pw[ends].sum()), "total_W": float(pw.sum())}


def axisym_mesh(f, shell_factor):
    """2D (r, z) half disk of radius R_AIR: tube with skin shells of 4 delta
    meshed at delta / shell_factor, a 0.5 mm square coil, air."""
    from netgen.occ import Glue, MoveTo, OCCGeometry, WorkPlane, X
    from ngsolve import Mesh
    d = skin_depth(f)
    sk = min(4 * d, (R - BORE) / 4)
    hs = d / shell_factor

    def rect(r0, z0, w, h, mh, name="wp"):
        fc = MoveTo(r0, z0).Rectangle(w, h).Face()
        fc.faces.name = name
        fc.maxh = mh
        return fc
    core = rect(BORE + sk, -H / 2 + sk, R - BORE - 2 * sk, H - 2 * sk,
                max(4 * hs, 5e-4))
    shells = [rect(BORE, -H / 2, sk, H, hs), rect(R - sk, -H / 2, sk, H, hs),
              rect(BORE + sk, -H / 2, R - BORE - 2 * sk, sk, hs),
              rect(BORE + sk, H / 2 - sk, R - BORE - 2 * sk, sk, hs)]
    cw = 5e-4
    coil = rect(LOOP_R - cw / 2, -cw / 2, cw, cw, cw / 2, name="coil")
    air = WorkPlane().Circle(0, 0, R_AIR).Face() * \
        MoveTo(0, -R_AIR).Rectangle(R_AIR, 2 * R_AIR).Face()
    air.faces.name = "air"
    air.edges.name = "outer"
    air.edges.Min(X).name = "axis"
    for p in [core, *shells, coil]:
        air = air - p
    return Mesh(OCCGeometry(Glue([air, core, *shells, coil]), dim=2)
                .GenerateMesh(maxh=0.01))


def axisym_solve(f, shell_factor, current):
    from radia import ih_axisym_coupled as C
    from ngsolve import Integrate, x as rr, y as zz
    mesh = axisym_mesh(f, shell_factor)
    em = C.AxisymEddyCurrent(mesh, frequency=f, workpiece="wp",
                             coils={"coil": current}, dirichlet="axis|outer")
    n = len(em.workpiece_elements)
    d = skin_depth(f)
    ratio, _ = em.skin_resolution(np.full(n, d))
    rec = em.solve(np.full(n, SIGMA), np.full(n, 1.0))
    # attribute each element's loss to the nearest face
    Q = em.heat_density()
    pe = np.asarray(Integrate(Q * 2 * math.pi * rr, mesh, element_wise=True,
                              order=8), float)
    cr = np.asarray(Integrate(rr, mesh, element_wise=True), float)
    cz = np.asarray(Integrate(zz, mesh, element_wise=True), float)
    ar = np.asarray(Integrate(1.0, mesh, element_wise=True), float)
    el = np.asarray(em.workpiece_elements)
    cr, cz = cr[el] / ar[el], cz[el] / ar[el]
    dist = np.stack([R - cr, cr - BORE, H / 2 - np.abs(cz)], axis=1)
    near = np.argmin(dist, axis=1)
    pw = pe[el]
    return {"P_joule_W": rec["P_joule_W"],
            "power_balance_relative_error": rec["power_balance_relative_error"],
            "outer_W": float(pw[near == 0].sum()),
            "bore_W": float(pw[near == 1].sum()),
            "ends_W": float(pw[near == 2].sum()),
            "skin_h_over_delta": ratio, "elements": int(mesh.ne)}


def part1(work, freqs, meshes, current):
    out = {}
    for f in freqs:
        row = {"skin_depth_m": skin_depth(f), "fem_sibc": {}, "axisym": {}}
        for label, (mp, mr) in meshes.items():
            vol = os.path.join(work, f"tube_{label}.vol")
            if not os.path.isfile(vol):
                G.air_mesh(lambda: G.cylinder_part(R, H, bore=BORE),
                           r_air=R_AIR, ring_r=RING_R, ring_h=RING_H, maxh_part=mp,
                           maxh_ring=mr, maxh_far=0.04).ngmesh.Save(vol)
            res = fem_solve(vol, f"tube_{label}_{int(f)}", f,
                            G.loop_paths((0, 0, 0), (0, 0, 1), LOOP_R),
                            current, work)
            row["fem_sibc"][label] = {
                "P_total_W": res["P_total"], "ndof": res.get("ndof"),
                "runtime_s": res["_runtime_s"], **wall_split_3d(res)}
            print(f"f={f:g} FEM-SIBC {label}: P={res['P_total']:.6e} W "
                  f"({res['_runtime_s']:.0f} s)", flush=True)
        for sf in (3, 6):
            a = axisym_solve(f, sf, current)
            row["axisym"][f"shell_delta_over_{sf}"] = a
            print(f"f={f:g} axisym delta/{sf}: P={a['P_joule_W']:.6e} W",
                  flush=True)
        ref = row["axisym"]["shell_delta_over_6"]
        row["relative_to_axisym"] = {
            label: {k: v[k2] / ref[k] - 1 for k, k2 in (
                ("P_joule_W", "P_total_W"), ("outer_W", "outer_W"),
                ("bore_W", "bore_W"), ("ends_W", "ends_W"))}
            for label, v in row["fem_sibc"].items()}
        out[f"{f:g}"] = row
    return out


def part2(work, f, ns, current, maxh, hole_maxh, tol, thermal_maxh, heat):
    from rotating_cross_hole_rotor_states import average_source, heat_run, rel
    def part():
        return G.cylinder_part(R, H, cross_hole=CROSS, hole_maxh=hole_maxh)
    sets, out = {}, {"frequency_Hz": f, "runs": {},
                     "em_maxh_m": {"part": maxh[0], "ring": maxh[1],
                                   "hole_walls": hole_maxh}}
    body_vol = os.path.join(work, "cross_body.vol")
    for frame in ("world", "body"):
        for n in ns:
            entries, recs = [], []
            for k in range(n):
                a = PERIOD * k / n
                if frame == "world":
                    vol = os.path.join(work, f"cross_world_N{n}_{k}.vol")
                    G.air_mesh(part, angle=a, r_air=R_AIR, ring_r=RING_R,
                               ring_h=RING_H, maxh_part=maxh[0],
                               maxh_ring=maxh[1],
                               maxh_far=0.04).ngmesh.Save(vol)
                    coil = G.loop_paths((SIDE_X, 0, 0), (1, 0, 0), SIDE_R)
                else:
                    vol = body_vol
                    if not os.path.isfile(vol):
                        G.air_mesh(part, r_air=R_AIR, ring_r=RING_R,
                                   ring_h=RING_H, maxh_part=maxh[0],
                                   maxh_ring=maxh[1],
                                   maxh_far=0.04).ngmesh.Save(vol)
                    coil = G.loop_paths((SIDE_X, 0, 0), (1, 0, 0), SIDE_R,
                                        angle=-a)
                res = fem_solve(vol, f"cross_{frame}_N{n}_{k}", f, coil,
                                current, work)
                entries.append({"angle_rad": a,
                                "qsurf_sol": os.path.relpath(res["qsurf_sol"],
                                                             work),
                                "em_vol": os.path.relpath(res["qsurf_em_vol"],
                                                          work)})
                recs.append({"angle_rad": a, "P_total_W": res["P_total"],
                             "runtime_s": res["_runtime_s"]})
            man = os.path.join(work, f"rotor_fem_{frame}_N{n}.json")
            with open(man, "w", encoding="utf-8") as fh:
                json.dump({"schema": ih_thermal.ROTOR_STATES_SCHEMA,
                           "axis": "z", "period_rad": PERIOD, "frame": frame,
                           "states": entries}, fh, indent=1)
            sets[(frame, n)] = man
            P = np.array([r["P_total_W"] for r in recs])
            print(f"{frame} N={n}: P {P.min():.4e}..{P.max():.4e} W", flush=True)
            out["runs"][f"{frame}_N{n}"] = {"em_states": recs}
    from ngsolve import Mesh
    from netgen.occ import OCCGeometry
    tpart = G.cylinder_part(R, H, cross_hole=CROSS, hole_maxh=hole_maxh)
    tpart.faces.name = "surf"
    tpart.solids.name = "wp"
    thermal = Mesh(OCCGeometry(tpart).GenerateMesh(maxh=thermal_maxh))
    avgs, snaps = {}, {}
    for (frame, n), man in sets.items():
        entry = out["runs"][f"{frame}_N{n}"]
        try:
            _, avg, snap, rec = average_source(man, thermal, tol)
            avgs[(frame, n)], snaps[(frame, n)] = avg, snap
            entry["average"] = rec
        except ValueError as exc:
            entry["refused"] = str(exc)
        print(frame, n, entry.get("refused", "ok")[:200], flush=True)
    cmp = {}
    for n in ns:
        if ("world", n) in avgs and ("body", n) in avgs:
            cmp[f"world_vs_body_N{n}_average"] = rel(avgs[("body", n)],
                                                     avgs[("world", n)])
            for th in snaps[("world", n)]:
                cmp[f"world_vs_body_N{n}_at_{th}"] = rel(
                    snaps[("body", n)][th], snaps[("world", n)][th])
    good = [n for n in ns if ("world", n) in avgs]
    for n0, n1 in zip(good, good[1:]):
        cmp[f"world_N{n0}_vs_N{n1}_average"] = rel(avgs[("world", n0)],
                                                   avgs[("world", n1)])
    out["comparison"] = cmp
    if not good:
        out["heat"] = {"skipped": "no rotor-state set passed the gate"}
        return out
    finest = max(good)
    out["heat"] = {}
    for frame in ("world", "body"):
        if (frame, finest) not in avgs:
            continue
        h = heat_run(sets[(frame, finest)], thermal, heat["rpm"], heat["dt"],
                     heat["t_end"], tol)
        out["heat"][f"{frame}_N{finest}"] = h
    out["heat_settings"] = heat
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--parts", default="axisym,rotation")
    ap.add_argument("--frequencies", default="10000,150000")
    ap.add_argument("--axisym-meshes", default="coarse:0.004:0.006,"
                                               "fine:0.0025:0.004")
    ap.add_argument("--current", type=float, default=100.0)
    ap.add_argument("--rotation-frequency", type=float, default=150000.0)
    ap.add_argument("--n-states", default="16,32")
    ap.add_argument("--rotation-maxh", default="0.0025:0.004")
    ap.add_argument("--rotation-current", type=float, default=500.0)
    ap.add_argument("--hole-maxh", type=float, default=0.001,
                    help="mesh size on the cross-hole walls, whose edges "
                         "carry the field singularity")
    ap.add_argument("--angle-step-tolerance", type=float, default=None)
    ap.add_argument("--thermal-maxh", type=float, default=0.002)
    ap.add_argument("--rpm", type=float, default=600.0)
    ap.add_argument("--dt", type=float, default=0.025)
    ap.add_argument("--t-end", type=float, default=0.5)
    ap.add_argument("--work", default=os.path.join(
        os.environ.get("TEMP", "C:\\temp"), "fem_sibc_through_hole"))
    ap.add_argument("--out", default=os.path.join(
        HERE, "results", "fem_sibc_through_hole.json"))
    ap.add_argument("--source-commit", default="")
    a = ap.parse_args()
    from _provenance import provenance
    prov = provenance(sys.argv, a.source_commit, [
        os.path.abspath(__file__), os.path.join(HERE, "fem_sibc_geometry.py"),
        os.path.join(HERE, "rotating_cross_hole_rotor_states.py"),
        os.path.join(ROOT, "src", "radia", "panels", "calc_fem_kelvin.py"),
        os.path.join(ROOT, "src", "radia", "ih_axisym_coupled.py"),
        os.path.join(ROOT, "src", "radia", "ih_thermal.py"),
        os.path.join(ROOT, "src", "radia", "panels", "calc_heat.py")])
    os.makedirs(a.work, exist_ok=True)
    parts = a.parts.split(",")
    out = {"case": "3D FEM-SIBC: axisymmetric brute-force check and rotor "
                   "states of a part with a through cross hole",
           "geometry_m": {"R": R, "H": H, "bore": BORE, "loop_radius": LOOP_R,
                          "cross_hole": CROSS, "side_loop_radius": SIDE_R,
                          "side_loop_plane_x": SIDE_X, "air_radius": R_AIR},
           "sigma_S_m": SIGMA, "provenance": prov}
    partial = a.out + ".partial.json"
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    unknown = set(parts) - {"axisym", "rotation"}
    if unknown:
        raise SystemExit(f"unknown parts {sorted(unknown)}")
    for name in parts:              # in the order given, each saved at once
        if name == "axisym":
            meshes = {m.split(":")[0]: (float(m.split(":")[1]),
                                        float(m.split(":")[2]))
                      for m in a.axisym_meshes.split(",")}
            out["axisym"] = part1(a.work, [float(v) for v in
                                           a.frequencies.split(",")],
                                  meshes, a.current)
        else:
            mh = [float(v) for v in a.rotation_maxh.split(":")]
            out["rotation"] = part2(
                a.work, a.rotation_frequency,
                [int(v) for v in a.n_states.split(",")], a.rotation_current,
                mh, a.hole_maxh, a.angle_step_tolerance, a.thermal_maxh,
                {"rpm": a.rpm, "dt": a.dt, "t_end": a.t_end})
        with open(partial, "w", encoding="utf-8") as fh:
            json.dump(out, fh, indent=1)
    with open(a.out, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=1)
    if os.path.exists(partial):
        os.remove(partial)
    print("wrote", a.out)


if __name__ == "__main__":
    from ngsolve import TaskManager
    with TaskManager():
        main()
