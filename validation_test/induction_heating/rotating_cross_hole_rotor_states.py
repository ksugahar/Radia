"""A rotating part with cross holes, heated by a side coil: rotor states
from the production SIBC eddy-current solve.

A copper cylinder with two opposite blind radial holes (two-fold) turns about its
axis in front of a circular coil facing its side wall.  The heat source seen
by the part changes with its angle, so the thermal solve needs one EM
solution per rotor angle (``--rotor-states``).  Two independent routes build
them:

``world``  the part is re-meshed turned by each angle; the coil stays put;
``body``   the part mesh stays put; the coil is turned by minus the angle.

Both describe the same physics with different meshes and coil
discretisations, so their revolution-averaged body-frame heat flux must
agree.  The script also records:

* convergence of the revolution average in the number of states N per
  half turn (the symmetry period),
* the angle-resolution gate on the coarsest set,
* the identity between the power of the averaged field and the mean knot
  power, and the transfer error of every knot against its EM power,
* an end-to-end rotating heat solve for both routes (energy balance and the
  peak temperature).

Output: results/rotating_cross_hole_rotor_states.json

    python validation_test/induction_heating/rotating_cross_hole_rotor_states.py
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

from radia import ih_thermal  # noqa: E402

R, H = 0.025, 0.025                 # cylinder radius, height
R_HOLE, DEPTH = 0.005, 0.010       # two opposite blind radial holes
COIL_R, COIL_X = 0.012, 0.036       # coil loop in the plane x = COIL_X
FREQ, CURRENT, SIGMA = 10000.0, 2000.0, 5.8e7
PERIOD = math.pi                    # two-fold part


def _part():
    """Cylinder with two opposite blind radial holes along x (two-fold).

    Blind holes keep the part simply connected; a through cross hole makes
    it genus 1 about x, which the SIBC loop extension (genus 1 about z only)
    refuses."""
    from netgen.occ import Axes, Cylinder, Pnt, X, Z
    body = Cylinder(Axes(Pnt(0, 0, -H / 2), Z), r=R, h=H)
    length = DEPTH + 0.002
    for x0 in (R - DEPTH, -R - 0.002):
        body = body - Cylinder(Pnt(x0, 0, 0), X, r=R_HOLE, h=length)
    return body


def part_surface(angle, maxh):
    """Surface mesh of the part turned by ``angle`` about z (faces 'sibc')."""
    from netgen.occ import Axis, Glue, OCCGeometry, Z
    from ngsolve import Mesh
    part = _part()
    if angle:
        part = part.Rotate(Axis((0, 0, 0), Z), math.degrees(angle))
    for f in part.faces:
        f.name = "sibc"
    return Mesh(OCCGeometry(Glue(list(part.faces))).GenerateMesh(maxh=maxh))


def part_volume(maxh):
    """Thermal volume mesh of the unturned part (faces 'surf')."""
    from netgen.occ import OCCGeometry
    from ngsolve import Mesh
    part = _part()
    part.faces.name = "surf"
    part.solids.name = "wp"
    return Mesh(OCCGeometry(part).GenerateMesh(maxh=maxh))


def coil_paths(angle, n_seg=720):
    """Filament loop of radius COIL_R in the plane x = COIL_X, turned by
    ``angle`` about z."""
    s = np.linspace(0, 2 * np.pi, n_seg + 1)
    pts = np.stack([np.full_like(s, COIL_X), COIL_R * np.cos(s),
                    COIL_R * np.sin(s)], axis=1)
    c, sn = math.cos(angle), math.sin(angle)
    rot = np.array([[c, -sn, 0], [sn, c, 0], [0, 0, 1]])
    pts = pts @ rot.T
    return [np.stack([pts[:-1], pts[1:]], axis=1)]


def em_state(work, tag, angle, frame, maxh, backend):
    """One production SIBC solve; returns the manifest entry and a record."""
    import calc_inductance as calc
    part_angle = angle if frame == "world" else 0.0
    coil_angle = 0.0 if frame == "world" else -angle
    vol = os.path.join(work, f"{tag}.vol")
    if frame == "world" or not os.path.isfile(vol):
        part_surface(part_angle, maxh).ngmesh.Save(vol)
    stem = f"{tag}_{frame}_{angle:.6f}"
    args = calc.build_argparser().parse_args([
        "--coil-solver", "peec", "--coil-step", "not-used.step",
        "--vol", vol, "--wp-label", "sibc", "--frequency", str(FREQ),
        "--current", str(CURRENT), "--sigma", str(SIGMA), "--mu-r", "1",
        "--h1-order", "1", "--wp-bem-backend", backend,
        "--msh-output", os.path.join(work, stem + ".msh")])
    source = dict(source_type="filament", paths=coil_paths(coil_angle),
                  I_fil=np.array([CURRENT + 0j]))
    t0 = time.perf_counter()
    res = calc._solve_workpiece_weak_coupled(args, source)
    rec = {"angle_rad": angle, "P_wp_W": float(res["P_wp"]),
           "power_balance_relative_error":
               float(res["wp_power_balance_relative_error"]),
           "runtime_s": time.perf_counter() - t0}
    return {"angle_rad": angle,
            "qsurf_sol": os.path.relpath(res["qsurf_sol"], work),
            "em_vol": os.path.relpath(vol, work)}, rec


def build_set(work, frame, n, maxh, backend):
    tag = f"{frame}_N{n}" if frame == "world" else "body_part"
    entries, recs = [], []
    for k in range(n):
        e, r = em_state(work, tag if frame == "body" else f"{tag}_{k}",
                        PERIOD * k / n, frame, maxh, backend)
        entries.append(e)
        recs.append(r)
    man = os.path.join(work, f"rotor_{frame}_N{n}.json")
    with open(man, "w", encoding="utf-8") as fh:
        json.dump({"schema": ih_thermal.ROTOR_STATES_SCHEMA, "axis": "z",
                   "period_rad": PERIOD, "frame": frame,
                   "states": entries}, fh, indent=1)
    return man, recs


def average_source(man, thermal, tolerance):
    from ngsolve import BND, GridFunction, H1, Integrate
    states, meta = ih_thermal.load_rotor_states(man)
    gf = GridFunction(H1(thermal, order=1))
    rot = ih_thermal.RotatingSurfaceSource(
        thermal, ["surf"], gf, states=states, power_tolerance=0.03,
        axis="z", period=meta["period_rad"], frame=meta["frame"],
        angle_step_tolerance=tolerance)
    snap = {}
    for theta in (0.0, PERIOD / 4, PERIOD / 2):
        rot.at(theta)
        snap[f"{theta:.6f}"] = gf.vec.FV().NumPy()[rot.vnrs].copy()
    rot.average(0.0, 2 * math.pi)
    avg = gf.vec.FV().NumPy()[rot.vnrs].copy()
    p_avg = float(Integrate(gf, thermal, BND,
                            definedon=thermal.Boundaries("surf")))
    knot_p = []
    for a in rot.angles:
        rot.at(a)
        knot_p.append(float(Integrate(gf, thermal, BND,
                                      definedon=thermal.Boundaries("surf"))))
    kp = np.r_[knot_p, knot_p[0]]
    L = np.diff(np.r_[rot.angles, rot.angles[0] + 2 * math.pi])
    p_trap = float(np.sum(0.5 * L * (kp[:-1] + kp[1:])) / (2 * math.pi))
    return rot, avg, snap, {"P_average_field_W": p_avg,
                            "P_knot_mean_W": p_trap,
                            "identity_relative_error": p_avg / p_trap - 1,
                            "knot_P_W": knot_p, **rot.audit()}


def rel(a, b):
    return float(np.linalg.norm(a - b) / np.linalg.norm(b))


def heat_run(man, thermal, rpm, dt, t_end, tolerance):
    import calc_heat
    res = calc_heat.solve_heat(
        "<cross-hole>", material="custom", rho=8960.0, cp=385.0, k=400.0,
        h_conv=0.0, heat_flux_boundaries="surf", rotor_states=man,
        rotation_rpm=rpm, angle_step_tolerance=tolerance,
        dt=dt, t_end=t_end, fes_order=2, _wp_mesh=thermal,
        _write_solution=False, power_tolerance=0.03)
    if "error" in res:
        raise RuntimeError(res["error"])
    return {"T_max_C": res["T_max_C"], "T_mean_C": res["T_mean_C"],
            "Q_input_J": res["Q_input_J"],
            "rotation": res["qsurf_projection"]["rotation"]}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--n-states", default="16,32,64")
    ap.add_argument("--maxh", type=float, default=0.0025)
    ap.add_argument("--thermal-maxh", type=float, default=0.002)
    ap.add_argument("--backend", default="intree-dense")
    ap.add_argument("--rpm", type=float, default=600.0)
    ap.add_argument("--dt", type=float, default=0.025)
    ap.add_argument("--t-end", type=float, default=0.5)
    ap.add_argument("--work", default=os.path.join(
        os.environ.get("TEMP", "C:\\temp"), "rotating_cross_hole"))
    ap.add_argument("--out", default=os.path.join(
        HERE, "results", "rotating_cross_hole_rotor_states.json"))
    ap.add_argument("--source-commit", default="")
    ap.add_argument("--angle-step-tolerance", type=float, default=None,
                    help="angle-resolution gate (default that of the heat "
                         "solver); recorded")
    a = ap.parse_args()
    from _provenance import provenance
    prov = provenance(sys.argv, a.source_commit, [
        os.path.abspath(__file__),
        os.path.join(ROOT, "src", "radia", "ih_thermal.py"),
        os.path.join(ROOT, "src", "radia", "panels", "calc_heat.py"),
        os.path.join(ROOT, "src", "radia", "panels", "calc_inductance.py")])
    os.makedirs(a.work, exist_ok=True)
    import ngsolve
    ngsolve.SetNumThreads(8)
    ns = [int(v) for v in a.n_states.split(",")]
    thermal = part_volume(a.thermal_maxh)
    sets = {}
    with ngsolve.TaskManager():
        for frame in ("world", "body"):
            for n in ns:
                man, recs = build_set(a.work, frame, n, a.maxh, a.backend)
                P = np.array([r["P_wp_W"] for r in recs])
                print(f"{frame} N={n}: P_wp {P.min():.4e}..{P.max():.4e} W "
                      f"({sum(r['runtime_s'] for r in recs):.0f} s)",
                      flush=True)
                sets[(frame, n)] = {"manifest": man, "em_states": recs}
    out = {"case": "copper cylinder with two opposite blind radial holes, side coil",
           "geometry_m": {"R": R, "H": H, "R_hole": R_HOLE, "hole_depth": DEPTH,
                          "coil_radius": COIL_R, "coil_plane_x": COIL_X},
           "frequency_Hz": FREQ, "current_A": CURRENT, "sigma_S_m": SIGMA,
           "symmetry_period_rad": PERIOD, "em_maxh_m": a.maxh,
           "thermal_maxh_m": a.thermal_maxh, "backend": a.backend,
           "angle_step_tolerance": a.angle_step_tolerance,
           "runs": {}, "provenance": prov}
    avgs, snaps = {}, {}
    for (frame, n), s in sets.items():
        key = f"{frame}_N{n}"
        entry = {"em_states": s["em_states"]}
        try:
            _, avg, snap, rec = average_source(s["manifest"], thermal,
                                               a.angle_step_tolerance)
            avgs[(frame, n)], snaps[(frame, n)] = avg, snap
            entry["average"] = rec
        except ValueError as exc:
            entry["refused"] = str(exc)
        out["runs"][key] = entry
        print(key, ("refused: " + entry["refused"][:300]) if "refused" in entry else
              f"P_avg={entry['average']['P_average_field_W']:.4e} W", flush=True)
    # route agreement and angular convergence
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
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    with open(a.out + ".partial.json", "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=1)
    # end-to-end rotating heat solve on the finest set of both routes
    if not good:
        raise SystemExit("no rotor-state set passed the angle-resolution "
                         "gate; see runs[*].refused in the output")
    finest = max(good)
    out["heat"] = {}
    for frame in ("world", "body"):
        t0 = time.perf_counter()
        h = heat_run(sets[(frame, finest)]["manifest"], thermal, a.rpm, a.dt,
                     a.t_end, a.angle_step_tolerance)
        h["runtime_s"] = time.perf_counter() - t0
        out["heat"][f"{frame}_N{finest}"] = h
        print(f"heat {frame}: T_max={h['T_max_C']:.3f} C Q={h['Q_input_J']:.4e}"
              f" J", flush=True)
    out["heat_settings"] = {"rpm": a.rpm, "dt_s": a.dt, "t_end_s": a.t_end,
                            "material": "copper-like constant (8960, 385, "
                                        "400)"}
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    with open(a.out, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=1)
    os.remove(a.out + ".partial.json")
    for k, v in cmp.items():
        print(f"{k}: {v:.3e}")
    print("wrote", a.out)


if __name__ == "__main__":
    main()
