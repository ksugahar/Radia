#!/usr/bin/env python3
"""verify_coil_field_independent.py -- end-to-end SF coil validation (A).

Design a coil with the GENERAL FE-direct stream-function stack
(``calc_streamfunction.py``), discretise it into equal-current contour turns,
and verify the designed coil's ACTUAL field over a DSV with THREE INDEPENDENT
field engines:

  1. ``_segment_field_B``               numpy straight-segment Biot-Savart -- the
                                        DESIGN routine (validated vs the
                                        circular-loop analytic to 6 digits).
  2. ``kelvin_source.biot_savart_B_cf`` an independent closed-form Biot-Savart
                                        CoefficientFunction (NGSolve, written
                                        for the FEM-SIBC scattered formulation
                                        -- a DIFFERENT implementation).
  3. ``rad.ObjFlmCur`` + ``rad.Fld``    Radia's C++ filament Biot-Savart -- an
                                        INDEPENDENT CODEBASE.

If the engines agree on the designed coil's field AND that field meets the
target over the DSV, the SF stack is validated end to end: the design routine
is not merely self-consistent -- an independent C++ codebase confirms it.

MRI-gradient-coil scale (cylinder r=0.15 m, L=0.5 m, DSV sphere r=0.05 m):
  uniform   target Bz = 1   (axisymmetric -> full rings; the safe ObjFlmCur path
                             where the documented tilted-loop ObjFlmCur bug does
                             not apply -> a true cross-CODEBASE check)
  gradient  target Bz = x   (transverse Gx fingerprint -- the general saddle
                             case; ObjFlmCur may diverge here by the documented
                             tilted-loop bug, flagged, not a stack error)

Run:  python verify_coil_field_independent.py [--order 2] [--nlevels 8]
"""
import argparse
import json
import os
import sys
import types

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_SRC = os.path.abspath(os.path.join(_HERE, "..", "..", "src"))
sys.path.insert(0, _SRC)

from radia.panels.calc_streamfunction import (         # noqa: E402
    _build_problem, _solve_and_metrics, _surface_triangles, _char_len,
    _contour_levels, _contour_segments, _segs_to_polylines,
    _build_gradpsi_kdt, _gradpsi_signs, _close_loop, _segment_field_B,
    _best_fit_current,
)


def _make_meshes(outdir, a=0.15, L=0.50, dsv=0.05):
    """Cylinder lateral surface (coil) + DSV sphere volume -> two .vol files."""
    from netgen.occ import Cylinder, Sphere, Pnt, Z, OCCGeometry
    from ngsolve import TaskManager
    os.makedirs(outdir, exist_ok=True)
    with TaskManager():     # caller wraps mesh-gen (TaskManager-Only policy)
        cyl = Cylinder(Pnt(0, 0, -L / 2), Z, r=a, h=L)
        lateral = max(cyl.faces, key=lambda f: f.mass)
        coil_path = os.path.join(outdir, "vc_coil_cyl.vol")
        OCCGeometry(lateral).GenerateMesh(maxh=0.04).Save(coil_path)
        sph = Sphere(Pnt(0, 0, 0), dsv)
        eval_path = os.path.join(outdir, "vc_eval_dsv.vol")
        OCCGeometry(sph).GenerateMesh(maxh=0.025).Save(eval_path)
    return coil_path, eval_path


def _design_loops(coil_vol, eval_vol, target_cf, order, nlevels,
                  confine="off"):
    """Run the FE-direct SF design, contour psi into orientation-consistent
    equal-current turns, return (loops_closed, signs, design_homogeneity)."""
    args = types.SimpleNamespace(
        coil_vol=coil_vol, eval_vol=eval_vol, target_cf=target_cf,
        order=order, regularize="h1", alpha=0.0, eval_max=400, confine=confine)
    P = _build_problem(args)
    psi_f, _fit, homo = _solve_and_metrics(P, 0.0)
    psi = P["R"] @ psi_f          # independent-DOF -> all-DOF (Abe eq.6)
    coil = P["coil"]
    verts = np.array([list(v.point) for v in coil.vertices])
    psi_v = psi[:coil.nv]
    tris = _surface_triangles(coil)
    tol = max(1e-12, 1e-6 * _char_len(verts, tris))
    loops = []
    for lev in _contour_levels(psi_v, nlevels):
        loops.extend(_segs_to_polylines(
            _contour_segments(verts, psi_v, tris, lev), tol))
    _kdt, _tK = _build_gradpsi_kdt(verts, psi_v, tris)        # grad-psi winding
    _fields, signs = _gradpsi_signs(loops, _kdt, _tK, P["mpts"], P["vector_b"],
                                    P["Bm"])
    loops_closed = [_close_loop(p) for p in loops]
    n_open = sum(1 for p in loops
                 if np.linalg.norm(p[0] - p[-1]) > 1e-9)
    # obs = DSV sphere vertices pulled 5% toward the centre so they sit
    # STRICTLY inside the eval mesh (boundary vertices fail CF point-location).
    # NOTE: hold the Mesh reference -- a temporary Mesh(eval_vol) is GC'd while
    # the .vertices comprehension iterates, so v.point reads freed memory.
    from ngsolve import Mesh
    em = Mesh(eval_vol)
    sv = np.array([list(v.point) for v in em.vertices])
    obs = sv.mean(axis=0) + 0.95 * (sv - sv.mean(axis=0))
    return dict(loops=loops_closed, signs=signs, design_homo=homo,
                obs=obs, eval_vol=eval_vol, n_open=n_open,
                vector_b=P["vector_b"], target_cf=target_cf)


# ---- three independent field engines (B over the DSV obs points) ----
def _field_segment(loops, signs, I, obs):
    """Engine 1: numpy straight-segment Biot-Savart (the design routine)."""
    B = np.zeros((len(obs), 3))
    for p, s in zip(loops, signs):
        B += _segment_field_B(p, obs, I * s)
    return B


def _field_ngsolve(loops, signs, I, eval_vol, obs):
    """Engine 2: kelvin_source.biot_savart_B_cf (independent NGSolve impl);
    evaluated per loop (small CFs) and summed in numpy."""
    from ngsolve import Mesh
    from radia.kelvin_source import biot_savart_B_cf
    m = Mesh(eval_vol)
    B = np.zeros((len(obs), 3))
    for p, s in zip(loops, signs):
        segs = [(p[k], p[k + 1]) for k in range(len(p) - 1)]
        cf = biot_savart_B_cf([segs], [I * s])
        for i, o in enumerate(obs):
            B[i] += np.array(cf(m(o[0], o[1], o[2])), dtype=float)
    return B


def _field_radia(loops, signs, I, obs):
    """Engine 3: Radia C++ rad.ObjFlmCur + rad.Fld (independent codebase)."""
    import radia as rad
    rad.UtiDelAll()
    objs = []
    for p, s in zip(loops, signs):
        pts = [[float(q[0]), float(q[1]), float(q[2])] for q in p]
        objs.append(rad.ObjFlmCur(pts, float(I * s)))
    cont = rad.ObjCnt(objs)
    B = np.asarray(rad.Fld(cont, "b", obs.tolist()), dtype=float).reshape(-1, 3)
    rad.UtiDelAll()
    return B


def _rel_diff(Ba, Bb):
    """Vector-difference relative error ||Ba-Bb|| / ||Bb|| over the grid."""
    num = float(np.linalg.norm(Ba - Bb))
    den = float(np.linalg.norm(Bb)) + 1e-30
    return num / den


def _uniformity(Bz):
    return float((Bz.max() - Bz.min()) / (abs(Bz.mean()) + 1e-30))


def _gradient_linearity(obs, Bz):
    """Fit Bz = b + G*x; return (G, relative nonlinearity over the DSV)."""
    x = obs[:, 0]
    A = np.column_stack([np.ones_like(x), x])
    coef, *_ = np.linalg.lstsq(A, Bz, rcond=None)
    b, G = coef
    resid = Bz - (b + G * x)
    swing = abs(G) * (x.max() - x.min()) + 1e-30
    return float(G), float(np.linalg.norm(resid) / np.sqrt(len(x)) / swing)


def run_case(name, target_cf, coil_vol, eval_vol, order, nlevels,
             with_ngsolve=False, confine="off"):
    print(f"\n=== case '{name}': target Bz = {target_cf} "
          f"(confine={confine}) ===")
    d = _design_loops(coil_vol, eval_vol, target_cf, order, nlevels, confine)
    loops, signs, obs = d["loops"], d["signs"], d["obs"]
    print(f"  design homogeneity (continuous psi) = {d['design_homo']:.3e}")
    print(f"  turns = {len(loops)} (open contours closed: {d['n_open']})")

    # best-fit ONE series current from the design routine, then reuse for all
    Bseg_unit = _field_segment(loops, signs, 1.0, obs)
    tgt = np.array([eval(target_cf, {"__builtins__": {}},
                         {"x": o[0], "y": o[1], "z": o[2]}) for o in obs],
                   dtype=float)
    I = _best_fit_current(Bseg_unit[:, 2], tgt)

    Bseg = _field_segment(loops, signs, I, obs)
    Brad = _field_radia(loops, signs, I, obs)             # Radia C++ (codebase)
    r_seg_rad = _rel_diff(Bseg, Brad)
    print(f"  best-fit series current I = {I:.3f} A")
    print(f"  ENGINE AGREEMENT (vector-difference relative):")
    print(f"    segment(numpy, design) vs Radia C++  = {r_seg_rad:.3e}")
    res = dict(case=name, target_cf=target_cf, confine=confine,
               n_turns=len(loops), n_open=d["n_open"],
               design_homogeneity=d["design_homo"],
               best_fit_current_A=I, rel_segment_vs_radia=r_seg_rad)
    if with_ngsolve:                          # optional 3rd independent impl
        Bng = _field_ngsolve(loops, signs, I, eval_vol, obs)
        res["rel_segment_vs_ngsolve"] = _rel_diff(Bseg, Bng)
        res["rel_ngsolve_vs_radia"] = _rel_diff(Bng, Brad)
        print(f"    segment(numpy)         vs NGSolve CF = "
              f"{res['rel_segment_vs_ngsolve']:.3e}")
        print(f"    NGSolve CF             vs Radia C++  = "
              f"{res['rel_ngsolve_vs_radia']:.3e}")

    Bz = Bseg[:, 2]
    if "x" in target_cf and target_cf.strip() == "x":
        G, nl = _gradient_linearity(obs, Bz)
        print(f"  gradient G = {G:.4e} T/m, nonlinearity over DSV = {nl:.3e}")
        res.update(gradient_G=G, gradient_nonlinearity=nl)
    else:
        u = _uniformity(Bz)
        print(f"  Bz uniformity (max-min)/mean over DSV = {u:.3e}")
        res.update(bz_uniformity=u)
    return res


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--order", type=int, default=2)
    ap.add_argument("--nlevels", type=int, default=8)
    ap.add_argument("--outdir", default=os.path.join(_HERE, "_verify_vols"))
    ap.add_argument("--cases", default="uniform,gradient,gradient_confined")
    ap.add_argument("--with-ngsolve", action="store_true",
                    help="also run the NGSolve biot_savart_B_cf engine "
                         "(per-point CF eval; can fail point-location)")
    ap.add_argument("--json", default=os.path.join(
        _HERE, "verify_coil_field_independent.json"),
        help="output JSON path (default: next to this script)")
    args = ap.parse_args()

    coil_vol, eval_vol = _make_meshes(args.outdir)
    # name -> (target Bz expr, confine).  gradient_confined shows the
    # method-completing boundary condition: psi=0 on the former edge -> the
    # contours CLOSE on this short former (n_open 0), so the single-current
    # coil is manufacturable AND still cross-validated by Radia C++.
    cases = {"uniform": ("1", "off"),
             "gradient": ("x", "off"),
             "gradient_confined": ("x", "abe")}
    results = []
    for name in args.cases.split(","):
        name = name.strip()
        if name not in cases:
            continue
        tgt, conf = cases[name]
        results.append(run_case(name, tgt, coil_vol, eval_vol,
                                args.order, args.nlevels,
                                with_ngsolve=args.with_ngsolve, confine=conf))

    out = {"meshes": {"coil": os.path.basename(coil_vol),
                      "eval": os.path.basename(eval_vol),
                      "cylinder_r_m": 0.15, "cylinder_L_m": 0.50,
                      "dsv_r_m": 0.05},
           "order": args.order, "nlevels": args.nlevels, "cases": results}
    with open(args.json, "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nsaved {args.json}")


if __name__ == "__main__":
    main()
