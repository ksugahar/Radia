"""demo_regcoil_parity_deliverable.py -- (a) of "surpass NESCOIL/REGCOIL".

REGCOIL/NESCOIL/FOCUS stop at the winding-surface CURRENT POTENTIAL psi (a field
and its iso-contours).  This shows that on the SAME vacuum problem radia is a
strict SUPERSET of the deliverable: it reaches the REGCOIL B.n parity AND, from
the same run, emits a windable artifact the design codes do not --

    target B.n on the plasma boundary
      -> design psi on the winding surface   (B.n residual ~ machine precision
         on a PRODUCIBLE target == REGCOIL forward-map parity)
      -> iso-contours of psi  (= the modular coils)
      -> STEP CAD  (OCC WriteStep)           <- REGCOIL does not emit this
      -> PEEC L, R  (one closed contour)     <- nor this

So (a) is not a tie: design AT PARITY, deliverable BEYOND.  This is the
"design-to-manufacture vs design-only" claim, measured.  (For STELLARATOR modular
coils the deliverable is the contour-coils as STEP + their inductance, not a
single wire; the single-stroke wire is the win for single-conductor MRI/IH coils
-- see radia_mcp.fusion_reactor 'stellarator_coil_design' for the honest nuance.)

Reuses the tested REGCOIL design core from demo_regcoil_fusion.py
(_torus_surface_vol / _plasma_points_normals / _A_normal / _design / _contours)
and the production manufacture helpers from calc_streamfunction.py
(_write_step_polylines / _peec_inductance / _close_loop).

Run:  python validation_test/stream_function/demo_regcoil_parity_deliverable.py
Golden: tests/panels/test_regcoil_parity_deliverable_golden.py
"""
from __future__ import annotations

import argparse
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
DOCS_STREAM_FUNCTION = os.path.join(REPO, "docs", "stream_function")
sys.path.insert(0, os.path.join(REPO, "src", "radia"))
sys.path.insert(0, os.path.join(REPO, "src", "radia", "panels"))
sys.path.insert(0, DOCS_STREAM_FUNCTION)       # to import docs-local fusion demo


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--work-dir",
                    default=os.path.join(os.environ.get("TEMP", "/tmp"),
                                         "sf_regcoil_deliverable"))
    ap.add_argument("--eval-max", type=int, default=160,
                    help="plasma-boundary B.n sample points")
    ap.add_argument("--wind-maxh", type=float, default=0.05)
    ap.add_argument("--plasma-maxh", type=float, default=0.04)
    ap.add_argument("--n-levels", type=int, default=14,
                    help="current-potential contours = modular coils")
    ap.add_argument("--step-output", default="",
                    help="STEP path (default <work-dir>/regcoil_coils.step)")
    ap.add_argument("--wire-diam", type=float, default=2.0e-3)
    ap.add_argument("--peec-freq", type=float, default=1.0e3)
    args = ap.parse_args()

    os.makedirs(args.work_dir, exist_ok=True)
    step_path = args.step_output or os.path.join(args.work_dir,
                                                 "regcoil_coils.step")

    import calc_streamfunction as C
    import demo_regcoil_fusion as D
    from ngsolve import Mesh, H1, specialcf, TaskManager

    with TaskManager():
        coil = Mesh(D._torus_surface_vol(
            D.A_WIND, args.wind_maxh, os.path.join(args.work_dir, "wind.vol")))
        plasma = Mesh(D._torus_surface_vol(
            D.A_PLASMA, args.plasma_maxh,
            os.path.join(args.work_dir, "plasma.vol")))
        fes = H1(coil, order=1, definedon=coil.Boundaries(".*"))
        n_cf = specialcf.normal(3)

        pts, nrm, _theta, _phi = D._plasma_points_normals(plasma, args.eval_max)
        A_n = D._A_normal(C, fes, n_cf, pts, nrm)

        # PRODUCIBLE target == REGCOIL forward-map parity case: uniform vertical
        # B.n (a PF / equilibrium / vertical-field coil) reproduces to ~machine
        # precision, exactly as REGCOIL/NESCOIL do.
        Bn = nrm[:, 2]
        psi, res, peak_grad = D._design(C, fes, coil, A_n, Bn, "h1", 1.0e-7)
        loops = D._contours(C, coil, psi, args.n_levels)
        if not loops:
            raise RuntimeError("design produced no contours")

        # ---- DELIVERABLE that REGCOIL/NESCOIL/FOCUS do NOT emit ----
        # (1) STEP CAD of the coil contours (OCC WriteStep)
        C._write_step_polylines(loops, step_path)
        step_ok = os.path.exists(step_path) and os.path.getsize(step_path) > 0
        # (2) PEEC L, R of one closed contour (a modular-coil turn)
        biggest = max(loops, key=len)
        peec = C._peec_inductance(C._close_loop(biggest), args.wire_diam,
                                  5.8e7, args.peec_freq)

    result = {
        "demo": "regcoil_parity_deliverable",
        "target": "uniform_vertical_PF (producible)",
        "regularize": "h1", "alpha_rel": 1.0e-7,
        "design": {                                # == what REGCOIL also produces
            "bn_residual_rel": res,
            "peak_grad_psi": peak_grad,
            "n_contours": len(loops),
        },
        "deliverable": {                           # what REGCOIL does NOT produce
            "step_file": step_path,
            "step_ok": bool(step_ok),
            "step_bytes": (os.path.getsize(step_path) if step_ok else 0),
            "peec_L_H": peec["L_H"],
            "peec_R_ohm": peec["R_ohm"],
            "peec_freq_Hz": peec["freq_Hz"],
            "peec_n_nodes": peec["n_nodes"],
        },
        "note": "REGCOIL/NESCOIL stop at the current potential psi (the B.n "
                "design). radia reaches the SAME B.n parity AND, from the same "
                "run, emits STEP CAD + PEEC L,R -- design-to-manufacture, not "
                "design-only.",
    }
    print(json.dumps(result))


if __name__ == "__main__":
    main()
