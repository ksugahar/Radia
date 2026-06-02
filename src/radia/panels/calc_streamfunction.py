#!/usr/bin/env python3
"""calc_streamfunction.py -- Stream-Function (SF) coil design (Layer 4, headless).

Design a surface coil on a STANDALONE 2D surface mesh (.vol) that produces a
target field over an evaluation region (.vol, surface OR volume).  The target
and the computed field are NGSolve CoefficientFunctions (no .sol / .csv).

I/O (per the 2026-06 design decisions):
  --coil-vol   coil surface mesh (.vol, 2D-in-3D).  psi = H1 on it (FE-direct).
  --eval-vol   evaluation region (.vol).  Surface OR volume; the homogeneity
               adapts to the mesh dimension automatically.
  --target-cf  target field as a CoefficientFunction EXPRESSION of x, y, z.
               SCALAR expr -> Bz target (e.g. "x" = Gx, "1" = uniform Bz,
               "x*x-y*y" = C2 ellipse shim).  3-VECTOR expr "(Bx,By,Bz)" ->
               full vector-B target.

THREE MODES (--method):
  design       target -> A psi = B (folded-Tikhonov RegularizedTSVD) -> psi,
               field homogeneity over the eval region, peak surface current.
  pareto       sweep the Tikhonov alpha -> the (homogeneity, peak current
               density) Pareto front (the --pareto-lever {alpha} L-curve;
               linf / geometry / sheetmetal levers: interface defined, bodies
               are follow-ons).
  manufacture  single-stroke chain + sheet-metal wire distortion + CAD(STEP) +
               PEEC inductance.  (interface defined; body is a follow-on.)

Surface-FE convention (verified): a standalone surface .vol loads with ne=0
(surface elements as boundary), so psi uses
``H1(coil, definedon=coil.Boundaries('.*'))`` with ``grad(v).Trace()`` and
``* ds`` -- the sphere-demo machinery generalised to a loaded mesh.
"""
import os
import sys
import argparse
import time

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)                                   # calc_common
sys.path.insert(0, os.path.join(_HERE, "..", ".."))         # src/radia package root
sys.path.insert(0, os.path.join(_HERE, ".."))               # src/radia

from calc_common import calc_main  # noqa: E402

_MU0_4PI = 1.0e-7


def _parse_target_cf(expr):
    """Eval a CoefficientFunction expression of x, y, z.  Returns (cf, dim)."""
    from ngsolve import x, y, z, CoefficientFunction, sin, cos, exp, sqrt, log
    import math
    scope = {"x": x, "y": y, "z": z, "CF": CoefficientFunction,
             "sin": sin, "cos": cos, "exp": exp, "sqrt": sqrt, "log": log,
             "pi": math.pi}
    val = eval(expr, {"__builtins__": {}}, scope)            # noqa: S307
    cf = CoefficientFunction(val)
    return cf, cf.dim


def _target_values(target_cf, vector_b, mesh, pts):
    """Evaluate the target CF at points -> flat array (1 or 3 per point)."""
    out = []
    for p in pts:
        tv = target_cf(mesh(p[0], p[1], p[2]))
        out.extend(list(tv) if vector_b else [float(tv)])
    return np.array(out, dtype=float)


def _assemble_biot_savart(fes, n, pts, comps):
    """A[row, dof]: field component(s) at each eval point from the FE-direct
    surface current K = n x grad_s phi (Setup B: grad(v).Trace(), * ds)."""
    from ngsolve import grad, ds, LinearForm, sqrt as ng_sqrt, Cross, x, y, z
    ndof = fes.ndof
    u, v = fes.TnT()
    A = np.zeros((len(pts) * len(comps), ndof))
    row = 0
    for p in pts:
        dxt, dyt, dzt = p[0] - x, p[1] - y, p[2] - z
        r2 = dxt * dxt + dyt * dyt + dzt * dzt
        r3 = r2 * ng_sqrt(r2)
        Kv = Cross(n, grad(v).Trace())
        cross = (Kv[1] * dzt - Kv[2] * dyt,
                 Kv[2] * dxt - Kv[0] * dzt,
                 Kv[0] * dyt - Kv[1] * dxt)
        for c in comps:
            Lf = LinearForm(fes)
            Lf += (_MU0_4PI / r3) * cross[c] * ds
            Lf.Assemble()
            A[row, :] = Lf.vec.FV().NumPy()
            row += 1
    return A


def _seminorm(fes, fi, regularize):
    """SPD seminorm S (free block): l2 = identity, h1 = surface H1 stiffness."""
    if regularize == "l2":
        return np.eye(len(fi))
    from ngsolve import BilinearForm, grad, ds, TaskManager
    u, v = fes.TnT()
    a = BilinearForm(fes, symmetric=True)
    a += grad(u).Trace() * grad(v).Trace() * ds
    a += 1.0e-10 * u * v * ds
    with TaskManager():
        a.Assemble()
    Sd = np.array(a.mat.ToDense())
    return Sd[np.ix_(fi, fi)]


def _element_centroids(mesh):
    """Interior (non-vertex) measure points: element centroids of the eval mesh."""
    vpts = np.array([list(p.point) for p in mesh.vertices])
    cents = []
    for el in mesh.Elements():
        vs = [v.nr for v in el.vertices]
        cents.append(vpts[vs].mean(axis=0))
    return np.array(cents) if cents else vpts


def _peak_current_density(fes, coil, gfu):
    """max |K| = max |grad_s psi|, sampled via an L2 projection of the |grad|
    CF onto the boundary-H1 nodal values (the documented vertex-DOF workaround
    for boundary-defined H1 point eval)."""
    from ngsolve import GridFunction, grad, sqrt as ng_sqrt, InnerProduct
    Kmag = ng_sqrt(InnerProduct(grad(gfu).Trace(), grad(gfu).Trace()))
    gK = GridFunction(fes)
    gK.Set(Kmag, definedon=coil.Boundaries(".*"))
    return float(np.max(np.abs(gK.vec.FV().NumPy())))


def _build_problem(args):
    """Shared setup: coil FES, eval constraint + measure points, target, A, S."""
    from ngsolve import Mesh, H1, specialcf
    from radia.stream_function import aca_tsvd, RegularizedTSVD

    coil = Mesh(args.coil_vol)
    fes = H1(coil, order=args.order, definedon=coil.Boundaries(".*"))
    n = specialcf.normal(3)
    fi = np.where(np.array(fes.FreeDofs()))[0]

    target_cf, tdim = _parse_target_cf(args.target_cf)
    if tdim not in (1, 3):
        raise ValueError(f"target-cf dim {tdim}; expected 1 (Bz) or 3 (B)")
    vector_b = (tdim == 3)
    comps = [0, 1, 2] if vector_b else [2]

    evalm = Mesh(args.eval_vol)
    cpts = np.array([list(p.point) for p in evalm.vertices])    # constraints
    if args.eval_max and len(cpts) > args.eval_max:
        cpts = cpts[np.linspace(0, len(cpts) - 1, args.eval_max).astype(int)]
    mpts = _element_centroids(evalm)                            # interior measure
    if args.eval_max and len(mpts) > args.eval_max:
        mpts = mpts[np.linspace(0, len(mpts) - 1, args.eval_max).astype(int)]

    Bc = _target_values(target_cf, vector_b, evalm, cpts)
    Bm = _target_values(target_cf, vector_b, evalm, mpts)
    Ac = _assemble_biot_savart(fes, n, cpts, comps)
    Am = _assemble_biot_savart(fes, n, mpts, comps)
    S = _seminorm(fes, fi, args.regularize)
    Af = Ac[:, fi]
    base = aca_tsvd(Ac.shape[0], len(fi), lambda i, j: float(Af[i, j]),
                    modes=Ac.shape[0], kmax=min(Ac.shape[0], len(fi)),
                    aca_eps=1e-10, method=3)
    reg = RegularizedTSVD.from_stiffness(base, S)
    md = float(np.mean(np.diag(Af.T @ Af)))
    return dict(coil=coil, fes=fes, fi=fi, n=n, vector_b=vector_b,
                comps=comps, Ac=Ac, Af=Af, Am=Am, Bc=Bc, Bm=Bm, reg=reg,
                md=md, n_constraint=len(cpts), n_measure=len(mpts))


def _solve_and_metrics(P, alpha):
    """Solve psi at a given alpha; return (psi_full, fit_rms, homogeneity, peak)."""
    Bc, Bm, Af, Am, fi, reg, md = (P["Bc"], P["Bm"], P["Af"], P["Am"],
                                   P["fi"], P["reg"], P["md"])
    psi_f = reg.solve(Bc, alpha=alpha * md) if alpha > 0 else reg.solve(Bc)
    fit_rms = float(np.linalg.norm(Af @ psi_f - Bc) / (np.linalg.norm(Bc) + 1e-30))
    # true homogeneity = field error at the INTERIOR measure points (not
    # constraints) -> captures the between-constraint deviation
    homo = float(np.linalg.norm(Am[:, fi] @ psi_f - Bm)
                 / (np.linalg.norm(Bm) + 1e-30))
    return psi_f, fit_rms, homo


def run_design(args):
    from ngsolve import GridFunction, TaskManager
    t0 = time.perf_counter()
    with TaskManager():
        P = _build_problem(args)
        t1 = time.perf_counter()
        psi_f, fit_rms, homo = _solve_and_metrics(P, args.alpha)
        psi = np.zeros(P["fes"].ndof)
        psi[P["fi"]] = psi_f
        gfu = GridFunction(P["fes"])
        gfu.vec.FV().NumPy()[:] = psi
        peak_J = _peak_current_density(P["fes"], P["coil"], gfu)
        t2 = time.perf_counter()
        result = {
            "method": "design",
            "target_cf": args.target_cf,
            "target_kind": "vector_B" if P["vector_b"] else "Bz",
            "coil_vol": os.path.basename(args.coil_vol),
            "eval_vol": os.path.basename(args.eval_vol),
            "order": args.order, "regularize": args.regularize,
            "alpha": args.alpha,
            "ndof": int(P["fes"].ndof), "ndof_free": int(len(P["fi"])),
            "n_constraint": int(P["n_constraint"]),
            "n_measure": int(P["n_measure"]),
            "n_constraints": int(P["Ac"].shape[0]),
            "fit_residual_rms": fit_rms,
            "homogeneity_rms": homo,
            "rms": homo,                       # primary metric = true homogeneity
            "peak_J": peak_J,
            "t_mesh_s": round(t1 - t0, 3),
            "t_solve_s": round(t2 - t1, 3),
            "t_total_s": round(time.perf_counter() - t0, 3),
        }
        if args.msh_output:
            try:
                from radia.gmsh_post_export import GmshPostExport
                post = GmshPostExport(P["coil"], boundary=True)
                post.add_field("psi", psi[:P["coil"].nv], ncomp=1)
                post.write(args.msh_output)
                result["msh"] = args.msh_output
            except Exception as e:                 # noqa: BLE001
                result["msh_error"] = str(e)
    return result


def run_pareto(args):
    """alpha-sweep -> (homogeneity, peak current density) Pareto front."""
    from ngsolve import GridFunction, TaskManager
    if args.pareto_lever != "alpha":
        return {"error": f"pareto-lever '{args.pareto_lever}' interface defined "
                         f"but body is a follow-on (task #52); use 'alpha'"}
    t0 = time.perf_counter()
    with TaskManager():
        P = _build_problem(args)
        lams = np.concatenate([[0.0], np.logspace(
            np.log10(args.alpha_min), np.log10(args.alpha_max),
            max(1, args.n_alpha - 1))])
        front = []
        for lam in lams:
            psi_f, fit_rms, homo = _solve_and_metrics(P, lam)
            psi = np.zeros(P["fes"].ndof); psi[P["fi"]] = psi_f
            gfu = GridFunction(P["fes"]); gfu.vec.FV().NumPy()[:] = psi
            peak = _peak_current_density(P["fes"], P["coil"], gfu)
            front.append({"alpha": float(lam), "homogeneity_rms": homo,
                          "peak_J": peak})
    return {
        "method": "pareto", "pareto_lever": args.pareto_lever,
        "target_cf": args.target_cf,
        "target_kind": "vector_B" if P["vector_b"] else "Bz",
        "coil_vol": os.path.basename(args.coil_vol),
        "eval_vol": os.path.basename(args.eval_vol),
        "order": args.order, "regularize": args.regularize,
        "ndof": int(P["fes"].ndof), "n_measure": int(P["n_measure"]),
        "n_alpha": len(front), "front": front,
        "t_total_s": round(time.perf_counter() - t0, 3),
    }


def run_manufacture(args):
    """single-stroke chain + sheet-metal distortion + CAD + PEEC.  Interface
    defined; body is a follow-on (task #52)."""
    return {
        "method": "manufacture",
        "error": "manufacture mode interface is defined but the body is a "
                 "follow-on (task #52): single-stroke chain + sheet-metal "
                 "distortion + STEP + PEEC.",
        "nlevels": args.nlevels, "chain_method": args.chain_method,
        "distort": args.distort, "step_output": args.step_output,
        "peec": args.peec,
    }


def run(args):
    dispatch = {"design": run_design, "pareto": run_pareto,
                "manufacture": run_manufacture}
    return dispatch[args.method](args)


def build_argparser():
    ap = argparse.ArgumentParser(
        description="Stream-function coil design (design / pareto / manufacture).")
    # ---- shared inputs ----
    ap.add_argument("--coil-vol", required=True,
                    help="coil surface mesh (.vol, 2D-in-3D)")
    ap.add_argument("--eval-vol", required=True,
                    help="evaluation region mesh (.vol; surface or volume)")
    ap.add_argument("--target-cf", required=True,
                    help="target field CoefficientFunction expr of x,y,z "
                         "(scalar -> Bz; 3-vector -> B)")
    ap.add_argument("--method", choices=["design", "pareto", "manufacture"],
                    default="design")
    # ---- solver (design + pareto) ----
    ap.add_argument("--order", type=int, default=3, help="psi FE order")
    ap.add_argument("--regularize", choices=["l2", "h1"], default="h1",
                    help="seminorm: l2 (min |psi|) or h1 (min surface-current "
                         "energy)")
    ap.add_argument("--alpha", type=float, default=0.0,
                    help="Tikhonov weight (design; 0 = exact-fit min-seminorm)")
    ap.add_argument("--eval-max", type=int, default=400,
                    help="cap on constraint / measure points (subsample)")
    # ---- pareto mode ----
    ap.add_argument("--pareto-lever",
                    choices=["alpha", "linf", "geometry", "sheetmetal"],
                    default="alpha",
                    help="front-pushing lever (alpha = L-curve sweep; "
                         "linf/geometry/sheetmetal: interface defined, "
                         "follow-on bodies)")
    ap.add_argument("--alpha-min", type=float, default=1e-4,
                    help="pareto alpha-sweep lower bound (relative)")
    ap.add_argument("--alpha-max", type=float, default=3e1,
                    help="pareto alpha-sweep upper bound (relative)")
    ap.add_argument("--n-alpha", type=int, default=12,
                    help="pareto alpha-sweep number of points")
    # ---- manufacture mode ----
    ap.add_argument("--nlevels", type=int, default=12,
                    help="number of iso-contours = wire turns (manufacture)")
    ap.add_argument("--chain-method",
                    choices=["field_aware", "kuijpers", "lobe", "greedy"],
                    default="field_aware",
                    help="single-stroke chain method (manufacture)")
    ap.add_argument("--distort", action="store_true",
                    help="single-current sheet-metal wire distortion (manufacture)")
    ap.add_argument("--distort-comps", default="rsz",
                    help="sheet-metal distortion components (manufacture)")
    ap.add_argument("--step-output", default="",
                    help="STEP CAD output path (manufacture)")
    ap.add_argument("--peec", action="store_true",
                    help="compute PEEC coil inductance (manufacture)")
    # ---- output ----
    ap.add_argument("--msh-output", default="",
                    help="optional GMSH .msh of psi on the coil surface")
    ap.add_argument("--output", default="", help="JSON output file")
    return ap


def main():
    return calc_main(run, build_argparser())


if __name__ == "__main__":
    main()
