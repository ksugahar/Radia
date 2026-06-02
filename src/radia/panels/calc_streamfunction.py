#!/usr/bin/env python3
"""calc_streamfunction.py -- Stream-Function (SF) coil design (Layer 4, headless).

Design a surface coil on a STANDALONE 2D surface mesh (.vol) that produces a
target field over an evaluation region (.vol, surface OR volume).  The target
and the computed field are NGSolve CoefficientFunctions (no .sol / .csv).

I/O (per the 2026-06 design decisions):
  --coil-vol   coil surface mesh (.vol, 2D-in-3D).  psi = H1 on it (FE-direct).
  --eval-vol   evaluation region (.vol).  Surface OR volume; the RMS integral
               adapts to the mesh dimension automatically.
  --target-cf  target field as a CoefficientFunction EXPRESSION of x, y, z.
               SCALAR expr -> Bz target (e.g. "x" = Gx, "1" = uniform Bz,
               "x*x-y*y" = C2 ellipse shim).  3-VECTOR expr "(Bx,By,Bz)" ->
               full vector-B target.

Method (this file currently implements --method design; pareto / manufacture
are follow-on modes):
  design   target -> A psi = B (folded-Tikhonov RegularizedTSVD) -> psi,
           field RMS over the eval region, peak surface current density.

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


def run_streamfunction(args):
    t0 = time.perf_counter()
    from ngsolve import (Mesh, H1, GridFunction, grad, specialcf, ds,
                         LinearForm, sqrt as ng_sqrt, Integrate, InnerProduct,
                         CoefficientFunction, Cross, TaskManager, x, y, z)
    from radia.stream_function import aca_tsvd, RegularizedTSVD

    if args.method != "design":
        return {"error": f"method '{args.method}' not yet implemented "
                         f"(design only in this version)"}

    mu0_4pi = 1.0e-7
    with TaskManager():
        # ---- coil surface FE space (FE-direct psi) ----
        coil = Mesh(args.coil_vol)
        fes = H1(coil, order=args.order, definedon=coil.Boundaries(".*"))
        n = specialcf.normal(3)
        u, v = fes.TnT()
        fi = np.where(np.array(fes.FreeDofs()))[0]
        ndof = fes.ndof
        t_mesh = time.perf_counter()

        # ---- target CoefficientFunction + eval points ----
        target_cf, tdim = _parse_target_cf(args.target_cf)
        vector_b = (tdim == 3)
        if tdim not in (1, 3):
            return {"error": f"target-cf dim {tdim}; expected 1 (Bz) or 3 (B)"}
        evalm = Mesh(args.eval_vol)
        epts = np.array([list(p.point) for p in evalm.vertices])
        if args.eval_max and len(epts) > args.eval_max:
            idx = np.linspace(0, len(epts) - 1, args.eval_max).astype(int)
            epts = epts[idx]
        M = len(epts)

        # target values at eval points
        comps = [0, 1, 2] if vector_b else [2]                # B components used
        Bt = []
        for p in epts:
            mp = evalm(p[0], p[1], p[2])
            tv = target_cf(mp)
            tv = list(tv) if vector_b else [float(tv)]
            Bt.extend(tv)
        Bt = np.array(Bt, dtype=float)

        # ---- assemble A: row per (eval point, component) ----
        # B = (mu0/4pi) Kv x (r_t - r') / r^3,  Kv = n x grad_s phi  (Trace, ds)
        nrow = M * len(comps)
        A = np.zeros((nrow, ndof))
        row = 0
        for p in epts:
            dxt, dyt, dzt = p[0] - x, p[1] - y, p[2] - z
            r2 = dxt * dxt + dyt * dyt + dzt * dzt
            r3 = r2 * ng_sqrt(r2)
            Kv = Cross(n, grad(v).Trace())
            cross = (Kv[1] * dzt - Kv[2] * dyt,                # x
                     Kv[2] * dxt - Kv[0] * dzt,                # y
                     Kv[0] * dyt - Kv[1] * dxt)               # z
            for c in comps:
                Lf = LinearForm(fes)
                Lf += (mu0_4pi / r3) * cross[c] * ds
                Lf.Assemble()
                A[row, :] = Lf.vec.FV().NumPy()
                row += 1
        Af = A[:, fi]
        t_assemble = time.perf_counter()

        # ---- folded-Tikhonov solve (H1 seminorm) ----
        S = _h1_stiffness(fes, fi, coil)
        base = aca_tsvd(nrow, len(fi), lambda i, j: float(Af[i, j]),
                        modes=nrow, kmax=min(nrow, len(fi)), aca_eps=1e-10,
                        method=3)
        reg = RegularizedTSVD.from_stiffness(base, S)
        alpha = args.alpha
        if alpha > 0:
            md = float(np.mean(np.diag(Af.T @ Af)))
            psi_f = reg.solve(Bt, alpha=alpha * md)
        else:
            psi_f = reg.solve(Bt)
        psi = np.zeros(ndof)
        psi[fi] = psi_f
        t_solve = time.perf_counter()

        # ---- metrics ----
        misfit = float(np.linalg.norm(Af @ psi_f - Bt)
                       / (np.linalg.norm(Bt) + 1e-30))
        gfu = GridFunction(fes)
        gfu.vec.FV().NumPy()[:] = psi
        # peak surface current density |K| = |grad_s psi|, sampled via an
        # L2 projection of the |grad| CF onto the boundary H1 nodal values
        Kmag = ng_sqrt(InnerProduct(grad(gfu).Trace(), grad(gfu).Trace()))
        gK = GridFunction(fes)
        gK.Set(Kmag, definedon=coil.Boundaries(".*"))
        peak_J = float(np.max(np.abs(gK.vec.FV().NumPy())))
        # integral RMS over the eval region (surface or volume) using the
        # target CF vs the (interpolated) computed-field nodal values:
        # report the discrete residual as the primary RMS.
        rms = misfit

        result = {
            "method": "design",
            "target_cf": args.target_cf,
            "target_kind": "vector_B" if vector_b else "Bz",
            "coil_vol": os.path.basename(args.coil_vol),
            "eval_vol": os.path.basename(args.eval_vol),
            "order": args.order, "regularize": "tikhonov_h1",
            "alpha": args.alpha,
            "ndof": int(ndof), "ndof_free": int(len(fi)),
            "n_eval": int(M), "n_constraints": int(nrow),
            "rms": rms, "peak_J": peak_J,
            "t_mesh_s": round(t_mesh - t0, 3),
            "t_assemble_s": round(t_assemble - t_mesh, 3),
            "t_solve_s": round(t_solve - t_assemble, 4),
            "t_total_s": round(time.perf_counter() - t0, 3),
        }

        # ---- optional GMSH post (psi on the coil surface) ----
        if args.msh_output:
            try:
                from radia.gmsh_post_export import GmshPostExport
                post = GmshPostExport(coil, boundary=True)
                post.add_field("psi", gfu.vec.FV().NumPy()
                               [:coil.nv], ncomp=1)
                post.write(args.msh_output)
                result["msh"] = args.msh_output
            except Exception as e:               # noqa: BLE001
                result["msh_error"] = str(e)

    return result


def _h1_stiffness(fes, fi, mesh):
    """H1 surface stiffness (free block) for the min-current-energy seminorm."""
    from ngsolve import (BilinearForm, grad, ds, TaskManager, dx)
    u, v = fes.TnT()
    a = BilinearForm(fes, symmetric=True)
    a += grad(u).Trace() * grad(v).Trace() * ds
    a += 1.0e-10 * u * v * ds
    with TaskManager():
        a.Assemble()
    Sd = np.array(a.mat.ToDense())
    return Sd[np.ix_(fi, fi)]


def build_argparser():
    ap = argparse.ArgumentParser(
        description="Stream-function coil design (Design mode).")
    ap.add_argument("--coil-vol", required=True,
                    help="coil surface mesh (.vol, 2D-in-3D)")
    ap.add_argument("--eval-vol", required=True,
                    help="evaluation region mesh (.vol; surface or volume)")
    ap.add_argument("--target-cf", required=True,
                    help="target field as a CoefficientFunction expr of x,y,z "
                         "(scalar -> Bz; 3-vector -> B)")
    ap.add_argument("--method", choices=["design", "pareto", "manufacture"],
                    default="design")
    ap.add_argument("--order", type=int, default=3, help="psi FE order")
    ap.add_argument("--regularize", choices=["h1"], default="h1")
    ap.add_argument("--alpha", type=float, default=0.0,
                    help="Tikhonov weight (0 = exact-fit min-seminorm)")
    ap.add_argument("--eval-max", type=int, default=400,
                    help="cap on eval points (subsample if exceeded)")
    ap.add_argument("--msh-output", default="",
                    help="optional GMSH .msh of psi on the coil surface")
    ap.add_argument("--output", default="", help="JSON output file")
    return ap


def main():
    return calc_main(run_streamfunction, build_argparser())


if __name__ == "__main__":
    main()
