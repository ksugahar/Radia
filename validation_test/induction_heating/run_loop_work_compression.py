"""Galerkin P0 and bordered body compression on a self-authored ring.

Run P0 dense/fmm or body dense/hacapk separately with identical mesh and
thread controls. HACApK retains the existing construction guard. The full
case set measures two carriers, scalar/jump impedance and two smooth lift
changes. --compact-cases keeps one carrier and the unperturbed lift for cost
runs. Fields are saved beside the path-free JSON, not committed by this script.
This measures implementation error at the same mesh, not physical accuracy.
The ACA-dependent agreement envelope is a validation target, not a theorem
converting local ACA truncation tolerance into global field error.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
import time
from pathlib import Path

import ngsolve as ng
import numpy as np
from radia import bem_loop_extension as extension
from radia.bem_sibc_solver import ScalarBIESIBCSolver
from radia.surface_impedance import PanelSurfaceImpedance
from run_loop_work_ring import ring_mesh
from run_loop_work_performance import memory


def run(level, backend, compact, repeats, threads, body_backend="dense", aca_eps=1e-10, gmres_tol=1e-10):
    mu0, omega = 4e-7*np.pi, 2*np.pi*5e4
    z = (1+1j)*np.sqrt(omega*mu0/(2*5.8e7))

    def incident(points):
        return .5*mu0*np.c_[-points[:, 1], points[:, 0], np.zeros(len(points))]

    records, fields, phases = [], {}, []
    for repeat in range(repeats):
        started = time.perf_counter()
        def phase(_tag,message):
            if message.startswith("HACApK compress:") or message.startswith("HACApK compress done") or message=="HACApK release source entries":
                phases.append(dict(repeat=repeat,phase=("before-compression" if ":" in message else "before-entry-release" if message=="HACApK release source entries" else "after-entry-release"),**memory()))
        with ng.TaskManager():
            mesh, points, triangles = ring_mesh(level)
            solver = ScalarBIESIBCSolver(mesh, order=1, assemble_dense=body_backend=="dense",
                use_intree_hacapk=body_backend=="hacapk", hacapk_aca_eps=aca_eps,
                use_intree_bem=True, intree_geom_order=1,
                intree_singular_n_q=6, intree_regular_quad_degree=7,
                loop_work_backend=backend,log_fn=phase)
            body_seconds = time.perf_counter()-started
            print("BIE assembled", len(points), "vertices", body_seconds, flush=True)
            original_geometry = extension._prepare_loop_geometry
            try:
                for offset in ([0.] if compact else [0., .00015]):
                    angle = np.arange(96)*2*np.pi/96
                    carrier = (.03+offset)*np.c_[np.cos(angle), np.sin(angle), np.zeros(96)]
                    for law in ["scalar", "jump"]:
                        impedance = z if law == "scalar" else PanelSurfaceImpedance(
                            np.where(points[triangles].mean(axis=1)[:, 0] >= 0, 3*z, .25*z))
                        baseline = None
                        for lift in ([0] if compact else [0, 1, 2]):
                            def geometry(*args, **kwargs):
                                result = original_geometry(*args, **kwargs).copy()
                                if lift:
                                    w = (.05*points[:, 0]/.03 if lift == 1 else .03*points[:, 2]/.003)
                                    mass = np.asarray(solver.M.sum(axis=1)).reshape(-1)
                                    w -= mass @ w/mass.sum()
                                    result["Theta"] = result["Theta"]+w[result["Tmap"]]
                                    if body_backend == "hacapk":
                                        from scipy.sparse.linalg import gmres
                                        from radia.bem_loop_work import _body_layer
                                        sl,dl = _body_layer(solver,"SL"),_body_layer(solver,"DL")
                                        rhs = .5*(solver.M @ w)-dl @ w
                                        correction,info = gmres(sl,rhs,rtol=1e-12,atol=0.,restart=100,maxiter=500)
                                        if info or np.linalg.norm(sl @ correction-rhs)/np.linalg.norm(rhs)>1e-6:
                                            raise RuntimeError("Lift-transform SL residual failed")
                                    else:
                                        correction=np.linalg.solve(solver.SL,(.5*solver.M-solver.DL) @ w)
                                    result["qT"] = result["qT"]-correction
                                return result

                            extension._prepare_loop_geometry = geometry
                            print("Loop case", repeat, offset, law, lift, flush=True)
                            case_started = time.perf_counter()
                            out = extension.solve_loop_extended(solver, -points[:, 2].astype(complex),
                                impedance, omega, incident, section_anchor=(.03, 0), carrier_ring=carrier,
                                _reuse_prepared=compact, hacapk=body_backend=="hacapk",
                                gmres=dict(tol=gmres_tol,restart=100,maxiter=500))
                            if baseline is None:
                                baseline = out
                            key = f"r{repeat}-o{offset}-{law}-g{lift}"
                            fields[key] = out["H_t_tri"]
                            record = dict(key=key, repeat=repeat, carrier_offset_m=offset,
                                impedance_law=law, lift=lift, alpha=[out["alpha"].real, out["alpha"].imag],
                                heat_power_w=out["P_total"], reaction_power_w=out["P_reaction"],
                                frozen_heat_power_w=out["P_frozen"],
                                power_balance_relative=abs(out["P_reaction"]/out["P_total"]-1),
                                linear_residual=out["linear_residual_rel"],
                                faraday_residual=out["faraday_residual_rel"],
                                seam_jump_deviation=out["theta_jump_max_deviation"],
                                lift_field_relative=np.linalg.norm(out["H_t_tri"]-baseline["H_t_tri"])/np.linalg.norm(baseline["H_t_tri"]),
                                lift_heat_relative=abs(out["P_total"]/baseline["P_total"]-1),
                                seconds=time.perf_counter()-case_started,
                                work=out["loop_work_diagnostics"],
                                gmres_iterations=out["loop_gmres_iterations"],
                                gmres_controls=dict(tol=gmres_tol,restart=100,maxiter=500) if body_backend=="hacapk" else None,
                                frozen_residual=out["loop_frozen_residual_rel"], **memory())
                            if record["linear_residual"] > 1e-6 or record["faraday_residual"] > 1e-6:
                                raise RuntimeError("Loop solution residual gate failed")
                            if record["power_balance_relative"] > .1:
                                raise RuntimeError("Reaction/heat gate failed")
                            records.append(record)
                            print(key, record["heat_power_w"], record["seconds"], flush=True)
            finally:
                extension._prepare_loop_geometry = original_geometry
        for row in records:
            if row["repeat"] == repeat:
                row["body_assembly_seconds"] = body_seconds
    return dict(schema="radia.loop-work-compression.v1", level=level,
        vertices=len(points), faces=len(triangles), backend=backend, platform_class=platform.system(),
        body_backend=body_backend, aca_eps=aca_eps if body_backend=="hacapk" else None,
        body_hmatrix_controls=(solver.hacapk_controls.copy() if body_backend=="hacapk" else None),
        agreement_bound=max(10*np.sqrt(aca_eps),1e-10) if body_backend=="hacapk" else .001,
        body_hmatrix_stats=([solver._SL_hacapk.GetStats(),solver._DL_hacapk.GetStats()]
                            if body_backend=="hacapk" else None),
        threads=threads, ngsolve_version=ng.__version__,
        geometry_sha256=hashlib.sha256(points.tobytes()+triangles.tobytes()).hexdigest(),
        model=dict(major_radius_m=.03, minor_radius_m=.003, frequency_hz=5e4,
                   conductivity_s_per_m=5.8e7, axial_field_a_per_m=1),
        records=records, body_storage_phases=phases, **memory()), fields


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--level", type=int, required=True)
    parser.add_argument("--backend", choices=["auto", "dense", "fmm"], required=True)
    parser.add_argument("--body-backend",choices=["dense","hacapk"],default="dense")
    parser.add_argument("--aca-eps",type=float,default=1e-10)
    parser.add_argument("--gmres-tol",type=float,default=1e-10)
    parser.add_argument("--compact-cases", action="store_true")
    parser.add_argument("--repeats", type=int, default=1)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--reference", type=Path)
    args = parser.parse_args()
    if args.threads < 1 or args.repeats < 1:
        parser.error("threads and repeats must be positive")
    ng.SetNumThreads(args.threads)
    result, fields = run(args.level, args.backend, args.compact_cases, args.repeats, args.threads,args.body_backend,args.aca_eps,args.gmres_tol)
    if args.reference:
        reference = json.loads(args.reference.read_text(encoding="utf-8"))
        if (reference["geometry_sha256"] != result["geometry_sha256"]
                or reference["threads"] != result["threads"]):
            raise ValueError("Reference mesh/thread controls differ")
        reference_fields = np.load(args.reference.with_suffix(".npz"))
        by_key = {row["key"]: row for row in reference["records"]}
        for row in result["records"]:
            prior = by_key[row["key"]]
            row["reference_heat_relative"] = abs(row["heat_power_w"]/prior["heat_power_w"]-1)
            row["reference_alpha_relative"] = abs(complex(*row["alpha"])/complex(*prior["alpha"])-1)
            row["reference_reaction_relative"] = abs(row["reaction_power_w"]/prior["reaction_power_w"]-1)
            row["reference_field_relative"] = np.linalg.norm(fields[row["key"]]-reference_fields[row["key"]])/np.linalg.norm(reference_fields[row["key"]])
            if any(row[key] > result["agreement_bound"]
                   for key in ("reference_heat_relative","reference_alpha_relative","reference_reaction_relative")):
                raise RuntimeError("Compression heat/alpha/reaction agreement bound exceeded")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    np.savez(args.output.with_suffix(".npz"), **fields)
    args.output.write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8")


if __name__ == "__main__":
    main()
