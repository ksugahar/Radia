"""Coarse direct-solver migration check; not acceptance of the p3 ROM campaign."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform

os.environ.update(TEAM_CASE="team10", TEAM_QUARTER="1", TEAM_COIL2="1")

import numpy as np
import ngsolve as ng
from aphi_model import ModelAphi, run_fom
from rom_aphi import ExtenderAphi, bits
from transient import waveform


def residual(mat, inverse, rhs, indices):
    sol = rhs.CreateVector()
    sol.data = inverse * rhs
    defect = rhs.CreateVector()
    defect.data = mat * sol - rhs
    r = defect.FV().NumPy()[indices]
    b = rhs.FV().NumPy()[indices]
    return float(np.linalg.norm(r) / max(np.linalg.norm(b), 1e-300))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--maxh", type=float, default=0.02)
    args = parser.parse_args()
    ng.SetNumThreads(2)
    with ng.TaskManager():
        model = ModelAphi(args.maxh, order=1, dt=0.002)
        ext = ExtenderAphi(model)
        free = np.flatnonzero(np.asarray(list(model.fes.FreeDofs()), dtype=bool))
        model.gfu.vec[:] = 0
        model.a.AssembleLinearization(model.gfu.vec)
        rhs = model.f.vec.CreateVector()
        rhs.data = model.fs[0].vec + model.fs[1].vec
        checks = {
            "tangent": residual(model.a.mat, model.a.mat.Inverse(
                model.fes.FreeDofs(), inverse="sparsecholesky"), rhs, free),
            "source_air": residual(model.k_rest.mat, ext.air_inv, rhs, model.idx_air),
        }
        trace = model.gfu.vec.CreateVector()
        trace[:] = 0
        trace.FV().NumPy()[model.idx_gamma] = np.random.default_rng(19).normal(
            size=len(model.idx_gamma))
        rhs.data = model.k_rest.mat * trace
        checks["trace_air"] = residual(model.k_rest.mat, ext.air_inv, rhs, model.idx_air)
        rhs.data = ext.k_ref.mat * trace
        checks["trace_interior"] = residual(ext.k_ref.mat, ext.int_inv, rhs, model.idx_int)
        steel = np.sort(np.r_[model.idx_gamma, model.idx_int])
        checks["indicator_steel"] = residual(ext.k_ref.mat, ext.k_ref.mat.Inverse(
            bits(model.fes.ndof, steel), inverse="sparsecholesky"), rhs, steel)
        snapshots, history, info = run_fom(model, waveform("mix_add"), 3)
    ok = (all(np.isfinite(v) and v <= 1e-7 for v in checks.values())
          and np.isfinite(snapshots).all() and len(history) == 3)
    root = Path(__file__).parent
    result = dict(schema="radia.surface-rom.direct-migration.v1", ok=bool(ok),
                  scope="coarse residual and three-step nonlinear smoke only",
                  inverse="sparsecholesky", maxh=args.maxh, order=1,
                  ndof=model.fes.ndof, relative_true_residuals=checks,
                  residual_limit=1e-7, history=history, platform_class=platform.system(),
                  python=platform.python_version(), ngsolve=ng.__version__,
                  linear_solves=info["linear_solves"],
                  sources={p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                           for p in sorted(root.glob("*.py"))})
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result, indent=2), flush=True)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
