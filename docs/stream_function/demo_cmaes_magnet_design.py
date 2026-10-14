"""CMA-ES magnetization-direction design (the nonlinear counterpart of TSVD).

The (ACA+)+TSVD solver in this folder handles the LINEAR design problem:
magnetization *amplitudes* with fixed directions (see demo_magnet_array.py).
When the design variables enter the field NONLINEARLY -- e.g. the magnetization
*directions* (angles) -- the problem is no longer a linear least-norm solve and
a black-box optimizer is the right tool.  This demo uses CMA-ES (Optuna's
`CmaEsSampler`, not re-implemented) to tune the in-plane magnetization angle of
every pixel of a fixed magnet array so the field over a target region is a
uniform transverse field.

This is the "+ CMA-ES" half of the "ACA stream function + CMA-ES" workflow
(IEEJ SA-25-020): the linear inner design uses (ACA+)+TSVD; the nonlinear outer
design uses CMA-ES.  Forward field is Radia's fixed-magnet field (`ObjRecMag` + `radia.Fld`).

Run:  python demo_cmaes_magnet_design.py [--trials N]
"""
import argparse
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent

import radia as rad

# ----------------------------------------------------------------------------
# Fixed geometry: planar array of magnet pixels with FREE in-plane magnetization
# angles; target = uniform transverse (Bx) field over a region above the array.
# ----------------------------------------------------------------------------
N_SIDE = 4                 # N = 16 pixels -> 16 angle design variables
N = N_SIDE * N_SIDE
ARRAY_HALF = 0.04          # array half-width [m]
CUBE = 0.012               # pixel size [m]
M0 = 9.5e5                 # fixed remanent magnetisation magnitude [A/m]
TARGET_Z = 0.03            # target plane height [m]
TARGET_HALF = 0.015        # target region half-width [m]

_xs = np.linspace(-ARRAY_HALF, ARRAY_HALF, N_SIDE)
_gx, _gy = np.meshgrid(_xs, _xs, indexing="ij")
POS = np.column_stack([_gx.ravel(), _gy.ravel(), np.zeros(N)])

_ts = np.linspace(-TARGET_HALF, TARGET_HALF, 5)
_tx, _ty = np.meshgrid(_ts, _ts, indexing="ij")
TARGET_PTS = np.column_stack([_tx.ravel(), _ty.ravel(),
                              np.full(_tx.size, TARGET_Z)])


def field_at_target(angles):
    """Bx,By,Bz over the target region for per-pixel in-plane angles [rad]."""
    rad.UtiDelAll()
    handles = []
    for k in range(N):
        th = float(angles[k])
        mag = [M0 * np.cos(th), M0 * np.sin(th), 0.0]
        handles.append(rad.magnet_box(POS[k].tolist(), [CUBE, CUBE, CUBE], mag))
    cont = rad.ObjCnt(handles)
    B = np.asarray(rad.Fld(cont, "b", TARGET_PTS.tolist()),
                   dtype=float).reshape(-1, 3)
    rad.UtiDelAll()
    return B


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--trials", type=int, default=200)
    args = ap.parse_args()

    try:
        import optuna
        from optuna.samplers import CmaEsSampler
    except ImportError:
        print("optuna not installed -- 'pip install optuna' to run this demo.")
        return
    optuna.logging.set_verbosity(optuna.logging.WARNING)

    # Reference field scale: all pixels magnetised along +x.
    B_aligned = field_at_target(np.zeros(N))
    B0 = 0.5 * float(np.mean(B_aligned[:, 0]))   # reachable uniform-Bx target [T]
    target = np.array([B0, 0.0, 0.0])

    def cost(angles):
        B = field_at_target(angles)
        return float(np.sqrt(np.mean(np.sum((B - target) ** 2, axis=1))) / abs(B0))

    def objective(trial):
        angles = np.array([trial.suggest_float(f"theta_{k}", -np.pi, np.pi)
                           for k in range(N)])
        return cost(angles)

    history = []
    study = optuna.create_study(direction="minimize",
                                sampler=CmaEsSampler(seed=42))
    study.optimize(objective, n_trials=args.trials,
                   callbacks=[lambda st, tr: history.append(
                       (tr.number, st.best_value))])

    best_angles = np.array([study.best_params[f"theta_{k}"] for k in range(N)])
    B_best = field_at_target(best_angles)
    first = study.trials[0].value
    best = study.best_value

    print("=== CMA-ES magnetization-direction design (uniform transverse field) ===")
    print(f"N pixels = {N}  ->  {N} angle design variables (CMA-ES, continuous)")
    print(f"target = uniform Bx = {B0*1e3:.3f} mT over a "
          f"{2*TARGET_HALF*1e3:.0f} mm region at z = {TARGET_Z*1e3:.0f} mm")
    print()
    print(f"objective (rel. RMS deviation from target field):")
    print(f"  first CMA-ES trial = {first:.4e}")
    print(f"  best  CMA-ES trial = {best:.4e}   ({first / best:.1f}x better)")
    print(f"achieved over target region:")
    print(f"  <Bx> = {np.mean(B_best[:,0])*1e3:.4f} mT  "
          f"(ripple pk-pk = {(B_best[:,0].max()-B_best[:,0].min())*1e3:.4f} mT)")
    print(f"  <|By|> = {np.mean(np.abs(B_best[:,1]))*1e3:.4f} mT, "
          f"<|Bz|> = {np.mean(np.abs(B_best[:,2]))*1e3:.4f} mT  (want ~0)")
    print(f"CMA-ES trials = {args.trials}")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        plt.rcParams["pdf.fonttype"] = 42
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.plot([h[0] for h in history], [h[1] for h in history], "o-", ms=2)
        ax.set_xlabel("trial"); ax.set_ylabel("best objective so far")
        ax.set_yscale("log")
        ax.grid(True, which="both", alpha=0.3)
        out = HERE / "demo_cmaes_magnet_design.png"
        fig.tight_layout(); fig.savefig(out, dpi=110)
        print(f"\nSaved plot: {out.name}")
    except ImportError:
        print("\n(matplotlib not installed; skipped plot)")


if __name__ == "__main__":
    main()
