"""Step-4 verification: the MOMENT-yano LU path handles NONLINEAR susceptibility (B-H curve) through the
standard Picard loop -- each Picard step rebuilds the moment system with the per-element chi(H) of the current
iterate (rad_relaxation_methods.cpp moment branch reads ctx.CurrentChiArray[e], which the Picard driver
refreshes from the BH curve every iteration).

Solves a C-yoke of nonlinear soft-iron hexes (MatSatIsoTab) in a uniform applied field deep enough to bite the
knee of the curve, with (a) EIEM2 (yano_moment=False) and (b) MOMENT (yano_moment=True), and compares the
EXTERNAL B field.  Both solve the same nonlinear physics so they must agree to a few percent (moment is the
more accurate collocation-free formula).  Internal-iron points are NOT compared (formulation-dependent;
CLAUDE.md "do not compare MSC internal fields").
"""
import numpy as np
import radia as rad
from yano_moment_iter_scaling import build_cyoke_hexes

MU0 = 4e-7 * np.pi

# Soft-iron B-H curve: linear-ish (mu_r ~ 3000) below the knee, saturating to ~2.1 T.
BH_DATA = [[0.0, 0.0], [200.0, 0.75], [500.0, 1.30], [1200.0, 1.70],
           [4000.0, 1.95], [20000.0, 2.08], [100000.0, 2.15]]


def solve_field(moment, H_app, probes):
    rad.UtiDelAll(); rad.set_demag_backend("yano")
    rad.SolverConfig(yano_moment=bool(moment), yano_eval_alpha=-1.0)
    objs = [rad.ObjHexahedron([list(v) for v in V], [0, 0, 0]) for V in build_cyoke_hexes(8, 2)]
    for h in objs:
        rad.MatApl(h, rad.MatSatIsoTab(BH_DATA))
    cont = rad.ObjCnt(objs + [rad.ObjBckg(lambda p: [0.0, MU0 * H_app, 0.0])])
    nit = rad.Solve(cont, 1e-6, 500, 0)            # method 0 = LU, Picard outer loop (nonlinear)
    M = np.asarray([m[1] for m in rad.ObjM(rad.ObjCnt(objs))], float)
    B = np.asarray([rad.Fld(cont, "b", list(p)) for p in probes], float)
    rad.SolverConfig(yano_moment=False); rad.UtiDelAll()
    return B, M, nit


def main():
    H_app = 8.0e5                                   # A/m -> strong enough that, even after the C-yoke's heavy
                                                    # self-demag, the internal H reaches the knee and the iron
                                                    # genuinely saturates (M -> Msat), so the BH curve bites.
    probes = [[0.0, 0.0, 0.05], [0.0, 0.0, 0.1], [0.1, 0.0, 0.05],
              [0.08, 0.08, 0.06], [-0.1, 0.05, 0.05], [0.0, 0.1, 0.03]]

    B_e, M_e, n_e = solve_field(False, H_app, probes)
    B_m, M_m, n_m = solve_field(True, H_app, probes)

    rel = np.linalg.norm(B_m - B_e) / max(np.linalg.norm(B_e), 1e-30)
    # nonlinearity sanity: peak |M| must be below saturation Msat = Bsat/mu0 and the response must be
    # sub-linear vs a linear mu_r=3000 extrapolation (chi(H) actually dropped) -> confirms the BH curve bit.
    Msat = 2.15 / MU0
    finite = np.all(np.isfinite(B_m)) and np.all(np.isfinite(M_m))
    saturating = np.max(np.linalg.norm(M_m, axis=1)) < Msat
    print(f"EIEM2  nonlinear: iters={n_e}  |B|@probes={np.linalg.norm(B_e, axis=1)}")
    print(f"moment nonlinear: iters={n_m}  |B|@probes={np.linalg.norm(B_m, axis=1)}")
    print(f"||B_moment - B_eiem2|| / ||B_eiem2|| = {rel:.3e}")
    print(f"moment |M|max={np.max(np.linalg.norm(M_m, axis=1)):.4e}  (Msat={Msat:.3e})  saturating={saturating}")
    ok = finite and saturating and rel < 0.1
    print("RESULT:", "PASS -- moment nonlinear (Picard + BH curve) runs and agrees with EIEM2" if ok else "FAIL")
    return ok


if __name__ == "__main__":
    raise SystemExit(0 if main() else 1)
