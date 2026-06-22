"""Step-2 verification: rad.Solve(backend='yano') with SolverConfig(yano_moment=True) runs the parameter-free
MOMENT system end-to-end (assemble -> LU -> write M) and gives a physically correct field.

Solves the same C-yoke soft-iron block in a uniform applied field with (a) EIEM2 (yano_moment=False) and
(b) MOMENT (yano_moment=True), and compares the external B field.  Both solve the same physics, so they must
agree to a few percent (moment is the more accurate one).  Also sanity-checks the moment M is finite and
nonzero.  (Temp scratch; promote once green.)
"""
import numpy as np
import radia as rad
from yano_moment_iter_scaling import build_cyoke_hexes

MU0 = 4e-7 * np.pi


def solve_field(hexes, mu_r, moment, probes):
    rad.UtiDelAll(); rad.set_demag_backend("yano")
    rad.SolverConfig(yano_moment=bool(moment), yano_eval_alpha=-1.0)
    objs = [rad.ObjHexahedron([list(v) for v in V], [0, 0, 0]) for V in hexes]
    for h in objs:
        rad.MatApl(h, rad.MatLin(mu_r))
    # uniform applied H = (0, 1e3, 0) A/m  ->  B_background = mu0 * H
    cont = rad.ObjCnt(objs + [rad.ObjBckg(lambda p: [0.0, MU0 * 1e3, 0.0])])
    rad.Solve(cont, 1e-6, 200, 0)                 # method 0 = LU
    M = np.asarray([m[1] for m in rad.ObjM(rad.ObjCnt(objs))], float)
    B = np.asarray([rad.Fld(cont, "b", list(p)) for p in probes], float)
    rad.SolverConfig(yano_moment=False)
    rad.UtiDelAll()
    return B, M


def main():
    nxy, nz, mu_r = 8, 2, 1000.0
    hexes = build_cyoke_hexes(nxy, nz)
    # EXTERNAL probes only (clearly outside the yoke bbox |x|,|y|<=0.06, |z|<=0.02); internal-iron points are
    # formulation-dependent (CLAUDE.md: do not compare MSC internal fields).
    probes = [[0.0, 0.0, 0.05], [0.0, 0.0, 0.1], [0.1, 0.0, 0.05],
              [0.08, 0.08, 0.06], [-0.1, 0.05, 0.05], [0.0, 0.1, 0.03]]

    B_e, M_e = solve_field(hexes, mu_r, False, probes)
    B_m, M_m = solve_field(hexes, mu_r, True, probes)

    rel = np.linalg.norm(B_m - B_e) / max(np.linalg.norm(B_e), 1e-30)
    finite = np.all(np.isfinite(B_m)) and np.all(np.isfinite(M_m))
    nonzero = np.linalg.norm(M_m) > 1e-6
    print(f"EIEM2  |B| at probes: {np.linalg.norm(B_e, axis=1)}")
    print(f"moment |B| at probes: {np.linalg.norm(B_m, axis=1)}")
    print(f"||B_moment - B_eiem2|| / ||B_eiem2|| = {rel:.3e}")
    print(f"moment M: finite={finite}, ||M||={np.linalg.norm(M_m):.3e} (EIEM2 ||M||={np.linalg.norm(M_e):.3e})")
    ok = finite and nonzero and rel < 0.1          # same physics -> within ~10%
    print("RESULT:", "PASS -- moment solve runs end-to-end and agrees with EIEM2" if ok else "FAIL")
    return ok


if __name__ == "__main__":
    raise SystemExit(0 if main() else 1)
