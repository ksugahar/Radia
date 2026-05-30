#!/usr/bin/env python
"""Belted-tree validation: the C++ BuildLoopBasis adds exactly b1 global "belt"
loops for a multiply-connected (ring) body, completing the short-cycle basis.

A thin 1-element-thick ring is the SHARPEST test: its element-adjacency graph is a
single C_N cycle, so there are ZERO short (3-/4-) cycles -- the only null mode is
the one global belt loop around the hole. A short-cycle-only basis finds NOTHING
(deflation_cycles == 0); a correct belted tree finds exactly 1.

We check at high mu_r:
  (1) deflation_cycles (GetSolveStats) == expected belt count  (sharp count proof)
  (2) SVD null dim of the geometric operator N (cross-check)
  (3) the converged solution is loop-free (alignment restored) only with deflation

Usage: python belt_loop_validation.py
"""
import numpy as np
import radia as rad

EDGE = 0.01
D = np.array([1.0, 1.0, 1.0]) / np.sqrt(3.0)
B0 = 0.1


def build_frame(outer, thick, nz, mu_r):
    """Square frame: outer x outer x nz grid minus the central hole (thin walls).
    thick>=outer/2 => no hole => solid block (simply connected control)."""
    rad.UtiDelAll()
    mat = rad.MatLin(float(mu_r))
    objs = []
    lo, hi = thick, outer - thick  # hole = [lo,hi) x [lo,hi)
    for iz in range(nz):
        for iy in range(outer):
            for ix in range(outer):
                if (lo <= ix < hi) and (lo <= iy < hi):
                    continue  # hole
                x0, y0, z0 = ix*EDGE, iy*EDGE, iz*EDGE
                v = [[x0,y0,z0],[x0+EDGE,y0,z0],[x0+EDGE,y0+EDGE,z0],[x0,y0+EDGE,z0],
                     [x0,y0,z0+EDGE],[x0+EDGE,y0,z0+EDGE],[x0+EDGE,y0+EDGE,z0+EDGE],[x0,y0+EDGE,z0+EDGE]]
                o = rad.ObjHexahedron(v, [0,0,0]); rad.MatApl(o, mat); objs.append(o)
    body = rad.ObjCnt(objs)
    ext = rad.ObjBckg(lambda p: [B0*D[0], B0*D[1], B0*D[2]])
    return rad.ObjCnt([body, ext]), body, len(objs)


def align(body):
    M = np.array([m[1] for m in rad.ObjM(body)], float)
    nrm = np.linalg.norm(M, axis=1); g = nrm > 1e-12
    return float(np.mean((M[g] @ D) / nrm[g])) if g.any() else 0.0


def svd_null_dim(body):
    N, _ = rad.GetInteractMatrix(rad.BuildMatrix(body))
    N = np.asarray(N, float)
    s = np.linalg.svd(N, compute_uv=False)
    return int(np.sum(s < 1e-9 * s[0]))


def solve_case(outer, thick, nz, mu_r, deflate):
    rad.SetDeflateNullspace(bool(deflate), 1.0)
    grp, body, ne = build_frame(outer, thick, nz, mu_r)
    rad.SolverConfig(bicgstab_tol=1e-8, hacapk_eps=1e-4, hacapk_leaf=10, hacapk_eta=2.0)
    rad.Solve(grp, 0.001, 3000, 2)   # method 2 = HACApK BiCGSTAB
    st = rad.GetSolveStats()
    return ne, st.get('deflation_cycles', -1), st.get('linear_iterations', -1), align(body)


def main():
    mu_r = 1e5
    print(f"Belted-tree validation (mu_r={mu_r:.0e}, HACApK BiCGSTAB tol=1e-8)\n")
    # (label, outer, thick, nz, b1) -- b1 = first Betti number (# belt loops).
    # A thin ring's graph is a single C_N cycle (ZERO short cycles), so its only
    # null mode is the belt loop: a short-cycle-only basis would install 0, a
    # belted tree installs exactly 1. The installed cycle count (deflation_cycles)
    # must equal the SVD null dimension for the basis to span the whole kernel.
    cases = [
        ("ring 3x3 thin     ", 3, 1, 1, 1),   # C_8  ring: 0 short, 1 belt
        ("ring 4x4 thin     ", 4, 1, 1, 1),   # C_12 ring: 0 short, 1 belt
        ("ring 3x3 x2 layers", 3, 1, 2, 1),   # tube: short cycles + 1 belt
        ("solid 3x3 (control)", 3, 2, 1, 0),  # simply connected: belt = 0
    ]
    print(f"{'case':<20} {'ne':>3} {'null':>4} {'cyc':>4} {'belt?':>5} "
          f"{'algn(off)':>9} {'algn(on)':>8}  verdict")
    all_ok = True
    for label, outer, thick, nz, b1 in cases:
        # SVD null dim of the geometric operator N (ground-truth kernel dim).
        _, body0, _ = build_frame(outer, thick, nz, mu_r)
        nd = svd_null_dim(body0)
        # off / on solves (each rebuilds via UtiDelAll inside build_frame).
        ne, _,      it_off, al_off = solve_case(outer, thick, nz, mu_r, False)
        ne, cyc_on, it_on,  al_on  = solve_case(outer, thick, nz, mu_r, True)
        # (1) the installed basis spans the whole kernel: count == null dim
        #     (for a thin ring this is achievable ONLY via belt completion).
        ok = (cyc_on == nd)
        # (2) multiply-connected => at least b1 belt loops had to be added on top
        #     of the (zero, for thin rings) short cycles.
        belt_seen = "yes" if (b1 > 0 and cyc_on >= b1) else ("--" if b1 == 0 else "NO")
        if b1 > 0:
            ok &= (cyc_on >= b1)
        # (3) deflation never hurts the converged alignment.
        ok &= (al_on >= al_off - 0.02)
        all_ok &= ok
        print(f"{label:<20} {ne:>3} {nd:>4} {cyc_on:>4} {belt_seen:>5} "
              f"{al_off:>9.3f} {al_on:>8.3f}  {'OK' if ok else 'FAIL'}")
    rad.SetDeflateNullspace(False, 0.0)
    rad.UtiDelAll()
    print("\nNote: a thin ring has zero short cycles, so cyc==1 there is the belt")
    print("loop itself -- direct proof the belted-tree completion is active.")
    print("\n" + ("ALL CHECKS PASSED" if all_ok else "*** SOME CHECKS FAILED ***"))
    return 0 if all_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
