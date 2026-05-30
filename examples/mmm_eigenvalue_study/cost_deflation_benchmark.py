#!/usr/bin/env python
"""Cost study: loop-basis deflation (shifted BiCGSTAB via L) vs plain HACApK.

For a range of cube sizes (mu_r high so the loop modes are active), compare:
  (P) plain HACApK BiCGSTAB  : SetDeflateNullspace(False)
  (D) deflated (loop shift)  : SetDeflateNullspace(True, alpha)
measuring per solve: linear iterations, wall time, and the solution's mean
alignment (quality: ~0.1 = ugly/loop-dominated, ~0.9 = clean).

This answers "is the loop-basis machinery worth it?": does deflation reduce
iterations (offsetting the L-build + L L^T overhead), and does it actually
clean the solution?
"""
import sys, time
import numpy as np
import radia as rad

EDGE = 0.01
D = np.array([1.0, 1.0, 1.0]) / np.sqrt(3.0)
B0 = 0.1


def build(n, mu_r):
    rad.UtiDelAll()
    mat = rad.MatLin(float(mu_r))
    objs = []
    for iz in range(n):
        for iy in range(n):
            for ix in range(n):
                x0, y0, z0 = ix*EDGE, iy*EDGE, iz*EDGE
                v = [[x0,y0,z0],[x0+EDGE,y0,z0],[x0+EDGE,y0+EDGE,z0],[x0,y0+EDGE,z0],
                     [x0,y0,z0+EDGE],[x0+EDGE,y0,z0+EDGE],[x0+EDGE,y0+EDGE,z0+EDGE],[x0,y0+EDGE,z0+EDGE]]
                o = rad.ObjHexahedron(v, [0,0,0]); rad.MatApl(o, mat); objs.append(o)
    cube = rad.ObjCnt(objs)
    ext = rad.ObjBckg(lambda p: [B0*D[0], B0*D[1], B0*D[2]])
    return rad.ObjCnt([cube, ext]), cube


def align(cube):
    M = np.array([m[1] for m in rad.ObjM(cube)], float)
    nrm = np.linalg.norm(M, axis=1); g = nrm > 1e-12
    return float(np.mean((M[g] @ D) / nrm[g]))


def solve_one(n, mu_r, deflate, alpha):
    rad.SetDeflateNullspace(bool(deflate), alpha)
    grp, cube = build(n, mu_r)
    rad.SolverConfig(bicgstab_tol=1e-8, hacapk_eps=1e-4, hacapk_leaf=10, hacapk_eta=2.0)
    t0 = time.perf_counter()
    rad.Solve(grp, 0.001, 3000, 2)
    dt = time.perf_counter() - t0
    st = rad.GetSolveStats()
    return dt, st.get('linear_iterations', -1), align(cube)


def main():
    mu_r = 1e5
    sizes = [4, 6, 8, 10, 12]
    a = sys.argv[1:]
    if a: sizes = [int(x) for x in a]
    print(f"mu_r={mu_r:.0e}, BiCGSTAB tol=1e-8")
    print(f"{'n':>3} {'dof':>7} | {'P_iter':>6} {'P_t(s)':>7} {'P_algn':>7} | "
          f"{'D_iter':>6} {'D_t(s)':>7} {'D_algn':>7} | {'t_ratio':>7}")
    for n in sizes:
        dof = n*n*n*6
        pt, pit, pal = solve_one(n, mu_r, False, 0.0)
        dt, dit, dal = solve_one(n, mu_r, True, 1.0)
        ratio = dt/pt if pt > 0 else float('nan')
        print(f"{n:>3} {dof:>7} | {str(pit):>6} {pt:>7.2f} {pal:>7.3f} | "
              f"{str(dit):>6} {dt:>7.2f} {dal:>7.3f} | {ratio:>7.2f}")
    rad.SetDeflateNullspace(False, 0.0)
    rad.UtiDelAll()


if __name__ == "__main__":
    main()
