#!/usr/bin/env python
"""Test the GENERAL C++ auto-deflation: rad.SetDeflateNullspace(True) makes the
HACApK solve build the local cycle basis from mesh topology (in C++) and apply
(A + alpha L L^T). Verify it cleans the solution on a STRUCTURED grid AND a
SHEARED (non-axis-aligned) grid -- the generalization goal.

Usage: python test_deflate_auto.py [nx ny nz] [mu_r] [shear] [alpha]
"""
import sys
import numpy as np
import radia as rad

EDGE = 0.01
D = np.array([1.0, 1.0, 1.0]) / np.sqrt(3.0)
B0 = 0.1


def build(nx, ny, nz, mu_r, shear=0.0):
    rad.UtiDelAll()
    mat = rad.MatLin(float(mu_r))
    objs = []
    for iz in range(nz):
        for iy in range(ny):
            for ix in range(nx):
                x0, y0, z0 = ix*EDGE, iy*EDGE, iz*EDGE
                base = [(x0,y0,z0),(x0+EDGE,y0,z0),(x0+EDGE,y0+EDGE,z0),(x0,y0+EDGE,z0),
                        (x0,y0,z0+EDGE),(x0+EDGE,y0,z0+EDGE),(x0+EDGE,y0+EDGE,z0+EDGE),(x0,y0+EDGE,z0+EDGE)]
                v = [[vx + shear*vz, vy, vz] for (vx, vy, vz) in base]
                o = rad.ObjHexahedron(v, [0,0,0]); rad.MatApl(o, mat); objs.append(o)
    cube = rad.ObjCnt(objs)
    ext = rad.ObjBckg(lambda p: [B0*D[0], B0*D[1], B0*D[2]])
    return rad.ObjCnt([cube, ext]), cube


def align(cube):
    M = np.array([m[1] for m in rad.ObjM(cube)], float)
    nrm = np.linalg.norm(M, axis=1); g = nrm > 1e-12
    return float(np.mean((M[g] @ D) / nrm[g]))


def main():
    nx, ny, nz = 4, 4, 2
    mu_r, shear, alpha = 1e5, 0.0, 1.0
    a = sys.argv[1:]
    if len(a) >= 3: nx, ny, nz = int(a[0]), int(a[1]), int(a[2])
    if len(a) >= 4: mu_r = float(a[3])
    if len(a) >= 5: shear = float(a[4])
    if len(a) >= 6: alpha = float(a[5])
    print("has SetDeflateNullspace:", hasattr(rad, "SetDeflateNullspace"))

    rad.SolverConfig(bicgstab_tol=1e-10, hacapk_eps=1e-4, hacapk_leaf=10, hacapk_eta=2.0)

    rad.SetDeflateNullspace(False, 0.0)
    grp, cube = build(nx, ny, nz, mu_r, shear); rad.Solve(grp, 0.001, 2000, 2)
    a_undef = align(cube)

    rad.SetDeflateNullspace(True, alpha)
    grp, cube = build(nx, ny, nz, mu_r, shear); rad.Solve(grp, 0.001, 2000, 2)
    a_def = align(cube)

    print(f"{nx}x{ny}x{nz} mu_r={mu_r:.0e} shear={shear}: "
          f"undeflated align={a_undef:.3f}  ->  auto-deflated align={a_def:.3f}  "
          f"{'OK' if a_def > a_undef + 0.2 else 'CHECK'}")
    rad.SetDeflateNullspace(False, 0.0)
    rad.UtiDelAll()


if __name__ == "__main__":
    main()
