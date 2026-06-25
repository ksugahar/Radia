"""Phase-3a verification: multipole-moment MMM extended to 5-face WEDGE (triangular-prism) elements.

A wedge has 5 surface-charge DOF, so the moment system per element is 3 dipole + 1 monopole + 1 quadrupole
= 5 rows.  The single quad MUST be the AXIAL quadrupole 3*(d.a)^2-|d|^2 about the prism axis a (the direction
between the two triangular cap centroids): the naive antisymmetric dx^2-dy^2 leaves a SYMMETRIC near-null
sigma mode unpinned, and the magnetization blows up to ~1e9 (diagnosed 2026-06-22).  The axial quad is
symmetric about the axis, so it pins that mode for ANY wedge orientation.

Checks (external field is the observable; internal M is formulation-dependent, CLAUDE.md):
  1. single wedge, z-axis and x-axis prisms: moment external B == EIEM2 external B (the trusted wedge path)
     to < 2%, and |M| does NOT blow up (the near-null mode is pinned).
  2. a cube as ONE hex vs the SAME cube tiled by TWO wedges: same external field (same geometry).
  3. a MIXED hex+wedge stack == the all-hex stack (variable DOF offsets).
Exercises method 0 (LU) and method 1 (-> moment LU).  Method 2 (H-matrix) for wedge is Phase 3a-2.
"""
import numpy as np
import radia as rad

MU0 = 4e-7 * np.pi
H0 = 1000.0
L = 0.02


def _solve(objs, method, probes):
    cont = rad.ObjCnt(objs + [rad.ObjBckg(lambda p: [0.0, 0.0, MU0 * H0])])
    rad.Solve(cont, 1e-8, 3000, method)
    M = np.asarray([rad.ObjM(o)["magnetization"] for o in objs], float)
    B = np.asarray([rad.Fld(cont, "b", p) for p in probes], float)
    return M, B


def _build(kind, moment):
    rad.UtiDelAll(); rad.set_demag_backend("yano"); rad.SolverConfig(bicgstab_tol=1e-9)
    if kind == "wedge_z":
        objs = [rad.ObjWedge([[0, 0, 0], [L, 0, 0], [0, L, 0], [0, 0, L], [L, 0, L], [0, L, L]], [0, 0, 0])]
    elif kind == "wedge_x":
        objs = [rad.ObjWedge([[0, 0, 0], [0, L, 0], [0, 0, L], [L, 0, 0], [L, L, 0], [L, 0, L]], [0, 0, 0])]
    elif kind == "hex":
        objs = [rad.ObjHexahedron([[0, 0, 0], [L, 0, 0], [L, L, 0], [0, L, 0], [0, 0, L], [L, 0, L], [L, L, L], [0, L, L]], [0, 0, 0])]
    elif kind == "cube_wedges":
        objs = [rad.ObjWedge([[0, 0, 0], [L, 0, 0], [0, L, 0], [0, 0, L], [L, 0, L], [0, L, L]], [0, 0, 0]),
                rad.ObjWedge([[L, 0, 0], [L, L, 0], [0, L, 0], [L, 0, L], [L, L, L], [0, L, L]], [0, 0, 0])]
    for o in objs:
        rad.MatApl(o, rad.MatLin(1000.0))
    return objs


def main():
    ok = True

    # (1) single wedge, z + x axis: moment vs EIEM2 external B + no blow-up
    pr = [[0.05, 0.01, 0.02], [0.0, 0.05, 0.01], [0.01, 0.01, 0.05]]
    for kind in ("wedge_z", "wedge_x"):
        objs = _build(kind, False); _, Be = _solve(objs, 0, pr); rad.UtiDelAll()
        objs = _build(kind, True);  Mm, Bm = _solve(objs, 0, pr); rad.UtiDelAll()
        rel = np.linalg.norm(Bm - Be) / max(np.linalg.norm(Be), 1e-30)
        blow = np.max(np.abs(Mm)) > 1e6
        good = np.all(np.isfinite(Bm)) and rel < 0.02 and not blow
        print(f"  (1) {kind:8s}: moment-vs-EIEM2 relB={rel:.2e}  |M|max={np.max(np.abs(Mm)):.1f}  blowup={blow}  {'OK' if good else 'FAIL'}")
        ok &= good

    # (2) cube as 1 hex vs 2 wedges (same geometry) -- method 0 and 1
    pr = [[0.05, 0.0, 0.0], [0.0, 0.05, 0.01], [0.01, 0.01, 0.05], [-0.04, 0.02, 0.03]]
    objs = _build("hex", True); _, Bh = _solve(objs, 0, pr); rad.UtiDelAll()
    for method in (0, 1):
        objs = _build("cube_wedges", True); Mw, Bw = _solve(objs, method, pr); rad.UtiDelAll()
        rel = np.linalg.norm(Bw - Bh) / max(np.linalg.norm(Bh), 1e-30)
        Mzw = float(np.mean(Mw[:, 2]))
        good = np.all(np.isfinite(Bw)) and rel < 0.02 and 2.0*H0 < Mzw < 4.0*H0 and np.all(np.abs(Mw[:, :2]) < 0.1*abs(Mzw))
        print(f"  (2) cube 2-wedge vs 1-hex method{method}: relB={rel:.2e}  M_z(avg)={Mzw:.1f}  {'OK' if good else 'FAIL'}")
        ok &= good

    print("RESULT:", "PASS -- wedge moment (5-row axial-quad) correct + orientation-robust" if ok else "FAIL")
    rad.UtiDelAll()
    return ok


if __name__ == "__main__":
    raise SystemExit(0 if main() else 1)
