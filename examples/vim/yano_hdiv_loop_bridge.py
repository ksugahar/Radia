"""yano-MSC <-> HDiv-VIM loop duality: the SAME loop space, two backends.

Radia keeps BOTH soft-iron demag backends (decision 2026-06-19): the collocation surface-charge
yano-type MSC (6 sigma DOF / hex, the scalable HACApK demag) and the FEEC HDiv-VIM (RT0 face flux,
loop-free by construction).  This example shows -- on one regular and one affine-sheared hex grid --
that they share the SAME loop structure, and that a star projection re-conditions the collocation
matrix the way the HDiv-VIM is conditioned by construction.

The structural fact (verified below, all cases exact):

  cell-adjacency cycle space  ==  HDiv-VIM ker(B)  ==  yano-MSC collocation near-null
  cycle = n_internal_faces - (n_cells - 1)   (first Betti number of the cell graph)

  (1) HDiv-VIM  N = B^T G B : the loops are ker(B), so N.loop = 0 EXACTLY (loop_res ~ 1e-16).
      This is why the HDiv-VIM nonlinear (Newton) iteration is mesh/mu_r-independent: the loop
      modes never enter the operator.
  (2) yano-MSC  collocation : the 6-sigma/hex matrix carries those same loops as a LATENT near-null
      subspace -> condition number ~1e16-1e17.  The loop SOURCE field cancels (element-common
      compensation), but the loop COEFFICIENTS remain unknowns, so an iterative solve can grow them.
  (3) star projection : restrict the unknowns to the non-loop (internal-flux) subspace -- the dense
      research form of the sparse q_internal = B^T a.  Drops the yano-MSC condition number from
      ~1e16 to ~40-65, i.e. into the HDiv-VIM regime, with an O(n_face) sparse projector (NOT a heavy
      Galerkin re-assembly).

Note: the HDiv mesh (MakeStructured3DMesh) is regular even for the sheared yano cells; the loop COUNT
is a topological invariant (shear-independent), so n_loop == cycle holds regardless -- the shear only
exercises the yano collocation kernel on non-axis-aligned hexes.

Reference: H. Yano and K. Sugahara, "Magnetic Moment Method with the Idea of Magnetic Surface Charge
Method", J. Magn. Soc. Jpn. 47(2), 52-56, 2023 (the yano-type MSC).  The loop / de Rham structure is
the FEEC RT0 reading of the same surface-charge discretization.
"""
import os
import sys
import json

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../../src/radia"))
import radia as rad
import ngsolve as ng
from ngsolve.meshes import MakeStructured3DMesh
from radia.vim import _core

HERE = os.path.dirname(os.path.abspath(__file__))


def hex_grid_cells(nx, ny, nz, shear=0.0):
    """Regular [0,1]^3 grid of nx*ny*nz hexes, optionally globally AFFINE-sheared.  A linear map
    preserves planarity, so the hex faces stay planar (ObjHexahedron requires planar faces; a
    per-node jitter does NOT -- it warps faces out of plane -- so we shear the whole grid instead).
    Returns (cells: list of 8x3 vertex arrays in VTK hex order, shared-internal-face count)."""
    nxp, nyp, nzp = nx + 1, ny + 1, nz + 1
    S = np.array([[1.0, 0.4 * shear, shear], [0.0, 1.0, 0.6 * shear], [0.0, 0.0, 1.0]])
    P = np.zeros((nxp, nyp, nzp, 3))
    for i in range(nxp):
        for j in range(nyp):
            for k in range(nzp):
                P[i, j, k] = S @ np.array([i / nx, j / ny, k / nz])
    cells = []
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                v = [P[i, j, k], P[i + 1, j, k], P[i + 1, j + 1, k], P[i, j + 1, k],
                     P[i, j, k + 1], P[i + 1, j, k + 1], P[i + 1, j + 1, k + 1], P[i, j + 1, k + 1]]
                cells.append(np.array(v))
    n_internal = (nx - 1) * ny * nz + nx * (ny - 1) * nz + nx * ny * (nz - 1)
    return cells, n_internal


def yano_matrix(cells, mu_r=1000.0):
    """Radia yano-MSC collocation interaction matrix (6 surface-charge DOF / hex) via the public
    BuildMatrix / GetInteractMatrix accessors."""
    rad.UtiDelAll()
    objs = []
    for V in cells:
        h = rad.ObjHexahedron(V.tolist(), [0.0, 0.0, 0.0])
        rad.MatApl(h, rad.MatLin(mu_r))
        objs.append(h)
    cont = rad.ObjCnt(objs)
    handle = rad.BuildMatrix(cont)
    A, dof = rad.GetInteractMatrix(handle)
    return np.asarray(A).reshape(dof, dof), dof


def report(nx, ny, nz, shear=0.0):
    cells, n_internal = hex_grid_cells(nx, ny, nz, shear)
    n_cells = len(cells)
    cycle = n_internal - (n_cells - 1)        # first Betti of the cell-adjacency graph (the loops)
    tag = f"{nx}x{ny}x{nz} shear={shear}"

    # (1) HDiv-VIM: the loop-free RT0 internal-flux operator.  N = B^T G B is applied via the production
    # C++ polytope charge-Gram H-matrix (the dense Python Gram path was removed); on these small structured
    # grids we densify N from the C++ H-matvec and take the loops = ker(B) of the sparse charge map.
    import radia._radia_pybind as _rp
    with ng.TaskManager():
        mesh = MakeStructured3DMesh(hexes=True, nx=nx, ny=ny, nz=nz)
        d = _core.build_demag(mesh)
        B = d["B_csr"]; ndof = int(d["ndof"])
        p = d["poly"]
        Hg = _rp._ChargeGramHMatrix(
            cell_tris=list(p["cell_tris"]), cell_troff=list(p["cell_troff"]),
            cell_cent=list(p["cell_cent"]), cell_meas=list(p["cell_meas"]),
            face_tris=list(p["face_tris"]), face_troff=list(p["face_troff"]),
            face_cent=list(p["face_cent"]), face_meas=list(p["face_meas"]),
            n_el=int(d["n_el"]), eps=1e-7, leaf=64, eta=2.0)
        N = np.column_stack([B.T @ np.asarray(Hg.matvec((B @ np.eye(ndof)[:, k]).tolist()), float)
                             for k in range(ndof)])
        Bd = B.toarray()
        sv_B = np.linalg.svd(Bd, compute_uv=False)
        rankB = int(np.sum(sv_B > 1e-9 * sv_B.max()))
        loops = np.linalg.svd(Bd)[2][rankB:, :]
    n_loop = ndof - rankB
    Nn = np.linalg.norm(N) or 1.0
    loop_res = (max(np.linalg.norm(N @ loops[k]) for k in range(n_loop)) / Nn) if n_loop else 0.0

    # (2) yano-MSC: collocation matrix + latent near-null (the loops) + conditioning
    A, dof = yano_matrix(cells)
    U, sv, Vt = np.linalg.svd(A)
    smax, smin = sv[0], sv[-1]
    cond = smax / smin if smin > 0 else float("inf")
    near_null = int(np.sum(sv < 1e-9 * smax))

    # (3) STAR PROJECTION (dense-SVD research form of q = B^T a): drop the `cycle` loop directions
    n_keep = dof - cycle
    P = Vt[:n_keep].T
    sv_star = np.linalg.svd(P.T @ A @ P, compute_uv=False)
    cond_star = sv_star[0] / sv_star[-1] if sv_star[-1] > 0 else float("inf")

    n_boundary = 6 * n_cells - 2 * n_internal
    print(f"\n=== {tag} : n_cells={n_cells}  n_internal_faces={n_internal}  cycle(loops)={cycle} ===")
    print(f"  [HDiv-VIM ] ndof(RT0)={d['ndof']:4d}  n_loop={n_loop:4d}  "
          f"loop_res(max||N loop||/||N||)={loop_res:.2e}")
    print(f"  [yano-MSC ] dof(6/hex)={dof:4d}  cond2={cond:.3e}  near_null(<1e-9 smax)={near_null}")
    print(f"  [star-proj] keep={n_keep} (drop {cycle} loops)  cond2_star={cond_star:.3e}  "
          f"-> {cond / cond_star:.1e}x better")
    print(f"  [bridge   ] local_face_dof=6*ncell={6 * n_cells}=2*int+bnd ({2 * n_internal}+{n_boundary})")
    ok = (n_loop == cycle == near_null)
    print(f"  [match    ] HDiv n_loop {n_loop} == cycle {cycle} == yano near_null {near_null}: "
          f"{'YES' if ok else 'NO'};  star fixes cond: {'YES' if cond_star < 1e6 else 'NO'}")
    return dict(tag=tag, n_cells=n_cells, n_internal=n_internal, cycle=cycle,
                hdiv_ndof=d["ndof"], hdiv_n_loop=n_loop, hdiv_loop_res=loop_res,
                yano_dof=dof, yano_cond=cond, yano_near_null=near_null,
                yano_cond_star=cond_star, n_keep=n_keep, match=ok)


def mur_sweep(nx, ny, nz, mur_list=(10, 100, 1000, 1e4, 1e5)):
    """The latent loops are REGULARIZED by the 1/chi diagonal of the system matrix A = -N + (1/chi)I
    (N is geometry-only; chi = mu_r - 1).  So they bite only at HIGH mu_r (1/chi -> 0): cond(A) grows
    ~mu_r and the iterative solve slows -- exactly the yano-MSC mesh/mu_r-dependence (the HDiv-VIM is
    mu_r-INDEPENDENT because the loops are in ker(B)).  Star-projecting the loops out of A makes the
    yano-MSC condition number mu_r-INDEPENDENT too."""
    cells, n_internal = hex_grid_cells(nx, ny, nz, shear=0.5)
    n_cells = len(cells)
    cycle = n_internal - (n_cells - 1)
    N, dof = yano_matrix(cells)               # raw geometry-only interaction (+N)
    _, _, Vt = np.linalg.svd(N)
    P = Vt[: dof - cycle].T                    # non-loop subspace (mu_r-independent: geometry only)
    print(f"\n=== mu_r sweep on {nx}x{ny}x{nz} shear=0.5 (cycle={cycle} loops) ===")
    print(f"  {'mu_r':>8} | {'cond(A=-N+I/chi)':>18} | {'cond(star A)':>14}")
    rows = []
    for mur in mur_list:
        chi = mur - 1.0
        A = -N + (1.0 / chi) * np.eye(dof)
        cA = np.linalg.cond(A)
        cS = np.linalg.cond(P.T @ A @ P)
        print(f"  {mur:>8.0f} | {cA:>18.3e} | {cS:>14.3e}")
        rows.append(dict(mu_r=float(mur), cond_full=float(cA), cond_star=float(cS)))
    return dict(grid=f"{nx}x{ny}x{nz}", cycle=cycle, sweep=rows)


def main():
    out = [report(2, 2, 2, shear=0.0),
           report(3, 3, 2, shear=0.0),
           report(3, 3, 2, shear=0.5),
           report(3, 3, 3, shear=0.5)]
    sweep = mur_sweep(3, 3, 2)
    all_ok = all(r["match"] for r in out) and all(r["yano_cond_star"] < 1e6 for r in out)
    with open(os.path.join(HERE, "yano_hdiv_loop_bridge.json"), "w") as f:
        json.dump({"all_ok": all_ok, "cases": out, "mur_sweep": sweep}, f, indent=2, default=float)
    print(f"\nall cases match + star-conditioned: {all_ok}")
    rad.UtiDelAll()
    return all_ok


if __name__ == "__main__":
    main()
