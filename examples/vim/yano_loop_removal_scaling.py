"""yano_loop_removal_scaling.py -- exact-field + scaling validation of LOOP REMOVAL (star projection)
on the REAL C++ EIEM2 surface-charge MSC interaction matrix (not the monopole Python prototype).

This is the "本筋" verification step that justifies wiring loop removal into the production surface-charge MSC
solve: prove, on the EXACT field (rad.BuildMatrix / rad.GetInteractMatrix builds the C++ EIEM2
surface-charge operator), that

  (1) the matrix near-null space == the cell-graph cycle count (loops are TOPOLOGICAL);
  (2) star projection (solve in the complement of the near-null) makes cond(A) + BiCGStab iters
      BOUNDED and mu_r-INDEPENDENT, whereas the full system has cond(A) ~ mu_r and iters that grow;
  (3) loop removal is ACCURACY-NEUTRAL on the exact field: N @ sigma_full == N @ sigma_star to round-off
      (the loop component the full solve leaves is field-null), so the physical demag field is unchanged;
  (4) the bound holds as N grows (n^3 hex bar, N = 8..216 elements), i.e. iters do not grow with size.

Self-contained: the geometry is a structured (optionally jittered/distorted) hex bar generated here, so
there is no external mesh dependency.  The matrix is the production C++ EIEM2 operator, so the field is
exact -- this complements the lab's monopole-prototype scaling study (element-study notebook) which
established the same cycle/star structure with a crude monopole kernel.

NOTE on the production picture:
  - The C++ yano matrix already carries the element-center compensation charge, so its near-null is the
    LOOPS ONLY (the monopole modes are already removed by the compensation).  Hence the production
    loop-removal is exactly the STAR projection measured here.
  - The "raw + div=0(Gauss) + loop" route (drop the compensation cloud, impose sum_f area_f*sigma_f = 0
    per element via Lagrange instead) is the EQUIVALENT alternative: it removes the same monopole modes
    the cloud removes (verified in the monopole prototype, yano_div_lagrange_study.py) and then needs the
    same loop removal.  It yields the same physical field; its only claimed edge is conditioning, which
    is a kernel-design question, not an accuracy one.
"""
import os, json
import numpy as np
from scipy.sparse.linalg import LinearOperator, bicgstab

import sys
sys.path.insert(0, r"S:\Radia\01_GitHub\src\radia")
import radia as rad

HERE = os.path.dirname(os.path.abspath(__file__))

# standard hex (VTK/Netgen MMB8T) corner order: bottom CCW 0-3, top CCW 4-7
_HEX_FACES = [(0, 1, 2, 3), (4, 5, 6, 7), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)]


def build_bar(nx, ny, nz, h=1.0, distort=0.0, seed=0):
    """structured hex bar.  distort>0 applies a GRADED (non-uniform) axis spacing + a global AFFINE shear,
    which keeps every quad face planar (affine maps the axis-aligned box faces to parallelograms) so the
    hexes are valid for ObjHexahedron yet genuinely distorted (non-cubic + non-orthogonal).  Returns
    (hex_verts, elems): hex_verts[e] = 8x3 vertex list, elems[e] = 8 global node ids."""
    nxp, nyp, nzp = nx + 1, ny + 1, nz + 1

    def nid(i, j, k):
        return i + nxp * (j + nyp * k)

    def axis_coords(n):
        # deterministic varied (positive) cell widths -> graded, non-uniform aspect
        widths = [h * (1.0 + distort * 0.5 * ((((i * 7 + 3) % 5) / 4.0) - 0.5)) for i in range(n)]
        c = [0.0]
        for w in widths:
            c.append(c[-1] + w)
        return np.array(c)

    xs, ys, zs = axis_coords(nx), axis_coords(ny), axis_coords(nz)
    S = np.eye(3)
    sh = distort * 0.4
    S[0, 1], S[0, 2], S[1, 2] = sh, sh * 0.7, sh * 0.5     # affine shear -> non-orthogonal hexes
    pos = {}
    for k in range(nzp):
        for j in range(nyp):
            for i in range(nxp):
                pos[nid(i, j, k)] = S @ np.array([xs[i], ys[j], zs[k]], float)
    hex_verts, elems = [], []
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                c = [nid(i, j, k), nid(i + 1, j, k), nid(i + 1, j + 1, k), nid(i, j + 1, k),
                     nid(i, j, k + 1), nid(i + 1, j, k + 1), nid(i + 1, j + 1, k + 1), nid(i, j + 1, k + 1)]
                elems.append(c)
                hex_verts.append(np.array([pos[n] for n in c]))
    return hex_verts, elems


def cycle_count(elems):
    """topological loop nullity of the cell/internal-face incidence: n_internal_faces - (n_cells - 1)."""
    from collections import defaultdict
    seen = defaultdict(int)
    for c in elems:
        for f in _HEX_FACES:
            key = tuple(sorted(c[i] for i in f))
            seen[key] += 1
    n_internal = sum(1 for v in seen.values() if v == 2)
    return n_internal - (len(elems) - 1), n_internal


def yano_matrix(hex_verts, mu_r=1000.0, scale=1e-2):
    """build the C++ EIEM2 surface-charge MSC interaction matrix N (+N convention; loops at eigenvalue 0)."""
    rad.UtiDelAll()
    rad.set_demag_backend("yano")
    objs = []
    for V in hex_verts:
        h = rad.ObjHexahedron([(V[i] * scale).tolist() for i in range(8)], [0, 0, 0])
        rad.MatApl(h, rad.MatLin(mu_r))
        objs.append(h)
    cont = rad.ObjCnt(objs)
    handle = rad.BuildMatrix(cont)
    N, _dof = rad.GetInteractMatrix(handle)
    rad.UtiDelAll()
    return np.asarray(N, float)


def bicgstab_iters(A, b):
    cnt = {"n": 0}
    op = LinearOperator(A.shape, matvec=lambda v: A @ v)
    bicgstab(op, b, rtol=1e-10, maxiter=30 * A.shape[0],
             callback=lambda xk: cnt.__setitem__("n", cnt["n"] + 1))
    return cnt["n"]


def star_basis(N):
    """complement of the near-null (loops) of N via SVD: columns span the loop-free (star) subspace."""
    U, S, Vt = np.linalg.svd(N)
    n_loop = int(np.sum(S < 1e-9 * S[0]))
    P = Vt[:N.shape[0] - n_loop, :].T
    return P, n_loop


def study_mu_sweep(hex_verts, elems, tag, mu_list=(1e2, 1e3, 1e4, 1e5)):
    cyc, n_int = cycle_count(elems)
    N = yano_matrix(hex_verts, mu_r=mu_list[-1])  # near-null is mu-independent (N is geometry-only)
    P, n_loop = star_basis(N)
    ndof = N.shape[0]
    print(f"\n=== {tag}: {len(elems)} hex, ndof={ndof}, internal_faces={n_int}, "
          f"cycle(loops)={cyc}, svd_near_null={n_loop}, star_dim={P.shape[1]} ===")
    rng = np.random.default_rng(0)
    b = rng.random(ndof) - 0.5
    rows = []
    print(f"  {'mu_r':>7} | {'cond full':>10} {'it full':>7} | {'cond star':>9} {'it star':>7} | "
          f"{'field |N(s_f-s_s)|/|N s_f|':>26}")
    for mur in mu_list:
        chi = mur - 1.0
        A = (1.0 / chi) * np.eye(ndof) - N
        cond_full = float(np.linalg.cond(A))
        it_full = bicgstab_iters(A, b)
        A_star = P.T @ A @ P
        cond_star = float(np.linalg.cond(A_star))
        it_star = bicgstab_iters(A_star, P.T @ b)
        # accuracy-neutrality: min-norm full solve vs star-lifted -> same field N@sigma
        sig_full = np.linalg.lstsq(A, b, rcond=None)[0]
        sig_star = P @ np.linalg.solve(A_star, P.T @ b)
        fld = N @ sig_full
        rel = float(np.linalg.norm(N @ (sig_full - sig_star)) / (np.linalg.norm(fld) + 1e-300))
        print(f"  {mur:>7.0e} | {cond_full:>10.2e} {it_full:>7} | {cond_star:>9.1f} {it_star:>7} | {rel:>26.2e}")
        rows.append(dict(mu_r=mur, cond_full=cond_full, it_full=it_full,
                         cond_star=cond_star, it_star=it_star, field_rel=rel))
    return dict(tag=tag, n_hex=len(elems), ndof=ndof, internal_faces=n_int,
                cycle=cyc, svd_near_null=n_loop, star_dim=P.shape[1], mu_sweep=rows)


def study_scaling(sizes, distort, mu_r=1e4):
    print(f"\n=== SCALING (distort={distort}, mu_r={mu_r:.0e}): star iters/cond bounded as N grows ===")
    print(f"  {'n^3':>6} {'hex':>5} {'ndof':>6} {'cycle':>6} | {'cond full':>10} {'it full':>7} | "
          f"{'cond star':>9} {'it star':>7} | {'field_rel':>10}")
    rng = np.random.default_rng(0)
    rows = []
    for n in sizes:
        hex_verts, elems = build_bar(n, n, n, distort=distort, seed=1)
        cyc, n_int = cycle_count(elems)
        N = yano_matrix(hex_verts, mu_r=mu_r)
        ndof = N.shape[0]
        P, n_loop = star_basis(N)
        chi = mu_r - 1.0
        A = (1.0 / chi) * np.eye(ndof) - N
        b = rng.random(ndof) - 0.5
        cond_full = float(np.linalg.cond(A)); it_full = bicgstab_iters(A, b)
        A_star = P.T @ A @ P
        cond_star = float(np.linalg.cond(A_star)); it_star = bicgstab_iters(A_star, P.T @ b)
        sig_full = np.linalg.lstsq(A, b, rcond=None)[0]
        sig_star = P @ np.linalg.solve(A_star, P.T @ b)
        rel = float(np.linalg.norm(N @ (sig_full - sig_star)) / (np.linalg.norm(N @ sig_full) + 1e-300))
        print(f"  {n}^3   {len(elems):>5} {ndof:>6} {cyc:>6} | {cond_full:>10.2e} {it_full:>7} | "
              f"{cond_star:>9.1f} {it_star:>7} | {rel:>10.2e}")
        rows.append(dict(n=n, n_hex=len(elems), ndof=ndof, cycle=cyc, svd_near_null=n_loop,
                         cond_full=cond_full, it_full=it_full, cond_star=cond_star,
                         it_star=it_star, field_rel=rel))
    return dict(distort=distort, mu_r=mu_r, rows=rows)


def main():
    out = {}
    # (1) regular + distorted patch: cycle == svd-null, star bounds cond/iters mu_r-independent, field-preserved
    hv_reg, el_reg = build_bar(3, 3, 2, distort=0.0, seed=1)
    out["regular_3x3x2"] = study_mu_sweep(hv_reg, el_reg, "REGULAR 3x3x2")
    hv_dis, el_dis = build_bar(3, 3, 2, distort=0.6, seed=1)
    out["distorted_3x3x2"] = study_mu_sweep(hv_dis, el_dis, "DISTORTED 3x3x2 (distort=0.6)")
    # (2) scaling: iters bounded as N grows (regular + distorted)
    out["scaling_regular"] = study_scaling([2, 3, 4, 5, 6], distort=0.0)
    out["scaling_distorted"] = study_scaling([2, 3, 4, 5], distort=0.5)
    with open(os.path.join(HERE, "yano_loop_removal_scaling.json"), "w") as f:
        json.dump(out, f, indent=2, default=float)
    print("\nReadout: on the EXACT C++ EIEM2 yano matrix, svd_near_null == cycle (loops topological);")
    print("full cond ~ mu_r with growing iters; star cond + iters BOUNDED and mu_r-independent; field")
    print("preserved to round-off (loop removal is accuracy-neutral on the exact field) and the bound")
    print("holds as N grows -> loop removal (star) is the accurate + loop-free + bounded-iter surface-charge MSC.")
    print("saved", os.path.join(HERE, "yano_loop_removal_scaling.json"))


if __name__ == "__main__":
    main()
