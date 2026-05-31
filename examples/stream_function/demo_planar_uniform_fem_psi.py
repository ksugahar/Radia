"""SF coil design with psi as direct H1 FE unknown on a surface mesh.

Same problem as demo_planar_uniform_coil.py (uniform Bz=B0 at z=h over a
square target region above a square source plane), but psi is now a
CONTINUOUS H1 GridFunction on a triangulated mesh of the source plane.

Why this matters:
  - psi lives in a true function space (H1) -- contour family deforms
    continuously under psi changes, no grid-spacing artefacts.
  - Essential BC psi=0 on the plane edge is built into the FES -> any
    contour automatically stays inside the plane (= manufacturable).
  - Optional Tikhonov H1-seminorm regularisation (--alpha A) controls
    smoothness vs accuracy.
  - Surface FES extends naturally to non-planar OCC surfaces (cylinder,
    sphere, conformal shielded geometry) by re-using the same code
    with a different mesh.
  - The Path-A compensated iteration should be smoother here than on
    the basis-loop demo because the level-set topology cannot jitter
    on a fixed mesh.

Pipeline:
  [1] 2D Netgen mesh of the source plane [-plane_half, plane_half]^2
  [2] H1(mesh, order=p, dirichlet="boundary") for psi
  [3] For each target j: assemble LinearForm L_j whose vector entries
      are A[j, i] = -mu0/(4*pi*r^3) * integral over mesh of
      grad(phi_i) . d_xy.  Stack into M x ndof matrix A.
  [4] Solve A . psi = B_target via SVD / Tikhonov-regularised lstsq.
  [5] Sample psi on a regular grid for marching-squares contour
      extraction; build single-stroke spiral (re-using the helpers
      from demo_planar_uniform_coil.py).
  [6][7] Field verification.

Reference for the SF-with-direct-FE-psi formulation: Liu, Hennig,
Korvink, "High-order-smooth discretised stream function coil design",
IEEE TM 48 (2012) 1179 -- they regularise via H2 seminorm; we offer
H1 here as a lighter knob.

Run: python demo_planar_uniform_fem_psi.py [--order 2] [--maxh 0.025]
                                            [--alpha 0.0] [--compensated-iter N]
"""
import argparse
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "src"))

from radia.biot_savart import h_segments_batch, MU0

# Re-use the planar contour extraction + spiral chain + field eval from
# the basis-loop demo; only psi solve differs.
from demo_planar_uniform_coil import (
    contour_polylines_xy, single_stroke_spiral_xy, chain_xy_to_3d,
    bz_at,
)


def build_fem_matrix(plane_half, maxh, order, dirichlet_bc, targets):
    """Mesh the plane and build M x ndof Biot-Savart matrix A.

    A[j, i] = Bz at targets[j] from psi = the i-th H1 basis function on
    the source plane (zero on the boundary).
    """
    import netgen.geom2d as g2d
    from ngsolve import (Mesh, H1, LinearForm, grad,
                          x as x_cf, y as y_cf, sqrt, dx, TaskManager)

    geo = g2d.SplineGeometry()
    geo.AddRectangle((-plane_half, -plane_half), (plane_half, plane_half),
                     bc=dirichlet_bc)
    ngmesh = geo.GenerateMesh(maxh=maxh)
    mesh = Mesh(ngmesh)
    fes = H1(mesh, order=order, dirichlet=dirichlet_bc)
    n_dof = fes.ndof
    n_free = int(sum(fes.FreeDofs()))
    print(f"  mesh: {mesh.ne} triangles, {mesh.nv} vertices, "
          f"H1 order={order}, ndof={n_dof}, free={n_free}")

    u, v = fes.TnT()
    M = len(targets)
    A = np.zeros((M, n_dof))
    mu0_4pi = 1.0e-7

    t0 = time.time()
    with TaskManager():
        for j, target in enumerate(targets):
            xt, yt, zt = float(target[0]), float(target[1]), float(target[2])
            dx_t = xt - x_cf
            dy_t = yt - y_cf
            r2 = dx_t * dx_t + dy_t * dy_t + zt * zt
            r3 = r2 * sqrt(r2)
            weight = -mu0_4pi / r3
            # Surface current K = z_hat x grad(psi) = (-d psi/dy, d psi/dx, 0)
            # Bz from K at target: (K x d_xy)_z = -grad(psi) . d_xy
            # Linear in psi -> weighted gradient integrated against test fn.
            L = LinearForm(fes)
            L += weight * (dx_t * grad(v)[0] + dy_t * grad(v)[1]) * dx
            L.Assemble()
            A[j, :] = L.vec.FV().NumPy()
    print(f"  matrix A: M={M} x ndof={n_dof}, assembled in {time.time()-t0:.2f} s")
    return A, fes, mesh, n_free


def build_h1_stiffness(fes, free_idx):
    """Assemble the H1 stiffness ``S = int grad u . grad v dA`` and
    return the dense free-DOF block ``S_free``."""
    from ngsolve import BilinearForm, grad, dx, TaskManager
    u, v = fes.TnT()
    a_stiff = BilinearForm(fes, symmetric=True)
    a_stiff += grad(u) * grad(v) * dx
    a_stiff += 1.0e-12 * u * v * dx       # tiny mass for SPD
    with TaskManager():
        a_stiff.Assemble()
    S = np.array(a_stiff.mat.ToDense())
    return S[np.ix_(free_idx, free_idx)]


def solve_tikhonov(A, B, fes, alpha=0.0, regularize="h1"):
    """Solve A . psi = B with one of three regularisation flavours.

    - regularize="l2" (alpha=0 default): min ||psi||_2 s.t. A psi = B
      via numpy lstsq.  Plain Euclidean min-norm; may have oscillations.
    - regularize="h1" (recommended): min ||grad psi||_S = psi^T S psi
      s.t. A psi = B exactly, where S is the H1 stiffness matrix.
      Lagrangian gives ``psi = S^-1 A^T (A S^-1 A^T)^-1 B`` -- the
      smoothest psi that hits the target.  No alpha to tune.
    - regularize="tikhonov" (alpha > 0): solve (A^T A + alpha S) psi
      = A^T B (the classical Tikhonov form).  Note: works only if
      alpha is below the smallest non-zero singular value squared of
      A_free (typically ~1e-13 for our scale); larger alpha makes the
      smoothness term dominate and the field fit collapses.
    """
    free = np.array(fes.FreeDofs())
    free_idx = np.where(free)[0]
    A_free = A[:, free_idx]
    if regularize == "l2" or (regularize == "tikhonov" and alpha <= 0):
        psi_free, _, rank, _ = np.linalg.lstsq(A_free, B, rcond=None)
    else:
        S_free = build_h1_stiffness(fes, free_idx)
        if regularize == "h1":
            # min psi^T S psi s.t. A_free psi = B
            S_inv_AT = np.linalg.solve(S_free, A_free.T)
            small_sys = A_free @ S_inv_AT
            lam = np.linalg.solve(small_sys, B)
            psi_free = S_inv_AT @ lam
        else:  # tikhonov
            lhs = A_free.T @ A_free + alpha * S_free
            rhs = A_free.T @ B
            psi_free = np.linalg.solve(lhs, rhs)
        rank = len(free_idx)
    psi = np.zeros(A.shape[1])
    psi[free_idx] = psi_free
    return psi, rank


# --------------------------------------------------------------------------
# Cached ACA+TSVD + regularised pseudo-inverse (Path-A inner solve)
# --------------------------------------------------------------------------
def build_regularized_cache(A, fes, regularize, alpha=0.0, aca_eps=1.0e-10,
                              verbose=False):
    """ACA+TSVD-factorise ``A_free`` and fold the regularisation stiffness.

    Returns the cached ``RegularizedTSVD`` plus the free-DOF index and
    the full-DOF length so subsequent solves can pad back to the full
    FES.  See ``radia.stream_function.RegularizedTSVD`` for the math.

    Supported ``regularize``:

      - ``"l2"``: S = identity -> reduces to standard min-Euclidean
        pseudo-inverse but ROUTED THROUGH the cache (same code path as
        h1, useful for Path-A timing comparisons).
      - ``"h1"``: S = int grad u . grad v dA (H1 seminorm + tiny mass).
      - ``"tikhonov"``: not folded into RegularizedTSVD (the AtA + alpha S
        normal-equations form has no clean (S^-1 V) factor) -- caller
        must drop back to ``solve_tikhonov`` for that mode.
    """
    from radia.stream_function import aca_tsvd, RegularizedTSVD
    free = np.array(fes.FreeDofs())
    free_idx = np.where(free)[0]
    A_free = A[:, free_idx]
    M_dim, N_dim = A_free.shape

    def entry(i, j):
        return float(A_free[i, j])

    t0 = time.time()
    res = aca_tsvd(M_dim, N_dim, entry, modes=M_dim,
                   kmax=min(M_dim, N_dim), aca_eps=aca_eps, method=3)
    t_aca = time.time() - t0

    if regularize == "l2":
        S_free = np.eye(N_dim)
    elif regularize == "h1":
        S_free = build_h1_stiffness(fes, free_idx)
    else:
        raise ValueError(
            f"regularize={regularize!r} not supported by RegularizedTSVD "
            f"cache; use solve_tikhonov directly")

    t0 = time.time()
    reg = RegularizedTSVD.from_stiffness(res, S_free)
    t_fold = time.time() - t0

    if verbose:
        print(f"     [cache] ACA+TSVD: rank={res.k_aca}/{res.modes}, "
              f"{t_aca*1000:.1f} ms")
        print(f"     [cache] S^-1 V + W_inv folded: {t_fold*1000:.1f} ms")

    return reg, free_idx, A.shape[1]


def solve_cached(reg, free_idx, n_full, B):
    """Pad the free-DOF cached solve back to full DOF (Dirichlet zero)."""
    psi_free = reg.solve(B)
    psi = np.zeros(n_full)
    psi[free_idx] = psi_free
    return psi


def sample_psi_grid(psi_vec, fes, mesh, plane_half, n_sample):
    """Sample a GridFunction psi on a regular grid for marching-squares."""
    from ngsolve import GridFunction
    gf = GridFunction(fes)
    gf.vec.FV().NumPy()[:] = psi_vec
    g = np.linspace(-plane_half, plane_half, n_sample)
    psi_grid = np.zeros((n_sample, n_sample))
    for k in range(n_sample):
        for l in range(n_sample):
            try:
                mp = mesh(g[k], g[l])
                psi_grid[k, l] = float(gf(mp))
            except Exception:
                psi_grid[k, l] = 0.0
    return psi_grid, g


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--plane-half", type=float, default=0.25)
    ap.add_argument("--maxh", type=float, default=0.025,
                    help="max triangle edge length on the source plane [m]")
    ap.add_argument("--order", type=int, default=2,
                    help="H1 polynomial order on the source plane")
    ap.add_argument("--regularize", choices=["l2", "h1", "tikhonov"],
                    default="h1",
                    help="psi regularisation: l2 = numpy min Euclidean norm "
                         "(may oscillate); h1 (default) = smoothest psi "
                         "exactly hitting B; tikhonov = mixed (alpha req)")
    ap.add_argument("--alpha", type=float, default=0.0,
                    help="Tikhonov mix weight (only used with "
                         "--regularize tikhonov)")
    ap.add_argument("--target-half", type=float, default=0.05)
    ap.add_argument("--target-z", type=float, default=0.10)
    ap.add_argument("--n-target", type=int, default=5)
    ap.add_argument("--B0", type=float, default=1.0e-3)
    ap.add_argument("--n-sample", type=int, default=81,
                    help="sample grid for psi (n_sample x n_sample)")
    ap.add_argument("--nlevels", type=int, default=12)
    ap.add_argument("--compensated-iter", type=int, default=0,
                    help="Path-A compensated iteration count")
    ap.add_argument("--compensated-step", type=float, default=1.0)
    args = ap.parse_args()

    print("=== SF -> planar uniform-Bz coil (H1 FE psi) ===")

    gt = np.linspace(-args.target_half, args.target_half, args.n_target)
    Xt, Yt = np.meshgrid(gt, gt, indexing="ij")
    targets = np.column_stack([Xt.ravel(), Yt.ravel(),
                                args.target_z * np.ones(Xt.size)])
    M = targets.shape[0]
    B_target = args.B0 * np.ones(M)

    print(f"[1] target: {M} obs at z={args.target_z*1e3:.0f}mm, "
          f"square half={args.target_half*1e3:.0f}mm, B0={args.B0*1e3:.2f} mT")

    A, fes, mesh, n_free = build_fem_matrix(
        args.plane_half, args.maxh, args.order, "boundary", targets)

    # Path-A iteration RE-USES the inner solve k times; for {l2, h1} we
    # build a RegularizedTSVD cache (ACA+TSVD factorisation + S^-1 V +
    # W_inv all precomputed) so each iter costs only O(k * (M + N)).
    use_cached = (args.compensated_iter > 0
                  and args.regularize in ("l2", "h1"))
    if use_cached:
        print(f"[2a] building RegularizedTSVD cache (regularize={args.regularize})")
        reg, free_idx, n_full = build_regularized_cache(
            A, fes, args.regularize, aca_eps=1.0e-10, verbose=True)
        psi_vec = solve_cached(reg, free_idx, n_full, B_target)
        rank = reg.base.modes
    else:
        psi_vec, rank = solve_tikhonov(A, B_target, fes, alpha=args.alpha,
                                        regularize=args.regularize)
        reg = None
        free_idx = None
        n_full = A.shape[1]
    A_psi = A @ psi_vec
    cont_res = float(np.linalg.norm(A_psi - B_target) /
                     (np.linalg.norm(B_target) + 1e-30))
    print(f"[2] continuous SF: ||A psi - B0||/||B0|| = {cont_res:.3e}"
          f" (alpha={args.alpha:.2g}, free-dof rank used={rank})")

    psi_grid, g_sample = sample_psi_grid(
        psi_vec, fes, mesh, args.plane_half, args.n_sample)

    polylines, dI = contour_polylines_xy(psi_grid, g_sample, g_sample,
                                          args.nlevels)
    n_wire = len(polylines)
    if n_wire == 0:
        print("[3] no contours; aborting")
        return
    print(f"[3] contours: {n_wire} closed curves, dI = {dI:.4g}")

    # ------- [3a] Path-A compensated iteration on the FE psi -----------
    if args.compensated_iter > 0:
        cache_tag = "cached" if use_cached else "direct"
        print(f"[3a] Path-A compensated iteration "
              f"({args.compensated_iter} iters, step={args.compensated_step},"
              f" inner={cache_tag})")
        best_psi = psi_vec.copy()
        best_polylines = polylines
        best_res_norm = float("inf")
        t_iter_sum = 0.0
        for it in range(args.compensated_iter):
            chain_it = single_stroke_spiral_xy(polylines)
            path_it = chain_xy_to_3d(chain_it)
            Bz_unit_it = bz_at(path_it, 1.0, targets)
            denom_it = float(np.dot(Bz_unit_it, Bz_unit_it))
            I_w_it = float(np.dot(Bz_unit_it, B_target) / denom_it) \
                if denom_it > 0 else 0.0
            residual = B_target - I_w_it * Bz_unit_it
            res_norm = float(np.linalg.norm(residual) /
                             (np.linalg.norm(B_target) + 1e-30))
            tag = ""
            if res_norm < best_res_norm:
                best_res_norm = res_norm
                best_psi = psi_vec.copy()
                best_polylines = polylines
                tag = " <-- best"
            print(f"     iter {it+1}: I_w = {I_w_it:.3e},"
                  f" residual = {res_norm:.3e}{tag}")
            t_inner0 = time.time()
            if use_cached:
                delta_psi = solve_cached(reg, free_idx, n_full, residual)
            else:
                delta_psi, _ = solve_tikhonov(A, residual, fes,
                                               alpha=args.alpha,
                                               regularize=args.regularize)
            t_iter_sum += time.time() - t_inner0
            psi_vec = psi_vec + args.compensated_step * delta_psi
            psi_grid, g_sample = sample_psi_grid(
                psi_vec, fes, mesh, args.plane_half, args.n_sample)
            polylines, dI = contour_polylines_xy(
                psi_grid, g_sample, g_sample, args.nlevels)
            if not polylines:
                print("     iter aborted: psi became flat")
                break
        psi_vec = best_psi
        polylines = best_polylines
        psi_grid, g_sample = sample_psi_grid(
            psi_vec, fes, mesh, args.plane_half, args.n_sample)
        n_done = max(1, it + 1)
        print(f"     final: best residual = {best_res_norm:.3e}, "
              f"inner-solve mean = {t_iter_sum / n_done * 1000:.2f} ms/iter")

    # ------- [4] single-stroke spiral ----------------------------------
    chain_xy = single_stroke_spiral_xy(polylines)
    path = chain_xy_to_3d(chain_xy)
    seg_lens = np.linalg.norm(np.diff(path, axis=0), axis=1)
    length = float(seg_lens.sum())
    print(f"[4] spiral chain: {len(path)} pts, {len(path)-1} segs,"
          f" length = {length:.3f} m")

    # ------- [6][7] field eval over target plane and on-axis ----------
    Bz_unit = bz_at(path, 1.0, targets)
    denom = float(np.dot(Bz_unit, Bz_unit))
    I_w = float(np.dot(Bz_unit, B_target) / denom) if denom > 0 else 0.0
    Bz_dsv = I_w * Bz_unit
    rms = float(np.linalg.norm(Bz_dsv - B_target) /
                (np.linalg.norm(B_target) + 1e-30))
    Bz_mean = float(np.mean(Bz_dsv))
    p2p = (float(np.max(Bz_dsv)) - float(np.min(Bz_dsv))) / (abs(Bz_mean) + 1e-30)
    print(f"[6] best-fit single current I_w = {I_w:.4g} A,"
          f" target-plane RMS = {rms:.3e}")
    print(f"    Bz over target: mean = {Bz_mean:.3e} T,"
          f" peak-to-peak / mean = {p2p:.3e}")
    print(f"    (target B0 = {args.B0:.3e} T)")

    zv = np.linspace(0.5 * args.target_z, 1.5 * args.target_z, 41)
    obs_axis = np.column_stack([np.zeros_like(zv), np.zeros_like(zv), zv])
    Bz_axis = bz_at(path, I_w, obs_axis)
    print(f"[7] on-axis Bz at z={args.target_z*1e3:.0f}mm: {Bz_axis[20]:.3e} T")

    # ------- plot -------------------------------------------------------
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig = plt.figure(figsize=(13, 4))
        ax1 = fig.add_subplot(1, 3, 1)
        cf = ax1.contourf(g_sample * 1e3, g_sample * 1e3, psi_grid.T,
                          levels=20, cmap="RdBu_r")
        lo, hi = psi_grid.min(), psi_grid.max()
        if hi > lo:
            lv = lo + (np.arange(args.nlevels) + 0.5) * (hi - lo) / args.nlevels
            ax1.contour(g_sample * 1e3, g_sample * 1e3, psi_grid.T,
                        levels=lv, colors="k", linewidths=0.5)
        ax1.set_title(f"H1 psi(x,y)  order={args.order}, free dofs={n_free}")
        ax1.set_xlabel("x [mm]"); ax1.set_ylabel("y [mm]")
        ax1.set_aspect("equal")
        fig.colorbar(cf, ax=ax1, fraction=0.045)

        ax2 = fig.add_subplot(1, 3, 2)
        ax2.plot(path[:, 0] * 1e3, path[:, 1] * 1e3, "b-", lw=0.5)
        ax2.set_title(f"Spiral chain ({len(path)} pts, {length:.2f} m)")
        ax2.set_xlabel("x [mm]"); ax2.set_ylabel("y [mm]")
        ax2.set_aspect("equal")
        ax2.grid(alpha=0.3)

        ax3 = fig.add_subplot(1, 3, 3)
        ax3.plot(zv * 1e3, Bz_axis * 1e3, "o-", ms=3, label="single-stroke coil")
        ax3.axhline(args.B0 * 1e3, color="k", ls="--",
                    label=f"target {args.B0*1e3:.2f} mT")
        ax3.axvline(args.target_z * 1e3, color="g", ls=":",
                    alpha=0.5, label=f"target z={args.target_z*1e3:.0f}mm")
        ax3.set_xlabel("z [mm]"); ax3.set_ylabel("Bz [mT]")
        ax3.set_title(f"On-axis Bz (RMS {rms:.1e}, p2p/mean {p2p:.1e})")
        ax3.legend(loc="best", fontsize=8); ax3.grid(alpha=0.3)

        out = HERE / "demo_planar_uniform_fem_psi.png"
        fig.tight_layout()
        fig.savefig(out, dpi=110)
        print(f"Saved plot: {out.name}")
    except ImportError:
        print("(matplotlib not installed; skipped plot)")


if __name__ == "__main__":
    main()
