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

INV_4PI = 1.0 / (4.0 * np.pi)


def bz_fast(path, current, obs, chunk=4000):
    """Fully-vectorised Bz from a polyline (broadcast over segments x obs).

    ``bz_at`` (-> ``h_segments_batch``) loops over segments in Python (~6k
    iterations per call); a Gauss-Newton coil-distortion Jacobian calls the
    field 100x per iteration, so the Python loop is ~100x too slow.  This
    keeps only the z-component of H and broadcasts (segments x obs) in one
    NumPy expression, chunked over segments to cap transient memory.  Matches
    ``bz_at`` to machine precision.
    """
    path = np.asarray(path, float)
    p1, p2 = path[:-1], path[1:]
    dl = p2 - p1
    L = np.linalg.norm(dl, axis=1)
    ok = L > 1e-30
    el = np.zeros_like(dl)
    el[ok] = dl[ok] / L[ok, None]
    obs = np.asarray(obs, float)
    Bz = np.zeros(len(obs))
    for s in range(0, len(p1), chunk):
        e = el[s:s + chunk]
        a = p1[s:s + chunk]
        b = p2[s:s + chunk]
        r1 = obs[None, :, :] - a[:, None, :]    # (n, M, 3)
        r2 = obs[None, :, :] - b[:, None, :]
        cz = e[:, None, 0] * r1[:, :, 1] - e[:, None, 1] * r1[:, :, 0]
        cx = e[:, None, 1] * r1[:, :, 2] - e[:, None, 2] * r1[:, :, 1]
        cy = e[:, None, 2] * r1[:, :, 0] - e[:, None, 0] * r1[:, :, 2]
        d = np.sqrt(cx * cx + cy * cy + cz * cz)
        r1m = np.linalg.norm(r1, axis=2)
        r2m = np.linalg.norm(r2, axis=2)
        good = (d > 1e-30) & (r1m > 1e-30) & (r2m > 1e-30)
        dsafe = np.where(d > 0, d, 1.0)
        c1 = (r1 * e[:, None, :]).sum(axis=2) / np.where(r1m > 0, r1m, 1.0)
        c2 = (r2 * e[:, None, :]).sum(axis=2) / np.where(r2m > 0, r2m, 1.0)
        scale = np.where(good, current * INV_4PI / dsafe * (c1 - c2), 0.0)
        Bz += MU0 * (scale * (cz / dsafe)).sum(axis=0)
    return Bz


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


# --------------------------------------------------------------------------
# Shim-loop compensation of the single-stroke residual (the iterative method
# that pushes the EASY-tier uniform-Bz single stroke into the ~100-200 ppm
# class).  Pipeline: high-order FE psi -> clean contours -> single stroke ->
# Path-A (re-solve FE, monotone) -> a FEW independent shim loops via LS-OMP.
# --------------------------------------------------------------------------
def make_plane_shim_loops(plane_half, n_grid, frac=0.92, size_factor=1.6):
    """A grid of square current loops on the source plane -- the shim
    correction basis.  ``size_factor`` scales each loop's half-side relative
    to half the grid spacing: 1.0 = touching, >1 = OVERLAPPING (bigger,
    stronger field at the target plane -> the LS-OMP shim currents stay
    small).  Tiny non-overlapping loops are weak at the DSV and force huge
    shim currents, so the default overlaps.  Returns (5, 3) closed loops at
    z = 0."""
    cg = np.linspace(-plane_half * frac, plane_half * frac, n_grid)
    half = (cg[1] - cg[0]) / 2.0 if n_grid > 1 else plane_half * frac
    s = half * size_factor
    loops = []
    for cx in cg:
        for cy in cg:
            loops.append(np.array([
                [cx - s, cy - s, 0.0], [cx + s, cy - s, 0.0],
                [cx + s, cy + s, 0.0], [cx - s, cy + s, 0.0],
                [cx - s, cy - s, 0.0]], dtype=float))
    return loops


def ls_omp_shim(A_fit, B_fit, B_chain_fit, k_max, mae_stop_ppm=0.0,
               max_current_ratio=None, I_w=1.0):
    """Order-Recursive / LS-OMP shim selection (monotone).

    Cancels the single-stroke residual ``r = B_fit - B_chain_fit`` with a few
    independent loops: at each step orthogonalise every candidate column
    against the current support, add the one whose normalised orthogonal
    residual reduction is largest, re-solve the whole-support least squares,
    repeat.

    Stopping (in priority order):
      - ``mae_stop_ppm`` > 0: stop once the fit-grid MAE (in ppm relative to
        the target) drops to it -- the design target.
      - ``max_current_ratio``: stop BEFORE the largest shim current exceeds
        this multiple of ``I_w`` (guards the overfit blow-up that happens
        when too many near-collinear loops are added).
      - else ``k_max`` loops.

    Returns ``(support, I_shim, mae_curve_ppm)``.
    """
    r0 = B_fit - B_chain_fit
    denom = float(np.sum(np.abs(B_fit))) + 1e-30
    support, I_S, mae_curve = [], np.zeros(0), []
    for _ in range(int(k_max)):
        resid = r0 - (A_fit[:, support] @ I_S if support else 0.0)
        if support:
            Q, _ = np.linalg.qr(A_fit[:, support])
            A_perp = A_fit - Q @ (Q.T @ A_fit)
        else:
            A_perp = A_fit
        score = np.abs(A_perp.T @ resid) / (np.linalg.norm(A_perp, axis=0)
                                            + 1e-30)
        if support:
            score[support] = -1.0
        cand = int(np.argmax(score))
        trial = support + [cand]
        I_trial, _, _, _ = np.linalg.lstsq(A_fit[:, trial], r0, rcond=None)
        # guard against the overfit current blow-up
        if (max_current_ratio is not None
                and float(np.max(np.abs(I_trial))) / (abs(I_w) + 1e-30)
                > max_current_ratio and support):
            break
        support, I_S = trial, I_trial
        mae = 1e6 * float(np.sum(np.abs(r0 - A_fit[:, support] @ I_S))) / denom
        mae_curve.append(mae)
        if mae_stop_ppm > 0.0 and mae <= mae_stop_ppm:
            break
    return support, I_S, mae_curve


def coil_distort_3d(path0, fit_obs, B_fit, eval_obs, B_eval, plane_half,
                    comps=(0, 1, 2), n_grid=6, n_iter=8, lam_disp_rel=0.3,
                    step=0.8, fd=1.0e-4):
    """Single-current "sheet-metal" (bankin-ho) coil distortion.

    The separate-feed LS-OMP shim (``ls_omp_shim``) cancels the single-stroke
    residual by adding INDEPENDENT loops with their own currents.  This
    routine instead keeps ONE series current and the stream-function contour
    LEVELS fixed, and BENDS the manufactured single-stroke wire ``path0`` with
    a smooth low-dimensional deformation field

        d(x, y) -> (delta_x, delta_y, delta_z)   (an ``n_grid x n_grid``
        bilinear control grid per active component -- the discrete realisation
        of an NGSolve VectorH1 mesh deformation)

    to minimise the Bz error over the DSV.  Gauss-Newton with a displacement
    penalty RELATIVE to ``mean(diag(JtJ))`` keeps the bend manufacturable:

        ``lam_disp_rel`` large  -> small bends (conservative, 3D-printable)
        ``lam_disp_rel`` small  -> aggressive bends, lower MAE

    ``comps`` selects the deformed axes: ``(0, 1, 2)`` = full sheet-metal,
    ``(2,)`` = pure out-of-plane bend, ``(0, 1)`` = in-plane reroute only.
    The series current is re-fit (optimal single current) at every step.

    Returns ``(path_best, I_w_best, mae0_ppm, mae_best_ppm, max_disp_mm,
    mae_curve_ppm)``.
    """
    comps = tuple(comps)
    cg = np.linspace(-plane_half, plane_half, n_grid)
    hc = cg[1] - cg[0] if n_grid > 1 else 2 * plane_half
    nc = n_grid * n_grid
    ndof = nc * len(comps)

    def interp(C, xy):
        x, y = xy[:, 0], xy[:, 1]
        fx = np.clip((x - cg[0]) / hc, 0, n_grid - 1.001)
        fy = np.clip((y - cg[0]) / hc, 0, n_grid - 1.001)
        i = fx.astype(int)
        j = fy.astype(int)
        a = fx - i
        b = fy - j
        return (C[i, j] * (1 - a) * (1 - b) + C[i + 1, j] * a * (1 - b)
                + C[i, j + 1] * (1 - a) * b + C[i + 1, j + 1] * a * b)

    def deform(dofs):
        w = path0.copy()
        xy = path0[:, :2]
        for k, ax in enumerate(comps):
            w[:, ax] += interp(dofs[k * nc:(k + 1) * nc].reshape(n_grid,
                                                                 n_grid), xy)
        return w

    def fit_Iw_from(Bu):
        d = float(np.dot(Bu, Bu))
        return float(np.dot(Bu, B_fit) / d) if d > 0 else 0.0

    def fit_Iw(w):
        return fit_Iw_from(bz_fast(w, 1.0, fit_obs))

    den_e = float(np.sum(np.abs(B_eval))) + 1e-30

    def mae_eval(w, Iw):
        return 1e6 * float(np.sum(np.abs(Iw * bz_fast(w, 1.0, eval_obs)
                                          - B_eval))) / den_e

    Iw0 = fit_Iw(path0)
    mae0 = mae_eval(path0, Iw0)
    dofs = np.zeros(ndof)
    best = (mae0, path0.copy(), Iw0, 0.0)
    mae_curve = [mae0]
    for _ in range(int(n_iter)):
        w = deform(dofs)
        Iw = fit_Iw(w)
        base = Iw * bz_fast(w, 1.0, fit_obs)
        r = B_fit - base
        J = np.empty((len(fit_obs), ndof))
        for dd in range(ndof):
            d2 = dofs.copy()
            d2[dd] += fd
            Bu = bz_fast(deform(d2), 1.0, fit_obs)   # one field eval per DOF
            J[:, dd] = (fit_Iw_from(Bu) * Bu - base) / fd
        JTJ = J.T @ J
        md = float(np.mean(np.diag(JTJ))) + 1e-30
        lam = 0.03 * md                       # LM step damping
        lam_disp = lam_disp_rel * md          # displacement penalty
        step_dofs = np.linalg.solve(
            JTJ + (lam + lam_disp) * np.eye(ndof),
            J.T @ r - lam_disp * dofs)
        dofs = dofs + step * step_dofs
        w = deform(dofs)
        Iw = fit_Iw(w)
        m = mae_eval(w, Iw)
        mae_curve.append(m)
        max_disp = float(np.max(np.abs(dofs))) * 1e3
        if m < best[0]:
            best = (m, w.copy(), Iw, max_disp)
    mae_best, path_best, Iw_best, disp_best = best
    return path_best, Iw_best, mae0, mae_best, disp_best, mae_curve


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
    # shim-loop compensation + honest dense-grid ppm evaluation
    ap.add_argument("--shim-loops", type=int, default=0,
                    help="after the single stroke (+ Path-A), add up to this "
                         "many independent shim correction loops via LS-OMP "
                         "to drive the residual into the ~100-200 ppm class. "
                         "0 = off.")
    ap.add_argument("--shim-tol-ppm", type=float, default=0.0,
                    help="stop adding shim loops once the dense-grid MAE drops "
                         "to this many ppm (e.g. 200). 0 = use --shim-loops.")
    ap.add_argument("--shim-grid", type=int, default=15,
                    help="shim-loop basis = shim_grid x shim_grid small square "
                         "loops tiling the source plane")
    ap.add_argument("--eval-n", type=int, default=21,
                    help="dense eval grid (eval_n x eval_n over the target "
                         "patch) for HONEST ppm reporting (the n_target grid "
                         "is the fit grid and over-reports accuracy)")
    # single-current "sheet-metal" (bankin-ho) coil distortion -- an ALTERNATIVE
    # to --shim-loops that keeps ONE series current (no extra feeds) and bends
    # the manufactured wire to cancel the single-stroke residual.
    ap.add_argument("--distort", action="store_true",
                    help="bend the single-stroke wire with a smooth 3D "
                         "deformation field (ONE current, contour levels "
                         "fixed) to cancel the single-stroke residual -- the "
                         "single-current alternative to separate-feed shims")
    ap.add_argument("--distort-comps", choices=["xyz", "xy", "z"],
                    default="xyz",
                    help="deformed axes: xyz = full sheet-metal (best), "
                         "z = pure out-of-plane bend, xy = in-plane reroute")
    ap.add_argument("--distort-grid", type=int, default=6,
                    help="deformation control grid (distort_grid x "
                         "distort_grid bilinear field per active component)")
    ap.add_argument("--distort-iter", type=int, default=10,
                    help="Gauss-Newton iterations for the distortion")
    ap.add_argument("--distort-penalty", type=float, default=0.3,
                    help="displacement penalty (relative to mean diag(JtJ)): "
                         "larger = smaller bends (3D-printable), smaller = "
                         "more aggressive bends + lower MAE. On the default "
                         "benchmark (flat ~12k ppm): 1.0 -> ~17mm/2000ppm, "
                         "0.3 -> ~24mm/1000ppm, 0.1 -> ~29mm/600ppm")
    ap.add_argument("--distort-fit-n", type=int, default=9,
                    help="distortion fit grid (distort_fit_n^2 obs); the "
                         "displacement penalty regularises the underdetermined "
                         "fit so a modest grid is fine and ~10x faster than a "
                         "dense one")
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

    # ------- [8] HONEST dense-grid evaluation + shim-loop compensation ----
    # The [6] RMS is on the n_target FIT grid (over-reports -- the SF solve is
    # massively underdetermined and over-fits those points).  Evaluate the
    # single-stroke field on a DENSE eval grid for the true degradation, and
    # optionally compensate the residual with a few LS-OMP shim loops.
    ge = np.linspace(-args.target_half, args.target_half, args.eval_n)
    Xe, Ye = np.meshgrid(ge, ge, indexing="ij")
    eval_obs = np.column_stack([Xe.ravel(), Ye.ravel(),
                                args.target_z * np.ones(Xe.size)])
    B_eval = args.B0 * np.ones(eval_obs.shape[0])

    def mae_ppm(field):
        return 1e6 * float(np.sum(np.abs(field - B_eval))
                           / (np.sum(np.abs(B_eval)) + 1e-30))

    Bz_eval = I_w * bz_at(path, 1.0, eval_obs)
    mae0 = mae_ppm(Bz_eval)
    print(f"[8] HONEST eval ({args.eval_n}x{args.eval_n} dense grid):"
          f" single-stroke MAE = {mae0:.0f} ppm  (fit-grid RMS {rms*1e3:.1f} mppm"
          f" over-reports)")

    if args.shim_loops > 0 or args.shim_tol_ppm > 0.0:
        # fit shims on a moderately dense grid (between fit and eval), then
        # report on the dense eval grid (no overfit illusion).
        gf = np.linspace(-args.target_half, args.target_half,
                         max(11, args.eval_n - 6))
        Xf, Yf = np.meshgrid(gf, gf, indexing="ij")
        fit_obs = np.column_stack([Xf.ravel(), Yf.ravel(),
                                   args.target_z * np.ones(Xf.size)])
        B_fit = args.B0 * np.ones(fit_obs.shape[0])
        loops = make_plane_shim_loops(args.plane_half, args.shim_grid)
        A_fit = np.column_stack([bz_at(lp, 1.0, fit_obs) for lp in loops])
        A_eval = np.column_stack([bz_at(lp, 1.0, eval_obs) for lp in loops])
        B_chain_fit = I_w * bz_at(path, 1.0, fit_obs)
        # cap k_max generously; the real stop is the ppm target or the
        # current-blow-up guard (keep shim currents <= 10x I_w).
        k_max = args.shim_loops if args.shim_loops > 0 else 40
        support, I_shim, mae_curve = ls_omp_shim(
            A_fit, B_fit, B_chain_fit, k_max,
            mae_stop_ppm=args.shim_tol_ppm, max_current_ratio=10.0, I_w=I_w)
        Bz_eval_corr = Bz_eval + A_eval[:, support] @ I_shim
        stop = (f" (stopped at <= {args.shim_tol_ppm:.0f} ppm)"
                if args.shim_tol_ppm > 0 else "")
        print(f"    + {len(support)} LS-OMP shim loops"
              f" ({len(support)+1} total current feeds){stop}:"
              f" eval MAE = {mae0:.0f} -> {mae_ppm(Bz_eval_corr):.0f} ppm"
              f"  (shim |I|/I_w max = "
              f"{np.max(np.abs(I_shim))/(abs(I_w)+1e-30):.2f})")
        # convergence curve (fit-grid MAE per added shim)
        marks = sorted(set([0, 1, 4, 9, len(mae_curve) - 1]))
        pts = ", ".join(f"+{m+1}:{mae_curve[m]:.0f}p"
                        for m in marks if 0 <= m < len(mae_curve))
        print(f"    shim convergence (fit MAE): {pts}")

    # ------- [9] single-current sheet-metal coil distortion -------------
    # The ALTERNATIVE to [8]'s separate-feed shims: keep ONE series current
    # and the contour LEVELS fixed, BEND the manufactured wire (NGSolve-style
    # mesh deformation) to cancel the single-stroke residual.  One feed,
    # 3D-printable.
    distorted_path = None
    if args.distort:
        # fit on a modest grid (the displacement penalty + coarse control grid
        # regularise the underdetermined fit; the honest dense eval grid is
        # what's tracked for "best"); report on the dense eval grid.
        gd = np.linspace(-args.target_half, args.target_half, args.distort_fit_n)
        Xd, Yd = np.meshgrid(gd, gd, indexing="ij")
        distort_fit_obs = np.column_stack([Xd.ravel(), Yd.ravel(),
                                           args.target_z * np.ones(Xd.size)])
        distort_B_fit = args.B0 * np.ones(distort_fit_obs.shape[0])
        comps = {"xyz": (0, 1, 2), "xy": (0, 1), "z": (2,)}[args.distort_comps]
        (distorted_path, I_w_d, mae0_d, mae_best_d, max_disp_mm,
         mae_curve_d) = coil_distort_3d(
            path, distort_fit_obs, distort_B_fit, eval_obs, B_eval,
            args.plane_half, comps=comps, n_grid=args.distort_grid,
            n_iter=args.distort_iter, lam_disp_rel=args.distort_penalty)
        print(f"[9] sheet-metal coil distortion ({args.distort_comps}, "
              f"1 current, {args.distort_grid}x{args.distort_grid} ctrl grid, "
              f"penalty={args.distort_penalty:g}):")
        print(f"    eval MAE = {mae0_d:.0f} -> {mae_best_d:.0f} ppm"
              f"  (SINGLE current I_w = {I_w_d:.4g} A, no extra feeds;"
              f" max bend = {max_disp_mm:.1f} mm)")
        marks = sorted(set([0, 1, 2, len(mae_curve_d) // 2,
                            len(mae_curve_d) - 1]))
        pts = ", ".join(f"it{m}:{mae_curve_d[m]:.0f}p"
                        for m in marks if 0 <= m < len(mae_curve_d))
        print(f"    distortion convergence (eval MAE): {pts}")

    # ------- plot -------------------------------------------------------
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        plt.rcParams["pdf.fonttype"] = 42
        fig = plt.figure(figsize=(13, 4))
        ax1 = fig.add_subplot(1, 3, 1)
        cf = ax1.contourf(g_sample * 1e3, g_sample * 1e3, psi_grid.T,
                          levels=20, cmap="RdBu_r")
        lo, hi = psi_grid.min(), psi_grid.max()
        if hi > lo:
            lv = lo + (np.arange(args.nlevels) + 0.5) * (hi - lo) / args.nlevels
            ax1.contour(g_sample * 1e3, g_sample * 1e3, psi_grid.T,
                        levels=lv, colors="k", linewidths=0.5)
        ax1.set_xlabel("$x$ (mm)"); ax1.set_ylabel("$y$ (mm)")
        ax1.set_aspect("equal")
        ax1.text(0.02, 0.96, "(a)", transform=ax1.transAxes, fontsize=10)
        fig.colorbar(cf, ax=ax1, fraction=0.045)

        ax2 = fig.add_subplot(1, 3, 2)
        if distorted_path is not None:
            # draw the distorted wire as ONE CONTINUOUS LINE coloured by z.
            # It is still a SINGLE STROKE -- the deformation only displaces
            # vertices, never splits the path; the green/black markers are the
            # two terminals (current enters at one, exits at the other, spiral
            # through every loop via the radial rungs).  Colour = the
            # out-of-plane sheet-metal bend.
            from matplotlib.collections import LineCollection
            P = distorted_path
            pts = (P[:, :2] * 1e3).reshape(-1, 1, 2)
            segs = np.concatenate([pts[:-1], pts[1:]], axis=1)
            zc = P[:, 2] * 1e3
            lc = LineCollection(segs, cmap="coolwarm", linewidths=0.7)
            lc.set_array(0.5 * (zc[:-1] + zc[1:]))
            ax2.add_collection(lc)
            ax2.plot(P[0, 0] * 1e3, P[0, 1] * 1e3, "go", ms=6,
                     label="terminal (in)")
            ax2.plot(P[-1, 0] * 1e3, P[-1, 1] * 1e3, "ks", ms=6,
                     label="terminal (out)")
            lim = args.plane_half * 1.05e3
            ax2.set_xlim(-lim, lim); ax2.set_ylim(-lim, lim)
            fig.colorbar(lc, ax=ax2, fraction=0.045, label="bend $z$ (mm)")
            zmax = float(np.max(np.abs(P[:, 2]))) * 1e3
            ax2.legend(loc="upper right", fontsize=7)
        else:
            ax2.plot(path[:, 0] * 1e3, path[:, 1] * 1e3, "b-", lw=0.5)
        ax2.set_xlabel("$x$ (mm)"); ax2.set_ylabel("$y$ (mm)")
        ax2.set_aspect("equal")
        ax2.text(0.02, 0.96, "(b)", transform=ax2.transAxes, fontsize=10)
        ax2.grid(alpha=0.3)

        ax3 = fig.add_subplot(1, 3, 3)
        ax3.plot(zv * 1e3, Bz_axis * 1e3, "o-", ms=3, label="single-stroke coil")
        ax3.axhline(args.B0 * 1e3, color="k", ls="--",
                    label=f"target {args.B0*1e3:.2f} mT")
        ax3.axvline(args.target_z * 1e3, color="g", ls=":",
                    alpha=0.5, label=f"target z={args.target_z*1e3:.0f}mm")
        ax3.set_xlabel("$z$ (mm)"); ax3.set_ylabel("$B_z$ (mT)")
        ax3.text(0.02, 0.96, "(c)", transform=ax3.transAxes, fontsize=10)
        ax3.legend(loc="best", fontsize=8); ax3.grid(alpha=0.3)

        out = HERE / "demo_planar_uniform_fem_psi.png"
        fig.tight_layout()
        fig.savefig(out, dpi=110)
        print(f"Saved plot: {out.name}")
    except ImportError:
        print("(matplotlib not installed; skipped plot)")


if __name__ == "__main__":
    main()
