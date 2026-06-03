"""SF coil design: H1 FE psi + advanced regularization + surface deformation.

Extends ``demo_planar_uniform_fem_psi.py`` with four research-mode knobs the
user requested 2026-05-30:

  1. ``--regularize h1_sigma --sigma-cf EXPR``
     Conductivity-weighted H1 seminorm:  min INT (1/sigma(x)) |grad psi|^2
     subject to A psi = B (exact).  sigma(x) is a CoefficientFunction
     expression in x, y (e.g. ``"1.0 + 100*((x*x+y*y) > 0.04)"`` to forbid
     current outside the inner disc of radius 0.2 m).  Physical: TRUE ohmic
     dissipation for non-uniform conductivity (litz wire, mask).

  2. ``--regularize inductance_diag``
     LUMPED self-inductance proxy:  min sum_i (mu0/4pi) * |grad phi_i|^2
     * length_i, where length_i ~ sqrt(element area / pi) is a per-DOF
     characteristic scale.  This is a *diagonal* approximation to the true
     self-inductance form L = (mu0/4pi) INT INT grad psi(x).grad psi(y)/|x-y|;
     the full off-diagonal form needs surface BEM (ngsolve.bem MaxwellSL on
     the source plane embedded as a thin 3D shell -- deferred, see TODO).

  3. ``--regularize linf --jmax VAL``
     Sequential quadratic programming (scipy SLSQP) with explicit per-vertex
     constraint ``|grad psi|(vertex) <= jmax``.  Caps peak surface current
     density -> avoids hot spots.  Slower than the closed-form L2/H1 solves;
     scales as O(ndof^2) per outer SLSQP iteration.

  4. ``--deform --deform-trials N --deform-params SET``
     Outer Optuna CMA-ES loop over surface-deformation parameters.  Each
     trial re-meshes a parameterised surface (planar + Gaussian bump or
     z-offset), runs the inner SF solve + single-stroke spiral + field
     evaluation, returns the cost (default = target-plane RMS).
     ``deform-params`` chooses what to vary:
       - ``zoff``: 1 param, source plane z-offset (translate base plane)
       - ``bump``: 3 params, Gaussian bump amplitude + center (x, y)
       - ``zoff+bump``: 4 params combined

Combine freely: ``--regularize h1_sigma --deform --deform-params bump``
optimises the bump SHAPE while solving the inner SF design with
sigma-weighted regularisation.

Run examples:

    # 1/sigma weighting: forbid current in the outer ring
    python demo_planar_uniform_fem_psi_advanced.py \\
        --regularize h1_sigma --sigma-cf "1.0 / (1.0 + 100*((x*x+y*y) > 0.04))"

    # Inductance-diagonal proxy
    python demo_planar_uniform_fem_psi_advanced.py --regularize inductance_diag

    # L-infinity peak limit
    python demo_planar_uniform_fem_psi_advanced.py --regularize linf --jmax 50.0

    # Surface deformation with Gaussian bump (3-param CMA-ES search)
    python demo_planar_uniform_fem_psi_advanced.py --deform \\
        --deform-params bump --deform-trials 20
"""
import argparse
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "src"))

from radia.biot_savart import h_segments_batch, MU0

from demo_planar_uniform_coil import (
    contour_polylines_xy, single_stroke_spiral_xy, chain_xy_to_3d,
)
from demo_planar_uniform_fem_psi import sample_psi_grid


# --------------------------------------------------------------------------
# Biot-Savart Bz evaluation -- DEFORMED surface aware.  The source surface
# can be at z = surface_z(x, y), not just z = 0.  All references in the
# advanced demo go through this function so deformation is transparent.
# --------------------------------------------------------------------------
def bz_at_deformed(path_3d, current, obs):
    """Bz from a 3D wire path at given observation points (Tesla)."""
    segs = np.stack([path_3d[:-1], path_3d[1:]], axis=1)
    H = h_segments_batch(segs, obs, current=current)
    return MU0 * H[:, 2]


# --------------------------------------------------------------------------
# 1. Build FE matrix for arbitrary surface_z(x, y) -- the deformation knob
# --------------------------------------------------------------------------
def build_fem_matrix_deformed(plane_half, maxh, order, dirichlet_bc, targets,
                              surface_z_fn=None):
    """Mesh the (x, y) plane and build M x ndof Biot-Savart matrix A.

    surface_z_fn : callable(x, y) -> z  (default: z=0 i.e. flat plane).
        The source surface is parameterised by (x, y) in the meshed plane,
        with physical z = surface_z_fn(x, y).  Gradients are taken w.r.t.
        the parametric (x, y); for small deformation this is a good
        approximation of the true surface gradient.
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

    u, v = fes.TnT()
    M = len(targets)
    A = np.zeros((M, n_dof))
    mu0_4pi = 1.0e-7

    # surface_z(x, y) as CoefficientFunction
    if surface_z_fn is None:
        zs_cf = 0.0      # flat plane (constant CF)
    else:
        zs_cf = surface_z_fn(x_cf, y_cf)

    with TaskManager():
        for j, target in enumerate(targets):
            xt, yt, zt = float(target[0]), float(target[1]), float(target[2])
            dx_t = xt - x_cf
            dy_t = yt - y_cf
            dz_t = zt - zs_cf            # account for deformed source z
            r2 = dx_t * dx_t + dy_t * dy_t + dz_t * dz_t
            r3 = r2 * sqrt(r2)
            weight = -mu0_4pi / r3
            L = LinearForm(fes)
            L += weight * (dx_t * grad(v)[0] + dy_t * grad(v)[1]) * dx
            L.Assemble()
            A[j, :] = L.vec.FV().NumPy()
    return A, fes, mesh, n_free


# --------------------------------------------------------------------------
# 2. Advanced regularisation solvers
# --------------------------------------------------------------------------
def solve_l2_aca(A, B, fes, aca_eps=1.0e-10):
    """Min ||psi||_2 s.t. A psi = B via HACApK (ACA+)+TSVD.

    Equivalent to ``np.linalg.lstsq(A_free, B)`` (Euclidean min-Euclidean-
    norm solution) but routed through the kernel-agnostic ACA+ pipeline
    in ``radia.stream_function``.  At M = 25 / N_free ~ 1773 the direct
    lstsq is faster (~10 ms vs ~50 ms for ACA+), but this path:

      1. Validates that ACA+TSVD works on FE-direct matrices (callback
         contract carries over from basis-loop case).
      2. Caches the TSVD factorisation for Path-A iteration re-use.
      3. Is the foundation for material-kernel acceleration (Radia MMM,
         shielded coil) where each matrix entry is expensive.

    Note: this is the L2 / Euclidean min-norm path.  The H1 min-seminorm
    path would need an A . S^(-1/2) preconditioner (Cholesky factor of
    the stiffness matrix), which is O(ndof^3) and only useful for
    very large M.
    """
    from radia.stream_function import aca_tsvd, pseudo_inverse_solve
    free = np.array(fes.FreeDofs())
    free_idx = np.where(free)[0]
    A_free = A[:, free_idx]
    M, N = A_free.shape

    def entry(i, j):
        return float(A_free[i, j])

    res = aca_tsvd(M, N, entry, modes=M, kmax=min(M, N),
                   aca_eps=aca_eps, method=3)
    psi_free = pseudo_inverse_solve(res, B, k_mode=res.modes)
    psi = np.zeros(A.shape[1])
    psi[free_idx] = psi_free
    return psi, res


def solve_h1(A, B, fes, sigma_cf_expr=None):
    """Min INT (1/sigma) |grad psi|^2 s.t. A psi = B (exact).

    sigma_cf_expr : Python expression string in (x, y) for sigma(x,y);
        None -> uniform sigma=1 (= standard H1 seminorm).
    """
    from ngsolve import (BilinearForm, grad, dx, TaskManager,
                          x as x_cf, y as y_cf, IfPos, exp, sqrt as ngs_sqrt,
                          sin, cos, log)
    free = np.array(fes.FreeDofs())
    free_idx = np.where(free)[0]
    A_free = A[:, free_idx]

    u, v = fes.TnT()
    a_stiff = BilinearForm(fes, symmetric=True)
    if sigma_cf_expr is None:
        a_stiff += grad(u) * grad(v) * dx
    else:
        sigma_cf = eval(sigma_cf_expr,
                        {"x": x_cf, "y": y_cf, "IfPos": IfPos,
                         "exp": exp, "sqrt": ngs_sqrt, "sin": sin, "cos": cos,
                         "log": log, "__builtins__": {}})
        # 1/sigma weighting
        a_stiff += (1.0 / sigma_cf) * grad(u) * grad(v) * dx
    a_stiff += 1.0e-12 * u * v * dx     # tiny mass for SPD
    with TaskManager():
        a_stiff.Assemble()
    S = np.array(a_stiff.mat.ToDense())
    S_free = S[np.ix_(free_idx, free_idx)]

    # Lagrangian: psi = S^-1 A^T (A S^-1 A^T)^-1 B
    S_inv_AT = np.linalg.solve(S_free, A_free.T)
    small = A_free @ S_inv_AT
    lam = np.linalg.solve(small, B)
    psi_free = S_inv_AT @ lam
    psi = np.zeros(A.shape[1])
    psi[free_idx] = psi_free
    return psi


def solve_inductance_diagonal(A, B, fes):
    """Lumped self-inductance proxy:
       min sum_i w_i (grad psi_i)^2, with w_i ~ sqrt(area_i/pi).

    The TRUE inductance form is L = (mu0/4pi) INT INT grad psi(x) . grad psi(y)
    / |x-y| dA(x) dA(y), a dense N x N coupling matrix.  Full assembly needs
    surface BEM (ngsolve.bem MaxwellSL on a thin 3D shell, deferred).  This
    diagonal proxy weights each gradient DOF by a local characteristic length
    sqrt(area/pi) -- captures the leading-order self-energy contribution and is
    cheap (same cost as H1).  Use as a "physics-flavoured" smoothing rather
    than a true inductance optimum.
    """
    from ngsolve import (BilinearForm, grad, dx, TaskManager, sqrt as ngs_sqrt,
                          CoefficientFunction)
    free = np.array(fes.FreeDofs())
    free_idx = np.where(free)[0]
    A_free = A[:, free_idx]

    # Build a piecewise-constant CF that returns sqrt(element area / pi).
    # For triangle area we can use the element measure: integrate(1) per element.
    # Simple proxy: use mesh.h (= per-element characteristic length).
    mesh = fes.mesh
    # Element-wise characteristic length: take sqrt of triangle area
    from ngsolve import specialcf
    h_cf = specialcf.mesh_size       # element diameter (~ characteristic length)
    weight_cf = h_cf                  # length scale per element

    u, v = fes.TnT()
    a_stiff = BilinearForm(fes, symmetric=True)
    a_stiff += weight_cf * grad(u) * grad(v) * dx
    a_stiff += 1.0e-12 * u * v * dx
    with TaskManager():
        a_stiff.Assemble()
    S = np.array(a_stiff.mat.ToDense())
    S_free = S[np.ix_(free_idx, free_idx)]

    S_inv_AT = np.linalg.solve(S_free, A_free.T)
    small = A_free @ S_inv_AT
    lam = np.linalg.solve(small, B)
    psi_free = S_inv_AT @ lam
    psi = np.zeros(A.shape[1])
    psi[free_idx] = psi_free
    return psi


def solve_linf(A, B, fes, jmax, init_psi=None, tol=1.0e-6, max_iter=80):
    """Solve A psi = B (approximately) s.t. max_vertex |grad psi| <= jmax.

    Uses scipy SLSQP on the free-DOF subspace.  Objective = least-squares
    residual ||A_free psi - B||^2 + tiny H1 smoothing.  Constraint =
    per-vertex |grad psi|^2 <= jmax^2 (one constraint per mesh vertex).

    Slow (O(ndof^2) per SLSQP iteration) -- only meaningful for moderate
    ndof (a few hundred).  Use H1 order=1 for L-infinity studies.

    Prefer ``solve_linf_irls`` (Lawson IRLS on the cached
    ``RegularizedTSVD``): faster and scales to full FE meshes.
    """
    from scipy.optimize import minimize
    from ngsolve import GridFunction, grad, x as x_cf, y as y_cf

    free = np.array(fes.FreeDofs())
    free_idx = np.where(free)[0]
    A_free = A[:, free_idx]
    n_free = len(free_idx)

    # Precompute basis-function gradient at each vertex of the mesh.
    # For each vertex v, store row e_v such that
    #   |grad psi(v)|^2 = (e_v_x . psi)^2 + (e_v_y . psi)^2
    # where e_v_x, e_v_y are length-n_free row vectors.
    mesh = fes.mesh
    n_verts = mesh.nv
    Gx = np.zeros((n_verts, n_free))
    Gy = np.zeros((n_verts, n_free))
    gf = GridFunction(fes)
    grad_gf = grad(gf)
    for vk in range(n_verts):
        # set ψ = e_i (i-th free DOF) and read grad(ψ)(vertex k)
        pass
    # Above loop is too slow per-vertex.  Build G via finite-difference on
    # each basis function: set GridFunction = unit basis, sample grad at
    # vertex coordinates.
    for ii, gi in enumerate(free_idx):
        gf.vec.FV().NumPy()[:] = 0.0
        gf.vec.FV().NumPy()[int(gi)] = 1.0
        for vk in range(n_verts):
            vert = mesh.vertices[vk]
            mp = mesh(*vert.point)
            gx, gy = grad_gf(mp)[:2]
            Gx[vk, ii] = float(gx)
            Gy[vk, ii] = float(gy)

    # objective: ||A_free psi - B||^2 + small smoothing
    AtA = A_free.T @ A_free
    AtB = A_free.T @ B
    eps_smooth = 1.0e-6
    Smooth = eps_smooth * (Gx.T @ Gx + Gy.T @ Gy)

    def fun(psi_free):
        return float(psi_free @ AtA @ psi_free
                      - 2.0 * AtB @ psi_free
                      + psi_free @ Smooth @ psi_free)

    def fun_grad(psi_free):
        return 2.0 * (AtA @ psi_free - AtB + Smooth @ psi_free)

    # constraint: jmax^2 - |grad psi(vk)|^2 >= 0 for each vk
    jmax2 = jmax * jmax

    def cstr_fn(psi_free):
        gx_v = Gx @ psi_free
        gy_v = Gy @ psi_free
        return jmax2 - (gx_v ** 2 + gy_v ** 2)

    def cstr_jac(psi_free):
        gx_v = Gx @ psi_free
        gy_v = Gy @ psi_free
        return -2.0 * (gx_v[:, None] * Gx + gy_v[:, None] * Gy)

    constraints = [{"type": "ineq", "fun": cstr_fn, "jac": cstr_jac}]
    if init_psi is None:
        x0 = np.zeros(n_free)
    else:
        x0 = init_psi[free_idx].copy()
    result = minimize(fun, x0, jac=fun_grad, method="SLSQP",
                      constraints=constraints,
                      options={"maxiter": max_iter, "ftol": tol, "disp": False})

    psi_free = result.x
    psi = np.zeros(A.shape[1])
    psi[free_idx] = psi_free
    return psi, result


def solve_linf_irls(A, B, fes, jmax=None, n_iter=12, weight_ratio_max=20.0,
                    damping=0.5, hot_quantile=0.85, hot_boost=3.0,
                    verbose=False):
    """L-infinity gradient cap via bounded-ratio active-set IRLS.

    Iteratively reweighted least squares on top of the cached
    ``RegularizedTSVD``.  Each iter solves the min-weighted-H1 problem
    subject to ``A psi = B`` exactly via the closed-form
    ``psi = (S^(k))^-1 V . W^-1 . Sigma^-1 . U^T B`` (the K x K matrix
    W = V^T (S^(k))^-1 V is rebuilt each iter; ACA+ runs ONCE).

    The reweighting is **bounded** to keep the linear system
    well-conditioned even when a few vertices reach near-zero gradient
    (which makes pure Lawson IRLS oscillate):

        iter k:
          peak  = max_v |grad psi^(k-1)|
          ratio = |grad psi^(k-1)|_v / peak     (clipped to [1/r, 1])
          base  = mix(uniform, ratio^2, damping)
          if |grad psi^(k-1)|_v > hot_quantile * peak:
              base[v] *= hot_boost              (active-set push)
          S^(k) = int base(x) |grad u . grad v| dA

    The clipping `weight_ratio_max` (default 20) bounds the
    conditioning of ``S^(k)`` relative to the uniform-weight S; the
    `damping` parameter mixes the new weight with the previous one
    (0 = no update, 1 = full Lawson); the `hot_boost` factor adds an
    active-set push so the vertices currently above ``hot_quantile``
    of the peak get extra penalty.

    Parameters
    ----------
    A : ndarray
        M x ndof Biot-Savart matrix (full FES DOF, Dirichlet not yet
        eliminated).
    B : ndarray, shape (M,)
        Target field.
    fes : NGSolve FESpace
        H1 source-plane space.
    jmax : float, optional
        Optional cap.  If set, ``psi`` is RESCALED at the end so
        ``max_v |grad psi| <= jmax`` exactly (the field at targets is
        scaled correspondingly; recover via the I_w best-fit-current
        step downstream).
    n_iter : int, optional
        IRLS iterations (default 12).
    weight_ratio_max : float, optional
        Maximum allowed weight ratio between any two vertices.  Caps
        conditioning of ``S^(k)`` (default 20).
    damping : float, optional
        New-weight mixing factor (default 0.5 = average of old and new).
    hot_quantile : float, optional
        Vertices with |grad| above ``hot_quantile * peak`` are flagged
        active and get ``hot_boost`` extra weight (default 0.85).
    hot_boost : float, optional
        Extra penalty factor on active vertices (default 3.0).
    verbose : bool, optional
        Print per-iter peak vs mean gradient.

    Returns
    -------
    (psi, info) : (ndarray shape (ndof,), dict)
        info has keys ``"peak_grad"``, ``"peak_history"``,
        ``"n_iter_done"``, ``"jmax_scale"``.
    """
    from ngsolve import (BilinearForm, GridFunction, grad, dx, TaskManager,
                         H1 as H1_space)
    from radia.stream_function import aca_tsvd, RegularizedTSVD

    free = np.array(fes.FreeDofs())
    free_idx = np.where(free)[0]
    A_free = A[:, free_idx]
    M_dim, N_dim = A_free.shape
    mesh = fes.mesh
    n_verts = mesh.nv

    # ACA+TSVD factorise A_free ONCE -- shared across all IRLS iters.
    def entry(i, j):
        return float(A_free[i, j])

    res = aca_tsvd(M_dim, N_dim, entry, modes=M_dim,
                   kmax=min(M_dim, N_dim), aca_eps=1.0e-10, method=3)

    vert_xy = np.array([list(v.point) for v in mesh.vertices])
    u_t, v_t = fes.TnT()
    gf = GridFunction(fes)
    fes_w = H1_space(mesh, order=1)        # weight CF lives here
    gf_w = GridFunction(fes_w)

    def assemble_weighted_stiffness(w_per_vertex):
        """``int w(x) |grad u| . |grad v| dA`` with w from vertex values."""
        gf_w.vec.FV().NumPy()[:n_verts] = w_per_vertex
        a_stiff = BilinearForm(fes, symmetric=True)
        a_stiff += gf_w * grad(u_t) * grad(v_t) * dx
        a_stiff += 1.0e-12 * u_t * v_t * dx
        with TaskManager():
            a_stiff.Assemble()
        S = np.array(a_stiff.mat.ToDense())
        return S[np.ix_(free_idx, free_idx)]

    def peak_gradient(psi_free):
        gf.vec.FV().NumPy()[:] = 0.0
        gf.vec.FV().NumPy()[free_idx] = psi_free
        grad_gf = grad(gf)
        gx = np.zeros(n_verts); gy = np.zeros(n_verts)
        for vk in range(n_verts):
            mp = mesh(*vert_xy[vk])
            gv = grad_gf(mp)
            gx[vk] = float(gv[0]); gy[vk] = float(gv[1])
        return np.sqrt(gx * gx + gy * gy)

    # --- IRLS with bounded-ratio + active-set push -----------------------
    w = np.ones(n_verts)                          # uniform start
    peak_hist = []
    info = {"peak_history": peak_hist, "n_iter_done": 0,
            "peak_grad": float("nan"), "jmax_scale": 1.0}
    psi_free = None

    for k in range(n_iter):
        S_free = assemble_weighted_stiffness(w)
        reg = RegularizedTSVD.from_stiffness(res, S_free)
        psi_free = reg.solve(B)
        grad_norm = peak_gradient(psi_free)
        peak = float(np.max(grad_norm))
        mean = float(np.mean(grad_norm))
        peak_hist.append(peak)
        info["n_iter_done"] = k + 1
        if verbose:
            print(f"     IRLS iter {k+1}: peak |grad| = {peak:.3e},"
                  f" mean = {mean:.3e}, peak/mean = {peak / (mean + 1e-30):.2f}")

        # New weight: bounded ratio + active-set boost
        rel = np.clip(grad_norm / (peak + 1e-30),
                      1.0 / weight_ratio_max, 1.0)
        w_new = rel ** 2.0
        # Active-set push on the top hot_quantile vertices
        thresh = float(np.quantile(grad_norm, hot_quantile))
        active = grad_norm >= thresh
        w_new = np.where(active, w_new * hot_boost, w_new)
        # Renormalise so the mean weight stays ~1 (prevents scale drift)
        w_new = w_new / np.mean(w_new)
        # Damped update
        w = (1.0 - damping) * w + damping * w_new

    info["peak_grad"] = peak_hist[-1] if peak_hist else float("nan")

    # Optional cap: scalar rescale so max |grad psi| <= jmax exactly.
    if jmax is not None and psi_free is not None and info["peak_grad"] > 0.0:
        if info["peak_grad"] > jmax:
            scale = float(jmax) / info["peak_grad"]
            psi_free = psi_free * scale
            info["jmax_scale"] = scale
            info["peak_grad"] = float(jmax)
            if verbose:
                print(f"     IRLS cap: scaled psi by {scale:.3e}"
                      f" -> max |grad| = {jmax:.3e}")
        else:
            if verbose:
                print(f"     IRLS cap not active: peak {info['peak_grad']:.3e}"
                      f" already <= jmax {jmax:.3e}")

    psi = np.zeros(A.shape[1])
    if psi_free is not None:
        psi[free_idx] = psi_free
    return psi, info


# --------------------------------------------------------------------------
# 3. Inner SF solve + chain + cost (used by --deform outer loop)
# --------------------------------------------------------------------------
def _compute_reg_norm(psi, A, fes, regularize, sigma_cf_expr=None):
    """Return ``psi_free^T S_free psi_free`` for the active regularisation.

    For ``l2*`` modes returns ``||psi_free||^2``; for ``linf*`` modes
    returns the H1 seminorm ``int |grad psi|^2 dA`` as a proxy (the
    L-inf cap is not a quadratic form).  Used by the constrained
    Optuna loop to evaluate the objective ``min psi^T S psi s.t.
    RMS <= eps``.
    """
    free = np.array(fes.FreeDofs())
    free_idx = np.where(free)[0]
    psi_free = psi[free_idx]
    if regularize in ("l2", "l2_aca"):
        return float(psi_free @ psi_free)
    # h1 / h1_sigma / inductance_diag / linf* all use a stiffness form
    from demo_planar_uniform_fem_psi import build_h1_stiffness
    if regularize == "h1_sigma" and sigma_cf_expr is not None:
        # Use the same 1/sigma weighting that was passed to solve_h1.
        # (Re-assembling here costs one BilinearForm + dense materialise;
        # ~50 ms at our N_free ~ 1773.)
        from ngsolve import (BilinearForm, grad, dx, TaskManager,
                              x as x_cf, y as y_cf, IfPos, exp,
                              sqrt as ngs_sqrt, sin, cos, log)
        u, v = fes.TnT()
        sigma_cf = eval(sigma_cf_expr,
                        {"x": x_cf, "y": y_cf, "IfPos": IfPos,
                         "exp": exp, "sqrt": ngs_sqrt, "sin": sin,
                         "cos": cos, "log": log, "__builtins__": {}})
        a_stiff = BilinearForm(fes, symmetric=True)
        a_stiff += (1.0 / sigma_cf) * grad(u) * grad(v) * dx
        a_stiff += 1.0e-12 * u * v * dx
        with TaskManager():
            a_stiff.Assemble()
        S = np.array(a_stiff.mat.ToDense())
        S_free = S[np.ix_(free_idx, free_idx)]
    elif regularize == "inductance_diag":
        from ngsolve import (BilinearForm, grad, dx, TaskManager, specialcf)
        u, v = fes.TnT()
        h_cf = specialcf.mesh_size
        a_stiff = BilinearForm(fes, symmetric=True)
        a_stiff += h_cf * grad(u) * grad(v) * dx
        a_stiff += 1.0e-12 * u * v * dx
        with TaskManager():
            a_stiff.Assemble()
        S = np.array(a_stiff.mat.ToDense())
        S_free = S[np.ix_(free_idx, free_idx)]
    else:
        # h1, linf, linf_irls -> plain H1 stiffness
        S_free = build_h1_stiffness(fes, free_idx)
    return float(psi_free @ S_free @ psi_free)


def solve_and_evaluate(plane_half, maxh, order, target_args, B0,
                       regularize="h1", sigma_cf_expr=None, jmax=None,
                       surface_z_fn=None, n_sample=81, nlevels=12,
                       compute_reg_norm=False):
    """Full inner-loop: mesh -> FE matrix -> solve -> chain -> Bz cost.

    Returns dict with cost (target-plane RMS), I_w, chain length, etc.
    If ``compute_reg_norm=True``, also returns ``reg_norm = psi^T S psi``
    in the dict (used by the constrained-Optuna ``--minimize-reg`` loop).
    """
    gt = np.linspace(-target_args["half"], target_args["half"],
                     target_args["n"])
    Xt, Yt = np.meshgrid(gt, gt, indexing="ij")
    targets = np.column_stack([Xt.ravel(), Yt.ravel(),
                                target_args["z"] * np.ones(Xt.size)])
    M = targets.shape[0]
    B_target = B0 * np.ones(M)

    A, fes, mesh, n_free = build_fem_matrix_deformed(
        plane_half, maxh, order, "boundary", targets, surface_z_fn)

    if regularize == "h1":
        psi = solve_h1(A, B_target, fes, sigma_cf_expr=None)
    elif regularize == "h1_sigma":
        psi = solve_h1(A, B_target, fes, sigma_cf_expr=sigma_cf_expr)
    elif regularize == "inductance_diag":
        psi = solve_inductance_diagonal(A, B_target, fes)
    elif regularize == "linf":
        psi_init = solve_h1(A, B_target, fes, sigma_cf_expr=None)
        psi, _ = solve_linf(A, B_target, fes, jmax, init_psi=psi_init)
    elif regularize == "linf_irls":
        # Bounded-ratio active-set IRLS on cached ACA+TSVD -- replaces
        # SLSQP for L-inf cap (faster, scales to full FE meshes)
        psi, _ = solve_linf_irls(A, B_target, fes, jmax=jmax,
                                 n_iter=12, verbose=False)
    elif regularize == "l2_aca":
        # ACA+TSVD-routed L2 min-norm.  Validates ACA+ on FE-direct matrix.
        psi, _ = solve_l2_aca(A, B_target, fes)
    else:
        # l2 = numpy lstsq (Euclidean min-norm)
        free_idx = np.where(np.array(fes.FreeDofs()))[0]
        psi_free, _, _, _ = np.linalg.lstsq(
            A[:, free_idx], B_target, rcond=None)
        psi = np.zeros(A.shape[1])
        psi[free_idx] = psi_free

    psi_grid, g_sample = sample_psi_grid(
        psi, fes, mesh, plane_half, n_sample)
    polylines, dI = contour_polylines_xy(psi_grid, g_sample, g_sample,
                                          nlevels)
    if not polylines:
        return {"cost": float("inf"), "rms": float("inf"), "I_w": 0.0,
                "length": 0.0, "n_wire": 0,
                "Bz_mean": 0.0, "p2p_mean": float("inf"),
                "reg_norm": float("inf"),
                "psi": psi, "fes": fes, "mesh": mesh,
                "polylines": [], "chain_xy": None, "path_xy": None}

    chain_xy = single_stroke_spiral_xy(polylines)
    path_xy = chain_xy_to_3d(chain_xy, z=0.0)
    # Lift chain to deformed surface: replace z = surface_z_fn(x, y)
    if surface_z_fn is not None:
        # Evaluate surface_z_fn numerically on the chain xy
        # (CF -> mesh point; for analytic fn we just call directly)
        for k in range(len(path_xy)):
            path_xy[k, 2] = _numeric_surface_z_fn(
                surface_z_fn, path_xy[k, 0], path_xy[k, 1])

    Bz_unit = bz_at_deformed(path_xy, 1.0, targets)
    denom = float(np.dot(Bz_unit, Bz_unit))
    I_w = float(np.dot(Bz_unit, B_target) / denom) if denom > 0 else 0.0
    Bz_dsv = I_w * Bz_unit
    rms = float(np.linalg.norm(Bz_dsv - B_target) /
                (np.linalg.norm(B_target) + 1e-30))
    Bz_mean = float(np.mean(Bz_dsv))
    p2p = (float(np.max(Bz_dsv)) - float(np.min(Bz_dsv))) / (abs(Bz_mean) + 1e-30)
    length = float(np.sum(np.linalg.norm(np.diff(path_xy, axis=0), axis=1)))
    reg_norm = (_compute_reg_norm(psi, A, fes, regularize, sigma_cf_expr)
                if compute_reg_norm else float("nan"))
    return {
        "cost": rms,
        "I_w": I_w,
        "length": length,
        "n_wire": len(polylines),
        "Bz_mean": Bz_mean,
        "p2p_mean": p2p,
        "rms": rms,
        "reg_norm": reg_norm,
        "psi": psi,
        "fes": fes,
        "mesh": mesh,
        "polylines": polylines,
        "chain_xy": chain_xy,
        "path_xy": path_xy,
    }


def _numeric_surface_z_fn(cf_or_fn, x, y):
    """Evaluate surface_z(x, y) numerically; the CF version is built later."""
    # If it's a Python callable returning a scalar, call directly.
    if callable(cf_or_fn):
        try:
            return float(cf_or_fn(x, y))
        except Exception:
            pass
    return 0.0


# --------------------------------------------------------------------------
# 4. Outer Optuna CMA-ES loop -- surface deformation search
# --------------------------------------------------------------------------
def run_deformation_search(args, trials=20):
    """Optuna outer loop over surface deformation parameters.

    Three modes, selected from ``args``:

      1. **Default** (no flag): single-objective CMA-ES minimising
         target-plane RMS.

      2. ``args.minimize_reg=True``: single-objective CMA-ES
         minimising the regularisation norm ``psi^T S psi`` subject
         to ``RMS <= args.eps_rms`` (penalty form).  This is the
         standard MRI-shim formulation -- "the most energy-efficient
         coil that meets the field-uniformity tolerance".  The
         constraint enters as

             cost = reg_norm + reg_scale * penalty_weight *
                    max(0, RMS/eps - 1)^2

         where ``reg_scale`` is set from the FIRST trial's reg_norm so
         the penalty has comparable magnitude to the reg term.

      3. ``args.pareto=True``: multi-objective NSGA-II returning
         ``(rms, reg_norm)`` tuples.  Optuna traces the Pareto front
         directly; access via ``study.best_trials`` (non-dominated
         set).  No penalty weights to tune.
    """
    import optuna
    optuna.logging.set_verbosity(optuna.logging.WARNING)
    target_args = {"half": args.target_half, "z": args.target_z,
                   "n": args.n_target}
    base_z = 0.0    # base plane z = 0

    minimize_reg = getattr(args, "minimize_reg", False)
    pareto_mode = getattr(args, "pareto", False)
    eps_rms = float(getattr(args, "eps_rms", 0.02))
    penalty_weight = float(getattr(args, "reg_penalty", 100.0))
    state = {"reg_scale": None,
             "best_feasible": (float("inf"), float("inf"))}  # (reg, rms)

    def surface_z_fn_from(deform_params):
        def surface_z_fn(x, y):
            from ngsolve import exp as ngs_exp
            is_cf = "CoefficientFunction" in type(x).__name__
            exp_fn = ngs_exp if is_cf else np.exp
            z_val = base_z
            if "zoff" in deform_params:
                z_val = z_val + deform_params["zoff"]
            if "bump" in deform_params:
                amp, cx, cy = deform_params["bump"]
                z_val = z_val + amp * exp_fn(
                    -((x - cx) ** 2 + (y - cy) ** 2) / 0.005)
            return z_val
        return surface_z_fn

    def evaluate_trial(trial):
        """Sample deformation params, run inner solve, return (rms, reg)."""
        deform_params = {}
        if "zoff" in args.deform_params:
            zoff = trial.suggest_float("zoff", -0.05, 0.05)
            deform_params["zoff"] = zoff
        if "bump" in args.deform_params:
            amp = trial.suggest_float("bump_amp", -0.05, 0.05)
            cx = trial.suggest_float("bump_cx", -0.15, 0.15)
            cy = trial.suggest_float("bump_cy", -0.15, 0.15)
            deform_params["bump"] = (amp, cx, cy)

        result = solve_and_evaluate(
            args.plane_half, args.maxh, args.order, target_args, args.B0,
            regularize=args.regularize, sigma_cf_expr=args.sigma_cf,
            jmax=args.jmax, surface_z_fn=surface_z_fn_from(deform_params),
            n_sample=args.n_sample, nlevels=args.nlevels,
            compute_reg_norm=(minimize_reg or pareto_mode))
        rms = float(result["rms"])
        reg = float(result.get("reg_norm", float("nan")))
        trial.set_user_attr("rms", rms)
        trial.set_user_attr("reg_norm", reg)
        return rms, reg

    if pareto_mode:
        # ---- Multi-objective NSGA-II ----
        def objective_mo(trial):
            rms, reg = evaluate_trial(trial)
            if not np.isfinite(reg):
                # Reject infeasible trials by penalising both objectives
                return float("inf"), float("inf")
            return rms, reg

        sampler = optuna.samplers.NSGAIISampler(
            seed=42, population_size=max(8, min(trials // 4, 24)))
        study = optuna.create_study(directions=["minimize", "minimize"],
                                     sampler=sampler)
        study.optimize(objective_mo, n_trials=trials)
        return study

    # ---- Single-objective CMA-ES (default or --minimize-reg) ----
    def objective_so(trial):
        rms, reg = evaluate_trial(trial)
        if not minimize_reg:
            return rms
        if not np.isfinite(reg):
            return float("inf")
        if state["reg_scale"] is None:
            state["reg_scale"] = max(reg, 1.0e-30)
        viol = max(0.0, rms / eps_rms - 1.0)
        cost = reg + state["reg_scale"] * penalty_weight * viol ** 2
        if rms <= eps_rms and reg < state["best_feasible"][0]:
            state["best_feasible"] = (reg, rms)
        trial.set_user_attr("feasible", rms <= eps_rms)
        return cost

    sampler = optuna.samplers.CmaEsSampler(seed=42)
    study = optuna.create_study(direction="minimize", sampler=sampler)
    study.optimize(objective_so, n_trials=trials)
    if minimize_reg:
        study.set_user_attr("best_feasible", state["best_feasible"])
        study.set_user_attr("eps_rms", eps_rms)
        study.set_user_attr("reg_scale", state["reg_scale"])
        study.set_user_attr("penalty_weight", penalty_weight)
    return study


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--plane-half", type=float, default=0.25)
    ap.add_argument("--maxh", type=float, default=0.025)
    ap.add_argument("--order", type=int, default=2)
    ap.add_argument("--target-half", type=float, default=0.05)
    ap.add_argument("--target-z", type=float, default=0.10)
    ap.add_argument("--n-target", type=int, default=5)
    ap.add_argument("--B0", type=float, default=1.0e-3)
    ap.add_argument("--n-sample", type=int, default=81)
    ap.add_argument("--nlevels", type=int, default=12)
    ap.add_argument("--regularize",
                    choices=["l2", "l2_aca", "h1", "h1_sigma",
                             "inductance_diag", "linf", "linf_irls"],
                    default="h1",
                    help="solver for the inner SF problem ('linf_irls' uses"
                         " Lawson IRLS on cached ACA+TSVD; faster than 'linf'"
                         " which uses scipy SLSQP)")
    ap.add_argument("--sigma-cf", type=str, default=None,
                    help="Python CF expression for sigma(x, y); used with"
                         " --regularize h1_sigma")
    ap.add_argument("--jmax", type=float, default=50.0,
                    help="L-inf cap on |grad psi| (use with --regularize linf)")
    ap.add_argument("--deform", action="store_true",
                    help="run Optuna CMA-ES on surface deformation params")
    ap.add_argument("--deform-params", type=str, default="bump",
                    help="comma-separated subset of {zoff, bump}")
    ap.add_argument("--deform-trials", type=int, default=20)
    # Reg-aware constrained optimisation: min psi^T S psi s.t. RMS <= eps
    ap.add_argument("--minimize-reg", action="store_true",
                    help="Optuna minimises the regularisation norm "
                         "psi^T S psi (= inductance proxy / surface "
                         "current L2 / ohmic dissipation, depending on "
                         "--regularize) subject to RMS <= --eps-rms. "
                         "Requires --deform.")
    ap.add_argument("--eps-rms", type=float, default=0.02,
                    help="target-plane RMS tolerance for --minimize-reg "
                         "(default 0.02 = 2 percent)")
    ap.add_argument("--reg-penalty", type=float, default=100.0,
                    help="penalty weight in the constrained-Optuna "
                         "objective (default 100). Scaled by first-trial "
                         "reg_norm so the units are comparable.")
    ap.add_argument("--pareto", action="store_true",
                    help="multi-objective NSGA-II returning (RMS, "
                         "reg_norm) tuples. Traces the Pareto front "
                         "directly; no penalty weights. Requires --deform.")
    args = ap.parse_args()
    if args.minimize_reg and not args.deform:
        ap.error("--minimize-reg requires --deform (it's an outer-loop "
                 "optimisation over deformation parameters)")
    if args.pareto and not args.deform:
        ap.error("--pareto requires --deform (NSGA-II searches over "
                 "deformation parameters)")
    if args.pareto and args.minimize_reg:
        ap.error("--pareto and --minimize-reg are mutually exclusive "
                 "(pareto is multi-objective, minimize-reg is "
                 "single-objective with penalty)")

    print(f"=== Advanced SF demo: regularize={args.regularize}, "
          f"deform={args.deform} ===")

    if args.deform:
        if args.pareto:
            print(f"[*] Outer NSGA-II: {args.deform_trials} trials on"
                  f" deform_params={args.deform_params!r},"
                  f" mode=Pareto(RMS, psi^T S psi)")
        elif args.minimize_reg:
            print(f"[*] Outer CMA-ES: {args.deform_trials} trials on"
                  f" deform_params={args.deform_params!r},"
                  f" mode=min(psi^T S psi) s.t. RMS <= {args.eps_rms:.2e}")
        else:
            print(f"[*] Outer CMA-ES: {args.deform_trials} trials on"
                  f" deform_params={args.deform_params!r}, mode=min(RMS)")
        t0 = time.time()
        study = run_deformation_search(args, trials=args.deform_trials)
        elapsed = time.time() - t0

        if args.pareto:
            print(f"\n=== NSGA-II finished in {elapsed:.1f} s ===")
            # Dedupe pareto trials: NSGA-II can evaluate the same params
            # across generations.  Group by (rms, reg) rounded to 4 sig figs.
            seen = set()
            pareto = []
            for t in sorted(study.best_trials, key=lambda t: t.values[0]):
                key = (round(t.values[0], 6), round(t.values[1], 2))
                if key in seen:
                    continue
                seen.add(key)
                pareto.append(t)
            n_dups = len(study.best_trials) - len(pareto)
            dup_tag = f" ({n_dups} duplicates removed)" if n_dups else ""
            print(f"  {len(pareto)} Pareto-optimal trials"
                  f" (of {len(study.trials)} total){dup_tag}:")
            print(f"  {'RMS':<12}{'reg_norm':<14}{'params'}")
            print(f"  {'-' * 60}")
            for t in pareto:
                rms_v, reg_v = t.values
                p_str = ", ".join(f"{k}={v:+.4f}" for k, v in t.params.items())
                print(f"  {rms_v*100:<10.2f}% {reg_v:<14.4e} {p_str}")

            # Save Pareto JSON + plot
            import json
            pareto_payload = {
                "n_trials": len(study.trials),
                "n_pareto": len(pareto),
                "elapsed_s": elapsed,
                "pareto": [
                    {"rms": t.values[0], "reg_norm": t.values[1],
                     "params": dict(t.params)} for t in pareto
                ],
                "all_trials": [
                    {"rms": t.user_attrs.get("rms"),
                     "reg_norm": t.user_attrs.get("reg_norm"),
                     "params": dict(t.params)}
                    for t in study.trials
                    if t.user_attrs.get("rms") is not None
                ],
            }
            out_json = HERE / "demo_pareto_results.json"
            with open(out_json, "w") as f:
                json.dump(pareto_payload, f, indent=2)
            print(f"\n  Saved: {out_json.name}")

            # Plot
            try:
                import matplotlib
                matplotlib.use("Agg")
                import matplotlib.pyplot as plt
                plt.rcParams["pdf.fonttype"] = 42
                fig, ax = plt.subplots(figsize=(7.0, 5.0))
                # All trials (dominated + Pareto)
                rms_all = [t.user_attrs.get("rms") for t in study.trials
                            if t.user_attrs.get("rms") is not None
                            and np.isfinite(t.user_attrs.get("rms"))]
                reg_all = [t.user_attrs.get("reg_norm") for t in study.trials
                            if t.user_attrs.get("rms") is not None
                            and np.isfinite(t.user_attrs.get("reg_norm"))]
                ax.scatter(np.array(rms_all) * 100, reg_all, s=18,
                           c="lightgray", label="dominated trials")
                # Pareto front (sorted by RMS)
                rms_p = np.array([t.values[0] for t in pareto]) * 100
                reg_p = np.array([t.values[1] for t in pareto])
                ax.plot(rms_p, reg_p, "o-", c="C0", lw=1.5, ms=7,
                        label="Pareto front")
                ax.set_xlabel("target-plane RMS  (%)")
                ax.set_ylabel(r"$\psi^T S \psi$  (regularisation norm)")
                ax.set_yscale("log" if reg_p.max() / max(reg_p.min(), 1e-30)
                              > 10 else "linear")
                ax.legend(loc="best", fontsize=9)
                ax.grid(alpha=0.3)
                out_png = HERE / "demo_pareto_plot.png"
                fig.tight_layout()
                fig.savefig(out_png, dpi=110)
                print(f"  Saved: {out_png.name}")
                plt.close(fig)
            except ImportError:
                print("  (matplotlib not installed; skipped plot)")

            # Pick the median Pareto trial for the downstream field/chain
            # evaluation (so the [6] block reports something interpretable)
            median_t = pareto[len(pareto) // 2]
            best = median_t.params
            deform_params = {}
            if "zoff" in best:
                deform_params["zoff"] = best["zoff"]
            if "bump_amp" in best:
                deform_params["bump"] = (best["bump_amp"], best["bump_cx"],
                                           best["bump_cy"])
            print(f"\n  (showing median Pareto trial: RMS={median_t.values[0]*100:.2f}%,"
                  f" reg={median_t.values[1]:.4e})")
        else:
            print(f"\n=== CMA-ES finished in {elapsed:.1f} s ===")
            print(f"  best cost = {study.best_value:.4e}")
            print(f"  best params = {study.best_params}")
            if args.minimize_reg:
                best_feasible = study.user_attrs.get(
                    "best_feasible", (float("inf"), float("inf")))
                n_feasible = sum(1 for t in study.trials
                                  if t.user_attrs.get("feasible", False))
                print(f"  feasible trials (RMS <= eps): {n_feasible} /"
                      f" {len(study.trials)}")
                if best_feasible[0] < float("inf"):
                    print(f"  best feasible: reg_norm = {best_feasible[0]:.4e},"
                          f" RMS = {best_feasible[1]:.3e}")
                else:
                    print(f"  WARNING: no feasible trial found "
                          f"(no trial reached RMS <= {args.eps_rms:.2e})")

            # Evaluate best params for full output
            best = study.best_params
            deform_params = {}
            if "zoff" in best:
                deform_params["zoff"] = best["zoff"]
            if "bump_amp" in best:
                deform_params["bump"] = (best["bump_amp"], best["bump_cx"],
                                           best["bump_cy"])

        def surface_z_fn(x, y):
            from ngsolve import exp as ngs_exp
            is_cf = "CoefficientFunction" in type(x).__name__
            exp_fn = ngs_exp if is_cf else np.exp
            z_val = 0.0
            if "zoff" in deform_params:
                z_val = z_val + deform_params["zoff"]
            if "bump" in deform_params:
                amp, cx, cy = deform_params["bump"]
                z_val = z_val + amp * exp_fn(
                    -((x - cx) ** 2 + (y - cy) ** 2) / 0.005)
            return z_val

        target_args = {"half": args.target_half, "z": args.target_z,
                       "n": args.n_target}
        result = solve_and_evaluate(
            args.plane_half, args.maxh, args.order, target_args, args.B0,
            regularize=args.regularize, sigma_cf_expr=args.sigma_cf,
            jmax=args.jmax, surface_z_fn=surface_z_fn,
            n_sample=args.n_sample, nlevels=args.nlevels,
            compute_reg_norm=args.minimize_reg)
    else:
        # single inner solve
        target_args = {"half": args.target_half, "z": args.target_z,
                       "n": args.n_target}
        t0 = time.time()
        result = solve_and_evaluate(
            args.plane_half, args.maxh, args.order, target_args, args.B0,
            regularize=args.regularize, sigma_cf_expr=args.sigma_cf,
            jmax=args.jmax, surface_z_fn=None,
            n_sample=args.n_sample, nlevels=args.nlevels,
            compute_reg_norm=args.minimize_reg)
        elapsed = time.time() - t0
        print(f"  solver: {args.regularize} (elapsed {elapsed:.2f} s)")

    print(f"[6] I_w = {result['I_w']:.4g} A,"
          f" target-plane RMS = {result['rms']:.3e}")
    print(f"    Bz mean = {result['Bz_mean']:.3e} T (target {args.B0:.3e} T),"
          f" p2p/mean = {result['p2p_mean']:.3e}")
    print(f"    chain: {result['n_wire']} contours,"
          f" wire length = {result['length']:.3f} m")
    if args.minimize_reg and "reg_norm" in result \
            and np.isfinite(result["reg_norm"]):
        print(f"    regularisation norm psi^T S psi = {result['reg_norm']:.4e}"
              f" (mode: {args.regularize})")


if __name__ == "__main__":
    main()
