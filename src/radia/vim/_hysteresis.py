"""B-input hysteresis stepping for the HDiv-VIM (RT1) charge-Gram demag solve.

Quasi-static hysteresis = MANY solves on ONE fixed geometry with an evolving
per-element material state.  The charge Gram G (and hence the demag operator
N = B^T G B) is chi-free geometry, so the HACApK H-matrix is built ONCE and
every step / nonlinear iteration reuses it -- the per-step cost is the W-CG
solve chain only.  This module wires the lab's B-input hysteresis models
(the C++ radTHysteresisMaterial Play/Energy family) into that loop.

Material protocol (duck-typed and deliberately FUNCTIONAL, so lab-local
research models -- e.g. a numpy vector-stop model -- can plug in without
touching radia):

    state0()            -> opaque committed state for ONE evaluation point
    forward(B, state)   -> H (A/m) from flux density B (T); PURE w.r.t. the
                           committed state (repeated calls with varying B
                           never advance the magnetic history)
    commit(B, state)    -> the new committed state after accepting B as the
                           converged flux density of this quasi-static step
    nu_B0()             -> small-signal |dH/dB| at the virgin state

The B-input constitutive relation is inverted POINTWISE per element: solve

    B = mu0 * (forward(B, state) + M)

which is a contraction with rate ~ mu0*dH/dB << 1 (play/stop branch slopes
are positive and bounded below mu0^-1), giving the material field
H_mat(M) = forward(B(M)) as a full VECTOR (no collinearity assumption).

Nonlinear outer iteration = the HANTILA POLARIZATION method (Hantila 1975;
the classical B-input demag iteration): split M = chi0 * H + M_p with a FIXED
chi0 = 1/nu0, keep the polarization source lagged.  Each outer iteration
solves the SPD system

    ( W(nu0) + N ) m^{k+1} = M_mass h_ext + INT (nu0 M^k - H_mat(M^k)) . v dx

whose LHS is CONSTANT across iterations AND steps (assembled once).  At the
fixed point the true constitutive equation INT H_mat(M).v + N m = M_mass h_ext
holds with H_mat unrestricted in direction -- in particular ANTI-PARALLEL
H/M on recoil branches beyond remanence, which a scalar secant nu = |H|/|M|
cannot represent (that structural failure is exactly the descending-branch
Picard divergence the retired moment path saw).  Contraction: the update map
has rate |1 - (dH/dM)/nu0| per mode, so nu0 must be an UPPER BOUND on the
differential reluctivity dH/dM.  For the lab's Play models with the
Hane-Sugahara (negative irreversible shape) convention the maximum slope is
the small-signal one, so nu0 defaults to the virgin value; override `nu0`
for materials whose branch slopes exceed it.

State discipline (the C++ two-buffer contract, rad_material_def.h): a forward
evaluation plays from the COMMITTED baseline (m_pk_prev) into scratch
(m_pk_current); only CommitState promotes scratch -> baseline.  This adapter
still RESTORES the committed state before EVERY evaluation, (a) to multiplex
ONE C++ handle across all elements and (b) to keep B |-> H referentially
transparent regardless of scratch (pinning / last-B warm-start) semantics.

The CALLER opens `with ngsolve.TaskManager():` (same contract as vim.Solve).
"""

import time

import numpy as np
import scipy.sparse as sp
import ngsolve as ng

from ._solve import _i32, _f64, _h_solve_mass_riesz, _resolve_gram_params
from ._vim import build_charge_gram

MU0 = 4.0e-7 * np.pi


class PlayHysteresisMaterial:
    """rad.MatPlayHysteresis-backed material implementing the duck-typed protocol.

    ONE C++ handle serves every element: the committed state lives Python-side
    (one flat array per element, from MatHysSaveState) and is restored into the
    handle before each evaluation, so the handle is a pure evaluator.
    """

    def __init__(self, K, eta, f_k_tables):
        import radia as rad
        self._rad = rad
        self._h = rad.MatPlayHysteresis(int(K), np.asarray(eta, float), f_k_tables)
        self._nu_rev = float(rad.MatHysGetNuRev(self._h))
        self._virgin = np.asarray(rad.MatHysSaveState(self._h), float).copy()

    @property
    def nu_rev(self):
        return self._nu_rev

    def state0(self):
        return self._virgin.copy()

    def forward(self, B, state):
        rad = self._rad
        rad.MatHysRestoreState(self._h, np.asarray(state, float))
        B = np.asarray(B, float)
        H_irr = np.asarray(rad.MatHysIrreversible(self._h, B), float).ravel()[:3]
        return self._nu_rev * B + H_irr

    def commit(self, B, state):
        rad = self._rad
        rad.MatHysRestoreState(self._h, np.asarray(state, float))
        rad.MatHysIrreversible(self._h, np.asarray(B, float))   # play from committed -> scratch
        rad.MatHysCommitState(self._h)                          # scratch -> new committed
        return np.asarray(rad.MatHysSaveState(self._h), float).copy()

    def nu_B0(self, eps=1e-9):
        H = self.forward(np.array([0.0, 0.0, eps]), self.state0())
        return float(np.linalg.norm(H) / eps)


def _solve_pointwise_B(material, state, M, B0, tol=1e-12, maxit=80):
    """Solve B = mu0*(forward(B, state) + M): fixed point, contraction ~mu0*dH/dB."""
    B = np.asarray(B0, float).copy()
    floor = MU0 * (float(np.linalg.norm(M)) + 1.0)
    for _ in range(maxit):
        Bn = MU0 * (material.forward(B, state) + np.asarray(M, float))
        d = float(np.linalg.norm(Bn - B))
        B = Bn
        if d <= tol * max(float(np.linalg.norm(B)), floor):
            return B
    raise RuntimeError(
        "vim.SolveHysteresis: pointwise B-inversion did not converge in %d fixed-point "
        "iterations (|M|=%.3e A/m).  B -> mu0*(H(B)+M) contracts at ~mu0*dH/dB; "
        "non-convergence means a pathological shape function (mu0*dH/dB >= 1, i.e. "
        "differential mu_r <= 1)." % (maxit, float(np.linalg.norm(M))))


def SolveHysteresis(mesh, h_steps, play=None, material=None, *,
                    nu0=None, gram_eps=None, leaf=32, eta=2.0, far_quad=None,
                    ho_far_factor=None, tol=1e-8, maxit=4000,
                    nl_maxit=200, nl_tol=1e-3, relax=1.0):
    """Quasi-static B-input hysteresis stepping on the RT1 HDiv-VIM charge Gram.

    The charge-Gram H-matrix is built ONCE (chi-free geometry) and reused by
    every step and every nonlinear iteration; each step runs the Hantila
    polarization iteration (constant SPD LHS W(nu0) + N, lagged vector
    polarization source) with the hysteresis material evaluated from the
    per-element COMMITTED state, then commits the converged flux density.

    Parameters
    ----------
    mesh    : NGSolve 3D mesh, pure TET / HEX / WEDGE (same scope as vim.Solve).
    h_steps : (n_steps, 3) applied uniform field H_ext per quasi-static step
              (A/m).  Steps are HISTORY: reversals between consecutive steps
              create the hysteresis branches.
    play    : (K, eta_thresholds_T, f_k_tables) -> builds PlayHysteresisMaterial.
    material: a duck-typed material (state0/forward/commit/nu_B0) -- exactly
              one of play / material.
    nu0     : polarization reluctivity (in H = nu0*M terms).  Must upper-bound
              the material's differential dH/dM for guaranteed contraction.
              Default: derived from material.nu_B0() (correct for Play models
              with the Hane-Sugahara negative-irreversible-shape convention,
              whose maximum branch slope is the small-signal one).
    nl_tol  : relative outer-step tolerance on m (engineering default 1e-3).
    relax   : under-relaxation on the polarization-source update (default 1.0;
              the iteration is a contraction, damping is a safety knob only).

    Returns dict with `steps` = per-step records (M (n_el,3), B, H, M_avg,
    B_avg, H_avg, iters, cg_iters, rel_step, t_step_s) + build info (ndof,
    n_el, n_charge, charge_gram_wall_s, t_setup_s, hmat_stats).  The CALLER
    opens `with ngsolve.TaskManager():`.
    """
    if mesh.dim != 3:
        raise ValueError("vim.SolveHysteresis: 3D meshes only (the 2D planar hysteresis "
                         "coupling is a separate increment)")
    _vtx = {len(el.vertices) for el in mesh.Elements(ng.VOL)}
    if _vtx not in ({4}, {8}, {6}):
        raise ValueError(
            "vim.SolveHysteresis supports a pure-TET (4-vertex), pure-HEX (8-vertex), or "
            "pure-WEDGE (6-vertex) mesh; got vertex counts %s." % sorted(_vtx))
    if (play is None) == (material is None):
        raise ValueError("vim.SolveHysteresis: provide EXACTLY ONE of play=(K, eta, f_k_tables) "
                         "or material=<duck-typed B-input material>")
    if play is not None:
        material = PlayHysteresisMaterial(*play)
    h_steps = np.asarray(h_steps, float)
    if h_steps.ndim != 2 or h_steps.shape[1] != 3 or h_steps.shape[0] < 1:
        raise ValueError("vim.SolveHysteresis: h_steps must be (n_steps, 3) applied H (A/m)")

    # polarization reluctivity nu0 (H = nu0*M terms): from the small-signal slope
    # H = nu_B B (linear), B = mu0 (H + M)  =>  H = [mu0 nu_B / (1 - mu0 nu_B)] M.
    if nu0 is None:
        x = MU0 * material.nu_B0()
        if not (0.0 < x < 1.0):
            raise ValueError("vim.SolveHysteresis: the material's small-signal slope gives "
                             "mu0*dH/dB = %.3e; a soft-iron hysteresis model needs "
                             "0 < mu0*dH/dB < 1 (differential mu_r > 1)" % x)
        nu0 = x / (1.0 - x)
    nu0 = float(nu0)
    if nu0 <= 0.0:
        raise ValueError("vim.SolveHysteresis: nu0 must be positive (got %r)" % nu0)

    # ---- ONE-TIME setup: fes + charge-Gram H-matrix (chi-free -> reused by every step) ----
    t_total = time.perf_counter()
    _gp = _resolve_gram_params(order=1, gram_backend="analytic", linear_solver="auto",
                               uniform_linear=False, gram_eps=gram_eps,
                               near_factor=None, far_quad=far_quad, ho_far_factor=ho_far_factor)
    fes = ng.HDiv(mesh, order=1)
    t_before_gram = time.perf_counter()
    B, H, M_mass = build_charge_gram(fes, eps=_gp["eps"], leafsize=leaf, eta=eta,
                                     far_quad=_gp["far_quad"], ho_far_factor=_gp["ho_far_factor"],
                                     nonlinear=True)
    charge_gram_wall_s = time.perf_counter() - t_before_gram
    Mm = sp.csr_matrix(M_mass)
    Bc = sp.csr_matrix(B)
    n_face = fes.ndof
    n_el = mesh.GetNE(ng.VOL)
    Bptr = _i32(Bc.indptr); Bidx = _i32(Bc.indices); Bdat = _f64(Bc.data)

    uf, vf = fes.TnT()
    l2 = ng.L2(mesh, order=0)
    gfM = ng.GridFunction(fes)
    gfHext = ng.GridFunction(fes)
    gfS = [ng.GridFunction(l2) for _ in range(3)]           # lagged polarization source (per component)
    cfS = ng.CoefficientFunction(tuple(gfS))
    vol_el = np.asarray(ng.Integrate(ng.CoefficientFunction(1.0), mesh, element_wise=True), float)
    w_el = vol_el / float(np.sum(vol_el))

    # CONSTANT SPD LHS: W(nu0) + N -- assembled once, reused by every iteration of every step.
    gfNu0 = ng.GridFunction(l2)
    gfNu0.vec.FV().NumPy()[:] = nu0
    a0 = ng.BilinearForm(fes)
    a0 += gfNu0 * uf * vf * ng.dx
    a0.Assemble()
    _r, _c, _v = a0.mat.COO()
    W0 = sp.coo_matrix((np.array(_v), (np.array(_r), np.array(_c))), shape=(n_face, n_face))

    def _solve_W0(rhs, x0=None):
        res = _h_solve_mass_riesz(H, Bptr, Bidx, Bdat, int(n_face),
                                  W0.row, W0.col, W0.data, 1.0, rhs, tol, int(maxit), x0=x0)
        it = int(res["iters"])
        if it >= int(maxit):
            raise RuntimeError(
                "vim.SolveHysteresis (inner W-CG): did NOT converge in %d iters (n_face=%d); "
                "tighten gram_eps or raise maxit." % (maxit, n_face))
        return np.asarray(res["m"], float), it

    def _M_el(m):
        gfM.vec.FV().NumPy()[:] = m
        return np.vstack([
            np.asarray(ng.Integrate(gfM[i], mesh, element_wise=True), float) / vol_el
            for i in range(3)
        ]).T.copy()

    def _polarization_rhs(rhs_src, s_el):
        for i in range(3):
            gfS[i].vec.FV().NumPy()[:] = s_el[:, i]
        lf = ng.LinearForm(fes)
        lf += ng.InnerProduct(cfS, vf) * ng.dx
        lf.Assemble()
        return rhs_src + lf.vec.FV().NumPy()

    states = [material.state0() for _ in range(n_el)]
    B_cache = np.zeros((n_el, 3))
    s_el = np.zeros((n_el, 3))          # lagged polarization source nu0*M - H_mat(M)
    m = np.zeros(n_face)
    t_setup_s = time.perf_counter() - t_total

    steps_out = []
    for istep in range(h_steps.shape[0]):
        hv = h_steps[istep]
        t_step = time.perf_counter()
        gfHext.Set(ng.CoefficientFunction(tuple(hv)))
        h_ext = gfHext.vec.FV().NumPy().copy()
        rhs_src = np.asarray(Mm @ h_ext).ravel()

        cg_total = 0
        rel = float("inf")
        M_el = None
        H_el = np.zeros((n_el, 3))
        nit = 0
        warm = (istep > 0)
        for it in range(int(nl_maxit)):
            m_new, cg_it = _solve_W0(_polarization_rhs(rhs_src, s_el),
                                     x0=m if (warm or it > 0) else None)
            cg_total += cg_it
            rel = float(np.linalg.norm(m_new - m)) / (float(np.linalg.norm(m_new)) + 1e-30)
            m = m_new
            M_el = _M_el(m)
            # material update: pointwise B-inversion from the COMMITTED states -> full-vector
            # H_mat (no collinearity assumption; recoil anti-parallel H/M is representable).
            s_new = np.empty((n_el, 3))
            for e in range(n_el):
                B0e = B_cache[e] if np.any(B_cache[e]) else MU0 * (M_el[e] + hv)
                Be = _solve_pointwise_B(material, states[e], M_el[e], B0e)
                B_cache[e] = Be
                H_el[e] = material.forward(Be, states[e])
                s_new[e] = nu0 * M_el[e] - H_el[e]
            nit = it + 1
            if rel < nl_tol and it > 0:
                break
            s_el = relax * s_new + (1.0 - relax) * s_el
        else:
            raise RuntimeError(
                "vim.SolveHysteresis: step %d (H_ext=%s) did NOT converge -- rel step %.2e > "
                "nl_tol %.1e after %d polarization iters.  The Hantila iteration contracts "
                "only when nu0 upper-bounds the material's dH/dM -- raise nu0, or reduce the "
                "field-step size." % (istep, np.array2string(hv, precision=3), rel, nl_tol, nit))

        # ---- COMMIT: advance every element's play state at the converged flux density ----
        for e in range(n_el):
            states[e] = material.commit(B_cache[e], states[e])

        steps_out.append(dict(
            h_applied=hv.copy(),
            M=M_el.copy(), B=B_cache.copy(), H=H_el.copy(),
            M_avg=np.asarray((w_el[:, None] * M_el).sum(axis=0), float),
            B_avg=np.asarray((w_el[:, None] * B_cache).sum(axis=0), float),
            H_avg=np.asarray((w_el[:, None] * H_el).sum(axis=0), float),
            iters=int(nit), cg_iters=int(cg_total), rel_step=float(rel),
            t_step_s=time.perf_counter() - t_step,
        ))

    out = dict(
        steps=steps_out,
        ndof=int(n_face), n_el=int(n_el), n_charge=int(Bc.shape[0]),
        nu0=float(nu0),
        t_setup_s=float(t_setup_s),
        charge_gram_wall_s=float(charge_gram_wall_s),
        t_steps_s=float(sum(s["t_step_s"] for s in steps_out)),
        total_wall_s_internal=float(time.perf_counter() - t_total),
    )
    if hasattr(H, "stats"):
        out["hmat_stats"] = dict(H.stats())
    return out
