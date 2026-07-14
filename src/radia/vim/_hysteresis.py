"""B-input hysteresis stepping for the HDiv-VIM (RT1) charge-Gram demag solve.

Quasi-static hysteresis = MANY solves on ONE fixed geometry with an evolving
per-element material state.  The charge Gram G (and hence the demag operator
N = B^T G B) is chi-free geometry, so the HACApK H-matrix is built ONCE and
every step / nonlinear iteration reuses it -- the per-step cost is the W-CG
solve chain only.  This module wires the lab's B-input hysteresis models
(the C++ radTHysteresisMaterial Play/Energy family) into that loop.

Material protocol (duck-typed, FUNCTIONAL, and BATCHED -- one call evaluates
every element, so numpy research models vectorize across elements):

    state0()             -> flat committed state (S,) for ONE evaluation point
    forward(B, states)   -> H (n,3) from flux densities B (n,3) and committed
                            states (n,S); PURE w.r.t. the states (repeated
                            calls with varying B never advance the history)
    commit(B, states)    -> new committed states (n,S) after accepting B as
                            the converged flux densities of this step
    nu_bound()           -> an UPPER BOUND on the material's dH/dB over the
                            operating range (NOT merely the origin slope --
                            the Hantila contraction requirement derives from
                            this bound)

The B-input constitutive relation is inverted POINTWISE per element (batched
across elements): solve

    B = mu0 * (forward(B, states) + M)

which is a contraction with rate ~ mu0*dH/dB << 1 (play/stop branch slopes
are positive and bounded below mu0^-1), giving the material field
H_mat(M) as a full VECTOR field (no collinearity assumption).

Nonlinear outer iteration = the HANTILA POLARIZATION method (Hantila 1975;
the classical B-input demag iteration): split M = chi0 * H + M_p with a FIXED
chi0 = 1/nu0, keep the polarization source lagged.  Each outer iteration
solves the SPD system

    ( W(nu0) + N ) m^{k+1} = M_mass h_ext + INT (nu0 M^k - H_mat(M^k)) . v dx

whose LHS is CONSTANT across iterations AND steps (assembled once).  The
material law is evaluated as a ONE POINT PER ELEMENT closure: M^k and
H_mat(M^k) are the ELEMENT-AVERAGED fields (elementwise constant), so the
fixed point satisfies the constitutive equation in the element-averaged sense

    INT H_mat(M_avg).v + N m = M_mass h_ext - nu0 INT (m - M_avg).v dx

-- intra-element RT1 fluctuations are closed LINEARLY at nu0, a single-point
discretization choice whose closure term vanishes under mesh refinement.
H_mat is unrestricted in direction -- in particular ANTI-PARALLEL H/M on
recoil branches beyond remanence, which a scalar secant nu = |H|/|M|
cannot represent (that structural failure is exactly the descending-branch
Picard divergence the retired moment path saw).  Contraction: the update map
has rate |1 - (dH/dM)/nu0| per mode, so nu0 must be an UPPER BOUND on
dH/dM (the inverse differential susceptibility).  The default derives from
material.nu_bound(), the material's own certified sup of dH/dB; override
`nu0` only when that bound is not trustworthy.

State discipline (the C++ two-buffer contract, rad_material_def.h): a forward
evaluation plays from the COMMITTED baseline (m_pk_prev) into scratch
(m_pk_current); only CommitState promotes scratch -> baseline.  The shipped
Play adapter still RESTORES the committed state before EVERY evaluation, (a)
to multiplex ONE C++ handle across all elements and (b) to keep B |-> H
referentially transparent regardless of scratch (pinning / last-B warm-start)
semantics.

The CALLER opens `with ngsolve.TaskManager():` (same contract as vim.Solve).
"""

import time

import numpy as np
import scipy.sparse as sp
import ngsolve as ng

from . import _solve as _solvemod
from ._solve import (_i32, _f64, _h_solve_mass_riesz, _resolve_gram_params,
                     _clear_cpp_solve_timings, _capture_cpp_solve_timings)
from ._vim import build_charge_gram

MU0 = 4.0e-7 * np.pi


class PlayHysteresisMaterial:
    """rad.MatPlayHysteresis-backed material implementing the duck-typed protocol.

    ONE C++ handle serves every element: the committed states live Python-side
    (one flat array per element, from MatHysSaveState) and are restored into
    the handle row-by-row INSIDE the C++ batch entries (MatHysForwardBatch /
    MatHysCommitBatch), so a batched call costs two Python<->C++ crossings in
    total instead of two per element.  The handle stays a pure evaluator:
    forward never commits; commit advances every row exactly once.
    """

    def __init__(self, K, eta, f_k_tables):
        import radia as rad
        self._rad = rad
        self._eta = np.asarray(eta, float).copy()
        # Valid |B| range = the smallest shape-function table extent: the eta=0 operator sees
        # |p_0| = |B| directly, and every f_k is EXTRAPOLATED beyond its table's largest r.
        self._r_max = float(min(np.max(np.asarray(r, float)) for r, _ in f_k_tables))
        self._h = rad.MatPlayHysteresis(int(K), self._eta, f_k_tables)
        self._nu_rev = float(rad.MatHysGetNuRev(self._h))
        self._virgin = np.asarray(rad.MatHysSaveState(self._h), float).copy()

    @property
    def nu_rev(self):
        return self._nu_rev

    def state0(self):
        return self._virgin.copy()

    def forward(self, B, states):
        return np.asarray(self._rad.MatHysForwardBatch(
            self._h, _f64(B), _f64(states)), float)

    def commit(self, B, states):
        return np.asarray(self._rad.MatHysCommitBatch(
            self._h, _f64(B), _f64(states)), float)

    def nu_bound(self):
        """Certified upper bound on dH/dB: the C++ ComputeNuRev scan of the
        virgin curve's maximum total slope.  Forward decomposes as
        H = nu_rev*B + H_irr with every H_irr slope <= 0 by construction, so
        nu_rev bounds dH/dB on every branch, not just at the origin."""
        return self._nu_rev

    def b_max(self):
        """Largest |B| (Tesla) the model was identified for -- the extent of the
        shape-function tables (the eta=0 operator evaluates f_0 at |p_0|=|B|).
        Beyond it every f_k EXTRAPOLATES, carrying no identified information, so
        SolveHysteresis raises rather than trusting an out-of-domain evaluation."""
        return self._r_max


def _solve_pointwise_B(material, states, M, B0, tol=1e-12, maxit=200):
    """Solve B = mu0*(forward(B, states) + M) for ALL elements at once.

    Each row is an independent contraction with rate ~mu0*dH/dB; the batch
    iterates until EVERY row's update is below tolerance.  Deep saturation
    (differential mu_r -> 1) pushes the rate toward 1 and legitimately needs
    more iterations, hence the generous default budget."""
    B = np.asarray(B0, float).copy()
    M = np.asarray(M, float)
    floor = MU0 * (np.linalg.norm(M, axis=1) + 1.0)
    for _ in range(maxit):
        H = material.forward(B, states)
        Bn = MU0 * (H + M)
        d = np.linalg.norm(Bn - B, axis=1)
        B = Bn
        if np.all(d <= tol * np.maximum(np.linalg.norm(B, axis=1), floor)):
            return B, H
    worst = int(np.argmax(d))
    raise RuntimeError(
        "vim.SolveHysteresis: pointwise B-inversion did not converge in %d fixed-point "
        "iterations (worst element %d, |M|=%.3e A/m).  B -> mu0*(H(B)+M) contracts at "
        "rate ~mu0*dH/dB, so deep saturation (differential mu_r -> 1) legitimately "
        "needs many iterations -- raise maxit for such drives; mu0*dH/dB >= 1 "
        "(differential mu_r <= 1) cannot converge at all."
        % (maxit, worst, float(np.linalg.norm(M[worst]))))


def SolveHysteresis(mesh, h_steps, play=None, material=None, *,
                    nu0=None, gram_eps=None, leaf=32, eta=2.0, far_quad=None,
                    ho_far_factor=None, tol=1e-8, maxit=4000,
                    nl_maxit=200, nl_tol=1e-3):
    """Quasi-static B-input hysteresis stepping on the RT1 HDiv-VIM charge Gram.

    The charge-Gram H-matrix is built ONCE (chi-free geometry) and reused by
    every step and every nonlinear iteration; each step runs the Hantila
    polarization iteration (constant SPD LHS nu0*M_mass + N, lagged vector
    polarization source) with the hysteresis material evaluated from the
    per-element COMMITTED states, then commits the converged flux density.
    The constant system's mass-Riesz PARDISO factor is warmed up during setup
    (booked in t_setup_s), so per-step wall times measure the reuse regime.

    Parameters
    ----------
    mesh    : NGSolve 3D mesh, pure TET / HEX / WEDGE (same scope as vim.Solve).
    h_steps : (n_steps, 3) applied uniform field H_ext per quasi-static step
              (A/m).  Steps are HISTORY: reversals between consecutive steps
              create the hysteresis branches.
    play    : (K, eta_thresholds_T, f_k_tables) -> builds PlayHysteresisMaterial.
    material: a duck-typed BATCHED material (state0/forward/commit/nu_bound,
              see the module docstring) -- exactly one of play / material.
    nu0     : polarization reluctivity (in H = nu0*M terms).  Must upper-bound
              the material's differential dH/dM for guaranteed contraction.
              Default: derived from material.nu_bound(), the material's own
              certified upper bound on dH/dB.
    nl_tol  : outer-step tolerance on m, measured against the running maximum
              of ||m|| over the loop (a uniform ABSOLUTE accuracy across the
              cycle) with a contraction-corrected acceptance (see the module
              docstring).  Engineering default 1e-3.

    Returns dict with `steps` = per-step records (M (n_el,3), B, H, M_avg,
    B_avg, H_avg, iters, cg_iters, rel_step, t_step_s) + build info (ndof,
    n_el, n_charge, charge_gram_wall_s, t_setup_s, hmat_stats) + the summed
    C++ inner-solve timing breakdown (cpp_solve_timings).  The CALLER opens
    `with ngsolve.TaskManager():`.
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

    # polarization reluctivity nu0 (H = nu0*M terms) from the material's certified
    # UPPER BOUND on dH/dB: with H = nu_B B and B = mu0 (H + M),
    # H = [mu0 nu_B / (1 - mu0 nu_B)] M, and nu_B -> dH/dM is monotone increasing,
    # so a sup bound on dH/dB maps to the sup bound on dH/dM that contraction needs.
    if nu0 is None:
        x = MU0 * float(material.nu_bound())
        if not (0.0 < x < 1.0):
            raise ValueError("vim.SolveHysteresis: material.nu_bound() gives "
                             "mu0*sup(dH/dB) = %.3e; the Hantila iteration needs "
                             "0 < mu0*dH/dB < 1 everywhere (differential mu_r > 1)" % x)
        nu0 = x / (1.0 - x)
    nu0 = float(nu0)
    if nu0 <= 0.0:
        raise ValueError("vim.SolveHysteresis: nu0 must be positive (got %r)" % nu0)

    # Optional out-of-range guard: a material identified only up to some |B|_max (e.g. the
    # play model's largest threshold) EXTRAPOLATES beyond it, and the demag-limited coupled
    # solve can be driven there (M runs away, unphysical B).  If the material exposes b_max(),
    # fail loud when a converged step exceeds it rather than trusting an out-of-domain law.
    b_max = None
    if hasattr(material, "b_max"):
        try:
            b_max = float(material.b_max())
        except Exception:
            b_max = None
        if b_max is not None and not (b_max > 0.0):
            b_max = None

    # ---- ONE-TIME setup: fes + charge-Gram H-matrix (chi-free -> reused by every step) ----
    t_total = time.perf_counter()
    _clear_cpp_solve_timings()
    _gp = _resolve_gram_params(gram_eps=gram_eps, far_quad=far_quad, ho_far_factor=ho_far_factor)
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

    uf = fes.TrialFunction()
    l2 = ng.L2(mesh, order=0)
    gfHext = ng.GridFunction(fes)
    vol_el = np.asarray(ng.Integrate(ng.CoefficientFunction(1.0), mesh, element_wise=True), float)
    w_el = vol_el / float(np.sum(vol_el))

    # Mixed element-integral matrix P (built once): row c*n_el + e holds INT_e u_c dx, so both
    # per-iteration material couplings are sparse mat-vecs instead of NGSolve assembly calls --
    # the element average M_el = (P m)/vol and the polarization load P^T s (exact, since the
    # lagged source is constant per element).  bonus_intorder=4 sets the quadrature to order
    # 1+0+4 = 5, the ng.Integrate default this replaces (matches to <=2e-14 on warped hexes);
    # L2(order=0) dof numbering is the element numbering, which the vstack layout relies on.
    assert l2.ndof == n_el
    wl2 = l2.TestFunction()
    _P_blocks = []
    for _c in range(3):
        _blf = ng.BilinearForm(trialspace=fes, testspace=l2)
        _blf += uf[_c] * wl2 * ng.dx(bonus_intorder=4)
        _blf.Assemble()
        _pr, _pc, _pv = _blf.mat.COO()
        _P_blocks.append(sp.csr_matrix((_f64(_pv), (np.asarray(_pr), np.asarray(_pc))),
                                       shape=(l2.ndof, n_face)))
    P = sp.vstack(_P_blocks).tocsr()
    PT = P.T.tocsr()

    # CONSTANT SPD LHS: nu0*M_mass + N.  The mass is uniform, so instead of assembling a
    # separate nu0-weighted BilinearForm we pass the ALREADY-BUILT M_mass with the scalar
    # inv_chi = nu0 (the same system; PCG is invariant under the preconditioner's scale,
    # and factoring M_mass itself shares the persistent factor with any other M_mass user).
    Mm_coo = Mm.tocoo()
    Wrow = _i32(Mm_coo.row); Wcol = _i32(Mm_coo.col); Wdat = _f64(Mm_coo.data)

    def _solve_W0(rhs, x0=None):
        res = _h_solve_mass_riesz(H, Bptr, Bidx, Bdat, int(n_face),
                                  Wrow, Wcol, Wdat, float(nu0), rhs, tol, int(maxit), x0=x0)
        _capture_cpp_solve_timings(res)
        it = int(res["iters"])
        if it >= int(maxit):
            raise RuntimeError(
                "vim.SolveHysteresis (inner W-CG): did NOT converge in %d iters (n_face=%d); "
                "tighten gram_eps or raise maxit." % (maxit, n_face))
        return np.asarray(res["m"], float), it

    # FACTOR WARMUP: the persistent mass-Riesz PARDISO factor of the constant mass is
    # built here (a zero-RHS solve converges in 0 CG iterations), so the one-time
    # analyze+factor lands in t_setup_s instead of the first step's t_step_s.
    _solve_W0(np.zeros(n_face))

    def _M_el(m):
        return (P @ m).reshape(3, n_el).T / vol_el[:, None]

    def _polarization_rhs(rhs_src, s_el):
        return rhs_src + PT @ s_el.T.ravel()

    states = np.tile(material.state0()[None, :], (n_el, 1))
    B_cache = None                      # previous converged B per element (None until the first solve)
    s_el = np.zeros((n_el, 3))          # lagged polarization source nu0*M - H_mat(M)
    m = np.zeros(n_face)
    m_scale = 0.0                       # running max ||m|| -> uniform absolute stop across the cycle
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
        H_el = None
        nit = 0
        d_prev = None
        for it in range(int(nl_maxit)):
            m_new, cg_it = _solve_W0(_polarization_rhs(rhs_src, s_el), x0=m)
            cg_total += cg_it
            d_now = float(np.linalg.norm(m_new - m))
            m = m_new
            m_scale = max(m_scale, float(np.linalg.norm(m_new)))
            rel = d_now / (m_scale + 1e-30)
            M_el = _M_el(m)
            # material update (BATCHED): pointwise B-inversion from the COMMITTED states ->
            # full-vector H_mat (recoil anti-parallel H/M is representable).
            B0 = MU0 * (M_el + hv[None, :]) if B_cache is None else B_cache
            B_cache, H_el = _solve_pointwise_B(material, states, M_el, B0)
            s_el = nu0 * M_el - H_el
            nit = it + 1
            if it > 0 and rel < nl_tol:
                # Contraction-corrected acceptance: the Cauchy increment under-estimates the
                # distance to the fixed point by ~q/(1-q).  Estimate q from successive
                # increments and require the corrected error below nl_tol as well -- never
                # LOOSER than the raw criterion (the correction factor is clamped >= 1), and
                # binding only for slow contractions q > 0.5.
                q = (d_now / d_prev) if (d_prev is not None and d_prev > 0.0) else 0.0
                if q < 1.0 and rel * max(1.0, q / (1.0 - q)) < nl_tol:
                    break
            d_prev = d_now
        else:
            raise RuntimeError(
                "vim.SolveHysteresis: step %d (H_ext=%s) did NOT converge -- rel step %.2e > "
                "nl_tol %.1e after %d polarization iters.  The Hantila iteration contracts "
                "only when nu0 upper-bounds the material's dH/dM -- raise nu0 (is "
                "material.nu_bound() a true upper bound for this drive?), or reduce the "
                "field-step size." % (istep, np.array2string(hv, precision=3), rel, nl_tol, nit))

        # ---- out-of-range guard: the converged flux density must stay within the material's
        # identified range, else the material law was extrapolated (No-Fallbacks: fail loud). ----
        if b_max is not None:
            b_peak = float(np.max(np.linalg.norm(B_cache, axis=1)))
            if b_peak > b_max * (1.0 + 1e-6):
                raise RuntimeError(
                    "vim.SolveHysteresis: step %d (H_ext=%s) drove |B|=%.3f T past the material's "
                    "identified range b_max=%.3f T -- the material law would be EXTRAPOLATED there. "
                    "Reduce the applied field so peak |B| stays <= b_max, or supply a material "
                    "identified to higher B."
                    % (istep, np.array2string(hv, precision=3), b_peak, b_max))

        # ---- COMMIT: advance every element's state at the converged flux density ----
        states = material.commit(B_cache, states)

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
        cpp_solve_timings=dict(_solvemod._LAST_CPP_SOLVE_TIMINGS),
        hmat_stats=dict(H.stats()),
    )
    return out
