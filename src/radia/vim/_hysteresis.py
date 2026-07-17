"""B-input hysteresis stepping for the HDiv-VIM BDM charge-Gram demag solve.

Quasi-static hysteresis = MANY solves on ONE fixed geometry with an evolving
material state.  BDM1 uses one state per element; BDM2 uses physical volume
quadrature points.  The charge Gram G (and hence the demag operator
N = B^T G B) is chi-free geometry, so the HACApK H-matrix is built ONCE and
every step / nonlinear iteration reuses it -- the per-step cost is the W-CG
solve chain only.  This module wires B-input hysteresis models into that
loop.  The production hard-magnet path is EnergyStopMaterial: an isotropic
vector Stop model whose shape tables are monotone, whose state lies in fixed
balls |s_k| <= eta_k, and whose batched trial/commit kernels run in C++.

Material protocol (duck-typed, FUNCTIONAL, and BATCHED -- one call evaluates
every constitutive state point):

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

The B-input constitutive relation is inverted POINTWISE at the constitutive
state points: solve

    B = mu0 * (forward(B, states) + M)

which is a contraction with rate ~ mu0*dH/dB << 1 (play/stop branch slopes
are positive and bounded below mu0^-1), giving the material field
H_mat(M) as a full VECTOR field (no collinearity assumption).

Nonlinear outer iteration = the HANTILA POLARIZATION method (Hantila 1975;
the classical B-input demag iteration): split M = chi0 * H + M_p with a FIXED
chi0 = 1/nu0, keep the polarization source lagged.  Each outer iteration
solves the SPD system

    ( W(nu0) + N ) m^{k+1} = M_mass h_ext + INT (nu0 M^k - H_mat(M^k)) . v dx

whose LHS is CONSTANT across iterations AND steps (assembled once).  BDM1 uses
the historical ONE POINT PER ELEMENT closure.  BDM2 evaluates the material at
volume quadrature points and returns it through the matching weighted transpose
of the physical Piola evaluation map.  The BDM1 fixed point satisfies the
constitutive equation in the element-averaged sense

    INT H_mat(M_avg).v + N m = M_mass h_ext - nu0 INT (m - M_avg).v dx

-- intra-element BDM1 fluctuations are closed LINEARLY at nu0, a single-point
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
from ._solve import (_f64, _h_solve_mass_riesz, _resolve_gram_params,
                     _clear_cpp_solve_timings, _capture_cpp_solve_timings,
                     _configure_cpp_mass)
from ._vim import build_charge_gram
from ._capabilities import validate_hdiv_configuration

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

    permanent_magnet_model = "simplified-play"
    permanent_magnet_level = 3

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


class EnergyStopMaterial:
    """C++ isotropic vector B-input energy Stop material.

    ``eta`` contains the Stop radii in tesla.  ``g_tables`` contains one
    ``(r, g)`` pair per branch, where r is in tesla and g is in A/m.  Every
    table must start at (0, 0), cover its eta, and be non-negative and
    non-decreasing.  These generic tables are the public constitutive input;
    fitted magnet-grade data are intentionally not bundled with Radia.

    ``alpha`` is the positive reversible reluctivity floor in A/(m T).
    ``gamma`` is a scalar or one value per branch in T m/A.  gamma=0 gives
    the rate-independent hard ball projection; gamma>0 uses the convex
    radial proximal update.  The explicit state is restartable: ``forward``
    is pure and ``commit`` is the only history-advancing operation.
    """

    permanent_magnet_model = "b-input-energy-stop"
    permanent_magnet_level = 4

    def __init__(self, eta, g_tables, *, alpha=5.0, gamma=0.0, b_max=None):
        import radia._radia_pybind as _rp

        eta = np.asarray(eta, dtype=float)
        if eta.ndim != 1 or eta.size == 0:
            raise ValueError("vim.EnergyStopMaterial: eta must be a non-empty 1D array")
        tables = list(g_tables)
        if len(tables) != eta.size:
            raise ValueError("vim.EnergyStopMaterial: g_tables must contain one table per eta")

        gamma = np.asarray(gamma, dtype=float)
        if gamma.ndim == 0:
            gamma = np.full(eta.size, float(gamma))
        if gamma.shape != eta.shape:
            raise ValueError("vim.EnergyStopMaterial: gamma must be scalar or have shape eta.shape")

        table_r = []
        table_g = []
        offsets = [0]
        for index, table in enumerate(tables):
            if len(table) != 2:
                raise ValueError("vim.EnergyStopMaterial: table %d must be an (r, g) pair" % index)
            r = np.asarray(table[0], dtype=float)
            g = np.asarray(table[1], dtype=float)
            if r.ndim != 1 or g.ndim != 1 or r.shape != g.shape:
                raise ValueError("vim.EnergyStopMaterial: table %d arrays must be matching 1D arrays" % index)
            table_r.extend(r)
            table_g.extend(g)
            offsets.append(len(table_r))

        self._core = _rp._EnergyStopMaterial(
            np.ascontiguousarray(eta),
            np.ascontiguousarray(table_r, dtype=float),
            np.ascontiguousarray(table_g, dtype=float),
            np.ascontiguousarray(offsets, dtype=np.int32),
            np.ascontiguousarray(gamma), float(alpha),
            0.0 if b_max is None else float(b_max),
        )

    @property
    def alpha(self):
        return float(self._core.alpha)

    @property
    def eta(self):
        return np.asarray(self._core.eta, dtype=float)

    @property
    def gamma(self):
        return np.asarray(self._core.gamma, dtype=float)

    @property
    def state_size(self):
        return int(self._core.state_size)

    def state0(self):
        return np.asarray(self._core.state0(), dtype=float)

    def forward(self, B, states):
        return np.asarray(self._core.forward_batch(_f64(B), _f64(states)), dtype=float)

    def commit(self, B, states):
        return np.asarray(self._core.commit_batch(_f64(B), _f64(states)), dtype=float)

    def stored_energy(self, B, states):
        """Return trial stored-energy density in J/m^3 for every row."""
        return np.asarray(self._core.stored_energy_batch(_f64(B), _f64(states)), dtype=float)

    def nu_bound(self):
        return float(self._core.nu_bound)

    def b_max(self):
        return float(self._core.b_max)


def _normalize_h_steps(h_steps):
    """Return (CoefficientFunction, uniform_vector_or_None, label) records."""
    try:
        numeric = np.asarray(h_steps, dtype=float)
    except (TypeError, ValueError):
        numeric = None
    if numeric is not None and numeric.ndim == 2 and numeric.shape[1] == 3:
        raw_steps = [row for row in numeric]
    else:
        try:
            raw_steps = list(h_steps)
        except TypeError as exc:
            raise ValueError(
                "vim.SolveHysteresis: h_steps must be a sequence of 3-vectors or "
                "3D NGSolve CoefficientFunctions") from exc
    if not raw_steps:
        raise ValueError("vim.SolveHysteresis: h_steps must contain at least one step")

    normalized = []
    for index, value in enumerate(raw_steps):
        try:
            vector = np.asarray(value, dtype=float)
        except (TypeError, ValueError):
            vector = None
        if vector is not None and vector.shape == (3,):
            vector = np.ascontiguousarray(vector)
            normalized.append((ng.CoefficientFunction(tuple(vector)), vector,
                               np.array2string(vector, precision=3)))
            continue
        if getattr(value, "dim", None) != 3:
            raise ValueError(
                "vim.SolveHysteresis: step %d must be a 3-vector or a 3D "
                "NGSolve CoefficientFunction" % index)
        normalized.append((value, None, "CoefficientFunction(step=%d)" % index))
    return normalized


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


class _QuadratureCoupling:
    """NGSolve-native BDM2 state evaluation and weak-load projection."""

    def __init__(self, fes, integration_order):
        from ngsolve.comp import IntegrationRuleSpace

        self.fes = fes
        self.mesh = fes.mesh
        self.integration_order = int(integration_order)
        self.m_field = ng.GridFunction(fes)
        self.scalar_space = IntegrationRuleSpace(
            self.mesh, order=self.integration_order)
        self.source_space = ng.VectorValued(self.scalar_space, dim=3)
        self.measure = ng.dx(intrules=self.scalar_space.GetIntegrationRules())
        self.source_field = ng.GridFunction(self.source_space)
        self.sample_field = ng.GridFunction(self.source_space)
        self.npoints = int(self.scalar_space.ndof)
        if self.npoints == 0:
            raise RuntimeError("vim.SolveHysteresis: empty BDM2 integration-rule space")

        trial = self.source_space.TrialFunction()
        test = fes.TestFunction()
        self.source_to_hdiv = ng.BilinearForm(
            trialspace=self.source_space, testspace=fes)
        self.source_to_hdiv += ng.InnerProduct(trial, test) * self.measure
        self.source_to_hdiv.Assemble()

        l2 = ng.L2(self.mesh, order=0)
        l2_test = l2.TestFunction()
        self.to_element = []
        for component in range(3):
            form = ng.BilinearForm(trialspace=self.source_space, testspace=l2)
            form += trial[component] * l2_test * self.measure
            form.Assemble()
            self.to_element.append(form.mat)

    def evaluate_magnetization(self, coefficients):
        self.m_field.vec.FV().NumPy()[:] = np.asarray(coefficients, dtype=float)
        self.sample_field.Interpolate(self.m_field)
        return self._values(self.sample_field)

    def evaluate_coefficient(self, coefficient):
        self.sample_field.Interpolate(coefficient)
        return self._values(self.sample_field)

    def _values(self, field):
        return np.asarray(field.vec.FV().NumPy(), dtype=float).reshape(3, self.npoints).T.copy()

    def _set_values(self, field, values):
        values = np.asarray(values, dtype=float).reshape(self.npoints, 3)
        field.vec.FV().NumPy()[:] = values.T.ravel()

    def weak_load(self, values):
        self._set_values(self.source_field, values)
        load = self.source_to_hdiv.mat.CreateColVector()
        self.source_to_hdiv.mat.Mult(self.source_field.vec, load)
        return np.asarray(load.FV().NumPy(), dtype=float).copy()

    def element_average(self, values, volumes):
        self._set_values(self.source_field, values)
        averaged = np.empty((len(volumes), 3), dtype=float)
        for component, matrix in enumerate(self.to_element):
            integral = matrix.CreateColVector()
            matrix.Mult(self.source_field.vec, integral)
            averaged[:, component] = np.asarray(
                integral.FV().NumPy(), dtype=float) / volumes
        return averaged


def SolveHysteresis(mesh, h_steps, play=None, material=None, *,
                    nu0=None, gram_eps=None, leaf=32, eta=2.0, far_quad=None,
                    ho_far_factor=None, tol=1e-8, maxit=4000,
                    nl_maxit=200, nl_tol=1e-3,
                    initial_b_path=None, initial_state=None, order=1,
                    state_quadrature_order=None, _operator_cache=None):
    """Quasi-static B-input hysteresis stepping on the BDM HDiv-VIM charge Gram.

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
    h_steps : sequence of applied fields per quasi-static step.  Each entry is
              either a uniform 3-vector in A/m or an NGSolve 3D
              CoefficientFunction.  Steps are HISTORY: reversals between
              consecutive steps create the hysteresis branches.
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
    initial_b_path : optional (n_init, 3) uniform flux-density history in T.
              This initializes the constitutive state before the coupled
              solve, e.g. a manufactured hard-magnet state.  It is not a
              replacement for explicitly modelling the magnetizing fixture.
    initial_state : optional state dict returned by an earlier call.  Restarts
              material history, RT coefficients, B cache, and convergence
              scale without hidden mutable material state.
    order   : HDiv order 1 or 2.  BDM1 retains one constitutive state per element;
              BDM2 stores state at physical volume quadrature points and returns
              its nonlinear source through the matching weighted weak load.
    state_quadrature_order : NGSolve ``IntegrationRuleSpace`` order for BDM2.
              The production default is 3 and is recorded in the restart state.

    Returns dict with `steps` = per-step records (M (n_el,3), B, H, M_avg,
    B_avg, H_avg, iters, cg_iters, rel_step, t_step_s) + build info (ndof,
    n_el, n_charge, charge_gram_wall_s, t_setup_s, hmat_stats) + the summed
    C++ inner-solve timing breakdown (cpp_solve_timings).  The CALLER opens
    `with ngsolve.TaskManager():`.  The returned `state` can be passed back as
    `initial_state` for a continuation run.  The final `gfM` and persistent
    C++ field evaluator satisfy `vim.FieldFromSolution(result, points)`.
    """
    if mesh.dim != 3:
        raise ValueError("vim.SolveHysteresis: 3D meshes only (the 2D planar hysteresis "
                         "coupling is a separate increment)")
    _vtx = {len(el.vertices) for el in mesh.Elements(ng.VOL)}
    if _vtx not in ({4}, {8}, {6}):
        raise ValueError(
            "vim.SolveHysteresis supports a pure-TET (4-vertex), pure-HEX (8-vertex), or "
            "pure-WEDGE (6-vertex) mesh; got vertex counts %s." % sorted(_vtx))
    order = int(order)
    validate_hdiv_configuration(mesh.dim, _vtx, order, mesh.GetCurveOrder())
    if state_quadrature_order is not None and order != 2:
        raise ValueError(
            "vim.SolveHysteresis: state_quadrature_order is a BDM2-only setting")
    if (play is None) == (material is None):
        raise ValueError("vim.SolveHysteresis: provide EXACTLY ONE of play=(K, eta, f_k_tables) "
                         "or material=<duck-typed B-input material>")
    if play is not None:
        material = PlayHysteresisMaterial(*play)
    h_step_records = _normalize_h_steps(h_steps)
    if initial_b_path is not None and initial_state is not None:
        raise ValueError("vim.SolveHysteresis: initial_b_path and initial_state are mutually exclusive")

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
    curve_order = int(mesh.GetCurveOrder())
    cache_contract = dict(
        mesh=mesh, order=order, curve_order=curve_order,
        vertex_counts=frozenset(_vtx), nu0=float(nu0),
        gram_eps=float(_gp["eps"]), leaf=int(leaf), eta=float(eta),
        far_quad=int(_gp["far_quad"]),
        ho_far_factor=float(_gp["ho_far_factor"]),
        operator_kind="hysteresis",
    )
    operator_reused = _operator_cache is not None
    if operator_reused:
        if not isinstance(_operator_cache, dict):
            raise TypeError(
                "vim.SolveHysteresis: internal operator cache must be owned by vim.HDivSolver")
        for key, value in cache_contract.items():
            cached = _operator_cache.get(key)
            matches = (cached is value) if key == "mesh" else (cached == value)
            if not matches:
                raise ValueError(
                    "vim.SolveHysteresis: prepared operator mismatch for %s "
                    "(cached=%r, requested=%r)" % (key, cached, value))
        fes = _operator_cache["fes"]
        H = _operator_cache["charge_gram"]
        n_charge = int(_operator_cache["n_charge"])
        charge_gram_wall_s = 0.0
    else:
        fes = ng.HDiv(mesh, order=order)
        t_before_gram = time.perf_counter()
        B, H, M_mass = build_charge_gram(
            fes, eps=_gp["eps"], leafsize=leaf, eta=eta,
            far_quad=_gp["far_quad"], ho_far_factor=_gp["ho_far_factor"],
            nonlinear=True)
        charge_gram_wall_s = time.perf_counter() - t_before_gram
        n_charge = int(B.shape[0])
        Mm = sp.csr_matrix(M_mass)
    n_face = fes.ndof
    n_el = mesh.GetNE(ng.VOL)

    uf = fes.TrialFunction()
    vf = fes.TestFunction()
    vol_el = np.asarray(ng.Integrate(ng.CoefficientFunction(1.0), mesh, element_wise=True), float)
    w_el = vol_el / float(np.sum(vol_el))
    if order == 1:
        l2 = ng.L2(mesh, order=0)
        assert l2.ndof == n_el
        wl2 = l2.TestFunction()
        blocks = []
        for component in range(3):
            form = ng.BilinearForm(trialspace=fes, testspace=l2)
            form += uf[component] * wl2 * ng.dx(bonus_intorder=4)
            form.Assemble()
            row, col, value = form.mat.COO()
            blocks.append(sp.csr_matrix(
                (_f64(value), (np.asarray(row), np.asarray(col))),
                shape=(l2.ndof, n_face)))
        P = sp.vstack(blocks).tocsr()
        PT = P.T.tocsr()
        n_state_points = n_el
        state_layout = "element-average"
        quadrature_order = None
        state_points = None

        def _M_state(m_coefficients):
            return (P @ m_coefficients).reshape(3, n_el).T / vol_el[:, None]

        def _state_to_element(values):
            return np.asarray(values, dtype=float)

        def _polarization_load(values):
            return PT @ np.asarray(values, dtype=float).T.ravel()
    else:
        quadrature_order = (3 if state_quadrature_order is None
                            else int(state_quadrature_order))
        if quadrature_order < 3:
            raise ValueError(
                "vim.SolveHysteresis: BDM2 state_quadrature_order must be >= 3")
        coupling = _QuadratureCoupling(fes, quadrature_order)
        n_state_points = coupling.npoints
        constant_average = coupling.element_average(
            np.ones((n_state_points, 3), dtype=float), vol_el)
        if not np.allclose(constant_average, 1.0, rtol=1.0e-12, atol=1.0e-12):
            raise RuntimeError(
                "vim.SolveHysteresis: BDM2 integration-rule measure does not match the mesh volume")
        state_layout = "quadrature"

        def _M_state(m_coefficients):
            return coupling.evaluate_magnetization(m_coefficients)

        def _state_to_element(values):
            return coupling.element_average(values, vol_el)

        def _polarization_load(values):
            return coupling.weak_load(values)

    # CONSTANT SPD LHS: nu0*M_mass + N.  The mass is uniform, so instead of assembling a
    # separate nu0-weighted BilinearForm we pass the ALREADY-BUILT M_mass with the scalar
    # inv_chi = nu0 (the same system; PCG is invariant under the preconditioner's scale,
    # and factoring M_mass itself shares the persistent factor with any other M_mass user).
    if not operator_reused:
        _configure_cpp_mass(H, Mm, int(n_face))
        _operator_cache = dict(
            cache_contract, fes=fes, charge_gram=H,
            n_face=int(n_face), n_el=int(n_el), n_charge=int(n_charge))

    def _solve_W0(rhs, x0=None):
        res = _h_solve_mass_riesz(
            H, None, int(n_face), float(nu0), rhs, tol, int(maxit), x0=x0)
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
    if not operator_reused:
        _solve_W0(np.zeros(n_face))

    def _M_el(m):
        return _state_to_element(_M_state(m))

    def _polarization_rhs(rhs_src, s_state):
        return rhs_src + _polarization_load(s_state)

    state0 = np.asarray(material.state0(), dtype=float)
    if state0.ndim != 1 or state0.size == 0 or not np.all(np.isfinite(state0)):
        raise ValueError("vim.SolveHysteresis: material.state0() must be a finite non-empty 1D array")
    states = np.tile(state0[None, :], (n_state_points, 1))
    B_cache = None                      # previous converged B at each material state point
    s_state = np.zeros((n_state_points, 3))
    m = np.zeros(n_face)
    m_scale = 0.0                       # running max ||m|| -> uniform absolute stop across the cycle

    if initial_b_path is not None:
        b_path = np.asarray(initial_b_path, dtype=float)
        if b_path.ndim != 2 or b_path.shape[1] != 3 or b_path.shape[0] < 1 or \
                not np.all(np.isfinite(b_path)):
            raise ValueError("vim.SolveHysteresis: initial_b_path must be a finite (n_init,3) array in T")
        for b_init in b_path:
            B_trial = np.repeat(b_init[None, :], n_state_points, axis=0)
            states = np.asarray(material.commit(B_trial, states), dtype=float)
        B_cache = np.repeat(b_path[-1][None, :], n_state_points, axis=0)
        s_state = -np.asarray(material.forward(B_cache, states), dtype=float)
    elif initial_state is not None:
        if not isinstance(initial_state, dict):
            raise ValueError("vim.SolveHysteresis: initial_state must be the returned state dict")
        try:
            states = np.asarray(initial_state["material_states"], dtype=float).copy()
            m = np.asarray(initial_state["m_coefficients"], dtype=float).copy()
            B_cache = np.asarray(initial_state["B"], dtype=float).copy()
        except KeyError as exc:
            raise ValueError("vim.SolveHysteresis: initial_state is missing %s" % exc) from exc
        expected_layout = initial_state.get("state_layout", "element-average")
        if expected_layout != state_layout:
            raise ValueError(
                "vim.SolveHysteresis: initial state layout %r does not match %r"
                % (expected_layout, state_layout))
        if states.shape != (n_state_points, state0.size):
            raise ValueError("vim.SolveHysteresis: initial material_states shape %s, expected %s" %
                             (states.shape, (n_state_points, state0.size)))
        if m.shape != (n_face,):
            raise ValueError("vim.SolveHysteresis: initial m_coefficients shape %s, expected %s" %
                             (m.shape, (n_face,)))
        if B_cache.shape != (n_state_points, 3):
            raise ValueError("vim.SolveHysteresis: initial B shape %s, expected %s" %
                             (B_cache.shape, (n_state_points, 3)))
        if not (np.all(np.isfinite(states)) and np.all(np.isfinite(m)) and
                np.all(np.isfinite(B_cache))):
            raise ValueError("vim.SolveHysteresis: initial_state arrays must be finite")
        m_scale = max(float(initial_state.get("m_scale", 0.0)), float(np.linalg.norm(m)))
        M_initial = _M_state(m)
        H_initial = np.asarray(material.forward(B_cache, states), dtype=float)
        s_state = nu0 * M_initial - H_initial
    t_setup_s = time.perf_counter() - t_total

    steps_out = []
    for istep, (h_cf, h_uniform, h_label) in enumerate(h_step_records):
        t_step = time.perf_counter()
        source = ng.LinearForm(fes)
        source += ng.InnerProduct(h_cf, vf) * ng.dx
        source.Assemble()
        rhs_src = np.asarray(source.vec.FV().NumPy(), dtype=float).copy()
        if h_uniform is not None:
            Hext_state = np.repeat(h_uniform[None, :], n_state_points, axis=0)
        else:
            if state_layout == "quadrature":
                Hext_state = coupling.evaluate_coefficient(h_cf)
            else:
                Hext_state = np.column_stack([
                    np.asarray(ng.Integrate(h_cf[c], mesh, element_wise=True), dtype=float)
                    for c in range(3)
                ]) / vol_el[:, None]
        Hext_el = _state_to_element(Hext_state)
        Hext_avg = np.asarray((w_el[:, None] * Hext_el).sum(axis=0), dtype=float)

        cg_total = 0
        rel = float("inf")
        M_el = None
        H_el = None
        nit = 0
        d_prev = None
        for it in range(int(nl_maxit)):
            m_new, cg_it = _solve_W0(_polarization_rhs(rhs_src, s_state), x0=m)
            cg_total += cg_it
            d_now = float(np.linalg.norm(m_new - m))
            m = m_new
            m_scale = max(m_scale, float(np.linalg.norm(m_new)))
            rel = d_now / (m_scale + 1e-30)
            M_state = _M_state(m)
            M_el = _state_to_element(M_state)
            # material update (BATCHED): pointwise B-inversion from the COMMITTED states ->
            # full-vector H_mat (recoil anti-parallel H/M is representable).
            B0 = MU0 * (M_state + Hext_state) if B_cache is None else B_cache
            B_cache, H_state = _solve_pointwise_B(material, states, M_state, B0)
            s_state = nu0 * M_state - H_state
            B_el = _state_to_element(B_cache)
            H_el = _state_to_element(H_state)
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
                "field-step size." % (istep, h_label, rel, nl_tol, nit))

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
                    % (istep, h_label, b_peak, b_max))

        # ---- COMMIT: advance every element's state at the converged flux density ----
        states = material.commit(B_cache, states)

        steps_out.append(dict(
            h_applied=Hext_avg.copy(), h_applied_avg=Hext_avg.copy(),
            h_applied_is_uniform=bool(h_uniform is not None),
            M=M_el.copy(), B=B_el.copy(), H=H_el.copy(),
            M_avg=np.asarray((w_el[:, None] * M_el).sum(axis=0), float),
            B_avg=np.asarray((w_el[:, None] * B_el).sum(axis=0), float),
            H_avg=np.asarray((w_el[:, None] * H_el).sum(axis=0), float),
            iters=int(nit), cg_iters=int(cg_total), rel_step=float(rel),
            t_step_s=time.perf_counter() - t_step,
        ))

    gfM = ng.GridFunction(fes)
    gfM.vec.FV().NumPy()[:] = np.asarray(m, dtype=float)
    out = _solvemod._OperatorBackedResult(
        steps=steps_out,
        ndof=int(n_face), n_el=int(n_el), n_charge=int(n_charge),
        gfM=gfM, order=order, curve_order=(2 if curve_order >= 2 else None),
        state_layout=state_layout, state_points=int(n_state_points),
        state_quadrature_order=quadrature_order,
        operator_reused=bool(operator_reused),
        nu0=float(nu0),
        t_setup_s=float(t_setup_s),
        charge_gram_wall_s=float(charge_gram_wall_s),
        t_steps_s=float(sum(s["t_step_s"] for s in steps_out)),
        total_wall_s_internal=float(time.perf_counter() - t_total),
        cpp_solve_timings=dict(_solvemod._LAST_CPP_SOLVE_TIMINGS),
        hmat_stats=dict(H.stats()),
        permanent_magnet_model=getattr(material, "permanent_magnet_model", "custom-b-input"),
        permanent_magnet_level=getattr(material, "permanent_magnet_level", None),
        state=dict(
            material_states=np.asarray(states, dtype=float).copy(),
            m_coefficients=np.asarray(m, dtype=float).copy(),
            B=np.asarray(B_cache, dtype=float).copy(),
            m_scale=float(m_scale),
            state_layout=state_layout,
            state_quadrature_order=quadrature_order,
        ),
        operator_cache=_operator_cache,
    )
    out["_charge_gram"] = H
    out["_m_coefficients"] = np.ascontiguousarray(m, dtype=np.float64)
    from ._field_batch import _materialize_field_evaluator
    _materialize_field_evaluator(out)
    out["total_wall_s_internal"] = float(time.perf_counter() - t_total)
    return out
