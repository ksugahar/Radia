"""
Applications: motor optimization, magnet design, induction heater
shape, accelerator pole face optimization.
"""

MOTOR_OPTIMIZATION = """
# IPM motor shape optimization (Gangl 2015 case study)

## Problem
Interior Permanent Magnet (IPM) brushless motor:
- Stator (outer) with coil slots
- Rotor (inner) with permanent magnets in iron cavities
- Thin air gap between them
- 2D cross-section model (typical for rotating-machine analysis)

Goal: smooth out the **cogging torque** by tuning the geometry of
a specific iron region of the rotor.

The cogging torque comes from interaction between PM stray field
and the variable reluctance as the rotor rotates.  Reducing it
improves motor smoothness and reduces audible noise.

## Mathematical formulation

Cost functional:
    J(Ω, u) = ∫_Γ₀ |B_r(u) - B_d|² ds

where Γ₀ is a circle inside the air gap, B_r is the radial
component of B = curl A (= ∇u × ẑ in 2D), and B_d is a TARGET
sinusoidal profile that gives zero cogging.

State equation:
    -div(ν(|∇u|²) ∇u) = f^PM    in D
    u = 0                        on ∂D

with f^PM the magnetization source term from the permanent magnets,
and ν the reluctivity (highly nonlinear in iron, constant in air).

Design domain Ω ⊂ Ω^ref ⊂ Ω_f^ref:  a small piece of the rotor
iron core (the user can only modify SOME of it; the rest is fixed
by manufacturing constraints).

## Constraint set

    O = {Ω ⊂ Ω^ref open, Lipschitz with bounded Lipschitz constant L_O}

Bounded Lipschitz constant prevents oscillation / chattering during
optimization.

## Existence of optimal shape

Gangl 2015 §3 proves: under Assumption 1 (β monotone + Lipschitz),
the minimization problem has a solution Ω* ∈ O.

The proof uses Hausdorff convergence of characteristic functions
and strong-H¹ convergence of u_n → u* (the latter requires a small
ε-regularity gain over H¹).

## Numerical results (Gangl 2015 §6)

Starting from an initial rotor design:
- ~100 iterations of gradient-based shape opt
- Cost J reduced by ~10× (1 → 0.1)
- Resulting Ω: smooth boundary with a characteristic "bump" that
  cancels the cogging Fourier modes

## Practical NGSolve workflow

```python
from ngsolve import *
mesh = generate_motor_mesh(...)

# Nonlinear reluctivity (B-H curve)
def nu_fn(absB):
    # Frohlich-type interpolation of B-H data
    return mur_air if absB < B_threshold else nu_curve(absB)
nu = CoefficientFunction(...)

# State equation
fes = H1(mesh, order=2, dirichlet="outer")
u = fes.TrialFunction(); v = fes.TestFunction()
a = BilinearForm(fes)
a += nu(absB) * grad(u) * grad(v) * dx
# (Newton's method for the nonlinear solve)

# Adjoint state
# (linear in p, after fixing u from the nonlinear state equation)

# Shape gradient (volume form)
# (loop over elements, accumulate ∇V contributions)

# Riesz-lift shape gradient to a smooth deformation field
W = VectorH1(mesh, order=2)
# ... solve auxiliary problem

# Mesh deformation
mesh.SetDeformation(W * step_size)
```

## Connection to Radia practice

Radia is used for:
- Accelerator magnets (where the cost is field uniformity or
  multipole content along beam axis — analogous to motor cost)
- Insertion devices (wigglers/undulators where the cost is
  spectral content of trajectory)
- Iron pole face design (cost = field flatness over a region)

For all these, the Gangl-Sturm framework applies directly.  Radia
itself doesn't run the optimization — that happens in NGSolve or
similar — but Radia gives the field-evaluation backbone for
validation of the optimized geometry.

## Connection to Radia + Mathematica + radia-mcp

Workflow for a PhD student doing magnet design:
1. Theory layer (radia_mcp.topology_optimization) — read this
   knowledge module
2. Symbolic layer (radia_mcp.mathematica) — derive shape-derivative
   for the specific cost functional in closed form
3. Discretization layer (radia_mcp.radia_ngsolve) — implement in
   NGSolve, run the optimization
4. Validation layer (Radia C++) — verify the optimized field
   against Radia's exact analytical formulas

This is the 4-layer stack from the killer demo in the README.
"""


FIELD_SYNTHESIS = r"""
# Field synthesis — the ANALYTIC / LINEAR-INVERSE branch of EM design opt

The rest of this server teaches the DENSITY / GRADIENT-SHAPE branch of EM
optimization (SIMP, level set, ON/OFF, Gangl-Sturm shape derivative,
Sokolowski topological derivative): iterate a GEOMETRY against a nonlinear
PDE + adjoint.  This topic covers the COMPLEMENTARY branch — design the
SOURCE (magnetization M, or surface current K) so that a TARGET field is
reproduced, exploiting the LINEARITY of magnetostatics in the source.  No
mesh-iteration, no adjoint: in the canonical cases the source is obtained in
CLOSED FORM (analytic inverse-source) or by ONE least-norm solve.

Two members of this branch:

  (1) COIL field synthesis (surface current K = n_hat x grad psi).
      target B over a DSV  ->  A psi = B  ((ACA+)+TSVD least-norm)
      -> iso-contours of psi = the wires.  This is the stream-function /
      target-field / TSVD-eigenmode lineage (Turner; Tomasi; M. Abe).
      ==> ALREADY SHIPPED + golden-locked in the `streamfunction` server.
          Use it; do NOT reimplement.  (streamfunction("overview"),
          ("harmonics"), ("fusion"); the Abe edge-equipotential BC =
          M. Abe's current-potential MRI-coil method.)

  (2) MAGNET field synthesis (permanent-magnet magnetization M).
      "arbitrary field generation by permanent magnets" — given a target
      field, find the radial-magnetization distribution M_r(r,theta).  This
      is the analytic inverse-source method documented below.


## Analytic PM-multipole inverse (2D polar) — VERIFIED

For a permanent magnet the bound-current source of the magnetostatic
vector potential A = A_z(r,theta) z_hat is curl M:

    -Laplacian(A_z) = mu0 (curl M)_z ,   (curl M)_z = -(1/r) dM_r/dtheta

Expand a purely RADIAL magnetization in azimuthal multipoles.  For a single
mode n with an r-INDEPENDENT amplitude (the lab "Hidaka form"):

    M_r        = (M_rn / n) sin(n (theta + theta_0))
    (curl M)_z = -(M_rn / r) cos(n (theta + theta_0))

The PARTICULAR solution is

    A_z,n = mu0 * M_rn * r * cos(n (theta + theta_0)) / (n^2 - 1)        (*)

i.e. a 1/(n^2 - 1) modal coefficient and a LINEAR-in-r potential, which
gives an r-INDEPENDENT (spatially uniform) 2n-pole field component inside
the magnet region:  B_r ~ sin(n(theta+theta_0)), B_theta ~ cos(...).

Because magnetostatics is LINEAR in M, a general radial magnetization is the
SUPERPOSITION of (*) over modes n — so an arbitrary target interior field
expanded in multipoles is inverted mode-by-mode in closed form.  The same
holds for a radially-VARYING magnetization M_r ~ r^{+/-n} (the lab
"Sugahara general form"); its particular solution carries r^{n+1}/(2n+1)
and r^{-n+1}/(2n-1) coefficients.

VERIFICATION (symbolic, golden-locked by
tests/mcp_server/test_topology_field_synthesis.py):
  * Hidaka single-mode:  Laplacian(A_z) - (curl M)_z = 0  exactly.
  * Sugahara general r^{+/-n}:  Laplacian(A_z) - (curl M)_z = 0  exactly.
(Convention: the lab scripts set mu0 = 1 and write the identity as
lapAz = rotMr; the test reproduces that identity symbolically.)

### The n = 1 degeneracy (do not miss this)
The coefficient 1/(n^2 - 1) DIVERGES at n = 1.  n = 1 radial magnetization
M_r ~ sin(theta + theta_0) IS a UNIFORMLY magnetized cylinder — a resonant
/ degenerate forcing whose particular solution is NOT of the r*cos form.
It is handled separately by the classical result: a uniformly magnetized 2D
cylinder has a UNIFORM interior field B_in = mu0 M / 2 (transverse demag
factor 1/2).  So the dipole (n=1) term is the special case, exactly as the
1/(n^2-1) pole signals.  (Cf. the ellipsoid/cylinder demag factors in
`analytical_formulas("ellipsoid")`.)


## Worked objective — air-gap harmonic content (SPM / Halbach ring)

A surface-permanent-magnet (SPM) ring with a chosen segmentation /
magnetization pattern is evaluated by the FOURIER HARMONIC CONTENT of A_z
on an air-gap circle r = R:

    a_0 = (1/2pi) integral A_z dtheta ,
    a_k = (1/pi) integral A_z cos(k theta) dtheta ,
    b_k = (1/pi) integral A_z sin(k theta) dtheta .

The design objective is to shape the magnet arrangement so the spectrum {a_k,
b_k} matches a target (e.g. a pure fundamental for low cogging / clean
back-EMF, or a prescribed gradient/shim profile).  The forward A_z(R,theta)
comes from Radia (rad.Fld) or the analytic multipole sum above; the harmonic
decomposition is a 1-D quadrature.  This pairs the analytic inverse (pick the
modes you want) with a Radia forward check (confirm the realised spectrum).
Cross-validation against an independent magnetostatic solver is kept as an
internal regression reference.


## Outer loop when the inverse is NOT closed-form
When the geometry is constrained (manufacturable segment angles, discrete
magnetization directions, soft-iron return paths) the linear inverse becomes
a bounded nonlinear fit — drive it with a DERIVATIVE-FREE optimizer:
Nelder-Mead simplex / swarm (see the `evolutionary` server and official
`optuna/optuna-mcp`).
The analytic multipole map makes a cheap, exactly-differentiable surrogate
for the inner objective.


## Honest scope
  * (*) is the INTERIOR particular solution for an r-independent single-mode
    radial M; a real finite magnet adds homogeneous r^{+/-n} terms fixed by
    the magnet inner/outer radii + soft-iron BCs (the Sugahara general form
    supplies those radial pieces; match coefficients to the boundary).
  * n=1 is degenerate (uniform-M cylinder, B=mu0 M/2) — never apply (*) there.
  * The COIL branch (stream function / TSVD / target field) is the shipped,
    tested `streamfunction` server; this topic only POINTS to it and adds the
    MAGNET (magnetization) analytic inverse + the SPM harmonic objective.

## References (field-synthesis lineage)
  * Turner, J. Phys. D 19, L147 (1986) — target-field method (gradient coils).
  * Tomasi, Magn. Reson. Med. (2001) — stream-function gradient-coil design.
  * M. Abe et al. — current-potential / TSVD-eigenmode MRI magnet & shim
    design (node current potentials + triangular FE; the edge-equipotential
    BC in the streamfunction server is named after this work).
  * Permanent-magnet "arbitrary field generation" + magnetic-body position
    control (Sugahara Lab analytic multipole note; the (*) derivation above).
See also the density/gradient-shape PM-design entries in
`topology_opt_applications` and the bibliography index.
"""


OUTER_LOOP_OPTIMIZERS = r"""
# Outer-loop optimizers for EM field / shape design (derivative-free)

When the inverse problem is NOT closed-form -- manufacturable bounds on
segment angles, discrete magnetization directions, soft-iron return paths, a
black-box FE/Radia forward solve in the loop -- the design is a bounded,
possibly non-smooth, gradient-free optimization.  The Sugahara Lab optimizer
toolbox (W:\...\MATLAB\30_Optimization) is the practical reference for this
outer loop; the canonical members:

  * NELDER-MEAD SIMPLEX (direct search, Nelder-Mead 1965) -- the workhorse
    derivative-free LOCAL optimizer.  A simplex of n+1 points crawls downhill
    by reflect / expand / contract / shrink.  No gradient, tolerates a
    NON-SMOOTH objective (the lab's CF = sum|A x|, an L1 cost, has no
    derivative at kinks -- Nelder-Mead handles it where gradient methods
    stall).  Local, so multistart / a global stage for multimodal fields.
  * BOUND CONSTRAINTS BY TRANSFORMATION (`fminsearchbnd`, D'Errico) -- the
    lab's actual tool: keep the unconstrained simplex but map the search
    variable into the feasible box so NO evaluation ever leaves it:
        both bounds [LB,UB]:  x = LB + (UB-LB)*(sin(t)+1)/2     (sin transform)
        lower only  [LB, inf): x = LB + t^2                     (quadratic)
        upper only  (-inf,UB]: x = UB - t^2
        equal bounds:          variable fixed, problem reduced in size
    The optimizer sees the unconstrained t; the objective sees the feasible x.
  * POPULATION / GLOBAL members for multimodal or discrete design -- PSO
    (Kennedy-Eberhart 1995; lab EQC_PSO / PSO dirs), GA (lab GA dirs), CMA-ES.
    See the `evolutionary` server (ga_de / pso / cma_es / immune_nsga) and
    official `optuna/optuna-mcp` (TPE / CMA-ES sampler, pruning) -- THIS topic is the
    direct-search LOCAL sibling those servers do not cover.

VERIFICATION (golden-locked by tests/mcp_server/test_topology_field_synthesis.py):
  * Nelder-Mead on Rosenbrock (1-x)^2 + 105(y-x^2)^2 from (3,3) -> (1,1), f=0.
  * fminsearchbnd transform with LB=(2,2), no UB -> the constrained minimum
    (2,4), f=1 (the boundary optimum; matches the fminsearchbnd doc example) --
    the transform is reproduced exactly.

WHERE IT PLUGS IN: the analytic PM-multipole / stream-function inverse
(topic `field_synthesis`) gives a CHEAP, exactly-evaluable surrogate for the
inner objective; this outer loop then searches the manufacturable parameters
(segment count/angles, former size, discrete easy-axis) against the Radia/FE
forward field.  Local direct search for fine tuning; population/global
(`evolutionary`, official `optuna/optuna-mcp`) when the field landscape is multimodal.
"""


LINEAR_INVERSE = r"""
TOPOLOGY / FIELD-SYNTHESIS APPLICATIONS -- REGULARIZED LINEAR INVERSION (TSVD / Tikhonov)
=========================================================================================
The field-synthesis / source-design inverse problem (topic `field_synthesis`) is LINEAR:
a forward matrix A maps a source distribution x (element magnetizations, coil
stream-function weights, shim moments) to the target field samples b = A x.  Build A by
columns -- the field at every observation point produced by each unit source mode (e.g.
each magnet element with Mx/My/Mz, weighted by its volume) -- then SOLVE for x given a
desired b.  A is intrinsically ILL-CONDITIONED (the field is smooth; high source modes
radiate weakly, so their singular values collapse), and the naive least-squares solution
amplifies measurement / discretization noise without bound.

The cure is a SPECTRAL FILTER of the SVD A = U diag(s) V^T:

    minimum-norm LS  x = sum_n (u_n^T b / s_n) v_n                 (= pinv(A) b)
    TSVD (rank k)    x = sum_{n<=k} (u_n^T b / s_n) v_n           (hard cut at mode k)
    Tikhonov (lam)   x = sum_n phi_n (u_n^T b / s_n) v_n,  phi_n = s_n^2/(s_n^2+lam^2)

TSVD and Tikhonov differ ONLY in the filter factors phi_n: a 1/0 step at mode k (TSVD)
vs a smooth roll-off s^2/(s^2+lam^2) (Tikhonov).  Both trade the residual ||A x - b||
(fit) against the solution norm ||x|| (noise amplification) -- the L-CURVE, whose corner
picks the regularization strength.  This is the standard planar gradient / shim-coil
design method and the linear core under the `field_synthesis` PM-multipole / stream-
function inverse: the analytic multipole identity gives the columns of A; this solver
inverts it stably.

IMPLEMENTATION (radia_mcp.topology_optimization.linear_inverse):
  tsvd_solve(A, b, k)        TSVD solution keeping k modes (k=rank -> pinv).
  tikhonov_solve(A, b, lam)  Tikhonov solution, argmin ||A x - b||^2 + lam^2 ||x||^2.
  filter_factors(s, lam)     phi_n = s_n^2/(s_n^2+lam^2).
  lcurve(A, b, lams)         (residual_norm, solution_norm) trade-off points.

VERIFICATION (golden-locked by tests/mcp_server/test_topology_linear_inverse.py, on a
controlled 8x5 rank-3 system with singular values [5, 1, 0.05]):
  * TSVD at k=rank == minimum-norm least-squares pinv(A) b to 1.2e-15.
  * Tikhonov lam->0 -> pinv (2.8e-8 at lam=1e-7).
  * filter factors phi = s^2/(s^2+lam^2) (the 0.05 mode damped to 0.059 at lam=0.2).
  * TSVD residual non-increasing and solution norm non-decreasing in k; Tikhonov L-curve
    residual increasing and solution norm decreasing as lam grows.

WHERE IT PLUGS IN: `field_synthesis` supplies the (analytic) forward map A; this topic
inverts it stably; `outer_loop` then tunes the manufacturable parameters around it.
"""


KRYLOV_SOLVERS = r"""
## Krylov inner solves -- Linear Conjugate Gradient (the `pcg` member)

`krylov.linear_conjugate_gradient(operator, b, x0=None, ...)` solves SPD systems
`A x = b` using only matrix-vector products.  It is the quadratic-minimization
workhorse behind PDE-constrained optimization, Newton/LM inner systems, and
regularized inverse problems when `A` is too large to factor:

    alpha = (r^T z)/(p^T A p),   x <- x + alpha p,   r <- r - alpha A p,
    beta  = (r_new^T z_new)/(r_old^T z_old),         p <- z_new + beta p.

With no preconditioner, `z=r`; with a preconditioner, `z=M^{-1}r`.  For an
n-dimensional SPD matrix, exact arithmetic terminates in <= n steps; in practice
the relative residual criterion is the stopping gate.  This is the readable,
matrix-free counterpart of dense `numpy.linalg.solve` / MATLAB `pcg`.

Verified (test_topology_krylov): dense SPD solve matches `numpy.linalg.solve` to
roundoff and converges within the matrix dimension; callable `matvec` path gives
the same answer; Jacobi-exact diagonal preconditioning converges in one step; and
non-SPD / invalid inputs fail loudly.
"""


NONLINEAR_LSQ = r"""
## Nonlinear least squares -- Levenberg-Marquardt (the `lsqnonlin` member)

`nonlinear_lsq.levenberg_marquardt(residual, x0, jac=None)` minimises 0.5||r(x)||^2 by
interpolating Gauss-Newton and gradient descent:

    (J^T J + lam I) delta = -J^T r ,   x <- x + delta ,

damping lam DOWN on a successful step, UP on a rejected one; the Jacobian is supplied or
formed by forward differences. It is the NONLINEAR complement of the linear inverse solvers
(`linear_inverse`: TSVD / Tikhonov, for A x = b) and of the derivative-free `outer_loop`
(Nelder-Mead): use LM when the model is nonlinear in the parameters but smooth and the
residual is a sum of squares (curve/field fitting, calibration, inverse design).

Magnet/coil uses: recover excitation/geometry from measured field samples, B-H curve fits,
over-determined nonlinear calibration -- the field-fitting cousin of the analytic/linear
`field_synthesis` inverse.

SCALING CAVEAT (learned): the gradient stop ||J^T r||_inf < gtol is ABSOLUTE, so a residual in
tiny physical units (B ~ 1e-5 T) can halt in a flat valley far from the optimum. NORMALISE the
residual to O(1) (divide by a characteristic scale); then noiseless data is recovered exactly.

Verified (test_topology_nonlinear_lsq): exponential-model parameter recovery to ~1e-15;
circular-loop on-axis field B_z=mu0 I a^2/(2(a^2+z^2)^{3/2}) -> (I,a) recovered to ~1e-14 after
normalisation; Rosenbrock residuals [1-x, 10(y-x^2)] -> (1,1), cost 0; and agreement with
scipy.optimize.least_squares(method='lm') to ~1e-15.
"""


GLOBAL_OPTIMIZERS = r"""
## Global optimization -- Differential Evolution (the multimodal/global member)

`global_optimizers.differential_evolution(func, bounds, seed=...)` is a GLOBAL, derivative-free
population optimizer (Storn-Price DE/rand/1/bin): mutate v = a + F(b-c) from three random members,
binomially cross with the target, keep the trial only if it lowers the objective. It is the global
complement to the LOCAL optimizers (Nelder-Mead in `outer_loop`, Levenberg-Marquardt in
`nonlinear_lsq`) and the LINEAR inverse (`linear_inverse`): reach for it when the objective is
MULTIMODAL / non-convex (competing harmonics in a magnet/coil layout, a boxed design space) where
a descent method started anywhere stalls in a local minimum.

SCALING CAVEAT (learned): the budget grows with dimension -- popsize~15/dim & a few hundred gens
solve 2-D multimodal problems, but a 5-D Rastrigin needs popsize~25/dim and ~2000 gens to reliably
hit the global optimum (smaller budgets stall one Rastrigin "ring" away). Seed the RNG for a
deterministic run.

Verified (test_topology_global_opt): global optimum of Rastrigin (2-D and 5-D), Ackley (3-D) and
Rosenbrock (4-D) to f < 1e-6 at the known minimizer, deterministic for a fixed seed, and matching
scipy.optimize.differential_evolution.
"""


def get_applications_documentation(topic: str = "all") -> str:
    """Dispatch by topic.

    Topics:
      "all"            - All applications
      "motor"          - IPM motor shape optimization (Gangl 2015 case study)
      "field_synthesis"- Analytic / linear-inverse branch: PM-multipole
                         magnetization inverse (verified) + air-gap harmonic
                         objective + pointer to the streamfunction coil branch
      "linear_inverse" - Regularized linear inversion (TSVD / Tikhonov filter
                         factors + L-curve), the stable solver behind
                         field_synthesis, verified vs the SVD pseudo-inverse
      "krylov"         - Linear Conjugate Gradient / PCG for matrix-free SPD
                         inner solves in optimization and inverse design
      "outer_loop"     - Derivative-free outer-loop optimizers (Nelder-Mead +
                         fminsearchbnd bound transform; pointer to evolutionary
                         / official optuna-mcp for population/global), verified on Rosenbrock
      "nonlinear_lsq"  - Levenberg-Marquardt nonlinear least squares (lsqnonlin):
                         the nonlinear sum-of-squares solver for field/curve fitting
                         and inverse design, verified vs known optima + scipy
      "global_optimizers" - Differential Evolution: global, derivative-free population
                         optimizer for MULTIMODAL/non-convex objectives where local
                         methods stall, verified on Rastrigin/Ackley/Rosenbrock + scipy
    """
    t = topic.lower().strip()
    if t == "all":
        return (MOTOR_OPTIMIZATION + "\n\n" + FIELD_SYNTHESIS
                + "\n\n" + LINEAR_INVERSE + "\n\n" + OUTER_LOOP_OPTIMIZERS
                + "\n\n" + KRYLOV_SOLVERS + "\n\n" + NONLINEAR_LSQ
                + "\n\n" + GLOBAL_OPTIMIZERS)
    if t in ("motor", "ipm"):
        return MOTOR_OPTIMIZATION
    if t in ("field_synthesis", "field-synthesis", "fieldsynthesis",
             "synthesis", "inverse_source", "inverse-source", "magnet",
             "magnet_design", "pm", "multipole", "arbitrary_field",
             "stream_function", "streamfunction", "target_field",
             "spm", "halbach", "harmonics"):
        return FIELD_SYNTHESIS
    if t in ("linear_inverse", "linear-inverse", "linearinverse", "tsvd",
             "truncated_svd", "tikhonov", "regularization", "regularisation",
             "svd", "pinv", "pseudoinverse", "l_curve", "lcurve", "ill_posed"):
        return LINEAR_INVERSE
    if t in ("krylov", "linear_cg", "conjugate_gradient", "cg", "pcg",
             "hestenes_stiefel", "hestenes-stiefel", "spd_solver",
             "matrix_free", "matrix-free"):
        return KRYLOV_SOLVERS
    if t in ("outer_loop", "outer-loop", "outerloop", "derivative_free",
             "derivative-free", "nelder_mead", "nelder-mead", "neldermead",
             "simplex", "fminsearch", "fminsearchbnd", "direct_search",
             "optimizer", "optimizers"):
        return OUTER_LOOP_OPTIMIZERS
    if t in ("nonlinear_lsq", "nonlinear-lsq", "nonlinearlsq", "lsqnonlin",
             "levenberg_marquardt", "levenberg-marquardt", "levenberg", "lm",
             "least_squares", "least-squares", "gauss_newton", "gauss-newton",
             "nonlinear_least_squares", "curve_fit", "fit"):
        return NONLINEAR_LSQ
    if t in ("global_optimizers", "global", "global_optimization", "differential_evolution",
             "differential-evolution", "de", "evolutionary", "population", "swarm",
             "pso", "stochastic", "multimodal", "global_search"):
        return GLOBAL_OPTIMIZERS
    return ("Unknown topic '%s'. Available: all, motor, field_synthesis, "
            "linear_inverse, krylov, outer_loop, nonlinear_lsq, global_optimizers." % topic)
