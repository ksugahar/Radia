"""Kanamori et al. (2016) "Continuous Optimization for Machine Learning"
textbook companion -- chapter-by-chapter knowledge for radia_mcp.optuna.

Source PDF:
    W:/03_文献・論文/01_機械学習と最適化/00_教科書/20_機械学習_基礎/
    機械学習のための連続最適化.pdf
    EasyOCR ja+en, 300 dpi, 354 pages (lab-OCR'd 2026).
    Authors: Kanamori Takafumi, Suzuki Taiji, Takeuchi Ichiro, Sato Issei
    Publisher: 講談社サイエンティフィク (Kodansha Scientific) MLP Series, 2016.

Why this is in radia_mcp.optuna:
    radia_mcp.optuna covers black-box optimization (Optuna, BBO).
    This textbook covers gradient-based + 2nd-order + duality + KKT --
    the THEORETICAL FOUNDATION for when gradient methods apply (and
    therefore when BBO is NOT the right choice).  Use this module as
    a reference for: deciding gradient vs BBO, understanding the
    inner-loop optimization inside many BBO-targeted problems
    (e.g., FEM inner solves), and the SVM/sparse/matrix-optimization
    chapters that connect to lab ML work.

Chapter map (from the OCR'd TOC):
    Part I  Intro
      1: Machine-learning inference + optimization-problem taxonomy
      2: 基礎事項 -- Taylor, implicit function thm, sym matrix
         eigenvalues, projection, rank-1 updates, norms,
         convex sets/functions/derivatives + subdifferential
    Part II Unconstrained optimization
      3: Optimality conditions + termination criteria
      4: Gradient descent (Armijo / Wolfe line search)
      5: Newton's method (mod. Newton, Levenberg-Marquardt)
      6: Conjugate gradient (Fletcher-Reeves, Polak-Ribiere)
      7: Quasi-Newton (BFGS, L-BFGS, DFP)
      8: Trust-region methods (Cauchy, dogleg, Steihaug)
    Part III Constrained optimization
      9: Equality-constrained optimality + sensitivity analysis
      10: Inequality-constrained optimality (KKT)
      11: Primal methods (barrier / interior-point)
      12: Dual methods (Lagrange duality, dual ascent)
    Part IV ML applications
      13: Iterative reweighted / MM / EM
      14: SVM (linear + kernel, dual problem, SMO)
      15: Sparse learning (LASSO, ISTA/FISTA, ADMM)
      16: Matrix-space optimization (nuclear norm, low-rank, SDP)
"""


_CHAPTERS = {

# ============================================================
# Part I -- Introduction
# ============================================================

"chapter01": r"""
# Chapter 1 -- Machine-learning inference and optimization

The book opens by classifying ML problems as one of:

| Problem | Variable | Typical loss |
|---|---|---|
| Regression | continuous y | squared, Huber |
| Classification | discrete y | logistic, hinge |
| Density estimation | distribution | log-likelihood |
| Structured prediction | combinatorial | structured hinge |
| Reinforcement learning | policy | -E[return] |

Optimization problem class:

    min_x f(x)            unconstrained
    min_x f(x)  s.t. g_i(x) <= 0,  h_j(x) = 0   constrained
    min_x f(x)  s.t. x ∈ C                       set-constrained

Convergence-rate vocabulary used throughout the book:
    linear:        ||x_{k+1} - x*|| <= rho ||x_k - x*||,  rho < 1
    superlinear:   rho_k -> 0
    quadratic:     ||x_{k+1} - x*|| <= C ||x_k - x*||^2
    sublinear:     1/k, 1/sqrt(k)   (typical of pure GD on smooth fns)

Decision rule for lab usage:
- Gradient available + convex/smooth          -> Ch 4-8 (gradient methods)
- Gradient unavailable / FEM black-box        -> radia_mcp.optuna (BBO)
- Discrete combinatorial                      -> radia_mcp.evolutionary
- Posterior sampling                          -> radia_mcp.mcmc
""",

"chapter02": r"""
# Chapter 2 -- Mathematical preliminaries

§2.1 Calculus / Linear algebra:
- Taylor expansion (1st, 2nd order); used in every method derivation.
- Implicit function theorem -- foundation of sensitivity analysis
  (Ch 9): given F(x,y) = 0 with ∂F/∂y invertible, can solve y(x).
- Symmetric matrix eigendecomposition; positive (semi)definite tests.
- Projection onto a subspace: P = A(A^T A)^{-1} A^T  (the hat matrix
  in least squares); P^2 = P, P^T = P.
- Rank-1 update: (A + uv^T)^{-1} = A^{-1} - A^{-1} u v^T A^{-1} /
  (1 + v^T A^{-1} u)  (Sherman-Morrison).  Used to derive BFGS (Ch 7).
- Norms: vector L^p (p = 1, 2, ∞), induced matrix norms, Frobenius,
  nuclear norm (= sum of singular values, Ch 16).

§2.2 Convex analysis:
- Convex set: tx + (1-t)y ∈ S for t ∈ [0,1].
- Convex function: f(tx + (1-t)y) <= t f(x) + (1-t) f(y).
- Equivalent (smooth case): Hessian ∇²f >= 0 PSD.
- Convex => every local min is global.  Why convex problems are nice.
- Subdifferential ∂f(x) at non-smooth points (e.g. |x| at 0):
  g ∈ ∂f(x) iff f(y) >= f(x) + g^T(y - x) for all y.
  Used in proximal methods (Ch 15) and SVM duality (Ch 14).
- Convex conjugate f*(p) = sup_x (p^T x - f(x)) -- the Legendre-
  Fenchel transform.  Fundamental tool for Fenchel duality (Ch 12).
""",

# ============================================================
# Part II -- Unconstrained optimization
# ============================================================

"chapter03": r"""
# Chapter 3 -- Optimality conditions + termination

Unconstrained smooth case:
- Necessary 1st-order:  ∇f(x*) = 0  (stationary point)
- Necessary 2nd-order:  ∇²f(x*) >= 0  (PSD)
- Sufficient:           ∇f(x*) = 0  AND  ∇²f(x*) > 0  (PD)
                        => x* is a strict local minimizer

Set constraint x ∈ C  (C closed convex):
- KKT-like: -∇f(x*) ∈ N_C(x*) (normal cone) iff x* is a local min
- Equivalent: <∇f(x*), x - x*> >= 0 for all x ∈ C

Termination criteria in practice (any of):
    ||∇f(x_k)|| < eps_grad                # gradient tolerance
    |f(x_{k+1}) - f(x_k)| < eps_func * (1 + |f(x_k)|)
    ||x_{k+1} - x_k|| < eps_step * (1 + ||x_k||)
    k > k_max                              # iteration cap

Lab pitfall: when comparing convergence across solvers, ALWAYS use the
SAME termination criterion.  A method that quotes its own ||∇f|| < 1e-6
may have walked 100x further than one that quotes ||delta x|| < 1e-6.
""",

"chapter04": r"""
# Chapter 4 -- Gradient descent

Update rule:    x_{k+1} = x_k - alpha_k ∇f(x_k)

Step-size strategies:
- Fixed (alpha_k = const): works iff f is L-smooth and alpha <= 1/L
- Diminishing (alpha_k = c/k): SGD-style
- Backtracking with Armijo condition:
      Start alpha = 1; while f(x_k - alpha g_k) > f(x_k) - c1*alpha*||g_k||^2:
         alpha *= rho                       # typical c1 = 1e-4, rho = 0.5
- Wolfe conditions (sufficient decrease + curvature):
      f(x + alpha p) <= f(x) + c1 alpha ∇f^T p     # Armijo
      |∇f(x + alpha p)^T p| <= c2 |∇f^T p|         # curvature; c2 ≈ 0.9

Convergence rates (convex L-smooth):
- Fixed step alpha = 1/L:   f(x_k) - f* = O(1/k)         (sublinear)
- Strongly-convex (μ > 0):  exp decay (1 - μ/L)^k

Accelerated gradient (Nesterov, not in this chapter but cross-ref):
- O(1/k^2) for smooth convex; O((1 - sqrt(μ/L))^k) for strongly convex.

Lab usage:
- Gradient available + cheap Hessian unavailable + 1e3-1e6 dim -> good
- Hessian cheap to factor -> use Newton (Ch 5)
- Quasi-1e7 dim with structure -> L-BFGS (Ch 7) almost always better
""",

"chapter05": r"""
# Chapter 5 -- Newton's method

Update rule:    x_{k+1} = x_k - [∇²f(x_k)]^{-1} ∇f(x_k)

Why:  Quadratic Taylor model at x_k, exact minimization of that model.
Convergence: locally quadratic (||x_{k+1} - x*|| <= C ||x_k - x*||^2).

Practical problems:
- ∇²f(x_k) may not be PD -> step direction can be uphill
- Globally divergent if started far from optimum
- O(n^3) factor per iteration -> infeasible for n > 1e4 unless sparse

Fix 1: Modified Newton
    H_k = ∇²f(x_k) + lambda I       choose lambda so H_k >> 0
    Solve H_k d_k = -∇f(x_k); add line search.

Fix 2: Levenberg-Marquardt (for least-squares problems)
    f(x) = 1/2 ||r(x)||^2
    J = ∂r/∂x   (Jacobian; cheap from autodiff)
    (J^T J + mu I) d_k = -J^T r        # mu adapts: increase on bad step
    Interpolates between Gauss-Newton (mu -> 0) and gradient descent
    (mu large).  De facto standard for nonlinear least-squares.

Cross-ref to lab:
- Inverse problems (BH-curve fit, Play hysteresis fit): LM is the right
  starting point.  See radia.hantila_solver for the Hantila method
  (Newton-flavored fixed-point with line search) as a related pattern.
- Optuna's BoTorch GP-EI inner loop uses L-BFGS, not Newton; the
  surrogate Hessian is too expensive at high dim.
""",

"chapter06": r"""
# Chapter 6 -- Conjugate gradient (CG)

CG for quadratic min: f(x) = (1/2) x^T A x - b^T x  with A symmetric PD.
Solves Ax = b in <= n steps exactly (in exact arithmetic).

Key idea: search directions p_k are A-conjugate: p_i^T A p_j = 0 for i != j.
Each new direction is built from current gradient + previous direction:
    p_{k+1} = -g_{k+1} + beta_k p_k

Nonlinear CG variants (for general convex f) differ in how beta_k is set:
- Fletcher-Reeves (FR):     beta = ||g_{k+1}||^2 / ||g_k||^2
- Polak-Ribiere (PR):       beta = g_{k+1}^T (g_{k+1} - g_k) / ||g_k||^2
- Polak-Ribiere+ (PR+):     beta = max(PR, 0)        # global convergence

Restart: every n iterations (or when beta < 0 in PR+), set beta = 0
(i.e., revert to steepest descent).

Memory: O(n) -- key advantage over BFGS O(n^2).
Convergence (quadratic case): k iterations achieve
    f(x_k) - f* <= 2 ((sqrt(kappa) - 1)/(sqrt(kappa) + 1))^{2k} (f(x_0) - f*)
where kappa = lambda_max / lambda_min  (condition number).

Lab usage:
- Inner solver for FEM (NGSolve uses CG with preconditioner) -- the
  bread-and-butter use of CG is solving K u = f, not nonlinear opt
- Compact-AMS / Compact-AMG preconditioners in radia.sparsesolv_ngsolve
  are designed to feed CG well -- see radia_mcp.matrix_solvers for the
  preconditioner side.
""",

"chapter07": r"""
# Chapter 7 -- Quasi-Newton (BFGS, L-BFGS, DFP)

Idea: approximate the Hessian H_k from gradient differences (no
second-derivative eval).  Update H_k -> H_{k+1} using rank-1 or rank-2
update that satisfies the secant equation:
    H_{k+1} s_k = y_k        where s_k = x_{k+1} - x_k,
                                    y_k = g_{k+1} - g_k

Common updates:
- DFP (Davidon-Fletcher-Powell): updates H directly
- BFGS (Broyden-Fletcher-Goldfarb-Shanno): updates H^{-1} = B directly:
      B_{k+1} = (I - rho s y^T) B_k (I - rho y s^T) + rho s s^T
      rho = 1/(y^T s)
- Symmetric rank-1 (SR1): cheap but can lose PD

BFGS is the standard.  Convergence: superlinear locally.

L-BFGS (limited memory): instead of storing n x n matrix, keep last m
(s_k, y_k) pairs (m = 5-20 typical).  Two-loop recursion computes the
search direction in O(m n).  Standard for large-scale (n = 1e4-1e7).

Implementation:
    scipy.optimize.minimize(f, x0, method='L-BFGS-B', jac=grad_f,
                            bounds=...)

Lab usage:
- Inverse problems with ~1e3-1e4 parameters (e.g., per-element
  magnetization fit): L-BFGS is the default.
- Most "general gradient minimizer" calls in the lab go to L-BFGS-B
  because it handles box bounds naturally.
- For 1e4+ dim, use ngsolve's BlockJacobi + CG -- L-BFGS still memory-
  bound on the (s_k, y_k) pairs.
""",

"chapter08": r"""
# Chapter 8 -- Trust-region (TR) methods

Alternative to line search: pick a TRUST RADIUS Delta_k, build a
quadratic model m_k(p) around x_k, minimize m_k SUBJECT TO ||p|| <= Delta_k:

    min_p   m_k(p) = f(x_k) + g_k^T p + (1/2) p^T B_k p
    s.t.   ||p|| <= Delta_k

Then evaluate actual decrease vs predicted:
    rho_k = (f(x_k) - f(x_k + p_k)) / (m_k(0) - m_k(p_k))
- rho_k > 0.75 and ||p_k|| ~ Delta_k:  Delta_{k+1} = 2 Delta_k   (expand)
- rho_k < 0.25:                        Delta_{k+1} = (1/4) Delta_k  (shrink)
- rho_k <= 0:                          reject step, keep x_{k+1} = x_k

Inner subproblem solvers:
- Cauchy point: steepest descent inside the trust region (cheap)
- Dogleg: combine Cauchy + Newton step on the dogleg path (PD B_k only)
- Steihaug (truncated CG): apply CG to model, terminate at TR boundary

Why TR > line search for non-convex:
- TR can REJECT a bad step (line search is committed once direction is
  picked)
- Naturally handles non-PD B_k (Cauchy point is always defined)
- More robust for ill-conditioned problems

Lab usage:
- scipy.optimize.minimize(method='trust-ncg') / 'trust-krylov' /
  'trust-exact' all use TR with different inner solvers.
- For nonlinear least-squares with bounds, scipy's 'trf'
  (trust-region-reflective) is essentially LM + TR.
""",

# ============================================================
# Part III -- Constrained optimization
# ============================================================

"chapter09": r"""
# Chapter 9 -- Equality-constrained optimality + sensitivity

Problem:    min_x f(x)    s.t.   h(x) = 0    (h: R^n -> R^m, m < n)

Lagrangian:  L(x, lambda) = f(x) + lambda^T h(x)

Necessary 1st-order (KKT for equality):
    ∇_x L = ∇f(x*) + sum_j lambda_j ∇h_j(x*) = 0
    h(x*) = 0
With LICQ (linear independence constraint qualification:
  {∇h_j(x*)} linearly indep) the lambda* is unique.

Necessary 2nd-order:
    d^T ∇²_x L d >= 0    for all d ∈ ker(∇h(x*))    (tangent space)

Sensitivity analysis (§9.4):
- Perturbed problem:  min f(x)  s.t.  h(x) = b
- f*(b) is the value function.  Then:  ∇_b f*(b) = -lambda*
- I.e., the Lagrange multiplier IS the marginal sensitivity of optimal
  value to constraint right-hand side.

This connects to:
- Adjoint method in PDE-constrained optimization (Ch 11): exactly the
  same identity, used at every step of shape optimization.
- Economic interpretation in SVM: dual variables = "support" weight.
""",

"chapter10": r"""
# Chapter 10 -- Inequality-constrained optimality (KKT)

Problem:  min f(x)  s.t.  g_i(x) <= 0 (i=1..m), h_j(x) = 0 (j=1..p)

Lagrangian:  L = f + sum_i mu_i g_i + sum_j lambda_j h_j

KKT conditions (necessary for local min under constraint qualif.):
    Stationarity:        ∇_x L = 0
    Primal feasibility:  g(x*) <= 0,  h(x*) = 0
    Dual feasibility:    mu* >= 0
    Complementary slack: mu_i* g_i(x*) = 0  (each i)

Constraint qualifications (CQs) that make KKT necessary:
- LICQ: active constraints have linearly independent gradients
- MFCQ (Mangasarian-Fromovitz): weaker, often used
- Slater's: for convex problems, exists strictly feasible point
  (g_i(x) < 0 for nonlinear inequality g_i)

If problem is convex AND Slater holds, KKT is necessary AND sufficient.

The KKT system is what every interior-point solver (Ch 11) tries to
satisfy numerically.

Lab usage:
- Box constraints (a <= x <= b) are inequality constraints -> KKT gives
  the projection rule used by L-BFGS-B and TRF.
- Nonlinear inequality constraints in topology optimization (volume
  bound) -> handled by either Lagrangian or barrier (Ch 11).
""",

"chapter11": r"""
# Chapter 11 -- Primal methods: barrier / interior-point

Problem:    min f(x)   s.t.   g_i(x) <= 0

Barrier function method (replaces inequality with smooth barrier):
    phi_mu(x) = f(x) - mu sum_i log(-g_i(x))     # log barrier
    or
    phi_mu(x) = f(x) + mu sum_i 1/(-g_i(x))      # inverse barrier

Properties:
- phi_mu has a unique min x*(mu) in the strict interior {g_i < 0}
- x*(mu) is "the central path"
- As mu -> 0+, x*(mu) -> x*, the constrained optimum
- ||x*(mu) - x*|| = O(mu) (linear central path tracking)

Algorithm sketch:
1. Pick mu_0 (large enough that x*(mu_0) is interior)
2. Solve min phi_mu(x) by Newton -- inner loop
3. Reduce mu (e.g. mu := mu / 10)
4. Use previous x* as warm start; go to 2
5. Stop when ||grad of barrier|| < tol AND mu < eps_mu

Primal-dual interior-point (modern variant):
- Treat (x, mu_i, lambda_j) as joint primal-dual variables
- Newton step on entire KKT system
- Polynomially convergent for LP/QP/SOCP/SDP (Karmarkar 1984 etc.)

Lab usage:
- IPOPT (large-scale NLP) uses interior-point with filter line search.
- CVXPY / MOSEK / SCS use interior-point for convex problems
  (LP/QP/SOCP/SDP).
- For topology optimization with volume constraint, MMA
  (Method of Moving Asymptotes, Svanberg 1987) is more common than
  pure IP; see radia_mcp.topology_optimization for details.
""",

"chapter12": r"""
# Chapter 12 -- Dual methods (Lagrange duality)

For the inequality-constrained problem:
    primal:  min_x  f(x)   s.t. g(x) <= 0, h(x) = 0
    dual:    max_{mu>=0, lambda}  q(mu, lambda) := inf_x L(x, mu, lambda)

Weak duality:  d* := max q  <=  min f =: p*   ALWAYS
Strong duality: d* = p*    HOLDS for convex problems with Slater's
                            constraint qualification.

Why care about dual:
- Sometimes the dual is much easier (e.g., SVM dual is QP over n
  variables = #training points, simpler than primal which is over
  feature space.  Ch 14 elaborates).
- Bounds on optimal value: q(mu, lambda) <= p* always.
- Sensitivity = dual variable (Ch 9 §9.4 in equality case).

Dual ascent algorithm:
1. Start with mu_0 >= 0, lambda_0
2. Solve primal: x_k = arg min_x L(x, mu_k, lambda_k)
3. Update dual: mu_{k+1} = max(0, mu_k + alpha g(x_k))
                lambda_{k+1} = lambda_k + alpha h(x_k)
4. Repeat.

Augmented Lagrangian (Powell-Hestenes-Rockafellar):
    L_aug = L + (rho/2) ||h(x)||^2 + (rho/2) ||max(0, g(x) + mu/rho)||^2

Penalty parameter rho gets updated, multipliers updated by gradient
ascent.  Robust replacement for pure dual ascent.

Lab usage:
- ADMM (Ch 15) is augmented-Lagrangian-based.
- SVM SMO (Ch 14) is coordinate dual ascent over (alpha_i, alpha_j).
- For radia-mcp.optuna, the connection: "constraints_func" returns
  inequality violations; Optuna's Deb-domination rule is a primal
  approach (no Lagrangian), but for expensive constraints, an outer
  augmented-Lagrangian loop with Optuna as inner solver works.
""",

# ============================================================
# Part IV -- Machine-learning applications
# ============================================================

"chapter13": r"""
# Chapter 13 -- Surrogate / MM / EM family

Majorization-Minimization (MM) framework:
- Replace hard f(x) with a tractable surrogate g(x | x_k) such that:
    g(x_k | x_k) = f(x_k)        (tight at current iterate)
    g(x | x_k) >= f(x)           (majorizes f everywhere)
- Iterate:  x_{k+1} = arg min_x g(x | x_k)
- Guaranteed monotone decrease: f(x_{k+1}) <= g(x_{k+1} | x_k) <=
  g(x_k | x_k) = f(x_k).

Special cases:
- EM (Expectation-Maximization) for latent-variable models:
    g(theta | theta_k) = E_{z|x, theta_k}[ log p(x, z | theta) ]
  This is the standard derivation of EM as MM with the lower bound from
  Jensen's inequality.
- Iteratively reweighted least squares (IRLS) for robust regression:
    g_k(x) = sum_i w_i(x_k) r_i(x)^2     (re-weight residuals each step)
- Proximal gradient (next chapter): linearize smooth part + keep
  non-smooth term.

Why MM is appealing:
- Each iteration solves an EASIER problem (often closed-form).
- Guaranteed monotone (no line search needed if surrogate is tight).
- Convergence analysis tied to majorization quality.

Lab usage:
- ESIM cell-problem in radia.esim is essentially MM: at each Karl
  iteration, the surface impedance is re-evaluated and a new linear
  problem is solved.
- Magnetic hysteresis Play/Energy models inside FEM Newton: each Newton
  step is an MM iteration with a quadratic surrogate.
""",

"chapter14": r"""
# Chapter 14 -- Support Vector Machines (SVM)

Linear SVM primal:
    min_w  (1/2) ||w||^2 + C sum_i max(0, 1 - y_i (w^T x_i + b))
    (hinge loss + L2 regularization)

Dual form (key insight: kernel trick lives here):
    max_alpha  sum_i alpha_i - (1/2) sum_{i,j} alpha_i alpha_j y_i y_j
                                                  K(x_i, x_j)
    s.t.       0 <= alpha_i <= C    (box constraint from C-SVM)
               sum_i alpha_i y_i = 0    (from b)
where K(x_i, x_j) = <phi(x_i), phi(x_j)> is the kernel.
Predict: f(x) = sum_i alpha_i y_i K(x_i, x) + b.

Sequential Minimal Optimization (SMO):
- Coordinate dual ascent on PAIRS (alpha_i, alpha_j) at a time.
- Constraint sum_i alpha_i y_i = 0 makes 1-variable steps infeasible,
  so update 2 at a time keeping the constraint.
- Sub-problem at each step is 1D quadratic; closed-form solution.
- Heuristic pair selection (maximum KKT violation).  This is what
  LIBSVM implements -- still state of the art for moderate n.

Kernel choices:
- Linear: K(x,y) = x^T y (when n_features < n_samples, primal is faster)
- RBF (Gaussian): K(x,y) = exp(-gamma ||x-y||^2)
- Polynomial: K(x,y) = (gamma x^T y + r)^d
- Custom (string kernels, graph kernels, ...)

Lab usage:
- Most lab ML work uses sklearn.svm.SVC; under the hood it's LIBSVM
  + SMO.
- For data > ~100k samples, switch to linear SVM (sklearn.svm.LinearSVC
  uses LIBLINEAR with OWL-QN / cordinate descent, ALL in the primal).
- Kernel SVM is rarely the right tool for FEM-data regression
  (smoothness is wrong); use Gaussian processes (radia_mcp.bayesian_opt)
  or neural surrogate (radia_mcp.pinn) instead.
""",

"chapter15": r"""
# Chapter 15 -- Sparse learning (LASSO, ISTA, FISTA, ADMM)

LASSO (Tibshirani 1996):
    min_w  (1/2) ||y - X w||^2 + lambda ||w||_1
    (L1 regularization promotes EXACT zeros in w)

Why L1 induces sparsity: subdifferential at 0 is the interval [-1, 1]
(non-singleton).  The optimum w_i* = 0 whenever |X_i^T r| <= lambda.

Proximal gradient (ISTA = Iterative Soft-Thresholding Algorithm):
For min f(w) + g(w) with f smooth, g non-smooth (here g = lambda ||.||_1):
    w_{k+1} = prox_{alpha g}( w_k - alpha ∇f(w_k) )
where the prox operator for L1 is the soft-threshold:
    soft(z, t) = sign(z) * max(|z| - t, 0)            # elementwise

FISTA (Fast ISTA, Beck-Teboulle 2009):
- Add Nesterov-style momentum to ISTA -> O(1/k^2) convergence (vs
  O(1/k) for plain ISTA).
- The default first-order method for L1 problems.

ADMM (Alternating Direction Method of Multipliers, Boyd et al. 2011):
    min f(x) + g(z)    s.t.   Ax + Bz = c
- Updates alternate x, z, mu (augmented Lagrangian dual).
- For LASSO: split into smooth (least squares on x) + non-smooth
  (soft-threshold on z) + multiplier update.
- Excellent for distributed / large-scale problems (each split solves
  independently).

Other sparse penalties:
- Elastic net: (alpha) ||w||_1 + ((1-alpha)/2) ||w||^2
- Group lasso: sum_g ||w_g||_2 (sparsity at group level)
- Nuclear norm (Ch 16): sum sigma_i of matrix W (matrix analog of L1)
- SCAD, MCP: non-convex penalties for less bias

Lab usage:
- sklearn.linear_model.Lasso / ElasticNet wrap coordinate descent.
- For "which input parameters matter" in design optimization: L1-
  regularized response surface fit is a cheap sensitivity proxy.
- Cross-ref radia_mcp.bayesian_opt for similar "feature relevance"
  via GP ARD lengthscales.
""",

"chapter16": r"""
# Chapter 16 -- Matrix-space optimization

Variables are matrices X ∈ R^{m x n}; objectives involve singular values.

Nuclear-norm minimization (matrix analog of LASSO):
    min_X  ||X||_*    s.t.   measurement constraints
where ||X||_* = sum_i sigma_i(X) (sum of singular values).
- Promotes LOW-RANK solutions (matrix analog of sparsity).
- Convex relaxation of rank(X) which is NP-hard.

Matrix completion (Netflix-style):
    min_X  ||X||_*    s.t.  X_{ij} = M_{ij}  for (i,j) ∈ Omega
Recovers full M from a few observed entries IF M is low-rank and
Omega is "good" (incoherence + enough samples).

Algorithms:
- Singular Value Thresholding (SVT):  soft-threshold sigma_i  values.
  This is the prox operator for nuclear norm.
- Proximal gradient using SVT -- ISTA/FISTA in matrix space.
- ALS (Alternating Least Squares) for fixed-rank factorization
  X = U V^T -- non-convex but fast.

Semidefinite programming (SDP):
    min  <C, X>    s.t.  <A_i, X> = b_i,  X >= 0  (PSD)
- Variables X are symmetric PSD matrices.
- Convex; solvable by interior-point methods (Ch 11).
- Many ML problems reduce to SDP (e.g., kernel learning, max-cut
  relaxation, sum-of-squares polynomial optimization).

Lab usage:
- Low-rank surrogates: when 50-D parameter space has effective dim
  ~5, a low-rank GP fit (BoTorch ARD) recovers this automatically.
- POD / reduced-basis methods for radia simulations: SVD of snapshot
  matrix; truncation = nuclear-norm penalty in spirit.
- Cross-ref radia_mcp.mor for Model Order Reduction lab context
  (Cauer Ladder Network is a structured low-rank approximation).
""",

}


_LAB_TIE_INS = r"""
# Lab tie-ins -- mapping textbook chapters to radia / NGSolve workflows

| Textbook chapter | radia / lab analog | Cross-ref MCP server |
|---|---|---|
| Ch 4 Gradient descent | scipy.optimize 'CG', topology MMA inner | radia_mcp.topology_optimization |
| Ch 5 Newton / LM | scipy 'trf' for inverse fits; Hantila solver | radia_mcp.magnetic_materials |
| Ch 6 Conjugate gradient | NGSolve CG preconditioned by Compact-AMS | radia_mcp.matrix_solvers |
| Ch 7 BFGS / L-BFGS | scipy 'L-BFGS-B' default for ~1e3-1e4 dim fits | -- |
| Ch 8 Trust region | scipy 'trust-ncg' / 'trust-krylov' | -- |
| Ch 9-10 KKT / sensitivity | adjoint method in shape derivative | radia_mcp.topology_optimization |
| Ch 11 Barrier / IPM | IPOPT, CVXPY/MOSEK for convex constrained | -- |
| Ch 12 Duality | SVM SMO; augmented Lagrangian for constrained BBO | -- |
| Ch 13 MM / EM | ESIM Karl iteration; Play hysteresis Newton inner | radia_mcp.ih |
| Ch 14 SVM | sklearn surrogate models (rare in lab) | -- |
| Ch 15 LASSO / ADMM | Sensitivity-as-LASSO feature selection | radia_mcp.bayesian_opt (ARD) |
| Ch 16 Nuclear norm / SDP | POD / RB surrogates; CLN MOR is structured low-rank | radia_mcp.mor |

When to USE the textbook vs Optuna:
- Gradient available + convex/smooth + dim < 1e6 -> textbook (Ch 4-8)
- FEM black-box (no gradient) -> Optuna (this server's recipes)
- Constrained + convex -> Ch 10-12 + CVXPY
- Constrained + non-convex -> Optuna constraints_func OR
  augmented-Lagrangian outer loop with Optuna inner
- Inverse problem with LM-friendly structure -> Ch 5 LM, not Optuna
"""


def get_kanamori2016_documentation(topic: str = "all") -> str:
    """Return curated chapter summary from the Kanamori (2016) textbook.

    Accepts: 'all', 'chapter01' .. 'chapter16', 'overview', 'lab',
    or aliases 'newton'/'cg'/'bfgs'/'trust'/'kkt'/'duality'/'svm'/
    'lasso'/'sparse'/'matrix'/'sdp' etc.
    """
    t = (topic or "all").strip().lower()
    # alias resolution
    aliases = {
        "overview": "_overview",
        "lab": "_lab",
        "intro": "chapter01",
        "preliminaries": "chapter02",
        "basics": "chapter02",
        "optimality": "chapter03",
        "termination": "chapter03",
        "gradient": "chapter04",
        "gd": "chapter04",
        "newton": "chapter05",
        "lm": "chapter05",
        "levenberg": "chapter05",
        "cg": "chapter06",
        "conjugate_gradient": "chapter06",
        "bfgs": "chapter07",
        "quasi_newton": "chapter07",
        "lbfgs": "chapter07",
        "trust": "chapter08",
        "trust_region": "chapter08",
        "equality": "chapter09",
        "sensitivity": "chapter09",
        "kkt": "chapter10",
        "inequality": "chapter10",
        "barrier": "chapter11",
        "interior_point": "chapter11",
        "ipm": "chapter11",
        "duality": "chapter12",
        "dual": "chapter12",
        "lagrangian": "chapter12",
        "mm": "chapter13",
        "em": "chapter13",
        "surrogate": "chapter13",
        "svm": "chapter14",
        "smo": "chapter14",
        "lasso": "chapter15",
        "sparse": "chapter15",
        "ista": "chapter15",
        "fista": "chapter15",
        "admm": "chapter15",
        "proximal": "chapter15",
        "matrix": "chapter16",
        "nuclear": "chapter16",
        "sdp": "chapter16",
        "low_rank": "chapter16",
    }
    t = aliases.get(t, t)

    if t == "all":
        parts = ["# Kanamori et al. (2016) Continuous Optimization for ML",
                 "",
                 _LAB_TIE_INS]
        for k in sorted(_CHAPTERS):
            parts.append(_CHAPTERS[k])
        return "\n\n".join(parts)
    if t in ("_overview", "overview"):
        return ("# Kanamori et al. (2016) -- chapter index\n\n"
                "Use topic='chapter01'..'chapter16' for full content, "
                "or aliases (newton/cg/bfgs/trust/kkt/duality/svm/"
                "lasso/sparse/matrix). Use 'lab' for the lab tie-ins "
                "table.\n\n"
                "Source PDF: W:/03_文献・論文/01_機械学習と最適化/"
                "00_教科書/20_機械学習_基礎/"
                "機械学習のための連続最適化.pdf "
                "(EasyOCR ja+en, 300 dpi, 354 pages).\n\n"
                + "\n".join(f"  - {k}" for k in sorted(_CHAPTERS)))
    if t in ("_lab", "lab"):
        return _LAB_TIE_INS
    if t in _CHAPTERS:
        return _CHAPTERS[t]
    return (f"Unknown topic: {topic!r}. Use 'all', 'overview', 'lab', "
            f"or 'chapter01' .. 'chapter16'. "
            f"Aliases: newton/cg/bfgs/trust/kkt/duality/svm/lasso/sparse/matrix.")
