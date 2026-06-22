"""Topology optimization knowledge — Nishiwaki-Kondoh-Yachi textbook +
Gangl-Sturm nonlinear-magnetostatics bundle, mirrored for the NGSolve
implementation track.

Source bundles (2026-05-26):
  1. 近藤継男・矢地謙太郎・西脇眞二『トポロジー最適化の基礎 — 弾性体ならびに
     熱流体関連工学諸問題への応用のために』(Corona Publishing 2024, 267 p.).
  2. Gangl-Sturm 3-paper bundle:
     - Part I sensitivity analysis (TU Wien/Graz, 2019)
     - Part II numerical realisation in NGSolve+ngsxfem (Strobl 2019)
     - IPM motor case study (SIAM SISC 2015, 37, B1002-B1025)

This module is the NGSolve-implementation mirror of the same topology-
optimization content, distilled from the textbook + Gangl-Sturm bundle above
into a concise, RAG-indexable form.  Continuous-learning pattern: each new
textbook / paper absorbed gets reflected into this NGSolve-side module.

Companion modules in radia_ngsolve.knowledge:
  - kelvin.py — needed for Gangl-Sturm's W/Q exterior cell problem
                (truncated nonlinear cell problem on truncated R^3)
  - mmm_core.py — MMM tangent reluctivity tensor for averaged adjoint
  - basis_functions.py — Nedelec HCurl edge elements (state space)
  - esim.py — SIBC / nonlinear surface impedance (related sensitivity
              framework for fluid-conductor TopOpt extensions)
"""


# ============================================================
# 1. Framework overview — the 6-section template + 6 dial numbers
# ============================================================

FRAMEWORK_OVERVIEW = """
# Topology Optimization — Sugahara Lab Framework Overview

## The 6-section template (Nishiwaki textbook ch.4-9)

Every chapter from ch.4 onwards follows the SAME six-section pattern.
This is the lab contribution and what makes the method generalise
cleanly to any new physics (heat, elasticity, fluid, magnetostatics):

```
§X.1  Design variable + nondim density
      Scalar field phi(x) on Omega (unrestricted) -> rho(x) in [0,1]
      via clip / Heaviside.  rho=0 (void), rho=1 (solid), 0<rho<1
      (interface band Omega*).

§X.2  Primal problem
      (a) State variable          T(x) / u(x) / (u, p, T)
      (b) Non-dimensionalisation   reference length L, velocity U
      (c) State equation           PDE (heat / equilibrium / NS / curl-curl)
      (d) Constitutive law         Fourier / Hooke / Newton / B-H
      (e) Boundary conditions      Dirichlet / Neumann / Robin
      (f) SIMP material interpolation:
            kappa(rho) = eps + (1-eps) rho^p     (heat / elasticity)
            alpha(rho) = alpha_max * q*(1-rho)/(q+rho)  (fluid Brinkman)

§X.3  Optimisation problem
      Lagrangian F = f + lambda*g + integral(adjoint*residual) dx
      KKT: lambda >= 0, g <= 0, lambda*g = 0

§X.4  Adjoint problem
      Pick adjoint variable so delta_state F = 0 identically.
      The adjoint PDE has the SAME stiffness operator as the primal
      (transposed for non-self-adjoint cases), with cost-gradient as
      source and modified BCs.  This is the universal pattern.

§X.5  Sensitivity + evolution equation
      With adjoint satisfied:
          S(x) = delta f / delta rho
               = -(d coeff/d rho) * (primal . adjoint) + lambda * dg/drho

      Update phi by reaction-diffusion PDE in pseudo-time:
          m * d phi / dt = C_diff * Laplacian(phi) - S(x, t)   in Omega
          d phi / dn = 0                                       on boundary

      Converged steady state d phi/dt = 0 IS a KKT stationary point of
      the Lagrangian.  Monotone descent guaranteed:
          d f / d t = -(1/m) integral S^2 dx <= 0.

§X.6  Application (self-adjoint examples)
      Self-adjoint = adjoint problem reduces to the primal -> u* = u
      (or T* = T).  Single primal solve gives both -> NO extra adjoint
      solve.  Compliance (elasticity), thermal-resistance (heat),
      dissipation (Stokes) all enjoy this property.  Non-self-adjoint
      examples: drag, lift, conjugate heat transfer, nonlinear magneto-
      statics with tracking objective.
```

## The 6 dial numbers

These ARE the entire dial set for any TopOpt problem in this textbook's
framework:

| Symbol  | Typical  | Role                                                  |
|---------|----------|-------------------------------------------------------|
| m       | 1.0      | Pseudo-time mass / damping                            |
| C_diff  | 1e-3     | Diffusion = filter radius R_f = sqrt(C_diff dt / m)   |
| p_SIMP  | 3        | Penalty exponent (drives binary 0/1)                  |
| eps_E   | 1e-5     | Ersatz floor (avoid singular K in void)               |
| s_H     | 0.5      | Heaviside half-width (level-set -> density mapping)   |
| C_dw    | from delta_if | Double-well coefficient (phase-field interface) |

Set these six numbers + write §X.2 (physics PDE) = entire dial set per
problem.

## Connection to NGSolve / Radia stack

For an NGSolve user:
  - §X.1 = scalar FES (`H1`, lin) for phi; physical rho via CF projection
  - §X.2 = standard FES + BilinearForm assembly with rho-dependent CF
  - §X.4 = same BilinearForm transposed (or symmetric for self-adjoint)
            with cost-gradient as LinearForm RHS
  - §X.5 = post-solve CF assembly + Helmholtz filter via BilinearForm
  - Helmholtz filter (App C) = standard regulariser, reuse the same FES

NGSolve's `BilinearForm`+`AssembleLinearization` machinery handles the
Newton tangent for nonlinear primal solves; the adjoint linear solve
reuses the same factored matrix.

## Bibliography quick-list

- Bendsoe & Kikuchi 1988 CMAME (homogenisation method)
- Bendsoe 1989 Struct Optim (SIMP)
- Borrvall & Petersson 2003 IJNMF (Brinkman penalisation for fluids)
- Wang et al 2003 CMAME / Allaire et al 2004 JCP (level set)
- Yamada-Izui-Nishiwaki-Takezawa 2010 CMAME (RD level set, lab signature)
- Sokolowski-Zochowski 1999 SIAM JCO (topological derivative)
- Sturm 2015 SIAM JCO (averaged adjoint, minmax Lagrangian)
- Gangl et al 2015 SIAM SISC B1002 (IPM motor case)
- Andreassen et al 2011 SMO (88-line MATLAB)
- Yaji et al 2016 JCP 307 (world-first LBM+TopOpt)
"""


# ============================================================
# 2. Textbook ch.1-2 — foundations (4 method families, OC, SCP)
# ============================================================

TEXTBOOK_CH1_2 = """
# Ch.1-2 — Foundations (西脇眞二 / Nishiwaki, Kyoto Univ)

## Ch.1 History of structural optimisation

Pre-1960  Fully-stressed design (FSD): recurrence-formula update; correct
          only in statically determinate problems but the recurrence
          idea survives as the modern "optimality criteria" method.

1960     Schmit: first true sizing optimisation, member dimensions as
          design variables, mathematical programming.  Canonical start.

1973     Zienkiewicz-Campbell: shape optimisation with boundary nodes
          as DV.  Suffered "波打ち現象" (sawtooth boundary) — Imam (1981)
          proposed (i) design element technique, (ii) super-curves,
          (iii) shape superposition.

2000s    Implicit-boundary shape opt: VOF (Abe-Kikuchi) and level-set
          (Wang, Allaire) move boundary via Hamilton-Jacobi
              d phi/dt + V * |grad phi| = 0
          Holes can VANISH but creating new ones requires special tricks
          (topological derivative).

Topology optimisation:  recast as material-distribution problem ->
boundary motion + hole nucleation/annihilation arise AUTOMATICALLY.

## Ch.2.1 Fixed design domain + characteristic function (Bendsoe-Kikuchi 1988)

Recast structural optimisation as material distribution.  Define a
fixed design domain D containing the unknown optimal Omega_a.

Characteristic function chi_Omega(x):  = 1 in Omega_a, = 0 in D \\ Omega_a.

Linear elastic equilibrium extended to fixed D:
    integral_D chi_Omega * E * eps(u) : eps(v) dOmega
        = integral_Gamma_t t.v dGamma

Integration domain D is FIXED -> sensitivities w.r.t. DV easy to derive.
Fundamental difficulty: chi_Omega in L^inf admits pathological
discontinuities (Cheng-Olhoff plate ribs, Haber bicycle-spoke) ->
optimal microstructure can be infinitely-fine perforated.  Classical
derivatives fail.  The "design space" must be RELAXED -> §2.2
homogenisation.

## Ch.2.2 Homogenisation design method (Bendsoe-Kikuchi 1988)

Replace discontinuous chi_Omega*E with a smoothed, homogenised tensor
E^H derived from a periodic MICROSTRUCTURE (unit cell) at infinitesimal
scale.  Two scales: macroscopic x and microscopic y = x/eps.

Cell-problem characteristic displacement chi(x, y) on unit cell Y:
    integral_Y E * (eps^y(chi_kl) - eps^y(.)) dY = 0

Homogenised elasticity tensor:
    E^H_ijkl = (1/|Y|) integral_Y {E_ijkl - E_ijpq * eps^y_pq(chi_kl)} dY

Microstructures used:  layered rank-n materials.  Rank-2 has orthogonal
layer families with analytic E^H.  Design vars: layer thicknesses a, b
+ rotation angle theta.

Result: optimum is LARGE GREY-SCALE REGIONS (perforated / composite).
Historically not manufacturable (1988), but ADDITIVE MANUFACTURING (post
2010) is reviving the method.

## Ch.2.3 Density method = SIMP (Bendsoe 1989, Sigmund SIMP)

Drop the microstructure.  Macro tensor:
    E^G(rho) = rho^p * E_0    (SIMP — Solid Isotropic Material with Penalisation)

  p = penalty exponent.  Larger p -> intermediate rho energetically
  expensive -> binary 0/1 driven.  p = 3 recommended (Sigmund-Bendsoe,
  Hashin-Shtrikman justification — SIMP curve is achievable by an
  isotropic microstructure for p >= 3).

  Fig: p=1 -> grey-scale; p=2, 3 -> crisp black/white.

RAMP variant (Stolpe-Svanberg 2001):
    E^G(rho) = rho / (1 + c*(1 - rho)) * E_0
  Avoids zero-slope at rho=0 that SIMP has (helps gradient solvers).

Pros: no microstructure assumption; trivially implementable
("88-line MATLAB", Andreassen 2011); extends to fluids/EM/heat (no
homogenisation needed).
Cons: E^G(rho) has no rigorous physical basis; assumes isotropic
response -> search space narrower than rank-n homogenisation.

For FLUID problems: alpha(rho_f) = alpha_max * q*(1-rho_f)/(q+rho_f)
(Borrvall-Petersson 2003; RAMP-style penalisation of inverse
permeability), since NS homogenisation is essentially impossible.

## Ch.2.4 Level-set method (Wang 2003, Allaire 2004, Yamada 2010)

Level-set function phi(x) represents structure boundary as its ZERO
LEVEL SET:
    phi > 0 in Omega_a,  phi = 0 on partial Omega_a,  phi < 0 in void

Yamada-school formulation uses PIECEWISE-CONSTANT phi (not signed
distance), because the design variable is updated by a REACTION-
DIFFUSION equation (Ch. 3) rather than Hamilton-Jacobi.

Regularisation by gradient-squared:
    f -> f + tau * integral_D |grad phi|^2 dOmega    (tau = diffusion coef)

The diffusion smooths phi and regularises the ill-posed problem.
Geometrical complexity of optimum is CONTROLLABLE via tau — large tau
-> simple, large-feature shapes (manufacturable); small tau -> fine
features, higher performance.  This trade-off knob is unique to the
RD/level-set route.

Fig: optimum is CRISPLY BLACK-WHITE, no grey-scale — boundary is the
analytic zero level set, ideal for CAD export.

## Ch.2.5 Stiffness-max formulation

Why not stress?  Stress is LOCAL: reducing it at one point can raise
it elsewhere; geometric corners give singularities (-> infinity), no
finite allowable.

Mean compliance (global stiffness measure):
    l_mc = integral_Gamma_t t.u dGamma

Since t is fixed during opt, minimising l_mc minimises work done by t
on boundary -> maximises stiffness.  l_mc = 2 * strain energy ->
always finite.  In practice the compliance optimum approximates a
fully-stressed design.

Volume constraint:  V = integral_D rho dOmega

Two dual formulations:
  (A) min compliance s.t. volume cap (classical academic)
  (B) min volume s.t. compliance cap (modern automotive light-weighting)

## Ch.2.6.1 Optimality Criteria (OC) — Bendsoe-Kikuchi 1988

Discretise: n density DV rho_i in [rho^L, rho^U], scalar constraint
g(rho) <= 0.  Lagrangian L = f + Lambda g + ...

When both side constraints inactive and dg/drho != 0 everywhere, KKT
gives  (df/drho_i) / (-Lambda dg/drho_i) = 1.  Iterated as fixed-point
recurrence (FSD-style):

    rho_i^(k+1) = rho_i^(k) * (B_i^(k))^eta
    B_i^(k)    = (df/drho_i) / (-Lambda^(k) dg/drho_i)

with damping eta ~ 0.5.  Lambda^(k) found by BISECTION so volume
constraint holds.

Move limit zeta ~ 0.1-0.2 keeps updates monotone-stable:
    rho^L(k) = max((1-zeta) rho_i^(k), rho^L)
    rho^U(k) = min((1+zeta) rho_i^(k), rho^U)

Applicability: requires dg/drho_i != 0.  Works extremely well for
compliance + volume (canonical problem) — fast + robust.  Fails for
fluid TopOpt and multiple-constraint problems -> use SCP.

## Ch.2.6.2 SCP — Sequential Convex Programming (CONLIN, MMA)

General TopOpt workhorse.  Structure: (i) convex approximation of f, g
around current iterate, (ii) separation into single-variable sub-
problems, (iii) solve via dual.

CONLIN (Fleury-Braibant 1986) intermediate variable:
    y_i = rho_i        if df/drho_i > 0  (linear part)
    y_i = 1/rho_i      if df/drho_i < 0  (reciprocal part — most
                                          conservative, matches the
                                          true f propto 1/rho for
                                          compliance)

For compliance, df/drho < 0 everywhere -> reciprocal approximation
-> excellent convergence.

MMA (Svanberg 1987) — moving asymptotes U_i, L_i:
    y_i = 1/(U_i - rho_i)   if df/drho_i > 0
    y_i = 1/(rho_i - L_i)   if df/drho_i < 0

By tuning U_i, L_i each iteration (further to damp oscillation, closer
to accelerate), MMA achieves robust convergence on essentially every
TopOpt problem.  TODAY MMA IS THE DEFAULT OPTIMISER in academic TO
codes.  88-line MATLAB code uses an MMA-derivative OC; `top88` ships
mmasub.m companion.
"""


# ============================================================
# 3. Textbook ch.3 — Evolution-equation framework (近藤継男)
# ============================================================

TEXTBOOK_CH3_EVOLUTION_EQ = """
# Ch.3 — Evolution-equation engine (近藤継男 / Kondoh)

The CORE chapter.  Every later chapter (ch.4-8) re-instantiates the
engine here verbatim — only the physics PDE in §X.2 and the
application objective in §X.6 change.

## 3.1 Three shape representations (density / phase-field / level-set)

(1) DENSITY METHOD: physical density rho(x) in [0,1].  Carry an
unrestricted auxiliary "pseudo-density" phi(x) (no bound), recover
rho via clip:  rho = max(0, min(phi, 1)).  Material interpolation
(SIMP, RAMP) drives intermediate values toward 0/1 at convergence;
interface band Omega* becomes thin.

(2) PHASE-FIELD METHOD: phase-field function varphi(x) in [-1, 1].
Density:  rho = (varphi + 1) / 2.  Double-well potential
    f_dw(varphi) = (C_dw / 4) * (varphi+1)^2 * (varphi-1)^2
enforces binarisation with thin interface.

(3) LEVEL-SET METHOD: phi(x) unrestricted.  Body boundary = zero
isosurface phi=0; phi>0 material, phi<0 void.  Density via smoothed
Heaviside  rho = H_PS(phi; s_H)  where H_PS is C^2 piecewise-quintic
in [-s_H, s_H].  Yamada's variant additionally clips to [-1, 1]
(§3.5.2).

In all three, the unrestricted scalar phi(x) is the DESIGN VARIABLE;
physical density rho(x) is a clipped / Heaviside-mapped projection of
phi(x).  The SAME engine (evolution equation, §3.2) updates phi for
all three families — only the concrete PDE differs.

## 3.2.1-3.2.2 Why use an evolution equation?

Classical alternatives: OC (§2.6.1) and SCP/MMA (§2.6.2).  Here the
textbook adopts a THIRD paradigm: solve the optimisation as the
CONVERGED STEADY STATE of an EVOLUTION EQUATION (発展方程式).

Formal sensitivities:
  - DESIGN sens     S_dsgn(x) = delta f[phi] / delta phi(x)
  - DENSITY sens    S_dens(x) = delta f[rho] / delta rho(x)
  - TOPOLOGY sens   S_topo(x) = (f[rho; Omega\\B_eps] - f[rho; Omega]) /
                                [delta rho(x) * (-4 pi eps^3 / 3)]  in 3D

  Related: S_dsgn = (drho/dphi) S_dens; both share sign because
  drho/dphi >= 0.

Pseudo-time evolution under simplest gradient flow:
    m * d phi/dt = -S(x, t),  phi(x, 0) = phi_init(x)

Monotone descent guarantee:
    df/dt = integral_Omega S * (d phi/dt) dx = -(1/m) integral_Omega S^2 dx <= 0

So f decreases monotonically and steady-state d phi/dt = 0 is a
stationary point of f.  Same descent argument works with S_dens
(because drho/dphi >= 0).  This is the unifying engine: every method
in ch.3 wraps a different PDE around -S as a reaction term, but
descent is preserved by construction.

## 3.2.3 Reaction-diffusion equation (the platform PDE)

For DENSITY method and LEVEL-SET method:
    m * d phi/dt = C_dif * Laplacian(phi) + Q(x, t)         in Omega
                   [diffusion]              [reaction]
    d phi/dn = 0                                            on boundary
    m = 1, C_dif ~ 1e-3 (typical)

DIFFUSION term acts as a spatial filter / regulariser; smooths
checkerboards and removes mesh-dependence.  Characteristic filter
radius R_filter ~ sqrt(C_dif * dt / m).

REACTION term is the OPTIMISATION DRIVER:
    Q(x, t) = -S_dens(x, t)    (density method, RD level-set)
    Q(x, t) = -S_topo(x, t)    (Yamada level-set)

Decoupled: "sensitivity drives change, diffusion regularises it".  No
manual sensitivity filtering needed (the diffusion takes care of it).

## 3.2.4 Allen-Cahn equation (phase-field platform)

Extends RD with a DOUBLE-WELL POTENTIAL TERM that drives binarisation
to phi = +-1:

Ginzburg-Landau free-energy functional:
    f_GL[phi] = integral_Omega {(C_dif/2) |grad phi|^2 + f_dw(phi)} dx
    f_dw(phi) = (C_dw/4) (phi+1)^2 (phi-1)^2

f_dw has two minima at phi = +-1 ("two wells").  Non-stationary driven
Allen-Cahn:
    m * d phi/dt = C_dif * Laplacian(phi)
                   - C_dw * (phi+1)*phi*(phi-1)
                   - S_dens(x, t)

The reaction term (phi+1)*phi*(phi-1) vanishes at phi in {-1, 0, +1},
sign flips across these zeros, so any deviation from +-1 is pushed
back to a well.

Interface thickness scaling:
    delta_if ~ tanh^(-1)(0.9) * 2 * sqrt(C_dif / C_dw)

Set C_dw = (2 tanh^(-1)(0.9) sqrt(C_dif) / delta_if)^2 to target a
specific interface band width.

## 3.5.2 Yamada's RD-based level-set method (LAB SIGNATURE)

Uses TOPOLOGY SENSITIVITY S_topo as the reaction term + clips
level-set to [-1, 1].  Semi-implicit time discretisation: treat
DIFFUSION IMPLICITLY (at t) and REACTION EXPLICITLY (at t-dt):

    m * [phi(x,t) - phi(x,t-dt)] / dt = C_dif * Laplacian(phi(x,t))
                                        - S_topo(x, t-dt)

Rearranging gives a HELMHOLTZ-TYPE EQUATION at each time step:

    R^2 * Laplacian(phi(x,t)) - [phi(x,t) - phi_tentative(x,t)] = 0
    R = sqrt(C_dif * dt / m)

Each step = one linear Helmholtz solve.  Solution interpreted as
local average of phi_tentative on radius R -> R is the EXPLICIT
filter radius.

Distinguishing features:
  (a) Topology sensitivity as driver permits HOLE NUCLEATION inside
      the bulk — true topology optimisation, not just boundary motion.
  (b) [-1, 1] clip keeps signed-distance-like behaviour without
      separate reinitialisation.
  (c) Semi-implicit Helmholtz step is robust, linear, filter radius
      R is explicit and tunable through dt and C_dif.

Reference: Yamada-Izui-Nishiwaki-Takezawa, CMAME 2010.

## Framework summary

Four reusable pieces every later chapter re-instantiates verbatim:

  1. STATE PROBLEM: physics PDE parameterised by rho(x) via material
     interpolation (SIMP, RAMP).  State u (or T, p, ...) is functional
     of rho.

  2. ADJOINT PROBLEM: solve once per opt step to obtain Lagrange-
     multiplier lambda(x).  Density sensitivity is an EXPLICIT volume
     integral (no finite-difference needed).

  3. SENSITIVITY: choose S_dens (density / phase-field / RD level-set)
     or S_topo (Yamada).

  4. EVOLUTION PDE for phi:
       DENSITY      : m d phi/dt = C_dif Lap(phi) - S_dens
       PHASE-FIELD  : m d phi/dt = C_dif Lap(phi) - C_dw (phi+1)*phi*(phi-1)
                                                  - S_dens
       RD LEVEL-SET : m d phi/dt = C_dif Lap(phi) - S_dens
       YAMADA LEVEL : m d phi/dt = C_dif Lap(phi) - S_topo

     Notice the COMMON BACKBONE "m d/dt = C_dif Laplacian + reaction".
     Differences: which sensitivity + whether double-well (phase-field)
     or Heaviside clip (level-set) is applied.

Computational loop:
  (a) Solve state.
  (b) Solve adjoint.
  (c) Assemble sensitivity.
  (d) Advance design PDE one time step (usually ONE Helmholtz solve).
  (e) Recompute rho = clip(phi) or rho = H_PS(phi).
  (f) Repeat until ||d phi/dt|| < tol.
"""


# ============================================================
# 4. Textbook ch.4-5 — heat conduction + elasticity (self-adjoint)
# ============================================================

TEXTBOOK_CH4_5_HEAT_ELASTICITY = """
# Ch.4-5 — Heat conduction + elasticity (近藤継男 / Kondoh)

Both chapters follow the 6-section template (§X.1 design var, §X.2
primal, §X.3 opt, §X.4 adjoint, §X.5 sens + evolution, §X.6 app).
Both ch.4 thermal compliance and ch.5 mean compliance are
SELF-ADJOINT — single primal solve per iter, no adjoint solve needed.

## Ch.4 Heat conduction (self-adjoint)

State PDE (steady, fixed design domain Omega):
    -div(kappa(rho) grad T) = Q             in Omega
    T = T_bar                                on Gamma_T
    q = q_bar                                on Gamma_q
    kappa(rho) = eps_k + (1-eps_k) rho^p     (SIMP, eps_k ~ 1e-5, p = 3)

Thermal compliance (= thermal resistance minimisation):
    f = integral_Omega Q T dV
        + integral_Gamma_q q_bar T dS
        - integral_Gamma_T q T_bar dS

Self-adjoint when Q does not depend on T (eq. 4.45 dQ/dT = 0):
adjoint heat equation reduces to primal -> T*(x) = T(x).  Sensitivity:

    S(x) = -p (1-eps_k) rho^(p-1) * |grad T|^2 + lambda / V_bar    in Omega
    S    = 0                                                         on boundary

Wherever |grad T|^2 is large, S is most negative -> evolution adds
material there.  Tree-like "heat-conduction trees" emerge automatically
for 2D square domain with T-T / q-T / q-h / h-T / h-h boundary combos.

## Ch.5 Linear elasticity (self-adjoint, classical Bendsoe-Kikuchi-Sigmund)

State PDE (linear elastic equilibrium):
    div(sigma) + b = 0                       in Omega
    u = u_bar                                on Gamma_u
    t = t_bar                                on Gamma_t
    sigma = C(rho) : eps,  eps = (1/2)(grad u + (grad u)^T)
    E(rho) = eps_E + (1-eps_E) rho^p         (SIMP, eps_E ~ 1e-5, p=3)
    C(rho) = (rho^p) C_0                     (linear in E since nu fixed)

Mean compliance:
    f = integral_Gamma_t u.t_bar dS
        + integral_Omega u.b dV
        - integral_Gamma_u t.u_bar dS

Self-adjoint when body force b does not depend on u (eq. 5.48
db/du = 0): u*(x) = u(x).  Sensitivity:

    S(x) = -p (1-eps_E) rho^(p-1) * sigma:eps / E_solid + lambda / V_bar
         = -p (1-eps_E) rho^(p-1) * 2 * W_strain(u)    + lambda / V_bar

(twice the strain-energy density per element, modulated by SIMP slope)

Wherever strain energy is highest, S is most negative -> add material.
Wherever strain energy is near zero, S ~ lambda/V_bar > 0 -> remove
material.  Classical truss-like topologies emerge mirroring the
principal-stress trajectories.

This SELF-ADJOINT structure is the reason BKS density method became
the workhorse of TopOpt: each design iteration costs ONE primal FEM
solve.

## NGSolve implementation pattern

```python
from ngsolve import *

mesh = Mesh(...)
fes_state = H1(mesh, order=1) if scalar else VectorH1(mesh, order=1)
fes_rho   = L2(mesh, order=0)   # piecewise constant per element

gfu = GridFunction(fes_state)
rho = GridFunction(fes_rho)
rho.vec[:] = volfrac    # initial uniform density

# SIMP interpolation as a CoefficientFunction:
p_SIMP = 3
eps    = 1e-5
kappa_simp = eps + (1 - eps) * rho**p_SIMP   # heat case
# or:  E_simp = eps_E + (1 - eps_E) * rho**p_SIMP  for elasticity

# State problem
a_state = BilinearForm(fes_state, symmetric=True)
a_state += kappa_simp * grad(u) * grad(v) * dx   # heat
# or:    += SymbolicBFI on C(rho):eps(u):eps(v) for elasticity
f_state = LinearForm(fes_state)
f_state += Q * v * dx  # plus boundary terms

# Solve:
a_state.Assemble()
f_state.Assemble()
gfu.vec.data = a_state.mat.Inverse() * f_state.vec

# Self-adjoint sensitivity (heat, single primal solve)
sens_cf = -p_SIMP * (1 - eps) * rho**(p_SIMP - 1) * grad(gfu)*grad(gfu)
sens_gf = GridFunction(fes_rho)
sens_gf.Set(sens_cf)
# Now use sens_gf + Helmholtz filter + RD evolution for phi update.
```
"""


# ============================================================
# 5. Textbook ch.6-8 — fluid TopOpt (Toyota CRDL line)
# ============================================================

TEXTBOOK_CH6_8_FLUID = """
# Ch.6-8 — Fluid TopOpt (近藤継男 / Kondoh, Toyota CRDL)

Borrvall-Petersson 2003 + Toyota CRDL (近藤継男, 川本敦史, 松森唯益)
established this line; Kondoh is the textbook author.

## Brinkman penalisation — the foundational trick

Solid region modelled as a porous medium with extremely high inverse
permeability alpha (Brinkman resistance coefficient).  Body force:
    b = -alpha(eta) * u
The no-slip wall condition is simulated by an enormous Darcy-type
sink that suppresses fluid velocity wherever solid material is placed.

RAMP interpolation (Rational Approximation of Material Properties):
    alpha(eta) = alpha_max * q * (1 - eta) / (q + eta)
    alpha_max ~ 1e4,  q ~ 0.1 or 0.01

eta = 0 (solid) -> alpha = alpha_max (no-slip)
eta = 1 (fluid) -> alpha = 0          (pure Stokes / NS)

Because alpha enters momentum as a velocity-proportional sink, NO
body-fitted remesh is needed — the wall moves as eta evolves on a
FIXED Eulerian grid.  THIS is why fluid TopOpt works.

For laminar problems scale by (1 + 1/Re) to keep alpha*u comparable
to diffusion grad^2 u / Re at all Reynolds numbers.

## Ch.6 Stokes flow

State (steady, vanishing-Re creeping flow):
    div(u) = 0
    -div(sigma) + alpha(eta) * u = 0
    sigma = -p I + tau,  tau = 2 D,  D = (grad u + (grad u)^T)/2
    u = u_bar  on Gamma_u,  t = t_bar  on Gamma_t

LINEAR in u and p -> adjoint problem is self-adjoint for flow-resistance
minimisation:
    p' = -p,  u' = +u    (closed-form self-adjoint solution)

Sensitivity:
    S(x) = -(d alpha / d eta) |u|^2 + lambda / V_bar    in Omega

NO second linear solve required.  Computational examples: width-6
channel -> two parallel pipes, width-10 -> single merged pipe.
Stokes TopOpt automatically discovers branching topology.

Ch.6.6.2 Drag minimisation (NOT self-adjoint):
    D = integral_Omega [alpha(eta) * (u . e_x)] dOmega
        (domain-integral form via Brinkman trick — avoids surface
         integral over moving body surface)

Adjoint has non-zero source alpha(eta)*e_x in momentum equation +
homogeneous u' = 0 on Gamma_u.  Full adjoint Stokes solve required.

## Ch.7 Laminar Navier-Stokes (1 <= Re <= 1000)

State adds convective term:
    div(u) = 0
    rho (u . grad) u - div(sigma) + alpha(eta) u = 0
    sigma = -p I + (2/Re) D

(u . grad) u makes the PDE NONLINEAR.  Newton/Picard iteration per
design step.  Reynolds Re = U L / nu enters via diffusion coefficient.

Convective term BREAKS self-adjointness for almost all objectives.
The adjoint PDE has SIGN-FLIPPED convective term + extra cross-coupling
-rho * (grad u)^T . u'.  Per iteration:
  (1) NONLINEAR primal NS solve (Newton, frozen-Jacobian inner loops)
  (2) LINEAR adjoint solve with same Jacobian + cost-gradient source
  (3) Sensitivity = -(d alpha/d eta) (u . u') + ...

Application portfolio (ch.7.6.1-7.6.5):
  7.6.1  Energy dissipation min in channel (self-adjoint at low Re)
  7.6.2  Mechanical energy loss min (boundary = volume form by Bernoulli)
  7.6.3  Drag min (domain-integral via Brinkman; Newton + adjoint)
  7.6.4  Lift max (drag-upper-bound constraint, airfoil-like shapes)
  7.6.5  Lift-drag ratio max (quotient rule combines both adjoints
                              into ONE adjoint NS solve via linearity)

## Ch.8 Conjugate heat transfer (Toyota CRDL heatsink topology)

State: mass + momentum + energy
    div(u) = 0
    rho (u.grad)u - div(sigma) + alpha(eta) u = 0
    rho cp (u.grad T) - div(kappa(eta) grad T) = Q

The SAME eta(x) field controls BOTH flow blockage (Brinkman alpha)
AND thermal conductivity (SIMP-like kappa(eta)).  This is what makes
the problem CONJUGATE.

Adjoint has structure inherited from BOTH the NS adjoint and the
heat-conduction adjoint.  CROSS-COUPLING term -rho cp (grad T) T' in
adjoint momentum: dT/du nonzero so adjoint u-equation picks up
T'-driven source.

Sensitivity (combining both physics):
    S(x) = -(d alpha/d eta) (u . u')
           + (d kappa/d eta) (grad T . grad T')
           + delta_eta f + lambda delta_eta g

The d alpha/d eta term adds solid (eta -> 0) wherever u . u' is positive
(drag reduction); d kappa/d eta places high-conductivity solid wherever
grad T . grad T' is positive (heat removal).

Heatsink topology application (ch.8.6):
  Minimise integral_Gamma_heat T dGamma subject to volume cap.
  Fractal / tree-like / dendritic fin patterns emerge AUTOMATICALLY.
  These are Bejan's constructal-theory shapes, arrived at without any
  geometric prior.  Toyota Prius hybrid power module cooling layout
  was designed exactly this way.

## Borrvall-Petersson 2003 + Toyota CRDL history

Borrvall & Petersson (Int. J. Numer. Meth. Fluids, 2003) introduced
the Brinkman penalisation trick to enable density-based TopOpt for
fluid problems on a fixed Eulerian grid.  Before this, fluid shape
optimisation required body-fitted remeshing per iteration -> gradient-
based optimisation impractical for topology changes.

Toyota CRDL (Nagakute, Aichi) team led by 近藤継男 (Kondoh — the
textbook author) and including 川本敦史 (Kawamoto), 松森唯益
(Matsumori) pioneered:
  - Laminar NS extension (Kawamoto et al 2011)
  - Conjugate heat transfer (Matsumori et al 2013)
  - Manufacturing-aware variants used in Toyota production:
    cold plates, oil galleries, exhaust gas recirculation passages
  - Density filter + Heaviside projection for crisp 0/1 designs
  - Continuation on alpha_max to avoid local minima at high Re
  - Natural convection (Boussinesq buoyancy) coupling
"""


# ============================================================
# 6. Textbook ch.9 — Lattice Boltzmann TopOpt (Yaji)
# ============================================================

TEXTBOOK_CH9_LBM = """
# Ch.9 — Lattice Boltzmann TopOpt (矢地謙太郎 / Yaji, Osaka Univ)

WORLD-FIRST combination of LBM + TopOpt (Yaji, Yamada, Yoshino,
Matsumoto, Izui, Nishiwaki — JCP 307 (2016) 355-377).

## LBM basics

Discrete velocity distribution functions f_i(x, t) — density of
fictitious particles moving at velocity c_i — instead of solving NS
directly.  Lattice gas: D2Q9 (2D, 9 velocities) or D3Q15 (3D).

BGK (Bhatnagar-Gross-Krook) discrete Boltzmann:
    Sh * d f_i / dt + c_i . grad f_i = -(1/eps) (f_i - f_i^eq)

After discretising in space (Delta x) and time (Delta t):
    f_i(x + c_i Delta_t, t + Delta_t) = f_i(x, t)
                                         - (1/tau) [f_i(x, t) - f_i^eq(x, t)]

with tau = eps / Delta_t.  FULLY EXPLICIT update using only nearest-
lattice-neighbour values -> embarrassingly parallel, GPU-friendly,
natural for transient flow.

Macroscopic recovery: rho = sum f_i, rho u = sum c_i f_i, p = rho/3.
Equilibrium (so Chapman-Enskog -> NS):
    f_i^eq = E_i rho [1 + 3 c_i.u + (9/2)(c_i.u)^2 - (3/2) u^2]

with weights {4/9, 1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36, 1/36} for D2Q9.

Kinematic viscosity: nu = (1/3)(tau - 1/2).  Require tau > 0.5.

## TopOpt under LBM

Define eta : Omega -> [0, 1].  Solid is modelled by injecting Brinkman-
type body force:
    g = -alpha(eta) u
    alpha(eta) = chi [1 + (1/Re)] (q+1)(1-eta) / (q+eta)
                  chi ~ 1e4, q in [0.01, 1]

Two-step fractional-step inclusion in LBM:
  Step 1 — standard stream-collide (no force) -> f_i^*
  Step 2 — correct: f_i = f_i^* + 3 E_i rho c_i.g Delta_t

Generic functional over SPACE AND TIME (key feature):
    F = integral_I integral_Omega f_Omega(u, p, eta) dOmega dt
        + integral_I integral_Gamma f_Gamma dGamma dt

The time-integration window I = [t_0, t_1] makes LBM the natural engine
for UNSTEADY-flow TopOpt.

## ALBM (adjoint LBM) — continuous-adjoint route

The continuous-adjoint route (Krause 2013 + Yaji 2017) is preferred
over discrete-adjoint (Tekitek 2006) because the bounce-back / velocity
/ pressure boundary conditions become tractable.

Extended functional with adjoint distribution hat-f_i:
    F-hat = F + sum_i integral integral hat-f_i R_i dOmega dt
    R_i = Sh d f_i/dt + c_i.grad f_i + (1/eps)(f_i - f_i^eq) + 3 E_i rho c_i.g

Integration by parts in t and x:
    integral hat-f d_t f dt = [hat-f f]_t0^t1 - integral d_t hat-f * f dt
    ... + divergence theorem on c_i . grad f_i

yields the ADJOINT BOLTZMANN EQUATION (sign-REVERSED time + advection):
    -Sh d hat-f_i/dt - c_i.grad hat-f_i
        = -(1/eps)(hat-f_i - hat-f_i^eq) - 3 E_i B_i c_i.g + (dF/df_i)

with adjoint equilibrium hat-f_i^eq = hat-rho + 3 c_i.hat-u where
    hat-rho = sum E_i hat-f_i [1 + 3 c_i.u + ...]
    hat-u   = sum c_i hat-f_i

Sign-reversed time and advection: ADJOINT PROPAGATES BACKWARDS in time
and along -c_i.  Same parallel structure as forward LBM -> same scalability.

Initial condition: hat-f_i(t_1) = 0 (terminal-time condition).
Boundary conditions inherit from forward BC family with UNKNOWN
distributions reversed.

Sensitivity (post-adjoint):
    dF/d eta = integral_I integral_Omega [sum_i (3 E_i rho c_i.hat-f_i) *
                (-d alpha/d eta) u + d_eta f_Omega] dOmega dt + boundary

## Why LBM matters for unsteady TopOpt

Pre-2016, fluid TopOpt was almost exclusively STEADY Stokes/laminar.
Unsteady NS-TopOpt was considered impractical because:
  (a) Backward-in-time adjoint requires full forward history ->
      memory blow-up.
  (b) Implicit time-integrators bury sensitivity in nonlinear iter.

LBM solves both:
  (a) Explicit stream-collide makes per-step memory predictable ->
      checkpointing easy.
  (b) Adjoint LB is structurally identical (just sign-reversed time/
      advection) -> inherits same parallel scalability.

Local-in-time refinement (Yaji 2017 / 2018) further bounds memory by
windowing -> industrially-relevant transient flows (turbomachinery,
biological pumping, automotive cooling) tractable.

## Application examples (ch.9.5)

(1) STEADY pressure-loss min (sec 9.5.1): 2D bent-pipe with 50% volume
cap.  Converges to smooth S-curve.

(2) UNSTEADY pressure-loss min (sec 9.5.2): time-averaged objective
over I = [0, 2 pi].  Optimised 2D dual-pipe develops a VERTICAL
CONNECTING CHANNEL between the two horizontal pipes — synchronised
suction-and-blowing across phase.  THIS GEOMETRY IS ABSENT IN STEADY
TopOpt — the canonical demonstration justifying LBM-TopOpt.

(3) STEADY forced convection (sec 9.5.3): branched channel maximising
solid-fluid interface area -> high heat exchange.

(4) UNSTEADY forced convection (sec 9.5.4): 3D, 200x100x40 lattice,
4900 Delta_t.  Complex branching channel with internal void cavities
-> captures unsteady mass exchange that 3D steady solvers cannot.

## Reference

Yaji-Yamada-Yoshino-Matsumoto-Izui-Nishiwaki (2016), JCP 307, 355-377.
Follow-up: local-in-time adjoint for unsteady flows (MEJ 2017);
large-scale local-in-time unsteady thermal-fluid (Struct Multi Optim 2018).
"""


# ============================================================
# 7. Textbook appendices — App A cookbook + App C filter + App D KKT
# ============================================================

TEXTBOOK_APPENDIX_A_OBJECTIVES = """
# App A — 5-physics objective catalogue

Canonical objective forms (volume + boundary integrals) for the five
physics covered in ch.4-8.  Each comes with the explicit functional-
derivative-w.r.t.-rho formula -> sensitivity assembly is one
substitution.

## A.1 Heat conduction TopOpt

State: -div(k grad T) - Q = 0, T = T* on Gamma_T, q = q* on Gamma_q.
Property: SIMP k(rho) = k_void + (k_solid - k_void) rho^p, p ~ 3.

Self-adjoint THERMAL COMPLIANCE (three equivalent forms):
  f_comp1 = integral_Omega Q T dOmega - integral_Gamma_q q* T dGamma
            + integral_Gamma_T q T dGamma
  f_comp2 = integral_Omega Q T dOmega - integral_Gamma_q q* T dGamma
            - integral_Omega (k/2) |grad T|^2 dOmega
            (when Q=0 and q*=0, MAXIMIZE thermal-field strain energy)
  f_comp3 = -integral_Omega (k/2) |grad T|^2 dOmega + 2 integral_Gamma_T q T dGamma
            (when T=0 on Gamma_T, MINIMIZE thermal strain energy)

Self-adjoint: T-hat = T, -T, or 2T.  Sensitivity:
    d f / d rho = -(dk/d rho) |grad T|^2
    < 0 ⇒ add material where temperature gradient is high.

## A.2 Elasticity TopOpt

State: -div(sigma) - b = 0, sigma = C(rho):epsilon.
Property: SIMP E(rho) = E_void + (E_solid - E_void) rho^p.

MEAN COMPLIANCE (most common):
  f_comp1 = integral_Omega b.u dOmega + integral_Gamma_t t*.u dGamma
            - integral_Gamma_u t.u dGamma
  f_comp2 = integral_Omega b.u + integral_Gamma_t t*.u dGamma
            - 2 integral_Omega W dOmega,  W = (1/2) sigma:epsilon
  f_comp3 = -2 integral W dOmega + 2 integral_Gamma_u t.u dGamma

Self-adjoint: u-hat = u, -u, or 2u.  Sensitivity:
    d f / d rho = -p rho^(p-1) sigma:epsilon / E_solid

This is the formula behind every textbook "MBB beam -> triangulated
truss" demo (Sigmund 88-line MATLAB code).

## A.3 Stokes TopOpt

State: NS without inertia + Brinkman g = -alpha(eta) u.
Property: alpha(eta) = chi*[1 + (q+1)(1-eta)/(q+eta)] / Re.

FLOWABILITY MAX = FLOW-RESISTANCE MIN (three forms, SIGN OPPOSITE
elasticity — flowability MAX not MIN):
  f_comp1 = integral b.u + integral_Gamma_t t.u + integral_Gamma_u t.u
  f_comp2 = ... - 2 integral_Omega (Phi + alpha u.u) dOmega + 2 integral_Gamma_u t.u
            Phi = (1/(2 Re)) sum (d u_i/d x_j + d u_j/d x_i)^2  (viscous dissipation)

Self-adjoint: u-hat = -u, 0, or -2u.  Sensitivity:
    d f / d eta = -(d alpha / d eta) u.u

Positive where fluid velocity is high -> reduce alpha (make more fluid).

## A.4 Laminar TopOpt

Same eta and Brinkman g, but full NS adds (u.grad)u -> breaks self-
adjointness; real adjoint solve required.

Catalog:
(1) FLOW COMPLIANCE MAX (A.4.2): Stokes form + boundary kinetic-energy
    outflux integral_Gamma (|u|^2/2)(u.n) dGamma.
(2) VISCOUS DISSIPATION MIN (A.4.3): f_diss = integral (Phi + alpha u.u) dOmega.
(3) MECHANICAL ENERGY LOSS MIN (A.4.3): f_loss = integral_Gamma p_t(u.n) dGamma
    where p_t = p + |u|^2/2 (Bernoulli total pressure).  By Bernoulli on
    inviscid streamlines, f_diss ≈ f_loss; with uniform inlet/outlet
    velocity, both reduce to Delta_p * V_dot.
(4) DRAG MIN (A.4.4):
      D = integral_Omega alpha u_x dOmega                 (Brinkman integral)
      D = integral_Gamma_infty [-u_x(u.n) - p n_x + (T_v n)_x] dGamma  (far-field)
      D = -u_infty . integral g dOmega ≈ integral (Phi + alpha u.u) / u_infty dOmega
    Two-or-three equivalent expressions.
(5) LIFT MAX (A.4.4): Y = integral alpha u_y dOmega.
(6) LIFT-DRAG RATIO MAX (A.4.4 [7]): f = -Y/D, quotient rule combines
    both adjoint sources into ONE adjoint NS solve.

## A.5 Heat transfer TopOpt

State: NS + thermal energy div(T u) - div(q) + Q - g.u = 0,
q = -(k/Pe) grad T.

(A.5.2) FIXED-T SOURCE (CPU die at T=1):
  HEAT-TRANSFER MAX (five equivalent forms):
    H_1 = -integral_Gamma_heat q dGamma       (direct heat flux integral)
    H_2 = integral_Gamma_non-heat (T(u.n) + q) dGamma   (energy balance)
    H_3 = -integral_Gamma T q dGamma          (whole-boundary T-weighted flux
                                                = "thermal compliance" analog)
    H_4 = integral_Omega (k/Pe) |grad T|^2 dOmega
          + integral_Omega T (u.grad T) dOmega  (volume form via Gauss)
    H_5 ≈ integral_Gamma_adiabatic T (u.n) dGamma  (convective wake approx)
  Maximising H ≡ maximising dimensionless mean Nu.

(A.5.3) FIXED-FLUX SOURCE (heater dissipating prescribed power):
  Since total heat is fixed, MINIMIZE heater-surface mean TEMPERATURE M.
"""


TEXTBOOK_APPENDIX_C_FILTER = """
# App C — Helmholtz density filter (the canonical regulariser)

PROBLEM: density-method TopOpt generically produces a "checkerboard"
pattern at the mesh scale because the single-point sensitivity
delta rho -> delta f is not coupled to neighbours.

FIX: filter the design field (or sensitivity) at length scale > mesh.

## The Helmholtz-type filter PDE

    -R^2 Laplacian(rho-tilde) + rho-tilde = rho      in Omega
    R^2 (d rho-tilde / d n) = 0                       on Gamma

  rho:        raw (noisy / pixelated) design field
  rho-tilde:  smoothed output
  R > 0:      FILTER PARAMETER

Effective FILTER RADIUS r ≈ 3.5 R (verified numerically in 1D step-
function test: rho-tilde reaches {0.02, 0.98} at x = +-r).

## Derivation

Euler-Lagrange of the master functional:
    f[rho-tilde] = integral_Omega [(1/2)(rho-tilde - rho)^2
                                    + (1/2) d (grad rho-tilde)^2] dOmega
                                    , d = R^2 > 0

First term -> rho-tilde approximates rho.  Second -> rho-tilde smoothed.
Setting d f / d rho-tilde = 0 yields the Helmholtz PDE.

## Physical intuition

3D Green's function of Helmholtz operator:
    G(r; R) = exp(-r/R) / (4 pi R^2 r)

decays EXPONENTIALLY (vs 1/r for Laplace) -> R sets a hard locality
scale.  integral_{R^3} G dV = 1 confirms it is properly normalised
averaging kernel.

## Filter-radius guidance

  R = (1-3) Delta_x  for normal use
  R = (5-10) Delta_x  for aggressive smoothing (early iterations)
  R = laser_spot      for additive manufacturing min feature
  R = milling_tool    for subtractive machining min feature

Continuation: start with large R, decrease as opt converges (similar
to SIMP penalty continuation p = 1 -> 3).

## NGSolve implementation

```python
from ngsolve import *

fes_rho = L2(mesh, order=0)
fes_rho_tilde = H1(mesh, order=1)

rho       = GridFunction(fes_rho)         # raw design field
rho_tilde = GridFunction(fes_rho_tilde)   # filtered field

a_filter = BilinearForm(fes_rho_tilde, symmetric=True)
a_filter += R**2 * grad(u_tilde) * grad(v_tilde) * dx
a_filter += u_tilde * v_tilde * dx
a_filter.Assemble()

f_filter = LinearForm(fes_rho_tilde)
f_filter += rho * v_tilde * dx
f_filter.Assemble()

# Each opt iter:
rho_tilde.vec.data = a_filter.mat.Inverse() * f_filter.vec
# Use rho_tilde for state solve + sensitivity.
```

This is THE workhorse regulariser for density-method TopOpt and is
mesh-independent (R in physical units, not element count).
"""


TEXTBOOK_APPENDIX_D_KKT = """
# App D — Constrained TopOpt via KKT + dual evolution

The classical KKT conditions are realised via TWO COUPLED PSEUDO-TIME
EVOLUTION PDEs (one for phi, one for lambda) instead of an inner-loop
dual solve.

## KKT conditions

For min f[phi] subject to g[phi] <= 0:

  L[phi; lambda] = f[phi] + lambda * g[phi]

KKT (necessary for optimality):
  (1) Stationarity:        S(x) = S_f(x) + lambda S_g(x) = 0 for all x
  (2) Primal feasibility:  g[phi] <= 0
  (3) Dual feasibility:    lambda >= 0
  (4) Complementary slackness:  lambda * g[phi] = 0

## Design-variable update (App D.2)

Target-tracking form (the production-quality update):
    phi_tentative(x, t) = phi(x, t) - C_S * S(x, t)
    phi_target(x, t)    = max(-1, min(phi_tentative, 1))
    m * d phi/dt        = C_diff Laplacian(phi)
                          - C_dw (-phi^3 + phi)         (double-well)
                          + (phi_target - phi)          (target tracking)

C_diff << 1 (e.g. 1e-5); C_dw set to make interface ≈ 2 Delta_x:
    C_dw = -(2 C_diff / delta_if) * tanh^(-1)(0.9),  delta_if ≈ 2 Delta_x

Steady state d phi/dt = 0 ⇔ S = 0 ⇔ KKT (1).

## Lagrange-multiplier update (App D.3)

Target-tracking form (production-quality):
    lambda_tentative(t) = lambda(t)
                          + (C_lg / m_l) lambda(t) sinh(g(t)/g_ref1)
                          + (C_g / m_l) dg/dt
    lambda_target(t)    = max(eps_l, lambda_tentative)
    m_l * d lambda/dt   = lambda_target - lambda(t)

The sinh nonlinearity makes lambda respond FASTER to g far from
feasibility but cap'd near feasibility -> smoother numerically than
pure exponential growth (D.16 basic form).

Equivalent to KKT (4) asymptotically:
  if final g < 0  -> lambda -> 0     (inactive constraint)
  if final g = 0  -> lambda -> finite (active constraint)
  if final g > 0  -> lambda grows until S drives phi to reduce g

## Complete coupled algorithm

For opt step t = 0, 1, 2, ...:
  1. State solve   : evaluate f[phi_t], g[phi_t], sensitivities S_f, S_g
  2. Update lambda : lambda_{t+1} via D.27 (target-tracking sinh form)
  3. Total sens    : S = S_f + lambda_{t+1} S_g
  4. Update phi    : phi_{t+1} via D.11 (target-tracking phase-field + regulariser)
  5. Filter phi    : Helmholtz filter (App C) -> rho-tilde before next state solve
  6. Check KKT     : |S|_max < eps_KKT and |g| < eps_g -> optimum found, else loop

The dual evolution-equation pattern decouples primal-dual interaction
into two parallel pseudo-time integrations -- both implementable as
simple explicit Euler updates, no internal KKT linear-algebra solve.
"""


# ============================================================
# 8. Gangl-Sturm — nonlinear magnetostatics (3-paper bundle)
# ============================================================

GANGL_PART_I_SENSITIVITY = """
# Gangl-Sturm Part I — Nonlinear magnetostatics sensitivity analysis

Sturm K. (TU Wien) & Gangl P. (TU Graz).  AANMPDE Strobl 2019.
Companion to Part II (numerical realisation on 3D IPM-style electric
motor).

Builds on Sturm's parametrised-Lagrangian / averaged-adjoint framework
(Delfour-Sturm theorem) and the topological-derivative concept of
Sokolowski-Zochowski (1999) extended from linear elliptic to the
nonlinear curl-curl magnetostatics PDE.

## Setup — 3D quasilinear magnetostatic problem

Find u in V = H_0(D, curl) ∩ H(D, div=0) such that
    integral_D nu_Omega(x, |curl u|) curl u . curl v dx = <F, v>  for all v in V

where the source <F, v> = integral_Omega_c J.v dx + integral_Omega_m nu_m M . curl v dx
captures the coil current density J on conductor sub-region Omega_c
plus the permanent-magnet term M on Omega_m.

Reluctivity is state-dependent (BH curve on each material subdomain):
    nu_Omega(x, |curl u|) = chi_Omega * nu_2(|curl u|) + chi_{D\\Omega} * nu_1(|curl u|)

Reluctivities nu_1, nu_2: R+ -> R+ are assumed Lipschitz continuous,
strongly monotone (Zarantonello), and C^2 with bounded derivatives.

Cost is boundary tracking on Gamma_0 (air-gap circle in motor design):
    J(Omega) = integral_Gamma_0 |curl u . n - B_d^n|^2 dS

i.e. matching the normal component of B = curl u to a desired profile
B_d^n on a target surface (cogging-torque surrogate for IPM motor).

## Why linear-case formulas fail in nonlinear setting

Sokolowski-Zochowski 1999 TD was for LINEAR elliptic PDEs.  Naively
extending fails because:

(a) FRECHET DERIVATIVE IS THE TANGENT TENSOR, not the chord:
    b(x, v, w) = [nu(|v|) * I + (nu'(|v|)/|v|) v ⊗ v] * w
    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    Exactly the Newton-Raphson tangent in the forward solve.
    Linear-case adjoint is WRONG by O(1) factor wherever the B-H curve
    has saturated at z.

(b) Even with correct tangent adjoint, the cell-problem corrector W
    now solves a NONLINEAR exterior PDE (depending on the full B-H
    curve, not just its linearisation at curl u_0(z)).  No closed-form
    polarisation tensor exists.

## Averaged adjoint (the key fix)

Sturm 2015 averaged-adjoint approach.  Instead of classical adjoint
linearised at u_0, define the AVERAGED ADJOINT q_eps in X by:

    integral_0^1 d_u L(eps, s*u_eps + (1-s)*u_0, q_eps)(phi) ds = 0
                                                          for all phi in X

Replacing d_u L|_{u_0} with the s-average over the line segment
between u_0 and u_eps gives the EXACT identity:
    L(eps, u_eps, psi) = L(eps, u_0, q_eps)
                        (fundamental theorem of calculus)

Concretely, with F_eps(x, v) := nu_eps(x, |v|) v and the tangent
tensor b_eps(x, v, w) above, the averaged adjoint equation reads:

    integral_D integral_0^1 b_eps(x, curl(s u_eps + (1-s) u_0), curl v) ds .
                                     curl q_eps dx
    = -integral_Gamma_0 (curl(u_eps + u_0).n - 2 B_d^n) curl v.n dS

This averaging is THE substitute for Frechet-differentiating the
nonlinear constitutive law along an a-priori-unknown limit direction.
It preserves a clean MINIMAX identity even when nu(|.|) is only Lipschitz
(not C^1).

## Delfour-Sturm theorem

Given ell : R -> R with ell(0) = 0:
  (H1) perturbed-state and averaged-adjoint sets are non-empty for
       all eps in [0, tau]
  (H2) d_ell L(0, u_0, q_0) := lim (L(eps, u_0, q_0) - L(0, u_0, q_0))/ell(eps) exists
  (H3) R := lim (L(eps, u_0, q_eps) - L(eps, u_0, q_0))/ell(eps) exists

Then: d_ell g(0) = d_ell L(0, u_0, q_0) + R

For the magnetostatic TD application:  ell(eps) = |omega_eps| ~ eps^d
(d = 3 here).  First term is EASY (explicit change-of-variables limit
at insertion point z).  R is the corrector — vanishes only in linear
case.

## Topological derivative formula

At insertion point z in D \\ Omega:

    TD(z) = [nu_2(|curl u_0(z)|) - nu_1(|curl u_0(z)|)]
              * curl u_0(z) . curl q_0(z)
            + R(u_0, q_0)

First term = easy "polarisation" part.  After eliminating Q by testing
the W-equation with Q, R becomes:

    R = (1/|omega|) integral_{R^3} [F_omega(x, curl W + curl u_0(z))
                                     - F_omega(x, curl u_0(z))
                                     - d_y F_omega(x, curl u_0(z)) curl W] .
                                     curl q_0(z) dx
        + (1/|omega|) integral_omega [d_y F_2 - d_y F_1](curl u_0(z))
                                      curl W . curl q_0(z) dx

W in BL(R^3) (Beppo-Levi space) solves the EXTERIOR NONLINEAR CELL
PROBLEM:

    integral_{R^3} [F_omega(x, curl W + curl w_0(z))
                    - F_omega(x, curl w_0(z))] . curl v dx
    = -integral_omega [nu_1(|curl w_0(z)|) - nu_2(|curl w_0(z)|)]
                       curl w_0(z) . curl v dx

In the LINEAR limit (nu_i constant), the correctors W and Q reduce to
the Polya-Szego polarisation tensor; R collapses to the classical
Sokolowski-Zochowski polarisation-tensor formula.  The NONLINEAR R
term is the HEADLINE contribution of this paper.

## Perturbation a-priori estimates

Lemma 1:  ||u_eps - u_0||_L2 + ||curl(u_eps - u_0)||_L2 <= C eps^(d/2)
Lemma 2:  ||q_eps - q_0||_H(D,curl) <= C eps^(d/2)

The eps^(d/2) rate comes from the inclusion having volume eps^d and
the energy estimate scaling like sqrt(perturbation volume) — standard
for elliptic problems with small inclusions.

## Practical implementation notes

  (1) The full corrector pair (W, Q) on R^3 is "numerically infeasible"
      as written — cannot mesh R^3.  Q-elimination via fundamental-
      theorem-of-calculus testing converts R into an integral over
      R^3 that, for practical inclusion shapes (sphere, ellipsoid), can
      be evaluated semi-analytically OR via a TRUNCATED-DOMAIN FE solve
      on a ball B_R(0) with R >> diam(omega).  KELVIN TRANSFORMATION
      (see radia_mcp.radia_ngsolve.kelvin) is the natural lab tool for
      this truncation.

  (2) The tangent tensor b_eps = nu*I + (nu'/|v|) v⊗v is the SAME tensor
      that appears in NEWTON-RAPHSON linearisation of the forward solve.
      A code that already does Newton on the magnetostatic state has
      the operator for the averaged adjoint essentially FOR FREE —
      just plug in s*u_eps + (1-s)*u_0 instead of u_0.

  (3) Helmholtz decomposition (slide 14) suggests gauge-fixing with
      H_0(D, curl) ∩ H(D, div=0) — in practice realised either via
      Nedelec edge elements + tree-cotree gauge, or via grad-div
      regularisation (radia_mcp.radia_ngsolve.basis_functions).

  (4) The TD field z -> TD(z) is computed on a DESIGN GRID covering
      D\\Omega and D ∩ Omega.  Descent inserts material where TD < 0
      and removes where TD > 0, in a level-set or characteristic-
      function update.

## Lab relevance

Directly relevant to motor / accelerator-magnet topology optimisation.
For the lab's production stack:
  (i) Averaged-adjoint equation can be assembled by calling
      radia.MatSatIsoTab (BH curve) tangent-reluctivity routines along
      the line segment between u_0 and u_eps.
  (ii) Exterior cell problem for W can be solved on a truncated NGSolve
       mesh with Kelvin transformation (matching the existing
       accelerator-magnet kelvin pipeline).
  (iii) TD field can drive a level-set update over a fixed Cubit hex
       mesh in the iron-yoke design domain.
"""


GANGL_PART_II_NUMERICS = """
# Gangl-Sturm Part II — NGSolve numerical realisation

Gangl-Sturm-Schoeberl AANMPDE Strobl 2019 (51-slide deck).
Demonstrates shape + topology optimisation of a 3D electric motor on
top of NGSolve + ngsxfem.

## NGSolve implementation overview

State PDE — quasilinear curl-curl with eps-regulariser:

    integral_D nu_Omega(x, |curl u|) curl u . curl v + eps * u . v dx = <F, v>

on design domain D, with u in V := H_0(D, curl), discretised using
NEDELEC EDGE ELEMENTS (NGSolve `HCurl` FE space).

KEY POINT: the eps * u . v REGULARISER IS ESSENTIAL.  Without it the
H(curl) bilinear form has a huge |grad H1| nullspace and is singular.
Use eps ~ 1e-6 (in SI units).

```python
from ngsolve import *

mesh = Mesh(...)
fes = HCurl(mesh, order=2)
u, v = fes.TnT()

# Reluctivity as a CoefficientFunction (BH-curve spline)
def nu_simp(rho, curl_u_norm):
    # ...interpolate the BH curve
    pass

a = BilinearForm(fes, symmetric=True)
a += nu_simp(rho, Norm(curl(u))) * curl(u) * curl(v) * dx
a += eps * u * v * dx                # regulariser — DO NOT OMIT
```

## Newton-Raphson primal solve

The reluctivity nu_i(s) * s being Lipschitz + strongly monotone makes
the global energy CONVEX.  Standard Newton with Armijo backtracking
on energy functional:

    W(u) = integral_D Nu(|curl u|) - <F, u>

Tangent bilinear form (= the same b(x, v, w) of Part I):

    K'(u)[du, v] = integral_D [nu(|curl u|) curl du . curl v
                                + nu'(|curl u|)/|curl u| (curl u . curl du)
                                                          (curl u . curl v)]
                                + eps du . v dx

SPD on H(curl) modulo gradients (the eps u.v term restores
invertibility).

NGSolve idiom:
```python
a = BilinearForm(fes)
a += nu(Norm(curl(u))) * curl(u) * curl(v) * dx
a += eps * u * v * dx
a.AssembleLinearization(gfu.vec)     # assembles tangent operator
```

Typical convergence: 4-8 Newton iterations to ||R|| / ||F|| < 1e-8 for
lab-class BH curves.  Damped Newton (line search on residual) required
near saturation knees where nu'(s) varies by orders of magnitude.

## Adjoint assembly in NGSolve

Adjoint p solves a LINEAR problem with the FROZEN Newton tangent K'(u)
and cost-gradient RHS:

    Find p in V s.t. K'(u)[w, p] = -dJ/du[w] for all w in V

SAME SPD-on-quotient operator as the last Newton step -> one extra
back-substitution or CG solve (reuse the factorisation).

NGSolve symbolic shape-derivative engine (slide 31):

```python
# Symbolic Lagrangian
G = Cost(gfu) + Equation(gfu, gfp)

# Forward + adjoint share the same physics-BilinearForm
a += Equation(u, v)

# Shape-derivative LinearForm on the velocity space
dJOmega = LinearForm(VEC)
dJOmega += G.Derive(xi, div(X))           # determinant perturbation
dJOmega += G.Derive(F, Grad(X))           # Jacobian perturbation
dJOmega += G.Derive(x, X[0]) + G.Derive(y, X[1])   # coordinate perturbations
```

G.Derive(symbol, direction) applies the chain rule through
CoefficientFunction expressions automatically — user never assembles
the Eulerian dJ formula by hand.

For the H(curl) version (slide 33) the COVARIANT TRANSFORMATION
    (u o T_t)(x) = (DT_t)^(-T) u-hat(x)
    curl(u o T_t) = (1/det DT_t) DT_t curl u-hat
is baked into Equation(u, v) via F^(-T) and (1/xi)*F*curl, and G.Derive
differentiates symbolically through these.

PITFALL: USING A LAGRANGIAN PULL-BACK FOR VECTOR POTENTIALS GIVES A
WRONG SHAPE DERIVATIVE.  The covariant transformation is mandatory for
H(curl).

## Level-set evolution (Amstutz-Andra fixed-point)

NOT a Hamilton-Jacobi PDE solver.  Convention:
    psi(x) > 0 <-> material 1,  psi(x) < 0 <-> material 2

Generalised topological derivative g_psi(x) defined as +T^(1->2)(x)
inside Omega_1 and -T^(2->1)(x) inside Omega_2 — one signed scalar
field over the whole design domain.

OPTIMALITY LEMMA: at a local optimum, psi is parallel to g_psi (i.e.
psi = g_psi up to positive normalisation).

UPDATE:
    psi_{k+1} = (1 - s) psi_k + s * g_psi

with step size s in [0, 1] chosen as large as possible such that
J(psi_{k+1}) < J(psi_k) (line search / backtracking).  No CFL
condition, no reinitialisation PDE — the L^2 normalisation of g_psi
handles stability.

In NGSolve, psi and g_psi live as L^2 GridFunctions on the background
tetrahedral mesh; update is one `.vec.data` assignment.

## TD hole nucleation in 3D

For the LINEAR-transmission model problem (slide 42):

    T^(Omega -> D\\Omega)(x) = 2 alpha_0 (alpha_1 - alpha_0)/(alpha_1 + alpha_0)
                                * pi * grad u(x) . grad p(x)   for x in Omega
    T^(D\\Omega -> Omega)(x) = 2 alpha_1 (alpha_0 - alpha_1)/(alpha_0 + alpha_1)
                                * pi * grad u(x) . grad p(x)   for x in D\\Omega

i.e. inner product of state and adjoint gradients scaled by Polya-Szego-
like ellipsoid polarisation factor.

T(x_0) < 0 -> J(Omega_eps) < J(Omega) for small hole radius eps, so
the NEGATIVE regions of the TD field are candidates for material swap.

Nucleation in 3D: evaluate g_psi on every mesh node (one NGSolve
grad(gfu) * grad(gfp) CF evaluated at vertices) -> feed into the
fixed-point update.  No explicit small-ball insertion, no remeshing,
no XFEM cut-cell creation.

The material discontinuity at the zero level set of psi is resolved by
NGSXFEM (Lehrenfeld), not by remeshing.

## Numerical pitfalls (Strobl 2019 slides + standard ecosystem)

  (1) NEVER drop the eps * u . v regulariser in the H(curl) state form.
      Without it the curl-curl operator has a |grad H1| nullspace and
      linearised Newton system is singular.  But eps too large biases
      the solution — keep eps ~ O(mu_0 * 1e-6) in SI.

  (2) Covariant pull-back (DT_t)^(-T) for H(curl) is NOT optional.
      Lagrangian pull-back for vector potential yields WRONG shape
      derivative.

  (3) G.Derive symbolic differentiation only sees what is in the
      BilinearForm.  If you compute B = curl u outside the form, the
      shape derivative will miss boundary terms.  Define cost J(u)
      inside Cost(u) as a CF expression of grad/curl gfu.

  (4) Mesh.Curve(p) for electric-motor air-gap target boundary Omega_g.
      Without it, surface normal n in |curl u . n - B_d^n|^2 is
      piecewise constant -> gradient field becomes facet-discretisation
      dominated rather than physics dominated.

  (5) Nitsche stabilisation: pick lambda = O(alpha_max) * 10-100.
      Too small -> cut-FEM bilinear form loses coercivity (Newton
      diverges).  Too large -> conditioning blows up.

  (6) Cut-cells from ngsxfem produce ill-shaped sub-tetrahedra near the
      zero level set.  Use CG with diagonal / BDDC preconditioner, not
      direct factorisation, when level set sits very close to a vertex.

  (7) Fixed-point step s in psi_{k+1} = (1-s) psi_k + s g_psi MUST be
      picked by backtracking on J.  Fixed s diverges as soon as g_psi
      changes sign.  Start with s = 1 and halve on failure (Armijo).
"""


GANGL_IPM_MOTOR_CASE = """
# Gangl IPM motor case study (SIAM SISC 2015)

Gangl-Langer-Laurain-Meftahi-Sturm (2015) SIAM J. Sci. Comput. 37,
B1002-B1025.  The canonical demonstration the Gangl group built around
the Part-I/II sensitivity theory.

## Problem setup

Interior Permanent Magnet (IPM) brushless electric motor, modelled as
a 2D cross-section (planar magnetostatic, justified for switched-
reluctance motor 2D vs 3D comparison [25, 42]).

Geometry: outer stator with ferromagnetic iron core, inner rotor
containing PMs, separated by thin air gap.  Eight PM regions Omega_ma
(M = (Mx, My) prescribed magnetisation) embedded in rotor iron
Omega_f_ref.

Coil regions Omega_c on inner stator face but modelled as air because
(a) coils switched OFF in this scenario (Ji = 0), (b) copper permeability
~ air.

Whole domain D = smooth bounded 2D disk with Dirichlet BC u = 0 on
outer boundary -> B.n = 0, no flux leaves D.

Governing PDE — nonlinear magnetostatics:
    -div(nu(|grad u|) grad u) = f
    u = z-component of magnetic vector potential A so B = curl((0,0,u))
    f = f0 + div(M) — PM magnetisation source

Reluctivity nu-hat(|B|) supplied as measured BH curve (monotone
decreasing from ~10^7 down to roughly 10^2-10^3 across |B| = 0 to 3 T),
C^1 spline-interpolated per Juttler-Pechstein.

## Objective — air-gap Br tracking (NOT direct cogging-torque)

NOT a directly-evaluated torque integral.  Tracking-type cost on the
RADIAL COMPONENT Br of the magnetic flux density along a circle Gamma_0
inside the air gap:

    J(Omega) = integral_Gamma_0 |Br(u_Omega) - Bd|^2 ds

Desired profile = pure sinusoid Bd(theta) = 0.5 * sin(4 theta) — 4-pole
pattern matching the 8 magnet regions (4 pole pairs).

Matching this single-harmonic target MINIMISES the Total Harmonic
Distortion (THD) of the air-gap flux, which IN TURN smooths rotation
and reduces cogging torque.

So 'cogging-torque minimisation' is solved INDIRECTLY through a THD /
L2-tracking surrogate.  THE PAPER DOES NOT compute torque vs rotor-
angle and DOES NOT sweep the rotor through a cogging period.

Only ONE rotor position evaluated (the coils-off equilibrium).
Rotation smoothness comes from the assumption that a sinusoidal
air-gap Br at this position implies low harmonic content at
neighbouring positions.

Initial L2 distance: 1.3033e-3
Final L2 distance:   0.94997e-3

## Design region

Design subregion Omega is a STRICT SUBSET of the rotor iron core
Omega_f_ref.  It is NOT the magnet region (PMs are FIXED: paper
enforces T_t(Omega_ma) = Omega_ma in velocity transformation), and NOT
the whole rotor.

Specifically, Omega is partitioned into 8 connected components
Omega_ref_k (k=1..8) — one per pole pair, each a small sub-patch of
rotor iron immediately ADJACENT to the magnet pocket.

Velocity-field constraints:
  V.n = 0 on south side Gamma_S of each Omega_ref_k (abuts the magnet)
  V.n <= 0 on N/E/W sides (design region cannot grow beyond Omega_ref)

Total: 8 independent design patches with shared 4-fold symmetry.

## Numerical results — IMPORTANT corrections to lab knowledge

After 35 iterations of gradient-descent shape update + line search:

    Cost reduction: **27%** (1.3033e-3 -> 0.94997e-3)
    NOT 87% as previously recorded in some lab knowledge bases.

Br peak amplitude reduced from ~+-0.45 T to ~+-0.40 T.  High-frequency
ripples near each pole-transition substantially flattened.

Wall time: ~26 minutes on a SINGLE CORE of a laptop.
Mesh: 44 810 DOFs, 89 454 elements, 53 488 of which are inside design
region.

COGGING TORQUE ITSELF IS NOT REPORTED NUMERICALLY (no mNm value in the
paper) — only the surrogate L2 cost.

## Lessons learned

  (1) THE LAGRANGIAN-WITHOUT-MATERIAL-DERIVATIVE TRICK is the paper's
      core contribution: for nonlinear PDE constraints, computing the
      material derivative u-dot is intractable.  The shape-Lagrangian
      formulation (Sturm 2015, Laurain-Sturm 2016) yields the shape
      derivative dJ(Omega; V) directly from state u, adjoint p, velocity V,
      WITHOUT u-dot — Theorem 5 / equation (27).

  (2) VOLUME (DOMAIN) EXPRESSION OF SHAPE DERIVATIVE is used, NOT
      boundary expression.  Volume form is more robust on coarse meshes
      and integrates the descent direction over a tubular neighbourhood.

  (3) BACKTRACKING LINE SEARCH on step size tau is essential — without
      it, large initial steps cross the iron/air interface in physically
      unrealistic ways.  Termination: tau too small for ANY element to
      switch state.

  (4) BINARY ELEMENT SWITCHING (centroid-in-polygon test) avoids mesh
      distortion at the cost of interface jaggedness.  XFEM /
      immersed-FEM / Nitsche-DG would give sharper interface.

  (5) TOPOLOGICAL DERIVATIVE was studied for the same motor earlier
      [11] but produces JAGGED boundaries.  Shape derivative gives
      smooth interfaces — the explicit motivation here.

  (6) WEAK FORMULATION OF div(M) source: integration by parts moves
      divergence onto v, so f = Ji + div(M) becomes integral Ji*v + M.grad(v)
      (eq. 40).  Standard PM source coupling in nonlinear magnetostatics.

  (7) 27% REDUCTION IS NOT DRAMATIC — binary switching limits achievable
      resolution.  Paper's emphasis is on RIGOR of sensitivity analysis
      for nonlinear problems, NOT on achieving a record cogging-torque
      reduction.

## Mapping to the lab NGSolve + Radia stack

(a) STATE SOLVE:
    -div(nu(|grad u|) grad u) = Ji + div(M) with measured BH curve maps
    directly to NGSolve H1(mesh, dirichlet='outer') + Newton-coupled
    GridFunction.  BH curve goes through radia.MatSatIsoTab or a C^1 spline.
    2D approximation (Bz alone) = FEMM/ELF 2D-planar regime.  In NGSolve
    this is the SCALAR A_z formulation on a triangular mesh.

(b) ADJOINT SOLVE:
    Equation (22) is a LINEAR elliptic problem with frozen coefficient:
        A(grad u) = nu(|grad u|) I_2 + 2 d_zeta(nu)(|grad u|^2) grad u (X) grad u
    (anisotropic 2x2 stiffness from BH-curve nonlinearity) and surface
    load 2*(Br(u) - Bd) on Gamma_0.  Assemble as BilinearForm with
    NGSolve CoefficientFunction.  radia.RadiaField is NOT needed since
    this is bounded-domain (motor in fixed disk, no Kelvin / open
    boundary).

(c) SHAPE DERIVATIVE:
    Eq. (27) is a sum of domain integrals involving grad u, grad p,
    div(V), DV, BH-curve derivatives.  Assemble as another BilinearForm
    on a velocity space (NGSolve VectorH1 over rotor sub-domain), then
    solve H1-Riesz problem (eq. 41-42) with alpha-weighted bilinear form
        b(V, W) = integral_Drot (alpha DV:DW + V.W) dx
    alpha = {1 in Omega, 10 in tube, 100 elsewhere}.

(d) UPDATER:
    Paper uses pure GRADIENT-DESCENT + BACKTRACKING LINE SEARCH on
    step size tau.  NOT OC (Optimality Criteria).  NOT MMA.  No volume
    constraint requires MMA's multiplier loop.  For our stack: replace
    binary element-switch with density (rho in [0, 1]) field + SIMP
    penalisation -> can use OC, or radia_mcp.optuna black-box if FEM
    solve wrapped as objective.

(e) ROTATION DISCRETISATION:
    Paper uses ZERO — only one rotor position.  If real cogging torque
    dT/dtheta minimisation is wanted, sweep N_theta ~ 8-16 positions
    per pole-pitch and sum |T(theta_k)|^2 or peak-to-peak T.  The
    shape-Lagrangian framework still applies (sum of independent terms).

(f) MESH STRATEGY:
    Paper uses ONE GLOBAL MESH with 8 design sub-patches.  In NGSolve:
    OCC -> Mesh() for whole motor cross-section, mark 8 design patches
    with Material('design_k'), let Mesh.Curve(2) handle curved magnet
    boundaries.

(g) PARALLELISM:
    PARDISO direct solver.  Substitute NGSolve's PARDISO interface
    (MKL-shipped).  For larger 3D problems switch to Compact HX
    preconditioner per lab POLICY on HCurl/H1 problems.

(h) RADIA ROLE:
    Marginal — this is a bounded-domain FEM problem with no open-
    boundary need, so radia.RadiaField / MMM / MSC do not apply.
"""


# ============================================================
# 9. NGSolve recipes for the full topology-optimization workflow
# ============================================================

LAB_NGSOLVE_RECIPES = """
# Lab NGSolve recipes for TopOpt

Practical recipes for implementing the textbook framework in NGSolve.

## Recipe A — Linear self-adjoint problem (heat / elasticity / Stokes)

The easy case.  Per iteration:
  1. Solve primal physics with material_property(rho) = eps + (1-eps) rho^p
  2. Compute sensitivity from primal alone (self-adjoint): one-line CF
     `S(x) = -p (1-eps) rho^(p-1) (primal_quantity^2)`
  3. Update phi via Helmholtz filter + reaction-diffusion evolution
  4. Loop

```python
from ngsolve import *
import numpy as np

mesh = Mesh(...)
fes_state    = H1(mesh, order=1, dirichlet='T_bdy')
fes_rho      = L2(mesh, order=0)            # piecewise constant
fes_phi      = H1(mesh, order=1)            # design variable space
fes_filtered = H1(mesh, order=1)            # filtered density

gfu   = GridFunction(fes_state)
rho   = GridFunction(fes_rho)
phi   = GridFunction(fes_phi)
rho_t = GridFunction(fes_filtered)

# Initial uniform density
rho.vec[:] = volfrac

# SIMP CoefficientFunction
p_SIMP, eps = 3, 1e-5
kappa_simp = eps + (1 - eps) * rho_t**p_SIMP

# 1. State solve
u, v = fes_state.TnT()
a = BilinearForm(fes_state, symmetric=True)
a += kappa_simp * grad(u) * grad(v) * dx
f = LinearForm(fes_state)
f += Q * v * dx
for opt_iter in range(MAX_ITER):
    a.Assemble()
    f.Assemble()
    gfu.vec.data = a.mat.Inverse() * f.vec

    # 2. Self-adjoint sensitivity (heat)
    sens_cf = -p_SIMP * (1 - eps) * rho_t**(p_SIMP - 1) * grad(gfu)*grad(gfu)
    sens = GridFunction(fes_rho)
    sens.Set(sens_cf)

    # 3. Helmholtz filter (App C)
    u_t, v_t = fes_filtered.TnT()
    a_filter = BilinearForm(fes_filtered, symmetric=True)
    a_filter += R_filter**2 * grad(u_t) * grad(v_t) * dx
    a_filter += u_t * v_t * dx
    f_filter = LinearForm(fes_filtered)
    f_filter += phi * v_t * dx
    a_filter.Assemble()
    f_filter.Assemble()
    rho_t.vec.data = a_filter.mat.Inverse() * f_filter.vec

    # 4. RD evolution for phi (one-step explicit Euler)
    # m d phi/dt = C_diff Lap(phi) - sens(x)
    # ...one Helmholtz-like solve per opt iter (Yamada semi-implicit)
```

## Recipe B — Nonlinear non-self-adjoint (magnetostatics + general NS)

Per iteration:
  1. PRIMAL SOLVE — Newton-Raphson (4-8 inner steps typical)
  2. ADJOINT ASSEMBLY with same Jacobian — reuse K'(u)
  3. ADJOINT SOLVE — one linear back-substitution
  4. SENSITIVITY = primal x adjoint inner product x material-slope
  5. DENSITY FILTER (Helmholtz, App C)
  6. PHI UPDATE via RD evolution PDE
  7. LAMBDA UPDATE via App D dual evolution
  8. Loop

```python
fes = HCurl(mesh, order=2)   # Nedelec edge elements
u, v = fes.TnT()

eps_reg = 1e-6   # H(curl) regulariser — DO NOT OMIT

# Define the nonlinear reluctivity (BH curve)
def nu_BH(rho, B_norm):
    # SIMP-style interpolation of BH curve
    nu_void = 1 / mu_0
    nu_iron = 1 / (mu_0 * mu_r_iron(B_norm))   # nonlinear
    return nu_void + (nu_iron - nu_void) * rho**p_SIMP

# Bilinear form (nonlinear)
a = BilinearForm(fes, symmetric=True)
a += nu_BH(rho_t, Norm(curl(u))) * curl(u) * curl(v) * dx
a += eps_reg * u * v * dx

# Newton-Raphson loop on the primal
gfu = GridFunction(fes)
gfu.vec[:] = 0
for newton_step in range(MAX_NEWTON):
    a.AssembleLinearization(gfu.vec)
    res = a.mat.CreateRowVector()
    a.Apply(gfu.vec, res)
    res.data -= f.vec
    if Norm(res) < tol_newton:
        break
    inv = a.mat.Inverse(fes.FreeDofs())
    du = -inv * res
    # Armijo line search on energy
    alpha = 1.0
    while True:
        gfu_trial = gfu + alpha * du
        if energy(gfu_trial) < energy(gfu) + 1e-4 * alpha * res.dot(du):
            break
        alpha *= 0.5
    gfu.vec.data += alpha * du

# Adjoint: linear with frozen Newton-tangent (already assembled)
gfp = GridFunction(fes)
f_adj = LinearForm(fes)
f_adj += dJ_du * v * dx   # cost-gradient RHS
f_adj.Assemble()
gfp.vec.data = a.mat.Inverse(fes.FreeDofs()) * f_adj.vec

# Sensitivity (averaged-adjoint for nonlinear case, see Gangl Part I)
# S = -(d nu / d rho) * (curl u . curl p) + lambda dg/drho
# ...
```

## Recipe D — Hybrid TD + level-set (Gangl line)

When the initial design has poor topology (no obvious starting candidate)
and you need to NUCLEATE NEW HOLES:

  1. TD FIELD evaluation: compute g_psi(x) = curl u . curl q everywhere
     on design grid (Gangl §3 polarisation expression for linear case;
     for nonlinear see Gangl Part I averaged-adjoint corrector).
  2. INSERT material wherever g_psi < 0 (or remove wherever g_psi > 0).
  3. LEVEL-SET UPDATE via Amstutz-Andra fixed-point:
        psi_{k+1} = (1 - s) psi_k + s g_psi
     Armijo backtracking on s.
  4. Optionally post-process with shape-derivative descent (smoother
     boundary).

Pure level-set methods cannot create new holes; pure TD methods give
jagged boundaries.  The hybrid is the lab default.

## Connection to other radia_ngsolve modules

  - `kelvin` — truncated cell-problem solve for Gangl-Sturm W/Q
                exterior corrector
  - `mmm_core` — MMM tangent reluctivity for averaged adjoint
  - `basis_functions` — Nedelec HCurl edge elements (state space)
  - `esim` — SIBC / nonlinear surface impedance (related sensitivity
              framework for fluid-conductor TopOpt extensions)
  - `sparsesolv` — PARDISO + Compact HX preconditioner for large 3D
"""


# ============================================================
# Dispatcher
# ============================================================

TOPICS = {
    "overview":            FRAMEWORK_OVERVIEW,
    "framework":           FRAMEWORK_OVERVIEW,    # alias
    "ch1_2":               TEXTBOOK_CH1_2,
    "foundations":         TEXTBOOK_CH1_2,        # alias
    "ch3":                 TEXTBOOK_CH3_EVOLUTION_EQ,
    "evolution_equation":  TEXTBOOK_CH3_EVOLUTION_EQ,  # alias
    "ch4_5":               TEXTBOOK_CH4_5_HEAT_ELASTICITY,
    "heat_elasticity":     TEXTBOOK_CH4_5_HEAT_ELASTICITY,  # alias
    "ch6_8":               TEXTBOOK_CH6_8_FLUID,
    "fluid":               TEXTBOOK_CH6_8_FLUID,  # alias
    "ch9":                 TEXTBOOK_CH9_LBM,
    "lbm":                 TEXTBOOK_CH9_LBM,      # alias
    "appendix_a":          TEXTBOOK_APPENDIX_A_OBJECTIVES,
    "objectives":          TEXTBOOK_APPENDIX_A_OBJECTIVES,   # alias
    "appendix_c":          TEXTBOOK_APPENDIX_C_FILTER,
    "helmholtz_filter":    TEXTBOOK_APPENDIX_C_FILTER,       # alias
    "appendix_d":          TEXTBOOK_APPENDIX_D_KKT,
    "kkt":                 TEXTBOOK_APPENDIX_D_KKT,          # alias
    "gangl_part1":         GANGL_PART_I_SENSITIVITY,
    "gangl_sensitivity":   GANGL_PART_I_SENSITIVITY,         # alias
    "averaged_adjoint":    GANGL_PART_I_SENSITIVITY,         # alias
    "gangl_part2":         GANGL_PART_II_NUMERICS,
    "ngsolve_implementation": GANGL_PART_II_NUMERICS,         # alias
    "gangl_motor":         GANGL_IPM_MOTOR_CASE,
    "ipm_motor":           GANGL_IPM_MOTOR_CASE,             # alias
    "ngsolve_recipes":     LAB_NGSOLVE_RECIPES,
    "recipes":             LAB_NGSOLVE_RECIPES,              # alias
}


def get_topology_optimization_documentation(topic: str = "all") -> str:
    """Return topology-optimization documentation for the requested topic.

    Args:
        topic: one of TOPICS keys, or 'all' for full document.

    Available canonical topics (aliases also accepted):
        overview / framework        — 6-section template + 6 dial numbers
        ch1_2 / foundations         — history + 4 method families + OC + SCP
        ch3 / evolution_equation    — reaction-diffusion engine (Kondoh / Yamada)
        ch4_5 / heat_elasticity     — self-adjoint compliance benchmarks
        ch6_8 / fluid               — Stokes + laminar NS + conjugate HT (Toyota CRDL)
        ch9 / lbm                   — lattice Boltzmann TopOpt (Yaji, world-first)
        appendix_a / objectives     — 5-physics objective cookbook
        appendix_c / helmholtz_filter — canonical density-filter regulariser
        appendix_d / kkt            — KKT dual-evolution scheme
        gangl_part1 / gangl_sensitivity / averaged_adjoint
                                    — nonlinear-magnetostatics sensitivity analysis
        gangl_part2 / ngsolve_implementation
                                    — Newton + adjoint + level-set in NGSolve
        gangl_motor / ipm_motor     — IPM motor case study (27% reduction)
        ngsolve_recipes / recipes   — practical NGSolve code patterns
    """
    key = topic.strip().lower()
    if key == "all":
        return "\n\n---\n\n".join([
            FRAMEWORK_OVERVIEW,
            TEXTBOOK_CH1_2,
            TEXTBOOK_CH3_EVOLUTION_EQ,
            TEXTBOOK_CH4_5_HEAT_ELASTICITY,
            TEXTBOOK_CH6_8_FLUID,
            TEXTBOOK_CH9_LBM,
            TEXTBOOK_APPENDIX_A_OBJECTIVES,
            TEXTBOOK_APPENDIX_C_FILTER,
            TEXTBOOK_APPENDIX_D_KKT,
            GANGL_PART_I_SENSITIVITY,
            GANGL_PART_II_NUMERICS,
            GANGL_IPM_MOTOR_CASE,
            LAB_NGSOLVE_RECIPES,
        ])
    if key not in TOPICS:
        return (f"Unknown topic '{topic}'. Available: "
                f"{', '.join(sorted(set(TOPICS)))}, or 'all'.")
    return TOPICS[key]
