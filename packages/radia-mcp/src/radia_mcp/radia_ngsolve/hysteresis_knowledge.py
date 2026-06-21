# -*- coding: utf-8 -*-
r"""Magnetic-hysteresis modelling knowledge for the radia-ngsolve core.

Synthesised from the magnetic-properties / hysteresis-model literature -- the model landscape, the
authoritative play-model formulation (Matsuo / Bobbio), the lab's own B-input congruency results, shape
identification, and the FE / iron-loss path.  Companion to the implementation in
:mod:`radia_mcp.radia_ngsolve.hysteresis` (the B-input vector Play operator) and the energy-based
lineage in :mod:`radia_mcp.motor.henrotte_lineage_knowledge`.
"""
from __future__ import annotations


LANDSCAPE = """\
## The hysteresis-model landscape

Rate-independent ("DC") scalar/vector hysteresis-model families, with the lab's practical view:

| Family | Idea | Pros / cons (for motor FE) |
|--------|------|----------------------------|
| Jiles-Atherton | anhysteretic + pinning ODE, ~5 params | few params, but limited flexibility, scalar only, minor-loop trouble |
| Preisach | weighted integral of rectangular hysterons (Everett fn) | very expressive; congruency + wiping-out built in; costly |
| **Play** | sum of play hysterons + shape functions | **equivalent to static Preisach but cheap; the practical FE choice** |
| Stop | sum of stop hysterons (the dual of play) | H-input dual; play/stop are complementary (Matsuo 2003) |
| Energy-based | friction-pinning, energy-consistent, vector | thermodynamic, vector-intrinsic (Henrotte 2006 / Jacques 2018) |
| Chua / Chan / Bouc-Wen / LLG | circuit / phenomenological / micromagnetic | niche; not the iron-loss workhorse |

The lab's working choice for electric-machine iron loss is the **play model** -- mathematically
equivalent to the static Preisach model (so equally expressive) but far cheaper, and widely used in
commercial motor-design hysteresis modelling.  Energy-based (Henrotte/Jacques) is the same friction
idea recast to be intrinsically vector and energy-consistent; the play and energy models coincide in
the unified Radia Type 5/6 formulation (energy ``chi`` == play threshold ``eta``).
"""

PLAY_FORMULATION = """\
## Play model -- authoritative formulation (Matsuo-Shimasaki 2003; Bobbio 1997)

Input x, output y (x,y = (B,H) for B-input, (H,B) for H-input):

    y(x) = sum_{k=1}^{2N} f_k( p_k(x) )                                   (1)

with the PLAY HYSTERON (width zeta_k), a clamp / backlash operator:

    p_k(x) = max( x - zeta_k , min( x + zeta_k , p_k^0 ) )                (2)

(p_k^0 = previous-step value; symmetric about x because B-H behaviour is odd).  Equal-interval widths
for deterministic Everett identification:

    zeta_k = (k - 1/2) * Xmax / N                                         (3)

The SHAPE FUNCTION f_k is piecewise-linear in p_k; its slope m_{k,j} is a RELUCTIVITY for B-input
(a permeability for H-input).  So a B-input play cell is "reluctivity a_k acting on the played flux".

radia_ngsolve.hysteresis implements exactly this (B-input): the play update
``p_k <- B - eta_k (B - p_k)/max(eta_k,|B-p_k|)`` is algebraically (2), and linear cells ``f_k=a_k r``
give the verification hook with analytic loop area ``4 a_k eta_k (Bm-eta_k)``.
"""

BINPUT_CONGRUENCY = """\
## B-input vs H-input + congruency (Hane-Sugahara-Morishita-Ahagon, IEEE TMAG 2026)

CONGRUENCY = all minor loops have identical SHAPE regardless of operating point, within the same input
range.  A B-input congruency-based model (Preisach or play) is congruent on the **H-axis**; an H-input
model is congruent on the **B-axis** (they cannot both hold).

LAB MEASUREMENT (non-oriented silicon steel 35A360, ring core, minor loops under several excitations):
- The material has **H-axis CONGRUENCY** (winding-current H is external to the material) and
  **B-axis INCONGRUENCY** (B is the magnetisation-influenced average).
- Therefore the **B-input play model reproduces the measured minor loops -- shapes, operating points,
  areas, and losses** -- while the **H-input play model fails** (it forces a B-axis congruency that the
  material does not have: wrong minor-loop shapes, even overshoot at loop closure).
- Measured vs B-input-play vs H-input-play minor-loop losses (J/m^3) agree well only for B-input
  (e.g. 17.2 / 15.6 vs 20.2; 4.5 / 3.9 vs 10.0).

CONSEQUENCE FOR radia-ngsolve:  the AGE / A_z core has B = curl A as primary, so the B-input play model
is BOTH the natural FE fit (A-formulation <-> B-input; T-Omega <-> H-input) AND the more accurate model
for real NOES.  This is why ``radia_ngsolve.hysteresis`` is B-input.  (The lab's congruency test is the
property ``test_hysteresis_play.test_minor_loop_congruency`` checks: H-axis congruency of the B-input
operator.)
"""

IDENTIFICATION = """\
## Shape-function identification (from measured DC loops)

Two routes in the lab, both deterministic:

1. EVERETT FUNCTION (Preisach-equivalent, the canonical route).  From the ascending Y+ and descending
   Y- branches of measured symmetric loops of amplitude a, the Everett function
       E(alpha,beta) = [ (Y-(a,alpha)-Y-(a,beta)) - (Y+(a,alpha)-Y+(a,beta)) ] / 2
   gives the piecewise-linear shape-function slopes m_{k,j} directly (Hane 2026 Eq. 5-7).  Needs ~30
   B-input loops (vs ~80 H-input) at equal Bmax intervals -> B-input also needs FEWER measurements.

2. STAIRCASE / s(i,j) RECURSION (ShapeFunction.m; lab ``radia.hysteresis_io.build_shape_functions``):
   the shape-function increments are differences of measured loop points across nested loops
   (e.g. s(0,2)=A(1,2)-B(1,1)+s(0,1), s(1,1)=B(1,1)-D(0,0)).  Equivalent bookkeeping of the same data.

LOOP INTERPOLATION (Morishita-Sugahara-Hane-Ahagon, Compumag 2025): fit each measured loop with an odd
5th/even-4th RATIONAL function, normalise by (Bmax,Hmax), interpolate the coefficients vs Bmax -> get
equally-spaced loops cheaply (so you do not have to measure dozens at strict equal intervals).

ROADMAP for radia_ngsolve.hysteresis: add an Everett-function identifier ``from_measured_loops(loops)``
returning piecewise-linear shape functions (then ``PlayHysteresis(eta, f_k=...)`` is a real material).
``radia.hysteresis_io`` already reads .hys / .mat loop files and builds shape functions on the C++ side.
"""

FE_AND_LOSS = """\
## FE coupling and iron loss

IRON LOSS (the immediate payoff -- replaces the fitted coreloss k_h).  Take the per-element flux-density
waveform B(theta) from a rotor sweep (the committed cogging / 3-phase AGE sweeps), run it through the
B-input play operator per element, and the hysteresis loss is the REAL loop area f * |∮ H dB| -- no
fitted k_h.  Classical-eddy + excess (Bertotti, coreloss.py) add on top; minor loops are handled
natively (the play operator tracks every reversal), which is exactly the iGSE refinement coreloss.py
flags as a TODO.

FE-COUPLED SOLVE (the core increment).  Make the constitutive H(B) history-dependent at each Gauss
point: per Newton/Picard step, given B = curl A, the play operator returns H and its differential
reluctivity dH/dB.  STATUS in radia_ngsolve.hysteresis:
- The DH/DB TENSOR is implemented and verified -- ``PlayHysteresis.play_tangent(Bx,By,px,py)`` returns
  the consistent  dH/dB = a_0 I + sum_{moving} a_k [I - (eta_k/|s_k|)(I - s^_k (x) s^_k)]  (s_k=B-p_k;
  pinned cells contribute 0), SPD by construction, matching a finite-difference of H to ~1e-11 incl. the
  off-diagonal vector coupling (Mitsuoka-Mifune-Matsuo, IEEE TMAG 2013 -- the Newton-Raphson ingredient).
- The FIXED-POINT (M-source) solve is a ROBUST working solver (test_hysteresis_fe_coupled): a
  current-driven ring core, per-element play states on the mesh, FE loop area ~ the model (~11%), ~25
  M-source sweeps/step.
- The VARIATIONAL (energy-line-search) NEWTON is the principled globally-convergent solver and is now
  implemented & verified (test_hysteresis_fe_variational).  ``PlayHysteresis.play_coenergy_density`` is
  the CONVEX incremental co-energy psi*(B) whose B-gradient is H and whose Hessian is the SPD
  play_tangent (Francois-Lavet-Henrotte 2013; the rate-independent-plasticity structure: p_prev <->
  back-stress, eta_k <-> yield radius, the play update <-> radial return).  So Pi(A)=INT psi*(B(A)) -
  INT J.A is CONVEX and a Newton step damped by an ARMIJO LINE SEARCH ON Pi (not the residual norm) is
  a guaranteed descent direction -> converges in a HANDFUL of iters/step (vs ~25 fixed-point sweeps),
  energy monotone, steady loop area matching the play model TIGHTER than fixed-point (~1.4% vs ~11%).
  NOTE: an earlier Newton attempt "diverged"; the cause was traced to an INEXACT element-average B
  (gfB.Set's L2-projection at order 1 inflates |B|, amplified by nu_air) -- NOT an intrinsic kink
  pathology.  Taking B exactly via Integrate(.,element_wise) on order-1 H1 (B=curl A is then P0 per
  element), BOTH plain and line-searched Newton converge; the energy line search is the rigorous
  globalisation/safety-net (guaranteed descent under deeper saturation / worse starts).
- This all extends the AGE nonlinear Picard (``airgap_motor_workflow.age_motor_nonlinear_solve``) from
  single-valued nu(|B|) to a hysteretic, history-carrying constitutive -- the play state p_k per element
  is the internal variable.
"""

SECTIONS = {
    "landscape": LANDSCAPE,
    "play_formulation": PLAY_FORMULATION,
    "binput_congruency": BINPUT_CONGRUENCY,
    "identification": IDENTIFICATION,
    "fe_and_loss": FE_AND_LOSS,
}


def get_hysteresis_knowledge(topic: str = "landscape") -> str:
    """Return magnetic-hysteresis modelling knowledge.  ``topic="all"`` returns everything."""
    t = topic.strip().lower()
    if t == "all":
        return "\n\n---\n\n".join(SECTIONS[k] for k in SECTIONS)
    if t not in SECTIONS:
        return f"Unknown topic {t!r}. Valid: {', '.join(sorted(SECTIONS))}\n"
    return SECTIONS[t]
