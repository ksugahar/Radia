r"""Time-domain IABC (Improvised/Improved Asymptotic Boundary Conditions /
spherical absorbing conditions based on electrical images) -- methodology and
verified design know-how for open-boundary truncation of FEM/FETD models.

An IABC truncates an unbounded exterior by concentric virtual material shells
whose (generally complex, frequency-dependent) eps, mu act as ELECTRICAL IMAGES
that cancel the reflection of the radiating multipole at the truncation surface,
backed by a perfect-conductor (PEC) outer wall.  It is designed PER FREQUENCY by
nulling the physical reflection coefficient of each multipole order, and it
degenerates to the classical static absorbing BC (Asymptotic BC / ballooning) in
the omega->0 limit.  This module documents how to extend that frequency-domain
design to the TIME DOMAIN (transient FETD / Newmark-beta) via a causal relaxation
network (see the `urn` tool) and an auxiliary-differential-equation (ADE)
material update, and records the verified design lessons.

References (open):
  K. Sugahara, "Spherical Absorbing Conditions Based on Electrical Images."
  K. Sugahara & Y. Sato, "Universal Relaxation Network", IEEE Access 2026
    (see the `urn` / `urn_fit` tools).
  J.-M. Jin, "The Finite Element Method in Electromagnetics", 3rd ed., Ch. 12
    (Newmark-beta time-domain FEM).
  Context: Bayliss-Turkel radiation BCs; Berenger PML; Dirichlet-to-Neumann maps.
"""

IABC_OVERVIEW = r"""
# IABC -- improvised/improved asymptotic (absorbing) boundary conditions

Goal: truncate an open exterior with a THIN region near the boundary instead of
meshing far into free space (or wrapping the domain in a PML).

Construction:
- Surround the truncation surface (circle in 2D, sphere in 3D) with one or more
  concentric virtual shells of material (eps_k, mu_k), and close with a PEC wall.
- The shells act as ELECTRICAL IMAGES: their (generally complex) eps, mu are
  chosen so the multipole field radiated from inside sees (near) zero reflection
  at the truncation surface.
- Design is PER MULTIPOLE ORDER and PER FREQUENCY: for each order l (or n in 2D)
  pick the shell material that nulls the reflection coefficient Gamma_l(omega).

Three regimes of the SAME operator (exterior wavenumber sets the character):
- omega -> 0  (static):     shells are real; IABC reduces to the classical static
                            asymptotic/ballooning BC; exterior multipole ~ r^{-(l+1)}.
- lossy/diffusive exterior: complex eps,mu; ties to a conductor DtN symbol.
- wave exterior (Helmholtz): spherical Hankel exterior; "improved" shells add a
                            Bayliss-Turkel-like curvature hierarchy.

Why time-domain is hard, in one line: a per-frequency-optimised (eps(omega),
mu(omega)) is generally NON-CAUSAL and the shell's reflection null is extremely
high-Q (a ~0.006% material error flips |Gamma| from 0 to 1), so it cannot be
marched in time as-is -- see topic `causality`.  The clean fix is to drop the
shell and realise the IABC's EFFECTIVE TERMINATION IMPEDANCE (the exact exterior
DtN) directly as a Robin BC + auxiliary ODEs -- exact, unconditionally stable,
seamless DC->radiation: see topic `exact_dtn` (the recommended route).

Relation to other truncations:
- DtN / Bayliss-Turkel: an IABC shell+PEC is a LOCAL material realisation of an
  absorbing operator; the exact nonlocal map is a Dirichlet-to-Neumann symbol.
- PML: both are absorbers near the boundary; a graded multilayer IABC optimised
  for broadband absorption emerges as a tapered impedance-matched stack, i.e. a
  discrete PML-like profile (topic `broadband_design`).
"""

IABC_FREQUENCY_DOMAIN = r"""
# Frequency-domain design (physical reflection of the multipole)

Spherical case (dipole l, normalised eps0=mu0=c=1, omega=k0, e^{+j w t}):
- Radial solutions per order l are spherical Bessel/Hankel; build them from
  half-integer cylindrical functions:  sj_l(z)=sqrt(pi/2z) J_{l+1/2}(z),
  sy_l(z)=sqrt(pi/2z) Y_{l+1/2}(z)  (scipy jv/yv accept complex z, so a complex
  shell wavenumber k=omega*sqrt(eps*mu) is handled directly).
- Transfer-matrix through each shell with the 2x2 layer matrix whose rows are the
  field and its conormal flux  (R, (1/p) dR/dr),  p = mu for the "mu-system" and
  p = eps for the "eps-system" (the two vector-potential polarisations).
- Outer PEC closes the stack: Dirichlet (field=0) for the mu-system, Neumann
  ((1/p)dR/dr=0) for the eps-system.
- Reflection of the OUTGOING multipole at the truncation radius:
      h2 = sj - i sy (outgoing),  h1 = sj + i sy (incoming);
      Gamma_l = -(row . T . v_out)/(row . T . v_in),
  with row=[1,0] (Dirichlet) or [0,1] (Neumann).  |Gamma|<=1 for a passive stack.

Per-frequency design: at each omega solve the small nonlinear system that drives
Gamma_l(omega)=0 for BOTH polarisations (4 real residuals -> 4 real unknowns =
Re/Im of one shell eps and mu).  Warm-start from DC upward for a smooth branch.
The omega->0 limit reproduces the static absorbing-BC material (the static IABC).

VERIFIED (open numpy/scipy):
- A single 1-shell dipole design nulls both polarisations to residual ~2e-9 at
  every omega (exact per-frequency null).
- The resulting ideal mu(omega) is NON-MONOTONIC (rises above its DC value to a
  resonant peak, then falls); ideal eps(omega) is small with a SIGN-CHANGING
  imaginary part (i.e. the per-frequency ideal is partly non-passive).  Both
  facts drive the causality discussion below.

Vectorise over omega (batched Bessel + batched 2x2 inverse, air-side Hankel
precomputed once): the scalar per-omega loop is ~100x slower and makes any
multistart optimisation run away.
"""

IABC_CAUSALITY = r"""
# Causality: why naive time-domain realisation fails, and what works

A per-frequency IABC gives (eps(omega), mu(omega)) sampled independently at each
omega.  Two failure modes when you try to march it in time:

1. NON-CAUSAL POINT VALUES.  Dropping a complex eps/mu sampled at one frequency
   into the time domain violates Kramers-Kronig (non-causal, unstable).

2. ILL-CONDITIONED MATERIAL->REFLECTION MAP.  Even fitting the ideal eps(omega),
   mu(omega) with a causal model and THEN evaluating the reflection FAILS:
   VERIFIED -- a URN fit of the ideal shell materials at <1% NRMSE
   (NRMSE_mu~0.8%, NRMSE_eps~3%) still produces a broadband reflection that
   BLOWS UP to max|Gamma| ~ 1.3 (>1).  The per-frequency design sits on a sharp,
   high-Q reflection null, so a sub-1% material error moves it off the null to
   O(1) reflection; and forcing a passive fit onto the partly non-passive ideal
   eps injects gain (|Gamma|>1).
   => Lesson: do NOT fit the ideal material then hope the reflection follows.
      Optimise the CAUSAL material DIRECTLY against the reflection objective.

3. NON-MONOTONIC TARGET.  The ideal mu(omega) overshoots its DC value (resonant
   peak).  A sum of Debye relaxations is monotonic and CANNOT reproduce it; a
   causal fit therefore needs RESONANT (Lorentzian / RLC) terms in addition to
   Debye, or one abandons material-matching for direct reflection optimisation.

Correct approach (passive-by-construction + direct objective):
- Parametrise the shell as a PASSIVE relaxation network -- a sum of Debye
  (Delta>=0, tau>0) and/or Lorentzian/RLC terms, optionally with a conductive
  (telegrapher) part eps'-j sigma/omega.  This is causal and |Gamma|<=1 BY
  CONSTRUCTION.
- Optimise its parameters (and the shell thickness) to minimise the broadband
  reflection directly (topic `broadband_design`).
- Realise the resulting poles in time with ADE / Newmark-beta (topic
  `timedomain`).  The poles map 1:1 to auxiliary ODEs; conductivity maps to a
  damping matrix with no auxiliary field at all.

Fundamental limit (state it, don't hide it):
- Broadband causal absorption only works where the truncation surface sees a
  PROPAGATING field, i.e. k*R0 >~ 1.  In the low-frequency EVANESCENT regime
  (k*R0 < 1) the multipole at the boundary is quasi-static (no propagating wave
  to absorb); no causal broadband absorber can null it -- only the exact
  per-frequency static image does, and the causal material's DC limit recovers
  that static IABC as omega->0.  The "price of causality" is concentrated in the
  evanescent band.
"""

IABC_TIMEDOMAIN = r"""
# Time-domain realisation (ADE + Newmark-beta)

Once the shell material is a causal relaxation network, each mechanism is a
LOCAL-in-time auxiliary ODE -- no history convolution:

- Conductive / telegrapher  eps(omega)=eps' - j sigma/omega :
    a DAMPING matrix R = int sigma N N ; NO auxiliary field.
- Debye  Delta/(1+j w tau) :
    one auxiliary field P with  tau dP/dt + P = Delta * u ,  plus a
    negative-stiffness term -- exactly one first-order ODE per pole.
- Lorentzian / RLC (resonant, needed for a non-monotonic mu) :
    one second-order auxiliary ODE (or two first-order) per pole.
- Fractional (Cole-Cole/CPE) :
    a short RC/RL ladder (a few auxiliary ODEs) -- emitted directly by the URN
    SPICE synthesis (the `urn_fit` tool).

Newmark-beta (Jin Ch.12; gamma=1/2, beta=1/4 unconditionally stable):
    K0 u_{n+1} = K1 u_n - K2 u_{n-1} + b_bar,
    K0 = T/dt^2 + R/(2 dt) + beta S,
    K1 = 2T/dt^2 - (1-2 beta) S,
    K2 = T/dt^2 - R/(2 dt) + beta S,
with T (eps mass), R (sigma damping), S (1/mu stiffness).  K0 is time-invariant
=> factor once, reuse every step.  Passivity of the relaxation network guarantees
stable stepping.

Reduced RADIAL form for one multipole order l (substitute u=rR; field R_l(r,t)):
    T  = int eps'(r) r^2 N N
    R  = int sigma(r) r^2 N N
    S  = int (1/mu'(r)) [ r^2 N' N' + l(l+1) N N ]     <- centrifugal l(l+1)/r^2
PEC outer => Dirichlet (mu-system).  This is the spherical analogue of a 1D
telegrapher FETD; a graded conductive+real-mu absorber needs only the damping
matrix R (no auxiliary ODE), which is the simplest stable broadband realisation.

VERIFIED (open): a single-pole Debye realised by the P-field ADE reproduces the
analytic Fresnel reflection of the dispersive medium in the time domain (in-band
|Gamma| within a few percent); a planar graded telegrapher absorber driven by a
band-limited pulse absorbs ~96% (time-domain reflection ~4%).
"""

IABC_BROADBAND_DESIGN = r"""
# Broadband design recipe (the optimisation that replaces single-omega fsolve)

Single-frequency design (solve Gamma=0 at one omega) is exact at that omega but
the reflection spikes to ~1 elsewhere ("the frequency response looks wrong").
Make it broadband and time-domain-realisable:

1. PASSIVE parametrisation (guarantees |Gamma|<=1 and causality):
   per layer  eps_j(omega)=eps'_j - j sigma_e_j/omega ,  mu'_j real
   (sigma_m=0 keeps 1/mu non-dispersive => cleanest time domain).  Use softplus
   for eps'_j, sigma_e_j, mu'_j >= 0.
2. STACK the reflection over a frequency GRID covering the band, for BOTH
   polarisations, into one residual vector.
3. THICKNESS as a free variable (single thin shell cannot be broadband); allow a
   GRADED MULTILAYER (N concentric layers, conductivity increasing outward).
4. GLOBAL search: multistart trf least-squares with graded physical inits beats
   a single local solve (and beat differential-evolution minimax in tests).

What emerges: a TAPERED IMPEDANCE-MATCHED ABSORBER -- inner layers eps'~mu'~1
(eta=sqrt(mu/eps)~1, matched to vacuum) with conductivity ramping up outward
(a discrete PML-like profile).

VERIFIED (open numpy/scipy, dipole l=1, truncation radius R0=1):
- Single thin causal shell, full band incl. evanescent: max|Gamma| ~ 1
  (limited by the un-absorbable low-omega evanescent regime; see `causality`).
- Restricting to the PROPAGATING band: single thin shell already gives
  max|Gamma| ~ 0.12 (kR0>=2), ~0.31 (kR0>=1).
- GRADED MULTILAYER on kR0 in [1,6]: N=3 -> max|Gamma| 0.092, N=4 -> 0.082
  (mean ~0.045), thickness H ~ 4-5 (~ one wavelength at the low-band edge).
- Planar analogue (normal incidence): 9 graded telegrapher layers ->
  max|Gamma| ~ 0.09, |Gamma| ~ 1-6% flat across the band.

The optimised (eps'_j, sigma_e_j, mu'_j) ARE the Newmark telegrapher
coefficients, so the design drops straight into the FETD of topic `timedomain`.
"""

IABC_APPLICATION = r"""
# Putting it together: a transient open boundary

Pipeline (frequency design -> causal network -> time stepping):
  1. Frequency-domain: design the absorbing shell(s) by direct broadband
     reflection optimisation over the band of interest (topic
     `broadband_design`); keep the material passive so it is causal.
  2. If the design is expressed as a target eps(omega)/mu(omega) or a boundary
     admittance / DtN symbol, identify its causal relaxation network with the
     `urn_fit` tool -> Debye/Lorentzian/fractional mechanisms + SPICE/ADE.
  3. Time domain: assemble T (eps), R (sigma), S (1/mu) with the centrifugal
     term per multipole order; add one auxiliary ODE per Debye/Lorentzian pole
     (conductivity needs none); march with Newmark-beta (topic `timedomain`).
  4. Validate: drive a band-limited multipole pulse, compare the time-domain
     reflection (reflected/incident spectrum vs an air-only far-wall reference)
     to the frequency-domain |Gamma(omega)| the absorber was designed for, and
     to an equal-thickness PML.

When to use which absorber:
  - Pure broadband absorption on a propagating band: a graded conductive
    (telegrapher) multilayer -- simplest, needs only the damping matrix.
  - A target with genuine dispersion / a non-monotonic mu (resonant peak), or a
    DtN symbol to match: fit it with URN (Debye + Lorentzian/RLC + fractional)
    and realise each mechanism as an auxiliary ODE.
  - Quasi-static / evanescent content (kR0<1): handled by the static-IABC limit
    (the causal material's DC value); broadband absorption is not available
    there and should not be claimed.
"""

IABC_EXACT_DTN = r"""
# The exact-impedance route (RECOMMENDED for a transient open boundary)

The cleanest time-domain open boundary uses NO absorbing shell at all.  An IABC
shell+PEC only presents an effective surface impedance = the exact exterior
Dirichlet-to-Neumann (DtN) symbol G_l(omega) at the truncation sphere.  Realise
THAT impedance directly as a Robin BC + auxiliary ODEs, instead of via a
fictitious lossy material (which is high-Q and does not causalise -- see
`causality`).

Air (wave) exterior: G_l is EXACTLY RATIONAL in s (outgoing h_l^{(2)} =
e^{-jx}/x * polynomial(1/x), so e^{-jx} cancels in the log-derivative):
    R0*G_l(s) = -s - 1 - N_l(s)/Q_l(s),   Q_l = reverse Bessel polynomial,
    deg l, ALL roots Re<0.
      l=1:  -s - 1 - 1/(s+1)                          pole s=-1
      l=2:  -s - 1 - 3(s+2)/(s^2+3s+3)                poles -1.5 +- 0.866 j
      l=3:  -s - 1 - 3(2s^2+10s+15)/(s^3+6s^2+15s+15) poles -2.32, -1.84+-1.75 j
    s->0 gives the static -(l+1)/R0; s->inf gives the Sommerfeld -s/c.
    => seamless DC -> evanescent -> radiation BY CONSTRUCTION (no fitting).

Time-domain realisation: each mode l = one l-state companion auxiliary-ODE
system (controllable-canonical realisation of N_l/Q_l) + a boundary time
derivative (the -s) + a constant (the -1), folded into Newmark-beta.  Poles
Re<0 + dissipative boundary => UNCONDITIONALLY STABLE, reflection zero up to
discretisation.  VERIFIED by a 1D radial FETD solve (Newmark-beta, R0=c=1,
demo_uu2_exact_dtn_fetd.py): realising the exact DtN as a Robin BC + l companion
auxiliary ODEs, the spurious reflection (truncated solve vs a free-space
reference on the SAME interior mesh) falls as O(h^2) for l=1,2,3 -- ~5.7e-4,
~1.4e-4, ~3.6e-5 on 100/200/400-element radial meshes (clean x4 per h-halving)
-- i.e. it is reflectionless in the continuum (discretisation-limited only).  It
beats a 1st-order Sommerfeld boundary on the same mesh by ~2000x, and the
interior energy drains to machine precision (~1e-16 of peak) without ever
growing => passive.  This is the right vehicle for "seamless, stable, transient
open boundary": realise the exact PASSIVE boundary response, do not causalise a
material.  (The DtN symbol / pole realisation itself is verified separately in
demo_uu_iabc_time_domain.py.)

IMPORTANT (prior art, state it): for a SEPARABLE boundary (sphere / cylinder /
half-space) this exact rational radiation BC is the classical Grote-Keller /
Bayliss-Turkel / Hagstrom-Warburton exact / high-order radiation condition --
reflectionless, and it does outperform a PML there, but that is established and
it does NOT generalise to arbitrary geometry (which is exactly where PML is
used).  So this is NOT "beating PML" in general.  The IABC-specific value is the
bridge "IABC shell == this exact termination impedance" and its PASSIVE
EQUIVALENT-CIRCUIT (relaxation-network / TLM) realisation -- fit/synthesise with
the `urn` / `urn_fit` tools.

Lossy / conductive (eddy-current / diffusion) exterior: G_l is NOT rational in s
(a sqrt(s) branch cut) -- BUT it IS EXACTLY RATIONAL in q=sqrt(s), with the SAME
reverse Bessel poles as the wave case above (gamma=ik -> gamma=sqrt(s)):
    G_n(s) = -s a^2 theta_{n-1}(a sqrt(s)) / theta_n(a sqrt(s)) - (n+1).
So a Cauer / CLN continued fraction IN q realises it EXACTLY with n+1 stages and
WELL-CONDITIONED (coeff spread ~1-30) -- decisively beating a Foster fit in s
(which floors ~1e-3 at 32 states and ill-conditions to ~1e5) and the URN (which
plateaus ~1.7e-2; URN's edge is measured/multi-mechanism data, not a clean
analytic DtN).  This promotes the lab CLN (Kameari-Sugahara eddy-current MOR) to
the open BOUNDARY itself; the sqrt(s) diffusion-memory element is realised by the
finite eddy-current ladder (real negative poles -> stable, finite auxiliary ODEs).
VERIFIED two ways: demo_xx3_cln_dtn_cauer.py (the EXACT Cauer in q=sqrt(s), n+1
stages -- the structural/unification result; wave realisation = demo_uu/uu2) and
demo_xx4_cln_mor_radial_eddy.py (the GENUINE lab CLN: a Lanczos/PVL MOR of a
~700-DOF radial eddy-current FEM down to a COMPACT ~16-stage INTEGER-ORDER Cauer
ladder reproducing G_n, monotone convergence, SPD => real negative poles =>
directly time-domain + stable -- the practical eddy-current open boundary).  And
demo_xx5_cln_fetd_reflection.py USES it as a TRANSIENT open boundary: a Crank-
Nicolson eddy-current diffusion FETD with the exterior Krylov-substructured to
~16 DOFs gives ~1e-6 spurious reflection (vs ~6-11% for Dirichlet/Neumann
truncation), monotone in the stage count, across n=1,2,3 -- the diffusive analog
of demo_uu2's wave reflection test.  So the reverse-Bessel/CLN open boundary is
reflectionless in time for BOTH the wave and the diffusion regime.
"""

_TOPICS = {
    "overview": IABC_OVERVIEW,
    "frequency_domain": IABC_FREQUENCY_DOMAIN,
    "causality": IABC_CAUSALITY,
    "exact_dtn": IABC_EXACT_DTN,
    "timedomain": IABC_TIMEDOMAIN,
    "broadband_design": IABC_BROADBAND_DESIGN,
    "application": IABC_APPLICATION,
}


def get_iabc_documentation(topic: str = "all") -> str:
    """Return IABC knowledge text.  topic in
    {all, overview, frequency_domain, causality, exact_dtn, timedomain,
     broadband_design, application}."""
    t = (topic or "all").strip().lower()
    if t == "all":
        return "\n".join(_TOPICS[k] for k in
                         ["overview", "frequency_domain", "causality",
                          "exact_dtn", "timedomain", "broadband_design",
                          "application"])
    if t in _TOPICS:
        return _TOPICS[t]
    return (f"Unknown topic '{topic}'. Options: all, " + ", ".join(_TOPICS) + ".")
