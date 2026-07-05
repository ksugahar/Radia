"""Sommerfeld integral / layered-medium Green's function knowledge module.

Absorbed 2026-05-26 from
``public-safe curated corpus``
(3 PDFs):

  1. ``Sommerfeld_Integral_Weyl_Identity.pdf`` -- Chew lecture notes (lecture
     35, pp. 349-370), pedagogical derivation of the Weyl + Sommerfeld
     identities from the 3D scalar Helmholtz equation.  Establishes the
     plane-wave spectrum representation of a point source and the
     Sommerfeld Integration Path (SIP).  Time convention exp(-i omega t).

  2. ``Exact_ClosedForm_Sommerfeld_Impedance_Plane.pdf`` -- Koh & Yook,
     IEEE Trans. Antennas Propagat. 54(9), 2568-2576, Sept 2006.
     DOI: 10.1109/TAP.2006.880747.  EXACT closed-form Sommerfeld integral
     for a dipole over an impedance half-plane via incomplete Weber
     integrals / incomplete Lipschitz-Hankel integrals / Lommel functions.
     Valid for ANY complex normalised impedance and ANY (rho, z+z'),
     unlike the patchwork Norton / Banos / Wait asymptotic developments.

  3. ``Spectral_Expansions_Sommerfeld.pdf`` -- Chew lecture 35 (same source
     as #1 from a slightly different revision), same Weyl + Sommerfeld
     derivation + stationary-phase asymptotic.  Does NOT cover DCIM /
     GPOF / VECTFIT / etc. (those are separate knowledge gaps; see
     NUMERICAL_EVALUATION below for the canonical references).

This module documents the THEORY and the SPECTRAL FORMULATION.  It does
NOT ship runnable code -- Sommerfeld integrals appear in lab production
only at the dipole-over-impedance-plane MOM kernel (rare, mostly cited
for ground-penetrating-radar and antenna-over-earth use cases).  For the
actual layered-medium MoM workflow the canonical implementation is
NGSolve's external ``ngsolve.bem`` Helmholtz solver coupled to a
Sommerfeld-integral Green's function table, generated offline by one of
the acceleration methods listed below.

POLICY: Radia core (MMM/MSC) is Laplace-kernel only.  This module is
the THEORY POINTER for users who need the FULL layered-medium Green's
function (Sommerfeld integrals) -- the production answer is usually:
either go to ngsolve.bem with a Helmholtz kernel + DCIM-accelerated
spectral domain, OR drop down to the impedance-plane closed form
(Koh-Yook 2006) if a single layer is enough.
"""


OVERVIEW = r"""
# Sommerfeld integral / layered-medium Green's function -- overview

## When do Sommerfeld integrals appear?

A Sommerfeld integral is the 1D spectral representation of a point-source
field above a planarly stratified medium.  It appears whenever:

  * a dipole (electric or magnetic) sits above OR inside a layered
    background (earth, multilayer dielectric, microstrip substrate,
    impedance ground plane, ...), AND
  * the time-harmonic Maxwell solution is needed AT another point in
    the same stack (microstrip-MoM source/observation Green's function,
    HF radio ground wave, ground-penetrating radar, subsurface comms,
    log-periodic antenna over earth).

## What does it solve?

The Sommerfeld integral evaluates the Maxwell Green's function for a
planarly layered medium.  The integral represents the field as a
continuous superposition of plane waves (propagating + evanescent),
each reflected/transmitted through the stack by analytic plane-wave
reflection / transmission coefficients computed via the standard
recursion (TM and TE polarisations evolve independently).

## Two equivalent identities (Chew lecture 35, exp(-i omega t) convention)

WEYL IDENTITY (35.1.11) -- spherical wave as 2D plane-wave integral:

    exp(i k_0 r) / r = (i / (2 pi)) *
        integral_(-inf)^(+inf) integral_(-inf)^(+inf)
            exp(i kx x + i ky y + i kz |z|) / kz   dkx dky

    with kx^2 + ky^2 + kz^2 = k_0^2,  Im(kz) > 0 (physical Riemann sheet).

The (kx, ky) disk |k_rho| <= k_0 carries propagating modes; the annulus
|k_rho| > k_0 carries evanescent modes that reconstruct the near-source
singularity.  The pole locus is the EWALD SPHERE (Ewald 1921).

SOMMERFELD IDENTITY (35.1.14 / 35.1.15) -- 1D Bessel/Hankel integral:

    exp(i k_0 r) / r = i *
        integral_0^inf  (k_rho / kz) J_0(k_rho rho) exp(i kz |z|)  dk_rho

    equivalently with H_0^(1) and integration over the SIP from
    -infinity to +infinity:

    exp(i k_0 r) / r = (i/2) *
        integral_(-inf)^(+inf)
            (k_rho / kz) H_0^(1)(k_rho rho) exp(i kz |z|)  dk_rho

These convert the 2D (kx, ky) integral to a single radial integration
over k_rho = sqrt(kx^2 + ky^2), at the price of a Bessel J_0 (or Hankel
H_0^(1)) kernel.

## TIME CONVENTION WARNING

Chew uses exp(-i omega t) (physics convention).  Stoll, Radia C++,
Sugahara lab CLN papers, and most engineering EM literature use
exp(+j omega t) (engineering convention).  TO CONVERT:

  * imaginary unit i in physics ↔ -j in engineering
  * "outgoing" spherical wave exp(+i k r) / r  ↔  exp(-j k r) / r
  * Hankel H_0^(1) (physics outgoing) ↔ H_0^(2) (engineering outgoing)
  * Conjugate any closed-form expression copied from physics texts.

The lab default is exp(+j omega t).  Apply the conjugate rule when
porting from Chew, Brekhovskikh, Felsen-Marcuvitz, etc.
"""


WEYL_DERIVATION = r"""
# Weyl identity -- derivation in detail

## Starting point: scalar Helmholtz Green's function

The 3D scalar Helmholtz equation with a point source

    ( d^2/dx^2 + d^2/dy^2 + d^2/dz^2 + k_0^2 ) phi(x,y,z)
        = -delta(x) delta(y) delta(z)

has the unique outgoing-wave solution

    phi(r) = exp(i k_0 r) / (4 pi r)        (Eq. 35.1.2 in Chew)

## Spectral representation by 3D Fourier transform

Apply the 3D Fourier transform.  The transformed delta is 1, the
Laplacian becomes -(kx^2 + ky^2 + kz^2), so

    tilde phi(kx, ky, kz) = 1 / (kx^2 + ky^2 + kz^2 - k_0^2)

Inverse transform:

    phi(x,y,z) = -1/(2 pi)^3 *
        triple integral over (kx, ky, kz) of
            exp(i kx x + i ky y + i kz z)
            / (k_0^2 - kx^2 - ky^2 - kz^2)  dkx dky dkz

The integrand has POLES at  kz = +/- sqrt(k_0^2 - kx^2 - ky^2)  which
lie on the real kz axis for a lossless k_0.  REGULARISATION: add a
small positive imaginary part to k_0,

    k_0 = k_0' + i k_0''   with   k_0'' > 0+

which physically corresponds to infinitesimal medium loss and
mathematically enforces the radiation condition (Sommerfeld 1909).

## Contour deformation: closing in the upper-half kz plane

For z > 0, the integrand decays in the upper half kz plane, so by
Jordan's lemma the semi-circular contour at infinity contributes zero.
By Cauchy's theorem, the real-axis kz integral equals the residue at
the UHP pole  kz = +sqrt(k_0^2 - kx^2 - ky^2).  Evaluating the residue
collapses the triple (kx, ky, kz) integral to a DOUBLE (kx, ky)
integral.  The z < 0 case is handled identically in the LHP.  Combining
both with |z|:

    phi(x,y,z) = i / (2 (2 pi)^2) *
        double integral over (kx, ky) of
            exp(i kx x + i ky y + i kz |z|) / kz  dkx dky    (35.1.10)

Equating with the closed form (35.1.2) yields the WEYL IDENTITY (35.1.11):

    exp(i k_0 r) / r = (i / (2 pi)) *
        double integral over (kx, ky) of
            exp(i kx x + i ky y + i kz |z|) / kz  dkx dky

with the constraint  kx^2 + ky^2 + kz^2 = k_0^2  and the Riemann sheet
rule  Im(kz) > 0  AND  Re(kz) > 0  for every (kx, ky).

## Physical interpretation

The Weyl identity expands a SPHERICAL wave (LHS) as a continuous integral
sum of PLANE WAVES (RHS) propagating in every direction, including
evanescent ones.  The EWALD SPHERE  kx^2 + ky^2 + kz^2 = k_0^2  is the
locus of poles in spectral space.  The disk  |k_rho| < k_0  carries
propagating modes; outside it the modes are evanescent in z.
"""


SOMMERFELD_DERIVATION = r"""
# Sommerfeld identity -- 1D radial form

## Change to polar coordinates

In the (kx, ky) plane use polar coordinates  k_rho, alpha  with
kx = k_rho cos(alpha), ky = k_rho sin(alpha).  In the (x, y) plane use
rho, phi  with  x = rho cos(phi), y = rho sin(phi).  Then

    dkx dky = k_rho dk_rho d(alpha)
    kx x + ky y = k_rho rho cos(alpha - phi)
    kz = sqrt(k_0^2 - k_rho^2)

## Inner alpha integral via Bessel J_0

The Bessel integral representation (Whittaker-Watson 1927, Sec. 17.21):

    J_0(k_rho rho) = (1 / 2 pi) *
        integral_0^(2 pi)  exp(i k_rho rho cos(alpha - phi))  d(alpha)

Substituting into the Weyl identity and performing the alpha integral
gives the SOMMERFELD IDENTITY (35.1.14):

    exp(i k_0 r) / r = i *
        integral_0^inf  (k_rho / kz) J_0(k_rho rho) exp(i kz |z|)  dk_rho

## Hankel-form symmetric extension (35.1.15)

Using J_0 = (1/2)(H_0^(1) + H_0^(2))  with the Hankel reflection
identity  H_0^(1)(exp(i pi) x) = -H_0^(2)(x)  the integral extends to
the full real line:

    exp(i k_0 r) / r = (i/2) *
        integral_(-inf)^(+inf)
            (k_rho / kz) H_0^(1)(k_rho rho) exp(i kz |z|)  dk_rho

## Physical interpretation

A spherical wave decomposes into CONICAL or CYLINDRICAL waves propagating
radially in rho with a plane wave in z.  For k_rho > k_0, kz is purely
imaginary and exp(i kz |z|) is exponentially decaying in z (evanescent).
"""


LAYERED_MEDIUM_DIPOLES = r"""
# Dipole over a planarly stratified medium

## TM / TE decomposition

Plane waves reflecting from a planarly stratified medium decouple into:

  * TM (Transverse Magnetic): Hz = 0, Ez != 0, reflection coefficient
    tilde R_{i,i+1}^{TM}  (generalised through the stack via recursion).
  * TE (Transverse Electric): Ez = 0, Hz != 0, reflection coefficient
    tilde R_{i,i+1}^{TE}.

Strategy: expand the primary dipole field via Weyl / Sommerfeld into
TM + TE plane-wave components, multiply each by its plane-wave
reflection / transmission coefficient (which is the same for every
k_rho due to phase matching across the stack), and re-integrate.

## Vertical Electric Dipole (VED) -- TM only

For alpha-hat = z-hat:  Hz = 0  identically, so only TM is excited.
The z-component of E in spectral form (35.2.5):

    Ez = -(I ell) / (8 pi omega eps) *
         integral_(-inf)^(+inf)
            (k_rho^3 / kz) H_0^(1)(k_rho rho) exp(i kz |z|)  dk_rho

Above a stratified medium (interface at z = -d_1), region 1 carries
incident + reflected wave (35.2.6):

    E_{1z} = -(I ell)/(4 pi omega eps_1)
        * integral_0^inf  dk_rho * (k_rho^3 / k_{1z}) J_0(k_rho rho)
        * [  exp(i k_{1z} |z|)
           + tilde R_{12}^{TM} exp(i k_{1z} z + 2 i k_{1z} d_1) ]

## Horizontal Electric Dipole (HED) -- TM + TE

For alpha-hat = x-hat, both Ez and Hz are non-zero.  After applying
the Sommerfeld identity inside d/dx and d/dy, only H_1^(1) and
cos(phi), sin(phi) azimuthal factors survive (35.2.12 / 35.2.13):

    E_{1z} = +(i I ell)/(4 pi omega eps_1) cos(phi)
        * integral_0^inf  dk_rho * k_rho^2 * J_1(k_rho rho)
        * [ +/- exp(i k_{1z} |z|)
            - tilde R_{12}^{TM} exp(i k_{1z} (z + 2 d_1)) ]

    H_{1z} = +(i I ell)/(4 pi) sin(phi)
        * integral_0^inf  dk_rho * (k_rho^2 / k_{1z}) * J_1(k_rho rho)
        * [  exp(i k_{1z} |z|)
            + tilde R_{12}^{TE} exp(i k_{1z} (z + 2 d_1)) ]

The +/- on the TM piece encodes the odd parity of Ez about z = 0.

## The reflection coefficients tilde R_{i,i+1}^{TM/TE}

These are computed by a STANDARD RECURSION through the stack (Chew
1995, ch. 2 of "Waves and Fields in Inhomogeneous Media"):

  * start at the deepest interface with the Fresnel coefficient,
  * propagate up layer-by-layer with the transmission-line analogy,
  * the result depends only on k_rho (phase-matched across all
    interfaces) and on the per-layer (mu_i, eps_i, d_i).

They are NOT part of the Sommerfeld identity itself -- treat them as
a black box called once per k_rho integrand evaluation.

## Sanity check (lab convention)

In the single-layer limit  tilde R^{TM} = tilde R^{TE} = 0,
the spectral integral (35.2.6) MUST reproduce the free-space dipole
field within discretisation tolerance.  This is the cheapest unit
test for any layered-medium Green's function code.
"""


BRANCH_CUTS_AND_SIP = r"""
# Branch cuts, Riemann sheets, and the Sommerfeld Integration Path

## Singularities of the Sommerfeld integrand

The radial integrand
    (k_rho / kz) H_0^(1)(k_rho rho) exp(i kz |z|)
has THREE types of singularity:

  1. LOGARITHMIC branch point of H_0^(1)(k_rho rho) at k_rho = 0
     (from H_0^(1)(x) ~ (2 i / pi) ln(x) as x -> 0).

  2. ALGEBRAIC branch points of kz = sqrt(k_0^2 - k_rho^2) at
     k_rho = +/- k_0.  The double-valued sqrt requires a Riemann
     sheet assignment.  The PHYSICAL sheet has Im(kz) > 0 (decaying
     away from the source).

  3. POLES of the layered-medium reflection coefficients tilde R^{TM/TE}
     at the SURFACE-WAVE k_rho values of the stack (these are the
     guided modes; their residues give the surface-wave contribution
     to the total field).

## Sommerfeld Integration Path (SIP)

The SIP is the deformed contour in the complex k_rho plane chosen so
that the integral (35.1.15) is well-defined EVEN in the lossless
medium where the algebraic branch points sit on the real axis.

  * Pass ABOVE the logarithmic branch point at k_rho = 0
    (consistent with the Hankel reflection convention).
  * Detour into the LHP just before k_rho = -k_0 and into the UHP
    just past k_rho = +k_0 so the integrand stays on the
    Im(kz) > 0 Riemann sheet.

In a LOSSY medium (k_0' + i k_0''), the branch points move into the
complex plane (-k_0 to LHP, +k_0 to UHP), and the SIP can be slid
back onto the real axis.  The lab's default numerical strategy is to
add an artificial loss (k_0'' = 1e-6 * k_0') for the analytic SIP
deformation, and remove it post-integration if the medium is supposed
to be lossless.

## Lossless-limit regularisation

Without an infinitesimal Im(k_0) > 0, the integrand of (35.1.10) is
strictly UNDEFINED (poles on the real axis, integral not absolutely
convergent).  Adding +i 0 to k_0 is the standard remedy and is
equivalent to imposing Sommerfeld's outgoing-radiation condition
(Schot 1992, "Eighty years of Sommerfeld's radiation condition").

## Tail behaviour (large k_rho)

For k_rho >> k_0:  kz ~ i k_rho, so  exp(i kz |z|) ~ exp(-k_rho |z|),
a Laplace-type decay.  This is the regime where evanescent modes
reconstruct the near-source 1/r singularity.  For numerical work:

  * Truncating k_rho at any finite cap UNDER-resolves the near-source
    singularity -- you cannot stop at k_0.
  * The standard preprocessing for accelerated quadrature is to
    SUBTRACT the quasi-static (k_rho -> infty) Laplace asymptote
    analytically and integrate only the (small) difference.
"""


KOH_YOOK_2006_CLOSED_FORM = r"""
# Koh-Yook 2006 -- exact closed-form Sommerfeld integral for the
# impedance plane

## Source

  Il-Suek Koh (Inha University) and Jong-Gwan Yook (Yonsei University),
  "Exact Closed-Form Expression of a Sommerfeld Integral for the
  Impedance Plane Problem," IEEE Trans. Antennas Propagat. 54(9),
  2568-2576, Sept 2006.  DOI: 10.1109/TAP.2006.880747

## Problem

Infinitesimal dipole at height z' above an impedance HALF-PLANE
characterised by a normalised surface impedance eta_s.  The dyadic
Green's function for ALL dipole orientations (HED, HMD, VED, VMD)
reduces to ONE primitive Sommerfeld integral

    I(p) = integral_0^inf  J_0(k_rho rho)
           * [exp(-u (z + z')) / (u + eta_s k_0)]  k_rho dk_rho

    u = sqrt(k_rho^2 - k_0^2),  radiation-condition branch.

All other entries of the dyadic Green's function are spatial
derivatives of I(p), so closing this one integral closes the
full half-space Green's function.

## What was previously open

  * Sommerfeld 1909: integral form, no closed-form evaluation.
  * Banos 1966: patchwork ASYMPTOTIC developments per region of
    (rho, z, z') space (far zone vs near zone vs grazing zone).
  * Norton ground-wave approximations: small-height asymptotics
    via Fresnel integrals; fail when both heights are non-trivial
    relative to wavelength.
  * Wait series: same intrinsic limitation.
  * Steepest descent (Felsen-Marcuvitz, Collin): geometric-optics
    term + lateral wave + cumbersome higher-order corrections.
  * Lindell-Alanen 1984 exact image theory: faster numerics but
    still requires numerical quadrature.
  * Chow et al. 1991 complex-image (Prony): discrete images;
    number of images grows when source & observation approach
    the surface, with accuracy degradation.

## Koh-Yook 2006 closed form

The paper builds the result in THREE layers:

  Layer A (Section II, eq. 8): EXPONENTIAL-INTEGRAL SERIES

    I(p) ~ exp(-w_0)/sqrt(rho^2 + (z+z')^2) * f_0(...)
           + sum_(n=1)^inf  a_n * E_n(arg_n)

    where E_n(x) is the exponential integral of order n.  Series is
    ABSOLUTELY CONVERGENT for any rho, z+z', and eta_s (real or
    complex, even non-physical).  Numerically the most stable.

  Layer B (Section III, eq. 10): LOMMEL-FUNCTION SERIES

    I(p) ~ sum_(n=0)^inf  c_n * U_n(omega, zeta)

    where U_n is the Lommel function of two variables.  Also
    absolutely convergent everywhere but numerically less stable
    because U_n grows fast.

  Layer C (Section IV, eq. 17): EXACT CLOSED FORM

    I(p) = (j pi / 2) H_0^(2)(k_0 r)
           + (closed-form impedance correction)

    r = sqrt(rho^2 + (z+z')^2)

    The correction is a finite combination of:
      * exp(-j k_0 r) / r        (image-source term)
      * H_0^(2)(k_0 r)           (Hankel function)
      * incomplete Weber integrals Q_nu(alpha, x), S_nu(alpha, x)
      * incomplete cylindrical function of Poisson form
      * incomplete Lipschitz-Hankel integrals Je_n(a, z), Ye_n(a, z)
      * exponential integral E_1 with a complex argument involving
        eta_s and k_0 r
      * polynomial coefficients in eta_s

    The "+/-" inside the closed form switches at  rho^2 - (z+z')^2 = 0
    (45 deg from the image).

## Special functions needed

  * Bessel J_0, J_n, Y_n
  * Hankel H_0^(2)
  * Exponential integral E_n (especially E_1)
  * Incomplete Weber integrals Q_nu, S_nu (Watson notation)
  * Incomplete cylindrical function of Poisson form
  * Incomplete Lipschitz-Hankel integrals (Dvorak 1995-96)
  * Lommel function U_n (intermediate only)
  * Legendre P_nu (asymptotic series, Sec. V)
  * Generalised hypergeometric 1F2 (intermediate only)

NOT needed: Fresnel integrals, error functions, modified Bessel.
(Norton's ground-wave family uses Fresnel; the impedance-plane
closed form does NOT.)

## Domain of validity

  * ANY complex eta_s, including the non-physical "eta_s = -10 + j"
    case the paper verifies explicitly (Fig. 6).
  * ANY rho >= 0, ANY z + z' > 0.  Near-field AND far-field both
    exact -- unlike steepest-descent, no near-field breakdown
    when source/observation are close to the surface.
  * Lossy ground (Re eta_s > 0) and PEC/PMC limits (eta_s -> 0,
    eta_s -> inf) both recovered analytically.

## Validation (paper figures)

All at f = 300 MHz unless stated; rho ~ O(wavelength); z + z' ~ 0.1 m.
Both eq. (8) and eq. (10) series and the closed form (eq. 17) match
direct numerical integration of (3) to many digits across all tested
parameter sweeps.  See Figs. 2-8 of the paper.

## Practical recipe

  * For RAPID evaluation, use eq. (8) (exponential-integral series).
  * For BENCHMARKING / closed-form reference, use eq. (17).
  * For LARGE k_0 r where (17) becomes expensive, use the Legendre
    P_nu asymptotic series of eq. (19) (Poincare divergent, truncate
    at minimum term).

## Application reach

  * HF / UHF ground-wave propagation (radio over earth).
  * Antennas near lossy ground (monopole, dipole, array; mutual
    coupling + ground reactance both matter).
  * Microstrip / printed circuits over impedance ground (Itoh 1989).
  * Subsurface comms, ground-penetrating radar surface response.
  * Wireless channel modeling at HF / UHF over forested or rough
    terrain (Sarabandi-Koh canopy-air interface work, cited).
  * Fast multipole acceleration for half-space Green's functions
    (companion FMM paper: Sarabandi-Koh, IEEE TAP 52(1), Jan 2004).

## Why this matters for lab

Lab production paths (Radia MMM/MSC, NGSolve FEM-Kelvin) handle
open-boundary problems by FAR-FIELD vanishing assumptions OR
Kelvin inversion.  They do NOT handle a dipole over a SEMI-INFINITE
LOSSY GROUND PLANE -- that requires a layered-medium Green's
function.  The Koh-Yook 2006 closed form is the RIGHT tool when:

  * a single layered interface is enough (semi-infinite earth model);
  * the lab needs an analytic kernel for MoM (lab will rarely build
    layered-medium MoM, but if asked: this is the kernel);
  * a reviewer asks "what about ground reflection at the antenna?"
    and a quick closed-form answer is needed.

For MULTI-LAYER stacks (microstrip, multi-layer PCB, geophysical
borehole), the Koh-Yook closed form is NOT enough -- use the
spectral integral (35.2.6 / 35.2.12 / 35.2.13) with numerical
quadrature accelerated by DCIM or GPOF (see NUMERICAL_EVALUATION).
"""


NUMERICAL_EVALUATION = r"""
# Numerical evaluation of Sommerfeld integrals -- canonical methods

The Chew lecture notes (paper 3 in this bundle) are PEDAGOGICAL --
they derive the Sommerfeld integral but do NOT survey numerical
acceleration.  The lab knowledge base therefore records the canonical
acceleration references separately.  This section is the POINTER LIST
for users who need to actually compute Sommerfeld integrals in code.

## Direct numerical quadrature

  * SIP integration: deform onto the standard Sommerfeld contour;
    Gauss-Patterson nested quadrature (already in
    `radia_mcp.radia_ngsolve.analytical_formulas`'s
    "gauss_patterson" topic) handles the oscillatory tail well.
  * For the tail (large k_rho), Gauss-Laguerre after the exponential
    substitution u = k_rho |z| works when |z| is not tiny.
  * Double-exponential (tanh-sinh) quadrature (Takahashi-Mori 1974)
    is robust for end-point singularities (k_rho = 0 logarithmic
    branch); modern tanh-sinh implementations (mpmath, scipy.integrate)
    handle the H_0^(1) branch cleanly.

## DCIM -- Discrete Complex Image Method (Chow et al. 1991)

  Y.L. Chow, J.J. Yang, D.G. Fang, G.E. Howard, IEEE Trans. Microw.
  Theory Tech. 39(3), 588-592, March 1991.

Idea: approximate the spectral kernel  tilde G(k_rho)  by a sum of
COMPLEX exponentials in k_z,

    tilde G(k_rho) ~ sum_n  a_n  exp(-b_n k_z)

each of which Sommerfeld-integrates ANALYTICALLY via the Weyl
identity to a complex spherical wave

    sum_n  a_n  exp(-j k r_n) / r_n

with  r_n = sqrt(rho^2 + (z + j b_n / k_0)^2)  (complex image
position).  Reduces a layered Green's function evaluation to a
finite sum of N ~ 5..20 image terms.  HUGE acceleration for MoM
loops where each Green's function call would otherwise need a SIP
integration.

The exponentials are fitted by GPOF (Generalised Pencil of Function,
Hua-Sarkar 1989, IEEE Trans. Antennas Propagat. 37(2), 229-234) or
by Prony's method.

## Aksun two-level DCIM (1996)

  M.I. Aksun, IEEE Trans. Microw. Theory Tech. 44(5), 651-658,
  May 1996.

Refines the original DCIM by extracting the quasi-static asymptote
AND the surface-wave poles BEFORE the exponential fit.  Produces a
much wider range of validity (near field AND far field, all
positions of z, z') than original DCIM.  Standard reference for
layered-medium MoM Green's function tables.

## Surface-wave pole extraction

For a multilayer stack with guided surface modes, tilde R^{TM/TE}
has POLES on the real k_rho axis (lossless) or just below (lossy).
These contribute LARGE RESIDUES that the exponential fit cannot
capture well.

Recipe:
  1. Locate poles by Cauchy-integral search or by direct contour
     deformation (Hua's method, also based on GPOF).
  2. Subtract pole residues from tilde G(k_rho) analytically.
  3. Fit the SMOOTH remainder with DCIM exponentials.
  4. Re-add the surface-wave residues as analytic terms in the
     final Green's function.

## Rational Function Fitting (RFF / VECTFIT)

  B. Gustavsen and A. Semlyen, IEEE Trans. Power Delivery 14(3),
  1052-1061, July 1999.

Alternative to GPOF: fit  tilde G(k_rho)  with a rational function
in k_rho directly.  Better-conditioned than Prony when the spectral
kernel has near-degenerate poles; standard in power-system EMTP
codes.  Less common in MoM literature but mathematically equivalent
for the smooth-tail part.

## Mixed-potential MoM (Michalski-Mosig 1997)

  K.A. Michalski, J.R. Mosig, IEEE Trans. Antennas Propagat. 45(3),
  508-519, March 1997.

The standard FORMULATION for layered-medium MoM: use the
mixed-potential integral equation (MPIE) instead of EFIE/MFIE, which
makes the Sommerfeld kernels DIRECTLY evaluable as scalar and vector
potentials with milder spectral behaviour.  Almost every modern
layered-MoM code uses MPIE + DCIM as the kernel evaluation strategy.

## Open-source / commercial implementations

  * Strata (Okhmatovski group, U. Manitoba): open-source layered-medium
    Green's function library with DCIM acceleration.
  * Mosig's various EPFL codes: research-grade FORTRAN/Matlab.
  * IE3D (Zeland Software, now Mentor Graphics): commercial MPIE+DCIM
    microstrip MoM (legacy, has been the de facto industry standard).
  * scikit-rf: transmission-line analogues (not Sommerfeld but uses
    the same spectral-domain plane-wave reflection coefficients).

## Workflow recommendation for radia-mcp users

If a user asks for Sommerfeld-integral evaluation in code, the
decision is:

  1. SINGLE INTERFACE, impedance plane?  -> Koh-Yook 2006 closed form
     (eq. 17), implement via the Lipschitz-Hankel + Weber library.
     This is exact for any complex eta_s and any geometry.

  2. SINGLE LAYER above PEC?  -> Direct SIP quadrature is feasible
     because the kernel is well-behaved; use Gauss-Patterson + tanh-sinh.

  3. MULTI-LAYER STACK?  -> Numerical quadrature alone is slow.
     Use ngsolve.bem with a Helmholtz kernel; if a tabulated layered
     Green's function is needed, run DCIM (Aksun 1996) offline once
     and call the resulting analytic image-sum at every MoM kernel
     evaluation.

  4. POWER-SYSTEM / TRANSMISSION-LINE analogue?  -> VECTFIT
     (Gustavsen-Semlyen).

Pure Sommerfeld-integral evaluation is RARELY a lab production need;
the more common production paths (Radia MMM/MSC, Kelvin-FEM, PEEC)
sidestep the layered-medium Green's function entirely by using
volume integrals + open-boundary truncation OR Kelvin inversion.
"""


LAB_USAGE_NOTES = r"""
# Lab-specific notes on Sommerfeld integrals

## Where Sommerfeld DOES appear in lab work (rare)

  * Antenna over earth in NDT setups (subsurface flaw detection):
    `radia_mcp.ndt` references this when the source-receiver pair
    sits above a half-space conductor.
  * Underground HF / VLF transmission line over earth: a sibling of
    the lab's litz / transmission topic, see
    `radia_mcp.litz_transmission`.
  * GPR (ground-penetrating radar) surface response: cited but not
    implemented in lab code.

## Where Sommerfeld DOES NOT appear

  * Radia MMM / MSC: volume / surface integral with Laplace kernel
    only.  No layered background; the lab uses an UNBOUNDED
    homogeneous space and lets the integrals decay naturally.
  * FEM-Kelvin: bounded computational domain via Kelvin inversion.
    Open boundary handled by inversion, not by Sommerfeld.
  * PEEC: filament + panel circuit extraction in homogeneous space.
    Layered MoM PEEC variants exist (microstrip PEEC) but are not
    lab production.
  * ngsolve.bem at HF: uses Helmholtz Green's function with free-space
    kernel, not Sommerfeld.  For high-frequency scattering off
    finite scatterers in homogeneous air this is the right answer.

## When a Sommerfeld answer IS needed -- decision

  * The question is "what does an HED at 50 m height radiate at
    1 km range over real lossy ground?"
    -> Koh-Yook 2006 closed form (SOMMERFELD_LAYERED tool).
  * The question is "what is the input impedance of a monopole over
    PEC ground"?
    -> Image theory + closed-form (Radia MMM with image symmetry,
       not Sommerfeld).
  * The question is "what does the field in layer 3 look like when
    a dipole is in layer 1 of a 5-layer dielectric stack"?
    -> ngsolve.bem with multilayer Green's function, accelerated
       by DCIM.  Lab will NOT typically build this from scratch.

## Cross-references inside radia-mcp

  * `bem_low_freq('loop_star')` -- the low-frequency stabilisation
    that PEEC uses for the homogeneous-space analogue.
  * `bem_mmm_msc('mmm_overview')` -- Radia's volume-integral
    alternative that sidesteps the Sommerfeld integral entirely.
  * `bem_h_matrix('aca')` -- the acceleration the lab uses for
    homogeneous-space MoM (NOT directly applicable to Sommerfeld
    kernels, which need DCIM instead).
  * `radia_mcp.radia_ngsolve.analytical_formulas('gauss_patterson')`
    -- nested quadrature suitable for SIP tail integration.
  * `radia_mcp.matrix_solvers.preconditioners_ams` -- companion
    krylov / preconditioner theory for the dense MoM matrices
    that arise once Sommerfeld evaluations populate them.
"""


REFERENCES = r"""
# Key references

## Source bundle absorbed (2026-05-26)

  W.C. Chew, "Sommerfeld Integral, Weyl Identity" (Lecture 35),
    in "Lectures on Electromagnetic Field Theory" lecture notes
    (Purdue University, ECE 604/618 series; revisions post-2018).
    -> Equations (35.1.1) through (35.5).  Pedagogical foundation.

  I.-S. Koh and J.-G. Yook, "Exact Closed-Form Expression of a
    Sommerfeld Integral for the Impedance Plane Problem,"
    IEEE Trans. Antennas Propagat. 54(9), 2568-2576, Sept 2006.
    DOI: 10.1109/TAP.2006.880747.
    -> Exact closed form via incomplete Weber / ILHI / ICF Poisson.

## Companion classical references (cited via Chew bibliography)

  A. Sommerfeld, "Ueber die Ausbreitung der Wellen in der drahtlosen
    Telegraphie," 1909.  Original dipole-over-half-space problem.

  W.C. Chew, "Waves and Fields in Inhomogeneous Media," IEEE Press
    1995 (also 1990).  Canonical modern monograph on layered-medium
    Green's functions; cited at almost every formula in the lecture
    notes above.

  J.A. Kong, "Electromagnetic fields due to dipole antennas over
    stratified anisotropic media," Geophysics 37(6), 985-996, 1972.
    The Ez / Hz reformulation used in the lecture notes.

  L. Brekhovskikh, "Waves in Layered Media," Academic Press 1980.

  L.B. Felsen, N. Marcuvitz, "Radiation and Scattering of Waves,"
    Wiley 1994 (also 1973).  Steepest-descent reference.

  V.A. Banos, "Dipole Radiation in the Presence of a Conducting
    Half-Space," Pergamon 1966.  The patchwork asymptotic
    developments that Koh-Yook 2006 supersedes.

## Numerical acceleration references (gap-filled by this module)

  Y.L. Chow, J.J. Yang, D.G. Fang, G.E. Howard, "A closed-form
    spatial Green's function for the thick microstrip substrate,"
    IEEE Trans. Microw. Theory Tech. 39(3), 588-592, March 1991.
    Original DCIM paper.

  M.I. Aksun, "A robust approach for the derivation of closed-form
    Green's functions," IEEE Trans. Microw. Theory Tech. 44(5),
    651-658, May 1996.  Two-level DCIM.

  K.A. Michalski, J.R. Mosig, "Multilayered media Green's functions
    in integral equation formulations," IEEE Trans. Antennas Propagat.
    45(3), 508-519, March 1997.  MPIE survey.

  Y. Hua, T.K. Sarkar, "Generalized pencil-of-function method for
    extracting poles of an EM system from its transient response,"
    IEEE Trans. Antennas Propagat. 37(2), 229-234, Feb 1989.  GPOF.

  B. Gustavsen, A. Semlyen, "Rational approximation of frequency
    domain responses by vector fitting," IEEE Trans. Power Delivery
    14(3), 1052-1061, July 1999.  VECTFIT.

  I.V. Lindell, E. Alanen, "Exact image theory for the Sommerfeld
    half-space problem," IEEE Trans. Antennas Propagat. AP-32,
    Parts I-III (Jan, Aug, Aug 1984).
"""


TOPICS = {
    "overview":          OVERVIEW,
    "weyl":              WEYL_DERIVATION,
    "weyl_derivation":   WEYL_DERIVATION,
    "sommerfeld":        SOMMERFELD_DERIVATION,
    "sommerfeld_derivation": SOMMERFELD_DERIVATION,
    "layered":           LAYERED_MEDIUM_DIPOLES,
    "layered_dipoles":   LAYERED_MEDIUM_DIPOLES,
    "ved_hed":           LAYERED_MEDIUM_DIPOLES,
    "branch_cuts":       BRANCH_CUTS_AND_SIP,
    "sip":               BRANCH_CUTS_AND_SIP,
    "riemann":           BRANCH_CUTS_AND_SIP,
    "koh_yook":          KOH_YOOK_2006_CLOSED_FORM,
    "impedance_plane":   KOH_YOOK_2006_CLOSED_FORM,
    "closed_form":       KOH_YOOK_2006_CLOSED_FORM,
    "numerical":         NUMERICAL_EVALUATION,
    "dcim":              NUMERICAL_EVALUATION,
    "gpof":              NUMERICAL_EVALUATION,
    "acceleration":      NUMERICAL_EVALUATION,
    "lab_usage":         LAB_USAGE_NOTES,
    "lab_notes":         LAB_USAGE_NOTES,
    "references":        REFERENCES,
    "bibliography":      REFERENCES,
}


def get_sommerfeld_layered_knowledge(topic: str = "all") -> str:
    """Return Sommerfeld layered-medium knowledge for the requested topic.

    Args:
        topic: one of TOPICS keys, or 'all' for full document.

    Available canonical topics (aliases also accepted):
        overview                  -- when to use, identity summary, conventions
        weyl / weyl_derivation    -- 2D plane-wave integral derivation
        sommerfeld / sommerfeld_derivation
                                  -- 1D Bessel/Hankel integral derivation
        layered / layered_dipoles / ved_hed
                                  -- VED + HED over a stratified medium
        branch_cuts / sip / riemann
                                  -- Riemann sheets + Sommerfeld contour
        koh_yook / impedance_plane / closed_form
                                  -- Koh-Yook 2006 exact closed form
        numerical / dcim / gpof / acceleration
                                  -- DCIM, GPOF, MPIE, RFF gap-filler
        lab_usage / lab_notes     -- where lab does + does NOT need this
        references / bibliography -- key papers
    """
    key = topic.strip().lower()
    if key == "all":
        return "\n\n---\n\n".join([
            OVERVIEW,
            WEYL_DERIVATION,
            SOMMERFELD_DERIVATION,
            LAYERED_MEDIUM_DIPOLES,
            BRANCH_CUTS_AND_SIP,
            KOH_YOOK_2006_CLOSED_FORM,
            NUMERICAL_EVALUATION,
            LAB_USAGE_NOTES,
            REFERENCES,
        ])
    if key not in TOPICS:
        return (f"Unknown topic '{topic}'. Available: "
                f"{', '.join(sorted(set(TOPICS)))}, or 'all'.")
    return TOPICS[key]
