"""Waveguide / cavity cutoff modes for radia-ngsolve -- the 2D Helmholtz
(transverse Laplacian) EIGENVALUE problem, the wave-physics counterpart of the
static BVP solvers.

A hollow metallic waveguide supports propagating modes whose transverse pattern
solves the scalar Helmholtz eigenproblem on the cross-section:

    -nabla_t^2 psi = k_c^2 psi ,

with the boundary condition set by the mode family (PEC walls):

    * TM modes : E_z = 0 on the wall  ->  DIRICHLET Laplacian
    * TE modes : dH_z/dn = 0 on the wall  ->  NEUMANN  Laplacian

The eigenvalues k_c^2 are the squared cutoff wavenumbers; below the cutoff
frequency f_c = c k_c/(2 pi) the mode is evanescent (does not propagate). The
SAME eigenproblem gives the TM resonances of a 2D cavity and the modes of a
vibrating membrane (drum) -- only the physical reading of k_c changes.

The Neumann problem always has a trivial constant (k_c = 0) eigenmode (no field);
:func:`helmholtz_cutoff_wavenumbers_2d` discards it and returns the lowest
PHYSICAL cutoff wavenumbers. Validated against the exact rectangular-guide
spectrum (examples/comsol_class/waveguide_cutoff.py).
"""
import math

from ngsolve import (H1, HCurl, BilinearForm, GridFunction, grad, curl, dx, ArnoldiSolver)

C0 = 299792458.0          # speed of light in vacuum [m/s]
MU0 = 4.0e-7 * math.pi     # vacuum permeability [H/m]  (eta0 = MU0*C0 = 376.73 Ohm)


def rectangular_waveguide_cutoff(a, b, m, n):
    """Exact cutoff FREQUENCY [Hz] of the TE_mn / TM_mn mode of a rectangular metallic
    waveguide of inner width ``a`` and height ``b`` [m]:

        f_c,mn = (c/2) sqrt((m/a)^2 + (n/b)^2),

    i.e. k_c = sqrt((m pi/a)^2 + (n pi/b)^2). Allowed indices: TE needs m,n >= 0 (not both 0,
    dominant TE10 -> f_c = c/(2a)); TM needs m,n >= 1 (lowest TM11). The empty band between the
    dominant TE10 and the next mode sets the single-mode bandwidth (e.g. WR-90: 6.56 -> 13.1 GHz).
    """
    return 0.5 * C0 * math.hypot(m / a, n / b)


def cutoff_frequency(kc, c=C0):
    """Cutoff frequency f_c = c k_c/(2 pi) [Hz] from a cutoff wavenumber k_c [1/m]."""
    return c * kc / (2.0 * math.pi)


def circular_waveguide_cutoff(radius, mode, m, n):
    """Exact cutoff FREQUENCY [Hz] of the TE_mn / TM_mn mode of a CIRCULAR metallic waveguide of
    inner ``radius`` a -- the Bessel-function counterpart of the rectangular
    :func:`rectangular_waveguide_cutoff`:

        TM_mn: k_c = j_mn / a    (j_mn  = n-th positive zero of J_m),
        TE_mn: k_c = j'_mn / a   (j'_mn = n-th positive zero of J'_m),    f_c = c k_c/(2 pi).

    The dominant mode is TE11 (j'_11 = 1.8412 -> f_c = 0.293 c/a); the next is TM01 (j_01 = 2.4048).
    Circular guides have the famous degeneracy TE_0n / TM_1n (because j'_0n = j_1n). ``mode`` is
    'TE' or 'TM'; m = azimuthal index (>=0), n = radial index (>=1)."""
    from scipy.special import jn_zeros, jnp_zeros
    if mode.upper() == "TE":
        z = jnp_zeros(m, n)[n - 1]
    elif mode.upper() == "TM":
        z = jn_zeros(m, n)[n - 1]
    else:
        raise ValueError("mode must be 'TE' or 'TM'")
    return C0 * z / (2.0 * math.pi * radius)


def waveguide_dispersion(frequency, fc, c=C0):
    """Dispersion of a hollow metallic waveguide mode ABOVE cutoff -- the propagation sequel to the
    cutoff eigenvalue (:func:`rectangular_waveguide_cutoff` / :func:`circular_waveguide_cutoff` /
    :func:`helmholtz_cutoff_wavenumbers_2d`). A mode of cutoff frequency ``fc`` [Hz] driven at
    ``frequency`` f > fc propagates with axial constant

        beta = (2 pi/c) sqrt(f^2 - fc^2) = k0 sqrt(1 - (fc/f)^2),    k0 = 2 pi f/c,

    so the GUIDE WAVELENGTH lambda_g = 2 pi/beta = lambda0 / sqrt(1-(fc/f)^2) exceeds the free-space
    lambda0 = c/f and obeys  1/lambda_g^2 = 1/lambda0^2 - 1/lambda_c^2. The PHASE velocity
    v_p = omega/beta = c/sqrt(1-(fc/f)^2) > c carries no energy; the GROUP (energy/signal) velocity
    v_g = d omega/d beta = c sqrt(1-(fc/f)^2) < c does, and the two satisfy the reciprocal identity

        v_p * v_g = c^2 .

    At cutoff (f -> fc+) beta -> 0, lambda_g -> inf, v_p -> inf, v_g -> 0 (standing wave, no axial
    transport); for f >> fc both velocities -> c (TEM-like). Below cutoff the mode is evanescent --
    see :func:`waveguide_evanescent_attenuation`.

    Returns a dict: ``beta`` [1/m], ``lambda_g`` [m], ``v_phase`` & ``v_group`` [m/s],
    ``fc_over_f``. Raises ValueError at or below cutoff (use the evanescent branch there)."""
    if frequency <= fc:
        raise ValueError("frequency must exceed the cutoff fc (mode is evanescent below it)")
    ratio = fc / frequency
    s = math.sqrt(1.0 - ratio * ratio)
    beta = (2.0 * math.pi * frequency / c) * s
    return {"beta": beta, "lambda_g": 2.0 * math.pi / beta,
            "v_phase": c / s, "v_group": c * s, "fc_over_f": ratio}


def guide_wavelength(frequency, fc, c=C0):
    """Guide wavelength lambda_g = lambda0/sqrt(1-(fc/f)^2) [m] (always > the free-space lambda0
    = c/f). The transverse standing-wave pattern stretches the axial period; component lengths
    (slots, irises, lambda_g/4 transformers) are set by lambda_g, not lambda0."""
    return waveguide_dispersion(frequency, fc, c)["lambda_g"]


def waveguide_evanescent_attenuation(frequency, fc, c=C0):
    """BELOW cutoff (f < fc) the mode does NOT propagate; its amplitude decays as exp(-alpha z) with

        alpha = (2 pi/c) sqrt(fc^2 - f^2)   [Np/m]  (purely reactive -- lossless 'cutoff' attenuator).

    The waveguide is a high-pass filter: alpha -> 0 as f -> fc-, and alpha -> 2 pi fc/c = k_c (the
    cutoff wavenumber) as f -> 0. Used for cutoff attenuators and below-cutoff isolation. Raises
    ValueError at or above fc (use :func:`waveguide_dispersion` there)."""
    if frequency >= fc:
        raise ValueError("frequency must be below the cutoff fc (use waveguide_dispersion above it)")
    return (2.0 * math.pi / c) * math.sqrt(fc * fc - frequency * frequency)


def helmholtz_cutoff_wavenumbers_2d(mesh, n_modes, bc="neumann", wall="wall",
                                    order=3, shift=-1.0):
    """Lowest cutoff wavenumbers k_c [1/m] of a waveguide cross-section by solving the 2D
    Helmholtz eigenproblem  -nabla_t^2 psi = k_c^2 psi  on ``mesh``.

    Generalised-eigenvalue solve (stiffness ``grad.grad`` vs mass) via NGSolve's
    :func:`ArnoldiSolver` with a small NEGATIVE shift, so the shifted operator (A - shift*M =
    A + |shift|*M) is SPD and never singular -- the eigenvalues nearest the shift are then the
    algebraically smallest. Works for ANY cross-section (rectangular, ridged, circular,
    L-shaped); only the geometry/mesh changes.

    Args:
        mesh    : 2D NGSolve mesh of the guide cross-section.
        n_modes : number of (physical) cutoff wavenumbers to return, ascending.
        bc      : 'neumann' for TE modes (dH_z/dn = 0), 'dirichlet' for TM modes (E_z = 0).
        wall    : boundary name of the metal wall (used as Dirichlet for the TM case).
        order   : H1 polynomial order (default 3 -- cutoffs converge fast).
        shift    : ArnoldiSolver shift (default -1.0; keep negative/below the spectrum).

    Returns the ascending list of k_c [1/m]; the trivial Neumann constant mode (k_c ~ 0) is
    dropped. Convert to frequency with :func:`cutoff_frequency`.
    """
    fes = H1(mesh, order=order, dirichlet=(wall if bc == "dirichlet" else ""))
    u, v = fes.TnT()
    A = BilinearForm(grad(u) * grad(v) * dx).Assemble()
    M = BilinearForm(u * v * dx).Assemble()
    nev = n_modes + (4 if bc == "neumann" else 2)          # pad for the dropped ~0 mode + spares
    gf = GridFunction(fes, multidim=nev)
    lams = ArnoldiSolver(A.mat, M.mat, fes.FreeDofs(), list(gf.vecs), shift=shift)
    kc2 = sorted(l.real for l in lams)
    kc2 = [x for x in kc2 if x > 1e-6 * max(kc2)]          # drop the constant (Neumann) null mode
    return [math.sqrt(x) for x in kc2[:n_modes]]


def rectangular_cavity_frequency(a, b, d, m, n, p):
    """Exact resonant FREQUENCY [Hz] of the TE_mnp / TM_mnp mode of a rectangular metallic CAVITY
    (a closed PEC box a x b x d):

        f_mnp = (c/2) sqrt((m/a)^2 + (n/b)^2 + (p/d)^2).

    The 3-D / closed-box counterpart of the (open-ended) waveguide cutoff
    :func:`rectangular_waveguide_cutoff` -- a third index p appears for the standing wave along the
    box length. The dominant mode (for a > d > b) is TE101. Used for microwave cavity filters,
    resonators, and accelerator cavities."""
    return 0.5 * C0 * math.sqrt((m / a) ** 2 + (n / b) ** 2 + (p / d) ** 2)


def maxwell_cavity_modes_3d(mesh, n_modes, shift, order=2, pec="pec"):
    """Lowest resonant frequencies [Hz] of a 3-D PEC cavity by the FULL-WAVE MAXWELL eigenproblem

        curl curl E = (omega/c)^2 E ,   n x E = 0 on the PEC walls,

    on an ``HCurl`` (edge / Nedelec) space -- the genuine vector-electromagnetic cavity (not the
    scalar Helmholtz box). The curl-curl operator has a LARGE gradient kernel (curl grad = 0, the
    spurious zero modes); they are suppressed by giving ``ArnoldiSolver`` a ``shift`` near the
    target k^2 (the shift-invert then amplifies the PHYSICAL resonances and starves the far-away
    null space). Estimate ``shift`` from the expected fundamental, e.g.
    ``(pi/a)**2 + (pi/d)**2`` for TE101.

    Args:
        mesh    : 3-D NGSolve mesh of the cavity; PEC walls named ``pec``.
        n_modes : number of resonant frequencies to return, ascending.
        shift    : ArnoldiSolver shift ~ the squared wavenumber k^2 of interest [1/m^2].
        order   : HCurl order (default 2).
        pec     : boundary name of the metal walls (tangential E = 0).

    Returns the ascending list of resonant frequencies [Hz] (the residual near-zero gradient modes
    are dropped). Validated against the exact box spectrum (examples/comsol_class/cavity_resonance.py)."""
    fes = HCurl(mesh, order=order, dirichlet=pec)
    u, v = fes.TnT()
    A = BilinearForm(curl(u) * curl(v) * dx).Assemble()
    M = BilinearForm(u * v * dx).Assemble()
    gf = GridFunction(fes, multidim=n_modes + 4)
    lams = ArnoldiSolver(A.mat, M.mat, fes.FreeDofs(), list(gf.vecs), shift=shift)
    k2 = sorted(l.real for l in lams if l.real > 1e-3 * shift)     # drop the gradient null space
    return [C0 * math.sqrt(x) / (2.0 * math.pi) for x in k2[:n_modes]]


def skin_depth(frequency, sigma, mu_r=1.0):
    """Skin depth delta = sqrt(2/(omega mu sigma)) [m] -- the 1/e penetration of an AC field into a
    conductor (mu = mu0 mu_r). Sets the surface resistance of cavity / waveguide walls and lines; at
    10 GHz copper delta ~ 0.66 um."""
    omega = 2.0 * math.pi * frequency
    return math.sqrt(2.0 / (omega * MU0 * mu_r * sigma))


def surface_resistance(frequency, sigma, mu_r=1.0):
    """Surface resistance R_s = sqrt(omega mu/(2 sigma)) = 1/(sigma*delta) [Ohm] of a good conductor --
    the real part of the surface impedance. The time-average wall loss per unit area is
    (R_s/2)|H_tan|^2, so R_s sets the cavity/line loss; it grows as sqrt(frequency)."""
    omega = 2.0 * math.pi * frequency
    return math.sqrt(omega * MU0 * mu_r / (2.0 * sigma))


def rectangular_cavity_q(a, b, d, sigma, mu_r=1.0):
    """Wall-loss (unloaded) quality factor Q of the dominant TE101 mode of a rectangular cavity
    (a x b x d; a along x, d along z, b the height) with conducting walls -- the LOSS sequel to the
    lossless resonance :func:`rectangular_cavity_frequency` / :func:`maxwell_cavity_modes_3d` (#59).
    Q = omega * U / P_wall, with U the stored energy and P_wall the wall dissipation; for TE101

        Q = (k a d)^3 b eta / (2 pi^2 R_s (2 a^3 b + 2 b d^3 + a^3 d + a d^3)),

    k = omega/c, eta = mu0 c the free-space impedance, R_s = :func:`surface_resistance`. Q ~
    (volume/surface)/delta ~ sqrt(sigma): a better conductor (smaller surface resistance / skin depth)
    gives a sharper resonance -- the bridge from the wave family to the skin depth (#45/#55).
    Validated against the numerically-integrated TE101 field energy and wall loss
    (examples/comsol_class/cavity_quality_factor.py)."""
    f = rectangular_cavity_frequency(a, b, d, 1, 0, 1)
    k = 2.0 * math.pi * f / C0
    Rs = surface_resistance(f, sigma, mu_r)
    eta = MU0 * C0
    num = (k * a * d) ** 3 * b * eta
    den = 2.0 * math.pi ** 2 * Rs * (2 * a ** 3 * b + 2 * b * d ** 3 + a ** 3 * d + a * d ** 3)
    return num / den
