"""Analytical formulas knowledge (mcp-server-radia-ngsolve, public tool).

Closed-form analytical expressions in :mod:`radia.analytical_formulas`,
sourced from the IEE Japan review

    若尾真治, 五十嵐一, 藤原耕二, 野口聡, 松尾哲司, 亀有昭久,
    "Useful Formulas of Analytical Integration in Electromagnetic Field
    Computations (Part 1..5)", SA / RM technical meeting, 2002-2004.

The package collects nine reference formula modules. They are the
trusted-baseline layer that the rest of Radia (HDiv-VIM, PEEC, FEM
panels, ngsolve.bem) is checked against. **For new analyses the first
question to ask is "is there a closed-form here that I can sanity-check
against?".** Cross-validation patterns are listed in the
``validation_use_cases`` topic.
"""


_OVERVIEW = """
# radia.analytical_formulas overview

Source: ``src/radia/analytical_formulas/`` (9 modules, ~30 functions).
Tests: ``tests/analytical_formulas/`` (145 + ; pytest in < 1 s).
Examples: ``docs/analytical_formulas/analytical_formulas.ipynb`` (12 runnable demos).
Docs: ``docs/analytical_formulas.md`` (PDF -> code cross-reference).

## Module index

Group B -- fields and global quantities

| Module | Purpose | PDF reference |
|--------|---------|---------------|
| ellipsoid       | Rotational ellipsoid demag factor + internal field + torque | Part 5 §5, eq 38-44 |
| ac_locus        | B_max / B_min of an AC vector phasor's time locus | Part 5 §4, eq 29-37 |
| shielding       | Static shielding factor of a magnetic cyl/sphere shell | Part 1 §5, eq 23-24 |
| rect_magnet_2d  | 2D rectangular bar A_z, B_x, B_y for arbitrary M | Part 2 §2, eq 2-3 |
| plate_eddy      | Thin rectangular-plate eddy current | Part 1 §6.1, eq 26-27 |

Group C -- coil geometries and numerical utilities

| Module | Purpose | PDF reference |
|--------|---------|---------------|
| solenoid_central   | Fabri F(alpha,beta) central + closed-form axial field | Part 4 §4, eq 26-27 |
| three_phase_line   | Triangle / planar / helical line field (multipole + far-field) | Part 4 §5, Part 5 §3 |
| elliptic_integrals | K(k), E(k) Hastings polynomial approximation | Part 3 §3, Tables 1-2 |
| gauss_legendre     | Gauss-Legendre nodes / weights to n=24 + integrate helpers | Part 3 §4, Table 3 |

Group D -- Part 6 / 8 / 9 extensions (eddy / AC / numerics / averaging)

| Module | Purpose | PDF reference |
|--------|---------|---------------|
| plate_eddy.plate_eddy_dissipation | Total Joule dissipation in thin-plate eddy current (numerical integration of analytic J series) | Part 6 §3 |
| shielding.shielding_factor_*_thin_ac | Frequency-dependent thin-shell AC shielding factor (complex S, sphere/cylinder) | Part 8 §2, eq 4-5 |
| shielding.spherical_shell_internal_field | All four uniform fields (interior, shell wall, M, image dipole) of magnetic spherical shell | Part 6 §2, eq 5 |
| conductor_impedance | Skin depth, planar surface impedance Z_s=(1+j)/(sigma delta), full Bessel cylindrical-conductor AC impedance | Part 6 §4-§5 |
| adaptive_quadrature | Gauss-Patterson nested rules n=0..3 + adaptive integrate with node-reuse caching | Part 9 §2, Table 1 |
| cuboid_average_field.average_B_in_box, average_demag_tensor | Spatial average of B over a target box from uniform M source box. Closed-form 64-corner C++ kernel (G1/G2 antiderivatives, sympy-derived in-house, ~40 us/call, 817x faster than Gauss-Legendre); method="numerical" fallback kept for cross-checks. v4.22.0+. | Part 6 §7 |
| induction_heating.cylinder_axial_eddy_loss (+ small_ka, thin_skin asymptotes) | Time-averaged Joule loss per unit length of a conductive cylinder inside a long axial-current solenoid (canonical induction-furnace workpiece). Closed-form via Kelvin functions ber/bei. Two limits: small-ka (Faraday eddy P=πσω²μ²H₀²a⁴/16) and thin-skin (P=πa·H₀²·Re(Z_s)). Validated in-house; replaces aborted Stafl 1967 §4.4 OCR transcription. v4.23.0+. | Smythe §11.07, Landau-Lifshitz §59, Jackson §5.18 |

## Conventions

* All quantities in **SI (m, A, A/m, A/m**2, T)**.
* Magnetisation ``M`` is in A/m (Radia convention).
  ``J_polarisation = mu_0 M`` in Tesla.
* Three-phase phasors use **peak** amplitude (eq 33 of Part 4):
  ``I_peak = sqrt(2) * I_rms`` for a sinusoidal RMS input.
* Time convention: ``b(t) = Re(B_phasor exp(j omega t))``.

## Quick API surface

```python
from radia.analytical_formulas import (
    # ellipsoid
    demag_factor_rotational, ellipsoid_internal_field, ellipsoid_torque,
    # ac vector locus
    ac_locus_axes, ac_locus_axes_batch,
    # shielding
    shielding_factor_cylinder, shielding_factor_sphere,
    # 2D rectangular bar
    rect_magnet_2d_A, rect_magnet_2d_B,
    # thin-plate eddy
    plate_eddy_T, plate_eddy_J,
    # solenoid central + axial
    fabri_F, solenoid_central_field, solenoid_axial_field,
    # three-phase line
    triangle_far_field_amplitude, planar_far_field_amplitude,
    helical_far_field_amplitude, field_xy, vector_potential_z,
    # elliptic integrals
    K_hastings_4, E_hastings_4,
    # quadrature
    gauss_legendre_integrate, gauss_legendre_integrate_2d,
)
```
"""


_ELLIPSOID = """
# ellipsoid — rotational ellipsoid demag factor + torque

Geometry: semi-axes ``(a, a, c)``, polar (symmetry) axis ``c``.

* prolate (cigar):  c > a
* oblate  (disk):   c < a
* sphere:           c = a

## Polar demag factor (PDF eq 39-41)

```python
N_z = demag_factor_prolate(c, a)   # for c >= a
N_z = demag_factor_oblate(c, a)    # for c <= a
N_x, N_y, N_z = demag_factor_rotational(c, a)
# N_x + N_y + N_z = 1, exactly
```

* prolate: ``N_z = (1 - e**2)/e**3 (atanh e - e)``,
  ``e = sqrt(1 - (a/c)**2)``
* oblate:  ``N_z = (1 + e**2)/e**3 (e - atan e)``,
  ``e = sqrt((a/c)**2 - 1)``
* near-sphere: Taylor branch at ``|c/a - 1| < 1e-4`` :
  ``N_z ~ 1/3 - 2 eps/15 + 8 eps**2/105``,
  ``eps = (c/a)**2 - 1``.

Reference values, validated against Osborn 1945 to 1e-13:

| c / a | N_z |
|-------|-----|
| 0.01  | 0.9844897069128692 |
| 0.1   | 0.8608042765278006 |
| 0.5   | 0.5272002825625699 |
| 1.0   | 0.3333333333333333 |
| 2.0   | 0.17356399753396423 |
| 5.0   | 0.0558209698024552 |
| 10.0  | 0.02028588030156383 |

## Internal field in a uniform external field (eq 38)

```python
H_in = ellipsoid_internal_field(H_0, chi_r, N)
# H_in = H_0 / (1 + chi_r * N)
```

## Torque on a tilted ellipsoid (eq 44)

```python
T_z = ellipsoid_torque(H_0, alpha, chi_r, a, c)
```

* ``T_z < 0`` is restoring (prolate aligns with field).
* ``T_z > 0`` is repulsive (oblate aligns perpendicular).
* ``T_z = 0`` at ``alpha = 0, pi/2, pi`` (sin(2 alpha) = 0).
* Scales as ``H_0**2 * V``.

## Use it for

* Validate HDiv-VIM magnetisation of a meshed ellipsoid against
  the analytical internal field (mu_r effective vs N_z prediction).
* Quick sanity check of permanent-magnet far-field demag corrections.
* Predict the torque on a soft-iron insert before running a 3D solve.
"""


_AC_LOCUS = """
# ac_locus — major / minor axis of an AC vector locus

For a complex 3-vector phasor ``B`` whose time-domain version is
``b(t) = Re(B exp(j omega t))``, the tip of ``b(t)`` traces a planar
ellipse. Closed form (PDF eq 37):

    s2 = |B_x|**2 + |B_y|**2 + |B_z|**2          (real)
    z2 = B_x**2 + B_y**2 + B_z**2                (complex; B . B)
    |B|_max = sqrt((s2 + |z2|) / 2)
    |B|_min = sqrt((s2 - |z2|) / 2)

```python
B_max, B_min = ac_locus_axes(B)                 # B shape (3,)
B_max, B_min = ac_locus_axes_batch(B_grid)      # B shape (..., 3)
```

## Special cases

* Linear:  ``B_min = 0``    (the locus collapses to a line).
* Circular: ``B_max = B_min`` (z2 = 0).
* Elliptical: anything in between.

## Use it for

* Peak ``B`` over an AC cycle for hysteresis-loss correlations
  (peak ``|B|`` is what core-loss models like Steinmetz consume,
  not RMS).
* Peak ``J`` for induction heating uniformity constraints.
* Distinguishing linear from rotating polarisation in eddy-current
  post-processing without loop-time-stepping.
"""


_SHIELDING = """
# shielding — static shielding factor of a magnetic shell

Shielding factor ``S = |H_inside| / |H_applied|`` (small = good shield).

```python
S_cyl = shielding_factor_cylinder(a, b, mu_r)   # eq 23
S_sph = shielding_factor_sphere(a, b, mu_r)     # eq 24
```

* a = inner radius, b = outer radius.
* mu_r = relative permeability of the shell material.
* ``S = 1`` at ``mu_r = 1`` (no shielding); decreases monotonically
  with ``mu_r`` and with shell thickness ``t = b - a``.

## High-mu_r limit

```
S_cyl ~ 4 b / (mu_r (b - a))     for b - a << a, mu_r >> 1
S_sph ~ S_asym / (1 + S_asym),
S_sph,asym = 9 / (2 mu_r (1 - (a/b)**3))
```

## Use it for

* Ballpark the wall thickness needed for a magnetic shield around an
  experiment.
* Validate a 3D HDiv-VIM solve of a high-mu cylinder / sphere shell
  by comparing the interior field to the analytical S.
"""


_RECT_MAGNET_2D = """
# rect_magnet_2d — 2D rectangular bar field

Cross-section ``2a x 2b``, infinitely long in z, uniform 2D
magnetisation ``M = (M_x, M_y)`` in A/m.

```python
A_z = rect_magnet_2d_A(point, a, b, M=(Mx, My))
B   = rect_magnet_2d_B(point, a, b, M=(Mx, My))   # shape (..., 2)
```

The returned ``A_z`` is in T m, ``B`` in Tesla.

## Useful limits

* Long thin slab (a << b), M along the long axis: ``B_x ~ mu_0 M_x``
  (no demag along the long direction).
* Same slab, M perpendicular: ``B ~ 0`` inside (demag = 1).
* Far field: 2D dipole law ``|B| ~ mu_0 m / (2 pi r**2)`` with
  ``m = M_x * 4 a b`` per unit length.

## Use it for

* 2D quick-check of 3D Radia ``CuboidMagnet`` results in the limit
  of a long bar (verified by example
  the cross_validation_3d_vs_2d section of ``docs/analytical_formulas/analytical_formulas.ipynb``).
* Sanity-check a 2D NGSolve magnetostatic solve of a rectangular
  permanent magnet.
"""


_PLATE_EDDY = """
# plate_eddy — thin rectangular plate eddy current

Plate of in-plane half-extents ``a x b``, thickness ``d``, in a
spatially-uniform perpendicular ``B_z(t)``. Slow-field limit:

    d / (mu sigma B_dot * min(a, b) / B) << 1   (PDF eq 25)

Then the in-plane current is divergence-free and is the curl of a
single z-component vector potential ``T_z`` (PDF eq 26-27):

```python
T = plate_eddy_T(point, a, b, sigma, Bdot, n_terms=200)   # A/m
J = plate_eddy_J(point, a, b, sigma, Bdot, n_terms=200)   # (..., 2) A/m**2
```

* T_z = 0 on x = +-a exactly (cos((2n+1) pi/2) = 0).
* J . n = 0 on the plate edges (the closed-form ``-sigma Bdot x / 4``
  in J_y cancels the cosh-cosh series at y = +-b).
* J = curl(T e_z): J_x = dT/dy, J_y = -dT/dx (verified to FD
  precision in the test suite).

## Use it for

* TEAM-style benchmark for thin-plate eddy current solvers.
* Quick estimate of eddy-current losses in a thin metallic shield
  exposed to an AC perpendicular field.
* Sanity-check a 3D shell-element FEM eddy current run.
"""


_SOLENOID = """
# solenoid_central — rectangular cross-section solenoid

Geometry: inner radius ``a_1``, outer radius ``a_2``, half-length ``b``,
azimuthal current density ``J``.

```python
F   = fabri_F(alpha, beta)                           # alpha = a2/a1, beta = b/a1
B_0 = solenoid_central_field(a1, a2, b, J)           # B_z(0, 0)
B_z = solenoid_axial_field(z, a1, a2, b, J)          # B_z(0, z), closed form
```

## Limits

* Long coil (``b -> infinity``): ``B_0 -> mu_0 J (a_2 - a_1)``.
* Thin pancake (``b << a_1``):    ``B_0 ~ mu_0 a_1 J beta ln(alpha)``.

## Use it for

* Quick design estimates for finite solenoids before running 3D Radia.
* Cross-validate a 3D coil-builder geometry: the on-axis field
  produced by a stacked set of ``radia.analytical_magnet.CurrentLoop``
  rings must converge to ``solenoid_axial_field`` at midpoint-rule
  ``O(h**2)`` rate (verified in
  the cross_validation_solenoid_currentloop section of ``docs/analytical_formulas/analytical_formulas.ipynb``;
  4x error reduction per grid halving).
"""


_THREE_PHASE_LINE = """
# three_phase_line — three-phase line magnetic field

Phasor convention: peak amplitude (eq 33 of Part 4)
``I_peak,k = I_peak exp(- j 2 k pi / 3)``. For RMS input ``I_rms``,
pass ``I_peak = sqrt(2) I_rms``.

## Direct evaluator (any geometry)

```python
A_z = vector_potential_z(positions, currents, observation)
B   = field_xy(positions, currents, observation)   # (..., 2)
```

* ``positions``: shape ``(n_lines, 2)`` Cartesian (x, y).
* ``currents``: shape ``(n_lines,)`` complex phasors.
* Use this for arbitrary line configurations, then post-process via
  :func:`ac_locus_axes` to get ``|B|_max, |B|_min``.

## Far-field amplitudes (specific arrangements)

| arrangement | spacing param | leading order | polarisation |
|-------------|----------------|---------------|--------------|
| triangle (3 lines on circle a) | a | 3 mu_0 a I_peak / (4 pi r**2) | circular |
| planar   (3 lines, spacing d)  | d | sqrt(3) mu_0 d I_peak / (2 pi r**2) | linear |
| helical (a << p << r)          | a, p | (3 pi mu_0 a I_peak / (4 sqrt(p r))) exp(-2 pi r / p) | rotating |

```python
B_tri = triangle_far_field_amplitude(r, I_peak, a)
B_pln = planar_far_field_amplitude(r, I_peak, d)
B_hel = helical_far_field_amplitude(r, I_peak, a, p)
```

## Hexagon caveat

The 1/r**3 hexagon asymptote in PDF eq 44 requires a SPECIFIC
phase-on-vertex assignment. The simple builder pair
``hexagon_positions + balanced_six_phase_currents`` is two
anti-phase 3-phase sets at opposite vertices, which gives a
1/r**2 dipole-dominant far-field (not 1/r**3).  The 1/r**3
formula is intentionally not exposed; use ``field_xy`` with the
arrangement appropriate to your problem.

## Use it for

* Predict the magnetic field at a sensitive instrument location near
  a power line.
* Sanity-check the polarisation (circular for triangle, linear for
  planar) when interpreting AC measurement data.
* Reference for far-field decay rates of helical superconducting
  cables (exponential).
"""


_ELLIPTIC_INTEGRALS = """
# elliptic_integrals — Hastings polynomial K(k), E(k)

Polynomial approximations of the complete elliptic integrals of the
first and second kind (PDF eq 13 + Tables 1-2):

    K(k) ~ sum (a_i + b_i ln(1/x)) x**i,  x = 1 - k**2
    E(k) ~ sum (c_i + d_i ln(1/x)) x**i

```python
K = K_hastings_2(k)   # < 4e-5 absolute
K = K_hastings_4(k)   # < 2e-8 absolute
E = E_hastings_2(k), E_hastings_4(k)
```

For machine-precision values, use ``scipy.special.ellipk(k**2)`` (note
scipy parameter is ``m = k**2``). The polynomial form is preserved
here as a documented reference and a no-scipy fallback.

## Use it for

* Closed-form ring-current field via :func:`fabri_F` etc., when full
  scipy is unavailable.
* Reference for verifying a custom K(k)/E(k) implementation.
"""


_GAUSS_LEGENDRE = """
# gauss_legendre — Gauss-Legendre quadrature up to n=24

```python
x, w = gauss_legendre_nodes_weights(n)        # 1 <= n <= 24
val = gauss_legendre_integrate(f, a, b, n=16)
val = gauss_legendre_integrate_2d(f, a, b, c, d, nx=16, ny=16)
```

* Nodes and weights computed via ``numpy.polynomial.legendre.leggauss``
  to machine precision.
* Exact for polynomials of degree ``2 n - 1``.
* Spectral convergence on smooth integrands (sin, exp, ...).

## Use it for

* Numerical integration in Radia panel / utility code where scipy.integrate
  is overkill.
* Reference rule for cross-checking custom quadrature implementations.
"""


_VALIDATION_USE_CASES = """
# Validation use cases — "I have X, sanity-check it with Y"

The analytical_formulas layer is the FIRST QUESTION to ask before a
numerical solve: "is there a closed-form here?"  This topic is the
practical mapping from "what kind of analysis am I running" to
"which analytical formula is the trusted reference".

## "I solved a permanent-magnet ellipsoid in 3D Radia"

* Compare ``M`` inside the ellipsoid to
  ``H_in = H_0 / (1 + chi_r * N_z)`` from
  :func:`ellipsoid_internal_field`.
* For the torque, compare to :func:`ellipsoid_torque`.

## "I solved an AC eddy current problem and want peak |B|"

* Take the complex phasor ``B = (Bx, By, Bz)`` at the point of interest
  and call :func:`ac_locus_axes` -> ``B_max, B_min``.
* Use ``B_max`` for hysteresis loss; ``B_max == B_min`` indicates
  circular polarisation.

## "I solved a magnetic shield in 3D and want to confirm S"

* Compare the interior |H| / |H_applied| to
  :func:`shielding_factor_cylinder` or :func:`shielding_factor_sphere`.

## "I solved a 2D infinite rectangular permanent magnet"

* Compare the field to :func:`rect_magnet_2d_B` (analytical exact).

## "I solved a thin-plate eddy current problem"

* Compare J(x, y, t = 0) to :func:`plate_eddy_J` and the streamlines
  to :func:`plate_eddy_T` contours.

## "I built a 3D solenoid coil model"

* Compare on-axis B_z to :func:`solenoid_axial_field`.
* Sample reference, ``a_1 = 50 mm, a_2 = 100 mm, b = 100 mm,
  J = 1e7 A/m**2``: ``B_0 = 0.502858 T`` (within 1e-5 of a 20 x 80
  CurrentLoop stack).

## "I have a 3-phase line and want field at a sensitive point"

* If the conductor cross-section is small relative to ``r``, use
  :func:`triangle_far_field_amplitude` or
  :func:`planar_far_field_amplitude` (1/r**2 decay).
* For arbitrary positions / phasors, use :func:`field_xy` and
  :func:`ac_locus_axes`.

## "I want to compute K(k), E(k) without scipy"

* Use :func:`K_hastings_4`, :func:`E_hastings_4` (~1e-8 absolute).

## "I need quadrature in pure Python without scipy.integrate"

* Use :func:`gauss_legendre_integrate` (n=16 default; spectral
  convergence on smooth integrands).

## Cross-validation tests already in the repo

* ``tests/analytical_formulas/test_cross_validation.py``:
  - 2D rectangular bar -> 3D CuboidMagnet long-bar limit (O(1/L_z**2))
  - high-mu_r spherical shielding asymptote
  - solenoid_axial_field -> stacked CurrentLoop midpoint rule (O(h**2),
    convergence ratio = 4 per halving)
* ``docs/analytical_formulas/analytical_formulas.ipynb`` (cross_validation_3d_vs_2d + cross_validation_solenoid_currentloop sections)

## NGSolve (FEM) cross-validation suite -- validation/force/validate_force_xval.py

The radia-ngsolve A-form / A-Phi solvers are checked against these
closed-forms (``python validation/force/validate_force_xval.py``; each is a real 3D solve). Besides the
TEAM problems (6/7/13/20/21) and COMSOL-recorded values, the analytic
cases include:
  - magnetised / permeable sphere demag: B_in = 3 mu_r/(mu_r+2) B0
    (mu_r=10 AND mu_r=1000 extreme-contrast, scattered-field A-form, ~1%)
  - permeable cylinder transverse: B_in = 2 mu_r/(mu_r+1) B0
  - coaxial-loop MUTUAL INDUCTANCE via the vector-potential flux linkage
    M = 2 pi a2 A_phi(a2,z2)/I1  vs Maxwell elliptic M (~2.4%)
  - coaxial-loop AXIAL FORCE via the Lorentz relation
    F_z = -2 pi a2 I2 B_r1(a2,z2)  vs I1 I2 dM/dz (elliptic) (~1%)
  - eddy diffusion closed-forms: slab B/cosh(kh), solid sphere B0 p/sinh(p),
    solid cylinder B0/I0(kR)  (kappa = sqrt(j w sigma mu0))
  - Helmholtz centre 8/(5 sqrt5) mu0 NI/R; annular-solenoid Biot-Savart;
    PM-sphere exterior dipole mu0/(4pi) 2m/z^3.
  - anti-Helmholtz central axial GRADIENT
    G = 3 mu0 NI a^2 s/(a^2+s^2)^(5/2) (opposed coils, finite-diff of B_z)
  - current-loop OFF-AXIS field (both B_r and B_z) vs the exact
    elliptic-integral loop field (Simpson/Jackson), not just on-axis.
These reuse the shielding / elliptic_integrals / demag formulas above as
the trusted reference, extending the cross-validation from radia's own
analytic+CurrentLoop layer to the FEM solver.

These two patterns -- 2D bar limit and stacked-CurrentLoop reference --
are reproducible templates. When you add a new closed-form to
:mod:`radia.analytical_formulas`, mirror them in a new
``cross_validation_*.py`` example so the cross-check is itself a
runnable artefact, not a one-off in a chat.
"""


_TOPICS = {
    "overview": _OVERVIEW,
    "ellipsoid": _ELLIPSOID,
    "ac_locus": _AC_LOCUS,
    "shielding": _SHIELDING,
    "rect_magnet_2d": _RECT_MAGNET_2D,
    "plate_eddy": _PLATE_EDDY,
    "solenoid_central": _SOLENOID,
    "three_phase_line": _THREE_PHASE_LINE,
    "elliptic_integrals": _ELLIPTIC_INTEGRALS,
    "gauss_legendre": _GAUSS_LEGENDRE,
    "validation_use_cases": _VALIDATION_USE_CASES,
}


def get_analytical_formulas_documentation(topic: str = "all") -> str:
    """Return analytical_formulas knowledge for a single topic or ``all``.

    Parameters
    ----------
    topic : str
        One of: ``"all"``, ``"overview"``, ``"ellipsoid"``,
        ``"ac_locus"``, ``"shielding"``, ``"rect_magnet_2d"``,
        ``"plate_eddy"``, ``"solenoid_central"``, ``"three_phase_line"``,
        ``"elliptic_integrals"``, ``"gauss_legendre"``,
        ``"validation_use_cases"``.

    Returns
    -------
    str
        Markdown-formatted documentation. ``topic = "all"`` concatenates
        every topic in declaration order.
    """
    if topic == "all":
        return "\n\n---\n\n".join(_TOPICS.values())
    if topic not in _TOPICS:
        return (
            f"Unknown topic '{topic}'. "
            f"Available: {', '.join(['all'] + list(_TOPICS.keys()))}"
        )
    return _TOPICS[topic]
