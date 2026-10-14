"""Uniformly-magnetized bodies: closed-form fields, demagnetizing factors and
magnet-magnet forces for radia-ngsolve -- the analytic counterpart of the
H(div) magnetization solvers.

A body of uniform magnetization ``M`` (or remanence ``B_rem = mu0 M``) carries an
equivalent surface "magnetic charge" sigma_m = M.n, and the resulting field is
the gradient of a scalar potential.  For a handful of canonical shapes the
potential integrates in closed form, giving textbook references against which a
finite-element / boundary-element magnetostatic solve can be checked without any
meshing.  This module collects those closed forms:

    * uniformly magnetized SPHERE      -- interior uniform field + exterior dipole
                                          (Jackson Sec. 5.10; Griffiths Ex. 6.1)
    * uniformly magnetized CYLINDER    -- on-axis field of a "bar magnet"
                                          (Furlani, *Permanent Magnet ...* Ch. 3)
    * uniformly magnetized ANNULUS     -- ring magnet = outer minus inner cylinder
    * rectangular-prism DEMAG factors  -- magnetometric N_x,N_y,N_z, sum = 1
                                          (Aharoni, J. Appl. Phys. 83, 3432 (1998))
    * finite-cylinder central DEMAG    -- ballistic factor, L/D -> inf -> 0
                                          (Joseph, J. Appl. Phys. 37, 4639 (1966))
    * magnetic POLARIZABILITY of sphere-- induced moment of a permeable sphere
                                          (Jackson Sec. 5.11)
    * two coaxial magnets FORCE        -- far-field dipole-dipole attraction
    * Halbach SEGMENTATION factor      -- discrete-segment bore-field reduction
                                          (Halbach, NIM 169, 1 (1980); Mallinson)

Everything is a pure closed form (``import math`` only); no permeability except
the vacuum constant ``MU0``.  SI units throughout (metre, ampere/metre, tesla,
newton).
"""
import math

MU0 = 4.0e-7 * math.pi          # vacuum permeability [H/m]


def uniform_magnetized_sphere_field(magnetization_M, radius_a, x, y, z):
    """Field of a uniformly magnetized sphere of radius ``a``, magnetization
    ``M`` along +z (Jackson, *Classical Electrodynamics* Sec. 5.10).

    Inside (r < a) the field is **uniform**::

        H_in = -M/3  (z-hat),     B_in = mu0 (H_in + M) = (2/3) mu0 M  (z-hat).

    Outside (r > a) the sphere looks like a point dipole of moment::

        m = (4/3) pi a^3 M ,

    with B(r) = (mu0/4pi)[ 3(m.r^)r^ - m ] / r^3.

    Returns a dict with the region, the B components [T], the interior H_z [A/m]
    and the equivalent dipole moment [A m^2].

    Raises ``ValueError`` for non-positive radius.
    """
    if radius_a <= 0.0:
        raise ValueError("radius_a must be positive")
    moment = (4.0 / 3.0) * math.pi * radius_a ** 3 * magnetization_M
    r = math.sqrt(x * x + y * y + z * z)
    if r < radius_a:
        return {
            "region": "interior",
            "Bx": 0.0,
            "By": 0.0,
            "Bz": (2.0 / 3.0) * MU0 * magnetization_M,
            "H_in_z": -magnetization_M / 3.0,
            "moment": moment,
        }
    # exterior point dipole m = moment z-hat
    pref = MU0 / (4.0 * math.pi)
    mdotr = moment * z                      # m . r  (m along z)
    Bx = pref * (3.0 * mdotr * x / r ** 5)
    By = pref * (3.0 * mdotr * y / r ** 5)
    Bz = pref * (3.0 * mdotr * z / r ** 5 - moment / r ** 3)
    return {
        "region": "exterior",
        "Bx": Bx,
        "By": By,
        "Bz": Bz,
        "H_in_z": -magnetization_M / 3.0,
        "moment": moment,
    }


def uniform_magnetized_cylinder_axial_field(magnetization_M, radius_a, length_L, z):
    """On-axis axial flux density of a uniformly (axially) magnetized cylinder
    ("bar magnet") of radius ``a`` and length ``L``, centred at the origin
    (Furlani, *Permanent Magnet and Electromechanical Devices* Ch. 3).

    With magnetization ``M`` along +z, at axial coordinate ``z`` (on the axis)::

        Bz(z) = (mu0 M / 2) * [ (z+L/2)/sqrt(a^2+(z+L/2)^2)
                              - (z-L/2)/sqrt(a^2+(z-L/2)^2) ].

    The centre value (z=0) is mu0 M L / sqrt(L^2 + 4 a^2); for a long rod
    (L >> a) this tends to mu0 M.  Far away (z >> L,a) it decays as a dipole of
    moment m = pi a^2 L M.

    Returns ``{Bz, centre, moment}`` (B in tesla, moment in A m^2).

    Raises ``ValueError`` for non-positive radius or length.
    """
    if radius_a <= 0.0:
        raise ValueError("radius_a must be positive")
    if length_L <= 0.0:
        raise ValueError("length_L must be positive")
    zp = z + 0.5 * length_L
    zm = z - 0.5 * length_L
    Bz = 0.5 * MU0 * magnetization_M * (
        zp / math.sqrt(radius_a * radius_a + zp * zp)
        - zm / math.sqrt(radius_a * radius_a + zm * zm)
    )
    centre = MU0 * magnetization_M * length_L / math.sqrt(
        length_L * length_L + 4.0 * radius_a * radius_a
    )
    moment = math.pi * radius_a * radius_a * length_L * magnetization_M
    return {"Bz": Bz, "centre": centre, "moment": moment}


def uniform_magnetized_annulus_axial_field(magnetization_M, r_inner, r_outer,
                                           length_L, z):
    """On-axis axial flux density of a uniformly (axially) magnetized RING
    magnet (annulus) of inner radius ``r_inner``, outer radius ``r_outer`` and
    length ``L`` (Furlani Ch. 3).

    The ring is the superposition of a solid outer cylinder minus a coaxial
    inner cylinder of the SAME magnetization, so::

        Bz_ring(z) = Bz_cyl(r_outer, z) - Bz_cyl(r_inner, z),

    each cylinder term being the on-axis form of
    :func:`uniform_magnetized_cylinder_axial_field`.

    Returns ``{Bz}`` (tesla).

    Raises ``ValueError`` unless ``r_outer > r_inner > 0``.
    """
    if r_inner <= 0.0:
        raise ValueError("r_inner must be positive")
    if r_outer <= r_inner:
        raise ValueError("r_outer must exceed r_inner")
    outer = uniform_magnetized_cylinder_axial_field(
        magnetization_M, r_outer, length_L, z)["Bz"]
    inner = uniform_magnetized_cylinder_axial_field(
        magnetization_M, r_inner, length_L, z)["Bz"]
    return {"Bz": outer - inner}


def rectangular_prism_demag_factors(a, b, c):
    """Magnetometric demagnetizing factors (N_x, N_y, N_z) of a rectangular
    prism with full side lengths ``a`` x ``b`` x ``c`` along x, y, z
    (Aharoni, J. Appl. Phys. 83, 3432 (1998), Eq. 1).

    Aharoni's closed form gives N_z for a prism of *semi-axes* (a/2,b/2,c/2);
    using semi-axes (p,q,r) and the abbreviations below,

        pi N_z = (b^2-c^2)/(2bc) ln[(sqrt(a^2+b^2+c^2)-a)/(sqrt(a^2+b^2+c^2)+a)]
               + (a^2-c^2)/(2ac) ln[(sqrt(a^2+b^2+c^2)-b)/(sqrt(a^2+b^2+c^2)+b)]
               + (b/2c) ln[(sqrt(a^2+b^2)+a)/(sqrt(a^2+b^2)-a)]
               + (a/2c) ln[(sqrt(a^2+b^2)+b)/(sqrt(a^2+b^2)-b)]
               + (c/2a) ln[(sqrt(b^2+c^2)-b)/(sqrt(b^2+c^2)+b)]
               + (c/2b) ln[(sqrt(a^2+c^2)-a)/(sqrt(a^2+c^2)+a)]
               + 2 atan(ab/(c sqrt(a^2+b^2+c^2)))
               + (a^3+b^3-2c^3)/(3abc)
               + (a^2+b^2-2c^2)/(3abc) sqrt(a^2+b^2+c^2)
               + (c/ab)(sqrt(a^2+c^2)+sqrt(b^2+c^2))
               - ((a^2+b^2)^1.5 + (b^2+c^2)^1.5 + (c^2+a^2)^1.5)/(3abc),

    where (a,b,c) here are the full side lengths used by Aharoni's normalized
    formula (his factors depend only on the ratios).  N_x and N_y follow by
    cyclic permutation, and N_x + N_y + N_z = 1 exactly.  A cube gives
    N_x = N_y = N_z = 1/3.

    Returns ``{Nx, Ny, Nz, sum}``.

    Raises ``ValueError`` for non-positive side lengths.
    """
    if a <= 0.0 or b <= 0.0 or c <= 0.0:
        raise ValueError("side lengths a, b, c must be positive")

    def _nz(a, b, c):
        # Aharoni (1998), Eq. (1); a,b,c are full side lengths along x,y,z.
        a2, b2, c2 = a * a, b * b, c * c
        abc = math.sqrt(a2 + b2 + c2)
        ab = math.sqrt(a2 + b2)
        bc = math.sqrt(b2 + c2)
        ca = math.sqrt(c2 + a2)
        t = (
            (b2 - c2) / (2.0 * b * c) * math.log((abc - a) / (abc + a))
            + (a2 - c2) / (2.0 * a * c) * math.log((abc - b) / (abc + b))
            + b / (2.0 * c) * math.log((ab + a) / (ab - a))
            + a / (2.0 * c) * math.log((ab + b) / (ab - b))
            + c / (2.0 * a) * math.log((bc - b) / (bc + b))
            + c / (2.0 * b) * math.log((ca - a) / (ca + a))
            + 2.0 * math.atan2(a * b, c * abc)
            + (a2 * a + b2 * b - 2.0 * c2 * c) / (3.0 * a * b * c)
            + (a2 + b2 - 2.0 * c2) / (3.0 * a * b * c) * abc
            + c / (a * b) * (ca + bc)
            - (ab ** 3 + bc ** 3 + ca ** 3) / (3.0 * a * b * c)
        )
        return t / math.pi

    Nz = _nz(a, b, c)
    Nx = _nz(b, c, a)
    Ny = _nz(c, a, b)
    return {"Nx": Nx, "Ny": Ny, "Nz": Nz, "sum": Nx + Ny + Nz}


def cylinder_central_demag_factor(length_over_diameter):
    """Central (ballistic) demagnetizing factor N of a finite circular cylinder
    magnetized along its axis, as a function of the aspect ratio gamma = L/D
    (Joseph, J. Appl. Phys. 37, 4639 (1966)).

    For an axially magnetized rod of length L and diameter D the *central*
    fluxmetric factor follows from the on-axis self-field; with the equivalent
    surface-pole model the relation reduces to the simple closed form::

        N_central = 1 - gamma / sqrt(gamma^2 + 1),       gamma = L/D,

    which is exactly the demagnetizing factor read off the centre of the on-axis
    B(z) of :func:`uniform_magnetized_cylinder_axial_field`, since
    B_centre = mu0 M L/sqrt(L^2+4a^2) = mu0 M (1 - N) with gamma = L/(2a).
    Limits: gamma -> 0 (thin disc) -> N = 1; gamma -> inf (long needle) -> N = 0.

    Returns ``{N_central}`` (dimensionless, in (0,1)).

    Raises ``ValueError`` for non-positive aspect ratio.
    """
    if length_over_diameter <= 0.0:
        raise ValueError("length_over_diameter must be positive")
    g = length_over_diameter
    N = 1.0 - g / math.sqrt(g * g + 1.0)
    return {"N_central": N}


def magnetic_polarizability_sphere(mu_r, radius_a):
    """Magnetic polarizability of a linear permeable sphere of relative
    permeability ``mu_r`` and radius ``a`` in a uniform applied field H0
    (Jackson, *Classical Electrodynamics* Sec. 5.11).

    The sphere acquires a uniform interior magnetization and an exterior dipole
    moment m = beta H0, with::

        beta = 4 pi a^3 (mu_r - 1)/(mu_r + 2)   [m^3].

    The volume-averaged magnetization per applied field is
    M_avg / H0 = 3 (mu_r - 1)/(mu_r + 2) = beta / ((4/3) pi a^3).  For a vacuum
    sphere (mu_r = 1) beta = 0; for a perfect soft magnet (mu_r -> inf)
    M_avg/H0 -> 3 (and beta -> 4 pi a^3).

    Returns ``{beta, M_avg_per_H0}``.

    Raises ``ValueError`` for non-positive radius or mu_r < 0.
    """
    if radius_a <= 0.0:
        raise ValueError("radius_a must be positive")
    if mu_r < 0.0:
        raise ValueError("mu_r must be non-negative")
    factor = (mu_r - 1.0) / (mu_r + 2.0)
    beta = 4.0 * math.pi * radius_a ** 3 * factor
    return {"beta": beta, "M_avg_per_H0": 3.0 * factor}


def coaxial_magnets_axial_force(b_rem, radius_a, length_L, gap, mu_r=1.0):
    """Far-field axial attraction/repulsion between two identical coaxial
    cylindrical magnets, modelled as point dipoles (maglev dipole-dipole limit).

    Each magnet of remanence ``B_rem``, radius ``a`` and length ``L`` carries an
    equivalent moment m = (B_rem/mu0) pi a^2 L.  Two collinear (head-to-tail)
    dipoles separated centre-to-centre by s attract with::

        F = 3 mu0 m^2 / (2 pi s^4),      s = gap + L,

    i.e. the force falls off as 1/s^4 (the gradient of the 1/s^3 dipole field).
    ``mu_r`` (>= 1) rescales the effective remanence by 1/mu_r for a soft return
    path; the default mu_r = 1 is the bare permanent-magnet pair.

    Returns ``{force_N, moment, sep_s}``.

    Raises ``ValueError`` for non-positive geometry or mu_r < 1.
    """
    if radius_a <= 0.0:
        raise ValueError("radius_a must be positive")
    if length_L <= 0.0:
        raise ValueError("length_L must be positive")
    if gap < 0.0:
        raise ValueError("gap must be non-negative")
    if mu_r < 1.0:
        raise ValueError("mu_r must be >= 1")
    moment = (b_rem / MU0) * math.pi * radius_a * radius_a * length_L / mu_r
    s = gap + length_L
    force = 3.0 * MU0 * moment * moment / (2.0 * math.pi * s ** 4)
    return {"force_N": force, "moment": moment, "sep_s": s}


def coaxial_magnets_force_gap_sweep_summary(b_rem, radius_a, length_L, gaps, mu_r=1.0, tol=1.0e-12):
    """Summarize the dipole-limit force contract over a gap sweep.

    ELF/MAGIC, radia-ngsolve, and notebooks may report axial magnet force from
    different APIs or tables.  In the far-field dipole limit a same-orientation
    coaxial pair must satisfy ``F*s^4 = constant`` with
    ``s = gap + length_L``.  This helper turns that identity into a compact
    solver-independent gate before reading product force tables.
    """

    gap_values = [float(value) for value in gaps]
    if not gap_values:
        raise ValueError("gaps must not be empty")
    if any(value < 0.0 for value in gap_values):
        raise ValueError("gaps must be non-negative")
    rows = []
    for gap in gap_values:
        row = coaxial_magnets_axial_force(b_rem, radius_a, length_L, gap, mu_r=mu_r)
        invariant = row["force_N"] * row["sep_s"] ** 4
        rows.append({
            "gap_m": gap,
            "sep_s": row["sep_s"],
            "force_N": row["force_N"],
            "moment": row["moment"],
            "force_s4_invariant": invariant,
        })

    reference = rows[0]["force_s4_invariant"]
    invariant_errors = [
        abs(row["force_s4_invariant"] - reference) / max(abs(reference), 1.0e-300)
        for row in rows
    ]
    checks = {
        "gaps_non_decreasing": all(
            rows[i]["gap_m"] <= rows[i + 1]["gap_m"] + float(tol)
            for i in range(len(rows) - 1)
        ),
        "force_decreases_with_gap": all(
            rows[i]["force_N"] + float(tol) >= rows[i + 1]["force_N"]
            for i in range(len(rows) - 1)
        ),
        "force_s4_invariant_ok": max(invariant_errors) <= float(tol),
        "force_positive_for_head_to_tail_pair": all(row["force_N"] > 0.0 for row in rows),
    }
    return {
        "schema": "radia-ngsolve.coaxial-magnets-force-gap-sweep.v1",
        "b_rem_T": float(b_rem),
        "radius_a_m": float(radius_a),
        "length_L_m": float(length_L),
        "mu_r": float(mu_r),
        "tol": float(tol),
        "rows": rows,
        "reference_force_s4_invariant": reference,
        "max_force_s4_invariant_rel_error": max(invariant_errors),
        "force_ratio_first_last": rows[0]["force_N"] / rows[-1]["force_N"],
        "checks": checks,
        "status": "ok" if all(checks.values()) else "needs_attention",
    }


def halbach_segmentation_factor(n_segments):
    """Bore-field reduction factor of an N-segment Halbach dipole ring relative
    to the ideal continuously-rotating magnetization (Halbach, Nucl. Instrum.
    Methods 169, 1 (1980); Mallinson, IEEE Trans. Magn. 9, 678 (1973)).

    Replacing the ideal sinusoidal magnetization by N equal piecewise-uniform
    segments multiplies the fundamental bore field by::

        K_seg = sin(2 pi / N) / (2 pi / N) = sinc(2/N),

    the average of the (rotated) magnetization over one segment.  K_seg -> 1 as
    N -> inf (e.g. N=8 -> 0.9003, N=16 -> 0.9745).

    Returns ``{seg_factor}`` (dimensionless, in (0,1)).

    Raises ``ValueError`` for fewer than 3 segments.
    """
    if n_segments < 3:
        raise ValueError("n_segments must be >= 3")
    x = 2.0 * math.pi / n_segments
    return {"seg_factor": math.sin(x) / x}
