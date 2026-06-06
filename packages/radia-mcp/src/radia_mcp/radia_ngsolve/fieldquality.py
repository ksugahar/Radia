"""Magnetic FIELD-QUALITY (multipole) analysis -- the post-processing an
accelerator-magnet designer actually runs on a 2D field: decompose B on a
reference circle into normal/skew multipoles (b_n, a_n), then read off the main
harmonic and the field errors (allowed + unallowed harmonics, in 1e-4 units).

Convention (CERN / "European", the one used in beam dynamics)::

    B_y(z) + i B_x(z) = sum_{n>=1} (b_n + i a_n) (z / R_ref)^{n-1},   z = x + i y

so n = 1 is the dipole, n = 2 the quadrupole, n = 3 the sextupole, ...  ``b_n``
are the NORMAL, ``a_n`` the SKEW coefficients, each with the units of B [T] (the
field the n-th pole produces at the reference radius ``R_ref``). The
US/"American" convention shifts the index by one (their n-1 is our n); we stick
to the European index throughout.

The unambiguous validation reference is a single straight line current: a current
``I`` at complex position ``z0`` produces ``B_y + i B_x = (mu0 I)/(2 pi) / (z - z0)``,
whose expansion about the origin (|z| < |z0|) is exactly

    b_n + i a_n = -(mu0 I)/(2 pi) * R_ref^{n-1} / z0^n            (:func:`line_current_multipoles`)

-- a closed form (no FEM, no convention ambiguity) that :func:`multipole_coefficients`
must reproduce. Superposing it over several wires gives the exact free-space
multipole content of any current arrangement (e.g. a 4-wire quadrupole), the
analytic check the FEM example is held to.
"""
import cmath
import math

MU0 = 4.0e-7 * math.pi


def _sample_circle(B, R_ref, center, n_samples):
    """Sample the planar field ``B`` at ``n_samples`` points on the circle of
    radius ``R_ref`` about ``center``. ``B`` is either a python callable
    ``phi -> (Bx, By)`` OR a 2-component NGSolve CoefficientFunction (then
    ``mesh`` must be passed to :func:`multipole_coefficients`). Returns the list
    of complex samples  f_k = By + i Bx  at phi_k = 2 pi k / n_samples.
    """
    cx, cy = center
    out = []
    for k in range(n_samples):
        phi = 2.0 * math.pi * k / n_samples
        bx, by = B(phi)
        out.append(complex(by, bx))           # f = By + i Bx
    return out


def multipole_coefficients(B, R_ref, center=(0.0, 0.0), n_max=12,
                           n_samples=None, mesh=None):
    """Decompose a 2D magnetic field into normal/skew multipoles on the circle
    ``|z - center| = R_ref``.

    Parameters
    ----------
    B : callable or NGSolve CoefficientFunction
        * ``callable``  -- ``phi -> (Bx, By)`` giving the field at the circle
          point ``center + R_ref (cos phi, sin phi)`` (radians). Lets you feed a
          closed-form field with NO mesh (used by the tests).
        * 2-component CF -- the FEM flux density ``CF((Bx, By))``; pass the
          ``mesh`` it lives on. Sampled via ``B(mesh(px, py))`` (real part taken,
          so a magnetostatic real field is fine).
    R_ref : reference radius [m] (must lie in the current-free bore; harmonics
        are only meaningful where there are no sources between ``center`` and it).
    center : (x, y) expansion centre [m].
    n_max : highest harmonic to return (n = 1 .. n_max).
    n_samples : circle samples (DFT length); default ``max(64, 8*n_max)`` and
        bumped to even. Must exceed ``2*n_max`` to resolve n_max without aliasing.
    mesh : NGSolve mesh, REQUIRED when ``B`` is a CoefficientFunction.

    Returns
    -------
    dict with 1-indexed lists (index ``n`` is harmonic ``n``; index 0 is unused
    and the constant term is dropped -- a magnetic field has no n=0):
        ``b``  : normal coefficients b_n [T] (len n_max+1)
        ``a``  : skew   coefficients a_n [T]
        ``C``  : magnitudes |b_n + i a_n| [T]
        ``phase`` : arg(b_n + i a_n) [rad]
        ``R_ref``, ``n_max``, ``n_samples``, ``center``.

    The k-th complex sample is ``f_k = By + i Bx``; since
    ``f(phi) = sum_n (b_n + i a_n) e^{i (n-1) phi}`` the coefficient is the
    ``(n-1)``-th DFT bin: ``b_n + i a_n = (1/N) sum_k f_k e^{-i (n-1) phi_k}``.
    """
    if n_samples is None:
        n_samples = max(64, 8 * n_max)
    if n_samples % 2:
        n_samples += 1
    if n_samples <= 2 * n_max:
        raise ValueError(f"n_samples={n_samples} too small for n_max={n_max}; "
                         f"need n_samples > 2*n_max for an alias-free DFT.")

    if mesh is not None:                       # B is an NGSolve CF on `mesh`
        cx, cy = center

        def sampler(phi, _B=B, _m=mesh, _cx=cx, _cy=cy, _R=R_ref):
            mp = _m(_cx + _R * math.cos(phi), _cy + _R * math.sin(phi))
            val = _B(mp)
            bx = getattr(val[0], "real", val[0])
            by = getattr(val[1], "real", val[1])
            return float(bx), float(by)
        f = _sample_circle(sampler, R_ref, center, n_samples)
    else:
        if not callable(B):
            raise TypeError("B must be a callable phi->(Bx,By) or a CF with mesh=")
        f = _sample_circle(B, R_ref, center, n_samples)

    b = [0.0] * (n_max + 1)
    a = [0.0] * (n_max + 1)
    C = [0.0] * (n_max + 1)
    phase = [0.0] * (n_max + 1)
    N = n_samples
    for n in range(1, n_max + 1):
        m = n - 1                              # DFT bin
        acc = 0j
        for k in range(N):
            phi = 2.0 * math.pi * k / N
            acc += f[k] * cmath.exp(-1j * m * phi)
        coeff = acc / N
        b[n] = coeff.real
        a[n] = coeff.imag
        C[n] = abs(coeff)
        phase[n] = cmath.phase(coeff)
    return {"b": b, "a": a, "C": C, "phase": phase, "R_ref": R_ref,
            "n_max": n_max, "n_samples": N, "center": tuple(center)}


def line_current_multipoles(I, r0, phi0, R_ref, n_max=12):
    """Exact normal/skew multipoles of a single straight line current ``I`` [A]
    at radius ``r0`` [m], angle ``phi0`` [rad], expanded about the origin
    (valid for ``R_ref < r0``):

        b_n + i a_n = -(mu0 I)/(2 pi) * R_ref^{n-1} / z0^n,   z0 = r0 e^{i phi0}.

    Returns the same dict shape as :func:`multipole_coefficients` (the reference
    it is validated against). Superpose (sum the b/a lists) for many wires.
    """
    z0 = cmath.rect(r0, phi0)
    b = [0.0] * (n_max + 1)
    a = [0.0] * (n_max + 1)
    C = [0.0] * (n_max + 1)
    phase = [0.0] * (n_max + 1)
    for n in range(1, n_max + 1):
        coeff = -(MU0 * I) / (2.0 * math.pi) * (R_ref ** (n - 1)) / (z0 ** n)
        b[n] = coeff.real
        a[n] = coeff.imag
        C[n] = abs(coeff)
        phase[n] = cmath.phase(coeff)
    return {"b": b, "a": a, "C": C, "phase": phase, "R_ref": R_ref,
            "n_max": n_max, "center": (0.0, 0.0)}


def superpose_multipoles(*coeff_dicts):
    """Sum several multipole dicts (same ``R_ref``/``n_max``) coefficient-wise --
    the field-quality of a multi-wire / multi-source magnet is the superposition
    of each source's :func:`line_current_multipoles`."""
    ref = coeff_dicts[0]
    n_max = ref["n_max"]
    b = [0.0] * (n_max + 1)
    a = [0.0] * (n_max + 1)
    for d in coeff_dicts:
        if d["n_max"] != n_max or abs(d["R_ref"] - ref["R_ref"]) > 1e-15:
            raise ValueError("superpose_multipoles: mismatched R_ref / n_max")
        for n in range(1, n_max + 1):
            b[n] += d["b"][n]
            a[n] += d["a"][n]
    C = [abs(complex(b[n], a[n])) for n in range(n_max + 1)]
    phase = [cmath.phase(complex(b[n], a[n])) if C[n] else 0.0
             for n in range(n_max + 1)]
    return {"b": b, "a": a, "C": C, "phase": phase, "R_ref": ref["R_ref"],
            "n_max": n_max, "center": ref.get("center", (0.0, 0.0))}


def field_errors(coeffs, main_n=None, units=1e4):
    """Field-quality table: each harmonic's magnitude RELATIVE to the main one,
    in ``units`` (default 1e4 -> the accelerator "units of 1e-4" / 'parts in
    10^4'). ``main_n`` defaults to the largest |C_n|.

    Returns ``dict(main_n, main_C, rel_b, rel_a, rel_C)`` where ``rel_*[n]`` =
    ``coeff_n / |C_main| * units``. By construction ``rel_C[main_n] == units``.
    A well-built magnet has |rel_C[n]| of a few units (a few 1e-4) for n != main.
    """
    C = coeffs["C"]
    n_max = coeffs["n_max"]
    if main_n is None:
        main_n = max(range(1, n_max + 1), key=lambda n: C[n])
    main_C = C[main_n]
    if main_C == 0:
        raise ValueError("main harmonic magnitude is zero")
    rel_b = [coeffs["b"][n] / main_C * units for n in range(n_max + 1)]
    rel_a = [coeffs["a"][n] / main_C * units for n in range(n_max + 1)]
    rel_C = [C[n] / main_C * units for n in range(n_max + 1)]
    return {"main_n": main_n, "main_C": main_C, "rel_b": rel_b,
            "rel_a": rel_a, "rel_C": rel_C, "units": units}
