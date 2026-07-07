"""Analytic sphere scattering (partial-wave series).

Ported from the readable MATLAB acoustic FEM/BEM teaching lane
(matlab-acoustic-fembem, matlab_api/acoustic/*SphereScattering.m) to Python.
e^{+ikr} radiation convention; the fluid host has sound speed c = 1 and density
rho = 1, so speeds and densities passed here are ratios relative to the fluid.

  soft_sphere_scattering     <- softSphereScattering.m     (sound-soft, p = 0)
  rigid_sphere_scattering    <- rigidSphereScattering.m    (sound-hard, dp/dn = 0)
  fluid_sphere_scattering    <- fluidSphereScattering.m    (Anderson 1950 transmission)
  elastic_sphere_scattering  <- elasticSphereScattering.m  (Faran 1951 solid sphere)

Spherical functions use scipy: j_l = spherical_jn(l, x), y_l = spherical_yn(l, x),
h_l^(1) = j_l + i y_l; first derivatives via scipy's
f_l'(x) = f_{l-1}(x) - (l+1)/x f_l(x) recurrence (identical to the MATLAB
half-integer form), and the second spherical Bessel derivative from the ODE
f'' = -(2/x) f' - (1 - l(l+1)/x^2) f.

Validation (validation_test/acoustics/): the returned fields match the MATLAB
partial-wave references to ~1e-14 (machine precision) and an independent
ngsolve.bem numerical Helmholtz BEM solve (sound-soft sphere) to ~2e-5.
"""

from __future__ import annotations

import numpy as np
from scipy.special import spherical_jn, spherical_yn, eval_legendre


# --- spherical special functions (e^{+ikr} outgoing Hankel) ------------------- #
def _jn(l, x):
    return spherical_jn(l, x)


def _jn_d(l, x):
    return spherical_jn(l, x, derivative=True)


def _jn_dd(l, x):
    # spherical Bessel ODE: f'' = -(2/x) f' - (1 - l(l+1)/x^2) f
    x = np.asarray(x, dtype=float)
    return -(2.0 / x) * _jn_d(l, x) - (1.0 - l * (l + 1) / x**2) * _jn(l, x)


def _h1(l, x):
    return spherical_jn(l, x) + 1j * spherical_yn(l, x)


def _h1_d(l, x):
    return spherical_jn(l, x, derivative=True) + 1j * spherical_yn(l, x, derivative=True)


def _prep(points, radius, allow_interior=False):
    pts = np.asarray(points, dtype=float).reshape(-1, 3)
    r = np.sqrt((pts**2).sum(axis=1))
    if not allow_interior and np.any(r < radius * (1 - 1e-9)):
        raise ValueError("evaluation points must lie on or outside the sphere r >= R")
    r_safe = np.maximum(r, 1e-30)
    costh = pts[:, 2] / r_safe
    return pts, r, r_safe, costh


# --- sound-soft (p = 0) ------------------------------------------------------- #
def soft_sphere_scattering(wavenumber, radius, points, terms=None):
    """Plane-wave exp(ikz) scattering by a sound-soft (p = 0) sphere.

    Returns a dict with ``scattered``, ``incident`` and ``total`` complex
    pressures (length N) at ``points`` (N x 3, exterior r >= R), plus the mode
    count and the last-mode magnitude ``truncation_tail``.
    """
    k = float(wavenumber)
    R = float(radius)
    pts, r, _, costh = _prep(points, R)
    if terms is None:
        terms = int(np.ceil(k * R)) + 12
    scattered = np.zeros(r.shape, dtype=complex)
    last = np.zeros(r.shape, dtype=complex)
    for l in range(terms + 1):
        Pl = eval_legendre(l, costh)
        a = -(1j**l) * (2 * l + 1) * _jn(l, k * R) / _h1(l, k * R)
        last = a * _h1(l, k * r) * Pl
        scattered += last
    incident = np.exp(1j * k * pts[:, 2])
    return {
        "kind": "soft_sphere_plane_wave_scattering_series",
        "wavenumber": k, "radius": R, "terms": int(terms),
        "truncation_tail": float(np.max(np.abs(last))),
        "scattered": scattered, "incident": incident, "total": incident + scattered,
    }


# --- sound-hard (dp/dn = 0) --------------------------------------------------- #
def rigid_sphere_scattering(wavenumber, radius, points, terms=1):
    """Plane-wave exp(ikz) scattering by a sound-hard (dp/dn = 0) sphere."""
    k = float(wavenumber)
    R = float(radius)
    pts, r, _, costh = _prep(points, R)
    x0 = k * R
    L = max(int(terms), int(np.ceil(k * max(R, float(r.max())))) + 12)
    scattered = np.zeros(r.shape, dtype=complex)
    last = np.zeros(r.shape, dtype=complex)
    for l in range(L + 1):
        Pl = eval_legendre(l, costh)
        A = -(1j**l) * (2 * l + 1) * _jn_d(l, x0) / _h1_d(l, x0)
        last = A * _h1(l, k * r) * Pl
        scattered += last
    incident = np.exp(1j * k * pts[:, 2])
    return {
        "kind": "rigid_sphere_plane_wave_scattering_series",
        "wavenumber": k, "radius": R, "terms": int(L),
        "truncation_tail": float(np.max(np.abs(last))),
        "scattered": scattered, "incident": incident, "total": incident + scattered,
    }


# --- penetrable fluid sphere (Anderson 1950 transmission) --------------------- #
def fluid_sphere_scattering(wavenumber, radius, points,
                            interior_wavenumber=None, density_ratio=1.0, terms=None):
    """Anderson (1950) plane-wave transmission by a penetrable fluid sphere.

    ``interior_wavenumber`` k1 and ``density_ratio`` rho1/rho0 set the contrast;
    k1 = k0 and density_ratio = 1 make the sphere acoustically invisible
    (total == incident).  The returned ``total`` is the interior series for
    r <= R and incident + scattered outside; ``inside_mask`` flags interior rows.
    """
    k0 = float(wavenumber)
    k1 = float(interior_wavenumber) if interior_wavenumber is not None else k0
    rhor = float(density_ratio)
    R = float(radius)
    pts, r, r_safe, costh = _prep(points, R, allow_interior=True)
    r_max = max(R, float(r.max()))
    L = max(0 if terms is None else int(terms),
            int(np.ceil(max(k0 * r_max, k1 * R))) + 12)
    inside = r <= R * (1 + 1e-12)
    x0, x1 = k0 * R, k1 * R
    total = np.zeros(r.shape, dtype=complex)
    last = np.zeros(r.shape, dtype=complex)
    for l in range(L + 1):
        Pl = eval_legendre(l, costh)
        a_inc = (1j**l) * (2 * l + 1)
        j0R, h0R, j1R = _jn(l, x0), _h1(l, x0), _jn(l, x1)
        dj0R, dh0R, dj1R = _jn_d(l, x0), _h1_d(l, x0), _jn_d(l, x1)
        # analytic elimination through the interior log-derivative beta (the naive
        # 2x2 solve is ill-conditioned at high l and pollutes the invisible case)
        beta = (k1 / rhor) * dj1R / j1R
        A = -a_inc * (k0 * dj0R - beta * j0R) / (k0 * dh0R - beta * h0R)
        B = (a_inc * j0R + A * h0R) / j1R
        mode = np.zeros(r.shape, dtype=complex)
        mode[inside] = B * _jn(l, k1 * r_safe[inside]) * Pl[inside]
        mode[~inside] = (a_inc * _jn(l, k0 * r_safe[~inside])
                         + A * _h1(l, k0 * r_safe[~inside])) * Pl[~inside]
        total += mode
        last = mode
    incident = np.exp(1j * k0 * pts[:, 2])
    return {
        "kind": "fluid_sphere_transmission_scattering_series",
        "wavenumber": k0, "interior_wavenumber": k1, "density_ratio": rhor,
        "radius": R, "terms": int(L), "truncation_tail": float(np.max(np.abs(last))),
        "incident": incident, "total": total, "inside_mask": inside,
    }


# --- elastic solid sphere (Faran 1951) ---------------------------------------- #
def _elastic_coeff(n, k, a, rhoF, omega, kL, cT, rhoS, lam, mu):
    """Scattered coefficient c_l from the r = a boundary conditions (Faran)."""
    x, xl = k * a, kL * a
    fluid_fac = k / (rhoF * omega**2)
    if cT == 0:
        # fluid interior (no shear): 2x2 for [c; A], sigma_rr = -lam kL^2 A j_l(xl)
        M = np.array([
            [fluid_fac * _h1_d(n, x), -kL * _jn_d(n, xl)],
            [_h1(n, x),               -lam * kL**2 * _jn(n, xl)],
        ], dtype=complex)
        rhs = np.array([-fluid_fac * _jn_d(n, x), -_jn(n, x)], dtype=complex)
        return np.linalg.solve(M, rhs)[0]

    kT = omega / cT
    xt = kT * a
    Ur_A = kL * _jn_d(n, xl)
    Ur_B = n * (n + 1) / a * _jn(n, xt)
    dUr_A = kL**2 * _jn_dd(n, xl)
    dUr_B = n * (n + 1) * (kT * _jn_d(n, xt) / a - _jn(n, xt) / a**2)
    # radial normal stress sigma_rr = -lam kL^2 phi + 2 mu d(u_r)/dr  (sign: sigma_rr = -p)
    Srr_A = -lam * kL**2 * _jn(n, xl) + 2 * mu * dUr_A
    Srr_B = 2 * mu * dUr_B
    Va_A = _jn(n, xl) / a
    Va_B = _jn(n, xt) / a + kT * _jn_d(n, xt)
    Vp_A = -_jn(n, xl) / a**2 + kL * _jn_d(n, xl) / a
    Vp_B = -_jn(n, xt) / a**2 + kT * _jn_d(n, xt) / a + kT**2 * _jn_dd(n, xt)
    # zero shear stress for an inviscid fluid: sigma_rtheta = mu(u_r/a + V' - V/a) = 0
    Srt_A = mu * (Ur_A / a + Vp_A - Va_A / a)
    Srt_B = mu * (Ur_B / a + Vp_B - Va_B / a)
    M = np.array([
        [fluid_fac * _h1_d(n, x), -Ur_A, -Ur_B],
        [_h1(n, x),                Srr_A, Srr_B],
        [0.0,                      Srt_A, Srt_B],
    ], dtype=complex)
    rhs = np.array([-fluid_fac * _jn_d(n, x), -_jn(n, x), 0.0], dtype=complex)
    return np.linalg.solve(M, rhs)[0]


def elastic_sphere_scattering(wavenumber, radius, points,
                              longitudinal_speed=2.0, shear_speed=1.0,
                              density_ratio=1.5, terms=0):
    """Faran (1951) plane-wave scattering by a solid elastic sphere in a fluid.

    Speeds are ratios to the fluid (c_fluid = 1): ``longitudinal_speed`` cL/c,
    ``shear_speed`` cT/c (0 = fluid), ``density_ratio`` rho_solid/rho_fluid.
    Exterior scattered-field reference (h_l(kr) is singular inside a solid), so
    ``points`` must be exterior (r >= R).  Limits: shear_speed -> 0 recovers the
    Anderson fluid sphere, very stiff recovers the rigid sphere.
    """
    k = float(wavenumber)
    a = float(radius)
    omega = k                       # fluid c = 1
    cL, cT = float(longitudinal_speed), float(shear_speed)
    rhoS, rhoF = float(density_ratio), 1.0
    mu = rhoS * cT**2
    lam = rhoS * (cL**2 - 2 * cT**2)
    kL = omega / cL
    L = int(terms) if terms and terms > 0 else int(np.ceil(k * a)) + 10
    # exterior scattered-field series (h_l(kr) singular inside): require r >= R
    pts, r, _, costh = _prep(points, a)
    coeff = [_elastic_coeff(l, k, a, rhoF, omega, kL, cT, rhoS, lam, mu) for l in range(L + 1)]
    scattered = np.zeros(r.shape, dtype=complex)
    last = np.zeros(r.shape, dtype=complex)
    for l in range(L + 1):
        Pl = eval_legendre(l, costh)
        last = (1j**l) * (2 * l + 1) * coeff[l] * _h1(l, k * r) * Pl
        scattered += last
    incident = np.exp(1j * k * pts[:, 2])
    return {
        "kind": "elastic_solid_sphere_faran_scattering_series",
        "wavenumber": k, "radius": a, "longitudinal_speed": cL, "shear_speed": cT,
        "density_ratio": rhoS, "terms": int(L), "truncation_tail": float(np.max(np.abs(last))),
        "incident": incident, "scattered": scattered, "total": incident + scattered,
    }
