"""Closed-form force between two uniformly magnetised cuboids (Akoun-Yonnet).

Two rectangular magnets with edges parallel to the axes, both magnetised
along z (polarisations J1, J2 in tesla, sign gives the direction).  With the
surface-charge model the force on magnet 2 is exact:

    F = J1 J2 / (4 pi mu0) * sum_{i,j,k,l,p,q in {0,1}} (-1)^(i+j+k+l+p+q) phi(u, v, w)

    u = alpha + (-1)^j A - (-1)^i a,  v = beta + (-1)^l B - (-1)^k b,
    w = gamma + (-1)^q C - (-1)^p c,  r = sqrt(u^2 + v^2 + w^2)

    phi_x = (v^2 - w^2)/2 ln(r - u) + u v ln(r - v) + v w atan(u v / (r w)) + r u / 2
    phi_y = (u^2 - w^2)/2 ln(r - v) + u v ln(r - u) + u w atan(u v / (r w)) + r v / 2
    phi_z = -u w ln(r - u) - v w ln(r - v) + u v atan(u v / (r w)) - r w

where (a, b, c) and (A, B, C) are the half sizes of magnets 1 and 2 and
(alpha, beta, gamma) is the centre of magnet 2 relative to magnet 1.

Reference: G. Akoun and J.-P. Yonnet, "3D analytical calculation of the forces
exerted between two cuboidal magnets", IEEE Trans. Magn. 20(5):1962-1964, 1984.

Rigid, uniform magnetisation (mu_r = 1) is assumed.  The implementation was
checked against direct Gauss-Legendre integration of the same charge model and
against the far-field dipole force -3 mu0 m1 m2 / (2 pi d^4).
"""
from __future__ import annotations

import math

import numpy as np

MU0 = 4.0e-7 * math.pi

MAX_REL_ROUNDING = 1.0e-6

__all__ = ["cuboid_magnet_force", "MAX_REL_ROUNDING"]


def _xlog(coef: float, arg: float) -> float:
    # c ln(x) with the limit c ln(x) -> 0 where the coefficient vanishes on an edge.
    return 0.0 if coef == 0.0 else coef * math.log(arg)


def _r_minus(r: float, a: float, b2: float) -> float:
    """r - a without cancellation: equals b2 / (r + a) when a > 0 (b2 = r^2 - a^2)."""
    return b2 / (r + a) if a > 0.0 else r - a


def cuboid_magnet_force(polarization_1_t: float, size_1_m, polarization_2_t: float,
                        size_2_m, center_offset_m) -> np.ndarray:
    """Force [N] on magnet 2 from magnet 1; both z-magnetised cuboids.

    Args:
        polarization_1_t, polarization_2_t: J = mu0 M along +z (negative = -z), tesla.
        size_1_m, size_2_m: full edge lengths (x, y, z), metres.
        center_offset_m: centre of magnet 2 minus centre of magnet 1, metres.

    Raises:
        ValueError: non-positive sizes, non-finite input, overlapping volumes, or a
            separation so large that double-precision cancellation in the corner
            sum could exceed MAX_REL_ROUNDING of |F| (a-posteriori bound).
    """
    s1 = np.asarray(size_1_m, float); s2 = np.asarray(size_2_m, float)
    d = np.asarray(center_offset_m, float)
    if s1.shape != (3,) or s2.shape != (3,) or d.shape != (3,):
        raise ValueError("sizes and offset must each have three components")
    vals = np.concatenate([s1, s2, d, [polarization_1_t, polarization_2_t]])
    if not np.all(np.isfinite(vals)):
        raise ValueError("inputs must be finite")
    if np.any(s1 <= 0) or np.any(s2 <= 0):
        raise ValueError("edge lengths must be positive")
    if np.all(np.abs(d) < 0.5 * (s1 + s2)):
        raise ValueError("magnet volumes overlap; the closed form needs disjoint magnets")
    a, b, c = 0.5 * s1
    A, B, C = 0.5 * s2
    al, be, ga = d
    terms = ([], [], [])
    for i in (0, 1):
        for j in (0, 1):
            u = al + (-1) ** j * A - (-1) ** i * a
            for k in (0, 1):
                for l in (0, 1):
                    v = be + (-1) ** l * B - (-1) ** k * b
                    for p in (0, 1):
                        for q in (0, 1):
                            w = ga + (-1) ** q * C - (-1) ** p * c
                            r = math.sqrt(u * u + v * v + w * w)
                            sgn = (-1) ** (i + j + k + l + p + q)
                            at = math.atan(u * v / (r * w)) if w != 0.0 else 0.0
                            ru = _r_minus(r, u, v * v + w * w)
                            rv = _r_minus(r, v, u * u + w * w)
                            fx = (_xlog(0.5 * (v * v - w * w), ru) + _xlog(u * v, rv)
                                  + v * w * at + 0.5 * r * u)
                            fy = (_xlog(0.5 * (u * u - w * w), rv) + _xlog(u * v, ru)
                                  + u * w * at + 0.5 * r * v)
                            fz = -_xlog(u * w, ru) - _xlog(v * w, rv) + u * v * at - r * w
                            terms[0].append(sgn * fx); terms[1].append(sgn * fy); terms[2].append(sgn * fz)
    F = np.array([math.fsum(t) for t in terms])
    # The alternating corner sum cancels as the separation grows (terms ~ d^2,
    # result ~ d^-4).  Bound the rounding error from the largest term and refuse
    # results whose bound exceeds MAX_REL_ROUNDING of |F| (use a dipole model there).
    bound = 64.0 * 8.0 * np.finfo(float).eps * max(abs(x) for t in terms for x in t)
    norm = float(np.linalg.norm(F))
    if not bound <= MAX_REL_ROUNDING * norm:
        raise ValueError(f"separation too large for the closed form in double precision "
                         f"(rounding bound {bound / max(norm, 1e-300):.1e} of |F| > {MAX_REL_ROUNDING:.0e}); "
                         f"use the far-field dipole force instead")
    return polarization_1_t * polarization_2_t / (4.0 * math.pi * MU0) * F
