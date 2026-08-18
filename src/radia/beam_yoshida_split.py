"""Yoshida split-operator tracking on a CanonicalHCurl chain.

The original EarlyTimes integrator (Ishi 2016 deck, Yoshida factorization
of Tokyo Astronomical Observatory lineage) revived with the field passed
ENTIRELY through :class:`radia.beam_canonical_hcurl.CanonicalHCurlChain`:
``a_y``, covariant ``a_s``, the curvature ``htilde(s)``, and the exact-H2
generator integral ``G = int_0^y a_y d eta`` are all closed-form chain
polynomial evaluations.

Paraxial Hamiltonian (deck convention, on-momentum, eB-rho normalized,
``a_x = 0`` axial gauge; ``atilde = a / (B rho)``):

    H = 1/2 px^2 + 1/2 (py - atilde_y)^2 + U(x, y, s) [+ ps]
    U = -htilde(s) x - atilde_s_cov(x, y, s)

For a pure dipole ``atilde_s_cov = -htilde x - htilde^2 x^2 / 2`` so
``U = htilde^2 x^2 / 2`` -- the deck's harmonic term emerges exactly.  The
split is the deck's palindrome (its ``p_s`` end caps make every inner
factor act at the STEP MIDPOINT ``s`` automatically):

    S2(sigma) = ps/2 . H1/4 . H3/2 . H1/4 . H2 . H1/4 . H3/2 . H1/4 . ps/2

with H1 = 1/2 px^2 + 1/2 htilde^2 x^2 (exact rotation, htilde frozen at
the step's midpoint s), H3 = U - 1/2 htilde^2 x^2 (exact coordinate kick),
and H2 = 1/2 (py - atilde_y)^2 solved EXACTLY: ``v = py - atilde_y`` is
conserved (x, s frozen), ``y`` drifts linearly, and the momentum kicks are
the x-derivative of the y-antiderivative ``G`` between the end heights --
the closed form of the deck's ``f = -int a_y dy'`` conjugation.

Every factor is an exact symplectic map at ANY amplitude and step; the
only approximations are the symmetric splitting (2nd order in the step;
``order=4`` composes Yoshida's w-coefficients for 4th order with a
negative middle step) and the paraxial Hamiltonian itself -- comparing
against the exact-sqrt canonical A-RK on the same chain isolates and
quantifies that paraxial gap for the first time.

On-momentum 4D state ``(x, px, y, py)`` (deck convention has no delta);
the ``ps`` ledger is bookkeeping only (nothing depends on ``ps``) and is
not tracked.
"""

from __future__ import annotations

import numpy as np

__all__ = ["track_yoshida_split", "YOSHIDA4_W1", "YOSHIDA4_W0"]

YOSHIDA4_W1 = 1.0 / (2.0 - 2.0 ** (1.0 / 3.0))
YOSHIDA4_W0 = 1.0 - 2.0 * YOSHIDA4_W1


def _h1_rotation(state, htilde, xi):
    """Exact flow of ``H1 = px^2/2 + htilde^2 x^2/2`` over parameter xi."""
    x, px = state[0], state[1]
    angle = xi * htilde
    cos = np.cos(angle)
    # sin(xi h)/h and h sin(xi h) are even in h; both are h->0 safe.
    sinc = xi * np.sinc(angle / np.pi)
    state[0] = cos * x + sinc * px
    state[1] = -htilde * htilde * sinc * x + cos * px


def _h3_kick(chain, rigidity, state, htilde, s_value, tau):
    """Exact kick of ``H3 = U - htilde^2 x^2/2`` (coordinates frozen)."""
    x = np.array([state[0]])
    y = np.array([state[2]])
    s = np.array([s_value])
    a, gradient = chain.vector_potential_and_gradient_frame(x, y, s)
    das_dx = gradient[0, 2, 0] / rigidity
    das_dy = gradient[0, 2, 1] / rigidity
    du3_dx = -htilde - das_dx - htilde * htilde * state[0]
    du3_dy = -das_dy
    state[1] -= tau * du3_dx
    state[3] -= tau * du3_dy


def _h2_exact(chain, rigidity, state, s_value, sigma):
    """Exact flow of ``H2 = (py - atilde_y)^2/2`` (x, s frozen).

    ``v = py - atilde_y`` is conserved, ``y`` drifts by ``sigma v``, and
    ``px`` receives ``d/dx`` of the y-antiderivative G between the end
    heights; ``py`` follows from the conservation law.
    """
    x = np.array([state[0]])
    s = np.array([s_value])
    y0 = state[2]
    ay0 = chain.vector_potential_frame(x, np.array([y0]), s)[0, 1] / rigidity
    v = state[3] - ay0
    y1 = y0 + sigma * v
    g0, g0_x = chain.ay_y_antiderivative_frame(x, np.array([y0]), s)
    g1, g1_x = chain.ay_y_antiderivative_frame(x, np.array([y1]), s)
    ay1 = chain.vector_potential_frame(x, np.array([y1]), s)[0, 1] / rigidity
    state[1] += (g1_x[0] - g0_x[0]) / rigidity
    state[2] = y1
    state[3] = v + ay1


def _s2_step(chain, rigidity, state, s_value, sigma):
    """One symmetric second-order step; returns the advanced ``s``."""
    s_mid = s_value + 0.5 * sigma
    htilde = float(chain.curvature_frame(np.array([s_mid]))[0])
    quarter = 0.25 * sigma
    half = 0.5 * sigma
    _h1_rotation(state, htilde, quarter)
    _h3_kick(chain, rigidity, state, htilde, s_mid, half)
    _h1_rotation(state, htilde, quarter)
    _h2_exact(chain, rigidity, state, s_mid, sigma)
    _h1_rotation(state, htilde, quarter)
    _h3_kick(chain, rigidity, state, htilde, s_mid, half)
    _h1_rotation(state, htilde, quarter)
    return s_value + sigma


def track_yoshida_split(chain, rigidity, state, s_span, *, step_count,
                        order=2):
    """Track the on-momentum canonical 4D state through the chain field.

    ``state = (x, px, y, py)`` with the same canonical convention as the
    chain A-RK (``py`` includes ``atilde_y``); ``s_span = (s0, s1)``;
    ``step_count`` uniform steps of the order-2 palindrome, or the
    order-4 Yoshida triple ``S2(w1 h) S2(w0 h) S2(w1 h)`` per step.
    """
    order = int(order)
    if order not in (2, 4):
        raise ValueError("order must be 2 or 4")
    count = int(step_count)
    if count < 1:
        raise ValueError("step_count must be positive")
    value = np.asarray(state, dtype=float).reshape(-1).copy()
    if value.size != 4 or not np.all(np.isfinite(value)):
        raise ValueError("state must be four finite canonical values")
    s0, s1 = (float(v) for v in s_span)
    if not s1 > s0:
        raise ValueError("s_span must advance")
    sigma = (s1 - s0) / count
    rigidity = float(rigidity)
    s_value = s0
    for _ in range(count):
        if order == 2:
            s_value = _s2_step(chain, rigidity, value, s_value, sigma)
        else:
            s_value = _s2_step(chain, rigidity, value, s_value,
                               YOSHIDA4_W1 * sigma)
            s_value = _s2_step(chain, rigidity, value, s_value,
                               YOSHIDA4_W0 * sigma)
            s_value = _s2_step(chain, rigidity, value, s_value,
                               YOSHIDA4_W1 * sigma)
    return value
