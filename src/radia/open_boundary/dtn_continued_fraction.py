# -*- coding: utf-8 -*-
"""Continued-fraction form of the separable eddy DtN.

The exact eddy/diffusion DtN of multipole n is a rational function of
q = R0*sqrt(s*mu_sigma) (`dtn_exact.eddy_dtn_rational_q`).  Euclid's algorithm
on its numerator and denominator gives the regular continued fraction

    G_n(q) = c_0(q) + 1/(c_1(q) + 1/(c_2(q) + ...)),

which terminates after exactly n+1 partial quotients, so the expansion
reproduces the symbol exactly.  This is the classical continued-fraction
representation used by exact and high-order local non-reflecting boundary
conditions (Grote-Keller; Hagstrom-Warburton); it is evaluated here as an
analytic non-reflecting boundary operator.
"""
import numpy as np
import numpy.polynomial.polynomial as _P

from .dtn_exact import eddy_dtn_rational_q

__all__ = ["continued_fraction_stages", "eval_continued_fraction"]


def continued_fraction_stages(n):
    """Partial quotients of the eddy DtN of multipole n in q = R0 sqrt(s*mu_sigma).

    Exact after n+1 partial quotients.  Returns a list of polynomials in q
    (short ascending-power ndarrays); evaluate with eval_continued_fraction.
    """
    A, den = eddy_dtn_rational_q(n)
    num = np.trim_zeros(np.asarray(A, float), 'b')
    den = np.trim_zeros(np.asarray(den, float), 'b')
    stages = []
    while len(num) and np.any(np.abs(den) > 1e-13) and len(stages) < 40:
        quo, rem = _P.polydiv(num, den)
        stages.append(quo)
        num, den = den, np.trim_zeros(rem, 'b')
        if len(den) == 0:
            break
    return stages


def eval_continued_fraction(stages, s, R0=1.0, mu_sigma=1.0):
    """Evaluate the continued fraction from continued_fraction_stages at s.

    Reproduces eddy_dtn to machine precision, because the expansion of the
    rational symbol terminates.
    """
    q = R0 * np.sqrt(complex(s) * mu_sigma)
    val = None
    for stage in reversed(stages):
        sv = sum(stage[k] * q ** k for k in range(len(stage)))
        val = sv if val is None else sv + 1.0 / val
    return val
