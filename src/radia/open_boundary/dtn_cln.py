# -*- coding: utf-8 -*-
"""Cauer continued-fraction realisation of the separable eddy DtN.

Scheduled for removal with the rest of the repository's Cauer-ladder code.  The
exact DtN symbols, their poles and the passive sqrt(s) realisation live in
`radia.open_boundary.dtn_exact`, which does not depend on this module.
"""
import numpy as np
import numpy.polynomial.polynomial as _P

from .dtn_exact import eddy_dtn_rational_q

__all__ = ["cauer_ladder", "eval_ladder"]


def cauer_ladder(n):
    """Cauer continued-fraction stages of the eddy DtN in q = R0 sqrt(s*mu_sigma).
    EXACT at n+1 stages.  Returns a list of stage polynomials (each a short
    ascending-power ndarray in q); evaluate with eval_ladder."""
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


def eval_ladder(stages, s, R0=1.0, mu_sigma=1.0):
    """Evaluate a Cauer ladder (from cauer_ladder) at the Laplace variable s.
    Reproduces eddy_dtn to machine precision (the ladder IS the exact operator)."""
    q = R0 * np.sqrt(complex(s) * mu_sigma)
    val = None
    for stage in reversed(stages):
        sv = sum(stage[k] * q ** k for k in range(len(stage)))
        val = sv if val is None else sv + 1.0 / val
    return val
