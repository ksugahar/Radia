r"""Composite (laminate) residual thermal stress / rule-of-mixtures CTE -- #38.

A FREE bonded sandwich B|A|B under uniform dT develops internal axial stress from CTE
mismatch alone: sigma_i = E_i (alpha_eff - alpha_i) dT, alpha_eff the stiffness-weighted
rule of mixtures, self-equilibrated (sum A_i sigma_i = 0). The radia elasticity solver
(region-wise E, alpha thermal eigenstrain) reproduces sigma_xx per layer and the common
central strain alpha_eff*dT. Tool-independent (closed form is the reference). The FREE-bar
counterpart of the externally blocked bar (#20) / clamped plate (#35); symmetric -> no bend.
"""
import math
import os
import sys

import pytest

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.radia_ngsolve.elasticity import (rule_of_mixtures_cte,
                                                laminate_residual_stress)

EA, alA = 70e9, 23e-6
EB, alB = 200e9, 12e-6
NU, DT, L, A, B = 0.3, 100.0, 40e-3, 2e-3, 1e-3   # core width A, each face B (2B=A)


def test_rule_of_mixtures_self_equilibrium():
    """Pure closed-form identities: alpha_eff bounded by the layer CTEs, and the
    free composite is self-equilibrated sum(A_i sigma_i) = 0 (no net axial force)."""
    a_eff = rule_of_mixtures_cte([EA, EB], [A, 2 * B], [alA, alB])
    assert alB < a_eff < alA                      # stiffness-weighted, between the two
    sA = laminate_residual_stress(EA, alA, a_eff, DT)
    sB = laminate_residual_stress(EB, alB, a_eff, DT)
    assert sA < 0.0 < sB                          # high-alpha core compressed, faces in tension
    assert abs(A * sA + 2 * B * sB) < 1e-6 * abs(A * sA)   # force balance
