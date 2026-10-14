"""Polynomial exactness of the curved BDM2 symmetric outer rules."""

import math

import numpy as np

from radia.vim._vim import _SYM10_TET, _SYM10_TRI


def test_degree10_tetrahedron_rule_integrates_all_monomials():
    points, weights = _SYM10_TET
    for i in range(11):
        for j in range(11-i):
            for k in range(11-i-j):
                actual = np.sum(weights*points[:, 0]**i*points[:, 1]**j*points[:, 2]**k)
                exact = math.factorial(i)*math.factorial(j)*math.factorial(k)
                exact /= math.factorial(i+j+k+3)
                assert abs(actual-exact) < 2e-16


def test_degree10_triangle_rule_integrates_all_monomials():
    points, weights = _SYM10_TRI
    for i in range(11):
        for j in range(11-i):
            actual = np.sum(weights*points[:, 0]**i*points[:, 1]**j)
            exact = math.factorial(i)*math.factorial(j)/math.factorial(i+j+2)
            assert abs(actual-exact) < 2e-16
