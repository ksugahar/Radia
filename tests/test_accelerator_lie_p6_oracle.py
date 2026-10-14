"""Lock the arbitrary-jet (p=6) symbolic Lie-map reference.

``validation_test/ffag_topopt/lie_map_symbolic_oracle_p6.wls`` derives, for a
committed median-plane-symmetric transverse jet through total degree six, the
homogeneous Hamiltonians ``H2..H6`` and the flow-map contractions ``R.u``,
``T[u,u]``, ``U[u^3]``, ``V[u^4]``, and the fifth-order ``W[u^5]`` around the
design orbit, plus a one-degree-of-freedom fifth-order Dragt--Finn composition
with known generators.  Wolfram Language is a derivation/golden tool only; this
test consumes the committed JSON without invoking it.

Two independent locks:

* analytic spot checks of tensor entries that close in one line
  (``H2_y_py = -cy01``, ``H6_x^6 = -720 cs60``, ``H2_dd = 1-beta^2``);
* an exact-Hamiltonian scipy flow: derivatives of the flow through fifth
  order involve only ``H2..H6``, so
  ``flow(eps*u) - sum_k eps^k/k! contraction_k`` must close at ``O(eps^6)``.
"""

import json
import math
from pathlib import Path

import numpy as np
import pytest

scipy_integrate = pytest.importorskip("scipy.integrate")

REFERENCE = (
    Path(__file__).resolve().parents[1]
    / "validation_test"
    / "ffag_topopt"
    / "lie_map_symbolic_reference_p6.json"
)


@pytest.fixture(scope="module")
def reference():
    return json.loads(REFERENCE.read_text(encoding="utf-8"))


def _jet(coefficients):
    terms = []
    for key, value in coefficients.items():
        i_str, j_str = key[1:].split("_")
        terms.append((int(i_str), int(j_str), float(value)))
    return terms


def test_tensor_entry_spot_checks(reference):
    entries = reference["hamiltonian_tensor_entries"]
    as_terms = dict(
        ((i, j), c) for i, j, c in _jet(reference["as_jet_coefficients"])
    )
    ay_terms = dict(
        ((i, j), c) for i, j, c in _jet(reference["ay_jet_coefficients"])
    )
    beta = reference["parameters"]["reference_beta"]
    hh = reference["parameters"]["hh"]
    # Median-plane parity: a_s even in y, a_y odd in y.
    assert all(j % 2 == 0 for (_, j) in as_terms)
    assert all(j % 2 == 1 for (_, j) in ay_terms)
    # Design-orbit condition H1=0.
    assert as_terms[(1, 0)] == pytest.approx(-hh, rel=1.0e-15)
    # -py*a_y cross term at second order.
    assert entries["H2_y_py"] == pytest.approx(-ay_terms[(0, 1)], rel=1.0e-13)
    # Pure-x sixth order comes only from a_s (a_y vanishes on y=0).
    assert entries["H6_x_x_x_x_x_x"] == pytest.approx(
        -math.factorial(6) * as_terms[(6, 0)], rel=1.0e-13
    )
    # Chromatic closure of the slip term.
    assert entries["H2_delta_delta"] == pytest.approx(
        1.0 - beta**2, rel=1.0e-13
    )
    residual = reference["linear_map"]["symplectic_residual"]
    assert abs(residual) < 1.0e-20


def _exact_rhs(reference):
    hh = reference["parameters"]["hh"]
    beta = reference["parameters"]["reference_beta"]
    mass_square = 1.0 / beta**2 - 1.0
    as_terms = _jet(reference["as_jet_coefficients"])
    ay_terms = _jet(reference["ay_jet_coefficients"])

    def value(terms, x, y):
        return sum(c * x**i * y**j for i, j, c in terms)

    def d_dx(terms, x, y):
        return sum(c * i * x ** (i - 1) * y**j for i, j, c in terms if i > 0)

    def d_dy(terms, x, y):
        return sum(c * j * x**i * y ** (j - 1) for i, j, c in terms if j > 0)

    def rhs(_s, z):
        x, px, y, py, _ell, delta = z
        a_y = value(ay_terms, x, y)
        root = math.sqrt((1.0 + delta) ** 2 - px**2 - (py - a_y) ** 2)
        metric = 1.0 + hh * x
        dH_dx = (
            -hh * root
            - metric * (py - a_y) * d_dx(ay_terms, x, y) / root
            - d_dx(as_terms, x, y)
        )
        dH_dy = (
            -metric * (py - a_y) * d_dy(ay_terms, x, y) / root
            - d_dy(as_terms, x, y)
        )
        dH_ddelta = -metric * (1.0 + delta) / root + (1.0 + delta) / (
            beta * math.sqrt((1.0 + delta) ** 2 + mass_square)
        )
        return np.array(
            [
                metric * px / root,
                -dH_dx,
                metric * (py - a_y) / root,
                -dH_dy,
                -dH_ddelta,
                0.0,
            ]
        )

    return rhs


@pytest.mark.parametrize(
    "probe_key,tensor_keys",
    [
        ("probe_u", ("Ru", "Tuu", "Uuuu", "Vuuuu", "Wuuuuu")),
        ("probe_v", ("Rv", "Tvv", "Uvvv", "Vvvvv", "Wvvvvv")),
    ],
)
def test_flow_contractions_close_at_sixth_order(
    reference, probe_key, tensor_keys
):
    rhs = _exact_rhs(reference)
    length = reference["parameters"]["segment_length"]
    contract = reference["variational_map_contractions"]
    u = np.asarray(contract[probe_key], dtype=float)
    tensors = [np.asarray(contract[key], dtype=float) for key in tensor_keys]

    def defect(eps):
        prediction = sum(
            eps**k / math.factorial(k) * tensor
            for k, tensor in enumerate(tensors, start=1)
        )
        solution = scipy_integrate.solve_ivp(
            rhs,
            (0.0, length),
            eps * u,
            method="DOP853",
            rtol=1.0e-13,
            atol=1.0e-16,
            max_step=length / 50.0,
        )
        assert solution.success
        return float(np.max(np.abs(solution.y[:, -1] - prediction)))

    coarse = defect(1.0e-2)
    fine = defect(5.0e-3)
    # O(eps^6) closure: exactly 64x per halving until the float64 floor.
    assert coarse < 5.0e-14
    assert fine < max(2.0e-15, coarse / 30.0)


def test_dragt_finn_fifth_order_sample_is_committed(reference):
    sample = reference["dragt_finn_fifth_order_sample"]
    # Known generators a factorize-fifth-order implementation must recover.
    assert sample["f3"]["q2p"] == pytest.approx(0.2, rel=1.0e-15)
    assert sample["f6"]["q6"] == pytest.approx(7.0 / 78.0, rel=1.0e-13)
    q_map = sample["q_map_coefficients"]
    p_map = sample["p_map_coefficients"]
    # The linear part is the identity (R is factored out of the sample).
    assert q_map["q1p0"] == pytest.approx(1.0, rel=1.0e-15)
    assert p_map["q0p1"] == pytest.approx(1.0, rel=1.0e-15)
    # Degree-5 rows exist and are nonzero: the f6 extraction has signal.
    assert abs(p_map["q5p0"]) > 1.0e-3
    assert any(
        abs(value) > 1.0e-12
        for key, value in q_map.items()
        if sum(int(part) for part in key.replace("q", "").split("p")) == 5
    )
