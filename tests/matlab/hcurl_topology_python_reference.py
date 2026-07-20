"""Generate deterministic Python references for the HCurl topology MEX test."""

from __future__ import annotations

from pathlib import Path
import sys

import numpy as np
from scipy.io import savemat


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

import radia._radia_pybind as _rp  # noqa: E402
from radia.hcurl_topology_optimization import (  # noqa: E402
    HCurlConductivityInterpolation,
    HCurlJouleLoadCase,
    linearize_hcurl_multifrequency_activation_joule_loss,
    linearize_hcurl_multifrequency_joule_loss,
)
from radia.vim._hcurl_tet_interaction import HCurlHMatrixOperator  # noqa: E402
from radia.vim._vim import _f64_buffer, _i32_buffer, _outer_tet  # noqa: E402


def _case_arrays(result):
    return {
        "state": np.column_stack([case.state for case in result.cases]),
        "adjoint": np.column_stack([case.adjoint for case in result.cases]),
        "case_objective": np.asarray([case.objective for case in result.cases]),
        "case_gradient": np.column_stack(
            [case.gradient for case in result.cases]
        ),
        "objective": result.objective,
        "gradient": result.gradient,
        "weights": result.weights,
    }


def build_reference():
    cell_vertices = np.array([
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0],
         [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        [[1.4, 0.2, 0.1], [2.2, 0.3, 0.1],
         [1.4, 1.1, 0.2], [1.5, 0.3, 0.8]],
    ])
    charge_hosts = np.array([0, 0, 1, 1], dtype=np.int32)
    host_parents = np.array([0, 1], dtype=np.int32)
    polynomial_coefficients = np.array([
        [1.0, 0.1, 0.0, 0.0],
        [0.0, 1.0, 0.2, 0.0],
        [1.0, -0.1, 0.3, 0.0],
        [0.2, 0.0, 0.5, 1.0],
    ])
    polynomial_exponents = np.array(
        [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]],
        dtype=np.int32,
    )
    quadrature_points, quadrature_weights = _outer_tet(4)
    aca_eps, leaf_size, eta = 1.0e-12, 4, 2.0
    gram = _rp._ChargeGramHMatrix.from_local_polynomials(
        cell_verts=_f64_buffer(cell_vertices),
        n_el=2,
        charge_host=_i32_buffer(charge_hosts),
        polynomial_coefficients=_f64_buffer(polynomial_coefficients),
        polynomial_exponents=_i32_buffer(polynomial_exponents),
        ref_tet_pts=_f64_buffer(quadrature_points),
        ref_tet_w=_f64_buffer(quadrature_weights),
        eps=aca_eps,
        leaf=leaf_size,
        eta=eta,
        build=True,
    )
    charge_maps = np.array([
        [[1.0, 0.2, -0.1], [0.1, 0.4, 0.3],
         [0.2, -0.3, 0.6], [0.4, 0.1, 0.2]],
        [[0.3, -0.1, 0.2], [0.2, 0.7, -0.2],
         [-0.4, 0.2, 0.1], [0.1, 0.5, 0.3]],
        [[-0.2, 0.5, 0.1], [0.6, 0.1, -0.3],
         [0.3, 0.4, -0.2], [-0.1, 0.2, 0.7]],
    ])
    mu = 1.3
    operator = HCurlHMatrixOperator(
        gram, charge_maps, mu=mu, cell_vertices=cell_vertices,
        charge_hosts=charge_hosts, host_parents=host_parents,
    )

    x = np.array([
        [1.0+0.2j, -0.3+0.4j],
        [-0.4+0.3j, 0.7-0.1j],
        [0.2-0.5j, 0.1+0.6j],
    ])
    left = np.array([0.2-0.1j, 0.7+0.3j, -0.4+0.2j])
    right = np.array([0.8+0.2j, -0.3+0.5j, 0.1-0.6j])
    velocity = np.empty((2, 2, 4, 3))
    affine_gradients = (
        np.eye(3),
        np.array([[0.2, 0.1, -0.05], [-0.1, 0.3, 0.2],
                  [0.15, -0.2, 0.1]]),
    )
    for direction, gradient in enumerate(affine_gradients):
        velocity[direction] = cell_vertices@gradient.T

    activation = np.array([0.58, 0.81])
    activation_power = 1.4
    resistance = np.array([
        [2.1, 0.2, -0.1], [0.2, 1.7, 0.15], [-0.1, 0.15, 1.4]
    ])
    resistance_jacobian = np.array([
        [[-0.8, 0.10, 0.03], [0.10, -0.5, -0.04],
         [0.03, -0.04, -0.3]],
        [[0.15, -0.02, 0.05], [-0.02, -0.1, 0.06],
         [0.05, 0.06, 0.2]],
    ])
    frequencies = np.array([73.0, 410.0])
    weights = np.array([0.4, 1.7])
    rhs = np.array([
        [1.0+0.2j, 0.6-0.1j],
        [-0.3+0.1j, 0.2+0.4j],
        [0.25-0.2j, -0.5+0.3j],
    ])
    rhs_jacobian = np.array([
        [[0.02+0.01j, -0.01+0.03j],
         [-0.03+0.02j, 0.04-0.01j],
         [0.01-0.02j, 0.02+0.01j]],
        [[-0.01+0.04j, 0.03+0.02j],
         [0.02-0.01j, -0.02+0.05j],
         [0.03+0.01j, -0.01-0.03j]],
    ])
    cases = tuple(HCurlJouleLoadCase(
        frequencies[index], rhs[:, index], weights[index],
        rhs_jacobian[:, :, index],
    ) for index in range(len(frequencies)))
    shape_result = linearize_hcurl_multifrequency_joule_loss(
        inductance=operator, resistance=resistance,
        resistance_jacobian=resistance_jacobian,
        cell_vertex_velocities=velocity, load_cases=cases,
    )

    cell_curl_grams = np.array([
        [[1.8, 0.15, -0.05], [0.15, 1.1, 0.08],
         [-0.05, 0.08, 0.9]],
        [[0.7, -0.04, 0.12], [-0.04, 1.4, 0.05],
         [0.12, 0.05, 1.2]],
    ])
    conductivity = HCurlConductivityInterpolation(5.0, 0.4, 2.3)
    activation_result = linearize_hcurl_multifrequency_activation_joule_loss(
        inductance=operator, cell_curl_grams=cell_curl_grams,
        activation=activation, load_cases=cases, conductivity=conductivity,
        inductance_power=activation_power,
    )

    reference = {
        "cell_vertices": cell_vertices,
        "gram_cell_vertices": cell_vertices.reshape(-1, 3),
        "charge_hosts": charge_hosts.reshape(-1, 1),
        "host_parents": host_parents.reshape(-1, 1),
        "polynomial_coefficients": polynomial_coefficients,
        "polynomial_exponents": polynomial_exponents,
        "quadrature_points": quadrature_points,
        "quadrature_weights": quadrature_weights.reshape(-1, 1),
        "aca_eps": aca_eps,
        "leaf_size": leaf_size,
        "eta": eta,
        "charge_maps": charge_maps,
        "mu": mu,
        "x": x,
        "left": left.reshape(-1, 1),
        "right": right.reshape(-1, 1),
        "velocity": velocity,
        "activation": activation.reshape(-1, 1),
        "activation_power": activation_power,
        "resistance": resistance,
        "resistance_jacobian": resistance_jacobian,
        "frequencies": frequencies.reshape(-1, 1),
        "weights": weights.reshape(-1, 1),
        "rhs": rhs,
        "rhs_jacobian": rhs_jacobian,
        "cell_curl_grams": cell_curl_grams,
        "conductivity_solid": conductivity.solid,
        "conductivity_void": conductivity.void,
        "conductivity_power": conductivity.power,
        "python_dense": operator.to_dense(),
        "python_matvec": operator.matmat(x),
        "python_directional": operator.directional_contractions(
            velocity, left, right
        ).reshape(-1, 1),
        "python_activation_dense": operator.activation_to_dense(
            activation, power=activation_power
        ),
        "python_activation_matvec": np.column_stack([
            operator.activation_matvec(
                activation, x[:, column], power=activation_power
            ) for column in range(x.shape[1])
        ]),
        "python_activation_contractions": operator.activation_contractions(
            activation, left, right, power=activation_power
        ).reshape(-1, 1),
        "python_shape_resistance": shape_result.cases[0].resistance,
        "python_shape_resistance_jacobian": (
            shape_result.cases[0].resistance_jacobian
        ),
        "python_shape_inductance": operator.to_dense(),
        "python_activation_resistance": activation_result.cases[0].resistance,
        "python_activation_inductance": operator.activation_to_dense(
            activation, power=activation_power
        ),
    }
    reference.update({f"python_shape_{key}": value
                      for key, value in _case_arrays(shape_result).items()})
    reference.update({f"python_activation_{key}": value
                      for key, value in _case_arrays(activation_result).items()})
    return reference


def main():
    if len(sys.argv) != 2:
        raise SystemExit("usage: hcurl_topology_python_reference.py OUTPUT.mat")
    output = Path(sys.argv[1]).resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    savemat(output, build_reference(), do_compression=False, oned_as="column")


if __name__ == "__main__":
    main()
