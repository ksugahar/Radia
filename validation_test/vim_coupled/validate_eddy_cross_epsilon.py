"""Separate sampled cross-kernel errors using the spherical mean-value theorem.

These are constant vector-density kernel fixtures, not admissible closed eddy
currents or tangential SIBC modes. They validate the native integration route,
not electromagnetic boundary conditions. The disjoint supports have no
singular pairs; touching/self-panel integration is a separate acceptance gate.
"""

import argparse
import hashlib
import json
import platform
import time
from pathlib import Path

import ngsolve as ng
import numpy as np

from radia import _radia_pybind as native
from radia.vim import _eddy_hybrid as eddy


def spherical_basis(order, kind, center, radius=0.5):
    """Tensor Gauss radial/polar rule and periodic azimuthal rule."""
    z, wz = np.polynomial.legendre.leggauss(order)
    phi = np.arange(2 * order) * np.pi / order
    zz, pp = np.meshgrid(z, phi, indexing="ij")
    directions = np.column_stack((
        (np.sqrt(1 - zz**2) * np.cos(pp)).ravel(),
        (np.sqrt(1 - zz**2) * np.sin(pp)).ravel(), zz.ravel(),
    ))
    angular_weights = np.repeat(wz, len(phi)) * np.pi / order
    if kind == "volume":
        r = radius * (z + 1) / 2
        wr = wz * radius / 2 * r**2
    elif kind == "surface":
        r, wr = np.array([radius]), np.array([radius**2])
    else:
        raise ValueError("kind must be volume or surface")
    points = (r[:, None, None] * directions).reshape(-1, 3) + center
    weights = (wr[:, None] * angular_weights).ravel()
    modes = np.zeros((1, len(points), 3))
    modes[0, :, 0] = 1.0
    return eddy.SampledCurrentBasis(
        points=points, weights=weights, modes=modes, kind=kind,
        names=(kind,),
    )


def run(orders=(3, 5), epsilons=(0.3, 0.1, 0.03), aca_tolerances=(1e-6, 1e-10)):
    records = []
    radius, separation = 0.5, 3.0
    measures = {"volume": 4 * np.pi * radius**3 / 3,
                "surface": 4 * np.pi * radius**2}
    started = time.perf_counter()
    for kinds in (("volume", "volume"), ("volume", "surface"),
                  ("surface", "surface")):
        exact = measures[kinds[0]] * measures[kinds[1]] / separation
        for order in orders:
            left = spherical_basis(order, kinds[0], np.zeros(3), radius)
            right = spherical_basis(order, kinds[1], np.array([separation, 0, 0]), radius)
            delta = left.points[:, None, :] - right.points[None, :, :]
            distance2 = np.sum(delta**2, axis=2)
            weights = left.weights[:, None] * right.weights[None, :]
            unregularized = float(np.sum(weights / np.sqrt(distance2)))
            for epsilon in epsilons:
                sampled = float(np.sum(weights / np.sqrt(distance2 + epsilon**2)))
                settings = [("compressed", aca, 8) for aca in aca_tolerances]
                settings.append(("dense-leaf", min(aca_tolerances),
                                 left.n_samples + right.n_samples + 1))
                for compression, aca, leaf_size in settings:
                    # mu=4*pi removes the physical prefactor for this kernel test.
                    backend = eddy.HACApKSampledLaplaceInteraction(
                        mu=4 * np.pi, kernel_epsilon=epsilon, aca_eps=aca,
                        leaf_size=leaf_size, cross_only=False,
                    )
                    with ng.TaskManager():
                        operator = backend.build_operator((left, right))
                        matrix = operator.to_dense()
                    stats = operator.stats()["hmatrix"]
                    observed = float(matrix[0, 1])
                    records.append({
                        "pair": "-".join(kinds), "quadrature_order": int(order),
                        "epsilon_m": float(epsilon), "aca_tolerance": float(aca),
                        "compression": compression, "leaf_size": leaf_size,
                        "lowrank_leaf_count": int(stats["n_lowrank"]),
                        "dense_leaf_count": int(stats["n_dense"]),
                        "sample_count": left.n_samples + right.n_samples,
                        "exact_unregularized_cross": float(exact),
                        "native_cross": observed,
                        "quadrature_relative_error": (unregularized - exact) / exact,
                        "epsilon_shift_relative": (sampled - unregularized) / exact,
                        "native_minus_sampled_relative": (observed - sampled) / exact,
                        "total_relative_error": (observed - exact) / exact,
                        "reciprocity_relative_error": float(abs(matrix[0, 1] - matrix[1, 0]) / exact),
                    })
    files = {"validator": Path(__file__), "python_operator": Path(eddy.__file__),
             "native_extension": Path(native.__file__)}
    return {
        "schema": "radia.eddy-cross-epsilon-validation.v1",
        "scope": "disjoint constant-density Laplace cross blocks only",
        "certifies_hdiv_sibc_coupled_solver": False,
        "runtime": {"hostname": platform.node(), "python": platform.python_version(),
                    "ngsolve": ng.__version__},
        "provenance": {key: {"path": str(path.resolve()),
                             "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
                       for key, path in files.items()},
        "records": records, "elapsed_seconds": time.perf_counter() - started,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = run()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(f"Recorded {len(result['records'])} native cross-block evaluations in "
          f"{result['elapsed_seconds']:.3f} s: {args.output}")


if __name__ == "__main__":
    main()
