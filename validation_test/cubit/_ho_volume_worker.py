"""NGSolve worker for the Cubit high-order volume convergence validation.

Run in a process separate from Cubit.  The Cubit command plugin embeds a
compact Netgen build, while this worker must use the installed NGSolve mapping
and quadrature implementation exactly as a solver would.
"""

from __future__ import annotations

import argparse
import json
import platform

import ngsolve
import radia
from radia.vol_check import check_mesh_quality
from ngsolve import BND, CF, Integrate, Mesh, TaskManager


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("vol_file")
    parser.add_argument("--integration-order", type=int, required=True)
    args = parser.parse_args()

    mesh = Mesh(args.vol_file)
    # The default quadrature order under-integrates curved elements from p = 3
    # on and made the convergence stall there.  Order 3p integrates an order-p
    # volume map exactly (radia.vol_check._measure_orders); areas integrate a
    # square root and take at least order 16.
    p = max(1, int(mesh.GetCurveOrder()))
    with TaskManager():
        volume = float(Integrate(CF(1), mesh, order=3 * p))
        area = float(Integrate(CF(1), mesh, BND, order=max(16, 3 * p)))
    quality = check_mesh_quality(
        mesh,
        min_scaled_jacobian=1.0e-6,
        integration_order=args.integration_order,
    )
    print(json.dumps({
        "versions": {
            "ngsolve": ngsolve.__version__,
            "radia": radia.__version__,
            "python": platform.python_version(),
        },
        "ngsolve_volume": volume,
        "ngsolve_area": area,
        "quality": quality,
    }))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
