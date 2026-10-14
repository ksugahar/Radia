"""NGSolve worker for the Cubit high-order volume convergence validation.

Run in a process separate from Cubit.  The Cubit command plugin embeds a
compact Netgen build, while this worker must use the installed NGSolve mapping
and quadrature implementation exactly as a solver would.
"""

from __future__ import annotations

import argparse
import json
import platform

import cubit_mesh_export
import ngsolve
from cubit_mesh_export.check import check_mesh_quality
from ngsolve import BND, CF, Integrate, Mesh, TaskManager


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("vol_file")
    parser.add_argument("--integration-order", type=int, required=True)
    args = parser.parse_args()

    mesh = Mesh(args.vol_file)
    with TaskManager():
        volume = float(Integrate(CF(1), mesh))
        area = float(Integrate(CF(1), mesh, BND))
    quality = check_mesh_quality(
        mesh,
        min_scaled_jacobian=1.0e-6,
        integration_order=args.integration_order,
    )
    print(json.dumps({
        "versions": {
            "ngsolve": ngsolve.__version__,
            "cubit_mesh_export": cubit_mesh_export.__version__,
            "python": platform.python_version(),
        },
        "ngsolve_volume": volume,
        "ngsolve_area": area,
        "quality": quality,
    }))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
