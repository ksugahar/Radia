"""Write deterministic SciPy acoustic references for MATLAB validation."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy.io import savemat

REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPOSITORY_ROOT / "src"))

import radia.acoustics as acoustic
from radia.acoustics.cq import (
    _cq_grid,
    bdf_delta,
    soft_sphere_scattering_complex_k,
)


def build_reference():
    points = np.array(
        [[0.0, 0.0, 1.4], [1.3, 0.0, -0.2], [-0.4, 1.1, 0.7]],
        dtype=float,
    )
    mixed_points = np.vstack(([[0.0, 0.0, 0.35]], points))
    soft = acoustic.soft_sphere_scattering(2.3, 0.9, points, terms=18)
    rigid = acoustic.rigid_sphere_scattering(2.3, 0.9, points, terms=18)
    fluid = acoustic.fluid_sphere_scattering(
        2.3,
        0.9,
        mixed_points,
        interior_wavenumber=1.4,
        density_ratio=1.75,
        terms=18,
    )
    elastic = acoustic.elastic_sphere_scattering(
        2.3,
        0.9,
        points,
        longitudinal_speed=2.4,
        shear_speed=1.1,
        density_ratio=1.8,
        terms=15,
    )
    complex_wavenumber = 0.75 + 0.45j
    complex_scattered = soft_sphere_scattering_complex_k(
        complex_wavenumber, 0.9, points, terms=18
    )
    zeta = np.array([[0.1 + 0.2j, -0.4 + 0.3j, 0.25 - 0.15j]])
    grids = {
        "bdf1": _cq_grid(15, 0.08, 1.2, "BDF1"),
        "bdf2": _cq_grid(16, 0.08, 1.2, "BDF2"),
    }
    result = {
        "points": points,
        "mixed_points": mixed_points,
        "soft_scattered": soft["scattered"].reshape(-1, 1),
        "soft_total": soft["total"].reshape(-1, 1),
        "rigid_scattered": rigid["scattered"].reshape(-1, 1),
        "fluid_total": fluid["total"].reshape(-1, 1),
        "fluid_inside_mask": fluid["inside_mask"].reshape(-1, 1),
        "elastic_scattered": elastic["scattered"].reshape(-1, 1),
        "complex_wavenumber": complex_wavenumber,
        "complex_scattered": complex_scattered.reshape(-1, 1),
        "zeta": zeta,
        "bdf1": np.asarray(bdf_delta(zeta, "BDF1")),
        "bdf2": np.asarray(bdf_delta(zeta, "BDF2")),
    }
    for name, grid in grids.items():
        result[f"{name}_cq_radius"] = grid["cq_radius"]
        result[f"{name}_cq_zeta"] = np.asarray(grid["zeta"]).reshape(1, -1)
        result[f"{name}_cq_nodes"] = np.asarray(grid["cq_nodes"]).reshape(1, -1)
        result[f"{name}_cq_wavenumbers"] = np.asarray(
            grid["cq_wavenumbers"]
        ).reshape(1, -1)
    return result


def main():
    if len(sys.argv) != 2:
        raise SystemExit("usage: acoustic_python_reference.py OUTPUT.mat")
    output = Path(sys.argv[1]).resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    savemat(output, build_reference(), do_compression=False, oned_as="column")


if __name__ == "__main__":
    main()
