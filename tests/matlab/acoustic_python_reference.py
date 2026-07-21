"""Write deterministic Python acoustic results for MATLAB/MEX adapter parity."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy.io import savemat

import radia.acoustics as acoustic
from radia.acoustics.cq import (
    bdf_delta,
    soft_sphere_scattering_complex_k,
)
from radia import _radia_pybind as native


def build_reference():
    points = np.array(
        [[0.0, 0.0, 1.4], [1.3, 0.0, -0.2], [-0.4, 1.1, 0.7]],
        dtype=float,
    )
    mixed_points = np.vstack(([[0.0, 0.0, 0.35]], points))
    soft = acoustic.soft_sphere_scattering(2.3, 0.9, points, terms=18)
    rigid = acoustic.rigid_sphere_scattering(2.3, 0.9, points, terms=18)
    fluid = acoustic.fluid_sphere_scattering(
        2.3, 0.9, mixed_points, interior_wavenumber=1.4,
        density_ratio=1.75, terms=18,
    )
    elastic = acoustic.elastic_sphere_scattering(
        2.3, 0.9, points, longitudinal_speed=2.4, shear_speed=1.1,
        density_ratio=1.8, terms=15,
    )
    complex_wavenumber = 0.75 + 0.45j
    complex_scattered = soft_sphere_scattering_complex_k(
        complex_wavenumber, 0.9, points, terms=18
    )
    zeta = np.array([[0.1 + 0.2j, -0.4 + 0.3j, 0.25 - 0.15j]])
    grids = {
        "bdf1": native._AcousticCQGrid(15, 0.08, 1.2, "BDF1"),
        "bdf2": native._AcousticCQGrid(16, 0.08, 1.2, "BDF2"),
    }
    return {
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
        **{
            f"{name}_cq_radius": grid["cq_radius"]
            for name, grid in grids.items()
        },
        **{
            f"{name}_cq_zeta": np.asarray(grid["zeta"]).reshape(1, -1)
            for name, grid in grids.items()
        },
        **{
            f"{name}_cq_nodes": np.asarray(grid["cq_nodes"]).reshape(1, -1)
            for name, grid in grids.items()
        },
        **{
            f"{name}_cq_wavenumbers": np.asarray(
                grid["cq_wavenumbers"]
            ).reshape(1, -1)
            for name, grid in grids.items()
        },
    }


def main():
    if len(sys.argv) != 2:
        raise SystemExit("usage: acoustic_python_reference.py OUTPUT.mat")
    output = Path(sys.argv[1]).resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    savemat(output, build_reference(), do_compression=False, oned_as="column")


if __name__ == "__main__":
    main()
