"""Native analytic sphere-scattering references.

The numerical partial-wave kernel lives in ``src/acoustic`` and is shared by
the Python pybind11 and MATLAB MEX adapters.  The convention is outgoing
``exp(+i k r)``; fluid sound speed and density are normalized to one.
"""

from __future__ import annotations

import numpy as np

from radia import _radia_pybind as _native


def _points(values):
    return np.ascontiguousarray(np.asarray(values, dtype=float).reshape(-1, 3))


def soft_sphere_scattering(wavenumber, radius, points, terms=None):
    """Plane-wave scattering by a sound-soft sphere."""
    return _native._AcousticSoftSphere(
        float(wavenumber), float(radius), _points(points),
        -1 if terms is None else int(terms),
    )


def rigid_sphere_scattering(wavenumber, radius, points, terms=1):
    """Plane-wave scattering by a sound-hard sphere."""
    return _native._AcousticRigidSphere(
        float(wavenumber), float(radius), _points(points), int(terms),
    )


def fluid_sphere_scattering(
    wavenumber,
    radius,
    points,
    interior_wavenumber=None,
    density_ratio=1.0,
    terms=None,
):
    """Anderson plane-wave transmission by a penetrable fluid sphere."""
    interior = float(wavenumber) if interior_wavenumber is None else float(
        interior_wavenumber
    )
    return _native._AcousticFluidSphere(
        float(wavenumber), float(radius), _points(points), interior,
        float(density_ratio), -1 if terms is None else int(terms),
    )


def elastic_sphere_scattering(
    wavenumber,
    radius,
    points,
    longitudinal_speed=2.0,
    shear_speed=1.0,
    density_ratio=1.5,
    terms=0,
):
    """Faran plane-wave scattering by an elastic solid sphere."""
    return _native._AcousticElasticSphere(
        float(wavenumber), float(radius), _points(points),
        float(longitudinal_speed), float(shear_speed), float(density_ratio),
        int(terms),
    )
