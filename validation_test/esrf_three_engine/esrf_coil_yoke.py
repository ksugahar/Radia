"""Shared physical contracts for the coil-driven ESRF yoke validations.

Examples 6 and 7 are nonlinear quadrupoles.  Their response mesh is the
iron-only Cubit mesh used by HDiv-MMM, while reduced-A and mixed
total/reduced-Omega use a second, conforming physical-air plus Kelvin mesh.
The current source must not be meshed: all three formulations receive the
same solid-current ``CoilBuilder`` geometry through one Radia object tree.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class ESRFCoilYokeCase:
    """Geometry-independent source and observation contract for one case."""

    number: int
    slug: str
    beam_axis: int
    kelvin_radius_m: float
    core_half_length_m: float
    transverse_offsets_m: tuple[float, ...]
    axial_stations_m: tuple[float, ...]
    gap_refinement_radius_m: float
    gap_refinement_half_length_m: float
    iron_size_m: float
    outer_air_size_m: float
    kelvin_size_m: float


_CASES = {
    6: ESRFCoilYokeCase(
        number=6,
        slug="quadrupole",
        # The Example-6 pole is an x extrusion; x is the beam direction.
        beam_axis=0,
        kelvin_radius_m=0.16,
        core_half_length_m=0.020,
        transverse_offsets_m=(-0.006, 0.0, 0.006),
        axial_stations_m=(-0.024, -0.012, 0.0, 0.012, 0.024),
        gap_refinement_radius_m=0.026,
        gap_refinement_half_length_m=0.032,
        iron_size_m=0.006,
        outer_air_size_m=0.016,
        kelvin_size_m=0.035,
    ),
    7: ESRFCoilYokeCase(
        number=7,
        slug="esrf_storage_ring_quadrupole",
        # The Example-7 yoke is extruded in z; z is the beam direction.
        beam_axis=2,
        kelvin_radius_m=0.40,
        core_half_length_m=0.100,
        transverse_offsets_m=(-0.008, 0.0, 0.008),
        axial_stations_m=(-0.150, -0.075, 0.0, 0.075, 0.150),
        gap_refinement_radius_m=0.043,
        gap_refinement_half_length_m=0.180,
        iron_size_m=0.010,
        outer_air_size_m=0.040,
        kelvin_size_m=0.080,
    ),
}


def get_case(number: int) -> ESRFCoilYokeCase:
    """Return the fixed three-engine contract for coil-driven ESRF cases."""
    try:
        return _CASES[int(number)]
    except KeyError as exc:
        raise ValueError("only ESRF coil-yoke examples 6 and 7 are supported") from exc


def observation_points(number: int) -> np.ndarray:
    """Return a symmetric gap stencil expressed in physical metres.

    Every point lies well inside the bore rather than on a material interface.
    The shared stencil measures both the central quadrupole gradient and its
    axial variation without treating a one-sided FE face trace as a field
    value.
    """
    case = get_case(number)
    transverse_axes = tuple(axis for axis in range(3) if axis != case.beam_axis)
    points: list[list[float]] = []
    for station in case.axial_stations_m:
        for first in case.transverse_offsets_m:
            for second in case.transverse_offsets_m:
                point = [0.0, 0.0, 0.0]
                point[case.beam_axis] = float(station)
                point[transverse_axes[0]] = float(first)
                point[transverse_axes[1]] = float(second)
                points.append(point)
    return np.asarray(points, dtype=float)


def core_selector(number: int, points: np.ndarray) -> np.ndarray:
    """Select the central magnetic-length portion of a shared stencil."""
    case = get_case(number)
    values = np.asarray(points, dtype=float)
    if values.ndim != 2 or values.shape[1] != 3:
        raise ValueError("points must have shape (n, 3)")
    return np.abs(values[:, case.beam_axis]) <= case.core_half_length_m + 1.0e-14


def build_radia_coil_source(number: int) -> tuple[int, dict[str, object]]:
    """Materialize the authoritative CoilBuilder sources once.

    The returned Radia container is consumed unchanged by HDiv, HCurl, and
    mixed-Omega.  It deliberately has no volume mesh counterpart.
    """
    import radia as rad

    from radia.esrf_examples import build_esrf_coils

    case = get_case(number)
    builders = build_esrf_coils(case.number)
    if not builders:
        raise RuntimeError(f"ESRF example {case.number} has no CoilBuilder source")
    objects = []
    entries = []
    for index, builder in enumerate(builders):
        if not builder.is_closed or builder.gap > 1.0e-12:
            raise RuntimeError(
                f"ESRF example {case.number} CoilBuilder {index} is open by "
                f"{builder.gap:.6e} m"
            )
        created = builder.to_radia(arc_max_segment_length=0.004)
        if not created:
            raise RuntimeError(
                f"ESRF example {case.number} CoilBuilder {index} created no "
                "Radia current objects"
            )
        objects.extend(created)
        entries.append(
            {
                "index": int(index),
                "current_A": float(builder.current),
                "segment_count": int(len(builder.segments)),
                "closed": True,
                "closure_gap_m": float(builder.gap),
                "radia_object_count": int(len(created)),
            }
        )
    return rad.ObjCnt(objects), {
        "authority": f"ESRF Example {case.number} CoilBuilder solid-current paths",
        "case": int(case.number),
        "source_is_meshed": False,
        "coil_count": int(len(builders)),
        "radia_object_count": int(len(objects)),
        "coils": entries,
    }


__all__ = [
    "ESRFCoilYokeCase",
    "build_radia_coil_source",
    "core_selector",
    "get_case",
    "observation_points",
]
