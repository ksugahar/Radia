"""Eddy-current brake (ECB) force computation + frequency response."""

from .lorentz import (
    compute_lorentz_force_via_foster,
    pm_field_z_dipole,
    pm_field_xz_dipole,
)
from .plate_response import (
    sweep_alpha_response,
    find_drag_peak,
    find_lift_crossover,
)

__all__ = [
    "compute_lorentz_force_via_foster",
    "pm_field_z_dipole",
    "pm_field_xz_dipole",
    "sweep_alpha_response",
    "find_drag_peak",
    "find_lift_crossover",
]
