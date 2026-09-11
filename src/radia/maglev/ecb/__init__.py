"""Eddy-current brake (ECB) force computation + frequency response."""

from .lorentz import (
    compute_lorentz_force_result_via_foster,
    compute_lorentz_force_result_via_hcurl_vim,
    compute_lorentz_force_torque_via_hcurl_vim,
    compute_lorentz_force_via_direct_scalar,
    compute_lorentz_force_via_foster,
    compute_lorentz_force_via_foster_verified,
    compute_lorentz_force_via_hcurl_vim,
    pm_field_xz_dipole,
    pm_field_z_dipole,
)
from .plate_response import (
    find_drag_peak,
    find_lift_crossover,
    sweep_alpha_response,
)

__all__ = [
    "compute_lorentz_force_result_via_foster",
    "compute_lorentz_force_result_via_hcurl_vim",
    "compute_lorentz_force_torque_via_hcurl_vim",
    "compute_lorentz_force_via_direct_scalar",
    "compute_lorentz_force_via_foster",
    "compute_lorentz_force_via_foster_verified",
    "compute_lorentz_force_via_hcurl_vim",
    "find_drag_peak",
    "find_lift_crossover",
    "pm_field_xz_dipole",
    "pm_field_z_dipole",
    "sweep_alpha_response",
]
