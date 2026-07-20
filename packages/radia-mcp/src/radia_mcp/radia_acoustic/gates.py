"""Cheap pre-solve gates for FSI and convolution quadrature."""

from __future__ import annotations

import numpy as np


def fsi_preflight_gate(*, wavenumber: float, c_longitudinal: float = 2.0,
                       c_transverse: float = 1.0, solid_density: float = 1.5,
                       fluid_density: float = 1.0, order: int = 1,
                       boundary: str = "gamma", radius_deviation: float = 0.0,
                       radius_tolerance: float = 0.03) -> dict:
    checks = {
        "positive_wavenumber": wavenumber > 0,
        "positive_wave_speeds": c_longitudinal > 0 and c_transverse > 0,
        "positive_densities": solid_density > 0 and fluid_density > 0,
        "positive_lame_lambda": c_longitudinal**2 > 2 * c_transverse**2,
        "valid_order": isinstance(order, int) and order >= 1,
        "boundary_named": bool(boundary.strip()),
        "spherical_dtn_geometry": 0 <= radius_deviation <= radius_tolerance,
    }
    return {
        "schema": "radia-mcp.radia-acoustic-fsi-preflight/v1",
        "ok": all(checks.values()), "checks": checks,
        "formulation": "NGSolve VectorH1 elasticity + spherical Helmholtz DtN",
        "next": "call radia.acoustics.fsi.fsi_dtn_solve only when ok is true",
    }


def cq_grid_gate(*, num_time: int, time_step: float, sound_speed: float = 1.0,
                 method: str = "BDF2") -> dict:
    from radia.acoustics.cq import bdf_delta
    method = method.upper()
    basic = isinstance(num_time, int) and num_time >= 4 and time_step > 0 and sound_speed > 0 and method in {"BDF1", "BDF2"}
    if not basic:
        return {"schema": "radia-mcp.radia-acoustic-cq-grid/v1", "ok": False,
                "errors": ["require num_time>=4, time_step>0, sound_speed>0, method BDF1/BDF2"]}
    n = np.arange(num_time)
    rho = (np.finfo(float).eps ** 0.5) ** (1.0 / num_time)
    zeta = rho * np.exp(-2j * np.pi * n / num_time)
    s = bdf_delta(zeta, method) / time_step
    kappa = 1j * s / sound_speed
    checks = {
        "laplace_nodes_right_half_plane": bool(np.all(s.real > 0)),
        "finite_nodes": bool(np.all(np.isfinite(s)) and np.all(np.isfinite(kappa))),
        "conjugate_pair_symmetry": bool(np.max(np.abs(s[1:] - np.conj(s[:0:-1]))) < 1e-10),
    }
    return {
        "schema": "radia-mcp.radia-acoustic-cq-grid/v1", "ok": all(checks.values()),
        "checks": checks, "method": method, "num_time": num_time,
        "time_step": time_step, "sound_speed": sound_speed, "cq_radius": float(rho),
        "min_real_s": float(np.min(s.real)), "max_abs_kappa": float(np.max(np.abs(kappa))),
        "convention": "s=delta(zeta)/dt; kappa=i*s/c",
    }
